# Make predictions of Systole and Diastole for all datasets.
# Uses "segments-classified.csv" to create a model to help find LV segments, then filters and
# aggregates the segment data to have one set of meta-info per ID. Combines this with the 
# train set to create a model for the Systole and Diastole volumes.

# TODO: with image dimension x/y added we could scale the areas with
# sum(area all segments)/area(whole image). Linear attributes could be
# scaled with the square of this. Poor man's approximation of 
# Voronoi segmentation. Image dimension x/y might be part of slice meta data.

source("util.R")

library(pROC)

# Threshold for LV segment probability
pSegmentThreshold <- 0.5

# Confidence level for predictions
confidence <- 0.95


# Identify the LV segments in the whole dataset by creating a model
# from the (manually) identified ones
classifiedSegments <- fread('segments-predict.csv') 

# Read full image segmentation
allSegments <- NULL
for (dataset in datasetFolders) {
  fname <- getSegmentFile(dataset)
  if (file.exists(fname)) {
    segInfo <- fread(fname)
  }
  if (is.null(allSegments)) {
    allSegments <- segInfo
  } else {
    allSegments <- rbind(allSegments, segInfo)
  }
}
setkey(allSegments, Id, Slice, Time, UUID)

# NOTE: UUID is not really unique

# First cut at 2-D normalization of the image data
allSegments <- mutate(allSegments,
                      areaMultiplier = pixelSpacing.x * pixelSpacing.y,
                      lengthMultiplier = sqrt(areaMultiplier),
                      
                      area = s.area*areaMultiplier,
                      area.ellipse = pi*s.radius.min*s.radius.max*areaMultiplier,
                      
                      perimeter = s.perimeter*lengthMultiplier,
                      radius.mean = s.radius.mean*lengthMultiplier,
                      radius.min = s.radius.min*lengthMultiplier,
                      radius.max = s.radius.max*lengthMultiplier,
                      radius.var = sqrt(s.radius.sd)*lengthMultiplier,
                      
                      majoraxis = m.majoraxis*lengthMultiplier,
                      roundness = 4*pi*area/(perimeter^2)) %>%
  rename(
    # "m.eccentricity" and "m.theta" are scale independent
    eccentricity = m.eccentricity,
    theta = m.theta) %>%
  select(-pixelSpacing.x,
         -pixelSpacing.y,
         -areaMultiplier,
         -lengthMultiplier,
         -m.majoraxis,
         -starts_with("s."))

# First, match the ones with the same UUID
segClassificationSet <- left_join(select(classifiedSegments, Id, Slice, Time, UUID, isLV), 
                                  allSegments, 
                                  by=c("Id", "Slice", "Time", "UUID"))

# For the others, find the best segment match per each Id/Slice/Time
# TODO: not sure what happens if none of the UUID's match
if (nrow(filter(segClassificationSet, is.na(Dataset))) > 0) {
  classifiedSegmentsOlder <- left_join(select(filter(segClassificationSet, is.na(Dataset)), Id, Slice, Time, UUID, isLV),
                                       select(classifiedSegments, Id, Slice, Time, UUID, isLV, m.cx, m.cy))
  segClassificationSet <- filter(segClassificationSet, !is.na(Dataset))
  s <- unique(select(classifiedSegmentsOlder, Id, Slice, Time))
  allCandidateSegments <- left_join(s, allSegments, by=c("Id","Slice","Time"))
  prevId <- -1
  prevSlice <- -1
  for (i in seq(nrow(s))) {
    candidateSegments <- filter(allCandidateSegments, Id==s$Id[i], Slice==s$Slice[i], Time==s$Time[i]) # all segs in this image
    lv <- left_join(s[i], classifiedSegmentsOlder, by=c("Id","Slice","Time")) %>% filter(isLV) # classified LV in this image
    
    if (!(prevId == s$Id[i] & prevSlice == s$Slice[i])) {
      cat("Matching existing classification to new segmentation for Id", 
          s$Id[i], "Slice", s$Slice[i], "# images", nrow(filter(s, Id == s$Id[i], Slice == s$Slice[i])), 
          "candidates:", nrow(candidateSegments), "classified:", nrow(lv), fill=T)
      prevId <- s$Id[i]
      prevSlice <- s$Slice[i]
    }
    
    candidateSegments$isLV <- NA
    candidateSegments <- candidateSegments[,names(segClassificationSet),with=F] # make sure cols have the same order
    
    if (nrow(lv) == 1) {
      # Set LV to the segment closest to the identified one. Note: similar code in classify.R
      distToLVSeg <- sqrt((candidateSegments$m.cx - lv$m.cx)^2 + (candidateSegments$m.cy - lv$m.cy)^2)
      segLV <- candidateSegments$segIndex[which.min(distToLVSeg)]
      candidateSegments$isLV <- (candidateSegments$segIndex == segLV & distToLVSeg < 5) # abs distance threshold like in classify.R
      segClassificationSet <- rbind(segClassificationSet, candidateSegments)
    } else {
      cat("WARN:",nrow(lv),"LV segments in image", s$Time[i], fill=T)
    }
    
    segClassificationSet <- rbind(segClassificationSet, candidateSegments)
  }
}

segClassificationSet <- filter(segClassificationSet, !is.na(isLV)) # just to be sure

# quick plots to verify results

# TODO: (trigonometric) interpolation & outlier detection
# NOTE: the graphs of area vs time have missing data points as well
# as outliers (typically near zero). These should be fixed somehow through
# interpolation or other means.
# ID 403/Slice 10 shows near-zero's for time 9..16
# ID 334/Slice 15 shows missing values for time 10..15
# ID 504/Slice 17 shows a near perfect graph
# See for example:
#  http://stats.stackexchange.com/questions/63233/fourier-trigonometric-interpolation

l <- unique(select(segClassificationSet, Id, Slice))
for (i in seq(nrow(l))) {
  slice <- filter(segClassificationSet, Id==l$Id[i], Slice==l$Slice[i])
  
  plotSlice(filter(slice, isLV))
}

# Build up the data set for training and classification

# TODO maybe join in some slice info to get nr of images per slice etc and
# create some predictors out of that. For now, just dropping all meta info.

segClassificationSet <- select(segClassificationSet, 
                               -m.cx, -m.cy, -distToROI,
                               -Id, -Slice, -Time, -Dataset, -ImgType, -Offset, -segIndex)
# We keep the meta-attributes as these are useful for the further roll-up
allSegments <- select(allSegments, 
                      -m.cx, -m.cy, -distToROI)

trainData <- segClassificationSet
valSet <- sample.int(nrow(trainData), 0.20*nrow(trainData))
trainDataPredictorsOnly <- select(trainData, -UUID, -isLV)

cat("Building segment model with",length(names(trainDataPredictorsOnly)),"predictors",fill=T)
uniVariateAnalysis <- data.frame(Predictor=names(trainDataPredictorsOnly),
                                 validation=sapply(seq(length(trainDataPredictorsOnly)), 
                                                   function(i) {auc(trainData$isLV[valSet], trainDataPredictorsOnly[[i]][valSet])}),
                                 train=sapply(seq(length(trainDataPredictorsOnly)), 
                                              function(i) {auc(trainData$isLV[-valSet], trainDataPredictorsOnly[[i]][-valSet])}))
uniVariateAnalysis <- gather(uniVariateAnalysis, dataset, auc, -Predictor)
print(ggplot(uniVariateAnalysis, aes(x=Predictor, y=auc, fill=dataset))+
        geom_bar(stat="identity",position="dodge")+
        theme(axis.text.x = element_text(angle = 45, hjust=1))+
        geom_hline(yintercept=0.52,linetype="dashed")+
        ggtitle("AUC of individual predictors for segmentation model"))

leftVentricleSegmentModel <- xgboost(data = as.matrix(trainDataPredictorsOnly[-valSet]), 
                                     label = trainData$isLV[-valSet], 
                                     max.depth = 2, eta = 0.1, nround = 100,
                                     objective = "binary:logistic", missing=NaN, verbose=0)
imp_matrix <- xgb.importance(feature_names = names(trainDataPredictorsOnly), model = leftVentricleSegmentModel)
print(xgb.plot.importance(importance_matrix = imp_matrix))

# Get an idea of the accuracy. Note, it seems very high always.
probLV <- predict(leftVentricleSegmentModel, as.matrix(trainDataPredictorsOnly))
cat("AUC for validation set:", auc(trainData$isLV[valSet], probLV[valSet]), fill=T)

# Distribution of probabilities
plotSet <- group_by(data.frame(predictedProbability = cut2(probLV, g=20), # equi-weight
                               isLV = trainData$isLV,
                               isVal = (seq(nrow(trainData)) %in% valSet)), predictedProbability) %>% 
  summarise(validation = sum(isLV & isVal)/sum(isVal),
            train = sum(isLV & !isVal)/(n()-sum(isVal)),
            count = n()) %>%
  gather(dataset, probability, -count, -predictedProbability)
print(ggplot(plotSet, aes(x=predictedProbability, y=probability, fill=dataset)) + 
        geom_bar(stat="identity",position="dodge") + 
        ggtitle("Segment Prediction") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Keep data for analysis elsewhere
write.csv(select(trainData, -UUID), "segmentTrainSet.csv", row.names=F)

# Apply on full dataset
cat("Apply segment model to", nrow(allSegments), "segments", fill=T)
allSegments$pLV <- predict(leftVentricleSegmentModel, 
                           as.matrix(select(allSegments, -Id, -Slice, -Time, -Dataset, -ImgType, -Offset, -segIndex, -UUID)))

# Aggregation to Image level
imageData <- left_join(group_by(allSegments, Id, Slice, Time) %>%
                         summarise(segIndex = segIndex[which.max(pLV)]),
                       allSegments)
# plotSlice( filter(imageData, Id == 657, Slice == 9))

sliceList <- fread("slicelist.csv")
imageData <- left_join(imageData, sliceList)
pLeftVentricle <- cut(imageData$pLV,10)
# This should show that lower slice order have a higher probabilities for the left ventricle (better segmentation)
print(ggplot(imageData, aes(x=pLeftVentricle, fill=factor(SliceOrder))) + geom_histogram()+
        ggtitle("LV Probability vs Slice Order") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))

idIncomplete <- unique(imageData$Id[which(!complete.cases(imageData))])
cat(length(idIncomplete),"incomplete cases for image data, Id's:",idIncomplete,fill=T)

# Filter on segmentation threshold. This will cause some images to be
# dropped but will increase the likelihood that they're segmented 
# correctly.
sliceData <- left_join(sliceList, group_by(filter(imageData, pLV >= pSegmentThreshold), Id, Slice) %>%
                         summarise(maxArea = max(area,na.rm=T),
                                   minArea = min(area,na.rm=T),
                                   meanArea = mean(area,na.rm=T),
                                   maxPerimeter = max(perimeter,na.rm=T),
                                   minPerimeter = min(perimeter,na.rm=T),
                                   meanPerimeter = mean(perimeter,na.rm=T))) %>%
  filter(SliceOrder <= 2)

idIncomplete <- unique(sliceData$Id[which(!complete.cases(sliceData))])
cat(length(idIncomplete),"incomplete slices for slice data, Id's:",idIncomplete,fill=T)

# Plot the areas for a couple of slices
for (n in seq(20)) {
  idData <- filter(sliceData, Id==n) %>% gather(metric, area, ends_with("Area"))
  if (!any(is.na(idData))) {
    print(ggplot(idData, aes(x=Slice, y=area, colour=metric))+geom_line()+geom_point()+ggtitle(paste("Id",n)))
  }
}

# Aggregates by Id...
# TODO get more meta-data in, and what happened to slice thickness etc??
# now only two predictors...
caseData <- group_by(sliceData, Id) %>%
  summarise(maxVolume = sum(maxArea, na.rm=T),
            minVolume = sum(minArea, na.rm=T),
            maxPerim  = mean(maxPerimeter, na.rm=T),
            minPerim  = mean(minPerimeter, na.rm=T),
            SlicePct  = n()/first(SliceCount))

# Train data
trainVolumes <- fread('data/train.csv') 
caseData <- left_join(caseData, trainVolumes)

# TODO deal with missing Id's - see slicelist

# Develop model (TODO: on -validation set)
# TODO: glm or svm?
systole_m <- lm(caseData$Systole ~ ., 
                data = select(caseData, -Systole, -Diastole))
diastole_m <- lm(caseData$Diastole ~ ., 
                 data = select(caseData, -Systole, -Diastole))

resultData <- data.frame(caseData, 
                         Systole = predict(systole_m, 
                                           newdata = caseData,
                                           interval="confidence", level=confidence),
                         Diastole = predict(diastole_m, 
                                            newdata = caseData,
                                            interval="confidence", level=confidence))

cat("Missing cases:",sum(is.na(resultData$Systole.fit) | is.na(resultData$Diastole.fit)),"out of",nrow(resultData),fill=T)

clipUpper <- function(v) {
  result <- quantile(v, (confidence+1)/2)
  if (is.na(result) | result < min(v)) {
    result <- max(v)
  }
  return(result)
}
clipLower <- function(v) {
  result <- quantile(v, (1-confidence)/2)
  if (is.na(result) | result > max(v)) {
    result <- max(v)
  }
  return(result)
}
resultData$Systole.fit  <- ifelse(is.na(resultData$Systole.fit),  mean(trainVolumes$Systole), resultData$Systole.fit)
resultData$Diastole.fit <- ifelse(is.na(resultData$Diastole.fit), mean(trainVolumes$Diastole), resultData$Diastole.fit)
resultData$Systole.lwr  <- ifelse(is.na(resultData$Systole.lwr),  clipLower(trainVolumes$Systole),  resultData$Systole.lwr)
resultData$Diastole.lwr <- ifelse(is.na(resultData$Diastole.lwr), clipLower(trainVolumes$Diastole), resultData$Diastole.lwr)
resultData$Systole.upr  <- ifelse(is.na(resultData$Systole.upr),  clipUpper(trainVolumes$Systole), resultData$Systole.upr)
resultData$Diastole.upr <- ifelse(is.na(resultData$Diastole.upr), clipUpper(trainVolumes$Diastole), resultData$Diastole.upr)

# Plot the results...
plotData <- select(resultData[which(complete.cases(resultData)),], contains("stole"), Id) %>% 
  gather(Phase, Actual, -Id, -contains(".")) %>%
  mutate(Predicted = ifelse(Phase == "Diastole", Diastole.fit, Systole.fit))
print(ggplot(plotData, 
             aes(x=Actual,y=Predicted,colour=Phase))+
        geom_point()+stat_smooth(method = "lm")+
        xlim(0, 600)+ylim(0, 600)+geom_abline(slope=1,linetype="dashed")+
        ggtitle("Actual vs Predicted..."))

print("Correlations:")
print(cor(resultData$Systole, resultData$Systole.fit, use="complete.obs"))
print(cor(resultData$Diastole, resultData$Diastole.fit, use="complete.obs"))

NPROBS <- 600

# Translate a predicted fit with certain low and high bounds and a confidence interval
# to a range of probabilities according to a logit sigmoid functions.
translateToProbs <- function(n.fit, n.low, n.high) 
{
  p.low <- (1-confidence)/2
  p.high <- (confidence+1)/2
  
  # alpha, beta and mu are solutions to fitting the logit to the two data points
  # defined by n.low and n.high vs the probabilities from the confidence interval.
  mu <- log(1/p.low - 1) / log(1/p.high - 1)
  alpha <- (mu*n.high - n.low)/(mu - 1)
  beta <- log(1/p.low - 1)/(alpha - n.low)
  
  p <- round(1/(1 + exp(0 - beta*(seq(NPROBS) - alpha))),3)
  
  return(p)
}

# Report on results
results_Train <- filter(resultData, !is.na(Systole) & !is.na(Diastole))
results_Test <- filter(resultData, !(Id %in% results_Train$Id))
results_Summary <- as.data.frame(t(data.frame( train_diastole = as.vector(summary(results_Train$Diastole.fit)),
                                               test_diastole = as.vector(summary(results_Test$Diastole.fit)),
                                               train_systole  = as.vector(summary(results_Train$Systole.fit)),
                                               test_systole  = as.vector(summary(results_Test$Systole.fit)))))
names(results_Summary) <- names(summary(results_Train$Diastole.fit))
print("Predictions summary:")
print(results_Summary)

# Create submission file
createProbabilities <- function(ds)
{
  probs <- matrix(nrow = 2*nrow(ds), ncol = NPROBS)
  for (i in 1:nrow(ds)) {
    Id <- ds$Id[i]
    probs[2*i-1,] <- translateToProbs(ds$Diastole.fit[ds$Id==Id], 
                                      ds$Diastole.lwr[ds$Id==Id], 
                                      ds$Diastole.upr[ds$Id==Id])
    probs[2*i,]   <- translateToProbs(ds$Systole.fit[ds$Id==Id], 
                                      ds$Systole.lwr[ds$Id==Id], 
                                      ds$Systole.upr[ds$Id==Id])
  }
  probs <- data.frame(Id = as.vector(sapply(ds$Id,function(n){paste(n,c('Diastole','Systole'),sep="_")})), probs)
  names(probs) <- c('Id', paste("P",0:(NPROBS-1),sep=""))

  return(probs)
}

testSubmission <- createProbabilities(results_Test)
write.csv(testSubmission, "submission.csv", row.names=F)

# Plot a few of the results
ds_plot <- gather(testSubmission[1:10,],Volume,Density,-Id)
ds_plot$Volume <- as.integer(gsub("P(.*)","\\1",ds_plot$Volume))
# Id <- factor(gsub("(.*)_(.*)","\\1",ds_plot$Id))
Phase <- factor(gsub("(.*)_(.*)","\\2",ds_plot$Id))
print(ggplot(data=ds_plot, aes(x=Volume, y=Density, colour=Phase))+
        geom_line(alpha=0.5)+
        ggtitle("Submissions"))

# Report on Kaggle's CRPS score
# Can probably be done much faster by operating on the whole matrix at once
trainProbabilities <- createProbabilities(results_Train)
crps <- 0
for (i in seq(nrow(results_Train))) {
  probs1 <- as.vector(as.matrix(trainProbabilities[2*i-1,2:ncol(trainProbabilities)]))
  truth1 <- ifelse(seq(NPROBS) >= results_Train$Diastole[i], 1, 0)

  probs2 <- as.vector(as.matrix(trainProbabilities[2*i,2:ncol(trainProbabilities)]))
  truth2 <- ifelse(seq(NPROBS) >= results_Train$Systole[i], 1, 0)
  
  crps <- crps + sum((probs1 - truth1)^2, na.rm=T) + sum((probs2 - truth2)^2, na.rm=T)
}
crps <- crps/nrow(trainProbabilities)/600
cat("CRPS score on train set:", crps,fill=T)

# TODO perhaps to cross validation 
