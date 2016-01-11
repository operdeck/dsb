# Make predictions of Systole and Diastole for all datasets.
# Uses "segments-classified.csv" to create a model to help find LV segments, then filters and
# aggregates the segment data to have one set of meta-info per ID. Combines this with the 
# train set to create a model for the Systole and Diastole volumes.

source("util.R")

library(pROC)

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

# TODO maybe filter on pLV threshold (0.5 or so)

sliceData <- group_by(imageData, Id, Slice) %>%
  summarise(maxArea = max(area),
            minArea = min(area),
            meanArea = mean(area),
            maxPerimeter = max(perimeter),
            minPerimeter = min(perimeter),
            meanPerimeter = mean(perimeter),
            # what happened to slice thickness etc?
            SliceCount = first(SliceCount),
            SliceRelIndex = first(SliceRelIndex),
            SliceOrder = first(SliceOrder)) %>%
  filter(SliceOrder <= 2) # Clearly, only the middle slices make some sense

# Plot the areas for a couple of slices
for (n in seq(20)) {
  idData <- filter(sliceData, Id==n) %>% gather(metric, area, ends_with("Area"))
  print(ggplot(idData, aes(x=Slice, y=area, colour=metric))+geom_line()+geom_point()+ggtitle(paste("Id",n)))
  idData <- filter(sliceData, Id==n) %>% gather(metric, perimeter, ends_with("Perimeter"))
  print(ggplot(idData, aes(x=Slice, y=perimeter, colour=metric))+geom_line()+geom_point()+ggtitle(paste("Id",n)))
}

# Aggregates by Id...
# TODO get more meta-data in, and what happened to slice thickness etc??
# now only two predictors...
caseData <- group_by(sliceData, Id) %>%
  summarise(maxVolume = sum(maxArea),
            minVolume = sum(minArea),
            SliceCount = first(SliceCount))

# Train data
trainVolumes <- fread('data/train.csv') 
caseData <- left_join(caseData, trainVolumes)

# TODO deal with missing Id's - see slicelist

# Develop model (TODO: on -validation set)
systole_m <- lm(caseData$Systole ~ ., data = select(caseData, -Systole, -Diastole))
diastole_m <- lm(caseData$Diastole ~ ., data = select(caseData, -Systole, -Diastole))

resultData <- data.frame(caseData, 
                         Systole_pred = predict(systole_m, 
                                                newdata = caseData),
                         Diastole_pred = predict(diastole_m, 
                                                 newdata = caseData))

# Plot the results...
plotData <- select(resultData, contains("stole"), Id) %>% 
  gather(Phase, Actual, -Id, -Diastole_pred, -Systole_pred) %>%
  mutate(Predicted = ifelse(Phase == "Diastole", Diastole_pred, Systole_pred))
print(ggplot(plotData, aes(x=Actual,y=Predicted,colour=Phase))+geom_point()+stat_smooth(method = "lm")+ggtitle("Actual vs Predicted..."))

print("Correlations:")
print(cor(resultData$Systole, resultData$Systole_pred, use="complete.obs"))
print(cor(resultData$Diastole, resultData$Diastole_pred, use="complete.obs"))

#TODO: report accuracy (RMSD?) on validation set

# now, translate this to probabilities by volume
# note, predictions can be NA, negative or otherwise out of bounds,
# also make sure Systole << Diastole - otherwise just fall back to default
# also make sure all IDs are covered



