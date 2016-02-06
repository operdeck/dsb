# Make predictions of Systole and Diastole for all datasets.
# Uses "segments-classified.csv" to create a model to help find LV segments, then filters and
# aggregates the segment data to have one set of meta-info per ID. Combines this with the 
# train set to create a model for the Systole and Diastole volumes.

# Volume = sum over slices *for a certain Time*
# Group by time:
#   Sum area(t) over slices (take scale and distance into account)
# This gives a view of volume over time by ID
# Min and Max volumes are the systole/diastole values

# Additional Meta data
# patient data contains age and sex information in DICOM tags
# seems useful in the final models!

source("util.R")

library(caret)
library(pROC)
library(lattice)
require(caret)
library(DMwR) # outlier detection

# source("https://bioconductor.org/biocLite.R")
# biocLite("impute"
         
# Threshold for LV segment probability
pSegmentThreshold <- 0.2

# Confidence level for predictions
confidence <- 0.95

# Used both in segment and case prediction
validationPercentage <- 0.20

# Identify the LV segments in the whole dataset by creating a model
# from the (manually) identified ones
classifiedSegments <- fread('segments-predict.csv') 

print("Reading image meta data")
imageList <- getImageList()

print("Reading segmentation")

imagePredictFile <- "allSegments-segmentsPredicted.csv"
segmentPredictFile <- "segmentTrainSet.csv"

skipSegmentPrediction <- F

if (skipSegmentPrediction & file.exists(imagePredictFile)) {
  # Keep data if we want to skip the segmentation predict phase
  print("!! Skipping segment prediction")
  allSegments <- fread(imagePredictFile)
} else {
  
  allSegments <- NULL
  for (dataset in unique(imageList$Dataset)) {
    if (file.exists(getSegmentFile(dataset))) {
      segmentsPerDataset <- fread(getSegmentFile(dataset))
      if (is.null(allSegments)) {
        allSegments <- segmentsPerDataset
      } else {
        removedSet <- setdiff(names(allSegments), names(segmentsPerDataset))
        addedSet <- setdiff(names(segmentsPerDataset), names(allSegments))
        diffSet <- paste(c(paste("-",removedSet), paste("+",addedSet)),collapse=", ")
        if (length(removedSet) + length(addedSet) > 0) {
          print(names(allSegments))
          print(names(segmentsPerDataset))
          print(diffSet)
          stop("Datasets do not match up. Please consider removing files.")
        }
        allSegments <- rbind(allSegments, segmentsPerDataset)
      }
    }
  }
  setkey(allSegments, Id, Slice, Time, UUID)
  
  # First, match the ones with the same UUID
  segClassificationSet <- left_join(select(classifiedSegments, Id, Slice, Time, UUID, isLV), 
                                    allSegments, 
                                    by=c("Id", "Slice", "Time", "UUID"))
  
  # For the others, find the best segment match per each Id/Slice/Time
  # TODO: not sure what happens if none of the UUID's match
  newSegmentationOldClassification <- filter(segClassificationSet, is.na(m.cx))
  segClassificationSet <- filter(segClassificationSet, !is.na(m.cx))
  
  if (nrow(newSegmentationOldClassification) > 0) {
    classifiedSegmentsOlder <- left_join(select(newSegmentationOldClassification, Id, Slice, Time, UUID, isLV),
                                         select(classifiedSegments, Id, Slice, Time, UUID, isLV, m.cx, m.cy))
    classifiedImagesOlder <- unique(select(classifiedSegmentsOlder, Id, Slice, Time))
    cat(nrow(classifiedImagesOlder), "images have new segmentation to match with existing classification",fill=T)
    allCandidateSegments <- left_join(classifiedImagesOlder, allSegments, by=c("Id","Slice","Time"))
    
    prevId <- -1
    prevSlice <- -1
    for (i in seq(nrow(classifiedImagesOlder))) {
      candidateSegments <- filter(allCandidateSegments, 
                                  Id==classifiedImagesOlder$Id[i], 
                                  Slice==classifiedImagesOlder$Slice[i], 
                                  Time==classifiedImagesOlder$Time[i]) # all segs in this image
      lv <- left_join(classifiedImagesOlder[i], 
                      classifiedSegmentsOlder, 
                      by=c("Id","Slice","Time")) %>% filter(isLV) # classified LV in this image
      if (!(prevId == classifiedImagesOlder$Id[i] & prevSlice == classifiedImagesOlder$Slice[i])) {
        cat("Matching existing classification to new segmentation for Id", classifiedImagesOlder$Id[i], 
            "Slice", classifiedImagesOlder$Slice[i], 
            "#img", nrow(filter(classifiedImagesOlder, Id == classifiedImagesOlder$Id[i], Slice == classifiedImagesOlder$Slice[i])), 
            "#segs:", nrow(candidateSegments), 
            "==>", nrow(lv), fill=T)
        prevId <- classifiedImagesOlder$Id[i]
        prevSlice <- classifiedImagesOlder$Slice[i]
      }
      
      candidateSegments$isLV <- NA
      candidateSegments <- candidateSegments[,names(newSegmentationOldClassification),with=F] # make sure cols have the same order
      
      if (nrow(lv) == 1) {
        # Set LV to the segment closest to the identified one. Note: similar code in classify.R
        distToLVSeg <- sqrt((candidateSegments$m.cx - lv$m.cx)^2 + (candidateSegments$m.cy - lv$m.cy)^2)
        segLV <- candidateSegments$segIndex[which.min(distToLVSeg)]
        candidateSegments$isLV <- (candidateSegments$segIndex == segLV & distToLVSeg < 5) # abs distance threshold like in classify.R
        segClassificationSet <- rbind(segClassificationSet, candidateSegments)
      } else {
        cat("WARN:",nrow(lv),"LV segments in image", classifiedImagesOlder$Time[i], fill=T)
      }
      
      segClassificationSet <- rbind(segClassificationSet, candidateSegments)
    }
  }
  
  # Quick report on the segmentation prediction data set.
  cat("Segment predict set has",nrow(segClassificationSet),"observations with a pos rate of",sum(segClassificationSet$isLV,na.rm=T)/nrow(segClassificationSet),fill=T)
  cat("   number of Ids   :",nrow(unique(select(segClassificationSet,Id))),"with identified LV",nrow(unique(select(filter(segClassificationSet,isLV),Id))),fill=T)
  cat("   number of Slices:",nrow(unique(select(segClassificationSet,Id,Slice))),"with identified LV",nrow(unique(select(filter(segClassificationSet,isLV),Id,Slice))),fill=T)
  cat("   number of Images:",nrow(unique(select(segClassificationSet,Id,Slice,Time))),"with identified LV",nrow(unique(select(filter(segClassificationSet,isLV),Id,Slice,Time))),fill=T)
  
  # Build up prediction set by creating derived variables and dropping non-predictors
  segClassificationSetIDs <- select(segClassificationSet, Id, Slice, Time)
  segClassificationSet <- createSegmentPredictSet(filter(left_join(segClassificationSet, 
                                                                   imageList, by=c("Id","Slice","Time")), 
                                                         !is.na(isLV)))
  
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
  
  # l <- head(unique(select(segClassificationSet, Id, Slice)),5)
  # for (i in seq(nrow(l))) {
  #   slice <- filter(segClassificationSet, Id==l$Id[i], Slice==l$Slice[i])
  #   
  #   plotSlice(filter(slice, isLV))
  # }
  
  # Build up the data set for training and classification
  
  # We keep the meta-attributes as these are useful for the further roll-up
  # allSegments <- select(allSegments, 
  #                       -m.cx, -m.cy, -distToROI)
  
  valSet <- sample.int(nrow(segClassificationSet), validationPercentage*nrow(segClassificationSet))
  trainDataPredictorsOnly <- select(segClassificationSet, -isLV)
  
  cat("Building segment model with",length(names(trainDataPredictorsOnly)),"predictors",fill=T)
  
  uniVariateAnalysis <- data.frame(Predictor=names(trainDataPredictorsOnly),
                                   validation=sapply(trainDataPredictorsOnly, 
                                                     function(p) {auc(segClassificationSet$isLV[valSet], p[valSet])}),
                                   train=sapply(trainDataPredictorsOnly, 
                                                function(p) {auc(segClassificationSet$isLV[-valSet], p[-valSet])}))
  uniVariateAnalysis <- gather(uniVariateAnalysis, dataset, auc, -Predictor)
  uniVariateAnalysis$Predictor <- factor(uniVariateAnalysis$Predictor, levels=arrange(uniVariateAnalysis,-auc)$Predictor)
  print(ggplot(uniVariateAnalysis, aes(x=Predictor, y=auc, fill=dataset))+
          geom_bar(stat="identity",position="dodge")+
          theme(axis.text.x = element_text(angle = 45, hjust=1))+
          geom_hline(yintercept=0.52,linetype="dashed")+
          ggtitle("AUC of individual predictors for segmentation model"))
  
  leftVentricleSegmentModel <- xgboost(data = data.matrix(trainDataPredictorsOnly[-valSet]), 
                                       label = segClassificationSet$isLV[-valSet], 
                                       max.depth = 2, eta = 0.1, nround = 100,
                                       objective = "binary:logistic", 
                                       missing=NaN, verbose=0)
  imp_matrix <- xgb.importance(feature_names = names(trainDataPredictorsOnly), model = leftVentricleSegmentModel)
  print(xgb.plot.importance(importance_matrix = imp_matrix))
  
  # Get an idea of the accuracy. Note, it seems very high always.
  probLV <- predict(leftVentricleSegmentModel, data.matrix(trainDataPredictorsOnly))
  cat("AUC for validation set:", auc(segClassificationSet$isLV[valSet], probLV[valSet]), fill=T)
  
  # Distribution of probabilities
  plotSet <- group_by(data.frame(predictedProbability = cut(probLV, 10), #cut2(probLV, g=20), # equi-weight
                                 isLV = segClassificationSet$isLV,
                                 isVal = (seq(nrow(segClassificationSet)) %in% valSet)), predictedProbability) %>% 
    summarise(validation = sum(isLV & isVal)/sum(isVal),
              train = sum(isLV & !isVal)/sum(!isVal),
              count = n()) %>%
    gather(dataset, probability, -count, -predictedProbability)
  print(ggplot(plotSet, aes(x=predictedProbability, y=probability, fill=dataset)) + 
          geom_bar(stat="identity",position="dodge") + 
          ggtitle("Segment Prediction") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  
  # Keep data for analysis elsewhere
  write.csv(segClassificationSet, segmentPredictFile, row.names=F)
  
  # Apply on full dataset
  cat("Apply segment model to", nrow(allSegments), "segments", fill=T)
  allSegments$pLV <- predict(leftVentricleSegmentModel, 
                             data.matrix(createSegmentPredictSet(
                               left_join(allSegments, 
                                         imageList, by=c("Id","Slice","Time")))),
                             missing=NaN)
  
  # Keep data if we want to skip the segmentation predict phase
  allSegments <- select(allSegments, -UUID)
  write.csv(allSegments, imagePredictFile, row.names=F)
}

# Remove segments with pLV < threshold
allSegments <- filter(allSegments, pLV > pSegmentThreshold)
print(qplot(allSegments$pLV))

# Keep only the segments with max pLV for each image
allSegments[, isLV := segIndex == segIndex[which.max(pLV)], by=c("Id","Slice","Time")]
imageData <- createImagePredictSet(left_join(imageList, 
                                             filter(allSegments, isLV),
                                             by=c("Id", "Slice", "Time")))
cat("Total",nrow(imageData),"images, of which",nrow(filter(imageData,isLV)),"have a detected LV",fill=T)

# This should show that lower slice order have a higher probabilities for the left ventricle (better segmentation)
pLeftVentricle <- cut(imageData$pLV,breaks=seq(0,1,by=0.1))
print(ggplot(imageData, aes(x=pLeftVentricle, fill=factor(SliceOrder))) + geom_bar(position="dodge")+
        ggtitle("LV Probability vs Slice Order") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Random IDs for plotting
sampleIds <- sample(unique(imageData$Id), 10)

# Data cleansing: outlier detection and missing value imputation
zRange <- range(0: quantile(imageData$area, na.rm=T, p=0.90)) 
for (plotId in unique(imageData$Id)) {
  cat("Outlier detection and missing value imputation for",plotId,fill=T)
  data3D <- imageData[Id == plotId, sapply(imageData, is.numeric), with=F]
  
  # View raw data - with missing and outliers
  print(wireframe(area ~ Time*SliceLocation, data=data3D,
                  shade=T, col.regions = terrain.colors(100), 
                  main=paste("Raw Id=",plotId),
                  zlim=zRange))
  
  # Outliers removal by time and slice (TODO cant we do this at once somehow?)
  similarityData <- select(data3D, Slice, Time, area)
  nOutliers <- 0
  for (t in unique(similarityData$Time)) {
    d <- similarityData[Time==t]$area
    if (sum(!is.na(d)) > 3) {
      outlier.scores <- lofactor(na.omit(d), k=3)
      #plot(density(outlier.scores))
      outliers <- which(!is.na(d))[which(outlier.scores > 3.0)] # somewhat arbitrary
      nOutliers <- nOutliers + length(outliers)
      outlierSlices <- similarityData[Time==t]$Slice[outliers]
      data3D[Time==t & Slice %in% outlierSlices, area := NA]
    }
  }
  for (s in unique(similarityData$Slice)) {
    d <- similarityData[Slice==s]$area
    if (sum(!is.na(d)) > 3) {
      outlier.scores <- lofactor(na.omit(d), k=3)
      #plot(density(outlier.scores))
      outliers <- which(!is.na(d))[which(outlier.scores > 3.0)] # somewhat arbitrary
      nOutliers <- nOutliers + length(outliers)
      outlierTimes <- similarityData[Slice==s]$Time[outliers]
      data3D[Slice==t & Time %in% outlierTimes, area := NA]
    }
  }
  cat(nOutliers, "outliers for Id", plotId, fill=T)
  
  #   outlier.scores <- lofactor(scale(similarityData[which(complete.cases(similarityData))]), k=5)
  #   #plot(density(outlier.scores))
  #   outliers <- which(complete.cases(similarityData))[which(outlier.scores > 1.2)]
  #   cat("Outliers:",outliers,fill=T)
  #   data3D$area[outliers] <- NA
  
  print(wireframe(area ~ Time*SliceLocation, data=data3D,
                  shade=T, col.regions = terrain.colors(100), 
                  main=paste("Outliers Removed Id=",plotId), 
                  zlim=zRange))
  
  # Missing imputation
  
  # this effectively smoothed and imputes
  # should we help bagging and add some 'shifted' rows of area data?
  # only set the missing values, not all
  # kNN doesnt work when there are too many missings
  preds <- predict(preProcess(data3D, method="bagImpute"), newdata=data3D) # [,!sapply(data3D, anyNA),with=F]
  data3D$area <- ifelse(is.na(data3D$area), preds$area, data3D$area)

  print(wireframe(area ~ Time*SliceLocation, data=data3D,
                  shade=T, col.regions = terrain.colors(100), 
                  main=paste("Missings Imputed Id=",plotId), 
                  zlim=zRange))
  
  # copy imputed and smoothed data back to imageData
  for (i in seq(nrow(data3D))) {
    imageData[Id == data3D[i]$Id & Slice == data3D[i]$Slice & Time == data3D[i]$Time, area := data3D[i]$area ]
  }
}

#
# Plot some graphs
#

print(ggplot(filter(imageData, !is.na(pLV) & Id %in% sampleIds), aes(x=SliceLocation, y=area, colour=factor(Id)))+geom_line()+
        ggtitle("Segment area over Slice"))
print(ggplot(filter(imageData, !is.na(pLV) & Id %in% sampleIds), aes(x=Time, y=area, colour=factor(Id)))+geom_boxplot()+
        ggtitle("Segment area over Time"))

# Aggregate up to Time level
timeData <- group_by(imageData, Id, Time) %>%
  summarise(volume     = sum(SliceThickness*area, na.rm=T),
            #volumeEllipse = sum(SliceThickness*area.ellipse, na.rm=T),
            #volumeMax  = sum(SliceThickness*pi*radius.max^2, na.rm=T),
            #volumeMin  = sum(SliceThickness*pi*radius.min^2, na.rm=T),
            #volumeMean = sum(SliceThickness*pi*radius.mean^2, na.rm=T),
            isLV = any(isLV))
print(ggplot(filter(timeData, Id %in% sampleIds), aes(x=Time, y=volume, colour=factor(Id)))+geom_line()+geom_point()+
        ggtitle("Volume over Time"))

cat("Total",nrow(timeData),"images, of which",nrow(filter(timeData,isLV)),"have a detected LV",fill=T)
# poorly segmented slices: filter(timeData, is.na(isLV), with=T)
timeData <- filter(timeData, !is.na(isLV)) # keep only the ones with an LV

# Now, aggregate up to Id
caseList <- getIdList(playlist=imageList)
caseData <- left_join(caseList, group_by(timeData, Id) %>%
                        summarise(
                          max_volume = max(volume, na.rm=T),
                          min_volume = min(volume, na.rm=T),
                          sd_volume  = sd(volume, na.rm=T),
                          
                          #max_volumeEllipse = max(volumeEllipse, na.rm=T),
                          #min_volumeEllipse = min(volumeEllipse, na.rm=T),
                          
                          #max_volumeMax = max(volumeMax, na.rm=T),
                          #min_volumeMax = min(volumeMax, na.rm=T),
                          
                          #max_volumeMin = max(volumeMin, na.rm=T),
                          #min_volumeMin = min(volumeMin, na.rm=T),
                          
                          #max_volumeMean = max(volumeMean, na.rm=T),
                          #min_volumeMean = min(volumeMean, na.rm=T),
                          
                          isLV = any(isLV)), 
                      by="Id")
cat("Total",nrow(caseData),"cases, of which",nrow(filter(caseData,isLV)),"have a detected LV",fill=T)

# Train data
trainVolumes <- fread('data/train.csv') 
caseData <- left_join(caseData, trainVolumes, by="Id")

# All cases (dev/train/test) with engineered extra features
casePredictSet <- select(caseData, 
                         -Id, -Dataset, -ImgType, -isLV)
# casePredictSet <- mutate(casePredictSet, # extra vars
#                          maxVolume2 = max_volumeEllipse^2,
#                          minVolume2 = max_volumeEllipse^2,
#                          maxVolumeSlice = max_volumeEllipse*SliceCount,
#                          minVolumeSlice = max_volumeEllipse*SliceCount)

### Build model (caret)

casePredictSetTrain <- casePredictSet[!is.na(casePredictSet$Systole) & !is.na(casePredictSet$Diastole)] # set with known outcomes

# IDEA: run this repeatedly, get the distribution of the outcomes to get confidence interval
casePredictValidationRows <- sample.int(nrow(casePredictSetTrain), 
                                        validationPercentage*nrow(casePredictSetTrain))

# TEMP not using dev/val split
casePredictSetTrainVal <- casePredictSetTrain[casePredictValidationRows,]  # validation % of this, for reporting error
casePredictSetTrainDev <- casePredictSetTrain[-casePredictValidationRows,] # to be used for training the model

fitControl <- trainControl(method = "cv",number = 10)
systole_m <- train(Systole ~ ., data = select(casePredictSetTrainDev, -Diastole), 
                   method = "lm", trControl = fitControl, verbose=F)
diastole_m <- train(Diastole ~ ., data = select(casePredictSetTrainDev, -Systole),
                    method = "lm", trControl = fitControl, verbose=F)
print(plot(varImp(systole_m)))
print(plot(varImp(diastole_m)))

# Caret is nice however it does not give me the confidence intervals
# maybe to repeated with 1-out sampling to get a 700 x 700 matrix etc then find mean and so
# for now - just falling back to direct lm

systole_m <- lm(Systole ~ ., data = select(casePredictSetTrainDev, -Diastole))
diastole_m <- lm(Diastole ~ ., data = select(casePredictSetTrainDev, -Systole))

# Predict ALL cases 
resultData <- data.frame(caseData, 
                         Systole = predict(systole_m, 
                                           newdata = casePredictSet,
                                           interval="confidence", level=confidence),
                         Diastole = predict(diastole_m, 
                                            newdata = casePredictSet,
                                            interval="confidence", level=confidence))

print("Result data with predictions:")
print(summary(resultData))

# Fill in missing predictions
cat("Missing cases:",sum(is.na(resultData$Systole.fit) | is.na(resultData$Diastole.fit)),"out of",nrow(resultData),fill=T)
cat("   ",resultData$Id[which(is.na(resultData$Systole.fit) | is.na(resultData$Diastole.fit))],fill=T)
cat("   will be replaced by mean",fill=T)
predictions <- c("Systole.fit", "Systole.lwr", "Systole.upr", "Diastole.fit", "Diastole.lwr", "Diastole.upr")
for (pName in predictions) {
  resultData[[pName]] <- ifelse(is.na(resultData[[pName]]), mean(resultData[[pName]],na.rm=T), resultData[[pName]])
}


# Plot the results... [which(complete.cases(resultData)),]
plotData <- select(resultData, contains("stole"), Id) %>% 
  gather(Phase, Actual, -Id, -contains(".")) %>%
  mutate(Predicted = ifelse(Phase == "Diastole", Diastole.fit, Systole.fit))
print(ggplot(na.omit(plotData), 
             aes(x=Actual,y=Predicted,colour=Phase))+
        geom_point()+stat_smooth(method = "lm")+
        xlim(0, 600)+ylim(0, 600)+geom_abline(slope=1,linetype="dashed")+
        ggtitle("Actual vs Predicted..."))

print("Correlations:")
print(cor(resultData$Systole, resultData$Systole.fit, use="complete.obs"))
print(cor(resultData$Diastole, resultData$Diastole.fit, use="complete.obs"))

# Max Volume (mL) - fixed value used in submissions
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