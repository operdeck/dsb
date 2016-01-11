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

# First cut at 2-D normalization of the image data
allSegments <- mutate(allSegments,
                      areaMultiplier = pixelSpacing.x * pixelSpacing.y,
                      lengthMultiplier = sqrt(areaMultiplier),
                      
                      area = s.area*areaMultiplier,
                      
                      perimeter = s.perimeter*lengthMultiplier,
                      radius.mean = s.radius.mean*lengthMultiplier,
                      radius.min = s.radius.min*lengthMultiplier,
                      radius.max = s.radius.max*lengthMultiplier,
                      
                      # TODO read into "m.majoraxis"    "m.eccentricity"  "m.theta"
                      # TODO figure out if s.radius.sd should be scaled as well
                      
                      roundness = 4*pi*area/(perimeter^2)) %>%
  select(-pixelSpacing.x,
         -pixelSpacing.y,
         -areaMultiplier,
         -lengthMultiplier)

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
  for (i in seq(nrow(s))) {
    candidateSegments <- filter(allCandidateSegments, Id==s$Id[i], Slice==s$Slice[i], Time==s$Time[i]) # all segs in this image
    lv <- left_join(s[i], classifiedSegmentsOlder, by=c("Id","Slice","Time")) %>% filter(isLV) # classified LV in this image
    
    cat("Matching new classification to old segmentation for Id", s$Id[i], "Slice", s$Slice[i], "Image", s$Time[i], 
        "candidates:", nrow(candidateSegments), "classified:", nrow(lv), fill=T)
    
    candidateSegments$isLV <- NA
    candidateSegments <- candidateSegments[,names(segClassificationSet),with=F] # make sure cols have the same order
    
    if (nrow(lv) == 1) {
      # Set LV to the segment closest to the identified one. Note: similar code in classify.R
      distToLVSeg <- sqrt((candidateSegments$m.cx - lv$m.cx)^2 + (candidateSegments$m.cy - lv$m.cy)^2)
      segLV <- candidateSegments$segIndex[which.min(distToLVSeg)]
      candidateSegments$isLV <- (candidateSegments$segIndex == segLV & distToLVSeg < 5) # abs distance threshold like in classify.R
      segClassificationSet <- rbind(segClassificationSet, candidateSegments)
    } else {
      cat("WARN:",nrow(lv),"LV segments in one image", fill=T)
      print(s[i,])
    }
    
    segClassificationSet <- rbind(segClassificationSet, candidateSegments)
  }
}

segClassificationSet <- filter(segClassificationSet, !is.na(isLV)) # just to be sure

# quick plot to verify results
l <- unique(select(segClassificationSet, Id, Slice))
for (i in seq(nrow(l))) {
  slice <- filter(segClassificationSet, Id==l$Id[i], Slice==l$Slice[i])
  
  plotData <- mutate(filter(slice, isLV), 
                     area.radius.mean = pi*radius.mean^2,
                     area.radius.max = pi*radius.max^2,
                     area.radius.min = pi*radius.min^2) %>% 
    gather(metric, area, starts_with("area"))
  if (nrow(filter(slice, isLV)) > 1) {
    print(ggplot(plotData, aes(x=Time, y=area, colour=metric))+geom_line()+
            ggtitle(paste("Segment area over Time for ID",
                          unique(slice$Id),"Slice",unique(slice$Slice))))
  } else {
    cat("No image with identified LV at all for slice:", l$Id[i], l$Slice[i], fill=T)
  }
}

# Build up the data set for training and classification

# TODO maybe join in some slice info to get nr of images per slice etc and
# create some predictors out of that. For now, just dropping all meta info.

segClassificationSet <- select(segClassificationSet, 
                               -m.cx, -m.cy, -distToROI,
                               -Id, -Slice, -Time, -Dataset, -ImgType, -Offset, -segIndex)
allSegments <- select(allSegments, 
                      -m.cx, -m.cy, -distToROI,
                      -Id, -Slice, -Time, -Dataset, -ImgType, -Offset, -segIndex)

trainData <- segClassificationSet
valSet <- sample.int(nrow(trainData), 0.20*nrow(trainData))
cat("Building segment model with",length(names(trainData)),"predictors")
leftVentricleSegmentModel <- xgboost(data = as.matrix(select(trainData[-valSet], -UUID, -isLV)), 
                                     label = trainData[-valSet]$isLV, 
                                     max.depth = 2, eta = 0.1, nround = 100,
                                     objective = "binary:logistic", missing=NaN, verbose=0)
imp_matrix <- xgb.importance(feature_names = names(select(trainData, -UUID, -isLV)), model = leftVentricleSegmentModel)
print(xgb.plot.importance(importance_matrix = imp_matrix))

# Get an idea of the accuracy. Note, it seems very high always.
train_probLV <- predict(leftVentricleSegmentModel, as.matrix(select(trainData[valSet], -UUID, -isLV)))
cat("AUC for validation set:", auc(trainData[valSet]$isLV, train_probLV), fill=T)
plotSet <- group_by(data.frame(predictedProbability = cut(train_probLV, 10), 
                               isLV = trainData[valSet]$isLV), predictedProbability) %>% 
  summarise(percentageLV = sum(isLV)/n())
print(ggplot(plotSet, aes(x=predictedProbability, y=percentageLV)) + geom_bar(stat="identity") + 
  ggtitle("LV predicted vs actual") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Apply on full dataset
# p <- predict(leftVentricleSegmentModel, as.matrix(allSegments))

stop()

# Read all segment info from all datasets
for (dataset in c('train','validate','test')) { 
  predictSet <- NULL
  fname <- paste("segments-",dataset,".csv",sep="")
  if (file.exists(fname)) {
    segInfo <- fread(fname)
    cat("Id's:",unique(segInfo$Id),fill=T)
    
    # Select only the middle segment for each Id (TODO: consider broader selection - see 'classify.R')
    slicesToConsider <- group_by(segInfo, Id) %>% 
      summarise(midSlice = sort(unique(Slice))[ceiling(length(unique(Slice))/2)]) 
    segInfo <- left_join(segInfo, slicesToConsider) %>% filter(Slice == midSlice)
    
    newData <- segInfo[, -which(sapply(segInfo, class) == "character"), with=F]
    newData <- newData[, -which(names(newData) %in% dropCols),with=F]
    setkey(newData, Id, Slice, Time, segIndex)
    setkeyv(segInfo, key(newData))
    segInfo$pLeftVentricle <- predict(leftVentricleSegmentModel, as.matrix(newData))
    segInfo <- left_join(segInfo, group_by(segInfo, Id, Slice, Time) %>% 
                           summarise(lvSeg = segIndex[which.max(pLeftVentricle)]) )
    # show lv segments:
    # print(unique(select(segInfo, Id, Slice, Time, lvSeg)))
    
    cat("AUC for train data set:",
        auc(segData$isLV, predict(leftVentricleSegmentModel, as.matrix(trainData))),
        fill=T)
    
    lvSegsOnly <- group_by(segInfo, Id, Slice, Time) %>% 
      summarise(UUID = UUID[which.max(pLeftVentricle)]) %>% left_join(segInfo)
    
    # TODO limit to only first 10 id's (or so)
    print(ggplot(filter(lvSegsOnly), aes(x=Time, y=sliceVolume, colour=factor(Id)))+geom_line()+
            ggtitle(paste("Volume over Time")))
    
    # Roll up over Time and Slice dimensions
    
    # TODO: roll up time dim. For all selected attributes, add sd, mean, min and max or p10 and p90
    timeRollUpConstCols <- c("Offset","sliceLocation","sliceThickness")
    timeRollUpCols <- c("m.majoraxis","m.eccentricity","m.theta",
                        "s.area","s.perimeter",
                        "s.radius.mean","s.radius.sd","s.radius.min","s.radius.max",
                        "roundness","dist","sliceArea","sliceVolume")
    allSlices <- unique(select(lvSegsOnly, Id, Slice))
    for (s in 1:nrow(allSlices)) {
      rollUpData <- filter(lvSegsOnly, Id == allSlices[s]$Id, Slice == allSlices[s]$Slice)
      rollUpData <- rollUpData[, which(names(rollUpData) %in% timeRollUpCols),with=F]
      
      data_mean <- sapply(rollUpData, mean)
      names(data_mean) <- paste(names(rollUpData), "mean", sep="_")
      
      data_sd <- sapply(rollUpData, sd)
      names(data_sd) <- paste(names(rollUpData), "sd", sep="_")
      
      data_quantiles <- sapply(rollUpData, quantile, probs=c(0.1,0.5,0.9))
      data_p10 <- data_quantiles[1,]
      names(data_p10) <- paste(names(rollUpData), "p10", sep="_")
      data_p50 <- data_quantiles[2,]
      names(data_p50) <- paste(names(rollUpData), "p50", sep="_")
      data_p90 <- data_quantiles[3,]
      names(data_p90) <- paste(names(rollUpData), "p90", sep="_")
      
      # TODO add 'first' for values relative constant to Time (Slice factors)
      data_aggregated <- c(Id = allSlices[s]$Id, Slice = allSlices[s]$Slice, 
                           data_mean, data_sd, data_p10, data_p50, data_p90)
      
      if (is.null(predictSet)) {
        predictSet <- as.data.frame(t(data_aggregated))
      } else {
        predictSet <- rbind(predictSet,t(data_aggregated))
      }
    }
    
    if (dataset == 'train') {
      train_predictSet <- predictSet
    } else if (dataset == 'validate') {
      validate_predictSet <- predictSet
    } else if (dataset == 'test') {
      test_predictSet <- predictSet
    } else {
      stop("Unknown dataset...")
    }
    
    # TODO roll up by Slice (later...)
    
    
  }
}

train <- fread('data/train.csv') 
print(ggplot(data=filter(train, Id %in% train_predictSet$Id) %>% gather(Phase, Volume, -Id), 
             aes(x=factor(Id, levels=sort(unique(Id))), y=Volume))+geom_line(size=5)+
        ggtitle("Actual volume in train set"))

# Prepare for predictions

# Select only the N most correlated fields with N = nrow / 2 (only because train set is limited sometimes)
predictData <- left_join(train_predictSet, train, "Id")
predictDataPredictorsOnly <- select(predictData, -contains("stole"), -Id, -Slice)
nPreds <- min(ceiling(sqrt(nrow(predictData))), ncol(predictDataPredictorsOnly))
predNames <- names(-sort(-sapply(predictDataPredictorsOnly, 
                                 function(x) {return (mean(abs(cor(predictData$Diastole, x)),
                                                           abs(cor(predictData$Systole, x))))}))[1:nPreds])
predictDataPredictorsOnly <- predictDataPredictorsOnly[, which(names(predictDataPredictorsOnly) %in% predNames)]

# Develop model (TODO: on -validation set)
# TODO of course we will create more powerful models when there is more data
systole_m <- lm(predictData$Systole ~ ., data = predictDataPredictorsOnly)
diastole_m <- lm(predictData$Diastole ~ ., data = predictDataPredictorsOnly)

resultData <- data.frame(predictData, 
                         Systole_pred = predict(systole_m, 
                                                newdata = predictDataPredictorsOnly),
                         Diastole_pred = predict(diastole_m, 
                                                 newdata = predictDataPredictorsOnly))

# Plot actual vs predicted
plotData <- select(resultData, contains("Diastole"), Id) %>% gather(Source, Volume, -Id)
plotData$Source <- factor(plotData$Source, levels=sort(unique(as.character(plotData$Source))))
print(ggplot(plotData, 
             aes(x=Id, y=Volume, fill=Source)) + geom_bar(stat="identity",position="dodge")+
        ggtitle("Diastole predictions"))

plotData <- select(resultData, contains("Systole"), Id) %>% gather(Source, Volume, -Id)
plotData$Source <- factor(plotData$Source, levels=sort(unique(as.character(plotData$Source))))
print(ggplot(plotData, 
             aes(x=Id, y=Volume, fill=Source)) + geom_bar(stat="identity",position="dodge")+
        ggtitle("Systole predictions"))

print("Correlations:")
print(cor(resultData$Systole, resultData$Systole_pred))
print(cor(resultData$Diastole, resultData$Diastole_pred))

#TODO: report accuracy (RMSD?) on validation set

#TODO: translate predicted values to distributions of the volumes

plotData <- select(resultData, contains("stole"), Id) %>% 
  gather(Phase, Actual, -Id, -Diastole_pred, -Systole_pred) %>%
  mutate(Predicted = ifelse(Phase == "Diastole", Diastole_pred, Systole_pred))
print(ggplot(plotData, aes(x=Actual,y=Predicted,colour=Phase))+geom_point()+stat_smooth(method = "lm")+ggtitle("Actual vs Predicted..."))

# predict on validate and test sets

validate_resultData <- data.frame(select(validate_predictSet, Id), 
                                  Systole_pred = predict(systole_m, newdata = validate_predictSet),
                                  Diastole_pred = predict(diastole_m, newdata = validate_predictSet))

print(head(validate_resultData))

# now, translate this to probabilities by volume
# note, predictions can be NA, negative or otherwise out of bounds,
# also make sure Systole << Diastole - otherwise just fall back to default



