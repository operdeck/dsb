# Make predictions of Systole and Diastole for all datasets.
# Uses "segments-classified.csv" to create a model to help find LV segments, then filters and
# aggregates the segment data to have one set of meta-info per ID. Combines this with the 
# train set to create a model for the Systole and Diastole volumes.

source("util.R")

library(pROC)

# Identify the LV segments in the whole dataset by creating a model
# from the (manually) identified ones
dropCols <- c('m.cx',  'm.cy')  # drop Id, segIndex
segData <- fread('segments-predict.csv', drop=dropCols)
trainData <- select(segData, -isLV)
trainData <- trainData[, -which(sapply(trainData, class) == "character"), with=F]
setkey(trainData, Id, Slice, Time, segIndex)
setkeyv(segData, key(trainData))
leftVentricleSegmentModel <- xgboost(data = as.matrix(trainData), 
                                     label = segData$isLV, 
                                     max.depth = 2, eta = 0.1, nround = 100,
                                     objective = "binary:logistic", missing=NaN, verbose=0)
imp_matrix <- xgb.importance(feature_names = names(trainData), model = leftVentricleSegmentModel)
print(xgb.plot.importance(importance_matrix = imp_matrix))

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



