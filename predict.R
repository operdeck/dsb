library(tidyr)
library(ggplot2)
library(dplyr)
library(data.table)

imageData <- fread('imageData.csv')

imageData <- group_by(imageData, case, series)
# add some meta info: index of image in series, nr of differently segmented images
# TODO: do this in extract.R - should all not be necessary anymore
for (n in 1:nrow(imageData)) {
  if (n == 1 || (imageData$series[n] != imageData$series[n-1])) {
    imageData$in_series_rank[n] <- 1
  } else {
    imageData$in_series_rank[n] <- 1+imageData$in_series_rank[n-1]
  }
  if (n == 1 || (imageData$case[n] != imageData$case[n-1])) {
    imageData$series_idx[n] <- 1
  } else {
    imageData$series_idx[n] <- 1+imageData$series_idx[n-1]
  }
}
imageData <- ungroup(imageData) %>% 
  mutate(Id=case) %>%
  select(-case)

imageData$Id <- factor(imageData$Id,levels=sort(unique(imageData$Id)))
print(ggplot(data=imageData, aes(x=series_idx, y=vol_max, colour=Id))+geom_point()+geom_line()+
  ggtitle("Volume vs image series per case"))

sumImageData <- select(imageData, vol_min, vol_max, Id) %>% gather(Phase, Volume, -Id)
print(ggplot(data=sumImageData, aes(y=Volume, x=Id, fill=Phase))+geom_boxplot()+
  ggtitle("Measured volume per case"))

train <- fread('data/train.csv') 
train$Id <- factor(train$Id,levels=sort(unique(train$Id)))
print(ggplot(data=filter(train, Id %in% imageData$Id) %>% gather(Phase, Volume, -Id), 
             aes(x=Id, y=Volume))+geom_line(size=5)+
        ggtitle("Actual volume in train set"))

# Prepare for predictions

predictData <- group_by(imageData, Id) %>%
  summarise(vol_min_median = median(vol_min),
            vol_max_median = median(vol_max)) # etc for all other possibly relevant predictors
predictData <- left_join(predictData, train, "Id")

#TODO: develop model on -validation set
systole_m <- lm(Systole ~ ., data = select(predictData, -Diastole, -Id))
diastole_m <- lm(Diastole ~ ., data = select(predictData, -Systole, -Id))

resultData <- data.frame(predictData, 
                         Systole_pred = predict(systole_m, newdata = select(predictData, -Diastole, -Id)),
                         Diastole_pred = predict(diastole_m, newdata = select(predictData, -Systole, -Id)))

plotData <- select(resultData, contains("stole"), Id) %>% gather(Source, Volume, -Id)
plotData$Source <- factor(plotData$Source, levels=sort(unique(as.character(plotData$Source))))
print(ggplot(plotData, 
       aes(x=Id, y=Volume, fill=Source)) + geom_bar(stat="identity",position="dodge"))

print("Correlations:")
print(cor(resultData$Systole, resultData$Systole_pred))
print(cor(resultData$Diastole, resultData$Diastole_pred))
#TODO: report accuracy (RMSD?) on validation set

#TODO: translate predicted values to distributions of the volumes

