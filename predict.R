library(tidyr)
library(ggplot2)
library(dplyr)
library(data.table)

imageData <- fread('imageData.csv')

imageData <- group_by(imageData, case, series)
# add some meta info: index of image in series, nr of differently segmented images
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

# scale measured to actual, combine in one plot
measured_med <- sapply(select(imageData, vol_min, vol_max), median)
actual_med <- sapply(select(filter(train, Id %in% imageData$Id), -Id), median)
scale <- actual_med[2] / measured_med[2]
print(scale)
# TODO combine both in one plot now

