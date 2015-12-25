library(tidyr)
library(ggplot2)
library(dplyr)
library(data.table)

imageData <- fread('imageData.csv')
imageData$Id <- factor(imageData$Id,levels=sort(unique(imageData$Id)))

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

plotData <- select(resultData, contains("stole"), Id) %>% 
  gather(Phase, Actual, -Id, -Diastole_pred, -Systole_pred) %>%
  mutate(Predicted = ifelse(Phase == "Diastole", Diastole_pred, Systole_pred))
print(ggplot(plotData, aes(x=Actual,y=Predicted,colour=Phase))+geom_point()+stat_smooth(method = "lm")+ggtitle("Actual vs Predicted..."))


