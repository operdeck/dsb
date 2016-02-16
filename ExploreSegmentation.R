source("Util.R")

# Exercise the feature engineering in 'createSegmentPredictSet' and validate the XGB model

classifiedSegments <- fread('segments-predict.csv') 
imageList <- getImageList()
segmentData <- fread('segments-train.csv')

d <- select(classifiedSegments, UUID, isLV) %>%
  left_join(segmentData, "UUID") %>%
  filter(!is.na(Id)) %>%
  left_join(imageList, by=c("Id","Slice","Time"))

train <- createSegmentPredictSet(d)

# NB below should be in above too - now no chars left
for (f in setdiff(names(train),"isLV")) {
  if (class(train[[f]])=="character") {
    levels <- unique(c(train[[f]], test[[f]]))
    train[[f]] <- as.integer(factor(train[[f]], levels=levels))
  }
}

train[is.na(train)] <- 0
train$isLV <- factor(train$isLV, c(T,F), c('True', 'False'))
nearZV <- nearZeroVar(train, saveMetrics = T)
print(nearZV)

inTrain <- sample(nrow(train), 0.8*nrow(train))
x.train <- train[inTrain, ]
x.test <- train[-inTrain, ]

xgb.grid <- expand.grid(nrounds=1:6 * 10,
                        max_depth=c(4, 6, 8, 10),
                        eta=1:3 * 0.1,
                        gamma=1,
                        colsample_bytree=0.9,
                        min_child_weight=5)

ctrl <- trainControl(
  method = "repeatedcv", # 10-fold repeated CV
  number = 5,
  repeats = 2,
  classProbs = T,
  verboseIter = T)

fit <- train(isLV ~ ., data=train[inTrain,], method='xgbTree', metric='AUC',
             trControl=ctrl, tuneGrid=xgb.grid)
print(fit)
print(plot(fit))
print(ggplot(varImp(fit)))


