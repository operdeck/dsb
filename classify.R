# Manual identification of LV segments with the purpose of creating a prediction dataset
# to do this automatically on the full dataset. 
# Creates "segments-predict.csv" with same data for only the (correctly) identified segments of
# each file / slice / Id.

# TODO if segment picture is available, show that next to the base one,
# overlay with original

source("util.R")

gc()
segPredictFile <- "segments-predict.csv"

# Id's of images that need re-classification because of manual error. All slices will be discarded.
redoClassificationIdList <- c() # list the ID's here

# TODO maybe matching UUID's (re-segmented) should take priority 
# TODO auto detect classified images with no LV segments

# Problematic images - need re-segmentation
# Too much light around edges:
#  train 6, slice 10
#  train 8, slice 58
#  train 8, slice 59
# Very dark, some images OK but some lack segments:
#  train 3, slice 46&47
# Near miss
#  train 7, slice 49 : almost all images OK but segment 13 (or so) is incorrect


showSegmentLabels <- function(imageDS) {
  for (i in 1:nrow(imageDS)) {
    text(x = imageDS$m.cx[i], y = imageDS$m.cy[i], 
         label = imageDS$segIndex[i], col = "red")
  }
}

showSingleImage <- function(segments) {
  f <- getImageFile(segments)
  if (!file.exists(f)) {
    print(f)
    stop("File doesnt exist on this system")
  }
  dicom <- readDICOMFile(f)
  img <- Image(normalize(dicom$img))
  if (dim(img)[1] > dim(img)[2]) {
    img <- rotate(img,-90)
  }
  for (j in 1:nrow(segments)) {
    radii <- c(segments$s.radius.mean[j], segments$s.radius.min[j], segments$s.radius.max[j])
    for (k in seq(length(radii))) {
      if (radii[k] >= 1) {
        img <- drawCircle(toRGB(img), 
                          x=segments$m.cx[j], y=segments$m.cy[j], radii[k], col=ifelse(k==1,"red","orange"))
      }
    }
  }
  EBImage::display(img,all=T,method="raster")
  text(10,20,getImageFile(segments),col="yellow",pos=4)
  showSegmentLabels(segments)
}

# Show all images of this slice
showAllSliceImages <- function(slice) {
  imgz <- list()
  #print(select(slice, Time, segIndex, starts_with("m."), starts_with("s.")))
  for (t in unique(slice$Time)) {
    f <- getImageFile(slice[which(slice$Time == t)[1],])
    dicom <- readDICOMFile(f)
    img <- Image(normalize(dicom$img))
    if (dim(img)[1] > dim(img)[2]) {
      img <- rotate(img,-90)
    }
    seg <- which(slice$Time == t & slice$isLV)
    #print(slice[seg,])
    
    radii <- c(slice$s.radius.mean[seg], slice$s.radius.min[seg], slice$s.radius.max[seg])
    for (k in seq(length(radii))) {
      if (radii[k] >= 1) {
        img <- drawCircle(toRGB(img), 
                          x=slice$m.cx[seg], y=slice$m.cy[seg], radii[k], col=ifelse(k==1,"red","orange"))
      }
    }
    imgz[[f]] <- img
  }
  EBImage::display(EBImage::combine(imgz),all=T,method="raster")
  text(500,40,paste("ID",unique(slice$Id),"Slice",unique(slice$Slice)),col="yellow",pos=4)
}

# Build (simple) segment prediction model and apply to test set
createSegmentModelAndApply <- function(train, test)
{
  cat("LV model built out of", sum(complete.cases(train)), 
      "complete cases from total",nrow(train),fill=T)
  train <- train[complete.cases(train),]
  preds <- rep(1.0/nrow(test), nrow(test))
  if(is.null(train) || nrow(train) < 10000) {
    return (preds)
  }

  segments <- test$segIndex
  
  train <- left_join(filter(train, !is.na(isLV)), imageList, by=c("Id", "Slice", "Time"))
  test <- test[, names(test) %in% setdiff(names(train), "isLV"), with=F]

  train <- createSegmentPredictSet(train)
  test <- createSegmentPredictSet(test)
  
  leftVentricleSegmentModel <- xgboost(data = data.matrix(select(train, -isLV)),
                                       missing=NaN,
                                       label = train$isLV, 
                                       max.depth = 4, eta = 0.1, nround = 100,
                                       objective = "binary:logistic", verbose=0)

  imp_matrix <- xgb.importance(feature_names = names(select(train, -isLV)), model = leftVentricleSegmentModel)
  print(xgb.plot.importance(importance_matrix = imp_matrix))
  
  probLV <- predict(leftVentricleSegmentModel, data.matrix(test), missing=NaN) #, missing=NaN)
  names(probLV) <- segments
  
  return (probLV)
}


# Read all segment info from the train set only and join with meta data
imageList <- getImageList()
allSegments <- 
  left_join(fread(getSegmentFile("train")), 
            imageList,
            by=c("Id","Slice","Time"))

# Read previous segment classification
segPredictSet <- NULL
if (file.exists(segPredictFile)) {
  segPredictSet <- fread(segPredictFile)
  
  # Only use a small set of fairly static attributes that are enough to identify a segment
  segPredictSet <- select(segPredictSet, 
                          Id, Slice, Time,   # uniquely identifies , image via join to playlist
                          UUID,              # unique segment identifier
                          starts_with("m."), # standard moments & shape attributes 
                          starts_with("s."),
                          isLV)
} else {
  print("No classification data yet")
}

# Select slices to prompt for
promptSlices <- unique(select(allSegments, Id, Slice)) 
promptSlices$ReDo <- promptSlices$Id %in% redoClassificationIdList
promptSlices$Random <- runif(nrow(promptSlices))
if (!is.null(segPredictSet)) {
  promptSlices <- left_join(promptSlices, 
                            mutate(unique(select(segPredictSet, Id, Slice)), classifiedAlready=T)) %>%
    mutate(classifiedAlready = !is.na(classifiedAlready)) %>%
    arrange(!ReDo, classifiedAlready, Random) # previously: Id, Slice instead of Random
}

if (nrow(promptSlices) < 1) {
  stop("Nothing to do, no slices to prompt for.")
}

for (sliceIndex in seq(nrow(promptSlices))) {
  # all segments of all images of one slice
  slice <- filter(allSegments, 
                  Id == promptSlices$Id[sliceIndex], 
                  Slice == promptSlices$Slice[sliceIndex])
  if (!all(file.exists(getImageFile(slice)))) {
    next
  }
  # Create simple segment model and predict pLV on whole slice here
  slice$pLV <- createSegmentModelAndApply(segPredictSet, slice)

  # list of images in this slice that need processing
  unProcessedImageIndices <- unique(slice$Time)
  while (length(unProcessedImageIndices) > 0) {
    cat("Id=",promptSlices$Id[sliceIndex],"Slice=",promptSlices$Slice[sliceIndex],
        "#Images:",length(unique(slice$Time)),
        "#Segments:",nrow(slice),fill=T)
    cat(length(unProcessedImageIndices),"images:",unProcessedImageIndices,fill=T)
    firstImageIndex <- unProcessedImageIndices[1] # all segments of first image
    if (length(unProcessedImageIndices) > 1) { # all segments of rest of images
      restOfUnprocessedImageIndices <- unProcessedImageIndices[2:length(unProcessedImageIndices)]
    } else {
      restOfUnprocessedImageIndices <- c()
    }
    
    # show firstImage
    segmentsFirstImage <- filter(slice, Time == firstImageIndex)
    segmentProbs <- round(segmentsFirstImage$pLV, 3)
    names(segmentProbs) <- segmentsFirstImage$segIndex
    print("Segments with prediction of being the Left Ventricle:")
    print(head(sort(segmentProbs, decreasing=T), 5))
    showSingleImage(segmentsFirstImage) 
    indexOfLVSegment <- readInteger(paste("Identify Left Ventricle in image",
                                          firstImageIndex,
                                          "(0=none, -1=skip image, -2=skip slice): "))
    
    if (indexOfLVSegment < 1) {
      if (indexOfLVSegment == -2) {
        break 
      }
      # set isLV for all segments of first image only
      if (indexOfLVSegment == 0) {
        slice[Time == firstImageIndex, isLV := FALSE] # no LV segment
      } else {
        slice[Time == firstImageIndex, isLV := NA]    # don't know / can't see
      }
      unProcessedImageIndices <- restOfUnprocessedImageIndices
    } else {
      # set isLV for all segments of first image (T/F)
      slice[Time == firstImageIndex, isLV := (segIndex == indexOfLVSegment)]
      
      # tentatively set isLV for all segments of the rest of the images (T/F/NA)
      segmentLV <- segmentsFirstImage[segmentsFirstImage$segIndex == indexOfLVSegment,]
      pixelAreaScale <- segmentLV$PixelSpacing.x * segmentLV$PixelSpacing.y
      pixelLengthScale <- sqrt(pixelAreaScale)
      distThreshold <- 3.75 / pixelLengthScale # 5 for 'normal' images
      
      slice[Time %in% restOfUnprocessedImageIndices, isLV := NA]
      slice[Time %in% restOfUnprocessedImageIndices, distToLVSeg := sqrt((m.cx - segmentLV$m.cx)^2 + (m.cy - segmentLV$m.cy)^2)]
      slice[Time %in% restOfUnprocessedImageIndices, segLV := segIndex[which.min(distToLVSeg)], by=Time]
      slice[Time %in% restOfUnprocessedImageIndices, isLV := (segIndex == segLV) & (distToLVSeg < distThreshold)]
      slice[Time %in% restOfUnprocessedImageIndices, anyLV := any(isLV, na.rm=T), by=Time]
      
      segmentsWithTentativeLV <- filter(slice, Time %in% restOfUnprocessedImageIndices, anyLV)
      imagesIndicesWithTentativeLV <- unique(segmentsWithTentativeLV$Time)
      cat(length(imagesIndicesWithTentativeLV),"images with tentative Left Ventricle assignment:",imagesIndicesWithTentativeLV,fill=T)
      
      # Show the (tentative) results with segment probabilities for review and prompt
      if (length(imagesIndicesWithTentativeLV) > 0) {
        plotSlice(rename(filter(segmentsWithTentativeLV, isLV), 
                         radius.mean = s.radius.mean,
                         radius.max = s.radius.max,
                         radius.min = s.radius.min))
        print("Prediction of being the Left Ventricle in plot:")
        segmentProbs <- round(filter(segmentsWithTentativeLV, isLV)$pLV, 3)
        ncol <- ceiling(sqrt(length(segmentProbs)))
        nrow <- ceiling(length(segmentProbs)/ncol)
        print(matrix(segmentProbs[1:(ncol*nrow)], ncol=ncol, nrow=nrow, byrow=T))
        showAllSliceImages(segmentsWithTentativeLV)
        identifiedCorrectly <- readline("Are all segments identified correctly (y/n): ")
        if (identifiedCorrectly == "y") {
          # keep isLV for the rest of the images for which isLV has been set
          unProcessedImageIndices <- setdiff(restOfUnprocessedImageIndices, imagesIndicesWithTentativeLV)
        } else {
          # ignore the tentative classification for the rest of the images
          unProcessedImageIndices <- restOfUnprocessedImageIndices
        }
      } else {
        # ignore the tentative classification for the rest of the images
        unProcessedImageIndices <- restOfUnprocessedImageIndices
      }      
    }
    View(slice)
    print("Left Ventricle segment summary:")
    print(table(slice$isLV, useNA="always"))
  }
  
  if (indexOfLVSegment == -2) {
    print("User triggered skip of this slice")
    next
  }
  # save data
  newData <- select(slice,
                    Id, Slice, Time,   # uniquely identifies image via join to playlist
                    UUID,              # unique segment identifier
                    starts_with("m."), # standard moments & shape attributes 
                    starts_with("s."),
                    isLV)
  
  if(nrow(newData) > 0) {
    if (is.null(segPredictSet)) {
      segPredictSet <- newData
    } else {
      # Check both are compatible before saving
      removedSet <- setdiff(names(segPredictSet), names(newData))
      addedSet <- setdiff(names(newData), names(segPredictSet))
      diffSet <- paste(c(paste("-",removedSet), paste("+",addedSet)),collapse=", ")
      if (length(removedSet) + length(addedSet) > 0) {
        print("Existing:")
        print(names(segPredictSet))
        print("New:")
        print(names(newData))
        print(paste(diffSet, collapse=","))
        stop("Datasets do not match up. Please consider removing or fixing segment classification.")
      }
      nDrop <- nrow(filter(segPredictSet, 
                           Slice %in% newData$Slice,
                           Id %in% newData$Id))
      if (nDrop > 0) {
        cat("...dropping",nDrop,"rows from previous results",fill=T)
      }
      segPredictSet <- rbind(filter(segPredictSet, 
                                    !(Slice %in% newData$Slice & 
                                        Id %in% newData$Id)), 
                             newData)
    }
    
    print("Progress:")
    segmentedByID <- left_join(getIdList(playlist=imageList), 
                               unique(select(segPredictSet, Id, isLV)), by=c("Id")) %>% 
      group_by(Dataset, Id) %>% summarise(hasSegs = any(isLV))
    segmentedBySlice <- left_join(getSliceList(playlist=imageList), 
                                  unique(select(segPredictSet, Id, Slice, isLV)), by=c("Id", "Slice")) %>% 
      group_by(Dataset, Id, Slice) %>% summarise(hasSegs = any(isLV))
    segmentedByImage <- left_join(getImageList(playlist=imageList), 
                                  unique(select(segPredictSet, Id, Slice, Time, isLV)), by=c("Id", "Slice", "Time")) %>% 
      group_by(Dataset, Id, Slice, Time) %>% summarise(hasSegs = any(isLV))
    
    classificationProgress <- data.frame(segmented = c(sum(segmentedByID$hasSegs,na.rm=T),
                                                       sum(segmentedBySlice$hasSegs,na.rm=T),
                                                       sum(segmentedByImage$hasSegs,na.rm=T)),
                                         total = c(nrow(segmentedByID),
                                                   nrow(segmentedBySlice),
                                                   nrow(segmentedByImage)))
    classificationProgress <- mutate(classificationProgress,
                                     percentage = paste(round(100*segmented/total,2),"%",sep=""))
    rownames(classificationProgress) <- c('By ID','By Slice','By Image')
    print(classificationProgress)
    
    write.csv(segPredictSet, segPredictFile, row.names=F)
  }
}

print("Done!")

