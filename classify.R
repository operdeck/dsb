# Manual identification of LV segments with the purpose of creating a prediction dataset
# to do this automatically on the full dataset. 
# Creates "segments-predict.csv" with same data for only the (correctly) identified segments of
# each file / slice / Id.

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
  for (t in unique(slice$Time)) {
    f <- getImageFile(slice[which(slice$Time == t)[1],])
    dicom <- readDICOMFile(f)
    img <- Image(normalize(dicom$img))
    if (dim(img)[1] > dim(img)[2]) {
      img <- rotate(img,-90)
    }
    seg <- which(slice$Time == t & slice$isLV)
    imgz[[f]] <- drawCircle(toRGB(img), 
                            x=slice$m.cx[seg], y=slice$m.cy[seg], radius=slice$s.radius.mean[seg], col="red")
  }
  EBImage::display(EBImage::combine(imgz),all=T,method="raster")
  text(500,40,paste("ID",unique(slice$Id),"Slice",unique(slice$Slice)),col="yellow",pos=4)
}


# Read all segment info from the train set only
# NOTE currently some slice meta info seems to be part of the segmentation data, but
# this shouldnt be the case. Regardless, drop it and join in from the official meta data.
allImages <- getImageList()
allSegmentationInfo <- 
  left_join(select(fread(getSegmentFile("train")), 
                   -FileName, -Offset, -SliceCount, -SliceIndex, -SliceOrder, 
                   -PixelSpacing.x, -PixelSpacing.y, -SliceLocation, -SliceThickness), 
            allImages,
            by=c("Dataset","Id","ImgType","Slice","Time"))

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
promptSlices <- unique(select(allSegmentationInfo, Id, Slice)) 
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
  slice <- filter(allSegmentationInfo, 
                  Id == promptSlices$Id[sliceIndex], 
                  Slice == promptSlices$Slice[sliceIndex])
  if (!all(file.exists(getImageFile(slice)))) {
    next
  }
  
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
    showSingleImage(segmentsFirstImage) 
    indexOfLVSegment <- readInteger(paste("Identify segment of left ventricle in image",
                                          firstImageIndex,
                                          "(0=none, -1=can't see): "))
    
    if (indexOfLVSegment < 1) {
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
      distThreshold <- 2 # 3.75 / pixelLengthScale # 5 for 'normal' images
      
      slice[Time %in% restOfUnprocessedImageIndices, isLV := NA]
      slice[Time %in% restOfUnprocessedImageIndices, distToLVSeg := sqrt((m.cx - segmentLV$m.cx)^2 + (m.cy - segmentLV$m.cy)^2)]
      slice[Time %in% restOfUnprocessedImageIndices, segLV := segIndex[which.min(distToLVSeg)], by=Time]
      slice[Time %in% restOfUnprocessedImageIndices, isLV := (segIndex == segLV) & (distToLVSeg < distThreshold)]
      slice[Time %in% restOfUnprocessedImageIndices, anyLV := any(isLV, na.rm=T), by=Time]
      
      segmentsWithTentativeLV <- filter(slice, Time %in% restOfUnprocessedImageIndices, anyLV)
      imagesIndicesWithTentativeLV <- unique(segmentsWithTentativeLV$Time)
      cat(length(imagesIndicesWithTentativeLV),"images with tentative Left Ventricle assignment:",imagesIndicesWithTentativeLV,fill=T)
      
      # Show the (tentative) results for review and prompt
      if (length(imagesIndicesWithTentativeLV) > 0) {
        plotSlice(rename(filter(segmentsWithTentativeLV, isLV), 
                         radius.mean = s.radius.mean,
                         radius.max = s.radius.max,
                         radius.min = s.radius.min))
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
    segmentedByID <- left_join(getIdList(playlist=allImages), 
                               unique(select(segPredictSet, Id, isLV)), by=c("Id")) %>% 
      group_by(Dataset, Id) %>% summarise(hasSegs = any(isLV))
    segmentedBySlice <- left_join(getSliceList(playlist=allImages), 
                                  unique(select(segPredictSet, Id, Slice, isLV)), by=c("Id", "Slice")) %>% 
      group_by(Dataset, Id, Slice) %>% summarise(hasSegs = any(isLV))
    segmentedByImage <- left_join(getImageList(playlist=allImages), 
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
  
  