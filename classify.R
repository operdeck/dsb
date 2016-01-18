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

showSingleImage <- function(ds) {
  f <- getImageFile(ds[1,])
  if (!file.exists(f)) {
    print(f)
    stop("File doesnt exist on this system")
  }
  dicom <- readDICOMFile(f)
  img <- Image(normalize(dicom$img))
  if (dim(img)[1] > dim(img)[2]) {
    img <- rotate(img,-90)
  }
  for (j in 1:nrow(ds)) {
    radii <- c(ds$s.radius.mean[j], ds$s.radius.min[j], ds$s.radius.max[j])
    for (k in seq(length(radii))) {
      if (radii[k] >= 1) {
        img <- drawCircle(toRGB(img), 
                          x=ds$m.cx[j], y=ds$m.cy[j], radii[k], col=ifelse(k==1,"red","orange"))
      }
    }
  }
  EBImage::display(img,all=T,method="raster")
  text(10,20,getImageFile(ds),col="yellow",pos=4)
  showSegmentLabels(ds)
}

# Show all images of this slice
showAllSliceImages <- function(ds) {
  imgz <- list()
  for (t in unique(ds$Time)) {
    f <- getImageFile(ds[which(ds$Time == t)[1],])
    dicom <- readDICOMFile(f)
    img <- Image(normalize(dicom$img))
    if (dim(img)[1] > dim(img)[2]) {
      img <- rotate(img,-90)
    }
    seg <- which(ds$Time == t & ds$isLV)
    imgz[[f]] <- drawCircle(toRGB(img), 
                            x=ds$m.cx[seg], y=ds$m.cy[seg], radius=ds$s.radius.mean[seg], col="red")
  }
  EBImage::display(EBImage::combine(imgz),all=T,method="raster")
  text(500,40,paste("ID",unique(ds$Id),"Slice",unique(ds$Slice)),col="yellow",pos=4)
}


# Read all segment info from the train set only
allSegmentationInfo <- fread(getSegmentFile("train"))

# TODO - join with image info here for the meta attributes and file name etc.

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
}

# Read all slices with meta info
sliceInfo <- getSliceList() # TODO should this be ImageList?

# TODO: don't
# Select only the middle slices (usually 2 per Id)
# TODO consider a wider selection
allSegmentationInfo <- left_join(allSegmentationInfo, select(sliceInfo, -ImgType, -Dataset), c("Id", "Slice")) %>%
  group_by(Id) %>% 
  filter(SliceOrder == min(SliceOrder))

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

if (nrow(promptSlices) > 0) {
  for (s in 1:nrow(promptSlices)) {
    slice <- filter(allSegmentationInfo, 
                    Id == promptSlices$Id[s], 
                    Slice == promptSlices$Slice[s])
    if (nrow(slice > 0)) {
      slice$isLV <- NA
      firstImage <- filter(slice, Time == slice$Time[1])
      if (all(file.exists(unique(getImageFile(firstImage))))) {
        showSingleImage(firstImage)
        
        cat("Id=",promptSlices$Id[s],"Slice=",promptSlices$Slice[s],
            "#Images=",length(unique(slice$Time)),
            "#Segment=",nrow(slice),fill=T)
        
        identifiedLVSegment <- readInteger("Identify segment of left ventricle in first image (0=none, -1=can't see): ")
        if (identifiedLVSegment == -1) {
          # Keep all at NA (we don't know)
        } else {
          if (identifiedLVSegment == 0) {
            # Set only segments of this image
            slice <- slice[Time == slice$Time[1], isLV := FALSE]
          } else {
            # Set segments of other images based to distance to segment in identified image
            # filter out any images with all segments that are too distant from the identified LV in the manually classified image
            lv <- firstImage[firstImage$segIndex == identifiedLVSegment,]
            
            slice$distToLVSeg <- sqrt((slice$m.cx - lv$m.cx)^2 + (slice$m.cy - lv$m.cy)^2)
            identifiedLVSegments <- group_by(slice, Time) %>% summarise(segLV = segIndex[which.min(distToLVSeg)],
                                                                        anyWithinThreshold = any(distToLVSeg < 5))
            slice <- left_join(slice, identifiedLVSegments, by=c("Time"))
            slice$isLV <- (slice$segIndex == slice$segLV)
            slice <- filter(slice, anyWithinThreshold)
            
            # graph of area vs image
            plotSlice(filter(slice, isLV) %>% rename(radius.mean = s.radius.mean,
                                                     radius.max = s.radius.max,
                                                     radius.min = s.radius.min))
            # all images
            showAllSliceImages(slice)
            
            identifiedCorrectly <- readline("Are all segments identified correctly (y/n): ")
            if (identifiedCorrectly != "y") {
              # Set only segments of this image
              slice$isLV <- NA
              slice <- slice[Time == slice$Time[1], isLV := (segIndex == identifiedLVSegment)]
            }
          }
        }
      } else {
        print(unique(getImageFile(firstImage)))
        print("WARN: Not all image files exist on this system")
      }
      View(slice)
      
      # save data
      newData <- filter(slice, !is.na(isLV)) %>% select( 
        Id, Slice, Time,   # uniquely identifies image via join to playlist
        UUID,              # unique segment identifier
        starts_with("m."), # standard moments & shape attributes 
        starts_with("s."),
        isLV)
      
      if(nrow(newData) > 0) {
        if (is.null(segPredictSet)) {
          segPredictSet <- newData
        } else {
          # TODO check both are compatible
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
      } else {
        print("No new data added")
      }
    }
    
    segmentedByID <- left_join(unique(select(sliceInfo, Dataset, Id)), unique(select(segPredictSet, Id, isLV))) %>% group_by(Dataset, Id) %>% summarise(hasSegs = any(isLV)) %>% filter(hasSegs)
    segmentedBySlice <- left_join(unique(select(sliceInfo, Dataset, Id, Slice)), unique(select(segPredictSet, Id, Slice, isLV))) %>% group_by(Dataset, Id, Slice) %>% summarise(hasSegs = any(isLV)) %>% filter(hasSegs)
    
    cat("Classified", nrow(segPredictSet), "segments in", nrow(filter(segPredictSet, isLV)), "Images", fill=T)
    cat("Classified", nrow(unique(select(segmentedByID, Id))), 
        "of", nrow(unique(select(sliceInfo, Id))), "Ids", fill=T)
    cat("Classified", nrow(unique(select(segmentedBySlice, Id, Slice))), 
        "of", nrow(unique(select(sliceInfo, Id, Slice))), "Slices", fill=T)
    
    write.csv(segPredictSet, segPredictFile, row.names=F)
  }
} else {
  print("No new segments to identify!")
}

