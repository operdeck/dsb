# Manual identification of LV segments with the purpose of creating a prediction dataset
# to do this automatically on the full dataset. 
# Creates "segments-predict.csv" with same data for only the (correctly) identified segments of
# each file / slice / Id.

source("util.R")

# TODO
#
# Manually incorrectly identified - add some 'shitlist' of Id/Slices to be re-done:
#  train 7, slice 49 (Time 13 incorrect)
#  train 11, slice 14 (Couple of Times, so small)
#
# TODO maybe also a re-confirm after the graph
# maybe use additional criterion of dist to LV < 2 (should really be very close), filter
# out the slice images that do not match this criterion (isLV is just unknown for them)
#
# Incorrect segmentations by segmentation code
#  train 6, slice 10 - too much light around edges
#  train 8, slice 58 - same
#  train 8, slice 59 - same
#  train 3, slice 46&47 - very dark, some images ok, but some lack segments



segPredictFile <- "segments-predict.csv"

showSegmentLabels <- function(imageDS) {
  for (i in 1:nrow(imageDS)) {
    text(x = imageDS$m.cx[i], y = imageDS$m.cy[i], 
         label = imageDS$segIndex[i], col = "red")
  }
}

showSingleImage <- function(ds) {
  f <- getImageFile(ds[1,])
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


# Read all segment info from all datasets
allSegments <- NULL
for (dataset in datasetFoldersForSegmentDetection) {
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

# Read playlist, which is a full list of all slices for all cases of all datasets
# which provides additional (meta) info about the slices
sliceInfo <- fread("slicelist.csv")

# Select only the middle slices (usually 2 per Id)
# TODO consider a wider selection
allSegments <- left_join(allSegments, select(sliceInfo, -ImgType, -Dataset), c("Id", "Slice")) %>%
  group_by(Id) %>% 
  filter(SliceOrder == min(SliceOrder))

# Select slices to prompt for
promptSlices <- unique(select(allSegments, Id, Slice)) 
if (!is.null(segPredictSet)) {
  promptSlices <- left_join(promptSlices, 
                            mutate(unique(select(segPredictSet, Id, Slice)), classifiedAlready=T)) %>%
    mutate(classifiedAlready = !is.na(classifiedAlready)) %>%
    arrange(classifiedAlready, Id, Slice)
}

if (nrow(promptSlices) > 0) {
  for (s in 1:nrow(promptSlices)) {
    slice <- filter(allSegments, Id == promptSlices$Id[s], Slice == promptSlices$Slice[s])
    if (nrow(slice > 0)) {
      slice$isLV <- NA
      firstImage <- filter(slice, Time == slice$Time[1])
      showSingleImage(firstImage)
      identifiedLVSegment <- readInteger("Identify segment of left ventricle in first image (0=none): ")
      if (identifiedLVSegment == 0) {
        # Set only segments of this image
        slice <- slice[Time == slice$Time[1], isLV := FALSE]
      } else {
        # Set segments of other images based to distance to segment in identified image
        lv <- firstImage[firstImage$segIndex == identifiedLVSegment,]
        
        slice$distToLVSeg <- sqrt((slice$m.cx - lv$m.cx)^2 + (slice$m.cy - lv$m.cy)^2)
        identifiedLVSegments <- group_by(slice, Time) %>% summarise(segLV = segIndex[which.min(distToLVSeg)])
        slice <- left_join(slice, identifiedLVSegments, by=c("Time"))
        slice$isLV <- (slice$segIndex == slice$segLV)
        
        showAllSliceImages(slice)
        identifiedCorrectly <- readline("Are all segments identified correctly (y/n): ")
        if (identifiedCorrectly != "y") {
          # Set only segments of this image
          slice$isLV <- NA
          slice <- slice[Time == slice$Time[1], isLV := (segIndex == identifiedLVSegment)]
        }
      }
    }
    View(slice)
    
    # visualize volume for this slice
    if (nrow(filter(slice, isLV))>0) {
      print(ggplot(filter(slice, isLV), aes(x=Time, y=s.area, colour=segIndex))+geom_line()+
              ggtitle(paste("Segment area over Time for ID",
                            unique(slice$Id),"Slice",unique(slice$Slice))))
    }
    
    # save data
    newData <- filter(slice, !is.na(isLV)) %>% select( 
      Id, Slice, Time,   # uniquely identifies image via join to playlist
      UUID,              # unique segment identifier
      starts_with("m."), # standard moments & shape attributes 
      starts_with("s."),
      isLV)
    
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
      
      segPredictSet <- rbind(segPredictSet, newData)
    }
    
    segmentedByID <- left_join(unique(select(playlist, Dataset, Id)), unique(select(segPredictSet, Id, isLV))) %>% group_by(Dataset, Id) %>% summarise(hasSegs = any(isLV)) %>% filter(hasSegs)
    segmentedBySlice <- left_join(unique(select(playlist, Dataset, Id, Slice)), unique(select(segPredictSet, Id, Slice, isLV))) %>% group_by(Dataset, Id, Slice) %>% summarise(hasSegs = any(isLV)) %>% filter(hasSegs)
    
    cat("Segmented by case : ", nrow(unique(select(segmentedByID, Id))), "of", nrow(unique(select(playlist, Id))), fill=T)
    cat("Segmented by slice: ", nrow(unique(select(segmentedBySlice, Id, Slice))), "of", nrow(unique(select(playlist, Id, Slice))), fill=T)
    
    write.csv(segPredictSet, segPredictFile, row.names=F)
  }
} else {
  print("No new segments to identify!")
}

