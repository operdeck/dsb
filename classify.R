# Manual identification of LV segments with the purpose of creating a prediction dataset
# to do this automatically on the full dataset. 
# Creates "segments-predict.csv" with same data for only the (correctly) identified segments of
# each file / slice / Id.

source("util.R")

segPredictFile <- "segments-predict.csv"

showSegmentLabels <- function(imageDS) {
  for (i in 1:nrow(imageDS)) {
    text(x = imageDS$m.cx[i], y = imageDS$m.cy[i], 
         label = imageDS$segIndex[i], col = "red")
  }
}

# Display and prompt for one slice
showSlice <- function(ds, highlight=F) {
  #   filez <-  ds$file
  #   filez <-  unique(ds$file)
  imgz <- list()
  for (i in 1:nrow(ds)) {
    f <- getImageFile(ds[i,])
    dicom <- readDICOMFile(f)
    img <- Image(normalize(dicom$img))
    if (dim(img)[1] > dim(img)[2]) {
      img <- rotate(img,-90)
    }
    if (!highlight & i == 1) {
      firstImageSegments <- filter(ds, Time == ds$Time[i])
      for (j in 1:nrow(firstImageSegments)) {
        if (firstImageSegments$s.radius.mean[j] >= 1) {
          img <- drawCircle(toRGB(img), 
                            x=firstImageSegments$m.cx[j], y=firstImageSegments$m.cy[j], radius=firstImageSegments$s.radius.mean[j], col="red")
        }
        if (firstImageSegments$s.radius.min[j] >= 1) {
          img <- drawCircle(toRGB(img), 
                            x=firstImageSegments$m.cx[j], y=firstImageSegments$m.cy[j], radius=firstImageSegments$s.radius.min[j], col="orange")
        }
        if (firstImageSegments$s.radius.max[j] >= 1) {
          img <- drawCircle(toRGB(img), 
                            x=firstImageSegments$m.cx[j], y=firstImageSegments$m.cy[j], radius=firstImageSegments$s.radius.max[j], col="orange")
        }
      }
      EBImage::display(img,all=T,method="raster")
      text(dim(img)[1]/2,40,paste("ID",unique(ds$Id),"Slice",unique(ds$Slice)),col="yellow",pos=4)
      showSegmentLabels(firstImageSegments)
    }
    if (highlight) {
      imgz[[f]] <- drawCircle(toRGB(img), 
                              x=ds$m.cx[i], y=ds$m.cy[i], radius=ds$s.radius.mean[i], col="red")
    } else {
      imgz[[f]] <- img
    }
  }
  EBImage::display(EBImage::combine(imgz),all=T,method="raster")
  text(500,40,paste("ID",unique(ds$Id),"Slice",unique(ds$Slice)),col="yellow",pos=4)
  showSegmentLabels(filter(ds, Time == ds$Time[1]))
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
                          Id, Slice, Time,   # uniquely identifies image via join to playlist
                          UUID,              # unique segment identifier
                          starts_with("m."), # standard moments & shape attributes 
                          starts_with("s."))
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

# TODO
# for each slice
#   1. show first image
#      prompt for segment
#      if valid: 
#           show all images of slice, with extrapolated segment (same pos)
#           prompt for ok
#           if ok - classify all segs of all imgs as T/F
#           if not ok - just classify the segs of the first img
#      if no valid seg:
#           classify all segs of first img as F


if (nrow(promptSlices) > 0) {
  for (s in 1:nrow(promptSlices)) {
    slice <- filter(allSegments, Id == promptSlices$Id[s], Slice == promptSlices$Slice[s])
    setkey(slice) # drop keys otherwise unique fails
    slice <- unique(slice) # TODO not sure why there are duplicates
    
    showSlice(slice)
    stop()
    # TODO create model and show predictions
    
    identifiedLVSegment <- readInteger("Identify segment of left ventricle in first image (0=none): ")
    filez <- unique(slice$file)
    slice <- slice[file == filez[1], isLV := (segIndex == identifiedLVSegment)]
    
    # position of identified segment
    identifiedLVSegment.x <- slice$m.cx[slice$file == filez[1] & slice$isLV]
    identifiedLVSegment.y <- slice$m.cy[slice$file == filez[1] & slice$isLV]
    
    # distance of all other segments to the identified one and identify closest one in each image
    slice$distToLVSeg <- sqrt((slice$m.cx - identifiedLVSegment.x)^2 + (slice$m.cy - identifiedLVSegment.y)^2)
    identifiedLVSegments <- group_by(slice, file) %>% summarise(extrapolatedLVSegment = segIndex[which.min(distToLVSeg)])
    slice <- left_join(slice, identifiedLVSegments)
    showSlice(select(filter(slice, segIndex == extrapolatedLVSegment), file, m.cx, m.cy, s.radius.mean, segIndex), 
              highlight=T)
    
    # TODO : prompt if all are ok and fix incorrect ones
    
    # set isLV for all other images in the same slice
    slice <- slice[file != filez[1], isLV := (segIndex == extrapolatedLVSegment)]
    slice <- select(slice, -distToLVSeg, -extrapolatedLVSegment)
    
    # visualize volume for this slice
    print(ggplot(filter(slice, isLV), aes(x=Time, y=sliceVolume, colour=segIndex))+geom_line()+
            ggtitle(paste("Volume over Time for ID",
                          unique(slice$Id),"Slice",unique(slice$Slice))))
    
    # save data
    if (is.null(segPredictSet)) {
      segPredictSet <- filter(slice, !is.na(isLV))
    } else {
      segPredictSet <- rbind(segPredictSet, filter(slice, !is.na(isLV)))
    }
    write.csv(segPredictSet, segPredictFile, row.names=F)
  }
} else {
  print("No new segments to identify!")
}

