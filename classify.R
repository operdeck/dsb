# Manual identification of LV segments with the purpose of creating a prediction dataset
# to do this automatically on the full dataset. 
# Creates "segments-predict.csv" with same data for only the (correctly) identified segments of
# each file / slice / Id.

source("util.R")

showSegmentLabels <- function(imageDS) {
  for (i in 1:nrow(imageDS)) {
    text(x = imageDS$m.cx[i], y = imageDS$m.cy[i], 
         label = imageDS$segIndex[i], col = "red")
  }
}

showSlice <- function(ds, highlight=F) {
  #   filez <-  ds$file
  #   filez <-  unique(ds$file)
  imgz <- list()
  for (i in 1:nrow(ds)) {
    f <- ds$file[i]
    dicom <- readDICOMFile(f)
    img <- Image(normalize(dicom$img))
    if (dim(img)[1] > dim(img)[2]) {
      img <- rotate(img,-90)
    }
    
    if (!highlight & i == 1) {
      firstImage <- filter(ds, file==ds$file[1])
      # TODO highlight actual segmentation, not just the circles
      for (i in 1:nrow(firstImage)) {
        if (firstImage$s.radius.mean[i] >= 1) {
          img <- drawCircle(toRGB(img), 
                            x=firstImage$m.cx[i], y=firstImage$m.cy[i], radius=firstImage$s.radius.mean[i], col="red")
        }
        if (firstImage$s.radius.min[i] >= 1) {
          img <- drawCircle(toRGB(img), 
                            x=firstImage$m.cx[i], y=firstImage$m.cy[i], radius=firstImage$s.radius.min[i], col="orange")
        }
        if (firstImage$s.radius.max[i] >= 1) {
          img <- drawCircle(toRGB(img), 
                            x=firstImage$m.cx[i], y=firstImage$m.cy[i], radius=firstImage$s.radius.max[i], col="orange")
        }
      }
      EBImage::display(img,all=T,method="raster")
      text(dim(img)[1]/2,40,paste("ID",unique(ds$Id),"Slice",unique(ds$Slice)),col="yellow",pos=4)
      showSegmentLabels(firstImage)
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
  showSegmentLabels(filter(ds, file==ds$file[1]))
}


# Read all segment info from all datasets
allSegments <- NULL
for (dataset in c('train','validate','test')) {
  fname <- paste("segments-",dataset,".csv",sep="")
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
segPredictFile <- "segments-predict.csv"
if (file.exists(segPredictFile)) {
  segPredictSet <- fread(segPredictFile)
  if (length(setdiff(names(allSegments), names(segPredictSet))) != 0 |
        length(setdiff(setdiff(names(segPredictSet), names(allSegments)), 
                       c("isLV", "midSlice"))) != 0) {
    print("Names have changed - existing segment prediction set is void")
    cat("New fields:", setdiff(names(allSegments), names(segPredictSet)), fill=T)
    cat("Prediction set adds:", setdiff(setdiff(names(segPredictSet), names(allSegments)), 
                                        c("isLV", "midSlice")), fill=T)
    segPredictSet <- NULL
  }
}

# Select only the middle segment for each Id (TODO: consider broader selection)
slicesToConsider <- group_by(allSegments, Id) %>% 
  summarise(midSlice = sort(unique(Slice))[ceiling(length(unique(Slice))/2)]) 
if (! is.null(segPredictSet)) {
  allSegments <- left_join(filter(allSegments, !(UUID %in% segPredictSet$UUID)), slicesToConsider) %>%
    filter(Slice == midSlice)
} else {
  allSegments <- left_join(allSegments, slicesToConsider) %>%
    filter(Slice == midSlice)
}


# Get input for these slices
promptSlices <- unique(select(allSegments, Id, Slice)) 
if (nrow(promptSlices) > 0) {
  for (s in 1:nrow(promptSlices)) {
    slice <- filter(allSegments, Id == promptSlices$Id[s], Slice == promptSlices$Slice[s])
    setkey(slice) # drop keys otherwise unique fails
    slice <- unique(slice) # TODO not sure why there are duplicates
    
    showSlice(slice)
    
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

