# Extract all potentially relevant info from segments from all images
# and write these with all meta info into "segments-<dataset>.csv"
# for later consumption by prediction code.

# Segmentation
# - use ROI of middle slices (x 1.5 or 2 perhaps) for the outer slices
# - for ideas on how to color EBimage segments: http://rpackages.ianhowson.com/bioc/EBImage/man/bwlabel.html
# - adaptive thresholding: http://rpackages.ianhowson.com/bioc/EBImage/man/thresh.html
# - segm: maybe use watershed (voronoi not a good idea)

# 80-8 is right scale
# 500-20 has same scale for img processing constants

# Id 80 - large scale
# No segmentation (too dark)
#  train 66 - 15
#  train  6 - 7
# TODO: get this automatically from the segmentation (processed but no segs)
# TODO: save segmented images

source("util.R")

# see bio image detection stuff
# http://bioconductor.wustl.edu/bioc/vignettes/EBImage/inst/doc/AnalysisWithEBImage.pdf

segmentImagesForOneSlice <- function(imgMetaData) {
  allImages <- vector(mode = "list", length = nrow(imgMetaData))
  for (i in seq(nrow(imgMetaData))) {
    f <- paste("data",imgMetaData$FileName[i],sep="/")
    if (!file.exists(f)) {
      print(imgMetaData[i,])
      print(f)
      stop("File doesnt exist - formatting issues or playlist out of sync with local filesystem?")
    }
    dicom <- readDICOMFile(f)
    img <- Image(normalize(dicom$img))
    if (dim(img)[1] > dim(img)[2]) {
      img <- rotate(img,-90)
    }
    allImages[[i]] <- img
  }
  # Calculate 'average' image
  isFirst <- T
  for (i in allImages) {
    if (isFirst) {
      avgImg <- i
      isFirst <- F
    } else {
      avgImg <- avgImg + i
    }
  }
  avgImg <- avgImg / length(allImages)
  
  sliceSegmentation <- NULL
  for (i in 1:nrow(imgMetaData)) {
    pixelAreaScale <- imgMetaData$PixelSpacing.x[i] * imgMetaData$PixelSpacing.y[i]
    pixelLengthScale <- sqrt(pixelAreaScale)
    
    img_original <- allImages[[i]] 
    img_subtracted <- normalize(abs(img_original - avgImg))
    img_roi_thresholded <- (img_subtracted > otsu(img_subtracted))
    # 
    img_roi <- medianFilter(img_roi_thresholded, min(1,round(1.5/pixelLengthScale))) # 2 for "normal images" scale 0.75
    
    # Find the ROI
    roi = as.data.frame(computeFeatures.moment(img_roi))
    coords <- as.data.frame(pos2coord(pos=which(img_roi>0),dim.mat=dim(img_roi)))
    names(coords) <- c('x','y')
    coords$distToROI <- floor(sqrt((coords$x - roi$m.cx)^2+(coords$y - roi$m.cy)^2))
    roi$radius <- quantile(coords$distToROI, probs=c(0.90)) # important constant: what to consider not moving
    if (roi$radius > 1) {
      roi_mask <- matrix(0, nrow=dim(img_roi)[1], ncol=dim(img_roi)[2])
      roi_mask <- drawCircle(roi_mask, x=roi$m.cx, y=roi$m.cy, roi$radius, col=1, fill=T)
      img_masked <- img_original*roi_mask
      scale <- 1/max(img_masked)
      img_original <- img_original*scale # scale original image wrt ROI
      
      kern <- makeBrush(min(1,round(3.75/pixelLengthScale)), shape= 'disc') # 5 for "normal images" scale 0.75
      
      img_blurred <- normalize(opening(img_original, kern))
      img_thresholded <- img_blurred > otsu(img_blurred) # Otsu??s threshold 
      img_segmented <- fillHull(bwlabel(img_thresholded))
      
      # Get segment meta data and select only the segments within the ROI
      segmentInfo <- cbind(imgMetaData[i,],
                           data.frame(computeFeatures.moment(img_segmented)),
                           data.frame(computeFeatures.shape(img_segmented)))
      segmentInfo$distToROI <- sqrt((segmentInfo$m.cx - roi$m.cx)^2 + (segmentInfo$m.cy - roi$m.cy)^2)
      segmentInfo$segIndex <- seq(nrow(segmentInfo))
      segmentInfo$UUID <- sapply(1:nrow(segmentInfo),UUIDgenerate) # unique ID for each segment
      
      # Only keep segments inside the ROI
      segmentInfo <- filter(segmentInfo, distToROI < roi$radius)
      img_colourSegs <- colorLabels(rmObjects(img_segmented, 
                                              setdiff(seq(max(img_segmented)),segmentInfo$segIndex)))
      
      img_comb <- EBImage::combine(drawCircle(img_colourSegs, x=roi$m.cx, y=roi$m.cy, roi$radius, "yellow", fill=FALSE, z=1),
                                   drawCircle(toRGB(img_original), x=roi$m.cx, y=roi$m.cy, roi$radius, "yellow", fill=FALSE, z=1),
                                   toRGB(img_subtracted),
                                   #toRGB(img_masked),
                                   toRGB(img_blurred),
                                   toRGB(img_thresholded),
                                   toRGB(img_roi))
      display(img_comb,all=T,method="raster")
      text(10,20,imgMetaData$FileName[i],col="yellow",pos=4)
      
      if (is.null(sliceSegmentation)) {
        sliceSegmentation <- segmentInfo
      } else {
        sliceSegmentation <- rbind(sliceSegmentation, segmentInfo)
      }
    } else {
      print("WARN: invalid ROI for segment (too small perhaps?)")
    }
  } # all images of one slice
  
  return(sliceSegmentation)
}

print("Reading image meta data")
imageList <- getImageList()

print("Reading previous segmentation")
allSegments <- NULL
for (dataset in unique(imageList$Dataset)) {
  if (file.exists(getSegmentFile(dataset))) {
    segmentsPerDataset <- fread(getSegmentFile(dataset))
    if (is.null(allSegments)) {
      allSegments <- segmentsPerDataset
    } else {
      removedSet <- setdiff(names(allSegments), names(segmentsPerDataset))
      addedSet <- setdiff(names(segmentsPerDataset), names(allSegments))
      diffSet <- paste(c(paste("-",removedSet), paste("+",addedSet)),collapse=", ")
      if (length(removedSet) + length(addedSet) > 0) {
        print(names(allSegments))
        print(names(segmentsPerDataset))
        print(diffSet)
        stop("Datasets do not match up. Please consider removing files.")
      }
      allSegments <- rbind(allSegments, segmentsPerDataset)
    }
  }
}


print("Starting to process images")
# Here we order the imageList:
# - first incomplete ID's, then unprocessed slices by their special 'mid order'
if (!is.null(allSegments)) {
  imageList <- left_join(imageList, 
                         mutate(unique(select(allSegments, Id, Dataset, ImgType, Slice, Time)), isProcessed=T),
                         by=c("Dataset", "Id", "ImgType", "Slice", "Time"))
  imageList$isProcessed <- ifelse(is.na(imageList$isProcessed), F, T)
} else {
  imageList$isProcessed <- F
}
imageList$Random <- runif(max(imageList$Id))[imageList$Id] # keep ID's together
imageList <- arrange(imageList, isProcessed, SliceOrder, Dataset, Random) %>% select(-Random)

# Process images per slice. Image of the same slice (usually) have same dimensions, location etc
sliceList <- unique(select(imageList, Dataset, Id, ImgType, starts_with("Slice")))
for (nextSlice in seq(nrow(sliceList))) {
  slice <- sliceList[nextSlice,]
  sliceImgMetaData <- left_join(select(slice, Dataset, Id, ImgType, Slice), 
                                imageList,
                                by=c("Dataset", "Id", "ImgType", "Slice")) %>% arrange(Time)
  
  if (all(file.exists(paste("data",sliceImgMetaData$FileName,sep="/")))) { # just in case not all files are on this filesystem
    cat("Processing",slice$Dataset,slice$Id,"slice",slice$Slice,
        paste("(",slice$SliceIndex, "/", slice$SliceCount, " order:", slice$SliceOrder, ")", sep=""),
        fill=T)
    
    sliceSegmentationInfo <- segmentImagesForOneSlice(sliceImgMetaData)
    cat("...",nrow(sliceSegmentationInfo),"segments in",nrow(sliceImgMetaData),"images",fill=T)
    
    if (is.null(allSegments)) {
      allSegments <- sliceSegmentationInfo
    } else {
      removedSet <- setdiff(names(allSegments), names(sliceSegmentationInfo))
      addedSet <- setdiff(names(sliceSegmentationInfo), names(allSegments))
      diffSet <- paste(c(paste("-",removedSet), paste("+",addedSet)),collapse=", ")
      if (length(removedSet) + length(addedSet) > 0) {
        print("Existing:")
        print(names(allSegments))
        print("New:")
        print(names(sliceSegmentationInfo))
        print(paste(diffSet, collapse=","))
        stop("Datasets do not match up. Please consider removing segment files.")
      }
      
      nDrop <- nrow(filter(allSegments, 
                           Slice %in% sliceSegmentationInfo$Slice,
                           Id %in% sliceSegmentationInfo$Id))
      if (nDrop > 0) {
        cat("...dropping",nDrop,"rows from previous results",fill=T)
      }
      
      allSegments <- rbind(filter(allSegments, 
                                          !(Slice %in% sliceSegmentationInfo$Slice & 
                                              Id %in% sliceSegmentationInfo$Id)), 
                                   sliceSegmentationInfo)
    }
    
    # Add segmentation to results and write file (once in a while)
    if ((nextSlice %% 10 == 0) | 
          (nextSlice == nrow(sliceList)) | 
          ((nextSlice < nrow(sliceList)) & (sliceList$Dataset[nextSlice+1] != sliceList$Dataset[nextSlice]))) {
      ds <- sliceList$Dataset[nextSlice]
      cat("...write", getSegmentFile(ds), "id's:",length(unique(filter(allSegments, Dataset==ds)$Id)), fill=T)
      write.csv(filter(allSegments, Dataset==ds), getSegmentFile(ds), row.names=F)
      
      print("...completeness:")
      print(group_by(left_join(sliceList, 
                               mutate(unique(select(allSegments, Id, Slice)), isSegmented=TRUE), 
                               by=c("Id", "Slice")), Dataset) %>% 
              summarise(complete = paste(round(100*sum(isSegmented, na.rm=T)/n(),2),"%",sep="" )))
      
      # Show progress
      ds <- group_by(left_join(sliceList, 
                               mutate(unique(select(allSegments, Id, Slice)), isSegmented=TRUE), 
                               by=c("Id", "Slice")), Dataset, Id) %>% 
        summarise(complete=round(100*sum(isSegmented, na.rm=T)/n())) %>%
        group_by(Dataset, complete) %>% summarise(n = n())
      print(ggplot(ds, aes(x=complete, y=n, fill=Dataset))+geom_bar(stat="identity")+
              ggtitle("Completeness of segmentation per Id"))
    }
  }
}

# write final results (again, just to be sure)
for (ds in unique(imageList$Dataset)) {
  cat("Write final", getSegmentFile(ds), "id's:",length(unique(filter(allSegments, Dataset==ds)$Id)), fill=T)
  write.csv(filter(allSegments, Dataset==ds), getSegmentFile(ds), row.names=F)
}




