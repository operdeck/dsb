# Extract all potentially relevant info from segments from all images
# and write these with all meta info into "segments-<dataset>.csv"
# for later consumption by prediction code.

# update
# source("https://bioconductor.org/biocLite.R")
# biocLite()

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

segmentImagesForOneSlice <- function(imgMetaData, prevROI=NULL) {
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
  
  if (is.null(prevROI)) {
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
  }
  
  sliceSegmentation <- NULL
  for (i in 1:nrow(imgMetaData)) {
    pixelAreaScale <- imgMetaData$PixelSpacing.x[i] * imgMetaData$PixelSpacing.y[i]
    pixelLengthScale <- sqrt(pixelAreaScale)
    
    img_original <- allImages[[i]] 
    
    if (is.null(prevROI)) {    
      img_roi_subtracted <- normalize(abs(img_original - avgImg))
      img_roi_thresholded <- (img_roi_subtracted > otsu(img_roi_subtracted))
      img_roi <- medianFilter(img_roi_thresholded, min(1,round(1.5/pixelLengthScale))) # 2 for "normal images" scale 0.75
      
      # Find the ROI
      roi <- as.data.frame(computeFeatures.moment(img_roi))
      coords <- as.data.frame(pos2coord(pos=which(img_roi>0),dim.mat=dim(img_roi)))
      names(coords) <- c('x','y')
      coords$distToROI <- floor(sqrt((coords$x - roi$m.cx)^2+(coords$y - roi$m.cy)^2))
      roi$radius <- quantile(coords$distToROI, probs=c(0.80)) # important constant: what to consider not moving
    } else {
      roi <- data.frame(m.cx = prevROI$x, m.cy = prevROI$y, radius = prevROI$r)
    }    
    
    if (roi$radius > 1) {
      # Use ROI to normalize image intensity. Using only half radius because often high intensity around borders.
      roi_mask <- matrix(0, nrow=dim(img_original)[1], ncol=dim(img_original)[2])
      roi_mask <- drawCircle(roi_mask, x=roi$m.cx, y=roi$m.cy, roi$radius/2, col=1, fill=T)
      img_masked <- img_original*roi_mask
      scale <- 1/max(img_masked)
      
      #img_original <- img_original*scale # scale original image wrt ROI
      
      maskSize <- 6  #3.75 # Size of erode/dilate mask in millimeter
      kern <- makeBrush(max(1,round(maskSize/pixelLengthScale)), shape= 'Gaussian', sigma=10) # 5 for "normal images" scale 0.75
      
      img_denoised <- normalize(opening(img_original*scale, kern))
      
      otsu <- otsu(img_denoised)
      thresholdStepSize <- 0.2
      thresholds <- seq(otsu, 1, by=thresholdStepSize) # hopefully the first one is good, but never now...
      if (otsu > thresholdStepSize) { thresholds <- c(thresholds, seq(otsu-thresholdStepSize, 0, by=-thresholdStepSize)) }
      
      thresholdIndex <- 0
      doneThresholding <- F
      while (!doneThresholding & (thresholdIndex < length(thresholds))) {
        thresholdIndex <- thresholdIndex + 1
        currentThreshold <- thresholds[thresholdIndex]
        if (thresholdIndex > 1) {
          cat("Re-thresholding",imgMetaData$FileName[i],"#",thresholdIndex,"at",currentThreshold,fill=T)
        }
        
        img_thresholded <- img_denoised > currentThreshold
        img_segmented <- fillHull(bwlabel(img_thresholded))
        
        # Get segment meta data and select only the segments within the ROI
        segmentInfo <- cbind(imgMetaData[i,],
                             data.frame(computeFeatures.moment(img_segmented)),
                             data.frame(computeFeatures.shape(img_segmented)))
        distToROI <- sqrt((segmentInfo$m.cx - roi$m.cx)^2 + (segmentInfo$m.cy - roi$m.cy)^2)
        
        segmentInfo$segIndex <- seq(nrow(segmentInfo))
        segmentInfo$UUID  <- sapply(1:nrow(segmentInfo),UUIDgenerate) # unique ID for each segment
        
        segmentInfo$ROI.x <- round(roi$m.cx)
        segmentInfo$ROI.y <- round(roi$m.cy)
        segmentInfo$ROI.r <- round(roi$radius)
        
        # Only keep segments inside the ROI
        minSegmentArea <- 10 # Minimum segment size in square mm
        segmentInfo <- filter(segmentInfo, distToROI < roi$radius, s.area > minSegmentArea/pixelAreaScale)
        img_colourSegs <- colorLabels(rmObjects(img_segmented, 
                                                setdiff(seq(max(img_segmented)),segmentInfo$segIndex)))
        
        # assume a good segmentation has at least 1 segment
        if (nrow(segmentInfo) >= 1) { doneThresholding <- T}
        
        if (is.null(prevROI)) {
          img_comb <- EBImage::combine(drawCircle((img_colourSegs+toRGB(scale*img_original))/2, x=roi$m.cx, y=roi$m.cy, roi$radius, "yellow", fill=FALSE, z=1),
                                       drawCircle(toRGB(img_original), x=roi$m.cx, y=roi$m.cy, roi$radius, "yellow", fill=FALSE, z=1),
                                       toRGB(img_roi_subtracted),
                                       #toRGB(normalize(img_masked)),
                                       toRGB(img_denoised),
                                       toRGB(img_thresholded),
                                       toRGB(img_roi))
        } else {
          img_comb <- EBImage::combine(drawCircle((img_colourSegs+toRGB(scale*img_original))/2, x=roi$m.cx, y=roi$m.cy, roi$radius, "yellow", fill=FALSE, z=1),
                                       drawCircle(toRGB(img_original), x=roi$m.cx, y=roi$m.cy, roi$radius, "yellow", fill=FALSE, z=1),
                                       toRGB(img_denoised),
                                       toRGB(img_thresholded))
        }
        display(img_comb,all=T,method="raster")
        text(10,20,imgMetaData$FileName[i],col="yellow",pos=4)
      
      }
      
      if (nrow(segmentInfo) > 0) {
        dirName <- getSegmentedImageDir(segmentInfo[1,])
        dir.create(dirName, showWarnings = F, recursive = T)
        fName <- getSegmentedImageFile(segmentInfo[1,], dirName)
        writeImage(img_colourSegs, fName)
      }
      
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
                         mutate(unique(select(allSegments, Id, Slice, Time)), isProcessed=T),
                         by=c("Id", "Slice", "Time"))
  imageList$isProcessed <- ifelse(is.na(imageList$isProcessed), F, T)
} else {
  imageList$isProcessed <- F
}
imageList$Random <- runif(max(imageList$Id))[imageList$Id] # keep ID's together

# imageList$SpecialOrder <- ifelse(imageList$SliceOrder > 1, imageList$SliceOrder, 10) # first last - temporary
# before random: -SliceOrder, Dataset
imageList <- arrange(imageList, isProcessed, Random) %>% select(-Random)

# Process images per slice. Image of the same slice (usually) have same dimensions, location etc

sliceList <- unique(select(imageList, Dataset, Id, ImgType, starts_with("Slice")))

for (nextSlice in seq(nrow(sliceList))) {
  slice <- sliceList[nextSlice,]
  sliceImgMetaData <- left_join(select(slice, Dataset, Id, ImgType, Slice), 
                                imageList,
                                by=c("Dataset", "Id", "ImgType", "Slice")) %>% arrange(Time)
  
  if (!all(file.exists(paste("data",sliceImgMetaData$FileName,sep="/")))) { 
    # just in case not all files are on this filesystem
    next
  }
  cat("Processing",slice$Dataset,slice$Id,"slice",slice$Slice,
      paste("(",slice$SliceIndex, "/", slice$SliceCount, " order:", slice$SliceOrder, ")", sep=""),
      fill=T)

  # ROI of all slices of the same Id with lower order (closer to the middle)
  ROIs <- filter(left_join(select(filter(allSegments, Id==slice$Id, Slice!=slice$Slice), 
                                  Id, Slice, starts_with("ROI.")),
                           imageList, by=c("Id", "Slice")), 
                 SliceOrder < slice$SliceOrder)
  if (nrow(ROIs)>0) {
    roi <- list(x = mean(ROIs$ROI.x), y = mean(ROIs$ROI.y), r = mean(ROIs$ROI.r))
  } else {
    roi <- NULL
  }
  
  sliceSegmentationInfo <- segmentImagesForOneSlice(sliceImgMetaData, roi)
  # drop meta-data as this will be joined in via image list upon read
  sliceSegmentationInfo <- select(sliceSegmentationInfo, 
                                  Id, Slice, Time, segIndex, UUID,
                                  starts_with("m."), starts_with("s."),
                                  starts_with("ROI."))
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

    # temp add Dataset to the segments for reporting and progress info
    allSegments <- left_join(allSegments, unique(select(imageList, Id, Slice, Time, Dataset)), by=c("Id","Slice","Time"))
    cat("...write", getSegmentFile(ds), "id's:",length(unique(filter(allSegments, Dataset==ds)$Id)), fill=T)
    
    write.csv(select(filter(allSegments, Dataset==ds), -Dataset), getSegmentFile(ds), row.names=F)
    
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
    
    # drop Dataset again - was only temp
    allSegments <- select(allSegments, -Dataset)
  }
}

# write final results (again, just to be sure)
for (ds in unique(imageList$Dataset)) {
  cat("Write final", getSegmentFile(ds), "id's:",length(unique(filter(allSegments, Dataset==ds)$Id)), fill=T)
  write.csv(filter(allSegments, Dataset==ds), getSegmentFile(ds), row.names=F)
}




