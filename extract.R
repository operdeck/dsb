# Extract all potentially relevant info from segments from all images
# and write these with all meta info into "segments-<dataset>.csv"
# for later consumption by prediction code.

# update
# source("https://bioconductor.org/biocLite.R")
# biocLite()

# Segmentation
# - for ideas on how to color EBimage segments: http://rpackages.ianhowson.com/bioc/EBImage/man/bwlabel.html
# - adaptive thresholding: http://rpackages.ianhowson.com/bioc/EBImage/man/thresh.html
# - segm: maybe use watershed (voronoi not a good idea)

source("util.R")

if (file.exists('lv.xgb.model')) {
  leftVentricleSegmentModel <- xgb.load('lv.xgb.model')
} else {
  leftVentricleSegmentModel <- NULL
}
thresholdHighLV <- 0.8 # used to decide to settle on a threshold and not try further, plus
                       # to get the segments from previous slices to determine the ROI
thresholdLowLV <- 0.4 # used to decide to re-segment image with full ROI

# see bio image detection stuff
# http://bioconductor.wustl.edu/bioc/vignettes/EBImage/inst/doc/AnalysisWithEBImage.pdf

getAverageImage <- function(allImages) {
  avgImg <- NULL
  for (i in allImages) {
    if (is.null(avgImg)) {
      avgImg <- i
    } else {
      avgImg <- avgImg + i
    }
  }
  avgImg <- avgImg / length(allImages)
  return(avgImg)  
}

readAllImages <- function(metaData) {
  allImages <- vector(mode = "list", length = nrow(metaData))
  for (i in seq(nrow(metaData))) {
    f <- paste("data",metaData$FileName[i],sep="/")
    if (!file.exists(f)) {
      print(metaData[i,])
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
  return(allImages)
}

findROI <- function (allImages, imgMetaData) {
  avgImg <- getAverageImage(allImages)
  
  display(avgImg,all=T,method="raster")
  text(10,20,paste("ID",unique(imgMetaData$Id),"slice",unique(imgMetaData$Slice)),col="yellow",pos=4)
  
  rois <- data.frame(m.cx = rep(0,nrow(imgMetaData)), m.cy = rep(0,nrow(imgMetaData)), radius = rep(0,nrow(imgMetaData)))
  for (i in 1:nrow(imgMetaData)) {
    pixelAreaScale <- imgMetaData$PixelSpacing.x[i] * imgMetaData$PixelSpacing.y[i]
    pixelLengthScale <- sqrt(pixelAreaScale)
    img_roi_subtracted <- normalize(abs(allImages[[i]] - avgImg))
    img_roi_thresholded <- (img_roi_subtracted > otsu(img_roi_subtracted))
    img_roi <- medianFilter(img_roi_thresholded, min(1,round(1.5/pixelLengthScale))) # 2 for "normal images" scale 0.75
    
    # Find the ROI
    moments <- computeFeatures.moment(img_roi)
    rois$m.cx[i] <- moments[1]
    rois$m.cy[i] <- moments[2]
    
    coords <- as.data.frame(pos2coord(pos=which(img_roi>0),dim.mat=dim(img_roi)))
    names(coords) <- c('x','y')
    coords$distToROI <- floor(sqrt((coords$x - moments[1])^2+(coords$y - moments[2])^2))
    rois$radius[i] <- quantile(coords$distToROI, probs=c(0.80)) # important constant: what to consider not moving
    
    symbols(x=moments[1],y=moments[2],circles=rois$radius[i],fg="blue",add=T)
  }
  # TODO - perhaps min, or overlap of rois is better measure, or a .95 value or something,
  # maybe drop the outliers (lofactor)
  roi <- list(x = mean(rois$m.cx), y = mean(rois$m.cy), r = mean(rois$radius))
}

segmentImagesForOneSlice <- function(imgMetaData, roi=NULL, roiAlt=NULL) {
  imgMetaData$best_probLV <- -Inf
  allImages <- readAllImages(imgMetaData)  
  hasPredefinedROI <- !is.null(roi)
  if (!hasPredefinedROI) {    
    roi <- findROI(allImages, imgMetaData)
  }    
  sliceSegmentation <- NULL
  for (i in 1:nrow(imgMetaData)) {
    pixelAreaScale <- imgMetaData$PixelSpacing.x[i] * imgMetaData$PixelSpacing.y[i]
    pixelLengthScale <- sqrt(pixelAreaScale)
    img_original <- allImages[[i]] 
    if (roi$r > 1) {
      roi.r.scale <- roi$r/2
      # Normalize intensity using ROI. Using only half radius because often high intensity around borders.
      intensity_mask <- drawCircle(matrix(0, nrow=dim(img_original)[1], ncol=dim(img_original)[2]), 
                                   x=roi$x, y=roi$y, roi.r.scale, col=1, fill=T)
      img_masked_to_roi <- img_original*intensity_mask # keeps only pixels inside ROI
      scale <- 1/max(img_masked_to_roi)
      
      # Clip image using ROI. Using a bit larger area than ROI itself.
      if (!is.null(roiAlt)) roi <- roiAlt
      roi.r.clip <- 1.5*roi$r
      clip_mask <- matrix(0, nrow=dim(img_original)[1], ncol=dim(img_original)[2])
      clip_mask[ max(1, round(roi$x - roi.r.clip)) : min(nrow(clip_mask), round(roi$x + roi.r.clip)),
                 max(1, round(roi$y - roi.r.clip)) : min(ncol(clip_mask), round(roi$y + roi.r.clip))] <- 1
      #clip_mask2 <- drawCircle(matrix(0, nrow=dim(img_original)[1], ncol=dim(img_original)[2]), 
      #                        x=roi$x, y=roi$y, roi.r.clip, col=1, fill=T)
      img_clipped <- img_original*scale*clip_mask
      
      # Filter image
      blurBrushSize <- 3.75 # Size of erode/dilate mask in millimeter
      blurBrush <- makeBrush(max(1,round(blurBrushSize/pixelLengthScale)), shape= 'Gaussian', sigma=10) # 5 for "normal images" scale 0.75
      img_filtered <- normalize(opening(img_clipped),blurBrush) #normalize(opening(img_original*scale, blurBrush))
      # Threshold and segment image. Will do multiple iterations if necessary.
      #plot(hist(img_filtered,xlim = c(0, 1),box = T))
      otsu_threshold <- otsu(img_filtered)
      thresholds <- c(seq(otsu_threshold, 1, by=0.1), seq(otsu_threshold, 0, by=-0.05)) # otsu first, then the rest
      best_img_thresholded <- NULL
      best_img_colourSegs <- NULL
      best_segmentInfo <- NULL
      best_probLV <- -Inf
      best_segIndex <- NULL
      for (thresholdIndex in seq(length(thresholds))) {
        currentThreshold <- thresholds[thresholdIndex]
        #cat("Threshold=",currentThreshold,"[",thresholdIndex, "]",fill=T)
        img_thresholded <- img_filtered > currentThreshold
        img_segmented <- fillHull(bwlabel(img_thresholded))
        if (max(img_segmented) > 0) {
          # Get segment meta data and select only the segments within the ROI
          segmentInfo <- cbind(imgMetaData[i,],
                               data.frame(computeFeatures.moment(img_segmented)),
                               data.frame(computeFeatures.shape(img_segmented)))
          distToROI <- sqrt((segmentInfo$m.cx - roi$x)^2 + (segmentInfo$m.cy - roi$y)^2)

          segmentInfo$segIndex <- seq(nrow(segmentInfo))
          segmentInfo$UUID  <- sapply(1:nrow(segmentInfo),UUIDgenerate) # unique ID for each segment
          segmentInfo$ROI.x <- round(roi$x)
          segmentInfo$ROI.y <- round(roi$y)
          segmentInfo$ROI.r <- round(roi$r)
  
          # Only keep segments inside the ROI
          segmentInfo <- filter(segmentInfo, distToROI < roi$r)
          img_colourSegs <- colorLabels(rmObjects(img_segmented, 
                                                  setdiff(seq(max(img_segmented)),segmentInfo$segIndex)))
          # Use LV model (if available) to evaluate quality of segmentation
          if (!is.null(leftVentricleSegmentModel)) {
            predictors <- createSegmentPredictSet(select(segmentInfo, -isProcessed, -best_probLV))
            #cat("*** Predictors:",names(predictors),fill=T)
            probLV <- round(predict(leftVentricleSegmentModel, data.matrix(predictors), missing=NaN),3)
            #names(probLV) <- segmentInfo$segIndex
            #print(head(sort(probLV, decreasing=T), 5))
            
            if (max(probLV) > best_probLV) {
              best_segIndex <- segmentInfo$segIndex[which.max(probLV)]
              best_img_thresholded <- img_thresholded
              best_img_colourSegs <- img_colourSegs
              best_segmentInfo <- segmentInfo
              best_probLV <- max(probLV)
            }
            if (max(probLV) > thresholdHighLV) break
          } else {
            best_img_thresholded <- img_thresholded
            best_img_colourSegs <- img_colourSegs
            best_segmentInfo <- segmentInfo
            if (nrow(segmentInfo) >= 2) break
          }
        }
      }
      #cat("Best pLV:", best_probLV, "[#", thresholdIndex, ": @", thresholds[thresholdIndex], "]", fill=T)
      imgMetaData$best_probLV[i] <- best_probLV
      
      # display result of best thresholding
      annotatedOriginal <- drawCircle((best_img_colourSegs+toRGB(scale*img_original))/2, 
                                      x=roi$x, y=roi$y, roi$r, "blue", fill=FALSE, z=1)
      if (!is.null(roiAlt)) {
        annotatedOriginal <- drawCircle(annotatedOriginal, 
                                        x=roiAlt$x, y=roiAlt$y, roiAlt$r, "green", fill=FALSE, z=1)
      }
      if (!hasPredefinedROI) {
        img_comb <- EBImage::combine(toRGB(best_img_thresholded), # will get segment indices overlayed
                                     annotatedOriginal,
                                     drawCircle(drawCircle(drawCircle(toRGB(img_original),
                                                                      x=roi$x, y=roi$y, roi$r, "blue", fill=FALSE, z=1),
                                                           x=roi$x, y=roi$y, roi.r.clip, "red", fill=FALSE, z=1),
                                                x=roi$x, y=roi$y, roi.r.scale, "yellow", fill=FALSE, z=1),
                                     best_img_colourSegs)
      } else {
        img_comb <- EBImage::combine(toRGB(best_img_thresholded), # will get segment indices overlayed
                                     annotatedOriginal)
      }
      display(img_comb,all=T,method="raster")
      text(10,20,imgMetaData$FileName[i],col="yellow",pos=4)
      showSegmentLabels(best_segmentInfo)
      if (!is.null(best_segIndex) & nrow(best_segmentInfo) > 0) {
        showSegmentLabels(best_segmentInfo[best_segmentInfo$segIndex==best_segIndex,],textColour="green")
      }
      
      if (!is.null(best_segmentInfo) && nrow(best_segmentInfo) > 0) {
        dirName <- getSegmentedImageDir(best_segmentInfo[1,])
        dir.create(dirName, showWarnings = F, recursive = T)
        fName <- getSegmentedImageFile(best_segmentInfo[1,], dirName)
        writeImage(img_colourSegs, fName)
      }
      
      if (!is.null(best_segmentInfo)) {
        if (is.null(sliceSegmentation)) {
          sliceSegmentation <- best_segmentInfo
        } else {
          sliceSegmentation <- rbind(sliceSegmentation, best_segmentInfo)
        }
      }
    } else {
      print("WARN: invalid ROI for segment (too small perhaps?)")
    }
  } # all images of one slice

  cat("pLV summary for this slice (best pLV for every image):",fill=T)
  print(summary(imgMetaData$best_probLV))
  if (!is.null(roiAlt) & mean(imgMetaData$best_probLV) < thresholdLowLV) {
    print("pLV below threshold - re-segment with ROI discovery")
    sliceSegmentation <- segmentImagesForOneSlice(imgMetaData)
  }
  
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
imageList$Random <- runif(nrow(imageList))
# imageList$SpecialOrder <- ifelse(T, imageList$SliceOrder, 1+max(imageList$SliceOrder))
# imageList <- arrange(imageList, isProcessed, SpecialOrder, Random) %>% select(-Random, -SpecialOrder)
imageList <- arrange(imageList, isProcessed, Random) %>% select(-Random)

# Process images per slice. Image of the same slice (usually) have same dimensions, location etc
sliceList <- unique(select(imageList, Dataset, Id, ImgType, starts_with("Slice")))
print(head(paste(sliceList$Id,sliceList$Slice,sep="/"),20))

for (nextSlice in seq(nrow(sliceList))) {
  slice <- sliceList[nextSlice,]
  sliceImgMetaData <- left_join(select(slice, Dataset, Id, ImgType, Slice), 
                                imageList,
                                by=c("Dataset", "Id", "ImgType", "Slice")) %>% arrange(Time)
  
  if (!all(file.exists(paste("data",sliceImgMetaData$FileName,sep="/")))) { 
    # just in case not all files are on this filesystem
    next
  }
  cat("Processing",nextSlice,"/",nrow(sliceList),":",slice$Dataset,slice$Id,"slice",slice$Slice,
      paste("(",slice$SliceIndex, "/", slice$SliceCount, " order:", slice$SliceOrder, ")", sep=""),
      ifelse(any(sliceImgMetaData$isProcessed), "Re-segmenting", "Segmenting"),
      fill=T)

  # Calculate ROI either from all slices of the same Id with lower order (closer to the middle), or, if available and
  # returning a high probability, the area of the predicted LV segments in the previous slices.
  previousImages <- filter(left_join(filter(allSegments, Id==slice$Id, Slice!=slice$Slice), 
                                  #Id, Slice, starts_with("ROI.")),
                           imageList, by=c("Id", "Slice", "Time")), 
                 SliceOrder < slice$SliceOrder, # closer to middle
                 ((SliceOrder == 1) | (SliceOrder%%2 == slice$SliceOrder%%2))) # same side of middle; maybe only last few
  roi <- NULL
  roiPrevLVs <- NULL
  if (nrow(previousImages)>0) {
    roi <- list(x = mean(previousImages$ROI.x), y = mean(previousImages$ROI.y), r = mean(previousImages$ROI.r))
    if (!is.null(leftVentricleSegmentModel)) {
      predictors <- createSegmentPredictSet(select(previousImages, -isProcessed))
      previousImages$pLV <- round(predict(leftVentricleSegmentModel, data.matrix(predictors), missing=NaN),3)
      previousImages[,isBestGuess := pLV == max(pLV), by="Time"]
      previousImages <- filter(previousImages, isBestGuess, pLV > thresholdHighLV) # high threshold
      if (nrow(previousImages) > 5) {
        roiPrevLVs <- list(x = median(previousImages$m.cx), 
                           y = median(previousImages$m.cy), 
                           r = max(previousImages$s.radius.max))
      }
    }
  }
  
  sliceSegmentationInfo <- segmentImagesForOneSlice(sliceImgMetaData, roi, roiPrevLVs)
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
  if ((nextSlice %% 100 == 0) | (nextSlice == nrow(sliceList))) {
    # temp add Dataset to the segments for reporting and progress info
    allSegments <- left_join(allSegments, unique(select(imageList, Id, Slice, Time, Dataset)), by=c("Id","Slice","Time"))
    for (ds in unique(sliceList$Dataset)) {
      cat("...write", getSegmentFile(ds), "id's:",length(unique(filter(allSegments, Dataset==ds)$Id)), fill=T)
      write.csv(select(filter(allSegments, Dataset==ds), -Dataset), getSegmentFile(ds), row.names=F)
    }
    print("...completeness:")
    print(group_by(left_join(sliceList, 
                             mutate(unique(select(allSegments, Id, Slice)), isSegmented=TRUE), 
                             by=c("Id", "Slice")), Dataset) %>% 
            dplyr::summarise(complete = paste(round(100*sum(isSegmented, na.rm=T)/n(),2),"%",sep="" )))
    # Plot progress
    ds <- group_by(left_join(sliceList, 
                             mutate(unique(select(allSegments, Id, Slice)), isSegmented=TRUE), 
                             by=c("Id", "Slice")), Dataset, Id) %>% 
      dplyr::summarise(complete=round(100*sum(isSegmented, na.rm=T)/n())) %>%
      group_by(Dataset, complete) %>% 
      dplyr::summarise(n = n())
    print(ggplot(ds, aes(x=complete, y=n, fill=Dataset))+geom_bar(stat="identity")+
            ggtitle("Completeness of segmentation per Id"))
    # drop Dataset again - was only temp
    allSegments <- select(allSegments, -Dataset)
  }
}

# write final results (again, just to be sure)
allSegments <- left_join(allSegments, unique(select(imageList, Id, Slice, Time, Dataset)), by=c("Id","Slice","Time"))
for (ds in unique(imageList$Dataset)) {
  cat("Write final", getSegmentFile(ds), "id's:",length(unique(filter(allSegments, Dataset==ds)$Id)), fill=T)
  write.csv(select(filter(allSegments, Dataset==ds), -Dataset), getSegmentFile(ds), row.names=F)
}

print("Done!")



