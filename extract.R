# Extract all potentially relevant info from segments from all images
# and write these with all meta info into "segments-<dataset>.csv"
# for later consumption by prediction code.

# TODO: pick up images first that have not yet been processed and append
# (unless names have changed)

source("util.R")

# after systole = min
# after diastole = max
#  ejection fraction = 100 * (Vd - Vs)/Vd

# model with 600 classes??
# caret can do multi-class see http://stackoverflow.com/questions/15585501/usage-of-caret-with-gbm-method-for-multiclass-classification
# probs are cumulative?
# see https://www.kaggle.com/c/second-annual-data-science-bowl/details/evaluation

# see bio image detection stuff
# http://bioconductor.wustl.edu/bioc/vignettes/EBImage/inst/doc/AnalysisWithEBImage.pdf

listSliceImages <- function(sliceInfo) {
  imgFolder <- getImageFolder(sliceInfo)
  imgFiles <- list.files(imgFolder, pattern=".*\\.dcm$")
  allImgMetaData <- data.frame(Id = sliceInfo$Id, 
                               Dataset = sliceInfo$Dataset,
                               ImgType = sliceInfo$ImgType,
                               Slice = sliceInfo$Slice,
                               Offset = as.integer(gsub("^IM-(\\d{4,})-(\\d{4})\\.dcm$", "\\1", imgFiles)),
                               Time = as.integer(gsub("^IM-(\\d{4,})-(\\d{4})\\.dcm$", "\\2", imgFiles)),
                               pixelSpacing.x = rep(NA, length(imgFiles)),
                               pixelSpacing.y = rep(NA, length(imgFiles)),
                               sliceLocation = rep(NA, length(imgFiles)),
                               sliceThickness = rep(NA, length(imgFiles)),
                               stringsAsFactors = T)
  
  if (any(is.na(allImgMetaData$Offset)) | any(is.na(allImgMetaData$Time))) {
    #     print(allImgMetaData)
    #     print(imgFiles)
    print("WARN: Unexpected formats in image file")
    return(NULL)
  }
  if (nrow(allImgMetaData) == 0) {
    print("WARN: No image files")
    return(NULL)
  }
  # any na's for Time/Offset -> return NULL and flag unexpected file format
  
  # read first image and set meta info - assume it's the same for all image files
  dicomHeader <- readDICOMFile(getImageFile(allImgMetaData[1,]), pixelData = FALSE)
  pixelSpacing <- extractHeader(dicomHeader$hdr, "PixelSpacing", numeric=FALSE)
  
  allImgMetaData$pixelSpacing.x <- as.numeric(gsub("^(.*) (.*)$", "\\1", pixelSpacing))
  allImgMetaData$pixelSpacing.y <- as.numeric(gsub("^(.*) (.*)$", "\\2", pixelSpacing))
  allImgMetaData$sliceLocation <- extractHeader(dicomHeader$hdr, "SliceLocation", numeric=TRUE)
  allImgMetaData$sliceThickness <- extractHeader(dicomHeader$hdr, "SliceThickness", numeric=TRUE)
  
  return(allImgMetaData)
}

segmentImagesForOneSlice <- function(imgMetaData) {
  allImages <- vector(mode = "list", length = nrow(imgMetaData))
  for (i in seq(nrow(imgMetaData))) {
    f <- getImageFile(imgMetaData[i,])
    if (!file.exists(f)) {
      print(imgMetaData[i,])
      print(f)
      stop("File doesnt exist - formatting issues?")
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
    img_original <- allImages[[i]] 
    img_subtracted <- normalize(abs(img_original - avgImg))
    img_roi_thresholded <- (img_subtracted > otsu(img_subtracted))
    img_roi <- medianFilter(img_roi_thresholded, 2)
    
    # Find the ROI
    roi = as.data.frame(computeFeatures.moment(img_roi))
    coords <- as.data.frame(pos2coord(pos=which(img_roi>0),dim.mat=dim(img_roi)))
    names(coords) <- c('x','y')
    coords$distToROI <- floor(sqrt((coords$x - roi$m.cx)^2+(coords$y - roi$m.cy)^2))
    roi$radius <- quantile(coords$distToROI, probs=c(0.95))
    if (roi$radius > 1) {
      #     print(qplot(coords$distToROI, stat = "ecdf", geom = "step")+
      #             geom_vline(xintercept=roi$radius))
      
      roi_mask <- matrix(0, nrow=dim(img_roi)[1], ncol=dim(img_roi)[2])
      roi_mask <- drawCircle(roi_mask, x=roi$m.cx, y=roi$m.cy, roi$radius, col=1, fill=T)
      img_masked <- img_original*roi_mask
      scale <- 1/max(img_masked)
      img_original <- img_original*scale # scale original image wrt ROI
      
      kern <- makeBrush(5, shape= 'disc')
      #img_masked <- opening(img_original*roi_mask, kern)
      
      # Segmentation of the original image
      #img_thresholded <- (img_masked > 2*otsu(img_masked)) # bit iffy here
      #img_eroded <- opening(img_thresholded, kern)
      #img_segmented <- fillHull(bwlabel(img_thresholded))
      
      img_blurred <- normalize(opening(img_original, kern))
      img_thresholded <- img_blurred > otsu(img_blurred) # Otsu€™s threshold 
      img_segmented <- fillHull(bwlabel(img_thresholded))
      
      img_comb <- EBImage::combine(drawCircle(colorLabels(img_segmented), x=roi$m.cx, y=roi$m.cy, roi$radius, "yellow", fill=FALSE, z=1),
                                   drawCircle(toRGB(img_original), x=roi$m.cx, y=roi$m.cy, roi$radius, "yellow", fill=FALSE, z=1),
                                   toRGB(img_subtracted),
                                   #toRGB(img_masked),
                                   toRGB(img_blurred),
                                   toRGB(img_thresholded),
                                   toRGB(img_roi))
      display(img_comb,all=T,method="raster")
      text(10,20,getImageFile(imgMetaData[i,]),col="yellow",pos=4)
      
      # Get segment meta data and select only the segments within the ROI
      segmentInfo <- cbind(imgMetaData[i,],
                           data.frame(computeFeatures.moment(img_segmented)),
                           data.frame(computeFeatures.shape(img_segmented)))
      #segmentInfo$roundness <- 4*pi*segmentInfo$s.area/(segmentInfo$s.perimeter^2)
      segmentInfo$distToROI <- sqrt((segmentInfo$m.cx - roi$m.cx)^2 + (segmentInfo$m.cy - roi$m.cy)^2)
      #     segmentInfo$imgIndex <- i # same as Time
      segmentInfo$segIndex <- 1:nrow(segmentInfo)
      segmentInfo$UUID <- sapply(1:nrow(segmentInfo),UUIDgenerate) # unique ID for each segment
      
      # Only keep segments inside the ROI
      segmentInfo <- filter(segmentInfo, distToROI < roi$radius)
      
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

# print(segmentImagesForOneSlice(filter(imgMetaData, Slice == sort(unique(imgMetaData$Slice))[2])))
# stop()

# processSegmentedImagesForOneSlice <- function(segmentedImages)
# {
#   allSegments <- NULL
#   for (img_original in allImages) { 
#     segments <- cbind(as.data.frame(computeFeatures.shape(img_segmented)), 
#                       as.data.frame(computeFeatures.moment(img_segmented)))
#     segments$roundness <- 4*pi*segments$s.area/(segments$s.perimeter^2)
#     segments$segmentNumber <- 1:nrow(segments)
#     
#     if (is.null(allSegments)) {
#       allSegments <- segments
#       maxSegmentNumber <- max(allSegments$segmentNumber)
#     } else {
#       for (i in 1:nrow(segments)) {
#         # find previous segments at approx the same position
#         distance <- sqrt((segments$m.cx[i] - allSegments$m.cx)^2 + (segments$m.cy[i] - allSegments$m.cy)^2)
#         # TODO eliminate doubles (segNr already in segments)
#         a <- which.min(distance)
#         if (distance[a] < 10) {
#           segments$segmentNumber[i] <- allSegments$segmentNumber[a]
#         } else {
#           maxSegmentNumber <- maxSegmentNumber+1
#           segments$segmentNumber[i] <- maxSegmentNumber
#         }
#       }
#       allSegments <- rbind(allSegments, segments)
#       for (i in 1:nrow(segments)) {
#         text(x = segments$m.cx[i], y = segments$m.cy[i], 
#              label = segments$segmentNumber[i], col = "orange")
#       }
#     } # segments of one image
#     allSegments$segmentNumber <- as.integer(allSegments$segmentNumber)
#   }
# }

playlist <- NULL
for (dataset in datasetFolders) {
  datasetDir <- paste("data",dataset,sep="/")
  if (dir.exists(datasetDir)) {
    Ids <- sort(as.integer(list.dirs(datasetDir, recursive=F, full.names=F)))
    for (Id in Ids) {
      IdFolder <- paste(datasetDir, Id, "study", sep="/")
      subFolders <- list.dirs(IdFolder, recursive=F, full.names=F)
      # NB excluding the 2 chamber and 4 chamber views here, only the short axis track
      saxFolders <- subFolders[ grepl("sax_[[:digit:]]+$", subFolders) ]
      slices <- sort(as.integer(gsub("sax_([[:digit:]]+)$", "\\1", saxFolders)))
      imgMetaData <- data.frame(Id = Id, 
                                Dataset = dataset,
                                ImgType = "sax",
                                Slice = slices,
                                SliceCount = length(slices),
                                SliceRelIndex = seq(length(slices)),
                                SliceOrder = midOrderSeqKeepSibblingsTogether(length(slices)),
                                stringsAsFactors = F)
      if (is.null(playlist)) {
        playlist <- imgMetaData
      } else {
        playlist <- rbind(playlist, imgMetaData)
      }
    }
  }
}
# playlist <- left_join(playlist, group_by(playlist, Id) %>% summarise(nSlices=n()))

# Show quick summary of the datasets & save for downstream use
print(select(playlist, Dataset, Id) %>% unique() %>% group_by(Dataset) %>% summarise(nIds = n()))
write.csv(playlist, "slicelist.csv", row.names=F)

allSegmentationInfo <- NULL
for (dataset in datasetFolders) {
  if (file.exists(getSegmentFile(dataset))) {
    segmentsPerDataset <- fread(getSegmentFile(dataset))
    if (is.null(allSegmentationInfo)) {
      allSegmentationInfo <- segmentsPerDataset
    } else {
      removedSet <- setdiff(names(allSegmentationInfo), names(segmentsPerDataset))
      addedSet <- setdiff(names(segmentsPerDataset), names(allSegmentationInfo))
      diffSet <- paste(c(paste("-",removedSet), paste("+",addedSet)),collapse=", ")
      if (length(removedSet) + length(addedSet) > 0) {
        print(names(allSegmentationInfo))
        print(names(segmentsPerDataset))
        print(diffSet)
        stop("Datasets do not match up. Please consider removing files.")
      }
      allSegmentationInfo <- rbind(allSegmentationInfo, segmentsPerDataset)
    }
  }
}

# Figure out which slices have been processed and have one of their two sibblings processed also (these are complete).
# This drives the order of processing of the slices

po <- left_join(select(playlist, Id, Slice), 
                select(allSegmentationInfo, Id, Slice, Dataset), by = c("Id", "Slice")) %>% 
  unique() %>%
  mutate(isProcessed = !is.na(Dataset)) %>% 
  select(-Dataset) %>%
  left_join(playlist) %>%
  group_by(Id) %>%
  mutate(sib1 = SliceRelIndex - 1,
         sib2 = SliceRelIndex + 1)
po <- as.data.frame(po)
playlist <- left_join(left_join(po, 
                                left_join( select(po, Id, sib1), 
                                           select(po, Id, SliceRelIndex, Slice, isProcessed), by = c("Id", "sib1"="SliceRelIndex")),
                                by=c("Id","sib1")),
                      left_join( select(po, Id, sib2), 
                                 select(po, Id, SliceRelIndex, Slice, isProcessed), by = c("Id", "sib2"="SliceRelIndex")),
                      by=c("Id","sib2")) %>%
  select(-sib1, -sib2, -Slice, -Slice.y) %>%
  rename(Slice = Slice.x, sib1Processed = isProcessed.y, sib2Processed = isProcessed, isProcessed = isProcessed.x) %>%
  mutate(sliceComplete = isProcessed & (sib1Processed | sib2Processed)) %>%
  group_by(Id) %>% mutate(anySliceCompletePerId = any(sliceComplete)) %>% ungroup()

# Here we order the playlist:
# - first incomplete ID's, then unprocessed slices by their special 'mid order'
playlist <- arrange(playlist, anySliceCompletePerId, isProcessed, SliceOrder, Dataset)

for (nextSlice in seq(nrow(playlist))) {
  cat("Processing",playlist$Dataset[nextSlice],
      playlist$Id[nextSlice],"slice",
      playlist$Slice[nextSlice],
      paste("(",playlist$SliceRelIndex[nextSlice], "/", playlist$SliceCount[nextSlice], " order:", playlist$SliceOrder[nextSlice], ")", sep=""),
      fill=T)
  
  # Add slice meta data and all images
  imgMetaData <- listSliceImages(select(playlist[nextSlice,], Id, Dataset, ImgType, Slice, SliceRelIndex))
  
  # NB slice thickness etc to be derived in predict, not here
  
  if (!is.null(imgMetaData)) {
    # Segment all images in this slice
    sliceSegmentationInfo <- segmentImagesForOneSlice(imgMetaData)
    cat("...",nrow(sliceSegmentationInfo),"segments in",nrow(imgMetaData),"images",fill=T)
    
    if (is.null(allSegmentationInfo)) {
      allSegmentationInfo <- sliceSegmentationInfo
    } else {
      removedSet <- setdiff(names(allSegmentationInfo), names(sliceSegmentationInfo))
      addedSet <- setdiff(names(sliceSegmentationInfo), names(allSegmentationInfo))
      diffSet <- paste(c(paste("-",removedSet), paste("+",addedSet)),collapse=", ")
      if (length(removedSet) + length(addedSet) > 0) {
        print("Existing:")
        print(names(allSegmentationInfo))
        print("New:")
        print(names(sliceSegmentationInfo))
        print(paste(diffSet, collapse=","))
        stop("Datasets do not match up. Please consider removing segment files.")
      }
      
      nDrop <- nrow(filter(allSegmentationInfo, 
                           Slice %in% sliceSegmentationInfo$Slice,
                           Id %in% sliceSegmentationInfo$Id))
      if (nDrop > 0) {
        cat("...dropping",nDrop,"rows from previous results",fill=T)
      }
      
      allSegmentationInfo <- rbind(filter(allSegmentationInfo, 
                                          !(Slice %in% sliceSegmentationInfo$Slice & 
                                              Id %in% sliceSegmentationInfo$Id)), 
                                   sliceSegmentationInfo)
    }
    
    # Add segmentation to results and write file (once in a while)
    if ((nextSlice %% 10 == 0) |
          ((nextSlice < nrow(playlist)) & (playlist$Dataset[nextSlice+1] != playlist$Dataset[nextSlice]))) {
      ds <- playlist$Dataset[nextSlice]
      cat("...write", getSegmentFile(ds), "id's:",length(unique(filter(allSegmentationInfo, Dataset==ds)$Id)), fill=T)
      write.csv(filter(allSegmentationInfo, Dataset==ds), getSegmentFile(ds), row.names=F)
      
      print("...completeness:")
      print(group_by(left_join(playlist, 
                               mutate(unique(select(allSegmentationInfo, Id, Slice)), isSegmented=TRUE), 
                               by=c("Id", "Slice")), Dataset) %>% 
              summarise(percentComplete = paste(round(100*sum(isSegmented, na.rm=T)/n(),2),"%",sep="" )))
      
      # Show progress
      ds <- group_by(left_join(playlist, 
                               mutate(unique(select(allSegmentationInfo, Id, Slice)), isSegmented=TRUE), 
                               by=c("Id", "Slice")), Dataset, Id) %>% 
        summarise(complete=round(100*sum(isSegmented, na.rm=T)/n())) %>%
        group_by(Dataset, complete) %>% summarise(n = n())
      print(ggplot(ds, aes(x=complete, y=n, fill=Dataset))+geom_bar(stat="identity")+ggtitle("Completeness of segmentation"))
    }
  } else {
    print("...no valid image data!")
  }
}

# write final results (again, just to be sure)
for (ds in datasetFolders) {
  cat("Write final", getSegmentFile(ds), "id's:",length(unique(filter(allSegmentationInfo, Dataset==ds)$Id)), fill=T)
  write.csv(filter(allSegmentationInfo, Dataset==ds), getSegmentFile(ds), row.names=F)
}

print("Segmented:")
print(group_by(playlist, Dataset) %>% 
        summarise(percentComplete = paste(round(100*sum(sliceComplete, na.rm=T)/n(),2),"%",sep="" )))





