# Extract all potentially relevant info from segments from all images
# and write these with all meta info into "segments-<dataset>.csv"
# for later consumption by prediction code.

# TODO: pick up images first that have not yet been processed and append
# (unless names have changed)

source("util.R")

datasetFolders <- c('train','validate','test')

# create mid-ordered sequence
midOrderSeq <- function(n) {
  abs(n %/% 2 + 1 - seq(n))*2 - ((1+sign(n %/% 2 + 1 - seq(n)))%/%2) + 1
}

# Extension of above so that there are always pairs of slices present
midOrderSeqKeepSibblingsTogether <- function(n) {
  s <- as.integer(factor(abs(n %/% 2 + 1 - 2*(seq(n) %/% 2))*2 + 1))
  if (length(s) > 1) {
    s[1] <- s[2]
    s[length(s)] <- s[length(s)-1]
  }
  return(s)
}

segmentFile <- function(ds) {
  paste("segments-",ds,".csv",sep="")
}

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
  imgFolder <- paste("data",sliceInfo$Dataset,sliceInfo$Id,"study",
                     paste(sliceInfo$ImgType, sliceInfo$Slice, sep="_"),sep="/")
  imgFiles <- list.files(imgFolder, pattern=".*\\.dcm$")
  isFirstImage <- T
  allImgMetaData <- NULL
  for (imgFile in imgFiles) {
    imgOffset <- as.integer(gsub("^IM-(\\d{4,})-(\\d{4})\\.dcm$", "\\1", imgFile))
    imgTime <- as.integer(gsub("^IM-(\\d{4,})-(\\d{4})\\.dcm$", "\\2", imgFile))
    fname <- paste(imgFolder, imgFile, sep="/")
    if (isFirstImage) {
      dicomHeader <- readDICOMFile(fname, pixelData = FALSE)
      pixelSpacing <- extractHeader(dicomHeader$hdr, "PixelSpacing", numeric=FALSE)
      isFirstImage <- F
    }
    imgMetaData <- data.frame(Id = sliceInfo$Id, 
                              Dataset = sliceInfo$Dataset,
                              ImgType = sliceInfo$ImgType,
                              Slice = sliceInfo$Slice,
                              SliceRelIndex = sliceInfo$SliceRelIndex,
                              Offset = imgOffset,
                              Time = imgTime,
                              pixelSpacing.x = as.numeric(gsub("^(.*) (.*)$", "\\1", pixelSpacing)),
                              pixelSpacing.y = as.numeric(gsub("^(.*) (.*)$", "\\2", pixelSpacing)),
                              sliceLocation = extractHeader(dicomHeader$hdr, "SliceLocation", numeric=TRUE),
                              sliceThickness = extractHeader(dicomHeader$hdr, "SliceThickness", numeric=TRUE),
                              file = fname,
                              stringsAsFactors = F)
    if (is.null(allImgMetaData)) {
      allImgMetaData <- imgMetaData
    } else {
      allImgMetaData <- rbind(allImgMetaData, imgMetaData)
    }
  }
  
  return(allImgMetaData)
}

segmentImagesForOneSlice <- function(imgMetaData) {
  allImages <- vector(mode = "list", length = nrow(imgMetaData))
  for (i in 1:nrow(imgMetaData)) {
    f <- imgMetaData$file[i]
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
    
    # Get segment meta data and select only the segments within the ROI
    segmentInfo <- cbind(imgMetaData[i,],
                         data.frame(computeFeatures.moment(img_segmented)),
                         data.frame(computeFeatures.shape(img_segmented)))
    segmentInfo$roundness <- 4*pi*segmentInfo$s.area/(segmentInfo$s.perimeter^2)
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
  } # all images of one slice
  
  return(select(sliceSegmentation, -file))
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
playlist <- left_join(playlist, group_by(playlist, Id) %>% summarise(nSlices=n()))

# Show quick summary of the datasets
print( group_by(playlist, Dataset, Id) %>% summarise(nSlices = n()) )
print( group_by(playlist, Dataset) %>% summarise(nCases = n()) )

previouslyProcessedIds <- NULL
allSegmentationInfo <- NULL
for (dataset in datasetFolders) {
  if (file.exists(segmentFile(dataset))) {
    segmentsPerDataset <- fread(segmentFile(dataset))

    # Only consider an ID processed if each slice has also it's sibling slice processed. This means
    # we may re-process a full ID instead of just a couple of slices.
    # If there is only 1 slice, it's ok too
    hasSiblings <- group_by(segmentsPerDataset, Id) %>% 
      select(Id, Slice, SliceRelIndex) %>% arrange(Id) %>%
      unique() %>%
      summarise(nProcessedSlices = n(),
                allSlicesHaveSibblings = all(((SliceRelIndex - 1) %in% SliceRelIndex) | ((SliceRelIndex + 1) %in% SliceRelIndex)))
    hasSiblings <- left_join(hasSiblings, unique(select(playlist, Id, nSlices)), copy=T)
    hasSiblings$ConsiderProcessed <- (hasSiblings$nSlices < 2) | (hasSiblings$allSlicesHaveSibblings)

    if (is.null(previouslyProcessedIds)) {
      previouslyProcessedIds <- sort(unique(hasSiblings$Id[hasSiblings$ConsiderProcessed]))
    } else {
      previouslyProcessedIds <- sort(unique(c(previouslyProcessedIds,hasSiblings$Id[hasSiblings$ConsiderProcessed])))
    }
    
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

playlist$isProcessed <- playlist$Id %in% previouslyProcessedIds
playlist <- arrange(playlist, isProcessed, SliceOrder, Dataset)

for (nextSlice in seq(nrow(playlist))) {
  cat("Processing",playlist$Dataset[nextSlice],
      playlist$Id[nextSlice],"slice",
      playlist$Slice[nextSlice],fill=T)
  
  # Add slice meta data and all images
  imgMetaData <- listSliceImages(select(playlist[nextSlice,], Id, Dataset, ImgType, Slice, SliceRelIndex))
  
  # NB slice thickness etc to be derived in predict, not here
  
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
  
  # Add segmentation to results and write file (or write once in a while)
  for (ds in unique(allSegmentationInfo$Dataset)) {
    cat("...write", segmentFile(ds), "id's:",length(unique(filter(allSegmentationInfo, Dataset==ds)$Id)), fill=T)
    
    write.csv(filter(allSegmentationInfo, Dataset==ds), segmentFile(ds), row.names=F)
  }
}


# for (dataset in datasetFolders) {
#   datasetDir <- paste("data",dataset,sep="/")
#   if (dir.exists(datasetDir)) {
#     Ids <- as.integer(list.dirs(datasetDir, recursive=F, full.names=F))
#     allSegInfo <- NULL
#     for (Id in Ids) {
#       IdFolder <- paste(datasetDir, Id, "study", sep="/")
#       imgMetaData <- processSingleId(Id, IdFolder)
#       
#       # figure out distance and area multiplier
#       # assuming there are 2 slices at least, also assuming constant distance and area
#       # but see named exceptions ID 123, 148 ... (TODO deal with those - or ignore)
#       sliceData <- group_by(filter(imgMetaData, Time == min(imgMetaData$Time)), Slice)
#       sliceDistance <- abs(sliceData$sliceLocation[2] - sliceData$sliceLocation[1])
#       if (is.na(sliceDistance)) { # data not available in all image meta data
#         sliceDistance <- sliceData$sliceThickness[1]
#       }
#       if (is.na(sliceDistance)) {
#         sliceDistance <- 8 # fall back scenario
#       }
#       sliceAreaMultiplier <- sliceData$pixelSpacing.x[1]*sliceData$pixelSpacing.y[1]
#       # TODO: slide area multiplier not always available - eg validation case 516
#       # fall back to average of others or to 2 hardcoded
#       cat("Slice distance: ", sliceDistance, fill=T)
#       cat("Slice area multiplier: ", sliceAreaMultiplier, fill=T)
#       
#       segInfo <- processImagesForOneId(Id, imgMetaData, sliceDistance, sliceAreaMultiplier)
#       segInfo$sliceArea <- segInfo$s.area*sliceAreaMultiplier
#       segInfo$sliceVolume <- segInfo$sliceArea*sliceDistance
#       
#       if (is.null(allSegInfo)) {
#         allSegInfo <- segInfo
#       } else {
#         allSegInfo <- rbind(allSegInfo, segInfo)
#       }
#       write.csv(allSegInfo, segmentFile(dataset), row.names=F)
#     }
#   }
# }
