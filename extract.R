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

processSingleId <- function(Id, folder) {
  cat("Processing meta info for case", Id, fill=T)
  allImgMetaData <- NULL
  subFolders <- list.dirs(folder, recursive=F, full.names=F)
  # NB excluding the 2 chamber and 4 chamber views here, only the short axis track
  print(folder)
  saxFolders <- subFolders[ grepl("sax_[[:digit:]]+$", subFolders) ]
  slices <- sort(as.integer(gsub("sax_([[:digit:]]+)$", "\\1", saxFolders)))
  for (slice in slices) { 
    imgFolder <- paste(folder, "/", "sax_", slice, sep="")
    imgFiles <- list.files(imgFolder, pattern=".*\\.dcm$")
    for (imgFile in imgFiles) {
      imgOffset <- as.integer(gsub("^IM-(\\d{4,})-(\\d{4})\\.dcm$", "\\1", imgFile))
      imgTime <- as.integer(gsub("^IM-(\\d{4,})-(\\d{4})\\.dcm$", "\\2", imgFile))
      fname <- paste(imgFolder, imgFile, sep="/")
      dicomHeader <- readDICOMFile(fname, pixelData = FALSE)
      pixelSpacing <- extractHeader(dicomHeader$hdr, "PixelSpacing", numeric=FALSE)
      imgMetaData <- data.frame(Id = Id, 
                                ImgType = "sax",
                                Slice = slice,
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
    
    #roi_mask <- matrix(0, nrow=dim(img_roi)[1], ncol=dim(img_roi)[2])
    #roi_mask <- drawCircle(roi_mask, x=roi$m.cx, y=roi$m.cy, roi$radius, col=1, fill=T)
    
    kern <- makeBrush(5, shape= 'disc')
    #img_masked <- opening(img_original*roi_mask, kern)
    
    # Segmentation of the original image
    #img_thresholded <- (img_masked > 2*otsu(img_masked)) # bit iffy here
    #img_eroded <- opening(img_thresholded, kern)
    #img_segmented <- fillHull(bwlabel(img_thresholded))
    img_blurred <- normalize(opening(img_original, kern))
    img_thresholded <- img_blurred > otsu(img_blurred) # Otsu’s threshold 
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
  
  return(sliceSegmentation)
}

# print(segmentImagesForOneSlice(filter(imgMetaData, Slice == sort(unique(imgMetaData$Slice))[2])))
# stop()

processSegmentedImagesForOneSlice <- function(segmentedImages)
{
  allSegments <- NULL
  for (img_original in allImages) { 
    segments <- cbind(as.data.frame(computeFeatures.shape(img_segmented)), 
                      as.data.frame(computeFeatures.moment(img_segmented)))
    segments$roundness <- 4*pi*segments$s.area/(segments$s.perimeter^2)
    segments$segmentNumber <- 1:nrow(segments)
    
    if (is.null(allSegments)) {
      allSegments <- segments
      maxSegmentNumber <- max(allSegments$segmentNumber)
    } else {
      for (i in 1:nrow(segments)) {
        # find previous segments at approx the same position
        distance <- sqrt((segments$m.cx[i] - allSegments$m.cx)^2 + (segments$m.cy[i] - allSegments$m.cy)^2)
        # TODO eliminate doubles (segNr already in segments)
        a <- which.min(distance)
        if (distance[a] < 10) {
          segments$segmentNumber[i] <- allSegments$segmentNumber[a]
        } else {
          maxSegmentNumber <- maxSegmentNumber+1
          segments$segmentNumber[i] <- maxSegmentNumber
        }
      }
      allSegments <- rbind(allSegments, segments)
      for (i in 1:nrow(segments)) {
        text(x = segments$m.cx[i], y = segments$m.cy[i], 
             label = segments$segmentNumber[i], col = "orange")
      }
    } # segments of one image
    allSegments$segmentNumber <- as.integer(allSegments$segmentNumber)
  }
}

processImagesForOneId <- function(Id, imgMetaData, sliceDistance, sliceAreaMultiplier) {
  allSegmentationInfo <- NULL
  for (slice in unique(imgMetaData$Slice)) {
    cat("Processing images for case", Id, "slice",slice,fill=T)
    sliceSegmentationInfo <- segmentImagesForOneSlice(filter(imgMetaData, Slice == slice))
    
    #print(sliceSegmentationInfo)
    if (is.null(allSegmentationInfo)) {
      allSegmentationInfo <- sliceSegmentationInfo
    } else {
      allSegmentationInfo <- rbind(allSegmentationInfo, sliceSegmentationInfo)
    }
  }
  return(allSegmentationInfo)
}

for (dataset in c('train','validate','test')) {
  datasetDir <- paste("data",dataset,sep="/")
  if (dir.exists(datasetDir)) {
    Ids <- as.integer(list.dirs(datasetDir, recursive=F, full.names=F))
    allSegInfo <- NULL
    for (Id in Ids) {
      IdFolder <- paste(datasetDir, Id, "study", sep="/")
      imgMetaData <- processSingleId(Id, IdFolder)
      
      # figure out distance and area multiplier
      # assuming there are 2 slices at least, also assuming constant distance and area
      # but see named exceptions ID 123, 148 ... (TODO deal with those - or ignore)
      sliceData <- group_by(filter(imgMetaData, Time == min(imgMetaData$Time)), Slice)
      sliceDistance <- abs(sliceData$sliceLocation[2] - sliceData$sliceLocation[1])
      if (is.na(sliceDistance)) { # data not available in all image meta data
        sliceDistance <- sliceData$sliceThickness[1]
      }
      if (is.na(sliceDistance)) {
        sliceDistance <- 8 # fall back scenario
      }
      sliceAreaMultiplier <- sliceData$pixelSpacing.x[1]*sliceData$pixelSpacing.y[1]
      # TODO: slide area multiplier not always available - eg validation case 516
      # fall back to average of others or to 2 hardcoded
      cat("Slice distance: ", sliceDistance, fill=T)
      cat("Slice area multiplier: ", sliceAreaMultiplier, fill=T)
      
      segInfo <- processImagesForOneId(Id, imgMetaData, sliceDistance, sliceAreaMultiplier)
      segInfo$sliceArea <- segInfo$s.area*sliceAreaMultiplier
      segInfo$sliceVolume <- segInfo$sliceArea*sliceDistance
      
      if (is.null(allSegInfo)) {
        allSegInfo <- segInfo
      } else {
        allSegInfo <- rbind(allSegInfo, segInfo)
      }
      write.csv(allSegInfo, paste("segments-",dataset,".csv",sep=""), row.names=F)
    }
  }
}

# imageData <- NULL
# 
# for (dataset in c('train','validate','test')) {
#   datasetDir <- paste("data",dataset,sep="/")
#   if (dir.exists(datasetDir)) {
#     IdFolders <- list.dirs(datasetDir, recursive=F, full.names=F)
#     for (Id in sample(IdFolders, length(IdFolders))) {
#       #   print(IdFolder) # = nr 1 .. 500
#       IdFolder <- paste(datasetDir, Id, "study", sep="/")
#       
#       serieFolders <- list.dirs(IdFolder, recursive=F, full.names=F)
#       # NB excluding the 2 chamber and 4 chamber views here, only the short axis track
#       serieFolders <- serieFolders[ grepl("sax_[[:digit:]]+$", serieFolders) ]
#       
#       for (serieFolder in sample(serieFolders, length(serieFolders))) { 
#         # serieFolder <- 'data/train/500/study/sax_22'
#         imgFolder <- paste(IdFolder, serieFolder, sep="/")
#         
#         print("Processing:")
#         print(imgFolder)
#         
#         seriesImageFiles <- list.files(imgFolder, pattern=".*dcm$")
#         
#         allSegments <- NULL
#         for (imgNumber in 1:length(seriesImageFiles)) {
#           f <- seriesImageFiles[imgNumber]
#           fname <- paste(imgFolder, f, sep="/")
#           
#           print(fname)
#           
#           dicom <- readDICOMFile(fname)
#           
#           
#           # read meta info
#           # see _read_all_dicom_images
#           # in https://www.kaggle.com/c/second-annual-data-science-bowl/details/fourier-based-tutorial
#           #           x <- readDICOMFile(fname)
#           #           pixelSpacing <- extractHeader(x$hdr, "PixelSpacing", numeric=FALSE)
#           #           pSmat <- header2matrix(pixelSpacing, ncol=2)
#           #           IOP <- extractHeader(x$hdr, "ImageOrientationPatient", numeric=FALSE)
#           #           IOPmat <- header2matrix(IOP, ncol=6)
#           
#           
#           
#           
#           img <- Image(normalize(dicom$img))
#           img <- rotate(img,-90)
#           # display(img, method = "raster")
#           
#           img_thresholded <- img > otsu(img) # Otsu’s threshold 
#           img_denoised <- medianFilter(img_thresholded, 4)
#           img_segmented <- fillHull(bwlabel(img_denoised))
#           img_comb = EBImage::combine(
#             colorLabels(img_segmented), # TODO try preserve colours
#             toRGB(img),
#             toRGB(img_thresholded),
#             toRGB(img_denoised)
#           )
#           
#           display(img_comb, title=imgFolder, all=T,method="raster")
#           #           for (i in 1:max(img_segmented)) {
#           #             coords <- as.data.frame(pos2coord(pos=which(img_segmented==i),dim.mat=dim(img_segmented)))
#           #             names(coords) <- c('x','y')
#           #             text(x = mean(coords$x), y = mean(coords$y), label = i, col = "orange")
#           #           }
#           
#           # Calculate image features. EBImage does a good job at this, see doc
#           # https://www.bioconductor.org/packages/3.3/bioc/manuals/EBImage/man/EBImage.pdf
#           # shape factors https://en.wikipedia.org/wiki/Shape_factor_(image_analysis_and_microscopy)
#           # http://www.cs.uu.nl/docs/vakken/ibv/reader/chapter8.pdf
#           segments <- cbind(as.data.frame(computeFeatures.shape(img_segmented)), 
#                             as.data.frame(computeFeatures.moment(img_segmented)))
#           # TODO basic, haralick, ..?
#           segments <- dplyr::filter(segments, s.area > 10) # early filter to remove small patches
#           segments$roundness <- 4*pi*segments$s.area/(segments$s.perimeter^2)
#           #           segments$compactness <- 1/segments$roundness
#           segments$segmentNumber <- 1:nrow(segments)
#           segments$imgNumber <- imgNumber
#           
#           if (is.null(allSegments)) {
#             allSegments <- segments
#             maxSegmentNumber <- max(allSegments$segmentNumber)
#           } else {
#             for (i in 1:nrow(segments)) {
#               # find previous segments at approx the same position
#               distance <- sqrt((segments$m.cx[i] - allSegments$m.cx)^2 + (segments$m.cy[i] - allSegments$m.cy)^2)
#               #               candidates <- which(distance < 10) # etc.
#               a <- which.min(distance)
#               if (distance[a] < 10 ){
#                 segments$segmentNumber[i] <- allSegments$segmentNumber[a]
#               } else {
#                 maxSegmentNumber <- maxSegmentNumber+1
#                 segments$segmentNumber[i] <- maxSegmentNumber
#               }
#             }
#             allSegments <- rbind(allSegments, segments)
#             for (i in 1:nrow(segments)) {
#               text(x = segments$m.cx[i], y = segments$m.cy[i], 
#                    label = segments$segmentNumber[i], col = "orange")
#             }
#           } # segments of one image
#         } # images
#         allSegments$segmentNumber <- as.integer(allSegments$segmentNumber)
#         
#         # Now try figure out which image segment is the one of interest
#         # must be large and small, maybe also filter on size or roundness
#         #         largeSegments <- unique(filter(allSegments, s.area>400)$segmentNumber)
#         #         smallSegments <- filter(allSegments, s.area<300)
#         #         candidateSegments <- intersect(largeSegments, smallSegments)
#         #         candidateSegments <- (group_by(allSegments, segmentNumber) %>% 
#         #                                 summarise(size_variance = sd(s.area)) %>% arrange(desc(size_variance)) %>% ungroup() %>% 
#         #                                 mutate(size_variance_scaled=scale(size_variance)) %>% filter(size_variance_scaled > 1))$segmentNumber
#         #         candidateSegments <- unique(allSegments$segmentNumber)
#         
#         data_aggregatedAllCandidateSegments <- NULL
#         for (s in unique(allSegments$segmentNumber)) {
#           data_candidateSegment <- filter(allSegments, segmentNumber == s) %>% 
#             select(-segmentNumber, -imgNumber, -m.cx, -m.cy)
#           
#           data_mean <- sapply(data_candidateSegment, mean)
#           names(data_mean) <- paste(names(data_candidateSegment), "mean", sep="_")
#           
#           data_sd <- sapply(data_candidateSegment, sd)
#           names(data_sd) <- paste(names(data_candidateSegment), "sd", sep="_")
#           
#           data_quantiles <- sapply(data_candidateSegment, quantile, probs=c(0.1,0.5,0.9))
#           data_p10 <- data_quantiles[1,]
#           names(data_p10) <- paste(names(data_candidateSegment), "p10", sep="_")
#           data_p50 <- data_quantiles[2,]
#           names(data_p50) <- paste(names(data_candidateSegment), "p50", sep="_")
#           data_p90 <- data_quantiles[3,]
#           names(data_p90) <- paste(names(data_candidateSegment), "p90", sep="_")
#           
#           data_aggregatedCandidateSegment <- c(segmentNumber=s, 
#                                                data_mean, data_sd, data_p10, data_p50, data_p90,
#                                                img_present = nrow(data_candidateSegment))
#           
#           if (is.null(data_aggregatedAllCandidateSegments)) {
#             data_aggregatedAllCandidateSegments <- as.data.frame(t(data_aggregatedCandidateSegment))
#           } else {
#             # only happens when there are multiple candidates
#             data_aggregatedAllCandidateSegments <- rbind(data_aggregatedAllCandidateSegments,t(data_aggregatedCandidateSegment))
#           }
#         }
#         # Add meta-info about the segments
#         if (nrow(data_aggregatedAllCandidateSegments) >= 1) {
#           row.names(data_aggregatedAllCandidateSegments) <- 1:nrow(data_aggregatedAllCandidateSegments)
#           data_aggregatedAllCandidateSegments$dataset <- dataset
#           data_aggregatedAllCandidateSegments$Id <- as.integer(Id)
#           data_aggregatedAllCandidateSegments$series <- serieFolder
#           data_aggregatedAllCandidateSegments$series_idx <- which(serieFolder == serieFolders)
#           data_aggregatedAllCandidateSegments$series_cnt <- length(serieFolders)
#           data_aggregatedAllCandidateSegments$img_cnt <- length(seriesImageFiles)
#         }
#         
#         # Perhaps threshold prediciton at 0.50 - if so, then assume yes
#         # or, 10 x next one up
#         
#         # Make prediction of Left Ventricle segment using previous data
#         segmentsInBestOrder <- data_aggregatedAllCandidateSegments$segmentNumber
#         if (file.exists('imageData.csv')) {
#           segData <- fread('imageData.csv', 
#                            drop=c('segmentNumber', 'dataset', 'Id', 'series', 'series_idx', 'series_cnt'))
#           if ('isLeftVentricle' %in% names(segData)) {
#             trainData <- select(segData, -isLeftVentricle)
#             bst <- xgboost(data = as.matrix(trainData), 
#                            label = segData$isLeftVentricle, 
#                            max.depth = 8, eta = 0.1, nround = 100,
#                            objective = "binary:logistic", missing=NaN, verbose=0)
#             imp_matrix <- xgb.importance(feature_names = names(trainData), model = bst)
#             print(xgb.plot.importance(importance_matrix = imp_matrix))
#             newData <- select(data_aggregatedAllCandidateSegments, -segmentNumber, -dataset, -Id, -series, -series_idx, -series_cnt)
#             pred <- predict(bst, as.matrix(newData), missing=NaN)
#             orderedPreds <- arrange(data.frame(segment=data_aggregatedAllCandidateSegments$segmentNumber, 
#                                                prob=pred), desc(prob))
#             print(head(orderedPreds))
#             segmentsInBestOrder <- orderedPreds$segment
#           }
#         }
#         print(segmentsInBestOrder)
#         leftVentricleSegment <- readInteger("Enter the segment number of the left ventricle (0 = none):")
#         cat("Segment = ", leftVentricleSegment, fill=T)
#         
#         # TODO create a model using the aggregates (min/max/mean/sd etc) of:
#         # "s.area"         "s.perimeter"    "s.radius.mean"  "s.radius.sd"    "s.radius.min"   "s.radius.max"  
#         # "m.eccentricity" "roundness"
#         
#         
#         data_aggregatedAllCandidateSegments$isLeftVentricle <- 
#           data_aggregatedAllCandidateSegments$segmentNumber == leftVentricleSegment        
#         
#         if (F) { # nrow(data_aggregatedAllCandidateSegments) > 1
#           cat("*** Identified multiple segments:",candidateSegments,fill=T)
#           
#           data_scaled <- as.data.frame(cbind(scale(allSegments[,sapply(allSegments,class) == "numeric"]),
#                                              allSegments[,sapply(allSegments,class) != "numeric"]))
#           plotData <- filter(data_scaled, segmentNumber %in% candidateSegments) %>% 
#             select(m.eccentricity, roundness, s.area, segmentNumber, imgNumber) %>%
#             gather(Feature, Value, -segmentNumber, -imgNumber)
#           print(ggplot(data=plotData, aes(x=imgNumber, y=Value, linetype=Feature, color=factor(segmentNumber)))+geom_line()+
#                   ggtitle(paste("Segments in",imgFolder)))
#           
#           bestIndex <- which.min(data_aggregatedAllCandidateSegments$roundness_sd)
#           bestCandidateSegment <- data_aggregatedAllCandidateSegments$segmentNumber[bestIndex]
#           cat("*** Selecting segment with smallest variation in roundness:",bestCandidateSegment,fill=T)
#           #           print(data_aggregatedAllCandidateSegments)
#           
#           data_aggregatedAllCandidateSegments <- data_aggregatedAllCandidateSegments[bestIndex,]
#           candidateSegments <- bestCandidateSegment
#           
#         }
#         
#         if (T) { # length(candidateSegments) == 1
#           #           cat("*** Identified segment:",candidateSegments,fill=T)
#           if (is.null(imageData)) {
#             imageData <- data_aggregatedAllCandidateSegments
#           } else {
#             imageData <- rbind(imageData, data_aggregatedAllCandidateSegments)
#           }
#         } else {
#           print(ggplot(data=allSegments, aes(x=imgNumber, y=s.area, colour=factor(segmentNumber)))+geom_line()+
#                   ggtitle(paste("All segments in", imgFolder)))
#           
#           cat("*** No segment identified in", imgFolder, fill=T)
#         }
#         write.csv(imageData, 'imageData.csv', row.names=F, quote=F)
#       } # for serie (sax_nnn etc) 
#       # TODO consider merge with existing data?
#     } # for Id
#   } else {
#     cat("No dataset in", datasetDir, fill=T)
#   }
# } # for dataset (train/validate/test)

