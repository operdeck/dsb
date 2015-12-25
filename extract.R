# http://www.r-bloggers.com/r-image-analysis-using-ebimage/
# http://www.bioconductor.org/packages/release/bioc/vignettes/EBImage/inst/doc/EBImage-introduction.html

# source("http://bioconductor.org/biocLite.R")
# biocLite("EBImage")

require('EBImage')
require('oro.dicom')
library(ggplot2)
library(dplyr)

# pos2coord<-function(pos=NULL, coord=NULL, dim.mat=NULL){
#   if(is.null(pos) & is.null(coord) | is.null(dim.mat)){
#     stop("must supply either 'pos' or 'coord', and 'dim.mat'")
#   }
#   if(is.null(pos) & !is.null(coord) & !is.null(dim.mat)){
#     pos <- ((coord[,2]-1)*dim.mat[1])+coord[,1] 
#     return(pos)
#   }
#   if(!is.null(pos) & is.null(coord) & !is.null(dim.mat)){
#     coord <- matrix(NA, nrow=length(pos), ncol=2)
#     coord[,1] <- ((pos-1) %% dim.mat[1]) +1
#     coord[,2] <- ((pos-1) %/% dim.mat[1]) +1
#     return(coord)
#   }
# }

# DICOM has all kind of meta info which we're not using right now

# after systole = min
# after diastole = max
#  ejection fraction = 100 * (Vd - Vs)/Vd


# Image processing creates 'image_features.csv'
#
# Id   - view     - slice    <image measurements>
# 1..500   ~10 - 20   ~1..30   [eg median]
#

# from this, predict Vd and Vs

# aggregate up - by Id+view min/max of measurements
# by Id: mean/stddev of aggregations

# model with 600 classes??
# caret can do multi-class see http://stackoverflow.com/questions/15585501/usage-of-caret-with-gbm-method-for-multiclass-classification
# probs are cumulative?
# see https://www.kaggle.com/c/second-annual-data-science-bowl/details/evaluation

# see bio image detection stuff
# http://bioconductor.wustl.edu/bioc/vignettes/EBImage/inst/doc/AnalysisWithEBImage.pdf

imageData <- NULL

stopAtImg <- ""
stopAtImg <- "data/train/1/study/sax_10/IM-4562-0027.dcm"

for (dataset in c('train','validate','test')) {
  datasetDir <- paste("data",dataset,sep="/")
  if (dir.exists(datasetDir)) {
    IdFolders <- list.dirs(datasetDir, recursive=F, full.names=F)
    for (Id in IdFolders) {
      #   print(IdFolder) # = nr 1 .. 500
      IdFolder <- paste(datasetDir, Id, "study", sep="/")
      
      serieFolders <- list.dirs(IdFolder, recursive=F, full.names=F)
      # NB excluding the 2 chamber and 4 chamber views here, only the short axis track
      serieFolders <- serieFolders[ grepl("sax_[[:digit:]]+$", serieFolders) ]
      
      for (serieFolder in serieFolders) { 
        # serieFolder <- 'data/train/500/study/sax_22'
        imgFolder <- paste(IdFolder, serieFolder, sep="/")
        
        print("Processing:")
        print(imgFolder)
        
        seriesImageFiles <- list.files(imgFolder, pattern=".*dcm$")
        
        allSegments <- NULL
        for (imgNumber in 1:length(seriesImageFiles)) {
          f <- seriesImageFiles[imgNumber]
          fname <- paste(imgFolder, f, sep="/")
          
          print(fname)
          
          dicom <- readDICOMFile(fname)
          img <- Image(normalize(dicom$img))
          img <- rotate(img,-90)
          # display(img, method = "raster")
          
          img_thresholded <- img > otsu(img) # Otsu’s threshold 
          img_denoised <- medianFilter(img_thresholded, 4)
          img_segmented <- bwlabel(img_denoised)
          img_comb = EBImage::combine(
            colorLabels(img_segmented), # TODO try preserve colours
            toRGB(img),
            toRGB(img_thresholded),
            toRGB(img_denoised)
          )
          
          display(img_comb, title=imgFolder, all=T,method="raster")
          #           for (i in 1:max(img_segmented)) {
          #             coords <- as.data.frame(pos2coord(pos=which(img_segmented==i),dim.mat=dim(img_segmented)))
          #             names(coords) <- c('x','y')
          #             text(x = mean(coords$x), y = mean(coords$y), label = i, col = "orange")
          #           }
          
          # Calculate image features. EBImage does a good job at this, see doc
          # https://www.bioconductor.org/packages/3.3/bioc/manuals/EBImage/man/EBImage.pdf
          # shape factors https://en.wikipedia.org/wiki/Shape_factor_(image_analysis_and_microscopy)
          # http://www.cs.uu.nl/docs/vakken/ibv/reader/chapter8.pdf
          segments <- cbind(as.data.frame(computeFeatures.shape(img_segmented)), 
                            as.data.frame(computeFeatures.moment(img_segmented)))
          segments <- dplyr::filter(segments, s.area > 10) # early filter to remove small patches
          segments$roundness <- 4*pi*segments$s.area/(segments$s.perimeter^2)
#           segments$compactness <- 1/segments$roundness
          segments$segmentNumber <- 1:nrow(segments)
          segments$imgNumber <- imgNumber
          #           for (i in 1:max(img_segmented)) {
          #             # all coordinates of the pixels in the segment:
          #             coords <- as.data.frame(pos2coord(pos=which(img_segmented==i),dim.mat=dim(img_segmented)))
          #             names(coords) <- c('x','y')
          #           }
          
          #           segments <- dplyr::filter(segments, RSD < 20) # remove non roundish things (TODO: rel to size)
          
          if (is.null(allSegments)) {
            allSegments <- segments
            maxSegmentNumber <- max(allSegments$segmentNumber)
          } else {
            for (i in 1:nrow(segments)) {
              # find previous segments at approx the same position
              distance <- sqrt((segments$m.cx[i] - allSegments$m.cx)^2 + (segments$m.cy[i] - allSegments$m.cy)^2)
              #               candidates <- which(distance < 10) # etc.
              a <- which.min(distance)
              if (distance[a] < 10 ){
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
        } # images
        #         allSegments$segmentNumber <- factor(allSegments$segmentNumber,levels=unique(allSegments$segmentNumber))
        allSegments$segmentNumber <- as.integer(allSegments$segmentNumber)
        
        # Now try figure out which image segment is the one of interest
        # must be large and small, maybe also filter on size or roundness
#         largeSegments <- unique(filter(allSegments, s.area>400)$segmentNumber)
        smallSegments <- filter(allSegments, s.area<300)
#         candidateSegments <- intersect(largeSegments, smallSegments)
        candidateSegments <- (group_by(smallSegments, segmentNumber) %>% 
                                summarise(size_variance = sd(s.area)) %>% arrange(desc(size_variance)) %>% ungroup() %>% 
                                mutate(size_variance_scaled=scale(size_variance)) %>% filter(size_variance_scaled > 1))$segmentNumber
        
        data_scaled <- as.data.frame(cbind(scale(allSegments[,sapply(allSegments,class) == "numeric"]),
                                           allSegments[,sapply(allSegments,class) != "numeric"]))
        data_aggregatedAllCandidateSegments <- NULL
        for (s in candidateSegments) {
          data_candidateSegment <- filter(data_scaled, segmentNumber == s) %>% select(-segmentNumber, -imgNumber)
          data_mean <- sapply(data_candidateSegment, mean)
          names(data_mean) <- paste(names(data_candidateSegment), "mean", sep="_")
          data_sd <- sapply(data_candidateSegment, sd)
          names(data_sd) <- paste(names(data_candidateSegment), "sd", sep="_")
          data_aggregatedCandidateSegment <- c(segmentNumber=s, data_mean, data_sd, 
                                               img_present = nrow(data_candidateSegment))
          #           s.area_p25 = quantile(s.area)[2],
          #           s.area_p50 = quantile(s.area)[3],
          #           s.area_p75 = quantile(s.area)[4],
          if (is.null(data_aggregatedAllCandidateSegments)) {
            data_aggregatedAllCandidateSegments <- as.data.frame(t(data_aggregatedCandidateSegment))
          } else {
            # only happens when there are multiple candidates
            data_aggregatedAllCandidateSegments <- rbind(data_aggregatedAllCandidateSegments,t(data_aggregatedCandidateSegment))
          }
        }
        
        if (length(candidateSegments) >= 1) {
          row.names(data_aggregatedAllCandidateSegments) <- 1:nrow(data_aggregatedAllCandidateSegments)
          data_aggregatedAllCandidateSegments$dataset <- dataset
          data_aggregatedAllCandidateSegments$Id <- as.integer(Id)
          data_aggregatedAllCandidateSegments$series <- serieFolder
          data_aggregatedAllCandidateSegments$series_idx <- which(serieFolder == serieFolders)
          data_aggregatedAllCandidateSegments$series_cnt <- length(serieFolders)
          data_aggregatedAllCandidateSegments$img_cnt <- length(seriesImageFiles)
        }
        
        if (length(candidateSegments) > 1) {
          cat("*** Identified multiple segments:",candidateSegments,fill=T)
          
          plotData <- filter(data_scaled, segmentNumber %in% candidateSegments) %>% 
            select(m.eccentricity, roundness, s.area, segmentNumber, imgNumber) %>%
            gather(Feature, Value, -segmentNumber, -imgNumber)
          print(ggplot(data=plotData, aes(x=imgNumber, y=Value, linetype=Feature, color=factor(segmentNumber)))+geom_line()+
                  ggtitle(paste("Candidate segments",imgFolder)))
          
          bestIndex <- which.min(data_aggregatedAllCandidateSegments$roundness_sd)
          bestCandidateSegment <- data_aggregatedAllCandidateSegments$segmentNumber[bestIndex]
          cat("*** Selecting segment with smallest variation in roundness:",bestCandidateSegment,fill=T)
#           print(data_aggregatedAllCandidateSegments)
          
          data_aggregatedAllCandidateSegments <- data_aggregatedAllCandidateSegments[bestIndex,]
          candidateSegments <- bestCandidateSegment

        }
        
        if (length(candidateSegments) == 1) {
          cat("*** Identified segment:",candidateSegments,fill=T)
          if (is.null(imageData)) {
            imageData <- data_aggregatedAllCandidateSegments
          } else {
            imageData <- rbind(imageData, data_aggregatedAllCandidateSegments)
          }
        } else {
          print(ggplot(data=allSegments, aes(x=imgNumber, y=s.area, colour=factor(segmentNumber)))+geom_line()+
                  ggtitle(paste("All segments in", imgFolder)))
          
          cat("*** No segment identified in", imgFolder, fill=T)
        }
      } # for serie (sax_nnn etc) 
    } # for Id
  } else {
    cat("No dataset in", datasetDir, fill=T)
  }
} # for dataset (train/validate/test)

write.csv(imageData, 'imageData.csv', row.names=F, quote=F)
