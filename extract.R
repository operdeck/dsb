# http://www.r-bloggers.com/r-image-analysis-using-ebimage/
# http://www.bioconductor.org/packages/release/bioc/vignettes/EBImage/inst/doc/EBImage-introduction.html

# source("http://bioconductor.org/biocLite.R")
# biocLite("EBImage")

require('EBImage')
require('oro.dicom')
library(ggplot2)
library(dplyr)

pos2coord<-function(pos=NULL, coord=NULL, dim.mat=NULL){
  if(is.null(pos) & is.null(coord) | is.null(dim.mat)){
    stop("must supply either 'pos' or 'coord', and 'dim.mat'")
  }
  if(is.null(pos) & !is.null(coord) & !is.null(dim.mat)){
    pos <- ((coord[,2]-1)*dim.mat[1])+coord[,1] 
    return(pos)
  }
  if(!is.null(pos) & is.null(coord) & !is.null(dim.mat)){
    coord <- matrix(NA, nrow=length(pos), ncol=2)
    coord[,1] <- ((pos-1) %% dim.mat[1]) +1
    coord[,2] <- ((pos-1) %/% dim.mat[1]) +1
    return(coord)
  }
}

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
        
        filez <- list.files(imgFolder, pattern=".*dcm$")
        
        allSegments <- NULL
        for (imgNumber in 1:length(filez)) {
          f <- filez[imgNumber]
          fname <- paste(imgFolder, f, sep="/")
          
          dicom <- readDICOMFile(fname)
          img <- Image(normalize(dicom$img))
          img <- rotate(img,-90)
          # display(img, method = "raster")
          
          img_thresholded <- img > otsu(img) # Otsu’s threshold 
          img_denoised <- medianFilter(img_thresholded, 4)
          img_segmented <- bwlabel(img_denoised)
          img_comb = EBImage::combine(
            normalize(img_segmented),
            img,
            img_thresholded,
            img_denoised
          )
          # 
          
          display(img_comb, title=imgFolder, all=T,method="raster")
          for (i in 1:max(img_segmented)) {
            coords <- as.data.frame(pos2coord(pos=which(img_segmented==i),dim.mat=dim(img_segmented)))
            names(coords) <- c('x','y')
            text(x = mean(coords$x), y = mean(coords$y), label = i, col = "orange")
          }
          
          #     print(table(segmented))
          # now find similar objects in consecutive images so we know what to look for...
          
          # get mean x, mean y and size for all segments...
          # then try match those in consecutive images, assign same color
          segments <- data.frame(segmentNumber=1:max(img_segmented)
                                 ,imgNumber=imgNumber
                                 #                          ,img=fname
          )
          
          for (i in 1:max(img_segmented)) {
            # TODO maybe discard too small segments (size < 10 or so)
            
            # all coordinates of the pixels in the segment:
            coords <- as.data.frame(pos2coord(pos=which(img_segmented==i),dim.mat=dim(img_segmented)))
            names(coords) <- c('x','y')
            
            # some aggregates that characterize the segment:
            segments$x[i] <- mean(coords$x)
            segments$y[i] <- mean(coords$y)
            segments$sd_x[i] <- sd(coords$x)
            segments$sd_y[i] <- sd(coords$y)
            segments$size[i] <- nrow(coords)
            segments$RSD[i] <- sqrt(segments$sd_x[i]^2 + segments$sd_y[i]^2) # indicator for roundness??
            segments$volume[i] <- segments$size[i]*sqrt(segments$size[i]) # ideally
            # TODO find something that characterizes 'roundness'/'circleness'
          }
          
          segments <- dplyr::filter(segments, size > 10) # remove small patches
          segments <- dplyr::filter(segments, RSD < 20) # remove non roundish things (TODO: rel to size)
          
          #   print(segments)
          
          if (is.null(allSegments)) {
            allSegments <- segments
            maxSegmentNumber <- max(allSegments$segmentNumber)
          } else {
            for (i in 1:nrow(segments)) {
              # find 'a' in allSegments such that it matches 'i' close enough
              # TODO make a function
              
              distance <- sqrt((segments$x[i] - allSegments$x)^2 + (segments$y[i] - allSegments$y)^2)
              
              candidates <- which(distance < 10) # etc.
              
              # apply criteria
              # then select the best one
              
              
              a <- which.min(distance)
              if (distance[a] < 10 ){
                #           & # absolute distance
                #             abs(1-allSegments$volume[a]/segments$volume[i]) < 0.20) { # relative volume tolerance
                #           cat(segments$segmentNumber[i],"Identical to",allSegments$segmentNumber[a],fill=T)
                
                segments$segmentNumber[i] <- allSegments$segmentNumber[a]
                #           cat(segments$segmentNumber[i],"Identical to",allSegments$segmentNumber[a],fill=T)
                
              } else {
                
                maxSegmentNumber <- maxSegmentNumber+1
                segments$segmentNumber[i] <- maxSegmentNumber
                
                #           cat(i,"NOT Identical to",allSegments$segmentNumber[a],fill=T)
                #           print(allSegments[a,])
                #           print(segments[i,])
              }
            }
            allSegments <- rbind(allSegments, segments)
            
            for (i in 1:nrow(segments)) {
              text(x = segments$x[i], y = segments$y[i], label = segments$segmentNumber[i], col = "orange")
            }
          } # segments of one image
        } # images
        
        # Now try figure out which image segment is the one of interest
        # this is based on heuristics - or we should do semi-supervised learning
        # todo: apply a 'roundness' criterion as well
        allSegmentsGrouped <- group_by(allSegments, segmentNumber) %>% 
          dplyr::summarize(size_sd = sd(size),
                    size_mean = mean(size),
                    size_p25 = quantile(size)[2],
                    size_p50 = quantile(size)[3],
                    size_p75 = quantile(size)[4],
                    vol_min = min(volume),
                    vol_max = max(volume),
                    vol_p25 = quantile(volume)[2],
                    vol_p50 = quantile(volume)[3],
                    vol_p75 = quantile(volume)[4],
                    img_used = n(),
                    img_cnt = length(filez))
        allSegmentsGrouped$dataset <- dataset
        allSegmentsGrouped$Id <- as.integer(Id)
        allSegmentsGrouped$series <- serieFolder
        allSegmentsGrouped$series_idx <- which(serieFolder == serieFolders)
        allSegmentsGrouped$series_cnt <- length(serieFolders)

        candidateSegments <- filter(allSegmentsGrouped, 
                                    size_mean > 200, 
                                    size_sd > 50)
        allSegments$segmentNumber <- factor(allSegments$segmentNumber,levels=unique(allSegments$segmentNumber))
        print(ggplot(data=filter(allSegments, segmentNumber %in% candidateSegments$segmentNumber), 
                     aes(x=imgNumber, y=size, colour=segmentNumber))+geom_line()+ggtitle(imgFolder))
        
        finalSegment <- intersect(unique(filter(allSegments, size>400) %>% select(segmentNumber)),
                                  unique(filter(allSegments, size<250) %>% select(segmentNumber)))
        
        if (nrow(finalSegment) >= 1) {
          print("Identified one or more segments:")
          seriesData <- filter(allSegmentsGrouped, segmentNumber %in% finalSegment$segmentNumber)
          seriesData$in_series_cnt <- nrow(seriesData) # number of candidate segments in this series
          seriesData$in_series_rank <- 1:nrow(seriesData) # assume an ordering (TODO we dont have that yet)
          print(seriesData)
          if (is.null(imageData)) {
            imageData <- seriesData
          } else {
            imageData <- rbind(imageData, seriesData)
          }
          # keep 'allSegments' info for this one segment (across the images)
        } else {
          print("No segment identified.")
        }
      } # for serie (sax_nnn etc) 
    } # for Id
  } else {
    cat("No dataset in", datasetDir, fill=T)
  }
} # for dataset (train/validate/test)

write.csv(imageData, 'imageData.csv', row.names=F, quote=F)
