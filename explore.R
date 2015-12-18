# http://www.r-bloggers.com/r-image-analysis-using-ebimage/
# http://www.bioconductor.org/packages/release/bioc/vignettes/EBImage/inst/doc/EBImage-introduction.html

# source("http://bioconductor.org/biocLite.R")
# biocLite("EBImage")

require('EBImage')
require('oro.dicom')

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

for (dataset in c('train','validate','test')) {
  datasetDir <- paste("data",dataset,sep="/")
  if (dir.exists(datasetDir)) {
    caseFolders <- list.dirs(datasetDir, recursive=F, full.names=F)
    for (caseFolder in caseFolders) {
    #   print(caseFolder) # = nr 1 .. 500
      caseFolder <- paste(datasetDir, caseFolder, "study", sep="/")
      
      serieFolders <- list.dirs(caseFolder, recursive=F, full.names=F)
      # NB excluding the 2 chamber and 4 chamber views here, only the short axis track
      serieFolders <- serieFolders[ grepl("sax_[[:digit:]]+$", serieFolders) ]
    #   print(imgFolders)
      
      for (serieFolder in serieFolders) {
        sliceFolder <- paste(caseFolder, serieFolder, sep="/")
    
        print(sliceFolder) # inside here are the .dcm frames
      }
    }
  } else {
    cat("No dataset in", datasetDir, fill=T)
  }
}

# Image processing creates 'image_features.csv'
#
# case   - view     - slice    <image measurements>
# 1..500   ~10 - 20   ~1..30   [eg median]
#

# from this, predict Vd and Vs

# aggregate up - by case+view min/max of measurements
# by case: mean/stddev of aggregations

# model with 600 classes??
# caret can do multi-class see http://stackoverflow.com/questions/15585501/usage-of-caret-with-gbm-method-for-multiclass-classification
# probs are cumulative?
# see https://www.kaggle.com/c/second-annual-data-science-bowl/details/evaluation

folder <- 'data/train/10/study/sax_11'
filez <- list.files(folder, pattern=".*dcm$")

allSegments <- NULL
for (f in filez) {
  fname <- paste(folder, f, sep="/")
 
  dicom <- readDICOMFile(fname)
  img <- Image(normalize(dicom$img))
  img <- rotate(img,-90)
  # display(img, method = "raster")
  
  img_comb = EBImage::combine(
    img,
    img > otsu(img) # Otsu’s threshold 
  )
  # display(img_comb, all=TRUE, method="raster")
  
  img_median = medianFilter(img_comb, 4)
  # display(img_median, all=T,method="raster")
  
  # segmentation
  segmented <- bwlabel(img_median[,,2])
  display(colorLabels( segmented), method="raster")
  
  print(table(segmented))
  # now find similar objects in consecutive images so we know what to look for...
  
  # get mean x, mean y and size for all segments...
  # then try match those in consecutive images, assign same color
  segments <- data.frame(segmentNumber=1:max(segmented),
                         img=fname)
  for (i in 1:max(segmented)) {
    # TODO maybe discard too small segments (size < 10 or so)
    
    # all coordinates of the pixels in the segment:
    coords <- as.data.frame(pos2coord(pos=which(segmented==i),dim.mat=dim(segmented)))
    names(coords) <- c('x','y')
    
    # some aggregates that characterize the segment:
    segments$x[i] <- mean(coords$x)
    segments$y[i] <- mean(coords$y)
    segments$sd_x[i] <- sd(coords$x)
    segments$sd_y[i] <- sd(coords$y)
    segments$size[i] <- nrow(coords)
    segments$volume[i] <- segments$size[i]*sqrt(segments$size[i]) # ideally
    # TODO find something that characterizes 'roundness'/'circleness'
  }
  segments <- filter(segments, size > 10) # remove small patches
  
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
      if (distance[a] < 10 & # absolute distance
            abs(1-allSegments$volume[a]/segments$volume[i]) < 0.20) { # relative volume tolerance
        segments$segmentNumber[i] <- allSegments$segmentNumber[a]
        
        cat(i,"Identical to",a,fill=T)
      } else {
        maxSegmentNumber <- maxSegmentNumber+1
        segments$segmentNumber[i] <- maxSegmentNumber
        
        cat(i,"NOT Identical to",a,fill=T)
        print(allSegments[a,])
        print(segments[i,])
      }
    }
    allSegments <- rbind(allSegments, segments)
    
    for (i in 1:nrow(segments)) {
      text(x = segments$x[i], y = segments$y[i], label = segments$segmentNumber[i], col = "orange")
    }
  }
}


