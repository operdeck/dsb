# http://www.r-bloggers.com/r-image-analysis-using-ebimage/
# http://www.bioconductor.org/packages/release/bioc/vignettes/EBImage/inst/doc/EBImage-introduction.html

# source("http://bioconductor.org/biocLite.R")
# biocLite("EBImage")

library(Hmisc)
require('EBImage')
require('oro.dicom')
library(tidyr)
library(data.table)
library(ggplot2)
library(dplyr)
library(xgboost)
library(uuid)

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

readInteger <- function(p = "Enter an integer: ")
{ 
  n <- readline(prompt=p)
  if(!grepl("^[-+]?[0-9]+$",n))
  {
    return(readInteger(p))
  }
  
  return(as.integer(n))
}

# create mid-ordered sequence
midOrderSeq <- function(n) {
  abs(n %/% 2 + 1 - seq(n))*2 - ((1+sign(n %/% 2 + 1 - seq(n)))%/%2) + 1
}

# Extension of above so that there are always pairs of slices present
midOrderSeqKeepSibblingsTogether <- function(n) {
  s <- abs((seq(n)-(n%/%2))%/%2)+1
  if (length(s) > 1) {
    s[1] <- s[2]
    s[length(s)] <- s[length(s)-1]
  }
  return(s)
}

getSegmentFile <- function(ds) {
  paste("segments-",ds,".csv",sep="")
}

getImageFile <- function(ds) {
  paste("data",unique(ds$FileName),sep="/")
}

# Plot interesting attributes over time for one slice
plotSlice <- function(slice) {
  plotData <- mutate(slice, 
                     area.radius.mean = pi*radius.mean^2,
                     area.radius.max = pi*radius.max^2,
                     area.radius.min = pi*radius.min^2) %>% 
    gather(metric, area, starts_with("area"))
  if (nrow(slice) > 1 & any(!is.na(plotData$area))) {
    print(ggplot(plotData, aes(x=Time, y=area, colour=metric))+geom_line()+geom_point()+
            ggtitle(paste("Segment area over Time for ID",
                          unique(slice$Id),"Slice",unique(slice$Slice))))
  } else {
    cat("No image with identified LV at all for slice:", unique(slice$Id), unique(slice$Slice), fill=T)
  }
}

showSegmentLabels <- function(segments, textColour="red") {
  if (nrow(segments) > 0) {
    for (i in 1:nrow(segments)) {
      text(x = segments$m.cx[i], y = segments$m.cy[i], 
           label = segments$segIndex[i], col = textColour)
    }
  }
}



getImageList <- function(type = "sax", playlist=NULL)
{
  if (is.null(playlist)) {
    playlist <- fread("imagelist.csv")
  }
  r <- filter(playlist, ImgType==type)
  setkey(r, Id, Dataset, ImgType, Slice, FileName)
  return(r)
}

getSliceList <- function(type = "sax", playlist=NULL)
{
  oneSliceGroup <- select(getImageList(type, playlist), 
                          Id, Dataset, ImgType, Offset, 
                          starts_with("Slice"),
                          starts_with("Patient")) %>% 
    arrange(Dataset, Id, ImgType, Slice) %>% unique()
  setkey(oneSliceGroup, Dataset, Id, ImgType, Slice)
  return(oneSliceGroup)
}

getIdList <- function(type = "sax", playlist=NULL)
{
  ids <- select(getSliceList(type, playlist), 
                Id, Dataset, ImgType, 
                SliceCount,
                starts_with("Patient")) %>% unique()
  setkey(ids, Id, Dataset, ImgType)
  return(ids)
}

# # Given a list of images with segments for one slice, return a list with every image
# # and only 0 or 1 segments per image, such that the overal distance to the Left Ventricle
# # is minimized.
# TODO - not tested, not finished
# getImagesWithLVSegments <- function(slice)
# {
#   # for ALL segments within the slice
#   for (lvCandidateIndex in seq(nrow(slice))) {
#     # for ALL other segments (at once)
#     # calculate distance to lvCandidate
#     slice$distToLVCandidate <- sqrt((slice$m.cx - slice$m.cx[lvCandidateIndex])^2 + 
#                                       (slice$m.cy - slice$m.cy[lvCandidateIndex])^2)
#     # identify the segment with smallest distance per Image
#     identifiedLVSegments <- group_by(slice, Time) %>% summarise(segLV = segIndex[which.min(distToLVCandidate)])
#     slice <- left_join(slice, identifiedLVSegments, by=c("Time"))
#     slice$isLV <- (slice$segIndex == slice$segLV)
#     
#     # calculate sum of pLV or something for the slices...
#     # keep track of the lowest pLV sum
#   }
# }

# TODO: add variation over time
# add rank of area (e.g.) vs other segments, maybe also over time

# Create segment prediction dataset
# ds could be just one Id/Slice, even be just the segments for a single image
createSegmentPredictSet <- function(ds)
{
  ds[, areaRank := frankv(s.area, order=-1L, ties.method="dense"), by=c("Id","Slice","Time")] # fast rank (data.table)
  ds[, areaRelSize := s.area/sum(s.area,na.rm=T), by=c("Id","Slice","Time")] # relative size
  
  dsByTime <- unique(select(ds,Id,Slice,Time))
  cat("DS has",nrow(dsByTime),"unique times of avg size",nrow(ds)/nrow(dsByTime),fill=T)
  dsBySlice <- unique(select(ds,Id,Slice))
  cat("DS has",nrow(dsBySlice),"unique slices of avg size",nrow(ds)/nrow(dsBySlice),fill=T)
  #for (i in nrow(dsBySlice)) {
  #  segments <- which(ds[])
  
  # take top-5 largest neighbours and add, for example,
  # top5.weighted.area (= sum distance*area; scaled) etc
  
  ds <- mutate(ds,
               areaMultiplier = PixelSpacing.x * PixelSpacing.y * 0.01, # 1 mL = 1000 mm3
               lengthMultiplier = sqrt(areaMultiplier),
               
               area = s.area*areaMultiplier,
               area.ellipse = pi*s.radius.min*s.radius.max*areaMultiplier,
               
               perimeter = s.perimeter*lengthMultiplier,
               radius.mean = s.radius.mean*lengthMultiplier,
               radius.min = s.radius.min*lengthMultiplier,
               radius.max = s.radius.max*lengthMultiplier,
               radius.var = sqrt(s.radius.sd)*lengthMultiplier,
               
               roi.size = pi*ROI.r*ROI.r*areaMultiplier,
               area.rel = area/roi.size,
               
               majoraxis = m.majoraxis*lengthMultiplier,
               roundness = 4*pi*area/(perimeter^2),
               
               slicePct = SliceIndex/SliceCount) %>%
    dplyr::rename(eccentricity = m.eccentricity,
           theta = m.theta) %>%
    select(-FileName,               # is an ID
           -segIndex,               # is (part of) ID
           -PixelSpacing.x,         # used to build another predictor
           -PixelSpacing.y,         # used to build another predictor
           -areaMultiplier,         # is temporary variable
           -lengthMultiplier,       # is temporary variable
           -m.majoraxis,            # is orientation dependent
           -starts_with("s."),      # renamed and scaled
           -starts_with("ROI."),    # radius used, others are position dependent
           -m.cx, -m.cy,            # is position dependent
           -UUID,                   # is an ID
           -Id, -Slice, -Time,      # is (part of) ID 
           -Dataset, -ImgType,      # is (part of) ID
           -SliceIndex, -SliceCount,# used to build another predictor
           #-distToROI,              # position dependent - not available always
           -Offset)                 # slice location is the better predictor
  
  # Keep only numerics plus the outcome (logic)
  ds <- ds[,sapply(ds, function(c) { return (is.numeric(c) | is.logical(c)) }),with=F] # keep only numerics

  # ID slice img   seg x/y <segment details>
  #
  #  1     1   1   0.1 0.2   0.5 0.6 
  #  1     1   1   0.2 0.3   0.3 0.7
  #  1     1   1   0.3 0.4   0.4 0.2 
  #
  #  1     1   2   0.1 0.2   0.5 0.6 
  #  1     1   2   0.2 0.3   0.3 0.7
  #  1     1   2   0.3 0.4   0.4 0.2 
  
  # TODO: add difficult aggregates over other images for the same slice
  
  return (ds[,sort(names(ds)),with=F])
}

# Create image prediction dataset
createImagePredictSet <- function(ds)
{
  ds <- mutate(ds,
               areaMultiplier = PixelSpacing.x * PixelSpacing.y * 0.01,  # 1 mL = 1000 mm3
               lengthMultiplier = sqrt(areaMultiplier),
               
               area = s.area*areaMultiplier,
               area.ellipse = pi*s.radius.min*s.radius.max*areaMultiplier,
               
               perimeter = s.perimeter*lengthMultiplier,
               radius.mean = s.radius.mean*lengthMultiplier,
               radius.min = s.radius.min*lengthMultiplier,
               radius.max = s.radius.max*lengthMultiplier,
               radius.var = sqrt(s.radius.sd)*lengthMultiplier,
               
               majoraxis = m.majoraxis*lengthMultiplier,
               roundness = 4*pi*area/(perimeter^2),
               
               slicePct = SliceIndex/SliceCount) %>%
    select(-starts_with("PixelSpacing."),
           -segIndex,               # segment attribute
           -starts_with("m."),      # segment attribute
           -m.majoraxis,            # is orientation dependent
           -areaMultiplier,         # is temporary variable
           -lengthMultiplier,       # is temporary variable
           -starts_with("s."),      # renamed and scaled
           -starts_with("ROI."),    # is position dependent
           -m.cx, -m.cy)            # is position dependent
  
  return (ds[,sort(names(ds)),with=F])
}

getSegmentedImageDir <- function(ds)
{
  paste("segmented",ds$Dataset, ds$Id, ds$Slice,sep="/")
}

getSegmentedImageFile <- function(ds, dirName = getSegmentedImageDir(ds))
{
  paste(dirName,paste(paste("SEG", ds$Offset, ds$Time, sep="-"), ".png", sep=""), sep="/")
}

