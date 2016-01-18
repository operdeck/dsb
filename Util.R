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

# getAllDatasetFolders <- function() {
#   f <- c('train','validate','test')
#   f[sapply(paste("data",f,sep="/"),dir.exists)]
# }

## TODO: clean up both of these - should be obsolete
# datasetFolders <- c('train','validate','test')
# datasetFoldersForSegmentDetection <- c('train','validate')

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

# getImageFolder <- function(entry) {
#   paste("data",
#         entry$Dataset,
#         entry$Id,
#         "study",
#         paste(entry$ImgType, entry$Slice, sep="_"),
#         sep="/")  
# }
# 
# getImageFile <- function(entry) {
#   paste(getImageFolder(entry), 
#         paste("IM-",sprintf("%04d",entry$Offset),"-",sprintf("%04d",entry$Time),".dcm",sep=""),
#         sep="/")
# }

# Plot interesting attributes over time for one slice
plotSlice <- function(slice) {
  plotData <- mutate(slice, 
                     area.radius.mean = pi*radius.mean^2,
                     area.radius.max = pi*radius.max^2,
                     area.radius.min = pi*radius.min^2) %>% 
    gather(metric, area, starts_with("area"))
  if (nrow(slice) > 1) {
    print(ggplot(plotData, aes(x=Time, y=area, colour=metric))+geom_line()+geom_point()+
            ggtitle(paste("Segment area over Time for ID",
                          unique(slice$Id),"Slice",unique(slice$Slice))))
  } else {
    cat("No image with identified LV at all for slice:", unique(slice$Id), unique(slice$Slice), fill=T)
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
  oneSliceGroup <- select(getImageList(type, playlist), Id, Dataset, ImgType, Offset, starts_with("Slice")) %>% 
    arrange(Dataset, Id, ImgType, Slice) %>% unique()
  setkey(oneSliceGroup, Dataset, Id, ImgType, Slice)
  return(oneSliceGroup)
}

getIdList <- function(type = "sax", playlist=NULL)
{
  ids <- select(getSliceList(type, playlist), Id, Dataset, ImgType, SliceCount) %>% unique()
  setkey(ids, Id, Dataset, ImgType)
  return(ids)
}

# # Given a list of images with segments for one slice, return a list with every image
# # and only 0 or 1 segments per image, such that the overal distance to the Left Ventricle
# # is minimized.
# TODO - not tested, not finished
getImagesWithLVSegments <- function(slice)
{
  # for ALL segments within the slice
  for (lvCandidateIndex in seq(nrow(slice))) {
    # for ALL other segments (at once)
    # calculate distance to lvCandidate
    slice$distToLVCandidate <- sqrt((slice$m.cx - slice$m.cx[lvCandidateIndex])^2 + 
                                      (slice$m.cy - slice$m.cy[lvCandidateIndex])^2)
    # identify the segment with smallest distance per Image
    identifiedLVSegments <- group_by(slice, Time) %>% summarise(segLV = segIndex[which.min(distToLVCandidate)])
    slice <- left_join(slice, identifiedLVSegments, by=c("Time"))
    slice$isLV <- (slice$segIndex == slice$segLV)
    
    # calculate sum of pLV or something for the slices...
    # keep track of the lowest pLV sum
  }
}
