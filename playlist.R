# List all images and save in a .csv file. List contains DICOM meta info
# and should be run prior to any of the other processing.

source("Util.R")

# contains all images
# with image DICOM meta info
# maybe also add image dimensions
# using SliceIndex instead of SliceRelIndex
# removed slicelist.csv from git, add imagelist.csv
# FileName added
# for all image types

listSliceImages <- function(sliceInfo) {
  imgFolder <- getImageFolder(sliceInfo)
  imgFiles <- list.files(imgFolder, pattern=".*\\.dcm$")
  imageMetaData <- data.table(sliceInfo,
  
                              FileName = imgFiles,
                              
                              # this will not work for some images that have
                              # a different format
                              Offset = as.integer(gsub("^IM-(\\d{4,})-(\\d{4})\\.dcm$", "\\1", imgFiles)),
                              Time = as.integer(gsub("^IM-(\\d{4,})-(\\d{4})\\.dcm$", "\\2", imgFiles)),
                              
                              # meta info
                              PixelSpacing.x = rep(NA, length(imgFiles)),
                              PixelSpacing.y = rep(NA, length(imgFiles)),
                              SliceLocation = rep(NA, length(imgFiles)),
                              SliceThickness = rep(NA, length(imgFiles)))

  for (i in seq(nrow(imageMetaData))) {
    dicomHeader <- readDICOMFile(paste(getImageFolder(imageMetaData[i,]),imageMetaData$FileName[i],sep="/"), 
                                 pixelData = TRUE)
    pixelSpacing <- extractHeader(dicomHeader$hdr, "PixelSpacing", numeric=FALSE)
    
    imageMetaData$PixelSpacing.x[i] <- as.numeric(gsub("^(.*) (.*)$", "\\1", pixelSpacing))
    imageMetaData$PixelSpacing.y[i] <- as.numeric(gsub("^(.*) (.*)$", "\\2", pixelSpacing))
    imageMetaData$SliceLocation[i]  <- extractHeader(dicomHeader$hdr, "SliceLocation", numeric=TRUE)
    imageMetaData$SliceThickness[i] <- extractHeader(dicomHeader$hdr, "SliceThickness", numeric=TRUE)
    img <- Image(dicomHeader$img)
    if (dim(img)[1] > dim(img)[2]) {
      imageMetaData$Dim.x <- dim(img)[1]
      imageMetaData$Dim.y <- dim(img)[2]
    } else {
      imageMetaData$Dim.x <- dim(img)[2]
      imageMetaData$Dim.y <- dim(img)[1]
    }
  }
  
  return(imageMetaData)
}

playlist <- NULL
for (dataset in getAllDatasetFolders()) {
  datasetDir <- paste("data",dataset,sep="/")
  Ids <- sort(as.integer(list.dirs(datasetDir, recursive=F, full.names=F)))
  for (Id in Ids) {
    IdFolder <- paste(datasetDir, Id, "study", sep="/")
    subFolders <- list.dirs(IdFolder, recursive=F, full.names=F)
    sliceMetaData <- data.table(Id = Id, 
                                Dataset = dataset,
                                ImgType = gsub("(.+)_([[:digit:]]+)$", "\\1", subFolders),
                                Slice = as.integer(gsub("(.+)_([[:digit:]]+)$", "\\2", subFolders))) %>%
      arrange(Dataset, ImgType, Slice)
    for (t in unique(sliceMetaData$ImgType)) {
      sliceMetaDataByImgType <- filter(sliceMetaData, ImgType==t)
      
      sliceMetaDataByImgType$SliceCount = nrow(sliceMetaDataByImgType)
      sliceMetaDataByImgType$SliceIndex = seq(nrow(sliceMetaDataByImgType))
      sliceMetaDataByImgType$SliceOrder = midOrderSeqKeepSibblingsTogether(nrow(sliceMetaDataByImgType))
      
      print(sliceMetaDataByImgType)
      for (s in seq(nrow(sliceMetaDataByImgType))) {
        imageMetaData <- listSliceImages(sliceMetaDataByImgType[s,])
        
        if (is.null(playlist)) {
          playlist <- imageMetaData
        } else {
          playlist <- rbind(playlist, imageMetaData)
        }
      }
    }
    write.csv(playlist, "imagelist.csv", row.names=F)
  }
}

getImageList <- function(type = "sax")
{
  playlist <- fread("imagelist.csv")
  r <- filter(playlist, ImgType==type)
  setkey(r, Id, Dataset, ImgType, Slice)
  return(r)
}

getSliceList <- function(type = "sax")
{
  slices <- select(getImageList(type), Id, Dataset, ImgType, Slice, Offset, starts_with("slice")) %>% unique()
  setkey(slices, Id, Dataset, ImgType, Slice)
  return(slices)
}

getIdList <- function(type = "sax")
{
  ids <- select(getSliceList(type), Id, Dataset, ImgType, SliceCount) %>% unique()
  setkey(ids, Id, Dataset, ImgType)
  return(ids)
}

# Show quick summary of the datasets & save for downstream use
print("Results:")
print(select(getIdList(), Dataset, Id) %>% group_by(Dataset) %>% summarise(n_Ids = n()))
print(select(getSliceList(), Dataset, Id, Slice) %>% group_by(Dataset) %>% summarise(n_Slices = n()))
print(select(getImageList(), Dataset, Id, Slice, FileName) %>% group_by(Dataset) %>% summarise(n_Images = n()))

