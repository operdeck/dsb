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
# Slice itself not so relevant?

print("Listing image files")
filez <- list.files("data", pattern=".*\\.dcm", recursive=T)
filePattern <- "^(.*)/([[:digit:]]+)/study/(.+)_([[:digit:]]+)/IM-(\\d{1,5})-(\\d{1,5})-{0,1}(\\d{0,5})\\.dcm$"
playlist <- data.table(FileName = filez,
                       Dataset = gsub(filePattern, "\\1", filez),
                       Id = as.numeric(gsub(filePattern, "\\2", filez)),
                       ImgType = gsub(filePattern, "\\3", filez),
                       Slice = as.numeric(gsub(filePattern, "\\4", filez)),
                       Offset = as.numeric(gsub(filePattern, "\\5", filez)),
                       Time = as.numeric(gsub(filePattern, "\\6", filez)),
                       SubSlice = as.numeric(gsub(filePattern, "\\7", filez)))
playlist$Slice <- ifelse(is.na(playlist$SubSlice), playlist$Slice, playlist$Slice*1000 + playlist$SubSlice)
playlist <- select(playlist, -SubSlice) %>% arrange(Dataset, Id, ImgType, Slice)
setkey(playlist, Dataset, Id, ImgType, Slice)

print("Creating slice groups")
sliceMetaData <- select(playlist, Dataset, Id, ImgType, Slice) %>% unique() %>% arrange(Dataset, Id, ImgType, Slice)
sliceMetaData$One <- 1 # trick to do a count within the := operator of data.table
sliceMetaData[, 
            c("SliceCount", "SliceIndex", "SliceOrder") := list(sum(One), seq(sum(One)), midOrderSeqKeepSibblingsTogether(sum(One))),
            by=c("Dataset", "Id", "ImgType")]
playlist <- left_join(playlist, sliceMetaData, by=c("Dataset", "Id", "ImgType", "Slice"))

print("Reading image data")
for (n in seq(nrow(playlist))) {
  img <- playlist[n,]
  print(img)
  dicomHeader <- readDICOMFile(paste("data", img$FileName, sep="/"), pixelData = FALSE)
  pixelSpacing <- extractHeader(dicomHeader$hdr, "PixelSpacing", numeric=FALSE)
  
  playlist$PixelSpacing.x[n] <- as.numeric(gsub("^(.*) (.*)$", "\\1", pixelSpacing))
  playlist$PixelSpacing.y[n] <- as.numeric(gsub("^(.*) (.*)$", "\\2", pixelSpacing))
  playlist$SliceLocation[n]  <- extractHeader(dicomHeader$hdr, "SliceLocation", numeric=TRUE)
  playlist$SliceThickness[n] <- extractHeader(dicomHeader$hdr, "SliceThickness", numeric=TRUE)
}
# TODO - fill in missings
write.csv(playlist, "imagelist.csv", row.names=F)

getImageList <- function(type = "sax")
{
  playlist <- fread("imagelist.csv")
  r <- filter(playlist, ImgType==type)
  setkey(r, Id, Dataset, ImgType, Slice, FileName)
  return(r)
}

getSliceList <- function(type = "sax")
{
  oneSliceGroup <- select(getImageList(type), Id, Dataset, ImgType, Offset, starts_with("Slice")) %>% 
    arrange(Dataset, Id, ImgType, Slice, SliceLocation) %>% unique()
  setkey(oneSliceGroup, Dataset, Id, ImgType, Slice, SliceLocation)
  return(oneSliceGroup)
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

