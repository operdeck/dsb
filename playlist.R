# List all images and save in a .csv file. List contains DICOM meta info
# and should be run prior to any of the other processing.
# Running it can take up ca 6 hours for 700 id's

source("Util.R")

# Now contains all images, also the 3 pair .dcm. Default filters on 'sax'.
# * using SliceIndex instead of SliceRelIndex
# * adds FileName
# * removed slicelist.csv from git, added imagelist.csv
# Fills in slice thickness etc if they're missing by mean
# TODO consider extracting patient info from header - possibly interesting for predictions
#  (0010,0040) Patient's Sex [M]
#  (0010,1010) Patient's Age [012Y]

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
            c("SliceCount", "SliceIndex", "SliceOrder") := list(sum(One), 
                                                                seq(sum(One)), 
                                                                midOrderSeqKeepSibblingsTogether(sum(One))),
            by=c("Dataset", "Id", "ImgType")]
playlist <- left_join(playlist, sliceMetaData, by=c("Dataset", "Id", "ImgType", "Slice")) %>% select(-One)

print("Reading image data")
for (n in seq(nrow(playlist))) {
  img <- playlist[n,]
  print(img$FileName)
  dicomHeader <- readDICOMFile(paste("data", img$FileName, sep="/"), pixelData = FALSE)
  pixelSpacing <- extractHeader(dicomHeader$hdr, "PixelSpacing", numeric=FALSE)

  playlist[n, "PixelSpacing.x" := as.numeric(gsub("^(.*) (.*)$", "\\1", pixelSpacing))]
  playlist[n, "PixelSpacing.y" := as.numeric(gsub("^(.*) (.*)$", "\\2", pixelSpacing))]
  playlist[n, "SliceLocation"  := extractHeader(dicomHeader$hdr, "SliceLocation", numeric=TRUE)]
  playlist[n, "SliceThickness" := extractHeader(dicomHeader$hdr, "SliceThickness", numeric=TRUE)]
}
print("Fill missing")
cat("Incomplete cases:",sum(!complete.cases(playlist)),"of",nrow(playlist),fill=T)
# Fill in missings. Maybe should do a (tree) model instead of just mean.
playlist[is.na(PixelSpacing.x), PixelSpacing.x := mean(PixelSpacing.x, na.rm=T)]
playlist[is.na(PixelSpacing.y), PixelSpacing.y := mean(PixelSpacing.y, na.rm=T)]
playlist[is.na(SliceLocation),  SliceLocation  := mean(SliceLocation, na.rm=T)]
playlist[is.na(SliceThickness), SliceThickness := mean(SliceThickness, na.rm=T)]

write.csv(playlist, "imagelist.csv", row.names=F)

# Show quick summary of the datasets & save for downstream use
print("Results:")
print(select(getIdList(), Dataset, Id) %>% group_by(Dataset) %>% summarise(n_Ids = n()))
print(select(getSliceList(), Dataset, Id, Slice) %>% group_by(Dataset) %>% summarise(n_Slices = n()))
print(select(getImageList(), Dataset, Id, Slice, FileName) %>% group_by(Dataset) %>% summarise(n_Images = n()))

