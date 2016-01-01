# Manual identification of LV segments. Creates "segments-classified-<dataset>.csv" files with
# same columns as "segments-<dataset>.csv" but with additional "isLeftVentricle" column that
# can be used for prediciton of LV segment in full dataset.

source("util.R")


showSlice <- function(ds, highlight=F) {
#   filez <-  ds$file
#   filez <-  unique(ds$file)
  imgz <- list()
  for (i in 1:nrow(ds)) {
    f <- ds$file[i]
    dicom <- readDICOMFile(f)
    img <- Image(normalize(dicom$img))
    if (dim(img)[1] > dim(img)[2]) {
      img <- rotate(img,-90)
    }
    if (highlight) {
      imgz[[f]] <- drawCircle(toRGB(img), 
                              x=ds$m.cx[i], y=ds$m.cy[i], radius=ds$s.radius.mean[i], col="red")
    } else {
      imgz[[f]] <- img
    }
  }
  EBImage::display(EBImage::combine(imgz),all=T,method="raster")
  text(500,40,paste("Slice",ds$Slice[1]),col="yellow",pos=4) # TODO show Id + slice
  
  # show segments in 1st image
  firstImage <- filter(ds, file==ds$file[1])
  for (i in 1:nrow(firstImage)) {
    text(x = firstImage$m.cx[i], y = firstImage$m.cy[i], 
         label = firstImage$segIndex[i], col = "red")
  }
}


# Read all segment info from all datasets
allSegInfo <- NULL
for (dataset in c('train','validate','test')) {
  fname <- paste("segments-",dataset,".csv",sep="")
  if (file.exists(fname)) {
    segInfo <- fread(fname)
  }
  if (is.null(allSegInfo)) {
    allSegInfo <- segInfo
  } else {
    allSegInfo <- rbind(allSegInfo, segInfo)
  }
}

# Read previous segment classification
segPredictSet <- NULL
fname <- "segments-predict.csv"
if (file.exists(fname)) {
  segPredictSet <- fread(fname)
}

# Find UUID's not classified yet
if (! is.null(segPredictSet)) {
  unclassified <- setdiff(allSegInfo$UUID, segPredictSet$UUID)
} else {
  unclassified <- allSegInfo$UUID
}

# Select only the middle segment for each Id (TODO: consider broader selection)
allSegInfo <- left_join(filter(allSegInfo, UUID %in% unclassified), 
                        group_by(allSegInfo, Id) %>% summarise(midSlice = sort(unique(Slice))[trunc(length(unique(Slice))/2)])) %>%
  filter(Slice == midSlice)

# Get input for these slices
inputSlices <- unique(select(allSegInfo, Id, Slice)) 
for (s in 1:nrow(inputSlices)) {
  slice <- filter(allSegInfo, Id == inputSlices$Id[s], Slice == inputSlices$Slice[s])
  setkey(slice) # drop keys
  slice <- unique(slice) # TODO not sure why there are duplicates
  
  showSlice(slice)
  
  lvSeg <- readInteger("Segment of left ventricle in first image (0=none): ")
  filez <- unique(slice$file)
  slice <- slice[file == filez[1], isLV := (segIndex == lvSeg)]
  lvSeg.x <- slice$m.cx[slice$file == filez[1] & slice$isLV]
  lvSeg.y <- slice$m.cy[slice$file == filez[1] & slice$isLV]
  slice$distToLVSeg <- sqrt((slice$m.cx - lvSeg.x)^2 + (slice$m.cy - lvSeg.y)^2)
  lvSegs <- group_by(slice, file) %>% summarise(lvSegGuessedIndex = segIndex[which.min(distToLVSeg)])
  slice <- left_join(slice, lvSegs)
  
  guessedSegsInAllImages <- 
    select(filter(slice, segIndex == lvSegGuessedIndex), file, m.cx, m.cy, s.radius.mean, segIndex)
  showSlice(guessedSegsInAllImages, highlight=T)
  
  # TODO : prompt if all are ok
  # fix incorrect ones
  # set isLV
  # save data
}


