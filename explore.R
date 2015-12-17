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

caseFolders <- list.dirs("data", recursive=F, full.names=F)
for (caseFolder in caseFolders) {
#   print(caseFolder) # = nr 1 .. 500
  caseFolder <- paste("data", caseFolder, "study", sep="/")
  
  serieFolders <- list.dirs(caseFolder, recursive=F, full.names=F)
  # NB excluding the 2 chamber and 4 chamber views here, only the short axis track
  serieFolders <- serieFolders[ grepl("sax_[[:digit:]]+$", serieFolders) ]
#   print(imgFolders)
  
  for (serieFolder in serieFolders) {
    sliceFolder <- paste(caseFolder, serieFolder, sep="/")

    print(sliceFolder) # inside here are the .dcm frames
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

folder <- 'data/10/study/sax_11'
filez <- list.files(folder, pattern=".*dcm$")

for (f in filez) {
  fname <- paste(folder, f, sep="/")
  
  dicom <- readDICOMFile(fname)
  img <- Image(normalize(dicom$img))
  img <- rotate(img,-90)
  # display(img, method = "raster")
  
  img_comb = combine(
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
  for (i in 1:max(segmented)) {
    
    coords <- as.data.frame(pos2coord(pos=which(segmented==i),dim.mat=dim(segmented)))
    names(coords) <- c('x','y')
    
    text(x = mean(coords$x), y = mean(coords$y), label = i, col = "orange")
    
    # for each segment some key characteristics
    # mid point, size, ?
    # then try match them from image to image
  }
}


