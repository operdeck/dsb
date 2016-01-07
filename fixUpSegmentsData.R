# Ad-hoc script to fix up segment data sets when formats
# change

# derived:
#       segmentInfo$roundness <- 4*pi*segmentInfo$s.area/(segmentInfo$s.perimeter^2)
# estimate of LV is area of ellipsoid with radius min and max of center (pi*a*b)
# slice thickness etc see Python code

source("util.R")

for (ds in datasetFolders) {
  f <- getSegmentFile(ds)
  if (file.exists(f)) {
    print(f)
    
    dataOld <- fread(f)
    
    ### FIX UP CODE HERE
    
    dataNew <- select(dataOld, -roundness)
    
    #####
    removedSet <- setdiff(names(dataOld), names(dataNew))
    addedSet <- setdiff(names(dataNew), names(dataOld))
    diffSet <- paste(c(paste("-",removedSet), paste("+",addedSet)),collapse=", ")
    
    if ((length(removedSet) + length(addedSet) > 0)) {
      print("Existing:")
      print(names(dataOld))
      print("New:")
      print(names(dataNew))
      print(paste(diffSet, collapse=","))
      
      # Be careful...
      write.csv(dataNew, f, row.names=F)
      
    } else {
      print("Identical")
    }    
    #   cat("Write final", getSegmentFile(ds), "id's:",length(unique(filter(allSegmentationInfo, Dataset==ds)$Id)), fill=T)
    #   write.csv(filter(allSegmentationInfo, Dataset==ds), getSegmentFile(ds), row.names=F)
  }
}
