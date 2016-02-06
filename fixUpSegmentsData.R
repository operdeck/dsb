# Ad-hoc script to fix up segment data sets when formats
# change

#TODO
# segmentation data:
#   =================
#   
#   ** Keep
# 
# id's
# "Id"             
# "Slice"          
# "Time" 
# "segIndex"       
# "UUID"          
# 
# actual segment data
# 
# m.* s.*
# 
# ** Drop
# 
# redundant - all derivable from ids's via slice or imagelist
# "Dataset"         
# "ImgType"        
# "FileName"     
# 
# reduntant - all coming from imagelist  
# "Offset"        
# "SliceCount"     
# "SliceIndex"     
# "SliceOrder"     
# "PixelSpacing.x" 
# "PixelSpacing.y" 
# "SliceLocation"  
# "SliceThickness"
# 
# "isProcessed"
# "distToROI"      
# 
# ** Add (derive from distToROI, m.cx and m.cy per image)
# ROI.x, ROI.y, ROI.r (int)
# 
# ** Change
# "m.cx"  (int)         
# "m.cy"  (int)
# m.majoraxis (int)
# s.area (int)
# s.perimeter (int)


source("util.R")
imageList <- getImageList()

fixupROI <- function(imageSegments) {
  n <- nrow(imageSegments)
  if (n > 1) {
    rA <- imageSegments$distToROI[1:(nrow(imageSegments)-1)]
    rB <- imageSegments$distToROI[2:nrow(imageSegments)]
    xA <- imageSegments$m.cx[1:(nrow(imageSegments)-1)]
    xB <- imageSegments$m.cx[2:nrow(imageSegments)]
    yA <- imageSegments$m.cy[1:(nrow(imageSegments)-1)]
    yB <- imageSegments$m.cy[2:nrow(imageSegments)]
    
    d2 <- (xB-xA)^2 + (yB-yA)^2 
    K  <- sqrt(((rA+rB)^2-d2)*(d2-(rA-rB)^2))/4
    f  <- (rA^2-rB^2)/d2/2
    
    x1 <- round((xB+xA)/2 + (xB-xA)*f + 2*(yB-yA)*K/d2,2)
    y1 <- round((yB+yA)/2 + (yB-yA)*f - 2*(xB-xA)*K/d2,2)
    x2 <- round((xB+xA)/2 + (xB-xA)*f - 2*(yB-yA)*K/d2,2)
    y2 <- round((yB+yA)/2 + (yB-yA)*f + 2*(xB-xA)*K/d2,2)
    
    if (n > 2) {
      # get most frequent X/Y value
      roi.x <- round(as.double(names(sort(table(c(x1,x2)),decreasing=T))[1]))
      roi.y <- round(as.double(names(sort(table(c(y1,y2)),decreasing=T))[1]))
      roi.r <- round(max(imageSegments$distToROI) + 5 )
    } else {
      roi.x <- round(mean(c(x1,x2)))
      roi.y <- round(mean(c(y1,y2)))
      roi.r <- round(max(imageSegments$distToROI) + 5 )
    }
  } else {
    roi.x <- round(imageSegments$m.cx[1])
    roi.y <- round(imageSegments$m.cy[1])
    roi.r <- round(imageSegments$distToROI[1] + 5)
  }  
  cat("ROI=",roi.x,",",roi.y," r=",roi.r,"n=",nrow(imageSegments),fill=T)
  
  dataNew[Id == imageSegments$Id & Slice == imageSegments$Slice & Time == imageSegments$Time, ROI.x := roi.x]
  dataNew[Id == imageSegments$Id & Slice == imageSegments$Slice & Time == imageSegments$Time, ROI.y := roi.y]
  dataNew[Id == imageSegments$Id & Slice == imageSegments$Slice & Time == imageSegments$Time, ROI.r := roi.r]
  
}

for (ds in unique(imageList$Dataset)) {
  f <- getSegmentFile(ds)
  if (file.exists(f)) {
    print(f)
    
    dataOld <- fread(f)
    
    ### FIX UP CODE HERE

    dataNew <- select(dataOld, -Dataset)
    
    #     dataNew <- select(dataOld, 
    #                       -Dataset, -ImgType, -FileName,
    #                       -Offset, -SliceCount, -SliceIndex, -SliceOrder, -PixelSpacing.x, -PixelSpacing.y,
    #                       -SliceLocation, -SliceThickness,
    #                       -isProcessed)
    #     dataNew$m.cx <- round(dataOld$m.cx)
    #     dataNew$m.cy <- round(dataOld$m.cy)
    #     dataNew$m.majoraxis <- round(dataOld$m.majoraxis)
    #     dataNew$s.area <- round(dataOld$s.area)
    #     dataNew$s.perimeter <- round(dataOld$s.perimeter)
    
    #     imgList <- unique(select(dataOld, Id, Slice, Time))
    #     for (i in 1:nrow(imgList)) {
    #       cat(i, "of", nrow(imgList), " %", round(100*i/nrow(imgList),2), " ")
    #       ds <- filter(dataOld, Id == imgList$Id[i], Slice == imgList$Slice[i], Time == imgList$Time[i])
    #       fixupROI(ds)
    #     }
    
    #dataNew <- select(dataOld, -Dataset.x, -Dataset.y, -Dataset)
    
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
