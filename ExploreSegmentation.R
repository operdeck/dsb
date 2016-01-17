source("Util.R")

for (t in seq(from=1, to=30, by=2)) {
f <- getImageFile(data.frame(Dataset="train", Id=200, ImgType="sax", Slice=44, Offset=6466, Time=t))
dicom <- readDICOMFile(f)
img <- Image(normalize(dicom$img))
if (dim(img)[1] > dim(img)[2]) {
  img <- rotate(img,-90)
}
kern <- makeBrush(5, shape= 'disc')
img_blurred <- normalize(opening(img, kern))
img_thresholded <- img_blurred > otsu(img_blurred) # Otsu??s threshold 
img_segmented <- fillHull(bwlabel(img_thresholded))

# img_comb <- EBImage::combine(img, img_blurred, img_thresholded, 
#                              thresh(img),fillHull(opening(thresh(img))), 
#                              img_segmented)
# display(img_comb)

j <- closing(fillHull(opening(thresh(img))))
display(EBImage::combine(toRGB(img), 
                         toRGB(j), 
                         colorLabels(bwlabel(j)), 
                         toRGB(img_segmented),
                         colorLabels(propagate(img, bwlabel(j), lambda=1))), 
        method="raster", all=T)
}
