# ROI x/y from points A, B and dist to ROI

# x = (1/2)(xB+xA) + (1/2)(xB-xA)(A$distToROI^2-B$distToROI^2)/d2 ± 2(yB-yA)K/d2 
# y = (1/2)(yB+yA) + (1/2)(yB-yA)(A$distToROI^2-B$distToROI^2)/d2 ± -2(xB-xA)K/d2 
# d = distance A-B

# for one image (vectorized):
# see http://2000clicks.com/mathhelp/GeometryConicSectionCircleIntersection.aspx

rA <- imageSegments$distToROI[1:(nrow(imageSegments)-1)]
rB <- imageSegments$distToROI[2:nrow(imageSegments)]
xA <- imageSegments$m.cx[1:(nrow(imageSegments)-1)]
xB <- imageSegments$m.cx[2:nrow(imageSegments)]
yA <- imageSegments$m.cy[1:(nrow(imageSegments)-1)]
yB <- imageSegments$m.cy[2:nrow(imageSegments)]
  
d2 <- (xB-xA)^2 + (yB-yA)^2 
K <- sqrt(((rA+rB)^2-d2)*(d2-(rA-rB)^2))/4
  
x1 <- round((xB+xA)/2 + (xB-xA)*(rA^2-rB^2)/d2/2 + 2*(yB-yA)*K/d2,2)
y1 <- round((yB+yA)/2 + (yB-yA)*(rA^2-rB^2)/d2/2 - 2*(xB-xA)*K/d2,2)
x2 <- round((xB+xA)/2 + (xB-xA)*(rA^2-rB^2)/d2/2 - 2*(yB-yA)*K/d2,2)
y2 <- round((yB+yA)/2 + (yB-yA)*(rA^2-rB^2)/d2/2 + 2*(xB-xA)*K/d2,2)
  
# get most frequent X/Y value
roi.x <- as.double(names(sort(table(c(x1,x2)),decreasing=T))[1])
roi.y <- as.double(names(sort(table(c(y1,y2)),decreasing=T))[1])
roi.r <- max(imageSegments$distToROI)  

slice$roi.x <- roi.x
slice$roi.y <- roi.y
slice$roi.r <- roi.r

slice$roi.dist <- sqrt((slice$m.cx-roi.x)^2 + (slice$m.cy-roi.y)^2)


