library(deldir)
library(spatstat)
library(magrittr)
library(dplyr)
library(igraph)
library(scales)
library(httr)
library(tidyverse)
library(ggnetwork)
library(ggplot2)
library(poweRlaw)
library(imager)
library(viridis)
library(plotrix)
library(openxlsx)
library(tidyr)
library(spdep)
library(maptools)
library(tmap)
library(OpenImageR)
library(dismo)
library(xROI)
library(sp)
library(raster)
library(gpclib)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

folder_path = "D:/Spring 2022/Research/Quantifying similarities between spatial arrangements of the axons in peripheral nerves/Data/Non-reported in prev paper/"
file_name = "11328_5_AB-2V_2_montage_pred_mi"
extension = ".csv"

axon_locations = unique(read.csv(paste(folder_path, file_name, extension, sep="")))

#fascicle_contour = ripras(axon_locations$X, axon_locations$Y)
#axon_marks = unclass(axon_locations$Name) %>% as.factor

axon_pp = ppp(x = axon_locations$X, y = axon_locations$Y)
plot(axon_pp, main=file_name, cex=0.8, pch=21, bg="black")

####this part will need changes depending on the fascicle
if(interactive()){
  # r1= drawPolygon(col='black', lwd=2)
  # r2= drawPolygon(col='black', lwd=2)
  # 
  # r3= drawPolygon(col='black', lwd=2)
  # r4= drawPolygon(col='black', lwd=2)
  # r5= drawPolygon(col='black', lwd=2)
  # r6= drawPolygon(col='black', lwd=2)
  # 
  # r7= drawPolygon(col='black', lwd=2)
  r8= drawPolygon(col='black', lwd=2)
}

####polygons
# p1 = Polygon(r1, hole = TRUE)
# ps1 = Polygons(list(p1),ID = "a")
# ps1@Polygons[[1]]@hole = TRUE
# 
# p2 = Polygon(r2, hole = TRUE)
# ps2 = Polygons(list(p2),ID = "b")
# ps2@Polygons[[1]]@hole = TRUE
# 
# p3 = Polygon(r3, hole = TRUE)
# ps3 = Polygons(list(p3),ID = "c")
# ps3@Polygons[[1]]@hole = TRUE
# 
# p4 = Polygon(r4, hole = TRUE)
# ps4 = Polygons(list(p4),ID = "d")
# ps4@Polygons[[1]]@hole = TRUE
# 
# p5 = Polygon(r5, hole = TRUE)
# ps5 = Polygons(list(p5),ID = "e")
# ps5@Polygons[[1]]@hole = TRUE
# 
# p6 = Polygon(r6, hole = TRUE)
# ps6 = Polygons(list(p6),ID = "f")
# ps6@Polygons[[1]]@hole = TRUE
# 
# p7 = Polygon(r7, hole = TRUE)
# ps7 = Polygons(list(p7),ID = "g")
# ps7@Polygons[[1]]@hole = TRUE

p8 = Polygon(r8, hole = FALSE)
ps8 = Polygons(list(p8),ID = "h")
ps8@Polygons[[1]]@hole = FALSE

####spatial polygons
sps_m = SpatialPolygons(list(ps8))


fascicle_contour = as.owin(sps_m)
####save and retrieve contour
saveRDS(fascicle_contour, paste(folder_path, file_name, ".rds", sep=""))
retrieved_contour = readRDS(paste(folder_path, file_name, ".rds", sep=""))

axon_pp = ppp(x = axon_locations$X, y = axon_locations$Y, window = retrieved_contour)
plot(axon_pp, main=file_name, cex=0.8, pch=21, bg=c("black"))

# axon_tess = dirichlet(axon_pp)
# plot(axon_tess, main=file_name)
# 
# #check if discontinued contour works for point generation
# test_pp = rpoispp(lambda = 0.03, win = retrieved_contour)
# plot(test_pp)