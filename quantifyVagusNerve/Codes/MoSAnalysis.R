#### loading required libraries (there might be more libraries loaded than required)
load.lib = c("deldir", "spatstat", "spatstat.utils", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "oce", "OneR", "RandomFieldsUtils", "RandomFields", "Cairo", "knitr", "scorepeak", "Rcpp", "emdist",
             "RImageJROI", "svglite", "transport", "Barycenter", "T4transport", "wvtool", "adimpro", "reshape2", "proxy", "RColorBrewer", "tictoc",
             "ggrepel", "scatterplot3d", "car", "e1071", "rgl", "this.path")

install.lib = load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(pkgs=lib, dependencies=TRUE)
sapply(load.lib, require, character=TRUE)


version = "complete"


#### extracting parent directory information for accessing input and output location
dir = this.dir()
folder = strsplit(dir, "/")
folder = folder[[1]][length(folder[[1]])]
parent = strsplit(dir, folder)

#### output plot location
plot_dest_path = paste(parent, "Plots/Moment of Shape/", sep="")
if (!dir.exists(plot_dest_path)) {dir.create(plot_dest_path, recursive=TRUE)}


#### source the inhomogeneous sector functions
source("LocalKSector.R")
source("Kinhomsector.R")


#### input location
if(version == "complete"){
  folder_path = paste(parent, "Data/Inputs/", sep = "")
}else if(version == "demo"){
  folder_path = paste(parent, "Data/Demo Inputs/", sep = "")
}else{
  print("Invalid version...\n")
  return()
}
data_files = list.files(path=folder_path, pattern = ".csv", full.names = FALSE)
extension = ".csv"

data_files = data_files[c(15, 18, 23, 28)]
#data_files = data_files[c(15, 18)]

df = data.frame()

#### point pattern and plots, comment out the while loop block if images are not required
i = 1
while (i <= length(data_files)) {
  
  #construct the i-th point pattern
  file_name = strsplit(data_files[i], ".csv")[[1]]
  print(file_name)
  
  axon_locations = unique(read.csv(paste(folder_path, file_name, extension, sep="")))
  retrieved_contour = readRDS(paste(folder_path, file_name, ".rds", sep=""))
  
  axon_marks = data.frame(Circularity = axon_locations$Circ., 
                          Extension = axon_locations$Ext., Dispersion = axon_locations$Disp., Elongation = axon_locations$Elong.)
  
  axon_pp = ppp(x=axon_locations$X, y=axon_locations$Y, marks = axon_marks,
                checkdup=F, window = retrieved_contour)
  
  axon_pp = rescale.ppp(axon_pp, s=320153.4) 
  axon_pp = shift.ppp(axon_pp, origin = "centroid")
  
  axon_pp = as.ppp(axon_pp)
  svglite(paste(plot_dest_path, file_name, "_pp.svg", sep=""), width = 10, height = 10)
  plot(axon_pp)
  dev.off()
  
  # Sample_Name = rep(file_name, axon_pp$n)
  # Sample_ID = rep(i, axon_pp$n)
  # temp = cbind(axon_pp$marks, Sample_Name, Sample_ID)
  # 
  # df = rbind(df, temp)
  
  dens = Smooth(axon_pp, sigma=bw.diggle(axon_pp), adjust=2, diggle=TRUE, eps=0.0025)
    
  jpeg(paste(plot_dest_path, file_name,"_Moment_of_Shape.jpeg", sep=""),
    width = 10,
    height = 10,
    units = "in",
    res=900)
  plot(dens, col=colorRampPalette(brewer.pal(9, "YlOrRd"))(9), box=FALSE)
  dev.off()
  
  ####priciple component analysis of the marks of the point pattern
  mark_pca = prcomp(axon_pp$marks, center = TRUE, scale. = TRUE, retx = TRUE)
  print(summary(mark_pca))
  
  pca_1 = mark_pca$x[, 1]
  pca_2 = mark_pca$x[, 2]

  ####retrieve the first and second principle component
  # pca_1 = multiply_by_matrix(as.matrix(axon_pp$marks, byrow=TRUE), as.vector(mark_pca$rotation[, 1]))
  # pca_2 = multiply_by_matrix(as.matrix(axon_pp$marks, byrow=TRUE), as.vector(mark_pca$rotation[, 2]))
  
  Sample_Name = rep(file_name, axon_pp$n)
  Sample_ID = rep(i, axon_pp$n)
  temp = data.frame(PC1=pca_1, PC2=pca_2, Sample_Name=Sample_Name, Sample_ID=Sample_ID)

  df = rbind(df, temp)
  
  # pca1.binned <- bin(pca1, nbins=3, labels= c("low", "medium", "high"), method="cluster")
  # 
  # ####create new version of point pattern with a specific principle component as feature
  # ganglia.ppp.2 <- ppp(x=ganglia$xC, y=ganglia$yC, window=win, marks = pca1.binned)
  # plot(Smooth.ppp(ganglia.ppp.2))
  # 
  # 
  # # colorset <- colourmap(c("yellow", "orange", "dark red"), breaks = c(-25000, -10000, -5000, 0))
  # # plot(ganglia.ppp.2, pch=21, bg=colorset, cex=1.5)
  # 
  # #colorset <- colourmap(c("yellow", "orange", "dark red"), breaks = c(-2.58e+04, -1.52e+04, -5.89e+03, -438))
  # 
  # plot(ganglia.ppp.2)
  
  
  i = i + 1
}


#### visualize
svglite(paste(plot_dest_path, "PCA.svg", sep=""), width = 9, height = 7)

ggplot(data = df, aes(x = PC1, y = PC2)) +
    geom_point(size = 1.5, shape=21, colour="grey40", aes(fill = Sample_Name)) +
    
    # geom_text_repel(size=2.5, 
    #                 segment.size=0.2, segment.color="grey40",
    #                 arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "first"),
    #                 max.overlaps = Inf) +
    
    #stat_ellipse(type="t", aes(colour=Nerve)) +
    
    theme_bw()+
    theme(legend.position="right",  legend.title=element_blank(), legend.text=element_text(size=10), 
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=12),
          plot.subtitle = element_text(hjust = 0.5, size=10),
          axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", size=0.25, linetype=2)) +
    
    xlab(expression(paste("PC 1"))) + ylab("PC 2")+
    
    labs(title = "PCA of the Moments of Shape",    # the titles needs changing for different runs
         subtitle = "Circularity, Extension, Dispersion, Elongation")

dev.off()


