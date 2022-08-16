#### loading required libraries (there might be more libraries loaded than required)
load.lib = c("deldir", "spatstat", "spatstat.utils", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "oce", "OneR", "RandomFieldsUtils", "RandomFields", "Cairo", "knitr", "scorepeak", "Rcpp", "emdist",
             "RImageJROI", "svglite", "transport", "Barycenter", "T4transport", "wvtool", "adimpro", "reshape2", "proxy", "RColorBrewer", "tictoc",
             "ggrepel", "scatterplot3d", "car", "e1071", "rgl", "this.path")

install.lib = load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(pkgs=lib, dependencies=TRUE)
sapply(load.lib, require, character=TRUE)


#### command line arguments contain analysis type, scaling information and entropy regularizer lambda
args = commandArgs()
print(args)

#### the number of arguments to check for should be changed if required
#### 11 arguments are used while using Rscript command from RStudio terminal
#### 7 arguments are required if run with slurm job
if(length(args) != 11){
  print("Invalid arguments...\n")
  return()
} else{
  analysis_type = args[7]
  scaling = as.numeric(args[8])
  lambda  = as.numeric(args[9])
  version = args[10]
  intr_dist = as.numeric(args[11])
}

print(analysis_type)
print(scaling)
print(lambda)
print(version)
print(intr_dist)


#### extracting parent directory information for accessing input and output location
dir = this.dir()
folder = strsplit(dir, "/")
folder = folder[[1]][length(folder[[1]])]
parent = strsplit(dir, folder)

#### output plot location
plot_dest_path = paste(parent, "Plots/T4transport/", version, "_", analysis_type, scaling, "_", lambda, "/", sep="")
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


#### point pattern and plots, comment out the while loop block if images are not required
i = 1
while (i <= length(data_files)) {
  
  #construct the i-th point pattern
  file_name = strsplit(data_files[i], ".csv")[[1]]
  print(file_name)
  
  axon_locations = unique(read.csv(paste(folder_path, file_name, extension, sep="")))
  retrieved_contour = readRDS(paste(folder_path, file_name, ".rds", sep=""))
  
  axon_pp = ppp(x=axon_locations$X, y=axon_locations$Y, checkdup=F, window = retrieved_contour)
  
  if(scaling == 1){
    print("scaled")
    axon_pp = rescale.ppp(axon_pp, s=axon_pp$window$yrange[2]) 
  }
  axon_pp = shift.ppp(axon_pp, origin = "centroid")
  
  axon_pp = as.ppp(axon_pp)
  svglite(paste(plot_dest_path, file_name, "_pp.svg", sep=""), width = 10, height = 10)
  plot(axon_pp, pch=21, bg="black", cex=0.7, main="")
  dev.off()
  
  if(analysis_type=="basic_density"){
    #basic density
    dens = density(axon_pp, sigma=bw.diggle(axon_pp), adjust = 2, diggle=TRUE, eps=0.0025)
    dens_mat = t(dens$v)
    
    jpeg(paste(plot_dest_path, file_name,"_Density.jpeg", sep=""),
         width = 10,
         height = 10,
         units = "in",
         res=900)
    plot(dens, main="", col=colorRampPalette(brewer.pal(9, "YlOrRd"))(9), box=FALSE)
    dev.off()
    
  }else if(analysis_type=="local_linhom"){
    pp_l = localLinhom(axon_pp, rvalue = intr_dist, verbose=FALSE)
    marks(axon_pp) = pp_l
    dens = Smooth(axon_pp, sigma=bw.diggle(axon_pp), adjust=2, diggle=TRUE, eps=0.0025)
    
    jpeg(paste(plot_dest_path, file_name,"_Local_Linhom.jpeg", sep=""),
         width = 10,
         height = 10,
         units = "in",
         res=900)
    plot(dens, main="", col=colorRampPalette(brewer.pal(9, "YlOrRd"))(9), box=FALSE)
    dev.off()
    
  }else if(analysis_type=="local_linhom_sector_horizontal"){
    pp_l = localLinhomsector(axon_pp, rvalue = intr_dist, begin=-15, end=15, units="degrees", verbose=FALSE)
    marks(axon_pp) = pp_l
    dens = Smooth(axon_pp, sigma=bw.diggle(axon_pp), adjust=2, diggle=TRUE, eps=0.0025)
    
    jpeg(paste(plot_dest_path, file_name,"_Local_Linhom_H.jpeg", sep=""),
         width = 10,
         height = 10,
         units = "in",
         res=900)
    plot(dens, main="", col=colorRampPalette(brewer.pal(9, "YlOrRd"))(9), box=FALSE)
    dev.off()
    
  }else if(analysis_type=="local_linhom_sector_vertical"){
    pp_l = localLinhomsector(axon_pp, rvalue = intr_dist, begin=90-15, end=90+15, units="degrees", verbose=FALSE)
    marks(axon_pp) = pp_l
    dens = Smooth(axon_pp, sigma=bw.diggle(axon_pp), adjust=2, diggle=TRUE, eps=0.0025)
    
    jpeg(paste(plot_dest_path, file_name,"_Local_Linhom_V.jpeg", sep=""),
         width = 10,
         height = 10,
         units = "in",
         res=900)
    plot(dens, main="", col=colorRampPalette(brewer.pal(9, "YlOrRd"))(9), box=FALSE)
    dev.off()
  }else{
    print("Invalid analysis type!")
  }
  
  i = i + 1
}


#### creating workbook to save results
wb = openxlsx::createWorkbook()
distance_matrix = matrix(0, nrow = length(data_files), ncol = length(data_files))
min_angle_matrix = matrix(0, nrow = length(data_files), ncol = length(data_files))

#### read the rotation of angle matrix previously determined by basic spatial density analysis
m = openxlsx::read.xlsx(xlsxFile = paste(parent, "Data/Supporting Files/Angle of Rotation/Angle_of_rotation_T4transport_", version, "_", scaling, "_", lambda,".xlsx", sep=""), 
                        sheet = "Min angle matrix", colNames = FALSE)
m = as.matrix(m)


#### computing pair-wise Sinkhorn distance
i = 1
while (i <= length(data_files)) {

  j = i+1
  while (j <= length(data_files)) {
    cat("Fascicle pair: ", i, " ", j, "\n")
    
    #construct the i-th point pattern
    file_name_1 = strsplit(data_files[i], ".csv")[[1]]
    print(file_name_1)
    
    axon_locations_1 = unique(read.csv(paste(folder_path, file_name_1, extension, sep="")))
    retrieved_contour_1 = readRDS(paste(folder_path, file_name_1, ".rds", sep=""))
    
    axon_pp_1 = ppp(x=axon_locations_1$X, y=axon_locations_1$Y, checkdup=F, window = retrieved_contour_1)
    if(scaling == 1){
    	print("scaled")
    	axon_pp_1 = rescale.ppp(axon_pp_1, s=axon_pp_1$window$yrange[2])
    }
    axon_pp_1 = shift.ppp(axon_pp_1, origin = "centroid") # making translation invariant
    
    axon_pp_1 = as.ppp(axon_pp_1)
    
    #construct the j-th point pattern
    file_name_2 = strsplit(data_files[j], ".csv")[[1]]
    print(file_name_2)
    
    axon_locations_2 = unique(read.csv(paste(folder_path, file_name_2, extension, sep="")))
    retrieved_contour_2 = readRDS(paste(folder_path, file_name_2, ".rds", sep=""))
    
    axon_pp_2 = ppp(x=axon_locations_2$X, y=axon_locations_2$Y, checkdup=F, window = retrieved_contour_2)
    
    if(scaling == 1){
    	print("scaled")
    	axon_pp_2 = rescale.ppp(axon_pp_2, s=axon_pp_2$window$yrange[2])
    }
    axon_pp_2 = shift.ppp(axon_pp_2, origin = "centroid") # making translation invariant
    
    # angle of rotation from the matrix
    theta = m[i, j]
    cat("theta: ", theta, "\n")
      
  	axon_pp_1_rotated = rotate.ppp(axon_pp_1, angle = theta*pi/180, centre = "centroid")
  	svglite(paste(plot_dest_path, paste(file_name_1, theta, "deg rotated", sep = " "), "_pp.svg", sep=""), width = 10, height = 10)
  	plot(axon_pp_1_rotated, pch=21, bg="black", cex=0.7, main=paste(file_name_1, theta, "deg rotated", sep = " "))
  	dev.off()

  	X = as.matrix(data.frame(axon_pp_1_rotated$x, axon_pp_1_rotated$y))
  	Y = as.matrix(data.frame(axon_pp_2$x, axon_pp_2$y))

  	if(analysis_type=="basic_density"){
  		WX = rep(1/axon_pp_1_rotated$n, axon_pp_1_rotated$n)
  		WY = rep(1/axon_pp_2$n, axon_pp_2$n)
  
  	}else if(analysis_type=="local_linhom"){
  		pp1_l = localLinhom(axon_pp_1_rotated, rvalue = intr_dist, verbose=FALSE)
  		pp2_l = localLinhom(axon_pp_2, rvalue = intr_dist, verbose=FALSE)
  
  		WX = pp1_l / sum(pp1_l)
  		WY = pp2_l / sum(pp2_l)
  
  	}else if(analysis_type=="local_linhom_sector_horizontal"){
  		pp1_l = localLinhomsector(axon_pp_1_rotated, rvalue = intr_dist, begin=-15, end=15, units="degrees", verbose=FALSE)
  		pp2_l = localLinhomsector(axon_pp_2, rvalue = intr_dist, begin=-15, end=15, units="degrees", verbose=FALSE)
  
  		WX = pp1_l / sum(pp1_l)
  		WY = pp2_l / sum(pp2_l)
  	 
  	}else if(analysis_type=="local_linhom_sector_vertical"){
  		pp1_l = localLinhomsector(axon_pp_1_rotated, rvalue = intr_dist, begin=90-15, end=90+15, units="degrees", verbose=FALSE)
  		pp2_l = localLinhomsector(axon_pp_2, rvalue = intr_dist, begin=90-15, end=90+15, units="degrees", verbose=FALSE)
  
  		WX = pp1_l / sum(pp1_l)
  		WY = pp2_l / sum(pp2_l)
  
  	}else{
  		print("Invalid analysis type!")
  	}

  	#cat("wx:\n", WX, "\nsum wx: ", sum(WX), "\n")
  	#cat("wy:\n", WY, "\nsum wy: ", sum(WY), "\n")

  	#using the function from T4transport package
  	if(lambda == -1){
  		print("Computing wasserstein distance...")
  		tic("Wasserstein")
  		dist_t4 = T4transport::wasserstein(X, Y)
  		toc()
  		cat(dist_t4$distance, "\n")
  	}
  	else{
  		print("Computing sinkhorn distance...")
  		tic(paste("Sinkhorn lambda = ", lambda, sep=""))
  		dist_t4 = T4transport::sinkhorn(X, Y, lambda = lambda, wx=WX, wy=WY) #default regularization param (lambda) is 0.1.
  		toc()
  		cat(dist_t4$distance, "\n")
  	}

	  distance_matrix[i, j] = dist_t4$distance
    min_angle_matrix[i, j] = theta
    
    j = j+1
  }
  
  cat("\n")
  i=i+1
}

addWorksheet(wb, sheetName = "Distance matrix")
writeData(wb, sheet = "Distance matrix", x=distance_matrix, rowNames = FALSE, colNames = FALSE)

addWorksheet(wb, sheetName = "Min angle matrix")
writeData(wb, sheet = "Min angle matrix", x=min_angle_matrix, rowNames = FALSE, colNames = FALSE)

openxlsx::saveWorkbook(wb, file =paste(plot_dest_path, version, "_", analysis_type, scaling, "_", lambda,"_all_result.xlsx", sep=""), overwrite = TRUE)
  
  
  

