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
plot_dest_path = paste(parent, "Plots/Barycenter/", version, "_", analysis_type, scaling, "_", lambda, "/", sep="")
if (!dir.exists(plot_dest_path)) {dir.create(plot_dest_path, recursive=TRUE)}


#### source the inhomogeneous sector functions
source("LocalKSector.R")
source("Kinhomsector.R")


#### auxiliary function, takes a matrix as input, replaces any NA values and normalizes the matrix to sum to one.
processMap1 <- function(X){
  X[is.na(X)] = 10^(-20) #very small number close to 0, but not 0.
  X = X / sum(X)
  
  M = nrow(X)
  P = ncol(X)
  cat(M, "*", P, "\n")
  
  return(X)
}


#### rvalue=0.0373 was for 29 fascicles, one of them has quite different spatial arrangement, 
#### hence different Linhom characteristics. we eliminated that one in rvalue
#### computation and came up with another rvalue=0.0139.

#### analysis of variance gives 0.032 (29 fascicles), 0.0352 (28 fascicles). we are using this one currently.
constructMap <- function(pp1, pp2, type){
  cat(type, "\n")
  # the eps value determines the number of columns in the density map matrix
  
  if(type=="basic_density"){
    pp_density_map_1 = density(pp1, sigma=bw.diggle(pp1), adjust=2, diggle=TRUE, eps=0.01) 
    X_1 = t(pp_density_map_1$v)
    
    pp_density_map_2 = density(pp2, sigma=bw.diggle(pp2), adjust=2, diggle=TRUE, eps=0.01)
    X_2 = t(pp_density_map_2$v)
    
    D_1 = pp_density_map_1
    D_2 = pp_density_map_2
    
  }else if(type=="local_linhom"){
    pp1_l = localLinhom(pp1, rvalue = intr_dist, verbose=FALSE)
    marks(pp1) = pp1_l
    
    D_1 = Smooth(pp1, sigma=bw.diggle(pp1), adjust=2, diggle=TRUE, eps=0.01)
    X_1 = t((D_1)$v)
    
    pp2_l = localLinhom(pp2, rvalue = intr_dist, verbose=FALSE)
    marks(pp2) = pp2_l
    
    D_2 = Smooth(pp2, sigma=bw.diggle(pp2), adjust=2, diggle=TRUE, eps=0.01)
    X_2 = t((D_2)$v)
    
  }else if(type=="local_linhom_sector_horizontal"){
    pp1_l = localLinhomsector(pp1, rvalue = intr_dist, begin=-15, end=15, units="degrees", verbose=FALSE)
    marks(pp1) = pp1_l
    
    D_1 = Smooth(pp1, sigma=bw.diggle(pp1), adjust=2, diggle=TRUE, eps=0.01)
    X_1 = t((D_1)$v)
    
    pp2_l = localLinhomsector(pp2, rvalue = intr_dist, begin=-15, end=15, units="degrees", verbose=FALSE)
    marks(pp2) = pp2_l
    
    D_2 = Smooth(pp2, sigma=bw.diggle(pp2),  adjust=2, diggle=TRUE, eps=0.01)
    X_2 = t((D_2)$v)
    
  }else if(type=="local_linhom_sector_vertical"){
    pp1_l = localLinhomsector(pp1, rvalue = intr_dist, begin=90-15, end=90+15, units="degrees", verbose=FALSE)
    marks(pp1) = pp1_l
    
    D_1 = Smooth(pp1, sigma=bw.diggle(pp1),adjust=2,diggle=TRUE, eps=0.01)
    X_1 = t((D_1)$v)
    
    pp2_l = localLinhomsector(pp2, rvalue = intr_dist, begin=90-15, end=90+15, units="degrees", verbose=FALSE)
    marks(pp2) = pp2_l
    
    D_2 = Smooth(pp2, sigma=bw.diggle(pp2), adjust=2, diggle=TRUE, eps=0.01)
    X_2 = t((D_2)$v)
    
  }else{
    print("Invalid argument!")
    return(list(NULL, NULL, NULL, NULL))
  }
  
  X_1 = processMap1(X_1)
  X_2 = processMap1(X_2)
  
  r_1 = nrow(X_1)
  c_1 = ncol(X_1)
  
  r_2 = nrow(X_2)
  c_2 = ncol(X_2)
  
  # 0-padding for making matrix equal size
  if(r_1 > r_2){
    X_2 = (padding(as.matrix(X_2), new_rows = r_1, new_cols = ncol(X_2), fill_value = 10^(-20)))$data
  }else{
    X_1 = (padding(as.matrix(X_1), new_rows = r_2, new_cols = ncol(X_1), fill_value = 10^(-20)))$data
  }
  
  if(c_1 > c_2){
    X_2 = (padding(as.matrix(X_2), new_rows = nrow(X_2), new_cols = c_1, fill_value = 10^(-20)))$data
  }else{
    X_1 = (padding(as.matrix(X_1), new_rows = nrow(X_1), new_cols = c_2, fill_value = 10^(-20)))$data
  }
  
  return(list(X_1, X_2, D_1, D_2))
  
}


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


#### creating workbook to save results
wb = openxlsx::createWorkbook()
distance_matrix = matrix(0, nrow = length(data_files), ncol = length(data_files))
min_angle_matrix = matrix(0, nrow = length(data_files), ncol = length(data_files))


#### computing pair-wise Sinkhorn distance
chunk = 1
i = chunk
while (i <= (chunk + 9)) {
  sinkhorn_info = data.frame(matrix(0, nrow = 8, ncol = length(data_files)))
  
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
      axon_pp_1 = rescale.ppp(axon_pp_1, s=320153.4)
    }
    axon_pp_1 = shift.ppp(axon_pp_1, origin = "centroid") 
    axon_pp_1 = as.ppp(axon_pp_1)
    
    #construct the j-th point pattern
    file_name_2 = strsplit(data_files[j], ".csv")[[1]]
    print(file_name_2)
    
    axon_locations_2 = unique(read.csv(paste(folder_path, file_name_2, extension, sep="")))
    retrieved_contour_2 = readRDS(paste(folder_path, file_name_2, ".rds", sep=""))
    
    axon_pp_2 = ppp(x=axon_locations_2$X, y=axon_locations_2$Y, checkdup=F, window = retrieved_contour_2)  
    if(scaling == 1){
      print("scaled")
      axon_pp_2 = rescale.ppp(axon_pp_2, s=320153.4)
    }
    axon_pp_2 = shift.ppp(axon_pp_2, origin = "centroid") 
    axon_pp_2 = as.ppp(axon_pp_2)
    
    col_info = c()
    row_info = c()
    
    #loop for rotating the i-th point pattern
    angle_of_rotation = 45 # in degree
    theta = 0
    while(theta < 360){
      cat("theta: ", theta, "\n")
      row_info[length(row_info) + 1] = theta
      
      axon_pp_1_rotated = rotate.ppp(axon_pp_1, angle = theta*pi/180, centre = "centroid")
      svglite(paste(plot_dest_path, paste(file_name_1, theta, "deg rotated", sep = " "), "_pp.svg", sep=""), width = 10, height = 10)
      plot(axon_pp_1_rotated, pch=21, bg="black", cex=0.7, main=paste(file_name_1, theta, "deg rotated", sep = " "))
      dev.off()
      
      result_map = constructMap(pp1=axon_pp_1_rotated, pp2=axon_pp_2, type=analysis_type)
      map_1 = result_map[[1]]
      map_2 = result_map[[2]]
      dens_1 = result_map[[3]]
      dens_2 = result_map[[4]]
      
      jpeg(paste(plot_dest_path, file_name_1,"_", theta, "_", analysis_type, ".jpeg", sep=""),
           width = 10,
           height = 10,
           units = "in",
           res=900)
      image(dens_1, main="", col=colorRampPalette(brewer.pal(9, "YlOrRd"))(9), box=FALSE)
      dev.off()
      
      jpeg(paste(plot_dest_path, file_name_2,"_", analysis_type, ".jpeg", sep=""),
           width = 10,
           height = 10,
           units = "in",
           res=900)
      image(dens_2, main="", col=colorRampPalette(brewer.pal(9, "YlOrRd"))(9), box=FALSE)
      dev.off()
      
      #===============================================
      r_1 = nrow(map_1)
      c_1 = ncol(map_1)
      
      r_2 = nrow(map_2)
      c_2 = ncol(map_2)
      
      #constructing cost/ distance matrix manually
      m1 = seq(0,1,length.out=c_1)
      m2 = seq(0,1,length.out=r_1)
      
      n1 = seq(0,1,length.out=c_2)
      n2 = seq(0,1,length.out=r_2)
      
      expand_m = data.frame(expand.grid(m1,rev(m2)))
      expand_n = data.frame(expand.grid(n1,rev(n2)))
      
      costm = proxy::dist(expand_m, expand_n, method = "Euclidean")
      costm = `dim<-`(c(costm), dim(costm))
      
      #flatten the images
      a = matrix(map_1, r_1*c_1, 1, byrow = T)
      b = matrix(map_2, r_2*c_2, 1, byrow = T)
      
      cat(r_1, "*", c_1, " ", r_2, "*", c_2, "\n")
      cat(dim(a), ", ", dim(b), ", ", dim(costm), "\n")
      
      #using the function from Barycenter package
      tic("Barycenter::Sinkhorn")
      skh_bary = Barycenter::Sinkhorn(a, b, costm, lambda = as.numeric(lambda)) #default regularization param (lambda) is 1.
      toc()
      cat(skh_bary$Distance, "\n")
      #===============================================
      
      col_info[length(col_info) + 1] = skh_bary$Distance
      theta = theta + angle_of_rotation
    }
    
    sinkhorn_info[, j] = col_info
    j = j+1
  }
  
  #include min, max, avg, variance and min index for all the sinkhorn distance per fascicle
  sinkhorn_info[nrow(sinkhorn_info) + 1, ] = apply(sinkhorn_info[1:8, ], 2, min)
  sinkhorn_info[nrow(sinkhorn_info) + 1, ] = apply(sinkhorn_info[1:8, ], 2, max)
  sinkhorn_info[nrow(sinkhorn_info) + 1, ] = apply(sinkhorn_info[1:8, ], 2, mean)
  sinkhorn_info[nrow(sinkhorn_info) + 1, ] = apply(sinkhorn_info[1:8, ], 2, var)
  sinkhorn_info[nrow(sinkhorn_info) + 1, ] = row_info[apply(sinkhorn_info[1:8, ], 2, which.min)]
  
  rownames(sinkhorn_info) = c(row_info, "min", "max", "avg", "var", "argmin")
  colnames(sinkhorn_info) = c(1:length(sinkhorn_info))
  
  addWorksheet(wb, sheetName = paste("Sheet", i, sep=""))
  writeData(wb, sheet = paste("Sheet", i, sep=""), x=sinkhorn_info, rowNames = TRUE)
  
  distance_matrix[i, ] = apply(sinkhorn_info[1:8, ], 2, min)
  min_angle_matrix[i, ] = row_info[apply(sinkhorn_info[1:8, ], 2, which.min)]
  
  cat("\n")
  i=i+1
}

addWorksheet(wb, sheetName = "Distance matrix")
writeData(wb, sheet = "Distance matrix", x=distance_matrix, rowNames = FALSE, colNames = FALSE)

addWorksheet(wb, sheetName = "Min angle matrix")
writeData(wb, sheet = "Min angle matrix", x=min_angle_matrix, rowNames = FALSE, colNames = FALSE)

openxlsx::saveWorkbook(wb, file =paste(plot_dest_path, version, "_", analysis_type, scaling, "_", lambda, "_", chunk,  "_all_result.xlsx", sep=""), overwrite = TRUE)


#### creating workbook to save the min angle matrix separately
wb2 = openxlsx::createWorkbook()

addWorksheet(wb2, sheetName = "Min angle matrix")
writeData(wb2, sheet = "Min angle matrix", x=min_angle_matrix, rowNames = FALSE, colNames = FALSE)

openxlsx::saveWorkbook(wb2, file =paste(parent, "Data/Supporting Files/Angle of Rotation/Angle_of_rotation_Barycenter_", version, "_", scaling, "_", lambda, "_", chunk,  ".xlsx", sep=""), 
                       overwrite = TRUE)


