#### loading required libraries (there might be more libraries loaded than required)
load.lib = c("deldir", "spatstat", "spatstat.utils", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "oce", "OneR", "RandomFieldsUtils", "RandomFields", "Cairo", "knitr", "scorepeak", "Rcpp", "emdist",
             "RImageJROI", "svglite", "transport", "Barycenter", "T4transport", "wvtool", "adimpro", "reshape2", "proxy", "RColorBrewer", "tictoc",
             "ggrepel", "scatterplot3d", "car", "e1071", "rgl", "this.path")

install.lib = load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib, dependencies=TRUE)
sapply(load.lib, require, character=TRUE)


#### extracting parent directory information for accessing input and output location
dir = this.dir()
folder = strsplit(dir, "/")
folder = folder[[1]][length(folder[[1]])]
parent = strsplit(dir, folder)


#### input location
folder_path = paste(parent, "Data/Inputs/", sep = "")
data_files = list.files(path=folder_path, pattern = ".csv", full.names = FALSE)
extension = ".csv"


pp_list = list()
y_ranges = c()

for (file_name in data_files) {
  file_name = strsplit(file_name, ".csv")[[1]]
  print(file_name)
  
  axon_locations = unique(read.csv(paste(folder_path, file_name, extension, sep="")))
  retrieved_contour = readRDS(paste(folder_path, file_name, ".rds", sep=""))
  
  axon_pp = ppp(x=axon_locations$X, y=axon_locations$Y, checkdup=F, window = retrieved_contour)
  
  pp_list[[length(pp_list) + 1]] = axon_pp
  
  cat("y-range: ", axon_pp$window$yrange[2], "\n\n")
  
  y_ranges[length(y_ranges) + 1] = axon_pp$window$yrange[2]
}

max(y_ranges)