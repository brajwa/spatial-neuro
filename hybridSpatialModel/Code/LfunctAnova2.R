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


#### removing the sample with quite different L-function behavior for r-value computation
#### sample: 11328_5_AB-2V_2_montage_pred
data_files = data_files[data_files != "11328_5_AB-2V_2_montage_pred_mi.csv"]


pp_list = list()

### computing the max r values in inhomogeneous L function for each sample
max_inhom_l_r = c()

for (file_name in data_files) {
  file_name = strsplit(file_name, ".csv")[[1]]
  print(file_name)
  
  axon_locations = unique(read.csv(paste(folder_path, file_name, extension, sep="")))
  retrieved_contour = readRDS(paste(folder_path, file_name, ".rds", sep=""))
  
  axon_pp = ppp(x=axon_locations$X, y=axon_locations$Y, checkdup=F, window = retrieved_contour)
  
  #axon_pp = rescale.ppp(axon_pp, s=axon_pp$window$yrange[2])
  axon_pp = rescale.ppp(axon_pp, s=320153.4)  #pre-computed maximum y-range of all the fascicles
  axon_pp = shift.ppp(axon_pp, origin = "centroid")
  
  pp_list[[length(pp_list) + 1]] = axon_pp
  
  inhom_l = Linhom(axon_pp, correction = "border")
  cat("Max r in Linhom: ", max(inhom_l$r), "\n\n")
  
  max_inhom_l_r[length(max_inhom_l_r) + 1] = max(inhom_l$r)
}


#### creating a common range of r-values
r_vect_size = 100
inhom_l_r_vect = seq(0, min(max_inhom_l_r), by=(min(max_inhom_l_r)/r_vect_size))


inhom_l_all = data.frame(sample=data_files, group=c(rep(1, 17), rep(2, 11))) #1=vagus,2=pelvic


#### computing inhomogeneous L function for the common range of r values (just computed above) for all the samples
#### to figure out the most varying r. This r will be used in local computations.
for(i in c(1:length(pp_list))){
  axon_pp = pp_list[[i]]
  inhom_l = Linhom(axon_pp, r=inhom_l_r_vect, correction = "border")
  
  inhom_l_all[i, 3:103] = inhom_l$border
}


f_val_vector = c()

for(i in c(3:103)){
  anova_res = summary(aov(inhom_l_all[, i]~inhom_l_all$group))
  f_val_vector[i-2] = anova_res[[1]]$`F value`[1]
}

f_ratio = data.frame(r=inhom_l_r_vect, f_ratio=f_val_vector)
f_ratio = f_ratio[order(f_ratio$f_ratio, decreasing = TRUE),]

#### top 3 r values with top variances
top_r = f_ratio[1:3,]$r

top_r[1] # the r-value to be used