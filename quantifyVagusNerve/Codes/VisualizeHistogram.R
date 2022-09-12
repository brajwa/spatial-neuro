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


#### input file containing the Sinkhorn distance matrix from previous computation
chosen_file = file.choose()

m = openxlsx::read.xlsx(xlsxFile = chosen_file, sheet = "Distance matrix", colNames = FALSE)
m = as.matrix(m)

vv = m[1:18, 1:18]
vv = vv[upper.tri(vv, diag = TRUE)]
vv_label = rep("Vagus-Vagus", length(vv))

pp = m[19:29, 19:29]
pp = pp[upper.tri(pp, diag = TRUE)]
pp_label = rep("Pelvic-Pelvic", length(pp))

vp = c(m[1:18, 19:29])
vp_label = rep("Vagus-Pelvic", length(vp))

d1 = data.frame(distance=vv, group=vv_label)
d2 = data.frame(distance=pp, group=pp_label)
d3 = data.frame(distance=vp, group=vp_label)

d = rbind(d1, d2, d3)

ggplot(d, aes(x = distance, colour = group)) +
  geom_density()