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
r = nrow(m)
c = ncol(m)

#### creating the symmetric matrix
d1 = m
d2 = t(m)
d = matrix(rep(0, r*c), nrow = r)

d[upper.tri(d)] = d1[upper.tri(d1)]
d[lower.tri(d)] = d2[lower.tri(d2)]

isSymmetric(d)
diag(d)

#### multi-dimensional scaling to map entities into the Sinkhorn space of spatial features
dd= cmdscale(d)

#### creating data frame with all data
df=data.frame(x=dd[,1], y=dd[,2])

#### the sample id, nerve, nerve location and sex info of the fascicles are saved in the following xlsx file,
#### in sequence of the the files processed
if(grepl("demo", chosen_file, fixed = TRUE)){
  sample_info = openxlsx::read.xlsx(xlsxFile = paste(parent, "Data/Supporting Files/Biological Information/Bio_info_demo.xlsx", sep=""), 
                                    rowNames = FALSE)
}else{
  sample_info = openxlsx::read.xlsx(xlsxFile = paste(parent, "Data/Supporting Files/Biological Information/Bio_info.xlsx", sep=""), 
                                    rowNames = FALSE)
}


#### including the biological information
df = cbind(df, sample_info)

#### visualize
ggplot(data = df, aes(x = x, y = y, shape=Nerve, colour = Nerve, label=Sample_ID_Avr)) +
  geom_point(size = 4) +
  geom_text_repel(size=3.5, max.overlaps = Inf, show.legend = F) +
  
  theme_bw()+
  theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=28), legend.title=element_blank(),
        legend.text=element_text(size=28),
        axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24),
        axis.title.x = element_text(size = 28), axis.title.y = element_text(size = 28),
        panel.grid.major = element_line(color="gray"),
        panel.grid.minor = element_line(color="gray")) +
  xlab(expression(paste("MDS dimension 1"))) + ylab("MDS dimension 2")+
  ggtitle("Visualizing the Sinkhorn space")




