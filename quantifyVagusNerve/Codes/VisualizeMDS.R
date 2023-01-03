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
svglite(paste(chosen_file, "_MDS_Poster.svg", sep=""), width = 9, height = 6)

ggplot(data = df, aes(x = x, y = y, label=Image_ID)) +
  #geom_point(size = 1.5, shape=21, colour="grey40", aes(fill = Nerve)) +
    
  geom_point(size = 3, aes(shape = Sex, colour=Nerve, fill= Nerve), show.legend = TRUE) +
  geom_point(size = 3, colour="grey40", aes(shape = Sex, fill= Nerve), show.legend = FALSE) +
  scale_shape_manual(values = c(21, 24)) + 
    
  geom_text_repel(size=5, 
                  segment.size=0.2, segment.color="grey40",
                  arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "first"),
                  max.overlaps = Inf) +
  
  #stat_ellipse(type="t", aes(colour=Nerve.Location.Avr)) +
  
  theme_bw()+
  theme(legend.position="top", legend.text=element_text(size=20),
        legend.box.margin=margin(-10,-10,-10,-10),
        plot.title = element_text(hjust = 0.5, size=24),
        legend.title = element_text(size=20),
        plot.subtitle = element_text(hjust = 0.5, size=20),
        axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(color = "grey", size=0.25, linetype=2)) +
  
  xlab(expression(paste("Dimension 1"))) + ylab("Dimension 2") +
  
  # labs(title = expression("Sinkhorn Space of Spatial Intensity, "*lambda*" = 0.01"),    # the titles needs changing for different runs
  #   subtitle = "Analysis: Kernel-smoothed Spatial Map")
    
  labs(title = expression("Sinkhorn Space of Spatial Satistics, "*lambda*" = 0.01"),    # the titles needs changing for different runs
    subtitle = "Analysis: Kernel-smoothed Spatial Map\nLocal Inhomogeneous L-Function")
  
  # labs(subtitle = "no sector")
  # labs(subtitle = expression("vertical sector (90 ± 15) "*degree))
    
  # labs(title = "Local Inhomogeneous L-Function",
  #      subtitle = expression("vertical sector (90 ± 15) "*degree))

dev.off()


