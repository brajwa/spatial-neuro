#### loading required libraries (there might be more libraries loaded than required)
load.lib = c("deldir", "spatstat", "spatstat.utils", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "oce", "OneR", "RandomFieldsUtils", "RandomFields", "Cairo", "knitr", "scorepeak", "Rcpp", "emdist",
             "RImageJROI", "svglite", "transport", "Barycenter", "T4transport", "wvtool", "adimpro", "reshape2", "proxy", "RColorBrewer", "tictoc",
             "ggrepel", "scatterplot3d", "car", "e1071", "rgl", "this.path", "ggbeeswarm")

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

#simple vagus vs. pelvic
all_vagus = d[1:18, 1:18]
all_vagus = all_vagus[upper.tri(all_vagus)]

all_pelvic = d[19:29, 19:29]
all_pelvic = all_pelvic[upper.tri(all_pelvic)]

all_vagus_pelvic = c(d[1:18, 19:29])

all = data.frame(Sinkhorn_distance = c(all_vagus, all_pelvic, all_vagus_pelvic),
                 Category = c(rep("Vag vs. Vag", 153),
                              rep("Pelv vs. Pelv", 55),
                              rep("Vag vs. Pelv", 198)))

write.csv(all, file = paste(chosen_file, "_all_dist.csv", sep=""))

svglite(paste(chosen_file, "_AllPair_Boxplot.svg", sep=""), width = 6, height = 4)

ggplot(all, aes(x = Category, y = Sinkhorn_distance, fill = Category)) +
    geom_boxplot(notch=TRUE, outlier.size = 1) +
    geom_quasirandom(cex=0.5, shape = 21, colour = "grey40", aes(fill=Category)) +
    
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="top", legend.text=element_text(size=16), legend.title = element_blank(),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=18),
          plot.subtitle = element_text(hjust = 0.5, size=16),
          axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", size=0.25, linetype=2)) +
    
    xlab(expression(paste("Sample categories"))) + ylab("Sinkhorn distance")+
    
    labs(title = "Pairwise Sinkhorn Distance (Spatial Intensity)",
         subtitle = "Analysis: Kernel-smoothed Spatial Map")    # the titles needs changing for different runs

dev.off()


#sub-categories
abd_vagus = c(6:11, 14:18)
cerv_vagus = c(1:5, 12, 13)

d_abd_vagus = d[abd_vagus, abd_vagus]
d_abd_vagus = d_abd_vagus[upper.tri(d_abd_vagus)]

d_cerv_vagus = d[cerv_vagus, cerv_vagus]
d_cerv_vagus = d_cerv_vagus[upper.tri(d_cerv_vagus)]

d_abd_cerv_vagus = c(d[abd_vagus, cerv_vagus])

vagus = data.frame(Sinkhorn_distance = c(d_abd_vagus, d_cerv_vagus, d_abd_cerv_vagus), 
                   Category = c(rep("Abd vs. Abd", 55), 
                                rep("Cerv vs. Cerv", 21), 
                                rep("Abd vs. Cerv", 77)))

write.csv(vagus, file = paste(chosen_file, "_intra_vagus_dist.csv", sep=""))

svglite(paste(chosen_file, "_Vagus_Boxplot.svg", sep=""), width = 6, height = 4)

ggplot(vagus, aes(x = Category, y = Sinkhorn_distance, fill = Category)) +
    geom_boxplot(notch=TRUE, outlier.size = 1) +
    geom_quasirandom(cex=0.5, shape = 21, colour = "grey40", aes(fill=Category)) +
    
    scale_fill_brewer(palette="Pastel1") +
    theme(legend.position="top", legend.text=element_text(size=16), legend.title = element_blank(),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=18),
          plot.subtitle = element_text(hjust = 0.5, size=16),
          axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", size=0.25, linetype=2)) +
    
    xlab(expression(paste("Sample categories"))) + ylab("Sinkhorn distance")+
    
    labs(title = "",
         subtitle = "(within vagus samples)")    # the titles needs changing for different runs

dev.off()


pelvic = c(19:29)

d_pelvic_abd = c(d[pelvic, abd_vagus])

d_pelvic_cerv = c(d[pelvic, cerv_vagus])

pelv = data.frame(Sinkhorn_distance = c(d_pelvic_abd, d_pelvic_cerv), 
                   Category = c(rep("Pelvic vs. Abd Vagus", 121), 
                                rep("Pelvic vs. Cerv Vagus", 77)))

write.csv(pelv, file = paste(chosen_file, "_inter_vagus_pelvic_dist.csv", sep=""))

svglite(paste(chosen_file, "_Pelvic_Boxplot.svg", sep=""), width = 6, height = 4)

ggplot(pelv, aes(x = Category, y = Sinkhorn_distance, fill = Category)) +
    geom_boxplot(notch=TRUE, outlier.size = 1) +
    geom_quasirandom(cex=0.5, shape = 21, colour = "grey40", aes(fill=Category)) +
    
    scale_fill_brewer(palette="Pastel2") +
    theme(legend.position="top", legend.text=element_text(size=16), legend.title = element_blank(),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=18),
          plot.subtitle = element_text(hjust = 0.5, size=16),
          axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", size=0.25, linetype=2)) +
    
    xlab(expression(paste("Sample categories"))) + ylab("Sinkhorn distance")+
    
    # labs(title = "",
    #      subtitle = "(between pelvic and vagus samples)")    # the titles needs changing for different runs
    
    labs(title = "",
         subtitle = "(between pelvic and vagus samples)")

dev.off()