#### loading required libraries (there might be more libraries loaded than required)
load.lib = c("deldir", "spatstat", "spatstat.utils", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "oce", "OneR", "RandomFieldsUtils", "RandomFields", "Cairo", "knitr", "scorepeak", "Rcpp", "emdist",
             "RImageJROI", "svglite", "transport", "Barycenter", "T4transport", "wvtool", "adimpro", "reshape2", "proxy", "RColorBrewer", "tictoc",
             "ggrepel", "scatterplot3d", "car", "e1071", "rgl", "this.path", "fpc", "GGally")

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
dd= cmdscale(d, k=4)

#### creating data frame with all data
df=data.frame(x=dd[,1], y=dd[,2], z=dd[,3], u=dd[, 4])

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

lda = fpc::discrcoord(xd = df[, 1:4], clvecd = as.numeric(factor(df[, 7])))

lda_df = data.frame(lda$proj)
names(lda_df) = c("LD1","LD2","LD3", "LD4")
plot(lda_df, col=as.numeric(factor(df[, 7])))

#### visualize
svglite(paste(chosen_file, "_LDA_4.svg", sep=""), width = 7, height = 7)

ggpairs(lda_df, aes(colour=df[, 7]))+
  theme(legend.position="top",  legend.title=element_blank(), legend.text=element_text(size=10), 
        legend.box.margin=margin(-10,-10,-10,-10),
        plot.title = element_text(hjust = 0.5, size=10),
        plot.subtitle = element_text(hjust = 0.5, size=9),
        axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(color = "grey", size=0.25, linetype=2))
  # labs(title = "",    # the titles needs changing for different runs
  #      subtitle = )

dev.off()


