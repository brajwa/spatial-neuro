#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)


setwd("~/GitHub/spatial-neuro/modelGanglionicNetwork/Codes")

#### source the functions from other files
source("AnalyzeGanglionicNetwork.R")
source("GenerateGangliaCenters.R")
source("GenerateNetwork.R")


#### extracting parent directory information for accessing input and output location
dir = this.dir()
folder = strsplit(dir, "/")
folder = folder[[1]][length(folder[[1]])]
parent = strsplit(dir, folder)

moment_folder = paste(parent, "Data/ENSMouse Moment Information (in um)/", sep="")
moment_files = list.files(moment_folder, recursive = TRUE, pattern = "\\.csv", full.names = TRUE)

face_feature_list = c("Area", "Perim.", "Ext.", "Disp.", "Elong.", "Eccentr.", "Orient.") 
select_feature = 5

for (i in c(1:11)) {
    moment_info = read.csv(moment_files[i])
    
    feature_info = moment_info[[face_feature_list[select_feature]]]
    #check if NaN present
    
    print(ggplot(data.frame(feature_info=feature_info), aes(x = "Distal", y = feature_info, fill = "Distal")) +
        geom_boxplot(notch=TRUE, outlier.size = 1) +
        #geom_quasirandom(cex=0.5, shape = 21, colour = "grey40", aes(fill=ENSLocation)) +
        
        scale_fill_brewer(palette="Set3") +
        theme(legend.position="top", legend.text=element_text(size=16), legend.title = element_blank(),
              legend.box.margin=margin(-10,-10,-10,-10),
              plot.title = element_text(hjust = 0.5, size=18),
              plot.subtitle = element_text(hjust = 0.5, size=16),
              axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
              axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
              panel.background = element_rect(fill='white', colour='black'),
              panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
        
        xlab(expression(paste("ENS location"))) + ylab(paste("Network Face ", face_feature_list[select_feature]))+
        
        labs(title = "") )   # the titles needs changing for different runs
    
    
}
