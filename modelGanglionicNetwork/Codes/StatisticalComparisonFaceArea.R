#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim", "RBGL", "ggbeeswarm")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)

#### face area
all_face_area_distal = unlist(combined_f_a_l[c(1:11)])
all_face_area_middle = unlist(combined_f_a_l[c(12:20)])
all_face_area_proximal = unlist(combined_f_a_l[c(21:35)])


ens_location = c(rep("distal", length(all_face_area_distal)),
                 rep("middle", length(all_face_area_middle)),
                 rep("proximal", length(all_face_area_proximal)))

all_face_area = data.frame(FaceArea=c(all_face_area_distal, all_face_area_middle, all_face_area_proximal),
                           ENSLocation= ens_location)

ggplot(all_face_area, aes(x = ENSLocation, y = FaceArea, fill = ENSLocation)) +
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
    
    xlab(expression(paste("ENS location"))) + ylab("Face area of the network")+
    
    labs(title = "Statistical comparison of face area from different part of ENS")    # the titles needs changing for different runs
