#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim", "RBGL", "svglite", "ggrepel", "devtools", "geosphere")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)

devtools::install_github("swarm-lab/Rvision")
require("Rvision")

#### from constructed network
setwd("~/GitHub/spatial-neuro/modelGanglionicNetwork/Codes")

#### source the functions from other files
source("AnalyzeGanglionicNetwork.R")


inhomLPlot <- function(sample_id, u_branch.ppp, min_max_r){
    #### creating a common range of r-values
    r_vect_size = 100
    inhom_l_r_vect = seq(0, min_max_r, by=min_max_r/r_vect_size)
    
    lohboot_linhom = lohboot(u_branch.ppp, fun = "Linhom", correction = "isotropic")
    env_linhom = envelope(u_branch.ppp, fun = "Linhom", correction = "isotropic", nsim = 39)
    
    my_color_map = c("CSR"="darkgrey", "Observed"="cyan3")
    g = ggplot()+
        geom_ribbon(aes(x=env_linhom$r, ymin=env_linhom$lo-env_linhom$r, ymax=env_linhom$hi-env_linhom$r, fill="CSR"), alpha=0.3) +
        
        geom_ribbon(aes(x=lohboot_linhom$r, ymin=lohboot_linhom$lo-lohboot_linhom$r, ymax=lohboot_linhom$hi-lohboot_linhom$r, fill="Observed"), alpha=0.3, show.legend = FALSE) +
        geom_line(aes(x=lohboot_linhom$r, y=lohboot_linhom$iso-lohboot_linhom$r, colour="Observed"), show.legend = FALSE) +
        
        theme(legend.position="top",  legend.title=element_blank(), legend.text=element_text(size=9),
              legend.box.margin=margin(0,-10,-10,-10),
              plot.title = element_text(hjust = 0.5, size=10),
              plot.subtitle = element_text(hjust = 0.5, size=8),
              axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9),
              axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
              panel.background = element_rect(fill='white', colour='black'),
              panel.grid.major = element_line(color = "grey", size=0.25, linetype=2)) +
        scale_color_manual(values = my_color_map) + 
        scale_fill_manual(values = my_color_map) +
        xlab(expression(paste("Interaction distance (r)"))) + ylab("L(r)") +
        labs(title="Besag's centered inhomogeneous L-function", subtitle = sample_id)
    
    print(g)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/proposal spp/linhom_", sample_id, ".svg", sep=""), width = 3.5, height = 3)
    par(mar = c(0, 0, 0, 0))  
    print(g)
    dev.off()
    
}


computeCommonMaxR <- function(branch_info_files, parent, output_folder_path, max_y){
    cat("Computing common max R...\n")
    max_inhom_l_r = c()
    
    for (i in c(1:length(branch_info_files))) {
        ens_location = strsplit(branch_info_files[i], "/")[[1]][11]
        sample_id = strsplit(strsplit(branch_info_files[i], "/")[[1]][11], "\\.")[[1]][1]

        data_struct_list = constructDataStruct(sample_id, parent, branch_info_files[i], output_folder_path, max_y)
        
        #### the returned values
        branch.all = data_struct_list[[1]]
        branch.ppp = data_struct_list[[2]]
        branch.lpp = data_struct_list[[3]]
        g1 = data_struct_list[[4]]
        hardcoreStrauss_model_param = data_struct_list[[5]]
        
        u_branch.ppp = unmark(branch.ppp)
        l_inhom = Linhom(u_branch.ppp, correction = "border")
        max_inhom_l_r[length(max_inhom_l_r) + 1] = max(l_inhom$r)
    }
    
    return(min(max_inhom_l_r))
}


#### extracting parent directory information for accessing input and output location
dir = this.dir()
folder = strsplit(dir, "/")
folder = folder[[1]][length(folder[[1]])]
parent = strsplit(dir, folder)

#### the TIF images of the ganglionic networks are preprocessed in Fiji (ImageJ) and the network information is extracted as .csv files
#### choose one/all of the ganglionic network samples with file chooser below
branch_info_folder = paste(parent, "Data/ENSMouse Branch Information (in um) v2.0/", sep="")
branch_info_files = list.files(branch_info_folder, recursive = TRUE, pattern = "\\.csv", full.names = TRUE)

max_y = 1 # 4539.812 found by computation; right now keeping everything unscaled as the moments can not be computed otherwise

min_max_r = computeCommonMaxR(branch_info_files, parent, output_folder_path, max_y)

pp_info = data.frame()
for (i in c(1:length(branch_info_files))) { # 2,13,21
    ens_location = strsplit(branch_info_files[i], "/")[[1]][11]
    sample_id = strsplit(strsplit(branch_info_files[i], "/")[[1]][11], "\\.")[[1]][1]
    cat("\n(", i, ") Location: ", ens_location, "\nSample Id: ", sample_id, "\n")
    
    data_struct_list = constructDataStruct(sample_id, parent, branch_info_files[i], output_folder_path, max_y)
    
    #### the returned values
    branch.all = data_struct_list[[1]]
    branch.ppp = data_struct_list[[2]]
    branch.lpp = data_struct_list[[3]]
    g1 = data_struct_list[[4]]
    hardcoreStrauss_model_param = data_struct_list[[5]]
    
    u_branch.ppp = unmark(branch.ppp)
    
    plot(u_branch.ppp, main="", pch=19, cex=1.2)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/proposal spp/pp_", sample_id, ".svg", sep=""), width = 2, height = 1.5)
    par(mar = c(0, 0, 0, 0))    
    plot(u_branch.ppp, main="", pch=19, cex=0.4)
    dev.off()
    
    inhomLPlot(sample_id, u_branch.ppp, min_max_r)
    
    pp_nn_dist = nndist(u_branch.ppp)
    pp_metrics = c(u_branch.ppp$n, summary(u_branch.ppp)$window$area, intensity.ppp(u_branch.ppp),
                   min(pp_nn_dist), max(pp_nn_dist), mean(pp_nn_dist), sd(pp_nn_dist),
                   max(pairdist(u_branch.ppp)), mean(branch.all$euclid), sd(branch.all$euclid))
    pp_info = rbind(pp_info, c(sample_id, pp_metrics))
    
}# loop ends for each sample

col_head = c("sample_id", "point_count", "window_area", "point_intensity", "min_nn_dist", "max_nn_dist", "avg_nn_dist",
             "sd_nn_dist", "max_pair_dist", "avg_branch_len", "sd_branch_len")
colnames(pp_info) = col_head

