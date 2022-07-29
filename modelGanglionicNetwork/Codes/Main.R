#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)


main <- function(){
    #### Step:1
    
    #### extracting parent directory information for accessing input and output location
    dir = this.dir()
    folder = strsplit(dir, "/")
    folder = folder[[1]][length(folder[[1]])]
    parent = strsplit(dir, folder)
    
    #### the TIF images of the ganglionic networks are prepocessed in Fiji (ImageJ) and the network information is extracted as .csv files
    #### choose one of the ganglionic network samples with file chooser below
    setwd(paste(parent, "Data/Branch Information (in um)/", sep=""))
    path_to_branch_info = file.choose()
        
    #2
    ganglia_info_list = analyzeGanglia(path_to_branch_info)
    
    Beta = ganglia_info_list[[1]]
    Gamma = ganglia_info_list[[2]]
    R = ganglia_info_list[[3]]
    H = ganglia_info_list[[4]]
    window = ganglia_info_list[[5]]
    ganglia_fitted_model = ganglia_info_list[[6]]
    
    print(ganglia_fitted_model)
    
    #3
    set.seed(Sys.time())
    
    g = generateGangliaCenters(Beta, Gamma, R, H, window=window, process_type=3, with_model=1, fitted_model=ganglia_fitted_model)
    plotGeneratedGanglia(g)
    
    ltest= envelope(g, Lest)
    ggplot(data = ltest, aes(r)) + 
    
    geom_ribbon(aes(ymin=lo-r, ymax=hi-r), alpha = 0.3)+
    geom_line(aes(y=lo-r, colour="CSR lower bound"), size=2) +
    geom_line(aes(y=hi-r, colour="CSR upper bound"), size=2) +
    #geom_point(aes(y=theo, colour="csr")) + 
    
    geom_line(aes(y=obs-r, colour="simulated"), size=2) + 
    #geom_point(aes(y=obs, colour="observed")) +  
    
    theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
          legend.text=element_text(size=30),
          axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
    xlab("Distance [pixels]") + ylab("L Function")
    
    #4
    branch_info_list = analyzeBranch(path_to_branch_info)
    
    branch_all = branch_info_list[[1]]
    orgKDE_angle = branch_info_list[[2]]
    orgKDE_length = branch_info_list[[3]]
    orgKDE_both = branch_info_list[[4]]
    meshedness = branch_info_list[[5]]
    network_density = branch_info_list[[6]]
    compactness = branch_info_list[[7]]
    
    network_info_list = generateNetworkEdges(g, branch_all, orgKDE_angle, orgKDE_length, orgKDE_both,
                                           meshedness, network_density, compactness)
    
    n = network_info_list[[1]] 
    g2_lin = network_info_list[[2]]
    EMD = network_info_list[[3]]
    
    plot(g2_lin)
    print(EMD)
    
    #5
    "The intensity profile for the neurons is generated manually in Fiji and provided to the next step."
    
    #6
    
    #7
    "Done in the generateNeuronCenters function"
    
    #8
}


main()

