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
    
    #### extract which sample we will be working on from the file name (the data files are to be saved in a particular format)
    sample_id = strsplit(path_to_branch_info, "\\\\")[[1]]
    sample_id = sample_id[length(sample_id)]
    sample_id = (strsplit(sample_id, "_")[[1]])[1]
    
    #### creating an output folder for the working sample
    output_folder_path = paste(parent, "Outputs/", sample_id, "/", sep="")
    if (!dir.exists(output_folder_path)) {dir.create(output_folder_path, recursive=TRUE)}
        
    #### Step:2
    
    #### analyzing the extracted information for generating a spatial model (hardcore-Strauss process) for the ganglia
    ganglia_info_list = analyzeGanglia(sample_id, parent, path_to_branch_info, output_folder_path)
    
    #### returned values
    Beta = ganglia_info_list[[1]]   #stationary intensity parameter (not required in our case)
    Gamma = ganglia_info_list[[2]]  #interaction parameter
    R = ganglia_info_list[[3]]      #interaction distance
    H = ganglia_info_list[[4]]      #hardcore distance
    window = ganglia_info_list[[5]]
    ganglia_fitted_model = ganglia_info_list[[6]]
    
    print(ganglia_fitted_model)
    
    #### Step:3
    
    #### generate ganglia with the fitted spatial model
    set.seed(Sys.time())
    
    g = generateGangliaCenters(Beta, Gamma, R, H, window=window, process_type=3, with_model=1, fitted_model=ganglia_model)
    plotGeneratedGanglia(g)
    
    ganglia_ppp = g
    
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
    
    stats = createWorkbook() 
    doc = read_pptx()
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(code = plot(ganglia_ppp, main=file_name, cex=1, pch=20, bg=1)), location = ph_location_fullsize())
    
    ltest= envelope(ganglia_ppp, Lest)
    ggobj = ggplot(data = ltest, aes(r)) + 
        
        geom_ribbon(aes(ymin=lo-r, ymax=hi-r), alpha = 0.3)+
        geom_line(aes(y=lo-r, colour="csr lower bound"), size=2) +
        geom_line(aes(y=hi-r, colour="csr upper bound"), size=2) +
        #geom_point(aes(y=theo, colour="csr")) + 
        
        geom_line(aes(y=obs-r, colour="simulated"), size=2) + 
        #geom_point(aes(y=obs, colour="observed")) +  
        
        theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
              legend.text=element_text(size=30),
              axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
              axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
        xlab("Distance") + ylab("L Function")+
        ggtitle(file_name)
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
    
    print(doc, target = paste(output_folder_path, "Simulated_PPP.pptx", sep=""))
    
    # write.csv(data.frame(x=ganglia_ppp$x, y=ganglia_ppp$y), 
    #         "D://Summer 2019//R codes//Research1.0InterganglionicNetwork2021//Outputs//Simulated thingys//ganglia_coord.csv", 
    #         row.names = F)
    
    
    
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

