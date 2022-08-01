#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)


####################################
generateNeuronCenters <- function(){
    #### extracting parent directory information for accessing input and output location
    dir = this.dir()
    folder = strsplit(dir, "/")
    folder = folder[[1]][length(folder[[1]])]
    parent = strsplit(dir, folder)
    
    #### the images of the neuron intensity profiles of the simulated ganglionic network are stored
    #### choose one of the intensity profiles with file chooser below
    setwd(paste(parent, "Data/Neuron Intensity Profiles/", sep=""))
    path_to_intensity_image = file.choose()
    
    #### extract which sample we will be working on from the file name (the data files are to be saved in a particular format)
    image_name = strsplit(path_to_intensity_image, "\\\\")[[1]]
    image_name = image_name[length(image_name)]
    image_name = (strsplit(image_name, "_")[[1]])[1]
    
    #### creating an output folder for the working sample
    output_folder_path = paste(parent, "Outputs/", image_name, "/", sep="")
    if (!dir.exists(output_folder_path)) {dir.create(output_folder_path, recursive=TRUE)}
    
    #### the spatial model parameters to generate neuron centers
    b_val = 0.003
    h_val = 0.01
    g_val =0.78
    r_val = 0.0285
    
    intensity_profile = cimg2im(load.image(path_to_intensity_image))
    
    #### simulate neuron centers with hardcore-Strauss process
    simulated_neuron_1 = rStraussHard(beta = b_val, gamma = g_val, R = r_val, H = h_val, 
                                      W=owin(intensity_profile$xrange, intensity_profile$yrange))
    
    plot(simulated_neuron_1, bg='magenta', pch=21, cex=0.9)
    
    #### thin the simulated neuron centers based on the loaded intensity profile
    simulated_neuron = rthin(X=simulated_neuron_1, P=intensity_profile)
    plot(simulated_neuron, bg='magenta', pch=21, cex=0.9)
    
    #### adjust the window of the simulated neurons
    simulated_neuron$y = simulated_neuron$y - simulated_neuron$window$yrange[2]
    simulated_neuron$x = simulated_neuron$x - simulated_neuron$window$xrange[1]
    simulated_neuron$window$xrange = simulated_neuron$window$xrange- simulated_neuron$window$xrange[1]
    simulated_neuron$window$yrange = simulated_neuron$window$yrange- simulated_neuron$window$yrange[2]
    
    simulated_neuron$y = simulated_neuron$y / abs(simulated_neuron$window$yrange[1]) * abs(ganglia_ppp$window$yrange[1])
    simulated_neuron$x = simulated_neuron$x /  abs(simulated_neuron$window$xrange[2]) * abs(ganglia_ppp$window$xrange[2])
    simulated_neuron$window$xrange = c(0, abs(ganglia_ppp$window$xrange[2]))
    simulated_neuron$window$yrange = c(ganglia_ppp$window$yrange[1], 0)
    
    #### creating a directory to save the simulated neuron related files
    sim_neuron_path = paste(output_folder_path, "Simulated Neuron/", sep="")
    if (!dir.exists(sim_neuron_path)) {dir.create(sim_neuron_path, recursive=TRUE)}
    
    write.csv(data.frame(x=simulated_neuron$x, y=simulated_neuron$y),
              paste(sim_neuron_path, image_name, "_Simulated_Neuron_coord.csv", sep=""),
              row.names = F)
}