load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)


selectPathToImageJ <- function(){
  path_to_imagej_exe = ""
  if(file.exists("ImageJExePath.txt")){
    content = readLines("ImageJExePath.txt")
    
    if(file.exists(content)){
      path_to_imagej_exe = content
    }
    else{
      print("Choose the path to ImageJ exe in your machine...\n")
      path_to_imagej_exe = file.choose(new = TRUE)
      writeLines(path_to_imagej_exe, "ImageJExePath.txt")
    }
  }
  else{
    print("Choose the path to ImageJ exe in your machine...\n")
    path_to_imagej_exe = file.choose(new = TRUE)
    writeLines(path_to_imagej_exe, "ImageJExePath.txt")
  }
  
  return(path_to_imagej_exe)
}


main <- function(){
  #1
  path_to_imagej_exe = selectPathToImageJ()
  
  "At this point, the automated image processing task should be done to segment an image and extract the branch information in xlsx 
  file, but the process could not be completed as of now because one of the Fiji plugins has GUI. So we segment the image manually and
  provide the path to the information file to the next step directly."
  # path_to_sample_network_image = ""
  # path_to_segment_macro = ""
  # path_to_branch_info = segmentNetworkImage(path_to_imagej_exe, path_to_segment_macro, path_to_sample_network_image)
  
  path_to_branch_info = "D:\\Summer 2019\\R codes\\Research1.0InterganglionicNetwork2021\\Inputs\\Branch Information\\ROHCHOPp01.NewbornLeftColon_Branches_Trimmed.csv"
  
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

