#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim", "RBGL", "svglite")

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

#### the TIF images of the ganglionic networks are preprocessed in Fiji (ImageJ) and the network information is extracted as .csv files
#### choose one/all of the ganglionic network samples with file chooser below
branch_info_folder = paste(parent, "Data/ENSMouse Branch Information (in um) v2.0/", sep="")
branch_info_files = list.files(branch_info_folder, recursive = TRUE, pattern = "\\.csv", full.names = TRUE)

network_feature_list = c("Node degree", "Edge angle", "Edge length", "Face area", "Face node count") 
select_feature = 4
cat("Network feature under consideration: ", network_feature_list[select_feature], "\n")

columns = c("feature_value","ens_location","sample_id") 
feature_info_combined = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(feature_info_combined) = columns

for (i in c(1: length(branch_info_files))) { #2,13,21
  ens_location = strsplit(branch_info_files[i], "/")[[1]][11]
  sample_id = strsplit(strsplit(branch_info_files[i], "/")[[1]][12], "\\.")[[1]][1]
  cat("\n(", i, ") Location: ", ens_location, "\nSample Id: ", sample_id, "\n")

  #### creating an output folder for the working sample
  output_folder_path = paste(parent, "Outputs/ENSMouse/", sample_id, "/", sep="")
  if (!dir.exists(output_folder_path)) {dir.create(output_folder_path, recursive=TRUE)}
  
  
  data_struct_list = constructDataStruct(sample_id, parent, branch_info_files[i], output_folder_path)
  
  #### the returned values
  branch.all = data_struct_list[[1]]
  branch.ppp = data_struct_list[[2]]
  branch.lpp = data_struct_list[[3]]
  g1 = data_struct_list[[4]]
  hardcoreStrauss_model_param = data_struct_list[[5]]
  
  plot(branch.lpp, main="", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", "dodgerblue", "white", "maroon1",
                                                         "mediumpurple"))
  #save the network plot
  pdf(paste(output_folder_path, "Networkv2.0.pdf", sep=""), width = 6, height = 6)
  par(mar = c(0,0,0,0))
  plot(branch.lpp, main="", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", "dodgerblue", "white", "maroon1",
                                                         "mediumpurple"))
  dev.off()
                                                        
  g1$layout = as.matrix(data.frame(x=branch.ppp$x, y=branch.ppp$y))
  plot(g1,
       vertex.size=2*igraph::degree(g1),
       vertex.label=NA,
       vertex.shape="circle",
       vertex.color=igraph::degree(g1))
  
  
  #### spatial point pattern metrics
  num_points_pp = branch.ppp$n
  area_pp = summary(branch.ppp)$window$area
  ints_pp = intensity.ppp(branch.ppp)
  
  closest_point_dist = min(nndist(branch.ppp))
  farthest_point_dist = max(pairdist(branch.ppp))
  
  cat("Number of points: ", num_points_pp, "\n",
      "Window area: ", area_pp, " um squared\n",
      "Homogeneous intensity of points: ", ints_pp, "\n",
      "Closest point distance: ", closest_point_dist, "\n",
      "Farthest point distance: ", farthest_point_dist, "\n")
  
  
  #### spatial network metrics
  
  # degree distribution
  table_degree = table(igraph::degree(g1))
  
  degree_frame = as.data.frame(igraph::degree(g1))
  colnames(degree_frame)[colnames(degree_frame) == 'igraph::degree(g1)'] = 'deg'
  
  ggplot(degree_frame, aes(x=deg)) + 
    geom_histogram(aes(y=after_stat(density)), colour="grey", fill="grey", binwidth = 0.5)+ 
    geom_density(alpha=1, colour="black", linewidth=1.5) +
    scale_x_continuous(limits=c(0, 10), breaks = seq(0,10, by=1))+
    labs(x = "Degree of the vertices", y = "Density", color = "")
  
  cat("Degree distribution:", table_degree, "\n")
  
  #### calculating the edge probability p from the degree distribution
  #### to use while generating ER random graph G(n, p)
  num_nodes = 0
  sum = 0
  for(i in 1:length(table_degree)){
    num_nodes = num_nodes + table_degree[[i]]
    sum = sum + (i * table_degree[[i]])
  }
  edge_probability = divide_by(sum, num_nodes^2)
  cat("Edge probability: ", edge_probability, "\n")
  
  #### Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
  #### This is sometimes also called the clustering coefficient.
  cluster_coeff = igraph::transitivity(g1, type = "global")
  cat("CC: ", cluster_coeff, "\n")
  
  #### alpha, gamma, psi
  N = num_points_pp
  E = length(branch.all$x1)
  A = area_pp
  L = sum(branch.all$euclid)
  
  meshedness = (E-N+1)/((2*N)-5)
  network_density = E/((3*N)-6)
  compactness = 1- ((4*A)/(L-(2*sqrt(A)))^2)
  
  cat("Meshedness: ", meshedness, "\n",
      "Network density: ", network_density, "\n",
      "Compactness: ", compactness, "\n")
  
  
  #### calculate the edge angle (in degree) and plot the distribution
  #### not scaled
  branch.all$angle = (apply(branch.all, 1, function(x) calcAngle(x)))
  
  ggplot(branch.all, aes(x=angle)) + 
    geom_histogram(aes(y=after_stat(density)), colour="grey", fill="grey", binwidth = 3)+
    geom_density(alpha=1, colour="black", linewidth=1.5) +
    labs(x = "Edge angle", y = "Density", color = "")
  
  
  #### calculate the edge length and plot the distribution
  #### not scaled
  branch.all$euclid = (apply(branch.all, 1, function(x) calcDist(x)))
  
  ggplot(branch.all, aes(x=euclid)) + 
    geom_histogram(aes(y=after_stat(density)), colour="grey", fill="grey", binwidth = 7)+
    geom_density(alpha=1, colour="black", linewidth=1.5) +
    labs(x = "Edge length (Euclidean)", y = "Density", color = "")
  
  
  #### netwrok faces
  g <- as_graphnel(g1) ## Convert igraph object to graphNEL object for planarity testing
  boyerMyrvoldPlanarityTest(g)
  
  face_list = planarFaceTraversal(g)
  face_node_count = sapply(face_list, length)
  
  #shoelace formula
  faceArea <- function(face, branch.ppp){
    n = length(face)
    area = 0
    
    for(i in c(1:(n-1))){
      area = area + ( branch.ppp$x[as.integer(face[i])] * branch.ppp$y[as.integer(face[i+1])] )
    }
    
    area = area + ( branch.ppp$x[as.integer(face[n])] * branch.ppp$y[as.integer(face[1])] )
    
    for(i in c(1:(n-1))){
      area = area - ( branch.ppp$x[as.integer(face[i+1])] * branch.ppp$y[as.integer(face[i])] )
    }
    
    area = area - ( branch.ppp$x[as.integer(face[1])] * branch.ppp$y[as.integer(face[n])] )
    
    area = abs(area) / 2
    
    return(area)
  }
  
  face_area_list = sapply(face_list, function(x) faceArea(x, branch.ppp))
  
  face_node_count = face_node_count[-which.max(face_area_list)]
  face_list = face_list[-which.max(face_area_list)]
  face_area_list = face_area_list[-which.max(face_area_list)]
  
  hist(face_node_count, breaks=100)
  #plot(density(face_area_list))
  
  ggplot(data.frame(area=face_area_list)) + 
    geom_density(aes(x=area, y=after_stat(density)), alpha=1, colour="black", linewidth=1.5) +
    labs(x = "Area of face", y = "Density", color = "")
  
  #combine the selected feature 
  if(select_feature == 1){
    
  }else if(select_feature == 2){
    edge_count = E
    feature_info_combined = rbind(feature_info_combined, data.frame(feature_value=branch.all$angle, 
                                                                    ens_location=rep(ens_location, edge_count), 
                                                                    sample_id=rep(sample_id, edge_count)) )
    
  }else if(select_feature == 3){
    edge_count = E
    feature_info_combined = rbind(feature_info_combined, data.frame(feature_value=branch.all$euclid, 
                                                                    ens_location=rep(ens_location, edge_count), 
                                                                    sample_id=rep(sample_id, edge_count)) )
    
  }else if(select_feature==4){
    face_count = length(face_area_list)
    feature_info_combined = rbind(feature_info_combined, data.frame(feature_value=face_area_list, 
                                                                    ens_location=rep(ens_location, face_count), 
                                                                    sample_id=rep(sample_id, face_count)) )
  }else if(select_feature == 5){
    face_count = length(face_area_list)
    feature_info_combined = rbind(feature_info_combined, data.frame(feature_value=face_node_count, 
                                                                    ens_location=rep(ens_location, face_count), 
                                                                    sample_id=rep(sample_id, face_count)) )
  }
}

# boxplot of select feature
ggplot(feature_info_combined, aes(x = ens_location, y = (feature_value), fill = ens_location)) +
  geom_boxplot(notch=TRUE, outlier.size = 1) +
  #geom_quasirandom(cex=0.5, shape = 21, colour = "grey40", aes(fill=ens_location)) +
  
  scale_fill_brewer(palette="Set3") +
  theme(legend.position="top", legend.text=element_text(size=16), legend.title = element_blank(),
        legend.box.margin=margin(-10,-10,-10,-10),
        plot.title = element_text(hjust = 0.5, size=18),
        plot.subtitle = element_text(hjust = 0.5, size=16),
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
  
  xlab(expression(paste("ENS location"))) + ylab(paste("(", network_feature_list[select_feature], ")",sep = "") )+
  
  labs(title = paste("Statistical comparison of Network ", network_feature_list[select_feature] ," from different part of ENS", sep="") )   # the titles needs changing for different runs

