#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim", "RBGL", "svglite", "ggrepel", "devtools")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)

devtools::install_github("swarm-lab/Rvision")
require("Rvision")

contourNodeCoord <- function(face, pp){
    df1 = data.frame(x = pp$x[as.integer(face)], y = pp$y[as.integer(face)])
    df1 = rbind(df1, df1[1, ])
    df2 = data.frame(df1[order(nrow(df1):1),])
    return(df2)
}


#### from watershed lines
# watershed_line_file = "C:/Users/sanja/Desktop/MAX_File_85_01-31-2019_CYM_1_4_Calb2_3B_YFP_DS_GFP-g_Hu-b.lif - TileScan_001_Merging001_ProjMax001_AdjustClr001-watershed-lines.tif"
# 
# w_lines = Rvision::image(watershed_line_file, colorspace = "GRAY")
# w_lines = Rvision::resize(w_lines, fx = 1/0.6602, fy = 1/0.6602) # to micron scale
# plot(w_lines)
# 
# dim_x = dim(w_lines)[2]
# dim_y = dim(w_lines)[1]
# 
# #w_lines_gray = Rvision::changeColorSpace(w_lines, "GRAY")
# w_lines_bin = w_lines < 200
# plot(w_lines_bin)
# 
# face_contours = Rvision::findContours(w_lines_bin)
# num_of_faces = max(face_contours$contours[, 1]) + 1 #0-indexed


#### from constructed network
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

max_y = 1 #4539.812 # found by computation and converting to nm scale from microns

for (i in c(2)) { # 2,13,21
    ens_location = strsplit(branch_info_files[i], "/")[[1]][11]
    sample_id = strsplit(strsplit(branch_info_files[i], "/")[[1]][12], "\\.")[[1]][1]
    cat("\n(", i, ") Location: ", ens_location, "\nSample Id: ", sample_id, "\n")
    
    data_struct_list = constructDataStruct(sample_id, parent, branch_info_files[i], output_folder_path, max_y)
    
    #### the returned values
    branch.all = data_struct_list[[1]]
    branch.ppp = data_struct_list[[2]]
    branch.lpp = data_struct_list[[3]]
    g1 = data_struct_list[[4]]
    hardcoreStrauss_model_param = data_struct_list[[5]]
    
    plot(branch.lpp, main="", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", "dodgerblue", "white", "maroon1",
                                                           "mediumpurple"))
                                                           #save the network plot
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
    tab_degree = data.frame(table(igraph::degree(g1)))
    tab_degree$Var1 = as.numeric(tab_degree$Var1)
    num_nodes = 0
    sum = 0
    for(j in 1:length(tab_degree$Var1)){
       num_nodes = num_nodes + tab_degree[j, ]$Freq
       sum = sum + (tab_degree[j, ]$Freq * tab_degree[j, ]$Var1)
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
    
    print(ggplot(branch.all, aes(x=angle)) + 
             geom_histogram(aes(y=after_stat(density)), colour="grey", fill="grey", binwidth = 3)+
             geom_density(alpha=1, colour="black", linewidth=1.5) +
             labs(x = "Edge angle", y = "Density", color = "") )
    
    
    #### calculate the edge length and plot the distribution
    #### not scaled
    branch.all$euclid = (apply(branch.all, 1, function(x) calcDist(x)))
    
    print(ggplot(branch.all, aes(x=euclid)) + 
             geom_histogram(aes(y=after_stat(density)), colour="grey", fill="grey", binwidth = 7)+
             geom_density(alpha=1, colour="black", linewidth=1.5) +
             labs(x = "Edge length (Euclidean)", y = "Density", color = "") )
    
    
    #### network faces
    #add additional edges for faces that are cut off at the boundary
    pp_nodes = data.frame(x=branch.ppp$x, y=branch.ppp$y)
    
    boundary_1 = pp_nodes[pp_nodes$x <= branch.ppp$window$xrange[1], ]
    boundary_1 = boundary_1[order(boundary_1$y), ]
    boundary_1 = as.integer(rownames(boundary_1))
    
    boundary_2 = pp_nodes[pp_nodes$y >= branch.ppp$window$yrange[2], ]
    boundary_2 = boundary_2[order(boundary_2$x), ]
    boundary_2 = as.integer(rownames(boundary_2))
    
    boundary_3 = pp_nodes[pp_nodes$x >= branch.ppp$window$xrange[2], ]
    boundary_3 = boundary_3[order(boundary_3$y, decreasing = TRUE), ]
    boundary_3 = as.integer(rownames(boundary_3))
    
    boundary_4 = pp_nodes[pp_nodes$y <= branch.ppp$window$yrange[1], ]
    boundary_4 = boundary_4[order(boundary_4$x, decreasing = TRUE), ]
    boundary_4 = as.integer(rownames(boundary_4))
    
    new_edges = data.frame(n1 = c(boundary_1, boundary_2, boundary_3, boundary_4),
                           n2 = c(boundary_1[2:length(boundary_1)], boundary_2, boundary_3, boundary_4, boundary_1[1]))
    
    g_p =  graph_from_data_frame(rbind(branch.all[, 5:6], new_edges), directed = FALSE)
    
    g <- as_graphnel(g_p) ## Convert igraph object to graphNEL object for planarity testing
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
}


#### face features computation
columns = c("Area", "Ext.", "Disp.", "Elong.", "Eccentr.", "Orient.") 
face_features = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(face_features) = columns

for(f in c(1: length(face_list))){
    cat("face id: ", f, "\n")
    
    #f_contour = face_contours$contours[face_contours$contours[, 1] == f, 2:3]
    f_contour = as.matrix(contourNodeCoord(face_list[[f]], branch.ppp))
    lines(f_contour, col="red", type="l", lwd=2)
    
    #f_contour[, 2] = f_contour[, 2] - dim_y
    #plot(f_contour, col="red", type="l", lwd=2)
    
    area = Rvision::contourArea(f_contour[,1], f_contour[,2])
    
    moments = Rvision::moments(f_contour)
    
    #### rotational invariants
    phi1 = moments$value[moments$moment == "nu02"] + moments$value[moments$moment == "nu20"]
    phi2 = ((moments$value[moments$moment == "nu02"] - moments$value[moments$moment == "nu20"]) * (moments$value[moments$moment == "nu02"] - moments$value[moments$moment == "nu20"])) 
        + (4 * moments$value[moments$moment == "nu11"] * moments$value[moments$moment == "nu11"])
    lambda1 = 2 * pi * (phi1 + sqrt(phi2))
    lambda2 = 2 * pi * (phi1 - sqrt(phi2))
        
    #### ext, disp, elong
    ext = log2(lambda1)
    disp = log2(sqrt(lambda1 * lambda2))
    elong = log2(sqrt(lambda1 / lambda2))
    
    #### orient, eccentr
    orient = 0.5 * atan2((2 * moments$value[moments$moment == "m11"]) , (moments$value[moments$moment == "m20"] - moments$value[moments$moment == "m02"]))
    orient = orient * 180 / pi
    eccentr = (((moments$value[moments$moment == "m02"] - moments$value[moments$moment == "m20"]) * (moments$value[moments$moment == "m02"] - moments$value[moments$moment == "m20"])) 
               + (4 * moments$value[moments$moment == "m11"] * moments$value[moments$moment == "m11"])) / moments$value[moments$moment == "m00"]
    
    face_features = rbind(face_features, data.frame(area, ext, disp, elong, eccentr, orient))
}



