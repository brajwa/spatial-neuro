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

####
computeFaceConvexity <- function(face, pp){
    face_node_coords = contourNodeCoord(face, pp)
    face_contour = as.matrix(face_node_coords)
    perim = contourPerimeter(face_contour)
    #cat("perimeter: ", perim, "\n")
    
    c_hull = chull(x=face_node_coords$x, y=face_node_coords$y)
    
    c_hull_coordinates = face_node_coords[c_hull,]
    hull_contour = as.matrix(c_hull_coordinates)
    hull_perim = contourPerimeter(hull_contour)
    #cat("convex hull perimeter: ", hull_perim, "\n")
    
    return(hull_perim / perim)
}


#### shoelace formula for computing the face area
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


#### Given a face of the network as a sequence of vertex id and the point pattern object of the vertices, 
#### this function returns a 2-column dataframe of the vertex coordinates
contourNodeCoord <- function(face, pp){
    df1 = data.frame(x = pp$x[as.integer(face)], y = pp$y[as.integer(face)])
    
    return(df1)
}


#### Given a face contour computes the perimeter
contourPerimeter <- function(f_c){
    f_c = rbind(f_c, f_c[1, ])  # add the first to the last to make a loop
    r = dim(f_c)[1]
    p = 0
    for(e in c(1: (r-1))){
        p = p + calcDist(c(f_c[e, ], f_c[e+1, ]))   # uses Euclidean distance routine
    }
    
    return(as.numeric(p))
}


#### Given the list of all faces (as sequence of vertices) and a particular face id,
#### this function computes all the face features under consideration of the given face id
computeFacefeatures <- function(f, face_list, u_branch.ppp, corner.ppp){
    #cat("face id: ", f, "\n")
    
    #f_contour = face_contours$contours[face_contours$contours[, 1] == 0, 2:3]  # this line was used when contours were computed from watershed lines
    if(!is.null(corner.ppp)){
        f_contour = as.matrix(contourNodeCoord(face_list[[f]], superimpose.ppp(u_branch.ppp, corner.ppp))) # in this case the contour is not a loop, as per example in documentation
    }else{
        f_contour = as.matrix(contourNodeCoord(face_list[[f]], u_branch.ppp))
    }
    
    #lines(f_contour, col="red", type="l", lwd=2) # draws each face on the actual network for ease of verification
    
    area = Rvision::contourArea(f_contour[,1], f_contour[,2])
    
    perim = contourPerimeter(f_contour)
    
    moments = Rvision::moments(f_contour)
    
    #### rotational invariants; use normalized central moments
    phi1 = moments$value[moments$moment == "nu02"] + moments$value[moments$moment == "nu20"]
    phi2 = ((moments$value[moments$moment == "nu02"] - moments$value[moments$moment == "nu20"]) * (moments$value[moments$moment == "nu02"] - moments$value[moments$moment == "nu20"])) 
    + (4 * moments$value[moments$moment == "nu11"] * moments$value[moments$moment == "nu11"])
    lambda1 = 2 * pi * (phi1 + sqrt(phi2))
    lambda2 = 2 * pi * (phi1 - sqrt(phi2))
    
    #### ext, disp, elong
    ext = log2(lambda1)
    disp = log2(sqrt(lambda1 * lambda2))
    elong = log2(sqrt(lambda1 / lambda2))
    
    #### orient, eccentr, use direct spatial moments
    orient = 0.5 * atan2((2 * moments$value[moments$moment == "mu11"]) , (moments$value[moments$moment == "mu20"] - moments$value[moments$moment == "mu02"]))
    orient = orient * 180 / pi
    eccentr = (((moments$value[moments$moment == "mu02"] - moments$value[moments$moment == "mu20"]) * (moments$value[moments$moment == "mu02"] - moments$value[moments$moment == "mu20"])) 
               + (4 * moments$value[moments$moment == "mu11"] * moments$value[moments$moment == "mu11"])) / moments$value[moments$moment == "m00"]
    
    return(data.frame(area, perim, ext, disp, elong, eccentr, orient))
}


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

max_y = 1 # 4539.812 found by computation; right now keeping everything unscaled as the moments can not be computed otherwise

#### a combined dataframe for all the face features of all the samples
columns_combined = c("edge_angle", "edge_len","ens_location", "sample_id") 
edge_features_combined = data.frame(matrix(nrow = 0, ncol = length(columns_combined)))
colnames(edge_features_combined) = columns_combined


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
    
    plot(branch.lpp, main="", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", "dodgerblue", "white", "maroon1",
                                                    "mediumpurple"))
    
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
    
    # ggplot(degree_frame, aes(x=deg)) + 
    #    geom_histogram(aes(y=after_stat(density)), colour="grey", fill="grey", binwidth = 0.5)+ 
    #    geom_density(alpha=1, colour="black", linewidth=1.5) +
    #    scale_x_continuous(limits=c(0, 10), breaks = seq(0,10, by=1))+
    #    labs(x = "Degree of the vertices", y = "Density", color = "")
    
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
    
    # print(ggplot(branch.all, aes(x=angle)) + 
    #          geom_histogram(aes(y=after_stat(density)), colour="grey", fill="grey", binwidth = 3)+
    #          geom_density(alpha=1, colour="black", linewidth=1.5) +
    #          labs(x = "Edge angle", y = "Density", color = "") )
    # 
    
    #### calculate the edge length and plot the distribution
    #### not scaled
    branch.all$euclid = (apply(branch.all, 1, function(x) calcDist(x)))
    
    # print(ggplot(branch.all, aes(x=euclid)) + 
    #          geom_histogram(aes(y=after_stat(density)), colour="grey", fill="grey", binwidth = 7)+
    #          geom_density(alpha=1, colour="black", linewidth=1.5) +
    #          labs(x = "Edge length (Euclidean)", y = "Density", color = "") )
    # 
    
    new_branch_info = data.frame(branch.all$angle, branch.all$euclid)
    if(grepl("dist", ens_location)){
        new_branch_info$ens_location = "distal"
    }else if(grepl("mid", ens_location)){
        new_branch_info$ens_location = "middle"
    }else if(grepl("prox", ens_location)){
        new_branch_info$ens_location = "proximal"
    }
    new_branch_info$sample_id = sample_id
    
    edge_features_combined = rbind(edge_features_combined, new_branch_info)
    
}# loop ends for each sample

sample_id = c("prox-sub-CYM1-4-sam-1-4-3", "mid-sub-CYM1-4-sam-1-4-2", "dist-sub-CYM1-4-sam-1-4-1")
edge_feature = edge_features_combined[edge_features_combined$sample_id %in% sample_id, ]

g = ggplot(edge_feature, aes(x = ens_location, y = (branch.all.angle), fill = ens_location)) +
    #geom_boxplot(notch=TRUE, outlier.size = 1) +
    #geom_quasirandom(cex=0.5, shape = 21, colour = "grey40", aes(fill=ens_location)) +
    
    geom_violinhalf(alpha=0.5, linewidth=0.2, outlier.size = 0.08) +
    coord_flip() +
    geom_boxplot(notch=TRUE, outlier.size = 0.08, width=0.1, linewidth=0.2) +
    
    scale_fill_brewer(palette="Accent") +
    theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(), 
          legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=8),
          axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),
          axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", linewidth=0.2, linetype=2)) +
    
    xlab(expression(paste("ENS location"))) + ylab("Edge angle")
plot(g)

svglite(paste("D:/Fall 2023/Research/Prelim/figures/feat_comp/edge_angle.svg", sep=""), width = 2.5, height = 2)
par(mar=c(0,0,0,0))
plot(g)
dev.off()