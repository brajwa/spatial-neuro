#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim", "RBGL", "svglite", "ggrepel", "devtools", "geosphere", "philentropy",
             "collections")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)

devtools::install_github("swarm-lab/Rvision")
require("Rvision")

setwd("~/GitHub/spatial-neuro/modelGanglionicNetwork/Codes")
#### source the functions from other files
source("AnalyzeGanglionicNetwork.R")


#### given a vector of values, in our case vertex degree (consider them as class labels)
#### outputs the softmax values for each of them as a dictionary structure
mySoftMax <- function(y){
    num = exp(y)
    denom = sum(num)
    soft_max = num / denom
    
    deg_prob = dict()
    
    for (i in c(1:length(y))) {
        deg_prob$set(as.character(y[i]), soft_max[i])
    }
    
    return(deg_prob)
}


#### computes the face convexity (perimeter based) from its convex hull
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


computeFaceList <- function(gen.ppp, network_extra2){
    temp_graph_obj = make_empty_graph() %>% add_vertices(gen.ppp$n)
    temp_graph_obj = add_edges(as.undirected(temp_graph_obj),
                               as.vector(t(as.matrix(network_extra2[,5:6]))))
    
    temp_g_o <- as_graphnel(temp_graph_obj)
    boyerMyrvoldPlanarityTest(temp_g_o)
    
    temp_face_list = planarFaceTraversal(temp_g_o)
    temp_face_node_count = sapply(temp_face_list, length)
    
    #### applying the shoe lace formula
    temp_face_area_list = sapply(temp_face_list, function(x) faceArea(x, gen.ppp))
    
    #### eliminating the outer face, it has the largest face area
    temp_face_node_count = temp_face_node_count[-which.max(temp_face_area_list)]
    temp_face_list = temp_face_list[-which.max(temp_face_area_list)]
    temp_face_area_list = temp_face_area_list[-which.max(temp_face_area_list)]
    
    return(temp_face_list)
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
    #f_contour = face_contours$contours[face_contours$contours[, 1] == 0, 2:3]  # this line was used when contours were computed from watershed lines
    
    if(!is.null(corner.ppp)){
        #cat("computing face features with corner ppp\n")
        f_contour = as.matrix(contourNodeCoord(face_list[[f]], superimpose.ppp(u_branch.ppp, corner.ppp))) # in this case the contour is not a loop, as per example in documentation
    }else{
        #cat("computing face features without corner ppp\n")
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
    
    #cat(area, perim, ext, disp, elong, eccentr, orient, "\n")
    return(data.frame(area, perim, ext, disp, elong, eccentr, orient))
}


#### plot the face features of two given networks for comparison 
#### can be original vs. initial DT or original vs. final network
comparePlotOrgSim <- function(org_face_feature, face_features, branch.all, network_extra1){
    cat("Number of faces in the original network: ", length(org_face_feature$X))
    cat("\nNumber of faces in the triangulated network: ", length(face_features$area), "\n")
    
    myColors = c("blue", "red")
    names(myColors) = c("real", "triangulated")
    custom_colors = scale_colour_manual(name = "Network type", values = myColors)
    custom_colors_2 = scale_fill_manual(name = "Network type", values = myColors)
    
    # node count
    den_org_nc = density(org_face_feature$Node_Count)
    den_org_nc = data.frame(x=den_org_nc$x, y=den_org_nc$y)
    
    den_sim_nc = density(face_features$Node_Count)
    den_sim_nc = data.frame(x=den_sim_nc$x, y=den_sim_nc$y)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/init tri figures/init_tri_node_count.svg", sep=""), width = 1.85, height = 1.35)
    print(ggplot() +
              geom_line(data = den_org_nc, aes(x=x, y=y, color = "real"), linewidth=0.1) +
              geom_area(data = den_org_nc, aes(x=x, y=y, fill = "real"), alpha=0.3) +
              geom_line(data = den_sim_nc, aes(x=x, y=y, color="triangulated"), linewidth=0.1) +
              geom_area(data = den_sim_nc, aes(x=x, y=y, fill="triangulated"), alpha=0.3) +
              
              theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(), 
                    legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
                    legend.box.margin=margin(-10,-10,-10,-10),
                    plot.title = element_text(hjust = 0.5, size=10),
                    plot.subtitle = element_text(hjust = 0.5, size=8),
                    axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),
                    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                    panel.background = element_rect(fill='white', colour='black'),
                    panel.grid.major = element_line(color = "grey", linewidth=0.2, linetype=2)) +
              custom_colors + custom_colors_2 +
              xlab(expression(paste("Face Node Count"))) + ylab("Density")  )# the titles needs changing for different runs
    dev.off()
    
    
    # area
    den_org_area = density(org_face_feature$Area_SL)
    den_org_area = data.frame(x=den_org_area$x, y=den_org_area$y)
    
    den_sim_area = density(face_features$area)
    den_sim_area = data.frame(x=den_sim_area$x, y=den_sim_area$y)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/init tri figures/init_tri_face_area.svg", sep=""), width = 1.85, height = 1.35)
    print(ggplot() +
              geom_line(data = den_org_area, aes(x=x, y=y, color = "real"), linewidth=0.1) +
              geom_area(data = den_org_area, aes(x=x, y=y, fill = "real"), alpha=0.3) +
              geom_line(data = den_sim_area, aes(x=x, y=y, color="triangulated"), linewidth=0.1) +
              geom_area(data = den_sim_area, aes(x=x, y=y, fill="triangulated"), alpha=0.3) +
              
              theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(), 
                    legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
                    legend.box.margin=margin(-10,-10,-10,-10),
                    plot.title = element_text(hjust = 0.5, size=10),
                    plot.subtitle = element_text(hjust = 0.5, size=8),
                    axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),
                    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                    panel.background = element_rect(fill='white', colour='black'),
                    panel.grid.major = element_line(color = "grey", linewidth=0.2, linetype=2)) +
              scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
              custom_colors + custom_colors_2 +
              xlab(expression(paste("Face Area"))) + ylab("Density"))  # the titles needs changing for different runs
    dev.off()
    
    
    # elongation
    den_org_elong = density(org_face_feature$Elong.)
    den_org_elong = data.frame(x=den_org_elong$x, y=den_org_elong$y)
    
    den_sim_elong = density(face_features$elong)
    den_sim_elong = data.frame(x=den_sim_elong$x, y= den_sim_elong$y)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/init tri figures/init_tri_face_elong.svg", sep=""), width = 1.85, height = 1.35)
    print(ggplot() +
              geom_line(data = den_org_elong, aes(x=x, y=y, color = "real"), linewidth=0.1) +
              geom_area(data = den_org_elong, aes(x=x, y=y, fill = "real"), alpha=0.3) +
              geom_line(data = den_sim_elong, aes(x=x, y=y, color="triangulated"), linewidth=0.1) +
              geom_area(data = den_sim_elong, aes(x=x, y=y, fill="triangulated"), alpha=0.3) +
              
              theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(), 
                    legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
                    legend.box.margin=margin(-10,-10,-10,-10),
                    plot.title = element_text(hjust = 0.5, size=10),
                    plot.subtitle = element_text(hjust = 0.5, size=8),
                    axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),
                    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                    panel.background = element_rect(fill='white', colour='black'),
                    panel.grid.major = element_line(color = "grey", linewidth=0.2, linetype=2)) +
              custom_colors + custom_colors_2 +
              xlab(expression(paste("Face Elongation"))) + ylab("Density")  )  # the titles needs changing for different runs
    dev.off()
    
    
    # orientation
    den_org_orient = density(org_face_feature$Orient.)
    den_org_orient = data.frame(x=den_org_orient$x, y=den_org_orient$y)
    
    den_sim_orient = density(face_features$orient)
    den_sim_orient = data.frame(x=den_sim_orient$x, y= den_sim_orient$y)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/init tri figures/init_tri_face_orient.svg", sep=""), width = 1.85, height = 1.35)
    print(ggplot() +
              geom_line(data = den_org_orient, aes(x=x, y=y, color = "real"), linewidth=0.1) +
              geom_area(data = den_org_orient, aes(x=x, y=y, fill = "real"), alpha=0.3) +
              geom_line(data = den_sim_orient, aes(x=x, y=y, color="triangulated"), linewidth=0.1) +
              geom_area(data = den_sim_orient, aes(x=x, y=y, fill="triangulated"), alpha=0.3) +
              
              theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(), 
                    legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
                    legend.box.margin=margin(-10,-10,-10,-10),
                    plot.title = element_text(hjust = 0.5, size=10),
                    plot.subtitle = element_text(hjust = 0.5, size=8),
                    axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),
                    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                    panel.background = element_rect(fill='white', colour='black'),
                    panel.grid.major = element_line(color = "grey", linewidth=0.2, linetype=2)) +
              custom_colors + custom_colors_2 +
              xlab(expression(paste("Face Orientation"))) + ylab("Density")  )  # the titles needs changing for different runs
    dev.off()
    
    
    # edge length
    den_org_e_len = density(branch.all$euclid)
    den_org_e_len = data.frame(x=den_org_e_len$x, y=den_org_e_len$y)
    
    den_sim_e_len = density(network_extra1$euclidDist)
    den_sim_e_len = data.frame(x=den_sim_e_len$x, y= den_sim_e_len$y)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/init tri figures/init_tri_edge_len.svg", sep=""), width = 1.85, height = 1.35)
    print(ggplot() +
              geom_line(data = den_org_e_len, aes(x=x, y=y, color = "real"), linewidth=0.1) +
              geom_area(data = den_org_e_len, aes(x=x, y=y, fill = "real"), alpha=0.3) +
              geom_line(data = den_sim_e_len, aes(x=x, y=y, color="triangulated"), linewidth=0.1) +
              geom_area(data = den_sim_e_len, aes(x=x, y=y, fill="triangulated"), alpha=0.3) +
              
              theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(), 
                    legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
                    legend.box.margin=margin(-10,-10,-10,-10),
                    plot.title = element_text(hjust = 0.5, size=10),
                    plot.subtitle = element_text(hjust = 0.5, size=8),
                    axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),
                    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                    panel.background = element_rect(fill='white', colour='black'),
                    panel.grid.major = element_line(color = "grey", linewidth=0.2, linetype=2)) +
              custom_colors + custom_colors_2 +
              xlab(expression(paste("Edge Length"))) + ylab("Density")  )  # the titles needs changing for different runs
    dev.off()
    
    
    # edge angle
    den_org_e_angle = density(apply(branch.all, 1, function(x) calcAngle(x)))
    den_org_e_angle = data.frame(x=den_org_e_angle$x, y=den_org_e_angle$y)
    
    den_sim_e_angle = density(network_extra1$anglecomp)
    den_sim_e_angle = data.frame(x=den_sim_e_angle$x, y= den_sim_e_angle$y)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/init tri figures/init_tri_edge_angle.svg", sep=""), width = 1.85, height = 1.35)
    print(ggplot() +
              geom_line(data = den_org_e_angle, aes(x=x, y=y, color = "real"), linewidth=0.1) +
              geom_area(data = den_org_e_angle, aes(x=x, y=y, fill = "real"), alpha=0.3)+
              geom_line(data = den_sim_e_angle, aes(x=x, y=y, color="triangulated"), linewidth=0.1) +
              geom_area(data = den_sim_e_angle, aes(x=x, y=y, fill="triangulated"), alpha=0.3) +
              
              theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(), 
                    legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
                    legend.box.margin=margin(-10,-10,-10,-10),
                    plot.title = element_text(hjust = 0.5, size=10),
                    plot.subtitle = element_text(hjust = 0.5, size=8),
                    axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),
                    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                    panel.background = element_rect(fill='white', colour='black'),
                    panel.grid.major = element_line(color = "grey", linewidth=0.2, linetype=2)) +
              custom_colors + custom_colors_2 +
              xlab(expression(paste("Edge Angle"))) + ylab("Density") ) # the titles needs changing for different runs
    dev.off()
    
}

comparePlotOrgSim2 <- function(org_face_feature, face_features, branch.all, network_extra1){
    cat("Number of faces in the original network: ", length(org_face_feature$X))
    cat("\nNumber of faces in the simulated network: ", length(face_features$area), "\n")
    
    myColors = c("blue", "green")
    names(myColors) = c("real", "simulated")
    custom_colors = scale_colour_manual(name = "Network type", values = myColors)
    custom_colors_2 = scale_fill_manual(name = "Network type", values = myColors)
    
    # node count
    den_org_nc = density(org_face_feature$Node_Count)
    den_org_nc = data.frame(x=den_org_nc$x, y=den_org_nc$y)
    
    den_sim_nc = density(face_features$Node_Count)
    den_sim_nc = data.frame(x=den_sim_nc$x, y=den_sim_nc$y)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/sim figures/sim_node_count.svg", sep=""), width = 1.85, height = 1.35)
    print(ggplot() +
              geom_line(data = den_org_nc, aes(x=x, y=y, color = "real"), linewidth=0.1) +
              geom_area(data = den_org_nc, aes(x=x, y=y, fill = "real"), alpha=0.3) +
              geom_line(data = den_sim_nc, aes(x=x, y=y, color="simulated"), linewidth=0.1) +
              geom_area(data = den_sim_nc, aes(x=x, y=y, fill="simulated"), alpha=0.3) +
              
              theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(), 
                    legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
                    legend.box.margin=margin(-10,-10,-10,-10),
                    plot.title = element_text(hjust = 0.5, size=10),
                    plot.subtitle = element_text(hjust = 0.5, size=8),
                    axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),
                    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                    panel.background = element_rect(fill='white', colour='black'),
                    panel.grid.major = element_line(color = "grey", linewidth=0.2, linetype=2)) +
              custom_colors + custom_colors_2 +
              xlab(expression(paste("Face Node Count"))) + ylab("Density"))# the titles needs changing for different runs
    dev.off()
    
    
    # area
    den_org_area = density(org_face_feature$Area_SL)
    den_org_area = data.frame(x=den_org_area$x, y=den_org_area$y)
    
    den_sim_area = density(face_features$area)
    den_sim_area = data.frame(x=den_sim_area$x, y=den_sim_area$y)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/sim figures/sim_face_area.svg", sep=""), width = 1.85, height = 1.35)
    print(ggplot() +
              geom_line(data = den_org_area, aes(x=x, y=y, color = "real"), linewidth=0.1) +
              geom_area(data = den_org_area, aes(x=x, y=y, fill = "real"), alpha=0.3) +
              geom_line(data = den_sim_area, aes(x=x, y=y, color="simulated"), linewidth=0.1) +
              geom_area(data = den_sim_area, aes(x=x, y=y, fill="simulated"), alpha=0.3) +
              
              theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(), 
                    legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
                    legend.box.margin=margin(-10,-10,-10,-10),
                    plot.title = element_text(hjust = 0.5, size=10),
                    plot.subtitle = element_text(hjust = 0.5, size=8),
                    axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),
                    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                    panel.background = element_rect(fill='white', colour='black'),
                    panel.grid.major = element_line(color = "grey", linewidth=0.2, linetype=2)) +
              scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
              custom_colors + custom_colors_2 +
              xlab(expression(paste("Face Area"))) + ylab("Density"))  # the titles needs changing for different runs
    dev.off()
    
    
    # elongation
    den_org_elong = density(org_face_feature$Elong.)
    den_org_elong = data.frame(x=den_org_elong$x, y=den_org_elong$y)
    
    den_sim_elong = density(face_features$elong)
    den_sim_elong = data.frame(x=den_sim_elong$x, y= den_sim_elong$y)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/sim figures/sim_face_elong.svg", sep=""), width = 1.85, height = 1.35)
    print(ggplot() +
              geom_line(data = den_org_elong, aes(x=x, y=y, color = "real"), linewidth=0.1) +
              geom_area(data = den_org_elong, aes(x=x, y=y, fill = "real"), alpha=0.3) +
              geom_line(data = den_sim_elong, aes(x=x, y=y, color="simulated"), linewidth=0.1) +
              geom_area(data = den_sim_elong, aes(x=x, y=y, fill="simulated"), alpha=0.3) +
              
              theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(), 
                    legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
                    legend.box.margin=margin(-10,-10,-10,-10),
                    plot.title = element_text(hjust = 0.5, size=10),
                    plot.subtitle = element_text(hjust = 0.5, size=8),
                    axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),
                    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                    panel.background = element_rect(fill='white', colour='black'),
                    panel.grid.major = element_line(color = "grey", linewidth=0.2, linetype=2)) +
              custom_colors + custom_colors_2 +
              xlab(expression(paste("Face Elongation"))) + ylab("Density"))  # the titles needs changing for different runs
    dev.off()
    
    
    # orientation
    den_org_orient = density(org_face_feature$Orient.)
    den_org_orient = data.frame(x=den_org_orient$x, y=den_org_orient$y)
    
    den_sim_orient = density(face_features$orient)
    den_sim_orient = data.frame(x=den_sim_orient$x, y= den_sim_orient$y)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/sim figures/sim_face_orient.svg", sep=""), width = 1.85, height = 1.35)
    print(ggplot() +
              geom_line(data = den_org_orient, aes(x=x, y=y, color = "real"), linewidth=0.1) +
              geom_area(data = den_org_orient, aes(x=x, y=y, fill = "real"), alpha=0.3) +
              geom_line(data = den_sim_orient, aes(x=x, y=y, color="simulated"), linewidth=0.1) +
              geom_area(data = den_sim_orient, aes(x=x, y=y, fill="simulated"), alpha=0.3) +
              
              theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(), 
                    legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
                    legend.box.margin=margin(-10,-10,-10,-10),
                    plot.title = element_text(hjust = 0.5, size=10),
                    plot.subtitle = element_text(hjust = 0.5, size=8),
                    axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),
                    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                    panel.background = element_rect(fill='white', colour='black'),
                    panel.grid.major = element_line(color = "grey", linewidth=0.2, linetype=2)) +
              custom_colors + custom_colors_2 +
              xlab(expression(paste("Face Orientation"))) + ylab("Density"))  # the titles needs changing for different runs
    dev.off()
    
    
    # edge length
    den_org_e_len = density(branch.all$euclid)
    den_org_e_len = data.frame(x=den_org_e_len$x, y=den_org_e_len$y)
    
    den_sim_e_len = density(network_extra1$euclidDist)
    den_sim_e_len = data.frame(x=den_sim_e_len$x, y= den_sim_e_len$y)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/sim figures/sim_edge_len.svg", sep=""), width = 1.85, height = 1.35)
    print(ggplot() +
              geom_line(data = den_org_e_len, aes(x=x, y=y, color = "real"), linewidth=0.1) +
              geom_area(data = den_org_e_len, aes(x=x, y=y, fill = "real"), alpha=0.3) +
              geom_line(data = den_sim_e_len, aes(x=x, y=y, color="simulated"), linewidth=0.1) +
              geom_area(data = den_sim_e_len, aes(x=x, y=y, fill="simulated"), alpha=0.3) +
              
              theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(), 
                    legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
                    legend.box.margin=margin(-10,-10,-10,-10),
                    plot.title = element_text(hjust = 0.5, size=10),
                    plot.subtitle = element_text(hjust = 0.5, size=8),
                    axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),
                    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                    panel.background = element_rect(fill='white', colour='black'),
                    panel.grid.major = element_line(color = "grey", linewidth=0.2, linetype=2)) +
              custom_colors + custom_colors_2 +
              xlab(expression(paste("Edge Length"))) + ylab("Density"))  # the titles needs changing for different runs
    dev.off()
    
    
    # edge angle
    den_org_e_angle = density(apply(branch.all, 1, function(x) calcAngle(x)))
    den_org_e_angle = data.frame(x=den_org_e_angle$x, y=den_org_e_angle$y)
    
    den_sim_e_angle = density(network_extra1$anglecomp)
    den_sim_e_angle = data.frame(x=den_sim_e_angle$x, y= den_sim_e_angle$y)
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/sim figures/sim_edge_angle.svg", sep=""), width = 1.85, height = 1.35)
    print(ggplot() +
              geom_line(data = den_org_e_angle, aes(x=x, y=y, color = "real"), linewidth=0.1) +
              geom_area(data = den_org_e_angle, aes(x=x, y=y, fill = "real"), alpha=0.3)+
              geom_line(data = den_sim_e_angle, aes(x=x, y=y, color="simulated"), linewidth=0.1) +
              geom_area(data = den_sim_e_angle, aes(x=x, y=y, fill="simulated"), alpha=0.3) +
              
              theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(), 
                    legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'),
                    legend.box.margin=margin(-10,-10,-10,-10),
                    plot.title = element_text(hjust = 0.5, size=10),
                    plot.subtitle = element_text(hjust = 0.5, size=8),
                    axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),
                    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
                    panel.background = element_rect(fill='white', colour='black'),
                    panel.grid.major = element_line(color = "grey", linewidth=0.2, linetype=2)) +
              custom_colors + custom_colors_2 +
              xlab(expression(paste("Edge Angle"))) + ylab("Density")) # the titles needs changing for different runs
    dev.off()
    
}


########################################################################
#### computes an edge weight based on the degree of its two end vertices
#### g2_degree is the list of degree of each vertex in the graph under consideration
#### x is the edge whose weight is being computed, the 5th and 6th index provides the indices of the end vertices
computeEdgeWeight <- function(g2_degree, x){
    return(g2_degree[x[5]] + g2_degree[x[6]])
}


####Given the degree of vertices of the current network and the max degree of the original network,
#### compute the probability of each vertex
computeVertexProb <- function(org_max_deg, g2_degree, network_extra1){
    y = c(min(g2_degree):max(g2_degree))
    deg_prob = mySoftMax(y)
    
    v_prob = c()
    for (d in g2_degree) {
        v_prob = c(v_prob, deg_prob$get(as.character(d)))
    }
    return(v_prob)
}

computeVertexProb2 <- function(org_max_deg, g2_degree, network_extra1){
    g2_degree_n = c()
    for (i in c(1:length(g2_degree))) {
        adj_vertices_from_df = c(network_extra1$ind1[network_extra1$ind2==i],
                                 network_extra1$ind2[network_extra1$ind1==i])
        g2_degree_n = c(g2_degree_n, g2_degree[i] + sum(g2_degree[adj_vertices_from_df]) )
    }
    
    y = c(min(g2_degree_n):max(g2_degree_n))
    deg_prob = mySoftMax(y)
    
    v_prob = c()
    for (d in g2_degree_n) {
        v_prob = c(v_prob, deg_prob$get(as.character(d)))
    }
    return(v_prob)
}


####given a vertex id and the network structure, finds its neighbors in a clockwise sequence
antiClockwiseNeighbors <- function(hd_v, gen.ppp, network_extra1){
    adj_vertices_from_df = c(network_extra1$ind1[network_extra1$ind2==hd_v],
                             network_extra1$ind2[network_extra1$ind1==hd_v])
    
    angles = c()
    for(v in adj_vertices_from_df){
        nom = gen.ppp$y[v]-gen.ppp$y[hd_v]
        denom = gen.ppp$x[v]-gen.ppp$x[hd_v]
        angle = (as.numeric(atan2(nom, denom)) *180/pi) + 180
        
        angles = c(angles, angle)
    }
    
    anti_clockwise_neighbors = adj_vertices_from_df[sort(angles, index.return=TRUE)$ix]
    return(anti_clockwise_neighbors)
}


#### Given a face and an edge find if that edge is present in the face
isEdgeOnFace <- function(face, edge){
    face = as.numeric(face)
    matched_index = which(face %in% edge)
    if(length(matched_index) < 2){
        return(FALSE)
    }else{
        if((abs(matched_index[1] - matched_index[2]) == 1) | (abs(matched_index[1] - matched_index[2]) == (length(face)-1))){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }
}


####Given a point pattern object, this function computes the boundary edges 
####to get closed boundary faces
####the additional edges this function computes only depend on the point pattern
####as it only connects the corner points of the point pattern window and 
####the boundary points of the point pattern sequentially
computeBoundaryEdges <- function(branch.ppp){
    pp_nodes = data.frame(x=branch.ppp$x, y=branch.ppp$y)
    
    ####include the window corner nodes to the point pattern
    pp_nodes = rbind(pp_nodes, data.frame(x=c(branch.ppp$window$xrange[1], branch.ppp$window$xrange[2], branch.ppp$window$xrange[2], branch.ppp$window$xrange[1]),
                                          y=c(branch.ppp$window$yrange[2], branch.ppp$window$yrange[2], branch.ppp$window$yrange[1], branch.ppp$window$yrange[1])))
    
    #### nodes at the left side of the network, ordered from bottom to top
    boundary_1 = pp_nodes[pp_nodes$x <= branch.ppp$window$xrange[1], ]
    boundary_1 = boundary_1[order(boundary_1$y), ]
    boundary_1 = as.integer(rownames(boundary_1))   # vertex ids
    
    #### nodes at the top, ordered from left to right
    boundary_2 = pp_nodes[pp_nodes$y >= branch.ppp$window$yrange[2], ]
    boundary_2 = boundary_2[order(boundary_2$x), ]
    boundary_2 = as.integer(rownames(boundary_2))
    
    #### nodes at the right side, ordered from top to bottom
    boundary_3 = pp_nodes[pp_nodes$x >= branch.ppp$window$xrange[2], ]
    boundary_3 = boundary_3[order(boundary_3$y, decreasing = TRUE), ]
    boundary_3 = as.integer(rownames(boundary_3))
    
    #### nodes at the bottom, ordered from right to left
    boundary_4 = pp_nodes[pp_nodes$y <= branch.ppp$window$yrange[1], ]
    boundary_4 = boundary_4[order(boundary_4$x, decreasing = TRUE), ]
    boundary_4 = as.integer(rownames(boundary_4))
    
    #### the new edges will connect the boundary nodes computed above sequentially and 
    ####create a loop by connecting the last node to the first one.
    #### to avoid shifting error
    if(length(boundary_1) == 1){
        b_1 = c()
    }else{
        b_1 = boundary_1[2:length(boundary_1)]
    }
    
    new_edges = data.frame(ind1 = c(boundary_1, boundary_2, boundary_3, boundary_4),
                           ind2 = c(b_1, boundary_2, boundary_3, boundary_4, boundary_1[1]))
    
    new_edges = new_edges[new_edges$ind1 != new_edges$ind2, ]  # removing self loops at the corner nodes
    
    return(new_edges)
    
}


#### given a face (as sequence of vertex ids) and the network edges
#### this function returns the edge ids that are on the face
compuateFaceEdges <- function(face, network_extra, gen.ppp){
    n = length(face)
    edges = c()
    for (i in c(1:(n-1))) {
        for (j in c((i+1):n)) {
            v1 = as.numeric(face[i])
            v2 = as.numeric(face[j])
            lines(c(gen.ppp$x[v1], gen.ppp$x[v2]), c(gen.ppp$y[v1], gen.ppp$y[v2]), col="blue", lwd=3)
            edges = c(edges, which((network_extra$ind1 == v1 & network_extra$ind2 == v2) |
                      (network_extra$ind1 == v2 & network_extra$ind2 == v1)) )
        }
    }
    return(edges)
}


#### Given the pont pattern and the data structure of the network, computes the clustering coefficient by constructing
#### the corresponding graph object
ccFromDataframe <- function(branch.ppp, network_extra){
    graph_obj =  make_empty_graph() %>% add_vertices(branch.ppp$n)
    graph_obj = add_edges(as.undirected(graph_obj), 
                          as.vector(t(as.matrix(network_extra[,5:6]))))
    
    #### Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
    #### This is sometimes also called the clustering coefficient.
    cluster_coeff = igraph::transitivity(graph_obj, type = "global")
    
    return(cluster_coeff)
}

####Given the pont pattern and the data structure of the network, computes the other network metrics
netMetrics <- function(gen.ppp, network_extra){
    #### alpha, gamma, psi (meshedness, network density and compactness parameters)
    N = gen.ppp$n
    E = length(network_extra$x1)
    A = summary(gen.ppp)$window$area
    L = sum(network_extra$euclidDist)
    
    meshedness = (E-N+1)/((2*N)-5)
    network_density = E/((3*N)-6)
    compactness = 1- ((4*A)/(L-(2*sqrt(A)))^2)
    
    cat("Meshedness: ", meshedness, ", Network density: ", network_density, ", Compactness: ", compactness, "\n")
    return(list(meshedness, network_density, compactness))
}


#### given a single or list of edge index and the current network structure,
#### eliminated the given edges from the network and 
#### recomputes all network features
eliminateEdges <- function(gen.ppp, network_extra1, edges_to_eliminate){
    if(is.null(edges_to_eliminate)){
        temp_network_extra1 = network_extra1
    }else{
        temp_network_extra1 = network_extra1[-c(edges_to_eliminate), ]
    }
    
    temp_graph_obj = make_empty_graph() %>% add_vertices(gen.ppp$n)
    temp_graph_obj = add_edges(as.undirected(temp_graph_obj), 
                               as.vector(t(as.matrix(temp_network_extra1[,5:6]))))
    
    temp_g_o <- as_graphnel(temp_graph_obj) 
    boyerMyrvoldPlanarityTest(temp_g_o)
    
    temp_face_list = planarFaceTraversal(temp_g_o)
    temp_face_node_count = sapply(temp_face_list, length)
    
    #### applying the shoe lace formula
    temp_face_area_list = sapply(temp_face_list, function(x) faceArea(x, gen.ppp))
    
    #### eliminating the outer face, it has the largest face area
    temp_face_node_count = temp_face_node_count[-which.max(temp_face_area_list)]
    temp_face_list = temp_face_list[-which.max(temp_face_area_list)]
    temp_face_area_list = temp_face_area_list[-which.max(temp_face_area_list)]
    
    #### face features computation
    temp_columns = c("Area_CF", "Perim.", "Ext.", "Disp.", "Elong.", "Eccentr.", "Orient.") # Area_CF: from contour function
    temp_face_features = data.frame(matrix(nrow = 0, ncol = length(temp_columns)))
    colnames(temp_face_features) = temp_columns
    
    temp_face_convexity_list = c()
    
    for(f in c(1: length(temp_face_list))){
        temp_f_feat = computeFacefeatures(f, temp_face_list, gen.ppp, NULL)
        temp_face_features = rbind(temp_face_features, temp_f_feat)
        
        temp_face_convexity_list = c(temp_face_convexity_list, computeFaceConvexity(temp_face_list[[f]], gen.ppp))
        
    }# loop ends for each face of the network
    temp_face_features$Node_Count = temp_face_node_count
    
    temp_triKDE_face_feat_1 = kde(as.matrix(data.frame(temp_face_area_list, temp_face_features$orient)))
    temp_triKDE_face_feat_2 = kde(as.matrix(data.frame((temp_face_node_count))), 
                                  h=density((temp_face_node_count))$bw)
    temp_triKDE_edge_feat = kde(as.matrix(data.frame((network_extra1$anglecomp),
                                                     (network_extra1$euclidDist))))
    
    #### remove the edges, there is a change
    cat("Edge(s) deleted\n")
    noChange = 0
    
    #### make necessary changes permanent
    network_extra1 = temp_network_extra1
    g2_degree = igraph::degree(temp_graph_obj, mode="total")
    print(table(g2_degree))
    
    face_list = temp_face_list
    face_area_list = temp_face_area_list
    face_node_count = temp_face_node_count
    triKDE_face_feat_1 = temp_triKDE_face_feat_1
    triKDE_face_feat_2 = temp_triKDE_face_feat_2
    triKDE_edge_feat = temp_triKDE_edge_feat
    tri_face_features = temp_face_features
    
    tri_face_convexity_mean = mean(temp_face_convexity_list)
    
    if(!is_connected(temp_graph_obj)){
        cat("NETWORK GOT DISCONNECTED\n")
    }
    
    return(list(noChange, network_extra1, g2_degree, face_list, face_area_list, face_node_count,
                triKDE_face_feat_1, triKDE_face_feat_2, triKDE_edge_feat, tri_face_features, tri_face_convexity_mean))
}


####rewire edges
rewireEdges <- function(gen.ppp, network_extra1, high_deg_vertices, g2_degree, face_list, org_face_convexity_mean,
                        orgKDE_edge_feat, triKDE_edge_feat){
    #### distance of each vertex from the point pattern boundary
    vertex_dist_boundary = bdist.points(gen.ppp)
    
    for(hd_v in high_deg_vertices){
        
        cat("\nhd_v: ", hd_v, "\n")
        
        #### detect the anti clockwise neighbors of a given vertex, wrapped the first one to the last
        anti_clockwise_neighbors = antiClockwiseNeighbors(hd_v, gen.ppp, network_extra1)
        anti_clockwise_neighbors[length(anti_clockwise_neighbors) + 1] = anti_clockwise_neighbors[1]
        
        for (n in c(1: (length(anti_clockwise_neighbors)-1))) {
            #### check for cc?
           
            p = anti_clockwise_neighbors[n]
            q = anti_clockwise_neighbors[n+1]
            
            pq_index = which((network_extra1$ind1==p & network_extra1$ind2==q) |
                           (network_extra1$ind2==p & network_extra1$ind1==q) )
            
            if(length(pq_index)==1){
                cat("Rewiring edge already exists\n")
                next
            }
            
            points(gen.ppp$x[p], gen.ppp$y[p], cex=1.5, pch=21, bg="red")
            points(gen.ppp$x[q], gen.ppp$y[q], cex=1.5, pch=21, bg="red")
            
            pq =  c(gen.ppp$x[p], gen.ppp$y[p],
                    gen.ppp$x[q], gen.ppp$y[q],
                    p, q,
                    sqrt( ((gen.ppp$x[p]-gen.ppp$x[q])^2) + ((gen.ppp$y[p]-gen.ppp$y[q])^2) ),
                    calcAngle(c(gen.ppp$x[p], gen.ppp$y[p], gen.ppp$x[q], gen.ppp$y[q])),
                    0)
            
            el1 = which( (network_extra1$ind1==hd_v & network_extra1$ind2==p) |
                             (network_extra1$ind2==hd_v & network_extra1$ind1==p))
            
            el2 = which( (network_extra1$ind1==hd_v & network_extra1$ind2==q) |
                             (network_extra1$ind2==hd_v & network_extra1$ind1==q))
            
            if(length(el1)==0 | length(el2)==0){
                next
            }
            
            l = selectEdge(c(el1, el2), network_extra1, orgKDE_edge_feat, triKDE_edge_feat)
            
            if(l == el1){
                v1 = hd_v
                v2 = p
                el = el1
            }else if(l == el2){
                v1 = hd_v
                v2 = q
                el = el2
            }
            
            #### if any of the end vertices of the selected edge has degree 3,
            #### and if that is on the boundary we won't allow it
            #### cause in the end after deleting the boundary edges it will disconnect the network
            if((vertex_dist_boundary[v1]==0 & g2_degree[v1]<=3 & !isCornerV(v1, gen.ppp)) |
               (vertex_dist_boundary[v2]==0 & g2_degree[v2]<=3  & !isCornerV(v2, gen.ppp))){
                cat("Rewiring edge kept [Boundary degree constraint]\n")
                next
            }
            
            #### if v1 and v2 both are boundary vertices, we keep the edge for now
            if((vertex_dist_boundary[v1]==0) & (vertex_dist_boundary[v2]==0) & 
               ((gen.ppp$x[v1]==gen.ppp$x[v2])| (gen.ppp$y[v1]==gen.ppp$y[v2]) ) ){
                cat("Rewiring edge kept [Boundary edge constraint]\n")
                next
            }
            
            #### just to see what happens
            if(g2_degree[hd_v] <= 3){
                cat("Rewiring edge internal degree constraint\n")
                next
            }
            
            network_extra2 = rbind(network_extra1, pq)
            #### recompute face list
            face_list = computeFaceList(gen.ppp, network_extra2)
            
            temp_network_extra1 = network_extra2[-c(el), ]
            
            temp_graph_obj = make_empty_graph() %>% add_vertices(gen.ppp$n)
            temp_graph_obj = add_edges(as.undirected(temp_graph_obj),
                                       as.vector(t(as.matrix(temp_network_extra1[,5:6]))))
            
            if(is_connected(temp_graph_obj)){
                #### finding out the faces the chosen edge is a part of
                face_index = which(unlist(lapply(face_list,
                                                 function(x) isEdgeOnFace(x, c(network_extra2[el, ]$ind1,
                                                                               network_extra2[el, ]$ind2)))))
                
                #### recompute the face features and KDE assuming the edge has been removed
                #### and there's a new face now
                temp_g_o <- as_graphnel(temp_graph_obj)
                boyerMyrvoldPlanarityTest(temp_g_o)
                
                temp_face_list = planarFaceTraversal(temp_g_o)
                temp_face_node_count = sapply(temp_face_list, length)
                
                #### applying the shoe lace formula
                temp_face_area_list = sapply(temp_face_list, function(x) faceArea(x, gen.ppp))
                
                #### eliminating the outer face, it has the largest face area
                temp_face_node_count = temp_face_node_count[-which.max(temp_face_area_list)]
                temp_face_list = temp_face_list[-which.max(temp_face_area_list)]
                temp_face_area_list = temp_face_area_list[-which.max(temp_face_area_list)]
                
                #### face features computation
                temp_columns = c("Area_CF", "Perim.", "Ext.", "Disp.", "Elong.", "Eccentr.", "Orient.") # Area_CF: from contour function
                temp_face_features = data.frame(matrix(nrow = 0, ncol = length(temp_columns)))
                colnames(temp_face_features) = temp_columns
                
                temp_face_convexity_list = c()
                
                for(f in c(1: length(temp_face_list))){
                    temp_f_feat = computeFacefeatures(f, temp_face_list, gen.ppp, NULL)
                    temp_face_features = rbind(temp_face_features, temp_f_feat)
                    
                    temp_face_convexity_list = c(temp_face_convexity_list, computeFaceConvexity(temp_face_list[[f]], gen.ppp))
                    
                }# loop ends for each face of the temp network
                temp_face_features$Node_Count = temp_face_node_count
                
                temp_triKDE_face_feat_1 = kde(as.matrix(data.frame(temp_face_area_list, temp_face_features$orient)))
                temp_triKDE_face_feat_2 = kde(as.matrix(data.frame((temp_face_node_count))), 
                                              h=density((temp_face_node_count))$bw)
                temp_triKDE_edge_feat = kde(as.matrix(data.frame((temp_network_extra1$anglecomp),
                                                                 (temp_network_extra1$euclidDist))))
                
                ####finding the new face
                face_p = c()
                for(f_i in face_index){
                    face_p = c(face_p, as.numeric(unlist(face_list[f_i])))
                }
                face_p = unique(face_p)
                face_p_index = which(sapply(lapply(temp_face_list, function(x) sort(as.numeric(unlist(x)))),
                                            identical, sort(face_p)))
                
                if(length(face_p_index)==0){
                    cat("\nError in face identification 2\n")
                    next
                }
                
                if(temp_face_convexity_list[face_p_index] < org_face_convexity_mean){ # another option: mean(temp_face_convexity_list) < org_face_convexity_mean 
                    cat("\nEdge kept [Face convexity constraint]\n")
                    next
                }else{
                    after_elim_0 = eliminateEdges(gen.ppp, network_extra2, el)
                }
                
                
            }else{
                after_elim_0 = eliminateEdges(gen.ppp, network_extra2, NULL)
            }
            
            noChange = after_elim_0[[1]]
            network_extra1 = after_elim_0[[2]]
            g2_degree = after_elim_0[[3]]
            face_list = after_elim_0[[4]]
            face_area_list = after_elim_0[[5]]
            face_node_count = after_elim_0[[6]]
            triKDE_face_feat_1 = after_elim_0[[7]]
            triKDE_face_feat_2 = after_elim_0[[8]]
            triKDE_edge_feat = after_elim_0[[9]]
            tri_face_features = after_elim_0[[10]]
            tri_face_convexity_mean = after_elim_0[[11]]
            
            graph_obj =  make_empty_graph() %>% add_vertices(gen.ppp$n)
            graph_obj = add_edges(as.undirected(graph_obj),
                                  as.vector(t(as.matrix(network_extra1[,5:6]))))
            
            #### Transitivity measures the probability that the adjacent vertices of a vertex are connected.
            #### This is sometimes also called the clustering coefficient.
            cluster_coeff_s = igraph::transitivity(graph_obj, type = "global")
            cat("\nCC Sim: ", cluster_coeff_s, "\n")
            
            #### construct and display as corresponding ppp and linnet
            degs = igraph::degree(graph_obj, mode="total")
            # ord = order(as.numeric(names(degs)))
            # degs = degs[ord]
            
            #### attach the degree information to the point pattern for proper visualization
            marks(gen.ppp) = factor(degs)
            gen.ppp$markformat = "factor"
            g_o_lin = linnet(gen.ppp, edges=as.matrix(network_extra1[,5:6]))
            branch.lpp_s = lpp(gen.ppp, g_o_lin )
            
            plot(branch.lpp_s, main="Rewired", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange",
                                                                     "dodgerblue", "white", "maroon1",
                                                                     "mediumpurple", "yellow", "cyan"))
            
            break
        }
    }
    
    after_elim_0 = eliminateEdges(gen.ppp, network_extra1, NULL)
    noChange = after_elim_0[[1]]
    network_extra1 = after_elim_0[[2]]
    g2_degree = after_elim_0[[3]]
    face_list = after_elim_0[[4]]
    face_area_list = after_elim_0[[5]]
    face_node_count = after_elim_0[[6]]
    triKDE_face_feat_1 = after_elim_0[[7]]
    triKDE_face_feat_2 = after_elim_0[[8]]
    triKDE_edge_feat = after_elim_0[[9]]
    tri_face_features = after_elim_0[[10]]
    tri_face_convexity_mean = after_elim_0[[11]]
        
    return(list(noChange, network_extra1, g2_degree, face_list, face_area_list, face_node_count,
                triKDE_face_feat_1, triKDE_face_feat_2, triKDE_edge_feat, tri_face_features, tri_face_convexity_mean))
    
}


#### given a list of tentative edges to reject, and the edge feature estimations of the
#### original and the current networks, selects an edge to reject
selectEdge <- function(tent_edges, network_extra1, orgKDE_edge_feat, triKDE_edge_feat){
    set.seed(Sys.time())
    tent_edges = sample(tent_edges) #shuffle
    
    for (te in tent_edges) {
        org_edge_est = predict(orgKDE_edge_feat, x=c((network_extra1$anglecomp)[te], (network_extra1$euclidDist)[te]))
        tri_edge_est = predict(triKDE_edge_feat, x=c((network_extra1$anglecomp)[te], (network_extra1$euclidDist)[te]))

        if(tri_edge_est  > org_edge_est){
            return(te)
        }
    }
    
    #### if none of them satisfies the criterion, return a random one of them
    cat("Random edge selected\n")
    
    set.seed(Sys.time())
    return(sample(tent_edges, 1))
}


#### given a list of tentative edges to reject, and the edge feature estimations of the
#### original and the current networks, selects many edge to reject
selectMultEdges <- function(tent_edges, network_extra1, orgKDE_edge_feat, triKDE_edge_feat){
    set.seed(Sys.time())
    tent_edges = sample(tent_edges) #shuffle
    
    e_list = c()
    for (te in tent_edges) {
        org_edge_est = predict(orgKDE_edge_feat, x=c((network_extra1$anglecomp)[te], (network_extra1$euclidDist)[te]))
        tri_edge_est = predict(triKDE_edge_feat, x=c((network_extra1$anglecomp)[te], (network_extra1$euclidDist)[te]))

        if(tri_edge_est  > org_edge_est){
            e_list = c(e_list, te)
        }
    }
    
    if(length(e_list) == 0){
        #### if none of them satisfies the criterion, return a random one of them
        cat("Random edge selected\n")
        return(tent_edges)
    }else{
        return(e_list)
    }
}

selectMultEdges2 <- function(tent_edges, network_extra1, orgKDE_edge_feat, triKDE_edge_feat){
    set.seed(Sys.time())
    return(sample(tent_edges)) #shuffle and return
    
}


#### detects if a vertex is a corner of the pp boundary
isCornerV <- function(v, gen.ppp){
    #cat("!!!!!!!!!!!!!!!!!!!!!!!!Inside isCornerV 1!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    v_x = gen.ppp$x[v]
    v_y = gen.ppp$y[v]
    
    cond1 = (v_x %in% gen.ppp$window$xrange)
    cond2 = (v_y %in% gen.ppp$window$yrange)
    
    #cat("!!!!!!!!!!!!!!!!!!!!!!!!Inside isCornerV 2!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    
    if(is.na(cond1) | is.na(cond2)){
        return(FALSE)
    }else{
        if(cond1 & cond2){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }
}


deterministicEdges_3 <- function(gen.ppp, branch.ppp, branch.all, org_face_feature, sample_id, org_face_convexity_mean, org_face_convexity_sd){
    #### construct the Delaunay triangulation on the parent points as a starter network
    #### ord_point_list: to maintain the order of the points
    ord_point_list = data.frame(x = gen.ppp$x, y = gen.ppp$y)
    
    network_triangulation = deldir(ord_point_list[, 1:2])
    
    #### construct a convenient data structure to keep the triangulation information
    network_extra1 = data.frame(x1=network_triangulation$delsgs$x1, y1=network_triangulation$delsgs$y1,
                                x2=network_triangulation$delsgs$x2, y2=network_triangulation$delsgs$y2,
                                ind1=network_triangulation$delsgs$ind1, ind2=network_triangulation$delsgs$ind2)
    
    network_extra1$euclidDist = apply(network_extra1, 1, function(x) sqrt( ((x[1]-x[3])^2) + ((x[2]-x[4])^2) ) ) 
    network_extra1$anglecomp = apply(network_extra1, 1, function(x) calcAngle(x))
    
    ###################################################################################
    #### distance of each vertex from the point pattern boundary
    vertex_dist_boundary = bdist.points(gen.ppp)
    
    ####index of the boundary edges
    bb_edges = which((vertex_dist_boundary[network_extra1$ind1]==0) & 
                         (vertex_dist_boundary[network_extra1$ind2]==0) &
                         ((gen.ppp$x[network_extra1$ind1]==gen.ppp$x[network_extra1$ind2])|
                              (gen.ppp$y[network_extra1$ind1]==gen.ppp$y[network_extra1$ind2]) ) )
    
    #### eliminate the edges from the triangulation whose length is larger/smaller than the max/min edge length
    #### in the original network
    #### then recompute the face features and KDE required
    org_max_edge_length = max(branch.all$euclid)
    org_min_edge_length = min(branch.all$euclid)
    edges_to_eliminate = c(which(network_extra1$euclidDist > org_max_edge_length),
                           which(network_extra1$euclidDist < org_min_edge_length))
    
    #### but we do not eliminate any boundary edge [for now]
    edges_to_eliminate = setdiff(edges_to_eliminate, bb_edges) 
    
    cat("Eliminating larger and smaller edges\n")
    after_elim = eliminateEdges(gen.ppp, network_extra1, edges_to_eliminate)
    
    noChange = after_elim[[1]]
    network_extra1 = after_elim[[2]]
    g2_degree = after_elim[[3]]
    face_list = after_elim[[4]]
    face_area_list = after_elim[[5]]
    face_node_count = after_elim[[6]]
    triKDE_face_feat_1 = after_elim[[7]]
    triKDE_face_feat_2 = after_elim[[8]]
    triKDE_edge_feat = after_elim[[9]]
    face_features = after_elim[[10]]
    face_convexity_mean = after_elim[[11]]
    
    cat("tri avg face convexity:", face_convexity_mean, "\n")
    ###################################################################################
    
    #### new
    #### Reminder: The initial triangulation already has closed faces at  the boundary, 
    #### no need to add extra edges before computing the faces.
    #### compute the face feature of the triangulation
    graph_obj =  graph_from_data_frame(unique(network_extra1[, 5:6]), directed = FALSE) 
    
    #### Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
    #### This is sometimes also called the clustering coefficient.
    cluster_coeff_t = igraph::transitivity(graph_obj, type = "global")
    cat("CC DT: ", cluster_coeff_t, "\n")
    
    #### create another copy of the data structure
    network_extra = rbind(network_extra1)
    
    #### construct and display as corresponding ppp and linnet
    degs = (igraph::degree(graph_obj, mode="total"))
    ord = order(as.numeric(names(degs)))
    degs = degs[ord]
    
    #### attach the degree information to the point pattern for proper visualization
    marks(gen.ppp) = factor(degs)
    gen.ppp$markformat = "factor"
    g_o_lin = linnet(gen.ppp, edges=as.matrix(network_extra[,5:6]))
    branch.lpp_dt = lpp(gen.ppp, g_o_lin )
    
    ####new
    ####plot the feature densities of the initial triangulation vs the original network for comparison
    comparePlotOrgSim(org_face_feature, face_features, branch.all, network_extra)
    ####new end
    
    plot(branch.lpp_dt, main="Initial DT", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", 
                                                                        "dodgerblue", "white", "maroon1", 
                                                                        "mediumpurple", "yellow", "cyan"))
    
    
    svglite(paste("D:/Fall 2023/Research/Prelim/figures/init tri figures/init_dt.svg", sep=""), width = 8, height = 6)
    plot(branch.lpp_dt, main="", pch=21, cex=1, bg=c("black", "red3", "green3", "orange", 
                                                                        "dodgerblue", "white", "maroon1", 
                                                                        "mediumpurple", "yellow", "cyan"))
    dev.off()                                                                    
    
    #### degree list of the constructed Delaunay graph again for edge weight calculation
    g_o_degree = degs
    
    network_extra$weight = apply(network_extra, 1, function(x) computeEdgeWeight(g_o_degree, x))
    network_extra$weight = range01(network_extra$weight)
    
    return(list(network_extra, face_list, face_area_list, face_node_count, triKDE_face_feat_1, triKDE_face_feat_2, triKDE_edge_feat, g_o_degree, face_features))
}


rejectionSampling_3 <- function(gen.ppp, branch.ppp, branch.all, org_face_feature, network_extra1, face_list, face_area_list, face_node_count, 
                    g2_degree, orgKDE_face_feat_1, orgKDE_face_feat_2, triKDE_face_feat_1, triKDE_face_feat_2, orgKDE_edge_feat, triKDE_edge_feat,
                    meshedness, network_density, compactness, cluster_coeff, org_max_deg,
                    sample_id, tri_face_features, org_face_convexity_mean, org_face_convexity_sd){
    set.seed(Sys.time())
    
    #### distance of each vertex from the point pattern boundary
    vertex_dist_boundary = bdist.points(gen.ppp)
    
    org_max_edge_length = max(branch.all$euclid)
    org_min_edge_length = min(branch.all$euclid)
    
    #### main loop
    noChange = 0
    loop_count = -1
    while (TRUE) {
        # q = readline()
        # if(q=="q"){
        #     break
        # }
        
        cat("\n--------------------\nnoChange value: ", noChange, "\n")
        cat("Loop count: ", loop_count, "\n")
        if(noChange >= 500){    # if the network has not been changed for certain iterations
            cat("\nNo edges rejected for 500 iterations.\n")
            break
        }
        
        # if((loop_count!=0) & ((loop_count %% 100)==0)){
        #     cat("...Rewiring...\n")
        # 
        #     lst = sort(g2_degree, index.return=TRUE, decreasing=TRUE)
        #     high_deg_vertices = (lapply(lst, `[`, lst$x %in% head(unique(lst$x), 1)))$ix
        #     if(length(high_deg_vertices) > 1){
        #         high_deg_vertices = sample(high_deg_vertices, 1)
        #     }
        # 
        #     after_rewire = rewireEdges(gen.ppp, network_extra1, high_deg_vertices, g2_degree, face_list, org_face_convexity_mean,
        #                                orgKDE_edge_feat, triKDE_edge_feat)
        # 
        #     network_extra1 = after_rewire[[2]]
        #     g2_degree = after_rewire[[3]]
        #     face_list = after_rewire[[4]]
        #     face_area_list = after_rewire[[5]]
        #     face_node_count = after_rewire[[6]]
        #     triKDE_face_feat_1 = after_rewire[[7]]
        #     triKDE_face_feat_2 = after_rewire[[8]]
        #     triKDE_edge_feat = after_rewire[[9]]
        #     tri_face_features = after_rewire[[10]]
        #     face_convexity_mean = after_rewire[[11]]
        #     
        #     comparePlotOrgSim2(org_face_feature, tri_face_features, branch.all, network_extra1)
        #     pause = readline()
        # 
        # }
        
        #### cc of the current network
        cc_cur = ccFromDataframe(gen.ppp, network_extra1)
        if(cc_cur <= cluster_coeff){
            cat("\nRejection sampling ended [CC reached target]\n")
            break
        }

        #### other network metrics of the current network
        metrics = netMetrics(gen.ppp, network_extra1)
        mesh = metrics[[1]]
        net_den = metrics[[2]]
        compact = metrics[[3]]
        if((mesh <= meshedness) | (net_den <= network_density) | (compact <= compactness)){
            cat("\nRejection sampling ended [Net metrics reached target]\n")
            break
        }
        
        if(loop_count <= 250){
            #### prepare the vertex probability from degree values
            prob_vertex = computeVertexProb(org_max_deg, g2_degree, network_extra1)
    
        }else if(loop_count > 250){
            #### prepare the vertex probability from degree values
            prob_vertex = computeVertexProb2(org_max_deg, g2_degree, network_extra1)
            
        }
        #### select a vertex at random or based on high degree
        selected_vertex = sample.int(gen.ppp$n, 1, prob = prob_vertex)
        cat("\nSelected vertex ID: ", selected_vertex, ", Degree of the selected vertex: ", g2_degree[selected_vertex], "\n")
        
        #### detect the neighbors of a given vertex
        adj_vertices_from_df = c(network_extra1$ind1[network_extra1$ind2==selected_vertex],
                                 network_extra1$ind2[network_extra1$ind1==selected_vertex])
        
        #### list of the edges that are between those neighboring vertices (if any)
        #### the indices are of network_extra1
        edge_bet_adj_vertices = which((network_extra1$ind1 %in% adj_vertices_from_df) &
                                          (network_extra1$ind2 %in% adj_vertices_from_df))
        
        tentative_edges = which( (network_extra1$ind1==selected_vertex) | (network_extra1$ind2==selected_vertex) )
        tentative_edges = unique( c(tentative_edges, edge_bet_adj_vertices) )
        
        selected_edges = selectEdge(tentative_edges, network_extra1, orgKDE_edge_feat, triKDE_edge_feat)
        
        if(length(selected_edges) == 0){
            cat("No edges selected based on edge length and/ or angle\n")
            quit()
        }
          
        for(selected_edge in selected_edges){
            loop_count = loop_count + 1
            
            cat("Selected edge ID: ", selected_edge, "\n\n")
            #print(network_extra1[selected_edge, ])
            
            if(is.na(network_extra1$ind1[selected_edge])){
                cat("Error in finding edge\n")
                next
            }
    
            #### check for the degree of the end vertices
            #### degree-one vertices only at the boundary
            v1 = network_extra1$ind1[selected_edge]
            v2 = network_extra1$ind2[selected_edge]
            lines(c(gen.ppp$x[v1], gen.ppp$x[v2]), c(gen.ppp$y[v1], gen.ppp$y[v2]), col="red", lwd=2.6)
    
            #### if any of the end vertices of the selected edge has degree 3,
            #### and if that is on the boundary we won't allow it
            #### cause in the end after deleting the boundary edges it will disconnect the network
            if((vertex_dist_boundary[v1]==0 & g2_degree[v1]<=3 & !isCornerV(v1, gen.ppp)) |
               (vertex_dist_boundary[v2]==0 & g2_degree[v2]<=3  & !isCornerV(v2, gen.ppp))){
                noChange = noChange + 1
                cat("\nEdge kept [Boundary degree constraint]\n")
                next
            }
    
            #### if v1 and v2 both are boundary vertices, we keep the edge for now
            if((vertex_dist_boundary[v1]==0) & (vertex_dist_boundary[v2]==0) & 
               ((gen.ppp$x[v1]==gen.ppp$x[v2])| (gen.ppp$y[v1]==gen.ppp$y[v2]) ) ){
                noChange = noChange + 1
                cat("\nEdge kept [Boundary edge constraint]\n")
                next
            }
    
            #### just to see what happens
            if(g2_degree[v1]<=2 | g2_degree[v2]<=2){
                noChange = noChange + 1
                cat("\nEdge kept [Internal degree constraint]\n")
                next
            }
    
            #### check if removing that edge disconnects the network
            temp_network_extra1 = network_extra1[-c(selected_edge), ]

            temp_graph_obj = make_empty_graph() %>% add_vertices(gen.ppp$n)
            temp_graph_obj = add_edges(as.undirected(temp_graph_obj),
                                       as.vector(t(as.matrix(temp_network_extra1[,5:6]))))
    
            if(is_connected(temp_graph_obj)){
                #### ensured that removing the edge will not disconnect the network
                #### now check for distribution of face features
    
                #### finding out the faces the chosen edge is a part of
                face_index = which(unlist(lapply(face_list,
                                                 function(x) isEdgeOnFace(x, c(network_extra1[selected_edge, ]$ind1,
                                                                               network_extra1[selected_edge, ]$ind2)))))
    
                #### recompute the face features and KDE assuming the edge has been removed
                #### and there's a new face now
                temp_g_o <- as_graphnel(temp_graph_obj)
                boyerMyrvoldPlanarityTest(temp_g_o)
    
                temp_face_list = planarFaceTraversal(temp_g_o)
                temp_face_node_count = sapply(temp_face_list, length)
    
                #### applying the shoe lace formula
                temp_face_area_list = sapply(temp_face_list, function(x) faceArea(x, gen.ppp))
    
                #### eliminating the outer face, it has the largest face area
                temp_face_node_count = temp_face_node_count[-which.max(temp_face_area_list)]
                temp_face_list = temp_face_list[-which.max(temp_face_area_list)]
                temp_face_area_list = temp_face_area_list[-which.max(temp_face_area_list)]
    
                #### face features computation
                temp_columns = c("Area_CF", "Perim.", "Ext.", "Disp.", "Elong.", "Eccentr.", "Orient.") # Area_CF: from contour function
                temp_face_features = data.frame(matrix(nrow = 0, ncol = length(temp_columns)))
                colnames(temp_face_features) = temp_columns
                
                temp_face_convexity_list = c()
    
                for(f in c(1: length(temp_face_list))){
                    temp_f_feat = computeFacefeatures(f, temp_face_list, gen.ppp, NULL)
                    temp_face_features = rbind(temp_face_features, temp_f_feat)
                    
                    temp_face_convexity_list = c(temp_face_convexity_list, computeFaceConvexity(temp_face_list[[f]], gen.ppp))
    
                }# loop ends for each face of the temp network
                temp_face_features$Node_Count = temp_face_node_count
    
                temp_triKDE_face_feat_1 = kde(as.matrix(data.frame(temp_face_area_list, temp_face_features$orient)))
                temp_triKDE_face_feat_2 = kde(as.matrix(data.frame((temp_face_node_count))), 
                                              h=density((temp_face_node_count))$bw)
                temp_triKDE_edge_feat = kde(as.matrix(data.frame((temp_network_extra1$anglecomp),
                                                                 (temp_network_extra1$euclidDist))))
    
                ####finding the new face
                face_p = c()
                for(f_i in face_index){
                    face_p = c(face_p, as.numeric(unlist(face_list[f_i])))
                }
                face_p = unique(face_p)
                face_p_index = which(sapply(lapply(temp_face_list, function(x) sort(as.numeric(unlist(x)))),
                                            identical, sort(face_p)))
    
                if(length(face_p_index)==0){
                    cat("\nError in face identification 2\n")
                    next
                }
                
                edge_reject = FALSE
                epsilon_f = 0
                epsilon_e = 0
                convexity_param = org_face_convexity_mean
                
                if(loop_count < 100){
                    epsilon_f = 1e-09
                    convexity_param = org_face_convexity_mean

                }else if(loop_count >= 100 & loop_count < 250){
                    cat("Epsilon increment 1\n")
                    epsilon_f = 1.5e-08
                    convexity_param = org_face_convexity_mean - (org_face_convexity_sd/2)

                }else if(loop_count >= 250){
                    cat("Epsilon increment 2\n")
                    epsilon_f = 2e-07
                    convexity_param = org_face_convexity_mean - org_face_convexity_sd
                }
                
                #### prediction
                org_est_1 = predict(orgKDE_face_feat_1, x=c(temp_face_area_list[face_p_index], temp_face_features$orient[face_p_index]))
                org_est_2 = predict(orgKDE_face_feat_2, x=c((temp_face_node_count)[face_p_index]))
                
                temp_tri_est_1 = predict(temp_triKDE_face_feat_1, x=c(temp_face_area_list[face_p_index], temp_face_features$orient[face_p_index]))
                temp_tri_est_2 = predict(temp_triKDE_face_feat_2, x=c((temp_face_node_count)[face_p_index]))
                
                cat("\nFace est. 1 diff: ", (org_est_1 - temp_tri_est_1), "\n")
                # cat("\nFace est. 2 diff: ", (org_est_2 - temp_tri_est_2), "\n")
                 
                cat("Face convexity of the new face: ", temp_face_convexity_list[face_p_index], "\n")
                 
                # org_edge_est = predict(orgKDE_edge_feat, x=network_extra1$anglecomp[selected_edge])
                # tri_edge_est = predict(temp_triKDE_edge_feat, x=network_extra1$anglecomp[selected_edge])
                # 
                # cat("\nEdge est. diff: ", (tri_edge_est - org_edge_est), "\n")
    
                if(temp_face_convexity_list[face_p_index] < convexity_param){ # another option: mean(temp_face_convexity_list) < org_face_convexity_mean 
                    #### keep the edge, no change
                    noChange = noChange + 1
                    cat("\nEdge kept [Face convexity constraint]\n")
                }else{
                    if(length(face_index)==2){
                        f1 = face_index[1]
                        f2 = face_index[2]
        
                        if((org_est_1+epsilon_f >= temp_tri_est_1) ){
                            edge_reject = TRUE
                        }
        
                    }else if(length(face_index)==1){
                        f1 = face_index[1]
        
                        if((org_est_1+epsilon_f >= temp_tri_est_1) ){
                            edge_reject = TRUE
                        }
        
                    }else{
                        cat("\nError in face identification 1\n")
                        quit()
                    }
                }
                
    
                if(edge_reject){
                    #### remove the edge, there is a change
                    noChange = 0
                    cat("\nEdge deleted\n")
    
                    #### make necessary changes permanent
                    network_extra1 = temp_network_extra1
                    g2_degree = igraph::degree(temp_graph_obj, mode="total")
                    print(table(g2_degree))
    
                    face_list = temp_face_list
                    face_area_list = temp_face_area_list
                    face_node_count = temp_face_node_count
                    triKDE_face_feat_1 = temp_triKDE_face_feat_1
                    triKDE_face_feat_2 = temp_triKDE_face_feat_2
                    triKDE_edge_feat =  temp_triKDE_edge_feat
                    tri_face_features = temp_face_features
    
                }else{
                    #### keep the edge, no change
                    noChange = noChange + 1
                    cat("\nEdge kept [Face feature and/or edge estimation constraint]\n")
                }
    
            }else{
                #### keep the edge, no change
                noChange = noChange + 1
                cat("\nEdge kept [Connectivity constraint]\n")
            }
    
            #### temporarily plotting each iteration
            graph_obj =  make_empty_graph() %>% add_vertices(gen.ppp$n)
            graph_obj = add_edges(as.undirected(graph_obj),
                                  as.vector(t(as.matrix(network_extra1[,5:6]))))
    
            #### Transitivity measures the probability that the adjacent vertices of a vertex are connected.
            #### This is sometimes also called the clustering coefficient.
            cluster_coeff_s = igraph::transitivity(graph_obj, type = "global")
            cat("\nCC Sim: ", cluster_coeff_s, "\n")
    
            #### construct and display as corresponding ppp and linnet
            degs = igraph::degree(graph_obj, mode="total")
            # ord = order(as.numeric(names(degs)))
            # degs = degs[ord]
    
            #### attach the degree information to the point pattern for proper visualization
            marks(gen.ppp) = factor(degs)
            gen.ppp$markformat = "factor"
            g_o_lin = linnet(gen.ppp, edges=as.matrix(network_extra1[,5:6]))
            branch.lpp_s = lpp(gen.ppp, g_o_lin )
    
            plot(branch.lpp_s, main="Sim", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange",
                                                                        "dodgerblue", "white", "maroon1",
                                                                        "mediumpurple", "yellow", "cyan"))
        }
    }
    
    #### compute index of the boundary edges again
    bb_edges_2 = which((vertex_dist_boundary[network_extra1$ind1]==0) & 
                         (vertex_dist_boundary[network_extra1$ind2]==0) & 
                           ((gen.ppp$x[network_extra1$ind1]==gen.ppp$x[network_extra1$ind2])|
                                 (gen.ppp$y[network_extra1$ind1]==gen.ppp$y[network_extra1$ind2]) ) )
    
    #### eliminate boundary-boundary edges
    cat("Eliminating boundary edges\n")
    after_elim_0 = eliminateEdges(gen.ppp, network_extra1, bb_edges_2)

    noChange = after_elim_0[[1]]
    network_extra1 = after_elim_0[[2]]
    g2_degree = after_elim_0[[3]]
    face_list = after_elim_0[[4]]
    face_area_list = after_elim_0[[5]]
    face_node_count = after_elim_0[[6]]
    triKDE_face_feat_1 = after_elim_0[[7]]
    triKDE_face_feat_2 = after_elim_0[[8]]
    triKDE_edge_feat = after_elim_0[[9]]
    tri_face_features = after_elim_0[[10]]
    face_convexity_mean = after_elim_0[[11]]
    
    #### figure out the vertex id of the corner points
    c_v = c()
    for (i in c(1:gen.ppp$n)) {
        if(isCornerV(i, gen.ppp)){
            c_v = c(c_v, i)
        }
    }
    
    #### compute index of the corner edges 
    c_e = which(network_extra1$ind1 %in% c_v | network_extra1$ind2 %in% c_v)
    
    #### remove corner points any edge incident on them
    gen.ppp_2 = subset.ppp(gen.ppp, !((x==gen.ppp$window$xrange[1] | x==gen.ppp$window$xrange[2]) 
                                     & (y==gen.ppp$window$yrange[1] | y==gen.ppp$window$yrange[2])))
    
    #### eliminate corner edges
    cat("Eliminating corner edges\n")
    after_elim_1 = eliminateEdges(gen.ppp, network_extra1, c_e)
    
    noChange = after_elim_1[[1]]
    network_extra1 = after_elim_1[[2]]
    g2_degree = after_elim_1[[3]]
    face_list = after_elim_1[[4]]
    face_area_list = after_elim_1[[5]]
    face_node_count = after_elim_1[[6]]
    triKDE_face_feat_1 = after_elim_1[[7]]
    triKDE_face_feat_2 = after_elim_1[[8]]
    triKDE_edge_feat = after_elim_1[[9]]
    tri_face_features = after_elim_1[[10]]
    face_convexity_mean = after_elim_1[[11]]
    
    #### final simulated network
    graph_obj =  make_empty_graph() %>% add_vertices(gen.ppp_2$n)
    graph_obj = add_edges(as.undirected(graph_obj), 
              as.vector(t(as.matrix(network_extra1[,5:6]))))
    
    #### Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
    #### This is sometimes also called the clustering coefficient.
    cluster_coeff_s = igraph::transitivity(graph_obj, type = "global")
    cat("CC Sim: ", cluster_coeff_s, "\n")
    
    #### construct and display as corresponding ppp and linnet
    degs = igraph::degree(graph_obj, mode="total")
    # ord = order(as.numeric(names(degs)))
    # degs = degs[ord]
    
    #### attach the degree information to the point pattern for proper visualization
    marks(gen.ppp_2) = factor(degs)
    gen.ppp_2$markformat = "factor"
    g_o_lin = linnet(gen.ppp_2, edges=as.matrix(network_extra1[,5:6]))
    branch.lpp_s = lpp(gen.ppp_2, g_o_lin )
    
    plot(branch.lpp_s, main="Sim", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", 
                                                                        "dodgerblue", "white", "maroon1", 
                                                                        "mediumpurple", "yellow", "cyan"))
    return(list(gen.ppp_2, network_extra1))                                    
}


generateNetworkEdges_3 <- function(gen.ppp, branch.ppp, branch_all, org_face_feature, orgKDE_face_feat_1, orgKDE_face_feat_2, orgKDE_edge_feat,
                                   meshedness, network_density, compactness, cluster_coeff, org_max_deg,
                                   sample_id, org_face_convexity_mean, org_face_convexity_sd){
    
    #### constructing the deterministic Delaunay triangulation as the initial ganglionic network
    triangulation_info_list = deterministicEdges_3(gen.ppp, branch.ppp, branch.all, org_face_feature, sample_id, org_face_convexity_mean, org_face_convexity_sd)
    
    #### returned values
    network_extra1 = triangulation_info_list[[1]]
    
    face_list = triangulation_info_list[[2]]
    face_area_list = triangulation_info_list[[3]]
    face_node_count = triangulation_info_list[[4]]
    
    triKDE_face_feat_1 = triangulation_info_list[[5]]
    triKDE_face_feat_2 = triangulation_info_list[[6]]
    triKDE_edge_feat = triangulation_info_list[[7]]
    
    g2_degree = triangulation_info_list[[8]]
    tri_face_features = triangulation_info_list[[9]]
    
    #### remove edges from the initial triangulation by rejection sampling
    sampled_net = rejectionSampling_3(gen.ppp, branch.ppp, branch.all, org_face_feature, network_extra1, face_list, face_area_list, face_node_count, 
                                        g2_degree, orgKDE_face_feat_1, orgKDE_face_feat_2, triKDE_face_feat_1, triKDE_face_feat_2, orgKDE_edge_feat, triKDE_edge_feat,
                                        meshedness, network_density, compactness, cluster_coeff, org_max_deg,
                                        sample_id, tri_face_features, org_face_convexity_mean, org_face_convexity_sd)
    
    gen.ppp = sampled_net[[1]]
    network_extra = sampled_net[[2]]
        
    #### create a graph from sampled triangulation
    g2 = make_empty_graph() %>% add_vertices(gen.ppp$n)
    g2 = add_edges(as.undirected(g2), as.vector(t(as.matrix(network_extra[,5:6]))))
    
    #### display as corresponding ppp and linnet
    g2_lin = linnet(gen.ppp, edges=as.matrix(network_extra[, 5:6]))
    
    return(list(gen.ppp, network_extra, g2_lin))
}


####main 
#### extracting parent directory information for accessing input and output location
dir = this.dir()
folder = strsplit(dir, "/")
folder = folder[[1]][length(folder[[1]])]
parent = strsplit(dir, folder)

face_folder = paste(parent, "Outputs/ENSMouse/FaceFeature/", sep="")
face_features_combined = read.csv(paste(face_folder, "FaceFeatures_3.csv", sep = ""))


#### the TIF images of the ganglionic networks are preprocessed in Fiji (ImageJ) and 
#### the network information is extracted as .csv files
branch_info_folder = paste(parent, "Data/ENSMouse Branch Information (in um) v2.0/", sep="")
branch_info_files = list.files(branch_info_folder, recursive = TRUE, pattern = "\\.csv", full.names = TRUE)

i = 2 # index of the ENS network we want to work on

ens_location = strsplit(branch_info_files[i], "/")[[1]][11]
sample_id = strsplit(strsplit(branch_info_files[i], "/")[[1]][11], "\\.")[[1]][1]
cat("\n(", i, ") Location: ", ens_location, "\nSample Id: ", sample_id, "\n")

max_y = 1 # 4539.812 found by computation; right now keeping everything unscaled as the moments can not be computed otherwise

data_struct_list = constructDataStruct(sample_id, parent, branch_info_files[i], output_folder_path, max_y)

#### the returned values
branch.all = data_struct_list[[1]]
branch.ppp = data_struct_list[[2]]  # marked point pattern, degree of the points as marks
branch.lpp = data_struct_list[[3]]
g1 = data_struct_list[[4]]
hardcoreStrauss_model_param = data_struct_list[[5]]

plot(branch.lpp, main="original", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", "dodgerblue", 
                                                               "white", "maroon1", "mediumpurple"))

svglite(paste("D:/Fall 2023/Research/Prelim/figures/org figures/org_net.svg", sep=""), width = 8, height = 6)
plot(branch.lpp, main="", pch=21, cex=1, bg=c("black", "red3", "green3", "orange", "dodgerblue", 
                                                        "white", "maroon1", "mediumpurple"))
dev.off()

#### alpha, gamma, psi (meshedness, network density and compactness parameters)
N = branch.ppp$n
E = length(branch.all$x1)
A = summary(branch.ppp)$window$area
L = sum(branch.all$euclid)

meshedness = (E-N+1)/((2*N)-5)
network_density = E/((3*N)-6)
compactness = 1- ((4*A)/(L-(2*sqrt(A)))^2)

cat("Meshedness: ", meshedness, ", Network density: ", network_density, ", Compactness: ", compactness, "\n")

#### Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
#### This is sometimes also called the clustering coefficient.
cluster_coeff = igraph::transitivity(g1, type = "global")
cat("CC original: ", cluster_coeff, "\n")

org_max_deg = max(igraph::degree(g1))

#### filter out the face face features of the sample under consideration
#### Reminder: the face features were computed assuming additional edges were computed 
#### to close the open faces at the boundary
#### similar thing should be done for the new network too (where needed).
face_feature = face_features_combined[face_features_combined$sample_id == sample_id, ]

org_face_convexity_mean = mean(face_feature$Convexity)
org_face_convexity_sd = sd(face_feature$Convexity)
cat("original avg face convexity: ", org_face_convexity_mean, " +- ", org_face_convexity_sd, "\n")

orgKDE_face_feat_1 = kde(as.matrix(data.frame(face_feature$Area_SL, face_feature$Orient.)))
orgKDE_face_feat_2 = kde(as.matrix(data.frame((face_feature$Node_Count) )),
                         h=density((face_feature$Node_Count))$bw)
orgKDE_edge_feat = kde(as.matrix(data.frame((apply(branch.all, 1, function(x) calcAngle(x))), 
                                            (branch.all$euclid) )))

####At this point, new point pattern will be simulated. 
####For now we are generating networks on the original point pattern.
####branch.ppp will be replaced by some new ppp object
gen.ppp = unmark(branch.ppp)
gen_corner.ppp = ppp(x=c(gen.ppp$window$xrange[1], gen.ppp$window$xrange[2], gen.ppp$window$xrange[2], gen.ppp$window$xrange[1]), 
                 y=c(gen.ppp$window$yrange[2], gen.ppp$window$yrange[2], gen.ppp$window$yrange[1], gen.ppp$window$yrange[1]),
                 window = gen.ppp$window)
gen.ppp = superimpose(gen.ppp, gen_corner.ppp)

#### call the network generation functions
network_info_list = generateNetworkEdges_3(gen.ppp, branch.ppp, branch_all, face_feature, orgKDE_face_feat_1, orgKDE_face_feat_2, orgKDE_edge_feat,
                                           meshedness, network_density, compactness, cluster_coeff, org_max_deg,
                                           sample_id, org_face_convexity_mean, org_face_convexity_sd)

#### returned values
gen.ppp_wo_c = network_info_list[[1]]
net_data_struct = network_info_list[[2]]
linnet_obj = network_info_list[[3]]

#### creating a list representing degree for every node an use it in spatstat pattern
graph_obj_0 = graph_from_data_frame(net_data_struct[, 5:6], directed = FALSE)
degs = (igraph::degree(graph_obj_0, mode="total"))
ord = order(as.numeric(names(degs)))
degs = degs[ord]

#### attach the degree information to the point pattern for proper visualization
marks(gen.ppp_wo_c) = factor(degs)
gen.ppp_wo_c$markformat = "factor"
branch.lpp_2 = lpp(gen.ppp_wo_c, linnet_obj )

svglite(paste("D:/Fall 2023/Research/Prelim/figures/sim figures/sim_net.svg", sep=""), width = 8, height = 6)
plot(branch.lpp_2, main="", pch=21, cex=1, bg=c("black", "red3", "green3", "orange", "dodgerblue", 
                                                                  "white", "maroon1", "mediumpurple", "yellow", "cyan"))
dev.off()                                                              

#### compute the face features of the newly constructed network to compare the density with the original face feature
#### Reminder: for the simulated network there are many open boundary faces
#### additional edges are required
new_edges = computeBoundaryEdges(gen.ppp_wo_c)

graph_obj =  graph_from_data_frame(unique(rbind(net_data_struct[, 5:6], new_edges)), directed = FALSE) # new graph object that combines the actual network and the new edges
graph_obj = delete_edges(graph_obj, which(which_multiple(graph_obj)))

degs = (igraph::degree(graph_obj, mode="total"))
ord = order(as.numeric(names(degs)))
degs = degs[ord]

#### attach the degree information to the point pattern for proper visualization
marks(gen.ppp) = factor(degs)
gen.ppp$markformat = "factor"
branch.lpp_3 = lpp(gen.ppp,  linnet(gen.ppp, edges=as.matrix(unique(rbind(net_data_struct[, 5:6], new_edges)))))

plot(branch.lpp_3, main="simulated(2)", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", "dodgerblue", 
                                                                  "white", "maroon1", "mediumpurple",  "yellow", "cyan"))

#### Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
#### This is sometimes also called the clustering coefficient.
cluster_coeff_s = igraph::transitivity(graph_obj, type = "global")
cat("CC simulated: ", cluster_coeff_s, "\n")

g_o <- as_graphnel(graph_obj) ## Convert igraph object to graphNEL object for planarity testing
boyerMyrvoldPlanarityTest(g_o)

face_list = planarFaceTraversal(g_o)
face_node_count = sapply(face_list, length)

#### applying the shoe lace formula
face_area_list = sapply(face_list, function(x) faceArea(x, gen.ppp))

#### eliminating the outer face, it has the largest face area
face_node_count = face_node_count[-which.max(face_area_list)]
face_list = face_list[-which.max(face_area_list)]
face_area_list = face_area_list[-which.max(face_area_list)]

#### face features computation
columns = c("Area_CF", "Perim.", "Ext.", "Disp.", "Elong.", "Eccentr.", "Orient.") # Area_CF: from contour function
face_features_sim = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(face_features_sim) = columns

face_convexity_list = c()

for(f in c(1: length(face_list))){
    f_feat = computeFacefeatures(f, face_list, gen.ppp, NULL)
    face_features_sim = rbind(face_features_sim, f_feat)
    
    face_convexity_list = c(face_convexity_list, computeFaceConvexity(face_list[[f]], gen.ppp))
    
}# loop ends for each face of the simulated network
face_features_sim$Node_Count = face_node_count

comparePlotOrgSim2(face_feature, face_features_sim, branch.all, net_data_struct)

#########################################################################
#(org_est1 < tri_est1) & (org_est2 < tri_est2) & (org_est >= temp_tri_est)