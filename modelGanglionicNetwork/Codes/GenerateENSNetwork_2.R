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

setwd("~/GitHub/spatial-neuro/modelGanglionicNetwork/Codes")
#### source the functions from other files
source("AnalyzeGanglionicNetwork.R")


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
    cat("face id: ", f, "\n")
    
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
    orient = 0.5 * atan2((2 * moments$value[moments$moment == "m11"]) , (moments$value[moments$moment == "m20"] - moments$value[moments$moment == "m02"]))
    orient = orient * 180 / pi
    eccentr = (((moments$value[moments$moment == "m02"] - moments$value[moments$moment == "m20"]) * (moments$value[moments$moment == "m02"] - moments$value[moments$moment == "m20"])) 
               + (4 * moments$value[moments$moment == "m11"] * moments$value[moments$moment == "m11"])) / moments$value[moments$moment == "m00"]
    
    return(data.frame(area, perim, ext, disp, elong, eccentr, orient))
}


#### plot the face features of two given networks for comparison 
#### can be original vs. initial DT or original vs. final network
comparePlotOrgSim <- function(org_face_feature, face_features){
    cat("Number of faces in the original network: ", length(org_face_feature$X))
    cat("\nNumber of faces in the simulated network: ", length(face_features$area))
    
    # area
    den_org_area = density(org_face_feature$Area_CF)
    den_org_area = data.frame(x=den_org_area$x, y=den_org_area$y)
    
    den_sim_area = density(face_features$area)
    den_sim_area = data.frame(x=den_sim_area$x, y=den_sim_area$y)
    
    print(ggplot(data = den_org_area, aes(y)) +
        geom_density(fill = "blue", alpha = 0.3) +
        geom_density(data = den_sim_area, aes(y), fill="red", alpha =0.3) +
        
        theme(legend.position="top", legend.text=element_text(size=16), legend.title = element_blank(),
              legend.box.margin=margin(-10,-10,-10,-10),
              plot.title = element_text(hjust = 0.5, size=18),
              plot.subtitle = element_text(hjust = 0.5, size=16),
              axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
              axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
              panel.background = element_rect(fill='white', colour='black'),
              panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
        xlab(expression(paste("Face Area"))) + ylab("Density")+
        labs(title = "Comparison of face feature of original and simulated networks") )   # the titles needs changing for different runs
    
    # elongation
    den_org_elong = density(org_face_feature$Elong.)
    den_org_elong = data.frame(x=den_org_elong$x, y=den_org_elong$y)
    
    den_sim_elong = density(face_features$elong)
    den_sim_elong = data.frame(x=den_sim_elong$x, y= den_sim_elong$y)
    
    print(ggplot(data = den_org_elong, aes(y)) +
        geom_density(fill = "blue", alpha = 0.3) +
        geom_density(data = den_sim_elong, aes(y), fill="red", alpha =0.3) +
        
        theme(legend.position="top", legend.text=element_text(size=16), legend.title = element_blank(),
              legend.box.margin=margin(-10,-10,-10,-10),
              plot.title = element_text(hjust = 0.5, size=18),
              plot.subtitle = element_text(hjust = 0.5, size=16),
              axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
              axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
              panel.background = element_rect(fill='white', colour='black'),
              panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
        xlab(expression(paste("Face Elongation"))) + ylab("Density")+
        labs(title = "Comparison of face feature of original and simulated networks")  )  # the titles needs changing for different runs
    
    # orientation
    den_org_orient = density(org_face_feature$Orient.)
    den_org_orient = data.frame(x=den_org_orient$x, y=den_org_orient$y)
    
    den_sim_orient = density(face_features$orient)
    den_sim_orient = data.frame(x=den_sim_orient$x, y= den_sim_orient$y)
    
    print(ggplot(data = den_org_orient, aes(y)) +
        geom_density(fill = "blue", alpha = 0.3) +
        geom_density(data = den_sim_orient, aes(y), fill="red", alpha =0.3) +
        
        theme(legend.position="top", legend.text=element_text(size=16), legend.title = element_blank(),
              legend.box.margin=margin(-10,-10,-10,-10),
              plot.title = element_text(hjust = 0.5, size=18),
              plot.subtitle = element_text(hjust = 0.5, size=16),
              axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
              axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
              panel.background = element_rect(fill='white', colour='black'),
              panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
        xlab(expression(paste("Face Orientation"))) + ylab("Density")+
        labs(title = "Comparison of face feature of original and simulated networks")  )  # the titles needs changing for different runs
}


########################################################################
#### computes an edge weight based on the degree of its two end vertices
#### g2_degree is the list of degree of each vertex in the graph under consideration
#### x is the edge whose weight is being computed, the 5th and 6th index provides the indices of the end vertices
computeEdgeWeight <- function(g2_degree, x){
    return(g2_degree[x[5]] + g2_degree[x[6]])
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


deterministicEdges_2 <- function(branch.ppp, branch.all, org_face_feature, sample_id){
    #### construct the Delaunay triangulation on the parent points as a starter network
    #### ord_point_list: to maintain the order of the points
    set.seed(Sys.time())
    ord_point_list = data.frame(x = branch.ppp$x, y = branch.ppp$y)
    
    network_triangulation = deldir(ord_point_list[, 1:2])
    
    #### construct a convenient data structure to keep the triangulation information
    network_extra1 = data.frame(x1=network_triangulation$delsgs$x1, y1=network_triangulation$delsgs$y1,
                                x2=network_triangulation$delsgs$x2, y2=network_triangulation$delsgs$y2,
                                ind1=network_triangulation$delsgs$ind1, ind2=network_triangulation$delsgs$ind2)
    
    network_extra1$euclidDist = apply(network_extra1, 1, function(x) sqrt( ((x[1]-x[3])^2) + ((x[2]-x[4])^2) ) ) 
    network_extra1$anglecomp = apply(network_extra1, 1, function(x) calcAngle(x))
    
    #### new
    #### reminder: The initial triangulation already has closed faces at  the boundary, no need to add extra edges before computing the faces.
    #### compute the face area of the triangulation
    graph_obj =  graph_from_data_frame(unique(network_extra1[, 5:6]), directed = FALSE) 
    
    #### Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
    #### This is sometimes also called the clustering coefficient.
    cluster_coeff_t = igraph::transitivity(graph_obj, type = "global")
    cat("CC DT: ", cluster_coeff_t, "\n")
    
    g_o <- as_graphnel(graph_obj) ## Convert igraph object to graphNEL object for planarity testing
    boyerMyrvoldPlanarityTest(g_o)
    
    face_list = planarFaceTraversal(g_o)
    face_node_count = sapply(face_list, length)
    
    #### applying the shoe lace formula
    #### the original point pattern is being unmarked 
    #### as we only need the coordinates for the shoelace formula
    u_branch.ppp = unmark(branch.ppp)
    face_area_list = sapply(face_list, function(x) faceArea(x, u_branch.ppp))
    
    #### eliminating the outer face, it has the largest face area
    face_node_count = face_node_count[-which.max(face_area_list)]
    face_list = face_list[-which.max(face_area_list)]
    face_area_list = face_area_list[-which.max(face_area_list)]
    
    #### face features computation
    columns = c("Area_CF", "Perim.", "Ext.", "Disp.", "Elong.", "Eccentr.", "Orient.") # Area_CF: from contour function
    face_features = data.frame(matrix(nrow = 0, ncol = length(columns)))
    colnames(face_features) = columns
    
    for(f in c(1: length(face_list))){
        f_feat = computeFacefeatures(f, face_list, u_branch.ppp, NULL)
        face_features = rbind(face_features, f_feat)
        
    }# loop ends for each face of the triangulation
    #### new end
    
    #### create another copy of the data structure
    network_extra = rbind(network_extra1)
    
    #### construct and display as corresponding ppp and linnet
    g_o_lin = linnet(branch.ppp, edges=as.matrix(network_extra[,5:6]))
    
    plot(g_o_lin, main="Initial DT", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", "dodgerblue", "white", "maroon1",
                                                                    "mediumpurple"))
    
    ####new
    ####plot the feature densities of the initial triangulation vs the original network for comparison
    comparePlotOrgSim(org_face_feature, face_features)
    ####new end
    
    #### degree list of the newly constructed Delaunay graph
    #### create a graph from Delaunay triangulation
    g_o_degree = igraph::degree(graph_obj)
    
    network_extra$weight = apply(network_extra, 1, function(x) computeEdgeWeight(g_o_degree, x))
    network_extra$weight = range01(network_extra$weight)
    
    triKDE_face_area = kde(as.matrix(face_area_list))
    
    return(list(network_extra, face_list, face_area_list, face_node_count, triKDE_face_area, g_o_degree))
    
}


rejectionSampling_2 <- function(branch.ppp, network_extra, face_list, face_area_list, face_node_count, 
                                g2_degree, orgKDE_face_area, triKDE_face_area, 
                                meshedness, network_density, compactness, sample_id){
    #### some initialization
    rejected1 = 0
    noChange = 0
    
    set.seed(Sys.time())
    
    network_extra$accepted = 0
    dist_to_boundary = bdist.points(branch.ppp)
    
    while (TRUE) {
        cat("noChange: ", noChange, "\n")  
        if(noChange == 200){    # if the network has not been changed for 200 iterations
            break
        }
        
        N = branch.ppp$n
        E = length(network_extra$x1)
        A = summary(branch.ppp)$window$area
        L = sum(network_extra$euclidDist)
        
        mesh = (E-N+1)/((2*N)-5)
        n_density = E/((3*N)-6)
        compact = 1- ((4*A)/(L-(2*sqrt(A)))^2)
        
        cat(meshedness, " ", network_density, " ", compactness, "\n")
        cat(mesh, " ", n_density, " ", compact, "\n\n")
        
        if(mesh <= meshedness || n_density <= network_density || compact <= compactness){   # if network is getting too sparse
            break
        }
        
        #### sample an edge index based on the degree-based weights assigned to the edges
        #i = sample.int(length(network_extra[, 1]), 1, prob = network_extra$weight)
        i = sample.int(length(network_extra[, 1]), 1)
        
        # if(network_extra$accepted[i] == 1){     # if the samples edge has already been accepted skip to the next sampling
        #     noChange = noChange + 1
        #     next
        # }
        
        cat("index: ", i, " current num of edges: ", length(network_extra[, 1]), "\n")
        
        #### check the degree of the nodes of the selected edge
        #### if removing the edge lowers the connectivity significantly or disconnects the network - skip
        fromDeg = g2_degree[network_extra[i, ]$ind1]
        toDeg = g2_degree[network_extra[i, ]$ind2]
        
        cat("from deg: ", fromDeg, " to degree: ", toDeg, "\n")
        
        # if(fromDeg <= 3 || toDeg <= 3){
        #     if(dist_to_boundary[network_extra[i, ]$ind1] > 100 && dist_to_boundary[network_extra[i, ]$ind2] > 100){     # we allow the degree of the nodes close to the boundary to be lower than the rest
        #         cat("Edge not removed [degree constraint]...\n\n")
        #         network_extra$accepted[i] = 1
        #         noChange = noChange + 1
        #         next
        #     }
        # }
        
        #### finding out the faces the chosen edge is a part of
        face_index = which(unlist(lapply(face_list, function(x) isEdgeOnFace(x, c(network_extra[i, ]$ind1, network_extra[i, ]$ind2)))))
        
        accept = FALSE
        for (f in face_index) {
            org_est = predict(orgKDE_face_area, x=face_area_list[f])
            tri_est = predict(triKDE_face_area, x=face_area_list[f])
            
            if(tri_est <= org_est){
                accept = accept | TRUE 
            }else{
                accept = accept | FALSE
            }
        }
        
        #### accept 
        if(accept){
            cat("edge accepted...\n")
            network_extra$accepted[i] = 1
        }
        else{
            #### check for network connectivity
            network_temp = network_extra[-c(i), ]
            g2 = make_empty_graph() %>% add_vertices(branch.ppp$n)
            g2 = add_edges(as.undirected(g2), as.vector(t(as.matrix(network_temp[,5:6]))))
            
            if(is_connected(g2)){
                #### adjust the degree of the nodes
                g2_degree[network_extra[i, ]$ind1] = g2_degree[network_extra[i, ]$ind1] - 1
                g2_degree[network_extra[i, ]$ind2] = g2_degree[network_extra[i, ]$ind2] - 1
                
                network_extra = network_temp
                
                cat("edge rejected...\n")
                rejected1 = rejected1 + 1
                noChange = 0
                
                #### recompute the degree based edge weights
                network_extra$weight = apply(network_extra, 1, function(x) computeEdgeWeight(g2_degree, x))
                network_extra$weight = range01(network_extra$weight)
                
                #### recompute the density estimation of the triangulation
                graph_obj =  graph_from_data_frame(unique(network_extra[, 5:6]), directed = FALSE) 
                
                g_o <- as_graphnel(graph_obj) ## Convert igraph object to graphNEL object for planarity testing
                boyerMyrvoldPlanarityTest(g_o)
                
                face_list = planarFaceTraversal(g_o)
                face_node_count = sapply(face_list, length)
                
                #### applying the shoe lace formula
                #### the original point pattern is being unmarked 
                #### as we only need the coordinates for the shoelace formula
                u_branch.ppp = unmark(branch.ppp)
                face_area_list = sapply(face_list, function(x) faceArea(x, u_branch.ppp))
                
                #### eliminating the outer face, it has the largest face area
                face_node_count = face_node_count[-which.max(face_area_list)]
                face_list = face_list[-which.max(face_area_list)]
                face_area_list = face_area_list[-which.max(face_area_list)]
                
                triKDE_face_area = kde(as.matrix(face_area_list))
            }
            else{
                cat("Edge not removed [connectivity constraint]...\n\n")
                network_extra$accepted[i] = 1
            }
        }
        noChange = noChange + 1
        cat("\n")
        
        # ####
        # g2_lin = linnet(branch.ppp, edges=as.matrix(network_extra[, 5:6]))
        # 
        # plot(branch.ppp, cex=1, pch=20, main="", bg=1)
        # plot(g2_lin, add=T)
        # Sys.sleep(2)
        # ####
    }
    
    return(network_extra)
    
}


generateNetworkEdges_2 <- function(branch.ppp, branch_all, org_face_feature, orgKDE_face_area,
                       meshedness, network_density, compactness,
                       sample_id){
    
    #### constructing the deterministic Delaunay triangulation as the initial ganglionic network
    triangulation_info_list = deterministicEdges_2(branch.ppp, branch.all, org_face_feature, sample_id)
    
    #### returned values
    network_extra1 = triangulation_info_list[[1]]
    face_list = triangulation_info_list[[2]]
    face_area_list = triangulation_info_list[[3]]
    face_node_count = triangulation_info_list[[4]]
    triKDE_face_area = triangulation_info_list[[5]]
    g2_degree = triangulation_info_list[[6]]
    
    #### remove edges from the initial triangulation by rejection sampling
    network_extra = rejectionSampling_2(branch.ppp, network_extra1, face_list, face_area_list, face_node_count, 
                                        g2_degree, orgKDE_face_area, triKDE_face_area, 
                                        meshedness, network_density, compactness, sample_id)
    
    #### observe if the distribution under consideration covers the targeted distribution
    # par(mar=c(5, 2, 1, 1))
    # plot(density(branch.all$angle), col="red", lty=2, xlab="Edge angle", ylab="Density", main="", xlim=c(-100,250), ylim=c(0, 0.02))
    # lines(density(network_extra1$anglecomp), col="blue")
    # lines(density(network_extra1$anglecomp), col="black", lty=3)
    # legend(x=0, y=0.02, legend=c("ENS angle", "Triangulation angle", "Simulated angle"), 
    #        col=c("red", "blue", "black"), 
    #        lty=c(2, 1, 3))
    # 
    # plot(density(branch.all$euclid), col="red", lty=2, xlab="Edge Length", ylab="Density", main="")
    # lines(density(network_extra1$euclidDist), col="blue")
    # lines(density(network_extra1$euclidDist), col="black", lty=3)
    # legend(x=300, y=0.006, legend=c("ENS edge len", "Triangulation edge len", "Simulated edge len"), 
    #        col=c("red", "blue", "black"), 
    #        lty=c(2, 1, 3))
    # 
    # #### visualize the bivariate distribution of the trimmed triangulation
    # den3d = kde2d(network_extra$anglecomp, network_extra$euclidDist)
    # persp(den3d, theta = -45, phi = 30, xlab="angle", ylab="edge len",
    #       ticktype = "detailed", shade = 0.75, col="lightblue")
    
    #### create a graph from sampled Delaunay triangulation
    g2 = make_empty_graph() %>% add_vertices(branch.ppp$n)
    g2 = add_edges(as.undirected(g2), as.vector(t(as.matrix(network_extra[,5:6]))))
    
    #### display as corresponding ppp and linnet
    g2_lin = linnet(branch.ppp, edges=as.matrix(network_extra[, 5:6]))
    
    plot(branch.ppp, cex=1, pch=20, main="", bg=1)
    plot(g2_lin, add=T)
    
    
    return(list(network_extra, g2_lin))
}


#### extracting parent directory information for accessing input and output location
dir = this.dir()
folder = strsplit(dir, "/")
folder = folder[[1]][length(folder[[1]])]
parent = strsplit(dir, folder)

face_folder = paste(parent, "Outputs/ENSMouse/FaceFeature/", sep="")
face_features_combined = read.csv(paste(face_folder, "FaceFeatures_2.csv", sep = ""))


#### the TIF images of the ganglionic networks are preprocessed in Fiji (ImageJ) and the network information is extracted as .csv files
#### choose one/all of the ganglionic network samples with file chooser below
branch_info_folder = paste(parent, "Data/ENSMouse Branch Information (in um) v2.0/", sep="")
branch_info_files = list.files(branch_info_folder, recursive = TRUE, pattern = "\\.csv", full.names = TRUE)

i = 2

ens_location = strsplit(branch_info_files[i], "/")[[1]][11]
sample_id = strsplit(strsplit(branch_info_files[i], "/")[[1]][12], "\\.")[[1]][1]
cat("\n(", i, ") Location: ", ens_location, "\nSample Id: ", sample_id, "\n")

max_y = 1 # 4539.812 found by computation; right now keeping everything unscaled as the moments can not be computed otherwise

data_struct_list = constructDataStruct(sample_id, parent, branch_info_files[i], output_folder_path, max_y)

#### the returned values
branch.all = data_struct_list[[1]]
branch.ppp = data_struct_list[[2]]
branch.lpp = data_struct_list[[3]]
g1 = data_struct_list[[4]]
hardcoreStrauss_model_param = data_struct_list[[5]]

plot(branch.lpp, main="original", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", "dodgerblue", "white", "maroon1",
                                                       "mediumpurple"))
                                                       #save the network plot

#### alpha, gamma, psi
N = branch.ppp$n
E = length(branch.all$x1)
A = summary(branch.ppp)$window$area
L = sum(branch.all$euclid)

meshedness = (E-N+1)/((2*N)-5)
network_density = E/((3*N)-6)
compactness = 1- ((4*A)/(L-(2*sqrt(A)))^2)

#### Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
#### This is sometimes also called the clustering coefficient.
cluster_coeff = igraph::transitivity(g1, type = "global")
cat("CC original: ", cluster_coeff, "\n")

#### filter out the face face features of the sample under consideration
#### reminder: the face features were computed assuming additional edges were computed to close the open faces at the boundary
#### similar thing should be done for the new network too if needed.
face_feature = face_features_combined[face_features_combined$sample_id == sample_id, ]
plot(density(face_feature$Area_SL), main="Feature density in original network")

orgKDE_face_area = kde(as.matrix(face_feature$Area_SL))

network_info_list = generateNetworkEdges_2(branch.ppp, branch_all, face_feature, orgKDE_face_area,
                                         meshedness, network_density, compactness,
                                         sample_id)

net_data_struct = network_info_list[[1]]
linnet_obj = network_info_list[[2]]

branch.lpp_2 = lpp(branch.ppp, linnet_obj )
plot(branch.lpp_2, main="simulated", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", "dodgerblue", "white", "maroon1",
                                                       "mediumpurple"))
                                                       

#### compute the face area of the newly constructed network to compare the density with the original face area
#### Reminder: for the simulated network there are many open boundary faces
#### additional edges are required
pp_nodes = data.frame(x=branch.ppp$x, y=branch.ppp$y)
#### corner nodes
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

#### the new edges will connect the boundary nodes computed above sequentially and create a loop by connecting the last node to the first one.
#### to avoid shifting error
if(length(boundary_1) == 1){
    b_1 = c()
}else{
    b_1 = boundary_1[2:length(boundary_1)]
}
new_edges = data.frame(ind1 = c(boundary_1, boundary_2, boundary_3, boundary_4),
                       ind2 = c(b_1, boundary_2, boundary_3, boundary_4, boundary_1[1]))
new_edges = new_edges[new_edges$ind1 != new_edges$ind2, ]  # removing self loops at the corner nodes

graph_obj =  graph_from_data_frame(unique(rbind(net_data_struct[, 5:6], new_edges)), directed = FALSE) # new graph object that combines the actual network and the new edges
graph_obj = delete_edges(graph_obj, which(which_multiple(graph_obj)))

#### Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
#### This is sometimes also called the clustering coefficient.
cluster_coeff_s = igraph::transitivity(graph_obj, type = "global")
cat("CC simulated: ", cluster_coeff_s, "\n")

g_o <- as_graphnel(graph_obj) ## Convert igraph object to graphNEL object for planarity testing
boyerMyrvoldPlanarityTest(g_o)

face_list = planarFaceTraversal(g_o)
face_node_count = sapply(face_list, length)

#### applying the shoe lace formula
#### include the corner nodes to a temporary point pattern, the original point pattern is being unmarked 
#### as we only need the coordinates for the shoelace formula
corner.ppp = ppp(x=c(branch.ppp$window$xrange[1], branch.ppp$window$xrange[2], branch.ppp$window$xrange[2], branch.ppp$window$xrange[1]), 
                 y=c(branch.ppp$window$yrange[2], branch.ppp$window$yrange[2], branch.ppp$window$yrange[1], branch.ppp$window$yrange[1]),
                 window = branch.ppp$window)
u_branch.ppp = unmark(branch.ppp)
face_area_list = sapply(face_list, function(x) faceArea(x, superimpose.ppp(u_branch.ppp, corner.ppp)))

#### eliminating the outer face, it has the largest face area
face_node_count = face_node_count[-which.max(face_area_list)]
face_list = face_list[-which.max(face_area_list)]
face_area_list = face_area_list[-which.max(face_area_list)]

#### face features computation
columns = c("Area_CF", "Perim.", "Ext.", "Disp.", "Elong.", "Eccentr.", "Orient.") # Area_CF: from contour function
face_features_sim = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(face_features_sim) = columns

for(f in c(1: length(face_list))){
    f_feat = computeFacefeatures(f, face_list, u_branch.ppp, corner.ppp)
    face_features_sim = rbind(face_features_sim, f_feat)
    
}# loop ends for each face of the triangulation

comparePlotOrgSim(face_feature, face_features_sim)