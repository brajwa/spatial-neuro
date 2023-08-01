#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim", "RBGL", "svglite", "ggrepel", "devtools", "geosphere", "philentropy")

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
    #cat("face id: ", f, "\n")
    
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
    orient = 0.5 * atan2((2 * moments$value[moments$moment == "m11"]) , (moments$value[moments$moment == "m20"] - moments$value[moments$moment == "m02"]))
    orient = orient * 180 / pi
    eccentr = (((moments$value[moments$moment == "m02"] - moments$value[moments$moment == "m20"]) * (moments$value[moments$moment == "m02"] - moments$value[moments$moment == "m20"])) 
               + (4 * moments$value[moments$moment == "m11"] * moments$value[moments$moment == "m11"])) / moments$value[moments$moment == "m00"]
    
    #cat(area, perim, ext, disp, elong, eccentr, orient, "\n")
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


deterministicEdges_3 <- function(branch.ppp, branch.all, org_face_feature, sample_id){
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
    #### Reminder: The initial triangulation already has closed faces at  the boundary, 
    #### no need to add extra edges before computing the faces.
    #### compute the face feature of the triangulation
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
    degs = (igraph::degree(graph_obj, mode="total"))
    ord = order(as.numeric(names(degs)))
    degs = degs[ord]
    
    #### attach the degree information to the point pattern for proper visualization
    marks(branch.ppp) = factor(degs)
    branch.ppp$markformat = "factor"
    g_o_lin = linnet(branch.ppp, edges=as.matrix(network_extra[,5:6]))
    branch.lpp_dt = lpp(branch.ppp, g_o_lin )
    
    plot(branch.lpp_dt, main="Initial DT", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", 
                                                                        "dodgerblue", "white", "maroon1", 
                                                                        "mediumpurple"))
                                                                  
    ####new
    ####plot the feature densities of the initial triangulation vs the original network for comparison
    comparePlotOrgSim(org_face_feature, face_features)
    ####new end
    
    #### degree list of the constructed Delaunay graph again for edge weight calculation
    g_o_degree = degs
    
    network_extra$weight = apply(network_extra, 1, function(x) computeEdgeWeight(g_o_degree, x))
    network_extra$weight = range01(network_extra$weight)
    
    triKDE_face_feat = kde(as.matrix(data.frame(face_area_list, 
                                                face_features$elong, 
                                                face_features$orient)))
    
    return(list(network_extra, face_list, face_area_list, face_node_count, triKDE_face_feat, g_o_degree, face_features))
}


rejectionSampling_3(branch.ppp, network_extra1, face_list, face_area_list, face_node_count, 
                    g2_degree, orgKDE_face_feat, triKDE_face_feat, 
                    meshedness, network_density, compactness, cluster_coeff, sample_id, tri_face_features){
    
    vertex_dist_boundary = bdist.points(branch.ppp)
    
    noChange = 0
    while (TRUE) {
        if(noChange == 200){    # if the network has not been changed for 200 iterations
            break
        }
        
        #### cc of the current network
        cc_cur = ccFromDataframe(branch.ppp, network_extra1)
        if(cc_cur <= cluster_coeff){
            cat("Rejection sampling ended [CC reached target]\n")
            break
        }
        
        #### select a vertex at random or based on high degree
        set.seed(Sys.time())
        
        #### prepare the vertex probability from degree values
        prob_vertex = g2_degree / sum(g2_degree)
        selected_vertex = sample.int(branch.ppp$n, 1, prob = prob_vertex)
        
        cat("\nSelected vertex ID: ", selected_vertex, ", Degree of the selected vertex: ", g2_degree[selected_vertex], "\n")
        #points(branch.ppp$x[selected_vertex], branch.ppp$y[selected_vertex], col="red", cex=2, pch=19)
        
        #### detect the neighbors of a given vertex
        adj_vertices_from_df = c(network_extra1$ind1[network_extra1$ind2==selected_vertex],
                                 network_extra1$ind2[network_extra1$ind1==selected_vertex])
        
        #### list of the edges that are between those neighboring vertices (if any)
        #### the indices are of network_extra1
        edge_bet_adj_vertices = which((network_extra1$ind1 %in% adj_vertices_from_df) & 
                                          (network_extra1$ind2 %in% adj_vertices_from_df))
        
        #### for each bet adj edge, there are two adj edges
        for (e in edge_bet_adj_vertices) {
            e1 = which(((network_extra1$ind1==network_extra1$ind1[e]) & (network_extra1$ind2==selected_vertex))|
                                                        ((network_extra1$ind2==network_extra1$ind1[e]) & (network_extra1$ind1==selected_vertex)) )
            
            e2 = which(((network_extra1$ind1==network_extra1$ind2[e]) & (network_extra1$ind2==selected_vertex))|
                                                        ((network_extra1$ind2==network_extra1$ind2[e]) & (network_extra1$ind1==selected_vertex)) )
            cat(e, e1, e2, "\n")
            
            #### pick an edge random from these three ones
            selected_edge = sample(c(e, e1, e2), 1)
            cat("Selected edge ID: ", selected_edge, "\n")
            
            #### check for the degree of the end vertices
            #### degree-one vertices only at the boundary
            v1 = network_extra1$ind1[selected_edge]
            v2 = network_extra1$ind2[selected_edge]
            #lines(c(branch.ppp$x[v1], branch.ppp$x[v2]), c(branch.ppp$y[v1], branch.ppp$y[v2]), col="red", lwd=2.5)
            
            #### if any of the end vertices of the selected edge has degree 2, deleting it will create a dangling vertex
            #### and if that is not on the boundary we won't allow it
            if((g2_degree[v1]==2 | g2_degree[v2]==2) & (vertex_dist_boundary[v1]!=0 | vertex_dist_boundary[v2]!=0)){
                #noChange = noChange + 1
                cat("Edge kept [Boundary degree constraint]\n")
                next
            }
            
            #### if v1 and v2 both are boundary vertices, eliminating the edge will increase the outer face
            #### but no other comparison will be possible for a new face
            #### so we reject the edge without any more comparison
            if((vertex_dist_boundary[v1]==0 & vertex_dist_boundary[v2]==0)){
                temp_network_extra1 = network_extra1[-c(selected_edge), ]
                temp_graph_obj = make_empty_graph() %>% add_vertices(branch.ppp$n)
                temp_graph_obj = add_edges(as.undirected(temp_graph_obj), 
                                           as.vector(t(as.matrix(temp_network_extra1[,5:6]))))
                temp_g_o <- as_graphnel(temp_graph_obj) 
                boyerMyrvoldPlanarityTest(temp_g_o)
                
                temp_face_list = planarFaceTraversal(temp_g_o)
                temp_face_node_count = sapply(temp_face_list, length)
                
                #### applying the shoe lace formula
                u_branch.ppp = unmark(branch.ppp)
                temp_face_area_list = sapply(temp_face_list, function(x) faceArea(x, u_branch.ppp))
                
                #### eliminating the outer face, it has the largest face area
                temp_face_node_count = temp_face_node_count[-which.max(temp_face_area_list)]
                temp_face_list = temp_face_list[-which.max(temp_face_area_list)]
                temp_face_area_list = temp_face_area_list[-which.max(temp_face_area_list)]
                
                #### face features computation
                temp_columns = c("Area_CF", "Perim.", "Ext.", "Disp.", "Elong.", "Eccentr.", "Orient.") # Area_CF: from contour function
                temp_face_features = data.frame(matrix(nrow = 0, ncol = length(temp_columns)))
                colnames(temp_face_features) = temp_columns
                
                for(f in c(1: length(temp_face_list))){
                    temp_f_feat = computeFacefeatures(f, temp_face_list, u_branch.ppp, NULL)
                    temp_face_features = rbind(temp_face_features, temp_f_feat)
                    
                }# loop ends for each face of the temp network
                
                temp_triKDE_face_feat = kde(as.matrix(data.frame(temp_face_area_list, 
                                                                 temp_face_features$elong, 
                                                                 temp_face_features$orient)))
                
                #### remove the edge, there is a change
                noChange = 0
                cat("Boundary edge deleted\n")
                
                #### make necessary changes permanent
                network_extra1 = temp_network_extra1
                g2_degree = igraph::degree(temp_graph_obj, mode="total")
                
                face_list = temp_face_list
                face_area_list = temp_face_area_list
                face_node_count = temp_face_node_count
                triKDE_face_feat = temp_triKDE_face_feat
                tri_face_features = temp_face_features
                
                next
            }
            
            #### check if removing that edge disconnects the network
            temp_network_extra1 = network_extra1[-c(selected_edge), ]
            #rownames(temp_network_extra1) = NULL # to reset the row index after row deletion
            
            temp_graph_obj = make_empty_graph() %>% add_vertices(branch.ppp$n)
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
                u_branch.ppp = unmark(branch.ppp)
                temp_face_area_list = sapply(temp_face_list, function(x) faceArea(x, u_branch.ppp))
                
                #### eliminating the outer face, it has the largest face area
                temp_face_node_count = temp_face_node_count[-which.max(temp_face_area_list)]
                temp_face_list = temp_face_list[-which.max(temp_face_area_list)]
                temp_face_area_list = temp_face_area_list[-which.max(temp_face_area_list)]
                
                #### face features computation
                temp_columns = c("Area_CF", "Perim.", "Ext.", "Disp.", "Elong.", "Eccentr.", "Orient.") # Area_CF: from contour function
                temp_face_features = data.frame(matrix(nrow = 0, ncol = length(temp_columns)))
                colnames(temp_face_features) = temp_columns
                
                for(f in c(1: length(temp_face_list))){
                    temp_f_feat = computeFacefeatures(f, temp_face_list, u_branch.ppp, NULL)
                    temp_face_features = rbind(temp_face_features, temp_f_feat)
                    
                }# loop ends for each face of the temp network
                
                temp_triKDE_face_feat = kde(as.matrix(data.frame(temp_face_area_list, 
                                                            temp_face_features$elong, 
                                                            temp_face_features$orient)))
                
                ####finding the new face
                face_p = c()
                for(f_i in face_index){
                    face_p = c(face_p, as.numeric(unlist(face_list[f_i])))
                }
                face_p = unique(face_p)
                face_p_index = which(sapply(lapply(temp_face_list, function(x) sort(as.numeric(unlist(x)))), 
                                            identical, sort(face_p)))
                
                edge_reject = FALSE
                
                #### prediction
                org_est = predict(orgKDE_face_feat, x=c(temp_face_area_list[face_p_index], 
                                                        temp_face_features$elong[face_p_index], 
                                                        temp_face_features$orient[face_p_index]))
                temp_tri_est = predict(temp_triKDE_face_feat, x=c(temp_face_area_list[face_p_index], 
                                                                  temp_face_features$elong[face_p_index], 
                                                                  temp_face_features$orient[face_p_index]))
                
                if(length(face_index)==2){
                    f1 = face_index[1]
                    f2 = face_index[2]
                    
                    org_est1 = predict(orgKDE_face_feat, x=c(face_area_list[f1], tri_face_features$elong[f1], tri_face_features$orient[f1]))
                    tri_est1 = predict(triKDE_face_feat, x=c(face_area_list[f1], tri_face_features$elong[f1], tri_face_features$orient[f1]))
                    
                    org_est2 = predict(orgKDE_face_feat, x=c(face_area_list[f2], tri_face_features$elong[f2], tri_face_features$orient[f2]))
                    tri_est2 = predict(triKDE_face_feat, x=c(face_area_list[f2], tri_face_features$elong[f2], tri_face_features$orient[f2]))
                    
                    if((org_est1 < tri_est1) & (org_est2 < tri_est2) & (org_est >= temp_tri_est)){
                        edge_reject = TRUE
                    }
                    
                }else if(length(face_index)==1){
                    f1 = face_index[1]

                    org_est1 = predict(orgKDE_face_feat, x=c(face_area_list[f1], tri_face_features$elong[f1], tri_face_features$orient[f1]))
                    tri_est1 = predict(triKDE_face_feat, x=c(face_area_list[f1], tri_face_features$elong[f1], tri_face_features$orient[f1]))
                    
                    if((org_est1 < tri_est1) & (org_est >= temp_tri_est)){
                        edge_reject = TRUE
                    }
                    
                }else{
                    cat("Error in face identification\n")
                    quit()
                }
                
                if(edge_reject){
                    #### remove the edge, there is a change
                    noChange = 0
                    cat("Edge deleted\n")
                    
                    #### make necessary changes permanent
                    network_extra1 = temp_network_extra1
                    g2_degree = igraph::degree(temp_graph_obj, mode="total")
                    
                    face_list = temp_face_list
                    face_area_list = temp_face_area_list
                    face_node_count = temp_face_node_count
                    triKDE_face_feat = temp_triKDE_face_feat
                    tri_face_features = temp_face_features
                    
                }else{
                    #### keep the edge, no change
                    noChange = noChange + 1
                    cat("Edge kept [Face feature estimation constraint]\n")
                }
                
            }else{
                #### keep the edge, no change
                noChange = noChange + 1
                cat("Edge kept [Connectivity constraint]\n")
            }
        }
        
        #### temporarily plotting each iteration
        graph_obj =  make_empty_graph() %>% add_vertices(branch.ppp$n)
        graph_obj = add_edges(as.undirected(graph_obj), 
                              as.vector(t(as.matrix(temp_network_extra1[,5:6]))))
        
        #### Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
        #### This is sometimes also called the clustering coefficient.
        cluster_coeff_s = igraph::transitivity(graph_obj, type = "global")
        cat("CC Sim: ", cluster_coeff_s, "\n")
        
        #### construct and display as corresponding ppp and linnet
        degs = igraph::degree(graph_obj, mode="total")
        # ord = order(as.numeric(names(degs)))
        # degs = degs[ord]
        
        #### attach the degree information to the point pattern for proper visualization
        marks(branch.ppp) = factor(degs)
        branch.ppp$markformat = "factor"
        g_o_lin = linnet(branch.ppp, edges=as.matrix(network_extra1[,5:6]))
        branch.lpp_s = lpp(branch.ppp, g_o_lin )
        
        plot(branch.lpp_s, main="Sim", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", 
                                                                    "dodgerblue", "white", "maroon1", 
                                                                    "mediumpurple", "yellow", "cyan"))
    }
    
    ####
    graph_obj =  make_empty_graph() %>% add_vertices(branch.ppp$n)
    graph_obj = add_edges(as.undirected(graph_obj), 
              as.vector(t(as.matrix(temp_network_extra1[,5:6]))))
    
    #### Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
    #### This is sometimes also called the clustering coefficient.
    cluster_coeff_s = igraph::transitivity(graph_obj, type = "global")
    cat("CC Sim: ", cluster_coeff_s, "\n")
    
    #### construct and display as corresponding ppp and linnet
    degs = igraph::degree(graph_obj, mode="total")
    # ord = order(as.numeric(names(degs)))
    # degs = degs[ord]
    
    #### attach the degree information to the point pattern for proper visualization
    marks(branch.ppp) = factor(degs)
    branch.ppp$markformat = "factor"
    g_o_lin = linnet(branch.ppp, edges=as.matrix(network_extra1[,5:6]))
    branch.lpp_s = lpp(branch.ppp, g_o_lin )
    
    plot(branch.lpp_s, main="Sim", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", 
                                                                        "dodgerblue", "white", "maroon1", 
                                                                        "mediumpurple", "yellow", "cyan"))
    return(network_extra1)                                        
}


generateNetworkEdges_3 <- function(branch.ppp, branch_all, org_face_feature, orgKDE_face_feat,
                                   meshedness, network_density, compactness, cluster_coeff,
                                   sample_id){
    
    #### constructing the deterministic Delaunay triangulation as the initial ganglionic network
    triangulation_info_list = deterministicEdges_3(branch.ppp, branch.all, org_face_feature, sample_id)
    
    #### returned values
    network_extra1 = triangulation_info_list[[1]]
    
    face_list = triangulation_info_list[[2]]
    face_area_list = triangulation_info_list[[3]]
    face_node_count = triangulation_info_list[[4]]
    
    triKDE_face_feat = triangulation_info_list[[5]]
    
    g2_degree = triangulation_info_list[[6]]
    tri_face_features = triangulation_info_list[[7]]
    
    #### remove edges from the initial triangulation by rejection sampling
    network_extra = rejectionSampling_3(branch.ppp, network_extra1, face_list, face_area_list, face_node_count, 
                                        g2_degree, orgKDE_face_feat, triKDE_face_feat, 
                                        meshedness, network_density, compactness, cluster_coeff, sample_id, tri_face_features)
    
    #### create a graph from sampled triangulation
    g2 = make_empty_graph() %>% add_vertices(branch.ppp$n)
    g2 = add_edges(as.undirected(g2), as.vector(t(as.matrix(network_extra[,5:6]))))
    
    #### display as corresponding ppp and linnet
    g2_lin = linnet(branch.ppp, edges=as.matrix(network_extra[, 5:6]))
    
    return(list(network_extra, g2_lin))
}


####main 
#### extracting parent directory information for accessing input and output location
dir = this.dir()
folder = strsplit(dir, "/")
folder = folder[[1]][length(folder[[1]])]
parent = strsplit(dir, folder)

face_folder = paste(parent, "Outputs/ENSMouse/FaceFeature/", sep="")
face_features_combined = read.csv(paste(face_folder, "FaceFeatures_2.csv", sep = ""))


#### the TIF images of the ganglionic networks are preprocessed in Fiji (ImageJ) and 
#### the network information is extracted as .csv files
branch_info_folder = paste(parent, "Data/ENSMouse Branch Information (in um) v2.0/", sep="")
branch_info_files = list.files(branch_info_folder, recursive = TRUE, pattern = "\\.csv", full.names = TRUE)

i = 2  # index of the ENS network we want to work on

ens_location = strsplit(branch_info_files[i], "/")[[1]][11]
sample_id = strsplit(strsplit(branch_info_files[i], "/")[[1]][12], "\\.")[[1]][1]
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

#### filter out the face face features of the sample under consideration
#### Reminder: the face features were computed assuming additional edges were computed 
#### to close the open faces at the boundary
#### similar thing should be done for the new network too (where needed).
face_feature = face_features_combined[face_features_combined$sample_id == sample_id, ]

orgKDE_face_feat = kde(as.matrix(data.frame(face_feature$Area_SL, 
                                            face_feature$Elong., 
                                            face_feature$Orient.))) # this will be also tested by replacing by PC1 of all the face features we want

####At this point, new point pattern will be simulated. 
####For now we are generating networks on the original point pattern.
####branch.ppp will be replaced by some new ppp object

#### call the network generation functions
network_info_list = generateNetworkEdges_3(branch.ppp, branch_all, face_feature, orgKDE_face_feat,
                                           meshedness, network_density, compactness, cluster_coeff,
                                           sample_id)

#### returned values
net_data_struct = network_info_list[[1]]
linnet_obj = network_info_list[[2]]

#### creating a list representing degree for every node an use it in spatstat pattern
graph_obj_0 = graph_from_data_frame(net_data_struct[, 5:6], directed = FALSE)
degs = (igraph::degree(graph_obj_0, mode="total"))
ord = order(as.numeric(names(degs)))
degs = degs[ord]

#### attach the degree information to the point pattern for proper visualization
marks(branch.ppp) = factor(degs)
branch.ppp$markformat = "factor"
branch.lpp_2 = lpp(branch.ppp, linnet_obj )

plot(branch.lpp_2, main="simulated", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", "dodgerblue", 
                                                                  "white", "maroon1", "mediumpurple"))
                                                                  

#### compute the face features of the newly constructed network to compare the density with the original face feature
#### Reminder: for the simulated network there are many open boundary faces
#### additional edges are required
new_edges = computeBoundaryEdges(branch.ppp)

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
    
}# loop ends for each face of the simulated network

comparePlotOrgSim(face_feature, face_features_sim)