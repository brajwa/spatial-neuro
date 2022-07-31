########################################################################
#### computes an edge weight based on the degree of its two end vertices
#### g2_degree is the list of degree of each vertex in the graph under consideration
#### x is the edge whose weight is being computed, the 5th and 6th index provides the indices of the end vertices
computeEdgeWeight <- function(g2_degree, x){
  return(g2_degree[x[5]] + g2_degree[x[6]])
}


##############################################################################
#### predicts the pdf of angle and/or length for original or synthetic network
#### parameter: kde is the given density estimation
#### parameter: data is the data point for which prediction to be done
#### parameter: mode is an integer; either 1/2/3; 1 for angle, 2 for length; 3 for both
predictKDE <- function(kde, data, mode){
  if(mode == 1){
    return(predict(kde, x=data$anglecomp))      # edge angle
  }
  if(mode == 2){
    return(predict(kde, x=data$euclidDist))     # edge length
  }
  if(mode == 3){
    return(predict(kde, x=c(data$anglecomp, data$euclidDist)))      # both angle and length
  }
  if(mode == 4){
    return(predict(kde, x=data$weight))     # degree-based weight
  }
  if(mode == 5){
    return(predict(kde, x=c(data$anglecomp, data$euclidDist, data$weight)))     # all
  }
}


#######################################################################################
deterministicEdges <- function(ganglia_ppp, branch.all, sample_id, output_folder_path){
  
    #### construct the Delaunay triangulation on the parent points (simulated point pattern) as a starter network
    #### temp0: to maintain the order of the points
    set.seed(Sys.time())
    temp0 = data.frame(x = ganglia_ppp$x, y = ganglia_ppp$y)
    
    network_triangulation = deldir(temp0[, 1:2])
    
    #### construct a convenient data structure to keep the triangulation information
    network_extra1 = data.frame(x1=network_triangulation$delsgs$x1, y1=network_triangulation$delsgs$y1,
                               x2=network_triangulation$delsgs$x2, y2=network_triangulation$delsgs$y2,
                               ind1=network_triangulation$delsgs$ind1, ind2=network_triangulation$delsgs$ind2)
    
    network_extra1$euclidDist = apply(network_extra1, 1, function(x) sqrt( ((x[1]-x[3])^2) + ((x[2]-x[4])^2) ) ) 
    network_extra1$anglecomp = apply(network_extra1, 1, function(x) calcAngle(x))
    
    #### filter too long edges, comparing with the original network's maximum branch length
    network_extra1 = network_extra1[network_extra1$euclidDist <= max(branch.all$euclid), ]
    
    #### create another copy of the data structure
    network_extra = rbind(network_extra1)
    
    #### create a graph from Delaunay triangulation
    g2 = make_empty_graph() %>% add_vertices(ganglia_ppp$n)
    g2 = add_edges(as.undirected(g2), as.vector(t(as.matrix(network_extra[,5:6]))))
    
    #### construct and display as corresponding ppp and linnet
    g2_lin = linnet(ganglia_ppp, edges=as.matrix(network_extra[,5:6]))
    
    par(mar=c(0,0,0,0), oma=c(0,0,0,0))
    plot(ganglia_ppp, cex=1, pch=20, main="", bg=1)
    plot(g2_lin, add=T)
    
    #### degree list of the newl contructed Delaunay graph
    g2_degree = degree(g2)
    
    network_extra$weight = apply(network_extra, 1, function(x) computeEdgeWeight(g2_degree, x))
    network_extra$weight = range01(network_extra$weight)
    
    N = ganglia_ppp$n
    E = length(network_extra$x1)
    A = (ganglia_ppp$window$xrange[2]-ganglia_ppp$window$xrange[1])*(ganglia_ppp$window$yrange[2]-ganglia_ppp$window$yrange[1])
    L = sum(network_extra$euclidDist)
    
    #### compute the bivariate density distribution of the initial triangulation; 
    #### triKDE holds the kernel density estimation;
    triData = data.frame(angle=network_extra$anglecomp, len=network_extra$euclidDist, weight=network_extra$weight)
    
    triKDE_angle = computeKDE(triData, 1) # 1: angle, 2: length, 3: both
    triKDE_length = computeKDE(triData, 2) # 1: angle, 2: length, 3: both
    triKDE_both = computeKDE(triData, 3) # 1: angle, 2: length, 3: both
    
    #### visualize the bivariate distribution of angle and length of the initial triangulation
    den3d = kde2d(network_extra$anglecomp, network_extra$euclidDist)
    persp(den3d, theta = -45, phi = 30, xlab="angle", ylab="edge len",
        ticktype = "detailed", shade = 0.75, col="lightblue")
    
    #### observe if the distribution under consideration covers the targeted distribution
    par(mar=c(5, 2, 1, 1))
    plot(density(branch.all$angle), col="red", lty=2, xlab="Edge angle", ylab="Density", main="", xlim=c(-100,250), ylim=c(0, 0.02))
    lines(density(network_extra1$anglecomp), col="blue")
    legend(x=0, y=0.02, legend=c("ENS angle", "Triangulation angle"), 
         col=c("red", "blue"), 
         lty=c(2, 1))
    
    plot(density(branch.all$euclid), col="red", lty=2, xlab="Edge Length", ylab="Density", main="")
    lines(density(network_extra1$euclidDist), col="blue")
    legend(x=300, y=0.006, legend=c("ENS edge len", "Triangulation edge len"), 
         col=c("red", "blue"), 
         lty=c(2, 1))
    
    #### creating .pptx file to store the plots
    doc = read_pptx()
    
    #### another graph for the Delaunay triangulation for better visualization
    marked_g2  = graph_from_data_frame(network_extra1[,5:6], directed=FALSE)
    degs = degree(marked_g2)
    ord = order(as.numeric(names(degs)))
    degs = degs[ord]
    marked_ganglia_ppp = ganglia_ppp
    marks(marked_ganglia_ppp) = factor(degs)
    marked_ganglia_ppp$markformat = "factor"
    
    marked_g2_lin = linnet(marked_ganglia_ppp, edges=as.matrix(network_extra1[,5:6]))
    
    marked_branch.lpp = lpp(marked_ganglia_ppp, marked_g2_lin )
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(code = plot(marked_branch.lpp, main="triangulation", cex=1, pch=21, bg=c(3,4,5,6,7,8,"white","forestgreen","tan3"))), location = ph_location_fullsize())
    
    #### creating a directory to save the simulated network related files
    sim_net_path = paste(output_folder_path, "Simulated Network/", sep="")
    if (!dir.exists(sim_net_path)) {dir.create(sim_net_path, recursive=TRUE)}
    
    print(doc, target = paste(sim_net_path, sample_id, "_Delaunay_Tri_Simulated.pptx", sep=""))
    
    return(list(network_extra, triKDE_angle, triKDE_length, triKDE_both, g2_degree))

}


##################################################################################################################
rejectionSampling <- function(ganglia_ppp, network_extra, g2_degree, orgKDE_angle, orgKDE_length, orgKDE_both,
                              triKDE_angle, triKDE_length, triKDE_both, meshedness, network_density, compactness,
                              sample_id, output_folder_path){
  
  #### some initialization
  rejected1 = 0
  noChange = 0
  
  set.seed(Sys.time())
  
  network_extra$accepted = 0
  dist_to_boundary = bdist.points(ganglia_ppp)
  
  while (TRUE) {
    if(noChange == 200){    # if the network has not been changed for 200 iterations
      break
    }
    
    N = ganglia_ppp$n
    E = length(network_extra$x1)
    A = (ganglia_ppp$window$xrange[2]-ganglia_ppp$window$xrange[1])*(ganglia_ppp$window$yrange[2]-ganglia_ppp$window$yrange[1])
    L = sum(network_extra$euclidDist)
    
    mesh = (E-N+1)/((2*N)-5)
    n_density = E/((3*N)-6)
    compact = 1- ((4*A)/(L-(2*sqrt(A)))^2)
    cat(mesh, " ", n_density, " ", compact, "\n")
    
    if(mesh <= meshedness || n_density <= network_density || compact <= compactness){   # if network is getting too sparse
      break
    }
    
    #### sample an edge index based on the degree-based weights assigned to the edges
    i = sample.int(length(network_extra[, 1]), 1, prob = network_extra$weight)
    
    if(network_extra$accepted[i] == 1){     # if the samples edge has already been accepted skip to the next sampling
      next
    }
    
    cat("index: ", i, " current num of edges: ", length(network_extra[, 1]), "\n")
    
    #### check the degree of the nodes of the selected edge
    #### if removing the edge lowers the connectivity significantly or disconnects the network - skip
    fromDeg = g2_degree[network_extra[i, ]$ind1]
    toDeg = g2_degree[network_extra[i, ]$ind2]
    
    cat("from deg: ", fromDeg, " to degree: ", toDeg, "\n")
    
    if(fromDeg <= 3 || toDeg <= 3){
      if(dist_to_boundary[network_extra[i, ]$ind1] > 35 && dist_to_boundary[network_extra[i, ]$ind2] > 35){     # we allow the degree of the nodes close to the boundary to be lower than the rest
        cat("Edge not removed [degree constraint]...\n\n")
        network_extra$accepted[i] = 1
        next
      }
    }
    
    #### compute the probability of acceptance or rejection for the chosen edge
    predictedProb = predictKDE(orgKDE_both, network_extra[i, ], 3)
    uniRandom = predictKDE(triKDE_both, network_extra[i, ], 3)

    cat("prediction: ", predictedProb, " uniform random: ", uniRandom, "\n")
    
    #### accept the edge if the probability of its being in the original network is higher than the probability of its being in the triangulation 
    if(uniRandom <= predictedProb){
      cat("edge accepted...\n")
      network_extra$accepted[i] = 1
    }
    else{
      #### check for network connectivity
      network_temp = network_extra[-c(i), ]
      g2 = make_empty_graph() %>% add_vertices(ganglia_ppp$n)
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
        triData = data.frame(angle=network_extra$anglecomp, len=network_extra$euclidDist, weight=network_extra$weight)
        
        #triKDE_angle = computeKDE(triData, 1) # 1: angle, 2: length, 3: both
        #triKDE_length = computeKDE(triData, 2) # 1: angle, 2: length, 3: both
        triKDE_both = computeKDE(triData, 3) # 1: angle, 2: length, 3: both
      }
      else{
        cat("Edge not removed [connectivity constraint]...\n\n")
        network_extra$accepted[i] = 1
      }
    }
    noChange = noChange + 1
    cat("\n")
  }
  
  return(network_extra)
}


###################################################################################################
generateNetworkEdges <- function(ganglia_ppp, branch.all, orgKDE_angle, orgKDE_length, orgKDE_both,
                                 meshedness, network_density, compactness,
                                 sample_id, output_folder_path){
  
    #### constructing the deterministic Delaunay triangulation as the initial ganglionic network
    triangulation_info_list = deterministicEdges(ganglia_ppp, branch.all, sample_id, output_folder_path)
    
    #### returned values
    network_extra1 = triangulation_info_list[[1]]
    triKDE_angle = triangulation_info_list[[2]]
    triKDE_length = triangulation_info_list[[3]]
    triKDE_both = triangulation_info_list[[4]]
    g2_degree = triangulation_info_list[[5]]
    
    #### remove edges from the initial triangulation by rejection sampling
    network_extra = rejectionSampling(ganglia_ppp, network_extra1, g2_degree, orgKDE_angle, orgKDE_length, orgKDE_both,
                                    triKDE_angle, triKDE_length, triKDE_both, meshedness, network_density, compactness,
                                    sample_id, output_folder_path)
    
    #### observe if the distribution under consideration covers the targeted distribution
    par(mar=c(5, 2, 1, 1))
    plot(density(branch.all$angle), col="red", lty=2, xlab="Edge angle", ylab="Density", main="", xlim=c(-100,250), ylim=c(0, 0.02))
    lines(density(network_extra1$anglecomp), col="blue")
    lines(density(network_extra1$anglecomp), col="black", lty=3)
    legend(x=0, y=0.02, legend=c("ENS angle", "Triangulation angle", "Simulated angle"), 
         col=c("red", "blue", "black"), 
         lty=c(2, 1, 3))
    
    plot(density(branch.all$euclid), col="red", lty=2, xlab="Edge Length", ylab="Density", main="")
    lines(density(network_extra1$euclidDist), col="blue")
    lines(density(network_extra1$euclidDist), col="black", lty=3)
    legend(x=300, y=0.006, legend=c("ENS edge len", "Triangulation edge len", "Simulated edge len"), 
         col=c("red", "blue", "black"), 
         lty=c(2, 1, 3))
    
    #### visualize the bivariate distribution of the trimmed triangulation
    den3d = kde2d(network_extra$anglecomp, network_extra$euclidDist)
    persp(den3d, theta = -45, phi = 30, xlab="angle", ylab="edge len",
        ticktype = "detailed", shade = 0.75, col="lightblue")
    
    #### create a graph from sampled Delaunay triangulation
    g2 = make_empty_graph() %>% add_vertices(ganglia_ppp$n)
    g2 = add_edges(as.undirected(g2), as.vector(t(as.matrix(network_extra[,5:6]))))
    
    #### display as corresponding ppp and linnet
    g2_lin = linnet(ganglia_ppp, edges=as.matrix(network_extra[, 5:6]))
    
    par(mar=c(0,0,0,0), oma=c(0,0,0,0))
    plot(ganglia_ppp, cex=1, pch=20, main="", bg=1)
    plot(g2_lin, add=T)
    
    #### creating .pptx file to store the plots
    doc = read_pptx()
    
    degree_frame_2 = as.data.frame(degree(g2))
    colnames(degree_frame_2)[colnames(degree_frame_2) == 'degree(g2)'] = 'deg'
    
    ggplot(degree_frame_2, aes(x=deg)) +
        geom_histogram(aes(y=..density..), colour="grey", fill="grey", binwidth = 0.5)+
        geom_density(alpha=1, colour="black", size=1.5) +
        scale_x_continuous(limits=c(0, 10), breaks = seq(0,10, by=1))+
        labs(x = "Degree of the vertices", y = "Density", color = "")
    
    ggplot(network_extra, aes(x=weight)) +
        geom_density(alpha=1, colour="black", size=1.5) +
        labs(x = "weight of the edges", y = "Density", color = "")
    
    marked_g2  = graph_from_data_frame(network_extra[,5:6], directed=FALSE)
    degs = degree(marked_g2)
    ord = order(as.numeric(names(degs)))
    degs = degs[ord]
    marked_ganglia_ppp = ganglia_ppp
    marks(marked_ganglia_ppp) = factor(degs)
    marked_ganglia_ppp$markformat = "factor"
    
    marked_g2_lin = linnet(marked_ganglia_ppp, edges=as.matrix(network_extra[, 5:6]))
    
    branch.lpp = lpp(marked_ganglia_ppp, marked_g2_lin )
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(code = plot(branch.lpp, main="Simulated Network", cex=1, pch=21, bg=c(1,2,3,4,5,6,7,8,"white","forestgreen","tan3"))), location = ph_location_fullsize())
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(code = ({plot(marked_g2_lin)})), location = ph_location_fullsize())
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(code =({plot(density(branch.all$angle), col="red", lty=2, lwd=2, xlab="Edge angle (degree)", ylab="Density", main="",
                                      cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, xlim=c(-50, 250), ylim=c(0,0.015))
                                    lines(density(network_extra1$anglecomp), lwd=2, col="blue")
                                    lines(density(network_extra$anglecomp), lwd=2, col="black", lty=6)
                                    legend(x=80, y=0.015, legend=c("Real ENS", "Delaunay triangulation", "Simulated network"),
                                           cex=2, col=c("red", "blue", "black"),
                                           lty=c(2, 1, 6))})), location = ph_location_fullsize())
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(code = ({plot(density(branch.all$euclid), col="red", lty=2, lwd=2, xlab="Edge length", ylab="Density", main="",
                                       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
                                    lines(density(network_extra1$euclidDist), lwd=2, col="blue")
                                    lines(density(network_extra$euclidDist), lwd=2, col="black", lty=6)
                                    legend(x=240, y=0.006, legend=c("Real ENS", "Delaunay triangulation", "Simulated network"),
                                        cex=2, col=c("red", "blue", "black"),
                                        lty=c(2, 1, 6))})), location = ph_location_fullsize())
    
    #### creating a directory to save the simulated network related files
    sim_net_path = paste(output_folder_path, "Simulated Network/", sep="")
    if (!dir.exists(sim_net_path)) {dir.create(sim_net_path, recursive=TRUE)}
    
    print(doc, target = paste(sim_net_path, sample_id, "_Final_Network_Simulated.pptx", sep=""))
    
    write.csv(data.frame(x1=network_extra$x1, y1=network_extra$y1, x2=network_extra$x2, y2=network_extra$y2), 
            paste(sim_net_path, sample_id, "_Simulated_edge_coord.csv", sep=""), 
            row.names = F)
    
    #### compute EMD between the original ENS and the simulated network
    
    #### orgKDE holds the kernel density estimation;
    #### orgKDE2 holds the same info in a convenient way to be used while calculating the Earth Mover's Distance
    data = branch.all[, 7:8]
    orgKDE = computeKDE(data, 3)
    predictionOrg = predict(orgKDE, x=data) 
    orgKDE2 = as.matrix(data.frame(z=predictionOrg, angle=orgKDE$x[, 1], edgelen=orgKDE$x[, 2]))
    
    #### triKDE holds the kernel density estimation;
    #### triKDE2 holds the same info in a convenient way to be used if calculating the Earth Mover's Distance
    triData = data.frame(angle=network_extra$anglecomp, edgelen=network_extra$euclidDist)
    triKDE = computeKDE(triData, 3)
    predictionTri = predict(triKDE, x=triData)
    triKDE2 = as.matrix(data.frame(z=predictionTri, angle=triData[, 1], edgelen=triData[, 2]))
    
    EMD = emd(A=orgKDE2, B=triKDE2, dist="euclidean")
    cat("Earth Mover's Distance: ", EMD, "\n")
    
    return(list(network_extra, g2_lin, EMD))
}