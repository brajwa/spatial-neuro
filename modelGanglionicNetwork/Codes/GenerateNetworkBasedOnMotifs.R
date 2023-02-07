########################################################################
#### computes an edge weight based on the degree of its two end vertices
#### g2_degree is the list of degree of each vertex in the graph under consideration
#### x is the edge whose weight is being computed, the 5th and 6th index provides the indices of the end vertices
computeEdgeWeight <- function(g2_degree, x){
  return(g2_degree[x[5]] + g2_degree[x[6]])
}


#### for now, use the original ENS point pattern, later it will be replaced by simulated point pattern from fitted point process
ganglia_ppp = branch.ppp


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

#### degree list of the newly constructed Delaunay graph
g2_degree = degree(g2)
table(g2_degree)

network_extra$weight = apply(network_extra, 1, function(x) computeEdgeWeight(g2_degree, x))
network_extra$weight = network_extra$weight / sum(network_extra$weight)

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
plot(marked_branch.lpp, main="triangulation", cex=1, pch=21, bg=c(3,4,5,6,7,8,"white","forestgreen","tan3", "black", "red3", "mediumpurple"))

#### graphlet frequency count in the Delaunay graph for size 3, 4 and 5
del_graphlet_count = c(na.omit(c(motifs(g2, size=3), motifs(g2, size=4), motifs(g2, size=5)))) # using the updated library function
rel_del_graphlet_count = del_graphlet_count / sum(del_graphlet_count)

plot(rel_del_graphlet_count, type="l")

########################################################################
#### some initialization
rejected = 0
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
  
  cat(meshedness, " ", network_density, " ", compactness, "\n")
  cat(mesh, " ", n_density, " ", compact, "\n\n")
  
  if(mesh <= meshedness || n_density <= network_density || compact <= compactness){   # if network is getting too sparse
    cat("Network got sparse\n")
    break
  }
  
  #### sample an edge index based on the degree-based weights assigned to the edges
  i = sample.int(length(network_extra[, 1]), 1, prob = network_extra$weight)
  
  if(network_extra$accepted[i] == 1){     # if the samples edge has already been accepted skip to the next sampling
    noChange = noChange + 1
    next
  }
  
  cat("index: ", i, " current num of edges: ", length(network_extra[, 1]), "\n")
  
  #### check for network connectivity
  network_temp = network_extra[-c(i), ]
  g3 = make_empty_graph() %>% add_vertices(ganglia_ppp$n)
  g3 = add_edges(as.undirected(g3), as.vector(t(as.matrix(network_temp[,5:6]))))
  
  if(is_connected(g3)){
    #### graphlet frequency count in the new graph for size 3, 4 and 5
    new_graphlet_count = c(na.omit(c(motifs(g3, size=3), motifs(g3, size=4), motifs(g3, size=5)))) # using the updated library function
    rel_new_graphlet_count = new_graphlet_count / sum(new_graphlet_count)
    
    if(l2norm(rel_new_graphlet_count, rel_graphlet_count) <= l2norm(rel_del_graphlet_count, rel_graphlet_count)){ # this checking can be weighted
      #### adjust the degree of the nodes
      g2_degree[network_extra[i, ]$ind1] = g2_degree[network_extra[i, ]$ind1] - 1
      g2_degree[network_extra[i, ]$ind2] = g2_degree[network_extra[i, ]$ind2] - 1
      
      network_extra = network_temp
      
      cat("edge rejected...\n")
      rejected = rejected + 1
      noChange = 0
      
      #### recompute the degree based edge weights
      network_extra$weight = apply(network_extra, 1, function(x) computeEdgeWeight(g2_degree, x))
      network_extra$weight = network_extra$weight / sum(network_extra$weight)
      
      #### update the graphlet profile of Delaunay
      del_graphlet_count = new_graphlet_count
      rel_del_graphlet_count = rel_new_graphlet_count
      
      marked_g2  = graph_from_data_frame(network_extra[,5:6], directed=FALSE)
      degs = degree(marked_g2)
      ord = order(as.numeric(names(degs)))
      degs = degs[ord]
      marked_ganglia_ppp = ganglia_ppp
      marks(marked_ganglia_ppp) = factor(degs)
      marked_ganglia_ppp$markformat = "factor"
      
      marked_g2_lin = linnet(marked_ganglia_ppp, edges=as.matrix(network_extra[,5:6]))
      
      marked_branch.lpp = lpp(marked_ganglia_ppp, marked_g2_lin )
      plot(marked_branch.lpp, main="triangulation", cex=1, pch=21, bg=c(3,4,5,6,7,8,"white","forestgreen","tan3", "black", "red3", "mediumpurple"))
      
    }else{
      noChange = noChange + 1
    }
    
  }else{
    network_extra$accepted[i] = 1
    noChange = noChange + 1
  }
  cat("\n")
}

plot(rel_del_graphlet_count, type="l")
