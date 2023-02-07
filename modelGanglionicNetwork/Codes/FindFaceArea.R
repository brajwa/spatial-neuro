if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RBGL")

#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim", "RBGL")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)

# pp = rpoispp(3, win=owin(c(0,1), c(0,1)))
# plot(pp)
# 
# del_tri = deldir(data.frame(x=pp$x, y=pp$y))
# g2 = make_empty_graph() %>% add_vertices(pp$n)
# g2 = add_edges(as.undirected(g2), as.vector(t(as.matrix(del_tri$delsgs[,5:6]))))
# 
# plot(g2)

# load ENS from data

g <- as_graphnel(g1) ## Convert igraph object to graphNEL object for planarity testing
boyerMyrvoldPlanarityTest(g)

face_list = planarFaceTraversal(g)
face_node_count = sapply(face_list, length)

face_list = face_list[-which.max(face_node_count)]
face_node_count = face_node_count[-which.max(face_node_count)]

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