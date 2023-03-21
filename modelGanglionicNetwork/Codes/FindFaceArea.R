# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("RBGL")

#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim", "RBGL")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)

# x = c(0,4,1,3,4,2)
# y = c(3,2,1,0,5,3)
# 
# pp = ppp(x, y, window = owin(xrange = c(0,5), yrange = c(0,5)))
# plot(pp)
# 
# del_tri = deldir(data.frame(x=pp$x, y=pp$y))
# 
# g2 = make_empty_graph() %>% add_vertices(pp$n) 
# g2$layout = as.matrix(data.frame(x=pp$x, y=pp$y))
# g2 = add_edges(as.undirected(g2), as.vector(t(as.matrix(del_tri$delsgs[1:8,5:6]))))
# 
# plot(g2)

# load ENS from data

g <- as_graphnel(g1) ## Convert igraph object to graphNEL object for planarity testing
boyerMyrvoldPlanarityTest(g)

face_list = planarFaceTraversal(g)
face_node_count = sapply(face_list, length)

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

combined_f_a_l[[length(combined_f_a_l) + 1]] = face_area_list

#### save list of list as xlsx, each list in separate sheet
#### 11 distal, 9 middle, 15 proximal
# sheet_name_list = as.character(c(1:length(combined_f_a_l)))
# combined_f_a_file_name = "C:/Users/sanja/Documents/GitHub/spatial-neuro/modelGanglionicNetwork/Outputs/ENSMouse/combined_face_area.xlsx"
# 
# write.xlsx(setNames(as.list(lapply(combined_f_a_l, data.frame)), names(sheet_name_list)), file=combined_f_a_file_name)


#############################################################################################33
# ggplot() +
#   geom_density(aes(x=combined_f_a_l[[1]], y=after_stat(density)), alpha=1, colour="black", linewidth=1) + #black=distal
#   geom_density(aes(x=combined_f_a_l[[2]], y=after_stat(density)), alpha=1, colour="black", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[3]], y=after_stat(density)), alpha=1, colour="black", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[4]], y=after_stat(density)), alpha=1, colour="black", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[5]], y=after_stat(density)), alpha=1, colour="black", linewidth=1) +
# 
#   geom_density(aes(x=combined_f_a_l[[6]], y=after_stat(density)), alpha=1, colour="black", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[7]], y=after_stat(density)), alpha=1, colour="black", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[8]], y=after_stat(density)), alpha=1, colour="black", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[9]], y=after_stat(density)), alpha=1, colour="black", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[10]], y=after_stat(density)), alpha=1, colour="black", linewidth=1) +
# 
#   geom_density(aes(x=combined_f_a_l[[11]], y=after_stat(density)), alpha=1, colour="black", linewidth=1) +
# 
#   geom_density(aes(x=combined_f_a_l[[12]], y=after_stat(density)), alpha=1, colour="red", linewidth=1) + #red=middle
#   geom_density(aes(x=combined_f_a_l[[13]], y=after_stat(density)), alpha=1, colour="red", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[14]], y=after_stat(density)), alpha=1, colour="red", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[15]], y=after_stat(density)), alpha=1, colour="red", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[16]], y=after_stat(density)), alpha=1, colour="red", linewidth=1) +
# 
#   geom_density(aes(x=combined_f_a_l[[17]], y=after_stat(density)), alpha=1, colour="red", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[18]], y=after_stat(density)), alpha=1, colour="red", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[19]], y=after_stat(density)), alpha=1, colour="red", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[20]], y=after_stat(density)), alpha=1, colour="red", linewidth=1) +
#   
#   geom_density(aes(x=combined_f_a_l[[21]], y=after_stat(density)), alpha=1, colour="dodgerblue", linewidth=1) + #dodgerblue=proximal
#   geom_density(aes(x=combined_f_a_l[[22]], y=after_stat(density)), alpha=1, colour="dodgerblue", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[23]], y=after_stat(density)), alpha=1, colour="dodgerblue", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[24]], y=after_stat(density)), alpha=1, colour="dodgerblue", linewidth=1) +
#   geom_density(aes(x=combined_f_a_l[[25]], y=after_stat(density)), alpha=1, colour="dodgerblue", linewidth=1) +
#   labs(x = "Area of face", y = "Density", color = "")
