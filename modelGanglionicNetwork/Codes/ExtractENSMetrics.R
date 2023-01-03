#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)


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
#### choose one of the ganglionic network samples with file chooser below
setwd(paste(parent, "Data/Branch Information (in um)/", sep=""))
path_to_branch_info = file.choose()

#### extract which sample we will be working on from the file name (the data files are to be saved in a particular format)
sample_id = strsplit(path_to_branch_info, "\\\\")[[1]]
sample_id = sample_id[length(sample_id)]
sample_id = (strsplit(sample_id, "_")[[1]])[1]

#### creating an output folder for the working sample
output_folder_path = paste(parent, "Outputs/", sample_id, "/", sep="")
if (!dir.exists(output_folder_path)) {dir.create(output_folder_path, recursive=TRUE)}


data_struct_list = constructDataStruct(sample_id, parent, path_to_branch_info, output_folder_path)

#### the returned values
branch.all = data_struct_list[[1]]
branch.ppp = data_struct_list[[2]]
branch.lpp = data_struct_list[[3]]
g1 = data_struct_list[[4]]
hardcoreStrauss_model_param = data_struct_list[[5]]

#### spatial point pattern metrics
num_points_pp = branch.ppp$n
area_pp = summary(branch.ppp)$window$area
ints_pp = intensity(branch.ppp)

closest_point_dist = min(nndist(branch.ppp))
farthest_point_dist = max(pairdist(branch.ppp))


#### spatial network metrics

# degree distribution
table_degree = table(degree(g1))

#### calculating the edge probability p from the degree distribution
#### to use while generating ER random graph G(n, p)
num_nodes = 0
sum = 0
for(i in 1:length(table_degree)){
  num_nodes = num_nodes + table_degree[[i]]
  sum = sum + (i * table_degree[[i]])
}
edge_probability = divide_by(sum, num_nodes^2)
print(edge_probability)

#### Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
#### This is sometimes also called the clustering coefficient.
cluster_coeff = transitivity(g1, type = "global")
print(cluster_coeff)

#### alpha, gamma, psi
N = num_points_pp
E = length(branch.all$x1)
A = area_pp
L = sum(branch.all$euclid)

meshedness = (E-N+1)/((2*N)-5)
network_density = E/((3*N)-6)
compactness = 1- ((4*A)/(L-(2*sqrt(A)))^2)