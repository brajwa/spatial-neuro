#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path",
             "resample")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)


####==================================================================================
####this R code manually generates all the 29 graphlets mentioned in the paper: 
####https://watermark.silverchair.com/bth436.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAnYwggJyBgkqhkiG9w0BBwagggJjMIICXwIBADCCAlgGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMDNnCZZRBjuvsadSlAgEQgIICKTooU99bAe9Gb_5Hk-AOqCPUvf4sqI_lWTXhmmkMHpELKopjGJTZe_5XtraZkX8VnwCGsp89Nt6nR-_Npp30MsT4BwZiHGWBk3k6H0pdP_0Nvo2iPw55e6S4NvcAQZt_W853wuXBsEmHuWx0mAaHvjZpf7SN_tyI6BqN9ZEYrDKY-ArP7APpux1D1j7QVfB2xl5WInibdKjzV0cnni2xcsgEeUGTv2iY5xcrznTLm_wjiI36hG2FlvpBbu0lru9BQmG5Ko6W7oDfJwZ8oG12aFdRJTjr13zS17QREKGHt50QtMxsvyTsfcyVDZXQfeHk3g_XEvtF09nOjB28RkJQ0VRJIS6ZMbYJXVaaaQtyUPUkRQm3vvfuYsdVyaJUEknOzW9Sb02jU-rgnU5v_KdOS5f00jeDGofW30rECpWsHxHIxnhBRgHPX5hwZLTMqhTTDes_VgTgwkt5feIIqdNdrQVfeEbd9mnctCVUBRnAr6vP6mp4vPaAu7KQM54Q1jdag3PoxKzAPtLCXcDjqy3k-_x6KxuzRNZszJQtwiozLWIfrKKDFOuCsOrthcJXvlgyYsgP7cOFvqrPm2Wc9KlSzlljWEvsBZSRNtsxFrBAOovhKulZQ5mEpMeoc02pkZZgGucNeytVMICl8TR1t-uhvaOaDtaKT4vNuo0q-HSVqoTyA7nL2ardaMtBMTY2G-imr_6tgLv_zg_hwMFf8Ap9nHdm3fKtt2V29no

####to store all the 29 graphlets, 3-5 node connected networks
graphlet = c()
graphlet_iso_class = c()

####3-node graphlets
graphlet_instance = make_ring(n=3, circular = FALSE)
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n = 3)
graphlet[[length(graphlet) + 1]] = graphlet_instance

####4-node graphlets
graphlet_instance = make_ring(n=4, circular = FALSE)
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_star(n = 4, mode = "undirected")
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=4)
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=3)%>%
  add_vertices(1)%>%
  add_edges(c(3,4))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=4)%>%
  add_edges(c(2,4))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=3)%>%
  add_vertices(1)%>%
  add_edges(c(1,4, 2,4, 3,4))
graphlet[[length(graphlet) + 1]] = graphlet_instance

####5-node graphlets
graphlet_instance = make_ring(n=5, circular = FALSE)
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=3, circular = FALSE)%>%
  add_vertices(2)%>%
  add_edges(c(3,4, 3,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_star(n = 5, mode = "undirected")
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=3)%>%
  add_vertices(2)%>%
  add_edges(c(2,4, 3,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=3, circular = FALSE)%>%
  add_vertices(2)%>%
  add_edges(c(3,4, 3,5, 4,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=3)%>%
  add_vertices(2)%>%
  add_edges(c(3,4, 3,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=5)
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=4)%>%
  add_vertices(1)%>%
  add_edges(c(5,4))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=4)%>%
  add_vertices(1)%>%
  add_edges(c(2,4, 5,4))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=3)%>%
  add_vertices(2)%>%
  add_edges(c(3,4, 3,5, 4,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=4)%>%
  add_vertices(1)%>%
  add_edges(c(1,3, 5,4))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=4)%>%
  add_vertices(1)%>%
  add_edges(c(1,5, 3,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=4)%>%
  add_vertices(1)%>%
  add_edges(c(3,5, 4,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=4)%>%
  add_vertices(1)%>%
  add_edges(c(2,4, 2,5, 4,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=3)%>%
  add_vertices(2)%>%
  add_edges(c(1,4, 2,4, 3,4, 2,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=4)%>%
  add_vertices(1)%>%
  add_edges(c(2,4, 2,5, 3,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=4)%>%
  add_vertices(1)%>%
  add_edges(c(1,5, 3,5, 2,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=4)%>%
  add_vertices(1)%>%
  add_edges(c(2,4, 2,5, 3,5, 4,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=4)%>%
  add_vertices(1)%>%
  add_edges(c(1,5, 2,5, 3,5, 4,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_ring(n=4)%>%
  add_vertices(1)%>%
  add_edges(c(1,3, 2,4, 2,5, 3,5, 4,5))
graphlet[[length(graphlet) + 1]] = graphlet_instance

graphlet_instance = make_full_graph(n=5)
graphlet[[length(graphlet) + 1]] = graphlet_instance

#### computing the isomorphism class of the graphlets; because ordering them in the increasing order of isomorphism class is important
for (i in c(1:29)) {
  graphlet_iso_class[length(graphlet_iso_class) + 1] = isomorphism_class(graphlet[[i]])
}

graphlet_ord = graphlet[order(graphlet_iso_class)]

#### graphlet frequency count in the given ENS network for size 3, 4 and 5
graphlet_count = c(na.omit(c(motifs(g1, size=3), motifs(g1, size=4), motifs(g1, size=5)))) # using the updated library function
rel_graphlet_count = graphlet_count / sum(graphlet_count)

graphlet_count
plot(rel_graphlet_count, type="l")

# #### manually counting the graphlet frequency
# pattern = graphlet_ord[[1]]
# plot(pattern)
# 
# iso = subgraph_isomorphisms(pattern, g1)      # takes a while
# iso_2 = lapply(iso, sort)
# iso_3 = lapply(iso_2, as.vector)
# iso_4 = unique(iso_3)
# 
# motifs = lapply(iso_4, function (x) { induced_subgraph(g1, vids = x, impl = "copy_and_delete") })
# length(motifs)

####generate ER random graphs with the same number of nodes and edges as the ENS network
how_many_ER = 50
ER_graphs = vector(mode = "list", length = how_many_ER)
ER_graphlet_count = matrix(nrow = how_many_ER, ncol = 29)    # a matrix

for (j in 1:how_many_ER) {
  ER_graphs[[j]] = erdos.renyi.game(N, E, type = "gnm")
  #plot(ER_graphs[[j]])
  
  ER_graphlet_count_j = c(na.omit(c(motifs(ER_graphs[[j]], size=3), motifs(ER_graphs[[j]], size=4), motifs(ER_graphs[[j]], size=5))))
  print(ER_graphlet_count_j)
  
  ER_graphlet_count[j, ] = ER_graphlet_count_j
}

ER_graphlet_count_mean = colMeans(ER_graphlet_count)
ER_graphlet_count_var =  colVars(ER_graphlet_count)

ER_graphlet_count_mean

z_score = (graphlet_count - ER_graphlet_count_mean)/sqrt(ER_graphlet_count_var)
z_score