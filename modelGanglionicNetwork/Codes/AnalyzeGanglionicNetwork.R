#################################################
range01 <- function(x){(x-min(x))/(max(x)-min(x))}


###############################################################################################################
# Computes the angle of an edge with the horizontal x axis in degree; the input parameter x has four components;
# the coordinates of the end points of the edge: x1, y1, x2, y2.
calcAngle <- function(x){
  nom = x[4]-x[2]
  denom = x[3]-x[1]
  angle = as.numeric(atan(nom/denom))*180/pi
  
  return((angle) %% 180)
}


#######################################################################################
# Computes the Euclidean Distance between the two nodes of an edge aka the edge length; 
# the input parameter x has four components; the coordinates of the end points of the edge: x1, y1, x2, y2.
calcDist <- function(x){
  return(sqrt( ((x[1]-x[3])^2) + ((x[2]-x[4])^2) ) )
}


####################################################################################
# computes the distribution of angle and/or length for original or synthetic network
# parameter: data is a matrix with 3 columns; column 1 has the angle; column 2 has the length; column 3 has the degree-based weight
# parameter: mode is an integer deciding which data columns to use
computeKDE <- function(data, mode){
  if(mode == 1){
    data = as.matrix(data[, 1])     # edge angle
    return(kde(x=data))
  }
  if(mode == 2){
    data = as.matrix(data[, 2])     # edge length
    return(kde(x=data))
  }
  if(mode == 3){
    data = as.matrix(data[, 1:2])   # both angle and length
    return(kde(x=data))
  }
  if(mode == 4){
    data = as.matrix(data[, 3])     # degree-based weight
    return(kde(x=data, h=0.01))
  }
  if(mode == 5){
    data = as.matrix(data[, 1:3])   # all
    return(kde(x=data))
  }
}


#########################################################################################
constructDataStruct <- function(sample_id, parent, branch_info_path, output_folder_path){
    #### extracting branch information file type. usually saved as .csv, checking for .xlsx in case.
    path_tokens = strsplit(branch_info_path, "\\.")
    file_extension = path_tokens[[1]][length(path_tokens[[1]])]
    
    branch = data.frame()
    
    if(file_extension == "csv"){
        branch = read.csv(branch_info_path)
    }
    else if(file_extension == "xlsx"){
        branch = read.xlsx(branch_info_path)
    }
    else{
        stop("Invalid input file type! Has to be csv or xlsx.\n")
    }
    
    #### checking Branch.length column: if the length is 1, it is not really  a branch, removing them
    branch.adj = branch[which(branch$Branch.length!=1),]
    
    hardcoreStrauss_model_param = mean(branch.adj$Branch.length)
    
    #### constructing a dataframe with the branch nodes, consider y coordinates are negated beforehand.
    branch.sim = data.frame(x1=branch.adj$V1.x, y1=branch.adj$V1.y, x2=branch.adj$V2.x, y2=branch.adj$V2.y)
    
    #### removing duplicate branches
    temp = which(branch.sim[,1]!=branch.sim[,3] &  branch.sim[,2]!=branch.sim[,4])
    branch.sim = branch.sim[temp,]
    
    #### create a unique id for each of the branches so, we can find the vertices
    branch.id = matrix(data=factor(as.numeric(as.matrix(branch.sim))), ncol=4)
    branch.id = data.frame(P1=paste(branch.id[,1],"-",branch.id[,2], sep=""), P2=paste(branch.id[,3],"-",branch.id[,4], sep=""))
    branch.id = c(as.matrix(branch.id)) %>% factor() #%>% unclass()
    levels(branch.id)  = 1:length(levels(branch.id))
    branch.id = data.frame(matrix(ncol=2, data=unclass(branch.id)))
    
    #### put together all info and extract position of all the vertices
    branch.all = data.frame(branch.sim, branch.id); colnames(branch.all) <- c("x1","y1","x2","y2","n1","n2")
    
    branch.xy = data.frame(x=c(branch.all$x1, branch.all$x2), y=c(branch.all$y1, branch.all$y2), n =c(branch.all$n1, branch.all$n2) )
    branch.xy = distinct(branch.xy)
    
    #### uncomment if required: normalize the branch length to [0,1] for convenience only
    #branch.xy[,1:2] = range01(as.vector(branch.xy[,1:2]))
    
    #### ordered by the vertex number
    branch.xy.ord = branch.xy[order(branch.xy$n),]
    
    #### creating point and line patterns using spatstat library
    #### creating a point pattern
    branch.ppp = ppp(x=branch.xy.ord$x, y=branch.xy.ord$y, window=owin(xrange=c(min(branch.xy.ord$x), max(branch.xy.ord$x)), yrange=c(min(branch.xy.ord$y), max(branch.xy.ord$y))))
    branch.ppp = unique.ppp(branch.ppp)
    
    #### plot the nodes
    #plot(branch.ppp, cex=1, pch=20, main="", bg=1)
    
    #### creating a line pattern
    branch.lin = linnet(branch.ppp, edges=as.matrix(branch.all[,5:6]))
    
    #### plot the branches
    #plot(branch.lin, add=T)
    
    #### creating the graph representation, nodes' color and shape represent their degree
    g1  = graph_from_data_frame(branch.all[,5:6], directed=FALSE)
    
    #### Location information is taken from the processed branches
    #plot(g1, vertex.size=2*degree(g1), vertex.label=NA, vertex.shape="circle", layout=as.matrix(branch.xy[,1:2]), vertex.color=degree(g1))
    
    #### creating a list representing degree for every node an use it in spatstat pattern
    degs = (degree(g1))
    ord = order(as.numeric(names(degs)))
    degs = degs[ord]
    
    #### required sometimes
    #marks(branch.ppp) = factor(degs)
    #branch.ppp$markformat = "factor"
    
    #### Create point pattern on linear network to store the information
    #### There are some problems with this approach. 
    #### represent the data as point patterns on linear network BEWARE: if there are any unconnected
    #### nodes, they will be removed! Consequently, the numbering of vertices changes and the colors are
    #### going to be wrong It is technically not a bug: the lpp format is used to represent points on the
    #### network ONLY
    branch.lpp = lpp(branch.ppp, branch.lin )
    
    return(list(branch.all, branch.ppp, branch.lpp, g1, hardcoreStrauss_model_param))
}


##########################################
summaryStat <- function(sample_id, parent, branch_info_path, output_folder_path){
    #### constructing spatial data structures for subsequent analysis
    data_struct_list = constructDataStruct(sample_id, parent, branch_info_path, output_folder_path)
    
    #### the returned values
    branch.all = data_struct_list[[1]]
    branch.ppp = data_struct_list[[2]]
    branch.lpp = data_struct_list[[3]]
    g1 = data_struct_list[[4]]
    hardcoreStrauss_model_param = data_struct_list[[5]]
    
    #### creating .pptx file to store the plots and workbook for the statistics
    stats = createWorkbook() 
    doc = read_pptx()
    
    #### slide for point pattern
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(code = plot(branch.ppp, main=sample_id, cex=1, pch=20)), location = ph_location_fullsize())
    
    #### slide for network
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(code = plot(branch.lpp, main=sample_id, cex=0.5, pch=21)), location = ph_location_fullsize())
    
    #### computing summary statistics and adding slides for them
    
    #### G-function
    gtest = envelope(branch.ppp, Gest)
    
    ggobj = ggplot(data = gtest, aes(r)) + 
                geom_line(aes(y=theo, colour="csr"), size=2) +
                geom_ribbon(aes(ymin=lo, ymax=hi), alpha = 0.3)+
                geom_line(aes(y=obs, colour="observed"), size=2) + 
                theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
                      legend.text=element_text(size=30),
                      axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
                      axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
                xlab("Interaction distance") + ylab("Nearest Neighbor Function G") +
                ggtitle(sample_id)
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
    
    addWorksheet(stats, "Gest")
    writeData(stats, sheet = "Gest", x = gtest)
    
    #### F-function
    ftest = envelope(branch.ppp, Fest)
    
    ggobj = ggplot(data = ftest, aes(r)) + 
                geom_line(aes(y=theo, colour="csr"), size=2) +
                geom_ribbon(aes(ymin=lo, ymax=hi), alpha = 0.3)+
                geom_line(aes(y=obs, colour="observed"), size=2) + 
                theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
                      legend.text=element_text(size=30),
                      axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
                      axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
                xlab("Interaction distance") + ylab("Empty-space Function F")+
                ggtitle(sample_id)
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
    
    addWorksheet(stats, "Fest")
    writeData(stats, sheet = "Fest", x = ftest)
    
    #### J-function
    jtest = envelope(branch.ppp, Jest)
    
    ggobj = ggplot(data = jtest, aes(r)) + 
                geom_line(aes(y=theo, colour="csr"), size=2) +
                geom_ribbon(aes(ymin=lo, ymax=hi), alpha = 0.3)+
                geom_line(aes(y=obs, colour="observed"), size=2) + 
                theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
                      legend.text=element_text(size=30),
                      axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
                      axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
                xlab("Interaction distance") + ylab("Summary Function J")+
                ggtitle(sample_id)
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
    
    addWorksheet(stats, "Jest")
    writeData(stats, sheet = "Jest", x = jtest)
    
    #### Ripley's K-function
    ktest = envelope(branch.ppp, Kest)
    
    ggobj = ggplot(data = ktest, aes(r)) + 
                geom_line(aes(y=theo, colour="csr"), size=2) +
                geom_ribbon(aes(ymin=lo, ymax=hi), alpha = 0.3)+
                geom_line(aes(y=obs, colour="observed"), size=2) + 
                theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
                      legend.text=element_text(size=30),
                      axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
                      axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
                xlab("Interaction distance") + ylab("K Function")+
                ggtitle(sample_id)
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
    
    addWorksheet(stats, "Kest")
    writeData(stats, sheet = "Kest", x = ktest)
    
    #### Besag's centered inhomogeneous L-function
    ltest= envelope(branch.ppp, Linhom)
    
    ggobj = ggplot(data = ltest, aes(r)) + 
                geom_ribbon(aes(ymin=lo-r, ymax=hi-r), alpha = 0.3)+
                geom_line(aes(y=lo-r, colour="csr lower bound"), size=2) +
                geom_line(aes(y=hi-r, colour="csr upper bound"), size=2) +
                geom_line(aes(y=obs-r, colour="observed"), size=2) + 
                theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
                      legend.text=element_text(size=30),
                      axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
                      axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
                xlab("Interaction distance") + ylab(" Besag's centered L Function")+
                ggtitle(sample_id)
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
    
    addWorksheet(stats, "Lest")
    writeData(stats, sheet = "Lest", x = ltest)
    
    #### creating a directory to save the summary statistics related files
    sum_stat_path = paste(output_folder_path, "Summary Statistics/", sep="")
    if (!dir.exists(sum_stat_path)) {dir.create(sum_stat_path, recursive=TRUE)}
    
    saveWorkbook(stats, paste(sum_stat_path, sample_id, "_STAT_Values.xlsx", sep=""), overwrite = T)
    print(doc, target = paste(sum_stat_path, sample_id, "_STAT_Figures.pptx", sep=""))

}


####################################################################################
analyzeGanglia <- function(sample_id, parent, branch_info_path, output_folder_path){
    #### constructing spatial data structures for subsequent analysis
    data_struct_list = constructDataStruct(sample_id, parent, branch_info_path, output_folder_path)
    
    #### the returned values
    branch.all = data_struct_list[[1]]  #data frame with all information
    branch.ppp = data_struct_list[[2]]  #point pattern object (ganglia)
    branch.lpp = data_struct_list[[3]]  #line pattern object (ganglionic network)
    g1 = data_struct_list[[4]]          #graph of the ganglionic network
    hardcoreStrauss_model_param = data_struct_list[[5]] #model parameter
    
    #### summary of the constructed point pattern
    print(summary(branch.ppp))
    
    #### inspecting the L-function behavior of the point pattern, change Linhom to Lest if homogeneous version is preferred
    ltest= envelope(branch.ppp, Linhom)
    
    #### plot the L-function
    ggplot(data = ltest, aes(r)) + 
        geom_ribbon(aes(ymin=lo-r, ymax=hi-r), alpha = 0.3)+
        
        geom_line(aes(y=lo-r, colour="CSR lower bound"), size=2) +
        geom_line(aes(y=hi-r, colour="CSR upper bound"), size=2) +
        
        geom_line(aes(y=obs-r, colour="observed"), size=2) + 

        theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
              legend.text=element_text(size=30),
              axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
              axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
        
        xlab("Interaction distance") + ylab("Besag's centered L-Function")
    
    #### plot the line pattern
    par(mar=c(0,0,0,0), oma=c(0,0,0,0))
    plot(branch.lpp, main="", pch=20, cex=1)

    #### spatial model of the point pattern (ganglia)
    ganglia_model = ppm(branch.ppp ~ x+y, interaction = StraussHard(hardcoreStrauss_model_param))  #mean(branch.adj$Branch.length)
    #print(ganglia_model)

    #### extracting the model parameters
    param_model = rmhmodel(ganglia_model)
    
    Beta = param_model$par$beta     #stationary intensity parameter (not required in our case)
    Gamma = param_model$par$gamma   #interaction parameter
    R = param_model$par$r           #interaction distance
    H = param_model$par$hc          #hardcore distance
    window = branch.ppp$window
    
    return(list(Beta, Gamma, R, H, window, ganglia_model))
}


#############################################
analyzeBranch <- function(sample_id, parent, branch_info_path, output_folder_path){
    #### constructing spatial data structures for subsequent analysis
    data_struct_list = constructDataStruct(sample_id, parent, branch_info_path, output_folder_path)
    
    #### the returned values
    branch.all = data_struct_list[[1]]  #data frame with all information
    branch.ppp = data_struct_list[[2]]  #point pattern object (ganglia)
    branch.lpp = data_struct_list[[3]]  #line pattern object (ganglionic network)
    g1 = data_struct_list[[4]]          #graph of the ganglionic network
    hardcoreStrauss_model_param = data_struct_list[[5]] #model parameter
  
    #### creating .pptx file to store the plots and workbook for the data
    stats = createWorkbook() 
    doc = read_pptx()
    
    par(mar=c(0,0,0,0), oma=c(0,0,0,0))
    plot(branch.lpp, main="", pch=21, cex=1)

    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(code = plot(branch.lpp, main="", pch=20, cex=1)), location = ph_location_fullsize())
    
    #### extracting information about connectivity
    table_degree = table(degree(g1))
    cat("Degree table:\n", table_degree, "\n")
    
    degree_frame = as.data.frame(degree(g1))
    colnames(degree_frame)[colnames(degree_frame) == 'degree(g1)'] = 'deg'
    
    ggobj = ggplot(degree_frame, aes(x=deg)) + 
                geom_histogram(aes(y=..density..), colour="grey", fill="grey", binwidth = 0.5)+ 
                geom_density(alpha=1, colour="black", size=1.5) +
                scale_x_continuous(limits=c(0, 10), breaks = seq(0,10, by=1))+
                labs(x = "Degree of the vertices", y = "Density", color = "")
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
    
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
    
    #### calculate the edge angle (in degree) and plot the distribution
    #### not scaled
    branch.all$angle = (apply(branch.all, 1, function(x) calcAngle(x)))
    
    ggobj = ggplot(branch.all, aes(x=angle)) + 
                geom_histogram(aes(y=..density..), colour="grey", fill="grey", binwidth = 3)+
                geom_density(alpha=1, colour="black", size=1.5) +
                labs(x = "Edge angle", y = "Density", color = "")
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
    
    #### calculate the edge length and plot the distribution
    #### not scaled
    branch.all$euclid = (apply(branch.all, 1, function(x) calcDist(x)))
    
    ggobj = ggplot(branch.all, aes(x=euclid)) + 
                geom_histogram(aes(y=..density..), colour="grey", fill="grey", binwidth = 7)+
                geom_density(alpha=1, colour="black", size=1.5) +
                labs(x = "Edge length (Euclidean)", y = "Density", color = "")
    
    doc = add_slide(doc, "Blank", "Office Theme")
    doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
    
    #### check if the edge angle and length are correlated; they both are standardized value
    plot(branch.all$angle, branch.all$euclid)
    
    #### visualize the bivariate distribution of the network's edge angle and length
    den3d = kde2d(branch.all$angle, branch.all$euclid)
    persp(den3d, theta = -45, phi = 30, xlab="Edge angle", ylab="Edge length",
        ticktype = "detailed", shade = 0.75, col="lightblue")
    
    #### creating a directory to save the branch characteristics
    branch_char_path = paste(output_folder_path, "Branch Characteristics/", sep="")
    if (!dir.exists(branch_char_path)) {dir.create(branch_char_path, recursive=TRUE)}
    
    print(doc, target = paste(branch_char_path, sample_id, "_Branch_Figures.pptx", sep=""))
    
    #### alpha, gamma, psi
    N = branch.ppp$n
    E = length(branch.all$x1)
    A = (branch.ppp$window$xrange[2]-branch.ppp$window$xrange[1])*(branch.ppp$window$yrange[2]-branch.ppp$window$yrange[1])
    L = sum(branch.all$euclid)
    
    meshedness = (E-N+1)/((2*N)-5)
    network_density = E/((3*N)-6)
    compactness = 1- ((4*A)/(L-(2*sqrt(A)))^2)
    
    cat("meshedness: ", meshedness, "\n")
    cat("network_density: ", network_density, "\n")
    cat("compactness: ", compactness, "\n")
    
    #### assigning a weight to each edge based on the degree of its end nodes
    branch.all$weight = apply(branch.all, 1, function(x) computeEdgeWeight(degree(g1), x))
    branch.all$weight = range01(branch.all$weight)
    
    #### compute the univariate or bivariate density distribution of the edge angle and/ or edge length 
    #### orgKDE holds the kernel density estimation
    data = branch.all[, 7:9] # column 7 has angle, 8 has length, 9 has weight
    
    orgKDE_angle = computeKDE(data, 1) # 1: angle, 2: length, 3: both
    orgKDE_length = computeKDE(data, 2) # 1: angle, 2: length, 3: both
    orgKDE_both = computeKDE(data, 3) # 1: angle, 2: length, 3: both
    
    return(list(branch.all, orgKDE_angle, orgKDE_length, orgKDE_both, meshedness, network_density, compactness))

}