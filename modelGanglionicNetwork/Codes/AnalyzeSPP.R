#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim", "RBGL", "svglite", "ggrepel")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)


setwd("~/GitHub/spatial-neuro/modelGanglionicNetwork/Codes")

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
#### choose one/all of the ganglionic network samples with file chooser below
branch_info_folder = paste(parent, "Data/ENSMouse Branch Information (in um) v2.0/", sep="")
branch_info_files = list.files(branch_info_folder, recursive = TRUE, pattern = "\\.csv", full.names = TRUE)

figure_folder = paste(parent, "Outputs/ENSMouse/SPPFeature/", sep="")

#### before processing anything we need the max range along the y-axis for rescaling
# max_y = -Inf
# for (i in c(1: length(branch_info_files))) {
#     branch_data = read.csv(branch_info_files[i])
#     new_max_y = max(max(branch_data$V1.y), max(branch_data$V2.y))
#     
#     cat("i: ", i, "; current max y: ", max_y, "; new max y: ", new_max_y, "\n")
#     
#     if(max_y <= new_max_y){
#         cat("max y updated\n")
#         max_y = new_max_y
#     }
# }
max_y = 4539.812 # found by computation


pp_list = list()
### computing the max r values in inhomogeneous L function for each sample
max_inhom_l_r = c()
#g = ggplot()

spp_feature_list = c("R", "HC", "Gamma") 
select_feature = 1
cat("Spatial feature under consideration: ", spp_feature_list[select_feature], "\n")

columns = c("spp_r", "spp_hc", "spp_gamma", "ens_location","sample_id") 
feature_info_combined = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(feature_info_combined) = columns

for (i in c(1: length(branch_info_files))) { # 2,13,21: 3 samples of the 3 categories but same subject
    
    ens_location = strsplit(branch_info_files[i], "/")[[1]][11]
    sample_id = strsplit(strsplit(branch_info_files[i], "/")[[1]][12], "\\.")[[1]][1]
    cat("\n(", i, ") Location: ", ens_location, "\nSample Id: ", sample_id, "\n")
    
    #### creating an output folder for the working sample
    output_folder_path = paste(parent, "Outputs/ENSMouse/", sample_id, "/", sep="")
    if (!dir.exists(output_folder_path)) {dir.create(output_folder_path, recursive=TRUE)}
    
    
    data_struct_list = constructDataStruct(sample_id, parent, branch_info_files[i], output_folder_path, max_y)
    
    #### the returned values
    branch.all = data_struct_list[[1]]
    branch.ppp = data_struct_list[[2]] #scale the coordinate and re-center
    branch.lpp = data_struct_list[[3]]
    g1 = data_struct_list[[4]]
    hardcoreStrauss_model_param = data_struct_list[[5]]
    
    u_branch.ppp = unmark(branch.ppp)
    pp_list[[length(pp_list) + 1]] = u_branch.ppp
    
    plot(u_branch.ppp, main=sample_id, pch=21, cex=1, bg="black")
    
    #plot(branch.lpp, main="", pch=21, cex=1.2, bg=c("black", "red3", "green3", "orange", "dodgerblue", "white", "maroon1",
    #                                                       "mediumpurple"))
                                                    
    
    #### spatial point pattern metrics
    num_points_pp = u_branch.ppp$n
    area_pp = summary(u_branch.ppp)$window$area
    ints_pp = intensity.ppp(u_branch.ppp)
    
    closest_point_dist = min(nndist(u_branch.ppp))
    farthest_point_dist = max(pairdist(u_branch.ppp))
    
    cat("Number of points: ", num_points_pp, "\n",
        "Window area: ", area_pp, " um squared\n",
        "Homogeneous intensity of points: ", ints_pp, "\n",
        "Closest point distance: ", closest_point_dist, "\n",
        "Farthest point distance: ", farthest_point_dist, "\n")
    
    #L-function
    l_inhom = Linhom(u_branch.ppp, correction = "border")
    
    max_inhom_l_r[length(max_inhom_l_r) + 1] = max(l_inhom$r)

    print(ggplot(l_inhom)+
        geom_abline(slope=0, intercept = 0)+
        geom_vline(xintercept = 0)+

        geom_line(aes(x=r, y=border-r, colour=sample_id)) +

        theme(legend.position="top",  legend.title=element_blank(), legend.text=element_text(size=8),
              legend.box.margin=margin(0,-10,-10,-10),
              plot.title = element_text(hjust = 0.5, size=12),
              axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
              axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
              panel.background = element_rect(fill='white', colour='black'),
              panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +

        xlab(expression(paste("Interaction distance (r)"))) + ylab("L(r)") +
        ggtitle("Besag's centered inhomogeneous L-function") )
    
    #g <- g + geom_line(data=l_inhom, aes(x=r, y=(border-r)))
    
    range_R = data.frame(r=seq(closest_point_dist+0.001, max(l_inhom$r), by=max(l_inhom$r)/100))
    
    p = profilepl(s=range_R, f=StraussHard, u_branch.ppp~x+y, 
                  aic=TRUE, rbord = closest_point_dist+0.001)
    
    cat("Interaction distance: ", parameters(p)$r, "\n",
        "Hardcore distance: ", parameters(p)$hc, "\n",
        "Interaction coeff: ", parameters(p)$gamma, "\n")
    
    # -AIC and Gamma values
    p_info = data.frame(r=p$param$r, aic=p$prof, gamma=p$allcoef$Interaction)
    
    print(ggplot(data = p_info) + 
        geom_point(aes(x=r, y=aic), size=1.5) +
        geom_line(aes(x=r, y=aic), size=1) +
        
        geom_vline(xintercept = parameters(p)$r, size=1, lty=2, colour="red") +
        
            theme(legend.position="top",  legend.title=element_blank(), legend.text=element_text(size=8),
                  legend.box.margin=margin(0,-10,-10,-10),
                  plot.title = element_text(hjust = 0.5, size=12),
                  axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
                  axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
                  panel.background = element_rect(fill='white', colour='black'),
                  panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
        xlab(expression(paste("Distance, r (", mu, "m)"))) + ylab("-AIC") )
    
    # if(select_feature == 1){
    #     feature_info_combined = rbind(feature_info_combined, data.frame(feature_value=parameters(p)$r, 
    #                                                                     ens_location=ens_location, 
    #                                                                     sample_id=sample_id) )
    # }else if(select_feature == 2){
    #     feature_info_combined = rbind(feature_info_combined, data.frame(feature_value=parameters(p)$hc, 
    #                                                                     ens_location=ens_location, 
    #                                                                     sample_id=sample_id) )
    # }else if(select_feature == 3){
    #     feature_info_combined = rbind(feature_info_combined, data.frame(feature_value=parameters(p)$gamma, 
    #                                                                     ens_location=ens_location, 
    #                                                                     sample_id=sample_id) )
    # }
    
    feature_info_combined = rbind(feature_info_combined, data.frame(spp_r=parameters(p)$r, 
                                                                    spp_hc=parameters(p)$hc, 
                                                                    spp_gamma=parameters(p)$gamma, 
                                                                    ens_location=ens_location, 
                                                                    sample_id=sample_id) )
}

write.csv(feature_info_combined, paste(figure_folder, "SppFeaturesperSample.csv", sep = ""))

# g = g + geom_abline(slope=0, intercept = 0)+
#         geom_vline(xintercept = 0)+
#         theme(legend.position="top",  legend.title=element_blank(), legend.text=element_text(size=8),
#           legend.box.margin=margin(0,-10,-10,-10),
#           plot.title = element_text(hjust = 0.5, size=12),
#           axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
#           axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
#           panel.background = element_rect(fill='white', colour='black'),
#           panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
#         xlab(expression(paste("Interaction distance (r)"))) + ylab("L(r)") +
#         ggtitle("Besag's centered inhomogeneous L-function") 
# g


#### creating a common range of r-values
# r_vect_size = 100
# inhom_l_r_vect = seq(0, min(max_inhom_l_r), by=(min(max_inhom_l_r)/r_vect_size))
# 
# 
# inhom_l_all = data.frame(sample=branch_info_files, 
#                          group=c(rep(1, 11), rep(2, 9), rep(3, 16))) # 1: distal, 2: middle, 3: proximal
# 
# 
# #### computing inhomogeneous L function for the common range of r values (just computed above) for all the samples
# #### to figure out the most varying r. This r will be used in local computations.
# for(i in c(1:length(pp_list))){
#     u_branch.ppp = pp_list[[i]]
#     inhom_l = Linhom(u_branch.ppp, r=inhom_l_r_vect, correction = "border")
#     
#     inhom_l_all[i, 3:103] = inhom_l$border
#     
#     print(ggplot(inhom_l)+
#               geom_abline(slope=0, intercept = 0)+
#               geom_vline(xintercept = 0)+
#               
#               geom_line(aes(x=r, y=border-r, colour=i)) +
#               
#               theme(legend.position="top",  legend.title=element_blank(), legend.text=element_text(size=8),
#                     legend.box.margin=margin(0,-10,-10,-10),
#                     plot.title = element_text(hjust = 0.5, size=12),
#                     axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
#                     axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
#                     panel.background = element_rect(fill='white', colour='black'),
#                     panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
#               
#               xlab(expression(paste("Interaction distance (r) [cropped]"))) + ylab("L(r)") +
#               ggtitle("Besag's centered inhomogeneous L-function") )
# }
# 
# 
# f_val_vector = c()
# 
# for(i in c(3:103)){
#     anova_res = summary(aov(inhom_l_all[, i]~inhom_l_all$group))
#     f_val_vector[i-2] = anova_res[[1]]$`F value`[1]
# }
# 
# f_ratio = data.frame(r=inhom_l_r_vect, f_ratio=f_val_vector)
# f_ratio = f_ratio[order(f_ratio$f_ratio, decreasing = TRUE),]
# 
# #### top 3 r values with top variances
# top_r = f_ratio[1:3,]$r
# 
# top_r[1] # the r-value to be used

# dist_pair = as.matrix(dist(inhom_l_all[, 3:103]))
# 
# isSymmetric(dist_pair)
# diag(dist_pair)
# 
# #### multi-dimensional scaling to map entities into the Euclidean space of Linhom
# dd= cmdscale(dist_pair)
# 
# #### creating data frame with all data
# df=data.frame(x=dd[,1], y=dd[,2])
# sample_info=c(rep("Distal", 11), rep("Middle", 9), rep("Proximal", 16))
# sample_seq = c(1: length(branch_info_files))
# 
# #### including the biological information
# df = cbind(sample_seq, df, sample_info)
# 
# ggplot(data = df, aes(x = x, y = y, label=sample_seq)) +
#     geom_point(size = 3, shape=21, colour="grey40", aes(fill = sample_info)) +
#     
#     geom_text_repel(size=5, 
#                     segment.size=0.2, segment.color="grey40",
#                     arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "first"),
#                     max.overlaps = Inf) +
#     
#     #stat_ellipse(type="t", aes(colour=Nerve.Location.Avr)) +
#     
#     theme_bw()+
#     theme(legend.position="top", legend.text=element_text(size=20),
#           legend.box.margin=margin(-10,-10,-10,-10),
#           plot.title = element_text(hjust = 0.5, size=24),
#           legend.title = element_text(size=20),
#           plot.subtitle = element_text(hjust = 0.5, size=20),
#           axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
#           axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
#           panel.background = element_rect(fill='white', colour='black'),
#           panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
#     
#     xlab(expression(paste("Dimension 1"))) + ylab("Dimension 2") +
#     
#     labs(title = "Euclidean Space of Linhom")
# 

ggplot(feature_info_combined, aes(x = ens_location, y = log(feature_value), fill = ens_location)) +
    geom_boxplot(notch=FALSE, outlier.size = 1) +
    geom_quasirandom(cex=0.5, shape = 21, colour = "grey40", aes(fill=ens_location)) +
    
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="top", legend.text=element_text(size=16), legend.title = element_blank(),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=18),
          plot.subtitle = element_text(hjust = 0.5, size=16),
          axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
    
    xlab(expression(paste("ENS location"))) + ylab(paste("log(", spp_feature_list[select_feature], ")",sep = "") )+
    
    #labs(title = paste("Statistical comparison of Face ", face_feature_list[select_feature] ," from different part of ENS", sep="") )   
    labs(title = paste("Statistical comparison of SPP ", spp_feature_list[select_feature], sep="")) # the titles needs changing for different runs

#### statistical tests
#### transformed values
# feature_info_combined$log_f = log10(feature_info_combined$feature_value)
# feature_info_combined$sqrt_f = sqrt(feature_info_combined$feature_value)
# feature_info_combined$cbrt_f = (feature_info_combined$feature_value)^(1/3)

ggdensity(feature_info_combined[feature_info_combined$ens_location == "distal", ]$feature_value)
ggdensity(feature_info_combined[feature_info_combined$ens_location == "middle", ]$feature_value)
ggdensity(feature_info_combined[feature_info_combined$ens_location == "proximal", ]$feature_value)

ggqqplot(feature_info_combined[feature_info_combined$ens_location == "distal", ]$feature_value)
ggqqplot(feature_info_combined[feature_info_combined$ens_location == "middle", ]$feature_value)
ggqqplot(feature_info_combined[feature_info_combined$ens_location == "proximal", ]$feature_value)


shapiro.test(feature_info_combined[feature_info_combined$ens_location == "distal", ]$feature_value)
shapiro.test(feature_info_combined[feature_info_combined$ens_location == "middle", ]$feature_value)
shapiro.test(feature_info_combined[feature_info_combined$ens_location == "proximal", ]$feature_value)

#### conclusion: face feature not normally distributed, needs non-parametric test for 3 groups
#### Kruskal-Wallis test

feature_info_combined$ens_location = as.factor(feature_info_combined$ens_location)
levels(feature_info_combined$ens_location)

group_by(feature_info_combined, ens_location) %>%
    summarise(
        count = n(),
        mean = mean(feature_value, na.rm = TRUE),
        sd = sd(feature_value, na.rm = TRUE),
        median = median(feature_value, na.rm = TRUE),
        IQR = IQR(feature_value, na.rm = TRUE)
    )

kruskal.test(feature_value ~ ens_location, data = feature_info_combined)

pairwise.wilcox.test(feature_info_combined$feature_value, feature_info_combined$ens_location,
                     p.adjust.method = "BH")