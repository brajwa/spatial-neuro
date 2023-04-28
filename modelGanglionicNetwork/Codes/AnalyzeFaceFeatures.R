#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim", "ggbeeswarm", "ggpubr")

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

moment_folder = paste(parent, "Data/ENSMouse Moment Information (in um)/", sep="")
moment_files = list.files(moment_folder, recursive = TRUE, pattern = "\\.csv", full.names = TRUE)

figure_folder = paste(parent, "Outputs/ENSMouse/FaceFeature/", sep="")

face_feature_list = c("Area", "Perim.", "Ext.", "Disp.", "Elong.", "Eccentr.", "Orient.") 
select_feature = 7
cat("Face feature under consideration: ", face_feature_list[select_feature], "\n")

columns = c("feature_value","ens_location","sample_id") 
feature_info_combined = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(feature_info_combined) = columns

for (i in c(1:length(moment_files))) { #2,13,21
    ens_location = strsplit(moment_files[i], "/")[[1]][11]
    sample_id = strsplit(strsplit(moment_files[i], "/")[[1]][12], "\\.")[[1]][1]
    cat("\nLocation: ", ens_location, "\nSample Id: ", sample_id, "\n")
    
    moment_info = read.csv(moment_files[i])
    feature_info = moment_info[[face_feature_list[select_feature]]]
    
    face_count = length(feature_info)
    cat("Original face count: ", face_count, "\n")
    
    #check if NaN or Inf present
    feature_info = feature_info[!is.na(feature_info) & !is.infinite(feature_info)]
    face_count = length(feature_info)
    cat("Face count (NaN removed): ", face_count, "\n")
    
    feature_info_combined = rbind(feature_info_combined, data.frame(feature_value=feature_info, 
                                                                    ens_location=rep(ens_location, face_count), 
                                                                    sample_id=rep(sample_id, face_count)) )
    
}

#write.csv(feature_info_combined, paste(figure_folder, face_feature_list[select_feature], ".csv", sep = ""))

#make the face orientation angle between 0 and 180
#feature_info_combined$feature_value = (360 + feature_info_combined$feature_value) %% 180

svglite(paste(figure_folder, "All_", face_feature_list[select_feature], ".svg",  sep=""), width = 6, height = 4)
ggplot(feature_info_combined, aes(x = ens_location, y = log(feature_value), fill = ens_location)) +
    geom_boxplot(notch=TRUE, outlier.size = 1) +
    #geom_quasirandom(cex=0.5, shape = 21, colour = "grey40", aes(fill=ens_location)) +
    
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="top", legend.text=element_text(size=16), legend.title = element_blank(),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=18),
          plot.subtitle = element_text(hjust = 0.5, size=16),
          axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
    
    xlab(expression(paste("ENS location"))) + ylab(paste("log(", face_feature_list[select_feature], ")",sep = "") )+
    
    #labs(title = paste("Statistical comparison of Face ", face_feature_list[select_feature] ," from different part of ENS", sep="") )   
    labs(title = paste("Statistical comparison of Face ", face_feature_list[select_feature], sep="")) # the titles needs changing for different runs
dev.off()

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