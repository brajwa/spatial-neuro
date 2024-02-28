#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim", "RBGL", "svglite", "qdap", "lme4", "emmeans")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)


#### load previously computed and stored features
face_feature_watershed = read.csv("C:/Users/sanja/Documents/GitHub/spatial-neuro/modelGanglionicNetwork/Outputs/ENSMouse/FaceFeature/Facefeatures.csv")
face_feature_spn = read.csv("C:/Users/sanja/Documents/GitHub/spatial-neuro/modelGanglionicNetwork/Outputs/ENSMouse/FaceFeature/Facefeatures_3.csv")

#face_feature_watershed = na.omit(face_feature_watershed)
face_feature_spn = na.omit(face_feature_spn)

#watershed_feature_index = 2
spn_feature_index = 11

#transformed_watershed_feature = asinh(face_feature_watershed[, watershed_feature_index])
transformed_spn_feature = (face_feature_spn[, spn_feature_index])

#watershed_animal = sapply(face_feature_watershed$sample_id, function(x) genXtract(x, "sub-", "-sam"))
#spn_animal = sapply(face_feature_spn$sample_id, function(x) genXtract(x, "sub-", "-sam"))

#watershed_animal = face_feature_watershed$sample_id
spn_animal = face_feature_spn$sample_id

#subset_watershed = data.frame(transformed_watershed_feature, face_feature_watershed$ens_location, watershed_animal)
subset_spn = data.frame(transformed_spn_feature, face_feature_spn$ens_location, spn_animal)

#subset_watershed = subset_watershed[is.finite(subset_watershed$transformed_watershed_feature), ]
subset_spn = subset_spn[is.finite(subset_spn$transformed_spn_feature), ]

# watershed_model = lmer(transformed_watershed_feature ~ face_feature_watershed.ens_location + (1 | watershed_animal), data = subset_watershed)
# emmeans_watershed = emmeans(watershed_model, ~ face_feature_watershed.ens_location)
# 
# svglite(paste("C:/Users/sanja/Documents/GitHub/spatial-neuro/modelGanglionicNetwork/Outputs/ENSMouse/FaceFeature/", 
#               "ParamLinModel_Watershed_Area.svg",  sep=""), width = 6, height = 4)
# plot(emmeans_watershed)
# dev.off()

spn_model = lmer(transformed_spn_feature ~ face_feature_spn.ens_location + (1 | spn_animal), data = subset_spn)
emmeans_spn = emmeans(spn_model, ~ face_feature_spn.ens_location)

# svglite(paste("C:/Users/sanja/Documents/GitHub/spatial-neuro/modelGanglionicNetwork/Outputs/ENSMouse/FaceFeature/", 
#               "ParamLinModel_SPN_Area_CF.svg",  sep=""), width = 6, height = 4)
plot(emmeans_spn)
#dev.off()

