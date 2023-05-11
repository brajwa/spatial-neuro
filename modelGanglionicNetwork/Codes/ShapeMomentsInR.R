#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim", "RBGL", "svglite", "ggrepel", "devtools")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)

devtools::install_github("swarm-lab/Rvision")

watershed_line_file = "C:/Users/sanja/Desktop/MAX_File_85_01-31-2019_CYM_1_4_Calb2_3B_YFP_DS_GFP-g_Hu-b.lif - TileScan_001_Merging001_ProjMax001_AdjustClr001-watershed-lines.jpg"

w_lines = Rvision::image((watershed_line_file))
plot(w_lines)

w_lines_gray = Rvision::changeColorSpace(w_lines, "GRAY")
w_lines_bin = w_lines_gray < 200
plot(w_lines_bin)

face_contours = Rvision::findContours(w_lines_bin)
num_of_faces = max(face_contours$contours[, 1]) + 1 #0-indexed

columns = c("Area", "Ext.", "Disp.", "Elong.", "Eccentr.", "Orient.") 
face_features = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(face_features) = columns

for(f in c(0: (num_of_faces-1))){
    cat("face id: ", f, "\n")
    
    f_contour = face_contours$contours[face_contours$contours[, 1] == f, 2:3]
    lines(f_contour, col="red", type="l", lwd=2)
    
    area = Rvision::contourArea(f_contour[,1], f_contour[,2])
    
    moments = Rvision::moments(f_contour)
    
    #### rotational invariants
    phi1 = moments$value[moments$moment == "nu02"] + moments$value[moments$moment == "nu20"]
    phi2 = ((moments$value[moments$moment == "nu02"] - moments$value[moments$moment == "nu20"]) * (moments$value[moments$moment == "nu02"] - moments$value[moments$moment == "nu20"])) 
        + (4 * moments$value[moments$moment == "nu11"] * moments$value[moments$moment == "nu11"])
    lambda1 = 2 * pi * (phi1 + sqrt(phi2))
    lambda2 = 2 * pi * (phi1 - sqrt(phi2))
        
    #### ext, disp, elong
    ext = log2(lambda1)
    disp = log2(sqrt(lambda1 * lambda2))
    elong = log2(sqrt(lambda1 / lambda2))
    
    #### orient, eccentr
    orient = 0.5 * atan2((2 * moments$value[moments$moment == "m11"]) , (moments$value[moments$moment == "m20"] - moments$value[moments$moment == "m02"]))
    orient = orient * 180 / pi
    eccentr = (((moments$value[moments$moment == "m02"] - moments$value[moments$moment == "m20"]) * (moments$value[moments$moment == "m02"] - moments$value[moments$moment == "m20"])) 
               + (4 * moments$value[moments$moment == "m11"] * moments$value[moments$moment == "m11"])) / moments$value[moments$moment == "m00"]
    
    face_features = rbind(face_features, data.frame(area, ext, disp, elong, eccentr, orient))
}



