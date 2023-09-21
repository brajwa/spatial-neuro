#### loading required libraries (there might be more libraries loaded than required)
load_lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
             "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
             "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
             "causaloptim", "RBGL", "svglite", "ggrepel", "devtools", "geosphere", "philentropy",
             "collections")

install_lib = load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)

#######################################################################################################################

Y <- rnoise(runif, square(1), max=100, dimyx=32)
Z <- Smooth(Y, sigma=0.07, normalise=TRUE, bleed=FALSE)
X <- rpoispp(Z)

plot(Y)

svglite(paste("D:/Fall 2023/Research/Prelim/figures/example pp/cox_intensity.svg", sep=""), width = 3, height = 1.8)
par(mar = c(0, 0, 0, 0))
plot(Z, main="", col=grey(seq(0, 1, length = 256)))
dev.off()

svglite(paste("D:/Fall 2023/Research/Prelim/figures/example pp/cox_pp.svg", sep=""), width = 1.8, height = 1.8)
par(mar = c(0, 0, 0, 0))
plot(X, main="", pch=21, cex=0.6, bg="black")
dev.off()


mu <- as.im(function(x,y){ exp(2 * x + 1) }, owin())
X <- rMatClust(10, 0.05, mu)
plot(X, main="", pch=21, cex=0.6, bg="black")

svglite(paste("D:/Fall 2023/Research/Prelim/figures/example pp/matern_pp.svg", sep=""), width = 1.8, height = 1.8)
par(mar = c(0, 0, 0, 0))
plot(X, main="", pch=21, cex=0.6, bg="black")
dev.off()

#######################################################################################################################

P <- rHardcore(100, 0.1)
svglite(paste("D:/Fall 2023/Research/Prelim/figures/example pp/hardcore_pp.svg", sep=""), width = 1.8, height = 1.8)
par(mar = c(0, 0, 0, 0))
plot(P, main="", cex=2.3, lwd=0.3)
plot(P, main="", pch=21, cex=0.6, bg="black", add=TRUE)
dev.off()

Q <- rStrauss(100, 0.5, 0.1)
svglite(paste("D:/Fall 2023/Research/Prelim/figures/example pp/strauss_pp.svg", sep=""), width = 1.8, height = 1.8)
par(mar = c(0, 0, 0, 0))
plot(Q, main="", cex=2.3, lwd=0.3)
plot(Q, main="", pch=21, cex=0.6, bg="black", add=TRUE)
dev.off()

#######################################################################################################################

A <- rpoispp(100)
svglite(paste("D:/Fall 2023/Research/Prelim/figures/example pp/pois_pp.svg", sep=""), width = 1.8, height = 1.8)
par(mar = c(0, 0, 0, 0))
plot(A, main="", pch=21, cex=0.6, bg="black")
dev.off()

i <- function(x,y){ 10 + 20*x + 30*y }
B <- rpoispp(i)
svglite(paste("D:/Fall 2023/Research/Prelim/figures/example pp/i_pois_pp.svg", sep=""), width = 1.8, height = 1.8)
par(mar = c(0, 0, 0, 0))
plot(B, main="", pch=21, cex=0.6, bg="black")
dev.off()