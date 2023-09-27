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
M <- rMatClust(10, 0.05, mu)
plot(M, main="", pch=21, cex=0.6, bg="black")

svglite(paste("D:/Fall 2023/Research/Prelim/figures/example pp/matern_pp.svg", sep=""), width = 1.8, height = 1.8)
par(mar = c(0, 0, 0, 0))
plot(M, main="", pch=21, cex=0.6, bg="black")
dev.off()

X_env_linhom = envelope(X, fun = "Linhom", correction = "isotropic", nsim = 39)
M_env_linhom = envelope(M, fun = "Linhom", correction = "isotropic", nsim = 39)

my_color_map = c("CSR"="darkgrey", "Cox"="cyan3", "Matern"="indianred1")
g1 = ggplot()+
    geom_ribbon(aes(x=X_env_linhom$r, ymin=X_env_linhom$lo-X_env_linhom$r, ymax=X_env_linhom$hi-X_env_linhom$r, fill="Cox"), alpha=0.2) +
    geom_ribbon(aes(x=M_env_linhom$r, ymin=M_env_linhom$lo-M_env_linhom$r, ymax=M_env_linhom$hi-M_env_linhom$r, fill="Matern"), alpha=0.2) +
    
    geom_line(aes(x=X_env_linhom$r, y=X_env_linhom$obs-X_env_linhom$r, colour="Cox"), show.legend = F) +
    geom_line(aes(x=M_env_linhom$r, y=M_env_linhom$obs-M_env_linhom$r, colour="Matern"), show.legend = F) +
    
    theme(legend.position="top",  legend.title=element_blank(), legend.text=element_text(size=8),
          legend.box.margin=margin(0,-10,-10,-10),
          legend.key.size = unit(0.5, "cm"),
          plot.title = element_text(hjust = 0.3, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=8),
          axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 8),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", size=0.25, linetype=2)) +
    scale_color_manual(values = my_color_map) + 
    scale_fill_manual(values = my_color_map) +
    xlab(expression(paste("Interaction distance (r)"))) + ylab("L(r)")
#labs(title="Besag's centered inhomogeneous L-function")

print(g1)

svglite(paste("D:/Fall 2023/Research/Prelim/figures/example pp/cluster_linhom.svg", sep=""), width = 2.5, height = 2)
par(mar = c(0, 0, 0, 0))  
print(g1)
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

P_env_linhom = envelope(P, fun = "Linhom", correction = "isotropic", nsim = 39)
Q_env_linhom = envelope(Q, fun = "Linhom", correction = "isotropic", nsim = 39)

my_color_map = c("CSR"="darkgrey", "Hardcore"="cyan3", "Strauss"="indianred1")
g2 = ggplot()+
    geom_ribbon(aes(x=P_env_linhom$r, ymin=P_env_linhom$lo-P_env_linhom$r, ymax=P_env_linhom$hi-P_env_linhom$r, fill="Hardcore"), alpha=0.2) +
    geom_ribbon(aes(x=Q_env_linhom$r, ymin=Q_env_linhom$lo-Q_env_linhom$r, ymax=Q_env_linhom$hi-Q_env_linhom$r, fill="Strauss"), alpha=0.2) +
    
    geom_line(aes(x=P_env_linhom$r, y=P_env_linhom$obs-P_env_linhom$r, colour="Hardcore"), show.legend = F) +
    geom_line(aes(x=Q_env_linhom$r, y=Q_env_linhom$obs-Q_env_linhom$r, colour="Strauss"), show.legend = F) +
    
    theme(legend.position="top",  legend.title=element_blank(), legend.text=element_text(size=8),
          legend.box.margin=margin(0,-10,-10,-10),
          legend.key.size = unit(0.5, "cm"),
          plot.title = element_text(hjust = 0.3, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=8),
          axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 8),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", size=0.25, linetype=2)) +
    scale_color_manual(values = my_color_map) + 
    scale_fill_manual(values = my_color_map) +
    xlab(expression(paste("Interaction distance (r)"))) + ylab("L(r)")
#labs(title="Besag's centered inhomogeneous L-function")

print(g2)

svglite(paste("D:/Fall 2023/Research/Prelim/figures/example pp/inhibition_linhom.svg", sep=""), width = 2.5, height = 2)
par(mar = c(0, 0, 0, 0))  
print(g2)
dev.off()


#######################################################################################################################

A <- rpoispp(100)
svglite(paste("D:/Fall 2023/Research/Prelim/figures/example pp/pois_pp.svg", sep=""), width = 1.5, height = 1.5)
par(mar = c(0, 0, 0, 0))
plot(A, main="", pch=21, cex=0.6, bg="black")
dev.off()

i <- function(x,y){ 10 + 20*x + 30*y }
B <- rpoispp(i)
svglite(paste("D:/Fall 2023/Research/Prelim/figures/example pp/i_pois_pp.svg", sep=""), width = 1.5, height = 1.5)
par(mar = c(0, 0, 0, 0))
plot(B, main="", pch=21, cex=0.6, bg="black")
dev.off()

A_env_linhom = envelope(A, fun = "Linhom", correction = "isotropic", nsim = 39)
B_env_linhom = envelope(B, fun = "Linhom", correction = "isotropic", nsim = 39)


my_color_map = c("CSR"="darkgrey", "Hom"="cyan3", "Inhom"="indianred1")
g3 = ggplot()+
    geom_ribbon(aes(x=A_env_linhom$r, ymin=A_env_linhom$lo-A_env_linhom$r, ymax=A_env_linhom$hi-A_env_linhom$r, fill="Hom"), alpha=0.2) +
    geom_ribbon(aes(x=B_env_linhom$r, ymin=B_env_linhom$lo-B_env_linhom$r, ymax=B_env_linhom$hi-B_env_linhom$r, fill="Inhom"), alpha=0.2) +
    
    geom_line(aes(x=A_env_linhom$r, y=A_env_linhom$obs-A_env_linhom$r, colour="Hom"), show.legend = F) +
    geom_line(aes(x=B_env_linhom$r, y=B_env_linhom$obs-B_env_linhom$r, colour="Inhom"), show.legend = F) +
    
    theme(legend.position="top",  legend.title=element_blank(), legend.text=element_text(size=8),
          legend.box.margin=margin(0,-10,-10,-10),
          legend.key.size = unit(0.5, "cm"),
          plot.title = element_text(hjust = 0.3, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=8),
          axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 8),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", size=0.25, linetype=2)) +
    scale_color_manual(values = my_color_map) + 
    scale_fill_manual(values = my_color_map) +
    xlab(expression(paste("Interaction distance (r)"))) + ylab("L(r)")
    #labs(title="Besag's centered inhomogeneous L-function")

print(g3)

svglite(paste("D:/Fall 2023/Research/Prelim/figures/example pp/pois_linhom.svg", sep=""), width = 2.5, height = 2)
par(mar = c(0, 0, 0, 0))  
print(g3)
dev.off()
