#### loading required libraries (there might be more libraries loaded than required)
load.lib = c("this.path", "spatstat", "ggplot2")

install.lib = load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib, dependencies=TRUE)
sapply(load.lib, require, character=TRUE)


#### extracting parent directory information for accessing input and output location
dir = this.dir()
folder = strsplit(dir, "/")
folder = folder[[1]][length(folder[[1]])]
parent = strsplit(dir, folder)


#### input location
folder_path = paste(parent, "Data/Peripheral Axon/Inputs/", sep = "")
data_files = list.files(path=folder_path, pattern = ".csv", full.names = FALSE)
extension = ".csv"


pp_list = list()

for (file_name in data_files) {
    file_name = strsplit(file_name, ".csv")[[1]]
    print(file_name)
    
    axon_locations = unique(read.csv(paste(folder_path, file_name, extension, sep="")))
    retrieved_contour = readRDS(paste(folder_path, file_name, ".rds", sep=""))
    
    axon_pp = ppp(x=axon_locations$X, y=axon_locations$Y, checkdup=F, window = retrieved_contour)
    axon_pp = rescale.ppp(axon_pp, s=320153.4)  #pre-computed maximum y-range of all the fascicles (ScalingFactorComputation.R)
    axon_pp = shift.ppp(axon_pp, origin = "centroid")
    
    pp_list[[length(pp_list) + 1]] = axon_pp
}

axon_pp = pp_list[[15]]
summary(axon_pp)
plot(axon_pp, pch=19, cex=0.6, bg="black")

epsilon = 0.01
init_r = min(nndist(axon_pp)) + epsilon

inhom_l = Linhom(axon_pp, correction = "border")
limit_r = max(inhom_l$r)

inhom_l = Linhom(axon_pp, correction = "border")
ggobj = ggplot(data = inhom_l) + 
  geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) + 
  geom_line(aes(x=inhom_l$r, y=inhom_l$border-inhom_l$r), linewidth=1) +
  theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(),
        #legend.box.margin=margin(-10,-10,-10,-10),
        plot.title = element_text(hjust = 0.5, size=10),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) + 
  xlab(expression(paste("Distance, r (", mu, "m scaled)"))) + ylab("Besag's centered L-function")
plot(ggobj)

####figuring out R
range_R = data.frame(r=seq(init_r, limit_r, by=0.002))
p = profilepl(s=range_R, f=StraussHard, axon_pp~x+y, aic=TRUE, rbord = init_r)
p

# -AIC and Gamma values
p_info = data.frame(r=p$param$r, aic=p$prof, gamma=exp(p$allcoef$Interaction))

ggobj = ggplot(data = p_info) + 
    geom_point(aes(x=r, y=aic), size=1.5) +
    geom_line(aes(x=r, y=aic), linewidth=1) +
    geom_vline(xintercept = parameters(p)$r, linewidth=1, lty=2, colour="red") +
    theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=10),
          axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) + 
    xlab(expression(paste("Distance, r (", mu, "m)"))) + ylab("-AIC")
plot(ggobj)

ggobj = ggplot(data = p_info) + 
    geom_point(aes(x=r, y=gamma), size=1.5) +
    geom_line(aes(x=r, y=gamma), size=1) +
    
    geom_vline(xintercept = parameters(p)$r, size=1, lty=2, colour="red") +
    geom_hline(yintercept = 1, size=0.5, colour="black") +
    
    theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=10),
          axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
    xlab(expression(paste("Distance, r (", mu, "m)"))) + ylab(expression(paste("Interaction parameter, ", gamma)))
plot(ggobj)