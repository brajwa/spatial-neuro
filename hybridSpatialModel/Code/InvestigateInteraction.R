#### loading required libraries (there might be more libraries loaded than required)
load.lib = c("this.path", "spatstat", "ggplot2", "tidyverse", "cowplot")

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


#### investigate the pelvic samples
pelvic_pp = pp_list[c(19:29)]
max_inhom_l_r = c()

for(axon_pp in pelvic_pp){
  summary(axon_pp)
  plot(axon_pp, pch=19, cex=0.3, bg="black")
  
  inhom_l = Linhom(axon_pp, correction = "border")
  max_inhom_l_r[length(max_inhom_l_r) + 1] = max(inhom_l$r)
  
  inhom_pcf = pcfinhom(axon_pp)
  
  par(mfrow = c(1,2))
  plot(inhom_l, .-r~r)
  plot(inhom_pcf)
  
  # ggobj = ggplot(data = inhom_l) + 
  #   geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) + 
  #   geom_line(aes(x=inhom_l$r, y=inhom_l$border-inhom_l$r), linewidth=1) +
  #   geom_vline(xintercept = min(nndist(axon_pp)), linewidth=1, lty=2, colour="blue") +
  #   theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(),
  #         #legend.box.margin=margin(-10,-10,-10,-10),
  #         plot.title = element_text(hjust = 0.5, size=10),
  #         plot.subtitle = element_text(hjust = 0.5, size=10),
  #         axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
  #         axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
  #         panel.background = element_rect(fill='white', colour='black'),
  #         panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) + 
  #   xlab(expression(paste("Distance, r (", mu, "m scaled)"))) + ylab("Besag's centered L-function")
  # plot(ggobj)
}

limit_r = min(max_inhom_l_r)

for(axon_pp in pelvic_pp){
  epsilon = 0.01
  init_r = min(nndist(axon_pp))

  ####figuring out R
  range_R = data.frame(r=seq(init_r, limit_r, by=0.001))
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
    geom_line(aes(x=r, y=gamma), linewidth=1) +
    
    geom_vline(xintercept = parameters(p)$r, linewidth=1, lty=2, colour="red") +
    geom_hline(yintercept = 1, linewidth=0.5, colour="black") +
    
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
  
  opt_r = as.numeric(p$fit$interaction$par[1])
  opt_hc = as.numeric(p$fit$interaction$par[2])
  opt_gamma = as.numeric(exp(p$fit$coef[4]))

  org_beta = intensity(axon_pp)
  org_win = axon_pp$window
  
  model_1 = list(cif="straush", par=list(beta=org_beta, gamma=opt_gamma, r=opt_r, hc=opt_hc), w=org_win)
  
  sim_1 = rmh(model=model_1)
  summary(sim_1)
  plot(sim_1, pch=19, cex=0.3, bg="black")
  
  d = density(axon_pp, sigma = 0.008)
  plot(d)
  
  Window(sim_1) = Window(d)
  d = d / max(d)
  d$v[is.na(d$v)] = 0
  sim_1_thinned = rthin(sim_1, d)
  Window(sim_1_thinned) = Window(axon_pp)
  plot(sim_1_thinned, pch=19, cex=0.3, bg="black")
  
  sim_inhom_l = Linhom(sim_1, correction = "border")
  ggobj = ggplot(data = sim_inhom_l) + 
    geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) + 
    geom_line(aes(x=sim_inhom_l$r, y=sim_inhom_l$border-sim_inhom_l$r), linewidth=1) +
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
  
  sim_inhom_l = Linhom(sim_1_thinned, correction = "border")
  ggobj = ggplot(data = sim_inhom_l) + 
    geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) + 
    geom_line(aes(x=sim_inhom_l$r, y=sim_inhom_l$border-sim_inhom_l$r), linewidth=1) +
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
  
}
