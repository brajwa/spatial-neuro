#### loading required libraries (there might be more libraries loaded than required)
load.lib = c("this.path", "spatstat", "ggplot2", "tidyverse", "cowplot", "svglite", "ggpmisc", "dplyr")

install.lib = load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib, dependencies=TRUE)
sapply(load.lib, require, character=TRUE)


#### computes and returns the x-coordinates at which the given function crosses the horizontal axis
#### parameter: points is a dataframe with two columns, x and y, defining the function points
computeCrossOverPoints <- function(points){
  points = points[is.finite(rowSums(points)),]
  x_intercepts = c()
  i = 1
  j = 2
  while(j <= length(points$x)){
    if( (points$y[i] <= 0 & points$y[j] >= 0) | (points$y[i] >= 0 & points$y[j] <= 0) ){
      x1 = points$x[i]
      x2 = points$x[j]
      
      y1 = points$y[i]
      y2 = points$y[j]
      
      x_cross = x1 + ((x2-x1)*y1)/(y1-y2)
      x_intercepts = c(x_intercepts, x_cross)
    }
    
    i = i + 1
    j = j + 1
  }
  return(x_intercepts)
}


#### computes and returns the interaction distances for the given point pattern that reflect significant
#### interaction based on its pair correlation function
detectLocalPeaks <- function(axon_pp, pp_id, peak_threshold, parent_path){
  inhom_pcf = pcfinhom(axon_pp)
  
  #### translate the pcf to y=0, to detect the crossover points
  translated_pcf = data.frame(x=inhom_pcf$r, y=inhom_pcf$trans-1)
  x_cross = computeCrossOverPoints(translated_pcf)
  x_cross = c(0, x_cross, max(inhom_pcf$r)) # the cross over points are now range endpoints
  
  #### takes absolute value of the pcf to detect the peaks within the computed ranges
  modified_pcf = data.frame(r=inhom_pcf$r, f=abs(inhom_pcf$trans-1))
  
  intr_info = data.frame(matrix(nrow=0, ncol = 4))
  i = 1
  j = 2
  while(j <= length(x_cross)){
    begin = x_cross[i]
    end = x_cross[j]
    
    filtered_pcf = modified_pcf[(modified_pcf$r>=begin & modified_pcf$r<=end), ]
    filtered_pcf = filtered_pcf[is.finite(rowSums(filtered_pcf)),]
    
    if(max(filtered_pcf$f) >= peak_threshold){
      intr_dist = filtered_pcf$r[which.max(filtered_pcf$f)]
      
      filtered_pcf_2 = translated_pcf[(translated_pcf$x>=begin & translated_pcf$x<=end), ]
      filtered_pcf_2 = filtered_pcf_2[is.finite(rowSums(filtered_pcf_2)),]
      if(sum(filtered_pcf_2$y) < 0){
        intr_type = "inhibition"
      }else{
        intr_type = "attraction"
      }
    }else{
      intr_dist = NA
      intr_type = "random"
    }
    
    intr_info = rbind(intr_info, c(begin, end, intr_dist, intr_type))
    
    i = i + 1
    j = j + 1
  }
  colnames(intr_info) = c("begin", "end", "intr_dist", "intr_type")
  
  ggobj = ggplot(data = inhom_pcf) + 
    geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) + 
    geom_hline(aes(yintercept=1)) +  
    geom_line(aes(x=r, y=trans), linewidth=1) +
    geom_vline(xintercept=as.numeric(intr_info$intr_dist[!is.na(intr_info$intr_dist)]), lty=2, colour="blue")
    theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=10),
          axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) + 
    xlab(expression(paste("Distance, r (", mu, "m scaled)"))) + ylab("Pair Correlation Function")
  plot(ggobj)
  # svglite(paste(parent_path, "Output/Peripheral Axon/Interaction Information/pelvic_pp_pcf_intr_", pp_id, ".svg", sep=""), width = 4.5, height = 3)
  # par(mar = c(0, 0, 0, 0))
  # plot(ggobj)
  # dev.off()
  
  return(intr_info)
  
  # p_r = modified_pcf$r[ggpmisc::find_peaks(modified_pcf$f)]
  # ggplot(data = modified_pcf, aes(x = r, y = f)) + 
  #   geom_line() + 
  #   stat_peaks(col = "red") +
  #   geom_vline(xintercept = p_r, linewidth=0.3, lty=2, colour="blue")
  
}


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
  
  # svglite(paste(parent, "Output/Peripheral Axon/PP/", file_name, ".svg", sep=""), width = 3, height = 3)
  # par(mar = c(0, 0, 0, 0))
  # plot(axon_pp, pch=19, cex=0.3, bg="black", main="")
  # dev.off()
}


#### investigate the pelvic samples
pelvic_pp = pp_list[c(19:29)]
max_inhom_l_r = c()

i = 1
for(axon_pp in pelvic_pp){
  summary(axon_pp)
  plot(axon_pp, pch=19, cex=0.3, bg="black")
  
  inhom_l = Linhom(axon_pp)
  max_inhom_l_r[length(max_inhom_l_r) + 1] = max(inhom_l$r)
  ggobj = ggplot(data = inhom_l) + 
    geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) + 
    geom_line(aes(x=inhom_l$r, y=inhom_l$border-inhom_l$r), linewidth=1) +
    geom_vline(xintercept = min(nndist(axon_pp)), linewidth=0.3, lty=2, colour="blue") +
    theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=10),
          axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) + 
    xlab(expression(paste("Distance, r (", mu, "m scaled)"))) + ylab("Besag's centered L-function")
  #plot(ggobj)
  # svglite(paste(parent, "Output/Peripheral Axon/Inhom-L/pelvic_pp_inhom_l_", i, ".svg", sep=""), width = 4.5, height = 3)
  # par(mar = c(0, 0, 0, 0))
  # plot(ggobj)
  # dev.off()
  
  inhom_pcf = pcfinhom(axon_pp)
  ggobj = ggplot(data = inhom_pcf) + 
    geom_hline(aes(yintercept=1)) + geom_vline(aes(xintercept=0)) + 
    geom_line(aes(x=inhom_pcf$r, y=inhom_pcf$trans), linewidth=1) +
    geom_vline(xintercept = min(nndist(axon_pp)), linewidth=0.3, lty=2, colour="blue") +
    theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=10),
          axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) + 
    xlab(expression(paste("Distance, r (", mu, "m scaled)"))) + ylab("Pair Correlation Function")
  #plot(ggobj)
  # svglite(paste(parent, "Output/Peripheral Axon/PCF/pelvic_pp_inhom_pcf_", i, ".svg", sep=""), width = 4.5, height = 3)
  # par(mar = c(0, 0, 0, 0))
  # plot(ggobj)
  # dev.off()
  
  #### compute envelope of pcf to get the CSR band
  #### this will act as peak threshold 
  env_pcf = envelope(axon_pp, "pcf")
  env_pcf = env_pcf[is.finite(rowSums(env_pcf)),]
  p_t = mean(env_pcf$hi-env_pcf$lo)
  cat(p_t, "\n")
  print(detectLocalPeaks(axon_pp, i, peak_threshold = p_t, parent))
  
  i = i + 1
}

limit_r = min(max_inhom_l_r)

i = 1
for(axon_pp in pelvic_pp){
  epsilon = 0.01
  init_r = min(nndist(axon_pp))

  ####figuring out R
  range_R = data.frame(r=seq(init_r, limit_r, by=0.001))
  p = profilepl(s=range_R, f=AreaInter, axon_pp~polynom(x, y, 2), aic=TRUE, rbord = init_r)
  p
  
  # -AIC and Gamma values
  p_info = data.frame(r=p$param$r, aic=p$prof, gamma=exp(p$allcoef$Interaction))
  
  ggobj = ggplot(data = p_info) + 
    geom_point(aes(x=r, y=aic), size=1.5) +
    geom_line(aes(x=r, y=aic), linewidth=1) +
    geom_vline(xintercept = parameters(p)$r, linewidth=0.5, lty=2, colour="#69b3a2") +
    geom_label(label=paste("Opt r: ", round(parameters(p)$r, 4), "\nOpt intr param: ", round(parameters(p)$eta, 2), sep = ""), x=parameters(p)$r, y=mean(p_info$aic),label.padding = unit(0.55, "lines"), # Rectangle size around label
               label.size = 0.25,color = "black", fill="#69b3a2", size=2)+
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
  svglite(paste(parent, "Output/Peripheral Axon/AIC/pelvic_pp_aic_", i, ".svg", sep=""), width = 4, height = 3)
  par(mar = c(0, 0, 0, 0))
  plot(ggobj)
  dev.off()
  
  # opt_r = as.numeric(p$fit$interaction$par[1])
  # opt_hc = as.numeric(p$fit$interaction$par[2])
  # opt_gamma = as.numeric(exp(p$fit$coef[4]))
  # org_beta = intensity(axon_pp)
  # org_win = axon_pp$window
  # model_1 = list(cif="straush", par=list(beta=org_beta, gamma=opt_gamma, r=opt_r, hc=opt_hc), w=org_win)
  # sim_1 = rmh(model=model_1)
  
  fitM = ppm(axon_pp~polynom(x, y, 2), AreaInter(r=parameters(p)$r))
  # fitM = ppm(axon_pp~polynom(x, y, 2), Hybrid(Hardcore(hc=init_r), AreaInter(r=parameters(p)$r)) )
  
  sim_1 = simulate(fitM)
  sim_1 = sim_1$`Simulation 1`
  summary(sim_1)
  plot(sim_1, pch=19, cex=0.3, bg="black")
  
  svglite(paste(parent, "Output/Peripheral Axon/Simulated PP/pelvic_pp_", i, "_sim2.svg", sep=""), width = 3, height = 3)
  par(mar = c(0, 0, 0, 0))
  plot(sim_1, pch=19, cex=0.3, bg="black", main="")
  dev.off()
  
  # d = density(axon_pp, sigma = 0.008)
  # plot(d)
  # 
  # Window(sim_1) = Window(d)
  # d = d / max(d)
  # d$v[is.na(d$v)] = 0
  # sim_1_thinned = rthin(sim_1, d)
  # Window(sim_1_thinned) = Window(axon_pp)
  # plot(sim_1_thinned, pch=19, cex=0.3, bg="black")
  
  sim_inhom_l = Linhom(sim_1)
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
  svglite(paste(parent, "Output/Peripheral Axon/Simulated PP/pelvic_pp_", i, "_sim2_inhom_l.svg", sep=""), width = 4.5, height = 3)
  par(mar = c(0, 0, 0, 0))
  plot(ggobj)
  dev.off()
  
  sim_inhom_pcf = pcfinhom(sim_1)
  ggobj = ggplot(data = sim_inhom_pcf) + 
    geom_hline(aes(yintercept=1)) + geom_vline(aes(xintercept=0)) + 
    geom_line(aes(x=sim_inhom_pcf$r, y=sim_inhom_pcf$trans), linewidth=1) +
    geom_vline(xintercept = min(nndist(axon_pp)), linewidth=0.3, lty=2, colour="blue") +
    theme(legend.position="top", legend.text=element_text(size=8), legend.title = element_blank(),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=10),
          axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) + 
    xlab(expression(paste("Distance, r (", mu, "m scaled)"))) + ylab("Pair Correlation Function")
  plot(ggobj)
  svglite(paste(parent, "Output/Peripheral Axon/Simulated PP/pelvic_pp_", i, "_sim2_inhom_pcf.svg", sep=""), width = 4.5, height = 3)
  par(mar = c(0, 0, 0, 0))
  plot(ggobj)
  dev.off()
  
  i = i + 1
  
}
