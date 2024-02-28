####A normalizing function
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

detectNegAICPeaks <-function(p_info, folder_path, file_name, extension){
  doc = read_pptx()
  
  peaks = detect_localmaxima(p_info$aic, w=5, boundary = "discard") # tune w here
  
  ggobj = ggplot() + 
    geom_line(aes(x=p_info$r, y=p_info$aic), size=1) +
    geom_point(aes(x=p_info$r[peaks], y=p_info$aic[peaks]), size=2, col="red") +
    
    theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
          legend.text=element_text(size=30),
          axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
    xlab(expression(paste("Index"))) + ylab("-AIC")
  
  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())

  tentative_w = length(p_info$r) %/% 2
  if(tentative_w %% 2 == 0) {
    tentative_w = tentative_w + 1
  }
  
  score = score_type1(p_info$aic, w=tentative_w) # tune w and score type here
  plot(seq(length(score)), score, type = "l", col = "red")

  true_peaks = score > quantile(score, probs=c(0.95)) & peaks # tune score threshold here
  #true_peaks = peaks

  ggobj = ggplot() + 
    geom_line(aes(x=p_info$r, y=p_info$aic), size=1) +
    geom_point(aes(x=p_info$r[peaks], y=p_info$aic[peaks]), size=2, col="red") +
    geom_point(aes(x=p_info$r[true_peaks], y=p_info$aic[true_peaks]), size=2, col="blue") +
    
    theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
          legend.text=element_text(size=30),
          axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
    xlab(expression(paste("Index"))) + ylab("-AIC")
  
  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
  
  print(doc, target = paste(folder_path, file_name, "_AIC_Local_Maxima.pptx", sep=""))
  
  true_peaks_r = p_info$r[true_peaks]
  
  return(list(true_peaks, true_peaks_r))
}

generateHybridModel <- function(folder_path, file_name, extension, fascicle_contour){
  doc = read_pptx()
  
  ####Reading data from file; using unique() to get rid of duplicate information; creating a new data frame to avoid nonsequential 
  ####data caused by the unique() function
  axon_locations = unique(read.xlsx(paste(folder_path, file_name, extension, sep="")))
  axon_locations = data.frame(Name=axon_locations$Name, CentroidX=axon_locations$CentroidX, CentroidY=axon_locations$CentroidY,
                              CentroidZ=axon_locations$CentroidZ, ShapeAdjustedEllipse_Fiber= axon_locations$ShapeAdjustedEllipse_Fiber)
  
  axon_pp_both = ppp(x = axon_locations$CentroidX, y = axon_locations$CentroidY, window = fascicle_contour, marks = axon_locations$Name)
  plot(axon_pp_both, main=file_name, cex=0.8, pch=21, bg=c("black", "cyan"))
  
  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(code = plot(axon_pp_both, main=file_name, cex=0.8, pch=21, bg=c("black", "cyan"))), 
                location = ph_location_fullsize())
  
  ###################################
  # filtering only Unmyelinated axons
  ###################################
  
  axon_locations = axon_locations[axon_locations$Name == 'Unmyelinated Axon',]
  axon_marks = data.frame(Name=axon_locations$Name, SAE=axon_locations$ShapeAdjustedEllipse_Fiber)
  
  ####bin the SAE values into categories to be included as marks
  SAE_binned = bin(axon_locations$ShapeAdjustedEllipse_Fiber, nbins=3, labels= c("small", "med", "large"), method="length")
  
  ####new marked point pattern
  axon_pp = ppp(x = axon_locations$CentroidX, y = axon_locations$CentroidY, window = fascicle_contour, marks = SAE_binned)
  plot(axon_pp, main=file_name, cex=0.7)
  
  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(code = plot(axon_pp, main=file_name, cex=0.7)), 
                location = ph_location_fullsize())
  
  ltest = envelope(axon_pp, Lest, nsim=39)
  ggobj = ggplot(data = ltest, aes(r)) + 
    
    geom_ribbon(aes(ymin=lo-r, ymax=hi-r), alpha = 0.3)+
    geom_line(aes(y=lo-r, colour="CSR lower bound"), size=2) +
    geom_line(aes(y=hi-r, colour="CSR upper bound"), size=2) +
    #geom_point(aes(y=theo, colour="csr")) + 
    
    geom_line(aes(y=obs-r, colour="observed"), size=2) + 
    #geom_point(aes(y=obs, colour="observed")) +  
    
    theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
          legend.text=element_text(size=30),
          axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
    xlab(expression(paste("Distance, r (", mu, "m)"))) + ylab("L-Function")
  
  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
  
  mean(axon_locations$ShapeAdjustedEllipse_Fiber)
  min(nndist(axon_pp))
  max(pairdist(axon_pp))
  
  epsilon = 0.01
  init_r = min(nndist(axon_pp)) + epsilon
  
  ####figuring out R
  range_R = data.frame(r=seq(init_r, max(ltest$r), by=0.05))
  p = profilepl(s=range_R, f=StraussHard, axon_pp~x+y+marks, aic=TRUE, rbord = init_r)
  p
  
  # -AIC and Gamma values
  p_info = data.frame(r=p$param$r, aic=p$prof, gamma=exp(p$allcoef$Interaction))
  
  ggobj = ggplot(data = p_info) + 
    geom_point(aes(x=r, y=aic), size=1.5) +
    geom_line(aes(x=r, y=aic), size=1) +
    
    geom_vline(xintercept = parameters(p)$r, size=1, lty=2, colour="red") +
    
    theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
          legend.text=element_text(size=30),
          axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
    xlab(expression(paste("Distance, r (", mu, "m)"))) + ylab("-AIC")
  
  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
  
  ggobj = ggplot(data = p_info) + 
    geom_point(aes(x=r, y=gamma), size=1.5) +
    geom_line(aes(x=r, y=gamma), size=1) +
    
    geom_vline(xintercept = parameters(p)$r, size=1, lty=2, colour="red") +
    geom_hline(yintercept = 1, size=0.5, colour="black") +
    
    theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
          legend.text=element_text(size=30),
          axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
    xlab(expression(paste("Distance, r (", mu, "m)"))) + ylab(expression(paste("Interaction parameter, ", gamma)))
  
  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
  
  print(doc, target = paste(folder_path, file_name, "_Gamma_Analysis.pptx", sep=""))
  
  ####residual of the initial model
  axon_model = as.ppm(p)
  doc = read_pptx()

  residual_K_init_model = Kres(axon_model, correction="best")
  ggobj = ggplot(data = residual_K_init_model, aes(r)) +

    geom_ribbon(aes(ymin=ilo, ymax=ihi), alpha = 0.3)+

    geom_line(aes(y=theo), size=1, linetype="dashed") +

    geom_line(aes(y=ires), size=1.2) +

    theme_bw()+
    theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(),
          legend.text=element_text(size=30),
          axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30),
          panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab(expression(paste("r"))) + ylab("R K(r)")

  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())

  residual_G_init_model = Gres(axon_model, correction="best")
  ggobj = ggplot(data = residual_G_init_model, aes(r)) +

    geom_ribbon(aes(ymin=lo, ymax=hi), alpha = 0.3)+

    geom_line(aes(y=theo), size=1, linetype="dashed") +

    geom_line(aes(y=hres), size=1.2) +

    theme_bw()+
    theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(),
          legend.text=element_text(size=30),
          axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30),
          panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab(expression(paste("r"))) + ylab("R G(r)")

  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())

  print(doc, target = paste(folder_path, file_name, "_Res_Sum_Functs.pptx", sep=""))
  # 
  # #--------------------------------------------------------------------------------------------------------------------------
  # ####hybrid model construction manually
  # manual_hybrid_interaction = Hybrid(StraussHard(r=parameters(p)$r), Strauss(r=0.5))
  # manual_hybrid_axon_model = ppm(axon_pp ~ x+y+marks, interaction = manual_hybrid_interaction)
  # 
  # manual_hybrid_axon_simulated = simulate(manual_hybrid_axon_model, nsim = 1)
  # plot(manual_hybrid_axon_simulated$`Simulation 1`, main=paste("Manual Hybrid Simulated\n", file_name, sep = " "), cex=0.7)
  # 
  # doc = add_slide(doc, "Blank", "Office Theme")
  # doc = ph_with(doc, dml(code = plot(manual_hybrid_axon_simulated$`Simulation 1`, main=paste("Manual Simulated Hybrid", file_name, sep = " "), cex=0.7)), 
  #               location = ph_location_fullsize())
  # 
  # manual_simulated_ltest_hybrid = envelope(manual_hybrid_axon_simulated$`Simulation 1`, Lest, nsim=39)
  # 
  # ggobj = ggplot(data = manual_simulated_ltest_hybrid, aes(r)) + 
  #   
  #   geom_ribbon(aes(ymin=lo-r, ymax=hi-r), alpha = 0.3)+
  #   geom_line(aes(y=lo-r, colour="CSR lower bound"), size=2) +
  #   geom_line(aes(y=hi-r, colour="CSR upper bound"), size=2) +
  #   #geom_point(aes(y=theo, colour="csr")) + 
  #   
  #   geom_line(aes(y=obs-r, colour="manual simulated hybrid"), size=2) + 
  #   #geom_point(aes(y=obs, colour="observed")) +  
  #   
  #   theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
  #         legend.text=element_text(size=30),
  #         axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
  #         axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
  #   xlab(expression(paste("Distance, r (", mu, "m)"))) + ylab("L-Function")
  # 
  # doc = add_slide(doc, "Blank", "Office Theme")
  # doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
  # 
  # residual_K_model_2 = Kres(manual_hybrid_axon_model, correction="best")
  # ggobj = ggplot(data = residual_K_model_2, aes(r)) + 
  #   
  #   geom_ribbon(aes(ymin=ilo, ymax=ihi), alpha = 0.3)+
  #   
  #   geom_line(aes(y=theo), size=1, linetype="dashed") + 
  #   
  #   geom_line(aes(y=ires), size=1.2) + 
  #   
  #   theme_bw()+
  #   theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
  #         legend.text=element_text(size=30),
  #         axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
  #         axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank()) + 
  #   xlab(expression(paste("r"))) + ylab("R K(r)")
  # 
  # doc = add_slide(doc, "Blank", "Office Theme")
  # doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
  # 
  # residual_G_model_2 = Gres(manual_hybrid_axon_model, correction="best")
  # ggobj = ggplot(data = residual_G_model_2, aes(r)) + 
  #   
  #   geom_ribbon(aes(ymin=lo, ymax=hi), alpha = 0.3)+
  #   
  #   geom_line(aes(y=theo), size=1, linetype="dashed") + 
  #   
  #   geom_line(aes(y=hres), size=1.2) + 
  #   
  #   theme_bw()+
  #   theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
  #         legend.text=element_text(size=30),
  #         axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
  #         axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30),
  #         panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank()) + 
  #   xlab(expression(paste("r"))) + ylab("R G(r)")
  # 
  # doc = add_slide(doc, "Blank", "Office Theme")
  # doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
  # 
  # print(doc, target = paste(folder_path, file_name, "_Res_Sum_Functs.pptx", sep=""))
  # 
  #--------------------------------------------------------------------------------------------------------------------------
  
  ####Finding local maxima
  peak_data = detectNegAICPeaks(p_info, folder_path, file_name, extension)
  
  true_peaks = peak_data[[1]]
  true_peaks_r = peak_data[[2]]
  
  hybrid_call_args = list()
  
  for (i in 1:length(true_peaks_r)) {
    if(true_peaks_r[i] == parameters(p)$r){
      hybrid_call_args[[i]] = StraussHard(r=true_peaks_r[i], h=parameters(p)$hc)
    }
    else{
      hybrid_call_args[[i]] = Strauss(r=true_peaks_r[i])
    }
  }
  
  hybrid_interact = do.call("Hybrid", hybrid_call_args)
  hybrid_axon_model = ppm(axon_pp ~ x+y+marks, interaction = hybrid_interact)
   
  return(list(hybrid_axon_model, ltest))
  
}


main <- function(){
  ####Library installation
  load.lib = c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", "poweRlaw",
               "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools",
               "officer", "rvg", "oce", "OneR", "RandomFieldsUtils", "RandomFields", "Cairo", "knitr", "scorepeak", "Rcpp", "emdist")
  
  install.lib = load.lib[!load.lib %in% installed.packages()]
  for(lib in install.lib) install.packages(lib, dependencies=TRUE)
  sapply(load.lib, require, character=TRUE)
  
  ####Data locations
  folder_path = "D:/Summer 2019/R codes/VagusNerve/Data/Rat vagus nerve_EM morphometric analysis_10-12-2020/AbdominalCervicalAxonLocations/"
  file_name = "AxonLocations_898_AbdVagAntGast_FASCICLE2"
  extension = ".xlsx"
  
  doc = read_pptx()
  
  fascicle_contour = readRDS(paste(folder_path, file_name, "_test.rds", sep=""))
  
  hybrid_ret_val = generateHybridModel(folder_path, file_name, extension, fascicle_contour)
  hybrid_axon_model = hybrid_ret_val[[1]]
  ltest = hybrid_ret_val[[2]]
  
  axon_simulated_hybrid = simulate(hybrid_axon_model, nsim = 1)
  plot(axon_simulated_hybrid$`Simulation 1`, main=paste("Simulated Hybrid", file_name, sep = " "), cex=0.7)
  
  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(code = plot(axon_simulated_hybrid$`Simulation 1`, main=paste("Simulated Hybrid", file_name, sep = " "), cex=0.7)), 
                location = ph_location_fullsize())
  
  simulated_ltest_hybrid = envelope(axon_simulated_hybrid$`Simulation 1`, Lest, nsim=39)
  
  ggobj = ggplot(data = simulated_ltest_hybrid, aes(r)) + 
    
    geom_ribbon(aes(ymin=lo-r, ymax=hi-r), alpha = 0.3)+
    geom_line(aes(y=lo-r, colour="CSR lower bound"), size=2) +
    geom_line(aes(y=hi-r, colour="CSR upper bound"), size=2) +
    #geom_point(aes(y=theo, colour="csr")) + 
    
    geom_line(aes(y=obs-r, colour="simulated hybrid"), size=2) + 
    #geom_point(aes(y=obs, colour="observed")) +  
    
    theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(), 
          legend.text=element_text(size=30),
          axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30)) + 
    xlab(expression(paste("Distance, r (", mu, "m)"))) + ylab("L-Function")
  
  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())
  
  print(doc, target = paste(folder_path, file_name, "_Simulated_Hybrid_Local_Maxima.pptx", sep=""))
  
  ####EMD between real and simulated
  simulated_hybrid_emd_mat = as.matrix(data.frame(r=simulated_ltest_hybrid$r, obs=simulated_ltest_hybrid$obs))
  original_emd_mat = as.matrix(data.frame(r=ltest$r, obs=ltest$obs))
  emd_val = emd(A=original_emd_mat, B=simulated_hybrid_emd_mat, dist="euclidean")
  # cat(emd_val)
  
  ####checking confidence interval for CSR
  original_emd_mat_low = as.matrix(data.frame(r=ltest$r, obs=ltest$lo))
  original_emd_mat_high = as.matrix(data.frame(r=ltest$r, obs=ltest$hi))
  emd_interval = emd(A=original_emd_mat_low, B=original_emd_mat_high, dist="euclidean")
  # cat(emd_interval)
  
  cat(emd_val, " ", emd_interval, "\n")
  
  doc = read_pptx()
  residual_K_model_2 = Kres(hybrid_axon_model, correction="best")
  ggobj = ggplot(data = residual_K_model_2, aes(r)) +

    geom_ribbon(aes(ymin=ilo, ymax=ihi), alpha = 0.3)+

    geom_line(aes(y=theo), size=1, linetype="dashed") +

    geom_line(aes(y=ires), size=1.2) +

    theme_bw()+
    theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(),
          legend.text=element_text(size=30),
          axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30),
          panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab(expression(paste("r"))) + ylab("R K(r)")

  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())

  residual_G_model_2 = Gres(hybrid_axon_model, correction="best")
  ggobj = ggplot(data = residual_G_model_2, aes(r)) +

    geom_ribbon(aes(ymin=lo, ymax=hi), alpha = 0.3)+

    geom_line(aes(y=theo), size=1, linetype="dashed") +

    geom_line(aes(y=hres), size=1.2) +

    theme_bw()+
    theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=20), legend.title=element_blank(),
          legend.text=element_text(size=30),
          axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30),
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30),
          panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab(expression(paste("r"))) + ylab("R G(r)")

  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(ggobj = ggobj), location = ph_location_fullsize())

  print(doc, target = paste(folder_path, file_name, "_Res_Sum_Functs_Hybrid.pptx", sep=""))

  
}
