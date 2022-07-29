generateGangliaCenters <- function(intensity, interaction, interact_radius, hardcore_dist, window, process_type=3, with_model=0, fitted_model){
  
  if(with_model == 0){
    if(process_type == 1){
      ganglia_centers = rHardcore(beta=intensity, R=hardcore_dist, W=window)
      return(ganglia_centers)
    }
    else if(process_type == 2){
      ganglia_centers = rStrauss(beta=intensity, gamma=interaction, R=interact_radius, W=window)
      return(ganglia_centers)
    }
    else if(process_type == 3){
      ganglia_centers = rStraussHard(beta=intensity, gamma=interaction, R=interact_radius, H=hardcore_dist, W=window)
      return(ganglia_centers)
    }
    else{
      stop("Invalid input for point process type!\n")
    }
  }
  
  if(with_model == 1){
    ganglia_centers = rmh.ppm(model = fitted_model, w=window)
    return(ganglia_centers)
  }
}

plotGeneratedGanglia <- function(ganglia_centers){
  plot(ganglia_centers, main="Point Pattern of Ganglia Centers", cex=0.75, pch=20)
}

savePlotGeneratedGanglia <- function(ganglia_centers){
  
}
