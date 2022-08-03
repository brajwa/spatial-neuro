###############################################################################################################################################
#### generates spatial point pattern according to the given parameters
#### with_model: 0 indicates the other parameter fitted_model is not provided or not to be used, 1 indicated the provided fitted_model is to be used
#### fitted_model: a spatial model object
#### process_type: to be used while generating point pattern without provided model. can take value 1 or 2 or 3.
####                1 indicates hardcore process, 2 indicates Strauss process, 3 indicates hardcore-Strauss process.
####                when with_model=1, the value of process_type doesn't make any difference.
generateGangliaCenters <- function(intensity, interaction, interact_radius, hardcore_dist, window, process_type=3, 
                                   with_model=0, fitted_model){
  
  if(with_model == 0){
        print("generating ganglia centers with parameters...")
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
        print("generating ganglia centers with model...")
        ganglia_centers = rmh.ppm(model = fitted_model, w=window)
        return(ganglia_centers)
  }
}


##################################################
plotGeneratedGanglia <- function(ganglia_centers){
  plot(ganglia_centers, main="Point Pattern of Ganglia Centers", cex=0.75, pch=20)
}

