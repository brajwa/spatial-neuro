generateNeuronCenters <- function(path_to_intensity_image, image_name, b_val, g_val, r_val, h_val){
  
  path_to_intensity_image = "D:\\Summer 2019\\R codes\\Research1.0InterganglionicNetwork2021\\Inputs\\Intensity Profile\\"
  image_name = "Thick_mask_simulation_3.jpg"
  
  b_val = 0.003
  h_val = 0.01
  
  g_val =0.78
  r_val = 0.0285
  
  # intensity_profile = cimg2im(load.image(paste(image_path, image_name, sep="")), W=simulated.ppp$window)
  intensity_profile = cimg2im(load.image(paste(path_to_intensity_image, image_name, sep="")))
  
  #simulated_neuron_1 = rHardcore(beta = b_val, R = h_val, W=owin(intensity_profile$xrange, intensity_profile$yrange))
  simulated_neuron_1 = rStraussHard(beta = b_val, gamma = g_val, R = r_val, H = h_val, W=owin(intensity_profile$xrange, intensity_profile$yrange))
  
  plot(simulated_neuron_1, bg='magenta', pch=21, cex=0.9)
  
  simulated_neuron = rthin(X=simulated_neuron_1, P=intensity_profile)
  
  plot(simulated_neuron, bg='magenta', pch=21, cex=0.9)
  
  simulated_neuron$y = simulated_neuron$y - simulated_neuron$window$yrange[2]
  simulated_neuron$x = simulated_neuron$x - simulated_neuron$window$xrange[1]
  simulated_neuron$window$xrange = simulated_neuron$window$xrange- simulated_neuron$window$xrange[1]
  simulated_neuron$window$yrange = simulated_neuron$window$yrange- simulated_neuron$window$yrange[2]
  
  simulated_neuron$y = simulated_neuron$y / abs(simulated_neuron$window$yrange[1]) * abs(ganglia_ppp$window$yrange[1])
  simulated_neuron$x = simulated_neuron$x /  abs(simulated_neuron$window$xrange[2]) * abs(ganglia_ppp$window$xrange[2])
  simulated_neuron$window$xrange = c(0, abs(ganglia_ppp$window$xrange[2]))
  simulated_neuron$window$yrange = c(ganglia_ppp$window$yrange[1], 0)
  
  #simulated_neuron$y = simulated_neuron$y - 3735
  
  plot(g2_lin, lwd=2)
  plot(simulated_neuron, bg='magenta', pch=21, cex=0.6, add=T)
  plot(ganglia_ppp, bg='black', pch=21, cex=0.6, add=T)
  
  write.csv(data.frame(x=simulated_neuron$x, y=simulated_neuron$y), 
            "D://Summer 2019//R codes//Research1.0InterganglionicNetwork2021//Outputs//Simulated thingys//neuron_coord_new.csv", 
            row.names = F)
  
  doc = read_pptx()
  
  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(code = plot(simulated_neuron, bg='magenta', pch=21, cex=0.8)), location = ph_location_fullsize())
  
  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(code = ({plot(g2_lin)
    plot(ganglia_ppp, bg='black', pch=21, cex=0.9, add=T)})), location = ph_location_fullsize())
  
  doc = add_slide(doc, "Blank", "Office Theme")
  doc = ph_with(doc, dml(code = ({plot(g2_lin, lwd=2)
    plot(ganglia_ppp, bg='magenta', pch=21, cex=0.9, add=T)
    plot(simulated_neuron, bg='magenta', pch=21, cex=0.9, add=T)})), location = ph_location_fullsize())
  
  
  print(doc, target = paste(path_to_intensity_image, "Neuron_Simulated_new.pptx", sep=""))
  
}