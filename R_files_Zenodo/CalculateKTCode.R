

# Takes 4-bar dimensions including diagonal and proposed rotation
# Returns standardized shape to Fixed = 1, kt, and xy coordinates of starting and ending configurations
# diagonal <- sqrt(Fixed^2 + In^2 - 2*Fixed*In*cos(theta_start*2*pi/360))
# Rotation in degrees

calc_KT_operc <- function(Fixed, In, Coupler, Out, Diagonal, rotation = 5){
  Lf <- 1
  Li <- In/Fixed
  Lc <- Coupler/Fixed
  Lo <- Out/Fixed
  Ld <- Diagonal/Fixed
  std.shape <- data.frame(Lf, Li, Lc, Lo, Ld)
  # All angles in radians
  rotation_rad <- rotation*pi/180
  
  # Test for a valid 4bar:
  if(Ld < Lc + Lo & Ld < Lf + Li) valid.start <- T else valid.start <- F
  
  if(valid.start){
    theta_in_start <- acos((Li^2 + 1 - Ld^2)/(2*Li*1))
    theta_1_start <-acos((1 + Ld^2 - Li^2)/(2*1*Ld))
    theta_2_start <-acos((Lo^2 + Ld^2 - Lc^2)/(2*Lo*Ld))
    theta_out_start <- theta_1_start + theta_2_start
  } else theta_out_start <- NA
  
  # Calculations for ending output angle
  if(is.na(theta_out_start) == F ){
  theta_in_end <- theta_in_start + rotation_rad
  Ld_end <- sqrt(1 + Li^2 - 2*1*Li*cos(theta_in_end))
  if(Ld_end < Lc + Lo & Ld_end < Lf + Li) valid.end <- T else valid.end <- F
  if(valid.end){
    theta_1_end <-acos((1 + Ld_end^2 - Li^2)/(2*1*Ld_end))
    theta_2_end <-acos((Lo^2 + Ld_end^2 - Lc^2)/(2*Lo*Ld_end))
    theta_out_end <- theta_1_end + theta_2_end
  } else theta_out_end <- NA
  
  # calculate kt
  kt <- abs(theta_out_end - theta_out_start)/rotation_rad
  
  # calculate starting configuration xy coordinates using fixed link as a horizontal line, fish facing right
  Input_Fixed_joint <- c(0,0) # same for starting and ending configuration
  Fixed_Output_joint <- c(Lf, 0) # same for starting and ending configuration
  theta_3_start <- pi/2 - theta_in_start
  theta_3_end <- pi/2 -theta_in_end
  Coupler_Input_joint_start <- c(Li*sin(theta_3_start),  -Li*cos(theta_3_start))
  Coupler_Input_joint_end <- c(Li*sin(theta_3_end),  -Li*cos(theta_3_end))
  theta_4_start <- pi - theta_out_start
  theta_4_end <- pi - theta_out_end
  Out_Coupler_joint_start <- c(Lf + Lo *cos(theta_4_start),  -Lo*sin(theta_4_start))
  Out_Coupler_joint_end <-  c(Lf + Lo *cos(theta_4_end),  -Lo*sin(theta_4_end))
  start.config <- rbind(Input_Fixed_joint, Fixed_Output_joint, Coupler_Input_joint_start, Out_Coupler_joint_start)
  colnames(start.config) <- c("x", "y")
  end.config <- rbind(Input_Fixed_joint, Fixed_Output_joint, Coupler_Input_joint_end, Out_Coupler_joint_end)
  colnames(end.config) <- c("x", "y")
  return(list(KT = kt , shape = std.shape, start.coord = start.config, end.coord = end.config))
  } else  return(list(KT = NA , shape = NA, start.coord = NA, end.coord = NA))
}


# Test example
Fixed <- 7.309
Coupler <- 6.91
In <- 4.7918
Out <- 0.895
Diagonal <- 7.2
rotation <- 5  # degrees
output <- calc_KT_operc(Fixed, In, Coupler, Out, Diagonal , rotation=rotation )
#






plot4bar <- function(coords){
  coords <- as.data.frame(coords)
  scaling <- 1.5
  edge <- 0.24
  xlims = range(coords[,1])*scaling + c(-edge, edge)
  ylims = range(coords[,2])*scaling  + c(-edge, edge)
  par(mar = c(0,0,0,0))
  plot(0~0, xlim = xlims, ylim = ylims, axes = F, xlab = "", ylab = "", cex = 0)
  line.dat.Fixed <- coords[c(1,2),] #Fixed
  lines(x = line.dat.Fixed$x, y = line.dat.Fixed$y, lwd = 4)
  line.dat.Out <- coords[c(2,4),] #Out
  lines(x = line.dat.Out$x, y = line.dat.Out$y, lwd = 4)
  line.dat.Coupler <- coords[c(4,3),] #Coupler
  lines(x = line.dat.Coupler$x, y = line.dat.Coupler$y, lwd = 4)
  line.dat.In <- coords[c(3,1),] #In
  lines(x = line.dat.In$x, y = line.dat.In$y, lwd = 4)
  text(x = mean(line.dat.Fixed[,1]), y  = mean(line.dat.Fixed[,2]) + 0.02, "Fixed")
  text(x = mean(line.dat.Out[,1]) + 0.25, y  = mean(line.dat.Out[,2]) + 0.02, "Output")
  text(x = mean(line.dat.Coupler[,1])+0.25, y  = mean(line.dat.Coupler[,2]) - 0.02, "Coupler")
  text(x = mean(line.dat.In[,1])-0.18, y  = mean(line.dat.In[,2]) + 0.02, "Input")
}
  
  # Test example
plot4bar(output$end.coord)





plot4bar_compare2 <- function(coords1, coords2, color1 = "black", color2 = "blue", lty1 = 1, lty2 = 2, lwd1 = 4, lwd2 = 4){
  coords1 <- as.data.frame(coords1)
  coords2 <- as.data.frame(coords2)
  coords.both <- rbind(coords1, coords2)
  scaling <- 1.5
  edge <- 0.24
  xlims = range(coords.both[,1])*scaling + c(-edge, 0)
  ylims = range(coords.both[,2])*scaling  + c(-0, edge)
  par(mar = c(0,0,0,0))
  plot(0~0, xlim = xlims, ylim = ylims, axes = F, xlab = "", ylab = "", cex = 0)
  
  #Plot First one
  line.dat.Fixed <- coords1[c(1,2),] #Fixed
  lines(x = line.dat.Fixed$x, y = line.dat.Fixed$y, lwd = lwd1, col = color1, lty = lty1)
  line.dat.Out <- coords1[c(2,4),] #Out
  lines(x = line.dat.Out$x, y = line.dat.Out$y, lwd = lwd1, col = color1, lty = lty1)
  line.dat.Coupler <- coords1[c(4,3),] #Coupler
  lines(x = line.dat.Coupler$x, y = line.dat.Coupler$y, lwd = lwd1, col = color1, lty = lty1)
  line.dat.In <- coords1[c(3,1),] #In
  lines(x = line.dat.In$x, y = line.dat.In$y, lwd = lwd1, col = color1, lty = lty1)
  text(x = mean(line.dat.Fixed[,1]), y  = mean(line.dat.Fixed[,2]) + 0.02, "Fixed")
  text(x = mean(line.dat.Out[,1]) + 0.25, y  = mean(line.dat.Out[,2]) + 0.02, "Output")
  text(x = mean(line.dat.Coupler[,1])+0.25, y  = mean(line.dat.Coupler[,2]) - 0.02, "Coupler")
  text(x = mean(line.dat.In[,1])-0.18, y  = mean(line.dat.In[,2]) + 0.02, "Input")
  
  #Plot Second one
  line.dat.Fixed <- coords2[c(1,2),] #Fixed
  lines(x = line.dat.Fixed$x, y = line.dat.Fixed$y, lwd = lwd2, col = color2, lty = lty2)
  line.dat.Out <- coords2[c(2,4),] #Out
  lines(x = line.dat.Out$x, y = line.dat.Out$y, lwd =  lwd2, col = color2, lty = lty2)
  line.dat.Coupler <- coords2[c(4,3),] #Coupler
  lines(x = line.dat.Coupler$x, y = line.dat.Coupler$y, lwd = lwd2, col = color2, lty = lty2)
  line.dat.In <- coords2[c(3,1),] #In
  lines(x = line.dat.In$x, y = line.dat.In$y, lwd =  lwd2, col = color2, lty = lty2)
}

# Example to test
plot4bar_compare2( output$start.coord, output$end.coord)





# Need to build a list with many elements, each element is a output$start.coord or equivalent

plot4bar_compare_Many <- function(coords_list, color1 = "black", lty1 = 1,  lwd1 = 2, plot.consensus = T, consensus.color = "red", consensus.lwd = 6){
 xlims <- c(-0.25, 1.5)
 ylims <- c(-1.5, 0.25)
  par(mar = c(0,0,0,0))
  plot(0~0, xlim = xlims, ylim = ylims, axes = F, xlab = "", ylab = "", cex = 0)
  n_to_plot <- length(coords_list)
  
  for(i in 1:n_to_plot){
    coords1 <- as.data.frame(coords_list[[i]])
    line.dat.Fixed <- coords1[c(1,2),] #Fixed
    lines(x = line.dat.Fixed$x, y = line.dat.Fixed$y, lwd = lwd1, col = color1, lty = lty1)
    line.dat.Out <- coords1[c(2,4),] #Out
    lines(x = line.dat.Out$x, y = line.dat.Out$y, lwd = lwd1, col = color1, lty = lty1)
    line.dat.Coupler <- coords1[c(4,3),] #Coupler
    lines(x = line.dat.Coupler$x, y = line.dat.Coupler$y, lwd = lwd1, col = color1, lty = lty1)
    line.dat.In <- coords1[c(3,1),] #In
    lines(x = line.dat.In$x, y = line.dat.In$y, lwd = lwd1, col = color1, lty = lty1)
    if(i == 1){
      text(x = mean(line.dat.Fixed[,1]), y  = mean(line.dat.Fixed[,2]) + 0.02, "Fixed")
      text(x = mean(line.dat.Out[,1]) + 0.25, y  = mean(line.dat.Out[,2]) + 0.02, "Output")
      text(x = mean(line.dat.Coupler[,1])+0.25, y  = mean(line.dat.Coupler[,2]) - 0.02, "Coupler")
      text(x = mean(line.dat.In[,1])-0.18, y  = mean(line.dat.In[,2]) + 0.02, "Input")
    }
  }
  if(plot.consensus){
    mean.coords <- matrix(0, nrow = 4, ncol = 2)
    for(i in 1:n_to_plot){ mean.coords <- mean.coords + coords_list[[i]]}
    mean.coords <- as.data.frame(mean.coords/n_to_plot)
    line.dat.Fixed <- mean.coords[c(1,2),] #Fixed
    lines(x = line.dat.Fixed$x, y = line.dat.Fixed$y, lwd = consensus.lwd, col = consensus.color, lty = lty1)
    line.dat.Out <- mean.coords[c(2,4),] #Out
    lines(x = line.dat.Out$x, y = line.dat.Out$y, lwd = consensus.lwd, col = consensus.color, lty = lty1)
    line.dat.Coupler <- mean.coords[c(4,3),] #Coupler
    lines(x = line.dat.Coupler$x, y = line.dat.Coupler$y, lwd = consensus.lwd, col = consensus.color, lty = lty1)
    line.dat.In <- mean.coords[c(3,1),] #In
    lines(x = line.dat.In$x, y = line.dat.In$y, lwd = consensus.lwd, col = consensus.color, lty = lty1)
    
  }
}

# Example 
coords_list <- list(output$start.coord, output$end.coord)
plot4bar_compare_Many(coords_list, lwd1 = 0.8)
