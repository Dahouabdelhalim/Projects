# =============================================================================
#
# This file contains general functions to generate figures used throughout this project
#
# figure functions overview
# 1) p.plot4bar           Plots a single, standardized Opercular 4-bar, with the Fixed Link set to one.
# 2) p.plot4bar_compare2  Plots two, standardized Opercular 4-bars.
# 3) p.plot4bar_compare_Many
# 4) p.plotPCxFunct             Plots lake-stream PC distances in lake-stream function distance 
# 5) p.plot.FunctVectorsxPC     Plots mean value for a given component-trait PC axis against the appropriate Function.
# =============================================================================
# 1) p.plot4bar
# -----
p.plot4bar <- function(coords){
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
#p.plot4bar(output$end.coord)
# -----
#############################################################
# 2) p.plot4bar_compare2
# -----
p.plot4bar_compare2 <- function(coords1, coords2, color1 = "black", color2 = "blue", lty1 = 1, lty2 = 2, lwd1 = 4, lwd2 = 4){
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
#p.plot4bar_compare2(output$start.coord, output$end.coord)
# -----
#############################################################
# 3) p.plot4bar_compare_Many
# -----
# Overlay multiple 4-bar plots on one another, with a mean 4-bar plotted.
# Need to build a list with many elements, each element is a output$start.coord or equivalent
p.plot4bar_compare_Many <- function(coords_list, color1 = "black", lty1 = 1,  lwd1 = 2, plot.consensus = T, consensus.color = "red", consensus.lwd = 6){
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
#coords_list <- list(output$start.coord, output$end.coord)
#p.plot4bar_compare_Many(coords_list, lwd1 = 0.8)
# -----
#############################################################
# 4) p.plotPCxFunct
# -----
#New plot of meanPC1.Lake - meanPC1.Stream for SI or KT components on x axis, meanSI/KT.Lake - meanSI/KT.Stream on y axis
p.plotPCxFunct <- function(PCaxisdata, Site, TrophicFunctionData, TrophicFunction, PCnumber, percent.explained){
  PC.means <- tapply(PCaxisdata, Site, mean, na.rm = TRUE)
  #names(PC.means)
  PC.distance <- vector()
  for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)) {
    j = i+1
    abcd <- PC.means[i]-PC.means[j]
    PC.distance <- cbind(PC.distance, abcd)
  }
  rownames(PC.distance) <- paste(TrophicFunction, PCnumber, "means.distance", sep = ".")
  x <- substr(names(PC.means), 1, 4)
  colnames(PC.distance) <- x[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)]

  
  Funct.means <- tapply(TrophicFunctionData, Site, mean, na.rm = TRUE)
  #names(Funct.means)
  Funct.distance <- vector()
  for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)) {
    j = i+1
    abcd <- Funct.means[i]-Funct.means[j]
    Funct.distance <- cbind(Funct.distance, abcd)
  }
  rownames(Funct.distance) <- paste(TrophicFunction, "SI.means.distance", sep = ".")
  y <- substr(names(Funct.means), 1, 4)
  colnames(Funct.distance) <- x[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)]

  #
  model <- lm(as.vector(Funct.distance)~as.vector(PC.distance))
  summary(model)
  par(las = 1, mar = c(6, 5, 2, 1))
  plot(Funct.distance~PC.distance, xlab = paste(TrophicFunction, "PC", PCnumber, ": L.mean - S.mean", sep = ""), ylab = paste(TrophicFunction, " L.mean - S.mean", sep = "."), main = paste("PC", PCnumber, " explains ", round(percent.explained*100, 1), "% of the variance in ", TrophicFunction, " components", sep = ""))
  abline(model)
  abline(h = 0, v = 0, lty = "dotted")
  legend.placement <- c(max(PC.distance), max(Funct.distance), min(PC.distance), min(Funct.distance))
    if (summary(model)$coefficients[2,1] > 1) positions <- c(legend.placement[1], legend.placement[2], 1, 1) else positions <- c(legend.placement[3], legend.placement[2], 0,1)
  text(positions[1], positions[2], paste("Linear Model\\nbeta = ",round(summary(model)$coefficients[2,1], 2), "\\nt = ", round(summary(model)$coefficients[2,3], 2),"\\ndf = ", model$df.residual, "\\np = ", round(summary(model)$coefficients[2,4], 2), "\\nAdj. R-squared = ", round(summary(model)$adj.r.squared, 2), sep = ""), adj = c(positions[3], positions[4]))
}

#Example
#PCaxisdata <- cleandatafinal$PC1.SIcomps
#Site <- cleandatafinal$site
#TrophicFunctionData <- cleandatafinal$SI
#TrophicFunction <- "Suction Index"
#PCnumber <- "1"

# -----

#############################################################
# 5) p.plot.FunctVectorsxPC
# -----
p.plot.FunctVectorsxPC <- function(PCaxisdata, Site, TrophicFunctionData, TrophicFunction, PCnumber){
  PC.means <- as.data.frame(tapply(PCaxisdata, Site, mean, na.rm = TRUE))
  colnames(PC.means) <- paste(TrophicFunction, "PC", PCnumber, sep = ".")
  PC.means$filename <- rownames(PC.means)
  Funct.means <- as.data.frame(tapply(TrophicFunctionData, Site, mean, na.rm = TRUE))
  colnames(Funct.means) <- paste(TrophicFunction, "coeff", sep = ".")  
  Funct.means$filename <- rownames(Funct.means)
  abcd <- merge(PC.means, Funct.means, "filename", all = TRUE)
  colour <- c("blue", "green")
  plot(abcd[,3]~(abcd[,2]), xlab = colnames(abcd[2]), ylab = colnames(abcd[3]), col = colour, pch = 16, main = paste("(non)Parallelism in ", TrophicFunction, " X PC", PCnumber, " space", sep = ""))
  for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)) {
    j = i + 1
    lines(matrix(data = c(abcd[i,2], abcd[j,2], abcd[i,3], abcd[j,3]), nrow = 2, ncol = 2), type = "l")
  }
}
# -----
#############################################################
# 6) p.plot.function.ordered
# -----
p.plot.function.ordered <- function(FunctionalData, Site, FunctionName, PCData, PCNumber){
  par(mar = c(4,4,1,1), mfrow = c(2,1))
  mean.site <- as.data.frame(tapply(FunctionalData, Site, mean, na.rm = TRUE))
    mean.site$sitenames <- rownames(mean.site)
    mean.site$color <- rep(c("blue", "green"), 16)
  sd.site <- as.data.frame(tapply(FunctionalData, Site, sd, na.rm = TRUE))
    sd.site$sitenames <- rownames(sd.site)
      mnPC <- as.data.frame(tapply(PCData, Site, mean, na.rm = TRUE))
      sdPC <- as.data.frame(tapply(PCData, Site, var, na.rm = TRUE))
        mnPC$sitenames <- rownames(mnPC)
        sdPC$sitenames <- rownames(sdPC)
  mnsd.site <- merge(mean.site, sd.site, by = "sitenames")
    mnsd.site <- merge(mnsd.site, mnPC, by = "sitenames")
    mnsd.site <- merge(mnsd.site, sdPC, by = "sitenames")
  colnames(mnsd.site) <- c("site", "mean", "color", "sd", "PCmean", "PCsd")
  
  #Order by site based on functional mean.
  mnsd.site <- mnsd.site[order(mnsd.site$mean),]
  
  #Plot Function
  ylimmin <- vector()
  ylimmax <- vector()
  for (i in 1:length(mnsd.site$mean)){
    ylimmin <- c(ylimmin, mnsd.site$mean[i] - mnsd.site$sd[i])
    ylimmax <- c(ylimmax, mnsd.site$mean[i] + mnsd.site$sd[i])
  }   
  plot(mnsd.site$mean, ylim = c(min(ylimmin)-0.1*abs(min(ylimmin)), max(ylimmax)+0.1*abs(max(ylimmax))), pch = 16, col = mnsd.site$color, ylab = FunctionName)
  for (i in 1:length(mnsd.site$mean)){
    lines(matrix(data = c(i, i, mnsd.site$mean[i], mnsd.site$mean[i] + mnsd.site$sd[i]), nrow = 2, ncol = 2))
    lines(matrix(data = c(i, i, mnsd.site$mean[i], mnsd.site$mean[i] - mnsd.site$sd[i]), nrow = 2, ncol = 2))
  }
  
  #Plot component PC, ordered by Functional means
  ylimmin <- vector()
  ylimmax <- vector()
  for (i in 1:length(mnsd.site$PCmean)){
    ylimmin <- c(ylimmin, mnsd.site$PCmean[i] - mnsd.site$PCsd[i])
    ylimmax <- c(ylimmax, mnsd.site$PCmean[i] + mnsd.site$PCsd[i])
  }   
  plot(mnsd.site$PCmean, ylim = c(min(ylimmin)-0.2*abs(min(ylimmin)), max(ylimmax)+0.2*abs(max(ylimmax))), pch = 16, col = mnsd.site$color, xlab = "PC means plotted in same order as Functional means", ylab = paste("PC", PCNumber, FunctionName, sep = " "))
  for (i in 1:length(mnsd.site$PCmean)){
    lines(matrix(data = c(i, i, mnsd.site$PCmean[i], mnsd.site$PCmean[i] + mnsd.site$PCsd[i]), nrow = 2, ncol = 2))
    lines(matrix(data = c(i, i, mnsd.site$PCmean[i], mnsd.site$PCmean[i] - mnsd.site$PCsd[i]), nrow = 2, ncol = 2))
  }
}


plot.Fig2.Col1 <- function(y,watershed, habitat, ylabel){
  group <- paste(habitat, watershed, sep = "_")
  y.means <- tapply(y, INDEX =group, mean, na.rm = T)
  y.sd <- tapply(y, INDEX =group, sd, na.rm = T)
  N <- summary(group)
  xvals <- c(rep(0, 16), rep(1, 16))
  plot(y.means ~ xvals, xlim = c(-0.25, 1.25), pch = 16, axes = F, xlab = "", ylab = ylabel, cex.lab = 1.1)
  box()
  axis(2)
  axis(1, at = c(0,1), labels = c("Lake", "Stream"))
  watershednames <- unique(watershed)
  for(i in 1:16){
    watershed.y <- y[watershed == watershednames[i]]
    watershed.habitat <- habitat[watershed == watershednames[i]]
    signif <- t.test(watershed.y ~ watershed.habitat)$p.value < 0.05
    lines(x = c(0,1), y = c(y.means[i] , y.means[i+16]), lty = 2 - signif, lwd = 1 + 2*signif)
  }
}
# plot.Fig2.Col1(y, watershed, habitat, ylabel = "Buccal suction index morphology PC1")

