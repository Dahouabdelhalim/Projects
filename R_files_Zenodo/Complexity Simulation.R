####### Define KT function

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






##################
# Sample morphospace
Nreps <- 100000
Lf <- rep(1, Nreps)
Li <- runif(n = Nreps, min = 0.44, max = 0.8) #empirical min max: 0.44, 0.80
Lc <- runif(n = Nreps, min = 0.77, max = 1.06)#empirical min max:  0.77, 1.06
Lo <- runif(n = Nreps, min = 0.04, max = 0.23)#empirical min max: 0.04, 0.23
Ld <- runif(n = Nreps, min = 0.75, max = 1.04)#empirical min max:  0.75, 1.04
KT.shape <- data.frame(Lf, Li, Lc, Lo, Ld)#empirical min max:

# Calculate KT for every valid 4-bar within the random morphospace
KTresult <- matrix(nrow = Nreps, ncol = 1)
for(i in 1:Nreps){
  Lf <- KT.shape$Lf[i]
  Li <- KT.shape$Li[i]
  Lc <- KT.shape$Lc[i]
  Lo <- KT.shape$Lo[i]
  Ld <- KT.shape$Ld[i]
  KTresult[i,1] <- calc_KT_operc(Lf, Li, Lc, Lo, Ld)$KT
}
KT.shape$KT <- KTresult

Nreps.actual <- (Nreps - sum(is.na(KT.shape$KT)))

# Test for correlation between KT and morphospace
modelKT <- lm(KT ~ Li + Lc + Lo + Ld, data = KT.shape, na.action = "na.omit")
summary(modelKT)



# Make an alfaro-type plot
library(fields)

hist(KT.shape$KT)
KT.shape.tr <- na.omit(KT.shape)
quantile(KT.shape.tr$KT, c(0.2, 0.5, 0.8))

# panel A
KT.min <- 3.22
KT.max <- 3.24
KT.surface <- KT.shape.tr[KT.shape.tr$KT > KT.min & KT.shape.tr$KT < KT.max, ]
dim(KT.surface)
head(KT.surface)
XY <- KT.surface[,2:3]
Z <- KT.surface$Lo
surface.fit <- Tps(XY, Z, m = 5)
surface <- predictSurface(surface.fit)
par(mar = c(5,5,2,1), oma = c(0,0,0,2))
plot.surface(surface, type = "p", xlab = "Li", ylab = "Lc" , zlab = "Lo", cex.lab = 1.5, main = "", add.legend = F, theta = 220, phi = 20, col = tim.colors(20))  # type = p, l, c, C, b

# Panel B
KT.min <- 4.4
KT.max <- 4.41
KT.surface <- KT.shape.tr[KT.shape.tr$KT > KT.min & KT.shape.tr$KT < KT.max, ]
dim(KT.surface)
head(KT.surface)
XY <- KT.surface[,2:3]
Z <- KT.surface$Lo
surface.fit <- Tps(XY, Z, m = 5)
surface <- predictSurface(surface.fit)
par(mar = c(5,5,2,1), oma = c(0,0,0,2))
plot.surface(surface, type = "p", xlab = "Li", ylab = "Lc" , zlab = "Lo", cex.lab = 1.5, main = "", add.legend = F, theta = 220, phi = 20)  # type = p, l, c, C, b


KT.min <- 6.73
KT.max <- 6.75
KT.surface <- KT.shape.tr[KT.shape.tr$KT > KT.min & KT.shape.tr$KT < KT.max, ]
dim(KT.surface)
head(KT.surface)
XY <- KT.surface[,2:3]
Z <- KT.surface$Lo
surface.fit <- Tps(XY, Z, m = 5)
surface <- predictSurface(surface.fit)
par(mar = c(5,5,2,1), oma = c(0,0,0,2))
plot.surface(surface, type = "p", xlab = "Li", ylab = "Lc" , zlab = "Lo", cex.lab = 1.5, main = "", add.legend = F, theta = 220, phi = 20, col = tim.colors(20), ticks = T)  # type = p, l, c, C, b



##################
# Repeat for SI
# Define function
f.SI <- function(CSAepaxial, Lin.momentarm.epaxial, Lout.momentarm.buccal, gape.width, buccal.length) {
  SI <- (CSAepaxial * (Lin.momentarm.epaxial/Lout.momentarm.buccal)) / (gape.width * buccal.length)
  return(SI)
}

Nreps.actual <- 100000
# Using empirically observed min & max
CSAepaxial <-  runif(n = Nreps.actual, min = 6.20, max = 18.69)
Lin.momentarm.epaxial <- runif(n = Nreps.actual, min = 1.06, max = 5.23) 
Lout.momentarm.buccal <- runif(n = Nreps.actual, min = 9.03, max = 25.20) 
gape.width <- runif(n = Nreps.actual, min = 1.64, max = 5.68) 
buccal.length <- runif(n = Nreps.actual, min = 6.37, max = 16.30) 

SIshape <- data.frame(CSAepaxial, Lin.momentarm.epaxial, Lout.momentarm.buccal, gape.width, buccal.length)


# Calculate SI within the random morphospace
SIresult <- matrix(nrow = Nreps.actual, ncol = 1)
for(i in 1:Nreps.actual){
  CSAepaxial <- SIshape$CSAepaxial[i]
  Lin.momentarm.epaxial <- SIshape$Lin.momentarm.epaxial[i]
  Lout.momentarm.buccal <- SIshape$Lout.momentarm.buccal[i]
  gape.width <- SIshape$gape.width[i]
  buccal.length <- SIshape$buccal.length[i]
  SIresult[i,1] <- f.SI(CSAepaxial, Lin.momentarm.epaxial , Lout.momentarm.buccal , gape.width, buccal.length)
}
SIshape$SI <- SIresult
# Test linear fit
modelSI <- lm(SI ~ CSAepaxial +Lin.momentarm.epaxial + Lout.momentarm.buccal + gape.width + buccal.length,data =  SIshape, na.action = "na.omit")

summary(modelKT)
summary(modelSI)




# Make an alfaro-type plot
library(fields)

hist(SIshape$SI)
SI.shape.tr <- na.omit(SIshape)
quantile(SI.shape.tr$SI, c(0.2, 0.5, 0.8))

# panel A
SI.min <- 0.03
SI.max <- 0.031
SI.surface <- SI.shape.tr[SI.shape.tr$SI > SI.min & SI.shape.tr$SI < SI.max, ]
dim(SI.surface)
head(SI.surface)
SI.surface.std <- SI.surface
SI.surface.std[,2:5] <- SI.surface[,2:5]/SI.surface[,1]
XY <- SI.surface.std[,2:3]
Z <- SI.surface.std$buccal.length
surface.fit <- Tps(XY, Z, m = 5)
surface <- predictSurface(surface.fit)
par(mar = c(5,5,2,2), oma = c(0,0,0,2))
plot.surface(surface, type = "p", xlab = "Epaxial input moment arm", ylab = "Buccal output moment arm" , zlab = "Buccal length", cex.lab = 1.5, main = "", add.legend = F, theta = 220, phi = 20, col = tim.colors(20))  # type = p, l, c, C, b

# Panel B
SI.min <- 0.056
SI.max <- 0.057
SI.surface <- SI.shape.tr[SI.shape.tr$SI > SI.min & SI.shape.tr$SI < SI.max, ]
dim(SI.surface)
head(SI.surface)
SI.surface.std <- SI.surface
SI.surface.std[,2:5] <- SI.surface[,2:5]/SI.surface[,1]
XY <- SI.surface.std[,2:3]
Z <- SI.surface.std$buccal.length
surface.fit <- Tps(XY, Z, m = 5)
surface <- predictSurface(surface.fit)
par(mar = c(5,5,2,2), oma = c(0,0,0,2))
plot.surface(surface, type = "p", xlab = "Epaxial input moment arm", ylab = "Buccal output moment arm" , zlab = "buccal.length", cex.lab = 1.5, main = "", add.legend = F, theta = 220, phi = 20, col = tim.colors(20))  # type = p, l, c, C, b


SI.min <- 0.11
SI.max <- 0.114
SI.surface <- SI.shape.tr[SI.shape.tr$SI > SI.min & SI.shape.tr$SI < SI.max, ]
dim(SI.surface)
head(SI.surface)
SI.surface.std <- SI.surface
SI.surface.std[,2:5] <- SI.surface[,2:5]/SI.surface[,1]
XY <- SI.surface.std[,2:3]
Z <- SI.surface.std$buccal.length
surface.fit <- Tps(XY, Z, m = 5)
surface <- predictSurface(surface.fit)
par(mar = c(5,5,2,2), oma = c(0,0,0,2))
plot.surface(surface, type = "p", xlab = "Epaxial input moment arm", ylab = "Buccal output moment arm" , zlab = "buccal.length", cex.lab = 1.5, main = "", add.legend = F, theta = 220, phi = 20, col = tim.colors(20))  # type = p, l, c, C, b










##################
# Repeat for LR
# Define function
f.jaw.leverratio.SI <- function(ClosingInLever, OutLever){
  LeverRatio <- OutLever/ClosingInLever
}

min(data.1a$jaw.opening.in.lever.length.mm)

# Using empirically observed min & max
InL <-  runif(n = Nreps.actual, min = 0.47, max = 1.58)
OutL <- runif(n = Nreps.actual, min = 1.74, max = 6.61) 

LRshape <- data.frame(InL, OutL)


# Calculate SI within the random morphospace
LRresult <- matrix(nrow = Nreps.actual, ncol = 1)
for(i in 1:Nreps.actual){
  ClosingInLever <- LRshape$InL[i]
  OutLever <- LRshape$OutL[i]
  LRresult[i,1] <- f.jaw.leverratio.kt(ClosingInLever, OutLever)
}
LRshape$LR <- LRresult
# Test linear fit
modelLR <- lm(LR ~ InL + OutL, data =  LRshape, na.action = "na.omit")
summary(modelLR)
summary(modelKT)
summary(modelSI)



# PCAs
KT.shape <- na.omit(KT.shape)
dim(KT.shape)
KTPCA <- prcomp(as.matrix(KT.shape[,2:5]), scale = T)
summary(KTPCA)

#SIPCA <- prcomp(as.matrix(SIshape[, c(1,3:5)]), scale = T)
SIPCA <- prcomp(as.matrix(SIshape[, c(1:5)]), scale = T)
summary(SIPCA)

summary(lm(KT.shape$KT ~ KTPCA$x[,1:3]))
summary(lm(SIshape$SI ~ SIPCA$x[,1:3]))
