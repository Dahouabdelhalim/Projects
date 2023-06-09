library(foreach) #version 1.5.2
library(nlme) #version 3.1-157

source('prep_data.R')

######The main function that separates within-taxa and across-taxa thermal sensitivities
#following Eq. 9 and 14 in Supplement 1 (see also Eq. 1 and 2 in the main text).
#Caution: Ea is used throughout in the following code to represent intra-specific activation energy
#which should be Eintra in the main manuscript and Supplement 1.
decomp <- function(oridat, remove.nonpositive = T, epsilon=0.01){
  
  #Weighed Covariance function
  wcov <- function(x, y, w){
    stopifnot(length(x) == length(y) && length(x) == length(w))
    stopifnot(round(sum(w),0) == 1L)   
    xbar <- sum(x*w)
    ybar <- sum(y*w)
    
    return(sum(w*(x-xbar)*(y-ybar)))
  }
  
  #How to deal with zero growth rate; the default is to remove all non-positive values
  if (remove.nonpositive) {
    oridat <- oridat %>% subset(Growth > 0)
  }else{
    oridat[oridat$Growth <= 0, 'Growth'] <- epsilon
  }
  
  Eapp       <- numeric(9L)  #Includes both terms of Einter equation and EL equation
  oridat$eps <- NA #Residual in Eqn. 1 (Supplement 1)
  oridat$xi  <- NA #Residual in Eqn. 12 (Supplement 1)
  
  #Extract each taxon (XoptL is the optimal temperature)
  unidat     <- plyr::ddply(oridat, .(ID), summarize,
                            ID       = ID[1],
                            Species  = Sp[1],
                            XoptL    = XoptL[1] )

  #grandmean of x
  GrandX <- mean(oridat$X)

  for (i in 1:nrow(unidat)){
    code              <- unidat$ID[i]
    
    #Pick up the data for the ith taxon
    tmp               <- oridat[oridat$ID == code,]
    unidat[i, 'mj'  ] <- nrow(tmp)
    unidat[i, 'Tbar'] <- mean(tmp$X)
    unidat[i, 'var1'] <- sum((tmp$X - GrandX)**2)
    
    #Calculate beta_{j} and E_{intra,j} 
    LM1               <- lm(log(Growth) ~ X, tmp)
    
    #LM2 is for computing y_{m,j} (see Eqn. 12 in Supplement 1) 
    tmp$X1            <- tmp$X - tmp$XoptL[1]
    LM2               <- lm(log(Growth) ~ X1, tmp)
    
    #Obtain beta_{j} in Eqn. 2 (Supplement 1)
    unidat[i, 'bj']   <- as.numeric(coef(LM1)[1])
    unidat[i, 'Ea']   <- as.numeric(coef(LM2)[2]) #Eintra (Ea) is the same for LM1 and LM2
    
    #Obtain y_{m,j} in Eqn. 12 (Supplement 1)
    unidat[i, 'ym']   <- as.numeric(coef(LM2)[1])

    #Add epsilon to oridat (Eqn. 1)
    oridat[oridat$ID == code, 'eps'] <- LM1$residuals
    
    #Add residuals (xi)  to oridat (Eqn. 12)
    oridat[oridat$ID == code, 'xi']  <- LM2$residuals
  }
 
  #Median Ea (Eintra)
  Ea_median <- median(unidat$Ea)
  
  #Total number of observations of the original data
  M <- sum(unidat$mj)
  
  #Number of taxa
  n <- nrow(unidat)
  
  #Calculate probabilities for each taxon
  p <- unidat$mj/M

  #variance of x in the original data
  VARx  <- var(oridat$X)
  
  #variance of x_bar with unequal probabilities
  VARxbar     <- sum(p*(unidat$Tbar - GrandX)^2)
  
  #variance of optimal temperature  with unequal probabilities
  #First, calculate mean Xm
  Mean_xm     <- sum(p*unidat$XoptL)
  VARxm       <- sum(p*(unidat$XoptL - Mean_xm)^2)
  
  #covariance between x_{m} and \\overline{x}
  COV_xm_xbar <- wcov(unidat$XoptL, unidat$Tbar,p)
  
  #covariance between x_{m}E_{a} and \\overline{x}
  unidat$EaXm   <- unidat$Ea * unidat$XoptL
  COV_Eaxm_xbar <- wcov(unidat$EaXm, unidat$Tbar, p)
  
  #First term of Eqn. 9 and Eqn. 14
  Eapp[1] <- sum(unidat$Ea * unidat$var1)/M/VARx

  #3rd term of Eqn. 14
  Eapp[2] <- -COV_Eaxm_xbar/VARx
  
  #Calculate Einter for phytoplankton
  Einter <- lm(bj ~ Tbar,  unidat, weights = p)  #weighed linear regression
  
  #Calculate EL for phytoplankton (Eqn. 14)
  EL <- lm(ym ~ XoptL, unidat, weights = p)  #weighed linear regression
  
  #Obtain residuals (nu) for EL (Eqn. 11)
  nu <- EL$residuals
  
  #Obtain residuals (beta_{j}) for Einter (Eqn. 2)
  beta <- Einter$residuals

  #check whether mean of weighed residuals is 0
  #sum(p*nu)

  EL     <- as.numeric(coef(EL)[2]) #Obtain EL
  Einter <- as.numeric(coef(Einter)[2]) #Obtain Einter
  
  # weighed variance of x
  # VARxbar1 <- wcov(unidat$Tbar, unidat$Tbar, p) confirms the calculation above
  
  #2nd term of Eqn. 14
  Eapp[3]  <- EL * COV_xm_xbar/VARx
  
  #4th term of Eqn. 14
  COVEXbar <-  wcov(unidat$Ea, unidat$Tbar, p) #weighed covariance between Ea (Eintra) and xbar 
  Eapp[4]  <-  GrandX  * COVEXbar  / VARx

  #5th term of Eqn. 14
  COVnuxbar   <- wcov(nu, unidat$Tbar, p)
  Eapp[5]     <- COVnuxbar/VARx
  
  #6th term of Eqn. 14
  COVXi_X  <- cov(oridat$X, oridat$xi)
  Eapp[6]  <- COVXi_X / VARx
  
  #2nd term of Eqn. 9 (EinterVar(Xbar)/Var(x))
  Eapp[7]  <- Einter * VARxbar /VARx
    
  #4th term of Eqn. 9 (Cov(beta, xbar)/Var(x))
  COVbetaxbar <- wcov(beta, unidat$Tbar, p)
  Eapp[8]     <- COVbetaxbar/VARx
  
  #5th term of Eqn. 9 (Cov(eps, X)/Var(x))
  COVepsX  <- cov(oridat$X, oridat$eps)
  Eapp[9]  <- COVepsX / VARx
  
  #run OLS regression for full data
  fullLM <- lm(log(Growth) ~ X, oridat) 
  
  return(list(n          = n,
              M          = M,
              VARx       = VARx,
              VARxm      = VARxm,
              VARxbar    = VARxbar,
              COVXmXbar  = COV_xm_xbar,
              COVEaxmxbar= COV_Eaxm_xbar,
              Einter     = Einter,
              EL         = EL, 
              XMEAN      = GrandX,
              COVEXbar   = COVEXbar,
              Ea_med     = Ea_median,
              EappLM     = as.numeric(coef(fullLM)[2]), 
              EappCal    = Eapp))
} 

#Use bootstrapping to estimate standard error of the above estimates
boot <- function(dat, Nrep = 100){
  result <- foreach(i=1:Nrep, .combine='rbind') %do% {
    
    #Randomly sample 60% of the taxa
    IDs    <- dat %>% select(ID) %>% unique()
    x      <- sample(IDs[,1], 0.6*nrow(IDs))
    x      <- subset(dat, ID %in% x)   #subsample
    x      <- decomp(x)
    
    c(x$EappLM, x$Ea_med, x$Einter, x$EL, sum(x$EappCal[c(1,4,7,8,9)]), 
      x$EappCal[1], x$EappCal[2], x$EappCal[3], x$EappCal[7])
  }
  MEAN <- apply(result,2,mean)
  SE   <- apply(result,2,sd)
  
  return(list(mean=MEAN, se=SE))
}

#OLS regression by considering cell size to generate Table S3
OLSsize <- function(dat){
  sdat <- dat %>% subset(Growth > 0)
  Z <- lm(log(Growth) ~ X + log(Volume), sdat)
  Z <- summary(Z)
  return(list(N   = sum(Z$df[1:2]),
             R2   = Z$r.squared, 
         Eapp     = coefficients(Z)[2,1],
         SE_Eapp  = coefficients(Z)[2,2],
         alpha    = coefficients(Z)[3,1],
         SE_alpha = coefficients(Z)[3,2]))
}
##################

PEuk_decomp       <- decomp(PEuk2) #Autotrophic protists
PEuk_decomp_boot  <-   boot(PEuk2) #Estimate the SE of each estimate by bootstrapping

#Check if Einter removing polar autotrophic protists
PEuk3             <- PEuk2 %>% subset(XoptL >= min(zoodat2$XoptL))
PEuk3_decomp      <- decomp(PEuk3)
OLSsize(PEuk2) #Check the size effect on Eapp of eukaryotic autotrophic protists

#Cyanobacteria
Cyn_decomp       <- decomp(Cyn2)
Cyn_decomp_boot  <-   boot(Cyn2) #Estimate the SE of each estimate by bootstrapping

#MicroZooplankton
mzoo_decomp      <- decomp(zoodat2)
mzoo_decomp_boot <-   boot(zoodat2)
OLSsize(zoodat2)

#Insect data in Rezende and Bozinovic (2019)
Ea_insect      <- decomp(Insect2)
Ea_insect_boot <- boot(Insect2)

#Heterotrophic bacteria in Smith et al. (2019)
Ea_HBac      <- decomp(HBac2)
Ea_HBac_boot <- boot(HBac2)

##Plot the temperature responses of each taxon(Fig. S1 and S2)
phypdffile <- 'AProtist_bytaxon_linear.pdf'
zoopdffile <- 'HProtist_bytaxon_linear.pdf'

plot_taxon <- function(dat, filename, caption = ''){

  dat$ID <- factor(dat$ID)
  
  sdat   <- dat %>% subset(Growth > 0)
  
  Ntaxon <- length(levels(dat$ID))
    
  YLAB <- expression(paste("Growth rate (" * d ^ -1 * ")"))
  
  pdf(filename,
      width = 7,
      height = 9,
      paper = 'a4')
  par(font.lab  = 1,
      family    = "serif",
      mgp       = c(2.2, 1, 0),
      mfrow     = c(4, 3),
      mar       = c(4, 4, 3.5, .2),
      oma       = c(4,4,0,0))
  
  for (i in 1:Ntaxon) {
 
      id   <- levels(dat$ID)[i]
      tmp  <- dat %>% subset(ID == id)
      
      #Linear fixed-effect model 
      LFE  <- lm(log(Growth) ~ X, data = tmp[tmp$Growth > 0,])
      
      #Essential to reset the newx
      minx <- min(tmp$Temperature)
      maxx <- max(tmp$Temperature)
      newx <- data.frame(Temperature = seq(minx, maxx, 0.01))
      newy <- data.frame('LFE' =numeric(nrow(newx)))
      
      #Find index for this taxon
      newy$LFE  <- exp(coef(LFE)[1] + coef(LFE)[2] * T.K(newx$Temperature))
      
      #Plotting original data points
      umax <- max(tmp$Growth, na.rm = T)

      plot(tmp$Temperature, tmp$Growth,
           las = 1,
           xlim = c(-2,   35),
           ylim = c(-0.05, umax + .2),
           xlab = 'Temperature (ÂºC)',
           ylab = YLAB,
           cex.lab = 1.1)
      lines(newx$Temperature, newy$LFE,  col='blue', lwd=1.3)
      txt  <- paste0(tmp$Sp[1])
      mtext(txt, adj = 0, cex=.7)  #Add caption
      text(34.4, umax*0.3, paste('Eintra = ', round( coef(LFE)[2], 2)), pos=2, col='blue')
  }
  #Add figure caption
  mtext(caption, side = 1, line = 0, outer = T, adj = 0)
  dev.off()
}

plot_taxon(PEuk2,   phypdffile, 'Fig. S1. Growth rate versus temperature plots for each autotrophic taxon')
plot_taxon(zoodat2, zoopdffile, 'Fig. S2. Growth rate versus temperature plots for each heterotrophic taxon')
