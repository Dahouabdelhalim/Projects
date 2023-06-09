######
###### Load libraries and define functions
###### 

library(tidyverse)
library(lubridate)
filter <- dplyr::filter
lag <- dplyr::lag

## Data
world <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")

## ggplot default
theme_set(theme_bw()+theme(strip.background = element_blank(),panel.grid = element_blank()))

## Utility functions
# Transform correlation into guassian variable
logitr <- function(r) {
  p <- (r+1)/2
  return(log(p/(1-p)))
}
# Inverse transformation
logitr_inv <- function(l){
  p <- 1/(1+exp(-l))
  return(p*2-1)
}

# Standardization
std <- function(x, na.rm = T) (x-mean(x, na.rm = na.rm))/sd(x, na.rm = T)
# Find the linear trend in a temporal series
trd <- function(x){
  if(sum(!is.na(x))>1){
    time <- 1:length(x)
    m <- lm(x~time, na.action = na.exclude)
    trd <- predict(m)
    return(trd)
  } else {
    return(x)
  }
}
# Remove linear trend from a temporal series
dtrd <- function(x){
  if(sum(!is.na(x))>1){
    time <- 1:length(x)
    m <- lm(x~time, na.action = na.exclude)
    # res <- x-time*m$coef[2]
    res <- residuals(m)
    return(res)
  } else {
    return(x)
  }
}

# Transform p-value into significance symbols
signif <- function(x) if(x>.1){""}else if(x>.05){""}else if(x>.01){"*"}else if(x>.001){"**"}else{"***"}

# 
max12 <- function(x) if(x>12){x-12}else{x}

# Calculate the variance inflation factor
VIF <- function(x){
  vars <- colnames(x)
  m <- map(seq_along(vars), function(i) lm(paste0(vars[i],"~",paste(vars[-i],collapse = "+")), data = x))
  rsq <- map(m, function(x) summary(x)$r.squared) %>% unlist
  vif <- map(rsq, function(x) 1/(1-x)) %>% unlist
  res <- data.frame(Var = vars, rsq = rsq, VIF = vif)
  return(res)
}

# Extract coefficients, associated standard errors and calculates p values from an lmer object
coefs.lmer <- function(x, norm.approx = T){
  estimate <- summary(x)$coefficients[,1]
  se <- summary(x)$coefficients[,2]
  tval <- summary(x)$coefficients[,3]
  if(!norm.approx){
    # get the Kenward-Roger-approximated degrees of freedom
    ddf <- pbkrtest::get_Lb_ddf(x, fixef(x))
    pval <- 2*(1-pt(abs(tval),ddf))
  }else{
    # assuming large ddf, the student distribution tends towards the normal distribution
    pval <- 2*(1-pnorm(abs(tval)))
  }
  return(data.frame(name = names(tval), estimate = estimate, se = se, tval = tval, pval = pval, signif = map(pval, signif) %>% unlist, row.names = NULL))
}

# Same for lm object
coefs.lm <- function(x, norm.approx = T){
  names <- rownames(summary(x)$coefficients)
  slope <- summary(x)$coefficients[,1]
  pval <- summary(x)$coefficients[,4]
  return(data.frame(name = names, slope = slope, pval = pval))
}

# Get the precision of a number
ndec <- function(x) {
  xchar <- as.character(x)
  dec <- str_extract(xchar, "(?<=\\\\.)[[:digit:]]+")
  res <- nchar(dec)
  res[is.na(res)] <- 0
  return(res)
}

# Find Whittaker biome for given MAP and MAT
Whittaker_find <- function(MAT, MAP){
  df <- data.frame(MAT = MAT, MAP = MAP/10)
  points <- sp::SpatialPoints(as.matrix(df))
  biomes <- sp::over(points, plotbiomes::Whittaker_biomes_poly)
  biomes <- factor(biomes$biome)
  return(biomes)
}

# Find Koppen climatic zones given coordinates
climatezones <- kgc::climatezones
KoppenGeiger_find <- function(Lon, Lat){
  kgc::LookupCZ(data = data.frame(Site = seq_along(Lon), 
                                  rndCoord.lon = kgc::RoundCoordinates(Lon), 
                                  rndCoord.lat = kgc::RoundCoordinates(Lat)))
}

