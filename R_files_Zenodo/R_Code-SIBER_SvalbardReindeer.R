##title: "SIBER ON SVALBARD REINDEER DATA"
##authors: "Tamara Hiltunen"

### GETTING STARTED

#Remove previously loaded items from the current environment and remove previous graphics.
rm(list=ls())
graphics.off()

#Set the seed each time so that the results are comparable. This is useful as 
#it means that anyone that runs this code, *should* get the same results as you, 
#although random number generators change from time to time.
set.seed(1)

#Load the Packages : SIBER_2.1.4, viridis_0.5.1 and readr_1.3.1

library(SIBER)
library(tidyverse)
library(dplyr)
library(viridis)
library(readr)
library(siar)

#Load in the data and wrangle into the required structure of:
#Columns: "iso1"      "iso2"      "group"     "community"

#For this data set iso 1 is d13C and iso2 is d15N
#group: Year (1995 - 2012, excl 2003 & 2020 thus 1-16)
#community = 1 as only from one community

#Years as per the following
#1 = 1995, 2 = 1996,  3 = 1997,  4 = 1998,  5 = 1999, 6 = 2000
#7 = 2001, 8 = 2002,  9 = 2004,  10 = 2005, 11 = 2006, 12 = 2007
#13 = 2008, 14 = 2009, 15 = 2011, 16 = 2012
#NB there is no data for 2003 and 2010

#With Viridris create a new palette with three colors for this data set
#TheN create the SIBER_OBJECT

palette(viridis(16))

mydata <- read.csv("Reindeer_data_schubert_correction.csv", 
                   header = TRUE, 
                   stringsAsFactors = FALSE)

names(mydata)

mydatasider <- mydata %>% dplyr::select(d13C, d15N, year)

mydatay <- mydatasider %>% rename(iso2 = d15N, iso1 = d13C, group = year)
mydatay <- mydatay %>% add_column(community = 1)
#mydatay$community <- as.integer(mydatay$group)

names(mydatay)
str(mydatay)

siber.year <- createSiberObject(mydatay)

### PLOTTING THE RAW DATA

community.hulls.args <- list(col = 1, lty = 1, lwd = 3)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

legend <- c("1995", "1996","1997","1998","1999","2000", "2001","2002","2004","2005","2006", "2007","2008","2009","2011","2012")
colours <- palette(viridis(16))
pchs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)

par(mfrow=c(1,1))
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))

plotSiberObject(siber.year,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                cex = 5,
                cex.lab = 1.5,
                cex.axis = 1.5,
                mai = c(1,1,1,1),
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\\u2030'),
                ylab = expression({delta}^15*N~'\\u2030'),
                x.limits = c(-26.5,-23.5),
                y.limits = c(1.5,6.5),
                main = "SIBER ellipses for each year")

#add legend to the plot
legend("topright", legend, col = colours, pch = pchs, bty = "n", cex =1, y.intersp=0.5)

##Summary statistics for each group: **TA, SEA and SEAc**

group.ML.y <- groupMetricsML(siber.year)
print(group.ML.y)

#RESULTS
#1.1995    1.1996    1.1997    1.1998    1.1999    1.2000    1.2001
#TA   0.6641000 0.9110000 0.9270000 1.0632500 0.6442500 0.6141000 0.7054500
#SEA  0.2444301 0.4625143 0.4177700 0.4602018 0.3362794 0.2604133 0.3610049
#SEAc 0.2632325 0.5203286 0.4557491 0.4956019 0.3699074 0.2804451 0.4011166
#1.2002    1.2004    1.2005    1.2006    1.2007    1.2008   1.2009
#TA   1.0026500 1.4149500 1.0017000 0.4200000 1.2626500 2.3055000 0.625650
#SEA  0.3671054 0.4885762 0.4082176 0.1825576 0.4231047 0.7832371 0.261066
#SEAc 0.3953443 0.5261590 0.4396190 0.1977707 0.4583634 0.8434861 0.281148
#1.2011    1.2012
#TA   0.8188000 0.8400500
#SEA  0.3055910 0.3150710
#SEAc 0.3259637 0.3316537


## FITTING BAYESIAN MODELS TO THE DATA

#Run parameters
parms <- list()
parms$n.iter <- 1 * 10^5   # number of iterations to run the model for
parms$n.burnin <- 2 * 10^3 # discard the first set of values
parms$n.thin <- 20     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains
parms$output = TRUE
parms$save.output = TRUE

parms$save.dir = getwd()

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior.y <- siberMVN(siber.year, parms, priors)

### COMPARING AMONG GROUPS USING STANDARD ELLIPSE AREA 

#The posterior estimates of the ellipses for each group can be used to
#calculate the SEA.B for each group.

SEA.B.y <- siberEllipses(ellipses.posterior.y)

siberDensityPlot(SEA.B.y, 
                 xticklabels = c("1995", "1996","1997","1998","1999","2000", "2001","2002","2004","2005","2006", "2007","2008","2009","2011","2012"), 
                 xlab = c("Year"), las = 2,                  
                 ylab = expression("Standard Ellipse Area " ('\\u2030' ^2) ),
                 bty = "o",
                 cex.axis =0.85,
                 cex.lab =5,
                 font.lab=2,
                 clr = matrix(palette(viridis(3)), nrow = 3, ncol = 16),
                 ylims = c(0,1.3))
#RUN above again to get the colours correct
# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B.y), group.ML.y[3,], col="red", pch = "x", lwd = 3, cex=1, las=3)

##### NO ROS YEARS COMPARED TO ROS YEARS AT BEGIN AND END STUDY PERIOD
##TOO many ellipses to do meaningful calculations of overlap between the different ellipses 
#Thus focusing on the years related to NO ROS and Heavy ROS

#Subset so that the years that had no ROS (1995, 2011) and were followed by a year with extreme ROS (1996, 2012),
names(mydatay)

mydatas <- filter (mydatay, group %in% c("1995", "1996", "2011", "2012"))

siber.sub <- createSiberObject(mydatas)

# BASIC PLOT

par(mfrow=c(1,1))
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))

legend2 <- c("1995", "1996","2011","2012")
colours2 <- palette(viridis(4))
pchs2 <- c(1,2,3,4)


community.hulls.args <- list(col = 1, lty = 1, lwd = 3)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

plotSiberObject(siber.sub,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                cex = 5,
                cex.lab = 1.5,
                cex.axis = 1.5,
                mai = c(1,1,1,1),
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\\u2030'),
                ylab=expression({delta}^15*N~'\\u2030'),
                x.limits = c(-26.5,-23.5),
                y.limits = c(1.5,6.5),
                main = "SIBER ellipses: 1995, 1996, 2011 & 2012")


legend("topright", legend2, col = colours2, pch = pchs2, bty = "n",cex =1, y.intersp=0.5)


# SUMMARY STATISTICS: **TA, SEA and SEAc**

group.ML.s <- groupMetricsML(siber.sub)
print(group.ML.s)

#1.1995    1.1996    1.2011    1.2012
#TA   0.6641000 0.9110000 0.8188000 0.8400500
#SEA  0.2444301 0.4625143 0.3055910 0.3150710
#SEAc 0.2632325 0.5203286 0.3259637 0.3316537

##CALCULATING THE AREA OF OVERLAP BETWEEN TWO ELLIPSES_VERY LONG PIECE OF CODE!!!CHUNK 1 detailed there after simplified into one chunk as numerous overlap possibilities with 16 years of data

#Define the ellipses
#Ellipse 1 (1995)= community 1, group 1995
ellipse1 <- "1.1995" 

# Ellipse 2 (1996): community 1, group 1996
ellipse2 <- "1.1996"

# Ellipse 3 (2011): community 1, group 2011
ellipse3 <- "1.2011"

# Ellipse 4 (2012): community 1, group 2012
ellipse4 <- "1.2012"


##COMPARING 1995 AND 1996
#To calculate the overlap between ellipses for groups 1 and 2 in community 1 and 2 respectively.
#The first ellipse is referenced using a character string representation where 
#in "x.y", "x" is the community, and "y" is the group within that community

# The overlap of the maximum likelihood fitted standard ellipses are estimated using
sea.overlap1.2 <- maxLikOverlap(ellipse1, ellipse2, siber.sub, 
                                p.interval = NULL, n = 100)

# the overlap betweeen the corresponding 95% prediction ellipses is given by:
ellipse95.overlap1.2 <- maxLikOverlap(ellipse1, ellipse2, siber.sub, 
                                      p.interval = 0.95, n = 100)

# so in this case, the overlap as a proportion of the non-overlapping area of 
# the two ellipses, would be
prop.95.over1.2 <- ellipse95.overlap1.2[3] / (ellipse95.overlap1.2[2] + 
                                                ellipse95.overlap1.2[1] -
                                                ellipse95.overlap1.2[3])
##COMPARING the 1995 to 2011

#calculate the overlaps

sea.overlap1.3 <- maxLikOverlap(ellipse1, ellipse3, siber.sub, 
                                 p.interval = NULL, n = 100)
ellipse95.overlap1.3 <- maxLikOverlap(ellipse1, ellipse3, siber.sub, 
                                       p.interval = 0.95, n = 100)
prop.95.over1.3 <- ellipse95.overlap1.3[3] / (ellipse95.overlap1.3[2] + 
                                                  ellipse95.overlap1.3[1] -
                                                  ellipse95.overlap1.3[3])

##COMPARING the 1995 to 2012

#calculate the overlaps

sea.overlap1.4 <- maxLikOverlap(ellipse1, ellipse4, siber.sub, 
                                p.interval = NULL, n = 100)
ellipse95.overlap1.4 <- maxLikOverlap(ellipse1, ellipse4, siber.sub, 
                                      p.interval = 0.95, n = 100)
prop.95.over1.4 <- ellipse95.overlap1.4[3] / (ellipse95.overlap1.4[2] + 
                                                ellipse95.overlap1.4[1] -
                                                ellipse95.overlap1.4[3])

##COMPARING the 1996 to 2011

#calculate the overlaps

sea.overlap2.3 <- maxLikOverlap(ellipse2, ellipse3, siber.sub, 
                                p.interval = NULL, n = 100)
ellipse95.overlap2.3 <- maxLikOverlap(ellipse2, ellipse3, siber.sub, 
                                      p.interval = 0.95, n = 100)
prop.95.over2.3 <- ellipse95.overlap2.3[3] / (ellipse95.overlap2.3[2] + 
                                                ellipse95.overlap2.3[1] -
                                                ellipse95.overlap2.3[3])

##COMPARING the 1996 to 2012

#calculate the overlaps

sea.overlap2.4 <- maxLikOverlap(ellipse2, ellipse4, siber.sub, 
                                p.interval = NULL, n = 100)
ellipse95.overlap2.4 <- maxLikOverlap(ellipse2, ellipse4, siber.sub, 
                                      p.interval = 0.95, n = 100)
prop.95.over2.4 <- ellipse95.overlap2.4[3] / (ellipse95.overlap2.4[2] + 
                                                ellipse95.overlap2.4[1] -
                                                ellipse95.overlap2.4[3])

##COMPARING the 2011 to 2012

#calculate the overlaps

sea.overlap3.4 <- maxLikOverlap(ellipse3, ellipse4, siber.sub, 
                                p.interval = NULL, n = 100)
ellipse95.overlap3.4 <- maxLikOverlap(ellipse3, ellipse4, siber.sub, 
                                      p.interval = 0.95, n = 100)
prop.95.over3.4 <- ellipse95.overlap3.4[3] / (ellipse95.overlap3.4[2] + 
                                                ellipse95.overlap3.4[1] -
                                                ellipse95.overlap3.4[3])

#All overlaps:

prop.95.over1.2 #0.002220698  = 0.2% 
prop.95.over1.3 #0.004224775  = 0.4%
prop.95.over1.4 #0.005236293  = 0.5%
prop.95.over2.3 #0.6027929    = 60%
prop.95.over2.4 #0.5728845    = 57%
prop.95.over3.4 #0.7086506    = 71%

### FITTING BAYESIAN MODELS TO THE DATA

ellipses.posterior.s <- siberMVN(siber.sub, parms, priors)

### COMPARING AMONG GROUPS USING STANDARD ELLIPSE AREA 

#The posterior estimates of the ellipses for each group can be used to
#calculate the SEA.B for each group.

SEA.B.s <- siberEllipses(ellipses.posterior.s)


siberDensityPlot(SEA.B.s, 
                 xticklabels = c("1995", "1996","2011","2012"), 
                 xlab = c("Year"), las = 3,
                 ylab = expression("Standard Ellipse Area " ('\\u2030' ^2) ),
                 bty = "L",
                 cex.axis =0.75,
                 clr = matrix(palette(viridis(3)), nrow = 3, ncol = 16),
                 ylims = c(0,1.5),
                 main = "SEA Bayseian Ellipse Area of Svalbard Reindeer: 1995, 1996, 2011, 2012")
#Run again to get correct colours
# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B.y), group.ML.y[3,], col="red", pch = "x", lwd = 3, cex=1, las=3)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.s.credibles <- lapply(
  as.data.frame(SEA.B.s), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.s.modes <- lapply(
  as.data.frame(SEA.B.s), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

SEA.B.s.modes
#$V1 0.2349247
#$V2 0.4186403
#$V3 0.2900127
#$V4 0.303053

summary(SEA.B.s[,1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1103  0.2147  0.2547  0.2662  0.3043  0.7711 
summary(SEA.B.s[,2])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1827  0.3950  0.4856  0.5208  0.6079  2.4792 
summary(SEA.B.s[,3])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1417  0.2666  0.3126  0.3265  0.3710  1.0530 
summary(SEA.B.s[,4])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1503  0.2778  0.3207  0.3317  0.3735  0.7557 

##Compare two ellipses for significant differences in SEA_B

#Relevant differences in SEAB are expected to be reflected by a PP ??? 0.95

# to test whether Group 1 SEA is larger than Group 2...
# you need to calculate the proportion of G1 ellipses that are greater
# than G2

#prop G1 1995 > 1996, 2011, 2012

Pg1.gt.g2 <- sum( SEA.B.s[,1] > SEA.B.s[,2] ) / nrow(SEA.B.s)
Pg1.gt.g3 <- sum( SEA.B.s[,1] > SEA.B.s[,3] ) / nrow(SEA.B.s)
Pg1.gt.g4 <- sum( SEA.B.s[,1] > SEA.B.s[,4] ) / nrow(SEA.B.s)


Pg1.gt.g2 # not significant although only just above 0.05
Pg1.gt.g3 # not significant 
Pg1.gt.g4 # not significant 


#prop G2 1996 > 1995, 2011,2012

Pg2.gt.g1 <- sum( SEA.B.s[,2] > SEA.B.s[,1] ) / nrow(SEA.B.s)
Pg2.gt.g3 <- sum( SEA.B.s[,2] > SEA.B.s[,3] ) / nrow(SEA.B.s)
Pg2.gt.g4 <- sum( SEA.B.s[,2] > SEA.B.s[,4] ) / nrow(SEA.B.s)

Pg2.gt.g1 # not significant although only just below 0.95
Pg2.gt.g3 # not significant 
Pg2.gt.g4 # not significant 


#prop G3 2011 > 1995, 1996, 2012

Pg3.gt.g1 <- sum( SEA.B.s[,3] > SEA.B.s[,1] ) / nrow(SEA.B.s)
Pg3.gt.g2 <- sum( SEA.B.s[,3] > SEA.B.s[,2] ) / nrow(SEA.B.s)
Pg3.gt.g4 <- sum( SEA.B.s[,3] > SEA.B.s[,4] ) / nrow(SEA.B.s)


Pg3.gt.g1 # not significant 
Pg3.gt.g2 # not significant 
Pg3.gt.g4 # not significant 


#prop G4 2012 > 1995, 1996, 2011

Pg4.gt.g1 <- sum( SEA.B.s[,4] > SEA.B.s[,1] ) / nrow(SEA.B.s)
Pg4.gt.g2 <- sum( SEA.B.s[,4] > SEA.B.s[,2] ) / nrow(SEA.B.s)
Pg4.gt.g3 <- sum( SEA.B.s[,4] > SEA.B.s[,3] ) / nrow(SEA.B.s)


Pg4.gt.g1 # not significant
Pg4.gt.g2 # not significant
Pg4.gt.g3 # not significant

