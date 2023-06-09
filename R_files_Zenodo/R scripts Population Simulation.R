########################################################################
## Van Irsel et al. 2021. State-dependent environmental sensitivity of reproductive success and survival in a shorebird. Ibis.
## 
## Authors: Jurrian van Irsel & Andrew M. Allen
## Date: 13-12-2021
##
## Title: R-Script for running the population simulation and averaging the multi-state mark-capture-recapture models
##
######################################################################
##
## 1. POPULATION SIMULATION
##  1.1. Averaging model parameter estimates
##  1.2. Run population simulation
##
######################################################################
library(data.table) #Model averaging, version 1.14.0
library(gdata) #Reduce function, version 2.18.0
library(MuMIn) #Version 1.43.17
library(popbio) #Population simulation version 2.7
library(RMark) #Run model.table() and view the results of the top performing multi-state live-dead recovery capture-mark-recapture models

######################################################################
##
## 1. POPULATION SIMULATION
##
######################################################################
##
## 1.1. Averaging model parameter estimates
##
#####################################################################
## Get data
covars <- read.csv('Covariates.csv')
load(TopModelResults.Rdata) ##This is saved .RData object containing the results of the RMark model comparisons
model.table(TopModelResults) ##Check results

#####################################################################
## Extract the model estimates of the top performing models
results1 <- TopModelResults[[1]]
results2 <- TopModelResults[[2]]
results3 <- TopModelResults[[3]]
results4 <- TopModelResults[[4]]
results5 <- TopModelResults[[5]]
results6 <- TopModelResults[[6]]
results7 <- TopModelResults[[7]]
results8 <- TopModelResults[[8]]
results9 <- TopModelResults[[9]]
results10 <- TopModelResults[[10]]
results11 <- TopModelResults[[11]]
results12 <- TopModelResults[[12]]
results13 <- TopModelResults[[13]]
results14 <- TopModelResults[[14]]

## Create function to calculate the real estimates
invlogit <- function(x){exp(x)/(1+exp(x))}
invlogit(4)

## Retrieve beta estimates of each of the top performing models and store into vector
Betas1 <- cbind(data.frame(Row.Names = rownames(results1$results$beta)), do.call(cbind,results1$results$beta))
Betas2 <- cbind(data.frame(Row.Names = rownames(results2$results$beta)), do.call(cbind,results2$results$beta))
Betas3 <- cbind(data.frame(Row.Names = rownames(results3$results$beta)), do.call(cbind,results3$results$beta))
Betas4 <- cbind(data.frame(Row.Names = rownames(results4$results$beta)), do.call(cbind,results4$results$beta))
Betas5 <- cbind(data.frame(Row.Names = rownames(results5$results$beta)), do.call(cbind,results5$results$beta))
Betas6 <- cbind(data.frame(Row.Names = rownames(results6$results$beta)), do.call(cbind,results6$results$beta))
Betas7 <- cbind(data.frame(Row.Names = rownames(results7$results$beta)), do.call(cbind,results7$results$beta))
Betas8 <- cbind(data.frame(Row.Names = rownames(results8$results$beta)), do.call(cbind,results8$results$beta))
Betas9 <- cbind(data.frame(Row.Names = rownames(results9$results$beta)), do.call(cbind,results9$results$beta))
Betas10 <- cbind(data.frame(Row.Names = rownames(results10$results$beta)), do.call(cbind,results10$results$beta))
Betas11 <- cbind(data.frame(Row.Names = rownames(results11$results$beta)), do.call(cbind,results11$results$beta))
Betas12 <- cbind(data.frame(Row.Names = rownames(results12$results$beta)), do.call(cbind,results12$results$beta))
Betas13 <- cbind(data.frame(Row.Names = rownames(results13$results$beta)), do.call(cbind,results13$results$beta))
Betas14 <- cbind(data.frame(Row.Names = rownames(results14$results$beta)), do.call(cbind,results14$results$beta))

head(Betas13,30) #Check results

####################
## Get beta estimates and standard errors of the estimates of each model. These will be used to average the beta estimates.
##Retrieve only the beta estimates of each model
BetasEst <- Reduce(function(x, y) merge(x, y, by = "Row.Names", all=TRUE), list(
  Betas1[1:2],
  Betas2[1:2],
  Betas3[1:2],
  Betas4[1:2],
  Betas5[1:2],
  Betas6[1:2],
  Betas7[1:2],
  Betas8[1:2],
  Betas9[1:2],
  Betas10[1:2],
  Betas11[1:2],
  Betas12[1:2],
  Betas13[1:2],
  Betas14[1:2]
))

## Separately store the beta estimates for the survival and transition parameter
BetasEstS <- BetasEst[substr(BetasEst$Row.Names, 1, 1)=="S",]
BetasEstS
BetasEstPsi <- BetasEst[substr(BetasEst$Row.Names, 1, 3)=="Psi",]
BetasEstPsi

## Retrieve only the standard errors of the beta estimates of each model
BetasSE <- Reduce(function(x, y) merge(x, y, by = "Row.Names", all=TRUE), list(
  Betas1[,c(1,3)],
  Betas2[,c(1,3)],
  Betas3[,c(1,3)],
  Betas4[,c(1,3)],
  Betas5[,c(1,3)],
  Betas6[,c(1,3)],
  Betas7[,c(1,3)],
  Betas8[,c(1,3)],
  Betas9[,c(1,3)],
  Betas10[,c(1,3)],
  Betas11[,c(1,3)],
  Betas12[,c(1,3)],
  Betas13[,c(1,3)],
  Betas14[,c(1,3)]
))
######################################################################
######################################################################

## Separately store the beta standard errors for the survival and transition parameter
BetasSES <- BetasSE[substr(BetasSE$Row.Names, 1, 1)=="S",]
BetasSEPsi <- BetasSE[substr(BetasSE$Row.Names, 1, 3)=="Psi",]

##Save all the files for the beta estimates and standard errors
write.csv(BetasEstS, "Survival BetaEst.csv")
write.csv(BetaSES, "Survival BetasSE.csv")
write.csv(BetasEstPsi, "Transition BetaEst.csv")
write.csv(BetasSEPsi, "Transitionl BetaSE.csv")

## We also need to save the model weights, obtained from the model.table. 
## Note that the above betas are sorted by model number, whereas the model table is sorted by performance. 
## It however includes the model number though so can be resorted
write.csv(model.table(TopModelResults), "Survival Beta weights.csv") 

##The next step was performed in Excel. The reason is that a lot of code needs to be written that can easily be done in excel using formulas
##We also use full-model averaging, that means all potential parameters (from the top models) need to be included in the model averaging. 
##Zeros were included if a parameter was not in the specified model (conditional averaging would include NAs here)
##The above files (Est, SE, weights) were combined for survival and transition
##And then using a type of concatenate formula, we joined all estimates into a single cell for the betas, standard errors and model weights
##These are then supplied to the par.avg function, whereby par.avg(c(beta estimates), c(standard errors), c(model weights))
##As seen below, there are 14 models to consider, which is why this is not written out by hand for each parameter. 
##Copies of the completed transition and survival files are included with the data archive to show how the par.avg() functions were prepared

## The below is based on modifications to the survival beta file. Betas, SEs and weights were combined to generate each line
SurvParams <- data.frame(t(data.frame(
  IntF= par.avg(c(2.6968468,2.6875071,2.689359,2.6395667,2.6915029,2.6898035,2.6896847,2.6792522,2.6401154,2.6458811,2.6809465,2.6882345,2.6806129,2.6793965),c(0.1100669,0.1092434,0.1107137,0.1083874,0.1108186,0.1112973,0.1112649,0.1098002,0.1099947,0.109075,0.1090669,0.1106366,0.1089588,0.1096103),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  StrN= par.avg(c(-0.4627831,-0.3338978,-0.3914192,-0.2855307,-0.4508018,-0.3940252,-0.383331,-0.349794,-0.4440229,-0.299601,-0.2971632,-0.4360327,-0.2919099,-0.4139326),c(0.1490305,0.1569841,0.1597976,0.1596052,0.1499101,0.1537888,0.1604784,0.1613995,0.148032,0.1601675,0.1615599,0.1515177,0.1604635,0.1506276),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  StrS= par.avg(c(0.3674262,0.3852633,0.4026704,0.4136475,0.3878667,0.3918634,0.4050135,0.4049438,0.4203465,0.4074912,0.3808558,0.4247344,0.3853947,0.38974),c(0.1908318,0.1935688,0.1988517,0.1919671,0.1936867,0.195788,0.1998386,0.1980644,0.1925322,0.1925065,0.1931207,0.1997429,0.1929969,0.1923777),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Fcoc= par.avg(c(0,0.1410297,0.0788918,0.1698476,0,0,0.0473743,0.1091324,0,0.19643,0.1721142,0,0.1692964,0),c(0,0.11158,0.1248658,0.1152935,0,0,0.1476252,0.1293403,0,0.1204227,0.1184742,0,0.1160781,0),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Ncoc= par.avg(c(0,-0.2621593,-0.1257465,-0.4014368,0,0,-0.022534,-0.1376857,0,-0.406794,-0.353471,0,-0.3032299,0),c(0,0.0878594,0.1152552,0.1099319,0,0,0.1276833,0.1179105,0,0.108238,0.0954279,0,0.0948087,0),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Scoc= par.avg(c(0,0.4372602,0.4134883,0.4002495,0,0,0.4215742,0.4375174,0,0.3975633,0.4543306,0,0.4617379,0),c(0,0.2054061,0.2204701,0.2015704,0,0,0.27,0.2277054,0,0.2118302,0.2156705,0,0.2117392,0),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Fmus= par.avg(c(0,0,-0.2115198,0,-0.2367677,-0.2772703,-0.2553145,-0.196281,-0.2378118,0,0,-0.2399344,0,-0.2344563),c(0,0,0.1141246,0,0.1083195,0.1187481,0.1349149,0.1179114,0.1072816,0,0,0.1176423,0,0.1107057),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Nmus= par.avg(c(0,0,0.2089502,0,0.2956289,0.3223845,0.3123415,0.2593165,0.3850355,0,0,0.315009,0,0.3396095),c(0,0,0.1295104,0,0.0993716,0.0937773,0.1296162,0.1244964,0.1084307,0,0,0.1016386,0,0.0994139),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Smus= par.avg(c(0,0,-0.1540241,0,-0.2463273,-0.3234036,-0.1681671,-0.1555037,-0.2213313,0,0,-0.3687773,0,-0.2642859),c(0,0,0.1808381,0,0.1613294,0.1860875,0.2095144,0.183402,0.1597802,0,0,0.1976125,0,0.1653075),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Fwind= par.avg(c(0,0,0,0.1353133,0,0,0,0,0.1416501,0.1413779,0,0,0,0),c(0,0,0,0.1303597,0,0,0,0,0.1280177,0.1331029,0,0,0,0),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Nwind= par.avg(c(0,0,0,0.418024,0,0,0,0,0.3674065,0.4194918,0,0,0,0),c(0,0,0,0.1328643,0,0,0,0,0.1276594,0.1301699,0,0,0,0),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Swind= par.avg(c(0,0,0,-0.1558765,0,0,0,0,-0.1952866,-0.1543033,0,0,0,0),c(0,0,0,0.181145,0,0,0,0,0.1743686,0.1818016,0,0,0,0),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Fbwl= par.avg(c(0,0,0,0,0,-0.0889491,-0.0651119,0,0,0,0,0,0,0),c(0,0,0,0,0,0.1198372,0.1405371,0,0,0,0,0,0,0),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Nbwl= par.avg(c(0,0,0,0,0,0.2381677,0.2335923,0,0,0,0,0,0,0),c(0,0,0,0,0,0.1028158,0.1136517,0,0,0,0,0,0,0),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Sbwl= par.avg(c(0,0,0,0,0,-0.1898713,-0.0057548,0,0,0,0,0,0,0),c(0,0,0,0,0,0.1627345,0.2075434,0,0,0,0,0,0,0),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Fwsi= par.avg(c(0,0,0,0,0,0,0,0.0891763,0,0,0.1281908,0,0.1482563,0.0616769),c(0,0,0,0,0,0,0,0.1207644,0,0,0.152795,0,0.1175184,0.1158346),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Nwsi= par.avg(c(0,0,0,0,0,0,0,0.2564384,0,0,0.4977023,0,0.1970188,0.2469941),c(0,0,0,0,0,0,0,0.1101816,0,0,0.1655756,0,0.1048928,0.1097391),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Swsi= par.avg(c(0,0,0,0,0,0,0,0.0296589,0,0,0.0866421,0,0.0582492,-0.0560517),c(0,0,0,0,0,0,0,0.169442,0,0,0.2138599,0,0.1681667,0.1525045),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Fprecip= par.avg(c(0,0,0,0,0,0,0,0,0,-0.1099073,-0.0263461,0.0046332,0,0),c(0,0,0,0,0,0,0,0,0,0.1016672,0.1357655,0.1117495,0,0),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Nprecip= par.avg(c(0,0,0,0,0,0,0,0,0,0.0855759,0.4027834,-0.0657579,0,0),c(0,0,0,0,0,0,0,0,0,0.1096771,0.1538181,0.1076518,0,0),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965)),
  Sprecip= par.avg(c(0,0,0,0,0,0,0,0,0,0.0051449,0.0434204,0.2289351,0,0),c(0,0,0,0,0,0,0,0,0,0.1514301,0.1945763,0.1696588,0,0),c(0.000233565799173199,0.0191964365104738,0.0255923292930968,0.314911822037176,0.0402546153107005,0.0649938512714137,0.0110773551862382,0.0316235770390362,0.374387991310169,0.0361581709179138,0.024444444009491,0.0055140129964886,0.0145238578985325,0.0370879704200965))
)))
SurvParams$row <- seq(1,nrow(SurvParams),1)

##Check that the relationships are roughly correct (there will always be deviations when averaging betas (what we do) and averaging reals (what RMark does)
PlotMUScheck <- data.frame(Cov = rep(seq(min(covariates$musselschier), max(covariates$musselschier), 0.1),3),State = rep(rep(c("F", "S", "N"),each = 39)))
PlotMUScheck$Surv <- NA
PlotMUScheck$Surv[1:39] <- invlogit(SurvParams[1,1] + SurvParams[7,1]*PlotMUScheck$Cov[1:39])
PlotMUScheck$Surv[40:78] <- invlogit(SurvParams[1,1] + SurvParams[3,1] + SurvParams[9,1]*PlotMUScheck$Cov[40:78])
PlotMUScheck$Surv[79:117] <- invlogit(SurvParams[1,1] + SurvParams[2,1] +SurvParams[8,1]*PlotMUScheck$Cov[79:117])

##The below is based on modifications to the Transition beta file. Betas, SEs and weights were combined to generate each line
TransitionParams <- data.frame(t(data.frame(
  FN= par.avg(c(-0.3906047,-0.3947338,-0.3937969,-0.4015675,-0.3921278,-0.3935666,-0.3940421,-0.3973544,-0.398596,-0.4014913,-0.3976457,-0.3929639,-0.3980877,-0.3954008),c(0.0609476,0.0609614,0.0609685,0.0610315,0.0609567,0.0609603,0.0609681,0.0609771,0.0610618,0.0610358,0.0609927,0.0609809,0.0609801,0.060971),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  FNrag= par.avg(c(-0.0147261,-0.0134376,-0.0149682,-0.0125992,-0.0164289,-0.0178984,-0.0175776,-0.0170793,-0.0149304,-0.0116463,-0.0161092,-0.0165739,-0.0142855,-0.0185761),c(0.0624014,0.0624419,0.0624412,0.0624604,0.0624129,0.062434,0.0624535,0.0624724,0.0624443,0.0624474,0.0624632,0.0624354,0.0624633,0.0624483),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  FNtemp= par.avg(c(0.3806204,0.3832492,0.3808878,0.3803413,0.3790523,0.3807376,0.3807352,0.3797198,0.3746173,0.3786592,0.3770827,0.3800761,0.3808446,0.3785103),c(0.0713035,0.0713694,0.0713604,0.071578,0.0713211,0.0713386,0.0713605,0.0714365,0.0715492,0.0715904,0.0715195,0.071345,0.0714681,0.0713959),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  FNtide= par.avg(c(-0.1752332,-0.1808855,-0.1835488,-0.1805908,-0.1818787,-0.1863423,-0.1864301,-0.1859502,-0.1805691,-0.1805788,-0.1826587,-0.1828292,-0.1825562,-0.1843037),c(0.0638122,0.0638833,0.0638902,0.0639766,0.0638607,0.0639021,0.063909,0.0639473,0.0639563,0.0639935,0.0639694,0.0638798,0.0639424,0.0639144),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  FS= par.avg(c(-0.5959613,-0.5959551,-0.5959396,-0.5959666,-0.5959354,-0.5959205,-0.5959246,-0.595925,-0.5959448,-0.5959765,-0.5959561,-0.5959293,-0.5959473,-0.5959185),c(0.0647739,0.0647716,0.0647708,0.064772,0.0647718,0.0647704,0.06477,0.064769,0.0647726,0.064772,0.0647698,0.0647714,0.0647699,0.0647704),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  FSrag= par.avg(c(0.2317742,0.2320025,0.2317127,0.2319606,0.2315763,0.2314853,0.2315383,0.2316781,0.231511,0.2320455,0.2318322,0.2315814,0.2320632,0.2315017),c(0.0621957,0.0622192,0.0621838,0.0622153,0.0621687,0.0621554,0.0621617,0.0621771,0.0621621,0.0622262,0.0621985,0.0621686,0.0622242,0.0621573),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  FStemp= par.avg(c(0.0154291,0.0155136,0.0154771,0.0153994,0.0154583,0.015482,0.0154788,0.0155115,0.0153554,0.0153879,0.0154233,0.0154826,0.0155537,0.0154924),c(0.0690478,0.0690057,0.0690054,0.0689876,0.0690271,0.069009,0.0690005,0.0689815,0.069018,0.0689855,0.0689737,0.0690232,0.0689816,0.0690107),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  FStide= par.avg(c(0.0213283,0.0210253,0.0211742,0.0208248,0.021315,0.021269,0.0212196,0.0210323,0.0211859,0.020798,0.0209858,0.0212722,0.0208567,0.0212172),c(0.0653877,0.0653487,0.0653703,0.0653065,0.0653872,0.0653849,0.0653802,0.065362,0.0653478,0.0653158,0.0653325,0.0653858,0.0653446,0.0653834),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  NF= par.avg(c(-1.3959847,-1.3980755,-1.3962776,-1.386081,-1.3944662,-1.3956458,-1.3958131,-1.3934282,-1.3825607,-1.3858231,-1.3932934,-1.3946987,-1.3956944,-1.3918244),c(0.0679296,0.0679249,0.0679485,0.0681626,0.0679338,0.0678951,0.0679004,0.0679256,0.0681574,0.0681721,0.0679578,0.0679251,0.067904,0.0679187),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  NFrag= par.avg(c(-0.1927582,-0.1927459,-0.1908844,-0.1875422,-0.1898518,-0.1859128,-0.1861459,-0.1853241,-0.1850663,-0.1883575,-0.1856808,-0.1891799,-0.1887848,-0.1846498),c(0.076502,0.0765997,0.0765731,0.0763103,0.0765226,0.076569,0.0765933,0.0765414,0.0761886,0.0762983,0.0765003,0.0765404,0.076591,0.0764842),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  NFtemp= par.avg(c(-0.1442515,-0.1475401,-0.1443429,-0.1277072,-0.1415131,-0.1398656,-0.1400449,-0.1388839,-0.120867,-0.1265958,-0.1338552,-0.1424886,-0.1439156,-0.1364685),c(0.0670483,0.0670456,0.0671575,0.0679244,0.0671591,0.0670785,0.0671338,0.0672185,0.0680664,0.0679613,0.0673486,0.0671494,0.0670842,0.0672175),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  NFtide= par.avg(c(0.1244626,0.1283775,0.12541,0.121066,0.1224192,0.1261942,0.1264878,0.1255883,0.1136105,0.1202765,0.1261799,0.123231,0.1287526,0.1227661),c(0.0634438,0.0634333,0.0634569,0.0635167,0.0633989,0.0634263,0.0634427,0.0634116,0.0634825,0.0635257,0.0633951,0.0634074,0.0634131,0.0633722),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  NS= par.avg(c(-2.2427763,-2.2447661,-2.2431396,-2.234937,-2.2415324,-2.2426213,-2.2427542,-2.2406642,-2.2318959,-2.2347812,-2.2410032,-2.2416927,-2.2426126,-2.2392388),c(0.0996446,0.0996386,0.0996581,0.09988,0.0996603,0.0996342,0.0996351,0.0996655,0.0998979,0.0998903,0.099712,0.0996515,0.0996397,0.0996701),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  NSrag= par.avg(c(0.2856828,0.2866804,0.2882403,0.2900807,0.2889166,0.293294,0.2931064,0.2935992,0.2914843,0.2890982,0.2930736,0.2896487,0.2906388,0.2938119),c(0.0862556,0.0863185,0.0863135,0.0862014,0.08629,0.0863194,0.0863315,0.0863086,0.086141,0.0862029,0.0862929,0.0863035,0.0863251,0.0862787),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  NStemp= par.avg(c(-0.4258637,-0.4288654,-0.4259049,-0.4129552,-0.4234627,-0.4216787,-0.4217964,-0.4211374,-0.4069637,-0.4119804,-0.416896,-0.4243226,-0.4256721,-0.4190645),c(0.0878373,0.0878434,0.0879131,0.0883423,0.0879184,0.0878505,0.0878863,0.0879418,0.0884277,0.088359,0.0879892,0.0879138,0.0878643,0.0879425),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  NStide= par.avg(c(0.0840924,0.0876261,0.0853278,0.0857763,0.0828868,0.0862596,0.0864539,0.086332,0.0798557,0.0853992,0.0882144,0.0834771,0.0884135,0.0840482),c(0.0857956,0.085764,0.0857362,0.0855961,0.0856951,0.0857468,0.0857523,0.0856755,0.0855237,0.0855704,0.0856079,0.0857103,0.0857315,0.085651),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  SF= par.avg(c(0.0211091,0.0210916,0.021065,0.0210686,0.0210807,0.0210494,0.021047,0.0210426,0.0210493,0.0210689,0.0210614,0.0210371,0.0210797,0.0210518),c(0.0703232,0.0703217,0.0703211,0.0703209,0.0703222,0.0703212,0.0703207,0.0703203,0.0703211,0.0703208,0.0703206,0.0703209,0.0703212,0.0703213),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  SFrag= par.avg(c(-0.3130129,-0.3129124,-0.3127454,-0.3127107,-0.3128312,-0.3126152,-0.3126364,-0.3126183,-0.3125627,-0.3126991,-0.3126596,-0.3124507,-0.3128614,-0.3126229),c(0.0755547,0.07554,0.0755198,0.0755164,0.0755318,0.0755069,0.0755061,0.0755027,0.0755014,0.075515,0.075508,0.0754878,0.0755319,0.0755068),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  SFtemp= par.avg(c(0.0905312,0.0903355,0.0903262,0.0902609,0.0904558,0.090329,0.0903037,0.0902628,0.090353,0.0902718,0.090305,0.090372,0.0902967,0.0903782),c(0.0742099,0.0741524,0.0741442,0.0741378,0.0741886,0.0741665,0.074133,0.0741261,0.0741695,0.0741375,0.0741227,0.0741704,0.0741377,0.0741729),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  SFtide= par.avg(c(0.0387709,0.0387296,0.038747,0.0386317,0.0387786,0.0387201,0.0387535,0.0387576,0.0386447,0.038611,0.0386447,0.0386191,0.0387634,0.0387323),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  SN= par.avg(c(-0.8480706,-0.8533916,-0.8529526,-0.860686,-0.8507211,-0.8537043,-0.8536142,-0.8568194,-0.8587199,-0.8603503,-0.8577856,-0.8514542,-0.856568,-0.8548683),c(0.0953584,0.0953585,0.0953485,0.0954144,0.0953727,0.0953827,0.0953194,0.0953467,0.0954575,0.0954171,0.0953676,0.0953588,0.0953778,0.095393),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  SNrag= par.avg(c(-0.2443454,-0.2461538,-0.2463897,-0.2462579,-0.2459378,-0.247483,-0.2479159,-0.2482753,-0.2458263,-0.2459696,-0.2485297,-0.2473614,-0.2474948,-0.2479255),c(0.0868563,0.0868972,0.086901,0.0869726,0.0868885,0.0869172,0.0869122,0.0869486,0.0869867,0.0869712,0.0869616,0.0869135,0.0869374,0.086947),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  SNtemp= par.avg(c(0.4256774,0.4266242,0.4241805,0.4283042,0.4233947,0.4247969,0.4232214,0.4242847,0.4254609,0.4275036,0.4225374,0.4232594,0.4265523,0.4245531),c(0.1063045,0.1065245,0.1065225,0.1068552,0.1064037,0.1065161,0.1064912,0.1065998,0.1067761,0.106868,0.106697,0.1065092,0.1066126,0.1065298),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461)),
  SNtide= par.avg(c(-0.4504373,-0.4585458,-0.4604058,-0.4603689,-0.4569414,-0.4619484,-0.4626398,-0.463929,-0.4587357,-0.4600314,-0.4630441,-0.459999,-0.4612621,-0.4606461),c(0.0931772,0.0934224,0.0933959,0.0935558,0.0932457,0.0933457,0.0933965,0.0934888,0.093396,0.0935627,0.0935777,0.0933378,0.0935045,0.0933397),c(0.0674594,0.067417,0.0674258,0.0674001,0.0674546,0.0674429,0.0674301,0.067417,0.0674333,0.0673982,0.0673967,0.0674365,0.0674115,0.0674461))
)))
TransitionParams$row <- seq(1,nrow(TransitionParams),1)

##Double check transitions are roughly correct. Some deviations expected due to averaging betas versus predicted reals. Remember it is a mlogit parameter
PlotRAGcheck <- data.frame(Cov = rep(seq(min(covariates$ragwormnext), max(covariates$ragwormnext), 0.1),2),State = rep(rep(c("FS", "NS"),each = 40)))
PlotRAGcheck$Transition <- NA
##When estimating the real estimates of transitions, remember it is a mlogit transformation
##that means the logit of transitions are confined to sum to 1, for example SS + SF + SN = 1
##and hence the below formula to obtain the reals taking into account other estimates
for(i in 1:40){
  PlotRAGcheck$Transition[i] <- (exp(TransitionParams["FS",1]+TransitionParams["FSrag",1]*PlotRAGcheck$Cov[i])	/
                                   (1+ sum(exp(TransitionParams["FS",1]+TransitionParams["FSrag",1]*PlotRAGcheck$Cov[i]),
                                           exp(TransitionParams["FN",1]+TransitionParams["FNrag",1]*PlotRAGcheck$Cov[i])))
  )
}
for(i in 41:80){
  PlotRAGcheck$Transition[i] <- (exp(TransitionParams["NS",1]+TransitionParams["NSrag",1]*PlotRAGcheck$Cov[i])	/
                                   (1+ sum(exp(TransitionParams["NS",1]+TransitionParams["NSrag",1]*PlotRAGcheck$Cov[i]),
                                           exp(TransitionParams["NF",1]+TransitionParams["NFrag",1]*PlotRAGcheck$Cov[i])))
  )}

##Combine the averaged survival and transition beta estimates
params<-rbind(SurvParams, TransitionParams)

##Prepare functions to estimate the average survival probability for each state (Successful breeder (S), Failed breeder (F) and Non-breeder (N) dependendt upon the environmental covariates. 
##The beta estimates of the effects of each environmental conditions acting upon each state are averaged over the top performing models.
surv <- function(p,windchillnext=0,cockle=0,musselschier=0,ragwormnext=0,avgtempincubationnext=0,tidalmaxnext=0) {
  survN<-plogis(p["IntF",1]+p["StrN",1]+p["Nwind",1]*windchillnext+p["Ncoc",1]*cockle+p["Nmus",1]*musselschier)
  survF<-plogis(p["IntF",1]+p["Fwind",1]*windchillnext+p["Fcoc",1]*cockle+p["Fmus",1]*musselschier)
  survS<-plogis(p["IntF",1]+p["StrS",1]+p["Swind",1]*windchillnext+p["Scoc",1]*cockle+p["Smus",1]*musselschier)
  return(c(survN,survF,survS))
}

##Prepare functions to estimate the average transition probability of each state transitioning into another state (Successful breeder (S), Failed breeder (F) and Non-breeder (N) dependendt upon the environmental covariates. 
##The beta estimates of the effects of each environmental conditions acting upon each transition are averaged over the top performing models. The average transition estimates are calculating by dividing the expontial of each environmental condition acting upon e.g. Non-breeder to Failed breeder, divided by the sum of the exponential of the environmental conditions acting upon the transition of Non-breeder to all other states.
trans <- function(p,windchillnext=0,cockle=0,musselschier=0,ragwormnext=0,avgtempincubationnext=0,tidalmaxnext=0) {
  psi<-matrix(NA,nrow=3,ncol=3)
  dimnames(psi)<-list(c("N","F","S"),c("N","F","S"))
  psi["F","N"]<-(
    exp(	p["NF",1]+p["NFrag",1]*ragwormnext+p["NFtemp",1]*avgtempincubationnext+
           p["NFtide",1]*tidalmaxnext)/
      (1+ sum(exp(p["NF",1]+p["NFrag",1]*ragwormnext+p["NFtemp",1]*avgtempincubationnext+
                    p["NFtide",1]*tidalmaxnext),
              exp(p["NS",1]+p["NSrag",1]*ragwormnext+p["NStemp",1]*avgtempincubationnext+
                    p["NStide",1]*tidalmaxnext)))
  )
  psi["S","N"]<-(
    exp(	p["NS",1]+p["NSrag",1]*ragwormnext+p["NStemp",1]*avgtempincubationnext+
           p["NStide",1]*tidalmaxnext)/
      (1+ sum(exp(p["NS",1]+p["NSrag",1]*ragwormnext+p["NStemp",1]*avgtempincubationnext+
                    p["NStide",1]*tidalmaxnext),
              exp(p["NF",1]+p["NFrag",1]*ragwormnext+p["NFtemp",1]*avgtempincubationnext+
                    p["NFtide",1]*tidalmaxnext)))
  )
  psi["N","N"]<-1-psi["F","N"]-psi["S","N"]
  
  psi["N","S"]<-(
    exp(	p["SN",1]+p["SNrag",1]*ragwormnext+p["SNtemp",1]*avgtempincubationnext+
           p["SNtide",1]*tidalmaxnext)/
      (1+ sum(exp(p["SN",1]+p["SNrag",1]*ragwormnext+p["SNtemp",1]*avgtempincubationnext+
                    p["SNtide",1]*tidalmaxnext),
              exp(p["SF",1]+p["SFrag",1]*ragwormnext+p["SFtemp",1]*avgtempincubationnext+
                    p["SFtide",1]*tidalmaxnext)))
  )
  psi["F","S"]<-(
    exp(	p["SF",1]+p["SFrag",1]*ragwormnext+p["SFtemp",1]*avgtempincubationnext+
           p["SFtide",1]*tidalmaxnext)/
      (1+ sum(exp(p["SF",1]+p["SFrag",1]*ragwormnext+p["SFtemp",1]*avgtempincubationnext+
                    p["SFtide",1]*tidalmaxnext),
              exp(p["SN",1]+p["SNrag",1]*ragwormnext+p["SNtemp",1]*avgtempincubationnext+
                    p["SNtide",1]*tidalmaxnext)))
  )
  psi["S","S"]<-1-psi["N","S"]-psi["F","S"]
  psi["S","F"]<-(
    exp(	p["FS",1]+p["FSrag",1]*ragwormnext+p["FStemp",1]*avgtempincubationnext+
           p["FStide",1]*tidalmaxnext)/
      (1+ sum(exp(p["FS",1]+p["FSrag",1]*ragwormnext+p["FStemp",1]*avgtempincubationnext+
                    p["FStide",1]*tidalmaxnext),
              exp(p["FN",1]+p["FNrag",1]*ragwormnext+p["FNtemp",1]*avgtempincubationnext+
                    p["FNtide",1]*tidalmaxnext)))
  )
  psi["N","F"]<-(
    exp(	p["FN",1]+p["FNrag",1]*ragwormnext+p["FNtemp",1]*avgtempincubationnext+
           p["FNtide",1]*tidalmaxnext)/
      (1+ sum(exp(p["FN",1]+p["FNrag",1]*ragwormnext+p["FNtemp",1]*avgtempincubationnext+
                    p["FNtide",1]*tidalmaxnext),
              exp(p["FS",1]+p["FSrag",1]*ragwormnext+p["FStemp",1]*avgtempincubationnext+
                    p["FStide",1]*tidalmaxnext)))
  )
  psi["F","F"]<-1-psi["S","F"]-psi["N","F"]
  
  return(psi)
}

s0 <- surv(p = params) #Average survival probabilities for Non-breeder, Failed- and Successful breeder respectively.
t0 <- trans(p = params) #Average transition probabilities of Non-breeder, Failed- and Successful breeder to all other states.
u0 <- t0*matrix(s0,nrow=3,byrow=T,ncol=3) #Average transition probability of Non-breeder, Failed- and Successful breeders conditional upon the survival per state.

############################################################################################################################################
##
## 1.2. RUN POPULATION SIMULATION WITH ENVIRONMENTAL COVARIATES
##
############################################################################################################################################
ssd <- round(stable.stage(u0)*500,0) # estimated population size 2002-2015

#Create function that runs 1000 simulations with fixed transition and survival probability per state using the estimates retrieved from the multi-state live-dead recovery capture-mark-recapture model
sim <- function(t, s, nsim=1000) {
  y <- 100 #number of time steps in each simulation
  book <- matrix(0,nrow=3,ncol=y) ##matrix to store population developments
  book[,1] <- round(stable.stage(t*matrix(s,nrow=3,byrow=T,ncol=3))*500,0) #first time step
  survivors <- successes <- numeric(nsim) ##objects to store survivors and successes in each sim
  for (i in 1:nsim) { ##start simulations
    nn <- book ##store results in nn
    for (j in 2:y) { ##first year is already done (i.e. book[,1]), so now 2:y which is 2:100
      for (k in 1:3) { ##number of breeding states
        if (nn[k,j-1]>0) {
          nn[,j] <- nn[,j]+rmultinom(size=rbinom(size=nn[k,j-1],prob=s[k],n=1),prob=t[,k],n=1) #survival and transitions based on model probabilities
        }
      }
    }
    survivors[i] <- sum(nn[,y]) ##assemble results
    successes[i] <- sum(nn[3,2:y]) ##assemble results
  }
  ##return average results (from the simulations)
  return(data.frame(meanSuccess=mean(successes),sdSuccesses=sd(successes),meanSurvivors=mean(survivors),sdSurvivors=sd(survivors))/sum(ssd))
}

#Set initial survival and transition probabilities for the simulation
sim(t = t0, s = s0)

#Run simulation for each covariate using the average beta estimates for survival or transition by setting the other covariates to zero.
#Remember the data is standardised, hence why 0 represent average (mean-centered to zero with 1 SD)
# Environmental covariate Common Cockle biomass
baseline<-sim(t=trans(p=params),s=surv(p=params),nsim=10000)[1]
dum<-21 #number of data points for environemntal covariates
dumCockle<-seq(from=-2,to=2,length=dum) #covariates values, from -2 SD to 2 SD
xCockle<-dumCockle
meanCockleSurv<-sdCockleSurv<-rep(NA,dum) #empty objects to fill with results

##For each simulation we calculate average effects for the parameters not of interest (i.e. transition/survival and only single environmental covariates

for (p in 1:dum) {
  print(p)
  z<-sim(t=trans(p=params),s=surv(p=params,cockle=dumCockle[p])) ##average transition with only cockle effects on survival
  meanCockleSurv[p]<-unlist(z[1])
  sdCockleSurv[p]<-unlist(z[2])
}
# Environmental covariate Winchill index
meanwindchillnextSurv<-sdwindchillnextSurv<-rep(NA,dum)
dumwindchillnext<-seq(from=-2,to=2,length=dum)
xwindchillnext<-dumwindchillnext
for (p in 1:dum) {
  print(p)
  z<-sim(t=trans(p=params),s=surv(p=params,windchillnext=dumwindchillnext[p]))
  meanwindchillnextSurv[p]<-unlist(z[1])
  sdwindchillnextSurv[p]<-unlist(z[2])
}

# Environmental covariate Blue Mussel biomass
meanmusselschierSurv<-sdmusselschierSurv<-rep(NA,dum)
dummusselschier<-seq(from=-2, to=2,length=dum)
xmusselschier<-dummusselschier
for (p in 1:dum) {
  print(p)
  z<-sim(t=trans(p=params),s=surv(p=params,musselschier=dummusselschier[p]))
  meanmusselschierSurv[p]<-unlist(z[1])
  sdmusselschierSurv[p]<-unlist(z[2])
}
# Environmental covariate Ragworm biomass
meanragwormnext<-sdragwormnext<-rep(NA,dum)
dumragwormnext<-seq(from=-2,to=2,length=dum)
xragwormnext<-dumragwormnext
for (p in 1:dum) {
  print(p)
  x<-sim(t=trans(p=params,ragwormnext=dumragwormnext[p]),s=surv(p=params))
  meanragwormnext[p]<-unlist(x[1])
  sdragwormnext[p]<-unlist(x[2])
}

# Environmental covariate Tidal Heigth
meantidalmaxnext<-sdtidalmaxnext<-rep(NA,dum)
dumtidalmaxnext<-seq(from=-2,to=2,length=dum)
xtidalmaxnext<-dumtidalmaxnext
for (p in 1:dum) {
  print(p)
  x<-sim(t=trans(p=params,tidalmaxnext=dumtidalmaxnext[p]),s=surv(p=params))
  meantidalmaxnext[p]<-unlist(x[1])
  sdtidalmaxnext[p]<-unlist(x[2])
}

# Environmental covariate average temperature during May and June
meanavgtempincubationnext<-sdavgtempincubationnext<-rep(NA,dum)
dumavgtempincubationnext<-seq(from=-2,to=2,length=dum)
xavgtempincubationnext<-dumavgtempincubationnext
for (p in 1:dum) {
  print(p)
  x<-sim(t=trans(p=params,avgtempincubationnext=dumavgtempincubationnext[p]),s=surv(p=params))
  meanavgtempincubationnext[p]<-unlist(x[1])
  sdavgtempincubationnext[p]<-unlist(x[2])
}

#Merge the simulation dataframe into a single dataframe for the environmental effects acting upon the survival probability
dataS <- melt(data.frame(meanCockleSurv, meanwindchillnextSurv, meanmusselschierSurv))
dataSD <- melt(data.frame(sdCockleSurv, sdwindchillnextSurv, sdmusselschierSurv))

dataC$MeanSuccess <- dataS$value
dataC$MeanSuccessSD <- dataSD$value

#Merge the simulation dataframe into a single dataframe for the environmental effects acting upon the transition probability
dataTr <- melt(data.frame(meanragwormnext, meanavgtempincubationnext, meantidalmaxnext))
dataTrSD <- melt(data.frame(sdragwormnext, sdavgtempincubationnext, sdtidalmaxnext))

dataR$MeanSuccess <- dataTr$value
dataR$MeanSuccessSD <- dataTrSD$value

write.csv(dataC, "LTRS winter survival.csv")
write.csv(dataR, "LTRS incubation.csv")