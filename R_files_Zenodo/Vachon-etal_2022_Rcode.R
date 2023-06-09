## Code for Vachon et al.2022: Distinctive, fine-scale distribution of Eastern Caribbean sperm whale vocal clans reflects island fidelity rather than environmental variables
 # Based on code from Pirotta et al.2011 (Modelling sperm whale habitat preference: a novel approach combining transect and follow data) and Eguiguren et al.2019 (Habitat use of culturally distinct Galápagos sperm whale Physeter macrocephalus clans)

## Uploading the data
setwd("")
dat <- read.table("", header=TRUE, sep=",")
head(dat); dim(dat)  # to visualize the structure of the data

## Upload the required libraries
library(geepack)         # for the GEEs (Wald's hypothesis tests allowed)
library(MuMIn)           # for model selection  (this also provides QIC)
library(MASS)            # for matrix
library(gee)
library(gam)
library(ROCR)            # to build the ROC curve
library(PresenceAbsence) # to build the confusion matrix
library(ggplot2)         # to build the partial residual plots
library(MRSea)           # to run ACF
library(car)             # for GVIF multicolinearity tests
library(stringr)         # for backward selection
library(mgcv)            # for qqGAM
library (beepr)          # Ping when analysis has finished running
library("plotrix")       # For stepwise cross validation
library (effects)        # For effect plots

#----------------STEP 1: Prepare variables---------------------------

# Look at variable distribution to determine if they should be logged

hist(dat$Lat)       
hist(dat$Long)      
hist(dat$Depth)
hist(dat$Slope)     #Needs to be log
hist(dat$Canyon)
hist(dat$Escarp)    
hist(dat$Abyss)
hist(dat$Shelf)
hist(dat$Ecurr)
hist(dat$Ncurr)
hist(dat$Zvelv)
hist(dat$Mvelv)
hist(dat$Inflow)
hist(dat$Channel)
hist(dat$Chla)

## Transform the variables (using LN)
# For variables that had to be log, create an additional column in the csv with log transformation

hist(log(dat$Slope))
hist(dat$logSlope)

## Center the variables
# scale subtracts the mean and divides by the standard dev (standardize the data)

Lat             <- as.vector(scale(dat$Lat))

Long       <- as.vector(scale(dat$Long))

Depth   <- as.vector(scale(dat$Depth))

logSlope          <- as.vector(scale(dat$logSlope))

Canyon               <- as.vector(scale(dat$Canyon))

Escarp               <- as.vector(scale(dat$Escarp))

Abyss               <- as.vector(scale(dat$Abyss))

Shelf               <- as.vector(scale(dat$Shelf))

Ecurr               <- as.vector(scale(dat$Ecurr))

Ncurr               <- as.vector(scale(dat$Ncurr))

Zvelv               <- as.vector(scale(dat$Zvelv))

Mvelv               <- as.vector(scale(dat$Mvelv))

Inflow               <- as.vector(scale(dat$Inflow))

Channel               <- as.vector(scale(dat$Channel))

Chla               <- as.vector(scale(dat$Chla))

Island               <- (dat$Island)

Windward               <- (dat$Windward)


#----------------STEP 2: Remove correlation----------------------------

## Calculate correlation coefficients for continuous variables
cor(Lat, Long) 
cor(Lat, Depth) 
cor(Lat, logSlope) 
cor(Lat, Canyon) 
cor(Lat, Escarp) 
cor(Lat, Abyss) 
cor(Lat, Shelf) 
cor(Lat, Ecurr) 
cor(Lat, Ncurr) 
cor(Lat, Zvelv) 
cor(Lat, Mvelv) 
cor(Lat, Inflow) 
cor(Lat, Channel) 
cor(Lat, Chla) 
cor(Long, Depth) 
cor(Long, logSlope) 
cor(Long, Canyon) 
cor(Long, Escarp) 
cor(Long, Abyss) 
cor(Long, Shelf) 
cor(Long, Ecurr) 
cor(Long, Ncurr) 
cor(Long, Zvelv) 
cor(Long, Mvelv) 
cor(Long, Inflow) 
cor(Long, Channel) 
cor(Long, Chla) 
cor(Depth, logSlope) 
cor(Depth, Canyon) 
cor(Depth, Escarp) 
cor(Depth, Abyss) 
cor(Depth, Shelf) 
cor(Depth, Ecurr) 
cor(Depth, Ncurr) 
cor(Depth, Zvelv) 
cor(Depth, Mvelv) 
cor(Depth, Inflow) 
cor(Depth, Channel) 
cor(Depth, Chla) 
cor(logSlope, Canyon) 
cor(logSlope, Escarp) 
cor(logSlope, Abyss) 
cor(logSlope, Shelf) 
cor(logSlope, Ecurr) 
cor(logSlope, Ncurr) 
cor(logSlope, Zvelv) 
cor(logSlope, Mvelv) 
cor(logSlope, Inflow) 
cor(logSlope, Channel) 
cor(logSlope, Chla) 
cor(Canyon, Escarp) 
cor(Canyon, Abyss) 
cor(Canyon, Shelf) 
cor(Canyon, Ecurr) 
cor(Canyon, Ncurr) 
cor(Canyon, Zvelv) 
cor(Canyon, Mvelv) 
cor(Canyon, Inflow) 
cor(Canyon, Channel) 
cor(Canyon, Chla) 
cor(Escarp, Abyss) 
cor(Escarp, Shelf) 
cor(Escarp, Ecurr) 
cor(Escarp, Ncurr) 
cor(Escarp, Zvelv) 
cor(Escarp, Mvelv) 
cor(Escarp, Inflow) 
cor(Escarp, Channel) 
cor(Escarp, Chla) 
cor(Abyss, Shelf) 
cor(Abyss, Ecurr) 
cor(Abyss, Ncurr) 
cor(Abyss, Zvelv) 
cor(Abyss, Mvelv) 
cor(Abyss, Inflow) 
cor(Abyss, Channel) 
cor(Abyss, Chla) 
cor(Shelf, Ecurr) 
cor(Shelf, Ncurr) 
cor(Shelf, Zvelv) 
cor(Shelf, Mvelv) 
cor(Shelf, Inflow) 
cor(Shelf, Channel) 
cor(Shelf, Chla) 
cor(Ecurr, Ncurr) 
cor(Ecurr, Zvelv) 
cor(Ecurr, Mvelv) 
cor(Ecurr, Inflow) 
cor(Ecurr, Channel) 
cor(Ecurr, Chla) 
cor(Ncurr, Zvelv) 
cor(Ncurr, Mvelv) 
cor(Ncurr, Inflow) 
cor(Ncurr, Channel) 
cor(Ncurr, Chla) 
cor(Zvelv, Mvelv) 
cor(Zvelv, Inflow) 
cor(Zvelv, Channel) 
cor(Zvelv, Chla) 
cor(Mvelv, Inflow) 
cor(Mvelv, Channel) 
cor(Mvelv, Chla) 
cor(Inflow, Channel) 
cor(Inflow, Chla) 
cor(Channel, Chla) 

#Create a matrix of variable correlations
#From this came up with potential models (all combinations of uncorrelated variables) --> done manually

##Calculate GVIF, make sure no multicollinearity in potential models
#Do for each model (discard models with GVIF > 3)

#GVIF calculation example                     
modA <-glm(Pres~Lat +
             Escarp +
             Abyss +
             Zvelv +
             Channel +
             Windward
           , family = "binomial", data=dat)

vifa <-as.data.frame(vif(modA))
vifa


#----------------STEP 3: Model selection-----------------------------

#Manual backward selection for each potential model

##Activate beepr (optional)
beep(sound = 1, expr = NULL) #Put after each run to know when they are finished

#Example (MODEL A)

#model A step 1
comb_A_s1_null  <- geeglm(Pres~ Lat + Escarp + Abyss + Zvelv + Channel + Windward, family = "binomial", id=Line_Id,corstr="independence", data=dat, scale.fix=TRUE)
beep(sound = 1, expr = NULL)

#Removing one term
comb_A_s1_lat  <- geeglm(Pres~ Escarp + Abyss + Zvelv + Channel + Windward, family = "binomial", id=Line_Id,corstr="independence", data=dat, scale.fix=TRUE)
beep(sound = 1, expr = NULL)

comb_A_s1_escarp  <- geeglm(Pres~ Lat + Abyss + Zvelv + Channel + Windward, family = "binomial", id=Line_Id,corstr="independence", data=dat, scale.fix=TRUE)
beep(sound = 1, expr = NULL)

comb_A_s1_abyss  <- geeglm(Pres~ Lat + Escarp + Zvelv + Channel + Windward, family = "binomial", id=Line_Id,corstr="independence", data=dat, scale.fix=TRUE)
beep(sound = 1, expr = NULL)

comb_A_s1_zvelv  <- geeglm(Pres~ Lat + Escarp + Abyss + Channel + Windward, family = "binomial", id=Line_Id,corstr="independence", data=dat, scale.fix=TRUE)
beep(sound = 1, expr = NULL)

comb_A_s1_channel  <- geeglm(Pres~ Lat + Escarp + Abyss + Zvelv + Windward, family = "binomial", id=Line_Id,corstr="independence", data=dat, scale.fix=TRUE)
beep(sound = 1, expr = NULL)

comb_A_s1_windward  <- geeglm(Pres~ Lat + Escarp + Abyss + Zvelv + Channel, family = "binomial", id=Line_Id,corstr="independence", data=dat, scale.fix=TRUE)
beep(sound = 1, expr = NULL)

QIC(comb_A_s1_null) #30822.41 LOWEST (stop)
QIC(comb_A_s1_lat) #30995.69
QIC(comb_A_s1_escarp) #32088.88
QIC(comb_A_s1_abyss) #30844.21
QIC(comb_A_s1_zvelv) #31184.73
QIC(comb_A_s1_channel) #31351.94
QIC(comb_A_s1_windward) #34449.72

#Keep removing variables in turn until the model with the lowest QIC is the "null" model

#Rearrange terms 
#based on how much removing them increases QIC
comb_A  <- geeglm(Pres~ Windward + Escarp + Channel + Zvelv + Lat + Abyss, family = "binomial", id=Line_Id,corstr="independence", data=dat, scale.fix=TRUE)

#Do for all potential models. 


#----------------STEP 4: Model validation------------------------------

##Stepwise cross-validation
# From Eguiguren et al. 2019

#Best model (that you are testing) --> remove variables from this model in turn to compare
comb_1 <- geeglm(Pres~ , family = "binomial", id=Line_Id,corstr="independence", data=dat, scale.fix=TRUE)
beep(sound = 1, expr = NULL)

model_candidate <- comb_1 

#Remove one encounter at a time
dat$Line_Id<-factor(dat$Line_Id) 
encounters <- levels(dat$Line_Id) 

#Remove from data
datFull$Line_Id <-as.character(datFull$Line_Id)

remove_index<-vector("list", length(encounters))
for(i in seq(1:length(encounters))){
  
  remove_index[[i]]<- which(datFull$Line_Id==encounters[[i]])
}

#Make lists for storing
sub_data1 <-vector("list",length(encounters))
sub_data_small <-vector("list",length(encounters))

accuracy_vals     <-vector("list",length(encounters))
qic_vals<-vector("list",length(encounters))
auc_vals<-vector("list",length(encounters))

lenc<- vector("list",length(encounters))
pres_abs <- vector("list",length(encounters)) 


#LOOP STARTS
for(i in seq(1:length(encounters))){
  
  
  #subset data
  
  
  sub_data1[[i]]      <-datFull[-c(remove_index[[i]]),]
  sub_data_small[[i]]  <-datFull[c(remove_index[[i]]),]
  lenc[[i]] <-dim(sub_data_small[[i]])[1]
  pres_abs[[i]] <-sub_data_small[[i]]$pres[1]
  
  sub_data1[[i]]$Line_Id <- as.factor(sub_data1[[i]]$Line_Id)
  sub_data_small[[i]]$Line_Id <- as.factor(sub_data_small[[i]]$Line_Id)
  
  #center and scale variables
  
  ##Variables for big dataset
  sub_data1[[i]]$Lat             <- as.vector(scale(sub_data1[[i]]$Lat))
  
  sub_data1[[i]]$Long       <- as.vector(scale(sub_data1[[i]]$Long))
  
  sub_data1[[i]]$Depth   <- as.vector(scale(sub_data1[[i]]$Depth))
  
  sub_data1[[i]]$logSlope          <- as.vector(scale(sub_data1[[i]]$logSlope))
  
  sub_data1[[i]]$Canyon               <- as.vector(scale(sub_data1[[i]]$Canyon))
  
  sub_data1[[i]]$Escarp               <- as.vector(scale(sub_data1[[i]]$Escarp))
  
  sub_data1[[i]]$Abyss               <- as.vector(scale(sub_data1[[i]]$Abyss))
  
  sub_data1[[i]]$Shelf               <- as.vector(scale(sub_data1[[i]]$Shelf))
  
  sub_data1[[i]]$Ecurr               <- as.vector(scale(sub_data1[[i]]$Ecurr))
  
  sub_data1[[i]]$Ncurr               <- as.vector(scale(sub_data1[[i]]$Ncurr))
  
  sub_data1[[i]]$Zvelv               <- as.vector(scale(sub_data1[[i]]$Zvelv))
  
  sub_data1[[i]]$Mvelv               <- as.vector(scale(sub_data1[[i]]$Mvelv))
  
  sub_data1[[i]]$Inflow               <- as.vector(scale(sub_data1[[i]]$Inflow))
  
  sub_data1[[i]]$Channel               <- as.vector(scale(sub_data1[[i]]$Channel))
  
  sub_data1[[i]]$Chla               <- as.vector(scale(sub_data1[[i]]$Chla))
  
  sub_data1[[i]]$Island               <- (sub_data1[[i]]$Island)
  
  sub_data1[[i]]$Windward               <- (sub_data1[[i]]$Windward)
  
  #sub data small
  sub_data_small[[i]]$Lat             <- as.vector(scale(sub_data_small[[i]]$Lat))
  
  sub_data_small[[i]]$Long       <- as.vector(scale(sub_data_small[[i]]$Long))
  
  sub_data_small[[i]]$Depth   <- as.vector(scale(sub_data_small[[i]]$Depth))
  
  sub_data_small[[i]]$logSlope          <- as.vector(scale(sub_data_small[[i]]$logSlope))
  
  sub_data_small[[i]]$Canyon               <- as.vector(scale(sub_data_small[[i]]$Canyon))
  
  sub_data_small[[i]]$Escarp               <- as.vector(scale(sub_data_small[[i]]$Escarp))
  
  sub_data_small[[i]]$Abyss               <- as.vector(scale(sub_data_small[[i]]$Abyss))
  
  sub_data_small[[i]]$Shelf               <- as.vector(scale(sub_data_small[[i]]$Shelf))
  
  sub_data_small[[i]]$Ecurr               <- as.vector(scale(sub_data_small[[i]]$Ecurr))
  
  sub_data_small[[i]]$Ncurr               <- as.vector(scale(sub_data_small[[i]]$Ncurr))
  
  sub_data_small[[i]]$Zvelv               <- as.vector(scale(sub_data_small[[i]]$Zvelv))
  
  sub_data_small[[i]]$Mvelv               <- as.vector(scale(sub_data_small[[i]]$Mvelv))
  
  sub_data_small[[i]]$Inflow               <- as.vector(scale(sub_data_small[[i]]$Inflow))
  
  sub_data_small[[i]]$Channel               <- as.vector(scale(sub_data_small[[i]]$Channel))
  
  sub_data_small[[i]]$Chla               <- as.vector(scale(sub_data_small[[i]]$Chla))
  
  sub_data_small[[i]]$Island               <- (sub_data_small[[i]]$Island)
  
  sub_data_small[[i]]$Windward               <- (sub_data_small[[i]]$Windward)
  
  
  #model fitting----
  
  ####model----
  formula_null<- as.formula(unlist(model_candidate$formula))#Formula stuff, fits the model inputted before
  
  
  #fit model with big subset 
  fit_null  <- geeglm(as.formula(formula_null), family = "binomial", id=Line_Id,corstr="independence", data=sub_data1[[i]], scale.fix=TRUE)
  #predict response for small subset
  pr_null   <-predict(fit_null, sub_data1[[i]],type="response")##
  pred_null <- prediction(pr_null,sub_data1[[i]]$Pres)
  perf_null <- performance(pred_null, measure = "tpr", x.measure = "fpr")
  
  auc_vals[[i]] <- performance(pred_null, measure ="auc")@y.values
  
  
  #choose cutoff
  
  y_null  <-as.data.frame(perf_null@y.values)
  x_null  <-as.data.frame(perf_null@x.values)
  fi_null <- atan(y_null/x_null) - pi/4                                         # to calculate the angle between the 45? line and the line joining the origin with the point (x;y) on the ROC curve
  L_null  <- sqrt(x_null^2+y_null^2)                                            # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
  d_null  <- L_null*sin(fi_null) 
  names(d_null)<-"vals"
  d_null$vals[which(d_null$vals=="NaN")]<-NA
  maxd_null<-max(d_null$vals, na.rm = T)
  maxpos_null<-which(d_null$vals==maxd_null)
  
  
  alpha_null <-as.data.frame(perf_null@alpha.values)
  thresh_null<-alpha_null[maxpos_null,] #AUC cutoff
  
  #build confusion matrix 
  datlength<-dim(sub_data_small[[i]])
  DATA_null           <- matrix(0,datlength[1],3)
  DATA_null           <- as.data.frame(DATA_null)
  names(DATA_null)    <-c("plotID", "Observed", "Predicted")
  DATA_null$plotID    <-1:datlength[1]
  DATA_null$Observed  <-sub_data_small[[i]]$Pres
  DATA_null$Predicted <-predict(fit_null, sub_data_small[[i]], type="response")
  conmat_null<-cmx(DATA_null, threshold = thresh_null)#this is the cutoff, prints confusion matrix
  
  sub_data_small[[i]]$predicted <-predict(fit_null, sub_data_small[[i]], type="response")#
  sub_data_small[[i]]$predicted_bin<-NA
  sub_data_small[[i]]$predicted_bin[which(sub_data_small[[i]]$predicted>thresh_null)]<-1
  sub_data_small[[i]]$predicted_bin[which(sub_data_small[[i]]$predicted<=thresh_null)]<-0
  
  accuracy_vals[[i]]<-length(which(sub_data_small[[i]]$predicted_bin==sub_data_small[[i]]$Pres))/dim(sub_data_small[[i]])[1]
  
  
  
  
  
}

#Stepwise cross-validation Output
formula_null
acc_vals<-unlist(accuracy_vals) #Accuracy_vals tells you the number of datapoints assigned correctly 
mean(acc_vals) #That is mean percentage right for each time remove an encounter
std.error(acc_vals)
#Then you can compare between models of different variable combinations and choose the one with highest acc_vals

##Calculate AUC and predictive accuracy (from Pirotta et al. 2011)

#Model
comb_1  <- geeglm(Pres~ , family = "binomial", id=Line_Id,corstr="independence", data=dat, scale.fix=TRUE)
beep(sound = 1, expr = NULL)

#QIC
QIC(comb_1) 

# Construction of the ROC curve
pr <- predict(comb_1,dat, type="response")                       # the final model is used to predict the data on the response scale (i.e. a value between 0 and 1)
pred <- prediction(pr,dat$Pres)                                    # to specify the vector of predictions (pr) and the vector of labels (i.e. the observed values "Pres")
perf <- performance(pred, measure="tpr", x.measure="fpr")          # to assess model performance in the form of the true positive rate and the false positive rate
plot(perf, colorize=TRUE, print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5)) # to plot the ROC curve

# Choice of the best cut-off probability

y<-as.data.frame(perf@y.values)
x<-as.data.frame(perf@x.values)
fi <- atan(y/x) - pi/4                                             # to calculate the angle between the 45° line and the line joining the origin with the point (x;y) on the ROC curve
L <- sqrt(x^2+y^2)                                                 # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
d <- L*sin(fi)                                                     # to calculate the distance between the 45° line and the ROC curve
write.table(d,"")                                                  # to write a table with the computed distances

# The table should then be opened in Microsoft Excel to find the maximum distance with the command "Sort", and the relative position (i.e. the number of the corresponding record)
#Use data>Text to column in Excel if all pasted in the same column
# MAX d= 0.229386843755222 --> position 16007

alpha<-as.data.frame(perf@alpha.values)                            # the alpha values represent the corresponding cut-offs
alpha[16007,]                                                       # to identify the alpha value (i.e. the cut-off) that corresponds to the maximum distance between the 45° line and the curve

# Best cutoff:  0.3400672
# This value can now be used to build the confusion matrix:

DATA<-matrix(0,26776,3)                                             # to build a matrix with 3 columns and n rows, where n is the dimension of the data set (here 8319 - the number of rows can be checked with dim(dat)) 
DATA<-as.data.frame(DATA)
names(DATA)<-c("plotID","Observed","Predicted")
DATA$plotID<-1:26776                                                # the first column is filled with an ID value that is unique for each row
DATA$Observed<-dat$Pres                                            # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(comb_1,dat,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = 0.3400672)                                   # the identified cut-off must be used here

#AUC
auc <- performance(pred, measure="auc")
auc@y.values

#Predictive accuracy
#all correct / all --> using data from confusion matrix
(3626+15889)/26776


#----------------STEP 5: Prediction maps-------------------------------

#From Pirotta et al. 2011

finalMOD  <- geeglm(Pres~ , family = "binomial", id=Line_Id,corstr="independence", data=dat, scale.fix=TRUE)
beep(sound = 1, expr = NULL)

pred<-read.table("", header=TRUE)            # import a grid of points created in a GIS. Each point should be associated with the values of the covariates that are retained in the final model
pred<-cbind(pred,predict(finalMOD,pred, type="response"))              # predict animal occurrence using the final model
write.table(pred, "") 

#Import this new dataset into arcGIS for visualization


#----------------OTHER---------------------------------------------------------

##ACF (make sure encounter is a good blocking variable)
#Need to do at the start, and again with the final models to make sure they are not auto correlated

gee.original <-geeglm(Pres~ , family = "binomial", id=Line_Id,corstr="independence", data=dat, scale.fix=TRUE)

ACF.original <-runACF(dat$Line_Id, gee.original, store=F)

##Effect plots
finalFULL  <- geeglm(Pres~ , family = "binomial", id=Line_Id,corstr="independence", data=dat, scale.fix=TRUE)

plotFULL  <- plot(allEffects(finalFULL), main="", ylab="")

