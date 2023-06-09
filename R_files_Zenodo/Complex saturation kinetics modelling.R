#--------------------------------------------------------------------------------------
#
#  July 2022
#  Code associated with publication:
#  "Kinetic modulation of bacterial hydrolases by microbial community structure in coastal waters" 
#   Authors: N. Abad, A. Uranga, B. Ayo, J.M. Arrieta,Z. Baña, I. Artolozaga, I. Azúa, J. Iriberri, Santos J. González-Rojí and M. Unanue
#
#
#  Creative Commons Licence: Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)
#--------------------------------------------------------------------------------------

# July 2022
# Script by N. Abad

#################### NOTES before starting the analysis #####################

# 1) It is required to previously install "nlstools" and "AICcmodavg" packages. Run the following
# code line if needed:

install.packages(c("nlstools","AICcmodavg"))

# 2) Input data is a file with n-rows (samples) and 2-columns: column 1 with header "S" corresponds to
# the different concentrations of substrate used in the experiment and column 2 with header "V" corresponds 
# to the hydrolysis rates.

#################################################

## 1_Load the file, omit the rows where there is no data (NA) and plot Hydrolysis rates vs Substrate concentration.

rm(list=ls())
dev.off()

library(nlstools)
library(AICcmodavg)

file<-file.choose()
setwd<-dirname(file)
setwd(setwd)
data<-read.delim(file)
data.with.na<-data
data<-na.omit(data.with.na)

# set the units in brackets to the enzymatic experiment
plot(data,xlab="Substrate concentration (µM)",ylab="Hydrolysis rate (nM·h-1)",
ylim=c(min(data$V),max(data$V)),xlim=c(min(data$S),max(data$S)))

#################################################

## 2_Perform the fit of different models of increasing complexity to the dataset

# MODEL 1: the simplest, is derived from MM kinetics assuming that all 
# of the data collected are in the region of first order kinetics (S<<Km)
# Parameters Km and Vmax cannot be estimated separately under such circumstances; 
# ratio Vmax/Km trated as one parameter can be determined reliably

### The code that is commented allows a simulation to better estimate 
### the starting values of model 1

#formula<-as.formula(V~(S/Tt))
#preview(formula, data=data, start= list(Tt=5)) 

model1 <- nls(V~(S/Tt), data=data,
            start =list(Tt = 2),
            alg = "port", trace = TRUE,
            control=list(maxiter=200, tol=1e-15))

RSS1 <- sum(residuals(model1)^2)
TSS1 <- sum((data$V- mean(data$V))^2)
summary(model1)
print(paste("R-squared:",format(1 - (RSS1/TSS1),digits=2)))
RSS1

abline(0,1/coef(model1)["Tt"], col="grey", lty=2,lwd=2) 

#####################################################

# MODEL 2: the Henri-Michaelis-Menten equation

# set the units in brackets to the enzymatic experiment
plot(data,xlab="Substrate concentration (µM)",ylab="Hydrolysis rate (nM·h-1)",
ylim=c(min(data$V),max(data$V)),
xlim=c(min(data$S),max(data$S)))

model2 <- nls(V~(Vmax*S)/(Km+S), data=data,
            start =list( Km = 10 , Vmax = max(data$V,na.rm=TRUE)),
            lower =list( Km = 1 , Vmax = min(data$V,na.rm=TRUE)),
            upper =list( Km = max(data$S,na.rm=TRUE),Vmax= 1000),
            alg = "port", trace = TRUE,
            control=list(maxiter=200, tol=1e-15))

RSS2 <- sum(residuals(model2)^2)
TSS2 <- sum((data$V- mean(data$V))^2)
summary(model2)
print(paste("R-squared:  ",format(1 - (RSS2/TSS2),digits=2)))
RSS2
            
if(substr(model2[1],1,5)!="Error")
	{
    	Km<-coef(model2)["Km"] 
   		Vmax<-coef(model2)["Vmax"]
		  smin<-min(data$S)   
		  smax<-max(data$S)
		  ds<-1
		  s <- seq(smin,smax,ds) 
		  n <- length(s) 
		  N <- min(data$V) 
		  for (j in 1:n) 
			{ 
			  N[j] <- (Vmax*s[j])/(Km+s[j]) 
      			}
		  lines(s,N, col="red") 
 		}

## Choose the model best representing experimental data on the basis of the corrected Akaike's 
## Information Criteria (AICc). 

## NOTE: is important to maintain the order in which the code below is written because the models 
## are in  hierarchical order.

anova(model1,model2)
AIC(model1,model2) 

AICc(model1)
AICc(model2) 

##########################################

# MODEL 3: assumes a system of two independent enzymes whose 
# kinetics are described by models 1 and 2

plot(data,xlab="Substrate concentration (µM)",ylab="Hydrolysis rate (nM·h-1)",
     ylim=c(min(data$V),max(data$V)),
     xlim=c(min(data$S),max(data$S)))

### The code that is commented allows a simulation to better estimate 
### the starting values of model 3. 

#formula<-as.formula(V~(S<= "Value")*((Vmax*S)/(Km+S)) + ((S> "Value")* (S/Tt)))
#preview(formula, data=data, start= list(Km= 1, Vmax= 7, Tt=3)) 

# NOTE: in the code line V~(S<= "Value")*((Vmax*S)/(Km+S)) + ((S> "Value")* (S/Tt),
# the "Value" item has to be substituted according to the concentration of substrates
# used in the enzymatic assay

model3 <- nls(V~(S<= 25)*((Vmax*S)/(Km+S)) + ((S> 25)* (S/Tt)), 
	data=data, start =list( Km = 1 , Vmax = 0.5, Tt= 50),
	alg = "port", trace = TRUE,
	control=list(maxiter=200, tol=1e-15))

RSS3 <- sum(residuals(model3)^2)
TSS3 <- sum(( data$V- mean(data$V))^2)
summary(model3)
print(paste("R-squared:  ",format(1 - (RSS3/TSS3),digits=2)))
RSS3

if(substr(model3[1],1,5)!="Error")
{
model.fit<-as.list(coef(model3))
Km<-coef(model3)["Km"] 
Vmax<-coef(model3)["Vmax"]
Tt<-coef(model3)["Tt"]
smin <- min(data$S)
smax <- max(data$S)
ds <- 1 
s <- seq(smin,smax,ds) 
n <- length(s) 
N1 <- min(data$V) 
N2<-  min(data$V)

for (j in 1:n) 
{ 
  N1[j] <- (Vmax*s[j])/(Km+s[j])
  N2[j] <- (s[j])/(Tt)
}

lines(s[1:n],N1[1:n], col="blue", lwd=2)
lines(s[1:n],N2[1:n],col="darkblue", lwd=2)

}

## Choose the model best representing experimental data on the basis of the corrected Akaike's 
## Information Criteria (AICc). 

## NOTE: is important to maintain the order in which the code below is written because the models 
## are in  hierarchical order.

## Depending on the result obtained in the previous AICc test, replace "modelX" by "model1" or "model 2"
anova(modelX,model3)
AIC(modelX,model3) 

AICc(modelX)
AICc(model3) 

#############################################################

# MODEL 4: the most complex model assumes a system of two independent 
# enzymes whose kinetics follow the Henri-Michaelis-Menten 
# equation (high- and low-affinity)

plot(data,xlab="Substrate concentration (µM)",ylab="Hydrolysis rate (nM·h-1)",
     	ylim=c(min(data$V),max(data$V)),
     	xlim=c(min(data$S),max(data$S)))

model4 <- nls(V~(VmaxL*S)/(KmL+S)+(VmaxH*S)/(KmH+S),data=data, 
              start =list(KmH=1, VmaxH= 0.4,
                          KmL=max(data$S,na.rm=TRUE),
                          VmaxL=max(data$V,na.rm=TRUE)),
              alg = "port", trace = TRUE,control=list(maxiter=200, tol=1e-10))

RSS4 <- sum(residuals(model4)^2)
TSS4 <- sum(( data$V- mean(data$V))^2)
summary(model4)
print(paste("R-squared:  ",format(1 - (RSS4/TSS4),digits=2)))
RSS4

if(substr(model4[1],1,5)!="Error")
    	{
model.fit<-as.list(coef(model4))
KmH<-coef(model4)["KmH"] 
VmaxH<-coef(model4)["VmaxH"]
KmL<-coef(model4)["KmL"]
VmaxL<-coef(model4)["VmaxL"]
smin <- min(data$S)
smax <- max(data$S)
ds <- 1 
s <- seq(smin,smax,ds) 
n <- length(s) 
N1 <- min(data$V) 
N2<-  min(data$V)
for (j in 1:n) 
{ 
  N1[j] <- (VmaxH*s[j])/(KmH+s[j])
  N2[j] <- (VmaxL*s[j])/(KmL+s[j])
}

lines(s[1:n],N1[1:n], col="darkmagenta", lwd=2)
lines(s[1:n],N2[1:n],col="darkorchid1", lwd=2)
}

## Choose the model best representing experimental data on the basis of the corrected Akaike's 
## Information Criteria (AICc). 

## NOTE: is important to maintain the order in which the code below is written because the models 
## are in  hierarchical order.

## Depending on the result obtained in the previous AICc test (a.k.a. modelX vs model3), replace "modelY" 

anova(modelY,model4)
AIC(modelY,model4) 
AICc(modelY) 
AICc(model4)

#################### THE END ####################