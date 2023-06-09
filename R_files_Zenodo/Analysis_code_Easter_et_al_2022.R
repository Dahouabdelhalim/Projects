#Set working directory to where data files are located on your computer
setwd()

#Load packages
library(nlme)
library(MuMIn)
library(ggplot2)
library(emmeans)
library(ggbeeswarm)

#########################################################
#Effects of behavioral variation on triangle transitivity
#Interaction directionality: From active
#########################################################

#Import data
Ttri_Data<-read.csv("Default values_Interactions from actives_Triangle transitivity.csv", header=TRUE)

Ttri_Data$ConditionF<-factor(Ttri_Data$Condition, levels=c("Uniform","ActVar","TurnVar","Uncorr","Corr"))
Ttri_Data$DensityStd<-(Ttri_Data$Density-mean(Ttri_Data$Density))/sd(Ttri_Data$Density)

#Run global model
M1<-lm(Ttri ~ ConditionF * DensityStd, data=Ttri_Data)

#Check residual plots
#Severe non-normality present
#Heteroscedasticity apparent for both behavioral condition and network density
R1<-resid(M1)
F1<-fitted(M1)
qqnorm(R1)
qqline(R1)
plot(R1~Ttri_Data$ConditionF)
plot(R1~Ttri_Data$DensityStd)
plot(R1~F1)

#Fit a number of variance structure to attempt to model this heteroscedasticity
M1<-gls(Ttri ~ ConditionF * DensityStd, data=Ttri_Data)
M2<-gls(Ttri ~ ConditionF * DensityStd, weights=varIdent(form=~1|ConditionF), data=Ttri_Data)
M3<-gls(Ttri ~ ConditionF * DensityStd, weights=varFixed(~DensityStd), data=Ttri_Data)
M4<-gls(Ttri ~ ConditionF * DensityStd, weights=varPower(form=~DensityStd), data=Ttri_Data)
M5<-gls(Ttri ~ ConditionF * DensityStd, weights=varExp(form=~DensityStd), data=Ttri_Data)
M6<-gls(Ttri ~ ConditionF * DensityStd, weights=varConstPower(form=~DensityStd), data=Ttri_Data)
M7<-gls(Ttri ~ ConditionF * DensityStd, weights=varPower(form=~DensityStd|ConditionF), data=Ttri_Data)
M8<-gls(Ttri ~ ConditionF * DensityStd, weights=varExp(form=~DensityStd|ConditionF), data=Ttri_Data)
M9<-gls(Ttri ~ ConditionF * DensityStd, weights=varConstPower(form=~DensityStd|ConditionF), data=Ttri_Data)
M10<-gls(Ttri ~ ConditionF * DensityStd, weights=varComb(varFixed(~DensityStd), varIdent(form=~1|ConditionF)), data=Ttri_Data)
M11<-gls(Ttri ~ ConditionF * DensityStd, weights=varComb(varPower(form=~DensityStd), varIdent(form=~1|ConditionF)), data=Ttri_Data)
M12<-gls(Ttri ~ ConditionF * DensityStd, weights=varComb(varExp(form=~DensityStd), varIdent(form=~1|ConditionF)), data=Ttri_Data)
M13<-gls(Ttri ~ ConditionF * DensityStd, weights=varComb(varConstPower(form=~DensityStd), varIdent(form=~1|ConditionF)), data=Ttri_Data)

#Select best-supported variant based on AICc
#In this case, M8 received the greatest support
AICc(M8)

#Inspection of residuals suggests adequate model specification
R1<-resid(M8, type="normalized")
F1<-fitted(M8)
qqnorm(R1)
qqline(R1)
plot(R1~Ttri_Data$ConditionF)
plot(R1~Ttri_Data$DensityStd)
plot(R1~F1)

M8<-gls(Ttri ~ ConditionF * DensityStd, weights=varExp(form=~DensityStd|ConditionF), data=Ttri_Data)


#Create the candidate model set for multimodel inference from all models nested within the global model
#As models with different fixed effects structures are being compared, ML estimation is used
dM8<-dredge(M8, REML=FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dM8.Reduced<-subset(dM8, !nested(.))

#Obtain model-averaged estimates from 95% confidence set
M8modAvg<-model.avg(dM8.Reduced, fit = TRUE, subset = cumsum(weight) <= 1)

#Get MAEs, USEs, and 95% CIs reported in Table 2
cbind(coefTable(M8modAvg, revised.var=TRUE, full=TRUE), confint(M8modAvg, revised.var=TRUE, full=TRUE))[,-3]

#######################
#To reproduce Figure 4a
#######################

predgls<-predict(M8modAvg)

ggplot(Ttri_Data, aes(x = Density, y = Ttri, color = Condition))+
	ylab("Triangle transitivity") + xlab("Network density") + 
	geom_point(show.legend=FALSE, position = position_jitter(w = 75, h = 0)) + 
	geom_line(aes(y = predgls), size = 1.5)+
	scale_color_discrete(name="Condition", labels = c("Activity: Variable", "Correlated", "TI: Variable", "Uncorrelated", "Uniform")) +
	scale_x_continuous(limits=c(70,1580),breaks=c(150,450,750,1050,1350)) + 
	scale_y_continuous(limits=c(-1,1), breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)) + 
	theme(legend.position = c(0.7, 0.2), legend.box="horizontal") +
	theme(axis.text=element_text(size=18), axis.title=element_text(size=18), legend.text=element_text(size=18), legend.title=element_text(size=18))+
	theme(axis.title.x=element_text(margin=margin(t=12)), axis.title.y=element_text(margin=margin(r=12)))+
	theme(panel.background=element_rect(fill="white"), legend.key=element_rect(fill="white"), axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1))

#########################################################
#Effects of behavioral variation on triangle transitivity
#Interaction directionality: To active
#########################################################

#Import data
Ttri_Data<-read.csv("Default values_Interactions to actives_Triangle transitivity.csv", header=TRUE)

Ttri_Data$ConditionF<-factor(Ttri_Data$Condition, levels=c("Uniform","ActVar","TurnVar","Uncorr","Corr"))
Ttri_Data$DensityStd<-(Ttri_Data$Density-mean(Ttri_Data$Density))/sd(Ttri_Data$Density)

#Run global model
M1<-lm(Ttri ~ ConditionF * DensityStd, data=Ttri_Data)

#Check residual plots
#Severe non-normality present
#Heteroscedasticity apparent for both behavioral condition and network density
R1<-resid(M1)
F1<-fitted(M1)
qqnorm(R1)
qqline(R1)
plot(R1~Ttri_Data$ConditionF)
plot(R1~Ttri_Data$DensityStd)
plot(R1~F1)

#Fit a number of variance structure to attempt to model this heteroscedasticity
M1<-gls(Ttri ~ ConditionF * DensityStd, data=Ttri_Data)
M2<-gls(Ttri ~ ConditionF * DensityStd, weights=varIdent(form=~1|ConditionF), data=Ttri_Data)
M3<-gls(Ttri ~ ConditionF * DensityStd, weights=varFixed(~DensityStd), data=Ttri_Data)
M4<-gls(Ttri ~ ConditionF * DensityStd, weights=varPower(form=~DensityStd), data=Ttri_Data)
M5<-gls(Ttri ~ ConditionF * DensityStd, weights=varExp(form=~DensityStd), data=Ttri_Data)
M6<-gls(Ttri ~ ConditionF * DensityStd, weights=varConstPower(form=~DensityStd), data=Ttri_Data)
M7<-gls(Ttri ~ ConditionF * DensityStd, weights=varPower(form=~DensityStd|ConditionF), data=Ttri_Data)
M8<-gls(Ttri ~ ConditionF * DensityStd, weights=varExp(form=~DensityStd|ConditionF), data=Ttri_Data)
M9<-gls(Ttri ~ ConditionF * DensityStd, weights=varConstPower(form=~DensityStd|ConditionF), data=Ttri_Data)
M10<-gls(Ttri ~ ConditionF * DensityStd, weights=varComb(varFixed(~DensityStd), varIdent(form=~1|ConditionF)), data=Ttri_Data)
M11<-gls(Ttri ~ ConditionF * DensityStd, weights=varComb(varPower(form=~DensityStd), varIdent(form=~1|ConditionF)), data=Ttri_Data)
M12<-gls(Ttri ~ ConditionF * DensityStd, weights=varComb(varExp(form=~DensityStd), varIdent(form=~1|ConditionF)), data=Ttri_Data)
M13<-gls(Ttri ~ ConditionF * DensityStd, weights=varComb(varConstPower(form=~DensityStd), varIdent(form=~1|ConditionF)), data=Ttri_Data)

#Select best-supported variant based on AICc
#In this case, M8 received the greatest support
AICc(M8)

#Inspection of residuals indicates no remaining issues with model validation
R1<-resid(M8, type="normalized")
F1<-fitted(M8)
qqnorm(R1)
qqline(R1)
plot(R1~Ttri_Data$ConditionF)
plot(R1~Ttri_Data$DensityStd)
plot(R1~F1)

M8<-gls(Ttri ~ ConditionF * DensityStd, weights=varExp(form=~DensityStd|ConditionF), data=Ttri_Data)

#Create the candidate model set for multimodel inference from all models nested within the global model
#As models with different fixed effects structures are being compared, ML estimation is used
dM8<-dredge(M8, REML=FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dM8.Reduced<-subset(dM8, !nested(.))

#Obtain model-averaged estimates from 95% confidence set
M8modAvg<-model.avg(dM8.Reduced, fit = TRUE, subset = cumsum(weight) <= 0.97)

#Get MAEs, USEs, and 95% CIs
cbind(coefTable(M8modAvg, revised.var=TRUE, full=TRUE), confint(M8modAvg, revised.var=TRUE, full=TRUE))[,-3]

#######################
#To reproduce Figure S1
#######################

predgls<-predict(M8modAvg)

ggplot(Ttri_Data, aes(x = Density, y = Ttri, color = Condition))+
	ylab("Triangle transitivity") + xlab("Network density") + 
	geom_point(show.legend=FALSE, position = position_jitter(w = 75, h = 0)) + 
	geom_line(aes(y = predgls), size = 1.5)+
	scale_color_discrete(name="Condition", labels = c("Activity: Variable", "Correlated", "TI: Variable", "Uncorrelated", "Uniform")) +
	scale_x_continuous(limits=c(70,1580),breaks=c(150,450,750,1050,1350)) + 
	scale_y_continuous(limits=c(-1,1), breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)) + 
	theme(legend.position = c(0.7, 0.2), legend.box="horizontal") +
	theme(axis.text=element_text(size=18), axis.title=element_text(size=18), legend.text=element_text(size=18), legend.title=element_text(size=18))+
	theme(axis.title.x=element_text(margin=margin(t=12)), axis.title.y=element_text(margin=margin(r=12)))+
	theme(panel.background=element_rect(fill="white"), legend.key=element_rect(fill="white"), axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1))

###########################################################
#Activity effects on triangle transitivity
#Interaction directionality (where applicable): From active
###########################################################

#Load data
AE_Data<-read.csv("Activity effects_Interactions from actives_Triangle transitivity.csv", header=TRUE)

AE_Data$ConditionF<-factor(AE_Data$Condition, levels=c("Uniform","ActVar"))

#Run global model
A1<-lm(Ttri~ConditionF*Initiation*Direction, data=AE_Data)

#Residuals suggest potential heteroscedasticity for Condition and Direction
R1<-resid(A1)
F1<-fitted(A1)
qqnorm(R1)
qqline(R1)
plot(R1~as.factor(AE_Data$ConditionF))
plot(R1~as.factor(AE_Data$Initiation))
plot(R1~as.factor(AE_Data$Direction))
plot(R1~F1)

#Incorporate variance structures to attempt to model this heteroscedasticity
A1<-gls(Ttri~ConditionF*Initiation*Direction, data=AE_Data)
A2<-gls(Ttri~ConditionF*Initiation*Direction, weights=varIdent(form=~1|ConditionF), data=AE_Data)
A3<-gls(Ttri~ConditionF*Initiation*Direction, weights=varIdent(form=~1|Initiation), data=AE_Data)
A4<-gls(Ttri~ConditionF*Initiation*Direction, weights=varIdent(form=~1|Direction), data=AE_Data)
A5<-gls(Ttri~ConditionF*Initiation*Direction, weights=varComb(varIdent(form=~1|ConditionF),varIdent(form=~1|Initiation)), data=AE_Data)
A6<-gls(Ttri~ConditionF*Initiation*Direction, weights=varComb(varIdent(form=~1|ConditionF),varIdent(form=~1|Direction)), data=AE_Data)
A7<-gls(Ttri~ConditionF*Initiation*Direction, weights=varComb(varIdent(form=~1|Initiation),varIdent(form=~1|Direction)), data=AE_Data)
A8<-gls(Ttri~ConditionF*Initiation*Direction, weights=varComb(varIdent(form=~1|ConditionF),varIdent(form=~1|Initiation),varIdent(form=~1|Direction)), data=AE_Data)

#On the basis of AICc, A6 is preferred

#Residuals appear reasonable
R6<-resid(A6, type="normalized")
F6<-fitted(A6)
qqnorm(R6)
qqline(R6)
plot(R6~as.factor(AE_Data$ConditionF))
plot(R6~as.factor(AE_Data$Initiation))
plot(R6~as.factor(AE_Data$Direction))
plot(R6~F6)


#Create the candidate model set for multimodel inference from all models nested within the global model
#As models with different fixed effects structures are being compared, ML estimation is used
dA6<-dredge(A6, REML=FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dA6.Reduced<-subset(dA6, !nested(.))

#Best-supported model received nearly all the Akaike weights (0.999)
#Includes the three-way interaction between behavioral condition, initiation variant, and directionality variant
#Get estimates reported in Table 3
summary(A6)
confint(A6)

######################
#To reproduce Figure 5
######################

ggplot(AE_Data, aes(x = Condition, y = Ttri, fill=Direction, color=Initiation))+
	ylab("Triangle transitivity") + xlab("Condition") + 
	geom_boxplot(alpha=1, lwd=1.3)+
	scale_fill_manual(values=c("#E69F00","#56B4E9"),name="Directedness",labels=c("From active","Random"))+
	scale_x_discrete(labels=c("Activity: Variable","Uniform"))+
	theme(legend.position = c(0.75, 0.75), legend.text=element_text(size=18), legend.title=element_text(size=18)) +
	theme(axis.text=element_text(size=18), axis.title=element_text(size=18))+
	scale_color_manual(values=c("black", "firebrick"), name="Initiation", labels=c("Active agents","Random"))+
	theme(axis.title.x=element_text(margin=margin(t=12)), axis.title.y=element_text(margin=margin(r=12)))+
	theme(panel.background=element_rect(fill="white"), legend.key=element_rect(fill="white"), axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1))

########################################
#T50
#Interaction directionality: From active
########################################

#Load data
Diff_Data<-read.csv("Diffusion data_Interactions from actives.csv", header=TRUE)

Diff_Data$ConditionF<-factor(Diff_Data$Condition, levels=c("Uniform","ActVar","TurnVar","Uncorrelated","Correlated"))

#Fit model
D1<-gls(T50~ConditionF, data=Diff_Data)

#Residuals indicate non-normality and heteroscedasticity
R1<-resid(D1)
F1<-fitted(D1)
qqnorm(R1)
qqline(R1)
plot(R1~as.factor(Diff_Data$ConditionF))
plot(R1~F1)

D1<-gls(T50~ConditionF, data=Diff_Data)
D2<-gls(T50~ConditionF, weights=varIdent(form=~1|ConditionF), data=Diff_Data)

#Heteroscedasticity is improved, but non-normality persists
R2<-resid(D2, type="normalized")
F2<-fitted(D2)
qqnorm(R2)
qqline(R2)
plot(R2~as.factor(Diff_Data$ConditionF))
plot(R2~F2)

#Fit data with a log transformation
D3<-gls(log(T50)~ConditionF, data=Diff_Data)
D4<-gls(log(T50)~ConditionF, weights=varIdent(form=~1|ConditionF), data=Diff_Data)

#Residuals are reasonable
R3<-resid(D3)
F3<-fitted(D3)
qqnorm(R3)
qqline(R3)
plot(R3~as.factor(Diff_Data$ConditionF))
plot(R3~F3)

#Create the candidate model set for multimodel inference from all models nested within the global model
#As models with different fixed effects structures are being compared, ML estimation is used
dD3<-dredge(D3, REML=FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dD3.Reduced<-subset(dD3, !nested(.))

#Best-supported model received essentially all the Akaike weights (1)
#Includes the effect of behavioral condition
#Get results reported in Table 4
summary(D3)
confint(D3)

#######################
#To reproduce Figure 6a
#######################

ggplot(data = Diff_Data, aes(y=T50, x=Condition, fill=Condition)) + 
	ylab("T50") + xlab("Condition") + 
	geom_boxplot(show.legend=FALSE, lwd=1.25)+
	geom_beeswarm(dodge.width = 0.75, aes(x=Condition, colour=Condition), show.legend=FALSE)+
	scale_x_discrete(labels = c("Act: Var", "Corr", "TI: Var", "Uncorr", "Uniform")) +
	scale_fill_manual(values=c("pink1", "khaki1", "darkseagreen1", "lightcyan1", "plum1"),name="Condition", labels = c("Activity: Variable", "Correlated", "TI: Variable", "Uncorrelated", "Uniform")) +
	scale_y_continuous(limits=c(0,1800)) + 
	theme(legend.position = c(0.85, 0.85)) + 
	theme(axis.text=element_text(size=18),axis.title=element_text(size=20))+
	theme(panel.background=element_rect(fill="white"), legend.key=element_rect(fill="white"), axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1))

########################################
#Transmission efficiency
#Interaction directionality: From active
########################################

#Load data
Diff_Data<-read.csv("Diffusion data_Interactions from actives.csv", header=TRUE)

Diff_Data$ConditionF<-factor(Diff_Data$Condition, levels=c("Uniform","ActVar","TurnVar","Uncorrelated","Correlated"))

#Fit model
E1<-gls(Efficiency~ConditionF, data=Diff_Data)

#Residuals indicate non-normality and heteroscedasticity
R1<-resid(E1)
F1<-fitted(E1)
qqnorm(R1)
qqline(R1)
plot(R1~as.factor(Diff_Data$Condition))
plot(R1~F1)

E1<-gls(Efficiency~ConditionF, data=Diff_Data)
E2<-gls(Efficiency~ConditionF, weights=varIdent(form=~1|ConditionF), data=Diff_Data)

#Heteroscedasticity is improved, but non-normality persists
R2<-resid(E2, type="normalized")
F2<-fitted(E2)
qqnorm(R2)
qqline(R2)
plot(R2~as.factor(Diff_Data$Condition))
plot(R2~F2)

#Fit data with a log transformation
E3<-gls(log(Efficiency)~ConditionF, data=Diff_Data)
E4<-gls(log(Efficiency)~ConditionF, weights=varIdent(form=~1|ConditionF), data=Diff_Data)

#Residuals are reasonable
R3<-resid(E3)
F3<-fitted(E3)
qqnorm(R3)
qqline(R3)
plot(R3~as.factor(Diff_Data$Condition))
plot(R3~F3)

#Create the candidate model set for multimodel inference from all models nested within the global model
#As models with different fixed effects structures are being compared, ML estimation is used
dE3<-dredge(E3, REML=FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dE3.Reduced<-subset(dE3, !nested(.))

#Best-supported model received essentially all the Akaike weights (1)
#Includes the effect of behavioral condition
#Get results reported in Table 5
summary(E3)
confint(E3)

#######################
#To reproduce Figure 4b
#######################

ggplot(data = Diff_Data, aes(y=Efficiency, x=Condition, fill=Condition)) + 
	ylab("Cumulative outgoing edges at T50") + xlab("Condition") + 
	geom_boxplot(show.legend=FALSE, lwd=1.25)+
	geom_beeswarm(dodge.width = 0.75, aes(x=Condition, colour=Condition), show.legend=FALSE)+
	scale_x_discrete(labels = c("Act: Var", "Corr", "TI: Var", "Uncorr", "Uniform")) +
	scale_fill_manual(values=c("pink1", "khaki1", "darkseagreen1", "lightcyan1", "plum1"),name="Condition", labels = c("Activity: Variable", "Correlated", "TI: Variable", "Uncorrelated", "Uniform")) +
	scale_y_continuous(limits=c(0,750)) + 
	theme(legend.position = c(0.85, 0.85)) + 
	theme(axis.text=element_text(size=18),axis.title=element_text(size=20))+
	theme(panel.background=element_rect(fill="white"), legend.key=element_rect(fill="white"), axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1))

########################################
#T50
#Interaction directionality: To active
########################################

#Load data
Diff_Data<-read.csv("Diffusion data_Interactions to actives.csv", header=TRUE)

Diff_Data$ConditionF<-factor(Diff_Data$Condition, levels=c("Uniform","ActVar","TurnVar","Uncorrelated","Correlated"))

#Fit model
D1<-lm(T50~ConditionF, data=Diff_Data)

#Residuals indicate non-normality and heteroscedasticity
R1<-resid(D1)
F1<-fitted(D1)
qqnorm(R1)
qqline(R1)
plot(R1~Diff_Data$ConditionF)
plot(R1~F1)

D1<-gls(T50~ConditionF, data=Diff_Data)
D2<-gls(T50~ConditionF, weights=varIdent(form=~1|ConditionF), data=Diff_Data)

#Heteroscedasticity is improved, but non-normality persists
R2<-resid(D2, type="normalized")
F2<-fitted(D2)
qqnorm(R2)
qqline(R2)
plot(R2~Diff_Data$ConditionF)
plot(R2~F2)

#Fit data with a log transformation
D3<-gls(log(T50)~ConditionF, data=Diff_Data)
D4<-gls(log(T50)~ConditionF, weights=varIdent(form=~1|ConditionF), data=Diff_Data)

#Residuals are reasonable
R3<-resid(D3)
F3<-fitted(D3)
qqnorm(R3)
qqline(R3)
plot(R3~Diff_Data$ConditionF)
plot(R3~F3)

#Create the candidate model set for multimodel inference from all models nested within the global model
#As models with different fixed effects structures are being compared, ML estimation is used
dD3<-dredge(D3, REML=FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dD3.Reduced<-subset(dD3, !nested(.))

#Best-supported model received essentially all the Akaike weights (1)
#Includes the effect of behavioral condition
summary(D3)
confint(D3)

#######################
#To reproduce Figure S2a
#######################

ggplot(data = Diff_Data, aes(y=T50, x=Condition, fill=Condition)) + 
	ylab("T50") + xlab("Condition") + 
	geom_boxplot(show.legend=FALSE, lwd=1.25)+
	geom_beeswarm(dodge.width = 0.75, aes(x=Condition, colour=Condition), show.legend=FALSE)+
	scale_x_discrete(labels = c("Act: Var", "Corr", "TI: Var", "Uncorr", "Uniform")) +
	scale_fill_manual(values=c("pink1", "khaki1", "darkseagreen1", "lightcyan1", "plum1"),name="Condition", labels = c("Activity: Variable", "Correlated", "TI: Variable", "Uncorrelated", "Uniform")) +
	scale_y_continuous(limits=c(0,1800)) + 
	theme(legend.position = c(0.85, 0.85)) + 
	theme(axis.text=element_text(size=18),axis.title=element_text(size=20))+
	theme(panel.background=element_rect(fill="white"), legend.key=element_rect(fill="white"), axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1))

######################################
#Transmission efficiency
#Interaction directionality: To active
######################################

#Load data
Diff_Data<-read.csv("Diffusion data_Interactions to actives.csv", header=TRUE)

Diff_Data$ConditionF<-factor(Diff_Data$Condition, levels=c("Uniform","ActVar","TurnVar","Uncorrelated","Correlated"))

#Fit model
E1<-lm(Efficiency~ConditionF, data=Diff_Data)

#Residuals indicate non-normality and heteroscedasticity
R1<-resid(E1)
F1<-fitted(E1)
qqnorm(R1)
qqline(R1)
plot(R1~Diff_Data$ConditionF)
plot(R1~F1)

E1<-gls(Efficiency~ConditionF, data=Diff_Data)
E2<-gls(Efficiency~ConditionF, weights=varIdent(form=~1|ConditionF), data=Diff_Data)

#Heteroscedasticity is improved, but non-normality persists
R2<-resid(E2, type="normalized")
F2<-fitted(E2)
qqnorm(R2)
qqline(R2)
plot(R2~Diff_Data$ConditionF)
plot(R2~F2)

#Fit data with a log transformation
E3<-gls(log(Efficiency)~ConditionF, data=Diff_Data)
E4<-gls(log(Efficiency)~ConditionF, weights=varIdent(form=~1|ConditionF), data=Diff_Data)

#Residuals are reasonable
R3<-resid(E3)
F3<-fitted(E3)
qqnorm(R3)
qqline(R3)
plot(R3~Diff_Data$ConditionF)
plot(R3~F3)

#Create the candidate model set for multimodel inference from all models nested within the global model
#As models with different fixed effects structures are being compared, ML estimation is used
dE3<-dredge(E3, REML=FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dE3.Reduced<-subset(dE3, !nested(.))

#Best-supported model received essentially all the Akaike weights (1)
#Includes the effect of behavioral condition
summary(E3)
confint(E3)

########################
#To reproduce Figure S2b
########################

ggplot(data = Diff_Data, aes(y=Efficiency, x=Condition, fill=Condition)) + 
	ylab("Efficiency") + xlab("Condition") + 
	geom_boxplot(show.legend=FALSE, lwd=1.25)+
	geom_beeswarm(dodge.width = 0.75, aes(x=Condition, colour=Condition), show.legend=FALSE)+
	scale_x_discrete(labels = c("Act: Var", "Corr", "TI: Var", "Uncorr", "Uniform")) +
	scale_fill_manual(values=c("pink1", "khaki1", "darkseagreen1", "lightcyan1", "plum1"),name="Condition", labels = c("Activity: Variable", "Correlated", "TI: Variable", "Uncorrelated", "Uniform")) +
	scale_y_continuous(limits=c(0,750)) + 
	theme(legend.position = c(0.85, 0.85)) + 
	theme(axis.text=element_text(size=18),axis.title=element_text(size=20))+
	theme(panel.background=element_rect(fill="white"), legend.key=element_rect(fill="white"), axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1))

##########################################
#Sensitivity analysis: Mean activity level
#Interaction directionality: From actives
##########################################

#Load data
SA_Act_Data<-read.csv("Sensitivity analysis_Mean activity_From active_Triangle transitivity.csv", header=TRUE)

SA_Act_Data$ConditionF<-factor(SA_Act_Data$Condition, levels=c("Uniform","ActVar","TurnVar","Uncorrelated","Correlated"))

#Run global model
SA1<-lm(Ttri ~ ConditionF * Am, data=SA_Act_Data)

#Check residual plots
#Heteroscedasticity apparent for behavioral condition
R1<-resid(SA1)
F1<-fitted(SA1)
qqnorm(R1)
qqline(R1)
plot(R1~SA_Act_Data$ConditionF)
plot(R1~SA_Act_Data$Am)
plot(R1~F1)

#Fit a variance structure to attempt to model this heteroscedasticity
SA1<-gls(Ttri ~ ConditionF * Am, data=SA_Act_Data)
SA2<-gls(Ttri ~ ConditionF * Am, weights=varIdent(form=~1|ConditionF), data=SA_Act_Data)

#Residuals appear reasonable
R2<-resid(SA2, type="normalized")
F2<-fitted(SA2)
qqnorm(R2)
qqline(R2)
plot(R2~SA_Act_Data$ConditionF)
plot(R2~SA_Act_Data$Am)
plot(R2~F2)

#Create the candidate model set for multimodel inference from all models nested within the global model
#As models with different fixed effects structures are being compared, ML estimation is used
dSA2<-dredge(SA2, REML=FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dSA2.Reduced<-subset(dSA2, !nested(.))

#Global model received essentially all the Akaike weights (1)
#Obtain estimates reported in Table S2
summary(SA2)
confint(SA2)

#Investigate interactions between mean activity and condition
#Estimates reported in Table S3
emtrends(SA2, pairwise~ConditionF, var="Am", mode="df.error")

########################
#To reproduce Figure S5a
########################

predgls<-predict(SA2)

ggplot(SA_Act_Data, aes(x = Am, y = Ttri, color = Condition))+
	ylab("Triangle transitivity") + xlab("Mean activity level") + 
	geom_point(show.legend=FALSE, position = position_jitter(w = 0.125, h = 0)) + 
	geom_line(aes(y = predgls), size = 1.5)+
	scale_color_discrete(name="Condition", labels = c("Activity: Variable", "Correlated", "TI: Variable", "Uncorrelated", "Uniform")) +
	scale_x_continuous(limits=c(0.35,2.65),breaks=c(0.5,1,1.5,2,2.5)) + 
	scale_y_continuous(limits=c(-1,1), breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)) + 
	theme(legend.position = c(0.7, 0.2), legend.box="horizontal") +
	theme(axis.text=element_text(size=18), axis.title=element_text(size=18), legend.text=element_text(size=18), legend.title=element_text(size=18))+
	theme(axis.title.x=element_text(margin=margin(t=12)), axis.title.y=element_text(margin=margin(r=12)))+
	theme(panel.background=element_rect(fill="white"), legend.key=element_rect(fill="white"), axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1))

#########################################
#Sensitivity analysis: Mean turning index
#Interaction directionality: From actives
#########################################

#Load data
SA_TI_Data<-read.csv("Sensitivity analysis_Mean turning index_From active_Triangle transitivity.csv", header=TRUE)

SA_TI_Data$ConditionF<-factor(SA_TI_Data$Condition, levels=c("Uniform","ActVar","TurnVar","Uncorrelated","Correlated"))

#Run global model
ST1<-lm(Ttri ~ ConditionF * Tm, data=SA_TI_Data)

#Check residual plots
#Heteroscedasticity apparent for behavioral condition
R1<-resid(ST1)
F1<-fitted(ST1)
qqnorm(R1)
qqline(R1)
plot(R1~SA_TI_Data$ConditionF)
plot(R1~SA_TI_Data$Tm)
plot(R1~F1)

#Fit a variance structure to attempt to model this heteroscedasticity
ST1<-gls(Ttri ~ ConditionF * Tm, data=SA_TI_Data)
ST2<-gls(Ttri ~ ConditionF * Tm, weights=varIdent(form=~1|ConditionF), data=SA_TI_Data)

#Residuals appear reasonable
R2<-resid(ST2, type="normalized")
F2<-fitted(ST2)
qqnorm(R2)
qqline(R2)
plot(R2~SA_TI_Data$ConditionF)
plot(R2~SA_TI_Data$Tm)
plot(R2~F2)

#Create the candidate model set for multimodel inference from all models nested within the global model
#As models with different fixed effects structures are being compared, ML estimation is used
dST2<-dredge(ST2, REML=FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dST2.Reduced<-subset(dST2, !nested(.))

#Global model received essentially all the Akaike weights (1)
#Obtain estimates reported in Table S4
summary(ST2)
confint(ST2)

#Investigate interactions between mean turning index and condition
#Estimates reported in Table S5
emtrends(ST2, pairwise~ConditionF, var="Tm", mode="df.error")

########################
#To reproduce Figure S5b
########################

predgls<-predict(ST2)

ggplot(SA_TI_Data, aes(x = Tm, y = Ttri, color = Condition))+
	ylab("Triangle transitivity") + xlab("Mean turning index") + 
	geom_point(show.legend=FALSE, position = position_jitter(w = 7.5, h = 0)) + 
	geom_line(aes(y = predgls), size = 1.5)+
	scale_color_discrete(name="Condition", labels = c("Activity: Variable", "Correlated", "TI: Variable", "Uncorrelated", "Uniform")) +
	scale_x_continuous(limits=c(7,163),breaks=c(30,60,90,120,150)) + 
	scale_y_continuous(limits=c(-1,1), breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)) + 
	theme(legend.position = c(0.7, 0.2), legend.box="horizontal") +
	theme(axis.text=element_text(size=18), axis.title=element_text(size=18), legend.text=element_text(size=18), legend.title=element_text(size=18))+
	theme(axis.title.x=element_text(margin=margin(t=12)), axis.title.y=element_text(margin=margin(r=12)))+
	theme(panel.background=element_rect(fill="white"), legend.key=element_rect(fill="white"), axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1))

#########################################
#T50
#Static Network Diffusion Analysis
#Interaction directionality: From actives
#########################################

#Load data
diffData<-read.csv("Diffusion data_Static networks.csv", header=TRUE)

diffData$ConditionF<-factor(diffData$Condition, levels=c("Uniform","ActVar","TurnVar","Uncorrelated","Correlated"))

M1<-gls(T50 ~ ConditionF, data=diffData)

#Non-normality and heteroscedasticity present
R1<-resid(M1)
F1<-fitted(M1)
qqnorm(R1)
qqline(R1)
plot(R1~as.factor(diffData$ConditionF))
plot(R1~F1)

#Fit a variance structure to attempt to model this heteroscedasticity
M2<-gls(T50 ~ ConditionF, weights=varIdent(form=~1|ConditionF), data=diffData)

#Non-normality persists
R1<-resid(M2, type="normalized")
F1<-fitted(M2)
qqnorm(R1)
qqline(R1)
plot(R1~as.factor(diffData$ConditionF))
plot(R1~F1)

#Apply log transformation
M3<-gls(log(T50) ~ ConditionF, data=diffData)
M4<-gls(log(T50) ~ ConditionF, weights=varIdent(form=~1|ConditionF), data=diffData)

#Residuals appear adequate, though some heteroscedasticity persists
R1<-resid(M4, type="normalized")
F1<-fitted(M4)
qqnorm(R1)
qqline(R1)
plot(R1~as.factor(diffData$ConditionF))
plot(R1~F1)

#Create the candidate model set for multimodel inference from all models nested within the global model
dM4<-dredge(M4, REML=FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dM4.Reduced<-subset(dM4, !nested(.))

#Obtain model-averaged estimates from 95% confidence set
M4modAvg<-model.avg(dM4.Reduced, fit = TRUE)

#Get MAEs, USEs, and 95% CIs reported in Table S6
cbind(coefTable(M4modAvg, revised.var=TRUE, full=TRUE), confint(M4modAvg, revised.var=TRUE, full=TRUE))[,-3]


########################
#To reproduce Figure S6a
########################

ggplot(data = diffData, aes(y=T50, x=Condition, fill=Condition)) + 
	ylab("T50") + xlab("Condition") + 
	geom_boxplot(show.legend=FALSE, lwd=1.25)+
	#geom_beeswarm(dodge.width = 0.75, aes(x=Condition, colour=Condition), show.legend=FALSE)+
	scale_x_discrete(labels = c("Act: Var", "Corr", "TI: Var", "Uncorr", "Uniform")) +
	scale_fill_manual(values=c("pink1", "khaki1", "darkseagreen1", "lightcyan1", "plum1"),name="Condition", labels = c("Activity: Variable", "Correlated", "TI: Variable", "Uncorrelated", "Uniform")) +
	scale_y_continuous(limits=c(0,25)) + 
	theme(legend.position = c(0.85, 0.85)) + 
	theme(axis.text=element_text(size=18),axis.title=element_text(size=20))+
	theme(panel.background=element_rect(fill="white"), legend.key=element_rect(fill="white"), axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1))

##########################################
#Transmission Efficiency
#Static Network Diffusion Analysis
#Interaction directionality: From actives
#########################################

#Load data
diffData<-read.csv("Diffusion data_Static networks.csv", header=TRUE)

diffData$ConditionF<-factor(diffData$Condition, levels=c("Uniform","ActVar","TurnVar","Uncorrelated","Correlated"))

M1<-gls(InteractionsPerInformed ~ ConditionF, data=diffData)

#Non-normality and heteroscedasticity present
R1<-resid(M1)
F1<-fitted(M1)
qqnorm(R1)
qqline(R1)
plot(R1~as.factor(diffData$ConditionF))
plot(R1~F1)

#Fit a variance structure to attempt to model this heteroscedasticity
M2<-gls(InteractionsPerInformed ~ ConditionF, weights=varIdent(form=~1|ConditionF), data=diffData)

#Non-normality persists
R1<-resid(M2, type="normalized")
F1<-fitted(M2)
qqnorm(R1)
qqline(R1)
plot(R1~as.factor(diffData$ConditionF))
plot(R1~F1)

#Apply log transformation
M3<-gls(log(InteractionsPerInformed) ~ ConditionF, data=diffData)
M4<-gls(log(InteractionsPerInformed) ~ ConditionF, weights=varIdent(form=~1|ConditionF), data=diffData)

#Residuals appear reasonable
R1<-resid(M4, type="normalized")
F1<-fitted(M4)
qqnorm(R1)
qqline(R1)
plot(R1~as.factor(diffData$ConditionF))
plot(R1~F1)

#Create the candidate model set for multimodel inference from all models nested within the global model
dM4<-dredge(M4, REML=FALSE)

#Models are removed if they are more complex versions of a model with better support (Richards SA, 2008, J Appl Ecol 45, 218-227)
dM4.Reduced<-subset(dM4, !nested(.))

#Global model received essentially all the Akaike weights (1)
#Estimates reported in Table S7
summary(M4)
confint(M4)

########################
#To reproduce Figure S6b
########################

ggplot(data = diffData, aes(y=InteractionsPerInformed, x=Condition, fill=Condition)) + 
	ylab("Efficiency") + xlab("Condition") + 
	geom_boxplot(show.legend=FALSE, lwd=1.25)+
	geom_beeswarm(dodge.width = 0.75, aes(x=Condition, colour=Condition), show.legend=FALSE)+
	scale_x_discrete(labels = c("Act: Var", "Corr", "TI: Var", "Uncorr", "Uniform")) +
	scale_fill_manual(values=c("pink1", "khaki1", "darkseagreen1", "lightcyan1", "plum1"),name="Condition", labels = c("Activity: Variable", "Correlated", "TI: Variable", "Uncorrelated", "Uniform")) +
	scale_y_continuous(limits=c(0,7)) + 
	theme(legend.position = c(0.85, 0.85)) + 
	theme(axis.text=element_text(size=18),axis.title=element_text(size=20))+
	theme(panel.background=element_rect(fill="white"), legend.key=element_rect(fill="white"), axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1))

##############################