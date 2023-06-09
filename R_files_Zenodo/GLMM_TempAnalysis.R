# Metadata:

#Mosquito_ID = unique identifier per mosquito
#Replicate: 1, 2 = biological replicate
#Batch: Mosquito cohort - we ignore this as this is mostly captured by variation across replicates
#Day: 3, 6, 9, 12, 15, 18, 21 - days post infection mosquitoes were destructively sampled for virus
#Temperature: 16, 20, 24, 28, 32, 34, 36, 38 C
#Infected: 0 = no, 1 = yes (positive bodies = 1)
#Disseminated: 0 = no, 1 = yes (positive head and legs = 1)
#Infectious = 0 = no, 1 = yes (positive saliva = 1)

#We are interested in the probability of infection, dissemination, and infectiousness out of all 
#mosquitoes exposed and the probability of infectiousness out of all of those mosquitoes successfully 
#infected.

##############################
# PRODUCTS

# 2 charts to copy from R Graphics Device

# 4 "best model" tables, one fore each response variable:
#	BEST_Infected
#	BEST_Disseminated
#	BEST_Infectious
#	BEST_Efficiency


# load packages, run the installs, if needed
library(readxl)	#install.packages("readxl")
library(lme4)	#install.packages("lme4")
library(MuMIn)	#install.packages("MuMIn")

antilogit <- function(x) { exp(x) / (1 + exp(x) ) }

############################################
#load data, center and scale continuous variables

setwd("C:/Users/Courtney Murdock/Documents/Dropbox/Research/Research/Papers&Posters&Talks/Tenure/Zika_NSFRAPID/Aim 2/ProcRoySoc/Dryad")
dat <- read_excel("InfectionData.xlsx")
dat$Tfactor <- as.factor(dat$Temperature)
dat$Dfactor <- as.factor(dat$Day)
dat$Tscale <- scale(dat$Temperature)
dat$Dscale <- scale(dat$Day)

summary(dat) # by-column summaries of variables
str(dat)	# structure of the data in each column

#infected only datafile
datinf <- subset(dat,Infected==1)
datinf$Tfactor <- as.factor(datinf$Temperature)
datinf$Dfactor <- as.factor(datinf$Day)
datinf$Tscale <- scale(datinf$Temperature)
datinf$Dscale <- scale(datinf$Day)

xtabs(~Infected+Day+Temperature, data =dat)
xtabs(~Disseminated+Day+Temperature, data =dat)
xtabs(~Infectious+Day+Temperature, data =dat)
xtabs(~Infectious+Day+Temperature, data =datinf)

xtabs(~Infected+Day, data =dat)
xtabs(~Disseminated+Day, data =dat)
xtabs(~Infectious+Day, data =dat)

xtabs(~Infected+Infectious, data =dat)
xtabs(~Infected+Infectious, data =dat)
xtabs(~Infected+Infectious, data =dat)

xtabs(~Infectious, data = subset(dat,Infected==1))
aggregate(Infectious~Temperature, data = subset(dat,Infected==1),mean)
aggregate(Infectious~Day, data = subset(dat,Infected==1),mean)


#######################################
# 3-panel figure: Infected, Disseminated, Infectious
df <- dat
windows(record=T)
par(mfrow = c(1,3))
days <- seq(min(dat$Day),max(dat$Day), length = 20)

## PANEL 1

	response.name <- "Infected"
	response<-as.matrix(df[,response.name])
	m1<-glmer(response ~ Dscale + Tscale + (1|Batch), data = df, family = binomial)
	m2<-glmer(response ~ Dscale*Tscale + (1|Batch), data = df, family = binomial)
	m3<-glmer(response ~ Dscale + I(Dscale^2) + Tscale + (1|Batch),	data = df, family = binomial)
	m4<-glmer(response ~ Tscale + I(Tscale^2) + Dscale + (1|Batch),	data = df, family = binomial)
	m5<-glmer(response ~ Tscale + I(Tscale^2) + Dscale + I(Dscale^2) + (1|Batch), data = df, family = binomial)
	m6<-glmer(response ~ Tscale*Dscale + I(Dscale^2)*I(Tscale^2) + (1|Batch), data = df, family = binomial)
	m7<-glmer(response ~ Tscale*Dscale + I(Dscale^2) + (1|Batch), data = df, family = binomial)
	m8<-glmer(response ~ Tscale*Dscale + I(Tscale^2) + (1|Batch), data = df, family = binomial)
	modlist <- list(m1,m2,m3,m4,m5,m6,m7,m8)
	(selection.table <- model.sel(modlist))
	
	BEST_Infected <- modlist[[as.numeric(rownames(selection.table[1,]))]]

	calc_y <- function(T = T){
		m <- modlist[[as.numeric(rownames(selection.table[1,]))]]
		#MODEL EXPECTATION ACROSS TEMPERATURE
		newdat<-data.frame(Dscale=seq(-1.42757,1.57742,length=20), 
			Tscale=rep(((T-attr(dat$Tscale, "scaled:center"))/attr(dat$Tscale, "scaled:scale")),20))
		mm<-model.matrix(~ Tscale * Dscale + I(Dscale^2) * I(Tscale^2),newdat)	
		yy<-mm%*%fixef(m)
		return(antilogit(yy))
	}
	plot(Infectious~Day,dat,pch="", ylab = response.name)
	lines(days, calc_y(T=16), col = "dark blue", lwd = 2)
	lines(days, calc_y(T=20), col = "light blue", lwd = 2)
	lines(days, calc_y(T=24), col = "dark green", lwd = 2)
	lines(days, calc_y(T=28), col = "light green", lwd = 2)
	lines(days, calc_y(T=32), col = "yellow", lwd = 2)
	lines(days, calc_y(T=34), col = "orange", lwd = 2)
	lines(days, calc_y(T=36), col = "magenta", lwd = 2)
	lines(days, calc_y(T=38), col = "red", lwd = 2)
	
## PANEL 2

	response.name <- "Disseminated"
	response<-as.matrix(df[,response.name])
	m1<-glmer(response ~ Dscale + Tscale + (1|Batch), data = df, family = binomial)
	m2<-glmer(response ~ Dscale*Tscale + (1|Batch), data = df, family = binomial)
	m3<-glmer(response ~ Dscale + I(Dscale^2) + Tscale + (1|Batch),	data = df, family = binomial)
	m4<-glmer(response ~ Tscale + I(Tscale^2) + Dscale + (1|Batch),	data = df, family = binomial)
	m5<-glmer(response ~ Tscale + I(Tscale^2) + Dscale + I(Dscale^2) + (1|Batch), data = df, family = binomial)
	m6<-glmer(response ~ Tscale*Dscale + I(Dscale^2)*I(Tscale^2) + (1|Batch), data = df, family = binomial)
	m7<-glmer(response ~ Tscale*Dscale + I(Dscale^2) + (1|Batch), data = df, family = binomial)
	m8<-glmer(response ~ Tscale*Dscale + I(Tscale^2) + (1|Batch), data = df, family = binomial)
	modlist <- list(m1,m2,m3,m4,m5,m6,m7,m8)
	(selection.table <- model.sel(modlist))
	
	BEST_Disseminated <- modlist[[as.numeric(rownames(selection.table[1,]))]]

	plot(Infectious~Day,dat,pch="", ylab = response.name)
	lines(days, calc_y(T=16), col = "dark blue", lwd = 2)
	lines(days, calc_y(T=20), col = "light blue", lwd = 2)
	lines(days, calc_y(T=24), col = "dark green", lwd = 2)
	lines(days, calc_y(T=28), col = "light green", lwd = 2)
	lines(days, calc_y(T=32), col = "yellow", lwd = 2)
	lines(days, calc_y(T=34), col = "orange", lwd = 2)
	lines(days, calc_y(T=36), col = "magenta", lwd = 2)
	lines(days, calc_y(T=38), col = "red", lwd = 2)

## PANEL 3

	response.name <- "Infectious"
	response<-as.matrix(df[,response.name])
	m1<-glmer(response ~ Dscale + Tscale + (1|Batch), data = df, family = binomial)
	m2<-glmer(response ~ Dscale*Tscale + (1|Batch), data = df, family = binomial)
	m3<-glmer(response ~ Dscale + I(Dscale^2) + Tscale + (1|Batch),	data = df, family = binomial)
	m4<-glmer(response ~ Tscale + I(Tscale^2) + Dscale + (1|Batch),	data = df, family = binomial)
	m5<-glmer(response ~ Tscale + I(Tscale^2) + Dscale + I(Dscale^2) + (1|Batch), data = df, family = binomial)
	m6<-glmer(response ~ Tscale*Dscale + I(Dscale^2)*I(Tscale^2) + (1|Batch), data = df, family = binomial)
	m7<-glmer(response ~ Tscale*Dscale + I(Dscale^2) + (1|Batch), data = df, family = binomial)
	m8<-glmer(response ~ Tscale*Dscale + I(Tscale^2) + (1|Batch), data = df, family = binomial)
	modlist <- list(m1,m2,m3,m4,m5,m6,m7,m8)
	(selection.table <- model.sel(modlist))
	
	BEST_Infectious <- modlist[[as.numeric(rownames(selection.table[1,]))]]

	plot(Infectious~Day,dat,pch="", ylab = response.name)
	lines(days, calc_y(T=16), col = "dark blue", lwd = 2)
	lines(days, calc_y(T=20), col = "light blue", lwd = 2)
	lines(days, calc_y(T=24), col = "dark green", lwd = 2)
	lines(days, calc_y(T=28), col = "light green", lwd = 2)
	lines(days, calc_y(T=32), col = "yellow", lwd = 2)
	lines(days, calc_y(T=34), col = "orange", lwd = 2)
	lines(days, calc_y(T=36), col = "magenta", lwd = 2)
	lines(days, calc_y(T=38), col = "red", lwd = 2)

	
#######################################
# 3-panel figure: Infected, Disseminated, Infectious
df <- datinf
windows(record=T)
par(mfrow = c(1,1))

## PANEL 1
	response.name <- "Infectious"
	response<-as.matrix(df[,response.name])
	m1<-glmer(response ~ Dscale + Tscale + (1|Batch), data = df, family = binomial)
	m2<-glmer(response ~ Dscale*Tscale + (1|Batch), data = df, family = binomial)
	m3<-glmer(response ~ Dscale + I(Dscale^2) + Tscale + (1|Batch),	data = df, family = binomial)
	m4<-glmer(response ~ Tscale + I(Tscale^2) + Dscale + (1|Batch),	data = df, family = binomial)
	m5<-glmer(response ~ Tscale + I(Tscale^2) + Dscale + I(Dscale^2) + (1|Batch), data = df, family = binomial)
	m6<-glmer(response ~ Tscale*Dscale + I(Dscale^2)*I(Tscale^2) + (1|Batch), data = df, family = binomial)
	m7<-glmer(response ~ Tscale*Dscale + I(Dscale^2) + (1|Batch), data = df, family = binomial)
	m8<-glmer(response ~ Tscale*Dscale + I(Tscale^2) + (1|Batch), data = df, family = binomial)
	modlist <- list(m1,m2,m3,m4,m5,m6,m7,m8)
	(selection.table <- model.sel(modlist))
	
	BEST_Efficiency <- modlist[[as.numeric(rownames(selection.table[1,]))]]

	plot(Infectious~Day,dat,pch="", ylab = "Efficiency")
	lines(days, calc_y(T=16), col = "dark blue", lwd = 2)
	lines(days, calc_y(T=20), col = "light blue", lwd = 2)
	lines(days, calc_y(T=24), col = "dark green", lwd = 2)
	lines(days, calc_y(T=28), col = "light green", lwd = 2)
	lines(days, calc_y(T=32), col = "yellow", lwd = 2)
	lines(days, calc_y(T=34), col = "orange", lwd = 2)
	lines(days, calc_y(T=36), col = "magenta", lwd = 2)
	lines(days, calc_y(T=38), col = "red", lwd = 2)
	
	
	
	