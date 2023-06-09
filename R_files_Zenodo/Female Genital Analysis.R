rm(list=ls())

library(ggplot2)
library(rcompanion)
library(scales)
library(nlme)
library(lme4)
library(car)
library(emmeans)
library(sjstats)
library(subselect)
library(HoRM)
library(plyr)
library(heplots)



data <- read.csv("Female Data.csv", header=TRUE, row.names=1)

names(data)

######################
#######  PCA  ########
######################
##Transformed data
pca<-princomp(data[,21:28])#, cor=TRUE)

summary(pca) # print variance accounted for by each axis
pca$loadings #loadings, not sure why it only prints the big ones
#Calculate relative contributions of each variable to a PC
loads <-with(pca, unclass(loadings))
loads
aload <- abs(loads) ## save absolute values
sweep(aload, 2, colSums(aload), "/")
colSums(sweep(aload, 2, colSums(aload), "/"))

#write.csv(loadings(pca),"loadings.csv") # export the pc loadings
ev <- eigen(cor(data[,21:28])) # get eigenvalues
ev$values # print eigenvalues for each axis

scores <- pca$scores
scores <- data.frame(scores)
data2 <- cbind(data, scores)

names(data2)

###############  MANCOVA on Significant PC Axes  ###############
################################################################
my.model <- lm(cbind(Comp.1, Comp.2, Comp.3) 
               ~ logMass + Sulf + Drainage + logMass*Sulf + logMass*Drainage + Sulf*Drainage + logMass*Sulf*Drainage, data=data2)
fit <- Manova(my.model, type="III")
summary(fit, test = "Wilks")
FemaleEtaSq <- etasq(my.model, test = "Wilks", anova = TRUE, partial = TRUE)
FemaleEtaSq
################################################################

###Extract EMMs for multivariate visualization
EMM_grid <- ref_grid(my.model, mult.name = "Comps")
EMM_grid
EMMs <- emmeans(my.model, ~ Sulf*Drainage | rep.meas)
EMMs
#write.csv(EMMs,"FemaleGenital_All_EMMs_SulfXDrainage.csv") # export the pc loadings

################################################################
###Post-hoc pairwise analyses
###The following code separates the data file into separate objects by Drainage
sites<-unique(data2$Drainage); n<-length(sites)
rhs<-paste("data2[data2$Drainage==sites[",1:n,"],1:ncol(data2)]",sep="")
eq<-paste(paste(sites,rhs,sep="<-"),collapse=";")
eval(parse(text=eq))

sites			#displays the names for each population/sample

#Taco
my.model <- lm(cbind(Comp.1, Comp.2, Comp.3) 
               ~ logMass + Sulf + logMass*Sulf, data=Taco)
fit <- Manova(my.model, type="III")
summary(fit, test = "Wilks")

#Puya
my.model <- lm(cbind(Comp.1, Comp.2, Comp.3) 
               ~ logMass + Sulf + logMass*Sulf, data=Puya)
fit <- Manova(my.model, type="III")
summary(fit, test = "Wilks")

#Ixta- non-sig intxn
my.model <- lm(cbind(Comp.1, Comp.2, Comp.3) 
               ~ logMass + Sulf, data=Ixta)
fit <- Manova(my.model, type="III")
summary(fit, test = "Wilks")

#Pich- non-sig intxn
my.model <- lm(cbind(Comp.1, Comp.2, Comp.3) 
               ~ logMass + Sulf, data=Pich)
fit <- Manova(my.model, type="III")
summary(fit, test = "Wilks")


################################################################
################################################################

###Plot Female EMMs for PC1 and PC2###
FemaleEMMs <- read.csv("FemaleGenital_RelativeMeasures_EMMs_PCA.csv", header=TRUE, row.names=1)
names(FemaleEMMs)

FemaleEMM_plot <- ggplot(data = FemaleEMMs,aes(x = emmean.FemComp1,y = emmean.FemComp2, shape=Drainage)) +
  geom_errorbar(aes(ymin = lower.SE.Comp2,ymax = upper.SE.Comp2, width=.005)) + 
  geom_errorbarh(aes(xmin = lower.SE.Comp1,xmax = upper.SE.Comp1, height=0.005)) + 
  scale_shape_manual(values=c(23, 22, 21, 24)) + geom_line(color = 'grey', linetype=2) + theme_classic() + theme(text = element_text(size=30), axis.text=element_text(colour = "black"), legend.position="none") +
  geom_point(data = subset(FemaleEMMs, Sulf != "S"), color = 'black', fill='dodgerblue3', size = 8) +
  geom_point(data = subset(FemaleEMMs, Sulf == "S"), color = 'black' , fill='yellow1', size = 8)
FemaleEMM_plot


####Calculate divergence vector scores in Excel then import for data visualization
FemaleData <- read.csv("FemaleGenitalDivergenceVectorScores.csv", header=TRUE, row.names=1)
FemaleData 

#BOXPLOTS
plot <- ggplot(FemaleData, aes(x=Drainage, y=DVS, fill=Sulf)) + 
  geom_boxplot(color="black", width= .9, position=position_dodge(.9))
plot <- plot + coord_flip()
plot <- plot + scale_x_discrete(limits=c("Pich", "Ixta", "Puya", "Taco")) +
  scale_fill_manual(values=c("dodgerblue3", "yellow1")) + theme(
    panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA, colour = "black", size = 1.5), legend.position="none") + theme(text = element_text(size=30), axis.text=element_text(colour = "black"))
plot

