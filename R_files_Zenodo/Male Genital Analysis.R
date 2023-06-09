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


data <- read.csv("Male Genitalia Data.csv", header=TRUE, row.names=1)

names(data)

######################
#######  PCA  ########
######################
##Shape Data 

pca<-princomp(data[,8:61])

summary(pca) # print variance accounted for by each axis
pca$loadings #loadings, not sure why it only prints the big ones
#write.csv(loadings(pca),"Male_Genital_PCA_loadings.csv") # export the pc loadings
plot(pca,type="lines") # scree plot
pca$scores # the principal components (scores)
biplot(pca) # print the scores and the loadings
#write.csv(pca$scores,"Male_Genital_PCA_Scores.csv") # export the pc loadings
ev <- eigen(cor(data[,8:61])) # get eigenvalues
ev$values # print eigenvalues for each axis

scores <- pca$scores
scores <- data.frame(scores)
data2 <- cbind(data, scores)

names(data2)

###############  MANCOVA on Significant PC Axes  ###############
################################################################
my.model <- lm(cbind(Comp.1, Comp.2, Comp.3,Comp.4, Comp.5, Comp.6,
                     Comp.7, Comp.8, Comp.9, Comp.10, Comp.11) 
               ~ logMass + Sulf + Drainage + logMass*Sulf + logMass*Drainage + Sulf*Drainage, data=data2)
fit <- Manova(my.model, type="III")
summary(fit, test = "Wilks")
MaleEtaSq <- etasq(my.model, test = "Wilks", anova = TRUE, partial = TRUE)
MaleEtaSq
################################################################

###pull out EMMs for Sulf*Drainage term
EMM_grid <- ref_grid(my.model, mult.name = "Comps")
EMM_grid
EMMs <- emmeans(my.model, ~ Sulf*Drainage | rep.meas)
EMMs
#write.csv(EMMs,"MaleGenital_All_EMMs_SulfXDrainage.csv") # export the EMMS


################################################################

###Plot male EMMs for PC1 and PC2###
MaleEMMs <- read.csv("MaleGenitalEMMs.csv", header=TRUE, row.names=1)
names(MaleEMMs)

MaleEMM_plot <- ggplot(data = MaleEMMs,aes(x = emmean.Comp1,y = emmean.Comp2, shape=Drainage)) +
  geom_errorbar(aes(ymin = lower.SE.Comp2,ymax = upper.SE.Comp2, width=.0025)) + 
  geom_errorbarh(aes(xmin = lower.SE.Comp1,xmax = upper.SE.Comp1, height=0.0025)) + 
  scale_shape_manual(values=c(23, 22, 21, 24)) + geom_line(color = 'grey', linetype=2) + theme_classic() + theme(text = element_text(size=30), axis.text=element_text(colour = "black"), legend.position="none") +
  geom_point(data = subset(MaleEMMs, Sulf != "S"), color = 'black', fill='dodgerblue3', size = 8) +
  geom_point(data = subset(MaleEMMs, Sulf == "S"), color = 'black' , fill='yellow1', size = 8)
MaleEMM_plot

################################################################



