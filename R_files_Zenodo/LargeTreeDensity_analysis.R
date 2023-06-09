#do the same thing I did for biomass and biomass turnover, but do it with large tree abundance

#open packages
#open packages
library(relaimpo)
library(dplyr)
library(ggplot2)
library(lme4)
library(MuMIn)
library(multcomp)
library(scales)
library(extrafont)
loadfonts(device="win")

#clear workspace
rm(list=ls())

#import tree abundance data from Slik et al. 2013 - file name is "LargeTreeDensity_data.csv"
combined_data #<- save "LargeTreeDensity_data" as "combined_data"
str(combined_data)

#look for correlations between lightning and environmental variables
cor.test(combined_data$log.CG.FRD,combined_data$temp.C)
cor.test(combined_data$log.CG.FRD,combined_data$precip.mm)

####################################################
#analyze tree density with one interaction (FRD:Continent) and then just other climatic variables

#run the full model
mod1<-lm(stems.ov70~log.CG.FRD + Continent + temp.C + precip.mm + CG.FRD:Continent, data = combined_data, na.action = "na.fail")
summary(mod1)
AIC(mod1)

#drop the interaction
mod2<-lm(stems.ov70~log.CG.FRD + Continent + temp.C + precip.mm, data = combined_data, na.action = "na.fail")
AIC(mod1,mod2)
#interaction is not significant

#drop temp
mod3<-lm(stems.ov70~log.CG.FRD + Continent + precip.mm, data = combined_data, na.action = "na.fail")
AIC(mod2,mod3)
#temp is not significant

#drop precip
mod4<-lm(stems.ov70~log.CG.FRD + Continent, data = combined_data, na.action = "na.fail")
AIC(mod3,mod4)
#precip is not significant (only marginally though)
summary(mod4)
plot(mod4)#model fit is great
hist(resid(mod4))#model residuals are normally distributed

#calculate relative contributions of different parameters to the full model
calc.relimp(mod4)

#drop continent
mod5<-lm(stems.ov70~log.CG.FRD, data = combined_data, na.action = "na.fail")
AIC(mod4,mod5)
#continent is very significant!!!

#drop FRD
mod6<-lm(stems.ov70~Continent, data = combined_data, na.action = "na.fail")
AIC(mod4,mod6)
#FRD is very significant!!!


###########################################################
#predict values from the model and compare to base formula
pred_bestmod<-confint(mod4)
pred_bestmod

#create a range of lightning-based values
lightning_range_Africa<-with(combined_data[combined_data$Continent=="Africa",],seq(min(log.CG.FRD),max(log.CG.FRD),.01))
lightning_range_Americas<-with(combined_data[combined_data$Continent=="America",],seq(min(log.CG.FRD),max(log.CG.FRD),.01))
lightning_range_Asia<-with(combined_data[combined_data$Continent=="Asia",],seq(min(log.CG.FRD),max(log.CG.FRD),.01))

#now add the best fit line
summary(mod4)$coefficients
predicted_Africa<-summary(mod4)$coefficients[1] + summary(mod4)$coefficients[2]*lightning_range_Africa
predicted_Americas<-(summary(mod4)$coefficients[1] + summary(mod4)$coefficients[3]) + summary(mod4)$coefficients[2]*lightning_range_Americas
predicted_Asia<-(summary(mod4)$coefficients[1] + summary(mod4)$coefficients[4]) + summary(mod4)$coefficients[2]*lightning_range_Asia

#then add the hi and low CI data
pred_bestmod
loCI_Africa<-pred_bestmod[1] + pred_bestmod[2]*lightning_range_Africa
hiCI_Africa<-pred_bestmod[5] + pred_bestmod[6]*lightning_range_Africa
loCI_Americas<-(pred_bestmod[1] + pred_bestmod[3]) + pred_bestmod[2]*lightning_range_Americas
hiCI_Americas<-(pred_bestmod[5] + pred_bestmod[7]) + pred_bestmod[6]*lightning_range_Americas
loCI_Asia<-(pred_bestmod[1] + pred_bestmod[4]) + pred_bestmod[2]*lightning_range_Asia
hiCI_Asia<-(pred_bestmod[5] + pred_bestmod[8]) + pred_bestmod[6]*lightning_range_Asia


#this might be easier with continent-specific curves
model_plotting_Africa<-data.frame(predicted_Africa,loCI_Africa,hiCI_Africa,lightning_range_Africa,rep("Africa",length(predicted_Africa)))
colnames(model_plotting_Africa)<-c("predicted_curve","loCI","hiCI","lightning_range","Continent")
model_plotting_Americas<-data.frame(predicted_Americas,loCI_Americas,hiCI_Americas,lightning_range_Americas,rep("Americas",length(predicted_Americas)))
colnames(model_plotting_Americas)<-c("predicted_curve","loCI","hiCI","lightning_range","Continent")
model_plotting_Asia<-data.frame(predicted_Asia,loCI_Asia,hiCI_Asia,lightning_range_Asia,rep("Asia",length(predicted_Asia)))
colnames(model_plotting_Asia)<-c("predicted_curve","loCI","hiCI","lightning_range","Continent")

#create x axis breaks and labels
breaks_normal<-c(seq(0.05,0.09,0.01),seq(.1,1,.1),seq(2,10,1),20,30,40)
breaks_log<-log(breaks_normal)
x_labels<-c(rep("",5),"0.1",rep("",8),"1.0",rep("",8),"10","20","30","40")

#create plot theme
theme_basis<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(axis.text.y=element_text(family = "Arial", 
                                 colour="black", face ="bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial", colour="black", face ="bold", size = 12))

################################################################
#now create the same figure with an expanded x-axis

#now plot the data using profile likelihood confidence intervals
treedensity_ENTLN_expandedx<-ggplot(data=combined_data,aes(x = combined_data$log.CG.FRD,y=combined_data$stems.ov70,color=combined_data$Continent))+
  geom_ribbon(data = model_plotting_Africa,color = "grey90", fill = "grey", alpha = .4, aes(x = lightning_range, y = predicted_curve, 
                                                                                            ymax = hiCI, ymin = loCI))+
  geom_ribbon(data = model_plotting_Americas,color = "grey90", fill = "grey", alpha = .4, aes(x = lightning_range, y = predicted_curve, 
                                                                                              ymax = hiCI, ymin = loCI))+
  geom_ribbon(data = model_plotting_Asia,color = "grey90", fill = "grey", alpha = .4, aes(x = lightning_range, y = predicted_curve, 
                                                                                          ymax = hiCI, ymin = loCI))+
  geom_line(data = model_plotting_Americas, size = 1,color = "gold2", aes(x = lightning_range,y = predicted_curve))+
  geom_line(data = model_plotting_Asia, size = 1,color = "#00BFC4", aes(x = lightning_range,y = predicted_curve))+
  geom_line(data = model_plotting_Africa, size = 1,color = "#F8766D", aes(x = lightning_range,y = predicted_curve))+
  geom_point(size=1.5, aes(color=combined_data$Continent))+
  scale_color_manual(values=c("#F8766D","gold2","#00BFC4"))+
  theme(legend.position="none")+
  theme_basis+
  scale_x_continuous(name= expression(bold("CG Lightning Frequency"~(fl~km^{-2}~yr^{-1}))),
                     breaks=breaks_log,labels = x_labels,expand = c(0,0.0), limits = c(log(0.04),log(25)))+
  scale_y_continuous(name = expression(bold("Large Tree Density"~(trees~ha^{-1}))),
                     breaks=seq(0,30,5),labels=c("0","","10","","20","","30"))+
  coord_cartesian(ylim = c(0,31))

treedensity_ENTLN_expandedx

ggsave("treedensity_ENTLN.tiff",treedensity_ENTLN_expandedx,dpi = 600, width = 3.25,height = 2.2,scale = 1.9, compression = "lzw")



###################################################################################
#lets now repeat the model averaging procedure of Slik et al. 2013 including all covariates
str(combined_data)

#ztransform all of the variables
combined_data$scaled_stems.ov70<-scale(combined_data$stems.ov70)
combined_data$scaled_log.CG.FRD<-scale(combined_data$log.CG.FRD)
combined_data$scaled_Soil1<-scale(combined_data$Soil1)
combined_data$scaled_Soil2<-scale(combined_data$Soil2)
combined_data$scaled_Soil3<-scale(combined_data$Soil3)
combined_data$scaled_Soil4<-scale(combined_data$Soil4)
combined_data$WD<-scale(combined_data$WD..basal.area.weighted.)
combined_data$ECM<-scale(combined_data$ECM..biomass.weighted..proportion.)
combined_data$wind<-scale(combined_data$Wind.dispersal..biomass.weighted.proportion.)
combined_data$warmestmonth<-scale(combined_data$Bio_05_max_temp_warmest_month)
combined_data$coldestmonth<-scale(combined_data$Bio_06_min_temp_coldest_month)
combined_data$wettestmonth<-scale(combined_data$Bio_13_prec_wettest_month)

#create full model
mod1_slik_vars<-lm(scaled_stems.ov70~ scaled_Soil1 + scaled_Soil2 + scaled_Soil3 + scaled_Soil4 + 
                     WD + ECM + wind+warmestmonth + coldestmonth + wettestmonth,data = combined_data, na.action = na.fail)
summary(mod1_slik_vars)
AICc(mod1_slik_vars)

#do model averaging
d1<-dredge(mod1_slik_vars)

#model averaging
modave_d1<-model.avg(d1)
summary(modave_d1)
modave_d1$coefficients#the coefficients exactly match those reported by Slik et al.
importance(modave_d1)#importance values also match exactly

#calculate relative contributions of different parameters to the full model
calc.relimp(mod1_slik_vars)#this indicates that the second most important parameter in the model is lightning frequency

######################################
#so let us run the same model, but with lightning
#create full model
mod2_slik_vars<-lm(scaled_stems.ov70~scaled_log.CG.FRD + scaled_Soil1 + scaled_Soil2 + scaled_Soil3 + scaled_Soil4 + 
                     WD + ECM + wind+warmestmonth + coldestmonth + wettestmonth,data = combined_data, na.action = na.fail)
summary(mod2_slik_vars)
AICc(mod2_slik_vars)

#do model averaging
d2<-dredge(mod2_slik_vars)
summary(d2)

#model averaging
modave_d2<-model.avg(d2)
summary(modave_d2)
modave_d2$coefficients
importance(modave_d2)

#calculate relative contributions of different parameters to the full model
calc.relimp(mod2_slik_vars)#this indicates that the second most important parameter in the model is lightning frequency
#relative importance values (percent of variation explained by each variable)

