#look at a associations between turnover and lightning flash rate density
#consider the effects of continent
##I would not be surprised if continents differ either because of negative feedbacks or adaptations
#open packages
library(relaimpo)
library(dplyr)
library(ggplot2)
library(lme4)
library(multcomp)
library(scales)
library(extrafont)
loadfonts(device="win")

#clear workspace
rm(list=ls())

#import dataset with Galbraith et al. 2013 data on woody biomass turnover rates
combined_data #<- import "Woodyturnover.csv" file here
str(combined_data)

###########################################

#construct the full model for analysis
ENTLNtrans_mod1<-lm(turnover_rate~log_ENTLN+Continent+temp.C+precip.mm+log_ENTLN:Continent,data=combined_data)
summary(ENTLNtrans_mod1)
AIC(ENTLNtrans_mod1)

#drop the lightning:continent interaction
ENTLNtrans_mod2<-lm(turnover_rate~log_ENTLN+Continent+temp.C+precip.mm,data=combined_data)
AIC(ENTLNtrans_mod1,ENTLNtrans_mod2)
#interaction is not significant

#drop temp
ENTLNtrans_mod3<-lm(turnover_rate~log_ENTLN+Continent+precip.mm,data=combined_data)
AIC(ENTLNtrans_mod2,ENTLNtrans_mod3)
#temp is not significant

#drop precipitation
ENTLNtrans_mod4<-lm(turnover_rate~log_ENTLN+Continent,data=combined_data)
AIC(ENTLNtrans_mod3,ENTLNtrans_mod4)
summary(ENTLNtrans_mod4)
#precip is not significant
plot(ENTLNtrans_mod4)
hist(resid(ENTLNtrans_mod4))
#some pretty good looking residuals

#calculate relative contributions of different parameters to the full model
calc.relimp(ENTLNtrans_mod4)

#drop continent
ENTLNtrans_mod5<-lm(turnover_rate~log_ENTLN,data=combined_data,na.action = na.fail)
AIC(ENTLNtrans_mod4,ENTLNtrans_mod5)
#continent is very significant

#drop lightning frequency
ENTLNtrans_mod6<-lm(turnover_rate~Continent,data=combined_data,na.action = na.fail)
AIC(ENTLNtrans_mod4,ENTLNtrans_mod6)
#lightning is very significant


########################################################
#make publishable figures for just the log-transformed ENTLN data
#create basic theme
theme_basis<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(axis.text.y=element_text(family = "Arial", 
                                 colour="black", face ="bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial", colour="black", face ="bold", size = 12))

###########################################################
#predict values from the model and compare to base formula
pred_bestmod<-confint(ENTLNtrans_mod4)
pred_bestmod

#create a range of lightning-based values
lightning_range_Africa<-with(combined_data[combined_data$Continent=="Africa",],seq(min(log_ENTLN),max(log_ENTLN),.01))
lightning_range_Americas<-with(combined_data[combined_data$Continent=="Americas",],seq(min(log_ENTLN),max(log_ENTLN),.01))
lightning_range_Asia<-with(combined_data[combined_data$Continent=="Asia",],seq(min(log_ENTLN),max(log_ENTLN),.01))

#now add the best fit line for each continent
summary(ENTLNtrans_mod4)$coefficients
predicted_Africa<-summary(ENTLNtrans_mod4)$coefficients[1] + summary(ENTLNtrans_mod4)$coefficients[2]*lightning_range_Africa
predicted_Americas<-(summary(ENTLNtrans_mod4)$coefficients[1] + summary(ENTLNtrans_mod4)$coefficients[3]) + summary(ENTLNtrans_mod4)$coefficients[2]*lightning_range_Americas
predicted_Asia<-(summary(ENTLNtrans_mod4)$coefficients[1] + summary(ENTLNtrans_mod4)$coefficients[4]) + summary(ENTLNtrans_mod4)$coefficients[2]*lightning_range_Asia

#then add the hi and low CIs
pred_bestmod
loCI_Africa<-pred_bestmod[1] + pred_bestmod[2]*lightning_range_Africa
hiCI_Africa<-pred_bestmod[5] + pred_bestmod[6]*lightning_range_Africa
loCI_Americas<-(pred_bestmod[1] + pred_bestmod[3]) + pred_bestmod[2]*lightning_range_Americas
hiCI_Americas<-(pred_bestmod[5] + pred_bestmod[7]) + pred_bestmod[6]*lightning_range_Americas
loCI_Asia<-(pred_bestmod[1] + pred_bestmod[4]) + pred_bestmod[2]*lightning_range_Asia
hiCI_Asia<-(pred_bestmod[5] + pred_bestmod[8]) + pred_bestmod[6]*lightning_range_Asia

#combine these vectors to make a dataframe for plotting
predicted_curve<-c(predicted_Africa,predicted_Americas,predicted_Asia)
loCI<-c(loCI_Africa,loCI_Americas,loCI_Asia)
hiCI<-c(hiCI_Africa,hiCI_Americas,hiCI_Asia)
Continent<-c(rep("Africa",length(predicted_Africa)),rep("Americas",length(predicted_Americas)),rep("Asia",length(predicted_Asia)))
lightning_range<-c(lightning_range_Africa,lightning_range_Americas,lightning_range_Asia)
model_plotting<-data.frame(predicted_curve,loCI,hiCI,lightning_range,Continent)
str(model_plotting)

#this might be easier with continent-specific curves
model_plotting_Africa<-data.frame(predicted_Africa,loCI_Africa,hiCI_Africa,lightning_range_Africa,rep("Africa",length(predicted_Africa)))
colnames(model_plotting_Africa)<-c("predicted_curve","loCI","hiCI","lightning_range","Continent")
model_plotting_Americas<-data.frame(predicted_Americas,loCI_Americas,hiCI_Americas,lightning_range_Americas,rep("Americas",length(predicted_Americas)))
colnames(model_plotting_Americas)<-c("predicted_curve","loCI","hiCI","lightning_range","Continent")
model_plotting_Asia<-data.frame(predicted_Asia,loCI_Asia,hiCI_Asia,lightning_range_Asia,rep("Asia",length(predicted_Asia)))
colnames(model_plotting_Asia)<-c("predicted_curve","loCI","hiCI","lightning_range","Continent")

#create x axis breaks and labels
breaks_normal<-c(seq(.05,.09,.01),seq(.1,1,.1),seq(2,10,1),20,30)
breaks_log<-log10(breaks_normal)
x_labels<-c(rep("",5),"0.1",rep("",8),"1.0",rep("",8),"10","20","30")

#############################################################################################
#create figure with a corrected axis

#now plot the data using profile likelihood confidence intervals
turnover_ENTLN_expandedx<-ggplot(data=combined_data,aes(x = combined_data$log_ENTLN,y=combined_data$turnover_rate,color=combined_data$Continent))+
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
  scale_x_continuous(name= expression(bold("CG Lightning Frequency"~(fl~km^{-2}~yr^{-1}))),breaks=breaks_log,
                     labels = x_labels, expand = c(0,0.0), limits = c(log(0.04,base = 10),log(25,base = 10)))+
  scale_y_continuous(name = expression(bold("Woody Biomass Turnover Rate"~('%'~yr^{-1}))),
                     breaks=seq(0.01,0.045,0.005),labels=c("1","","2","","3","","4",""),expand = c(0,0.001))

turnover_ENTLN_expandedx

ggsave("turnover_ENTLN.tiff",turnover_ENTLN_expandedx,dpi = 600, width = 3.25,height = 2.2,scale = 1.9, compression = "lzw")


