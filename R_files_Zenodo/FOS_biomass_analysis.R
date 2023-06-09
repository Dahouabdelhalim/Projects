#now compare these data to the biomass and height data from FOS dataset


#open packages
library(dplyr)
library(ggplot2)
library(lme4)
library(multcomp)
library(scales)
library(extrafont)
loadfonts(device="win")

#clear workspace
rm(list=ls())

#import data - import data
combined_data #import "FOS_biomass_analysis.csv" as "combined_data
str(combined_data)

#################################
#create a continent variable
levels(as.factor(combined_data$CountryName))
combined_data$Continent<-combined_data$CountryName
inter1<-gsub("Bolivia","Americas",combined_data$Continent)
inter2<-gsub("Brazil", "Americas",inter1)
inter3<-gsub("Cameroon", "Africa",inter2)
inter4<-gsub("Congo, Democratic republic", "Africa",inter3)
inter5<-gsub("Costa Rica", "Americas",inter4)
inter6<-gsub("French Guiana", "Americas",inter5)
inter7<-gsub("Gabon", "Africa",inter6)
inter8<-gsub("Ghana", "Africa",inter7)
inter9<-gsub("Liberia", "Africa",inter8)
inter10<-gsub("Malaysia", "Asia",inter9)
inter11<-gsub("Panama", "Americas",inter10)
inter12<-gsub("Peru", "Americas",inter11)
inter13<-gsub("Guyana", "Americas",inter12)
levels(as.factor(inter13))
combined_data$Continent<-as.factor(inter13)

#try averaging values for each "Plot_ID" value
agg_data<-combined_data%>%
  group_by(Continent,Plot_ID)%>%
  summarise(CG.FRD = mean(CG...10kA.FRD, na.rm = TRUE),
            log.CG.FRD = mean(log.CG.FRD, na.rm = TRUE),
            AGB_Chave.mean = mean(AGB_Chave, na.rm = TRUE),
            precip.mm.mean = mean(precip.mm, na.rm = TRUE),
            temp.C.mean = mean(temp.C, na.rm = TRUE))
str(agg_data)

#is lightning frequency correlated with climatic variables?
cor.test(agg_data$log.CG.FRD,agg_data$precip.mm.mean)
cor.test(agg_data$log.CG.FRD,agg_data$temp.C.mean)


#also scale all of the variables
agg_data$scaled_log.CG<-as.numeric(scale(agg_data$log.CG.FRD))
agg_data$scaled_AGB_Chave<-as.numeric(scale(agg_data$AGB_Chave.mean))
agg_data$scaled_precip.mm<-as.numeric(scale(agg_data$precip.mm.mean))
agg_data$scaled_temp.C<-as.numeric(scale(agg_data$temp.C.mean))

#no biomass data for Australia, so drop it here too (we don't use Australia in the other analyses anyway)
agg_data<-agg_data[agg_data$Continent!="Australia",]

#test the simple model we used for the other approaches
reg1<-lm(AGB_Chave.mean ~ scaled_log.CG + scaled_precip.mm + scaled_temp.C + Continent + Continent:log.CG.FRD, data = agg_data)
summary(reg1)

reg2<-lm(AGB_Chave.mean ~ scaled_log.CG + scaled_precip.mm + scaled_temp.C + Continent, data = agg_data)
AIC(reg1,reg2)#no effect of the interaction

reg3<-lm(AGB_Chave.mean ~ scaled_log.CG + scaled_precip.mm + Continent, data = agg_data)
AIC(reg2,reg3)#no effect of temperature

reg4<-lm(AGB_Chave.mean ~ scaled_log.CG + Continent, data = agg_data)
AIC(reg3,reg4)#no effect of precip

reg5<-lm(AGB_Chave.mean ~ scaled_log.CG , data = agg_data)
AIC(reg4,reg5)#no effect of continent
#best fit model here
summary(reg5)
plot(reg5)#imperfect, but overall it is it quite decent
hist(resid(reg5))

reg6<-lm(AGB_Chave.mean ~ 1, data = agg_data)
AIC(reg5,reg6)#moderately strong effect of lightning frequency

###############################
#create figure representing best fit model

#create data for 
#predict values from the model and compare to base formula
pred_bestmod<-confint(reg5)
pred_bestmod

#create a range of lightning-based values
lightning_range<-with(combined_data,seq(min(log.CG.FRD),max(log.CG.FRD),.01))

#create data for the line and the confidence intervals
#first create the best fit line
summary(reg5)$coefficients
predicted_curve<-summary(reg5)$coefficients[1] + summary(reg5)$coefficients[2]*lightning_range

#then add the hi and low CI data
pred_bestmod
loCI<-pred_bestmod[1] + pred_bestmod[2]*lightning_range
hiCI<-pred_bestmod[3] + pred_bestmod[4]*lightning_range


#combine these vectors into a data frame
model_plotting<-data.frame(predicted_curve,loCI,hiCI,lightning_range)
str(model_plotting)

#I generally use ggplot2 to visualize data
theme_basis<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(axis.text.y=element_text(family = "Arial", 
                                 colour="black", face ="bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial", colour="black", face ="bold", size = 12))

#create axis breaks on a log scale
breaks_normal<-c(seq(.05,.09,.01),seq(.1,.9,.1),seq(1,10,1),20)
log_breaks<-log(breaks_normal)
x_labels<-c(rep("",5),"0.1",rep("",8),"1.0",rep("",8),"10.0","20")


#now plot these relationships
biomass_ENTLN<-ggplot(data = agg_data, aes(x = log.CG.FRD, y = AGB_Chave.mean)) + theme_basis +
  geom_point(size = 2, aes(color = Continent))+
  scale_color_manual(values=c("#F8766D","gold2","#00BFC4"))+
  geom_smooth(method = "lm",color = "grey20", fill = "grey", alpha = .4)+
  theme(legend.position="none")+
  scale_x_continuous(name= expression(bold("CG Lightning Frequency"~(fl~km^{-2}~yr^{-1}))),breaks=log_breaks,
                     labels = x_labels, expand = c(0,0.0), limits = c(log(0.04),log(25)))+
  scale_y_continuous(name = expression(bold("Aboveground Biomass "~(Mg~ha^{-1}))),
                     breaks=seq(0,800,100),labels=c("0","","200","","400","","600","","800"))

biomass_ENTLN#the 20 missing rows have no biomass data

ggsave("biomass_ENTLN.tiff",biomass_ENTLN,dpi = 600, width = 3.25,height = 2.2,scale = 1.9, compression = "lzw")



