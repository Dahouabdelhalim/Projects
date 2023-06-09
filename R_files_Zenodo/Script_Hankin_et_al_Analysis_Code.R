# Script_Hankin_et_al_Analysis_Code.R
#
# This script will reproduce the results (tables and figures) presented in 
# Hankin et al. (2018).
#
# Hankin, L.E., Higuera, P.E., Davis, K.T., Dobrowski, S.Z. 2018. Accuracy
# of node and bud-scar counts for aging two dominant conifers in western
# North America. Forest Ecology and Management. In Press.
#
# Forest Ecology and Management DOI: https://doi.org/10.1016/j.foreco.2018.06.001
# Dryad Repository DOI: https://doi.org/10.5061/dryad.63m17n5
#
# Data File Requirements: 
# (1) Data_TreeAges.csv - all tree age data for running the code
# (2) Data_AgeTool.csv - data for running the tool alone
#
# Created by: L.E. Hankin
# Created: June 2018
# University of Montana, PaleoEcology and Fire Ecology Lab
#
# Contact:
# Lacey Hankin: lacey.hankin@gmail.com
# Philip Higuera: phiguera@umontana.edu
#
# Set working directory ---------------------------------------------------
setwd("L:/4_archivedData/Hankin_et_al_2018_FEM") # The directory should be changed  
#                                           to reflect your data file location

# Data Setup --------------------------------------------------------------

# Load required packages
# If packages have not been previously installed, uncomment the following code,
# or install these packages manually:
#install.packages("ggplot2","data.table","plyr","reshape2","car","dplyr",
#      "lme4","lmerTest","piecewiseSEM","merTools","foreach","scales",
#      "cowplot","nlme","effects")
require(ggplot2);require(data.table);require(plyr);require(reshape2);
require(car);require(dplyr);require(lme4);require(lmerTest);
require(piecewiseSEM);require(merTools);require(foreach);require(scales);
require(cowplot);require(nlme);require(effects)

# Load in the data
All<-read.csv("Data_TreeAges.csv")

# Remove NAs
All<-All[!is.na(All$AgeDif),]
All<-All[!is.na(All$Nodes),]
All<-All[!is.na(All$Age),]

# Remove count data that did not meet confidence criteria and are prior to fires
All<-All[!(All$Confidence==2),]
All<-All[!(All$TimeSinceFire<0),]

# Establish categories of precision to evaluate bias between nodes and ring counts
All$zero_yr<-ifelse(All$AgeDif==0,1,0)
All$one_yr<-ifelse(abs(All$AgeDif)<=1,1,0)
All$five_yr<-ifelse(abs(All$AgeDif)<=5,1,0)
All$two_yr<-ifelse(abs(All$AgeDif)<=2,1,0)

# Subset data per region
NRall<-subset(All,Region=="NR")
CAall<-subset(All,Region=="CA")
SWall<-subset(All,Region=="SW")
NRall_2017<-subset(NRall,Sample.year==2017)

# Create species-specific subsets of data
NRpipo<-subset(NRall,Species=="PIPO")
NRpsme<-subset(NRall,Species=="PSME")
SWpipo<-subset(SWall,Species=="PIPO")
SWpsme<-subset(SWall,Species=="PSME")
CApipo<-subset(CAall,Species=="PIPO")
CApsme<-subset(CAall,Species=="PSME")
Allpipo<-subset(All,Species=="PIPO")
Allpsme<-subset(All,Species=="PSME")
NRpipo_2017<-subset(NRall_2017,Species=="PIPO")
NRpsme_2017<-subset(NRall_2017,Species=="PSME")

# Establish categories of precision to evaluate bias in subset of data used
  # for bud-scar evaluation
NRall_2017$zero_yr_node<-ifelse(NRall_2017$AgeDif==0,1,0)
NRall_2017$one_yr_node<-ifelse(abs(NRall_2017$AgeDif)<=1,1,0)
NRall_2017$five_yr_node<-ifelse(abs(NRall_2017$AgeDif)<=5,1,0)
NRall_2017$two_yr_node<-ifelse(abs(NRall_2017$AgeDif)<=2,1,0)

NRall_2017$zero_yr_bud<-ifelse(NRall_2017$AgeDif_Buds==0,1,0)
NRall_2017$one_yr_bud<-ifelse(abs(NRall_2017$AgeDif_Buds)<=1,1,0)
NRall_2017$five_yr_bud<-ifelse(abs(NRall_2017$AgeDif_Buds)<=5,1,0)
NRall_2017$two_yr_bud<-ifelse(abs(NRall_2017$AgeDif_Buds)<=2,1,0)

## Age vs. Nodes analysis ---------------------------------------------------

# Fit linear mixed effects models WITHOUT non-constant variance structure 
  # to establish baseline for each region (Zuur et al. 2009)
NRme<-lmer(Nodes~Age*Species+(1|Site),data=NRall)
summary(NRme)
Anova(NRme)
confint(NRme)

SWme<-lmer(Nodes~Age*Species+(1|Site),data=SWall)
summary(SWme)
Anova(SWme)
confint(SWme)

CAme<-lmer(Nodes~Age*Species+(1|Site),data=CAall)
summary(CAme)
Anova(CAme)
confint(CAme)

Allme<-lmer(Nodes~Age*Species+(1|Site),data=All)
summary(Allme)
Anova(Allme)
confint(Allme)

# Add exponential/power variance structure using generalized least 
  # squares method (Zuur et al. 2009)
head(NRall)
sapply(NRall, function(x) sum(is.na(x)))
NRall<-NRall[is.na(NRall$Age)==F,] # Get rid of two rows with NA for age
NRall<-NRall[is.na(NRall$Nodes)==F,] # Get rid of 6 rows with NA for nodes (gls
  # doesn't work with missing data like lmer will)
NRme<-lmer(Nodes~Age*Species+(1|Site),data=NRall)

# Evaluate different variance structures for each region

# Northern Rockies
NRgls<-gls(Nodes~Age*Species,data=NRall) # Linear model no random effects or 
  # variance structure to compare to variance structures
vf1Fixed <- varFixed(~Age) # Specify that var increases with age (in a 
  # multiplicative way)
NRgls1<-gls(Nodes~Age*Species,weights = vf1Fixed,data=NRall) 
summary(NRgls1)
# plot(NRgls1)
vf3 <- varPower(form =~ Age) # Variance varies as a power of the covariate
NRgls3<-gls(Nodes~Age*Species,weights = vf3,data=NRall) 
vf4 <- varExp(form =~ Age) # Exponential of the variance covariate
NRgls4<-gls(Nodes~Age*Species,weights = vf4,data=NRall)
AIC(NRgls,NRgls1,NRgls3,NRgls4) # Variance structure in gls4 is best
summary(NRgls4)
# plot(NRgls)
# plot(NRgls4) 
E1 <- resid(NRgls4)
# plot(E1 ~ Age,
#         ylab = "Ordinary residuals", data = NRall) # These are still cone 
  # shaped, need to normalize by variance
E2 <- resid(NRgls4, type = "normalized")
# plot(E2 ~ Age,
#     ylab = "Normalized residuals", data = NRall) # This is where the pattern 
  # should be gone, which it is for the most part so this worked

# Now add nestedness (site random effect)
# Package nlme and function lme allows you to include the weights call/variance 
  # structure we defined above:
NRme1<-lme(Nodes~Age*Species,random=~1|Site,data=NRall)
summary(NRme1)
# plot(NRme1)
NRme2<-lme(Nodes~Age*Species,random=~1|Site,weights = vf4,data=NRall)
summary(NRme2)
# plot(NRme2)
summary(NRme)
Anova(NRme)
confint(NRme)
intervals(NRme2) 
NRsq<-sem.model.fits(NRme2)

# Evaluate R2 as the predicted values as a function of observed values
obsnr<-NRall$Nodes
prednr<-fitted(NRme2)
ponr<-lm(prednr~obsnr)
nrsq<-summary(ponr)$r.squared

# California
head(CAall)
sapply(CAall, function(x) sum(is.na(x)))
CAall<-CAall[is.na(CAall$Age)==F,] # Get rid of two rows with NA for age
CAall<-CAall[is.na(CAall$Nodes)==F,] # Get rid of 6 rows with NA for nodes (gls 
  # doesn't work with missing data like lmer will)
CAme<-lmer(Nodes~Age*Species+(1|Site),data=CAall)

CAgls<-gls(Nodes~Age*Species,data=CAall) # Linear model no random effects or 
  # variance structure to compare to variance structures
vf1Fixed <- varFixed(~Age) # Specify that var increases with age (in a 
  # multiplicative way)
CAgls1<-gls(Nodes~Age*Species,weights = vf1Fixed,data=CAall) 
summary(CAgls1)
# plot(CAgls1)
vf3 <- varPower(form =~ Age) # Variance varies as a power of the covariate
CAgls3<-gls(Nodes~Age*Species,weights = vf3,data=CAall) 
vf4 <- varExp(form =~ Age) # Exponential of the variance covariate
CAgls4<-gls(Nodes~Age*Species,weights = vf4,data=CAall)
AIC(CAgls,CAgls1,CAgls3,CAgls4) # Variance structure in gla2 and gls4 is best, 
  # stay consistent with exponential
summary(CAgls4)
# plot(CAgls)
# plot(CAgls4) 
E1 <- resid(CAgls4)
# plot(E1 ~ Age,
#     ylab = "Ordinary residuals", data = CAall) # These are still cone shaped, 
  # need to normalize by variance
E2 <- resid(CAgls4, type = "normalized")
# plot(E2 ~ Age,
#     ylab = "Normalized residuals", data = CAall) # This is where the pattern 
  # should be gone, which it is for the most part so this worked

# Now add nestedness (site random effect)
# Package nlme and function lme allows you to include the weights call/variance 
  # structure we defined above:
CAme2<-lme(Nodes~Age*Species,random=~1|Site,weights = varExp(form =~ Age),
           data=CAall)
summary(CAme2)
# plot(CAme2)
CAme3<-lme(Nodes~Age*Species,random=~1|Site,data=CAall)
summary(CAme3)
# plot(CAme3)
summary(CAme)
Anova(CAme)
confint(CAme)
intervals(CAme2) 
CAsq<-sem.model.fits(CAme2)
# plot(CAme2, resid(., type = "normalized") ~ Age, abline = 0)

obsca<-CAall$Nodes
predca<-fitted(CAme2)
poca<-lm(predca~obsca)
casq<-summary(poca)$r.squared

# Southwest
head(SWall)
sapply(SWall, function(x) sum(is.na(x)))
# doesn't work with missing data like lmer will)
SWme<-lmer(Nodes~Age*Species+(1|Site),data=SWall)

SWgls<-gls(Nodes~Age*Species,data=SWall) # Linear model no random effects or 
  # variance structure to compare to variance structures
vf1Fixed <- varFixed(~Age) # Specify that var increases with age (in a 
  # multiplicative way)
SWgls1<-gls(Nodes~Age*Species,weights = vf1Fixed,data=SWall) 
summary(SWgls1)
# plot(SWgls1)
vf3 <- varPower(form =~ Age) # Variance varies as a power of the covariate
SWgls3<-gls(Nodes~Age*Species,weights = vf3,data=SWall) 
vf4 <- varExp(form =~ Age) # Exponential of the variance covariate
SWgls4<-gls(Nodes~Age*Species,weights = vf4,data=SWall)
AIC(SWgls,SWgls1,SWgls3,SWgls4) # Variance structure in gls3 is best
summary(SWgls4)
# plot(SWgls)
# plot(SWgls4)  
E1 <- resid(SWgls4)
# plot(E1 ~ Age,
#     ylab = "Ordinary residuals", data = SWall) # These are still cone shaped, 
  # need to normalize by variance
E2 <- resid(SWgls4, type = "normalized")
# plot(E2 ~ Age,
#     ylab = "Normalized residuals", data = SWall) # This is where the pattern 
  # should be gone, which it is for the most part so this worked

# Now add nestedness (site random effect)
# Package nlme and function lme allows you to include the weights call/variance 
  # structure we defined above:
SWme1<-lme(Nodes~Age*Species,random=~1|Site,data=SWall)
summary(SWme1)
# plot(SWme1)
SWme2<-lme(Nodes~Age*Species,random=~1|Site,weights = vf3,data=SWall)
summary(SWme2)
# plot(SWme2)
summary(SWme)
Anova(SWme)
confint(SWme)
intervals(SWme2)  
SWsq<-sem.model.fits(SWme2)
# plot(SWme1, resid(., type = "normalized") ~ Age, abline = 0)

obssw<-SWall$Nodes
predsw<-fitted(SWme2)
posw<-lm(predsw~obssw)
swsq<-summary(posw)$r.squared

# All regions
head(All)
sapply(All, function(x) sum(is.na(x)))
All<-All[is.na(All$Age)==F,] # Get rid of two rows with NA for age
All<-All[is.na(All$Nodes)==F,] # Get rid of 6 rows with NA for nodes (gls 
  # doesn't work with missing data like lmer will)
Allme<-lmer(Nodes~Age*Species+(1|Site),data=All)

Allgls<-gls(Nodes~Age*Species,data=All) # Linear model no random effects or 
  # variance structure to compare to variance structures
vf1Fixed <- varFixed(~Age) # Specify that var increases with age (in a 
  # multiplicative way)
Allgls1<-gls(Nodes~Age*Species,weights = vf1Fixed,data=All) 
summary(Allgls1)
# plot(Allgls1)
vf3 <- varPower(form =~ Age) # Variance varies as a power of the covariate
Allgls3<-gls(Nodes~Age*Species,weights = vf3,data=All) 
vf4 <- varExp(form =~ Age) # Exponential of the variance covariate
Allgls4<-gls(Nodes~Age*Species,weights = vf4,data=All)
AIC(Allgls,Allgls1,Allgls3,Allgls4) # Variance structure in gls4 is best
summary(Allgls4)
# plot(Allgls)
# plot(Allgls4) 
E1 <- resid(Allgls4)
# plot(E1 ~ Age,
#     ylab = "Ordinary residuals", data = All) # These are still cone shaped, 
  # need to normalize by variance
E2 <- resid(Allgls4, type = "normalized")
# plot(E2 ~ Age,
#     ylab = "Normalized residuals", data = All) # This is where the pattern 
  # should be gone, which it is for the most part so this worked

# Now add nestedness (site random effect)
# Package nlme and function lme allows you to include the weights call/variance 
  # structure we defined above:
Allme1<-lme(Nodes~Age*Species,random=~1|Site,data=All)
summary(Allme1)
# plot(Allme1)
Allme2<-lme(Nodes~Age*Species,random=~1|Site,weights = vf4,data=All)
summary(Allme2)
# plot(Allme2)
summary(Allme)
Anova(Allme)
confint(Allme)
intervals(Allme2) 
arsq<-sem.model.fits(Allme2)

# plot(Allme2, resid(., type = "normalized") ~ Age, abline = 0)

obs<-All$Nodes
pred<-fitted(Allme2)
po<-lm(pred~obs)
arsq<-summary(po)$r.squared


### Make Figure 2 - Nodes vs. Age scatterplot ----------------------------------

# Northern Rockies

nr_r <- paste("R^2 == ", format(round(nrsq,2),nsmall=2))
nrnewdf<-expand.grid(Species=c("PIPO","PSME"),Age=seq(0,24,1))
#nrnewdf$pred<-predict(NRme,nrnewdf,re.form=NA)
nrnewdf$pred<-predict(NRme2,nrnewdf,level=0) # Predictions from lme need level=0
  # instead of re.form=NA to use only fixed effects

NR<-ggplot(NRall, aes(Age, Nodes, color=Species, shape=Species))+
  geom_point(position=position_jitter(width=.5,height=.5))+
  scale_shape_manual(values=c(17,1)) +
  geom_abline(color="grey")+
  geom_abline(slope=1.5,intercept=0,color="darkgrey",lty=2)+
  geom_abline(slope=.5,intercept=0,color="darkgrey",lty=2)+
  geom_line(data=nrnewdf,aes(Age,pred,color=Species),lwd=1)+
  scale_color_manual(values=c("black", "grey50"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  ylab(" ")+
  xlab(" ")+
  scale_x_continuous(expand = c(0, 0),limits=c(0,25)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(0,25))+
  annotate("text", x = 4, y = 24, label = nr_r, parse=TRUE, color="black")+
  ggtitle("Northern Rockies")
NR

# Southwest

sw_r <- paste("R^2 == ", format(round(swsq,2),nsmall=2))
swnewdf<-expand.grid(Species=c("PIPO","PSME"),Age=seq(0,max(SWall$Age),1))
#swnewdf$pred<-predict(SWme,swnewdf,re.form=NA)
swnewdf$pred<-predict(SWme2,swnewdf,level=0) # Predictions from lme need level=0
  # instead of re.form=NA to use only fixed effects

SW<-ggplot(SWall, aes(Age, Nodes, color=Species, shape=Species))+
  geom_point(position=position_jitter(width=.5,height=.5))+
  scale_shape_manual(values=c(17,1)) +
  geom_abline(color="grey")+
  geom_abline(slope=1.5,intercept=0,color="darkgrey",lty=2)+
  geom_abline(slope=.5,intercept=0,color="darkgrey",lty=2)+
  geom_line(data=swnewdf,aes(Age,pred,color=Species),lwd=1)+
  scale_color_manual(values=c("black", "grey50"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  ylab("Nodes (#)")+
  xlab("Age (yr)")+
  scale_x_continuous(expand = c(0, 0),limits=c(0,25)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(0,25))+
  annotate("text", x = 4, y = 24, label = sw_r, parse=TRUE, color="black")+
  ggtitle("Southwest")
SW

# California

ca_r <- paste("R^2 == ", format(round(casq,2),nsmall=2))
canewdf<-expand.grid(Species=c("PIPO","PSME"),Age=seq(0,max(CAall$Age),1))
#canewdf$pred<-predict(CAme,canewdf,re.form=NA)
canewdf$pred<-predict(CAme2,canewdf,level=0) # Predictions from lme need level=0
  # instead of re.form=NA to use only fixed effects

CA<-ggplot(CAall, aes(Age, Nodes, color=Species, shape=Species))+
  geom_point(position=position_jitter(width=.5,height=.5))+
  scale_shape_manual(values=c(17,1)) +
  geom_abline(color="grey")+
  geom_abline(slope=1.5,intercept=0,color="darkgrey",lty=2)+
  geom_abline(slope=.5,intercept=0,color="darkgrey",lty=2)+
  geom_line(data=canewdf,aes(Age,pred,color=Species),lwd=1)+
  scale_color_manual(values=c("black", "grey50"))+
  theme_bw()+
  theme(legend.position=c(.9,.2))+
  theme(legend.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  ylab("Nodes (#)")+
  xlab(" ")+
  scale_x_continuous(expand = c(0, 0),limits=c(0,25)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(0,25))+
  annotate("text", x = 4, y = 24, label = ca_r, parse=TRUE, color="black")+
  annotate("text", x = 13, y = 23, label = "paste(\\"+50%\\")", 
           parse=TRUE, color="dark grey")+
  annotate("text", x = 23, y = 10, label = "paste(\\"-50%\\")", 
           parse=TRUE, color="dark grey")+
  ggtitle("California")
CA

# All regions

all_r <- paste("R^2 == ", format(round(arsq,2),nsmall=2))
allnewdf<-expand.grid(Species=c("PIPO","PSME"),Age=seq(0,24,1))
#allnewdf$pred<-predict(Allme,allnewdf,re.form=NA)
allnewdf$pred<-predict(Allme2,allnewdf,level=0) # Predictions from lme need 
  # level=0 instead of re.form=NA to use only fixed effects

All_nodeage<-ggplot(All, aes(Age, Nodes, color=Species, shape=Species))+
  geom_point(position=position_jitter(width=.5,height=.5))+
  scale_shape_manual(values=c(17,1)) +
  geom_abline(color="grey")+
  geom_abline(slope=1.5,intercept=0,color="darkgrey",lty=2)+
  geom_abline(slope=.5,intercept=0,color="darkgrey",lty=2)+
  geom_line(data=allnewdf,aes(Age,pred,color=Species),lwd=1)+
  scale_color_manual(values=c("black", "grey50"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  ylab(" ")+
  xlab("Age (yr)")+
  scale_x_continuous(expand = c(0, 0),limits=c(0,25)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(0,25))+
  annotate("text", x = 4, y = 24, label = all_r, parse=TRUE, color="black")+
  ggtitle("All Regions")
All_nodeage

plot_grid(CA,NR,SW,All_nodeage,ncol=2,align="v",labels = c("a", "b", "c", "d"))
# ggsave("Fig2_bw.tiff",dpi=200)


#### Calculate accuracy and bias metrics ------------------------------------

# Northern Rockies

# Percent accuracy of node counts
(NR_zero_acc<-mean(NRall$zero_yr,na.rm=T))*100
(NR_one_acc<-mean(NRall$one_yr,na.rm=T))*100
(NR_two_acc<-mean(NRall$two_yr,na.rm=T))*100
(NR_five_acc<-mean(NRall$five_yr,na.rm=T))*100

# Mean error
(NR_mean_error<-mean(NRall$AgeDif,na.rm=T))
(NR_pipo_mean_error<-mean(NRpipo$AgeDif))
(NR_psme_mean_error<-mean(NRpsme$AgeDif))

# 2017 nodes
# Percent accuracy of node counts
(NR_zero_acc17<-mean(NRall_2017$zero_yr_node))*100
(NR_one_acc17<-mean(NRall_2017$one_yr_node))*100
(NR_five_acc17<-mean(NRall_2017$five_yr_node))*100
(NR_two_acc17<-mean(NRall_2017$two_yr_node))*100

# Mean error
(NR_mean_error17<-mean(NRall_2017$AgeDif))
(NR_pipo_mean_error17<-mean(NRpipo_2017$AgeDif))
(NR_psme_mean_error17<-mean(NRpsme_2017$AgeDif))

# Northern Rockies 2017 bud scars#

# Percent accuracy of bud scar counts
(NR_zero_acc_bud<-mean(NRall_2017$zero_yr_bud))*100
(NR_one_acc_bud<-mean(NRall_2017$one_yr_bud))*100
(NR_two_acc_bud<-mean(NRall_2017$two_yr_bud))*100
(NR_five_acc_bud<-mean(NRall_2017$five_yr_bud))*100

# Mean error
(NR_mean_error_bud<-mean(NRall_2017$AgeDif_Buds))
(NR_pipo_mean_error_bud<-mean(NRpipo_2017$AgeDif_Buds))
(NR_psme_mean_error_bud<-mean(NRpsme_2017$AgeDif_Buds))

mean(NRall_2017$Buds-NRall_2017$Nodes,na.rm=T)

# NRpipo
# Percent accuracy of node counts
(NRpipo_zero_acc<-mean(NRpipo$zero_yr))*100
(NRpipo_one_acc<-mean(NRpipo$one_yr))*100
(NRpipo_five_acc<-mean(NRpipo$five_yr))*100
(NRpipo_two_acc<-mean(NRpipo$two_yr))*100
# NRpsme
# Percent accuracy of node counts
(NRpsme_zero_acc<-mean(NRpsme$zero_yr))*100
(NRpsme_one_acc<-mean(NRpsme$one_yr))*100
(NRpsme_five_acc<-mean(NRpsme$five_yr))*100
(NRpsme_two_acc<-mean(NRpsme$two_yr))*100

# Southwest#

# Percent accuracy of node counts
(SW_zero_acc<-mean(SWall$zero_yr))*100
(SW_one_acc<-mean(SWall$one_yr))*100
(SW_five_acc<-mean(SWall$five_yr))*100
(SW_two_acc<-mean(SWall$two_yr))*100

# Mean error
(SW_mean_error<-mean(SWall$AgeDif))
(SW_pipo_mean_error<-mean(SWpipo$AgeDif))
(SW_psme_mean_error<-mean(SWpsme$AgeDif))

# Swpipo
# Percent accuracy of node counts
(SWpipo_zero_acc<-mean(SWpipo$zero_yr))*100
(SWpipo_one_acc<-mean(SWpipo$one_yr))*100
(SWpipo_five_acc<-mean(SWpipo$five_yr))*100
(SWpipo_two_acc<-mean(SWpipo$two_yr))*100
# SWpsme
# Percent accuracy of node counts
(SWpsme_zero_acc<-mean(SWpsme$zero_yr))*100
(SWpsme_one_acc<-mean(SWpsme$one_yr))*100
(SWpsme_five_acc<-mean(SWpsme$five_yr))*100
(SWpsme_two_acc<-mean(SWpsme$two_yr))*100

# California#

# Percent accuracy of node counts
(CA_zero_acc<-mean(CAall$zero_yr))*100
(CA_one_acc<-mean(CAall$one_yr))*100
(CA_five_acc<-mean(CAall$five_yr))*100
(CA_two_acc<-mean(CAall$two_yr))*100

# Mean error
(CA_mean_error<-mean(CAall$AgeDif))
(CA_pipo_mean_error<-mean(CApipo$AgeDif))
(CA_psme_mean_error<-mean(CApsme$AgeDif))

# CApipo
# Percent accuracy of node counts
(CApipo_zero_acc<-mean(CApipo$zero_yr))*100
(CApipo_one_acc<-mean(CApipo$one_yr))*100
(CApipo_five_acc<-mean(CApipo$five_yr))*100
(CApipo_two_acc<-mean(CApipo$two_yr))*100
# CApsme
# Percent accuracy of node counts
(CApsme_zero_acc<-mean(CApsme$zero_yr))*100
(CApsme_one_acc<-mean(CApsme$one_yr))*100
(CApsme_five_acc<-mean(CApsme$five_yr))*100
(CApsme_two_acc<-mean(CApsme$two_yr))*100

# All together#

# Percent accuracy of node counts
(All_zero_acc<-mean(All$zero_yr,na.rm=T))*100
(All_one_acc<-mean(All$one_yr,na.rm=T))*100
(All_five_acc<-mean(All$five_yr,na.rm=T))*100
(All_two_acc<-mean(All$two_yr,na.rm=T))*100

# Mean error
(All_mean_error<-mean(All$AgeDif,na.rm=T))

# PIPO 
# Percent accuracy#
(pipo_zero_acc<-mean(Allpipo$zero_yr))*100
(pipo_one_acc<-mean(Allpipo$one_yr))*100
(pipo_five_acc<-mean(Allpipo$five_yr))*100
(pipo_two_acc<-mean(Allpipo$two_yr))*100

# Mean error
(pipo_mean_error<-mean(Allpipo$AgeDif))

# PSME
# Percent accuracy#
(psme_zero_acc<-mean(Allpsme$zero_yr))*100
(psme_one_acc<-mean(Allpsme$one_yr))*100
(psme_five_acc<-mean(Allpsme$five_yr))*100
(psme_two_acc<-mean(Allpsme$two_yr))*100

# Mean error
(psme_mean_error<-mean(Allpsme$AgeDif))


##### Bias vs. age - bias mixed effects model ---------------------------------

# Establish growth rate variable
All$growth<-All$Height_cm/All$Age

# Baseline mixed effects model without variance structure
biasmefull<-lmer(AgeDif~(Age+Species+Region+growth)^2+(1|Site),REML=F,data=All)
summary(biasmefull)
Anova(biasmefull)
sem.model.fits(biasmefull)

biasme1<-lmer(AgeDif~Age+Species+Region+growth+Age:Species+Age:Region+Age:growth+
                Species:Region+Region:growth+(1|Site),data=All)
summary(biasme1)
Anova(biasme1)
biasrsq<-sem.model.fits(biasme1)
confint(biasme1)
# plot(allEffects(biasme1))

# Add in non-constant variance structure
head(All)
sapply(All, function(x) sum(is.na(x)))
All<-All[is.na(All$Age)==F,] # Get rid of two rows with NA for age
All<-All[is.na(All$Nodes)==F,] # Get rid of 6 rows with NA for nodes (gls 
  # doesn't work with missing data like lmer will)
Allbias<-lmer(AgeDif~Age+Species+Region+growth+Age:Species+Age:Region+Age:growth+
              Species:Region+Region:growth+(1|Site),data=All)

Allbiasgls<-gls(AgeDif~Age+Species+Region+growth+Age:Species+Age:Region+
                  Age:growth+ Species:Region+Region:growth,data=All) 
# Linear model no random effects or variance structure to compare to 
  # variance structures
vf1Fixed <- varFixed(~Age) # Specify that var increases with age (in a 
  # multiplicative way)
Allbiasgls1<-gls(AgeDif~Age+Species+Region+growth+Age:Species+Age:Region+
                   Age:growth+Species:Region+Region:growth,weights = vf1Fixed,
                 data=All) 
summary(Allbiasgls1)
# plot(Allbiasgls1)
vf3 <- varPower(form =~ Age) # Variance varies as a power of the covariate
Allbiasgls3<-gls(AgeDif~Age+Species+Region+growth+Age:Species+Age:Region+
                   Age:growth+
               Species:Region+Region:growth,weights = vf3,data=All) 
vf4 <- varExp(form =~ Age) # Exponential of the variance covariate
Allbiasgls4<-gls(AgeDif~Age+Species+Region+growth+Age:Species+Age:Region+
                   Age:growth+
               Species:Region+Region:growth,weights = vf4,data=All)
AIC(Allbiasgls,Allbiasgls1,Allbiasgls3,Allbiasgls4)
summary(Allgls4)
# plot(Allgls)
# plot(Allgls4)  
E1 <- resid(Allgls4)
# plot(E1 ~ Age,
#     ylab = "Ordinary residuals", data = All) # These are still cone shaped,
  # need to normalize by variance
E2 <- resid(Allgls4, type = "normalized")
# plot(E2 ~ Age,
#     ylab = "Normalized residuals", data = All) # This is where the pattern 
  # should be gone, which it is for the most part so this worked

# Now add nestedness (site random effect)
# Package nlme and function lme allows you to include the weights call/variance 
  # structure we defined above:
Allme2<-lme(AgeDif~Age+Species+Region+growth+Age:Species+Age:Region+Age:growth+
              Species:Region+Region:growth,random=~1|Site,weights = vf4,
            data=All)
summary(Allme2)
# plot(Allme2)
summary(Allme)
Anova(Allme)
confint(Allme)
intervals(Allme2)  
biasrq<-sem.model.fits(Allme2)

obsbias<-All$AgeDif
predbias<-fitted(Allme2)
pobias<-lm(predbias~obsbias)
biasrsq<-summary(pobias)$r.squared

# Calculate predictions

# Northern Rockies
brsq <- paste("r^2 == ", format(round(biasrsq,2),nsmall=2))
biasnewdf<-expand.grid(Species=c("PIPO","PSME"),Region=c("CA","NR","SW"),
                       Age=seq(0,24,1),growth=4.25)
# biasnewdf$pred<-predict(biasme1,biasnewdf,re.form=NA)
biasnewdf$pred<-predict(Allme2,biasnewdf,level=0) # Predictions from lme need 
  # level=0 instead of re.form=NA to use only fixed effects
head(biasnewdf)
biasnewnr<-subset(biasnewdf,Region=="NR")

# California
biasnewca<-expand.grid(Species=c("PIPO","PSME"),Region=c("CA","NR","SW"),
                       Age=seq(0,max(CAall$Age),1),growth=4.25)
biasnewca$pred<-predict(Allme2,biasnewca,level=0) # Predictions from lme need 
  # level=0 instead of re.form=NA to use only fixed effects
biasnewca<-subset(biasnewca,Region=="CA")
head(biasnewca)

# Southwest
biasnewsw<-expand.grid(Species=c("PIPO","PSME"),Region=c("CA","NR","SW"),
                       Age=seq(0,max(SWall$Age),1),growth=4.25)
biasnewsw$pred<-predict(Allme2,biasnewsw,level=0) # Predictions from lme need 
  # level=0 instead of re.form=NA to use only fixed effects
biasnewsw<-subset(biasnewsw,Region=="SW")
head(biasnewsw)

# All regions together
biasmeall<-lme(AgeDif~Age+Species+growth+Age:Species+Age:growth,random=~1|Site,
               weights = vf4,data=All)
biasnewall<-expand.grid(Species=c("PIPO","PSME"),Age=seq(0,24,1),growth=4.25)
biasnewall$pred<-predict(biasmeall,biasnewall,level=0) # Predictions from lme 
  # need level=0 instead of re.form=NA to use only fixed effects
head(biasnewall)


###### Figure 3 - Bias as a function of sample age -----------------------------

# Northern Rockies 

bias_r <- paste("R^2 == ", round(biasrsq,2))

NR_bias<-ggplot(NRall, aes(Age, AgeDif, color=Species, shape=Species)) +
  scale_shape_manual(values=c(17,1)) +
  geom_point(position=position_jitter(width=.5,height=.5))+
  geom_hline(yintercept=0,color="grey")+
  geom_line(data=biasnewnr,aes(Age,pred, colour=Species),lwd=1)+
  scale_color_manual(values=c("black", "grey50"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,25)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(-5,20))+
  ylab(" ")+
  xlab(" ")+
  ggtitle("Northern Rockies")

# Southwest

SW_bias<-ggplot(SWall, aes(Age, AgeDif, color=Species, shape=Species)) +
  scale_shape_manual(values=c(17,1)) +
  geom_point(position=position_jitter(width=.5,height=.5))+
  geom_hline(yintercept=0,color="grey")+
  geom_line(data=biasnewsw,aes(Age,pred, colour=Species),lwd=1)+
  scale_color_manual(values=c("black", "grey50"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,25)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(-5,20))+
  ylab("Bias (Rings - Nodes)")+
  xlab("Age (yr)")+
  ggtitle("Southwest")

# California

CA_bias<-ggplot(CAall, aes(Age, AgeDif, color=Species, shape=Species)) +
  scale_shape_manual(values=c(17,1)) +
  geom_point(position=position_jitter(width=.5,height=.5))+
  geom_hline(yintercept=0,color="grey")+
  geom_line(data=biasnewca,aes(Age,pred, colour=Species),lwd=1)+
  scale_color_manual(values=c("black", "grey50"))+
  theme_bw()+
  theme(legend.position=c(.8,.9))+
  theme(legend.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text.y = element_text(size=10))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,25)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(-5,20))+
  ylab("Bias (Rings - Nodes)")+
  xlab(" ")+
  ggtitle("California")

# All Regions

all_bias<-ggplot(All, aes(Age, AgeDif, color=Species, shape=Species)) +
  scale_shape_manual(values=c(17,1)) +
  geom_point(position=position_jitter(width=.5,height=.5))+
  geom_hline(yintercept=0,color="grey")+
  geom_line(data=biasnewall,aes(Age,pred, colour=Species),lwd=1)+
  scale_color_manual(values=c("black", "grey50"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,25)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(-5,20))+
  ylab(" ")+
  xlab("Age (yr)")+
  ggtitle("All Regions")

plot_grid(CA_bias,NR_bias,SW_bias,all_bias,ncol=2,align="v",
          labels = c("a", "b", "c", "d"))
# ggsave("Fig3_bw.tiff",dpi=200)


####### Figure S1 - Bias as a function of growth rate --------------------------

# I(Height_cm/Age)
# Prediction dataframes for plotting
biasnewdf<-expand.grid(Species=c("PIPO","PSME"),Region=c("CA","NR","SW"),
                       Age=c(5,10,15),growth=seq(0,41))
biasnewdf$pred<-predict(Allme2,biasnewdf,level=0)
head(biasnewdf)

# Northern rockies
biasnewnr<-expand.grid(Species=c("PIPO","PSME"),Region=c("CA","NR","SW"),
                       Age=c(5,10,15),growth=seq(0,38))
biasnewnr$pred<-predict(Allme2,biasnewnr,level=0)
biasnewnr<-subset(biasnewnr,Region=="NR")

# California
biasnewca<-expand.grid(Species=c("PIPO","PSME"),Region=c("CA","NR","SW"),
                       Age=c(5,10,15),growth=seq(0,41))
biasnewca$pred<-predict(Allme2,biasnewca,level=0)
biasnewca<-subset(biasnewca,Region=="CA")

# Southwest
biasnewsw<-expand.grid(Species=c("PIPO","PSME"),Region=c("CA","NR","SW"),
                       Age=c(5,10,15),growth=seq(0,30))
biasnewsw$pred<-predict(Allme2,biasnewsw,level=0)
biasnewsw<-subset(biasnewsw,Region=="SW")

# All regions together

biasmeall<-lme(AgeDif~Age+Species++growth+Age:Species+Age:growth,random=~1|Site,
               weights = vf4,data=All)
biasnewall<-expand.grid(Species=c("PIPO","PSME"),Age=c(5,10,15),
                        growth=seq(0,41))
biasnewall$pred<-predict(biasmeall,biasnewall,level=0) # Predictions from lme 
  # need level=0 instead of re.form=NA to use only fixed effects
head(biasnewall)

# Northern Rockies

bias_r <- paste("R^2 == ", round(biasrsq,2))

NR_bias_gr<-ggplot(NRall, aes(I(Height_cm/Age), AgeDif, color=Species, 
                              shape=Species)) +
  scale_shape_manual(values=c(1,2)) +
  geom_point(position=position_jitter(width=.5,height=.5))+
  geom_hline(yintercept=0,color="grey")+
  geom_line(data=subset(biasnewnr,Age==5),aes(growth,pred, colour=Species),
            lwd=1)+
  geom_line(data=subset(biasnewnr,Age==10),aes(growth,pred, colour=Species),
            lty=2,lwd=1)+
  geom_line(data=subset(biasnewnr,Age==15),aes(growth,pred, colour=Species),
            lty=3,lwd=1)+
  scale_color_manual(values=c("blue", "black"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,45)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(-5,20))+
  ylab(" ")+
  xlab(" ")+
  ggtitle("Northern Rockies")

# Southwest

SW_bias_gr<-ggplot(SWall, aes(I(Height_cm/Age), AgeDif, color=Species, 
                              shape=Species)) +
  scale_shape_manual(values=c(1,2)) +
  geom_point(position=position_jitter(width=.5,height=.5))+
  geom_hline(yintercept=0,color="grey")+
  geom_line(data=subset(biasnewsw,Age==5),aes(growth,pred, colour=Species),
            lwd=1)+
  geom_line(data=subset(biasnewsw,Age==10),aes(growth,pred, colour=Species),
            lty=2,lwd=1)+
  geom_line(data=subset(biasnewsw,Age==15),aes(growth,pred, colour=Species),
            lty=3,lwd=1)+
  scale_color_manual(values=c("blue", "black"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,45)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(-5,20))+
  ylab("Bias (Rings - Nodes)")+
  xlab("Growth rate (cm/yr)")+
  ggtitle("Southwest")

# California

CA_bias_gr<-ggplot(CAall, aes(I(Height_cm/Age), AgeDif, color=Species, 
                              shape=Species)) +
  scale_shape_manual(values=c(1,2)) +
  geom_point(position=position_jitter(width=.5,height=.5))+
  geom_hline(yintercept=0,color="grey")+
  geom_line(data=subset(biasnewca,Age==5),aes(growth,pred, colour=Species),
            lwd=1)+
  geom_line(data=subset(biasnewca,Age==10),aes(growth,pred, colour=Species),
            lty=2,lwd=1)+
  geom_line(data=subset(biasnewca,Age==15),aes(growth,pred, colour=Species),
            lty=3,lwd=1)+
  scale_color_manual(values=c("blue", "black"))+
  theme_bw()+
  theme(legend.position=c(.8,.9))+
  theme(legend.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text.y = element_text(size=10))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,45)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(-5,20))+
  ylab("Bias (Rings - Nodes)")+
  xlab(" ")+
  ggtitle("California")

# All Regions

all_bias_gr<-ggplot(All, aes(I(Height_cm/Age), AgeDif, color=Species, 
                             shape=Species)) +
  scale_shape_manual(values=c(1,2)) +
  geom_point(position=position_jitter(width=.5,height=.5))+
  geom_hline(yintercept=0,color="grey")+
  geom_line(data=subset(biasnewall,Age==5),aes(growth,pred, colour=Species),
            lwd=1)+
  geom_line(data=subset(biasnewall,Age==10),aes(growth,pred, colour=Species),
            lty=2,lwd=1)+
  geom_line(data=subset(biasnewall,Age==15),aes(growth,pred, colour=Species),
            lty=3,lwd=1)+
  scale_color_manual(values=c("blue", "black"))+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,45)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(-5,20))+
  ylab(" ")+
  xlab("Growth rate (cm/yr)")+
  ggtitle("All Regions")

plot_grid(CA_bias_gr,NR_bias_gr,SW_bias_gr,all_bias_gr,ncol=2,align="v",
          labels = c("a", "b", "c", "d"))
# ggsave("FigS1.tiff",dpi=200)


######## Test differences in bias -------------------------------------------

# Mean error tests
t.test(Allpipo$AgeDif,Allpsme$AgeDif,conf.level=0.95)

a<-aov(AgeDif~Region,data=All)
summary(a)
TukeyHSD(a)

# Accuracy tests
region_0<-aov(zero_yr~Region,data=All)
summary(region_0)
region_1<-aov(one_yr~Region,data=All)
summary(region_1)
region_2<-aov(two_yr~Region,data=All)
summary(region_2)
region_5<-aov(five_yr~Region,data=All)
summary(region_5)

pzero<-t.test(Allpipo$zero_yr,Allpsme$zero_yr,conf.level=0.95)[3]
pone<-t.test(Allpipo$one_yr,Allpsme$one_yr,conf.level=0.95)[3]
ptwo<-t.test(Allpipo$two_yr,Allpsme$two_yr,conf.level=0.95)[3]
pfive<-t.test(Allpipo$five_yr,Allpsme$five_yr,conf.level=0.95)[3]

(pval<-c(pzero,pone,ptwo,pfive))
(pval<-as.data.frame(pval))
p.adjust(pval, method = "bonferroni")


######### Predict age from nodes --------------------------------------------

# Fit linear mixed effects models without variance structure as baseline
  # for each region
predNRme<-lmer(Age~Nodes+Species+(1|Site),data=NRall)
summary(predNRme)
Anova(predNRme)
confint(predNRme)
prednrsq<-sem.model.fits(predNRme)

predSWme<-lmer(Age~Nodes*Species+(1|Site),data=SWall)
summary(predSWme)
Anova(predSWme)
confint(predSWme)
predsrsq<-sem.model.fits(predSWme)

predCAme<-lmer(Age~Nodes+Species+(1|Site),data=CAall)
summary(predCAme)
Anova(predCAme)
confint(predCAme)
predcrsq<-sem.model.fits(predCAme)

predAllme<-lmer(Age~Nodes+Species+(1|Site),data=All)
summary(predAllme)
Anova(predAllme)
confint(predAllme)
predarsq<-sem.model.fits(predAllme)

# qqnorm(residuals(predNRme))
# qqnorm(residuals(predCAme))
# qqnorm(residuals(predSWme))
# qqnorm(residuals(predAllme))

# Non-constant variance, need variance structure
# Add in exponential/power variance structure using gls (Zuur et al. 2009)

# All regions
head(All)
sapply(All, function(x) sum(is.na(x)))
All<-All[is.na(All$Age)==F,] # Get rid of two rows with NA for age
All<-All[is.na(All$Nodes)==F,] # Get rid of 6 rows with NA for nodes (gls 
# doesn't work with missing data like lmer will)
Allme<-lmer(Age~Nodes*Species+(1|Site),data=All)

Allgls<-gls(Age~Nodes*Species,data=All) # Linear model no random effects or 
  # variance structure to compare to variance structures
# plot(Allgls)
vf1Fixed <- varFixed(~Nodes) # Specify that var increases with age (in a 
  # multiplicative way)
Allgls1<-gls(Age~Nodes*Species,weights = vf1Fixed,data=All) 
summary(Allgls1)
vf3 <- varPower(form =~ Nodes) # Variance varies as a power of the covariate
Allgls3<-gls(Age~Nodes*Species,weights = vf3,data=All) 
vf4 <- varExp(form =~ Nodes) # Exponential of the variance covariate
Allgls4<-gls(Age~Nodes*Species,weights = vf4,data=All)
AIC(Allgls,Allgls1,Allgls3,Allgls4) # Variance structure in gls4 is best
summary(Allgls4)
sem.model.fits(Allgls4)
# Now add nestedness (site random effect)
# Package nlme and function lme allows you to include the weights call/variance 
  # structure we defined above:
predAllme2<-lme(Age~Nodes*Species,random=~1|Site,data=All)
summary(predAllme2)
# plot(predAllme2)
predAllme<-lme(Age~Nodes*Species,random=~1|Site,weights = vf4,data=All)
summary(predAllme)
# plot(predAllme)
intervals(predAllme2) 
predarsq<-sem.model.fits(predAllme)
# plot(predAllme, resid(., type = "normalized") ~ Age, abline = 0)

obs<-All$Age
pred<-fitted(predAllme)
po<-lm(pred~obs)
arsq<-summary(po)$r.squared


# California
head(CAall)
sapply(CAall, function(x) sum(is.na(x)))
CAall<-CAall[is.na(CAall$Age)==F,] # Get rid of two rows with NA for age
CAall<-CAall[is.na(CAall$Nodes)==F,] # Get rid of 6 rows with NA for nodes (gls 
  # doesn't work with missing data like lmer will)
CAallme<-lmer(Age~Nodes*Species+(1|Site),data=CAall)
summary(CAallme)
sem.model.fits(CAallme)

CAallgls<-gls(Age~Nodes*Species,data=CAall) # Linear model no random effects or 
  # variance structure to compare to variance structures
# plot(CAallgls)
vf1Fixed <- varFixed(~Nodes) # Specify that var increases with age (in a 
  # multiplicative way)
CAallgls1<-gls(Age~Nodes*Species,weights = vf1Fixed,data=CAall) 
summary(CAallgls1)
vf3 <- varPower(form =~ Nodes) # Variance varies as a power of the covariate
CAallgls3<-gls(Age~Nodes*Species,weights = vf3,data=CAall) 
vf4 <- varExp(form =~ Nodes) # Exponential of the variance covariate
CAallgls4<-gls(Age~Nodes*Species,weights = vf4,data=CAall)
AIC(CAallgls,CAallgls1,CAallgls3,CAallgls4) # Variance structure in gls4 is best
summary(CAallgls4)
predcrq<-sem.model.fits(CAallgls4)

# Now add nestedness (site random effect)
# Package nlme and function lme CAallows you to include the weights 
# CAall/variance structure we defined above:
predCAallme<-lme(Age~Nodes*Species,random=~1|Site,weights = vf4,data=CAall)
summary(predCAallme)
# plot(predCAallme)
intervals(predCAallme) #not too different 
predcrsq<-sem.model.fits(predCAallme)

predCAallme2<-lme(Age~Nodes*Species,random=~1|Site,data=CAall)
summary(predCAallme2)
# plot(predCAallme2)
intervals(predCAallme2)  
predcrsq<-sem.model.fits(predCAallme2)

# plot(CAall$Nodes,CAall$Age)
# plot(predCAallme, resid(., type = "normalized") ~ Nodes, abline = 0)

obsca<-CAall$Age
predca<-fitted(predCAallme)
poca<-lm(predca~obsca)
casq<-summary(poca)$r.squared

# Northern rockies
head(NRall)
sapply(NRall, function(x) sum(is.na(x)))
NRall<-NRall[is.na(NRall$Age)==F,] # Get rid of two rows with NA for age
NRall<-NRall[is.na(NRall$Nodes)==F,] # Get rid of 6 rows with NA for nodes (gls 
  # doesn't work with missing data like lmer will)
NRallme<-lmer(Age~Nodes*Species+(1|Site),data=NRall)
summary(NRallme)
sem.model.fits(NRallme)

NRallgls<-gls(Age~Nodes*Species,data=NRall) # Linear model no random effects or 
  # variance structure to compare to variance structures
# plot(NRallgls)
vf1Fixed <- varFixed(~Nodes) # Specify that var increases with age (in a 
  # multiplicative way )
NRallgls1<-gls(Age~Nodes*Species,weights = vf1Fixed,data=NRall) 
summary(NRallgls1)
vf3 <- varPower(form =~ Nodes) # Variance varies as a power of the covariate
NRallgls3<-gls(Age~Nodes*Species,weights = vf3,data=NRall) 
vf4 <- varExp(form =~ Nodes) # Exponential of the variance covariate
NRallgls4<-gls(Age~Nodes*Species,weights = vf4,data=NRall)
AIC(NRallgls,NRallgls1,NRallgls3,NRallgls4) # Variance structure in gls4 is best
summary(NRallgls4)
sem.model.fits(NRallgls4)

# Now add nestedness (site random effect)
# Package nlme and function lme NRallows you to include the weights 
# NRall/variance structure we defined above:
predNRallme1<-lme(Age~Nodes*Species,random=~1|Site,data=NRall)
# plot(predNRallme1)
predNRallme<-lme(Age~Nodes*Species,random=~1|Site,weights = vf4,data=NRall)
# plot(predNRallme)
summary(predNRallme)
intervals(predNRallme) #not too different 
prednrsq<-sem.model.fits(predNRallme)
# plot(predNRallme1, resid(., type = "normalized") ~ Nodes, abline = 0)

obsnr<-NRall$Age
prednr<-fitted(predNRallme)
ponr<-lm(prednr~obsnr)
nrsq<-summary(ponr)$r.squared

# Southwest
head(SWall)
sapply(SWall, function(x) sum(is.na(x)))
SWall<-SWall[is.na(SWall$Age)==F,] # Get rid of two rows with NA for age
SWall<-SWall[is.na(SWall$Nodes)==F,] # Get rid of 6 rows with NA for nodes (gls 
  # doesn't work with missing data like lmer will)
SWallme<-lmer(Age~Nodes*Species+(1|Site),data=SWall)
summary(SWallme)
sem.model.fits(SWallme)

SWallgls<-gls(Age~Nodes*Species,data=SWall) # Linear model no random effects or 
    # variance structure to compare to variance structures
# plot(SWallgls)
vf1Fixed <- varFixed(~Nodes) # Specify that var increases with age (in a 
  # multiplicative way)
SWallgls1<-gls(Age~Nodes*Species,weights = vf1Fixed,data=SWall) 
summary(SWallgls1)
vf3 <- varPower(form =~ Nodes) # Variance varies as a power of the covariate
SWallgls3<-gls(Age~Nodes*Species,weights = vf3,data=SWall) 
vf4 <- varExp(form =~ Nodes) # Exponential of the variance covariate
SWallgls4<-gls(Age~Nodes*Species,weights = vf4,data=SWall)
AIC(SWallgls,SWallgls1,SWallgls3,SWallgls4) # Variance structure in gls4 is best
summary(SWallgls4)
predsrsq<-sem.model.fits(SWallgls4)

# Now add nestedness (site random effect)
# Package nlme and function lme SWallows you to include the weights 
# SWall/variance structure we defined above:
predSWallme<-lme(Age~Nodes*Species,random=~1|Site,weights = vf4,data=SWall)
# plot(predSWallme)
summary(predSWallme)
intervals(predSWallme) 
# predsrsq<-sem.model.fits(predSWallme)

predSWallme2<-lme(Age~Nodes*Species,random=~1|Site,data=SWall)
# plot(predSWallme2)
summary(predSWallme2)
intervals(predSWallme2) 
# predsrsq<-sem.model.fits(predSWallme2)
# plot(predSWallme2, resid(., type = "normalized") ~ Nodes, abline = 0)

obssw<-SWall$Age
predsw<-fitted(predSWallme)
posw<-lm(predsw~obssw)
swsq<-summary(posw)$r.squared

# Calculate prediction intervals

# newdata = data.frame(x=c(1:20))
# predict(lm, newdata, interval="predict") 

nr_r_pi<-paste("R^2 == ", format(round(nrsq,2),nsmall=2))

ca_r_pi<-paste("R^2 == ", format(round(casq,2),nsmall=2))

sw_r_pi<-paste("R^2 == ", format(round(swsq,2),nsmall=2))

all_r_pi<-paste("R^2 == ", format(round(arsq,2),nsmall=2))

summary(predNRallme)

# Northern Rockies
nrdf2<-expand.grid(Species=c("PIPO","PSME"),Nodes=seq(0,22,1),
                   Site=unique(NRall$Site))
nrdf2$pred <- predict(predNRallme, nrdf2,level=0)
head(nrdf2)

# Calculate standard error for mixed models
# Need to account for non-constant variance structure
# [-2] drops response from formula
Designmat <- model.matrix(formula(predNRallme)[-2], nrdf2)
predvar <- diag(Designmat %*% vcov(predNRallme) %*% t(Designmat)) 
nrdf2$SE <- sqrt(predvar) 
nrdf2$SE2 <- sqrt(predvar+predNRallme$sigma^2)

cmult <- 1.96 

summary(predNRallme)

# Add variance structure to prediction intervals
nrdf2$SE.b <- sqrt(predvar*exp(-0.0617747*2*nrdf2$Nodes)) # Parameter estimate
nrdf2$SE2.b <- sqrt((predvar + 
                       predNRallme$sigma^2)*exp(-0.0617747*2*nrdf2$Nodes)) 
# Parameter estimate with variance structure
nrdf2$lwr<-nrdf2$pred-cmult*nrdf2$SE2.b
nrdf2$upr<-nrdf2$pred+cmult*nrdf2$SE2.b

head(nrdf2)

nrdf3<-ddply(nrdf2,c("Species","Nodes"), summarize,pred=mean(pred),
             upr=mean(upr),lwr=mean(lwr),SE=mean(SE2.b))
head(nrdf3)

# Test to see if prediction intervals increase with increasing nodes
nrdf3[nrdf3$Nodes==1,]$upr-nrdf3[nrdf3$Nodes==1,]$lwr
nrdf3[nrdf3$Nodes==20,]$upr-nrdf3[nrdf3$Nodes==20,]$lwr

# write.csv(nrdf3,"nr_pi.csv")

# Northern Rockies PI plot

(NR_pi<-ggplot(data=NRall, aes(x = Nodes,color=Species,shape=Species)) +
  geom_point(aes(y = Age),position=position_jitter(width=.5,height=.5)) +
  geom_ribbon(data=nrdf3,aes(x=Nodes,ymin = lwr, ymax = upr, 
                             fill = Species,col=Species), alpha = 0.2)+
  theme_bw()+
  scale_color_manual(values=c("blue", "black"))+
  scale_fill_manual(values=c("blue", "black"))+
  scale_shape_manual(values=c(1,2)) +
  geom_line(data=nrdf3,aes(Nodes,pred,col=Species))+
  theme(legend.position="none")+
  #theme(legend.position=c(.9,.2))+
  #theme(legend.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,23)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(-.5,30))+
  ylab("")+
  xlab("")+
  annotate("text", x = 17, y = 5, label = nr_r_pi, parse=TRUE, color="black")+
  ggtitle("Northern Rockies"))
NR_pi

# California

cadf2<-expand.grid(Species=c("PIPO","PSME"),Nodes=seq(0,16,1),
                   Site=unique(CAall$Site))
cadf2$pred <- predict(predCAallme, cadf2, level=0)
# Calculate standard error for mixed models
# [-2] drops response from formula
Designmat <- model.matrix(formula(predCAallme)[-2], cadf2)
predvar <- diag(Designmat %*% vcov(predCAallme) %*% t(Designmat)) 
cadf2$SE <- sqrt(predvar) 
cadf2$SE2 <- sqrt(predvar+predCAallme$sigma^2)

cmult <- 1.96  

summary(predCAallme)

# Add variance structure to prediction intervals
cadf2$SE.b <- sqrt(predvar)*exp(-0.1137666*2*cadf2$Nodes) # Parameter estimate
cadf2$SE2.b <- sqrt((predvar + 
                       predCAallme$sigma^2)*exp(-0.1137666*2*cadf2$Nodes))
cadf2$lwr<-cadf2$pred-cmult*cadf2$SE2.b
cadf2$upr<-cadf2$pred+cmult*cadf2$SE2.b

head(cadf2)

cadf3<-ddply(cadf2,c("Species","Nodes"), summarize,pred=mean(pred),upr=mean(upr)
             ,lwr=mean(lwr),SE=mean(SE2.b))
head(cadf3)

# Test to see if prediction intervals increase with increasing nodes
cadf3[cadf3$Nodes==1,]$upr-cadf3[cadf3$Nodes==1,]$lwr
cadf3[cadf3$Nodes==10,]$upr-cadf3[cadf3$Nodes==10,]$lwr

# write.csv(cadf3,"ca_pi.csv")

# California PI plot 

(CA_pi<-ggplot(data=CAall, aes(x = Nodes,color=Species,shape=Species)) +
  geom_point(aes(y = Age),position=position_jitter(width=.5,height=.5)) +
  geom_ribbon(data=cadf3,aes(x=Nodes,ymin = lwr, ymax = upr, 
                             fill = Species,col=Species), alpha = 0.2)+
  theme_bw()+
  scale_color_manual(values=c("blue", "black"))+
  scale_fill_manual(values=c("blue", "black"))+
  scale_shape_manual(values=c(1,2)) +
  geom_line(data=cadf3,aes(Nodes,pred,col=Species))+
  theme(legend.position=c(.9,.2))+
  theme(legend.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,23)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(0,30))+
  ylab("Age (yr)")+
  xlab("")+
  annotate("text", x = 15, y = 5, label = ca_r_pi, parse=TRUE, color="black")+
  ggtitle("California"))
CA_pi


# Southwest
swdf2<-expand.grid(Species=c("PIPO","PSME"),Nodes=seq(0,13,1),
                   Site=unique(SWall$Site))
swdf2$pred <- predict(predSWallme, swdf2, level=0)
# Calculate standard error for mixed models
# [-2] drops response from formula
Designmat <- model.matrix(formula(predSWallme)[-2], swdf2)
predvar <- diag(Designmat %*% vcov(predSWallme) %*% t(Designmat)) 
swdf2$SE <- sqrt(predvar) 
swdf2$SE2 <- sqrt(predvar+predSWallme$sigma^2)

cmult <- 1.96 

summary(predSWallme)

# Add variance structure to prediction intervals
swdf2$SE.b <- sqrt(predvar)*exp(-0.1637558*2*swdf2$Nodes) # Parameter estimate
swdf2$SE2.b <- sqrt((predvar + 
                       predSWallme$sigma^2)*exp(-0.1637558*2*swdf2$Nodes))
swdf2$lwr<-swdf2$pred-cmult*swdf2$SE2.b
swdf2$upr<-swdf2$pred+cmult*swdf2$SE2.b

head(swdf2)

swdf3<-ddply(swdf2,c("Species","Nodes"), summarize,pred=mean(pred),
             upr=mean(upr),lwr=mean(lwr),SE=mean(SE2.b))
head(swdf3)

# Test to see if prediction intervals increase with increasing nodes
swdf3[swdf3$Nodes==1,]$upr-swdf3[swdf3$Nodes==1,]$lwr
swdf3[swdf3$Nodes==10,]$upr-swdf3[swdf3$Nodes==10,]$lwr

# write.csv(swdf3,"sw_pi.csv")

# SW PI plot

(SW_pi<-ggplot(data=SWall, aes(x = Nodes,color=Species,shape=Species)) +
  geom_point(aes(y = Age),position=position_jitter(width=.5,height=.5)) +
  geom_ribbon(data=swdf3,aes(x=Nodes,ymin = lwr, ymax = upr, 
                             fill = Species,col=Species), alpha = 0.2)+
  theme_bw()+
  scale_color_manual(values=c("blue", "black"))+
  scale_fill_manual(values=c("blue", "black"))+
  scale_shape_manual(values=c(1,2)) +
  geom_line(data=swdf3,aes(Nodes,pred,col=Species))+
  theme(legend.position="none")+
  #theme(legend.position=c(.9,.2))+
  #theme(legend.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,23)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(-3,30))+
  ylab("Age (yr)")+
  xlab("Nodes (#)")+
  annotate("text", x = 17, y = 5, label = sw_r_pi, parse=TRUE, color="black")+
  ggtitle("Southwest"))
SW_pi


# All

adf2<-expand.grid(Species=c("PIPO","PSME"),Nodes=seq(0,22,1),
                  Site=unique(All$Site))
adf2$pred <- predict(predAllme, adf2, level=0)
# Calculate standard error for mixed models
# [-2] drops response from formula
Designmat <- model.matrix(formula(predAllme)[-2], adf2)
predvar <- diag(Designmat %*% vcov(predAllme) %*% t(Designmat)) 
adf2$SE <- sqrt(predvar) 
adf2$SE2 <- sqrt(predvar+predAllme$sigma^2)

cmult <- 1.96  

summary(predAllme)

# Add variance structure to prediction intervals
adf2$SE.b <- sqrt(predvar)*exp(-0.04369598*2*adf2$Nodes) # Parameter estimate
adf2$SE2.b <- sqrt((predvar + predAllme$sigma^2)*exp(-0.04369598*2*adf2$Nodes))
adf2$lwr<-adf2$pred-cmult*adf2$SE2.b
adf2$upr<-adf2$pred+cmult*adf2$SE2.b

head(adf2)

adf3<-ddply(adf2,c("Species","Nodes"), summarize,pred=mean(pred),upr=mean(upr),
            lwr=mean(lwr),SE=mean(SE2.b))
head(adf3)

# Test to see if prediction intervals increase with increasing nodes
adf3[adf3$Nodes==1,]$upr-adf3[adf3$Nodes==1,]$lwr
adf3[adf3$Nodes==10,]$upr-adf3[adf3$Nodes==10,]$lwr

# write.csv(adf3,"a_pi.csv")

ap<-subset(adf3,Species=="PIPO")
aps<-subset(adf3,Species=="PSME")

mean(ap$upr-ap$lwr)
mean(aps$upr-aps$lwr)

# All regions PI plot

(All_pi<-ggplot(data=All, aes(x = Nodes,color=Species,shape=Species)) +
  geom_point(aes(y = Age),position=position_jitter(width=.5,height=.5)) +
  geom_ribbon(data=adf3,aes(x=Nodes,ymin = lwr, ymax = upr, 
                            fill = Species,col=Species), alpha = 0.2)+
  theme_bw()+
  scale_color_manual(values=c("blue", "black"))+
  scale_fill_manual(values=c("blue", "black"))+
  scale_shape_manual(values=c(1,2)) +
  geom_line(data=adf3,aes(Nodes,pred,col=Species))+
  theme(legend.position="none")+
  #theme(legend.position=c(.9,.2))+
  #theme(legend.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,23)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(0,40))+
  ylab("")+
  xlab("Nodes (#)")+
  annotate("text", x = 17, y = 5, label = all_r_pi, parse=TRUE, color="black")+
  ggtitle("All Regions"))
All_pi

plot_grid(CA_pi,NR_pi,SW_pi,All_pi,ncol=2,align="v",
          labels = c("a", "b", "c", "d"))
# ggsave("Fig_pred_gls.tiff",dpi=200)


########## Prediction tool -----------------------------------------------------
# Also available independently in separate script provided

# Function for tree age predictions and 95% prediction intervals from node counts
# Regions = "All", "CA", "NR", or "SW"
# Species = "PIPO" (ponderosa pine); "PSME" (Douglas-fir)

tree_age <- function(species,node,region){
  y <- data.frame()
  node <- node[node<23]
  for(i in node){
    dat <- read.csv("Data_AgeTool.csv")
    dat <- dat[dat$Species==species,]
    dat <- dat[dat$Region==region,]
    dat <- dat[dat$Nodes==i,]
    Nodes <- i
    EstTreeAge <- dat$pred
    UprPredInt <- dat$upr
    LwrPredInt <- dat$lwr
    age <- data.frame(Nodes,EstTreeAge,UprPredInt,LwrPredInt)
    y <- rbind(y, age)
  }
  return(y)
}

# Example
mynode <- c(1,5,9,12)
tree_age("PIPO",mynode,"All") 

# Test tool with holdout dataset
# 80% of the sample size
smp_size <- floor(0.80 * nrow(All))

# Set the seed to make  partition reproducible
set.seed(150)
train_ind <- sample(seq_len(nrow(All)), size = smp_size)

train <- All[train_ind, ]
test <- All[-train_ind, ]
train <-train[!is.na(train$Age),]
train <-train[!is.na(train$Nodes),]
train <-train[!is.na(train$Species),]

test <-test[!is.na(test$Age),]
test <-test[!is.na(test$Nodes),]
test <-test[!is.na(test$Species),]

head(test)
trainlm<-lme(Age~Nodes*Species,random=~1|Site,weights = vf4,data=train)


summary(trainlm)
sem.model.fits(trainlm)

test$pred<-predict(trainlm, test, level=0)
test$pred<-round(test$pred, 0)
head(test)
# test_node<-predict(trainlm,test,level=0)
# test_node<-round(test_node,0)
# head(test_node)

# test$predicted<-test_node

# Calculate accuracy metrics for test dataset

test$AgeDif_pred<-test$Age-test$pred

head(test)

# Actual accuracy for test group
(test_zero_acc<-mean(test$zero_yr,na.rm=T))*100
(test_one_acc<-mean(test$one_yr,na.rm=T))*100
(test_five_acc<-mean(test$five_yr,na.rm=T))*100
(test_two_acc<-mean(test$two_yr,na.rm=T))*100

# Accuracy for test group using predictions
test$zero_yr_pred<-ifelse(test$AgeDif_pred==0,1,0)
test$one_yr_pred<-ifelse(abs(test$AgeDif_pred)<=1,1,0)
test$five_yr_pred<-ifelse(abs(test$AgeDif_pred)<=5,1,0)
test$two_yr_pred<-ifelse(abs(test$AgeDif_pred)<=2,1,0)

(test_zero_acc<-mean(test$zero_yr_pred,na.rm=T))*100
(test_one_acc<-mean(test$one_yr_pred,na.rm=T))*100
(test_five_acc<-mean(test$five_yr_pred,na.rm=T))*100
(test_two_acc<-mean(test$two_yr_pred,na.rm=T))*100

mean(test$AgeDif_pred,na.rm=T)

head(test)
test$predest<-test$Sample.year-test$pred+1


########### Figure S2 - prediction tool example --------------------------------
# Observed ring-based ages
obs<-ggplot(data=test)+
  geom_histogram(aes(Age),col="black",fill="grey",binwidth=1)+
  ylab("Count")+
  xlab("Rings (#)")+
  scale_y_continuous(limits=c(0, 130))+
  scale_x_continuous(breaks=seq(0,25,5),limits=c(0,25))+
  theme_bw()+
  theme(panel.grid.minor=element_blank())+
  ggtitle("")

# Predicted ages
pred<-ggplot(data=test)+
  geom_histogram(aes(pred),col="black",fill="grey",binwidth=1)+
  ylab("")+
  xlab("Predicted Tree Age (yr)")+
  scale_x_continuous(breaks=seq(0,25,5),limits=c(0,25))+
  theme_bw()+scale_y_continuous(limits=c(0,130))+
  theme(panel.grid.minor=element_blank())+
  ggtitle("")

# Observed field-based node counts
obsnode<-ggplot(data=test)+
  geom_histogram(aes(Nodes),col="black",fill="grey",binwidth=1)+
  ylab("")+
  xlab("Nodes (#)")+
  scale_y_continuous(limits=c(0, 130))+
  scale_x_continuous(breaks=seq(0,25,5),limits=c(0,25))+
  theme_bw()+
  theme(panel.grid.minor=element_blank())+
  ggtitle("")

plot_grid(obsnode,obs,pred,ncol=1)
# ggsave("Correction_Tool.tiff", dpi=200)
