
library(ape)
library(MuMIn)
library(caper)
library(dplyr)
library(lme4)
library(phytools)
library(merTools)
library(car)



setwd("")

rawmad<-read.csv(file = "MAD_Output_All_Filter_1.csv", header =TRUE)

rawmad <- rawmad %>% 
  filter(Common_Name != "Chestnut-sided Warbler")

#Function that extracts slopes for species specific LMM 
coefficients<- function (commonname){
        rawmad<-rawmad[rawmad$Common_Name == commonname ,]
        regression<-lmer(MAD~Year + (1|GridID), data = rawmad)
        MAD.Shift<-data.frame(fixef(regression))
        RMSE<-RMSE.merMod(regression)
        Common.name<-unique(rawmad$Common_Name)
        data.frame(Common.name,MAD.Shift[2,],RMSE)
}

#Create dataframe of slopes == MAD Shifts
  x<-unique(rawmad$Common_Name)
coef.df<-lapply(x, coefficients)   
  coef2.df<-data.frame(do.call(rbind, coef.df))
  
#Read in morphology parameters  
para<-read.csv(file="MorphologyParameters.csv", header= TRUE)
    para<-merge(coef2.df,para, by="Common.name", all=TRUE)  
    para<-na.omit(para)
 
#Read phylogenetic tree for our species list
    
    
tree<-read.nexus(file="50MRCtree.nex")

#check tree
plot(tree, cex=0.7)


#PRMA to generate residuals (mass adj)
    q<-setNames(log(para$WLI), para$Tree_name)
    g<-setNames((para$Body.mass..log.), para$Tree_name)
    correction2<-phyl.RMA(q, g, tree, method="lambda")
    Mass.adj.WLI<-data.frame(correction2$resid)

    #merge mass adjusted WLI to original parameter dataframe    
      row.names(para)<-para$Tree_name
        para<-merge(para, Mass.adj.WLI, by="row.names")
          para<-para[-c(1)]
            names(para)[11]<-"Mass.adj.WLI"
          
 #create object which contains comparative data (morphologies and mad shift) with tree
 #warning will show if any tips on the tree are dropped due to missing data (this should happen, no data for Chestnut Sided warbler)    
    
            
cdata<-comparative.data(tree, para, names.col = Tree_name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)   


{mod1<- pgls(MAD.Shift.2...~ log(WAR) + log(Migration.Distance.Centroid), data = cdata, lambda="ML")
mod2<- pgls(MAD.Shift.2... ~ log(HWI) + log(Migration.Distance.Centroid), data = cdata, lambda="ML")
mod3<- pgls(MAD.Shift.2... ~ Mass.adj.WLI + log(Migration.Distance.Centroid), data = cdata, lambda="ML")
mod4<- pgls(MAD.Shift.2... ~ Body.mass..log. +log(Migration.Distance.Centroid), data = cdata, lambda="ML")
mod5<- pgls(MAD.Shift.2... ~ log(WAR), data = cdata, lambda="ML")
mod6<- pgls(MAD.Shift.2...~ log(HWI), data = cdata, lambda="ML")
mod7<- pgls(MAD.Shift.2... ~ Mass.adj.WLI, data = cdata, lambda="ML")
mod8<- pgls(MAD.Shift.2...~ Body.mass..log. , data = cdata, lambda="ML")
mod9<- pgls(log(WAR) ~ log(Migration.Distance.Centroid), data = cdata, lambda="ML")
mod10<- pgls(log(HWI) ~ log(Migration.Distance.Centroid), data = cdata, lambda="ML")
mod11<- pgls(Body.mass..log. ~ log(Migration.Distance.Centroid), data = cdata, lambda="ML")
mod12<- pgls(MAD.Shift.2...~ log(Migration.Distance.Centroid), data = cdata, lambda="ML")
}
#Change mod# to view model statistics of desired model
summary(mod1)

#Model Selection
FULL <- pgls(MAD.Shift.2... ~ (log(WAR) + log(HWI) + Mass.adj.WLI + Body.mass..log. + log(Migration.Distance.Centroid)), data= cdata, lambda="ML")
D2 <- dredge(FULL, m.lim=c(1,1), trace=2, extra=c(R2=function(M) 1-M$RSSQ/M$NSSQ, L=function(M) M$param["lambda"]))
D2 


#Diagnostics
#plot.Linear.ModWar<-plot(resid(mod1), param$MAD.Shift.2...)
#plot.Linear.ModWar<-plot(resid(mod2), param$MAD.Shift.2...)
#plot.Linear.ModWar<-plot(resid(mod3), param$MAD.Shift.2...)
#plot.Linear.ModWar<-plot(resid(mod4), param$MAD.Shift.2...)

#Colinearity
para2<-subset(para, select= c(MAD.Shift.2...,WAR,HWI,Mass.adj.WLI,Body.mass..log.,Migration.Distance.Centroid))
para2$Migration.Distance.Centroid<-as.numeric(para2$Migration.Distance.Centroid)
pairs(para2)

cor(para2, method ="pearson")

colWAR<-lm(para$MAD.Shift.2... ~ log(para$WAR)  + log(para$Migration.Distance.Centroid))
colHWI<-lm(para$MAD.Shift.2... ~ log(para$HWI)  + log(para$Migration.Distance.Centroid))
colWLI<-lm(para$MAD.Shift.2... ~ log(para$WLI)  + log(para$Migration.Distance.Centroid))
colMss<-lm(para$MAD.Shift.2... ~ para$Body.mass..log.+ log(para$Migration.Distance.Centroid))


VIFWAR<-vif(colWAR)
VIFHWI<-vif(colHWI)
VIFWLI<-vif(colWLI)
VIFMSS<-vif(colMss)


#check WAR 
summary(lm(para$HWI ~ para$WAR))

#Code for grid specific mad
rawmad<-read.csv(file = "MAD_Output_All_Filter_1.csv", header =TRUE)


rawmad_groupedgrid <- rawmad %>% group_by(Year, GridID)
madsgroupedgrid <-rawmad_groupedgrid %>% summarise(meanmad=mean(MAD)) 

coefficientsGRID<- function (Grid){
  rawmad<-rawmad[rawmad$GridID == Grid,]
  regression<-lm(rawmad$MAD~rawmad$Year)
  MAD.Shift<-data.frame(coef(regression))
  GridID<-unique(rawmad$GridID)
  data.frame(GridID,MAD.Shift[2,])
} 

y<-unique(rawmad$GridID)
y

coefG.df<-lapply(y, coefficientsGRID)   
coefG2.df<-data.frame(do.call(rbind, coefG.df))




