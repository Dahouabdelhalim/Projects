
#Reset working directory to source Location
setwd("C:/Users/emchenoweth/Desktop/New folder/Step 2 Acoustics")
library("dplyr", lib.loc="~/R/win-library/3.5")

CS2014Chum<-read.csv("CS2014Chum.csv")
CS2014Coho<-read.csv("CS2014Coho.csv")
KB2014Herring<-read.csv("KB2014Herring.csv")
SC2014Herring<-read.csv("SC2014Herring.csv")
SC2014Krill<-read.csv("SC2014Krill.csv")
SS2012Krill<-read.csv("SS2012Krill.csv")
TI2013Herring_large<-read.csv("TI2013Herring_large.csv")
TI2013Herring_small<-read.csv("TI2013Herring_small.csv")

#Function
Energy_gain.nodf<-function(Vp = NA, Vp.sd = NA, Ddetect = NA, Ddetect.sd= NA, assim_eff = 0.84,
                           mass_ind,energy_g, Sv_mean,TS,lungetype = NA,var = var,Prey = NA, 
                           Data = 1,Vc_0 = NA, proportion_captured = "formula", grid = NA){
  if(proportion_captured == "formula"){
    df<-Create.df(var = var, Prey = Prey, Vc_0 = Vc_0)
    Capture.prop.stochastic(Vp=Vp, Ddetect=Ddetect, Vp.sd=Vp.sd, Ddetect.sd=Ddetect.sd, df=df, lungetype=lungetype, var=var, grid = grid)
  }
  if(proportion_captured != "formula"){prop_captured<-proportion_captured}
  engulf_vol_m3<-(8.0094*Lbody ^3.21)/1020 # equation from GOldbogen et al 2012 for kg.  / 1020 to accoutn for weight of seawater
  
  Data<-as.data.frame(Data)
  Data$.ind_m3<-10^(Sv_mean/10)/10^(TS/10)
  Data$.ind_captured<-10^(Sv_mean/10)/10^(TS/10)*prop_captured*engulf_vol_m3
  Data$.biomass_m3<-10^(Sv_mean/10)/10^(TS/10)*mass_ind
  Data$.biomass_captured<- Data$.ind_captured*mass_ind
  Data$.energy_gain<- Data$.biomass_captured*energy_g*assim_eff
  Data$.prop_captured<-prop_captured
  #Data$.Vp<-Vp
  Data$.mass_ind<-mass_ind
  Data$.energy_g<-energy_g
  Data$.lungetype<-lungetype
  Data$.prey<-Prey
  print(Data)
}

#No Whales

Energygain<-c()
#if(!is.na(EnergyParameters[n,]$CaptureProp_surface)){
  for(n in 1:nrow(EnergyParameters)){
    n
    assign(paste(EnergyParameters[n,]$Prey.Group),
           Energy_gain.nodf(mass_ind = EnergyParameters[n,]$Mean.Wet.mass,
                            energy_g = EnergyParameters[n,]$Energy,
                            Sv_mean = get(paste(EnergyParameters[n,]$Sv))$Sv_mean,
                            proportion_captured = EnergyParameters[n,]$CaptureProp_NoWhaleBehave,
                            TS = EnergyParameters[n,]$TS,
                            var = var,
                            lungetype = "subsurface",
                            Prey = EnergyParameters[n,]$Prey,
                            Data = get(paste(EnergyParameters[n,]$Data)) ))
  }
#}

HerringTI_small_shallow$group<-"HerringTI_small_shallow"
HerringTI_large_shallow$group<-"HerringTI_large_shallow"
HerringSC_shallow$group<-"HerringSC_shallow"
HerringKB_shallow$group<-"HerringKB_shallow"
KrillSS_subsurface$group<-"KrillSS_subsurface"
KrillSC_subsurface$group<-"KrillSC_subsurface"
ChumEarly_shallow$group<-"ChumEarly_shallow"
ChumLate_shallow$group<-"ChumLate_shallow"
CohoAtRelease_shallow$group<-"CohoAtRelease_shallow"
CohoAfterRelease_shallow$group<-"CohoAfterRelease_shallow"
HerringTI_small_subsurface$group<-"HerringTI_small_subsurface"
HerringTI_large_subsurface$group<-"HerringTI_large_subsurface"
HerringSC_subsurface$group<-"HerringSC_subsurface"
HerringKB_subsurface$group<-"HerringKB_subsurface"

EnergyAll<-rbind.fill(HerringTI_small_shallow, 
                      HerringTI_large_shallow,
                      HerringSC_shallow, 
                      KrillSS_subsurface, 
                      KrillSC_subsurface, 
                      ChumEarly_shallow,
                      ChumLate_shallow, 
                      CohoAtRelease_shallow, 
                      CohoAfterRelease_shallow, 
                      HerringKB_shallow,
                      HerringTI_small_subsurface,
                      HerringTI_large_subsurface,
                      HerringSC_subsurface,
                      HerringKB_subsurface
                      )

EnergyAll<-EnergyAll[log(EnergyAll$.energy_gain)>0,]
EnergyAll$group<-as.factor(EnergyAll$group)
#ggplot()+geom_boxplot(aes(group,.energy_gain), data = EnergyAll)+scale_y_log10("Energy Captured (kJ)",labels = scales::comma)
#ggplot()+geom_point(data = EnergyAll,aes(group,.prop_captured, color = group, size = 3))
#ggplot()+geom_point(data = EnergyAll,aes(.biomass_m3,.energy_gain, color = group))+geom_smooth(aes(.biomass_m3,.energy_gain),method = "lm", data = EnergyAll)
#ggplot()+geom_point(data = EnergyAll,aes(.energy_g,.energy_gain, color = group))+geom_smooth(aes(.energy_g,.energy_gain),method = "lm", data = EnergyAll)
#ggplot()+geom_point(data = EnergyAll,aes(.energy_g,log(.energy_gain), color = group))+geom_smooth(aes(.energy_g,log(.energy_gain)),method = "lm", data = EnergyAll)

EnergyAll$.energy_m3<-EnergyAll$.biomass_m3*EnergyAll$.energy_g
results.NoWhales<-as.data.frame(summaryBy(.energy_m3+.prop_captured+.energy_gain+Sv_mean+.ind_m3+.ind_captured+.biomass_m3+.biomass_captured ~group, FUN = c(mean,se), data = EnergyAll, na.rm = TRUE))
# if this doesn't run its because somewhere there is code introducing "mean" as a variable, which conflicts with its duties as a function.  Sometimes it happens sometimes
#not.  Just rm(mean) and the thing runs. Frustrating.
results.NoWhales<-results.NoWhales[order(-results.NoWhales$.energy_gain.mean),]
results.NoWhales$gross.energy_gainkJ.Lunge.mean<-results.NoWhales$.energy_gain.mean
results.NoWhales
#with CaptureProp_NoWhaleBehave


###Whale Tactics
EnergyParameters$CaptureProp_WhaleTactics<-NA
EnergyParameters[EnergyParameters$depth =="shallow",]$CaptureProp_WhaleTactics<-EnergyParameters[EnergyParameters$depth =="shallow",]$CaptureProp_surface
EnergyParameters[EnergyParameters$depth =="subsurface",]$CaptureProp_WhaleTactics<-EnergyParameters[EnergyParameters$depth =="subsurface",]$CaptureProp_subsurface



Energygain<-c()
for(n in 1:nrow(EnergyParameters)){
  n
  assign(paste0(EnergyParameters[n,]$Prey.Group,"_WT"),
         Energy_gain.nodf(mass_ind = EnergyParameters[n,]$Mean.Wet.mass,
                          energy_g = EnergyParameters[n,]$Energy,
                          Sv_mean = get(paste(EnergyParameters[n,]$Sv))$Sv_mean,
                          proportion_captured = EnergyParameters[n,]$CaptureProp_WhaleTactics,
                          TS = EnergyParameters[n,]$TS,
                          var = var,
                          lungetype = EnergyParameters[n,]$depth,
                          Prey = EnergyParameters[n,]$Prey,
                          Data = get(paste(EnergyParameters[n,]$Data)) ))
}
#}



HerringTI_small_shallow_WT$group<-"HerringTI_small_shallow"
HerringTI_large_shallow_WT$group<-"HerringTI_large_shallow"
HerringSC_shallow_WT$group<-"HerringSC_shallow"
HerringKB_shallow_WT$group<-"HerringKB_shallow"
KrillSS_subsurface_WT$group<-"KrillSS_subsurface"
KrillSC_subsurface_WT$group<-"KrillSC_subsurface"
ChumEarly_shallow_WT$group<-"ChumEarly_shallow"
ChumLate_shallow_WT$group<-"ChumLate_shallow"
CohoAtRelease_shallow_WT$group<-"CohoAtRelease_shallow"
CohoAfterRelease_shallow_WT$group<-"CohoAfterRelease_shallow"
HerringTI_small_subsurface_WT$group<-"HerringTI_small_subsurface"
HerringTI_large_subsurface_WT$group<-"HerringTI_large_subsurface"
HerringSC_subsurface_WT$group<-"HerringSC_subsurface"
HerringKB_subsurface_WT$group<-"HerringKB_subsurface"

EnergyAll<-rbind.fill(HerringTI_small_shallow_WT, 
                      HerringTI_large_shallow_WT,
                      HerringSC_shallow_WT, 
                      KrillSS_subsurface_WT, 
                      KrillSC_subsurface_WT, 
                      ChumEarly_shallow_WT,
                      ChumLate_shallow_WT, 
                      CohoAtRelease_shallow_WT, 
                      CohoAfterRelease_shallow_WT, 
                      HerringKB_shallow_WT,
                      HerringTI_small_subsurface_WT,
                      HerringTI_large_subsurface_WT,
                      HerringSC_subsurface_WT,
                      HerringKB_subsurface_WT
)

EnergyAll<-EnergyAll[log(EnergyAll$.energy_gain)>0,]
EnergyAll$group<-as.factor(EnergyAll$group)
#ggplot()+geom_boxplot(aes(group,.energy_gain), data = EnergyAll)+scale_y_log10("Energy Captured (kJ)",labels = scales::comma)
#ggplot()+geom_point(data = EnergyAll,aes(group,.prop_captured, color = group, size = 3))
#ggplot()+geom_point(data = EnergyAll,aes(.biomass_m3,.energy_gain, color = group))+geom_smooth(aes(.biomass_m3,.energy_gain),method = "lm", data = EnergyAll)
#ggplot()+geom_point(data = EnergyAll,aes(.energy_g,.energy_gain, color = group))+geom_smooth(aes(.energy_g,.energy_gain),method = "lm", data = EnergyAll)
#ggplot()+geom_point(data = EnergyAll,aes(.energy_g,log(.energy_gain), color = group))+geom_smooth(aes(.energy_g,log(.energy_gain)),method = "lm", data = EnergyAll)

EnergyAll$.energy_m3<-EnergyAll$.biomass_m3*EnergyAll$.energy_g
results.Tactics<-as.data.frame(summaryBy(.energy_m3+.prop_captured+.energy_gain+Sv_mean+.ind_m3+.ind_captured+.biomass_m3+.biomass_captured ~group, FUN = c(mean,se), data = EnergyAll, na.rm = TRUE))
results.Tactics<-results.Tactics[order(-results.Tactics$.energy_gain.mean),]
results.Tactics$gross.energy_gainkJ.Lunge.mean<-results.Tactics$.energy_gain.mean
results.Tactics
names(results.Tactics)[which(names(results.Tactics)=="group")]<-"Prey.Group"
names(results.NoWhales)[which(names(results.NoWhales)=="group")]<-"Prey.Group"
results.Tactics<-results.Tactics[order(results.Tactics$Prey.Group),]
results.NoWhales<-results.NoWhales[order(results.NoWhales$Prey.Group),]

#use to compare columns and find unique (FALSE) and duplicates (TRUE)
results.NoWhales==results.Tactics

KeepDup<-names(results.Tactics)[which(names(results.Tactics)%in% c(".energy_gain.mean",".ind_captured.mean",".biomass_captured.mean",".energy_gain.se",".ind_captured.se",".biomass_captured.se","gross.energy_gainkJ.Lunge.mean"))]
names(results.Tactics)[which(names(results.Tactics)%in% c(".energy_gain.mean",".ind_captured.mean",".biomass_captured.mean",".energy_gain.se",".ind_captured.se",".biomass_captured.se","gross.energy_gainkJ.Lunge.mean"))]<-paste0("WT",KeepDup)

KeepDup2<-names(results.NoWhales)[which(names(results.NoWhales)%in% c(".energy_gain.mean",".ind_captured.mean",".biomass_captured.mean",".energy_gain.se",".ind_captured.se",".biomass_captured.se","gross.energy_gainkJ.Lunge.mean"))]
names(results.NoWhales)[which(names(results.NoWhales)%in% c(".energy_gain.mean",".ind_captured.mean",".biomass_captured.mean",".energy_gain.se",".ind_captured.se",".biomass_captured.se","gross.energy_gainkJ.Lunge.mean"))]<-paste0("PC",KeepDup2)
results.NoWhales %>% select(starts_with("WT"))

InandOut<-left_join(EnergyParameters, results.NoWhales)
results.Both<-left_join(InandOut,results.Tactics %>% select(starts_with("WT"), matches("Prey.Group")), by = "Prey.Group")










nrow(EnergyAll[EnergyAll$tag== "mn13_111a" | EnergyAll$tag== "mn13_112a"| EnergyAll$tag== "mn13_113a"| EnergyAll$tag== "mn13_113d" & EnergyAll$p<22.5, ])
nrow(EnergyAll[EnergyAll$tag== "mn13_108a" | EnergyAll$tag== "mn13_109a" & EnergyAll$p<22.5, ])

nrow(EnergyAll[EnergyAll$tag== "mn13_111a" | EnergyAll$tag== "mn13_112a"| EnergyAll$tag== "mn13_113a"| EnergyAll$tag== "mn13_113d" & EnergyAll$p>=22.5, ])
nrow(EnergyAll[EnergyAll$tag== "mn13_108a" | EnergyAll$tag== "mn13_109a" & !is.na(EnergyAll$tag) & EnergyAll$p>=22.5, ])



####### test proportion captured####

a<-3
b<-10
Create.df(var = var,Prey = EnergyParameters[a,]$Prey, Vc_0 = EnergyParameters[a,]$Whale.initial.speed..subsurface.)
Capture.prop.stochastic(Vp = EnergyParameters[a,]$Prey.speed,
                        Vp.sd = 0,
                        Ddetect = EnergyParameters[a,]$Mean.Reaction.distance,
                        Ddetect.sd = 0,
                        df = df,
                        lungetype="subsurface",
                        Tclose=Tclose,
                        var = var)
prop_captured

##### simulate speed distribution ###

#HerringTI
s<-rnorm(1000,0.217,0.021)*rnorm(1000,5.746,1.4)
mean(s)
sd(s)
hist(s)

#HerringSC
s<-rnorm(1000,0.202,0.008)*rnorm(1000,5.746,1.4)
mean(s)
sd(s)
hist(s)

#HerringKB
s<-rnorm(1000,0.196,0.020)*rnorm(1000,5.746,1.4)
mean(s)
sd(s)
hist(s)

#KrillSS
s<-rnorm(1000,0.016,0.004)*rnorm(1000,9.5,1.5)
mean(s)
sd(s)
hist(s)

#KrillSC
s<-rnorm(1000,0.018,0.003)*rnorm(1000,9.5,1.5)
mean(s)
sd(s)
hist(s)

#ChumEarly
s<-rnorm(1000,0.058,0.005)*rnorm(1000,3.8,3.2)
mean(s)
sd(s)
hist(s)

#ChumEarly
s<-rnorm(1000,0.066,0.009)*rnorm(1000,3.8,3.2)
mean(s)
sd(s)
hist(s)

#CohoAtRelease
s<-rnorm(1000,0.121,0.006)*rnorm(1000,8.1,6.6)
mean(s)
sd(s)
hist(s)

#CohoAfterRelease
s<-rnorm(1000,0.128,0.011)rnorm(1000,8.1,6.6)
mean(s)
sd(s)
hist(s)



######
EnergyParameters$PreyEscape<-EnergyParameters$Mean.Reaction.distance*EnergyParameters$Prey.speed
EnergyParameters$PreyEscapeSD<-sqrt(EnergyParameters$Mean.Reaction.distance^2*EnergyParameters$Prey.speed^2*+EnergyParameters$Prey.speed^2*EnergyParameters$SD.Reaction.distance^2+EnergyParameters$SD.Prey.speed^2*EnergyParameters$Mean.Reaction.distance^2)
mod<-(CaptureProp_surface~PreyEscape+PreyEscapeSD+Whale.initial.speed..surface., data = EnergyParameters)

#####

Energy_gain.nodf<-function(Vp = NA, Vp.sd = NA, Ddetect = NA, Ddetect.sd= NA, assim_eff = 0.84,
                           mass_ind,energy_g, Sv_mean,TS,lungetype = NA,var = var,Prey = NA, 
                           Data = 1,Vc_0 = NA, proportion_captured = "formula", grid = NA)
  
  
  for (n in 1:nrow(EnergyParameters)){
    Create.df(var = var,Prey = EnergyParameters[a,]$Prey, Vc_0 = mean(EnergyParameters$Whale.initial.speed..subsurface., na.rm = TRUE))
    Capture.prop.stochastic(Vp = EnergyParameters[a,]$Prey.speed,
                            Vp.sd = EnergyParameters[a,]$SD.Prey.speed,
                            Ddetect = EnergyParameters[a,]$Mean.Reaction.distance,
                            Ddetect.sd = EnergyParameters[a,]$SD.Reaction.distance,
                            df = df,
                            lungetype="subsurface",
                            Tclose=Tclose,
                            var = var)
  }

for (n in 1:nrow(EnergyParameters)){
  Energy_gain.nodf(mass_ind = EnergyParameters[n,]$Mean.Wet.mass,
                   energy_g = EnergyParameters[n,]$Energy,
                   Sv_mean = get(EnergyParameters[n,]$Sv)$Sv_mean,
                   proportion_captured = EnergyParameters[n,]$CaptureProp_NoWhaleBehave,
                   TS = EnergyParameters[n,]$TS,
                   var = var,
                   lungetype = "subsurface",
                   Prey = EnergyParameters[n,]$Prey)}



