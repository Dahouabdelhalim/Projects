
setwd("C:/Users/emchenoweth/Desktop/New folder/Step 4 Combine Costs and Gains")

results.Both$lunge.rate<-NA
results.Both$lunge.rate.se<-NA
results.Both$metabolic.costs_min<-NA
results.Both$metabolic.costs_min.se<-NA
results.Both$metabolic.costs_lunge<-NA
results.Both$metabolic.costs_lunge.se<-NA
results.Both$net.energy.gain_lunge<-NA #KJ/lunge
results.Both$net.energy.gain_min<-NA #KJ/min

Swim<-nrow(results.Both)+1
results.Both[Swim,]$Prey.Group<-"Swimming"
results.Both[Swim,]$WTgross.energy_gainkJ.Lunge.mean<-0
results.Both[Swim,]$metabolic.costs_lunge<-3*BMR*12.2
results.Both[Swim,]$metabolic.costs_min<-(3*BMR*60)
results.Both[Swim,]$net.energy.gain_min<-(-3*BMR*60)
results.Both[Swim,]$.energy_m3.mean<-0
results.Both[Swim,]$.prop_captured.mean<-0
results.Both$Prey.species<-as.character(results.Both$Prey.species)
results.Both[Swim,]$Prey.species<-"None"
results.Both$Prey.Group<-as.factor(results.Both$Prey.Group)

results.Both$Null<-results.Both$.energy_m3.mean*24.75114*0.84
results.Both$B<-results.Both$PC.energy_gain.mean
results.Both$C<-results.Both$WT.energy_gain.mean
results.Both$metabolic.costs_min<-NA
results.Both$lunge.rate<-NA
results.Both$net.energy.gain_min<-NA


results.Both$.energy_m3.mean*24.75114*0.84
results.Both$PC.energy_gain.mean

results.Both$metabolic.costs_lunge<-NA
results.Both[results.Both$Prey.species=="Chum salmon",]$metabolic.costs_lunge<-EnergyCosts[EnergyCosts$group == "ChumTB",]$meta.costs.mean
results.Both[results.Both$Prey.species=="Coho salmon",]$metabolic.costs_lunge<-EnergyCosts[EnergyCosts$group == "CohoKB",]$meta.costs.mean
results.Both[results.Both$Prey.Group=="KrillSS_subsurface",]$metabolic.costs_lunge<-mean(EnergyCosts[EnergyCosts$group == "KrillSS",]$meta.costs.mean)
results.Both[results.Both$Prey.Group=="KrillSC_subsurface",]$metabolic.costs_lunge<-EnergyCosts[EnergyCosts$group == "KrillSC",]$meta.costs.mean
results.Both[results.Both$Prey.Group=="Swimming",]$metabolic.costs_lunge<-(3*BMR*12.2)

results.Both[results.Both$Prey.species=="Chum salmon",]$lunge.rate<-EnergyCosts[EnergyCosts$group == "ChumTB",]$lunge_min
results.Both[results.Both$Prey.species=="Coho salmon",]$lunge.rate<-EnergyCosts[EnergyCosts$group == "CohoKB",]$lunge_min
results.Both[results.Both$Prey.Group=="KrillSS_subsurface",]$lunge.rate<-mean(EnergyCosts[EnergyCosts$group == "KrillSS",]$lunge_min)
results.Both[results.Both$Prey.Group=="KrillSC_subsurface",]$lunge.rate<-EnergyCosts[EnergyCosts$group == "KrillSC",]$lunge_min
results.Both[results.Both$Prey.Group=="Swimming",]$lunge.rate<-0

results.Both[results.Both$Prey.species=="Chum salmon",]$metabolic.costs_min<-EnergyCosts[EnergyCosts$group == "ChumTB",]$meta.costs_min
results.Both[results.Both$Prey.species=="Coho salmon",]$metabolic.costs_min<-EnergyCosts[EnergyCosts$group == "CohoKB",]$meta.costs_min
results.Both[results.Both$Prey.Group=="KrillSS_subsurface",]$metabolic.costs_min<-mean(EnergyCosts[EnergyCosts$group == "KrillSS",]$meta.costs_min)
results.Both[results.Both$Prey.Group=="KrillSC_subsurface",]$metabolic.costs_min<-EnergyCosts[EnergyCosts$group == "KrillSC",]$meta.costs_min

AllLunges$tag<-as.factor(AllLunges$tag)
LungeSC<-AllLunges[AllLunges$tag == "mn14_110a"|AllLunges$tag == "mn14_113a"|AllLunges$tag == "mn14_113b",]
DivesSC<-AllDives[AllDives$tag == "mn14_110a"|AllDives$tag == "mn14_113a"|AllDives$tag == "mn14_113b",]
results.Both[results.Both$Prey.Group=="HerringSC_subsurface",]$metabolic.costs_lunge<-mean(LungeSC[LungeSC$max.p>=22.5,]$meta.cost)
results.Both[results.Both$Prey.Group=="HerringSC_shallow",]$metabolic.costs_lunge<-mean(LungeSC[LungeSC$max.p<22.5,]$meta.cost)

results.Both[results.Both$Prey.Group=="HerringSC_shallow",]$metabolic.costs_min<-mean(DivesSC[DivesSC$p.max<22.5,]$total.meta.cost/(DivesSC[DivesSC$p.max<22.5,]$dur/60), na.rm = TRUE)
results.Both[results.Both$Prey.Group=="HerringSC_subsurface",]$metabolic.costs_min<-mean(DivesSC[DivesSC$p.max>=22.5,]$total.meta.cost/(DivesSC[DivesSC$p.max>=22.5,]$dur/60), na.rm = TRUE)
results.Both[results.Both$Prey.Group=="HerringSC_shallow",]$lunge.rate<-mean(DivesSC[DivesSC$p.max<22.5,]$lunges_min, na.rm = TRUE)
results.Both[results.Both$Prey.Group=="HerringSC_subsurface",]$lunge.rate<-mean(DivesSC[DivesSC$p.max>=22.5,]$lunges_min, na.rm = TRUE)


LungeTI1<-AllLunges[AllLunges$tag == "mn13_108a"|AllLunges$tag == "mn13_109a",]
DivesTI1<-AllDives[AllDives$tag == "mn13_108a"|AllDives$tag == "mn13_109a",]
results.Both[results.Both$Prey.Group=="HerringTI_small_subsurface",]$metabolic.costs_lunge<-mean(LungeTI1[LungeTI1$max.p>=22.5,]$meta.cost)
results.Both[results.Both$Prey.Group=="HerringTI_small_shallow",]$metabolic.costs_lunge<-mean(LungeTI1[LungeTI1$max.p<22.5,]$meta.cost)

results.Both[results.Both$Prey.Group=="HerringTI_small_shallow",]$metabolic.costs_min<-mean(DivesTI1[DivesTI1$p.max<22.5,]$total.meta.cost/(DivesTI1[DivesTI1$p.max<22.5,]$dur/60), na.rm = TRUE)
results.Both[results.Both$Prey.Group=="HerringTI_small_subsurface",]$metabolic.costs_min<-mean(DivesTI1[DivesTI1$p.max>=22.5,]$total.meta.cost/(DivesTI1[DivesTI1$p.max>=22.5,]$dur/60), na.rm = TRUE)
results.Both[results.Both$Prey.Group=="HerringTI_small_shallow",]$lunge.rate<-mean(DivesTI1[DivesTI1$p.max<22.5,]$lunges_min, na.rm = TRUE)
results.Both[results.Both$Prey.Group=="HerringTI_small_subsurface",]$lunge.rate<-mean(DivesTI1[DivesTI1$p.max>=22.5,]$lunges_min, na.rm = TRUE)

LungeTI2<-AllLunges[AllLunges$tag == "mn13_113d"|AllLunges$tag == "mn13_112a"|AllLunges$tag == "mn13_113a"|AllLunges$tag == "mn13_111a",]
DivesTI2<-AllDives[AllDives$tag == "mn13_113d"|AllDives$tag == "mn13_112a"|AllDives$tag == "mn13_113a"|AllDives$tag == "mn13_111a",]
results.Both[results.Both$Prey.Group=="HerringTI_large_subsurface",]$metabolic.costs_lunge<-mean(LungeTI2[LungeTI2$max.p>=22.5,]$meta.cost)
results.Both[results.Both$Prey.Group=="HerringTI_large_shallow",]$metabolic.costs_lunge<-mean(LungeTI2[LungeTI2$max.p<22.5,]$meta.cost)

results.Both[results.Both$Prey.Group=="HerringTI_large_shallow",]$metabolic.costs_min<-mean(DivesTI2[DivesTI2$p.max<22.5,]$total.meta.cost/(DivesTI2[DivesTI2$p.max<22.5,]$dur/60), na.rm = TRUE)
results.Both[results.Both$Prey.Group=="HerringTI_large_subsurface",]$metabolic.costs_min<-mean(DivesTI2[DivesTI2$p.max>=22.5,]$total.meta.cost/(DivesTI2[DivesTI2$p.max>=22.5,]$dur/60), na.rm = TRUE)
results.Both[results.Both$Prey.Group=="HerringTI_large_shallow",]$lunge.rate<-mean(DivesTI2[DivesTI2$p.max<22.5,]$lunges_min, na.rm = TRUE)
results.Both[results.Both$Prey.Group=="HerringTI_large_subsurface",]$lunge.rate<-mean(DivesTI2[DivesTI2$p.max>=22.5,]$lunges_min, na.rm = TRUE)

LungeHerringAll<-AllLunges[AllLunges$tag == "mn13_108a"|AllLunges$tag == "mn13_109a"|AllLunges$tag == "mn14_110a"|AllLunges$tag == "mn14_113a"|AllLunges$tag == "mn14_113b",]
DivesHerringAll<-AllDives[AllDives$tag == "mn13_108a"|AllDives$tag == "mn13_109a"|AllDives$tag == "mn14_110a"|AllDives$tag == "mn14_113a"|AllDives$tag == "mn14_113b",]
results.Both[results.Both$Prey.Group=="HerringKB_subsurface",]$metabolic.costs_lunge<-mean(LungeHerringAll[LungeHerringAll$max.p>=22.5,]$meta.cost)
results.Both[results.Both$Prey.Group=="HerringKB_shallow",]$metabolic.costs_lunge<-mean(LungeHerringAll[LungeHerringAll$max.p<22.5,]$meta.cost)

results.Both[results.Both$Prey.Group=="HerringKB_shallow",]$metabolic.costs_min<-mean(DivesHerringAll[DivesHerringAll$p.max<22.5,]$total.meta.cost/(DivesHerringAll[DivesHerringAll$p.max<22.5,]$dur/60), na.rm = TRUE)
results.Both[results.Both$Prey.Group=="HerringKB_subsurface",]$metabolic.costs_min<-mean(DivesHerringAll[DivesHerringAll$p.max>=22.5,]$total.meta.cost/(DivesHerringAll[DivesHerringAll$p.max>=22.5,]$dur/60), na.rm = TRUE)
results.Both[results.Both$Prey.Group=="HerringKB_shallow",]$lunge.rate<-mean(DivesHerringAll[DivesHerringAll$p.max<22.5,]$lunges_min, na.rm = TRUE)
results.Both[results.Both$Prey.Group=="HerringKB_subsurface",]$lunge.rate<-mean(DivesHerringAll[DivesHerringAll$p.max>=22.5,]$lunges_min, na.rm = TRUE)

results.Both$D<-results.Both$WT.energy_gain.mean-results.Both$metabolic.costs_lunge
results.Both$net.energy.gain_min<-results.Both$WTgross.energy_gainkJ.Lunge.mean*results.Both$lunge.rate-results.Both$metabolic.costs_min
results.Both$E<-results.Both$net.energy.gain_min
results.Both
results.Both[order(-results.Both$E),]

