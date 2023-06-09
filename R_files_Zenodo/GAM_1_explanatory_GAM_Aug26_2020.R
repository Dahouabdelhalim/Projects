library("mgcv")
using <- readRDS("0.2 explanatory_bl_sighting_180s_15min_eff1 in1km June3 2020.rds")

using <- as.data.frame(sapply( using, as.numeric ))
using <- using[complete.cases(using), ]

env.1 <-gam(has_bird ~ s(whole_chl,k=4) + s(whole_sal,k=4) + s(whole_sst,k=4)  + s(max_cur_speed,k=4) 
                  + s(bathymetry,k=4) + s(aspect,k=4) +s(dist2coast,k=4)
                  +s(seastate,k=4),data = using,family=binomial)
env.1.1<-update(env.1,.~.- s(max_cur_speed,k=4)) 
env.1.2<-update(env.1,.~.- s(dist2coast,k=4)) 
env.1.3<-update(env.1,.~.- s(bathymetry,k=4)) 
env.1.4<-update(env.1,.~.- s(aspect,k=4))
env.1.5<-update(env.1,.~.- s(whole_sst,k=4))
env.1.6<-update(env.1,.~.- s(whole_sal,k=4)) 
env.1.7<-update(env.1,.~.- s(whole_chl,k=4)) 

env.2<-update(env.1,.~.- s(aspect,k=4)+ s(anomalies,k=4)) 

env.2.1<-update(env.1,.~.- s(aspect,k=4)+ s(anomalies,k=4)- s(max_cur_speed,k=4))
env.2.2<-update(env.1,.~.- s(aspect,k=4)+ s(anomalies,k=4)- s(dist2coast,k=4))
env.2.3<-update(env.1,.~.- s(aspect,k=4)+ s(anomalies,k=4)- s(bathymetry,k=4)) 
env.2.4<-update(env.1,.~.- s(aspect,k=4)+ s(anomalies,k=4)- s(anomalies,k=4))
env.2.5<-update(env.1,.~.- s(aspect,k=4)+ s(anomalies,k=4)- s(whole_sst,k=4))
env.2.6<-update(env.1,.~.- s(aspect,k=4)+ s(anomalies,k=4)- s(whole_sal,k=4)) 
env.2.7<-update(env.1,.~.- s(aspect,k=4)+ s(anomalies,k=4)- s(whole_chl,k=4))


env.3<-update(env.1,.~.- s(aspect,k=4)+ s(roughness,k=4))

env.3.1<-update(env.3,.~. - s(max_cur_speed,k=4))
env.3.2<-update(env.3,.~.- s(dist2coast,k=4))
env.3.3<-update(env.3,.~.- s(bathymetry,k=4)) 
env.3.4<-update(env.3,.~.- s(roughness,k=4))
env.3.5<-update(env.3,.~.- s(whole_sst,k=4))
env.3.6<-update(env.3,.~.- s(whole_sal,k=4)) 
env.3.7<-update(env.3,.~.- s(whole_chl,k=4)) 


env.4<-update(env.1,.~.- s(max_cur_speed,k=4)+ s(stratf_indx,k=4))

env.4.1<-update(env.4,.~.- s(stratf_indx,k=4)) 
env.4.2<-update(env.4,.~.- s(dist2coast,k=4)) 
env.4.3<-update(env.4,.~.- s(bathymetry,k=4)) 
env.4.4<-update(env.4,.~.- s(aspect,k=4))
env.4.5<-update(env.4,.~.- s(whole_sst,k=4))
env.4.6<-update(env.4,.~.- s(whole_sal,k=4)) 
env.4.7<-update(env.4,.~.- s(whole_chl,k=4))

env.5<-update(env.2,.~.- s(max_cur_speed,k=4)+ s(stratf_indx,k=4))

env.5.1<-update(env.2,.~.- s(max_cur_speed,k=4)+ s(stratf_indx,k=4)- s(stratf_indx,k=4))
env.5.2<-update(env.2,.~.- s(max_cur_speed,k=4)+ s(stratf_indx,k=4)- s(dist2coast,k=4))
env.5.3<-update(env.2,.~.- s(max_cur_speed,k=4)+ s(stratf_indx,k=4)- s(bathymetry,k=4))
env.5.4<-update(env.2,.~.- s(max_cur_speed,k=4)+ s(stratf_indx,k=4)- s(anomalies,k=4))
env.5.5<-update(env.2,.~.- s(max_cur_speed,k=4)+ s(stratf_indx,k=4)- s(whole_sst,k=4))
env.5.6<-update(env.2,.~.- s(max_cur_speed,k=4)+ s(stratf_indx,k=4)- s(whole_sal,k=4)) 
env.5.7<-update(env.2,.~.- s(max_cur_speed,k=4)+ s(stratf_indx,k=4)- s(whole_chl,k=4))

env.6<-update(env.3,.~.- s(max_cur_speed,k=4)+ s(stratf_indx,k=4))

env.6.1<-update(env.6,.~.- s(stratf_indx,k=4))
env.6.2<-update(env.6,.~.- s(dist2coast,k=4))
env.6.3<-update(env.6,.~.- s(bathymetry,k=4)) 
env.6.4<-update(env.6,.~.- s(roughness,k=4))
env.6.5<-update(env.6,.~.- s(whole_sst,k=4))
env.6.6<-update(env.6,.~.- s(whole_sal,k=4)) 
env.6.7<-update(env.6,.~.- s(whole_chl,k=4))

aics<-AIC(env.1,env.1.1,env.1.2,env.1.3,env.1.4,env.1.5,env.1.6,env.1.7,
          env.2,env.2.1,env.2.2,env.2.3,env.2.4,env.2.5,env.2.6,env.2.7,
          env.3,env.3.1,env.3.2,env.3.3,env.3.4,env.3.5,env.3.6,env.3.7,
          env.4,env.4.1,env.4.2,env.4.3,env.4.4,env.4.5,env.4.6,env.4.7,
          env.5,env.5.1,env.5.2,env.5.3,env.5.4,env.5.5,env.5.6,env.5.7,
          env.6,env.6.1,env.6.2,env.6.3,env.6.4,env.6.5,env.6.6,env.6.7)[2]
delta<-aics-min(aics) 
delta[order(delta$AIC),, drop = FALSE]


env.4.5.1<-update(env.4.5,.~.- s(whole_chl, k = 4)) 
env.4.5.2<-update(env.4.5,.~.- s(whole_sal, k = 4)) 
env.4.5.3<-update(env.4.5,.~.- s(bathymetry, k = 4)) 
env.4.5.4<-update(env.4.5,.~.- s(aspect, k = 4)) 
env.4.5.5<-update(env.4.5,.~.- s(dist2coast, k = 4)) 
env.4.5.6<-update(env.4.5,.~.- s(stratf_indx, k = 4)) 

env.5.5.1<-update(env.5.5,.~.- s(whole_chl, k = 4)) 
env.5.5.2<-update(env.5.5,.~.- s(whole_sal, k = 4)) 
env.5.5.3<-update(env.5.5,.~.- s(bathymetry, k = 4)) 
env.5.5.4<-update(env.5.5,.~.- s(anomalies, k = 4)) 
env.5.5.5<-update(env.5.5,.~.- s(dist2coast, k = 4)) 
env.5.5.6<-update(env.5.5,.~.- s(stratf_indx, k = 4)) 


aics<-AIC(env.1,env.1.1,env.1.2,env.1.3,env.1.4,env.1.5,env.1.6,env.1.7,
          env.2,env.2.1,env.2.2,env.2.3,env.2.4,env.2.5,env.2.6,env.2.7,
          env.3,env.3.1,env.3.2,env.3.3,env.3.4,env.3.5,env.3.6,env.3.7,
          env.4,env.4.1,env.4.2,env.4.3,env.4.4,env.4.5,env.4.6,env.4.7,
          env.5,env.5.1,env.5.2,env.5.3,env.5.4,env.5.5,env.5.6,env.5.7,
          env.6,env.6.1,env.6.2,env.6.3,env.6.4,env.6.5,env.6.6,env.6.7,
          env.4.5.1,env.4.5.2,env.4.5.3,env.4.5.4,env.4.5.5,env.4.5.6,
          env.5.5.1,env.5.5.2,env.5.5.3,env.5.5.4,env.5.5.5,env.5.5.6)[2]

delta<-aics-min(aics) 
delta[order(delta$AIC),, drop = FALSE]

fish.1 <-mgcv::gam(has_bird ~ s(MAC, k = 4) + s(SPR, k = 4) + s(ANE, k = 4) + s(PIL, k = 4) + s(HOM, k = 4) + s(HER, k = 4) +s(BOF, k = 4)
                   +s(seastate, k = 4),data = using,family=binomial)


fish.1.1<-update(fish.1,.~.- s(MAC, k = 4)) 
fish.1.2<-update(fish.1,.~.- s(SPR, k = 4)) 
fish.1.3<-update(fish.1,.~.- s(ANE, k = 4)) 
fish.1.4<-update(fish.1,.~.- s(PIL, k = 4)) 
fish.1.5<-update(fish.1,.~.- s(HOM, k = 4)) 
fish.1.6<-update(fish.1,.~.- s(HER, k = 4)) 
fish.1.7<-update(fish.1,.~.- s(BOF, k = 4)) 

#1.1
fish.1.1.2<-update(fish.1.1,.~.- s(SPR, k = 4)) 
fish.1.1.3<-update(fish.1.1,.~.- s(ANE, k = 4)) 
fish.1.1.4<-update(fish.1.1,.~.- s(PIL, k = 4)) 
fish.1.1.5<-update(fish.1.1,.~.- s(HOM, k = 4)) 
fish.1.1.6<-update(fish.1.1,.~.- s(HER, k = 4)) 
fish.1.1.7<-update(fish.1.1,.~.- s(BOF, k = 4)) 

#1.2
fish.1.2.2<-update(fish.1.2,.~.- s(MAC, k = 4)) 
fish.1.2.3<-update(fish.1.2,.~.- s(ANE, k = 4)) 
fish.1.2.4<-update(fish.1.2,.~.- s(PIL, k = 4)) 
fish.1.2.5<-update(fish.1.2,.~.- s(HOM, k = 4)) 
fish.1.2.6<-update(fish.1.2,.~.- s(HER, k = 4)) 
fish.1.2.7<-update(fish.1.2,.~.- s(BOF, k = 4)) 

#1.3
fish.1.3.2<-update(fish.1.3,.~.- s(MAC, k = 4)) 
fish.1.3.3<-update(fish.1.3,.~.- s(SPR, k = 4)) 
fish.1.3.4<-update(fish.1.3,.~.- s(PIL, k = 4)) 
fish.1.3.5<-update(fish.1.3,.~.- s(HOM, k = 4)) 
fish.1.3.6<-update(fish.1.3,.~.- s(HER, k = 4)) 
fish.1.3.7<-update(fish.1.3,.~.- s(BOF, k = 4))

#1.4
fish.1.4.2<-update(fish.1.4,.~.- s(MAC, k = 4)) 
fish.1.4.3<-update(fish.1.4,.~.- s(SPR, k = 4)) 
fish.1.4.4<-update(fish.1.4,.~.- s(ANE, k = 4)) 
fish.1.4.5<-update(fish.1.4,.~.- s(HOM, k = 4)) 
fish.1.4.6<-update(fish.1.4,.~.- s(HER, k = 4)) 
fish.1.4.7<-update(fish.1.4,.~.- s(BOF, k = 4))

#1.5
fish.1.5.2<-update(fish.1.5,.~.- s(MAC, k = 4)) 
fish.1.5.3<-update(fish.1.5,.~.- s(SPR, k = 4)) 
fish.1.5.4<-update(fish.1.5,.~.- s(ANE, k = 4)) 
fish.1.5.5<-update(fish.1.5,.~.- s(PIL, k = 4)) 
fish.1.5.6<-update(fish.1.5,.~.- s(HER, k = 4)) 
fish.1.5.7<-update(fish.1.5,.~.- s(BOF, k = 4))

#1.6
fish.1.6.2<-update(fish.1.6,.~.- s(MAC, k = 4)) 
fish.1.6.3<-update(fish.1.6,.~.- s(SPR, k = 4)) 
fish.1.6.4<-update(fish.1.6,.~.- s(ANE, k = 4)) 
fish.1.6.5<-update(fish.1.6,.~.- s(PIL, k = 4)) 
fish.1.6.6<-update(fish.1.6,.~.- s(HOM, k = 4)) 
fish.1.6.7<-update(fish.1.6,.~.- s(BOF, k = 4))

#1.7
fish.1.7.2<-update(fish.1.7,.~.- s(MAC, k = 4)) 
fish.1.7.3<-update(fish.1.7,.~.- s(SPR, k = 4)) 
fish.1.7.4<-update(fish.1.7,.~.- s(ANE, k = 4)) 
fish.1.7.5<-update(fish.1.7,.~.- s(PIL, k = 4)) 
fish.1.7.6<-update(fish.1.7,.~.- s(HOM, k = 4)) 
fish.1.7.7<-update(fish.1.7,.~.- s(HER, k = 4))

aics<-AIC(fish.1,fish.1.1,fish.1.2,fish.1.3,fish.1.4,fish.1.5,fish.1.6,fish.1.7,
          fish.1.1.2,fish.1.1.3,fish.1.1.4,fish.1.1.5,fish.1.1.6,fish.1.1.7,
          fish.1.2.2,fish.1.2.3,fish.1.2.4,fish.1.2.5,fish.1.2.6,fish.1.2.7,
          fish.1.3.2,fish.1.3.3,fish.1.3.4,fish.1.3.5,fish.1.3.6,fish.1.3.7,
          fish.1.4.2,fish.1.4.3,fish.1.4.4,fish.1.4.5,fish.1.4.6,fish.1.4.7,
          fish.1.5.2,fish.1.5.3,fish.1.5.4,fish.1.5.5,fish.1.5.6,fish.1.5.7,
          fish.1.6.2,fish.1.6.3,fish.1.6.4,fish.1.6.5,fish.1.6.6,fish.1.6.7,
          fish.1.7.2,fish.1.7.3,fish.1.7.4,fish.1.7.5,fish.1.7.6,fish.1.7.7)[2]
delta<-aics-min(aics) 
delta[order(delta$AIC),, drop = FALSE]

fish.1.2.4.1<-update(fish.1.2.4,.~.- s(MAC, k = 4))
fish.1.2.4.2<-update(fish.1.2.4,.~.- s(ANE, k = 4))
fish.1.2.4.3<-update(fish.1.2.4,.~.- s(HOM, k = 4))
fish.1.2.4.4<-update(fish.1.2.4,.~.- s(HER, k = 4))
fish.1.2.4.5<-update(fish.1.2.4,.~.- s(BOF, k = 4))

fish.1.4.3.1<-update(fish.1.4.3,.~.- s(MAC, k = 4))
fish.1.4.3.2<-update(fish.1.4.3,.~.- s(ANE, k = 4))
fish.1.4.3.3<-update(fish.1.4.3,.~.- s(HOM, k = 4))
fish.1.4.3.4<-update(fish.1.4.3,.~.- s(HER, k = 4))
fish.1.4.3.5<-update(fish.1.4.3,.~.- s(BOF, k = 4))

fish.1.4.6.1<-update(fish.1.4.6,.~.- s(MAC, k = 4))
fish.1.4.6.2<-update(fish.1.4.6,.~.- s(ANE, k = 4))
fish.1.4.6.3<-update(fish.1.4.6,.~.- s(HOM, k = 4))
fish.1.4.6.4<-update(fish.1.4.6,.~.- s(BOF, k = 4))
fish.1.4.6.5<-update(fish.1.4.6,.~.- s(SPR, k = 4))

aics<-AIC(fish.1,fish.1.1,fish.1.2,fish.1.3,fish.1.4,fish.1.5,fish.1.6,fish.1.7,
          fish.1.1.2,fish.1.1.3,fish.1.1.4,fish.1.1.5,fish.1.1.6,fish.1.1.7,
          fish.1.2.2,fish.1.2.3,fish.1.2.4,fish.1.2.5,fish.1.2.6,fish.1.2.7,
          fish.1.3.2,fish.1.3.3,fish.1.3.4,fish.1.3.5,fish.1.3.6,fish.1.3.7,
          fish.1.4.2,fish.1.4.3,fish.1.4.4,fish.1.4.5,fish.1.4.6,fish.1.4.7,
          fish.1.5.2,fish.1.5.3,fish.1.5.4,fish.1.5.5,fish.1.5.6,fish.1.5.7,
          fish.1.6.2,fish.1.6.3,fish.1.6.4,fish.1.6.5,fish.1.6.6,fish.1.6.7,
          fish.1.7.2,fish.1.7.3,fish.1.7.4,fish.1.7.5,fish.1.7.6,fish.1.7.7,
          fish.1.2.4.1,fish.1.2.4.2,fish.1.2.4.3,fish.1.2.4.4,fish.1.2.4.5,
          fish.1.4.6.1,fish.1.4.6.2,fish.1.4.6.3,fish.1.4.6.4,fish.1.4.6.5,
          fish.1.4.3.1,fish.1.4.3.2,fish.1.4.3.3,fish.1.4.3.4,fish.1.4.3.5)[2]

delta<-aics-min(aics) 
delta[order(delta$AIC),, drop = FALSE]

fish.1.2.4.5.1<-update(fish.1.2.4.5,.~.- s(MAC, k = 4))
fish.1.2.4.5.2<-update(fish.1.2.4.5,.~.- s(ANE, k = 4))
fish.1.2.4.5.3<-update(fish.1.2.4.5,.~.- s(HOM, k = 4))
fish.1.2.4.5.4<-update(fish.1.2.4.5,.~.- s(HER, k = 4))

aics<-AIC(fish.1,fish.1.1,fish.1.2,fish.1.3,fish.1.4,fish.1.5,fish.1.6,fish.1.7,
          fish.1.1.2,fish.1.1.3,fish.1.1.4,fish.1.1.5,fish.1.1.6,fish.1.1.7,
          fish.1.2.2,fish.1.2.3,fish.1.2.4,fish.1.2.5,fish.1.2.6,fish.1.2.7,
          fish.1.3.2,fish.1.3.3,fish.1.3.4,fish.1.3.5,fish.1.3.6,fish.1.3.7,
          fish.1.4.2,fish.1.4.3,fish.1.4.4,fish.1.4.5,fish.1.4.6,fish.1.4.7,
          fish.1.5.2,fish.1.5.3,fish.1.5.4,fish.1.5.5,fish.1.5.6,fish.1.5.7,
          fish.1.6.2,fish.1.6.3,fish.1.6.4,fish.1.6.5,fish.1.6.6,fish.1.6.7,
          fish.1.7.2,fish.1.7.3,fish.1.7.4,fish.1.7.5,fish.1.7.6,fish.1.7.7,
          fish.1.2.4.1,fish.1.2.4.2,fish.1.2.4.3,fish.1.2.4.4,fish.1.2.4.5,
          fish.1.4.6.1,fish.1.4.6.2,fish.1.4.6.3,fish.1.4.6.4,fish.1.4.6.5,
          fish.1.4.3.1,fish.1.4.3.2,fish.1.4.3.3,fish.1.4.3.4,fish.1.4.3.5,
          fish.1.2.4.5.1,fish.1.2.4.5.2,fish.1.2.4.5.3,fish.1.2.4.5.4)[2]

delta<-aics-min(aics) 
delta[order(delta$AIC),, drop = FALSE]

aics<-AIC(env.4,env.4.5,env.5,env.5.5,env.1,
          fish.1.2.4,fish.1.4,fish.1.2,fish.1.2.4.5,fish.1)[2]
delta<-aics-min(aics)
delta[order(delta$AIC),, drop = FALSE]
round(delta[order(delta$AIC),, drop = FALSE],1)

