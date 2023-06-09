library("mgcv")

effort_data<- readRDS(file = "0.2 1km_explain_bl_sight_180s_wlatlon_in1km June3 2020.rds")
effort_data <- effort_data[complete.cases(effort_data), ]


mg.2 <-mgcv::gam(has_bird ~ s(max_cur_speed)+s(bathymetry)+s(aspect) +s(whole_sst) +s(whole_sal) +s(whole_chl)+s(cell_lat)+s(cell_lon),data = effort_data,family=binomial)
mg.3<-update(mg.2,.~.- s(aspect)+ s(anomalies)) 
mg.4<-update(mg.2,.~.- s(aspect)+ s(roughness)) 

mg.2.1<-update(mg.2,.~.- s(max_cur_speed)+ s(stratf_indx)) 
mg.3.1<-update(mg.3,.~.- s(max_cur_speed)+ s(stratf_indx))
mg.4.1<-update(mg.4,.~.- s(max_cur_speed)+ s(stratf_indx)) 

aics<-AIC(mg.2,mg.3,mg.4,mg.2.1,mg.3.1,mg.4.1)[2]
delta<-aics-min(aics) 
delta[order(delta$AIC),, drop = FALSE]


mg.3.1.1<-update(mg.3.1,.~.- s(bathymetry))
mg.3.1.2<-update(mg.3.1,.~.- s(whole_sst))
mg.3.1.3<-update(mg.3.1,.~.- s(whole_sal))
mg.3.1.4<-update(mg.3.1,.~.- s(whole_chl))
mg.3.1.5<-update(mg.3.1,.~.- s(cell_lat))
mg.3.1.6<-update(mg.3.1,.~.- s(cell_lon))
mg.3.1.7<-update(mg.3.1,.~.- s(anomalies))
mg.3.1.8<-update(mg.3.1,.~.- s(stratf_indx))

mg.2.1.1<-update(mg.2.1,.~.- s(bathymetry))
mg.2.1.2<-update(mg.2.1,.~.- s(whole_sst))
mg.2.1.3<-update(mg.2.1,.~.- s(whole_sal))
mg.2.1.4<-update(mg.2.1,.~.- s(whole_chl))
mg.2.1.5<-update(mg.2.1,.~.- s(cell_lat))
mg.2.1.6<-update(mg.2.1,.~.- s(cell_lon))
mg.2.1.7<-update(mg.2.1,.~.- s(aspect))
mg.2.1.8<-update(mg.2.1,.~.- s(stratf_indx))

mg.4.1.1<-update(mg.4.1,.~.- s(bathymetry))
mg.4.1.2<-update(mg.4.1,.~.- s(whole_sst))
mg.4.1.3<-update(mg.4.1,.~.- s(whole_sal))
mg.4.1.4<-update(mg.4.1,.~.- s(whole_chl))
mg.4.1.5<-update(mg.4.1,.~.- s(cell_lat))
mg.4.1.6<-update(mg.4.1,.~.- s(cell_lon))
mg.4.1.7<-update(mg.4.1,.~.- s(roughness))
mg.4.1.8<-update(mg.4.1,.~.- s(stratf_indx))


aics<-AIC(mg.2,mg.3,mg.4,mg.2.1,mg.3.1,mg.4.1,
          mg.3.1.1,mg.3.1.2,mg.3.1.3,mg.3.1.4,mg.3.1.5,mg.3.1.6,mg.3.1.7,mg.3.1.8,
          mg.2.1.1,mg.2.1.2,mg.2.1.3,mg.2.1.4,mg.2.1.5,mg.2.1.6,mg.2.1.7,mg.2.1.8,
          mg.4.1.1,mg.4.1.2,mg.4.1.3,mg.4.1.4,mg.4.1.5,mg.4.1.6,mg.4.1.7,mg.4.1.8)[2]
delta<-aics-min(aics) 
delta[order(delta$AIC),, drop = FALSE]
