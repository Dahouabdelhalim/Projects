# KINEMATIC VECTOR COMPARISONS & FIGURE S9

# 1)Construction of dataframe with kinematic, population, and diet variables 
# 
kinem.df.res.pairwise1<-data.frame(kinem.df.res[,c("Population","Watershed","Habitat","Spec","LD1","simpson",
                                                   "Op4.2d.KT.res","Max4.2d.KT.res","Mand.2d.openingLR.res","Mand.2d.closingLR.res")], "Dim" = as.vector(rep("2",times=59)))
kinem.df.res.pairwise2<-data.frame(kinem.df.res[,c("Population","Watershed","Habitat","Spec","LD1","simpson",
                                                   "Op4.2d.KT.15.res","Max4.2d.KT.10.res","Mand.2d.openingLR.res","Mand.2d.closingLR.res")], "Dim" = as.vector(rep("2",times=59)))
kinem.df.res.pairwise3<-data.frame(kinem.df.res[,c("Population","Watershed","Habitat","Spec","LD1","simpson",
                                                   "Op4.2d.KT.30.res","Max4.2d.KT.20.res","Mand.2d.openingLR.res","Mand.2d.closingLR.res")], "Dim" = as.vector(rep("2",times=59)))
kinem.df.res.pairwise4<-data.frame(kinem.df.res[,c("Population","Watershed","Habitat","Spec","LD1","simpson",
                                                   "Op4.3d.KT.res","Max4.3d.KT.res","Mand.3d.openingLR.res","Mand.3d.closingLR.res")], "Dim" = as.vector(rep("3",times=59)))
kinem.df.res.pairwise5<-data.frame(kinem.df.res[,c("Population","Watershed","Habitat","Spec","LD1","simpson",
                                                   "Op4.3d.KT.15.res","Max4.3d.KT.10.res","Mand.3d.openingLR.res","Mand.3d.closingLR.res")], "Dim" = as.vector(rep("3",times=59)))
kinem.df.res.pairwise6<-data.frame(kinem.df.res[,c("Population","Watershed","Habitat","Spec","LD1","simpson",
                                                   "Op4.3d.KT.30.res","Max4.3d.KT.20.res","Mand.3d.openingLR.res","Mand.3d.closingLR.res")], "Dim" = as.vector(rep("3",times=59)))
colnames(kinem.df.res.pairwise1)<-c("Population","Watershed","Habitat","Spec","LD1","simpson","Op4","Max4","Mand.openingLR","Mand.closingLR","Dim")
colnames(kinem.df.res.pairwise2)<-c("Population","Watershed","Habitat","Spec","LD1","simpson","Op4","Max4","Mand.openingLR","Mand.closingLR","Dim")
colnames(kinem.df.res.pairwise3)<-c("Population","Watershed","Habitat","Spec","LD1","simpson","Op4","Max4","Mand.openingLR","Mand.closingLR","Dim")
colnames(kinem.df.res.pairwise4)<-c("Population","Watershed","Habitat","Spec","LD1","simpson","Op4","Max4","Mand.openingLR","Mand.closingLR","Dim")
colnames(kinem.df.res.pairwise5)<-c("Population","Watershed","Habitat","Spec","LD1","simpson","Op4","Max4","Mand.openingLR","Mand.closingLR","Dim")
colnames(kinem.df.res.pairwise6)<-c("Population","Watershed","Habitat","Spec","LD1","simpson","Op4","Max4","Mand.openingLR","Mand.closingLR","Dim")
kinem.df.res.pairwise<-cbind(rbind(kinem.df.res.pairwise1,kinem.df.res.pairwise2,kinem.df.res.pairwise3,
                                   kinem.df.res.pairwise4,kinem.df.res.pairwise5,kinem.df.res.pairwise6),
                             "tilt"= c(rep("2D-0", times=59),rep("2D-mid", times=59),rep("2D-max", times=59),
                                       rep("3D-0", times=59),rep("3D-mid", times=59),rep("3D-max", times=59)))
KT.LR.2Dx3D_0PW<-cbind(kinem.df.res.pairwise$Op4,kinem.df.res.pairwise$Max4,
                    kinem.df.res.pairwise$Mand.openingLR,kinem.df.res.pairwise$Mand.closingLR)

kinem.df.res.pairwise.short<-cbind(rbind(kinem.df.res.pairwise1,kinem.df.res.pairwise4),
                             "tilt"= c(rep("2D-0", times=59),rep("3D-0", times=59)))
KT.LR.2Dx3D_0PW.short<-cbind(kinem.df.res.pairwise.short$Op4,kinem.df.res.pairwise.short$Max4,
                             kinem.df.res.pairwise.short$Mand.openingLR,kinem.df.res.pairwise.short$Mand.closingLR)

# 1.1) Table S16
# Watershed*Dimensions*LD
kinem.pairwise.lm.rrpp.short<-lm.rrpp(data = kinem.df.res.pairwise.short,KT.LR.2Dx3D_0PW.short~Watershed*tilt*LD1, iter = 9999,SS.type = "II")
kinem.PW.WatxtiltxLD.lm.short<-pairwise(kinem.pairwise.lm.rrpp.short, groups = interaction(kinem.df.res.pairwise.short$Watershed,kinem.df.res.pairwise.short$tilt),
                                  covariate = kinem.df.res.pairwise.short$LD1)
anova.lm.rrpp(kinem.pairwise.lm.rrpp.short)
WatxtiltLD.angle.short<-summary(kinem.PW.WatxtiltxLD.lm.short, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
WatxtiltLD.angleP.short<-summary(kinem.PW.WatxtiltxLD.lm.short, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
WatxtiltLD.dist.short<-summary(kinem.PW.WatxtiltxLD.lm.short, confidence = 0.95, test.type = "dist")$pairwise.tables$D
WatxtiltLD.distP.short<-summary(kinem.PW.WatxtiltxLD.lm.short, confidence = 0.95, test.type = "dist")$pairwise.tables$P
# Dimensions*LD
kinem.tot.pairwise.lm.rrpp.short<-lm.rrpp(data = kinem.df.res.pairwise.short,KT.LR.2Dx3D_0PW.short~tilt*LD1, iter = 9999,SS.type = "II")
kinem.PW.tot.tiltxLD.lm.short<-pairwise(kinem.tot.pairwise.lm.rrpp.short, groups = kinem.df.res.pairwise.short$tilt,
                                        covariate = kinem.df.res.pairwise.short$LD1)
tot.tiltxLD.angle.short<-summary(kinem.PW.tot.tiltxLD.lm.short, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
tot.tiltxLD.angleP.short<-summary(kinem.PW.tot.tiltxLD.lm.short, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
tot.tiltxLD.dist.short<-summary(kinem.PW.tot.tiltxLD.lm.short, confidence = 0.95, test.type = "dist")$pairwise.tables$D
tot.tiltxLD.distP.short<-summary(kinem.PW.tot.tiltxLD.lm.short, confidence = 0.95, test.type = "dist")$pairwise.tables$P
# Watershed*Dimensions*Simpson
kinem.pairwise.simpson.lm.rrpp.short<-lm.rrpp(data = kinem.df.res.pairwise.short,KT.LR.2Dx3D_0PW.short~Watershed*tilt*simpson, iter = 9999,SS.type = "II")
kinem.PW.Watxtiltxsimpson.lm.short<-pairwise(kinem.pairwise.simpson.lm.rrpp.short, groups = interaction(kinem.df.res.pairwise.short$Watershed,kinem.df.res.pairwise.short$tilt),
                                        covariate = kinem.df.res.pairwise.short$simpson)
Watxtiltxsimpson.angle.short<-summary(kinem.PW.Watxtiltxsimpson.lm.short, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
Watxtiltxsimpson.angleP.short<-summary(kinem.PW.Watxtiltxsimpson.lm.short, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
Watxtiltxsimpson.dist.short<-summary(kinem.PW.Watxtiltxsimpson.lm.short, confidence = 0.95, test.type = "dist")$pairwise.tables$D
Watxtiltxsimpson.distP.short<-summary(kinem.PW.Watxtiltxsimpson.lm.short, confidence = 0.95, test.type = "dist")$pairwise.tables$P
# Dimensions*Simpson
kinem.tot.pairwise.simpson.lm.rrpp.short<-lm.rrpp(data = kinem.df.res.pairwise.short,KT.LR.2Dx3D_0PW.short~tilt*simpson, iter = 9999,SS.type = "II")
kinem.PW.tot.tiltxsimpson.lm.short<-pairwise(kinem.tot.pairwise.simpson.lm.rrpp.short, groups = kinem.df.res.pairwise.short$tilt,
                                        covariate = kinem.df.res.pairwise.short$LD1)
tot.tiltxsimpson.angle.short<-summary(kinem.PW.tot.tiltxsimpson.lm.short, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
tot.tiltxsimpson.angleP.short<-summary(kinem.PW.tot.tiltxsimpson.lm.short, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
tot.tiltxsimpson.dist.short<-summary(kinem.PW.tot.tiltxsimpson.lm.short, confidence = 0.95, test.type = "dist")$pairwise.tables$D
tot.tiltxsimpson.distP.short<-summary(kinem.PW.tot.tiltxsimpson.lm.short, confidence = 0.95, test.type = "dist")$pairwise.tables$P




# 2.0) Watershed * kinematic model versions * dietLD vector comparisons
kinem.pairwise.lm.rrpp<-lm.rrpp(data = kinem.df.res.pairwise,KT.LR.2Dx3D_0PW~Watershed*tilt*LD1, iter = 9999,SS.type = "II")
kinem.PW.WatxtiltxLD.lm<-pairwise(kinem.pairwise.lm.rrpp, groups = interaction(kinem.df.res.pairwise$Watershed,kinem.df.res.pairwise$tilt),
                                  covariate = kinem.df.res.pairwise$LD1)
anova.lm.rrpp(kinem.pairwise.lm.rrpp)
WatxtiltLD.angle<-summary(kinem.PW.WatxtiltxLD.lm, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
WatxtiltLD.angleP<-summary(kinem.PW.WatxtiltxLD.lm, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
WatxtiltLD.dist<-summary(kinem.PW.WatxtiltxLD.lm, confidence = 0.95, test.type = "dist")$pairwise.tables$D
WatxtiltLD.distP<-summary(kinem.PW.WatxtiltxLD.lm, confidence = 0.95, test.type = "dist")$pairwise.tables$P

# 2.1) Watershed * kinematic model versions * simpson vector comparisons
kinem.pairwise.simpson.lm.rrpp<-lm.rrpp(data = kinem.df.res.pairwise,KT.LR.2Dx3D_0PW~Watershed*tilt*simpson, iter = 9999,SS.type = "II")
kinem.PW.Watxtiltxsimp.lm<-pairwise(kinem.pairwise.simpson.lm.rrpp, groups = interaction(kinem.df.res.pairwise$Watershed,kinem.df.res.pairwise$tilt),
                                    covariate = kinem.df.res.pairwise$simpson)
anova.lm.rrpp(kinem.pairwise.simpson.lm.rrpp)
Watxtiltxsimp.angle<-summary(kinem.PW.Watxtiltxsimp.lm, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
Watxtiltxsimp.angleP<-summary(kinem.PW.Watxtiltxsimp.lm, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
Watxtiltxsimp.dist<-summary(kinem.PW.Watxtiltxsimp.lm, confidence = 0.95, test.type = "dist")$pairwise.tables$D
Watxtiltxsimp.distP<-summary(kinem.PW.Watxtiltxsimp.lm, confidence = 0.95, test.type = "dist")$pairwise.tables$P

# 2.2)kinematic model versions * dietLD vector comparisons (NO WATERSHED FACTOR)
kinem.pairwise.totLD.lm.rrpp<-lm.rrpp(data = kinem.df.res.pairwise,KT.LR.2Dx3D_0PW~tilt*LD1, iter = 9999,SS.type = "II")
kinem.PW.totLD.lm<-pairwise(kinem.pairwise.totLD.lm.rrpp, groups = kinem.df.res.pairwise$tilt,covariate = kinem.df.res.pairwise$LD1)
anova.lm.rrpp(kinem.pairwise.totLD.lm.rrpp)
totLD.angle<-summary(kinem.PW.totLD.lm, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
totLD.angleP<-summary(kinem.PW.totLD.lm, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
totLD.dist<-summary(kinem.PW.totLD.lm, confidence = 0.95, test.type = "dist")$pairwise.tables$D
totLD.distP<-summary(kinem.PW.totLD.lm, confidence = 0.95, test.type = "dist")$pairwise.tables$P

# 2.3)kinematic model versions * simpson vector comparisons (NO WATERSHED FACTOR)
kinem.pairwise.totsimpson.lm.rrpp<-lm.rrpp(data = kinem.df.res.pairwise,KT.LR.2Dx3D_0PW~tilt*simpson, iter = 9999,SS.type = "II")
kinem.PW.totsimpson.lm<-pairwise(kinem.pairwise.totsimpson.lm.rrpp, groups = kinem.df.res.pairwise$tilt,covariate = kinem.df.res.pairwise$simpson)
anova.lm.rrpp(kinem.pairwise.totsimpson.lm.rrpp)
totsimpson.angle<-summary(kinem.PW.totsimpson.lm, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
totsimpson.angleP<-summary(kinem.PW.totsimpson.lm, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
totsimpson.dist<-summary(kinem.PW.totsimpson.lm, confidence = 0.95, test.type = "dist")$pairwise.tables$D
totsimpson.distP<-summary(kinem.PW.totsimpson.lm, confidence = 0.95, test.type = "dist")$pairwise.tables$P



#
#
#
#
#
#
###    PLOTS


## 3.1) dataframes of factors and kinematic PC scores for kin. models w/ planar rotation, 
# and with maximum tilt of rotational axis in Op4 and Max4.
biomech.2d<-rda(kinem.df.res[c("Op4.2d.KT.res","Max4.2d.KT.res","Mand.2d.openingLR.res", "Mand.2d.closingLR.res")])
biomech.2d.df<-data.frame(summary(biomech.2d)["sites"],
                          kinem.df.res[c("Population","Watershed","Habitat","LD1","simpson")])
biomech.2d.maxtilt<-rda(kinem.df.res[c("Op4.2d.KT.30.res","Max4.2d.KT.20.res","Mand.2d.openingLR.res", "Mand.2d.closingLR.res")])
biomech.2d.maxtilt.df<-data.frame(summary(biomech.2d.maxtilt)["sites"],
                                  kinem.df.res[c("Population","Watershed","Habitat","LD1","simpson")])
biomech.3d<-rda(kinem.df.res[c("Op4.3d.KT.res","Max4.3d.KT.res", "Mand.3d.openingLR.res", "Mand.3d.closingLR.res")])
biomech.3d.df<-data.frame(summary(biomech.3d)["sites"], 
                          kinem.df.res[c("Population","Watershed","Habitat","LD1","simpson")])
biomech.3d.maxtilt<-rda(kinem.df.res[c("Op4.3d.KT.30.res","Max4.3d.KT.20.res", "Mand.3d.openingLR.res", "Mand.3d.closingLR.res")])
biomech.3d.maxtilt.df<-data.frame(summary(biomech.3d.maxtilt)["sites"],
                                  kinem.df.res[c("Population","Watershed","Habitat","LD1","simpson")])

## 3.2) Kinematic PC score means by population
biomech.2d.PCmeans<-cbind( biomech.2d.df %>% group_by(Population) %>%dplyr::summarize(PC1 = mean(sites.PC1)),
                           PC2=as.data.frame(biomech.2d.df %>% group_by(Population) %>%dplyr::summarize(PC2 = mean(sites.PC2)))[,"PC2"])
biomech.2d.PCmeans<-distinct(base::merge(biomech.2d.df[c("Habitat", "Watershed", "Population","LD1","simpson")],
                                         biomech.2d.PCmeans, by = "Population", all.x = F))

biomech.3d.PCmeans<-cbind( biomech.3d.df %>% group_by(Population) %>%dplyr::summarize(PC1 = mean(sites.PC1)),
                           PC2=as.data.frame(biomech.3d.df %>% group_by(Population) %>%dplyr::summarize(PC2 = mean(sites.PC2)))[,"PC2"])
biomech.3d.PCmeans<-distinct(base::merge(biomech.3d.df[c("Habitat", "Watershed", "Population","LD1","simpson")],
                                         biomech.3d.PCmeans, by = "Population", all.x = F))

biomech.3d.maxtilt.PCmeans<-cbind( biomech.3d.maxtilt.df %>% group_by(Population) %>%dplyr::summarize(PC1 = mean(sites.PC1)),
                                   PC2=as.data.frame(biomech.3d.maxtilt.df %>% group_by(Population) %>%dplyr::summarize(PC2 = mean(sites.PC2)))[,"PC2"])
biomech.3d.maxtilt.PCmeans<-distinct(base::merge(biomech.3d.maxtilt.df[c("Habitat", "Watershed", "Population","LD1","simpson")],
                                                 biomech.3d.maxtilt.PCmeans, by = "Population", all.x = F))


## 3) 2D mechanical function Plots
#  3.1) 2D version of Figure 5A
biomech.functions.2dxLD<- ggplot(data = biomech.2d.df, 
                                 aes(x=LD1, y=sites.PC1))+geom_point(alpha = .3, stroke = 2)+
  xlab("LD1")+ylab("PC1") #+ #geom_density2d(alpha=.10) +
#stat_ellipse(type = "norm", level = .95)
biomech.functions.2dxLD<- biomech.functions.2dxLD+ aes(colour=Habitat, shape=Watershed, stroke = 1)+ 
  scale_shape_manual(values = c(17,16,15,3,4,18)) + xlab("LD") +
  scale_color_manual(values = c("slategray2","red4"))
biomech.functions.2dxLD<-biomech.functions.2dxLD +  
  geom_line(data = biomech.2d.PCmeans, aes(x = LD1, y = PC1, group = Watershed), size = .5, color = "black") +
  geom_point(data = biomech.2d.PCmeans, aes(x = LD1, y = PC1), size = 3, stroke = 2.5)  
biomech.functions.2dxLD

#  3.2) 2D version of Figure 5B
biomech.functions.PC1.2dxsimpson<- ggplot(data = biomech.2d.df, 
                                          aes(x=simpson, y=sites.PC1))+geom_point(alpha = .3, stroke = 2)+
  xlab("simpson")+ylab("PC1") #+ #geom_density2d(alpha=.10) +
#stat_ellipse(type = "norm", level = .95)
biomech.functions.PC1.2dxsimpson<- biomech.functions.PC1.2dxsimpson+ aes(colour=Habitat, shape=Watershed, stroke = 1)+ 
  scale_shape_manual(values = c(17,16,15,3,4,18)) + xlab("simpson") +
  scale_color_manual(values = c("slategray2","red4"))
biomech.functions.PC1.2dxsimpson<-biomech.functions.PC1.2dxsimpson +  
  geom_line(data = biomech.2d.PCmeans, aes(x = simpson, y = PC1, group = Watershed), size = .5, color = "black") +
  geom_point(data = biomech.2d.PCmeans, aes(x = simpson, y = PC1), size = 3, stroke = 2.5)  
biomech.functions.PC1.2dxsimpson


## 4) 3D mechanical function Plots
#  4.1) No tilt on Op4 or Max4 version of Figure 5A
biomech.functions.3dxLD<- ggplot(data = biomech.3d.df, 
                                 aes(x=LD1, y=sites.PC1))+geom_point(alpha = .3, stroke = 2)+
  xlab("LD1")+ylab("PC1") #+ #geom_density2d(alpha=.10) +
#stat_ellipse(type = "norm", level = .95)
biomech.functions.3dxLD<- biomech.functions.3dxLD+ aes(colour=Habitat, shape=Watershed, stroke = 1)+ 
  scale_shape_manual(values = c(17,16,15,3,4,18)) + xlab("diet LD") +
  scale_color_manual(values = c("slategray2","red4"))
biomech.functions.3dxLD<-biomech.functions.3dxLD +  
  geom_line(data = biomech.3d.PCmeans, aes(x = LD1, y = PC1, group = Watershed), size = .5, color = "black") +
  geom_point(data = biomech.3d.PCmeans, aes(x = LD1, y = PC1), size = 3, stroke = 2.5)  
biomech.functions.3dxLD

#  4.2) No tilt on Op4 or Max4 version of Figure 5B
biomech.functions.PC1.3dxsimpson<- ggplot(data = biomech.3d.df, 
                                          aes(x=simpson, y=sites.PC1))+geom_point(alpha = .3, stroke = 2)+
  xlab("simpson")+ylab("PC1") #+ #geom_density2d(alpha=.10) +
#stat_ellipse(type = "norm", level = .95)
biomech.functions.PC1.3dxsimpson<- biomech.functions.PC1.3dxsimpson+ aes(colour=Habitat, shape=Watershed, stroke = 1)+ 
  scale_shape_manual(values = c(17,16,15,3,4,18)) + xlab("simpson") +
  scale_color_manual(values = c("slategray2","red4"))
biomech.functions.PC1.3dxsimpson<-biomech.functions.PC1.3dxsimpson +  
  geom_line(data = biomech.3d.PCmeans, aes(x = simpson, y = PC1, group = Watershed), size = .5, color = "black") +
  geom_point(data = biomech.3d.PCmeans, aes(x = simpson, y = PC1), size = 3, stroke = 2.5)  
biomech.functions.PC1.3dxsimpson


## 5) 3D tilted axis mechanical function Plots
#  5.1) Figure 5A
biomech.functions.3d.maxtiltxLD<- ggplot(data = biomech.3d.maxtilt.df, 
                                         aes(x=LD1, y=sites.PC1))+geom_point(alpha = .3, stroke = 2)+
  xlab("Diet LD")+ylab("Kin. PC1 (94.6%)") #+ #geom_density2d(alpha=.10) +
#stat_ellipse(type = "norm", level = .95)
biomech.functions.3d.maxtiltxLD<- biomech.functions.3d.maxtiltxLD+ aes(colour=Habitat, shape=Watershed, stroke = 1)+ 
  scale_shape_manual(values = c(17,16,15,3,4,18)) + xlab("Diet LD") +
  scale_color_manual(values = c("slategray2","red4"))
biomech.functions.3d.maxtiltxLD<-biomech.functions.3d.maxtiltxLD +  
  geom_line(data = biomech.3d.maxtilt.PCmeans, aes(x = LD1, y = PC1, group = Watershed), size = .5, color = "black") +
  geom_point(data = biomech.3d.maxtilt.PCmeans, aes(x = LD1, y = PC1), size = 3, stroke = 2.5)  
biomech.functions.3d.maxtiltxLD
#  5.2) Figure 5B
biomech.functions.PC1.3d.maxtiltxsimpson<- ggplot(data = biomech.3d.maxtilt.df, 
                                                  aes(x=simpson, y=sites.PC1))+geom_point(alpha = .3, stroke = 2)+
  xlab("Gini-Simpson Index")+ylab("")  #+ #geom_density2d(alpha=.10) +
#stat_ellipse(type = "norm", level = .95)
biomech.functions.PC1.3d.maxtiltxsimpson<- biomech.functions.PC1.3d.maxtiltxsimpson+ aes(colour=Habitat, shape=Watershed, stroke = 1)+ 
  scale_shape_manual(values = c(17,16,15,3,4,18)) + xlab("Gini-Simpson Index") +
  scale_color_manual(values = c("slategray2","red4"))
biomech.functions.PC1.3d.maxtiltxsimpson<-biomech.functions.PC1.3d.maxtiltxsimpson +  
  geom_line(data = biomech.3d.maxtilt.PCmeans, aes(x = simpson, y = PC1, group = Watershed), size = .5, color = "black") +
  geom_point(data = biomech.3d.maxtilt.PCmeans, aes(x = simpson, y = PC1), size = 3, stroke = 2.5) 
biomech.functions.PC1.3d.maxtiltxsimpson

# 6) combining Figure 5A and 5B to form Figure 5
ggarrange(biomech.functions.3d.maxtiltxLD, biomech.functions.PC1.3d.maxtiltxsimpson,
          labels = c("A", "B"),ncol = 2, nrow = 1, common.legend = T, legend = "right")



#  7) PCAs for figure S9
Op4.3d.kinem.df<-data.frame(plotTangentSpace(Op4.3d.sym)$pc.scores[,1:4], kinem.df[,-c(1:2)],kinem.df.res[c("LD1","Habitat", "Watershed", "Population")])
Max4.3d.kinem.df<-data.frame(plotTangentSpace(Max4.3d.sym)$pc.scores[,1:4], kinem.df[,-c(1:2)],kinem.df.res[c("LD1","Habitat", "Watershed", "Population")])
Mand.3d.kinem.df<-data.frame(plotTangentSpace(Mand.3d.sym)$pc.scores[,1:4], kinem.df[,-c(1:2)],kinem.df.res[c("LD1","Habitat", "Watershed", "Population")])
tot.3d.kinem.df<-data.frame(plotTangentSpace(spec.sym.3d.gpa$coords)$pc.scores[,1:4], biomech.3d.df[,1:2], kinem.df[,-c(1:2)],kinem.df.res[c("LD1","Habitat", "Watershed", "Population")])
tot.3d.kinem.maxtilt.df<-data.frame(plotTangentSpace(spec.sym.3d.gpa$coords)$pc.scores[,1:4], biomech.3d.maxtilt.df[,1:2], kinem.df[,-c(1:2)],kinem.df.res[c("LD1","Habitat", "Watershed", "Population")])
tot.3d.kinem.maxtilt.df<-cbind(tot.3d.kinem.maxtilt.df,maxtilt3D_2d.dispPC1=summary(rda(cbind(kinem.df.res$Op4.3d.KT.30.res-kinem.df.res$Op4.2d.KT.res,
                                                                                              kinem.df.res$Max4.3d.KT.20.res-kinem.df.res$Max4.2d.KT.res, 
                                                                                              kinem.df.res$Mand.3d.openingLR.res-kinem.df.res$Mand.2d.openingLR.res, 
                                                                                              kinem.df.res$Mand.3d.closingLR.res-kinem.df.res$Mand.2d.closingLR.res)))$sites[,1])
## PC1 explains 89.9% of variance in maxtilt3D_2d.disp
tot.3d.kinem.maxtilt.df<-cbind(tot.3d.kinem.maxtilt.df,maxtilt.dispPC1=summary(rda(cbind(kinem.df.res$Op4.3d.KT.30.res-kinem.df.res$Op4.2d.KT.30.res,
                                                                                         kinem.df.res$Max4.3d.KT.20.res-kinem.df.res$Max4.2d.KT.20.res, 
                                                                                         kinem.df.res$Mand.3d.openingLR.res-kinem.df.res$Mand.2d.openingLR.res, 
                                                                                         kinem.df.res$Mand.3d.closingLR.res-kinem.df.res$Mand.2d.closingLR.res)))$sites[,1])
## PC1 explains 67.3% of variance in maxtilt.disp

# 7.1) Figure S9 plots
ggarrange(
  ggplot(data = Op4.3d.kinem.df, aes(x=PC1, y=PC2, colour = Op4.3d.KT.30, stroke = 3, shape = Habitat))+geom_point()+scale_colour_gradient2(low = "red",mid = "white", high = "steelblue", midpoint = mean(Op4.2d.KT))+
    xlab("PC1 - 63.6% ")+ylab("PC2 - 10.2%")+labs(colour = "Op4 KT")+
    stat_ellipse(aes(group=Habitat, linetype = Habitat), level = 0.8, type = c("norm")),
  ggplot(data = Max4.3d.kinem.df, aes(x=PC1, y=PC2, colour = Max4.3d.KT.20, stroke = 3, shape = Habitat))+geom_point()+scale_colour_gradient2(low = "red",mid = "white", high = "steelblue", midpoint = mean(Max4.2d.KT))+
    xlab("PC1 - 34.3% ")+ylab("PC2 - 21.4%")+labs(colour = "Max4 KT")+
    stat_ellipse(aes(group=Habitat, linetype = Habitat), level = 0.8, type = c("norm")),
  ggplot(data = Mand.3d.kinem.df, aes(x=PC1, y=PC2, colour = Mand.3d.closingLR, stroke = 3, shape = Habitat))+geom_point()+scale_colour_gradient2(low = "red",mid = "white", high = "steelblue", midpoint = mean(Mand.2d.closingLR))+
    xlab("PC1 - 39.3% ")+ylab("PC2 - 27.8%")+ labs(colour = "Closing LR")+
    stat_ellipse(aes(group=Habitat, linetype = Habitat), level = 0.8, type = c("norm")),
  ggplot(data = Mand.3d.kinem.df, aes(x=PC1, y=PC2, colour = Mand.3d.openingLR, stroke = 3, shape = Habitat))+geom_point()+scale_colour_gradient2(low = "red",mid = "white", high = "steelblue", midpoint = mean(Mand.2d.openingLR))+
    xlab("PC1 - 39.3% ")+ylab("PC2 - 27.8%") + labs(colour = "Opening LR")+
    stat_ellipse(aes(group=Habitat, linetype = Habitat), level = 0.8, type = c("norm")),
  ggplot(data = tot.3d.kinem.df, aes(x=PC1, y=PC2, colour = tot.3d.kinem.df$sites.PC1, stroke = 3, shape = Habitat))+geom_point()+scale_colour_gradient2(low = "red",mid = "white", high = "steelblue", midpoint = mean(tot.2d.kinem.df$sites.PC1))+
    xlab("PC1 - 48.6% ")+ylab("PC2 - 12.9%") + labs(colour = "Kinematic PC1")+
    stat_ellipse(aes(group=Habitat, linetype = Habitat), level = 0.8, type = c("norm")),
  ggplot(data = Op4.3d.kinem.df, aes(x=PC1, y=PC2, colour = Op4.3d.KT.30-Op4.2d.KT, stroke = 3, shape = Habitat))+
    geom_point()+scale_colour_gradient2(low = "red",mid = "white", high = "steelblue", midpoint = mean(Op4.3d.KT.30-Op4.2d.KT))+
    xlab("PC1 - 63.6% ")+ylab("PC2 - 10.2%")+ labs(colour = "KT Disparity") +
    stat_ellipse(aes(group=Habitat, linetype = Habitat), level = 0.8, type = c("norm")),
  ggplot(data = Max4.3d.kinem.df, aes(x=PC1, y=PC2, colour = Max4.3d.KT.20-Max4.2d.KT, stroke = 3, shape = Habitat))+
    geom_point()+scale_colour_gradient2(low = "red",mid = "white", high = "steelblue", midpoint = mean(Max4.3d.KT.20-Max4.2d.KT))+
    xlab("PC1 - 34.3% ")+ylab("PC2 - 21.4%") +labs(colour = "KT Disparity") +
    stat_ellipse(aes(group=Habitat, linetype = Habitat), level = 0.8, type = "norm"),
  ggplot(data = Mand.3d.kinem.df, aes(x=PC1, y=PC2, colour = Mand.3d.closingLR-Mand.2d.closingLR, stroke = 3, shape = Habitat))+
    geom_point()+scale_colour_gradient2(low = "red",mid = "white", high = "steelblue", midpoint = mean(Mand.3d.closingLR-Mand.2d.closingLR))+
    xlab("PC1 - 39.3% ")+ylab("PC2 - 27.5%")+ labs(colour = "LR Disparity") +
    stat_ellipse(aes(group=Habitat, linetype = Habitat), level = 0.8, type = "norm") ,
  ggplot(data = Mand.3d.kinem.df, aes(x=PC1, y=PC2, colour = Mand.3d.openingLR-Mand.2d.openingLR, stroke = 3, shape = Habitat))+
    geom_point()+scale_colour_gradient2(low = "red",mid = "white", high = "steelblue", midpoint = mean(Mand.3d.openingLR-Mand.2d.openingLR))+
    xlab("PC1 - 39.3% ")+ylab("PC2 - 27.5%")+ labs(colour = "LR Disparity") +
    stat_ellipse(aes(group=Habitat,linetype = Habitat), level = 0.8, type = "norm") ,
  ggplot(data = tot.3d.kinem.maxtilt.df, aes(x=PC1, y=PC2, colour = maxtilt3D_2d.dispPC1, stroke = 3, shape = Habitat))+
    geom_point()+scale_colour_gradient2(low = "red",mid = "white", high = "steelblue", midpoint = mean(tot.3d.kinem.maxtilt.df$maxtilt3D_2d.dispPC1))+
    xlab("PC1 - 48.6% ")+ylab("PC2 - 12.9%")+ labs(colour = "Mech. PC1 Disp.") +
    stat_ellipse(aes(group=Habitat, linetype = Habitat), level = 0.8, type = "norm"),
  ncol =5 , nrow =2 , labels = c("A","B","C","D","E",  "F","G","H","I","J"))

# 8)
# wireframe function for figure S9 PC references
plot.coords <- function(A, W, points.col="black", points.cex=1, lines.col="black", lines.wd=2, bg.col=NULL, 
                        main=NULL, main.line=2, main.cex=2, legend=NULL, legend.pos="topright", legend.title="", 
                        legend.col=NULL, legend.cex=1.2, legend.lwd=2, legend.bty="n", params=NULL, add=FALSE) {
  if (!is.null(params)) {par3d(params)}
  points.col <- rep(points.col, length.out=nrow(A))
  points.cex <- rep(points.cex, length.out=nrow(A))
  lines.col <- rep(lines.col, length.out=nrow(W))
  lines.wd <- rep(lines.wd, length.out=nrow(W))
  if (!is.null(bg.col)) rgl.bg(sphere=TRUE, color=bg.col, lit=FALSE, back="fill")
  plot3d(A, type="s", col=points.col, xlab="", ylab="", zlab="", size=points.cex, aspect=FALSE, box=F, axes=F, add=add)
  if (!is.null(main) | !is.null(legend)) {
    if (!is.null(legend) & is.null(legend.col)) stop("must supply legend colors")
    bgplot3d({plot.new()
      if (!is.null(main)) title(main=main, line=main.line, cex.main=main.cex)
      if (!is.null(legend)) legend(legend.pos, title=legend.title, legend=legend, col=legend.col, lwd=legend.lwd, cex=legend.cex, bty=legend.bty)})}
  for (i in 1:nrow(W)) {
    segments3d(rbind(A[W[i,1],], A[W[i,2],]), lwd=lines.wd[i], col=lines.col[i])
  }
}

Op4.2d.sym<-spec.sym.2d.gpa$coords[1:4,,]
Max4.2d.sym<-spec.sym.2d.gpa$coords[c(4:5,7:8),,]
Mand.2d.sym<-spec.sym.2d.gpa$coords[c(3:6),,]
# 8.2) PCAs, with min and max coords for reference wireframes
Op4.3d.pc1min<-plotTangentSpace(spec.sym.3d.gpa$coords[c(1:4,9:12),,])$pc.shapes$PC1min
rownames(Op4.3d.pc1min)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")
Op4.3d.pc1max<-plotTangentSpace(spec.sym.3d.gpa$coords[c(1:4,9:12),,])$pc.shapes$PC1max
rownames(Op4.3d.pc1max)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")
Op4.3d.pc2min<-plotTangentSpace(spec.sym.3d.gpa$coords[c(1:4,9:12),,])$pc.shapes$PC2min
rownames(Op4.3d.pc2min)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")
Op4.3d.pc2max<-plotTangentSpace(spec.sym.3d.gpa$coords[c(1:4,9:12),,])$pc.shapes$PC2max
rownames(Op4.3d.pc2max)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")
Max4.3d.pc1min<-plotTangentSpace(spec.sym.3d.gpa$coords[c(4:5,7:8,12:13,15:16),,])$pc.shapes$PC1min
rownames(Max4.3d.pc1min)<-c("11r","6r","3r","10r","11l","6l","3l","10l")
Max4.3d.pc1max<-plotTangentSpace(spec.sym.3d.gpa$coords[c(4:5,7:8,12:13,15:16),,])$pc.shapes$PC1max
rownames(Max4.3d.pc1max)<-c("11r","6r","3r","10r","11l","6l","3l","10l")
Max4.3d.pc2min<-plotTangentSpace(spec.sym.3d.gpa$coords[c(4:5,7:8,12:13,15:16),,])$pc.shapes$PC2min
rownames(Max4.3d.pc2min)<-c("11r","6r","3r","10r","11l","6l","3l","10l")
Max4.3d.pc2max<-plotTangentSpace(spec.sym.3d.gpa$coords[c(4:5,7:8,12:13,15:16),,])$pc.shapes$PC2max
rownames(Max4.3d.pc2max)<-c("11r","6r","3r","10r","11l","6l","3l","10l")
Mand.3d.pc1min<-plotTangentSpace(spec.sym.3d.gpa$coords[c(3:6,11:14),,])$pc.shapes$PC1min
rownames(Mand.3d.pc1min)<-c("12r","11r","6r","2r","12l","11l","6l","2l")
Mand.3d.pc1max<-plotTangentSpace(spec.sym.3d.gpa$coords[c(3:6,11:14),,])$pc.shapes$PC1max
rownames(Mand.3d.pc1max)<-c("12r","11r","6r","2r","12l","11l","6l","2l")
Mand.3d.pc2min<-plotTangentSpace(spec.sym.3d.gpa$coords[c(3:6,11:14),,])$pc.shapes$PC2min
rownames(Mand.3d.pc2min)<-c("12r","11r","6r","2r","12l","11l","6l","2l")
Mand.3d.pc2max<-plotTangentSpace(spec.sym.3d.gpa$coords[c(3:6,11:14),,])$pc.shapes$PC2max
rownames(Mand.3d.pc2max)<-c("12r","11r","6r","2r","12l","11l","6l","2l")
sym.3d.pc1min<-plotTangentSpace(spec.sym.3d.gpa$coords)$pc.shapes$PC1min
rownames(sym.3d.pc1min)<-c("vOpr","Opr","12r","11r","6r","2r","3r","10r","vOpl","Opl","12l","11l","6l","2l","3l","10l")
sym.3d.pc1max<-plotTangentSpace(spec.sym.3d.gpa$coords)$pc.shapes$PC1max
rownames(sym.3d.pc1max)<-c("vOpr","Opr","12r","11r","6r","2r","3r","10r","vOpl","Opl","12l","11l","6l","2l","3l","10l")
sym.3d.pc2min<-plotTangentSpace(spec.sym.3d.gpa$coords)$pc.shapes$PC2min
rownames(sym.3d.pc2min)<-c("vOpr","Opr","12r","11r","6r","2r","3r","10r","vOpl","Opl","12l","11l","6l","2l","3l","10l")
sym.3d.pc2max<-plotTangentSpace(spec.sym.3d.gpa$coords)$pc.shapes$PC2max
rownames(sym.3d.pc2max)<-c("vOpr","Opr","12r","11r","6r","2r","3r","10r","vOpl","Opl","12l","11l","6l","2l","3l","10l")


W.sym.tot<-cbind(c("vOpr","vOpr","Opr","12r","11r","11r","11r","6r","10r","vOpl","vOpl","Opl","12l","11l","11l","11l","6l","10l","3r","2r","11r"),
                 c("Opr","12r","11r","11r","2r","6r","10r","3r","3r","Opl","12l","11l","11l","2l","6l","10l","3l","3l","3l","2l","11l"))
W.Op4<-rbind(W.sym.tot[1:4,],W.sym.tot[10:13,])
W.Max4<-rbind(W.sym.tot[6:9,], W.sym.tot[15:18,])
W.Dent<-rbind(W.sym.tot[4:6,], W.sym.tot[13:15,])

# 8.3) calibrate view
plot.coords(sym.3d.pc1min, W.sym.tot, points.col="black", lines.col="black", add = T)
view1<-par3d() # these 3 lines set views of wireframes
view2<-par3d()
view3<-par3d()

# 8.4) Lever system PC axis wireframes for Fig. S9 (PC1 row 1, PC2 row 2)
layout3d(matrix(1:12, nrow = 2), sharedMouse = T) #establishes grid for subsequent wireframes for lever systems
plot.coords(Op4.3d.pc1min, W.Op4, points.col="slategray4", lines.col="slategray4")
plot.coords(Op4.3d.pc1max, W.Op4, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view1$userMatrix)        
plot.coords(Op4.3d.pc2min, W.Op4, points.col="slategray4", lines.col="slategray4")
plot.coords(Op4.3d.pc2max, W.Op4, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view1$userMatrix) 
plot.coords(Op4.3d.pc1min, W.Op4, points.col="slategray4", lines.col="slategray4")
plot.coords(Op4.3d.pc1max, W.Op4, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view2$userMatrix)        
plot.coords(Op4.3d.pc2min, W.Op4, points.col="slategray4", lines.col="slategray4")
plot.coords(Op4.3d.pc2max, W.Op4, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view2$userMatrix)
plot.coords(Max4.3d.pc1min, W.Max4, points.col="slategray4", lines.col="slategray4")
plot.coords(Max4.3d.pc1max, W.Max4, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view1$userMatrix)        
plot.coords(Max4.3d.pc2min, W.Max4, points.col="slategray4", lines.col="slategray4")
plot.coords(Max4.3d.pc2max, W.Max4, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view1$userMatrix) 
plot.coords(Max4.3d.pc1min, W.Max4, points.col="slategray4", lines.col="slategray4")
plot.coords(Max4.3d.pc1max, W.Max4, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view2$userMatrix)        
plot.coords(Max4.3d.pc2min, W.Max4, points.col="slategray4", lines.col="slategray4")
plot.coords(Max4.3d.pc2max, W.Max4, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view2$userMatrix)
plot.coords(Mand.3d.pc1min, W.Dent, points.col="slategray4", lines.col="slategray4")
plot.coords(Mand.3d.pc1max, W.Dent, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view1$userMatrix)        
plot.coords(Mand.3d.pc2min, W.Dent, points.col="slategray4", lines.col="slategray4")
plot.coords(Mand.3d.pc2max, W.Dent, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view1$userMatrix) 
plot.coords(Mand.3d.pc1min, W.Dent, points.col="slategray4", lines.col="slategray4")
plot.coords(Mand.3d.pc1max, W.Dent, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view2$userMatrix)        
plot.coords(Mand.3d.pc2min, W.Dent, points.col="slategray4", lines.col="slategray4")
plot.coords(Mand.3d.pc2max, W.Dent, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view2$userMatrix)

# 8.5) Total PC axis wireframes for Fig. S9 (PC1 row 1, PC2 row 2)
layout3d(matrix(1:6, nrow=2), sharedMouse = TRUE) #grid for wireframes of total landmark set
plot.coords(sym.3d.pc1min, W.sym.tot, points.col="slategray4", lines.col="slategray4")
plot.coords(sym.3d.pc1max, W.sym.tot, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view1$userMatrix)        
plot.coords(sym.3d.pc2min, W.sym.tot, points.col="slategray4", lines.col="slategray4")
plot.coords(sym.3d.pc2max, W.sym.tot, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view1$userMatrix) 
plot.coords(sym.3d.pc1min, W.sym.tot, points.col="slategray4", lines.col="slategray4")
plot.coords(sym.3d.pc1max, W.sym.tot, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view2$userMatrix)        
plot.coords(sym.3d.pc2min, W.sym.tot, points.col="slategray4", lines.col="slategray4")
plot.coords(sym.3d.pc2max, W.sym.tot, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view2$userMatrix)   
plot.coords(sym.3d.pc1min, W.sym.tot, points.col="slategray4", lines.col="slategray4")
plot.coords(sym.3d.pc1max, W.sym.tot, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view3$userMatrix)  
plot.coords(sym.3d.pc2min, W.sym.tot, points.col="slategray4", lines.col="slategray4")
plot.coords(sym.3d.pc2max, W.sym.tot, points.col="sienna1", lines.col="sienna1", add = T)
rgl.viewpoint(userMatrix = view3$userMatrix)







#######
#
# Kinematic Vector matrix (Figure S6; here to end)
####

KT.LR.2D.res.0.0<-cbind(kinem.df.res$Op4.2d.KT.res,kinem.df.res$Max4.2d.KT.res,
                        kinem.df.res$Mand.2d.openingLR.res,kinem.df.res$Mand.2d.closingLR.res)
KT.LR.2D.res.15.10<-cbind(kinem.df.res$Op4.2d.KT.15.res,kinem.df.res$Max4.2d.KT.10.res,
                          kinem.df.res$Mand.2d.openingLR.res,kinem.df.res$Mand.2d.closingLR.res)
KT.LR.2D.res.30.20<-cbind(kinem.df.res$Op4.2d.KT.30.res,kinem.df.res$Max4.2d.KT.20.res,
                          kinem.df.res$Mand.2d.openingLR.res,kinem.df.res$Mand.2d.closingLR.res)
KT.LR.3D.res.0.0<-cbind(kinem.df.res$Op4.3d.KT.res,kinem.df.res$Max4.3d.KT.res,
                        kinem.df.res$Mand.3d.openingLR.res,kinem.df.res$Mand.3d.closingLR.res)
KT.LR.3D.res.15.10<-cbind(kinem.df.res$Op4.3d.KT.15.res,kinem.df.res$Max4.3d.KT.10.res,
                          kinem.df.res$Mand.3d.openingLR.res,kinem.df.res$Mand.3d.closingLR.res)
KT.LR.3D.res.30.20<-cbind(kinem.df.res$Op4.3d.KT.30.res,kinem.df.res$Max4.3d.KT.20.res,
                          kinem.df.res$Mand.3d.openingLR.res,kinem.df.res$Mand.3d.closingLR.res)

kinem.df.res.rrpp<-rrpp.data.frame(kinem.df.res)

summary(Mech.res.2D.lmrrpp<-lm.rrpp(KT.LR.2D.res.0.0~Watershed*LD1, 
                                    data = kinem.df.res.rrpp, print.progress = F))
Mech.2d.pairwise<-pairwise(Mech.res.2D.lmrrpp, groups = kinem.df.res.rrpp$Watershed, covariate = kinem.df.res.rrpp$LD1)

summary(Mech.res.3D.lmrrpp<-lm.rrpp(KT.LR.3D.res.0.0~Watershed*LD1, 
                                    data = kinem.df.res.rrpp, print.progress = F))
Mech.3d.pairwise<-pairwise(Mech.res.3D.lmrrpp, groups = kinem.df.res.rrpp$Watershed, covariate = kinem.df.res.rrpp$LD1)

summary(Mech.res.3D.half.lmrrpp<-lm.rrpp(KT.LR.3D.res.15.10~Watershed*LD1, 
                                         data = kinem.df.res.rrpp, print.progress = F))
Mech.3d.half.pairwise<-pairwise(Mech.res.3D.half.lmrrpp, groups = kinem.df.res.rrpp$Watershed, covariate = kinem.df.res.rrpp$LD1)

summary(Mech.res.3D.full.lmrrpp<-lm.rrpp(KT.LR.3D.res.30.20~Watershed*LD1, 
                                         data = kinem.df.res.rrpp, print.progress = F))
Mech.3d.full.pairwise<-pairwise(Mech.res.3D.full.lmrrpp, groups = kinem.df.res.rrpp$Watershed, covariate = kinem.df.res.rrpp$LD1)




Mech.2d.pairwise.angles<-summary(Mech.2d.pairwise, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
Mech.2d.pairwise.p<-summary(Mech.2d.pairwise, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
Mech.2d.pairwise.D<-summary(Mech.2d.pairwise, confidence = 0.95, test.type = "dist")$pairwise.tables$D
Mech.2d.pairwise.Dp<-summary(Mech.2d.pairwise, confidence = 0.95, test.type = "dist")$pairwise.tables$P


get_upper_tri<-function(Mech.2d.pairwise.angles){
  Mech.2d.pairwise.angles[upper.tri(Mech.2d.pairwise.angles)] <- NA
  return(Mech.2d.pairwise.angles)
}
Mech.2d.pairwise.angles.melted<-melt(get_upper_tri(Mech.2d.pairwise.angles))
Mech.2d.pairwise.angles.melted[Mech.2d.pairwise.angles.melted==0]<-NA
Mech.2d.pairwise.angles.matrix<-ggplot(data = Mech.2d.pairwise.angles.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 90, 
                                                    na.value = "transparent", limit = c(0,180), guide_colorbar(title = "Angle (deg)")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,1)), color = "black", size = 3)

Mech.2d.pairwise.p
get_lower_tri<-function(Mech.2d.pairwise.p){
  Mech.2d.pairwise.p[lower.tri(Mech.2d.pairwise.p)] <- NA
  return(Mech.2d.pairwise.p)
}
Mech.2d.pairwise.p.melted<-melt(get_lower_tri(Mech.2d.pairwise.p))
Mech.2d.pairwise.p.melted[Mech.2d.pairwise.p.melted==0]<-NA
Mech.2d.pairwise.p.matrix<-ggplot(data = Mech.2d.pairwise.p.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)


get_upper_tri<-function(Mech.2d.pairwise.D){
  Mech.2d.pairwise.D[upper.tri(Mech.2d.pairwise.D)] <- NA
  return(Mech.2d.pairwise.D)
}
Mech.2d.pairwise.D.melted<-melt(get_upper_tri(Mech.2d.pairwise.D))
Mech.2d.pairwise.D.melted[Mech.2d.pairwise.D.melted==0]<-NA
Mech.2d.pairwise.D.matrix<-ggplot(data = Mech.2d.pairwise.D.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 0.2, 
                                                    na.value = "transparent", limit = c(0,0.4), guide_colorbar(title = "Vector Length diff.")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)

Mech.2d.pairwise.Dp
get_lower_tri<-function(Mech.2d.pairwise.Dp){
  Mech.2d.pairwise.Dp[lower.tri(Mech.2d.pairwise.Dp)] <- NA
  return(Mech.2d.pairwise.Dp)
}
Mech.2d.pairwise.Dp.melted<-melt(get_lower_tri(Mech.2d.pairwise.Dp))
Mech.2d.pairwise.Dp.melted[Mech.2d.pairwise.Dp.melted==0]<-NA
Mech.2d.pairwise.Dp.matrix<-ggplot(data = Mech.2d.pairwise.Dp.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)



Mech.3d.pairwise.angles<-summary(Mech.3d.pairwise, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
Mech.3d.pairwise.p<-summary(Mech.3d.pairwise, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
Mech.3d.pairwise.D<-summary(Mech.3d.pairwise, confidence = 0.95, test.type = "dist")$pairwise.tables$D
Mech.3d.pairwise.Dp<-summary(Mech.3d.pairwise, confidence = 0.95, test.type = "dist")$pairwise.tables$P


get_upper_tri<-function(Mech.3d.pairwise.angles){
  Mech.3d.pairwise.angles[upper.tri(Mech.3d.pairwise.angles)] <- NA
  return(Mech.3d.pairwise.angles)
}
Mech.3d.pairwise.angles.melted<-melt(get_upper_tri(Mech.3d.pairwise.angles))
Mech.3d.pairwise.angles.melted[Mech.3d.pairwise.angles.melted==0]<-NA
Mech.3d.pairwise.angles.matrix<-ggplot(data = Mech.3d.pairwise.angles.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 90, 
                                                    na.value = "transparent", limit = c(0,180), guide_colorbar(title = "Angle (deg)")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,1)), color = "black", size = 3)

Mech.3d.pairwise.p
get_lower_tri<-function(Mech.3d.pairwise.p){
  Mech.3d.pairwise.p[lower.tri(Mech.3d.pairwise.p)] <- NA
  return(Mech.3d.pairwise.p)
}
Mech.3d.pairwise.p.melted<-melt(get_lower_tri(Mech.3d.pairwise.p))
Mech.3d.pairwise.p.melted[Mech.3d.pairwise.p.melted==0]<-NA
Mech.3d.pairwise.p.matrix<-ggplot(data = Mech.3d.pairwise.p.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)


get_upper_tri<-function(Mech.3d.pairwise.D){
  Mech.3d.pairwise.D[upper.tri(Mech.3d.pairwise.D)] <- NA
  return(Mech.3d.pairwise.D)
}
Mech.3d.pairwise.D.melted<-melt(get_upper_tri(Mech.3d.pairwise.D))
Mech.3d.pairwise.D.melted[Mech.3d.pairwise.D.melted==0]<-NA
Mech.3d.pairwise.D.matrix<-ggplot(data = Mech.3d.pairwise.D.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 0.2, 
                                                    na.value = "transparent", limit = c(0,0.4), guide_colorbar(title = "Vector Length diff.")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)

Mech.3d.pairwise.Dp
get_lower_tri<-function(Mech.3d.pairwise.Dp){
  Mech.3d.pairwise.Dp[lower.tri(Mech.3d.pairwise.Dp)] <- NA
  return(Mech.3d.pairwise.Dp)
}
Mech.3d.pairwise.Dp.melted<-melt(get_lower_tri(Mech.3d.pairwise.Dp))
Mech.3d.pairwise.Dp.melted[Mech.3d.pairwise.Dp.melted==0]<-NA
Mech.3d.pairwise.Dp.matrix<-ggplot(data = Mech.3d.pairwise.Dp.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)



Mech.3d.pairwise.half.angles<-summary(Mech.3d.half.pairwise, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
Mech.3d.pairwise.half.p<-summary(Mech.3d.half.pairwise, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
Mech.3d.pairwise.half.D<-summary(Mech.3d.half.pairwise, confidence = 0.95, test.type = "dist")$pairwise.tables$D
Mech.3d.pairwise.half.Dp<-summary(Mech.3d.half.pairwise, confidence = 0.95, test.type = "dist")$pairwise.tables$P


get_upper_tri<-function(Mech.3d.pairwise.half.angles){
  Mech.3d.pairwise.half.angles[upper.tri(Mech.3d.pairwise.half.angles)] <- NA
  return(Mech.3d.pairwise.half.angles)
}
Mech.3d.pairwise.half.angles.melted<-melt(get_upper_tri(Mech.3d.pairwise.half.angles))
Mech.3d.pairwise.half.angles.melted[Mech.3d.pairwise.half.angles.melted==0]<-NA
Mech.3d.pairwise.half.angles.matrix<-ggplot(data = Mech.3d.pairwise.half.angles.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 90, 
                                                    na.value = "transparent", limit = c(0,180), guide_colorbar(title = "Angle (deg)")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,1)), color = "black", size = 3)


get_lower_tri<-function(Mech.3d.pairwise.half.p){
  Mech.3d.pairwise.half.p[lower.tri(Mech.3d.pairwise.half.p)] <- NA
  return(Mech.3d.pairwise.half.p)
}
Mech.3d.pairwise.half.p.melted<-melt(get_lower_tri(Mech.3d.pairwise.half.p))
Mech.3d.pairwise.half.p.melted[Mech.3d.pairwise.half.p.melted==0]<-NA
Mech.3d.pairwise.half.p.matrix<-ggplot(data = Mech.3d.pairwise.half.p.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)


get_upper_tri<-function(Mech.3d.pairwise.half.D){
  Mech.3d.pairwise.half.D[upper.tri(Mech.3d.pairwise.half.D)] <- NA
  return(Mech.3d.pairwise.half.D)
}
Mech.3d.pairwise.half.D.melted<-melt(get_upper_tri(Mech.3d.pairwise.half.D))
Mech.3d.pairwise.half.D.melted[Mech.3d.pairwise.half.D.melted==0]<-NA
Mech.3d.pairwise.half.D.matrix<-ggplot(data = Mech.3d.pairwise.half.D.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 0.2, 
                                                    na.value = "transparent", limit = c(0,0.4), guide_colorbar(title = "Vector Length diff.")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)

Mech.3d.pairwise.half.Dp
get_lower_tri<-function(Mech.3d.pairwise.half.Dp){
  Mech.3d.pairwise.half.Dp[lower.tri(Mech.3d.pairwise.half.Dp)] <- NA
  return(Mech.3d.pairwise.half.Dp)
}
Mech.3d.pairwise.half.Dp.melted<-melt(get_lower_tri(Mech.3d.pairwise.half.Dp))
Mech.3d.pairwise.half.Dp.melted[Mech.3d.pairwise.half.Dp.melted==0]<-NA
Mech.3d.pairwise.half.Dp.matrix<-ggplot(data = Mech.3d.pairwise.half.Dp.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)



Mech.3d.pairwise.full.angles<-summary(Mech.3d.full.pairwise, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
Mech.3d.pairwise.full.p<-summary(Mech.3d.full.pairwise, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
Mech.3d.pairwise.full.D<-summary(Mech.3d.full.pairwise, confidence = 0.95, test.type = "dist")$pairwise.tables$D
Mech.3d.pairwise.full.Dp<-summary(Mech.3d.full.pairwise, confidence = 0.95, test.type = "dist")$pairwise.tables$P


get_upper_tri<-function(Mech.3d.pairwise.full.angles){
  Mech.3d.pairwise.full.angles[upper.tri(Mech.3d.pairwise.full.angles)] <- NA
  return(Mech.3d.pairwise.full.angles)
}
Mech.3d.pairwise.full.angles.melted<-melt(get_upper_tri(Mech.3d.pairwise.full.angles))
Mech.3d.pairwise.full.angles.melted[Mech.3d.pairwise.full.angles.melted==0]<-NA
Mech.3d.pairwise.full.angles.matrix<-ggplot(data = Mech.3d.pairwise.full.angles.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 90, 
                                                    na.value = "transparent", limit = c(0,180), guide_colorbar(title = "Angle (deg)")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,1)), color = "black", size = 3)


get_lower_tri<-function(Mech.3d.pairwise.full.p){
  Mech.3d.pairwise.full.p[lower.tri(Mech.3d.pairwise.full.p)] <- NA
  return(Mech.3d.pairwise.full.p)
}
Mech.3d.pairwise.full.p.melted<-melt(get_lower_tri(Mech.3d.pairwise.full.p))
Mech.3d.pairwise.full.p.melted[Mech.3d.pairwise.full.p.melted==0]<-NA
Mech.3d.pairwise.full.p.matrix<-ggplot(data = Mech.3d.pairwise.full.p.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)


get_upper_tri<-function(Mech.3d.pairwise.full.D){
  Mech.3d.pairwise.full.D[upper.tri(Mech.3d.pairwise.full.D)] <- NA
  return(Mech.3d.pairwise.full.D)
}
Mech.3d.pairwise.full.D.melted<-melt(get_upper_tri(Mech.3d.pairwise.full.D))
Mech.3d.pairwise.full.D.melted[Mech.3d.pairwise.full.D.melted==0]<-NA
Mech.3d.pairwise.full.D.matrix<-ggplot(data = Mech.3d.pairwise.full.D.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 0.2, 
                                                    na.value = "transparent", limit = c(0,0.4), guide_colorbar(title = "Vector Length diff.")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)

Mech.3d.pairwise.full.Dp
get_lower_tri<-function(Mech.3d.pairwise.full.Dp){
  Mech.3d.pairwise.full.Dp[lower.tri(Mech.3d.pairwise.full.Dp)] <- NA
  return(Mech.3d.pairwise.full.Dp)
}
Mech.3d.pairwise.full.Dp.melted<-melt(get_lower_tri(Mech.3d.pairwise.full.Dp))
Mech.3d.pairwise.full.Dp.melted[Mech.3d.pairwise.full.Dp.melted==0]<-NA
Mech.3d.pairwise.full.Dp.matrix<-ggplot(data = Mech.3d.pairwise.full.Dp.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)



kin.angle.matrix<-ggarrange(Mech.2d.pairwise.angles.matrix, Mech.3d.pairwise.angles.matrix, Mech.3d.pairwise.half.angles.matrix,Mech.3d.pairwise.full.angles.matrix,
                            ncol = 4, nrow = 1, common.legend = T, legend = "left")
kin.anglepval.matrix<-ggarrange(Mech.2d.pairwise.p.matrix, Mech.3d.pairwise.p.matrix, Mech.3d.pairwise.half.p.matrix,Mech.3d.pairwise.full.p.matrix,
                                ncol = 4, nrow = 1, common.legend = T, legend = "right")

kin.dist.matrix<-ggarrange(Mech.2d.pairwise.D.matrix, Mech.3d.pairwise.D.matrix, Mech.3d.pairwise.half.D.matrix,Mech.3d.pairwise.full.D.matrix,
                           ncol = 4, nrow = 1, common.legend = T, legend = "left")
kin.distpval.matrix<-ggarrange(Mech.2d.pairwise.Dp.matrix, Mech.3d.pairwise.Dp.matrix, Mech.3d.pairwise.half.Dp.matrix,Mech.3d.pairwise.full.Dp.matrix,
                               ncol = 4, nrow = 1, common.legend = T, legend = "right")
ggarrange(kin.anglepval.matrix,kin.distpval.matrix, ncol = 1, nrow = 2, widths = c(1,1,1,1))


# Diet LD covariate angle
mech.angle.matrix<-ggarrange(Mech.2d.pairwise.angles.matrix, Mech.3d.pairwise.angles.matrix, 
                             Mech.3d.pairwise.half.angles.matrix, Mech.3d.pairwise.full.angles.matrix,
                             ncol = 4, nrow = 1, common.legend = T, legend = "left")
# Diet LD covariate vector magnitude
mech.D.matrix<-ggarrange( Mech.2d.pairwise.D.matrix, Mech.3d.pairwise.D.matrix, 
                          Mech.3d.pairwise.half.D.matrix, Mech.3d.pairwise.full.D.matrix,
                          ncol = 4, nrow = 1, common.legend = T, legend = "left")
# Diet LD covariate angle p-val
mech.p.matrix<-ggarrange(Mech.2d.pairwise.p.matrix, Mech.3d.pairwise.p.matrix, 
                         Mech.3d.pairwise.half.p.matrix, Mech.3d.pairwise.full.p.matrix,
                         ncol = 4, nrow = 1, common.legend = T, legend = "right")
# Diet LD covariate vector magnitude p-val
mech.Dp.matrix<-ggarrange(Mech.2d.pairwise.Dp.matrix, Mech.3d.pairwise.Dp.matrix, 
                          Mech.3d.pairwise.half.Dp.matrix, Mech.3d.pairwise.full.Dp.matrix,
                          ncol = 4, nrow = 1, common.legend = T, legend = "right")

