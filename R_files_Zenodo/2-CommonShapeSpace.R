library(geomorph)
library(Morpho)
library(vegan)
library(ggplot2)
library(MASS)
library(shapes)
library(abind)
library(RRPP)
library(cowplot)


# ANALYSES IN COMMON 2D-3D SHAPE SPACE

# Inherited from "1-SeparateShapeSpaces.R"

# Op4.2dC.res, Op4bar.2d.gpa
# Op4.3dC.res, Op4bar.3d.gpa
# Dent.2dC.res, Dent.2d.gpa
# Dent.3dC.res, Dent.3d.gpa
# Max4.2dC.res, Max4bar.2d.gpa
# Max4.3dC.res, Max4bar.3d.gpa
# spec.sym.2d.res, spec.sym.2d.gpa
# spec.sym.3d.res, spec.sym.3d.gpa


#  1. (FOLLOWING CARDINI '14 AND  CARDINI & CHIAPPELLI '19) TO PRODUCE COMMON SHAPE SPACE  for Op4 
#   (uses both sides of symmetric landmarks)
Op4bar.2d.sym
Op4bar.2d3d.sym<-abind(abind(Op4bar.2d.sym, array(0, replace(dim(Op4bar.2d.sym), 2, 1)), along = 2),Op4bar.3d.sym, along=3)

Op4bar2d3d.gpa<-gpagen(Op4bar.2d3d.sym)
Op4.geodf.2d<-geomorph.data.frame(coords=Op4bar2d3d.gpa$coords[,,1:59],
                               Csize = Op4bar2d3d.gpa$Csize[1:59],
                               Dim = rep("2", times = 59),
                               spec = ldmk.factors$Spec)
Op4.geodf.3d<-geomorph.data.frame(coords=Op4bar2d3d.gpa$coords[,,60:118],
                                  Csize = Op4bar2d3d.gpa$Csize[60:118],
                                  Dim = rep("3", times = 59),
                                  spec = ldmk.factors$Spec)

Op4.geodf.2d.res<-procD.lm(coords~Csize,data =Op4.geodf.2d)$residuals
Op4.geodf.3d.res<-procD.lm(coords~Csize,data =Op4.geodf.3d)$residuals

Op4.geodf<-geomorph.data.frame(coords=abind(arrayspecs(Op4.geodf.2d.res,8,3),arrayspecs(Op4.geodf.3d.res,8,3), along=3),
                               pop = rep(ldmk.factors$Population, times = 2),
                               watershed = rep(ldmk.factors$Watershed, times = 2), 
                               Csize = Op4bar2d3d.gpa$Csize,
                               habitat = rep(ldmk.factors$Habitat, times = 2), 
                               shannon = rep(ldmk.factors$shannon, times = 2),
                               simpson = rep(ldmk.factors$simpson, times = 2), 
                               dietPC1 = rep(ldmk.factors$PC1, times = 2),
                               dietLD1 = rep(ldmk.factors$LD1, times = 2),
                               Dim = append(rep("2", times = 59),rep("3", times = 59)),
                               spec = rownames(as.matrix(Op4bar2d3d.gpa$Csize)))



# 1) Puts Op4 Csize residuals for each lever system into common shape space (does the same as above, 
#    but by binding 2D and 3D residuals created in code file 1). This method will be used subsequently, as it uses simpler code.
#    (uses only one side of symmetric landmarks)
Op4.geodf<-geomorph.data.frame(coords=abind(abind(Op4.2dC.res, array(0, replace(dim(Op4.2dC.res), 2, 1)), along = 2),Op4.3dC.res, along=3),
                               pop = rep(ldmk.factors$Population, times = 2),
                               watershed = rep(ldmk.factors$Watershed, times = 2), 
                               Csize = append(Op4bar.2d.gpa$Csize, Op4bar.3d.gpa$Csize, after = 59),
                               habitat = rep(ldmk.factors$Habitat, times = 2), 
                               shannon = rep(ldmk.factors$shannon, times = 2),
                               simpson = rep(ldmk.factors$simpson, times = 2), 
                               dietPC1 = rep(ldmk.factors$PC1, times = 2),
                               dietLD1 = rep(ldmk.factors$LD1, times = 2),
                               Dim = append(rep("2", times = 59),rep("3", times = 59)),
                               spec = rownames(as.matrix(append(Op4bar.2d.gpa$Csize, Op4bar.3d.gpa$Csize, after = 59))))

# 1.1) Fig 4A
Op4.2d3d.tanpca<-plotTangentSpace(Op4.geodf$coords)
Op4.2d3d.pcaplot<- ggplot(data = data.frame(Op4.2d3d.tanpca$pc.scores,Op4.geodf$Dim,Op4.geodf$spec, 
                                            Op4.geodf$habitat, Op4.geodf$watershed, Op4.geodf$dietLD1), 
                          aes(x=PC1, y=PC2))+geom_point()+
  xlab("PC1 - 57.7%")+ylab("PC2 - 18.5%") +
  geom_line(aes(group = Op4.geodf$spec)) + theme_cowplot()
Op4.2d3d.pcaplot<- Op4.2d3d.pcaplot+ aes(colour=Op4.geodf$dietLD1, shape=Op4.geodf$Dim, stroke = 1.5) +
  scale_colour_gradient2(low = "red",mid = "gray98", high = "steelblue", midpoint = .75, name = "Diet LD") +
  scale_shape_discrete(name = "Landmark Dim.")
Op4.2d3d.pcaplot
# 1.2) Fig 4E
Op4.2d3d.simp.pcaplot<- ggplot(data = data.frame(Op4.2d3d.tanpca$pc.scores,Op4.geodf$Dim,Op4.geodf$spec, 
                                                 Op4.geodf$habitat, Op4.geodf$watershed, Op4.geodf$simpson), 
                               aes(x=PC1, y=PC2))+geom_point()+
  xlab("PC1 - 57.7%")+ylab("PC2 - 18.5%") +
  geom_line(aes(group = Op4.geodf$spec))+ theme_cowplot()
Op4.2d3d.simp.pcaplot<- Op4.2d3d.simp.pcaplot+ aes(colour=Op4.geodf$simpson, shape=Op4.geodf$Dim, stroke = 1.5) +
  scale_colour_gradient2(low = "red",mid = "gray98", high = "steelblue", midpoint = .86, name = "Simpson Div.") +
  scale_shape_discrete(name = "Landmark Dim.")
Op4.2d3d.simp.pcaplot

# 1.3) Loop calculating Procrustes distances between 2D and 3D versions of same specimen
Op4.2d3d.procdist<-1:59
x <-1
repeat{
  Op4.2d3d.procdist[x]<-procdist(Op4.geodf$coords[,,x],Op4.geodf$coords[,,(x+59)], type = "full")
  x<- x+1
  if(x>59){break}
}
Op4.procdist.df<-data.frame(procdist = Op4.2d3d.procdist,ldmk.factors)
# 1.4) ANCOVA calculating effects of diet and watershed on Procrustes distance between 2D and 3D versions
Anova(lm(procdist ~ Watershed * LD1, data = Op4.procdist.df), type = "II")
Anova(lm(procdist ~ Watershed * simpson, data = Op4.procdist.df), type = "II")

# 1.5) Op4 vector comparisons in common shape space
summary(Op4.WatxLDxDim<-procD.lm(coords~watershed*dietLD1*Dim, SS.type = "II", data = Op4.geodf, logsz = T, print.progress = F))
Op4.WatxLDxDim.PW<-pairwise(Op4.WatxLDxDim, groups = interaction(Op4.geodf$watershed,Op4.geodf$Dim), covariate = Op4.geodf$dietLD1)
summary(Op4.WatxLDxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(Op4.WatxLDxDim.PW, confidence = 0.95, test.type = "dist")

summary(Op4.LDxDim<-procD.lm(coords~dietLD1*Dim, SS.type = "II", data = Op4.geodf, logsz = T, print.progress = F))
Op4.LDxDim.PW<-pairwise(Op4.LDxDim, groups = Op4.geodf$Dim, covariate = Op4.geodf$dietLD1)
summary(Op4.LDxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(Op4.LDxDim.PW, confidence = 0.95, test.type = "dist")

summary(Op4.WatxSimpxDim<-procD.lm(coords~watershed*simpson*Dim, SS.type = "II", data = Op4.geodf, logsz = T, print.progress = F))
Op4.WatxSimpxDim.PW<-pairwise(Op4.WatxSimpxDim, groups = interaction(Op4.geodf$watershed,Op4.geodf$Dim), covariate = Op4.geodf$simpson)
summary(Op4.WatxSimpxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(Op4.WatxSimpxDim.PW, confidence = 0.95, test.type = "dist")

summary(Op4.SimpxDim<-procD.lm(coords~simpson*Dim, SS.type = "II", data = Op4.geodf, logsz = T, print.progress = F))
Op4.SimpxDim.PW<-pairwise(Op4.SimpxDim, groups = Op4.geodf$Dim, covariate = Op4.geodf$simpson)
summary(Op4.SimpxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(Op4.SimpxDim.PW, confidence = 0.95, test.type = "dist")


#
#
# 2) Puts Mandible Csize residuals for each lever system into common shape space
Dent.geodf<-geomorph.data.frame(coords=abind(abind(Dent.2dC.res, array(0, replace(dim(Dent.2dC.res), 2, 1)), along = 2),Dent.3dC.res, along=3),
                                pop = rep(ldmk.factors$Population, times = 2),
                                watershed = rep(ldmk.factors$Watershed, times = 2), 
                                Csize = append(Dent.2d.gpa$Csize, Dent.3d.gpa$Csize, after = 59),
                                habitat = rep(ldmk.factors$Habitat, times = 2), 
                                shannon = rep(ldmk.factors$shannon, times = 2),
                                simpson = rep(ldmk.factors$simpson, times = 2), 
                                dietPC1 = rep(ldmk.factors$PC1, times = 2),
                                dietLD1 = rep(ldmk.factors$LD1, times = 2),
                                Dim = append(rep("2", times = 59),rep("3", times = 59)),
                                spec = rownames(as.matrix(append(Dent.2d.gpa$Csize, Dent.3d.gpa$Csize, after = 59))))

# 2.1) Fig 4B
Dent.2d3d.tanpca<-plotTangentSpace(Dent.geodf$coords)
Dent.2d3d.pcaplot<- ggplot(data = data.frame(Dent.2d3d.tanpca$pc.scores,Dent.geodf$Dim,Dent.geodf$spec, 
                                             Dent.geodf$habitat, Dent.geodf$watershed, Dent.geodf$dietLD1), 
                           aes(x=PC1, y=PC2))+geom_point()+
  xlab("PC1 - 60.4%")+ylab("PC2 - 14.5%") +
  geom_line(aes(group = Dent.geodf$spec))+ theme_cowplot()
Dent.2d3d.pcaplot<- Dent.2d3d.pcaplot+ aes(colour=Dent.geodf$dietLD1, shape=Dent.geodf$Dim, stroke = 1.5) +
  scale_colour_gradient2(low = "red",mid = "gray98", high = "steelblue", midpoint = .75, name = "Diet LD") +
  scale_shape_discrete(name = "Landmark Dim.")
Dent.2d3d.pcaplot
# 2.2) Fig 4F
Dent.2d3d.simp.pcaplot<- ggplot(data = data.frame(Dent.2d3d.tanpca$pc.scores,Dent.geodf$Dim,Dent.geodf$spec, 
                                                  Dent.geodf$habitat, Dent.geodf$watershed, Dent.geodf$simpson), 
                                aes(x=PC1, y=PC2))+geom_point()+
  xlab("PC1 - 60.4%")+ylab("PC2 - 14.5%") +
  geom_line(aes(group = Dent.geodf$spec))+ theme_cowplot()
Dent.2d3d.simp.pcaplot<- Dent.2d3d.simp.pcaplot+ aes(colour=Dent.geodf$simpson, shape=Dent.geodf$Dim, stroke = 1.5) +
  scale_colour_gradient2(low = "red",mid = "gray98", high = "steelblue", midpoint = .86, name = "Simpson Div.") +
  scale_shape_discrete(name = "Landmark Dim.")
Dent.2d3d.simp.pcaplot
plot(ldmk.factors$simpson, Dent.2d3d.procdist)
# 2.3) Loop calculating Procrustes distances between 2D and 3D versions of same specimen
Dent.2d3d.procdist<-1:59
x <-1
repeat{
  Dent.2d3d.procdist[x]<-procdist(Dent.geodf$coords[,,x],Dent.geodf$coords[,,(x+59)], type = "full")
  x<- x+1
  if(x>59){break}
}
plot(ldmk.factors$LD1,Dent.2d3d.procdist)
Dent.procdist.df<-data.frame(procdist = Dent.2d3d.procdist,ldmk.factors)
# 2.4) ANCOVAs calculating effects of diet and watershed on Procrustes distance between 2D and 3D versions
Anova(lm(procdist ~ Watershed * LD1, data = Dent.procdist.df), type = "II")
Anova(lm(procdist ~ Watershed * simpson, data = Dent.procdist.df), type = "II")
# 2.5) Mandibular lever vector comparisons in common shape space
summary(Dent.WatxLDxDim<-procD.lm(coords~watershed*dietLD1*Dim, SS.type = "II", data = Dent.geodf, logsz = T, print.progress = F))
Dent.WatxLDxDim.PW<-pairwise(Dent.WatxLDxDim, groups = interaction(Dent.geodf$watershed,Dent.geodf$Dim), covariate = Dent.geodf$dietLD1)
summary(Dent.WatxLDxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(Dent.WatxLDxDim.PW, confidence = 0.95, test.type = "dist")

summary(Dent.LDxDim<-procD.lm(coords~dietLD1*Dim, SS.type = "II", data = Dent.geodf, logsz = T, print.progress = F))
Dent.LDxDim.PW<-pairwise(Dent.LDxDim, groups = Dent.geodf$Dim, covariate = Dent.geodf$dietLD1)
summary(Dent.LDxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(Dent.LDxDim.PW, confidence = 0.95, test.type = "dist")

summary(Dent.WatxSimpxDim<-procD.lm(coords~watershed*simpson*Dim, SS.type = "II", data = Dent.geodf, logsz = T, print.progress = F))
Dent.WatxSimpxDim.PW<-pairwise(Dent.WatxSimpxDim, groups = interaction(Dent.geodf$watershed,Dent.geodf$Dim), covariate = Dent.geodf$simpson)
summary(Dent.WatxSimpxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(Dent.WatxSimpxDim.PW, confidence = 0.95, test.type = "dist")

summary(Dent.SimpxDim<-procD.lm(coords~simpson*Dim, SS.type = "II", data = Dent.geodf, logsz = T, print.progress = F))
Dent.SimpxDim.PW<-pairwise(Dent.SimpxDim, groups = Dent.geodf$Dim, covariate = Dent.geodf$simpson)
summary(Dent.SimpxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(Dent.SimpxDim.PW, confidence = 0.95, test.type = "dist")

#
#
# 3) Puts Max4 Csize residuals for each lever system into common shape space
Max4.geodf<-geomorph.data.frame(coords=abind(abind(Max4.2dC.res, array(0, replace(dim(Max4.2dC.res), 2, 1)), along = 2),Max4.3dC.res, along=3),
                                pop = rep(ldmk.factors$Population, times = 2),
                                watershed = rep(ldmk.factors$Watershed, times = 2), 
                                Csize = append(Max4bar.2d.gpa$Csize, Max4bar.3d.gpa$Csize, after = 59),
                                habitat = rep(ldmk.factors$Habitat, times = 2), 
                                shannon = rep(ldmk.factors$shannon, times = 2),
                                simpson = rep(ldmk.factors$simpson, times = 2), 
                                dietPC1 = rep(ldmk.factors$PC1, times = 2),
                                dietLD1 = rep(ldmk.factors$LD1, times = 2), 
                                Dim = append(rep("2", times = 59),rep("3", times = 59)),
                                spec = rownames(as.matrix(append(Max4bar.2d.gpa$Csize, Max4bar.3d.gpa$Csize, after = 59))))
# 3.1) Fig 4C
Max4.2d3d.tanpca<-plotTangentSpace(Max4.geodf$coords)
Max4.2d3d.pcaplot<- ggplot(data = data.frame(Max4.2d3d.tanpca$pc.scores,Max4.geodf$Dim,Max4.geodf$spec, 
                                             Max4.geodf$habitat, Max4.geodf$watershed, Max4.geodf$dietLD1), 
                           aes(x=PC1, y=PC2))+geom_point()+
  xlab("PC1 - 37.3%")+ylab("PC2 - 31.0%")+
  geom_line(aes(group = Max4.geodf$spec))+ theme_cowplot()
Max4.2d3d.pcaplot<- Max4.2d3d.pcaplot+ aes(colour=Max4.geodf$dietLD1, shape=Max4.geodf$Dim, stroke = 1.5) +
  scale_colour_gradient2(low = "red",mid = "gray98", high = "steelblue", midpoint = .75,name = "Diet LD") +
  scale_shape_discrete(name = "Landmark Dim.")
Max4.2d3d.pcaplot
# 3.2) Fig 4G
Max4.2d3d.simp.pcaplot<- ggplot(data = data.frame(Max4.2d3d.tanpca$pc.scores,Max4.geodf$Dim,Max4.geodf$spec, 
                                                  Max4.geodf$habitat, Max4.geodf$watershed, Max4.geodf$simpson), 
                                aes(x=PC1, y=PC2))+geom_point()+
  xlab("PC1 - 37.3%")+ylab("PC2 - 31.0%")+
  geom_line(aes(group = Max4.geodf$spec))+ theme_cowplot()
Max4.2d3d.simp.pcaplot<- Max4.2d3d.simp.pcaplot+ aes(colour=Max4.geodf$simpson, shape=Max4.geodf$Dim, stroke = 1.5) +
  scale_colour_gradient2(low = "red",mid = "gray98", high = "steelblue", midpoint = .86,name = "Simpson Div.") +
  scale_shape_discrete(name = "Landmark Dim.")
Max4.2d3d.simp.pcaplot
plot(ldmk.factors$simpson, Max4.2d3d.procdist)
# 3.3) Loop calculating Procrustes distances between 2D and 3D versions of same specimen
Max4.2d3d.procdist<-1:59
x <-1
repeat{
  Max4.2d3d.procdist[x]<-procdist(Max4.geodf$coords[,,x],Max4.geodf$coords[,,(x+59)], type = "full")
  x<- x+1
  if(x>59){break}
}
plot(ldmk.factors$LD1,Max4.2d3d.procdist)
Max4.procdist.df<-data.frame(procdist = Max4.2d3d.procdist,ldmk.factors)
# 3.4) ANCOVAs calculating effects of diet and watershed on Procrustes distance between 2D and 3D versions
Anova(lm(procdist ~ Watershed * LD1, data = Max4.procdist.df), type = "II")
Anova(lm(procdist ~ Watershed * simpson, data = Max4.procdist.df), type = "II")
# 3.5) Max4 lever vector comparisons in common shape space
summary(Max4.WatxLDxDim<-procD.lm(coords~watershed*dietLD1*Dim, SS.type = "II", data = Max4.geodf, logsz = T, print.progress = F))
Max4.WatxLDxDim.PW<-pairwise(Max4.WatxLDxDim, groups = interaction(Max4.geodf$watershed,Max4.geodf$Dim), covariate = Max4.geodf$dietLD1)
summary(Max4.WatxLDxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(Max4.WatxLDxDim.PW, confidence = 0.95, test.type = "dist")

summary(Max4.LDxDim<-procD.lm(coords~dietLD1*Dim, SS.type = "II", data = Max4.geodf, logsz = T, print.progress = F))
Max4.LDxDim.PW<-pairwise(Max4.LDxDim, groups = Max4.geodf$Dim, covariate = Max4.geodf$dietLD1)
summary(Max4.LDxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg", show.vectors = T)
summary(Max4.LDxDim.PW, confidence = 0.95, test.type = "dist")

summary(Max4.WatxSimpxDim<-procD.lm(coords~watershed*simpson*Dim, SS.type = "II", data = Max4.geodf, logsz = T, print.progress = F))
Max4.WatxSimpxDim.PW<-pairwise(Max4.WatxSimpxDim, groups = interaction(Max4.geodf$watershed,Max4.geodf$Dim), covariate = Max4.geodf$simpson)
summary(Max4.WatxSimpxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(Max4.WatxSimpxDim.PW, confidence = 0.95, test.type = "dist")

summary(Max4.SimpxDim<-procD.lm(coords~simpson*Dim, SS.type = "II", data = Max4.geodf, logsz = T, print.progress = F))
Max4.SimpxDim.PW<-pairwise(Max4.SimpxDim, groups = Max4.geodf$Dim, covariate = Max4.geodf$simpson)
summary(Max4.SimpxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(Max4.SimpxDim.PW, confidence = 0.95, test.type = "dist")

#
#
# 4) Puts Total landmark Csize residuals for each lever system into common shape space
tot.geodf<-geomorph.data.frame(coords=abind(abind(spec.sym.2d.res, array(0, replace(dim(spec.sym.2d.res), 2, 1)), along = 2),spec.sym.3d.res, along=3),
                               pop = rep(ldmk.factors$Population, times = 2),
                               watershed = rep(ldmk.factors$Watershed, times = 2), 
                               Csize = append(spec.sym.2d.gpa$Csize, spec.sym.3d.gpa$Csize, after = 59),
                               habitat = rep(ldmk.factors$Habitat, times = 2), 
                               shannon = rep(ldmk.factors$shannon, times = 2),
                               simpson = rep(ldmk.factors$simpson, times = 2), 
                               dietPC1 = rep(ldmk.factors$PC1, times = 2),
                               dietLD1 = rep(ldmk.factors$LD1, times = 2), 
                               Dim = append(rep("2", times = 59),rep("3", times = 59)),
                               spec = rownames(as.matrix(append(spec.sym.2d.gpa$Csize, spec.sym.3d.gpa$Csize, after = 59))))

# 4.1) Fig 4D
tot.2d3d.tanpca<-plotTangentSpace(tot.geodf$coords)
tot.2d3d.pcaplot<- ggplot(data = data.frame(tot.2d3d.tanpca$pc.scores,tot.geodf$Dim,tot.geodf$spec, 
                                            tot.geodf$habitat, tot.geodf$watershed, tot.geodf$dietLD1), 
                          aes(x=PC1, y=PC2))+geom_point()+
  xlab("PC1 - 32.2%")+ylab("PC2 - 19.4%")+
  geom_line(aes(group = Op4.geodf$spec)) + theme_cowplot()
tot.2d3d.pcaplot<- tot.2d3d.pcaplot+ aes(colour=tot.geodf$dietLD1, shape=tot.geodf$Dim, stroke = 1.5) +
  scale_colour_gradient2(low = "red",mid = "gray98", high = "steelblue", midpoint = .75,name = "Diet LD") +
  scale_shape_discrete(name = "Landmark Dim.")
tot.2d3d.pcaplot
# 4.2) Fig 4H
tot.2d3d.simp.pcaplot<- ggplot(data = data.frame(tot.2d3d.tanpca$pc.scores,tot.geodf$Dim,tot.geodf$spec, 
                                                 tot.geodf$habitat, tot.geodf$watershed, tot.geodf$simpson), 
                               aes(x=PC1, y=PC2))+geom_point()+
  xlab("PC1 - 32.2%")+ylab("PC2 - 19.4%")+
  geom_line(aes(group = Op4.geodf$spec))
tot.2d3d.simp.pcaplot<- tot.2d3d.simp.pcaplot+ aes(colour=tot.geodf$simpson, shape=tot.geodf$Dim, stroke = 1.5) +
  scale_colour_gradient2(low = "red",mid = "gray98", high = "steelblue", midpoint = .86,name = "Simpson Div.") +
  scale_shape_discrete(name = "Landmark Dim.")+ theme_cowplot()
tot.2d3d.simp.pcaplot
plot(ldmk.factors$simpson, tot.2d3d.procdist)
# 4.3) Loop calculating Procrustes distances between 2D and 3D versions of same specimen
tot.2d3d.procdist<-1:59
x <-1
repeat{
  tot.2d3d.procdist[x]<-procdist(tot.geodf$coords[,,x],tot.geodf$coords[,,(x+59)], type = "full")
  x<- x+1
  if(x>59){break}
}
tot.procdist.df<-data.frame(procdist = tot.2d3d.procdist,ldmk.factors)
# 4.4) ANCOVAs calculating effects of diet and watershed on Procrustes distance between 2D and 3D versions
Anova(lm(procdist ~ Watershed * LD1, data = tot.procdist.df), type = "II")
summary(lm(procdist ~ Watershed * simpson, data = tot.procdist.df), type = "II")

# 4.5) TOTAL lever vector comparisons in common shape space
summary(tot.WatxLDxDim<-procD.lm(coords~watershed*dietLD1*Dim, SS.type = "II", data = tot.geodf, logsz = T, print.progress = F))
tot.WatxLDxDim.PW<-pairwise(tot.WatxLDxDim, groups = interaction(tot.geodf$watershed,tot.geodf$Dim), covariate = tot.geodf$dietLD1)
summary(tot.WatxLDxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(tot.WatxLDxDim.PW, confidence = 0.95, test.type = "dist")

summary(tot.LDxDim<-procD.lm(coords~dietLD1*Dim, SS.type = "II", data = tot.geodf, logsz = T, print.progress = F))
tot.LDxDim.PW<-pairwise(tot.LDxDim, groups = tot.geodf$Dim, covariate = tot.geodf$dietLD1)
summary(tot.LDxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg", show.vectors = T)
summary(tot.LDxDim.PW, confidence = 0.95, test.type = "dist")

summary(tot.WatxSimpxDim<-procD.lm(coords~watershed*simpson*Dim, SS.type = "II", data = tot.geodf, logsz = T, print.progress = F))
tot.WatxSimpxDim.PW<-pairwise(tot.WatxSimpxDim, groups = interaction(tot.geodf$watershed,tot.geodf$Dim), covariate = tot.geodf$simpson)
summary(tot.WatxSimpxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(tot.WatxSimpxDim.PW, confidence = 0.95, test.type = "dist")

summary(tot.SimpxDim<-procD.lm(coords~simpson*Dim, SS.type = "II", data = tot.geodf, logsz = T, print.progress = F))
tot.SimpxDim.PW<-pairwise(tot.SimpxDim, groups = tot.geodf$Dim, covariate = tot.geodf$simpson)
summary(tot.SimpxDim.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")
summary(tot.SimpxDim.PW, confidence = 0.95, test.type = "dist")

#
#
#
#
# 5) Fig 4: COMBINED 2D3D SHAPE SPACE PLOTS
ggarrange(ggarrange(Op4.2d3d.pcaplot,Dent.2d3d.pcaplot,Max4.2d3d.pcaplot,tot.2d3d.pcaplot,ncol = 4, nrow = 1,
                    labels = c("A","B", "C","D"), hjust = -2,
                    align = "hv",legend = "right", common.legend = T),
          ggarrange(Op4.2d3d.simp.pcaplot,Dent.2d3d.simp.pcaplot,Max4.2d3d.simp.pcaplot,tot.2d3d.simp.pcaplot,
                    ncol = 4, nrow = 1,
                    labels = c("E","F","G","H"), hjust = -2,
                    align = "hv",legend = "right", common.legend = T),ncol = 1, nrow = 2)


#6) JEvoBio Graphical Abstract
graphAbsA<- ggplot(data = data.frame(tot.2d3d.tanpca$pc.scores,tot.geodf$Dim,tot.geodf$spec, 
                                            tot.geodf$habitat, tot.geodf$watershed, tot.geodf$dietLD1), 
                          aes(x=PC1, y=PC2))+geom_point()+
  xlab("PC1")+ylab("PC2")+
  geom_line(aes(group = Op4.geodf$spec)) + theme_cowplot()
graphAbsA<- graphAbsA+ aes(colour=tot.geodf$dietLD1, shape=tot.geodf$Dim, stroke = 1.5) +
  scale_colour_gradient2(low = "red",mid = "gray98", high = "steelblue", midpoint = .75,name = "Diet Type") +
  scale_shape_discrete(name = "Landmark Dim.")
graphAbsB<- ggplot(data = data.frame(tot.2d3d.tanpca$pc.scores,tot.geodf$Dim,tot.geodf$spec, 
                                                 tot.geodf$habitat, tot.geodf$watershed, tot.geodf$simpson), 
                               aes(x=PC1, y=PC2))+geom_point()+
  xlab("PC1")+ylab("PC2")+
  geom_line(aes(group = Op4.geodf$spec))
graphAbsB<- graphAbsB+ aes(colour=tot.geodf$simpson, shape=tot.geodf$Dim, stroke = 1.5) +
  scale_colour_gradient2(low = "red",mid = "gray98", high = "steelblue", midpoint = .86,name = "Diet Diversity") +
  scale_shape_discrete(name = "Landmark Dim.")+ theme_cowplot()
ggarrange(graphAbsA,graphAbsB,ncol = 1, nrow = 2,
          labels = c("A","B"), hjust = -2,
          align = "hv",legend = "right", common.legend = F)



# 1.4,2.4,3.4,4.4 #THESE TESTS ALSO PERFORMED ABOVE, just put in this block for convenience
Anova(lm(procdist ~ Watershed * LD1, data = Op4.procdist.df), type = "II")
Anova(lm(procdist ~ Watershed * simpson, data = Op4.procdist.df), type = "II")
Anova(lm(procdist ~ Watershed * LD1, data = Dent.procdist.df), type = "II")
Anova(lm(procdist ~ Watershed * simpson, data = Dent.procdist.df), type = "II")
Anova(lm(procdist ~ Watershed * LD1, data = Max4.procdist.df), type = "II")
Anova(lm(procdist ~ Watershed * simpson, data = Max4.procdist.df), type = "II")
Anova(lm(procdist ~ Watershed * LD1, data = tot.procdist.df), type = "II")
Anova(lm(procdist ~ Watershed * simpson, data = tot.procdist.df), type = "II")