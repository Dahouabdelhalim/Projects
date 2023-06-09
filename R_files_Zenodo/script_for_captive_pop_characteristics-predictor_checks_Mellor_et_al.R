#script for models making captive pop-predictor checks (see Table S3)
library(ape)
library(caper)
library(dplyr)

parrottree <- read.nexus("consensus_parrot_tree.nex")
parrot_data<-read.csv("parrot_comp_data.csv", header=TRUE)
parrot_comp_data <- comparative.data(phy=parrottree, data= parrot_data, 
                                     names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
attach(parrot_data)

#First set: running PGLS regressions with proportion adult as the outcome, 
#correlated against every  wild biology predictor in turn.
#starting with maximum foraging group size as the predictor

m1<-pgls(Prop_adult ~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m1)
par(mfrow=c(2,2))
#model diagnostics 
plot(m1)
#checking for normality of the phylogenetic residuals
res1<- residuals(m1, phylo = TRUE) #extracts phylogenetic residuals from the pgls model
shapiro.test(res1)
#now checking for outliers (phylogenetic residual +/- >3)
res1<- res1/sqrt(var(res1))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m1$residuals)#matches the residuals up with the species names 
rownames(res1)[(abs(res1)>3)]#gives the names of the outliers. 

#same again but with communal roosting as the predictor 
m2<-pgls(Prop_adult ~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m2)
par(mfrow=c(2,2))
plot(m2)
res2<- residuals(m2, phylo = TRUE) 
shapiro.test(res2)
res2<- res2/sqrt(var(res2))[1] 
rownames(res2)<-rownames(m2$residuals) 
rownames(res2)[(abs(res2)>3)]

#now with % diet needing extensive food search as the predictor. 
m3<-pgls(Prop_adult ~ Food_search, parrot_comp_data, lambda ='ML')
summary(m3)
par(mfrow=c(2,2))
#model diagnostics 
plot(m3)
res3<- residuals(m3, phylo = TRUE)
shapiro.test(res3)
res3<- res3/sqrt(var(res3))[1] 
rownames(res3)<-rownames(m3$residuals)
rownames(res3)[(abs(res3)>3)]

#now with food handling as the predictor
m4<-pgls(Prop_adult ~ sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m4)
par(mfrow=c(2,2))
plot(m4)
res4<- residuals(m4, phylo = TRUE) 
shapiro.test(res4)
res4<- res4/sqrt(var(res4))[1] 
rownames(res4)<-rownames(m4$residuals) 
rownames(res4)[(abs(res4)>3)]

#now with diet breadth as the predictor
m5<-pgls(Prop_adult ~ Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m5)
par(mfrow=c(2,2))
plot(m5)
res5<- residuals(m5, phylo = TRUE) 
shapiro.test(res5)
res5<- res5/sqrt(var(res5))[1] 
rownames(res5)<-rownames(m5$residuals) 
rownames(res5)[(abs(res5)>3)]

#now with habitat breadth as the predictor
m6<-pgls(Prop_adult ~ Habitat_breadth, parrot_comp_data,lambda ='ML')
summary(m6)
par(mfrow=c(2,2))
plot(m6)
res6<- residuals(m6, phylo = TRUE) 
shapiro.test(res6)
res6<- res6/sqrt(var(res6))[1] 
rownames(res6)<-rownames(m6$residuals) 
rownames(res6)[(abs(res6)>3)]

#now with innovation rate as the predictor. First have to filter the data so that only 
#species whose native range is within the geographical regions the innovation rate 
#database included. Note the model includes 'research effort' (n of published papers) 
#as a control (see MS)
parrot_data_inn<-filter(parrot_data, Inn_region==1)
#make a new comparative object using this filtered dataset
parrot_comp_data1 <- comparative.data(phy=parrottree, data= parrot_data_inn, 
                                      names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m7<-pgls(Prop_adult ~ Innovation+log(Research_effort), parrot_comp_data1, lambda ='ML')
summary(m7)
par(mfrow=c(2,2))
plot(m7)
res7<- residuals(m7, phylo = TRUE) 
shapiro.test(res7)
res7<- res7/sqrt(var(res7))[1] 
rownames(res7)<-rownames(m7$residuals) 
rownames(res7)[(abs(res7)>3)]

#now with brain volume as the predictor. Includes body mass as a control (for allometery; 
#see MS)
#Nandayus nenday was IDd as an outlier (phylogenetic residual +/->3) so removed from this model
new_parrot_data <-parrot_data[-c(125),]
new_parrot_comp_data <- comparative.data(phy=parrottree, data= new_parrot_data, 
                                         names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m8<-pgls(Prop_adult ~ log(Brain_vol)+log(Body_mass), new_parrot_comp_data, lambda ='ML')
summary(m8)
par(mfrow=c(2,2))
plot(m8)
res8<- residuals(m8, phylo = TRUE) 
shapiro.test(res8)
res8<- res8/sqrt(var(res8))[1] 
rownames(res8)<-rownames(m8$residuals) 
rownames(res8)[(abs(res8)>3)]

#now with IUCN status as the predictor
m9<-pgls(Prop_adult ~ IUCN_code, parrot_comp_data,lambda ='ML')
summary(m9)
par(mfrow=c(2,2))
plot(m9)
res9<- residuals(m9, phylo = TRUE) 
shapiro.test(res9)
res9<- res9/sqrt(var(res9))[1] 
rownames(res9)<-rownames(m9$residuals) 
rownames(res9)[(abs(res9)>3)]

#Second set: running PGLS regressions with proportion known sex as the outcome, correlated against
#every  wild biology predictor in turn.
#starting with maximum feeding group size as the predictor
m10<-pgls(log(Prop_known+1) ~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m10)
par(mfrow=c(2,2))
plot(m10)
res10<- residuals(m10, phylo = TRUE) 
shapiro.test(res10)
res10<- res10/sqrt(var(res10))[1] 
rownames(res10)<-rownames(m10$residuals) 
rownames(res10)[(abs(res10)>3)]

#same again but with communal roosting as the predictor 
m11<-pgls(log(Prop_known+1) ~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m11)
par(mfrow=c(2,2))
plot(m11)
res11<- residuals(m11, phylo = TRUE) 
shapiro.test(res11)
res11<- res11/sqrt(var(res11))[1] 
rownames(res11)<-rownames(m11$residuals) 
rownames(res11)[(abs(res11)>3)]

#now with % diet needing extensive food search as the predictor. 
m12<-pgls(Prop_known ~ Food_search, parrot_comp_data, lambda ='ML')
summary(m12)
par(mfrow=c(2,2))
#model diagnostics 
plot(m12)
res12<- residuals(m12, phylo = TRUE)
shapiro.test(res12)
res12<- res12/sqrt(var(res12))[1] 
rownames(res12)<-rownames(m12$residuals)
rownames(res12)[(abs(res12)>3)]

#now with food handling as the predictor
m13<-pgls(Prop_known ~ sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m13)
par(mfrow=c(2,2))
plot(m13)
res13<- residuals(m13, phylo = TRUE) 
shapiro.test(res13)
res13<- res13/sqrt(var(res13))[1] 
rownames(res13)<-rownames(m13$residuals) 
rownames(res13)[(abs(res13)>3)]

#now with diet breadth as the predictor
m14<-pgls(Prop_known ~ Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m14)
par(mfrow=c(2,2))
plot(m14)
res14<- residuals(m14, phylo = TRUE) 
shapiro.test(res14)
res14<- res14/sqrt(var(res14))[1] 
rownames(res14)<-rownames(m14$residuals) 
rownames(res14)[(abs(res14)>3)]

#now with habitat breadth as the predictor
m15<-pgls(Prop_known ~ Habitat_breadth, parrot_comp_data,lambda ='ML')
summary(m15)
par(mfrow=c(2,2))
plot(m15)
res15<- residuals(m15, phylo = TRUE) 
shapiro.test(res15)
res15<- res15/sqrt(var(res15))[1] 
rownames(res15)<-rownames(m15$residuals) 
rownames(res15)[(abs(res15)>3)]

#now with innovation rate as the predictor. First have to filter the data so that only 
#species whose native range is within the geographical regions the innovation rate 
#database included. Note the model includes 'research effort' (n of published papers) 
#as a control (see MS)
parrot_data_inn<-filter(parrot_data, Inn_region==1)
#make a new comparative object using this filtered dataset
parrot_comp_data1 <- comparative.data(phy=parrottree, data= parrot_data_inn, 
                                      names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m16<-pgls(Prop_known ~ Innovation+log(Research_effort), parrot_comp_data1, lambda ='ML')
summary(m16)
par(mfrow=c(2,2))
plot(m16)
res16<- residuals(m16, phylo = TRUE) 
shapiro.test(res16)
res16<- res16/sqrt(var(res16))[1] 
rownames(res16)<-rownames(m16$residuals) 
rownames(res16)[(abs(res16)>3)]

#now with brain volume as the predictor. Includes body mass as a control (for allometery; 
#see MS)
m17<-pgls(Prop_known ~ log(Brain_vol)+log(Body_mass), parrot_comp_data, lambda ='ML')
summary(m17)
par(mfrow=c(2,2))
plot(m17)
res17<- residuals(m17, phylo = TRUE) 
shapiro.test(res17)
res17<- res17/sqrt(var(res17))[1] 
rownames(res17)<-rownames(m17$residuals) 
rownames(res17)[(abs(res17)>3)]

#now with IUCN status as the predictor
m18<-pgls(Prop_known ~ IUCN_code, parrot_comp_data,lambda ='ML')
summary(m18)
par(mfrow=c(2,2))
plot(m18)
res18<- residuals(m18, phylo = TRUE) 
shapiro.test(res18)
res18<- res18/sqrt(var(res18))[1] 
rownames(res18)<-rownames(m18$residuals) 
rownames(res18)[(abs(res18)>3)]

#Third set: running PGLS regressions with proportion female as the outcome, correlated against
#every  wild biology predictor in turn.
#starting with maximum feeding group size as the predictor

m19<-pgls(sqrt(Prop_female) ~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m19)
par(mfrow=c(2,2))
plot(m19)
res19<- residuals(m19, phylo = TRUE) 
shapiro.test(res19)
res19<- res19/sqrt(var(res19))[1]
rownames(res19)<-rownames(m19$residuals)
rownames(res19)[(abs(res19)>3)]

#same again but with communal roosting as the predictor 
m20<-pgls(Prop_female ~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m20)
par(mfrow=c(2,2))
plot(m20)
res20<- residuals(m20, phylo = TRUE) 
shapiro.test(res20)
res20<- res20/sqrt(var(res20))[1] 
rownames(res20)<-rownames(m20$residuals) 
rownames(res20)[(abs(res20)>3)]

#now with % diet needing extensive food search as the predictor. 
m21<-pgls(Prop_female ~ Food_search, parrot_comp_data, lambda ='ML')
summary(m21)
par(mfrow=c(2,2))
#model diagnostics 
plot(m21)
res21<- residuals(m21, phylo = TRUE)
shapiro.test(res21)
res21<- res21/sqrt(var(res21))[1] 
rownames(res21)<-rownames(m21$residuals)
rownames(res21)[(abs(res21)>3)]

#now with food handling as the predictor
m22<-pgls(Prop_female~ sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m22)
par(mfrow=c(2,2))
plot(m22)
res22<- residuals(m22, phylo = TRUE) 
shapiro.test(res22)
res22<- res22/sqrt(var(res22))[1] 
rownames(res22)<-rownames(m22$residuals) 
rownames(res22)[(abs(res22)>3)]

#now with diet breadth as the predictor
m23<-pgls(log(Prop_female+1) ~ Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m23)
par(mfrow=c(2,2))
plot(m23)
res23<- residuals(m23, phylo = TRUE) 
shapiro.test(res23)
res23<- res23/sqrt(var(res23))[1] 
rownames(res23)<-rownames(m23$residuals) 
rownames(res23)[(abs(res23)>3)]

#now with habitat breadth as the predictor
m24<-pgls(Prop_female ~ Habitat_breadth, parrot_comp_data,lambda ='ML')
summary(m24)
par(mfrow=c(2,2))
plot(m24)
res24<- residuals(m24, phylo = TRUE) 
shapiro.test(res24)
res24<- res24/sqrt(var(res24))[1] 
rownames(res24)<-rownames(m24$residuals) 
rownames(res24)[(abs(res24)>3)]

#now with innovation rate as the predictor. First have to filter the data so that only 
#species whose native range is within the geographical regions the innovation rate 
#database included. Note the model includes 'research effort' (n of published papers) 
#as a control (see MS)
parrot_data_inn<-filter(parrot_data, Inn_region==1)
#make a new comparative object using this filtered dataset
parrot_comp_data1 <- comparative.data(phy=parrottree, data= parrot_data_inn, 
                                      names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m25<-pgls(Prop_female ~ Innovation+log(Research_effort), parrot_comp_data1, lambda ='ML')
summary(m25)
par(mfrow=c(2,2))
plot(m25)
res25<- residuals(m25, phylo = TRUE) 
shapiro.test(res25)
res25<- res25/sqrt(var(res25))[1] 
rownames(res25)<-rownames(m25$residuals) 
rownames(res25)[(abs(res25)>3)]

#now with brain volume as the predictor. Includes body mass as a control (for allometery; 
#see MS)
#Forpus coelestis was not formally IDd as an outlier, but its inclusion made residuals non-normal 
#so removed here
new_parrot_data <-parrot_data[-c(106),]
new_parrot_comp_data <- comparative.data(phy=parrottree, data= new_parrot_data, 
                                         names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m26<-pgls(Prop_female ~ log(Brain_vol)+log(Body_mass), new_parrot_comp_data, lambda ='ML')
summary(m26)
par(mfrow=c(2,2))
plot(m26)
res26<- residuals(m26, phylo = TRUE) 
shapiro.test(res26)
res26<- res26/sqrt(var(res26))[1] 
rownames(res26)<-rownames(m26$residuals) 
rownames(res26)[(abs(res26)>3)]

#now with IUCN status as the predictor
m27<-pgls(Prop_female ~ IUCN_code, parrot_comp_data,lambda ='ML')
summary(m27)
par(mfrow=c(2,2))
plot(m27)
res27<- residuals(m27, phylo = TRUE) 
shapiro.test(res27)
res27<- res27/sqrt(var(res27))[1] 
rownames(res27)<-rownames(m27$residuals) 
rownames(res27)[(abs(res27)>3)]

#Fourth set: running PGLS regressions with proportion hand-reared as the outcome, correlated against
#every  wild biology predictor in turn.
#starting with maximum feeding group size as the predictor

m28<-pgls(sqrt(Human_reared) ~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m28)
par(mfrow=c(2,2))
plot(m28)
res28<- residuals(m28, phylo = TRUE) 
shapiro.test(res28)
res28<- res28/sqrt(var(res28))[1]
rownames(res28)<-rownames(m28$residuals)
rownames(res28)[(abs(res28)>3)]

#same again but with communal roosting as the predictor 
m29<-pgls(sqrt(Human_reared) ~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m29)
par(mfrow=c(2,2))
plot(m29)
res29<- residuals(m29, phylo = TRUE) 
shapiro.test(res29)
res29<- res29/sqrt(var(res29))[1] 
rownames(res29)<-rownames(m29$residuals) 
rownames(res29)[(abs(res29)>3)]

#now with % diet needing extensive food search as the predictor. 
m30<-pgls(sqrt(Human_reared) ~ Food_search, parrot_comp_data, lambda ='ML')
summary(m30)
par(mfrow=c(2,2))
#model diagnostics 
plot(m30)
res30<- residuals(m30, phylo = TRUE)
shapiro.test(res30)
res30<- res30/sqrt(var(res30))[1] 
rownames(res30)<-rownames(m30$residuals)
rownames(res30)[(abs(res30)>3)]

#now with food handling as the predictor
m31<-pgls(sqrt(Human_reared)~ sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m31)
par(mfrow=c(2,2))
plot(m31)
res31<- residuals(m31, phylo = TRUE) 
shapiro.test(res31)
res31<- res31/sqrt(var(res31))[1] 
rownames(res31)<-rownames(m31$residuals) 
rownames(res31)[(abs(res31)>3)]

#now with diet breadth as the predictor
m32<-pgls(sqrt(Human_reared) ~ Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m32)
par(mfrow=c(2,2))
plot(m32)
res32<- residuals(m32, phylo = TRUE) 
shapiro.test(res32)
res32<- res32/sqrt(var(res32))[1] 
rownames(res32)<-rownames(m32$residuals) 
rownames(res32)[(abs(res32)>3)]

#now with habitat breadth as the predictor
m33<-pgls(sqrt(Human_reared) ~ Habitat_breadth, parrot_comp_data,lambda ='ML')
summary(m33)
par(mfrow=c(2,2))
plot(m33)
res33<- residuals(m33, phylo = TRUE) 
shapiro.test(res33)
res33<- res33/sqrt(var(res33))[1] 
rownames(res33)<-rownames(m33$residuals) 
rownames(res33)[(abs(res33)>3)]

#now with innovation rate as the predictor. First have to filter the data so that only 
#species whose native range is within the geographical regions the innovation rate 
#database included. Note the model includes 'research effort' (n of published papers) 
#as a control (see MS)
parrot_data_inn<-filter(parrot_data, Inn_region==1)
#make a new comparative object using this filtered dataset
parrot_comp_data1 <- comparative.data(phy=parrottree, data= parrot_data_inn, 
                                      names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m34<-pgls(sqrt(Human_reared) ~ Innovation+log(Research_effort), parrot_comp_data1, lambda ='ML')
summary(m34)
par(mfrow=c(2,2))
plot(m34)
res34<- residuals(m34, phylo = TRUE) 
shapiro.test(res34)
res34<- res34/sqrt(var(res34))[1] 
rownames(res34)<-rownames(m34$residuals) 
rownames(res34)[(abs(res34)>3)]

#now with brain volume as the predictor. Includes body mass as a control (for allometery; 
#see MS)
m35<-pgls(Human_reared ~ log(Brain_vol)+log(Body_mass), parrot_comp_data, lambda ='ML')
summary(m35)
par(mfrow=c(2,2))
plot(m35)
res35<- residuals(m35, phylo = TRUE) 
shapiro.test(res35)
res35<- res35/sqrt(var(res35))[1] 
rownames(res35)<-rownames(m35$residuals) 
rownames(res35)[(abs(res35)>3)]

#now with IUCN status as the predictor
m36<-pgls(sqrt(Human_reared) ~ IUCN_code, parrot_comp_data,lambda ='ML')
summary(m36)
par(mfrow=c(2,2))
plot(m36)
res36<- residuals(m36, phylo = TRUE) 
shapiro.test(res36)
res36<- res36/sqrt(var(res36))[1] 
rownames(res36)<-rownames(m36$residuals) 
rownames(res36)[(abs(res36)>3)]

#Fifth set: running PGLS regressions with proportion housed in a standard-sized
#cage as the outcome, correlated against every  wild biology predictor in turn.
#starting with maximum feeding group size as the predictor
m37<-pgls(sqrt(Stand_cage)~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m37)
par(mfrow=c(2,2))
plot(m37)
res37<- residuals(m37, phylo = TRUE) 
shapiro.test(res37)
res37<- res37/sqrt(var(res37))[1]
rownames(res37)<-rownames(m37$residuals)
rownames(res37)[(abs(res37)>3)]

#same again but with communal roosting as the predictor 
m38<-pgls(log(Stand_cage+1) ~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m38)
par(mfrow=c(2,2))
plot(m38)
res38<- residuals(m38, phylo = TRUE) 
shapiro.test(res38)
res38<- res38/sqrt(var(res38))[1] 
rownames(res38)<-rownames(m38$residuals) 
rownames(res38)[(abs(res38)>3)]

#now with % diet needing extensive food search as the predictor. 
#Diopsittaca nobilis an outlier, but its inclusion made residuals non-normal 
#so removed here
new_parrot_data <-parrot_data[-c(96),]
new_parrot_comp_data <- comparative.data(phy=parrottree, data= new_parrot_data, 
                                         names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m39<-pgls(sqrt(Stand_cage) ~ Food_search, new_parrot_comp_data, lambda ='ML')
summary(m39)
par(mfrow=c(2,2))
#model diagnostics 
plot(m39)
res39<- residuals(m39, phylo = TRUE)
shapiro.test(res39)
res39<- res39/sqrt(var(res39))[1] 
rownames(res39)<-rownames(m39$residuals)
rownames(res39)[(abs(res39)>3)]

#now with food handling as the predictor
m40<-pgls(sqrt(Stand_cage)~ sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m40)
par(mfrow=c(2,2))
plot(m40)
res40<- residuals(m40, phylo = TRUE) 
shapiro.test(res40)
res40<- res40/sqrt(var(res40))[1] 
rownames(res40)<-rownames(m40$residuals) 
rownames(res40)[(abs(res40)>3)]

#now with diet breadth as the predictor
m41<-pgls(sqrt(Stand_cage) ~ Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m41)
par(mfrow=c(2,2))
plot(m41)
res41<- residuals(m41, phylo = TRUE) 
shapiro.test(res41)
res41<- res41/sqrt(var(res41))[1] 
rownames(res41)<-rownames(m41$residuals) 
rownames(res41)[(abs(res41)>3)]

#now with habitat breadth as the predictor
m42<-pgls(Stand_cage ~ Habitat_breadth, parrot_comp_data,lambda ='ML')
summary(m42)
par(mfrow=c(2,2))
plot(m42)
res42<- residuals(m42, phylo = TRUE) 
shapiro.test(res42)
res42<- res42/sqrt(var(res42))[1] 
rownames(res42)<-rownames(m42$residuals) 
rownames(res42)[(abs(res42)>3)]

#now with innovation rate as the predictor. First have to filter the data so that only 
#species whose native range is within the geographical regions the innovation rate 
#database included. Note the model includes 'research effort' (n of published papers) 
#as a control (see MS)
parrot_data_inn<-filter(parrot_data, Inn_region==1)
#make a new comparative object using this filtered dataset
parrot_comp_data1 <- comparative.data(phy=parrottree, data= parrot_data_inn, 
                                      names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m43<-pgls(sqrt(Stand_cage) ~ Innovation+log(Research_effort), parrot_comp_data1, lambda ='ML')
summary(m43)
par(mfrow=c(2,2))
plot(m43)
res43<- residuals(m43, phylo = TRUE) 
shapiro.test(res43)
res43<- res43/sqrt(var(res43))[1] 
rownames(res43)<-rownames(m43$residuals) 
rownames(res43)[(abs(res43)>3)]

#now with brain volume as the predictor. Includes body mass as a control (for allometery; 
#see MS)
m44<-pgls(log(Stand_cage+1) ~ log(Brain_vol)+log(Body_mass), parrot_comp_data, lambda ='ML')
summary(m44)
par(mfrow=c(2,2))
plot(m44)
res44<- residuals(m44, phylo = TRUE) 
shapiro.test(res44)
res44<- res44/sqrt(var(res44))[1] 
rownames(res44)<-rownames(m44$residuals) 
rownames(res44)[(abs(res44)>3)]

#now with IUCN status as the predictor
m45<-pgls(sqrt(Stand_cage) ~ IUCN_code, parrot_comp_data,lambda ='ML')
summary(m45)
par(mfrow=c(2,2))
plot(m45)
res45<- residuals(m45, phylo = TRUE) 
shapiro.test(res45)
res45<- res45/sqrt(var(res45))[1] 
rownames(res45)<-rownames(m45$residuals) 
rownames(res45)[(abs(res45)>3)]

#Sixth set: running PGLS regressions with proportion socially isolated as the outcome, correlated against
#every  wild biology predictor in turn.
#starting with maximum feeding group size as the predictor

m46<-pgls(Prop_isolated~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m46)
par(mfrow=c(2,2))
plot(m46)
res46<- residuals(m46, phylo = TRUE) 
shapiro.test(res46)
res46<- res46/sqrt(var(res46))[1]
rownames(res46)<-rownames(m46$residuals)
rownames(res46)[(abs(res46)>3)]

#same again but with communal roosting as the predictor 
m47<-pgls(Prop_isolated ~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m47)
par(mfrow=c(2,2))
plot(m47)
res47<- residuals(m47, phylo = TRUE) 
shapiro.test(res47)
res47<- res47/sqrt(var(res47))[1] 
rownames(res47)<-rownames(m47$residuals) 
rownames(res47)[(abs(res47)>3)]

#now with food search as the predictor
m48<-pgls(Prop_isolated ~ Food_search, parrot_comp_data, lambda ='ML')
summary(m48)
par(mfrow=c(2,2))
#model diagnostics 
plot(m48)
res48<- residuals(m48, phylo = TRUE)
shapiro.test(res48)
res48<- res48/sqrt(var(res48))[1] 
rownames(res48)<-rownames(m48$residuals)
rownames(res48)[(abs(res48)>3)]

#now with food handling as the predictor
m49<-pgls(Prop_isolated~ sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m49)
par(mfrow=c(2,2))
plot(m49)
res49<- residuals(m49, phylo = TRUE) 
shapiro.test(res49)
res49<- res49/sqrt(var(res49))[1] 
rownames(res49)<-rownames(m49$residuals) 
rownames(res49)[(abs(res49)>3)]

#now with diet breadth as the predictor
m50<-pgls(Prop_isolated ~ Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m50)
par(mfrow=c(2,2))
plot(m50)
res50<- residuals(m50, phylo = TRUE) 
shapiro.test(res50)
res50<- res50/sqrt(var(res50))[1] 
rownames(res50)<-rownames(m50$residuals) 
rownames(res50)[(abs(res50)>3)]

#now with habitat breadth as the predictor
m51<-pgls(Prop_isolated ~ Habitat_breadth, parrot_comp_data,lambda ='ML')
summary(m51)
par(mfrow=c(2,2))
plot(m51)
res51<- residuals(m51, phylo = TRUE) 
shapiro.test(res51)
res51<- res51/sqrt(var(res51))[1] 
rownames(res51)<-rownames(m51$residuals) 
rownames(res51)[(abs(res51)>3)]

#now with innovation rate as the predictor. First have to filter the data so that only 
#species whose native range is within the geographical regions the innovation rate 
#database included. Note the model includes 'research effort' (n of published papers) 
#as a control (see MS)
parrot_data_inn<-filter(parrot_data, Inn_region==1)
#make a new comparative object using this filtered dataset
parrot_comp_data1 <- comparative.data(phy=parrottree, data= parrot_data_inn, 
                                      names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m52<-pgls(Prop_isolated ~ Innovation+log(Research_effort), parrot_comp_data1, lambda ='ML')
summary(m52)
par(mfrow=c(2,2))
plot(m52)
res52<- residuals(m52, phylo = TRUE) 
shapiro.test(res52)
res52<- res52/sqrt(var(res52))[1] 
rownames(res52)<-rownames(m52$residuals) 
rownames(res52)[(abs(res52)>3)]

#now with brain volume as the predictor. Includes body mass as a control (for allometery; 
#see MS)
m53<-pgls(Prop_isolated ~ log(Brain_vol)+log(Body_mass), parrot_comp_data, lambda ='ML')
summary(m53)
par(mfrow=c(2,2))
plot(m53)
res53<- residuals(m53, phylo = TRUE) 
shapiro.test(res53)
res53<- res53/sqrt(var(res53))[1] 
rownames(res53)<-rownames(m53$residuals) 
rownames(res53)[(abs(res53)>3)]

#now with IUCN status as the predictor
m54<-pgls(Prop_isolated ~ IUCN_code, parrot_comp_data,lambda ='ML')
summary(m54)
par(mfrow=c(2,2))
plot(m54)
res54<- residuals(m54, phylo = TRUE) 
shapiro.test(res54)
res54<- res54/sqrt(var(res54))[1] 
rownames(res54)<-rownames(m54$residuals) 
rownames(res54)[(abs(res54)>3)]

#Seventh set: running PGLS regressions with proportion with short feed times as the outcome, 
#correlated against every wild biology predictor in turn.
#starting with maximum feeding group size as the predictor

#Cacatua sanguinea IDd as an outlier so removed here
new_parrot_data <-parrot_data[-c(73),]
new_parrot_comp_data <- comparative.data(phy=parrottree, data= new_parrot_data, 
                                         names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m55<-pgls(Short_feed~ log(Max_feed_size), new_parrot_comp_data, lambda ='ML')
summary(m55)
par(mfrow=c(2,2))
plot(m55)
res55<- residuals(m55, phylo = TRUE) 
shapiro.test(res55)
res55<- res55/sqrt(var(res55))[1]
rownames(res55)<-rownames(m55$residuals)
rownames(res55)[(abs(res55)>3)]

#same again but with communal roosting as the predictor 
m56<-pgls(Short_feed ~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m56)
par(mfrow=c(2,2))
plot(m56)
res56<- residuals(m56, phylo = TRUE) 
shapiro.test(res56)
res56<- res56/sqrt(var(res56))[1] 
rownames(res56)<-rownames(m56$residuals) 
rownames(res56)[(abs(res56)>3)]

#now with food search as the predictor
m57<-pgls(Short_feed ~ Food_search, parrot_comp_data, lambda ='ML')
summary(m57)
par(mfrow=c(2,2))
#model diagnostics 
plot(m57)
res57<- residuals(m57, phylo = TRUE)
shapiro.test(res57)
res57<- res57/sqrt(var(res57))[1] 
rownames(res57)<-rownames(m57$residuals)
rownames(res57)[(abs(res57)>3)]

#now with food handling as the predictor
m58<-pgls(Short_feed~ sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m58)
par(mfrow=c(2,2))
plot(m58)
res58<- residuals(m58, phylo = TRUE) 
shapiro.test(res58)
res58<- res58/sqrt(var(res58))[1] 
rownames(res58)<-rownames(m58$residuals) 
rownames(res58)[(abs(res58)>3)]

#now with diet breadth as the predictor
m59<-pgls(Short_feed ~ Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m59)
par(mfrow=c(2,2))
plot(m59)
res59<- residuals(m59, phylo = TRUE) 
shapiro.test(res59)
res59<- res59/sqrt(var(res59))[1] 
rownames(res59)<-rownames(m59$residuals) 
rownames(res59)[(abs(res59)>3)]

#now with habitat breadth as the predictor
m60<-pgls(Short_feed ~ Habitat_breadth, parrot_comp_data,lambda ='ML')
summary(m60)
par(mfrow=c(2,2))
plot(m60)
res60<- residuals(m60, phylo = TRUE) 
shapiro.test(res60)
res60<- res60/sqrt(var(res60))[1] 
rownames(res60)<-rownames(m60$residuals) 
rownames(res60)[(abs(res60)>3)]

#now with innovation rate as the predictor. First have to filter the data so that only 
#species whose native range is within the geographical regions the innovation rate 
#database included. Note the model includes 'research effort' (n of published papers) 
#as a control (see MS)
parrot_data_inn<-filter(parrot_data, Inn_region==1)
#make a new comparative object using this filtered dataset
parrot_comp_data1 <- comparative.data(phy=parrottree, data= parrot_data_inn, 
                                      names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m61<-pgls(Short_feed ~ Innovation+log(Research_effort), parrot_comp_data1, lambda ='ML')
summary(m61)
par(mfrow=c(2,2))
plot(m61)
res61<- residuals(m61, phylo = TRUE) 
shapiro.test(res61)
res61<- res61/sqrt(var(res61))[1] 
rownames(res61)<-rownames(m61$residuals) 
rownames(res61)[(abs(res61)>3)]

#now with brain volume as the predictor. Includes body mass as a control (for allometery; 
#see MS)
m62<-pgls(Short_feed ~ log(Brain_vol)+log(Body_mass), parrot_comp_data, lambda ='ML')
summary(m62)
par(mfrow=c(2,2))
plot(m62)
res62<- residuals(m62, phylo = TRUE) 
shapiro.test(res62)
res62<- res62/sqrt(var(res62))[1] 
rownames(res62)<-rownames(m62$residuals) 
rownames(res62)[(abs(res62)>3)]

#now with IUCN status as the predictor
m63<-pgls(Short_feed ~ IUCN_code, parrot_comp_data,lambda ='ML')
summary(m63)
par(mfrow=c(2,2))
plot(m63)
res63<- residuals(m63, phylo = TRUE) 
shapiro.test(res63)
res63<- res63/sqrt(var(res63))[1] 
rownames(res63)<-rownames(m63$residuals) 
rownames(res63)[(abs(res63)>3)]

#Eighth set: running PGLS regressions with median diet diversity as the outcome, 
#correlated against every  wild biology predictor in turn.
#starting with maximum feeding group size as the predictor

m64<-pgls(Cap_diet_div~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m64)
par(mfrow=c(2,2))
plot(m64)
res64<- residuals(m64, phylo = TRUE) 
shapiro.test(res64)
res64<- res64/sqrt(var(res64))[1]
rownames(res64)<-rownames(m64$residuals)
rownames(res64)[(abs(res64)>3)]

#same again but with communal roosting as the predictor 
m65<-pgls(Cap_diet_div ~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m65)
par(mfrow=c(2,2))
plot(m65)
res65<- residuals(m65, phylo = TRUE) 
shapiro.test(res65)
res65<- res65/sqrt(var(res65))[1] 
rownames(res65)<-rownames(m65$residuals) 
rownames(res65)[(abs(res65)>3)]

#now with food search as the predictor
m66<-pgls(Cap_diet_div ~ Food_search, parrot_comp_data, lambda ='ML')
summary(m66)
par(mfrow=c(2,2))
#model diagnostics 
plot(m66)
res66<- residuals(m66, phylo = TRUE)
shapiro.test(res66)
res66<- res66/sqrt(var(res66))[1] 
rownames(res66)<-rownames(m66$residuals)
rownames(res66)[(abs(res66)>3)]

#now with food handling as the predictor
m67<-pgls(Cap_diet_div~ sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m67)
par(mfrow=c(2,2))
plot(m67)
res67<- residuals(m67, phylo = TRUE) 
shapiro.test(res67)
res67<- res67/sqrt(var(res67))[1] 
rownames(res67)<-rownames(m67$residuals) 
rownames(res67)[(abs(res67)>3)]

#now with diet breadth as the predictor
m68<-pgls(log(Cap_diet_div) ~ Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m68)
par(mfrow=c(2,2))
plot(m68)
res68<- residuals(m68, phylo = TRUE) 
shapiro.test(res68)
res68<- res68/sqrt(var(res68))[1] 
rownames(res68)<-rownames(m68$residuals) 
rownames(res68)[(abs(res68)>3)]

#now with habitat breadth as the predictor
m69<-pgls(Cap_diet_div ~ Habitat_breadth, parrot_comp_data,lambda ='ML')
summary(m69)
par(mfrow=c(2,2))
plot(m69)
res69<-residuals(m69, phylo = TRUE) 
shapiro.test(res69)
res69<- res69/sqrt(var(res69))[1] 
rownames(res69)<-rownames(m69$residuals) 
rownames(res69)[(abs(res69)>3)]

#now with innovation rate as the predictor. First have to filter the data so that only 
#species whose native range is within the geographical regions the innovation rate 
#database included. Note the model includes 'research effort' (n of published papers) 
#as a control (see MS)
parrot_data_inn<-filter(parrot_data, Inn_region==1)
#make a new comparative object using this filtered dataset
parrot_comp_data1 <- comparative.data(phy=parrottree, data= parrot_data_inn, 
                                      names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m70<-pgls(log(Cap_diet_div) ~ Innovation+log(Research_effort), parrot_comp_data1, lambda ='ML')
summary(m70)
par(mfrow=c(2,2))
plot(m70)
res70<-residuals(m70, phylo = TRUE) 
shapiro.test(res70)
res70<-res70/sqrt(var(res70))[1] 
rownames(res70)<-rownames(m70$residuals) 
rownames(res70)[(abs(res70)>3)]

#now with brain volume as the predictor. Includes body mass as a control (for allometery; 
#see MS)
m71<-pgls(Cap_diet_div ~ log(Brain_vol)+log(Body_mass), parrot_comp_data, lambda ='ML')
summary(m71)
par(mfrow=c(2,2))
plot(m71)
res71<- residuals(m71, phylo = TRUE) 
shapiro.test(res71)
res71<-res71/sqrt(var(res71))[1] 
rownames(res71)<-rownames(m71$residuals) 
rownames(res71)[(abs(res71)>3)]

#now with IUCN status as the predictor
m72<-pgls(Cap_diet_div ~ IUCN_code, parrot_comp_data,lambda ='ML')
summary(m72)
par(mfrow=c(2,2))
plot(m72)
res72<- residuals(m72, phylo = TRUE) 
shapiro.test(res72)
res72<- res72/sqrt(var(res72))[1] 
rownames(res72)<-rownames(m72$residuals) 
rownames(res72)[(abs(res72)>3)]

#Ninth set: running PGLS regressions with median early enrichment diversity as the outcome, 
#correlated against every  wild biology predictor in turn.
#starting with maximum feeding group size as the predictor

m73<-pgls(log(Early_EE+1)~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m73)
par(mfrow=c(2,2))
plot(m73)
res73<- residuals(m73, phylo = TRUE) 
shapiro.test(res73)
res73<- res73/sqrt(var(res73))[1]
rownames(res73)<-rownames(m73$residuals)
rownames(res73)[(abs(res73)>3)]

#same again but with communal roosting as the predictor. This one 
#only just passes a normality test; however, an outlier(Anodorhynchus_
#hyacinthinus pops up when transformations are applied, and further 
#ones pop up after its removal. In the interest of retaining sample 
#size (and these not being focal analyses), this one the best choice 
m74<-pgls(Early_EE~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m74)
par(mfrow=c(2,2))
plot(m74)
res74<- residuals(m74, phylo = TRUE) 
shapiro.test(res74)
res74<- res74/sqrt(var(res74))[1] 
rownames(res74)<-rownames(m74$residuals) 
rownames(res74)[(abs(res74)>3)]

#now with food search as the predictor
m75<-pgls(Early_EE ~ Food_search, parrot_comp_data, lambda ='ML')
summary(m75)
par(mfrow=c(2,2))
#model diagnostics 
plot(m75)
res75<- residuals(m75, phylo = TRUE)
shapiro.test(res75)
res75<- res75/sqrt(var(res75))[1] 
rownames(res75)<-rownames(m75$residuals)
rownames(res75)[(abs(res75)>3)]

#now with food handling as the predictor
m76<-pgls(Early_EE~ sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m76)
par(mfrow=c(2,2))
plot(m76)
res76<- residuals(m76, phylo = TRUE) 
shapiro.test(res76)
res76<- res76/sqrt(var(res76))[1] 
rownames(res76)<-rownames(m76$residuals) 
rownames(res76)[(abs(res76)>3)]

#now with diet breadth as the predictor
m77<-pgls(Early_EE ~ Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m77)
par(mfrow=c(2,2))
plot(m77)
res77<-residuals(m77, phylo = TRUE) 
shapiro.test(res77)
res77<-res77/sqrt(var(res77))[1] 
rownames(res77)<-rownames(m77$residuals) 
rownames(res77)[(abs(res77)>3)]

#now with habitat breadth as the predictor
m78<-pgls(Early_EE ~ Habitat_breadth, parrot_comp_data,lambda ='ML')
summary(m78)
par(mfrow=c(2,2))
plot(m78)
res78<-residuals(m78, phylo = TRUE) 
shapiro.test(res78)
res78<- res78/sqrt(var(res78))[1] 
rownames(res78)<-rownames(m78$residuals) 
rownames(res78)[(abs(res78)>3)]

#now with innovation rate as the predictor. First have to filter the data so that only 
#species whose native range is within the geographical regions the innovation rate 
#database included. Note the model includes 'research effort' (n of published papers) 
#as a control (see MS)
parrot_data_inn<-filter(parrot_data, Inn_region==1)
#make a new comparative object using this filtered dataset
parrot_comp_data1 <- comparative.data(phy=parrottree, data= parrot_data_inn, 
                                      names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m79<-pgls(Early_EE ~ Innovation+log(Research_effort), parrot_comp_data1, lambda ='ML')
summary(m79)
par(mfrow=c(2,2))
plot(m79)
res79<-residuals(m79, phylo = TRUE) 
shapiro.test(res79)
res79<-res79/sqrt(var(res79))[1] 
rownames(res79)<-rownames(m79$residuals) 
rownames(res79)[(abs(res79)>3)]

#now with brain volume as the predictor. Includes body mass as a control (for allometery; 
#see MS)
m80<-pgls(Early_EE ~ log(Brain_vol)+log(Body_mass), parrot_comp_data, lambda ='ML')
summary(m80)
par(mfrow=c(2,2))
plot(m80)
res80<- residuals(m80, phylo = TRUE) 
shapiro.test(res80)
res80<-res80/sqrt(var(res80))[1] 
rownames(res80)<-rownames(m80$residuals) 
rownames(res80)[(abs(res80)>3)]

#now with IUCN status as the predictor
m81<-pgls(log(Early_EE+1) ~ IUCN_code, parrot_comp_data,lambda ='ML')
summary(m81)
par(mfrow=c(2,2))
plot(m81)
res81<- residuals(m81, phylo = TRUE) 
shapiro.test(res81)
res81<- res81/sqrt(var(res81))[1] 
rownames(res81)<-rownames(m81$residuals) 
rownames(res81)[(abs(res81)>3)]

#Tenth set: running PGLS regressions with median current enrichment diversity as the outcome, 
#correlated against every  wild biology predictor in turn.
#starting with maximum feeding group size as the predictor

m82<-pgls(Current_EE~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m82)
par(mfrow=c(2,2))
plot(m82)
res82<-residuals(m82, phylo = TRUE) 
shapiro.test(res82)
res82<- res82/sqrt(var(res82))[1]
rownames(res82)<-rownames(m82$residuals)
rownames(res82)[(abs(res82)>3)]

#same again but with communal roosting as the predictor.
m83<-pgls(Current_EE~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m83)
par(mfrow=c(2,2))
plot(m83)
res83<- residuals(m83, phylo = TRUE) 
shapiro.test(res83)
res83<- res83/sqrt(var(res83))[1] 
rownames(res83)<-rownames(m83$residuals) 
rownames(res83)[(abs(res83)>3)]

#Myiopsitta_monachus was IDd as an outlier and residuals only just pass
#a normality test. However, its removal results in non-normal residuals
#and other outliers pop up. So, retaining this one in the model
#now with food search as the predictor
#Myiopsitta_monachus and Nandayus nenday was IDd as outliers (phylogenetic residual +/->3) 
#and their inclusion affected normality of the residuals, so removed from this model
new_parrot_data <-parrot_data[-c(124, 125),]
new_parrot_comp_data <- comparative.data(phy=parrottree, data= new_parrot_data, 
                                         names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m84<-pgls(Current_EE ~Food_search, new_parrot_comp_data, lambda ='ML')
summary(m84)
par(mfrow=c(2,2))
#model diagnostics 
plot(m84)
res84<- residuals(m84, phylo = TRUE)
shapiro.test(res84)
res84<- res84/sqrt(var(res84))[1] 
rownames(res84)<-rownames(m84$residuals)
rownames(res84)[(abs(res84)>3)]

#now with food handling as the predictor
#Myiopsitta_monachus and Nandayus nenday was IDd as outliers (phylogenetic residual +/->3) 
#and their inclusion affected normality of the residuals, so removed from this model

m85<-pgls(Current_EE~ sqrt(Food_handling), new_parrot_comp_data, lambda ='ML')
summary(m85)
par(mfrow=c(2,2))
plot(m85)
res85<- residuals(m85, phylo = TRUE) 
shapiro.test(res85)
res85<- res85/sqrt(var(res85))[1] 
rownames(res85)<-rownames(m85$residuals) 
rownames(res85)[(abs(res85)>3)]

#now with diet breadth as the predictor
#Myiopsitta_monachus and Nandayus nenday was IDd as outliers (phylogenetic residual +/->3) 
#and their inclusion affected normality of the residuals, so removed from this model
m86<-pgls(Current_EE ~ Diet_breadth, new_parrot_comp_data, lambda ='ML')
summary(m86)
par(mfrow=c(2,2))
plot(m86)
res86<-residuals(m86, phylo = TRUE) 
shapiro.test(res86)
res86<-res86/sqrt(var(res86))[1] 
rownames(res86)<-rownames(m86$residuals) 
rownames(res86)[(abs(res86)>3)]

#now with habitat breadth as the predictor
#Myiopsitta_monachus and Nandayus nenday was IDd as outliers (phylogenetic residual +/->3) 
#and their inclusion affected normality of the residuals, so removed from this model
m87<-pgls(Current_EE ~ Habitat_breadth, new_parrot_comp_data,lambda ='ML')
summary(m87)
par(mfrow=c(2,2))
plot(m87)
res87<-residuals(m87, phylo = TRUE) 
shapiro.test(res87)
res87<- res87/sqrt(var(res87))[1] 
rownames(res87)<-rownames(m87$residuals) 
rownames(res87)[(abs(res87)>3)]

#now with innovation rate as the predictor. First have to filter the data so that only 
#species whose native range is within the geographical regions the innovation rate 
#database included. Note the model includes 'research effort' (n of published papers) 
#as a control (see MS)
parrot_data_inn<-filter(parrot_data, Inn_region==1)
#make a new comparative object using this filtered dataset
parrot_comp_data1 <- comparative.data(phy=parrottree, data= parrot_data_inn, 
                                      names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m88<-pgls(Current_EE ~ Innovation+log(Research_effort), parrot_comp_data1, lambda ='ML')
summary(m88)
par(mfrow=c(2,2))
plot(m88)
res88<-residuals(m88, phylo = TRUE) 
shapiro.test(res88)
res88<-res88/sqrt(var(res88))[1] 
rownames(res88)<-rownames(m88$residuals) 
rownames(res88)[(abs(res88)>3)]

#now with brain volume as the predictor. Includes body mass as a control (for allometery; 
#see MS)
m89<-pgls(Current_EE ~ log(Brain_vol)+log(Body_mass), parrot_comp_data, lambda ='ML')
summary(m89)
par(mfrow=c(2,2))
plot(m89)
res89<- residuals(m89, phylo = TRUE) 
shapiro.test(res89)
res89<-res89/sqrt(var(res89))[1] 
rownames(res89)<-rownames(m89$residuals) 
rownames(res89)[(abs(res89)>3)]

#now with IUCN status as the predictor
#Myiopsitta_monachus and Nandayus nenday was IDd as outliers (phylogenetic residual +/->3) 
#and their inclusion affected normality of the residuals, so removed from this model
m90<-pgls(Current_EE ~ IUCN_code, new_parrot_comp_data,lambda ='ML')
summary(m90)
par(mfrow=c(2,2))
plot(m90)
res90<- residuals(m90, phylo = TRUE) 
shapiro.test(res90)
res90<-res90/sqrt(var(res90))[1] 
rownames(res90)<-rownames(m90$residuals) 
rownames(res90)[(abs(res90)>3)]