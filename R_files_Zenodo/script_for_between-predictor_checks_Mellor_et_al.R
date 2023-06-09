#script for models making between-predictor checks (see Table S5)
library(ape)
library(caper)
library(dplyr)
library(phylolm)
library(geiger)

parrottree <- read.nexus("consensus_parrot_tree.nex")
parrot_data<-read.csv("parrot_comp_data.csv", header=TRUE)
parrot_comp_data <- comparative.data(phy=parrottree, data= parrot_data, 
                                     names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
attach(parrot_data)

#First set: running PGLS regressions with maximum feeding group size as the outcome, 
#correlated against every other wild biology predictor in turn.
#starting with communal roosting as the predictor

m1<-pgls(log(Max_feed_size) ~ Communal_roost, parrot_comp_data, lambda ='ML')
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

#same again but with % diet needing extensive search as the predictor 
m2<-pgls(log(Max_feed_size) ~ Food_search, parrot_comp_data, lambda ='ML')
summary(m2)
par(mfrow=c(2,2))
plot(m2)
res2<- residuals(m2, phylo = TRUE) 
shapiro.test(res2)
res2<- res2/sqrt(var(res2))[1] 
rownames(res2)<-rownames(m2$residuals) 
rownames(res2)[(abs(res2)>3)]

#now with % diet needing extensive food handling as the predictor. 
#Poicephalus_cryptoxanthus was IDd as an outlier (phylogenetic residual +/->3) 
#so removed from this model
new_parrot_data <-parrot_data[-c(154),]
new_parrot_comp_data <- comparative.data(phy=parrottree, data= new_parrot_data, 
                                    names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m3<-pgls(log(Max_feed_size) ~ sqrt(Food_handling), new_parrot_comp_data, lambda ='ML')
summary(m3)
par(mfrow=c(2,2))
#model diagnostics 
plot(m3)
res3<- residuals(m3, phylo = TRUE)
shapiro.test(res3)
res3<- res3/sqrt(var(res3))[1] 
rownames(res3)<-rownames(m3$residuals)
rownames(res3)[(abs(res3)>3)]

#now with diet breadth as the predictor
m4<-pgls(log(Max_feed_size) ~ Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m4)
par(mfrow=c(2,2))
plot(m4)
res4<- residuals(m4, phylo = TRUE) 
shapiro.test(res4)
res4<- res4/sqrt(var(res4))[1] 
rownames(res4)<-rownames(m4$residuals) 
rownames(res4)[(abs(res4)>3)]

#now with habitat breadth as the predictor
m5<-pgls(log(Max_feed_size) ~ Habitat_breadth, parrot_comp_data, lambda ='ML')
summary(m5)
par(mfrow=c(2,2))
plot(m5)
res5<- residuals(m5, phylo = TRUE) 
shapiro.test(res5)
res5<- res5/sqrt(var(res5))[1] 
rownames(res5)<-rownames(m5$residuals) 
rownames(res5)[(abs(res5)>3)]

#now with innovation rate as the predictor. First have to filter the data so that only 
#species whose native range is within the geographical regions the innovation rate 
#database included. Note the model includes 'research effort' (n of published papers) 
#as a control (see MS)
parrot_data_inn<-filter(parrot_data, Inn_region==1)
#make a new comparative object using this filtered dataset
parrot_comp_data1 <- comparative.data(phy=parrottree, data= parrot_data_inn, 
                                  names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m6<-pgls(log(Max_feed_size) ~ log(Research_effort) + log(Innovation+1), parrot_comp_data1, 
         lambda ='ML')
summary(m6)
par(mfrow=c(2,2))
plot(m6)
res6<- residuals(m6, phylo = TRUE) 
shapiro.test(res6)
res6<- res6/sqrt(var(res6))[1] 
rownames(res6)<-rownames(m6$residuals) 
rownames(res6)[(abs(res6)>3)]

#now with brain volume as the predictor. Includes body mass as a control (for allometery; 
#see MS)

m7<-pgls(log(Max_feed_size) ~ log(Brain_vol) + Body_mass, parrot_comp_data, lambda ='ML')
summary(m7)
par(mfrow=c(2,2))
plot(m7)
res7<- residuals(m7, phylo = TRUE) 
shapiro.test(res7)
res7<- res7/sqrt(var(res7))[1] 
rownames(res7)<-rownames(m7$residuals) 
rownames(res7)[(abs(res7)>3)]

#now with IUCN as the predictor
m8<-pgls(log(Max_feed_size) ~ log(IUCN_code), parrot_comp_data, lambda ='ML')
summary(m8)
par(mfrow=c(2,2))
plot(m8)
res8<- residuals(m8, phylo = TRUE) 
shapiro.test(res8)
res8<- res8/sqrt(var(res8))[1] 
rownames(res8)<-rownames(m8$residuals) 
rownames(res8)[(abs(res8)>3)]

#Second set: running phylogenetic logistic regressions with communal roosting as the 
#outcome, correlated against every other wild biology predictor in turn. Being categorical 
#with two levels (yes=1, no=0) we used phylogenetic logistic regression to analyse 
#communal roosting as an outcome. Starting with maximum feeding group size as the predictor

#have to pare dataset down to complete cases first
parrot_data9<-parrot_data[!is.na(parrot_data$Communal_roost_code),]
parrot_data9<-parrot_data9[!is.na(parrot_data9$Max_feed_size),]
#match up the tip labels to the row names in the dataframe (so the model knows to match the spp names to the tree )
row.names(parrot_data9)<-parrot_data9$Species 
#see if any names missing in dataset but are on tree, create object listing these
miss_species<-name.check(parrottree, parrot_data9)
#pull these names out, and put into new object
species<-miss_species$tree_not_data
#prune tree to remove these species
pruned_parrottree<-drop.tip(phy = parrottree, species)
#check all match up now (should say OK)
name.check(pruned_parrottree, parrot_data9)

# Cecile Ane (one of package's authors) recommended changing logistic_IG10 in the method to logistic_MPLE (done here)
#she does not recommend changing the btol.
m9<-phyloglm(Communal_roost_code ~ log(Max_feed_size), parrot_data9, pruned_parrottree, 
      method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, 
      start.alpha=NULL, boot = 0, full.matrix = TRUE)
summary(m9)
#get n species for reporting
m9$n
plot(m9)
#extracting residuals to check for outliers (i.e., phylogenetic residual +/- >3)
res9<-residuals(m9, phylo=TRUE)
res9<- res9/sqrt(var(res9))[1] #standardises residuals by sqrt of their variance
#take a look to see how they're distributed
hist(res9, breaks=10)
rownames(res9)<-rownames(m9$residuals)#matches the residuals up with the species names 
res9<-as.data.frame(res9)
rownames(res9)[(abs(res9)>3)]#gives the names of the outliers. 

#same again but for food search as the predictor 
parrot_data10<-parrot_data[!is.na(parrot_data$Communal_roost_code),]
parrot_data10<-parrot_data10[!is.na(parrot_data10$Food_search),]
row.names(parrot_data10)<-parrot_data10$Species 
miss_species<-name.check(parrottree, parrot_data10)
species<-miss_species$tree_not_data
pruned_parrottree<-drop.tip(phy = parrottree, species)
name.check(pruned_parrottree, parrot_data10)

m10<-phyloglm(Communal_roost_code ~ Food_search, parrot_data10, pruned_parrottree, 
        method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, 
        start.alpha=NULL, boot = 0, full.matrix = TRUE)
summary(m10)
m10$n
plot(m10)
res10<-residuals(m10, phylo=TRUE)
res10<-res10/sqrt(var(res10))[1] 
hist(res10, breaks=10)
rownames(res10)<-rownames(m10$residuals) 
res10<-as.data.frame(res10)
rownames(res10)[(abs(res10)>3)]. 

#same again but for food handling as the predictor 
parrot_data11<-parrot_data[!is.na(parrot_data$Communal_roost_code),]
parrot_data11<-parrot_data11[!is.na(parrot_data11$Food_handling),]
row.names(parrot_data11)<-parrot_data11$Species 
miss_species<-name.check(parrottree, parrot_data11)
species<-miss_species$tree_not_data
pruned_parrottree<-drop.tip(phy = parrottree, species)
name.check(pruned_parrottree, parrot_data11)

m11<-phyloglm(Communal_roost_code ~ sqrt(Food_handling), parrot_data11, pruned_parrottree, 
    method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, 
    start.alpha=NULL, boot = 0, full.matrix = TRUE)
summary(m11)
m11$n
plot(m11)
res11<-residuals(m11, phylo=TRUE)
res11<-res11/sqrt(var(res11))[1] 
hist(res11, breaks=20)
rownames(res11)<-rownames(m11$residuals) 
res11<-as.data.frame(res11)
rownames(res11)[(abs(res11)>3)]. 

#same again but for diet breadth as the predictor 
parrot_data12<-parrot_data[!is.na(parrot_data$Communal_roost_code),]
parrot_data12<-parrot_data12[!is.na(parrot_data12$Diet_breadth),]
row.names(parrot_data12)<-parrot_data12$Species 
miss_species<-name.check(parrottree, parrot_data12)
species<-miss_species$tree_not_data
pruned_parrottree<-drop.tip(phy = parrottree, species)
name.check(pruned_parrottree, parrot_data12)

m12<-phyloglm(Communal_roost_code ~ Diet_breadth, parrot_data12, pruned_parrottree, 
              method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, 
              start.alpha=NULL, boot = 0, full.matrix = TRUE)
summary(m12)
m12$n
plot(m12)
res12<-residuals(m12, phylo=TRUE)
res12<-res12/sqrt(var(res12))[1] 
hist(res12, breaks=30)
rownames(res12)<-rownames(m12$residuals) 
res12<-as.data.frame(res12)
rownames(res12)[(abs(res12)>3)] 

#same again but for habitat breadth as the predictor 
parrot_data13<-parrot_data[!is.na(parrot_data$Communal_roost_code),]
parrot_data13<-parrot_data13[!is.na(parrot_data13$Habitat_breadth),]
row.names(parrot_data13)<-parrot_data13$Species 
miss_species<-name.check(parrottree, parrot_data13)
species<-miss_species$tree_not_data
pruned_parrottree<-drop.tip(phy = parrottree, species)
name.check(pruned_parrottree, parrot_data13)

m13<-phyloglm(Communal_roost_code ~ Habitat_breadth, parrot_data13, pruned_parrottree, 
              method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL,
              start.alpha=NULL, boot = 0, full.matrix = TRUE)
summary(m13)
m13$n
plot(m13)
res13<-residuals(m13, phylo=TRUE)
res13<-res13/sqrt(var(res13))[1] 
hist(res13, breaks=30)
rownames(res13)<-rownames(m13$residuals) 
res13<-as.data.frame(res13)
rownames(res13)[(abs(res13)>3)] 

#same again but for innovation rate as the predictor
parrot_data_inn<-filter(parrot_data, Inn_region==1)
parrot_data14<-parrot_data_inn[!is.na(parrot_data_inn$Communal_roost_code),]
parrot_data14<-parrot_data14[!is.na(parrot_data14$Innovation),]
parrot_data14<-parrot_data14[!is.na(parrot_data14$Research_effort),]
row.names(parrot_data14)<-parrot_data14$Species 
miss_species<-name.check(parrottree, parrot_data14)
species<-miss_species$tree_not_data
pruned_parrottree<-drop.tip(phy = parrottree, species)
name.check(pruned_parrottree, parrot_data14)

m14<-phyloglm(Communal_roost_code ~ log(Research_effort)+Innovation, parrot_data14, 
      pruned_parrottree, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, 
      start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE)
summary(m14)
m14$n
plot(m14)
res14<-residuals(m14, phylo=TRUE)
res14<-res14/sqrt(var(res14))[1] 
hist(res14, breaks=30)
rownames(res14)<-rownames(m14$residuals) 
res14<-as.data.frame(res14)
rownames(res14)[(abs(res14)>3)] 

#same again, but with relative brain volume as the predictor
parrot_data15<-parrot_data[!is.na(parrot_data$Communal_roost_code),]
parrot_data15<-parrot_data15[!is.na(parrot_data15$Brain_vol),]
parrot_data15<-parrot_data15[!is.na(parrot_data15$Body_mass),]
row.names(parrot_data15)<-parrot_data15$Species 
miss_species<-name.check(parrottree, parrot_data15)
species<-miss_species$tree_not_data
pruned_parrottree<-drop.tip(phy = parrottree, species)
name.check(pruned_parrottree, parrot_data15)

m15<-phyloglm(Communal_roost_code ~ log(Brain_vol)+log(Body_mass), parrot_data15, 
        pruned_parrottree, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, 
        start.beta=NULL, start.alpha=NULL, boot = 0, full.matrix = TRUE)
summary(m15)
m15$n
plot(m15)
res15<-residuals(m15)
res15<-res15/sqrt(var(res15))[1] 
hist(res15, breaks=20)
rownames(res15)<-rownames(m15$residuals) 
res15<-as.data.frame(res15)
rownames(res15)[(abs(res15)>3)] 

#same again, but with IUCN status as the predictor
parrot_data16<-parrot_data[!is.na(parrot_data$Communal_roost_code),]
parrot_data16<-parrot_data16[!is.na(parrot_data16$IUCN_code),]
row.names(parrot_data16)<-parrot_data16$Species 
miss_species<-name.check(parrottree, parrot_data16)
species<-miss_species$tree_not_data
pruned_parrottree<-drop.tip(phy = parrottree, species)
name.check(pruned_parrottree, parrot_data16)

m16<-phyloglm(Communal_roost_code ~ IUCN_code, parrot_data16, pruned_parrottree, 
        method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4, start.beta=NULL, 
        start.alpha=NULL, boot = 0, full.matrix = TRUE)
summary(m16)
m16$n
plot(m16)
res16<-residuals(m16)
res16<-res16/sqrt(var(res16))[1] 
hist(res16, breaks=30)
rownames(res16)<-rownames(m16$residuals) 
res16<-as.data.frame(res16)
rownames(res16)[(abs(res16)>3)] 

#Third set: running PGLS regressions with food search as the outcome, correlated against
#every other wild biology predictor in turn.
#starting with maximum feeding group size as the predictor

m17<-pgls(sqrt(Food_search) ~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m17)
par(mfrow=c(2,2))
plot(m17)
res17<-residuals(m17, phylo = TRUE) 
shapiro.test(res17)
res17<- res17/sqrt(var(res17))[1] 
rownames(res17)<-rownames(m17$residuals) 
rownames(res17)[(abs(res17)>3)]
par(mfrow=c(1,1))

#now with communal roosting as the predictor. Note that Brotogeris_versicolurus is IDd as 
#an outlier, but removing it results in several other species being IDd as outliers. In 
#the interest of retaining sample size, we decided to leave Brotogeris_versicolurus in.
m18<-pgls(Food_search ~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m18)
par(mfrow=c(2,2))
plot(m18)
res18<-residuals(m18, phylo = TRUE) 
shapiro.test(res18)
res18<- res18/sqrt(var(res18))[1] 
rownames(res18)<-rownames(m18$residuals) 
rownames(res18)[(abs(res18)>3)]

#now with food handling as the predictor
m19<-pgls(Food_search ~ sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m19)
par(mfrow=c(2,2))
plot(m19)
res19<-residuals(m19, phylo = TRUE) 
shapiro.test(res19)
res19<- res19/sqrt(var(res19))[1] 
rownames(res19)<-rownames(m19$residuals) 
rownames(res19)[(abs(res19)>3)]

#now with diet breadth as the predictor. Poicephalus_cryptoxanthus was IDd as an outlier
#(phylogenetic residual +/->3) so removed from this model
new_parrot_data <-parrot_data[-c(158),]
new_parrot_comp_data <- comparative.data(phy=parrottree, data= new_parrot_data, 
                              names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)

m20<-pgls(Food_search ~ Diet_breadth, new_parrot_comp_data, lambda ='ML')
summary(m20)
par(mfrow=c(2,2))
plot(m20)
res20<-residuals(m20, phylo = TRUE) 
shapiro.test(res20)
res20<- res20/sqrt(var(res20))[1] 
rownames(res20)<-rownames(m20$residuals) 
rownames(res20)[(abs(res20)>3)]

#now with habitat breadth as the predictor
m21<-pgls(Food_search ~ Habitat_breadth, parrot_comp_data, lambda ='ML')
summary(m21)
par(mfrow=c(2,2))
plot(m21)
res21<-residuals(m21, phylo = TRUE) 
shapiro.test(res21)
res21<- res21/sqrt(var(res21))[1] 
rownames(res21)<-rownames(m21$residuals) 
rownames(res21)[(abs(res21)>3)]

#now with innovation rate as the predictor
m22<-pgls(sqrt(Food_search) ~ log(Research_effort) + (Innovation), parrot_comp_data1, lambda ='ML')
summary(m22)
par(mfrow=c(2,2))
plot(m22)
res22<-residuals(m22, phylo = TRUE) 
shapiro.test(res22)
res22<-res22/sqrt(var(res22))[1] 
rownames(res22)<-rownames(m22$residuals) 
rownames(res22)[(abs(res22)>3)]

#now with brain volume as the predictor
m23<-pgls(Food_search ~ log(Body_mass) +log(Brain_vol), parrot_comp_data, lambda ='ML')
summary(m23)
par(mfrow=c(2,2))
plot(m23)
res23<-residuals(m23, phylo = TRUE) 
shapiro.test(res23)
res23<-res23/sqrt(var(res23))[1] 
rownames(res23)<-rownames(m23$residuals) 
rownames(res23)[(abs(res23)>3)]

#now with IUCN status as the predictor
m24<-pgls(Food_search ~ IUCN_code, parrot_comp_data, lambda ='ML')
summary(m24)
par(mfrow=c(2,2))
plot(m24)
res24<-residuals(m24, phylo = TRUE) 
shapiro.test(res24)
res24<-res24/sqrt(var(res24))[1] 
rownames(res24)<-rownames(m24$residuals) 
rownames(res24)[(abs(res24)>3)]

#Fourth set: running PGLS regressions with food handling as the outcome, correlated against
#every other wild biology predictor in turn.
#starting with maximum feeding group size as the predictor

m25<-pgls(sqrt(Food_handling) ~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m25)
par(mfrow=c(2,2))
plot(m25)
res25<-residuals(m25, phylo = TRUE) 
shapiro.test(res25)
res25<-res25/sqrt(var(res25))[1] 
rownames(res25)<-rownames(m25$residuals) 
rownames(res25)[(abs(res25)>3)]

#now with communal roosting as the predictor
m26<-pgls(Food_handling ~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m26)
par(mfrow=c(2,2))
plot(m26)
res26<-residuals(m26, phylo = TRUE) 
shapiro.test(res26)
res26<-res26/sqrt(var(res26))[1] 
rownames(res26)<-rownames(m26$residuals) 
rownames(res26)[(abs(res26)>3)]

#now with food search as the predictor. Alisterus amboinensis was an outlier with a 
#residual > -3. Removing it, though, resulted in the subsequent model's residuals then 
#becoming non-normally distributed with other species having residual values identifying 
#them as outliers. In the interest of sample size, we therefore decided to retain this 
#species 
m27<-pgls(sqrt(Food_handling) ~ Food_search, parrot_comp_data, lambda ='ML')
summary(m27)
par(mfrow=c(2,2))
plot(m27)
res27<-residuals(m27, phylo = TRUE) 
shapiro.test(res27)
res27<-res27/sqrt(var(res27))[1] 
rownames(res27)<-rownames(m27$residuals) 
rownames(res27)[(abs(res27)>3)]

#now with diet breadth as the predictor
m28<-pgls(sqrt(Food_handling) ~ Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m28)
par(mfrow=c(2,2))
plot(m28)
res28<-residuals(m28, phylo = TRUE) 
shapiro.test(res28)
res28<-res28/sqrt(var(res28))[1] 
rownames(res28)<-rownames(m28$residuals) 
rownames(res28)[(abs(res28)>3)]

#now with habitat breadth as the predictor
m29<-pgls(sqrt(Food_handling) ~ sqrt(Habitat_breadth), parrot_comp_data, lambda ='ML')
summary(m29)
par(mfrow=c(2,2))
plot(m29)
res29<-residuals(m29, phylo = TRUE) 
shapiro.test(res29)
res29<-res29/sqrt(var(res29))[1] 
rownames(res29)<-rownames(m29$residuals) 
rownames(res29)[(abs(res29)>3)]

#now with innovation rate as the predictor
m29<-pgls(log(Food_handling+1) ~ log(Research_effort) +Innovation, parrot_comp_data1, 
          lambda ='ML')
summary(m29)
par(mfrow=c(2,2))
plot(m29)
res29<-residuals(m29, phylo = TRUE) 
shapiro.test(res29)
res29<-res29/sqrt(var(res29))[1] 
rownames(res29)<-rownames(m29$residuals) 
rownames(res29)[(abs(res29)>3)]

#now with brain volume as a predictor
m30<-pgls(Food_handling ~ log(Brain_vol)+log(Body_mass), parrot_comp_data,lambda ='ML')
summary(m30)
par(mfrow=c(2,2))
plot(m30)
res30<-residuals(m30, phylo = TRUE) 
shapiro.test(res30)
res30<-res30/sqrt(var(res30))[1] 
rownames(res30)<-rownames(m30$residuals) 
rownames(res30)[(abs(res30)>3)]

#now with IUCN as a predictor
m31<-pgls(Food_handling ~ IUCN_code, parrot_comp_data,lambda ='ML')
summary(m31)
par(mfrow=c(2,2))
plot(m31)
res31<-residuals(m31, phylo = TRUE) 
shapiro.test(res31)
res31<-res31/sqrt(var(res31))[1] 
rownames(res31)<-rownames(m31$residuals) 
rownames(res31)[(abs(res31)>3)]

#Fifth set: running PGLS regressions with diet breadth as the outcome, correlated against
#every other wild biology predictor in turn.
#starting with maximum feeding group size as the predictor
m32<-pgls(log(Diet_breadth)~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m32)
par(mfrow=c(2,2))
plot(m32)
res32<-residuals(m32, phylo = TRUE) 
shapiro.test(res32)
res32<-res32/sqrt(var(res32))[1] 
rownames(res32)<-rownames(m32$residuals) 
rownames(res32)[(abs(res32)>3)]

#now with communal roosting as the predictor
m33<-pgls(sqrt(Diet_breadth)~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m33)
par(mfrow=c(2,2))
plot(m33)
res33<-residuals(m33, phylo = TRUE) 
shapiro.test(res33)
res33<-res33/sqrt(var(res33))[1] 
rownames(res33)<-rownames(m33$residuals) 
rownames(res33)[(abs(res33)>3)]

#now with food search as the predictor
m34<-pgls(Diet_breadth~ Food_search, parrot_comp_data, lambda ='ML')
summary(m34)
par(mfrow=c(2,2))
plot(m34)
res34<-residuals(m34, phylo = TRUE) 
shapiro.test(res34)
res34<-res34/sqrt(var(res34))[1] 
rownames(res34)<-rownames(m34$residuals) 
rownames(res34)[(abs(res34)>3)]

#now with food handling as the predictor
m35<-pgls(Diet_breadth~ sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m35)
par(mfrow=c(2,2))
plot(m35)
res35<-residuals(m35, phylo = TRUE) 
shapiro.test(res35)
res35<-res35/sqrt(var(res35))[1] 
rownames(res35)<-rownames(m35$residuals) 
rownames(res35)[(abs(res35)>3)]

#now with habitat breadth as the predictor
m36<-pgls(Diet_breadth~ Habitat_breadth, parrot_comp_data, lambda ='ML')
summary(m36)
par(mfrow=c(2,2))
plot(m36)
res36<-residuals(m36, phylo = TRUE) 
shapiro.test(res36)
res36<-res36/sqrt(var(res36))[1] 
rownames(res36)<-rownames(m36$residuals) 
rownames(res36)[(abs(res36)>3)]

#now with innovation rate as the predictor
m37<-pgls(Diet_breadth~ log(Innovation+1)+log(Research_effort), parrot_comp_data1, 
          lambda ='ML')
summary(m37)
par(mfrow=c(2,2))
plot(m37)
res37<-residuals(m37, phylo = TRUE) 
shapiro.test(res37)
res37<-res37/sqrt(var(res37))[1] 
rownames(res37)<-rownames(m37$residuals) 
rownames(res37)[(abs(res37)>3)]

#now with brain volume as the predictor
m38<-pgls(Diet_breadth~ log(Brain_vol)+log(Body_mass), parrot_comp_data, lambda ='ML')
summary(m38)
par(mfrow=c(2,2))
plot(m38)
res38<-residuals(m38, phylo = TRUE) 
shapiro.test(res38)
res38<-res38/sqrt(var(res38))[1] 
rownames(res38)<-rownames(m38$residuals) 
rownames(res38)[(abs(res38)>3)]

#now with IUCN as the predictor
m39<-pgls(log(Diet_breadth)~ log(IUCN_code), parrot_comp_data, lambda ='ML')
summary(m39)
par(mfrow=c(2,2))
plot(m39)
res39<-residuals(m39, phylo = TRUE) 
shapiro.test(res39)
res39<-res39/sqrt(var(res39))[1] 
rownames(res39)<-rownames(m39$residuals) 
rownames(res39)[(abs(res39)>3)]

#Sixth set: running PGLS regressions with habitat breadth as the outcome, correlated against
#every other wild biology predictor in turn.
#starting with maximum feeding group size as the predictor
m40<-pgls(Habitat_breadth~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m40)
par(mfrow=c(2,2))
plot(m40)
res40<-residuals(m40, phylo = TRUE) 
shapiro.test(res40)
res40<-res40/sqrt(var(res40))[1] 
rownames(res40)<-rownames(m40$residuals) 
rownames(res40)[(abs(res40)>3)]

#now with communal roosting as a predictor
m41<-pgls(Habitat_breadth~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m41)
par(mfrow=c(2,2))
plot(m41)
res41<-residuals(m41, phylo = TRUE) 
shapiro.test(res41)
res41<-res41/sqrt(var(res41))[1] 
rownames(res41)<-rownames(m41$residuals) 
rownames(res41)[(abs(res41)>3)]

#now with food search as a predictor
m42<-pgls(Habitat_breadth~ Food_search, parrot_comp_data, lambda ='ML')
summary(m42)
par(mfrow=c(2,2))
plot(m42)
res42<-residuals(m42, phylo = TRUE) 
shapiro.test(res42)
res42<-res42/sqrt(var(res42))[1] 
rownames(res42)<-rownames(m42$residuals) 
rownames(res42)[(abs(res42)>3)]

#now with food handling as a predictor
m43<-pgls(Habitat_breadth~ sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m43)
par(mfrow=c(2,2))
plot(m43)
res43<-residuals(m43, phylo = TRUE) 
shapiro.test(res43)
res43<-res43/sqrt(var(res43))[1] 
rownames(res43)<-rownames(m43$residuals) 
rownames(res43)[(abs(res43)>3)]

#now with diet breadth as a predictor
m44<-pgls(Habitat_breadth~ Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m44)
par(mfrow=c(2,2))
plot(m44)
res44<-residuals(m44, phylo = TRUE) 
shapiro.test(res44)
res44<-res44/sqrt(var(res44))[1] 
rownames(res44)<-rownames(m44$residuals) 
rownames(res44)[(abs(res44)>3)]

#now with innovation rate as a predictor
m45<-pgls(Habitat_breadth~ Innovation+log(Research_effort), parrot_comp_data1, 
          lambda ='ML')
summary(m45)
par(mfrow=c(2,2))
plot(m45)
res45<-residuals(m45, phylo = TRUE) 
shapiro.test(res45)
res45<-res45/sqrt(var(res45))[1] 
rownames(res45)<-rownames(m45$residuals) 
rownames(res45)[(abs(res45)>3)]

#now with brain volume as a predictor
m46<-pgls(Habitat_breadth~ log(Brain_vol)+log(Body_mass), parrot_comp_data, lambda ='ML')
summary(m46)
par(mfrow=c(2,2))
plot(m46)
res46<-residuals(m46, phylo = TRUE) 
shapiro.test(res46)
res46<-res46/sqrt(var(res46))[1] 
rownames(res46)<-rownames(m46$residuals) 
rownames(res46)[(abs(res46)>3)]

#now with IUCN as a predictor
m47<-pgls(Habitat_breadth~ IUCN_code, parrot_comp_data, lambda ='ML')
summary(m47)
par(mfrow=c(2,2))
plot(m47)
res47<-residuals(m47, phylo = TRUE) 
shapiro.test(res47)
res47<-res47/sqrt(var(res47))[1] 
rownames(res47)<-rownames(m47$residuals) 
rownames(res47)[(abs(res47)>3)]

#Seventh set: running PGLS regressions with (relative) brain volume as the outcome, 
#correlated against every other wild biology predictor in turn.
#starting with maximum feeding group size as the predictor
m48<-pgls(log(Brain_vol)~ log(Body_mass)+log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m48)
par(mfrow=c(2,2))
plot(m48)
res48<-residuals(m48, phylo = TRUE) 
shapiro.test(res48)
res48<-res48/sqrt(var(res48))[1] 
rownames(res48)<-rownames(m48$residuals) 
rownames(res48)[(abs(res48)>3)]

#now with communal roosting as the predictor
m49<-pgls(log(Brain_vol)~ log(Body_mass)+Communal_roost, parrot_comp_data, lambda ='ML')
summary(m49)
par(mfrow=c(2,2))
plot(m49)
res49<-residuals(m49, phylo = TRUE) 
shapiro.test(res49)
res49<-res49/sqrt(var(res49))[1] 
rownames(res49)<-rownames(m49$residuals) 
rownames(res49)[(abs(res49)>3)]

#now with food search as the predictor
m50<-pgls(log(Brain_vol)~ log(Body_mass)+Food_search, parrot_comp_data, lambda ='ML')
summary(m50)
par(mfrow=c(2,2))
plot(m50)
res50<-residuals(m50, phylo = TRUE) 
shapiro.test(res50)
res50<-res50/sqrt(var(res50))[1] 
rownames(res50)<-rownames(m50$residuals) 
rownames(res50)[(abs(res50)>3)]

#now with food handling as the predictor
m51<-pgls(log(Brain_vol)~ log(Body_mass)+sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m51)
par(mfrow=c(2,2))
plot(m51)
res51<-residuals(m51, phylo = TRUE) 
shapiro.test(res51)
res51<-res51/sqrt(var(res51))[1] 
rownames(res51)<-rownames(m51$residuals) 
rownames(res51)[(abs(res51)>3)]

#now with diet breadth as the predictor
m52<-pgls(log(Brain_vol)~ log(Body_mass)+Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m52)
par(mfrow=c(2,2))
plot(m52)
res52<-residuals(m52, phylo = TRUE) 
shapiro.test(res52)
res52<-res52/sqrt(var(res52))[1] 
rownames(res52)<-rownames(m52$residuals) 
rownames(res52)[(abs(res52)>3)]

#now with habitat breadth as the predictor
m53<-pgls(log(Brain_vol)~ log(Body_mass)+Habitat_breadth, parrot_comp_data, lambda ='ML')
summary(m53)
par(mfrow=c(2,2))
plot(m53)
res53<-residuals(m53, phylo = TRUE) 
shapiro.test(res53)
res53<-res53/sqrt(var(res53))[1] 
rownames(res53)<-rownames(m53$residuals) 
rownames(res53)[(abs(res53)>3)]

#now with innovation rate as the predictor
m54<-pgls(log(Brain_vol)~ log(Body_mass)+Innovation+Research_effort, parrot_comp_data1, 
          lambda ='ML')
summary(m54)
par(mfrow=c(2,2))
plot(m54)
res54<-residuals(m54, phylo = TRUE) 
shapiro.test(res54)
res54<-res54/sqrt(var(res54))[1] 
rownames(res54)<-rownames(m54$residuals) 
rownames(res54)[(abs(res54)>3)]

#now IUCN status as the predictor
m55<-pgls(log(Brain_vol)~ log(Body_mass)+IUCN_code, parrot_comp_data, lambda ='ML')
summary(m55)
par(mfrow=c(2,2))
plot(m55)
res55<-residuals(m55, phylo = TRUE) 
shapiro.test(res55)
res55<-res55/sqrt(var(res55))[1] 
rownames(res55)<-rownames(m55$residuals) 
rownames(res55)[(abs(res55)>3)]

#Eighth set: running PGLS regressions with IUCN status as the outcome, 
#correlated against every other wild biology predictor in turn.
#starting with maximum feeding group size as the predictor
m56<-pgls(IUCN_code~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(m56)
par(mfrow=c(2,2))
plot(m56)
res56<-residuals(m56, phylo = TRUE) 
shapiro.test(res56)
res56<-res56/sqrt(var(res56))[1] 
rownames(res56)<-rownames(m56$residuals) 
rownames(res56)[(abs(res56)>3)]

#now with communal roosting as the predictor
m57<-pgls(sqrt(IUCN_code)~ Communal_roost, parrot_comp_data, lambda ='ML')
summary(m57)
par(mfrow=c(2,2))
plot(m57)
res57<-residuals(m57, phylo = TRUE) 
shapiro.test(res57)
res57<-res57/sqrt(var(res57))[1] 
rownames(res57)<-rownames(m57$residuals) 
rownames(res57)[(abs(res57)>3)]

#now with food search as the predictor
m58<-pgls(IUCN_code~ Food_search, parrot_comp_data, lambda ='ML')
summary(m58)
par(mfrow=c(2,2))
plot(m58)
res58<-residuals(m58, phylo = TRUE) 
shapiro.test(res58)
res58<-res58/sqrt(var(res58))[1] 
rownames(res58)<-rownames(m58$residuals) 
rownames(res58)[(abs(res58)>3)]

#now with food handling as the predictor
m59<-pgls(IUCN_code~ sqrt(Food_handling), parrot_comp_data, lambda ='ML')
summary(m59)
par(mfrow=c(2,2))
plot(m59)
res59<-residuals(m59, phylo = TRUE) 
shapiro.test(res59)
res59<-res59/sqrt(var(res59))[1] 
rownames(res59)<-rownames(m59$residuals) 
rownames(res59)[(abs(res59)>3)]

#now with diet breadth as the predictor
m60<-pgls(IUCN_code~ Diet_breadth, parrot_comp_data, lambda ='ML')
summary(m60)
par(mfrow=c(2,2))
plot(m60)
res60<-residuals(m60, phylo = TRUE) 
shapiro.test(res60)
res60<-res60/sqrt(var(res60))[1] 
rownames(res60)<-rownames(m60$residuals) 
rownames(res60)[(abs(res60)>3)]

#now with habitat breadth as the predictor
m61<-pgls(IUCN_code~ Habitat_breadth, parrot_comp_data, lambda ='ML')
summary(m61)
par(mfrow=c(2,2))
plot(m61)
res61<-residuals(m61, phylo = TRUE) 
shapiro.test(res61)
res61<-res61/sqrt(var(res61))[1] 
rownames(res61)<-rownames(m61$residuals) 
rownames(res61)[(abs(res61)>3)]

#now with innovation rate as the predictor
m62<-pgls(IUCN_code~ Innovation+log(Research_effort), parrot_comp_data1, lambda ='ML')
summary(m62)
par(mfrow=c(2,2))
plot(m62)
res62<-residuals(m62, phylo = TRUE) 
shapiro.test(res62)
res62<-res62/sqrt(var(res62))[1] 
rownames(res62)<-rownames(m62$residuals) 
rownames(res62)[(abs(res62)>3)]

#now with brain volume as the predictor
m63<-pgls(log(IUCN_code)~ log(Brain_vol)+log(Body_mass), parrot_comp_data, lambda ='ML')
summary(m63)
par(mfrow=c(2,2))
plot(m63)
res63<-residuals(m63, phylo = TRUE) 
shapiro.test(res63)
res63<-res63/sqrt(var(res63))[1] 
rownames(res63)<-rownames(m63$residuals) 
rownames(res63)[(abs(res63)>3)]

