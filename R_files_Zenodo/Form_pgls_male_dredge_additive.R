#PGLS plumage versus vocal traits in males, ecological and behavior additive effects


Form.male<-read.csv("Form_total_male_ok.csv",row.names = 1)
head(Form.male)

#Contrast values were inverted so that small values now represent lower contrast and higher  values represent higher contrast. (See Analyses of plumage data section for more info)
Form.male$contrast.dorsal.inv<-Form.male$contrast.dorsal*-1 #inversion of dorsal contrast
Form.male$contrast.ventral.inv<-Form.male$contrast.ventral*-1 #inversion of ventral contrast
Form.male$contrast.wing.coverts.inv<-Form.male$contrast.wing.coverts*-1 #inversion of wing contrast
head(Form.male)


library(tidyr)
#remove NAs
Form.male<-drop_na(Form.male)
Form.male

#Load packages
library(ape)
library(geiger)
library(nlme)
library(phytools)

library(ggplot2)
library(dplyr)

library(MuMIn)

#Open tree 
Form.tree <- read.tree("formicivorini_UCEs_SAAC.tre")
plot(Form.tree)
Form.tree

#Check tree and data
obj<-name.check(Form.tree,Form.male)
obj

#Remove tips that we don't have data
Form.tree<- drop.tip(Form.tree, obj$tree_not_data)
plot(Form.tree)

#Check tree and data again
name.check(Form.tree, Form.male)

#Branch length transformation. Necessary to avoid error you were getting before. Does not affect results. 
Form.tree$edge.length <- Form.tree$edge.length * 1000

####################################################################################################################################################################################################################################################################################
#Models 

#1. maxPower.dorsal	(BM)
#versus Note.count (BM-->Lambda) + Note.type (BM) + Delta.Time..s. (OU) + Peak.Freq..Hz. (BM) +	Note.diversity (BM-->Lambda em 2°) + Note.rate	(BM) + Song.bandwidth (OU) + Loudsong.mod.rate (OU-->Lambda em 2°)

#BM model
mod.bro1<-gls(maxPower.dorsal ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro1<-dredge(global.model=mod.bro1) 
head(d.bro1)

write.table(d.bro1,"AICs_MPdorsalXvocaleco_BM_male.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro1<-coefTable(d.bro1) 
se.bro1

se.bro1.df<-do.call(rbind.data.frame,se.bro1 )
dim(se.bro1.df)
str(se.bro1.df)

write.table(se.bro1.df[1:7677,],"SE_MPdorsalXvocaleco_BM_male1.csv",sep = ",")
write.table(se.bro1.df[7678:15360,],"SE_MPdorsalXvocaleco_BM_male2.csv",sep = ",")



######
#OU model
mod.ou1<-gls(maxPower.dorsal ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou1<-dredge(global.model=mod.ou1)

write.table(d.ou1,"AICs_MPdorsalXvocaleco_OU_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.ou1<-coefTable(d.ou1) 
head(se.ou1)
se.ou1.df<-do.call(rbind.data.frame,se.ou1 )
dim(se.ou1.df)
str(se.ou1.df)

write.table(se.ou1.df[1:7682,],"SE_MPdorsalXvocaleco_OU_male1.csv",sep = ",")
write.table(se.ou1.df[7683:15360,],"SE_MPdorsalXvocaleco_OU_male2.csv",sep = ",")


######
#Pagel's Lambda   
mod.pl1<-gls(maxPower.dorsal ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl1<-dredge(global.model=mod.pl1)

write.table(d.pl1,"AICs_MPdorsalXvocaleco_PL_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl1<-coefTable(d.pl1) 
head(se.pl1)
se.pl1.df<-do.call(rbind.data.frame,se.pl1 )
dim(se.pl1.df)
str(se.pl1.df)

write.table(se.pl1.df[1:7682,],"SE_MPdorsalXvocaleco_pl_male1.csv",sep = ",")
write.table(se.pl1.df[7683:15360,],"SE_MPdorsalXvocaleco_pl_male2.csv",sep = ",")

#####################################
##################################
#Summary statistics of the best models 

#Dorsal maximum energy X loudsong duration >> all BM
#1a
mbest.bm1a<-gls(maxPower.dorsal ~ Delta.Time..s. + Habitat_exposure + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.bm1a)

#1b
mbest.bm1b<-gls(maxPower.dorsal ~ Delta.Time..s. + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.bm1b)


#Dorsal maximum energy X Loudsong.mod.rate >> BM e PL
#11a BM
mbest.bm11a<-gls(maxPower.dorsal ~ Loudsong.mod.rate + Habitat_exposure + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.bm11a)

#11b BM
mbest.bm11b<-gls(maxPower.dorsal ~ Loudsong.mod.rate + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.bm11b)

#11c PL
mbest.pl11c<-gls(maxPower.dorsal ~ Loudsong.mod.rate + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.pl11c)

#Dorsal maximum energy X note count >> BM 
#111a
mbest.bm111a<-gls(maxPower.dorsal ~ Note.count + Habitat_exposure + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.bm111a)

#111b
mbest.bm111b<-gls(maxPower.dorsal ~ Note.count + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.bm111b)

#Dorsal maximum energy X note diversity >> BM 
#1111a
mbest.bm1111a<-gls(maxPower.dorsal ~ Note.diversity + Habitat_exposure + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.bm1111a)

#1111b
mbest.bm1111b<-gls(maxPower.dorsal ~ Note.diversity + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.bm1111b)


###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#2.	LuminanceMean.dorsal  (OU)	
#versus Note.count (BM) + Note.type (BM) + Delta.Time..s. (OU) + Peak.Freq..Hz. (BM) +	Note.diversity (BM) + Note.rate	(BM) + Song.bandwidth (OU) + Loudsong.mod.rate (OU)

#BM model
mod.bro2<-gls(LuminanceMean.dorsal ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro2<-dredge(global.model=mod.bro2)

write.table(d.bro2,"AICs_LumdorsalXvocaleco_BM_male.csv",sep = ",",row.names = F)


#obtain SE (standard error)
se.bro2<-coefTable(d.bro2) 
se.bro2.df<-do.call(rbind.data.frame,se.bro2 )
dim(se.bro2.df)
str(se.bro2.df)

write.table(se.bro2.df[1:7675,],"SE_LumdorsalXvocaleco_BM_male1.csv",sep = ",")
write.table(se.bro2.df[7676:15360,],"SE_LumdorsalXvocaleco_BM_male2.csv",sep = ",")

#OU model
mod.ou2<-gls(LuminanceMean.dorsal ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou2<-dredge(global.model=mod.ou2)

write.table(d.ou2,"AICs_LumdorsalXvocaleco_OU_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.ou2<-coefTable(d.ou2) 
se.ou2.df<-do.call(rbind.data.frame,se.ou2 )
dim(se.ou2.df)
str(se.ou2.df)

write.table(se.ou2.df[1:7677,],"SE_LumdorsalXvocaleco_OU_male1.csv",sep = ",")
write.table(se.ou2.df[7678:15360,],"SE_LumdorsalXvocaleco_OU_male2.csv",sep = ",")


######
#Pagel's Lambda model
mod.pl2<-gls(LuminanceMean.dorsal ~ Note.count + Note.diversity + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.7, phy = Form.tree), data = Form.male, method = "ML")
mod.pl2

d.pl2<-dredge(global.model=mod.pl2)

write.table(d.pl2,"AICs_LumdorsalXvocaleco_pl_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl2<-coefTable(d.pl2) 
se.pl2.df<-do.call(rbind.data.frame,se.pl2 )
dim(se.pl2.df)
str(se.pl2.df)

write.table(se.pl2.df[1:7677,],"SE_LumdorsalXvocaleco_pl_male1.csv",sep = ",")
write.table(se.pl2.df[7678:15360,],"SE_LumdorsalXvocaleco_pl_male2.csv",sep = ",")

#####################################
##################################
#Summary statistics of the best models 

#Dorsal luminance X loudsong bandwidth >> all OU
mbest.ou2a<-gls(LuminanceMean.dorsal ~ Song.bandwidth + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

summary(mbest.ou2a)

mbest.ou2b<-gls(LuminanceMean.dorsal ~ Song.bandwidth + Habitat_exposure + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

summary(mbest.ou2b)

mbest.ou2c<-gls(LuminanceMean.dorsal ~ Song.bandwidth + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

summary(mbest.ou2c)

mbest.ou2d<-gls(LuminanceMean.dorsal ~ Song.bandwidth + Foraging_strata, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

summary(mbest.ou2d)

mbest.ou2e<-gls(LuminanceMean.dorsal ~ Song.bandwidth + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

summary(mbest.ou2e)


#Dorsal luminance X Loudsong.mod.rate >> all lambda
mbest.pl2a<-gls(LuminanceMean.dorsal ~ Loudsong.mod.rate + Mixed.Flocking, correlation  = corPagel(0.7, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.pl2a)

mbest.pl2b<-gls(LuminanceMean.dorsal ~ Loudsong.mod.rate + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(0.7, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.pl2b)

#Dorsal luminance X note type >> all OU
mbest.ou22a<-gls(LuminanceMean.dorsal ~ Note.type + Foraging_strata, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

summary(mbest.ou22a)

mbest.ou22b<-gls(LuminanceMean.dorsal ~ Note.type + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

summary(mbest.ou22b)

###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#3. contrast.dorsal	(BM)

#BM model
mod.bro3<-gls(contrast.dorsal.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro3<-dredge(global.model=mod.bro3)

write.table(d.bro3,"AICs_ContdorsalXvocaleco_BM_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro3<-coefTable(d.bro3) 
se.bro3.df<-do.call(rbind.data.frame,se.bro3 )
dim(se.bro3.df)
str(se.bro3.df)

write.table(se.bro3.df[1:7671,],"SE_ContdorsalXvocaleco_BM_male1.csv",sep = ",")
write.table(se.bro3.df[7672:15360,],"SE_ContdorsalXvocaleco_BM_male2.csv",sep = ",")


#OU model
mod.ou3<-gls(contrast.dorsal.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou3<-dredge(global.model=mod.ou3)

write.table(d.ou3,"AICs_ContdorsalXvocaleco_OU_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.ou3<-coefTable(d.ou3) 
se.ou3.df<-do.call(rbind.data.frame,se.ou3)
dim(se.ou3.df)
str(se.ou3.df)

write.table(se.ou3.df[1:7667,],"SE_ContdorsalXvocaleco_OU_male1.csv",sep = ",")
write.table(se.ou3.df[7668:15360,],"SE_ContdorsalXvocaleco_OU_male2.csv",sep = ",")

######
#Pagel's Lambda model
mod.pl3<-gls(contrast.dorsal.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl3<-dredge(global.model=mod.pl3)

write.table(d.pl3,"AICs_ContdorsalXvocaleco_pl_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl3<-coefTable(d.pl3) 
se.pl3.df<-do.call(rbind.data.frame,se.pl3)
dim(se.pl3.df)
str(se.pl3.df)

write.table(se.pl3.df[1:7667,],"SE_ContdorsalXvocaleco_pl_male1.csv",sep = ",")
write.table(se.pl3.df[7668:15360,],"SE_ContdorsalXvocaleco_pl_male2.csv",sep = ",")

###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#4. maxPower.ventral 		

#BM model
mod.bro4<-gls(maxPower.ventral ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro4<-dredge(global.model=mod.bro4)

write.table(d.bro4,"AICs_MPventralXvocaleco_BM_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro4<-coefTable(d.bro4) 
se.bro4.df<-do.call(rbind.data.frame,se.bro4 )
dim(se.bro4.df)
str(se.bro4.df)

write.table(se.bro4.df[1:7675,],"SE_MPventralXvocaleco_BM_male1.csv",sep = ",")
write.table(se.bro4.df[7676:15360,],"SE_MPventralXvocaleco_BM_male2.csv",sep = ",")


#OU model
mod.ou4<-gls(maxPower.ventral ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou4<-dredge(global.model=mod.ou4)

write.table(d.ou4,"AICs_MPventralXvocaleco_OU_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.ou4<-coefTable(d.ou4) 
se.ou4.df<-do.call(rbind.data.frame,se.ou4 )
dim(se.ou4.df)
str(se.ou4.df)

write.table(se.ou4.df[1:7669,],"SE_MPventralXvocaleco_OU_male1.csv",sep = ",")
write.table(se.ou4.df[7670:15360,],"SE_MPventralXvocaleco_OU_male2.csv",sep = ",")

######
#Pagel's Lambda model
mod.pl4<-gls(maxPower.ventral ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl4<-dredge(global.model=mod.pl4)

write.table(d.pl4,"AICs_MPventralXvocaleco_pl_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl4<-coefTable(d.pl4) 
se.pl4.df<-do.call(rbind.data.frame,se.pl4 )
dim(se.pl4.df)
str(se.pl4.df)

write.table(se.pl4.df[1:7669,],"SE_MPventralXvocaleco_pl_male1.csv",sep = ",")
write.table(se.pl4.df[7670:15360,],"SE_MPventralXvocaleco_pl_male2.csv",sep = ",")


###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#5. LuminanceMean.ventral

#BM model
mod.bro5<-gls(LuminanceMean.ventral ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro5<-dredge(global.model=mod.bro5)

write.table(d.bro5,"AICs_lumventralXvocaleco_BM_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro5<-coefTable(d.bro5) 
se.bro5.df<-do.call(rbind.data.frame,se.bro5 )
dim(se.bro5.df)
str(se.bro5.df)

write.table(se.bro5.df[1:7669,],"SE_lumventralXvocaleco_BM_male1.csv",sep = ",")
write.table(se.bro5.df[7670:15360,],"SE_lumventralXvocaleco_BM_male2.csv",sep = ",")


#OU model
mod.ou5<-gls(LuminanceMean.ventral ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou5<-dredge(global.model=mod.ou5)

write.table(d.ou5,"AICs_lumventralXvocaleco_OU_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.OU5<-coefTable(d.ou5) 
se.OU5.df<-do.call(rbind.data.frame,se.OU5 )
dim(se.OU5.df)
str(se.OU5.df)

write.table(se.OU5.df[1:7670,],"SE_lumventralXvocaleco_OU_male1.csv",sep = ",")
write.table(se.OU5.df[7671:15360,],"SE_lumventralXvocaleco_OU_male2.csv",sep = ",")



######
#Pagel's Lambda model
mod.pl5<-gls(LuminanceMean.ventral ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl5<-dredge(global.model=mod.pl5)

write.table(d.pl5,"AICs_lumventralXvocaleco_pl_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl5<-coefTable(d.pl5) 
se.pl5.df<-do.call(rbind.data.frame,se.pl5 )
dim(se.pl5.df)
str(se.pl5.df)

write.table(se.pl5.df[1:7670,],"SE_lumventralXvocaleco_pl_male1.csv",sep = ",")
write.table(se.pl5.df[7671:15360,],"SE_lumventralXvocaleco_pl_male2.csv",sep = ",")


###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#6. contrast.ventral 	

#BM model
mod.bro6<-gls(contrast.ventral.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro6<-dredge(global.model=mod.bro6)

write.table(d.bro6,"AICs_contventralXvocaleco_BM_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro6<-coefTable(d.bro6) 
se.bro6.df<-do.call(rbind.data.frame,se.bro6 )
dim(se.bro6.df)
str(se.bro6.df)

write.table(se.bro6.df[1:7669,],"SE_contventralXvocaleco_BM_male1.csv",sep = ",")
write.table(se.bro6.df[7670:15360,],"SE_contventralXvocaleco_BM_male2.csv",sep = ",")


#OU model
mod.ou6<-gls(contrast.ventral.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou6<-dredge(global.model=mod.ou6)

write.table(d.ou6,"AICs_contventralXvocaleco_OU_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.OU6<-coefTable(d.ou6) 
se.OU6.df<-do.call(rbind.data.frame,se.OU6 )
dim(se.OU6.df)
str(se.OU6.df)

write.table(se.OU6.df[1:7673,],"SE_contventralXvocaleco_OU_male1.csv",sep = ",")
write.table(se.OU6.df[7674:15360,],"SE_contventralXvocaleco_OU_male2.csv",sep = ",")


######
#Pagel's Lambda model
mod.pl6<-gls(contrast.ventral.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.7, phy = Form.tree), data = Form.male, method = "ML")

d.pl6<-dredge(global.model=mod.pl6)

write.table(d.pl6,"AICs_contventralXvocaleco_pl_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl6<-coefTable(d.pl6) 
se.pl6.df<-do.call(rbind.data.frame,se.pl6 )
dim(se.pl6.df)
str(se.pl6.df)

write.table(se.pl6.df[1:7673,],"SE_contventralXvocaleco_pl_male1.csv",sep = ",")
write.table(se.pl6.df[7674:15360,],"SE_contventralXvocaleco_pl_male2.csv",sep = ",")


###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#7. maxPower.wing.coverts	(BM) 

#BM model
mod.bro7<-gls(maxPower.wing.coverts ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro7<-dredge(global.model=mod.bro7)

write.table(d.bro7,"AICs_MPwingXvocaleco_BM_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro7<-coefTable(d.bro7) 
se.bro7.df<-do.call(rbind.data.frame,se.bro7 )
dim(se.bro7.df)
str(se.bro7.df)

write.table(se.bro7.df[1:7673,],"SE_MPwingXvocaleco_BM_male1.csv",sep = ",")
write.table(se.bro7.df[7674:15360,],"SE_MPwingXvocaleco_BM_male2.csv",sep = ",")


#OU model
mod.ou7<-gls(maxPower.wing.coverts ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou7<-dredge(global.model=mod.ou7)

write.table(d.ou7,"AICs_MPwingXvocaleco_OU_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.OU7<-coefTable(d.ou7) 
se.OU7.df<-do.call(rbind.data.frame,se.OU7 )
dim(se.OU7.df)
str(se.OU7.df)

write.table(se.OU7.df[1:7673,],"SE_MPwingXvocaleco_OU_male1.csv",sep = ",")
write.table(se.OU7.df[7674:15360,],"SE_MPwingXvocaleco_OU_male2.csv",sep = ",")


######
#Pagel's Lambda model
mod.pl7<-gls(maxPower.wing.coverts ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl7<-dredge(global.model=mod.pl7)

write.table(d.pl7,"AICs_MPwingXvocaleco_pl_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl7<-coefTable(d.pl7) 
se.pl7.df<-do.call(rbind.data.frame,se.pl7 )
dim(se.pl7.df)
str(se.pl7.df)

write.table(se.pl7.df[1:7673,],"SE_MPwingXvocaleco_pl_male1.csv",sep = ",")
write.table(se.pl7.df[7674:15360,],"SE_MPwingXvocaleco_pl_male2.csv",sep = ",")


###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#8. LuminanceMean.wing.coverts (BM)		

#BM model
mod.bro8<-gls(LuminanceMean.wing.coverts ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro8<-dredge(global.model=mod.bro8)

write.table(d.bro8,"AICs_lumwingXvocaleco_BM_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro8<-coefTable(d.bro8) 
se.bro8.df<-do.call(rbind.data.frame,se.bro8 )
dim(se.bro8.df)
str(se.bro8.df)

write.table(se.bro8.df[1:7668,],"SE_lumwingXvocaleco_BM_male1.csv",sep = ",")
write.table(se.bro8.df[7669:15360,],"SE_lumwingXvocaleco_BM_male2.csv",sep = ",")


#OU model
mod.ou8<-gls(LuminanceMean.wing.coverts ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou8<-dredge(global.model=mod.ou8)

write.table(d.ou8,"AICs_lumwingXvocaleco_OU_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.OU8<-coefTable(d.ou8) 
se.OU8.df<-do.call(rbind.data.frame,se.OU8 )
dim(se.OU8.df)
str(se.OU8.df)

write.table(se.OU8.df[1:7671,],"SE_lumwingXvocaleco_OU_male1.csv",sep = ",")
write.table(se.OU8.df[7672:15360,],"SE_lumwingXvocaleco_OU_male2.csv",sep = ",")


######
#Pagel's Lambda model
mod.pl8<-gls(LuminanceMean.wing.coverts ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.9, phy = Form.tree), data = Form.male, method = "ML")

d.pl8<-dredge(global.model=mod.pl8)

write.table(d.pl8,"AICs_lumwingXvocaleco_pl_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl8<-coefTable(d.pl8) 
se.pl8.df<-do.call(rbind.data.frame,se.pl8 )
dim(se.pl8.df)
str(se.pl8.df)

write.table(se.pl8.df[1:7671,],"SE_lumwingXvocaleco_pl_male1.csv",sep = ",")
write.table(se.pl8.df[7672:15360,],"SE_lumwingXvocaleco_pl_male2.csv",sep = ",")



#####################################
##################################
#Summary statistics of the best models 

#Wing coverts luminance X loudsong bandwidth >> all PL
#1a
mbest.pl8a<-gls(LuminanceMean.wing.coverts ~ Song.bandwidth + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.9, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.pl8a)

mbest.pl8b<-gls(LuminanceMean.wing.coverts ~ Song.bandwidth + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(0.9, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.pl8b)



###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#9. contrast.wing.coverts 	

#BM model
mod.bro9<-gls(contrast.wing.coverts.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro9<-dredge(global.model=mod.bro9)

write.table(d.bro9,"AICs_contwingXvocaleco_BM_male.csv",sep = ",",row.names = F)

se.bro9<-coefTable(d.bro9) 
se.bro9.df<-do.call(rbind.data.frame,se.bro9 )
dim(se.bro9.df)
str(se.bro9.df)

write.table(se.bro9.df[1:7669,],"SE_contwingXvocaleco_BM_male1.csv",sep = ",")
write.table(se.bro9.df[7670:15360,],"SE_contwingXvocaleco_BM_male2.csv",sep = ",")


#OU model
mod.ou9<-gls(contrast.wing.coverts.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou9<-dredge(global.model=mod.ou9)

write.table(d.ou9,"AICs_contwingXvocaleco_OU_male.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.OU9<-coefTable(d.ou9) 
se.OU9.df<-do.call(rbind.data.frame,se.OU9 )
dim(se.OU9.df)
str(se.OU9.df)

write.table(se.OU9.df[1:7675,],"SE_contwingXvocaleco_OU_male1.csv",sep = ",")
write.table(se.OU9.df[7676:15360,],"SE_contwingXvocaleco_OU_male2.csv",sep = ",")


######
#Pagel's Lambda model
mod.pl9<-gls(contrast.wing.coverts.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.5, phy = Form.tree), data = Form.male, method = "ML")

d.pl9<-dredge(global.model=mod.pl9)

write.table(d.pl9,"AICs_contwingXvocaleco_pl_male_1000.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl9<-coefTable(d.pl9) 
se.pl9.df<-do.call(rbind.data.frame,se.pl9 )
dim(se.pl9.df)
str(se.pl9.df)

write.table(se.pl9.df[1:7675,],"SE_contwingXvocaleco_pl_male1.csv",sep = ",")
write.table(se.pl9.df[7676:15360,],"SE_contwingXvocaleco_pl_male2.csv",sep = ",")


#####################################
##################################
#Summary statistics of the best models 

#Wing coverts contrast X loudsong duration >> all PL
#1a
mbest.pl9a<-gls(contrast.wing.coverts.inv ~ Delta.Time..s. + Habitat_exposure + Foraging_strata, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

summary(mbest.pl9a)

