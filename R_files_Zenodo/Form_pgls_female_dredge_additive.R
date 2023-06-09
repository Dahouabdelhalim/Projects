#PGLS plumage versus vocal traits in females, ecological and behavior additive effects


Form.female<-read.csv("Form_total_female_ok.csv",row.names = 1)
head(Form.female)


#Contrast values were inverted so that small values now represent lower contrast and higher  values represent higher contrast. (See Analyses of plumage data section for more info)
Form.female$contrast.dorsal.inv<-Form.female$contrast.dorsal*-1 #inversion of dorsal contrast
Form.female$contrast.ventral.inv<-Form.female$contrast.ventral*-1 #inversion of ventral contrast
Form.female$contrast.wing.coverts.inv<-Form.female$contrast.wing.coverts*-1 #inversion of wing contrast

head(Form.female)

library(tidyr)
#remove NAs
Form.female<-drop_na(Form.female)
head(Form.female)

#load packages
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
obj<-name.check(Form.tree,Form.female)
obj

#Remove tips that we don't have data
Form.tree<- drop.tip(Form.tree, obj$tree_not_data)
plot(Form.tree)

#Check tree and data again
name.check(Form.tree, Form.female)

#Branch length transformation. Necessary to avoid error you were getting before. Does not affect results. 

Form.tree$edge.length <- Form.tree$edge.length * 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000 


####################################################################################################################################################################################################################################################################################
#Models 

#1. maxPower.dorsal	(BM)
#versus Note.count (BM) + Note.type (BM) + Delta.Time..s. (OU) + Peak.Freq..Hz. (BM) +	Note.diversity (BM) + Note.rate	(BM) + Song.bandwidth (OU) + Loudsong.mod.rate (OU)

#BM model
mod.bro1<-gls(maxPower.dorsal ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro1<-dredge(global.model=mod.bro1)

write.table(d.bro1,"AICs_MPdorsalXvocaleco_BM_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro1<-coefTable(d.bro1) 
se.bro1.df<-do.call(rbind.data.frame,se.bro1 )
dim(se.bro1.df)
str(se.bro1.df)

write.table(se.bro1.df[1:7675,],"SE_MPdorsalXvocaleco_BM_female1.csv",sep = ",")
write.table(se.bro1.df[7676:15360,],"SE_MPdorsalXvocaleco_BM_female2.csv",sep = ",")

######
#OU model
mod.ou1<-gls(maxPower.dorsal ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou1<-dredge(global.model=mod.ou1)

write.table(d.ou1,"AICs_MPdorsalXvocaleco_OU_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.ou1<-coefTable(d.ou1) 
head(se.ou1)
se.ou1.df<-do.call(rbind.data.frame,se.ou1 )
dim(se.ou1.df)
str(se.ou1.df)

write.table(se.ou1.df[1:7683,],"SE_MPdorsalXvocaleco_OU_female1.csv",sep = ",")
write.table(se.ou1.df[7684:15360,],"SE_MPdorsalXvocaleco_OU_female2.csv",sep = ",")


######
#Pagel's Lambda model  
mod.pl1<-gls(maxPower.dorsal ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

d.pl1<-dredge(global.model=mod.pl1)
d.pl1


write.table(d.pl1,"AICs_MPdorsalXvocaleco_pl_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl1<-coefTable(d.pl1) 
head(se.pl1)
se.pl1.df<-do.call(rbind.data.frame,se.pl1 )
dim(se.pl1.df)
str(se.pl1.df)

write.table(se.pl1.df[1:7683,],"SE_MPdorsalXvocaleco_pl_female1.csv",sep = ",")
write.table(se.pl1.df[7684:15360,],"SE_MPdorsalXvocaleco_pl_female2.csv",sep = ",")


#####################################
##################################
#Summary statistics of the best models 

#Dorsal maximum energy X loudsong bandwidth >> all BM PL BM
#1a
f_mbest.bm1a<-gls(maxPower.dorsal ~ Song.bandwidth + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm1a)

#1b
f_mbest.pl1b<-gls(maxPower.dorsal ~ Song.bandwidth + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.pl1b)

#1c
f_mbest.bm1c<-gls(maxPower.dorsal ~ Song.bandwidth + Habitat_exposure + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm1c)

#Dorsal maximum energy X Loudsong.mod.rate >> all OU BM OU
#1a
f_mbest.ou11a<-gls(maxPower.dorsal ~ Loudsong.mod.rate + Habitat_exposure + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

summary(f_mbest.ou11a)

#11b
f_mbest.bm11b<-gls(maxPower.dorsal ~ Loudsong.mod.rate + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm11b)

##11c
f_mbest.ou11c<-gls(maxPower.dorsal ~ Loudsong.mod.rate + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

summary(f_mbest.ou11c)

#Dorsal maximum energy X note diversity >> all BM 
#111a
f_mbest.bm111a<-gls(maxPower.dorsal ~ Note.diversity + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm111a)

#111b
f_mbest.bm111b<-gls(maxPower.dorsal ~ Note.diversity + Habitat_exposure  + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm111b)

#111c
f_mbest.bm111c<-gls(maxPower.dorsal ~ Note.diversity + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm111c)


#Dorsal maximum energy X note type >> all BM 
#1111a
f_mbest.bm1111a<-gls(maxPower.dorsal ~ Note.type + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm1111a)

###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#2.	LuminanceMean.dorsal  
#versus Note.count (BM) + Note.type (BM) + Delta.Time..s. (OU) + Peak.Freq..Hz. (BM) +	Note.diversity (BM) + Note.rate	(BM) + Song.bandwidth (OU) + Loudsong.mod.rate (OU)

#BM model
mod.bro2<-gls(LuminanceMean.dorsal ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro2<-dredge(global.model=mod.bro2)

write.table(d.bro2,"AICs_LumdorsalXvocaleco_BM_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro2<-coefTable(d.bro2) 
se.bro2.df<-do.call(rbind.data.frame,se.bro2 )
dim(se.bro2.df)
str(se.bro2.df)

write.table(se.bro2.df[1:7676,],"SE_LumdorsalXvocaleco_BM_female1.csv",sep = ",")
write.table(se.bro2.df[7677:15360,],"SE_LumdorsalXvocaleco_BM_female2.csv",sep = ",")


#OU model
mod.ou2<-gls(LuminanceMean.dorsal ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou2<-dredge(global.model=mod.ou2)

write.table(d.ou2,"AICs_LumdorsalXvocaleco_OU_female.csv",sep = ",",row.names = F)

se.ou2<-coefTable(d.ou2) 
se.ou2.df<-do.call(rbind.data.frame,se.ou2 )
dim(se.ou2.df)
str(se.ou2.df)

write.table(se.ou2.df[1:7679,],"SE_LumdorsalXvocaleco_OU_female1.csv",sep = ",")
write.table(se.ou2.df[7680:15360,],"SE_LumdorsalXvocaleco_OU_female2.csv",sep = ",")


#####
#Pagel's Lambda model  
mod.pl2<-gls(LuminanceMean.dorsal ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

d.pl2<-dredge(global.model=mod.pl2)

write.table(d.pl2,"AICs_LumdorsalXvocaleco_pl_female10000.csv",sep = ",",row.names = F)


#obtain SE (standard error)
se.pl2<-coefTable(d.pl2) 
se.pl2.df<-do.call(rbind.data.frame,se.pl2 )
dim(se.pl2.df)
str(se.pl2.df)

write.table(se.pl2.df[1:7679,],"SE_LumdorsalXvocaleco_pl_female1.csv",sep = ",")
write.table(se.pl2.df[7680:15360,],"SE_LumdorsalXvocaleco_pl_female2.csv",sep = ",")

#######
#Run models separately that did not run in the dredge
mod.pl2<-gls(LuminanceMean.dorsal ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

#Only duration
mod.pl2_dur1<-gls(LuminanceMean.dorsal ~ Delta.Time..s., correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_dur1)
AICc(mod.pl2_dur1)

mod.pl2_dur2<-gls(LuminanceMean.dorsal ~ Delta.Time..s. + Habitat_exposure, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_dur2)


mod.pl2_dur3<-gls(LuminanceMean.dorsal ~ Delta.Time..s. + Foraging_strata, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_dur3)
AICc(mod.pl2_dur3)

mod.pl2_dur4<-gls(LuminanceMean.dorsal ~ Delta.Time..s. + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_dur4)


mod.pl2_dur5<-gls(LuminanceMean.dorsal ~ Delta.Time..s. + Habitat_exposure + Foraging_strata, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_dur5)
AICc(mod.pl2_dur5)

mod.pl2_dur6<-gls(LuminanceMean.dorsal ~ Delta.Time..s. + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_dur6)


mod.pl2_dur7<-gls(LuminanceMean.dorsal ~ Delta.Time..s. + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_dur7)
AICc(mod.pl2_dur7)

mod.pl2_dur8<-gls(LuminanceMean.dorsal ~ Delta.Time..s. + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_dur8)

#Only Loudsong.mod.rate >> FOR + MSF model
mod.pl2_fs7<-gls(LuminanceMean.dorsal ~ Loudsong.mod.rate + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_fs7)
AICc(mod.pl2_fs7)


#Only note count >>  HAB + FOR
mod.pl2_nc5<-gls(LuminanceMean.dorsal ~ Note.count + Habitat_exposure + Foraging_strata, correlation  = corPagel(0.8, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_nc5)
AICc(mod.pl2_nc5)

#Only note diversity >>  HAB + FOR
mod.pl2_nd5<-gls(LuminanceMean.dorsal ~ Note.diversity + Habitat_exposure + Foraging_strata, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_nd5)
AICc(mod.pl2_nd5)

#Only note rate >>  + HAB
mod.pl2_nr2<-gls(LuminanceMean.dorsal ~ Note.rate + Habitat_exposure, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_nr2)
AICc(mod.pl2_nr2)

#Only note type >> 
mod.pl2_nt1<-gls(LuminanceMean.dorsal ~ Note.type, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_nt1)
AICc(mod.pl2_nt1)

#Only loudsong bandwidth >> HAB + FOR
mod.pl2_sb5<-gls(LuminanceMean.dorsal ~ Song.bandwidth + Habitat_exposure + Foraging_strata, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl2_sb5)
AICc(mod.pl2_sb5)

#####################################
##################################
#Summary statistics of the best models 

#Dorsal luminance X loudsong duration >> all PL
#2a
f_mbest.pl2a<-gls(LuminanceMean.dorsal ~ Delta.Time..s., correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.pl2a)

#2b
f_mbest.pl2b<-gls(LuminanceMean.dorsal ~ Delta.Time..s. + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.pl2b)


#Dorsal luminance X loudsong bandwidth >> all PL
#22a
f_mbest.pl22a<-gls(LuminanceMean.dorsal ~ Song.bandwidth, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.pl22a)


#Dorsal luminance X Loudsong.mod.rate >> OU BM OU
#222a
f_mbest.ou222a<-gls(LuminanceMean.dorsal ~ Loudsong.mod.rate, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

summary(f_mbest.ou222a)

#222b
f_mbest.bm222b<-gls(LuminanceMean.dorsal ~ Loudsong.mod.rate + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm222b)

#222c
f_mbest.ou222c<-gls(LuminanceMean.dorsal ~ Loudsong.mod.rate + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

summary(f_mbest.ou222c)

#Dorsal luminance X note count >> all PL
#2222a
f_mbest.pl2222a<-gls(LuminanceMean.dorsal ~ Note.count, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.pl2222a)

#2222b
f_mbest.pl2222b<-gls(LuminanceMean.dorsal ~ Note.count + Habitat_exposure, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.pl2222b)



###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#3. contrast.dorsal	(BM)

#BM
mod.bro3<-gls(contrast.dorsal.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro3<-dredge(global.model=mod.bro3)

write.table(d.bro3,"AICs_ContdorsalXvocaleco_BM_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro3<-coefTable(d.bro3) 
se.bro3.df<-do.call(rbind.data.frame,se.bro3 )
dim(se.bro3.df)
str(se.bro3.df)

write.table(se.bro3.df[1:7669,],"SE_ContdorsalXvocaleco_BM_female1.csv",sep = ",")
write.table(se.bro3.df[7670:15360,],"SE_ContdorsalXvocaleco_BM_female2.csv",sep = ",")



#OU model
mod.ou3<-gls(contrast.dorsal.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou3<-dredge(global.model=mod.ou3)

write.table(d.ou3,"AICs_ContdorsalXvocaleco_OU_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.ou3<-coefTable(d.ou3) 
se.ou3.df<-do.call(rbind.data.frame,se.ou3)
dim(se.ou3.df)
str(se.ou3.df)

write.table(se.ou3.df[1:7673,],"SE_ContdorsalXvocaleco_OU_female1.csv",sep = ",")
write.table(se.ou3.df[7674:15360,],"SE_ContdorsalXvocaleco_OU_female2.csv",sep = ",")


######
#Pagel's Lambda model   
mod.pl3<-gls(contrast.dorsal.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

d.pl3<-dredge(global.model=mod.pl3)

write.table(d.pl3,"AICs_ContdorsalXvocaleco_pl_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl3<-coefTable(d.pl3) 
se.pl3.df<-do.call(rbind.data.frame,se.pl3)
dim(se.pl3.df)
str(se.pl3.df)

write.table(se.pl3.df[1:7673,],"SE_ContdorsalXvocaleco_pl_female1.csv",sep = ",")
write.table(se.pl3.df[7674:15360,],"SE_ContdorsalXvocaleco_pl_female2.csv",sep = ",")


#####################################
##################################
#Summary statistics of the best models 

#Dorsal contrast X peak frequency >> all BM
#3a
f_mbest.bm3a<-gls(contrast.dorsal.inv ~ Peak.Freq..Hz. + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm3a)


#Dorsal contrast X Loudsong.mod.rate >> all BM
#33a
f_mbest.bm33a<-gls(contrast.dorsal.inv ~ Loudsong.mod.rate + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm33a)


#Dorsal contrast X note count >> all BM
#333a
f_mbest.bm333a<-gls(contrast.dorsal.inv ~ Note.count + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm333a)

#333b
f_mbest.bm333b<-gls(contrast.dorsal.inv ~ Note.count + Habitat_exposure + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm333b)

#333c
f_mbest.bm333c<-gls(contrast.dorsal.inv ~ Note.count + Habitat_exposure + Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm333c)

#333d
f_mbest.bm333d<-gls(contrast.dorsal.inv ~ Note.count + Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm333d)


#Dorsal contrast X note diversity>> all BM
#nda
f_mbest.bm.nda<-gls(contrast.dorsal.inv ~ Note.diversity + Habitat_exposure + Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm.nda)

#ndb
f_mbest.bm.ndb<-gls(contrast.dorsal.inv ~ Note.diversity + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm.ndb)

#ndc
f_mbest.bm.ndc<-gls(contrast.dorsal.inv ~ Note.diversity + Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm.ndc)

#ndd
f_mbest.bm.ndd<-gls(contrast.dorsal.inv ~ Note.diversity + Habitat_exposure + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm.ndd)

###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#4. maxPower.ventral (BM)		

#BM model
mod.bro4<-gls(maxPower.ventral ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro4<-dredge(global.model=mod.bro4)

write.table(d.bro4,"AICs_MPventralXvocaleco_BM_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro4<-coefTable(d.bro4) 
se.bro4.df<-do.call(rbind.data.frame,se.bro4 )
dim(se.bro4.df)
str(se.bro4.df)

write.table(se.bro4.df[1:7674,],"SE_MPventralXvocaleco_BM_female1.csv",sep = ",")
write.table(se.bro4.df[7675:15360,],"SE_MPventralXvocaleco_BM_female2.csv",sep = ",")


#OU model
mod.ou4<-gls(maxPower.ventral ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou4<-dredge(global.model=mod.ou4)

write.table(d.ou4,"AICs_MPventralXvocaleco_OU_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.ou4<-coefTable(d.ou4) 
se.ou4.df<-do.call(rbind.data.frame,se.ou4 )
dim(se.ou4.df)
str(se.ou4.df)

write.table(se.ou4.df[1:7666,],"SE_MPventralXvocaleco_OU_female1.csv",sep = ",")
write.table(se.ou4.df[7667:15360,],"SE_MPventralXvocaleco_OU_female2.csv",sep = ",")


########
#Pagel's Lambda model  
mod.pl4<-gls(maxPower.ventral ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.6, phy = Form.tree), data = Form.female, method = "ML")


d.pl4<-dredge(global.model=mod.pl4)

write.table(d.pl4,"AICs_MPventralXvocaleco_pl_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl4<-coefTable(d.pl4) 
se.pl4.df<-do.call(rbind.data.frame,se.pl4 )
dim(se.pl4.df)
str(se.pl4.df)

write.table(se.pl4.df[1:7666,],"SE_MPventralXvocaleco_pl_female1.csv",sep = ",")
write.table(se.pl4.df[7667:15360,],"SE_MPventralXvocaleco_pl_female2.csv",sep = ",")

#######
#Run models separately that did not run in the dredge
#Only duration >> Hab + For and Hab + For + MSF
mod.pl4_dur5<-gls(maxPower.ventral ~ Delta.Time..s. + Habitat_exposure + Foraging_strata, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl4_dur5)
AICc(mod.pl4_dur5)

mod.pl4_dur8<-gls(maxPower.ventral ~ Delta.Time..s. + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl4_dur8)
AICc(mod.pl4_dur8)


#Note count >> Habitat exposure + Foraging strata ; Habitat exposure + MSF ; Habitat exposure + Foraging strata + MSF
mod.pl4_nc5<-gls(maxPower.ventral ~ Note.count + Habitat_exposure + Foraging_strata, correlation  = corPagel(0.5, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl4_nc5)
AICc(mod.pl4_nc5)

mod.pl4_nc6<-gls(maxPower.ventral ~ Note.count + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl4_nc6)
AICc(mod.pl4_nc6)

mod.pl4_nc8<-gls(maxPower.ventral ~ Note.count + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl4_nc8)
AICc(mod.pl4_nc8)


##Only note diversity >> 
mod.pl4_nt8<-gls(maxPower.ventral ~ Note.diversity + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl4_nt8)
AICc(mod.pl4_nt8)


#Only note rate >>  
mod.pl4_nr2<-gls(maxPower.ventral ~ Note.rate + Habitat_exposure, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl4_nr2)
AICc(mod.pl4_nr2)

mod.pl4_nr3<-gls(maxPower.ventral ~ Note.rate + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl4_nr3)
AICc(mod.pl4_nr3)

mod.pl4_nr4<-gls(maxPower.ventral ~ Note.rate + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl4_nr4)
AICc(mod.pl4_nr4)

mod.pl4_nr5<-gls(maxPower.ventral ~ Note.rate + Habitat_exposure + Foraging_strata, correlation  = corPagel(0, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl4_nr5)
AICc(mod.pl4_nr5)

mod.pl4_nr6<-gls(maxPower.ventral ~ Note.rate + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(0.3, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl4_nr6)
AICc(mod.pl4_nr6)



#####################################
##################################
#Summary statistics of the best models 

#Ventral mp X peak frequency >> all PL
#4a
f_mbest.pl4a<-gls(maxPower.ventral ~ Peak.Freq..Hz., correlation  = corPagel(0.6, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.pl4a)

#Ventral mp X note rate >> all PL
#44a
f_mbest.pl44a<-gls(maxPower.ventral ~ Note.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.6, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.pl44a)



###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#5. LuminanceMean.ventral	(BM)

#BM model
mod.bro5<-gls(LuminanceMean.ventral ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro5<-dredge(global.model=mod.bro5)

write.table(d.bro5,"AICs_lumventralXvocaleco_BM_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro5<-coefTable(d.bro5) 
se.bro5.df<-do.call(rbind.data.frame,se.bro5 )
dim(se.bro5.df)
str(se.bro5.df)

write.table(se.bro5.df[1:7668,],"SE_lumventralXvocaleco_BM_female1.csv",sep = ",")
write.table(se.bro5.df[7669:15360,],"SE_lumventralXvocaleco_BM_female2.csv",sep = ",")



#OU model
mod.ou5<-gls(LuminanceMean.ventral ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou5<-dredge(global.model=mod.ou5)

write.table(d.ou5,"AICs_lumventralXvocaleco_OU_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.OU5<-coefTable(d.ou5) 
se.OU5.df<-do.call(rbind.data.frame,se.OU5 )
dim(se.OU5.df)
str(se.OU5.df)

write.table(se.OU5.df[1:7669,],"SE_lumventralXvocaleco_OU_female1.csv",sep = ",")
write.table(se.OU5.df[7670:15360,],"SE_lumventralXvocaleco_OU_female2.csv",sep = ",")


######
#Pagel's Lambda model  
mod.pl5<-gls(LuminanceMean.ventral ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

d.pl5<-dredge(global.model=mod.pl5)

write.table(d.pl5,"AICs_lumventralXvocaleco_pl_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl5<-coefTable(d.pl5) 
se.pl5.df<-do.call(rbind.data.frame,se.pl5 )
dim(se.pl5.df)
str(se.pl5.df)

write.table(se.pl5.df[1:7669,],"SE_lumventralXvocaleco_pl_female1.csv",sep = ",")
write.table(se.pl5.df[7670:15360,],"SE_lumventralXvocaleco_pl_female2.csv",sep = ",")


#####################################
##################################
#Summary statistics of the best models 

#Ventral luminance X peak frequency >> all PL
#5a
f_mbest.bm5a<-gls(LuminanceMean.ventral ~ Peak.Freq..Hz., correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm5a)

#5b
f_mbest.bm5b<-gls(LuminanceMean.ventral ~ Peak.Freq..Hz. + Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm5b)

#5c
f_mbest.bm5c<-gls(LuminanceMean.ventral ~ Peak.Freq..Hz. + Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm5c)




###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#6. contrast.ventral (OU)  

#BM model
mod.bro6<-gls(contrast.ventral.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro6<-dredge(global.model=mod.bro6)

write.table(d.bro6,"AICs_contventralXvocaleco_BM_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro6<-coefTable(d.bro6) 
se.bro6.df<-do.call(rbind.data.frame,se.bro6 )
dim(se.bro6.df)
str(se.bro6.df)

write.table(se.bro6.df[1:7673,],"SE_contventralXvocaleco_BM_female1.csv",sep = ",")
write.table(se.bro6.df[7674:15360,],"SE_contventralXvocaleco_BM_female2.csv",sep = ",")


#OU model
mod.ou6<-gls(contrast.ventral.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou6<-dredge(global.model=mod.ou6)

write.table(d.ou6,"AICs_contventralXvocaleco_OU_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.OU6<-coefTable(d.ou6) 
se.OU6.df<-do.call(rbind.data.frame,se.OU6 )
dim(se.OU6.df)
str(se.OU6.df)

write.table(se.OU6.df[1:7673,],"SE_contventralXvocaleco_OU_female1.csv",sep = ",")
write.table(se.OU6.df[7674:15360,],"SE_contventralXvocaleco_OU_female2.csv",sep = ",")



#Pagel's Lambda   
mod.pl6<-gls(contrast.ventral.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")


d.pl6<-dredge(global.model=mod.pl6)

write.table(d.pl6,"AICs_contventralXvocaleco_pl_female_sb.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl6<-coefTable(d.pl6) 
se.pl6.df<-do.call(rbind.data.frame,se.pl6 )
dim(se.pl6.df)
str(se.pl6.df)

write.table(se.pl6.df[1:7673,],"SE_contventralXvocaleco_pl_female1.csv",sep = ",")
write.table(se.pl6.df[7674:15360,],"SE_contventralXvocaleco_pl_female2.csv",sep = ",")


#Run models separately that did not run in the dredge
#Only loudsongbandwidth >>> 
mod.pl6_sb4<-gls(contrast.ventral.inv ~ Song.bandwidth + Mixed.Flocking, correlation  = corPagel(0.8, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl6_sb4)
AICc(mod.pl6_sb4)

mod.pl6_sb6<-gls(contrast.ventral.inv ~ Song.bandwidth + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(0.8, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl6_sb6)
AICc(mod.pl6_sb6)

mod.pl6_sb7<-gls(contrast.ventral.inv ~ Song.bandwidth + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.8, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl6_sb7)
AICc(mod.pl6_sb7)

mod.pl6_sb8<-gls(contrast.ventral.inv ~ Song.bandwidth + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl6_sb8)
AICc(mod.pl6_sb8)


#####################################
##################################
#Summary statistics of the best models 

#Ventral contrast X Loudsong.mod.rate >> all BM
#6a
f_mbest.bm6a<-gls(contrast.ventral.inv ~  Loudsong.mod.rate, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm6a)

#Ventral contrast X note count >> all BM
#66a
f_mbest.bm66a<-gls(contrast.ventral.inv ~ Note.count, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm66a)

#Ventral contrast X note diversity >> BM OU
#666a
f_mbest.bm666a<-gls(contrast.ventral.inv ~ Note.diversity, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm666a)


#666b
f_mbest.ou666b<-gls(contrast.ventral.inv ~ Note.diversity, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

summary(f_mbest.ou666b)


###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#7. maxPower.wing.coverts	(BM) 

#BM model
mod.bro7<-gls(maxPower.wing.coverts ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro7<-dredge(global.model=mod.bro7)

write.table(d.bro7,"AICs_MPwingXvocaleco_BM_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro7<-coefTable(d.bro7) 
se.bro7.df<-do.call(rbind.data.frame,se.bro7 )
dim(se.bro7.df)
str(se.bro7.df)

write.table(se.bro7.df[1:7671,],"SE_MPwingXvocaleco_BM_female1.csv",sep = ",")
write.table(se.bro7.df[7672:15360,],"SE_MPwingXvocaleco_BM_female2.csv",sep = ",")


#OU model
mod.ou7<-gls(maxPower.wing.coverts ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou7<-dredge(global.model=mod.ou7)

write.table(d.ou7,"AICs_MPwingXvocaleco_OU_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.OU7<-coefTable(d.ou7) 
se.OU7.df<-do.call(rbind.data.frame,se.OU7 )
dim(se.OU7.df)
str(se.OU7.df)

write.table(se.OU7.df[1:7670,],"SE_MPwingXvocaleco_OU_female1.csv",sep = ",")
write.table(se.OU7.df[7671:15360,],"SE_MPwingXvocaleco_OU_female2.csv",sep = ",")


######
#Pagel's Lambda model  
mod.pl7<-gls(maxPower.wing.coverts ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

d.pl7<-dredge(global.model=mod.pl7)

write.table(d.pl7,"AICs_MPwingXvocaleco_pl_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl7<-coefTable(d.pl7) 
se.pl7.df<-do.call(rbind.data.frame,se.pl7 )
dim(se.pl7.df)
str(se.pl7.df)

write.table(se.pl7.df[1:7670,],"SE_MPwingXvocaleco_pl_female1.csv",sep = ",")
write.table(se.pl7.df[7671:15360,],"SE_MPwingXvocaleco_pl_female2.csv",sep = ",")


#####################################
##################################
#Summary statistics of the best models 

#Wing Maximum energy X peak frequency >> all BM
#7a
f_mbest.bm7a<-gls(maxPower.wing.coverts ~ Peak.Freq..Hz., correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm7a)



###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#8. LuminanceMean.wing.coverts 

#BM model
mod.bro8<-gls(LuminanceMean.wing.coverts ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro8<-dredge(global.model=mod.bro8)

write.table(d.bro8,"AICs_lumwingXvocaleco_BM_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro8<-coefTable(d.bro8) 
se.bro8.df<-do.call(rbind.data.frame,se.bro8 )
dim(se.bro8.df)
str(se.bro8.df)

write.table(se.bro8.df[1:7666,],"SE_lumwingXvocaleco_BM_female1.csv",sep = ",")
write.table(se.bro8.df[7667:15360,],"SE_lumwingXvocaleco_BM_female2.csv",sep = ",")



#OU model
mod.ou8<-gls(LuminanceMean.wing.coverts ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou8<-dredge(global.model=mod.ou8)

write.table(d.ou8,"AICs_lumwingXvocaleco_OU_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.OU8<-coefTable(d.ou8) 
se.OU8.df<-do.call(rbind.data.frame,se.OU8 )
dim(se.OU8.df)
str(se.OU8.df)

write.table(se.OU8.df[1:7669,],"SE_lumwingXvocaleco_OU_female1.csv",sep = ",")
write.table(se.OU8.df[7670:15360,],"SE_lumwingXvocaleco_OU_female2.csv",sep = ",")



#Pagel's Lambda lambda 
mod.pl8<-gls(LuminanceMean.wing.coverts ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML")

d.pl8<-dredge(global.model=mod.pl8)

write.table(d.pl8,"AICs_lumwingXvocaleco_pl_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl8<-coefTable(d.pl8) 
se.pl8.df<-do.call(rbind.data.frame,se.pl8 )
dim(se.pl8.df)
str(se.pl8.df)

write.table(se.pl8.df[1:7669,],"SE_lumwingXvocaleco_pl_female1.csv",sep = ",")
write.table(se.pl8.df[7670:15360,],"SE_lumwingXvocaleco_pl_female2.csv",sep = ",")

#Run models separately that did not run in the dredge
#Only duration >> 
mod.pl8_dur6<-gls(LuminanceMean.wing.coverts ~ Delta.Time..s. + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(0.8, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_dur6)
AICc(mod.pl8_dur6)

mod.pl8_dur7<-gls(LuminanceMean.wing.coverts ~ Delta.Time..s. + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_dur7)
AICc(mod.pl8_dur7)

mod.pl8_dur8<-gls(LuminanceMean.wing.coverts ~ Delta.Time..s. + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.4, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_dur8)
AICc(mod.pl8_dur8)

#Loudsong.mod.rate 
mod.pl8_fs3<-gls(LuminanceMean.wing.coverts ~ Loudsong.mod.rate + Foraging_strata, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_fs3)
AICc(mod.pl8_fs3)

mod.pl8_fs5<-gls(LuminanceMean.wing.coverts ~ Loudsong.mod.rate + Habitat_exposure + Foraging_strata, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_fs5)
AICc(mod.pl8_fs5)

mod.pl8_fs6<-gls(LuminanceMean.wing.coverts ~ Loudsong.mod.rate + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_fs6)
AICc(mod.pl8_fs6)

mod.pl8_fs7<-gls(LuminanceMean.wing.coverts ~ Loudsong.mod.rate + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_fs7)
AICc(mod.pl8_fs7)

mod.pl8_fs8<-gls(LuminanceMean.wing.coverts ~ Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_fs8)
AICc(mod.pl8_fs8)


#Only note diversity >> model all covs
mod.pl8_nt8<-gls(LuminanceMean.wing.coverts ~ Note.diversity + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_nt8)
AICc(mod.pl8_nt8)


#Only note rate
mod.pl8_nr5<-gls(LuminanceMean.wing.coverts ~ Note.rate + Habitat_exposure + Foraging_strata, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_nr5)
AICc(mod.pl8_nr5)

mod.pl8_nr6<-gls(LuminanceMean.wing.coverts ~ Note.rate + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_nr6)
AICc(mod.pl8_nr6)

mod.pl8_nr8<-gls(LuminanceMean.wing.coverts ~ Note.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_nr8)
AICc(mod.pl8_nr8)

#Only Note type >>
mod.pl8_nt7<-gls(LuminanceMean.wing.coverts ~ Note.type + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_nt7)
AICc(mod.pl8_nt7)

#Only peak freq
mod.pl8_pf6<-gls(LuminanceMean.wing.coverts ~ Peak.Freq..Hz. + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_pf6)
AICc(mod.pl8_pf6)

mod.pl8_pf8<-gls(LuminanceMean.wing.coverts ~ Peak.Freq..Hz. + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_pf8)
AICc(mod.pl8_pf8)


#Only loudsong bandwidth
mod.pl8_sb6<-gls(LuminanceMean.wing.coverts ~ Song.bandwidth + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_sb6)
AICc(mod.pl8_sb6)

mod.pl8_sb8<-gls(LuminanceMean.wing.coverts ~ Song.bandwidth + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.6, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl8_sb8)
AICc(mod.pl8_sb8)

#####################################
##################################
#Summary statistics of the best models 

#Wing luminance X note diversity >> all PL
#8a
f_mbest.bm8a<-gls(LuminanceMean.wing.coverts ~ Note.diversity + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm8a)

#8b
f_mbest.bm8b<-gls(LuminanceMean.wing.coverts ~ Note.diversity + Foraging_strata, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm8b)


#Wing luminance X note rate >> all PL
#88a
f_mbest.bm88a<-gls(LuminanceMean.wing.coverts ~ Note.rate + Foraging_strata, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML")

summary(f_mbest.bm88a)


###########################################################################################
#-----------------------------------------------------------------------------------------#
###########################################################################################

#9. contrast.wing.coverts 

#BM model
mod.bro9<-gls(contrast.wing.coverts.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro9<-dredge(global.model=mod.bro9)

write.table(d.bro9,"AICs_contwingXvocaleco_BM_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.bro9<-coefTable(d.bro9) 
se.bro9.df<-do.call(rbind.data.frame,se.bro9 )
dim(se.bro9.df)
str(se.bro9.df)

write.table(se.bro9.df[1:7670,],"SE_contwingXvocaleco_BM_female1.csv",sep = ",")
write.table(se.bro9.df[7671:15360,],"SE_contwingXvocaleco_BM_female2.csv",sep = ",")



#OU model
mod.ou9<-gls(contrast.wing.coverts.inv ~ Note.count + Note.type + Delta.Time..s. + Peak.Freq..Hz. + Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou9<-dredge(global.model=mod.ou9)

write.table(d.ou9,"AICs_contwingXvocaleco_OU_female.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.OU9<-coefTable(d.ou9) 
se.OU9.df<-do.call(rbind.data.frame,se.OU9 )
dim(se.OU9.df)
str(se.OU9.df)

write.table(se.OU9.df[1:7675,],"SE_contwingXvocaleco_OU_female1.csv",sep = ",")
write.table(se.OU9.df[7676:15360,],"SE_contwingXvocaleco_OU_female2.csv",sep = ",")



#Pagel's Lambda   

#Run models separately that did not run in the dredge

#Part a -->Note type , loudsong duration e peak frequency
mod.pl9a<-gls(contrast.wing.coverts.inv ~ Note.type + Delta.Time..s. + Peak.Freq..Hz. + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")


d.pl9a<-dredge(global.model=mod.pl9a)

write.table(d.pl9a,"AICs_contwingXvocaleco_pl_femalea.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl9a<-coefTable(d.pl9a) 
se.pl9a.df<-do.call(rbind.data.frame,se.pl9a )
dim(se.pl9a.df)
str(se.pl9a.df)

write.table(se.pl9a.df[1:7675,],"SE_contwingXvocaleco_pl_female1a.csv",sep = ",")
write.table(se.pl9a.df[7676:15360,],"SE_contwingXvocaleco_pl_female2a.csv",sep = ",")

###########
#Part b --> Note diversity, note rate, song bandwidth e Loudsong.mod.rate
mod.pl9b<-gls(contrast.wing.coverts.inv ~ Note.diversity + Note.rate + Song.bandwidth + Loudsong.mod.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML")

d.pl9b<-dredge(global.model=mod.pl9b)

write.table(d.pl9b,"AICs_contwingXvocaleco_pl_femaleb.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.pl9b<-coefTable(d.pl9b) 
se.pl9b.df<-do.call(rbind.data.frame,se.pl9b )
dim(se.pl9b.df)
str(se.pl9b.df)

write.table(se.pl9b.df[1:7675,],"SE_contwingXvocaleco_pl_female1b.csv",sep = ",")


#################
#Only note count

mod.pl9_nc1<-gls(contrast.wing.coverts.inv ~ Note.count, correlation  = corPagel(0.3, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl9_nc1)
AICc(mod.pl9_nc1)

mod.pl9_nc2<-gls(contrast.wing.coverts.inv ~ Note.count + Habitat_exposure, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl9_nc2)
AICc(mod.pl9_nc2)

mod.pl9_nc3<-gls(contrast.wing.coverts.inv ~ Note.count + Foraging_strata, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl9_nc3)
AICc(mod.pl9_nc3)

mod.pl9_nc4<-gls(contrast.wing.coverts.inv ~ Note.count + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl9_nc4)
AICc(mod.pl9_nc4)

mod.pl9_nc5<-gls(contrast.wing.coverts.inv ~ Note.count + Habitat_exposure + Foraging_strata, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl9_nc5)
AICc(mod.pl9_nc5)


mod.pl9_nc6<-gls(contrast.wing.coverts.inv ~ Note.count + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(0.2, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl9_nc6)
AICc(mod.pl9_nc6)


mod.pl9_nc7<-gls(contrast.wing.coverts.inv ~ Note.count + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.1, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl9_nc7)
AICc(mod.pl9_nc7)


mod.pl9_nc8<-gls(contrast.wing.coverts.inv ~ Note.count + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML") 
summary(mod.pl9_nc8)
AICc(mod.pl9_nc8)


###################
#Only note diversity 
mod.pl9_nd1<-gls(contrast.wing.coverts.inv ~ Note.diversity, correlation  = corPagel(0, phy = Form.tree), data = Form.female, method = "ML") 
summary(mod.pl9_nd1)
AICc(mod.pl9_nd1)

mod.pl9_nd2<-gls(contrast.wing.coverts.inv ~ Note.diversity + Habitat_exposure, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML") 
summary(mod.pl9_nd2)
AICc(mod.pl9_nd2)
mod.pl9_nd2<-gls(contrast.wing.coverts.inv ~ Note.diversity + Habitat_exposure, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML")


mod.pl9_nd3<-gls(contrast.wing.coverts.inv ~ Note.diversity + Foraging_strata, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML") 
summary(mod.pl9_nd3)
AICc(mod.pl9_nd3)


mod.pl9_nd4<-gls(contrast.wing.coverts.inv ~ Note.diversity + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML") 

summary(mod.pl9_nd4)
AICc(mod.pl9_nd4)


mod.pl9_nd5<-gls(contrast.wing.coverts.inv ~ Note.diversity + Habitat_exposure + Foraging_strata, correlation  = corPagel(0.2, phy = Form.tree), data = Form.female, method = "ML") 
summary(mod.pl9_nd5)
AICc(mod.pl9_nd5)


mod.pl9_nd6<-gls(contrast.wing.coverts.inv ~ Note.diversity + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(0.8, phy = Form.tree), data = Form.female, method = "ML") 
summary(mod.pl9_nd6)
AICc(mod.pl9_nd6)


mod.pl9_nd7<-gls(contrast.wing.coverts.inv ~ Note.diversity + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.4, phy = Form.tree), data = Form.female, method = "ML") 
summary(mod.pl9_nd7)
AICc(mod.pl9_nd7)


mod.pl9_nd8<-gls(contrast.wing.coverts.inv ~ Note.diversity + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.5, phy = Form.tree), data = Form.female, method = "ML")
summary(mod.pl9_nd8)
AICc(mod.pl9_nd8)



#Only note rate 
mod.pl9_nr5<-gls(contrast.wing.coverts.inv ~ Note.rate + Habitat_exposure + Foraging_strata, correlation  = corPagel(0.3, phy = Form.tree), data = Form.female, method = "ML") 
summary(mod.pl9_nr5)
AICc(mod.pl9_nr5)

mod.pl9_nr6<-gls(contrast.wing.coverts.inv ~ Note.rate + Habitat_exposure + Mixed.Flocking, correlation  = corPagel(0.4, phy = Form.tree), data = Form.female, method = "ML") 
summary(mod.pl9_nr6)
AICc(mod.pl9_nr6)

mod.pl9_nr7<-gls(contrast.wing.coverts.inv ~ Note.rate + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.4, phy = Form.tree), data = Form.female, method = "ML") 
summary(mod.pl9_nr7)
AICc(mod.pl9_nr7)

mod.pl9_nr8<-gls(contrast.wing.coverts.inv ~ Note.rate + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(0.4, phy = Form.tree), data = Form.female, method = "ML") 
summary(mod.pl9_nr8)
AICc(mod.pl9_nr8)

