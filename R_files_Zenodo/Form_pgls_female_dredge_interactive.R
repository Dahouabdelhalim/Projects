#PGLS plumage versus vocal traits in females, ecological and behavior interactive effects


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
#Models that had positive or negative associations to test interactions

#1. maxPower.dorsal	versus Loudsong.mod.rate (BM OU), Note.diversity (BM)  
#BM model
#a)interaction Habitat

mod.bro1a<-gls(maxPower.dorsal ~ Note.diversity * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro1a<-dredge(global.model=mod.bro1a)  
head(d.bro1a)
write.table(d.bro1a,"AICs_MPdorsalXvocaleco_BM_female_hab.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro1a<-coefTable(d.bro1a)  

se.bro1a

se.bro1a.df<-do.call(rbind.data.frame,se.bro1a )
dim(se.bro1a.df)
str(se.bro1a.df)

write.table(se.bro1a.df,"SE_MPdorsalXvocaleco_BM_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bro1b<-gls(maxPower.dorsal ~ Note.diversity * Foraging_strata  + Loudsong.mod.rate * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro1b<-dredge(global.model=mod.bro1b)  
head(d.bro1b)
write.table(d.bro1b,"AICs_MPdorsalXvocaleco_BM_female_for.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro1b<-coefTable(d.bro1b)  

se.bro1b

se.bro1b.df<-do.call(rbind.data.frame,se.bro1b )
dim(se.bro1b.df)
str(se.bro1b.df)

write.table(se.bro1b.df,"SE_MPdorsalXvocaleco_BM_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.bro1c<-gls(maxPower.dorsal ~ Note.diversity * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro1c<-dredge(global.model=mod.bro1c)  
head(d.bro1c)
write.table(d.bro1c,"AICs_MPdorsalXvocaleco_BM_female_msf.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro1c<-coefTable(d.bro1c)  

se.bro1c

se.bro1c.df<-do.call(rbind.data.frame,se.bro1c)
dim(se.bro1c.df)
str(se.bro1c.df)

write.table(se.bro1c.df,"SE_MPdorsalXvocaleco_BM_female_msf.csv",sep = ",")


####################################################################################################################################################################################################################################################################################################################################################################
#OU model

#a)interaction Habitat
mod.ou1a<-gls(maxPower.dorsal ~ Note.diversity * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou1a<-dredge(global.model=mod.ou1a)  
head(d.ou1a)
write.table(d.ou1a,"AICs_MPdorsalXvocaleco_OU_female_hab.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.ou1a<-coefTable(d.ou1a)  

se.ou1a

se.ou1a.df<-do.call(rbind.data.frame,se.ou1a )
dim(se.ou1a.df)
str(se.ou1a.df)

write.table(se.ou1a.df,"SE_MPdorsalXvocaleco_OU_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.ou1b<-gls(maxPower.dorsal ~ Note.diversity * Foraging_strata  + Loudsong.mod.rate * Foraging_strata, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou1b<-dredge(global.model=mod.ou1b)  
head(d.ou1b)
write.table(d.ou1b,"AICs_MPdorsalXvocaleco_OU_female_for.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.ou1b<-coefTable(d.ou1b)  

se.ou1b

se.ou1b.df<-do.call(rbind.data.frame,se.ou1b )
dim(se.ou1b.df)
str(se.ou1b.df)

write.table(se.ou1b.df,"SE_MPdorsalXvocaleco_OU_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.ou1c<-gls(maxPower.dorsal ~ Note.diversity * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou1c<-dredge(global.model=mod.ou1c)  
head(d.ou1c)
write.table(d.ou1c,"AICs_MPdorsalXvocaleco_OU_female_msf.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.ou1c<-coefTable(d.ou1c)  

se.ou1c

se.ou1c.df<-do.call(rbind.data.frame,se.ou1c)
dim(se.ou1c.df)
str(se.ou1c.df)

write.table(se.ou1c.df,"SE_MPdorsalXvocaleco_OU_female_msf.csv",sep = ",")

####
###########
##########################################################################################################################################################################################################################################################################################
##2.	LuminanceMean.dorsal versus Delta.Time..s. (BM PL OU) + Loudsong.mod.rate (BM PL OU)  + Note.count (BM PL OU) + Song.bandwidth (PL OU) 
#BM model
#a)interaction Habitat

mod.bro2a<-gls(LuminanceMean.dorsal ~ Note.count * Habitat_exposure + Delta.Time..s. * Habitat_exposure + Song.bandwidth * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro2a<-dredge(global.model=mod.bro2a)  
head(d.bro2a)
write.table(d.bro2a,"AICs_LumdorsalXvocaleco_BM_female_hab.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro2a<-coefTable(d.bro2a)  

se.bro2a

se.bro2a.df<-do.call(rbind.data.frame,se.bro2a )
dim(se.bro2a.df)
str(se.bro2a.df)

write.table(se.bro2a.df,"SE_LumdorsalXvocaleco_BM_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bro2b<-gls(LuminanceMean.dorsal ~ Note.count * Foraging_strata + Delta.Time..s. * Foraging_strata + Song.bandwidth * Foraging_strata + Loudsong.mod.rate * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro2b<-dredge(global.model=mod.bro2b)  
head(d.bro2b)
write.table(d.bro2b,"AICs_LumdorsalXvocaleco_BM_female_for.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro2b<-coefTable(d.bro2b)  

se.bro2b

se.bro2b.df<-do.call(rbind.data.frame,se.bro2b )
dim(se.bro2b.df)
str(se.bro2b.df)

write.table(se.bro2b.df,"SE_LumdorsalXvocaleco_BM_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.bro2c<-gls(LuminanceMean.dorsal ~ Note.count * Mixed.Flocking + Delta.Time..s. * Mixed.Flocking + Song.bandwidth * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro2c<-dredge(global.model=mod.bro2c)  
head(d.bro2c)
write.table(d.bro2c,"AICs_LumdorsalXvocaleco_BM_female_msf.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro2c<-coefTable(d.bro2c)  

se.bro2c

se.bro2c.df<-do.call(rbind.data.frame,se.bro2c)
dim(se.bro2c.df)
str(se.bro2c.df)

write.table(se.bro2c.df,"SE_LumdorsalXvocaleco_BM_female_msf.csv",sep = ",")

####################################################################################################################################################################################################################################################################################################################################################################
#OU model

#a)interaction Habitat
mod.ou2a<-gls(LuminanceMean.dorsal ~ Note.count * Habitat_exposure + Delta.Time..s. * Habitat_exposure + Song.bandwidth * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure
,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou2a<-dredge(global.model=mod.ou2a)  
head(d.ou2a)
write.table(d.ou2a,"AICs_LumdorsalXvocaleco_OU_female_hab.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.ou2a<-coefTable(d.ou2a)  

se.ou2a

se.ou2a.df<-do.call(rbind.data.frame,se.ou2a )
dim(se.ou2a.df)
str(se.ou2a.df)

write.table(se.ou2a.df,"SE_LumdorsalXvocaleco_OU_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.ou2b<-gls(LuminanceMean.dorsal ~ Note.count * Foraging_strata + Delta.Time..s. * Foraging_strata + Song.bandwidth * Foraging_strata + Loudsong.mod.rate * Foraging_strata
,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou2b<-dredge(global.model=mod.ou2b)  
head(d.ou2b)
write.table(d.ou2b,"AICs_LumdorsalXvocaleco_OU_female_for.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.ou2b<-coefTable(d.ou2b)  

se.ou2b

se.ou2b.df<-do.call(rbind.data.frame,se.ou2b )
dim(se.ou2b.df)
str(se.ou2b.df)

write.table(se.ou2b.df,"SE_LumdorsalXvocaleco_OU_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.ou2c<-gls(LuminanceMean.dorsal ~ Note.count * Mixed.Flocking + Delta.Time..s. * Mixed.Flocking + Song.bandwidth * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking
,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou2c<-dredge(global.model=mod.ou2c)  
head(d.ou2c)
write.table(d.ou2c,"AICs_LumdorsalXvocaleco_OU_female_msf.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.ou2c<-coefTable(d.ou2c)  

se.ou2c

se.ou2c.df<-do.call(rbind.data.frame,se.ou2c)
dim(se.ou2c.df)
str(se.ou2c.df)

write.table(se.ou2c.df,"SE_LumdorsalXvocaleco_OU_female_msf.csv",sep = ",")


####################################################################################################################################################################################################################################################################################################################################################################
#Pagel's Lambda model

#a)interaction Habitat
mod.pl2a<-gls(LuminanceMean.dorsal ~ Note.count * Habitat_exposure + Delta.Time..s. * Habitat_exposure + Song.bandwidth * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure
, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML")

d.pl2a<-dredge(global.model=mod.pl2a)  
head(d.pl2a)
write.table(d.pl2a,"AICs_LumdorsalXvocaleco_PL_female_hab.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.pl2a<-coefTable(d.pl2a)  

se.pl2a

se.pl2a.df<-do.call(rbind.data.frame,se.pl2a )
dim(se.pl2a.df)
str(se.pl2a.df)

write.table(se.pl2a.df,"SE_LumdorsalXvocaleco_PL_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata



#Run separately note count duration
mod.pl2bnc<-gls(LuminanceMean.dorsal ~ Note.count * Foraging_strata + Delta.Time..s. * Foraging_strata, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML")

d.pl2bnc<-dredge(global.model=mod.pl2bnc)

write.table(d.pl2bnc,"AICs_LumdorsalXvocaleco_PL_female_for_ncdur.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.pl2bnc<-coefTable(d.pl2bnc)  

se.pl2bnc.df<-do.call(rbind.data.frame,se.pl2bnc)

write.table(se.pl2bnc.df,"SE_LumdorsalXvocaleco_PL_female_for_ncdur.csv",sep = ",")

#Run separately Loudsong.mod.rate
mod.pl2bfs<-gls(LuminanceMean.dorsal ~ Loudsong.mod.rate * Foraging_strata, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML")

d.pl2bfs<-dredge(global.model=mod.pl2bfs)

write.table(d.pl2bfs,"AICs_LumdorsalXvocaleco_PL_female_for_fs.csv",sep = ",",row.names = F)  
#obtain SE (standard error)
se.pl2bfs<-coefTable(d.pl2bfs)  

se.pl2bfs.df<-do.call(rbind.data.frame,se.pl2bfs)

write.table(se.pl2bfs.df,"SE_LumdorsalXvocaleco_PL_female_for_fs.csv",sep = ",")

#run separately Song.bandwidth
mod.pl2bsb<-gls(LuminanceMean.dorsal ~ Song.bandwidth * Foraging_strata, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

d.pl2bsb<-dredge(global.model=mod.pl2bsb)

write.table(d.pl2bsb,"AICs_LumdorsalXvocaleco_PL_female_for_sb.csv",sep = ",",row.names = F)  
#obtain SE (standard error)
se.pl2bsb<-coefTable(d.pl2bsb)  

se.pl2bsb.df<-do.call(rbind.data.frame,se.pl2bsb)

write.table(se.pl2bsb.df,"SE_LumdorsalXvocaleco_PL_female_for_sb.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF


#Part 1
mod.pl2c1<-gls(LuminanceMean.dorsal ~ Note.count * Mixed.Flocking + Delta.Time..s. * Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

d.pl2c1<-dredge(global.model=mod.pl2c1) 

write.table(d.pl2c1,"AICs_LumdorsalXvocaleco_PL_female_msf_1.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl2c1<-coefTable(d.pl2c1) 

se.pl2c1.df<-do.call(rbind.data.frame,se.pl2c1)

write.table(se.pl2c1.df,"SE_LumdorsalXvocaleco_PL_female_msf_1.csv",sep = ",")

#Part 2
mod.pl2c2<-gls(LuminanceMean.dorsal ~ Song.bandwidth * Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.female, method = "ML")

d.pl2c2<-dredge(global.model=mod.pl2c2) 

write.table(d.pl2c2,"AICs_LumdorsalXvocaleco_PL_female_msf_2.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl2c2<-coefTable(d.pl2c2) 

se.pl2c2.df<-do.call(rbind.data.frame,se.pl2c2)

write.table(se.pl2c2.df,"SE_LumdorsalXvocaleco_PL_female_msf_2.csv",sep = ",")

#Part 3
mod.pl2c3<-gls(LuminanceMean.dorsal ~ Loudsong.mod.rate * Mixed.Flocking, correlation  = corPagel(0.9, phy = Form.tree), data = Form.female, method = "ML")

d.pl2c3<-dredge(global.model=mod.pl2c3) 

write.table(d.pl2c3,"AICs_LumdorsalXvocaleco_PL_female_msf_3.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl2c3<-coefTable(d.pl2c3) 

se.pl2c3.df<-do.call(rbind.data.frame,se.pl2c3)

write.table(se.pl2c3.df,"SE_LumdorsalXvocaleco_PL_female_msf_3.csv",sep = ",")


####
###########
##########################################################################################################################################################################################################################################################################################
##3. contrast.dorsal	versus Note.diversity (BM) + Loudsong.mod.rate (BM OU) 
#BM model
#a)interaction Habitat

mod.bro3a<-gls(contrast.dorsal.inv ~ Note.diversity * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro3a<-dredge(global.model=mod.bro3a) 

head(d.bro3a)
write.table(d.bro3a,"AICs_ContdorsalXvocaleco_BM_female_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro3a<-coefTable(d.bro3a) 

se.bro3a.df<-do.call(rbind.data.frame,se.bro3a )


write.table(se.bro3a.df,"SE_ContdorsalXvocaleco_BM_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bro3b<-gls(contrast.dorsal.inv ~ Note.diversity * Foraging_strata + Loudsong.mod.rate * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro3b<-dredge(global.model=mod.bro3b) 
head(d.bro3b)

write.table(d.bro3b,"AICs_ContdorsalXvocaleco_BM_female_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro3b<-coefTable(d.bro3b) 

se.bro3b.df<-do.call(rbind.data.frame,se.bro3b )


write.table(se.bro3b.df,"SE_ContdorsalXvocaleco_BM_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.bro3c<-gls(contrast.dorsal.inv ~ Note.diversity * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro3c<-dredge(global.model=mod.bro3c) 
head(d.bro3c)

write.table(d.bro3c,"AICs_ContdorsalXvocaleco_BM_female_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro3c<-coefTable(d.bro3c) 

se.bro3c.df<-do.call(rbind.data.frame,se.bro3c)


write.table(se.bro3c.df,"SE_ContdorsalXvocaleco_BM_female_msf.csv",sep = ",")


##################################
############################
#OU model

#a)interaction Habitat

mod.ou3a<-gls(contrast.dorsal.inv ~ Note.diversity * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou3a<-dredge(global.model=mod.ou3a) 
head(d.ou3a)

write.table(d.ou3a,"AICs_ContdorsalXvocaleco_OU_female_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou3a<-coefTable(d.ou3a) 

se.ou3a.df<-do.call(rbind.data.frame,se.ou3a )


write.table(se.ou3a.df,"SE_ContdorsalXvocaleco_OU_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.ou3b<-gls(contrast.dorsal.inv ~ Note.diversity * Foraging_strata + Loudsong.mod.rate * Foraging_strata, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou3b<-dredge(global.model=mod.ou3b) 
head(d.ou3b)

write.table(d.ou3b,"AICs_ContdorsalXvocaleco_OU_female_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou3b<-coefTable(d.ou3b) 

se.ou3b.df<-do.call(rbind.data.frame,se.ou3b )


write.table(se.ou3b.df,"SE_ContdorsalXvocaleco_OU_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.ou3c<-gls(contrast.dorsal.inv ~ Note.diversity * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou3c<-dredge(global.model=mod.ou3c) 
head(d.ou3c)

write.table(d.ou3c,"AICs_ContdorsalXvocaleco_OU_female_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou3c<-coefTable(d.ou3c) 

se.ou3c.df<-do.call(rbind.data.frame,se.ou3c)


write.table(se.ou3c.df,"SE_ContdorsalXvocaleco_OU_female_msf.csv",sep = ",")


####
###########
##########################################################################################################################################################################################################################################################################################
##4. maxPower.ventral versus Peak.Freq..Hz. (BM PL OU) + Song.bandwidth  (BM PL OU)
#BM model
#a)interaction Habitat
mod.bro4a<-gls(maxPower.ventral ~ Song.bandwidth * Habitat_exposure + Peak.Freq..Hz. * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro4a<-dredge(global.model=mod.bro4a) 
head(d.bro4a)

write.table(d.bro4a,"AICs_MPventralXvocaleco_BM_female_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro4a<-coefTable(d.bro4a) 

se.bro4a.df<-do.call(rbind.data.frame,se.bro4a )


write.table(se.bro4a.df,"SE_MPventralXvocaleco_BM_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bro4b<-gls(maxPower.ventral ~ Song.bandwidth * Foraging_strata + Peak.Freq..Hz. * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro4b<-dredge(global.model=mod.bro4b) 
head(d.bro4b)

write.table(d.bro4b,"AICs_MPventralXvocaleco_BM_female_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro4b<-coefTable(d.bro4b) 

se.bro4b.df<-do.call(rbind.data.frame,se.bro4b )


write.table(se.bro4b.df,"SE_MPventralXvocaleco_BM_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.bro4c<-gls(maxPower.ventral ~ Song.bandwidth * Mixed.Flocking + Peak.Freq..Hz. * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro4c<-dredge(global.model=mod.bro4c) 
head(d.bro4c)

write.table(d.bro4c,"AICs_MPventralXvocaleco_BM_female_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro4c<-coefTable(d.bro4c) 

se.bro4c

se.bro4c.df<-do.call(rbind.data.frame,se.bro4c)
dim(se.bro4c.df)
str(se.bro4c.df)

write.table(se.bro4c.df,"SE_MPventralXvocaleco_BM_female_msf.csv",sep = ",")

####################################################################################################################################################################################################################################################################################################################################################################
#OU model

#a)interaction Habitat
mod.ou4a<-gls(maxPower.ventral ~ Song.bandwidth * Habitat_exposure + Peak.Freq..Hz. * Habitat_exposure,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou4a<-dredge(global.model=mod.ou4a) 
head(d.ou4a)

write.table(d.ou4a,"AICs_MPventralXvocaleco_OU_female_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou4a<-coefTable(d.ou4a) 

se.ou4a

se.ou4a.df<-do.call(rbind.data.frame,se.ou4a )
dim(se.ou4a.df)
str(se.ou4a.df)

write.table(se.ou4a.df,"SE_MPventralXvocaleco_OU_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.ou4b<-gls(maxPower.ventral ~ Song.bandwidth * Foraging_strata + Peak.Freq..Hz. * Foraging_strata,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou4b<-dredge(global.model=mod.ou4b) 
head(d.ou4b)

write.table(d.ou4b,"AICs_MPventralXvocaleco_OU_female_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou4b<-coefTable(d.ou4b) 

se.ou4b

se.ou4b.df<-do.call(rbind.data.frame,se.ou4b )
dim(se.ou4b.df)
str(se.ou4b.df)

write.table(se.ou4b.df,"SE_MPventralXvocaleco_OU_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.ou4c<-gls(maxPower.ventral ~ Song.bandwidth * Mixed.Flocking + Peak.Freq..Hz. * Mixed.Flocking,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou4c<-dredge(global.model=mod.ou4c) 
head(d.ou4c)
write.table(d.ou4c,"AICs_MPventralXvocaleco_OU_female_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou4c<-coefTable(d.ou4c) 

se.ou4c

se.ou4c.df<-do.call(rbind.data.frame,se.ou4c)
dim(se.ou4c.df)
str(se.ou4c.df)

write.table(se.ou4c.df,"SE_MPventralXvocaleco_OU_female_msf.csv",sep = ",")


####################################################################################################################################################################################################################################################################################################################################################################
#Pagel's Lambda model

#a)interaction Habitat
mod.pl4a<-gls(maxPower.ventral ~ Song.bandwidth * Habitat_exposure + Peak.Freq..Hz. * Habitat_exposure, correlation  = corPagel(0, phy = Form.tree), data = Form.female, method = "ML")

d.pl4a<-dredge(global.model=mod.pl4a)  
head(d.pl4a)
write.table(d.pl4a,"AICs_MPventralXvocaleco_PL_female_hab.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.pl4a<-coefTable(d.pl4a)  

se.pl4a

se.pl4a.df<-do.call(rbind.data.frame,se.pl4a )
dim(se.pl4a.df)
str(se.pl4a.df)

write.table(se.pl4a.df,"SE_MPventralXvocaleco_PL_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.pl4b<-gls(maxPower.ventral ~ Song.bandwidth * Foraging_strata + Peak.Freq..Hz. * Foraging_strata, correlation  = corPagel(0.6, phy = Form.tree), data = Form.female, method = "ML")

d.pl4b<-dredge(global.model=mod.pl4b)  
head(d.pl4b)
write.table(d.pl4b,"AICs_MPventralXvocaleco_PL_female_for.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.pl4b<-coefTable(d.pl4b)  

se.pl4b

se.pl4b.df<-do.call(rbind.data.frame,se.pl4b )
dim(se.pl4b.df)
str(se.pl4b.df)

write.table(se.pl4b.df,"SE_MPventralXvocaleco_PL_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.pl4c<-gls(maxPower.ventral ~ Song.bandwidth * Mixed.Flocking + Peak.Freq..Hz. * Mixed.Flocking, correlation  = corPagel(0, phy = Form.tree), data = Form.female, method = "ML")

d.pl4c<-dredge(global.model=mod.pl4c)  
head(d.pl4c)
write.table(d.pl4c,"AICs_MPventralXvocaleco_PL_female_msf.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.pl4c<-coefTable(d.pl4c)  

se.pl4c

se.pl4c.df<-do.call(rbind.data.frame,se.pl4c)
dim(se.pl4c.df)
str(se.pl4c.df)

write.table(se.pl4c.df,"SE_MPventralXvocaleco_PL_female_msf.csv",sep = ",")


####
###########
##########################################################################################################################################################################################################################################################################################
##5. LuminanceMean.ventral versus Peak.Freq..Hz. (BM) 
#BM model
#a)interaction Habitat

mod.bro5a<-gls(LuminanceMean.ventral ~ Peak.Freq..Hz. * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro5a<-dredge(global.model=mod.bro5a)  
head(d.bro5a)
write.table(d.bro5a,"AICs_LumventralXvocaleco_BM_female_hab.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro5a<-coefTable(d.bro5a)  

se.bro5a

se.bro5a.df<-do.call(rbind.data.frame,se.bro5a )
dim(se.bro5a.df)
str(se.bro5a.df)

write.table(se.bro5a.df,"SE_LumventralXvocaleco_BM_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bro5b<-gls(LuminanceMean.ventral ~ Peak.Freq..Hz. * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro5b<-dredge(global.model=mod.bro5b)  
head(d.bro5b)
write.table(d.bro5b,"AICs_LumventralXvocaleco_BM_female_for.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro5b<-coefTable(d.bro5b)  

se.bro5b

se.bro5b.df<-do.call(rbind.data.frame,se.bro5b )
dim(se.bro5b.df)
str(se.bro5b.df)

write.table(se.bro5b.df,"SE_LumventralXvocaleco_BM_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.bro5c<-gls(LuminanceMean.ventral ~ Peak.Freq..Hz. * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro5c<-dredge(global.model=mod.bro5c)  
head(d.bro5c)
write.table(d.bro5c,"AICs_LumventralXvocaleco_BM_female_msf.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro5c<-coefTable(d.bro5c)  

se.bro5c

se.bro5c.df<-do.call(rbind.data.frame,se.bro5c)
dim(se.bro5c.df)
str(se.bro5c.df)

write.table(se.bro5c.df,"SE_LumventralXvocaleco_BM_female_msf.csv",sep = ",")



########
###################
####
###########
##########################################################################################################################################################################################################################################################################################
###6. contrast.ventral versus Note.count (BM OU) + Note diversity (BM OU) + Loudsong.mod.rate (BM OU)

#BM model
#a)interaction Habitat

mod.bro6a<-gls(contrast.ventral.inv ~ Note.count * Habitat_exposure + Note.diversity * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro6a<-dredge(global.model=mod.bro6a)  
head(d.bro6a)
write.table(d.bro6a,"AICs_contventralXvocaleco_BM_female_hab.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro6a<-coefTable(d.bro6a)  

se.bro6a

se.bro6a.df<-do.call(rbind.data.frame,se.bro6a )
dim(se.bro6a.df)
str(se.bro6a.df)

write.table(se.bro6a.df,"SE_contventralXvocaleco_BM_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bro6b<-gls(contrast.ventral.inv ~ Note.count * Foraging_strata + Note.diversity * Foraging_strata + Loudsong.mod.rate * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro6b<-dredge(global.model=mod.bro6b)  
head(d.bro6b)
write.table(d.bro6b,"AICs_contventralXvocaleco_BM_female_for.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro6b<-coefTable(d.bro6b)  

se.bro6b

se.bro6b.df<-do.call(rbind.data.frame,se.bro6b )
dim(se.bro6b.df)
str(se.bro6b.df)

write.table(se.bro6b.df,"SE_contventralXvocaleco_BM_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.bro6c<-gls(contrast.ventral.inv ~ Note.count * Mixed.Flocking  + Note.diversity * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro6c<-dredge(global.model=mod.bro6c)  
head(d.bro6c)
write.table(d.bro6c,"AICs_contventralXvocaleco_BM_female_msf.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro6c<-coefTable(d.bro6c)  

se.bro6c

se.bro6c.df<-do.call(rbind.data.frame,se.bro6c)
dim(se.bro6c.df)
str(se.bro6c.df)

write.table(se.bro6c.df,"SE_contventralXvocaleco_BM_female_msf.csv",sep = ",")

####################################################################################################################################################################################################################################################################################################################################################################
#OU model

#a)interaction Habitat
mod.ou6a<-gls(contrast.ventral.inv ~ Note.count * Habitat_exposure + Note.diversity * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou6a<-dredge(global.model=mod.ou6a)  
head(d.ou6a)
write.table(d.ou6a,"AICs_contventralXvocaleco_OU_female_hab.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.ou6a<-coefTable(d.ou6a)  

se.ou6a

se.ou6a.df<-do.call(rbind.data.frame,se.ou6a )
dim(se.ou6a.df)
str(se.ou6a.df)

write.table(se.ou6a.df,"SE_contventralXvocaleco_OU_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.ou6b<-gls(contrast.ventral.inv ~ Note.count * Foraging_strata + Note.diversity * Foraging_strata + Loudsong.mod.rate * Foraging_strata, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou6b<-dredge(global.model=mod.ou6b)  
head(d.ou6b)
write.table(d.ou6b,"AICs_contventralXvocaleco_OU_female_for.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.ou6b<-coefTable(d.ou6b)  

se.ou6b

se.ou6b.df<-do.call(rbind.data.frame,se.ou6b )
dim(se.ou6b.df)
str(se.ou6b.df)

write.table(se.ou6b.df,"SE_contventralXvocaleco_OU_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.ou6c<-gls(contrast.ventral.inv ~ Note.count * Mixed.Flocking  + Note.diversity * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou6c<-dredge(global.model=mod.ou6c)  
head(d.ou6c)
write.table(d.ou6c,"AICs_contventralXvocaleco_OU_female_msf.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.ou6c<-coefTable(d.ou6c)  

se.ou6c.df<-do.call(rbind.data.frame,se.ou6c)

write.table(se.ou6c.df,"SE_contventralXvocaleco_OU_female_msf.csv",sep = ",")

#
#
#
##########################################################################################################################################################################################################################################################################################
###7. maxPower.wing.coverts	versus Peak.Freq..Hz.(BM) 
#BM model
#a)interaction Habitat

mod.bro7a<-gls(maxPower.wing.coverts ~ Peak.Freq..Hz. * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro7a<-dredge(global.model=mod.bro7a)  
head(d.bro7a)
write.table(d.bro7a,"AICs_MPwingXvocaleco_BM_female_hab.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro7a<-coefTable(d.bro7a)  

se.bro7a

se.bro7a.df<-do.call(rbind.data.frame,se.bro7a )
dim(se.bro7a.df)
str(se.bro7a.df)

write.table(se.bro7a.df,"SE_MPwingXvocaleco_BM_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bro7b<-gls(maxPower.wing.coverts ~ Peak.Freq..Hz. * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro7b<-dredge(global.model=mod.bro7b)  
head(d.bro7b)
write.table(d.bro7b,"AICs_MPwingXvocaleco_BM_female_for.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro7b<-coefTable(d.bro7b)  

se.bro7b

se.bro7b.df<-do.call(rbind.data.frame,se.bro7b )
dim(se.bro7b.df)
str(se.bro7b.df)

write.table(se.bro7b.df,"SE_MPwingXvocaleco_BM_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.bro7c<-gls(maxPower.wing.coverts ~ Peak.Freq..Hz. * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro7c<-dredge(global.model=mod.bro7c)  
head(d.bro7c)
write.table(d.bro7c,"AICs_MPwingXvocaleco_BM_female_msf.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro7c<-coefTable(d.bro7c)  

se.bro7c

se.bro7c.df<-do.call(rbind.data.frame,se.bro7c)
dim(se.bro7c.df)
str(se.bro7c.df)

write.table(se.bro7c.df,"SE_MPwingXvocaleco_BM_female_msf.csv",sep = ",")


#
#
#
##########################################################################################################################################################################################################################################################################################
###9. contrast.wing.coverts versus Note.diversity (BM PL OU)
#BM model
#a)interaction Habitat
mod.bro9a<-gls(contrast.wing.coverts.inv ~ Note.diversity * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro9a<-dredge(global.model=mod.bro9a)  
head(d.bro9a)
write.table(d.bro9a,"AICs_contwingXvocaleco_BM_female_hab.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro9a<-coefTable(d.bro9a)  

se.bro9a.df<-do.call(rbind.data.frame,se.bro9a )

write.table(se.bro9a.df,"SE_contwingXvocaleco_BM_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bro9b<-gls(contrast.wing.coverts.inv ~ Note.diversity * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro9b<-dredge(global.model=mod.bro9b)  
head(d.bro9b)
write.table(d.bro9b,"AICs_contwingXvocaleco_BM_female_for.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro9b<-coefTable(d.bro9b)  

se.bro9b.df<-do.call(rbind.data.frame,se.bro9b )

write.table(se.bro9b.df,"SE_contwingXvocaleco_BM_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.bro9c<-gls(contrast.wing.coverts.inv ~ Note.diversity * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.female, method = "ML")

d.bro9c<-dredge(global.model=mod.bro9c)  
head(d.bro9c)
write.table(d.bro9c,"AICs_contwingXvocaleco_BM_female_msf.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.bro9c<-coefTable(d.bro9c)  

se.bro9c.df<-do.call(rbind.data.frame,se.bro9c)

write.table(se.bro9c.df,"SE_contwingXvocaleco_BM_female_msf.csv",sep = ",")


#OU model
#a)interaction Habitat
mod.ou9a<-gls(contrast.wing.coverts.inv ~ Note.diversity * Habitat_exposure,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou9a<-dredge(global.model=mod.ou9a)  
head(d.ou9a)
write.table(d.ou9a,"AICs_contwingXvocaleco_OU_female_hab.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.ou9a<-coefTable(d.ou9a)  

se.ou9a.df<-do.call(rbind.data.frame,se.ou9a )

write.table(se.ou9a.df,"SE_contwingXvocaleco_OU_female_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.ou9b<-gls(contrast.wing.coverts.inv ~ Note.diversity * Foraging_strata,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou9b<-dredge(global.model=mod.ou9b)  
head(d.ou9b)
write.table(d.ou9b,"AICs_contwingXvocaleco_OU_female_for.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.ou9b<-coefTable(d.ou9b)  

se.ou9b.df<-do.call(rbind.data.frame,se.ou9b )

write.table(se.ou9b.df,"SE_contwingXvocaleco_OU_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.ou9c<-gls(contrast.wing.coverts.inv ~ Note.diversity * Mixed.Flocking,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.female, method = "ML")

d.ou9c<-dredge(global.model=mod.ou9c)  
head(d.ou9c)
write.table(d.ou9c,"AICs_contwingXvocaleco_OU_female_msf.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.ou9c<-coefTable(d.ou9c)  

se.ou9c.df<-do.call(rbind.data.frame,se.ou9c)

write.table(se.ou9c.df,"SE_contwingXvocaleco_OU_female_msf.csv",sep = ",")


####################################################################################################################################################################################################################################################################################################################################################################
#Pagel's Lambda model

#a)interaction Habitat
mod.pl9a<-gls(contrast.wing.coverts.inv ~ Note.diversity * Habitat_exposure, correlation  = corPagel(0.6, phy = Form.tree), data = Form.female, method = "ML")

d.pl9a<-dredge(global.model=mod.pl9a)  
head(d.pl9a)
write.table(d.pl9a,"AICs_contwingXvocaleco_PL_female_hab.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.pl9a<-coefTable(d.pl9a)  

se.pl9a.df<-do.call(rbind.data.frame,se.pl9a )


write.table(se.pl9a.df,"SE_contwingXvocaleco_PL_female_hab.csv",sep = ",")

##################################################################################################################################################################################
#b)interaction Foraging strata

mod.pl9b<-gls(contrast.wing.coverts.inv ~ Note.diversity * Foraging_strata, correlation  = corPagel(0.6, phy = Form.tree), data = Form.female, method = "ML")

d.pl9b<-dredge(global.model=mod.pl9b)  
head(d.pl9b)
write.table(d.pl9b,"AICs_contwingXvocaleco_PL_female_for.csv",sep = ",",row.names = F)  

#obtain SE (standard error)
se.pl9b<-coefTable(d.pl9b)  

se.pl9b.df<-do.call(rbind.data.frame,se.pl9b )


write.table(se.pl9b.df,"SE_contwingXvocaleco_PL_female_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF


#run separately
#Only note diversity 
#1
mod.pl9c.nd1<-gls(contrast.wing.coverts.inv ~ Note.diversity, correlation  = corPagel(0, phy = Form.tree), data = Form.female, method = "ML") 

summary(mod.pl9c.nd1)
AICc(mod.pl9c.nd1)

#2
mod.pl9c.nd2<-gls(contrast.wing.coverts.inv ~ Note.diversity * Mixed.Flocking, correlation  = corPagel(0.7, phy = Form.tree), data = Form.female, method = "ML")

summary(mod.pl9c.nd2)
AICc(mod.pl9c.nd2)

