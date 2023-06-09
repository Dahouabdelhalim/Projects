#PGLS plumage versus vocal traits in males, ecological and behavior interactive effects

#abrir df
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
#Models that had positive or negative associations to test interactions

#1. maxPower.dorsal	versus Delta.Time..s. (BM OU), Loudsong.mod.rate (BM PL OU), Note.diversity (BM PL OU), Peak.Freq..Hz. (BM)  
#BM model
#a)interaction Habitat

mod.bro1a<-gls(maxPower.dorsal ~ Delta.Time..s. * Habitat_exposure + Peak.Freq..Hz. * Habitat_exposure + Note.diversity * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro1a<-dredge(global.model=mod.bro1a) 
head(d.bro1a)
write.table(d.bro1a,"AICs_MPdorsalXvocaleco_BM_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro1a<-coefTable(d.bro1a) 

se.bro1a

se.bro1a.df<-do.call(rbind.data.frame,se.bro1a )
dim(se.bro1a.df)
str(se.bro1a.df)

write.table(se.bro1a.df,"SE_MPdorsalXvocaleco_BM_male_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bro1b<-gls(maxPower.dorsal ~ Delta.Time..s. * Foraging_strata + Peak.Freq..Hz. * Foraging_strata + Note.diversity * Foraging_strata + Loudsong.mod.rate * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro1b<-dredge(global.model=mod.bro1b) 
head(d.bro1b)
write.table(d.bro1b,"AICs_MPdorsalXvocaleco_BM_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro1b<-coefTable(d.bro1b) 

se.bro1b

se.bro1b.df<-do.call(rbind.data.frame,se.bro1b )
dim(se.bro1b.df)
str(se.bro1b.df)

write.table(se.bro1b.df,"SE_MPdorsalXvocaleco_BM_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.bro1c<-gls(maxPower.dorsal ~ Delta.Time..s. * Mixed.Flocking + Peak.Freq..Hz. * Mixed.Flocking + Note.diversity * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro1c<-dredge(global.model=mod.bro1c) 
head(d.bro1c)
write.table(d.bro1c,"AICs_MPdorsalXvocaleco_BM_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro1c<-coefTable(d.bro1c) 

se.bro1c

se.bro1c.df<-do.call(rbind.data.frame,se.bro1c)
dim(se.bro1c.df)
str(se.bro1c.df)

write.table(se.bro1c.df,"SE_MPdorsalXvocaleco_BM_male_msf.csv",sep = ",")

####################################################################################################################################################################################################################################################################################################################################################################
#OU model

#a)interaction Habitat
mod.ou1a<-gls(maxPower.dorsal ~ Delta.Time..s. * Habitat_exposure + Peak.Freq..Hz. * Habitat_exposure + Note.diversity * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou1a<-dredge(global.model=mod.ou1a) 
head(d.ou1a)
write.table(d.ou1a,"AICs_MPdorsalXvocaleco_OU_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou1a<-coefTable(d.ou1a) 

se.ou1a

se.ou1a.df<-do.call(rbind.data.frame,se.ou1a )
dim(se.ou1a.df)
str(se.ou1a.df)

write.table(se.ou1a.df,"SE_MPdorsalXvocaleco_OU_male_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.ou1b<-gls(maxPower.dorsal ~ Delta.Time..s. * Foraging_strata + Peak.Freq..Hz. * Foraging_strata + Note.diversity * Foraging_strata + Loudsong.mod.rate * Foraging_strata,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou1b<-dredge(global.model=mod.ou1b) 
head(d.ou1b)
write.table(d.ou1b,"AICs_MPdorsalXvocaleco_OU_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou1b<-coefTable(d.ou1b) 

se.ou1b

se.ou1b.df<-do.call(rbind.data.frame,se.ou1b )
dim(se.ou1b.df)
str(se.ou1b.df)

write.table(se.ou1b.df,"SE_MPdorsalXvocaleco_OU_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.ou1c<-gls(maxPower.dorsal ~ Delta.Time..s. * Mixed.Flocking + Peak.Freq..Hz. * Mixed.Flocking + Note.diversity * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou1c<-dredge(global.model=mod.ou1c) 
head(d.ou1c)
write.table(d.ou1c,"AICs_MPdorsalXvocaleco_OU_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou1c<-coefTable(d.ou1c) 

se.ou1c

se.ou1c.df<-do.call(rbind.data.frame,se.ou1c)
dim(se.ou1c.df)
str(se.ou1c.df)

write.table(se.ou1c.df,"SE_MPdorsalXvocaleco_OU_male_msf.csv",sep = ",")


####################################################################################################################################################################################################################################################################################################################################################################
#Pagel's Lambda model

#a)interaction Habitat
mod.pl1a<-gls(maxPower.dorsal ~ Delta.Time..s. * Habitat_exposure + Peak.Freq..Hz. * Habitat_exposure + Note.diversity * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl1a<-dredge(global.model=mod.pl1a) 
head(d.pl1a)
write.table(d.pl1a,"AICs_MPdorsalXvocaleco_PL_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl1a<-coefTable(d.pl1a) 

se.pl1a

se.pl1a.df<-do.call(rbind.data.frame,se.pl1a )
dim(se.pl1a.df)
str(se.pl1a.df)

write.table(se.pl1a.df,"SE_MPdorsalXvocaleco_PL_male_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.pl1b<-gls(maxPower.dorsal ~ Delta.Time..s. * Foraging_strata + Peak.Freq..Hz. * Foraging_strata + Note.diversity * Foraging_strata + Loudsong.mod.rate * Foraging_strata, correlation  = corPagel(0.9, phy = Form.tree), data = Form.male, method = "ML")

d.pl1b<-dredge(global.model=mod.pl1b) 
head(d.pl1b)
write.table(d.pl1b,"AICs_MPdorsalXvocaleco_PL_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl1b<-coefTable(d.pl1b) 

se.pl1b

se.pl1b.df<-do.call(rbind.data.frame,se.pl1b )
dim(se.pl1b.df)
str(se.pl1b.df)

write.table(se.pl1b.df,"SE_MPdorsalXvocaleco_PL_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.pl1c<-gls(maxPower.dorsal ~ Delta.Time..s. * Mixed.Flocking + Peak.Freq..Hz. * Mixed.Flocking + Note.diversity * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl1c<-dredge(global.model=mod.pl1c) 
head(d.pl1c)
write.table(d.pl1c,"AICs_MPdorsalXvocaleco_PL_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl1c<-coefTable(d.pl1c) 

se.pl1c

se.pl1c.df<-do.call(rbind.data.frame,se.pl1c)
dim(se.pl1c.df)
str(se.pl1c.df)

write.table(se.pl1c.df,"SE_MPdorsalXvocaleco_PL_male_msf.csv",sep = ",")


####
###########
##########################################################################################################################################################################################################################################################################################
##2.	LuminanceMean.dorsal versus Loudsong.mod.rate (PL OU)  + Song.bandwidth (BM OU) + 
#BM model
#a)interaction Habitat

mod.bro2a<-gls(LuminanceMean.dorsal ~ Song.bandwidth * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro2a<-dredge(global.model=mod.bro2a) 
head(d.bro2a)
write.table(d.bro2a,"AICs_LumdorsalXvocaleco_BM_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro2a<-coefTable(d.bro2a) 

se.bro2a

se.bro2a.df<-do.call(rbind.data.frame,se.bro2a )
dim(se.bro2a.df)
str(se.bro2a.df)

write.table(se.bro2a.df,"SE_LumdorsalXvocaleco_BM_male_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bro2b<-gls(LuminanceMean.dorsal ~ Song.bandwidth * Foraging_strata + Loudsong.mod.rate * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro2b<-dredge(global.model=mod.bro2b) 
head(d.bro2b)
write.table(d.bro2b,"AICs_LumdorsalXvocaleco_BM_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro2b<-coefTable(d.bro2b) 

se.bro2b

se.bro2b.df<-do.call(rbind.data.frame,se.bro2b )
dim(se.bro2b.df)
str(se.bro2b.df)

write.table(se.bro2b.df,"SE_LumdorsalXvocaleco_BM_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.bro2c<-gls(LuminanceMean.dorsal ~ Song.bandwidth * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro2c<-dredge(global.model=mod.bro2c) 
head(d.bro2c)
write.table(d.bro2c,"AICs_LumdorsalXvocaleco_BM_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro2c<-coefTable(d.bro2c) 

se.bro2c

se.bro2c.df<-do.call(rbind.data.frame,se.bro2c)
dim(se.bro2c.df)
str(se.bro2c.df)

write.table(se.bro2c.df,"SE_LumdorsalXvocaleco_BM_male_msf.csv",sep = ",")

####################################################################################################################################################################################################################################################################################################################################################################
#OU model

#a)interaction Habitat
mod.ou2a<-gls(LuminanceMean.dorsal ~ Song.bandwidth * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou2a<-dredge(global.model=mod.ou2a) 
head(d.ou2a)
write.table(d.ou2a,"AICs_LumdorsalXvocaleco_OU_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou2a<-coefTable(d.ou2a) 

se.ou2a

se.ou2a.df<-do.call(rbind.data.frame,se.ou2a )
dim(se.ou2a.df)
str(se.ou2a.df)

write.table(se.ou2a.df,"SE_LumdorsalXvocaleco_OU_male_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.ou2b<-gls(LuminanceMean.dorsal ~ Song.bandwidth * Foraging_strata + Loudsong.mod.rate * Foraging_strata,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou2b<-dredge(global.model=mod.ou2b) 
head(d.ou2b)
write.table(d.ou2b,"AICs_LumdorsalXvocaleco_OU_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou2b<-coefTable(d.ou2b) 

se.ou2b

se.ou2b.df<-do.call(rbind.data.frame,se.ou2b )
dim(se.ou2b.df)
str(se.ou2b.df)

write.table(se.ou2b.df,"SE_LumdorsalXvocaleco_OU_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.ou2c<-gls(LuminanceMean.dorsal ~ Song.bandwidth * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou2c<-dredge(global.model=mod.ou2c) 
head(d.ou2c)
write.table(d.ou2c,"AICs_LumdorsalXvocaleco_OU_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou2c<-coefTable(d.ou2c) 

se.ou2c

se.ou2c.df<-do.call(rbind.data.frame,se.ou2c)
dim(se.ou2c.df)
str(se.ou2c.df)

write.table(se.ou2c.df,"SE_LumdorsalXvocaleco_OU_male_msf.csv",sep = ",")


####################################################################################################################################################################################################################################################################################################################################################################
#Pagel's Lambda model

#a)interaction Habitat
mod.pl2a<-gls(LuminanceMean.dorsal ~ Song.bandwidth * Habitat_exposure + Loudsong.mod.rate * Habitat_exposure, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl2a<-dredge(global.model=mod.pl2a) 
head(d.pl2a)
write.table(d.pl2a,"AICs_LumdorsalXvocaleco_PL_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl2a<-coefTable(d.pl2a) 

se.pl2a

se.pl2a.df<-do.call(rbind.data.frame,se.pl2a )
dim(se.pl2a.df)
str(se.pl2a.df)

write.table(se.pl2a.df,"SE_LumdorsalXvocaleco_PL_male_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.pl2b<-gls(LuminanceMean.dorsal ~ Song.bandwidth * Foraging_strata + Loudsong.mod.rate * Foraging_strata, correlation  = corPagel(0.1, phy = Form.tree), data = Form.male, method = "ML")

d.pl2b<-dredge(global.model=mod.pl2b) 
head(d.pl2b)
write.table(d.pl2b,"AICs_LumdorsalXvocaleco_PL_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl2b<-coefTable(d.pl2b) 

se.pl2b

se.pl2b.df<-do.call(rbind.data.frame,se.pl2b )
dim(se.pl2b.df)
str(se.pl2b.df)

write.table(se.pl2b.df,"SE_LumdorsalXvocaleco_PL_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.pl2c<-gls(LuminanceMean.dorsal ~ Song.bandwidth * Mixed.Flocking + Loudsong.mod.rate * Mixed.Flocking, correlation  = corPagel(0.1, phy = Form.tree), data = Form.male, method = "ML")

d.pl2c<-dredge(global.model=mod.pl2c) 
head(d.pl2c)
write.table(d.pl2c,"AICs_LumdorsalXvocaleco_PL_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl2c<-coefTable(d.pl2c) 

se.pl2c

se.pl2c.df<-do.call(rbind.data.frame,se.pl2c)
dim(se.pl2c.df)
str(se.pl2c.df)

write.table(se.pl2c.df,"SE_LumdorsalXvocaleco_PL_male_msf.csv",sep = ",")


####
###########
##########################################################################################################################################################################################################################################################################################
##3. contrast.dorsal	versus Peak frequency (BM) 
#BM model
#a)interaction Habitat

mod.bro3a<-gls(contrast.dorsal.inv ~ Peak.Freq..Hz. * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro3a<-dredge(global.model=mod.bro3a) 
head(d.bro3a)
write.table(d.bro3a,"AICs_ContdorsalXvocaleco_BM_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro3a<-coefTable(d.bro3a) 

se.bro3a.df<-do.call(rbind.data.frame,se.bro3a )


write.table(se.bro3a.df,"SE_ContdorsalXvocaleco_BM_male_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bro3b<-gls(contrast.dorsal.inv ~ Peak.Freq..Hz. * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro3b<-dredge(global.model=mod.bro3b) 
head(d.bro3b)
write.table(d.bro3b,"AICs_ContdorsalXvocaleco_BM_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro3b<-coefTable(d.bro3b) 

se.bro3b.df<-do.call(rbind.data.frame,se.bro3b )


write.table(se.bro3b.df,"SE_ContdorsalXvocaleco_BM_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.bro3c<-gls(contrast.dorsal.inv ~ Peak.Freq..Hz. * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro3c<-dredge(global.model=mod.bro3c) 
head(d.bro3c)
write.table(d.bro3c,"AICs_ContdorsalXvocaleco_BM_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro3c<-coefTable(d.bro3c) 

se.bro3c.df<-do.call(rbind.data.frame,se.bro3c)


write.table(se.bro3c.df,"SE_ContdorsalXvocaleco_BM_male_msf.csv",sep = ",")




####
###########
##########################################################################################################################################################################################################################################################################################
##5. LuminanceMean.ventral versus Peak frequency (PL BM) 
#BM model
#a)interaction Habitat

mod.bro5a<-gls(LuminanceMean.ventral ~ Peak.Freq..Hz. * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro5a<-dredge(global.model=mod.bro5a) 
head(d.bro5a)
write.table(d.bro5a,"AICs_LumventralXvocaleco_BM_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro5a<-coefTable(d.bro5a) 

se.bro5a

se.bro5a.df<-do.call(rbind.data.frame,se.bro5a )
dim(se.bro5a.df)
str(se.bro5a.df)

write.table(se.bro5a.df,"SE_LumventralXvocaleco_BM_male_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bro5b<-gls(LuminanceMean.ventral ~ Peak.Freq..Hz. * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro5b<-dredge(global.model=mod.bro5b) 
head(d.bro5b)
write.table(d.bro5b,"AICs_LumventralXvocaleco_BM_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro5b<-coefTable(d.bro5b) 

se.bro5b

se.bro5b.df<-do.call(rbind.data.frame,se.bro5b )
dim(se.bro5b.df)
str(se.bro5b.df)

write.table(se.bro5b.df,"SE_LumventralXvocaleco_BM_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.bro5c<-gls(LuminanceMean.ventral ~ Peak.Freq..Hz. * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro5c<-dredge(global.model=mod.bro5c) 
head(d.bro5c)
write.table(d.bro5c,"AICs_LumventralXvocaleco_BM_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro5c<-coefTable(d.bro5c) 

se.bro5c

se.bro5c.df<-do.call(rbind.data.frame,se.bro5c)
dim(se.bro5c.df)
str(se.bro5c.df)

write.table(se.bro5c.df,"SE_LumventralXvocaleco_BM_male_msf.csv",sep = ",")


####################################################################################################################################################################################################################################################################################################################################################################
#Pagel's Lambda model

#a)interaction Habitat
mod.pl5a<-gls(LuminanceMean.ventral ~ Peak.Freq..Hz. * Habitat_exposure, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl5a<-dredge(global.model=mod.pl5a) 
head(d.pl5a)
write.table(d.pl5a,"AICs_LumventralXvocaleco_PL_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl5a<-coefTable(d.pl5a) 

se.pl5a

se.pl5a.df<-do.call(rbind.data.frame,se.pl5a )
dim(se.pl5a.df)
str(se.pl5a.df)

write.table(se.pl5a.df,"SE_LumventralXvocaleco_PL_male_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.pl5b<-gls(LuminanceMean.ventral ~ Peak.Freq..Hz. * Foraging_strata, correlation  = corPagel(0.9, phy = Form.tree), data = Form.male, method = "ML")

d.pl5b<-dredge(global.model=mod.pl5b) 
head(d.pl5b)
write.table(d.pl5b,"AICs_LumventralXvocaleco_PL_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl5b<-coefTable(d.pl5b) 

se.pl5b

se.pl5b.df<-do.call(rbind.data.frame,se.pl5b )
dim(se.pl5b.df)
str(se.pl5b.df)

write.table(se.pl5b.df,"SE_LumventralXvocaleco_PL_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.pl5c<-gls(LuminanceMean.ventral ~ Peak.Freq..Hz. * Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl5c<-dredge(global.model=mod.pl5c) 
head(d.pl5c)
write.table(d.pl5c,"AICs_LumventralXvocaleco_PL_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl5c<-coefTable(d.pl5c) 

se.pl5c

se.pl5c.df<-do.call(rbind.data.frame,se.pl5c)
dim(se.pl5c.df)
str(se.pl5c.df)

write.table(se.pl5c.df,"SE_LumventralXvocaleco_PL_male_msf.csv",sep = ",")


########
###################
####
###########
##########################################################################################################################################################################################################################################################################################
###6. contrast.ventral versus Note diversity (BM PL OU)

#BM model
#a)interaction Habitat

mod.bro6a<-gls(contrast.ventral.inv ~ Note.diversity * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro6a<-dredge(global.model=mod.bro6a) 
head(d.bro6a)
write.table(d.bro6a,"AICs_contventralXvocaleco_BM_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro6a<-coefTable(d.bro6a) 

se.bro6a

se.bro6a.df<-do.call(rbind.data.frame,se.bro6a )
dim(se.bro6a.df)
str(se.bro6a.df)

write.table(se.bro6a.df,"SE_contventralXvocaleco_BM_male_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bro6b<-gls(contrast.ventral.inv ~ Note.diversity * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro6b<-dredge(global.model=mod.bro6b) 
head(d.bro6b)
write.table(d.bro6b,"AICs_contventralXvocaleco_BM_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro6b<-coefTable(d.bro6b) 

se.bro6b

se.bro6b.df<-do.call(rbind.data.frame,se.bro6b )
dim(se.bro6b.df)
str(se.bro6b.df)

write.table(se.bro6b.df,"SE_contventralXvocaleco_BM_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.bro6c<-gls(contrast.ventral.inv ~ Note.diversity * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.male, method = "ML")

d.bro6c<-dredge(global.model=mod.bro6c) 
head(d.bro6c)
write.table(d.bro6c,"AICs_contventralXvocaleco_BM_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bro6c<-coefTable(d.bro6c) 

se.bro6c

se.bro6c.df<-do.call(rbind.data.frame,se.bro6c)
dim(se.bro6c.df)
str(se.bro6c.df)

write.table(se.bro6c.df,"SE_contventralXvocaleco_BM_male_msf.csv",sep = ",")

####################################################################################################################################################################################################################################################################################################################################################################
#OU model

#a)interaction Habitat
mod.ou6a<-gls(contrast.ventral.inv ~ Note.diversity * Habitat_exposure,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou6a<-dredge(global.model=mod.ou6a) 
head(d.ou6a)
write.table(d.ou6a,"AICs_contventralXvocaleco_OU_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou6a<-coefTable(d.ou6a) 

se.ou6a

se.ou6a.df<-do.call(rbind.data.frame,se.ou6a )
dim(se.ou6a.df)
str(se.ou6a.df)

write.table(se.ou6a.df,"SE_contventralXvocaleco_OU_male_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.ou6b<-gls(contrast.ventral.inv ~ Note.diversity * Foraging_strata,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou6b<-dredge(global.model=mod.ou6b) 
head(d.ou6b)
write.table(d.ou6b,"AICs_contventralXvocaleco_OU_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou6b<-coefTable(d.ou6b) 

se.ou6b

se.ou6b.df<-do.call(rbind.data.frame,se.ou6b )
dim(se.ou6b.df)
str(se.ou6b.df)

write.table(se.ou6b.df,"SE_contventralXvocaleco_OU_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.ou6c<-gls(contrast.ventral.inv ~ Note.diversity * Mixed.Flocking,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou6c<-dredge(global.model=mod.ou6c) 
head(d.ou6c)
write.table(d.ou6c,"AICs_contventralXvocaleco_OU_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou6c<-coefTable(d.ou6c) 

se.ou6c.df<-do.call(rbind.data.frame,se.ou6c)

write.table(se.ou6c.df,"SE_contventralXvocaleco_OU_male_msf.csv",sep = ",")


####################################################################################################################################################################################################################################################################################################################################################################
#Pagel's Lambda model

#a)interaction Habitat
mod.pl6a<-gls(contrast.ventral.inv ~ Note.diversity * Habitat_exposure, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl6a<-dredge(global.model=mod.pl6a) 
head(d.pl6a)
write.table(d.pl6a,"AICs_contventralXvocaleco_PL_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl6a<-coefTable(d.pl6a) 

se.pl6a.df<-do.call(rbind.data.frame,se.pl6a )


write.table(se.pl6a.df,"SE_contventralXvocaleco_PL_male_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.pl6b<-gls(contrast.ventral.inv ~ Note.diversity * Foraging_strata, correlation  = corPagel(0.9, phy = Form.tree), data = Form.male, method = "ML")

d.pl6b<-dredge(global.model=mod.pl6b) 
head(d.pl6b)
write.table(d.pl6b,"AICs_contventralXvocaleco_PL_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl6b<-coefTable(d.pl6b) 

se.pl6b

se.pl6b.df<-do.call(rbind.data.frame,se.pl6b )
dim(se.pl6b.df)
str(se.pl6b.df)

write.table(se.pl6b.df,"SE_contventralXvocaleco_PL_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.pl6c<-gls(contrast.ventral.inv ~ Note.diversity * Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl6c<-dredge(global.model=mod.pl6c) 
head(d.pl6c)
write.table(d.pl6c,"AICs_contventralXvocaleco_PL_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl6c<-coefTable(d.pl6c) 

se.pl6c

se.pl6c.df<-do.call(rbind.data.frame,se.pl6c)
dim(se.pl6c.df)
str(se.pl6c.df)

write.table(se.pl6c.df,"SE_contventralXvocaleco_PL_male_msf.csv",sep = ",")



#
#
#
##########################################################################################################################################################################################################################################################################################
###9. contrast.wing.coverts versus Delta.Time..s. (PL OU)

#OU model

#a)interaction Habitat
mod.ou9a<-gls(contrast.wing.coverts.inv ~ Delta.Time..s. * Habitat_exposure,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou9a<-dredge(global.model=mod.ou9a) 
head(d.ou9a)
write.table(d.ou9a,"AICs_contwingXvocaleco_OU_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou9a<-coefTable(d.ou9a) 

se.ou9a.df<-do.call(rbind.data.frame,se.ou9a )

write.table(se.ou9a.df,"SE_contwingXvocaleco_OU_male_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.ou9b<-gls(contrast.wing.coverts.inv ~ Delta.Time..s. * Foraging_strata,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou9b<-dredge(global.model=mod.ou9b) 
head(d.ou9b)
write.table(d.ou9b,"AICs_contwingXvocaleco_OU_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou9b<-coefTable(d.ou9b) 

se.ou9b.df<-do.call(rbind.data.frame,se.ou9b )

write.table(se.ou9b.df,"SE_contwingXvocaleco_OU_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.ou9c<-gls(contrast.wing.coverts.inv ~ Delta.Time..s. * Mixed.Flocking,correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.male, method = "ML")

d.ou9c<-dredge(global.model=mod.ou9c) 
head(d.ou9c)
write.table(d.ou9c,"AICs_contwingXvocaleco_OU_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou9c<-coefTable(d.ou9c) 

se.ou9c.df<-do.call(rbind.data.frame,se.ou9c)

write.table(se.ou9c.df,"SE_contwingXvocaleco_OU_male_msf.csv",sep = ",")


####################################################################################################################################################################################################################################################################################################################################################################
#Pagel's Lambda model

#a)interaction Habitat
mod.pl9a<-gls(contrast.wing.coverts.inv ~ Delta.Time..s. * Habitat_exposure, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl9a<-dredge(global.model=mod.pl9a) 
head(d.pl9a)
write.table(d.pl9a,"AICs_contwingXvocaleco_PL_male_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl9a<-coefTable(d.pl9a) 

se.pl9a.df<-do.call(rbind.data.frame,se.pl9a )


write.table(se.pl9a.df,"SE_contwingXvocaleco_PL_male_hab.csv",sep = ",")

##################################################################################################################################################################################
#b)interaction Foraging strata

mod.pl9b<-gls(contrast.wing.coverts.inv ~ Delta.Time..s. * Foraging_strata, correlation  = corPagel(0.9, phy = Form.tree), data = Form.male, method = "ML")

d.pl9b<-dredge(global.model=mod.pl9b) 
head(d.pl9b)
write.table(d.pl9b,"AICs_contwingXvocaleco_PL_male_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl9b<-coefTable(d.pl9b) 

se.pl9b.df<-do.call(rbind.data.frame,se.pl9b )


write.table(se.pl9b.df,"SE_contwingXvocaleco_PL_male_for.csv",sep = ",")


#################################################################################################################################################################################
#c)interaction MSF

mod.pl9c<-gls(contrast.wing.coverts.inv ~ Delta.Time..s. * Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.male, method = "ML")

d.pl9c<-dredge(global.model=mod.pl9c) 
head(d.pl9c)
write.table(d.pl9c,"AICs_contwingXvocaleco_PL_male_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl9c<-coefTable(d.pl9c) 

se.pl9c.df<-do.call(rbind.data.frame,se.pl9c)


write.table(se.pl9c.df,"SE_contwingXvocaleco_PL_male_msf.csv",sep = ",")

