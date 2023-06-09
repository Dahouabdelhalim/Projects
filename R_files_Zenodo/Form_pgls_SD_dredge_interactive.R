#PGLS sexual dimorphism, ecological and behavior interactive effects

Form.SDeco<-read.csv("Form_SD_eco_9traits.csv",row.names = 1)
head(Form.SDeco)

library(tidyr)
#remove NAs
Form.SDeco<-drop_na(Form.SDeco)
Form.SDeco

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
obj<-name.check(Form.tree,Form.SDeco)
obj

#Remove tips that we don't have data
Form.tree<- drop.tip(Form.tree, obj$tree_not_data)
plot(Form.tree)

#Check tree and data again
name.check(Form.tree, Form.SDeco)

#Branch length transformation. Necessary to avoid error you were getting before. Does not affect results. 
Form.tree$edge.length <- Form.tree$edge.length * 100

####################################################################################################################################################################################################################################################################################
#1. ED_plum versus ED_vocal (BM model)

#a)interaction Habitat

mod.bm.a<-gls(ED_plum ~ ED_vocal * Habitat_exposure, correlation  = corBrownian(1, phy = Form.tree), data = Form.SDeco, method = "ML")

d.bm.a<-dredge(global.model=mod.bm.a) 
head(d.bm.a)

write.table(d.bm.a,"AICs_SD_plumXSD_vocal_BM_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bm.a<-coefTable(d.bm.a) 

se.bm.a

se.bm.a.df<-do.call(rbind.data.frame,se.bm.a )
dim(se.bm.a.df)
str(se.bm.a.df)

write.table(se.bm.a.df,"SE_SD_plumXSD_vocaleco_BM_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.bm.b<-gls(ED_plum ~ ED_vocal * Foraging_strata, correlation  = corBrownian(1, phy = Form.tree), data = Form.SDeco, method = "ML")

d.bm.b<-dredge(global.model=mod.bm.b) 

head(d.bm.b)
write.table(d.bm.b,"AICs_SD_plumXSD_vocal_BM_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bm.b<-coefTable(d.bm.b) 

se.bm.b

se.bm.b.df<-do.call(rbind.data.frame,se.bm.b)
dim(se.bm.b.df)
str(se.bm.b.df)

write.table(se.bm.b.df,"SE_SD_plumXSD_vocaleco_BM_for.csv",sep = ",")

##################################################################################################################################################################################
#c)interaction MSF

mod.bm.c<-gls(ED_plum ~ ED_vocal * Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.SDeco, method = "ML")

d.bm.c<-dredge(global.model=mod.bm.c) 
head(d.bm.c)

write.table(d.bm.c,"AICs_SD_plumXSD_vocal_BM_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.bm.c<-coefTable(d.bm.c) 

se.bm.c

se.bm.c.df<-do.call(rbind.data.frame,se.bm.c)
dim(se.bm.c.df)
str(se.bm.c.df)

write.table(se.bm.c.df,"SE_SD_plumXSD_vocaleco_BM_msf.csv",sep = ",")


####################################################################################################################################################################################################################################################################################################################################################################
#1. ED_plum versus ED_vocal (OU model)

#a)interaction Habitat

mod.ou.a<-gls(ED_plum ~ ED_vocal * Habitat_exposure, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.SDeco, method = "ML")

d.ou.a<-dredge(global.model=mod.ou.a) 
head(d.ou.a)

write.table(d.ou.a,"AICs_SD_plumXSD_vocal_OU_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou.a<-coefTable(d.ou.a) 

se.ou.a

se.ou.a.df<-do.call(rbind.data.frame,se.ou.a )
dim(se.ou.a.df)
str(se.ou.a.df)

write.table(se.ou.a.df,"SE_SD_plumXSD_vocaleco_OU_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.ou.b<-gls(ED_plum ~ ED_vocal * Foraging_strata, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.SDeco, method = "ML")
summary(mod.ou.b)

d.ou.b<-dredge(global.model=mod.ou.b) 
head(d.ou.b)

write.table(d.ou.b,"AICs_SD_plumXSD_vocal_OU_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou.b<-coefTable(d.ou.b) 
se.ou.b

se.ou.b.df<-do.call(rbind.data.frame,se.ou.b)
dim(se.ou.b.df)
str(se.ou.b.df)

write.table(se.ou.b.df,"SE_SD_plumXSD_vocaleco_OU_for.csv",sep = ",")

##################################################################################################################################################################################
#c)interaction MSF

mod.ou.c<-gls(ED_plum ~ ED_vocal * Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.SDeco, method = "ML")

d.ou.c<-dredge(global.model=mod.ou.c) 
head(d.ou.c)

write.table(d.ou.c,"AICs_SD_plumXSD_vocal_OU_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.ou.c<-coefTable(d.ou.c) 
se.ou.c

se.ou.c.df<-do.call(rbind.data.frame,se.ou.c)
dim(se.ou.c.df)
str(se.ou.c.df)

write.table(se.ou.c.df,"SE_SD_plumXSD_vocaleco_OU_msf.csv",sep = ",")

####################################################################################################################################################################################################################################################################################################################################################################
#1. ED_plum versus ED_vocal (Pagel's lambda)

#a)interaction Habitat

mod.pl.a<-gls(ED_plum ~ ED_vocal * Habitat_exposure, correlation  = corPagel(1, phy = Form.tree), data = Form.SDeco, method = "ML")

d.pl.a<-dredge(global.model=mod.pl.a) 
head(d.pl.a)

write.table(d.pl.a,"AICs_SD_plumXSD_vocal_PL_hab.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl.a<-coefTable(d.pl.a) 
se.pl.a

se.pl.a.df<-do.call(rbind.data.frame,se.pl.a )
dim(se.pl.a.df)
str(se.pl.a.df)

write.table(se.pl.a.df,"SE_SD_plumXSD_vocaleco_PL_hab.csv",sep = ",")


##################################################################################################################################################################################
#b)interaction Foraging strata

mod.pl.b<-gls(ED_plum ~ ED_vocal * Foraging_strata, correlation  = corPagel(1, phy = Form.tree), data = Form.SDeco, method = "ML")

d.pl.b<-dredge(global.model=mod.pl.b) 
head(d.pl.b)

write.table(d.pl.b,"AICs_SD_plumXSD_vocal_PL_for.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl.b<-coefTable(d.pl.b) 

se.pl.b

se.pl.b.df<-do.call(rbind.data.frame,se.pl.b)
dim(se.pl.b.df)
str(se.pl.b.df)

write.table(se.pl.b.df,"SE_SD_plumXSD_vocaleco_PL_for.csv",sep = ",")

##################################################################################################################################################################################
#c)interaction MSF

mod.pl.c<-gls(ED_plum ~ ED_vocal * Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.SDeco, method = "ML")

d.pl.c<-dredge(global.model=mod.pl.c) 
head(d.pl.c)

write.table(d.pl.c,"AICs_SD_plumXSD_vocal_PL_msf.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pl.c<-coefTable(d.pl.c) 

se.pl.c

se.pl.c.df<-do.call(rbind.data.frame,se.pl.c)
dim(se.pl.c.df)
str(se.pl.c.df)

write.table(se.pl.c.df,"SE_SD_plumXSD_vocaleco_PL_msf.csv",sep = ",")


###########################################################################################

