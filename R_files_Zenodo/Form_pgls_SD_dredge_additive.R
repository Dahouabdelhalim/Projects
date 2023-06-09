#PGLS sexual dimorphism, ecological and behavior additive effects

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
#Models

#1. ED_plum	versus ED_vocal (BM model)
mod.bro1<-gls(ED_plum ~ ED_vocal + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corBrownian(1, phy = Form.tree), data = Form.SDeco, method = "ML")

d.bro1<-dredge(global.model=mod.bro1) #run all associations at once
head(d.bro1)
write.table(d.bro1,"AICs_SD_plumXSD_vocal_BM.csv",sep = ",",row.names = F) 

se.bro1<-coefTable(d.bro1) #obtain SE (standard error) off all tested associations

se.bro1

se.bro1.df<-do.call(rbind.data.frame,se.bro1 )
dim(se.bro1.df)
str(se.bro1.df)

write.table(se.bro1.df,"SE_SD_plumXSD_vocaleco_BM.csv",sep = ",")

  

########################################################################################################################################################################################
######
#1. ED_plum	versus ED_vocal (OU model)
mod.ou1<-gls(ED_plum ~ ED_vocal + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corMartins(1, phy = Form.tree, fixed = F), data = Form.SDeco, method = "ML")

d.ou1<-dredge(global.model=mod.ou1)

write.table(d.ou1,"AICs_SD_plumXSD_vocal_OU.csv",sep = ",",row.names = F)

#obtain SE (standard error)
se.ou1<-coefTable(d.ou1) 
head(se.ou1)
se.ou1.df<-do.call(rbind.data.frame,se.ou1 )
dim(se.ou1.df)
str(se.ou1.df)

write.table(se.ou1.df,"SE_SD_plumXSD_vocaleco_OU.csv",sep = ",")


####################################################################################################################################################################################################
##################################################################################################

#1. ED_plum	versus ED_vocal (Pagel's lambda)

mod.pag<-gls(ED_plum ~ ED_vocal + Habitat_exposure + Foraging_strata + Mixed.Flocking, correlation  = corPagel(1, phy = Form.tree), data = Form.SDeco, method = "ML")

d.pag<-dredge(global.model=mod.pag) 
head(d.pag)
write.table(d.pag,"AICs_SD_plumXSD_vocal_pag.csv",sep = ",",row.names = F) 

#obtain SE (standard error)
se.pag<-coefTable(d.pag) 

se.pag

se.pag.df<-do.call(rbind.data.frame,se.pag )
dim(se.pag.df)
str(se.pag.df)

write.table(se.pag.df,"SE_SD_plumXSD_vocaleco_pag.csv",sep = ",")



#############
#Summary statistics of the best model

mod.bm1<-gls(ED_plum ~ ED_vocal, correlation  = corBrownian(1, phy = Form.tree), data = Form.SDeco, method = "ML")

summary(mod.bm1)

