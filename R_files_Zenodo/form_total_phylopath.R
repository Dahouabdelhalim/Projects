
#open df
Form<-read.csv("Form_total_ok.csv",row.names = 1)
head(Form)

#Contrast values were inverted so that small values now represent lower contrast and higher  values represent higher contrast. (See Analyses of plumage data section for more info)
#Male contrast values
Form$CT.dor_m.iv<-Form$CT.dor_m*-1 #inversion of dorsal contrast
Form$CT.vent_m.iv<-Form$CT.vent_m*-1 #inversion of ventral contrast
Form$CT.win_m.iv<-Form$CT.win_m*-1 #inversion of wing contrast
#Female contrast values
Form$CT.dor_f.iv<-Form$CT.dor_f*-1 #inversion of dorsal contrast
Form$CT.vent_f.iv<-Form$CT.vent_f*-1 #inversion of ventral contrast
Form$CT.win_f.iv<-Form$CT.win_f*-1 #inversion of wing contrast


head(Form)


library(tidyr)
#remove NAs
Form<-drop_na(Form)
Form

#load packages
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(ggplot2)
library(dplyr)

#Open tree
Form.tree <- read.tree("formicivorini_UCEs_SAAC.tre")
plot(Form.tree)
Form.tree
#Check tree and data
obj<-name.check(Form.tree,Form)
obj
#Remove tips that we don't have data
Form.tree<- drop.tip(Form.tree, obj$tree_not_data)
plot(Form.tree)

name.check(Form.tree, Form)

library(phylopath)


#Models tested
m1.m<-define_model_set(null=c(),
                       palla=c(ME.dor_m ~ HE,Lum.dor_m ~ HE, CT.dor_m.iv ~ HE, ME.vent_m ~ HE, Lum.vent_m ~ HE, CT.vent_m.iv ~ HE, ME.win_m ~ HE, Lum.win_m ~ HE, CT.win_m.iv ~ HE, 
                               ME.dor_f ~ HE, Lum.dor_f ~ HE, CT.dor_f.iv ~ HE, ME.vent_f ~ HE, Lum.vent_f ~ HE, CT.vent_f.iv ~ HE, ME.win_f ~ HE, Lum.win_f ~ HE, CT.win_f.iv ~ HE, 
                               plum.SD ~ HE),
                       
                       pallb=c(ME.dor_m ~ FE, Lum.dor_m ~ FE, CT.dor_m.iv ~ FE, ME.vent_m ~ FE, Lum.vent_m ~ FE, CT.vent_m.iv ~ FE, ME.win_m ~ FE, Lum.win_m ~ FE, CT.win_m.iv ~ FE, 
                               ME.dor_f ~ FE, Lum.dor_f ~ FE, CT.dor_f.iv ~ FE, ME.vent_f ~ FE, Lum.vent_f ~ FE, CT.vent_f.iv ~ FE, ME.win_f ~ FE, Lum.win_f ~ FE, CT.win_f.iv ~ FE,
                               plum.SD ~ FE),
                       
                       pallc=c(ME.dor_m ~ MF, Lum.dor_m ~ MF, CT.dor_m.iv ~ MF, ME.vent_m ~ MF, Lum.vent_m ~ MF, CT.vent_m.iv ~ MF, ME.win_m ~ MF, Lum.win_m ~ MF, CT.win_m.iv ~ MF,
                               ME.dor_f ~ MF, Lum.dor_f ~ MF, CT.dor_f.iv ~ MF, ME.vent_f ~ MF, Lum.vent_f ~ MF, CT.vent_f.iv ~ MF, ME.win_f ~ MF, Lum.win_f ~ MF, CT.win_f.iv ~ MF,
                               plum.SD ~ MF),
                       
                       pallabc=c(ME.dor_m ~ HE + FE + MF,Lum.dor_m ~ HE + FE + MF, CT.dor_m.iv ~ HE + FE + MF,ME.vent_m ~ HE + FE + MF, Lum.vent_m ~ HE + FE + MF,CT.vent_m.iv ~ HE + FE + MF, ME.win_m ~ HE + FE + MF,Lum.win_m ~ HE + FE + MF, CT.win_m.iv ~ HE + FE + MF, 
                                 ME.dor_f ~ HE + FE + MF, Lum.dor_f ~ HE + FE + MF, CT.dor_f.iv ~ HE + FE + MF, ME.vent_f ~ HE + FE + MF, Lum.vent_f ~ HE + FE + MF, CT.vent_f.iv ~ HE + FE + MF, ME.win_f ~ HE + FE + MF,Lum.win_m ~ HE + FE + MF, CT.win_f.iv ~ HE + FE + MF,
                                 plum.SD ~ HE + FE + MF),
                       
                       valla=c(N.cou_m ~ HE,N.typ_m ~ HE,L.dur_m ~ HE,Peak.F_m ~ HE, N.div_m ~ HE,N.rat_m ~ HE,L.BW_m ~ HE,L.mr_m ~ HE, 
                               N.cou_f ~ HE, N.typ_f ~ HE, L.dur_f ~ HE, Peak.F_f ~ HE, N.div_f ~ HE, N.rat_f ~ HE, L.BW_f ~ HE, L.mr_f ~ HE,
                               vocal.SD ~ HE),
                       
                       vallb=c(N.cou_m ~ FE,N.typ_m ~ FE,L.dur_m ~ FE,Peak.F_m ~ FE, N.div_m ~ FE,N.rat_m ~ FE,L.BW_m ~ FE,L.mr_m ~ FE,
                               N.cou_f ~ FE, N.typ_f ~ FE, L.dur_f ~ FE, Peak.F_f ~ FE, N.div_f ~ FE, N.rat_f ~ FE, L.BW_f ~ FE, L.mr_f ~ FE,
                               vocal.SD ~ FE),
                       
                       vallc=c(N.cou_m ~ MF, N.typ_m ~ MF, L.dur_m ~ MF, Peak.F_m ~ MF, N.div_m ~ MF,N.rat_m ~ MF,L.BW_m ~ MF,L.mr_m ~ MF, 
                               N.cou_f ~ MF, N.typ_f ~ MF, L.dur_f ~ MF, Peak.F_f ~ MF, N.div_f ~ MF, N.rat_f ~ MF, L.BW_f ~ MF, L.mr_f ~ MF,
                               vocal.SD ~ MF),
                       
                       vallabc=c(N.cou_m ~ HE + FE + MF, N.typ_m ~ HE + FE + MF, L.dur_m ~ HE + FE + MF, Peak.F_m ~ HE + FE + MF, N.div_m ~ HE + FE + MF, N.rat_m ~ HE + FE + MF, L.BW_m ~ HE + FE + MF, L.mr_m ~ HE + FE + MF, 
                                 N.cou_f ~ HE + FE + MF, N.typ_f ~ HE + FE + MF, L.dur_f ~ HE + FE + MF, Peak.F_f ~ HE + FE + MF, N.div_f ~ HE + FE + MF, N.rat_f ~ HE + FE + MF, L.BW_f ~ HE + FE + MF, L.mr_f ~ HE + FE + MF,
                                 vocal.SD ~ HE + FE + MF),
                       
                       pvalla=c(ME.dor_m ~ HE, Lum.dor_m ~ HE, CT.dor_m.iv ~ HE, ME.vent_m ~ HE, Lum.vent_m ~ HE, CT.vent_m.iv ~ HE, ME.win_m ~ HE, Lum.win_m ~ HE, CT.win_m.iv ~ HE,
                                N.cou_m ~ HE, N.typ_m ~ HE,L.dur_m ~ HE, Peak.F_m ~ HE, N.div_m ~ HE, N.rat_m ~ HE, L.BW_m ~ HE, L.mr_m ~ HE, 
                                ME.dor_f ~ HE, Lum.dor_f ~ HE, CT.dor_f.iv ~ HE, ME.vent_f ~ HE, Lum.vent_f ~ HE, CT.vent_f.iv ~ HE, ME.win_f ~ HE, Lum.win_f ~ HE, CT.win_f.iv ~ HE, 
                                N.cou_f ~ HE, N.typ_f ~ HE, L.dur_f ~ HE, Peak.F_f ~ HE, N.div_f ~ HE, N.rat_f ~ HE, L.BW_f ~ HE, L.mr_f ~ HE,
                                plum.SD ~ HE, 
                                vocal.SD ~ HE),
                       
                       pvallb=c(ME.dor_m ~ FE, Lum.dor_m ~ FE, CT.dor_m.iv ~ FE, ME.vent_m ~ FE, Lum.vent_m ~ FE, CT.vent_m.iv ~ FE, ME.win_m ~ FE, Lum.win_m ~ FE, CT.win_m.iv ~ FE,
                                N.cou_m ~ FE, N.typ_m ~ FE,L.dur_m ~ FE, Peak.F_m ~ FE, N.div_m ~ FE, N.rat_m ~ FE,L.BW_m ~ FE, L.mr_m ~ FE, 
                                ME.dor_f ~ FE, Lum.dor_f ~ FE, CT.dor_f.iv ~ FE, ME.vent_f ~ FE, Lum.vent_f ~ FE, CT.vent_f.iv ~ FE, ME.win_f ~ FE, Lum.win_f ~ FE, CT.win_f.iv ~ FE, 
                                N.cou_f ~ FE, N.typ_f ~ FE, L.dur_f ~ FE, Peak.F_f ~ FE, N.div_f ~ FE, N.rat_f ~ FE, L.BW_f ~ FE, L.mr_f ~ FE,
                                plum.SD ~ FE, 
                                vocal.SD ~ FE),
                       
                       pvallc=c(ME.dor_m ~ MF, Lum.dor_m ~ MF, CT.dor_m.iv ~ MF, ME.vent_m ~ MF, Lum.vent_m ~ MF, CT.vent_m.iv ~ MF, ME.win_m ~ MF, Lum.win_m ~ MF, CT.win_m.iv ~ MF,
                                N.cou_m ~ MF, N.typ_m ~ MF, L.dur_m ~ MF, Peak.F_m ~ MF, N.div_m ~ MF, N.rat_m ~ MF, L.BW_m ~ MF, L.mr_m ~ MF, 
                                ME.dor_f ~ MF, Lum.dor_f ~ MF, CT.dor_f.iv ~ MF, ME.vent_f ~ MF, Lum.vent_f ~ MF, CT.vent_f.iv ~ MF, ME.win_f ~ MF, Lum.win_f ~ MF, CT.win_f.iv ~ MF, 
                                N.cou_f ~ MF, N.typ_f ~ MF, L.dur_f ~ MF, Peak.F_m ~ MF, N.div_f ~ MF, N.rat_f ~ MF, L.BW_f ~ MF, L.mr_f ~ MF,
                                plum.SD ~ MF, 
                                vocal.SD ~ MF),
                       
                       pvallabc=c(ME.dor_m ~ HE + FE + MF,Lum.dor_m ~ HE + FE + MF, CT.dor_m.iv ~ HE + FE + MF,ME.vent_m ~ HE + FE + MF, Lum.vent_m ~ HE + FE + MF,CT.vent_m.iv ~ HE + FE + MF, ME.win_m ~ HE + FE + MF,Lum.win_m ~ HE + FE + MF, CT.win_m.iv ~ HE + FE + MF,
                                  N.cou_m ~ HE + FE + MF,N.typ_m ~ HE + FE + MF, L.dur_m ~ HE + FE + MF,Peak.F_m ~ HE + FE + MF, N.div_m ~ HE + FE + MF,N.rat_m ~ HE + FE + MF, L.BW_m ~ HE + FE + MF,L.mr_m ~ HE + FE + MF, 
                                  ME.dor_f ~ HE + FE + MF, Lum.dor_f ~ HE + FE + MF, CT.dor_f.iv ~ HE + FE + MF, ME.vent_f ~ HE + FE + MF, Lum.vent_f ~ HE + FE + MF, CT.vent_f.iv ~ HE + FE + MF, ME.win_f ~ HE + FE + MF, Lum.win_f ~ HE + FE + MF, CT.win_f.iv ~ HE + FE + MF, 
                                  N.cou_f ~ HE + FE + MF, N.typ_f ~ HE + FE + MF, L.dur_f ~ HE + FE + MF, Peak.F_f ~ HE + FE + MF, N.div_f ~ HE + FE + MF, N.rat_f ~ HE + FE + MF, L.BW_f ~ HE + FE + MF, L.mr_f ~ HE + FE + MF,
                                  plum.SD ~ HE + FE + MF,
                                  vocal.SD ~ HE + FE + MF))

#All models plotted
plot_model_set(m1.m,text_size=3)


####################################################################################################################################################################################

#########################
#Test models under different evolutionary models:

#1. BM
p.BM <- phylo_path(m1.m, Form, Form.tree,model = 'BM')
p.BM

s.BM <- summary(p.BM)
s.BM

write.table(s.BM,"phylopath_m&f&SD-BM.csv",sep = ",")


plot(s.BM)

#Plot best model
b.BM <- best(p.BM)
plot(b.BM)

#Coeficient values
b.BM$coef
write.table(b.BM$coef,"coef_m&f&SD-BM.csv",sep = ",")

#Standard errors
b.BM$se
write.table(b.BM$se,"SE_m&f&SD-BM.csv",sep = ",")



###################################################################################################
#2. OU
p.OU <- phylo_path(m1.m, Form, Form.tree,model = 'OUfixedRoot')
p.OU

s.OU <- summary(p.OU)
s.OU

write.table(s.OU,"phylopath_m&f&SD-OU.csv",sep = ",")

plot(s.OU)

#Plot best model
b.OU <- best(p.OU)
plot(b.OU)

