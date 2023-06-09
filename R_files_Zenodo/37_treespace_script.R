###############################################################################
##plotting treespace and coding for trees to test the position of KUMIP314087##
###############################################################################

library(ape)
library(phangorn)

#read in the sets of trees from each phylogenetic analysis
trees_max <-read.nexus("13_maxinfo.trprobs.nex",force.multi=TRUE) #read in trees from maxinfo
trees_min <-read.nexus("6_minassumpt.trprobs.nex",force.multi=TRUE) #read in trees from minassumpt
treesBI<-c(trees_max,trees_min) #joint list of trees inferred from maxinfo and minassumpt models
treesMPr <-read.nexus("17_mp_trees.nex",force.multi=TRUE) #read in trees from MP analysis
treesMP <-unroot(treesMPr)
treesBI_MP<-c(treesBI,treesMP) #joint list of trees inferred from BI and MP


#calculate RF distances
rf.KUMIP<-RF.dist(treesBI_MP,rooted = FALSE) 
#run MDS to obtain a treespace
mds<-cmdscale(rf.KUMIP,k=2,add = TRUE,eig = TRUE) 
mds_points_BI_MP <- write (mds$points, "mds_points_BI_MP.csv")


#plot by model (Fig. 4C)
pdf("fig_4c.pdf")
plot(mds$points[,1],mds$points[,2],xlab = 'MDS 1',ylab = 'MDS 2')
points(mds$points[(length(trees_max)+1):length(treesBI),],pch=19,col='grey') #color trees_min
points(mds$points[(length(treesBI)+1):length(treesBI_MP),],cex=2,pch=19,col='black') #color trees_MP
points(mds$points[1:length(trees_max),],pch=1, col='black')  #color trees_max
dev.off()


#vectors can be created to test monophyly of any clade you are interested in
#the following are clades with KUMIP that have >0 trees resultant and tallied in Table 1
KUMIP_op<-c('Opabinia_regalis','KUMIP_314087') #KUMIP with opabinia
pabs_deut<-c('Opabinia_regalis','KUMIP_314087','Fuxianhuia_xiaoshibaensis', 'Chengjiangocaris_kunmingensis', 'Kuamaia_lata', 'Misszhouia_longicaudata', 'Alalcomenaeus_sp', 'Leanchoilia_superlata', 'Limulus_polyphemus', 'Triops_cancriformis') #KUMIP with opabinia + deuteropods
KUMIP_deut<-c('KUMIP_314087', 'Fuxianhuia_xiaoshibaensis', 'Chengjiangocaris_kunmingensis', 'Kuamaia_lata', 'Misszhouia_longicaudata', 'Alalcomenaeus_sp', 'Leanchoilia_superlata', 'Limulus_polyphemus', 'Triops_cancriformis') #KUMIP with deuteropods
KUMIP_radio<-c('KUMIP_314087', 'Aegirocassis_benmoulae', 'Hurdia_victoria', 'Hurdia_triangulata', 'Cambroraster_falcatus', 'Schinderhannes_bartelsi', 'Peytoia_nathorsti', 'Anomalocaris_saron', 'Anomalocaris_canadensis', 'Amplectobelua_symbrachiata', 'Lyrarapax_unguispinus') #KUMIP with radiodonts
pabs_radio_deut<-c('KUMIP_314087','Opabinia_regalis', 'Aegirocassis_benmoulae', 'Hurdia_victoria', 'Hurdia_triangulata', 'Cambroraster_falcatus', 'Schinderhannes_bartelsi', 'Peytoia_nathorsti', 'Anomalocaris_saron', 'Anomalocaris_canadensis', 'Amplectobelua_symbrachiata', 'Lyrarapax_unguispinus', 'Fuxianhuia_xiaoshibaensis', 'Chengjiangocaris_kunmingensis', 'Kuamaia_lata', 'Misszhouia_longicaudata', 'Alalcomenaeus_sp', 'Leanchoilia_superlata', 'Limulus_polyphemus', 'Triops_cancriformis') #KUMIP with opabinia + radiodonts + deuteropods 
KUMIP_radio_deut<-c('KUMIP_314087', 'Aegirocassis_benmoulae', 'Hurdia_victoria', 'Hurdia_triangulata', 'Cambroraster_falcatus', 'Schinderhannes_bartelsi', 'Peytoia_nathorsti', 'Anomalocaris_saron', 'Anomalocaris_canadensis', 'Amplectobelua_symbrachiata', 'Lyrarapax_unguispinus', 'Fuxianhuia_xiaoshibaensis', 'Chengjiangocaris_kunmingensis', 'Kuamaia_lata', 'Misszhouia_longicaudata', 'Alalcomenaeus_sp', 'Leanchoilia_superlata', 'Limulus_polyphemus', 'Triops_cancriformis') #KUMIP with radiodonts + deuteropods
KUMIP_hurdiid<-c('KUMIP_314087', 'Aegirocassis_benmoulae', 'Hurdia_victoria', 'Hurdia_triangulata', 'Cambroraster_falcatus', 'Schinderhannes_bartelsi', 'Peytoia_nathorsti') #KUMIP with hurdiids
KUMIP_pamb<-c('KUMIP_314087','Pambdelurion_whittingtoni') #KUMIP with pambdelurion
#the following are other clades with KUMIP we tested that result in 0 trees
KUMIP_radio_no_schin<-c('KUMIP_314087', 'Aegirocassis_benmoulae', 'Hurdia_victoria', 'Hurdia_triangulata', 'Cambroraster_falcatus', 'Peytoia_nathorsti', 'Anomalocaris_saron', 'Anomalocaris_canadensis', 'Amplectobelua_symbrachiata', 'Lyrarapax_unguispinus') #KUMIP with radiodonts except schinderhannes
KUMIP_hurdiid_no_schin<-c('KUMIP_314087', 'Aegirocassis_benmoulae', 'Hurdia_victoria', 'Hurdia_triangulata', 'Cambroraster_falcatus', 'Peytoia_nathorsti') #KUMIP with hurdiids except schinderhannes
KUMIP_amp_anom<-c('KUMIP_314087', 'Anomalocaris_saron', 'Anomalocaris_canadensis', 'Amplectobelua_symbrachiata', 'Lyrarapax_unguispinus') #KUMIP with amplectobeluids and anomalocaridids
KUMIP_amp<-c('KUMIP_314087', 'Amplectobelua_symbrachiata', 'Lyrarapax_unguispinus') #KUMIP with amplectobeluids
KUMIP_nomms<-c('KUMIP_314087', 'Anomalocaris_saron', 'Anomalocaris_canadensis') #KUMIP with anomalocaridids 
KUMIP_schin<-c('KUMIP_314087', 'Schinderhannes_bartelsi') #KUMIP with schinderhannes
KUMIP_kery<-c('KUMIP_314087','Kerygmachela_kierkegaardi') #KUMIP with kerygmachela
#the following are clades related to radiodonts, not to KUMIP, but for deuteropod and radiodont figures
radio<-c('Aegirocassis_benmoulae', 'Hurdia_victoria', 'Hurdia_triangulata', 'Cambroraster_falcatus', 'Schinderhannes_bartelsi', 'Peytoia_nathorsti', 'Anomalocaris_saron', 'Anomalocaris_canadensis', 'Amplectobelua_symbrachiata', 'Lyrarapax_unguispinus') #monophyletic radiodonts
radio_deut<-c('Aegirocassis_benmoulae', 'Hurdia_victoria', 'Hurdia_triangulata', 'Cambroraster_falcatus', 'Schinderhannes_bartelsi', 'Peytoia_nathorsti', 'Anomalocaris_saron', 'Anomalocaris_canadensis', 'Amplectobelua_symbrachiata', 'Lyrarapax_unguispinus', 'Fuxianhuia_xiaoshibaensis', 'Chengjiangocaris_kunmingensis', 'Kuamaia_lata', 'Misszhouia_longicaudata', 'Alalcomenaeus_sp', 'Leanchoilia_superlata', 'Limulus_polyphemus', 'Triops_cancriformis') #radiodonts with deuteropods
amp_anom_deut<-c('Anomalocaris_saron', 'Anomalocaris_canadensis', 'Amplectobelua_symbrachiata', 'Lyrarapax_unguispinus', 'Fuxianhuia_xiaoshibaensis', 'Chengjiangocaris_kunmingensis', 'Kuamaia_lata', 'Misszhouia_longicaudata', 'Alalcomenaeus_sp', 'Leanchoilia_superlata', 'Limulus_polyphemus', 'Triops_cancriformis') #amplectobeluids + anomalocaridids with deuteropods
hurdiid_deut<-c('Aegirocassis_benmoulae', 'Hurdia_victoria', 'Hurdia_triangulata', 'Cambroraster_falcatus', 'Schinderhannes_bartelsi', 'Peytoia_nathorsti',  'Fuxianhuia_xiaoshibaensis', 'Chengjiangocaris_kunmingensis', 'Kuamaia_lata', 'Misszhouia_longicaudata', 'Alalcomenaeus_sp', 'Leanchoilia_superlata', 'Limulus_polyphemus', 'Triops_cancriformis') #hurdiids with deuteropods
hurdiids<-c('Aegirocassis_benmoulae', 'Hurdia_victoria', 'Hurdia_triangulata', 'Cambroraster_falcatus', 'Schinderhannes_bartelsi', 'Peytoia_nathorsti') #monophyletic hurdiids
amp_anom<-c('Anomalocaris_saron', 'Anomalocaris_canadensis', 'Amplectobelua_symbrachiata', 'Lyrarapax_unguispinus') #amplectobeluids + anomalocaridids monophyletic
#the following are more clades related to radiodonts, not to KUMIP
dino<-c('Opabinia_regalis','KUMIP_314087','Aegirocassis_benmoulae', 'Hurdia_victoria', 'Hurdia_triangulata', 'Cambroraster_falcatus', 'Schinderhannes_bartelsi', 'Peytoia_nathorsti', 'Anomalocaris_saron', 'Anomalocaris_canadensis', 'Amplectobelua_symbrachiata', 'Lyrarapax_unguispinus') #dinocaridids, represents opabiniids with radiodonts
radio_no_schin<-c('Aegirocassis_benmoulae', 'Hurdia_victoria', 'Hurdia_triangulata', 'Cambroraster_falcatus', 'Peytoia_nathorsti', 'Anomalocaris_saron', 'Anomalocaris_canadensis', 'Amplectobelua_symbrachiata', 'Lyrarapax_unguispinus') #monophyletic radiodonts excluding schinderhannes
radio_no_schin_deut<-c('Aegirocassis_benmoulae', 'Hurdia_victoria', 'Hurdia_triangulata', 'Cambroraster_falcatus', 'Peytoia_nathorsti', 'Anomalocaris_saron', 'Anomalocaris_canadensis', 'Amplectobelua_symbrachiata', 'Lyrarapax_unguispinus', 'Fuxianhuia_xiaoshibaensis', 'Chengjiangocaris_kunmingensis', 'Kuamaia_lata', 'Misszhouia_longicaudata', 'Alalcomenaeus_sp', 'Leanchoilia_superlata', 'Limulus_polyphemus', 'Triops_cancriformis') #radiodonts excluding schinderhannes, with deuteropods

#query trees to find monophyly of the clades named above
mono.KUMIP_op<-vector() #KUMIP with opabinia
for (i in 1:length(treesBI_MP)){
  mono.KUMIP_op[i]<-is.monophyletic(treesBI_MP[[i]],tips = KUMIP_op,reroot = TRUE)
}
mono.pabs_deut<-vector() #KUMIP with opabinia + deuteropods
for (i in 1:length(treesBI_MP)){
  mono.pabs_deut[i]<-is.monophyletic(treesBI_MP[[i]],tips = pabs_deut,reroot = TRUE)
}
mono.KUMIP_deut<-vector() #KUMIP with deuteropods
for (i in 1:length(treesBI_MP)){
  mono.KUMIP_deut[i]<-is.monophyletic(treesBI_MP[[i]],tips = KUMIP_deut,reroot = TRUE)
}
mono.KUMIP_radio<-vector() #KUMIP with radiodonts
for (i in 1:length(treesBI_MP)){
  mono.KUMIP_radio[i]<-is.monophyletic(treesBI_MP[[i]],tips = KUMIP_radio,reroot = TRUE)
}
mono.pabs_radio_deut<-vector() #KUMIP with opabinia + radiodonts + deuteropods
for (i in 1:length(treesBI_MP)){
  mono.pabs_radio_deut[i]<-is.monophyletic(treesBI_MP[[i]],tips = pabs_radio_deut,reroot = TRUE)
}
mono.KUMIP_radio_deut<-vector() #KUMIP with radiodonts + deuteropods
for (i in 1:length(treesBI_MP)){
  mono.KUMIP_radio_deut[i]<-is.monophyletic(treesBI_MP[[i]],tips = KUMIP_radio_deut,reroot = TRUE)
}
mono.KUMIP_hurdiid<-vector() #KUMIP with hurdiids
for (i in 1:length(treesBI_MP)){
  mono.KUMIP_hurdiid[i]<-is.monophyletic(treesBI_MP[[i]],tips = KUMIP_hurdiid,reroot = TRUE)
}
mono.KUMIP_pamb<-vector() #KUMIP with pambdelurion
for (i in 1:length(treesBI_MP)){
  mono.KUMIP_pamb[i]<-is.monophyletic(treesBI_MP[[i]],tips = KUMIP_pamb,reroot = TRUE)
}
mono.KUMIP_radio_no_schin<-vector() #KUMIP with radiodonts except schinderhannes
for (i in 1:length(treesBI_MP)){
  mono.KUMIP_radio_no_schin[i]<-is.monophyletic(treesBI_MP[[i]],tips = KUMIP_radio_no_schin,reroot = TRUE)
}
mono.KUMIP_hurdiid_no_schin<-vector()  #KUMIP with hurdiids except schinderhannes
for (i in 1:length(treesBI_MP)){
  mono.KUMIP_hurdiid_no_schin[i]<-is.monophyletic(treesBI_MP[[i]],tips = KUMIP_hurdiid_no_schin,reroot = TRUE)
}
mono.KUMIP_amp_anom<-vector() #KUMIP with amplectobeluids and anomalocaridids
for (i in 1:length(treesBI_MP)){
  mono.KUMIP_amp_anom[i]<-is.monophyletic(treesBI_MP[[i]],tips = KUMIP_amp_anom,reroot = TRUE)
}
mono.KUMIP_amp<-vector() #KUMIP with amplectobeluids 
for (i in 1:length(treesBI_MP)){
  mono.KUMIP_amp[i]<-is.monophyletic(treesBI_MP[[i]],tips = KUMIP_amp,reroot = TRUE)
}
mono.KUMIP_nomms<-vector() #KUMIP with anomalocaridids
for (i in 1:length(treesBI_MP)){
  mono.KUMIP_nomms[i]<-is.monophyletic(treesBI_MP[[i]],tips = KUMIP_nomms,reroot = TRUE)
}
mono.KUMIP_schin<-vector()  #KUMIP with schinderhannes
for (i in 1:length(treesBI_MP)){
  mono.KUMIP_schin[i]<-is.monophyletic(treesBI_MP[[i]],tips = KUMIP_schin,reroot = TRUE)
}
mono.KUMIP_kery<-vector() #KUMIP with kerygmachela
for (i in 1:length(treesBI_MP)){
  mono.KUMIP_kery[i]<-is.monophyletic(treesBI_MP[[i]],tips = KUMIP_kery,reroot = TRUE)
}
mono.radio<-vector() #monophyletic radiodonts
for (i in 1:length(treesBI_MP)){
  mono.radio[i]<-is.monophyletic(treesBI_MP[[i]],tips = radio,reroot = TRUE)
}
mono.radio_deut<-vector() #radiodonts with deuteropods
for (i in 1:length(treesBI_MP)){
  mono.radio_deut[i]<-is.monophyletic(treesBI_MP[[i]],tips = radio_deut,reroot = TRUE)
}
mono.amp_anom_deut<-vector() #amplectobeluids + anomalocaridids with deuteropods
for (i in 1:length(treesBI_MP)){
  mono.amp_anom_deut[i]<-is.monophyletic(treesBI_MP[[i]],tips = amp_anom_deut,reroot = TRUE)
}
mono.hurdiid_deut<-vector() #hurdiids with deuteropods
for (i in 1:length(treesBI_MP)){
  mono.hurdiid_deut[i]<-is.monophyletic(treesBI_MP[[i]],tips = hurdiid_deut,reroot = TRUE)
}
mono.hurdiids<-vector() #hurdiids monophyletic
for (i in 1:length(treesBI_MP)){
  mono.hurdiids[i]<-is.monophyletic(treesBI_MP[[i]],tips = hurdiids,reroot = TRUE)
}
mono.amp_anom<-vector() #amplectobeluids + anomalocaridids monophyletic
for (i in 1:length(treesBI_MP)){
  mono.amp_anom[i]<-is.monophyletic(treesBI_MP[[i]],tips = amp_anom,reroot = TRUE)
}
mono.dino<-vector() #dinocaridids
for (i in 1:length(treesBI_MP)){
  mono.dino[i]<-is.monophyletic(treesBI_MP[[i]],tips = dino,reroot = TRUE)
}
mono.radio_no_schin<-vector() #monophyletic radiodonts excluding schinderhannes
for (i in 1:length(treesBI_MP)){
  mono.radio_no_schin[i]<-is.monophyletic(treesBI_MP[[i]],tips = radio_no_schin,reroot = TRUE)
}
mono.radio_no_schin_deut<-vector() #radiodonts excluding schinderhannes, with deuteropods
for (i in 1:length(treesBI_MP)){
  mono.radio_no_schin_deut[i]<-is.monophyletic(treesBI_MP[[i]],tips = radio_no_schin_deut,reroot = TRUE)
}


#count how many trees had monophyly for each clade tested (Table 1 - in the table broken down by model)
sum(mono.KUMIP_op)
sum(mono.pabs_deut)
sum(mono.KUMIP_deut)
sum(mono.KUMIP_radio)
sum(mono.pabs_radio_deut)
sum(mono.KUMIP_radio_deut)
sum(mono.KUMIP_hurdiid)
sum(mono.KUMIP_pamb)
sum(mono.KUMIP_schin)
#count how many trees had monophyly for each clade tested (other clades)
sum(mono.KUMIP_radio_no_schin)
sum(mono.KUMIP_hurdiid_no_schin)
sum(mono.KUMIP_amp_anom)
sum(mono.KUMIP_amp)
sum(mono.KUMIP_nomms)
sum(mono.KUMIP_kery)
sum(mono.radio)
sum(mono.radio_deut)
sum(mono.amp_anom_deut)
sum(mono.hurdiid_deut)
sum(mono.amp_anom)
sum(mono.hurdiids)
sum(mono.dino)
sum(mono.radio_no_schin)
sum(mono.radio_no_schin_deut)


#to see an individual tree, in this case tree 200
is.monophyletic(treesBI_MP[[200]],tips = KUMIP_op,reroot = TRUE,plot = TRUE)
plot.phylo(root(treesBI_MP[[200]],15))


#plot by bipartition for KUMIP relationships with color palette (Fig. 4B)
pdf("fig_4b.pdf")
plot(mds$points[,1],mds$points[,2],xlab = 'MDS 1',ylab = 'MDS 2')
points(mds$points[,1][which(mono.pabs_deut==TRUE)],mds$points[,2][which(mono.pabs_deut==TRUE)],cex=1.5,pch=19,col="#b66dff") #color KUMIP with opabinia + deuteropods
points(mds$points[,1][which(mono.KUMIP_deut==TRUE)],mds$points[,2][which(mono.KUMIP_deut==TRUE)],pch=19,col="#490092")  #color KUMIP with deuteropods
points(mds$points[,1][which(mono.KUMIP_op==TRUE)],mds$points[,2][which(mono.KUMIP_op==TRUE)],cex=1.5,pch=19,col="#ff6db6") #color KUMIP with opabinia
points(mds$points[,1][which(mono.KUMIP_radio_deut==TRUE)],mds$points[,2][which(mono.KUMIP_radio_deut==TRUE)],pch=19,col="#006ddb") #color KUMIP with radiodonts + deuteropods
points(mds$points[,1][which(mono.KUMIP_radio==TRUE)],mds$points[,2][which(mono.KUMIP_radio==TRUE)],cex=1.5,pch=19,col="#24ff24") #color KUMIP with radiodonts
points(mds$points[,1][which(mono.KUMIP_pamb==TRUE)],mds$points[,2][which(mono.KUMIP_pamb==TRUE)],cex=1.5,pch=19,col="#004949") #color KUMIP with pambdelurion
dev.off()


#plot by bipartition for deuteropod relationships (Supp Fig. 5)
pdf("fig_s5c.pdf")
plot(mds$points[,1],mds$points[,2],xlab = 'MDS 1',ylab = 'MDS 2')
points(mds$points[,1][which(mono.pabs_deut==TRUE)],mds$points[,2][which(mono.pabs_deut==TRUE)],pch=19,col="#b66dff") #color KUMIP with opabinia + deuteropods
points(mds$points[,1][which(mono.radio_deut==TRUE)],mds$points[,2][which(mono.radio_deut==TRUE)],cex=2,pch=19,col="#920000") #color radiodonts + deuteropods
points(mds$points[,1][which(mono.amp_anom_deut==TRUE)],mds$points[,2][which(mono.amp_anom_deut==TRUE)],pch=19,col="#ffff6d") #color amplectobeluids + anomalocaridids + deuteropods
points(mds$points[,1][which(mono.hurdiid_deut==TRUE)],mds$points[,2][which(mono.hurdiid_deut==TRUE)],pch=19,col="#ffb6db") #color hurdiids + deuteropods
dev.off()


#plot by bipartition for radiodont relationships (Supp Fig. 7)
pdf("fig_s7.pdf")
plot(mds$points[,1],mds$points[,2],xlab = 'MDS 1',ylab = 'MDS 2')
points(mds$points[,1][which(mono.radio==TRUE)],mds$points[,2][which(mono.radio==TRUE)],cex=2,pch=19,col="#009292") #color monophyletic radiodonts
points(mds$points[,1][which(mono.amp_anom==TRUE)],mds$points[,2][which(mono.amp_anom==TRUE)],pch=19,col="#b6dbff") #color monophyletic amplectobeluids + anomalocaridids
points(mds$points[,1][which(mono.hurdiids==TRUE)],mds$points[,2][which(mono.hurdiids==TRUE)],pch=19,col="#db6d00") #color monophyletic hurdiids
dev.off()


##################################################################################
##analyses with proboscis coded as uncertain, Table 1 right side and Supp Fig. 4##
##################################################################################


#read in the sets of trees with proboscis coded as uncertain
trees_max_np <-read.nexus("31_maxinfo_np.trprobs.nex",force.multi=TRUE) #read in trees from maxinfo with proboscis uncertain
trees_min_np <-read.nexus("24_minassumpt_np.trprobs.nex",force.multi=TRUE) #read in trees from minassumpt with proboscis uncertain
treesBI_np <- c(trees_max_np,trees_min_np) #joint list of trees inferred from maxinfo and minassumpt models with proboscis uncertain
treesMPr_np <-read.nexus("35_mp_np_trees.nex",force.multi=TRUE) #read in trees from MP analysis with proboscis uncertain
treesMP_np <-unroot(treesMPr_np)
treesBI_MP_np <-c(treesBI_np,treesMP_np) #joint list of trees inferred from BI and MP with proboscis uncertain


#calculate RF distances with proboscis coded as uncertain
rf.KUMIP_np<-RF.dist(treesBI_MP_np,rooted = FALSE) 
#run MDS to obtain a treespace with proboscis coded as uncertain
mds_np<-cmdscale(rf.KUMIP_np,k=2,add = TRUE,eig = TRUE) 
mds_points_BI_MP_np <- write (mds_np$points, "mds_points_BI_MP_np.csv")


#plot by model with proboscis coded as uncertain (Supp Fig. 4B)
pdf("fig_s4b.pdf")
plot(mds_np$points[,1],mds_np$points[,2],xlab = 'MDS 1',ylab = 'MDS 2')
points(mds_np$points[(length(trees_max_np)+1):length(treesBI_np),],pch=19,col='grey') #color trees_min
points(mds_np$points[(length(treesBI_np)+1):length(treesBI_MP_np),],cex=2,pch=19,col='black') #color trees_MP
points(mds_np$points[1:length(trees_max_np),],pch=1, col='black')  #color trees_max
dev.off()


#use the original vectors to test monophyly of clades
#query trees to find monophyly of the clades in trees with proboscis coded as uncertain
mono.KUMIP_op_np<-vector() #KUMIP with opabinia with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.KUMIP_op_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = KUMIP_op,reroot = TRUE)
}
mono.pabs_deut_np<-vector() #KUMIP with opabinia + deuteropods with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.pabs_deut_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = pabs_deut,reroot = TRUE)
}
mono.KUMIP_radio_deut_np<-vector() #KUMIP with opabinia + radiodonts + deuteropods
for (i in 1:length(treesBI_MP_np)){
  mono.KUMIP_radio_deut_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = KUMIP_radio_deut,reroot = TRUE)
}
mono.KUMIP_deut_np<-vector() #KUMIP with deuteropods with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.KUMIP_deut_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = KUMIP_deut,reroot = TRUE)
}
mono.KUMIP_radio_np<-vector() #KUMIP with radiodonts with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.KUMIP_radio_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = KUMIP_radio,reroot = TRUE)
}
mono.KUMIP_hurdiid_np<-vector() #KUMIP with hurdiids with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.KUMIP_hurdiid_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = KUMIP_hurdiid,reroot = TRUE)
}
mono.KUMIP_pamb_np<-vector() #KUMIP with pambdelurion with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.KUMIP_pamb_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = KUMIP_pamb,reroot = TRUE)
}
mono.KUMIP_kery_np<-vector() #KUMIP with kerygmachela with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.KUMIP_kery_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = KUMIP_kery,reroot = TRUE)
}
mono.KUMIP_radio_no_schin_np<-vector() #KUMIP with radiodonts except schinderhannes with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.KUMIP_radio_no_schin_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = KUMIP_radio_no_schin,reroot = TRUE)
}
mono.KUMIP_hurdiid_no_schin_np<-vector()  #KUMIP with hurdiids except schinderhannes with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.KUMIP_hurdiid_no_schin_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = KUMIP_hurdiid_no_schin,reroot = TRUE)
}
mono.KUMIP_amp_anom_np<-vector() #KUMIP with amplectobeluids and anomalocaridids with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.KUMIP_amp_anom_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = KUMIP_amp_anom,reroot = TRUE)
}
mono.KUMIP_amp<-vector_np() #KUMIP with amplectobeluids with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.KUMIP_amp_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = KUMIP_amp,reroot = TRUE)
}
mono.KUMIP_nomms_np<-vector() #KUMIP with anomalocaridids with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.KUMIP_nomms_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = KUMIP_nomms,reroot = TRUE)
}
mono.KUMIP_schin_np<-vector()  #KUMIP with schinderhannes with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.KUMIP_schin_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = KUMIP_schin,reroot = TRUE)
}
mono.radio_np<-vector() #monophyletic radiodonts with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.radio_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = radio,reroot = TRUE)
}
mono.radio_deut_np<-vector() #radiodonts with deuteropods with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.radio_deut_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = radio_deut,reroot = TRUE)
}
mono.amp_anom_deut_np<-vector() #amplectobeluids + anomalocaridids with deuteropods with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.amp_anom_deut_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = amp_anom_deut,reroot = TRUE)
}
mono.hurdiid_deut_np<-vector() #hurdiids with deuteropods with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.hurdiid_deut_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = hurdiid_deut,reroot = TRUE)
}
mono.hurdiids_np<-vector() #hurdiids monophyletic with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.hurdiids_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = hurdiids,reroot = TRUE)
}
mono.amp_anom_np<-vector() #amplectobeluids + anomalocaridids monophyletic with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.amp_anom_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = amp_anom,reroot = TRUE)
}
mono.dino_np<-vector() #dinocaridids with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.dino_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = dino,reroot = TRUE)
}
mono.radio_no_schin_np<-vector() #monophyletic radiodonts excluding schinderhannes with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.radio_no_schin_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = radio_no_schin,reroot = TRUE)
}
mono.radio_no_schin_deut_np<-vector() #radiodonts excluding schinderhannes, with deuteropods with proboscis coded as uncertain
for (i in 1:length(treesBI_MP_np)){
  mono.radio_no_schin_deut_np[i]<-is.monophyletic(treesBI_MP_np[[i]],tips = radio_no_schin_deut,reroot = TRUE)
}


#plot by bipartition for KUMIP relationships with proboscis coded as uncertain (Supp Fig. 4A)
pdf("fig_s4a.pdf")
plot(mds_np$points[,1],mds_np$points[,2],xlab = 'MDS 1',ylab = 'MDS 2')
points(mds_np$points[,1][which(mono.KUMIP_op_np==TRUE)],mds_np$points[,2][which(mono.KUMIP_op_np==TRUE)],cex=1.5,pch=19,col="#ff6db6") #color KUMIP with opabinia with proboscis coded as uncertain 
points(mds_np$points[,1][which(mono.KUMIP_radio_deut_np==TRUE)],mds_np$points[,2][which(mono.KUMIP_radio_deut_np==TRUE)],pch=19,col="#006ddb") #color KUMIP with radiodonts + deuteropods with proboscis coded as uncertain 
points(mds_np$points[,1][which(mono.pabs_deut_np==TRUE)],mds_np$points[,2][which(mono.pabs_deut_np==TRUE)],cex=1.5,pch=19,col="#b66dff") #color KUMIP with opabinia + deuteropods with proboscis coded as uncertain 
points(mds_np$points[,1][which(mono.KUMIP_deut_np==TRUE)],mds_np$points[,2][which(mono.KUMIP_deut_np==TRUE)],pch=19,col="#490092")  #color KUMIP with deuteropods with proboscis coded as uncertain  
points(mds_np$points[,1][which(mono.KUMIP_radio_np==TRUE)],mds_np$points[,2][which(mono.KUMIP_radio_np==TRUE)],cex=1.5,pch=19,col="#24ff24") #color KUMIP with radiodonts with proboscis coded as uncertain 
points(mds_np$points[,1][which(mono.KUMIP_pamb_np==TRUE)],mds_np$points[,2][which(mono.KUMIP_pamb_np==TRUE)],cex=1.5,pch=19,col="#004949") #color KUMIP with pambdelurion with proboscis coded as uncertain 
points(mds_np$points[,1][which(mono.KUMIP_hurdiid_np==TRUE)],mds_np$points[,2][which(mono.KUMIP_hurdiid_np==TRUE)],cex=1.5,pch=19,col="#924900") #color KUMIP with hurdiids with proboscis coded as uncertain 
dev.off()


