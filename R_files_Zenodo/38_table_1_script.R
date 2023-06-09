################################################################
##tallying trees with each position of KUMIP314087 for Table 1##
################################################################

library(ape)

#read in the sets of trees from each phylogenetic analysis
trees_max <-read.nexus("13_maxinfo.trprobs.nex",force.multi=TRUE) #read in trees from maxinfo
trees_min <-read.nexus("6_minassumpt.trprobs.nex",force.multi=TRUE) #read in trees from minassumpt
treesBI<-c(trees_max,trees_min) #joint list of trees inferred from maxinfo and minassumpt models
treesMPr <-read.nexus("17_mp_trees.nex",force.multi=TRUE) #read in trees from MP analysis
treesMP <-unroot(treesMPr)
treesBI_MP<-c(treesBI,treesMP) #joint list of trees inferred from BI and MP

#read in the sets of trees with proboscis coded as uncertain
trees_max_np <-read.nexus("31_maxinfo_np.trprobs.nex",force.multi=TRUE) #read in trees from maxinfo with proboscis uncertain
trees_min_np <-read.nexus("24_minassumpt_np.trprobs.nex",force.multi=TRUE) #read in trees from minassumpt with proboscis uncertain
treesBI_np <- c(trees_max_np,trees_min_np) #joint list of trees inferred from maxinfo and minassumpt models with proboscis uncertain
treesMPr_np <-read.nexus("35_mp_np_trees.nex",force.multi=TRUE) #read in trees from MP analysis with proboscis uncertain
treesMP_np <-unroot(treesMPr_np)
treesBI_MP_np <-c(treesBI_np,treesMP_np) #joint list of trees inferred from BI and MP with proboscis uncertain


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


#query trees to find monophyly of the clades named above

#Table 1, column 1: trees_min
mono.KUMIP_op<-vector() #KUMIP with opabinia
for (i in 1:length(trees_min)){
  mono.KUMIP_op[i]<-is.monophyletic(trees_min[[i]],tips = KUMIP_op,reroot = TRUE)
}
mono.pabs_deut<-vector() #KUMIP with opabinia + deuteropods
for (i in 1:length(trees_min)){
  mono.pabs_deut[i]<-is.monophyletic(trees_min[[i]],tips = pabs_deut,reroot = TRUE)
}
mono.KUMIP_deut<-vector() #KUMIP with deuteropods
for (i in 1:length(trees_min)){
  mono.KUMIP_deut[i]<-is.monophyletic(trees_min[[i]],tips = KUMIP_deut,reroot = TRUE)
}
mono.KUMIP_radio<-vector() #KUMIP with radiodonts
for (i in 1:length(trees_min)){
  mono.KUMIP_radio[i]<-is.monophyletic(trees_min[[i]],tips = KUMIP_radio,reroot = TRUE)
}
mono.pabs_radio_deut<-vector() #KUMIP with opabinia + radiodonts + deuteropods
for (i in 1:length(trees_min)){
  mono.pabs_radio_deut[i]<-is.monophyletic(trees_min[[i]],tips = pabs_radio_deut,reroot = TRUE)
}
mono.KUMIP_radio_deut<-vector() #KUMIP with radiodonts + deuteropods
for (i in 1:length(trees_min)){
  mono.KUMIP_radio_deut[i]<-is.monophyletic(trees_min[[i]],tips = KUMIP_radio_deut,reroot = TRUE)
}
mono.KUMIP_hurdiid<-vector() #KUMIP with hurdiids
for (i in 1:length(trees_min)){
  mono.KUMIP_hurdiid[i]<-is.monophyletic(trees_min[[i]],tips = KUMIP_hurdiid,reroot = TRUE)
}
mono.KUMIP_pamb<-vector() #KUMIP with pambdelurion
for (i in 1:length(trees_min)){
  mono.KUMIP_pamb[i]<-is.monophyletic(trees_min[[i]],tips = KUMIP_pamb,reroot = TRUE)
}
mono.KUMIP_radio_no_schin<-vector() #KUMIP with radiodonts except schinderhannes
for (i in 1:length(trees_min)){
  mono.KUMIP_radio_no_schin[i]<-is.monophyletic(trees_min[[i]],tips = KUMIP_radio_no_schin,reroot = TRUE)
}
mono.KUMIP_hurdiid_no_schin<-vector()  #KUMIP with hurdiids except schinderhannes
for (i in 1:length(trees_min)){
  mono.KUMIP_hurdiid_no_schin[i]<-is.monophyletic(trees_min[[i]],tips = KUMIP_hurdiid_no_schin,reroot = TRUE)
}
mono.KUMIP_amp_anom<-vector() #KUMIP with amplectobeluids and anomalocaridids
for (i in 1:length(trees_min)){
  mono.KUMIP_amp_anom[i]<-is.monophyletic(trees_min[[i]],tips = KUMIP_amp_anom,reroot = TRUE)
}
mono.KUMIP_amp<-vector() #KUMIP with amplectobeluids 
for (i in 1:length(trees_min)){
  mono.KUMIP_amp[i]<-is.monophyletic(trees_min[[i]],tips = KUMIP_amp,reroot = TRUE)
}
mono.KUMIP_nomms<-vector() #KUMIP with anomalocaridids
for (i in 1:length(trees_min)){
  mono.KUMIP_nomms[i]<-is.monophyletic(trees_min[[i]],tips = KUMIP_nomms,reroot = TRUE)
}
mono.KUMIP_schin<-vector()  #KUMIP with schinderhannes
for (i in 1:length(trees_min)){
  mono.KUMIP_schin[i]<-is.monophyletic(trees_min[[i]],tips = KUMIP_schin,reroot = TRUE)
}
mono.KUMIP_kery<-vector() #KUMIP with kerygmachela
for (i in 1:length(trees_min)){
  mono.KUMIP_kery[i]<-is.monophyletic(trees_min[[i]],tips = KUMIP_kery,reroot = TRUE)
}

#queries to get "none of the above" 
mono.net_pabs_deut<-vector() #KUMIP with opabinia + deuteropods, exclusively 
for (i in 1:length(trees_min)){
  mono.net_pabs_deut[i]<-(is.monophyletic(trees_min[[i]],tips = pabs_deut,reroot = TRUE) -
(is.monophyletic(trees_min[[i]],tips = pabs_deut,reroot = TRUE) & is.monophyletic(trees_min[[i]],tips = KUMIP_deut,reroot = TRUE)) - 
(is.monophyletic(trees_min[[i]],tips = pabs_deut,reroot = TRUE) & is.monophyletic(trees_min[[i]],tips = KUMIP_op,reroot = TRUE)))
}

#count how many trees had monophyly for each clade tested 
sum(mono.KUMIP_op)
sum(mono.KUMIP_radio_deut)
sum(mono.pabs_deut)
sum(mono.KUMIP_deut)
sum(mono.KUMIP_radio)
sum(mono.KUMIP_pamb)
sum(mono.KUMIP_kery)
sum(mono.KUMIP_hurdiid)
sum(mono.KUMIP_schin)

#sum for none of the above
length(trees_min) - (sum(mono.KUMIP_op) + sum(mono.net_pabs_deut) + sum(mono.KUMIP_deut) + sum(mono.KUMIP_radio) + sum(mono.KUMIP_radio_deut) + sum(mono.KUMIP_hurdiid) + sum(mono.KUMIP_pamb) + sum(mono.KUMIP_schin))



#use the original vectors to test monophyly of clades
#Table 1, column 2: trees_max
mono.KUMIP_op<-vector() #KUMIP with opabinia
for (i in 1:length(trees_max)){
  mono.KUMIP_op[i]<-is.monophyletic(trees_max[[i]],tips = KUMIP_op,reroot = TRUE)
}
mono.pabs_deut<-vector() #KUMIP with opabinia + deuteropods
for (i in 1:length(trees_max)){
  mono.pabs_deut[i]<-is.monophyletic(trees_max[[i]],tips = pabs_deut,reroot = TRUE)
}
mono.KUMIP_deut<-vector() #KUMIP with deuteropods
for (i in 1:length(trees_max)){
  mono.KUMIP_deut[i]<-is.monophyletic(trees_max[[i]],tips = KUMIP_deut,reroot = TRUE)
}
mono.KUMIP_radio<-vector() #KUMIP with radiodonts
for (i in 1:length(trees_max)){
  mono.KUMIP_radio[i]<-is.monophyletic(trees_max[[i]],tips = KUMIP_radio,reroot = TRUE)
}
mono.pabs_radio_deut<-vector() #KUMIP with opabinia + radiodonts + deuteropods
for (i in 1:length(trees_max)){
  mono.pabs_radio_deut[i]<-is.monophyletic(trees_max[[i]],tips = pabs_radio_deut,reroot = TRUE)
}
mono.KUMIP_radio_deut<-vector() #KUMIP with radiodonts + deuteropods
for (i in 1:length(trees_max)){
  mono.KUMIP_radio_deut[i]<-is.monophyletic(trees_max[[i]],tips = KUMIP_radio_deut,reroot = TRUE)
}
mono.KUMIP_hurdiid<-vector() #KUMIP with hurdiids
for (i in 1:length(trees_max)){
  mono.KUMIP_hurdiid[i]<-is.monophyletic(trees_max[[i]],tips = KUMIP_hurdiid,reroot = TRUE)
}
mono.KUMIP_pamb<-vector() #KUMIP with pambdelurion
for (i in 1:length(trees_max)){
  mono.KUMIP_pamb[i]<-is.monophyletic(trees_max[[i]],tips = KUMIP_pamb,reroot = TRUE)
}
mono.KUMIP_radio_no_schin<-vector() #KUMIP with radiodonts except schinderhannes
for (i in 1:length(trees_max)){
  mono.KUMIP_radio_no_schin[i]<-is.monophyletic(trees_max[[i]],tips = KUMIP_radio_no_schin,reroot = TRUE)
}
mono.KUMIP_hurdiid_no_schin<-vector()  #KUMIP with hurdiids except schinderhannes
for (i in 1:length(trees_max)){
  mono.KUMIP_hurdiid_no_schin[i]<-is.monophyletic(trees_max[[i]],tips = KUMIP_hurdiid_no_schin,reroot = TRUE)
}
mono.KUMIP_amp_anom<-vector() #KUMIP with amplectobeluids and anomalocaridids
for (i in 1:length(trees_max)){
  mono.KUMIP_amp_anom[i]<-is.monophyletic(trees_max[[i]],tips = KUMIP_amp_anom,reroot = TRUE)
}
mono.KUMIP_amp<-vector() #KUMIP with amplectobeluids 
for (i in 1:length(trees_max)){
  mono.KUMIP_amp[i]<-is.monophyletic(trees_max[[i]],tips = KUMIP_amp,reroot = TRUE)
}
mono.KUMIP_nomms<-vector() #KUMIP with anomalocaridids
for (i in 1:length(trees_max)){
  mono.KUMIP_nomms[i]<-is.monophyletic(trees_max[[i]],tips = KUMIP_nomms,reroot = TRUE)
}
mono.KUMIP_schin<-vector()  #KUMIP with schinderhannes
for (i in 1:length(trees_max)){
  mono.KUMIP_schin[i]<-is.monophyletic(trees_max[[i]],tips = KUMIP_schin,reroot = TRUE)
}
mono.KUMIP_kery<-vector() #KUMIP with kerygmachela
for (i in 1:length(trees_max)){
  mono.KUMIP_kery[i]<-is.monophyletic(trees_max[[i]],tips = KUMIP_kery,reroot = TRUE)
}

#queries to get "none of the above" for Table 1
mono.net_pabs_deut<-vector() #KUMIP with opabinia + deuteropods, exclusively 
for (i in 1:length(trees_max)){
  mono.net_pabs_deut[i]<-(is.monophyletic(trees_max[[i]],tips = pabs_deut,reroot = TRUE) -
(is.monophyletic(trees_max[[i]],tips = pabs_deut,reroot = TRUE) & is.monophyletic(trees_max[[i]],tips = KUMIP_deut,reroot = TRUE)) - 
(is.monophyletic(trees_max[[i]],tips = pabs_deut,reroot = TRUE) & is.monophyletic(trees_max[[i]],tips = KUMIP_op,reroot = TRUE)))
}

#count how many trees had monophyly for each clade tested 
sum(mono.KUMIP_op)
sum(mono.KUMIP_radio_deut)
sum(mono.pabs_deut)
sum(mono.KUMIP_deut)
sum(mono.KUMIP_radio)
sum(mono.KUMIP_pamb)
sum(mono.KUMIP_kery)
sum(mono.KUMIP_hurdiid)
sum(mono.KUMIP_schin)

#sum for none of the above
length(trees_max) - (sum(mono.KUMIP_op) + sum(mono.net_pabs_deut) + sum(mono.KUMIP_deut) + sum(mono.KUMIP_radio) + sum(mono.KUMIP_radio_deut) + sum(mono.KUMIP_hurdiid) + sum(mono.KUMIP_pamb) + sum(mono.KUMIP_schin))



#use the original vectors to test monophyly of clades
#Table 1, column 3: treesMP
mono.KUMIP_op<-vector() #KUMIP with opabinia
for (i in 1:length(treesMP)){
  mono.KUMIP_op[i]<-is.monophyletic(treesMP[[i]],tips = KUMIP_op,reroot = TRUE)
}
mono.pabs_deut<-vector() #KUMIP with opabinia + deuteropods
for (i in 1:length(treesMP)){
  mono.pabs_deut[i]<-is.monophyletic(treesMP[[i]],tips = pabs_deut,reroot = TRUE)
}
mono.KUMIP_deut<-vector() #KUMIP with deuteropods
for (i in 1:length(treesMP)){
  mono.KUMIP_deut[i]<-is.monophyletic(treesMP[[i]],tips = KUMIP_deut,reroot = TRUE)
}
mono.KUMIP_radio<-vector() #KUMIP with radiodonts
for (i in 1:length(treesMP)){
  mono.KUMIP_radio[i]<-is.monophyletic(treesMP[[i]],tips = KUMIP_radio,reroot = TRUE)
}
mono.pabs_radio_deut<-vector() #KUMIP with opabinia + radiodonts + deuteropods
for (i in 1:length(treesMP)){
  mono.pabs_radio_deut[i]<-is.monophyletic(treesMP[[i]],tips = pabs_radio_deut,reroot = TRUE)
}
mono.KUMIP_radio_deut<-vector() #KUMIP with radiodonts + deuteropods
for (i in 1:length(treesMP)){
  mono.KUMIP_radio_deut[i]<-is.monophyletic(treesMP[[i]],tips = KUMIP_radio_deut,reroot = TRUE)
}
mono.KUMIP_hurdiid<-vector() #KUMIP with hurdiids
for (i in 1:length(treesMP)){
  mono.KUMIP_hurdiid[i]<-is.monophyletic(treesMP[[i]],tips = KUMIP_hurdiid,reroot = TRUE)
}
mono.KUMIP_pamb<-vector() #KUMIP with pambdelurion
for (i in 1:length(treesMP)){
  mono.KUMIP_pamb[i]<-is.monophyletic(treesMP[[i]],tips = KUMIP_pamb,reroot = TRUE)
}
mono.KUMIP_radio_no_schin<-vector() #KUMIP with radiodonts except schinderhannes
for (i in 1:length(treesMP)){
  mono.KUMIP_radio_no_schin[i]<-is.monophyletic(treesMP[[i]],tips = KUMIP_radio_no_schin,reroot = TRUE)
}
mono.KUMIP_hurdiid_no_schin<-vector()  #KUMIP with hurdiids except schinderhannes
for (i in 1:length(treesMP)){
  mono.KUMIP_hurdiid_no_schin[i]<-is.monophyletic(treesMP[[i]],tips = KUMIP_hurdiid_no_schin,reroot = TRUE)
}
mono.KUMIP_amp_anom<-vector() #KUMIP with amplectobeluids and anomalocaridids
for (i in 1:length(treesMP)){
  mono.KUMIP_amp_anom[i]<-is.monophyletic(treesMP[[i]],tips = KUMIP_amp_anom,reroot = TRUE)
}
mono.KUMIP_amp<-vector() #KUMIP with amplectobeluids 
for (i in 1:length(treesMP)){
  mono.KUMIP_amp[i]<-is.monophyletic(treesMP[[i]],tips = KUMIP_amp,reroot = TRUE)
}
mono.KUMIP_nomms<-vector() #KUMIP with anomalocaridids
for (i in 1:length(treesMP)){
  mono.KUMIP_nomms[i]<-is.monophyletic(treesMP[[i]],tips = KUMIP_nomms,reroot = TRUE)
}
mono.KUMIP_schin<-vector()  #KUMIP with schinderhannes
for (i in 1:length(treesMP)){
  mono.KUMIP_schin[i]<-is.monophyletic(treesMP[[i]],tips = KUMIP_schin,reroot = TRUE)
}
mono.KUMIP_kery<-vector() #KUMIP with kerygmachela
for (i in 1:length(treesMP)){
  mono.KUMIP_kery[i]<-is.monophyletic(treesMP[[i]],tips = KUMIP_kery,reroot = TRUE)
}

#queries to get "none of the above" for Table 1
mono.net_pabs_deut<-vector() #KUMIP with opabinia + deuteropods, exclusively 
for (i in 1:length(treesMP)){
  mono.net_pabs_deut[i]<-(is.monophyletic(treesMP[[i]],tips = pabs_deut,reroot = TRUE) -
(is.monophyletic(treesMP[[i]],tips = pabs_deut,reroot = TRUE) & is.monophyletic(treesMP[[i]],tips = KUMIP_deut,reroot = TRUE)) - 
(is.monophyletic(treesMP[[i]],tips = pabs_deut,reroot = TRUE) & is.monophyletic(treesMP[[i]],tips = KUMIP_op,reroot = TRUE)))
}

#count how many trees had monophyly for each clade tested 
sum(mono.KUMIP_op)
sum(mono.KUMIP_radio_deut)
sum(mono.pabs_deut)
sum(mono.KUMIP_deut)
sum(mono.KUMIP_radio)
sum(mono.KUMIP_pamb)
sum(mono.KUMIP_kery)
sum(mono.KUMIP_hurdiid)
sum(mono.KUMIP_schin)

#sum for none of the above
length(treesMP) - (sum(mono.KUMIP_op) + sum(mono.net_pabs_deut) + sum(mono.KUMIP_deut) + sum(mono.KUMIP_radio) + sum(mono.KUMIP_radio_deut) + sum(mono.KUMIP_hurdiid) + sum(mono.KUMIP_pamb) + sum(mono.KUMIP_schin))



#use the original vectors to test monophyly of clades
#query trees to find monophyly of the clades in trees with proboscis coded as uncertain
#Table 1, column 4: trees_min_np
mono.KUMIP_op<-vector() #KUMIP with opabinia
for (i in 1:length(trees_min_np)){
  mono.KUMIP_op[i]<-is.monophyletic(trees_min_np[[i]],tips = KUMIP_op,reroot = TRUE)
}
mono.pabs_deut<-vector() #KUMIP with opabinia + deuteropods
for (i in 1:length(trees_min_np)){
  mono.pabs_deut[i]<-is.monophyletic(trees_min_np[[i]],tips = pabs_deut,reroot = TRUE)
}
mono.KUMIP_deut<-vector() #KUMIP with deuteropods
for (i in 1:length(trees_min_np)){
  mono.KUMIP_deut[i]<-is.monophyletic(trees_min_np[[i]],tips = KUMIP_deut,reroot = TRUE)
}
mono.KUMIP_radio<-vector() #KUMIP with radiodonts
for (i in 1:length(trees_min_np)){
  mono.KUMIP_radio[i]<-is.monophyletic(trees_min_np[[i]],tips = KUMIP_radio,reroot = TRUE)
}
mono.pabs_radio_deut<-vector() #KUMIP with opabinia + radiodonts + deuteropods
for (i in 1:length(trees_min_np)){
  mono.pabs_radio_deut[i]<-is.monophyletic(trees_min_np[[i]],tips = pabs_radio_deut,reroot = TRUE)
}
mono.KUMIP_radio_deut<-vector() #KUMIP with radiodonts + deuteropods
for (i in 1:length(trees_min_np)){
  mono.KUMIP_radio_deut[i]<-is.monophyletic(trees_min_np[[i]],tips = KUMIP_radio_deut,reroot = TRUE)
}
mono.KUMIP_hurdiid<-vector() #KUMIP with hurdiids
for (i in 1:length(trees_min_np)){
  mono.KUMIP_hurdiid[i]<-is.monophyletic(trees_min_np[[i]],tips = KUMIP_hurdiid,reroot = TRUE)
}
mono.KUMIP_pamb<-vector() #KUMIP with pambdelurion
for (i in 1:length(trees_min_np)){
  mono.KUMIP_pamb[i]<-is.monophyletic(trees_min_np[[i]],tips = KUMIP_pamb,reroot = TRUE)
}
mono.KUMIP_radio_no_schin<-vector() #KUMIP with radiodonts except schinderhannes
for (i in 1:length(trees_min_np)){
  mono.KUMIP_radio_no_schin[i]<-is.monophyletic(trees_min_np[[i]],tips = KUMIP_radio_no_schin,reroot = TRUE)
}
mono.KUMIP_hurdiid_no_schin<-vector()  #KUMIP with hurdiids except schinderhannes
for (i in 1:length(trees_min_np)){
  mono.KUMIP_hurdiid_no_schin[i]<-is.monophyletic(trees_min_np[[i]],tips = KUMIP_hurdiid_no_schin,reroot = TRUE)
}
mono.KUMIP_amp_anom<-vector() #KUMIP with amplectobeluids and anomalocaridids
for (i in 1:length(trees_min_np)){
  mono.KUMIP_amp_anom[i]<-is.monophyletic(trees_min_np[[i]],tips = KUMIP_amp_anom,reroot = TRUE)
}
mono.KUMIP_amp<-vector() #KUMIP with amplectobeluids 
for (i in 1:length(trees_min_np)){
  mono.KUMIP_amp[i]<-is.monophyletic(trees_min_np[[i]],tips = KUMIP_amp,reroot = TRUE)
}
mono.KUMIP_nomms<-vector() #KUMIP with anomalocaridids
for (i in 1:length(trees_min_np)){
  mono.KUMIP_nomms[i]<-is.monophyletic(trees_min_np[[i]],tips = KUMIP_nomms,reroot = TRUE)
}
mono.KUMIP_schin<-vector()  #KUMIP with schinderhannes
for (i in 1:length(trees_min_np)){
  mono.KUMIP_schin[i]<-is.monophyletic(trees_min_np[[i]],tips = KUMIP_schin,reroot = TRUE)
}
mono.KUMIP_kery<-vector() #KUMIP with kerygmachela
for (i in 1:length(trees_min_np)){
  mono.KUMIP_kery[i]<-is.monophyletic(trees_min_np[[i]],tips = KUMIP_kery,reroot = TRUE)
}

#queries to get "none of the above" for Table 1
mono.net_pabs_deut<-vector() #KUMIP with opabinia + deuteropods, exclusively 
for (i in 1:length(trees_min_np)){
  mono.net_pabs_deut[i]<-(is.monophyletic(trees_min_np[[i]],tips = pabs_deut,reroot = TRUE) -
(is.monophyletic(trees_min_np[[i]],tips = pabs_deut,reroot = TRUE) & is.monophyletic(trees_min_np[[i]],tips = KUMIP_deut,reroot = TRUE)) - 
(is.monophyletic(trees_min_np[[i]],tips = pabs_deut,reroot = TRUE) & is.monophyletic(trees_min_np[[i]],tips = KUMIP_op,reroot = TRUE)))
}

#count how many trees had monophyly for each clade tested 
sum(mono.KUMIP_op)
sum(mono.KUMIP_radio_deut)
sum(mono.pabs_deut)
sum(mono.KUMIP_deut)
sum(mono.KUMIP_radio)
sum(mono.KUMIP_pamb)
sum(mono.KUMIP_kery)
sum(mono.KUMIP_hurdiid)
sum(mono.KUMIP_schin)

#sum for none of the above
length(trees_min_np) - (sum(mono.KUMIP_op) + sum(mono.net_pabs_deut) + sum(mono.KUMIP_deut) + sum(mono.KUMIP_radio) + sum(mono.KUMIP_radio_deut) + sum(mono.KUMIP_hurdiid) + sum(mono.KUMIP_pamb) + sum(mono.KUMIP_schin))



#use the original vectors to test monophyly of clades
#Table 1, column 5: trees_max_np
mono.KUMIP_op<-vector() #KUMIP with opabinia
for (i in 1:length(trees_max_np)){
  mono.KUMIP_op[i]<-is.monophyletic(trees_max_np[[i]],tips = KUMIP_op,reroot = TRUE)
}
mono.pabs_deut<-vector() #KUMIP with opabinia + deuteropods
for (i in 1:length(trees_max_np)){
  mono.pabs_deut[i]<-is.monophyletic(trees_max_np[[i]],tips = pabs_deut,reroot = TRUE)
}
mono.KUMIP_deut<-vector() #KUMIP with deuteropods
for (i in 1:length(trees_max_np)){
  mono.KUMIP_deut[i]<-is.monophyletic(trees_max_np[[i]],tips = KUMIP_deut,reroot = TRUE)
}
mono.KUMIP_radio<-vector() #KUMIP with radiodonts
for (i in 1:length(trees_max_np)){
  mono.KUMIP_radio[i]<-is.monophyletic(trees_max_np[[i]],tips = KUMIP_radio,reroot = TRUE)
}
mono.pabs_radio_deut<-vector() #KUMIP with opabinia + radiodonts + deuteropods
for (i in 1:length(trees_max_np)){
  mono.pabs_radio_deut[i]<-is.monophyletic(trees_max_np[[i]],tips = pabs_radio_deut,reroot = TRUE)
}
mono.KUMIP_radio_deut<-vector() #KUMIP with radiodonts + deuteropods
for (i in 1:length(trees_max_np)){
  mono.KUMIP_radio_deut[i]<-is.monophyletic(trees_max_np[[i]],tips = KUMIP_radio_deut,reroot = TRUE)
}
mono.KUMIP_hurdiid<-vector() #KUMIP with hurdiids
for (i in 1:length(trees_max_np)){
  mono.KUMIP_hurdiid[i]<-is.monophyletic(trees_max_np[[i]],tips = KUMIP_hurdiid,reroot = TRUE)
}
mono.KUMIP_pamb<-vector() #KUMIP with pambdelurion
for (i in 1:length(trees_max_np)){
  mono.KUMIP_pamb[i]<-is.monophyletic(trees_max_np[[i]],tips = KUMIP_pamb,reroot = TRUE)
}
mono.KUMIP_radio_no_schin<-vector() #KUMIP with radiodonts except schinderhannes
for (i in 1:length(trees_max_np)){
  mono.KUMIP_radio_no_schin[i]<-is.monophyletic(trees_max_np[[i]],tips = KUMIP_radio_no_schin,reroot = TRUE)
}
mono.KUMIP_hurdiid_no_schin<-vector()  #KUMIP with hurdiids except schinderhannes
for (i in 1:length(trees_max_np)){
  mono.KUMIP_hurdiid_no_schin[i]<-is.monophyletic(trees_max_np[[i]],tips = KUMIP_hurdiid_no_schin,reroot = TRUE)
}
mono.KUMIP_amp_anom<-vector() #KUMIP with amplectobeluids and anomalocaridids
for (i in 1:length(trees_max_np)){
  mono.KUMIP_amp_anom[i]<-is.monophyletic(trees_max_np[[i]],tips = KUMIP_amp_anom,reroot = TRUE)
}
mono.KUMIP_amp<-vector() #KUMIP with amplectobeluids 
for (i in 1:length(trees_max_np)){
  mono.KUMIP_amp[i]<-is.monophyletic(trees_max_np[[i]],tips = KUMIP_amp,reroot = TRUE)
}
mono.KUMIP_nomms<-vector() #KUMIP with anomalocaridids
for (i in 1:length(trees_max_np)){
  mono.KUMIP_nomms[i]<-is.monophyletic(trees_max_np[[i]],tips = KUMIP_nomms,reroot = TRUE)
}
mono.KUMIP_schin<-vector()  #KUMIP with schinderhannes
for (i in 1:length(trees_max_np)){
  mono.KUMIP_schin[i]<-is.monophyletic(trees_max_np[[i]],tips = KUMIP_schin,reroot = TRUE)
}
mono.KUMIP_kery<-vector() #KUMIP with kerygmachela
for (i in 1:length(trees_max_np)){
  mono.KUMIP_kery[i]<-is.monophyletic(trees_max_np[[i]],tips = KUMIP_kery,reroot = TRUE)
}

#queries to get "none of the above" for Table 1
mono.net_pabs_deut<-vector() #KUMIP with opabinia + deuteropods, exclusively 
for (i in 1:length(trees_max_np)){
  mono.net_pabs_deut[i]<-(is.monophyletic(trees_max_np[[i]],tips = pabs_deut,reroot = TRUE) -
(is.monophyletic(trees_max_np[[i]],tips = pabs_deut,reroot = TRUE) & is.monophyletic(trees_max_np[[i]],tips = KUMIP_deut,reroot = TRUE)) - 
(is.monophyletic(trees_max_np[[i]],tips = pabs_deut,reroot = TRUE) & is.monophyletic(trees_max_np[[i]],tips = KUMIP_op,reroot = TRUE)))
}

#count how many trees had monophyly for each clade tested 
sum(mono.KUMIP_op)
sum(mono.KUMIP_radio_deut)
sum(mono.pabs_deut)
sum(mono.KUMIP_deut)
sum(mono.KUMIP_radio)
sum(mono.KUMIP_pamb)
sum(mono.KUMIP_kery)
sum(mono.KUMIP_hurdiid)
sum(mono.KUMIP_schin)

#sum for none of the above
length(trees_max_np) - (sum(mono.KUMIP_op) + sum(mono.net_pabs_deut) + sum(mono.KUMIP_deut) + sum(mono.KUMIP_radio) + sum(mono.KUMIP_radio_deut) + sum(mono.KUMIP_hurdiid) + sum(mono.KUMIP_pamb) + sum(mono.KUMIP_schin))



#use the original vectors to test monophyly of clades
#Table 1, column 6: treesMP_np
mono.KUMIP_op<-vector() #KUMIP with opabinia
for (i in 1:length(treesMP_np)){
  mono.KUMIP_op[i]<-is.monophyletic(treesMP_np[[i]],tips = KUMIP_op,reroot = TRUE)
}
mono.pabs_deut<-vector() #KUMIP with opabinia + deuteropods
for (i in 1:length(treesMP_np)){
  mono.pabs_deut[i]<-is.monophyletic(treesMP_np[[i]],tips = pabs_deut,reroot = TRUE)
}
mono.KUMIP_deut<-vector() #KUMIP with deuteropods
for (i in 1:length(treesMP_np)){
  mono.KUMIP_deut[i]<-is.monophyletic(treesMP_np[[i]],tips = KUMIP_deut,reroot = TRUE)
}
mono.KUMIP_radio<-vector() #KUMIP with radiodonts
for (i in 1:length(treesMP_np)){
  mono.KUMIP_radio[i]<-is.monophyletic(treesMP_np[[i]],tips = KUMIP_radio,reroot = TRUE)
}
mono.pabs_radio_deut<-vector() #KUMIP with opabinia + radiodonts + deuteropods
for (i in 1:length(treesMP_np)){
  mono.pabs_radio_deut[i]<-is.monophyletic(treesMP_np[[i]],tips = pabs_radio_deut,reroot = TRUE)
}
mono.KUMIP_radio_deut<-vector() #KUMIP with radiodonts + deuteropods
for (i in 1:length(treesMP_np)){
  mono.KUMIP_radio_deut[i]<-is.monophyletic(treesMP_np[[i]],tips = KUMIP_radio_deut,reroot = TRUE)
}
mono.KUMIP_hurdiid<-vector() #KUMIP with hurdiids
for (i in 1:length(treesMP_np)){
  mono.KUMIP_hurdiid[i]<-is.monophyletic(treesMP_np[[i]],tips = KUMIP_hurdiid,reroot = TRUE)
}
mono.KUMIP_pamb<-vector() #KUMIP with pambdelurion
for (i in 1:length(treesMP_np)){
  mono.KUMIP_pamb[i]<-is.monophyletic(treesMP_np[[i]],tips = KUMIP_pamb,reroot = TRUE)
}
mono.KUMIP_radio_no_schin<-vector() #KUMIP with radiodonts except schinderhannes
for (i in 1:length(treesMP_np)){
  mono.KUMIP_radio_no_schin[i]<-is.monophyletic(treesMP_np[[i]],tips = KUMIP_radio_no_schin,reroot = TRUE)
}
mono.KUMIP_hurdiid_no_schin<-vector()  #KUMIP with hurdiids except schinderhannes
for (i in 1:length(treesMP_np)){
  mono.KUMIP_hurdiid_no_schin[i]<-is.monophyletic(treesMP_np[[i]],tips = KUMIP_hurdiid_no_schin,reroot = TRUE)
}
mono.KUMIP_amp_anom<-vector() #KUMIP with amplectobeluids and anomalocaridids
for (i in 1:length(treesMP_np)){
  mono.KUMIP_amp_anom[i]<-is.monophyletic(treesMP_np[[i]],tips = KUMIP_amp_anom,reroot = TRUE)
}
mono.KUMIP_amp<-vector() #KUMIP with amplectobeluids 
for (i in 1:length(treesMP_np)){
  mono.KUMIP_amp[i]<-is.monophyletic(treesMP_np[[i]],tips = KUMIP_amp,reroot = TRUE)
}
mono.KUMIP_nomms<-vector() #KUMIP with anomalocaridids
for (i in 1:length(treesMP_np)){
  mono.KUMIP_nomms[i]<-is.monophyletic(treesMP_np[[i]],tips = KUMIP_nomms,reroot = TRUE)
}
mono.KUMIP_schin<-vector()  #KUMIP with schinderhannes
for (i in 1:length(treesMP_np)){
  mono.KUMIP_schin[i]<-is.monophyletic(treesMP_np[[i]],tips = KUMIP_schin,reroot = TRUE)
}
mono.KUMIP_kery<-vector() #KUMIP with kerygmachela
for (i in 1:length(treesMP_np)){
  mono.KUMIP_kery[i]<-is.monophyletic(treesMP_np[[i]],tips = KUMIP_kery,reroot = TRUE)
}

#queries to get "none of the above" for Table 1
mono.net_pabs_deut<-vector() #KUMIP with opabinia + deuteropods, exclusively 
for (i in 1:length(treesMP_np)){
  mono.net_pabs_deut[i]<-(is.monophyletic(treesMP_np[[i]],tips = pabs_deut,reroot = TRUE) -
(is.monophyletic(treesMP_np[[i]],tips = pabs_deut,reroot = TRUE) & is.monophyletic(treesMP_np[[i]],tips = KUMIP_deut,reroot = TRUE)) - 
(is.monophyletic(treesMP_np[[i]],tips = pabs_deut,reroot = TRUE) & is.monophyletic(treesMP_np[[i]],tips = KUMIP_op,reroot = TRUE)))
}

#count how many trees had monophyly for each clade tested 
sum(mono.KUMIP_op)
sum(mono.KUMIP_radio_deut)
sum(mono.pabs_deut)
sum(mono.KUMIP_deut)
sum(mono.KUMIP_radio)
sum(mono.KUMIP_pamb)
sum(mono.KUMIP_kery)
sum(mono.KUMIP_hurdiid)
sum(mono.KUMIP_schin)

#sum for none of the above
length(treesMP_np) - (sum(mono.KUMIP_op) + sum(mono.net_pabs_deut) + sum(mono.KUMIP_deut) + sum(mono.KUMIP_radio) + sum(mono.KUMIP_radio_deut) + sum(mono.KUMIP_hurdiid) + sum(mono.KUMIP_pamb) + sum(mono.KUMIP_schin))



