
#####################################################################################
#	Analysis for Schumacher and Carlson 2022										#
#																					#
#                Combines new otophysan data with published                         #
#				 osetoglossiform data from Sukhum et al 2018						#
#																					#
#				 Requires: SchumacherAndCarlson_SupplementaryFile1.csv				#
#             			   ConvMosEvol_tree.nex										#
#						   actinopt_12k_treePL.tre									#
#						   31Calibrations_Chronogram_newick.tre						#
#						   phylo.fda.v0.2.R											#
#				 Revised May 2022													#
#																					#
#####################################################################################


## wd=eLife Submission
#setwd("PATH");

library(plyr)
library(tidyverse)
library(car)
library(phytools)
library(ape)
library(nlme)
library(emmeans)
library(MuMIn)
library(bayou)
`%!in%` = Negate(`%in%`)  ##code for making the operator not in used below

## reading in the data
comb_log_data<-read.csv(file="SchumacherAndCarlson_SupplementaryFile1.csv", header=TRUE, sep=",");

## Make tbv-roi and append all of that to comb_log_data 
tb.ob<-comb_log_data$Total.Brain.Volume - comb_log_data$OB
tb.tel<-comb_log_data$Total.Brain.Volume - comb_log_data$TEL
tb.hb<-comb_log_data$Total.Brain.Volume - comb_log_data$HB
tb.ot<-comb_log_data$Total.Brain.Volume - comb_log_data$OT
tb.ts<-comb_log_data$Total.Brain.Volume - comb_log_data$TS
tb.cb<-comb_log_data$Total.Brain.Volume - comb_log_data$CB
tb.rob<-comb_log_data$Total.Brain.Volume - comb_log_data$RoB
comb_log_data<-cbind(comb_log_data,tb.ob,tb.tel,tb.hb,tb.ot,tb.ts,tb.cb,tb.rob)

comb_log_data[,4:12]<-log10(comb_log_data[,4:12]);
comb_log_data[,17:23]<-log10(comb_log_data[,17:23]);
comb_log_data$Species<-as.factor(comb_log_data$Species);
comb_log_data$ID<-as.factor(comb_log_data$ID);
comb_log_data$Family<-as.factor(comb_log_data$Family);

spp_shapes <- c("Apteronotus_albifrons" = 21, "Brachyhypopomus_gauderio" = 21, "Brevimyrus_niger" = 21, "Brienomyrus_brachyistius" = 22, "Campylomormyrus_spp" = 23, "Chitala_ornata" = 21, "Corydoras_sterbai" = 22, "Danio_rerio" = 21, "Devario_aequipinnatus" = 22, "Eigenmannia_limbata" = 22, "Eigenmannia_virescens" = 23, "Electrophorus_spp" = 22, "Gnathonemus_petersii" = 24, "Gymnarchus_niloticus" = 25, "Gymnocorymbus_ternetzi" = 23, "Gymnotus_carapo" = 23, "Gymnotus_javari" = 24, "Kryptopterus_vitreolus" = 21, "Microglanis_iheringi" = 23, "Mormyrus_tapirus" = 8, "Pantodon_buchholzi" = 22, "Petrocephalus_tenuicauda" = 4, "Phenacogrammus_interruptus" = 24, "Steatogenys_elegans" = 25, "Sternarchella_calhamazon" = 24, "Sternarchella_orthos" = 25, "Sternopygus_macrurus" = 8, "Synodontis_multipunctatus" = 24, "Synodontis_ocellifer" = 4, "Synodontis_petricola" = 25, "Synodontis_soloni" = 8, "Xenomystus_nigri" = 23)   
spp_fills <- c("Apteronotus_albifrons" = "blue", "Brachyhypopomus_gauderio" = "seagreen2", "Brevimyrus_niger" = "mediumorchid1", "Brienomyrus_brachyistius" = "mediumorchid1", "Campylomormyrus_spp" = "mediumorchid1", "Chitala_ornata" = "white", "Corydoras_sterbai" = "gray", "Danio_rerio" = "black", "Devario_aequipinnatus" = "black", "Eigenmannia_limbata" = "blue", "Eigenmannia_virescens" = "blue", "Electrophorus_spp" = "seagreen2", "Gnathonemus_petersii" = "mediumorchid1", "Gymnarchus_niloticus" = "salmon", "Gymnocorymbus_ternetzi" = "black", "Gymnotus_carapo" = "seagreen2", "Gymnotus_javari" = "seagreen2", "Kryptopterus_vitreolus" = "gray", "Microglanis_iheringi" = "gray", "Mormyrus_tapirus" = "mediumorchid1", "Pantodon_buchholzi" = "white", "Petrocephalus_tenuicauda" = "mediumorchid1", "Phenacogrammus_interruptus" = "black", "Steatogenys_elegans" = "seagreen2", "Sternarchella_calhamazon" = "blue", "Sternarchella_orthos" = "blue", "Sternopygus_macrurus" = "blue", "Synodontis_multipunctatus" = "burlywood4", "Synodontis_ocellifer" = "burlywood4", "Synodontis_petricola" = "burlywood4", "Synodontis_soloni" = "burlywood4", "Xenomystus_nigri" = "white")   
sim_colors <- c("Apteronotus_albifrons" = "Black", "Brachyhypopomus_gauderio" = "Black", "Brevimyrus_niger" = "Black", "Brienomyrus_brachyistius" = "Black", "Campylomormyrus_spp" = "Black", "Chitala_ornata" = "Black", "Danio_rerio" = "Black", "Devario_aequipinnatus" = "Black", "Eigenmannia_limbata" = "Black", "Eigenmannia_virescens" = "Black", "Electrophorus_spp" = "Black", "Gnathonemus_petersii" = "Black", "Gymnarchus_niloticus" = "Black", "Gymnocorymbus_ternetzi" = "Black", "Gymnotus_carapo" = "Black", "Gymnotus_javari" = "Black", "Kryptopterus_vitreolus" = "Black", "Mormyrus_tapirus" = "mediumorchid1", "Pantodon_buchholzi" = "Black", "Petrocephalus_tenuicauda" = "mediumorchid1", "Phenacogrammus_interruptus" = "Black", "Steatogenys_elegans" = "Black", "Sternarchella_calhamazon" = "Black", "Sternarchella_orthos" = "Black", "Sternopygus_macrurus" = "Blue", "Synodontis_multipunctatus" = "Black", "Synodontis_petricola" = "Black", "Synodontis_soloni" = "burlywood4", "Xenomystus_nigri" = "Black", "Corydoras_sterbai" = "Black", "Microglanis_iheringi" = "Black", "Synodontis_ocellifer" = "burlywood4")   


######################################
#   species avg matrix for ancovas   #
######################################

## do ancovas on species averages
# caluclating species averages and sd
OB_M <- tapply(comb_log_data$OB, INDEX = comb_log_data$Species, FUN = mean)
OB_sd <- tapply(comb_log_data$OB, INDEX = comb_log_data$Species, FUN = sd)
tmp1<-as.vector(OB_M)
tmp21<-as.vector(OB_sd)

TEL_M <- tapply(comb_log_data$TEL, INDEX = comb_log_data$Species, FUN = mean)
TEL_sd <- tapply(comb_log_data$TEL, INDEX = comb_log_data$Species, FUN = sd)
tmp2<-as.vector(TEL_M)
tmp22<-as.vector(TEL_sd)

HB_M <- tapply(comb_log_data$HB, INDEX = comb_log_data$Species, FUN = mean)
HB_sd <- tapply(comb_log_data$HB, INDEX = comb_log_data$Species, FUN = sd)
tmp4<-as.vector(HB_M)
tmp24<-as.vector(HB_sd)

OT_M <- tapply(comb_log_data$OT, INDEX = comb_log_data$Species, FUN = mean)
OT_sd <- tapply(comb_log_data$OT, INDEX = comb_log_data$Species, FUN = sd)
tmp5<-as.vector(OT_M)
tmp25<-as.vector(OT_sd)

TS_M <- tapply(comb_log_data$TS, INDEX = comb_log_data$Species, FUN = mean)
TS_sd <- tapply(comb_log_data$TS, INDEX = comb_log_data$Species, FUN = sd)
tmp6<-as.vector(TS_M)
tmp26<-as.vector(TS_sd)

CB_M <- tapply(comb_log_data$CB, INDEX = comb_log_data$Species, FUN = mean)
CB_sd <- tapply(comb_log_data$CB, INDEX = comb_log_data$Species, FUN = sd)
tmp7<-as.vector(CB_M)
tmp27<-as.vector(CB_sd)

RoB_M <- tapply(comb_log_data$RoB, INDEX = comb_log_data$Species, FUN = mean)
RoB_sd <- tapply(comb_log_data$RoB, INDEX = comb_log_data$Species, FUN = sd)
tmp8<-as.vector(RoB_M)
tmp28<-as.vector(RoB_sd)

TBV_M <- tapply(comb_log_data$Total.Brain.Volume, INDEX = comb_log_data$Species, FUN = mean)
TBV_sd <- tapply(comb_log_data$Total.Brain.Volume, INDEX = comb_log_data$Species, FUN = sd)
tmp9<-as.vector(TBV_M)
tmp29<-as.vector(TBV_sd)

BM_M <- tapply(comb_log_data$Mass.g, INDEX = comb_log_data$Species, FUN = mean)
BM_sd <- tapply(comb_log_data$Mass.g, INDEX = comb_log_data$Species, FUN = sd)
tmp0<-as.vector(BM_M)
tmp20<-as.vector(BM_sd)

Species <- as.vector(sort(unique(comb_log_data$Species)))
Code_M<-c('Apteronotidae','Hypopomidae','Mormyridae','Mormyridae','Mormyridae','Notopteridae','Callichthyidae','Cyprinidae','Cyprinidae','Sternopygidae','Sternopygidae','Gymnotidae','Mormyridae','Gymnarchidae','Characidae','Gymnotidae','Gymnotidae','Siluridae','Pseudopimelodidae','Mormyridae','Pantodontidae','Mormyridae','Alestidae','Rhamphichthyidae','Apteronotidae','Apteronotidae','Sternopygidae','Mochokidae','Mochokidae','Mochokidae','Mochokidae','Notopteridae')
Type_M<-c('tubegen','tubegen','tubegen','tubegen','tubegen','out2','ampno','out2','out2','tubegen','tubegen','tubegen','tubegen','tubegen','out2','tubegen','tubegen','ampno','ampno','tubegen','out2','tubegen','out2','tubegen','tubegen','tubegen','tubegen','ampegen','ampegen','ampegen','ampegen','ampno')

tmpmatrix<-cbind(Species, tmp0, tmp1, tmp2, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, Code_M, Type_M)
mode(tmpmatrix) = "numeric"
avg_comb_data<-data.frame(tmpmatrix)
avg_comb_data$Species<-Species
avg_comb_data$Code_M<-Code_M
avg_comb_data$Type_M<-Type_M
colnames(avg_comb_data)<-c("Species","BM_M","OB_M","TEL_M","HB_M","OT_M","TS_M","CB_M","RoB_M","TBV_M","Code_M","Type_M")
rownames(avg_comb_data)<-Species

tmpmatrix2<-cbind(Species, tmp20, tmp21, tmp22, tmp24, tmp25, tmp26, tmp27, tmp28, tmp29, Code_M)
mode(tmpmatrix2) = "numeric"
sd_comb_data<-data.frame(tmpmatrix2)
sd_comb_data$Species<-Species
sd_comb_data$Code_M<-Code_M
colnames(sd_comb_data)<-c("Species","BM_sd","OB_sd","TEL_sd","HB_sd","OT_sd","TS_sd","CB_sd","RoB_sd","TBV_sd","Code_M")
rownames(sd_comb_data)<-Species

TBV_OB <- tapply(comb_log_data$tb.ob, INDEX = comb_log_data$Species, FUN = mean)
TBV_TEL <- tapply(comb_log_data$tb.tel, INDEX = comb_log_data$Species, FUN = mean)
TBV_HB <- tapply(comb_log_data$tb.hb, INDEX = comb_log_data$Species, FUN = mean)
TBV_OT <- tapply(comb_log_data$tb.ot, INDEX = comb_log_data$Species, FUN = mean)
TBV_TS <- tapply(comb_log_data$tb.ts, INDEX = comb_log_data$Species, FUN = mean)
TBV_CB <- tapply(comb_log_data$tb.cb, INDEX = comb_log_data$Species, FUN = mean)
TBV_ROB <- tapply(comb_log_data$tb.rob, INDEX = comb_log_data$Species, FUN = mean)

avg_comb_data <- cbind(avg_comb_data,TBV_OB,TBV_TEL,TBV_HB,TBV_OT,TBV_TS,TBV_CB,TBV_ROB)

########################################
#          phylogenetic trees          #
########################################

## subset for electrosensory phenotype
tubegen_sub<-subset(comb_log_data, comb_log_data$Electrogenic=='Yes' & comb_log_data$Electroreceptor=='Tuberous')
ampegen_sub<-subset(comb_log_data, comb_log_data$Electrogenic=='Yes' & comb_log_data$Electroreceptor=='Only_ampullary')
out_sub<-subset(comb_log_data, comb_log_data$Electrogenic=='No')

## subset for lineage
gym_sub<-subset(comb_log_data, comb_log_data$Electrogenic=='Yes' & comb_log_data$Electroreceptor=='Tuberous' & comb_log_data$Lineage=='Otophysan')
morm_sub<-subset(comb_log_data, comb_log_data$Electrogenic=='Yes' & comb_log_data$Electroreceptor=='Tuberous' & comb_log_data$Lineage=='Osteoglossiform')
allgout_sub <- subset(comb_log_data, comb_log_data$Electrogenic=='No' & comb_log_data$Lineage=='Otophysan')
allmout_sub <- subset(comb_log_data, comb_log_data$Electrogenic=='No' & comb_log_data$Lineage=='Osteoglossiform')

tubegen_names<-unique(tubegen_sub$Species)
ampegen_names<-unique(ampegen_sub$Species)
out_names<-unique(out_sub$Species)
gym_names<-unique(gym_sub$Species)
morm_names<-unique(morm_sub$Species)
allgout_names<-unique(allgout_sub$Species)
allmout_names<-unique(allmout_sub$Species)
gwave_names<-c("Apteronotus_albifrons","Eigenmannia_limbata","Eigenmannia_virescens","Sternarchella_calhamazon","Sternarchella_orthos","Sternopygus_macrurus")
gpulse_names<-c("Brachyhypopomus_gauderio","Electrophorus_spp","Gymnotus_carapo","Gymnotus_javari","Steatogenys_elegans")

syno_names<-ampegen_names
cat_names<-c("Corydoras_sterbai","Kryptopterus_vitreolus","Microglanis_iheringi")
gout_names<-subset(allgout_names, allgout_names %!in% cat_names)
mout_names<-subset(allmout_names, allmout_names !=  "Xenomystus_nigri")

## phylogenetic trees
all_names<-unique(comb_log_data$Species)
osteo_names<-comb_log_data %>% subset(Lineage=='Osteoglossiform') %>% .$Species %>% unique()
ostario_names<-comb_log_data %>% subset(Lineage=='Otophysan') %>% .$Species %>% unique()
full_clado_names<-c(as.character(all_names),"Anguilla_anguilla","Chanos_chanos","Clupea_harengus","Elops_saurus")

## trees where siluri are sister to gymnoti
tel_tree_sg<-read.nexus("ConvMosEvol_tree.nex") ## whole tree w 189 species from elopamorpha to ostariophysi
all_tree_sg<-drop.tip(tel_tree_sg, setdiff(tel_tree_sg$tip.label, all_names))  #tree with all study species
osteo_tree_sg<-drop.tip(tel_tree_sg, setdiff(tel_tree_sg$tip.label, osteo_names))
ostario_tree_sg<-drop.tip(tel_tree_sg, setdiff(tel_tree_sg$tip.label, ostario_names))
full_clado_tree_sg<-drop.tip(tel_tree_sg, setdiff(tel_tree_sg$tip.label, full_clado_names))


## actinopterygii phylogeny from Rabosky et al 2018
actino_rab_tree<-read.tree("~/Desktop/Dropbox/Washington University, St. Louis/Carlson Lab/Phylogenetics/Rabosky et al 2018 dryad data/dataFiles/trees/ultrametric_12K/actinopt_12k_treePL.tre")
## correct genus for B niger,
idx<-grep("Brienomyrus_niger",actino_rab_tree$tip.label)
actino_rab_tree$tip.label[idx]<-"Brevimyrus_niger"

## use the krypto with the shortest branch length from the genus node as k vitreolus
idx<-grep("Kryptopterus_limpok",actino_rab_tree$tip.label)
actino_rab_tree$tip.label[idx]<-"Kryptopterus_vitreolus"

## trim tree for my species
all_tree_rab<-drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, all_names))  #tree with all study species

## which of my species are missing from the Rabosky phylogeny?
## 22 species are included, note that electrophorus and campylomormyrus are going to be problematic as they are spp
## Missing: Brevimyrus_niger, Campylomormyrus_spp, Eigenmannia limbata, Electrophorus_spp, Kryptopterus_vitreolus,
##			Pantodon_buchholzi, Petrocephalus_tenuicauda, Steatogenys_elegans, Sternarchella_calhamazon,
##			Sternarchella_orthos
##	B niger is listed as Brienomyrus_niger
##	Tree has 9 Campy species
##	Tree has 5 Kryptopterus species
##	Tree has 21 Petrocephalus species
##	Tree has no Electrophorus species, 3 Eigenmannia but not good for limbata, no Steatogenys, no Sternarchella
##	Other gymnote spp in tree and Tsuboi data: Sternarchorhamphus_muelleri, Apteronotus_leptorhynchus, Orthosternarchus_tamandua, Brachyhypopomus_brevirostris
##	Additional individuals for: A albifrons, G carapo
##	Other characiform spp in tree and Tsuboi data: Parodon_nasus, Pygocentrus_nattereri
##	Other cypriniform spp in tree and Tsuboi data: Carassius_auratus, Carassius_carassius, Cyprinus_carpio, Gobio_gobio
##	Other siluriform spp in tree and Tsuboi data: Diplomystes_nahuelbutaensis, Ictalurus_punctatus, Noturus_flavus, Hypostomus_plecostomoides, Pimelodus_pictus, Plotosus_lineatus
## To check for other species in the full tree: 
# idx<-grep("Brevimyrus",actino_rab_tree$tip.label)
# actino_rab_tree$tip.label[idx]

## what species should be used for Kryptopterus?
#idx<-grep("Kryptopterus",actino_rab_tree$tip.label)
#krypto_names<-actino_rab_tree$tip.label[idx]
#kyrpto_tree_rab<-drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, c(as.character(all_names),krypto_names)))  #tree with all study species
## Answer: Kryptopterus_limpok has the shortest branch lengths of all Kryptos

library(ggpubr)
library(ggplot2)
library(ggtree)
library(cowplot)

## visualizing the tree
ggtree(all_tree_sg) + geom_tiplab(size=3,offset=.2) + xlim(0, 250) + geom_treescale(offset=-1)

## cladogram of ensen phenos
clado_type<-c("tubegen","tubegen","tubegen","tubegen","tubegen","out","oamp","out","out","tubegen","tubegen","tubegen","tubegen","tubegen","out","tubegen","tubegen","oamp","oamp","tubegen","out","tubegen","out","tubegen","tubegen","tubegen","tubegen","ampegen","ampegen","ampegen","ampegen","oamp","out","out","out","out")
names(clado_type)<-full_clado_names

cls <- list(c1=morm_names,
            c2=gym_names,
            c3=syno_names,
            c4=cat_names,
            c5="Xenomystus_nigri",
            c6=gout_names,
            c7=mout_names,
            c8=c("Anguilla_anguilla","Elops_saurus","Clupea_harengus","Chanos_chanos"))
full_clado_tree_sg <- groupOTU(full_clado_tree_sg, cls)
clado_group_vals<-c(c1="darkolivegreen3", c2="darkolivegreen3", c3="darkolivegreen3", c4="darkolivegreen3", c5="darkolivegreen3", c6="black", c7="black", c8="black")
clado_group_outs<-c(c1="magenta", c2="magenta", c3="black", c4="darkolivegreen3", c5="darkolivegreen3", c6="black", c7="black", c8="black")

library(ggnewscale)

paper_clado<-ggtree(full_clado_tree_sg, branch.length="none", size=1.8, aes(colour=group)) + scale_colour_manual(values=c("black",clado_group_outs)) +
	new_scale_color() +
	geom_tree(aes(color=group), size=1) + 
	geom_tiplab(size=3,offset=.2) + xlim(0, 30) + scale_color_manual(values=c("black",clado_group_vals)) +
	geom_cladelabel(node=63, label="Osteoglossiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5) +
	geom_cladelabel(node=66, label="Mormyroids", angle=270,  offset = 12, offset.text=.7, barsize=1.5, hjust=.5) +
	geom_cladelabel(node=42, label="Otophysans", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5) +
	geom_cladelabel(node=51, label="Gymnotiformes", angle=270,  offset = 12, offset.text=.7, barsize=1.5, hjust=.5) +
	geom_cladelabel(node=45, label="Siluriformes", angle=270,  offset = 12, offset.text=.7, barsize=1.5, hjust=.5)
	
paper_clado2<-paper_clado %>% rotate(37) %>% rotate(38) %>% rotate(61) %>% rotate(50) %>% rotate(55) 


## tree of just taxa with brain data
cls2 <- list(c1=morm_names,
            c2=gym_names,
            c3=syno_names,
            c4=cat_names,
            c5="Xenomystus_nigri",
            c6=gout_names,
            c7=mout_names)
all_tree_sg2 <- groupOTU(all_tree_sg, cls2)


paper_clado1.1<-ggtree(all_tree_sg2, branch.length="none", size=1.8, aes(colour=group)) + scale_colour_manual(values=c("black",clado_group_outs)) +
	new_scale_color() +
	geom_tree(aes(color=group), size=1) + 
	geom_tiplab(size=3,offset=.2) + xlim(0, 30) + scale_color_manual(values=c("black",clado_group_vals)) 
	
paper_clado3<-paper_clado1.1 %>% rotate(33) %>% rotate(42) %>% rotate(47) %>% rotate(53) 


############################
#   visualizing rab tree   #
############################

pdf("rab_phylo_ultrametric.pdf", width = 8.5, height = 11)
ggtree(tel_tree_rab, size=.1) + geom_tiplab(size=.5,offset=0,vjust=0)
dev.off()


## rab tel phylogeny colored by order
rab_cols<-c("Acipenseriformes"="antiquewhite3", "Amiiformes" = "aquamarine2", "Anguilliformes" = "aquamarine4", "Atheriniformes" = "azure4", 
	"Aulopiformes" = "black", "Aulopiformes " = "black", "Beloniformes" = "blueviolet", "Beryciformes" = "cadetblue", "Characiformes" = "chartreuse4",
	"Clupeiformes" = "chocolate1", "Cypriniformes" = "darkolivegreen3", "Cyprinodontiformes" = "darksalmon", "Esociformes" = "darkseagreen",
	"Gadiformes" = "darkslateblue", "Gasterosteiformes" = "deeppink", "Gobiesociformes" = "deeppink4", "Gymnotiformes" = "firebrick1",
	"Lepisosteiformes" = "goldenrod1", "Lophiiformes" = "khaki3", "Mugiliformes" = "lavenderblush2", "Myctophiformes" = "lightgreen", 
	"Notacanthiformes" = "lightpink3", "Ophidiiformes" = "lightseagreen", "Osmeriformes" = "magenta", "Osmeriformes " = "magenta",
	"Osteoglossiformes" = "mediumspringgreen", "Perciformes" = "rosybrown3", "Polypteriformes" = "plum", "Saccopharyngiformes" = "tan",
	"Salmoniformes" = "blue", "Scorpaeniformes" = "yellowgreen", "Siluriformes" = "yellow2", "Stephanoberyciformes" = "wheat", "Stomiiformes" = "sienna",
	"Syngnathiformes" = "red4", "Tetraodontiformes" = "lightskyblue", "Zeiformes" = "mediumpurple4")


Pairs<-unique(phylo_tel_data[,c(1,6)])
Ord<-split(Pairs$Genus_Species,Pairs$Order)

## w updated order
Pairs<-unique(phylo_tel_data[,c(1,7)])
Ord<-split(Pairs$Genus_Species,Pairs$Updated_Order)

colorings = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
set.seed(102)
rab_cols1<- sample(colorings, length(Ord))
rab_cols<-c("black",rab_cols1)

tel_tree_rab <- groupOTU(tel_tree_rab, Ord)
rab_tree_order<-ggtree(tel_tree_rab, size=.2, aes(colour=group)) + scale_color_manual(values=rab_cols) + geom_tiplab(size=.5,offset=0,vjust=0)

rab_col_expl<-cbind(summary(Ord),rab_cols1)

pdf("rab_phylo_updated_order_color_ultrametric.pdf", width = 8.5, height = 11)
rab_tree_order
dev.off()


###################################
#  Visualizing Hughes et al tree  #
###################################

hughes_tree<-read.tree("~/Desktop/Dropbox/Washington University, St. Louis/Carlson Lab/Phylogenetics/Hughes et al 2019 dryad data/TreeFiles/31Calibrations_Chronogram_newick.tre")

fig1_names<-c("Polypteridae_Polypterus_bichir","Acipenseridae_Acipenser_sinensis","Amiidae_Amia_calva","Lepisosteidae_Lepisosteus_oculatus",
	"Megalopidae_Megalops_cyprinoides","Anguillidae_Anguilla_anguilla","Morymridae_Gnathonemus_petersii","Clupeidae_Clupea_harengus","Chanidae_Chanos_chanos",
	"Danionidae_Danio_rerio","Serrasalmidae_Pygocentrus_nattereri","Apteronotidae_Apteronotus_albifrons","Callichthyidae_Corydoras_julii","Lepidogalaxiidae_Lepidogalaxias_salamandroides",
	"Esocidae_Esox_lucius","Salmonidae_Salmo_salar","Argentinidae_Argentina_sp","Galaxiidae_Galaxias_maculatus","Stomiidae_Borostomias_antarcticus","Osmeridae_Osmerus_eperlanus",
	"Synodontidae_Synodus_intermedius","Ateleopodidae_Guentherus_altivela","Myctophidae_Benthosema_glaciale","Polymixiidae_Polymixia_japonica","Percopsidae_Percopsis_transmontana",
	"Zeidae_Zeus_faber","Stylephoridae_Stylephorus_chordatus","Gadidae_Gadus_morhua","Lampridae_Lampris_guttatus","Monocentridae_Monocentris_japonica",
	"Berycidae_Beryx_splendens","Holocentridae_Neoniphon_sammara","Ophidiidae_Brotula_barbata")

fig1_tree<-drop.tip(hughes_tree, setdiff(hughes_tree$tip.label, fig1_names))

fig1_tree$tip.label[grep("Polypteridae_Polypterus_bichir",fig1_tree$tip.label)]<-"Polypteriformes"
fig1_tree$tip.label[grep("Acipenseridae_Acipenser_sinensis",fig1_tree$tip.label)]<-"Acipenseriformes"
fig1_tree$tip.label[grep("Amiidae_Amia_calva",fig1_tree$tip.label)]<-"Amiiformes"
fig1_tree$tip.label[grep("Lepisosteidae_Lepisosteus_oculatus",fig1_tree$tip.label)]<-"Lepisosteiformes"
fig1_tree$tip.label[grep("Megalopidae_Megalops_cyprinoides",fig1_tree$tip.label)]<-"Elopiformes"
fig1_tree$tip.label[grep("Anguillidae_Anguilla_anguilla",fig1_tree$tip.label)]<-"Anguilliformes"
fig1_tree$tip.label[grep("Morymridae_Gnathonemus_petersii",fig1_tree$tip.label)]<-"Osteoglossiformes"
fig1_tree$tip.label[grep("Clupeidae_Clupea_harengus",fig1_tree$tip.label)]<-"Clupeiformes"
fig1_tree$tip.label[grep("Chanidae_Chanos_chanos",fig1_tree$tip.label)]<-"Gonorynchiformes"
fig1_tree$tip.label[grep("Danionidae_Danio_rerio",fig1_tree$tip.label)]<-"Cypriniformes"
fig1_tree$tip.label[grep("Serrasalmidae_Pygocentrus_nattereri",fig1_tree$tip.label)]<-"Characiformes"
fig1_tree$tip.label[grep("Apteronotidae_Apteronotus_albifrons",fig1_tree$tip.label)]<-"Gymnotiformes"
fig1_tree$tip.label[grep("Callichthyidae_Corydoras_julii",fig1_tree$tip.label)]<-"Siluriformes"
fig1_tree$tip.label[grep("Lepidogalaxiidae_Lepidogalaxias_salamandroides",fig1_tree$tip.label)]<-"Lepidogalaxiiformes"
fig1_tree$tip.label[grep("Esocidae_Esox_lucius",fig1_tree$tip.label)]<-"Esociformes"
fig1_tree$tip.label[grep("Salmonidae_Salmo_salar",fig1_tree$tip.label)]<-"Salmoniformes"
fig1_tree$tip.label[grep("Argentinidae_Argentina_sp",fig1_tree$tip.label)]<-"Argentiniformes"
fig1_tree$tip.label[grep("Galaxiidae_Galaxias_maculatus",fig1_tree$tip.label)]<-"Galaxiiformes"
fig1_tree$tip.label[grep("Stomiidae_Borostomias_antarcticus",fig1_tree$tip.label)]<-"Stomiatiformes"
fig1_tree$tip.label[grep("Osmeridae_Osmerus_eperlanus",fig1_tree$tip.label)]<-"Osmeriformes"
fig1_tree$tip.label[grep("Synodontidae_Synodus_intermedius",fig1_tree$tip.label)]<-"Aulopiformes"
fig1_tree$tip.label[grep("Ateleopodidae_Guentherus_altivela",fig1_tree$tip.label)]<-"Ateleopodiformes"
fig1_tree$tip.label[grep("Myctophidae_Benthosema_glaciale",fig1_tree$tip.label)]<-"Myctophiformes"
fig1_tree$tip.label[grep("Polymixiidae_Polymixia_japonica",fig1_tree$tip.label)]<-"Polymixiiformes"
fig1_tree$tip.label[grep("Percopsidae_Percopsis_transmontana",fig1_tree$tip.label)]<-"Percopsiformes"
fig1_tree$tip.label[grep("Zeidae_Zeus_faber",fig1_tree$tip.label)]<-"Zeiformes"
fig1_tree$tip.label[grep("Stylephoridae_Stylephorus_chordatus",fig1_tree$tip.label)]<-"Stylephoriformes"
fig1_tree$tip.label[grep("Gadidae_Gadus_morhua",fig1_tree$tip.label)]<-"Gadiformes"
fig1_tree$tip.label[grep("Lampridae_Lampris_guttatus",fig1_tree$tip.label)]<-"Lampriformes"
fig1_tree$tip.label[grep("Monocentridae_Monocentris_japonica",fig1_tree$tip.label)]<-"Trachichthyiformes"
fig1_tree$tip.label[grep("Berycidae_Beryx_splendens",fig1_tree$tip.label)]<-"Beryciformes"
fig1_tree$tip.label[grep("Holocentridae_Neoniphon_sammara",fig1_tree$tip.label)]<-"Holocentriformes"
fig1_tree$tip.label[grep("Ophidiidae_Brotula_barbata",fig1_tree$tip.label)]<-"Percomorphaceae"


paper_fig1_bl<-ggtree(fig1_tree, size=1) +
   	geom_tiplab(size=3,offset=.04, vjust=.4) + 
	theme_tree2(plot.margin=margin(2, 100, 2, 2)) + coord_cartesian(clip='off') 
#	+scale_x_continuous(breaks=c(0:4))
paper_fig1_bl2 <- revts(paper_fig1_bl) +
	scale_x_continuous(breaks=c(-4:0),labels = abs)
paper_fig1_chrono <- paper_fig1_bl2 %>% rotate(34) %>% rotate(39) %>% rotate(46)

## Adding phylopics 
fish_pics <- c("Polypteriformes" = "59c88fc6-9586-4884-81dc-63d392a7944b",
	"Acipenseriformes" = "5ce55e7c-827d-4049-b128-ae0bc8fb2981",
	"Amiiformes" = "c957bd69-9ca4-465e-b86e-f28824ab739d",
	"Lepisosteiformes" = "fd00131d-7916-44bb-9420-2aaac7cf6798",
	"Elopiformes" = "3491012a-91eb-49ea-bd20-479e6b58935c", #made myself - megalops
	"Anguilliformes" = "629e1135-d0c5-4f02-bb13-8cbf63e48692",
	"Osteoglossiformes" = "c2ee7f75-56c3-457e-b55e-56028f105315",
	"Clupeiformes" = "0b89df58-7eae-40c5-9676-d35e3449afb2",
	"Gonorynchiformes" = "7e282f99-579d-4115-a128-f83af138bfab", #made myself - chanos chanos
	"Cypriniformes" = "aa67de30-9e97-45e3-bfcb-bbc0b0d7ef39", 
	"Characiformes" = "56971d2a-799c-4993-9ffc-7455417efd90",
	"Gymnotiformes" = "a6e1b574-cf8d-4714-8f93-9959b649d4f7",
	"Siluriformes" = "ae1d3ca5-f6d9-49dd-a431-e1bea5a3bdf3",
	"Lepidogalaxiiformes" = "debd0fa6-7481-42f9-b100-42eae4a136ce",
	"Esociformes" = "f8369dec-bdf6-432b-a0c4-41ee5d75286d",
	"Salmoniformes" = "9150be88-6910-4374-aa54-a7f8f3d79fb6",
	"Argentiniformes" = "e9f4e72d-377e-44ed-8722-0362fbd589ad",
	"Galaxiiformes" = "57ad745e-81c2-4ba9-a4e2-6d216cfd7325",
	"Stomiatiformes" = "5400612d-8088-470b-a7b0-36ecde18bbaf",
	"Osmeriformes" = "f1f91d08-b850-4600-ad64-622ce87f0199",
	"Aulopiformes" = "f2396149-3019-45a6-b87c-9ac5bfeba8b1", #made myself - synodus
	"Ateleopodiformes" = "0d1d0e12-82ef-447f-8e3d-85155f8f727c", #made myself - ijimaia
	"Myctophiformes" = "93b9ea15-3319-4b74-8003-e4a9f25ec680",
	"Polymixiiformes" = "95214cf9-2647-4dae-8ad3-5cd325818966",
	"Percopsiformes" = "e88bb918-37fe-42a1-8cac-3381f7066035",
	"Zeiformes" = "1cc605fb-013f-4997-b88c-6d2c71952e8d",
	"Stylephoriformes" = "dd082061-8753-48a2-be7c-a2fdc2b45e22", #made myself - Stylephorus
	"Gadiformes" = "dac37960-8810-4e23-8395-758cb9962b47",
	"Lampriformes" = "4dd00c92-a2b7-4169-80df-c0971724b4a3",
	"Trachichthyiformes" = "ddcfb5b4-c1a8-41d3-82e2-201011174d58", #made myself - Anoplaster
	"Beryciformes" = "820b392c-570d-440c-a75b-570c6a4d4dbd", #made myself - Beryx
	"Holocentriformes" = "893f7dc2-d270-44b5-9195-0d9f49c0f670",
	"Percomorphaceae" = "84c7e672-2593-44a6-a807-cffbd3156cc5")


d <- cbind(fig1_tree$tip.label,fish_pics)
colnames(d)<-c("tip.label","fish_pics")
d <- data.frame(d)

phylopic_fig1_chrono <- paper_fig1_chrono %<+% d + geom_tiplab(aes(image=fish_pics), geom="phylopic", offset=.95)

## saving plots to pdf
# Customizing the output
pdf("fig1_phylo.pdf",         # File name
    width = 7, height = 7)	  # Size
phylopic_fig1_chrono
dev.off() 



#############################################################
#          phylogenetic pca w individual rotations          #
#############################################################

## the pca of species means
## creating a species avg matrix
Avg_CB <- tapply(comb_log_data$CB, INDEX = comb_log_data$Species, FUN = mean)
Avg_TEL <- tapply(comb_log_data$TEL, INDEX = comb_log_data$Species, FUN = mean)
Avg_OB <- tapply(comb_log_data$OB, INDEX = comb_log_data$Species, FUN = mean)
Avg_OT <- tapply(comb_log_data$OT, INDEX = comb_log_data$Species, FUN = mean)
Avg_RoB <- tapply(comb_log_data$RoB, INDEX = comb_log_data$Species, FUN = mean)
Avg_HB <- tapply(comb_log_data$HB, INDEX = comb_log_data$Species, FUN = mean)
Avg_TS <- tapply(comb_log_data$TS, INDEX = comb_log_data$Species, FUN = mean)
avg_region_data <- cbind(Avg_OB, Avg_TEL, Avg_HB, Avg_OT, Avg_TS, Avg_CB, Avg_RoB)

zscore_region_data <- scale(avg_region_data, center=TRUE, scale=TRUE)

tree<-all_tree_sg 			## change this to represent which tree I want to use in the pPCA

avg_pca<-phyl.pca(tree,avg_comb_data[,c(3,4,5,6,7,8,9)],method="lambda",mode="corr")

## obtain among species & among trait phylogenetic covariance matrices
## with our fitted lambda
obj<-phyl.vcv(avg_region_data[tree$tip.label,],vcv(tree),avg_pca$lambda)
V<-obj$R
C<-obj$C
n<-nrow(avg_region_data)
m<-ncol(avg_region_data)
## transform the data such that the columns have phylogenetic 
## variances of 1.0
Y<-avg_region_data/matrix(rep(sqrt(diag(V)),n),n,m,byrow=TRUE)
## check (should be 1.0s)
diag(phyl.vcv(as.matrix(Y)[tree$tip.label,],C,1)$R)
## invert among species covariance matrix
invC<-solve(C)[rownames(Y),rownames(Y)]
## compute ancestral states to recenter
a<-matrix(colSums(invC%*%as.matrix(Y))/sum(invC),m,1)
A<-matrix(rep(a,n),n,m,byrow=TRUE)
Y<-Y-A

## compute scores using the phyl.pca rotation on individual data
indiv_dat<-comb_log_data[,c(1,2,3,5,6,7,8,9,10,11)]
indiv_dat2<-as.matrix(comb_log_data[,c(5,6,7,8,9,10,11)])
n<-nrow(indiv_dat2)
Yi<-indiv_dat2/matrix(rep(sqrt(diag(V)),n),n,m,byrow=TRUE)
A<-matrix(rep(a,n),n,m,byrow=TRUE)
Yi<-as.matrix(Yi)-A
Si<-Yi%*%avg_pca$Evec

indiv_ppca<-data.frame(cbind(as.factor(comb_log_data$Species)),Si)
indiv_ppca2<-as.data.frame(indiv_ppca)
indiv_ppca2<-cbind(comb_log_data$Species,comb_log_data$Mass.g,comb_log_data$Total.Brain.Volume,indiv_ppca2[,c(2,3,4,5,6,7,8)],comb_log_data$Electroreceptor,comb_log_data$Electrogenic)
colnames(indiv_ppca2)<-c("Species","Mass.g","Total.Brain.Volume","PC1","PC2","PC3","PC4","PC5","PC6","PC7","recpt","eo")

pers<-summary(avg_pca)
eigs<-pers$sdev^2
percent <- round(eigs / sum(eigs) * 100, 2)
percentage <- paste( colnames(Si), "(", paste( as.character(percent), "%", ")", sep="") )


## Calculating convex polygons
## stat fxn to calculate the polygons
StatChull <- ggproto("StatChull", Stat,
  compute_group = function(data, scales) {
    data[chull(data$x, data$y), , drop = FALSE]
  },
  
  required_aes = c("x", "y")
)
## layer fxn to add the polygons to the ggplot object
stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

## plotting ppca with convex polygons for lineage
ppca_plot_all<-ggplot(indiv_ppca, aes(x=PC1, y=PC2, shape=comb_log_data$Species, fill=comb_log_data$Species, color=comb_log_data$Species, group=comb_log_data$Group)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + ylab(percentage[2]) + xlab(percentage[1]) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	stat_chull(fill=NA, colour="black")
ppca_plot_all_pc23<-ggplot(indiv_ppca, aes(x=PC2, y=PC3, shape=comb_log_data$Species, fill=comb_log_data$Species, color=comb_log_data$Species, group=comb_log_data$Group)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + xlab(percentage[2]) + ylab(percentage[3]) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	stat_chull(fill=NA, colour="black")
ppca_plot_all_pc34<-ggplot(indiv_ppca, aes(x=PC3, y=PC4, shape=comb_log_data$Species, fill=comb_log_data$Species, color=comb_log_data$Species, group=comb_log_data$Group)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + xlab(percentage[3]) + ylab(percentage[4]) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	stat_chull(fill=NA, colour="black")

## plotting ppca with convex hulls for egen+recpt
hyp_eorecpt1<-paste(comb_log_data$Electrogenic,comb_log_data$Electroreceptor)
ppca_plot_all_eorecpt<-ggplot(indiv_ppca, aes(x=PC1, y=PC2, shape=comb_log_data$Species, fill=comb_log_data$Species, color=comb_log_data$Species, group=hyp_eorecpt1)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + ylab(percentage[2]) + xlab(percentage[1]) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	stat_chull(fill=NA, colour="black", aes(lty=hyp_eorecpt1)) + scale_linetype_manual(values=c(2,4,3,1))
ppca_plot_all_pc23_eorecpt<-ggplot(indiv_ppca, aes(x=PC2, y=PC3, shape=comb_log_data$Species, fill=comb_log_data$Species, color=comb_log_data$Species, group=hyp_eorecpt1)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + xlab(percentage[2]) + ylab(percentage[3]) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	stat_chull(fill=NA, colour="black", aes(lty=hyp_eorecpt1)) + scale_linetype_manual(values=c(2,4,3,1))
ppca_plot_all_pc34_eorecpt<-ggplot(indiv_ppca, aes(x=PC3, y=PC4, shape=comb_log_data$Species, fill=comb_log_data$Species, color=comb_log_data$Species, group=hyp_eorecpt1)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + xlab(percentage[3]) + ylab(percentage[4]) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	stat_chull(fill=NA, colour="black", aes(lty=hyp_eorecpt1)) + scale_linetype_manual(values=c(2,4,3,1))


## plotting the loadings
eig_vecs<-avg_pca$Evec
rownames(eig_vecs)<-c("OB","TEL","HB","OT","TS","CB","RoB")
eig_vecs<-as.data.frame(eig_vecs)
ppca_load<-ggplot() + theme_classic() + geom_segment(data=eig_vecs, mapping=aes(x=0, y=0, xend=PC1, yend=PC2), arrow=arrow(type="open", length = unit(0.12, "cm")), size=0.35, color="black") + geom_text(data = eig_vecs, aes(x=PC1, y=PC2, label = row.names(eig_vecs)), nudge_x = c(-.02,-.08,-.06,-.055,-.055,-.02,-.08), nudge_y = c(.05,0,-.02,.02,0,-.05,0), size = 2.5) + coord_cartesian(xlim = c(-.6, .1), ylim = c(-.6, .7)) + xlab("PC1 Eigenvectors") + ylab("PC2 Eigenvectors")
ppca_load_pc23<-ggplot() + theme_classic() + geom_segment(data=eig_vecs, mapping=aes(x=0, y=0, xend=PC2, yend=PC3), arrow=arrow(length = unit(0.12, "cm")), size=0.35, color="black") + geom_text(data = eig_vecs, aes(x=PC2, y=PC3, label = row.names(eig_vecs)), nudge_x = c(.12,.16,-.11,.05,-.09,-.1,.15), nudge_y = c(.05,-.04,.07,-.12,-.07,.07,-.05), size = 3) + coord_cartesian(xlim = c(-.6, .8), ylim = c(-1, .7)) + xlab("PC2 Eigenvectors") + ylab("PC3 Eigenvectors")
ppca_load_pc34<-ggplot() + theme_classic() + geom_segment(data=eig_vecs, mapping=aes(x=0, y=0, xend=PC3, yend=PC4), arrow=arrow(length = unit(0.12, "cm")), size=0.35, color="black") + geom_text(data = eig_vecs, aes(x=PC3, y=PC4, label = row.names(eig_vecs)), nudge_x = c(.17,.17,.16,-.15,0,0,-.21), nudge_y = c(0,.07,.01,0,.07,-.08,.01), size = 2.7) + coord_cartesian(xlim = c(-1, .9), ylim = c(-.9, .7)) + xlab("PC3 Eigenvectors") + ylab("PC4 Eigenvectors")

## insetting the loadings 
ppca_inset<- ppca_plot_all + annotation_custom(ggplotGrob(ppca_load), xmin = 16, xmax = 60, ymin = -24, ymax = -5)
ppca_inset_pc34<- ppca_plot_all_pc34 + annotation_custom(ggplotGrob(ppca_load_pc34+theme_cowplot()+theme(axis.title=element_text(size=10),axis.text=element_text(size=8))), xmin = 6, xmax = 16.5, ymin = 0, ymax = 6.5)
ppca_inset_eorecpt<- ppca_plot_all_eorecpt + annotation_custom(ggplotGrob(ppca_load), xmin = 16, xmax = 60, ymin = -24, ymax = -5)
ppca_inset_pc34_eorecpt<- ppca_plot_all_pc34_eorecpt + annotation_custom(ggplotGrob(ppca_load_pc34+theme_cowplot()+theme(axis.title=element_text(size=10),axis.text=element_text(size=8))), xmin = 6, xmax = 16.5, ymin = 0, ymax = 6.5)

## pc1 & pc2 v total brain volume
morm_pc1_plot<-ggplot(comb_log_data, aes(y=indiv_ppca$PC1, x=Total.Brain.Volume, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('PC1')) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "mediumorchid1", "black", "mediumorchid1", "black", "black", "black", "black", "Blue", "black", "gray", "black", "gray", "black"));
morm_pc2_plot<-ggplot(comb_log_data, aes(y=indiv_ppca$PC2, x=Total.Brain.Volume, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('PC2')) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "mediumorchid1", "black", "mediumorchid1", "black", "black", "black", "black", "Blue", "black", "gray", "black", "gray", "black"));

pc1tbv.lm<-lm(indiv_ppca$PC1 ~ Total.Brain.Volume, data=comb_log_data)
morm_pc1_plot2 <- ggplot(comb_log_data, aes(y=indiv_ppca$PC1, x=Total.Brain.Volume, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('PC1')) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "mediumorchid1", "black", "mediumorchid1", "black", "black", "black", "black", "Blue", "black", "gray", "black", "gray", "black")) + geom_abline(slope = pc1tbv.lm$coefficients[2], intercept = pc1tbv.lm$coefficients[1])

## do a correlation test
shapiro.test(indiv_ppca$PC1)
shapiro.test(comb_log_data$Total.Brain.Volume)
## normal distribution = "pearson", otherwise do "spearman"
cor.test(comb_log_data$Total.Brain.Volume, indiv_ppca$PC1 , method="spearman")


###################################
#     z-score normalized pPCA     #
###################################

zscore_region_data <- scale(avg_region_data, center=TRUE, scale=TRUE)

tree<-all_tree_sg 			## change this to represent which tree I want to use in the pPCA

avg_pca_z<-phyl.pca(tree,zscore_region_data,method="lambda",mode="corr")

## obtain among species & among trait phylogenetic covariance matrices
## with our fitted lambda
obj<-phyl.vcv(zscore_region_data[tree$tip.label,],vcv(tree),avg_pca_z$lambda)
V<-obj$R
C<-obj$C
n<-nrow(zscore_region_data)
m<-ncol(zscore_region_data)
## transform the data such that the columns have phylogenetic 
## variances of 1.0
Y<-zscore_region_data/matrix(rep(sqrt(diag(V)),n),n,m,byrow=TRUE)
## check (should be 1.0s)
diag(phyl.vcv(as.matrix(Y)[tree$tip.label,],C,1)$R)
## invert among species covariance matrix
invC<-solve(C)[rownames(Y),rownames(Y)]
## compute ancestral states to recenter
a<-matrix(colSums(invC%*%as.matrix(Y))/sum(invC),m,1)
A<-matrix(rep(a,n),n,m,byrow=TRUE)
Y<-Y-A

## compute scores using the phyl.pca rotation on individual data
#indiv_dat_z<-comb_log_data[,c(1,2,3,5,6,7,8,9,10,11)]
indiv_dat_z2<-as.matrix(scale(comb_log_data[,c(5,6,7,8,9,10,11)]), center=TRUE, scale=TRUE)
n<-nrow(indiv_dat_z2)
Yi<-indiv_dat_z2/matrix(rep(sqrt(diag(V)),n),n,m,byrow=TRUE)
A<-matrix(rep(a,n),n,m,byrow=TRUE)
Yi<-as.matrix(Yi)-A
Si<-Yi%*%avg_pca_z$Evec

indiv_ppca_z<-data.frame(cbind(as.factor(comb_log_data$Species)),Si)
indiv_ppca_z2<-as.data.frame(indiv_ppca_z)
indiv_ppca_z2<-cbind(comb_log_data$Species,comb_log_data$Mass.g,comb_log_data$Total.Brain.Volume,indiv_ppca_z2[,c(2,3,4,5,6,7,8)],comb_log_data$Electroreceptor,comb_log_data$Electrogenic)
colnames(indiv_ppca_z2)<-c("Species","Mass.g","Total.Brain.Volume","PC1","PC2","PC3","PC4","PC5","PC6","PC7","recpt","eo")

pers_z<-summary(avg_pca_z)
eigs_z<-pers_z$sdev^2
percent_z <- round(eigs_z / sum(eigs_z) * 100, 2)
percentage_z <- paste( colnames(Si), "(", paste( as.character(percent_z), "%", ")", sep="") )



## plotting ppca with convex hulls for egen+recpt
hyp_eorecpt1<-paste(comb_log_data$Electrogenic,comb_log_data$Electroreceptor)
ppca_plot_z_all_eorecpt<-ggplot(indiv_ppca_z, aes(x=PC1, y=PC2, shape=comb_log_data$Species, fill=comb_log_data$Species, color=comb_log_data$Species, group=hyp_eorecpt1)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + ylab(percentage_z[2]) + xlab(percentage_z[1]) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	new_scale_color() +
	stat_chull(fill=NA, aes(colour=hyp_eorecpt1,lty=hyp_eorecpt1)) + scale_linetype_manual(values=c(2,2,3,1)) + scale_color_manual(values=c("grey","black","black","black"))
ppca_plot_z_all_pc23_eorecpt<-ggplot(indiv_ppca_z, aes(x=PC2, y=PC3, shape=comb_log_data$Species, fill=comb_log_data$Species, color=comb_log_data$Species, group=hyp_eorecpt1)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + xlab(percentage_z[2]) + ylab(percentage_z[3]) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	new_scale_color() +
	stat_chull(fill=NA, aes(colour=hyp_eorecpt1,lty=hyp_eorecpt1)) + scale_linetype_manual(values=c(2,2,3,1)) + scale_color_manual(values=c("grey","black","black","black"))
ppca_plot_z_all_pc34_eorecpt<-ggplot(indiv_ppca_z, aes(x=PC3, y=PC4, shape=comb_log_data$Species, fill=comb_log_data$Species, color=comb_log_data$Species, group=hyp_eorecpt1)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + xlab(percentage_z[3]) + ylab(percentage_z[4]) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	new_scale_color() +
	stat_chull(fill=NA, aes(colour=hyp_eorecpt1,lty=hyp_eorecpt1)) + scale_linetype_manual(values=c(2,2,3,1)) + scale_color_manual(values=c("grey","black","black","black"))


## plotting the loadings
eig_vecs_z<-avg_pca_z$Evec
rownames(eig_vecs_z)<-c("OB","TEL","HB","OT","TS","CB","RoB")
eig_vecs_z<-as.data.frame(eig_vecs_z)
ppca_load_z<-ggplot() + theme_classic() + geom_segment(data=eig_vecs_z, mapping=aes(x=0, y=0, xend=PC1, yend=PC2), arrow=arrow(type="open", length = unit(0.12, "cm")), size=0.35, color="black") + geom_text(data = eig_vecs_z, aes(x=PC1, y=PC2, label = row.names(eig_vecs_z)), nudge_x = c(-.02,-.08,-.06,-.055,-.055,-.02,-.08), nudge_y = c(.05,0,-.02,.02,0,-.05,0), size = 2.5) + coord_cartesian(xlim = c(-.6, .1), ylim = c(-.6, .7)) + xlab("PC1 Eigenvectors") + ylab("PC2 Eigenvectors")
ppca_load_z_pc23<-ggplot() + theme_classic() + geom_segment(data=eig_vecs_z, mapping=aes(x=0, y=0, xend=PC2, yend=PC3), arrow=arrow(length = unit(0.12, "cm")), size=0.35, color="black") + geom_text(data = eig_vecs_z, aes(x=PC2, y=PC3, label = row.names(eig_vecs_z)), nudge_x = c(.12,.16,-.11,.05,-.09,-.1,.15), nudge_y = c(.05,-.04,.07,-.12,-.07,.07,-.05), size = 3) + coord_cartesian(xlim = c(-.6, .8), ylim = c(-1, .7)) + xlab("PC2 Eigenvectors") + ylab("PC3 Eigenvectors")
ppca_load_z_pc34<-ggplot() + theme_classic() + geom_segment(data=eig_vecs_z, mapping=aes(x=0, y=0, xend=PC3, yend=PC4), arrow=arrow(length = unit(0.12, "cm")), size=0.35, color="black") + geom_text(data = eig_vecs_z, aes(x=PC3, y=PC4, label = row.names(eig_vecs_z)), nudge_x = c(.17,.17,.16,-.15,0,0,-.21), nudge_y = c(0,.07,.01,0,.07,-.08,.01), size = 2.7) + coord_cartesian(xlim = c(-1, .9), ylim = c(-.9, .7)) + xlab("PC3 Eigenvectors") + ylab("PC4 Eigenvectors")

## insetting the loadings 
ppca_inset_z_eorecpt<- ppca_plot_z_all_eorecpt + annotation_custom(ggplotGrob(ppca_load_z+theme_cowplot()+theme(axis.title=element_text(size=10),axis.text=element_text(size=8))), xmin = 6, xmax = 45, ymin = -21, ymax = -5)
ppca_inset_z_pc34_eorecpt<- ppca_plot_z_all_pc34_eorecpt + annotation_custom(ggplotGrob(ppca_load_z_pc34+theme_cowplot()+theme(axis.title=element_text(size=10),axis.text=element_text(size=8))), xmin = 6, xmax = 15.5, ymin = 0, ymax = 5.5)




###################################################################
#       gls model comparison on pca data (like Mull et al)        #
###################################################################

## need these for avg mass and tbv for each spp
avg_comb_data_sub<-subset(avg_comb_data, rownames(avg_comb_data) != 'Gymnarchus_niloticus')
sd_comb_data_sub<-subset(sd_comb_data, rownames(sd_comb_data) != 'Gymnarchus_niloticus')

## trees w out gymnarchus
no_gym_tree_sg<-drop.tip(all_tree_sg, "Gymnarchus_niloticus")


tree_sub<-no_gym_tree_sg 			## change this to represent which tree I want to use in the pPCA

## dataset is indiv_ppca2, which has the relevant info in a dataframe
sub_indiv_ppca<-subset(indiv_ppca2, indiv_ppca2$Species != 'Gymnarchus_niloticus')

## mean, stdev, se matrices for sub_indiv_ppca
Avg_PC1 <- ddply(sub_indiv_ppca, "Species", summarise,
			   N    = length(PC1),
               mean = mean(PC1),
               sd   = sd(PC1),
               se   = sd / sqrt(N))
if (any(is.na(Avg_PC1$se))) {
            Avg_PC1$se[which(is.na(Avg_PC1$se))] <- mean(Avg_PC1$se, na.rm = TRUE)}
Avg_PC2 <- ddply(sub_indiv_ppca, "Species", summarise,
			   N    = length(PC2),
               mean = mean(PC2),
               sd   = sd(PC2),
               se   = sd / sqrt(N))
if (any(is.na(Avg_PC2$se))) {
            Avg_PC2$se[which(is.na(Avg_PC2$se))] <- mean(Avg_PC2$se, na.rm = TRUE)}
Avg_PC3 <- ddply(sub_indiv_ppca, "Species", summarise,
			   N    = length(PC3),
               mean = mean(PC3),
               sd   = sd(PC3),
               se   = sd / sqrt(N))
if (any(is.na(Avg_PC3$se))) {
            Avg_PC3$se[which(is.na(Avg_PC3$se))] <- mean(Avg_PC3$se, na.rm = TRUE)}
Avg_PC4 <- ddply(sub_indiv_ppca, "Species", summarise,
			   N    = length(PC4),
               mean = mean(PC4),
               sd   = sd(PC4),
               se   = sd / sqrt(N))
if (any(is.na(Avg_PC4$se))) {
            Avg_PC4$se[which(is.na(Avg_PC4$se))] <- mean(Avg_PC4$se, na.rm = TRUE)}
Avg_PC5 <- ddply(sub_indiv_ppca, "Species", summarise,
			   N    = length(PC5),
               mean = mean(PC5),
               sd   = sd(PC5),
               se   = sd / sqrt(N))
if (any(is.na(Avg_PC5$se))) {
            Avg_PC5$se[which(is.na(Avg_PC5$se))] <- mean(Avg_PC5$se, na.rm = TRUE)}

avg_pc_data <- data.frame(Avg_PC1$mean, Avg_PC2$mean, Avg_PC3$mean, Avg_PC4$mean, Avg_PC5$mean)
colnames(avg_pc_data) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
rownames(avg_pc_data) <- Avg_PC1$Species

PC1_SE <- Avg_PC1$se
names(PC1_SE) <- Avg_PC1$Species
PC2_SE <- Avg_PC2$se
names(PC2_SE) <- Avg_PC2$Species
PC3_SE <- Avg_PC3$se
names(PC3_SE) <- Avg_PC3$Species
PC4_SE <- Avg_PC4$se
names(PC4_SE) <- Avg_PC4$Species
PC5_SE <- Avg_PC5$se
names(PC5_SE) <- Avg_PC5$Species


tb <- avg_comb_data_sub$TBV_M
bm <- avg_comb_data_sub$BM_M
#pc <- avg_pc_data[,1]	## change for which pc I am testing
#pc <- avg_pc_data[,2]
#pc <- avg_pc_data[,3]
pc <- avg_pc_data[,4]
#pc <- avg_pc_data[,5]

eo <- c("yes","yes","yes","yes","yes","no","no","no","no","yes","yes","yes","yes","no","yes","yes","no","no","yes","no","yes","no","yes","yes","yes","yes","yes","yes","yes","yes","no")
recpt <- c("tub","tub","tub","tub","tub","none","oamp","none","none","tub","tub","tub","tub","none","tub","tub","oamp","oamp","tub","none","tub","none","tub","tub","tub","tub","oamp","oamp","oamp","oamp","oamp")
data <- data.frame(avg_comb_data_sub$Species,tb,bm,pc,row.names=avg_comb_data_sub$Species)		## w phylogeny
colnames(data)<-c("Species","tb","bm","pc")

## for pc1
correlation=corPagel(.9,tree_sub)
pc.brbm <- gls(pc ~ bm, correlation = correlation, data=data, na.action = na.omit, method = "ML")
pc.eo <- gls(pc ~ bm + eo, correlation = correlation, data=data, na.action = na.omit, method = "ML")
pc.recpt <- gls(pc ~ bm + recpt, correlation = correlation, data=data, na.action = na.omit, method = "ML")
pc.both <- gls(pc ~ bm + eo + recpt, correlation = correlation, data=data, na.action = na.omit, method = "ML")

# AICc
AICc(pc.brbm)
AICc(pc.eo)
AICc(pc.recpt)
AICc(pc.both)

## for everyone else
correlation=corPagel(.9,tree_sub, form=~Species)
pc.brbm <- gls(pc ~ bm + tb, correlation = correlation, data=data, na.action = na.omit, method = "ML")
pc.eo <- gls(pc ~ bm + tb + eo, correlation = correlation, data=data, na.action = na.omit, method = "ML")
pc.recpt <- gls(pc ~ bm + tb + recpt, correlation = correlation, data=data, na.action = na.omit, method = "ML")
pc.both <- gls(pc ~ bm + tb + eo + recpt, correlation = correlation, data=data, na.action = na.omit, method = "ML")

# AICc
AICc(pc.brbm)
AICc(pc.eo)
AICc(pc.recpt)
AICc(pc.both)


###################################################
#   phylogenetic flexible discriminant analysis   #
###################################################

## note some of these libraries may already be loaded
#library(class)
#library(lattice)
#library(mda)
#library(nnet)
#source("phylo.fda.v0.2.R")
source('~/Desktop/Dropbox/Washington University, St. Louis/Carlson Lab/uCT/phylo.fda.v0.2.R')

tree<-all_tree_sg 

eo <- c("yes","yes","yes","yes","yes","no","no","no","no","yes","yes","yes","yes","yes","no","yes","yes","no","no","yes","no","yes","no","yes","yes","yes","yes","yes","yes","yes","yes","no")
recpt <- c("tub","tub","tub","tub","tub","none","oamp","none","none","tub","tub","tub","tub","tub","none","tub","tub","oamp","oamp","tub","none","tub","none","tub","tub","tub","tub","oamp","oamp","oamp","oamp","oamp")
eorecpt <- paste(eo,recpt)
egenrecpt <- c("tubegen","tubegen","tubegen","tubegen","tubegen","out","out","out","out","tubegen","tubegen","tubegen","tubegen","tubegen","out","tubegen","tubegen","out","out","tubegen","out","tubegen","out","tubegen","tubegen","tubegen","tubegen","ampegen","ampegen","ampegen","ampegen","out")
lin <- c("oto","oto","osteo","osteo","osteo","osteo","oto","oto","oto","oto","oto","oto","osteo","osteo","oto","oto","oto","oto","oto","osteo","osteo","osteo","oto","oto","oto","oto","oto","oto","oto","oto","oto","osteo")
egenrecptlin <- paste(egenrecpt,lin)

discrim_dat<-cbind(avg_region_data,as.numeric(as.factor(eo))-1,as.numeric(as.factor(recpt))-1,as.numeric(as.factor(eorecpt))-1,as.numeric(as.factor(egenrecpt))-1,as.numeric(as.factor(egenrecptlin))-1)		## recodes the text as a number
colnames(discrim_dat)<-c("OB","TEL","HB","OT","TS","CB","RoB","eo","recpt","eorecpt","egenrecpt","egenrecptlin")
discrim_dat<-as.data.frame(discrim_dat)

## Ordering data to match tip labels.
pfda_dat <- discrim_dat[,c(1:7,10)]
pfda_dat <- pfda_dat[match(tree$tip.label,rownames(pfda_dat)),]

## Defining groups and taxa.
gA <- pfda_dat$eorecpt # contains data on electrosensory phenotypes
names(gA) <- rownames(pfda_dat)
taxaA <- rownames(pfda_dat) # species names
#rownames(pfda_dat) <- taxaA # assigning species names to rows   

# Identifying optimal lambda: where is the strongest correlation between brain region data and electrosensory phenotype among taxa?
X <- pfda_dat[,1:7]	# variables we want from pfda_dat, i.e. just the region data
filename_stem <- "ecology" # A plot will appear in separate window; a pdf is saved in your working directory.
ol1 <- optLambda(X,gA,tree,idc=filename_stem )
ol1$optlambda # displaying the optimal lambda value

## optimal lambda value == 0, maybe this is bc of the convergent evolution

## pFDA, val = lambda value, tried 1, 0, and the fit lambda from the pca
#pfda_1 <- phylo.fda(X,gA,tree,val=1,treetrans=rescale)
pfda_0 <- phylo.fda(X,gA,tree,val=0,treetrans=rescale)
#pfda_pcl <- phylo.fda(X,gA,tree,val=avg_pca$lambda,treetrans=rescale)

#pfda_1_values <- predict(pfda_1,type="variates")
#pfda_1_class <- predict(pfda_1,type="class")
#pfda_1_prob <- predict(pfda_1,type="posterior")	# prob of each spp belonging to each class
#pfda_1_df <- as.data.frame(pfda_1_values)
#pfda_1_df <- cbind(rownames(pfda_dat),pfda_1_df,pfda_dat$eorecpt,pfda_1_class)
#colnames(pfda_1_df) <- c("Species","DA1","DA2","DA3","eorecpt","pred.class")

pfda_0_values <- predict(pfda_0,type="variates")
pfda_0_class <- predict(pfda_0,type="class")
pfda_0_prob <- predict(pfda_0,type="posterior")	# prob of each spp belonging to each class
pfda_0_df <- as.data.frame(pfda_0_values)
pfda_0_df <- cbind(rownames(pfda_dat),pfda_0_df,pfda_dat$eorecpt,pfda_0_class)
colnames(pfda_0_df) <- c("Species","DA1","DA2","DA3","eorecpt","pred.class")

#pfda_pcl_values <- predict(pfda_pcl,type="variates")
#pfda_pcl_class <- predict(pfda_pcl,type="class")
#pfda_pcl_prob <- predict(pfda_pcl,type="posterior")	# prob of each spp belonging to each class
#pfda_pcl_df <- as.data.frame(pfda_pcl_values)
#pfda_pcl_df <- cbind(rownames(pfda_dat),pfda_pcl_df,pfda_dat$eorecpt,pfda_pcl_class)
#colnames(pfda_pcl_df) <- c("Species","DA1","DA2","DA3","eorecpt","pred.class")

## plotting pfda results
#pfda_1_plot <- ggplot(pfda_1_df, aes(x=DA1, y=DA2, shape=Species, fill=Species, color=Species, group=eorecpt)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
#	new_scale_color() +
#	stat_chull(fill=NA, aes(lty=eorecpt, colour=eorecpt)) + scale_linetype_manual(values=c(2,2,3,1)) + scale_colour_manual(values=c("grey","black","black","black"))
#pfda_1_plot_ld23 <- ggplot(pfda_1_df, aes(x=DA2, y=DA3, shape=Species, fill=Species, color=Species, group=eorecpt)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
#	new_scale_color() +
#	stat_chull(fill=NA, aes(lty=eorecpt, colour=eorecpt)) + scale_linetype_manual(values=c(2,2,3,1)) + scale_colour_manual(values=c("grey","black","black","black"))

pfda_0_plot <- ggplot(pfda_0_df, aes(x=DA1, y=DA2, shape=Species, fill=Species, color=Species, group=eorecpt)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	new_scale_color() +
	stat_chull(fill=NA, aes(lty=eorecpt, colour=eorecpt)) + scale_linetype_manual(values=c(2,2,3,1)) + scale_colour_manual(values=c("grey","black","black","black"))
pfda_0_plot_ld23 <- ggplot(pfda_0_df, aes(x=DA2, y=DA3, shape=Species, fill=Species, color=Species, group=eorecpt)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	new_scale_color() +
	stat_chull(fill=NA, aes(lty=eorecpt, colour=eorecpt)) + scale_linetype_manual(values=c(2,2,3,1)) + scale_colour_manual(values=c("grey","black","black","black"))

#pfda_pcl_plot <- ggplot(pfda_pcl_df, aes(x=DA1, y=DA2, shape=Species, fill=Species, color=Species, group=eorecpt)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
#	new_scale_color() +
#	stat_chull(fill=NA, aes(lty=eorecpt, colour=eorecpt)) + scale_linetype_manual(values=c(2,2,3,1)) + scale_colour_manual(values=c("grey","black","black","black"))
#pfda_pcl_plot_ld23 <- ggplot(pfda_pcl_df, aes(x=DA2, y=DA3, shape=Species, fill=Species, color=Species, group=eorecpt)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="none") + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
#	new_scale_color() +
#	stat_chull(fill=NA, aes(lty=eorecpt, colour=eorecpt)) + scale_linetype_manual(values=c(2,2,3,1)) + scale_colour_manual(values=c("grey","black","black","black"))


#####################################################################
#           removing G niloticus which lacks mass data              #
#####################################################################

sub_comb_log_data<-subset(comb_log_data, Species !="Gymnarchus_niloticus")
morm_mass_plot<-ggplot(sub_comb_log_data, aes(x=Mass.g, y=Total.Brain.Volume, shape=Species, fill=Species, color=Species, group=ID)) + geom_line() + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,23,23,24,21,23,7,22,9,24,25,24,25,7,24,4,25,7,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(y=bquote('Log Total Brain Volume' ~(mm^3)), x = bquote('Log Body Mass ' (g))) + scale_fill_manual(values = c("Blue", "seagreen2", "mediumorchid1", "mediumorchid1", "mediumorchid1", "white", "gray", "Black", "Black", "Blue", "Blue", "seagreen2", "mediumorchid1", "Black", "seagreen2", "seagreen2", "gray", "gray", "mediumorchid1", "white", "mediumorchid1", "Black", "seagreen2", "Blue", "Blue", "Blue", "gray", "gray", "gray", "gray", "white")) + scale_colour_manual(values = c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "mediumorchid1", "black", "mediumorchid1", "black", "black", "black", "black", "Blue", "black", "gray", "black", "gray", "black"));

tree_sub<-no_gym_tree_sg 			## change this to represent which tree I want to use

## linear regression of brain vs body w phylogeny
sub_total <- sub_comb_log_data$Total.Brain.Volume
names(sub_total) <- sub_comb_log_data$Species
sub_mass <- sub_comb_log_data$Mass.g
names(sub_mass) <- sub_comb_log_data$Species

## Ives for just y
TBV_suberror <- ddply(sub_comb_log_data, "Species", summarise,
			   N    = length(Total.Brain.Volume),
               mean = mean(Total.Brain.Volume),
               sd   = sd(Total.Brain.Volume),
               se   = sd / sqrt(N))
if (any(is.na(TBV_suberror$se))) {
            TBV_suberror$se[which(is.na(TBV_suberror$se))] <- mean(TBV_suberror$se, na.rm = TRUE)}
Mass_suberror <- ddply(sub_comb_log_data, "Species", summarise,
			   N    = length(Mass.g),
               mean = mean(Mass.g),
               sd   = sd(Mass.g),
               se   = sd / sqrt(N))
if (any(is.na(Mass_suberror$se))) {
            Mass_suberror$se[which(is.na(Mass_suberror$se))] <- mean(Mass_suberror$se, na.rm = TRUE)}
avg_mass_data<-data.frame(TBV_suberror$Species,TBV_suberror$mean,Mass_suberror$mean)
colnames(avg_mass_data)<-c("Species","TBV","Mass")
rownames(avg_mass_data)<-avg_mass_data$Species

TBVsub_SE <- TBV_suberror$se
names(TBVsub_SE) <- TBV_suberror$Species

text <- paste("Pulse mormyroids",
              "Pulse gymnotiforms",
              "Wave gymnotiforms",
              "Synodontis siluriforms",
             "Non-electric siluriforms",
             "Outgroup osteoglossiforms",
             "Outgroup otophysans", sep = " ")
text <- paste("Pulse mormyroids            ",
              "Non-electric siluriforms ",
              "Pulse gymnotiforms         ",
              "Outgroup osteoglossiforms",
              "Wave gymnotiforms         ",
              "Outgroup otophysans      ",
              "Synodontis siluriforms   ")

text.p <- ggparagraph(text = text, size = 9, color = "black")

#################################################################
#    brain body allometry evolution across teleosts analysis    #
#################################################################

## Combines data from Tsuboi et al 2018, Sukhum et al 2016 (reffed in Tsuboi 2021), and this study
## Brain mass data from this study was calculated via a volume conversion

tel_data <- read.csv("eLife Revisions/SchumacherAndCarlson_SupplementaryFile2.csv")
tel_data[,11:12] <- log10(tel_data[,11:12])
tel_tree_rab<-drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, tel_data$Genus_Species))  #tree with all study species

phylo_tel_data <- subset(tel_data, Genus_Species %in% tel_tree_rab$tip.label)		## 3673/4068 individs and 881/1017 spps remain

## comp raw, log, ln dat
#raw_plot<-ggplot(raw_tel_data, aes(x=Body.weight..g., y=Brain.weight..g.)) + geom_point() + theme_classic()
#log10_plot<-ggplot(tel_data, aes(x=Body.weight..g., y=Brain.weight..g.)) + geom_point() + theme_classic()
#ln_plot<-ggplot(ln_tel_data, aes(x=Body.weight..g., y=Brain.weight..g.)) + geom_point() + theme_classic()


## First, is my conversion okay to use with Tsuboi and Sukhum et al's data?
## Ancova comparing brain allometries from my data to their data for spp with 3+ individs in both groups
## Qualifying spps: X nigri, C ornata, G petersii, B brachyistius, B niger, P buchholzi
xn_dat <- subset(phylo_tel_data, Genus_Species == "Xenomystus_nigri")
co_dat <- subset(phylo_tel_data, Genus_Species == "Chitala_ornata")
gp_dat <- subset(phylo_tel_data, Genus_Species == "Gnathonemus_petersii")
bb_dat <- subset(phylo_tel_data, Genus_Species == "Brienomyrus_brachyistius")
bn_dat <- subset(phylo_tel_data, Genus_Species == "Brevimyrus_niger")
pb_dat <- subset(tel_data, Genus_Species == "Pantodon_buchholzi")
mt_dat <- subset(tel_data, Genus_Species == "Mormyrus_tapirus")
gc_dat <- subset(tel_data, Genus_Species == "Gymnotus_carapo")
aa_dat <- subset(tel_data, Genus_Species == "Apteronotus_albifrons")
so_dat <- subset(tel_data, Genus_Species == "Sternarchella_orthos")
leg_dat <- rbind(gp_dat,gc_dat)
overlap_dat <- rbind(xn_dat,co_dat,gp_dat,bb_dat,bn_dat,pb_dat,mt_dat,gc_dat,aa_dat,so_dat)

xn.gls <- gls(Brain.weight..g. ~ Body.weight..g. + Method + Body.weight..g.:Method, data=xn_dat, na.action = na.omit)
xn.gls2 <- gls(Brain.weight..g. ~ Body.weight..g. + Method, data=xn_dat, na.action = na.omit)
## Ancovas
## anova = type I, where the order of the variables in the model matters w/ first variable explaining the most variation
anova(xn.gls) #p value of x:type shows slope
anova(xn.gls2)		## xn is still sig diff
xn.mass <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(xn_dat, Method == "Mass"))
xn.vol <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(xn_dat, Method == "Volume Conversion"))
xn.all <- gls(Brain.weight..g. ~ Body.weight..g., data=xn_dat)

co.gls <- gls(Brain.weight..g. ~ Body.weight..g. + Method + Body.weight..g.:Method, data=co_dat, na.action = na.omit)
co.gls2 <- gls(Brain.weight..g. ~ Body.weight..g. + Method, data=co_dat, na.action = na.omit)
## Ancovas
## anova = type I, where the order of the variables in the model matters w/ first variable explaining the most variation
anova(co.gls) #p value of x:type shows slope
anova(co.gls2)
co.mass <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(co_dat, Method == "Mass"))
co.vol <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(co_dat, Method == "Volume Conversion"))
co.all <- gls(Brain.weight..g. ~ Body.weight..g., data=co_dat)

gp.gls <- gls(Brain.weight..g. ~ Body.weight..g. + Method + Body.weight..g.:Method, data=gp_dat, na.action = na.omit)
gp.gls2 <- gls(Brain.weight..g. ~ Body.weight..g. + Method, data=gp_dat, na.action = na.omit)
## Ancovas
## anova = type I, where the order of the variables in the model matters w/ first variable explaining the most variation
anova(gp.gls) #p value of x:type shows slope
anova(gp.gls2)		## gp is also still sig diff
gp.mass <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(gp_dat, Method == "Mass"))
gp.vol <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(gp_dat, Method == "Volume Conversion"))
gp.all <- gls(Brain.weight..g. ~ Body.weight..g., data=gp_dat)

bb.gls <- gls(Brain.weight..g. ~ Body.weight..g. + Method + Body.weight..g.:Method, data=bb_dat, na.action = na.omit)
bb.gls2 <- gls(Brain.weight..g. ~ Body.weight..g. + Method, data=bb_dat, na.action = na.omit)
## Ancovas
## anova = type I, where the order of the variables in the model matters w/ first variable explaining the most variation
anova(bb.gls) #p value of x:type shows slope
anova(bb.gls2)
bb.mass <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(bb_dat, Method == "Mass"))
bb.vol <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(bb_dat, Method == "Volume Conversion"))
bb.all <- gls(Brain.weight..g. ~ Body.weight..g., data=bb_dat)

bn.gls <- gls(Brain.weight..g. ~ Body.weight..g. + Method + Body.weight..g.:Method, data=bn_dat, na.action = na.omit)
bn.gls2 <- gls(Brain.weight..g. ~ Body.weight..g. + Method, data=bn_dat, na.action = na.omit)
## Ancovas
## anova = type I, where the order of the variables in the model matters w/ first variable explaining the most variation
anova(bn.gls) #p value of x:type shows slope
anova(bn.gls2)
bn.mass <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(bn_dat, Method == "Mass"))
bn.vol <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(bn_dat, Method == "Volume Conversion"))
bn.all <- gls(Brain.weight..g. ~ Body.weight..g., data=bn_dat)

pb.gls <- gls(Brain.weight..g. ~ Body.weight..g. + Method + Body.weight..g.:Method, data=pb_dat, na.action = na.omit)
pb.gls2 <- gls(Brain.weight..g. ~ Body.weight..g. + Method, data=pb_dat, na.action = na.omit)
## Ancovas
## anova = type I, where the order of the variables in the model matters w/ first variable explaining the most variation
anova(pb.gls) #p value of x:type shows slope
anova(pb.gls2)		## pb is still sig diff
pb.mass <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(pb_dat, Method == "Mass"))
pb.vol <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(pb_dat, Method == "Volume Conversion"))
pb.all <- gls(Brain.weight..g. ~ Body.weight..g., data=pb_dat)

xn.plot <- ggplot(xn_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="X nigri", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) +
	geom_abline(slope = xn.mass$coefficients[2], intercept = xn.mass$coefficients[1],color="magenta") + geom_abline(slope = xn.vol$coefficients[2], intercept = xn.vol$coefficients[1],color="cyan")
co.plot <- ggplot(co_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="C ornata", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) +
	geom_abline(slope = co.mass$coefficients[2], intercept = co.mass$coefficients[1],color="magenta") + geom_abline(slope = co.vol$coefficients[2], intercept = co.vol$coefficients[1],color="cyan")
gp.plot <- ggplot(gp_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="G petersii", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) +
	geom_abline(slope = gp.mass$coefficients[2], intercept = gp.mass$coefficients[1],color="magenta") + geom_abline(slope = gp.vol$coefficients[2], intercept = gp.vol$coefficients[1],color="cyan")
bb.plot <- ggplot(bb_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="B brachyistius", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) +
	geom_abline(slope = bb.mass$coefficients[2], intercept = bb.mass$coefficients[1],color="magenta") + geom_abline(slope = bb.vol$coefficients[2], intercept = bb.vol$coefficients[1],color="cyan")
bn.plot <- ggplot(bn_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="B niger", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) +
	geom_abline(slope = bn.mass$coefficients[2], intercept = bn.mass$coefficients[1],color="magenta") + geom_abline(slope = bn.vol$coefficients[2], intercept = bn.vol$coefficients[1],color="cyan")
pb.plot <- ggplot(pb_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="P buchholzi", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) +
	geom_abline(slope = pb.mass$coefficients[2], intercept = pb.mass$coefficients[1],color="magenta") + geom_abline(slope = pb.vol$coefficients[2], intercept = pb.vol$coefficients[1],color="cyan")
aa.plot <- ggplot(subset(phylo_tel_data, Genus_Species == "Apteronotus_albifrons"), aes(x=Body.weight..g., y=Brain.weight..g., color=Method)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="A albifrons", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan"))
gc.plot <- ggplot(subset(phylo_tel_data, Genus_Species == "Gymnotus_carapo"), aes(x=Body.weight..g., y=Brain.weight..g., color=Method)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="G carapo", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan"))
mt.plot <- ggplot(subset(phylo_tel_data, Genus_Species == "Mormyrus_tapirus"), aes(x=Body.weight..g., y=Brain.weight..g., color=Method)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="M tapirus", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan"))

vol.mass.plots <- ggarrange(xn.plot,co.plot,gp.plot,bb.plot,bn.plot,pb.plot,mt.plot,aa.plot,gc.plot, nrow=3, ncol=3, common.legend = TRUE)

## plots w study indicated by shape, also gls across all individs in black
aa.all <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(phylo_tel_data, Genus_Species == "Apteronotus_albifrons"))
gc.all <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(phylo_tel_data, Genus_Species == "Gymnotus_carapo"))
mt.all <- gls(Brain.weight..g. ~ Body.weight..g., data=subset(phylo_tel_data, Genus_Species == "Mormyrus_tapirus"))

ref_shapes <- c("Sukhum et al (2016) Proc.B. 283: 20162157"=1, "This study"=2, "Kaufman JA (2003) Curr. Anthropol. 44:705-707."=3, "Nilsson GE (1996) J. Exp. Biol. 199:603-607"=4, "Albert JS et al. (1999) Proceedings of the 5th Indo-Pacific Fish Conference (Seret B, Sire J-Y, eds), pp 647-656"=5, "Bauchot R, Diagne M, Platel R, Ridet J-M, Bauchot M-L. in Fishbase. www.fishbase.org. date of access: Apr.-13-17"=6)

xn.plot <- ggplot(xn_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method, shape=Reference)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="X nigri", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) + scale_shape_manual(values=ref_shapes) +
	geom_abline(slope = xn.mass$coefficients[2], intercept = xn.mass$coefficients[1],color="magenta") + geom_abline(slope = xn.vol$coefficients[2], intercept = xn.vol$coefficients[1],color="cyan") + geom_abline(slope = xn.all$coefficients[2], intercept = xn.all$coefficients[1],color="black")
co.plot <- ggplot(co_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method, shape=Reference)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="C ornata", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) + scale_shape_manual(values=ref_shapes) +
	geom_abline(slope = co.mass$coefficients[2], intercept = co.mass$coefficients[1],color="magenta") + geom_abline(slope = co.vol$coefficients[2], intercept = co.vol$coefficients[1],color="cyan") + geom_abline(slope = co.all$coefficients[2], intercept = co.all$coefficients[1],color="black")
gp.plot <- ggplot(gp_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method, shape=Reference)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="G petersii", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) + scale_shape_manual(values=ref_shapes) +
	geom_abline(slope = gp.mass$coefficients[2], intercept = gp.mass$coefficients[1],color="magenta") + geom_abline(slope = gp.vol$coefficients[2], intercept = gp.vol$coefficients[1],color="cyan") + geom_abline(slope = gp.all$coefficients[2], intercept = gp.all$coefficients[1],color="black")
bb.plot <- ggplot(bb_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method, shape=Reference)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="B brachyistius", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) + scale_shape_manual(values=ref_shapes) +
	geom_abline(slope = bb.mass$coefficients[2], intercept = bb.mass$coefficients[1],color="magenta") + geom_abline(slope = bb.vol$coefficients[2], intercept = bb.vol$coefficients[1],color="cyan") + geom_abline(slope = bb.all$coefficients[2], intercept = bb.all$coefficients[1],color="black")
bn.plot <- ggplot(bn_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method, shape=Reference)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="B niger", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) + scale_shape_manual(values=ref_shapes) +
	geom_abline(slope = bn.mass$coefficients[2], intercept = bn.mass$coefficients[1],color="magenta") + geom_abline(slope = bn.vol$coefficients[2], intercept = bn.vol$coefficients[1],color="cyan") + geom_abline(slope = bn.all$coefficients[2], intercept = bn.all$coefficients[1],color="black")
pb.plot <- ggplot(pb_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method, shape=Reference)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="P buchholzi", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) + scale_shape_manual(values=ref_shapes) +
	geom_abline(slope = pb.mass$coefficients[2], intercept = pb.mass$coefficients[1],color="magenta") + geom_abline(slope = pb.vol$coefficients[2], intercept = pb.vol$coefficients[1],color="cyan") + geom_abline(slope = pb.all$coefficients[2], intercept = pb.all$coefficients[1],color="black")
aa.plot <- ggplot(subset(phylo_tel_data, Genus_Species == "Apteronotus_albifrons"), aes(x=Body.weight..g., y=Brain.weight..g., color=Method, shape=Reference)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="A albifrons", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) + scale_shape_manual(values=ref_shapes) +
	geom_abline(slope = aa.all$coefficients[2], intercept = aa.all$coefficients[1],color="black")
gc.plot <- ggplot(subset(phylo_tel_data, Genus_Species == "Gymnotus_carapo"), aes(x=Body.weight..g., y=Brain.weight..g., color=Method, shape=Reference)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="G carapo", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) + scale_shape_manual(values=ref_shapes) +
	geom_abline(slope = gc.all$coefficients[2], intercept = gc.all$coefficients[1],color="black")
mt.plot <- ggplot(subset(phylo_tel_data, Genus_Species == "Mormyrus_tapirus"), aes(x=Body.weight..g., y=Brain.weight..g., color=Method, shape=Reference)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="M tapirus", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) + scale_shape_manual(values=ref_shapes) +
	geom_abline(slope = mt.all$coefficients[2], intercept = mt.all$coefficients[1],color="black")
for.legend <- ggplot(leg_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method, shape=Reference)) + geom_point() + theme_classic() + theme(legend.position="right") + labs(title="M tapirus", x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) + scale_shape_manual(values=ref_shapes)
leg <- get_legend(for.legend)
as_ggplot(leg)

vol.mass.plots <- ggarrange(xn.plot,co.plot,gp.plot,bb.plot,bn.plot,pb.plot,mt.plot,aa.plot,gc.plot, nrow=3, ncol=3, common.legend = TRUE, legend = "none")

## plotting all of the overlap data
overlap.plot <- ggplot(overlap_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Method, shape=Genus_Species)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_color_manual(values=c("magenta","cyan")) + scale_shape_manual(values=c(1:10))


## plotting all data from osteo and oto
ost_bm_dat <- subset(phylo_tel_data, Order == "Osteoglossiformes")
gym_bm_dat <- subset(phylo_tel_data, Order == "Gymnotiformes")
sil_bm_dat <- subset(phylo_tel_data, Order == "Siluriformes")
cha_bm_dat <- subset(phylo_tel_data, Order == "Characiformes")
cyp_bm_dat <- subset(phylo_tel_data, Order == "Cypriniformes")
rel_bm_dat <- rbind(ost_bm_dat, gym_bm_dat, sil_bm_dat, cha_bm_dat, cyp_bm_dat)

rel.bm.plot <- ggplot(rel_bm_dat, aes(x=Body.weight..g., y=Brain.weight..g., color=Genus_Species, shape=Order)) + geom_point() + theme_classic() + theme(legend.position="right") + labs(x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) 
leg2 <- get_legend(rel.bm.plot)
as_ggplot(leg2)


###################################
#   bayesian allometry modeling   #
###################################

## see README and other files for code for this analysis


###################################################
#						   						  #
#      GLS Analysis and Ancovas, w phylogeny      #
#												  #
###################################################

##############################################
#    pgls regression lines for each group    #
##############################################

## first get std error for each species and region, use mean(std err) for spp with one individual
OB_error <- ddply(comb_log_data, "Species", summarise,
			   N    = length(OB),
               mean = mean(OB),
               sd   = sd(OB),
               se   = sd / sqrt(N))
if (any(is.na(OB_error$se))) {
            OB_error$se[which(is.na(OB_error$se))] <- mean(OB_error$se, na.rm = TRUE)}
TEL_error <- ddply(comb_log_data, "Species", summarise,
			   N    = length(TEL),
               mean = mean(TEL),
               sd   = sd(TEL),
               se   = sd / sqrt(N))
if (any(is.na(TEL_error$se))) {
            TEL_error$se[which(is.na(TEL_error$se))] <- mean(TEL_error$se, na.rm = TRUE)}
HB_error <- ddply(comb_log_data, "Species", summarise,
			   N    = length(HB),
               mean = mean(HB),
               sd   = sd(HB),
               se   = sd / sqrt(N))
if (any(is.na(HB_error$se))) {
            HB_error$se[which(is.na(HB_error$se))] <- mean(HB_error$se, na.rm = TRUE)}
OT_error <- ddply(comb_log_data, "Species", summarise,
			   N    = length(OT),
               mean = mean(OT),
               sd   = sd(OT),
               se   = sd / sqrt(N))
if (any(is.na(OT_error$se))) {
            OT_error$se[which(is.na(OT_error$se))] <- mean(OT_error$se, na.rm = TRUE)}
TS_error <- ddply(comb_log_data, "Species", summarise,
			   N    = length(TS),
               mean = mean(TS),
               sd   = sd(TS),
               se   = sd / sqrt(N))
if (any(is.na(TS_error$se))) {
            TS_error$se[which(is.na(TS_error$se))] <- mean(TS_error$se, na.rm = TRUE)}
CB_error <- ddply(comb_log_data, "Species", summarise,
			   N    = length(CB),
               mean = mean(CB),
               sd   = sd(CB),
               se   = sd / sqrt(N))
if (any(is.na(CB_error$se))) {
            CB_error$se[which(is.na(CB_error$se))] <- mean(CB_error$se, na.rm = TRUE)}
ROB_error <- ddply(comb_log_data, "Species", summarise,
			   N    = length(RoB),
               mean = mean(RoB),
               sd   = sd(RoB),
               se   = sd / sqrt(N))
if (any(is.na(ROB_error$se))) {
            ROB_error$se[which(is.na(ROB_error$se))] <- mean(ROB_error$se, na.rm = TRUE)}

Total.Brain.Volume_error <- ddply(comb_log_data, "Species", summarise,
			   N    = length(Total.Brain.Volume),
               mean = mean(Total.Brain.Volume),
               sd   = sd(Total.Brain.Volume),
               se   = sd / sqrt(N))
if (any(is.na(Total.Brain.Volume_error$se))) {
            Total.Brain.Volume_error$se[which(is.na(Total.Brain.Volume_error$se))] <- mean(Total.Brain.Volume_error$se, na.rm = TRUE)}

avg_log_data <- data.frame(OB_error$Species, OB_error$mean, TEL_error$mean, HB_error$mean, OT_error$mean, TS_error$mean, CB_error$mean, ROB_error$mean, Total.Brain.Volume_error$mean)
colnames(avg_log_data) <- c("Species", "OB", "TEL", "HB", "OT", "TS", "CB", "ROB", "TBV")
rownames(avg_log_data) <- avg_log_data$Species


tree<-all_tree_sg 			## change this to represent which tree I want to use in the pgls

tubegen.tree<-drop.tip(tree, setdiff(tree$tip.label, tubegen_names))
ampegen.tree<-drop.tip(tree, setdiff(tree$tip.label, ampegen_names))
out.tree<-drop.tip(tree, setdiff(tree$tip.label, out_names))
gym.tree<-drop.tip(tree, setdiff(tree$tip.label, gym_names))
morm.tree<-drop.tip(tree, setdiff(tree$tip.label, morm_names))
allgout.tree<-drop.tip(tree, setdiff(tree$tip.label, allgout_names))
allmout.tree<-drop.tip(tree, setdiff(tree$tip.label, allmout_names))
gwave.tree<-drop.tip(tree, setdiff(tree$tip.label, gwave_names))
gpulse.tree<-drop.tip(tree, setdiff(tree$tip.label, gpulse_names))

## subset for eo + recept type
tubegen_avg<-subset(avg_log_data, Species %in% tubegen_names)
ampegen_avg<-subset(avg_log_data, Species %in% ampegen_names)
out_avg<-subset(avg_log_data, Species %in% out_names)
gym_avg<-subset(avg_log_data, Species %in% gym_names)
morm_avg<-subset(avg_log_data, Species %in% morm_names)
allgout_avg<-subset(avg_log_data, Species %in% allgout_names)
allmout_avg<-subset(avg_log_data, Species %in% allmout_names)
gwave_avg<-subset(avg_log_data, Species %in% gwave_names)
gpulse_avg<-subset(avg_log_data, Species %in% gpulse_names)

hyp_egenrecpt <- c("tubegen","tubegen","tubegen","tubegen","tubegen","out","out","out","out","tubegen","tubegen","tubegen","tubegen","tubegen","out","tubegen","tubegen","out","out","tubegen","out","tubegen","out","tubegen","tubegen","tubegen","tubegen","ampegen","ampegen","ampegen","ampegen","out")
hyp_eolin <- c("gym","gym","morm","morm","morm","mout","gout","gout","gout","gym","gym","gym","morm","morm","gout","gym","gym","gout","gout","morm","mout","morm","gout","gym","gym","gym","gym","syno","syno","syno","syno","mout")


## pgls lines
OB_all<-gls(OB~TBV, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
OB_egenrecpt<-gls(OB~TBV+hyp_egenrecpt, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
OB_eolin<-gls(OB~TBV+hyp_eolin, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
OB_lam2<-as.numeric(OB_eolin$modelStruct)
OB_lam<-as.numeric(OB_egenrecpt$modelStruct)
if (OB_lam > 1)
	OB_lam = 1
if (OB_lam < 0)
	OB_lam = 0
if (OB_lam2 > 1)
	OB_lam2 = 1
if (OB_lam2 < 0)
	OB_lam2 = 0
OB_tubegen<-gls(OB~TBV, data=tubegen_avg, correlation=corPagel(OB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OB_ampegen<-gls(OB~TBV, data=ampegen_avg, correlation=corPagel(OB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OB_out<-gls(OB~TBV, data=out_avg, correlation=corPagel(OB_lam,phy=out.tree,form=~Species, fixed=TRUE))
OB_gym<-gls(OB~TBV, data=gym_avg, correlation=corPagel(OB_lam2,phy=gym.tree,form=~Species, fixed=TRUE))
OB_morm<-gls(OB~TBV, data=morm_avg, correlation=corPagel(OB_lam2,phy=morm.tree,form=~Species, fixed=TRUE))
OB_allgout<-gls(OB~TBV, data=allgout_avg, correlation=corPagel(OB_lam2,phy=allgout.tree,form=~Species, fixed=TRUE))
OB_allmout<-gls(OB~TBV, data=allmout_avg, correlation=corPagel(OB_lam2,phy=allmout.tree,form=~Species, fixed=TRUE))
OB_gwave<-gls(OB~TBV, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OB_gpulse<-gls(OB~TBV, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

TEL_all<-gls(TEL~TBV, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
TEL_egenrecpt<-gls(TEL~TBV+hyp_egenrecpt, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
TEL_eolin<-gls(TEL~TBV+hyp_eolin, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
TEL_lam2<-as.numeric(TEL_eolin$modelStruct)
TEL_lam<-as.numeric(TEL_egenrecpt$modelStruct)
if (TEL_lam > 1)
	TEL_lam = 1
if (TEL_lam < 0)
	TEL_lam = 0
if (TEL_lam2 > 1)
	TEL_lam2 = 1
if (TEL_lam2 < 0)
	TEL_lam2 = 0
TEL_tubegen<-gls(TEL~TBV, data=tubegen_avg, correlation=corPagel(TEL_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TEL_ampegen<-gls(TEL~TBV, data=ampegen_avg, correlation=corPagel(TEL_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TEL_out<-gls(TEL~TBV, data=out_avg, correlation=corPagel(TEL_lam,phy=out.tree,form=~Species, fixed=TRUE))
TEL_gym<-gls(TEL~TBV, data=gym_avg, correlation=corPagel(TEL_lam2,phy=gym.tree,form=~Species, fixed=TRUE))
TEL_morm<-gls(TEL~TBV, data=morm_avg, correlation=corPagel(TEL_lam2,phy=morm.tree,form=~Species, fixed=TRUE))
TEL_allgout<-gls(TEL~TBV, data=allgout_avg, correlation=corPagel(TEL_lam2,phy=allgout.tree,form=~Species, fixed=TRUE))
TEL_allmout<-gls(TEL~TBV, data=allmout_avg, correlation=corPagel(TEL_lam2,phy=allmout.tree,form=~Species, fixed=TRUE))
TEL_gwave<-gls(TEL~TBV, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
TEL_gpulse<-gls(TEL~TBV, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

HB_all<-gls(HB~TBV, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
HB_egenrecpt<-gls(HB~TBV+hyp_egenrecpt, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
HB_eolin<-gls(HB~TBV+hyp_eolin, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
HB_lam2<-as.numeric(HB_eolin$modelStruct)
HB_lam<-as.numeric(HB_egenrecpt$modelStruct)
if (HB_lam > 1)
	HB_lam = 1
if (HB_lam < 0)
	HB_lam = 0
if (HB_lam2 > 1)
	HB_lam2 = 1
if (HB_lam2 < 0)
	HB_lam2 = 0
HB_tubegen<-gls(HB~TBV, data=tubegen_avg, correlation=corPagel(HB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
HB_ampegen<-gls(HB~TBV, data=ampegen_avg, correlation=corPagel(HB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
HB_out<-gls(HB~TBV, data=out_avg, correlation=corPagel(HB_lam,phy=out.tree,form=~Species, fixed=TRUE))
HB_gym<-gls(HB~TBV, data=gym_avg, correlation=corPagel(HB_lam2,phy=gym.tree,form=~Species, fixed=TRUE))
HB_morm<-gls(HB~TBV, data=morm_avg, correlation=corPagel(HB_lam2,phy=morm.tree,form=~Species, fixed=TRUE))
HB_allgout<-gls(HB~TBV, data=allgout_avg, correlation=corPagel(HB_lam2,phy=allgout.tree,form=~Species, fixed=TRUE))
HB_allmout<-gls(HB~TBV, data=allmout_avg, correlation=corPagel(HB_lam2,phy=allmout.tree,form=~Species, fixed=TRUE))
HB_gwave<-gls(HB~TBV, data=gwave_avg, correlation=corPagel(0.8756104,phy=gwave.tree,form=~Species, fixed=TRUE))
HB_gpulse<-gls(HB~TBV, data=gpulse_avg, correlation=corPagel(0.8756104,phy=gpulse.tree,form=~Species, fixed=TRUE))

OT_all<-gls(OT~TBV, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
OT_egenrecpt<-gls(OT~TBV+hyp_egenrecpt, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
OT_eolin<-gls(OT~TBV+hyp_eolin, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
OT_lam2<-as.numeric(OT_eolin$modelStruct)
OT_lam<-as.numeric(OT_egenrecpt$modelStruct)
if (OT_lam > 1)
	OT_lam = 1
if (OT_lam < 0)
	OT_lam = 0
if (OT_lam2 > 1)
	OT_lam2 = 1
if (OT_lam2 < 0)
	OT_lam2 = 0
OT_tubegen<-gls(OT~TBV, data=tubegen_avg, correlation=corPagel(OT_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OT_ampegen<-gls(OT~TBV, data=ampegen_avg, correlation=corPagel(OT_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OT_out<-gls(OT~TBV, data=out_avg, correlation=corPagel(OT_lam,phy=out.tree,form=~Species, fixed=TRUE))
OT_gym<-gls(OT~TBV, data=gym_avg, correlation=corPagel(OT_lam2,phy=gym.tree,form=~Species, fixed=TRUE))
OT_morm<-gls(OT~TBV, data=morm_avg, correlation=corPagel(OT_lam2,phy=morm.tree,form=~Species, fixed=TRUE))
OT_allgout<-gls(OT~TBV, data=allgout_avg, correlation=corPagel(OT_lam2,phy=allgout.tree,form=~Species, fixed=TRUE))
OT_allmout<-gls(OT~TBV, data=allmout_avg, correlation=corPagel(OT_lam2,phy=allmout.tree,form=~Species, fixed=TRUE))
OT_gwave<-gls(OT~TBV, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OT_gpulse<-gls(OT~TBV, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

TS_all<-gls(TS~TBV, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
TS_egenrecpt<-gls(TS~TBV+hyp_egenrecpt, data=avg_log_data, correlation=corPagel(.5,phy=tree,form=~Species))
TS_eolin<-gls(TS~TBV+hyp_eolin, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
TS_lam2<-as.numeric(TS_eolin$modelStruct)
TS_lam<-as.numeric(TS_egenrecpt$modelStruct)
if (TS_lam > 1)
	TS_lam = 1
if (TS_lam < 0)
	TS_lam = 0
if (TS_lam2 > 1)
	TS_lam2 = 1
if (TS_lam2 < 0)
	TS_lam2 = 0
TS_tubegen<-gls(TS~TBV, data=tubegen_avg, correlation=corPagel(TS_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TS_ampegen<-gls(TS~TBV, data=ampegen_avg, correlation=corPagel(TS_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TS_out<-gls(TS~TBV, data=out_avg, correlation=corPagel(TS_lam,phy=out.tree,form=~Species, fixed=TRUE))
TS_gym<-gls(TS~TBV, data=gym_avg, correlation=corPagel(TS_lam2,phy=gym.tree,form=~Species, fixed=TRUE))
TS_morm<-gls(TS~TBV, data=morm_avg, correlation=corPagel(TS_lam2,phy=morm.tree,form=~Species, fixed=TRUE))
TS_allgout<-gls(TS~TBV, data=allgout_avg, correlation=corPagel(TS_lam2,phy=allgout.tree,form=~Species, fixed=TRUE))
TS_allmout<-gls(TS~TBV, data=allmout_avg, correlation=corPagel(TS_lam2,phy=allmout.tree,form=~Species, fixed=TRUE))
TS_gwave<-gls(TS~TBV, data=gwave_avg, correlation=corPagel(0.6012342,phy=gwave.tree,form=~Species, fixed=TRUE))
TS_gpulse<-gls(TS~TBV, data=gpulse_avg, correlation=corPagel(0.6012342,phy=gpulse.tree,form=~Species, fixed=TRUE))

CB_all<-gls(CB~TBV, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
CB_egenrecpt<-gls(CB~TBV+hyp_egenrecpt, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
CB_eolin<-gls(CB~TBV+hyp_eolin, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
CB_lam2<-as.numeric(CB_eolin$modelStruct)
CB_lam<-as.numeric(CB_egenrecpt$modelStruct)
if (CB_lam > 1)
	CB_lam = 1
if (CB_lam < 0)
	CB_lam = 0
if (CB_lam2 > 1)
	CB_lam2 = 1
if (CB_lam2 < 0)
	CB_lam2 = 0
CB_tubegen<-gls(CB~TBV, data=tubegen_avg, correlation=corPagel(CB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
CB_ampegen<-gls(CB~TBV, data=ampegen_avg, correlation=corPagel(CB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
CB_out<-gls(CB~TBV, data=out_avg, correlation=corPagel(CB_lam,phy=out.tree,form=~Species, fixed=TRUE))
CB_gym<-gls(CB~TBV, data=gym_avg, correlation=corPagel(CB_lam2,phy=gym.tree,form=~Species, fixed=TRUE))
CB_morm<-gls(CB~TBV, data=morm_avg, correlation=corPagel(CB_lam2,phy=morm.tree,form=~Species, fixed=TRUE))
CB_allgout<-gls(CB~TBV, data=allgout_avg, correlation=corPagel(CB_lam2,phy=allgout.tree,form=~Species, fixed=TRUE))
CB_allmout<-gls(CB~TBV, data=allmout_avg, correlation=corPagel(CB_lam2,phy=allmout.tree,form=~Species, fixed=TRUE))
CB_gwave<-gls(CB~TBV, data=gwave_avg, correlation=corPagel(0.7507515,phy=gwave.tree,form=~Species, fixed=TRUE))
CB_gpulse<-gls(CB~TBV, data=gpulse_avg, correlation=corPagel(0.7507515,phy=gpulse.tree,form=~Species, fixed=TRUE))

ROB_all<-gls(ROB~TBV, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
ROB_egenrecpt<-gls(ROB~TBV+hyp_egenrecpt, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
ROB_eolin<-gls(ROB~TBV+hyp_eolin, data=avg_log_data, correlation=corPagel(.8,phy=tree,form=~Species))
ROB_lam2<-as.numeric(ROB_eolin$modelStruct)
ROB_lam<-as.numeric(ROB_egenrecpt$modelStruct)
if (ROB_lam > 1)
	ROB_lam = 1
if (ROB_lam < 0)
	ROB_lam = 0
if (ROB_lam2 > 1)
	ROB_lam2 = 1
if (ROB_lam2 < 0)
	ROB_lam2 = 0
ROB_tubegen<-gls(ROB~TBV, data=tubegen_avg, correlation=corPagel(ROB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
ROB_ampegen<-gls(ROB~TBV, data=ampegen_avg, correlation=corPagel(ROB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
ROB_out<-gls(ROB~TBV, data=out_avg, correlation=corPagel(ROB_lam,phy=out.tree,form=~Species, fixed=TRUE))
ROB_gym<-gls(ROB~TBV, data=gym_avg, correlation=corPagel(ROB_lam2,phy=gym.tree,form=~Species, fixed=TRUE))
ROB_morm<-gls(ROB~TBV, data=morm_avg, correlation=corPagel(ROB_lam2,phy=morm.tree,form=~Species, fixed=TRUE))
ROB_allgout<-gls(ROB~TBV, data=allgout_avg, correlation=corPagel(ROB_lam2,phy=allgout.tree,form=~Species, fixed=TRUE))
ROB_allmout<-gls(ROB~TBV, data=allmout_avg, correlation=corPagel(ROB_lam2,phy=allmout.tree,form=~Species, fixed=TRUE))
ROB_gwave<-gls(ROB~TBV, data=gwave_avg, correlation=corPagel(0,phy=gwave.tree,form=~Species, fixed=TRUE))
ROB_gpulse<-gls(ROB~TBV, data=gpulse_avg, correlation=corPagel(0,phy=gpulse.tree,form=~Species, fixed=TRUE))

## morm line excl gymnarchus
nogymmorm_avg<-subset(morm_avg, Species != "Gymnarchus_niloticus")
nogymmorm.tree<-drop.tip(tree, setdiff(tree$tip.label, as.character(nogymmorm_avg$Species)))
OB_nogymmorm<-gls(OB~TBV, data=nogymmorm_avg, correlation=corPagel(0.8476816,phy=nogymmorm.tree,form=~Species, fixed=TRUE))
TEL_nogymmorm<-gls(TEL~TBV, data=nogymmorm_avg, correlation=corPagel(0.8626002,phy=nogymmorm.tree,form=~Species, fixed=TRUE))
HB_nogymmorm<-gls(HB~TBV, data=nogymmorm_avg, correlation=corPagel(0.7417891,phy=nogymmorm.tree,form=~Species, fixed=TRUE))
OT_nogymmorm<-gls(OT~TBV, data=nogymmorm_avg, correlation=corPagel(1,phy=nogymmorm.tree,form=~Species, fixed=TRUE))
TS_nogymmorm<-gls(TS~TBV, data=nogymmorm_avg, correlation=corPagel(0.6713374,phy=nogymmorm.tree,form=~Species, fixed=TRUE))
CB_nogymmorm<-gls(CB~TBV, data=nogymmorm_avg, correlation=corPagel(0.9435256,phy=nogymmorm.tree,form=~Species, fixed=TRUE))
ROB_nogymmorm<-gls(ROB~TBV, data=nogymmorm_avg, correlation=corPagel(0,phy=nogymmorm.tree,form=~Species, fixed=TRUE))


###########
#  plots  #
###########

## plotting the lines for eo + tuberous vs eo + ampullary vs not electric
morm_ob_plot_egenrecpt<-ggplot(comb_log_data, aes(x=Total.Brain.Volume, y=OB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Olfactory Bulbs", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) + 
	geom_abline(slope = OB_tubegen$coefficients[2], intercept = OB_tubegen$coefficients[1],color="Black")+geom_abline(slope = OB_ampegen$coefficients[2], intercept = OB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OB_out$coefficients[2], intercept = OB_out$coefficients[1],color="Black",lty=2)
morm_tel_plot_egenrecpt<-ggplot(comb_log_data, aes(x=Total.Brain.Volume, y=TEL, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Telencephalon", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TEL_tubegen$coefficients[2], intercept = TEL_tubegen$coefficients[1],color="Black")+geom_abline(slope = TEL_ampegen$coefficients[2], intercept = TEL_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TEL_out$coefficients[2], intercept = TEL_out$coefficients[1],color="Black",lty=2)
morm_ellhb_plot_egenrecpt<-ggplot(comb_log_data, aes(x=Total.Brain.Volume, y=HB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Hind Brain", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = HB_tubegen$coefficients[2], intercept = HB_tubegen$coefficients[1],color="Black")+geom_abline(slope = HB_ampegen$coefficients[2], intercept = HB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = HB_out$coefficients[2], intercept = HB_out$coefficients[1],color="Black",lty=2)
morm_ot_plot_egenrecpt<-ggplot(comb_log_data, aes(x=Total.Brain.Volume, y=OT, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Optic Tectum", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = OT_tubegen$coefficients[2], intercept = OT_tubegen$coefficients[1],color="Black")+geom_abline(slope = OT_ampegen$coefficients[2], intercept = OT_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OT_out$coefficients[2], intercept = OT_out$coefficients[1],color="Black",lty=2)
morm_ts_plot_egenrecpt<-ggplot(comb_log_data, aes(x=Total.Brain.Volume, y=TS, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Torus Semicircularis", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TS_tubegen$coefficients[2], intercept = TS_tubegen$coefficients[1],color="Black")+geom_abline(slope = TS_ampegen$coefficients[2], intercept = TS_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TS_out$coefficients[2], intercept = TS_out$coefficients[1],color="Black",lty=2)
morm_cb_plot_egenrecpt<-ggplot(comb_log_data, aes(x=Total.Brain.Volume, y=CB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Cerebellum", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = CB_tubegen$coefficients[2], intercept = CB_tubegen$coefficients[1],color="Black")+geom_abline(slope = CB_ampegen$coefficients[2], intercept = CB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = CB_out$coefficients[2], intercept = CB_out$coefficients[1],color="Black",lty=2)
morm_rob_plot_egenrecpt<-ggplot(comb_log_data, aes(x=Total.Brain.Volume, y=RoB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Rest of Brain", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = ROB_tubegen$coefficients[2], intercept = ROB_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROB_ampegen$coefficients[2], intercept = ROB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROB_out$coefficients[2], intercept = ROB_out$coefficients[1],color="Black",lty=2)

morm_all_plot_egenrecpt<-ggarrange(morm_ob_plot_egenrecpt,morm_tel_plot_egenrecpt,morm_ellhb_plot_egenrecpt,morm_ot_plot_egenrecpt,morm_ts_plot_egenrecpt,morm_cb_plot_egenrecpt,morm_rob_plot_egenrecpt, ncol=3, nrow=3, common.legend = TRUE, legend = "none");
morm_all_plot_egenrecpt2<-ggarrange(morm_ob_plot_egenrecpt,morm_tel_plot_egenrecpt,morm_ellhb_plot_egenrecpt,morm_ot_plot_egenrecpt,morm_ts_plot_egenrecpt,morm_cb_plot_egenrecpt,morm_rob_plot_egenrecpt, ncol=4, nrow=2, common.legend = TRUE, legend = "none");


## plotting morms vs gyms
morm_ob_plot_mormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=OB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Olfactory Bulbs", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) + 
	geom_abline(slope = OB_gym$coefficients[2], intercept = OB_gym$coefficients[1],color="Black")+geom_abline(slope = OB_morm$coefficients[2], intercept = OB_morm$coefficients[1],color="Black",lty=2)
morm_tel_plot_mormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=TEL, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Telencephalon", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TEL_gym$coefficients[2], intercept = TEL_gym$coefficients[1],color="Black")+geom_abline(slope = TEL_morm$coefficients[2], intercept = TEL_morm$coefficients[1],color="Black",lty=2)
morm_ellhb_plot_mormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=HB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Hindbrain", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = HB_gym$coefficients[2], intercept = HB_gym$coefficients[1],color="Black")+geom_abline(slope = HB_morm$coefficients[2], intercept = HB_morm$coefficients[1],color="Black",lty=2)
morm_ot_plot_mormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=OT, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Optic Tectum", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = OT_gym$coefficients[2], intercept = OT_gym$coefficients[1],color="Black")+geom_abline(slope = OT_morm$coefficients[2], intercept = OT_morm$coefficients[1],color="Black",lty=2)
morm_ts_plot_mormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=TS, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Torus Semicircularis", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TS_gym$coefficients[2], intercept = TS_gym$coefficients[1],color="Black")+geom_abline(slope = TS_morm$coefficients[2], intercept = TS_morm$coefficients[1],color="Black",lty=2)
morm_cb_plot_mormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=CB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Cerebellum", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = CB_gym$coefficients[2], intercept = CB_gym$coefficients[1],color="Black")+geom_abline(slope = CB_morm$coefficients[2], intercept = CB_morm$coefficients[1],color="Black",lty=2)
morm_rob_plot_mormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=RoB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Rest of Brain", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = ROB_gym$coefficients[2], intercept = ROB_gym$coefficients[1],color="Black")+geom_abline(slope = ROB_morm$coefficients[2], intercept = ROB_morm$coefficients[1],color="Black",lty=2)

morm_all_plot_mormgym<-ggarrange(morm_ob_plot_mormgym,morm_tel_plot_mormgym,morm_ellhb_plot_mormgym,morm_ot_plot_mormgym,morm_ts_plot_mormgym,morm_cb_plot_mormgym,morm_rob_plot_mormgym, ncol=3, nrow=3, common.legend = TRUE, legend = "none");
morm_all_plot_mormgym2<-ggarrange(morm_ob_plot_mormgym,morm_tel_plot_mormgym,morm_ellhb_plot_mormgym,morm_ot_plot_mormgym,morm_ts_plot_mormgym,morm_cb_plot_mormgym,morm_rob_plot_mormgym, ncol=4, nrow=2, common.legend = TRUE, legend = "none");

## plotting morms vs gyms + morm line excl gymnarchus
morm_ob_plot_nogymmormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=OB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Olfactory Bulbs", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) + 
	geom_abline(slope = OB_gym$coefficients[2], intercept = OB_gym$coefficients[1],color="Black")+geom_abline(slope = OB_morm$coefficients[2], intercept = OB_morm$coefficients[1],color="Black",lty=2) +
	geom_abline(slope = OB_nogymmorm$coefficients[2], intercept = OB_nogymmorm$coefficients[1],color="gray",lty=2)
morm_tel_plot_nogymmormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=TEL, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Telencephalon", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TEL_gym$coefficients[2], intercept = TEL_gym$coefficients[1],color="Black")+geom_abline(slope = TEL_morm$coefficients[2], intercept = TEL_morm$coefficients[1],color="Black",lty=2) +
	geom_abline(slope = TEL_nogymmorm$coefficients[2], intercept = TEL_nogymmorm$coefficients[1],color="gray",lty=2)
morm_ellhb_plot_nogymmormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=HB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Hindbrain", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = HB_gym$coefficients[2], intercept = HB_gym$coefficients[1],color="Black")+geom_abline(slope = HB_morm$coefficients[2], intercept = HB_morm$coefficients[1],color="Black",lty=2) +
	geom_abline(slope = HB_nogymmorm$coefficients[2], intercept = HB_nogymmorm$coefficients[1],color="gray",lty=2)
morm_ot_plot_nogymmormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=OT, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Optic Tectum", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = OT_gym$coefficients[2], intercept = OT_gym$coefficients[1],color="Black")+geom_abline(slope = OT_morm$coefficients[2], intercept = OT_morm$coefficients[1],color="Black",lty=2) +
	geom_abline(slope = OT_nogymmorm$coefficients[2], intercept = OT_nogymmorm$coefficients[1],color="gray",lty=2)
morm_ts_plot_nogymmormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=TS, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Torus Semicircularis", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TS_gym$coefficients[2], intercept = TS_gym$coefficients[1],color="Black")+geom_abline(slope = TS_morm$coefficients[2], intercept = TS_morm$coefficients[1],color="Black",lty=2) +
	geom_abline(slope = TS_nogymmorm$coefficients[2], intercept = TS_nogymmorm$coefficients[1],color="gray",lty=2)
morm_cb_plot_nogymmormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=CB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Cerebellum", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = CB_gym$coefficients[2], intercept = CB_gym$coefficients[1],color="Black")+geom_abline(slope = CB_morm$coefficients[2], intercept = CB_morm$coefficients[1],color="Black",lty=2) +
	geom_abline(slope = CB_nogymmorm$coefficients[2], intercept = CB_nogymmorm$coefficients[1],color="gray",lty=2)
morm_rob_plot_nogymmormgym<-ggplot(tubegen_sub, aes(x=Total.Brain.Volume, y=RoB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Rest of Brain", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = ROB_gym$coefficients[2], intercept = ROB_gym$coefficients[1],color="Black")+geom_abline(slope = ROB_morm$coefficients[2], intercept = ROB_morm$coefficients[1],color="Black",lty=2) +
	geom_abline(slope = ROB_nogymmorm$coefficients[2], intercept = ROB_nogymmorm$coefficients[1],color="gray",lty=2)

morm_all_plot_nogymmormgym<-ggarrange(morm_ob_plot_nogymmormgym,morm_tel_plot_nogymmormgym,morm_ellhb_plot_nogymmormgym,morm_ot_plot_nogymmormgym,morm_ts_plot_nogymmormgym,morm_cb_plot_nogymmormgym,morm_rob_plot_nogymmormgym, ncol=3, nrow=3, common.legend = TRUE, legend = "none");
morm_all_plot_nogymmormgym2<-ggarrange(morm_ob_plot_nogymmormgym,morm_tel_plot_nogymmormgym,morm_ellhb_plot_nogymmormgym,morm_ot_plot_nogymmormgym,morm_ts_plot_nogymmormgym,morm_cb_plot_nogymmormgym,morm_rob_plot_nogymmormgym, ncol=4, nrow=2, common.legend = TRUE, legend = "none");


## plotting wave vs pulse gyms
g_ob_plot_gpulsegwave<-ggplot(gym_sub, aes(x=Total.Brain.Volume, y=OB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Olfactory Bulbs", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) + 
	geom_abline(slope = OB_gwave$coefficients[2], intercept = OB_gwave$coefficients[1],color="Blue")+geom_abline(slope = OB_gpulse$coefficients[2], intercept = OB_gpulse$coefficients[1],color="seagreen2",lty=2)
g_tel_plot_gpulsegwave<-ggplot(gym_sub, aes(x=Total.Brain.Volume, y=TEL, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Telencephalon", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TEL_gwave$coefficients[2], intercept = TEL_gwave$coefficients[1],color="Blue")+geom_abline(slope = TEL_gpulse$coefficients[2], intercept = TEL_gpulse$coefficients[1],color="seagreen2",lty=2)
g_ellhb_plot_gpulsegwave<-ggplot(gym_sub, aes(x=Total.Brain.Volume, y=HB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Hindbrain", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = HB_gwave$coefficients[2], intercept = HB_gwave$coefficients[1],color="Blue")+geom_abline(slope = HB_gpulse$coefficients[2], intercept = HB_gpulse$coefficients[1],color="seagreen2",lty=2)
g_ot_plot_gpulsegwave<-ggplot(gym_sub, aes(x=Total.Brain.Volume, y=OT, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Optic Tectum", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = OT_gwave$coefficients[2], intercept = OT_gwave$coefficients[1],color="Blue")+geom_abline(slope = OT_gpulse$coefficients[2], intercept = OT_gpulse$coefficients[1],color="seagreen2",lty=2)
g_ts_plot_gpulsegwave<-ggplot(gym_sub, aes(x=Total.Brain.Volume, y=TS, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Torus Semicircularis", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TS_gwave$coefficients[2], intercept = TS_gwave$coefficients[1],color="Blue")+geom_abline(slope = TS_gpulse$coefficients[2], intercept = TS_gpulse$coefficients[1],color="seagreen2",lty=2)
g_cb_plot_gpulsegwave<-ggplot(gym_sub, aes(x=Total.Brain.Volume, y=CB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Cerebellum", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = CB_gwave$coefficients[2], intercept = CB_gwave$coefficients[1],color="Blue")+geom_abline(slope = CB_gpulse$coefficients[2], intercept = CB_gpulse$coefficients[1],color="seagreen2",lty=2)
g_rob_plot_gpulsegwave<-ggplot(gym_sub, aes(x=Total.Brain.Volume, y=RoB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Rest of Brain", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = ROB_gwave$coefficients[2], intercept = ROB_gwave$coefficients[1],color="Blue")+geom_abline(slope = ROB_gpulse$coefficients[2], intercept = ROB_gpulse$coefficients[1],color="seagreen2",lty=2)

g_all_plot_gpulsegwave<-ggarrange(g_ob_plot_gpulsegwave,g_tel_plot_gpulsegwave,g_ellhb_plot_gpulsegwave,g_ot_plot_gpulsegwave,g_ts_plot_gpulsegwave,g_cb_plot_gpulsegwave,g_rob_plot_gpulsegwave, ncol=3, nrow=3, common.legend = TRUE, legend = "none");
g_all_plot_gpulsegwave2<-ggarrange(g_ob_plot_gpulsegwave,g_tel_plot_gpulsegwave,g_ellhb_plot_gpulsegwave,g_ot_plot_gpulsegwave,g_ts_plot_gpulsegwave,g_cb_plot_gpulsegwave,g_rob_plot_gpulsegwave, ncol=4, nrow=2, common.legend = TRUE, legend = "none");


## plotting allmouts vs allgouts
allmout_ob_plot_moutgout<-ggplot(out_sub, aes(x=Total.Brain.Volume, y=OB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Olfactory Bulbs", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) + 
	geom_abline(slope = OB_allgout$coefficients[2], intercept = OB_allgout$coefficients[1],color="gray")+geom_abline(slope = OB_allmout$coefficients[2], intercept = OB_allmout$coefficients[1],color="gray",lty=2)
allmout_tel_plot_moutgout<-ggplot(out_sub, aes(x=Total.Brain.Volume, y=TEL, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Telencephalon", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TEL_allgout$coefficients[2], intercept = TEL_allgout$coefficients[1],color="gray")+geom_abline(slope = TEL_allmout$coefficients[2], intercept = TEL_allmout$coefficients[1],color="gray",lty=2)
allmout_ellhb_plot_moutgout<-ggplot(out_sub, aes(x=Total.Brain.Volume, y=HB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Hind Brain", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = HB_allgout$coefficients[2], intercept = HB_allgout$coefficients[1],color="gray")+geom_abline(slope = HB_allmout$coefficients[2], intercept = HB_allmout$coefficients[1],color="gray",lty=2)
allmout_ot_plot_moutgout<-ggplot(out_sub, aes(x=Total.Brain.Volume, y=OT, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Optic Tectum", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = OT_allgout$coefficients[2], intercept = OT_allgout$coefficients[1],color="gray")+geom_abline(slope = OT_allmout$coefficients[2], intercept = OT_allmout$coefficients[1],color="gray",lty=2)
allmout_ts_plot_moutgout<-ggplot(out_sub, aes(x=Total.Brain.Volume, y=TS, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Torus Semicircularis", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TS_allgout$coefficients[2], intercept = TS_allgout$coefficients[1],color="gray")+geom_abline(slope = TS_allmout$coefficients[2], intercept = TS_allmout$coefficients[1],color="gray",lty=2)
allmout_cb_plot_moutgout<-ggplot(out_sub, aes(x=Total.Brain.Volume, y=CB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Cerebellum", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = CB_allgout$coefficients[2], intercept = CB_allgout$coefficients[1],color="gray")+geom_abline(slope = CB_allmout$coefficients[2], intercept = CB_allmout$coefficients[1],color="gray",lty=2)
allmout_rob_plot_moutgout<-ggplot(out_sub, aes(x=Total.Brain.Volume, y=RoB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Rest of Brain", x=bquote('Log Total Brain Volume' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = ROB_allgout$coefficients[2], intercept = ROB_allgout$coefficients[1],color="gray")+geom_abline(slope = ROB_allmout$coefficients[2], intercept = ROB_allmout$coefficients[1],color="gray",lty=2)

allmout_all_plot_moutgout<-ggarrange(allmout_ob_plot_moutgout,allmout_tel_plot_moutgout,allmout_ellhb_plot_moutgout,allmout_ot_plot_moutgout,allmout_ts_plot_moutgout,allmout_cb_plot_moutgout,allmout_rob_plot_moutgout, ncol=3, nrow=3, common.legend = TRUE, legend = "none");
allmout_all_plot_moutgout2<-ggarrange(allmout_ob_plot_moutgout,allmout_tel_plot_moutgout,allmout_ellhb_plot_moutgout,allmout_ot_plot_moutgout,allmout_ts_plot_moutgout,allmout_cb_plot_moutgout,allmout_rob_plot_moutgout, ncol=4, nrow=2, common.legend = TRUE, legend = "none");


## all plots together
all_3_comp_plots<-ggarrange(morm_ob_plot_mormgym+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),morm_tel_plot_mormgym+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),morm_ellhb_plot_mormgym+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),morm_ot_plot_mormgym+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),morm_ts_plot_mormgym+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),morm_cb_plot_mormgym+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),morm_rob_plot_mormgym+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),
	allmout_ob_plot_moutgout+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),allmout_tel_plot_moutgout+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),allmout_ellhb_plot_moutgout+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),allmout_ot_plot_moutgout+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),allmout_ts_plot_moutgout+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),allmout_cb_plot_moutgout+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),allmout_rob_plot_moutgout+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),
	g_ob_plot_gpulsegwave+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),g_tel_plot_gpulsegwave+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),g_ellhb_plot_gpulsegwave+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),g_ot_plot_gpulsegwave+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),g_ts_plot_gpulsegwave+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),g_cb_plot_gpulsegwave+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),g_rob_plot_gpulsegwave+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()), 
	ncol=7, nrow=3, common.legend = TRUE, legend = "none", labels=c("A","","","","","","","B","","","","","","","C","","","","","",""), heights=c(1,.9,.9));
all_3_comp_plots2<-annotate_figure(all_3_comp_plots, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90, size = 11),
                    bottom = text_grob(bquote('Log Total Brain Volume' ~(mm^3)), size = 11))
all_3_comp_plots3<-ggarrange(all_3_comp_plots2,text.p, ncol=1, nrow=2, heights=c(1,.1))



#########################
# pgls lines, tbv - roi #
#########################

tree<-all_tree_sg 			## change this to represent which tree I want to use in the pgls

## subset for eo + recept type
tubegen_avg<-subset(avg_comb_data, Species %in% tubegen_names)
ampegen_avg<-subset(avg_comb_data, Species %in% ampegen_names)
out_avg<-subset(avg_comb_data, Species %in% out_names)
gym_avg<-subset(avg_comb_data, Species %in% gym_names)
morm_avg<-subset(avg_comb_data, Species %in% morm_names)
allgout_avg<-subset(avg_comb_data, Species %in% allgout_names)
allmout_avg<-subset(avg_comb_data, Species %in% allmout_names)
gwave_avg<-subset(avg_comb_data, Species %in% gwave_names)
gpulse_avg<-subset(avg_comb_data, Species %in% gpulse_names)

## pgls lines
OB_all<-gls(OB_M~TBV_OB+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
OB_lam<-as.numeric(OB_all$modelStruct)
if (OB_lam > 1)
	OB_lam = 1
OB_tubegen<-gls(OB_M~TBV_OB, data=tubegen_avg, correlation=corPagel(OB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OB_ampegen<-gls(OB_M~TBV_OB, data=ampegen_avg, correlation=corPagel(OB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OB_out<-gls(OB_M~TBV_OB, data=out_avg, correlation=corPagel(OB_lam,phy=out.tree,form=~Species, fixed=TRUE))
OB_gym<-gls(OB_M~TBV_OB, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
OB_morm<-gls(OB_M~TBV_OB, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
OB_allgout<-gls(OB_M~TBV_OB, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
OB_allmout<-gls(OB_M~TBV_OB, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
OB_gwave<-gls(OB_M~TBV_OB, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OB_gpulse<-gls(OB_M~TBV_OB, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

TEL_all<-gls(TEL_M~TBV_TEL+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
TEL_lam<-as.numeric(TEL_all$modelStruct)
if (TEL_lam > 1)
	TEL_lam = 1
TEL_tubegen<-gls(TEL_M~TBV_TEL, data=tubegen_avg, correlation=corPagel(TEL_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TEL_ampegen<-gls(TEL_M~TBV_TEL, data=ampegen_avg, correlation=corPagel(TEL_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TEL_out<-gls(TEL_M~TBV_TEL, data=out_avg, correlation=corPagel(TEL_lam,phy=out.tree,form=~Species, fixed=TRUE))
TEL_gym<-gls(TEL_M~TBV_TEL, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
TEL_morm<-gls(TEL_M~TBV_TEL, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
TEL_allgout<-gls(TEL_M~TBV_TEL, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
TEL_allmout<-gls(TEL_M~TBV_TEL, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
TEL_gwave<-gls(TEL_M~TBV_TEL, data=gwave_avg, correlation=corPagel(0.6797043,phy=gwave.tree,form=~Species, fixed=TRUE))
TEL_gpulse<-gls(TEL_M~TBV_TEL, data=gpulse_avg, correlation=corPagel(0.6797043,phy=gpulse.tree,form=~Species, fixed=TRUE))

HB_all<-gls(HB_M~TBV_HB+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
HB_lam<-as.numeric(HB_all$modelStruct)
if (HB_lam > 1)
	HB_lam = 1
HB_tubegen<-gls(HB_M~TBV_HB, data=tubegen_avg, correlation=corPagel(HB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
HB_ampegen<-gls(HB_M~TBV_HB, data=ampegen_avg, correlation=corPagel(HB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
HB_out<-gls(HB_M~TBV_HB, data=out_avg, correlation=corPagel(HB_lam,phy=out.tree,form=~Species, fixed=TRUE))
HB_gym<-gls(HB_M~TBV_HB, data=gym_avg, correlation=corPagel(0.838,phy=gym.tree,form=~Species, fixed=TRUE))
HB_morm<-gls(HB_M~TBV_HB, data=morm_avg, correlation=corPagel(0.838,phy=morm.tree,form=~Species, fixed=TRUE))
HB_allgout<-gls(HB_M~TBV_HB, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
HB_allmout<-gls(HB_M~TBV_HB, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
HB_gwave<-gls(HB_M~TBV_HB, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
HB_gpulse<-gls(HB_M~TBV_HB, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

OT_all<-gls(OT_M~TBV_OT+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
OT_lam<-as.numeric(OT_all$modelStruct)
if (OT_lam > 1)
	OT_lam = 1
OT_tubegen<-gls(OT_M~TBV_OT, data=tubegen_avg, correlation=corPagel(OT_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OT_ampegen<-gls(OT_M~TBV_OT, data=ampegen_avg, correlation=corPagel(OT_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OT_out<-gls(OT_M~TBV_OT, data=out_avg, correlation=corPagel(OT_lam,phy=out.tree,form=~Species, fixed=TRUE))
OT_gym<-gls(OT_M~TBV_OT, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
OT_morm<-gls(OT_M~TBV_OT, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
OT_allgout<-gls(OT_M~TBV_OT, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
OT_allmout<-gls(OT_M~TBV_OT, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
OT_gwave<-gls(OT_M~TBV_OT, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OT_gpulse<-gls(OT_M~TBV_OT, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

TS_all<-gls(TS_M~TBV_TS+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
TS_lam<-as.numeric(TS_all$modelStruct)
if (TS_lam > 1)
	TS_lam = 1
if (TS_lam < 0)
	TS_lam = 0
TS_tubegen<-gls(TS_M~TBV_TS, data=tubegen_avg, correlation=corPagel(TS_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TS_ampegen<-gls(TS_M~TBV_TS, data=ampegen_avg, correlation=corPagel(TS_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TS_out<-gls(TS_M~TBV_TS, data=out_avg, correlation=corPagel(TS_lam,phy=out.tree,form=~Species, fixed=TRUE))
TS_gym<-gls(TS_M~TBV_TS, data=gym_avg, correlation=corPagel(0.056,phy=gym.tree,form=~Species, fixed=TRUE))
TS_morm<-gls(TS_M~TBV_TS, data=morm_avg, correlation=corPagel(0.056,phy=morm.tree,form=~Species, fixed=TRUE))
TS_allgout<-gls(TS_M~TBV_TS, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
TS_allmout<-gls(TS_M~TBV_TS, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
TS_gwave<-gls(TS_M~TBV_TS, data=gwave_avg, correlation=corPagel(0.5631124,phy=gwave.tree,form=~Species, fixed=TRUE))
TS_gpulse<-gls(TS_M~TBV_TS, data=gpulse_avg, correlation=corPagel(0.5631124,phy=gpulse.tree,form=~Species, fixed=TRUE))

CB_all<-gls(CB_M~TBV_CB+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
CB_lam<-as.numeric(CB_all$modelStruct)
if (CB_lam > 1)
	CB_lam = 1
CB_tubegen<-gls(CB_M~TBV_CB, data=tubegen_avg, correlation=corPagel(CB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
CB_ampegen<-gls(CB_M~TBV_CB, data=ampegen_avg, correlation=corPagel(CB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
CB_out<-gls(CB_M~TBV_CB, data=out_avg, correlation=corPagel(CB_lam,phy=out.tree,form=~Species, fixed=TRUE))
CB_gym<-gls(CB_M~TBV_CB, data=gym_avg, correlation=corPagel(0.8520273,phy=gym.tree,form=~Species, fixed=TRUE))
CB_morm<-gls(CB_M~TBV_CB, data=morm_avg, correlation=corPagel(0.8520273,phy=morm.tree,form=~Species, fixed=TRUE))
CB_allgout<-gls(CB_M~TBV_CB, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
CB_allmout<-gls(CB_M~TBV_CB, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
CB_gwave<-gls(CB_M~TBV_CB, data=gwave_avg, correlation=corPagel(0.7886864,phy=gwave.tree,form=~Species, fixed=TRUE))
CB_gpulse<-gls(CB_M~TBV_CB, data=gpulse_avg, correlation=corPagel(0.7886864,phy=gpulse.tree,form=~Species, fixed=TRUE))

ROB_all<-gls(RoB_M~TBV_ROB+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
ROB_lam<-as.numeric(ROB_all$modelStruct)
if (ROB_lam > 1)
	ROB_lam = 1
ROB_tubegen<-gls(RoB_M~TBV_ROB, data=tubegen_avg, correlation=corPagel(ROB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
ROB_ampegen<-gls(RoB_M~TBV_ROB, data=ampegen_avg, correlation=corPagel(ROB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
ROB_out<-gls(RoB_M~TBV_ROB, data=out_avg, correlation=corPagel(ROB_lam,phy=out.tree,form=~Species, fixed=TRUE))
ROB_gym<-gls(RoB_M~TBV_ROB, data=gym_avg, correlation=corPagel(0.8200716,phy=gym.tree,form=~Species, fixed=TRUE))
ROB_morm<-gls(RoB_M~TBV_ROB, data=morm_avg, correlation=corPagel(0.8200716,phy=morm.tree,form=~Species, fixed=TRUE))
ROB_allgout<-gls(RoB_M~TBV_ROB, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
ROB_allmout<-gls(RoB_M~TBV_ROB, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
ROB_gwave<-gls(RoB_M~TBV_ROB, data=gwave_avg, correlation=corPagel(0,phy=gwave.tree,form=~Species, fixed=TRUE))
ROB_gpulse<-gls(RoB_M~TBV_ROB, data=gpulse_avg, correlation=corPagel(0,phy=gpulse.tree,form=~Species, fixed=TRUE))

######################
#  plots, tbv - roi  #
######################

## plotting the lines for eo + tuberous vs eo + ampullary vs not electric, tbv - roi
morm_ob_plot_egenrecpt_roi<-ggplot(comb_log_data, aes(x=tb.ob, y=OB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Olfactory Bulbs", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) + 
	geom_abline(slope = OB_tubegen$coefficients[2], intercept = OB_tubegen$coefficients[1],color="Black")+geom_abline(slope = OB_ampegen$coefficients[2], intercept = OB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OB_out$coefficients[2], intercept = OB_out$coefficients[1],color="Black",lty=2)
morm_tel_plot_egenrecpt_roi<-ggplot(comb_log_data, aes(x=tb.tel, y=TEL, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Telencephalon", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TEL_tubegen$coefficients[2], intercept = TEL_tubegen$coefficients[1],color="Black")+geom_abline(slope = TEL_ampegen$coefficients[2], intercept = TEL_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TEL_out$coefficients[2], intercept = TEL_out$coefficients[1],color="Black",lty=2)
morm_ellhb_plot_egenrecpt_roi<-ggplot(comb_log_data, aes(x=tb.hb, y=HB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Hind Brain", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = HB_tubegen$coefficients[2], intercept = HB_tubegen$coefficients[1],color="Black")+geom_abline(slope = HB_ampegen$coefficients[2], intercept = HB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = HB_out$coefficients[2], intercept = HB_out$coefficients[1],color="Black",lty=2)
morm_ot_plot_egenrecpt_roi<-ggplot(comb_log_data, aes(x=tb.ot, y=OT, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Optic Tectum", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = OT_tubegen$coefficients[2], intercept = OT_tubegen$coefficients[1],color="Black")+geom_abline(slope = OT_ampegen$coefficients[2], intercept = OT_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OT_out$coefficients[2], intercept = OT_out$coefficients[1],color="Black",lty=2)
morm_ts_plot_egenrecpt_roi<-ggplot(comb_log_data, aes(x=tb.ts, y=TS, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Torus Semicircularis", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TS_tubegen$coefficients[2], intercept = TS_tubegen$coefficients[1],color="Black")+geom_abline(slope = TS_ampegen$coefficients[2], intercept = TS_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TS_out$coefficients[2], intercept = TS_out$coefficients[1],color="Black",lty=2)
morm_cb_plot_egenrecpt_roi<-ggplot(comb_log_data, aes(x=tb.cb, y=CB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Cerebellum", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = CB_tubegen$coefficients[2], intercept = CB_tubegen$coefficients[1],color="Black")+geom_abline(slope = CB_ampegen$coefficients[2], intercept = CB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = CB_out$coefficients[2], intercept = CB_out$coefficients[1],color="Black",lty=2)
morm_rob_plot_egenrecpt_roi<-ggplot(comb_log_data, aes(x=tb.rob, y=RoB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=c(21,21,21,22,23,21,22,21,22,22,23,22,24,25,23,23,24,21,23,8,22,4,24,25,24,25,8,24,4,25,8,23)) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Rest of Brain", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = ROB_tubegen$coefficients[2], intercept = ROB_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROB_ampegen$coefficients[2], intercept = ROB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROB_out$coefficients[2], intercept = ROB_out$coefficients[1],color="Black",lty=2)

morm_all_plot_egenrecpt_roi<-ggarrange(morm_ob_plot_egenrecpt_roi,morm_tel_plot_egenrecpt_roi,morm_ellhb_plot_egenrecpt_roi,morm_ot_plot_egenrecpt_roi,morm_ts_plot_egenrecpt_roi,morm_cb_plot_egenrecpt_roi,morm_rob_plot_egenrecpt_roi, ncol=3, nrow=3, common.legend = TRUE, legend = "none");
morm_all_plot_egenrecpt_roi2<-ggarrange(morm_ob_plot_egenrecpt_roi+theme(axis.title=element_text(size=9)),morm_tel_plot_egenrecpt_roi+theme(axis.title=element_text(size=9)),morm_ellhb_plot_egenrecpt_roi+theme(axis.title=element_text(size=9)),morm_ot_plot_egenrecpt_roi+theme(axis.title=element_text(size=9)),morm_ts_plot_egenrecpt_roi+theme(axis.title=element_text(size=9)),morm_cb_plot_egenrecpt_roi+theme(axis.title=element_text(size=9)),morm_rob_plot_egenrecpt_roi+theme(axis.title=element_text(size=9)), ncol=4, nrow=2, common.legend = TRUE, legend = "none");

## plotting morms vs gyms, tbv - roi
morm_ob_plot_mormgym_roi<-ggplot(tubegen_sub, aes(x=tb.ob, y=OB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Olfactory Bulbs", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) + 
	geom_abline(slope = OB_gym$coefficients[2], intercept = OB_gym$coefficients[1],color="Black")+geom_abline(slope = OB_morm$coefficients[2], intercept = OB_morm$coefficients[1],color="Black",lty=2)
morm_tel_plot_mormgym_roi<-ggplot(tubegen_sub, aes(x=tb.tel, y=TEL, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Telencephalon", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TEL_gym$coefficients[2], intercept = TEL_gym$coefficients[1],color="Black")+geom_abline(slope = TEL_morm$coefficients[2], intercept = TEL_morm$coefficients[1],color="Black",lty=2)
morm_ellhb_plot_mormgym_roi<-ggplot(tubegen_sub, aes(x=tb.hb, y=HB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Hindbrain", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = HB_gym$coefficients[2], intercept = HB_gym$coefficients[1],color="Black")+geom_abline(slope = HB_morm$coefficients[2], intercept = HB_morm$coefficients[1],color="Black",lty=2)
morm_ot_plot_mormgym_roi<-ggplot(tubegen_sub, aes(x=tb.ot, y=OT, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Optic Tectum", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = OT_gym$coefficients[2], intercept = OT_gym$coefficients[1],color="Black")+geom_abline(slope = OT_morm$coefficients[2], intercept = OT_morm$coefficients[1],color="Black",lty=2)
morm_ts_plot_mormgym_roi<-ggplot(tubegen_sub, aes(x=tb.ts, y=TS, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Torus Semicircularis", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TS_gym$coefficients[2], intercept = TS_gym$coefficients[1],color="Black")+geom_abline(slope = TS_morm$coefficients[2], intercept = TS_morm$coefficients[1],color="Black",lty=2)
morm_cb_plot_mormgym_roi<-ggplot(tubegen_sub, aes(x=tb.cb, y=CB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Cerebellum", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = CB_gym$coefficients[2], intercept = CB_gym$coefficients[1],color="Black")+geom_abline(slope = CB_morm$coefficients[2], intercept = CB_morm$coefficients[1],color="Black",lty=2)
morm_rob_plot_mormgym_roi<-ggplot(tubegen_sub, aes(x=tb.rob, y=RoB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Rest of Brain", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = ROB_gym$coefficients[2], intercept = ROB_gym$coefficients[1],color="Black")+geom_abline(slope = ROB_morm$coefficients[2], intercept = ROB_morm$coefficients[1],color="Black",lty=2)

morm_all_plot_mormgym_roi<-ggarrange(morm_ob_plot_mormgym_roi,morm_tel_plot_mormgym_roi,morm_ellhb_plot_mormgym_roi,morm_ot_plot_mormgym_roi,morm_ts_plot_mormgym_roi,morm_cb_plot_mormgym_roi,morm_rob_plot_mormgym_roi, ncol=3, nrow=3, common.legend = TRUE, legend = "none");
morm_all_plot_mormgym_roi2<-ggarrange(morm_ob_plot_mormgym_roi,morm_tel_plot_mormgym_roi,morm_ellhb_plot_mormgym_roi,morm_ot_plot_mormgym_roi,morm_ts_plot_mormgym_roi,morm_cb_plot_mormgym_roi,morm_rob_plot_mormgym_roi, ncol=4, nrow=2, common.legend = TRUE, legend = "none");


## plotting wave vs pulse gyms, tbv - roi
g_ob_plot_gpulsegwave_roi<-ggplot(gym_sub, aes(x=tb.ob, y=OB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Olfactory Bulbs", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) + 
	geom_abline(slope = OB_gwave$coefficients[2], intercept = OB_gwave$coefficients[1],color="Blue")+geom_abline(slope = OB_gpulse$coefficients[2], intercept = OB_gpulse$coefficients[1],color="seagreen2",lty=2)
g_tel_plot_gpulsegwave_roi<-ggplot(gym_sub, aes(x=tb.tel, y=TEL, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Telencephalon", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TEL_gwave$coefficients[2], intercept = TEL_gwave$coefficients[1],color="Blue")+geom_abline(slope = TEL_gpulse$coefficients[2], intercept = TEL_gpulse$coefficients[1],color="seagreen2",lty=2)
g_ellhb_plot_gpulsegwave_roi<-ggplot(gym_sub, aes(x=tb.hb, y=HB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Hindbrain", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = HB_gwave$coefficients[2], intercept = HB_gwave$coefficients[1],color="Blue")+geom_abline(slope = HB_gpulse$coefficients[2], intercept = HB_gpulse$coefficients[1],color="seagreen2",lty=2)
g_ot_plot_gpulsegwave_roi<-ggplot(gym_sub, aes(x=tb.ot, y=OT, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Optic Tectum", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = OT_gwave$coefficients[2], intercept = OT_gwave$coefficients[1],color="Blue")+geom_abline(slope = OT_gpulse$coefficients[2], intercept = OT_gpulse$coefficients[1],color="seagreen2",lty=2)
g_ts_plot_gpulsegwave_roi<-ggplot(gym_sub, aes(x=tb.ts, y=TS, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Torus Semicircularis", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TS_gwave$coefficients[2], intercept = TS_gwave$coefficients[1],color="Blue")+geom_abline(slope = TS_gpulse$coefficients[2], intercept = TS_gpulse$coefficients[1],color="seagreen2",lty=2)
g_cb_plot_gpulsegwave_roi<-ggplot(gym_sub, aes(x=tb.cb, y=CB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Cerebellum", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = CB_gwave$coefficients[2], intercept = CB_gwave$coefficients[1],color="Blue")+geom_abline(slope = CB_gpulse$coefficients[2], intercept = CB_gpulse$coefficients[1],color="seagreen2",lty=2)
g_rob_plot_gpulsegwave_roi<-ggplot(gym_sub, aes(x=tb.rob, y=RoB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Rest of Brain", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = spp_fills) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = ROB_gwave$coefficients[2], intercept = ROB_gwave$coefficients[1],color="Blue")+geom_abline(slope = ROB_gpulse$coefficients[2], intercept = ROB_gpulse$coefficients[1],color="seagreen2",lty=2)

g_all_plot_gpulsegwave_roi<-ggarrange(g_ob_plot_gpulsegwave_roi,g_tel_plot_gpulsegwave_roi,g_ellhb_plot_gpulsegwave_roi,g_ot_plot_gpulsegwave_roi,g_ts_plot_gpulsegwave_roi,g_cb_plot_gpulsegwave_roi,g_rob_plot_gpulsegwave_roi, ncol=3, nrow=3, common.legend = TRUE, legend = "none");
g_all_plot_gpulsegwave_roi2<-ggarrange(g_ob_plot_gpulsegwave_roi,g_tel_plot_gpulsegwave_roi,g_ellhb_plot_gpulsegwave_roi,g_ot_plot_gpulsegwave_roi,g_ts_plot_gpulsegwave_roi,g_cb_plot_gpulsegwave_roi,g_rob_plot_gpulsegwave_roi, ncol=4, nrow=2, common.legend = TRUE, legend = "none");


## plotting allmouts vs allgouts, tbv - roi
allmout_ob_plot_moutgout_roi<-ggplot(out_sub, aes(x=tb.ob, y=OB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Olfactory Bulbs", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) + 
	geom_abline(slope = OB_allgout$coefficients[2], intercept = OB_allgout$coefficients[1],color="gray")+geom_abline(slope = OB_allmout$coefficients[2], intercept = OB_allmout$coefficients[1],color="gray",lty=2)
allmout_tel_plot_moutgout_roi<-ggplot(out_sub, aes(x=tb.tel, y=TEL, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Telencephalon", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TEL_allgout$coefficients[2], intercept = TEL_allgout$coefficients[1],color="gray")+geom_abline(slope = TEL_allmout$coefficients[2], intercept = TEL_allmout$coefficients[1],color="gray",lty=2)
allmout_ellhb_plot_moutgout_roi<-ggplot(out_sub, aes(x=tb.hb, y=HB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Hind Brain", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = HB_allgout$coefficients[2], intercept = HB_allgout$coefficients[1],color="gray")+geom_abline(slope = HB_allmout$coefficients[2], intercept = HB_allmout$coefficients[1],color="gray",lty=2)
allmout_ot_plot_moutgout_roi<-ggplot(out_sub, aes(x=tb.ot, y=OT, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Optic Tectum", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = OT_allgout$coefficients[2], intercept = OT_allgout$coefficients[1],color="gray")+geom_abline(slope = OT_allmout$coefficients[2], intercept = OT_allmout$coefficients[1],color="gray",lty=2)
allmout_ts_plot_moutgout_roi<-ggplot(out_sub, aes(x=tb.ts, y=TS, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Torus Semicircularis", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = TS_allgout$coefficients[2], intercept = TS_allgout$coefficients[1],color="gray")+geom_abline(slope = TS_allmout$coefficients[2], intercept = TS_allmout$coefficients[1],color="gray",lty=2)
allmout_cb_plot_moutgout_roi<-ggplot(out_sub, aes(x=tb.cb, y=CB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Cerebellum", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = CB_allgout$coefficients[2], intercept = CB_allgout$coefficients[1],color="gray")+geom_abline(slope = CB_allmout$coefficients[2], intercept = CB_allmout$coefficients[1],color="gray",lty=2)
allmout_rob_plot_moutgout_roi<-ggplot(out_sub, aes(x=tb.rob, y=RoB, shape=Species, fill=Species, color=Species)) + scale_shape_manual(values=spp_shapes) + geom_point() + theme_classic() + theme(legend.position="bottom") + labs(title="Rest of Brain", x=bquote('Log (Total Brain - Region Volume)' ~(mm^3)), y = bquote('Log Region Volume' ~(mm^3))) + scale_fill_manual(values = c(     "white", "gray", "Black", "Black",      "Black",   "gray", "gray",  "white",  "Black",      "white")) + scale_colour_manual(values = sim_colors) +
	geom_abline(slope = ROB_allgout$coefficients[2], intercept = ROB_allgout$coefficients[1],color="gray")+geom_abline(slope = ROB_allmout$coefficients[2], intercept = ROB_allmout$coefficients[1],color="gray",lty=2)

allmout_all_plot_moutgout_roi<-ggarrange(allmout_ob_plot_moutgout_roi,allmout_tel_plot_moutgout_roi,allmout_ellhb_plot_moutgout_roi,allmout_ot_plot_moutgout_roi,allmout_ts_plot_moutgout_roi,allmout_cb_plot_moutgout_roi,allmout_rob_plot_moutgout_roi, ncol=3, nrow=3, common.legend = TRUE, legend = "none");
allmout_all_plot_moutgout_roi2<-ggarrange(allmout_ob_plot_moutgout_roi,allmout_tel_plot_moutgout_roi,allmout_ellhb_plot_moutgout_roi,allmout_ot_plot_moutgout_roi,allmout_ts_plot_moutgout_roi,allmout_cb_plot_moutgout_roi,allmout_rob_plot_moutgout_roi, ncol=4, nrow=2, common.legend = TRUE, legend = "none");


## all plots together, tbv - roi
all_3_comp_plots_roi<-ggarrange(morm_ob_plot_mormgym_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),morm_tel_plot_mormgym_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),morm_ellhb_plot_mormgym_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),morm_ot_plot_mormgym_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),morm_ts_plot_mormgym_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),morm_cb_plot_mormgym_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),morm_rob_plot_mormgym_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_text(size=11)),
	allmout_ob_plot_moutgout_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),allmout_tel_plot_moutgout_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),allmout_ellhb_plot_moutgout_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),allmout_ot_plot_moutgout_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),allmout_ts_plot_moutgout_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),allmout_cb_plot_moutgout_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),allmout_rob_plot_moutgout_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),
	g_ob_plot_gpulsegwave_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),g_tel_plot_gpulsegwave_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),g_ellhb_plot_gpulsegwave_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),g_ot_plot_gpulsegwave_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),g_ts_plot_gpulsegwave_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),g_cb_plot_gpulsegwave_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()),g_rob_plot_gpulsegwave_roi+theme(axis.text=element_text(size=8),axis.title=element_blank(),plot.title=element_blank()), 
	ncol=7, nrow=3, common.legend = TRUE, legend = "none", labels=c("A","","","","","","","B","","","","","","","C","","","","","",""), heights=c(1,.9,.9));
all_3_comp_plots_roi2<-annotate_figure(all_3_comp_plots_roi, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90, size = 11),
                    bottom = text_grob(bquote('Log (Total Brain - Region Volume)' ~(mm^3)), size = 11))
all_3_comp_plots_roi3<-ggarrange(all_3_comp_plots_roi2,text.p, ncol=1, nrow=2, heights=c(1,.1))

## make all 3 comp plots width = 7.08661, height = 4.41898

################################
# pgls lines, region by region #
################################

tree<-all_tree_sg 			## change this to represent which tree I want to use in the pgls

## subset for eo + recept type
tubegen_avg<-subset(avg_comb_data, Species %in% tubegen_names)
ampegen_avg<-subset(avg_comb_data, Species %in% ampegen_names)
out_avg<-subset(avg_comb_data, Species %in% out_names)

gym_avg<-subset(avg_comb_data, Species %in% gym_names)
morm_avg<-subset(avg_comb_data, Species %in% morm_names)
allgout_avg<-subset(avg_comb_data, Species %in% allgout_names)
allmout_avg<-subset(avg_comb_data, Species %in% allmout_names)
gwave_avg<-subset(avg_comb_data, Species %in% gwave_names)
gpulse_avg<-subset(avg_comb_data, Species %in% gpulse_names)

## pgls lines
## OB as x
OBTEL_all<-gls(TEL_M~OB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
OB_lam<-as.numeric(OBTEL_all$modelStruct)
if (OB_lam > 1)
	OB_lam = 1
OBTEL_tubegen<-gls(TEL_M~OB_M, data=tubegen_avg, correlation=corPagel(OB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OBTEL_ampegen<-gls(TEL_M~OB_M, data=ampegen_avg, correlation=corPagel(OB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OBTEL_out<-gls(TEL_M~OB_M, data=out_avg, correlation=corPagel(OB_lam,phy=out.tree,form=~Species, fixed=TRUE))
OBTEL_gym<-gls(TEL_M~OB_M, data=gym_avg, correlation=corPagel(0,phy=gym.tree,form=~Species, fixed=TRUE))
OBTEL_morm<-gls(TEL_M~OB_M, data=morm_avg, correlation=corPagel(0,phy=morm.tree,form=~Species, fixed=TRUE))
OBTEL_allgout<-gls(TEL_M~OB_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
OBTEL_allmout<-gls(TEL_M~OB_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
OBTEL_gwave<-gls(TEL_M~OB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OBTEL_gpulse<-gls(TEL_M~OB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

OBHB_all<-gls(HB_M~OB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(1,phy=tree,form=~Species))
OB_lam<-as.numeric(OBHB_all$modelStruct)
if (OB_lam > 1)
	OB_lam = 1
OBHB_tubegen<-gls(HB_M~OB_M, data=tubegen_avg, correlation=corPagel(OB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OBHB_ampegen<-gls(HB_M~OB_M, data=ampegen_avg, correlation=corPagel(OB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OBHB_out<-gls(HB_M~OB_M, data=out_avg, correlation=corPagel(OB_lam,phy=out.tree,form=~Species, fixed=TRUE))
OBHB_gym<-gls(HB_M~OB_M, data=gym_avg, correlation=corPagel(0,phy=gym.tree,form=~Species, fixed=TRUE))
OBHB_morm<-gls(HB_M~OB_M, data=morm_avg, correlation=corPagel(0,phy=morm.tree,form=~Species, fixed=TRUE))
OBHB_allgout<-gls(HB_M~OB_M, data=allgout_avg, correlation=corPagel(0.323637,phy=allgout.tree,form=~Species, fixed=TRUE))
OBHB_allmout<-gls(HB_M~OB_M, data=allmout_avg, correlation=corPagel(0.323637,phy=allmout.tree,form=~Species, fixed=TRUE))
OBHB_gwave<-gls(HB_M~OB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OBHB_gpulse<-gls(HB_M~OB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

OBOT_all<-gls(OT_M~OB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.9,phy=tree,form=~Species))
OB_lam<-as.numeric(OBOT_all$modelStruct)
if (OB_lam > 1)
	OB_lam = 1
OBOT_tubegen<-gls(OT_M~OB_M, data=tubegen_avg, correlation=corPagel(OB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OBOT_ampegen<-gls(OT_M~OB_M, data=ampegen_avg, correlation=corPagel(OB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OBOT_out<-gls(OT_M~OB_M, data=out_avg, correlation=corPagel(OB_lam,phy=out.tree,form=~Species, fixed=TRUE))
OBOT_gym<-gls(OT_M~OB_M, data=gym_avg, correlation=corPagel(0.8676742,phy=gym.tree,form=~Species, fixed=TRUE))
OBOT_morm<-gls(OT_M~OB_M, data=morm_avg, correlation=corPagel(0.8676742,phy=morm.tree,form=~Species, fixed=TRUE))
OBOT_allgout<-gls(OT_M~OB_M, data=allgout_avg, correlation=corPagel(0.8538995,phy=allgout.tree,form=~Species, fixed=TRUE))
OBOT_allmout<-gls(OT_M~OB_M, data=allmout_avg, correlation=corPagel(0.8538995,phy=allmout.tree,form=~Species, fixed=TRUE))
OBOT_gwave<-gls(OT_M~OB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OBOT_gpulse<-gls(OT_M~OB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

OBTS_all<-gls(TS_M~OB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
OB_lam<-as.numeric(OBTS_all$modelStruct)
if (OB_lam > 1)
	OB_lam = 1
OBTS_tubegen<-gls(TS_M~OB_M, data=tubegen_avg, correlation=corPagel(OB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OBTS_ampegen<-gls(TS_M~OB_M, data=ampegen_avg, correlation=corPagel(OB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OBTS_out<-gls(TS_M~OB_M, data=out_avg, correlation=corPagel(OB_lam,phy=out.tree,form=~Species, fixed=TRUE))
OBTS_gym<-gls(TS_M~OB_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
OBTS_morm<-gls(TS_M~OB_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
OBTS_allgout<-gls(TS_M~OB_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
OBTS_allmout<-gls(TS_M~OB_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
OBTS_gwave<-gls(TS_M~OB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OBTS_gpulse<-gls(TS_M~OB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

OBCB_all<-gls(CB_M~OB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
OB_lam<-as.numeric(OBCB_all$modelStruct)
if (OB_lam > 1)
	OB_lam = 1
OBCB_tubegen<-gls(CB_M~OB_M, data=tubegen_avg, correlation=corPagel(OB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OBCB_ampegen<-gls(CB_M~OB_M, data=ampegen_avg, correlation=corPagel(OB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OBCB_out<-gls(CB_M~OB_M, data=out_avg, correlation=corPagel(OB_lam,phy=out.tree,form=~Species, fixed=TRUE))
OBCB_gym<-gls(CB_M~OB_M, data=gym_avg, correlation=corPagel(0,phy=gym.tree,form=~Species, fixed=TRUE))
OBCB_morm<-gls(CB_M~OB_M, data=morm_avg, correlation=corPagel(0,phy=morm.tree,form=~Species, fixed=TRUE))
OBCB_allgout<-gls(CB_M~OB_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
OBCB_allmout<-gls(CB_M~OB_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
OBCB_gwave<-gls(CB_M~OB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OBCB_gpulse<-gls(CB_M~OB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

OBROB_all<-gls(RoB_M~OB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(1,phy=tree,form=~Species))
OB_lam<-as.numeric(OBROB_all$modelStruct)
if (OB_lam > 1)
	OB_lam = 1
OBROB_tubegen<-gls(RoB_M~OB_M, data=tubegen_avg, correlation=corPagel(OB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OBROB_ampegen<-gls(RoB_M~OB_M, data=ampegen_avg, correlation=corPagel(OB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OBROB_out<-gls(RoB_M~OB_M, data=out_avg, correlation=corPagel(OB_lam,phy=out.tree,form=~Species, fixed=TRUE))
OBROB_gym<-gls(RoB_M~OB_M, data=gym_avg, correlation=corPagel(0.4686382,phy=gym.tree,form=~Species, fixed=TRUE))
OBROB_morm<-gls(RoB_M~OB_M, data=morm_avg, correlation=corPagel(0.4686382,phy=morm.tree,form=~Species, fixed=TRUE))
OBROB_allgout<-gls(RoB_M~OB_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
OBROB_allmout<-gls(RoB_M~OB_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
OBROB_gwave<-gls(RoB_M~OB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OBROB_gpulse<-gls(RoB_M~OB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

## tel as x
TELOB_all<-gls(OB_M~TEL_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
TEL_lam<-as.numeric(TELOB_all$modelStruct)
if (TEL_lam > 1)
	TEL_lam = 1
TELOB_tubegen<-gls(OB_M~TEL_M, data=tubegen_avg, correlation=corPagel(TEL_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TELOB_ampegen<-gls(OB_M~TEL_M, data=ampegen_avg, correlation=corPagel(TEL_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TELOB_out<-gls(OB_M~TEL_M, data=out_avg, correlation=corPagel(TEL_lam,phy=out.tree,form=~Species, fixed=TRUE))
TELOB_gym<-gls(OB_M~TEL_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
TELOB_morm<-gls(OB_M~TEL_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
TELOB_allgout<-gls(OB_M~TEL_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
TELOB_allmout<-gls(OB_M~TEL_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
TELOB_gwave<-gls(OB_M~TEL_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
TELOB_gpulse<-gls(OB_M~TEL_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

TELHB_all<-gls(HB_M~TEL_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.9,phy=tree,form=~Species))
TEL_lam<-as.numeric(TELHB_all$modelStruct)
if (TEL_lam > 1)
	TEL_lam = 1
TELHB_tubegen<-gls(HB_M~TEL_M, data=tubegen_avg, correlation=corPagel(TEL_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TELHB_ampegen<-gls(HB_M~TEL_M, data=ampegen_avg, correlation=corPagel(TEL_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TELHB_out<-gls(HB_M~TEL_M, data=out_avg, correlation=corPagel(TEL_lam,phy=out.tree,form=~Species, fixed=TRUE))
TELHB_gym<-gls(HB_M~TEL_M, data=gym_avg, correlation=corPagel(0.9969506,phy=gym.tree,form=~Species, fixed=TRUE))
TELHB_morm<-gls(HB_M~TEL_M, data=morm_avg, correlation=corPagel(0.9969506,phy=morm.tree,form=~Species, fixed=TRUE))
TELHB_allgout<-gls(HB_M~TEL_M, data=allgout_avg, correlation=corPagel(0.9267392,phy=allgout.tree,form=~Species, fixed=TRUE))
TELHB_allmout<-gls(HB_M~TEL_M, data=allmout_avg, correlation=corPagel(0.9267392,phy=allmout.tree,form=~Species, fixed=TRUE))
TELHB_gwave<-gls(HB_M~TEL_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
TELHB_gpulse<-gls(HB_M~TEL_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

TELOT_all<-gls(OT_M~TEL_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.5,phy=tree,form=~Species))
TEL_lam<-as.numeric(TELOT_all$modelStruct)
if (TEL_lam > 1)
	TEL_lam = 1
TELOT_tubegen<-gls(OT_M~TEL_M, data=tubegen_avg, correlation=corPagel(TEL_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TELOT_ampegen<-gls(OT_M~TEL_M, data=ampegen_avg, correlation=corPagel(TEL_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TELOT_out<-gls(OT_M~TEL_M, data=out_avg, correlation=corPagel(TEL_lam,phy=out.tree,form=~Species, fixed=TRUE))
TELOT_gym<-gls(OT_M~TEL_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
TELOT_morm<-gls(OT_M~TEL_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
TELOT_allgout<-gls(OT_M~TEL_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
TELOT_allmout<-gls(OT_M~TEL_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
TELOT_gwave<-gls(OT_M~TEL_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
TELOT_gpulse<-gls(OT_M~TEL_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

TELTS_all<-gls(TS_M~TEL_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
TEL_lam<-as.numeric(TELTS_all$modelStruct)
if (TEL_lam > 1)
	TEL_lam = 1
TELTS_tubegen<-gls(TS_M~TEL_M, data=tubegen_avg, correlation=corPagel(TEL_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TELTS_ampegen<-gls(TS_M~TEL_M, data=ampegen_avg, correlation=corPagel(TEL_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TELTS_out<-gls(TS_M~TEL_M, data=out_avg, correlation=corPagel(TEL_lam,phy=out.tree,form=~Species, fixed=TRUE))
TELTS_gym<-gls(TS_M~TEL_M, data=gym_avg, correlation=corPagel(0.7487472,phy=gym.tree,form=~Species, fixed=TRUE))
TELTS_morm<-gls(TS_M~TEL_M, data=morm_avg, correlation=corPagel(0.7487472,phy=morm.tree,form=~Species, fixed=TRUE))
TELTS_allgout<-gls(TS_M~TEL_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
TELTS_allmout<-gls(TS_M~TEL_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
TELTS_gwave<-gls(TS_M~TEL_M, data=gwave_avg, correlation=corPagel(0.4636762,phy=gwave.tree,form=~Species, fixed=TRUE))
TELTS_gpulse<-gls(TS_M~TEL_M, data=gpulse_avg, correlation=corPagel(0.4636762,phy=gpulse.tree,form=~Species, fixed=TRUE))

TELCB_all<-gls(CB_M~TEL_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
TEL_lam<-as.numeric(TELCB_all$modelStruct)
if (TEL_lam > 1)
	TEL_lam = 1
TELCB_tubegen<-gls(CB_M~TEL_M, data=tubegen_avg, correlation=corPagel(TEL_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TELCB_ampegen<-gls(CB_M~TEL_M, data=ampegen_avg, correlation=corPagel(TEL_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TELCB_out<-gls(CB_M~TEL_M, data=out_avg, correlation=corPagel(TEL_lam,phy=out.tree,form=~Species, fixed=TRUE))
TELCB_gym<-gls(CB_M~TEL_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
TELCB_morm<-gls(CB_M~TEL_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
TELCB_allgout<-gls(CB_M~TEL_M, data=allgout_avg, correlation=corPagel(0.01939753,phy=allgout.tree,form=~Species, fixed=TRUE))
TELCB_allmout<-gls(CB_M~TEL_M, data=allmout_avg, correlation=corPagel(0.01939753,phy=allmout.tree,form=~Species, fixed=TRUE))
TELCB_gwave<-gls(CB_M~TEL_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
TELCB_gpulse<-gls(CB_M~TEL_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

TELROB_all<-gls(RoB_M~TEL_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
TEL_lam<-as.numeric(TELROB_all$modelStruct)
if (TEL_lam > 1)
	TEL_lam = 1
TELROB_tubegen<-gls(RoB_M~TEL_M, data=tubegen_avg, correlation=corPagel(TEL_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TELROB_ampegen<-gls(RoB_M~TEL_M, data=ampegen_avg, correlation=corPagel(TEL_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TELROB_out<-gls(RoB_M~TEL_M, data=out_avg, correlation=corPagel(TEL_lam,phy=out.tree,form=~Species, fixed=TRUE))
TELROB_gym<-gls(RoB_M~TEL_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
TELROB_morm<-gls(RoB_M~TEL_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
TELROB_allgout<-gls(RoB_M~TEL_M, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
TELROB_allmout<-gls(RoB_M~TEL_M, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
TELROB_gwave<-gls(RoB_M~TEL_M, data=gwave_avg, correlation=corPagel(0.3419171,phy=gwave.tree,form=~Species, fixed=TRUE))
TELROB_gpulse<-gls(RoB_M~TEL_M, data=gpulse_avg, correlation=corPagel(0.3419171,phy=gpulse.tree,form=~Species, fixed=TRUE))

## HB as x
HBOB_all<-gls(OB_M~HB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
HB_lam<-as.numeric(HBOB_all$modelStruct)
if (HB_lam > 1)
	HB_lam = 1
HBOB_tubegen<-gls(OB_M~HB_M, data=tubegen_avg, correlation=corPagel(HB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
HBOB_ampegen<-gls(OB_M~HB_M, data=ampegen_avg, correlation=corPagel(HB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
HBOB_out<-gls(OB_M~HB_M, data=out_avg, correlation=corPagel(HB_lam,phy=out.tree,form=~Species, fixed=TRUE))
HBOB_gym<-gls(OB_M~HB_M, data=gym_avg, correlation=corPagel(0.8590364,phy=gym.tree,form=~Species, fixed=TRUE))
HBOB_morm<-gls(OB_M~HB_M, data=morm_avg, correlation=corPagel(0.8590364,phy=morm.tree,form=~Species, fixed=TRUE))
HBOB_allgout<-gls(OB_M~HB_M, data=allgout_avg, correlation=corPagel(0.3555335,phy=allgout.tree,form=~Species, fixed=TRUE))
HBOB_allmout<-gls(OB_M~HB_M, data=allmout_avg, correlation=corPagel(0.3555335,phy=allmout.tree,form=~Species, fixed=TRUE))
HBOB_gwave<-gls(OB_M~HB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
HBOB_gpulse<-gls(OB_M~HB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

HBTEL_all<-gls(TEL_M~HB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
HB_lam<-as.numeric(HBTEL_all$modelStruct)
if (HB_lam > 1)
	HB_lam = 1
HBTEL_tubegen<-gls(TEL_M~HB_M, data=tubegen_avg, correlation=corPagel(HB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
HBTEL_ampegen<-gls(TEL_M~HB_M, data=ampegen_avg, correlation=corPagel(HB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
HBTEL_out<-gls(TEL_M~HB_M, data=out_avg, correlation=corPagel(HB_lam,phy=out.tree,form=~Species, fixed=TRUE))
HBTEL_gym<-gls(TEL_M~HB_M, data=gym_avg, correlation=corPagel(0.8744299,phy=gym.tree,form=~Species, fixed=TRUE))
HBTEL_morm<-gls(TEL_M~HB_M, data=morm_avg, correlation=corPagel(0.8744299,phy=morm.tree,form=~Species, fixed=TRUE))
HBTEL_allgout<-gls(TEL_M~HB_M, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
HBTEL_allmout<-gls(TEL_M~HB_M, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
HBTEL_gwave<-gls(TEL_M~HB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
HBTEL_gpulse<-gls(TEL_M~HB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

HBOT_all<-gls(OT_M~HB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
HB_lam<-as.numeric(HBOT_all$modelStruct)
if (HB_lam > 1)
	HB_lam = 1
HBOT_tubegen<-gls(OT_M~HB_M, data=tubegen_avg, correlation=corPagel(HB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
HBOT_ampegen<-gls(OT_M~HB_M, data=ampegen_avg, correlation=corPagel(HB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
HBOT_out<-gls(OT_M~HB_M, data=out_avg, correlation=corPagel(HB_lam,phy=out.tree,form=~Species, fixed=TRUE))
HBOT_gym<-gls(OT_M~HB_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
HBOT_morm<-gls(OT_M~HB_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
HBOT_allgout<-gls(OT_M~HB_M, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
HBOT_allmout<-gls(OT_M~HB_M, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
HBOT_gwave<-gls(OT_M~HB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
HBOT_gpulse<-gls(OT_M~HB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

HBTS_all<-gls(TS_M~HB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
HB_lam<-as.numeric(HBTS_all$modelStruct)
if (HB_lam > 1)
	HB_lam = 1
HBTS_tubegen<-gls(TS_M~HB_M, data=tubegen_avg, correlation=corPagel(HB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
HBTS_ampegen<-gls(TS_M~HB_M, data=ampegen_avg, correlation=corPagel(HB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
HBTS_out<-gls(TS_M~HB_M, data=out_avg, correlation=corPagel(HB_lam,phy=out.tree,form=~Species, fixed=TRUE))
HBTS_gym<-gls(TS_M~HB_M, data=gym_avg, correlation=corPagel(0.557723,phy=gym.tree,form=~Species, fixed=TRUE))
HBTS_morm<-gls(TS_M~HB_M, data=morm_avg, correlation=corPagel(0.557723,phy=morm.tree,form=~Species, fixed=TRUE))
HBTS_allgout<-gls(TS_M~HB_M, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
HBTS_allmout<-gls(TS_M~HB_M, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
HBTS_gwave<-gls(TS_M~HB_M, data=gwave_avg, correlation=corPagel(0.5832345,phy=gwave.tree,form=~Species, fixed=TRUE))
HBTS_gpulse<-gls(TS_M~HB_M, data=gpulse_avg, correlation=corPagel(0.5832345,phy=gpulse.tree,form=~Species, fixed=TRUE))

HBCB_all<-gls(CB_M~HB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
HB_lam<-as.numeric(HBCB_all$modelStruct)
if (HB_lam > 1)
	HB_lam = 1
HBCB_tubegen<-gls(CB_M~HB_M, data=tubegen_avg, correlation=corPagel(HB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
HBCB_ampegen<-gls(CB_M~HB_M, data=ampegen_avg, correlation=corPagel(HB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
HBCB_out<-gls(CB_M~HB_M, data=out_avg, correlation=corPagel(HB_lam,phy=out.tree,form=~Species, fixed=TRUE))
HBCB_gym<-gls(CB_M~HB_M, data=gym_avg, correlation=corPagel(0.7930674,phy=gym.tree,form=~Species, fixed=TRUE))
HBCB_morm<-gls(CB_M~HB_M, data=morm_avg, correlation=corPagel(0.7930674,phy=morm.tree,form=~Species, fixed=TRUE))
HBCB_allgout<-gls(CB_M~HB_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
HBCB_allmout<-gls(CB_M~HB_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
HBCB_gwave<-gls(CB_M~HB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
HBCB_gpulse<-gls(CB_M~HB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

HBROB_all<-gls(RoB_M~HB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
HB_lam<-as.numeric(HBROB_all$modelStruct)
if (HB_lam > 1)
	HB_lam = 1
HBROB_tubegen<-gls(RoB_M~HB_M, data=tubegen_avg, correlation=corPagel(HB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
HBROB_ampegen<-gls(RoB_M~HB_M, data=ampegen_avg, correlation=corPagel(HB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
HBROB_out<-gls(RoB_M~HB_M, data=out_avg, correlation=corPagel(HB_lam,phy=out.tree,form=~Species, fixed=TRUE))
HBROB_gym<-gls(RoB_M~HB_M, data=gym_avg, correlation=corPagel(0.5776069,phy=gym.tree,form=~Species, fixed=TRUE))
HBROB_morm<-gls(RoB_M~HB_M, data=morm_avg, correlation=corPagel(0.5776069,phy=morm.tree,form=~Species, fixed=TRUE))
HBROB_allgout<-gls(RoB_M~HB_M, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
HBROB_allmout<-gls(RoB_M~HB_M, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
HBROB_gwave<-gls(RoB_M~HB_M, data=gwave_avg, correlation=corPagel(0,phy=gwave.tree,form=~Species, fixed=TRUE))
HBROB_gpulse<-gls(RoB_M~HB_M, data=gpulse_avg, correlation=corPagel(0,phy=gpulse.tree,form=~Species, fixed=TRUE))

## ot as x
OTOB_all<-gls(OB_M~OT_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
OT_lam<-as.numeric(OTOB_all$modelStruct)
if (OT_lam > 1)
	OT_lam = 1
OTOB_tubegen<-gls(OB_M~OT_M, data=tubegen_avg, correlation=corPagel(OT_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OTOB_ampegen<-gls(OB_M~OT_M, data=ampegen_avg, correlation=corPagel(OT_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OTOB_out<-gls(OB_M~OT_M, data=out_avg, correlation=corPagel(OT_lam,phy=out.tree,form=~Species, fixed=TRUE))
OTOB_gym<-gls(OB_M~OT_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
OTOB_morm<-gls(OB_M~OT_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
OTOB_allgout<-gls(OB_M~OT_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
OTOB_allmout<-gls(OB_M~OT_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
OTOB_gwave<-gls(OB_M~OT_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OTOB_gpulse<-gls(OB_M~OT_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

OTTEL_all<-gls(TEL_M~OT_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
OT_lam<-as.numeric(OTTEL_all$modelStruct)
if (OT_lam > 1)
	OT_lam = 1
OTTEL_tubegen<-gls(TEL_M~OT_M, data=tubegen_avg, correlation=corPagel(OT_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OTTEL_ampegen<-gls(TEL_M~OT_M, data=ampegen_avg, correlation=corPagel(OT_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OTTEL_out<-gls(TEL_M~OT_M, data=out_avg, correlation=corPagel(OT_lam,phy=out.tree,form=~Species, fixed=TRUE))
OTTEL_gym<-gls(TEL_M~OT_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
OTTEL_morm<-gls(TEL_M~OT_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
OTTEL_allgout<-gls(TEL_M~OT_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
OTTEL_allmout<-gls(TEL_M~OT_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
OTTEL_gwave<-gls(TEL_M~OT_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OTTEL_gpulse<-gls(TEL_M~OT_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

OTHB_all<-gls(HB_M~OT_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
OT_lam<-as.numeric(OTHB_all$modelStruct)
if (OT_lam > 1)
	OT_lam = 1
OTHB_tubegen<-gls(HB_M~OT_M, data=tubegen_avg, correlation=corPagel(OT_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OTHB_ampegen<-gls(HB_M~OT_M, data=ampegen_avg, correlation=corPagel(OT_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OTHB_out<-gls(HB_M~OT_M, data=out_avg, correlation=corPagel(OT_lam,phy=out.tree,form=~Species, fixed=TRUE))
OTHB_gym<-gls(HB_M~OT_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
OTHB_morm<-gls(HB_M~OT_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
OTHB_allgout<-gls(HB_M~OT_M, data=allgout_avg, correlation=corPagel(0.0765538,phy=allgout.tree,form=~Species, fixed=TRUE))
OTHB_allmout<-gls(HB_M~OT_M, data=allmout_avg, correlation=corPagel(0.0765538,phy=allmout.tree,form=~Species, fixed=TRUE))
OTHB_gwave<-gls(HB_M~OT_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OTHB_gpulse<-gls(HB_M~OT_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

OTTS_all<-gls(TS_M~OT_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
OT_lam<-as.numeric(OTTS_all$modelStruct)
if (OT_lam > 1)
	OT_lam = 1
OTTS_tubegen<-gls(TS_M~OT_M, data=tubegen_avg, correlation=corPagel(OT_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OTTS_ampegen<-gls(TS_M~OT_M, data=ampegen_avg, correlation=corPagel(OT_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OTTS_out<-gls(TS_M~OT_M, data=out_avg, correlation=corPagel(OT_lam,phy=out.tree,form=~Species, fixed=TRUE))
OTTS_gym<-gls(TS_M~OT_M, data=gym_avg, correlation=corPagel(0.9111959,phy=gym.tree,form=~Species, fixed=TRUE))
OTTS_morm<-gls(TS_M~OT_M, data=morm_avg, correlation=corPagel(0.9111959,phy=morm.tree,form=~Species, fixed=TRUE))
OTTS_allgout<-gls(TS_M~OT_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
OTTS_allmout<-gls(TS_M~OT_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
OTTS_gwave<-gls(TS_M~OT_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OTTS_gpulse<-gls(TS_M~OT_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

OTCB_all<-gls(CB_M~OT_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
OT_lam<-as.numeric(OTCB_all$modelStruct)
if (OT_lam > 1)
	OT_lam = 1
OTCB_tubegen<-gls(CB_M~OT_M, data=tubegen_avg, correlation=corPagel(OT_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OTCB_ampegen<-gls(CB_M~OT_M, data=ampegen_avg, correlation=corPagel(OT_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OTCB_out<-gls(CB_M~OT_M, data=out_avg, correlation=corPagel(OT_lam,phy=out.tree,form=~Species, fixed=TRUE))
OTCB_gym<-gls(CB_M~OT_M, data=gym_avg, correlation=corPagel(0.9951343,phy=gym.tree,form=~Species, fixed=TRUE))
OTCB_morm<-gls(CB_M~OT_M, data=morm_avg, correlation=corPagel(0.9951343,phy=morm.tree,form=~Species, fixed=TRUE))
OTCB_allgout<-gls(CB_M~OT_M, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
OTCB_allmout<-gls(CB_M~OT_M, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
OTCB_gwave<-gls(CB_M~OT_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OTCB_gpulse<-gls(CB_M~OT_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

OTROB_all<-gls(RoB_M~OT_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
OT_lam<-as.numeric(OTROB_all$modelStruct)
if (OT_lam > 1)
	OT_lam = 1
OTROB_tubegen<-gls(RoB_M~OT_M, data=tubegen_avg, correlation=corPagel(OT_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
OTROB_ampegen<-gls(RoB_M~OT_M, data=ampegen_avg, correlation=corPagel(OT_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
OTROB_out<-gls(RoB_M~OT_M, data=out_avg, correlation=corPagel(OT_lam,phy=out.tree,form=~Species, fixed=TRUE))
OTROB_gym<-gls(RoB_M~OT_M, data=gym_avg, correlation=corPagel(0.9269965,phy=gym.tree,form=~Species, fixed=TRUE))
OTROB_morm<-gls(RoB_M~OT_M, data=morm_avg, correlation=corPagel(0.9269965,phy=morm.tree,form=~Species, fixed=TRUE))
OTROB_allgout<-gls(RoB_M~OT_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
OTROB_allmout<-gls(RoB_M~OT_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
OTROB_gwave<-gls(RoB_M~OT_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
OTROB_gpulse<-gls(RoB_M~OT_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

## ts as x
TSOB_all<-gls(OB_M~TS_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
TS_lam<-as.numeric(TSOB_all$modelStruct)
if (TS_lam > 1)
	TS_lam = 1
TSOB_tubegen<-gls(OB_M~TS_M, data=tubegen_avg, correlation=corPagel(TS_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TSOB_ampegen<-gls(OB_M~TS_M, data=ampegen_avg, correlation=corPagel(TS_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TSOB_out<-gls(OB_M~TS_M, data=out_avg, correlation=corPagel(TS_lam,phy=out.tree,form=~Species, fixed=TRUE))
TSOB_gym<-gls(OB_M~TS_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
TSOB_morm<-gls(OB_M~TS_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
TSOB_allgout<-gls(OB_M~TS_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
TSOB_allmout<-gls(OB_M~TS_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
TSOB_gwave<-gls(OB_M~TS_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
TSOB_gpulse<-gls(OB_M~TS_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

TSTEL_all<-gls(TEL_M~TS_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
TS_lam<-as.numeric(TSTEL_all$modelStruct)
if (TS_lam > 1)
	TS_lam = 1
TSTEL_tubegen<-gls(TEL_M~TS_M, data=tubegen_avg, correlation=corPagel(TS_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TSTEL_ampegen<-gls(TEL_M~TS_M, data=ampegen_avg, correlation=corPagel(TS_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TSTEL_out<-gls(TEL_M~TS_M, data=out_avg, correlation=corPagel(TS_lam,phy=out.tree,form=~Species, fixed=TRUE))
TSTEL_gym<-gls(TEL_M~TS_M, data=gym_avg, correlation=corPagel(0.8081982,phy=gym.tree,form=~Species, fixed=TRUE))
TSTEL_morm<-gls(TEL_M~TS_M, data=morm_avg, correlation=corPagel(0.8081982,phy=morm.tree,form=~Species, fixed=TRUE))
TSTEL_allgout<-gls(TEL_M~TS_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
TSTEL_allmout<-gls(TEL_M~TS_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
TSTEL_gwave<-gls(TEL_M~TS_M, data=gwave_avg, correlation=corPagel(0.6656835,phy=gwave.tree,form=~Species, fixed=TRUE))
TSTEL_gpulse<-gls(TEL_M~TS_M, data=gpulse_avg, correlation=corPagel(0.6656835,phy=gpulse.tree,form=~Species, fixed=TRUE))

TSHB_all<-gls(HB_M~TS_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
TS_lam<-as.numeric(TSHB_all$modelStruct)
if (TS_lam > 1)
	TS_lam = 1
TSHB_tubegen<-gls(HB_M~TS_M, data=tubegen_avg, correlation=corPagel(TS_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TSHB_ampegen<-gls(HB_M~TS_M, data=ampegen_avg, correlation=corPagel(TS_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TSHB_out<-gls(HB_M~TS_M, data=out_avg, correlation=corPagel(TS_lam,phy=out.tree,form=~Species, fixed=TRUE))
TSHB_gym<-gls(HB_M~TS_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
TSHB_morm<-gls(HB_M~TS_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
TSHB_allgout<-gls(HB_M~TS_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
TSHB_allmout<-gls(HB_M~TS_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
TSHB_gwave<-gls(HB_M~TS_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
TSHB_gpulse<-gls(HB_M~TS_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

TSOT_all<-gls(OT_M~TS_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
TS_lam<-as.numeric(TSOT_all$modelStruct)
if (TS_lam > 1)
	TS_lam = 1
TSOT_tubegen<-gls(OT_M~TS_M, data=tubegen_avg, correlation=corPagel(TS_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TSOT_ampegen<-gls(OT_M~TS_M, data=ampegen_avg, correlation=corPagel(TS_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TSOT_out<-gls(OT_M~TS_M, data=out_avg, correlation=corPagel(TS_lam,phy=out.tree,form=~Species, fixed=TRUE))
TSOT_gym<-gls(OT_M~TS_M, data=gym_avg, correlation=corPagel(0.8933667,phy=gym.tree,form=~Species, fixed=TRUE))
TSOT_morm<-gls(OT_M~TS_M, data=morm_avg, correlation=corPagel(0.8933667,phy=morm.tree,form=~Species, fixed=TRUE))
TSOT_allgout<-gls(OT_M~TS_M, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
TSOT_allmout<-gls(OT_M~TS_M, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
TSOT_gwave<-gls(OT_M~TS_M, data=gwave_avg, correlation=corPagel(0.9685348,phy=gwave.tree,form=~Species, fixed=TRUE))
TSOT_gpulse<-gls(OT_M~TS_M, data=gpulse_avg, correlation=corPagel(0.9685348,phy=gpulse.tree,form=~Species, fixed=TRUE))

TSCB_all<-gls(CB_M~TS_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
TS_lam<-as.numeric(TSCB_all$modelStruct)
if (TS_lam > 1)
	TS_lam = 1
TSCB_tubegen<-gls(CB_M~TS_M, data=tubegen_avg, correlation=corPagel(TS_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TSCB_ampegen<-gls(CB_M~TS_M, data=ampegen_avg, correlation=corPagel(TS_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TSCB_out<-gls(CB_M~TS_M, data=out_avg, correlation=corPagel(TS_lam,phy=out.tree,form=~Species, fixed=TRUE))
TSCB_gym<-gls(CB_M~TS_M, data=gym_avg, correlation=corPagel(0.6803758,phy=gym.tree,form=~Species, fixed=TRUE))
TSCB_morm<-gls(CB_M~TS_M, data=morm_avg, correlation=corPagel(0.6803758,phy=morm.tree,form=~Species, fixed=TRUE))
TSCB_allgout<-gls(CB_M~TS_M, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
TSCB_allmout<-gls(CB_M~TS_M, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
TSCB_gwave<-gls(CB_M~TS_M, data=gwave_avg, correlation=corPagel(0.3071866,phy=gwave.tree,form=~Species, fixed=TRUE))
TSCB_gpulse<-gls(CB_M~TS_M, data=gpulse_avg, correlation=corPagel(0.3071866,phy=gpulse.tree,form=~Species, fixed=TRUE))

TSROB_all<-gls(RoB_M~TS_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
TS_lam<-as.numeric(TSROB_all$modelStruct)
if (TS_lam > 1)
	TS_lam = 1
TSROB_tubegen<-gls(RoB_M~TS_M, data=tubegen_avg, correlation=corPagel(TS_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
TSROB_ampegen<-gls(RoB_M~TS_M, data=ampegen_avg, correlation=corPagel(TS_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
TSROB_out<-gls(RoB_M~TS_M, data=out_avg, correlation=corPagel(TS_lam,phy=out.tree,form=~Species, fixed=TRUE))
TSROB_gym<-gls(RoB_M~TS_M, data=gym_avg, correlation=corPagel(0.8475078,phy=gym.tree,form=~Species, fixed=TRUE))
TSROB_morm<-gls(RoB_M~TS_M, data=morm_avg, correlation=corPagel(0.8475078,phy=morm.tree,form=~Species, fixed=TRUE))
TSROB_allgout<-gls(RoB_M~TS_M, data=allgout_avg, correlation=corPagel(0.9951045,phy=allgout.tree,form=~Species, fixed=TRUE))
TSROB_allmout<-gls(RoB_M~TS_M, data=allmout_avg, correlation=corPagel(0.9951045,phy=allmout.tree,form=~Species, fixed=TRUE))
TSROB_gwave<-gls(RoB_M~TS_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
TSROB_gpulse<-gls(RoB_M~TS_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

## cb as x
CBOB_all<-gls(OB_M~CB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
CB_lam<-as.numeric(CBOB_all$modelStruct)
if (CB_lam > 1)
	CB_lam = 1
CBOB_tubegen<-gls(OB_M~CB_M, data=tubegen_avg, correlation=corPagel(CB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
CBOB_ampegen<-gls(OB_M~CB_M, data=ampegen_avg, correlation=corPagel(CB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
CBOB_out<-gls(OB_M~CB_M, data=out_avg, correlation=corPagel(CB_lam,phy=out.tree,form=~Species, fixed=TRUE))
CBOB_gym<-gls(OB_M~CB_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
CBOB_morm<-gls(OB_M~CB_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
CBOB_allgout<-gls(OB_M~CB_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
CBOB_allmout<-gls(OB_M~CB_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
CBOB_gwave<-gls(OB_M~CB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
CBOB_gpulse<-gls(OB_M~CB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

CBTEL_all<-gls(TEL_M~CB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
CB_lam<-as.numeric(CBTEL_all$modelStruct)
if (CB_lam > 1)
	CB_lam = 1
CBTEL_tubegen<-gls(TEL_M~CB_M, data=tubegen_avg, correlation=corPagel(CB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
CBTEL_ampegen<-gls(TEL_M~CB_M, data=ampegen_avg, correlation=corPagel(CB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
CBTEL_out<-gls(TEL_M~CB_M, data=out_avg, correlation=corPagel(CB_lam,phy=out.tree,form=~Species, fixed=TRUE))
CBTEL_gym<-gls(TEL_M~CB_M, data=gym_avg, correlation=corPagel(0.9407809,phy=gym.tree,form=~Species, fixed=TRUE))
CBTEL_morm<-gls(TEL_M~CB_M, data=morm_avg, correlation=corPagel(0.9407809,phy=morm.tree,form=~Species, fixed=TRUE))
CBTEL_allgout<-gls(TEL_M~CB_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
CBTEL_allmout<-gls(TEL_M~CB_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
CBTEL_gwave<-gls(TEL_M~CB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
CBTEL_gpulse<-gls(TEL_M~CB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

CBHB_all<-gls(HB_M~CB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
CB_lam<-as.numeric(CBHB_all$modelStruct)
if (CB_lam > 1)
	CB_lam = 1
CBHB_tubegen<-gls(HB_M~CB_M, data=tubegen_avg, correlation=corPagel(CB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
CBHB_ampegen<-gls(HB_M~CB_M, data=ampegen_avg, correlation=corPagel(CB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
CBHB_out<-gls(HB_M~CB_M, data=out_avg, correlation=corPagel(CB_lam,phy=out.tree,form=~Species, fixed=TRUE))
CBHB_gym<-gls(HB_M~CB_M, data=gym_avg, correlation=corPagel(0.9838753,phy=gym.tree,form=~Species, fixed=TRUE))
CBHB_morm<-gls(HB_M~CB_M, data=morm_avg, correlation=corPagel(0.9838753,phy=morm.tree,form=~Species, fixed=TRUE))
CBHB_allgout<-gls(HB_M~CB_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
CBHB_allmout<-gls(HB_M~CB_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
CBHB_gwave<-gls(HB_M~CB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
CBHB_gpulse<-gls(HB_M~CB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

CBOT_all<-gls(OT_M~CB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
CB_lam<-as.numeric(CBOT_all$modelStruct)
if (CB_lam > 1)
	CB_lam = 1
CBOT_tubegen<-gls(OT_M~CB_M, data=tubegen_avg, correlation=corPagel(CB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
CBOT_ampegen<-gls(OT_M~CB_M, data=ampegen_avg, correlation=corPagel(CB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
CBOT_out<-gls(OT_M~CB_M, data=out_avg, correlation=corPagel(CB_lam,phy=out.tree,form=~Species, fixed=TRUE))
CBOT_gym<-gls(OT_M~CB_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
CBOT_morm<-gls(OT_M~CB_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
CBOT_allgout<-gls(OT_M~CB_M, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
CBOT_allmout<-gls(OT_M~CB_M, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
CBOT_gwave<-gls(OT_M~CB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
CBOT_gpulse<-gls(OT_M~CB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

CBTS_all<-gls(TS_M~CB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
CB_lam<-as.numeric(CBTS_all$modelStruct)
if (CB_lam > 1)
	CB_lam = 1
CBTS_tubegen<-gls(TS_M~CB_M, data=tubegen_avg, correlation=corPagel(CB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
CBTS_ampegen<-gls(TS_M~CB_M, data=ampegen_avg, correlation=corPagel(CB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
CBTS_out<-gls(TS_M~CB_M, data=out_avg, correlation=corPagel(CB_lam,phy=out.tree,form=~Species, fixed=TRUE))
CBTS_gym<-gls(TS_M~CB_M, data=gym_avg, correlation=corPagel(0.7379803,phy=gym.tree,form=~Species, fixed=TRUE))
CBTS_morm<-gls(TS_M~CB_M, data=morm_avg, correlation=corPagel(0.7379803,phy=morm.tree,form=~Species, fixed=TRUE))
CBTS_allgout<-gls(TS_M~CB_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
CBTS_allmout<-gls(TS_M~CB_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
CBTS_gwave<-gls(TS_M~CB_M, data=gwave_avg, correlation=corPagel(0,phy=gwave.tree,form=~Species, fixed=TRUE))
CBTS_gpulse<-gls(TS_M~CB_M, data=gpulse_avg, correlation=corPagel(0,phy=gpulse.tree,form=~Species, fixed=TRUE))

CBROB_all<-gls(RoB_M~CB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
CB_lam<-as.numeric(CBROB_all$modelStruct)
if (CB_lam > 1)
	CB_lam = 1
CBROB_tubegen<-gls(RoB_M~CB_M, data=tubegen_avg, correlation=corPagel(CB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
CBROB_ampegen<-gls(RoB_M~CB_M, data=ampegen_avg, correlation=corPagel(CB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
CBROB_out<-gls(RoB_M~CB_M, data=out_avg, correlation=corPagel(CB_lam,phy=out.tree,form=~Species, fixed=TRUE))
CBROB_gym<-gls(RoB_M~CB_M, data=gym_avg, correlation=corPagel(0.833309,phy=gym.tree,form=~Species, fixed=TRUE))
CBROB_morm<-gls(RoB_M~CB_M, data=morm_avg, correlation=corPagel(0.833309,phy=morm.tree,form=~Species, fixed=TRUE))
CBROB_allgout<-gls(RoB_M~CB_M, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
CBROB_allmout<-gls(RoB_M~CB_M, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
CBROB_gwave<-gls(RoB_M~CB_M, data=gwave_avg, correlation=corPagel(0,phy=gwave.tree,form=~Species, fixed=TRUE))
CBROB_gpulse<-gls(RoB_M~CB_M, data=gpulse_avg, correlation=corPagel(0,phy=gpulse.tree,form=~Species, fixed=TRUE))

## rob as x
ROBTEL_all<-gls(TEL_M~RoB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
ROB_lam<-as.numeric(ROBTEL_all$modelStruct)
if (ROB_lam > 1)
	ROB_lam = 1
ROBTEL_tubegen<-gls(TEL_M~RoB_M, data=tubegen_avg, correlation=corPagel(ROB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
ROBTEL_ampegen<-gls(TEL_M~RoB_M, data=ampegen_avg, correlation=corPagel(ROB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
ROBTEL_out<-gls(TEL_M~RoB_M, data=out_avg, correlation=corPagel(ROB_lam,phy=out.tree,form=~Species, fixed=TRUE))
ROBTEL_gym<-gls(TEL_M~RoB_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
ROBTEL_morm<-gls(TEL_M~RoB_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
ROBTEL_allgout<-gls(TEL_M~RoB_M, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
ROBTEL_allmout<-gls(TEL_M~RoB_M, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
ROBTEL_gwave<-gls(TEL_M~RoB_M, data=gwave_avg, correlation=corPagel(0.3546719,phy=gwave.tree,form=~Species, fixed=TRUE))
ROBTEL_gpulse<-gls(TEL_M~RoB_M, data=gpulse_avg, correlation=corPagel(0.3546719,phy=gpulse.tree,form=~Species, fixed=TRUE))

ROBOB_all<-gls(OB_M~RoB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
ROB_lam<-as.numeric(ROBOB_all$modelStruct)
if (ROB_lam > 1)
	ROB_lam = 1
ROBOB_tubegen<-gls(OB_M~RoB_M, data=tubegen_avg, correlation=corPagel(ROB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
ROBOB_ampegen<-gls(OB_M~RoB_M, data=ampegen_avg, correlation=corPagel(ROB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
ROBOB_out<-gls(OB_M~RoB_M, data=out_avg, correlation=corPagel(ROB_lam,phy=out.tree,form=~Species, fixed=TRUE))
ROBOB_gym<-gls(OB_M~RoB_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
ROBOB_morm<-gls(OB_M~RoB_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
ROBOB_allgout<-gls(OB_M~RoB_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
ROBOB_allmout<-gls(OB_M~RoB_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
ROBOB_gwave<-gls(OB_M~RoB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
ROBOB_gpulse<-gls(OB_M~RoB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

ROBHB_all<-gls(HB_M~RoB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
ROB_lam<-as.numeric(ROBHB_all$modelStruct)
if (ROB_lam > 1)
	ROB_lam = 1
ROBHB_tubegen<-gls(HB_M~RoB_M, data=tubegen_avg, correlation=corPagel(ROB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
ROBHB_ampegen<-gls(HB_M~RoB_M, data=ampegen_avg, correlation=corPagel(ROB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
ROBHB_out<-gls(HB_M~RoB_M, data=out_avg, correlation=corPagel(ROB_lam,phy=out.tree,form=~Species, fixed=TRUE))
ROBHB_gym<-gls(HB_M~RoB_M, data=gym_avg, correlation=corPagel(0.7079042,phy=gym.tree,form=~Species, fixed=TRUE))
ROBHB_morm<-gls(HB_M~RoB_M, data=morm_avg, correlation=corPagel(0.7079042,phy=morm.tree,form=~Species, fixed=TRUE))
ROBHB_allgout<-gls(HB_M~RoB_M, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
ROBHB_allmout<-gls(HB_M~RoB_M, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
ROBHB_gwave<-gls(HB_M~RoB_M, data=gwave_avg, correlation=corPagel(0.3768407,phy=gwave.tree,form=~Species, fixed=TRUE))
ROBHB_gpulse<-gls(HB_M~RoB_M, data=gpulse_avg, correlation=corPagel(0.3768407,phy=gpulse.tree,form=~Species, fixed=TRUE))

ROBOT_all<-gls(OT_M~RoB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.9,phy=tree,form=~Species))
ROB_lam<-as.numeric(ROBOT_all$modelStruct)
if (ROB_lam > 1)
	ROB_lam = 1
ROBOT_tubegen<-gls(OT_M~RoB_M, data=tubegen_avg, correlation=corPagel(ROB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
ROBOT_ampegen<-gls(OT_M~RoB_M, data=ampegen_avg, correlation=corPagel(ROB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
ROBOT_out<-gls(OT_M~RoB_M, data=out_avg, correlation=corPagel(ROB_lam,phy=out.tree,form=~Species, fixed=TRUE))
ROBOT_gym<-gls(OT_M~RoB_M, data=gym_avg, correlation=corPagel(1,phy=gym.tree,form=~Species, fixed=TRUE))
ROBOT_morm<-gls(OT_M~RoB_M, data=morm_avg, correlation=corPagel(1,phy=morm.tree,form=~Species, fixed=TRUE))
ROBOT_allgout<-gls(OT_M~RoB_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
ROBOT_allmout<-gls(OT_M~RoB_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
ROBOT_gwave<-gls(OT_M~RoB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
ROBOT_gpulse<-gls(OT_M~RoB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

ROBTS_all<-gls(TS_M~RoB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
ROB_lam<-as.numeric(ROBTS_all$modelStruct)
if (ROB_lam > 1)
	ROB_lam = 1
ROBTS_tubegen<-gls(TS_M~RoB_M, data=tubegen_avg, correlation=corPagel(ROB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
ROBTS_ampegen<-gls(TS_M~RoB_M, data=ampegen_avg, correlation=corPagel(ROB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
ROBTS_out<-gls(TS_M~RoB_M, data=out_avg, correlation=corPagel(ROB_lam,phy=out.tree,form=~Species, fixed=TRUE))
ROBTS_gym<-gls(TS_M~RoB_M, data=gym_avg, correlation=corPagel(0.7898533,phy=gym.tree,form=~Species, fixed=TRUE))
ROBTS_morm<-gls(TS_M~RoB_M, data=morm_avg, correlation=corPagel(0.7898533,phy=morm.tree,form=~Species, fixed=TRUE))
ROBTS_allgout<-gls(TS_M~RoB_M, data=allgout_avg, correlation=corPagel(0,phy=allgout.tree,form=~Species, fixed=TRUE))
ROBTS_allmout<-gls(TS_M~RoB_M, data=allmout_avg, correlation=corPagel(0,phy=allmout.tree,form=~Species, fixed=TRUE))
ROBTS_gwave<-gls(TS_M~RoB_M, data=gwave_avg, correlation=corPagel(1,phy=gwave.tree,form=~Species, fixed=TRUE))
ROBTS_gpulse<-gls(TS_M~RoB_M, data=gpulse_avg, correlation=corPagel(1,phy=gpulse.tree,form=~Species, fixed=TRUE))

ROBCB_all<-gls(CB_M~RoB_M+hyp_egenrecpt, data=avg_comb_data, correlation=corPagel(.8,phy=tree,form=~Species))
ROB_lam<-as.numeric(ROBCB_all$modelStruct)
if (ROB_lam > 1)
	ROB_lam = 1
ROBCB_tubegen<-gls(CB_M~RoB_M, data=tubegen_avg, correlation=corPagel(ROB_lam,phy=tubegen.tree,form=~Species, fixed=TRUE))
ROBCB_ampegen<-gls(CB_M~RoB_M, data=ampegen_avg, correlation=corPagel(ROB_lam,phy=ampegen.tree,form=~Species, fixed=TRUE))
ROBCB_out<-gls(CB_M~RoB_M, data=out_avg, correlation=corPagel(ROB_lam,phy=out.tree,form=~Species, fixed=TRUE))
ROBCB_gym<-gls(CB_M~RoB_M, data=gym_avg, correlation=corPagel(0.9278729,phy=gym.tree,form=~Species, fixed=TRUE))
ROBCB_morm<-gls(CB_M~RoB_M, data=morm_avg, correlation=corPagel(0.9278729,phy=morm.tree,form=~Species, fixed=TRUE))
ROBCB_allgout<-gls(CB_M~RoB_M, data=allgout_avg, correlation=corPagel(1,phy=allgout.tree,form=~Species, fixed=TRUE))
ROBCB_allmout<-gls(CB_M~RoB_M, data=allmout_avg, correlation=corPagel(1,phy=allmout.tree,form=~Species, fixed=TRUE))
ROBCB_gwave<-gls(CB_M~RoB_M, data=gwave_avg, correlation=corPagel(0.2610847,phy=gwave.tree,form=~Species, fixed=TRUE))
ROBCB_gpulse<-gls(CB_M~RoB_M, data=gpulse_avg, correlation=corPagel(0.2610847,phy=gpulse.tree,form=~Species, fixed=TRUE))

######################################################
#  plotting uCT results, region x region, spp means  #
######################################################

ob_blank<-ggplot(avg_log_data, aes(x=OB, y=OB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(-.85, .1), ylim = c(-.85, .1));
ob_tel<-ggplot(avg_log_data, aes(x=OB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBTEL_tubegen$coefficients[2], intercept = OBTEL_tubegen$coefficients[1],color="Black")+geom_abline(slope = OBTEL_ampegen$coefficients[2], intercept = OBTEL_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OBTEL_out$coefficients[2], intercept = OBTEL_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_ellhb<-ggplot(avg_log_data, aes(x=OB, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBHB_tubegen$coefficients[2], intercept = OBHB_tubegen$coefficients[1],color="Black")+geom_abline(slope = OBHB_ampegen$coefficients[2], intercept = OBHB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OBHB_out$coefficients[2], intercept = OBHB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_ot<-ggplot(avg_log_data, aes(x=OB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBOT_tubegen$coefficients[2], intercept = OBOT_tubegen$coefficients[1],color="Black")+geom_abline(slope = OBOT_ampegen$coefficients[2], intercept = OBOT_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OBOT_out$coefficients[2], intercept = OBOT_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_ts<-ggplot(avg_log_data, aes(x=OB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBTS_tubegen$coefficients[2], intercept = OBTS_tubegen$coefficients[1],color="Black")+geom_abline(slope = OBTS_ampegen$coefficients[2], intercept = OBTS_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OBTS_out$coefficients[2], intercept = OBTS_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_cb<-ggplot(avg_log_data, aes(x=OB, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBCB_tubegen$coefficients[2], intercept = OBCB_tubegen$coefficients[1],color="Black")+geom_abline(slope = OBCB_ampegen$coefficients[2], intercept = OBCB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OBCB_out$coefficients[2], intercept = OBCB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_rob<-ggplot(avg_log_data, aes(x=OB, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBROB_tubegen$coefficients[2], intercept = OBROB_tubegen$coefficients[1],color="Black")+geom_abline(slope = OBROB_ampegen$coefficients[2], intercept = OBROB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OBROB_out$coefficients[2], intercept = OBROB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));

tel_ob<-ggplot(avg_log_data, aes(x=TEL, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELOB_tubegen$coefficients[2], intercept = TELOB_tubegen$coefficients[1],color="Black")+geom_abline(slope = TELOB_ampegen$coefficients[2], intercept = TELOB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TELOB_out$coefficients[2], intercept = TELOB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
tel_blank<-ggplot(avg_log_data, aes(x=TEL, y=TEL)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_ellhb<-ggplot(avg_log_data, aes(x=TEL, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELHB_tubegen$coefficients[2], intercept = TELHB_tubegen$coefficients[1],color="Black")+geom_abline(slope = TELHB_ampegen$coefficients[2], intercept = TELHB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TELHB_out$coefficients[2], intercept = TELHB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_ot<-ggplot(avg_log_data, aes(x=TEL, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELOT_tubegen$coefficients[2], intercept = TELOT_tubegen$coefficients[1],color="Black")+geom_abline(slope = TELOT_ampegen$coefficients[2], intercept = TELOT_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TELOT_out$coefficients[2], intercept = TELOT_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_ts<-ggplot(avg_log_data, aes(x=TEL, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELTS_tubegen$coefficients[2], intercept = TELTS_tubegen$coefficients[1],color="Black")+geom_abline(slope = TELTS_ampegen$coefficients[2], intercept = TELTS_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TELTS_out$coefficients[2], intercept = TELTS_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_cb<-ggplot(avg_log_data, aes(x=TEL, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELCB_tubegen$coefficients[2], intercept = TELCB_tubegen$coefficients[1],color="Black")+geom_abline(slope = TELCB_ampegen$coefficients[2], intercept = TELCB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TELCB_out$coefficients[2], intercept = TELCB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_rob<-ggplot(avg_log_data, aes(x=TEL, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELROB_tubegen$coefficients[2], intercept = TELROB_tubegen$coefficients[1],color="Black")+geom_abline(slope = TELROB_ampegen$coefficients[2], intercept = TELROB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TELROB_out$coefficients[2], intercept = TELROB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

ellhb_ob<-ggplot(avg_log_data, aes(x=HB, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBOB_tubegen$coefficients[2], intercept = HBOB_tubegen$coefficients[1],color="Black")+geom_abline(slope = HBOB_ampegen$coefficients[2], intercept = HBOB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = HBOB_out$coefficients[2], intercept = HBOB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
ellhb_tel<-ggplot(avg_log_data, aes(x=HB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBTEL_tubegen$coefficients[2], intercept = HBTEL_tubegen$coefficients[1],color="Black")+geom_abline(slope = HBTEL_ampegen$coefficients[2], intercept = HBTEL_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = HBTEL_out$coefficients[2], intercept = HBTEL_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_blank<-ggplot(avg_log_data, aes(x=HB, y=HB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_ot<-ggplot(avg_log_data, aes(x=HB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBOT_tubegen$coefficients[2], intercept = HBOT_tubegen$coefficients[1],color="Black")+geom_abline(slope = HBOT_ampegen$coefficients[2], intercept = HBOT_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = HBOT_out$coefficients[2], intercept = HBOT_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_ts<-ggplot(avg_log_data, aes(x=HB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBTS_tubegen$coefficients[2], intercept = HBTS_tubegen$coefficients[1],color="Black")+geom_abline(slope = HBTS_ampegen$coefficients[2], intercept = HBTS_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = HBTS_out$coefficients[2], intercept = HBTS_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_cb<-ggplot(avg_log_data, aes(x=HB, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBCB_tubegen$coefficients[2], intercept = HBCB_tubegen$coefficients[1],color="Black")+geom_abline(slope = HBCB_ampegen$coefficients[2], intercept = HBCB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = HBCB_out$coefficients[2], intercept = HBCB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_rob<-ggplot(avg_log_data, aes(x=HB, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBROB_tubegen$coefficients[2], intercept = HBROB_tubegen$coefficients[1],color="Black")+geom_abline(slope = HBROB_ampegen$coefficients[2], intercept = HBROB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = HBROB_out$coefficients[2], intercept = HBROB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

ot_ob<-ggplot(avg_log_data, aes(x=OT, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTOB_tubegen$coefficients[2], intercept = OTOB_tubegen$coefficients[1],color="Black")+geom_abline(slope = OTOB_ampegen$coefficients[2], intercept = OTOB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OTOB_out$coefficients[2], intercept = OTOB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
ot_tel<-ggplot(avg_log_data, aes(x=OT, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTTEL_tubegen$coefficients[2], intercept = OTTEL_tubegen$coefficients[1],color="Black")+geom_abline(slope = OTTEL_ampegen$coefficients[2], intercept = OTTEL_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OTTEL_out$coefficients[2], intercept = OTTEL_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_ellhb<-ggplot(avg_log_data, aes(x=OT, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTHB_tubegen$coefficients[2], intercept = OTHB_tubegen$coefficients[1],color="Black")+geom_abline(slope = OTHB_ampegen$coefficients[2], intercept = OTHB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OTHB_out$coefficients[2], intercept = OTHB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_blank<-ggplot(avg_log_data, aes(x=OT, y=OT)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_ts<-ggplot(avg_log_data, aes(x=OT, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTTS_tubegen$coefficients[2], intercept = OTTS_tubegen$coefficients[1],color="Black")+geom_abline(slope = OTTS_ampegen$coefficients[2], intercept = OTTS_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OTTS_out$coefficients[2], intercept = OTTS_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_cb<-ggplot(avg_log_data, aes(x=OT, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTCB_tubegen$coefficients[2], intercept = OTCB_tubegen$coefficients[1],color="Black")+geom_abline(slope = OTCB_ampegen$coefficients[2], intercept = OTCB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OTCB_out$coefficients[2], intercept = OTCB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_rob<-ggplot(avg_log_data, aes(x=OT, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTROB_tubegen$coefficients[2], intercept = OTROB_tubegen$coefficients[1],color="Black")+geom_abline(slope = OTROB_ampegen$coefficients[2], intercept = OTROB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = OTROB_out$coefficients[2], intercept = OTROB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

ts_ob<-ggplot(avg_log_data, aes(x=TS, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSOB_tubegen$coefficients[2], intercept = TSOB_tubegen$coefficients[1],color="Black")+geom_abline(slope = TSOB_ampegen$coefficients[2], intercept = TSOB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TSOB_out$coefficients[2], intercept = TSOB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
ts_tel<-ggplot(avg_log_data, aes(x=TS, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSTEL_tubegen$coefficients[2], intercept = TSTEL_tubegen$coefficients[1],color="Black")+geom_abline(slope = TSTEL_ampegen$coefficients[2], intercept = TSTEL_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TSTEL_out$coefficients[2], intercept = TSTEL_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_ellhb<-ggplot(avg_log_data, aes(x=TS, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSHB_tubegen$coefficients[2], intercept = TSHB_tubegen$coefficients[1],color="Black")+geom_abline(slope = TSHB_ampegen$coefficients[2], intercept = TSHB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TSHB_out$coefficients[2], intercept = TSHB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_ot<-ggplot(avg_log_data, aes(x=TS, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSOT_tubegen$coefficients[2], intercept = TSOT_tubegen$coefficients[1],color="Black")+geom_abline(slope = TSOT_ampegen$coefficients[2], intercept = TSOT_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TSOT_out$coefficients[2], intercept = TSOT_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_blank<-ggplot(avg_log_data, aes(x=TS, y=TS)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_cb<-ggplot(avg_log_data, aes(x=TS, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSCB_tubegen$coefficients[2], intercept = TSCB_tubegen$coefficients[1],color="Black")+geom_abline(slope = TSCB_ampegen$coefficients[2], intercept = TSCB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TSCB_out$coefficients[2], intercept = TSCB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_rob<-ggplot(avg_log_data, aes(x=TS, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSROB_tubegen$coefficients[2], intercept = TSROB_tubegen$coefficients[1],color="Black")+geom_abline(slope = TSROB_ampegen$coefficients[2], intercept = TSROB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = TSROB_out$coefficients[2], intercept = TSROB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

cb_ob<-ggplot(avg_log_data, aes(x=CB, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBOB_tubegen$coefficients[2], intercept = CBOB_tubegen$coefficients[1],color="Black")+geom_abline(slope = CBOB_ampegen$coefficients[2], intercept = CBOB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = CBOB_out$coefficients[2], intercept = CBOB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
cb_tel<-ggplot(avg_log_data, aes(x=CB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBTEL_tubegen$coefficients[2], intercept = CBTEL_tubegen$coefficients[1],color="Black")+geom_abline(slope = CBTEL_ampegen$coefficients[2], intercept = CBTEL_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = CBTEL_out$coefficients[2], intercept = CBTEL_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_ellhb<-ggplot(avg_log_data, aes(x=CB, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBHB_tubegen$coefficients[2], intercept = CBHB_tubegen$coefficients[1],color="Black")+geom_abline(slope = CBHB_ampegen$coefficients[2], intercept = CBHB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = CBHB_out$coefficients[2], intercept = CBHB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_ot<-ggplot(avg_log_data, aes(x=CB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBOT_tubegen$coefficients[2], intercept = CBOT_tubegen$coefficients[1],color="Black")+geom_abline(slope = CBOT_ampegen$coefficients[2], intercept = CBOT_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = CBOT_out$coefficients[2], intercept = CBOT_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_ts<-ggplot(avg_log_data, aes(x=CB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBTS_tubegen$coefficients[2], intercept = CBTS_tubegen$coefficients[1],color="Black")+geom_abline(slope = CBTS_ampegen$coefficients[2], intercept = CBTS_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = CBTS_out$coefficients[2], intercept = CBTS_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_blank<-ggplot(avg_log_data, aes(x=CB, y=CB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_rob<-ggplot(avg_log_data, aes(x=CB, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBROB_tubegen$coefficients[2], intercept = CBROB_tubegen$coefficients[1],color="Black")+geom_abline(slope = CBROB_ampegen$coefficients[2], intercept = CBROB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = CBROB_out$coefficients[2], intercept = CBROB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

rob_ob<-ggplot(avg_log_data, aes(x=ROB, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBOB_tubegen$coefficients[2], intercept = ROBOB_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROBOB_ampegen$coefficients[2], intercept = ROBOB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROBOB_out$coefficients[2], intercept = ROBOB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
rob_tel<-ggplot(avg_log_data, aes(x=ROB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBTEL_tubegen$coefficients[2], intercept = ROBTEL_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROBTEL_ampegen$coefficients[2], intercept = ROBTEL_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROBTEL_out$coefficients[2], intercept = ROBTEL_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_ellhb<-ggplot(avg_log_data, aes(x=ROB, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBHB_tubegen$coefficients[2], intercept = ROBHB_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROBHB_ampegen$coefficients[2], intercept = ROBHB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROBHB_out$coefficients[2], intercept = ROBHB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_ot<-ggplot(avg_log_data, aes(x=ROB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBOT_tubegen$coefficients[2], intercept = ROBOT_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROBOT_ampegen$coefficients[2], intercept = ROBOT_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROBOT_out$coefficients[2], intercept = ROBOT_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_ts<-ggplot(avg_log_data, aes(x=ROB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBTS_tubegen$coefficients[2], intercept = ROBTS_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROBTS_ampegen$coefficients[2], intercept = ROBTS_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROBTS_out$coefficients[2], intercept = ROBTS_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_cb<-ggplot(avg_log_data, aes(x=ROB, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBCB_tubegen$coefficients[2], intercept = ROBCB_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROBCB_ampegen$coefficients[2], intercept = ROBCB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROBCB_out$coefficients[2], intercept = ROBCB_out$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_blank<-ggplot(avg_log_data, aes(x=ROB, y=ROB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));

## ob, tel, hb, ot, ts, cb, rob
all_regionxregion_plot<-ggarrange(ob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()), ncol=7, nrow=7, common.legend = TRUE, legend = "none");
annotate_figure(all_regionxregion_plot, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

## ob, ot, tel, rob, hb, ts, cb
all_regionxregion_plot2<-ggarrange(ob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()), ncol=7, nrow=7, common.legend = TRUE, legend = "none");
annotate_figure(all_regionxregion_plot2, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

## ob, ot, rob, tel, hb, ts, cb
all_regionxregion_plot3<-ggarrange(ob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()), ncol=7, nrow=7, common.legend = TRUE, legend = "none");
annotate_figure(all_regionxregion_plot3, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

## saving plots to pdf
# Customizing the output
pdf("regionxregion_esens_means_Feb2.pdf",         # File name
    width = 7.874, height = 7.874) # Width and height in inches

#all_regionxregion_plot3
annotate_figure(all_regionxregion_plot3, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

# Closing the graphical device
dev.off() 


##################################################################
#  plotting uCT results, region x region, spp means, morm v gym  #
##################################################################

mormgym_avg<-subset(avg_log_data, Species %in% tubegen_names)

ob_blank<-ggplot(mormgym_avg, aes(x=OB, y=OB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(-.85, .1), ylim = c(-.85, .1));
ob_tel<-ggplot(mormgym_avg, aes(x=OB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBTEL_gym$coefficients[2], intercept = OBTEL_gym$coefficients[1],color="Black")+geom_abline(slope = OBTEL_morm$coefficients[2], intercept = OBTEL_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_ellhb<-ggplot(mormgym_avg, aes(x=OB, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBHB_gym$coefficients[2], intercept = OBHB_gym$coefficients[1],color="Black")+geom_abline(slope = OBHB_morm$coefficients[2], intercept = OBHB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_ot<-ggplot(mormgym_avg, aes(x=OB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBOT_gym$coefficients[2], intercept = OBOT_gym$coefficients[1],color="Black")+geom_abline(slope = OBOT_morm$coefficients[2], intercept = OBOT_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_ts<-ggplot(mormgym_avg, aes(x=OB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBTS_gym$coefficients[2], intercept = OBTS_gym$coefficients[1],color="Black")+geom_abline(slope = OBTS_morm$coefficients[2], intercept = OBTS_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_cb<-ggplot(mormgym_avg, aes(x=OB, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBCB_gym$coefficients[2], intercept = OBCB_gym$coefficients[1],color="Black")+geom_abline(slope = OBCB_morm$coefficients[2], intercept = OBCB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_rob<-ggplot(mormgym_avg, aes(x=OB, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBROB_gym$coefficients[2], intercept = OBROB_gym$coefficients[1],color="Black")+geom_abline(slope = OBROB_morm$coefficients[2], intercept = OBROB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));

tel_ob<-ggplot(mormgym_avg, aes(x=TEL, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELOB_gym$coefficients[2], intercept = TELOB_gym$coefficients[1],color="Black")+geom_abline(slope = TELOB_morm$coefficients[2], intercept = TELOB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
tel_blank<-ggplot(mormgym_avg, aes(x=TEL, y=TEL)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_ellhb<-ggplot(mormgym_avg, aes(x=TEL, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELHB_gym$coefficients[2], intercept = TELHB_gym$coefficients[1],color="Black")+geom_abline(slope = TELHB_morm$coefficients[2], intercept = TELHB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_ot<-ggplot(mormgym_avg, aes(x=TEL, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELOT_gym$coefficients[2], intercept = TELOT_gym$coefficients[1],color="Black")+geom_abline(slope = TELOT_morm$coefficients[2], intercept = TELOT_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_ts<-ggplot(mormgym_avg, aes(x=TEL, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELTS_gym$coefficients[2], intercept = TELTS_gym$coefficients[1],color="Black")+geom_abline(slope = TELTS_morm$coefficients[2], intercept = TELTS_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_cb<-ggplot(mormgym_avg, aes(x=TEL, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELCB_gym$coefficients[2], intercept = TELCB_gym$coefficients[1],color="Black")+geom_abline(slope = TELCB_morm$coefficients[2], intercept = TELCB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_rob<-ggplot(mormgym_avg, aes(x=TEL, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELROB_gym$coefficients[2], intercept = TELROB_gym$coefficients[1],color="Black")+geom_abline(slope = TELROB_morm$coefficients[2], intercept = TELROB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

ellhb_ob<-ggplot(mormgym_avg, aes(x=HB, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBOB_gym$coefficients[2], intercept = HBOB_gym$coefficients[1],color="Black")+geom_abline(slope = HBOB_morm$coefficients[2], intercept = HBOB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
ellhb_tel<-ggplot(mormgym_avg, aes(x=HB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBTEL_gym$coefficients[2], intercept = HBTEL_gym$coefficients[1],color="Black")+geom_abline(slope = HBTEL_morm$coefficients[2], intercept = HBTEL_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_blank<-ggplot(mormgym_avg, aes(x=HB, y=HB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_ot<-ggplot(mormgym_avg, aes(x=HB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBOT_gym$coefficients[2], intercept = HBOT_gym$coefficients[1],color="Black")+geom_abline(slope = HBOT_morm$coefficients[2], intercept = HBOT_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_ts<-ggplot(mormgym_avg, aes(x=HB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBTS_gym$coefficients[2], intercept = HBTS_gym$coefficients[1],color="Black")+geom_abline(slope = HBTS_morm$coefficients[2], intercept = HBTS_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_cb<-ggplot(mormgym_avg, aes(x=HB, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBCB_gym$coefficients[2], intercept = HBCB_gym$coefficients[1],color="Black")+geom_abline(slope = HBCB_morm$coefficients[2], intercept = HBCB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_rob<-ggplot(mormgym_avg, aes(x=HB, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBROB_gym$coefficients[2], intercept = HBROB_gym$coefficients[1],color="Black")+geom_abline(slope = HBROB_morm$coefficients[2], intercept = HBROB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

ot_ob<-ggplot(mormgym_avg, aes(x=OT, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTOB_gym$coefficients[2], intercept = OTOB_gym$coefficients[1],color="Black")+geom_abline(slope = OTOB_morm$coefficients[2], intercept = OTOB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
ot_tel<-ggplot(mormgym_avg, aes(x=OT, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTTEL_gym$coefficients[2], intercept = OTTEL_gym$coefficients[1],color="Black")+geom_abline(slope = OTTEL_morm$coefficients[2], intercept = OTTEL_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_ellhb<-ggplot(mormgym_avg, aes(x=OT, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTHB_gym$coefficients[2], intercept = OTHB_gym$coefficients[1],color="Black")+geom_abline(slope = OTHB_morm$coefficients[2], intercept = OTHB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_blank<-ggplot(mormgym_avg, aes(x=OT, y=OT)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_ts<-ggplot(mormgym_avg, aes(x=OT, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTTS_gym$coefficients[2], intercept = OTTS_gym$coefficients[1],color="Black")+geom_abline(slope = OTTS_morm$coefficients[2], intercept = OTTS_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_cb<-ggplot(mormgym_avg, aes(x=OT, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTCB_gym$coefficients[2], intercept = OTCB_gym$coefficients[1],color="Black")+geom_abline(slope = OTCB_morm$coefficients[2], intercept = OTCB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_rob<-ggplot(mormgym_avg, aes(x=OT, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTROB_gym$coefficients[2], intercept = OTROB_gym$coefficients[1],color="Black")+geom_abline(slope = OTROB_morm$coefficients[2], intercept = OTROB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

ts_ob<-ggplot(mormgym_avg, aes(x=TS, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSOB_gym$coefficients[2], intercept = TSOB_gym$coefficients[1],color="Black")+geom_abline(slope = TSOB_morm$coefficients[2], intercept = TSOB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
ts_tel<-ggplot(mormgym_avg, aes(x=TS, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSTEL_gym$coefficients[2], intercept = TSTEL_gym$coefficients[1],color="Black")+geom_abline(slope = TSTEL_morm$coefficients[2], intercept = TSTEL_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_ellhb<-ggplot(mormgym_avg, aes(x=TS, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSHB_gym$coefficients[2], intercept = TSHB_gym$coefficients[1],color="Black")+geom_abline(slope = TSHB_morm$coefficients[2], intercept = TSHB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_ot<-ggplot(mormgym_avg, aes(x=TS, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSOT_gym$coefficients[2], intercept = TSOT_gym$coefficients[1],color="Black")+geom_abline(slope = TSOT_morm$coefficients[2], intercept = TSOT_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_blank<-ggplot(mormgym_avg, aes(x=TS, y=TS)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_cb<-ggplot(mormgym_avg, aes(x=TS, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSCB_gym$coefficients[2], intercept = TSCB_gym$coefficients[1],color="Black")+geom_abline(slope = TSCB_morm$coefficients[2], intercept = TSCB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_rob<-ggplot(mormgym_avg, aes(x=TS, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSROB_gym$coefficients[2], intercept = TSROB_gym$coefficients[1],color="Black")+geom_abline(slope = TSROB_morm$coefficients[2], intercept = TSROB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

cb_ob<-ggplot(mormgym_avg, aes(x=CB, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBOB_gym$coefficients[2], intercept = CBOB_gym$coefficients[1],color="Black")+geom_abline(slope = CBOB_morm$coefficients[2], intercept = CBOB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
cb_tel<-ggplot(mormgym_avg, aes(x=CB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBTEL_gym$coefficients[2], intercept = CBTEL_gym$coefficients[1],color="Black")+geom_abline(slope = CBTEL_morm$coefficients[2], intercept = CBTEL_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_ellhb<-ggplot(mormgym_avg, aes(x=CB, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBHB_gym$coefficients[2], intercept = CBHB_gym$coefficients[1],color="Black")+geom_abline(slope = CBHB_morm$coefficients[2], intercept = CBHB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_ot<-ggplot(mormgym_avg, aes(x=CB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBOT_gym$coefficients[2], intercept = CBOT_gym$coefficients[1],color="Black")+geom_abline(slope = CBOT_morm$coefficients[2], intercept = CBOT_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_ts<-ggplot(mormgym_avg, aes(x=CB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBTS_gym$coefficients[2], intercept = CBTS_gym$coefficients[1],color="Black")+geom_abline(slope = CBTS_morm$coefficients[2], intercept = CBTS_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_blank<-ggplot(mormgym_avg, aes(x=CB, y=CB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_rob<-ggplot(mormgym_avg, aes(x=CB, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBROB_gym$coefficients[2], intercept = CBROB_gym$coefficients[1],color="Black")+geom_abline(slope = CBROB_morm$coefficients[2], intercept = CBROB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

rob_ob<-ggplot(mormgym_avg, aes(x=ROB, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBOB_gym$coefficients[2], intercept = ROBOB_gym$coefficients[1],color="Black")+geom_abline(slope = ROBOB_morm$coefficients[2], intercept = ROBOB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
rob_tel<-ggplot(mormgym_avg, aes(x=ROB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBTEL_gym$coefficients[2], intercept = ROBTEL_gym$coefficients[1],color="Black")+geom_abline(slope = ROBTEL_morm$coefficients[2], intercept = ROBTEL_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_ellhb<-ggplot(mormgym_avg, aes(x=ROB, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBHB_gym$coefficients[2], intercept = ROBHB_gym$coefficients[1],color="Black")+geom_abline(slope = ROBHB_morm$coefficients[2], intercept = ROBHB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_ot<-ggplot(mormgym_avg, aes(x=ROB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBOT_gym$coefficients[2], intercept = ROBOT_gym$coefficients[1],color="Black")+geom_abline(slope = ROBOT_morm$coefficients[2], intercept = ROBOT_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_ts<-ggplot(mormgym_avg, aes(x=ROB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBTS_gym$coefficients[2], intercept = ROBTS_gym$coefficients[1],color="Black")+geom_abline(slope = ROBTS_morm$coefficients[2], intercept = ROBTS_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_cb<-ggplot(mormgym_avg, aes(x=ROB, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBCB_gym$coefficients[2], intercept = ROBCB_gym$coefficients[1],color="Black")+geom_abline(slope = ROBCB_morm$coefficients[2], intercept = ROBCB_morm$coefficients[1],color="Black",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_blank<-ggplot(mormgym_avg, aes(x=ROB, y=ROB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));

## ob, tel, hb, ot, ts, cb, rob
all_regionxregion_plot_mormgym<-ggarrange(ob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()), ncol=7, nrow=7, common.legend = TRUE, legend = "none");
annotate_figure(all_regionxregion_plot_mormgym, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

## ob, ot, tel, rob, hb, ts, cb
all_regionxregion_plot_mormgym2<-ggarrange(ob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()), ncol=7, nrow=7, common.legend = TRUE, legend = "none");
annotate_figure(all_regionxregion_plot_mormgym2, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

## ob, ot, rob, tel, hb, ts, cb
all_regionxregion_plot_mormgym3<-ggarrange(ob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()), ncol=7, nrow=7, common.legend = TRUE, legend = "none");
annotate_figure(all_regionxregion_plot_mormgym3, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

## saving plots to pdf
# Customizing the mormput
pdf("regionxregion_mormgym_means_Feb3.pdf",         # File name
    width = 7.874, height = 7.874) # Width and height in inches

annotate_figure(all_regionxregion_plot_mormgym3, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

# Closing the graphical device
dev.off() 


###################################################################
#  plotting uCT results, region x region, spp means, mout v gout  #
###################################################################

allmoutallgout_avg<-subset(avg_log_data, Species %in% out_names)

ob_blank<-ggplot(allmoutallgout_avg, aes(x=OB, y=OB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(-.85, .1), ylim = c(-.85, .1));
ob_tel<-ggplot(allmoutallgout_avg, aes(x=OB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBTEL_allgout$coefficients[2], intercept = OBTEL_allgout$coefficients[1],color="grey")+geom_abline(slope = OBTEL_allmout$coefficients[2], intercept = OBTEL_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_ellhb<-ggplot(allmoutallgout_avg, aes(x=OB, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBHB_allgout$coefficients[2], intercept = OBHB_allgout$coefficients[1],color="grey")+geom_abline(slope = OBHB_allmout$coefficients[2], intercept = OBHB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_ot<-ggplot(allmoutallgout_avg, aes(x=OB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBOT_allgout$coefficients[2], intercept = OBOT_allgout$coefficients[1],color="grey")+geom_abline(slope = OBOT_allmout$coefficients[2], intercept = OBOT_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_ts<-ggplot(allmoutallgout_avg, aes(x=OB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBTS_allgout$coefficients[2], intercept = OBTS_allgout$coefficients[1],color="grey")+geom_abline(slope = OBTS_allmout$coefficients[2], intercept = OBTS_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_cb<-ggplot(allmoutallgout_avg, aes(x=OB, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBCB_allgout$coefficients[2], intercept = OBCB_allgout$coefficients[1],color="grey")+geom_abline(slope = OBCB_allmout$coefficients[2], intercept = OBCB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_rob<-ggplot(allmoutallgout_avg, aes(x=OB, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBROB_allgout$coefficients[2], intercept = OBROB_allgout$coefficients[1],color="grey")+geom_abline(slope = OBROB_allmout$coefficients[2], intercept = OBROB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));

tel_ob<-ggplot(allmoutallgout_avg, aes(x=TEL, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELOB_allgout$coefficients[2], intercept = TELOB_allgout$coefficients[1],color="grey")+geom_abline(slope = TELOB_allmout$coefficients[2], intercept = TELOB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
tel_blank<-ggplot(allmoutallgout_avg, aes(x=TEL, y=TEL)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_ellhb<-ggplot(allmoutallgout_avg, aes(x=TEL, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELHB_allgout$coefficients[2], intercept = TELHB_allgout$coefficients[1],color="grey")+geom_abline(slope = TELHB_allmout$coefficients[2], intercept = TELHB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_ot<-ggplot(allmoutallgout_avg, aes(x=TEL, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELOT_allgout$coefficients[2], intercept = TELOT_allgout$coefficients[1],color="grey")+geom_abline(slope = TELOT_allmout$coefficients[2], intercept = TELOT_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_ts<-ggplot(allmoutallgout_avg, aes(x=TEL, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELTS_allgout$coefficients[2], intercept = TELTS_allgout$coefficients[1],color="grey")+geom_abline(slope = TELTS_allmout$coefficients[2], intercept = TELTS_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_cb<-ggplot(allmoutallgout_avg, aes(x=TEL, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELCB_allgout$coefficients[2], intercept = TELCB_allgout$coefficients[1],color="grey")+geom_abline(slope = TELCB_allmout$coefficients[2], intercept = TELCB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_rob<-ggplot(allmoutallgout_avg, aes(x=TEL, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELROB_allgout$coefficients[2], intercept = TELROB_allgout$coefficients[1],color="grey")+geom_abline(slope = TELROB_allmout$coefficients[2], intercept = TELROB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

ellhb_ob<-ggplot(allmoutallgout_avg, aes(x=HB, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBOB_allgout$coefficients[2], intercept = HBOB_allgout$coefficients[1],color="grey")+geom_abline(slope = HBOB_allmout$coefficients[2], intercept = HBOB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
ellhb_tel<-ggplot(allmoutallgout_avg, aes(x=HB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBTEL_allgout$coefficients[2], intercept = HBTEL_allgout$coefficients[1],color="grey")+geom_abline(slope = HBTEL_allmout$coefficients[2], intercept = HBTEL_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_blank<-ggplot(allmoutallgout_avg, aes(x=HB, y=HB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_ot<-ggplot(allmoutallgout_avg, aes(x=HB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBOT_allgout$coefficients[2], intercept = HBOT_allgout$coefficients[1],color="grey")+geom_abline(slope = HBOT_allmout$coefficients[2], intercept = HBOT_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_ts<-ggplot(allmoutallgout_avg, aes(x=HB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBTS_allgout$coefficients[2], intercept = HBTS_allgout$coefficients[1],color="grey")+geom_abline(slope = HBTS_allmout$coefficients[2], intercept = HBTS_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_cb<-ggplot(allmoutallgout_avg, aes(x=HB, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBCB_allgout$coefficients[2], intercept = HBCB_allgout$coefficients[1],color="grey")+geom_abline(slope = HBCB_allmout$coefficients[2], intercept = HBCB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_rob<-ggplot(allmoutallgout_avg, aes(x=HB, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBROB_allgout$coefficients[2], intercept = HBROB_allgout$coefficients[1],color="grey")+geom_abline(slope = HBROB_allmout$coefficients[2], intercept = HBROB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

ot_ob<-ggplot(allmoutallgout_avg, aes(x=OT, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTOB_allgout$coefficients[2], intercept = OTOB_allgout$coefficients[1],color="grey")+geom_abline(slope = OTOB_allmout$coefficients[2], intercept = OTOB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
ot_tel<-ggplot(allmoutallgout_avg, aes(x=OT, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTTEL_allgout$coefficients[2], intercept = OTTEL_allgout$coefficients[1],color="grey")+geom_abline(slope = OTTEL_allmout$coefficients[2], intercept = OTTEL_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_ellhb<-ggplot(allmoutallgout_avg, aes(x=OT, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTHB_allgout$coefficients[2], intercept = OTHB_allgout$coefficients[1],color="grey")+geom_abline(slope = OTHB_allmout$coefficients[2], intercept = OTHB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_blank<-ggplot(allmoutallgout_avg, aes(x=OT, y=OT)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_ts<-ggplot(allmoutallgout_avg, aes(x=OT, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTTS_allgout$coefficients[2], intercept = OTTS_allgout$coefficients[1],color="grey")+geom_abline(slope = OTTS_allmout$coefficients[2], intercept = OTTS_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_cb<-ggplot(allmoutallgout_avg, aes(x=OT, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTCB_allgout$coefficients[2], intercept = OTCB_allgout$coefficients[1],color="grey")+geom_abline(slope = OTCB_allmout$coefficients[2], intercept = OTCB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_rob<-ggplot(allmoutallgout_avg, aes(x=OT, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTROB_allgout$coefficients[2], intercept = OTROB_allgout$coefficients[1],color="grey")+geom_abline(slope = OTROB_allmout$coefficients[2], intercept = OTROB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

ts_ob<-ggplot(allmoutallgout_avg, aes(x=TS, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSOB_allgout$coefficients[2], intercept = TSOB_allgout$coefficients[1],color="grey")+geom_abline(slope = TSOB_allmout$coefficients[2], intercept = TSOB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
ts_tel<-ggplot(allmoutallgout_avg, aes(x=TS, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSTEL_allgout$coefficients[2], intercept = TSTEL_allgout$coefficients[1],color="grey")+geom_abline(slope = TSTEL_allmout$coefficients[2], intercept = TSTEL_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_ellhb<-ggplot(allmoutallgout_avg, aes(x=TS, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSHB_allgout$coefficients[2], intercept = TSHB_allgout$coefficients[1],color="grey")+geom_abline(slope = TSHB_allmout$coefficients[2], intercept = TSHB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_ot<-ggplot(allmoutallgout_avg, aes(x=TS, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSOT_allgout$coefficients[2], intercept = TSOT_allgout$coefficients[1],color="grey")+geom_abline(slope = TSOT_allmout$coefficients[2], intercept = TSOT_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_blank<-ggplot(allmoutallgout_avg, aes(x=TS, y=TS)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_cb<-ggplot(allmoutallgout_avg, aes(x=TS, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSCB_allgout$coefficients[2], intercept = TSCB_allgout$coefficients[1],color="grey")+geom_abline(slope = TSCB_allmout$coefficients[2], intercept = TSCB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_rob<-ggplot(allmoutallgout_avg, aes(x=TS, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSROB_allgout$coefficients[2], intercept = TSROB_allgout$coefficients[1],color="grey")+geom_abline(slope = TSROB_allmout$coefficients[2], intercept = TSROB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

cb_ob<-ggplot(allmoutallgout_avg, aes(x=CB, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBOB_allgout$coefficients[2], intercept = CBOB_allgout$coefficients[1],color="grey")+geom_abline(slope = CBOB_allmout$coefficients[2], intercept = CBOB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
cb_tel<-ggplot(allmoutallgout_avg, aes(x=CB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBTEL_allgout$coefficients[2], intercept = CBTEL_allgout$coefficients[1],color="grey")+geom_abline(slope = CBTEL_allmout$coefficients[2], intercept = CBTEL_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_ellhb<-ggplot(allmoutallgout_avg, aes(x=CB, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBHB_allgout$coefficients[2], intercept = CBHB_allgout$coefficients[1],color="grey")+geom_abline(slope = CBHB_allmout$coefficients[2], intercept = CBHB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_ot<-ggplot(allmoutallgout_avg, aes(x=CB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBOT_allgout$coefficients[2], intercept = CBOT_allgout$coefficients[1],color="grey")+geom_abline(slope = CBOT_allmout$coefficients[2], intercept = CBOT_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_ts<-ggplot(allmoutallgout_avg, aes(x=CB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBTS_allgout$coefficients[2], intercept = CBTS_allgout$coefficients[1],color="grey")+geom_abline(slope = CBTS_allmout$coefficients[2], intercept = CBTS_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_blank<-ggplot(allmoutallgout_avg, aes(x=CB, y=CB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_rob<-ggplot(allmoutallgout_avg, aes(x=CB, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBROB_allgout$coefficients[2], intercept = CBROB_allgout$coefficients[1],color="grey")+geom_abline(slope = CBROB_allmout$coefficients[2], intercept = CBROB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

rob_ob<-ggplot(allmoutallgout_avg, aes(x=ROB, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBOB_allgout$coefficients[2], intercept = ROBOB_allgout$coefficients[1],color="grey")+geom_abline(slope = ROBOB_allmout$coefficients[2], intercept = ROBOB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
rob_tel<-ggplot(allmoutallgout_avg, aes(x=ROB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBTEL_allgout$coefficients[2], intercept = ROBTEL_allgout$coefficients[1],color="grey")+geom_abline(slope = ROBTEL_allmout$coefficients[2], intercept = ROBTEL_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_ellhb<-ggplot(allmoutallgout_avg, aes(x=ROB, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBHB_allgout$coefficients[2], intercept = ROBHB_allgout$coefficients[1],color="grey")+geom_abline(slope = ROBHB_allmout$coefficients[2], intercept = ROBHB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_ot<-ggplot(allmoutallgout_avg, aes(x=ROB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBOT_allgout$coefficients[2], intercept = ROBOT_allgout$coefficients[1],color="grey")+geom_abline(slope = ROBOT_allmout$coefficients[2], intercept = ROBOT_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_ts<-ggplot(allmoutallgout_avg, aes(x=ROB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBTS_allgout$coefficients[2], intercept = ROBTS_allgout$coefficients[1],color="grey")+geom_abline(slope = ROBTS_allmout$coefficients[2], intercept = ROBTS_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_cb<-ggplot(allmoutallgout_avg, aes(x=ROB, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBCB_allgout$coefficients[2], intercept = ROBCB_allgout$coefficients[1],color="grey")+geom_abline(slope = ROBCB_allmout$coefficients[2], intercept = ROBCB_allmout$coefficients[1],color="grey",lty=2)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_blank<-ggplot(allmoutallgout_avg, aes(x=ROB, y=ROB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));

## ob, tel, hb, ot, ts, cb, rob
all_regionxregion_plot_allmoutallgout<-ggarrange(ob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()), ncol=7, nrow=7, common.legend = TRUE, legend = "none");
annotate_figure(all_regionxregion_plot_allmoutallgout, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

## ob, ot, tel, rob, hb, ts, cb
all_regionxregion_plot_allmoutallgout2<-ggarrange(ob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()), ncol=7, nrow=7, common.legend = TRUE, legend = "none");
annotate_figure(all_regionxregion_plot_allmoutallgout2, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

## ob, ot, rob, tel, hb, ts, cb
all_regionxregion_plot_allmoutallgout3<-ggarrange(ob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()), ncol=7, nrow=7, common.legend = TRUE, legend = "none");
annotate_figure(all_regionxregion_plot_allmoutallgout3, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

## saving plots to pdf
# Customizing the allmoutput
pdf("regionxregion_allmoutallgout_means_Feb3.pdf",         # File name
    width = 7.874, height = 7.874) # Width and height in inches

annotate_figure(all_regionxregion_plot_allmoutallgout3, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

# Closing the graphical device
dev.off() 


######################################################################
#  plotting uCT results, region x region, spp means, gwave v gpulse  #
######################################################################

gwavegpulse_avg<-subset(avg_log_data, Species %in% gym_names)

ob_blank<-ggplot(gwavegpulse_avg, aes(x=OB, y=OB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(-.85, .1), ylim = c(-.85, .1));
ob_tel<-ggplot(gwavegpulse_avg, aes(x=OB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBTEL_gwave$coefficients[2], intercept = OBTEL_gwave$coefficients[1],color="blue")+geom_abline(slope = OBTEL_gpulse$coefficients[2], intercept = OBTEL_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_ellhb<-ggplot(gwavegpulse_avg, aes(x=OB, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBHB_gwave$coefficients[2], intercept = OBHB_gwave$coefficients[1],color="blue")+geom_abline(slope = OBHB_gpulse$coefficients[2], intercept = OBHB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_ot<-ggplot(gwavegpulse_avg, aes(x=OB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBOT_gwave$coefficients[2], intercept = OBOT_gwave$coefficients[1],color="blue")+geom_abline(slope = OBOT_gpulse$coefficients[2], intercept = OBOT_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_ts<-ggplot(gwavegpulse_avg, aes(x=OB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBTS_gwave$coefficients[2], intercept = OBTS_gwave$coefficients[1],color="blue")+geom_abline(slope = OBTS_gpulse$coefficients[2], intercept = OBTS_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_cb<-ggplot(gwavegpulse_avg, aes(x=OB, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBCB_gwave$coefficients[2], intercept = OBCB_gwave$coefficients[1],color="blue")+geom_abline(slope = OBCB_gpulse$coefficients[2], intercept = OBCB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));
ob_rob<-ggplot(gwavegpulse_avg, aes(x=OB, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OBROB_gwave$coefficients[2], intercept = OBROB_gwave$coefficients[1],color="blue")+geom_abline(slope = OBROB_gpulse$coefficients[2], intercept = OBROB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(-.85, .2), ylim = c(.4, 1.8));

tel_ob<-ggplot(gwavegpulse_avg, aes(x=TEL, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELOB_gwave$coefficients[2], intercept = TELOB_gwave$coefficients[1],color="blue")+geom_abline(slope = TELOB_gpulse$coefficients[2], intercept = TELOB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
tel_blank<-ggplot(gwavegpulse_avg, aes(x=TEL, y=TEL)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_ellhb<-ggplot(gwavegpulse_avg, aes(x=TEL, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELHB_gwave$coefficients[2], intercept = TELHB_gwave$coefficients[1],color="blue")+geom_abline(slope = TELHB_gpulse$coefficients[2], intercept = TELHB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_ot<-ggplot(gwavegpulse_avg, aes(x=TEL, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELOT_gwave$coefficients[2], intercept = TELOT_gwave$coefficients[1],color="blue")+geom_abline(slope = TELOT_gpulse$coefficients[2], intercept = TELOT_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_ts<-ggplot(gwavegpulse_avg, aes(x=TEL, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELTS_gwave$coefficients[2], intercept = TELTS_gwave$coefficients[1],color="blue")+geom_abline(slope = TELTS_gpulse$coefficients[2], intercept = TELTS_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_cb<-ggplot(gwavegpulse_avg, aes(x=TEL, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELCB_gwave$coefficients[2], intercept = TELCB_gwave$coefficients[1],color="blue")+geom_abline(slope = TELCB_gpulse$coefficients[2], intercept = TELCB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
tel_rob<-ggplot(gwavegpulse_avg, aes(x=TEL, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TELROB_gwave$coefficients[2], intercept = TELROB_gwave$coefficients[1],color="blue")+geom_abline(slope = TELROB_gpulse$coefficients[2], intercept = TELROB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

ellhb_ob<-ggplot(gwavegpulse_avg, aes(x=HB, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBOB_gwave$coefficients[2], intercept = HBOB_gwave$coefficients[1],color="blue")+geom_abline(slope = HBOB_gpulse$coefficients[2], intercept = HBOB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
ellhb_tel<-ggplot(gwavegpulse_avg, aes(x=HB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBTEL_gwave$coefficients[2], intercept = HBTEL_gwave$coefficients[1],color="blue")+geom_abline(slope = HBTEL_gpulse$coefficients[2], intercept = HBTEL_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_blank<-ggplot(gwavegpulse_avg, aes(x=HB, y=HB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_ot<-ggplot(gwavegpulse_avg, aes(x=HB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBOT_gwave$coefficients[2], intercept = HBOT_gwave$coefficients[1],color="blue")+geom_abline(slope = HBOT_gpulse$coefficients[2], intercept = HBOT_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_ts<-ggplot(gwavegpulse_avg, aes(x=HB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBTS_gwave$coefficients[2], intercept = HBTS_gwave$coefficients[1],color="blue")+geom_abline(slope = HBTS_gpulse$coefficients[2], intercept = HBTS_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_cb<-ggplot(gwavegpulse_avg, aes(x=HB, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBCB_gwave$coefficients[2], intercept = HBCB_gwave$coefficients[1],color="blue")+geom_abline(slope = HBCB_gpulse$coefficients[2], intercept = HBCB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ellhb_rob<-ggplot(gwavegpulse_avg, aes(x=HB, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = HBROB_gwave$coefficients[2], intercept = HBROB_gwave$coefficients[1],color="blue")+geom_abline(slope = HBROB_gpulse$coefficients[2], intercept = HBROB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

ot_ob<-ggplot(gwavegpulse_avg, aes(x=OT, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTOB_gwave$coefficients[2], intercept = OTOB_gwave$coefficients[1],color="blue")+geom_abline(slope = OTOB_gpulse$coefficients[2], intercept = OTOB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
ot_tel<-ggplot(gwavegpulse_avg, aes(x=OT, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTTEL_gwave$coefficients[2], intercept = OTTEL_gwave$coefficients[1],color="blue")+geom_abline(slope = OTTEL_gpulse$coefficients[2], intercept = OTTEL_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_ellhb<-ggplot(gwavegpulse_avg, aes(x=OT, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTHB_gwave$coefficients[2], intercept = OTHB_gwave$coefficients[1],color="blue")+geom_abline(slope = OTHB_gpulse$coefficients[2], intercept = OTHB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_blank<-ggplot(gwavegpulse_avg, aes(x=OT, y=OT)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_ts<-ggplot(gwavegpulse_avg, aes(x=OT, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTTS_gwave$coefficients[2], intercept = OTTS_gwave$coefficients[1],color="blue")+geom_abline(slope = OTTS_gpulse$coefficients[2], intercept = OTTS_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_cb<-ggplot(gwavegpulse_avg, aes(x=OT, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTCB_gwave$coefficients[2], intercept = OTCB_gwave$coefficients[1],color="blue")+geom_abline(slope = OTCB_gpulse$coefficients[2], intercept = OTCB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ot_rob<-ggplot(gwavegpulse_avg, aes(x=OT, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = OTROB_gwave$coefficients[2], intercept = OTROB_gwave$coefficients[1],color="blue")+geom_abline(slope = OTROB_gpulse$coefficients[2], intercept = OTROB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

ts_ob<-ggplot(gwavegpulse_avg, aes(x=TS, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSOB_gwave$coefficients[2], intercept = TSOB_gwave$coefficients[1],color="blue")+geom_abline(slope = TSOB_gpulse$coefficients[2], intercept = TSOB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
ts_tel<-ggplot(gwavegpulse_avg, aes(x=TS, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSTEL_gwave$coefficients[2], intercept = TSTEL_gwave$coefficients[1],color="blue")+geom_abline(slope = TSTEL_gpulse$coefficients[2], intercept = TSTEL_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_ellhb<-ggplot(gwavegpulse_avg, aes(x=TS, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSHB_gwave$coefficients[2], intercept = TSHB_gwave$coefficients[1],color="blue")+geom_abline(slope = TSHB_gpulse$coefficients[2], intercept = TSHB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_ot<-ggplot(gwavegpulse_avg, aes(x=TS, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSOT_gwave$coefficients[2], intercept = TSOT_gwave$coefficients[1],color="blue")+geom_abline(slope = TSOT_gpulse$coefficients[2], intercept = TSOT_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_blank<-ggplot(gwavegpulse_avg, aes(x=TS, y=TS)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_cb<-ggplot(gwavegpulse_avg, aes(x=TS, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSCB_gwave$coefficients[2], intercept = TSCB_gwave$coefficients[1],color="blue")+geom_abline(slope = TSCB_gpulse$coefficients[2], intercept = TSCB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
ts_rob<-ggplot(gwavegpulse_avg, aes(x=TS, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = TSROB_gwave$coefficients[2], intercept = TSROB_gwave$coefficients[1],color="blue")+geom_abline(slope = TSROB_gpulse$coefficients[2], intercept = TSROB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

cb_ob<-ggplot(gwavegpulse_avg, aes(x=CB, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBOB_gwave$coefficients[2], intercept = CBOB_gwave$coefficients[1],color="blue")+geom_abline(slope = CBOB_gpulse$coefficients[2], intercept = CBOB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
cb_tel<-ggplot(gwavegpulse_avg, aes(x=CB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBTEL_gwave$coefficients[2], intercept = CBTEL_gwave$coefficients[1],color="blue")+geom_abline(slope = CBTEL_gpulse$coefficients[2], intercept = CBTEL_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_ellhb<-ggplot(gwavegpulse_avg, aes(x=CB, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBHB_gwave$coefficients[2], intercept = CBHB_gwave$coefficients[1],color="blue")+geom_abline(slope = CBHB_gpulse$coefficients[2], intercept = CBHB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_ot<-ggplot(gwavegpulse_avg, aes(x=CB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBOT_gwave$coefficients[2], intercept = CBOT_gwave$coefficients[1],color="blue")+geom_abline(slope = CBOT_gpulse$coefficients[2], intercept = CBOT_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_ts<-ggplot(gwavegpulse_avg, aes(x=CB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBTS_gwave$coefficients[2], intercept = CBTS_gwave$coefficients[1],color="blue")+geom_abline(slope = CBTS_gpulse$coefficients[2], intercept = CBTS_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_blank<-ggplot(gwavegpulse_avg, aes(x=CB, y=CB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));
cb_rob<-ggplot(gwavegpulse_avg, aes(x=CB, y=ROB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = CBROB_gwave$coefficients[2], intercept = CBROB_gwave$coefficients[1],color="blue")+geom_abline(slope = CBROB_gpulse$coefficients[2], intercept = CBROB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(.4, 1.8));

rob_ob<-ggplot(gwavegpulse_avg, aes(x=ROB, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBOB_gwave$coefficients[2], intercept = ROBOB_gwave$coefficients[1],color="blue")+geom_abline(slope = ROBOB_gpulse$coefficients[2], intercept = ROBOB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(.4, 1.8), ylim = c(-.85, .2));
rob_tel<-ggplot(gwavegpulse_avg, aes(x=ROB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBTEL_gwave$coefficients[2], intercept = ROBTEL_gwave$coefficients[1],color="blue")+geom_abline(slope = ROBTEL_gpulse$coefficients[2], intercept = ROBTEL_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_ellhb<-ggplot(gwavegpulse_avg, aes(x=ROB, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBHB_gwave$coefficients[2], intercept = ROBHB_gwave$coefficients[1],color="blue")+geom_abline(slope = ROBHB_gpulse$coefficients[2], intercept = ROBHB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_ot<-ggplot(gwavegpulse_avg, aes(x=ROB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBOT_gwave$coefficients[2], intercept = ROBOT_gwave$coefficients[1],color="blue")+geom_abline(slope = ROBOT_gpulse$coefficients[2], intercept = ROBOT_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_ts<-ggplot(gwavegpulse_avg, aes(x=ROB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBTS_gwave$coefficients[2], intercept = ROBTS_gwave$coefficients[1],color="blue")+geom_abline(slope = ROBTS_gpulse$coefficients[2], intercept = ROBTS_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_cb<-ggplot(gwavegpulse_avg, aes(x=ROB, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBCB_gwave$coefficients[2], intercept = ROBCB_gwave$coefficients[1],color="blue")+geom_abline(slope = ROBCB_gpulse$coefficients[2], intercept = ROBCB_gpulse$coefficients[1],color="seagreen2",lty=1)	
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));
rob_blank<-ggplot(gwavegpulse_avg, aes(x=ROB, y=ROB)) + geom_point(color='white') + theme_classic(base_size=7) + labs(x="Log Region Volume", y = "Log Region Volume") 
#+ 
#	coord_cartesian(xlim = c(0.4, 1.8), ylim = c(0.4, 1.8));

## ob, tel, hb, ot, ts, cb, rob
all_regionxregion_plot_gwavegpulse<-ggarrange(ob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()), ncol=7, nrow=7, common.legend = TRUE, legend = "none");
annotate_figure(all_regionxregion_plot_gwavegpulse, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

## ob, ot, tel, rob, hb, ts, cb
all_regionxregion_plot_gwavegpulse2<-ggarrange(ob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()), ncol=7, nrow=7, common.legend = TRUE, legend = "none");
annotate_figure(all_regionxregion_plot_gwavegpulse2, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

## ob, ot, rob, tel, hb, ts, cb
all_regionxregion_plot_gwavegpulse3<-ggarrange(ob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ot_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),rob_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),tel_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ellhb_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),ts_cb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ot + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_rob + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_tel + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ellhb + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_ts + theme(axis.title.y = element_blank(),axis.title.x = element_blank()),cb_blank + theme(axis.title.y = element_blank(),axis.title.x = element_blank()), ncol=7, nrow=7, common.legend = TRUE, legend = "none");
annotate_figure(all_regionxregion_plot_gwavegpulse3, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

## saving plots to pdf
# Customizing the gpulseput
pdf("regionxregion_gwavegpulse_means_Feb3.pdf",         # File name
    width = 7.874, height = 7.874) # Width and height in inches

annotate_figure(all_regionxregion_plot_gwavegpulse3, left = text_grob(bquote('Log Region Volume' ~(mm^3)), rot = 90),
                    bottom = text_grob(bquote('Log Region Volume' ~(mm^3))))

# Closing the graphical device
dev.off() 


####################################
#   plots of regions against ROB   #
####################################

rob_ob2<-ggplot(comb_log_data, aes(x=RoB, y=OB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(title="Olfactory Bulbs", x="Log RoB Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBOB_tubegen$coefficients[2], intercept = ROBOB_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROBOB_ampegen$coefficients[2], intercept = ROBOB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROBOB_out$coefficients[2], intercept = ROBOB_out$coefficients[1],color="Black",lty=2)	
rob_tel2<-ggplot(comb_log_data, aes(x=RoB, y=TEL, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(title="Telencephalon", x="Log RoB Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBTEL_tubegen$coefficients[2], intercept = ROBTEL_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROBTEL_ampegen$coefficients[2], intercept = ROBTEL_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROBTEL_out$coefficients[2], intercept = ROBTEL_out$coefficients[1],color="Black",lty=2)	
rob_ellhb2<-ggplot(comb_log_data, aes(x=RoB, y=HB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(title="Hindbrain", x="Log RoB Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBHB_tubegen$coefficients[2], intercept = ROBHB_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROBHB_ampegen$coefficients[2], intercept = ROBHB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROBHB_out$coefficients[2], intercept = ROBHB_out$coefficients[1],color="Black",lty=2)	
rob_ot2<-ggplot(comb_log_data, aes(x=RoB, y=OT, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(title="Optic Tectum", x="Log RoB Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBOT_tubegen$coefficients[2], intercept = ROBOT_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROBOT_ampegen$coefficients[2], intercept = ROBOT_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROBOT_out$coefficients[2], intercept = ROBOT_out$coefficients[1],color="Black",lty=2)	
rob_ts2<-ggplot(comb_log_data, aes(x=RoB, y=TS, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(title="Torus Semicircularis", x="Log RoB Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBTS_tubegen$coefficients[2], intercept = ROBTS_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROBTS_ampegen$coefficients[2], intercept = ROBTS_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROBTS_out$coefficients[2], intercept = ROBTS_out$coefficients[1],color="Black",lty=2)	
rob_cb2<-ggplot(comb_log_data, aes(x=RoB, y=CB, shape=Species, color=Species, fill=Species)) + scale_shape_manual(values=spp_shapes) + scale_fill_manual(values=spp_fills) + scale_color_manual(values=sim_colors) + geom_point(size=1) + theme_classic(base_size=7) + labs(title="Cerebellum", x="Log RoB Volume", y = "Log Region Volume") + 
	geom_abline(slope = ROBCB_tubegen$coefficients[2], intercept = ROBCB_tubegen$coefficients[1],color="Black")+geom_abline(slope = ROBCB_ampegen$coefficients[2], intercept = ROBCB_ampegen$coefficients[1],color="Black",lty=3)+geom_abline(slope = ROBCB_out$coefficients[2], intercept = ROBCB_out$coefficients[1],color="Black",lty=2)	

region_rob_plot<-ggarrange(rob_ob2, rob_tel2, rob_ellhb2, rob_ot2, rob_ts2, rob_cb2, ncol=3, nrow=2, common.legend = TRUE, legend = "none")


#####################################################################
#																	#
#	           Meaningful hypothesis testing ancovas				#
#																	#
#####################################################################

##############################################################################################
# Hypothesis Electrogenic + Receptor Type -- tuberous + eo vs ampullary + eo vs not electric #
##############################################################################################

## Ancova
tree <- all_tree_sg

#tb <- TBV_M
#br <- OB_M 			## change this to the diff regions
#br <- TEL_M
#br <- HB_M
#br <- OT_M
#br <- TS_M
#br <- CB_M
#br <- RoB_M

## region x region - roi
#tb <- TBV_OB
#tb <- TBV_TEL
#tb <- TBV_HB
#tb <- TBV_OT
#tb <- TBV_TS
#tb <- TBV_CB
#tb <- TBV_ROB

## region x region, tb = x, br = y
#tb <- OB_M 			
#tb <- TEL_M
#tb <- HB_M
#tb <- OT_M
#tb <- TS_M
#tb <- CB_M
tb <- RoB_M

br <- OB_M 			
#br <- TEL_M
#br <- HB_M
#br <- OT_M
#br <- TS_M
#br <- CB_M
#br <- RoB_M

##phylogenetically comparing eorecpt groups
hyp_c <- c("tubegen","tubegen","tubegen","tubegen","tubegen","out","out","out","out","tubegen","tubegen","tubegen","tubegen","tubegen","out","tubegen","tubegen","out","out","tubegen","out","tubegen","out","tubegen","tubegen","tubegen","tubegen","ampegen","ampegen","ampegen","ampegen","out")
data <- data.frame(avg_comb_data$Species,tb,br,row.names=avg_comb_data$Species)
colnames(data)<-c("Species","tb","br")
## with phylogeny, 1st gls is slope, 2nd gls is intercept
correlation <- corPagel(.8,tree)
bm.gls <- gls(br ~ tb + hyp_c + tb:hyp_c, correlation = correlation, data=data, na.action = na.omit)
bm.gls2 <- gls(br ~ tb + hyp_c, correlation = correlation, data=data, na.action = na.omit)

## test for normality
# Get residuals for gls
reg.res<-resid(bm.gls)
# shaprio wilk normality test
shapiro.test(reg.res)
# levene test for homogeneity of variances, looks at the residuals by group
leveneTest(reg.res ~ hyp_c)
## test for heteroskedasticy of gls by looking at standardized residuals vs fitted values

## Ancovas
## anova = type I, where the order of the variables in the model matters w/ first variable explaining the most variation
anova(bm.gls) #p value of x:type shows slope
anova(bm.gls2)

## posthoc test
pairs(emmeans(bm.gls, "hyp_c"))
pairs(emmeans(bm.gls2, "hyp_c")) ## does a tukey pairwise comparison!

## effect size
eff_size(emmeans(bm.gls2, "hyp_c"), sigma=sigma(bm.gls2), edf=28)
eff_size(emmeans(bm.gls, "hyp_c"), sigma=sigma(bm.gls), edf=26)

## for pairwise slope 
pairs(emtrends(bm.gls,"hyp_c",var="tb",mode = "df.error"))


##############################################
# Hypothesis Egen + Tub -- morms vs gymnotes #
##############################################

## Ancova
tree <- tubegen.tree

#tb <- TBV_M[tubegen_names]
#br <- OB_M[tubegen_names] 			## change this to the diff regions
#br <- TEL_M[tubegen_names]
#br <- HB_M[tubegen_names]
#br <- OT_M[tubegen_names]
#br <- TS_M[tubegen_names]
#br <- CB_M[tubegen_names]
br <- RoB_M[tubegen_names]

## region x region - roi
#tb <- TBV_OB[tubegen_names]
#tb <- TBV_TEL[tubegen_names]
#tb <- TBV_HB[tubegen_names]
#tb <- TBV_OT[tubegen_names]
#tb <- TBV_TS[tubegen_names]
#tb <- TBV_CB[tubegen_names]
tb <- TBV_ROB[tubegen_names]

## region x region, tb = x, br = y
#tb <- OB_M[tubegen_names] 			
#tb <- TEL_M[tubegen_names]
#tb <- HB_M[tubegen_names]
#tb <- OT_M[tubegen_names]
#tb <- TS_M[tubegen_names]
#tb <- CB_M[tubegen_names]
#tb <- RoB_M[tubegen_names]

#br <- OB_M[tubegen_names]
#br <- TEL_M[tubegen_names]
#br <- HB_M[tubegen_names]
#br <- OT_M[tubegen_names]
#br <- TS_M[tubegen_names]
#br <- CB_M[tubegen_names]
#br <- RoB_M[tubegen_names]

##phylogenetically comparing eorecpt groups
hyp_c <- c("gym","gym","morm","morm","morm","gym","gym","gym","morm","morm","gym","gym","morm","morm","gym","gym","gym","gym")
data <- data.frame(tubegen_avg$Species,tb,br,row.names=tubegen_avg$Species)
colnames(data)<-c("Species","tb","br")
## with phylogeny, 1st gls is slope, 2nd gls is intercept
correlation <- corPagel(.8,tree)
bm.gls <- gls(br ~ tb + hyp_c + tb:hyp_c, correlation = correlation, data=data, na.action = na.omit)
bm.gls2 <- gls(br ~ tb + hyp_c, correlation = correlation, data=data, na.action = na.omit)

## test for normality
# Get residuals for gls
reg.res<-resid(bm.gls)
# shaprio wilk normality test
shapiro.test(reg.res)
# levene test for homogeneity of variances, looks at the residuals by group
leveneTest(reg.res ~ hyp_c)
## test for heteroskedasticy of gls by looking at standardized residuals vs fitted values

## Ancovas
## do type III ancova if there is a sig interaction, do type II ancova if there is not sig interaction
## anova = type I, where the order of the variables in the model matters w/ first variable explaining the most variation
anova(bm.gls) #p value of x:type shows slope
anova(bm.gls2)

## effect size
eff_size(emmeans(bm.gls2, "hyp_c"), sigma=sigma(bm.gls2), edf=15)
eff_size(emmeans(bm.gls, "hyp_c"), sigma=sigma(bm.gls), edf=14)


###############################################################
# Hypothesis Egen + Tub -- morms vs gymnotes, excl gymnarchus #
###############################################################

## Ancova
nogymtubegen_avg<-subset(tubegen_avg, Species != "Gymnarchus_niloticus")
nogymtubegen_names<-as.character(nogymtubegen_avg$Species)
nogymtubegen.tree<-drop.tip(all_tree_sg, setdiff(all_tree_sg$tip.label, nogymtubegen_names))
tree <- nogymtubegen.tree

#tb <- TBV_M[nogymtubegen_names]
#br <- OB_M[nogymtubegen_names] 			## change this to the diff regions
#br <- TEL_M[nogymtubegen_names]
#br <- HB_M[nogymtubegen_names]
#br <- OT_M[nogymtubegen_names]
#br <- TS_M[nogymtubegen_names]
#br <- CB_M[nogymtubegen_names]
#br <- RoB_M[nogymtubegen_names]

## region x region - roi
#tb <- TBV_OB[nogymtubegen_names]
#tb <- TBV_TEL[nogymtubegen_names]
#tb <- TBV_HB[nogymtubegen_names]
#tb <- TBV_OT[nogymtubegen_names]
#tb <- TBV_TS[nogymtubegen_names]
#tb <- TBV_CB[nogymtubegen_names]
#tb <- TBV_ROB[nogymtubegen_names]

## region x region, tb = x, br = y
#tb <- OB_M[nogymtubegen_names]
#tb <- TEL_M[nogymtubegen_names]
#tb <- HB_M[nogymtubegen_names]
#tb <- OT_M[nogymtubegen_names]
#tb <- TS_M[nogymtubegen_names]
#tb <- CB_M[nogymtubegen_names]
tb <- RoB_M[nogymtubegen_names]

#br <- OB_M[nogymtubegen_names]
#br <- TEL_M[nogymtubegen_names]
#br <- HB_M[nogymtubegen_names]
#br <- OT_M[nogymtubegen_names]
#br <- TS_M[nogymtubegen_names]
br <- CB_M[nogymtubegen_names]
#br <- RoB_M[nogymtubegen_names]

##phylogenetically comparing eorecpt groups
hyp_c <- c("gym","gym","morm","morm","morm","gym","gym","gym","morm","gym","gym","morm","morm","gym","gym","gym","gym")
data <- data.frame(nogymtubegen_avg$Species,tb,br,row.names=nogymtubegen_avg$Species)
colnames(data)<-c("Species","tb","br")
## with phylogeny, 1st gls is slope, 2nd gls is intercept
correlation <- corPagel(.8,tree)
bm.gls <- gls(br ~ tb + hyp_c + tb:hyp_c, correlation = correlation, data=data, na.action = na.omit)
bm.gls2 <- gls(br ~ tb + hyp_c, correlation = correlation, data=data, na.action = na.omit)

## test for normality
# Get residuals for gls
reg.res<-resid(bm.gls)
# shaprio wilk normality test
shapiro.test(reg.res)
# levene test for homogeneity of variances, looks at the residuals by group
leveneTest(reg.res ~ hyp_c)
## test for heteroskedasticy of gls by looking at standardized residuals vs fitted values

## Ancovas
## anova = type I, where the order of the variables in the model matters w/ first variable explaining the most variation
anova(bm.gls) #p value of x:type shows slope
anova(bm.gls2)

## effect size
eff_size(emmeans(bm.gls2, "hyp_c"), sigma=sigma(bm.gls2), edf=15)
eff_size(emmeans(bm.gls, "hyp_c"), sigma=sigma(bm.gls), edf=14)


########################################
# Hypothesis Wave vs Pulse -- gymnotes #
########################################

## Ancova
tree <- gym.tree

#tb <- TBV_M[gym_names]
#br <- OB_M[gym_names] 			## change this to the diff regions
#br <- TEL_M[gym_names]
#br <- HB_M[gym_names]
#br <- OT_M[gym_names]
#br <- TS_M[gym_names]
#br <- CB_M[gym_names]
br <- RoB_M[gym_names]

## region x region - roi
#tb <- TBV_OB[gym_names]
#tb <- TBV_TEL[gym_names]
#tb <- TBV_HB[gym_names]
#tb <- TBV_OT[gym_names]
#tb <- TBV_TS[gym_names]
#tb <- TBV_CB[gym_names]
tb <- TBV_ROB[gym_names]

## region x region, tb = x, br = y
#tb <- OB_M[gym_names]
#tb <- TEL_M[gym_names]
#tb <- HB_M[gym_names]
#tb <- OT_M[gym_names]
#tb <- TS_M[gym_names]
#tb <- CB_M[gym_names]
#tb <- RoB_M[gym_names]

#br <- OB_M[gym_names]
#br <- TEL_M[gym_names]
#br <- HB_M[gym_names]
#br <- OT_M[gym_names]
#br <- TS_M[gym_names]
#br <- CB_M[gym_names]
#br <- RoB_M[gym_names]

##phylogenetically comparing eorecpt groups
hyp_c <- c("wave","pulse","wave","wave","pulse","pulse","pulse","pulse","wave","wave","wave")
data <- data.frame(gym_avg$Species,tb,br,row.names=gym_avg$Species)
colnames(data)<-c("Species","tb","br")
## with phylogeny, 1st gls is slope, 2nd gls is intercept
correlation <- corPagel(.8,tree)
bm.gls <- gls(br ~ tb + hyp_c + tb:hyp_c, correlation = correlation, data=data, na.action = na.omit)
bm.gls2 <- gls(br ~ tb + hyp_c, correlation = correlation, data=data, na.action = na.omit)

## test for normality
# Get residuals for gls
reg.res<-resid(bm.gls)
# shaprio wilk normality test
shapiro.test(reg.res)
# levene test for homogeneity of variances, looks at the residuals by group
leveneTest(reg.res ~ hyp_c)
## test for heteroskedasticy of gls by looking at standardized residuals vs fitted values

## Ancovas
## anova = type I, where the order of the variables in the model matters w/ first variable explaining the most variation
anova(bm.gls) #p value of x:type shows slope
anova(bm.gls2)

## effect size
eff_size(emmeans(bm.gls2, "hyp_c"), sigma=sigma(bm.gls2), edf=8)
eff_size(emmeans(bm.gls, "hyp_c"), sigma=sigma(bm.gls), edf=7)


######################################################################
# Hypothesis All Out -- morm outs vs gym outs (includ not egen cats) #
######################################################################

## Ancova
tree <- out.tree

#tb <- TBV_M[out_names]
#br <- OB_M[out_names] 			## change this to the diff regions
#br <- TEL_M[out_names]
#br <- HB_M[out_names]
#br <- OT_M[out_names]
#br <- TS_M[out_names]
#br <- CB_M[out_names]
br <- RoB_M[out_names]

## region x region - roi
#tb <- TBV_OB[out_names]
#tb <- TBV_TEL[out_names]
#tb <- TBV_HB[out_names]
#tb <- TBV_OT[out_names]
#tb <- TBV_TS[out_names]
#tb <- TBV_CB[out_names]
tb <- TBV_ROB[out_names]

## region x region, tb = x, br = y
#tb <- OB_M[out_names]
#tb <- TEL_M[out_names]
#tb <- HB_M[out_names]
#tb <- OT_M[out_names]
#tb <- TS_M[out_names]
#tb <- CB_M[out_names]
#tb <- RoB_M[out_names]

#br <- OB_M[out_names]
#br <- TEL_M[out_names]
#br <- HB_M[out_names]
#br <- OT_M[out_names]
#br <- TS_M[out_names]
#br <- CB_M[out_names]
#br <- RoB_M[out_names]

##phylogenetically comparing eorecpt groups
hyp_c <- c("mout","gout","gout","gout","gout","gout","gout","mout","gout","mout")
data <- data.frame(out_avg$Species,tb,br,row.names=out_avg$Species)
colnames(data)<-c("Species","tb","br")
## with phylogeny, 1st gls is slope, 2nd gls is intercept
correlation <- corPagel(.8,tree)
bm.gls <- gls(br ~ tb + hyp_c + tb:hyp_c, correlation = correlation, data=data, na.action = na.omit)
bm.gls2 <- gls(br ~ tb + hyp_c, correlation = correlation, data=data, na.action = na.omit)

## test for normality
# Get residuals for gls
reg.res<-resid(bm.gls)
# shaprio wilk normality test
shapiro.test(reg.res)
# levene test for homogeneity of variances, looks at the residuals by group
leveneTest(reg.res ~ hyp_c)
## test for heteroskedasticy of gls by looking at standardized residuals vs fitted values

## Ancovas
## anova = type I, where the order of the variables in the model matters w/ first variable explaining the most variation
anova(bm.gls) #p value of x:type shows slope
anova(bm.gls2)

## effect size
eff_size(emmeans(bm.gls2, "hyp_c"), sigma=sigma(bm.gls2), edf=7)
eff_size(emmeans(bm.gls, "hyp_c"), sigma=sigma(bm.gls), edf=6)
