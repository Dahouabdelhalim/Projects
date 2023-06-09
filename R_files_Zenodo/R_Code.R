## LOAD PACKAGES ##
#General
library(tidyverse)
library(magrittr)
library(stringr)
library(readr)
library(reshape2)
## Biological
library(Rphylip)
library(arbutus)
library(vegan)
library(ape)
library(phangorn)
library(phytools)
library(OUwie)
library(picante)
library(geiger)
library(phylobase)
library(fields)
library(phylosignal)
library(geomorph)
library(phylopath)
library(phylolm)
library(adegenet)
## Graphics
library(FactoMineR)
library(factoextra)
library(corrplot)
library(BAMMtools)
library(gridExtra)
library(colorRamps)
library(PerformanceAnalytics)

# Set paths to input data
#setwd("~/tentilla_morph/Supplementary_Materials/Dryad/")

#Load raw data
read.csv("raw_morphology_data.csv") -> numbers
numbers$Species = as.character(numbers$Species)
categorical <- read.csv("raw_categorical_data.csv")[,-2]
rownames(categorical) = categorical$Species

#Correct species spellings
numbers$Species[which(numbers$Species == "Agalma okeni")] <- "Agalma okenii"
numbers$Species[which(numbers$Species == "Forskalia edwardsi")] <- "Forskalia edwardsii"
numbers$Species[which(numbers$Species == "Rhizophysa eysenhardti")] <- "Rhizophysa eysenhardtii"

#Merge Nanomias
numbers$Species[which(numbers$Species == "Nanomia cara")] <- "Nanomia sp"
numbers$Species[which(numbers$Species == "Nanomia bijuga")] <- "Nanomia sp"

#Polish fields
castnumbers = numbers
castnumbers[castnumbers=="needDIC"] <- NA
castnumbers[castnumbers=="needConfocal"] <- NA
castnumbers[castnumbers=="needMature"] <- NA
castnumbers[castnumbers=="needReslide"] <- NA
castnumbers[castnumbers==-1] <- NA
castnumbers[,c(-1,-2,-3)] <- sapply(castnumbers[,c(-1,-2,-3)],as.character)
castnumbers[,c(-1,-2,-3)] <- sapply(castnumbers[,c(-1,-2,-3)],as.numeric)

#Define functions to deal with NAs and decimals
mean.na <- function(a){mean(a, na.rm = TRUE)}
var.na <- function(a){var(a, na.rm = TRUE)}
round3 <- function(a){round(a, 3)}

#Compute morphometric ratios
coiledness = castnumbers$Cnidoband.free.length..um./castnumbers$Cnidoband.length..um.
heteroneme_elongation = castnumbers$Heteroneme.free.length..um./castnumbers$Heteroneme.width..um.
haploneme_elongation = castnumbers$Haploneme.free.length..um./castnumbers$Haploneme.width..um.
desmoneme_elongation = castnumbers$Desmoneme.length..um./castnumbers$Desmoneme.width..um.
rhopaloneme_elongation = castnumbers$Rhopaloneme.length..um./castnumbers$Rhopaloneme.width..um.
heteroneme_shaft_extension = castnumbers$Heteroneme.free.length..um./castnumbers$Heteroneme.shaft.free.length..um.
Heteroneme_to_CB = castnumbers$Heteroneme.free.length..um./(castnumbers$Cnidoband.free.length..um.+0.001)
total_heteroneme_volume = castnumbers$Heteroneme.volume..um3.*castnumbers$Heteroneme.number
total_haploneme_volume = (castnumbers$Cnidoband.free.length..um./castnumbers$Haploneme.width..um.)*((4*pi/3)*0.5*castnumbers$Haploneme.free.length..um.*((0.5*castnumbers$Haploneme.width..um.)^2))
cnidomic_index = log(total_heteroneme_volume + total_haploneme_volume)
cnidomic_index[which(is.na(cnidomic_index))] <- log(total_haploneme_volume[which(is.na(cnidomic_index))])
cnidomic_index[which(is.na(cnidomic_index))] <- log(castnumbers[which(is.na(cnidomic_index)), "Heteroneme.volume..um3."])
cnidomic_index[cnidomic_index==-Inf]<-NA

#Define surface and volume functions
surface_ellipsoid <- function(L, W){
  a <- L/2
  b <- W/2
  c = b
  S <- 4*pi*((((a*b)^1.6)+((a*c)^1.6)+((c*b)^1.6))/3)^(1/1.6)
  return(S)
}
volume_ellipsoid <- function(L, W){
  a <- L/2
  b <- W/2
  c = b
  V <- (4/3)*pi*a*b*c
  return(V)
}

SAV_haploneme = surface_ellipsoid(castnumbers$Haploneme.free.length..um., castnumbers$Haploneme.width..um.)/volume_ellipsoid(castnumbers$Haploneme.free.length..um., castnumbers$Haploneme.width..um.)
names(SAV_haploneme) = castnumbers$Species
SAV_heteroneme = surface_ellipsoid(castnumbers$Heteroneme.free.length..um., castnumbers$Heteroneme.width..um.)/volume_ellipsoid(castnumbers$Heteroneme.free.length..um., castnumbers$Heteroneme.width..um.)
names(SAV_heteroneme) = castnumbers$Species

#Combine all morphometric variables
morphometrics = data.frame(castnumbers$slide_id, castnumbers$Species,coiledness, heteroneme_elongation, haploneme_elongation, desmoneme_elongation, rhopaloneme_elongation,heteroneme_shaft_extension, Heteroneme_to_CB, total_heteroneme_volume, total_haploneme_volume, cnidomic_index, SAV_haploneme)
names(morphometrics)[1:2] = c("slide_id"," Species")

#Combine morphometrics with regular characters
castnumbers = data.frame(castnumbers,morphometrics[,c(-1,-2)])

#Logtransform variables which are not normal
castlogs = castnumbers
non_normal = apply(castnumbers[,c(-1,-2,-3)], 2, shapiro.test) %>% lapply(function(x){x$p.value}) %>% unlist() 
non_normal = which(as.vector(non_normal) < 0.05)+3
castlogs[,non_normal] <- sapply(castnumbers[,non_normal], log)
castlogs[castlogs==-Inf]<-NA

#Means, varinces, standard errors
castmeans <- aggregate(. ~  Species, data = castnumbers[,c(-1,-3)], mean.na, na.action = na.pass)
castmean_logs <- aggregate(. ~  Species, data = castlogs[,c(-1,-3)], mean.na, na.action = na.pass)
castvariances <- aggregate(. ~  Species, data = castnumbers[,c(-1,-3)], FUN = var.na, na.action = na.pass) #[,-2]
castvariance_logs <- aggregate(. ~  Species, data = castlogs[,c(-1,-3)], FUN = var.na, na.action = na.pass)
castvariances[is.na(castvariances)] <- 0
castvariance_logs[is.na(castvariance_logs)] <- 0
castN <- aggregate(. ~  Species, data = castnumbers[,c(-1,-3)], length, na.action = na.pass) 
castN = data.frame(castN$Species, rowMeans(castN[,-1])) #[,-2]
names(castN) = c("Species", "N_specimens")
castSE <- cbind(castvariances$Species, sqrt(castvariances[,-1])/sqrt(castN$N_specimens))
castSE_logs <- cbind(castvariance_logs$Species, sqrt(castvariance_logs[,-1])/sqrt(castN[,-1]))
names(castSE)[1] <- "Species"
names(castSE_logs)[1] <- "Species"

#Make overview of data and heatmap
Ses = castSE[,-1]
Ses[Ses == 0] = NA
morphdata = cbind(castN, castmeans[,-1], Ses)
morphdata = morphdata[,c(1:3,33,4,34,5,35,6,36,7,37,8,38,9,39,10,40,11,41,12,42,13,43,14,44,15,45,16,46,17,47,18,48,19,49,20,50,21,51,22,52,23,53,24,54,25,55,26,56,27,57,28,58,29,59,30,60,31,61,32,62)]
names(morphdata) = str_replace_all(names(morphdata),".1","_SE")
#write.csv(morphdata, "characterdata.csv")

heatdata = as.matrix(castmean_logs[,-1])
rownames(heatdata) = castmean_logs$Species
heatdata[is.nan(heatdata)]<- -1
#hcolors = grDevices::terrain.colors(20)
library(colorspace)
sequential_hcl()
hcolors <- colorRampPalette(c("#e0ecf4","#9ebcda","#8856a7"), bias=1)
hcolors <- colorRampPalette(c("#ffeda0","#feb24c","#f03b20"), bias=1)
hcolors <- colorRampPalette(c("yellow","purple"), bias=1)
#hcolors <- sequential_hcl(17, h = 245, c = c(40, 75, 0), l = c(30, 95), power = 1)
hcolors <- c(rep("#000000FF",10), hcolors(30))
heatmap(heatdata, scale = "column", cexCol = 0.2, col=hcolors, keep.dendro = T)

#Load phylogenetic tree
consensus = read.nexus("../TreeBase/RB_constrained_timetree/TimeTree_siphs_mcmc_MAP.tre") %>% drop.tip(56:61)
consensus$tip.label = str_replace_all(consensus$tip.label,"_"," ")
#Switch Nanomia bijuga for Nanaomia sp
consensus$tip.label[which(consensus$tip.label == "Nanomia bijuga")] <- "Nanomia sp"

#Prune quant matrix to tree species
matrix = castnumbers[which(!is.na(castnumbers$Species)),]
matrix_logs = castlogs[which(!is.na(castlogs$Species)),]
sharedspp = matrix$Species[which(matrix$Species %in% consensus$tip.label)]
sharedmatrix = matrix[which(matrix$Species %in% sharedspp),]
sharedlogs = matrix_logs[which(matrix$Species %in% sharedspp),]
sharedmeans = castmeans[which(castmeans$Species %in% sharedspp),]
sharedmean_logs = castmean_logs[which(castmean_logs$Species %in% sharedspp),]
sharedvars = castvariances[which(castvariances$Species %in% sharedspp),]
sharedvar_logs = castvariance_logs[which(castvariance_logs$Species %in% sharedspp),]

#Prune categorical matrix to tree species
cat_rowNAs = apply(categorical[,-1], 1, function(x) sum(is.na(x)))
catmatrix = categorical[which(cat_rowNAs<1),] %>% .[,-1]
sharedspp_cat = rownames(catmatrix)[which(rownames(catmatrix) %in% consensus$tip.label)]
sharedcategorical = catmatrix[which(rownames(catmatrix) %in% sharedspp_cat),]

#Prune tree to quant matrix species
nodatatipnames = consensus$tip.label[which(!(consensus$tip.label %in% sharedmatrix$Species))]
nodatatips = c(1:length(consensus$tip.label))[which(consensus$tip.label %in% nodatatipnames)]
prunedtree = drop.tip(consensus, nodatatips)
prunedtree$tip.label
plot(prunedtree)
#ultram = chronos(prunedtree)
ultram = prunedtree
plot(ultram)

#Prune tree to categorical matrix species
nodatatipnames_cat = consensus$tip.label[which(!(consensus$tip.label %in% sharedspp_cat))]
nodatatips_cat = c(1:length(consensus$tip.label))[which(consensus$tip.label %in% nodatatipnames_cat)]
cat_tree = drop.tip(consensus, nodatatips_cat)
cat_tree$tip.label
plot(cat_tree)
#ultram_cat = chronos(cat_tree)
ultram_cat = cat_tree
plot(ultram_cat)
sharedcategorical = sharedcategorical[match(ultram_cat$tip.label, rownames(sharedcategorical)),]
cprunedmatrix = sharedcategorical[which(rownames(sharedcategorical)%in%rownames(sharedmatrix)),] %>% .[which(sapply(.,function(x) length(unique(x)))>1)]
sharedbinary = sharedcategorical
sharedbinary$Haploneme.type = as.character(sharedbinary$Haploneme.type)
sharedbinary[sharedbinary=="Isorhizas"] = 0
sharedbinary[sharedbinary=="Anisorhizas"] = 1
sharedbinary$Haploneme.type = as.numeric(sharedbinary$Haploneme.type)
sharedbinary$Heteroneme.type = as.character(sharedbinary$Heteroneme.type)
sharedbinary[sharedbinary=="Stenotele"] = 0
sharedbinary[sharedbinary=="Microbasic mastigophore"] = 1
sharedbinary[sharedbinary=="Eurytele"] = NA
sharedbinary$Heteroneme.type = as.numeric(sharedbinary$Heteroneme.type)

#Purely numerical character sets without slide or species
Q_sharedmatrix = sharedmatrix[,c(-1,-2,-3)]
rownames(sharedmeans) = sharedmeans$Species
Q_sharedmeans = sharedmeans[,-1]
rownames(sharedvars) = sharedvars$Species
Q_sharedvars = sharedvars[,-1]
Q_sharedmean_logs = sharedmean_logs[,-1]
rownames(Q_sharedmean_logs) = sharedmean_logs$Species

### COMPARATIVE ANALYSES ###

phylosignals = as.data.frame(matrix(ncol=3, nrow=ncol(sharedmean_logs[,-1])))
#phylosignals = as.data.frame(matrix(ncol=3, nrow=ncol(sharedmeans[,-1])))
names(phylosignals) = c("K", "P", "Ntaxa")
for(i in 2:ncol(sharedmean_logs)){
  CH_I=as.numeric(sharedmean_logs[,i])
  #CH_I=as.numeric(sharedmeans[,i])
  names(CH_I) = sharedmean_logs$Species
  #names(CH_I) = sharedmeans$Species
  SE_I = as.numeric(castSE[,i]) %>% log()
  SE_I[which(SE_I==-Inf)] = 0
  names(SE_I) = castSE$Species
  CH_I= CH_I[!is.na(CH_I)]
  SE_I= SE_I[!is.na(SE_I)]
  SE_I = SE_I[which(names(SE_I)%in%names(CH_I))]
  CH_I = CH_I[which(names(CH_I)%in%names(SE_I))]
  treeI = drop.tip(ultram,which(!(ultram$tip.label %in% names(CH_I))))
  class(treeI) = "phylo"
  rownames(phylosignals)[i-1] <- names(sharedmean_logs)[i]
  PSIG <- phylosig(treeI, CH_I, se = SE_I, test=T)
  phylosignals[i-1,1] <- PSIG$K
  phylosignals[i-1,2] <- PSIG$P
  phylosignals[i-1,3] <- length(CH_I)
  phylosignals$K = round(phylosignals$K,3)
}
#write.csv(phylosignals, "PhylosignalsWSE_log.csv")

#Model support
AICdf = as.data.frame(matrix(ncol=6,nrow=ncol(Q_sharedmean_logs)))
colnames(AICdf) = c("Variable", "white_noise", "starBM", "BM", "EB", "OU")
for(c in 1:ncol(Q_sharedmean_logs)){
  startree <- rescale(ultram, "lambda", 0)
  C = Q_sharedmean_logs[,c]
  names(C) = rownames(Q_sharedmean_logs)
  C = C[!is.na(C)]
  Ctree = drop.tip(ultram, which(!(ultram$tip.label %in% names(C))))
  startree = drop.tip(startree, which(!(startree$tip.label %in% names(C))))
  Cse = castSE[,c+1] %>% log() %>% abs()
  Cse[which(Cse == Inf)] <- 0
  names(Cse) = castSE$Species
  Cse = Cse[which(names(Cse) %in% names(C))]
  model_matrix = matrix("NA", nrow = 5, ncol = 3)
  colnames(model_matrix) = c("aicc","aicc_best","dAICc")
  row.names(model_matrix) = c("white", "starBM", "BM", "EB", "OU")
  for(j in 1:dim(model_matrix)[1]){
    if(j==2){
      temp_model = fitContinuous(startree, C, model="BM", SE = Cse)$opt
    }
    else{
      temp_model = fitContinuous(Ctree, C, model=row.names(model_matrix)[j], SE = Cse)$opt
    }
    model_matrix = apply(model_matrix,2, as.numeric)
    row.names(model_matrix) = c("white", "starBM", "BM", "EB", "OU")
    model_matrix[j, "aicc"] <- temp_model$aicc
  }
  model_matrix[,"aicc_best"] <- min(model_matrix[,"aicc"])
  model_matrix[,"dAICc"] <- model_matrix[, "aicc"] - model_matrix[j, "aicc_best"]
  print(names(Q_sharedmean_logs)[c])
  string_c <- c(names(Q_sharedmean_logs)[c], model_matrix[,3])
  names(string_c) = colnames(AICdf)
  AICdf[c,] <- string_c
}
AICdf[,2:6] = apply(AICdf[,2:6], 2, as.numeric) %>% apply(2, round3)
#write.csv(AICdf, "log_model_support.csv")

#Model adequacy
worthy_models = AICdf[which(AICdf$white_noise != 0 & AICdf$starBM != 0),]
MAD = as.data.frame(matrix(ncol = 6, nrow = nrow(worthy_models)))
names(MAD) = c("msig", "cvar", "svar", "sasr", "shgt", "dcfd")
rownames(MAD) = worthy_models$Variable
for(m in 1:nrow(worthy_models)){
  C = Q_sharedmean_logs[,which(names(Q_sharedmean_logs) == worthy_models$Variable[m])]
  names(C) = rownames(Q_sharedmean_logs)
  C = C[!is.na(C)]
  Ctree = drop.tip(ultram, which(!(ultram$tip.label %in% names(C))))
  C = C[match(Ctree$tip.label, names(C))]
  class(Ctree)="phylo"
  FC <- fitContinuous(Ctree,C,model=worthy_models$Best_model[m])
  UTC <- make_unit_tree(FC)
  picstat_data <- calculate_pic_stat(UTC)
  sim <- simulate_char_unit(UTC)
  picstat_sim <- calculate_pic_stat(sim)
  compare_pic_stat(picstat_data, picstat_sim) %>% .$p.values -> MAD[m,]
}
phy_models = data.frame(worthy_models,MAD)
phy_models[,7:12] = apply(phy_models[,7:12], 2,round3)
#write.csv(phy_models, "model_adequacy.csv")

#Nematocyst shape evolution - phylomorphospace figure
het_el = sharedmean_logs$heteroneme_elongation
names(het_el) = sharedmean_logs$Species
hap_el = sharedmean_logs$haploneme_elongation
names(hap_el) = sharedmean_logs$Species
het_el = het_el[which(!is.na(het_el))]
hap_el = hap_el[which(!is.na(hap_el))]
hap_el = hap_el[which(names(hap_el) %in% names(het_el))]
het_el = het_el[which(names(het_el) %in% names(hap_el))]
het_el = het_el[match(names(hap_el), names(het_el))]
elon_data = cbind(het_el, hap_el) %>% as.data.frame()
names(elon_data) = c("Heteroneme elongation (log um)", "Haploneme elongation (log um)")
elon_tree = drop.tip(ultram, which(!(ultram$tip.label %in% names(hap_el))))
phylomorphospace(elon_tree, elon_data, label="horizontal")

## Character correlations ##
C = sharedmatrix[,c(-1,-3)] %>% .[which(!is.na(rowSums(.[,-1]))),]
Cspecies = C$Species %>% as.character()
C = as.matrix(C[,-1])
rownames(C) = Cspecies
Ctree = drop.tip(ultram,which(!(ultram$tip.label %in% rownames(C))))
class(Ctree) = "phylo"
PICi = Rcontrast(tree=Ctree, X=C, path="phylip-3.695/exe", cleanup=TRUE)

#Correlation visualizations R2
phy_corr = PICi$VarA.Correlations
intra_corr = PICi$VarE.Correlations
phy_corr[upper.tri(phy_corr)] <- NA
combicorr = phy_corr
combicorr[upper.tri(combicorr)] = intra_corr[upper.tri(intra_corr)]
rownames(combicorr) = colnames(C)
colnames(combicorr) = colnames(C)
phylointramatrix <-corrplot(combicorr, diag=F, tl.cex = 0.4, tl.col="black") #phylo correlations vs intraspecific correlations
PICOLS = combicorr
PICOLS[upper.tri(PICOLS)]=cor(C)[upper.tri(cor(C))]
phyloregmatrix <- corrplot(PICOLS, diag=F, tl.cex = 0.4, tl.col="black") #phylo correlations vs regular correlations for figure

#Scatterplot phylo vs regular correlations
cbind(as.vector(phy_corr[lower.tri(phy_corr)]), as.vector(cor(C)[lower.tri(cor(C))])) %>% as.data.frame()->phyreg
names(phyreg)<-c("Phylo", "Reg")
ggplot(phyreg, aes(x=Reg, y=Phylo, color=(Phylo+Reg)/2)) + geom_point() + geom_hline(yintercept = 0)  + geom_vline(xintercept = 0) + theme_bw()
abline(h=0)
abline(v=0)

phylomatrix = phyloregmatrix
phylomatrix[upper.tri(phylomatrix)]<-NA
melt(phylomatrix) %>% .[which(!is.na(.[,3])),] -> phylocolumn
names(phylocolumn)<-c("CH1","CH2","R2_phylo")
regmatrix = phyloregmatrix
regmatrix[lower.tri(regmatrix)]<-NA
melt(regmatrix) %>% .[which(!is.na(.[,3])),] -> regcolumn
regcolumn = regcolumn[,c(2,1,3)]
names(regcolumn)<-c("CH1","CH2","R2_regular")
phyreg_full = left_join(regcolumn,phylocolumn)
opposites = phyreg_full[which(phyreg_full[,3]*phyreg_full[,4] < 0),]
opposites[which(opposites$R2_phylo<0),] %>% .[order(.$R2_regular-.$R2_phylo),]
opposites[which(opposites$R2_regular<0),] %>% .[order(.$R2_regular-.$R2_phylo),]
#phyreg_full[which(phyreg_full$R2_regular>0.5 & phyreg_full$R2_phylo>0 & phyreg_full$R2_phylo<0.3),] %>% .[order(.$R2_regular-.$R2_phylo),]

## Retrieve diet data ##
#Retrieve diet info from literature BINARY
GC = read.csv("literature_diet_data.tsv", header = T, sep='\\t')
GC$Prey.type = factor(GC$Prey.type, levels=unique(GC$Prey.type))
GC$Siphonophore.species = as.character(GC$Siphonophore.species)
#Fix typos#
GC$Siphonophore.species[which(GC$Siphonophore.species == "Nanomia bijuga")] <- "Nanomia sp"
GC$Siphonophore.species[which(GC$Siphonophore.species == "Rhizophysa eyesenhardti")] <- "Rhizophysa eysenhardtii"
GC$Siphonophore.species[which(GC$Siphonophore.species == "Agalma okeni")] <- "Agalma okenii"
GC = GC[,1:2]
names(GC)<-c("species", "character")
#write.csv(GC, "gutcontentliteraturereview.csv")

#Prune to tree species#
GC = GC[which(GC$species %in% ultram$tip.label),]
GC = split(GC,GC$character)
nrowGC = purrr::map(GC,nrow) %>% as.numeric()
GC = GC[which(nrowGC>2)]
GC = purrr::map(GC,unique)
diet= matrix(ncol=length(GC),nrow=length(unique(sharedmeans$Species))) %>% as.data.frame()
names(diet) = names(GC)
rownames(diet) = unique(sharedmeans$Species)
for(E in GC){
  print(E$species)
  for(S in E$species){
    diet[which(rownames(diet) == S),as.character(unique(E$character))] = 1
  }
}
diet[is.na(diet)] <- 0
diet = diet[which(rowSums(diet)>0),which(colSums(diet)<nrow(diet))]

#Add personal observations of the authors
cladeB <- matrix(rep(c(0,1,0,0,0,0,0,0,0,0),3),nrow=10, ncol = 3) %>% t() %>% as.data.frame()
rownames(cladeB) = c("Erenna richardi", "Erenna sirena", "Stephanomia amphytridis")
krilleaters = matrix(rep(c(0,0,0,1,0,0,1,0,0,0),4),nrow=10, ncol = 4) %>% t() %>% as.data.frame()
rownames(krilleaters) = c("Praya dubia", "Resomia ornicephala", "Lychnagalma utricularia", "Bargmannia amoena")
Gelatinous_diet = rep(0,nrow(diet))
diet = cbind(diet, Gelatinous_diet)
names(diet)[ncol(diet)] = "Gelatinous diet"
gelateaters = matrix(rep(c(0,0,0,0,0,0,0,0,0,1),1),nrow=10, ncol = 1) %>% t() %>% as.data.frame()
names(gelateaters) = names(diet)
rownames(gelateaters) = c("Apolemia rubriversa")
names(gelateaters)=names(diet)
names(krilleaters)=names(diet)
names(cladeB)=names(diet)
diet = rbind(diet, cladeB, krilleaters, gelateaters)
#Decapod diet column actually encompasses decapods, krill, and mysids (large crustaceans/shrimp like animals)

# Retrive ROV annotation data
VARS <- read.csv("raw_ROV_data.csv")
VARS_curated = VARS[which(VARS$Siphonophore.concept %in% ultram$tip.label | VARS$Siphonophore.concept=="Nanomia bijuga"),]
VARS_cast = acast(VARS_curated, Siphonophore.concept~Prey.taxonomy, fun.aggregate = length)

#Prune morphological matrix to species in diet
dprunedmatrix = sharedmeans[which(sharedmeans$Species%in%rownames(diet)),]
dprunedmatrix_logs = sharedmean_logs[which(sharedmean_logs$Species%in%rownames(diet)),]
#Prune tree to diet species
dprunedTree = drop.tip(ultram, which(!(ultram$tip.label %in% rownames(diet))))

## Retrieve prey selectivity ##
selectivity = read.csv("Purcell1981_selectivity.csv", header=T, sep=",")[,c(1,31:39)]
selectivity$Species <- as.character(selectivity$Species)
selectivity$Species[which(selectivity$Species == "Nanomia bijuga")] <- "Nanomia sp"
selectivity$Species[which(selectivity$Species == "Rhizophysa eyesenhardti")] <- "Rhizophysa eysenhardtii"
selectivity = aggregate(. ~ Species, data=selectivity, FUN=mean)
rownames(selectivity) = selectivity$Species
selectivity=selectivity[,-1]
# Prune morphological data to selectivity tree
Sprunedmatrix = sharedmeans[which(sharedmeans$Species%in%rownames(selectivity)),]
Sprunedmatrix_logs = sharedmean_logs[which(sharedmean_logs$Species%in%rownames(selectivity)),]
selectivity = selectivity[which(rownames(selectivity) %in% sharedmeans$Species),]
# Prune tree to selectivity species
Sprunedtree = drop.tip(ultram, which(!(ultram$tip.label %in% rownames(selectivity))))

#Retrieve copepod prey length from literature
preylength = read.csv("Purcell1984_preylength.csv", header = T, sep=',')[,c(2,6)]
names(preylength) = c("Species","Copepod prey length (mm)")
PL_pruned_matrix = sharedmean_logs[which(sharedmean_logs$Species %in% preylength$Species),] %>% data.frame(preylength[which(preylength$Species %in% sharedmean_logs$Species),2])
PL_pruned_SEs = sharedvar_logs[which(sharedvar_logs$Species %in% preylength$Species),] %>% data.frame(preylength[which(preylength$Species %in% sharedmean_logs$Species),2])
names(PL_pruned_matrix)[ncol(PL_pruned_matrix)] <- "Copepod prey length (mm)"
names(PL_pruned_SEs)[ncol(PL_pruned_matrix)] <- "Copepod prey length (mm)"

#PGLS of characters vs copepod prey length
pGLSp_preylength = as.data.frame(matrix(nrow=ncol(PL_pruned_matrix[,-1]), ncol=1))
colnames(pGLSp_preylength) = "Copepod prey length"
row.names(pGLSp_preylength) = names(PL_pruned_matrix[,-1])
pGLS_PL_sign = pGLSp_preylength
pGLS_sign = pGLSp_preylength
for(i in 2:(ncol(PL_pruned_matrix)-1)){
  CH_I=as.numeric(PL_pruned_matrix[,i])
  names(CH_I) = PL_pruned_matrix$Species
  CH_I = CH_I[!is.na(CH_I)]
  CH_J=PL_pruned_matrix[,ncol(PL_pruned_matrix)]
  names(CH_J) = PL_pruned_matrix$Species
  CH_I = CH_I[!is.na(CH_I)]
  CH_I = CH_I[which(names(CH_I) %in% names(CH_J))]
  CH_J = CH_J[which(names(CH_J) %in% names(CH_I))]
  SE_I = PL_pruned_SEs[which(PL_pruned_SEs$Species %in% names(CH_I)),i]
  names(SE_I) = PL_pruned_SEs[which(PL_pruned_SEs$Species %in% names(CH_I)),1]
  SE_I = SE_I[match(names(CH_I), names(SE_I))]
  pgls_tree = drop.tip(ultram, which(!(ultram$tip.label %in% names(CH_I))))
  pg_data = as.data.frame(cbind(CH_I,CH_J)) %>% .[which(!is.na(rowSums(.))),]
  pgls.SEy(CH_I~CH_J, data=pg_data, tree = pgls_tree, se = SE_I) -> TTABLE
  pgls.SEy(CH_I~sample(CH_J, length(CH_J)), data=pg_data, tree = pgls_tree, se = SE_I) -> TTABLE_null
  TTABLE <- gls(CH_I ~ CH_J, correlation =  corBrownian(phy = pgls_tree), data = pg_data, method = "ML") %>% summary() %>% .$tTable %>% as.data.frame()
  TTABLE$`p-value`[2] -> pGLSp_preylength[i-1,1]
  if(pGLSp_preylength[i-1,1]<0.05){
    print(paste(rownames(pGLSp_preylength)[i-1],"Copepod prey length",pGLSp_preylength[(i-1),1],sep=" "))
    ValueIJ = TTABLE %>% .$Value %>% .[2]
    if(ValueIJ>0){pGLS_PL_sign[i-1,1] = "+"}
    else pGLS_PL_sign[i,1] = "-"
  }
}
write.csv(pGLSp_preylength, "PGLS_Pval_preylength.csv")
pGLS_PL_sign <- pGLS_PL_sign[which(rowSums(is.na(pGLS_PL_sign))!=ncol(pGLS_PL_sign)),which(colSums(is.na(pGLS_PL_sign))!=nrow(pGLS_PL_sign))]
write.csv(pGLS_PL_sign, "PGLS_SIGN_preylength.csv")

## Soft bodied / Hard bodied prey ##
#Soft vs Hard
soft_hard = cbind((diet$`Copepod diet`+diet$`Crustacean diet` + diet$`Decapod diet`+ diet$`Amphipod diet` + diet$`Ostracod diet`), (diet$`Fish diet`+diet$`Chaetognath diet`+diet$`Gelatinous diet`+diet$`Mollusc diet`+diet$`Polychaete diet`))
colnames(soft_hard) = c("Hard", "Soft")
rownames(soft_hard) = rownames(diet)
soft_hard <- as.data.frame(soft_hard)
soft_hard[soft_hard>0]<-1
softORhard = as.data.frame(soft_hard$Hard+soft_hard$Soft)
rownames(softORhard)=rownames(soft_hard)
names(softORhard) = "Type"
softORhard[softORhard==2]<-"Both"
softORhard[which(soft_hard$Hard==1 & soft_hard$Soft==0),]<-"Hard"
softORhard[which(soft_hard$Hard==0 & soft_hard$Soft==1),]<-"Soft"

##Test literature hypotheses ##

#Nematocyst type number vs soft/hard specialization test
nematocyst_number=as.numeric(categorical$Heteroneme.type!="")[which(!is.na(categorical$Desmonemes))]+as.numeric(categorical$Haploneme.type!="")[which(!is.na(categorical$Desmonemes))]+as.numeric(categorical$Rhopalonemes!="0")[which(!is.na(categorical$Desmonemes))]+as.numeric(categorical$Desmonemes!="0")[which(!is.na(categorical$Desmonemes))]
names(nematocyst_number)=categorical$Species[which(!is.na(categorical$Desmonemes))]
CH_J=softORhard[,1]
names(CH_J) = rownames(softORhard)
CH_J[CH_J!="Soft"]<-0
CH_J[CH_J=="Soft"]<-1
CH_I = nematocyst_number
CH_J = CH_J[which(names(CH_J)%in%names(CH_I))]
CH_I = CH_I[which(names(CH_I)%in%names(CH_J))]
treeIJ = drop.tip(ultram, which(!(ultram$tip.label %in% names(CH_J))))
CH_I = CH_I[match(treeIJ$tip.label, names(CH_I))]
CH_J = CH_J[match(treeIJ$tip.label, names(CH_J))]
datIJ = data.frame(preytype = CH_J, character = CH_I)
datIJ$preytype = as.numeric(as.character(datIJ$preytype))
logit = glm(preytype~character,data=datIJ, family="binomial")
print(logit)
phylogit = phyloglm(preytype~character,phy=treeIJ,data=datIJ,boot=10, btol=10)
print(phylogit)
c("Soft specilization", "Number of nematocyst types", phylogit$n, phylogit$aic, coef(summary(phylogit))[2,6], coef(summary(phylogit))[2,1], logit$aic, coef(summary(logit))[2,4], coef(summary(logit))[2,1])

#Desmoneme/Rhopaloneme vs crustacean
dprunedCat = binary[which(rownames(binary)%in% rownames(diet)),]
CH_I = dprunedCat$Rhopalonemes #same result with desmonemes, identical phylogenetic distribution
names(CH_I) = rownames(dprunedCat)
CH_J = diet$`Crustacean diet`
names(CH_J) = rownames(diet)
CH_J = CH_J[which(names(CH_J)%in%names(CH_I))]
CH_I = CH_I[which(names(CH_I)%in%names(CH_J))]
treeIJ = drop.tip(ultram, which(!(ultram$tip.label %in% names(CH_J))))
CH_I = CH_I[match(treeIJ$tip.label, names(CH_I))]
CH_J = CH_J[match(treeIJ$tip.label, names(CH_J))]
fitPagel(treeIJ, CH_I, CH_J)

#Differentiated cnidobands with heteronemes
CH_I = dprunedCat$Proximal.heteronemes
names(CH_I) = rownames(dprunedCat)
CH_J = diet$`Crustacean diet`
names(CH_J) = rownames(diet)
CH_J = CH_J[which(names(CH_J)%in%names(CH_I))]
CH_I = CH_I[which(names(CH_I)%in%names(CH_J))]
treeIJ = drop.tip(ultram, which(!(ultram$tip.label %in% names(CH_J))))
CH_I = CH_I[match(treeIJ$tip.label, names(CH_I))]
CH_J = CH_J[match(treeIJ$tip.label, names(CH_J))]
fitPagel(treeIJ, CH_I, CH_J) %>% .$P

## OUwie regime tree of feeding guilds ##

regimetree = dprunedTree
hypdiet = c("Small crustacean", "Small crustacean", "Small crustacean", "Small crustacean", "Small crustacean", "Large crustacean", "Mixed", "Mixed", "Mixed", "Large crustacean", "Large crustacean", "Mixed", "Small crustacean", "Large crustacean", "Fish", "Fish", "Fish", "Large crustacean", "Gelatinous", "Fish", "Fish", "Fish")
names(hypdiet) = regimetree$tip.label #Important thing here is that Forskalia tholoides was used as a proxy for Forskalia sp. in Purcell 1981b
DietAnc = ace(hypdiet, regimetree, type = "discrete")$lik.anc
HypDietAnc = 1:nrow(DietAnc)
for(row in 1:nrow(DietAnc)){
  SCR = scale(DietAnc[row,])
  print(SCR)
  HypDietAnc[row] <- rownames(SCR)[which(SCR[,1] == max(SCR[,1]))] %>% print()
}
regimetree$node.label = as.factor(HypDietAnc)
plotTree(regimetree)
nodelabels(text=regimetree$node.label,frame="none",adj=c(1.6,-0.45), cex=0.6);tiplabels(text=hypdiet, frame="none", cex=0.6, adj=c(0,2))
RTorder = regimetree$tip.label

##SIMMAP feeding guilds ##
set.seed(1111111)
make.simmap(regimetree, hypdiet, nsim = 100) -> feeding_sim
plotTree(regimetree, lwd = 4)
feeding_sim %>% plotSimmap(lwd = 4, add = T)
colors = c( "black","green3","blue","cyan","red")
names(colors) = c("Fish","Large crustacean","Mixed","Small crustacean","Gelatinous")
nodelabels(pie=(describe.simmap(feeding_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
par(ask=F)

#PGLS of characters vs Purcell Selectivity
pGLSp_sel = as.data.frame(matrix(nrow=ncol(Sprunedmatrix_logs[,-1]), ncol=ncol(selectivity)))
pGLS__sel_sign = pGLSp_sel
colnames(pGLSp_sel) = names(selectivity)
row.names(pGLSp_sel) = names(Sprunedmatrix_logs[,-1])
pGLS_sel_sign = pGLSp_sel
for(i in 2:ncol(Sprunedmatrix_logs)){
  CH_I=as.numeric(Sprunedmatrix_logs[,i])
  names(CH_I) = Sprunedmatrix_logs$Species
  CH_I = CH_I[!is.na(CH_I)]
  for(j in 1:ncol(selectivity)){
    CH_J=as.numeric(selectivity[,j])
    names(CH_J) = rownames(selectivity)
    CH_J = CH_J[!is.na(CH_J)]
    CH_I = CH_I[which(names(CH_I) %in% names(CH_J))]
    CH_J = CH_J[which(names(CH_J) %in% names(CH_I))]
    SE_I = castSE_logs[which(castSE$Species %in% names(CH_I)),i]
    names(SE_I) = castSE_logs[which(castSE$Species %in% names(CH_I)),1]
    SE_I = SE_I[match(names(CH_I), names(SE_I))]
    pgls_tree = drop.tip(Sprunedtree, which(!(Sprunedtree$tip.label %in% names(CH_I))))
    pg_data = as.data.frame(cbind(CH_I,CH_J)) %>% .[which(!is.na(rowSums(.))),]
    if(length(unique(CH_J))>1){
      pgls.SEy(CH_I~CH_J, data=pg_data, tree = pgls_tree, se = SE_I) -> TTABLE
      pgls.SEy(CH_I~sample(CH_J, length(CH_J)), data=pg_data, tree = pgls_tree, se = SE_I) -> TTABLE_null
      LikRat <- 2*(TTABLE$logLik[1] - TTABLE_null$logLik[1])
      pchisq(LikRat, df=1, lower.tail = F) -> pGLSp_sel[i-1,j]
      if(pGLSp_sel[i-1,j]<0.05){
        print(paste(rownames(pGLSp_sel)[i-1],colnames(pGLSp_sel)[j],pGLSp_sel[(i-1),(j)],TTABLE$coefficients[2], sep=" "))
        ValueIJ = TTABLE %>% .$coefficients %>% .[2]
        if(ValueIJ>0){pGLS_sel_sign[i-1,j] = "+"}
        else pGLS_sel_sign[i,j] = "-"
      }
    }
    else TTABLE<-NA; TTABLE_null<-NA; LikRat<-NA; pGLSp_sel[i-1,j]<-NA
  }
}
write.csv(pGLSp_sel, "PGLS_Pval_Selectivity.csv")
pGLS_sel_sign <- pGLS_sel_sign[which(rowSums(is.na(pGLSp_sel))!=ncol(pGLS_sel_sign)),which(colSums(is.na(pGLS_sel_sign))!=nrow(pGLS_sel_sign))]
write.csv(pGLS_sel_sign, "PGLS_SIGN_Selectivity.csv")

## Logistic regressions ##
logits = data.frame(matrix(ncol = 9, nrow = 0))
names(logits) = c("Prey type", "Character", "Ntaxa", "phyloglmAIC", "phyloglmP","phyloglm_b","glmAIC", "glmP", "glm_b")
logit_diet = diet[,-7]
logits[1,] = rep(1,9)
for(prey in 1:ncol(logit_diet)){
  for(char in 2:ncol(dprunedmatrix_logs)){
    CH_J=logit_diet[,prey]
    names(CH_J) = rownames(logit_diet)
    CH_I = dprunedmatrix_logs[,char]
    names(CH_I) = dprunedmatrix_logs$Species
    CH_I = CH_I[!is.na(CH_I)]
    CH_J = CH_J[which(names(CH_J)%in%names(CH_I))]
    treeIJ = drop.tip(dprunedTree, which(!(dprunedTree$tip.label %in% names(CH_J))))
    CH_I = CH_I[match(treeIJ$tip.label, names(CH_I))]
    CH_J = CH_J[match(treeIJ$tip.label, names(CH_J))]
    datIJ = data.frame(preytype = CH_J, character = CH_I)
    if(length(unique(CH_I))>1 & length(unique(CH_J))>1){
      logit = glm(preytype~character,data=datIJ, family="binomial")
      print(logit)
      phylogit = phyloglm(preytype~character,phy=treeIJ,data=datIJ,boot=10, btol=30)
      print(phylogit)
      rowIJ = c(names(logit_diet)[prey], names(dprunedmatrix_logs)[char], phylogit$n, phylogit$aic, coef(summary(phylogit))[2,6], coef(summary(phylogit))[2,1], logit$aic, coef(summary(logit))[2,4], coef(summary(logit))[2,1])
    } else{
      rowIJ = c(names(logit_diet)[prey], names(dprunedmatrix_logs)[char], NA, NA, NA, NA, NA, NA, NA)
    }
    print(rowIJ)
    logits = rbind(logits, rowIJ)
  }
}
logits = logits[-1,]
logits[,4:9] = sapply(logits[,4:9], as.numeric)
logits[,4:9] = sapply(logits[,4:9], function(x){round(x,3)})
write.csv(logits, "diet_morph_logits.csv", row.names = F)

## OUWIE ##
OUwie_matrix = as.data.frame(matrix(ncol=5, nrow=nrow(dprunedmatrix_logs)))
names(OUwie_matrix) = c("Character", "Ntips", "dAICc_BM", "dAICc_OU1", "dAICc_OUm")
library(magrittr)
for(c in 2:ncol(dprunedmatrix_logs)){
  Cdata = data.frame(dprunedmatrix$Species,hypdiet,dprunedmatrix_logs[,c])
  names(Cdata) = c("species", "regime", "trait")
  Cdata$trait[which(Cdata$trait == "NaN")] <- NA
  Cdata = Cdata[!is.na(Cdata$trait),]
  Cdata$trait = as.numeric(Cdata$trait)
  Cdata$species = as.character(Cdata$species)
  Cdata$regime = as.character(Cdata$regime)
  Ctree = drop.tip(regimetree,which(!(regimetree$tip.label %in% Cdata$species)))
  BM <- OUwie(Ctree, Cdata, model="BM1")
  OU1 <- OUwie(Ctree, Cdata, model="OU1")
  OUm <- OUwie(Ctree, Cdata, model="OUMV")
  OUwie_matrix[c-1,1] = names(dprunedmatrix_logs)[c]
  OUwie_matrix[c-1,2] = nrow(Cdata)
  OUwie_matrix[c-1,3] = BM$AICc - min(c(BM$AICc, OU1$AICc, OUm$AICc))
  OUwie_matrix[c-1,4] = OU1$AICc - min(c(BM$AICc, OU1$AICc, OUm$AICc))
  OUwie_matrix[c-1,5] = OUm$AICc - min(c(BM$AICc, OU1$AICc, OUm$AICc))
}
OUwie_matrix[,3:5] <- apply(OUwie_matrix[,3:5], 2, round3)
write.csv(OUwie_matrix, "OU_diet_test.csv")

#OUwie model adequacy
#Model adequacy
worthy_models = OUwie_matrix
MAD = as.data.frame(matrix(ncol = 6, nrow = nrow(worthy_models)))
names(MAD) = c("msig", "cvar", "svar", "sasr", "shgt", "dcfd")
rownames(MAD) = worthy_models$Character
for(m in 1:nrow(worthy_models)){
  C = dprunedmatrix_logs[,which(names(dprunedmatrix_logs) == worthy_models$Character[m])]
  names(C) = dprunedmatrix_logs$Species
  C = C[!is.na(C)]
  Ctree = drop.tip(regimetree, which(!(regimetree$tip.label %in% names(C))))
  C = C[match(Ctree$tip.label, names(C))]
  class(Ctree)="phylo"
  FC <- fitContinuous(Ctree,C,model=worthy_models$Best_model[m])
  UTC <- make_unit_tree(FC)
  picstat_data <- calculate_pic_stat(UTC)
  sim <- simulate_char_unit(UTC)
  picstat_sim <- calculate_pic_stat(sim)
  compare_pic_stat(picstat_data, picstat_sim) %>% .$p.values -> MAD[m,]
}
phy_models = data.frame(worthy_models,MAD)
phy_models[,6:11] = apply(phy_models[,6:11], 2,round3)
#write.csv(phy_models, "OUwie_model_adequacy.csv")

## DAPC ##
#General prediction of feeding guild (transforming NAs to 0)
ldamtrix = sharedmean_logs
ldamtrix[is.na(ldamtrix)] <- 0
ldamtrix$Species = as.character(ldamtrix$Species)
ldamtrix = ldamtrix[which(ldamtrix$Species %in% names(hypdiet)),]
prunediets = hypdiet[which(names(hypdiet) %in% ldamtrix$Species)]
#prunediets["Hippopodius hippopus"]<-"Ostracod"
ldamtrix = cbind(ldamtrix, prunediets[match(ldamtrix$Species, names(prunediets))])
names(ldamtrix)[ncol(ldamtrix)]<-"Diet"
dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$Diet, n.pca = 29, n.da = 5) -> DAPClogs_test
DAPClogs_guilds <- dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$Diet, n.pca = optim.a.score(DAPClogs_test, n=100)$best, n.da = 5)
scatter(DAPClogs_guilds, clabel = 0.5, col=c("blue","dark green","orange","purple","red"))
contributions_guilds = rowSums(DAPClogs_guilds$var.contr %*% diag(DAPClogs_guilds$eig))
contributions_guilds = contributions_guilds[order(contributions_guilds,decreasing = T)]
top_guilds = contributions_guilds[which(contributions_guilds>summary(contributions_guilds)[5])]
top_guilds <-  as.matrix(top_guilds)
colnames(top_guilds)<-"Variable contribution"
summary(DAPClogs_guilds)$assign.prop

#Particular diet items that have enough variability (Copepods, fish, large crustaceans)
#Copepods
ldamtrix = sharedmean_logs
ldamtrix[is.na(ldamtrix)] <- 0
ldamtrix$Species = as.character(ldamtrix$Species)
ldamtrix = ldamtrix[which(ldamtrix$Species %in% rownames(diet)),]
prunediets = diet[which(rownames(diet) %in% ldamtrix$Species),]
ldamtrix = cbind(ldamtrix, prunediets[match(ldamtrix$Species, rownames(prunediets)), "Copepod diet"])
names(ldamtrix)[ncol(ldamtrix)]<-"CopepodDiet"
dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$CopepodDiet, n.pca = 29, n.da = 2) -> DAPClogs_test
dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$CopepodDiet, n.pca = optim.a.score(DAPClogs_test, n=100)$best, n.da = 2) -> DAPClogs_Cope
summary.dapc(DAPClogs_Cope)$assign.prop
scatter(DAPClogs_Cope, clabel = 0.5)
contributions_cope = rowSums(DAPClogs_Cope$var.contr %*% DAPClogs_Cope$eig)
contributions_cope = contributions_cope[order(contributions_cope,decreasing = T)]
top_cope = contributions_cope[which(contributions_cope>summary(contributions_cope)[5])]
top_cope <-  as.matrix(top_cope)
colnames(top_cope)<-"Variable contribution"

#Copepod GLM
top_cope %>% rownames() %>% paste(collapse = " + ")
glmatrix_cope = ldamtrix
glmatrix_cope[is.na(glmatrix_cope)]<-0
CopepodGLM <- glm(CopepodDiet~cnidomic_index + Tentacle.width..um. + haploneme_elongation + SAV_haploneme + Haploneme.row.number..um. + Cnidoband.length..um. + Cnidoband.width..um. + Cnidoband.free.length..um., data = glmatrix_cope)
summary(CopepodGLM)
1-(CopepodGLM$deviance/CopepodGLM$null.deviance)

#Copepod selectivity GLM
copselmatrix = cbind(Sprunedmatrix, selectivity$Selectivity.Copepods)
copselmatrix = copselmatrix[which(copselmatrix$Species %in% ldamtrix$Species),]
names(copselmatrix)[ncol(copselmatrix)]<-"CopepodSelectivity"
copselmatrix[is.na(copselmatrix)]<-0
CopepodSelGLM <- glm(CopepodSelectivity~cnidomic_index + Tentacle.width..um. + haploneme_elongation + SAV_haploneme + Haploneme.row.number..um. + Cnidoband.length..um. + Cnidoband.width..um. + Cnidoband.free.length..um., data = copselmatrix)
summary(CopepodSelGLM)
1-(CopepodSelGLM$deviance/CopepodSelGLM$null.deviance)

#Fish
ldamtrix = sharedmean_logs
ldamtrix[is.na(ldamtrix)] <- 0
ldamtrix$Species = as.character(ldamtrix$Species)
ldamtrix = ldamtrix[which(ldamtrix$Species %in% rownames(diet)),]
prunediets = diet[which(rownames(diet) %in% ldamtrix$Species),]
ldamtrix = cbind(ldamtrix, prunediets[match(ldamtrix$Species, rownames(prunediets)), "Fish diet"])
names(ldamtrix)[ncol(ldamtrix)]<-"FishDiet"
dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$FishDiet, n.pca = 29, n.da = 2) -> DAPClogs_test
dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$FishDiet, n.pca = optim.a.score(DAPClogs_test, n=100)$best, n.da = 2) -> DAPClogs_Fish
summary.dapc(DAPClogs_Fish)$assign.prop
scatter(DAPClogs_Fish, clabel = 0.5)
contributions_fish = rowSums(DAPClogs_Fish$var.contr %*% DAPClogs_Fish$eig)
contributions_fish = contributions_fish[order(contributions_fish,decreasing = T)]
top_fish = contributions_fish[which(contributions_fish>summary(contributions_fish)[5])]
top_fish <-  as.matrix(top_fish)
colnames(top_fish)<-"Variable contribution"

#Fish GLM
top_fish %>% rownames() %>% paste(collapse = " + ")
glmatrix_fish = ldamtrix
glmatrix_fish[is.na(glmatrix_fish)]<-0
FishGLM <- glm(FishDiet~total_haploneme_volume + Heteroneme.volume..um3. + cnidomic_index + total_heteroneme_volume + Cnidoband.length..um. + Cnidoband.free.length..um. + Involucrum.length..um. + Pedicle.width..um., data = glmatrix_fish)
summary(FishGLM)
1-(FishGLM$deviance/FishGLM$null.deviance)

#Fish selectivity GLM
fishselmatrix = cbind(Sprunedmatrix, selectivity$Selectivity.Fish)
fishselmatrix = fishselmatrix[which(fishselmatrix$Species %in% ldamtrix$Species),]
names(fishselmatrix)[ncol(fishselmatrix)]<-"FishSelectivity"
fishselmatrix[is.na(fishselmatrix)]<-0
FishSelGLM <- glm(FishSelectivity~total_haploneme_volume + Heteroneme.volume..um3. + cnidomic_index + total_heteroneme_volume + Cnidoband.length..um. + Cnidoband.free.length..um. + Involucrum.length..um. + Pedicle.width..um., data = fishselmatrix)
summary(FishSelGLM)
1-(FishSelGLM$deviance/FishSelGLM$null.deviance) 

#Large crustaceans
ldamtrix = sharedmean_logs
ldamtrix[is.na(ldamtrix)] <- 0
ldamtrix$Species = as.character(ldamtrix$Species)
ldamtrix = ldamtrix[which(ldamtrix$Species %in% rownames(diet)),]
prunediets = diet[which(rownames(diet) %in% ldamtrix$Species),]
ldamtrix = cbind(ldamtrix, prunediets[match(ldamtrix$Species, rownames(prunediets)), "Decapod diet"])
names(ldamtrix)[ncol(ldamtrix)]<-"ShrimpDiet"
dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$ShrimpDiet, n.pca = 29, n.da = 2) -> DAPClogs_test
dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$ShrimpDiet, n.pca = optim.a.score(DAPClogs_test, n=100)$best, n.da = 2) -> DAPClogs_Shrimp
summary.dapc(DAPClogs_Shrimp)$assign.prop
scatter(DAPClogs_Shrimp, clabel = 0.5)
contributions_shrimp = rowSums(DAPClogs_Shrimp$var.contr %*% DAPClogs_Shrimp$eig)
contributions_shrimp = contributions_shrimp[order(contributions_shrimp,decreasing = T)]
top_shrimp = contributions_shrimp[which(contributions_shrimp>summary(contributions_shrimp)[5])]
top_shrimp <-  as.matrix(top_shrimp)
colnames(top_shrimp)<-"Variable contribution"

#Shrimp GLM
top_shrimp %>% rownames() %>% paste(collapse = " + ")
glmatrix_shrimp = ldamtrix
glmatrix_shrimp[is.na(glmatrix_shrimp)]<-0
ShrimpGLM <- glm(ShrimpDiet~Involucrum.length..um. + total_heteroneme_volume + Elastic.strand.width..um. + Rhopaloneme.length..um. + Heteroneme.volume..um3. + haploneme_elongation + Desmoneme.length..um. + Tentacle.width..um., data = glmatrix_shrimp)
summary(ShrimpGLM)
1-(ShrimpGLM$deviance/ShrimpGLM$null.deviance)

#Large crustacean selectivity GLM
decselmatrix = cbind(Sprunedmatrix, (selectivity$Selectivity.Decapod.larvae+selectivity$Selectivity.Shrimp)/2)
decselmatrix = decselmatrix[which(decselmatrix$Species %in% ldamtrix$Species),]
names(decselmatrix)[ncol(decselmatrix)]<-"ShrimpSelectivity"
decselmatrix[is.na(decselmatrix)]<-0
DecapodSelGLM <- glm(ShrimpSelectivity~Involucrum.length..um. + total_heteroneme_volume + Elastic.strand.width..um. + Rhopaloneme.length..um. + Heteroneme.volume..um3. + haploneme_elongation + Desmoneme.length..um. + Tentacle.width..um., data = decselmatrix)
summary(DecapodSelGLM)
1-(DecapodSelGLM$deviance/DecapodSelGLM$null.deviance)

#Soft bodied vs hard bodied prey
prunediets = softORhard
ldamtrix = sharedmean_logs
ldamtrix[is.na(ldamtrix)]<-0
ldamtrix$Species = as.character(ldamtrix$Species)
ldamtrix = ldamtrix[which(ldamtrix$Species %in% rownames(prunediets)),]
ldamtrix = cbind(ldamtrix, prunediets[match(ldamtrix$Species, rownames(prunediets)), 1])
names(ldamtrix)[ncol(ldamtrix)]<-"SoftHard"
rownames(ldamtrix) = ldamtrix$Species
dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$SoftHard, n.pca = 29, n.da = 3) -> DAPClogs_test
dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$SoftHard, n.pca = optim.a.score(DAPClogs_test, n=100)$best, n.da = 2) -> DAPClogs_SoftHard
summary.dapc(DAPClogs_SoftHard)$assign.prop
scatter(DAPClogs_SoftHard, clabel = 0.5, col = c("purple", "red", "blue"))
contributions_softHard = rowSums(DAPClogs_SoftHard$var.contr %*% DAPClogs_SoftHard$eig)
contributions_softHard = contributions_softHard[order(contributions_softHard,decreasing = T)]
top_softHard = contributions_softHard[which(contributions_softHard>summary(contributions_softHard)[5])]
top_softHard <-  as.matrix(top_softHard)
colnames(top_softHard)<-"Variable contribution"
predictionset <- castmean_logs[,which(names(castmean_logs) %in% names(ldamtrix))]
predictionset[is.na(predictionset)]<-0
rownames(predictionset)<-predictionset$Species
predictionset=predictionset[which(!(predictionset$Species %in% c("Thermopalia taraxaca", "Nectadamas richardi", "Halistemma cupulifera", "Cardianecta parchelion"))),] #Remove taxa with many unmeasured NAs that create zero-biases
predictionset=predictionset[which(!(predictionset$Species %in% ldamtrix$Species)),]
preDIET <- predict(DAPClogs_SoftHard, predictionset[,-1])
preDIET$posterior %>% round(5) %>%  as.matrix() %>% heatmap(scale="col", cexCol=0.5, col=c("white", "white", gray.colors(10)[10:1],"black", "black"),Colv=NA)

##VCV analysis###

#For the whole taxa set
regimematrix = sharedmean_logs
rownames(regimematrix) = regimematrix$Species
#NAs to Zeroes
regimematrix[is.na(regimematrix)] <- 0
simple.vcv <- evol.vcv(tree=ultram, as.matrix(regimematrix[,2:22]))
simple.vcv.corrs <- cov2cor(simple.vcv$R.single)
corrplot(simple.vcv.corrs)
#NA taxa removed
regimematrix = sharedmean_logs
rownames(regimematrix) = regimematrix$Species
prunedmatrix <- as.matrix(regimematrix[which(!is.na(rowSums(regimematrix[,-1]))),c(2,3,5:22)])
prunedphylo <- drop.tip(ultram, which(!(ultram$tip.label %in% rownames(prunedmatrix))))
simple.vcv <- evol.vcv(tree=prunedphylo, prunedmatrix)
simple.vcv.corrs <- cov2cor(simple.vcv$R.single)
corrplot(simple.vcv.corrs)

#For the taxa with diet data
par(ask=F)
dietregimematrix = dprunedmatrix_logs
rownames(dietregimematrix) = dietregimematrix$Species
#Nas to zeroes
#dietregimematrix[is.na(dietregimematrix)] <- 0
diet.simple.vcv <- evol.vcv(tree=regimetree, as.matrix(dietregimematrix[,2:22]))
diet.simple.vcv.corrs <- cov2cor(diet.simple.vcv$R.single)
corrplot(diet.simple.vcv.corrs)
#NA taxa removed   ----  TOO MANY PAIRS ARE COMPUTATIONALLY SINGULAR -----
# dietregimematrix = dprunedmatrix_logs
# rownames(dietregimematrix) = dietregimematrix$Species
# pruneddietmatrix <- as.matrix(dietregimematrix[which(!is.na(rowSums(dietregimematrix[,-1]))),c(2,3,5,6,7,8,9,10,11,12)])
# pruneddietphylo <- drop.tip(regimetree, which(!(regimetree$tip.label %in% rownames(pruneddietmatrix))))
# diet.simple.vcv <- evol.vcv(tree=pruneddietphylo, as.matrix(pruneddietmatrix))
# diet.simple.vcv.corrs <- cov2cor(diet.simple.vcv$R.single)
# corrplot(diet.simple.vcv.corrs)

##VCVacross the Regimes in the tree
#Loop for evol.vcvlite with character pairs
VCVlist <- list()
it = 0
for(i in 2:ncol(dietregimematrix)){
  for(j in 2:ncol(dietregimematrix)){
    IJmatrix = dietregimematrix[which(!is.na(rowSums(dietregimematrix[,c(i,j)]))),]
    IJtree = drop.tip.simmap(feeding_sim[[19]], feeding_sim[[19]]$tip.label[which(!(feeding_sim[[19]]$tip.label %in% IJmatrix$Species))])
    #remaining_states = IJtree$node.states %>% as.vector() %>%  unique()
    #Qmod = IJtree$Q
    #Qmod = Qmod[which(rownames(Qmod) %in% remaining_states),which(colnames(Qmod) %in% remaining_states)]
    #IJtree$Q <- Qmod
    IJtree$node.label=NULL
    IJmatrix = IJmatrix[match(IJtree$tip.label, rownames(IJmatrix)),]
    it <- it+1
    if(i != j ){
      evcvIJ <- try(evolvcv.lite(IJtree, as.matrix(IJmatrix[,c(i,j)])), silent=T)
      VCVlist[[it]] <- evcvIJ
      names(VCVlist)[it] <- paste(names(IJmatrix)[i], names(IJmatrix)[j], collapse = " ~ ")
    }
  }
}

ntaxamatrix = matrix(ncol=ncol(dietregimematrix)-1,nrow=ncol(dietregimematrix)-1) %>% as.data.frame()
names(ntaxamatrix) <- names(dietregimematrix)[-1]
rownames(ntaxamatrix) <- names(dietregimematrix)[-1]

#Get the number of taxa per character pair
for(i in 2:ncol(dietregimematrix)){
  for(j in 2:ncol(dietregimematrix)){
    which(!is.na(rowSums(dietregimematrix[,c(i,j)]))) %>% length() -> ntaxamatrix[i-1,j-1] #%>% print()
  }
}
ntaxamelt <- cbind(rep(rownames(ntaxamatrix),ncol(ntaxamatrix)), melt(ntaxamatrix, id.var = NULL))
names(ntaxamelt) <- c("Character_1","Character_2","Ntaxa")
ggplot(ntaxamelt, aes(Character_1, Character_2)) + geom_tile(aes(fill = Ntaxa), colour = "black") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


VCVshortlist <- VCVlist[which(lapply(VCVlist, class)=="evolvcv.lite")]
bestmodel_number = lapply(VCVshortlist, function(y){lapply(y,function(x){x$AIC}) %>% unlist() %>% .[which(. == min(.))] %>% names()})
bestmodel_number <- bestmodel_number %>% unlist() %>% as.vector() %>% str_extract("[0-9]") %>% as.numeric()
bestmodel_descriptors <- list()
names(bestmodel_descriptors)=names(VCVshortlist)
for(i in 1:length(bestmodel_number)){
  VCVshortlist[[i]][bestmodel_number[i]][[1]] -> bestmodel_descriptors[[i]]
}
names(bestmodel_descriptors)=names(VCVshortlist)

##R.simple VCV##
generalmodel_descriptors <- list()
names(generalmodel_descriptors)=names(VCVshortlist)
for(i in 1:length(VCVshortlist)){
  VCVshortlist[[i]][1][[1]] -> generalmodel_descriptors[[i]]
}
names(generalmodel_descriptors)=names(VCVshortlist)
generalmodel_Simple <- generalmodel_descriptors %>% lapply(function(x){x$R %>% unlist() %>% matrix(nrow=2,ncol=2) %>% cov2cor() %>% as.data.frame()})
for(i in 1:length(generalmodel_Simple)){
  rownames(generalmodel_Simple[[i]])<- names(generalmodel_Simple)[i] %>% str_split(" ") %>% unlist()
  names(generalmodel_Simple[[i]])<- names(generalmodel_Simple)[i] %>% str_split(" ") %>% unlist()
  generalmodel_Simple[[i]] %<>% mutate(parameter=rownames(.))
}
generalmodel_Simple %>% reduce(full_join) -> Simplejoin
Simplejoin <- aggregate(. ~ Simplejoin$parameter, data = Simplejoin[,-3], mean.na, na.action = na.pass)
rownames(Simplejoin)<-Simplejoin$`Simplejoin$parameter`
Simplejoin <- Simplejoin[,-1]
Simplejoin = Simplejoin[match(sort(rownames(Simplejoin)), rownames(Simplejoin)), match(sort(names(Simplejoin)), names(Simplejoin))]
#Simplejoin[is.na(Simplejoin)] <- 0
corrplot(Simplejoin %>% as.matrix())

##FISH SPECIFIC VCV##
bestmodel_Fish <- bestmodel_descriptors %>% lapply(function(y){if(is.matrix(y$R["Fish"][[1]])){y$R["Fish"] %>% unlist() %>% matrix(nrow=2,ncol=2) %>% cov2cor() %>% as.data.frame()}})
bestmodel_Fish <- bestmodel_Fish[which(sapply(bestmodel_Fish, function(e) is.list(e)))]
for(i in 1:length(bestmodel_Fish)){
  rownames(bestmodel_Fish[[i]])<- names(bestmodel_Fish)[i] %>% str_split(" ") %>% unlist()
  names(bestmodel_Fish[[i]])<- names(bestmodel_Fish)[i] %>% str_split(" ") %>% unlist()
  bestmodel_Fish[[i]] %<>% mutate(parameter=rownames(.))
}
bestmodel_Fish %>% reduce(full_join) -> Fishjoin
Fishjoin <- aggregate(. ~ Fishjoin$parameter, data = Fishjoin[,-3], mean.na, na.action = na.pass)
rownames(Fishjoin)<-Fishjoin$`Fishjoin$parameter`
Fishjoin <- Fishjoin[,-1]
Fishjoin = Fishjoin[match(sort(rownames(Fishjoin)), rownames(Fishjoin)), match(sort(names(Fishjoin)), names(Fishjoin))]
#Fishjoin[is.na(Fishjoin)] <- 0
corrplot(Fishjoin %>% as.matrix())
#difference with whole diet tree simple VCV
subtractFish <- Fishjoin - Simplejoin[which(rownames(Simplejoin) %in% rownames(Fishjoin)),which(names(Simplejoin) %in% names(Fishjoin))]
subtractFish[is.na(subtractFish)] <- 0
(subtractFish/max(abs(subtractFish),na.rm = T)) %>% as.matrix() %>% corrplot()
#difference with preceding regime VCV
subtractFish <- Fishjoin - LargeCrustaceanjoin[which(rownames(Simplejoin) %in% rownames(Fishjoin)),which(names(Simplejoin) %in% names(Fishjoin))]
subtractFish[is.na(subtractFish)] <- 0
(subtractFish/max(abs(subtractFish),na.rm = T)) %>% as.matrix() %>% corrplot()

##GELATINOUS SPECIFIC VCV## not all items in list have this sub R matrix
bestmodel_Gelatinous <- bestmodel_descriptors %>% lapply(function(y){if(is.matrix(y$R["Gelatinous"][[1]])){y$R["Gelatinous"] %>% unlist() %>% matrix(nrow=2,ncol=2) %>% cov2cor() %>% as.data.frame()}})
bestmodel_Gelatinous <- bestmodel_Gelatinous[which(sapply(bestmodel_Gelatinous, function(e) is.list(e)))]
for(i in 1:length(bestmodel_Gelatinous)){
  rownames(bestmodel_Gelatinous[[i]])<- names(bestmodel_Gelatinous)[i] %>% str_split(" ") %>% unlist()
  names(bestmodel_Gelatinous[[i]])<- names(bestmodel_Gelatinous)[i] %>% str_split(" ") %>% unlist()
  bestmodel_Gelatinous[[i]] %<>% mutate(parameter=rownames(.))
}
bestmodel_Gelatinous %>% reduce(full_join) -> Gelatinousjoin
Gelatinousjoin <- aggregate(. ~ Gelatinousjoin$parameter, data = Gelatinousjoin[,-3], mean.na, na.action = na.pass)
rownames(Gelatinousjoin)<-Gelatinousjoin$`Gelatinousjoin$parameter`
Gelatinousjoin <- Gelatinousjoin[,-1]
Gelatinousjoin = Gelatinousjoin[match(sort(rownames(Gelatinousjoin)), rownames(Gelatinousjoin)), match(sort(names(Gelatinousjoin)), names(Gelatinousjoin))]
#Gelatinousjoin[is.na(Gelatinousjoin)] <- 0
corrplot(Gelatinousjoin %>% as.matrix())
#difference with whole diet tree simple VCV
subtractGelatinous <- Gelatinousjoin - Simplejoin[which(rownames(Simplejoin) %in% rownames(Gelatinousjoin)),which(names(Simplejoin) %in% names(Gelatinousjoin))]
subtractGelatinous[is.na(subtractGelatinous)] <- 0
(subtractGelatinous/max(abs(subtractGelatinous),na.rm = T)) %>% as.matrix() %>% corrplot()
#difference with preceding regime VCV
subtractGelatinous <- Gelatinousjoin - LargeCrustaceanjoin[which(rownames(Simplejoin) %in% rownames(Gelatinousjoin)),which(names(Simplejoin) %in% names(Gelatinousjoin))]
subtractGelatinous[is.na(subtractGelatinous)] <- 0
(subtractGelatinous/max(abs(subtractGelatinous),na.rm = T)) %>% as.matrix() %>% corrplot()

##LARGE CRUSTACEAN SPECIFIC VCV##
bestmodel_LargeCrustacean <- bestmodel_descriptors %>% lapply(function(x){x$R["Large crustacean"] %>% unlist() %>% matrix(nrow=2,ncol=2) %>% cov2cor() %>% as.data.frame()})
for(i in 1:length(bestmodel_LargeCrustacean)){
  rownames(bestmodel_LargeCrustacean[[i]])<- names(bestmodel_LargeCrustacean)[i] %>% str_split(" ") %>% unlist()
  names(bestmodel_LargeCrustacean[[i]])<- names(bestmodel_LargeCrustacean)[i] %>% str_split(" ") %>% unlist()
  bestmodel_LargeCrustacean[[i]] %<>% mutate(parameter=rownames(.))
}
bestmodel_LargeCrustacean %>% reduce(full_join) -> LargeCrustaceanjoin
LargeCrustaceanjoin <- aggregate(. ~ LargeCrustaceanjoin$parameter, data = LargeCrustaceanjoin[,-3], mean.na, na.action = na.pass)
rownames(LargeCrustaceanjoin)<-LargeCrustaceanjoin$`LargeCrustaceanjoin$parameter`
LargeCrustaceanjoin <- LargeCrustaceanjoin[,-1]
LargeCrustaceanjoin = LargeCrustaceanjoin[match(sort(rownames(LargeCrustaceanjoin)), rownames(LargeCrustaceanjoin)), match(sort(names(LargeCrustaceanjoin)), names(LargeCrustaceanjoin))]
#LargeCrustaceanjoin[is.na(LargeCrustaceanjoin)] <- 0
corrplot(LargeCrustaceanjoin %>% as.matrix())
#difference with whole diet tree simple VCV
subtractLargeCrustacean <- LargeCrustaceanjoin - Simplejoin
subtractLargeCrustacean[is.na(subtractLargeCrustacean)] <- 0
(subtractLargeCrustacean/max(abs(subtractLargeCrustacean),na.rm = T)) %>% as.matrix() %>% corrplot()

##GENERALIST SPECIFIC VCV##
bestmodel_Mixed <- bestmodel_descriptors %>% lapply(function(x){x$R["Mixed"] %>% unlist() %>% matrix(nrow=2,ncol=2) %>% cov2cor() %>% as.data.frame()})
for(i in 1:length(bestmodel_Mixed)){
  rownames(bestmodel_Mixed[[i]])<- names(bestmodel_Mixed)[i] %>% str_split(" ") %>% unlist()
  names(bestmodel_Mixed[[i]])<- names(bestmodel_Mixed)[i] %>% str_split(" ") %>% unlist()
  bestmodel_Mixed[[i]] %<>% mutate(parameter=rownames(.))
}
bestmodel_Mixed %>% reduce(full_join) -> Mixedjoin
Mixedjoin <- aggregate(. ~ Mixedjoin$parameter, data = Mixedjoin[,-3], mean.na, na.action = na.pass)
rownames(Mixedjoin)<-Mixedjoin$`Mixedjoin$parameter`
Mixedjoin <- Mixedjoin[,-1]
Mixedjoin = Mixedjoin[match(sort(rownames(Mixedjoin)), rownames(Mixedjoin)), match(sort(names(Mixedjoin)), names(Mixedjoin))]
#Mixedjoin[is.na(Mixedjoin)] <- 0
corrplot(Mixedjoin %>% as.matrix())
#difference with whole diet tree simple VCV
subtractMixed <- Mixedjoin - Simplejoin
subtractMixed[is.na(subtractMixed)] <- 0
(subtractMixed/max(abs(subtractMixed),na.rm = T)) %>% as.matrix() %>% corrplot()
subtractMixed <- Mixedjoin - LargeCrustaceanjoin
subtractMixed[is.na(subtractMixed)] <- 0
(subtractMixed/max(abs(subtractMixed),na.rm = T)) %>% as.matrix() %>% corrplot()

##SMALL CRUSTACEAN SPECIFIC VCV##
bestmodel_SmallCrustacean <- bestmodel_descriptors %>% lapply(function(y){if(is.matrix(y$R["Small crustacean"][[1]])){y$R["Small crustacean"] %>% unlist() %>% matrix(nrow=2,ncol=2) %>% cov2cor() %>% as.data.frame()}})
bestmodel_SmallCrustacean <- bestmodel_SmallCrustacean[which(sapply(bestmodel_SmallCrustacean, function(e) is.list(e)))]
for(i in 1:length(bestmodel_SmallCrustacean)){
  rownames(bestmodel_SmallCrustacean[[i]])<- names(bestmodel_SmallCrustacean)[i] %>% str_split(" ") %>% unlist()
  names(bestmodel_SmallCrustacean[[i]])<- names(bestmodel_SmallCrustacean)[i] %>% str_split(" ") %>% unlist()
  bestmodel_SmallCrustacean[[i]] %<>% mutate(parameter=rownames(.))
}
bestmodel_SmallCrustacean %>% reduce(full_join) -> SmallCrustaceanjoin
SmallCrustaceanjoin <- aggregate(. ~ SmallCrustaceanjoin$parameter, data = SmallCrustaceanjoin[,-3], mean.na, na.action = na.pass)
rownames(SmallCrustaceanjoin)<-SmallCrustaceanjoin$`SmallCrustaceanjoin$parameter`
SmallCrustaceanjoin <- SmallCrustaceanjoin[,-1]
SmallCrustaceanjoin = SmallCrustaceanjoin[match(sort(rownames(SmallCrustaceanjoin)), rownames(SmallCrustaceanjoin)), match(sort(names(SmallCrustaceanjoin)), names(SmallCrustaceanjoin))]
#SmallCrustaceanjoin[is.na(SmallCrustaceanjoin)] <- 0
corrplot(SmallCrustaceanjoin %>% as.matrix())
#difference with whole diet tree simple VCV
subtractSmallCrustacean <- SmallCrustaceanjoin - Simplejoin[which(rownames(Simplejoin) %in% rownames(SmallCrustaceanjoin)),which(names(Simplejoin) %in% names(SmallCrustaceanjoin))]
subtractSmallCrustacean[is.na(subtractSmallCrustacean)] <- 0
(subtractSmallCrustacean/max(abs(subtractSmallCrustacean),na.rm = T)) %>% as.matrix() %>% corrplot()
subtractSmallCrustacean <- SmallCrustaceanjoin - LargeCrustaceanjoin[which(rownames(LargeCrustaceanjoin) %in% rownames(SmallCrustaceanjoin)),which(names(LargeCrustaceanjoin) %in% names(SmallCrustaceanjoin))]
subtractSmallCrustacean[is.na(subtractSmallCrustacean)] <- 0
(subtractSmallCrustacean/max(abs(subtractSmallCrustacean),na.rm = T)) %>% as.matrix() %>% corrplot()

### Best model supported for each pair

bestmodel_types <- bestmodel_descriptors %>% lapply(function(x){x$description})
names(bestmodel_types)<-names(VCVshortlist)
VCVsolutions = data.frame(str_split_fixed(names(bestmodel_types), pattern=" ", n=Inf), unlist(bestmodel_types))
rownames(VCVsolutions)=1:nrow(VCVsolutions)
names(VCVsolutions) <- c("Character_1", "Character_2", "Best_model")
castSolutions = matrix(nrow=length(unique(VCVsolutions$Character_2)),ncol=length(unique(VCVsolutions$Character_2))) %>% as.data.frame()
names(castSolutions) <- sort(unique(VCVsolutions$Character_2))
rownames(castSolutions) <- sort(unique(VCVsolutions$Character_2))
VCVsolutions$Best_model <- as.character(VCVsolutions$Best_model)
for(i in 1:length(unique(VCVsolutions$Character_2))){
  for(j in 1:length(unique(VCVsolutions$Character_2))){
    try(castSolutions[i,j]<-VCVsolutions$Best_model[which(VCVsolutions$Character_1==names(castSolutions)[i] & VCVsolutions$Character_2==rownames(castSolutions)[j])],silent=T)
  }
}
castSolutions[is.na(castSolutions)]<-"Non-computable vcv"
ggplot(VCVsolutions, aes(Character_1, Character_2)) + geom_tile(aes(fill = Best_model), colour = "black") + scale_fill_manual(values=c("green", "blue", "orange", "red")) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

bestmodel_types %>% length()
bestmodel_types[which(sapply(bestmodel_types, function(e) e=="no common structure"))] %>% length()
bestmodel_types[which(sapply(bestmodel_types, function(e) e=="common rates, different correlation"))] %>% length()
