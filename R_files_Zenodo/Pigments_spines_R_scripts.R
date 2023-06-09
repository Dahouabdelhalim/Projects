require(brms)
require(ape)
require(phytools)
require(smatr)
require(tidyverse)

#Function to calculate standard error
ste <- function(x) {
  if (length(which(is.na(x)))>0){
    y<-x[-which(is.na(x))]
  } else {
    y<-x
  }
  sd(y)/sqrt(length(y))
}

#set directory for data
dir<-#My directory

# Load data and summarise ---------------------------------------------------------------

d<-read_csv(paste0(dir, "Publication data.csv"))

#Create hatch-level means of all variables
sumf<-dplyr::summarise(group_by(d, Species, Mother, Stage, Family, Uniquemother), RD=mean(RD_length, na.rm=T), RD.se=ste(RD_length), Dorsal=mean(Dorsal_spine, na.rm=T), Dorsal.se=ste(Dorsal_spine), Rostral=mean(Rostral_spine, na.rm=T), Rostral.se=ste(Rostral_spine), Total=mean(Total_spine, na.rm=T), Total.se=ste(Total_spine), Pigprop=mean(Pigment_proportion, na.rm=T), Pigprop.se=ste(Pigment_proportion), Pigarea=mean(Pigment_area, na.rm=T), Pigarea.se=ste(Pigment_area), Clength=mean(Carapace_length, na.rm=T), Clength.se=ste(Carapace_length), Cheight=mean(Carapace_height, na.rm=T), Cheight.se=ste(Carapace_height), Antenna=mean(Antenna_length, na.rm=T), Antenna.se=ste(Antenna_length))

#Summarise hatch-level means into species-level means (for statistical analysis)
sum<-dplyr::summarise(group_by(subset(sumf, Stage==1), Species, Family), RD_length=mean(RD, na.rm=T), RD_length.se=ste(RD), Dorsal_spine=mean(Dorsal, na.rm=T), Dorsal_spine.se=ste(Dorsal), Rostral_spine=mean(Rostral, na.rm=T), Rostral_spine.se=ste(Rostral), Total_spine=mean(Total, na.rm=T), Total_spine.se=ste(Total), Pigment_proportion=mean(Pigprop, na.rm=T), Pigment_proportion.se=ste(Pigprop), Pigment_area=mean(Pigarea, na.rm=T), Pigment_area.se=ste(Pigarea), Carapace_length=mean(Clength, na.rm=T), Carapace_length.se=ste(Clength), Carapace_height=mean(Cheight, na.rm=T), Carapace_height.se=ste(Cheight), Antenna_length=mean(Antenna, na.rm=T), Antenna_length.se=ste(Antenna))

#Create species-level means from raw data
sumtot<-dplyr::summarise(group_by(filter(d, Stage==1), Species, Stage, Family), RD=mean(RD_length, na.rm=T), RD.se=ste(RD_length), Dorsal=mean(Dorsal_spine, na.rm=T), Dorsal.se=ste(Dorsal_spine), Rostral=mean(Rostral_spine, na.rm=T), Rostral.se=ste(Rostral_spine), Total=mean(Total_spine, na.rm=T), Total.se=ste(Total_spine), Total.sd=sd(Total_spine, na.rm=T), Pigprop=mean(Pigment_proportion, na.rm=T), Pigprop.se=ste(Pigment_proportion), Pigprop.sd=sd(Pigment_proportion, na.rm=T), Pigarea=mean(Pigment_area, na.rm=T), Pigarea.se=ste(Pigment_area), Clength=mean(Carapace_length, na.rm=T), Clength.se=ste(Carapace_length), Cheight=mean(Carapace_height, na.rm=T), Cheight.se=ste(Carapace_height), Antenna=mean(Antenna_length, na.rm=T), Antenna.se=ste(Antenna_length))

# Load phylogenetic tree --------------------------------------------------

treeb<-read.nexus(paste0(dir, "infile.nex.con.tre"))
#Make phytools believe tree is ultrametric (even thought it already is. Phytools gets confused by very small rounding change to numbers)
treeb<-force.ultrametric(treeb)
#Change species names of placeholder taxa and update species names for those that have changed names recently
treeb$tip.label[which(treeb$tip.label=="Gecarcinus_lateralis")]<-"Gecarcinus_ruricola"
treeb$tip.label[which(treeb$tip.label=="Pitho_lherminieri")]<-"Pitho_laevigata"
treeb$tip.label[which(treeb$tip.label=="Microphrys_bicornutus")]<-"Omalacantha_bicornuta"
treeb$tip.label[which(treeb$tip.label=="Cancer_antennarius")]<-"Romaleon_antennarium"
treeb$tip.label[which(treeb$tip.label=="Eurypanopeus_abbreviatus")]<-"Eurypanopeus_sp."
#Add missing species to tree
treeb2<-add.species.to.genus(treeb, "Armases_americanum", where = "random")
#Prune tree to include only study species
treeb2.prun<-drop.tip(treeb2, setdiff(treeb2$tip.label, gsub(" ", "_", unique(d$Species))))



# Prepare data and tree for statistical models ----------------------------

#Create covariance matrix from tree
inv.phylo <- MCMCglmm::inverseA(treeb2.prun, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

#Prepare data for model
dstat<-d%>%
  filter(Stage==1)%>% #only include stage 1 larvae (only a few stage 2 larvae collected for 1 species)
  left_join(dplyr::select(ungroup(sum), RD_length, Total_spine, Species, Pigment_proportion, Carapace_length), by=c("Species"), suffix=c("", "_specmean"))%>%
  mutate(Species2=Species,
         Species=gsub(" ", "_", Species), #rename species names to match tree
         Species2=fct_recode(Species2,
                             `A. americanum`="Armases americanum",
                             `A. ricordi`="Armases ricordi",
                             `C. floridanus`="Cataleptodius floridanus",
                             `C. guanhumi`="Cardisoma guanhumi",
                             `C. integer`="Cyclograpsus integer",
                             `G. grapsus`="Grapsus grapsus",
                             `G. lividus`="Geograpsus lividus",
                             `G. ruricola`="Gecarcinus ruricola",
                             `O. bicornuta`="Omalacantha bicornuta",
                             `M. coryphe`="Mithraculus coryphe",
                             `M. forceps`="Mithraculus forceps",
                             `M. sculptus`="Mithraculus sculptus",
                             `P. laevigata`="Pitho laevigata",
                             `P. transversus`="Pachygrapsus transversus",
                             `M. mordax`="Minuca mordax",
                             `M. rapax`="Minuca rapax",
                             `H. oregonensis`="Hemigrapsus oregonensis",
                             `P. crassipes`="Pachygrapsus crassipes",
                             `P. producta`="Pugettia producta",
                             `R. antennarium`="Romaleon antennarium"))%>%
  arrange(Family, Species2)%>%
  mutate(Species2=factor(Species2, unique(Species2)))

#Create variables at scale of mother and species
##Add hatch means to dataset
dstat<-left_join(dstat, dplyr::select(ungroup(sumf), Uniquemother, RD, Total, Pigprop, Clength), by="Uniquemother", suffix=c("", "_mothmean"))%>%
  dplyr::rename(RD_length_mothmean=RD, Total_spine_mothmean=Total, Pigment_proportion_mothmean=Pigprop, Carapace_length_mothmean=Clength)
##Create offset variables at the hatch and individual level so that variation at those levels is decoupled from variation at any other level.
dstat$Pig_moth_spec<-dstat$Pigment_proportion_mothmean-dstat$Pigment_proportion_specmean
dstat$Pig_within_moth<-dstat$Pigment_proportion-dstat$Pigment_proportion_mothmean
#Create phylo variable for phylogenetic analysis
dstat$Phylo<-dstat$Species

#Remove data with NAs for analysis
dstat<-dstat%>%
  filter(!is.na(Total_spine) & !is.na(Pigment_proportion_specmean) & !is.na(Pig_moth_spec) & !is.na(Pig_within_moth) & !is.na(Carapace_length))

#Center and scale predictors, log transform responses
dstat$Carapace_length.s<-(dstat$Carapace_length-mean(dstat$Carapace_length, na.rm=T))/sd(dstat$Carapace_length, na.rm=T)
dstat$Pigment_proportion.s<-(dstat$Pigment_proportion-mean(dstat$Pigment_proportion, na.rm=T))/sd(dstat$Pigment_proportion, na.rm=T)
dstat$Total_spine.l<-log(dstat$Total_spine)
dstat$RD_length.l<-log(dstat$RD_length)
dstat$Pigment_proportion_specmean.s<-(dstat$Pigment_proportion_specmean-mean(dstat$Pigment_proportion_specmean, na.rm=T))/sd(dstat$Pigment_proportion_specmean, na.rm=T)
dstat$Pig_moth_spec.s<-(dstat$Pig_moth_spec)/sd(dstat$Pig_moth_spec, na.rm=T)
dstat$Pig_within_moth.s<-(dstat$Pig_within_moth)/sd(dstat$Pig_within_moth, na.rm=T)

# Fit pigmentation vs. spine length models --------------------------------------------------------------

#Total spine length model
TOTmodelphy<-brm(Total_spine.l ~ (Pigment_proportion_specmean.s+Pig_moth_spec.s+Pig_within_moth.s+Carapace_length.s)^3 + (1|Uniquemother) + (1|Phylo),
                 data=dplyr::select(dstat, Total_spine.l, Pigment_proportion_specmean.s, Pig_moth_spec.s, Pig_within_moth.s, Carapace_length.s, Phylo, Species, Uniquemother),
                 family=gaussian(), cov_ranef = list(Phylo=A), prior=c(prior(normal(0,5), "b"),
                                                                       prior(normal(0,10), "Intercept"),
                                                                       prior(cauchy(0,5), "sd"),
                                                                       prior(cauchy(0,5), "sigma")),
                 chains = 3, cores = 3, 
                 iter = 10000, warmup = 2500)

#Rostral-dorsal length model
RDmodelphy<-brm(RD_length.l ~ (Pigment_proportion_specmean.s+Pig_moth_spec.s+Pig_within_moth.s+Carapace_length.s)^3 + (1|Uniquemother) + (1|Phylo),
                data=dplyr::select(dstat, RD_length.l, Pigment_proportion_specmean.s, Pig_moth_spec.s, Pig_within_moth.s, Carapace_length.s, Phylo, Species, Uniquemother),
                family=gaussian(), cov_ranef = list(Phylo=A), prior=c(prior(normal(0,5), "b"),
                                                                      prior(normal(0,10), "Intercept"),
                                                                      prior(cauchy(0,5), "sd"),
                                                                      prior(cauchy(0,5), "sigma")),
                chains = 3, cores = 3, 
                iter = 10000, warmup = 2500)


# Fit allometric models ---------------------------------------------------

#Intraspecific models
TS<-sma(Total_spine ~ Carapace_length*Species, data=dplyr::select(dstat, Total_spine, Carapace_length, Species), slope.test=1)
RD<-sma(RD_length ~ Carapace_length*Species, data=dplyr::select(dstat, RD_length, Carapace_length, Species), slope.test=1)

#Interspecific models
YRD<-sumtot$RD
names(YRD)<-gsub(" ", "_", sumtot$Species)
RD_PIC<-pic(YRD, treeb2.prun)

YTS<-sumtot$Total
names(YTS)<-gsub(" ", "_", sumtot$Species)
Total_PIC<-pic(YTS, treeb2.prun)

XCL<-sumtot$Clength
names(XCL)<-gsub(" ", "_", sumtot$Species)
Total_PIC<-pic(XCL, treeb2.prun)

TSspecphyRMA<-phyl.RMA(x=XCL, y=YTS, treeb2.prun)
RDspecphyRMA<-phyl.RMA(x=XCL, y=YRD, treeb2.prun)


# Load field data and process for statistics ------------------------------
df<-read_csv(paste0(dir, "Publication field data.csv"))

dfstat<-filter(df, !is.na(Carapacelength) & !is.na(Pigprop) & !is.na(Pigpropav) & !is.na(TSL))

#Create individual offset for individual pigment proportion as we did with the lab data
dfstat$Pigpropavdiff<-dfstat$Pigprop-dfstat$Pigpropav
dfstat$SpecStage=paste(dfstat$Species, dfstat$Stage, sep="_")
dfstat<-dfstat%>%
  group_by(Species)%>%
  summarise(N=n())%>%
  right_join(dfstat)%>%
  mutate(Species2=fct_recode(Species, "Sesarmid1" = "Sesarma1"))

#Log transform response variables
dfstat$TSL.l<-log(dfstat$TSL+1)
dfstat$RostralDorsal.l<-log(dfstat$RostralDorsal)

#Center and scale predictors
dfstat$Pigpropav.s<-(dfstat$Pigpropav-mean(dfstat$Pigpropav))/sd(dfstat$Pigpropav)
dfstat$Pigpropavdiff.s<-(dfstat$Pigpropavdiff)/sd(dfstat$Pigpropavdiff)
dfstat$Carapacelength.s<-(dfstat$Carapacelength-mean(dfstat$Carapacelength))/sd(dfstat$Carapacelength)


# Fit field pigmentation vs. spine length models --------------------------------------------------------

TOTmodelphyfield<-brm(TSL.l ~ Pigpropav.s + Pigpropavdiff.s + Carapacelength.s + (1|Species/Stage), 
                      data=dfstat, family=gaussian(),
                      prior=c(prior(normal(0,5), "b"),
                              prior(normal(0,10), "Intercept"),
                              prior(cauchy(0,5), "sd"),
                              prior(cauchy(0,5), "sigma")),
                      chains=3, cores=3,
                      iter = 1e4, warmup = 2500, control=list(adapt_delta=0.9))

RDmodelphyfield<-brm(RostralDorsal.l ~ Pigpropav.s + Pigpropavdiff.s + Carapacelength.s + (1|Species/Stage), 
                     data=dfstat, family=gaussian(),
                     prior=c(prior(normal(0,5), "b"),
                             prior(normal(0,10), "Intercept"),
                             prior(cauchy(0,5), "sd"),
                             prior(cauchy(0,5), "sigma")),
                     chains = 3, cores = 3, 
                     iter = 1e4, warmup = 2500, control=list(adapt_delta=0.9))
