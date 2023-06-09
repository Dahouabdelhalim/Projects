# r script to the publication 'The costs and benefits of larger brains in fishes'

# loading required packages
library("rfishbase")
library("ape")
library("adephylo")
library("nlme")
library("phylobase")
library("caper")

# loading data (saved as tab-delimited .txt files from the .xlsx file 'FischerJungwirth2022_JEB_S7' as provided on dryad)
all<-read.delim("all.txt")
intraspeciesbrains<-read.delim("8speciesindbrains.txt")
publications<-read.delim("publicationdata.txt")

# re-arranging and managing data
all$RepInvestment<-c(all$MaxFecundity*all$MaxLarvalSize)
all$AsFactorClass<-as.factor(all$Class)
all$AsFactorOrder<-as.factor(all$Order)
all$AsFactorFamily<-as.factor(all$Family)
all$AsFactorGenus<-as.factor(all$Genus)
all$AsFactorSpecies<-as.factor(all$Species)
fish.phylo.formula<- ~AsFactorClass/AsFactorOrder/AsFactorFamily/AsFactorGenus/AsFactorSpecies
phylo.per.class.formula<- ~AsFactorOrder/AsFactorFamily/AsFactorGenus/AsFactorSpecies


# ways to check out what is available on fishbase ----
data(fishbase)
tables()
fishrev.description<-docs()
fishrev.species<-c("Astatotilapia burtoni", "Neolamprologus pulcher", "Nothobranchius furzeri", "Salmo salar")
fishrev.brains<-brains(fishrev.species)
fishrev.commonnames<-common_names(fishrev.species)
fishrev.countref<-countrysub(fishrev.species)
fishrev.country<-coun(fishrev.species)
fishrev.diet<-diet(fishrev.species)
fishrev.dietitems<-diet_items(fishrev.species)
fishrev.distribution<-distribution(fishrev.species)
fishrev.ecology<-ecology(fishrev.species)
fishrev.ecosystem<-ecosystem(fishrev.species)
fishrev.family<-family(fishrev.species)
fishrev.faoareas<-faoareas(fishrev.species)
fishrev.fecundity<-fecundity(fishrev.species)
fishrev.fooditems<-fooditems(fishrev.species)
fishrev.genetics<-genetics(fishrev.species)
fishrev.introductions<-introductions(fishrev.species)
fishrev.larvae<-larvae(fishrev.species)
fishrev.lengthfrequency<-length_freq(fishrev.species)
fishrev.lengthlength<-length_length(fishrev.species)
fishrev.lengthweight<-length_weight(fishrev.species)
fishrev.maturity<-maturity(fishrev.species)
fishrev.morphology<-morphology(fishrev.species)
fishrev.morphometrics<-morphometrics(fishrev.species)
fishrev.oxygen<-oxygen(fishrev.species)
fishrev.popchar<-popchar(fishrev.species)
fishrev.popgrowth<-popgrowth(fishrev.species)
fishrev.popqby<-popqb(fishrev.species)
fishrev.predators<-predators(fishrev.species)
fishrev.ration<-ration(fishrev.species)
fishrev.reproduction<-reproduction(fishrev.species)
fishrev.spawning<-spawning(fishrev.species)
fishrev.speciestable<-species(fishrev.species)
fishrev.speed<-speed(fishrev.species)
fishrev.stocks<-stocks(fishrev.species)
fishrev.swimming<-swimming(fishrev.species)
fishrev.synonyms<-synonyms(fishrev.species)


# getting species data----
fishspecies<-species(fields=c("SpecCode",
                              "Species",
                              "Genus",
                              "FamCode",
                              "Subfamily",
                              "Saltwater",
                              "DemersPelag",
                              "AnaCat",
                              "LongevityWild",
                              "LongevityCaptive",
                              "Length",
                              "LTypeMaxM",
                              "Weight"))
write.table(fishspecies, file="species.txt")
length(subset(fishspecies$SpecCode, !is.na(fishspecies$Subfamily)))
length(subset(fishspecies$SpecCode, !is.na(fishspecies$AnaCat)))
length(subset(fishspecies$SpecCode, !is.na(fishspecies$LongevityWild)))
length(subset(fishspecies$SpecCode, !is.na(fishspecies$LongevityCaptive)))
length(subset(fishspecies$SpecCode, !is.na(fishspecies$LongevityCaptive) & !is.na(fishspecies$LongevityWild)))
length(subset(fishspecies$SpecCode, !is.na(fishspecies$Length)))
length(subset(fishspecies$SpecCode, !is.na(fishspecies$LTypeMaxM)))
length(subset(fishspecies$SpecCode, !is.na(fishspecies$Weight)))

# getting length conversion data----
fishlengthconversions<-length_length(fields=c("SpecCode",
                                              "Species",
                                              "Length1", 
                                              "Length2",
                                              "a",
                                              "b"))
write.table(fishlengthconversions, file="conversions.txt")
length(subset(fishlengthconversions$SpecCode, !is.na(fishlengthconversions$a)))

# getting sociality data-----
fishsociality<-ecology(fields=c("SpecCode",
                                "Species",
                                "Solitary"))
write.table(fishsociality, file="solitary.txt")


# getting brain data----
fishbrains<-brains(fields=c("SpecCode",
                            "Species",
                            "BodyWeight",
                            "BrainWeight",
                            "EncCoeff",
                            "SL",
                            "TL"))
write.table(fishbrains, file="brains.txt")

# getting fecundity data----
fishfecundity<-fecundity(fields=c("SpecCode",
                                  "Species",
                                  "FecundityMin",
                                  "WeightMin",
                                  "FecundityMax",
                                  "WeightMax",
                                  "FecundityMean",
                                  "WeightMean"))
write.table(fishfecundity, file="fecundity.txt")

length(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMin)))
length(unique(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMin))))
length(subset(fishfecundity$SpecCode, !is.na(fishfecundity$WeightMin)))
length(unique(subset(fishfecundity$SpecCode, !is.na(fishfecundity$WeightMin))))

length(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMax)))
length(unique(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMax))))
length(subset(fishfecundity$SpecCode, !is.na(fishfecundity$WeightMax)))
length(unique(subset(fishfecundity$SpecCode, !is.na(fishfecundity$WeightMax))))

length(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMean)))
length(unique(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMean))))
length(subset(fishfecundity$SpecCode, !is.na(fishfecundity$WeightMean)))
length(unique(subset(fishfecundity$SpecCode, !is.na(fishfecundity$WeightMean))))

# getting FecundityMax from FecundityMin
plot(FecundityMax~FecundityMin, data=fishfecundity)
lm.fecundity.mintomax<-lm(FecundityMax~FecundityMin, data=fishfecundity)
summary(lm.fecundity.mintomax)
plot(FecundityMax~FecundityMin, data=fishfecundity, xlim=c(0, 3000000), ylim=c(0, 3000000))
abline(lm.fecundity.mintomax)
### FecunityMax=FecundityMin*3.167+224100

# plotting relationships between maternal weight and fecundity
plot(log(fishfecundity$FecundityMin)~log(fishfecundity$WeightMin), xlab="Minimum Clutch Size [log]", ylab="Absolute Individual Weight [log(g)]")
length(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMin) & !is.na(fishfecundity$WeightMin)))
text(x=12, y=5, c("N_Total: 246"))
length(unique(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMin) & !is.na(fishfecundity$WeightMin))))
text(x=12, y=4, c("N_Species: 174"))

plot(log(fishfecundity$FecundityMax)~log(fishfecundity$WeightMax), xlab="Maximum Clutch Size [log]", ylab="Absolute Individual Weight [log(g)]")
length(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMax) & !is.na(fishfecundity$WeightMax)))
text(x=12, y=5, c("N_Total: 212"))
length(unique(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMax) & !is.na(fishfecundity$WeightMax))))
text(x=12, y=4, c("N_Species: 145"))

plot(log(fishfecundity$FecundityMean)~log(fishfecundity$WeightMean))
length(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMean) & !is.na(fishfecundity$WeightMean)))
text(x=5, y=5, c("N_Total: 5"))
length(unique(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMean) & !is.na(fishfecundity$WeightMean))))
text(x=5, y=4, c("N_Species: 3"))

plot(log(fishfecundity$FecundityMax)~log(fishfecundity$FecundityMin), xlab="Minimum Clutch Size [log]", ylab="Maximum Clutch Size [log]")
length(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMin) & !is.na(fishfecundity$FecundityMax)))
length(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMin)))
length(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMax)))
text(x=15, y=5, c("N_Total: 1970 (min=2300; max=2607)"))
length(unique(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMin) & !is.na(fishfecundity$FecundityMax))))
length(unique(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMin))))
length(unique(subset(fishfecundity$SpecCode, !is.na(fishfecundity$FecundityMax))))
text(x=15, y=4, c("N_Species: 1204 (min=1374; max=1570)"))
abline(a=0, b=1)

# getting social/reproduction data----
fishreproduction<-reproduction(fields=c("SpecCode",
                                        "Species",
                                        "Fertilization",
                                        "MatingSystem",
                                        "RepGuild1",
                                        "RepGuild2",
                                        "ParentalCare"))
write.table(fishreproduction, file="reproduction.txt")

length(subset(fishreproduction$SpecCode, !is.na(fishreproduction$MatingSystem)))
length(subset(fishreproduction$SpecCode, !is.na(fishreproduction$RepGuild1)))
length(subset(fishreproduction$SpecCode, !is.na(fishreproduction$RepGuild2)))
length(subset(fishreproduction$SpecCode, !is.na(fishreproduction$ParentalCare)))

# getting larva data----
fishlarvae<-larvae(fields=c("SpecCode",
                            "Species",
                            "LhMax",
                            "LhMin",
                            "LhMid",
                            "LarvalDurationMin",
                            "LarvalDurationMax",
                            "LarvalDurationMod"))
write.table(fishlarvae, file="larvae.txt")


# larva size
length(subset(fishlarvae$SpecCode, !is.na(fishlarvae$LhMin)))
length(unique(subset(fishlarvae$SpecCode, !is.na(fishlarvae$LhMin))))

length(subset(fishlarvae$SpecCode, !is.na(fishlarvae$LhMax)))
length(unique(subset(fishlarvae$SpecCode, !is.na(fishlarvae$LhMax))))

length(subset(fishlarvae$SpecCode, !is.na(fishlarvae$LhMid)))
length(unique(subset(fishlarvae$SpecCode, !is.na(fishlarvae$LhMid))))

# getting LhMax from LhMin
plot(LhMax~LhMin, data= fishlarvae)
lm.larvae.mintomax<-lm(LhMax~LhMin, data= fishlarvae)
summary(lm.larvae.mintomax)
abline(lm.larvae.mintomax)
### LhMax=LhMin*1.248+0.085

# getting LhMax from LhMid
plot(LhMax~LhMid, data= fishlarvae)
lm.larvae.midtomax<-lm(LhMax~LhMid, data= fishlarvae)
summary(lm.larvae.midtomax)
abline(lm.larvae.midtomax)
### LhMax=LhMin*1.094+0.206

# larval duration
length(subset(fishlarvae$SpecCode, !is.na(fishlarvae$LarvalDurationMin)))
length(unique(subset(fishlarvae$SpecCode, !is.na(fishlarvae$LarvalDurationMin))))

length(subset(fishlarvae$SpecCode, !is.na(fishlarvae$LarvalDurationMax)))
length(unique(subset(fishlarvae$SpecCode, !is.na(fishlarvae$LarvalDurationMax))))

length(subset(fishlarvae$SpecCode, !is.na(fishlarvae$LarvalDurationMod)))
length(unique(subset(fishlarvae$SpecCode, !is.na(fishlarvae$LarvalDurationMod))))

# getting spawning data----
fishspawning<-spawning(fields=c("SpecCode",
                                "Species",
                                "SpawningGround"))
write.table(fishspawning, file="spawning.txt")
length(subset(fishspawning$SpecCode, !is.na(fishspawning$SpawningGround)))
length(unique(subset(fishspawning$SpecCode, !is.na(fishspawning$SpawningGround))))
summary(fishspawning$SpawningGround)


fishrev.allspecies.ecol<-ecology(fields=c("SpecCode",
                                          "Species", 
                                          "Solitary",
                                          "Symbiosis",
                                          "Commensalism",
                                          "Mutualism",
                                          "Schooling",
                                          "SchoolingFrequency",
                                          "SchoolingLifestage",
                                          "Shoaling",
                                          "ShoalingFrequency",
                                          "ShoalingLifestage",
                                          "AssociationsWith"))
write.table(fishrev.allspecies.ecol, "ecol.txt", sep="\\t")

fishrev.allspecies.spec<-species(fields=c("SpecCode",
                                          "Species", 
                                          "LongevityWild",
                                          "LongevityCaptive",
                                          "Length"))
write.table(fishrev.allspecies.spec, "spec.txt", sep="\\t")

fishrev.allspecies.coun<-country(fields=c("SpecCode",
                                          "Species",
                                          "C_Code",
                                          "Importance",
                                          "Aquaculture",
                                          "Threatened"))
write.table(fishrev.allspecies.coun, "coun.txt", sep="\\t")



# Appendix 2.1: Brain size versus body length - inter-specific----

#subset data to make sure no NAs are left
brains.vs.length.data<-subset(all, !is.na(all$BrainWeightMaxWeightIndWithBrain) & !is.na(all$Length))

# make the tree of the reduced (NA free) dataset, adding branch lengths in the second step, adding colours in the third step
brains.vs.length.tree<-as.phylo(fish.phylo.formula, data=brains.vs.length.data, collapse=TRUE)
brains.vs.length.tree.l<-compute.brlen(brains.vs.length.tree, method="Grafen")
plot.phylo(brains.vs.length.tree.l, show.tip.label=FALSE)
plot.phylo(brains.vs.length.tree.l, show.tip.label=FALSE, type="fan")
brains.vs.length.tree.l.actinopteri.species<-which.edge(brains.vs.length.tree.l, subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Actinopteri"))
brains.vs.length.tree.l.cladistii.species<-which.edge(brains.vs.length.tree.l, subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Cladistii"))
brains.vs.length.tree.l.coelacanthi.species<-which.edge(brains.vs.length.tree.l, subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Coelacanthi"))
brains.vs.length.tree.l.dipneusti.species<-which.edge(brains.vs.length.tree.l, subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Dipneusti"))
brains.vs.length.tree.l.elasmobranchii.species<-which.edge(brains.vs.length.tree.l, subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Elasmobranchii"))
brains.vs.length.tree.l.holocephali.species<-which.edge(brains.vs.length.tree.l, subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Holocephali"))
brains.vs.length.tree.l.myxini.species<-which.edge(brains.vs.length.tree.l, subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Myxini"))
brains.vs.length.tree.l.petromyzonti.species<-which.edge(brains.vs.length.tree.l, subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Petromyzonti"))
brains.vs.length.tree.l.classcolours<- rep("darkgrey", dim(brains.vs.length.tree.l$edge)[1]) 
brains.vs.length.tree.l.classcolours[brains.vs.length.tree.l.actinopteri.species] <- "black"
brains.vs.length.tree.l.classcolours[brains.vs.length.tree.l.cladistii.species] <- "red1"
brains.vs.length.tree.l.classcolours[brains.vs.length.tree.l.coelacanthi.species] <- "darkorange"
brains.vs.length.tree.l.classcolours[brains.vs.length.tree.l.dipneusti.species] <- "gold1"
brains.vs.length.tree.l.classcolours[brains.vs.length.tree.l.elasmobranchii.species] <- "green3"
brains.vs.length.tree.l.classcolours[brains.vs.length.tree.l.holocephali.species] <- "aquamarine1"
brains.vs.length.tree.l.classcolours[brains.vs.length.tree.l.myxini.species] <- "steelblue1"
brains.vs.length.tree.l.classcolours[brains.vs.length.tree.l.petromyzonti.species] <- "purple1"
par(mar=c(1, 1, 1, 1))
plot(brains.vs.length.tree.l, lwd=3, edge.color=brains.vs.length.tree.l.classcolours, show.tip.label=FALSE, type="fan")
legend("bottomright", c("Actinopteri", "Cladistii", "Coelacanthi", "Dipneusti", "Elasmobranchii", "Holocephali", "Myxini", "Petromyzonti"),
       col=c("black", "red1", "darkorange", "gold1", "green3", "aquamarine1", "steelblue1", "purple2"),
       lty=c(1, 1, 1, 1, 1, 1, 1, 1),
       lwd=c(2, 2, 2, 2, 2, 2, 2, 2),
       pt.lwd=c(.5, .5, .5, .5, .5, .5, .5, .5), cex=0.75)

# investigating phyologenetic autocorrelation in brain size with correlogram
par(mar=c(3, 5, 1, 1))
fish.phylo.formula.BrainSize<-BrainWeightMaxWeightIndWithBrain~Class/Order/Family/Genus
brains.vs.length.correl.brains<-correlogram.formula(fish.phylo.formula.BrainSize, data=brains.vs.length.data)
plot(brains.vs.length.correl.brains)
text(x=3.5, y=0.7, c("Brain Size"))

# investigating phyologenetic autocorrelation in body length with correlogram
par(mar=c(3, 5, 1, 1))
fish.phylo.formula.Length<-Length~Class/Order/Family/Genus
brains.vs.length.correl.length<-correlogram.formula(fish.phylo.formula.Length, data=brains.vs.length.data)
plot(brains.vs.length.correl.length)
text(x=3.5, y=0.45, c("Body Length"))

# plot both correlograms together
par(mfrow=c(1, 2), mar=c(5, 5, 1, 1))
plot(brains.vs.length.correl.brains)
text(x=3.5, y=0.65, c("Brain Size"))
plot(brains.vs.length.correl.length)
text(x=3.5, y=0.45, c("Body Length"))
par(mfrow=c(1, 1))

# investigating phyolegentic signals in brain size with orthonormal transform
brains.vs.length.tree.orthomatrix<-as.matrix(orthobasis.phylo(brains.vs.length.tree.l))
brains.vs.length.tree.orthomatrix.reduced<-brains.vs.length.tree.orthomatrix[, 1:2]
anova(lm(brains.vs.length.data$BrainWeightMaxWeightIndWithBrain~brains.vs.length.tree.orthomatrix.reduced))
orthogram(brains.vs.length.data$BrainWeightMaxWeightIndWithBrain, brains.vs.length.tree.l)

# investigating phyolegentic signals in body length with orthonormal transform
anova(lm(brains.vs.length.data$Length~brains.vs.length.tree.orthomatrix.reduced))
orthogram(brains.vs.length.data$Length, brains.vs.length.tree.l)

# investigating the phylogenetic correlation between brain size and body length 
brains.vs.length.ppca.matrix<-matrix(c(brains.vs.length.data$BrainWeightMaxWeightIndWithBrain,
                                       brains.vs.length.data$Length), ncol=2)
brains.vs.length.ppca.data<-phylo4d(brains.vs.length.tree.l, brains.vs.length.ppca.matrix)
brains.vs.length.ppca<-ppca(brains.vs.length.ppca.data)
plot(brains.vs.length.ppca)

# computing a phylogenetically controlled generalised least squares model for the correlation between brain size and body length
brains.vs.length.pGLS.brownian<-corBrownian(phy=brains.vs.length.tree.l)
brains.vs.length.pGLS.data<-data.frame(brains.vs.length.data$Species, brains.vs.length.data$BrainWeightMaxWeightIndWithBrain, brains.vs.length.data$Length)
brains.vs.length.pGLS<-gls(brains.vs.length.data.BrainWeightMaxWeightIndWithBrain~brains.vs.length.data.Length, data=brains.vs.length.pGLS.data, correlation=brains.vs.length.pGLS.brownian)
summary(brains.vs.length.pGLS)

# plotting relationships between brain size and body length
length(subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Actinopteri"))
length(subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Elasmobranchii"))
length(subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Petromyzonti"))
length(subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Holocephali"))
length(subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Dipneusti"))
length(subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Myxini"))
length(subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Cladistii"))
length(subset(brains.vs.length.data$Species, brains.vs.length.data$Class=="Coelacanthi"))

par(mfrow=c(1, 1), mar=c(5, 5, 1, 1))
clip(-100, 100, -100, 100)
plot(log(BrainWeightMaxWeightIndWithBrain)~log(Length), data=brains.vs.length.data, col=0, xlab="Body length [ln(cm)]", ylab="Brain weight [ln(mg)]", cex.lab=1.5)
points(log(BrainWeightMaxWeightIndWithBrain)~log(Length), pch=21, bg=0, data=subset(brains.vs.length.data, brains.vs.length.data$Class=="Actinopteri"))
points(log(BrainWeightMaxWeightIndWithBrain)~log(Length), pch=21, bg="green3", data=subset(brains.vs.length.data, brains.vs.length.data$Class=="Elasmobranchii"))
points(log(BrainWeightMaxWeightIndWithBrain)~log(Length), pch=21, bg="purple2", data=subset(brains.vs.length.data, brains.vs.length.data$Class=="Petromyzonti"))
points(log(BrainWeightMaxWeightIndWithBrain)~log(Length), pch=22, bg="aquamarine1", data=subset(brains.vs.length.data, brains.vs.length.data$Class=="Holocephali"))
points(log(BrainWeightMaxWeightIndWithBrain)~log(Length), pch=22, bg="gold1", data=subset(brains.vs.length.data, brains.vs.length.data$Class=="Dipneusti"))
points(log(BrainWeightMaxWeightIndWithBrain)~log(Length), pch=22, bg="steelblue1", data=subset(brains.vs.length.data, brains.vs.length.data$Class=="Myxini"))
points(log(BrainWeightMaxWeightIndWithBrain)~log(Length), pch=22, bg="red1", data=subset(brains.vs.length.data, brains.vs.length.data$Class=="Cladistii"))
points(log(BrainWeightMaxWeightIndWithBrain)~log(Length), pch=24, bg="darkorange", data=subset(brains.vs.length.data, brains.vs.length.data$Class=="Coelacanthi"))

lm.brains.vs.length.all<-lm(log(BrainWeightMaxWeightIndWithBrain)~log(Length), data=brains.vs.length.data)
clip(log(min(brains.vs.length.data$Length)), 
     log(max(brains.vs.length.data$Length)), -100, 100)
abline(lm.brains.vs.length.all, col="grey", lwd=2)
lm.brains.vs.length.actinopteri<-lm(log(BrainWeightMaxWeightIndWithBrain)~log(Length), data=subset(brains.vs.length.data, brains.vs.length.data$Class=="Actinopteri"))
clip(log(min(subset(brains.vs.length.data$Length, brains.vs.length.data$Class=="Actinopteri"))), 
     log(max(subset(brains.vs.length.data$Length, brains.vs.length.data$Class=="Actinopteri"))), -100, 100)
abline(lm.brains.vs.length.actinopteri, lwd=2)
lm.brains.vs.length.elasmobranchii<-lm(log(BrainWeightMaxWeightIndWithBrain)~log(Length), data=subset(brains.vs.length.data, brains.vs.length.data$Class=="Elasmobranchii"))
clip(log(min(subset(brains.vs.length.data$Length, brains.vs.length.data$Class=="Elasmobranchii"))), 
     log(max(subset(brains.vs.length.data$Length, brains.vs.length.data$Class=="Elasmobranchii"))), -100, 100)
abline(lm.brains.vs.length.elasmobranchii, lwd=2, col="green3")
lm.brains.vs.length.holocephali<-lm(log(BrainWeightMaxWeightIndWithBrain)~log(Length), data=subset(brains.vs.length.data, brains.vs.length.data$Class=="Holocephali"))
clip(log(min(subset(brains.vs.length.data$Length, brains.vs.length.data$Class=="Holocephali"))), 
     log(max(subset(brains.vs.length.data$Length, brains.vs.length.data$Class=="Holocephali"))), -100, 100)
abline(lm.brains.vs.length.holocephali, lwd=2, col="aquamarine1")
lm.brains.vs.length.petromyzonti<-lm(log(BrainWeightMaxWeightIndWithBrain)~log(Length), data=subset(brains.vs.length.data, brains.vs.length.data$Class=="Petromyzonti"))
clip(log(min(subset(brains.vs.length.data$Length, brains.vs.length.data$Class=="Petromyzonti"))), 
     log(max(subset(brains.vs.length.data$Length, brains.vs.length.data$Class=="Petromyzonti"))), -100, 100)
abline(lm.brains.vs.length.petromyzonti, lwd=2, col="purple2")

clip(-100, 100, -100, 100)
legend("bottomright", c("Actinopteri", "Cladistii", "Coelacanthi", "Dipneusti", "Elasmobranchii", "Holocephali", "Myxini", "Petromyzonti"),
       pt.bg=c("white", "red1", "darkorange", "gold1", "green3", "aquamarine1", "steelblue1", "purple2"),
       pch=c(21, 22, 24, 22, 21, 22, 22, 21),
       lty=c(0, 0, 0, 0, 0, 0, 0, 0),
       lwd=c(2, 2, 2, 2, 2, 2, 2, 2),
       text.col=c(0, 0, 0, 0, 0, 0, 0, 0),
       pt.lwd=c(.5, .5, .5, .5, .5, .5, .5, .5),
       col=c("black", "red1", "darkorange", "gold1", "green3", "aquamarine1", "steelblue1", "purple2"))
legend("bottomright", c("Actinopteri", "Cladistii", "Coelacanthi", "Dipneusti", "Elasmobranchii", "Holocephali", "Myxini", "Petromyzonti"),
       pt.bg=c("white", "red1", "darkorange", "gold1", "green3", "aquamarine1", "steelblue1", "purple2"),
       pch=c(21, 22, 24, 22, 21, 22, 22, 21),
       lty=c(0, 0, 0, 0, 0, 0, 0, 0),
       lwd=c(2, 2, 2, 2, 2, 2, 2, 2),
       pt.lwd=c(.5, .5, .5, .5, .5, .5, .5, .5),
       bg=NA)
par(mfrow=c(1, 1))

# Appendix 2.2: Brain size versus body length - intra-specific----
par(mfrow=c(1, 1), mar=c(5, 5, 1, 1))
clip(-100, 100, -100, 100)
plot(log(BrainWeight)~log(SL), data=subset(intraspeciesbrains, !is.na(intraspeciesbrains$SL)), col=0, xlab="Body length [ln(cm)]", ylab="Brain weight [ln(mg)]", cex.lab=1.50)
points(log(BrainWeight)~log(SL), data=subset(intraspeciesbrains, intraspeciesbrains$Species=="Poecilia mexicana"), pch=21, bg=0)



# Appendix 2.3: Brain size versus body weight - inter-specific----

#subset data to make sure no NAs are left
brains.vs.weight.data<-subset(all, !is.na(all$BrainWeightMaxWeightIndWithBrain) & !is.na(all$Weight))

# a quick and dirty plot
plot(log(BrainWeightMaxWeightIndWithBrain)~log(Weight), data=brains.vs.weight.data)

# make the tree of the reduced (NA free) dataset, adding branch lengths in the second step
brains.vs.weight.tree<-as.phylo(fish.phylo.formula, data=brains.vs.weight.data, collapse=TRUE)
brains.vs.weight.tree.l<-compute.brlen(brains.vs.weight.tree, method="Grafen")
plot.phylo(brains.vs.weight.tree.l, show.tip.label=FALSE)
plot.phylo(brains.vs.weight.tree.l, show.tip.label=FALSE, type="fan")
brains.vs.weight.tree.l.actinopteri.species<-which.edge(brains.vs.weight.tree.l, subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Actinopteri"))
brains.vs.weight.tree.l.cladistii.species<-which.edge(brains.vs.weight.tree.l, subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Cladistii"))
brains.vs.weight.tree.l.coelacanthi.species<-which.edge(brains.vs.weight.tree.l, subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Coelacanthi"))
brains.vs.weight.tree.l.dipneusti.species<-which.edge(brains.vs.weight.tree.l, subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Dipneusti"))
brains.vs.weight.tree.l.elasmobranchii.species<-which.edge(brains.vs.weight.tree.l, subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Elasmobranchii"))
brains.vs.weight.tree.l.holocephali.species<-which.edge(brains.vs.weight.tree.l, subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Holocephali"))
brains.vs.weight.tree.l.myxini.species<-which.edge(brains.vs.weight.tree.l, subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Myxini"))
brains.vs.weight.tree.l.petromyzonti.species<-which.edge(brains.vs.weight.tree.l, subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Petromyzonti"))
brains.vs.weight.tree.l.classcolours<- rep("darkgrey", dim(brains.vs.weight.tree.l$edge)[1]) 
brains.vs.weight.tree.l.classcolours[brains.vs.weight.tree.l.actinopteri.species] <- "black"
brains.vs.weight.tree.l.classcolours[brains.vs.weight.tree.l.cladistii.species] <- "red1"
brains.vs.weight.tree.l.classcolours[brains.vs.weight.tree.l.coelacanthi.species] <- "darkorange"
brains.vs.weight.tree.l.classcolours[brains.vs.weight.tree.l.dipneusti.species] <- "gold1"
brains.vs.weight.tree.l.classcolours[brains.vs.weight.tree.l.elasmobranchii.species] <- "green3"
brains.vs.weight.tree.l.classcolours[brains.vs.weight.tree.l.holocephali.species] <- "aquamarine1"
brains.vs.weight.tree.l.classcolours[brains.vs.weight.tree.l.myxini.species] <- "steelblue1"
brains.vs.weight.tree.l.classcolours[brains.vs.weight.tree.l.petromyzonti.species] <- "purple1"
par(mar=c(1, 1, 1, 1))
plot(brains.vs.weight.tree.l, lwd=3, edge.color=brains.vs.weight.tree.l.classcolours, show.tip.label=FALSE, type="fan")
legend("bottomright", c("Actinopteri", "Cladistii", "Coelacanthi", "Dipneusti", "Elasmobranchii", "Holocephali", "Myxini", "Petromyzonti"),
       col=c("black", "red1", "darkorange", "gold1", "green3", "aquamarine1", "steelblue1", "purple2"),
       lty=c(1, 1, 1, 1, 1, 1, 1, 1),
       lwd=c(2, 2, 2, 2, 2, 2, 2, 2),
       pt.lwd=c(.5, .5, .5, .5, .5, .5, .5, .5), cex=0.75)

# investigating phyologenetic autocorrelation in brain size with correlogram
fish.phylo.formula.BrainSize<-BrainWeightMaxWeightIndWithBrain~Class/Order/Family/Genus
brains.vs.weight.correl.brains<-correlogram.formula(fish.phylo.formula.BrainSize, data=brains.vs.weight.data)
plot(brains.vs.weight.correl.brains)

# investigating phyologenetic autocorrelation in body weight with correlogram
fish.phylo.formula.Weight<-MaxWeightSpecies~Class/Order/Family/Genus
brains.vs.weight.correl.weight<-correlogram.formula(fish.phylo.formula.Weight, data=brains.vs.weight.data)
plot(brains.vs.weight.correl.weight)

# investigating phyolegentic signals in brain size with orthonormal transform
brains.vs.weight.tree.orthomatrix<-as.matrix(orthobasis.phylo(brains.vs.weight.tree.l))
brains.vs.weight.tree.orthomatrix.reduced<-brains.vs.weight.tree.orthomatrix[, 1:2]
anova(lm(brains.vs.weight.data$BrainWeightMaxWeightIndWithBrain~brains.vs.weight.tree.orthomatrix.reduced))
orthogram(brains.vs.weight.data$BrainWeightMaxWeightIndWithBrain, brains.vs.weight.tree.l)

# investigating phyolegentic signals in body weight with orthonormal transform
anova(lm(brains.vs.weight.data$Length~brains.vs.weight.tree.orthomatrix.reduced))
orthogram(brains.vs.weight.data$Length, brains.vs.weight.tree.l)

# investigating the phylogenetic correlation between brain size and body weight 
brains.vs.weight.ppca.matrix<-matrix(c(brains.vs.weight.data$BrainWeightMaxWeightIndWithBrain,
                                       brains.vs.weight.data$MaxWeightSpecies), ncol=2)
brains.vs.weight.ppca.data<-phylo4d(brains.vs.weight.tree.l, brains.vs.weight.ppca.matrix)
brains.vs.weight.ppca<-ppca(brains.vs.weight.ppca.data)
plot(brains.vs.weight.ppca)


# computing a phylogenetically controlled generalised least squares model for the correlation between brain size and body weight
brains.vs.weight.pGLS.brownian<-corBrownian(phy=brains.vs.weight.tree.l)
brains.vs.weight.pGLS.data<-data.frame(brains.vs.weight.data$Species, brains.vs.weight.data$BrainWeightMaxWeightIndWithBrain, brains.vs.weight.data$MaxWeightSpecies)
brains.vs.weight.pGLS<-gls(brains.vs.weight.data.BrainWeightMaxWeightIndWithBrain~brains.vs.weight.data.MaxWeightSpecies, data=brains.vs.weight.pGLS.data, correlation=brains.vs.weight.pGLS.brownian)
summary(brains.vs.weight.pGLS)


# plotting relationships between brain size and body weight
length(subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Actinopteri"))
length(subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Elasmobranchii"))
length(subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Petromyzonti"))
length(subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Holocephali"))
length(subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Dipneusti"))
length(subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Myxini"))
length(subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Cladistii"))
length(subset(brains.vs.weight.data$Species, brains.vs.weight.data$Class=="Coelacanthi"))

par(mfrow=c(1, 1), mar=c(5, 5, 1, 1))
clip(-100, 100, -100, 100)
plot(logBrain~logWeight, data=brains.vs.weight.data, col=0, xlab="Body weight [ln(g)]", ylab="Brain weight [ln(mg)]", cex.lab=1.5, ylim=c(0, 12.5))
points(logBrain~logWeight, pch=21, bg=0, data=subset(brains.vs.weight.data, brains.vs.weight.data$Class=="Actinopteri"))
points(logBrain~logWeight, pch=21, bg="green3", data=subset(brains.vs.weight.data, brains.vs.weight.data$Class=="Elasmobranchii"))
points(logBrain~logWeight, pch=21, bg="purple2", data=subset(brains.vs.weight.data, brains.vs.weight.data$Class=="Petromyzonti"))
#points(logBrain~logWeight, pch=22, bg="aquamarine1", data=subset(brains.vs.weight.data, brains.vs.weight.data$Class=="Holocephali"))
#points(logBrain~logWeight, pch=22, bg="gold1", data=subset(brains.vs.weight.data, brains.vs.weight.data$Class=="Dipneusti"))
#points(logBrain~logWeight, pch=22, bg="steelblue1", data=subset(brains.vs.weight.data, brains.vs.weight.data$Class=="Myxini"))
#points(logBrain~logWeight, pch=22, bg="red1", data=subset(brains.vs.weight.data, brains.vs.weight.data$Class=="Cladistii"))
#points(logBrain~logWeight, pch=24, bg="darkorange", data=subset(brains.vs.weight.data, brains.vs.weight.data$Class=="Coelacanthi"))

lm.brains.vs.weight.all<-lm(logBrain~logWeight, data=brains.vs.weight.data)
clip(log(min(brains.vs.weight.data$Weight)), 
     log(max(brains.vs.weight.data$Weight)), -100, 100)
#abline(lm.brains.vs.weight.all, col="grey", lwd=2)
lm.brains.vs.weight.actinopteri<-lm(logBrain~logWeight, data=subset(brains.vs.weight.data, brains.vs.weight.data$Class=="Actinopteri"))
summary(lm.brains.vs.weight.actinopteri)
clip(log(min(subset(brains.vs.weight.data$Weight, brains.vs.weight.data$Class=="Actinopteri"))), 
     log(max(subset(brains.vs.weight.data$Weight, brains.vs.weight.data$Class=="Actinopteri"))), -100, 100)
abline(lm.brains.vs.weight.actinopteri, lwd=2)
lm.brains.vs.weight.elasmobranchii<-lm(logBrain~logWeight, data=subset(brains.vs.weight.data, brains.vs.weight.data$Class=="Elasmobranchii"))
summary(lm.brains.vs.weight.elasmobranchii)
clip(log(min(subset(brains.vs.weight.data$Weight, brains.vs.weight.data$Class=="Elasmobranchii"))), 
     log(max(subset(brains.vs.weight.data$Weight, brains.vs.weight.data$Class=="Elasmobranchii"))), -100, 100)
abline(lm.brains.vs.weight.elasmobranchii, lwd=2, col="green3")
lm.brains.vs.weight.holocephali<-lm(logBrain~logWeight, data=subset(brains.vs.weight.data, brains.vs.weight.data$Class=="Holocephali"))
clip(log(min(subset(brains.vs.weight.data$Weight, brains.vs.weight.data$Class=="Holocephali"))), 
     log(max(subset(brains.vs.weight.data$Weight, brains.vs.weight.data$Class=="Holocephali"))), -100, 100)
#abline(lm.brains.vs.weight.holocephali, lwd=2, col="aquamarine1")
lm.brains.vs.weight.petromyzonti<-lm(logBrain~logWeight, data=subset(brains.vs.weight.data, brains.vs.weight.data$Class=="Petromyzonti"))
summary(lm.brains.vs.weight.petromyzonti)
clip(log(min(subset(brains.vs.weight.data$Weight, brains.vs.weight.data$Class=="Petromyzonti"))), 
     log(max(subset(brains.vs.weight.data$Weight, brains.vs.weight.data$Class=="Petromyzonti"))), -100, 100)
abline(lm.brains.vs.weight.petromyzonti, lwd=2, col="purple2")

clip(-100, 100, -100, 100)
legend("bottomright", c("Actinopteri", "Elasmobranchii", "Petromyzonti"),
       pt.bg=c("white", "green3", "purple2"),
       pch=c(21, 21, 21),
       lty=c(1, 1, 1),
       lwd=c(2, 2, 2),
       text.col=c(0, 0, 0),
       pt.lwd=c(.5, .5, .5),
       col=c("black", "green3", "purple2"))
legend("bottomright", c("Actinopteri", "Elasmobranchii", "Petromyzonti"),
       pt.bg=c("white", "green3", "purple2"),
       pch=c(21, 21, 21),
       lty=c(0, 0, 0),
       lwd=c(2, 2, 2),
       pt.lwd=c(.5, .5, .5),
       bg=NA)

# Appendix 2.4: Brain size versus body weight - intra-specific----
par(mfrow=c(1, 1), mar=c(5, 5, 1, 1))
clip(-100, 100, -100, 100)
plot(log(BrainWeight)~log(BodyWeight), data=intraspeciesbrains, col=0, xlab="Body weight [ln(g)]", ylab="Brain weight [ln(mg)]", cex.lab=1.50)
points(log(BrainWeight)~log(BodyWeight), data=subset(intraspeciesbrains, intraspeciesbrains$Species=="Poecilia mexicana"), pch=21, bg=0)
points(log(BrainWeight)~log(BodyWeight), data=subset(intraspeciesbrains, intraspeciesbrains$Species=="Taurulus bubalis"), pch=22, bg=0)
points(log(BrainWeight)~log(BodyWeight), data=subset(intraspeciesbrains, intraspeciesbrains$Species=="Geotria australis"), pch=21, bg="purple2")
points(log(BrainWeight)~log(BodyWeight), data=subset(intraspeciesbrains, intraspeciesbrains$Species=="Torpedo marmorata"), pch=21, bg="green3")
points(log(BrainWeight)~log(BodyWeight), data=subset(intraspeciesbrains, intraspeciesbrains$Species=="Caspiomyzon wagneri"), pch=22, bg="purple2")
points(log(BrainWeight)~log(BodyWeight), data=subset(intraspeciesbrains, intraspeciesbrains$Species=="Polypterus palmas"), pch=21, bg="red1")
points(log(BrainWeight)~log(BodyWeight), data=subset(intraspeciesbrains, intraspeciesbrains$Species=="Latimeria chalumnae"), pch=21, bg="darkorange")
points(log(BrainWeight)~log(BodyWeight), data=subset(intraspeciesbrains, intraspeciesbrains$Species=="Carcharias taurus"), pch=22, bg="green3")
legend("bottomright", c("Poecilia mexicana", "Taurulus bubalis", "Polypterus palmas", "Latimeria chalumnae", "Carcharias taurus", "Torpedo marmorata", "Caspiomyzon wagneri", "Geotria australis"),
       pt.bg=c("white", "white", "red1", "darkorange", "green3", "green3", "purple2", "purple2"),
       pch=c(21, 22, 21, 21, 22, 21, 22, 21),
       lty=c(0, 0, 0, 0, 0, 0, 0, 0),
       lwd=c(2, 2, 2, 2, 2, 2, 2, 2),
       text.col=c(0, 0, 0, 0, 0, 0, 0, 0),
       pt.lwd=c(.5, .5, .5, .5, .5, .5, .5, .5),
       col=c("black", "red1", "darkorange", "gold1", "green3", "aquamarine1", "steelblue1", "purple2"))
legend("bottomright", c("Poecilia mexicana", "Taurulus bubalis", "Polypterus palmas", "Latimeria chalumnae", "Carcharias taurus", "Torpedo marmorata", "Caspiomyzon wagneri", "Geotria australis"),
       pt.bg=c("white", "white", "red1", "darkorange", "green3", "green3", "purple2", "purple2"),
       pch=c(21, 22, 21, 21, 22, 21, 22, 21),
       lty=c(0, 0, 0, 0, 0, 0, 0, 0),
       lwd=c(2, 2, 2, 2, 2, 2, 2, 2),
       pt.lwd=c(.5, .5, .5, .5, .5, .5, .5, .5),
       bg=NA)

# Appendix 3 - COSTS.BRAINS. investigating potential costs of having a large brain across each class separately----

# subsetting all data to remove any NAs
costs.brains.data<-subset(all, !is.na(all$Longevity) & !is.na(all$MaxFecundity) & !is.na(all$BrainWeightMaxWeightIndWithBrain))

# subsetting NA free data for the three classes in which regressions can be made
costs.brains.data.actinopteri<-subset(costs.brains.data, costs.brains.data$Class=="Actinopteri")
length(subset(costs.brains.data$Species, costs.brains.data$Class=="Actinopteri"))
costs.brains.data.elasmobranchii<-subset(costs.brains.data, costs.brains.data$Class=="Elasmobranchii")
length(subset(costs.brains.data$Species, costs.brains.data$Class=="Elasmobranchii"))
costs.brains.data.petromyzonti<-subset(costs.brains.data, costs.brains.data$Class=="Petromyzonti")
length(subset(costs.brains.data$Species, costs.brains.data$Class=="Petromyzonti"))

# making Z values of all the variables of interest for the three classes
costs.brains.data.actinopteri$BrainZ<-c((costs.brains.data.actinopteri$logBrain-mean(costs.brains.data.actinopteri$logBrain))/sd(costs.brains.data.actinopteri$logBrain))
costs.brains.data.elasmobranchii$BrainZ<-c((costs.brains.data.elasmobranchii$logBrain-mean(costs.brains.data.elasmobranchii$logBrain))/sd(costs.brains.data.elasmobranchii$logBrain))
costs.brains.data.petromyzonti$BrainZ<-c((costs.brains.data.petromyzonti$logBrain-mean(costs.brains.data.petromyzonti$logBrain))/sd(costs.brains.data.petromyzonti$logBrain))
costs.brains.data.actinopteri$WeightZ<-c((costs.brains.data.actinopteri$logWeight-mean(costs.brains.data.actinopteri$logWeight))/sd(costs.brains.data.actinopteri$logWeight))
costs.brains.data.elasmobranchii$WeightZ<-c((costs.brains.data.elasmobranchii$logWeight-mean(costs.brains.data.elasmobranchii$logWeight))/sd(costs.brains.data.elasmobranchii$logWeight))
costs.brains.data.petromyzonti$WeightZ<-c((costs.brains.data.petromyzonti$logWeight-mean(costs.brains.data.petromyzonti$logWeight))/sd(costs.brains.data.petromyzonti$logWeight))
costs.brains.data.actinopteri$LongevityZ<-c((costs.brains.data.actinopteri$logLongevity-mean(costs.brains.data.actinopteri$logLongevity))/sd(costs.brains.data.actinopteri$logLongevity))
costs.brains.data.elasmobranchii$LongevityZ<-c((costs.brains.data.elasmobranchii$logLongevity-mean(costs.brains.data.elasmobranchii$logLongevity))/sd(costs.brains.data.elasmobranchii$logLongevity))
costs.brains.data.petromyzonti$LongevityZ<-c((costs.brains.data.petromyzonti$logLongevity-mean(costs.brains.data.petromyzonti$logLongevity))/sd(costs.brains.data.petromyzonti$logLongevity))
costs.brains.data.actinopteri$FecundityZ<-c((costs.brains.data.actinopteri$logFecundity-mean(costs.brains.data.actinopteri$logFecundity))/sd(costs.brains.data.actinopteri$logFecundity))
costs.brains.data.elasmobranchii$FecundityZ<-c((costs.brains.data.elasmobranchii$logFecundity-mean(costs.brains.data.elasmobranchii$logFecundity))/sd(costs.brains.data.elasmobranchii$logFecundity))
costs.brains.data.petromyzonti$FecundityZ<-c((costs.brains.data.petromyzonti$logFecundity-mean(costs.brains.data.petromyzonti$logFecundity))/sd(costs.brains.data.petromyzonti$logFecundity))

# making trees for each class in which regressions can be made
# Actinopteri trees
costs.brains.tree.actinopteri<-as.phylo(phylo.per.class.formula, data=costs.brains.data.actinopteri, collapse=TRUE)
costs.brains.tree.actinopteri.l<-compute.brlen(costs.brains.tree.actinopteri, method="Grafen")
par(mar=c(1, 1, 1, 1), mfrow=c(1, 2))
plot.phylo(costs.brains.tree.actinopteri.l, show.tip.label=FALSE)
plot.phylo(costs.brains.tree.actinopteri.l, show.tip.label=FALSE, type="fan")
# Figure A3.1
par(mar=c(0, 12, 0, 4), mfrow=c(1, 1))
plot.phylo(costs.brains.tree.actinopteri.l, show.tip.label=TRUE, cex=0.3)

# Elasmobranchii trees
costs.brains.tree.elasmobranchii<-as.phylo(phylo.per.class.formula, data=costs.brains.data.elasmobranchii, collapse=TRUE)
costs.brains.tree.elasmobranchii.l<-compute.brlen(costs.brains.tree.elasmobranchii, method="Grafen")
par(mar=c(1, 1, 1, 1), mfrow=c(1, 2))
plot.phylo(costs.brains.tree.elasmobranchii.l, show.tip.label=FALSE, edge.color=c("green3"))
plot.phylo(costs.brains.tree.elasmobranchii.l, show.tip.label=FALSE, type="fan", edge.color=c("green3"))
# Figure A2.2
par(mar=c(0, 6, 0, 1), mfrow=c(1, 1))
plot.phylo(costs.brains.tree.elasmobranchii.l, show.tip.label=TRUE, edge.color=c("green3"))

#Petromyzonti trees
costs.brains.tree.petromyzonti<-as.phylo(phylo.per.class.formula, data=costs.brains.data.petromyzonti, collapse=TRUE)
costs.brains.tree.petromyzonti.l<-compute.brlen(costs.brains.tree.petromyzonti, method="Grafen")
par(mar=c(1, 1, 1, 1), mfrow=c(1, 2))
plot.phylo(costs.brains.tree.petromyzonti.l, show.tip.label=FALSE, edge.color=c("purple2"))
plot.phylo(costs.brains.tree.petromyzonti.l, show.tip.label=FALSE, type="fan", edge.color=c("purple2"))
# Figure A3.3
par(mar=c(0, 6, 0, 1), mfrow=c(1, 1))
plot.phylo(costs.brains.tree.petromyzonti.l, show.tip.label=TRUE, edge.color=c("purple2"))

# investigating phyologenetic autocorrelation in absolute brain size, body weight, fecundity and longevity
phylo.signal.formula.costs.brains.brain<-BrainWeightMaxWeightIndWithBrain~Order/Family/Genus
phylo.signal.costs.brains.brain.actinopteri<-correlogram.formula(phylo.signal.formula.costs.brains.brain, data=costs.brains.data.actinopteri)
phylo.signal.formula.costs.brains.weight<-Weight~Order/Family/Genus
phylo.signal.costs.brains.weight.actinopteri<-correlogram.formula(phylo.signal.formula.costs.brains.weight, data=costs.brains.data.actinopteri)
phylo.signal.formula.costs.brains.fecundity<-MaxFecundity~Order/Family/Genus
phylo.signal.costs.brains.fecundity.actinopteri<-correlogram.formula(phylo.signal.formula.costs.brains.fecundity, data=costs.brains.data.actinopteri)
phylo.signal.formula.costs.brains.longevity<-Longevity~Order/Family/Genus
phylo.signal.costs.brains.longevity.actinopteri<-correlogram.formula(phylo.signal.formula.costs.brains.longevity, data=costs.brains.data.actinopteri)
phylo.signal.formula.costs.brains.brain<-BrainWeightMaxWeightIndWithBrain~Order/Family/Genus
phylo.signal.costs.brains.brain.elasmobranchii<-correlogram.formula(phylo.signal.formula.costs.brains.brain, data=costs.brains.data.elasmobranchii)
phylo.signal.formula.costs.brains.weight<-Weight~Order/Family/Genus
phylo.signal.costs.brains.weight.elasmobranchii<-correlogram.formula(phylo.signal.formula.costs.brains.weight, data=costs.brains.data.elasmobranchii)
phylo.signal.formula.costs.brains.fecundity<-MaxFecundity~Order/Family/Genus
phylo.signal.costs.brains.fecundity.elasmobranchii<-correlogram.formula(phylo.signal.formula.costs.brains.fecundity, data=costs.brains.data.elasmobranchii)
phylo.signal.formula.costs.brains.longevity<-Longevity~Order/Family/Genus
phylo.signal.costs.brains.longevity.elasmobranchii<-correlogram.formula(phylo.signal.formula.costs.brains.longevity, data=costs.brains.data.elasmobranchii)
phylo.signal.formula.costs.brains.brain<-BrainWeightMaxWeightIndWithBrain~Genus
phylo.signal.costs.brains.brain.petromyzonti<-correlogram.formula(phylo.signal.formula.costs.brains.brain, data=costs.brains.data.petromyzonti)
phylo.signal.formula.costs.brains.weight<-Weight~Genus
phylo.signal.costs.brains.weight.petromyzonti<-correlogram.formula(phylo.signal.formula.costs.brains.weight, data=costs.brains.data.petromyzonti)
phylo.signal.formula.costs.brains.fecundity<-MaxFecundity~Genus
phylo.signal.costs.brains.fecundity.petromyzonti<-correlogram.formula(phylo.signal.formula.costs.brains.fecundity, data=costs.brains.data.petromyzonti)
phylo.signal.formula.costs.brains.longevity<-Longevity~Genus
phylo.signal.costs.brains.longevity.petromyzonti<-correlogram.formula(phylo.signal.formula.costs.brains.longevity, data=costs.brains.data.petromyzonti)

# Figure A3.4
par(mfrow=c(4, 3), mar=c(2, 5, 1, 1))
plot(phylo.signal.costs.brains.brain.actinopteri, ylim=c(0, 0.5))
legend("bottomleft", c("A"), bty="n")
plot(phylo.signal.costs.brains.brain.elasmobranchii)
legend("bottomleft", c("B"), bty="n", text.col="green3")
plot(phylo.signal.costs.brains.brain.petromyzonti)
legend("bottomleft", c("C"), bty="n", text.col="purple3")
plot(phylo.signal.costs.brains.weight.actinopteri)
legend("bottomleft", c("D"), bty="n")
plot(phylo.signal.costs.brains.weight.elasmobranchii, ylim=c(-0.1, 0.4))
legend("bottomleft", c("E"), bty="n", text.col="green3")
plot(phylo.signal.costs.brains.weight.petromyzonti)
legend("bottomleft", c("F"), bty="n", text.col="purple3")
plot(phylo.signal.costs.brains.fecundity.actinopteri)
legend("bottomleft", c("G"), bty="n")
plot(phylo.signal.costs.brains.fecundity.elasmobranchii)
legend("bottomleft", c("H"), bty="n", text.col="green3")
plot(phylo.signal.costs.brains.fecundity.petromyzonti)
legend("bottomleft", c("I"), bty="n", text.col="purple3")
plot(phylo.signal.costs.brains.longevity.actinopteri, ylim=c(-0.2, 0.7))
legend("bottomleft", c("J"), bty="n")
plot(phylo.signal.costs.brains.longevity.elasmobranchii)
legend("bottomleft", c("K"), bty="n", text.col="green3")
plot(phylo.signal.costs.brains.longevity.petromyzonti)
legend("bottomleft", c("L"), bty="n", text.col="purple3")
par(mfrow=c(1, 1), mar=c(5, 5, 1, 1))

# computing phylogenetically controlled generalised least squares models for each class in which regressions can be made

# actinopteri - ln(absolute)
costs.brains.pGLS.brownian.actinopteri<-corBrownian(phy=costs.brains.tree.actinopteri.l)
costs.brains.pGLS.model.actinopteri<-gls(logBrain~logFecundity+logLongevity+logWeight, data=costs.brains.data.actinopteri, correlation=costs.brains.pGLS.brownian.actinopteri)
# Table A3.5
summary(costs.brains.pGLS.model.actinopteri)
# Figure A3.5
par(mfrow=c(1, 3), mar=c(5, 5, 1, 1))
plot(resid(costs.brains.pGLS.model.actinopteri)~fitted(costs.brains.pGLS.model.actinopteri), xlab="Fitted values", ylab="Residuals")
legend("topleft", c("A"), bty="n", cex=1.5)
abline(h=0, lty=2)
hist(costs.brains.pGLS.model.actinopteri$residuals, main="Costly brain - Actinopteri", xlab="Residuals")
legend("topleft", c("B"), bty="n", cex=1.5)
qqnorm(costs.brains.pGLS.model.actinopteri$residuals, main="")
qqline(costs.brains.pGLS.model.actinopteri$residuals)
legend("topleft", c("C"), bty="n", cex=1.5)
costs.brains.pGLS.model.actinopteri.residuals<-resid(costs.brains.pGLS.model.actinopteri)
shapiro.test(costs.brains.pGLS.model.actinopteri.residuals)

# actinopteri - Z values of ln transformed data
costs.brains.pGLS.brownian.actinopteri<-corBrownian(phy=costs.brains.tree.actinopteri.l)
costs.brains.pGLS.model.actinopteriZ<-gls(BrainZ~FecundityZ+LongevityZ+WeightZ, data=costs.brains.data.actinopteri, correlation=costs.brains.pGLS.brownian.actinopteri)
# Table A3.6
summary(costs.brains.pGLS.model.actinopteriZ)
# Figure A3.6
par(mfrow=c(1, 3), mar=c(5, 5, 1, 1))
plot(resid(costs.brains.pGLS.model.actinopteriZ)~fitted(costs.brains.pGLS.model.actinopteriZ), xlab="Fitted values", ylab="Residuals")
abline(h=0, lty=2)
hist(costs.brains.pGLS.model.actinopteriZ$residuals, main="Costly brain - Actinopteri Z", xlab="Residuals")
qqnorm(costs.brains.pGLS.model.actinopteriZ$residuals, main="")
qqline(costs.brains.pGLS.model.actinopteriZ$residuals)
costs.brains.pGLS.model.actinopteriZ.residuals<-resid(costs.brains.pGLS.model.actinopteriZ)
shapiro.test(costs.brains.pGLS.model.actinopteriZ.residuals)

# elasmobranchii - ln(absolute values)
costs.brains.pGLS.brownian.elasmobranchii<-corBrownian(phy=costs.brains.tree.elasmobranchii.l)
costs.brains.pGLS.model.elasmobranchii<-gls(logBrain~logFecundity+logLongevity+logWeight, data=costs.brains.data.elasmobranchii, correlation=costs.brains.pGLS.brownian.elasmobranchii)
# Table A3.7
summary(costs.brains.pGLS.model.elasmobranchii)
# Figure A3.7
plot(resid(costs.brains.pGLS.model.elasmobranchii)~fitted(costs.brains.pGLS.model.elasmobranchii), xlab="Fitted values", ylab="Residuals")
abline(h=0, lty=2)
hist(costs.brains.pGLS.model.elasmobranchii$residuals, main="Costly brain - Elasmobranchii", xlab="Residuals")
qqnorm(costs.brains.pGLS.model.elasmobranchii$residuals, main="")
qqline(costs.brains.pGLS.model.elasmobranchii$residuals)
costs.brains.pGLS.model.elasmobranchii.residuals<-resid(costs.brains.pGLS.model.elasmobranchii)
shapiro.test(costs.brains.pGLS.model.elasmobranchii.residuals)

# elasmobranchii - Z values of ln transformed data
costs.brains.pGLS.brownian.elasmobranchii<-corBrownian(phy=costs.brains.tree.elasmobranchii.l)
costs.brains.pGLS.model.elasmobranchiiZ<-gls(BrainZ~FecundityZ+LongevityZ+WeightZ, data=costs.brains.data.elasmobranchii, correlation=costs.brains.pGLS.brownian.elasmobranchii)
# Table A3.8
summary(costs.brains.pGLS.model.elasmobranchiiZ)
# Figure A3.8
par(mfrow=c(1, 3), mar=c(5, 5, 1, 1))
plot(resid(costs.brains.pGLS.model.elasmobranchiiZ)~fitted(costs.brains.pGLS.model.elasmobranchiiZ), xlab="Fitted values", ylab="Residuals")
abline(h=0, lty=2)
hist(costs.brains.pGLS.model.elasmobranchiiZ$residuals, main="Costly brain - Elasmobranchii Z", xlab="Residuals")
qqnorm(costs.brains.pGLS.model.elasmobranchiiZ$residuals, main="")
qqline(costs.brains.pGLS.model.elasmobranchiiZ$residuals)
costs.brains.pGLS.model.elasmobranchiiZ.residuals<-resid(costs.brains.pGLS.model.elasmobranchiiZ)
shapiro.test(costs.brains.pGLS.model.elasmobranchiiZ.residuals)

# petromyzonti - ln(absolute values)
costs.brains.pGLS.brownian.petromyzonti<-corBrownian(phy=costs.brains.tree.petromyzonti.l)
costs.brains.pGLS.model.petromyzonti<-gls(logBrain~logFecundity+logLongevity+logWeight, data=costs.brains.data.petromyzonti, correlation=costs.brains.pGLS.brownian.petromyzonti)
# Table A3.9
summary(costs.brains.pGLS.model.petromyzonti)
# Figure A3.9
plot(resid(costs.brains.pGLS.model.petromyzonti)~fitted(costs.brains.pGLS.model.petromyzonti), xlab="Fitted values", ylab="Residuals")
abline(h=0, lty=2)
hist(costs.brains.pGLS.model.petromyzonti$residuals, main="Costly brain - Petromyzonti", xlab="Residuals")
qqnorm(costs.brains.pGLS.model.petromyzonti$residuals, main="")
qqline(costs.brains.pGLS.model.petromyzonti$residuals)
costs.brains.pGLS.model.petromyzonti.residuals<-resid(costs.brains.pGLS.model.petromyzonti)
shapiro.test(costs.brains.pGLS.model.petromyzonti.residuals)
par(mfrow=c(1, 1), mar=c(5, 5, 1, 1))

# petromyzonti - Z values of ln transformed data
costs.brains.pGLS.brownian.petromyzonti<-corBrownian(phy=costs.brains.tree.petromyzonti.l)
costs.brains.pGLS.model.petromyzontiZ<-gls(BrainZ~FecundityZ+LongevityZ+WeightZ, data=costs.brains.data.petromyzonti, correlation=costs.brains.pGLS.brownian.petromyzonti)
# Table A3.10
summary(costs.brains.pGLS.model.petromyzontiZ)
# Figure A3.10
par(mfrow=c(1, 3), mar=c(5, 5, 1, 1))
plot(resid(costs.brains.pGLS.model.petromyzontiZ)~fitted(costs.brains.pGLS.model.petromyzontiZ), xlab="Fitted values", ylab="Residuals")
abline(h=0, lty=2)
hist(costs.brains.pGLS.model.petromyzontiZ$residuals, main="Costly brain - Petromyzonti Z", xlab="Residuals")
qqnorm(costs.brains.pGLS.model.petromyzontiZ$residuals, main="")
qqline(costs.brains.pGLS.model.petromyzontiZ$residuals)
costs.brains.pGLS.model.petromyzontiZ.residuals<-resid(costs.brains.pGLS.model.petromyzontiZ)
shapiro.test(costs.brains.pGLS.model.petromyzontiZ.residuals)

# proper plot ALL - absolute (logged) values
layout(matrix(c(1, 1, 2, 3, 4, 5), nrow=3, ncol=2, byrow = T), heights=c(1, 4, 4, 4, 4))
layout.show(5)
par(mar=c(1, 5, 1, 1))
plot(logBrain~logFecundity, data=costs.brains.data, xlab="", ylab="", xaxt="n", yaxt="n", col=0)
legend("left", bty="n", c("Actinopteri", "Coelacanthi"),
       pch=c(21, 24), pt.bg=c(0, "darkorange"), cex=1)
legend("center", bty="n", c("Dipneusti", "Elasmobranchi"),
       pch=c(22, 21), pt.bg=c("gold1", "green3"), cex=1)
legend("right", bty="n", c("Petromyzonti"),
       pch=c(21), pt.bg=c("purple3"), cex=1)
par(mar=c(5, 5, 1, 1))
plot(logBrain~logFecundity, data=costs.brains.data, col=0, xlab="Fecundity [ln(eggs)]", ylab="Brain size [ln(mg)]", cex.lab=1.5)
points(logBrain~logFecundity, data=costs.brains.data.actinopteri)
#lines(x=c(5.5, 19), y=c((5.5*(-0.0045)+3.5343), (19*(-0.0045)+3.5343)), lwd=2, lty=2)
points(logBrain~logFecundity, data=costs.brains.data.elasmobranchii, pch=21, bg="green3")
#lines(x=c(0.5, 4.5), y=c((0.5*(-0.4139)+7.225), (4.5*(-0.4139)+7.225)), lwd=2, lty=1, col="green3")
points(logBrain~logFecundity, data=costs.brains.data.petromyzonti, pch=21, bg="purple3")
#lines(x=c(8.5, 13), y=c((8.5*(-0.2146)+3.8082), (13*(-0.2146)+3.8082)), lwd=2, lty=2, col="purple3")
points(logBrain~logFecundity, data=subset(costs.brains.data, costs.brains.data$Class=="Coelacanthi"), pch=24, bg="darkorange")
points(logBrain~logFecundity, data=subset(costs.brains.data, costs.brains.data$Class=="Dipneusti"), pch=22, bg="gold1")
legend("topright", bty="n", c("A"), cex=1.5)

plot(logBrain~logLongevity, data=costs.brains.data, col=0, xlab="Lifespan [ln(years)]", ylab="", cex.lab=1.5)
points(logBrain~logLongevity, data=costs.brains.data.actinopteri)
#lines(x=c(0.05, 4.5), y=c((0.05*(-0.0068)+3.5343), (4.5*(-0.0068)+3.5343)), lwd=2, lty=2)
points(logBrain~logLongevity, data=costs.brains.data.elasmobranchii, pch=21, bg="green3")
#lines(x=c(1.7, 6), y=c((1.7*(-0.0039)+7.225), (6*(-0.0039)+7.225)), lwd=2, lty=1, col="green3")
points(logBrain~logLongevity, data=costs.brains.data.petromyzonti, pch=21, bg="purple3")
#lines(x=c(1.5, 2.5), y=c((1.5*(-0.1214)+3.8082), (2.5*(-0.1214)+3.8082)), lwd=2, lty=2, col="purple3")
points(logBrain~logLongevity, data=subset(costs.brains.data, costs.brains.data$Class=="Coelacanthi"), pch=24, bg="darkorange")
points(logBrain~logLongevity, data=subset(costs.brains.data, costs.brains.data$Class=="Dipneusti"), pch=22, bg="gold1")
legend("topright", bty="n", c("B"), cex=1.5)

# proper plot - Z values
#layout(matrix(c(1, 1, 2, 3, 2, 3), nrow=2, ncol=2, byrow = T), heights=c(1, 5))
#par(mar=c(1, 5, 1, 1))
#plot(log(BrainWeightMaxWeightIndWithBrain)~log(MaxFecundity), data=costs.brains.data, xlab="", ylab="", xaxt="n", yaxt="n", col=0)
#legend("left", bty="n", c("Actinopteri", "Coelacanthi"), pch=c(21, 24), pt.bg=c(0, "darkorange"), cex=1.5)
#legend("center", bty="n", c("Dipneusti", "Elasmobranchi"), pch=c(22, 21), pt.bg=c("gold1", "green3"), cex=1.5)
#legend("right", bty="n", c("Petromyzonti"), pch=c(21), pt.bg=c("purple3"), cex=1.5)
#par(mar=c(5, 5, 1, 1))
plot(BrainZ~FecundityZ, data=costs.brains.data.actinopteri, col=0, xlab="Fecundity [Z]", ylab="Brain size [Z]", cex.lab=1.5, ylim=c(-2.75, 2.25), xlim=c(-2.75, 2.75))
points(BrainZ~FecundityZ, data=costs.brains.data.actinopteri)
clip(min(costs.brains.data.actinopteri$FecundityZ), max(costs.brains.data.actinopteri$FecundityZ), -5, 5)
abline(costs.brains.pGLS.model.actinopteriZ$coefficients[["(Intercept)"]], costs.brains.pGLS.model.actinopteriZ$coefficients[["FecundityZ"]], lty=2, lwd=2)
clip(-100, 100, -100, 100)
points(BrainZ~FecundityZ, data=costs.brains.data.elasmobranchii, pch=21, bg="green3")
clip(min(costs.brains.data.elasmobranchii$FecundityZ), max(costs.brains.data.elasmobranchii$FecundityZ), -5, 5)
abline(costs.brains.pGLS.model.elasmobranchiiZ$coefficients[["(Intercept)"]], costs.brains.pGLS.model.elasmobranchiiZ$coefficients[["FecundityZ"]], col="green3", lty=1, lwd=2)
clip(-100, 100, -100, 100)
points(BrainZ~FecundityZ, data=costs.brains.data.petromyzonti, pch=21, bg="purple3")
clip(min(costs.brains.data.petromyzonti$FecundityZ), max(costs.brains.data.petromyzonti$FecundityZ), -5, 5)
abline(costs.brains.pGLS.model.petromyzontiZ$coefficients[["(Intercept)"]], costs.brains.pGLS.model.petromyzontiZ$coefficients[["FecundityZ"]], col="purple3", lty=2, lwd=2)
clip(-100, 100, -100, 100)
legend("topright", bty="n", c("C"), cex=1.5)

plot(BrainZ~LongevityZ, data=costs.brains.data.actinopteri, col=0, xlab="Lifespan [Z]", ylab="", cex.lab=1.5, ylim=c(-2.75, 2.25), xlim=c(-2.75, 4))
points(BrainZ~LongevityZ, data=costs.brains.data.actinopteri)
clip(min(costs.brains.data.actinopteri$LongevityZ), max(costs.brains.data.actinopteri$LongevityZ), -5, 5)
abline(costs.brains.pGLS.model.actinopteriZ$coefficients[["(Intercept)"]], costs.brains.pGLS.model.actinopteriZ$coefficients[["LongevityZ"]], lty=2, lwd=2)
clip(-100, 100, -100, 100)
points(BrainZ~LongevityZ, data=costs.brains.data.elasmobranchii, pch=21, bg="green3")
clip(min(costs.brains.data.elasmobranchii$LongevityZ), max(costs.brains.data.elasmobranchii$LongevityZ), -5, 5)
abline(costs.brains.pGLS.model.elasmobranchiiZ$coefficients[["(Intercept)"]], costs.brains.pGLS.model.elasmobranchiiZ$coefficients[["LongevityZ"]], col="green3", lty=2, lwd=2)
clip(-100, 100, -100, 100)
points(BrainZ~LongevityZ, data=costs.brains.data.petromyzonti, pch=21, bg="purple3")
clip(min(costs.brains.data.petromyzonti$LongevityZ), max(costs.brains.data.petromyzonti$LongevityZ), -5, 5)
abline(costs.brains.pGLS.model.petromyzontiZ$coefficients[["(Intercept)"]], costs.brains.pGLS.model.petromyzontiZ$coefficients[["LongevityZ"]], col="purple3", lty=1, lwd=2)
clip(-100000, 100000, -100000, 100000)
legend("topright", bty="n", c("D"), cex=1.5)

# Fig A3.11
par(mar=c(5, 5, 1, 1), mfrow=c(1, 2))
plot(logBrain~logFecundity, data=costs.brains.data, col=0, xlab="Fecundity [ln(eggs)]", ylab="Brain size [ln(mg)]", cex.lab=1.5)
points(logBrain~logFecundity, data=costs.brains.data.actinopteri)
#lines(x=c(5.5, 19), y=c((5.5*(-0.0045)+3.5343), (19*(-0.0045)+3.5343)), lwd=2, lty=2)
points(logBrain~logFecundity, data=costs.brains.data.elasmobranchii, pch=21, bg="green3")
#lines(x=c(0.5, 4.5), y=c((0.5*(-0.4139)+7.225), (4.5*(-0.4139)+7.225)), lwd=2, lty=1, col="green3")
points(logBrain~logFecundity, data=costs.brains.data.petromyzonti, pch=21, bg="purple3")
#lines(x=c(8.5, 13), y=c((8.5*(-0.2146)+3.8082), (13*(-0.2146)+3.8082)), lwd=2, lty=2, col="purple3")
points(logBrain~logFecundity, data=subset(costs.brains.data, costs.brains.data$Class=="Coelacanthi"), pch=24, bg="darkorange")
points(logBrain~logFecundity, data=subset(costs.brains.data, costs.brains.data$Class=="Dipneusti"), pch=22, bg="gold1")
legend("topright", bty="n", c("A"), cex=1.5)

plot(logBrain~logLongevity, data=costs.brains.data, col=0, xlab="Lifespan [ln(years)]", ylab="", cex.lab=1.5, xlim=c(0, 6))
points(logBrain~logLongevity, data=costs.brains.data.actinopteri)
#lines(x=c(0.05, 4.5), y=c((0.05*(-0.0068)+3.5343), (4.5*(-0.0068)+3.5343)), lwd=2, lty=2)
points(logBrain~logLongevity, data=costs.brains.data.elasmobranchii, pch=21, bg="green3")
#lines(x=c(1.7, 6), y=c((1.7*(-0.0039)+7.225), (6*(-0.0039)+7.225)), lwd=2, lty=1, col="green3")
points(logBrain~logLongevity, data=costs.brains.data.petromyzonti, pch=21, bg="purple3")
#lines(x=c(1.5, 2.5), y=c((1.5*(-0.1214)+3.8082), (2.5*(-0.1214)+3.8082)), lwd=2, lty=2, col="purple3")
points(logBrain~logLongevity, data=subset(costs.brains.data, costs.brains.data$Class=="Coelacanthi"), pch=24, bg="darkorange")
points(logBrain~logLongevity, data=subset(costs.brains.data, costs.brains.data$Class=="Dipneusti"), pch=22, bg="gold1")
legend("topright", bty="n", c("B"), cex=1.5)
clip(-100, 100, -100, 100)
legend("bottomright", c("Actinopteri", "Coelacanthi", "Dipneusti", "Elasmobranchii", "Petromyzonti"),
       pt.bg=c("white", "darkorange", "gold1", "green3", "purple2"),
       pch=c(21, 24, 22, 21, 21),
       cex=.9)
par(mfrow=c(1, 1))

# Figure 2 COSTS OF BRAINS
par(mar=c(5, 5, 1, 1), mfrow=c(3, 2))
par(mar=c(3, 5, 2, 1))
plot(BrainZ~FecundityZ, data=costs.brains.data.actinopteri, col=0, xlab="", ylab="Brain size [Z]", cex.lab=1.5, ylim=c(-2.75, 2.25), xlim=c(-2.75, 2.75))
points(BrainZ~FecundityZ, data=costs.brains.data.actinopteri)
clip(min(costs.brains.data.actinopteri$FecundityZ), max(costs.brains.data.actinopteri$FecundityZ), -5, 5)
abline(costs.brains.pGLS.model.actinopteriZ$coefficients[["(Intercept)"]], costs.brains.pGLS.model.actinopteriZ$coefficients[["FecundityZ"]], lty=2, lwd=2)
clip(-100, 100, -100, 100)
text(x=-2.7, y=1.9, c("A"), cex=2)
#legend("topleft", bty="n", c("A"), cex=1.5)

par(mar=c(3, 5, 2, 1))
plot(BrainZ~LongevityZ, data=costs.brains.data.actinopteri, col=0, xlab="", ylab="", cex.lab=1.5, ylim=c(-2.75, 2.25), xlim=c(-2.75, 4))
points(BrainZ~LongevityZ, data=costs.brains.data.actinopteri)
clip(min(costs.brains.data.actinopteri$LongevityZ), max(costs.brains.data.actinopteri$LongevityZ), -5, 5)
abline(costs.brains.pGLS.model.actinopteriZ$coefficients[["(Intercept)"]], costs.brains.pGLS.model.actinopteriZ$coefficients[["LongevityZ"]], lty=2, lwd=2)
clip(-100, 100, -100, 100)
text(x=-2.7, y=1.9, c("B"), cex=2)
#legend("topleft", bty="n", c("B"), cex=1.5)

par(mar=c(4, 5, 1, 1))
plot(BrainZ~FecundityZ, data=costs.brains.data.elasmobranchii, pch=21, bg="green3", xlab="", ylab="Brain size [Z]", cex.lab=1.5, ylim=c(-2.75, 2.25), xlim=c(-2.75, 2.75))
clip(min(costs.brains.data.elasmobranchii$FecundityZ), max(costs.brains.data.elasmobranchii$FecundityZ), -5, 5)
abline(costs.brains.pGLS.model.elasmobranchiiZ$coefficients[["(Intercept)"]], costs.brains.pGLS.model.elasmobranchiiZ$coefficients[["FecundityZ"]], col="green3", lty=1, lwd=2)
clip(-100, 100, -100, 100)
text(x=-2.7, y=-2, c("C"), cex=2)
#legend("bottomleft", bty="n", c("C"), cex=1.5)

par(mar=c(4, 5, 1, 1))
plot(BrainZ~LongevityZ, data=costs.brains.data.elasmobranchii, pch=21, bg="green3", xlab="", ylab="", cex.lab=1.5, ylim=c(-2.75, 2.25), xlim=c(-2.75, 4))
clip(min(costs.brains.data.elasmobranchii$LongevityZ), max(costs.brains.data.elasmobranchii$LongevityZ), -5, 5)
abline(costs.brains.pGLS.model.elasmobranchiiZ$coefficients[["(Intercept)"]], costs.brains.pGLS.model.elasmobranchiiZ$coefficients[["LongevityZ"]], col="green3", lty=2, lwd=2)
clip(-100, 100, -100, 100)
text(x=-2.7, y=1.9, c("D"), cex=2)
#legend("topleft", bty="n", c("D"), cex=1.5)

par(mar=c(5, 5, 0, 1))
plot(BrainZ~FecundityZ, data=costs.brains.data.petromyzonti, pch=21, bg="purple3",xlab="Fecundity [Z]", ylab="Brain size [Z]", cex.lab=1.5, ylim=c(-2.75, 2.25), xlim=c(-2.75, 2.75))
clip(min(costs.brains.data.petromyzonti$FecundityZ), max(costs.brains.data.petromyzonti$FecundityZ), -5, 5)
abline(costs.brains.pGLS.model.petromyzontiZ$coefficients[["(Intercept)"]], costs.brains.pGLS.model.petromyzontiZ$coefficients[["FecundityZ"]], col="purple3", lty=2, lwd=2)
clip(-100, 100, -100, 100)
text(x=-2.7, y=1.9, c("E"), cex=2)
#legend("topleft", bty="n", c("E"), cex=1.5)

par(mar=c(5, 5, 0, 1))
plot(BrainZ~LongevityZ, data=costs.brains.data.petromyzonti, pch=21, bg="purple3", xlab="Lifespan [Z]", ylab="", cex.lab=1.5, ylim=c(-2.75, 2.25), xlim=c(-2.75, 4))
clip(min(costs.brains.data.petromyzonti$LongevityZ), max(costs.brains.data.petromyzonti$LongevityZ), -5, 5)
abline(costs.brains.pGLS.model.petromyzontiZ$coefficients[["(Intercept)"]], costs.brains.pGLS.model.petromyzontiZ$coefficients[["LongevityZ"]], col="purple3", lty=1, lwd=2)
clip(-100000, 100000, -100000, 100000)
text(x=-2.7, y=1.9, c("F"), cex=2)
#legend("topleft", bty="n", c("F"), cex=1.5)
#legend("bottomright", c("Actinopteri", "Elasmobranchii", "Petromyzonti"),
#       pt.bg=c("white", "green3", "purple2"),
#       pch=c(21, 21,21),
#       lty=c(1, 1, 1),
#       lwd=c(2, 2, 2),
#       text.col=c(0, 0, 0),
#       pt.lwd=c(.5, .5, .5),
#       col=c("black", "green3", "purple2"),
#       cex=0.9)
#legend("bottomright", c("Actinopteri", "Elasmobranchii", "Petromyzonti"),
#       pt.bg=c("white", "green3", "purple2"),
#       pch=c(21, 21, 21),
#       lty=c(0, 0, 0),
#       lwd=c(2, 2, 2),
#       pt.lwd=c(.5, .5, .5),
#       bg=NA,
#       cex=0.9)
par(mfrow=c(1, 1))

# Appendix 4 - REPCAR.BRAINS. investigating potential links between breeding system and brain size ----

#subset data to remove all NAs
repcar.brains.data<-subset(all, !is.na(all$logBrain) & !is.na(all$logWeight) & !is.na(all$Monogamy) & !is.na(all$ParentalCare))
length(subset(repcar.brains.data$Species, repcar.brains.data$Class=="Actinopteri"))
length(subset(repcar.brains.data$Species, repcar.brains.data$Class=="Elasmobranchii"))
length(subset(repcar.brains.data$Species, repcar.brains.data$Class=="Petromyzonti"))

# quick plots
par(mar=c(5, 5, 1, 1), mfrow=c(1, 2))
boxplot(logBrain~Monogamy*Class, data= repcar.brains.data)
boxplot(logBrain~ParentalCare*Class, data=repcar.brains.data)

# subset to only have data for actinopteri
repcar.brains.data.actinopteri<-subset(repcar.brains.data, repcar.brains.data$Class=="Actinopteri")
length(repcar.brains.data.actinopteri$Species)

# making Z values 
repcar.brains.data.actinopteri$BrainZ<-c((repcar.brains.data.actinopteri$logBrain-mean(repcar.brains.data.actinopteri$logBrain))/sd(repcar.brains.data.actinopteri$logBrain))
repcar.brains.data.actinopteri$WeightZ<-c((repcar.brains.data.actinopteri$logWeight-mean(repcar.brains.data.actinopteri$logWeight))/sd(repcar.brains.data.actinopteri$logWeight))

# make trees
repcar.brains.tree.actinopteri<-as.phylo(phylo.per.class.formula, data=repcar.brains.data.actinopteri, collapse=TRUE)
repcar.brains.tree.actinopteri.l<-compute.brlen(repcar.brains.tree.actinopteri, method="Grafen")
par(mar=c(1, 1, 1, 1), mfrow=c(1, 2))
plot.phylo(repcar.brains.tree.actinopteri.l, show.tip.label=FALSE)
plot.phylo(repcar.brains.tree.actinopteri.l, show.tip.label=FALSE, type="fan")
#Figure A4.1
par(mar=c(0, 12, 0, 4), mfrow=c(1, 1))
plot.phylo(repcar.brains.tree.actinopteri.l, show.tip.label=TRUE, cex=0.2)

# investigating phyologenetic autocorrelation in absolute brain size, body weight, mating system, and parental care
phylo.signal.formula.repcar.brains.brain<-BrainWeightMaxWeightIndWithBrain~Order/Family/Genus
phylo.signal.repcar.brains.brain.actinopteri<-correlogram.formula(phylo.signal.formula.repcar.brains.brain, data=repcar.brains.data.actinopteri)
phylo.signal.formula.repcar.brains.weight<-Weight~Order/Family/Genus
phylo.signal.repcar.brains.weight.actinopteri<-correlogram.formula(phylo.signal.formula.repcar.brains.weight, data=repcar.brains.data.actinopteri)
phylo.signal.formula.repcar.brains.monogamy<-Monogamy~Order/Family/Genus
phylo.signal.repcar.brains.monogamy.actinopteri<-correlogram.formula(phylo.signal.formula.repcar.brains.monogamy, data=repcar.brains.data.actinopteri)
phylo.signal.formula.repcar.brains.parental<-ParentalCare~Order/Family/Genus
phylo.signal.repcar.brains.parental.actinopteri<-correlogram.formula(phylo.signal.formula.repcar.brains.parental, data=repcar.brains.data.actinopteri)

# Figure A4.2
par(mfrow=c(4, 1), mar=c(2, 5, 1, 1))
plot(phylo.signal.repcar.brains.brain.actinopteri, ylim=c(0, .7))
legend("bottomleft", c("A"), bty="n")
plot(phylo.signal.repcar.brains.weight.actinopteri)
legend("bottomleft", c("B"), bty="n")
plot(phylo.signal.repcar.brains.monogamy.actinopteri)
legend("bottomleft", c("C"), bty="n")
plot(phylo.signal.repcar.brains.parental.actinopteri, ylim=c(.6, 1.1))
legend("bottomleft", c("D"), bty="n")

# computing phylogenetically controlled generalised least squares model - Actinopteri absolute values
repcar.brains.pGLS.brownian.actinopteri<-corBrownian(phy=repcar.brains.tree.actinopteri.l)
repcar.brains.pGLS.model.actinopteri<-gls(logBrain~Monogamy+ParentalCare+logWeight, data=repcar.brains.data.actinopteri, correlation=repcar.brains.pGLS.brownian.actinopteri)
# Table A4.3
summary(repcar.brains.pGLS.model.actinopteri)
# Figure A4.3
par(mfrow=c(1, 3), mar=c(5, 5, 1, 1))
plot(resid(repcar.brains.pGLS.model.actinopteri)~fitted(repcar.brains.pGLS.model.actinopteri), xlab="Fitted values", ylab="Residuals")
abline(h=0, lty=2)
hist(repcar.brains.pGLS.model.actinopteri$residuals, main="Breed sys - Actinopteri", xlab="Residuals")
qqnorm(repcar.brains.pGLS.model.actinopteri$residuals, main="")
qqline(repcar.brains.pGLS.model.actinopteri$residuals)
repcar.brains.pGLS.model.actinopteri.residuals<-resid(repcar.brains.pGLS.model.actinopteri)
shapiro.test(repcar.brains.pGLS.model.actinopteri.residuals)

# computing phylogenetically controlled generalised least squares model - Actinopteri Z values
repcar.brains.pGLS.brownian.actinopteri<-corBrownian(phy=repcar.brains.tree.actinopteri.l)
repcar.brains.pGLS.model.actinopteriZ<-gls(BrainZ~Monogamy+ParentalCare+WeightZ, data=repcar.brains.data.actinopteri, correlation=repcar.brains.pGLS.brownian.actinopteri)
# Table 4.4
summary(repcar.brains.pGLS.model.actinopteriZ)
# Figure A4.4
par(mfrow=c(1, 3), mar=c(5, 5, 1, 1))
plot(resid(repcar.brains.pGLS.model.actinopteriZ)~fitted(repcar.brains.pGLS.model.actinopteriZ), xlab="Fitted values", ylab="Residuals")
abline(h=0, lty=2)
hist(repcar.brains.pGLS.model.actinopteriZ$residuals, main="Breed sys Z - Actinopteri ", xlab="Residuals")
qqnorm(repcar.brains.pGLS.model.actinopteriZ$residuals, main="")
qqline(repcar.brains.pGLS.model.actinopteriZ$residuals)
repcar.brains.pGLS.model.actinopteriZ.residuals<-resid(repcar.brains.pGLS.model.actinopteriZ)
shapiro.test(repcar.brains.pGLS.model.actinopteriZ.residuals)

# Figure A4.5
par(mar=c(5, 5, 1, 0), mfrow=c(1, 2))
boxplot(logBrain~Monogamy, data=repcar.brains.data.actinopteri, ylab="Brain size [ln(mg)]", cex.lab=1.5, ylim=c(0, 8), xaxt="n", col=c("pink1", "pink4"))
axis(side=1, at=c(1, 2), mgp=c(3, 1, 0), c("no", "yes"))
text(x=c(1), y=c(0.5), c(length(subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$Monogamy=="0"))))
text(x=c(2), y=c(0.5), c(length(subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$Monogamy=="1"))))
text(x=0.6, y=7.75, c("A"), cex=1.5)
par(mar=c(5, 4, 1, 1))
boxplot(logBrain~ParentalCare, data=repcar.brains.data.actinopteri, ylab="", cex.lab=1.5, ylim=c(0, 8), xaxt="n", col=c("royalblue1", "royalblue4"), xlab="Parental care")
axis(side=1, at=c(1, 2), mgp=c(3, 1, 0), c("no", "yes"))
text(x=c(1), y=c(0.5), c(length(subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$ParentalCare=="0"))))
text(x=c(2), y=c(0.5), c(length(subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$ParentalCare=="1"))))
text(x=0.6, y=7.75, c("B"), cex=1.5)

# proper plot Z
par(mar=c(5, 5, 1, 0), mfrow=c(1, 2))
boxplot(BrainZ~Monogamy, data=repcar.brains.data.actinopteri, ylab="Brain size [Z]", cex.lab=1.5, ylim=c(-3.5, 1.5), xaxt="n", col=c("pink1", "pink4"))
axis(side=1, at=c(1, 2), mgp=c(3, 1, 0), c("no", "yes"))
text(x=c(1), y=c(-3.5), c(length(subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$Monogamy=="0"))))
text(x=c(2), y=c(-3.5), c(length(subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$Monogamy=="1"))))
text(x=0.6, y=1.5, c("A"), cex=1.5)
par(mar=c(5, 4, 1, 1))
boxplot(BrainZ~ParentalCare, data=repcar.brains.data.actinopteri, ylab="", cex.lab=1.5, ylim=c(-3.5, 1.5), xaxt="n", col=c("royalblue1", "royalblue4"), xlab="Parental care")
axis(side=1, at=c(1, 2), mgp=c(3, 1, 0), c("no", "yes"))
text(x=c(1), y=c(-3.5), c(length(subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$ParentalCare=="0"))))
text(x=c(2), y=c(-3.5), c(length(subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$ParentalCare=="1"))))
text(x=0.6, y=1.5, c("B"), cex=1.5)

# counts of breeding systems
length(subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$Monogamy=="1" & repcar.brains.data.actinopteri$ParentalCare=="1"))
length(subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$Monogamy=="1" & repcar.brains.data.actinopteri$ParentalCare=="0"))
length(subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$Monogamy=="0" & repcar.brains.data.actinopteri$ParentalCare=="1"))
length(subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$Monogamy=="0" & repcar.brains.data.actinopteri$ParentalCare=="0"))


# coloured trees
repcar.brains.tree.actinopteri.l.monogamous<-which.edge(repcar.brains.tree.actinopteri.l, subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$Monogamy=="1"))
repcar.brains.tree.actinopteri.l.monogamycolours<- rep("pink1", dim(repcar.brains.tree.actinopteri.l$edge)[1]) 
repcar.brains.tree.actinopteri.l.monogamycolours[repcar.brains.tree.actinopteri.l.monogamous] <- "pink4"
par(mar=c(0 ,2, 0, 0))
plot(repcar.brains.tree.actinopteri.l, lwd=3, edge.color=repcar.brains.tree.actinopteri.l.monogamycolours, show.tip.label=FALSE)

repcar.brains.tree.actinopteri.l.caring<-which.edge(repcar.brains.tree.actinopteri.l, subset(repcar.brains.data.actinopteri$Species, repcar.brains.data.actinopteri$ParentalCare=="1"))
repcar.brains.tree.actinopteri.l.carecolours<- rep("royalblue1", dim(repcar.brains.tree.actinopteri.l$edge)[1]) 
repcar.brains.tree.actinopteri.l.carecolours[repcar.brains.tree.actinopteri.l.caring] <- "royalblue4"
par(mar=c(0 ,0, 0, 0))
plot(repcar.brains.tree.actinopteri.l, lwd=3, edge.color=repcar.brains.tree.actinopteri.l.carecolours, show.tip.label=FALSE, direction = "leftwards")

# additional plot - weight effects
par(mfrow=c(1, 2), mar=c(5, 5, 1, 1))
plot(BrainZ~WeightZ, data=repcar.brains.data.actinopteri, col=0)
points(BrainZ~WeightZ, data=subset(repcar.brains.data.actinopteri, repcar.brains.data.actinopteri$Monogamy=="1"), bg=c("pink4"), pch=21)
points(BrainZ~WeightZ, data=subset(repcar.brains.data.actinopteri, repcar.brains.data.actinopteri$Monogamy=="0"), bg=c("pink1"), pch=21)
abline(repcar.brains.pGLS.model.actinopteriZ$coefficients[["(Intercept)"]], repcar.brains.pGLS.model.actinopteriZ$coefficients[["WeightZ"]], lty=1, lwd=2)
plot(BrainZ~WeightZ, data=repcar.brains.data.actinopteri, col=0)
points(BrainZ~WeightZ, data=subset(repcar.brains.data.actinopteri, repcar.brains.data.actinopteri$ParentalCare=="1"), bg=c("royalblue4"), pch=21)
points(BrainZ~WeightZ, data=subset(repcar.brains.data.actinopteri, repcar.brains.data.actinopteri$ParentalCare=="0"), bg=c("royalblue1"), pch=21)
abline(repcar.brains.pGLS.model.actinopteriZ$coefficients[["(Intercept)"]], repcar.brains.pGLS.model.actinopteriZ$coefficients[["WeightZ"]], lty=1, lwd=2)

# Appendix 5 - HABITAT.BRAINS. investigating potential links between environmental complexity and brain size----

# subsetting all data to remove any NAs
habitat.brains.data<-subset(all, !is.na(all$Migratory) & !is.na(all$Pelagic) & !is.na(all$BrainWeightMaxWeightIndWithBrain))

# quick plot
par(mfrow=c(1, 2), mar=c(5, 5, 1, 1))
boxplot(log(BrainWeightMaxWeightIndWithBrain)~Migratory*Class, data=habitat.brains.data)
boxplot(log(BrainWeightMaxWeightIndWithBrain)~Pelagic*Class, data=habitat.brains.data)

# checking sample sizes
length(subset(habitat.brains.data$Species, habitat.brains.data$Class=="Actinopteri"))
length(subset(habitat.brains.data$Species, habitat.brains.data$Class=="Cladistii"))
length(subset(habitat.brains.data$Species, habitat.brains.data$Class=="Coelacanthi"))
length(subset(habitat.brains.data$Species, habitat.brains.data$Class=="Dipneusti"))
length(subset(habitat.brains.data$Species, habitat.brains.data$Class=="Elasmobranchii"))
length(subset(habitat.brains.data$Species, habitat.brains.data$Class=="Holocephali"))
length(subset(habitat.brains.data$Species, habitat.brains.data$Class=="Myxini"))
length(subset(habitat.brains.data$Species, habitat.brains.data$Class=="Petromyzonti"))

# making data sets for those classes for which regressions can be made
habitat.brains.data.actinopteri<-subset(habitat.brains.data, habitat.brains.data$Class=="Actinopteri")
habitat.brains.data.elasmobranchii<-subset(habitat.brains.data, habitat.brains.data$Class=="Elasmobranchii")
habitat.brains.data.petromyzonti<-subset(habitat.brains.data, habitat.brains.data$Class=="Petromyzonti")

# making Z values for all the variables of interest
habitat.brains.data.actinopteri$BrainZ<-c((habitat.brains.data.actinopteri$logBrain-mean(habitat.brains.data.actinopteri$logBrain))/sd(habitat.brains.data.actinopteri$logBrain))
habitat.brains.data.actinopteri$WeightZ<-c((habitat.brains.data.actinopteri$logWeight-mean(habitat.brains.data.actinopteri$logWeight))/sd(habitat.brains.data.actinopteri$logWeight))
habitat.brains.data.elasmobranchii$BrainZ<-c((habitat.brains.data.elasmobranchii$logBrain-mean(habitat.brains.data.elasmobranchii$logBrain))/sd(habitat.brains.data.elasmobranchii$logBrain))
habitat.brains.data.elasmobranchii$WeightZ<-c((habitat.brains.data.elasmobranchii$logWeight-mean(habitat.brains.data.elasmobranchii$logWeight))/sd(habitat.brains.data.elasmobranchii$logWeight))
habitat.brains.data.petromyzonti$BrainZ<-c((habitat.brains.data.petromyzonti$logBrain-mean(habitat.brains.data.petromyzonti$logBrain))/sd(habitat.brains.data.petromyzonti$logBrain))
habitat.brains.data.petromyzonti$WeightZ<-c((habitat.brains.data.petromyzonti$logWeight-mean(habitat.brains.data.petromyzonti$logWeight))/sd(habitat.brains.data.petromyzonti$logWeight))

# making trees for each class in which regressions can be made
habitat.brains.tree.actinopteri<-as.phylo(phylo.per.class.formula, data=habitat.brains.data.actinopteri, collapse=TRUE)
habitat.brains.tree.actinopteri.l<-compute.brlen(habitat.brains.tree.actinopteri, method="Grafen")
par(mar=c(1, 1, 1, 1), mfrow=c(1, 2))
plot.phylo(habitat.brains.tree.actinopteri.l, show.tip.label=FALSE)
plot.phylo(habitat.brains.tree.actinopteri.l, show.tip.label=FALSE, type="fan")
habitat.brains.tree.elasmobranchii<-as.phylo(phylo.per.class.formula, data=habitat.brains.data.elasmobranchii, collapse=TRUE)
habitat.brains.tree.elasmobranchii.l<-compute.brlen(habitat.brains.tree.elasmobranchii, method="Grafen")
# Figure A5.1
par(mar=c(0, 5, 0, 5), mfrow=c(1, 1))
abitat.brains.tree.elasmobranchii<-as.phylo(phylo.per.class.formula, data=habitat.brains.data.elasmobranchii, collapse=TRUE)
plot.phylo(habitat.brains.tree.actinopteri.l, show.tip.label=TRUE, cex=0.1)

par(mar=c(1, 1, 1, 1), mfrow=c(1, 2))
plot.phylo(habitat.brains.tree.elasmobranchii.l, show.tip.label=FALSE, edge.color=c("green3"))
plot.phylo(habitat.brains.tree.elasmobranchii.l, show.tip.label=FALSE, type="fan", edge.color=c("green3"))
habitat.brains.tree.petromyzonti<-as.phylo(phylo.per.class.formula, data=habitat.brains.data.petromyzonti, collapse=TRUE)
habitat.brains.tree.petromyzonti.l<-compute.brlen(habitat.brains.tree.petromyzonti, method="Grafen")
par(mar=c(1, 1, 1, 1), mfrow=c(1, 2))
plot.phylo(habitat.brains.tree.petromyzonti.l, show.tip.label=FALSE, edge.color=c("purple2"))
plot.phylo(habitat.brains.tree.petromyzonti.l, show.tip.label=FALSE, type="fan", edge.color=c("purple2"))
par(mar=c(1, 1, 1, 1), mfrow=c(1, 3))
plot.phylo(habitat.brains.tree.actinopteri.l, show.tip.label=FALSE)
plot.phylo(habitat.brains.tree.elasmobranchii.l, show.tip.label=FALSE, edge.color=c("green3"))
plot.phylo(habitat.brains.tree.petromyzonti.l, show.tip.label=FALSE, edge.color=c("purple2"))
par(mar=c(5, 5, 1, 1), mfrow=c(1, 1))

# investigating phyologenetic autocorrelation in absolute brain size, body weight, pelagic care and migratory
phylo.signal.formula.habitat.brains.brain<-BrainWeightMaxWeightIndWithBrain~Order/Family/Genus
phylo.signal.habitat.brains.brain.actinopteri<-correlogram.formula(phylo.signal.formula.habitat.brains.brain, data=habitat.brains.data.actinopteri)
phylo.signal.formula.habitat.brains.weight<-Weight~Order/Family/Genus
phylo.signal.habitat.brains.weight.actinopteri<-correlogram.formula(phylo.signal.formula.habitat.brains.weight, data=habitat.brains.data.actinopteri)
phylo.signal.formula.habitat.brains.migratory<-Migratory~Order/Family/Genus
phylo.signal.habitat.brains.migratory.actinopteri<-correlogram.formula(phylo.signal.formula.habitat.brains.migratory, data=habitat.brains.data.actinopteri)
phylo.signal.formula.habitat.brains.pelagic<-Pelagic~Order/Family/Genus
phylo.signal.habitat.brains.pelagic.actinopteri<-correlogram.formula(phylo.signal.formula.habitat.brains.pelagic, data=habitat.brains.data.actinopteri)
phylo.signal.formula.habitat.brains.brain<-BrainWeightMaxWeightIndWithBrain~Order/Family/Genus
phylo.signal.habitat.brains.brain.elasmobranchii<-correlogram.formula(phylo.signal.formula.habitat.brains.brain, data=habitat.brains.data.elasmobranchii)
phylo.signal.formula.habitat.brains.weight<-Weight~Order/Family/Genus
phylo.signal.habitat.brains.weight.elasmobranchii<-correlogram.formula(phylo.signal.formula.habitat.brains.weight, data=habitat.brains.data.elasmobranchii)
phylo.signal.formula.habitat.brains.migratory<-Migratory~Order/Family/Genus
phylo.signal.habitat.brains.migratory.elasmobranchii<-correlogram.formula(phylo.signal.formula.habitat.brains.migratory, data=habitat.brains.data.elasmobranchii)
phylo.signal.formula.habitat.brains.pelagic<-Pelagic~Order/Family/Genus
phylo.signal.habitat.brains.pelagic.elasmobranchii<-correlogram.formula(phylo.signal.formula.habitat.brains.pelagic, data=habitat.brains.data.elasmobranchii)
phylo.signal.formula.habitat.brains.brain<-BrainWeightMaxWeightIndWithBrain~Family/Genus
phylo.signal.habitat.brains.brain.petromyzonti<-correlogram.formula(phylo.signal.formula.habitat.brains.brain, data=habitat.brains.data.petromyzonti)
phylo.signal.formula.habitat.brains.weight<-Weight~Family/Genus
phylo.signal.habitat.brains.weight.petromyzonti<-correlogram.formula(phylo.signal.formula.habitat.brains.weight, data=habitat.brains.data.petromyzonti)
phylo.signal.formula.habitat.brains.migratory<-Migratory~Family/Genus
phylo.signal.habitat.brains.migratory.petromyzonti<-correlogram.formula(phylo.signal.formula.habitat.brains.migratory, data=habitat.brains.data.petromyzonti)
phylo.signal.formula.habitat.brains.pelagic<-Pelagic~Family/Genus
phylo.signal.habitat.brains.pelagic.petromyzonti<-correlogram.formula(phylo.signal.formula.habitat.brains.pelagic, data=habitat.brains.data.petromyzonti)

par(mfrow=c(4, 3), mar=c(2, 5, 1, 1))
plot(phylo.signal.habitat.brains.brain.actinopteri, ylim=c(0.3, 0.65))
legend("bottomleft", c("A"), bty="n")
plot(phylo.signal.habitat.brains.brain.elasmobranchii, ylim=c(-0.15, 0.55))
legend("bottomleft", c("B"), bty="n", text.col="green3")
plot(phylo.signal.habitat.brains.brain.petromyzonti, ylim=c(-0.2, 0.35))
legend("bottomleft", c("C"), bty="n", text.col="purple3")
plot(phylo.signal.habitat.brains.weight.actinopteri, ylim=c(-0.05, 0.2))
legend("bottomleft", c("D"), bty="n")
plot(phylo.signal.habitat.brains.weight.elasmobranchii, ylim=c(-0.05, 0.3))
legend("bottomleft", c("E"), bty="n", text.col="green3")
plot(phylo.signal.habitat.brains.weight.petromyzonti, ylim=c(-0.15, 0.4))
legend("bottomleft", c("F"), bty="n", text.col="purple3")
plot(phylo.signal.habitat.brains.migratory.actinopteri, ylim=c(0.25, 0.8))
legend("bottomleft", c("G"), bty="n")
plot(phylo.signal.habitat.brains.weight.elasmobranchii, col=c(0), lty=0, pch=0, xaxt="n", yaxt="n", ylab="", xlab="")
legend("bottomleft", c("H"), bty="n", text.col="green3")
plot(phylo.signal.habitat.brains.migratory.petromyzonti, ylim=c(-0.2,  0.1))
legend("bottomleft", c("I"), bty="n", text.col="purple3")
plot(phylo.signal.habitat.brains.pelagic.actinopteri, ylim=c(0.25, 0.7))
legend("bottomleft", c("J"), bty="n")
plot(phylo.signal.habitat.brains.pelagic.elasmobranchii, ylim=c(-0.1, 0.4))
legend("bottomleft", c("K"), bty="n", text.col="green3")
plot(phylo.signal.habitat.brains.weight.petromyzonti, col=c(0), lty=0, pch=0, xaxt="n", yaxt="n", ylab="", xlab="")
legend("bottomleft", c("L"), bty="n", text.col="purple3")
par(mfrow=c(1, 1), mar=c(5, 5, 1, 1))

# phylogenetic signal Actinopteri only
par(mfrow=c(4, 1), mar=c(2, 5, 1, 1))
plot(phylo.signal.habitat.brains.brain.actinopteri, ylim=c(0.3, 0.65))
legend("bottomleft", c("A"), bty="n")
plot(phylo.signal.habitat.brains.weight.actinopteri, ylim=c(-0.05, 0.2))
legend("bottomleft", c("B"), bty="n")
plot(phylo.signal.habitat.brains.migratory.actinopteri, ylim=c(0.25, 0.8))
legend("bottomleft", c("C"), bty="n")
plot(phylo.signal.habitat.brains.pelagic.actinopteri, ylim=c(0.25, 0.7))
legend("bottomleft", c("D"), bty="n")

# computing phylogenetically controlled generalised least squares models for Actinopteri - absolute values
habitat.brains.pGLS.brownian.actinopteri<-corBrownian(phy=habitat.brains.tree.actinopteri.l)
habitat.brains.pGLS.model.actinopteri<-gls(logBrain~Migratory+Pelagic+logWeight, data=habitat.brains.data.actinopteri, correlation=habitat.brains.pGLS.brownian.actinopteri)
# Table A5.3
summary(habitat.brains.pGLS.model.actinopteri)
# Figure A5.3
par(mfrow=c(1, 3), mar=c(5, 5, 1, 1))
plot(resid(habitat.brains.pGLS.model.actinopteri)~fitted(habitat.brains.pGLS.model.actinopteri), xlab="Fitted values", ylab="Residuals")
abline(h=0, lty=2)
hist(habitat.brains.pGLS.model.actinopteri$residuals, main="Habitat - Actinopteri", xlab="Residuals")
qqnorm(habitat.brains.pGLS.model.actinopteri$residuals, main="")
qqline(habitat.brains.pGLS.model.actinopteri$residuals)
habitat.brains.pGLS.model.actinopteri.residuals<-resid(habitat.brains.pGLS.model.actinopteri)
shapiro.test(habitat.brains.pGLS.model.actinopteri.residuals)

# computing phylogenetically controlled generalised least squares models for Actinopteri - Z values
habitat.brains.pGLS.brownian.actinopteri<-corBrownian(phy=habitat.brains.tree.actinopteri.l)
habitat.brains.pGLS.model.actinopteriZ<-gls(BrainZ~Migratory+Pelagic+WeightZ, data=habitat.brains.data.actinopteri, correlation=habitat.brains.pGLS.brownian.actinopteri)
# Table A5.4
summary(habitat.brains.pGLS.model.actinopteriZ)
# Figure A5.4
par(mfrow=c(1, 3), mar=c(5, 5, 1, 1))
plot(resid(habitat.brains.pGLS.model.actinopteriZ)~fitted(habitat.brains.pGLS.model.actinopteriZ), xlab="Fitted values", ylab="Residuals")
abline(h=0, lty=2)
hist(habitat.brains.pGLS.model.actinopteriZ$residuals, main="Habitat - Actinopteri Z", xlab="Residuals")
qqnorm(habitat.brains.pGLS.model.actinopteriZ$residuals, main="")
qqline(habitat.brains.pGLS.model.actinopteriZ$residuals)
habitat.brains.pGLS.model.actinopteriZ.residuals<-resid(habitat.brains.pGLS.model.actinopteriZ)
shapiro.test(habitat.brains.pGLS.model.actinopteriZ.residuals)

# proper plot
par(mar=c(5, 5, 1, 1), mfrow=c(1, 2))
boxplot(log(BrainWeightMaxWeightIndWithBrain)~Migratory*Class, data=subset(habitat.brains.data, habitat.brains.data$Class=="Actinopteri" | habitat.brains.data$Class=="Petromyzonti"), col=c(0, 0, "purple3", "purple3"), ylab="Brain size [ln(mg)]", xaxt="n", xlab="", cex.lab=1.5, ylim=c(-0.75, 9))
axis(side=1, at=c(1, 2, 3, 4), c("no", "yes", "no", "yes"))
axis(side=1, at=c(2.5), "Migration", mgp=c(3, 3, 1), cex.axis=1.5, tick=F)
legend("topright", pch=c(22, 22), pt.bg=c(0, "purple3"), c("Actinopteri", "Petromyzotni"))
legend("topleft", bty="n", c(""), cex=1.5, pch="A")
text(x=1, y=-0.5, c(length(subset(habitat.brains.data.actinopteri$Species, habitat.brains.data.actinopteri$Migratory=="0"))))
text(x=2, y=-0.5, c(length(subset(habitat.brains.data.actinopteri$Species, habitat.brains.data.actinopteri$Migratory=="1"))))
text(x=3, y=-0.5, c(length(subset(habitat.brains.data.petromyzonti$Species, habitat.brains.data.petromyzonti$Migratory=="0"))))
text(x=4, y=-0.5, c(length(subset(habitat.brains.data.petromyzonti$Species, habitat.brains.data.petromyzonti$Migratory=="1"))))

boxplot(log(BrainWeightMaxWeightIndWithBrain)~Pelagic*Class, data=subset(habitat.brains.data, habitat.brains.data$Class=="Actinopteri" | habitat.brains.data$Class=="Elasmobranchii"), col=c(0, 0, "green3", "green3"), ylab="", xaxt="n", xlab="", cex.lab=1.5, ylim=c(0, 13))
axis(side=1, at=c(1, 2, 3, 4), c("no", "yes", "no", "yes"))
axis(side=1, at=c(2.5), "Pelagic", mgp=c(3, 3, 1), cex.axis=1.5, tick=F)
legend("bottomright", pch=c(22, 22), pt.bg=c(0, "green3"), c("Actinopteri", "Elasmobranchii"))
legend("topleft", bty="n", c(""), cex=1.5, pch="B")

# proper plot - Actinopteri only
par(mar=c(5, 5, 1, 1), mfrow=c(1, 2))
boxplot(logBrain~Migratory, data=subset(habitat.brains.data, habitat.brains.data$Class=="Actinopteri"), col=c("gold1", "gold4"), ylab="Brain size [ln(mg)]", xaxt="n", xlab="", cex.lab=1.5, ylim=c(-0.75, 9))
axis(side=1, at=c(1, 2), c("no", "yes"))
axis(side=1, at=c(1.5), "Migration", mgp=c(3, 3, 1), cex.axis=1.5, tick=F)
legend("topleft", bty="n", c(""), cex=1.5, pch="A")
text(x=1, y=-0.5, c(length(subset(habitat.brains.data.actinopteri$Species, habitat.brains.data.actinopteri$Migratory=="0"))))
text(x=2, y=-0.5, c(length(subset(habitat.brains.data.actinopteri$Species, habitat.brains.data.actinopteri$Migratory=="1"))))

boxplot(logBrain~Pelagic, data=subset(habitat.brains.data, habitat.brains.data$Class=="Actinopteri"), col=c("springgreen1", "springgreen4"), ylab="", xaxt="n", xlab="", cex.lab=1.5, ylim=c(-0.75, 9))
axis(side=1, at=c(1, 2, 3, 4), c("no", "yes", "no", "yes"))
axis(side=1, at=c(1.5), "Pelagic", mgp=c(3, 3, 1), cex.axis=1.5, tick=F)
legend("topleft", bty="n", c(""), cex=1.5, pch="B")
text(x=1, y=-0.5, c(length(subset(habitat.brains.data.actinopteri$Species, habitat.brains.data.actinopteri$Pelagic=="0"))))
text(x=2, y=-0.5, c(length(subset(habitat.brains.data.actinopteri$Species, habitat.brains.data.actinopteri$Pelagic=="1"))))

# Figure 3 - habitat.brains. proper plot - Actinopteri only Z
par(mar=c(5, 5, 1, 1), mfrow=c(1, 2))
boxplot(BrainZ~Migratory, data=habitat.brains.data.actinopteri, col=c("gold1", "gold4"), ylab="Brain size [Z]", xaxt="n", xlab="", cex.lab=1.5, ylim=c(-3.5, 2.5))
axis(side=1, at=c(1, 2), c("no", "yes"))
axis(side=1, at=c(1.5), "Migration", mgp=c(3, 3, 1), cex.axis=1.5, tick=F)
legend("topleft", bty="n", c(""), cex=1.5, pch="A")
text(x=1, y=-3.5, c(length(subset(habitat.brains.data.actinopteri$Species, habitat.brains.data.actinopteri$Migratory=="0"))))
text(x=2, y=-3.5, c(length(subset(habitat.brains.data.actinopteri$Species, habitat.brains.data.actinopteri$Migratory=="1"))))

boxplot(BrainZ~Pelagic, data=habitat.brains.data.actinopteri, col=c("springgreen1", "springgreen4"), ylab="", xaxt="n", xlab="", cex.lab=1.5, ylim=c(-3.5, 2.5))
axis(side=1, at=c(1, 2, 3, 4), c("no", "yes", "no", "yes"))
axis(side=1, at=c(1.5), "Pelagic", mgp=c(3, 3, 1), cex.axis=1.5, tick=F)
legend("topleft", bty="n", c(""), cex=1.5, pch="B")
text(x=1, y=-3.5, c(length(subset(habitat.brains.data.actinopteri$Species, habitat.brains.data.actinopteri$Pelagic=="0"))))
text(x=2, y=-3.5, c(length(subset(habitat.brains.data.actinopteri$Species, habitat.brains.data.actinopteri$Pelagic=="1"))))


# coloured trees
habitat.brains.tree.actinopteri.l.migrating<-which.edge(habitat.brains.tree.actinopteri.l, subset(habitat.brains.data.actinopteri$Species, habitat.brains.data.actinopteri$Migratory=="1"))
habitat.brains.tree.actinopteri.l.monogamycolours<- rep("gold1", dim(habitat.brains.tree.actinopteri.l$edge)[1]) 
habitat.brains.tree.actinopteri.l.monogamycolours[habitat.brains.tree.actinopteri.l.migrating] <- "gold4"
par(mar=c(0 ,0, 0, 0))
plot(habitat.brains.tree.actinopteri.l, lwd=3, edge.color=habitat.brains.tree.actinopteri.l.monogamycolours, show.tip.label=FALSE)

habitat.brains.tree.actinopteri.l.pelagic<-which.edge(habitat.brains.tree.actinopteri.l, subset(habitat.brains.data.actinopteri$Species, habitat.brains.data.actinopteri$Pelagic=="1"))
habitat.brains.tree.actinopteri.l.carecolours<- rep("springgreen1", dim(habitat.brains.tree.actinopteri.l$edge)[1]) 
habitat.brains.tree.actinopteri.l.carecolours[habitat.brains.tree.actinopteri.l.pelagic] <- "springgreen4"
par(mar=c(0 ,0, 0, 0))
plot(habitat.brains.tree.actinopteri.l, lwd=3, edge.color=habitat.brains.tree.actinopteri.l.carecolours, show.tip.label=FALSE, direction = "leftwards")

# Figure 1: Publication bias----
par(mar=c(5, 5, 1, 1), mfrow=c(1, 1))
plot(cum.Nr.publ~Year, data=publications, col=0, xlab="Year of publication", ylab="Cumulative number of articles", cex.lab=1.5)
lines(cum.Nr.publ~Year, data=subset(publications, publications$Taxon=="fish"), lty=1, lwd=3, col="deepskyblue3")
points(cum.Nr.publ~Year, data=subset(publications, publications$Taxon=="fish"), pch=1, col="deepskyblue3")
lines(cum.Nr.publ~Year, data=subset(publications, publications$Taxon=="mammal"), lty=1, lwd=3, col="chocolate4")
points(cum.Nr.publ~Year, data=subset(publications, publications$Taxon=="mammal"), pch=2, col="chocolate4")
lines(cum.Nr.publ~Year, data=subset(publications, publications$Taxon=="bird"), lty=1, lwd=3, col="chartreuse4")
points(cum.Nr.publ~Year, data=subset(publications, publications$Taxon=="bird"), pch=0, col="chartreuse4")
legend("topleft", c("Birds", "Mammals", "Fishes"), lwd=c(2, 2, 2), pch=c(0, 2, 1), col=c("chartreuse4","chocolate4", "deepskyblue3"))
