rm(list = ls())
dat <- read.csv(file.choose(), sep = "\\t") #load in ebd_data
head(dat)

#######################
# SUBSETTING THE DATA # 
#######################
# subsetting to only complete checklists and removing incidental records
dat <- droplevels(subset(dat, ALL.SPECIES.REPORTED == 1))
dat <- droplevels(subset(dat, PROTOCOL.TYPE != "Incidental"))

# stripping away unnecessary columns
kmng <- dat[, -which(names(dat) %in% c("GLOBAL.UNIQUE.IDENTIFIER", "LAST.EDITED.DATE", "TAXONOMIC.ORDER", "CATEGORY",
                                       "SUBSPECIES.COMMON.NAME", "SUBSPECIES.SCIENTIFIC.NAME", "BREEDING.BIRD.ATLAS.CODE",
                                       "AGE.SEX", "COUNTRY", "COUNTRY.CODE", "STATE", "STATE.CODE", "COUNTY", "COUNTY.CODE",
                                       "IBA.CODE", "BCR.CODE", "USFWS.CODE", "ATLAS.BLOCK", "PROTOCOL.CODE", "PROJECT.CODE",
                                       "EFFORT.AREA.HA", "GROUP.IDENTIFIER", "HAS.MEDIA", "APPROVED", "REVIEWED", "REASON", 
                                       "TRIP.COMMENTS", "SPECIES.COMMENTS", "BREEDING.BIRD.ATLAS.CATEGORY"))]
kmng <- droplevels(kmng)
head(kmng)                                     
sort(unique(kmng$LOCALITY))

hotspots <- c(# 'official eBird' titles for Eaglenest hotspots
  "Eaglenest WLS--Khellong", "Eaglenest WLS--Sessni", "Eaglenest WLS--Bompu Camp",
  "Eaglenest WLS--Eaglenest Pass", "Eaglenest WLS--Lama Camp",
  
  # Khellong cluster
  "Khellong",
  
  # Sessni cluster: 1250m
  "Sessni Camp", "Eaglenest Wildlife Sanctuary--Sessni",
  "ses-1.1250m", "ses-1.1300m", "ses-1.1350m", "ses-2.1150m", "ses-2.1200m", # transects by Umesh Srinivasan (Srinivasan et al. 2018 Proceedings B, 285, 20172593)
  
  # Bompu cluster
  "Eaglenest WLS, Bompu Camp", "Eaglenest WS---Bompu", "Eaglenest WS--Bompu Camp", 
  "Bhompu", "Bomphu", "bompu camp", "Bompu camp", "Bompu Camp", "Bompu Camp, Eaglenest WLS",
  "Bompu Camp, Eagles Nest WS", "Bompu Camp, West Kameng IN-Arunachal Pradesh", 
  "Bompu Eglenest Sanctuary.", "Bompu, Eaglenest WLS",
  "EWLS-Plot B", "EWLS-Plot C", "EWLS-Plot D", "EWLS-Plot J", # point counts by Binod Borah at Bompu camp elevation (Borah et al. 2018 J Applied Ecology, 55, 1637-1646)
  "bom-1.1950m", "bom-1.2000m", "bom-1.2050m","bom-2.1850m", "bom-2.1900m", # transects by Umesh Srinivasan (Srinivasan et al. 2018 Proceedings B, 285, 20172593)
  
  # Lama Camp cluster
  "lama camp", "Lama Camp", "Lama Camp (27.1525, 92.4675)", "Lama Camp (27.1569, 92.4606)", 
  "Lama Camp (27.1569, 92.4608)", "Lama Camp (27.1586, 92.4603)", "Lama Camp (27.1739, 92.4622)",
  "Lama Camp area", "Lama Camp surroundings", "Lama Camp, Eaglenest Wildlife Sanctuary", 
  "Lama Camp, Eaglenest Wildlife Sanctuary: Bugun Point", "Lama Camp, Eaglenest WLS", 
  "Lama Camp, Eaglenest WS", "Lama Camp, Eagles Nest WS","Eaglenest WS--Lama Camp",
  "lam-1.2250m", "lam-1.2300m", "lam-2.2350m", "lam-2.2400m", "lam-2.2450m", # transects by Umesh Srinivasan (Srinivasan et al. 2018 Proceedings B, 285, 20172593)
  
  # Eaglenest Pass cluster: 2780m
  "Eaglenest--Eaglenest Pass @ 27.1252x92.4559", "Eaglenest Pass", "Eaglenest Pass, Eaglenest WLS")

# subsetting by hotspots
ewls.htspt <- droplevels(subset(kmng, LOCALITY %in% hotspots))

# locations names in the data
table(ewls.htspt$LOCALITY)

# reclassifying location names into 5 consistently named locations
ewls.htspt$LOCALITY <- as.character(ewls.htspt$LOCALITY)
ewls.htspt$LOCALITY[which(ewls.htspt$LOCALITY %in% c("Bhompu", "Bomphu", "Bompu camp",
                                                     "Bompu Camp", "Bompu Camp. Eaglenest WLS",
                                                     "Bompu Camp, West Kameng IN-Arunachal Pradesh",
                                                     "Bompu Eglenest Sanctuary.", "Eaglenest WLS--Bompu Camp",
                                                     "Eaglenest WLS, Bompu Camp", "EWLS-Plot B", "EWLS-Plot C",
                                                     "EWLS-Plot D", "EWLS-Plot J", "Bompu Camp, Eaglenest WLS",
                                                     "bom-1.1950m", "bom-1.2000m", "bom-1.2050m","bom-2.1850m", "bom-2.1900m"))] <- "Bompu"
ewls.htspt$LOCALITY[which(ewls.htspt$LOCALITY %in% c("Eaglenest Pass", "Eaglenest--Eaglenest Pass @ 27.1252x92.4559",
                                                     "Eaglenest WLS--Eaglenest Pass"))] <- "Pass"
ewls.htspt$LOCALITY[which(ewls.htspt$LOCALITY %in% c("Eaglenest WLS--Lama Camp", "lama camp", "Lama Camp",
                                                     "Lama Camp area", "Lama Camp, Eaglenest Wildlife Sanctuary",
                                                     "Lama Camp, Eaglenest Wildlife Sanctuary: Bugun Point",
                                                     "Lama Camp, Eaglenest WLS", "Lama Camp, Eaglenest WS",
                                                     "lam-1.2250m", "lam-1.2300m", "lam-2.2350m", "lam-2.2400m", "lam-2.2450m"))] <- "Lama"
ewls.htspt$LOCALITY[which(ewls.htspt$LOCALITY %in% c("Eaglenest Wildlife Sanctuary--Sessni", "Eaglenest WLS--Sessni",
                                                     "ses-1.1250m", "ses-1.1300m", "ses-1.1350m", "ses-2.1150m", "ses-2.1200m"))] <- "Sessni"
ewls.htspt$LOCALITY[which(ewls.htspt$LOCALITY %in% c("Khellong", "Eaglenest WLS--Khellong"))] <- "Khellong"

sum(table(ewls.htspt$LOCALITY))

# subsetting to only breeding season
yy.mm.dd <- unlist(strsplit(as.character(ewls.htspt$OBSERVATION.DATE), "-"))
ewls.htspt$year <- as.numeric(yy.mm.dd[seq(1,length(yy.mm.dd),by=3)])
ewls.htspt$month <- as.numeric(yy.mm.dd[seq(2,length(yy.mm.dd),by=3)])
ewls.htspt <- droplevels(subset(ewls.htspt, month %in% 3:6))
#ewls.htspt <- droplevels(subset(ewls.htspt, month %in% 4:6)) 
#if march data is included, then keep as is: if excluded, then comment out line 85, and run line 86

# subsetting to checklists that are either stationary or have low effort in terms of distance
# stationary and traveling checklists:
dist.threshold <- 1 # threshold for inclusion of traveling checklists in the data (km)
ewls.htspt$EFFORT.DISTANCE.KM[which(is.na(ewls.htspt$EFFORT.DISTANCE.KM))] <- 0
ewls.summer <- subset(ewls.htspt, EFFORT.DISTANCE.KM <= dist.threshold)

##############################
# preliminary numbers/graphs #
##############################
length(unique(ewls.summer$SAMPLING.EVENT.IDENTIFIER)) # number of checklists
colSums(table(ewls.summer$SAMPLING.EVENT.IDENTIFIER, ewls.summer$LOCALITY) != 0) # number of checklists per location

# too few lists from Khellong, Sessni and Eaglenet Pass; removing these locations from the data
ewls.final <- droplevels(subset(ewls.summer, LOCALITY %in% c("Bompu", "Lama")))
length(unique(ewls.final$SAMPLING.EVENT.IDENTIFIER))

#a subset of ewls.final may be seen in the Table S2.

# temporal distribution of records
hist(ewls.final$year, main = "", br = 30,
     xlab = "Year", ylab = "Number of records") # clear separation of record into two year clusters:
# 1. 2006 to 2010
# 2. 2015 to 2019
# no records from the years 2011, 2012, 2013 and 2014

# splitting years into "early" and "later" periods
ewls.final$period <- as.character(nrow(ewls.final))
ewls.final$period[which(ewls.final$year %in% 2006:2010)] <- "early"
ewls.final$period[which(ewls.final$year %in% 2015:2019)] <- "later"

# checklist properties; using a loop; "duplicated()" is messing up the checklist properties
list.no <- unique(ewls.final$SAMPLING.EVENT.IDENTIFIER)
time <- loc <- character(length(list.no)); no.sp <- numeric(length(list.no))
for(i in 1:length(list.no)){
  temp <- droplevels(subset(ewls.final, SAMPLING.EVENT.IDENTIFIER == list.no[i]))
  time[i] <- temp$period[1]; loc[i] <- temp$LOCALITY[1]; no.sp[i] <- nrow(temp)
}
list.props <- data.frame(list.no, no.sp, loc, time)

# selecting lists that are longer than median list length
selected.lists <- droplevels(subset(list.props, no.sp >= median(list.props$no.sp)))
ewls.selected <- droplevels(subset(ewls.final, SAMPLING.EVENT.IDENTIFIER %in% selected.lists$list.no))
nrow(ewls.selected)

par(mfrow = c(2,1)); par(mai = c(0.6, 0.6, 0.1, 0.1)); par(oma = c(0,2,0,0))
boxplot(no.sp ~ time, data = selected.lists[which(selected.lists$loc == "Bompu"),], las = 1,
        names = c("2006 to 2010", "2015 to 2019"),
        ylab = "Number of species per checklist")
text(.5, 68, labels = "a", font = 2, cex = 1.5)
boxplot(no.sp ~ time, data = selected.lists[which(selected.lists$loc == "Lama"),], las = 1,
        names = c("2006 to 2010", "2015 to 2019"),
        ylab = "Number of species per checklist")
text(0.5, 59, labels = "b", cex = 1.5, font = 2)
mtext("Checklist length (number of species per checklist)", 2, line = 3, at = 70)

# earlier checklists have more species than later checklists

##########################
# SPECIES-LEVEL ANALYSES #
##########################
# logistic regressions to
# estimate proportion of checklists in which a species occurred over the two time periods, while
# accounting for differences in list length

# species selection
threshold <- 10 # selecting species that have been reported a reasonable number of times (no rare species) from each location
bom <- droplevels(subset(ewls.selected, LOCALITY == "Bompu"))
lam <- droplevels(subset(ewls.selected, LOCALITY == "Lama"))

sp.bom <- subset(data.frame(table(bom$SCIENTIFIC.NAME)), Freq >= threshold)[1]
sp.lam <- subset(data.frame(table(lam$SCIENTIFIC.NAME)), Freq >= threshold)[1]

species.both <- sort(union(sp.bom$Var1, sp.lam$Var1))
eagles <- c("Ictinaetus malaiensis", "Nisaetus nipalensis") #remove eagles: see paper
species.both <- setdiff(species.both, eagles)
species <- species.both

# 'newdata' for predictions from logistic regression
new <- data.frame(expand.grid(40, 
                              c("Bompu", "Lama"), 
                              c("early", "later"))) # estimate for checklist with 40 species
dimnames(new) <- list(1:nrow(new), c("no.sp", "locn", "time"))

# empty matrices for storing proportions
bompu <- lama <- matrix(0, length(species), 2)
dimnames(bompu) <- dimnames(lama) <- list(species, c("p.early", "p.later"))
n.bom <- n.lam <- numeric(length(species))

for(i in 1:length(species)){
  temp <- droplevels(subset(ewls.selected, SCIENTIFIC.NAME == species.both[i]))
  
  # response variable (binomial: present in a list or absent in a list)
  resp <- numeric(length(selected.lists$list.no))
  resp[which(selected.lists$list.no %in% temp$SAMPLING.EVENT.IDENTIFIER)] <- 1
  
  locn <- selected.lists$loc; time <- selected.lists$time; no.sp <- selected.lists$no.sp
  
  n.bom[i] <- length(which(temp$LOCALITY == "Bompu"))
  n.lam[i] <- length(which(temp$LOCALITY == "Lama"))
  
  mod.data <- data.frame(resp, locn, time, no.sp)
  mod.glm <- glm(resp ~ no.sp + locn + time + locn:time, family = "binomial")
  
  preds <- round(predict(mod.glm, newdata = new, type = "response"), 2)
  data.frame(new, preds)
  
  bompu[i,1] <- preds[1]; bompu[i,2] <- preds[3]
  lama[i,1] <- preds[2]; lama[i,2] <- preds[4]
}
bom.diff <- bompu[,2] - bompu[,1]; lam.diff <- ses.diff <- lama[,2] - lama[,1]

bbb <- cbind(bom.diff, species)
ccc <- cbind(lam.diff, species)

# although model estimates some near-zero probability of occurrence for species never recorded at a site; 
# removing these estimates (replacing with NA); then remove the eagles
bom.diff[which(!species %in% sp.bom$Var1)] <- NA
lam.diff[which(!species %in% sp.lam$Var1)] <- NA

par(mfrow = c(2,1)); par(oma = c(2, 2, 0, 0))
hist(bom.diff, col = "grey", las = 1, ylim = c(0,10), # xlim = c(-0.4,0.4),
     br = seq(-1, 1, by = 0.1), 
     xlab = "", ylab = "", main = "")
text(-0.95, 8, labels = "a", font = 2, cex = 1.5)
hist(lam.diff, col = "grey", las = 1, ylim = c(0,12), # xlim = c(-0.4,0.4), 
     br = seq(-1, 1, by = 0.1), 
     xlab = "Difference in probability of occurrence", ylab = "", main = "")
text(-0.95, 9, labels = "b", font = 2, cex = 1.5)
mtext("Number of species", 2, line = 3, at = 12)
mtext("Difference in probability of occurrence", 1, line = 3)

prop.data <- data.frame(species, bom.diff, lam.diff, n.bom, n.lam)
dimnames(prop.data) <- list(1:nrow(prop.data), c("species", "bom.diff", "lam.diff", "n.bom", "n.lam"))

# relationship between elevational range and change in probability of reporting
elev.data <- read.csv(file.choose()) #load in species_elevation_data
elev.data <- droplevels(subset(elev.data, species %in% species.both)) #elev.data initially has 94 species, because of Phylloscopus intermedius
elev.data <- elev.data[order(elev.data$species),]
rownames(elev.data)<-NULL
attach(elev.data)

range.mid <- apply(cbind(E.summer.low, E.summer.high), 1, mean)

bom.elev.diff <- range.mid - 1950
lam.elev.diff <- range.mid - 2350

hist(bom.elev.diff)
hist(lam.elev.diff)

prop.data <- data.frame(prop.data, bom.elev.diff, lam.elev.diff, mass, E.summer.low, E.summer.high)

# excluding species that are likely to be seasonal migrants, especially high elevation breeders likely to be on passage
bom.high.sp <- as.character(prop.data$species[which(prop.data$E.summer.low > 1950)])
bom.low.sp <- as.character(prop.data$species[which(prop.data$E.summer.high < 1950)])
prop.data$bom.diff[which(prop.data$species %in% c(bom.high.sp, bom.low.sp))] <- NA
bom.dat <- droplevels(subset(prop.data, !is.na(bom.diff)))

lam.high.sp <- as.character(prop.data$species[which(prop.data$E.summer.low > 2350)])
lam.low.sp <- as.character(prop.data$species[which(prop.data$E.summer.high < 2350)])
prop.data$lam.diff[which(prop.data$species %in% c(lam.high.sp, lam.low.sp))] <- NA
lam.dat <- droplevels(subset(prop.data, !is.na(lam.diff)))

common.sp <- intersect(lam.dat$species, bom.dat$species)
com <- prop.data[which(prop.data$species %in% common.sp),]

library(DescTools)
par(mfrow = c(1,2)); par(mai = c(0.7, 0.7, 0.1, 0.1)); par(oma = c(1,1,0,0))
plot(bom.diff ~ bom.elev.diff, data = bom.dat, las = 1, xlab = "", ylab = "",
     pch = ifelse(bom.dat$species %in% common.sp, 1, 2), cex = 1)
text(-800, 0.8, labels = "a", font = 2, cex = 1.5)
abline(v = 0, lty = 2); abline(h = 0, lty = 2)
mod.bom <- lm(bom.diff ~ bom.elev.diff, data = bom.dat)
lines.lm(mod.bom, col = "black", conf.level = 0.95)

# text(com$bom.elev.diff, com$bom.diff, labels = as.character(com$species), cex = 0.4)

plot(lam.diff ~ lam.elev.diff, data = lam.dat, las = 1, xlab = "", ylab = "",
     pch = ifelse(lam.dat$species %in% common.sp, 1, 2), cex = 1)
text(-800, 0.575, labels = "b", font = 2, cex = 1.5)
abline(v = 0, lty = 2); abline(h = 0, lty = 2)
mod.lam <- lm(lam.diff ~ lam.elev.diff, data = lam.dat)
lines.lm(mod.lam, col = "black", conf.level = 0.95)

summary(mod.bom)
summary(mod.lam)


# text(com$lam.elev.diff, com$lam.diff, labels = as.character(com$species), cex = 0.4)
par(xpd = F)
mtext("Difference in elevational range mid-point and eBird hotspot location (m)", 1, at = -1100, line = 2.5, cex = 0.8)
mtext("Change in probability of occurrence", 2, line = 24.5, cex = 0.8)

summary(mod.bom)
summary(mod.lam)

# testing for phylogenetic autocorrelation
# tree
tree <- ape::read.nexus(file.choose()) #load in the .nex file of 93 species
plot(tree[[1]])

# need consistency in the scientific names between the tree and eBird
tree.species <- sort(tree$tree_2041$tip.label)
gen.sp <- unlist(strsplit(tree.species, "_"))
genus <- gen.sp[seq(1, length(gen.sp), by = 2)]
spces <- gen.sp[seq(2, length(gen.sp), by = 2)]

tree.sp <- paste(genus, spces, sep = " ")
species[which(is.na(match(species, tree.sp)))]

species.renamed <- species
species.renamed[which(species == "Actinodura strigula")] <- "Minla strigula"
species.renamed[which(species == "Caprimulgus jotaka")] <- "Caprimulgus indicus"
species.renamed[which(species == "Carpodacus sipahi")] <- "Haematospiza sipahi"
species.renamed[which(species == "Cettia castaneocoronata")] <- "Tesia castaneocoronata"
species.renamed[which(species == "Chelidorhynx hypoxanthus")] <- "Chelidorhynx hypoxantha"
species.renamed[which(species == "Cyanoderma chrysaeum")] <- "Stachyris chrysaea"
species.renamed[which(species == "Cyanoderma ruficeps")] <- "Stachyris ruficeps"
species.renamed[which(species == "Ficedula hodgsoni")] <- "Muscicapella hodgsoni"
species.renamed[which(species == "Grammatoptila striata")] <- "Garrulax striatus"
species.renamed[which(species == "Hierococcyx sparverioides")] <- "Cuculus sparverioides"
species.renamed[which(species == "Horornis fortipes")] <- "Cettia fortipes"
species.renamed[which(species == "Ianthocincla caerulata")] <- "Garrulax caerulatus"
species.renamed[which(species == "Lalage melaschistos")] <- "Coracina melaschistos"
species.renamed[which(species == "Lioparus chrysotis")] <- "Alcippe chrysotis"
species.renamed[which(species == "Machlolophus spilonotus")] <- "Parus spilonotus"
species.renamed[which(species == "Phylloscopus castaniceps")] <- "Seicercus castaniceps"
species.renamed[which(species == "Phylloscopus poliogenys")] <- "Seicercus poliogenys"
species.renamed[which(species == "Psilopogon franklinii")] <- "Megalaima franklinii"
species.renamed[which(species == "Psilopogon virens")] <- "Megalaima virens"
species.renamed[which(species == "Psittiparus ruficeps")] <- "Paradoxornis ruficeps"
species.renamed[which(species == "Pteruthius aeralatus")] <- "Pteruthius flaviscapis"
species.renamed[which(species == "Schoeniparus castaneceps")] <- "Alcippe castaneceps"
species.renamed[which(species == "Schoeniparus cinereus")] <- "Alcippe cinereus"
species.renamed[which(species == "Tarsiger rufilatus")] <- "Tarsiger cyanurus"
species.renamed[which(species == "Trochalopteron affine")] <- "Garrulax affinis"
species.renamed[which(species == "Trochalopteron erythrocephalum")] <- "Garrulax erythrocephalus"
species.renamed[which(species == "Trochalopteron imbricatum")] <- "Garrulax imbricatus"

new.gen.sp <- unlist(strsplit(species.renamed, " "))
new.gen <- new.gen.sp[seq(1, length(new.gen.sp), by = 2)]
new.sps <- new.gen.sp[seq(2, length(new.gen.sp), by = 2)]

new.names <- paste(new.gen, new.sps, sep = "_")
new.names[1] <- "[Redacted_due_to_Vulnerability_3]"
new.names[2] <- "[Redacted_due_to_Vulnerability_9]"
new.names <- new.names[-c(3:5)]  
phyl.dat <- data.frame(new.names, bom.diff, lam.diff)

library(ape)
# bompu differences
moran.scores <- moran.pvals <- numeric(100)
for (i in 1:100) {
  treei <- tree[[i]]
  w <- 1/cophenetic(treei)
  diag(w) <- 0
  moran.scores[i] <- Moran.I(phyl.dat$bom.diff, w, na.rm = T)$observed
}
round(quantile(moran.scores, p = c(0.025, 0.975)), 2) # includes zero; no significant phylogenetic structure

# lama differences
moran.scores <- moran.pvals <- numeric(100)
for (i in 1:100) {
  treei <- tree[[i]]
  w <- 1/cophenetic(treei)
  diag(w) <- 0
  moran.scores[i] <- Moran.I(phyl.dat$lam.diff, w, na.rm = T)$observed
}
round(quantile(moran.scores, p = c(0.025, 0.975)), 2) # includes zero; no significant phylogenetic structure

# The 95% CIs for phylogenetic Moran's I for 100 trees for the difference in proportion of checklists in 
# which a species was reported ranged from were [-0.02, 0.01] for Bompu and [-0.04, 0.01] for Lama Camp,
# indicating that species'responses were not phylogenetically strucutred.