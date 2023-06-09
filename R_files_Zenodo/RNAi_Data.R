#2020_03_27 CR Fisher
require(dplyr)
require(tidyr)
##### Set wd to where you are keeping the datatables
setwd("C:/Users/cruth/Documents/ONCOPELTUSTEMP/DATA/")
#setwd("C:/Users/cruth/Google Drive/Treehoppers/MediaPosterPresentation/PublicationDraftFolder/OncopeltusRNAi_Manuscript/datatables/")
#####
options(stringsAsFactors = FALSE)
df <- read.csv("CombinedTable_ForR_CRF.txt", header=TRUE, sep="\\t")
dim(df)
colnames(df)
str(df)
## All variables except for Experiment are character, which is fine. 
## We can co-erce to number when we need to.

## NAs mean that the trait wasn't on the scoring sheet when this
## experiment was scored. NS means that the trait was on the sheet,
## but not scoreable for this individual. 

## Let's look at the How.Informative field values
df.inf <- arrange(df, How.Informative) %>% 
  select(Experiment, DPI, Individual, Gene.Control, Stage.Scored, How.Informative) %>%
  filter(How.Informative < 5)

### 169 specimens have less than or equal to 5 informative fields

count(df.inf, Gene.Control)

# # A tibble: 9 x 2
# Gene        n
# <chr>   <int>
#   1 control    49
# 2 hth        12
# 3 iro        68
# 4 mir         2
# 5 nub         3
# 6 srf         6
# 7 tio         1
# 8 tup        10
# 9 ush        18

## What about the ones that are just zeroes? 

df.inf <- filter(df.inf, How.Informative == 0)
dim(df.inf)
count(df.inf, Gene.Control)
# [1] 153   6
# > count(df.inf, Gene.Control)
# # A tibble: 8 x 2
# Gene.Control     n
# <chr>        <int>
#   1 control         49
# 2 hth             11
# 3 iro             68
# 4 mir              1
# 5 nub              2
# 6 srf              6
# 7 tup             10
# 8 ush              6
NumUnInformative <- count(df.inf, Gene.Control)

write.csv(NumUnInformative, "Number_Uninformative_by_Gene.txt")

### Let's make a table of just the ones that are uninformative, 
### In case we want to look at those specimens again later. 
df.inf <- filter(df, How.Informative==0)
write.table(df.inf, "Uninformative_Table.tab", sep="\\t")
getwd()

df.inf <- filter(df, !(Stage.Scored=="adult")) %>%
  filter(!(Stage.Scored=="L5-surgical")) %>% 
  filter(!(Stage.Scored=="adult-stuck.in.eclosion")) %>%
  filter(How.Informative==0)

### We'll keep a filtered file of just the uninformative nymphs, too, 
### In case we want to go look at those specimens later.
nymphdata <- df.inf
count(df.inf, Gene.Control)
write.table(nymphdata, "UninformativeNymphs.tab", sep="\\t")

####################################################
#
#     Data cleanup for publication
#
###################################################
##### Let's only keep the informative samples, then:
df.inf <- arrange(df, Gene.Control) %>% 
  filter(How.Informative > 0)

##### And actually let's please filter down to just adults. I've already 
##### manually curated this, calling it "adult" if things were able to be scored
## All I really want here is the population sizes of each group, so:
Exp.Groups <- df.inf %>% filter(Stage.Scored=="adult") %>%
  count(Gene.Control)

# Let's keep a table of our experimental group sizes -- we'll want that for 
# calculating penetrance
write.table(Exp.Groups, "GroupSizesOfScorableAdults.txt", sep="\\t", quote=FALSE)

str(df.inf)

#### And now we need to try to collapse down some of the data fields to
#### understand the results in terms of Our major phenotype categories: Collar,
#### posterior pronotum, scutellum, wings, pleural margins, supracoxal lobes,
#### and propleuron/pronotum junction

## Collar: column names = Collar
## Pronotum: Pronotum, D.T1.Posterior.margin.central, D.T1.Posterior.margin.corners
## Scutellum: Scutellar.lobe, Mesonotum
## Forewings: L.FW.hemelytra, L.FW.membranous, R.FW.hemelytra, R.FW.membranous
## Wing hinges: Wing.hinges, Wing.sclerites, T2.FW.L.Hinge.Area, T2.FW.L.Axillary.sclerites, T3.HW.L.Hinge.Area, T3.HW.L.Axillary.sclerites, T2.FW.R.Hinge.Area, T2.FW.R.Axillary.sclerites, T3.HW.R.Hinge.Area, T3.HW.R.Axillary.sclerites 
## Hindwings: L.HW, T3.HW.L.Wing.blade, R.HW, T3.HW.R.Wing.blade
## pleural margins: Pos.Pl1, T1.LL.posterior.pleural.margin, T1.RL.posterior.pleural.margin, Pos.Pl2, T2.LL.posterior.pleural.margin, T2.RL.posterior.pleural.margin, Pos.Pl.3, T3.LL.post.pleural.margin, T3.RL.post.pleural.margin, T2.L.wing.base, T2.R.wing.base, T3.LL.wing.base, T3.RL.wing.base   
## Supracoxal lobes:Pl.lobes1, T1.LL.supracoxal.lobes, T1.RL.supracoxal.lobes, Pl.lobes2, T2.LL.supracoxal.lobes, T2.RL.supracoxal.lobes
## Pronotal/propleuron junction: Pn.Pl.Junction, T1.LL.pronotum.pleural.junction, T1.RL.pronotum.pleural.junction

#scutellum.fields <- c("Scutellar.lobe", "Mesonotum")
#collar.fields <- c("Collar")
#pronotum.fields <- c("Pronotum", "D.T1.Posterior.margin.central", "D.T1.Posterior.margin.corners")
#forewings.fields <- c("L.FW.hemelytra", "L.FW.membranous", "R.FW.hemelytra", "R.FW.membranous")
#hindwings.fields <- c("L.HW", "T3.HW.L.Wing.blade", "R.HW", "T3.HW.R.Wing.blade")
# winghinges.fields <- c("Wing.hinges", "Wing.sclerites", "T2.FW.L.Hinge.Area", "T2.FW.L.Axillary.sclerites", 
#                        "T3.HW.L.Hinge.Area", "T3.HW.L.Axillary.sclerites", "T2.FW.R.Hinge.Area", "T2.FW.R.Axillary.sclerites",
#                        "T3.HW.R.Hinge.Area", "T3.HW.R.Axillary.sclerites")
# pleuralmargins.fields <- c("Pos.Pl1", "T1.LL.posterior.pleural.margin", "T1.RL.posterior.pleural.margin", "Pos.Pl2", 
#                            "T2.LL.posterior.pleural.margin", "T2.RL.posterior.pleural.margin", "Pos.Pl.3", 
#                            "T3.LL.post.pleural.margin", "T3.RL.post.pleural.margin", "T2.L.wing.base", "T2.R.wing.base", 
#                            "T3.LL.wing.base", "T3.RL.wing.base")
#supracoxallobes.fields <- c(Pl.lobes1, T1.LL.supracoxal.lobes, T1.RL.supracoxal.lobes, Pl.lobes2, T2.LL.supracoxal.lobes,
                           # T2.RL.supracoxal.lobes)
#propleuron.pronotum.junction.fields <- c("Pn.Pl.Junction", "T1.LL.pronotum.pleural.junction", "T1.RL.pronotum.pleural.junction")

#/--------------------------------------------\\#
#                                              #
# The rest of this script collapses fields     #
# that are synonymous. Duplicative/synonymous  #
# fields in the scoring sheet are the result of#
# multiple scorers working independently from  #
# each other over the course of several months.#
#                                              #
#\\--------------------------------------------/#

###################### SUPRACOXAL LOBES #######################################

df.sclobes <- df.inf %>% select(Experiment, DPI, Individual, Gene.Control, Pl.lobes1, T1.LL.supracoxal.lobes, T1.RL.supracoxal.lobes, Pl.lobes2, T2.LL.supracoxal.lobes,
                                T2.RL.supracoxal.lobes)

## I need to coerce all of these to integers...
df.sclobes$Pl.lobes1 <- as.integer(df.sclobes$Pl.lobes1) # NAs will be introduced by coercion and that's FINE
df.sclobes$T1.LL.supracoxal.lobes <- as.integer(df.sclobes$T1.LL.supracoxal.lobes)
df.sclobes$T1.RL.supracoxal.lobes <- as.integer(df.sclobes$T1.RL.supracoxal.lobes)
df.sclobes$Pl.lobes2 <- as.integer(df.sclobes$Pl.lobes2)
df.sclobes$T2.LL.supracoxal.lobes <- as.integer(df.sclobes$T2.LL.supracoxal.lobes)
df.sclobes$T2.RL.supracoxal.lobes <- as.integer(df.sclobes$T2.RL.supracoxal.lobes)

str(df.sclobes)

## Make a matrix that's just the data fields
dat.sclobes <- select(df.sclobes, Pl.lobes2, Pl.lobes1, T1.RL.supracoxal.lobes, T1.LL.supracoxal.lobes, T2.RL.supracoxal.lobes, T2.LL.supracoxal.lobes)
### and we'll add a column called "aggregate" to sum across the rows, removing NAs
dat.sclobes <- mutate(dat.sclobes, aggregate=apply(dat.sclobes, 1, sum, na.rm=TRUE))
### Now we'll put the dataframe back together
df.sclobes <- data.frame(c(df.sclobes[,c(1:4)], dat.sclobes))

###################################################
#        Synonymize duplicative phenotypes        #
###################################################
### In order to just pick the ones that had supracoxal lobe phenotypes recorded, 
### We just need to filter to the ones where aggregate > 0
df.sclobes <- filter(df.sclobes, aggregate>0)
### which turns out to be 228 of 949

########### Add to Frequency Table ##########
ScLobes.Phenotypes <- count(df.sclobes, Gene.Control) 

FreqTable <- full_join(Exp.Groups, ScLobes.Phenotypes, by="Gene.Control")
colnames(FreqTable) <- c("Gene.Control", "Total", "ScLobes")

write.table(FreqTable, "PhenotypeFrequencyTable.txt", sep="\\t", quote=FALSE)
##################### COLLAR #########################
df.collar <- df.inf %>% select(Experiment, DPI, Individual, Gene.Control, Collar)
df.collar$Collar <- as.integer(df.collar$Collar)
df.collar <- filter(df.collar, Collar > 0)
Collar.Phenotypes <- count(df.collar, Gene.Control)
########## Add to Frequency Table ##########
FreqTable <- full_join(FreqTable, Collar.Phenotypes, by="Gene.Control")
colnames(FreqTable)[4] <- "Collar"
write.table(FreqTable, "PhenotypeFrequencyTable.txt", sep="\\t", quote=FALSE)
############# PRONOTUM #####################
pronotum.fields <- c("Pronotum", "D.T1.Posterior.margin.central", "D.T1.Posterior.margin.corners")

df.pronotum <- df.collar <- df.inf %>% select(Experiment, DPI, Individual, Gene.Control, Pronotum, 
                                              D.T1.Posterior.margin.central, D.T1.Posterior.margin.corners)

df.pronotum$Pronotum <- as.integer(df.pronotum$Pronotum)
df.pronotum$D.T1.Posterior.margin.central <- as.integer(df.pronotum$D.T1.Posterior.margin.central)
df.pronotum$D.T1.Posterior.margin.corners <- as.integer(df.pronotum$D.T1.Posterior.margin.corners)
dat.pronotum <- df.pronotum %>% select(Pronotum, D.T1.Posterior.margin.central, D.T1.Posterior.margin.corners)
dat.pronotum <- mutate(dat.pronotum, aggregate=apply(dat.pronotum, 1, sum, na.rm=TRUE))
df.pronotum <- data.frame(c(df.pronotum[,1:4], dat.pronotum))
df.pronotum <- filter(df.pronotum, aggregate>1)
Pronotum.Phenotypes <- count(df.pronotum, Gene.Control)
############ Add to Frequency Table #########
FreqTable <- full_join(FreqTable, Pronotum.Phenotypes, by="Gene.Control")
colnames(FreqTable)[5] <- "Pronotum"
write.table(FreqTable, "PhenotypeFrequencyTable.txt", sep="\\t", quote=FALSE)
############ ProJunction ###################
#propleuron.pronotum.junction.fields <- c(Pn.Pl.Junction, T1.LL.pronotum.pleural.junction, T1.RL.pronotum.pleural.junction)

df.pronpljunction <- df.inf %>% select(Experiment, DPI, Individual, Gene.Control, Pn.Pl.Junction, T1.LL.pronotum.pleural.junction, T1.RL.pronotum.pleural.junction)

df.pronpljunction$Pn.Pl.Junction <- as.integer(df.pronpljunction$Pn.Pl.Junction)
df.pronpljunction$T1.LL.pronotum.pleural.junction <- as.integer(df.pronpljunction$T1.LL.pronotum.pleural.junction)
df.pronpljunction$T1.RL.pronotum.pleural.junction <- as.integer(df.pronpljunction$T1.RL.pronotum.pleural.junction)

dat.pronpljunction <- select(df.pronpljunction, Pn.Pl.Junction, T1.LL.pronotum.pleural.junction, T1.RL.pronotum.pleural.junction)
dat.pronpljunction <- mutate(dat.pronpljunction, aggregate=apply(dat.pronpljunction, 1, sum, na.rm=TRUE))

df.pronpljunction <- data.frame(c(df.pronpljunction[,1:4], dat.pronpljunction))
df.pronpljunction <- filter(df.pronpljunction, aggregate > 0)
## only 281 of these!
PronPlJunction.Phenotypes <- count(df.pronpljunction, Gene.Control)
############ Add to Frequency Table #########
FreqTable <- full_join(FreqTable, PronPlJunction.Phenotypes, by="Gene.Control")
colnames(FreqTable)[6] <- "PronPlJunction"
write.table(FreqTable, "PhenotypeFrequencyTable.txt", sep="\\t", quote=FALSE)
############### FOREWINGS ###############
#forewings.fields <- c(L.FW.hemelytra, L.FW.membranous, R.FW.hemelytra, R.FW.membranous)
df.forewings <- df.inf %>% select(Experiment, DPI, Individual, Gene.Control, L.FW.hemelytra, L.FW.membranous, R.FW.hemelytra, R.FW.membranous)

df.forewings$L.FW.hemelytra <- as.integer(df.forewings$L.FW.hemelytra)
df.forewings$L.FW.membranous <- as.integer(df.forewings$L.FW.membranous)
df.forewings$R.FW.hemelytra <- as.integer(df.forewings$R.FW.hemelytra)
df.forewings$R.FW.membranous <- as.integer(df.forewings$R.FW.membranous)

dat.forewings <- df.forewings %>% select(L.FW.hemelytra, L.FW.hemelytra, R.FW.hemelytra, R.FW.membranous)
dat.forewings <- mutate(dat.forewings, aggregate = apply(dat.forewings, 1, sum, na.rm=TRUE))
df.forewings <- data.frame(c(df.forewings[,1:4], dat.forewings))
df.forewings <- filter(df.forewings, aggregate>0)

Forewings.Phenotype <- count(df.forewings, Gene.Control)
############ Add to Frequency Table #########
FreqTable <- full_join(FreqTable, Forewings.Phenotype, by="Gene.Control")
colnames(FreqTable)[7] <- "Forewings"
write.table(FreqTable, "PhenotypeFrequencyTable.txt", sep="\\t", quote=FALSE)
########## HIND WINGS ###################
#hindwings.fields <- c(L.HW, T3.HW.L.Wing.blade, R.HW, T3.HW.R.Wing.blade)
df.hindwings <- df.inf %>% select(Experiment, DPI, Individual, Gene.Control, L.HW, T3.HW.L.Wing.blade, R.HW, T3.HW.R.Wing.blade)

df.hindwings$L.HW <- as.integer(df.hindwings$L.HW)
df.hindwings$T3.HW.L.Wing.blade <- as.integer(df.hindwings$T3.HW.L.Wing.blade)
df.hindwings$R.HW <- as.integer(df.hindwings$R.HW)
df.hindwings$T3.HW.R.Wing.blade <- as.integer(df.hindwings$T3.HW.R.Wing.blade)

dat.hindwings <- df.hindwings %>% select(L.HW, T3.HW.L.Wing.blade, R.HW, T3.HW.R.Wing.blade)
dat.hindwings <- mutate(dat.hindwings, aggregate = apply(dat.hindwings, 1, sum, na.rm=TRUE))

df.hindwings <- data.frame(c(df.hindwings[,1:4], dat.hindwings))
df.hindwings <- filter(df.hindwings, aggregate > 0)
Hindwings.Phenotype <- count(df.hindwings, Gene.Control)
############ Add to Frequency Table #########
FreqTable <- full_join(FreqTable, Hindwings.Phenotype, by="Gene.Control")
colnames(FreqTable)[8] <- "Hindwings"
write.table(FreqTable, "PhenotypeFrequencyTable.txt", sep="\\t", quote=FALSE)
############ WINGS OVERALL #################
dat.wings <- data.frame(c(dat.hindwings, dat.forewings))
dat.wings <- dat.wings[,c(1:4,6:8)]
dat.wings <- mutate(dat.wings, aggregate = apply(dat.wings, 1, sum, na.rm=TRUE))

df.wings <- data.frame(c(df.inf[,1:5], dat.wings))
df.wings <- filter(df.wings, aggregate>0)
Wings.Phenotype <- count(df.wings, Gene.Control)
############ Add to Frequency Table #########
FreqTable <- full_join(FreqTable, Wings.Phenotype, by="Gene.Control")
colnames(FreqTable)[9] <- "Wings"
write.table(FreqTable, "PhenotypeFrequencyTable.txt", sep="\\t", quote=FALSE)

######### SCUTELLUM ##################
#scutellum.fields <- c(Scutellar.lobe, Mesonotum)

df.scutellum <- df.inf %>% select(Experiment, DPI, Individual, Gene.Control, Scutellar.lobe, Mesonotum)
df.scutellum$Scutellar.lobe <- as.integer(df.scutellum$Scutellar.lobe)
df.scutellum$Mesonotum <- as.integer(df.scutellum$Mesonotum)

dat.scutellum <- df.scutellum %>% select(Scutellar.lobe, Mesonotum)
dat.scutellum <- mutate(dat.scutellum, aggregate = apply(dat.scutellum, 1, sum, na.rm=TRUE))
df.scutellum <- data.frame(c(df.scutellum[,1:4], dat.scutellum))
df.scutellum <- filter(df.scutellum, aggregate > 0)
Scutellum.Phenotype <- count(df.scutellum, Gene.Control)
############ Add to Frequency Table #########
FreqTable <- full_join(FreqTable, Scutellum.Phenotype, by="Gene.Control")
colnames(FreqTable)[10] <- "Scutellum"
write.table(FreqTable, "PhenotypeFrequencyTable.txt", sep="\\t", quote=FALSE)
########## WING HINGES ################
# winghinges.fields <- c(Wing.hinges, Wing.sclerites, T2.FW.L.Hinge.Area, T2.FW.L.Axillary.sclerites, 
#                        T3.HW.L.Hinge.Area, T3.HW.L.Axillary.sclerites, T2.FW.R.Hinge.Area, T2.FW.R.Axillary.sclerites,
#                        T3.HW.R.Hinge.Area, T3.HW.R.Axillary.sclerites)

df.winghinges <- df.inf %>% select(Experiment, DPI, Individual, Gene.Control, Wing.hinges, Wing.sclerites, T2.FW.L.Hinge.Area, T2.FW.L.Axillary.sclerites, 
                                                           T3.HW.L.Hinge.Area, T3.HW.L.Axillary.sclerites, T2.FW.R.Hinge.Area, T2.FW.R.Axillary.sclerites,
                                                           T3.HW.R.Hinge.Area, T3.HW.R.Axillary.sclerites)

df.winghinges$Wing.hinges <- as.integer(df.winghinges$Wing.hinges)
df.winghinges$Wing.sclerites <- as.integer(df.winghinges$Wing.sclerites)
df.winghinges$T2.FW.L.Hinge.Area <- as.integer(df.winghinges$T2.FW.L.Hinge.Area)
df.winghinges$T2.FW.L.Axillary.sclerites <- as.integer(df.winghinges$T2.FW.L.Axillary.sclerites)
df.winghinges$T3.HW.L.Hinge.Area <- as.integer(df.winghinges$T3.HW.L.Hinge.Area)
df.winghinges$T3.HW.L.Axillary.sclerites <- as.integer(df.winghinges$T3.HW.L.Axillary.sclerites)
df.winghinges$T2.FW.R.Hinge.Area <- as.integer(df.winghinges$T2.FW.R.Hinge.Area)
df.winghinges$T2.FW.R.Axillary.sclerites <- as.integer(df.winghinges$T2.FW.R.Axillary.sclerites)
df.winghinges$T3.HW.R.Hinge.Area <- as.integer(df.winghinges$T3.HW.R.Hinge.Area)
df.winghinges$T3.HW.R.Axillary.sclerites <- as.integer(df.winghinges$T3.HW.R.Axillary.sclerites)

dat.winghinges <- df.winghinges[,5:14]
dat.winghinges <- mutate(dat.winghinges, aggregate = apply(dat.winghinges, 1, sum, na.rm=TRUE))

df.winghinges <- data.frame(c(df.winghinges[,1:4], dat.winghinges))

df.winghinges <- filter(df.winghinges, aggregate>0)
Winghinges.Phenotype <- count(df.winghinges, Gene.Control)
############ Add to Frequency Table #########
FreqTable <- full_join(FreqTable, Winghinges.Phenotype, by="Gene.Control")
colnames(FreqTable)[11] <- "WingHinges"
write.table(FreqTable, "PhenotypeFrequencyTable.txt", sep="\\t", quote=FALSE)
################ PLEURAL MARGINS ############
# pleuralmargins.fields <- c(Pos.Pl1, T1.LL.posterior.pleural.margin, T1.RL.posterior.pleural.margin, Pos.Pl2, 
#                            T2.LL.posterior.pleural.margin, T2.RL.posterior.pleural.margin, Pos.Pl.3, 
#                            T3.LL.post.pleural.margin, T3.RL.post.pleural.margin, T2.L.wing.base, T2.R.wing.base, 
#                            T3.LL.wing.base, T3.RL.wing.base)

df.pleuralmargins <- df.inf %>% select(Experiment, DPI, Individual, Gene.Control, Pos.Pl1, T1.LL.posterior.pleural.margin, T1.RL.posterior.pleural.margin, Pos.Pl2, 
                                       T2.LL.posterior.pleural.margin, T2.RL.posterior.pleural.margin, Pos.Pl.3, 
                                       T3.LL.post.pleural.margin, T3.RL.post.pleural.margin, T2.L.wing.base, T2.R.wing.base, 
                                       T3.LL.wing.base, T3.RL.wing.base)

dat.pleuralmargins <- df.pleuralmargins[,5:17]
dat.pleuralmargins$Pos.Pl1 <- as.integer(dat.pleuralmargins$Pos.Pl1)
dat.pleuralmargins$T1.LL.posterior.pleural.margin <- as.integer(dat.pleuralmargins$T1.LL.posterior.pleural.margin)
dat.pleuralmargins$T1.RL.posterior.pleural.margin <- as.integer(dat.pleuralmargins$T1.RL.posterior.pleural.margin)
dat.pleuralmargins$Pos.Pl2 <- as.integer(dat.pleuralmargins$Pos.Pl2)
dat.pleuralmargins$T2.LL.posterior.pleural.margin <- as.integer(dat.pleuralmargins$T2.LL.posterior.pleural.margin)
dat.pleuralmargins$T2.RL.posterior.pleural.margin <- as.integer(dat.pleuralmargins$T2.RL.posterior.pleural.margin)
dat.pleuralmargins$Pos.Pl.3 <- as.integer(dat.pleuralmargins$Pos.Pl.3)
dat.pleuralmargins$T3.LL.post.pleural.margin <- as.integer(dat.pleuralmargins$T3.LL.post.pleural.margin)
dat.pleuralmargins$T3.RL.post.pleural.margin <- as.integer(dat.pleuralmargins$T3.RL.post.pleural.margin)
dat.pleuralmargins$T2.L.wing.base <- as.integer(dat.pleuralmargins$T2.L.wing.base)
dat.pleuralmargins$T2.R.wing.base <- as.integer(dat.pleuralmargins$T2.R.wing.base)
dat.pleuralmargins$T3.LL.wing.base <- as.integer(dat.pleuralmargins$T3.LL.wing.base)
dat.pleuralmargins$T3.RL.wing.base <- as.integer(dat.pleuralmargins$T3.RL.wing.base)

dat.pleuralmargins <- mutate(dat.pleuralmargins, aggregate = apply(dat.pleuralmargins, 1, sum, na.rm=TRUE))

df.pleuralmargins <- data.frame(c(df.pleuralmargins[,1:4], dat.pleuralmargins))
df.pleuralmargins <- filter(df.pleuralmargins, aggregate>0)

PleuralMargins.Phenotype <- count(df.pleuralmargins, Gene.Control)
############ Add to Frequency Table #########
FreqTable <- full_join(FreqTable, PleuralMargins.Phenotype, by="Gene.Control")
colnames(FreqTable)[12] <- "PleuralMargins"
write.table(FreqTable, "PhenotypeFrequencyTable.txt", sep="\\t", quote=FALSE)
############################################
# Get frequencies from FreqTable           #
############################################

FreqTable <- FreqTable %>% mutate(CollarFreq = (Collar/Total)) %>%
  mutate(ScLobesFreq = ScLobes/Total) %>%
  mutate(PronotumFreq = Pronotum/Total) %>%
  mutate(PronPlJunctionFreq = PronPlJunction/Total) %>%
  mutate(ForewingsFreq = Forewings/Total) %>%
  mutate(HindwingsFreq = Hindwings/Total) %>%
  mutate(ScutellumFreq = Scutellum/Total) %>%
  mutate(WingHingesFreq = WingHinges/Total) %>%
  mutate(PleuralMarginsFreq = PleuralMargins/Total) %>%
  mutate(WingsFreq = Wings/Total)

colnames(FreqTable)

FreqTable_2 <- data.frame(FreqTable$Gene.Control, FreqTable$Total, FreqTable$ScLobes, 
                          round(FreqTable$ScLobesFreq, 2), FreqTable$Collar, 
                          round(FreqTable$CollarFreq, 2), 
                          FreqTable$Pronotum, round(FreqTable$PronotumFreq, 2),  
                          FreqTable$PronPlJunction, round(FreqTable$PronPlJunctionFreq, 2),
                          FreqTable$Forewings, round(FreqTable$ForewingsFreq, 2),
                          FreqTable$Hindwings, round(FreqTable$HindwingsFreq, 2),
                          FreqTable$Wings, round(FreqTable$WingsFreq, 2),
                          FreqTable$Scutellum, round(FreqTable$ScutellumFreq, 2),
                          FreqTable$WingHinges, round(FreqTable$WingHingesFreq, 2),
                          FreqTable$PleuralMargins, round(FreqTable$PleuralMarginsFreq, 2))

###########################################
# Write some subset tables                #
###########################################

df.InfoCols <- select(df.inf, Experiment, DPI, Individual, Gene.Control, Treatment)
df.collar <- df.collar
df.forewings <- data.frame(c(df.InfoCols, dat.forewings))
df.hindwings <- data.frame(c(df.InfoCols, dat.hindwings))
df.wings <- data.frame(c(df.InfoCols, dat.wings))
df.pleuralmargins <- data.frame(c(df.InfoCols, dat.pleuralmargins))
df.pronotum <- data.frame(c(df.InfoCols, dat.pronotum))
df.pronpljunction <- data.frame(c(df.InfoCols, dat.pronpljunction))
df.sclobes <- data.frame(c(df.InfoCols, dat.sclobes))
df.winghinges <- data.frame(c(df.InfoCols, dat.winghinges))
df.scutellum <- data.frame(c(df.InfoCols, dat.scutellum))

setwd("DataSubsets/")
write.table(df.collar, "CollarPhenotypes.txt", sep="\\t", quote=FALSE)
write.table(df.forewings, "ForewingPhenotypes.txt", sep="\\t", quote=FALSE)
write.table(df.hindwings, "HindwingPhenotypes.txt", sep="\\t", quote=FALSE)
write.table(df.pleuralmargins, "PleuralMarginPhenotypes.txt", sep="\\t", quote=FALSE)
write.table(df.pronotum, "PronotumPhenotypes.txt", sep="\\t", quote=FALSE)
write.table(df.pronpljunction, "PronPlJunctionPhenotypes.txt", sep="\\t", quote=FALSE)
write.table(df.sclobes, "SupracoxalLobesPhenotypes.txt", sep="\\t", quote=FALSE)
write.table(df.scutellum, "ScutellumPhenotypes.txt", sep="\\t", quote=FALSE)
write.table(df.winghinges, "WinghingesPhenotypes.txt", sep="\\t", quote=FALSE)
write.table(df.wings, "WingsPhenotypes.txt", sep="\\t", quote=FALSE)

colnames(FreqTable_2) <- c(
  "Gene.or.Control", 
  "Total", 
  "ScLobes", 
  "ScLobes.Freq", 
  "Collar", 
  "Collar.Freq", 
  "Pronotum", 
  "Pronotum.Freq", 
  "PronPlJunction", 
  "PronPlJunction.Freq", 
  "Forewings", 
  "Forewings.Freq", 
  "Hindwings", 
  "Hindwings.Freq", 
  "Wings.either", 
  "Wings.either.Freq",
  "Scutellum", 
  "Scutellum.Freq", 
  "Wing.Hinges", 
  "Wing.Hinges.Freq", 
  "Pleural.Margins", 
  "Pleural.Margins.Freq"
)
write.table(FreqTable_2, "PhenotypeFrequencies.txt", sep="\\t", quote=FALSE)




Penetrance <- select(FreqTable, Gene.Control, Total, ScLobesFreq, CollarFreq, 
                     PronotumFreq, PronPlJunctionFreq, ForewingsFreq, 
                     HindwingsFreq, WingsFreq, ScutellumFreq, WingHingesFreq, 
                     PleuralMarginsFreq)
Penetrance <- data.frame(Penetrance$Gene.Control, Penetrance$Total, 
                         round(Penetrance$ScLobesFreq, 2), 
                         round(Penetrance$CollarFreq, 2), 
                         round(Penetrance$PronotumFreq, 2), 
                         round(Penetrance$PronPlJunctionFreq, 2),
                         round(Penetrance$ForewingsFreq, 2),
                         round(Penetrance$HindwingsFreq, 2), 
                         round(Penetrance$WingsFreq, 2), 
                         round(Penetrance$ScutellumFreq, 2), 
                         round(Penetrance$PleuralMarginsFreq, 2),
                         round(Penetrance$WingHingesFreq, 2))

colnames(Penetrance) <- c(
  "Gene.or.Control", 
  "Total", 
  "Penetrance.ScLobes", 
  "Penetrance.Collar", 
  "Penetrance.Pronotum", 
  "Penetrance.PronPlJunction", 
  "Penetrance.Forewings", 
  "Penetrance.Hindwings", 
  "Penetrance.Wings", 
  "Penetrance.Scutellum", 
  "Penetrance.PleuralMargins", 
  "Penetrance.WingHinges"
)
write.table(Penetrance, file="PhenotypePenetrance.txt", quote=FALSE, sep="\\t")

###############################################################################
# Save RData workspace 
save.image(file="../RNAi_Data.RData")
###############################################################################
# Save session info

writeLines(capture.output(sessionInfo()), "../sessionInfo.txt")
