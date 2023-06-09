#BREEDING SYSTEM ASSESSMENT
library(MASS)
#Seed Count Data for Breeding Systems Results
setwd("~/Documents/Preparation for Manuscript")
gc<-read.csv("GC_Crosses_December_2021.csv")
labels(gc)
#all gayleana data, to determine cross effect on total seeds
Gall <- gc[which(gc$Species == 'G'), ]
#Put cross types in order of hypothesis
Gall$Cross <- factor(Gall$Cross, levels=c("BS","AUT","GEIT","SIB","WIN","BW"))
levels(Gall$Cross)
aov.Gall = aov(TotSeeds ~ Cross, data = Gall)
summary(aov.Gall)
#significant differences when all cross types considered

#outcrossing gayleana data
Goutcross <- gc[which((gc$Species == 'G' & gc$Cross == 'SIB')
                      | (gc$Species == 'G' & gc$Cross == 'WIN')
                      | (gc$Species == 'G' & gc$Cross == 'BW')), ]
summary(Goutcross)
hist(Goutcross$TotSeeds)
aov.Goutcross = aov(TotSeeds ~ Cross, data = Goutcross)
summary(aov.Goutcross) #Cross not sig
TukeyHSD(aov.Goutcross) #No pairwise comparison is sig (BW/SIB, BW/WIN, SIB/WIN)
#compare by maternal population
boxplot(TotSeeds~MaternalPop, data = Goutcross, main = "G Outcross")
aov.GoutcrossPop = aov(TotSeeds ~ MaternalPop, data = Goutcross)
summary(aov.GoutcrossPop) #MatPop is sig for outcrossed data
TukeyHSD(aov.GoutcrossPop)
#compare by maternal line
aov.GoutcrossMatLine = aov(TotSeeds ~ Maternal, data = Goutcross)
summary(aov.GoutcrossMatLine)

#Self crosses gayleana
Gself <- gc[which((gc$Species == 'G' & gc$Cross == 'AUT' ) 
                  | (gc$Species == 'G' & gc$Cross == 'GEIT')), ]
summary(Gself)
hist(Gself$TotSeeds)
aov.GselfType <- aov(TotSeeds ~ Cross, data = Gself)
summary(aov.GselfType) #sig difference between hand-pollination autogamous and unmanipulated controls 
#Note, only 14 self-crosses produced seed.

#compare self and outcrossed treatments O.gayleana
aov.GCrossType <- (aov(TotSeeds ~ Cross_Type, data = Gall))
summary(aov.GCrossType)
summary(Gall) #very sig difference between self and outcrossed groups of crosses

#O hart filifolia
Fall <- gc[which(gc$Species == 'F' & !gc$Cross == 'GF'), ]
#outcrossed filifolia 
Foutcross <- gc[which((gc$Species == 'F' & gc$Cross == 'SIB') 
                      | (gc$Species == 'F' & gc$Cross == 'WIN') 
                      | (gc$Species == 'F' & gc$Cross == 'BW')), ]
summary(Foutcross)
aov.Foutcross = aov(TotSeeds ~ Cross, data = Foutcross)
summary(aov.Foutcross) #Cross type not sig for outcrossed fili
TukeyHSD(aov.Foutcross) #No pairwise comparison sig (BW/SIB, BW/WIN, SIB/WIN)
#compare by maternal population
aov.FoutcrossPop = aov(TotSeeds ~ MaternalPop, data = Foutcross)
boxplot(TotSeeds ~ MaternalPop, data = Foutcross)
summary(aov.FoutcrossPop) #MatPop is sig for outcrossed fili

#self crosses filifolia
Fself<-gc[which((gc$Species == 'F' & gc$Cross == 'AUT' ) 
                | (gc$Species == 'F' & gc$Cross == 'GEIT')), ]
summary(Fself)
aov.Fselftype <- aov(TotSeeds ~ Cross, data = Fself)
summary(aov.Fselftype) #sig difference between autogamous hand pollination crosses and unmanipulated controls
Fself
Fself2<-gc

#to determine difference between self and outcrosses O hart filifolia
aov.FCrossType <- (aov(TotSeeds ~ Cross_Type, data = Fall))
summary(aov.FCrossType)
summary(Fall) #very sig difference between self and outcrossed groups of crosses


#FLORAL MORPHOLOGY
#Figure 2 - PCA of both species, floral measurements including herkogamy, 
  #tube length, corolla diameter (avg), and floral flare
#PCA with only first flower, herkogamy is left blank if no style measurement
setwd("~/Documents/Preparation for Manuscript")
fm.2 <- read.csv("FloralMorphGC_noYA_2022.csv")
labels(fm.2)
fm.pca3 <- fm.2[which(fm.2$Flower == '1'), c("Species", "Site", "corolla_avg", "Floral_flare", "Tube_length", "Herkogamy")] 
fm.pca3<-na.omit(fm.pca3)
summary(fm.pca3)
pca3<-prcomp(fm.pca3[,3:6], center=T, scale=T) #This also standardizes data, note: going to column 6 adds Herkogamy
Species <- factor(fm.pca3$Species) #I don't actually think this is a necessary step.
Site <- factor(fm.pca3$Site)
library(ggplot2)
library(ggfortify)
library(vegan)
library(cluster)
summary(pca3)
pca3
head(pca3$rotation) #find loadings
autoplot((pca3), data = fm.pca3, colour = 'Site', shape = 'Species', loadings = TRUE, loadings.colour = 'black', loadings.label = TRUE,
         main = "Both species with herkogamy") #PCA Scatterplot. This one works w/ loadings
#shaded ellipses
autoplot(pca3, data = fm.pca3, colour = 'Site', shape = 'Species', loadings = TRUE, loadings.colour = 'black', loadings.label = TRUE,
         main = "Both species with herkogamy", frame = TRUE, frame.type = 'norm', frame.colour = 'Species') #somewhat unexpectedly, this works but has filled in circles.
#Transparent ellipses
autoplot(pca3, data = fm.pca3, colour = 'Site', shape = 'Species', loadings = TRUE, loadings.colour = 'black', loadings.label = TRUE,
         main = "Both species with herkogamy and first flower", frame = TRUE, frame.type = 'norm', frame.colour = 'Species', frame.alpha = 0) #somewhat unexpectedly, this works. frame.alpha changes the transparency.

#Simple Floral Stats
setwd("~/Documents/Preparation for Manuscript")
fm<-read.csv("FloralMorphGC_noYA_2021.csv")
#Are floral morph (fm) qualities different between O. gayleana and O. hart fili?
aov.corolla = aov(corolla_avg ~ Species, data = fm)
summary(aov.corolla)
aov.tube = aov(Tube_length ~ Species, data = fm)
summary(aov.tube)
aov.flare = aov(Floral_flare ~ Species, data = fm)
summary(aov.flare)
aov.herkogamy = aov(Herkogamy ~ Species, data = fm)
summary(aov.herkogamy)
boxplot(Herkogamy~Species, data = fm)
boxplot(Tube_length~Species, data = fm)
boxplot(corolla_avg~Species, data = fm)
boxplot(Floral_flare~Species, data = fm)
#Yes, highly significant floral morph characteristics for corolla, tube, and flare. Herk is not sig.

#Are floral morph qualities different between populations for O. gayleana? Growth chamber plants only.
Gmorph <- fm[which(fm$Species == 'G'), ]
aov.Gsite = aov(corolla_avg ~ Site, data = Gmorph)
summary(aov.Gsite)
boxplot(corolla_avg~Site, data = Gmorph, xlab = "O. gayleana Maternal Pop Site", ylab = "Corolla Diameter")
TukeyHSD(aov.Gsite)
#corolla = yes - Croton is sig different from all other sites
aov.Gsite2 = aov(Tube_length ~ Site, data = Gmorph)
summary(aov.Gsite2)
boxplot(Tube_length~Site, data = Gmorph, xlab = "O. gayleana Maternal Pop Site", ylab = "Floral Tube Length")
TukeyHSD(aov.Gsite2)
#tube length = yes, various pops are sig different from each other
aov.Gsite3 = aov(Floral_flare ~ Site, data = Gmorph)
summary(aov.Gsite3)
boxplot(Floral_flare~Site, data = Gmorph, xlab = "O. gayleana Maternal Pop Site", ylab = "Floral flare")
TukeyHSD(aov.Gsite3)
#floralflare = yes - Croton is sig different from all other sites
aov.Gsite4 = aov(Herkogamy~Site, data = Gmorph)
summary(aov.Gsite4)
boxplot(Herkogamy~Site, data = Gmorph, xlab = "O.gayleana Maternal Pop Site")
TukeyHSD(aov.Gsite4)
#Herkogamy not sig for O. gay

#O.hart fili stats
#Are floral morph qualities different between populations for O. hart. fili?
Fmorph <- fm[which(fm$Species == 'F'), ]
aov.Fsite = aov(corolla_avg ~ Site, data = Fmorph)
summary(aov.Fsite)
boxplot(corolla_avg~Site, data = Fmorph, xlab = "O.hart.fili Maternal Pop Site", ylab = "Corolla Diameter")
TukeyHSD(aov.Fsite)
#corolla = yes - Croton is sig different from three other sites
aov.Fsite2 = aov(Tube_length ~ Site, data = Fmorph)
summary(aov.Fsite2)
boxplot(Tube_length~Site, data = Fmorph, xlab = "O.hart.fili Maternal Pop Site", ylab = "Floral Tube Length")
TukeyHSD(aov.Fsite2)
#tube length = yes, YB is sig diff from PhD and YA
aov.Fsite3 = aov(Floral_flare ~ Site, data = Fmorph)
summary(aov.Fsite3)
boxplot(Floral_flare~Site, data = Fmorph, xlab = "O.hart.fili Maternal Pop Site", ylab = "Floral Flare")
TukeyHSD(aov.Fsite3)
#floralflare = No
aov.Fsite4 = aov(Herkogamy ~ Site, data = Fmorph)
summary(aov.Fsite4)
boxplot(Herkogamy~Site, data = Fmorph, xlab = "O.hart.fili Maternal Pop Site", ylab = "Herkogamy")
TukeyHSD(aov.Fsite4)
#Herkogamy = yes significant, SH is diff from PhD

#Are fm qualities different between wild plants and growth chamber plants? 
# May be put in supp materials, not for main manuscript
#Only including sites with data from both GC and wild (YA and YB for both, SH for O. gay)
gcwild<-read.csv("GCWild_FloralMorphMaster.csv")
Ggcwild <- gcwild[which(gcwild$Species == 'G'), ]
aov.Ggcwild = aov(corolla_avg ~ Location, data = Ggcwild)
summary(aov.Ggcwild)
boxplot(corolla_avg~Location, data = Ggcwild)
aov.Ggcwild2 = aov(Tube_length ~ Location, data = Ggcwild)
summary(aov.Ggcwild2)
boxplot(Tube_length~Location, data = Ggcwild)
aov.Ggcwild3 = aov(Floral_flare ~ Location, data = Ggcwild)
summary(aov.Ggcwild3)
#O gay - Yes, Corolla length and tube length are significantly larger for growth chamber O. gayleana plants
#Floral flare not significant
Fgcwild <- gcwild[which(gcwild$Species == 'F'), ]
aov.Fgcwild = aov(corolla_avg ~ Location, data = Fgcwild)
summary(aov.Fgcwild)
boxplot(corolla_avg~Location, data = Fgcwild)
aov.Fgcwild2 = aov(Tube_length ~ Location, data = Fgcwild)
summary(aov.Fgcwild2)
boxplot(Tube_length~Location, data = Fgcwild)
aov.Fgcwild3 = aov(Floral_flare ~ Location, data = Fgcwild)
summary(aov.Fgcwild3)
boxplot(Floral_flare~Location, data = Fgcwild)
#O hart fili - Yes, corolla length, tube length, and floral flare all sig larger for growth chamber plants

#POPULATION GENETIC PARAMETERS
#ANOVAs for genetic data. Need to differentiate landscape scale so that we are not oversampling Yeso metapops
setwd("~/Documents/Preparation for Manuscript")
caly<-read.csv("GeneticAnalysis.csv", header=T)
summary(caly)
names(caly)
library(ggplot2)
landscape <- caly[which(caly$Scale == 'landscape'), ]
summary(landscape) #note that these values include both species
local <- caly[which(caly$Scale == 'local'), ] 
aov.Na = aov(Na ~ Sp, data = landscape) #Na = number of alleles (A in paper)
summary(aov.Na)
boxplot(Na~Sp, data = landscape)
aov.Nalocal = aov(Na ~ Sp, data = local)
summary(aov.Nalocal)
boxplot(Na~Sp, data = local)

aov.Ne = aov(Ne ~ Sp, data = landscape) #Ne = number of effective alleles
summary(aov.Ne)
boxplot(Ne~Sp, data = landscape)

aov.He = aov(He ~ Sp, data = landscape) #He = expected heterozygosity
summary(aov.He)
boxplot(He~Sp, data = landscape)

aov.Ho = aov(Ho ~ Sp, data = landscape) #Ho = observed heterozygosity
summary(aov.Ho)
boxplot(Ho~Sp, data = landscape)

aov.Fis = aov(Fis ~ Sp, data = landscape) #Fis = inbreeding coefficient. 
summary(aov.Fis)
boxplot(Fis~Sp, data = landscape)

aov.Prv = aov(Prv ~ Sp, data = landscape)
summary(aov.Prv)
boxplot(Prv~Sp, data = landscape)

aov.R = aov(Relatedness ~ Sp, data = landscape)
summary(aov.R)
boxplot(Relatedness~Sp, data = landscape)

R<- as.numeric(caly$Relatedness)
FisNum<-as.numeric(caly$Fis)#maybe don't need to do that?
RF<-lm(formula = Relatedness ~ Fis, data = caly) #this is all scales and both species
summary(RF)
OhfRlocal <- local[which(local$Sp == 'Ohf'), ] 
OhfRlocal
RFlocalF <- lm(formula = Fis ~ Relatedness, data = OhfRlocal)
summary(RFlocalF)
OgRlocal <- local[which(local$Sp == 'G'), ] 
OgRlocal
RFlocalG <- lm(formula = Fis ~ Relatedness, data = OgRlocal)
summary(RFlocalG)
OhfRlandscape <- landscape[which(landscape$Sp == 'Ohf'), ]
OhfRlandscape
RFlandscapeF <- lm(formula = Fis ~ Relatedness, data = OhfRlandscape)
summary(RFlandscapeF)
OgRlandscape <- landscape[which(landscape$Sp == 'G'), ]
OgRlandscape
RFlandscapeG <- lm(formula = Fis ~ Relatedness, data = OgRlandscape)
summary(RFlandscapeG)
#inbreeding did not vary significantly by relatedness for either species at either scale.

#Figure 3
#FST vs Geographic Distance Graph - all populations
#With all pops, ANCOVA is significant and Mantel is significant for both species
install.packages("ggplot2")
library("ggplot2")
setwd("~/Documents/Preparation for Manuscript")
fst<-read.csv("FST_GeogDistance.csv")
labels(fst)
Species <- factor(fst$Species)
Distance <- factor(fst$Distance)
FST_Linearized <- as.numeric(fst$FST_Lin_Spag) #Note - need to make this numeric so the scale appears correctly on y-axis, can't be a factor
FST<- as.numeric(fst$FST)
plot.default(fst$Distance, fst$FST_Lin_Spag, data = fst) #scatterplot of both G and F with no differentiation
Fili <- fst[which(fst$Species == 'Ohf'), ]
summary(Fili)
Gay <- fst[which(fst$Species == 'G'), ]
summary(Gay)
ggplot(fst, aes(Distance,FST_Linearized)) + geom_point(aes(color = Species))
ggplot(fst, aes(Distance,FST)) + geom_point(aes(color = Species)) #just to see if FST values non-linearized looks the same
fitF<- lm(FST_Lin_Spag~Distance, data = Fili) #fit regression line for O. hart fili
fitF  #Fili equation: Fst = 0.0001279x + 0.02555
summary(fitF)  #Adjusted Rsquared value: 0.443, F-stat: 22.47 on 1 and 26 DF, pvalue = 6.67e-05
fitG <- lm(FST_Lin_Spag~Distance, data = Gay) #fit regression line for O. gay
fitG #O.gay equation: Fst = 0.0004016x + 0.1368
summary(fitG)  #Adjusted Rsquared value: 0.3633, F-stat: 16.4 on 1 and 26 DF, pvalue = 0.00041
#color graph
ggplot(fst, aes(Distance,FST_Linearized)) + geom_point(aes(color = Species)) + geom_abline(intercept = 0.02555, slope = 0.0001279) + geom_abline(intercept = 0.1368, slope = 0.0004016)
#black and white graph
ggplot(fst, aes(Distance,FST_Linearized)) + geom_point(aes(shape = Species)) + scale_shape_manual(values = c(2,19)) + geom_abline(intercept = 0.02555, slope = 0.0001279) + geom_abline(intercept = 0.1368, slope = 0.0004016)
#R shapes - 1 is open circle, 2 is open triangle, 17 is filled triangle

#ANCOVA to compare the regression lines to each other to determine if they are different in slope and intercept.
#categorical factor = species, x-var= distance, y-var = fst
model<-aov(FST_Lin_Spag~Distance*Species, data = fst)
summary(model)
#because the interaction of species and distance is significant, we can say that the slope of the regression lines are significantly different

#Mantel Test to check for Isolation by Distance - are genetic distances and geographical distances between populations correlated?
install.packages("ade4")
install.packages("vegan")
library(ade4)
library(vegan)
dist<-read.csv("Matrix_Distance.csv")
is.matrix(dist) #not a workable matrix
distmatrix <- as.dist(dist)
distmatrix
class(distmatrix) #this should be the class necessary for Mantel test

#O.gayleana with linearized FST and distance values from SPAGeDi
fstgaylin<-read.csv("Matrix_FSTGay_Lin.csv")
fstgaylinmatrix<-as.dist(fstgaylin)
fstgaylinmatrix
mantel.rtest(distmatrix, fstgaylinmatrix, nrepet = 10000) 
mantel(distmatrix, fstgaylinmatrix, method = "spearman", permutations = 10000, na.rm = TRUE)
#Both significant
#O. gay Mantel statistic r = 0.6634, p = 0.0052. This is a strong relationship between Fst and distance (0.6634)
#O.hart fili
fstfili<-read.csv("Matrix_FSTFili_Lin.csv")
fstfilimatrix<-as.dist(fstfili)
fstfilimatrix
mantel.rtest(distmatrix, fstfilimatrix, nrepet = 10000) #hmm. Not significant correlation between distance and Fst
mantel(distmatrix, fstfilimatrix, method = "spearman", permutations = 10000, na.rm = TRUE)
#Both tests significant
#Ohf Mantel statistic r = 0.7094, p = 0.0063. This is also a strong relationship.

#Chloroplast Data - Landscape
setwd("~/Documents/Preparation for Manuscript")
chloro<-read.csv("GeneticAnalysis_Chloro.csv", header=T)
summary(chloro)
names(chloro)
library(ggplot2)
landscape <- chloro[which(chloro$Scale == 'landscape'), ]
aov.NaChloro = aov(Na ~ Sp, data = landscape) #Na = number of alleles (A in paper)
summary(aov.NaChloro)
boxplot(Na~Sp, data = landscape)
aov.NeChloro = aov(Ne ~ Sp, data = landscape)
summary(aov.NeChloro)
boxplot(Ne~Sp, data = landscape)
aov.IChloro = aov(I ~ Sp, data = landscape)
summary(aov.IChloro)
aov.PrvChloro = aov(Unique ~ Sp, data= landscape)
summary(aov.PrvChloro)
boxplot(Unique~Sp, data = landscape)
aov.HaploTotal = aov(Haplo_Total ~ Sp, data = landscape)
summary(aov.HaploTotal)

#Yeso Chloro
local<- chloro[which(chloro$Scale == 'local'), ]
aov.NaChloroYeso = aov(Na ~ Sp, data = local)
summary(aov.NaChloroYeso)
aov.NeChloroYeso = aov(Ne ~ Sp, data = local)
summary(aov.NeChloroYeso)
aov.IChloroYeso = aov(I ~ Sp, data = local)
summary(aov.IChloroYeso)
aov.PrvChloroYeso = aov(Unique ~ Sp, data = local)
summary(aov.PrvChloroYeso)
boxplot(Unique~Sp, data = local)
aov.HaploTotalYeso = aov(Haplo_Total ~ Sp, data = local)
summary(aov.HaploTotalYeso)