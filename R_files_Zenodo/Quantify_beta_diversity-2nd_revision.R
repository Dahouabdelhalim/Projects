
install.packages("adespatial")
install.packages("adegraphics")
install.packages("goeveg")
install.packages("betapart")
install.packages("vegan")
install.packages("MASS")
install.packages("labdsv")
install.packages("car")
install.packages("vegan3d")
install.packages("mvoutlier")
install.packages("BiodiversityR")
install.packages("deldir")
install.packages("ade4")


library(adespatial)
library(adegraphics)
library(goeveg)
library(betapart)
library(vegan)
library(MASS)
library(labdsv)
library(car)
library(vegan3d)
library(mvoutlier)
library(BiodiversityR)
library(deldir)
library(ade4)

? beta.div.comp


###################################### Quantitative form ¦Â-Diversity

####### mangrove community
##beta-entire study area

mangrove=read.csv("mangrove_beta.csv", head=T, row.names = 1)
head(mangrove)
mangrove.bd.part <- beta.div.comp(mangrove, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
mangrove.bd.part
summary(mangrove.bd.part$D)# summary of Dissimilarity matrix

mangrove.bd.part1 <- beta.div.comp(mangrove, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
mangrove.bd.part1
summary(mangrove.bd.part1$D)# summary of Dissimilarity matrix


##### beta-reference site
mangrove_m=read.csv("mangrove_beta_m.csv", head=T, row.names = 1)
head(mangrove_m)
mangrove.bd.part <- beta.div.comp(mangrove_m, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
mangrove.bd.part
summary(mangrove.bd.part$D)# summary of Dissimilarity matrix

mangrove.bd.part1 <- beta.div.comp(mangrove_m, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
#coef="S","BJ": Baselga family, Jaccard-based indices
mangrove.bd.part1
summary(mangrove.bd.part1$D)# summary of Dissimilaritry matrix


######beta-passive restoration

mangrove_n=read.csv("mangrove_beta_P.csv", head=T, row.names = 1)
head(mangrove_n)
mangrove.bd.part <- beta.div.comp(mangrove_n, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
mangrove.bd.part
summary(mangrove.bd.part$D)# summary of Dissimilarity matrix

mangrove.bd.part1 <- beta.div.comp(mangrove_n, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
#coef="S","BJ": Baselga family, Jaccard-based indices
mangrove.bd.part1
summary(mangrove.bd.part1$D)# summary of Dissimilarity matrix


######beta-active restoration

mangrove_sn=read.csv("mangrove_beta_A.csv", head=T, row.names = 1)
head(mangrove_sn)
mangrove.bd.part <- beta.div.comp(mangrove_sn, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
mangrove.bd.part
summary(mangrove.bd.part$D)# summary of Dissimilarity matrix

mangrove.bd.part1 <- beta.div.comp(mangrove_sn, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
#coef="S","BJ": Baselga family, Jaccard-based indices
mangrove.bd.part1
summary(mangrove.bd.part1$D)# summary of Dissimilarity matrix



############################################### waterbirds
##beta-entire study area

bird=read.csv("bird_beta.csv", head=T, row.names = 1)
head(bird)
bird.bd.part <- beta.div.comp(bird, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
bird.bd.part
summary(bird.bd.part$D)# summary of Dissimilarity matrix

bird.bd.part1 <- beta.div.comp(bird, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
bird.bd.part1
summary(bird.bd.part1$D)# summary of Dissimilarity matrix


##### beta-reference site
bird_m=read.csv("bird_beta_m.csv", head=T, row.names = 1)
head(bird_m)
bird.bd.part <- beta.div.comp(bird_m, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
bird.bd.part
summary(bird.bd.part$D)# summary of Dissimilarity matrix

bird.bd.part1 <- beta.div.comp(bird_m, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
#coef="S","BJ": Baselga family, Jaccard-based indices
bird.bd.part1
summary(bird.bd.part1$D)# summary of Dissimilarity matrix


######beta-passive restoration site

bird_n=read.csv("bird_beta_P.csv", head=T, row.names = 1)
head(bird_n)
bird.bd.part <- beta.div.comp(bird_n, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
bird.bd.part
summary(bird.bd.part$D)# summary of Dissimilarity matrix

bird.bd.part1 <- beta.div.comp(bird_n, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
#coef="S","BJ": Baselga family, Jaccard-based indices
bird.bd.part1
summary(bird.bd.part1$D)# summary of Dissimilarity matrix


######beta-active restoration site

bird_sn=read.csv("bird_beta_A.csv", head=T, row.names = 1)
head(bird_sn)
bird.bd.part <- beta.div.comp(bird_sn, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
bird.bd.part
summary(bird.bd.part$D)# summary of Dissimilarity matrix

bird.bd.part1 <- beta.div.comp(bird_sn, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
#coef="S","BJ": Baselga family, Jaccard-based indices
bird.bd.part1
summary(bird.bd.part1$D)# summary of Dissimilarity matrix


############################################################### fish

#################beta-entire study area

fish=read.csv("fish_beta.csv", head=T, row.names = 1)
head(fish)
fish.bd.part <- beta.div.comp(fish, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
fish.bd.part
summary(fish.bd.part$D)# summary of Dissimilarity matrix

fish.bd.part1 <- beta.div.comp(fish, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
fish.bd.part1
summary(fish.bd.part1$D)# summary of Dissimilarity matrix


##### beta-reference site
fish_m=read.csv("fish_beta_m.csv", head=T, row.names = 1)
head(fish_m)
fish.bd.part <- beta.div.comp(fish_m, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
fish.bd.part
summary(fish.bd.part$D)# summary of Dissimilarity matrix

fish.bd.part1 <- beta.div.comp(fish_m, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
#coef="S","BJ": Baselga family, Jaccard-based indices
fish.bd.part1
summary(fish.bd.part1$D)# summary of Dissimilarity matrix


######beta-passive restoration site

fish_n=read.csv("fish_beta_P.csv", head=T, row.names = 1)
fish.bd.part <- beta.div.comp(fish_n, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
fish.bd.part
summary(fish.bd.part$D)# summary of Dissimilarity matrix

fish.bd.part1 <- beta.div.comp(fish_n, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
#coef="S","BJ": Baselga family, Jaccard-based indices
fish.bd.part1
summary(fish.bd.part1$D)# summary of Dissimilarity matrix


######beta-active restoration site

fish_sn=read.csv("fish_beta_A.csv", head=T, row.names = 1)
head(fish_sn)
fish.bd.part <- beta.div.comp(fish_sn, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
fish.bd.part
summary(fish.bd.part$D)# summary of Dissimilarity matrix

fish.bd.part1 <- beta.div.comp(fish_sn, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
#coef="S","BJ": Baselga family, Jaccard-based indices
fish.bd.part1
summary(fish.bd.part1$D)# summary of Dissimilarity matrix



############################################################################# macrobenthos

############# beta-entire study area

benthos=read.csv("benthos_beta.csv", head=T, row.names = 1)

benthos.bd.part <- beta.div.comp(benthos, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
benthos.bd.part
summary(benthos.bd.part$D)# summary of Dissimilarity matrix

benthos.bd.part1 <- beta.div.comp(benthos, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
benthos.bd.part1
summary(benthos.bd.part1$D)# summary of Dissimilarity matrix


##### beta-reference site
benthos_m=read.csv("benthos_beta_m.csv", head=T, row.names = 1)
benthos.bd.part <- beta.div.comp(benthos_m, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
benthos.bd.part
summary(benthos.bd.part$D)# summary of Dissimilarity matrix

benthos.bd.part1 <- beta.div.comp(benthos_m, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
#coef="S","BJ": Baselga family, Jaccard-based indices
benthos.bd.part1
summary(benthos.bd.part1$D)# summary of Dissimilarity matrix


######beta-passive restoration site

benthos_n=read.csv("benthos_beta_P.csv", head=T, row.names = 1)
benthos.bd.part <- beta.div.comp(benthos_n, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
benthos.bd.part
summary(benthos.bd.part$D)# summary of Dissimilarity matrix

benthos.bd.part1 <- beta.div.comp(benthos_n, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
#coef="S","BJ": Baselga family, Jaccard-based indices
benthos.bd.part1
summary(benthos.bd.part1$D)# summary of Dissimilarity matrix


######beta-active restoration site

benthos_sn=read.csv("benthos_beta_A.csv", head=T, row.names = 1)
benthos.bd.part <- beta.div.comp(benthos_sn, coef="BS", quant=T) #Decompose D in replacement and richness difference components
#coef="S","S" or "Sorensen": Podani family, S??rensen-based indices.
benthos.bd.part
summary(benthos.bd.part$D)# summary of Dissimilarity matrix

benthos.bd.part1 <- beta.div.comp(benthos_sn, coef="BJ", quant=T) #Decompose D in replacement and richness difference components
#coef="S","BJ": Baselga family, Jaccard-based indices
benthos.bd.part1
summary(benthos.bd.part1$D)# summary of Dissimilarity matrix
