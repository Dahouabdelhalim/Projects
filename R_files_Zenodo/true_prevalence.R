setwd("C:/Users/llf44/OneDrive/Documents/ACADEMIA/Graduate School/R Files/Cornell/bee_pathogens/2015-FINAL/science/verified")
pathogens<-read.csv("pathogens_2015R.csv")
flower<-read.csv("flowers.csv")

library(epiR)
library(car)

#for how many bee species do we have at least one pathogen amplifying?
x<-as.data.frame.matrix(table(pathogens$Gen.sp,pathogens$pathogen))
names(x)<-c("absent","present")
x$species<-row.names(x)
x<-droplevels(subset(x, !species %in% c("Unknown.unknown","Melissodes.Melissodes.sp","Megachile.Megachile.sp","Lasioglossum.Lasioglossum.sp","Hylaeus.Hylaeus.sp","Hoplitis.Hoplitis.sp","Ceratina.Ceratina.sp","Andrena.Andrena.sp")))
y<-droplevels(subset(x, !present==0))
((30)/46)*100 # percent species that have at least one pathogen

#for how many flower species do we have at least one pathogen amplifying? 
flower$total<-flower$C.b+flower$NosC+flower$Neog
flower$pathogen <- pmin(flower$total, 1)
table(flower$flower.sp,flower$pathogen) #9/12 = 75%

# Is the "host control" primer, which was developed for Apidae, biased by family?
table(pathogens$family,pathogens$ap)
contable <- array(c(0,33,2,15,2,4,312,28,165,12), dim=c(5,2)) # no differences across families
chisq.test(contable)

## In Bees
#prevalence of Neogregarines
table(pathogens$neo)
epi.prev(pos=34, tested=575, se=.95, sp=.95, method = "blaker", units = 100, conf.level = 0.95)$tp

#prevalence of N.ceranae
table(pathogens$N.c)
epi.prev(pos=90, tested=575, se=.95, sp=.95, method = "blaker", units = 100, conf.level = 0.95)$tp

#prevalence of Trypanosomes
table(pathogens$tryp)
epi.prev(pos=154, tested=575, se=.95, sp=.95, method = "blaker", units = 100, conf.level = 0.95)$tp

#prevalence of N.bombi
table(pathogens$N.b)
epi.prev(pos=5, tested=575, se=.95, sp=.95, method = "blaker", units = 100, conf.level = 0.95)$tp

#prevalence of any pathogen
table(pathogens$pathogen)
epi.prev(pos=233, tested=575, se=.95, sp=.95, method = "blaker", units = 100, conf.level = 0.95)$tp

##In Flowers
#prevalence of Neogregarines
table(flower$Neog)
epi.prev(pos=1, tested=81, se=.95, sp=.95, method = "blaker", units = 100, conf.level = 0.95)$tp

#prevalence of N.ceranae
table(flower$NosC)
epi.prev(pos=15, tested=81, se=.95, sp=.95, method = "blaker", units = 100, conf.level = 0.95)$tp

#prevalence of Trypanosomes
table(flower$C.b)
epi.prev(pos=14, tested=81, se=.95, sp=.95, method = "blaker", units = 100, conf.level = 0.95)$tp

#prevalence of any pathogen
table(flower$pathogen)
epi.prev(pos=27, tested=81, se=.95, sp=.95, method = "blaker", units = 100, conf.level = 0.95)$tp


#prevalence of N.bombi: 0