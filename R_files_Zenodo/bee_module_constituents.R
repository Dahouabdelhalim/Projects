setwd("C:/Users/llf44/OneDrive/Documents/ACADEMIA/Graduate School/R Files/Cornell/bee_pathogens/2015-FINAL/science/verified")

pathogensR<-read.csv("pathogens_2015R.csv")

library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)

pathogensR$Gen.sp <- paste(pathogensR$genus, pathogensR$species, sep='.')
pathogensR$Gen.sp[pathogensR$Gen.sp == "Lasioglossum.Lasioglossum.sp"] <- "Lasioglossum.sp"
pathogensR$Gen.sp[pathogensR$Gen.sp == "Andrena.Andrena.sp"] <- "Andrena.sp"
pathogensR$Gen.sp[pathogensR$Gen.sp == "Hoplitis.Hoplitis.sp"] <- "Hoplitis.sp"
pathogensR$Gen.sp[pathogensR$Gen.sp == "Hylaeus.Hylaeus.sp"] <- "Hylaeus.sp"
pathogensR$Gen.sp[pathogensR$Gen.sp == "Lasioglossum.Lasioglossum.sp"] <- "Lasioglossum.sp"
pathogensR$Gen.sp[pathogensR$Gen.sp == "Melissodes.Melissodes.sp"] <- "Melissodes.sp"
pathogensR$Gen.sp[pathogensR$Gen.sp == "Unknown.unknown"] <- "Unknown"
pathogensR$Gen.sp[pathogensR$Gen.sp == "Ceratina.Ceratina.sp"] <- "Ceratina.sp"
pathogensR$Gen.sp[pathogensR$Gen.sp == "Megachile.sp.Megachile.sp"] <- "Megachile.sp"
pathogensR$Gen.sp[pathogensR$Gen.sp == "Megachile.Megachile.sp"] <- "Megachile.sp"

#first will calculate pathogen prevalence in each module
PathSum <- pathogensR %>%
  group_by(mod.group) %>%
  summarise_at(.vars = vars(tryp:pathogen), .funs = funs(sum(.)))

PathCount <- pathogensR %>%
  group_by(mod.group) %>%
  summarise(count= n())

CombinedPath<-inner_join(PathCount, PathSum, 'mod.group') %>%
  group_by(mod.group) %>%
  mutate(Prev=pathogen/count)

#now will comine that information with species level summary
SpMod<- pathogensR %>%
  group_by(Gen.sp,mod.group) %>%
  summarise(count=n())

dat<-merge(SpMod,CombinedPath[,c(1,8)],by="mod.group")
data_wide <-as.data.frame(t( spread(dat, Gen.sp, count)))

#write.csv(data_wide,"module.constituents.csv")
