# species x traits PCA and related illustrations


library(SYNCSA)
library(tidyverse)


traits<-read_csv("traits_pca_ready.csv")

Tachet_traits <- traits %>% 
  select(species_id, matches("^[A-Z]{2}\\\\d")) 

trait.full_original <- Tachet_traits

# remove animals which have more than 7 NAs, and check that this is successful
remover <-c(4076, 4461, 4681, 5201, 5206, 5211, 5421, 6066, 6231, 6236, 6531, 6771, 6776, 6781, 6786,
            6791, 6801, 6926, 6936, 6941, 6946, 6951, 7251, 7316, 7886, 8026)

length(remover)

sort(as.numeric(remover))==sort(as.numeric(trait.full_original$species_id[apply(is.na(trait.full_original[,-c(1)]), 1, sum)>7]))

head(trait.full_original)
trait.full_complete_spp <- trait.full_original[-which(trait.full_original$species_id %in% remover),]
dim(trait.full_complete_spp)
head(trait.full_complete_spp)

# convert traits into ranks
Ttraits<-trait.full_complete_spp[,-c(1)]
Ttraits<-Ttraits[,]
Ttraits <- sapply(1:ncol(Ttraits), function(i) rank(Ttraits[,i], na.last = "keep")) # transform each trait in ranks
colnames(Ttraits) <- colnames(trait.full_complete_spp)[-c(1)]
rownames(Ttraits) <- rownames(trait.full_complete_spp)

# View(traits)
head(Ttraits)
dim(Ttraits)

#Run PCA, get results to interpret
RES_pca<-SYNCSA::pca(Ttraits)
RES_pca$eigenvalues
RES_pca$eigenvalues[1:4,] %>% view()
head(RES_pca$individuals)
RES_pca$variables %>% view

#identify traits associated positively and negatively with axes
RES_pca$variables %>%
  as.data.frame() %>% 
  select(Axis.1) %>% 
  filter(Axis.1>0.5|Axis.1<(-0.5)) %>% 
  rownames_to_column(var = "traits") %>% 
  view()

RES_pca$variables %>%
  as.data.frame() %>% 
  select(Axis.2) %>% 
  filter(Axis.2>0.5|Axis.2<(-0.5)) %>% 
  rownames_to_column(var = "traits") %>% 
  view()

RES_pca$variables %>%
  as.data.frame() %>% 
  select(Axis.3) %>% 
  filter(Axis.3>0.4|Axis.3<(-0.4)) %>% 
  rownames_to_column(var = "traits") %>% 
  view()

RES_pca$variables %>%
  as.data.frame() %>% 
  select(Axis.4) %>% 
  filter(Axis.4>0.5|Axis.4<(-0.5)) %>% 
  rownames_to_column(var = "traits") %>% 
  view()

#identify species associated with positive and negative ends of axes
minitraits<-traits %>% select(species_id, taxon_name) %>% mutate(species_id = as.character(species_id))
species.pca <- cbind(trait.full_complete_spp$species_id,  RES_pca$individuals[,1:4]) %>% 
  as.data.frame() %>% 
  mutate(species_id = as.character(V1)) %>% 
  left_join(minitraits)

trait.full<-traits

head(trait.full)
trait.full$ord

taxa.name<-trait.full$ord
taxa.name<-ifelse(is.na(taxa.name), trait.full$subclass, taxa.name)
taxa.name<-ifelse(is.na(taxa.name), trait.full$class, taxa.name)
taxa.name<-ifelse(is.na(taxa.name), trait.full$phylum, taxa.name)

taxa.name[which(trait.full$family=="Chironomidae")] <- "Chironomidae"
taxa.name[which(trait.full$family=="Culicidae")] <- "Culicidae"
taxa.name[which(taxa.name=="Diptera")] <- "Other_Diptera"

taxa.name <- ifelse(taxa.name == "Trombidiformes", "Acari", taxa.name)
taxa.name <- ifelse(taxa.name == "Rhynchobdellida", "Hirudinea", taxa.name)
taxa.name <- ifelse(taxa.name == "Turbellaria", "Platyhelminthes", taxa.name)
taxa.name <- ifelse(taxa.name == "Harpacticoida", "Copepoda", taxa.name)
taxa.name <- ifelse(taxa.name == "Cyclopoida", "Copepoda", taxa.name)
taxa.name <- ifelse(taxa.name == "Haplotaxida", "Oligochaeta", taxa.name)
taxa.name <- ifelse(taxa.name == "Nematomorpha", "Nematoda", taxa.name)

position.remove<-which(taxa.name=="Annelida") # remove Annelida
position.remove<-c(position.remove,which(is.na(taxa.name))) # remove unknown taxa 
position.remove

taxa.name<-taxa.name[-1*position.remove] # remove Annelida and unknown
taxa.name
table(taxa.name)

gra_options<-matrix(NA,length(table(taxa.name)),3)
rownames(gra_options)<-names(table(taxa.name))
colnames(gra_options)<-c("pch","col","freq")
gra_options[,3]<-table(taxa.name)

gra_options["Chironomidae",1:2]<-c(19,"black") #ok
gra_options["Coleoptera",1:2]<-c(19,"red") #ok
gra_options["Culicidae",1:2]<-c(19,"grey45")  #ok
gra_options["Ephemeroptera",1:2]<-c(15,"gold3") #ok
gra_options["Hemiptera",1:2]<-c(21,"red") # ok
gra_options["Hirudinea",1:2]<-c(17,"gold3") #ok
gra_options["Lepidoptera",1:2]<-c(19,"blue") #ok
gra_options["Megaloptera",1:2]<-c(19,"gold3") #ok
gra_options["Nematoda",1:2]<-c(17,"red") # ok
gra_options["Odonata",1:2]<-c(24,"blue") #ok
gra_options["Oligochaeta",1:2]<-c(24,"green4") # ok
gra_options["Opisthopora",1:2]<-c(15,"green4") #ok
gra_options["Other_Diptera",1:2]<-c(21,"black") #ok
gra_options["Panpulmonata",1:2]<-c(15,"red") #ok
gra_options["Platyhelminthes",1:2]<-c(17,"green4") #ok
gra_options["Podocopida",1:2]<-c(24,"red") # ok
gra_options["Trichoptera",1:2]<-c(22,"red") # ok

gra_options 
gra_options[order(as.numeric(gra_options[,3] )),]

pch_vec<-as.numeric(gra_options[,1])
names(pch_vec)<-rownames(gra_options)
pch_col<-gra_options[,2]

#====plotting PCA======

pdf("plots_axes_0.7.7_ranks.pdf", 8,12)

par(mfrow=c(4,2), oma = c(1, 1, 2.5, 1),mgp = c(2, 0.5, 0),font.lab=1,cex.lab=1, cex.axis=1)
layout(matrix(c(1:8),4,2,byrow=T), widths=c(1,1,1,1),heights=c(1,1,1,0.4))

axes <- 1:2 # select axes 

ind.scores <-RES_pca$individuals[-1*position.remove,axes] # remove unknown taxa and Annelida
ind.scores
per.inertia <- RES_pca$eigenvalues[axes,2]
per.inertia <- per.inertia*100
x.lim<-c(min(ind.scores[, 1]), max(ind.scores[, 1]))
y.lim<-c(min(ind.scores[, 2]), max(ind.scores[, 2]))
xlim<-c(min(y.lim,x.lim), max(y.lim,x.lim))
ylim<-xlim

par(mar = c(4, 4, 0.5, 0.5))
plot(ind.scores[,1], ind.scores[,2], xlab = paste("PC ", axes[1], " (", round(per.inertia[1], 2), "%)", sep = ""), ylab = paste("PC ", axes[2], " (", round(per.inertia[2], 2), "%)", sep = ""), type = "n", xlim = x.lim, ylim = y.lim, asp = 1, bty = "n", xaxt = "n", yaxt = "n",  mgp = c(1, 1, 0))
points(ind.scores[,1], ind.scores[, 2], pch = pch_vec[taxa.name], col = pch_col[taxa.name])
abline(h = 0, v = 0)

var.scores <- RES_pca$variables[,axes]
var.scores
cor.min<-0.5
var.scores<-var.scores[apply(abs(var.scores)>=cor.min, 1, sum)>0,] # filter correlation
var.scores
circle <- seq(0, 2*pi, length = 200)
par(mar=c(1,1,1,1))
graphics::plot(cos(circle), sin(circle), type = 'l', col = "black", xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 1), asp = 1, bty = "n", xaxt = "n", yaxt = "n")
graphics::abline(h = 0, v = 0, lty = 3, col = "black")
graphics::arrows(0 ,0, x1 = var.scores[,1]*1, y1 = var.scores[,2]*1, length = 0.05, col = "black", lwd =0.7)

# draw each label separately (adjust if necessary)
rownames(var.scores)
arrow <- 1
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 2
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 3
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 4
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1)+0.05, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 5
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 6
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 7
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1)-0.03, (var.scores[arrow, 2]*1.1)+0.05, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 8
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1)+0.05, (var.scores[arrow, 2]*1.1)-0.04, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 9
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1)-0.1, (var.scores[arrow, 2]*1.1)+0.02, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 10
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1)+0.04, (var.scores[arrow, 2]*1.1)-0.02, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 11
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 12
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 13
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1)+0.04, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 14
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1)-0.03, (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 15
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 16
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1)-0.06, (var.scores[arrow, 2]*1.1)-0.04, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 17
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 18
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 19
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 20
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1)-0.07, (var.scores[arrow, 2]*1.1)-0.04, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 21
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 22
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 23
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1)-0.05, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)

axes <- c(1,3) # select axes 
ind.scores <-RES_pca$individuals[-1*position.remove,axes] # remove unknown taxa and Annelida
ind.scores
per.inertia <- RES_pca$eigenvalues[axes,2]
per.inertia <- per.inertia*100
x.lim<-c(min(ind.scores[, 1]), max(ind.scores[, 1]))
y.lim<-c(min(ind.scores[, 2]), max(ind.scores[, 2]))
xlim<-c(min(y.lim,x.lim), max(y.lim,x.lim))
ylim<-xlim

par(mar = c(4, 4, 0.5, 0.5))
plot(ind.scores[,1], ind.scores[,2], xlab = paste("PC ", axes[1], " (", round(per.inertia[1], 2), "%)", sep = ""), ylab = paste("PC ", axes[2], " (", round(per.inertia[2], 2), "%)", sep = ""), type = "n", xlim = x.lim, ylim = y.lim, asp = 1, bty = "n", xaxt = "n", yaxt = "n",  mgp = c(1, 1, 0))
points(ind.scores[,1], ind.scores[, 2], pch = pch_vec[taxa.name], col = pch_col[taxa.name])
abline(h = 0, v = 0)

var.scores <- RES_pca$variables[,axes]
var.scores
cor.min<-0.5
var.scores<-var.scores[apply(abs(var.scores)>=cor.min, 1, sum)>0,] # filter correlation
var.scores

circle <- seq(0, 2*pi, length = 200)
par(mar=c(1,1,1,1))

graphics::plot(cos(circle), sin(circle), type = 'l', col = "black", xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 1), asp = 1, bty = "n", xaxt = "n", yaxt = "n")
graphics::abline(h = 0, v = 0, lty = 3, col = "black")
graphics::arrows(0 ,0, x1 = var.scores[,1]*1, y1 = var.scores[,2]*1, length = 0.05, col = "black", lwd =0.7)

# draw each label separately (adjust if necessary)
rownames(var.scores)
arrow <- 1
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 2
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 3
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 4
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 5
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1)-0.04, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 6
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1)+0.02, (var.scores[arrow, 2]*1.1)-0.02, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 7
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1)+0.04, (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 8
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 9
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 10
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 11
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 12
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 13
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 14
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 15
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1)-0.01, (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 16
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 17
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 18
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)


axes <- c(1,4) # select axes 

ind.scores <-RES_pca$individuals[-1*position.remove,axes] # remove unknown taxa and Annelida
ind.scores
per.inertia <- RES_pca$eigenvalues[axes,2]
per.inertia <- per.inertia*100
x.lim<-c(min(ind.scores[, 1]), max(ind.scores[, 1]))
y.lim<-c(min(ind.scores[, 2]), max(ind.scores[, 2]))
xlim<-c(min(y.lim,x.lim), max(y.lim,x.lim))
ylim<-xlim

par(mar = c(4, 4, 0.5, 0.5))
plot(ind.scores[,1], ind.scores[,2], xlab = paste("PC ", axes[1], " (", round(per.inertia[1], 2), "%)", sep = ""), ylab = paste("PC ", axes[2], " (", round(per.inertia[2], 2), "%)", sep = ""), type = "n", xlim = x.lim, ylim = y.lim, asp = 1, bty = "n", xaxt = "n", yaxt = "n",  mgp = c(1, 1, 0))
points(ind.scores[,1], ind.scores[, 2], pch = pch_vec[taxa.name], col = pch_col[taxa.name])
abline(h = 0, v = 0)

var.scores <- RES_pca$variables[,axes]
var.scores
cor.min<-0.5

var.scores<-var.scores[apply(abs(var.scores)>=cor.min, 1, sum)>0,] # filter correlation
var.scores

par(mar=c(1,1,1,1))
circle <- seq(0, 2*pi, length = 200)
graphics::plot(cos(circle), sin(circle), type = 'l', col = "black", xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 1), asp = 1, bty = "n", xaxt = "n", yaxt = "n")
graphics::abline(h = 0, v = 0, lty = 3, col = "black")
graphics::arrows(0 ,0, x1 = var.scores[,1]*1, y1 = var.scores[,2]*1, length = 0.05, col = "black", lwd =0.7)

# draw each label separately (adjust if necessary)
rownames(var.scores)
arrow <- 1
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 2
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 3
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 4
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 5
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 6
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 7
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 8
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 9
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 10
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 11
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 12
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1)-0.04, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 13
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 14
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1)-0.06, (var.scores[arrow, 2]*1.1)-0.06, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 15
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1)+0.03, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 16
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 17
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 18
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1), (var.scores[arrow, 2]*1.1), labels = paste(rownames(var.scores)[arrow]), cex = 0.8)
arrow <- 19
rownames(var.scores)[arrow]
text((var.scores[arrow,1]*1.1)+0.05, (var.scores[arrow, 2]*1.1)-0.05, labels = paste(rownames(var.scores)[arrow]), cex = 0.8)

par(mar=c(0,0,0,0))

legend_order<-order(as.numeric(gra_options[,3]),decreasing = TRUE)
gra_options[legend_order,]
legend_order<-c(legend_order[2:3],legend_order[1],legend_order[4:length(legend_order)]) # organize lengend for show Culicidae, Chironomidae and Other Diptera first all
legend_order

rownames(gra_options)[which(rownames(gra_options)=="Other_Diptera")] <- "Other Diptera"

plot(1,1, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
legend("topright",legend = rownames(gra_options)[legend_order[1:8]],pch = pch_vec[legend_order[1:8]],col = pch_col[legend_order[1:8]], bty = "n", ncol = 2)

plot(1,1, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
legend("topleft",legend = rownames(gra_options)[legend_order[9:16]],pch = pch_vec[legend_order[9:16]],col = pch_col[legend_order[9:16]], bty = "n", ncol = 2)

dev.off()

write.table(RES_pca$variables,sep="\\t", file="RES_pca_variables_0.7.7_ranks.txt")
write.table(cbind(species_id=trait.full_complete_spp$species_id,RES_pca$individuals),sep="\\t", file="RES_pca_individuals_0.7.7_ranks_cleaned.txt", row.names=FALSE)
write.table(RES_pca$eigenvalues,sep="\\t", file="eigenvalues_0.7.7_ranks.txt")


