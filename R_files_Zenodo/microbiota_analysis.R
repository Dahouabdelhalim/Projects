##Microbiota analysis

require("ape")
require("vegan")
require("reshape2")

#Import distance matrix of fish guts & water from Apoyo, Xiloa, Nicaragua & Managua

#Either weighted/unweighted unifrac or bray-curtis data

distmat_complete <- read.csv("distance-matrix_wuf.tsv", dec = ",", sep = "\\t", as.is = T)

#distmat_complete <- read.csv("distance-matrix_wuf.tsv", dec = ",", sep = "\\t", as.is = T)
#distmat_complete <- read.csv("distance-matrix_uwuf.tsv", dec = ",", sep = "\\t", as.is = T)
#distmat_complete <- read.csv("distance-matrix_bc.tsv", dec = ",", sep = "\\t", as.is = T)

rownames(distmat_complete) <- distmat_complete[,1]
distmat_complete <- distmat_complete[,-1]
names_vec <- rownames(distmat_complete)
colnames(distmat_complete) <- names_vec

#Import metadata file

ids_complete <- read.csv("complete-samples.txt", dec = ",", sep = "\\t", as.is = T)

#ids_complete <- read.csv("samples.txt", dec = ",", sep = "\\t", as.is = T)

ids_complete$pcolor[ids_complete$species=="cit_nic"] <- "limegreen"
ids_complete$pcolor[ids_complete$species=="cit_man"] <- "darkgreen"
ids_complete$pcolor[ids_complete$species=="amarillo"] <- "mediumblue"
ids_complete$pcolor[ids_complete$species=="xiloaensis"] <- "dodgerblue3"
ids_complete$pcolor[ids_complete$species=="sagittae"] <- "skyblue2"
ids_complete$pcolor[ids_complete$species=="astorquii"] <- "darkred"
ids_complete$pcolor[ids_complete$species=="chancho"] <- "firebrick3"
ids_complete$pcolor[ids_complete$species=="globosus"] <- "indianred2"
ids_complete$pcolor[ids_complete$species=="zaliosus"] <- "darkorange1"
ids_complete$pcolor[ids_complete$species=="water"] <- "black"

ids_complete$lcolor[ids_complete$lake=="nicaragua"] <- "darkgreen"
ids_complete$lcolor[ids_complete$lake=="managua"] <- "blue"
ids_complete$lcolor[ids_complete$lake=="xiloa"] <- "purple"
ids_complete$lcolor[ids_complete$lake=="apoyo"] <- "orange"

ids_complete$pch[ids_complete$lake=="nicaragua"] <- "16"
ids_complete$pch[ids_complete$lake=="managua"] <- "18"
ids_complete$pch[ids_complete$lake=="xiloa"] <- "15"
ids_complete$pch[ids_complete$lake=="apoyo"] <- "17"

ids_complete$pch <- as.numeric(ids_complete$pch)

pcoa_complete <- pcoa(distmat_complete)

pcoa_axes <- pcoa_complete$vectors
pcoa_axes <- as.data.frame(pcoa_axes)

pcoa_ids <- cbind(ids_complete, pcoa_axes)

plot(pcoa_axes$Axis.2 ~ pcoa_axes$Axis.1, pch = pcoa_ids$pch, col = pcoa_ids$pcolor)

#Test whether bacterial community differs between lake water and fish guts
adonis(distmat_complete ~ type, data = pcoa_ids)

#Taxa barplots

taxa <- read.csv("bacterial_phyla.csv", dec = ",", sep = ";", as.is = T)

rownames(taxa) <- taxa[,1]
taxa <- taxa[,-1]

barplot(as.matrix(taxa), legend = rownames(taxa))

#Alpha diversity

adiv <- read.csv("alpha_div.csv", dec = ",", sep = ";", as.is = T)
adiv$type <- as.factor(adiv$type)
adiv$habitat <- as.factor(adiv$habitat)
adiv$species <- factor(adiv$species,levels=c("water_nic","water_man","water_apo", "water_xil","cit_nic", "cit_man", "astorquii", "chancho", "globosus", "zaliosus", "amarillo", "xiloaensis", "sagittae"))

plot(adiv$observed_otus ~ adiv$species)

#Test whether alpha diversity differs between lake water and fish guts
wilcox.test(adiv$observed_otus ~ adiv$type)

adiv_water <- adiv[1:16,]
adiv_fish <- adiv[17:146,]
adiv_fish_nomanagua <- adiv[c(17:26,37:146),]

#Test whether alpha diversity differs between great lakes and crater lakes for lake water and fish guts
wilcox.test(adiv_water$observed_otus ~ adiv_water$habitat)
wilcox.test(adiv_fish$observed_otus ~ adiv_fish$habitat)
wilcox.test(adiv_fish_nomanagua$observed_otus ~ adiv_fish_nomanagua$habitat)

#Test whether alpha diversity differs between species in crater lakes
#Change species names accordingly
wilcox.test(adiv_fish$observed_otus[which(adiv_fish$species == "amarillo")], adiv_fish$observed_otus[which(adiv_fish$species == "xiloaensis")])

#Subset: only fish guts

distmat_fish <- distmat_complete[1:130,1:130]
ids_fish <- ids_complete[1:130,]

plot(ids_fish$n ~ ids_fish$c, pch = ids_fish$pch, col = ids_fish$pcolor)

adonis(distmat_fish ~ habitat*lake, data = ids_fish, strate = A:B)

##SIA data

#All locations

#Remove fish with missing SIA data
fish_sia <- subset(ids_fish, !is.na(ids_fish$n))
fish_sia$lake <- as.factor(fish_sia$lake)
fish_sia$habitat <- as.factor(fish_sia$habitat)

#Differences in Carbon/nitrogen among lakes
kruskal.test(fish_sia$c ~ fish_sia$lake)
kruskal.test(fish_sia$n ~ fish_sia$lake)

#Subset: only crater lake fish

distmat_crater <- distmat_complete[c(2:14,16,17,19:21,23:57,59,61:65,67:90,92:95,98,100:103,107,110,111,114:116,118:127,129,130),c(2:14,16,17,19:21,23:57,59,61:65,67:90,92:95,98,100:103,107,110,111,114:116,118:127,129,130)]
ids_crater <- ids_complete[c(2:14,16,17,19:21,23:57,59,61:65,67:90,92:95,98,100:103,107,110,111,114:116,118:127,129,130),]

#Subset: only crater lake fish with SIA data
distmat_crater_sia <- distmat_crater[c(1:11,13:66,68,70,71,73:80,82:94,96:98,100:110),c(1:11,13:66,68,70,71,73:80,82:94,96:98,100:110)]
ids_crater_sia <- ids_crater[c(1:11,13:66,68,70,71,73:80,82:94,96:98,100:110),]

pcoa_crater <- pcoa(distmat_crater)

pcoa_axes_crater <- pcoa_crater$vectors
pcoa_axes_crater <- as.data.frame(pcoa_axes_crater)

pcoa_ids_crater <- cbind(ids_crater, pcoa_axes_crater)

plot(pcoa_axes_crater$Axis.2 ~ pcoa_axes_crater$Axis.1, pch = pcoa_ids_crater$pch, col = pcoa_ids_crater$pcolor)

#Xiloa distmat

distmat_xiloa <- distmat_crater[c(1,5,8,10:12,14,17,24,26,27,29:31,33,35:37,40:45,47,48,50:52,54:56,58,60:63,65,67:69,71,72,75:77,80:82,89,91,95,98,101,106,110),c(1,5,8,10:12,14,17,24,26,27,29:31,33,35:37,40:45,47,48,50:52,54:56,58,60:63,65,67:69,71,72,75:77,80:82,89,91,95,98,101,106,110)]
ids_xiloa <- ids_crater[c(1,5,8,10:12,14,17,24,26,27,29:31,33,35:37,40:45,47,48,50:52,54:56,58,60:63,65,67:69,71,72,75:77,80:82,89,91,95,98,101,106,110),]

adonis(distmat_xiloa ~ species, data = ids_xiloa)

#Pairwise distances between individuals

D <- as.matrix(distmat_xiloa)
D[upper.tri(D)] <- NA
x <- subset(melt(D), value!=0)
x <- as.data.frame(x)

row.names(x) <- NULL

#Nitrogen
n_dist <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  dif <- abs(ids_xiloa$n[which(ids_xiloa$id==id1)]-ids_xiloa$n[which(ids_xiloa$id==id2)])
  n_dist[g] <- dif
}

x$n_dist <- n_dist 

plot(value~n_dist, data = x)
abline(lm(value~n_dist, data = x))
cor.test(x$value, x$n_dist, method = "pearson")

#Carbon
c_dist <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  dif <- abs(ids_xiloa$c[which(ids_xiloa$id==id1)]-ids_xiloa$c[which(ids_xiloa$id==id2)])
  c_dist[g] <- dif
}

x$c_dist <- c_dist 

plot(value~c_dist, data = x)
abline(lm(value~c_dist, data = x))
cor.test(x$value, x$c_dist, method = "pearson")

#Species comparisons
#Test for nitrogen and carbon separately, adjust code accordingly
eco <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  if  (ids_xiloa$species[which(ids_xiloa$id==id1)]==ids_xiloa$species[which(ids_xiloa$id==id2)]) {
    eco[g] <- "intra"
  }
  else {
    eco[g] <- "inter"
  }
}

x$eco <- eco

sp1 <- NA
sp2 <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  sp1[g] <- (ids_xiloa$species[which(ids_xiloa$id==id1)])
  sp2[g] <- (ids_xiloa$species[which(ids_xiloa$id==id2)])
}

x$sp1 <- sp1
x$sp2 <- sp2

dif <- abs(ids_xiloa$n[which(ids_xiloa$id==id1)]-ids_xiloa$n[which(ids_xiloa$id==id2)])

x_intra <- x[which(x$eco=="intra"),]
x_inter <- x[which(x$eco=="inter"),]

x <- rbind(x_intra, x_inter)
x$eco <- as.factor(x$eco)

boxplot(x$n_dist ~ x$eco)
kruskal.test(x$n_dist ~ x$eco)

boxplot(x$value ~ x$eco)
kruskal.test(x$value ~ x$eco)

#Apoyo distmat

distmat_apoyo <- distmat_crater[c(2:4,6,7,9,13,15,16,18:23,25,28,32,34,38,39,46,49,53,57,59,64,66,70,73,74,78,79,83:88,90,92:94,96,97,99,100,102:105,107:109),c(2:4,6,7,9,13,15,16,18:23,25,28,32,34,38,39,46,49,53,57,59,64,66,70,73,74,78,79,83:88,90,92:94,96,97,99,100,102:105,107:109)]
ids_apoyo <- ids_crater[c(2:4,6,7,9,13,15,16,18:23,25,28,32,34,38,39,46,49,53,57,59,64,66,70,73,74,78,79,83:88,90,92:94,96,97,99,100,102:105,107:109),]

adonis(distmat_apoyo ~ species, data = ids_apoyo)

#Pairwise distances between individuals

D <- as.matrix(distmat_apoyo)
D[upper.tri(D)] <- NA
x <- subset(melt(D), value!=0)
x <- as.data.frame(x)

row.names(x) <- NULL

#Nitrogen
n_dist <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  dif <- abs(ids_apoyo$n[which(ids_apoyo$id==id1)]-ids_apoyo$n[which(ids_apoyo$id==id2)])
  n_dist[g] <- dif
}

x$n_dist <- n_dist 

plot(value~n_dist, data = x)
abline(lm(value~n_dist, data = x))
cor.test(x$value, x$n_dist, method = "pearson")

#carbon
c_dist <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  dif <- abs(ids_apoyo$c[which(ids_apoyo$id==id1)]-ids_apoyo$c[which(ids_apoyo$id==id2)])
  c_dist[g] <- dif
}

x$c_dist <- c_dist 

plot(value~c_dist, data = x)
abline(lm(value~c_dist, data = x))
cor.test(x$value, x$c_dist, method = "pearson")

#Species comparisons
#Test for nitrogen and carbon separately, adjust code accordingly
eco <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  if  (ids_apoyo$species[which(ids_apoyo$id==id1)]==ids_apoyo$species[which(ids_apoyo$id==id2)]) {
    eco[g] <- "intra"
  }
  else {
    eco[g] <- "inter"
  }
}

x$eco <- eco

sp1 <- NA
sp2 <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  sp1[g] <- (ids_apoyo$species[which(ids_apoyo$id==id1)])
  sp2[g] <- (ids_apoyo$species[which(ids_apoyo$id==id2)])
}

x$sp1 <- sp1
x$sp2 <- sp2

dif <- abs(ids_apoyo$n[which(ids_apoyo$id==id1)]-ids_apoyo$n[which(ids_apoyo$id==id2)])

x_intra <- x[which(x$eco=="intra"),]
x_inter <- x[which(x$eco=="inter"),]

x <- rbind(x_intra, x_inter)
x$eco <- as.factor(x$eco)

boxplot(x$n_dist ~ x$eco)
kruskal.test(x$n_dist ~ x$eco)

boxplot(x$value ~ x$eco)
kruskal.test(x$value ~ x$eco)

#Z-normalization function

znorm <- function(ts){
  ts.mean <- mean(ts)
  ts.dev <- sd(ts)
  (ts - ts.mean)/ts.dev
}

ids_apo_sia <- ids_crater_sia[which(ids_crater_sia$lake=="apoyo"),]

ids_apo_sia$n_norm = znorm(ids_apo_sia$n)
ids_apo_sia$c_norm = znorm(ids_apo_sia$c)

#Test for differences in SIA among species within Apoyo
ids_apo_sia$species <- as.factor(ids_apo_sia$species)
kruskal.test(ids_apo_sia$n ~ ids_apo_sia$species)
kruskal.test(ids_apo_sia$c ~ ids_apo_sia$species)

ids_xil_sia <- ids_crater_sia[which(ids_crater_sia$lake=="xiloa"),]

#Test for differences in SIA among species within Xiloa
ids_xil_sia$species <- as.factor(ids_xil_sia$species)
kruskal.test(ids_xil_sia$n ~ ids_xil_sia$species)
kruskal.test(ids_xil_sia$c ~ ids_xil_sia$species)

ids_xil_sia$n_norm = znorm(ids_xil_sia$n)
ids_xil_sia$c_norm = znorm(ids_xil_sia$c)

ids_crater_sia_new <- rbind(ids_xil_sia, ids_apo_sia) 

#Arrange ids_crater_sia_new to match order in distmat_crater_sia
crater_sia_ids <- rownames(distmat_crater_sia) 

ids_crater_sia_new <- ids_crater_sia_new[match(crater_sia_ids, ids_crater_sia_new$id),]

adonis(distmat_crater_sia ~ lake*n_norm, data = ids_crater_sia_new)
adonis(distmat_crater_sia ~ lake*c_norm, data = ids_crater_sia_new)

plot(ids_crater_sia_new$n_norm ~ ids_crater_sia_new$c_norm, col = ids_crater_sia_new$pcolor, pch = ids_crater_sia_new$pch)

#Picrust data

#all fish
pwab <- read.csv("pathwayabundance.csv", sep = ";", dec = ",")
rownames(pwab) <- pwab$id
pwab <- pwab[,-1]
pwab_ids <- pwab[,1:6]

#crater lake subset
pwab_crater <- pwab[c(1:54,76:130),]
pwab_crater_ids <- pwab_crater[,1:6]

pwab_crater_dist <- dist(pwab_crater, method = "euclidean")
pwab_crater_dist <- as.matrix(pwab_crater_dist)
pwab_crater_dist <- as.data.frame(pwab_crater_dist)

#Apoyo subset
pwab_apo <- pwab[1:54,]
pwab_apo_ids <- pwab_apo[,1:6]

pwab_apo_dist <- dist(pwab_apo, method = "euclidean")
pwab_apo_dist <- as.matrix(pwab_apo_dist)
pwab_apo_dist <- as.data.frame(pwab_apo_dist)

adonis(pwab_apo_dist ~ species, data = pwab_apo_ids)

#Xiloa subset
pwab_xil <- pwab[76:130,]
pwab_xil_ids <- pwab_xil[,1:6]

pwab_xil_dist <- dist(pwab_xil, method = "euclidean")
pwab_xil_dist <- as.matrix(pwab_xil_dist)
pwab_xil_dist <- as.data.frame(pwab_xil_dist)

adonis(pwab_xil_dist ~ species, data = pwab_xil_ids)

pwab_xil_ids_sia <- subset(pwab_xil_ids, !is.na(pwab_xil_ids$n))
pwab_apo_ids_sia <- subset(pwab_apo_ids, !is.na(pwab_apo_ids$n))

pwab_xil_ids_sia$n_norm = znorm(pwab_xil_ids_sia$n)
pwab_xil_ids_sia$c_norm = znorm(pwab_xil_ids_sia$c)

pwab_apo_ids_sia$n_norm = znorm(pwab_apo_ids_sia$n)
pwab_apo_ids_sia$c_norm = znorm(pwab_apo_ids_sia$c)

pwab_crater_ids_new <- (rbind(pwab_apo_ids_sia, pwab_xil_ids_sia))

#Remove samples with no SIA data
pwab_crater_dist_sia <- pwab_crater_dist[-c(46,59,92,94,96,101,105), -c(46,59,92,94,96,101,105)]

pwab_crater_new <- pwab_crater[-c(46,59,92,94,96,101,105),]
pwab_crater_new <- cbind(pwab_crater_new, pwab_crater_ids_new$n_norm, pwab_crater_ids_new$c_norm)
colnames(pwab_crater_new)[466] <- "n_norm"
colnames(pwab_crater_new)[467] <- "c_norm"

adonis(pwab_crater_dist_sia ~ lake*n_norm, data = pwab_crater_ids_new)
adonis(pwab_crater_dist_sia ~ lake*c_norm, data = pwab_crater_ids_new)

##Intra- and interspecific distances

#Apoyo

#Pairwise distances between individuals

D <- as.matrix(pwab_apo_dist)
D[upper.tri(D)] <- NA
x <- subset(melt(D), value!=0)
x <- as.data.frame(x)

row.names(x) <- NULL

#nitrogen
n_dist <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  dif <- abs(pwab_apo_ids$n[which(rownames(pwab_apo_ids)==id1)]-pwab_apo_ids$n[which(rownames(pwab_apo_ids)==id2)])
  n_dist[g] <- dif
}

x$n_dist <- n_dist 

#carbon
c_dist <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  dif <- abs(pwab_apo_ids$c[which(rownames(pwab_apo_ids)==id1)]-pwab_apo_ids$c[which(rownames(pwab_apo_ids)==id2)])
  c_dist[g] <- dif
}

x$c_dist <- c_dist 

#Species comparisons

eco <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  if  (pwab_apo_ids$species[which(rownames(pwab_apo_ids)==id1)]==pwab_apo_ids$species[which(rownames(pwab_apo_ids)==id2)]) {
    eco[g] <- "intra"
  }
  else {
    eco[g] <- "inter"
  }
}

x$eco <- eco

sp1 <- NA
sp2 <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  sp1[g] <- pwab_apo_ids$species[which(rownames(pwab_apo_ids)==id1)]
  sp2[g] <- pwab_apo_ids$species[which(rownames(pwab_apo_ids)==id2)]
}

x$sp1 <- sp1
x$sp2 <- sp2

dif <- abs(pwab_apo_ids$n[which(rownames(pwab_apo_ids)==id1)]-pwab_apo_ids$n[which(rownames(pwab_apo_ids)==id2)])

x_intra <- x[which(x$eco=="intra"),]
x_inter <- x[which(x$eco=="inter"),]

x <- rbind(x_intra, x_inter)
x$eco <- as.factor(x$eco)

boxplot(x$value ~ x$eco)
kruskal.test(x$value ~ x$eco)

#Xiloa

#Pairwise distances between individuals

D <- as.matrix(pwab_xil_dist)
D[upper.tri(D)] <- NA
x <- subset(melt(D), value!=0)
x <- as.data.frame(x)

row.names(x) <- NULL

#nitrogen
n_dist <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  dif <- abs(pwab_xil_ids$n[which(rownames(pwab_xil_ids)==id1)]-pwab_xil_ids$n[which(rownames(pwab_xil_ids)==id2)])
  n_dist[g] <- dif
}

x$n_dist <- n_dist 

#carbon
c_dist <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  dif <- abs(pwab_xil_ids$c[which(rownames(pwab_xil_ids)==id1)]-pwab_xil_ids$c[which(rownames(pwab_xil_ids)==id2)])
  c_dist[g] <- dif
}

x$c_dist <- c_dist 

#Species comparisons

eco <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  if  (pwab_xil_ids$species[which(rownames(pwab_xil_ids)==id1)]==pwab_xil_ids$species[which(rownames(pwab_xil_ids)==id2)]) {
    eco[g] <- "intra"
  }
  else {
    eco[g] <- "inter"
  }
}

x$eco <- eco

sp1 <- NA
sp2 <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  sp1[g] <- pwab_xil_ids$species[which(rownames(pwab_xil_ids)==id1)]
  sp2[g] <- pwab_xil_ids$species[which(rownames(pwab_xil_ids)==id2)]
}

x$sp1 <- sp1
x$sp2 <- sp2

dif <- abs(pwab_xil_ids$n[which(rownames(pwab_xil_ids)==id1)]-pwab_xil_ids$n[which(rownames(pwab_xil_ids)==id2)])

x_intra <- x[which(x$eco=="intra"),]
x_inter <- x[which(x$eco=="inter"),]

x <- rbind(x_intra, x_inter)
x$eco <- as.factor(x$eco)

boxplot(x$value ~ x$eco)
kruskal.test(x$value ~ x$eco)

##Parallelism of OTUs and pathway abundance across crater lakes

#OTU abundance

#Import otu abundance file

otu <- read.csv("otu_table_fish_prop.csv", dec = ",", sep = ";")

rownames(otu) <- otu[,1]
otu <- otu[,-1]

otu_data <- otu[,9:484]
otu_ids <- otu[,1:8]
otu_data_filtered <- otu_data[which(otu_data[131,]>0.0001)]

otu_new <- cbind(otu_ids,otu_data_filtered) 
otu_new <- otu_new[1:130,]
otu_new$habitat <- as.factor(otu_new$habitat)
otu_new$ecotype <- as.factor(otu_new$ecotype)
otu_new$species <- as.factor(otu_new$species)

#Subset of crater lakes separately
#Apoyo
otu_apo <- otu_new[which(otu_new$lake=="apoyo"),]
otu_xil <- otu_new[which(otu_new$lake=="xiloa"),]

otu_apo_sia <- subset(otu_apo, !is.na(otu_apo$n))
otu_xil_sia <- subset(otu_xil, !is.na(otu_xil$n))

otu_apo_sia$n_norm = znorm(otu_apo_sia$n)
otu_apo_sia$c_norm = znorm(otu_apo_sia$c)

otu_xil_sia$n_norm = znorm(otu_xil_sia$n)
otu_xil_sia$c_norm = znorm(otu_xil_sia$c)

otu_crater_new <- rbind(otu_apo_sia, otu_xil_sia)

#Test whether OTU abundance is associated with SIA, both crater lakes combined

#linear model for nitrogen

p_values_n_norm <- NA
p_values_int <- NA

for (g in 1:323){
  lm_crater <- lm(otu_crater_new[,g+8] ~ otu_crater_new$lake*otu_crater_new$n_norm)
  summary <- summary(lm_crater)
  p_values_n_norm[g] <- summary$coefficients[3,4]
  p_values_int[g] <- summary$coefficients[4,4]
}

p_fdr_n_norm <- p.adjust(p_values_n_norm, method = "fdr", n = 323)
p_fdr_int <- p.adjust(p_values_int, method = "fdr", n = 323)

otu_crater_new[104,9:331] <- p_values_n_norm
otu_crater_new[105,9:331] <- p_fdr_n_norm
otu_crater_new[106,9:331] <- p_values_int
otu_crater_new[107,9:331] <- p_fdr_int

otu_crater_parallel_n <- otu_crater_new[which(otu_crater_new[104,]<0.05 & otu_crater_new[106,]>0.05)]

length(colnames(otu_crater_parallel_n))

#linear model for carbon

p_values_c_norm <- NA
p_values_int <- NA

for (g in 1:323){
  lm_crater <- lm(otu_crater_new[,g+8] ~ otu_crater_new$lake*otu_crater_new$c_norm)
  summary <- summary(lm_crater)
  p_values_c_norm[g] <- summary$coefficients[3,4]
  p_values_int[g] <- summary$coefficients[4,4]
}

p_fdr_c_norm <- p.adjust(p_values_c_norm, method = "fdr", n = 323)
p_fdr_int <- p.adjust(p_values_int, method = "fdr", n = 323)

otu_crater_new[104,9:331] <- p_values_c_norm
otu_crater_new[105,9:331] <- p_fdr_c_norm
otu_crater_new[106,9:331] <- p_values_int
otu_crater_new[107,9:331] <- p_fdr_int

otu_crater_parallel_c <- otu_crater_new[which(otu_crater_new[104,]<0.05 & otu_crater_new[106,]>0.05)]

length(colnames(otu_crater_parallel_c))

##metabolic pathway abundance

#Test whether metabolic pathways are associated with SIA, both crater lakes combined

#linear model

p_values_n_norm <- NA
p_values_int <- NA

for (g in 1:459){
  lm_crater <- lm(pwab_crater_new[,g+6] ~ pwab_crater_new$lake*pwab_crater_new$n_norm)
  summary <- summary(lm_crater)
  p_values_n_norm[g] <- summary$coefficients[3,4]
  p_values_int[g] <- summary$coefficients[4,4]
}

p_fdr_n_norm <- p.adjust(p_values_n_norm, method = "fdr", n = 459)
p_fdr_int <- p.adjust(p_values_int, method = "fdr", n = 459)

pwab_crater_new[103,7:465] <- p_values_n_norm
pwab_crater_new[104,7:465] <- p_fdr_n_norm
pwab_crater_new[105,7:465] <- p_values_int
pwab_crater_new[106,7:465] <- p_fdr_int

pwab_crater_parallel_n <- pwab_crater_new[which(pwab_crater_new[103,]<0.05 & pwab_crater_new[105,]>0.05)]

length(colnames(pwab_crater_parallel_n))

#linear model

p_values_c_norm <- NA
p_values_int <- NA

for (g in 1:459){
  lm_crater <- lm(pwab_crater_new[,g+6] ~ pwab_crater_new$lake*pwab_crater_new$c_norm)
  summary <- summary(lm_crater)
  p_values_c_norm[g] <- summary$coefficients[3,4]
  p_values_int[g] <- summary$coefficients[4,4]
}

p_fdr_c_norm <- p.adjust(p_values_c_norm, method = "fdr", n = 459)
p_fdr_int <- p.adjust(p_values_int, method = "fdr", n = 459)

pwab_crater_new[103,7:465] <- p_values_c_norm
pwab_crater_new[104,7:465] <- p_fdr_n_norm
pwab_crater_new[105,7:465] <- p_values_int
pwab_crater_new[106,7:465] <- p_fdr_int

pwab_crater_parallel_c <- pwab_crater_new[which(pwab_crater_new[103,]<0.05 & pwab_crater_new[105,]>0.05)]

length(colnames(pwab_crater_parallel_c))
