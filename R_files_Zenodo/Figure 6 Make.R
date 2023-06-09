## R code to produce Figure 6 in:
## Johnston EC, Cunning R, Burgess SC. Cophylogeny and specificity between cryptic coral species (Pocillopora spp.) at Moâ€™orea and their symbionts (Symbiodiniaceae).
## Code written by Erika Johnston. February 2022. Send comments or corrections to ejohnston@bio.fsu.edu
## R version 4.1.2 (2022-03-01)

# initialization ----------------------------------------------------------
rm(list=ls()) #clear all variables

#setwd("")

library(ggplot2)
library(ggstance)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(tidyr)
library(vegan)
library(BiodiversityR)
library(adonis)
library(devtools)

# Get pairwiseAdonis
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# Import data
dat <- read.csv("Figure 4 and 6 data.csv")

# Plot based on majority ITS2 type profile
dat <- dat %>%
  filter(Type_profile_Prop > 0.43)

# Remove columns that are all zeros
dat <- dat[, colSums(dat != 0) > 0]

----------#### PCoA analysis for all Pocillopora species and haplotypes ####--------------
# Calculate distance matrix
# Remove colonies that have C116, D, or A ITS2 type profiles
sym <- dat[
  dat$Coral.ID!="61A58"
  & dat$Coral.ID!="21D10"
  & dat$Coral.ID!="H15"
  & dat$Coral.ID!="22C17"
  & dat$Coral.ID!="G05"
  & dat$Coral.ID!="G44"
  & dat$Coral.ID!="G49"
  & dat$Coral.ID!="C68"
  & dat$Coral.ID!="J17",]

# Create distance matrix
sym_DIV <- sym[, -c(1:6)]
sym_matrix <- vegdist(sym_DIV, method = "bray")

# Running PCoA in vegan (cmdscale)
mod_sym <- cmdscale(sym_matrix, eig=T)

#Extract points
data.scores = as.data.frame(scores(mod_sym))
xy <- data.scores[,1:2]
plot(xy)
mod_sym 

#PERMANOVA
all_test <- adonis(sym_DIV ~ sym$Species.haplotype, method = "bray", permutations = 99999)
all_test

#Pairwise PERMANOVA for all Pocillopora species/haplotypes
pw <- pairwise.adonis2(sym[,7:487]~Species.haplotype, data=sym)
pw

#Order ITS2 type profiles
sym$Type_profile = factor(sym$Type_profile,levels=c(
  "C1.C42.2",
  "C1.C3.C42.2.C1b",
  "C1.C3.C42.2.C1b.C3cg",
  "C42.2.C1.C1b.C3.C1au",
  "C42.2.C1.C42i.C3.C1b.C1au",
  "C1.C42.2.C3.C1b.C1au.C115n",
  "C1.C42.2.C3cg.C1b.C115k.C45c",
  "C1.C42.2.C1b.C1au.C3.C3cg.C1j",
  "C1.C1ag.C1ah.C42.2.C3cg.C1b.C3",
  "C1.C42.2.C1b.C3cg.C3.C1au.C41p",
  "C1.C1m.C1ag.C42.2.C3cg.C1b.C3cw.C1au",
  "C1d.C42.2.C1.C1k.C1b",
  "C1d.C42.2.C1.C3cg.C1b.C1au",
  "C1d.C1.C42.2.C1b.C3cg.C3",
  "C1d.C1.C42.2.C1k.C1b.C3cg.C3cw",
  "C1.C42.2.C1d.C1b.C3.C3cg.C115n.C1au",
  "C1d.C1.C42.2.C3cg.C1b.C45c.C115k.C3.C1au",
  "C1.C1ah.C1ag.C42.2.C3cg.C1b.C1dj.C115k.C45c",
  "C1d.C1.C42.2.C3.C1b.C3cg.C115k.C45c.C1au.C41p",
  "C42a.C42.2.C1.C1b",
  "C42.2.C42a.C1b.C1.C1j",
  "C42a.C42.2.C1b.C1.C42b",
  "C42a.C42.2.C1.C3.C1b.C1au",
  "C42a.C42.2.C1.C1b.C3.C1j.C1au",
  "C42a.C42.2.C1b.C1.C42ao.C1au",
  "C42.2.C42a.C1.C1j.C1b.C1ci.C1au",
  "C42g.C1.C42.2.C42a.C1b.C42h.C3.C1au",
  "C1.C42.2.C42a.C1b.C1az.C1au.C115d.C3.C1j",
  "C42a.C1.C42.2.C1j.C1b.C1au.C3.C115l.C115k",
  "C1.C42.2.C42g.C42a.C1b.C3.C1az.C1au.C42az",
  "C1.C42.2.C42g.C42a.C1b.C1au.C1az.C3.C115k.C115d.C41p",
  "C3.C116.C3dc",
  "C116.C1.C116q",
  "C15.C15ev.C15dt",
  "C116.C116i.C116q",
  "C42.2.C116.C1.C1b",
  "A1.A1fi",
  "D1",
  "D1ch",
  "D1.D6",
  "D1.D6.D4.D2.2.D1h.D2c.D2f"))

#Set ITS2 profile colors
my_colors_profile <- c("C1.C42.2"="#ffedfe",
                       "C1.C3.C42.2.C1b"="#fbccff",
                       "C1.C3.C42.2.C1b.C3cg"="#f2b3ff",
                       "C42.2.C1.C1b.C3.C1au"="#fa7de1",
                       "C42.2.C1.C42i.C3.C1b.C1au"="#d152eb",
                       "C1.C42.2.C3.C1b.C1au.C115n"="#bf4e92",
                       "C1.C42.2.C3cg.C1b.C115k.C45c"="#e35675",
                       "C1.C42.2.C1b.C1au.C3.C3cg.C1j"="#bd2848",
                       "C1.C1ag.C1ah.C42.2.C3cg.C1b.C3"="#8a017a",
                       "C1.C42.2.C1b.C3cg.C3.C1au.C41p"="#4a0475",
                       "C1.C1m.C1ag.C42.2.C3cg.C1b.C3cw.C1au"="#230436",
                       "C1d.C42.2.C1.C1k.C1b"="#fff5b5",
                       "C1d.C42.2.C1.C3cg.C1b.C1au"="#fcf5ca",
                       "C1d.C1.C42.2.C1b.C3cg.C3"="#f5e171",
                       "C1d.C1.C42.2.C1k.C1b.C3cg.C3cw"="#fcbf74",
                       "C1.C42.2.C1d.C1b.C3.C3cg.C115n.C1au"="#ed9f18",
                       "C1d.C1.C42.2.C3cg.C1b.C45c.C115k.C3.C1au"="#de7309",
                       "C1.C1ah.C1ag.C42.2.C3cg.C1b.C1dj.C115k.C45c"="#e85702",
                       "C1d.C1.C42.2.C3.C1b.C3cg.C115k.C45c.C1au.C41p"="#b85707",
                       "C42a.C42.2.C1.C1b"="#c4fffc",
                       "C42.2.C42a.C1b.C1.C1j"="#9cfffa",
                       "C42a.C42.2.C1b.C1.C42b"="#5cf1ff",
                       "C42a.C42.2.C1.C3.C1b.C1au"="#74f1fc",
                       "C42a.C42.2.C1.C1b.C3.C1j.C1au"="#21dced",
                       "C42a.C42.2.C1b.C1.C42ao.C1au"="#259dd9",
                       "C42.2.C42a.C1.C1j.C1b.C1ci.C1au"="#52a0fa",
                       "C42g.C1.C42.2.C42a.C1b.C42h.C3.C1au"="#286bb8",
                       "C1.C42.2.C42a.C1b.C1az.C1au.C115d.C3.C1j"="#0353ad",
                       "C42a.C1.C42.2.C1j.C1b.C1au.C3.C115l.C115k"="#143ec9",
                       "C1.C42.2.C42g.C42a.C1b.C3.C1az.C1au.C42az"="#01068f",
                       "C1.C42.2.C42g.C42a.C1b.C1au.C1az.C3.C115k.C115d.C41p"="#05165c",
                       "C3.C116.C3dc"="#edebeb",
                       "C116.C1.C116q"="#b8b6b6",
                       "C15.C15ev.C15dt"="#adadad",
                       "C116.C116i.C116q"="#858585",
                       "C42.2.C116.C1.C1b"="#5e5d5d",
                       "A1.A1fi"="#ffee00",
                       "D1"="#eafc5d",
                       "D1ch"="#c0d428",
                       "D1.D6"="#a1c91e",
                       "D1.D6.D4.D2.2.D1h.D2c.D2f"="#018a26")



## Plot
# Percent of variance explained by axis
round(mod_sym$eig*100/sum(mod_sym$eig),1)

quartz(width=10,height=6)
par(mfrow=c(1,2),mar=c(5,4,4,0))
plot(xy,xlim=c(-0.25,0.26),ylim=c(-0.45,0.15),type="n", xlab="PCoA1 (30.6%)", ylab="PCoA2 (22.3%)")

# Move ITS2 profile colors into a dataframe
my_colors_profile  <- as.data.frame(my_colors_profile)
my_colors_to_plot <- my_colors_profile[rownames(my_colors_profile) %in% unique(sym$Type_profile),]
my_profiles_to_plot <- rownames(my_colors_profile)[rownames(my_colors_profile) %in% unique(sym$Type_profile)]

# Assign symbols to a species/haplotype.
my_coralID_symbols <- cbind.data.frame(Species.haplotype=unique(sym$Species.haplotype), symbol=c(1,2,3,4,5,6,7))

# Add points, colored by ITS2 Profile, with symbols based on Species.
for(i in 1:length(my_profiles_to_plot)){
  
  # Create a vector of pch symbols that match that desired for each species (or depth, as appropriate)
  tmp <- sym[which(sym$Type_profile==my_profiles_to_plot[i]),]
  pch.vec <- my_coralID_symbols[match(tmp$Species.haplotype,my_coralID_symbols$Species.haplotype),2]
  
  #
  points(x=xy[which(sym$Type_profile==my_profiles_to_plot[i]),1],
         y=xy[which(sym$Type_profile==my_profiles_to_plot[i]),2],
         col = my_colors_to_plot[i],
         pch=pch.vec,
         cex=1)
}
ordiellipse(mod_sym, sym$Species.haplotype, lty = 'dashed', label=F, cex = 0.9, font=2)

# Add legend
par(mar=c(5,1,4,1),xpd=T)
plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",axes=F,xlab="",ylab="")
legend(0,1.05, legend=my_profiles_to_plot, pch=19, col=my_colors_to_plot, cex=0.65, bty="n")
legend(0,0, legend= my_coralID_symbols$Species.haplotype, pch = my_coralID_symbols$symbol, cex=0.65, bty="n")


----------#### Haplotype 10 by depth ####--------------

## Prepare data
# Subset Haplotype 10 and remove colony 21D10 because it has ITS2 profile C116
H10_depth <- sym[sym$Species.haplotype== "Haplotype 10",]

#Rename depths
H10_depth$Depth.m <- ifelse(grepl("5",H10_depth$Depth.m),"5m", as.character(H10_depth$Depth.m))
H10_depth$Depth.m <- ifelse(grepl("10",H10_depth$Depth.m),"10m", as.character(H10_depth$Depth.m))
H10_depth$Depth.m <- ifelse(grepl("20",H10_depth$Depth.m),"20m", as.character(H10_depth$Depth.m))

#### PCoA analysis
# Create distance matrix
sym_H10_DIV <- H10_depth[, -c(1:7)]
sym_H10 <- vegdist(sym_H10_DIV, method = "bray")

# Running PCoA in vegan (cmdscale)
mod_sym_H10 <- cmdscale(sym_H10,eig=T)

#Extract points
data.scores = as.data.frame(scores(mod_sym_H10))
xy <- data.scores[,1:2]
plot(xy)
mod_sym_H10

#PERMANOVA
adonis(sym_H10_DIV ~ H10_depth$Depth.m, method = "bray", permutations = 99999)

#Pairwise PERMANOVA for depths
pw <- pairwise.adonis2(H10_depth[,7:487]~Depth.m, data=H10_depth, permutations = 99999)
pw

## Plot 
# Percent of variance explained by axis
round(mod_sym_H10$eig*100/sum(mod_sym_H10$eig),1)

quartz(width=10,height=6)
par(mfrow=c(1,2),mar=c(5,4,4,0))
plot(xy,xlim=c(-0.21,0.2),ylim=c(-0.2,0.28),type="n", xlab="PCoA1 (26.7%)", ylab="PCoA2 (24.6%)")

# Haplotype 10 colors and symbols
my_colors_to_plot <- my_colors_profile[rownames(my_colors_profile) %in% unique(H10_depth$Type_profile),]
my_profiles_to_plot <- rownames(my_colors_profile)[rownames(my_colors_profile) %in% unique(H10_depth$Type_profile)]
my_coralID_symbols <- cbind.data.frame(Depth=unique(H10_depth$Depth.m), symbol=c(15,16,17))
# Reorder depths
my_coralID_symbols$Depth = factor(my_coralID_symbols$Depth,levels=c("5m","10m","20m"))

# Add points, colored by ITS2 Profile, with symbols based on depth.
for(i in 1:length(my_profiles_to_plot)){
  
  # Create a vector of pch symbols that match that desired for each depth
  tmp <- H10_depth[which(H10_depth$Type_profile==my_profiles_to_plot[i]),]
  pch.vec <- my_coralID_symbols[match(tmp$Depth,my_coralID_symbols$Depth),2]
  
  points(x=xy[which(H10_depth$Type_profile==my_profiles_to_plot[i]),1],
         y=xy[which(H10_depth$Type_profile==my_profiles_to_plot[i]),2],
         col = my_colors_to_plot[i],
         pch=pch.vec,
         cex=1)
}
ordiellipse(mod_sym_H10, H10_depth$Depth.m, lty = 'dashed', label=F, cex = 0.95)

par(mar=c(5,1,4,1),xpd=T)
plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",axes=F,xlab="",ylab="")
legend(0,1.05, legend=my_profiles_to_plot, pch=19, col=my_colors_to_plot, cex=0.85, bty="n")
legend(0,0, legend= my_coralID_symbols$Depth, pch = my_coralID_symbols$symbol, cex=1.2, bty="n")

------------------- ##### P. meandrina and Haplotype 8a by depth ##### -------------------

# Remove Species/haplotypes and samples not of interest and sample C68, which hosts C42.2.C116.C1.C1b 
Pmean_depth <- sym[sym$Species.haplotype!= "Haplotype 10" 
                        & sym$Species.haplotype!= "Haplotype 11"
                        & sym$Species.haplotype!= "Haplotype 2"
                        & sym$Species.haplotype!= "P. verrucosa"
                        & sym$Species.haplotype!= "P. grandis"
                        & sym$Coral.ID!= "C68",]

#Rename depths
Pmean_depth$Depth.m <- ifelse(grepl("5",Pmean_depth$Depth.m),"5m", as.character(Pmean_depth$Depth.m))
Pmean_depth$Depth.m <- ifelse(grepl("10",Pmean_depth$Depth.m),"10m", as.character(Pmean_depth$Depth.m))
Pmean_depth$Depth.m <- ifelse(grepl("20",Pmean_depth$Depth.m),"20m", as.character(Pmean_depth$Depth.m))

#### PCoA analysis
# Create distance matrix
sym_Pmean_DIV <- Pmean_depth[, -c(1:6)]
sym_Pmean <- vegdist(sym_Pmean_DIV, method = "bray")

# Running PCoA in vegan (cmdscale)
mod_sym_Pmean <- cmdscale(sym_Pmean,eig=T)

#Extract points
data.scores = as.data.frame(scores(mod_sym_Pmean))
xy <- data.scores[,1:2]
plot(xy)
mod_sym_Pmean

#PERMANOVA
adonis(sym_Pmean_DIV ~ Pmean_depth$Depth.m, method = "bray", permutations = 99999)

#Pairwise PERMANOVA for depths
pw <- pairwise.adonis2(Pmean_depth[,7:487]~Depth.m, data=Pmean_depth, permutations = 99999)
pw

## Plot
# Percent of variance explained by axis
round(mod_sym_Pmean$eig*100/sum(mod_sym_Pmean$eig),1)

quartz(width=10,height=6)
par(mfrow=c(1,2),mar=c(5,4,4,0))
plot(xy,xlim=c(-0.32,0.13),ylim=c(-0.15,0.45),type="n", xlab="PCoA1 (30.7%)", ylab="PCoA2 (19.9%)")

# P. meandrina and Haplotype 8a colors and symbols
my_colors_to_plot <- my_colors_profile[rownames(my_colors_profile) %in% unique(Pmean_depth$Type_profile),]
my_profiles_to_plot <- rownames(my_colors_profile)[rownames(my_colors_profile) %in% unique(Pmean_depth$Type_profile)]
my_coralID_symbols <- cbind.data.frame(Depth=unique(Pmean_depth$Depth.m), symbol=c(15,16,17))
# Reorder depths
my_coralID_symbols$Depth = factor(my_coralID_symbols$Depth,levels=c("5m","10m","20m"))

# Add points, colored by ITS2 Profile, with symbols based on depth.
for(i in 1:length(my_profiles_to_plot)){
  
  # Create a vector of pch symbols that match that desired for each depth
  tmp <- Pmean_depth[which(Pmean_depth$Type_profile==my_profiles_to_plot[i]),]
  pch.vec <- my_coralID_symbols[match(tmp$Depth,my_coralID_symbols$Depth),2]
  
  points(x=xy[which(Pmean_depth$Type_profile==my_profiles_to_plot[i]),1],
         y=xy[which(Pmean_depth$Type_profile==my_profiles_to_plot[i]),2],
         col = my_colors_to_plot[i],
         pch=pch.vec,
         cex=1)
}
ordiellipse(mod_sym_Pmean, Pmean_depth$Depth.m, lty = 'dashed', label=F, cex = 0.95)
# Plot legend
par(mar=c(5,1,4,1),xpd=T)
plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",axes=F,xlab="",ylab="")
legend(0,1.05, legend=my_profiles_to_plot, pch=19, col=my_colors_to_plot, cex=0.85, bty="n")
legend(0,0, legend= my_coralID_symbols$Depth, pch = my_coralID_symbols$symbol, cex=1.2, bty="n")

------------------- ##### P. grandis by depth ##### -------------------
## Prepare data
# Remove Species/haplotypes and samples not of interest and colonies that have C116, D, or A type profiles
Peyd_depth <- sym[sym$Species.haplotype== "P. grandis"
                       & sym$Coral.ID!="22C17"
                       & sym$Coral.ID!="G05"
                       & sym$Coral.ID!="H15",]

#Rename depths
Peyd_depth$Depth.m <- ifelse(grepl("5",Peyd_depth$Depth.m),"5m", as.character(Peyd_depth$Depth.m))
Peyd_depth$Depth.m <- ifelse(grepl("10",Peyd_depth$Depth.m),"10m", as.character(Peyd_depth$Depth.m))
Peyd_depth$Depth.m <- ifelse(grepl("20",Peyd_depth$Depth.m),"20m", as.character(Peyd_depth$Depth.m))

#### PCoA analysis
# Create distance matrix
sym_Peyd_DIV <- Peyd_depth[, -c(1:6)]
sym_Peyd <- vegdist(sym_Peyd_DIV, method = "bray")

# Running PCoA in vegan (cmdscale)
mod_sym_Peyd <- cmdscale(sym_Peyd,eig=T)

#Extract points
data.scores = as.data.frame(scores(mod_sym_Peyd))
xy <- data.scores[,1:2]
plot(xy)
mod_sym_Peyd

#PERMANOVA
adonis(sym_Peyd_DIV ~ Peyd_depth$Depth.m, method = "bray", permutations = 99999)

#Pairwise PERMANOVA for depths
pw <- pairwise.adonis2(Peyd_depth[,7:487]~Depth.m, data=Peyd_depth, permutations = 99999)
pw


## Plot
# Percent of variance explained by axis
round(mod_sym_Peyd$eig*100/sum(mod_sym_Peyd$eig),1)

quartz(width=10,height=6)
par(mfrow=c(1,2),mar=c(5,4,4,0))
plot(xy,xlim=c(-0.25,0.35),ylim=c(-0.25,0.2),type="n", xlab="PCoA1 (60%)", ylab="PCoA2 (10.5%)")

# P. grandis colors and symbols
my_colors_to_plot <- my_colors_profile[rownames(my_colors_profile) %in% unique(Peyd_depth$Type_profile),]
my_profiles_to_plot <- rownames(my_colors_profile)[rownames(my_colors_profile) %in% unique(Peyd_depth$Type_profile)]
my_coralID_symbols <- cbind.data.frame(Depth=unique(Peyd_depth$Depth.m), symbol=c(15,16,17))
# Reorder depths
my_coralID_symbols$Depth = factor(my_coralID_symbols$Depth,levels=c("5m","10m","20m"))

# Add points, colored by ITS2 Profile, with symbols based on depth.
for(i in 1:length(my_profiles_to_plot)){
  
  # Create a vector of pch symbols that match that desired for each depth
  tmp <- Peyd_depth[which(Peyd_depth$Type_profile==my_profiles_to_plot[i]),]
  pch.vec <- my_coralID_symbols[match(tmp$Depth,my_coralID_symbols$Depth),2]
  
  points(x=xy[which(Peyd_depth$Type_profile==my_profiles_to_plot[i]),1],
         y=xy[which(Peyd_depth$Type_profile==my_profiles_to_plot[i]),2],
         col = my_colors_to_plot[i],
         pch=pch.vec,
         cex=1)
}
ordiellipse(mod_sym_Peyd, Peyd_depth$Depth.m, lty = 'dashed', label=F, cex = 0.95)
# Plot legend
par(mar=c(5,1,4,1),xpd=T)
plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",axes=F,xlab="",ylab="")
legend(0,1.05, legend=my_profiles_to_plot, pch=19, col=my_colors_to_plot, cex=0.85, bty="n")
legend(0,0, legend= my_coralID_symbols$Depth, pch = my_coralID_symbols$symbol, cex=1.2, bty="n")

------------------- ##### Haplotype 11 and 2 by depth ##### -------------------
## Prepare data
# Remove Species/haplotypes and samples not of interest and colonies with ITS2 type D profile
H11_depth <- sym[sym$Species.haplotype!= "Haplotype 10" 
                     & sym$Species.haplotype!= "P. verrucosa"
                     & sym$Species.haplotype!= "P. meandrina"
                     & sym$Species.haplotype!= "P. grandis"
                     & sym$Species.haplotype!= "Haplotype 8a"
                     & sym$Coral.ID!="61A58"
                     & sym$Coral.ID!="G44"
                     & sym$Coral.ID!="G49"
                     & sym$Coral.ID!="F04"
                     & sym$Coral.ID!="J17",]  

#Rename depths
H11_depth$Depth.m <- ifelse(grepl("5",H11_depth$Depth.m),"5m", as.character(H11_depth$Depth.m))
H11_depth$Depth.m <- ifelse(grepl("10",H11_depth$Depth.m),"10m", as.character(H11_depth$Depth.m))
H11_depth$Depth.m <- ifelse(grepl("20",H11_depth$Depth.m),"20m", as.character(H11_depth$Depth.m))

#### PCoA analysis
# Create distance matrix
sym_H11_DIV <- H11_depth[, -c(1:6)]
sym_H11 <- vegdist(sym_H11_DIV, method = "bray")

# Running PCoA in vegan (cmdscale)
mod_sym_H11 <- cmdscale(sym_H11,eig=T)

#Extract points
data.scores = as.data.frame(scores(mod_sym_H11))
xy <- data.scores[,1:2]
plot(xy)
mod_sym_H11

#PERMANOVA
adonis(sym_H11_DIV ~ H11_depth$Depth.m, method = "bray", permutations = 99999)

#Pairwise PERMANOVA for depths
pw <- pairwise.adonis2(H11_depth[,7:487]~Depth.m, data=H11_depth, permutations = 99999)
pw


## Plot
# Percent of variance explained by axis
round(mod_sym_H11$eig*100/sum(mod_sym_H11$eig),1)

quartz(width=10,height=6)
par(mfrow=c(1,2),mar=c(5,4,4,0))
plot(xy,xlim=c(-0.25,0.25),ylim=c(-0.25,0.07),type="n", xlab="PCoA1 (39.6%)", ylab="PCoA2 (20.3%)")

# Haplotype 11 and haplotype 2 colors and symbols
my_colors_to_plot <- my_colors_profile[rownames(my_colors_profile) %in% unique(H11_depth$Type_profile),]
my_profiles_to_plot <- rownames(my_colors_profile)[rownames(my_colors_profile) %in% unique(H11_depth$Type_profile)]
my_coralID_symbols <- cbind.data.frame(Depth=unique(H11_depth$Depth.m), symbol=c(15,16,17))
# Reorder depths
my_coralID_symbols$Depth = factor(my_coralID_symbols$Depth,levels=c("5m","10m","20m"))

# Add points, colored by ITS2 Profile, with symbols based on depth.
for(i in 1:length(my_profiles_to_plot)){
  
  # Create a vector of pch symbols that match that desired for each depth
  tmp <- H11_depth[which(H11_depth$Type_profile==my_profiles_to_plot[i]),]
  pch.vec <- my_coralID_symbols[match(tmp$Depth,my_coralID_symbols$Depth),2]
  
  points(x=xy[which(H11_depth$Type_profile==my_profiles_to_plot[i]),1],
         y=xy[which(H11_depth$Type_profile==my_profiles_to_plot[i]),2],
         col = my_colors_to_plot[i],
         pch=pch.vec,
         cex=1)
}
ordiellipse(mod_sym_H11, H11_depth$Depth.m, lty = 'dashed', label=F, cex = 0.95)
# Plot legend
par(mar=c(5,1,4,1),xpd=T)
plot(c(0,1),c(0,1),type="n",xaxt="n",yaxt="n",axes=F,xlab="",ylab="")
legend(0,1.05, legend=my_profiles_to_plot, pch=19, col=my_colors_to_plot, cex=0.85, bty="n")
legend(0,0, legend= my_coralID_symbols$Depth, pch = my_coralID_symbols$symbol, cex=1.2, bty="n")

