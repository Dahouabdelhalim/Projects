library(MASS)
library(vegan)
library(ggplot2)

#############################################  mangrove  NMDS 
mangrove <- read.delim('mangrove_table.txt', row.names = 1, sep = '\\t', stringsAsFactors = FALSE, check.names = FALSE)
head(mangrove)
mangrove <- data.frame(t(mangrove))
bray_dis <- vegdist(mangrove, method = 'bray')
write.table(as.matrix(bray_dis), 'bray_distance_mangrove.txt', sep = '\\t', col.names = NA, quote = FALSE)

nmds_dis <- metaMDS(bray_dis, k = 2)

#ordiplot() 
nmds_dis_species <- wascores(nmds_dis$points, mangrove)
species<- data.frame(nmds_dis_species)
species$name <- rownames(species)
nmds_dis_site <- data.frame(nmds_dis$points)
par(mfrow = c(2, 2))
ordiplot(nmds_dis, type = 'none', main = paste('Stress =', round(nmds_dis$stress, 4)))
points(nmds_dis_species, pch = 3, cex = 0.5, col = 'gray')
points(nmds_dis, pch = 19, cex = 0.7, col = c(rep('red', 5), rep('green3', 5), rep('blue', 5)))

################################ Permanova test for mangrove community

mangrove <- read.delim('mangrove_table.txt', row.names = 1, sep = '\\t', stringsAsFactors = FALSE, check.names = FALSE)
mangrove <- data.frame(t(mangrove))
group <- read.delim('mangrove_group.txt', sep = '\\t', stringsAsFactors = FALSE)
bray_dis <- vegdist(mangrove, method = 'bray')
adonis_result <- adonis(bray_dis~site, group, permutations = 999)
adonis_result
summary(adonis_result)



################################################ macrobenthos NMDS 
benthos <- read.delim('benthos_table.txt', row.names = 1, sep = '\\t', stringsAsFactors = FALSE, check.names = FALSE)
head(benthos)
benthos <- data.frame(t(benthos))
bray_dis <- vegdist(benthos, method = 'bray')
write.table(as.matrix(bray_dis), 'bray_distance_benthos.txt', sep = '\\t', col.names = NA, quote = FALSE)
#NMDS 
nmds_dis <- metaMDS(bray_dis, k = 2)

#ordiplot() 
nmds_dis_species <- wascores(nmds_dis$points, benthos)
species<- data.frame(nmds_dis_species)
species$name <- rownames(species)
nmds_dis_site <- data.frame(nmds_dis$points)
par(mfrow = c(2, 2))
ordiplot(nmds_dis, type = 'none', main = paste('Stress =', round(nmds_dis$stress, 4)))
points(nmds_dis_species, pch = 3, cex = 0.5, col = 'gray')
points(nmds_dis, pch = 19, cex = 0.7, col = c(rep('red', 4), rep('green3', 4), rep('blue', 4)))

################################ Permanova test for macrobenthos community

benthos <- read.delim('benthos_table.txt', row.names = 1, sep = '\\t', stringsAsFactors = FALSE, check.names = FALSE)
benthos <- data.frame(t(benthos))
group <- read.delim('benthos_group.txt', sep = '\\t', stringsAsFactors = FALSE)
bray_dis <- vegdist(benthos, method = 'bray')
adonis_result <- adonis(bray_dis~site, group, permutations = 999)
adonis_result


################################################ fish NMDS 
fish <- read.delim('fish_table.txt', row.names = 1, sep = '\\t', stringsAsFactors = FALSE, check.names = FALSE)
head(fish)
fish <- data.frame(t(fish))
bray_dis <- vegdist(fish, method = 'bray')
write.table(as.matrix(bray_dis), 'bray_distance_fish.txt', sep = '\\t', col.names = NA, quote = FALSE)
#NMDS 
nmds_dis <- metaMDS(bray_dis, k = 2)

#ordiplot() 
nmds_dis_species <- wascores(nmds_dis$points, fish)
species<- data.frame(nmds_dis_species)
species$name <- rownames(species)
nmds_dis_site <- data.frame(nmds_dis$points)
par(mfrow = c(2, 2))
ordiplot(nmds_dis, type = 'none', main = paste('Stress =', round(nmds_dis$stress, 4)))
points(nmds_dis_species, pch = 3, cex = 0.5, col = 'gray')
points(nmds_dis, pch = 19, cex = 0.7, col = c(rep('red', 4), rep('green3', 4), rep('blue', 4)))

################################ Permanova test for fish community

fish <- read.delim('fish_table_v2.txt', row.names = 1, sep = '\\t', stringsAsFactors = FALSE, check.names = FALSE)
fish <- data.frame(t(fish))
group <- read.delim('fish_group.txt', sep = '\\t', stringsAsFactors = FALSE)
bray_dis <- vegdist(fish, method = 'bray')
adonis_result <- adonis(bray_dis~site, group, permutations = 999)
adonis_result

###################################################  waterbird  NMDS 
bird <- read.delim('bird_table.txt', row.names = 1, sep = '\\t', stringsAsFactors = FALSE, check.names = FALSE)
head(bird)
bird <- data.frame(t(bird))
bray_dis <- vegdist(bird, method = 'bray')
write.table(as.matrix(bray_dis), 'bray_distance_bird.txt', sep = '\\t', col.names = NA, quote = FALSE)
#NMDS 
nmds_dis <- metaMDS(bray_dis, k = 2)

#ordiplot() 
nmds_dis_species <- wascores(nmds_dis$points, bird)
species<- data.frame(nmds_dis_species)
species$name <- rownames(species)
nmds_dis_site <- data.frame(nmds_dis$points)
par(mfrow = c(2, 2))
ordiplot(nmds_dis, type = 'none', main = paste('Stress =', round(nmds_dis$stress, 4)))
points(nmds_dis_species, pch = 3, cex = 0.5, col = 'gray')
points(nmds_dis, pch = 19, cex = 0.7, col = c(rep('red', 5), rep('green3', 5), rep('blue', 5)))

################################ Permanova test for waterbird community

bird <- read.delim('bird_table.txt', row.names = 1, sep = '\\t', stringsAsFactors = FALSE, check.names = FALSE)
bird <- data.frame(t(bird))
group <- read.delim('bird_group.txt', sep = '\\t', stringsAsFactors = FALSE)
bray_dis <- vegdist(bird, method = 'bray')
adonis_result <- adonis(bray_dis~site, group, permutations = 999)
adonis_result
summary(adonis_result)