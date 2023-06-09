library(geomorph)
library(abind)

### load landmarks and links connecting them
landmarks_rep1 <- readland.tps(paste(getwd(), "/Phrynus_group_B_species_replicate_1_curves_appended.tps", sep=""))
landmarks_rep2 <- readland.tps(paste(getwd(), "/Phrynus_group_B_species_replicate_2_curves_appended.tps", sep=""))
landmarks_links <- read.csv(paste(getwd(), "/links.csv", sep=""))
landmarks_links <- as.matrix(landmarks_links)
specimen_info <- read.csv(paste(getwd(), "/specimen_descriptions.csv", sep=""))

### plot landmarks before gpa
plotAllSpecimens(landmarks_rep1, mean=FALSE)                 
plotAllSpecimens(landmarks_rep2, mean=FALSE) 

### perform gpa              
landmarks_gpa_rep1 <- gpagen(landmarks_rep1, print.progress = FALSE)
summary(landmarks_gpa_rep1)
landmarks_gpa_rep2 <- gpagen(landmarks_rep2, print.progress = FALSE)
summary(landmarks_gpa_rep2)

### generate means of replicates

coordinates_rep1 <- landmarks_gpa_rep1$coords
coordinates_rep2 <- landmarks_gpa_rep2$coords

specimen_coordinates_rep1 <- list()
specimen_coordinates_rep2 <- list()
bound_specimen_replicates <- list()
mean_shape_specimens_coordinates <- list()

for (i in 1:length(coordinates_rep1[1,1,])) {

specimen_coordinates_rep1[[i]] <- landmarks_gpa_rep1$coords[,,i]
specimen_coordinates_rep2[[i]] <- landmarks_gpa_rep2$coords[,,i]
bound_specimen_replicates[[i]] <- abind(specimen_coordinates_rep1[[i]],specimen_coordinates_rep2[[i]], along=3)
mean_shape_specimens_coordinates[[i]] <- mshape(bound_specimen_replicates[[i]])

}

array_mean_coordinates <- array(as.numeric(unlist(mean_shape_specimens_coordinates)), dim = c(length(coordinates_rep1[,1,1]), length(coordinates_rep1[1,,1]), i))

### plot points of all specimens, the mean and the links connecting the means
plotAllSpecimens(landmarks_gpa_rep1$coords,links=landmarks_links, mean=TRUE, plot.param = list(pt.bg = "green", mean.cex=1, link.col="red",txt.pos=3, txt.cex=1))
plotAllSpecimens(landmarks_gpa_rep2$coords,links=landmarks_links, mean=TRUE, plot.param = list(pt.bg = "green", mean.cex=1, link.col="red",txt.pos=3, txt.cex=1))

###  Traditional PCA 
PCA_rep1 <- gm.prcomp(landmarks_gpa_rep1$coords)
summary(PCA_rep1)
PCA_rep2 <- gm.prcomp(landmarks_gpa_rep2$coords)
summary(PCA_rep2)
PCA_means <- gm.prcomp(array_mean_coordinates)
summary(PCA_means)

# plot PCA
pdf(file = paste(getwd(), '/results/PCA_rep1.pdf', sep = ''))
plot(PCA_rep1, main = "PCA_rep1", pch = specimen_info$Sex, col = specimen_info$color)
dev.off()
pdf(file = paste(getwd(), '/results/PCA_rep2.pdf', sep = ''))
plot(PCA_rep2, main = "PCA_rep2", pch = specimen_info$Sex, col = specimen_info$color)
dev.off()
pdf(file = paste(getwd(), '/results/PCA_means.pdf', sep = ''))
plot(PCA_means, main = "PCA_means", pch = specimen_info$Sex, col = specimen_info$color)
dev.off()

### plot species separately
P_calypso <- specimen_info$Species == "Phrynus_calypso"
pdf(file = paste(getwd(), '/results/P_calypso.pdf', sep = ''))
plotAllSpecimens(array_mean_coordinates[,,P_calypso],links=landmarks_links, mean=TRUE, plot.param = list(pt.bg = "green", mean.cex=2, link.col="red",txt.pos=3, txt.cex=1))
dev.off()

P_exsul <- specimen_info$Species == "Phrynus_exsul"
pdf(file = paste(getwd(), '/results/P_exsul.pdf', sep = ''))
plotAllSpecimens(array_mean_coordinates[,,P_exsul],links=landmarks_links, mean=TRUE, plot.param = list(pt.bg = "green", mean.cex=2, link.col="red",txt.pos=3, txt.cex=1))
dev.off()

P_goesii <- specimen_info$Species == "Phrynus_goesii"
pdf(file = paste(getwd(), '/results/P_goesii.pdf', sep = ''))
plotAllSpecimens(array_mean_coordinates[,,P_goesii],links=landmarks_links, mean=TRUE, plot.param = list(pt.bg = "green", mean.cex=2, link.col="red",txt.pos=3, txt.cex=1))
dev.off()

P_sp.nov <- specimen_info$Species == "Phrynus_sp.nov"
pdf(file = paste(getwd(), '/results/P_sp.nov.pdf', sep = ''))
plotAllSpecimens(array_mean_coordinates[,,P_sp.nov],links=landmarks_links, mean=TRUE, plot.param = list(pt.bg = "green", mean.cex=2, link.col="red",txt.pos=3, txt.cex=1))
dev.off()

P_longipes <- specimen_info$Species == "Phrynus_longipes"
pdf(file = paste(getwd(), '/results/P_longipes.pdf', sep = ''))
plotAllSpecimens(array_mean_coordinates[,,P_longipes],links=landmarks_links, mean=TRUE, plot.param = list(pt.bg = "green", mean.cex=2, link.col="red",txt.pos=3, txt.cex=1))
dev.off()

P_pinarensis <- specimen_info$Species == "Phrynus_pinarensis"
pdf(file = paste(getwd(), '/results/P_pinarensis.pdf', sep = ''))
plotAllSpecimens(array_mean_coordinates[,,P_pinarensis],links=landmarks_links, mean=TRUE, plot.param = list(pt.bg = "green", mean.cex=2, link.col="red",txt.pos=3, txt.cex=1))
dev.off()

P_tessellatus <- specimen_info$Species == "Phrynus_tessellatus"
pdf(file = paste(getwd(), '/results/P_tessellatus.pdf', sep = ''))
plotAllSpecimens(array_mean_coordinates[,,P_tessellatus],links=landmarks_links, mean=TRUE, plot.param = list(pt.bg = "green", mean.cex=2, link.col="red",txt.pos=3, txt.cex=1))
dev.off()