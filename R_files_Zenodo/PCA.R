#author: "Giorgia Auteri"
#Script modified from Auteri, Giorgia G. and L. Lacey Knowles. 2020. Decimated little brown bats
# show potential for adaptive change. Scientific reports: 10, 3023. doi https://doi.org/10.1038/s41598-020-59797-4

library(adegenet)
library(plyr)
library(vegan)

#Clear working directory
rm(list = ls())

# Import data
gind1 <- read.structure("Structure25Miss.stru", 
                        n.ind = 15,    # Number of individuals
                        n.loc = 1880, # Number of variant sites (SNPs)
                        col.lab = 1, 
                        row.marknames = 0, 
                        col.pop = 0, 
                        onerowperind = TRUE, 
                        NA.char=0, 
                        ask = FALSE)
gind <- scaleGen (gind1, NA.method = c("mean")) #default center & scale = TRUE, 

#PCA
pca<-dudi.pca(df = gind, center = T, scale = F, scannf = FALSE, nf = 2) #selects the first 2 axes (nf=2) 

gind[,0] #list all names and their position
pca$li #list the pca axis

pca$eig[1]/sum(pca$eig) #calculates % explanation of each pc axis
pca$eig[2]/sum(pca$eig)
pca$eig[3]/sum(pca$eig)

#####PCA 1 and 2
IndInfo <- read.table("IndInfo.txt",
                      header = TRUE,
                      sep = "\\t")

plotDat <- pca$li[,1:2]
plotDat$ID <- NA
plotDat$ID <- rownames(plotDat)
plotDat <- merge(x = plotDat,
                       y = IndInfo,
                       by.x = "ID",
                       by.y = "ID")

library(ggplot2)
library(RColorBrewer)
ggplot(plotDat, aes(x=Axis1, y=Axis2, group=(Pop))) +
  geom_point(aes(shape=IR, fill=as.character(Pop)), size = 3) + 
  scale_color_brewer(palette="Dark2")+
  scale_shape_manual(values = c(21, 23)) +
  guides(fill = guide_legend(override.aes=list(shape=21)))+ #legend corresponding colors
  theme_bw()
