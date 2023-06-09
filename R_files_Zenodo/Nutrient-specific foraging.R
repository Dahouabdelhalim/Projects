##### Libraries #####

library('mvabund')
library('ggtern')
library('vegan')
library("flashClust")
library("dendextend")
library("plyr")
library("dplyr")
library("ggplot2")
library("viridis")
library("RColorBrewer")
library("cluster")
library("cooccur")
library("ggrepel")
library("clValid")
library("econullnetr")
library("gplots")
library("ggalt")
library("car")
library("scales")
library("Rcpp")

##### Tropho-species cluster method determination #####

tropho <- read.csv("Mean macros per taxon.csv")
rownames(tropho) <- tropho[,1]
tropho_taxon <- tropho$Taxon
tropho <- tropho[2:4]
summary(tropho)


tropho_sc <- as.data.frame(scale(tropho))
summary(tropho_sc)

trophodist<- dist(tropho_sc, method = "euclidean")

# average

trophotreeAVG <- hclust(trophodist, method = "average")
plot(trophotreeAVG, main="")

x <- c(3:43)
for (i in x) {
  trophocut_avg <- cutree(trophotreeAVG, k = i )
  trophodunn <- dunn(distance= trophodist, clusters = trophocut_avg, method= 'euclidean') 
  print(trophodunn)
}

plot(trophotreeAVG, main="")
rect.hclust(trophotreeAVG, k = 41, border = 2:28)

trophocut_avg41 <- cutree(trophotreeAVG, k = 41)
trophodunn_avg41 <- dunn(distance= trophodist, clusters = trophocut_avg41, method= 'euclidean') 
trophodunn_avg41


# single

trophotreeSIN <- hclust(trophodist, method = "single")
plot(trophotreeSIN, main="")

x <- c(5:50)
for (i in x) {
  trophocut_sin <- cutree(trophotreeSIN, k = i )
  trophodunn <- dunn(distance= trophodist, clusters = trophocut_sin, method= 'euclidean') 
  print(trophodunn)
}

rect.hclust(trophotreeSIN, k = 5, border = 2:28)


trophocut_sin5 <- cutree(trophotreeSIN, k = 5)
trophodunn_sin5 <- dunn(distance= trophodist, clusters = trophocut_sin5, method= 'euclidean') 
trophodunn_sin5


# complete

trophotreeCOM <- hclust(trophodist, method = "complete")
plot(trophotreeCOM, main="")

x <- c(5:21)
for (i in x) {
  trophocut_com <- cutree(trophotreeCOM, k = i )
  trophodunn <- dunn(distance= trophodist, clusters = trophocut_com, method= 'euclidean') 
  print(trophodunn)
}

rect.hclust(trophotreeCOM, k = 20, border = 2:28)

trophocut_com20 <- cutree(trophotreeCOM, k = 20)
trophodunn_com20 <- dunn(distance= trophodist, clusters = trophocut_com20, method= 'euclidean') 
trophodunn_com20


# mcquitty

trophotreeMCQ <- hclust(trophodist, "mcquitty")
plot(trophotreeMCQ, main="")

x <- c(5:50)
for (i in x) {
  trophocut_mcq <- cutree(trophotreeMCQ, k = i )
  trophodunn <- dunn(distance= trophodist, clusters = trophocut_mcq, method= 'euclidean') 
  print(trophodunn)
}

rect.hclust(trophotreeMCQ, k = 29, border = 2:28)

trophocut_mcq29 <- cutree(trophotreeMCQ, k = 29)
trophodunn_mcq29 <- dunn(distance= trophodist, clusters = trophocut_mcq29, method= 'euclidean') 
trophodunn_mcq29


# median

trophotreeMED <- hclust(trophodist, "median")
plot(trophotreeMED, main="")

x <- c(5:40)
for (i in x) {
  trophocut_med <- cutree(trophotreeMED, k = i )
  trophodunn <- dunn(distance= trophodist, clusters = trophocut_med, method= 'euclidean') 
  print(trophodunn)
}

rect.hclust(trophotreeMED, k = 31, border = 2:28)

trophocut_med31 <- cutree(trophotreeMED, k = 31)
trophodunn_med31 <- dunn(distance= trophodist, clusters = trophocut_med31, method= 'euclidean') 
trophodunn_med31


# centroid

trophotreeCEN <- hclust(trophodist, "centroid")
plot(trophotreeCEN, main="")

x <- c(5:50)
for (i in x) {
  trophocut_cen <- cutree(trophotreeCEN, k = i )
  trophodunn <- dunn(distance= trophodist, clusters = trophocut_cen, method= 'euclidean') 
  print(trophodunn)
}

rect.hclust(trophotreeCEN, k = 43, border = 2:28)

trophocut_cen43 <- cutree(trophotreeAVG, k = 43)
trophodunn_cen43 <- dunn(distance= trophodist, clusters = trophocut_cen43, method= 'euclidean') 
trophodunn_cen43


##### Tropho-species clustering #####

trophotreeCOM <- hclust(trophodist, method = "complete")
plot(trophotreeCOM, main="")

x <- c(5:21)
for (i in x) {
  trophocut_com <- cutree(trophotreeCOM, k = i )
  trophodunn <- dunn(distance= trophodist, clusters = trophocut_com, method= 'euclidean') 
  print(trophodunn)
}

plot(trophotreeCOM, main="")
rect.hclust(trophotreeCOM, k = 20, border = 2:28)

trophocut_com20 <- cutree(trophotreeCOM, k = 20)
trophodunn_com20 <- dunn(distance= trophodist, clusters = trophocut_com20, method= 'euclidean') 
trophodunn_com20

tropho_cl <- mutate(tropho, cluster = trophocut_com20)
count(tropho_cl,cluster)

ggplot(tropho_cl, aes(x=Lipid, y = Protein, color = factor(cluster))) + geom_point()

TrophoClusters <- table(tropho_cl$cluster,tropho_taxon)

write.csv(TrophoClusters, "TrophoClusters.csv")

plot(trophotreeCOM, main="")
rect.hclust(trophotreeCOM, k = 20, border = 2:28)
abline(h = 0.74, col = 'red')

pdf("TS dendrogram.pdf", width = 12, height = 8) 
plot(trophotreeCOM, main="")
rect.hclust(trophotreeCOM, k = 20, border = 2:28)
dev.off()

tiff(file = "TS dendrogram.tiff", width = 24, height = 16, res = 600, units = "cm")
plot(trophotreeCOM, main="")
rect.hclust(trophotreeCOM, k = 20, border = 2:28)
dev.off()

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x,method="euclidean")
d <- distfunc(tropho_sc)
fit <- hclustfunc(d)
clusters <- cutree(fit, h=0.74)
nofclust.height <-  length(unique(as.vector(clusters)));

cl.row <- hclustfunc(distfunc(tropho_sc))
cl.col <- hclustfunc(distfunc(t(tropho_sc)))

hmcols <- rev(viridis(2750))
selcol <- colorRampPalette(brewer.pal(12,"Set3"))
selcol2 <- colorRampPalette(brewer.pal(8,"Accent"))
clustcol.height = selcol2(nofclust.height);

heatmap.2(as.matrix(tropho_sc), 
          trace='none', 
          dendrogram='both', 
          key=F,
          Colv=T, 
          scale='row',
          hclust=hclustfunc, distfun=distfunc, col=hmcols,
          symbreak=T,
          margins=c(7,10), keysize=0.1,
          lwid=c(5,0.5,3), lhei=c(0.05,0.5),
          lmat=rbind(c(5,0,4),c(3,1,2)),
          labRow=rownames(tropho_sc),
          RowSideColors=clustcol.height[clusters], cexRow = 0.8, cexCol = 1.5)

pdf("TS dendroheatmap.pdf", width = 8, height = 12) 
heatmap.2(as.matrix(tropho_sc), 
          trace='none', 
          dendrogram='both', 
          key=F,
          Colv=T, 
          scale='row',
          hclust=hclustfunc, distfun=distfunc, col=hmcols,
          symbreak=T,
          margins=c(7,10), keysize=0.1,
          lwid=c(5,0.5,3), lhei=c(0.05,0.5),
          lmat=rbind(c(5,0,4),c(3,1,2)),
          labRow=rownames(tropho_sc),
          RowSideColors=clustcol.height[clusters], cexRow = 0.8, cexCol = 1.5)
dev.off()

tiff(file = "TS dendroheatmap.tiff", width = 16, height = 24, res = 600, units = "cm")
heatmap.2(as.matrix(tropho_sc), 
          trace='none', 
          dendrogram='both', 
          key=F,
          Colv=T, 
          scale='row',
          hclust=hclustfunc, distfun=distfunc, col=hmcols,
          symbreak=T,
          margins=c(7,10), keysize=0.1,
          lwid=c(5,0.5,3), lhei=c(0.05,0.5),
          lmat=rbind(c(5,0,4),c(3,1,2)),
          labRow=rownames(tropho_sc),
          RowSideColors=clustcol.height[clusters], cexRow = 0.8, cexCol = 1.5)
dev.off()


Clustpal <- brewer.pal(8, "Accent")
Clustpal <- colorRampPalette(Clustpal)(20)

ggtern(tropho_cl, aes(x=Lipid,y=Carbohydrate, z=Protein))+
  geom_text(aes(label = as.factor(cluster), colour = as.factor(cluster)),
            vjust=-0.40) +
  scale_colour_manual(values=Clustpal) +
  theme_bw() +
  guides(col=FALSE) +
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.7))   + 
  scale_L_continuous(limits=c(.0,.7))  + 
  scale_R_continuous(limits=c(.3,1))

pdf("TS cluster ggtern.pdf", width = 8, height = 8) 
ggtern(tropho_cl, aes(x=Lipid,y=Carbohydrate, z=Protein))+
  geom_text(aes(label = as.factor(cluster), colour = as.factor(cluster)),
            vjust=-0.40) +
  scale_colour_manual(values=Clustpal) +
  theme_bw() +
  guides(col=FALSE) +
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.7))   + 
  scale_L_continuous(limits=c(.0,.7))  + 
  scale_R_continuous(limits=c(.3,1))
dev.off()

tiff(file = "TS cluster ggtern.tiff", width = 16, height = 16, res = 600, units = "cm")
ggtern(tropho_cl, aes(x=Lipid,y=Carbohydrate, z=Protein))+
  geom_text(aes(label = as.factor(cluster), colour = as.factor(cluster)),
            vjust=-0.40) +
  scale_colour_manual(values=Clustpal) +
  theme_bw() +
  guides(col=FALSE) +
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.7))   + 
  scale_L_continuous(limits=c(.0,.7))  + 
  scale_R_continuous(limits=c(.3,1))
dev.off()

##### Tropho-species macro cluster method determination #####

TSn <- read.csv("Tropho macro cluster.csv")
summary(TSn)

rownames(TSn) <- TSn[,1]
TSn_cluster <- TSn$Cluster
TSncarb <- TSn[2]
summary(TSncarb)

TSncarb_sc <- as.data.frame(scale(TSncarb))
summary(TSncarb_sc)

carbdist<- dist(TSncarb_sc, method = "euclidean")

# average

carbtreeSIN <- hclust(carbdist, method = "average")
plot(carbtreeSIN, main="")

x <- c(5:20)
for (i in x) {
  carbcut_sin <- cutree(carbtreeSIN, k = i )
  carbdunn <- dunn(distance= carbdist, clusters = carbcut_sin, method= 'euclidean') 
  print(carbdunn)
}

plot(carbtreeSIN, main="")
rect.hclust(carbtreeSIN, k = 11, border = 2:28)

carbcut_sin11 <- cutree(carbtreeSIN, k = 11)
carbdunn_sin11 <- dunn(distance= carbdist, clusters = carbcut_sin11, method= 'euclidean') 
carbdunn_sin11

# single

carbtreeSIN <- hclust(carbdist, method = "single")
plot(carbtreeSIN, main="")

x <- c(5:20)
for (i in x) {
  carbcut_sin <- cutree(carbtreeSIN, k = i )
  carbdunn <- dunn(distance= carbdist, clusters = carbcut_sin, method= 'euclidean') 
  print(carbdunn)
}

rect.hclust(carbtreeSIN, k = 10, border = 2:28)

carbcut_sin10 <- cutree(carbtreeSIN, k = 10)
carbdunn_sin10 <- dunn(distance= carbdist, clusters = carbcut_sin10, method= 'euclidean') 
carbdunn_sin10


# complete

carbtreeCOM <- hclust(carbdist, method = "complete")
plot(carbtreeCOM, main="")

x <- c(5:20)
for (i in x) {
  carbcut_com <- cutree(carbtreeCOM, k = i )
  carbdunn <- dunn(distance= carbdist, clusters = carbcut_com, method= 'euclidean') 
  print(carbdunn)
}

rect.hclust(carbtreeCOM, k = 13, border = 2:28)

carbcut_com13 <- cutree(carbtreeCOM, k = 13)
carbdunn_com13 <- dunn(distance= carbdist, clusters = carbcut_com13, method= 'euclidean') 
carbdunn_com13


# mcquitty

carbtreeMCQ <- hclust(carbdist, "mcquitty")
plot(carbtreeMCQ, main="")

x <- c(5:20)
for (i in x) {
  carbcut_mcq <- cutree(carbtreeMCQ, k = i )
  carbdunn <- dunn(distance= carbdist, clusters = carbcut_mcq, method= 'euclidean') 
  print(carbdunn)
}

rect.hclust(carbtreeMCQ, k = 13, border = 2:28)

carbcut_mcq13 <- cutree(carbtreeSIN, k = 13)
carbdunn_mcq13 <- dunn(distance= carbdist, clusters = carbcut_mcq13, method= 'euclidean') 
carbdunn_mcq13


# median

carbtreeMED <- hclust(carbdist, "median")
plot(carbtreeMED, main="")

x <- c(5:20)
for (i in x) {
  carbcut_med <- cutree(carbtreeMED, k = i )
  carbdunn <- dunn(distance= carbdist, clusters = carbcut_med, method= 'euclidean') 
  print(carbdunn)
}

rect.hclust(carbtreeMED, k = 10, border = 2:28)

carbcut_med10 <- cutree(carbtreeMED, k = 10)
carbdunn_med10 <- dunn(distance= carbdist, clusters = carbcut_med10, method= 'euclidean') 
carbdunn_med10


# centroid

carbtreeCEN <- hclust(carbdist, "centroid")
plot(carbtreeCEN, main="")

x <- c(5:20)
for (i in x) {
  carbcut_cen <- cutree(carbtreeCEN, k = i )
  carbdunn <- dunn(distance= carbdist, clusters = carbcut_cen, method= 'euclidean') 
  print(carbdunn)
}

rect.hclust(carbtreeCEN, k = 10, border = 2:28)

carbcut_cen10 <- cutree(carbtreeCEN, k = 10)
carbdunn_cen10 <- dunn(distance= carbdist, clusters = carbcut_cen10, method= 'euclidean') 
carbdunn_cen10

##### Tropho-species macro clustering #####

# Carbohydrate clustering

TSncarb <- TSn[2]
summary(TSncarb)
TSncarb_sc <- as.data.frame(scale(TSncarb))
summary(TSncarb_sc)
carbdist<- dist(TSncarb_sc, method = "euclidean")

carbtreeSIN <- hclust(carbdist, method = "single")
plot(carbtreeSIN, main="")

x <- c(5:20)
for (i in x) {
  carbcut_sin <- cutree(carbtreeSIN, k = i )
  carbdunn <- dunn(distance= carbdist, clusters = carbcut_sin, method= 'euclidean') 
  print(carbdunn)
}

rect.hclust(carbtreeSIN, k = 10, border = 2:28)


carbcut_sin10 <- cutree(carbtreeSIN, k = 10)
carbdunn_sin10 <- dunn(distance= carbdist, clusters = carbcut_sin10, method= 'euclidean') 
carbdunn_sin10

carb_cl <- mutate(TSncarb, cluster = carbcut_sin10)
count(carb_cl,cluster)

CarbClusters <- table(carb_cl$cluster,TSn_cluster)

write.csv(CarbClusters, "CarbClusters.csv")

pdf("Carb dendrogram.pdf", width = 12, height = 8) 
plot(carbtreeSIN, main="")
rect.hclust(carbtreeSIN, k = 10, border = 2:28)
dev.off()

tiff(file = "Carb dendrogram.tiff", width = 24, height = 16, res = 600, units = "cm")
plot(carbtreeSIN, main="")
rect.hclust(carbtreeSIN, k = 10, border = 2:28)
dev.off()

# Lipid clustering

TSnlip <- TSn[3]
summary(TSnlip)
TSnlip_sc <- as.data.frame(scale(TSnlip))
summary(TSnlip_sc)
lipdist<- dist(TSnlip_sc, method = "euclidean")

liptreeSIN <- hclust(lipdist, method = "single")
plot(liptreeSIN, main="")

x <- c(5:20)
for (i in x) {
  lipcut_sin <- cutree(liptreeSIN, k = i )
  lipdunn <- dunn(distance= lipdist, clusters = lipcut_sin, method= 'euclidean') 
  print(lipdunn)
}

rect.hclust(liptreeSIN, k = 7, border = 2:28)


lipcut_sin7 <- cutree(liptreeSIN, k = 7)
lipdunn_sin7 <- dunn(distance= lipdist, clusters = lipcut_sin7, method= 'euclidean') 
lipdunn_sin7

lip_cl <- mutate(TSnlip, cluster = lipcut_sin7)
count(lip_cl,cluster)

LipClusters <- table(lip_cl$cluster,TSn_cluster)

write.csv(LipClusters, "LipClusters.csv")

pdf("Lip dendrogram.pdf", width = 12, height = 8) 
plot(liptreeSIN, main="")
rect.hclust(liptreeSIN, k = 7, border = 2:28)
dev.off()

tiff(file = "Lip dendrogram.tiff", width = 24, height = 16, res = 600, units = "cm")
plot(liptreeSIN, main="")
rect.hclust(liptreeSIN, k = 7, border = 2:28)
dev.off()

# Protein clustering

TSnprot <- TSn[4]
summary(TSnprot)
TSnprot_sc <- as.data.frame(scale(TSnprot))
summary(TSnprot_sc)
protdist<- dist(TSnprot_sc, method = "euclidean")

prottreeSIN <- hclust(protdist, method = "single")
plot(prottreeSIN, main="")

x <- c(5:20)
for (i in x) {
  protcut_sin <- cutree(prottreeSIN, k = i )
  protdunn <- dunn(distance= protdist, clusters = protcut_sin, method= 'euclidean') 
  print(protdunn)
}

rect.hclust(prottreeSIN, k = 6, border = 2:28)


protcut_sin6 <- cutree(prottreeSIN, k = 6)
protdunn_sin6 <- dunn(distance= protdist, clusters = protcut_sin6, method= 'euclidean') 
protdunn_sin6

prot_cl <- mutate(TSnprot, cluster = protcut_sin6)
count(prot_cl,cluster)

ProtClusters <- table(prot_cl$cluster,TSn_cluster)

write.csv(ProtClusters, "ProtClusters.csv")

pdf("Prot dendrogram.pdf", width = 12, height = 8) 
plot(prottreeSIN, main="")
rect.hclust(prottreeSIN, k = 6, border = 2:28)
dev.off()

tiff(file = "Prot dendrogram.tiff", width = 24, height = 16, res = 600, units = "cm")
plot(prottreeSIN, main="")
rect.hclust(prottreeSIN, k = 6, border = 2:28)
dev.off()

##### Tropho-species comparison #####

TS <- read.csv("TS mean macros.csv")
rownames(TS) <- TS[,1]
TSmacros <- TS[2:4]
summary(TS)


TSpal <- viridis_pal(option = "D")(20)

ggtern(TS, aes(x=Lipid,y=Carbohydrate, z=Protein))+
  geom_text(aes(label = as.factor(TS), colour = as.factor(TS)), vjust=-0.40) +
  scale_colour_manual(values=TSpal, name = "Tropho-species") +
  theme_bw() +
  theme_legend_position('tr') + 
  guides(col=FALSE) +
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.7))   + 
  scale_L_continuous(limits=c(.0,.7))  + 
  scale_R_continuous(limits=c(.3,1))

pdf("TS ggtern.pdf", width = 8, height = 8) 
ggtern(TS, aes(x=Lipid,y=Carbohydrate, z=Protein))+
  geom_text(aes(label = as.factor(TS), colour = as.factor(TS)), vjust=-0.40) +
  scale_colour_manual(values=TSpal, name = "Tropho-species") +
  theme_bw() +
  theme_legend_position('tr') + 
  guides(col=FALSE) +
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.7))   + 
  scale_L_continuous(limits=c(.0,.7))  + 
  scale_R_continuous(limits=c(.3,1))
dev.off()

tiff(file = "TS ggtern.tiff", width = 16, height = 16, res = 600, units = "cm")
ggtern(TS, aes(x=Lipid,y=Carbohydrate, z=Protein))+
  geom_text(aes(label = as.factor(TS), colour = as.factor(TS)), vjust=-0.40) +
  scale_colour_manual(values=TSpal, name = "Tropho-species") +
  theme_bw() +
  theme_legend_position('tr') + 
  guides(col=FALSE) +
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.7))   + 
  scale_L_continuous(limits=c(.0,.7))  + 
  scale_R_continuous(limits=c(.3,1))
dev.off()

##### ENNR All #####

ennr <- read.csv("Fam_ENNR_Diet_All.csv")

invertsennr <- read.csv("Fam_ENNR_Inverts.csv")

all.null <- generate_null_net(ennr[,2:79], invertsennr[,2:78],
                              sims = 999, data.type = "names",
                              summary.type = "sum",
                              r.samples = invertsennr[,1],
                              c.samples = ennr[,1])

all.null.tab <- test_interactions(all.null)

all.null.tab

write.csv(all.null.tab, "ALL_SES_ENNR_out.csv")

ennrout <- read.csv("ALL_SES_ENNR_output.csv")

ennrout$SES2 <- ennrout$SES-min(ennrout$SES)
ennrout$SES2 <- scale(ennrout$SES2, center = FALSE)
summary(ennrout$SES2)

write.csv(ennrout, "ALL_SES_ENNR_scaled.csv")

##### ENNR Simulation Extraction and Macronutrient Comparison ####

mcn <- read.csv("Family level macronutrient data.csv")

null.1 <- generate_null_net(ennr[,2:79], invertsennr[,2:78],
                            sims = 999, data.type = "names",
                            summary.type = "none",
                            r.samples = invertsennr[,1],
                            c.samples = ennr[,1])

# Extract the outputs from the null model
rand.nut <- data.frame(null.1$rand.data[, 1:2], 
                       ID = rep(seq(1, nrow(ennr)),
                                length(unique(null.1$rand$Iteration))),
                       null.1$rand.data[, 3:ncol(null.1$rand.data)])

rand.nut <- reshape2::melt(rand.nut, 
                           id.vars = c("Consumer", "Iteration", "ID"))


# Add in the nutrient data
rand.nut <- rand.nut[rand.nut$value > 0, ]
rand.nut <- merge(rand.nut, mcn, by.x = "variable", by.y = "Family")


# Re-arrange/order to tidy up data frame
rand.nut <- rand.nut[order(rand.nut$Consumer, rand.nut$ID, rand.nut$Iteration), ]
rand.nut <- rand.nut[, c(2, 4, 3, 1, 6:8)]
colnames(rand.nut)[4] <- "Resource"


# Summarise for individuals
r.nut.sum <- aggregate(rand.nut[, 5:7], 
                       by = list(rand.nut$Consumer, rand.nut$ID,
                                 rand.nut$Iteration), mean)
colnames(r.nut.sum)[1:3] <- c("Consumer", "ID", "Iteration")



# +++++++++++++++++
# Generate equivalent data for the observed
obs.nut <- data.frame(Consumer = ennr[, 2], ID = seq(1, nrow(ennr)), 
                      ennr[, -c(1:2)])

obs.nut <- reshape2::melt(obs.nut, id.vars = c("Consumer", "ID"))
obs.nut <- obs.nut[obs.nut$value > 0, ]
obs.nut <- merge(obs.nut, mcn, by.x = "variable", by.y = "Family")

obs.nut <- obs.nut[order(obs.nut$ID), ]

obs.nut <- obs.nut[, c(2, 3, 1, 5:7)]
colnames(obs.nut)[3] <- "Resource"


obs.nut.sum <- aggregate(obs.nut[, 4:6], 
                         by = list(obs.nut$Consumer, obs.nut$ID), mean)
colnames(obs.nut.sum)[1:2] <- c("Consumer", "ID")
head(obs.nut.sum)
# +++++++++++++++++




# +++++++++++++++++
# Individual-level comparison
d1 <- aggregate(r.nut.sum[, 4:6], by = list(r.nut.sum$Consumer, r.nut.sum$ID),
                mean)
d2 <- aggregate(r.nut.sum[, 4:6], by = list(r.nut.sum$Consumer, r.nut.sum$ID), 
                quantile, probs = 0.025)
d3 <- aggregate(r.nut.sum[, 4:6], by = list(r.nut.sum$Consumer, r.nut.sum$ID),
                quantile, probs = 0.975)


# Check that the order is retained and combine
identical(d1$Group.2, d2$Group.2)
identical(d1$Group.2, d3$Group.2)

r.nut.indiv <- data.frame(Consumer = d1[, 1], ID = d1[, 2],
                          mean.C = d1[, 3], low.C = d2[, 3], upp.C = d3[, 3],
                          mean.L = d1[, 4], low.L = d2[, 4], upp.L = d3[, 4],
                          mean.P = d1[, 5], low.P = d2[, 5], upp.P = d3[, 5])


identical(as.character(obs.nut.sum$Consumer), 
          as.character(r.nut.indiv$Consumer))
identical(obs.nut.sum$ID, r.nut.indiv$ID)


# Combine observed mean nutrients with means and 95% CIs across the simulations
indiv <- cbind(obs.nut.sum, r.nut.indiv[, -c(1:2)])
head(indiv)

write.csv(indiv, "ALL_NewENNRSimOut.csv")

##### Expected vs Observed #####

evo <- read.csv("ALL_NewENNRSimStack.csv")

mvEVO <- mvabund(evo[,2:4])

EVO.lm<-manylm((mvEVO) ~ Consumer
                 , data=evo)

plot(EVO.lm)

anoEVO <- anova(EVO.lm, p.uni="adjusted")

anoEVO

NRspal <- brewer.pal(3, "Accent")
NRspal <- c("#7FC97F", "#BEAED4")

spinutvsexp = ggtern(evo[2:4], aes(x=Lip,y=Carb, z=Prot, fill=evo$Consumer,shape=evo$Consumer))+
  geom_point(aes(colour = as.factor(evo$Consumer))) +
  scale_shape_manual(values=c(21,24), name = "Consumer", guide = guide_legend(override.aes = list(linetype = c(0))))+
  scale_fill_manual(values=NRspal, name = "Consumer") +
  scale_colour_manual(values=NRspal, name = "Consumer") +
  theme_bw() +
  theme_legend_position('tr') + 
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.6))   + 
  scale_L_continuous(limits=c(.05,.65))  + 
  scale_R_continuous(limits=c(.35,.95)) +
  geom_encircle(alpha=0.2,size=1) 

Evocent = ddply(evo, .(Consumer), summarize, Carb = mean(Carb), Prot = mean(Prot), Lip = mean(Lip))

spinutvsexp = ggtern(evo[2:4], aes(x=Lip,y=Carb, z=Prot, fill=evo$Consumer,shape=evo$Consumer))+
  scale_shape_manual(values=c(25,21), name = "Consumer", guide = guide_legend(override.aes = list(linetype = c(0))))+
  scale_fill_manual(values=NRspal, name = "Consumer") +
  scale_colour_manual(values=NRspal, name = "Consumer") +
  theme_bw() +
  theme_legend_position('tr') +
  geom_encircle(alpha=0.2,size=1)+ 
  geom_point(aes(colour = as.factor(evo$Consumer)), alpha = 0.8) +
  geom_point(data = Evocent, size = 5, alpha = 0.8, shape = c(25,21), color = "black", fill = c("#7FC97F", "#BEAED4")) +
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.6))   + 
  scale_L_continuous(limits=c(.05,.65))  + 
  scale_R_continuous(limits=c(.35,.95))  

print(spinutvsexp)

pdf("Exp vs Obs tern.pdf", width = 10, height = 10) 
print(spinutvsexp)
dev.off()

tiff(file = "Exp vs Obs tern.tiff", width = 20, height = 20, res = 600, units = "cm")
print(spinutvsexp)
dev.off()

##### Comparison of spider nutrient intake #####

NRs <- read.csv("ALL_NewENNRSimOutput.csv")
NRsprey <- NRs[4:6]
summary(NRs)
NRs$Sex <- as.factor(NRs$Sex)
NRs$Life <- as.factor(NRs$Maturity)
NRs$Genus <- as.factor(NRs$Genus)

NRspal <- brewer.pal(5, "Accent")

mvNRsprey <- mvabund((NRs[,4:6]))
meanvar.plot(mvNRsprey)
NRsm1<-manylm(mvNRsprey ~ Sex + Life + Genus, data=NRs)
plot(NRsm1)
anoNRsm1 <- anova(NRsm1, p.uni="adjusted")
anoNRsm1

##### Spider nutrient intake comparison ternary plots #####

evogen <- droplevels(NRs[!NRs$Genus == 'Expected',])

evogen$Genus <- factor(evogen$Genus, levels=c("Bathyphantes", "Erigone",  "Tenuiphantes", "Microlinyphia", "Pardosa"))

NRsgenpal <- viridis_pal(option = "D", direction=1)(5)
NRsgenpal
NRsgenpal <- c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")

spinutvsgenus = ggtern(evogen[4:6], aes(x=Lip,y=Carb, z=Prot, fill=evogen$Genus,shape=evogen$Genus))+
  geom_point(aes(colour = as.factor(evogen$Genus))) +
  scale_shape_manual(values=c(6,15,16,17,18,20), name = "Genus", guide = guide_legend(override.aes = list(linetype = c(0))))+
  scale_fill_manual(values=NRsgenpal, name = "Genus") +
  scale_colour_manual(values=NRsgenpal, name = "Genus") +
  theme_bw() +
  theme_legend_position('tr') + 
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.6))   + 
  scale_L_continuous(limits=c(.05,.65))  + 
  scale_R_continuous(limits=c(.35,.95)) +
  geom_encircle(alpha=0.2,size=1) 

Gencent = ddply(evogen, .(Genus), summarize, Obs.Carb = mean(Obs.Carb), Obs.Prot = mean(Obs.Prot), Obs.Lip = mean(Obs.Lip))

spinutvsgenus = ggtern(evogen[4:6], aes(x=Obs.Lip,y=Obs.Carb, z=Obs.Prot, fill=evogen$Genus,shape=evogen$Genus))+
  scale_shape_manual(values=c(25,22,21,24,23), name = "Genus", guide = guide_legend(override.aes = list(linetype = c(0))))+
  scale_fill_manual(values=NRsgenpal, name = "Genus") +
  scale_colour_manual(values=NRsgenpal, name = "Genus") +
  theme_bw() +
  theme_legend_position('tr') +
  geom_encircle(alpha=0.2,size=1)+ 
  geom_point(aes(colour = as.factor(evogen$Genus)), alpha = 0.8) +
  geom_point(data = Gencent, size = 5, alpha = 0.8, shape = c(25,22,21,24,23), color = "black", fill = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) +
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.6))   + 
  scale_L_continuous(limits=c(.05,.65))  + 
  scale_R_continuous(limits=c(.35,.95))  

print(spinutvsgenus)

pdf("Genus vs nutrient tern.pdf", width = 10, height = 10) 
print(spinutvsgenus)
dev.off()

tiff(file = "Genus vs nutrient tern.tiff", width = 20, height = 20, res = 600, units = "cm")
print(spinutvsgenus)
dev.off()

evosex <- droplevels(NRs[!NRs$Sex == 'N/A',])

evosex$Sex <- factor(evosex$Sex, levels=c("Female", "Male"))

NRssexpal <- brewer.pal(4, "Set1")
NRssexpal<- c("#E41A1C", "#377EB8")

spinutvssex = ggtern(evosex[4:6], aes(x=Lip,y=Carb, z=Prot, fill=evosex$Sex,shape=evosex$Sex))+
  geom_point(aes(colour = as.factor(evosex$Sex))) +
  scale_shape_manual(values=c(15,16,17), name = "Sex", guide = guide_legend(override.aes = list(linetype = c(0))))+
  scale_fill_manual(values=NRssexpal, name = "Sex") +
  scale_colour_manual(values=NRssexpal, name = "Sex") +
  theme_bw() +
  theme_legend_position('tr') + 
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.6))   + 
  scale_L_continuous(limits=c(.05,.65))  + 
  scale_R_continuous(limits=c(.35,.95)) +
  geom_encircle(alpha=0.2,size=1) 

Sexcent = ddply(evosex, .(Sex), summarize, Obs.Carb = mean(Obs.Carb), Obs.Prot = mean(Obs.Prot), Obs.Lip = mean(Obs.Lip))

spinutvssex = ggtern(evosex[4:6], aes(x=Obs.Lip,y=Obs.Carb, z=Obs.Prot, fill=evosex$Sex,shape=evosex$Sex))+
  scale_shape_manual(values=c(21,24), name = "Sex", guide = guide_legend(override.aes = list(linetype = c(0))))+
  scale_fill_manual(values=NRssexpal, name = "Sex") +
  scale_colour_manual(values=NRssexpal, name = "Sex") +
  theme_bw() +
  theme_legend_position('tr') +
  geom_encircle(alpha=0.2,size=1)+ 
  geom_point(aes(colour = as.factor(evosex$Sex)), alpha = 0.8) +
  geom_point(data = Sexcent, size = 5, alpha = 0.8, shape = c(21,24), color = "black", fill = c("#E41A1C", "#377EB8")) +
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.6))   + 
  scale_L_continuous(limits=c(.05,.65))  + 
  scale_R_continuous(limits=c(.35,.95))  

print(spinutvssex)

pdf("Sex vs nutrient tern.pdf", width = 10, height = 10) 
print(spinutvssex)
dev.off()

tiff(file = "Sex vs nutrient tern.tiff", width = 20, height = 20, res = 600, units = "cm")
print(spinutvssex)
dev.off()

evogen$Maturity <- factor(evogen$Maturity, levels=c("Adult", "Juvenile"))

NRsmatpal <- brewer.pal(4, "Set1")
NRsmatpal <- c("#4DAF4A", "#984EA3")

spinutvsmat = ggtern(evogen[4:6], aes(x=Lip,y=Carb, z=Prot, fill=evogen$Maturity,shape=evogen$Maturity))+
  geom_point(aes(colour = as.factor(evogen$Maturity))) +
  scale_shape_manual(values=c(15,16), name = "Maturity", guide = guide_legend(override.aes = list(linetype = c(0))))+
  scale_fill_manual(values=NRsmatpal, name = "Maturity") +
  scale_colour_manual(values=NRsmatpal, name = "Maturity") +
  theme_bw() +
  theme_legend_position('tr') + 
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.6))   + 
  scale_L_continuous(limits=c(.05,.65))  + 
  scale_R_continuous(limits=c(.35,.95)) +
  geom_encircle(alpha=0.2,size=1) 

Matcent = ddply(evogen, .(Maturity), summarize, Obs.Carb = mean(Obs.Carb), Obs.Prot = mean(Obs.Prot), Obs.Lip = mean(Obs.Lip))

spinutvsmat = ggtern(evogen[4:6], aes(x=Obs.Lip,y=Obs.Carb, z=Obs.Prot, fill=evogen$Maturity,shape=evogen$Maturity))+
  scale_shape_manual(values=c(21,24), name = "Maturity", guide = guide_legend(override.aes = list(linetype = c(0))))+
  scale_fill_manual(values=NRsmatpal, name = "Maturity") +
  scale_colour_manual(values=NRsmatpal, name = "Maturity") +
  theme_bw() +
  theme_legend_position('tr') +
  geom_encircle(alpha=0.2,size=1)+ 
  geom_point(aes(colour = as.factor(evogen$Maturity)), alpha = 0.8) +
  geom_point(data = Matcent, size = 5, alpha = 0.8, shape = c(21,24), color = "black", fill = c("#4DAF4A", "#984EA3")) +
  xlab("Lipid") + ylab("Carbohydrate") + zlab("Protein") +
  scale_T_continuous(limits=c(.0,.6))   + 
  scale_L_continuous(limits=c(.05,.65))  + 
  scale_R_continuous(limits=c(.35,.95))  

print(spinutvsmat)

pdf("Maturity vs nutrient tern.pdf", width = 10, height = 10) 
print(spinutvsmat)
dev.off()

tiff(file = "Maturity vs nutrient tern.tiff", width = 20, height = 20, res = 600, units = "cm")
print(spinutvsmat)
dev.off()


##### Analysis of taxonomic or life history based deviation in nutrient-specific foraging from null models #####

mvNRennr <- mvabund(NRs[,13:15])

nrennrm1<-manylm((mvNRennr) ~ Genus + Maturity + Sex +
                 Genus:Maturity + Genus:Sex + Maturity:Sex
                 , data=NRs)

plot(nrennrm1)

anodiff <- anova(nrennrm1, p.uni="adjusted")

anodiff


##### TS ENNR Genus #####

tsennr <- read.csv("TS_ENNR_Diet_Genus.csv")

tsinvertsennr <- read.csv("TS_ENNR_Inverts.csv")

genus.null <- generate_null_net(tsennr[,2:22], tsinvertsennr[,2:21],
                                sims = 999, data.type = "names",
                                summary.type = "sum",
                                r.samples = tsinvertsennr[,1],
                                c.samples = tsennr[,1])

#par(mfrow = c(2,3))
par(mfrow = c(1,1))
plot_preferences(genus.null, "Bathyphantes", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)

pdf("Bathyphantes_TS_ENNR.pdf", width = 6, height = 8) 
plot_preferences(genus.null, "Bathyphantes", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

tiff(file = "Bathyphantes_TS_ENNR.tiff", width = 12, height = 16, res = 600, units = "cm")
plot_preferences(genus.null, "Bathyphantes", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

plot_preferences(genus.null, "Erigone", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)

pdf("Erigone_TS_ENNR.pdf", width = 6, height = 8) 
plot_preferences(genus.null, "Erigone", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

tiff(file = "Erigone_TS_ENNR.tiff", width = 12, height = 16, res = 600, units = "cm")
plot_preferences(genus.null, "Erigone", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

plot_preferences(genus.null, "Tenuiphantes", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)

pdf("Tenuiphantes_TS_ENNR.pdf", width = 6, height = 8) 
plot_preferences(genus.null, "Tenuiphantes", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

tiff(file = "Tenuiphantes_TS_ENNR.tiff", width = 12, height = 16, res = 600, units = "cm")
plot_preferences(genus.null, "Tenuiphantes", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

plot_preferences(genus.null, "Microlinyphia", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)

pdf("Microlinyphia_TS_ENNR.pdf", width = 6, height = 8) 
plot_preferences(genus.null, "Microlinyphia", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

tiff(file = "Microlinyphia_TS_ENNR.tiff", width = 12, height = 16, res = 600, units = "cm")
plot_preferences(genus.null, "Microlinyphia", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

plot_preferences(genus.null, "Pardosa", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)

pdf("Pardosa_TS_ENNR.pdf", width = 6, height = 8) 
plot_preferences(genus.null, "Pardosa", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

tiff(file = "Pardosa_TS_ENNR.tiff", width = 12, height = 16, res = 600, units = "cm")
plot_preferences(genus.null, "Pardosa", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

gen.links <- test_interactions(genus.null, signif.level = 0.95)

gen.links

write.table(gen.links, "Genus TS ENNR.csv")


# Call test_interactions
gbti <- test_interactions(genus.null, signif.level = 0.95)
gbti <- gbti[gbti$Consumer == "Bathyphantes", ]
gbti[, 3] <- ifelse(rowSums(gbti[, 3:6]) == 0, NA, gbti[, 3])
gbti[, 4] <- ifelse(rowSums(gbti[, 3:6]) == 0, NA, gbti[, 4])
gbti[, 5] <- ifelse(rowSums(gbti[, 3:6]) == 0, NA, gbti[, 5])
gbti[, 6] <- ifelse(rowSums(gbti[, 3:6]) == 0, NA, gbti[, 6])
gbti
gbti.t <- gbti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

gbti <- gbti[c(2,3,6,9,11,12,16,19),]


# Set up maximum x-axis value for xlim. Add an additional 5%
gbmin.x <- min(gbti[, 3:6], na.rm = TRUE)
gbmin.x <- max(0, gbmin.x, na.rm = TRUE)
gbmax.x <- max(gbti[, 3:6], na.rm = TRUE)
gbmax.x <- gbmax.x * 1.05
gbti$Setup <- seq(gbmin.x, gbmax.x, length.out = nrow(gbti))

# Plot built up in 2 stages: i) using min and max values to set the
#   y-axis range without having to use ylim (so this can be customised
#   by the user), ii) the main dbarplot and label, and iii) the error
graphics::dotchart(gbti$Setup, labels = paste(gbti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Bathyphantes")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(gbti)){
  eval(parse(text = paste("lines(x = c(gbti$Lower.", 0.95 * 100,
                          ".CL[i], gbti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(gbti$Test[i] == "Weaker") p.col <- res.col[1]
  if(gbti$Test[i] == "ns" | is.na(gbti$Test[i])) p.col <- res.col[2]
  if(gbti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(gbti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}

pdf("Bathyphantes_TS_ENNRsig.pdf", width = 6, height = 6) 
graphics::dotchart(gbti$Setup, labels = paste(gbti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Bathyphantes")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(gbti)){
  eval(parse(text = paste("lines(x = c(gbti$Lower.", 0.95 * 100,
                          ".CL[i], gbti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(gbti$Test[i] == "Weaker") p.col <- res.col[1]
  if(gbti$Test[i] == "ns" | is.na(gbti$Test[i])) p.col <- res.col[2]
  if(gbti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(gbti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

tiff(file = "Bathyphantes_TS_ENNRsig.tiff", width = 12, height = 12, res = 600, units = "cm")
graphics::dotchart(gbti$Setup, labels = paste(gbti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Bathyphantes")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(gbti)){
  eval(parse(text = paste("lines(x = c(gbti$Lower.", 0.95 * 100,
                          ".CL[i], gbti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(gbti$Test[i] == "Weaker") p.col <- res.col[1]
  if(gbti$Test[i] == "ns" | is.na(gbti$Test[i])) p.col <- res.col[2]
  if(gbti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(gbti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()


# Call test_interactions
geti <- test_interactions(genus.null, signif.level = 0.95)
geti <- geti[geti$Consumer == "Erigone", ]
geti[, 3] <- ifelse(rowSums(geti[, 3:6]) == 0, NA, geti[, 3])
geti[, 4] <- ifelse(rowSums(geti[, 3:6]) == 0, NA, geti[, 4])
geti[, 5] <- ifelse(rowSums(geti[, 3:6]) == 0, NA, geti[, 5])
geti[, 6] <- ifelse(rowSums(geti[, 3:6]) == 0, NA, geti[, 6])
geti
geti.t <- geti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

geti <- geti[c(1,2,4,11,12,15,16),]


# Set up maximum x-axis value for xlim. Add an additional 5%
gemin.x <- min(geti[, 3:6], na.rm = TRUE)
gemin.x <- max(0, gemin.x, na.rm = TRUE)
gemax.x <- max(geti[, 3:6], na.rm = TRUE)
gemax.x <- gemax.x * 1.05
geti$Setup <- seq(gemin.x, gemax.x, length.out = nrow(geti))

# Plot built up in 2 stages: i) using min and max values to set the
#   y-axis range without having to use ylim (so this can be customised
#   by the user), ii) the main dbarplot and label, and iii) the error
graphics::dotchart(geti$Setup, labels = paste(geti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Erigone")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(geti)){
  eval(parse(text = paste("lines(x = c(geti$Lower.", 0.95 * 100,
                          ".CL[i], geti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(geti$Test[i] == "Weaker") p.col <- res.col[1]
  if(geti$Test[i] == "ns" | is.na(geti$Test[i])) p.col <- res.col[2]
  if(geti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(geti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}

pdf("Erigone_TS_ENNRsig.pdf", width = 6, height = 6) 
graphics::dotchart(geti$Setup, labels = paste(geti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Erigone")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(geti)){
  eval(parse(text = paste("lines(x = c(geti$Lower.", 0.95 * 100,
                          ".CL[i], geti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(geti$Test[i] == "Weaker") p.col <- res.col[1]
  if(geti$Test[i] == "ns" | is.na(geti$Test[i])) p.col <- res.col[2]
  if(geti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(geti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

tiff(file = "Erigone_TS_ENNRsig.tiff", width = 12, height = 12, res = 600, units = "cm")
graphics::dotchart(geti$Setup, labels = paste(geti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Erigone")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(geti)){
  eval(parse(text = paste("lines(x = c(geti$Lower.", 0.95 * 100,
                          ".CL[i], geti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(geti$Test[i] == "Weaker") p.col <- res.col[1]
  if(geti$Test[i] == "ns" | is.na(geti$Test[i])) p.col <- res.col[2]
  if(geti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(geti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()



# Call test_interactions
gtti <- test_interactions(genus.null, signif.level = 0.95)
gtti <- gtti[gtti$Consumer == "Tenuiphantes", ]
gtti[, 3] <- ifelse(rowSums(gtti[, 3:6]) == 0, NA, gtti[, 3])
gtti[, 4] <- ifelse(rowSums(gtti[, 3:6]) == 0, NA, gtti[, 4])
gtti[, 5] <- ifelse(rowSums(gtti[, 3:6]) == 0, NA, gtti[, 5])
gtti[, 6] <- ifelse(rowSums(gtti[, 3:6]) == 0, NA, gtti[, 6])
gtti
gtti.t <- gtti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

gtti <- gtti[c(2,3,4,5,6,7,11,12,13,15,16,19),]


# Set up maximum x-axis value for xlim. Add an additional 5%
gtmin.x <- min(gtti[, 3:6], na.rm = TRUE)
gtmin.x <- max(0, gtmin.x, na.rm = TRUE)
gtmax.x <- max(gtti[, 3:6], na.rm = TRUE)
gtmax.x <- gtmax.x * 1.05
gtti$Setup <- seq(gtmin.x, gtmax.x, length.out = nrow(gtti))

# Plot built up in 2 stages: i) using min and max values to set the
#   y-axis range without having to use ylim (so this can be customised
#   by the user), ii) the main dbarplot and label, and iii) the error
graphics::dotchart(gtti$Setup, labels = paste(gtti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Tenuiphantes")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(gtti)){
  eval(parse(text = paste("lines(x = c(gtti$Lower.", 0.95 * 100,
                          ".CL[i], gtti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(gtti$Test[i] == "Weaker") p.col <- res.col[1]
  if(gtti$Test[i] == "ns" | is.na(gtti$Test[i])) p.col <- res.col[2]
  if(gtti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(gtti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}

pdf("Tenuiphantes_TS_ENNRsig.pdf", width = 6, height = 8) 
graphics::dotchart(gtti$Setup, labels = paste(gtti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Tenuiphantes")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(gtti)){
  eval(parse(text = paste("lines(x = c(gtti$Lower.", 0.95 * 100,
                          ".CL[i], gtti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(gtti$Test[i] == "Weaker") p.col <- res.col[1]
  if(gtti$Test[i] == "ns" | is.na(gtti$Test[i])) p.col <- res.col[2]
  if(gtti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(gtti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

tiff(file = "Tenuiphantes_TS_ENNRsig.tiff", width = 12, height = 16, res = 600, units = "cm")
graphics::dotchart(gtti$Setup, labels = paste(gtti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Tenuiphantes")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(gtti)){
  eval(parse(text = paste("lines(x = c(gtti$Lower.", 0.95 * 100,
                          ".CL[i], gtti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(gtti$Test[i] == "Weaker") p.col <- res.col[1]
  if(gtti$Test[i] == "ns" | is.na(gtti$Test[i])) p.col <- res.col[2]
  if(gtti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(gtti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

# Call test_interactions
gmti <- test_interactions(genus.null, signif.level = 0.95)
gmti <- gmti[gmti$Consumer == "Microlinyphia", ]
gmti[, 3] <- ifelse(rowSums(gmti[, 3:6]) == 0, NA, gmti[, 3])
gmti[, 4] <- ifelse(rowSums(gmti[, 3:6]) == 0, NA, gmti[, 4])
gmti[, 5] <- ifelse(rowSums(gmti[, 3:6]) == 0, NA, gmti[, 5])
gmti[, 6] <- ifelse(rowSums(gmti[, 3:6]) == 0, NA, gmti[, 6])
gmti
gmti.t <- gmti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

gmti <- gmti[c(1,3,4,5,16),]


# Set up maximum x-axis value for xlim. Add an additional 5%
gmmin.x <- min(gmti[, 3:6], na.rm = TRUE)
gmmin.x <- max(0, gmmin.x, na.rm = TRUE)
gmmax.x <- max(gmti[, 3:6], na.rm = TRUE)
gmmax.x <- gmmax.x * 1.05
gmti$Setup <- seq(gmmin.x, gmmax.x, length.out = nrow(gmti))

# Plot built up in 2 stages: i) using min and max values to set the
#   y-axis range without having to use ylim (so this can be customised
#   by the user), ii) the main dbarplot and label, and iii) the error
graphics::dotchart(gmti$Setup, labels = paste(gmti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Microlinyphia")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(gmti)){
  eval(parse(text = paste("lines(x = c(gmti$Lower.", 0.95 * 100,
                          ".CL[i], gmti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(gmti$Test[i] == "Weaker") p.col <- res.col[1]
  if(gmti$Test[i] == "ns" | is.na(gmti$Test[i])) p.col <- res.col[2]
  if(gmti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(gmti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}

pdf("Microlinyphia_TS_ENNRsig.pdf", width = 6, height = 4) 
graphics::dotchart(gmti$Setup, labels = paste(gmti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Microlinyphia")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(gmti)){
  eval(parse(text = paste("lines(x = c(gmti$Lower.", 0.95 * 100,
                          ".CL[i], gmti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(gmti$Test[i] == "Weaker") p.col <- res.col[1]
  if(gmti$Test[i] == "ns" | is.na(gmti$Test[i])) p.col <- res.col[2]
  if(gmti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(gmti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

tiff(file = "Microlinyphia_TS_ENNRsig.tiff", width = 12, height = 8, res = 600, units = "cm")
graphics::dotchart(gmti$Setup, labels = paste(gmti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Microlinyphia")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(gmti)){
  eval(parse(text = paste("lines(x = c(gmti$Lower.", 0.95 * 100,
                          ".CL[i], gmti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(gmti$Test[i] == "Weaker") p.col <- res.col[1]
  if(gmti$Test[i] == "ns" | is.na(gmti$Test[i])) p.col <- res.col[2]
  if(gmti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(gmti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()


# Call test_interactions
gpti <- test_interactions(genus.null, signif.level = 0.95)
gpti <- gpti[gpti$Consumer == "Pardosa", ]
gpti[, 3] <- ifelse(rowSums(gpti[, 3:6]) == 0, NA, gpti[, 3])
gpti[, 4] <- ifelse(rowSums(gpti[, 3:6]) == 0, NA, gpti[, 4])
gpti[, 5] <- ifelse(rowSums(gpti[, 3:6]) == 0, NA, gpti[, 5])
gpti[, 6] <- ifelse(rowSums(gpti[, 3:6]) == 0, NA, gpti[, 6])
gpti
gpti.t <- gpti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

gpti <- gpti[c(2,5,15,16),]


# Set up maximum x-axis value for xlim. Add an additional 5%
gpmin.x <- min(gpti[, 3:6], na.rm = TRUE)
gpmin.x <- max(0, gpmin.x, na.rm = TRUE)
gpmax.x <- max(gpti[, 3:6], na.rm = TRUE)
gpmax.x <- gpmax.x * 1.05
gpti$Setup <- seq(gpmin.x, gpmax.x, length.out = nrow(gpti))

# Plot built up in 2 stages: i) using min and max values to set the
#   y-axis range without having to use ylim (so this can be customised
#   by the user), ii) the main dbarplot and label, and iii) the error
graphics::dotchart(gpti$Setup, labels = paste(gpti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Pardosa")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(gpti)){
  eval(parse(text = paste("lines(x = c(gpti$Lower.", 0.95 * 100,
                          ".CL[i], gpti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(gpti$Test[i] == "Weaker") p.col <- res.col[1]
  if(gpti$Test[i] == "ns" | is.na(gpti$Test[i])) p.col <- res.col[2]
  if(gpti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(gpti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}

pdf("Pardosa_TS_ENNRsig.pdf", width = 6, height = 4) 
graphics::dotchart(gpti$Setup, labels = paste(gpti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Pardosa")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(gpti)){
  eval(parse(text = paste("lines(x = c(gpti$Lower.", 0.95 * 100,
                          ".CL[i], gpti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(gpti$Test[i] == "Weaker") p.col <- res.col[1]
  if(gpti$Test[i] == "ns" | is.na(gpti$Test[i])) p.col <- res.col[2]
  if(gpti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(gpti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

tiff(file = "Pardosa_TS_ENNRsig.tiff", width = 12, height = 8, res = 600, units = "cm")
graphics::dotchart(gpti$Setup, labels = paste(gpti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Pardosa")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(gpti)){
  eval(parse(text = paste("lines(x = c(gpti$Lower.", 0.95 * 100,
                          ".CL[i], gpti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(gpti$Test[i] == "Weaker") p.col <- res.col[1]
  if(gpti$Test[i] == "ns" | is.na(gpti$Test[i])) p.col <- res.col[2]
  if(gpti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(gpti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

##### Tropho-ternary Genus #####

TS <- read.csv("TS mean macros.csv")
TS$TS <- as.factor(TS$TS)
rownames(TS) <- TS[,1]
summary(TS)

nrow(gbti.t)

dgb <- merge(gbti.t, TS, by.x = "Resource", by.y = "TS")

dgb

bathyphantes_tschoicetern <- ggtern(dgb, aes(x = Lipid, y = Carbohydrate, z = Protein)) +
  geom_point(aes(size = abs(SES),colour = Test)) +
  scale_colour_manual(breaks = c("Stronger", "ns", "Weaker"),
                      values = c("#fc8d59", "#ffffbf", "#91bfdb")) +
  theme_dark()+
  scale_T_continuous(limits=c(.0,.8))   + 
  scale_L_continuous(limits=c(.0,.8))  + 
  scale_R_continuous(limits=c(.2,1))

bathyphantes_tschoicetern

pdf("Bathyphantes_TS_ENNR_Tern.pdf", width = 8, height = 6) 
bathyphantes_tschoicetern
dev.off()

tiff(file = "Bathyphantes_TS_ENNR_Tern.tiff", width = 16, height = 12, res = 600, units = "cm")
bathyphantes_tschoicetern
dev.off()

nrow(geti.t)

dge <- merge(geti.t, TS, by.x = "Resource", by.y = "TS")

dge

erigone_tschoicetern <- ggtern(dge, aes(x = Lipid, y = Carbohydrate, z = Protein)) +
  geom_point(aes(size = abs(SES),colour = Test)) +
  scale_colour_manual(breaks = c("Stronger", "ns", "Weaker"),
                      values = c("#fc8d59", "#ffffbf", "#91bfdb")) +
  theme_dark()+
  scale_T_continuous(limits=c(.0,.8))   + 
  scale_L_continuous(limits=c(.0,.8))  + 
  scale_R_continuous(limits=c(.2,1))

erigone_tschoicetern

pdf("Erigone_TS_ENNR_Tern.pdf", width = 8, height = 6) 
erigone_tschoicetern
dev.off()

tiff(file = "Erigone_TS_ENNR_Tern.tiff", width = 16, height = 12, res = 600, units = "cm")
erigone_tschoicetern
dev.off()

nrow(gtti.t)

dgt <- merge(gtti.t, TS, by.x = "Resource", by.y = "TS")

dgt

tenuiphantes_tschoicetern <- ggtern(dgt, aes(x = Lipid, y = Carbohydrate, z = Protein)) +
  geom_point(aes(size = abs(SES),colour = Test)) +
  scale_colour_manual(breaks = c("Stronger", "ns", "Weaker"),
                      values = c("#fc8d59", "#ffffbf", "#91bfdb")) +
  theme_dark()+
  scale_T_continuous(limits=c(.0,.8))   + 
  scale_L_continuous(limits=c(.0,.8))  + 
  scale_R_continuous(limits=c(.2,1))

tenuiphantes_tschoicetern

pdf("Tenuiphantes_TS_ENNR_Tern.pdf", width = 8, height = 6) 
tenuiphantes_tschoicetern
dev.off()

tiff(file = "Tenuiphantes_TS_ENNR_Tern.tiff", width = 16, height = 12, res = 600, units = "cm")
tenuiphantes_tschoicetern
dev.off()

nrow(gmti.t)

dgm <- merge(gmti.t, TS, by.x = "Resource", by.y = "TS")

dgm

microlinyphia_tschoicetern <- ggtern(dgm, aes(x = Lipid, y = Carbohydrate, z = Protein)) +
  geom_point(aes(size = abs(SES),colour = Test)) +
  scale_colour_manual(breaks = c("Stronger", "ns", "Weaker"),
                      values = c("#fc8d59", "#ffffbf", "#91bfdb")) +
  theme_dark()+
  scale_T_continuous(limits=c(.0,.8))   + 
  scale_L_continuous(limits=c(.0,.8))  + 
  scale_R_continuous(limits=c(.2,1))

microlinyphia_tschoicetern

pdf("Microlinyphia_TS_ENNR_Tern.pdf", width = 8, height = 6) 
microlinyphia_tschoicetern
dev.off()

tiff(file = "Microlinyphia_TS_ENNR_Tern.tiff", width = 16, height = 12, res = 600, units = "cm")
microlinyphia_tschoicetern
dev.off()

nrow(gpti.t)

dgp <- merge(gpti.t, TS, by.x = "Resource", by.y = "TS")

dgp

pardosa_tschoicetern <- ggtern(dgp, aes(x = Lipid, y = Carbohydrate, z = Protein)) +
  geom_point(aes(size = abs(SES),colour = Test)) +
  scale_colour_manual(breaks = c("Stronger", "ns", "Weaker"),
                      values = c("#fc8d59", "#ffffbf", "#91bfdb")) +
  theme_dark()+
  scale_T_continuous(limits=c(.0,.8))   + 
  scale_L_continuous(limits=c(.0,.8))  + 
  scale_R_continuous(limits=c(.2,1))

pardosa_tschoicetern

pdf("Pardosa_TS_ENNR_Tern.pdf", width = 8, height = 6) 
pardosa_tschoicetern
dev.off()

tiff(file = "Pardosa_TS_ENNR_Tern.tiff", width = 16, height = 12, res = 600, units = "cm")
pardosa_tschoicetern
dev.off()

##### TS ENNR Sex #####

sextsennr <- read.csv("TS_ENNR_Diet_Sex.csv")

tsinvertsennr <- read.csv("TS_ENNR_Inverts.csv")

sex.null <- generate_null_net(sextsennr[,2:22], tsinvertsennr[,2:21],
                              sims = 999, data.type = "names",
                              summary.type = "sum",
                              r.samples = tsinvertsennr[,1],
                              c.samples = sextsennr[,1])

#par(mfrow = c(1,2))
par(mfrow = c(1,1))
plot_preferences(sex.null, "Female", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)

pdf("Female_TS_ENNR.pdf", width = 6, height = 8) 
plot_preferences(sex.null, "Female", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

tiff(file = "Female_TS_ENNR.tiff", width = 12, height = 16, res = 600, units = "cm")
plot_preferences(sex.null, "Female", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

plot_preferences(sex.null, "Male", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)

pdf("Male_TS_ENNR.pdf", width = 6, height = 8) 
plot_preferences(sex.null, "Male", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

tiff(file = "Male_TS_ENNR.tiff", width = 12, height = 16, res = 600, units = "cm")
plot_preferences(sex.null, "Male", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()


sex.links <- test_interactions(sex.null, signif.level = 0.95)

sex.links

write.table(sex.links, "Sex TS ENNR.csv")

# Call test_interactions
sfti <- test_interactions(sex.null, signif.level = 0.95)
sfti <- sfti[sfti$Consumer == "Female", ]
sfti[, 3] <- ifelse(rowSums(sfti[, 3:6]) == 0, NA, sfti[, 3])
sfti[, 4] <- ifelse(rowSums(sfti[, 3:6]) == 0, NA, sfti[, 4])
sfti[, 5] <- ifelse(rowSums(sfti[, 3:6]) == 0, NA, sfti[, 5])
sfti[, 6] <- ifelse(rowSums(sfti[, 3:6]) == 0, NA, sfti[, 6])
sfti
sfti.t <- sfti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

sfti <- sfti[c(1,2,3,4,5,6,7,11,12,15,16),]


# Set up maximum x-axis value for xlim. Add an additional 5%
sfmin.x <- min(sfti[, 3:6], na.rm = TRUE)
sfmin.x <- max(0, sfmin.x, na.rm = TRUE)
sfmax.x <- max(sfti[, 3:6], na.rm = TRUE)
sfmax.x <- sfmax.x * 1.05
sfti$Setup <- seq(sfmin.x, sfmax.x, length.out = nrow(sfti))

# Plot built up in 2 stages: i) using min and max values to set the
#   y-axis range without having to use ylim (so this can be customised
#   by the user), ii) the main dbarplot and label, and iii) the error
graphics::dotchart(sfti$Setup, labels = paste(sfti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Female")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(sfti)){
  eval(parse(text = paste("lines(x = c(sfti$Lower.", 0.95 * 100,
                          ".CL[i], sfti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(sfti$Test[i] == "Weaker") p.col <- res.col[1]
  if(sfti$Test[i] == "ns" | is.na(sfti$Test[i])) p.col <- res.col[2]
  if(sfti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(sfti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}

pdf("Female_TS_ENNRsig.pdf", width = 6, height = 7) 
graphics::dotchart(sfti$Setup, labels = paste(sfti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Female")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(sfti)){
  eval(parse(text = paste("lines(x = c(sfti$Lower.", 0.95 * 100,
                          ".CL[i], sfti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(sfti$Test[i] == "Weaker") p.col <- res.col[1]
  if(sfti$Test[i] == "ns" | is.na(sfti$Test[i])) p.col <- res.col[2]
  if(sfti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(sfti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

tiff(file = "Female_TS_ENNRsig.tiff", width = 12, height = 14, res = 600, units = "cm")
graphics::dotchart(sfti$Setup, labels = paste(sfti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Female")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(sfti)){
  eval(parse(text = paste("lines(x = c(sfti$Lower.", 0.95 * 100,
                          ".CL[i], sfti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(sfti$Test[i] == "Weaker") p.col <- res.col[1]
  if(sfti$Test[i] == "ns" | is.na(sfti$Test[i])) p.col <- res.col[2]
  if(sfti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(sfti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

# Call test_interactions
smti <- test_interactions(sex.null, signif.level = 0.95)
smti <- smti[smti$Consumer == "Male", ]
smti[, 3] <- ifelse(rowSums(smti[, 3:6]) == 0, NA, smti[, 3])
smti[, 4] <- ifelse(rowSums(smti[, 3:6]) == 0, NA, smti[, 4])
smti[, 5] <- ifelse(rowSums(smti[, 3:6]) == 0, NA, smti[, 5])
smti[, 6] <- ifelse(rowSums(smti[, 3:6]) == 0, NA, smti[, 6])
smti
smti.t <- smti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

smti <- smti[c(2,4,5,6,7,11,12,13,15,16,19),]


# Set up maximum x-axis value for xlim. Add an additional 5%
smmin.x <- min(smti[, 3:6], na.rm = TRUE)
smmin.x <- max(0, smmin.x, na.rm = TRUE)
smmax.x <- max(smti[, 3:6], na.rm = TRUE)
smmax.x <- smmax.x * 1.05
smti$Setup <- seq(smmin.x, smmax.x, length.out = nrow(smti))

# Plot built up in 2 stages: i) using min and max values to set the
#   y-axis range without having to use ylim (so this can be customised
#   by the user), ii) the main dbarplot and label, and iii) the error
graphics::dotchart(smti$Setup, labels = paste(smti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Male")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(smti)){
  eval(parse(text = paste("lines(x = c(smti$Lower.", 0.95 * 100,
                          ".CL[i], smti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(smti$Test[i] == "Weaker") p.col <- res.col[1]
  if(smti$Test[i] == "ns" | is.na(smti$Test[i])) p.col <- res.col[2]
  if(smti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(smti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}


pdf("Male_TS_ENNRsig.pdf", width = 6, height = 7) 
graphics::dotchart(smti$Setup, labels = paste(smti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Male")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(smti)){
  eval(parse(text = paste("lines(x = c(smti$Lower.", 0.95 * 100,
                          ".CL[i], smti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(smti$Test[i] == "Weaker") p.col <- res.col[1]
  if(smti$Test[i] == "ns" | is.na(smti$Test[i])) p.col <- res.col[2]
  if(smti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(smti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

tiff(file = "Male_TS_ENNRsig.tiff", width = 12, height = 14, res = 600, units = "cm")
graphics::dotchart(smti$Setup, labels = paste(smti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Male")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(smti)){
  eval(parse(text = paste("lines(x = c(smti$Lower.", 0.95 * 100,
                          ".CL[i], smti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(smti$Test[i] == "Weaker") p.col <- res.col[1]
  if(smti$Test[i] == "ns" | is.na(smti$Test[i])) p.col <- res.col[2]
  if(smti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(smti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

##### Tropho-ternary Sex #####

nrow(sfti.t)

dsf <- merge(sfti.t, TS, by.x = "Resource", by.y = "TS")

dsf

female_tschoicetern <- ggtern(dsf, aes(x = Lipid, y = Carbohydrate, z = Protein)) +
  geom_point(aes(size = abs(SES),colour = Test)) +
  scale_colour_manual(breaks = c("Stronger", "ns", "Weaker"),
                      values = c("#fc8d59", "#ffffbf", "#91bfdb")) +
  theme_dark()+
  scale_T_continuous(limits=c(.0,.8))   + 
  scale_L_continuous(limits=c(.0,.8))  + 
  scale_R_continuous(limits=c(.2,1))

female_tschoicetern

pdf("Female_TS_ENNR_Tern.pdf", width = 8, height = 6) 
female_tschoicetern
dev.off()

tiff(file = "Female_TS_ENNR_Tern.tiff", width = 16, height = 12, res = 600, units = "cm")
female_tschoicetern
dev.off()

nrow(smti.t)

dsm <- merge(smti.t, TS, by.x = "Resource", by.y = "TS")

dsm

male_tschoicetern <- ggtern(dsm, aes(x = Lipid, y = Carbohydrate, z = Protein)) +
  geom_point(aes(size = abs(SES),colour = Test)) +
  scale_colour_manual(breaks = c("Stronger", "ns", "Weaker"),
                      values = c("#fc8d59", "#ffffbf", "#91bfdb")) +
  theme_dark()+
  scale_T_continuous(limits=c(.0,.8))   + 
  scale_L_continuous(limits=c(.0,.8))  + 
  scale_R_continuous(limits=c(.2,1))

male_tschoicetern

pdf("Male_TS_ENNR_Tern.pdf", width = 8, height = 6) 
male_tschoicetern
dev.off()

tiff(file = "Male_TS_ENNR_Tern.tiff", width = 16, height = 12, res = 600, units = "cm")
male_tschoicetern
dev.off()

##### TS ENNR Life #####

lifetsennr <- read.csv("TS_ENNR_Diet_Life.csv")

tsinvertsennr <- read.csv("TS_ENNR_Inverts.csv")

life.null <- generate_null_net(lifetsennr[,2:22], tsinvertsennr[,2:21],
                               sims = 999, data.type = "names",
                               summary.type = "sum",
                               r.samples = tsinvertsennr[,1],
                               c.samples = lifetsennr[,1])

#par(mfrow = c(1,2))
par(mfrow = c(1,1))
plot_preferences(life.null, "Adult", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)

pdf("Adult_TS_ENNR.pdf", width = 6, height = 8) 
plot_preferences(life.null, "Adult", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

tiff(file = "Adult_TS_ENNR.tiff", width = 12, height = 16, res = 600, units = "cm")
plot_preferences(life.null, "Adult", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

plot_preferences(life.null, "Juvenile", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)

pdf("Juvenile_TS_ENNR.pdf", width = 6, height = 8) 
plot_preferences(life.null, "Juvenile", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()

tiff(file = "Juvenile_TS_ENNR.tiff", width = 12, height = 16, res = 600, units = "cm")
plot_preferences(life.null, "Juvenile", signif.level = 0.95, type = "counts",
                 xlab = "Num. of prey detections", p.cex = 1.5, l.cex = 0.8, lwd = 2)
dev.off()


life.links <- test_interactions(life.null, signif.level = 0.95)

life.links

write.table(life.links, "Life TS ENNR.csv")


# Call test_interactions
lati <- test_interactions(life.null, signif.level = 0.95)
lati <- lati[lati$Consumer == "Adult", ]
lati[, 3] <- ifelse(rowSums(lati[, 3:6]) == 0, NA, lati[, 3])
lati[, 4] <- ifelse(rowSums(lati[, 3:6]) == 0, NA, lati[, 4])
lati[, 5] <- ifelse(rowSums(lati[, 3:6]) == 0, NA, lati[, 5])
lati[, 6] <- ifelse(rowSums(lati[, 3:6]) == 0, NA, lati[, 6])
lati
lati.t <- lati

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

lati <- lati[c(1,2,4,5,6,7,11,12,13,15,16,19),]


# Set up maximum x-axis value for xlim. Add an additional 5%
lamin.x <- min(lati[, 3:6], na.rm = TRUE)
lamin.x <- max(0, lamin.x, na.rm = TRUE)
lamax.x <- max(lati[, 3:6], na.rm = TRUE)
lamax.x <- lamax.x * 1.05
lati$Setup <- seq(lamin.x, lamax.x, length.out = nrow(lati))

# Plot built up in 2 stages: i) using min and max values to set the
#   y-axis range without having to use ylim (so this can be customised
#   by the user), ii) the main dbarplot and label, and iii) the error
graphics::dotchart(lati$Setup, labels = paste(lati$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Adult")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(lati)){
  eval(parse(text = paste("lines(x = c(lati$Lower.", 0.95 * 100,
                          ".CL[i], lati$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(lati$Test[i] == "Weaker") p.col <- res.col[1]
  if(lati$Test[i] == "ns" | is.na(lati$Test[i])) p.col <- res.col[2]
  if(lati$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(lati$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}

pdf("Adult_TS_ENNRsig.pdf", width = 6, height = 8) 
graphics::dotchart(lati$Setup, labels = paste(lati$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Adult")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(lati)){
  eval(parse(text = paste("lines(x = c(lati$Lower.", 0.95 * 100,
                          ".CL[i], lati$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(lati$Test[i] == "Weaker") p.col <- res.col[1]
  if(lati$Test[i] == "ns" | is.na(lati$Test[i])) p.col <- res.col[2]
  if(lati$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(lati$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

tiff(file = "Adult_TS_ENNRsig.tiff", width = 12, height = 16, res = 600, units = "cm")
graphics::dotchart(lati$Setup, labels = paste(lati$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Adult")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(lati)){
  eval(parse(text = paste("lines(x = c(lati$Lower.", 0.95 * 100,
                          ".CL[i], lati$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(lati$Test[i] == "Weaker") p.col <- res.col[1]
  if(lati$Test[i] == "ns" | is.na(lati$Test[i])) p.col <- res.col[2]
  if(lati$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(lati$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

# Call test_interactions
ljti <- test_interactions(life.null, signif.level = 0.95)
ljti <- ljti[ljti$Consumer == "Juvenile", ]
ljti[, 3] <- ifelse(rowSums(ljti[, 3:6]) == 0, NA, ljti[, 3])
ljti[, 4] <- ifelse(rowSums(ljti[, 3:6]) == 0, NA, ljti[, 4])
ljti[, 5] <- ifelse(rowSums(ljti[, 3:6]) == 0, NA, ljti[, 5])
ljti[, 6] <- ifelse(rowSums(ljti[, 3:6]) == 0, NA, ljti[, 6])
ljti
ljti.t <- ljti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

ljti <- ljti[c(4,5,6,11,12,13,16),]


# Set up maximum x-axis value for xlim. Add an additional 5%
ljmin.x <- min(ljti[, 3:6], na.rm = TRUE)
ljmin.x <- max(0, ljmin.x, na.rm = TRUE)
ljmax.x <- max(ljti[, 3:6], na.rm = TRUE)
ljmax.x <- ljmax.x * 1.05
ljti$Setup <- seq(ljmin.x, ljmax.x, length.out = nrow(ljti))

# Plot built up in 2 stages: i) using min and max values to set the
#   y-axis range without having to use ylim (so this can be customised
#   by the user), ii) the main dbarplot and label, and iii) the error
graphics::dotchart(ljti$Setup, labels = paste(ljti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Juvenile")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(ljti)){
  eval(parse(text = paste("lines(x = c(ljti$Lower.", 0.95 * 100,
                          ".CL[i], ljti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(ljti$Test[i] == "Weaker") p.col <- res.col[1]
  if(ljti$Test[i] == "ns" | is.na(ljti$Test[i])) p.col <- res.col[2]
  if(ljti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(ljti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}


pdf("Juvenile_TS_ENNRsig.pdf", width = 6, height = 6) 
graphics::dotchart(ljti$Setup, labels = paste(ljti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Juvenile")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(ljti)){
  eval(parse(text = paste("lines(x = c(ljti$Lower.", 0.95 * 100,
                          ".CL[i], ljti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(ljti$Test[i] == "Weaker") p.col <- res.col[1]
  if(ljti$Test[i] == "ns" | is.na(ljti$Test[i])) p.col <- res.col[2]
  if(ljti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(ljti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

tiff(file = "Juvenile_TS_ENNRsig.tiff", width = 12, height = 12, res = 600, units = "cm")
graphics::dotchart(ljti$Setup, labels = paste(ljti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Juvenile")
graphics::abline(v = 0, lty = 2, col = "dimgrey")


res.col <- c("#67A9CF", "#F7F7F7", "#EF8A62")

for (i in 1:nrow(ljti)){
  eval(parse(text = paste("lines(x = c(ljti$Lower.", 0.95 * 100,
                          ".CL[i], ljti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i))", sep = "")))
  if(ljti$Test[i] == "Weaker") p.col <- res.col[1]
  if(ljti$Test[i] == "ns" | is.na(ljti$Test[i])) p.col <- res.col[2]
  if(ljti$Test[i] == "Stronger") p.col <- res.col[3]
  graphics::points(ljti$Observed[i], i, pch = 21, col = "black",
                   bg = p.col, cex = 2)
}
dev.off()

##### Tropho-ternary Life #####

nrow(lati.t)

dla <- merge(lati.t, TS, by.x = "Resource", by.y = "TS")

dla

adult_tschoicetern <- ggtern(dla, aes(x = Lipid, y = Carbohydrate, z = Protein)) +
  geom_point(aes(size = abs(SES),colour = Test)) +
  scale_colour_manual(breaks = c("Stronger", "ns", "Weaker"),
                      values = c("#fc8d59", "#ffffbf", "#91bfdb")) +
  theme_dark()+
  scale_T_continuous(limits=c(.0,.8))   + 
  scale_L_continuous(limits=c(.0,.8))  + 
  scale_R_continuous(limits=c(.2,1))

adult_tschoicetern

pdf("Adult_TS_ENNR_Tern.pdf", width = 8, height = 6) 
adult_tschoicetern
dev.off()

tiff(file = "Adult_TS_ENNR_Tern.tiff", width = 16, height = 12, res = 600, units = "cm")
adult_tschoicetern
dev.off()

nrow(ljti.t)

dlj <- merge(ljti.t, TS, by.x = "Resource", by.y = "TS")

dlj

juv_tschoicetern <- ggtern(dlj, aes(x = Lipid, y = Carbohydrate, z = Protein)) +
  geom_point(aes(size = abs(SES),colour = Test)) +
  scale_colour_manual(breaks = c("Stronger", "ns", "Weaker"),
                      values = c("#fc8d59", "#ffffbf", "#91bfdb")) +
  theme_dark()+
  scale_T_continuous(limits=c(.0,.8))   + 
  scale_L_continuous(limits=c(.0,.8))  + 
  scale_R_continuous(limits=c(.2,1))

juv_tschoicetern

pdf("Juvenile_TS_ENNR_Tern.pdf", width = 8, height = 6) 
juv_tschoicetern
dev.off()

tiff(file = "Juvenile_TS_ENNR_Tern.tiff", width = 16, height = 12, res = 600, units = "cm")
juv_tschoicetern
dev.off()

##### TS ENNR All #####

tsall_ennr <- read.csv("TS_ENNR_Diet_ALL.csv")

invertsennr <- read.csv("TS_ENNR_Inverts.csv")

tsall.null <- generate_null_net(tsall_ennr[,2:22], invertsennr[,2:21],
                              sims = 999, data.type = "names",
                              summary.type = "sum",
                              r.samples = invertsennr[,1],
                              c.samples = ennr[,1])

tsall.null.tab <- test_interactions(tsall.null)

tsall.null.tab

write.csv(tsall.null.tab, "TS_ALL_SES_ENNR.csv")

##### TS ENNR PerMANOVA #####

tsennrperm <- read.csv("TS_ALL_SES_ENNRforperm.csv")

tsennrperm$Genus <- as.factor(tsennrperm$Genus)
tsennrperm$Maturity <- as.factor(tsennrperm$Maturity)
tsennrperm$Sex <- as.factor(tsennrperm$Sex)

tsennrvals <- (tsennrperm[,4:23])

tspermennr <- adonis(tsennrvals ~ Sex + Genus + Maturity  
                   , data=tsennrperm, permutations = 9999, method = "euclidean")

tspermennr

tssimgen <- simper(tsennrvals, tsennrperm$Genus, permutations=9999)

summary(tssimgen)            

tssimgen    

tssimmat <- simper(tsennrvals, tsennrperm$Maturity, permutations=9999)

summary(tssimmat)            

tssimmat 

tssimsex <- simper(tsennrvals, tsennrperm$Sex, permutations=9999)

summary(tssimsex)            

tssimsex  


##### Genus differences TS ENNR #####

# Call test_interactions
gbti <- test_interactions(genus.null, signif.level = 0.95)
gbti <- gbti[gbti$Consumer == "Bathyphantes", ]
gbti[, 3] <- ifelse(rowSums(gbti[, 3:6]) == 0, NA, gbti[, 3])
gbti[, 4] <- ifelse(rowSums(gbti[, 3:6]) == 0, NA, gbti[, 4])
gbti[, 5] <- ifelse(rowSums(gbti[, 3:6]) == 0, NA, gbti[, 5])
gbti[, 6] <- ifelse(rowSums(gbti[, 3:6]) == 0, NA, gbti[, 6])
gbti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

gbti <- gbti[c(1,2,6,10,11,13,14,15,17,18,19),]


# Set up maximum x-axis value for xlim. Add an additional 5%
gbmin.x <- min(gbti[, 3:6], na.rm = TRUE)
gbmin.x <- max(0, gbmin.x, na.rm = TRUE)
gbmax.x <- max(gbti[, 3:6], na.rm = TRUE)
gbmax.x <- gbmax.x * 1.05
gbti$Setup <- seq(gbmin.x, gbmax.x, length.out = nrow(gbti))

# Call test_interactions
geti <- test_interactions(genus.null, signif.level = 0.95)
geti <- geti[geti$Consumer == "Erigone", ]
geti[, 3] <- ifelse(rowSums(geti[, 3:6]) == 0, NA, geti[, 3])
geti[, 4] <- ifelse(rowSums(geti[, 3:6]) == 0, NA, geti[, 4])
geti[, 5] <- ifelse(rowSums(geti[, 3:6]) == 0, NA, geti[, 5])
geti[, 6] <- ifelse(rowSums(geti[, 3:6]) == 0, NA, geti[, 6])
geti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

geti <- geti[c(1,2,6,10,11,13,14,15,17,18,19),]


# Set up maximum x-axis value for xlim. Add an additional 5%
gemin.x <- min(geti[, 3:6], na.rm = TRUE)
gemin.x <- max(0, gemin.x, na.rm = TRUE)
gemax.x <- max(geti[, 3:6], na.rm = TRUE)
gemax.x <- gemax.x * 1.05
geti$Setup <- seq(gemin.x, gemax.x, length.out = nrow(geti))

# Call test_interactions
gtti <- test_interactions(genus.null, signif.level = 0.95)
gtti <- gtti[gtti$Consumer == "Tenuiphantes", ]
gtti[, 3] <- ifelse(rowSums(gtti[, 3:6]) == 0, NA, gtti[, 3])
gtti[, 4] <- ifelse(rowSums(gtti[, 3:6]) == 0, NA, gtti[, 4])
gtti[, 5] <- ifelse(rowSums(gtti[, 3:6]) == 0, NA, gtti[, 5])
gtti[, 6] <- ifelse(rowSums(gtti[, 3:6]) == 0, NA, gtti[, 6])
gtti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

gtti <- gtti[c(1,2,6,10,11,13,14,15,17,18,19),]


# Set up maximum x-axis value for xlim. Add an additional 5%
gtmin.x <- min(gtti[, 3:6], na.rm = TRUE)
gtmin.x <- max(0, gtmin.x, na.rm = TRUE)
gtmax.x <- max(gtti[, 3:6], na.rm = TRUE)
gtmax.x <- gtmax.x * 1.05
gtti$Setup <- seq(gtmin.x, gtmax.x, length.out = nrow(gtti))

# Call test_interactions
gmti <- test_interactions(genus.null, signif.level = 0.95)
gmti <- gmti[gmti$Consumer == "Microlinyphia", ]
gmti[, 3] <- ifelse(rowSums(gmti[, 3:6]) == 0, NA, gmti[, 3])
gmti[, 4] <- ifelse(rowSums(gmti[, 3:6]) == 0, NA, gmti[, 4])
gmti[, 5] <- ifelse(rowSums(gmti[, 3:6]) == 0, NA, gmti[, 5])
gmti[, 6] <- ifelse(rowSums(gmti[, 3:6]) == 0, NA, gmti[, 6])
gmti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

gmti <- gmti[c(1,2,6,10,11,13,14,15,17,18,19),]


# Set up maximum x-axis value for xlim. Add an additional 5%
gmmin.x <- min(gmti[, 3:6], na.rm = TRUE)
gmmin.x <- max(0, gmmin.x, na.rm = TRUE)
gmmax.x <- max(gmti[, 3:6], na.rm = TRUE)
gmmax.x <- gmmax.x * 1.05
gmti$Setup <- seq(gmmin.x, gmmax.x, length.out = nrow(gmti))

# Call test_interactions
gpti <- test_interactions(genus.null, signif.level = 0.95)
gpti <- gpti[gpti$Consumer == "Pardosa", ]
gpti[, 3] <- ifelse(rowSums(gpti[, 3:6]) == 0, NA, gpti[, 3])
gpti[, 4] <- ifelse(rowSums(gpti[, 3:6]) == 0, NA, gpti[, 4])
gpti[, 5] <- ifelse(rowSums(gpti[, 3:6]) == 0, NA, gpti[, 5])
gpti[, 6] <- ifelse(rowSums(gpti[, 3:6]) == 0, NA, gpti[, 6])
gpti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

gpti <- gpti[c(1,2,6,10,11,13,14,15,17,18,19),]


# Set up maximum x-axis value for xlim. Add an additional 5%
gpmin.x <- min(gpti[, 3:6], na.rm = TRUE)
gpmin.x <- max(0, gpmin.x, na.rm = TRUE)
gpmax.x <- max(gpti[, 3:6], na.rm = TRUE)
gpmax.x <- gpmax.x * 1.05
gpti$Setup <- seq(gpmin.x, gpmax.x, length.out = nrow(gpti))

# Plot built up in 2 stages: i) using min and max values to set the
#   y-axis range without having to use ylim (so this can be customised
#   by the user), ii) the main dbarplot and label, and iii) the error
graphics::dotchart(gtti$Setup, labels = paste(gtti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Choice differences between genera")
graphics::abline(v = 0, lty = 2, col = "dimgrey")

for (i in 1:nrow(gbti)){
  eval(parse(text = paste("lines(x = c(gbti$Lower.", 0.95 * 100,
                          ".CL[i], gbti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#440154FF', 0.5))", sep = "")))
  graphics::points(gbti$Observed[i], i, pch = if(gbti$Test[i] == "ns") {
    21
  } else {
    if(gbti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#440154FF",0.5), cex = 2)
}
for (i in 1:nrow(geti)){
  eval(parse(text = paste("lines(x = c(geti$Lower.", 0.95 * 100,
                          ".CL[i], geti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#3B528BFF', 0.5))", sep = "")))
  graphics::points(geti$Observed[i], i, pch = if(geti$Test[i] == "ns") {
    21
  } else {
    if(geti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#3B528BFF",0.5), cex = 2)
}
for (i in 1:nrow(gmti)){
  eval(parse(text = paste("lines(x = c(gmti$Lower.", 0.95 * 100,
                          ".CL[i], gmti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#21908CFF', 0.5))", sep = "")))
  graphics::points(gmti$Observed[i], i, pch = if(gmti$Test[i] == "ns") {
    21
  } else {
    if(gmti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#21908CFF",0.5), cex = 2)
}
for (i in 1:nrow(gpti)){
  eval(parse(text = paste("lines(x = c(gpti$Lower.", 0.95 * 100,
                          ".CL[i], gpti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#5DC863FF', 0.5))", sep = "")))
  graphics::points(gpti$Observed[i], i, pch = if(gpti$Test[i] == "ns") {
    21
  } else {
    if(gpti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#5DC863FF",0.5), cex = 2)
}
for (i in 1:nrow(gtti)){
  eval(parse(text = paste("lines(x = c(gtti$Lower.", 0.95 * 100,
                          ".CL[i], gtti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#FDE725FF', 0.5))", sep = "")))
  graphics::points(gtti$Observed[i], i, pch = if(gtti$Test[i] == "ns") {
    21
  } else {
    if(gtti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#FDE725FF",0.5), cex = 2)
}


pdf("Genus_TS_ENNRsigdiff.pdf", width = 6, height = 10) 
graphics::dotchart(gtti$Setup, labels = paste(gtti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Choice differences between genera")
graphics::abline(v = 0, lty = 2, col = "dimgrey")

for (i in 1:nrow(gbti)){
  eval(parse(text = paste("lines(x = c(gbti$Lower.", 0.95 * 100,
                          ".CL[i], gbti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#440154FF', 0.5))", sep = "")))
  graphics::points(gbti$Observed[i], i, pch = if(gbti$Test[i] == "ns") {
    21
  } else {
    if(gbti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#440154FF",0.5), cex = 2)
}
for (i in 1:nrow(geti)){
  eval(parse(text = paste("lines(x = c(geti$Lower.", 0.95 * 100,
                          ".CL[i], geti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#3B528BFF', 0.5))", sep = "")))
  graphics::points(geti$Observed[i], i, pch = if(geti$Test[i] == "ns") {
    21
  } else {
    if(geti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#3B528BFF",0.5), cex = 2)
}
for (i in 1:nrow(gmti)){
  eval(parse(text = paste("lines(x = c(gmti$Lower.", 0.95 * 100,
                          ".CL[i], gmti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#21908CFF', 0.5))", sep = "")))
  graphics::points(gmti$Observed[i], i, pch = if(gmti$Test[i] == "ns") {
    21
  } else {
    if(gmti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#21908CFF",0.5), cex = 2)
}
for (i in 1:nrow(gpti)){
  eval(parse(text = paste("lines(x = c(gpti$Lower.", 0.95 * 100,
                          ".CL[i], gpti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#5DC863FF', 0.5))", sep = "")))
  graphics::points(gpti$Observed[i], i, pch = if(gpti$Test[i] == "ns") {
    21
  } else {
    if(gpti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#5DC863FF",0.5), cex = 2)
}
for (i in 1:nrow(gtti)){
  eval(parse(text = paste("lines(x = c(gtti$Lower.", 0.95 * 100,
                          ".CL[i], gtti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#FDE725FF', 0.5))", sep = "")))
  graphics::points(gtti$Observed[i], i, pch = if(gtti$Test[i] == "ns") {
    21
  } else {
    if(gtti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#FDE725FF",0.5), cex = 2)
}
dev.off()

tiff(file = "Genus_TS_ENNRsigdiff.tiff", width = 12, height = 20, res = 600, units = "cm")
graphics::dotchart(gtti$Setup, labels = paste(gtti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "Choice differences between genera")
graphics::abline(v = 0, lty = 2, col = "dimgrey")

for (i in 1:nrow(gbti)){
  eval(parse(text = paste("lines(x = c(gbti$Lower.", 0.95 * 100,
                          ".CL[i], gbti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#440154FF', 0.5))", sep = "")))
  graphics::points(gbti$Observed[i], i, pch = if(gbti$Test[i] == "ns") {
    21
  } else {
    if(gbti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#440154FF",0.5), cex = 2)
}
for (i in 1:nrow(geti)){
  eval(parse(text = paste("lines(x = c(geti$Lower.", 0.95 * 100,
                          ".CL[i], geti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#3B528BFF', 0.5))", sep = "")))
  graphics::points(geti$Observed[i], i, pch = if(geti$Test[i] == "ns") {
    21
  } else {
    if(geti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#3B528BFF",0.5), cex = 2)
}
for (i in 1:nrow(gmti)){
  eval(parse(text = paste("lines(x = c(gmti$Lower.", 0.95 * 100,
                          ".CL[i], gmti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#21908CFF', 0.5))", sep = "")))
  graphics::points(gmti$Observed[i], i, pch = if(gmti$Test[i] == "ns") {
    21
  } else {
    if(gmti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#21908CFF",0.5), cex = 2)
}
for (i in 1:nrow(gpti)){
  eval(parse(text = paste("lines(x = c(gpti$Lower.", 0.95 * 100,
                          ".CL[i], gpti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#5DC863FF', 0.5))", sep = "")))
  graphics::points(gpti$Observed[i], i, pch = if(gpti$Test[i] == "ns") {
    21
  } else {
    if(gpti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#5DC863FF",0.5), cex = 2)
}
for (i in 1:nrow(gtti)){
  eval(parse(text = paste("lines(x = c(gtti$Lower.", 0.95 * 100,
                          ".CL[i], gtti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#FDE725FF', 0.5))", sep = "")))
  graphics::points(gtti$Observed[i], i, pch = if(gtti$Test[i] == "ns") {
    21
  } else {
    if(gtti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#FDE725FF",0.5), cex = 2)
}
dev.off()

##### Sex differences TS ENNR #####

# Call test_interactions
sfti <- test_interactions(sex.null, signif.level = 0.95)
sfti <- sfti[sfti$Consumer == "Female", ]
sfti[, 3] <- ifelse(rowSums(sfti[, 3:6]) == 0, NA, sfti[, 3])
sfti[, 4] <- ifelse(rowSums(sfti[, 3:6]) == 0, NA, sfti[, 4])
sfti[, 5] <- ifelse(rowSums(sfti[, 3:6]) == 0, NA, sfti[, 5])
sfti[, 6] <- ifelse(rowSums(sfti[, 3:6]) == 0, NA, sfti[, 6])
sfti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

sfti <- sfti[c(8,17),]


# Set up maximum x-axis value for xlim. Add an additional 5%
sfmin.x <- min(sfti[, 3:6], na.rm = TRUE)
sfmin.x <- max(0, sfmin.x, na.rm = TRUE)
sfmax.x <- max(sfti[, 3:6], na.rm = TRUE)
sfmax.x <- sfmax.x * 1.05
sfti$Setup <- seq(sfmin.x, sfmax.x, length.out = nrow(sfti))

# Call test_interactions
smti <- test_interactions(sex.null, signif.level = 0.95)
smti <- smti[smti$Consumer == "Male", ]
smti[, 3] <- ifelse(rowSums(smti[, 3:6]) == 0, NA, smti[, 3])
smti[, 4] <- ifelse(rowSums(smti[, 3:6]) == 0, NA, smti[, 4])
smti[, 5] <- ifelse(rowSums(smti[, 3:6]) == 0, NA, smti[, 5])
smti[, 6] <- ifelse(rowSums(smti[, 3:6]) == 0, NA, smti[, 6])
smti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

smti <- smti[c(8,17),]


# Set up maximum x-axis value for xlim. Add an additional 5%
smmin.x <- min(smti[, 3:6], na.rm = TRUE)
smmin.x <- max(0, smmin.x, na.rm = TRUE)
smmax.x <- max(smti[, 3:6], na.rm = TRUE)
smmax.x <- smmax.x * 1.05
smti$Setup <- seq(smmin.x, smmax.x, length.out = nrow(smti))

# Plot built up in 2 stages: i) using min and max values to set the
#   y-axis range without having to use ylim (so this can be customised
#   by the user), ii) the main dbarplot and label, and iii) the error
graphics::dotchart(sfti$Setup, labels = paste(sfti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "")
graphics::abline(v = 0, lty = 2, col = "dimgrey")

for (i in 1:nrow(sfti)){
  eval(parse(text = paste("lines(x = c(sfti$Lower.", 0.95 * 100,
                          ".CL[i], sfti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#440154FF', 0.5))", sep = "")))
  graphics::points(sfti$Observed[i], i, pch = if(sfti$Test[i] == "ns") {
    21
  } else {
    if(sfti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#440154FF", 0.5), cex = 2)
}


for (i in 1:nrow(smti)){
  eval(parse(text = paste("lines(x = c(smti$Lower.", 0.95 * 100,
                          ".CL[i], smti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#FDE725FF', 0.5))", sep = "")))
  graphics::points(smti$Observed[i], i, pch = if(smti$Test[i] == "ns") {
    21
  } else {
    if(smti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#FDE725FF", 0.5), cex = 2)
}



pdf("Sex_TS_ENNRsigdiff.pdf", width = 6, height = 4) 

graphics::dotchart(sfti$Setup, labels = paste(sfti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "")
graphics::abline(v = 0, lty = 2, col = "dimgrey")

for (i in 1:nrow(sfti)){
  eval(parse(text = paste("lines(x = c(sfti$Lower.", 0.95 * 100,
                          ".CL[i], sfti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#440154FF', 0.5))", sep = "")))
  graphics::points(sfti$Observed[i], i, pch = if(sfti$Test[i] == "ns") {
    21
  } else {
    if(sfti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#440154FF", 0.5), cex = 2)
}


for (i in 1:nrow(smti)){
  eval(parse(text = paste("lines(x = c(smti$Lower.", 0.95 * 100,
                          ".CL[i], smti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#FDE725FF', 0.5))", sep = "")))
  graphics::points(smti$Observed[i], i, pch = if(smti$Test[i] == "ns") {
    21
  } else {
    if(smti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#FDE725FF", 0.5), cex = 2)
}
dev.off()


tiff(file = "Sex_TS_ENNRsigdiff.tiff", width = 12, height = 8, res = 600, units = "cm")
graphics::dotchart(sfti$Setup, labels = paste(sfti$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "")
graphics::abline(v = 0, lty = 2, col = "dimgrey")

for (i in 1:nrow(sfti)){
  eval(parse(text = paste("lines(x = c(sfti$Lower.", 0.95 * 100,
                          ".CL[i], sfti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#440154FF', 0.5))", sep = "")))
  graphics::points(sfti$Observed[i], i, pch = if(sfti$Test[i] == "ns") {
    21
  } else {
    if(sfti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#440154FF", 0.5), cex = 2)
}


for (i in 1:nrow(smti)){
  eval(parse(text = paste("lines(x = c(smti$Lower.", 0.95 * 100,
                          ".CL[i], smti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#FDE725FF', 0.5))", sep = "")))
  graphics::points(smti$Observed[i], i, pch = if(smti$Test[i] == "ns") {
    21
  } else {
    if(smti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#FDE725FF", 0.5), cex = 2)
}
dev.off()

##### Life stage differences TS ENNR #####

# Call test_interactions
lati <- test_interactions(life.null, signif.level = 0.95)
lati <- lati[lati$Consumer == "Adult", ]
lati[, 3] <- ifelse(rowSums(lati[, 3:6]) == 0, NA, lati[, 3])
lati[, 4] <- ifelse(rowSums(lati[, 3:6]) == 0, NA, lati[, 4])
lati[, 5] <- ifelse(rowSums(lati[, 3:6]) == 0, NA, lati[, 5])
lati[, 6] <- ifelse(rowSums(lati[, 3:6]) == 0, NA, lati[, 6])
lati

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

lati <- lati[c(1,5,7,8,10,11,14,16,19),]


# Set up maximum x-axis value for xlim. Add an additional 5%
lamin.x <- min(lati[, 3:6], na.rm = TRUE)
lamin.x <- max(0, lamin.x, na.rm = TRUE)
lamax.x <- max(lati[, 3:6], na.rm = TRUE)
lamax.x <- lamax.x * 1.05
lati$Setup <- seq(lamin.x, lamax.x, length.out = nrow(lati))

# Call test_interactions
ljti <- test_interactions(life.null, signif.level = 0.95)
ljti <- ljti[ljti$Consumer == "Juvenile", ]
ljti[, 3] <- ifelse(rowSums(ljti[, 3:6]) == 0, NA, ljti[, 3])
ljti[, 4] <- ifelse(rowSums(ljti[, 3:6]) == 0, NA, ljti[, 4])
ljti[, 5] <- ifelse(rowSums(ljti[, 3:6]) == 0, NA, ljti[, 5])
ljti[, 6] <- ifelse(rowSums(ljti[, 3:6]) == 0, NA, ljti[, 6])
ljti

# EDIT 'ti' - to just the prey taxa that you want to show on the plot

ljti <- ljti[c(1,5,7,8,10,11,14,16,19),]


# Set up maximum x-axis value for xlim. Add an additional 5%
ljmin.x <- min(ljti[, 3:6], na.rm = TRUE)
ljmin.x <- max(0, ljmin.x, na.rm = TRUE)
ljmax.x <- max(ljti[, 3:6], na.rm = TRUE)
ljmax.x <- ljmax.x * 1.05
ljti$Setup <- seq(ljmin.x, ljmax.x, length.out = nrow(ljti))

# Plot built up in 2 stages: i) using min and max values to set the
#   y-axis range without having to use ylim (so this can be customised
#   by the user), ii) the main dbarplot and label, and iii) the error
graphics::dotchart(lati$Setup, labels = paste(lati$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "")
graphics::abline(v = 0, lty = 2, col = "dimgrey")

for (i in 1:nrow(lati)){
  eval(parse(text = paste("lines(x = c(lati$Lower.", 0.95 * 100,
                          ".CL[i], lati$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#440154FF', 0.5))", sep = "")))
  graphics::points(lati$Observed[i], i, pch = if(lati$Test[i] == "ns") {
    21
  } else {
    if(lati$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#440154FF",0.5), cex = 2)
}
for (i in 1:nrow(ljti)){
  eval(parse(text = paste("lines(x = c(ljti$Lower.", 0.95 * 100,
                          ".CL[i], ljti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#FDE725FF', 0.5))", sep = "")))
  graphics::points(ljti$Observed[i], i, pch = if(ljti$Test[i] == "ns") {
    21
  } else {
    if(ljti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#FDE725FF",0.5), cex = 2)
}

pdf("Life_TS_ENNRsigdiff.pdf", width = 6, height = 8) 
graphics::dotchart(lati$Setup, labels = paste(lati$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "")
graphics::abline(v = 0, lty = 2, col = "dimgrey")

for (i in 1:nrow(lati)){
  eval(parse(text = paste("lines(x = c(lati$Lower.", 0.95 * 100,
                          ".CL[i], lati$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#440154FF', 0.5))", sep = "")))
  graphics::points(lati$Observed[i], i, pch = if(lati$Test[i] == "ns") {
    21
  } else {
    if(lati$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#440154FF",0.5), cex = 2)
}
for (i in 1:nrow(ljti)){
  eval(parse(text = paste("lines(x = c(ljti$Lower.", 0.95 * 100,
                          ".CL[i], ljti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#FDE725FF', 0.5))", sep = "")))
  graphics::points(ljti$Observed[i], i, pch = if(ljti$Test[i] == "ns") {
    21
  } else {
    if(ljti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#FDE725FF",0.5), cex = 2)
}
dev.off()

tiff(file = "Life_TS_ENNRsigdiff.tiff", width = 12, height = 16, res = 600, units = "cm")
graphics::dotchart(lati$Setup, labels = paste(lati$Resource, " ", sep = ""),
                   col = 1, pt.cex = 0, cex = 1.5, main = "")
graphics::abline(v = 0, lty = 2, col = "dimgrey")

for (i in 1:nrow(lati)){
  eval(parse(text = paste("lines(x = c(lati$Lower.", 0.95 * 100,
                          ".CL[i], lati$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#440154FF', 0.5))", sep = "")))
  graphics::points(lati$Observed[i], i, pch = if(lati$Test[i] == "ns") {
    21
  } else {
    if(lati$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#440154FF",0.5), cex = 2)
}
for (i in 1:nrow(ljti)){
  eval(parse(text = paste("lines(x = c(ljti$Lower.", 0.95 * 100,
                          ".CL[i], ljti$Upper.", 0.95 * 100,
                          ".CL[i]), y = c(i, i), col = alpha('#FDE725FF', 0.5))", sep = "")))
  graphics::points(ljti$Observed[i], i, pch = if(ljti$Test[i] == "ns") {
    21
  } else {
    if(ljti$Test[i] == "Weaker") {
      25
    } else{
      24
    }
  }, col = "black",
  bg = alpha("#FDE725FF",0.5), cex = 2)
}
dev.off()

