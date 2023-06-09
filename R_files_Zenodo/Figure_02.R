rm(list=ls())

########
# Xing et al., Figure 02
########


##goal: Allometric scaling; plot of eye socket length against skull length in phylogenetic context

#files:
#"skull.orbit.avg.csv"
#"Stage2_MayrParSho_Ericson.nex"
#"Stage2_MayrParSho_Hackett.nex"

##### 1 Preliminaries

#please set working directory
setwd(path/to/your/working/directory)

#call libraries
require(ape)
require(geiger)
require(nlme)
require(phytools)
require(plyr)
require(ggplot2)


#load tree and data
phy.all.1 <- read.nexus("Stage2_MayrParSho_Ericson.nex") #reading in tree samples
phy1 <- multi2di(consensus(phy.all.1, p=0.5)) #majority rule consensus tree; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5837340/
phy1 <- consensus.edges(phy.all.1, method="mean.edge", consensus.tree=phy1, if.absent="zero") #this option looks better than "ignore"
write.tree(phy1, "phy1_Ericson_consensus.tre")
#phy1 <- read.tree("phy1_Ericson_consensus.tre") #that's for future runs (comment out the 4 preceding lines!)

phy.all.2 <- read.nexus("Stage2_MayrParSho_Hackett.nex") #reading in tree samples
phy2 <- multi2di(consensus(phy.all.2, p=0.5))
phy2 <- consensus.edges(phy.all.2, method="mean.edge", consensus.tree=phy2, if.absent="zero")
write.tree(phy2, "phy2_Hackett_consensus.tre")
#phy2 <- read.tree("phy2_Hackett_consensus.tre") #that's for future runs (comment out the 4 preceding lines!)

sockets <- read.csv("skull.orbit.avg.csv", header=T)
rownames(sockets) <- sockets$taxon #assigning rownames (used for matching purposes)
head(sockets)



##### with Ericson consensus tree, i.e. phy1 from above

#match data with tree
data <- sockets[phy1$tip.label,]
compare.data <- treedata(phy1, data, sort=T, warnings=T) # the treedata function compares data with tree and prunes the tree automatically
phy <- compare.data$phy
head(compare.data$data)

ordered.sockets <- data.frame(taxon=as.character(compare.data$data[,1]),
                              group1=as.character(compare.data$data[,5]), 
                              group2=as.character(compare.data$data[,6]), 
                              sk=as.numeric(compare.data$data[,3]),
                              ol=as.numeric(compare.data$data[,4]))
rownames(ordered.sockets) <- ordered.sockets$taxon
ordered.sockets <- ordered.sockets[phy$tip.label,]
head(ordered.sockets)
length(ordered.sockets$taxon)

#sort data
variables <- ordered.sockets[,4:5]
sk <- variables[,1]
names(sk) <- phy$tip.label
orb <- variables[,2]
names(orb) <- phy$tip.label

#round values
sk <- signif(log10(sk), 3)
orb <- signif(log10(orb),3)

#putting it all in a dataframe
M <- data.frame(taxon=phy$tip.label, group1=ordered.sockets$group1, group2=ordered.sockets$group2, sk, orb)
head(M)

#PGLS
fit1 <- gls(orb~sk, data=M, correlation=corPagel(0.5, phy, fixed=FALSE), method="ML")
fit1

#generate prediction
N <- cbind(M, Y=predict(fit1, M))
head(N)

#Not possible to do exact prediction belts for GLS :/

p <- ggplot(N, aes(x=sk, y=orb)) + xlim(0.7, 2.2) + ylim(0.6, 1.8) +
  geom_point(aes(shape=factor(group1)), alpha=1/4, size=4) +
  geom_line(aes(sk, Y)) +
  labs(x = "Log10 Skull Length, mm", y = "Log10 Orbit Length, mm") 
p

weenie <- data.frame(sk = log10(7.1), ol = log10(5.1))
R <- cbind(weenie, Y=predict(fit1, weenie))


pdf("orbit_size_PGLS_Ericson_cons_tree.pdf")
p + geom_point(x=R[[1]], y=R[[2]], colour="black", shape=15, size=5)+ scale_color_grey() + scale_fill_grey() + theme_classic(base_size = 22)+ coord_fixed() +theme(legend.position="none")
dev.off()



##### with Hackett consensus tree, i.e. phy2 from above

#match data with tree
data <- sockets[phy2$tip.label,]
compare.data <- treedata(phy2, data, sort=T, warnings=T) # the treedata function compares data with tree and prunes the tree automatically
phy <- compare.data$phy

ordered.sockets <- data.frame(taxon=as.character(compare.data$data[,1]),
                              group1=as.character(compare.data$data[,5]), 
                              group2=as.character(compare.data$data[,6]), 
                              sk=as.numeric(compare.data$data[,3]),
                              ol=as.numeric(compare.data$data[,4]))
rownames(ordered.sockets) <- ordered.sockets$taxon
ordered.sockets <- ordered.sockets[phy$tip.label,]
head(ordered.sockets)

#sort data
variables <- ordered.sockets[,4:5]
sk <- variables[,1]
names(sk) <- phy$tip.label
orb <- variables[,2]
names(orb) <- phy$tip.label

#round values
sk <- signif(log10(sk), 3)
orb <- signif(log10(orb),3)

#putting it all in a dataframe
M <- data.frame(taxon=phy$tip.label, group1=ordered.sockets$group1, group2=ordered.sockets$group2, sk, orb)
head(M)

#PGLS
fit2 <- gls(orb~sk, data=M, correlation=corPagel(0.5, phy, fixed=FALSE), method="ML")
fit2

#generate prediction
N <- cbind(M, Y=predict(fit2, M))
N 

head(N)

#Not possible to calculate exact prediction belts for GLS to the best of my knowledge :/

#plot
p <- ggplot(N, aes(x=sk, y=orb)) + xlim(0.7, 2.2) + ylim(0.6, 1.8) +
  geom_point(aes(shape=factor(group1)), alpha=1/4, size=4) +
  geom_line(aes(sk, Y)) +
  labs(x = "Log10 Skull Length, mm", y = "Log10 Orbit Length, mm") 
p

weenie <- data.frame(sk = log10(7.1), ol = log10(4.92))
R <- cbind(weenie, Y=predict(fit1, weenie))

pdf("orbit_size_PGLS_Hackett_cons_tree.pdf")
p + geom_point(x=R[[1]], y=R[[2]], colour="black", shape=15, size=5)+ scale_color_grey() + scale_fill_grey() + theme_classic(base_size = 18)+ coord_fixed() + theme(legend.position="none")
dev.off()


#adding taxon names
p + geom_point(x=R[[1]], y=R[[2]], colour="black", shape=15, size=5)+ 
  scale_color_grey() + scale_fill_grey() +
  theme_classic(base_size = 18)+ coord_fixed() + 
  geom_text(data = N, aes(x = sk, y = orb, label = taxon), 
            size = 2, vjust = 0, hjust = -0.1) +
  theme(legend.position="none")


#adding in some other birds to get a better idea how the amber bird sizes up
weenie <- data.frame(sk = log10(7.1), ol = log10(4.92))
gallus_gallus <- data.frame(sk = 1.6, ol = 1.35)
rhea_americana <- data.frame(sk = 1.9, ol = 1.68)
struthio_camelus <- data.frame(sk = 1.69, ol = 1.5)
mellisuga <- data.frame(sk = 0.946, ol = 0.743)
poi <- rbind(weenie, gallus_gallus, rhea_americana, struthio_camelus, mellisuga)
rownames(poi) <- c("Weenie", "Chicken", "Rhea", "Ostrich", "Mellisuga")
R2 <- cbind(poi, Y=predict(fit1, poi))

pdf("Figure_02.pdf")
p + geom_point(x=weenie[[1]], y=weenie[[2]], colour="black", shape=15, size=5) +
  geom_point(x=gallus_gallus[[1]], y=gallus_gallus[[2]], colour="red", shape=15, size=5) +
  geom_point(x=rhea_americana[[1]], y=rhea_americana[[2]], colour="blue", shape=15, size=5) +
  geom_point(x=struthio_camelus[[1]], y=struthio_camelus[[2]], colour="green", shape=15, size=5) +
  geom_point(x=mellisuga[[1]], y=mellisuga[[2]], colour="purple", shape=15, size=5) +
  scale_color_grey() + scale_fill_grey() +
  theme_classic(base_size = 18)+ coord_fixed() + 
  theme(legend.position="none")
dev.off()




















