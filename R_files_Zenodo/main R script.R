# The evolution of life history theory
# R script for analyses reported in paper

###Head#####
# Libraries used
library(bibliometrix)
library(ggplot2)
library(dplyr)
library(DescTools)
library(psych)
library(cowplot)

###Descriptive analyses of overall dataset including figure 1#####
# Set working directory now
refs_import <- readFiles("pre2010part1.bib", "pre2010part2.bib", "post2010part1.bib", "post2010part2.bib")
d<-convert2df(refs_import, format="bibtex", dbsource="isi")
r=biblioAnalysis(d)
print(r)
# Make figure 1
figdat=xtabs(~d$PY)
figdat=as.data.frame(figdat)
figdat$cumfreq=cumsum(figdat$Freq)
figdat=subset(figdat, d.PY!="2018")
fig1=ggplot(figdat, aes(x=as.numeric(as.character(d.PY)), y=Freq)) + 
  geom_point() + 
  geom_smooth(se=F) + 
  theme_bw() + 
  xlab("Year") + 
  ylab("Publications")  +
  guides(colour=guide_legend(title=NULL)) + 
  theme(legend.position=c(0.25, 0.75))
fig1
png("figure1.png", res=300, units="in", width=4, height=3)
fig1
dev.off()

###Connection indices for the clusters (figure 3)####
# Read in map and cluster files created by VOS viewer
# Pre-2010
network <- read.csv(file = "pre2010network.csv")
map<-read.csv("pre2010map.csv")
# Make matrix of connections
clusters=select(map, id, cluster)
d=merge(network, clusters, by.x="from.id", by.y="id")
d=merge(d, clusters, by.x="to.id", by.y="id")
parta=xtabs(~d$cluster.x + d$cluster.y) 
partb=t(parta)
link.matrix=parta + partb
xtabs(~map$cluster)
sum(xtabs(~map$cluster))
# Numbers of papers per cluster
n=c(120, 115, 109, 93, 63)
# Expected number of links
expected=n%*%t(n)
diag(expected)=diag(expected)-n
expected
f=link.matrix/expected
# Make data frame of the connection indices
x=c(rep("A1", times=5), rep("A2", times=5), rep("A3", times=5), rep("A4", times=5), rep("A5", times=5))
y=rep(c("A1", "A2", "A3", "A4", "A5"), times=5)
z=as.vector(f)
l=data.frame(x, y, z)
# Figure 3a
figurea=ggplot(l, aes(x=x, y=z)) + 
  geom_bar(stat="identity", fill="darkred", colour="black") + 
  theme_bw() + 
  facet_grid(~y) + 
  xlab("Cluster") + 
  ylab("Connection index")  
  coord_cartesian(ylim=c(0, 1))
figurea
# Overall connection probability
mean(l$z)
# Calculate within/between ratio
counter=0
ratios=NULL
for (i in levels(l$x)) {
  print(i)
  counter=counter+1
  current.subset=subset(l, x==i)
  current.within=mean(current.subset$z[current.subset$y==i])
  current.between=mean(current.subset$z[current.subset$y!=i])
  current.ratio=current.within/current.between
  ratios[counter]=current.ratio
}
ratios
# Gini coefficients
dplyr::summarise(group_by(l, x), g=Gini(z))

#Post 2010
network2 <- read.csv(file = "post2010network.csv")
map2<-read.csv("post2010map.csv")
xtabs(~map2$cluster)
sum(xtabs(~map2$cluster))
clusters2=select(map2, id, cluster)
d2=merge(network2, clusters2, by.x="cluster.from", by.y="id")
d2=merge(d2, clusters2, by.x="cluster.to", by.y="id")
parta2=xtabs(~d2$cluster.x + d2$cluster.y) 
partb2=t(parta2)
link.matrix2=parta2 + partb2
link.matrix2
xtabs(~map2$cluster)
# Numbers of papers per cluster
n2=c(184, 155, 90, 40, 29)
# Expected number of links
expected2=n2%*%t(n2)
diag(expected2)=diag(expected2)-n2
f2=link.matrix2/expected2
# Make data frame of connection indices
x=c(rep("B1", times=5), rep("B2", times=5), rep("B3", times=5), rep("B4", times=5), rep("B5", times=5))
y=rep(c("B1", "B2", "B3", "B4", "B5"), times=5)
z=as.vector(f2)
l2=data.frame(x, y, z)
# Mean connection probability
mean(l2$z)
# Within/between ratio
counter=0
ratios=NULL
for (i in levels(l2$x)) {
  print(i)
  counter=counter+1
  current.subset=subset(l2, x==i)
  current.within=mean(current.subset$z[current.subset$y==i])
  current.between=mean(current.subset$z[current.subset$y!=i])
  current.ratio=current.within/current.between
  ratios[counter]=current.ratio
}
ratios
# Gini coeffient
dplyr::summarise(group_by(l2, x), g=Gini(z))
# Figure 3b
figureb=ggplot(l2, aes(x=x, y=z)) + 
  geom_bar(stat="identity", fill="darkred", colour="black") + 
  theme_bw() + 
  facet_grid(~y) + 
  xlab("Cluster") + 
  ylab("Connection index") + 
  coord_cartesian(ylim=c(0, 1))
figureb
mean(l2$z)
dplyr::summarise(group_by(l2, x), g=Gini(z))

# Save figure 3
png("figure3.png", res=300, width=6, height=6, units="in")
plot_grid(figurea, figureb, nrow=2, labels=c("A", "B"))
dev.off()

###Working out most cited papers per cluster (for table 1)####
refs_import <- readFiles("pre2010part1.bib", "pre2010part2.bib", "post2010part1.bib", "post2010part2.bib")
d<-convert2df(refs_import, format="bibtex", dbsource="isi")
map1<-read.csv("pre2010map.csv")
map2<-read.csv("post2010map.csv")
map1$cluster=paste("A", map1$cluster)
map2$cluster=paste("B", map2$cluster)
map1$doi=substring(map1$url, first=17, last=1000)
map1$doi=toupper(map1$doi)
map2$doi=substring(map2$url, first=17, last=1000)
map2$doi=toupper(map2$doi)
colnames(map1)=colnames(map2)
maps=rbind(map1, map2)
f=merge(d, maps, by.x="DI", by.y="doi", all.x=TRUE, all.y=F)
xtabs(~f$cluster)
sum(xtabs(~f$cluster))
xtabs(~map1$cluster)
xtabs(~map2$cluster)
# Most cited sources by cluster
citationsA1<-citations(subset(f, cluster=="A 1"), field="article", sep=".  ")
citationsA1$Cited[1:20]
citationsA2<-citations(subset(f, cluster=="A 2"), field="article", sep=".  ")
citationsA2$Cited[1:20]
citationsA3<-citations(subset(f, cluster=="A 3"), field="article", sep=".  ")
citationsA3$Cited[1:20]
citationsA4<-citations(subset(f, cluster=="A 4"), field="article", sep=".  ")
citationsA4$Cited[1:20]
citationsA5<-citations(subset(f, cluster=="A 5"), field="article", sep=".  ")
citationsA5$Cited[1:20]

citationsB1<-citations(subset(f, cluster=="B 1"), field="article", sep=".  ")
citationsB1$Cited[1:20]
citationsB2<-citations(subset(f, cluster=="B 2"), field="article", sep=".  ")
citationsB2$Cited[1:20]
citationsB3<-citations(subset(f, cluster=="B 3"), field="article", sep=".  ")
citationsB3$Cited[1:20]
citationsB4<-citations(subset(f, cluster=="B 4"), field="article", sep=".  ")
citationsB4$Cited[1:20]
citationsB5<-citations(subset(f, cluster=="B 5"), field="article", sep=".  ")
# Data frame of most cited references
cited.refs=rbind(as.data.frame(citationsA1$Cited[1:10]), 
                 as.data.frame(citationsA2$Cited[1:10]), 
                 as.data.frame(citationsA3$Cited[1:10]), 
                 as.data.frame(citationsA4$Cited[1:10]), 
                 as.data.frame(citationsA5$Cited[1:10]), 
                 as.data.frame(citationsB1$Cited[1:10]), 
                 as.data.frame(citationsB2$Cited[1:10]), 
                 as.data.frame(citationsB3$Cited[1:10]), 
                 as.data.frame(citationsB4$Cited[1:10]), 
                 as.data.frame(citationsB5$Cited[1:10]))

cited.refs$cluster=c(rep("A1", times=10),
                     rep("A2", times=10),
                     rep("A3", times=10),
                     rep("A4", times=10),
                     rep("A5", times=10),
                     rep("B1", times=10),
                     rep("B2", times=10),
                     rep("B3", times=10),
                     rep("B4", times=10),
                     rep("B5", times=10)
)
# You can write this to a csv
write.csv(cited.refs, file="cited.refs.csv")

###Miscellaneous - relating life history theory to 'pace of life'#####
pol_import <- readFiles("savedpolrecs.bib")
p = convert2df(pol_import, format="bibtex", dbsource="isi")
citationspol<-citations(p, field="article", sep=".  ")
citationspol$Cited[1:20]

dfB2=as.data.frame(citationsB2$Cited)
grep("REALE", dfB2$CR)
dfB2[182, ]
grep("RICKLEFS", dfB2$CR)
dfB2[1193, ]
grep("BIRO", dfB2$CR)
dfB2[287, ]

dfB5=as.data.frame(citationsB5$Cited)
grep("REALE", dfB5$CR)
grep("RICKLEFS", dfB5$CR)
grep("BIRO", dfB5$CR)

