## R code to produce Figure 4, and Table 2, in: 
## Johnston EC, Wyatt ASJ, Leichter JJ, Burgess SC. Niche differences in co-occuring cryptic coral species (Pocillopora spp). Coral Reefs.
## Code written by Scott Burgess. March 2021. Send comments or corrections to sburgess@bio.fsu.edu
## R version 3.6.3 (2020-02-29)

rm(list=ls()) #clear all variables

library(vegan)
library(tidyr)
library(dplyr)
library(data.table)


#setwd("")
# Import relative abundance data
dat_long <- read.csv("Figure 3 and 4 Data.csv")

# Import data on temperature and light
PO_data <- read.csv("Figure 1 and 4 Physical Data.csv")

# Convert count data to proportion
names(dat_long)[names(dat_long) == "Depth.m"] <- "Depth"
dat_long <- table(Species.haplotype = dat_long$Species.haplotype, Site = dat_long$Site, Depth = dat_long$Depth)
dat_long <- as.data.frame(dat_long)
dat_long$sum <- ave(dat_long$Freq, dat_long$Depth, dat_long$Site, FUN=sum)
dat_long <- dat_long %>% group_by(.dots=c("Species.haplotype","Depth","Site")) %>% mutate(prop=Freq/sum)


# Prepare data
dat_wide <- spread(dat_long,Species.haplotype,prop)
dat_wide$Depth <- as.factor(dat_wide$Depth)
dat_wide$Site <- as.factor(dat_wide$Site)
dat_wide <- dat_wide[,-c(3,4)]

dat_wide <- dat_wide %>%
  mutate_all(~replace(., .=='NA', NA)) %>%
  summarize_all(~paste(unique(na.omit(.)), collapse = ','))

prop_mat <- dat_wide[,-c(1,2)]
prop_mat <- sapply( prop_mat, as.numeric )
rownames(prop_mat) <- interaction(dat_wide$Depth, dat_wide$Site)


PO_data <- PO_data[order(PO_data$Depth),]

### Distance-based Redundancy Analysis (db-RDA).
## See, for example:
## https://stackoverflow.com/questions/51715281/dbrda-in-r-how-to-with-abundance-data-and-missing-values-for-environmental-data#51767138
## https://archetypalecology.wordpress.com/2018/02/21/distance-based-redundancy-analysis-db-rda-in-r/


####### Model 1 - Constrained by Depth + Site
m1 <- dbrda(prop_mat ~ Depth + Site, distance="bray",sqrt.dist=F,data= dat_wide)
summary(m1) # The first dbRDA axis explained 79.72%, the second axis explained 6.67%


# Perform test of each factor using marginal effect
# which is the  unique effect of each variable conditional on the presence of the other variable in the model
anova(m1,by="margin",permutations=999999) # Table 2 in paper. Note that the p-values will change slightly each time this line is run, more so when permutations are reduced (quicker, but less accurate).


####### Model 2 - Constrained by Environmental variables
m2 <- dbrda(prop_mat ~ 
					Temp.MaxDailyMean.JunNov
					+Temp.MaxDailyMean.DecMay
					+ Temp.MinDailyMean.DecMay
					+ Temp.MeanDailyVariance.DecMay
					+ Temp.MaxDailyVariance.DecMay
					+ Light.Mean.Em2d.DecFeb
					+ Light.Min.Em2d.DecFeb,
					PO_data,
					distance="bray",sqrt.dist=F)

sppscores(m2) <- prop_mat # add species scores
summary(m2)

# What is the unique effect of each variable conditional on the presence of the other variable in the model?
anova(m2,by="margin",permutations=99999) # No one variable is significant, but...
anova(m2,permutations=99999) # ...but the overall model is significant.


############ Prepare plot for paper 
# Note that the analysis uses ALL haplotypes, but here we're only plotting species scores for the common haplotypes

# Extract Site scores (weighted sums of species scores)
xy <- summary(m2)$sites[,1:2]

# Extract species scores
haps <- summary(m2)$species[,1:2]

# Get focal haplotypes and names for plotting
focal.haplotypes <- c(
	"Haplotype 1a_Pm",
	"Haplotype 1a_Pe",
	"Haplotype 8a",
	"Haplotype 11",
	"Haplotype 3b",
	"Haplotype 10")
haplotype.index <- match(focal.haplotypes,rownames(haps))
haps <- haps[haplotype.index,]
hap.names <- c(
	"P. meandrina",
	"P. eydouxi",
	"Hap 8a",
	"Hap 11",
	"Hap 3b",
	"Hap 10")


# For plotting
foo <- data.frame(cbind(xy, Depth = as.character(dat_wide$Depth), Site=as.character(dat_wide$Site)))
foo5 <- foo[foo$Depth=="5",1:2]
foo10 <- foo[foo$Depth=="10",1:2]
foo20 <- foo[foo$Depth=="21",1:2]
fooS1 <- foo[foo$Site=="1",1:2]
fooS2 <- foo[foo$Site=="2",1:2]
fooS4 <- foo[foo$Site=="4",1:2]
fooS5 <- foo[foo$Site=="5",1:2]
polygon.color <- "grey"

SC <- summary(m2)$biplot[,1:2]
SCnames <- c(
	"Max Daily\\nTemp (Jun-Nov)",
	"Max Daily\\nTemp (Dec-May)",
	"Min Daily Temp",
	"Temp Variance (Mean)",
	"Temp Variance (Max)",
	"Mean\\nLight",
	"Min Light")


### Figure 4

dev.new(width=9,height=5)
par(mfrow=c(1,2),mar=c(1,2,2,1),oma=c(3,3,0,0))

# Depth
plot(xy,xlim=c(-1.6,1.6),ylim=c(-1,1.3),xaxt="n",yaxt="n",type="n",cex.lab=1.2)
legend('topleft',legend="a)",bty="n",cex=1.2,adj=0.3)
abline(v=0,h=0,lty=2)
# Add polygons by depth
ch <- chull(foo5)
coords <- as.matrix(foo5[c(ch,ch[1]),])
polygon(coords,col=adjustcolor(polygon.color,alpha.f=0.4),border=NA)
text(x=-0.8,y=0.55,labels="5m",col=polygon.color,cex=1.5)
ch <- chull(foo10)
coords <- as.matrix(foo10[c(ch,ch[1]),])
polygon(coords[c(2,1,3,4,2),],col=adjustcolor(polygon.color,alpha.f=0.4),border=NA)
text(x=0,y=0.55,labels="10m",col=polygon.color,cex=1.5)
ch <- chull(foo20)
coords <- as.matrix(foo20[c(ch,ch[1]),])
polygon(coords,col=adjustcolor(polygon.color,alpha.f=0.4),border=NA)
text(x=1.25,y=0.55,labels="20m",col=polygon.color,cex=1.5)

# Add location of each site.depth, but label by site only
points(xy,pch=19,col=adjustcolor("white",alpha.f=0.6),cex=2.5)
text(xy, as.character(dat_wide$Site),cex=1.2)

# Species names / scores
text(haps, hap.names,cex=1,col="blue")
axis(side=1,at=seq(-1,1,0.5),cex.axis=1.2)
axis(side=2,at=seq(-1,1,0.5),las=2,cex.axis=1.2)


plot(xy,xlim=c(-1.6,1.6),ylim=c(-1,1.3),xaxt="n",yaxt="n",type="n",cex.lab=1.2)
legend('topleft',legend="b)",bty="n",cex=1.2,adj=0.3)
abline(v=0,h=0,lty=2)
# Add polygons by depth
ch <- chull(foo5)
coords <- as.matrix(foo5[c(ch,ch[1]),])
polygon(coords,col=adjustcolor(polygon.color,alpha.f=0.4),border=NA)
text(x=-0.8,y=0.45,labels="5m",col=polygon.color,cex=1.5)
ch <- chull(foo10)
coords <- as.matrix(foo10[c(ch,ch[1]),])
polygon(coords[c(2,1,3,4,2),],col=adjustcolor(polygon.color,alpha.f=0.4),border=NA)
text(x=0,y=0.45,labels="10m",col=polygon.color,cex=1.5)
ch <- chull(foo20)
coords <- as.matrix(foo20[c(ch,ch[1]),])
polygon(coords,col=adjustcolor(polygon.color,alpha.f=0.4),border=NA)
text(x=1.25,y=0.45,labels="20m",col=polygon.color,cex=1.5)

# Add location of each site.depth, but label by site only
# Species names / scores
axis(side=1,at=seq(-1,1,0.5),cex.axis=1.2)
axis(side=2,at=seq(-1,1,0.5),las=2,cex.axis=1.2)

arrows(rep(0,length(SC[,1])),rep(0,length(SC[,1])),SC[,1],SC[,2],length=0.1,col="blue")
xoffsets <- c(-0.28,0,0,-0.05,-0.1,-0.2,-0.4)
yoffsets <- c(-0.15,-0.12,-0.1,0.1,0.1,0.1,0)
text(SC[,1]+xoffsets,SC[,2]+yoffsets,SCnames,col="blue")

# Add Proportion Explained for dbRDA1 and dbRDA2
summary(m2) # The first dbRDA axis explained 81.53%, the second axis explained 5.62%
mtext("dbRDA1 (82%)",side=1,line=1.5,outer=T,cex=1.2)
mtext("dbRDA2 (6%)",side=2,line=1,outer=T,cex=1.2)



