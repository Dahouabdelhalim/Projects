
### R code to conduct analyses in Hahn et al. "Population variation, environmental gradients, and 
      ### the evolutionary ecology of plant defense against herbivory" 
      ### published in the American Naturalist.
  ## Contact Phil Hahn (phil.hahn@mso.umt.edu, www.plant-herbivore-interactions.net) for questions or additional details.


### load necassary packages, install if necassary, set wd if necassary ##
library(lattice)
library(reshape2)
library(Hmisc)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)


##########################################################
#### Load and process Site info and data from natural population  ###########

site <- read.csv("Hahn_AmNat18_SiteInfo.csv", header=T)

popp <-read.csv("Hahn_AmNat18_NaturalTraits.csv", header=T)
popp
pop <- merge(popp,site,by="Pop")
pop$latex <- log(pop$latexmgf)  ### log transform latex
pop$No_ram <- log(pop$No_ram+1) ### log transform # ramets
pop$cards1 <- sqrt(pop$cards)   ### sqrt transform cardenolide concentration

### change column names
colnames(pop)[colnames(pop)=="Height"] <- "Height"
colnames(pop)[colnames(pop)=="slaf"] <- "SLA"
colnames(pop)[colnames(pop)=="cards1"] <- "Cards"
colnames(pop)[colnames(pop)=="d13Cf"] <- "dC13"
colnames(pop)[colnames(pop)=="Nf"] <- "LeafN"
colnames(pop)[colnames(pop)=="latex"] <- "Latex"
colnames(pop)[colnames(pop)=="No_fruit"] <- "Pods"



#############################################################################
##############################################################################
########### FULL ANALYSIS OF DATA FROM NATURAL POPULATIONS ##################

### herbivore (beetle) abundance and damage analysis 
b1 <- lm(log(beetles+1)~clim, data=pop)
summary(b1)
anova(b1)
b2 <- lm(Damage~log(beetles+1), data=pop)
summary(b2)
anova(b2)
b3 <- lm(Damage~clim, data=pop)
summary(b3)
anova(b3)

head(pop)
br <- colorRampPalette(c("darkblue","blue","cyan","green","yellow","orange","red"))
pop$Col <- br(17)[as.numeric(cut(pop$clim,breaks=17))]

pop1 <- pop[c("Pop","clim","soils","beetles","Damage","Height","No_ram","Pods","SLA","dC13","LeafN",
              "Trichome","Cards","Latex","Col")]
popgraph <- pop1  ### save un-scaled trait values to make figure later
pop1[6:14] <- scale(pop1[6:14])  ### center and scale pop-average traits in natural populations
pop1


#################################################################################################
##### NATURAL POPS CLINAL ANALYSIS #############################################################

##### reshape scaled trait for analysis 
pop2 <- melt(pop1, id.vars=c("Pop","clim","soils","beetles","Damage","Col"), 
             measure.vars=c("Height","No_ram","Pods","SLA","dC13","LeafN",
                            "Trichome","Cards","Latex"))
colnames(pop2)[7] <- "Trait"
head(pop2)

hist(log(pop2$Damage+1))

#### make plot of each trait ploted against bioclimPC to examine
xyplot(value~clim|Trait, data=pop2, 
       ylab=list("centered trait value",cex=1.5),xlab=list("Bioclim PC",cex=1.5),
       main="Clinal variation of traits",
       type = c("r", "p"), pch=16, col=pop2$Col, cex=1.25, lty=1, lwd=2,
       par.settings = list(strip.background=list(col="lightgrey"),cex=1.5))

#### contruct and examine model with climate, damage, soils predicting values of each trait
## results are in Table C2; can also swap "sqrt(Damage+1)" for beetles
ff <- lmer(value~clim*Trait+log(beetles+1)*Trait+soils*Trait+(1|Pop), data=pop2)
summary(ff)
anova(ff) ### Table is in Appendix: Table C2


###### construct reduced model -- eq. 1 in manuscript
f1 <- lmer(value~clim*Trait+(1|Pop), data=pop2)
summary(f1)
anova(f1)
f1a <- emtrends(f1, ~Trait, var="clim") ### examine slopes for each trait with bioclimPC
f1b <- as.data.frame(summary(f1a))  ## save slopes in dataframe to make figure later (fig. 3A)


##########################################################################################
#### GARDEN CLINAL ANALYSIS ###################################
#### load and process trait data from garden ###

s1 <- read.csv("Hahn_AmNat18_GardenTraits.csv", header=T)
head(s1)
summary(s1)

s1$cards1 <- sqrt(s1$cards)  ## square-root transform cardenolide concentrations
s1$slat <- sqrt(s1$sla)      ## square-root transform specific leaf area

hist(s1$cards)
hist(s1$cards1)
summary(s1$cards)
summary(s1$cards1)

colnames(s1)[colnames(s1)=="Hgt"] <- "Height"
colnames(s1)[colnames(s1)=="Lstems"] <- "No_ram"
colnames(s1)[colnames(s1)=="PodMass"] <- "Pods"
colnames(s1)[colnames(s1)=="slat"] <- "SLA"
colnames(s1)[colnames(s1)=="cards1"] <- "Cards"
colnames(s1)[colnames(s1)=="d13C"] <- "dC13"
colnames(s1)[colnames(s1)=="N"] <- "LeafN"
colnames(s1)[colnames(s1)=="log_latex_mg"] <- "Latex"
colnames(s1)[colnames(s1)=="Trics"] <- "Trichome"

#### use aggregate to make pop-averages 
g2m <- aggregate(cbind(clim,soils,beetles,Height,No_ram,Pods,SLA,dC13,LeafN,Trichome,Cards,Latex)~Pop, data=s1, mean)
g2 <- s1[c("Pop","Fam","Plot",'clim','soils','beetles','Height','No_ram','Pods','SLA','dC13','LeafN','Trichome','Cards','Latex')]
gargraph <- merge(g2m,pop1[c("Pop","Col")], by="Pop")  ## save un-scaled pop-average traits for figure later

g2m[5:13] <- scale(g2m[5:13]) ### center and scale pop-average traits
g2m
g2[7:15] <- scale(g2[7:15])  ### center and scale indivdual-plant traits
head(g2)

g2ma <- merge(g2m,pop1[c("Pop","Col")], by="Pop")  ### merge pop-average traits for natural pops and garden 

#### reshape garden trait data (non-aggregated) wide to long for analysis
g3 <- melt(g2, id.vars=c("Pop","Fam","Plot","clim","soils","beetles"),
            measure.vars=c("Height","No_ram","Pods","SLA","dC13","LeafN",
                           "Trichome","Cards","Latex"))

g3m <- melt(g2ma, id.vars=c("Pop","clim","soils","beetles","Col"),
            measure.vars=c("Height","No_ram","Pods","SLA","dC13","LeafN",
                           "Trichome","Cards","Latex"))

colnames(g3)[7] <- "Trait"
head(g3)
colnames(g3m)[6] <- "Trait"
head(g3m)

## examine clinal patterns for garden data
xyplot(value~clim|Trait, data=g3, 
       ylab=list("centered trait value",cex=1.5),xlab=list("Climate PC",cex=1.5),
       main="Clinal variation of traits",
       type = c("r", "p"), pch=1, cex=1, lty=1, lwd=2,
       par.settings = list(strip.background=list(col="lightgrey"),cex=1))

#### construct model to analyze garden data, including climate, beetles, and soils 
  ## results are in Appendix: Table C3
garf <- lmer(value~clim*Trait+log(beetles+1)*Trait+soils*Trait+(1|Pop)+(1|Pop:Trait)+(1|Pop:Trait:Fam), data=g3)
summary(garf)
anova(garf)

#### construct reduced model to analyze garden data -- eq. 2
gar1 <- lmer(value~clim*Trait+(1|Pop)+(1|Pop:Trait)+(1|Pop:Trait:Fam), data=g3)
summary(gar1)
anova(gar1)
gar1a <- emtrends(gar1, ~Trait, var="clim") ### examine slopes for each trait with bioclimPC
gar1b <- as.data.frame(summary(gar1a))     ### save slopes in dataframe for figure (fig. 3D)


########################################################################################################
##### RAW CORRELATIONS #############
#########################################################################################################

## NATURAL POP CORRS
pt <- pop1[6:14]

rcorr(as.matrix(pt))  #### Table E1 above diagonal
pt1 <- na.omit(pop1)
pt2 <- pt1[6:14]

### conduct PCA on traits in natural populations -- Table E4
ptpc <- prcomp(pt2, scale=TRUE, center=TRUE)
summary(ptpc)                               
ptpc
plot(ptpc)
biplot(ptpc, pc.biplot=TRUE)

pt1$pc1 <-  predict(ptpc)[,1]*-1  ### invert axis so that all traits load positively
pt1$pc2 <-  predict(ptpc)[,2]

## save pca output for making graph later
pcrot <- as.data.frame(ptpc$rotation[,1:2], stringsAsFactors=TRUE)
pcrot$Trait <- c("Hgt","Ram","Pods","SLA","d13C","Lf_N","Tric","Card","Ltx")
pcrot$PC1s <- pcrot$PC1*-5 #scale loadings for plotting as vector arrows
pcrot$PC2s <- pcrot$PC2*5 #scale loadings for plotting as vector arrows

##### garden corrs
rcorr(as.matrix(g2ma[5:13]))  ## table E1 below diagonal

cor.test(g2ma$Cards[g2ma$Cards<1.8],g2ma$Latex[g2ma$Cards<1.8])  ## remove outlier for cards~latex corr

gspc <- prcomp(g2ma[5:13], scale=FALSE, center=FALSE)
summary(gspc)
gspc
plot(gspc)
biplot(gspc, pc.biplot=TRUE)
g2ma$pc1 <-  predict(gspc)[,1]
g2ma$pc2 <-  predict(gspc)[,2]

### save pca output for making graph later
gcrot <- as.data.frame(gspc$rotation[,1:2], stringsAsFactors=TRUE)
gcrot$Trait <- c("Hgt","Ram","Pods","SLA","dC13","Lf_N","Tric","Card","Ltx")
gcrot$PC1s <- gcrot$PC1*5 #scale loadings for plotting as vector arrows
gcrot$PC2s <- gcrot$PC2*5 #scale loadings for plotting as vector arrows






####################################################################################
########### CREATE FIGURES #################################
############################################################

######################################
##### Beetle FIGURE - fig 2A-B

tiff("Fig2_beetlesdamage.tiff", width=3.5, height=5, units='in',
     res=300, compression='lzw')

par(mfrow=c(2,1), cex=1, oma=c(.5,.5,.5,.5),mar=c(4,3.5,.5,.5),mgp=c(2.5,1,0))
plot(log(beetles+1)~clim, data=pop, xlab="Bioclim PC",
     ylab="log # beetles", pch=21, bg=Col, cex=1.25, lwd=2, las=1)
abline(b1, lwd=2)
plot(Damage~log(beetles+1), data=pop, xlab="log # beetles",
     ylab="% leaf damage", pch=21, bg=Col, cex=1.25, lwd=2, las=1)
abline(b2, lwd=2)
legend("bottomright",c("-4","0","4"), col=c('blue','green','red'), 
       title="Bioclim PC", xpd = TRUE, pch=16, cex=.8,horiz=TRUE)

dev.off()
## end beetle plot

###################################
### fig. 3A natural pops dotchart
tiff("Fig3A_FieldDot.tiff", width=3.5, height=4, units='in',
     res=300, compression='lzw')
par(mfrow=c(1,1), cex=1, oma=c(.5,.5,.5,.5),mar=c(4,3.5,3,.5),mgp=c(2.5,1,0))
dotchart(f2b$clim.trend, labels=f2b$Trait, pch=16, pt.cex=1.25, xlim=range(c(-.2,.35)),
         color=c("red","red","black","black","red","black","black","red","red"),
         xlab="Coefficient with Bioclim PC", main="Natural")
abline(v=0, col="grey", lwd=2)
for(i in 1:10){
  lines(x=c(f2b$lower.CL[i],f2b$upper.CL[i]), y=c(i,i), lwd=2)
}

dev.off()
### end dotchart

######### GARDEN DOT PLOT #################################

tiff("Fig3b_GardenDot.tiff", width=3.5, height=4, units='in',
     res=300, compression='lzw')
par(mfrow=c(1,1), cex=1, oma=c(.5,.5,.5,.5),mar=c(4,3.5,3,.5),mgp=c(2.5,1,0))
dotchart(gar1b$clim.trend, labels=gar1b$Trait, pch=16, pt.cex=1.25, xlim=range(c(-.2,.35)),
         color=c("black","black","black","red","black",
                 "black","black","red","red"),
         xlab="Coefficient with Bioclim PC", main="Garden")
abline(v=0, col="grey", lwd=2)
for(i in 1:10){
  lines(x=c(gar1b$lower.CL[i],gar1b$upper.CL[i]), y=c(i,i), lwd=2)
}
dev.off()
### end dotchart

#########################################################################
###### CLINE PLOTS ############################
#########################################################################
############ PLOT CLINES for Height, d13C, Cards, latex #######################

#### mods - backtransform variables for plots
fl1 <- lm((exp(No_ram)-1)~clim, data=popgraph)
fl2 <- lm(Cards**2~clim, data=popgraph)
gr1 <- lm((exp(No_ram)-1)~clim, data=gargraph)
gr2 <- lm(Cards**2~clim, data=gargraph)

tiff("Fig3_FGCline.tiff", width=6.5, height=5, units='in',
     res=300, compression='lzw')

par(mfrow=c(2,2), cex=1, oma=c(2,.5,.5,.5),mar=c(3.5,3.5,1,.5),mgp=c(2.5,1,0))
plot((exp(No_ram)-1)~clim, data=popgraph, xlab="",xlim=c(-6,4.3),ylim=c(0,9),
     ylab="No. ramets", pch=21, bg=Col, cex=1.25, lwd=2, las=1)
abline(fl1, lwd=2)
text(-6,8.5,"B")
plot((exp(No_ram)-1)~clim, data=gargraph, xlab="",xlim=c(-6,4.3),ylim=c(0,9),
     ylab="No. ramets", pch=21, bg=Col, cex=1.25, lwd=2, las=1)
#abline(g1, lwd=2)
text(-6,8.5,"E")

plot(Cards**2~clim, data=popgraph, xlab="",xlim=c(-6,4.3),ylim=c(0,.85),
     ylab="Cardenolides (mg/g)", pch=21, bg=Col, cex=1.25, lwd=2, las=1)
abline(fl2, lwd=2)
text(-5.95,.8325,"C")
plot(Cards**2~clim, data=gargraph, xlab="",xlim=c(-6,4.3),ylim=c(0,.3),
     ylab="Cardenolides (mg/g)", pch=21, bg=Col, cex=1.25, lwd=2, yaxt="n")
axis(side=2, at=seq(0,.3,by=.1), las=1)
abline(gr2, lwd=2)
mtext("Bioclim PC", side=1, outer=TRUE, line=-.5, cex=1.25)
legend("topleft",c("-4","0","4"), col=c('blue','green','red'), 
       title="Bioclim PC", xpd = TRUE, pch=16, cex=.8,horiz=TRUE)
text(-6,.293,"F")

dev.off()



#############
##### TRIAT CORRELATION FIGURE ######
tiff("Fig4a_FieldTraits.tiff", width=3.5, height=10, units='in',
     res=300, compression='lzw')

par(mfrow=c(4,1), cex=1, oma=c(1,1,.5,.5),mar=c(3.5,3.5,1,.5),mgp=c(2.5,1,0), xpd=FALSE)
plot(pc2~pc1, data=pt1, pch=21, bg=Col, ylim=range(c(-3,3.1)),xlim=range(c(-3.5,3.3)), 
     xlab="Trait PC1 (36.0%)",ylab="Trait PC2 (21.5%)", main="Natural", las=1)
abline(a=0,b=0,lty=3)
abline(h=0,v=0,lty=3)
arrows(0,0,pcrot$PC1s,pcrot$PC2s, length=.1, lwd=2, col="grey")
text(pcrot$PC1s,pcrot$PC2s, labels=pcrot$Trait, pos=4, font=2, cex=.8, col="black")
text(-3.1,3,"A")

plot(Latex~Height, data=popgraph, xlab="Stem height (cm)", xlim=c(40,110), ylim=c(2,4.5),
     ylab="log(Latex (mg))", pch=21, bg=Col, cex=1.25, lwd=2, las=1)
text(41,4.4,"B")

plot(Cards**2~Height, data=popgraph, xlab="Stem height (cm)",xlim=c(40,110), ylim=c(0,.85),
     ylab="Cardenolides (mg/g)", pch=21, bg=Col, cex=1.25, lwd=2, las=1)
text(41,.83,"C")

plot(Cards**2~Latex, data=popgraph, xlab="log(Latex (mg))",xlim=c(2,4.5), ylim=c(0,.85),
     ylab="Cardenolides (mg/g)", pch=21, bg=Col, cex=1.25, lwd=2, las=1)
text(2.05,.83,"D")

dev.off()


#####################################################
#### GARDEN pop averages


##### TRIAT CORRELATION FIGURE ######
tiff("Fig4b_GardenTraits.tiff", width=3.5, height=10, units='in',
     res=300, compression='lzw')

par(mfrow=c(4,1), cex=1, oma=c(1,1,.5,.5),mar=c(3.5,4,1,.5),mgp=c(2.5,1,0), xpd=FALSE)
plot(pc2~pc1, data=g2ma, pch=21, bg=Col,cex=1, las=1, 
     ylim=range(c(-3,3.1)),xlim=range(c(-3.5,3.3)),
     xlab="Trait PC1 (33.8%)",ylab="Trait PC2 (24.0%)",main="Garden")
abline(a=0,b=0,lty=3)
abline(h=0,v=0,lty=3)
arrows(0,0,gcrot$PC1s,gcrot$PC2s, length=.1, lwd=2, col="grey")
text(gcrot$PC1s,gcrot$PC2s, labels=gcrot$Trait, pos=3,
     font=2, cex=.8, col="black")
text(-3.1,3,"E")
plot(Latex~(exp(No_ram)-1), data=gargraph, pch=21, bg=Col,cex=1.25,lwd=2, xlim=c(2.5,8), ylim=c(2.6,3.7),
     xlab="No. ramets",ylab="log(Latex (mg))",las=1)
text(2.5,3.65,"F")

plot(Cards**2~(exp(No_ram)-1), data=gargraph, pch=21, bg=Col,cex=1.25, lwd=2,xlim=c(2.5,8), ylim=c(0,.25), 
     xlab="No. ramets",ylab="Cardenolides (mg/g)",yaxt="n")
axis(side=2, at=seq(0,.3,by=.1), las=1)
text(2.5,.24,"G")

plot(Cards**2~Latex, data=gargraph, pch=21, bg=Col,cex=1.25,lwd=2,xlim=c(2.6,3.7), ylim=c(-.02,.25), 
     xlab="log(Latex (mg))",ylab="Cardenolides (mg/g)",yaxt="n")
points(Cards**2~Latex, data=gargraph, pch=4, subset=Cards>.45)
axis(side=2, at=seq(0,.3,by=.1), las=1)
text(2.61,.24,"H")
legend("bottomright",c("-4","0","4"), col=c('blue','green','red'), 
       title="Bioclim PC", pch=16, cex=.8,horiz=TRUE)
dev.off()

### end plot

rcorr(as.matrix(gargraph[5:13]))

lmtc <- lm(Cards~Latex, data=g2ma) #, subset=Cards<1.8
summary(lmtc)
outlierTest(lmtc)
cooks.distance(lmtc)

lmtc2 <- lm(Cards~Latex, data=g2ma, subset=Cards<1.8) #, subset=Cards<1.8
summary(lmtc2)
cor.test(g2ma$Cards[g2ma$Cards<1.8], g2ma$Latex[g2ma$Cards<1.8])
plot(Cards~Latex, data=g2ma, subset=Cards<1.8)


########################################################################################
###### APPENDIX ANALYSIS: CLIMATE AND SIZED CORRECTIONS  ###############################
########################################################################################

#############################################################
##### CLIMATE-CORRECTED CORRELATIONS #############

#### natural pops clim corrected pca
pop3 <- na.omit(pop2)
r1a <- lmer(value~clim*Trait+(1|Pop), data=pop3)  # construct model to control for climate (same model as f1 above)

pop3$resid <- resid(r1a)  ## extract residuals from climate-corrected model

pop4 <- dcast(pop3, Pop+clim+soils+beetles~Trait, value.var="resid")


rcorr(as.matrix(pop4[5:13]))  ##### Table E2 - above diagonal
pop5 <- as.matrix(na.omit(pop4[5:13]))

ccpc <- prcomp(pop5, scale=F, center=F)  ## Table E4: Natural pop trait correlations, after climate-correction 
summary(ccpc)
ccpc
plot(ccpc)
biplot(ccpc, pc.biplot=TRUE)
ccpc$pc1 <-  predict(ccpc)[,1]*-1
ccpc$pc2 <-  predict(ccpc)[,2]


#### garden climate-corrected
r1g <- lmer(value~clim*Trait+(1|Pop), data=g3m) ### construct model to correct for climate using population-average data

g3m$resid <- resid(r1g)

g4m <- dcast(g3m, Pop+clim+soils+beetles~Trait, value.var="resid")

rcorr(as.matrix(g4m[5:13]))  #### Table E2 - below diagonal

gcpc <- prcomp(g4m[5:13], scale=F, center=F)  ## Table E4: Garden trait correlations, after climate-correction
summary(gcpc)
gcpc
plot(gcpc)
biplot(gcpc, pc.biplot=TRUE)
gcpc$pc1 <-  predict(gcpc)[,1]*-1
gcpc$pc2 <-  predict(gcpc)[,2]

latcard <- lm(Cards~Latex, data=g4m, subset=Pop!='Arl')
summary(latcard)


#############################################################################
####### SIZE-CORRECTED CORRELATIONS ########################################

#### SIZE CORRECTION FIELD
### wide to long
size0 <- melt(pop1, id.vars=c("Pop","clim","soils","beetles","Col","Height","No_ram"), 
              measure.vars=c("Pods","SLA","dC13","LeafN",
                             "Trichome","Cards","Latex"))
colnames(size0)[8] <- "Trait"
head(size0)

### construct models to correct for size (size is stem height and no. of stems)
size0 <- na.omit(size0)
size1 <- lmer(value~Height*Trait+No_ram*Trait+(1|Pop), data=size0)
summary(size1)
anova(size1)

size0$sizeres <- resid(size1)

size2 <- dcast(size0, Pop+clim+soils+beetles~Trait, value.var="sizeres")


rcorr(as.matrix(size2[5:11])) ### Table E3 - above diagonal
size3 <- na.omit(size2[5:11])

spc <- prcomp(size3, scale=TRUE, center=TRUE)  ### Table E4
summary(spc)
spc
plot(spc)
biplot(spc, pc.biplot=TRUE)
spc$pc1 <-  predict(spc)[,1]*-1
spc$pc2 <-  predict(spc)[,2]



#### SIZE CORRECTION GARDEN
### wide to long using population averages
sizeg <- melt(g2ma, id.vars=c("Pop","clim","soils","beetles","Col","Height","No_ram"), 
              measure.vars=c("Pods","SLA","dC13","LeafN",
                             "Trichome","Cards","Latex"))
colnames(sizeg)[8] <- "Trait"
head(sizeg)

### construct models to correct for size (size is stem height and no. of stems)
sizeg1 <- lmer(value~Height*Trait+No_ram*Trait+(1|Pop), data=sizeg)
summary(sizeg1)
anova(sizeg1)

sizeg$sizeres <- resid(sizeg1)

sizeg2 <- dcast(sizeg, Pop+clim+soils+beetles~Trait, value.var="sizeres")


rcorr(as.matrix(sizeg2[5:11]))  ### Table E3 - below diagonal

sgpc <- prcomp(sizeg2[5:11])   ### Table E4
summary(sgpc)
sgpc
plot(sgpc)
biplot(sgpc, pc.biplot=TRUE)
sgpc$pc1 <-  predict(sgpc)[,1]*-1
sgpc$pc2 <-  predict(sgpc)[,2]




