library(ape); library(phytools); library(MCMCglmm)
library(metafor); library(cowplot);library(plotMCMC)
library(plyr); library(emmeans); library(dplyr);library(cowplot)

setwd("C:/Users/Dochtermann/Dropbox/Working/Projects/AGA/H2 Revision & Reanalysis")
Data <- read.csv("Compiled2.csv")

Tree<-read.newick("newickNEW.txt")

#### 1. Clean up tree & plotting the tree ####
Tree$node.label<-NULL; Tree<-collapse.singles(Tree) 
inv.phylo<-inverseA(Tree,nodes="TIPS",scale=TRUE)

tips<-Tree$tip.label; genera<-unique(sapply(strsplit(tips,"_"),function(x) x[1]))
ii<-sapply(genera,function(x,y) grep(x,y)[1],y=tips)
tree<-drop.tip(Tree,setdiff(Tree$tip.label,tips[ii]))
#plotTree(tree,ftype="i")
tree$tip.label<-sapply(strsplit(tree$tip.label,"_"),function(x) x[1])
#plotTree(tree,ftype="i")
#plot(tree,"radial")

#### 2. Global Mean, Variance Components & Heterogeneity (Intercept Only) ####
CorMatrix <- cov2cor(as.matrix(inv.phylo$Ainv))

Model1 <- rma(Zr, vi=1/(Data$N -3), data=Data,method="REML")
Model1

Model2.1 <- rma.mv(Zr, V=1/(Data$N -3), 
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="REML")

Model2.2 <- rma.mv(Zr, V=1/(Data$N -3), 
                   random = list(~ 1 | Study.ID, ~ 1 | phylo,~1|Estimate.ID),
                   R=list(phylo=CorMatrix),
                   data=Data, method="REML")

Model2.2b <- rma.mv(Zr, V=1/(Data$N.C -3), 
                   random = list(~ 1 | Study.ID, ~ 1 | phylo,~1|Estimate.ID),
                   R=list(phylo=CorMatrix),
                   data=Data, method="REML")

Model2.1

cbind(Model2.1$b,
      Model2.1$ci.lb,
      Model2.1$ci.ub)

cbind(transf.ztor(Model2.1$b),
      transf.ztor(Model2.1$ci.lb),
      transf.ztor(Model2.1$ci.ub))

s2m <- sum(Data$N-3)*41/(sum(Data$N-3)^2-sum((Data$N-3)^2)) # typical sampling error
s2t <- sum(Model2.2$sigma2)+s2m # total sampling error

I2t <- sum(Model2.2$sigma2)/s2t
I2t*100 # 99.79%

I2s <- Model2.2$sigma2[1]/s2t
I2s*100 # 41.43% varaince due to Study

I2p <- Model2.2$sigma2[2]/s2t
I2p*100 # 0% varaince due to Phylogeny

I2e <- Model2.2$sigma2[3]/s2t
I2e*100 # 58.36% - residuals against sampling error

I2se <- s2m/s2t
I2se*100

### 3. Tests of overall moderator effects ####
Model3.1 <- rma.mv(Zr, V=1/(Data$N -3), mod=~D.W+L.F+Behavior+Thermy+Vert,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="REML")

Model3.1ML <- rma.mv(Zr, V=1/(Data$N -3), mod=~D.W+L.F+Behavior+Thermy+Vert,
                     random = list(~ 1 | Study.ID, ~ 1 | phylo),
                     R=list(phylo=CorMatrix),
                     data=Data, method="ML")

anova(Model3.1) #Global ANOVA
Model3.1

#Significance tests by category (run as d-i-d, hence ML fitting)
Model3.2 <- rma.mv(Zr, V=1/(Data$N -3), mod=~L.F+Behavior+Thermy+Vert,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="ML")

Model3.3 <- rma.mv(Zr, V=1/(Data$N -3), mod=~D.W+Behavior+Thermy+Vert,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="ML")

Model3.4 <- rma.mv(Zr, V=1/(Data$N -3), mod=~D.W+L.F+Thermy+Vert,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="ML")

Model3.5 <- rma.mv(Zr, V=1/(Data$N -3), mod=~D.W+L.F+Behavior+Vert,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="ML")

Model3.6 <- rma.mv(Zr, V=1/(Data$N -3), mod=~D.W+L.F+Behavior+Thermy,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="ML")

anova(Model3.1ML,Model3.2) #Wild vs Domestic
anova(Model3.1ML,Model3.3) #Field vs Lab
anova(Model3.1ML,Model3.4)  #Behavior
anova(Model3.1ML,Model3.5) #Thermy
anova(Model3.1ML,Model3.6) #Vert

#### 3b. Sample Size Uncertainty ####

#These analyses were just done with a sample size (N.C) set to a small value (N=25).
#This was done to check whether inferences were sensitive to possible mis-estimation of 
#sampling variance. The inferences were all the same.

Model3.1bML <- rma.mv(Zr, V=1/(Data$N.C-3), mod=~D.W+L.F+Behavior+Thermy+Vert,
                     random = list(~ 1 | Study.ID, ~ 1 | phylo),
                     R=list(phylo=CorMatrix),
                     data=Data, method="ML")

#Significance tests by category (run as d-i-d, hence ML fitting)
Model3.2b <- rma.mv(Zr, V=1/(Data$N.C-3), mod=~L.F+Behavior+Thermy+Vert,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="ML")

Model3.3b <- rma.mv(Zr, V=1/(Data$N.C-3), mod=~D.W+Behavior+Thermy+Vert,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="ML")

Model3.4b <- rma.mv(Zr, V=1/(Data$N.C-3), mod=~D.W+L.F+Thermy+Vert,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="ML")

Model3.5b <- rma.mv(Zr, V=1/(Data$N.C-3), mod=~D.W+L.F+Behavior+Vert,
                    random = list(~ 1 | Study.ID, ~ 1 | phylo),
                    R=list(phylo=CorMatrix),
                    data=Data, method="ML")

Model3.6b <- rma.mv(Zr, V=1/(Data$N.C-3), mod=~D.W+L.F+Behavior+Thermy,
                    random = list(~ 1 | Study.ID, ~ 1 | phylo),
                    R=list(phylo=CorMatrix),
                    data=Data, method="ML")

anova(Model3.1bML,Model3.2b)
anova(Model3.1bML,Model3.3b)
anova(Model3.1bML,Model3.4b)
anova(Model3.1bML,Model3.5b)
anova(Model3.1bML,Model3.6b)

### 4. Category specific estimates ####
Model4.1 <- rma.mv(Zr, V=1/(Data$N -3), mod=~D.W-1,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="REML")

Model4.2 <- rma.mv(Zr, V=1/(Data$N -3), mod=~L.F-1,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="REML")

Model4.3 <- rma.mv(Zr, V=1/(Data$N -3), mod=~Behavior-1,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="REML")

Model4.4 <- rma.mv(Zr, V=1/(Data$N -3), mod=~Thermy-1,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="REML")

Model4.5 <- rma.mv(Zr, V=1/(Data$N -3), mod=~Vert-1,
                   random = list(~ 1 | Study.ID, ~ 1 | phylo),
                   R=list(phylo=CorMatrix),
                   data=Data, method="REML")

cbind(transf.ztor(Model4.1$b),
      transf.ztor(Model4.1$ci.lb),
      transf.ztor(Model4.1$ci.ub))

cbind(transf.ztor(Model4.2$b),
      transf.ztor(Model4.2$ci.lb),
      transf.ztor(Model4.2$ci.ub))

cbind(transf.ztor(Model4.3$b),
      transf.ztor(Model4.3$ci.lb),
      transf.ztor(Model4.3$ci.ub))

cbind(transf.ztor(Model4.4$b),
      transf.ztor(Model4.4$ci.lb),
      transf.ztor(Model4.4$ci.ub))

cbind(transf.ztor(Model4.5$b),
      transf.ztor(Model4.5$ci.lb),
      transf.ztor(Model4.5$ci.ub))

#### 5. Species specific estimates ####
Model5<- rma.mv(Zr, V=1/(Data$N -3), mod=~phylo-1,
                random = list(~ 1 | Study.ID), 
                data=Data, method="REML")

Meta.Species=summary(Model5)
Meta.DF=data.frame(Sp=substr(attr(Meta.Species$b,"dimnames")[[1]],6,36),
                   M = Meta.Species$b, CI.lb = Meta.Species$ci.lb, CI.ub = Meta.Species$ci.ub)

Meta.DF <- Meta.DF[match(Tree$tip.label, Meta.DF$Sp),]

#### 6. Eggert regression for publication bias ####

trimfill(Model1)

DATAnoNA<-c(na.omit(Data$N))
Resid.all=residuals(Model1)
Precision.all=sqrt(1/(DATAnoNA -3))
model <- rma(yi=Resid.all,sei=1/Precision.all)
regtest(model,model="lm")
#funnel(model,yaxi="seinv")

#### 7. Plotting ####
#By Moderator Category (Fig1)

ForFigs <- read_excel("ForFigs.xlsx")
ForFigs<-as.data.frame(ForFigs)

ForFigs$Category2<-factor(ForFigs$Category2,
         levels=rev(c(ForFigs$Category2)))
ForFigs$Category2


Fig2 <- ggplot(ForFigs,aes(ES,Category2)) + 
  geom_errorbarh(aes(xmax = upperci, xmin = lowerci), height = 0.15) +
  geom_point(aes(size=N,shape=Category2,fill=Category2)) + #scale_size_identity(trans="sqrt") +
  scale_shape_manual(values=c(23,rep(16,21))) + scale_fill_manual(values=c("white",rep("black",21))) +
  geom_vline(xintercept = 0, linetype = "longdash") + 
  scale_x_continuous(breaks = seq(-.5,.7,.1)) +
  labs(x=expression('h'^2),y="") +
  theme(axis.title.x  = element_text(size=20),
        axis.text.x  = element_text(size=14),
        axis.text.y  = element_text(size=16),
        axis.line.y=element_blank(),axis.ticks.y = element_blank(),
        legend.position="none")
Fig2 
#tiff(file="Fig2.tiff",width=7,height=11,units='in',res=600)
Fig2
#dev.off()


#By phylogeny (Fig2)
Meta.DF$M_r=transf.ztor(Meta.DF$M)
Meta.DF$CI.lb_r=transf.ztor(Meta.DF$CI.lb)
Meta.DF$CI.ub_r=transf.ztor(Meta.DF$CI.ub)
Meta.DF$CI.ub_r=transf.ztor(Meta.DF$CI.ub)
Meta.DF[which(Meta.DF$M_r<0),5]<-0

Meta.DF$Sp2=as.factor(c("Lobesia botrana","Grapholita molesta",
                        "Achroia grisella","Diabrotica undecimpunctata",
                        "Tribolium castaneum","Nicrophorus vespilloides",
                        "Musca domestica","Drosophila ananassae",
                        "Drosophila melanogaster","Sepsis cynipsea",
                        "Panorpa vulgaris","Cotesia glomerata",
                        "Nasonia vitripennis","Gryllus integer",
                        "Aquarius remigis","Enallagma hageni",
                        "Enallagma geminatum","Larinioides sclopetarius",
                        "Agelenopsis pennsylvanica","Hydroides dianthus",
                        "Argopecten purpuratus","Euprymna tasmanica",
                        "Rana sylvatica","Rhinella marina",
                        "Plethodon cinereus","Papio anubis×cynocephalus",
                        "Macaca mulatta","Pongo pygmaeus×abelii",
                        "Pan troglodytes","Tamiasciurus hudsonicus",
                        "Marmota flaviventris","Tamias striatus",
                        "Mus musculus","Rattus norvegicus",
                        "Phodopus sungorus","Peromyscus maniculatus",
                        "Ovis aries","Ovis canadensis",
                        "Bos taurus","Sus scrofa",
                        "Crocuta crocuta","Ficedula albicollis",
                        "Sialia mexicana","Serinus canaria",
                        "Passer domesticus","Taeniopygia guttata",
                        "Cyanistes caeruleus","Parus major",
                        "Acrocephalus sechellensis","Aegithalos caudatus",
                        "Hirundo rustica","Petrochelidon pyrrhonota",
                        "Falco naumanni","Tachymarptis melba",
                        "Picoides borealis","Gallus gallus",
                        "Sterna hirundo","Larus michahellis",
                        "Chrysemys picta","Zootoca vivipara",
                        "Eulamprus quoyii","Amatitlania siquia",
                        "Neolamprologus pulcher","Lucania goodei",
                        "Girardinichthys multiradiatus","Xiphophorus birchmanni",
                        "Poecilia reticulata","Gasterosteus aculeatus",
                        "Salmo salar","Salmo trutta",
                        "Danio rerio"))

#tiff(file="Fig1.tiff",width=7,height=11,units='in',res=600)
par(mfrow=c(1,2))
par(mar=c(4.5,2,1.35,1))
plot(Tree, font=1, cex=0.8, x.lim=90, show.tip.label = F)
par(mar=c(4,3,1,0))
plot(Meta.DF$M_r, 1:length(Meta.DF$M_r), ylab=NA, yaxt="n", bty="n", 
     xlim=c(-.6,1), ylim=c(0.25, length(Meta.DF$M_r)+.5), xlab=expression('h'^2),
     pch=16, cex=0.8, cex.axis=.9)
abline(v=0,lty=3)
mtext(Meta.DF$Sp2, 2, -1, at=1:length(Meta.DF$M), las=2, cex=.8, font=3)
segments(Meta.DF$CI.lb_r, 1:length(Meta.DF$M_r), 
         Meta.DF$CI.ub_r, 1:length(Meta.DF$M_r), lwd=1.25)
#dev.off()

#Publication Bias
Fig3=ggplot(Data, aes(x=Zr, y=sqrt(Data$N-3))) + geom_point(shape=1) +
  geom_vline(xintercept = 0) + geom_vline(xintercept = 0.2891878, linetype = "longdash") +
  labs(x="Zr", y="Precision (1/SE)") + xlim(-1.6,1.6) + 
  scale_y_continuous(breaks=seq(0,110,10), limits=c(0,110))
Fig3

#tiff(file="Fig3.tiff",width=6,height=6,units='in',res=600)
Fig3
#dev.off()