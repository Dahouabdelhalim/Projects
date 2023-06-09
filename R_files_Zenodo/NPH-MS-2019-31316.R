#============
#INPUT DATA
#============
p <- read.table("HB_Data_12-13-19.txt",sep="\\t",header=T)
str(p)

#Set desired variables to factors
p$Round <- as.factor(p$Round)
p$BoxNo <- as.factor(p$BoxNo)
p$Acc <- as.factor(p$Acc)
p$Pouch <- as.factor(p$Pouch)
str(p)

#================================================================
#SPATIAL AND TEMPORAL EFFECTS ON SHOOT BIOMASS AND NODULE TRAITS
#================================================================
#Caclulate trait mean of replicate pouches within accession, including controls
p2=ddply(p,c("Strain", "HostSpecies","Acc","Ploidy", "BoxNo", "Round"),
         summarise, Shoot.mean=mean(Shoot, na.rm = TRUE),
         Total.mean=mean(Total, na.rm=TRUE),
         QNC.mean=mean(QNodColor, na.rm=TRUE),
         QNC2.mean=mean(QNodColor2, na.rm = TRUE),
         NodNo.mean=mean(NodNo, na.rm=TRUE),
         NodBM.mean=mean(NodBMmg, na.rm=TRUE))
head(p2)

st1 <- lm(p2$Shoot.mean ~ p2$Ploidy * p2$Strain + p2$BoxNo)
summary(st1)
anova(st1)

st2 <- lm(p2$Shoot.mean ~ p2$Ploidy * p2$Strain + p2$Round)
summary(st2)
anova(st2)

st3 <- lm(p2$QNC.mean ~ p2$Ploidy * p2$Strain + p2$BoxNo)
summary(st3)
anova(st3)

st4 <- lm(p2$QNC.mean ~ p2$Ploidy * p2$Strain + p2$Round)
summary(st4)
anova(st4)

st5 <- lm(p2$NodNo.mean ~ p2$Ploidy * p2$Strain + p2$BoxNo)
summary(st5)
anova(st5)

st6 <- lm(p2$NodNo.mean ~ p2$Ploidy * p2$Strain + p2$Round)
summary(st6)
anova(st6)

st7 <- lm(p2$NodBM.mean ~ p2$Ploidy * p2$Strain + p2$BoxNO)
summary(st7)
anova(st7)

st8 <- lm(p2$NodBM.mean ~ p2$Ploidy * p2$Strain + p2$Round)
summary(st8)
anova(st8)

#=======================
#SHOOT BIOMASS
#=======================
#Calculate mean shoot biomass of replicate pouches within accession, including controls
p3=ddply(p,c("Strain", "HostSpecies","Acc","Ploidy"),
         summarise, Shoot.mean=mean(Shoot, na.rm = TRUE),
         Root.mean=mean(Root, na.rm=TRUE),
         Total.mean=mean(Total, na.rm=TRUE),
         QNC.mean=mean(QNodColor, na.rm=TRUE),
         QNC2.mean=mean(QNodColor2, na.rm=TRUE),
         NodNo.mean=mean(NodNo, na.rm=TRUE),
         NodBM.mean=mean(NodBMmg, na.rm=TRUE))
head(p3)

#Shoot Biomass Model
shoot <- lmer(Shoot.mean ~ Ploidy * Strain + (1|HostSpecies/Acc), data = p3)
plot(shoot)
E_shoot = resid(shoot, type = "deviance")
hist(E_shoot)
qqnorm(E_shoot)
qqline(E_shoot)
summary(shoot)
anova(shoot)
shootm <- emmeans(shoot, specs = pairwise ~ Ploidy:Strain)
shootm


#Shoot summary and plot
shoots <- ddply(p3, c("Ploidy", "Strain"), summarise, N = length(Shoot.mean), mean = mean(Shoot.mean),colors=mean(QNC.mean), sd = sd(Shoot.mean), se = sd/sqrt(N))
shoots

shootmean <- tapply(p3$Shoot.mean, p3$Ploidy, mean)
shootmean
shootmax <- tapply(p3$Shoot.mean, p3$Ploidy, max)
shootmax
shootmin <- tapply(p3$Shoot.mean, p3$Ploidy, min)
shootmin

#Shoot Biomass Figure (all 21 strains + controls)
shootsort <- p3[order(p3$Shoot.mean),]
shootsort$Strain=factor(shootsort$Strain,levels=c("Uninoc", "M30", "USDA205", "USDA4894", "USDA207", "KH53b", "M270", "USDA1021", "KH36b", "WSM419", "KH53a", "KH36c", "1021", "KH46c", "KH30a", "KH16b", "USDA1002", "KH36d", "M210", "KH48e", "KH35c", "USDA4893"))

shootp <- ggplot(shootsort, aes(x = Strain, y = Shoot.mean, group = Ploidy, colour = Ploidy)) + 
  geom_point(size=0.6, alpha=0.6) +
  stat_summary(fun.y=mean, geom="line", size=1) +
  stat_summary(fun.y=mean, geom="point", size=2) +
  xlab("Sinorhizobium Strain") + 
  ylab("Shoot Biomass (g)") + theme_classic() + scale_color_manual(values=c("gray50", "black")) + theme(axis.text.x=element_text(angle=90, hjust=1))
shootp

#Nodule Color Matrix - Shoot Biomass
shootsort2 <- shoots[order(shoots$colors),]
shootsort2

nod.col=matrix(ncol=nrow(shootsort2)/2,nrow=2, dimnames=list(c("2X", "4X")))
nod.col[1,]=shootsort2$colors[shootsort2$Ploidy=="2X"]
nod.col[2,]=shootsort2$colors[shootsort2$Ploidy=="4X"]

col2 <- colorRampPalette(brewer.pal(9,"YlOrRd"))(50)
levelplot(nod.col, col.regions= col2)
heatmap(nod.col, Colv=NA, Rowv=NA, col=col2)

#=========================================
#Test of Lineage Effects - Shoot Biomass
#=========================================
shoot2 <- lm(Shoot.mean ~ Ploidy * Strain, data = p3)
plot(shoot2)
E_shoot2 = resid(shoot2, type = "deviance")
hist(E_shoot2)
qqnorm(E_shoot2)
qqline(E_shoot2)
summary(shoot2)
anova(shoot2)

#Compare models
anova(shoot,shoot2)

#Variation explained by accession and host species
Var_Random_effect <- VarCorr(shoot)
Var_Random_effect
variances <- as.data.frame(Var_Random_effect)
variances

#accession
variances$vcov[1]/sum(variances$vcov)

#host species
variances$vcov[2]/sum(variances$vcov)

#====================
#TOTAL BIOMASS
#====================
total <- lmer(Total.mean ~ Ploidy * Strain + (1|HostSpecies/Acc), data = p3)
plot(total)
E_total = resid(total, type = "deviance")
hist(E_total)
qqnorm(E_total)
qqline(E_total)
summary(total)
anova(total)
totm <- emmeans(total, specs = pairwise ~ Ploidy:Strain)
totm

totals<- ddply(p3, c("Ploidy", "Strain"), summarise, N = length(Total.mean), mean = mean(Total.mean),colors=mean(QNC.mean), sd = sd(Total.mean), se = sd/sqrt(N))
totals

totalmean <- tapply(p3$Total.mean, p3$Ploidy, mean)
totalmean
totalmax <- tapply(p3$Total.mean, p3$Ploidy, max)
totalmax
totalmin <- tapply(p3$Total.mean, p3$Ploidy, min)
totalmin

#Total Biomass Figure (all 21 strains + controls)
totalsort <- p3[order(p3$Total.mean),]
totalsort$Strain=factor(totalsort$Strain,levels=c("Uninoc", "M30", "USDA205", "USDA4894", "USDA207", "KH53b", "M270", "USDA1021", "KH36b", "WSM419", "KH53a", "KH36c", "1021", "KH46c", "KH30a", "KH16b", "USDA1002", "KH36d", "M210", "KH48e", "KH35c", "USDA4893"))

totalp <- ggplot(totalsort, aes(x = Strain, y = Total.mean, group = Ploidy, colour = Ploidy)) + 
  geom_point(size=0.6, alpha=0.6) +
  stat_summary(fun.y=mean, geom="line", size=1) +
  stat_summary(fun.y=mean, geom="point", size=2) +
  xlab("Sinorhizobium Strain") + ylab("Total Shoot and Root Biomass (mg)") + theme_classic() + scale_color_manual(values=c("gray50", "black")) + theme(axis.text.x=element_text(angle=90, hjust=1))
totalp

#======================================
#HOST GROWTH RESPONSE OF SHOOT BIOMASS
#======================================
#Calculate Average Shoot Biomass for Controls within Round
contr=p[p$Strain=="Uninoc",]
head(contr)

p4<-ddply(contr,c("HostSpecies","Acc","Ploidy","Round"),
          summarise,SBM.Contr=mean(ShootBMg, na.rm = TRUE))
head(p4)

#Merge dataframes
p5 = merge(p,p4)
head(p5)

#Check that new dataframe has same number of rows as old one
nrow(p)==nrow(p5)

#Calculate Host Growth Response of Shoot Biomass (HGR)
p5$HGR.SBM=((p5$ShootBMg-p5$SBM.Contr)/p5$SBM.Contr)*100
head(p5)

#Remove uninoculated plants from dataframe
p6=p5[p5$Strain!="Uninoc",]
head(p6)

#Calculate mean HGR of replicate pouches within accession and strain (USE FOR ALL 21 STRAINS)
p7<-ddply(p6,c("Strain", "HostSpecies","Acc","Ploidy"),
          summarise,HGRS.mean=mean(HGR.SBM, na.rm = TRUE),
          Shoot.mean=mean(Shoot, na.rm=TRUE),
          Root.mean=mean(Shoot, na.rm=TRUE),
          QNC.mean=mean(QNodColor, na.rm=TRUE),
          QNC2.mean=mean(QNodColor2, na.rm=TRUE),
          NodNo.mean=mean(NodNo, na.rm=TRUE),
          NodBM.mean=mean(NodBMmg, na.rm=TRUE))
head(p7)

#HGR-SHOOT MODEL
hgr <- lmer(HGRS.mean ~ Ploidy * Strain + (1|HostSpecies/Acc), data = p7)
plot(hgr)
E_hgr = resid(hgr, type = "deviance")
hist(E_hgr)
qqnorm(E_hgr)
qqline(E_hgr)
summary(hgr)
anova(hgr)
hgrm <- emmeans(hgr, specs = pairwise ~ Ploidy:Strain)
hgrm

hgrs <- ddply(p7, c("Ploidy", "Strain"), summarise, N = length(HGRS.mean), mean = mean(HGRS.mean),colors=mean(QNC.mean), sd = sd(HGRS.mean), se = sd/sqrt(N))
hgrs

hgrmean <- tapply(p7$HGRS.mean, p7$Ploidy, mean)
hgrmean
hgrmax <- tapply(p7$HGRS.mean, p7$Ploidy, max)
hgrmax
hgrmin <- tapply(p7$HGRS.mean, p7$Ploidy, min)
hgrmin

#Host Growth Response Figure (all 21 strains) 
hgrsort <- p7[order(p7$HGRS.mean),]
hgrsort$Strain=factor(hgrsort$Strain,levels=c("Uninoc", "M30", "USDA205", "USDA4894", "USDA207", "KH53b", "M270", "USDA1021", "KH36b", "WSM419", "KH53a", "KH36c", "1021", "KH46c", "KH30a", "KH16b", "USDA1002", "KH36d", "M210", "KH48e", "KH35c", "USDA4893"))

hgrp <- ggplot(hgrsort, aes(x = Strain, y = HGRS.mean, group = Ploidy, colour = Ploidy)) + 
  geom_point(size=0.6, alpha=0.6) +
  stat_summary(fun.y=mean, geom="line", size=1) +
  stat_summary(fun.y=mean, geom="point", size=2) +
  geom_hline(yintercept=0, linetype=2, color="black") +
  xlab("Sinorhizobium Strain") + ylab("Host Growth Response - Shoot Biomass") + theme_classic() + scale_color_manual(values=c("gray50", "black")) + theme(axis.text.x=element_text(angle=90, hjust=1))
hgrp

#Nodule Color Matrix - HGR Shoot 
hgrsort2 <- hgrs[order(hgrs$colors),]
hgrsort2

nod.col2=matrix(ncol=nrow(hgrsort2)/2,nrow=2, dimnames=list(c("2X", "4X")))
nod.col2[1,]=hgrsort2$colors[hgrsort2$Ploidy=="2X"]
nod.col2[2,]=hgrsort2$colors[hgrsort2$Ploidy=="4X"]

col3 <- colorRampPalette(brewer.pal(9,"YlOrRd"))(50)
levelplot(nod.col2, col.regions= col3)
heatmap(nod.col2, Colv=NA, Rowv=NA, col=col3)

#===================================
#LINEAGE EFFECTS ON HGR
#===================================
hgr2 <- lm(HGRS.mean ~ Ploidy * Strain, data = p7)
plot(hgr2)
E_hgr2 = resid(hgr2, type = "deviance")
hist(E_hgr2)
qqnorm(E_hgr2)
qqline(E_hgr2)
summary(hgr2)
anova(hgr2)

#Compare models
anova(hgr,hgr2)

#Variation explained by accession and host species
Var_Random_effect <- VarCorr(hgr)
Var_Random_effect
variances <- as.data.frame(Var_Random_effect)
variances

#accession
variances$vcov[1]/sum(variances$vcov)

#host species
variances$vcov[2]/sum(variances$vcov)
#=========================================
#PLOIDY AVERAGES
#=========================================
p8 <- ddply(p6,c("Strain","Ploidy"),
            summarise,HGRS.mean=mean(HGR.SBM, na.rm = TRUE),
            Shoot.mean=mean(Shoot, na.rm=TRUE),
            Root.mean=mean(Shoot, na.rm=TRUE),
            QNC.mean=mean(QNodColor, na.rm=TRUE),
            QNC2.mean=mean(QNodColor2, na.rm=TRUE),
            NodNo.mean=mean(NodNo, na.rm=TRUE),
            NodBM.mean=mean(NodBMmg, na.rm=TRUE))
p8

hgrmeanc <- tapply(p8$HGRS.mean, p8$Ploidy, mean)
hgrmeanc
hgrmaxc <- tapply(p8$HGRS.mean, p8$Ploidy, max)
hgrmaxc
hgrminc <- tapply(p8$HGRS.mean, p8$Ploidy, min)
hgrminc

shootmeanc <- tapply(p8$Shoot.mean, p8$Ploidy, mean)
shootmeanc
shootmaxc <- tapply(p8$Shoot.mean, p8$Ploidy, max)
shootmaxc
shootminc <- tapply(p8$Shoot.mean, p8$Ploidy, min)
shootminc

#=========================================
#SPATIAL AND TEMPORTAL EFFECTS ON HGR 
#=========================================
p8<-ddply(p6,c("Strain", "HostSpecies","Acc","Ploidy", "BoxNo", "Round"),
          summarise,HGRS.mean=mean(HGR.SBM, na.rm = TRUE))
head(p8)

st9 <- lm(p8$HGRS.mean ~ p8$Ploidy * p8$Strain + p8$BoxNo)
summary(st9)
anova(st9)

st10 <- lm(p8$HGRS.mean ~ p8$Ploidy * p8$Strain + p8$Round)
summary(st10)
anova(st10)

#===================================
#NODULE TRAIT ANALYSES
#===================================
#Create dataset of 17 Nodulating Strains 
#Remove plants inoculated with non-nodulating strains from dataframe
p9=p7[p7$Strain!="USDA4894"&p7$Strain!="USDA205"&p7$Strain!="USDA207"&p7$Strain!="M30",]
head(p9)

#Trait Correlations
corfun<-function(x, y) {
  corr=(cor.test(x, y,
                 alternative="two.sided", method="pearson"))
}

#Nod Number and Nod Biomass
ddply(p9, .(Ploidy), summarise,z=corfun(NodBM.mean,NodNo.mean)$statistic,
      pval=corfun(NodBM.mean,NodNo.mean)$p.value,
      r.est=corfun(NodBM.mean,NodNo.mean)$estimate,
      alt=corfun(NodBM.mean,NodNo.mean)$alternative
) 

#Nod Color and Nod Number
ddply(p9, .(Ploidy), summarise,z=corfun(NodNo.mean,QNC.mean)$statistic,
      pval=corfun(NodNo.mean,QNC.mean)$p.value,
      r.est=corfun(NodNo.mean,QNC.mean)$estimate,
      alt=corfun(NodNo.mean,QNC.mean)$alternative
) 

#Nodule Number
nodno <- lmer(NodNo.mean ~ Ploidy * Strain + (1|HostSpecies/Acc), data = p9)
plot(nodno)
E_nodno = resid(nodno, type = "deviance")
hist(E_nodno)
qqnorm(E_nodno)
qqline(E_nodno)
summary(nodno)
anova(nodno)
nodnom <- emmeans(nodno, specs = pairwise ~ Ploidy:Strain)
nodnom

#Nodule Biomass
nodbm <- lmer(NodBM.mean ~ Ploidy * Strain + (1|HostSpecies/Acc), data = p9)
plot(nodbm)
E_nodbm = resid(nodbm, type = "deviance")
hist(E_nodbm)
qqnorm(E_nodbm)
qqline(E_nodbm)
summary(nodbm)
anova(nodbm)
nodbmm <- emmeans(nodbm, specs = pairwise ~ Ploidy:Strain)
nodbmm

nodbmmean <- tapply(p9$NodBM.mean, p9$Ploidy, mean)
nodbmmean

#Nodule Biomass + Root Biomass Covariate
nodbmr <- lmer(NodBM.mean ~ Ploidy * Strain + Root.mean + (1|HostSpecies/Acc), data = p9)
plot(nodbmr)
E_nodbmr = resid(nodbmr, type = "deviance")
hist(E_nodbmr)
qqnorm(E_nodbmr)
qqline(E_nodbmr)
summary(nodbmr)
anova(nodbmr)

#Nodule Number + Root Biomass Covariate
nodnor <- lmer(NodNo.mean ~ Ploidy * Strain + Root.mean + (1|HostSpecies/Acc), data = p9)
plot(nodnor)
E_nodnor = resid(nodnor, type = "deviance")
hist(E_nodnor)
qqnorm(E_nodnor)
qqline(E_nodnor)
summary(nodnor)
anova(nodnor)

#Nodule Color
qnc <- lmer(QNC.mean ~ Ploidy * Strain + (1|HostSpecies/Acc), data = p9)
plot(qnc)
E_qnc = resid(qnc, type = "deviance")
hist(E_qnc)
qqnorm(E_qnc)
qqline(E_qnc)
summary(qnc)
anova(qnc)
qncm <- emmeans(qnc, specs = pairwise ~ Ploidy:Strain)
qncm

#Nodule Color 2
qnc2 <- lmer(QNC2.mean ~ Ploidy * Strain + (1|HostSpecies/Acc), data = p9)
plot(qnc2)
E_qnc2 = resid(qnc2, type = "deviance")
hist(E_qnc2)
qqnorm(E_qnc2)
qqline(E_qnc2)
summary(qnc2)
anova(qnc2)

#3D Nodule Trait plot by ploidy
n1 <- as.data.frame(nodnom$emmeans) #NodNumber
n1
n2 <- as.data.frame(nodbmm$emmeans) #NodBM
n2
n3 <- as.data.frame(qncm$emmeans) #NodColor
n3

scatter3D(n1$emmean, n3$emmean, n2$emmean,
          bty = "u", colvar = as.integer(n1$Ploidy), pch=16, cex=1.5,
          colkey = list(at = c(1,2), labels= c("2X","4X")),col = c("gray50", "gray20"),
          ticktype = "detailed", xlab= "Nodule Number",
          ylab="Nodule Color",zlab="Total Nodule Biomass", type="h")

#=============================================
#PLANT PERFORMANCE REGRESSIONS - SHOOT BIOMASS
#=============================================
corfun<-function(x, y) {
  corr=(cor.test(x, y,
                 alternative="two.sided", method="pearson"))
}

#Nodule Biomass and Shoot Biomass
ddply(p9, .(Ploidy), summarise,z=corfun(Shoot.mean,NodBM.mean)$statistic,
      pval=corfun(Shoot.mean,NodBM.mean)$p.value,
      r.est=corfun(Shoot.mean,NodBM.mean)$estimate,
      alt=corfun(Shoot.mean,NodBM.mean)$alternative
) 

r.test(84, 0.6688938, 0.7039179)

#Nodule Color and Shoot Biomass
ddply(p9, .(Ploidy), summarise,z=corfun(Shoot.mean,QNC.mean)$statistic,
      pval=corfun(Shoot.mean,QNC.mean)$p.value,
      r.est=corfun(Shoot.mean,QNC.mean)$estimate,
      alt=corfun(Shoot.mean,QNC.mean)$alternative
) 

r.test(84, 0.5673017, 0.5841607)

#Correlation Plot
c1 <- ggscatter(p9, x = "NodBM.mean", y = "Shoot.mean", color = "Ploidy", 
                conf.int = F, 
                cor.coef = FALSE, cor.method = "pearson", palette = c("gray70", "gray20"),
                xlab = "Total Nodule Biomass", ylab = "Shoot Biomass")
c1 + stat_cor(aes(color = Ploidy), label.x=0.25) 


#=============================================
#PLANT PERFORMANCE REGRESSIONS - HGR
#=============================================

#Nodule Biomass and HGR
ddply(p9, .(Ploidy), summarise,z=corfun(HGRS.mean,NodBM.mean)$statistic,
      pval=corfun(HGRS.mean,NodBM.mean)$p.value,
      r.est=corfun(HGRS.mean,NodBM.mean)$estimate,
      alt=corfun(HGRS.mean,NodBM.mean)$alternative
) 

r.test(84, 0.4200151, 0.6417595)

#Nodule Color and HGR
ddply(p9, .(Ploidy), summarise,z=corfun(HGRS.mean,QNC.mean)$statistic,
      pval=corfun(HGRS.mean,QNC.mean)$p.value,
      r.est=corfun(HGRS.mean,QNC.mean)$estimate,
      alt=corfun(HGRS.mean,QNC.mean)$alternative
) 

r.test(84,  0.5494668, 0.6073102)

#Correlation Plot
c3 <- ggscatter(p9, x = "HGRS.mean", y = "NodBM.mean", color = "Ploidy", 
                conf.int = F, 
                cor.coef = FALSE, cor.method = "pearson", palette = c("gray70", "gray20"),
                xlab = "Host Growth Response - Shoot", ylab = "Total Nodule Biomass")
c3 + stat_cor(aes(color = Ploidy), label.x=0.25) 

#Correlation Plot
c4 <- ggscatter(p9, x = "NodBM.mean", y = "HGRS.mean", color = "Ploidy", 
                conf.int = F, 
                cor.coef = FALSE, cor.method = "pearson", palette = c("gray70", "gray20"),
                xlab = "Total Nodule Biomass", ylab = "Host Growth Response - Shoot")
c4 + stat_cor(aes(color = Ploidy), label.x=0.25) 

#======================================================
#SHOOT AND HGR ANALYSES WITH ONLY 17 NODULATING STRAINS 
#======================================================
shoot2 <- lmer(Shoot.mean ~ Ploidy * Strain + (1|HostSpecies/Acc), data = p9)
plot(shoot2)
E_shoot2 = resid(shoot2, type = "deviance")
hist(E_shoot2)
qqnorm(E_shoot2)
qqline(E_shoot2)
summary(shoot2)
anova(shoot2)


hgr2 <- lmer(HGRS.mean ~ Ploidy * Strain + (1|HostSpecies/Acc), data = p9)
plot(hgr2)
E_hgr2 = resid(hgr2, type = "deviance")
hist(E_hgr2)
qqnorm(E_hgr2)
qqline(E_hgr2)
summary(hgr2)
anova(hgr2)

#======================================================
#RDPI ANALYSES - HGR
#======================================================
#RDPI Analysis including all 21 strains
#Create list separating by accesion 
L=split(p7,p7$Acc)
L

#Function to calculate RDPI of HGR
rdpi=function(x){
  m=as.matrix(x$HGRS.mean)
  dist=vegdist(m,method="canberra")
  sum(dist)/length(dist)
}

#Make data frame for plot 
rdpi.data=data.frame(rdpi=unlist(lapply(L, rdpi)),ploidy=c("4X","2X","2X","2X","4X","4X","4X","2X","4X","2X"))
rdpi.data

#t.test
t.test(rdpi.data$rdpi[rdpi.data$ploidy=="2X"], 
       rdpi.data$rdpi[rdpi.data$ploidy=="4X"])

#RDPI range 
rdpir <- tapply(rdpi.data$rdpi, rdpi.data$ploidy, range)
rdpir

#RDPI variance
rdpiv <- tapply(rdpi.data$rdpi, rdpi.data$ploidy, var)
rdpiv

#RDPI plot
ggplot(rdpi.data,aes(x=ploidy,y=rdpi,fill=ploidy))+
  geom_point(shape=21)+
  stat_summary(fun.y = mean,geom = "point",shape=21,size=3)+
  theme_bw()

#--------
#RDPI - 17 Nodulating Strains
#Create list separating by accesion 
L2=split(p9,p9$Acc)
L2

#Make data frame for plot 
rdpi.data2=data.frame(rdpi=unlist(lapply(L2, rdpi)),ploidy=c("4X","2X","2X","2X","4X","4X","4X","2X","4X","2X"))
rdpi.data2

#t.test
t.test(rdpi.data2$rdpi[rdpi.data2$ploidy=="2X"], 
       rdpi.data2$rdpi[rdpi.data2$ploidy=="4X"])

#======================================
#COST OF SPECIALIZATION ANALYSES - HGR
#======================================
#21 Strains
#Calculate cost of specialization
cost=function(x){
  max_val=max(x$HGRS.mean)
  dist=max_val-x$HGRS.mean
  dist=dist[dist!=0]/max_val
  sum(dist)/length(dist)
}

#Make data frame for plot
cost.data=data.frame(cost=unlist(lapply(L, cost)),ploidy=c("4X","2X","2X","2X","4X","4X","4X","2X","4X","2X"))
cost.data

#t.test
t.test(cost.data$cost[cost.data$ploidy=="2X"], 
       cost.data$cost[cost.data$ploidy=="4X"])

#cost range 
costr <- tapply(cost.data$cost, cost.data$ploidy, range)
costr

#cost variance
costv <- tapply(cost.data$cost, cost.data$ploidy, var)
costv

#cost plot
ggplot(cost.data,aes(x=ploidy,y=cost,fill=ploidy))+
  geom_point(shape=21)+
  stat_summary(fun.y = mean,geom = "point",shape=21,size=3)+
  theme_bw()

#---------------
#17 Strains

#Make data frame for plot
cost.data2=data.frame(cost=unlist(lapply(L2, cost)),ploidy=c("4X","2X","2X","2X","4X","4X","4X","2X","4X","2X"))
cost.data2

#t.test
t.test(cost.data2$cost[cost.data2$ploidy=="2X"], 
       cost.data2$cost[cost.data2$ploidy=="4X"])

#--------------------------
#Correlation RDPI and Cost - 21 strains
cor.test(rdpi.data$rdpi, cost.data$cost)
