#### Example script that can be used to analyze GI-trade-off in Arabidopsis competition experiment. 
### Code by Samuel Wuest, with contributions from Pascal Niklaus 

######################################################################
### define functions used to calculate G-I trade off
### distance to fitted model
###
### called only by distTofit

.dfun <- function(x, ...) {
  a <- list(...)
  unname(sqrt((a$x0-x)^2 + (a$y0-predict(a$model, newdata = data.frame(indPer=x)))^2))
}

######################################################################
###
### distToFit function
###
### arguments:
###   (x0, y0)   coordinates of point
###   model      fitted regression curve, must have parameter 'indPer'
###
### result (list):
###   (x, y)     closest point on regression
###   d          distance (x0,y0)--(x,y)

distToFit <- function(x0, y0, model) {
  p <- x0+0.5*(y0-x0)  # <- better for near 1:1 models
  p <- x0
  for(i in 1:20) {
    ff <- nlm(.dfun, p=p, model=model, x0=x0,y0=y0)
    if(ff$code != 4)
      break;
    p <- ff$estimate
  }
  if(ff$code==4)
    stop("Did not converge, iteration limit exceeded!")
  x1 <- ff$estimate
  y1 <- unname(predict(model, newdata = data.frame(indPer=x1)))
  list(x=x1, y=y1, d=.dfun(x1,x0=x0,y0=y0,model=model))
}

### packages used: lsmeans; lme4; lmerTest (caution: lsmeans-Function both in lsmeans and lmerTest packages)
require(lsmeans)
require(lme4)
require(lmerTest)

## load data
d <- read.csv("competition.csv")
str(d)

##################################################################
######################### some pre-processing:

### remove pots in which only one individual survived
d2 <- subset(d, !is.na(Realized_PosB)&!is.na(Realized_PosA))
## accessionID has to be coded as factor
d2$AccID_Tester_posA <- as.factor(d2$AccID_Tester_posA); d2$AccID_Target_posB <- as.factor(d2$AccID_Target_posB)
### add genotype composition
d2$GenoComp <- paste(d2$AccID_Tester_posA, d2$AccID_Target_posB, sep=" ")
##################################################################

##################################################################
### now start calculating necessary coefficients for the analysis:

######################
### a first model to fit effects of target genotype on target biomass across all mixtures, and extract individual performances
mod1 <- lm(Biomass_mg_posB~Block+AccID_Target_posB,data=subset(d2, CommunityType=="mixture"))
## extract the coefficients:
### indPer: mean individual biomass across all mixtures
indPer <- as.data.frame(summary(lsmeans::lsmeans(mod1, specs=~AccID_Target_posB)))[,1:2]

######################
### also estimate tester individual biomass across all mixtures:
mod2 <- lm(Biomass_mg_posA~Block+AccID_Tester_posA,data=subset(d2, CommunityType=="mixture"))
testerPer <- as.data.frame(summary(lsmeans::lsmeans(mod2, specs=~AccID_Tester_posA)))[,1:2]
d2$meanTesterPer <- testerPer$lsmean[match(d2$AccID_Tester_posA, testerPer$AccID_Tester_posA)]
### effect of tester genotypic mean on neighbor biomass:
anova(lmer(Biomass_mg_posB~Block+meanTesterPer+(1|GenoComp), data=subset(d2,CommunityType=="mixture")))
# Analysis of Variance Table of type III  with  Satterthwaite approximation for degrees of freedom
# Sum Sq Mean Sq NumDF  DenDF F.value    Pr(>F)    
# Block         806496  806496     1 917.54  97.268 < 2.2e-16 ***
#  meanTesterPer 731518  731518     1 960.03  88.225 < 2.2e-16 ***


######################
### model estimating the performance of a given genotype composition: extract only monoculture performances
genoCompMod <- lm(PotBiomass_mg~Block+GenoComp, data=d2)
### extract all coefficients:
commMeans <- as.data.frame(summary(lsmeans::lsmeans(genoCompMod, specs=~GenoComp)))[,1:2]
commMeans$geno1 <- substr(commMeans$GenoComp,1,4)
commMeans$geno2 <-   substr(commMeans$GenoComp,6,10)
### extract monoculture performance (monoPer) estimates:
monoPer <- subset(commMeans, geno1==geno2)
######################

### put everything into a "mydat"-data frame that contains genotypic values for all 98 accessions 
mydat1 <- merge(indPer, monoPer, by.x="AccID_Target_posB", by.y="geno1", all.x=T)[,c(1,2,4)]
colnames(mydat1) <- c("accID", "indPer", "monoPer")
head(mydat1)

######################
### also add the monoculture washable rootmass: 
genoCompMod2 <- lm(RootMass_mg_corrected~Block+GenoComp, data=d2)
commMeans2 <- as.data.frame(summary(lsmeans::lsmeans(genoCompMod2, specs=~GenoComp)))[,1:2]
commMeans2$geno1 <- substr(commMeans2$GenoComp,1,4)
commMeans2$geno2 <-   substr(commMeans2$GenoComp,6,10)
monoRoot <- subset(commMeans2, geno1==geno2); colnames(monoRoot)[2] <- "monoRoot"; monoRoot <- monoRoot[,-c(1,4)]

#### merge all estimates into a single data frame (mydat)
mydat <- merge(mydat1, monoRoot, by.x="accID", by.y="geno1", all.x=T)
######################
### shoot to root ration in monocultures:
mydat$monoRootShoot <- mydat$monoRoot/mydat$monoPer

######################
### add info on allele at SNP Chr 3, position 15294955; this is done manually here to avoid loading the whole SNP-dataset
allele <- as.data.frame(matrix(data=c("5837" , "6009" , "6043" , "6046" , "6897" , "6898" , "6899" , "6900" , "6901" , "6903" , "6904" , "6905" , "6906" , "6907" , "6908" , "6909" , "6910" , "6911" , "6913" , "6914" , "6915" , "6916" , "6917" , "6918" , "6919" , "6920" , "6921" , "6922" , "6926" , "6927" , "6928" , "6929" , "6931" , "6932" , "6933" , "6936" , "6937" , "6939" , "6940" , "6942" , "6943" , "6944" , "6945" , "6946" , "6951" , "6956" , "6958" , "6959" , "6960" , "6961" , "6962" , "6963" , "6965" , "6966" , "6968" , "6969" , "6970" , "6971" , "6972" , "6973" , "6974" , "6975" , "6976" , "6977" , "6978" , "6979" , "6980" , "6981" , "6982" , "6983" , "6984" , "7186" , "7275" , "7328" , "7378" , "7413" , "7514" , "7515" , "7516" , "7517" , "7518" , "7519" , "7520" , "7521" , "7523" , "7524" , "7525" , "7526" , "8213" , "8214" , "8215" , "8264" , "8274" , "8297" , "8304" , "8387" , "8412" , "8424" , "A" , "A" , "A" , "A" , "A" , "A" , "C" , "A" , "A" , "C" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "C" , "C" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "C" , "A" , "A" , "A" , "C" , "C" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "C" , "A" , "A" , "C" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "C" , "A" , "A" , "A" , "C" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "A" , "C" , "C"), ncol=2, byrow=F))
mydat$geno <- allele[match(mydat$accID, allele[,1]),2]

### also add genotypic rosette diameter estimate across all compositions
mod1 <- lm(Rosette_PosB_22to24March~Block+AccID_Target_posB,data=d2)
tmp <- summary(lsmeans::lsmeans(mod1, ~AccID_Target_posB))
mydat$Rosette <- tmp$lsmean[match(mydat$accID, tmp$AccID_Target_posB)]



############################################################### 
###### now calculate G-I trade-off and run some tests:
###############################################################

### fit a quadratic model explaining mean monoculture performance 
d.lm <- lm(monoPer ~ indPer + I(indPer^2), data=mydat)
anova(d.lm)

### plot with alleles at cooperation-associated locus indicated by color; remove asp=1 if nicer plot required
plot(monoPer~indPer, data=mydat, pch=16, col=geno, cex=1.3, xlim=c(100, 800), ylim=c(200, 1500), asp=1)
predicted.intervals <- predict(d.lm,data.frame(indPer=seq(50,900,5)),interval='confidence',level=0.99)
lines(seq(50,900,5), predicted.intervals[,2],col='grey',lwd=1)
lines(seq(50,900,5), predicted.intervals[,1],col='black',lwd=1)
lines(seq(50,900,5), predicted.intervals[,3],col='grey',lwd=1)

### add an arbitrary point and check if G-I trade-off is correctly derived

x0 <- 500
y0 <- 1500

### calculate nearest point on curve and distance

r <- distToFit(x0, y0, d.lm)

points(c(x0,r$x),c(y0,r$y),pch=16,col="blue")

lines(c(x0,r$x),c(y0,r$y),pch=16,col="blue")

### calculate the distance to the regression line; points above the line will have positive distances, points below the line will have negative distances
### we will use monoculture performance per individual plant (i.e. divide monoculture values by two)
tmpMat <- mydat[,colnames(mydat)%in%c("indPer", "monoPer")]; tmpMat$monoPer <- 0.5*tmpMat$monoPer
d.lm2 <- lm(monoPer ~ indPer + I(indPer^2), data=tmpMat)
mydat$tradeoff  <- sign(resid(d.lm2))*apply(tmpMat, 1, function(m) return(distToFit(x0=m[1], y0=m[2], d.lm2)$d))
hist(mydat$tradeoff)

### one could also argue that the residual is the interesting trait (more productive in monoculture than expected based on indiv. perf. in mixture)
mydat$resid <- resid(d.lm2); plot(resid~tradeoff, data=mydat, pch=16, col=geno)


##### some simple models explaining G-I tradeoff 

### association with allele at Chr3:15294955
### please note: this is w/o adjustment for pop-structure; however, results correspond quite well with EMMAX-results from easygwas
anova(lm(tradeoff~geno, data=mydat))
summary(lm(tradeoff~geno, data=mydat))

### compare this to:
anova(lm(resid~geno, data=mydat))

#### write datat into a file for GWAS (easygwas.ethz.ch)
out <-mydat[,c(1,1,2:ncol(mydat))]; colnames(out)[1:2] <- c("FID", "IID"); out <- out[,-grep("geno", colnames(out))]
#write.table(out, file="phenotype", sep="\\t", row.names=F)


#### trait-based analyses
anova(lm(tradeoff~Rosette, data=mydat))
anova(lm(tradeoff~monoRootShoot, data=mydat))
anova(lm(monoRootShoot~geno, data=mydat))


plot(tradeoff~monoRootShoot, data=mydat, pch=16)
abline(lm(tradeoff~monoRootShoot, data=mydat))


noco2 <- c("r", "r", "s", NA,"r", "s", "i", "r", "s", "r", "r", "r", "r", "s", "r", "s", "r", "s", "r", "s",NA, "r", "r", "s", "s", "s", "s", "r", "s", "s", "s", "s", "s", "r", "s", "r", "s", "r", "r", "s","s", "r", "s", "r", "r", NA, "i", "r", "r", "s","r", "r", NA,"r", "r", "s", "s", "s", "s", "r","r","s", "r", "r", "r", "s", "r", "r", "r", "s", "r", NA,NA,NA,NA,NA,"r", "s", "s", "r","s", "s", "r", "r", "r", NA,"s", "s", "s", "r", "s", NA,NA,NA,NA,NA,NA,"i")
names(noco2) <- c("5837" , "6009" , "6043" , "6046" , "6897" , "6898" , "6899" , "6900" , "6901" , "6903" , "6904" , "6905" , "6906" , "6907" , "6908" , "6909" , "6910" , "6911" , "6913" , "6914" , "6915" , "6916" , "6917" , "6918" , "6919" , "6920" , "6921" , "6922" , "6926" , "6927" , "6928" , "6929" , "6931" , "6932" , "6933" , "6936" , "6937" , "6939" , "6940" , "6942" , "6943" , "6944" , "6945" , "6946" , "6951" , "6956" , "6958" , "6959" , "6960" , "6961" , "6962" , "6963" , "6965" , "6966" , "6968" , "6969" , "6970" , "6971" , "6972" , "6973" , "6974" , "6975" , "6976" , "6977" , "6978" , "6979" , "6980" , "6981" , "6982" , "6983" , "6984" , "7186" , "7275" , "7328" , "7378" , "7413" , "7514" , "7515" , "7516" , "7517" , "7518" , "7519" , "7520" , "7521" , "7523" , "7524" , "7525" , "7526" , "8213" , "8214" , "8215" , "8264" , "8274" , "8297" , "8304" , "8387" , "8412" , "8424" )

mydat$Noco2 <- noco2[match(mydat$accID, names(noco2))]
mydat$Noco2TwoClasses <- ifelse(mydat$Noco2=="i", "r", ifelse(mydat$Noco2=="s", "anf", mydat$Noco2))

anova(lm(tradeoff~Noco2, data=mydat))
anova(lm(tradeoff~Noco2TwoClasses, data=mydat))

fisher.test(table(mydat$Noco2, mydat$geno))
barplot(apply(table(mydat$Noco2, mydat$geno)[c(2,1,3),], 2, function(x) x/sum(x)), legend.text = levels(mydat$Noco2), col=c("blue", "lightblue2", "red4"))


### as an example, plot G-I trade-off value ~ allele 

mydat$col <- ifelse(mydat$geno=="A", "blue", "dark red")
### s.e.m for each group
myses <- tapply(mydat$tradeoff, mydat$geno, function(x) sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)])))
#### mean for each group
mymeans <- tapply(mydat$tradeoff, mydat$geno, mean, na.rm=T)

### plot
plot(main="",c(0.5,1.5), mymeans, type="n", xlim=c(0.2, 1.8), xlab="allele", ylab="Trade-off (mg/mg)", ylim=c(min(mydat$tradeoff, na.rm=T),max(mydat$tradeoff, na.rm=T)), xaxt="n", cex=2, las=1)
segments(x0 = c(0.5,1.5), mymeans-myses, c(0.5,1.5), mymeans+myses, lwd=4, col="dark red")
points(jitter(as.integer(as.factor(mydat$geno))-0.5), mydat$tradeoff, pch=16, col=mydat$col,cex=1.6)
segments(x0=0.5, y0=mymeans[1], x1=1.5, y1=mymeans[2], lwd=3, col="dark red")
axis(side=1,at=c(0.5, 1.5), labels=levels(mydat$geno))

