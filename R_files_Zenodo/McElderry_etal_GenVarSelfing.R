################################################################################
##                Incipient Evolution of Selfing Syndrome                     ##
##                 University of Tennessee, Knoxville                         ##
##                       Data from 1999 - 2002                                ##
## Created: October 2018                                                      ##
##      By: Robert McElderry | mcelderone@gmail.com                           ##
################################################################################
library(plotrix)
library(MCMCglmm)
library(coda)
# __| Overview |____________________________________________________________----
# A study of floral morphology conducted on Collinsia verna using 
# Optimas Image Analysis. 
# 
# Genotypes from 3 populations, Enlow Fork, Braddock Trail, and Ten Mile Creek, 
# all in Western Pennsylvania, were cloned and clones were designated to either 
# this morphology study, the autogamy study, or the inbreeding depression study.  
# 
# Genotypes 100-199 are from Enlow Fork, 200-299 are from Braddock Trail, and 
# 300-399 are from Ten Mile Creek.
# 
# At least 2 flowers from each stage were collected from each clone and stored 
# in ethanol for floral morphology measurements.Stage # = the number of anthers 
# dehisced.  Stage 0 = 0 anthers dehisced and Stage 4 = all four anthers 
# dehisced.
# 
# There was a separate macro and calibrations for Parts A, B, C, and D.  For 
# an explanantion of of how each measurement was taken, see the Image Analysis 
# instructions in C:datasheets and protocols/Protocol Sheets for Image Analysis
# A "." (RM converted to "NA") indicates no data was taken.
#_______________________________________________________________________________ 

# __| COLUMN DESCRIPTION |______________________________________________________
#         GENOTYPE  =	100-199 -> BT, 200-299 -> EF, 300-399 -> TMC
#             CLONE =	Clone number 
#             STAGE =	Floral stage, number of anthers dehisced (0-4)
#         REPLICATE =	Replicate number (1 to 3 flowers per clone)
#           SEPAL_L	=	length of 1 sepal
#       UPPERPTL_L	=	length of upper petal - trace and not a straight line
#       LOWERPTL_L	=	length of lower petal - trace and not a straight line
#           KEEL_L	=	Keel length - trace and not a straight line
#           KEEL_A	=	Keel areas was calculated by tracing around the 
#                     folded keel (so total petal area is acually about 
#                     twice this if you consider both sides of the fold) 
#         STG_ANT_L	=	distance between stigma and nearest dehisced anther -- 
#                     measured in macro -- adjusted where necessary for 
#                     incorrect calibrations.  In old protocol data (see old 
#                     protocol column) the anthers were manipulated before 
#                     this measurement was taken, so measurements are not 
#                     very accurate.
#   minimum STG-ANT =	this is an estimation of distance from stigma to 
#                     nearest dehisced anther.  It is the smallest 
#                     difference between the length of any dehisced filament 
#                     and the length of the style.
#             ZONE	=	(+) style above nearest dehisced anther, (0) stigma at 
#                     or very near nearest dehisced anther, (-) stigma 
#                     below nearest dehisced anther
#           OVARY_W	=	ovary width measured again at higher magnification, 
#                     this (and all below) were adjusted for some data that 
#                     was taken and calibrated incorrectly
#           FILM1_L	=	linear length of filament 1 (1st anther to dehisce).
#       FILM1CURV_L	=	trace of filament 1.  
#           FILM2_L	=	linear length of filament 2 (2nd anther to dehisce)
#       FILM2CURV_L	=	trace of filament 2.  
#           FILM3_L	=	linear length of filament 3 (3rd anther to dehisce).
#       FILM3CURV_L	=	trace of filament 3.  
#           FILM4_L	=	linear length of filament 4 (4th anther to dehisce)
#       FILM4CURV_L	=	trace of filament 4.
#         ANTHER_A	=	anther area of the next anther to dehisce.
#           STYLE_L	=	linear length of style.
#           STYLE_C	=	trace of style.
# old protocol set	=	a "y" in this column indicates that the flower  was 
#                     in the first batch of image analysis done.  The keel 
#                     was cut and filaments were adjusted so you could see 
#                     them all before stigma-anther distance was measured, 
#                     so stigma-anther distance is not as accurate for this 
#                     set of data.
#___________________________________________________________________________---- 
#_______________________________________________________________________________
# -- | 1. Floral trait data |---------------------------------------------------
rm(list=ls(all=T))  ##Clears R's short term memory
setwd("C:/Users/mceld/Dropbox/UTK_2018-19/Manuscripts/Spigler_etal_GenVarSelfing_2019/Data Analysis")
dat <- read.csv("Data/1_Spigler_etal_CollinsiaFloralTraits.csv", header=TRUE,stringsAsFactors=T)
summary(dat)
dim(dat) # 4813   25
# Modifications
dat$old.protocol.set <- as.character(dat$old.protocol.set)
dat$old.protocol.set[is.na(dat$old.protocol.set)==T] <- 'n'
dat$old.protocol.set <- as.factor(dat$old.protocol.set)
dat$POPULATION <- factor(NA,levels=1:3,labels=c('BT','EF','TMC'))
dat$POPULATION[dat$GENOTYPE>=100 & dat$GENOTYPE<200] <- 'EF'
dat$POPULATION[dat$GENOTYPE>=200 & dat$GENOTYPE<300] <- 'BT'
dat$POPULATION[dat$GENOTYPE>=300 & dat$GENOTYPE<400] <- 'TMC'
dat[!is.na(dat$FILM1_L) & dat$FILM1_L<=2,]
dat$plantID <- paste(dat$GENOTYPE,'_',dat$CLONE,sep='') # Helps link other datasets
# Caveat | There are some weird numbers in this data set. Filament lenghts of 
# 0.02 are just not reasonable. The data sorely need cleaning. I simply shift 
# decimal points below to put a few petal measurements in the right order of 
# magnitude, but these are just a few of the signs I have seen indicating 
# consistent unchecked issues.

# All filament and style lengths were measured for each stage.
# we are interested only in open anthers. So only compute SAD for those 
# filaments with dehisced antherw, i.e., SAD1 = distance to filament1.
# SPZD | minimum distance to the pollen zone (region between all dehisced anthers). 
# |SPZD| = minimum STG-ANT except in the pollen zone where SPZD = 0 
# mSAD | mean distance considering all dehisced anthers. 
# all of these metrics have issues, but each contributes an important piece. 
# What if we combined SPZD and mSAD in different stages/regions? SPZD, the 
# minimum distance to pollen zone seems appropriate until the stigma enters the 
# pollen zone. Then mSAD may make the most sense, representing the mean distance 
# to all dehisced anthers. Here mean distance may arguably represent probability 
# of contact or intensity of contact. Outside of the pollen zone its all about 
# getting there and less about whether or not two other dehisced anthers are 
# further away. a fourth (really fifth) composite variable may be needed. 

dat <- subset(dat,STAGE!=0) # Without STAGE, data are of no use.
dat$SPZD <- dat$mSAD <- dat$SAD1 <- dat$SAD2 <- dat$SAD3 <- dat$SAD4 <- NA
dat$CONTACT <- dat$PolZone <- dat$PolIn <- dat$PolOut <- NA 
for (i in 1:nrow(dat)) {
  if (dat$STAGE[i]==1) {
    dat$SAD1[i] <- dat$STYLE_L[i] - dat$ FILM1_L[i]
    dat$PolOut[i] <- max(dat[i,c('FILM1_L')]) 
    dat$PolIn[i] <- min(dat[i,c('FILM1_L')])
  }
  if (dat$STAGE[i]==2) {
    dat$SAD1[i] <- dat$STYLE_L[i] - dat$ FILM1_L[i]
    dat$SAD2[i] <- dat$STYLE_L[i] - dat$ FILM2_L[i]
    dat$PolOut[i] <- max(dat[i,c('FILM1_L','FILM2_L')]) 
    dat$PolIn[i] <- min(dat[i,c('FILM1_L','FILM2_L')])
  }
  if (dat$STAGE[i]==3) {
    dat$SAD1[i] <- dat$STYLE_L[i] - dat$ FILM1_L[i]
    dat$SAD2[i] <- dat$STYLE_L[i] - dat$ FILM2_L[i]
    dat$SAD3[i] <- dat$STYLE_L[i] - dat$ FILM3_L[i]
    dat$PolOut[i] <- max(dat[i,c('FILM1_L','FILM2_L','FILM3_L')])
    dat$PolIn[i] <- min(dat[i,c('FILM1_L','FILM2_L','FILM3_L')])
  }
  if (dat$STAGE[i]==4) {
    dat$SAD1[i] <- dat$STYLE_L[i] - dat$ FILM1_L[i]
    dat$SAD2[i] <- dat$STYLE_L[i] - dat$ FILM2_L[i]
    dat$SAD3[i] <- dat$STYLE_L[i] - dat$ FILM3_L[i]
    dat$SAD4[i] <- dat$STYLE_L[i] - dat$ FILM4_L[i]
    dat$PolOut[i] <- max(dat[i,c('FILM1_L','FILM2_L','FILM3_L','FILM4_L')])
    dat$PolIn[i] <- min(dat[i,c('FILM1_L','FILM2_L','FILM3_L','FILM4_L')])
  }
  dat$PolZone[i] <- dat$PolOut[i] - dat$PolIn[i]
  SADs <- c(dat$SAD1[i],dat$SAD2[i],dat$SAD3[i],dat$SAD4[i])
  n <- sum((is.na(SADs)==F)*1)
  if (n==0) dat$SPZD[i] <- NA
  else  { dat$SPZD[i] <- SADs[which.min(abs(SADs))]; dat$mSAD[i] <- mean(SADs,na.rm=T) }
  if(sum((SADs<0)*1,na.rm=T)==n) dat$CONTACT[i] <- 1 # inside of flower
  if(sum((SADs>0)*1,na.rm=T)==n) dat$CONTACT[i] <- 3 # outside of flower
  if(sum((SADs<0)*1,na.rm=T)<n & sum((SADs>0)*1,na.rm=T)<n) {
    dat$CONTACT[i] <- 2 # pollen zone
    dat$SPZD[i] <- 0 # SAD is zero if in the pollen zone.
  }
  rm(SADs)
}


# CONTACT | Indicate whether the style is inside of the pollen zone (shorter 
# than the inner pollen zone boundary), within the pollen zone, or outside of 
#  the pollen zone. 
dat$CONTACT <- factor(dat$CONTACT,levels=1:3,labels=c('in','pz','out'))
# SAD | SPZD when shorter or longer than the pollen zone; mSAD when within the 
# pollen zone.
dat$SAD <- dat$SPZD
dat$SAD[dat$CONTACT=='pz'] <- dat$mSAD[dat$CONTACT=='pz']   # Use SAD!!!!!!!!!!!
nrow(dat[is.infinite(dat$SAD)==T,]) # Check | Should equal zero

x11(width=12, height=8); par(mfrow=c(2,4))
plot(dat$SPZD,dat$SAD,xlab='SPZD',ylab='SAD')
plot(dat$SAD,dat$mSAD,xlab='SAD',ylab='mean SAD')
plot(dat$SAD,dat$minimum.STG.ANT,xlab='SAD',ylab='abs SPZD')
plot(dat$mSAD,dat$minimum.STG.ANT,xlab='mean SAD',ylab='abs SPZD')
boxplot(dat$SPZD~as.factor(dat$STAGE),xlab='Flower Stage',ylab='SPZD')
boxplot(dat$SAD~as.factor(dat$STAGE),xlab='Flower Stage',ylab='SAD')
boxplot(dat$mSAD~as.factor(dat$STAGE),xlab='Flower Stage',ylab='meanSAD')
boxplot(dat$minimum.STG.ANT ~as.factor(dat$STAGE),xlab='Flower Stage',ylab='abs SPZD')
# Recreate between stage correlations
# old SAD
gendat <- tapply(dat$minimum.STG.ANT,list(as.factor(dat$GENOTYPE),dat$STAGE),mean,na.rm=T)
st.cor <- array(NA,dim=c(4,3,3))
dimnames(st.cor) <- list(c('s1.s2','s2.s3','s3.s4','s1.s4'),c('est','lo','hi'),
                         c('old','SAD','mSAD'))
for (i in 1:3) {
  st.cor[i,1,1] <- cor.test(gendat[,i],gendat[,i+1],use='complete.obs')$estimate
  st.cor[i,2:3,1] <- cor.test(gendat[,i],gendat[,i+1],use='complete.obs')$conf.int[1:2]
}
st.cor[4,1,1] <- cor.test(gendat[,1],gendat[,4],use='complete.obs')$estimate
st.cor[4,2:3,1] <- cor.test(gendat[,1],gendat[,4],use='complete.obs')$conf.int[1:2]
# new SAD
gendat <- tapply(dat$SAD,list(as.factor(dat$GENOTYPE),dat$STAGE),mean,na.rm=T)
for (i in 1:3) {
  st.cor[i,1,2] <- cor.test(gendat[,i],gendat[,i+1],use='complete.obs')$estimate
  st.cor[i,2:3,2] <- cor.test(gendat[,i],gendat[,i+1],use='complete.obs')$conf.int[1:2]
}
st.cor[4,1,2] <- cor.test(gendat[,1],gendat[,4],use='complete.obs')$estimate
st.cor[4,2:3,2] <- cor.test(gendat[,1],gendat[,4],use='complete.obs')$conf.int[1:2]
# mean SAD
gendat <- tapply(dat$mSAD,list(as.factor(dat$GENOTYPE),dat$STAGE),mean,na.rm=T)
for (i in 1:3) {
  st.cor[i,1,3] <- cor.test(gendat[,i],gendat[,i+1],use='complete.obs')$estimate
  st.cor[i,2:3,3] <- cor.test(gendat[,i],gendat[,i+1],use='complete.obs')$conf.int[1:2]
}
st.cor[4,1,3] <- cor.test(gendat[,1],gendat[,4],use='complete.obs')$estimate
st.cor[4,2:3,3] <- cor.test(gendat[,1],gendat[,4],use='complete.obs')$conf.int[1:2]
st.cor # absolute value disrupts the correlation between stages. We want one stage to represent all stages. It matters if the stigma is outside or inside the flower. They should not be treated the same.

x11();
plot(1:3,st.cor[1:3,1,1],type='l',lty=1,col=1,lwd=2,xaxt='n',xlim=c(0.75,4.25),
     ylim=c(-0.1,1.0),xlab='Stage transition',ylab='Correlation')
axis(1,1:4,labels=dimnames(st.cor)[[1]])
lines(1:3,st.cor[1:3,1,2],lty=2,col=3,lwd=2)
lines(1:3,st.cor[1:3,1,3],lty=3,col=4,lwd=2)
plotCI(x=1:3,y=st.cor[1:3,1,1],liw=st.cor[1:3,1,1]-st.cor[1:3,2,1],
       uiw=st.cor[1:3,3,1]-st.cor[1:3,1,1],add=T) 
plotCI(x=1:3,y=st.cor[1:3,1,2],liw=st.cor[1:3,1,2]-st.cor[1:3,2,2],
       uiw=st.cor[1:3,3,2]-st.cor[1:3,1,2],add=T,col=3)
plotCI(x=1:3,y=st.cor[1:3,1,3],liw=st.cor[1:3,1,3]-st.cor[1:3,2,3],
       uiw=st.cor[1:3,3,3]-st.cor[1:3,1,3],add=T,col=4)
plotCI(x=c(4,4,4),y=st.cor[4,1,1:3],liw=st.cor[4,1,1:3]-st.cor[4,2,1:3],
       uiw=st.cor[4,3,1:3]-st.cor[4,1,1:3],add=T,col=c(1,3,4))
legend('topright',c('Old','SAD','mSAD'),col=c(1,3,4),lty=1:3,lwd=2,pch=1)
# Switch the sign on SAD. It's easier to discuss this way
dat$SAD <- -dat$SAD
dat$aSAD <- abs(dat$SAD)
#_______________________________________________________________________________

#-- | Petal Traits |-----------------------------------------------------------|
summary(dat[,c('UPPERPTL_L','KEEL_L')])
# There are a few outliers here that need attention.
dat[dat$UPPERPTL_L<5 & is.na(dat$UPPERPTL_L)==FALSE,c('UPPERPTL_L','LOWERPTL_L',
                                                      'KEEL_L','KEEL_A')]  # Def a decimal error
up <- dat[dat$UPPERPTL_L<5 & is.na(dat$UPPERPTL_L)==FALSE,c('UPPERPTL_L')]
up*10
dat[dat$UPPERPTL_L<5 & is.na(dat$UPPERPTL_L)==FALSE,c('UPPERPTL_L')] <- up*10
plot(dat$KEEL_L,dat$KEEL_A)
dat[dat$KEEL_L>14 & is.na(dat$KEEL_L)==FALSE,c('UPPERPTL_L','LOWERPTL_L',
                                               'KEEL_L','KEEL_A')]  # typo?
kl <- dat[dat$KEEL_L>14 & is.na(dat$KEEL_L)==FALSE,'KEEL_L']
kl-10
dat[dat$KEEL_L>14 & is.na(dat$KEEL_L)==FALSE,'KEEL_L'] <- kl-10
# better
summary(dat$ANTHER_A) # there are some unreal sizes here too
# set a cutoff of 2.5mm.
summary(dat$ANTHER_A[dat$ANTHER_A<=2.5])
dat$ANTHER_A[dat$ANTHER_A>2.5] <- NA

# Use petal correlations to determine a PCA score for each individual. 
# The average PCA score per clone is our measure of flower size
# Variables are all measured in mm and the variances are all similar.
# Therefore use an unscaled PCA.
cov(dat[,c('SEPAL_L','UPPERPTL_L','LOWERPTL_L','KEEL_L','KEEL_A')],use="pairwise.complete.obs")
cor(dat[,c('SEPAL_L','UPPERPTL_L','LOWERPTL_L','KEEL_L','KEEL_A')],use="pairwise.complete.obs")
# Exclude sepal length
# keel area is better correlated with upper and lower petal lenghts than is keel length.
X <- dat[,c('UPPERPTL_L','LOWERPTL_L','KEEL_L','KEEL_A')]
theserows <- which(rowSums((!is.na(X))*1)==4)
X <- as.matrix(X[theserows,])
prin.x<-prcomp(X,center=F, scale=F)
# fa.x <- factanal(X,1)
print(prin.x)
summary(prin.x) # PC1 explains 99.4% of the variation in all four petal metrics.
print(fa.x)
names(fa.x)
# fa.x$loadings

par(mfrow=c(2,2));
biplot(prin.x)
qqnorm(prin.x$x[,1]);qqline(prin.x$x[,1]);
qqnorm(prin.x$x[,2]);qqline(prin.x$x[,2]);
plot(-prin.x$x[,1],X[,4]) # very good.
# reverse the sign to properly scale with original measures
dat$Flower.Size[theserows] <- - prin.x$x[,1] 
pairs(dat[,c('UPPERPTL_L','LOWERPTL_L','KEEL_L','KEEL_A','Flower.Size')])
cor(dat[,c('UPPERPTL_L','LOWERPTL_L','KEEL_L','KEEL_A','Flower.Size')],use='complete.obs')
# PC scores are highly correlated with all three petal metrics, higher than 
# their pairwise correlations
# Condense stamen traits as well.
X <- cbind(dat[,c('PolOut','PolIn')],log1p(dat$PolZone))
theserows <- which(rowSums((!is.na(X))*1)==3)
X <- as.matrix(X[theserows,])
prin.x<-prcomp(X,center=F, scale=F)
summary(prin.x) # PC1 explains 99.7% of the variation in all three stamen metrics.

par(mfrow=c(2,2));
biplot(prin.x)
qqnorm(prin.x$x[,1]);qqline(prin.x$x[,1]);
qqnorm(prin.x$x[,2]);qqline(prin.x$x[,2]);
plot(-prin.x$x[,1],X[,3]) # very good.
# reverse the sign to properly scale with original measures
dat$Stamen[theserows] <- - prin.x$x[,1] # PC representing pollen distance and zone width
# Finally, calculate style curvature.
plot(dat$STYLE_C,(dat$STYLE_C-dat$STYLE_L)/dat$STYLE_C)
dat$StyleCurve <- (dat$STYLE_C-dat$STYLE_L)/dat$STYLE_C
ArcLength <- function(A,L) { # INPUT = Arc length (A) and Secand length (L)
  theta <- seq(0.001,2*pi,0.001) # uses a range of angles (0,2*pi) to find one angle 
  A.Length <- sqrt((L/2/tan(theta/2))^2 + L^2/4)*theta
  centralangle <- median(theta[which(abs(A.Length-A)<=0.01)])
  radius <- sqrt((L/2/tan(centralangle/2))^2 + L^2/4)
  return(list(centralangle,radius))
}
ArcLength(dat$STYLE_C[1],dat$STYLE_L[1])
ArcLength(8.9,8.3)[[1]]
dat$StyleAngle <- dat$StyleRadius <- NA
for(i in 1:nrow(dat)) {
  dat$StyleAngle[i] <- ArcLength(dat$STYLE_C[i],dat$STYLE_L[i])[[1]]
  dat$StyleRadius[i] <- log(ArcLength(dat$STYLE_C[i],dat$STYLE_L[i])[[2]])
}

X <- dat[,c('STYLE_L','StyleAngle','StyleRadius')]
theserows <- which(rowSums((!is.na(X))*1)==3)
X <- as.matrix(X[theserows,])
prin.x<-prcomp(X,center=F, scale=F)
summary(prin.x) # PC1 explains 99.3% of the variation in all style metrics.
prin.x$rotation
par(mfrow=c(2,2));
biplot(prin.x)
qqnorm(prin.x$x[,1]);qqline(prin.x$x[,1]);
qqnorm(prin.x$x[,2]);qqline(prin.x$x[,2]);
plot(-prin.x$x[,1],X[,1]) # very good.
# reverse the sign to properly scale with original measures
dat$Style <- NA
dat$Style[theserows] <- - prin.x$x[,1] # PC representing pollen distance and zone width

graphics.off(); x11(width=10, height=10)
grys <- c('gray80','gray60','gray40','gray20')
pairs(dat[,c('Flower.Size','STYLE_L','Style','StyleAngle','StyleRadius',
             'Stamen','SAD')],lower.panel = NULL,gap=0.5,
      col=c('gray80','gray60','gray40','gray20')[dat$STAGE],row1attop=FALSE,pch=19,
      c('Flower\\nsize (PC1)','Style\\nlength','Style\\n(PC1)','Style\\nangle','Style\\nradius',
        'Stamen\\n(PC1)','SAD'))
# These two to show metavariable correlations
# Simple function to add histograms to scatterplot matrix along the diagonal.
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = grys[3])
}
graphics.off();w=h=10;x11(width=w,height=h)
pairs(dat[,c('SAD','STYLE_L','Flower.Size','UPPERPTL_L','LOWERPTL_L','KEEL_L',
             'KEEL_A')],lower.panel = NULL,gap=0.5,pch=19,col=grys[dat$STAGE],
      row1attop=FALSE)
pairs(cbind(dat[,c('SAD','STYLE_L','Flower.Size','Stamen','PolIn','PolOut')],
            log1p(dat$PolZone)),lower.panel = NULL,gap=0.5,pch=19,
      col=grys[dat$STAGE],row1attop=FALSE)
# All together
graphics.off();w=h=10;x11(width=w,height=h)
pairs(cbind(dat[,c('UPPERPTL_L','LOWERPTL_L','KEEL_L','KEEL_A','SAD','STYLE_L',
                   'PolOut','PolIn')],log1p(dat$PolZone)),lower.panel = NULL,
      gap=0.5,pch=19,col=grys[dat$STAGE],row1attop=FALSE,text.panel=NULL,
      diag.panel=panel.hist )
# Save Figure
# dev.copy(pdf,'Figures/FloralTrait_SPlotM.pdf',width=w,height=h)
# dev.off()
# dev.copy(win.metafile,'Figures/FloralTrait_SPlotM.emf',width=w,height=h)
# dev.off()

x11()
pairs(dat[dat$POPULATION=='BT',c('Flower.Size','STYLE_L','Stamen','SAD')],
      lower.panel = NULL,gap=0.5, col=grys[dat$STAGE[dat$POPULATION=='BT']],
      row1attop=FALSE,main='BT',pch=19)
x11()
pairs(dat[dat$POPULATION=='EF',c('Flower.Size','STYLE_L','Stamen','SAD')],
      lower.panel = NULL,gap=0.5,col=grys[dat$STAGE[dat$POPULATION=='EF']],
      row1attop=FALSE,main='EF',pch=19)
x11()
pairs(dat[dat$POPULATION=='TMC',c('Flower.Size','STYLE_L','Stamen','SAD')],
      lower.panel = NULL,gap=0.5,col=grys[dat$STAGE[dat$POPULATION=='TMC']],
      row1attop=FALSE,main='TMC',pch=19)
#_______________________________________________________________________________

#_______________________________________________________________________________
# -- | 2. Stigmatic receptivity data |-------------------------------------------
# headers = "population" "genotype"   "clone"      "stage"      "score" 
# score 1-5 from no to highest bubbling of hydrogen peroxide on stigma
sdat<-read.csv("Data/2_Spigler_etal_CollinsiaStigmaRecep.csv",header=TRUE)
names(sdat)
str(sdat)
summary(sdat)
boxplot(log(sdat$score)~sdat$stage)
sdat$plantID <- paste(sdat$genotype,'_',sdat$clone,sep='') 
plantIDs <- levels(as.factor(sdat$plantID))
newdat <- data.frame('plantID'=plantIDs,'pop'=NA,'recep2'=NA,'recep3'=NA,'recep4'=NA)
for (i in 1:length(plantIDs)) {
  thisdat <- subset(sdat,plantID==plantIDs[i])
  newdat$recep2[i] <- mean(thisdat$score[thisdat$stage==2],na.rm=T)
  newdat$recep3[i] <- mean(thisdat$score[thisdat$stage==3],na.rm=T)
  newdat$recep4[i] <- mean(thisdat$score[thisdat$stage==4],na.rm=T)
  newdat$pop[i] <- as.character(thisdat$population[1])
  newdat$gen[i] <- thisdat$genotype[1]
}
# newdat will be shared with dat below
#_______________________________________________________________________________

#_______________________________________________________________________________
# -- | 3. Autonomous selfing ability | -----------------------------------------
# Load the third data file. This data contains the proportion of fruits set per 
# flowers initiated at the whole plant level. Three clones were used for floral 
# traits, and three clones were used for selfing ability
asdat <- read.csv('Data/3_Spigler_etal_CollinsiaAutoSelfing.csv',header=T)
names(asdat)
str(asdat)
summary(asdat)

asdat[asdat$GTYPE==100,] # clones 1, 2, 5
dat[dat$GENOTYPE==100,]  # clones 0, 3, 4 (different from AS clones)
#_______________________________________________________________________________
# -- | 4. Functions to be used below | -----------------------------------------
# It simplifies the coding below to define a few functions here.

# mSEpredict | computes the means and standard errors from an ANOVA model
mSEpredict <- function(model, data) {
  m.est <- tapply(predict(model,type='response',se.fit=T)$fit,
                  list(data$Population),median)
  m.se <- tapply(predict(model,type='response',se.fit=T)$se.fit,
                 list(data$Population),median)
  return(rbind('Means'=m.est,'SEs'=m.se))
}

# CI Plot Functions | plots means as points with error bars:
cexaxs <- 1.1; cexlbs <- 1.4; cexmains <- 0.8; cexpts <- 1.2
lwds=2; alwds=2
graphics.off(); w=3.15;h=4.6;m1=4;m2=3;m3=0;m4=1; x11(width=w,height=h); 
# graphics.off(); w=7.08;h=4.6; x11(width=w,height=h); 
# Panel A:
cifun <- function(means,sem,xlabname,ylabname,siglets,panlet,collist,
                  yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet) {
  plotCI(means,uiw=sem,xlab=xlabname, ylab='',axes=T,xaxt='n',
         ylim=c(min(means*(2-yfac)-sem),max(means*yfac+sem)),bty='n',
         col=collist,cex.lab=clab,cex.axis=caxis,sfrac=.012,lwd=1,
         gap=0.0,pch=pchs,cex=pcex, pt.bg=par("bg"),xlim=c(0.2,3.8),add=F)
  box(bty='l',lwd=alwds)
  title(ylab=ylabname,line=2.5,cex.lab=clab)
  for(i in 1:3) { text(i,(means*yfac2+sem)[i],siglets[i],cex=clet,pos=3)}
  axis(1,1:3,labels=c('BT','EF','TMC'),cex.axis=cxlab)
  mtext(panlet, side=3, line=.25, at=0.1, outer=F, cex=clet, font=1)
}
# Definitions: -------------------------------------------|
# argnames = group names
# siglets = significance letters for pairwise comparisions
# panlet = panel letter
# collist = Color palette
# yfac, yfac2 = Scaling
# clab,caxis,cnames,pcex,cxlab,clet = Scaling
# pchs = Point types
# --------------------------------------------------------|

# Plot parameter estimates | plots means as points with error bars:
# Receives output from MCMCglmm object
plot.estimates <- function(x) {
  if (class(x) != "summary.mcmc")
    x <- summary(x)
  n <- dim(x$solutions)[1]
  par(mar=c(2, 9, 4, 1))
  plot(x$solutions[,1], n:1,
       yaxt="n", ylab="",
       xlim=range(x$solutions[,2:3])*1.2,
       pch=19,
       main="Posterior means and 95% credible intervals")
  grid()
  axis(2, at=n:1, rownames(x$solutions), las=2)
  arrows(x$solutions[,2], n:1, x$solutions[,3], n:1, code=0)
  abline(v=0, lty=2)
}

# Different from zero??
sigdif <- function(CovMat) { # receives upper and lower 95% CI
  sig_r <- which({CovMat[,1]>0 & CovMat[,2]>0} | {CovMat[,1]<0 & CovMat[,2]<0})
  return(as.numeric(sig_r))
}

#_______________________________________________________________________________
# -- | 5. Aggregate genotypic dataset | ----------------------------------------
# To combine floral trait, stigmatic receptivity, and autonomous selfing data, we
# need to calculate genotypic means. Floral trait and stigmatic receptivity were 
# measured on the same clones, so create a clonal mean data set for them first. 
# Then calculate genotypic means before merging with selfing data, after also 
# calculating genotypic means for selfing. 
# 
# Identify columns in floral trait data to aggregate, and create empty objects
# * cdats is a list containing four data sets, each with clone means
# * gdats is also a list containing four data sets, each with genotype means

# thesecols <- (c("UPPERPTL_L","LOWERPTL_L","KEEL_L","KEEL_A","SAD",
#                 "Flower.Size",'STYLE_L','STYLE_C','PolOut','PolIn',
#                 'PolZone','Stamen','GENOTYPE'))
thesecols <- (c("UPPERPTL_L","LOWERPTL_L","KEEL_L","KEEL_A","SAD","Flower.Size",
                'STYLE_L','STYLE_C','FILM1_L','FILM2_L','FILM3_L','FILM4_L',
                'PolIn','PolOut','PolZone','Stamen','GENOTYPE','aSAD'))
dat[1:20,thesecols]
# Data organizational note: Each plant/clone has a developmental trajectory. To store this information, set up a list fo datasets. Each data set contains all floral trait values for each floral stage, i.e., 1 to 4 stamens dehisced. What we want is to estimate the genetic varition underlying six key traits: flower size, stamen length, style length, stigma-anther distance, stigmatic receptivity, and selfing ability. 
# Selfing ability does not have a developmental trajectory and therfore requires only a mean. 

# Stig Recep data already aggregated
newdat$recep1 <- 1 # set receptivity at stage 1 to 1, which is no activity
newdat$Population <- newdat$pop
receps <- c('recep1','recep2','recep3','recep4')
gdats <- cdats <- list(dat,dat,dat,dat) # place holders for 4 generated datasets

for(i in 1:4) {
  theserows <- which(dat$STAGE==i)
  # average within clones first for main dataset
  dat1 <- aggregate(dat[theserows,thesecols],list(as.factor(dat[theserows,
                    'plantID'])),mean,na.rm=T)
  dat1$plantID <- dat1$Group.1
  cdat <- merge(dat1,newdat)[,c('plantID','Population',thesecols,receps[i])]
  # Now average over clones within genotypes
  dat1 <- aggregate(cdat[,c(thesecols,receps[i])],list(as.factor(cdat[,
                          'GENOTYPE'])),mean,na.rm=T)
  # Aggregate Automatic Selfing Ability data
  dat2 <- aggregate(asdat[,'AA'],list(as.factor(asdat[,'GTYPE'])),mean,na.rm=T)
  names(dat2)[2] <- 'AS'
  # create Population identifier
  dat2$Population <- factor(NA,levels=1:3,labels=c('BT','EF','TMC'))
  thing <- as.numeric(as.character(dat2$Group.1))
  dat2$Population[thing>=100 & thing<200] <- 'EF'
  dat2$Population[thing>=200 & thing<300] <- 'BT'
  dat2$Population[thing>=300 & thing<400] <- 'TMC'
  gdat <- merge(dat1,dat2)
  gdats[[i]] <- gdat
  cdats[[i]] <- cdat
  rm(theserows,cdat,dat1,dat2,thing,gdat)
}
summary(cdats)
summary(cdats[[1]])
summary(gdats)
summary(gdats[[1]])
for (i in 1:4) {
  gdats[[i]]$Stage <- i
  names(gdats[[i]])[20] <- 'SR'
  gdats[[i]]$lnSR <- log(gdats[[i]]$SR) 
  cdats[[i]]$Stage <- i
  names(cdats[[i]])[21] <- 'SR'
  cdats[[i]]$lnSR <- log(cdats[[i]]$SR) 
}

Gdat <- rbind(gdats[[1]],gdats[[2]],gdats[[3]],gdats[[4]])
Cdat <- rbind(cdats[[1]],cdats[[2]],cdats[[3]],cdats[[4]])
names(Gdat)
names(Cdat)

Ntab <- matrix(NA,nrow=3,ncol=5)
dimnames(Ntab) <- list(c('BT','EF','TMC'),c('Genotype','AS Clones',
                         'FT Clones','SR Flowers','Morph Flowers'))
Ntab[1,1] <- length(levels(as.factor(Gdat$GENOTYPE[Gdat$Population=='BT'])))
Ntab[2,1] <- length(levels(as.factor(Gdat$GENOTYPE[Gdat$Population=='EF'])))
Ntab[3,1] <- length(levels(as.factor(Gdat$GENOTYPE[Gdat$Population=='TMC'])))
Ntab[1,2] <- nrow(asdat[asdat$POP =='BT',])
Ntab[2,2] <- nrow(asdat[asdat$POP =='EF',])
Ntab[3,2] <- nrow(asdat[asdat$POP =='TMC',])
Ntab[1,3] <- length(levels(droplevels(Cdat$plantID[Cdat$Population=='BT'])))
Ntab[2,3] <- length(levels(droplevels(Cdat$plantID[Cdat$Population=='EF'])))
Ntab[3,3] <- length(levels(droplevels(Cdat$plantID[Cdat$Population=='TMC'])))
Ntab[1,4] <- nrow(sdat[sdat$population =='BT',])
Ntab[2,4] <- nrow(sdat[sdat$population =='EF',])
Ntab[3,4] <- nrow(sdat[sdat$population =='TMC',])
Ntab[1,5] <- nrow(dat[dat$STAGE > 0 & dat$POPULATION=='BT',])
Ntab[2,5] <- nrow(dat[dat$STAGE > 0 & dat$POPULATION=='EF',])
Ntab[3,5] <- nrow(dat[dat$STAGE > 0 & dat$POPULATION=='TMC',])
Ntab

x11()
pairs(Gdat[Gdat$Population=='BT',c('Flower.Size','STYLE_L','Stamen','SAD',
                                         'SR','AS')],lower.panel = NULL,gap=0.5,
    col=grys[Gdat$Stage[Gdat$Population=='BT']],row1attop=FALSE,main='BT',pch=19)

x11()
pairs(Gdat[Gdat$Population=='EF',c('Flower.Size','STYLE_L','Stamen','SAD',
                                         'SR','AS')],lower.panel = NULL,gap=0.5,
      col=grys[Gdat$Stage[Gdat$Population=='EF']],row1attop=FALSE,main='EF',pch=19)
x11()
pairs(Gdat[Gdat$Population=='TMC',c('Flower.Size','STYLE_L','Stamen','SAD',
                                         'SR','AS')],lower.panel = NULL,gap=0.5,
      col=grys[Gdat$Stage[Gdat$Population=='TMC']],row1attop=FALSE,
      main='TMC',pch=19)
## ln(SR) doesn't really seem different from SR when looking at correlations 
#  within a single stage. 

#--Figure 3--------------------------------------------------------------------|
cexaxs <- 1.1; cexlbs <- 1.4; cexmains <- 0.8; cexpts <- 1.2
lwds=2; alwds=2
graphics.off(); w=3.5;h=4.6;m1=4;m2=3;m3=0;m4=1; x11(width=w,height=h); 
# graphics.off(); w=7.08;h=4.6; x11(width=w,height=h); 
# Panel A:
par(mfrow=c(3,2),oma=c(0.4,0.4,0,0),mar=c(m3,m2,m1,m4))
ylims1 <- range(Gdat$SAD); xlims1 <- range(Gdat$STYLE_L)
plot(Gdat[Gdat$Population=='BT','STYLE_L'],Gdat[Gdat$Population=='BT','SAD'],
     col=grys[Gdat$Stage[Gdat$Population=='BT']],cex=cexpts,pch=19,
     xlab='',ylab='',main='',axes=FALSE,ylim=ylims1, xlim=xlims1)
box(bty='l',lwd=alwds)
axis(1,at=seq(4,10,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=-1.25,
     labels=c('','','',''))
axis(2,at=seq(-2,4,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=0.75)
abline(h=0,lty=2,lwd=lwds,col='gray60')
# Stigmatic receptivity v SAD
ylims2 <- range(Gdat$SR,na.rm=T); xlims2 <- range(Gdat$SAD)
# Panel D:
par(mar=c(m3,m2,m1,m4))
plot(Gdat[Gdat$Population=='BT','SAD'],Gdat[Gdat$Population=='BT','SR'],
     col=grys[Gdat$Stage[Gdat$Population=='BT']],cex=cexpts,pch=19,
     xlab='',ylab='',main='',axes=FALSE,ylim=ylims2, xlim=xlims2)
box(bty='l',lwd=alwds)
axis(1,at=seq(-2,4,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=-1.25,
     labels=c('','','',''))
axis(2,at=seq(1,5,1),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=0.75)
abline(v=0,lty=2,lwd=lwds,col='gray60')
mtext('BT',cex=cexmains,side=3,at=4,line=-0.5)
# Panel B:
par(mar=c(m1/2,m2,m1/2,m4))
plot(Gdat[Gdat$Population=='EF','STYLE_L'],Gdat[Gdat$Population=='EF','SAD'],
     col=grys[Gdat$Stage[Gdat$Population=='EF']],cex=cexpts,pch=19,
     xlab='',ylab='',main='',axes=FALSE,ylim=ylims1, xlim=xlims1)
box(bty='l',lwd=alwds)
axis(1,at=seq(4,10,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=-1.25,
     labels=c('','','','',''))
axis(2,at=seq(-2,4,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=0.75)
abline(h=0,lty=2,lwd=lwds,col='gray60')
title(ylab='Stigma-anther distance (mm)',line=-1.15,cex.lab=cexlbs,outer=TRUE)
# Panel E:
par(mar=c(m1/2,m2,m1/2,m4))
plot(Gdat[Gdat$Population=='EF','SAD'],Gdat[Gdat$Population=='EF','SR'],
     col=grys[Gdat$Stage[Gdat$Population=='EF']],cex=cexpts,pch=19,
     xlab='',ylab='',main='',axes=FALSE,ylim=ylims2, xlim=xlims2)
box(bty='l',lwd=alwds)
axis(1,at=seq(-2,4,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=-1.25,
     labels=c('','','','',''))
axis(2,at=seq(1,5,1),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=0.75)
title(ylab='Stigmatic receptivity',line=1.7,cex.lab=cexlbs)
mtext('EF',cex=cexmains,side=3,at=4,line=-0.5)
abline(v=0,lty=2,lwd=lwds,col='gray60')
# Panel C:
par(mar=c(m1,m2,m3,m4))
plot(Gdat[Gdat$Population=='TMC','STYLE_L'],Gdat[Gdat$Population=='TMC','SAD'],
     col=grys[Gdat$Stage[Gdat$Population=='TMC']],cex=cexpts,pch=19,
     xlab='',ylab='',main='',axes=FALSE,ylim=ylims1, xlim=xlims1)
box(bty='l',lwd=alwds)
axis(1,at=seq(4,10,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=-0.7)
axis(2,at=seq(-2,4,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=0.75)
abline(h=0,lty=2,lwd=lwds,col='gray60')
title(xlab='Style length (mm)\\n',line=3,cex.lab=cexlbs)
# Panel F:
par(mar=c(m1,m2,m3,m4))
plot(Gdat[Gdat$Population=='TMC','SAD'],Gdat[Gdat$Population=='TMC','SR'],
     col=grys[Gdat$Stage[Gdat$Population=='TMC']],cex=cexpts,pch=19,
     xlab='',ylab='',main='',axes=FALSE,ylim=ylims2, xlim=xlims2)
box(bty='l',lwd=alwds)
axis(1,at=seq(-2,4,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=-0.7)
axis(2,at=seq(1,5,1),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=0.75)
abline(v=0,lty=2,lwd=lwds,col='gray60')
title(xlab='Stigma-anther\\ndistance (mm)',line=3,cex.lab=cexlbs)
mtext('TMC',cex=cexmains,side=3,at=4,line=-0.5)
# dev.copy(pdf,'Figures/SAD_SL_SR_Scatterplots.pdf',width=w,height=h)
# dev.off()
# dev.copy(win.metafile,'Figures/SAD_SL_SR_Scatterplots.emf',width=6.5,height=5)
# dev.off()

# #--Figure 3--------------------------------------------------------------------|
# cexaxs <- 1.0; cexlbs <- 1.2; cexmains <- 1.1; cexpts <- 1.2
# lwds=2; alwds=2
# graphics.off(); w=3.5;h=2.6;m1=2.45;m2=3.8;m3=1.5;m4=0; x11(width=w,height=h); 
# # graphics.off(); w=7.08;h=4.6; x11(width=w,height=h); 
# # Panel A:
# par(mfrow=c(2,3),oma=c(0.4,0.2,0.3,0),mar=c(m1,m2,m3,m4))
# ylims1 <- range(Gdat$SAD); xlims1 <- range(Gdat$STYLE_L)
# plot(Gdat[Gdat$Population=='BT','STYLE_L'],Gdat[Gdat$Population=='BT','SAD'],
#      col=grys[Gdat$Stage[Gdat$Population=='BT']],cex=cexpts,pch=19,
#      xlab='',ylab='',main='',axes=FALSE,ylim=ylim1s, xlim=xlims1)
# box(bty='l',lwd=alwds)
# axis(1,at=seq(4,10,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=-1.25)
# axis(2,at=seq(-2,4,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=0.95)
# title(main='BT',cex.lab=cexmains,line=0.15)
# abline(h=0,lty=2,lwd=lwds,col='gray60')
# # Panel B:
# par(mar=c(m1,m2/2,m3,m2/2))
# plot(Gdat[Gdat$Population=='EF','STYLE_L'],Gdat[Gdat$Population=='EF','SAD'],
#      col=grys[Gdat$Stage[Gdat$Population=='EF']],cex=cexpts,pch=19,
#      xlab='',ylab='',main='',axes=FALSE,ylim=ylims1, xlim=xlims1)
# box(bty='l',lwd=alwds)
# axis(1,at=seq(4,10,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=-1.25)
# axis(2,at=seq(-2,4,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=0.95,
#      labels=c('','','',''))
# title(main='EF',cex.lab=cexmains,line=0.15)
# abline(h=0,lty=2,lwd=lwds,col='gray60')
# # Panel C:
# par(mar=c(m1,m4,m3,m2))
# plot(Gdat[Gdat$Population=='TMC','STYLE_L'],Gdat[Gdat$Population=='TMC','SAD'],
#      col=grys[Gdat$Stage[Gdat$Population=='TMC']],cex=cexpts,pch=19,
#      xlab='',ylab='',main='',axes=FALSE,ylim=ylims1, xlim=xlims1)
# box(bty='l',lwd=alwds)
# axis(1,at=seq(4,10,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=-1.25)
# axis(2,at=seq(-2,4,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=0.95,
#      labels=c('','','',''))
# title(main='TMC',cex.lab=cexmains,line=0.15)
# title(ylab='Stigma-anther distance (mm)',line=1.6,cex.lab=cexlbs)
# abline(h=0,lty=2,lwd=lwds,col='gray60')
# title(xlab='Style length (mm)',line=1.4,cex.lab=cexlbs)
# # Stigmatic receptivity v SAD
# ylims2 <- range(Gdat$SR,na.rm=T); xlims2 <- range(Gdat$SAD)
# # Panel D:
# par(mar=c(m1,m2,m3,m4))
# plot(Gdat[Gdat$Population=='BT','SAD'],Gdat[Gdat$Population=='BT','SR'],
#      col=grys[Gdat$Stage[Gdat$Population=='BT']],cex=cexpts,pch=19,
#      xlab='',ylab='',main='',axes=FALSE,ylim=ylims2, xlim=xlims2)
# box(bty='l',lwd=alwds)
# axis(1,at=seq(-2,4,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=-1.25)
# axis(2,at=seq(1,5,1),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=0.95)
# abline(v=0,lty=2,lwd=lwds,col='gray60')
# # Panel E:
# par(mar=c(m1,m2/2,m3,m2/2))
# plot(Gdat[Gdat$Population=='EF','SAD'],Gdat[Gdat$Population=='EF','SR'],
#      col=grys[Gdat$Stage[Gdat$Population=='EF']],cex=cexpts,pch=19,
#      xlab='',ylab='',main='',axes=FALSE,ylim=ylims2, xlim=xlims2)
# box(bty='l',lwd=alwds)
# axis(1,at=seq(-2,4,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=-1.25)
# axis(2,at=seq(1,5,1),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=0.95,
#      labels=c('','','','',''))
# title(ylab='Stigmatic\\nreceptivity',line=1.6,cex.lab=cexlbs)
# abline(v=0,lty=2,lwd=lwds,col='gray60')
# # Panel F:
# par(mar=c(m1,m4,m3,m2))
# plot(Gdat[Gdat$Population=='TMC','SAD'],Gdat[Gdat$Population=='TMC','SR'],
#      col=grys[Gdat$Stage[Gdat$Population=='TMC']],cex=cexpts,pch=19,
#      xlab='',ylab='',main='',axes=FALSE,ylim=ylims2, xlim=xlims2)
# box(bty='l',lwd=alwds)
# axis(1,at=seq(-2,4,2),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=-1.25)
# axis(2,at=seq(1,5,1),cex.axis=cexaxs,lwd.ticks=alwds,tck=-0.05,padj=0.95,
#      labels=c('','','','',''))
# abline(v=0,lty=2,lwd=lwds,col='gray60')
# title(xlab='Stigma-anther\\ndistance (mm)',line=-0.8,cex.lab=cexlbs,outer=TRUE)
# 
# dev.copy(pdf,'Figures/SAD_SL_SR_Scatterplots.pdf',width=w,height=h)
# dev.off()
# # dev.copy(win.metafile,'Figures/SAD_SL_SR_Scatterplots.emf',width=6.5,height=5)
# # dev.off()


# Collect correlation matrices.
BTcors <- EFcors <- TMCcors <- array(NA,dim=c(6,6,4))
dimnames(BTcors) <- list(c('FS','Sty','Sta','SAD','SR','AS'),c('FS','Sty','Sta',
                          'SAD','SR','AS'),c('S1','S2','S3','S4'))
dimnames(EFcors) <- dimnames(TMCcors) <- dimnames(BTcors)
BTcors[-5,-5,1] <- round(cor(gdats[[1]][gdats[[1]]$Population=='BT',
                          c('Flower.Size','STYLE_L','Stamen','SAD','AS')],
                          use="pairwise.complete.obs"),3)
EFcors[-5,-5,1] <- round(cor(gdats[[1]][gdats[[1]]$Population=='EF',
                          c('Flower.Size','STYLE_L','Stamen','SAD','AS')],
                             use="pairwise.complete.obs"),3)
TMCcors[-5,-5,1] <- round(cor(gdats[[1]][gdats[[1]]$Population=='TMC',
                          c('Flower.Size','STYLE_L','Stamen','SAD','AS')],
                             use="pairwise.complete.obs"),3)
for (i in 2:4) {
  BTcors[,,i] <- round(cor(gdats[[i]][gdats[[i]]$Population=='BT',
                          c('Flower.Size','STYLE_L','Stamen','SAD','SR','AS')],
                               use="pairwise.complete.obs"),3)
  EFcors[,,i] <- round(cor(gdats[[i]][gdats[[i]]$Population=='EF',
                          c('Flower.Size','STYLE_L','Stamen','SAD','SR','AS')],
                               use="pairwise.complete.obs"),3)
  TMCcors[,,i] <- round(cor(gdats[[i]][gdats[[i]]$Population=='TMC',
                          c('Flower.Size','STYLE_L','Stamen','SAD','SR','AS')],
                                use="pairwise.complete.obs"),3)
}
BTcors
BTbigCor <- EFbigCor <- TMCbigCor <- matrix (NA,nrow=6*4,ncol=6)
dimnames(BTbigCor) <- list(c(paste0('AS.S',1:4),paste0('SR.S',1:4),
                             paste0('SAD.S',1:4),paste0('Sta.S',1:4),
                             paste0('Sty.S',1:4),paste0('FS.S',1:4)),
                           dimnames(BTcors)[[1]])
dimnames(EFbigCor) <- dimnames(TMCbigCor) <- dimnames(BTbigCor)
# Fill the lower diagonal of bigCors
for (i in 2:6) {
  BTbigCor[(1:4)+20,i] <- BTcors[i,1,1:4] # FS w/ i
  EFbigCor[(1:4)+20,i] <- EFcors[i,1,1:4] # FS w/ i
  TMCbigCor[(1:4)+20,i] <- TMCcors[i,1,1:4] # FS w/ i
}
for (i in 3:6) {
  BTbigCor[(1:4)+16,i] <- BTcors[i,2,1:4] # Sty w/ i
  EFbigCor[(1:4)+16,i] <- EFcors[i,2,1:4] # Sty w/ i
  TMCbigCor[(1:4)+16,i] <- TMCcors[i,2,1:4] # Sty w/ i
}
for (i in 4:6) {
  BTbigCor[(1:4)+12,i] <- BTcors[i,3,1:4] # Sta w/ i
  EFbigCor[(1:4)+12,i] <- EFcors[i,3,1:4] # Sta w/ i
  TMCbigCor[(1:4)+12,i] <- TMCcors[i,3,1:4] # Sta w/ i
}
for (i in 5:6) {
  BTbigCor[(1:4)+8,i] <- BTcors[i,4,1:4] # SAD w/ i
  EFbigCor[(1:4)+8,i] <- EFcors[i,4,1:4] # SAD w/ i
  TMCbigCor[(1:4)+8,i] <- TMCcors[i,4,1:4] # SAD w/ i
}
BTbigCor[(1:4)+4,6] <- BTcors[6,5,1:4] # SR w/ i
EFbigCor[(1:4)+4,6] <- EFcors[6,5,1:4] # SR w/ i
TMCbigCor[(1:4)+4,6] <- TMCcors[6,5,1:4] # SR w/ i
BTbigCor
EFbigCor
TMCbigCor
# Figure 
# BigCor <- cbind(BTbigCor[,1],EFbigCor[,1],TMCbigCor[,1],
#                BTbigCor[,2],EFbigCor[,2],TMCbigCor[,2],
#                BTbigCor[,3],EFbigCor[,3],TMCbigCor[,3],
#                BTbigCor[,4],EFbigCor[,4],TMCbigCor[,4],
#                BTbigCor[,5],EFbigCor[,5],TMCbigCor[,5],
#                BTbigCor[,6],EFbigCor[,6],TMCbigCor[,6])
BigCor <- cbind(TMCbigCor[,1],EFbigCor[,1],BTbigCor[,1],
                TMCbigCor[,2],EFbigCor[,2],BTbigCor[,2],
                TMCbigCor[,3],EFbigCor[,3],BTbigCor[,3],
                TMCbigCor[,4],EFbigCor[,4],BTbigCor[,4],
                TMCbigCor[,5],EFbigCor[,5],BTbigCor[,5],
                TMCbigCor[,6],EFbigCor[,6],BTbigCor[,6])

#--Figure 2--------------------------------------------------------------------|
cexaxs=0.8; cexlbs=1.1;cexcor=0.7
lwds=2; alwds=2
graphics.off(); w=7.25;h1=6.5; x11(width=w,height=h1); 
rgb.palette <- colorRampPalette(c('#0093F2','white','#FF8E00'), space = "rgb")
par(oma=c(0,2.5,1.5,0), mar=c(0,1,1,0))
# BigCor[c(21:24,17:20,13:16,9:12,5:8,1:4),]
image(x=1:nrow(BigCor),y=1:ncol(BigCor),
      z=BigCor[c(21:24,17:20,13:16,9:12,5:8,1:4),],col=rgb.palette(220),
      xlab='',ylab='',axes=F,xlim=c(0.5,24.5),ylim=c(0.5,18.5))
h <- c(3,seq(3,18,3))
v <- c(0,seq(4,20,4),20)
for (i in 1:length(v)) {
    lines(c(0.5,v[i]+0.5),c(h[i]+0.5,h[i]+0.5),lty=1,lwd=1,col='gray40')
    lines(c(v[i]+0.5,v[i]+0.5),c(18.5,h[i]+0.5),lty=1,lwd=1,col='gray40')
}
axis(3,at=1:(nrow(BigCor)-4),label=rep(1:4,5),tick=FALSE,padj=1.75,cex.axis=cexaxs)
axis(2,at=4:ncol(BigCor),label=rep(c('TMC','EF','BT'),5),tick=FALSE,
     las=2,hadj=0.37,cex.axis=cexaxs)
mtext('                 Population',side=2,line=1.2,outer=TRUE,cex=cexlbs)
mtext('Progression of anther dehiscence               ',side=3,line=0,
      outer=TRUE,cex=cexlbs)
# traitnames <- c('FS','SL','ST','SAD','SR','AS')
traitnames <- c('Flower size\\n(PC1)\\n','Style\\nlength\\n','Pollen\\nzone\\n(PC1)',
                'Stigma\\nanther\\ndistance','Stigmatic\\nreceptivity\\n','Auto.\\nselfing\\nability')
for(i in 1:6) {
  text(seq(4,24,4)[i]-1.5,seq(3,18,3)[i]-1,traitnames[i],cex=cexlbs)
}
for (i in 1:nrow(BigCor)) {
  for (j in 1:ncol(BigCor))
    if(!is.na(BigCor[i,j])) {
      if(BigCor[i,j] > 0.3 | BigCor[i,j] < -0.3) {
        text(c(21:24,17:20,13:16,9:12,5:8,1:4)[i],j,
             round(BigCor[i,j],1),cex=cexcor, col='gray40', font=2)
      }
    }
}
#--inset---
rho.lvls <- seq(0.8,-0.8,-0.1)
rho.lvls[-c(1,5,9,13,17)] <- ''
par(new=TRUE) # overlay existing plot
par(mar=c(0,0,0,0)) # strip out the margins for the inset plot
par(fig=c(0.75,0.95,0.05,0.575)) # fig shrinks and places relative to figure region
plot(1:5,1:5,type='n',axes=FALSE,xlab='',ylab='')
par(lend=1)
legend('topleft',legend=rho.lvls,col=rgb.palette(17)[17:1],lty=1,lwd=15,
       y.intersp=0.63,cex=1.1,bty='n',inset=0,title.adj=0.2) 
polygon(x=c(1.37,2.45,2.45,1.37),y=c(5.02,5.02,1.83,1.83),border='gray40',lwd=1)
# Save Figure
# dev.copy(pdf,'Figures/FloralTrait_Correlations.pdf',width=w,height=h1)
# dev.off()
# dev.copy(win.metafile,'Figures/FloralTrait_Correlations.emf',width=8,height=7)
# dev.off()
#___________________________________________________________________________----
#_______________________________________________________________________________
# -- | 6. Multivariate Analysis of Variance | ----------------------------------
# Estimate trait differences between populations
aggregate(gdats[[4]],list(gdats[[4]][,'Population']),mean,na.rm=T)[,c(-2,-23)]

# -- 6.A | Petal trait differences among populations (Stage 4 only)
thesecols <- c("UPPERPTL_L","LOWERPTL_L","KEEL_L","KEEL_A")
theserows <- which(rowSums((!is.na(gdats[[4]][,thesecols]))*1)==length(thesecols))
Yp <- as.matrix(gdats[[4]][theserows,thesecols])
mp <- manova(Yp~Population, data=gdats[[4]][theserows,])
summary(mp)
summary.aov(mp) # for univariate results

# Set up Tukey HSD tests for trait-specific pairwise comparsions
u1 <- aov(Yp[,1]~Population,data=gdats[[4]]) # upper petal length
up.diff <- TukeyHSD(u1)
u2 <- aov(Yp[,2]~Population,data=gdats[[4]]) # lower petal length
lp.diff <- TukeyHSD(u2)
u3 <- aov(Yp[,3]~Population,data=gdats[[4]]) # keel petal length
kl.diff <- TukeyHSD(u3)
u4 <- aov(Yp[,4]~Population,data=gdats[[4]]) # keel petal area
ka.diff <- TukeyHSD(u4)

# 1. Upper petals ARE different among populations: both EF & BT > TMC.
up.diff
# 2. Lower petals ARE different among populations: only BT > TMC.
lp.diff
# 3. Keel petal length is NOT different among populations.
kl.diff
# 4. Keel petal area also is NOT different among populations.
ka.diff

# Collect mean and standard errors from the ANOVAs
up.est <- mSEpredict(u1, gdats[[4]]) 
lp.est <- mSEpredict(u2, gdats[[4]]) 
kl.est <- mSEpredict(u3, gdats[[4]]) 
ka.est <- mSEpredict(u4, gdats[[4]]) 

# -- | Figures | Plot means and standard errors for each petal trait (stage 4 only)
# Set the Controls |
cols <- 'black' # Color palette
clab <- 1.8; pchs <- 19; pcex <- 2;  # Scaling factors, Point type
caxis <- cxlab <- 1.7; clet <- 1.5; # Scaling

graphics.off(); w<-8; h<-8
x11(width=w,height=h); par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(2.5,4.5,3,1.5))  
# Panel a: upper petal length
siglets <- c('a','a','b'); yfac <- 1.01; yfac2 <- 1.003
cifun(up.est[1,],up.est[2,],xlabname='',ylabname='Upper petal length',
      siglets,'a',cols,yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet)
# Panel b: keel length
siglets <- c('a','a','a'); yfac <- 1.01; yfac2 <- 1.003
cifun(kl.est[1,],kl.est[2,],xlabname='',ylabname='Keel length',
      siglets,'b',cols,yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet)
# Panel c: lower petal length
siglets <- c('a','ab','b'); yfac <- 1.01; yfac2 <- 1.003
cifun(lp.est[1,],lp.est[2,],xlabname='',ylabname='Lower petal length',
      siglets,'c',cols,yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet)
# Panel d: Keel area 
siglets <- c('a','a','a'); yfac <- 1.01; yfac2 <- 1.003
cifun(ka.est[1,],ka.est[2,],xlabname='',ylabname='Keel area',
      siglets,'d',cols,yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet)
# Save Figure
# dev.copy(pdf,'Figures/PetalTraitComparisons.pdf',width=w,height=h)
# dev.off()
# dev.copy(win.metafile,'Figures/PetalTraitComparisons.emf',width=w,height=h)
# dev.off()

# -- 6.B | Stamen length differences among populations (Stage 4 only)
thesecols <- c('FILM1_L','FILM2_L','FILM3_L','FILM4_L')
theserows <- which(rowSums((!is.na(gdats[[4]][,thesecols]))*1)==length(thesecols))
Ys <- as.matrix(gdats[[4]][theserows,thesecols])
ms <- manova(Ys~Population, data=gdats[[4]][theserows,])
summary(ms)
summary.aov(ms) # for univariate results

# Set up Tukey HSD tests for trait-specific pairwise comparsions
s1 <- aov(Ys[,1]~Population,data=gdats[[4]]) # Stamen 1 length
s1.diff <- TukeyHSD(s1)
s2 <- aov(Ys[,2]~Population,data=gdats[[4]]) # Stamen 2 length
s2.diff <- TukeyHSD(s2)
s3 <- aov(Ys[,3]~Population,data=gdats[[4]]) # Stamen 3 length
s3.diff <- TukeyHSD(s3)
s4 <- aov(Ys[,4]~Population,data=gdats[[4]]) # Stamen 4 length
s4.diff <- TukeyHSD(s4)

# 1. Stamen 1 length IS different among populations: only BT > TMC.
s1.diff
# 2. Stamen 2 length IS different among populations: only BT > TMC.
s2.diff
# 3. Stamen 3 length IS different among populations: only BT > TMC.
s3.diff
# 4. Stamen 4 length IS different among populations: only BT > TMC.
s4.diff

# Collect mean and standard errors from the ANOVAs
s1.est <- mSEpredict(s1, gdats[[4]]) 
s2.est <- mSEpredict(s2, gdats[[4]]) 
s3.est <- mSEpredict(s3, gdats[[4]]) 
s4.est <- mSEpredict(s4, gdats[[4]]) 

# -- | Figure | Plot means and standard errors for each stamen length at stage 4.
graphics.off(); w<-8; h<-8
x11(width=w,height=h); par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(2.5,4.5,3,1.5))  
# Panel a: Stamen 1 length
siglets <- c('a','ab','b'); yfac <- 1.01; yfac2 <- 1.003
cifun(s1.est[1,],s1.est[2,],xlabname='',ylabname='Stamen 1 length',
      siglets,'a',cols,yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet)
# Panel b: Stamen 2 length
siglets <- c('a','ab','b'); yfac <- 1.01; yfac2 <- 1.003
cifun(s2.est[1,],s2.est[2,],xlabname='',ylabname='Stamen 2 length',
      siglets,'b',cols,yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet)
# Panel c: Stamen 3 length
siglets <- c('a','ab','b'); yfac <- 1.01; yfac2 <- 1.003
cifun(s3.est[1,],s3.est[2,],xlabname='',ylabname='Stamen 3 length',
      siglets,'c',cols,yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet)
# Panel d: Stamen 4 length
siglets <- c('a','ab','b'); yfac <- 1.01; yfac2 <- 1.003
cifun(s4.est[1,],s4.est[2,],xlabname='',ylabname='Stamen 4 length',
      siglets,'d',cols,yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet)
# Save Figure
# dev.copy(pdf,'Figures/StamenLengthComparisons.pdf',width=w,height=h)
# dev.off()
# dev.copy(win.metafile,'Figures/StamenLengthComparisons.emf',width=w,height=h)
# dev.off()

# -- 6.C | Floral trait differences among populations (Stage 4 only)
thesecols <- c('Flower.Size','STYLE_L','Stamen','SAD','SR','AS')
theserows <- which(rowSums((!is.na(gdats[[4]][,thesecols]))*1)==length(thesecols))
Yf <- as.matrix(gdats[[4]][theserows,thesecols])
# Tests of Assumptions



# Graphical Assessment of Multivariate Normality
# If we have p x 1 multivariate normal random vector then the squared Mahalanobis distance between x and µ is going to be chi-square distributed with p degrees of freedom. We can use this fact to construct a Q-Q plot to assess multivariate normality.
x <- as.matrix(Yf) # n x p numeric matrix
center <- colMeans(x) # centroid
n <- nrow(x); p <- ncol(x); cov <- cov(x);
d <- mahalanobis(x,center,cov) # distances
cexaxs <- 0.9; cexlbs <- 1.1; cexmains <- 0.8; cexpts <- 1.2
lwds=2; alwds=2
graphics.off(); w=3.15;h=3.15; x11(width=w,height=h); 
par(oma=c(0.2,0.2,0.2,0.2),mar=c(3.5,3.5,1,1))
qqplot(qchisq(ppoints(n),df=p),d, main="",ylab="",
       xlab='',cex=cexpts,pch=19,bty='n',
       cex.axis=cexaxs,cex.lab=cexlbs)
title(ylab='Mahalanobis distance',line=2.3,cex.lab=cexlbs) 
title(xlab=expression(paste(Chi^2,'distance')),line=2.3,cex.lab=cexlbs)
box(bty='l',lwd=alwds)
# QQ Plot Assessing Multivariate Normality
abline(a=0,b=1,lty=3,col='gray40',lwd=lwds)
dev.copy(pdf,'Figures/MVTraitNormalQQplot.pdf',width=w,height=h)
dev.off()


# 1) All variables are reasonably Gaussian
cexaxs <- 0.9; cexlbs <- 1.4; cexmains <- 0.8; cexpts <- 1.2
lwds=2; alwds=2
graphics.off(); w=3.15;h=4.5;m1=1.5;m2=1.5;m3=1;m4=1; x11(width=w,height=h); 
par(mfrow=c(3,2),oma=c(2.6,2.6,0.2,0.2),mar=c(m1,m2,m3,m4))
for (j in 1:length(thesecols)) {
  qqnorm(Yf[,j],cex=cexpts,pch=19,bty='n',main='',xlab='',ylab='',
         cex.axis=cexaxs,cex.lab=cexlbs);
  qqline(Yf[,j],lty=3,col='gray40',lwd=lwds)
  if(j==3) title(ylab='Sample quantiles',line=1.3,cex.lab=cexlbs,outer=TRUE) 
  if(j==5) title(xlab='Theoretical quantiles',line=1.3,cex.lab=cexlbs,outer=TRUE)
  box(bty='l',lwd=alwds)
  mtext(c('Flower\\nsize','Style\\nlength','Pollen\\nzone',
          'Stigma-anther\\ndistance','Stigmatic\\nreceptivity',
          'Auto. selfing\\nability')[j],side=3,at=-2.7,line=-2,adj=0,cex=cexaxs-0.2)
}
# dev.copy(pdf,'Figures/TraitNormalQQplots.pdf',width=w,height=h)
# dev.off()

# 2) Variances are all reasonable homogeneous
graphics.off(); w=3.15;h=4;m1=1.5;m2=3.5;m3=1;m4=1; x11(width=w,height=h); 
par(mfrow=c(3,2),oma=c(2.6,0.2,0.2,0.2),mar=c(m1,m2,m3,m4))
for (j in 1:length(thesecols)) {
  plot(Yf[,j]~gdats[[4]][,'Population'],cex=cexpts,pch=19,bty='l',main='',
          xlab='',ylab='',cex.axis=cexaxs,cex.lab=cexlbs)
  title(ylab=c('Flower size','Style length','Pollen zone','Stig-ant. dist.',
               'Stigma. recep.','Auto. selfing')[j],line=2.3,cex.lab=cexlbs)
  if(j==5) title(xlab='Population',line=1.3,cex.lab=cexlbs,outer=TRUE)
}
dev.copy(pdf,'Figures/TraitBoxplots.pdf',width=w,height=h)
dev.off()

# MANOVA: The actual test
mf <- manova(Yf~Population, data=gdats[[4]][theserows,])
summary(mf)
summary.aov(mf) # for univariate results

# Set up Tukey HSD tests for trait-specific pairwise comparsions
f1 <- aov(Yf[,1]~Population,data=gdats[[4]]) # Flower size
fs.diff <- TukeyHSD(f1)
f2 <- aov(Yf[,2]~Population,data=gdats[[4]]) # Stigma length
sl.diff <- TukeyHSD(f2)
f3 <- aov(Yf[,3]~Population,data=gdats[[4]]) # Stamen traits
st.diff <- TukeyHSD(f3)
f4 <- aov(Yf[,4]~Population,data=gdats[[4]]) # Stigma-anther distance
sad.diff <- TukeyHSD(f4)
f5 <- aov(Yf[,5]~Population,data=gdats[[4]]) # Stigmatic receptivity
sr.diff <- TukeyHSD(f5)
f6 <- aov(Yf[,6]~Population,data=gdats[[4]]) # Autonomous selfing ability
as.diff <- TukeyHSD(f6)

# 1. Flower size is NOT different among populations.
fs.diff
# 2. Style length is NOT different among populations.
sl.diff
# 3. Stamen traits ARE different among populations: only BT > TMC.
st.diff
# 4. Stigma-anther distance is NOT different among populations.
sad.diff
# 5. Stigmatic receptivity IS different among populations: both BT and EF > TMC.
sr.diff
# 6. Autonomous selfing ability IS different among populations: both EF & BT > TMC.
as.diff
  
# Collect mean and standard errors from the ANOVAs
fs.est <- mSEpredict(f1, gdats[[4]]) 
sl.est <- mSEpredict(f2, gdats[[4]]) 
st.est <- mSEpredict(f3, gdats[[4]]) 
sad.est <- mSEpredict(f4, gdats[[4]])
sr.est <- mSEpredict(f5, gdats[[4]]) 
as.est <- mSEpredict(f6, gdats[[4]]) 

# -- | Figure | Plot means and standard errors for each stamen length at stage 4
cols='black' # Color palette
clab=1.5; pchs=19; pcex=2;  # Scaling factors, Point type
caxis=1.2; cxlab=1.3; clet=1.1; # Scaling
graphics.off(); w=7.5; h=4
x11(width=w,height=h); par(mfrow=c(2,3),oma=c(0.2,0.2,0.2,0.2),mar=c(2.3,4,2,1))  
# Panel a: flower size
siglets <- c('a','a','a'); yfac <- 1.01; yfac2 <- 1.003
cifun(fs.est[1,],fs.est[2,],xlabname='',ylabname='Flower size (PC1)',
      siglets,'a',cols,yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet)
# Panel c: Style length
siglets <- c('a','a','a'); yfac <- 1.01; yfac2 <- 1.003
cifun(sl.est[1,],sl.est[2,],xlabname='',ylabname='Style length (mm)',
      siglets,'c',cols,yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet)
# Panel e: stigmatic receptivity
siglets <- c('a','a','b'); yfac <- 1.05; yfac2 <- 1.01
cifun(sr.est[1,],sr.est[2,],xlabname='',ylabname='Stigmatic recep.',
      siglets,'e',cols,yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet)
# Panel b: Stamen traits
siglets <- c('a','ab','b'); yfac <- 1.01; yfac2 <- 1.003
cifun(st.est[1,],st.est[2,],xlabname='',ylabname='Pollen zone (PC1)',
      siglets,'b',cols,yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet)
# Panel d: Stigma-anther distance
siglets <- c('a','a','a'); yfac <- 1.01; yfac2 <- 1.003
plotCI(sad.est[1,],uiw=sad.est[2,],xlab='',ylab='',axes=T,
       xaxt='n',ylim=c(-0.07,0.2),col=cols,cex.lab=clab,cex.axis=caxis,
       sfrac=.012,lwd=1,gap=0.0,pch=pchs,cex=pcex, pt.bg=par("bg"),
       xlim=c(0.2,3.8),add=F,bty='n')
box(bty='l',lwd=alwds)
title(ylab='Stig.-anth. dist. (mm)',line=2.5,cex.lab=clab)
for(i in 1:3) { text(i,(sad.est[1,]*yfac2+sad.est[2,])[i],siglets[i],cex=clet,pos=3)}
axis(1,1:3,labels=c('BT','EF','TMC'),cex.axis=cxlab)
mtext('d', side=3, line=.25, at=0.1, outer=F, cex=clet, font=1)
# Panel f: automatic selfing ability
siglets <- c('a','a','b'); yfac <- 1.1;  yfac2 <- 1.01
cifun(as.est[1,],as.est[2,],xlabname='',ylabname='Auto.selfing rate',
      siglets,'f',cols,yfac,yfac2,clab,caxis,pchs,pcex,cxlab,clet)
# Save Figure
# dev.copy(pdf,'Figures/TraitComparisons_Fig1.pdf',width=w,height=h)
# dev.off()
# dev.copy(win.metafile,'Figures/TraitComparisons_Fig2.emf',width=5.9,height=7.8)
# dev.off()

#_______________________________________________________________________________
# -- | 7. Multivariate LMM using MCMCglmm | ------------------------------------
# Estimate Genetic correlations between the following traits
# to compare within and among genotype trait variation, estimate covariances 
# among clone means for each trait (except AS). Then calculate covariances among
# AS and other traits using geneotypic means, so AS can be included. Finally, 
# paste AS-trait covariances with other trait variance-covariance matrix to make 
# one comprehensive G matrix. For correlations, divide all covariances by the 
# square roots of the variance for each trait involved. 
# cdat contains the clonal means
summary (cdats)
# Center and standardize all data
ntraits <- 5 # FS, Sty, Sta, SAD, SR
thesecols <- c('GENOTYPE','Population','Flower.Size','STYLE_L','Stamen','SAD','SR')
thesecols <- c('GENOTYPE','Population','Flower.Size','STYLE_L','Stamen','aSAD','SR')
theserows <- which(rowSums((!is.na(cdats[[4]][,thesecols[-1:-2]]))*1)==5)
aggregate(cdats[[4]][theserows,],list(cdats[[4]][theserows,'Population']),
          mean,na.rm=T)
# Standardized traits not needed.

# The priors for the variance structures (R and G) are lists with the
# expected (co)variances (V) and degree of belief parameter (nu) for the inverse
# Wishart, and also the mean vector (alpha.mu) and covariance matrix (alpha.V)
# for the redundant working parameters. The defaults are nu=0, V=1, alpha.mu=0,
# and alpha.V=0. When alpha.V is non-zero, parameter expanded algorithms are
# used.

# Genetic Correlations | 
# Pc = Gc + Rc is not actually correct. See Roff 1997 p77.
# Pc = Gc*sqrt(hx2,hy2) + Rc*sqrt((1-hx2)(1-hy2))
# if the two heritabilites are very small, the phenotypic correlation is 
# determined primarily by the environmental correlation, whereas if they are 
# both high, it is the genetic correlation that is most important. Note that 
# the phenotypic corrleation does not by itself give any idea of the importance 
# of the genetic relationship between the two traits. (Roff 1997)

# Pearson product-moment correlations | Via (1984) in Roff 1997
# The estimation of genetic correlations by product-moment correlation between 
# family means allows an easy means for calculating confidence intervals. It is 
# an approximation because the variance and covariance terms contain a fraction 
# of the within-family error term: Cov_m = Cov_among + Cov_within / n
# n is family size.

# All variables are reasonably Gaussian
graphics.off()
Y <- cdats[[4]][,c('Flower.Size','STYLE_L','Stamen','SAD','SR')]
x11();par(mfrow=c(3,2),oma=c(1,1,1,1))
for (j in 1:ntraits) {
  qqnorm(Y[,j]);qqline(Y[,j])
}

# -----> 7A. Load MCMCglmm models and skip to line | --------------------------
load('Output/CvernaModelFits.Rdata')
# First fit a global model with all populations together. Evaluate convergence, 
# autocorrelation, and parameter identifiability. Determine the MCMC parameters 
# needed to obtain convergence and an effective sample size of >=1000.
burns <- 10000; thins <- 100; nitts <- 1.2*1000*thins + burns # for 1000 MCMC samples
nsam <- (nitts - burns) / thins; nsam; nitts
families <- rep("gaussian",ntraits); # trait means are gaussian
priors<-list(R=list(V=diag(ntraits), nu=1.002), # Slope and Intercept
              G=list(G1=list(V=diag(ntraits), nu=1.002)))
# start<-as.POSIXlt(Sys.time()) # Exceeds 4.0 minute ----------------------------|
# gmALL1 <- MCMCglmm(cbind(Flower.Size,STYLE_L,Stamen,SAD,SR) ~ trait - 1,
#                   random = ~us(trait):GENOTYPE, rcov = ~us(trait):units,
#                   data = cdats[[4]],family = families,
#                   verbose = FALSE,prior=priors,nitt=nitts,thin=thins,
#                   burnin=burns)
# end <-as.POSIXlt(Sys.time())
# duration1 = end-start; duration1 # Exceeds 4.0 minute -------------------------|
x11(width=8, height=8)
plot.estimates(gmALL1) # function defined above.

# Convergence diagnostics (copied from Titus von der Malsburg, MCMCglmm-intro on GitHub)
x11(8,16);par(mfrow=c(5,2), mar=c(2,2,1,0))
plot(gmALL1$Sol, auto.layout=F)

# In the trace, look for high autocorrelation, which is not good, and for 
# stationary traces.i.e., the sampling process remains in one part of the 
# parameter space

# Autocorrelation can be reduced through thinning. 
# Modify sample size and thinning together to obtain a large-enough set of usable samples.

autocorr.plot(gmALL1$Sol) # great drop off.
autocorr.plot(gmALL1$VCV) # great drop off
# VCV dimensions is 1000 rows (20K MCMC samples thinned by 20) and 50 columns
# Columns are for our G matrix (5 x 5) + R (5 x 5)
autocorr.diag(gmALL1$Sol)
autocorr.diag(gmALL1$VCV) # need to thin by more than 10 -> 20 is good
summary(effectiveSize(gmALL1$Sol)) # no reduction in ESS 
summary(effectiveSize(gmALL1$VCV)) # Reduction is 626

# Run the global model 3 times (3 chains) and check convergence diagnostics.
# start<-as.POSIXlt(Sys.time()) # Exceeds 8 minutes -----------------------------|
# gmALL2 <- MCMCglmm(cbind(Flower.Size,STYLE_L,Stamen,SAD,SR) ~ trait - 1,
#                    random = ~us(trait):GENOTYPE, rcov = ~us(trait):units,
#                    data = cdats[[4]],family = families,
#                    verbose = FALSE,prior=priors,nitt=nitts,thin=thins,
#                    burnin=burns)
# gmALL3 <- MCMCglmm(cbind(Flower.Size,STYLE_L,Stamen,SAD,SR) ~ trait - 1,
#                    random = ~us(trait):GENOTYPE, rcov = ~us(trait):units,
#                    data = cdats[[4]],family = families,
#                    verbose = FALSE,prior=priors,nitt=nitts,thin=thins,
#                    burnin=burns)
# end <-as.POSIXlt(Sys.time())
# end-start # Exceeds 8 minutes -------------------------------------------------|
gmALL <- lapply(list(gmALL1,gmALL2,gmALL3), function(m) m$Sol)
gmSol <- do.call(mcmc.list, gmALL)
gmALL <- lapply(list(gmALL1,gmALL2,gmALL3), function(m) m$VCV)
gmVCV <- do.call(mcmc.list, gmALL)
# The Gelman-Rubin test statistic is called the scale reduction factor. 
# The closer this factor is to 1, the better the convergence of our chains. 
# In practice, values below 1.1 can be acceptable and values below 1.02 are good. 

par(mfrow=c(5,2), mar=c(2,2,1,2))
gelman.plot(gmSol, auto.layout=F)
gelman.plot(gmVCV[,1:25], auto.layout=F)
# Compute scale reduction factors. 
gelman.diag(gmSol) # good
gelman.diag(gmVCV[,1:5]) # good
gelman.diag(gmVCV[,6:11]) # good
gelman.diag(gmVCV[,12:17]) # good
gelman.diag(gmVCV[,18:23]) # good
gelman.diag(gmVCV[,24:30]) # good
gelman.diag(gmVCV[,30:36]) # good
gelman.diag(gmVCV[,37:42]) # good
gelman.diag(gmVCV[,43:47]) # good
gelman.diag(gmVCV[,48:50]) # good

# Visually confirm adequate chain mixing. 
par(mfrow=c(5,2), mar=c(2,1,1,1))
plot(gmSol, ask=F, auto.layout=F)
plot(gmVCV, ask=F, auto.layout=F)
# Revisit Effect size | 3000 total MCMC samples
range(effectiveSize(gmSol)) #  ESS > 2909 /3000 = 97%
range(effectiveSize(gmVCV)) # ESS > 2750 /3000 = 86%
# save(gmALL1,gmALL2,gmALL3,file='Output/CvernaModelFits.Rdata')

gmALL1$DIC # DIC = 5260.169
summary(gmALL1) 
posterior.mode(gmALL1$Sol)
posterior.mode(gmALL1$VCV)[1:(ntraits*ntraits)]
HPDinterval(gmALL1$VCV)

x11(width=5, height=5)
plot.estimates(gmALL1) # function defined above.
ALLest <- posterior.mode(gmALL1$Sol)
#_______________________________________________________________________________

#_______________________________________________________________________________ 
# -----> 7B. Analyze by Population | -------------------------------------------
# Raw Data for Vg calculation
# start<-as.POSIXlt(Sys.time()) # Exceeds 4 minutes -----------------------------|
# gmBT <- MCMCglmm(cbind(Flower.Size,STYLE_L,Stamen,SAD,SR) ~ trait - 1,
#                  random = ~us(trait):GENOTYPE, rcov = ~us(trait):units,
#                  data = subset(cdats[[4]],Population=='BT'),family = families,
#                  verbose = FALSE,prior=priors,nitt=nitts,thin=thins,
#                  burnin=burns)
# gmEF <- MCMCglmm(cbind(Flower.Size,STYLE_L,Stamen,SAD,SR) ~ trait - 1,
#                  random = ~us(trait):GENOTYPE, rcov = ~us(trait):units,
#                  data = subset(cdats[[4]],Population=='EF'),family = families,
#                  verbose = FALSE,prior=priors,nitt=nitts,thin=thins,
#                  burnin=burns)
# gmTMC <- MCMCglmm(cbind(Flower.Size,STYLE_L,Stamen,SAD,SR) ~ trait - 1,
#                   random = ~us(trait):GENOTYPE, rcov = ~us(trait):units,
#                   data = subset(cdats[[4]],Population=='TMC'),family = families,
#                   verbose = FALSE,prior=priors,nitt=nitts,thin=thins,
#                   burnin=burns)
# end <-as.POSIXlt(Sys.time())
# end-start # Exceeds 4 minutes -------------------------------------------------|
range(effectiveSize(gmBT$Sol)) 
range(effectiveSize(gmEF$Sol))
range(effectiveSize(gmTMC$Sol))
range(effectiveSize(gmBT$VCV))
range(effectiveSize(gmEF$VCV))
range(effectiveSize(gmTMC$VCV)) #  ESS 695 /1000 = 69%
# save(gmBT,gmEF,gmTMC,gmALL1,gmALL2,gmALL3,file='Output/CvernaModelFits.Rdata')
## For calculation of 95% CI on correlations. Calculate correlations for each 
# MCMC sample, then compute quantiles.
nsam
dim(gmBT$VCV)
GvBTs <- GvEFs <- GvTMCs <- GcBTs <- GcEFs <- GcTMCs <- array(NA,dim=c(5,5,nsam))
HIBTs <- HIEFs <- HITMCs <- array(NA,dim=c(5,2,nsam))
for (i in 1:nsam) {
  # Genetic variance-covariance matrix
  GvBT <- matrix(gmBT$VCV[i,1:25],nrow=ntraits)
  GvEF <- matrix(gmEF$VCV[i,1:25],nrow=ntraits)
  GvTMC <- matrix(gmTMC$VCV[i,1:25],nrow=ntraits)
  # Residual variance-covariance
  RvBT <- matrix(gmBT$VCV[i,26:50],nrow=ntraits)
  RvEF <- matrix(gmEF$VCV[i,26:50],nrow=ntraits)
  RvTMC <- matrix(gmTMC$VCV[i,26:50],nrow=ntraits)
  # Total Phenotypic variance-covariance
  PvBT <- GvBT + RvBT;
  PvEF <- GvEF + RvEF;
  PvTMC <- GvTMC + RvTMC; 
  # Calculate correlations | r = CovXY / sqrt(VX VY)
  GcBT <- GvBT/sqrt(diag(ntraits)%*%diag(GvBT)%*%diag(GvBT))
  GcEF <- GvEF/sqrt(diag(ntraits)%*%diag(GvEF)%*%diag(GvEF))
  GcTMC <- GvTMC/sqrt(diag(ntraits)%*%diag(GvTMC)%*%diag(GvTMC))
  # Trait values
  ZBT <- matrix(gmBT$Sol[i,],nrow=ntraits,ncol=1) 
  ZEF <- matrix(gmEF$Sol[i,],nrow=ntraits,ncol=1) 
  ZTMC <- matrix(gmTMC$Sol[i,],nrow=ntraits,ncol=1) 
  # Calculate heritability and evolvability
  HIBT <- cbind('H2'=diag(GvBT)/diag(PvBT),'I'=diag(GvBT)/ZBT^2)
  HIEF <- cbind('H2'=diag(GvEF)/diag(PvEF),'I'=diag(GvEF)/ZEF^2)
  HITMC <- cbind('H2'=diag(GvTMC)/diag(PvTMC),'I'=diag(GvTMC)/ZTMC^2)
  # Store each VCV and Cor mat
  GvBTs[,,i] <- GvBT
  GvEFs[,,i] <- GvEF
  GvTMCs[,,i] <- GvTMC
  GcBTs[,,i] <- GcBT
  GcEFs[,,i] <- GcEF
  GcTMCs[,,i] <- GcTMC
  HIBTs[,,i] <- HIBT
  HIEFs[,,i] <- HIEF
  HITMCs[,,i] <- HITMC
  rm(GvBT,GvEF,GvTMC,RvBT,RvEF,RvTMC,PvBT,PvEF,PvTMC,GcBT,GcEF,GcTMC,ZBT,ZEF,
     ZTMC,HIBT,HIEF,HITMC)  
}
GvBTe <- GvEFe <- GvTMCe <- GcBTe <- GcEFe <- GcTMCe <- array(NA,dim=c(5,5,2))
HIBTe <- HIEFe <- HITMCe <- array(NA,dim=c(5,2,2))
for (i in 1:5) {
  for (j in 1:5) {
    GvBTe[i,j,] <- HPDinterval(as.mcmc.list(as.mcmc(GvBTs[i,j,])))[[1]][1:2]
    GvEFe[i,j,] <- HPDinterval(as.mcmc.list(as.mcmc(GvEFs[i,j,])))[[1]][1:2]
    GvTMCe[i,j,] <- HPDinterval(as.mcmc.list(as.mcmc(GvTMCs[i,j,])))[[1]][1:2]
    GcBTe[i,j,] <- HPDinterval(as.mcmc.list(as.mcmc(GcBTs[i,j,])))[[1]][1:2]
    GcEFe[i,j,] <- HPDinterval(as.mcmc.list(as.mcmc(GcEFs[i,j,])))[[1]][1:2]
    GcTMCe[i,j,] <- HPDinterval(as.mcmc.list(as.mcmc(GcTMCs[i,j,])))[[1]][1:2]
  }
  for (k in 1:2) {
    HIBTe[i,k,] <- HPDinterval(as.mcmc.list(as.mcmc(HIBTs[i,k,])))[[1]][1:2]
    HIEFe[i,k,] <- HPDinterval(as.mcmc.list(as.mcmc(HIEFs[i,k,])))[[1]][1:2]
    HITMCe[i,k,] <- HPDinterval(as.mcmc.list(as.mcmc(HITMCs[i,k,])))[[1]][1:2]
  }
}
HIBTe
# Collect the posterior modes for each variance
GvBT <- matrix(posterior.mode(gmBT$VCV[,1:25]),nrow=ntraits,ncol=ntraits)
GvEF <- matrix(posterior.mode(gmEF$VCV[,1:25]),nrow=ntraits,ncol=ntraits)
GvTMC <- matrix(posterior.mode(gmTMC$VCV[,1:25]),nrow=ntraits,ncol=ntraits)
# Residual variance-covariance
RvBT <- matrix(posterior.mode(gmBT$VCV[,26:50]),nrow=ntraits,ncol=ntraits)
RvEF <- matrix(posterior.mode(gmEF$VCV[,26:50]),nrow=ntraits,ncol=ntraits)
RvTMC <- matrix(posterior.mode(gmTMC$VCV[,26:50]),nrow=ntraits,ncol=ntraits)
# Total Phenotypic variance-covariance
PvBT <- GvBT + RvBT;
PvEF <- GvEF + RvEF;
PvTMC <- GvTMC + RvTMC; 
# Compute correlations
GcBT <- GvBT/sqrt(diag(ntraits)%*%diag(GvBT)%*%diag(GvBT))
GcEF <- GvEF/sqrt(diag(ntraits)%*%diag(GvEF)%*%diag(GvEF))
GcTMC <- GvTMC/sqrt(diag(ntraits)%*%diag(GvTMC)%*%diag(GvTMC))
# Trait values
ZBT <- posterior.mode(gmBT$Sol) 
ZEF <- posterior.mode(gmEF$Sol) 
ZTMC <- posterior.mode(gmTMC$Sol)
# Compute heritability, and evolvability using posterior modes   
HIBT <- cbind('H2'=diag(GvBT)/diag(PvBT),'I'=diag(GvBT)/ZBT^2)
HIEF <- cbind('H2'=diag(GvEF)/diag(PvEF),'I'=diag(GvEF)/ZEF^2)
HITMC <- cbind('H2'=diag(GvTMC)/diag(PvTMC),'I'=diag(GvTMC)/ZTMC^2)

# Different from one or negative one? Note: ones on the diagonal are expected
which(GcBTe[,,1]<=-1);which(GcBTe[,,2]>=1) # answer = integer(0) [1] 1 7 13 19 25
which(GcEFe[,,1]<=-1);which(GcEFe[,,2]>=1) # answer = integer(0) [1] 1 7 13 19 25
which(GcTMCe[,,1]<=-1);which(GcTMCe[,,2]>=1) # answer = integer(0) [1] 1 7 13 19 25
which(HIBTe[,1,2]>=1) # answer = integer(0) 
which(HIEFe[,1,2]>=1) # answer = integer(0) 
which(HITMCe[,1,2]>=1) # answer = integer(0) 
# --> No perfect correlations. AND No perfectly heritable traits.

# Different from zero?? uses sigdif defined above.
# Correlations
GcBT <- round(GcBT,3)
GcBT[sigdif(cbind(c(GcBTe[,,1]),c(GcBTe[,,2])))] <- paste(GcBT[sigdif(cbind(
  c(GcBTe[,,1]),c(GcBTe[,,2])))],'*',sep='')
GcEF <- round(GcEF,3)
GcEF[sigdif(cbind(c(GcEFe[,,1]),c(GcEFe[,,2])))] <- paste(GcEF[sigdif(cbind(
  c(GcEFe[,,1]),c(GcEFe[,,2])))],'*',sep='')
GcTMC <- round(GcTMC,3)
GcTMC[sigdif(cbind(c(GcTMCe[,,1]),c(GcTMCe[,,2])))] <- paste(GcTMC[sigdif(cbind(
  c(GcTMCe[,,1]),c(GcTMCe[,,2])))],'*',sep='')
# Covariances
GBT <- round(GvBT,3)
GBT[sigdif(cbind(c(GvBTe[,,1]),c(GvBTe[,,2])))] <- paste(GBT[sigdif(cbind(
  c(GvBTe[,,1]),c(GvBTe[,,2])))],'*',sep='')
GEF <- round(GvEF,3)
GEF[sigdif(cbind(c(GvEFe[,,1]),c(GvEFe[,,2])))] <- paste(GEF[sigdif(cbind(
  c(GvEFe[,,1]),c(GvEFe[,,2])))],'*',sep='')
GTMC <- round(GvTMC,3)
GTMC[sigdif(cbind(c(GvTMCe[,,1]),c(GvTMCe[,,2])))] <- paste(GTMC[sigdif(cbind(
  c(GvTMCe[,,1]),c(GvTMCe[,,2])))],'*',sep='')
# Heritability and Evolvability
# Compile the matrices with Heritability and Evolvability
GBT[upper.tri(GBT, diag = FALSE)] <- GcBT[upper.tri(GcBT, diag = FALSE)]
GBT <- cbind(round(HIBT ,3),GBT)
GEF[upper.tri(GEF, diag = FALSE)] <- GcEF[upper.tri(GcEF, diag = FALSE)]
GEF <- cbind(round(HIEF,3),GEF)
GTMC[upper.tri(GTMC, diag = FALSE)] <- GcTMC[upper.tri(GcTMC, diag = FALSE)]
GTMC <- cbind(round(HITMC,3),GTMC)
GBT

Zs <- array(NA,dim=c(ntraits+1,3,3)); 
dimnames(Zs) <- list(c('Flower.Size','STYLE_L','Stamen','SAD','SR','AS'),
                     c('BT','EF','TMC'),c('mode','lHPD','uHPD'))
Zs[1:5,,1] <- cbind(posterior.mode(gmBT$Sol),posterior.mode(gmEF$Sol),
                    posterior.mode(gmTMC$Sol))
Zs[1:5,,2] <- cbind(HPDinterval(gmBT$Sol),HPDinterval(gmEF$Sol),
                    HPDinterval(gmTMC$Sol))[,c(1,3,5)]
Zs[1:5,,3] <- cbind(HPDinterval(gmBT$Sol),HPDinterval(gmEF$Sol),
                    HPDinterval(gmTMC$Sol))[,c(2,4,6)]
dimnames(GBT) <- list(c('Flower','Style','Stamen','SAD','SR'),
                      c('H2','I','Flower','Style','Stamen','SAD','SR'))
dimnames(GEF) <- dimnames(GTMC) <- dimnames(GBT)
# Add a row and column for AS
GBT <- cbind(GBT,'AS' = NA)
GEF <- cbind(GEF,'AS' = NA)
GTMC <- cbind(GTMC,'AS' = NA)
GBT <- rbind(GBT,'AS' = NA)
GEF <- rbind(GEF,'AS' = NA)
GTMC <- rbind(GTMC,'AS' = NA)
#_______________________________________________________________________________

#_______________________________________________________________________________
# -----> 7C. Add in Automatic Selfing Ability (AS) | ---------------------------
# Different clones were used for floral traits and selfing ability, which 
# complicates multi trait genetic covariance calculations.  
# I can't look at within versus between variance components, but I can look at 
# mean genotype variance-covariance structure. 
# used for selfing ability.

# Add one to priors and distributions for addition of  AS
priorsa<-list(R=list(V=diag(ntraits+1)*1, nu=0.002))
families <- rep("gaussian",ntraits+1); 
# gdats is the data set of genotypic means
# start<-as.POSIXlt(Sys.time()) # Less than 1 minute ----------------------------|
# gmBTf <- MCMCglmm(cbind(Flower.Size,STYLE_L,Stamen,SAD,SR,AS) ~ trait - 1, 
#                    rcov = ~us(trait):units,family = families, verbose = FALSE,
#                    data = subset(gdats[[4]],Population=='BT'), prior=priorsa,
#                    nitt=nitts,thin=thins, burnin=burns)
# gmEFf <- MCMCglmm(cbind(Flower.Size,STYLE_L,Stamen,SAD,SR,AS) ~ trait - 1, 
#                    rcov = ~us(trait):units,family = families, verbose = FALSE,
#                    data = subset(gdats[[4]],Population=='EF'), prior=priorsa,
#                    nitt=nitts,thin=thins, burnin=burns)
# gmTMCf <- MCMCglmm(cbind(Flower.Size,STYLE_L,Stamen,SAD,SR,AS) ~ trait - 1, 
#                     rcov = ~us(trait):units,family = families, verbose = FALSE,
#                     data = subset(gdats[[4]],Population=='TMC'), prior=priorsa,
#                     nitt=nitts,thin=thins, burnin=burns)
# end <-as.POSIXlt(Sys.time())
# end-start # Less than 1 minute ------------------------------------------------|
# Raw Data for Vg calculation
range(effectiveSize(gmBTf$Sol)) 
range(effectiveSize(gmEFf$Sol))
range(effectiveSize(gmTMCf$Sol))
range(effectiveSize(gmBTf$VCV))
range(effectiveSize(gmEFf$VCV))
range(effectiveSize(gmTMCf$VCV)) #  ESS 695 /1000 = 69%
# save(gmBT,gmEF,gmTMC,gmBTf,gmEFf,gmTMCf,gmALL1,gmALL2,gmALL3,
#      file='Output/CvernaModelFits.Rdata')
## For calculation of 95% CI on correlations. Calculate correlations for each 
# MCMC sample, then compute quantiles.
nsam
dim(gmBTf$VCV) # 36 rows (6 x 6 traits)
GvBTfs <- GvEFfs <- GvTMCfs <- GcBTfs <- array(NA,dim=c(6,6,nsam))
GcEFfs <- GcTMCfs <- array(NA,dim=c(6,6,nsam))
for (i in 1:nsam) {
  # Variance-covariance matrix
  GvBTf <- matrix(gmBTf$VCV[i,],nrow=ntraits+1)
  GvEFf <- matrix(gmEFf$VCV[i,],nrow=ntraits+1)
  GvTMCf <- matrix(gmTMCf$VCV[i,],nrow=ntraits+1)
  # Calculate correlations | r = CovXY / sqrt(VX VY)
  GcBTf <- GvBTf/sqrt(diag(ntraits+1)%*%diag(GvBTf)%*%diag(GvBTf))
  GcEFf <- GvEFf/sqrt(diag(ntraits+1)%*%diag(GvEFf)%*%diag(GvEFf))
  GcTMCf <- GvTMCf/sqrt(diag(ntraits+1)%*%diag(GvTMCf)%*%diag(GvTMCf))
  # Store each VCV and Cor mat
  GvBTfs[,,i] <- GvBTf
  GvEFfs[,,i] <- GvEFf
  GvTMCfs[,,i] <- GvTMCf
  GcBTfs[,,i] <- GcBTf
  GcEFfs[,,i] <- GcEFf
  GcTMCfs[,,i] <- GcTMCf
  rm(GvBTf,GvEFf,GvTMCf,GcBTf,GcEFf,GcTMCf)  
}
GvBTfe <- GvEFfe <- GvTMCfe <- GcBTfe <- GcEFfe <- GcTMCfe <- array(NA,dim=c(6,6,2))
for (i in 1:6) {
  for (j in 1:6) {
    GvBTfe[i,j,] <- HPDinterval(as.mcmc.list(as.mcmc(GvBTfs[i,j,])))[[1]][1:2]
    GvEFfe[i,j,] <- HPDinterval(as.mcmc.list(as.mcmc(GvEFfs[i,j,])))[[1]][1:2]
    GvTMCfe[i,j,] <- HPDinterval(as.mcmc.list(as.mcmc(GvTMCfs[i,j,])))[[1]][1:2]
    GcBTfe[i,j,] <- HPDinterval(as.mcmc.list(as.mcmc(GcBTfs[i,j,])))[[1]][1:2]
    GcEFfe[i,j,] <- HPDinterval(as.mcmc.list(as.mcmc(GcEFfs[i,j,])))[[1]][1:2]
    GcTMCfe[i,j,] <- HPDinterval(as.mcmc.list(as.mcmc(GcTMCfs[i,j,])))[[1]][1:2]
  }
}
# Posterior modes
GvBTf <- matrix(posterior.mode(gmBTf$VCV[,]),nrow=ntraits+1,ncol=ntraits+1)
GvEFf <- matrix(posterior.mode(gmEFf$VCV[,]),nrow=ntraits+1,ncol=ntraits+1)
GvTMCf <- matrix(posterior.mode(gmTMCf$VCV[,]),nrow=ntraits+1,ncol=ntraits+1)
# Compute correlations
GcBTf <- GvBTf/sqrt(diag(ntraits+1)%*%diag(GvBTf)%*%diag(GvBTf))
GcEFf <- GvEFf/sqrt(diag(ntraits+1)%*%diag(GvEFf)%*%diag(GvEFf))
GcTMCf <- GvTMCf/sqrt(diag(ntraits+1)%*%diag(GvTMCf)%*%diag(GvTMCf))

# Different from one or negative one? Note: ones on the diagonal are expected
which(GcBTfe[,,1]<=-1);which(GcBTfe[,,2]>=1) # answer = integer(0) [1] 1 8 15 22 29 36
which(GcEFfe[,,1]<=-1);which(GcEFfe[,,2]>=1) # answer = integer(0) [1] 1 8 15 22 29 36
which(GcTMCfe[,,1]<=-1);which(GcTMCfe[,,2]>=1) # answer = integer(0) [1] 1 8 15 22 29 36
# --> No perfect correlations. AND No perfectly heritable traits.

# Different from zero?? uses sigdif defined above.
# Correlations
GcBTas <- round(GcBTf,3)
GcBTas[sigdif(cbind(c(GcBTfe[,,1]),c(GcBTfe[,,2])))] <- paste(GcBTas[sigdif(cbind(
  c(GcBTfe[,,1]),c(GcBTfe[,,2])))],'*',sep='')
GcEFas <- round(GcEFf,3)
GcEFas[sigdif(cbind(c(GcEFfe[,,1]),c(GcEFfe[,,2])))] <- paste(GcEFas[sigdif(cbind(
  c(GcEFfe[,,1]),c(GcEFfe[,,2])))],'*',sep='')
GcTMCas <- round(GcTMCf,3)
GcTMCas[sigdif(cbind(c(GcTMCfe[,,1]),c(GcTMCfe[,,2])))] <- paste(GcTMCas[sigdif(cbind(
  c(GcTMCfe[,,1]),c(GcTMCfe[,,2])))],'*',sep='')
# Covariances
GBTas <- round(GvBTf,3)
GBTas[sigdif(cbind(c(GvBTfe[,,1]),c(GvBTfe[,,2])))] <- paste(GBTas[sigdif(cbind(
  c(GvBTfe[,,1]),c(GvBTfe[,,2])))],'*',sep='')
GEFas <- round(GvEFf,3)
GEFas[sigdif(cbind(c(GvEFfe[,,1]),c(GvEFfe[,,2])))] <- paste(GEFas[sigdif(cbind(
  c(GvEFfe[,,1]),c(GvEFfe[,,2])))],'*',sep='')
GTMCas <- round(GvTMCf,3)
GTMCas[sigdif(cbind(c(GvTMCfe[,,1]),c(GvTMCfe[,,2])))] <- paste(GTMCas[sigdif(cbind(
  c(GvTMCfe[,,1]),c(GvTMCfe[,,2])))],'*',sep='')

GBTas
GEFas
GTMCas
# None of the covariances are significant

# Decompose AS variance
asBT <- MCMCglmm(AA ~ 1,random = ~GTYPE, 
                 rcov = ~units, data = subset(asdat,POP=='BT'),
                 family = "gaussian",verbose = FALSE,nitt=nitts,thin=thins, 
                 burnin=burns)
asEF <- MCMCglmm(AA ~ 1,random = ~GTYPE, 
                 rcov = ~units, data = subset(asdat,POP=='EF'),
                 family = "gaussian",verbose = FALSE,nitt=nitts,thin=thins, 
                 burnin=burns)
asTMC <- MCMCglmm(AA ~ 1,random = ~GTYPE, 
                  rcov = ~units, data = subset(asdat,POP=='TMC'),
                  family = "gaussian",verbose = FALSE,nitt=nitts,thin=thins, 
                  burnin=burns)

posterior.mode(asBT$VCV[,1])
HPDinterval(asBT$VCV[,1]) # Vg for AS is significantly different from zero
HPDinterval(asBT$VCV[,1] / (asBT$VCV[,1]+asBT$VCV[,2])) # H2 also significant
HPDinterval(asEF$VCV[,1]) # Vg for AS is significantly different from zero
HPDinterval(asEF$VCV[,1] / (asEF$VCV[,1]+asEF$VCV[,2])) # H2 also significant
HPDinterval(asTMC$VCV[,1]) # Vg for AS is significantly different from zero
HPDinterval(asTMC$VCV[,1] / (asTMC$VCV[,1]+asTMC$VCV[,2])) # H2 also significant
GvBTas <- posterior.mode(asBT$VCV[,1])
RvBTas <- posterior.mode(asBT$VCV[,2])
GvEFas <- posterior.mode(asEF$VCV[,1])
RvEFas <- posterior.mode(asEF$VCV[,2])
GvTMCas <- posterior.mode(asTMC$VCV[,1])
RvTMCas <- posterior.mode(asTMC$VCV[,2])
PvBTas <- GvBTas + RvBTas;
PvEFas <- GvEFas + RvEFas;
PvTMCas <- GvTMCas + RvTMCas;
PvTMCas; GvTMCas; RvTMCas;
var(asdat$AA) # Total variation matches. 
GvTMCas / PvTMCas # heritability 
  
# Calculate correlations | r = CovXY / sqrt(VX VY)
GvBTf[,6]/sqrt(c(diag(GvBT),GvBTas[1])*GvBTas[1]) # This overestimates Gen Corr and is inappropriate.
GvBTf[,6]/sqrt(diag(GvBTf)*diag(GvBTf)[6]) # these are correlations among genotypes, but perhaps is less biased than the calculation above.
GcBTf [,6] 
# Equal to correlations of family means
newTraits <- matrix(NA,nrow=ntraits+1,ncol=6)
ASvar <- matrix(NA,nrow=3,ncol=4)
dimnames(newTraits) <- list(c('FS','SL','Po','SAD','SR','AS'),
                            c('Bv','Ev','Tv','Bc','Ec','Tc'))
dimnames(ASvar) <- list(c('B','E','T'),c('Gv','Pv','H2','I'))
newTraits[,1] <- GvBTf[,6]
newTraits[,2] <- GvEFf[,6]
newTraits[,3] <- GvTMCf[,6]
newTraits[,4] <- GcBTf [,6]
newTraits[,5] <- GcEFf [,6]
newTraits[,6] <- GcTMCf [,6]
ASvar[,1] <- c(GvBTas[1],GvEFas[1],GvTMCas[1])
ASvar[,2] <- c(PvBTas[1],PvEFas[1],PvTMCas[1])
ASvar[,3] <- ASvar[,1] / ASvar[,2]
ASmeans <- c(posterior.mode(asBT$Sol),posterior.mode(asEF$Sol),
             posterior.mode(asTMC$Sol))
ASvar[,4] <- ASvar[,1] / ASmeans^2
# All three Gv's are significantly different from zero.
# None of the G cov's are significantly different from zero.
# None of the G cor's are significantly different from zero.

# Add AS info to previously compiled matrices
GBT[6,3:8] <- round(newTraits[,1],3)
GEF[6,3:8] <- round(newTraits[,2],3)
GTMC[6,3:8] <- round(newTraits[,3],3)
GBT[-6,8] <- round(newTraits[-6,4],3)
GEF[-6,8] <- round(newTraits[-6,5],3)
GTMC[-6,8] <- round(newTraits[-6,6],3)
GBT[6,1:2] <- round(ASvar[1,3:4],3)
GEF[6,1:2] <- round(ASvar[2,3:4],3)
GTMC[6,1:2] <- round(ASvar[3,3:4],3)
GBT[6,8] <- paste(round(ASvar[1,1],3),'*',sep='')
GEF[6,8] <- paste(round(ASvar[2,1],3),'*',sep='')
GTMC[6,8] <- paste(round(ASvar[3,1],3),'*',sep='')
GBT; GEF; GTMC
# -----> 7D. Save G matrices | -------------------------------------------------
# Save genetic variance/covariance matrices 
dimnames(GvBT) <- list(c('Flower','Style','Stamen','SAD','SR'),
                       c('Flower','Style','Stamen','SAD','SR'))
dimnames(GvEF) <- dimnames(GvTMC) <- dimnames(GvBT)
BTcov <- rbind(cbind(GvBT,'AS'=newTraits[1:5,1]),'AS'=newTraits[,1])
EFcov <- rbind(cbind(GvEF,'AS'=newTraits[1:5,2]),'AS'=newTraits[,2])
TMCcov <- rbind(cbind(GvTMC,'AS'=newTraits[1:5,3]),'AS'=newTraits[,3])
BTcov[6,6] <- ASvar[1,1]
EFcov[6,6] <- ASvar[2,1]
TMCcov[6,6] <- ASvar[3,1]
Zs[6,,1] <- ASmeans 
Zs[6,,2] <- c(HPDinterval(asBT$VCV[,1])[[1]],HPDinterval(asEF$VCV[,1])[[1]],
              HPDinterval(asTMC$VCV[,1])[[1]])
Zs[6,,3] <- c(HPDinterval(asBT$VCV[,1])[[2]],HPDinterval(asEF$VCV[,1])[[2]],
              HPDinterval(asTMC$VCV[,1])[[2]])
save(BTcov,EFcov,TMCcov,Zs,file='Output/CvernaGmatrices.Rdata')
# Sink a text file with all of the main findings
sink('Output/Spigler_etal_GenVarSelfing_Results_Output.txt')
print('Spigler_etal_GenVarSelfing_Results_Output')
print('')
print('Sample size breakdown')
Ntab
print('')
print('Correlations among genotypic mean values')
print('BT')
BTbigCor
print('EF')
EFbigCor
print('TMC')
TMCbigCor
print('')
print('MANOVA | Petal trait differences among populations (Stage 4 only)')
summary.aov(mp) # for univariate results
print('')
print('Upper petals ARE different among populations: both EF & BT > TMC.')
up.diff
print('Lower petals ARE different among populations: only BT > TMC.')
lp.diff
print('Keel petal length is NOT different among populations.')
kl.diff
print('Keel petal area also is NOT different among populations.')
ka.diff
print('')
print('MANOVA | Stamen length differences among populations (Stage 4 only)')
summary.aov(ms) # for univariate results
print('')
print('Stamen 1 length IS different among populations: only BT > TMC.')
s1.diff
print('Stamen 2 length IS different among populations: only BT > TMC.')
s2.diff
print('Stamen 3 length IS different among populations: only BT > TMC.')
s3.diff
print('Stamen 4 length IS different among populations: only BT > TMC.')
s4.diff
print('')
print('MANOVA | Floral trait differences among populations (Stage 4 only)')
summary(mf) # for multivariate results
summary.aov(mf) # for univariate results
print('')
print('Flower size is NOT different among populations.')
fs.diff
print('Style length is NOT different among populations.')
sl.diff
print('Stamen traits ARE different among populations: only BT > TMC.')
st.diff
print('Stigma-anther distance is NOT different among populations.')
sad.diff
print('Stigmatic receptivity IS different among populations: both BT and EF > TMC.')
sr.diff
print('Autonomous selfing ability IS different among populations: both EF & BT > TMC.')
as.diff
print('')
print('MCMC sampling parameters')
print('Burn in')
burns
print('Thinning')
thins
print('Number of iterations')
nitts
print('Number of MCMC samples analyzed')
nsam
print('')
print('Summaries for Genetic variance-covariance')
print('BT G and Genetic Correlation table')
GBT
print('EF G and Genetic Correlation table')
GEF
print('TMC G and Genetic Correlation table')
GTMC
print('')
print('G matrices')
print('BT G matrix')
BTcov
print('EF G matrix')
EFcov
print('TMC G matrix')
TMCcov
print('')
print('Trait values: Posterior mode and HPD interval')
Zs

sink()
#_______________________________________________________________________________
#  END SCRIPT ------------------------------------------------------------------
#_______________________________________________________________________________