library(maptools)
library(spatstat)
library(spatialEco)

#CLARK-EVANS Analysis----

#LOAD DATA - x,y coordinates
field<-read.csv("...Demographics_XYcoord.csv")

str(field)
onlyT<-subset(field, Species=='Tlong')


#The Clark and Evans (1954) aggregation index R is a crude measure of clustering or ordering of a
#point pattern. It is the ratio of the observed mean nearest neighbour distance in the pattern to that
#expected for a Poisson point process of the same intensity. A value R > 1 suggests ordering, while
#R < 1 suggests clustering.

window<-disc(radius=3, centre=c(0,0)) #designate size of plot

#Make data sets for each plot
A<-subset(onlyT, PlotID=='A')
A<-as.ppp(A, window)
clarkevans.test(A, alternative="clustered")
mean(nndist(A))
clarkevans.test(dataT,marks=PlotID, alternative="clustered")
B<-subset(onlyT, PlotID=='B')
B<-as.ppp(B, window)
clarkevans.test(B, alternative="clustered")
C<-subset(onlyT, PlotID=='C')
C<-as.ppp(C, window)
clarkevans.test(C, alternative="clustered")
D<-subset(onlyT, PlotID=='D')
D<-as.ppp(D, window)
clarkevans.test(D, alternative="clustered")
E<-subset(onlyT, PlotID=='E')
E<-as.ppp(E, window)
clarkevans.test(E, alternative="clustered")
F<-subset(onlyT, PlotID=='F')
F<-as.ppp(F, window)
clarkevans.test(F, alternative="clustered")
BM1<-subset(onlyT, PlotID=='BM1')
BM1<-as.ppp(BM1, window)
clarkevans.test(BM1, alternative="clustered")
BM2<-subset(onlyT, PlotID=='BM2')
BM2<-as.ppp(BM2, window)
clarkevans.test(BM2, alternative="clustered")
BR<-subset(onlyT, PlotID=='BR')
BR<-as.ppp(BR, window)
clarkevans.test(BR, alternative="clustered")

#Determine the nearest neighbor distance to a Pam - closest Tlong
Ap<-subset(field,PlotID=='A')
Ap<-as.ppp(Ap,window)
Ap$nnA<-nndist.ppp(Ap, method="C") #Pam nndist = 1.155
plot(Ap)
Bp<-subset(field,PlotID=='B')
Bp<-as.ppp(Bp,window)
Bp$nnB<-nndist.ppp(Bp, method="C") #Pam nndist = 1.086
Bp$nnB
plot(Bp)
BM1p<-subset(field,PlotID=='BM1')
BM1p<-as.ppp(BM1p,window)
BM1p$nnB<-nndist.ppp(BM1p, method="C") #Pam nndist = 0.593
BM1p$nnB

BM2p<-subset(field,PlotID=='BM2')
BM2p<-as.ppp(BM2p,window)
BM2p$nnB<-nndist.ppp(BM2p, method="C") #Pam nndist = 1.083
BM2p$nnB

BRp<-subset(field,PlotID=='BR')
BRp<-as.ppp(BRp,window)
BRp$nnB<-nndist.ppp(BRp, method="C") #Pam nndist = 0.822
BRp$nnB

Cp<-subset(field,PlotID=='C')
Cp<-as.ppp(Cp,window)
Cp$nnB<-nndist.ppp(Cp, method="C") #Pam nndist = 0.719
Cp$nnB

Dp<-subset(field,PlotID=='D')
Dp<-as.ppp(Dp,window)
Dp$nnB<-nndist.ppp(Dp, method="C") #Pam nndist = 0.868
Dp$nnB

Ep<-subset(field,PlotID=='E'& NestID !="E2") #Removed the other Pam colony, which was the clostest neighbor to the Pam; only considering closest Tlong colony
Ep<-as.ppp(Ep,window)
Ep$nnB<-nndist.ppp(Ep, method="C") #Pam nndist = 0.778
Ep$nnB
plot(Ep)

Fp<-subset(field,PlotID=='F')
Fp<-as.ppp(Fp,window)
Fp$nnB<-nndist.ppp(Fp, method="C") #Pam nndist =0.927
Fp$nnB

PamNN<-c(1.155,1.086,0.593,1.083,0.822,0.719,0.868,0.778,0.927)
summary(PamNN) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5930  0.7780  0.8680  0.8923  1.0830  1.1550 
var(PamNN)
sd(PamNN)

#LEVENE'S TEST----
library(car)

#LOAD DATA
FieldN<-read.csv("...Demographics - Total.csv")

#Create different datasubsets for each plot
A<-FieldN[which(FieldN$PlotID=='A'), ]
B<-FieldN[which(FieldN$PlotID=='B'), ]
C<-FieldN[which(FieldN$PlotID=='C'), ]
D<-FieldN[which(FieldN$PlotID=='D'), ]
E<-FieldN[which(FieldN$PlotID=='E'), ]
Fi<-FieldN[which(FieldN$PlotID=='F'), ]
BR<-FieldN[which(FieldN$PlotID=='BR'), ]
BM1<-FieldN[which(FieldN$PlotID=='BM1'), ]
BM2<-FieldN[which(FieldN$PlotID=='BM2'), ]

#Plot histograms
par(mfrow = c(3, 3))
hist(A$Tw,breaks=5,col="blue")
hist(B$Tw,breaks=5,col="blue")
hist(C$Tw,breaks=5,col="blue")
hist(D$Tw,breaks=5,col="blue")
hist(E$Tw,breaks=5,col="blue")
hist(Fi$Tw,breaks=5,col="blue")
hist(BR$Tw,breaks=5,col="blue")
hist(BM1$Tw,breaks=5,col="blue")
hist(BM2$Tw,breaks=5,col="blue")
hist(FieldN$Tw,breaks=20)

par(mfrow = c(3, 3))
hist(A$BroodTot,breaks=5,col="blue")
hist(B$BroodTot,breaks=5,col="blue")
hist(C$BroodTot,breaks=5,col="blue")
hist(D$BroodTot,breaks=5,col="blue")
hist(E$BroodTot,breaks=5,col="blue")
hist(Fi$BroodTot,breaks=5,col="blue")
hist(BR$BroodTot,breaks=5,col="blue")
hist(BM1$BroodTot,breaks=5,col="blue")
hist(BM2$BroodTot,breaks=5,col="blue")
hist(FieldN$BroodTot,breaks=20)


#Levene's Test for Homogeneity of variance
#Brood Content
fligner.test(BroodTot~PlotID,data=FieldN)
leveneTest(BroodTot~PlotID,data=FieldN)

#Worker content
fligner.test(Tw~PlotID,data=FieldN)
leveneTest(Tw~PlotID,data=FieldN)

n<-data.frame(table(FieldN$PlotID))
fligner.test(Freq~Var1,data=n)
leveneTest(Freq~Var1,data=n)

