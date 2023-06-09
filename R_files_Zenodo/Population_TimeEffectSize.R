library('compute.es')
library('pwr')

library(ggplot2)
library(gridExtra)
library(car)
library(lme4)
library(stargazer)
library(dotwhisker)
library(arsenal)
library(RColorBrewer)

SAVSmorph<-read.csv("~/Dropbox/NSF_collections_PostDoc/NSF_PRFB_Morphological_Analyses/FinalData_Code_BenhamBowie_2020_EvolApplications/MorphologicalData_BenhamBowie_EvolApplications_2020.csv",header=TRUE)

SAVSmorph<-SAVSmorph[!(SAVSmorph$Subspecies== "alaudinus" & SAVSmorph$wing > 70.0),]
SAVS.temp<-subset(SAVSmorph, age == "adult" & Subspecies == "alaudinus" & sex=="male" & Habitat=="marsh")
savs.sa<-SAVS.temp[,"SurfaceArea"]
SAVS.tidal.male<-SAVS.temp[complete.cases(savs.sa),]

tarsusSA.lm<-lm(SAVS.tidal.male$SurfaceArea~SAVS.tidal.male$wing)
BlengthTarsus.lm<-lm(SAVS.tidal.male$BillLength~SAVS.tidal.male$wing)
BwidthTarsus.lm<-lm(SAVS.tidal.male$BillWidth~SAVS.tidal.male$wing)
BdepthTarsus.lm<-lm(SAVS.tidal.male$BillDepth~SAVS.tidal.male$wing)

SAtarCorr<-residuals(tarsusSA.lm)
BlengthCorr<-residuals(BlengthTarsus.lm)
BwidthCorr<-residuals(BwidthTarsus.lm)
BdepthCorr<-residuals(BdepthTarsus.lm)

rows=rownames(SAVS.tidal.male)
print(length(rows))

SAVS.tidal.male<-cbind(SAVS.tidal.male, BlengthCorr, BwidthCorr, BdepthCorr, SAtarCorr)

SumTable<-tableby(region~Time, data = SAVS.tidal.male)
print(summary(SumTable))

#SAVS.tidal.m<-na.omit(SAVS.tidal.male)
#Clim.pca<-prcomp(SAVS.tidal.m[,c(15:20,25:27)], center=TRUE, scale=TRUE)
#print(summary(Clim.pca))
#print(Clim.pca)
#PC1<-Clim.pca$x[,1]
#PC2<-Clim.pca$x[,2]
#PC3<-Clim.pca$x[,3]

#Full.Dat<-cbind(SAVS.tidal.m, PC1,PC2,PC3)
#Full.Dat2<-as.data.frame(Full.Dat)


#PC1_PC2<-ggplot(Full.Dat2, aes(y=PC2, x=PC1, color=latitude, shape=factor(Time))) + #theme_bw() + geom_point(size=5) + labs(y="PC2 (20.3%)", x="PC1 (66.6%)") + #theme(axis.text=element_text(size=16),axis.title=element_text(size=18), #legend.title=element_text(size=16, color="black"), legend.text=element_text(color="black", #size=14)) 

#print(PC1_PC2)

SAVS.Humboldt<-subset(SAVS.tidal.male, Habitat=="marsh" & region=="HumboldtBay")
SAVS.Sfbay<-subset(SAVS.tidal.male, Habitat=="marsh" & region=="SFBay")
SAVS.Morro<-subset(SAVS.tidal.male, Habitat=="marsh" & region=="Morro")
SAVS.SanPablo<-subset(SAVS.tidal.male, Habitat=="marsh" & region=="SanPablo")
SAVS.Suisun<-subset(SAVS.tidal.male, Habitat=="marsh" & region=="Suisun")

#t.tests by Time period for each region and surface area
Humboldt.SA.t<-t.test(SAtarCorr ~Time, data=SAVS.Humboldt)
Sfbay.SA.t<-t.test(SAtarCorr ~Time, data=SAVS.Sfbay)
Morro.SA.t<-t.test(SAtarCorr ~Time, data=SAVS.Morro)
SAVS.SanPablo.t<-t.test(SAtarCorr ~Time, data=SAVS.SanPablo)
SAVS.Suisun.t<-t.test(SAtarCorr ~Time, data=SAVS.Suisun)


print(Humboldt.SA.t)
print(SAVS.SanPablo.t)
print(Sfbay.SA.t)
print(Morro.SA.t)
print(SAVS.Suisun.t)

print(mean(SAVS.Humboldt$Mean_Salinity))
print(mean(SAVS.Sfbay$Mean_Salinity))
print(mean(SAVS.Morro$Mean_Salinity))
print(mean(SAVS.SanPablo$Mean_Salinity))
print(mean(SAVS.Suisun$Mean_Salinity))

tvalue<-c(Morro.SA.t$statistic, Sfbay.SA.t$statistic, SAVS.Suisun.t$statistic, SAVS.SanPablo.t$statistic, Humboldt.SA.t$statistic)
n.t<-c(10,31,4,17,15)
n.c<-c(25,10,13,8,5)

ttest.result<-data.frame(id=1:5, t=tvalue, n.t=n.t, n.c=n.c)
print(ttest.result)

d<-tes(t,n.1=n.t,n.2=n.c,level=95, dig=2, id=id, data=ttest.result)
print(d)

CohenD<-d[,5]

VarD<-d[,6]
sd.D<-sqrt(VarD)
CI.D<-sd.D*1.96

ci.plot<-matrix(CI.D,1,5,byrow=TRUE)
print(ci.plot)

location=c("Morro Bay", "SF Bay", "Suisun Bay", "San Pablo Bay", "Humboldt Bay")

group.bar<-data.frame(location, CohenD)
print(group.bar)
error.bar <- function(x, y, upper, lower=upper, length=0.05,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) !=length(upper))
stop("vectors must be same length")
arrows(y+upper, x, y-lower, x, angle=90, code=3, lwd=0.5, length=length, ...)
 }

barz<-barplot(group.bar$CohenD, names.arg=location, xlim=c(-2,6), horiz=TRUE, col="blue", beside=TRUE, xlab="Effect size (d) of bill surface area", lwd=2)

error.bar(barz,group.bar$CohenD,ci.plot)
#title(main="Effect size of Time on bill width")
abline(v=0, col="black", lwd=2)
