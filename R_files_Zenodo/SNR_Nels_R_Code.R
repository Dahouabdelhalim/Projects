# Annotated R Code for Reid et al. 2020. Annual understory plant recovery
# dynamics in a temperate woodland mosaic during a decade of ecological
# restoration. Natural Areas Journal.

#Last updated by JLR on 2019-12-09

###############################################################################
# SET UP 
###############################################################################

#Libraries
library(BiodiversityR) #Diversity tools
library(plyr) #Data management tools
library(reshape2) #Data management tools
library(shape) #Draws nice arrows

#Data
veg<-read.csv('Data_Veg.csv',header=TRUE)
pheno<-read.csv('Data_Phenology.csv',header=TRUE)
traits<-read.csv('Data_SpecList.csv',header=TRUE)

#Format data
veg<-subset(veg,L=='H')
veg<-merge(veg,traits,by='Species')
veg<-merge(veg,pheno,by='Species')

#Unique years
Years<-unique(veg$Y)[order(unique(veg$Y))]

#define spring & summer phenology
pheno$PHENO<-ifelse(pheno$Start<=6 & pheno$End<=7,"Spring",'Summer')

###############################################################################
# FUNCTIONS
###############################################################################

## Restoration intervention timeline bars
RxLines<-function(){
  par(xpd=FALSE)
  abline(v=1+(328/365),lty=1,col='black',lwd=1) #prescribed burns
  abline(v=4+(57/365),lty=1,col='black',lwd=1)
  abline(v=7+(74/365),lty=1,col='black',lwd=1)
  abline(v=9+(319/365),lty=1,col='black',lwd=1)
  abline(v=12+(10/365),lty=1,col='black',lwd=1)
  abline(v=6,lty=2,col='black',lwd=1) #Juniper removal
  box()
  axis(side=1,at=seq(2001,2013,1),labels=NA)
}

## Rx toptext
RxToptext<-function(){
  mtext(side=3,line=0.5,text='Fire',at=1+(328/365),las=2)
  mtext(side=3,line=0.5,text='Fire',at=4+(57/365),las=2)
  mtext(side=3,line=0.5,text='Fire',at=7+(74/365),las=2)
  mtext(side=3,line=0.5,text='Fire',at=9+(319/365),las=2)
  mtext(side=3,line=0.5,text='Fire',at=12+(10/365),las=2)
  mtext(side=3,line=0.5,text='Thinning',at=6,las=2)
}


###############################################################################
###############################################################################
# SUMMARY STATISTICS
###############################################################################

#Total number of stems sampled
print(stemCount<-sum(veg$Stems))

#Percent of stems identified to species level
print(stemCompleteIDCount<-sum(veg[veg$idSpecLevel==1,]$Stems))
print(percentIDcomplete<-(stemCompleteIDCount/stemCount)*100)

#Total no. species encountered
print(noSpecies<-length(unique(veg[veg$idSpecLevel==1,]$Species)))

#Total no. genera encountered
print(noGenus<-length(unique(veg[veg$idSpecLevel==1,]$Genus)))

#Total no. families encountered  
print(noFamily<-length(unique(veg[veg$idSpecLevel==1,]$Family)))

#Top 3 families for stems
FamStems<-setNames(aggregate(veg$Stems,by=list(veg$Family),sum),c('Family','Stems'))
FamStems<-FamStems[order(-FamStems$Stems),]
FamStems$Percent<-(FamStems$Stems/sum(FamStems$Stems))*100

#Top 3 families for diversity
FamDiv<-with(veg[veg$idSpecLevel==1,],setNames(aggregate(Species,by=list(Family),function(x)length(unique(x))),c('Family','Species')))
FamDiv<-FamDiv[order(-FamDiv$Species),]

#Family summary table
FamTable<-merge(FamStems,FamDiv,by='Family',all.x=TRUE,all.y=TRUE)
FamTable$Percent<-sprintf("%.3f",FamTable$Percent)
write.csv(FamTable,"Table_FamTable.csv",row.names=FALSE)

#Stem count by species, 2001-2012
aggAll<-with(veg,setNames(aggregate(Stems,by=list(Species,Y),sum),c('Species','Y','Stems')))
aggAll<-dcast(aggAll,Species~Y,value.var='Stems',sum)
firstthree<-(aggAll$`2001`+aggAll$`2002`+aggAll$`2003`)/3
lastthree<-(aggAll$`2010`+aggAll$`2011`+aggAll$`2012`)/3
aggAll$threeYrDelta<-sprintf("%.0f",((lastthree-firstthree)/firstthree)*100)
aggAll$threeYrDelta[aggAll$threeYrDelta=='NaN']<-0
write.csv(aggAll,'Table_GroundStemsSummary.csv',row.names=FALSE)

#Most abundance species in 2001
DomSp<-function(Hab,Yr){
  v<-subset(veg,idSpecLevel==1 & H==Hab & Y==Yr)
  v1<-setNames(aggregate(v$Q,by=list(v$Species),length),c('Species','Q'))
  v1<-v1[with(v1,order(-Q)),]
  print(head(v1))
  plot(v1$Q~seq(1:length(rownames(v1))))
}

DomSp('J',2001)
DomSp('W',2001)
DomSp('R',2001)

###############################################################################
# CALCULATE FLORISTIC QUALITY INDICES
###############################################################################

#Remove incomplete species identities
vegC<-veg[veg$idSpecLevel==1,]

#Subset to a single habitat & season (or not)
#Vs<-vegC[vegC$H=='J',]
Vs<-vegC

#Create a dataframe to store Species Richness results
FQA<-data.frame(
  Y=as.integer(), #Year
  S.obs=as.integer(), #No. species observed
  S.chao1=as.numeric(), #Chao1 estimate
  se.chao1=as.numeric(), #Chao1 standard error
  S.ACE=as.numeric(), #ACE estimate
  se.ACE=as.numeric(), #ACE standard error
  SW=as.numeric(), #Shannon-Weiner diversity
  FA=as.numeric(), #Fisher's alpha
  PJ=as.numeric(), #Pielou's J
  meanC=as.numeric(), #Mean C score (omitting non-natives)
  nn=as.numeric(), #Native species richness
  percentNative=as.numeric(), #Percent native species
  I=as.numeric(), #FQI (omitting non-natives)
  AI=as.numeric(), #Adjusted FQI
  Iy<-as.numeric(), #Cover-weighted FQI
  Checker=as.numeric() #Should equal zero, if not, species richness calculated wrong somewhere
)

#Pretty FQA
pFQA<-data.frame(
  Year=as.integer(), #Year
  ObservedSpecies=as.integer(), #No. species observed
  NativeSpecies=as.numeric(), #Native species richness
  percentNative=as.numeric(), #Percent native species
  Chao=as.numeric(), #Chao1 estimate
  ACE=as.numeric(), #ACE estimate
  Shannon=as.numeric(), #Shannon-Weiner diversity
  Fisher=as.numeric(), #Fisher's alpha
  Pielou=as.numeric(), #Pielou's J
  I=as.numeric(), #FQI (omitting non-natives)
  meanC=as.numeric(), #Mean C score (omitting non-natives)
  AI=as.numeric(), #Adjusted FQI
  Iy<-as.numeric() #Cover-weighted FQI
)


for(i in 1:11){ #Loop to estimate species richness for each year
  Y=Years[i] #Select year_i
  
  #Create a vegan matrix for year i
  matrixi<-with(Vs[Vs$Y==Y,],setNames(aggregate(Stems,by=list(Species),sum),c('Species','Abun')))
  
  #Calculate estimated & observed species richness
  estR<-estimateR(matrixi[,2]) #Estimate species richness
  S.obs<-estR[1] #Extract observed species richness
  S.chao1<-estR[2] #Extract Chao1 estimate
  se.chao1<-estR[3] #Extract Chao1 standard error
  S.ACE<-estR[4] #Extract ACE estimate
  se.ACE<-estR[5] #Extract ACE standard error
  
  #Some diversity functions
  SW<-diversity(matrixi[,2]) #Calculate Shannon-Weiner diversity
  FA<-fisher.alpha(matrixi[,2]) #Calculate Fisher's Alpha
  PJ<-SW/log(S.obs) #Calculate Pielou's J
  
  #Some FQA functions
  syi<-merge(matrixi,traits,by='Species',all.x=TRUE,all.y=FALSE)
  syi$C[syi$C==99]<-0 #Change non-native C to zero
  syi<-syi[!is.na(syi$C),] #remove NAs
  meanC<-mean(syi$C) #Calculate meanC
  nn<-length(syi$Species[syi$NT=='Nt']) #Calculate native species richness (nn)
  nt<-length(syi$Species) #Calculate total species richness (nt)
  percentNative<-(nn/nt)*100 #Calculate percent native species
  I<-meanC*sqrt(nn) #Calculate FQI
  AI<-100*(meanC/10)*(sqrt(nn)/sqrt(nt)) #Calculate adjusted FQI
  
  #Cover-weighted FQI
  library(plyr)
  Vs$Q10<-round_any(Vs$Q,10,ceiling)   #Create a variable combining 5 & 10 m quadrat observations
  castSprI<-acast(Vs[Vs$Y==Y,],Species~T+Q10,value.var='Cov',max) #Take the max value for each species to prevent double-counting between SP & SU
  castSprI[is.infinite(castSprI)]<-0 #Fill in zeros
  meanCov<-setNames(as.data.frame(rowMeans(castSprI)),'meanCov')
  meanCov$Species<-rownames(meanCov)
  meanCov<-merge(meanCov,traits,by='Species',all.x=TRUE,all.y=FALSE)
  meanCov$C[meanCov$C==99]<-0 #Make non-native C values 0
  meanCov$MCXC<-meanCov$meanCov*meanCov$C
  Cy<-sum(meanCov$MCXC)
  Iy<-Cy/length(meanCov$Species)  
  
  #Add a checker - species richness calculated two different ways should be the same
  Checker<-nt-S.obs
  
  #FQA Pretty
  chaoPretty<-paste(sprintf('%.1f',S.chao1),sprintf('%.1f',se.chao1),sep= " ± ")
  acePretty<-paste(sprintf('%.1f',S.ACE),sprintf('%.1f',se.ACE),sep= " ± ")
  fqaPretty<-as.data.frame(cbind(Year=Y,
                                 S.obs,nn,percentNative=sprintf('%.1f',percentNative),
                                 chaoPretty,acePretty,
                                 Shannon=sprintf('%.2f',SW),Fisher=sprintf('%.2f',FA),
                                 Pielou=sprintf('%.2f',PJ),
                                 I=sprintf('%.2f',I),meanC=sprintf('%.2f',meanC),
                                 AI=sprintf('%.2f',AI),Iy=sprintf('%.2f',Iy)))  
  
  #FQA Stat
  fqaStat<-as.data.frame(cbind(Y,
                               S.obs,nn,percentNative,
                               S.chao1,S.ACE,
                               SW,FA,PJ,
                               I,meanC,AI,Iy))
  
  #Add results to Species Richness output
  FQA<-rbind(FQA,fqaStat)
  pFQA<-rbind(pFQA,fqaPretty)
}


#Format output
FQA<-as.data.frame(t(FQA[,2:(length(colnames(FQA)))]))
colnames(FQA)<-Years

#Simple linear regression for each variable
Time<-c(0,1,2,3,4,6,7,8,9,10,11)

#Create a dataframe to store Species Richness results
LMResults<-data.frame(
  adjR<-as.character(),
  Beta=as.character(),
  t=as.character()
)

Pnone<-data.frame(pnone=as.numeric())

for(i in 1:length(rownames(FQA))){
  M<-summary(lm(as.numeric(FQA[i,1:11])~Time))
  adjR<-sprintf("%.2f",M$adj.r.squared)
  Beta<-paste(sprintf("%.3f",M$coefficients[2,1]),sprintf("%.3f",M$coefficients[2,2]),sep=" ± ")
  t<-sprintf("%.2f",M$coefficients[2,3])
  p<-as.numeric(M$coefficients[2,4])
  results<-cbind(adjR,Beta,t)
  LMResults<-rbind(LMResults,results)
  Pnone=rbind(Pnone,p)
}

PBon<-p.adjust(Pnone[,1],method='bonferroni')

LMResults$Padj<-paste(sprintf('%.4f',PBon),
                      ifelse(PBon<=0.0001,'****',
                             ifelse(PBon<=0.001,'***',
                                    ifelse(PBon<=0.01,'**',
                                           ifelse(PBon<=0.05,'*',
                                                  ifelse(PBon<=0.1,'.',''))))),sep="")

PrettyFQA<-t(pFQA[,2:length(colnames(pFQA))])
colnames(PrettyFQA)<-Years
PrettyFQA<-cbind(PrettyFQA,LMResults)

#Save results
#write.csv(PrettyFQA,'Table_FQA_Forest.csv',row.names=TRUE)
#write.csv(t(FQA),'Table_FQA_DBW_Raw.csv',row.names=TRUE)

###############################################################################
# ANALYZE FLORISTIC QUALITY DYNAMICS
###############################################################################

#Read FQA data from exported tables
FS<-read.csv("Table_FQA_DBW_Raw.csv")
  FS$Y<-as.numeric(FS$X)+0.5
FG<-read.csv("Table_FQA_Glade_Raw.csv")
  FG$Y<-as.numeric(FG$X)+0.5
FW<-read.csv("Table_FQA_Woodland_Raw.csv")
  FW$Y<-as.numeric(FW$X)+0.5
FF<-read.csv("Table_FQA_Forest_Raw.csv")
  FF$Y<-as.numeric(FF$X)+0.5

#Simple years for analysis
Ys<-c(1,2,3,4,5,7,8,9,10,11,12)
  
#Did species density increase linearly or logistically?
plot(FS[,3]~Ys)
LinMod<-lm(FS[,3]~Ys);summary(LinMod)
plot(resid(LinMod))

#Linear model for glade species density
plot(FG[,3]~Ys)
LinMod<-lm(FG[,3]~Ys);summary(LinMod)
abline(lm(FG[,3]~Ys))
plot(resid(LinMod))

#Linear model for woodland species density
plot(FW[,3]~Ys)
LinMod<-lm(FW[,3]~Ys);summary(LinMod)
abline(lm(FW[,3]~Ys))
plot(resid(LinMod))

#Log model for woodland species density
plot(FW[,3]~Ys)
LogMod<-lm(FW[,3]~log(Ys));summary(LogMod)
lines(Ys,predict(LogMod),col='red')

#Linear model for forest species density
plot(FF[,3]~Ys)
LinMod<-lm(FF[,3]~Ys);summary(LinMod)
abline(lm(FW[,3]~Ys))
plot(resid(LinMod))

#Linear model for Site-level I
plot(FS$I~Ys)
LinMod<-lm(FS$I~Ys);summary(LinMod)
LogMod<-lm(FS$I~log(Ys));summary(LogMod)
lines(Ys,predict(LinMod),col='red')

#Linear model for glade I
plot(FW$I~Ys)
LinMod<-lm(FW$I~Ys);summary(LinMod)
LogMod<-lm(FW$I~log(Ys));summary(LogMod)

#model for DBW Iy
plot(FS$Iy~Ys)
LinMod<-lm(FS$Iy~Ys);summary(LinMod)
LogMod<-lm(FS$Iy~log(Ys));summary(LogMod)
lines(Ys,predict(LogMod),col='red')

#Model for glade Iy
plot(FG$Iy~Ys)
LinMod<-lm(FG$Iy~Ys);summary(LinMod)
LogMod<-lm(FG$Iy~log(Ys));summary(LogMod)
lines(Ys,predict(LinMod),col='red')

#Model for woodland Iy
plot(FW$Iy~Ys)
LinMod<-lm(FW$Iy~Ys);summary(LinMod)
LogMod<-lm(FW$Iy~log(Ys));summary(LogMod)
lines(Ys,predict(LinMod),col='red')

#Model for forest Iy
plot(FF$Iy~Ys)
LogMod<-lm(FF$Iy~log(Ys));summary(LogMod)
lines(Ys,predict(LogMod),col='red')

#Broken stick model for forest Iy
library(segmented)
FIY<-FF$Iy
LinMod2<-lm(FIY~Ys)
SegMod<-segmented(LinMod2,seg.Z=~Ys);summary(SegMod)
lines(Ys,predict(SegMod),col='red')

#Graphics function for lines and shaded area
fqaGraphics<-function(x){
  abline(v=2001+(328/365),lty=2,col='gray85') #prescribed burns
  abline(v=2004+(57/365),lty=2,col='gray85')
  abline(v=2007+(74/365),lty=2,col='gray85')
  abline(v=2009+(319/365),lty=2,col='gray85')
  abline(v=2012+(10/365),lty=2,col='gray85')
  polygon(c(2005.8,2005.8,2007,2007),c(-1000,3000,3000,-1000),col='gray85',border=NA) #Shade gray - shrubs
  box()
  par(xpd=NA)
  text(2002,x,pos=3,offset=1.3,'fire',srt=45,col='gray75',adj=0)
  text(2004,x,pos=3,offset=1.3,'fire',srt=45,col='gray75',adj=0)
  text(2007.7,x,pos=3,offset=1.3,'fire',srt=45,col='gray75',adj=0)
  text(2010,x,pos=3,offset=1.3,'fire',srt=45,col='gray75',adj=0)
  text(2012,x,pos=3,offset=1.3,'fire',srt=45,col='gray75',adj=0)
  text(2006.5,x,pos=3,offset=1.7,'juniper\\nremoval',srt=45,col='gray75',adj=0)
  par(new=TRUE)
  par(xpd=FALSE)
}

#Axis function to prevent overplotting
axisFQA<-function(x){
  axis(side=1,at=seq(2002,2010,4))
  axis(side=1,at=seq(2004,2012,4))
}

DBWpalCOL<-c('darkgreen','black','maroon')
DBWpalBG<-c('yellow','white','pink')

#FQA Scatterplot
tiff('Fig_FQA_Scatter_Site.tif',res=300,height=80,width=180,units='mm')
par(mfrow=c(1,3),mai=c(1,0.5,0.5,0.1),omi=c(0,0,0,0),mgp=c(2.5,0.7,0))

plot(FS$S.obs~FS$Y,type='n',ylim=c(150,210),las=1,ylab='Obs. native species density',xlab=NA,xaxt='n')
fqaGraphics(210)
points(FS$nn~FS$Y,pch=16)
abline(lm(FS$nn~FS$Y))
axisFQA()

par(mgp=c(2,0.7,0))
plot(FS$I~FS$Y,type='n',las=1,ylim=c(49,56),ylab=expression(paste('Floristic quality (',italic(I),')')),xlab=NA,xaxt='n')
fqaGraphics(56)
points(FS$I~FS$Y,pch=16)
abline(lm(FS$I~FS$Y))
axisFQA()

plot(FS$Iy~FS$Y,type='n',las=1,
     ylab=expression(paste('Cover-weighted floristic quality ( ',italic(I[gamma]),')')),
     xlab=NA,ylim=c(1,2.6),xaxt='n')
fqaGraphics(2.6)
points(FS$Iy~FS$Y,pch=16)
lines(Ys+2000.5,predict(lm(FS$Iy~log(Ys))))
axisFQA()

dev.off()

#FQA Scatterplot - Community level
tiff('Fig_FQA_Scatter_Comm.tif',res=300,height=80,width=180,units='mm')
par(mfrow=c(1,3),mai=c(1,0.5,0.5,0.1),omi=c(0,0,0,0),mgp=c(2.5,0.7,0))

plot(FS$S.obs~FS$Y,type='n',ylim=c(90,150),las=1,ylab='Obs. native species density',xlab=NA,xaxt='n')
fqaGraphics(150)
lines(FF$nn~FF$Y,pch=22,col=DBWpalCOL[3],bg=DBWpalBG[3],type='b',lty=2)
points(FG$nn~FG$Y,pch=21,bg=DBWpalBG[1],col=DBWpalCOL[1])
abline(lm(FG$nn~FG$Y),col=DBWpalCOL[1])
points(FW$nn~FW$Y,pch=24,bg=DBWpalBG[2],col=DBWpalCOL[2])
lines(Ys+2000.5,predict(lm(FW$nn~log(Ys))),lty=3,col=DBWpalCOL[2])
axisFQA()
legend('topleft',legend=c('Juniper','Woodland','Forest'),
       pch=c(21,24,22),pt.bg=DBWpalBG,col=DBWpalCOL,lty=c(1,3,2),bty='n')

plot(FS$I~FS$Y,type='n',las=1,ylim=c(35,45),ylab=expression(paste('Floristic quality (',italic(I),')')),xlab=NA,xaxt='n')
fqaGraphics(45)
lines(FW$I~FW$Y,pch=24,bg=DBWpalBG[2],col=DBWpalCOL[2],type='b',lty=3)
lines(FF$I~FF$Y,pch=22,bg=DBWpalBG[3],col=DBWpalCOL[3],type='b',lty=2)
lines(FG$I~FG$Y,pch=21,bg=DBWpalBG[1],col=DBWpalCOL[1],type='b')
axisFQA()

plot(FS$Iy~FS$Y,type='n',las=1,
     ylab=expression(paste('Cover-weighted floristic quality ( ',italic(I[gamma]),')')),xlab=NA,ylim=c(1,5),xaxt='n')
fqaGraphics(5)
points(FG$Iy~FG$Y,pch=21,bg=DBWpalBG[1],col=DBWpalCOL[1])
abline(lm(FG$Iy~FG$Y),col=DBWpalCOL[1])
points(FF$Iy~FF$Y,pch=22,bg=DBWpalBG[3],col=DBWpalCOL[3])
lines(FF$Y,predict(segmented(lm(FIY~Ys))),col=DBWpalCOL[3],lty=2)
points(FW$Iy~FW$Y,pch=24,bg=DBWpalBG[2],col=DBWpalCOL[2])
abline(lm(FW$Iy~FW$Y),lty=3)
axisFQA()

dev.off()
###############################################################################
#POLYGON PLOT - GUILD S + Cover ~ TIME
###############################################################################

#Unique habitat list
Habitats<-c('J','W','R')

#Subset to drop incomplete IDs
v<-droplevels(subset(veg,idSpecLevel==1))

#Calculate species richness for each guild per habitat per year
GuildRich<-function(x){
vx<-droplevels(subset(v,H==Habitats[x]))
vx2<-setNames(aggregate(vx$Species,by=list(vx$Y,vx$PHYS),function(x)length(unique(x))),c('Y','Guild','Species'))
vx3<-dcast(vx2,Y~Guild)
vx3$S<-rowSums(vx3[,2:length(colnames(vx3))])
assign('vx3',vx3,envir=.GlobalEnv)
}

#Draw polygons for species richness per guild per year
GuildPoly<-function(vx3){
polygon(c(vx3$Y,rev(vx3$Y)),c(rep(0,11),rev(vx3$S)),col='blue') #Forbs
polygon(c(vx3$Y,rev(vx3$Y)),c(vx3$FORB,rev(vx3$S)),col='pink') #TREE
polygon(c(vx3$Y,rev(vx3$Y)),c(vx3$FORB+vx3$TREE,rev(vx3$S)),col='green') #Grass
polygon(c(vx3$Y,rev(vx3$Y)),c(vx3$FORB+vx3$TREE+vx3$GRASS,rev(vx3$S)),col='seagreen') #Shrub
polygon(c(vx3$Y,rev(vx3$Y)),c(vx3$FORB+vx3$TREE+vx3$GRASS+vx3$SHRUB,rev(vx3$S)),col='yellow') #Sedge
polygon(c(vx3$Y,rev(vx3$Y)),c(vx3$FORB+vx3$TREE+vx3$GRASS+vx3$SHRUB+vx3$SEDGE,rev(vx3$S)),col='gray80') #Woody vine
polygon(c(vx3$Y,rev(vx3$Y)),c(vx3$FORB+vx3$TREE+vx3$GRASS+vx3$SHRUB+vx3$SEDGE+vx3$WVINE,rev(vx3$S)),col='red') #Fern
}

#Calculate mean percent cover per 0.5 m2 for each guild per habitat per year (taking the max of 5 & 10 quads)
GuildCov<-function(x){
  vx<-droplevels(subset(v,H==Habitats[x]))
  vx$Q10<-paste(vx$T,round_any(vx$Q,10,ceiling))
  vx2<-setNames(aggregate(vx$Cov,by=list(vx$Y,vx$PHYS,vx$Q10),max),c('Y','PHYS','Q10','Cov'))
  vx3<-dcast(vx2,Y+PHYS~Q10)
  vx3[is.na(vx3)]<-0
  vx3$Cov<-rowMeans(vx3[,3:length(colnames(vx3))])
  vx4<-dcast(vx3,Y~PHYS,value.var='Cov')
  vx4$Cov<-rowSums(vx4[,2:length(colnames(vx4))])
  assign('vx4',vx4,envir=.GlobalEnv)
}

#Draw polygons for mean percent cover per guild per year
GuildCovPoly<-function(vx4){
  polygon(c(vx4$Y,rev(vx4$Y)),c(rep(0,11),rev(vx4$Cov)),col='blue') #Forbs
  polygon(c(vx4$Y,rev(vx4$Y)),c(vx4$FORB,rev(vx4$Cov)),col='pink') #TREE
  polygon(c(vx4$Y,rev(vx4$Y)),c(vx4$FORB+vx4$TREE,rev(vx4$Cov)),col='green') #Grass
  polygon(c(vx4$Y,rev(vx4$Y)),c(vx4$FORB+vx4$TREE+vx4$GRASS,rev(vx4$Cov)),col='seagreen') #Shrub
  polygon(c(vx4$Y,rev(vx4$Y)),c(vx4$FORB+vx4$TREE+vx4$GRASS+vx4$SHRUB,rev(vx4$Cov)),col='yellow') #Sedge
  polygon(c(vx4$Y,rev(vx4$Y)),c(vx4$FORB+vx4$TREE+vx4$GRASS+vx4$SHRUB+vx4$SEDGE,rev(vx4$Cov)),col='gray80') #Woody vine
  polygon(c(vx4$Y,rev(vx4$Y)),c(vx4$FORB+vx4$TREE+vx4$GRASS+vx4$SHRUB+vx4$SEDGE+vx4$WVINE,rev(vx4$Cov)),col='red') #Fern
}

#Analyze guild richness by treatment
GuildRich(1)
vx3$ForbPC<-with(vx3,FORB/S)
with(vx3,plot(ForbPC~Ys))
with(vx3,summary(lm(ForbPC~Ys)))
vx3$BA<-c(rep('B',5),rep('A',6))
with(vx3,summary(lm(ForbPC~BA)))

#Assign C categories
v$catC<-ifelse(v$C<=3,'Ruderal',
               ifelse(v$C<=6,'Matrix',
                      ifelse(v$C<=10,'Conservative','Non-native')))

#Calculate species richness for each C per habitat per year
cRich<-function(x){
  vx<-droplevels(subset(v,H==Habitats[x]))
  vx2<-setNames(aggregate(vx$Species,by=list(vx$Y,vx$catC),function(x)length(unique(x))),c('Y','catC','Species'))
  vx3<-dcast(vx2,Y~catC)
  vx3$S<-rowSums(vx3[,2:length(colnames(vx3))])
  assign('vx3',vx3,envir=.GlobalEnv)
}


cols<-c('black','gray40','gray80','white')

#Draw polygons for species richness per C per year
cPoly<-function(vx3){
  polygon(c(vx3$Y,rev(vx3$Y)),c(rep(0,11),rev(vx3$S)),col=cols[1]) #Non-native
  polygon(c(vx3$Y,rev(vx3$Y)),c(vx3$`Non-native`,rev(vx3$S)),col=cols[2]) #Ruderal
  polygon(c(vx3$Y,rev(vx3$Y)),c(vx3$`Non-native`+vx3$Ruderal,rev(vx3$S)),col=cols[3]) #Matrix
  polygon(c(vx3$Y,rev(vx3$Y)),c(vx3$`Non-native`+vx3$Ruderal+vx3$Matrix,rev(vx3$S)),col=cols[4]) #Conservative
}

#Calculate mean percent cover per 0.5 m2 for each C per habitat per year (taking the max of 5 & 10 quads)
cCov<-function(x){
  vx<-droplevels(subset(v,H==Habitats[x]))
  vx$Q10<-paste(vx$T,round_any(vx$Q,10,ceiling))
  vx2<-setNames(aggregate(vx$Cov,by=list(vx$Y,vx$catC,vx$Q10),max),c('Y','catC','Q10','Cov'))
  vx3<-dcast(vx2,Y+catC~Q10)
  vx3[is.na(vx3)]<-0
  vx3$Cov<-rowMeans(vx3[,3:length(colnames(vx3))])
  vx4<-dcast(vx3,Y~catC,value.var='Cov')
  vx4$Cov<-rowSums(vx4[,2:length(colnames(vx4))])
  assign('vx4',vx4,envir=.GlobalEnv)
}

#Draw polygons for mean percent cover per C per year
cCovPoly<-function(vx4){
  polygon(c(vx4$Y,rev(vx4$Y)),c(rep(0,11),rev(vx4$Cov)),col=cols[1]) #Non-native
  polygon(c(vx4$Y,rev(vx4$Y)),c(vx4$`Non-native`,rev(vx4$Cov)),col=cols[2]) #Ruderal
  polygon(c(vx4$Y,rev(vx4$Y)),c(vx4$`Non-native`+vx4$Ruderal,rev(vx4$Cov)),col=cols[3]) #Matrix
  polygon(c(vx4$Y,rev(vx4$Y)),c(vx4$`Non-native`+vx4$Ruderal+vx4$Matrix,rev(vx4$Cov)),col=cols[4]) #Conservative
}

#Plot change in mean cover by (A) physiognomy and (B) conservatism
tiff('Fig_PolygonsMeanCov.tif',res=300,height=140,width=180,units='mm')
par(mfrow=c(2,3),mai=c(0.1,0.1,0,0),omi=c(0.8,0.5,0.5,1.0),mgp=c(2.5,0.7,0))

GuildCov(1)
plot(vx4$Cov~vx4$Y,type='l',ylim=c(0,100),las=1,xlab='',xaxt='n')
mtext(side=2,outer=TRUE,'Mean cover (%)',line=2)
fqaGraphics(100)
GuildCovPoly(vx4)
legend('bottomleft','Juniper',text.col='white',bty='n')

GuildCov(2)
plot(vx4$Cov~vx4$Y,type='l',ylim=c(0,100),las=1,xlab='',xaxt='n',ylab='',yaxt='n')
fqaGraphics(100)
GuildCovPoly(vx4)
legend('bottomleft','Woodland',text.col='white',bty='n')

GuildCov(3)
plot(vx4$Cov~vx4$Y,type='l',ylim=c(0,100),las=1,xlab='',xaxt='n',ylab='',yaxt='n')
fqaGraphics(100)
GuildCovPoly(vx4)
legend('bottomleft','Forest',text.col='white',bty='n')

par(xpd=NA)
legend(2012.3,100,legend=c('Forb','Tree','Grass','Shrub','Sedge','Woody vine','Fern'),
       fill=c('blue','pink','green','seagreen','yellow','gray80','red'),bty='n',title='Physiognomy')
par(xpd=FALSE)

cCov(1)
plot(vx4$Cov~vx4$Y,type='l',ylim=c(0,100),las=1,xlab='',ylab='Mean cover (%)',xaxt='n')
fqaGraphics(1000)
cCovPoly(vx4)
axisFQA()
legend('bottomleft','Juniper',text.col='white',bty='n')

cCov(2)
plot(vx4$Cov~vx4$Y,type='l',ylim=c(0,100),las=1,xlab='',xaxt='n',ylab='',yaxt='n')
fqaGraphics(1000)
cCovPoly(vx4)
axisFQA()
legend('bottomleft','Woodland',text.col='white',bty='n')

cCov(3)
plot(vx4$Cov~vx4$Y,type='l',ylim=c(0,100),las=1,xlab='',xaxt='n',ylab='',yaxt='n')
fqaGraphics(1000)
cCovPoly(vx4)
axisFQA()
legend('bottomleft','Forest',text.col='white',bty='n')

par(xpd=NA)
legend(2012.3,100,legend=c('Non-native','Ruderal','Matrix','Conservative'),fill=c(cols[1:4]),
       title="Conservatism",bty='n')
par(xpd=FALSE)

dev.off()

#Plot change in species density by (A) physiognomy and (B) conservatism
tiff('Fig_PolygonsSpecDens.tif',res=300,height=140,width=180,units='mm')
par(mfrow=c(2,3),mai=c(0.1,0.1,0,0),omi=c(0.8,0.5,0.5,1.0),mgp=c(2.5,0.7,0))

GuildRich(1)
plot(vx3$S~vx3$Y,type='l',ylim=c(0,160),las=1,xlab='',xaxt='n',ylab='')
fqaGraphics(160)
GuildPoly(vx3)
legend('bottomleft','Juniper',text.col='white',bty='n')

GuildRich(2)
plot(vx3$S~vx3$Y,type='l',ylim=c(0,160),las=1,xlab='',ylab='',xaxt='n',yaxt='n')
fqaGraphics(160)
GuildPoly(vx3)
legend('bottomleft','Woodland',text.col='white',bty='n')

GuildRich(3)
plot(vx3$S~vx3$Y,type='l',ylim=c(0,160),las=1,xlab='',ylab='',xaxt='n',yaxt='n')
fqaGraphics(160)
GuildPoly(vx3)
legend('bottomleft','Forest',text.col='white',bty='n')

par(xpd=NA)
legend(2012.3,160,legend=c('Forb','Tree','Grass','Shrub','Sedge','Woody vine','Fern'),
       fill=c('blue','pink','green','seagreen','yellow','gray80','red'),bty='n',title='Physiognomy')
par(xpd=FALSE)

cRich(1)
plot(vx3$S~vx3$Y,type='l',ylim=c(0,160),las=1,xlab='',ylab='',xaxt='n')
fqaGraphics(1600)
cPoly(vx3)
axisFQA()
legend('bottomleft','Juniper',text.col='white',bty='n')

cRich(2)
plot(vx3$S~vx3$Y,type='l',ylim=c(0,160),las=1,xlab='',ylab='',xaxt='n',yaxt='n')
fqaGraphics(1600)
cPoly(vx3)
axisFQA()
legend('bottomleft','Woodland',text.col='white',bty='n')

cRich(3)
plot(vx3$S~vx3$Y,type='l',ylim=c(0,160),las=1,xlab='',ylab='',xaxt='n',yaxt='n')
fqaGraphics(1600)
cPoly(vx3)
axisFQA()
legend('bottomleft','Forest',text.col='white',bty='n')

par(xpd=NA)
legend(2012.3,160,legend=c('Non-native','Ruderal','Matrix','Conservative'),fill=c(cols[1:4]),
       title="Conservatism",bty='n')
par(xpd=FALSE)

mtext(side=2,outer=TRUE,'Species density',line=2)

dev.off()

###############################################################################
#POLYGON PLOT - C ~ TIME
###############################################################################

#Subset to drop incomplete IDs
v<-droplevels(subset(veg,idSpecLevel==1))

#Assign C categories
v$catC<-ifelse(v$C<=3,'Ruderal',
               ifelse(v$C<=6,'Matrix',
                      ifelse(v$C<=10,'Conservative','Non-native')))

#Calculate species richness for each C per habitat per year
cRich<-function(x){
  vx<-droplevels(subset(v,H==Habitats[x]))
  vx2<-setNames(aggregate(vx$Species,by=list(vx$Y,vx$catC),function(x)length(unique(x))),c('Y','catC','Species'))
  vx3<-dcast(vx2,Y~catC)
  vx3$S<-rowSums(vx3[,2:length(colnames(vx3))])
  assign('vx3',vx3,envir=.GlobalEnv)
}


cols<-c('black','gray40','gray80','white')

#Draw polygons for species richness per C per year
cPoly<-function(vx3){
  polygon(c(vx3$Y,rev(vx3$Y)),c(rep(0,11),rev(vx3$S)),col=cols[1]) #Non-native
  polygon(c(vx3$Y,rev(vx3$Y)),c(vx3$`Non-native`,rev(vx3$S)),col=cols[2]) #Ruderal
  polygon(c(vx3$Y,rev(vx3$Y)),c(vx3$`Non-native`+vx3$Ruderal,rev(vx3$S)),col=cols[3]) #Matrix
  polygon(c(vx3$Y,rev(vx3$Y)),c(vx3$`Non-native`+vx3$Ruderal+vx3$Matrix,rev(vx3$S)),col=cols[4]) #Conservative
}

#Calculate mean percent cover per 0.5 m2 for each C per habitat per year (taking the max of 5 & 10 quads)
cCov<-function(x){
  vx<-droplevels(subset(v,H==Habitats[x]))
  vx$Q10<-paste(vx$T,round_any(vx$Q,10,ceiling))
  vx2<-setNames(aggregate(vx$Cov,by=list(vx$Y,vx$catC,vx$Q10),max),c('Y','catC','Q10','Cov'))
  vx3<-dcast(vx2,Y+catC~Q10)
  vx3[is.na(vx3)]<-0
  vx3$Cov<-rowMeans(vx3[,3:length(colnames(vx3))])
  vx4<-dcast(vx3,Y~catC,value.var='Cov')
  vx4$Cov<-rowSums(vx4[,2:length(colnames(vx4))])
  assign('vx4',vx4,envir=.GlobalEnv)
}

#Draw polygons for mean percent cover per C per year
cCovPoly<-function(vx4){
  polygon(c(vx4$Y,rev(vx4$Y)),c(rep(0,11),rev(vx4$Cov)),col=cols[1]) #Non-native
  polygon(c(vx4$Y,rev(vx4$Y)),c(vx4$`Non-native`,rev(vx4$Cov)),col=cols[2]) #Ruderal
  polygon(c(vx4$Y,rev(vx4$Y)),c(vx4$`Non-native`+vx4$Ruderal,rev(vx4$Cov)),col=cols[3]) #Matrix
  polygon(c(vx4$Y,rev(vx4$Y)),c(vx4$`Non-native`+vx4$Ruderal+vx4$Matrix,rev(vx4$Cov)),col=cols[4]) #Conservative
}

#Plot C frequency of S & Cover by habitat & year
tiff('Fig_CPoly.tif',res=300,height=160,width=180,units='mm')
par(mfrow=c(2,3),mai=c(0.3,0.5,0.5,0),omi=c(0.8,0,0,0.2),mgp=c(2.5,0.7,0))

cRich(1)
plot(vx3$S~vx3$Y,type='l',ylim=c(0,160),las=1,xlab='',ylab='Obs. species density',xaxt='n')
fqaGraphics(160)
cPoly(vx3)
axisFQA()
legend('bottomleft','Glade',text.col='white',bty='n')

cRich(3)
plot(vx3$S~vx3$Y,type='l',ylim=c(0,160),las=1,xlab='',ylab='',xaxt='n')
fqaGraphics(160)
cPoly(vx3)
axisFQA()
legend('bottomleft','Woodland',text.col='white',bty='n')

cRich(2)
plot(vx3$S~vx3$Y,type='l',ylim=c(0,160),las=1,xlab='',ylab='',xaxt='n')
fqaGraphics(160)
cPoly(vx3)
axisFQA()
legend('bottomleft','Forest',text.col='white',bty='n')

cCov(1)
plot(vx4$Cov~vx4$Y,type='l',ylim=c(0,100),las=1,xlab='',ylab='Mean cover (%)',xaxt='n')
fqaGraphics(100)
cCovPoly(vx4)
axisFQA()
legend('bottomleft','Glade',text.col='white',bty='n')

cCov(3)
plot(vx4$Cov~vx4$Y,type='l',ylim=c(0,100),las=1,xlab='',xaxt='n',ylab='')
fqaGraphics(100)
cCovPoly(vx4)
axisFQA()
legend('bottomleft','Woodland',text.col='white',bty='n')

par(xpd=NA)
legend(1990,-30,ncol=4,legend=c('Non-native','Ruderal species','Matrix species','Conservative species'),
       fill=c(cols[1:4]))
par(xpd=FALSE)

cCov(2)
plot(vx4$Cov~vx4$Y,type='l',ylim=c(0,100),las=1,xlab='',xaxt='n',ylab='')
fqaGraphics(100)
cCovPoly(vx4)
axisFQA()
legend('bottomleft','Forest',text.col='white',bty='n')

dev.off()

###############################################################################
# NMDS
###############################################################################

#Pairwise adonis Function - Paula Martins/Robert Danczak rendition of original function by Pedro Martinez Arbizu
pairwise.adonis <- function(x,factors, sim.method,binary=FALSE, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~ 
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method,binary=FALSE);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

#Subset veg to included only complete species identifications
V<-droplevels(subset(veg,idSpecLevel==1))

#Add a Year_Transect variable
V$YT<-with(V,paste(Y,T,sep='_'))

#Cast into a vegan-style matrix, calculating average percent cover per quadrat
VegMat<-acast(V,YT~ACRON,value.var='Cov',function(x)round_any(sum(x)/10,accuracy=1,f=ceiling))

#Remove species that account for <5% of records (i.e., those present in <16.5 transect years)
SpecCount<-apply(VegMat,2,function(x)length(x[x>0]))
VegMat<-VegMat[,SpecCount>17]

#Environmental dataframe
env<-setNames(as.data.frame(rownames(VegMat)),'YT')
env$H<-substr(env$YT,6,6)
env$HY<-paste(env$H,substr(env$YT,1,4),sep="_")

#Abbreviate species names for plotting
#colnames(DM)<-make.cepnames(colnames(DM)) #No need to do this since we already have abbreviations

#Relativize by species maximum percent cover
#VegMat<-apply(VegMat,2,function(x)x/max(x)) #Better fit w/o relativization

#Set species plotting prioritization
priSpp<-colSums(VegMat)

#Run NMDS
MD1<-metaMDS(VegMat,distance='bray',k=3,trymax=1000)

#Check stress
stressplot(MD1)

#PERMANOVA
adonis(VegMat~HY,data=env,method='bray',binary=FALSE,perm=999) #Regular permanova

#Pairwise PERMANOVA
ad1<-pairwise.adonis(x=VegMat,factors=env$HY,sim.method='bray',binary=FALSE,p.adjust.m='bonferroni');ad1 #Pairwise permanova

#PERMDISP (test for dispersion)
bd1<-betadisper(vegdist(VegMat,method='bray'),group=env$HY,type='centroid',bias.adjust=TRUE)
anova(bd1)
tk1<-as.data.frame(TukeyHSD(bd1)$group);tk1
tk1$pairsForward<-gsub("-"," vs ",rownames(tk1)) #Add a merging column
#Add a reversed merging column (for some reason, the pairs are backwards in PERMDISP vs. PERMANOVA)
tk1$pairsReversed<-paste(substr(tk1$pairsForward,11,16),substr(tk1$pairsForward,1,6),sep=' vs ')

#Merge PERMANOVA & PERMDISP results
permResults<-merge(ad1,tk1,by.x='pairs',by.y='pairsReversed',all.x=TRUE)

#Pairs I want
PairsIWant<-c('J_2001 vs R_2001',
              'J_2001 vs W_2001',
              'R_2001 vs W_2001',
              'J_2012 vs R_2012',
              'J_2012 vs W_2012',
              'R_2012 vs W_2012',
              'J_2001 vs J_2012',
              'R_2001 vs R_2012',
              'W_2001 vs W_2012')

#Subset to include only selected pairs
permResults<-permResults[permResults$pairs %in% PairsIWant,]

#Format PERMANOVA/PERMDISP results table
permResults<-subset(permResults,select=c('pairs','F.Model','R2','p.value','diff','lwr','upr','p adj'))
permResults$diff<-permResults$diff*-1
permResults$lwr<-permResults$lwr*-1
permResults$upr<-permResults$upr*-1
permResults$L<-pmin(permResults$lwr,permResults$upr)
permResults$U<-pmax(permResults$lwr,permResults$upr)
permResults$diff<-with(permResults,paste(sprintf("%.3f",diff)," (",sprintf("%.3f",L),"-",sprintf("%.3f",U),")",sep=""))
permResults<-subset(permResults,select=c('pairs','F.Model','R2','p.value','diff','p adj'))
colnames(permResults)<-c('Contrast','F','R2','P1','Dispersion','P2')
permResults$F<-round(permResults$F,2)
permResults$R2<-round(permResults$R2,2)
permResults$P1<-round(permResults$P1,4)
permResults$P2<-round(permResults$P2,4)
permResults

#Save pairwise PERMANOVA table
write.csv(permResults,'Table_GroundPERMANOVA.csv',row.names=FALSE)

#Export species scores
SpecScores<-as.data.frame(MD1$species)
SpecScores$ACRON<-rownames(SpecScores)
SpecScores<-merge(SpecScores,traits,all.x=TRUE,all.y=FALSE)
SpecScores$PCH<-with(SpecScores,
                     ifelse(PHYS=='FORB',21,
                            ifelse(PHYS=='FERN',21,
                                   ifelse(PHYS=='SEDGE',21,
                                          ifelse(PHYS=='GRASS',21,25)))))
SpecScores$BG<-with(SpecScores,
                     ifelse(NT=='Ad','black',
                            ifelse(C<=3,'gray40',
                                   ifelse(C>=7,'white','gray80'))))

SpecScores$COL<-with(SpecScores,
                     ifelse(PHYS=='FORB','blue',
                            ifelse(PHYS=='FERN','red',
                                   ifelse(PHYS=='SEDGE','yellow',
                                          ifelse(PHYS=='GRASS','green',
                                                 ifelse(PHYS=='TREE','pink',
                                                        ifelse(PHYS=='SHRUB','seagreen','gray80')))))))

#write.csv(SpecScores,"Table_NMDSGroundSpecScores.csv",row.names=TRUE)

#Extract site scores for plotting
MDsite<-as.data.frame(MD1$points)
MDsite$HY<-paste(substr(rownames(MDsite),6,6),substr(rownames(MDsite),3,4),sep='')

#Calculate site centroids for each habitat_x_year
MDcen<-setNames(aggregate(list(MDsite$MDS1,MDsite$MDS2,MDsite$MDS3),by=list(MDsite$HY),mean),c('HY','MDS1','MDS2','MDS3'))
MDcen$H<-substr(MDcen$HY,1,1)
MDcen$Y<-substr(MDcen$HY,2,3)

#Set plot controls
plotCols<-c('darkgreen','red','black')

#Plot NMDS result
tiff('Fig_GroundNMDSBray_Pub_v23.tif',res=300,height=80,width=160,units='mm')
par(mfrow=c(1,2),mai=c(0,0.4,0.3,0),omi=c(1,0.3,0,1.2),mgp=c(1.5,0.3,0),
    tcl=-0.2,las=1,cex.axis=0.8)
plot(MD1,choices=c(1,2),type='n',ylab="",xlab="",xlim=c(-1.1,1.4))
#abline(h=0,lty=2,col='gray70')
#abline(v=0,lty=2,col='gray70')
with(MDcen[MDcen$Y=='01',],arrows(MDS1[1],MDS2[1],MDS1[3],MDS2[3],length=0,lty=3,col='gray50'))
with(MDcen[MDcen$Y=='01',],arrows(MDS1[3],MDS2[3],MDS1[2],MDS2[2],length=0,lty=3,col='gray50'))
with(MDcen[MDcen$Y=='01',],arrows(MDS1[2],MDS2[2],1.0,MDS2[2],length=0,lty=3,col='gray50'))
text(1.0,MDcen[MDcen$Y=='01',]$MDS2[2],'2001',adj=0,col='gray50',cex=0.8)
with(MDcen[MDcen$Y=='05',],arrows(MDS1[1],MDS2[1],MDS1[3],MDS2[3],length=0,lty=3,col='gray50'))
with(MDcen[MDcen$Y=='05',],arrows(MDS1[3],MDS2[3],MDS1[2],MDS2[2],length=0,lty=3,col='gray50'))
with(MDcen[MDcen$Y=='05',],arrows(MDS1[2],MDS2[2],1.0,MDS2[2],length=0,lty=3,col='gray50'))
text(1.0,MDcen[MDcen$Y=='05',]$MDS2[2],'2005',adj=0,col='gray50',cex=0.8)
with(MDcen[MDcen$Y=='12',],arrows(MDS1[1],MDS2[1],MDS1[3],MDS2[3],length=0,lty=3,col='gray50'))
with(MDcen[MDcen$Y=='12',],arrows(MDS1[3],MDS2[3],MDS1[2],MDS2[2],length=0,lty=3,col='gray50'))
with(MDcen[MDcen$Y=='12',],arrows(MDS1[2],MDS2[2],1.0,MDS2[2],length=0,lty=3,col='gray50'))
text(1.0,MDcen[MDcen$Y=='12',]$MDS2[2],'2012',adj=0,col='gray50',cex=0.8)
for(i in 1:11){with(MDcen[MDcen$H=='J',],Arrows(MDS1[i],MDS2[i],MDS1[i+1],MDS2[i+1],arr.type='triangle',arr.length=0,col=DBWpalCOL[1]))}
for(i in 1:11){with(MDcen[MDcen$H=='J',],points(MDS2[i]~MDS1[i],pch=21,bg=DBWpalBG[1],col=DBWpalCOL[1]))}
for(i in 1:11){with(MDcen[MDcen$H=='W',],Arrows(MDS1[i],MDS2[i],MDS1[i+1],MDS2[i+1],arr.type='triangle',arr.length=0,col=DBWpalCOL[2]))}
for(i in 1:11){with(MDcen[MDcen$H=='W',],points(MDS2[i]~MDS1[i],pch=24,bg=DBWpalBG[2],col=DBWpalCOL[2]))}
for(i in 1:11){with(MDcen[MDcen$H=='R',],Arrows(MDS1[i],MDS2[i],MDS1[i+1],MDS2[i+1],arr.type='triangle',arr.length=0,col=DBWpalCOL[3]))}
for(i in 1:11){with(MDcen[MDcen$H=='R',],points(MDS2[i]~MDS1[i],pch=22,bg=DBWpalBG[3],col=DBWpalCOL[3]))}
text(-0.4,0.4,'Juniper',adj=1,cex=0.8,col=DBWpalCOL[1])
text(0.5,-0.45,'Forest',adj=0,cex=0.8,col=DBWpalCOL[3])
text(-0.2,-0.6,'Woodland',adj=0,cex=0.8,col=DBWpalCOL[2])
mtext(side=3,adj=0.03,line=-1.2,'A')
mtext(side=2,line=2,'NMDS Axis 2',las=3,cex=0.8)
mtext(side=1,line=2,'NMDS Axis 1',cex=0.8)

plot(MD1,choices=c(1,2),type='n',ylab="",xlab="",xlim=c(-1.1,1.4))
#abline(h=0,lty=2,col='gray70')
#abline(v=0,lty=2,col='gray70')
with(SpecScores,points(MDS2~MDS1,pch=PCH,bg=BG,col='black'))
mtext(side=3,adj=0.03,line=-1.2,'B')
mtext(side=1,line=2,'NMDS Axis 1',cex=0.8)

par(xpd=NA)
legend(1.5,1.2,legend=c('Woody','Herbaceous'),pch=c(25,21),bty='n',title='Physiognamy',
       cex=0.8,pt.cex=1)
legend(1.5,0.1,legend=c('Non-native','Ruderal','Matrix','Conservative'),
       fill=c('black','gray40','gray80','white'),bty='n',title='Conservatism',
       cex=0.8,pt.cex=1)
par(xpd=FALSE)

dev.off()


#Plot NMDS result - full figure for supplement
tiff('Fig_GroundNMDSBray_SM_v23.tif',res=300,height=150,width=160,units='mm')
par(mfrow=c(2,2),mai=c(0,0.4,0.3,0),omi=c(1,0.3,0,1.2),mgp=c(1.5,0.3,0),
    tcl=-0.2,las=1,cex.axis=0.8)

plot(MD1,choices=c(1,2),type='n',ylab="",xlab="",xlim=c(-1.1,1.4))
with(MDcen[MDcen$Y=='01',],arrows(MDS1[1],MDS2[1],MDS1[3],MDS2[3],length=0,lty=3,col='gray50'))
with(MDcen[MDcen$Y=='01',],arrows(MDS1[3],MDS2[3],MDS1[2],MDS2[2],length=0,lty=3,col='gray50'))
with(MDcen[MDcen$Y=='01',],arrows(MDS1[2],MDS2[2],1.0,MDS2[2],length=0,lty=3,col='gray50'))
text(1.0,MDcen[MDcen$Y=='01',]$MDS2[2],'2001',adj=0,col='gray50',cex=0.8)
with(MDcen[MDcen$Y=='05',],arrows(MDS1[1],MDS2[1],MDS1[3],MDS2[3],length=0,lty=3,col='gray50'))
with(MDcen[MDcen$Y=='05',],arrows(MDS1[3],MDS2[3],MDS1[2],MDS2[2],length=0,lty=3,col='gray50'))
with(MDcen[MDcen$Y=='05',],arrows(MDS1[2],MDS2[2],1.0,MDS2[2],length=0,lty=3,col='gray50'))
text(1.0,MDcen[MDcen$Y=='05',]$MDS2[2],'2005',adj=0,col='gray50',cex=0.8)
with(MDcen[MDcen$Y=='12',],arrows(MDS1[1],MDS2[1],MDS1[3],MDS2[3],length=0,lty=3,col='gray50'))
with(MDcen[MDcen$Y=='12',],arrows(MDS1[3],MDS2[3],MDS1[2],MDS2[2],length=0,lty=3,col='gray50'))
with(MDcen[MDcen$Y=='12',],arrows(MDS1[2],MDS2[2],1.0,MDS2[2],length=0,lty=3,col='gray50'))
text(1.0,MDcen[MDcen$Y=='12',]$MDS2[2],'2012',adj=0,col='gray50',cex=0.8)
for(i in 1:11){with(MDcen[MDcen$H=='J',],Arrows(MDS1[i],MDS2[i],MDS1[i+1],MDS2[i+1],arr.type='triangle',arr.length=0,col=DBWpalCOL[1]))}
for(i in 1:11){with(MDcen[MDcen$H=='J',],points(MDS2[i]~MDS1[i],pch=21,bg=DBWpalBG[1],col=DBWpalCOL[1]))}
for(i in 1:11){with(MDcen[MDcen$H=='W',],Arrows(MDS1[i],MDS2[i],MDS1[i+1],MDS2[i+1],arr.type='triangle',arr.length=0,col=DBWpalCOL[2]))}
for(i in 1:11){with(MDcen[MDcen$H=='W',],points(MDS2[i]~MDS1[i],pch=24,bg=DBWpalBG[2],col=DBWpalCOL[2]))}
for(i in 1:11){with(MDcen[MDcen$H=='R',],Arrows(MDS1[i],MDS2[i],MDS1[i+1],MDS2[i+1],arr.type='triangle',arr.length=0,col=DBWpalCOL[3]))}
for(i in 1:11){with(MDcen[MDcen$H=='R',],points(MDS2[i]~MDS1[i],pch=22,bg=DBWpalBG[3],col=DBWpalCOL[3]))}
text(-0.4,0.4,'Juniper',adj=1,cex=0.8,col=DBWpalCOL[1])
text(0.5,-0.45,'Forest',adj=0,cex=0.8,col=DBWpalCOL[3])
text(-0.2,-0.6,'Woodland',adj=0,cex=0.8,col=DBWpalCOL[2])
mtext(side=3,adj=0.03,line=-1.2,'A')
mtext(side=2,line=2,'NMDS Axis 2',las=3,cex=0.8)

plot(MD1,choices=c(1,2),type='n',ylab="",xlab="",xlim=c(-1.1,1.4))
orditorp(MD1, display = "species", choices=c(1,2), priority = priSpp,pcol='gray70',cex=0.8)
mtext(side=3,adj=0.01,line=-1.2,'B')

plot(MD1,choices=c(1,3),type='n',xlab="",ylab="",xlim=c(-1.1,1.4))
for(i in 1:11){with(MDcen[MDcen$H=='J',],Arrows(MDS1[i],MDS3[i],MDS1[i+1],MDS3[i+1],arr.type='triangle',arr.length=0,col=DBWpalCOL[1]))}
for(i in 1:11){with(MDcen[MDcen$H=='J',],points(MDS3[i]~MDS1[i],pch=21,bg=DBWpalBG[1],col=DBWpalCOL[1]))}
for(i in 1:11){with(MDcen[MDcen$H=='R',],Arrows(MDS1[i],MDS3[i],MDS1[i+1],MDS3[i+1],arr.type='triangle',arr.length=0,col=DBWpalCOL[3]))}
for(i in 1:11){with(MDcen[MDcen$H=='R',],points(MDS3[i]~MDS1[i],pch=22,bg=DBWpalBG[3],col=DBWpalCOL[3]))}
for(i in 1:11){with(MDcen[MDcen$H=='W',],Arrows(MDS1[i],MDS3[i],MDS1[i+1],MDS3[i+1],arr.type='triangle',arr.length=0.2,col=DBWpalCOL[2]))}
for(i in 1:11){with(MDcen[MDcen$H=='W',],points(MDS3[i]~MDS1[i],pch=24,bg=DBWpalBG[2],col=DBWpalCOL[2]))}
mtext(side=3,adj=0.03,line=-1.2,'C')
mtext(side=1,line=2,'NMDS Axis 1',cex=0.8)
mtext(side=2,line=2,'NMDS Axis 3',las=3,cex=0.8)

plot(MD1,choices=c(1,3),type='n',ylab="",xlim=c(-1.1,1.4),xlab='')
#abline(h=0,lty=2,col='gray70')
#abline(v=0,lty=2,col='gray70')
orditorp(MD1, display = "species", choices=c(1,3), priority = priSpp,pcol='gray70',cex=0.8)
mtext(side=3,adj=0.03,line=-1.2,'D')
mtext(side=1,line=2,'NMDS Axis 1',cex=0.8)

dev.off()

###############################################################################  
#TABLE SHOWING INCREASERS, DECREASERS, AND STABLE SPECIES
###############################################################################

vIDS<-veg[veg$idSpecLevel=='1',]#Select fully identified species
vIDS_aggStems<-setNames(aggregate(vIDS$Stems,by=list(vIDS$Species,vIDS$Y),sum),c('Species','Y','Stems')) #Aggregate stem count
vIDS_aggStems<-subset(vIDS_aggStems,Y=='2001' | Y>2009) #Select only 2001, 2010, 2011, and 2012
vIDS_aggStems<-setNames(dcast(vIDS_aggStems,Species~Y,sum),c('Species','Y01','Y10','Y11','Y12')) #Sum stems
vIDS_aggStems$lastThree<-round(with(vIDS_aggStems,(Y10+Y11+Y12)/3),1) #Take mean of last three years
vIDS_aggStems$Delta<-vIDS_aggStems$lastThree-vIDS_aggStems$Y01 #Get change over time 01-last three
vx<-vIDS_aggStems #Give a short name for ease below

#New species
newSpecies<-vx[vx$Y01==0,][with(vx[vx$Y01==0,],order(-lastThree)),]
newSpecies<-merge(newSpecies,traits,by='Species')
newSpecies<-merge(newSpecies,pheno,by='Species')
write.csv(newSpecies,'Table_newSpecies.csv',row.names=FALSE)

#Species that disappeared
disappeared<-vx[vx$Y01>0 & vx$Y10<1 & vx$Y11<1 & vx$Y12<1,]
disappeared<-disappeared[with(disappeared,order(Delta)),]
disappeared<-merge(disappeared,traits,by='Species')
disappeared<-merge(disappeared,pheno,by='Species')
write.csv(disappeared,'Table_disappeared.csv',row.names=FALSE)  

#Biggest increasers
increasers<-head(vx[with(vx,order(-Delta)),],20)
increasers<-merge(increasers,traits,by='Species')
increasers<-merge(increasers,pheno,by='Species')
write.csv(increasers,'Table_increasers.csv',row.names=FALSE)

#Biggest decreasers
decreasers<-head(vx[with(vx,order(Delta)),],20)
decreasers<-merge(decreasers,traits,by='Species')
decreasers<-merge(decreasers,pheno,by='Species')
write.csv(decreasers,'Table_decreasers.csv',row.names=FALSE)

#Most stable species
stable<-vx[with(vx,order(abs(Delta))),]
stable<-subset(stable,Y01>=10)
stable<-head(stable,20)
stable<-merge(stable,traits,by='Species')
stable<-merge(stable,pheno,by='Species')
write.csv(stable,'Table_stable.csv',row.names=FALSE)

vx2<-merge(vx,traits,by='Species')
vx2<-merge(vx2,pheno,by='Species')
write.csv(vx2,'Table_IDS.csv',row.names=FALSE)

