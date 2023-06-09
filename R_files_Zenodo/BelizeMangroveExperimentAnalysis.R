
# Edit Date: 4/20/22
# Author: A.Forde



# ---------------------------------------------------------------------------- #
# SETTING UP SESSION
# ---------------------------------------------------------------------------- #

Setup <- function(){}

# clear all existing data 
rm(list=ls())

# setwd("~/Desktop/WD")

# Load packages

require(DHARMa)

require(contrast)
require(effects)

require(VGAM)
require(lavaan)

require(plyr)
require(reshape2)

require(lme4)
require(nlme)
require(AICcmodavg)

require(MASS)



# ---------------------------------------------------------------------------- #
# PROCESSING DATA
# ---------------------------------------------------------------------------- #

ProcessingData <- function(){}


# --------------------------------------------------------- #
# Load main experimental dataset

# Create dataframe Ca (all data)
# Create dataframe Cb (only 54 weeks)

Ca<-read.csv("BPTreeData.csv")

Tr<-read.csv("BPTreeAttributes.csv")
Tr<-Tr[,1:4]
names(Tr)[2:4]<-c("S","P","N")
Tr$N<-as.factor(Tr$N)

Ca<-merge(Tr,Ca,by="Tree",all=TRUE)

names(Ca)[5]<-"Wks"

# Create variable for arthropod herbivore abundance (PE)

Ca$PE <- Ca$Lep.Borer + Ca$Aratus

Ca$PE <- Ca$PE + grepl("cricket",Ca$Notes)
Ca$PE <- Ca$PE + grepl("katydid",Ca$Notes)

Ca$PE <- Ca$PE + grepl("whitefly",Ca$Notes)
Ca$PE <- Ca$PE + grepl("2 whitefly",Ca$Notes)
Ca$PE <- Ca$PE + grepl("bagworm",Ca$Notes)
Ca$PE <- Ca$PE + grepl("cerambicid",Ca$Notes)

Ca$PE <- Ca$PE + grepl("thrips",Ca$Notes)
Ca$PE <- Ca$PE + grepl("2 thrips",Ca$Notes)

# Remove unused columns and sort

Ca<-subset(Ca,select=-c(Date,Cocoons,Littorina,Notes))

Ca<-Ca[order(Ca$Wks,Ca$Tree),]

# Create variable for number of live terminal shoots on plants

Ca$TotBud<-Ca$LiveBuds+Ca$DamagBuds
Ca$TotBud[Ca$Wks==0]<-Ca$LiveBuds[Ca$Wks==0]

Cb<-subset(Ca,subset=(Ca$Wks==54))


# --------------------------------------------------------- #
# Processing Floc

# Create df floc with an estimate of floc depth for each tree x date

Ca$MFloc<-NA

#Only one value recorded at time=0, so this is set to the "mean" for this time
#All other times have 5 replicate measurements per tree

Ca$MFloc[Ca$Wks==0]<-Ca$Fc1[Ca$Wks==0]

Ca$MFloc[Ca$Wks!=0]<-(Ca$Fc1[Ca$Wks!=0]+Ca$Fc2[Ca$Wks!=0]+Ca$Fc3[Ca$Wks!=0]
                      +Ca$Fc4[Ca$Wks!=0]+Ca$Fc5[Ca$Wks!=0])/5

d<-ddply(.data=Ca,.variables=c("Tree","S","P"),.fun="summarise",m=mean(MFloc,na.rm=TRUE))

names(d)[4]<-"Floc"

d$LF<-log(d$Floc)

floc <-d


# --------------------------------------------------------- #
# Deciding on start and end values for growth variables for each marked twig

#Measures other than leaf num
BB<-read.csv("BPBranchData.csv")

BB$TwigID<-paste(BB$Tree,BB$Twig,sep=".")

#Twig 7.2 and 38.3 do not have Jan (week0) values
BB[BB$Date=="JAN",]

#But they both have Mar (week10) values (so wk10 will be their startval)
BB[BB$Date=="MAR",]

#Leaf num was not measured in Jan (week 0)

#All twigs except 2.3 have March values
#2.3 leaf num was not accurately measured until Sept (week33), so it will 
#be excluded
BB[BB$Date=="MAR",]

#10.3, 12.1, 13.2 and all 3 40's do not have Feb (week 54) values
ddply(BB[BB$Date=="FEB",],.variables="TwigID",.fun="summarise",S=sum(Leaves))

#They all have Sept values though, except for 13.2 which died before 20 weeks
#into the experiment. Exclude 13.2 but use Sep endvalues for others
ddply(BB[BB$Date=="SEPT",],.variables="TwigID",.fun="summarise",S=sum(Leaves))


# --------------------------------------------------------- #
# Create derived dataframe with summary values for the earliest and latest measurements for each branch for total length, main axis length, total (new) shoots, shoots/mainlength, total leaves, nodes on main axis

#Create empty dataframe with rows corresponding to branches (3/tree)
Tree<-1:40
Twig<-1:3
x<-expand.grid(Twig,Tree)
G<-data.frame(x[2],x[1])
names(G)<-c("Tree","Twig")
G$TwigID<-paste(G$Tree,G$Twig,sep=".")

#Add data on difference in time between earliest and latest
#growth measures on each tree (not leaves)
G$Time1<-54
G$Time1[G$TwigID==7.2|G$TwigID==38.3]<-44
G$Time1[G$TwigID==10.3|G$TwigID==12.1|G$TwigID==40.1|G$TwigID==40.2|G$TwigID==40.3]<-33

#Add data on difference in time between earliest and latest
#LEAF number measures on each tree
G$Time2<-44
G$Time2[G$TwigID==10.3|G$TwigID==12.1|G$TwigID==40.1|G$TwigID==40.2|G$TwigID==40.3]<-23

#Add data on starting length of each twig which is both the 
#Starting total and "main axis" length
#(First measures are in Jan except for branches 7.2 and 38.3) 
#(col 6)
x1<-ddply(BB[BB$Date=="JAN",],.variables=c("TwigID"),.fun="summarise",S=sum(Length))

G<-merge(G,x1,by="TwigID",all=TRUE)
names(G)[6]<-"Tlen1"

G$Tlen1[G$TwigID==7.2]<-BB$Length[BB$Date=="MAR"&BB$TwigID==7.2]
G$Tlen1[G$TwigID==38.3]<-BB$Length[BB$Date=="MAR"&BB$TwigID==38.3]

#Add data on starting num leaves on each twig which is both the 
#Starting total and "main axis" leaf number
#(First measures are in Mar except for branches 7.2 and 38.3)
#(col 7) 
x1<-ddply(BB[BB$Date=="MAR",],.variables=c("TwigID"),.fun="summarise",S=sum(Leaves))

G<-merge(G,x1,by="TwigID",all=TRUE)
names(G)[7]<-"TLv1"


#Add data on ending total length (col 8)
x1<-ddply(BB[BB$Date=="FEB",],.variables=c("TwigID"),.fun="summarise",S=sum(Length))

G<-merge(G,x1,by="TwigID",all=TRUE)
names(G)[8]<-"Tlen2"

x1<-ddply(BB[BB$Date=="SEPT",],.variables=c("TwigID"),.fun="summarise",S=sum(Length))
names(x1)[2]<-"Length"

G$Tlen2[G$TwigID==10.3]<-x1$Length[x1$TwigID==10.3]
G$Tlen2[G$TwigID==12.1]<-x1$Length[x1$TwigID==12.1]
G$Tlen2[G$TwigID==40.1]<-x1$Length[x1$TwigID==40.1]
G$Tlen2[G$TwigID==40.2]<-x1$Length[x1$TwigID==40.2]
G$Tlen2[G$TwigID==40.3]<-x1$Length[x1$TwigID==40.3]


#Add data on ending total number of leaves (col 9)
x1<-ddply(BB[BB$Date=="FEB",],.variables=c("TwigID"),.fun="summarise",S=sum(Leaves))

G<-merge(G,x1,by="TwigID",all=TRUE)
names(G)[9]<-"TLv2"

x1<-ddply(BB[BB$Date=="SEPT",],.variables=c("TwigID"),.fun="summarise",S=sum(Leaves))
names(x1)[2]<-"Leaves"

G$TLv2[G$TwigID==10.3]<-x1$Leaves[x1$TwigID==10.3]
G$TLv2[G$TwigID==12.1]<-x1$Leaves[x1$TwigID==12.1]
G$TLv2[G$TwigID==40.1]<-x1$Leaves[x1$TwigID==40.1]
G$TLv2[G$TwigID==40.2]<-x1$Leaves[x1$TwigID==40.2]
G$TLv2[G$TwigID==40.3]<-x1$Leaves[x1$TwigID==40.3]


#Add data on ending number of nodes on main axis (col 10)
x1<-BB[BB$TwigType=="M"&BB$Date=="FEB",c("TwigID","Nodes")]

G<-merge(G,x1,by="TwigID",all=TRUE)
names(G)[10]<-"Mnod2"

x1<-subset(BB,subset=(Date==c("SEPT")&BB$TwigType==c("M")))

G$Mnod2[G$TwigID==10.3]<-x1$Nodes[x1$TwigID==10.3]
G$Mnod2[G$TwigID==12.1]<-x1$Nodes[x1$TwigID==12.1]
G$Mnod2[G$TwigID==40.1]<-x1$Nodes[x1$TwigID==40.1]
G$Mnod2[G$TwigID==40.2]<-x1$Nodes[x1$TwigID==40.2]
G$Mnod2[G$TwigID==40.3]<-x1$Nodes[x1$TwigID==40.3]


#Add data on ending main axis length (col 11)
x1<-BB[BB$TwigType=="M"&BB$Date=="FEB",c("TwigID","Length")]
G<-merge(G,x1,by="TwigID",all=TRUE)
names(G)[11]<-"Mlen2"

x1<-subset(BB,subset=(Date==c("SEPT")&BB$TwigType==c("M")))

G$Mlen2[G$TwigID==10.3]<-x1$Length[x1$TwigID==10.3]
G$Mlen2[G$TwigID==12.1]<-x1$Length[x1$TwigID==12.1]
G$Mlen2[G$TwigID==40.1]<-x1$Length[x1$TwigID==40.1]
G$Mlen2[G$TwigID==40.2]<-x1$Length[x1$TwigID==40.2]
G$Mlen2[G$TwigID==40.3]<-x1$Length[x1$TwigID==40.3]


#Add data on ending main axis leaf number (col 12)
x1<-BB[BB$TwigType=="M"&BB$Date=="FEB",c("TwigID","Leaves")]
G<-merge(G,x1,by="TwigID",all=TRUE)
names(G)[12]<-"MLv2"

x1<-subset(BB,subset=(Date==c("SEPT")&BB$TwigType==c("M")))

G$MLv2[G$TwigID==10.3]<-x1$Leaves[x1$TwigID==10.3]
G$MLv2[G$TwigID==12.1]<-x1$Leaves[x1$TwigID==12.1]
G$MLv2[G$TwigID==40.1]<-x1$Leaves[x1$TwigID==40.1]
G$MLv2[G$TwigID==40.2]<-x1$Leaves[x1$TwigID==40.2]
G$MLv2[G$TwigID==40.3]<-x1$Leaves[x1$TwigID==40.3]


#Add data on ending total number of shoots (including main) (col 13)
x1<-ddply(BB[BB$Date=="FEB",],.variables=c("TwigID"),.fun="summarise",L=length(Length))

G<-merge(G,x1,by="TwigID",all=TRUE)
names(G)[13]<-"Shoots"

x1<-ddply(BB[BB$Date=="SEPT",],.variables=c("TwigID"),.fun="summarise",L=length(Length))
names(x1)[2]<-"Shoots"

G$Shoots[G$TwigID==10.3]<-x1$Shoots[x1$TwigID==10.3]
G$Shoots[G$TwigID==12.1]<-x1$Shoots[x1$TwigID==12.1]
G$Shoots[G$TwigID==40.1]<-x1$Shoots[x1$TwigID==40.1]
G$Shoots[G$TwigID==40.2]<-x1$Shoots[x1$TwigID==40.2]
G$Shoots[G$TwigID==40.3]<-x1$Shoots[x1$TwigID==40.3]

x1<-ddply(BB[BB$Date=="FEB",],.variables=c("TwigID"),.fun="summarise",S=sum(Nodes))


#Add data on ending total number of "L" (notL1orL2) (primary) shoots (col 14)
x1<-ddply(BB[BB$Date=="FEB"&(BB$TwigType=="L"|BB$TwigType=="M"),],.variables=c("TwigID"),.fun="summarise",L=length(Length))

G<-merge(G,x1,by="TwigID",all=TRUE)
names(G)[14]<-"PrimShoot"

x1<-ddply(BB[BB$Date=="SEPT"&(BB$TwigType=="L"|BB$TwigType=="M")&!is.na(BB$TwigType),],.variables=c("TwigID"),.fun="summarise",L=length(Length))
names(x1)[2]<-"Shoots"

G$PrimShoot[G$TwigID==10.3]<-x1$Shoots[x1$TwigID==10.3]
G$PrimShoot[G$TwigID==12.1]<-x1$Shoots[x1$TwigID==12.1]
G$PrimShoot[G$TwigID==40.1]<-x1$Shoots[x1$TwigID==40.1]
G$PrimShoot[G$TwigID==40.2]<-x1$Shoots[x1$TwigID==40.2]
G$PrimShoot[G$TwigID==40.3]<-x1$Shoots[x1$TwigID==40.3]

G$PrimShoot<-G$PrimShoot-1

#Add data on whether main axis meristem is dead (col 15)

x1<-BB[BB$TwigType=="M"&BB$Date=="SEPT"&!is.na(BB$TwigType),c("TwigID","DeadMeri")]
G<-merge(G,x1,by="TwigID",all=TRUE)
names(G)[15]<-"Dead"

# branches 5.2,25.1,27.2 don't have data for wwk 33 but M were seen alive at 54 weeks
G$Dead[is.na(G$Dead)]<-0


#Add data on initial num nodes on main axis (col16)
x1<-BB[BB$TwigType=="M"&BB$Date=="JAN",c("TwigID","Nodes")]

G<-merge(G,x1,by="TwigID",all=TRUE)
names(G)[16]<-"Mnod1"

G$Mnod1[G$TwigID==7.2]<-BB$Nodes[BB$Date=="MAR"&BB$TwigID==7.2]
G$Mnod1[G$TwigID==38.3]<-BB$Nodes[BB$Date=="MAR"&BB$TwigID==38.3]

G<-subset(G,subset=!is.na(TwigID))

#-------------#
#Add treatment data!

CT<-read.csv("BPTreeAttributes.csv")
names(CT)[2:4]<-c("S","P","N")
CT<-CT[,1:4]

Dat<-merge(CT,G,by="Tree",all=TRUE)

Dat<-subset(Dat,subset=TwigID!=13.2)

#-------------#
#Dat contains growth measures at the level of branches
#Mdat contains growth measures at the level of trees
#both Dat and Mdat contain data on starting and ending conditions

#Calculate number of NEW shoots (i.e. shoots-1)
Dat$NewShoots<-Dat$Shoots-1

#Create variable for rate of total shoots produced per month
Dat$ShootRate<-(Dat$NewShoots)/((Dat$Time1*7)/30.43)

#Create a variable for rate of primary shoots (produced by main axis)
#per month. NA given for main axes that had a dead meristem in Sept (wk 33)
Dat$PShootRate<-(Dat$PrimShoot)/((Dat$Time1*7)/30.43)
Dat$PShootRate[Dat$Dead==1]<-NA

#Create variable for rate of change in total length of branches 
#(all segments) per month
Dat$Dtlen<-(Dat$Tlen2-Dat$Tlen1)/((Dat$Time1*7)/30.43)

#Create variable for rate of change in length of branch main axis 
#per month. NA given for main axes that had a dead meristem in Sept (wk 33)
Dat$Dmlen<-(Dat$Mlen2-Dat$Tlen1)/((Dat$Time1*7)/30.43)
Dat$Dmlen[Dat$Dead==1]<-NA

#Create variable for rate of change in the total num leaves per month
Dat$DtLv<-(Dat$TLv2-Dat$TLv1)/((Dat$Time2*7)/30.43)

#Create variable for rate of change in the num leaves on main axis per month
#NA given for main axes that had a dead meristem in Sept (wk 33)
Dat$DmLv<-(Dat$MLv2-Dat$TLv1)/((Dat$Time2*7)/30.43)
Dat$DmLv[Dat$Dead==1]<-NA

#Create variable for final internode length of main axis. NA given
#to branches with dead main axis meristems in Sept (wk 33)
Dat$MIN<-Dat$Mlen2/Dat$Mnod2
Dat$MIN[Dat$Dead==1]<-NA

#Create variable for total shoots per total end length
Dat$SPT<-Dat$Shoots/(Dat$Tlen2)

#Create variable for total shoots per total change in length
Dat$SPdT<-Dat$Shoots/(Dat$Dtlen)

#Create variable for shoots on main axis only per end length
Dat$PSPM<-(Dat$PrimShoot+1)/(Dat$Mlen2)

#-------------#

# Summarise variables at the level of trees by taking means
# Create Tree level dataframe Mdat

Mdat <- ddply(Dat,.variables=c("Tree","S","P","N"),.fun="summarise",
              Dtl=mean(Dtlen,na.rm=TRUE),Dml=mean(Dmlen,na.rm=TRUE),
              DtLv=mean(DtLv,na.rm=TRUE),DmLv=mean(DmLv,na.rm=TRUE),
              DtShu=mean(ShootRate,na.rm=TRUE),DmShu=mean(PShootRate,na.rm=TRUE),
              MIN=mean(MIN,na.rm=TRUE),NewShoots=mean(NewShoots,na.rm=TRUE),
              SPT=mean(SPT,na.rm=TRUE), SPdT=mean(SPdT,na.rm=TRUE), 
              PSPM=mean(PSPM,na.rm=TRUE))

#Add floc data to Mdat frame

f<-subset(floc,select=c("Tree","Floc","LF"))
Mdat<-merge(Mdat,f,by="Tree")

#Create log transformed stem elong (mm per 6 months)

Mdat$LS<-log(Mdat$Dtl*6)


# --------------------------------------------------------- #
# Calculate change in bud number and bud damage and add to Mdat

Ca$Battack<-Ca$DamagBuds+Ca$DeformedLeaves

T0 <- subset(Ca,subset=Wks==0)
T0<-T0[,c("Tree","TotBud")]

T1 <- subset(Ca,subset=Wks==54)
T1<-T1[,c("Tree","TotBud","DeadBuds","Battack","DamagBuds")]

T40 <- subset(Ca,subset=Wks==33)
T40<-T40[,c("Tree","TotBud","DeadBuds","Battack","DamagBuds")]

BU<-merge(T0,T1,by="Tree",all=TRUE)

BU$Time<-54

BU$Time[BU$Tree==40]<-33
BU$Time[BU$Tree==12]<-33

BU$TotBud.y[BU$Tree==40]<-T40$TotBud[T40$Tree==40]
BU$DeadBuds[BU$Tree==40]<-T40$DeadBuds[T40$Tree==40]
BU$Battack[BU$Tree==40]<-T40$Battack[T40$Tree==40]
BU$DamagBuds[BU$Tree==40]<-T40$DamagBuds[T40$Tree==40]
BU$TotBud.y[BU$Tree==12]<-T40$TotBud[T40$Tree==12]
BU$DeadBuds[BU$Tree==12]<-T40$DeadBuds[T40$Tree==12]
BU$Battack[BU$Tree==12]<-T40$Battack[T40$Tree==12]
BU$DamagBuds[BU$Tree==12]<-T40$DamagBuds[T40$Tree==12]

BU$Change<-BU$TotBud.y-BU$TotBud.x

BU$Change<-BU$Change/((BU$Time*7)/30.44)

names(BU)[c(2,3,4)]<-c("StartBuds","FinalBuds","Bdead")

BU<-BU[,c("Tree","StartBuds","FinalBuds","Change","Bdead","DamagBuds","Battack","Time")]

#-------------#
# Add bud damage and bud change data to tree level dataset

Mdat <- merge(Mdat,BU,by="Tree",all=TRUE)

Mdat$N<-as.factor(Mdat$N)


# --------------------------------------------------------- #
# Creating growth PCA and add to Mdat

#Looking at correlations between growth measures, for use in PCA
#Dtl (change in tot length of marked branches), DtLv (change in leaves on marked branches), MIN (final internode len of main axis of marked branches), and Change (change in total live and damaged buds on tree)
#And creating PC1

Grow<-Mdat[,c("Dtl","DtLv","MIN","Change")]

cor(Grow)

Grow<-Mdat[,c("Dtl","DtLv","Change")]

cor(Grow)

pr<-prcomp(Mdat[,c("Dtl","DtLv","Change")],scale=TRUE)

Mdat$PC1<-pr$x[,1]

Grow<-Mdat[,c("Dtl","DtLv","Change","PC1")]

cor(Grow)


# --------------------------------------------------------- #
# Add leaf elemental data to Mdat

x<-read.csv("BPNutrients.csv")
names(x)[c(2,3,4)]<-c("PercentN","PercentC","PercentP")
Mdat<-merge(Mdat,x,by="Tree",all=T)

x$C<-x$PercentC/12.011
x$N<-x$PercentN/14.007
x$P<-x$PercentP/30.974

x$CN<-x$C/x$N
x$CP<-x$C/x$P

x<-subset(x,select=c(Tree,CN,CP))

Mdat<-merge(Mdat,x,by="Tree",all=T)




# ---------------------------------------------------------------------------- #
# RESULTS: Comparison of plant characteristics
# ---------------------------------------------------------------------------- #

ResultsPlantCharacteristics <- function(){}


# --------------------------------------------------------- #
# Floc depth

ddply(.data=Mdat,.variables="P",.fun="summarise",m=mean(Floc),s=sd(Floc))
t.test(Mdat$Floc[Mdat$P=="fringe"], Mdat$Floc[Mdat$P=="dwarf"], var.equal=FALSE)
t.test(Mdat$Floc[Mdat$N==1], Mdat$Floc[Mdat$N==0], var.equal=FALSE)

# --------------------------------------------------------- #
# Starting terminal shoots

ddply(.data=Mdat,.variables="P",.fun="summarise",m=mean(StartBuds),s=sd(StartBuds))
t.test(Mdat$StartBuds[Mdat$P=="fringe"], Mdat$StartBuds[Mdat$P=="dwarf"], var.equal=FALSE)
t.test(Mdat$StartBuds[Mdat$N==1], Mdat$StartBuds[Mdat$N==0], var.equal=FALSE)

# --------------------------------------------------------- #
# Ending terminal shoots

ddply(.data=Mdat,.variables="P",.fun="summarise",m=mean(FinalBuds),s=sd(FinalBuds))
t.test(Mdat$FinalBuds[Mdat$P=="fringe"], Mdat$FinalBuds[Mdat$P=="dwarf"], var.equal=FALSE)
t.test(Mdat$FinalBuds[Mdat$N==1], Mdat$FinalBuds[Mdat$N==0], var.equal=FALSE)

# --------------------------------------------------------- #
# Shoot number and biomass relationship

BC<-read.csv("BranchClipData.csv")

CT<-read.csv("CageTreatments.csv")
names(CT)[1:4]<-c("Tree","S","P","N")

CT<-CT[,1:4]

D<-merge(CT,BC,by="Tree")

D<-subset(D,subset=Wks==54)

D$TotalMass<-D$StemMass+D$LeafMass

#par(mfrow=c(1,2))
#plot(D$LiveBuds,D$LeafMass)
#plot(D$LiveBuds,D$TotalMass)

#plot(D$LeafMass[D$P=="dwarf"],D$LiveBuds[D$P=="dwarf"])
#plot(D$LeafMass[D$P=="fringe"],D$TotalBuds[D$P=="fringe"])

g1<-lm(LeafMass~LiveBuds*P,data=D)
g2<-lm(LeafMass~LiveBuds+P,data=D)
g3<-lm(LeafMass~LiveBuds,data=D)
g4<-lm(LeafMass~P,data=D)
g5<-lm(LeafMass~1,data=D)

L<-list()

L[[1]] <- g1
L[[2]] <- g2
L[[3]] <- g3
L[[4]] <- g4
L[[5]] <- g5

n<-1:5

aictab(L,n)

summary(g3)

anova(g1,g2,test="Chisq")

# plot(g3)

g1<-lm(TotalMass~LiveBuds*P,data=D)
g2<-lm(TotalMass~LiveBuds+P,data=D)
g3<-lm(TotalMass~LiveBuds,data=D)
g4<-lm(TotalMass~P,data=D)
g5<-lm(TotalMass~1,data=D)

L<-list()

L[[1]] <- g1
L[[2]] <- g2
L[[3]] <- g3
L[[4]] <- g4
L[[5]] <- g5

n<-1:5

aictab(L,n)

anova(g1,g2,test="Chisq")

summary(g3)



# --------------------------------------------------------- #
# Stem elongation

ddply(.data=Mdat,.variables="P",.fun="summarise",m=mean(Dtl),s=sd(Dtl))
t.test(Mdat$Dtl[Mdat$P=="fringe"], Mdat$Dtl[Mdat$P=="dwarf"], var.equal=FALSE)
t.test(Mdat$Dtl[Mdat$N==1], Mdat$Dtl[Mdat$N==0], var.equal=FALSE)

Mdat$DtlCY <- Mdat$Dtl/10*12


# --------------------------------------------------------- #
# Change in leaf abundance

ddply(.data=Mdat,.variables="P",.fun="summarise",m=mean(DtLv),s=sd(DtLv))
t.test(Mdat$DtLv[Mdat$P=="fringe"], Mdat$DtLv[Mdat$P=="dwarf"], var.equal=FALSE)
t.test(Mdat$DtLv[Mdat$N==1], Mdat$DtLv[Mdat$N==0], var.equal=FALSE)

# --------------------------------------------------------- #
# Leaf thickness

#Leaf area and drymass measured from 10 fringe and 10 dwarf branches at 
#the boa site; 1 leaf per branch, non-experimental trees in March 2010

BSA<-read.csv("BelizeSLA.csv")

BSA$Area.cm<-BSA$Area.mm/100

BSA$LMA<-BSA$Mass.g/(BSA$Area.cm)

ddply(.data=BSA,.variables="TreeType",.fun="summarise",m=mean(LMA),s=sd(LMA))
t.test(BSA$LMA[BSA$TreeType=="fringe"], BSA$LMA[BSA$TreeType=="dwarf"], var.equal=FALSE)


# --------------------------------------------------------- #
# Leaf toughness

#Leaf Toughness averaged over 3 leaves per experimental tree
#measured in Sep 2010

#Load toughness data, collected Sept 2010
BLT<-read.csv("BelizeLeafToughness.csv")

A<-ddply(.data=BLT,.variables=c("Tree"),.fun="summarise",m=mean(Toughness.N))
names(A)[2]<-"Tough"

Mdat <- merge(Mdat,A,all=T,by="Tree")

ddply(.data=Mdat,.variables="P",.fun="summarise",m=mean(Tough),s=sd(Tough))
t.test(Mdat$Tough[Mdat$P=="fringe"], Mdat$Tough[Mdat$P=="dwarf"], var.equal=FALSE)
t.test(Mdat$Tough[Mdat$N==1], Mdat$Tough[Mdat$N==0], var.equal=FALSE)


# --------------------------------------------------------- #
# Leaf %N

ddply(.data=Mdat,.variables="P",.fun="summarise",m=mean(PercentN,na.rm=TRUE),s=sd(PercentN,na.rm=TRUE))
t.test(Mdat$PercentN[Mdat$P=="fringe"], Mdat$PercentN[Mdat$P=="dwarf"], var.equal=FALSE)
t.test(Mdat$PercentN[Mdat$N==1], Mdat$PercentN[Mdat$N==0], var.equal=FALSE)


# --------------------------------------------------------- #
# Leaf %P

dp<-ddply(.data=Mdat,.variables="P",.fun="summarise",m=mean(PercentP,na.rm=TRUE),s=sd(PercentP,na.rm=TRUE))
t.test(Mdat$PercentP[Mdat$P=="fringe"], Mdat$PercentP[Mdat$P=="dwarf"], var.equal=FALSE)
t.test(Mdat$PercentP[Mdat$N==1], Mdat$PercentP[Mdat$N==0], var.equal=FALSE)

(dp$s[1]/sqrt(20))*2
(dp$s[1]/sqrt(19))*2


# --------------------------------------------------------- #
# Leaf C:N

ddply(.data=Mdat,.variables="P",.fun="summarise",m=mean(CN,na.rm=TRUE),s=sd(CN,na.rm=TRUE))
t.test(Mdat$CN[Mdat$P=="fringe"], Mdat$CN[Mdat$P=="dwarf"], var.equal=FALSE)
t.test(Mdat$CN[Mdat$N==1], Mdat$CN[Mdat$N==0], var.equal=FALSE)

# --------------------------------------------------------- #
# Leaf C:P

ddply(.data=Mdat,.variables="P",.fun="summarise",m=mean(CP,na.rm=TRUE),s=sd(CP,na.rm=TRUE))
t.test(Mdat$CP[Mdat$P=="fringe"], Mdat$CP[Mdat$P=="dwarf"], var.equal=FALSE)
t.test(Mdat$CP[Mdat$N==1], Mdat$CP[Mdat$N==0], var.equal=FALSE)


# --------------------------------------------------------- #
# Leaf N:P

Mdat$NP <- Mdat$PercentN/Mdat$PercentP

ddply(.data=Mdat,.variables="P",.fun="summarise",m=mean(NP,na.rm=TRUE),s=sd(NP,na.rm=TRUE))
t.test(Mdat$NP[Mdat$P=="fringe"], Mdat$NP[Mdat$P=="dwarf"], var.equal=FALSE)



# ---------------------------------------------------------------------------- #
# Evidence of eased nutrient limitation in floc trees due to floc deposits
# ---------------------------------------------------------------------------- #

ResultsNutrientLimitation <- function(){}


# --------------------------------------------------------- #
# Substrate and leaves phosphorous comparison

CN.S<-read.csv("CNP-substrates.csv")
CN.L<-read.csv("CNP-leaves.csv")

#boxplot(CN.S$phos~CN.S$Substrate)
#boxplot(CN.S$CPratio~CN.S$Substrate)

vf1 <- varIdent(form=~1|Substrate)
g1 <- gls(phos ~ Substrate  , data=CN.S, weights=vf1,method="ML")
g2 <- gls(phos ~ Substrate  , data=CN.S, method="ML")

anova(g1,g2)

summary(g1)

g1 <- gls(phos ~ Substrate  , data=CN.S, weights=vf1,method="ML")
g3 <- gls(phos ~ 1  , data=CN.S, weights=vf1,method="ML")

anova(g1,g3)

t.test(CN.S$phos[CN.S$Substrate=="Floc"],CN.S$phos[CN.S$Substrate=="Peat"])



# --------------------------------------------------------- #
# Comparison of stem elong and % for our study vs. Feller 2003

ddply(.data=Mdat,.variables="P",.fun="summarise",m=mean(DtlCY),s=(sd(DtlCY)/sqrt(length(DtlCY)))*1.96)

ddply(.data=Mdat,.variables="P",.fun="summarise",m=mean(PercentP, na.rm=TRUE),s=(sd(PercentP, na.rm=TRUE)/sqrt(length(PercentP[!is.na(PercentP)])))*1.96)


# --------------------------------------------------------- #
# SEM

mod.3<-'PercentP~Floc + PC1
				PC1~Floc+P 
				'
mod.3.fit<-sem(mod.3, data=Mdat)

summary(mod.3.fit, standardized=TRUE, rsq=TRUE)



# ---------------------------------------------------------------------------- #
# Impacts of birds and tree types on arthropods, herbivory and growth
# ---------------------------------------------------------------------------- #

ResultsArthropodsHerbivoryGrowth <- function(){}

# --------------------------------------------------------- #
# Function for computing Mcfadden's pseudo R square used in following analyses

#Mcfadden's pseudo R square
#x = full model
#y = null model
#k = number of parameters

PRsq<- function (x, y, k ) {
a= 1- (logLik(x)/logLik(y))
b= 1- ((logLik(x)-k)/logLik(y))
cat ("Pseudo R Square Estimates", "\\n")
cat ("McFadden's:          ", a, "\\n")
cat ("adjusted McFadden's: ", b, "\\n")
}


# --------------------------------------------------------- #
# Arthropod density analysis 

Ca$Arth<- Ca$OtherInsects + Ca$Lep.Borer + Ca$Aratus + Ca$Ants + Ca$Spider

x1<-subset(Ca,subset=Wks==54)

x1$ADens<-x1$Arth/x1$TotBud

x1 <- subset(x1,select=c(Tree,TotBud,Arth,ADens))

Mdat <- merge(Mdat,x1,by="Tree",all=TRUE)

# log transform bud number

Mdat$Lbud <- log(Mdat$TotBud)

# Global model

g1<-glm(Arth~N*P+S+offset(Lbud),family="poisson",data=Mdat)

# Model Selection 

g2 <- update(g1, ~ N*P+offset(Lbud))
g3 <- update(g1, ~ N+P+S+offset(Lbud))
g4 <- update(g1, ~ N+P+offset(Lbud))
g5 <- update(g1, ~ N+S+offset(Lbud))
g6 <- update(g1, ~ P+S+offset(Lbud))
g7 <- update(g1, ~ N+offset(Lbud))
g8 <- update(g1, ~ P+offset(Lbud))
g9 <- update(g1, ~ S+offset(Lbud))
g10 <- update(g1, ~ 1+offset(Lbud))

L<-list()
L[[1]] <- g1
L[[2]] <- g2
L[[3]] <- g3
L[[4]] <- g4
L[[5]] <- g5
L[[6]] <- g6
L[[7]] <- g7
L[[8]] <- g8
L[[9]] <- g9
L[[10]] <- g10

names <- c("Birds x Planttype + Site","Birds x Planttype", "Birds + Planttype + Site", "Birds + Planttype", "Birds + Site", "Planttype + Site", "Birds", "Planttype", "Site", "Intercept")

aictab(L,names)

# Best Model

gbest<-glm(Arth~N+offset(Lbud),family="poisson",data=Mdat)

# Check assumptions

summary(gbest)
sum( resid(gbest, type="pearson")^2 ) / df.residual(gbest)
# plot(gbest)

# LRT hypothesis tests

g2<-glm(Arth~N+P+S+offset(Lbud),family="poisson",data=Mdat)

anova(g1,g2, test="Chisq")

g3<-glm(Arth~N+P+offset(Lbud),family="poisson",data=Mdat)
g4<-glm(Arth~N+S+offset(Lbud),family="poisson",data=Mdat)
g5<-glm(Arth~S+P+offset(Lbud),family="poisson",data=Mdat)

anova(g2,g3,test="Chisq")
anova(g2,g4,test="Chisq")
anova(g2,g5,test="Chisq")

gFull<-glm(Arth~N+offset(Lbud),family="poisson",data=Mdat)
gNull<-glm(Arth~offset(Lbud),family="poisson",data=Mdat)

PRsq(x=gFull,y=gNull,k=3)

# Differences between group means

Dat <- Mdat
Dat$Y <- Dat$ADens

d<-ddply(.data=Dat,.variables=c("N"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]


# --------------------------------------------------------- #
# Leaf damage analysis

# Read in leaf damage data

LD<-read.csv("BPLeafDmg.csv")

LD$LeafSize<-LD$Area.without.Damage

LD<-subset(LD,subset=!is.na(Damage))

S1<-Mdat[,c("Tree","P","N","PC1","LS","S")]

M<-merge(S1,LD,by="Tree")

head(M)

# Calculate mean leaf damage for each experimental unit

D<-ddply(.data=M,.variables="Tree",.fun="summarise",Dam=sum(Damage)/(length(Damage)))

summary(D$summarise)

table(M$Tree)

names(D)[2]<-"Dam"

D$DamPer <- D$Dam

# divide by 100 to change percentages to proportions

D$Dam<-D$Dam/100

# to calculate logit of damage proportion, add minimum value to all observations to remove zero's

D$DamLogit<-logitlink(D$Dam+min(D$Dam[D$Dam!=0],na.rm=TRUE))

Mdat<-merge(Mdat,D,by="Tree",all=TRUE)

# Global model

g1<-glm(DamLogit~N*P+S,family="gaussian",data=Mdat)

# plot(g1)

# Variances do not differ between bird treatments

vf1id<-varIdent(form= ~1 | N)
g1a<-gls(DamLogit~N*P+S,data=Mdat, method="ML")
g1b<-gls(DamLogit~N*P+S,data=Mdat, weights=vf1id, method="ML")

anova(g1a,g1b)

# Model Selection 

g2 <- update(g1, ~ N*P)
g3 <- update(g1, ~ N+P+S)
g4 <- update(g1, ~ N+P)
g5 <- update(g1, ~ N+S)
g6 <- update(g1, ~ P+S)
g7 <- update(g1, ~ N)
g8 <- update(g1, ~ P)
g9 <- update(g1, ~ S)
g10 <- update(g1, ~ 1)

L<-list()
L[[1]] <- g1
L[[2]] <- g2
L[[3]] <- g3
L[[4]] <- g4
L[[5]] <- g5
L[[6]] <- g6
L[[7]] <- g7
L[[8]] <- g8
L[[9]] <- g9
L[[10]] <- g10

names <- c("Birds x Planttype + Site","Birds x Planttype", "Birds + Planttype + Site", "Birds + Planttype", "Birds + Site", "Planttype + Site", "Birds", "Planttype", "Site", "Intercept")

aictab(L,names)

# Best Model

gbest<-glm(DamLogit~N+P,family="gaussian",data=Mdat)

summary(gbest)
# plot(gbest)

# LRT hypothesis tests

g2<-glm(DamLogit~N+P+S,family="gaussian",data=Mdat)
anova(g1,g2,test="Chisq")

g3<-glm(DamLogit~N+P,family="gaussian",data=Mdat)
g4<-glm(DamLogit~N+S,family="gaussian",data=Mdat)
g5<-glm(DamLogit~P+S,family="gaussian",data=Mdat)

anova(g2, g3, test="Chisq")
anova(g2, g4, test="Chisq")
anova(g2, g5, test="Chisq")

gFull<-glm(DamLogit~N+P,family="gaussian",data=Mdat)
gNull<-glm(DamLogit~1,family="gaussian",data=Mdat)

# Calculating McFadden R^2

PRsq(x=gFull,y=gNull, k=4)

# Differences between group means (Tree 10 is an outlier, according to diagnostics)

Dat <- Mdat
Dat$Y <- Dat$DamPer

d<-ddply(.data=Dat[Dat$Tree!=10,],.variables=c("N"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]

d<-ddply(.data=Dat[Dat$Tree!=10,],.variables=c("P"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]


# --------------------------------------------------------- #
# Bud damage analysis

# Global model

g1<-glm(Battack~N*P+S+offset(Lbud),family="poisson",data=Mdat)

# Model Selection 

g2 <- update(g1, ~ N*P+offset(Lbud))
g3 <- update(g1, ~ N+P+S+offset(Lbud))
g4 <- update(g1, ~ N+P+offset(Lbud))
g5 <- update(g1, ~ N+S+offset(Lbud))
g6 <- update(g1, ~ P+S+offset(Lbud))
g7 <- update(g1, ~ N+offset(Lbud))
g8 <- update(g1, ~ P+offset(Lbud))
g9 <- update(g1, ~ S+offset(Lbud))
g10 <- update(g1, ~ 1+offset(Lbud))

L<-list()
L[[1]] <- g1
L[[2]] <- g2
L[[3]] <- g3
L[[4]] <- g4
L[[5]] <- g5
L[[6]] <- g6
L[[7]] <- g7
L[[8]] <- g8
L[[9]] <- g9
L[[10]] <- g10

names <- c("Birds x Planttype + Site","Birds x Planttype", "Birds + Planttype + Site", "Birds + Planttype", "Birds + Site", "Planttype + Site", "Birds", "Planttype", "Site", "Intercept")

aictab(L,names)

# Best Model

gbest<-glm(Battack~N+P+S+offset(Lbud),family="poisson",data=Mdat)

# Best model overdispersion parameter is too high = 2.34

sum( resid(gbest, type="pearson")^2 ) / df.residual(gbest)

# Negative binomial alternative global model

g1<-glm.nb(Battack~N*P+S+offset(Lbud),data=Mdat)

# Model Selection 

g2 <- update(g1, ~ N*P+offset(Lbud))
g3 <- update(g1, ~ N+P+S+offset(Lbud))
g4 <- update(g1, ~ N+P+offset(Lbud))
g5 <- update(g1, ~ N+S+offset(Lbud))
g6 <- update(g1, ~ P+S+offset(Lbud))
g7 <- update(g1, ~ N+offset(Lbud))
g8 <- update(g1, ~ P+offset(Lbud))
g9 <- update(g1, ~ S+offset(Lbud))
g10 <- update(g1, ~ 1+offset(Lbud))

L<-list()
L[[1]] <- g1
L[[2]] <- g2
L[[3]] <- g3
L[[4]] <- g4
L[[5]] <- g5
L[[6]] <- g6
L[[7]] <- g7
L[[8]] <- g8
L[[9]] <- g9
L[[10]] <- g10

names <- c("Birds x Planttype + Site","Birds x Planttype", "Birds + Planttype + Site", "Birds + Planttype", "Birds + Site", "Planttype + Site", "Birds", "Planttype", "Site", "Intercept")

aictab(L,names)

# Best model overdispersion parameter better = 1.42

gbest <- g3

sum( resid(gbest, type="pearson")^2 ) / df.residual(gbest)
# plot(gbest)

# LRT hypothesis tests

g2 <- glm.nb(Battack~N+P+S+offset(Lbud),data=Mdat)
anova(g1,g2, test="Chisq")

g3 <- glm.nb(Battack~N+P+offset(Lbud),data=Mdat)
g4 <- glm.nb(Battack~N+S+offset(Lbud),data=Mdat)
g5 <- glm.nb(Battack~P+S+offset(Lbud),data=Mdat)

anova(g2,g3, test="Chisq")
anova(g2,g4, test="Chisq")
anova(g2,g5, test="Chisq")

gFULL <- glm.nb(Battack~N+P+S+offset(Lbud),data=Mdat)
gNULL <- glm.nb(Battack~1+offset(Lbud),data=Mdat)

# McFaddon R^2

PRsq(gFULL,gNULL,k=5)

# Differences between group means 

Dat <- Mdat
Dat$Y <- Dat$Battack/Dat$TotBud

d<-ddply(.data=Dat,.variables=c("N"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]

d<-ddply(.data=Dat,.variables=c("P"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]

d<-ddply(.data=Dat,.variables=c("S"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]


# --------------------------------------------------------- #
# PCA growth analysis

vf2id<-varIdent(form= ~1 | P)

# Global model (much more variation in Floc/Fringe plants)

g1<-gls(PC1~N*P+S,data=Mdat, weights=vf2id, method="ML")

# plot(g1)
# plot(predict(g1),resid(g1, type="pearson"))
# qqnorm(resid(g1, type="pearson"))
# qqline(resid(g1, type="pearson"))

# Model Selection 

g2 <- update(g1, ~ N*P)
g3 <- update(g1, ~ N+P+S)
g4 <- update(g1, ~ N+P)
g5 <- update(g1, ~ N+S)
g6 <- update(g1, ~ P+S)
g7 <- update(g1, ~ N)
g8 <- update(g1, ~ P)
g9 <- update(g1, ~ S)
g10 <- update(g1, ~ 1)

L<-list()
L[[1]] <- g1
L[[2]] <- g2
L[[3]] <- g3
L[[4]] <- g4
L[[5]] <- g5
L[[6]] <- g6
L[[7]] <- g7
L[[8]] <- g8
L[[9]] <- g9
L[[10]] <- g10

names <- c("Birds x Planttype + Site","Birds x Planttype", "Birds + Planttype + Site", "Birds + Planttype", "Birds + Site", "Planttype + Site", "Birds", "Planttype", "Site", "Intercept")

aictab(L,names)

# LRT hypothesis tests

g2<-gls(PC1~N+P+S,data=Mdat, weights=vf2id, method="ML")

anova(g1,g2)

g3<-gls(PC1~N+P,data=Mdat, weights=vf2id, method="ML")
g4<-gls(PC1~N+S,data=Mdat, weights=vf2id, method="ML")
g5<-gls(PC1~P+S,data=Mdat, weights=vf2id, method="ML")

anova(g2,g3)
anova(g2,g4)
anova(g2,g5)

gFull<-gls(PC1~P,data=Mdat, weights=vf2id, method="ML")
gNull<-gls(PC1~1,data=Mdat, weights=vf2id, method="ML")

# McFaddon's R^2

PRsq(gFull,gNull,k=3)

# Differences between group means 

Dat <- Mdat
Dat$Y <- Dat$PC1 + abs(min(Dat$PC1))

d<-ddply(.data=Dat,.variables=c("P"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]


# ---------------------------------------------------------------------------- #
# Supplemental analyses: Repeated arthropod surveys and Arthropod Herbivores
# ---------------------------------------------------------------------------- #

SupplementalArthAnalyses <- function(){}


# --------------------------------------------------------- #
# Repeated measures GLMM for all arthropods and all surveys

Dat <- subset(Ca,subset=Wks!=0,select=c("Tree","S","P","N","TotBud","Arth"))

Dat$ADens<-Dat$Arth/Dat$TotBud

# log transform bud number

Dat$Lbud <- log(Dat$TotBud)

# Global model

g1 <- glmer( Arth ~ N*P + S + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)

# Model selection

g2 <- glmer( Arth ~ N*P + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g3 <- glmer( Arth ~ N+P + S + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g4 <- glmer( Arth ~ N+P + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g5 <- glmer( Arth ~ N + S + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g6 <- glmer( Arth ~ P + S + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g7 <- glmer( Arth ~ N + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g8 <- glmer( Arth ~ P + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g9 <- glmer( Arth ~ S + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g10 <- glmer( Arth ~ 1 + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)

L<-list()
L[[1]] <- g1
L[[2]] <- g2
L[[3]] <- g3
L[[4]] <- g4
L[[5]] <- g5
L[[6]] <- g6
L[[7]] <- g7
L[[8]] <- g8
L[[9]] <- g9
L[[10]] <- g10
n<-1:10
aictab(L,n)

# Best model has poor general fit and overdispersion

so <- simulateResiduals(fittedModel = g4, n = 10000, refit=F)
testUniformity(simulationOutput = so, plot=FALSE)
testDispersion(simulationOutput = so, plot=FALSE, alternative="greater")

# Alternative negative binomial global model

g1 <- glmer.nb( Arth ~ N*P + S + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)

# Model selection - some models have singular fit because random effect doesn't explain any variation (term=0)

g2 <- glmer.nb( Arth ~ N*P + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g3 <- glmer.nb( Arth ~ N+P + S + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g4 <- glmer.nb( Arth ~ N+P + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g5 <- glmer.nb( Arth ~ N+ S + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g6 <- glmer.nb( Arth ~ P + S + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g7 <- glmer.nb( Arth ~ N + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g8 <- glmer.nb( Arth ~ P + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g9 <- glmer.nb( Arth ~ S + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g10 <- glmer.nb( Arth ~ 1 + offset(Lbud) + (1|Tree), family="poisson",data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)

L<-list()
L[[1]] <- g1
L[[2]] <- g2
L[[3]] <- g3
L[[4]] <- g4
L[[5]] <- g5
L[[6]] <- g6
L[[7]] <- g7
L[[8]] <- g8
L[[9]] <- g9
L[[10]] <- g10

names <- c("Birds x Planttype + Site","Birds x Planttype", "Birds + Planttype + Site", "Birds + Planttype", "Birds + Site", "Planttype + Site", "Birds", "Planttype", "Site", "Intercept")

aictab(L,names)

# Fit of best nb model is good

so <- simulateResiduals(fittedModel = g4, n = 10000, refit=F)
testUniformity(simulationOutput = so, plot=FALSE)
testDispersion(simulationOutput = so, plot=FALSE, alternative="greater")

# LRT hypothesis tests - some models have singular fit because random effect has no explanatory power

g0 <- glmer.nb( Arth ~ N*P + offset(Lbud) + (1|Tree),data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g1 <- glmer.nb( Arth ~ N+P + offset(Lbud) + (1|Tree),data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g2 <- glmer.nb( Arth ~ N + offset(Lbud) + (1|Tree),data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g3 <- glmer.nb( Arth ~ P + offset(Lbud) + (1|Tree),data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)
g4 <- glmer.nb( Arth ~ N + P + S + offset(Lbud) + (1|Tree),data=Dat, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),na.action=na.fail)

anova(g0,g1,test="Chisq")
anova(g1,g2,test="Chisq")
anova(g1,g3,test="Chisq") 
anova(g1,g4,test="Chisq") 

# Differences between group means 

Dat$Y <- Dat$ADens

d<-ddply(.data=Dat,.variables=c("N"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]

d<-ddply(.data=Dat,.variables=c("P"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]


# --------------------------------------------------------- #

# Herbivore analysis - final time point only 

Dat <- subset(Ca,subset=Wks==54)

Dat$PED<-Dat$PE/Dat$TotBud

# log transform bud number

Dat$Lbud <- log(Dat$TotBud)

# Global model 

g1<-glm(PE~N*P+S+offset(Lbud),family="poisson",data=Dat)

# Model selection

g1<-glm(PE~N*P+S+offset(Lbud),family="poisson",data=Dat)
g2<-glm(PE~N*P+offset(Lbud),family="poisson",data=Dat)
g3<-glm(PE~N+P+S+offset(Lbud),family="poisson",data=Dat)
g4<-glm(PE~N+P+offset(Lbud),family="poisson",data=Dat)
g5<-glm(PE~N+S+offset(Lbud),family="poisson",data=Dat)
g6<-glm(PE~P+S+offset(Lbud),family="poisson",data=Dat)
g7<-glm(PE~N+offset(Lbud),family="poisson",data=Dat)
g8<-glm(PE~P+offset(Lbud),family="poisson",data=Dat)
g9<-glm(PE~S+offset(Lbud),family="poisson",data=Dat)
g10<-glm(PE~1+offset(Lbud),family="poisson",data=Dat)

L<-list()
L[[1]] <- g1
L[[2]] <- g2
L[[3]] <- g3
L[[4]] <- g4
L[[5]] <- g5
L[[6]] <- g6
L[[7]] <- g7
L[[8]] <- g8
L[[9]] <- g9
L[[10]] <- g10

names <- c("Birds x Planttype + Site","Birds x Planttype", "Birds + Planttype + Site", "Birds + Planttype", "Birds + Site", "Planttype + Site", "Birds", "Planttype", "Site", "Intercept")

aictab(L,names)

# Overdispersion parameter of best model =1.42

sum( resid(g5, type="pearson")^2 ) / df.residual(g5)

# LRT hypothesis tests

g1<-glm(PE~N+offset(Lbud),family="poisson",data=Dat)
g2<-glm(PE~N+S +offset(Lbud),family="poisson",data=Dat)
g3<-glm(PE~N+P +offset(Lbud),family="poisson",data=Dat)
g4<-glm(PE~1+offset(Lbud),family="poisson",data=Dat)

anova(g1,g2,test="Chisq") 
anova(g1,g3,test="Chisq")
anova(g1,g4,test="Chisq")

# Differences between group means 

Dat$Y <- Dat$PED

d<-ddply(.data=Dat,.variables=c("N"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]

d<-ddply(.data=Dat,.variables=c("P"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]

d<-ddply(.data=Dat,.variables=c("S"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]



# ---------------------------------------------------------------------------- #
# Supplemental analyses: Analysis of individual growth measures
# ---------------------------------------------------------------------------- #


SupplementalGrowthAnalyses <- function(){}


# --------------------------------------------------------- #

# Analysis of individual growth measures: Change in terminal shoots
# gls used due to greater variance in floc trees

Dat <- Mdat

Dat$Y <- Dat$Change

g1 <- glm(Y~N*P+S,family="gaussian",data=Dat)

# boxplot(Dat$Y~Dat$P*Dat$S)
# plot(g1)
# plot(predict(g1),resid(g1, type="pearson"))
# qqnorm(resid(g1, type="pearson"))
# qqline(resid(g1, type="pearson"))

vf2id<-varIdent(form= ~1 | P)
g1<-gls(Y~N*P+S,data=Dat, weights=vf2id, method="ML")
g2<-gls(Y~N*P+S,data=Dat, method="ML")

anova(g1,g2)

# Global model (much more variation in Floc/Fringe plants)

g1<-gls(Y~N*P+S,data=Dat, weights=vf2id, method="ML")

# Model Selection 

g2 <- update(g1, ~ N*P)
g3 <- update(g1, ~ N+P+S)
g4 <- update(g1, ~ N+P)
g5 <- update(g1, ~ N+S)
g6 <- update(g1, ~ P+S)
g7 <- update(g1, ~ N)
g8 <- update(g1, ~ P)
g9 <- update(g1, ~ S)
g10 <- update(g1, ~ 1)

L<-list()
L[[1]] <- g1
L[[2]] <- g2
L[[3]] <- g3
L[[4]] <- g4
L[[5]] <- g5
L[[6]] <- g6
L[[7]] <- g7
L[[8]] <- g8
L[[9]] <- g9
L[[10]] <- g10

names <- c("Birds x Planttype + Site","Birds x Planttype", "Birds + Planttype + Site", "Birds + Planttype", "Birds + Site", "Planttype + Site", "Birds", "Planttype", "Site", "Intercept")

aictab(L,names)

# Best model

g1<-gls(Y~P,data=Dat, weights=vf2id, method="ML")

# LRT hypothesis tests

g2<-gls(Y~1,data=Dat, weights=vf2id, method="ML") # P: <0.001
g3<-gls(Y~P+N,data=Dat, weights=vf2id, method="ML")
g4<-gls(Y~P+S,data=Dat, weights=vf2id, method="ML")

anova(g1,g2) 
anova(g1,g3)
anova(g1,g4)

# Differences between group means 

d<-ddply(.data=Dat,.variables=c("P"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]


# --------------------------------------------------------- #

# Analysis of individual growth measures: Change in shoot length
# gls used due to greater variance in floc trees

Dat <- Mdat

Dat$Y <- Dat$Dtl

g1 <- glm(Y~N*P+S,family="gaussian",data=Dat)

# boxplot(Dat$Y~Dat$P*Dat$S)
# plot(g1)
# plot(predict(g1),resid(g1, type="pearson"))
# qqnorm(resid(g1, type="pearson"))
# qqline(resid(g1, type="pearson"))

vf2id<-varIdent(form= ~1 | P)
g1<-gls(Y~N*P+S,data=Dat, weights=vf2id, method="ML")
g2<-gls(Y~N*P+S,data=Dat, method="ML")

anova(g1,g2)

# Global model (much more variation in Floc/Fringe plants)

g1<-gls(Y~N*P+S,data=Dat, weights=vf2id, method="ML")

#plot(g1)

# Model Selection 

g2 <- update(g1, ~ N*P)
g3 <- update(g1, ~ N+P+S)
g4 <- update(g1, ~ N+P)
g5 <- update(g1, ~ N+S)
g6 <- update(g1, ~ P+S)
g7 <- update(g1, ~ N)
g8 <- update(g1, ~ P)
g9 <- update(g1, ~ S)
g10 <- update(g1, ~ 1)

L<-list()
L[[1]] <- g1
L[[2]] <- g2
L[[3]] <- g3
L[[4]] <- g4
L[[5]] <- g5
L[[6]] <- g6
L[[7]] <- g7
L[[8]] <- g8
L[[9]] <- g9
L[[10]] <- g10

names <- c("Birds x Planttype + Site","Birds x Planttype", "Birds + Planttype + Site", "Birds + Planttype", "Birds + Site", "Planttype + Site", "Birds", "Planttype", "Site", "Intercept")

aictab(L,names)

# Best model

g1<-gls(Y~P,data=Dat, weights=vf2id, method="ML")

# LRT hypothesis tests

g2<-gls(Y~1,data=Dat, weights=vf2id, method="ML") # P: <0.001
g3<-gls(Y~P+N,data=Dat, weights=vf2id, method="ML")
g4<-gls(Y~P+S,data=Dat, weights=vf2id, method="ML")

anova(g1,g2) 
anova(g1,g3)
anova(g1,g4)

# Differences between group means 

d<-ddply(.data=Dat,.variables=c("P"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]


# --------------------------------------------------------- #

# Analysis of individual growth measures: Change in leaf number
# gls used due to greater variance in floc trees

Dat <- Mdat

Dat$Y <- Dat$DtLv

g1 <- glm(Y~N*P+S,family="gaussian",data=Dat)

# boxplot(Dat$Y~Dat$P*Dat$S)
# plot(g1)
# plot(predict(g1),resid(g1, type="pearson"))
# qqnorm(resid(g1, type="pearson"))
# qqline(resid(g1, type="pearson"))

vf2id<-varIdent(form= ~1 | P)
g1<-gls(Y~N*P+S,data=Dat, weights=vf2id, method="ML")
g2<-gls(Y~N*P+S,data=Dat, method="ML")

anova(g1,g2)

# Global model (much more variation in Floc/Fringe plants)

g1<-gls(Y~N*P+S,data=Dat, weights=vf2id, method="ML")

#plot(g1)

# Model Selection 

g2 <- update(g1, ~ N*P)
g3 <- update(g1, ~ N+P+S)
g4 <- update(g1, ~ N+P)
g5 <- update(g1, ~ N+S)
g6 <- update(g1, ~ P+S)
g7 <- update(g1, ~ N)
g8 <- update(g1, ~ P)
g9 <- update(g1, ~ S)
g10 <- update(g1, ~ 1)

L<-list()
L[[1]] <- g1
L[[2]] <- g2
L[[3]] <- g3
L[[4]] <- g4
L[[5]] <- g5
L[[6]] <- g6
L[[7]] <- g7
L[[8]] <- g8
L[[9]] <- g9
L[[10]] <- g10

names <- c("Birds x Planttype + Site","Birds x Planttype", "Birds + Planttype + Site", "Birds + Planttype", "Birds + Site", "Planttype + Site", "Birds", "Planttype", "Site", "Intercept")

aictab(L,names)

# Best model

g1<-gls(Y~P,data=Dat, weights=vf2id, method="ML")

# LRT hypothesis tests

g2<-gls(Y~1,data=Dat, weights=vf2id, method="ML") # P: <0.001
g3<-gls(Y~P+N,data=Dat, weights=vf2id, method="ML")
g4<-gls(Y~P+S,data=Dat, weights=vf2id, method="ML")

anova(g1,g2) 
anova(g1,g3)
anova(g1,g4)

# Differences between group means 

d<-ddply(.data=Dat,.variables=c("P"),.fun="summarise",
m=mean(Y,na.rm=TRUE),
se=sd(Y,na.rm=TRUE)/sqrt(length(Y)))
d
(d$m[1]-d$m[2])/d$m[2] 
(d$m[2]-d$m[1])/d$m[1]
d$m[1]/d$m[2]
d$m[2]/d$m[1]

# END
