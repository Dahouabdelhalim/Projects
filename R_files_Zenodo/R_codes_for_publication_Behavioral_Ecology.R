library(tnet)
library(bipartite)
library(lmtest)
library(car)
library(rcompanion)

#CAMP1
#out of camp proximity networks camp 1
CAMP1gpsOUT_OF_CAMP=read.csv("~CAMP1gpsOUT_OF_CAMP.csv",header=TRUE,row.names=1,check.names=FALSE) 
CAMP1gpsOUT_OF_CAMP<-as.matrix(CAMP1gpsOUT_OF_CAMP)
#in of camp proximity networks camp 1
CAMP1gpsIN_CAMP=read.csv("~CAMP1gpsIN_CAMP.csv",header=TRUE,row.names=1,check.names=FALSE) 
CAMP1gpsIN_CAMP<-as.matrix(CAMP1gpsIN_CAMP)
#friendship  nomination camp 1
CAMP1friendship=read.csv("~CAMP1friendship.csv",header=TRUE,row.names=1,check.names=FALSE) 
CAMP1friendship<-as.matrix(CAMP1friendship)
nrow(CAMP1friendship)
#kinship matrix camp 1
CAMP1kinship<-read.csv("~CAMP1kinship.csv",header=TRUE,row.names=1,check.names=FALSE) 
CAMP1kinship<-as.matrix(CAMP1kinship)
#hunting nomination matrix camp 1
CAMP1hunter<-read.csv("~CAMP1hunter.csv",header=TRUE,row.names=1,check.names=FALSE) 
CAMP1hunter<-as.matrix(CAMP1hunter)

#CAMP2
#out os of camp proximity networks camp 2
CAMP2gpsOUT_OF_CAMP=read.csv("~CAMP2gpsOUT_OF_CAMP.csv",header=TRUE,row.names=1,check.names=FALSE) 
CAMP2gpsOUT_OF_CAMP<-as.matrix(CAMP2gpsOUT_OF_CAMP)
#in of camp proximity networks camp 2
CAMP2gpsIN_CAMP=read.csv("~CAMP2gpsIN_CAMP.csv",header=TRUE,row.names=1,check.names=FALSE) 
CAMP2gpsIN_CAMP<-as.matrix(CAMP2gpsIN_CAMP)
#firenship network camp 2
CAMP2friendship=read.csv("CAMP2friendship~csv",header=TRUE,row.names=1,check.names=FALSE) 
CAMP2friendship<-as.matrix(CAMP2friendship)
#kinship matrix camp 2
CAMP2kinship<-read.csv("~CAMP2kinship.csv",header=TRUE,row.names=1,check.names=FALSE)
CAMP2kinship<-as.matrix(CAMP2kinship)
#hunting nomination matrix camp 2
CAMP2hunter<-read.csv("~CAMP2hunter.csv",header=TRUE,row.names=1,check.names=FALSE)
CAMP2hunter<-as.matrix(CAMP2hunter)

#CAMP3
#out os of camp proximity networks camp 3
CAMP3gpsOUT_OF_CAMP=read.csv("~CAMP3gpsOUT_OF_CAMP.csv",header=TRUE,row.names=1,check.names=FALSE) 
CAMP3gpsOUT_OF_CAMP<-as.matrix(CAMP3gpsOUT_OF_CAMP)
#in of camp proximity networks camp 3
CAMP3gpsIN_CAMP=read.csv("~CAMP3gpsIN_CAMP.csv",header=TRUE,row.names=1,check.names=FALSE) 
CAMP3gpsIN_CAMP<-as.matrix(CAMP3gpsIN_CAMP)
#firenship network camp 3
CAMP3friendship=read.csv("~CAMP3friendship.csv",header=TRUE,row.names=1,check.names=FALSE) 
CAMP3friendship<-as.matrix(CAMP3friendship)
#kinship matrix camp 3
CAMP2kinship=read.csv("~CAMP2kinship.csv",header=TRUE,row.names=1,check.names=FALSE) 
CAMP2kinship<-as.matrix(CAMP2kinship)
#hunting nomination matrix camp 3
CAMP3hunter=read.csv("~CAMP3hunter.csv",header=TRUE,row.names=1,check.names=FALSE) 
CAMP3hunter<-as.matrix(CAMP3hunter)


#FRIENDSHIP POPULARITY

M<-as.matrix(CAMP1friendship)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
friend_strength1<-v


M<-as.matrix(CAMP2friendship)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
friend_strength2<-v



M<-as.matrix(CAMP3friendship)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
friend_strength3<-v


friend_strength<-c(friend_strength1,friend_strength2,friend_strength3)

#FORAGING REPUTATION

M<-as.matrix(CAMP1hunter)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
foraging_strength1<-v


M<-as.matrix(CAMP2hunter)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
foraging_strength2<-v



M<-as.matrix(CAMP3hunter)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
foraging_strength3<-v


foraging_strength<-c(foraging_strength1,foraging_strength2,foraging_strength3)


#gps_strengthOUT_OF_CAMP
M<-as.matrix(CAMP1gpsOUT_OF_CAMP)
M<-M/10000#returning the original proximity values, see README 
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
gps_strength1<-v
gps_strength1

M<-as.matrix(CAMP2gpsOUT_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
gps_strength2<-v




M<-as.matrix(CAMP3gpsOUT_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
gps_strength3<-v

gps_strength3
gps_strengthOUT<-c(gps_strength1,gps_strength2,gps_strength3)
length(gps_strengthOUT)
gps_strengthOUT

#gps_strengthIN_OF_CAMP
M<-as.matrix(CAMP1gpsIN_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
gps_strength1<-v
gps_strength1

M<-as.matrix(CAMP2gpsIN_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
gps_strength2<-v




M<-as.matrix(CAMP3gpsIN_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
gps_strength3<-v

gps_strength3
gps_strengthIN<-c(gps_strength1,gps_strength2,gps_strength3)
length(gps_strengthIN)
gps_strengthIN

#gps degreeOUT_OF_CAMP

M<-as.matrix(CAMP1gpsOUT_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
gps_strength1<-v


M<-as.matrix(CAMP2gpsOUT_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
v
gps_strength2<-v

M<-as.matrix(CAMP3gpsOUT_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
gps_strength3<-v
gps_strength3
gps_degreeOUT<-c(gps_strength1,gps_strength2,gps_strength3)
length(gps_degreeOUT)


#gps betweenness OUT OF CAMP 

M<-as.matrix(CAMP1gpsOUT_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-betweenness_w(b,alpha=0.5)  
v<-v[,-1]
gps_betweenness1<-v

M<-as.matrix(CAMP2gpsOUT_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-betweenness_w(b,alpha=0.5) 
v
v<-v[,-1]
gps_betweenness2<-v
gps_betweenness2



M<-as.matrix(CAMP3gpsOUT_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-betweenness_w(b,alpha=0.5)  
v<-v[,-1]
gps_betweenness3<-v
gps_betweenness3

gps_betweennessOUT<-c(gps_betweenness1,gps_betweenness2,gps_betweenness3)


#KINDEGREE
M<-as.matrix(CAMP1kinship)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
v<-c(v,0)
kin_strength1<-v
length(kin_strength1)

M<-as.matrix(CAMP2kinship)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
kin_strength2<-v

M<-as.matrix(CAMP3kinship)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
kin_strength3<-v


kin_degree<-c(kin_strength1,kin_strength2,kin_strength3)

#gps degree IN CAMP
M<-as.matrix(CAMP1gpsIN_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
gps_strength1<-v

M<-as.matrix(CAMP2gpsIN_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
v
gps_strength2<-v


M<-as.matrix(CAMP3gpsIN_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
gps_strength3<-v

gps_degreeIN<-c(gps_strength1,gps_strength2,gps_strength3)

#gps betweenness IN CAMP

M<-as.matrix(CAMP1gpsIN_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-betweenness_w(b,alpha=0.5)  
v<-v[,-1]
sd(v)
gps_betweenness1<-v


M<-as.matrix(CAMP2gpsIN_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-betweenness_w(b,alpha=0.5) 
v<-v[,-1]
sd(v,na.rm=T)
gps_betweenness2<-v

M<-as.matrix(CAMP3gpsIN_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-betweenness_w(b,alpha=0.5)  
v<-v[,-1]
sd(v,na.rm=T)
gps_betweenness3<-v


gps_betweennessIN<-c(gps_betweenness1,gps_betweenness2,gps_betweenness3)

# GPS CLUSTERING IM CAMP
M<-as.matrix(CAMP2gpsIN_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
#b<-symmetrise_w(b)
v<-clustering_local_w(b, measure=c("am"))
v<-v[,-1]
gps_clustering2<-v

M<-as.matrix(CAMP1gpsIN_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
#cvc<-quantile(b, probs = 0.75)
#cvc<-unname(cvc)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
#b<-symmetrise_w(b)
v<-clustering_local_w(b, measure=c("am"))
v<-v[,-1]
gps_clustering1<-v
#gps_betweenness1<-scale(gps_betweenness1)

M<-as.matrix(CAMP3gpsIN_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
#cvc<-quantile(b, probs = 0.75)
#cvc<-unname(cvc)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
#b<-symmetrise_w(b)
v<-clustering_local_w(b, measure=c("am")) 
v<-v[,-1]
v
gps_clustering3<-v

gps_clusteringIN<-c(gps_clustering1,gps_clustering2,gps_clustering3)

# GPS CLUSTERING OUT OF CAMP

M<-as.matrix(CAMP1gpsOUT_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
#b<-symmetrise_w(b)
v<-clustering_local_w(b, measure=c("am"))
v<-v[,-1]
v
sd(v,na.rm=T)
gps_clustering1<-v

M<-as.matrix(CAMP2gpsOUT_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
#b<-symmetrise_w(b)
v<-clustering_local_w(b, measure=c("am"))
v<-v[,-1]
v
sd(v,na.rm=T)
gps_clustering2<-v


M<-as.matrix(CAMP3gpsOUT_OF_CAMP)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
#b<-symmetrise_w(b)
v<-clustering_local_w(b, measure=c("am")) 
v<-v[,-1]
v
sd(v,na.rm=T)
gps_clustering3<-v

gps_clusteringOUT<-c(gps_clustering1,gps_clustering2,gps_clustering3)
gps_clusteringOUT


#library(simr)
#library(lme4)
all10<-read.csv("~FILE1.csv") #IMPORTING FILE WITH CAMP MEMBERSHIP, SEX AND AGE DATA
nrow(all10)
all10<-all10[1:81,]


all10<-all10[,c("NAME","sex","age2","CAMP")]#
#COMPILING NETWORK METRICS AND FILE 1 DATA
all10<-cbind(all10,gps_betweennessIN)
all10<-cbind(all10,gps_betweennessOUT)
all10<-cbind(all10,friend_strength)
all10<-cbind(all10,gps_strengthOUT)
all10<-cbind(all10,gps_strengthIN)
all10<-cbind(all10,foraging_strength)
all10<-cbind(all10,gps_degreeIN)
all10<-cbind(all10,gps_degreeOUT)
all10<-cbind(all10,gps_clusteringIN)
all10<-cbind(all10,gps_clusteringOUT)
all10<-cbind(all10,kin_degree)



#nrow(all10)
#all10$age2<-scale(all10$age2)
str(all10$CAMP)

all10<-all10[all10$sex=="1",]#man data extracted

mean(all10$age)
median(all10$age)
nrow(all10)

all10$age

# SCALING NETWORK METRICS DATA DEPENDING ON CAMP MEMBERSHIP

#betweenness in camp
a<-all10$gps_betweennessIN[1:14]#first camp
all10$NAME[1:14]
#mean(a)
#sd(a)
#range(a)
a<-scale(as.numeric(a))
#length(a) 
a<-scale(as.numeric(a))
#length(a)
b<-all10$gps_betweennessIN[15:31]#second camp
all10$NAME[15:31]
#mean(b)
#sd(b)
#range(b)
b<-scale(as.numeric(b))
#length(b)
c<-all10$gps_betweennessIN[32:35]
all10$NAME[32:35]
#mean(c)
#sd(c)
#range(c)
c<-scale(as.numeric(c))
#length(c)
all10$gps_betweennessIN<-c(a,b,c)

#betweenness out
a<-all10$gps_betweennessOUT[1:14]#first camp
#all10$NAME[1:14]
#mean(a)
#sd(a)
#range(a)
a<-scale(as.numeric(a))
#length(a)
b<-all10$gps_betweennessOUT[15:31]#second camp
#all10$NAME[15:31]
#mean(b)
#sd(b)
#range(b)
b<-scale(as.numeric(b))
#length(b)
c<-all10$gps_betweennessOUT[32:35]
#all10$NAME[32:35]
#mean(c)
#sd(c)
#range(c)
c<-scale(as.numeric(c))
#length(c)
all10$gps_betweennessOUT<-c(a,b,c)



#friend strength
a<-all10$friend_strength[1:14]#first camp
#all10$NAME[1:14]
a<-scale(as.numeric(a))
#length(a)
b<-all10$friend_strength[15:31]#second camp
#all10$NAME[15:31]
b<-scale(as.numeric(b))
#length(b)
c<-all10$friend_strength[32:35]
#all10$NAME[32:35]
c<-scale(as.numeric(c))
#length(c)
all10$friend_strength<-c(a,b,c)

#foraging strength
a<-all10$foraging_strength[1:14]#first camp
#all10$NAME[1:14]
a<-scale(as.numeric(a))
#length(a)
b<-all10$foraging_strength[15:31]#second camp
#all10$NAME[15:31]
b<-scale(as.numeric(b))
#length(b)
c<-all10$foraging_strength[32:35]
#all10$NAME[32:35]
c<-scale(as.numeric(c))
#length(c)
all10$foraging_strength<-c(a,b,c)

#GPS_DEGREE_IN

a<-all10$gps_degreeIN[1:14]#first camp
#all10$NAME[1:14]
#mean(a)
#sd(a)
#range(a)
a<-scale(as.numeric(a))
#length(a)
b<-all10$gps_degreeIN[15:31]#second camp
#all10$NAME[15:31]
#mean(b)
#sd(b)
#range(b)
b<-scale(as.numeric(b))
#length(b)
c<-all10$gps_degreeIN[32:35]
#all10$NAME[32:35]
#mean(c)
#sd(c)
#range(c)
c<-scale(as.numeric(c))
#length(c)
all10$gps_degreeIN<-c(a,b,c)

# gps_degree_out

a<-all10$gps_degreeOUT[1:14]#first camp
#all10$NAME[1:14]
#mean(a)
#sd(a)
#range(a)
#sd(a)
a<-scale(as.numeric(a))
#length(a)
#sd(a)
b<-all10$gps_degreeOUT[15:31]#second camp
all10$NAME[15:31]
#mean(b)
#sd(b)
#range(b)
b<-scale(as.numeric(b))
#length(b)
c<-all10$gps_degreeOUT[32:35]
#mean(c)
#sd(c)
#range(c)
c<-scale(as.numeric(c))
#length(c)
all10$gps_degreeOUT<-c(a,b,c)

# gps_strength_out

a<-all10$gps_strengthOUT[1:14]#first camp
#all10$NAME[1:14]
#mean(a)
#sd(a)
a<-scale(as.numeric(a))
#length(a)
b<-all10$gps_strengthOUT[15:31]#second camp
#all10$NAME[15:31]
b<-scale(as.numeric(b))
#length(b)
c<-all10$gps_strengthOUT[32:35]
c<-scale(as.numeric(c))
#length(c)
all10$gps_strengthOUT<-c(a,b,c)

# gps strength in

a<-all10$gps_strengthIN[1:14]#first camp
#all10$NAME[1:14]
a<-scale(as.numeric(a))
#length(a)
b<-all10$gps_strengthIN[15:31]#second camp
#all10$NAME[15:31]
b<-scale(as.numeric(b))
#length(b)
c<-all10$gps_strengthIN[32:35]
c<-scale(as.numeric(c))
#length(c)
all10$gps_strengthIN<-c(a,b,c)



#kin degree
a<-all10$kin_degree[1:14]#first camp
#all10$NAME[1:14]
a<-scale(as.numeric(a))
#length(a)
b<-all10$kin_degree[15:31]#second camp
#all10$NAME[15:31]
b<-scale(as.numeric(b))
#length(b)
c<-all10$kin_degree[32:35]
#all10$NAME[32:35]
c<-scale(as.numeric(c))
#length(c)
all10$kin_degree<-c(a,b,c)

#CLUSTERING in
#all10$gps_clusteringIN[is.nan(all10$gps_clusteringIN)] <- 0 
a<-all10$gps_clusteringIN[1:14]#first camp
#all10$NAME[1:14]
#mean(a,na.rm=T)
#sd(a,na.rm=T)
#range(a)
a<-scale(as.numeric(a))
#length(a)
b<-all10$gps_clusteringIN[15:31]#second camp
#all10$NAME[15:31]
#mean(b)
#sd(b)
#range(b)
b<-scale(as.numeric(b))
#length(b)
c<-all10$gps_clusteringIN[32:35]
#all10$NAME[32:35]
#length(c)
#mean(c,na.rm=T)
#sd(c,na.rm=T)
#range(c,na.rm=T)
c<-scale(as.numeric(c))
#all10$gps_clusteringIN<-c(a,b,c)

#CLUSTERING OUT
#all10$gps_clusteringOUT[is.nan(all10$gps_clusteringOUT)] <- 0 
a<-all10$gps_clusteringOUT[1:14]#first camp
#all10$NAME[1:14]
#mean(a,na.rm=T)
#sd(a,na.rm=T)
#range(a,na.rm=T)
a<-scale(as.numeric(a))
#length(a)
b<-all10$gps_clusteringOUT[15:31]#second camp
#all10$NAME[15:31]
#mean(b,na.rm=T)
#sd(b,na.rm=T)
#range(b,na.rm=T)
b<-scale(as.numeric(b))
#length(b)
c<-all10$gps_clusteringOUT[32:35]
#all10$NAME[32:35]
#length(c)
#mean(c,na.rm=T)
#sd(c)
#c<-scale(as.numeric(c))
all10$gps_clusteringOUT<-c(a,b,c)

cor.test(all10$friend_strength,all10$foraging_strength,method = "spearman")
cor.test(all10$gps_degreeOUT,all10$gps_betweennessOUT,method = "spearman")
cor.test(all10$gps_degreeOUT,all10$gps_clusteringOUT,method = "spearman",na.rm=T)
cor.test(all10$gps_clusteringOUT,all10$gps_betweennessOUT,method = "spearman",na.rm=T)

cor.test(all10$gps_degreeIN,all10$gps_betweennessIN,method = "spearman")
cor.test(all10$gps_degreeIN,all10$gps_clusteringIN,method = "spearman",na.rm=T)
cor.test(all10$gps_clusteringIN,all10$gps_betweennessIN,method = "spearman",na.rm=T)
cor.test(all10$foraging_strength,all10$age2,method = "spearman",na.rm=T)
cor.test(all10$friend_strength,all10$age2,method = "spearman",na.rm=T)




#BETWEENNESS OUT OF CAMP

library(MCMCpack)
nrow(all10)

all11<-all10[-20,]

range(all11$gps_betweennessOUT)


all11$gps_betweennessOUT<-all11$gps_betweennessOUT+1.161896
range(all11$gps_betweennessOUT)


#all11<-all10
library(rcompanion)
all11$gps_betweennessOUT =#DATA TRANSFORMATION
  transformTukey(all11$gps_betweennessOUT,
                 plotit=FALSE)


shapiro.test(all11$gps_betweennessOUT)

all11$BETWEENNESS<-all11$gps_betweennessOUT
all11$POPULARITY<-all11$friend_strength
all11$HUNTING_REPUTATION<-all11$foraging_strength
all11$AGE<-all11$age2
all11$NUMBER_OF_KIN<-all11$kin_degree

all11<-cbind(all11,all11$BETWEENNESS)
all11<-cbind(all11,all11$POPULARITY)
all11<-cbind(all11,all11$HUNTING_REPUTATION)
all11<-cbind(all11,all11$AGE)
all11<-cbind(all11,all11$NUMBER_OF_KIN)

log_fit2 <-MCMCregress(BETWEENNESS~AGE+POPULARITY+HUNTING_REPUTATION+NUMBER_OF_KIN, seed=123,thin=10,burnin = 10000, mcmc=100000, data=all11)
summary(log_fit2)
summary(log_fit2 , probs=c(.025, .975), digits=3)
dev.off()
plot(log_fit2)

#BETWEENNESS IN CAMP
range(all10$gps_betweennessIN)

all10$gps_betweennessIN<-all10$gps_betweennessIN+1.009963
all10$gps_betweennessIN =#DATA TRANSFORMATION
  transformTukey(all10$gps_betweennessIN,
                 plotit=FALSE)


shapiro.test(all10$gps_betweennessIN)

kruskal.test(gps_betweennessIN ~ CAMP, data = all10) 
all10$BETWEENNESS<-all10$gps_betweennessIN
all10$POPULARITY<-all10$friend_strength
all10$HUNTING_REPUTATION<-all10$foraging_strength
all10$AGE<-all10$age2
all10$NUMBER_OF_KIN<-all10$kin_degree

all10<-cbind(all10,all10$BETWEENNESS)
all10<-cbind(all10,all10$POPULARITY)
all10<-cbind(all10,all10$HUNTING_REPUTATION)
all10<-cbind(all10,all10$AGE)
all10<-cbind(all10,all10$NUMBER_OF_KIN)

log_fit2 <-MCMCregress(BETWEENNESS~AGE+POPULARITY+HUNTING_REPUTATION+NUMBER_OF_KIN, seed=123,thin=10,burnin = 10000, mcmc=100000, data=all10)
summary(log_fit2)
summary(log_fit2 , probs=c(.025, .975), digits=3)
dev.off()
plot(log_fit2)


#CLUSTERING OUR OF CAMP
#library(MCMCpack)
#library(rcompanion)

nrow(all10)
all14<-all10
all14<-all14[-5,]#removing individuals with NA values
all14<-all14[-7,]
all14<-all14[-14,]
all14<-all14[-16,]
all14<-all14[-26,]
all14<-all14[-26,]
all14$gps_clusteringOUT
#all10$gps_clusteringOUT<-(all10$gps_clusteringOUT)
range(all14$gps_clusteringOUT)
#all10$gps_clusteringOUT[is.nan(all10$gps_clusteringOUT)] <- 0 
#all14<-all14[1:26,]
#all14$gps_clusteringOUT<-all14$gps_clusteringOUT+1.439349#1.915734
#cc<-log(all14$gps_clusteringOUT)
#all14$gps_clusteringOUT =#DATA TRANSFORMATION
#  transformTukey(all14$gps_clusteringOUT)
#all10$gps_degreeOUT<-gps_degreeOUT

shapiro.test(all14$gps_clusteringOUT)

kruskal.test(gps_clusteringOUT ~ CAMP, data = all14) 
all14$CLUSTERING<-all14$gps_clusteringOUT
all14$POPULARITY<-all14$friend_strength
all14$HUNTING_REPUTATION<-all14$foraging_strength
all14$AGE<-all14$age2
all14$NUMBER_OF_KIN<-all14$kin_degree

all14<-cbind(all14,all14$CLUSTERING)
all14<-cbind(all14,all14$POPULARITY)
all14<-cbind(all14,all14$HUNTING_REPUTATION)
all14<-cbind(all14,all14$AGE)
all14<-cbind(all14,all14$NUMBER_OF_KIN)

log_fit2 <-MCMCregress(CLUSTERING~AGE+POPULARITY+HUNTING_REPUTATION+NUMBER_OF_KIN,sigma.mu = NA,sigma.var = NA, seed=123,thin=10,burnin = 10000, mcmc=100000, data=all14)
summary(log_fit2)
summary(log_fit2 , probs=c(.025, .975), digits=3)
dev.off()
plot(log_fit2)


#CLUSTERING IN CAMP
nrow(all10)


all17<-all10[-35,]#removing individuals with NA values

range(all17$gps_clusteringIN)
#all17$gps_clusteringIN<-all17$gps_clusteringIN+1.74#1.915734
#all17$gps_clusteringIN<-all17$gps_clusteringIN

#all17$gps_clusteringIN =#DATA TRANSFORMATION
  #  transformTukey(all17$gps_clusteringIN,
                 #                 plotit=FALSE)
#all10$gps_degreeOUT<-gps_degreeOUT

shapiro.test(all17$gps_clusteringIN)

kruskal.test(gps_clusteringIN ~ CAMP, data = all17) 

all17$CLUSTERING<-all17$gps_clusteringIN
all17$POPULARITY<-all17$friend_strength
all17$HUNTING_REPUTATION<-all17$foraging_strength
all17$AGE<-all17$age2
all17$NUMBER_OF_KIN<-all17$kin_degree

all17<-cbind(all17,all17$CLUSTERING)
all17<-cbind(all17,all17$POPULARITY)
all17<-cbind(all17,all17$HUNTING_REPUTATION)
all17<-cbind(all17,all17$AGE)
all17<-cbind(all17,all17$NUMBER_OF_KIN)

log_fit2 <-MCMCregress(CLUSTERING~AGE+POPULARITY+HUNTING_REPUTATION+NUMBER_OF_KIN, seed=123,thin=10,burnin = 10000, mcmc=100000, data=all17)
summary(log_fit2)
summary(log_fit2 , probs=c(.025, .975), digits=3)
dev.off()
plot(log_fit2)



#degree out OF CAMP
library(rcompanion)
nrow(all10)
range(all10$gps_degreeOUT)
#all10$gps_degreeOUT<-all10$gps_degreeOUT+1.860016#1.915734
#all10$gps_degreeOUT =#DATA TRANSFORMATION
#  transformTukey(all10$gps_degreeOUT,
#                 plotit=FALSE)

shapiro.test(all10$gps_degreeOUT)
all10$DEGREE<-all10$gps_degreeOUT
all10$POPULARITY<-all10$friend_strength
all10$HUNTING_REPUTATION<-all10$foraging_strength
all10$AGE<-all10$age2
all10$NUMBER_OF_KIN<-all10$kin_degree

all10<-cbind(all10,all10$DEGREE)
all10<-cbind(all10,all10$POPULARITY)
all10<-cbind(all10,all10$HUNTING_REPUTATION)
all10<-cbind(all10,all10$AGE)
all10<-cbind(all10,all10$NUMBER_OF_KIN)

log_fit2 <-MCMCregress(DEGREE~AGE+POPULARITY+HUNTING_REPUTATION+NUMBER_OF_KIN, seed=123,thin=10,burnin = 10000, mcmc=100000, data=all10)
summary(log_fit2)
summary(log_fit2 , probs=c(.025, .975), digits=3)
dev.off()
plot(log_fit2)


#degreeIN CAMP
range(all10$gps_degreeIN)

#all10$gps_degreeIN<-all10$gps_degreeIN+1.747669

#all10$gps_degreeIN =#DATA TRANSFORMATION
# transformTukey(all10$gps_degreeIN,
#                 plotit=FALSE)
#all10$gps_degreeOUT<-gps_degreeOUT

shapiro.test(all10$gps_degreeIN)

kruskal.test(gps_degreeIN ~ CAMP, data = all10) 

all10$DEGREE<-all10$gps_degreeIN
all10$POPULARITY<-all10$friend_strength
all10$HUNTING_REPUTATION<-all10$foraging_strength
all10$AGE<-all10$age2
all10$NUMBER_OF_KIN<-all10$kin_degree

all10<-cbind(all10,all10$DEGREE)
all10<-cbind(all10,all10$POPULARITY)
all10<-cbind(all10,all10$HUNTING_REPUTATION)
all10<-cbind(all10,all10$AGE)
all10<-cbind(all10,all10$NUMBER_OF_KIN)

log_fit2 <-MCMCregress(DEGREE~AGE+POPULARITY+HUNTING_REPUTATION+NUMBER_OF_KIN, seed=123,thin=10,burnin = 10000, mcmc=100000, data=all10)
summary(log_fit2)
summary(log_fit2 , probs=c(.025, .975), digits=3)
dev.off()
plot(log_fit2)

#STRENGTH IN CAMP
nrow(all10)
range(all10$gps_strengthIN)
#all10<-all10[1:31,]
#all10$gps_strengthIN<-all10$gps_strengthIN+3#1.915734
#all10$gps_strengthIN =#DATA TRANSFORMATION
# transformTukey(all10$gps_strengthIN,
#                plotit=FALSE)

shapiro.test(all10$gps_strengthIN)

kruskal.test(gps_strengthIN ~ CAMP, data = all10) 

all10$STRENGTH<-all10$gps_strengthIN
all10$POPULARITY<-all10$friend_strength
all10$HUNTING_REPUTATION<-all10$foraging_strength
all10$AGE<-all10$age2
all10$NUMBER_OF_KIN<-all10$kin_degree

all10<-cbind(all10,all10$STRENGTH)
all10<-cbind(all10,all10$POPULARITY)
all10<-cbind(all10,all10$HUNTING_REPUTATION)
all10<-cbind(all10,all10$AGE)
all10<-cbind(all10,all10$NUMBER_OF_KIN)

log_fit2 <-MCMCregress(STRENGTH~AGE+POPULARITY+HUNTING_REPUTATION+NUMBER_OF_KIN, seed=123,thin=10,burnin = 10000, mcmc=100000, data=all10)
summary(log_fit2)
summary(log_fit2 , probs=c(.025, .975), digits=3)
dev.off()
plot(log_fit2)
log_fit1 <-MCMCregress(gps_degreeIN~age2+friend_strength+foraging_strength+kin_degree, seed=123,thin=10,burnin = 10000, mcmc=100000, data=all11)
summary(log_fit1)
plot(log_fit1)




#STRENGTH OUT OF CAMP

nrow(all10)
range(all10$gps_strengthOUT)

#all10$gps_strengthOUT<-all10$gps_strengthOUT+1.88
#all10$gps_strengthOUT =#DATA TRANSFORMATION
#  transformTukey(all10$gps_strengthOUT,
#                plotit=FALSE)

shapiro.test(all10$gps_strengthOUT)

all10$STRENGTH<-all10$gps_strengthOUT
all10$POPULARITY<-all10$friend_strength
all10$HUNTING_REPUTATION<-all10$foraging_strength
all10$AGE<-all10$age2
all10$NUMBER_OF_KIN<-all10$kin_degree

all10<-cbind(all10,all10$STRENGTH)
all10<-cbind(all10,all10$POPULARITY)
all10<-cbind(all10,all10$HUNTING_REPUTATION)
all10<-cbind(all10,all10$AGE)
all10<-cbind(all10,all10$NUMBER_OF_KIN)

log_fit2 <-MCMCregress(STRENGTH~AGE+POPULARITY+HUNTING_REPUTATION+NUMBER_OF_KIN, seed=123,thin=10,burnin = 10000, mcmc=100000, data=all10)
summary(log_fit2)
summary(log_fit2 , probs=c(.025, .975), digits=3)
dev.off()
plot(log_fit2)


#permutation p-values (example)

all20<-all10
log_fit1 <-lm(gps_clusteringIN~age2+friend_strength+foraging_strength+kin_degree, data=all20)
summary(log_fit1)
OBS.VALUE<-summary(log_fit1)[["coefficients"]][, "Pr(>|t|)"] #[, "t value"]#[, "Estimate"]# [, "t value"]#
OBS.VALUE<-OBS.VALUE[5]
OBS.VALUE

set.seed(1234)
cluster_vect2 <- ""    
for (i in 1:1000)    
{
  cv<-all20$gps_clusteringIN[1:14]
  
  cv<-sample(cv, 14, replace = FALSE, prob = NULL)
  
  cv1<-all20$gps_clusteringIN[15:31]
  
  cv1<-sample(cv1, 17, replace = FALSE, prob = NULL)
  
  cv2<-all20$gps_clusteringIN[32:34]
  
  cv2<-sample(cv2, 3, replace = FALSE, prob = NULL)
  
  cv3<-c(cv,cv1,cv2)
  all20$gps_clusteringIN<-c(cv,cv1,cv2)
  
  log_fit1 <-lm(gps_clusteringIN~age2+friend_strength+foraging_strength+kin_degree, data=all20)
  
  VV<-summary(log_fit1)[["coefficients"]][, "Pr(>|t|)"] #[, "t value"]#[, "Estimate"]#[, "Pr(>|t|)"]
  VV<-VV[5]
  VV
  tmpval2<-VV    
  cluster_vect2 <- c(cluster_vect2, tmpval2)
  
}
gggg<-as.numeric(cluster_vect2)
gggg<-gggg[-1]
#gggg
sum(OBS.VALUE> gggg)/1000

OBS.VALUE




#GRAPHS

#clustering  out of camp
library(arm)
m1 <-lmer(gps_clusteringOUT~friend_strength+age2+foraging_strength+kin_degree+(1|CAMP), data=all14,REML = T, na.action=na.fail)
vif.mer(m1)
sims<-sim(m1,n.sims=200)

new_temp<-data.frame(friend_strength=seq(-2,4,length=60),foraging_strength=mean(all14$foraging_strength),kin_degree=mean(all14$kin_degree),age2=mean(all14$age2))
pred_temp<-predict(m1,newdata=new_temp,re.form=~0)
#compute bootstrapped confidence interval, see ?predict.merMod
ci_line<-bootMer(m1,FUN=function(.) predict(.,newdata=new_temp,re.form=~0),nsim=200)
ci_regT<-apply(ci_line$t,2,function(x) x[order(x)][c(5,195)])
par(mfrow=c(1,2))
plot(gps_clusteringOUT~friend_strength,all14,col="#999933",xlab="POPULARITY", ylab="OUF-OF-CAMP CLUSTERING")
lines(new_temp$friend_strength,pred_temp,lwd=2,col="#332288")
lines(new_temp$friend_strength,ci_regT[1,],lty=2,col="#332288")
lines(new_temp$friend_strength,ci_regT[2,],lty=2,col="#332288")
text(-1,0.4, "**",cex=2.5, col="#117733")

new_n<-data.frame(friend_strength=mean(all14$friend_strength),foraging_strength=seq(-2,4,length=20),kin_degree=mean(all14$kin_degree),age2=mean(all14$age2))
pred_n<-predict(m1,newdata=new_n,re.form=~0)
#compute bootstrapped confidence interval, see ?predict.merMod
ci_line<-bootMer(m1,FUN=function(.) predict(.,newdata=new_n,re.form=~0),nsim=200)
ci_regT<-apply(ci_line$t,2,function(x) x[order(x)][c(5,195)])
par(mfrow=c(1,2))
plot(gps_clusteringOUT~foraging_strength,all14,col="#999933",xlab="HUNTING REPUTATION", ylab="OUF-OF-CAMP CLUSTERING")
lines(new_n$foraging_strength,pred_n,lwd=2,col="#332288")
lines(new_n$foraging_strength,ci_regT[1,],lty=2,col="#332288")
lines(new_n$foraging_strength,ci_regT[2,],lty=2,col="#332288")
#model2<‐extend(model1, along="x", n=20)
#mtext(text="Regression curves with 95% confidence intervals",side=3,outer=FALSE,at=-3)
text(2.5,0.2, "**",cex=2.5, col="#117733")

#clustering  in camp
library(arm)
m1 <-lmer(gps_clusteringIN~friend_strength+age2+foraging_strength+kin_degree+(1|CAMP), data=all17,REML = T, na.action=na.fail)

sims<-sim(m1,n.sims=200)

new_temp<-data.frame(friend_strength=seq(-2,4,length=60),foraging_strength=mean(all17$foraging_strength),kin_degree=mean(all17$kin_degree),age2=mean(all17$age2))
pred_temp<-predict(m1,newdata=new_temp,re.form=~0)
#compute bootstrapped confidence interval, see ?predict.merMod
ci_line<-bootMer(m1,FUN=function(.) predict(.,newdata=new_temp,re.form=~0),nsim=200)
ci_regT<-apply(ci_line$t,2,function(x) x[order(x)][c(5,195)])
par(mfrow=c(1,2))
plot(gps_clusteringIN~friend_strength,all17,col="#999933",xlab="POPULARITY", ylab="IN-CAMP CLUSTERING")
lines(new_temp$friend_strength,pred_temp,lwd=2,col="#332288")
lines(new_temp$friend_strength,ci_regT[1,],lty=2,col="#332288")
lines(new_temp$friend_strength,ci_regT[2,],lty=2,col="#332288")
#text(2.5,-1.55, "*",cex=2.5, col="#117733")

new_n<-data.frame(friend_strength=mean(all17$friend_strength),foraging_strength=seq(-2,4,length=20),kin_degree=mean(all17$kin_degree),age2=mean(all17$age2))
pred_n<-predict(m1,newdata=new_n,re.form=~0)
#compute bootstrapped confidence interval, see ?predict.merMod
ci_line<-bootMer(m1,FUN=function(.) predict(.,newdata=new_n,re.form=~0),nsim=200)
ci_regT<-apply(ci_line$t,2,function(x) x[order(x)][c(5,195)])
par(mfrow=c(1,2))
plot(gps_clusteringIN~foraging_strength,all17,col="#999933",xlab="HUNTING REPUTATION", ylab="IN-CAMP CLUSTERING")
lines(new_n$foraging_strength,pred_n,lwd=2,col="#332288")
lines(new_n$foraging_strength,ci_regT[1,],lty=2,col="#332288")
lines(new_n$foraging_strength,ci_regT[2,],lty=2,col="#332288")
#model2<‐extend(model1, along="x", n=20)
#mtext(text="Regression curves with 95% confidence intervals",side=3,outer=FALSE,at=-3)
#text(2.5,0.30, "*",cex=2.5, col="#117733")


#betweenness  out of camp
library(arm)
m1 <-lmer(gps_betweennessOUT~friend_strength+age2+foraging_strength+kin_degree+(1|CAMP), data=all11,REML = T, na.action=na.fail)
sims<-sim(m1,n.sims=200)
new_temp<-data.frame(friend_strength=seq(-2,4,length=60),foraging_strength=mean(all11$foraging_strength),kin_degree=mean(all11$kin_degree),age2=mean(all11$age2))
pred_temp<-predict(m1,newdata=new_temp,re.form=~0)
#compute bootstrapped confidence interval, see ?predict.merMod
ci_line<-bootMer(m1,FUN=function(.) predict(.,newdata=new_temp,re.form=~0),nsim=200)
ci_regT<-apply(ci_line$t,2,function(x) x[order(x)][c(5,195)])
par(mfrow=c(1,2))
plot(gps_betweennessOUT~friend_strength,all11,col="#999933",xlab="POPULARITY", ylab="OUT-OF-CAMP BETWEENNESS")
lines(new_temp$friend_strength,pred_temp,lwd=2,col="#332288")
lines(new_temp$friend_strength,ci_regT[1,],lty=2,col="#332288")
lines(new_temp$friend_strength,ci_regT[2,],lty=2,col="#332288")
text(2.5,0.2, "***",cex=2.5, col="#117733")

new_n<-data.frame(friend_strength=mean(all11$friend_strength),foraging_strength=seq(-2,4,length=20),kin_degree=mean(all11$kin_degree),age2=mean(all11$age2))
pred_n<-predict(m1,newdata=new_n,re.form=~0)
#compute bootstrapped confidence interval, see ?predict.merMod
ci_line<-bootMer(m1,FUN=function(.) predict(.,newdata=new_n,re.form=~0),nsim=200)
ci_regT<-apply(ci_line$t,2,function(x) x[order(x)][c(5,195)])
par(mfrow=c(1,2))
plot(gps_betweennessOUT~foraging_strength,all11,col="#999933",xlab="HUNTING REPUTATION", ylab="OUT-OF-CAMP BETWEENNESS")
lines(new_n$foraging_strength,pred_n,lwd=2,col="#332288")
lines(new_n$foraging_strength,ci_regT[1,],lty=2,col="#332288")
lines(new_n$foraging_strength,ci_regT[2,],lty=2,col="#332288")
#model2<‐extend(model1, along="x", n=20)
#mtext(text="Regression curves with 95% confidence intervals",side=3,outer=FALSE,at=-3)
text(-0.7,0.17, "**",cex=2.5, col="#117733")



#betweeness in camp
m1 <-lmer(gps_betweennessIN~friend_strength+age2+foraging_strength+kin_degree+(1|CAMP), data=all10,REML = T, na.action=na.fail)
sims<-sim(m1,n.sims=200)
new_temp<-data.frame(friend_strength=seq(-2,4,length=60),foraging_strength=mean(all10$foraging_strength),kin_degree=mean(all10$kin_degree),age2=mean(all10$age2))
pred_temp<-predict(m1,newdata=new_temp,re.form=~0)
#compute bootstrapped confidence interval, see ?predict.merMod
ci_line<-bootMer(m1,FUN=function(.) predict(.,newdata=new_temp,re.form=~0),nsim=200)
ci_regT<-apply(ci_line$t,2,function(x) x[order(x)][c(5,195)])
par(mfrow=c(1,2))
plot(gps_betweennessIN~friend_strength,all10,col="#999933",xlab="POPULARITY", ylab="IN-CAMP BETWEENNESS")
lines(new_temp$friend_strength,pred_temp,lwd=2,col="#332288")
lines(new_temp$friend_strength,ci_regT[1,],lty=2,col="#332288")
lines(new_temp$friend_strength,ci_regT[2,],lty=2,col="#332288")
#text(2.5,-0.20, "ns",cex=1.5, col="#117733")

new_n<-data.frame(friend_strength=mean(all10$friend_strength),foraging_strength=seq(-2,4,length=20),kin_degree=mean(all10$kin_degree),age2=mean(all10$age2))
pred_n<-predict(m1,newdata=new_n,re.form=~0)
#compute bootstrapped confidence interval, see ?predict.merMod
ci_line<-bootMer(m1,FUN=function(.) predict(.,newdata=new_n,re.form=~0),nsim=200)
ci_regT<-apply(ci_line$t,2,function(x) x[order(x)][c(5,195)])
par(mfrow=c(1,2))
plot(gps_betweennessIN~foraging_strength,all10,col="#999933",xlab="HUNTING REPUTATION", ylab="IN-CAMP BETWEENNESS")
lines(new_n$foraging_strength,pred_n,lwd=2,col="#332288")
lines(new_n$foraging_strength,ci_regT[1,],lty=2,col="#332288")
lines(new_n$foraging_strength,ci_regT[2,],lty=2,col="#332288")
#model2<‐extend(model1, along="x", n=20)
#text(2.5,-0.21, "**",cex=2.5, col="#117733")

#DEGREE OUT of camp
m1 <-lmer(gps_degreeOUT~friend_strength+age2+foraging_strength+kin_degree+(1|CAMP), data=all10,REML = T, na.action=na.fail)

sims<-sim(m1,n.sims=200)


new_temp<-data.frame(friend_strength=seq(-2,4,length=60),foraging_strength=mean(all10$foraging_strength),kin_degree=mean(all10$kin_degree),age2=mean(all10$age2))
pred_temp<-predict(m1,newdata=new_temp,re.form=~0)
#compute bootstrapped confidence interval, see ?predict.merMod
ci_line<-bootMer(m1,FUN=function(.) predict(.,newdata=new_temp,re.form=~0),nsim=200)
ci_regT<-apply(ci_line$t,2,function(x) x[order(x)][c(5,195)])
par(mfrow=c(1,2))
plot(gps_degreeOUT~friend_strength,all10,col="#999933",xlab="POPULARITY", ylab="OUT-OF-CAMP DEGREE")
lines(new_temp$friend_strength,pred_temp,lwd=2,col="#332288")
lines(new_temp$friend_strength,ci_regT[1,],lty=2,col="#332288")
lines(new_temp$friend_strength,ci_regT[2,],lty=2,col="#332288")
text(2.4,-1.7, "**",cex=2.5, col="#117733")

new_n<-data.frame(friend_strength=mean(all10$friend_strength),foraging_strength=seq(-2,4,length=20),kin_degree=mean(all10$kin_degree),age2=mean(all10$age2))
pred_n<-predict(m1,newdata=new_n,re.form=~0)
#compute bootstrapped confidence interval, see ?predict.merMod
ci_line<-bootMer(m1,FUN=function(.) predict(.,newdata=new_n,re.form=~0),nsim=200)
ci_regT<-apply(ci_line$t,2,function(x) x[order(x)][c(5,195)])
par(mfrow=c(1,2))
plot(gps_degreeOUT~foraging_strength,all10,col="#999933",xlab="HUNTING REPUTATION", ylab="OUT-OF-CAMP DEGREE")
lines(new_n$foraging_strength,pred_n,lwd=2,col="#332288")
lines(new_n$foraging_strength,ci_regT[1,],lty=2,col="#332288")
lines(new_n$foraging_strength,ci_regT[2,],lty=2,col="#332288")
#model2<‐extend(model1, along="x", n=20)

#degree in camp
m1 <-lmer(gps_degreeIN~friend_strength+age2+foraging_strength+kin_degree+(1|CAMP), data=all10,REML = T, na.action=na.fail)

sims<-sim(m1,n.sims=200)

new_temp<-data.frame(friend_strength=seq(-2,4,length=60),foraging_strength=mean(all10$foraging_strength),kin_degree=mean(all10$kin_degree),age2=mean(all10$age2))
pred_temp<-predict(m1,newdata=new_temp,re.form=~0)
#compute bootstrapped confidence interval, see ?predict.merMod
ci_line<-bootMer(m1,FUN=function(.) predict(.,newdata=new_temp,re.form=~0),nsim=200)
ci_regT<-apply(ci_line$t,2,function(x) x[order(x)][c(5,195)])
par(mfrow=c(1,2))
plot(gps_degreeIN~friend_strength,all10,col="#999933",xlab="POPULARITY", ylab="IN-CAMP DEGREE")
lines(new_temp$friend_strength,pred_temp,lwd=2,col="#332288")
lines(new_temp$friend_strength,ci_regT[1,],lty=2,col="#332288")
lines(new_temp$friend_strength,ci_regT[2,],lty=2,col="#332288")
#text(2.5,5, "**",cex=2.5, col="#228B22")

new_n<-data.frame(friend_strength=mean(all10$friend_strength),foraging_strength=seq(-2,4,length=20),kin_degree=mean(all10$kin_degree),age2=mean(all10$age2))
pred_n<-predict(m1,newdata=new_n,re.form=~0)
#compute bootstrapped confidence interval, see ?predict.merMod
ci_line<-bootMer(m1,FUN=function(.) predict(.,newdata=new_n,re.form=~0),nsim=200)
ci_regT<-apply(ci_line$t,2,function(x) x[order(x)][c(5,195)])
par(mfrow=c(1,2))
plot(gps_degreeIN~foraging_strength,all10,col="#999933",xlab="HUNTING REPUTATION", ylab="IN-CAMP DEGREE")
lines(new_n$foraging_strength,pred_n,lwd=2,col="#332288")
lines(new_n$foraging_strength,ci_regT[1,],lty=2,col="#332288")
lines(new_n$foraging_strength,ci_regT[2,],lty=2,col="#332288")
#model2<‐extend(model1, along="x", n=20)
#text(2.5,5, "**",cex=2.5, col="#228B22")

# SUPPLEMENTARY MATERIAL
#calculating correlation between network metrics derived from proximity networks with 2m and 5m cut-off values

GPSkib=read.csv("~/Dropbox/aSTRESS HUMANS/ANALIZA/kibonola/GPS network Kibonola partial.csv",header=TRUE,row.names=1,check.names=FALSE) 
GPSkib<-as.matrix(GPSkib)
colnames(GPSkib)<-c("A1","A3","A3","A4","A5","A6","A7","A8","A9","A10")
rownames(GPSkib)<-c("A1","A3","A3","A4","A5","A6","A7","A8","A9","A10")
write.csv(GPSkib,"~/Dropbox/Piotr/chapter 4 data/for submission/for submission/Colette revisions/FOR SUBMISSION/Major revisions/files for publication/CAMP3gps2M.csv")  

GPSkib5M=read.csv("~/Dropbox/aSTRESS HUMANS/ANALIZA/kibonola/GPS network Kibonola partial5m.csv",header=TRUE,row.names=1,check.names=FALSE) 
GPSkib5M<-as.matrix(GPSkib5M)
colnames(GPSkib5M)<-c("A1","A3","A3","A4","A5","A6","A7","A8","A9","A10")
rownames(GPSkib5M)<-c("A1","A3","A3","A4","A5","A6","A7","A8","A9","A10")
write.csv(GPSkib5M,"~/Dropbox/Piotr/chapter 4 data/for submission/for submission/Colette revisions/FOR SUBMISSION/Major revisions/files for publication/CAMP3gps5M.csv") 

#gps_strength2
M<-as.matrix(CAMP1gps2M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
gps_strength1<-v
#gps_strength1<-scale(gps_strength1)
gps_strength1

M<-as.matrix(CAMP2gps2M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
gps_strength2<-v
#gps_strength2<-scale(gps_strength2)

M<-as.matrix(CAMP3gps2M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
gps_strength3<-v
#gps_strength3<-scale(gps_strength3)
gps_strength3
gps_strength<-c(gps_strength1,gps_strength2,gps_strength3)
length(gps_strength)


#gps_strength2
M<-as.matrix(CAMP1gps5M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
gps_strength1<-v
#gps_strength1<-scale(gps_strength1)
gps_strength1

M<-as.matrix(CAMP2gps5M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
gps_strength2<-v
#gps_strength2<-scale(gps_strength2)

M<-as.matrix(CAMP3gps5M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=1, type="in")
v<-v[,-1]
gps_strength3<-v
#gps_strength3<-scale(gps_strength3)
gps_strength3
gps_strength5<-c(gps_strength1,gps_strength2,gps_strength3)
length(gps_strength5)

#gps degree 2m

M<-as.matrix(CAMP1gps2M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
#b_endeshMen<-b_endeshMen[,3]
cvc<-mean(b)
cvc
#cvc<-mean(b_endeshMen)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
gps_strength1<-v
#gps_strength1<-scale(gps_strength1)

M<-as.matrix(CAMP2gps2M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
#b_endMen<-b_endMen[,3]
#cvc<-mean(b_endMen)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
gps_strength2<-v
#gps_strength2<-scale(gps_strength2)



M<-as.matrix(CAMP3gps2M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
#b_kibMen<-b_kibMen[,3]
cvc<-mean(b)
#cvc
#cvc<-mean(b_kibMen)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
gps_strength3<-v
#gps_strength3<-scale(gps_strength3)
gps_strength3
gps_degree<-c(gps_strength1,gps_strength2,gps_strength3)
length(gps_degree)

#gps degree 5m

M<-as.matrix(CAMP1gps5M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
#b_endeshMen<-b_endeshMen[,3]
cvc<-mean(b)
cvc
#cvc<-mean(b_endeshMen)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
gps_strength1<-v
#gps_strength1<-scale(gps_strength1)

M<-as.matrix(CAMP2gps5M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
#b_endMen<-b_endMen[,3]
#cvc<-mean(b_endMen)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
gps_strength2<-v
#gps_strength2<-scale(gps_strength2)

M<-as.matrix(CAMP3gps5M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
#b_kibMen<-b_kibMen[,3]
cvc<-mean(b)
#cvc
#cvc<-mean(b_kibMen)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-degree_w(b, measure=c("alpha"), alpha=0, type="in")
v<-v[,-1]
gps_strength3<-v
#gps_strength3<-scale(gps_strength3)
gps_strength3
gps_degree5<-c(gps_strength1,gps_strength2,gps_strength3)
length(gps_degree5)

#gps betweenness 2m

M<-as.matrix(CAMP1gps2M)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-betweenness_w(b,alpha=0.5)  
v<-v[,-1]
gps_betweenness1<-v
#gps_betweenness1<-scale(gps_betweenness1)

M<-as.matrix(CAMP2gps2M)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-betweenness_w(b,alpha=0.5) 
v
v<-v[,-1]
gps_betweenness2<-v
#gps_betweenness2<-scale(gps_betweenness2)
gps_betweenness2

M<-as.matrix(CAMP3gps2M)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-betweenness_w(b,alpha=0.5)  
v<-v[,-1]
gps_betweenness3<-v
#gps_betweenness3<-scale(gps_betweenness3)
gps_betweenness3

gps_betweenness<-c(gps_betweenness1,gps_betweenness2,gps_betweenness3)

#gps betweenness 2m

M<-as.matrix(CAMP1gps5M)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-betweenness_w(b,alpha=0.5)  
v<-v[,-1]
gps_betweenness1<-v
#gps_betweenness1<-scale(gps_betweenness1)

M<-as.matrix(CAMP2gps5M)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-betweenness_w(b,alpha=0.5) 
v
v<-v[,-1]
gps_betweenness2<-v
#gps_betweenness2<-scale(gps_betweenness2)
gps_betweenness2

M<-as.matrix(CAMP3gps5M)
M<-M/10000
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
v<-betweenness_w(b,alpha=0.5)  
v<-v[,-1]
gps_betweenness3<-v
#gps_betweenness3<-scale(gps_betweenness3)
gps_betweenness3

gps_betweenness5<-c(gps_betweenness1,gps_betweenness2,gps_betweenness3)

# GPS CLUSTERING 2m
M<-as.matrix(CAMP1gps2M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
b<-symmetrise_w(b)
v<-clustering_local_w(b, measure=c("am"))
#v[is.na(v)] <- 0
v<-v[,-1]
sd(v,na.rm=T)
gps_betweenness2<-v
v
#gps_betweenness2<-scale(gps_betweenness2)

M<-as.matrix(CAMP2gps2M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
b<-symmetrise_w(b)
v<-clustering_local_w(b, measure=c("am"))
v<-v[,-1]
v
sd(v,na.rm=T)
gps_betweenness1<-v
#gps_betweenness1<-scale(gps_betweenness1)

M<-as.matrix(CAMP3gps2M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
b<-symmetrise_w(b)
v<-clustering_local_w(b, measure=c("am")) 
v<-v[,-1]
sd(v,na.rm=T)
#v[is.na(v)] <- 0
gps_betweenness3<-v
#gps_betweenness3<-scale(gps_betweenness3)

gps_clustering<-c(gps_betweenness1,gps_betweenness2,gps_betweenness3)
gps_clustering

# GPS CLUSTERING 5m

M<-as.matrix(CAMP1gps5M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
b<-symmetrise_w(b)
v<-clustering_local_w(b, measure=c("am"))
v<-v[,-1]
sd(v,na.rm=T)
gps_betweenness1<-v
#gps_betweenness1<-scale(gps_betweenness1)

M<-as.matrix(CAMP2gps5M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
b<-symmetrise_w(b)
v<-clustering_local_w(b, measure=c("am"))
#v[is.na(v)] <- 0
v<-v[,-1]
sd(v,na.rm=T)
gps_betweenness2<-v

#gps_betweenness2<-scale(gps_betweenness2)
M<-as.matrix(CAMP3gps5M)
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<-b[,3]
cvc<-mean(b)
cvc
M[M<cvc]<-0
b<-web2edges(M, webName=NULL, weight.column=TRUE, both.directions=FALSE,is.one.mode=T, out.files=c("edges", "names", "groups")[1:2],return=T, verbose=T)
b<- as.tnet(b, type="weighted one-mode tnet")
b<-symmetrise_w(b)
v<-clustering_local_w(b, measure=c("am")) 
v<-v[,-1]
sd(v,na.rm=T)
#v[is.na(v)] <- 0
gps_betweenness3<-v
#gps_betweenness3<-scale(gps_betweenness3)

gps_clustering5<-c(gps_betweenness1,gps_betweenness2,gps_betweenness3)
gps_clustering5



cor.test(gps_strength,gps_strength5,method = "spearman")
cor.test(gps_degree,gps_degree5,method = "spearman")
cor.test(gps_betweenness,gps_betweenness5,method = "spearman")
cor.test(gps_clustering,gps_clustering5,method = "spearman")

#correlation between proximity matrices of Camp 3 derived from two different data collection periods:
library(vegan)
veg.dist <- vegdist(CAMP3gps_period2) # Bray-Curtis
env.dist <- vegdist(CAMP3gps_period2, "euclid")
mantel(veg.dist, env.dist)
mantel(veg.dist, env.dist, method="spear",permutations=100000)

mantel.test(CAMP3gps_period2,CAMP3gps_period2, nperm = 9999, graph = FALSE,
            alternative = "two.sided")
library(ade4)
gen <- quasieuclid(as.dist(CAMP3gps_period1))
geo <- quasieuclid(as.dist(CAMP3gps_period2))
r1<-mantel.rtest(geo,gen)
plot(r1)
r1
