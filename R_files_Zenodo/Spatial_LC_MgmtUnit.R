setwd("M:/Miller Lab/Projects/Dove Models/Spatial")
require(reshape)
require(R2WinBUGS)
require(spdep)
require(plyr)

DAT <- read.csv("Data/Wings.csv")    ###  Wing data sets
k <- DAT["MOLT"]
j = k=="+"|k=="*"
levels(k[,]) = c(levels(k[,]),"10")
k[j] = 10
DAT[,"MOLT"] <- droplevels(k)

#####  Crosswalk database has extra information that needs to be joined to the data file

LL <- read.csv("Data/Spatial_Link.csv")  ###  Latitude longitude for county centroids - based on FIPS ID
FIPS <- read.csv("Data/FIPS.csv")  ### FIPS IDs based on SNO and CNO (state and county numbers)
Region <- read.csv("Data/Region.csv")  ###  Assigns regions from latest Otis estimates based on states

#Load land-use covariates, "cov"
load("Data/AllFly_Prop.RData")

DAT <- subset(DAT, ST %in% Region[Region$REG %in% unique(Region$REG),2]) #Remove any states that aren't in hunting regions

#######  Gets rid of all wings collected after the first 2 weeks in September and cleans up the rest of the data
DAT <-  subset(DAT,MO == 9 & is.na(MOLT) == FALSE & MOLT != "U" & MOLT != "D" & MOLT != "u" & MOLT != "d" & MOLT != ""
               & is.na(AS) == FALSE & AS != "0" & AS != "+" & AS != "")
DAT <-  subset(DAT,DA < 15)

DAT<-droplevels(DAT)
DAT$MOLT<-as.numeric(as.character(DAT$MOLT))
DAT$AS<-as.numeric(as.character(DAT$AS))
####  hexagon centers

scalar = 1.5
xn <- ceiling(10/(scalar*sqrt(3)))
yn <- ceiling(11/scalar)

mn.x = -97
mn.y = 37

xys=matrix(ncol = 2)

for (a in c(-xn:xn)*sqrt(3)*scalar+mn.x){
  for (b in c(-yn:yn)*scalar+mn.y){
    xys = rbind(xys,matrix(c(b,a),1,2) )
}}

xys = xys[-1,]
for (a in (c((-xn-1):xn)+0.5)*sqrt(3)*scalar+mn.x){
  for (b in (c((-yn-1):yn)+0.5)*scalar+mn.y){
    xys = rbind(xys,matrix(c(b,a),1,2) )
}}

xys <- xys[order(xys[,1]),]
xys <- cbind(1:nrow(xys),xys)

LL <- data.frame(LL,cell = rep(0,nrow(LL)))

for (a in 1:nrow(LL)){
  LL[a,"cell"] <- xys[which.min((LL[a,1]-xys[,2])^2+(LL[a,2]-xys[,3])^2),1]
}

adj = dnearneigh(xys[,2:3],0,2.26)
num <- sapply(adj,length)
sumNeigh <- sum(num)
adj = unlist(adj)

######  Join additionaLLl fields to the data
#DAT <- join(DAT,FIPS)
DAT <- join(DAT,LL)
DAT <- join(DAT,Region,by = c("ST"))
DAT <-  subset(DAT,is.na(Lat)==0)
DAT[DAT[,"MOLT"]==8,"AS"] = 3 
DAT[DAT[,"MOLT"]==9,"AS"] = 3
DAT[DAT[,"MOLT"]==10,"AS"] = 3

####### Join land-use covariates
cov$CNO<-as.integer(as.character(cov$COUNTY)) #In dove dataset not all counties have names, have to use CNO
cov$YEAR<-as.integer(cov$Year)
test<- join(DAT, cov, by = c("ST","CNO","YEAR"))
test<- test[which(is.na(test$hurb)==FALSE),] #Missing some states, no covariates for 2007

#######  The casting formula creates a cross-tab querie that is used to generate estimates
#######  The values on the left hand side of the equation define the groups the estimates are calculated for

Age <- cast(test, FIPS + YEAR + cell + REG + lurb + hurb + bigag + agfood + natw + nato ~ AS, value = "ID", fun.aggregate = length)  ### when ST is used data is grouped by state and year
Age <- data.frame(Age[,c(1:10)], AD = Age[,11], TOT = Age[,11]+Age[,12]+Age[,13]) # Covariates, then add wings for all age groups
Age <- Age[Age[,11]>0,] #Remove any records with zero captured adults

Age$MU<-"x" #Assign Management Unit Based off of Region
Age$MU[which(as.numeric(Age$REG)==1| as.numeric(Age$REG)==2|as.numeric(Age$REG)==3|as.numeric(Age$REG)==4)]<-"CMU"
Age$MU[which(as.numeric(Age$REG)==5| as.numeric(Age$REG)==6|as.numeric(Age$REG)==7|as.numeric(Age$REG)==8
       |as.numeric(Age$REG)==9|as.numeric(Age$REG)==10)]<-"EMU"
Age$MU[which(as.numeric(Age$REG)==11| as.numeric(Age$REG)==12|as.numeric(Age$REG)==13)]<-"WMU"
Age$MU<-as.factor(Age$MU)

nobs <- nrow(Age)
ncell <- nrow(xys)
cell <- Age$cell
urb<- Age$lurb + Age$hurb
bigag <- Age$bigag
agfood <- Age$agfood
nato <- Age$nato
natw <- Age$natw
AD <- Age$AD
N <- Age$TOT
JV <- Age$TOT-Age$AD
region<-as.numeric(Age$MU)

ni = 10000
nb = 2000
nt = 2
nc = 3

#inits <- function() list (p = rbeta(ncell,10,6),m = 0.6)
inits <- function() list (m = 0, alpha=rep(0,ncell), beta_urb=rep(0,3), beta_bigag=rep(0,3),
                          beta_agfood=rep(0,3), beta_nato=rep(0,3),beta_natw=rep(0,3))

meanHab <- matrix(NA,3,6) #Containers for habitat values for plotting
maxHab <- matrix(NA,3,6)
minHab <- matrix(NA,3,6)
for (a in 1:3){
  meanHab[a,1:5] <- c(mean(urb[region==a]), mean(bigag[region==a]), #Take the mean habitat values for each mgmt unit
    mean(agfood[region==a]), mean(natw[region==a]), mean(nato[region==a]))
  maxHab[a,1:6] <- c(quantile(urb[region==a],0.99), quantile(bigag[region==a],0.99), #Take the max habitat values for each mgmt unit
                     quantile(agfood[region==a],0.99), quantile(natw[region==a],0.99), quantile(nato[region==a],0.99), 
                     quantile(1 - urb[region==a] - bigag[region==a] - agfood[region==a] - natw[region==a] - nato[region==a],0.99))
  minHab[a,1:6] <- c(quantile(urb[region==a],0.01), quantile(bigag[region==a],0.01), #Take the min habitat values for each mgmt unit
                     quantile(agfood[region==a],0.01), quantile(natw[region==a],0.01), quantile(nato[region==a],0.01), 
                     quantile(1 - urb[region==a] - bigag[region==a] - agfood[region==a] - natw[region==a] - nato[region==a],0.01))
}
meanHab[,6] <- 1- rowSums(meanHab[,1:5])
derHab <- array(NA,c(21,6,6,3))
#For plotting results, need increases in one habitat to come at decreases in others
for (a in 1:6){
  for (b in 1:3){
    out <- matrix(NA,21,6)
    out[,a] <- seq(0,1,0.05)
    out[,-a] <- matrix(1-seq(0,1,0.05),21,1)%*%matrix(meanHab[b,-a]/sum(meanHab[b,-a]),1,5)
    derHab[,,a,b] <- out #makes matrix of plotting proportions of the different habitats
  }
}

data <- list(AD = AD, N=N, nobs = nobs, cell = cell, ncell= ncell,sumNeigh = sumNeigh, 
             adj = adj, num=num, region=region, urb=urb, bAg = bigag, 
             food = agfood, wood = natw, nato = nato , derHab11 = derHab[,,1,1],
             derHab21 = derHab[,,2,1], derHab31 = derHab[,,3,1], derHab41 = derHab[,,4,1], derHab51 = derHab[,,5,1],
             derHab61 = derHab[,,6,1], derHab12 = derHab[,,1,2], derHab22 = derHab[,,2,2], derHab32 = derHab[,,3,2],
             derHab42 = derHab[,,4,2], derHab52 = derHab[,,5,2], derHab62 = derHab[,,6,2], derHab13 = derHab[,,1,3],
             derHab23 = derHab[,,2,3], derHab33 = derHab[,,3,3], derHab43 = derHab[,,4,3], derHab53 = derHab[,,5,3], 
             derHab63 = derHab[,,6,3])

params <- c("AR","tau.o","spacetau","sigma.o","spacesigma","alpha",
            "beta_urb","beta_bigag","beta_agfood","beta_nato",
            "beta_natw","beta_0","outUrb","outAgfood","outAg",
            "outNatw","outNato","outOth")


out <- bugs(data = data, inits = inits, params, model.file ="C:/Users/djm516/Documents/GitHub/Dove_Recruitment/Models/Spatial_LC_MgmtUnit.txt",
            n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=FALSE, bugs.directory="C:/Users/djm516/Downloads/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14")

sum.out <- data.frame(AR=out$summary[1:333,1],cell =xys[,1],xys[,2:3])
l <- cast(Age, cell ~ .,value = "TOT" ,fun.aggregate = sum)
sum.out <- join(sum.out,l)
              
save.image(file = "Results/MgmtUnit_Data&Results.Rdata")
save(out,sum.out,file = "Results/MgmtUnit_Results_Final.Rdata")