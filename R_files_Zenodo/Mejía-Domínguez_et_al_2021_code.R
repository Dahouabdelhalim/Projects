####################################################
###   Mejía-Domínguez, N.R. et al. 2021 JVS      ###
####################################################

library(likelihood)
library(nlme)

# read csv file 'data_dbh>2.5' page from Mejía-Domínguez_et_al_2021_Data.xlsx
neighbors<-read.csv("'data_dbh>2.5.csv", header = TRUE)
# selection of tree species neighborhood
neighbors$SP <- CODI

# read csv file 'Seed_rain', 'Seed_bank' or 'Saplings_2007_2008' pages from Mejía-Domínguez_et_al_2021_Data.xlsx
seedrain <- read.csv("Seed_rain.csv", header = TRUE)
targets<-subset(seedrain,x>=15 & x<=85 & y>=15 & y<=85)
# define targets numbers: seed rain = 20, seed bank = 47 and saplings = 81
ntarget <- 20

# selection target stage
TARGET<-"Srain"
# selection target species
targets$variable <- targets$CHPE

# define neighborhood radius (5,10 or 15)
sr<- 5
RADIO<-5

# initial values parameters
ainitial<- 5
BAMax_X <-10
DISTMax_X <-50
BAdistMax_X <-10
NIMax_X <- 100
# for weibull function
bWei <- 1
cWei <- -4
# for log-normal function
blogBA<-2
clogBA<-1 
blogDIST<-15
clogDIST<-1
blogBAdist<-5
clogBAdist<-1
blogNI<-2
clogNI<-1
# for gaussian function
bgauBA<-2
cgauBA<-1
bgauDIST<-1.55
cgauDIST<-7
bgauBAdist<-5
cgauBAdist<-0.9
bgauNI<-5
cgauNI<-1

##################################################
###### --------------------------------  #########
# function to calculate distances of neighbors to target
neighdist<-function(targetx, targety, neighborx, neighbory) {
  sqrt((targetx-neighborx)^2 + (targety-neighbory)^2)       }  

# determine the maximum number of neighbors for any target in a square
max.neighbors.l <- 0
max.neigh <- vector ("numeric",length=ntarget)
for (i in 1:ntarget)
{
  max.neighbors <- max(max.neighbors.l,
                       nrow(subset(neighbors,  SP_spp != 0 &
                                     x < targets$x[i] + sr &
                                     x > targets$x[i] - sr &
                                     y < targets$y[i] + sr &
                                     y > targets$y[i] - sr 
                       )))
  max.neigh[i]<-max.neighbors
}
maxmax.neigh<-max(max.neigh)
max.neighbors<-maxmax.neigh
max.neighbors
# initialize matrices for species, distance and dbh to neighbors
distances <- matrix(0, nrow=nrow(targets), ncol=max.neighbors)
dbh.m     <- matrix(0, nrow=nrow(targets), ncol=max.neighbors)
BA.m      <- matrix(0, nrow=nrow(targets), ncol=max.neighbors)
BA.mh     <- matrix(0, nrow=nrow(targets), ncol=max.neighbors)
species   <- matrix(0, nrow=nrow(targets), ncol=max.neighbors)
sp        <- matrix(0, nrow=nrow(targets), ncol=max.neighbors)
# Go back and populate matrix
for (i in 1:ntarget) {
  neighbors.for.point <- subset(neighbors,  SP_spp != 0 &
                                  x < targets$x[i] + sr &
                                  x > targets$x[i] - sr &
                                  y < targets$y[i] + sr &
                                  y > targets$y[i] - sr)
  if (nrow(neighbors.for.point) > 0) {
    distances[i ,1:nrow(neighbors.for.point)] <-
      neighdist(targets$x[i], targets$y[i],
                neighbors.for.point$x, neighbors.for.point$y)
    dbh.m  [i,1:nrow(neighbors.for.point)]<- neighbors.for.point$dbh 
    
    # use whatever variable name you have for dbh and species codes
    BA.m  [i,1:nrow(neighbors.for.point)]<- pi*(neighbors.for.point$dbh/2)^2
    BA.mh [i,1:nrow(neighbors.for.point)]<- ((pi*(neighbors.for.point$dbh/2)^2)*10000)/pi*sr^2
    species[i ,1:nrow(neighbors.for.point)]<- neighbors.for.point$SP_spp
    sp[i ,1:nrow(neighbors.for.point)]<- neighbors.for.point$SP
  }
}
# Replace 0 distances with NAs
distances <- ifelse(distances == 0, NA, distances)
# drop neighbors in the square that are outside the desired neighbor radius, for made a circle
radius <- sr
distances <- ifelse(distances > radius, NA, distances)
# Replace 0 dbh with NAs
dbh.m <- ifelse(dbh.m == 0, NA, dbh.m)
# drop neighbors in the square that are outside the desired neighbor radius, for made a circle
radius <- sr
dbh.m <- ifelse(dbh.m > radius, NA, dbh.m)

### Lambdas setup ###
nsp <- length(levels(as.factor(species)))-1 
top.neighbors <- sort(table(ifelse(species==0,NA,species)), decreasing=T)[1:nsp]
top.neighbors.sppcodes <- as.factor(names(top.neighbors))
neigh.lambda.indexes <- matrix(0, nrow=nrow(species), ncol=ncol(species))
num <- length(levels(as.factor(species)))
neigh.lambda.indexes <- apply(species, c(1,2),
                              function(x) {if (x %in% top.neighbors.sppcodes) which(top.neighbors.sppcodes
                                                                                    == x) else NA})
#
initial.lambdas <- rep(1,nsp)
lambda.lower.bounds<-rep(-1,nsp)
lambda.upper.bounds<-rep(1,nsp)
names(initial.lambdas)<-c(as.character(top.neighbors.sppcodes))
names(lambda.lower.bounds)<-c(as.character(top.neighbors.sppcodes))
names(lambda.upper.bounds)<-c(as.character(top.neighbors.sppcodes))

# define pdf error
mypdf<-dpois

##############################
# Models
##############################
###--- Null model
null_model <- function(a){
  a
}

var <- list(lambda="predicted",x="variable",log=TRUE)
par <- list(a=ainitial)
par_lo <- list(a=1)
par_hi <- list(a=10000)
targets.na <- subset(targets, !is.na(targets$variable))
anneal(null_model,par, var, targets.na, par_lo, par_hi, mypdf,"variable",hessian=TRUE,max_iter=10000, show_display = FALSE)


###--- Models with lambdas
###--- Gaussian NI
GauModel_NI <- function(a,b,c,alpha,beta,lambdas){
  lambda.vals <- lambdas[neigh.lambda.indexes]
  dim(lambda.vals) <- dim(neigh.lambda.indexes)
  XX <- rowSums(lambda.vals*(dbh.m^alpha)/(distances^beta), na.rm=T)
  a * exp (-0.5 * ((XX - b)/c)^2 ) }

var <- list(lambda="predicted",x="variable",log=TRUE)
par <- list(a=ainitial, b=bgauNI, c=cgauNI,alpha=0, beta=1,lambdas = initial.lambdas)
par_lo <- list(a=1, b=0.01, c=0.01,alpha=0, beta=0,lambdas = lambda.lower.bounds)
par_hi <- list(a=10000, b= bgauhiNI, c=cgauhiNI,alpha=4, beta=4,lambdas = lambda.upper.bounds)
targets.na <- subset(targets, !is.na(targets$variable))
anneal(GauModel_NI, par, var, targets, par_lo, par_hi, mypdf,"variable",hessian=TRUE,max_iter=25000, show_display = FALSE)


###--- Lognormal NI
LognorModel_NI <- function(a,b,c,alpha,beta,lambdas){
  lambda.vals <- lambdas[neigh.lambda.indexes]
  dim(lambda.vals) <- dim(neigh.lambda.indexes)
  XX <- rowSums(lambda.vals*(dbh.m^alpha)/(distances^beta), na.rm=T)
  a * exp ( -0.5 * ( log (XX/b) / c )^2 ) }

var <- list(lambda="predicted",x="variable",log=TRUE)
par <- list(a=ainitial, b=blogNI, c=clogNI,alpha=1, beta=1,lambdas = initial.lambdas)
par_lo <- list(a=1, b=0.01, c=0.01,alpha=0, beta=0,lambdas = lambda.lower.bounds)
par_hi <- list(a=10000, b= bloghiNI, c=cloghiNI,alpha=4, beta=4,lambdas = lambda.upper.bounds)
targets.na <- subset(targets, !is.na(targets$variable))
anneal(LognorModel_NI, par, var, targets, par_lo, par_hi, mypdf,"variable",hessian=TRUE,max_iter=25000, show_display = FALSE)

###--- Weibull NI
WeiModel_NI <- function(a,b,c,alpha,beta,lambdas){
  lambda.vals <- lambdas[neigh.lambda.indexes]
  dim(lambda.vals) <- dim(neigh.lambda.indexes)
  XX <- rowSums(lambda.vals*(dbh.m^alpha)/(distances^beta), na.rm=T)
  
  a  * exp ( - b * XX^c )
}

var <- list(lambda="predicted",x="variable",log=TRUE)
par <- list(a=ainitial, b=bWei, c=cWei,alpha=0, beta=1,lambdas = initial.lambdas)
par_lo <- list(a=1, b=0.01, c=-10,alpha=0, beta=0,lambdas = lambda.lower.bounds)
par_hi <- list(a=10000, b= 50, c=10,alpha=4, beta=4,lambdas = lambda.upper.bounds)
targets.na <- subset(targets, !is.na(targets$variable))
anneal(WeiModel_NI, par, var, targets, par_lo, par_hi, mypdf,"variable",hessian=TRUE,max_iter=25000, show_display = FALSE)