######------ PART A: CROSS-VALIDATION METHOD THAT PRODUCES MEASURES OF ACCURACY, COMMISSION, AND OMISSION FOR MODEL PREDICTIONS TO 20% OF DATA ------######

######------ PART A, SECTION 1: BUILDING PREDICTIVE MODEL WITH 80% OF DATA ------######

library(jagsUI)
library(psych)

setwd("C:/Users/lpetracca/Desktop/git/Jaguar")

#Import detection history
as.data.frame(detect<-read.csv("Jaguar_Max6_FP.csv"))
detect[detect=="-"] <-NA 
detect<-apply(detect, 2, function(x){as.numeric(x)})
as.data.frame(detect)
y <- detect
dim(y)

#code produces full detection histories for 1225 sampled cells
y_orig_full <- as.data.frame(y)
y_orig_full <- y_orig_full[!is.na(y_orig_full$S.1),]

#code produces max by row for 1225 sampled cells
y_orig_max <- apply(y, 1, max, na.rm=TRUE)
positions <- which(y_orig_max ==1 | y_orig_max == 0)
length(positions)
y_orig_max <- y_orig_max[positions]
length(y_orig_max)

#y now has 1225 rows instead of 1442
y <- y[positions,]

#Site Covariates, with rows extracted only for surveyed cells
as.data.frame(site_cov<-read.csv("SiteCovariates.csv", row.names=1))
site_cov <- site_cov[positions,]

#standardization of site covariates
Elev<-as.numeric(scale(site_cov$Elev))
DistStrictPA<-as.numeric(scale(site_cov$Dist_StrictPA))
DistSettle<-as.numeric(scale(site_cov$Dist_MajorSettlement))
Prey<-as.numeric(scale(site_cov$Prey))
PercCover<-as.numeric(site_cov$PercCover/100)
Agropast<-as.numeric(site_cov$Agropast)
CorridorTest<-as.numeric(site_cov$CorridorTest)
Size<-as.numeric(scale(site_cov$Size))

#Import corridor number for detection probability (six columns), with rows extracted only for surveyed cells
as.data.frame(corridor<-read.csv("Corridor_Max6_FP.csv", row.names=1))
corridor[corridor=="-"]<-NA
corridor[(is.na(corridor))] <- 6
corridor<-apply(corridor, 2, function(x){as.numeric(x)}) 
corridor <- corridor[,1]
Corridor <- as.numeric(corridor+1)
Corridor <- Corridor[positions]
ncorridors <- max(as.numeric(Corridor))

#Import effort data, with rows extracted only for surveyed cells
as.data.frame(effort<-read.csv("Effort_Max6_FP.csv", row.names=1))
effort[effort=="-"]<-NA
effort<-apply(effort, 2, function(x){as.numeric(x)}) 
#apply mean value to missing values
effort[(is.na(effort))] <- 0.545
as.data.frame(effort)
effort <- effort[positions,]

#Import data on whether or not there were four than four interviews collected in the cell, with rows extracted only for surveyed cells
as.data.frame((morethanfour<-read.csv("MoreThanFour_v3.csv", row.names=1, header=TRUE)))
morethanfour[morethanfour=="-"]<-NA
morethanfour<-apply(morethanfour, 2, function(x){as.numeric(x)}) 
#apply mean value to missing values
morethanfour[(is.na(morethanfour))] <- 0.67755
morethanfour <- morethanfour[positions]

#Import field season data for detection (1 column), with rows extracted only for surveyed cells
as.data.frame(FieldSeason<-read.csv("FieldSeason.csv", header=FALSE))
FieldSeason<-apply(FieldSeason, 2, function(x){as.numeric(x)})
as.data.frame(FieldSeason)
FieldSeason <- FieldSeason[positions]
NFieldSeason <- max(as.numeric(FieldSeason))

#CREATING INDEX CALLED "CLASS" SUCH THAT EACH SITE IS ASSIGNED VALUES 1-5
#IMPORTANT TO BE 1-5 FOLLOWED BY 1-5 ETC SUCH THAT EACH CORRIDOR IS REPRESENTED IN EACH CLASS
set.seed=100
N=245
class <- vector(mode = "logical", length = 1225)
length(class)
count=1:5
for (i in 1:N){
  class[count] <- sample.int(5, size=5, replace=FALSE)  
  count=count+5
}
length(class)



# Define model
sink("crossval_predicty.txt")
cat("
    model {
    
    #Priors
    #these are the random intercepts on psi for corridor
    for (i in 1:ncorridors){
    alpha.occ[i] ~ dnorm(mu.int, tau.int)
    }
    
    #these are the random intercepts on p for field season
    for (i in 1:NFieldSeason){
    alpha.p[i] ~ dnorm(mu.p, tau.p)
    }
    
    #hyperparameter on psi
    mu.int ~ dnorm(0, 0.001)
    tau.int <- 1 / (sigma.int * sigma.int)
    sigma.int ~ dunif(0, 100)
    
    #hyperparameter on p
    mu.p ~ dunif(0,10)  
    tau.p <- 1 / (sigma.p * sigma.p)
    sigma.p ~ dunif(0, 1)
    
    #priors for beta parameters
    agropast.beta ~ dnorm(0, 0.001)
    distsettle.beta ~ dnorm(0, 0.001)
    perccover.beta ~ dnorm(0, 0.001)
    prey.beta ~ dnorm(0, 0.001)
    strictpa.beta ~ dnorm(0, 0.001)
    elev.beta ~ dnorm(0, 0.001)
    
    effort.beta ~ dnorm(0, 0.001)
    
    #prior for probability of false positives
    alpha.fp ~ dunif(-10, 0)
    #prior for beta on false positive probability
    morethanfour.beta ~ dunif(0, 10)
    
    
    # Likelihood
    
    for (i in 1:R) { #start initial loop over the R sites
    # True state model for the partially observed true state
    logitpsi[i]<- alpha.occ[corridor[i]] + agropast.beta * Agropast[i] + distsettle.beta * DistSettle[i] + perccover.beta * PercCover[i] + prey.beta * Prey[i] + strictpa.beta * DistStrictPA[i] + elev.beta * Elev[i] 
    # Code to avoid issues with logit    
    logitpsitrun[i]<-min(999,max(-999,logitpsi[i]))
    psi[i]<-1/(1+exp(-logitpsitrun[i]))
    z[i] ~ dbern(psi[i])		# True occupancy z at site i
    
    # Model for false positive probability, based on whether there were >4 interviews conducted in the unit
    logitfp[i]<- alpha.fp + morethanfour.beta * morethanfour[i]
    # Code to avoid issues with logit
    logitfptrun[i]<-min(999,max(-999,logitfp[i]))
    fp[i]<-1/(1+exp(-logitfptrun[i]))
    
    # Model for detection probability, based on field season and effort
    for (j in 1:T) { # nested random effects
    logitp[i,j]<-alpha.p[FieldSeason[i]] + effort.beta * effort[i,j]
    # Code to avoid issues with logit
    logitptrun[i,j]<-min(999,max(-999,logitp[i,j]))
    p[i,j]<-1/(1+exp(-logitptrun[i,j]))
    }}
    
    #this is where we build a model using 80% (980/1225) of the surveyed sites
    for(i in 1:980){
    for(j in 1:T){
    #we have z and p for every site
    #this is where all the model fitting is happening
    #covariates on logit p estimated for 80% of sites
    y[i,j] ~ dbern(eff.py[i,j])
    eff.py[i,j] <- z[i] * p[i,j] + (1-z[i]) * fp[i]
    }
    }
    
    #this is where we cast out the model to the remaining 20% of sites
    for(i in 981:1225){
    #we only want eff.predy for each interview that was conducted in the site (anywhere from 1-6 per site)
    for(j in 1:V[i]){
    pred.y[i-980,j] ~ dbern(eff.predy[i-980,j])
    eff.predy[i-980,j]<- z[i] * p[i,j] + (1-z[i]) * fp[i]
    }
    X[i] <- (V[i])
    #this is where y_pred_max is calculated on the site level
    #a site will get a value of 1 if it has a 1 in any of the interviews
    y_pred_max[i-980] <- max(pred.y[i-980,1:X[i]])
    }
    
    #getting at accuracy components
    # a = true positive
    a <- y_orig_max + y_pred_max == 2
    # b = true negative
    b <- y_orig_max + y_pred_max == 0
    # c = false negative
    c <- y_orig_max - y_pred_max == 1
    # d = false positive
    d <- y_orig_max - y_pred_max == -1
    
    truepos <- sum(a)
    
    trueneg <- sum(b)
    
    falseneg <- sum(c)
    
    falsepos <- sum(d)
    
    denom <- trueneg + falsepos + falseneg + truepos

    accuracy_rate <- (trueneg+truepos)/denom
    sensitivity <- truepos / (truepos + falseneg)
    specificity <- trueneg / (trueneg + falsepos)

    }
    ",fill=TRUE)
sink()

######------ PART A, SECTION 2: APPLYING MODEL TO EACH OF FIVE FOLDS ------######

#MCMC settings
ni=400000
nb=300000
nt=10
nc=3

parameters <- c("kappa", "accuracy_rate", "sensitivity", "specificity", "truepos", "trueneg", "falsepos", "falseneg", "y_pred_max")


#NUMBER ONE
ord <- c(which(class!=1),which(class==1))
y_orig_max = y_orig_max[class==1]

#puts original detection histories in order of class (class 1 at end)
y_orig_full <- y_orig_full[ord,]
#gets number of interviews per cell, with class 1 at end
sum_int <- as.vector(apply(y_orig_full, 1, function(x) sum(!is.na(x))))
#gets number of interviews per cell for class 1
sum_int_class <- sum_int[981:1225]

#getting detection histories for class 1
# y_orig_test <- y_orig_full[981:1225,]
data <- list (y = y[class!=1,], yvalidate = y[class==1,], effort=effort[ord,], morethanfour=morethanfour[ord], Agropast=Agropast[ord], DistSettle=DistSettle[ord], PercCover=PercCover[ord], Prey=Prey[ord], DistStrictPA = DistStrictPA[ord], Elev=Elev[ord], corridor = (Corridor[ord]), FieldSeason=as.numeric(FieldSeason[ord]),
              ncorridors = 12, NFieldSeason = max(as.numeric(FieldSeason)), R = dim(y)[1], T = dim(y)[2], V=sum_int, y_orig_max=y_orig_max)
zst <- apply(y, 1, max)
inits <- function(){list(z = zst, mu.int = rnorm(1,0,2), sigma.int= rlnorm(1), alpha.occ = rnorm(ncorridors, 0, 2), agropast.beta = rnorm(1, 0, 2), distsettle.beta = rnorm(1, 0, 2), perccover.beta = rnorm(1, 0, 2), prey.beta = rnorm(1, 0, 2), 
                         strictpa.beta = rnorm(1, 0, 2), elev.beta = rnorm(1, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out1 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_predicty.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt)
save(out1,file = "E:/Jaguar/Results_PredictY/XVal1.Rdata")
write.csv(out1$summary, file = "E:/Jaguar/Results_PredictY/XVal1.csv")


#NUMBER TWO
ord <- c(which(class!=2),which(class==2))
y_orig_max = y_orig_max[class==2]

#puts original detection histories in order of class (class 2 at end)
y_orig_full <- y_orig_full[ord,]
#gets number of interviews per cell, with class 2 at end
sum_int <- as.vector(apply(y_orig_full, 1, function(x) sum(!is.na(x))))
#gets number of interviews per cell for class 2
sum_int_class <- sum_int[981:1225]

data <- list (y = y[class!=2,], yvalidate = y[class==2,], effort=effort[ord,], morethanfour=morethanfour[ord], Agropast=Agropast[ord], DistSettle=DistSettle[ord], PercCover=PercCover[ord], Prey=Prey[ord], DistStrictPA = DistStrictPA[ord], Elev=Elev[ord], Corridor = (Corridor[ord]), FieldSeason=as.numeric(FieldSeason[ord]),
              ncorridors = 12, NFieldSeason = max(as.numeric(FieldSeason)), R = dim(y)[1], T = dim(y)[2], V=sum_int, y_orig_max=y_orig_max)
zst <- apply(y, 1, max)
inits <- function(){list(z = zst, mu.int = rnorm(12,0,2), sigma.int= rlnorm(12), agropast.beta = rnorm(12, 0, 2), distsettle.beta = rnorm(12, 0, 2), perccover.beta = rnorm(12, 0, 2), prey.beta = rnorm(12, 0, 2), 
                         strictpa.beta = rnorm(12, 0, 2), elev.beta = rnorm(12, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out2 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_predicty.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out2,file = "E:/Jaguar/Results_PredictY/XVal2.Rdata")
write.csv(out2$summary, file = "E:/Jaguar/Results_PredictY/XVal2.csv")


#NUMBER THREE
ord <- c(which(class!=3),which(class==3))
y_orig_max = y_orig_max[class==3]

#puts original detection histories in order of class (class 3 at end)
y_orig_full <- y_orig_full[ord,]
#gets number of interviews per cell, with class 3 at end
sum_int <- as.vector(apply(y_orig_full, 1, function(x) sum(!is.na(x))))
#gets number of interviews per cell for class 3
sum_int_class <- sum_int[981:1225]


data <- list (y = y[class!=3,], yvalidate = y[class==3,], effort=effort[ord,], morethanfour=morethanfour[ord], Agropast=Agropast[ord], DistSettle=DistSettle[ord], PercCover=PercCover[ord], Prey=Prey[ord], DistStrictPA = DistStrictPA[ord], Elev=Elev[ord], Corridor = (Corridor[ord]), FieldSeason=as.numeric(FieldSeason[ord]),
              ncorridors = 12, NFieldSeason = max(as.numeric(FieldSeason)), R = dim(y)[1], T = dim(y)[2], V=sum_int, y_orig_max=y_orig_max)
zst <- apply(y, 1, max)
inits <- function(){list(z = zst, mu.int = rnorm(12,0,2), sigma.int= rlnorm(12), agropast.beta = rnorm(12, 0, 2), distsettle.beta = rnorm(12, 0, 2), perccover.beta = rnorm(12, 0, 2), prey.beta = rnorm(12, 0, 2), 
                         strictpa.beta = rnorm(12, 0, 2), elev.beta = rnorm(12, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out3 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_predicty.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out3,file = "E:/Jaguar/Results_PredictY/XVal3.Rdata")
write.csv(out3$summary, file = "E:/Jaguar/Results_PredictY/XVal3.csv")


#NUMBER FOUR
ord <- c(which(class!=4),which(class==4))
y_orig_max = y_orig_max[class==4]

#puts original detection histories in order of class (class 4 at end)
y_orig_full <- y_orig_full[ord,]
#gets number of interviews per cell, with class 4 at end
sum_int <- as.vector(apply(y_orig_full, 1, function(x) sum(!is.na(x))))
#gets number of interviews per cell for class 4
sum_int_class <- sum_int[981:1225]


data <- list (y = y[class!=4,], yvalidate = y[class==4,], effort=effort[ord,], morethanfour=morethanfour[ord], Agropast=Agropast[ord], DistSettle=DistSettle[ord], PercCover=PercCover[ord], Prey=Prey[ord], DistStrictPA = DistStrictPA[ord], Elev=Elev[ord], Corridor = (Corridor[ord]), FieldSeason=as.numeric(FieldSeason[ord]),
              ncorridors = 12, NFieldSeason = max(as.numeric(FieldSeason)), R = dim(y)[1], T = dim(y)[2], V=sum_int, y_orig_max=y_orig_max)
zst <- apply(y, 1, max)
inits <- function(){list(z = zst, mu.int = rnorm(12,0,2), sigma.int= rlnorm(12), agropast.beta = rnorm(12, 0, 2), distsettle.beta = rnorm(12, 0, 2), perccover.beta = rnorm(12, 0, 2), prey.beta = rnorm(12, 0, 2), 
                         strictpa.beta = rnorm(12, 0, 2), elev.beta = rnorm(12, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out4 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_predicty.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out4,file = "E:/Jaguar/Results_PredictY/XVal4.Rdata")
write.csv(out4$summary, file = "E:/Jaguar/Results_PredictY/XVal4.csv")


#NUMBER FIVE
ord <- c(which(class!=5),which(class==5))
y_orig_max = y_orig_max[class==5]

#puts original detection histories in order of class (class 5 at end)
y_orig_full <- y_orig_full[ord,]
#gets number of interviews per cell, with class 5 at end
sum_int <- as.vector(apply(y_orig_full, 1, function(x) sum(!is.na(x))))
#gets number of interviews per cell for class 5
sum_int_class <- sum_int[981:1225]


data <- list (y = y[class!=5,], yvalidate = y[class==5,], effort=effort[ord,], morethanfour=morethanfour[ord], Agropast=Agropast[ord], DistSettle=DistSettle[ord], PercCover=PercCover[ord], Prey=Prey[ord], DistStrictPA = DistStrictPA[ord], Elev=Elev[ord], Corridor = (Corridor[ord]), FieldSeason=as.numeric(FieldSeason[ord]),
              ncorridors = 12, NFieldSeason = max(as.numeric(FieldSeason)), R = dim(y)[1], T = dim(y)[2], V=sum_int, y_orig_max=y_orig_max)
zst <- apply(y, 1, max)
inits <- function(){list(z = zst, mu.int = rnorm(12,0,2), sigma.int= rlnorm(12), agropast.beta = rnorm(12, 0, 2), distsettle.beta = rnorm(12, 0, 2), perccover.beta = rnorm(12, 0, 2), prey.beta = rnorm(12, 0, 2), 
                         strictpa.beta = rnorm(12, 0, 2), elev.beta = rnorm(12, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out5 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_predicty.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out5,file = "E:/Jaguar/Results_PredictY/XVal5.Rdata")
write.csv(out5$summary, file = "E:/Jaguar/Results_PredictY/XVal5.csv")


######------ PART A, SECTION 3: CALCULATION OF ACCURACY STATISTICS ACROSS ALL FIVE FOLDS ------######

load("E:/Jaguar/Results_PredictY/XVal1.Rdata")
load("E:/Jaguar/Results_PredictY/XVal2.Rdata")
load("E:/Jaguar/Results_PredictY/XVal3.Rdata")
load("E:/Jaguar/Results_PredictY/XVal4.Rdata")
load("E:/Jaguar/Results_PredictY/XVal5.Rdata")

save(out1,out2,out3,out4,out5,file = "E:/Jaguar/Results_PredictY/XVal_ALL.Rdata")

load("E:/Jaguar/Results_PredictY/XVal_ALL.Rdata")

all_sens <- c(out1$sims.list$sensitivity, out2$sims.list$sensitivity, out3$sims.list$sensitivity, out4$sims.list$sensitivity, out5$sims.list$sensitivity)
all_spec <- c(out1$sims.list$specificity, out2$sims.list$specificity, out3$sims.list$specificity, out4$sims.list$specificity, out5$sims.list$specificity)
all_acc_rate <- c(out1$sims.list$accuracy_rate, out2$sims.list$accuracy_rate, out3$sims.list$accuracy_rate, out4$sims.list$accuracy_rate, out5$sims.list$accuracy_rate)

mean(all_sens)
mean(all_spec)
mean(all_acc_rate)

sens_quant <- quantile(all_sens, probs = c(0.025, 0.975), na.rm = FALSE, names=TRUE)
spec_quant <- quantile(all_spec, probs = c(0.025, 0.975), na.rm = FALSE, names=TRUE)
acc_rate_quant <- quantile(all_acc_rate, probs = c(0.025, 0.975), na.rm = FALSE, names=TRUE)





######------ PART B: CROSS-VALIDATION METHOD (based on Hooten & Hobbs 2015, adapted by Miller & Grant 2015) TESTING IF INCLUSION OF RANDOM EFFECTS IMPROVES MODEL FIT ------######

######------ PART B, Section 1: BUILDING BASE JAGUAR MODEL ------######

library(jagsUI)
setwd("C:/Users/lpetracca/Desktop/git/Jaguar")
as.data.frame(site_cov<-read.csv("SiteCovariates.csv", row.names=1))
dim(site_cov)
#only including 1440 sites (instead of 1442) because 1440 is divisible among five folds
site_cov <- site_cov[1:1440,]
dim(site_cov)

Elev<-as.numeric(scale(site_cov$Elev))
DistStrictPA<-as.numeric(scale(site_cov$Dist_StrictPA))
DistSettle<-as.numeric(scale(site_cov$Dist_MajorSettlement))
Prey<-as.numeric(scale(site_cov$Prey))
PercCover<-as.numeric(site_cov$PercCover/100)
Agropast<-as.numeric(site_cov$Agropast)
CorridorTest<-as.numeric(site_cov$CorridorTest)
Size<-as.numeric(scale(site_cov$Size))

#Import detection history
as.data.frame(detect<-read.csv("Jaguar_Max6_FP.csv"))
detect[detect=="-"] <-NA 
detect<-apply(detect, 2, function(x){as.numeric(x)})
as.data.frame(detect)
y <- detect
y <- y[1:1440,]

#Import corridor number for detection probability (six columns)
as.data.frame(corridor<-read.csv("Corridor_Max6_FP.csv", row.names=1))
corridor[corridor=="-"]<-NA
corridor[(is.na(corridor))] <- 6
corridor<-apply(corridor, 2, function(x){as.numeric(x)}) 
as.data.frame(corridor)
Corridor <- as.numeric(corridor+1)
ncorridors <- max(as.numeric(Corridor))
Corridor <- Corridor[1:1440]

#Import effort data
as.data.frame(effort<-read.csv("Effort_Max6_FP.csv", row.names=1))
effort[effort=="-"]<-NA
effort<-apply(effort, 2, function(x){as.numeric(x)}) 
#apply mean to missing values
effort[(is.na(effort))] <- 0.545
as.data.frame(effort)
effort <- effort[1:1440,]

#Import data on whether or not there were four than four interviews collected in the cell
as.data.frame((morethanfour<-read.csv("MoreThanFour_v3.csv", row.names=1, header=TRUE)))
morethanfour[morethanfour=="-"]<-NA
morethanfour<-apply(morethanfour, 2, function(x){as.numeric(x)}) 
#apply mean to missing values
morethanfour[(is.na(morethanfour))] <- 0.67755
morethanfour <- as.numeric(morethanfour[1:1440])

#Import field season data for detection (1 column)
as.data.frame(FieldSeason<-read.csv("FieldSeason.csv", header=FALSE))
FieldSeason<-apply(FieldSeason, 2, function(x){as.numeric(x)})
as.data.frame(FieldSeason)
FieldSeason <- FieldSeason[1:1440]
NFieldSeason <- max(as.numeric(FieldSeason))


#CREATING INDEX CALLED "CLASS" SUCH THAT EACH SITE IS ASSIGNED VALUES 1-5
#IMPORTANT TO BE 1-5 FOLLOWED BY 1-5 ETC SUCH THAT EACH CORRIDOR IS REPRESENTED IN EACH CLASS
set.seed=100
N=288
class <- vector(mode = "logical", length = 1440)
length(class)
count=1:5
for (i in 1:N){
  class[count] <- sample.int(5, size=5, replace=FALSE)  
  count=count+5
}



#TESTING ASSIGNMENT OF Y WITH CLASS 1 ON TOP
y- as.data.frame(y)
dim(y)
ord <- c(which(class!=1),which(class==1))
y <- y[ord,] 
y
tail(y)

######------ PART B, Section 2: CREATING MODEL FILE WITH CORRIDOR AS RANDOM EFFECT ------######

sink("crossval_randomeffect.txt")
cat("
    model {
    
    #Priors
    #these are the random intercepts on psi for corridor
    for (i in 1:ncorridors){
    alpha.occ[i] ~ dnorm(mu.int, tau.int)
    }
    
    #these are the random intercepts on p for field season
    for (i in 1:NFieldSeason){
    alpha.p[i] ~ dnorm(mu.p, tau.p)
    }
    
    #hyperparameter on psi
    mu.int ~ dnorm(0, 0.001)
    tau.int <- 1 / (sigma.int * sigma.int)
    sigma.int ~ dunif(0, 100)
    
    #hyperparameter on p
    mu.p ~ dunif(0,10)  
    tau.p <- 1 / (sigma.p * sigma.p)
    sigma.p ~ dunif(0, 1)
    
    #priors for beta parameters
    agropast.beta ~ dnorm(0, 0.001)
    distsettle.beta ~ dnorm(0, 0.001)
    perccover.beta ~ dnorm(0, 0.001)
    prey.beta ~ dnorm(0, 0.001)
    strictpa.beta ~ dnorm(0, 0.001)
    elev.beta ~ dnorm(0, 0.001)
    
    effort.beta ~ dnorm(0, 0.001)
    
    #prior for probability of false positives
    alpha.fp ~ dunif(-10, 0)
    #prior for beta on false positive probability
    morethanfour.beta ~ dunif(0, 10)
    
    # Likelihood
    for (i in 1:R) { #start initial loop over the R sites
    # True state model for the partially observed true state
    logitpsi[i]<- alpha.occ[corridor[i]] + agropast.beta * Agropast[i] + distsettle.beta * DistSettle[i] + perccover.beta * PercCover[i] + prey.beta * Prey[i] + strictpa.beta * DistStrictPA[i] + elev.beta * Elev[i] 
    # Code to avoid issues with logit
    logitpsitrun[i]<-min(999,max(-999,logitpsi[i]))
    psi[i]<-1/(1+exp(-logitpsitrun[i]))
    z[i] ~ dbern(psi[i])		# True occupancy z at site i
    
    # Model for false positive probability, based on whether there were >4 interviews conducted in the unit
    logitfp[i]<- alpha.fp + morethanfour.beta * morethanfour[i]
    logitfptrun[i]<-min(999,max(-999,logitfp[i]))
    fp[i]<-1/(1+exp(-logitfptrun[i]))
    
    # Model for detection probability, based on field season and effort  
    for (j in 1:T) { # nested random effects
    logitp[i,j]<-alpha.p[FieldSeason[i]] + effort.beta * effort[i,j]
    logitptrun[i,j]<-min(999,max(-999,logitp[i,j]))
    p[i,j]<-1/(1+exp(-logitptrun[i,j]))
    }}
    

#this is where a model is built from 80% of the data (1152 sites)
    for(i in 1:1152){
    for(j in 1:T){
    #we have z (from line 107) and p (from line 116) for every site
    #this is where all the model fitting is happening
    #covariates on logit p estimated from first 1152
    y[i,j] ~ dbern(eff.py[i,j])
    eff.py[i,j] <- z[i] * p[i,j] + (1-z[i]) * fp[i]
    
    }
    }
    
#now the model is cast out to the remaining 20% of sites
    for(i in 1153:1440){
    for(j in 1:T){
    eff.predicted[(i-1152),j]<- z[i] * p[i,j] + (1-z[i]) * fp[i]
    }
    }
    
    # Derived quantities
    occ.fs <- sum(z[])			# Number of occupied sites among 150
    }
    ",fill=TRUE)
sink()



######------ PART B, Section 3: CREATING MODEL FILE WITH CORRIDOR AS FIXED EFFECT ------######

sink("crossval_fixedeffect.txt")
cat("
    model {
    
    #Priors
    for (i in 1:ncorridors){
    
    #here, corridor is a fixed effect, so mu is not a hyperparameter
    mu.int[i] ~ dunif(-100,100)
    tau.int[i] <- 1 / (sigma.int[i] * sigma.int[i])
    sigma.int[i] ~ dunif(0, 100)
    
    #each beta is estimated separately for each corridor
    agropast.beta[i] ~ dnorm(0, 0.001)
    distsettle.beta[i] ~ dnorm(0, 0.001)
    perccover.beta[i] ~ dnorm(0, 0.001)
    prey.beta[i] ~ dnorm(0, 0.001)
    strictpa.beta[i] ~ dnorm(0, 0.001)
    elev.beta[i] ~ dnorm(0, 0.001)
    }
    
    #these are the random intercepts on p for field season
    for (i in 1:NFieldSeason){
    alpha.p[i] ~ dnorm(mu.p, tau.p)
    }
    
    #hyperparameter on p
    mu.p ~ dunif(0,10)  
    tau.p <- 1 / (sigma.p * sigma.p)
    sigma.p ~ dunif(0, 1)
    
    effort.beta ~ dnorm(0, 0.001)
    
    #prior for probability of false positives
    alpha.fp ~ dunif(-10, 0)
    #prior for beta on false positive probability
    morethanfour.beta ~ dunif(0, 10)
    
    
    # Likelihood
    
    for (i in 1:R) { #start initial loop over the R sites
    # True state model for the partially observed true state
    logitpsi[i]<- mu.int[Corridor[i]] + agropast.beta[Corridor[i]] * Agropast[i] + distsettle.beta[Corridor[i]] * DistSettle[i] + perccover.beta[Corridor[i]] * PercCover[i] + prey.beta[Corridor[i]] * Prey[i] + strictpa.beta[Corridor[i]] * DistStrictPA[i] + elev.beta[Corridor[i]] * Elev[i]
    # Code to avoid issues with logit
    logitpsitrun[i]<-min(999,max(-999,logitpsi[i]))
    psi[i]<-1/(1+exp(-logitpsitrun[i]))
    z[i] ~ dbern(psi[i])		# True occupancy z at site i
    
    # Model for false positive probability, based on whether there were >4 interviews conducted in the unit
    logitfp[i]<- alpha.fp + morethanfour.beta * morethanfour[i]
    # Code to avoid issues with logit
    logitfptrun[i]<-min(999,max(-999,logitfp[i]))
    fp[i]<-1/(1+exp(-logitfptrun[i]))
    
    # Model for detection probability, based on field season and effort
    for (j in 1:T) { # nested random effects
    logitp[i,j]<-alpha.p[FieldSeason[i]] + effort.beta * effort[i,j]
    # Code to avoid issues with logit
    logitptrun[i,j]<-min(999,max(-999,logitp[i,j]))
    p[i,j]<-1/(1+exp(-logitptrun[i,j]))
    }}
    
    
#this is where a model is built from 80% of the data (1152 sites)
    for(i in 1:1152){
    for(j in 1:T){
    #we have z (from line 107) and p (from line 116) for every site
    #this is where all the model fitting is happening
    #covariates on logit p estimated from first 1152
    y[i,j] ~ dbern(eff.py[i,j])
    eff.py[i,j] <- z[i] * p[i,j] + (1-z[i]) * fp[i]
    
    }
    }
    
#now the model is cast out to the remaining 20% of sites
    for(i in 1153:1440){
    for(j in 1:T){
    eff.predicted[(i-1152),j]<- z[i] * p[i,j] + (1-z[i]) * fp[i]
    }
    }
    
    # Derived quantities
    occ.fs <- sum(z[])			# Number of occupied sites among 150
    }
    ",fill=TRUE)
sink()


######------ PART B, Section 2: RUNNING MODELS WITH CORRIDOR AS RANDOM EFFECT ------######

parameters <- c("eff.predicted")

# MCMC settings
nc <- 3
nb <- 40000
ni <- 400000
nt <- 6

#NUMBER ONE
ord <- c(which(class!=1),which(class==1))
#you need to set up y and yval
data <- list (y = y[class!=1,], yvalidate = y[class==1,], effort=effort[ord,], morethanfour=morethanfour[ord], Agropast=Agropast[ord], DistSettle=DistSettle[ord], PercCover=PercCover[ord], Prey=Prey[ord], DistStrictPA = DistStrictPA[ord], Elev=Elev[ord], corridor=as.numeric(Corridor[ord]), FieldSeason=as.numeric(FieldSeason[ord]),
              ncorridors = max(as.numeric(Corridor)), NFieldSeason = max(as.numeric(FieldSeason)), R = dim(y)[1], T = dim(y)[2])
zst <- apply(y, 1, max)	
inits <- function(){list(z = zst, mu.int = rnorm(1,0,2), sigma.int= rlnorm(1), alpha.occ = rnorm(ncorridors, 0, 2), agropast.beta = rnorm(1, 0, 2), distsettle.beta = rnorm(1, 0, 2), perccover.beta = rnorm(1, 0, 2), prey.beta = rnorm(1, 0, 2), 
                         strictpa.beta = rnorm(1, 0, 2), elev.beta = rnorm(1, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7),alpha.p = runif(NFieldSeason,0,10), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out1 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_randomeffect.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out1,file = "E:/Jaguar/Results/XVal1_RandomEffect.Rdata")
write.csv(out1$summary, file = "XVal1_RandomEffect.csv")

#NUMBER TWO
ord <- c(which(class!=2),which(class==2))
data <- list (y = y[class!=2,], yvalidate = y[class==2,], effort=effort[ord,], morethanfour=morethanfour[ord], Agropast=Agropast[ord], DistSettle=DistSettle[ord], PercCover=PercCover[ord], Prey=Prey[ord], DistStrictPA = DistStrictPA[ord], Elev=Elev[ord], corridor=as.numeric(Corridor[ord]), FieldSeason=as.numeric(FieldSeason[ord]),
              ncorridors = max(as.numeric(Corridor)), NFieldSeason = max(as.numeric(FieldSeason)), R = dim(y)[1], T = dim(y)[2])
zst <- apply(y, 1, max)	
inits <- function(){list(z = zst, mu.int = rnorm(1,0,2), sigma.int= rlnorm(1), alpha.occ = rnorm(ncorridors, 0, 2), agropast.beta = rnorm(1, 0, 2), distsettle.beta = rnorm(1, 0, 2), perccover.beta = rnorm(1, 0, 2), prey.beta = rnorm(1, 0, 2), 
                         strictpa.beta = rnorm(1, 0, 2), elev.beta = rnorm(1, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7),alpha.p = runif(NFieldSeason,0,10), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out2 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_randomeffect.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out2,file = "E:/Jaguar/Results/XVal2_RandomEffect.Rdata")
write.csv(out2$summary, file = "XVal2_RandomEffect.csv")


#NUMBER THREE
ord <- c(which(class!=3),which(class==3))
data <- list (y = y[class!=3,], yvalidate = y[class==3,], effort=effort[ord,], morethanfour=morethanfour[ord], Agropast=Agropast[ord], DistSettle=DistSettle[ord], PercCover=PercCover[ord], Prey=Prey[ord], DistStrictPA = DistStrictPA[ord], Elev=Elev[ord], corridor=as.numeric(Corridor[ord]), FieldSeason=as.numeric(FieldSeason[ord]),
              ncorridors = max(as.numeric(Corridor)), NFieldSeason = max(as.numeric(FieldSeason)), R = dim(y)[1], T = dim(y)[2])
zst <- apply(y, 1, max)	
inits <- function(){list(z = zst, mu.int = rnorm(1,0,2), sigma.int= rlnorm(1), alpha.occ = rnorm(ncorridors, 0, 2), agropast.beta = rnorm(1, 0, 2), distsettle.beta = rnorm(1, 0, 2), perccover.beta = rnorm(1, 0, 2), prey.beta = rnorm(1, 0, 2), 
                         strictpa.beta = rnorm(1, 0, 2), elev.beta = rnorm(1, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7),alpha.p = runif(NFieldSeason,0,10), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out3 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_randomeffect.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out3,file = file="E:/Jaguar/Results/XVal3_RandomEffect.Rdata")
write.csv(out3$summary, file = "XVal3_RandomEffect.csv")


#NUMBER FOUR
ord <- c(which(class!=4),which(class==4))
data <- list (y = y[class!=4,], yvalidate = y[class==4,], effort=effort[ord,], morethanfour=morethanfour[ord], Agropast=Agropast[ord], DistSettle=DistSettle[ord], PercCover=PercCover[ord], Prey=Prey[ord], DistStrictPA = DistStrictPA[ord], Elev=Elev[ord], corridor=as.numeric(Corridor[ord]), FieldSeason=as.numeric(FieldSeason[ord]),
              ncorridors = max(as.numeric(Corridor)), NFieldSeason = max(as.numeric(FieldSeason)), R = dim(y)[1], T = dim(y)[2])
zst <- apply(y, 1, max)	
inits <- function(){list(z = zst, mu.int = rnorm(1,0,2), sigma.int= rlnorm(1), alpha.occ = rnorm(ncorridors, 0, 2), agropast.beta = rnorm(1, 0, 2), distsettle.beta = rnorm(1, 0, 2), perccover.beta = rnorm(1, 0, 2), prey.beta = rnorm(1, 0, 2), 
                         strictpa.beta = rnorm(1, 0, 2), elev.beta = rnorm(1, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7),alpha.p = runif(NFieldSeason,0,10), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out4 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_randomeffect.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out4,file = "Results/XVal4_RandomEffect.Rdata")
write.csv(out4$summary, file = "XVal4_RandomEffect.csv")


#NUMBER FIVE
ord <- c(which(class!=5),which(class==5))
data <- list (y = y[class!=5,], yvalidate = y[class==5,], effort=effort[ord,], morethanfour=morethanfour[ord], Agropast=Agropast[ord], DistSettle=DistSettle[ord], PercCover=PercCover[ord], Prey=Prey[ord], DistStrictPA = DistStrictPA[ord], Elev=Elev[ord], corridor=as.numeric(Corridor[ord]), FieldSeason=as.numeric(FieldSeason[ord]),
              ncorridors = max(as.numeric(Corridor)), NFieldSeason = max(as.numeric(FieldSeason)), R = dim(y)[1], T = dim(y)[2])
zst <- apply(y, 1, max)	
inits <- function(){list(z = zst, mu.int = rnorm(1,0,2), sigma.int= rlnorm(1), alpha.occ = rnorm(ncorridors, 0, 2), agropast.beta = rnorm(1, 0, 2), distsettle.beta = rnorm(1, 0, 2), perccover.beta = rnorm(1, 0, 2), prey.beta = rnorm(1, 0, 2), 
                         strictpa.beta = rnorm(1, 0, 2), elev.beta = rnorm(1, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7),alpha.p = runif(NFieldSeason,0,10), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out5 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_randomeffect.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out5,file = "Results/XVal5_RandomEffect.Rdata")      
write.csv(out5$summary, file = "XVal5_RandomEffect.csv")


save(out1,out2,out3,out4,out5,file = "Results/XVal_RandomEffect.Rdata")



######------ PART B, Section 3: RUNNING MODELS WITH CORRIDOR AS FIXED EFFECT ------######

parameters <- c("eff.predicted")

# MCMC settings
nc <- 3
nb <- 40000
ni <- 400000
nt <- 6

#NUMBER ONE
ord <- c(which(class!=1),which(class==1))
data <- list(y = y[class!=1,], yvalidate = y[class==1,], effort=effort, morethanfour=morethanfour, Agropast=Agropast, DistSettle=DistSettle, PercCover=PercCover, Prey=Prey, DistStrictPA = DistStrictPA, Elev=Elev,
             ncorridors = max(as.numeric(Corridor)), NFieldSeason = max(as.numeric(FieldSeason)), FieldSeason=as.numeric(FieldSeason), Corridor=Corridor, R = dim(y)[1], T = dim(y)[2])
zst <- apply(y, 1, max)
inits <- function(){list(z = zst, mu.int = rnorm(12,0,2), sigma.int= rlnorm(12), agropast.beta = rnorm(12, 0, 2), distsettle.beta = rnorm(12, 0, 2), perccover.beta = rnorm(12, 0, 2), prey.beta = rnorm(12, 0, 2), 
                         strictpa.beta = rnorm(12, 0, 2), elev.beta = rnorm(12, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out1 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_fixedeffect.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out1,file = "E:/Jaguar/Results/XVal1_FixedEffect.Rdata")
write.csv(out1$summary, file = "XVal1_FixedEffect.csv")

#NUMBER TWO
ord <- c(which(class!=2),which(class==2))
data <- list(y = y[class!=2,], yvalidate = y[class==2,], effort=effort, morethanfour=morethanfour, Agropast=Agropast, DistSettle=DistSettle, PercCover=PercCover, Prey=Prey, DistStrictPA = DistStrictPA, Elev=Elev,
             ncorridors = max(as.numeric(Corridor)), NFieldSeason = max(as.numeric(FieldSeason)), FieldSeason=as.numeric(FieldSeason), Corridor=Corridor, R = dim(y)[1], T = dim(y)[2])
zst <- apply(y, 1, max)
inits <- function(){list(z = zst, mu.int = rnorm(12,0,2), sigma.int= rlnorm(12), agropast.beta = rnorm(12, 0, 2), distsettle.beta = rnorm(12, 0, 2), perccover.beta = rnorm(12, 0, 2), prey.beta = rnorm(12, 0, 2), 
                         strictpa.beta = rnorm(12, 0, 2), elev.beta = rnorm(12, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out2 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_fixedeffect.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out2,file = "E:/Jaguar/Results/XVal2_FixedEffect.Rdata")
write.csv(out2$summary, file = "XVal2_FixedEffect.csv")

#NUMBER THREE
ord <- c(which(class!=3),which(class==3))
data <- list(y = y[class!=3,], yvalidate = y[class==3,], effort=effort, morethanfour=morethanfour, Agropast=Agropast, DistSettle=DistSettle, PercCover=PercCover, Prey=Prey, DistStrictPA = DistStrictPA, Elev=Elev,
             ncorridors = max(as.numeric(Corridor)), NFieldSeason = max(as.numeric(FieldSeason)), FieldSeason=as.numeric(FieldSeason), Corridor=Corridor, R = dim(y)[1], T = dim(y)[2])
zst <- apply(y, 1, max)
inits <- function(){list(z = zst, mu.int = rnorm(12,0,2), sigma.int= rlnorm(12), agropast.beta = rnorm(12, 0, 2), distsettle.beta = rnorm(12, 0, 2), perccover.beta = rnorm(12, 0, 2), prey.beta = rnorm(12, 0, 2), 
                         strictpa.beta = rnorm(12, 0, 2), elev.beta = rnorm(12, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out3 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_fixedeffect.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out3,file = "E:/Jaguar/Results/XVal3_FixedEffect.Rdata")
write.csv(out3$summary, file = "XVal3_FixedEffect.csv")
out3

#NUMBER FOUR
ord <- c(which(class!=4),which(class==4))
data <- list(y = y[class!=4,], yvalidate = y[class==4,], effort=effort, morethanfour=morethanfour, Agropast=Agropast, DistSettle=DistSettle, PercCover=PercCover, Prey=Prey, DistStrictPA = DistStrictPA, Elev=Elev,
             ncorridors = max(as.numeric(Corridor)), NFieldSeason = max(as.numeric(FieldSeason)), FieldSeason=as.numeric(FieldSeason), Corridor=Corridor, R = dim(y)[1], T = dim(y)[2])
zst <- apply(y, 1, max)
zst[1059] <- 1
inits <- function(){list(z = zst, mu.int = rnorm(12,0,2), sigma.int= rlnorm(12), agropast.beta = rnorm(12, 0, 2), distsettle.beta = rnorm(12, 0, 2), perccover.beta = rnorm(12, 0, 2), prey.beta = rnorm(12, 0, 2), 
                         strictpa.beta = rnorm(12, 0, 2), elev.beta = rnorm(12, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out4 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_fixedeffect.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out4,file = "E:/Jaguar/Results/XVal4_FixedEffect.Rdata")
write.csv(out4$summary, file = "XVal4_FixedEffect.csv")
out4

#NUMBER FIVE
ord <- c(which(class!=5),which(class==5))
data <- list(y = y[class!=5,], yvalidate = y[class==5,], effort=effort, morethanfour=morethanfour, Agropast=Agropast, DistSettle=DistSettle, PercCover=PercCover, Prey=Prey, DistStrictPA = DistStrictPA, Elev=Elev,
             ncorridors = max(as.numeric(Corridor)), NFieldSeason = max(as.numeric(FieldSeason)), FieldSeason=as.numeric(FieldSeason), Corridor=Corridor, R = dim(y)[1], T = dim(y)[2])
zst <- apply(y, 1, max)
inits <- function(){list(z = zst, mu.int = rnorm(12,0,2), sigma.int= rlnorm(12), agropast.beta = rnorm(12, 0, 2), distsettle.beta = rnorm(12, 0, 2), perccover.beta = rnorm(12, 0, 2), prey.beta = rnorm(12, 0, 2), 
                         strictpa.beta = rnorm(12, 0, 2), elev.beta = rnorm(12, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}
out5 <- jagsUI(data = data, inits = inits, parameters.to.save = parameters, model.file ="crossval_fixedeffect.txt",n.chains =nc, n.iter =ni,n.burnin =nb, n.thin =nt, parallel=TRUE)
save(out5,file = "E:/Jaguar/Results/XVal5_FixedEffect.Rdata")
write.csv(out5$summary, file = "XVal5_FixedEffect.csv")


save(out1,out2,out3,out4,out5,file = "Results/XVal_FixedEffect.Rdata")


load("E:/Jaguar/Results/XVal1_RandomEffect.Rdata")
load("E:/Jaguar/Results/XVal2_RandomEffect.Rdata")
load("E:/Jaguar/Results/XVal3_RandomEffect.Rdata")
load("E:/Jaguar/Results/XVal4_RandomEffect.Rdata")
load("E:/Jaguar/Results/XVal5_RandomEffect.Rdata")
save(out1,out2,out3,out4,out5,file = "E:/Jaguar/Results/XVal_RandomEffect_ALL.Rdata")

load("E:/Jaguar/Results/XVal1_FixedEffect.Rdata")
load("E:/Jaguar/Results/XVal2_FixedEffect.Rdata")
load("E:/Jaguar/Results/XVal3_FixedEffect.Rdata")
load("E:/Jaguar/Results/XVal4_FixedEffect.Rdata")
load("E:/Jaguar/Results/XVal5_FixedEffect.Rdata")
save(out1,out2,out3,out4,out5,file = "E:/Jaguar/Results/XVal_FixedEffect_ALL.Rdata")


# objects(out1$sims.list)
# (out1$sims.list$eff.predicted)
# out1$summary

######------ PART B, Section 4: USE OF MILLER & GRANT'S (2015) CROSS VALIDATION APPROACH COMPARING DEVIANCE WITH AND WITHOUT HIERARCHICAL STRUCTURE ------######

#random effect model

load(file = "E:/Jaguar/Results/XVal_RandomEffect_ALL.Rdata")

as.data.frame(detect<-read.csv("Jaguar_Max6_FP.csv"))
detect[detect=="-"] <-NA 
detect<-apply(detect, 2, function(x){as.numeric(x)})
as.data.frame(detect)
y <- detect
y <- y[1:1440,]


set.seed=100
N=288
class <- vector(mode = "logical", length = 1440)
length(class)
count=1:5
for (i in 1:N){
  class[count] <- sample.int(5, size=5, replace=FALSE)  
  count=count+5
}

SSER = 0
devRandom = 0
for (l in 1:5) {
  yvalidate = y[class==l,]
  yvalidate = array(as.numeric(yvalidate),c(288,6))
  dat = get(paste("out",l,sep = ""))$sims.list$eff.predicted
  mn <- apply(dat,c(2,3),mean)
  head(mn)
  dim(mn)
  #SSER = SSER + sum((mn-yval)^2,na.rm = TRUE)
  devRandom = devRandom + -2*sum(log(mn^yvalidate*(1-mn)^(1-yvalidate)),na.rm = TRUE)
}
devRandom

#fixed effect model

load(file = "E:/Jaguar/Results/XVal_FixedEffect_ALL.Rdata")

as.data.frame(detect<-read.csv("Jaguar_Max6_FP.csv"))
detect[detect=="-"] <-NA 
detect<-apply(detect, 2, function(x){as.numeric(x)})
as.data.frame(detect)
y <- detect
y <- y[1:1440,]


set.seed=100
N=288
class <- vector(mode = "logical", length = 1440)
length(class)
count=1:5
for (i in 1:N){
  class[count] <- sample.int(5, size=5, replace=FALSE)  
  count=count+5
}

SSER = 0
devFixed = 0
for (l in 1:5) {
  yvalidate = y[class==l,]
  yvalidate = array(as.numeric(yvalidate),c(288,6))
  dat = get(paste("out",l,sep = ""))$sims.list$eff.predicted
  mn <- apply(dat,c(2,3),mean)
  head(mn)
  dim(mn)
  #SSER = SSER + sum((mn-yval)^2,na.rm = TRUE)
  devFixed = devFixed + -2*sum(log(mn^yvalidate*(1-mn)^(1-yvalidate)),na.rm = TRUE)
}
devFixed



