library(phytools)
library(geiger)
library(HDInterval)

setwd("/scratch/osalehzi/phy/Final/")

#Use the new Bayesthresh function of phytools. Here's a totorial: http://www.phytools.org/eqg2015/threshold.html

#load the tree and comparative data
t <- read.tree('input/Aphidoidea-MCC.tre')
d <- read.csv('input/aphidtree.csv')

# modify mcc tip labels,
phy.names<-read.csv("input/phy-names.csv",header = T)
t$tip.label<-as.character(phy.names$tip_labels)

##do some formatting, drop species with missing data and get the intersection of our tree and trait datasets
rownames(d) <- d$id
keepers <- c('wings', 'woody')
dat <- d[keepers]
dat <- na.omit(dat)
td <- treedata(t, dat, warnings=F)
tree <- td$phy
comp.dat <- data.frame(td$data)
comp.dat$animal <- rownames(comp.dat)

foo <- c('wings', 'woody')
data.pair <- comp.dat[foo]
dp <- as.matrix(data.pair)


## set some parameters for the MCMC
sample<-500
ngen<-10000000
burnin<-0.5*ngen
# changed matrix dp to data.pair because matrix coerced factors as character
m1 <- threshBayes(tree, data.pair, types=c('discrete', 'continuous'), ngen=ngen, control=list(sample=sample))
#check for stationarity

par(mfrow=c(2,1))
plot(m1$par[,"logL"] ~ m1$par[,"gen"],type="l", xlab="generation",ylab="logL")
abline(v=burnin, lty=3, col='green', lwd=2)
#the correlation coefficent is in a row called 'r' in the m1$par matrix
#estimate the correlation
r.value <- mean(m1$par[(burnin/sample+1):nrow(m1$par),"r"])
r.HPD <- hdi(m1$par[(burnin/sample+1):nrow(m1$par),"r"])
print("woody")
print(r.value)
print(r.HPD)
#look at the posterior distribution
plot(density(m1$par[(burnin/sample+1):nrow(m1$par),"r"],bw=0.1),xlab="r",main="posterior density for r")
#lines(rep(cov2cor(true.Cov)[1,3],2),par()$usr[3:4],lty="dashed")

save.image("models/thresh/thresh-woody10M.RData")
