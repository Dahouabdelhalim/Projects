library(coda)
library(geiger)
library(MCMCglmm)
library(HDInterval)

setwd("/scratch/osalehzi/phy/Final/")

#load the tree and comparative data
t <- read.tree('input/Aphidoidea-MCC.tre')
d <- read.csv('input/aphidtree.csv')

# modify mcc tip labels,
phy.names<-read.csv("input/phy-names.csv",header = T)
t$tip.label<-as.character(phy.names$tip_labels)

# do some formatting, drop species with missing data and get the intersection of our tree and trait datasets
rownames(d) <- d$id
keepers <- c('wings', 'host', 'db', 'woody', 'ant')
dat <- d[keepers]
dat <- na.omit(dat)
td <- treedata(t, dat, warnings=F)
tree1 <- td$phy
comp.dat <- data.frame(td$data)
comp.dat$animal <- rownames(comp.dat)
rownames(comp.dat) <- NULL

comp.dat$wings<-as.factor(comp.dat$wings)
comp.dat$host<-as.factor(comp.dat$host)
comp.dat$woody<-as.factor(comp.dat$woody)
comp.dat$ant<-as.factor(comp.dat$ant)
comp.dat$db<-log(comp.dat$db)
comp.dat$db<-(comp.dat$db-mean(comp.dat$db))/sd(comp.dat$db)


INtree <- inverseA(tree1,nodes="TIPS")

#set the prior
prior <-list(R=list(V=diag(2),nu=2), G=list(G1=list(V=diag(2),nu=2)))
prior2<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))
prior3 = list(R = list(V = diag(2), n = 2),G = list(G1 = list(V = diag(2)/3, n = 2)))

##use MCMC to estimate the parameters of generalized linear mixed models
m1 <-MCMCglmm(wings ~ host + db + woody + ant, random=~animal, 
              data=comp.dat, ginverse=list(animal=INtree$Ainv),
              prior=prior2,family="categorical", nitt=10000000, thin=2000, burnin=100000)
summary(m1)

save.image("models/pglmm/pglmm10M.RData")