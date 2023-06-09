# This R script compares models of evolutionary rate shifts between 
# pollination syndromes: symmetric rates (b to h = h to b), all rates 
# different, or irreversible b to h. Then, using the preferred irreversible 
# model, it performs stochastic character mapping. It uses the phylogenetic 
# tree file "chronogram_v8_2.tre" and the floral trait factor analysis output 
# data file "famd_coords.csv"

# UPDATED IN REVIEW: Added Bisse to simultaneously account for impacts on 
# speciation and extinction when estimating irreversibility

library(geiger)
library(phytools)
library(diversitree)

setwd("~/Dropbox/Costus grant/Dena and Kathleen/FloralEvolution/Dryad_Deposit")

# chronogram
tree <- read.nexus("chronogram_v8_2.tre")
plot(tree, cex=0.4)

#################
## use syndromes from Random Forest predictions
################
# floral Dim characters
data <- read.csv("famd_coords.csv",row.names=1)

# match tree and data
(name.check(tree, data) -> phyOverlap)
drop.tip(tree, phyOverlap$tree_not_data) -> phyComparativeTree
data <- data[!(row.names(data)  %in% phyOverlap$data_not_tree), ]

#use RF syndrome predictions
data$...2[33] = "hummingbird" #Costus_sp_nov_19168
data$...2[39] = "hummingbird" #Costus_varzearum_19252
data$...2[51] = "hummingbird" #Costus_dirzoi_98079

# syndrome transition models
#prepare data
dat <- as.data.frame(data$...2)
rownames(dat) <- row.names(data)
synd_geiger <- treedata(phyComparativeTree, dat)
#symmetric rates
(synd_SYM_geiger <- fitDiscrete(synd_geiger$phy, synd_geiger$data [,1],
                                type="discrete", model = "SYM"))
#all rates
(synd_ARD_geiger <- fitDiscrete(synd_geiger$phy, synd_geiger$data [,1],
                                type="discrete", model = "ARD"))
# model: only B to H, but not reverse
model<-matrix(c(0,0,1,0),2,2)
colnames(model)<-rownames(model)<-c("bee","hummingbird")
model
(fitIrr<-fitDiscrete(synd_geiger$phy, synd_geiger$data [,1],model=model))
#aic weights
(aicc<-setNames(
  c(synd_SYM_geiger$opt$aicc,synd_ARD_geiger$opt$aicc,fitIrr$opt$aicc),
  c("ER","ARD","Irr")))
aic.w(aicc) #Irr model best
# ER       ARD       Irr 
# 0.3505889 0.1647360 0.4846752 

# stochastic character mapping of syndrome
fmode<-setNames(as.factor(data$...2),rownames(data))
x <- make.simmap(phyComparativeTree, x=fmode, model= model, nsim=100)
summary(x)
# 100 trees with a mapped discrete character with states:
#  bee, hummingbird 
#
# trees have 13.11 changes between states on average
#
# changes are of the following types:
#   bee,hummingbird hummingbird,bee
# x->y           13.11               0
#
# mean total time spent in each state is:
#  bee hummingbird     total
# raw  0.4986612   0.2781890 0.7768502
# prop 0.6419013   0.3580987 1.0000000

table(summary(x)$count) # range changes between states
#    0  13  14  15 
#  100 180  18   2 

obj <- summary(x)

pdf(file="Fig4phylo.pdf",width=11,height=11,paper='special') 
hh<-max(nodeHeights(phyComparativeTree))*0.02 # for label offset
plot(phyComparativeTree,no.margin=TRUE,label.offset=hh,edge.width=2)
nodelabels(pie=obj$ace,piecol=c("#1F639B","#ED553B"), cex=0.3)
tiplabels(pie=obj$tips, piecol=c("#1F639B","#ED553B"),cex=0.3)
dev.off()

###########
#BISSE
###########
library(tidyverse)

names(dat) <- c("syn")
dat$syn 

#format trait vector for bisse
dat.bisse <- dat %>% mutate(syn = case_when(syn=="bee" ~ "0", syn=="hummingbird" ~ "1"))
dat.bisse <- as.numeric(dat.bisse$syn)
names(dat.bisse) <- row.names(dat)

#bisse
lik <- make.bisse(phyComparativeTree, dat.bisse)
lik(pars) #-211.3333

#starting point for ML search
(p <- starting.point.bisse(phyComparativeTree))

#initiat ML search
fit <- find.mle(lik, p)
fit$lnLik #log-likelihood value
round(coef(fit), 3) #coefficients

# lambda: speciation
# mu: extinction
# q: transition x to y (0=bee, 1=bird)
# lambda0 lambda1     mu0     mu1     q01     q10 
# 73.954  46.780   0.000   0.000  28.197   0.000 

#test hypothesis that transition rates are different
lik.l <- constrain(lik, q01 ~ q10)
fit.l <- find.mle(lik.l, p[argnames(lik.l)])
fit.l$lnLik # 125.8509
round(rbind(full=coef(fit), equal.q=coef(fit.l, TRUE)), 3)
anova(fit, equal.l=fit.l) # not significantly better model
# although with this small a tree and 6 parameters that is not at all surprising

#test hypothesis that speciation rates are different
lik.l <- constrain(lik, lambda0 ~ lambda1)
fit.l <- find.mle(lik.l, p[argnames(lik.l)])
fit.l$lnLik # 125.8925
round(rbind(full=coef(fit), equal.l=coef(fit.l, TRUE)), 3)
anova(fit, equal.l=fit.l) # not significantly better model
# although with this small a tree and 6 parameters that is not at all surprising

#because we fitting 6 parameters with only 52 taxa
#priors needed (according to bisse tutorial) so that
#posterior distribution is proper
prior <- make.prior.exponential(1 / (2 * (p[1] - p[3])))
#establish parameter update method using a short chain run
set.seed(1)
tmp <- mcmc(lik, fit$par, nsteps=100, prior=prior,
            lower=0, w=rep(1, 6), print.every=0)
w <- diff(sapply(tmp[2:7], range))
w
#run chain 10K steps
samples <- mcmc(lik, fit$par, nsteps=10000, w=w, lower=0, prior=prior,
                print.every=0)

#make marginal plots showing rates: spec, ext, trans
pdf(file="~/Dropbox/Costus grant/Dena and Kathleen/FloralEvolution/Bisse/FigBisse.pdf",width=11,height=4,paper='special') 
par(mar=c(5,5,3,3), mfrow=c(1,3))
col <- c("#1F639B","#ED553B")
col.lines <- c("#ED553B", "#1F639B")

#speciation
profiles.plot(samples[c("lambda0", "lambda1")], col.line=col, las=1,
              xlab="Speciation rate",cex.lab=2)
legend("topright",legend=c("bee","bird"), col=col, pch=15, cex=1.1)
abline(v=c(coef(fit)[2], coef(fit)[1]), col=col.lines)

#extinction
profiles.plot(samples[c("mu0", "mu1")], col.line=col, las=1,
              xlab="Extinction rate",cex.lab=2, ylab="")
legend("topright",legend=c("bee","bird"), col=col, pch=15, cex=1.1)
abline(v=c(coef(fit)[4], coef(fit)[3]), col=col.lines)
mtext("Bisse model posterior probability distributions",line=1, cex=1.5)

#transition
profiles.plot(samples[c("q01", "q10")], col.line=col, las=1,
              xlab="Transition rate",cex.lab=2, ylab="")
legend("topright",legend=c("bee to bird","bird to bee"), col=col, pch=22, cex=1.1, pt.bg=col)
abline(v=c(coef(fit)[6], coef(fit)[5]), col=col.lines)

dev.off()

#################
## use syndromes from taxonomic treatments
################
#################
## Use only taxa with pollinator observations
#################
#pollinator observations (for subsetting dataset)
polobs <- read.csv("Pollinator_observations.csv")
polobs <- polobs[,c(1:3)]
polobs <- subset(polobs, observed=="yes")
rownames(polobs) <- polobs$tip

# match tree and data
(name.check(tree, polobs) -> phyOverlap)
drop.tip(tree, phyOverlap$tree_not_data) -> phyComparativeTree
data <- data[!(row.names(data)  %in% phyOverlap$data_not_tree), ]

#prepare data
dat <- as.data.frame(polobs$syndrome)
rownames(dat) <- row.names(polobs)
synd_geiger <- treedata(phyComparativeTree, dat)

#symmetric rates
(synd_SYM_geiger <- fitDiscrete(synd_geiger$phy, synd_geiger$data [,1],
                                type="discrete", model = "SYM"))
#all rates
(synd_ARD_geiger <- fitDiscrete(synd_geiger$phy, synd_geiger$data [,1],
                                type="discrete", model = "ARD"))
# model: only B to H, but not reverse
model<-matrix(c(0,0,1,0),2,2)
colnames(model)<-rownames(model)<-c("bee","hummingbird")
model
(fitIrr<-fitDiscrete(synd_geiger$phy, synd_geiger$data [,1],model=model))
#aic weights
(aicc<-setNames(
  c(synd_SYM_geiger$opt$aicc,synd_ARD_geiger$opt$aicc,fitIrr$opt$aicc),
  c("ER","ARD","Irr")))
aic.w(aicc) #Irr model best
#  ER       ARD       Irr 
# 0.1521715 0.2018753 0.6459532

# stochastic character mapping of syndrome
fmode<-setNames(as.factor(polobs$syndrome),rownames(polobs))
x <- make.simmap(phyComparativeTree, x=fmode, model= model, nsim=100)
summary(x) 
# 100 trees with a mapped discrete character with states:
#  bee, hummingbird 
#
# trees have 7.45 changes between states on average
#
# changes are of the following types:
#  bee,hummingbird hummingbird,bee
#  x->y            7.45               0

table(summary(x)$count) # range changes between states
#  0   7   8   9 
# 100 116  78   6 