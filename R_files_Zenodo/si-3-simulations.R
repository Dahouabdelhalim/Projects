# Code to repeat simulation analysis

#######################
#HEADERS###############
#######################

library(geiger)
library(pez)
source("si-2-functions.R")
library(picante)
library(parallel)
library(xtable)
library(phytools)
library(MuMIn)
library(mvMORPH)
library(OUwie)
library(RPANDA)

inv.logit <- function(x)
    exp(x) / (exp(x)+1)

find.clade <- function(tree, min.size, max.size, iter=100){
    clade <- rep(FALSE, length(tree$tip.label))
    for(i in seq_len(iter)){
        clade <- extract.clade(tree, sample(seq(from=length(tree$tip.label)+1,length.out=tree$Nnode-1),1))$tip.label
        clade <- tree$tip.label %in% clade
        if(sum(clade)>=min.size & sum(clade)<=max.size)
            return(clade)
    }
    return(NA)
}

null.var <- function(c.data, clade, n.reps){
    vars <- numeric(n.reps)
    for(i in seq_along(vars)){
        c.data$comm <- c.data$comm[,sample(ncol(c.data$comm))]
        t <- clade.var(c.data)
        vars[i] <- t$variance[clade]
    }
    return(vars)
}

make.comm <- function(tree, trait, n.sites, clade, rnd=FALSE){
    comm <- matrix(nrow=n.sites, ncol=length(trait), dimnames=list(seq_len(n.sites),tree$tip.label))
    for(i in seq_len(n.sites)){
        sp <- sample(length(tree$tip.label), 1)
        curr.trait <- abs(trait - trait[sp])
        curr.trait[curr.trait > 1] <- 1
        comm[i,] <- rbinom(length(tree$tip.label), 1, 1-curr.trait)
    }
    return(comparative.comm(tree, comm))
}

sim.trait <- function(tree, clade, tree.sd, clade.sd){
    trait <- sim.char(tree, clade.sd, model="BM")[,,1]
    null  <- sim.char(tree, tree.sd, model="BM")[,,1]
    # Prepare output and return
    return(ifelse(clade, trait, null))
}

calc.wrap <- function(i){
    tree <- sim.bdtree(n=data$n.spp[i])
    clade <- find.clade(tree, data$min.clade[i], data$max.clade[i], data$max.iter[i])
    if(any(is.na(clade)))
        return(setNames(c(NA,NA,NA,NA,NA,NA,NA), c("var","exp.var","mean.null","median.null","sd.null","p","n.clade")))
    clade <- setNames(clade, tree$tip.label)
    mrca.clade <- getMRCA(tree, which(clade))
    trait <- sim.trait(tree, clade, data$tree.sd[i], data$clade.sd[i])
    c.data <- make.comm(tree, trait, data$n.sites[i])
    var <- clade.var(c.data)
    null <- null.var(c.data, mrca.clade, 9999)

    return(setNames(
        c(var$variance[mrca.clade], var$expect.var[mrca.clade], mean(null), median(null), sd(null), rank(c(var$variance[mrca.clade],null))[1]/10000,sum(clade)),
        c("var","exp.var","mean.null","median.null","sd.null","p","n.clade")
    ))
}

#######################
# Simulate ############
#######################
data <- data.frame(expand.grid(n.spp=c(50,100), n.sites=c(50,100), min.clade=c(.05,.1), max.clade=c(.1,.2), tree.sd=seq(.5,2,.5), clade.sd=c(10^seq(-3,3,.25),rep(1,20)), iter=1, max.iter=10))
data <- data[data$max.clade > data$min.clade,]
data$min.clade <- round(data$min.clade * data$n.spp)
data$max.clade <- round(data$max.clade * data$n.spp)
data$clade.sd <- data$tree.sd * data$clade.sd

# This line assumes you have 12 cores on your computer
output <- mcMap(calc.wrap, seq_len(nrow(data)), mc.cores=12)

# This takes a while so you may want to save the workspace
save.image("si-3-simulations.RData")


#######################
# Model    ############
#######################
data <- cbind(data, do.call(rbind,output))
data <- na.omit(data)

data$ratio <- with(data, clade.sd/tree.sd)
data$result <- with(data, as.numeric(ifelse(ratio < 1, ifelse(p > .975, TRUE, FALSE), ifelse(p < .025, TRUE, FALSE))))
data$eff.size <- with(data, ifelse(ratio>1, clade.sd/tree.sd, tree.sd/clade.sd))
data$trt <- with(data, ifelse(ratio>1, "over", ifelse(ratio==1, "null", "under")))
cols <- setNames(c("orange","skyblue","grey70"), c("over","under","null"))
data$log10.eff.size <- log10(data$eff.size)

# Make summary plot
pdf("simulation-plot.pdf")
{
    # Plot basics
    par(mar=c(5.1,5.1,4.1,2.1))
    with(data, plot(p ~ eff.size, log="x", axes=FALSE, xlab=expression(paste("effect size (", sigma[clade]^2, ":", sigma[tree]^2, " or ", sigma[tree]^2, ":", sigma[clade]^2,")")), ylab="probability", cex.lab=1.5, type="n"))
    with(data[data$eff.size>1,], points(p ~ eff.size, pch=20, col=cols[trt]))
    with(data[data$eff.size==1,], points(p ~ jitter(eff.size), pch=20, col=cols[trt]))
    axis(1, cex=1.5, at=1, label="none")
    axis(1, cex=1.5, at=c(sort(unique(data$eff.size))[2], 5, 10, 50, 100, 500, 1000), label=c("",5,10,50,100,500,1000))
    axis(2, cex=1.5)

    # Overdispersion (fast trait)
    model <- glm(p ~ log10.eff.size, data=data[data$trt %in% c("over"),], family=quasibinomial)
    pred.dat <- data.frame(log10.eff.size=seq(0.25,3,.01))
    pred.dat$p <- predict(model, pred.dat, type="response")
    with(pred.dat, lines(p ~ I(10^log10.eff.size), col="red", lwd=3))

    # Underdispersion (slow trait)
    model <- glm(p ~ log10.eff.size, data=data[data$trt %in% c("under"),], family=quasibinomial)
    pred.dat <- data.frame(log10.eff.size=seq(0.25,3,.01))
    pred.dat$p <- predict(model, pred.dat, type="response")
    with(pred.dat, lines(p ~ I(10^log10.eff.size), col="blue", lwd=3))

    # Null (same rate)
    mean <- mean(data$p)
    lines(c(.9,1.1), rep(mean,2), col="black", lwd=3)
    under.null <- with(data[data$eff.size==1,], sum(p<.025)/length(p)) * 100
    lines(c(1.1,1.1), c(.025,0), col="black", lwd=3)
    text(1.15, 0, paste0(round(under.null,0),"%"), adj=0)
    over.null <- with(data[data$eff.size==1,], sum(p>.975)/length(p)) * 100
    lines(c(1.1,1.1), c(.975,1), col="black", lwd=3)
    text(1.15, 1, paste0(round(over.null,0),"%"), adj=0)

    # Legend
    text(rep(250,2), c(.7,.3), c("clustering", "overdispersion"), col=c("blue","red"), cex=1.5)
    text(1.3, .5, "random", col="grey40", cex=1.5, adj=.5, srt=90)
}
dev.off()

##########################
# Stats
# Run models
overdispersion <- glm(p ~ log10.eff.size + n.clade + tree.sd + factor(n.spp) + factor(n.sites), data=data[data$trt=="over",], family=quasibinomial)
clustering <- glm(p ~ log10.eff.size + n.clade + tree.sd + factor(n.spp) + factor(n.sites), data=data[data$trt=="under",], family=quasibinomial)
null <- lm(p ~ n.clade + tree.sd + factor(n.spp) + factor(n.sites), data=data[data$trt=="null",])

# Summaries for MS
xtable(summary(overdispersion))
summary(overdispersion)
xtable(summary(clustering))
summary(clustering)
xtable(summary(null))
summary(null)

# Text descriptions for MS
predict(overdispersion, data.frame(log10.eff.size=log10(10), n.clade=5, n.spp=50, n.sites=50, tree.sd=1), type="response")
predict(overdispersion, data.frame(log10.eff.size=log10(10), n.clade=5, n.spp=100, n.sites=100, tree.sd=1), type="response")

predict(null, data.frame(n.clade=10, n.spp=100, n.sites=50, tree.sd=1))
predict(null, data.frame(n.clade=5, n.spp=100, n.sites=50, tree.sd=1))
