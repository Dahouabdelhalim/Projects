############################
#HEADERS####################
############################
library(caper)
library(picante)
library(pez)
library(parallel)
library(phylolm)
library(phytools)
library(weights)

############################
# DISPERSION################
############################
dispersion <- readRDS("dispersion.RDS")

estimates <- sapply(dispersion, function(x) tryCatch(c(x$migration$DEstimate, x$breeding$DEstimate, x$refuge$DEstimate, x$tracking$DEstimate), error=function(z) rep(NA,4)))
p.0.1 <- apply(sapply(dispersion, function(x) tryCatch(c(x$migration$Pval0, x$breeding$Pval0, x$refuge$Pval0, x$tracking$Pval0), error=function(z) rep(NA,4))), 1, function(x) sum(x <= 0.01, na.rm=TRUE))
p.1.1 <- apply(sapply(dispersion, function(x) tryCatch(c(x$migration$Pval1, x$breeding$Pval1, x$refuge$Pval1, x$tracking$Pval1), error=function(z) rep(NA,4))), 1, function(x) sum(x <= 0.01, na.rm=TRUE))
p.0.5 <- apply(sapply(dispersion, function(x) tryCatch(c(x$migration$Pval0, x$breeding$Pval0, x$refuge$Pval0, x$tracking$Pval0), error=function(z) rep(NA,4))), 1, function(x) sum(x <= 0.05, na.rm=TRUE))
p.1.5 <- apply(sapply(dispersion, function(x) tryCatch(c(x$migration$Pval1, x$breeding$Pval1, x$refuge$Pval1, x$tracking$Pval1), error=function(z) rep(NA,4))), 1, function(x) sum(x <= 0.05, na.rm=TRUE))
coefs <- matrix(
    c(apply(estimates, 1, mean, na.rm=TRUE), apply(estimates, 1, sd, na.rm=TRUE), apply(estimates, 1, function(x) sum(!is.na(x)))),
    ncol=4, nrow=3, dimnames=list(c("mean","sd","count"), c("migration","breeding","refuge","tracking")),
    byrow=TRUE
)
coefs <- rbind(coefs, 100*p.0.5/coefs["count",], 100*p.0.1/coefs["count",], 100*p.1.5/coefs["count",], 100*p.1.1/coefs["count",])
rownames(coefs)[4:7] <- c("% P(D=0)<=0.05", "% P(D=0)<=0.01", "% P(D=1)<=0.05", "% P(D=1)<=0.01")

############################
#CONTINGENCY TABLE##########
############################
#Calculations
contingency.wrapper <- function(x, n.null=10){
    comm <- as.matrix(sample2matrix(data.frame(type=c(x$data$migration.type,x$data$movement), abund=rep(1, nrow(x$data)), species=rep(rownames(x$data),2))))
    x <- comparative.comm(x$phy, comm)
    phylosor <- as.matrix(.phylosor(x))
    null <- array(NA, c(dim(phylosor),n.null))
    for(i in seq_len(n.null)){
        x$comm <- randomizeMatrix(x$comm, "trialswap", prod(dim(x$comm),10))
        null[,,i] <- as.matrix(.phylosor(x))
    }
    upper <- matrix(rowSums(apply(null, 3, function(x) phylosor > x)), dim(phylosor))
    lower <- matrix(rowSums(apply(null, 3, function(x) phylosor < x)), dim(phylosor))
    return(list(obs=phylosor, lower=lower, upper=upper))
}
saveRDS(mcMap(contingency.wrapper, readRDS("basic_workspace.RDS"), mc.cores=24), "contingency.RDS")
#contingency <- readRDS("certain.contingency.RDS")

#Post-processing
contingency.summary <- function(sims){
    mean <- apply(array(unlist(lapply(sims, function(x) x$obs)), dim=c(dim(sims[[1]]$obs),length(sims))), 1:2, mean)
    sd <- apply(array(unlist(lapply(sims, function(x) x$obs)), dim=c(dim(sims[[1]]$obs),length(sims))), 1:2, sd)
    lower <- apply(array(unlist(lapply(sims, function(x) x$lower)), dim=c(dim(sims[[1]]$lower),length(sims))), 1:2, function(y) sum(y)/(length(y)*100))
    upper <- apply(array(unlist(lapply(sims, function(x) x$upper)), dim=c(dim(sims[[1]]$upper),length(sims))), 1:2, function(y) sum(y)/(length(y)*100))
    rownames(mean) <- colnames(mean) <- rownames(lower) <- colnames(lower) <- rownames(upper) <- colnames(upper) <- rownames(sims[[1]]$obs)
    return(list(mean=mean, sd=sd, lower.p=lower, upper.p=upper))
}
contingency <- contingency.summary(contingency)
write.csv(contingency, "neat/certain.contingency.csv")

############################
#REGRESSIONS################
############################
regression.summary <- function(sims, p.val=0.05){
    .summary <- function(sims, which, p.val){
        coef <- sapply(sims, function(x) tryCatch(coef(x[[which]]), error=function(y) rep(NA,length(coef(sims[[1]][[which]])))))
        se <- sapply(sims, function(x) tryCatch(summary(x[[which]])$coefficients[,2], error=function(y) rep(NA,length(coef(sims[[1]][[which]])))))
        p <- sapply(sims, function(x) tryCatch(summary(x[[which]])$coefficients[,4], error=function(y) rep(NA,length(coef(sims[[1]][[which]])))))
        output <- matrix(ncol=nrow(coef), nrow=4, dimnames=list(c("estimate","SE","p","total"),rownames(coef)))
        subset <- !is.na(coef[1,])
        output[4,] <- rep(sum(subset), ncol(output))
        for(i in seq_len(ncol(output))){
            output[1,i] <- paste(round(mean(coef[i,subset]),2), "±", round(sd(coef[i,subset]),4))
            output[2,i] <- paste(round(mean(se[i,subset]),2), "±", round(sd(se[i,subset]),4))
            output[3,i] <- paste(round(mean(p[i,subset]),4), "±", round(sd(p[i,subset]),4))
        }
        return(output)
    }

    return(list(pglm=.summary(sims,1,p.val), pgls=.summary(sims,2,p.val)))
}

data <- readRDS("regression_workspace.RDS")

regressions <- readRDS("regressions.RDS")
t <- regression.summary(regressions)
alphas <- sapply(regressions, function(x) tryCatch(x[[1]]$alpha, error=function(e) NA))
c(mean(alphas, na.rm=TRUE), sd(alphas, na.rm=TRUE))
write.csv(t$pglm, "~/Dropbox/2013 Mammal migration/analyses/certain.regressions.pglm.csv")
write.csv(t$pgls, "~/Dropbox/2013 Mammal migration/analyses/certain.regressions.pgls.csv")

regressions <- readRDS("regressions_treatened.RDS")
t <- regression.summary(regressions)
alphas <- sapply(regressions, function(x) tryCatch(x[[1]]$alpha, error=function(e) NA))
c(mean(alphas, na.rm=TRUE), sd(alphas, na.rm=TRUE))
write.csv(t$pglm, "~/Dropbox/2013 Mammal migration/analyses/certain.threatened.regressions.pglm.csv")
write.csv(t$pgls, "~/Dropbox/2013 Mammal migration/analyses/certain.threatened.regressions.pgls.csv")

subset.regression.summary <- function(sims, p.val=0.05){
    .summary <- function(sims, which, p.val){
        coef <- sapply(sims, function(x) tryCatch(coef(x[[which]]), error=function(y) rep(NA,length(coef(sims[[1]][[which]])))))
        se <- sapply(sims, function(x) tryCatch(summary(x[[which]])$coefficients[,2], error=function(y) rep(NA,length(coef(sims[[1]][[which]])))))
        p <- sapply(sims, function(x) tryCatch(summary(x[[which]])$coefficients[,4], error=function(y) rep(NA,length(coef(sims[[1]][[which]])))))
        output <- matrix(ncol=nrow(coef), nrow=4, dimnames=list(c("estimate","SE","p","total"),rownames(coef)))
        subset <- !is.na(coef[1,])
        output[4,] <- rep(sum(subset), ncol(output))
        for(i in seq_len(ncol(output))){
            output[1,i] <- paste(round(mean(coef[i,subset]),2), "±", round(sd(coef[i,subset]),4))
            output[2,i] <- paste(round(mean(se[i,subset]),2), "±", round(sd(se[i,subset]),4))
            output[3,i] <- paste(round(mean(p[i,subset]),4), "±", round(sd(p[i,subset]),4))
        }
        return(output)
    }

    return(list(r=.summary(sims,1,p.val), b=.summary(sims,2,p.val), t=.summary(sims,3,p.val)))
}

subset.regression <- readRDS("subset.threatened.regressions.RDS")
subset.regression.pglm <- subset.regression.summary(lapply(subset.regression, function(x) x[[1]]))
r.alphas <- sapply(subset.regression, function(x) tryCatch(x[[1]][[1]]$alpha, error=function(e) NA))
c(mean(r.alphas, na.rm=TRUE), sd(r.alphas, na.rm=TRUE))
b.alphas <- sapply(subset.regression, function(x) tryCatch(x[[1]][[2]]$alpha, error=function(e) NA))
c(mean(b.alphas, na.rm=TRUE), sd(b.alphas, na.rm=TRUE))
t.alphas <- sapply(subset.regression, function(x) tryCatch(x[[1]][[3]]$alpha, error=function(e) NA))
c(mean(t.alphas, na.rm=TRUE), sd(t.alphas, na.rm=TRUE))
write.csv(subset.regression.pglm, "subset.threatened.regressions.pglm.csv")


############################
#RECONSTRUCTIONS############
############################
saveRDS(mcMap(function(x) make.simmap(x$phy, setNames(x$data$migration, rownames(x$data)), model="ER", nsim=1000, pi="estimated"), readRDS("basic_workspace.RDS"), mc.cores=24), "reconstruction.RDS")
