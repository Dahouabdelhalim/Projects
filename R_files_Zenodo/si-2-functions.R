# Function to calculate variance among clades

#Headers
require(vegan)
require(pez)
require(caper)
require(geiger)
require(MASS)

#Basic function
clade.var <- function(data, alpha=0.05, clade=NULL){
    #Setup
    ci <- function(x, alpha=0.05){
        upper <- ((length(x)-1)*var(x)) / qchisq(alpha/2, length(x)-1)
        lower <- ((length(x)-1)*var(x)) / qchisq(1-(alpha/2), length(x)-1)
        return(cbind(upper, lower))
    }
    if(!inherits(data, "comparative.comm")) stop("Must use pez::comparative.comm object")
    clade.mat <- clade.matrix(data$phy)$clade.matrix

    variance <- pd <- n.spp <- numeric(nrow(clade.mat))
    cis <- matrix(NA, nrow=nrow(clade.mat), ncol=2)
    for(i in seq_len(nrow(clade.mat))){
        spp <- unname(clade.mat[i,]==1)
        n.spp[i] <- sum(spp)
        variance[i] <- var(rowSums(as.matrix(comm(data)[,spp])))
        cis[i,] <- ci(rowSums(as.matrix(comm(data)[,spp])), alpha)
        if(n.spp[i] == 1) pd[i] <- 0 else pd[i] <- sum(drop_tip(phy(data),setdiff(species(data),species(data)[spp]))$edge.length)/2
    }
    expect.var <- apply(clade.mat, 1, function(x) sum(variance[which(unname(x)==1)]))
    output <- data.frame(cbind(variance, pd, n.spp, expect.var, cis))
    names(output) <- c("variance", "pd", "n.spp", "expect.var", "upper.ci", "lower.ci")
    if(!is.null(clade))
        return(output[clade,])
    return(output)
}
