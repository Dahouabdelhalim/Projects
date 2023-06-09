# This script defines a suite of useful functions that will
# be used often adapted from Sarah Orlofske & Brett
# Melbourne's code by Max Joseph 6/8/14

############################################################################################### Negative log likelihood functions for evaluating the
############################################################################################### biological models above with data and incorporating
############################################################################################### stochasticity in the distribution generally, parameters to
############################################################################################### be estimated are packed into the object p fixed parameters
############################################################################################### are named and need not be unrolled our expected values
############################################################################################### (meta) are calculated with the mechanistic model observed
############################################################################################### numbers of metacercariae are called metacerc

# Const1 Function
NLL.const1 <- function(p, Np, H = NULL, v = NULL, t = NULL, metacerc = NULL) {
    meta <- const1(par = p, Np = Np, H = H, v = v, t = t)
    nll <- -sum(dbinom(metacerc, prob = meta/Np, size = Np, log = TRUE))
    return(nll)
}

# Const2 Function
NLL.const2 <- function(p, Np, H = NULL, v = NULL, t = NULL, metacerc = NULL) {
    meta <- const2(par = p, Np = Np, H = H, v = v, t = t)
    nll <- -sum(dbinom(metacerc, prob = meta/Np, size = Np, log = TRUE))
    return(nll)
}

# Density dependent 1 Function
NLL.den1 <- function(p, Np, H = NULL, v = NULL, t = NULL, metacerc = NULL) {
    meta <- densdep1(par = p, Np = Np, H = H, v = v, t = t)
    nll <- -sum(dbinom(metacerc, prob = meta/Np, size = Np, log = TRUE))
    return(nll)
}

# Density dependent 2 Function
NLL.den2 <- function(p, Np, H = NULL, v = NULL, t = NULL, metacerc = NULL) {
    meta <- densdep2(par = p, Np = Np, H = H, v = v, t = t)
    nll <- -sum(dbinom(metacerc, prob = meta/Np, size = Np, log = TRUE))
    return(nll)
}

# Frequency dependent 1 Function
NLL.freq1 <- function(p, Np, H = NULL, v = NULL, t = NULL, metacerc = NULL) {
    meta <- freqdep1(par = p, Np = Np, H = H, v = v, t = t)
    nll <- -sum(dbinom(metacerc, prob = meta/Np, size = Np, log = TRUE))
    return(nll)
}

# Frequency dependent 2 Function
NLL.freq2 <- function(p, Np, H = NULL, v = NULL, t = NULL, metacerc = NULL) {
    meta <- freqdep2(par = p, Np = Np, H = H, v = v, t = t)
    nll <- -sum(dbinom(metacerc, prob = meta/Np, size = Np, log = TRUE))
    return(nll)
}

NLL.freq2_ode <- function(p, Np, H = NULL, v = NULL, t = NULL, metacerc = NULL) {
    meta <- fd2_ode(par = p, Np = Np, H = H, v = v, t = t)
    nll <- -sum(dbinom(metacerc, prob = meta/Np, size = Np, log = TRUE))
    return(nll)
}

# Power C Relationship
NLL.powC <- function(p, Np, H = NULL, v = NULL, t = NULL, metacerc = NULL) {
    meta <- powerC(par = p, Np = Np, H = H, v = v, t = t)
    nll <- -sum(dbinom(metacerc, prob = meta/Np, size = Np, log = TRUE))
    return(nll)
}

# Power H Relationship
NLL.powH <- function(p, Np, H = NULL, v = NULL, t = NULL, metacerc = NULL) {
    meta <- powerH(par = p, Np = Np, H = H, v = v, t = t)
    nll <- -sum(dbinom(metacerc, prob = meta/Np, size = Np, log = TRUE))
    return(nll)
}

# Power CH Relationship
NLL.powCH <- function(p, Np, H = NULL, v = NULL, t = NULL, metacerc = NULL) {
    meta <- powerCH(par = p, Np = Np, H = H, v = v, t = t)
    nll <- -sum(dbinom(metacerc, prob = meta/Np, size = Np, log = TRUE))
    return(nll)
}

# Negative binomial model 1
NLL.NB1 <- function(p, Np, H = NULL, v = NULL, t = NULL, metacerc = NULL) {
    meta <- NB1(par = p, Np = Np, H = H, v = v, t = t)
    nll <- -sum(dbinom(metacerc, prob = meta/Np, size = Np, log = TRUE))
    return(nll)
}

# Negative binomial model 2
NLL.NB2 <- function(p, Np, H = NULL, v = NULL, t = NULL, metacerc = NULL) {
    meta <- NB2(par = p, Np = Np, H = H, v = v, t = t)
    nll <- -sum(dbinom(metacerc, prob = meta/Np, size = Np, log = TRUE))
    return(nll)
} 
