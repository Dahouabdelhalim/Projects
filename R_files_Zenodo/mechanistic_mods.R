# This script defines a suite of useful functions that will
# be used often adapted from Sarah Orlofske & Brett
# Melbourne's code by Max Joseph 6/8/14

# Biological models for trematode Ribeiroia ondatrae
# transmission of cercariae to tadpole second intermediate
# hosts Each model contains a separate 'biological model'
# which is then combined with the stochastic component in the
# negative log likelihood.

# Models are taken from the literature and have been
# converted from forms for microparasite transmission

# All models were initially developed as differential
# equations then solved analytically (in Mathematica) to the
# form used here.

# _____Parameter
# definitions_____________________________________________________
# Np = intial number of parasites added as cercariae Beta =
# the transmission parameter H = number of hosts (Note: T was
# used in the orginial model description) t = time 30 minutes
# (end of the experiment after which tadpoles were removed
# and no other cercariae could transmit) v = container volume
# q = parameter in the Power Relationship and Host-parasitoid
# interference model c = parameter in the Asymptotic
# transmission model

# Return M, the number of metacercariae in the tadpole (those
# cercariae that have successfully transmitted)

library(deSolve)
library(gsl)

### Constant risk 1 Per-capita risk of infection for parasites
const1 <- function(par, Np = NULL, H = NULL, v = NULL, t = NULL) {
    Beta <- par[1]
    Ct <- Np * exp(-Beta * t)
    # num. metacercariae in frog = initial # cercariae - # left
    # in pool
    M <- Np - Ct
    return(M)
}
attr(const1, "data") <- c('Np', 't')

### Constant risk 2 Constant risk of infection for hosts
const2 <- function(par, Np = NULL, H = NULL, v = NULL, t = NULL) {
    Beta <- par[1]
    Ct <- Np - Beta * H * t
    stopifnot(Beta * H * t <= Np)  # can't transmit more cercariae than we have
    M <- Np - Ct
    return(M)
}
attr(const2, "data") <- c('Np', 't', 'H')

### Density dependent 1 Density dependence in both parasites
### and hosts
densdep1 <- function(par, Np = NULL, H = NULL, v = NULL, t = NULL) {
    Beta <- par[1]
    Ct <- Np * exp(-Beta * H * t / v ^ 2)
    stopifnot(Np >= Ct)
    M <- Np - Ct
    return(M)
}
attr(densdep1, "data") <- c('Np', 't', 'H', 'v')

### Density dependent 2 Density dependence in parasites only
densdep2 <- function(par, Np = NULL, H = NULL, v = NULL, t = NULL) {
    Beta <- par[1]
    Ct <- Np * exp(-Beta * H * t / v)
    M <- Np - Ct
    return(M)
}
attr(densdep2, "data") <- c('Np', 't', 'H', 'v')


### Frequency dependent 1 Dependence on numbers independent of
### density
freqdep1 <- function(par, Np = NULL, H = NULL, v = NULL, t = NULL) {
    Beta <- par[1]
    Ct <- Np * exp(-Beta * H * t)
    M <- Np - Ct
    return(M)
}
attr(freqdep1, "data") <- c('Np', 't', 'H')

### Frequency dependent 2 Dependent on ratio of parasites to
### hosts
freqdep2 <- function(par, Np = NULL, H = NULL, v = NULL, t = NULL) {
    Beta <- par[1]
    require(gsl)  # for the Lambert W function
    Ct <- H * lambert_W0((Np * exp((Np / H) - (Beta * t / v))) / H)
    M <- Np - Ct
    return(M)
}
attr(freqdep2, "data") <- c('Np', 't', 'H', 'v')

# alternative fd2 ode formulation for when lambertW function
# fails
fd2_ode <- function(par, Np = NULL, H = NULL, v = NULL, t = NULL) {
    beta <- par[1]
    require(deSolve)
    # Numerical solution for negative binomial model 1 dC/dt = -
    # beta * ((c(t) * h) / ((c(t) + h) * v))
    fdmodel <- function(t, C, parms) {
        beta <- parms[1]
        H <- parms[2]
        v <- parms[3]
        dC <- -beta * ((C * H)/((C + H) * v))
        return(list(dC))
    }
    
    # First check for correct vector lengths
    veclengths <- c(length(beta), length(Np), length(H), length(v), 
        length(t))
    maxveclength <- max(veclengths)
    if (any(!(veclengths %in% c(1, maxveclength)))) {
        stop("Vectors must be same length")
    }
    # Expand parameters to vectors if necessary
    if (length(beta) < maxveclength) 
        beta <- rep(beta, maxveclength)
    if (length(Np) < maxveclength) 
        Np <- rep(Np, maxveclength)
    if (length(H) < maxveclength) 
        H <- rep(H, maxveclength)
    if (length(v) < maxveclength) 
        v <- rep(v, maxveclength)
    if (length(t) < maxveclength) 
        t <- rep(t, maxveclength)
    # Solve
    M <- rep(NA, maxveclength)
    for (i in 1:maxveclength) {
        out <- ode(y = Np[i], times = c(0, t[i]), func = fdmodel, 
            parms = c(beta[i], H[i], v[i]), method = "lsoda")
        # convert to number of metacercariae
        M[i] <- Np[i] - out[2, 2]
    }
    return(M)
}
attr(fd2_ode, "data") <- c('Np', 't', 'H', 'v')

### Power on C
powerC <- function(par, Np = NULL, H = NULL, v = NULL, t = NULL) {
    Beta <- par[1]
    q <- par[2]
    Ct <- ((-1 + q) * (((Np ^ (1 - q)) / (-1 + q)) + Beta*H*t)) ^ (1 / (1 - q))
    M <- Np - Ct
    return(M)
}
attr(powerC, "data") <- c('Np', 't', 'H')

### Power on H
powerH <- function(par, Np = NULL, H = NULL, v = NULL, t = NULL) {
    Beta <- par[1]
    p <- par[2]
    Ct <- Np * exp(-Beta * t * H ^ p)
    M <- Np - Ct
    return(M)
}
attr(powerH, "data") <- c('Np', 't', 'H')

### Power on C and H
powerCH <- function(par, Np = NULL, H = NULL, v = NULL, t = NULL) {
    Beta <- par[1]
    q <- par[2]
    p <- par[3]
    Ct <- ((-1 + q) * (((Np^(1 - q))/(-1 + q)) + Beta * (H^p) * t))^(1/(1 - q))
    M <- Np - Ct
    return(M)
}
attr(powerCH, "data") <- c('Np', 't', 'H')

### Negative binomial 1
NB1 <- function(par, Np = NULL, H = NULL, v = NULL, t = NULL) {
    require(deSolve)
    Beta <- par[1]
    k <- par[2]
    NBmodel <- function(t, C, parms) {
        beta <- parms[1]
        k <- parms[2]
        dC <- -k * log(1 + beta * C/k)
        return(list(dC))
    }
    
    # First check for correct vector lengths
    veclengths <- c(length(Beta), length(Np), length(k), length(t))
    maxveclength <- max(veclengths)
    if (any(!(veclengths %in% c(1, maxveclength)))) 
        stop("Vectors must be same length")
    # Expand parameters to vectors if necessary
    if (length(Beta) < maxveclength) 
        Beta <- rep(Beta, maxveclength)
    if (length(Np) < maxveclength) 
        Np <- rep(Np, maxveclength)
    if (length(k) < maxveclength) 
        k <- rep(k, maxveclength)
    if (length(t) < maxveclength) 
        t <- rep(t, maxveclength)
    # Solve
    M <- rep(NA, maxveclength)
    for (i in 1:maxveclength) {
        out <- ode(y = Np[i], times = c(0, t[i]), func = NBmodel, 
            parms = c(Beta[i], k[i]), method = "radau")
        # convert to number of metacercariae
        M[i] <- Np[i] - out[2, 2]
    }
    return(M)
}
attr(NB1, "data") <- c('Np', 't')

### Negative binomial 2
NB2 <- function(par, Np = NULL, H = NULL, v = NULL, t = NULL) {
    require(deSolve)
    Beta <- par[1]
    k <- par[2]
    NBmodel <- function(t, C, parms) {
        beta <- parms[1]
        k <- parms[2]
        H <- parms[3]
        dC <- -k * H * log(1 + beta * C/k)
        return(list(dC))
    }
    
    # First check for correct vector lengths
    veclengths <- c(length(Beta), length(Np), length(k), length(H), 
        length(t))
    maxveclength <- max(veclengths)
    if (any(!(veclengths %in% c(1, maxveclength)))) 
        stop("Vectors must be same length")
    # Expand parameters to vectors if necessary
    if (length(Beta) < maxveclength) 
        Beta <- rep(Beta, maxveclength)
    if (length(Np) < maxveclength) 
        Np <- rep(Np, maxveclength)
    if (length(k) < maxveclength) 
        k <- rep(k, maxveclength)
    if (length(H) < maxveclength) 
        H <- rep(H, maxveclength)
    if (length(t) < maxveclength) 
        t <- rep(t, maxveclength)
    # Solve
    M <- rep(NA, maxveclength)
    for (i in 1:maxveclength) {
        out <- ode(y = Np[i], times = c(0, t[i]), func = NBmodel, 
            parms = c(Beta[i], k[i], H[i]), method = "radau")
        # convert to number of metacercariae
        M[i] <- Np[i] - out[2, 2]
    }
    return(M)
} 
attr(NB2, "data") <- c('Np', 't', 'H')
