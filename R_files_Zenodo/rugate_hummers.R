# Rugated thin-film optical modeling
# Capabiltiites:
# - wavelength-dependent complex refractive index

# Default parameters for melanin from Stavenga BoP 2015 study
refrind <- function(A=1.648,B=23700,C=0.56,lambda=300:700,lambda_i=270) {
  cbind(wl=lambda,
        nker = 1.532 + 5890/lambda^2,  # bird keratin,
        nmel = A + B/lambda^2,
        kmel = C * exp(-lambda/lambda_i))
}
# plot(nmel~wl, data=refrind(A=1.648,B=23700,C=0.56,lambda_i=270), type='l', ylim=c(1.54, 2))
# lines(nmel~wl, data=refrind(A=1.548,B=23700,C=0.56,lambda_i=270), type='l', lty=2)
# lines(nmel~wl, data=refrind(A=1.348,B=23700,C=0.56,lambda_i=270), type='l', lty=3)


# HR(d=100,)
# ps = platelet spacing
# pt = platelet thickness
# phi = air space density
# air = air space diameter
# cor = cortex thickness

# Wrapper function to run analyses
# For birds (Stavenga et al. 2015): A = 1.648, B = 23700, C (a_m in their paper) = 0.56, lambda_i (b_m in their paper) = 270 nm
# For beetles: A=1.56, B=36000, C=1.62, lambda_i=142
HR <- function(pt, air, ps=NULL, cor=0, ptsmall=0, smallair, phi=1, layers, ag0=0, N=500, lim=c(300, 700), nwl=100, A=1.648, B=23700, C=0.56, lambda_i=270, nmel=NULL, nker=NULL, kmel=NULL, nend=1.56-0*1i, shape="sphere", ...) {
  lambda <- seq(lim[1], lim[2], length=nwl)
  pars <- cbind(nker = 1.532 + 5890/lambda^2,  # bird keratin,
                nmel = A + B/lambda^2,
                kmel = C * exp(-lambda/lambda_i))
  #
  if (!(is.null(nmel) & is.null(nker))) {
    pars <- cbind(nker = rep(nker, nwl), nmel=rep(nmel, nwl), kmel=rep(kmel, nwl))
  }
  # Make gradients for all wavelengths
  x <- lapply(1:nrow(pars), function(x) {gen_gradient(pt=pt, ptsmall=ptsmall, smallair=smallair, ps=ps, air=air, layers=layers, cor=cor, kmel=pars[x, 'kmel'], nmel=pars[x, 'nmel'], phi=phi, nker=pars[x, 'nker'], N=N, shape=shape)})
  N <- sapply(x, "[[", "n")  # RI as function of position and wl
  K <- sapply(x, "[[", "k")  # k as a function of position and wl
  D <- x[[1]]$d  # layer thickness
  rugated(N=N, K=K, D=D, ag0=ag0, lim=lim, nwl=nwl, nend=nend, ...)
}

# Convert n/k into dielectric variables (epsilons)
eps <- function(n, k) {
  ereal <- n^2 - k^2
  eimag <- 2 * n * k
  eps <- complex(real = ereal, imaginary = -eimag)
  eps
}

# Volume average dielectic constant mixing rule
navg2 <- function(v2, n1, n2, k1, k2) {
  v1 <- 1 - v2
  e1 <- eps(n1,k1)
  e2 <- eps(n2,k2)
  eavg <- e1*v1 + e2*v2
  navg <- sqrt((sqrt(Re(eavg)^2 + Im(eavg)^2) + Re(eavg)) / 2)
  kavg <- sqrt((sqrt(Re(eavg)^2 + Im(eavg)^2) - Re(eavg)) / 2)
  complex(real=navg, imaginary=-kavg)
}

# Refractive index function for air spheres in melanin platelets
# ps = platelet spacing
# pt = platelet thickness
# phi = air space density
# air = air space diameter
# cor = cortex thickness

# gen_gradient(plot=T, layers=10)

# gen_gradient(cor=100,air=100,smallair=0,ptsmall=0,pt=100,plot=T)

# cor=100
# air=150
# smallair=0
# ptsmall=0
# pt=100

gen_gradient <- function(cor=100, air=150, smallair=10, ptsmall=50, pt=200, ps=NULL, phi=1, nmel=2, kker=0, kmel=0.1, nker=1.56, nair=1, N=100, layers=3, shape=c("sphere", "block"), plot=FALSE) {
  shape <- match.arg(shape)
  if (is.null(ps)) {
    ps <- pt
  }
   # stack of platelets
  # melanin at top of first platelet
  z1 <- c(0, (pt - air)/2)
  n1 <- rep(nmel, 2)
  k1 <- rep(-kmel, 2)

  # air sphere layer
  if (shape=="block") {
    z2 <- c(0, air)
    n2 <- rep(nair, 2)
    k2 <- rep(0, 2)  
  }
  if (shape=="sphere") {
    a <- air/2
    z2 <- seq(-air, 0, length=N)
    n2 <- nmel + phi * (nair-nmel) * (1 - ((z2/a) + 1)^2)
    z2 <- rev(abs(z2))
    vf_mel <- (n2-nair)/(nmel-nair)
    n_avg <- sapply(vf_mel, FUN=navg2, n1=nair, n2=nmel, k1=0, k2=kmel)
    n2 <- Re(n_avg)
    k2 <- Im(n_avg)
  }
  # bottom melanin layer
  z3 <- c(0, ((pt-air)/2))
  n3 <- rep(nmel, 2)
  k3 <- rep(-kmel, 2)
  # combine all 3 parts of a melanosomes now:
  zstruc <- c(z1, z2 + max(z1), z3+max(z1)+max(z2))
  nstruc <- c(n1, n2, n3)
  kstruc <- c(k1, k2, k3)
  # keratin between platelets
  if ((ps-pt)!=0) {
    zstruc <- c(zstruc, c(max(zstruc), ps - pt + max(zstruc)))
    nstruc <- c(nstruc, rep(nker, 2))
    kstruc <- c(kstruc, rep(-kker, 2))
  }
  # make layers
  zstruc <- as.numeric(t(t(replicate(layers, zstruc)) + (((1:layers)*ps)-ps) ) )
  nstruc <- as.numeric(replicate(layers, nstruc))
  kstruc <- as.numeric(replicate(layers, kstruc))
  # small platelets
  if (ptsmall!=0) {
    # melanin at top of platelet
    z1 <- c(0, (ptsmall-smallair)/2)
    n1 <- rep(nmel, 2)
    k1 <- rep(-kmel, 2)
    # air sphere layer
    if (shape=="block") {
      z2 <- c(0, smallair)
      n2 <- rep(nair, 2)
      k2 <- rep(0, 2)
    }
    if (shape=="sphere") {
      a <- smallair/2
      z2 <- seq(-smallair, 0, length=N)
      n2 <- nmel + phi * (nair-nmel) * (1 - ((z2/a) + 1)^2)
      z2 <- rev(abs(z2))
      vf_mel <- (n2-nair)/(nmel-nair)
      n_avg <- sapply(vf_mel, FUN=navg2, n1=nair, n2=nmel, k1=0, k2=kmel)
      n2 <- Re(n_avg)
      k2 <- Im(n_avg)
    }
    # bottom melanin layer
    z3 <- c(0, ((ptsmall-smallair)/2))
    n3 <- rep(nmel, 2)
    k3 <- rep(-kmel, 2)
    z <- c(z1, z2+max(z1), z3+max(z1)+max(z2))
    n <- c(n1,n2,n3)
    k <- c(k1,k2,k3)
    # keratin spacing between
    if (ps != pt) {
      z4 <- c(0, ps-pt)
      n4 <- rep(nker, 2)
      k4 <- rep(kker, 2)
      z <- c(z, z4 + max(z))
      n <- c(n, n4)
      k <- c(k, k4)  
    }
    zstruc <- c(z, zstruc + max(z))
    nstruc <- c(n, nstruc)
    kstruc <- c(k, kstruc)
  }
  # add cortex
  if (cor!=0) {
    z <- c(0, cor)
    n <- rep(nker, 2)
    k <- rep(-kker, 2)
    zstruc <- c(z, zstruc + max(z))
    nstruc <- c(n, nstruc)
    kstruc <- c(k, kstruc)
  }
  res <- data.frame(pos = zstruc, n = nstruc, k = kstruc)
  if (plot) {
    par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(1.5,.5,0),ps=8,mex=.75)
    plot(n~pos,data=res,type='l', ylim=range(nair,nmel,nker))
    plot(k~pos,data=res,type='l', ylim=-range(kmel,kker))
  }
  res <- res[complete.cases(res), ]
  res <- res[!duplicated(res), ]
  res <- res[diff(res$pos)!=0, ]
  # calculate layer thicknesses
  newn <- res$n[-nrow(res)]
  newd <- diff(res$pos)
  newk <- res$k[-nrow(res)]
  res <- data.frame(d=newd, n=newn, k=newk)
  invisible(res)
}

# gen_gradient(cor=100, ps=120, pt=100, air=50, N=100, phi=1, plot=T, layers=2)
# gen_gradient(cor=200, ptsmall=0, ps=250, pt=150, air=0, N=100, plot=T)
# gen_gradient(cor=100, pt=100, air=50, N=100, phi=1, plot=T, layers=10)
# gen_gradient(cor=0, pt=100, air=50, N=100, phi=1, plot=T, layers=10)
# gen_gradient(cor=100, ptsmall=50, smallair=10, pt=200, air=50, N=100, phi=1, plot=T, layers=5)


# Function to perform transfer matrix calculations

rugated <- function(N, K, D, ag0=0, n0=1-0*1i, nend=1-0*1i, lim=c(300, 700), nwl=100, plot=TRUE) {
  ms <- NULL
  mp <- NULL
  wavel <- seq(lim[1], lim[2], length=nwl)
  ag0 <- ag0/180*pi
  for (i in seq_along(wavel)) {
    n <- complex(real = N[, i], imaginary = K[, i])
    #air - first layer
    agstart <- asin(n0*sin(ag0)/n[1])
    r1p <- (n0*cos(ag0)-n[1]*cos(agstart)) / (n0*cos(ag0)+n[1]*cos(agstart))
    I1p <- matrix(c(1, r1p, r1p, 1), nrow=2, ncol=2)
    r1s <- (n[1]*cos(ag0)-n0*cos(agstart)) / (n[1]*cos(ag0)+n0*cos(agstart))
    I1s <- matrix(c(1, r1s, r1s, 1), nrow=2, ncol=2)
    # b1 <- (2*pi*(d/length(n))*n[1]*cos(agstart)) / wavel[i]
    b1 <- (2*pi*D[1]*n[1]*cos(agstart)) / wavel[i]
    L1 <- matrix(c(exp(b1*1i), 0, 0, exp(b1*-1i)), nrow=2, ncol=2)
    SCAT1s <- I1s %*% L1
    SCAT1p <- I1p %*% L1
    #################################
    #START ITERATION FOR VARIABLE RI#
    #################################
      for (j in 2:length(n)){
      #COMPLEX ANGLE OF INCIDENCE (SNELL'S LAW)
      ag1 <- asin(n0*sin(ag0)/n[j-1]) 	#previous layer
      ag2 <- asin(n0*sin(ag0)/n[j])	#current layer

      #INTERFACE MATRICES

      #previous layer - current layer
      r2p <- (n[j-1] * cos(ag1) - n[j] * cos(ag2)) / (n[j-1] * cos(ag1) + n[j] * cos(ag2))
      I2p <- matrix(c(1,r2p,r2p,1), nrow=2, ncol=2)
      r2s <- (n[j] * cos(ag1)-n[j-1] * cos(ag2))/(n[j] * cos(ag1)+n[j-1] * cos(ag2))
      I2s <- matrix(c(1,r2s,r2s,1),nrow=2,ncol=2)

      #TRANSFER MATRICES

      #current layer
      # b2 <- (2*pi*(d/length(n))*n[j]*cos(ag2)) / wavel[i]
      b2 <- (2*pi*D[j]*n[j]*cos(ag2)) / wavel[i]
      L2 <- matrix(c(exp(b2*1i), 0, 0, exp(b2*-1i)), nrow=2, ncol=2)

      SCAT1s <- SCAT1s %*% I2s %*% L2
      SCAT1p <- SCAT1p %*% I2p %*% L2

      }

      ###############################
      #TESTING FOR MULTIPLE PLATELETS
      ###############################
      #Final passage to substrate (assuming same as first interface)
      aglast <- asin(n0*sin(ag0)/n[length(n)])
      agfinish <- asin(n0*sin(ag0)/nend)
      rFp <- (n[length(n)]*cos(aglast)-nend*cos(agfinish))/(n[length(n)]*cos(aglast)+nend*cos(agfinish))
      IFp <- matrix(c(1,rFp,rFp,1),nrow=2,ncol=2)
      rFs <- (nend*cos(aglast)-n[length(n)]*cos(agfinish))/(nend*cos(aglast)+n[length(n)]*cos(agfinish))
      IFs <- matrix(c(1, rFs, rFs, 1), nrow=2, ncol=2)
      SCAT1s <- SCAT1s %*% IFs
      SCAT1p <- SCAT1p %*% IFp
      #AMPLITUDE REFLECTIVITY
      refs <- SCAT1s[2, 1] / SCAT1s[1, 1]
      refp <- SCAT1p[2, 1] / SCAT1p[1, 1]
      #FINAL REFLECTANCE
      REFs <- refs %*% Conj(refs)
      REFp <- refp %*% Conj(refp)
      #MODEL
      ms <- c(ms, REFs * 100)
      mp <- c(mp, REFp * 100)
  }

  model.s <- as.numeric(ms)
  model.p <- as.numeric(mp)
  refl <- rowMeans(cbind(model.s, model.p))
  model <- (model.s + model.p) / 2

  #DATA TABLE
  wavelength <- wavel
  mfilm <- data.frame(wavelength, model.s, model.p, refl, model)

  #GRAPHIC AND RESULT
  if (plot) {
    plot(wavelength, model.s, xlab="Wavelength", ylab="Reflectance (%)", lwd=2, 
      col="red", type="l", ylim=c(0,100))
    #points(wavelength,smooth.spline(model.p,df=15)$y,typ="l", lwd=2, col="blue")
    points(wavelength, model.p, lwd=2, col="blue", type="l")
    points(wavelength, model, type="l", lwd=2, lty=3)
  }

  invisible(mfilm)
}
