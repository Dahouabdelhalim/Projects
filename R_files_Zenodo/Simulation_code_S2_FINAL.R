
sc <- function(x1, x2, precision = 10^-2) {
  flag <- F
  if (dim(x1)[1] == dim(x2)[1]) {
    n <- dim(x1)[1]
    flag <- T
    i <- 0
    while (i < n & flag == T) {
      flag <- F
      i <- i + 1
      for (j in 1:n) {
        if (mean((x1[i, ] - x2[j, ])^2) < precision^2) 
          flag <- T
      }
      
    }
  }
  return(flag)
}


etaN <- 0.1
etaP <- -0.1


n <- 10
p <- 10

##Number of traits
q <- 3

##Precision for distinguishing between similar species using sc function described above
precision <- .05



maxtime <- 3 * 10^6
starttime <- 100


int <- .01*maxtime

# parameters governing evolutionary rates (genetic variance and )

  sig2N <- 0.1
  sdV <- 0.1


  sig2P <- 0.1
  sdU <- 0.1


# parameter set removing r from P equations
r <- 0.05
a <- 2
f <- 0.1
b <- 0.05
cc <- 1
m <- 0.01
g <- 0.01

initN <- .4
initP <- .1
initV <- .1
initU <- 1

C <- matrix(etaN, nrow = q, ncol = q)
diag(C) <- 1
D <- matrix(-etaP, nrow = q, ncol = q)
diag(D) <- 1

N <- matrix(initN, nrow = n, ncol = 1)
P <- matrix(initP, nrow = p, ncol = 1)
V <- matrix(initV, nrow = n, ncol = q)
U <- matrix(initU, nrow = p, ncol = q)

Nlist <- array(0, c(n, maxtime))
Plist <- array(0, c(p, maxtime))
Vlist <- array(0, c(n, q, maxtime))
Ulist <- array(0, c(p, q, maxtime))

##Running model without trait perturbation to obtain equilibrium density and traits
for (time in 1:starttime) {
  
  A <- a + f * diag(V %*% C %*% t(V))
  B <- b * exp(-(V^2) %*% t(U^2))
  M <- m + g/diag(U %*% D %*% t(U))
  
  Nt <- N * exp(r * (1 - A * sum(N) - B %*% P))
  Pt <- P * exp(cc * r*t(B) %*% N - M)
  Vt <- V + sig2N * (-2 * r * f * sum(N) * V %*% C + 2 * r * b * V * exp(-(V^2) %*% t(U^2)) %*% (U^2 * 
                                                                                                   array(P, c(p, q))))
  Ut <- U + sig2P * (-2 *r* b * cc * U * exp(-(U^2) %*% t(V^2)) %*% (V^2 * array(N, c(n, q))) + 2 * g * 
                       array(1/diag(U %*% D %*% t(U))^2, c(p, q)) * (U %*% D))
  
  N <- Nt
  P <- Pt
  V <- abs(Vt)
  U <- abs(Ut)
}



Nlist[, 1] <- N
Plist[, 1] <- P
Vlist[, , 1] <- V
Ulist[, , 1] <- U

for (time in 2:maxtime) {

# Random trait perturbation
  
  if(time == maxtime/3){
    V <- matrix(exp(rnorm(n = n * q, mean = 0, sd = sdV)), nrow = n, ncol = q) * V
    U <- matrix(exp(rnorm(n = p * q, mean = 0, sd = sdU)), nrow = p, ncol = q) * U
  }
  
  A <- a + f * diag(V %*% C %*% t(V))
  B <- b * exp(-(V^2) %*% t(U^2))
  M <- m + g/diag(U %*% D %*% t(U))
  
  Nt <- N * exp(r * (1 - A * sum(N) - B %*% P))
  Pt <- P * exp(cc * t(B) %*% N - M)
  Vt <- V + sig2N * (-2 * r * f * sum(N) * V %*% C + 2 * r * b * V * exp(-(V^2) %*% t(U^2)) %*% (U^2 * 
                                                                                                   array(P, c(p, q))))
  Ut <- U + sig2P * (-2 * b * cc * U * exp(-(U^2) %*% t(V^2)) %*% (V^2 * array(N, c(n, q))) + 2 * g * 
                       array(1/diag(U %*% D %*% t(U))^2, c(p, q)) * (U %*% D))
  
  N <- Nt
  P <- Pt
  V <- abs(Vt)
  U <- abs(Ut)
  

  Nlist[, time] <- N
  Plist[, time] <- P
  Vlist[, , time] <- V
  Ulist[, , time] <- U
}


  par(mfcol = c(2, 1 + q))
  
  matplot((1:100)/10,t(Nlist[, int * (1 + 0:((dim(Nlist)[2]/int - 1)))]), type = "l", ylab = "Resource density", lwd=2, xlab = '')
  matplot((1:100)/10,t(Plist[, int * (1 + 0:((dim(Nlist)[2]/int - 1)))]), type = "l", ylab = "Consumer density", xlab = expression(paste('Time (3 x ', 10^5, ')')), lwd=2)
  
  for (i in 1:q) {
    matplot((1:100)/10, t(abs(Vlist[, i, int * (1 + 0:((dim(Nlist)[2]/int - 1)))])), type = "l", ylab = paste("Resource trait", i), lwd=2, xlab = '')
    matplot((1:100)/10, t(abs(Ulist[, i, int * (1 + 0:((dim(Nlist)[2]/int - 1)))])), type = "l", ylab = paste("Consumer trait", i), xlab = expression(paste('Time (2 x ', 10^5, ')')), lwd=2)
  }


# find unique species.
unique.prey <- 1:n
for(i in 1:(n - 1)) for(j in (i+1):n) {
  if(sc(matrix(V[i,], nrow=1), matrix(V[j,], nrow=1), precision = precision) == T) unique.prey[j] <- 0
}
unique.prey <- unique.prey[unique.prey > 10^(-8)]
unique.prey <- unique.prey[N[unique.prey] > 10^(-8)]

unique.pred <- 1:p
for(i in 1:(p - 1)) for(j in (i+1):p) {
  if(sc(matrix(U[i,], nrow=1), matrix(U[j,], nrow=1), precision = precision) == T) unique.pred[j] <- 0
}
unique.pred <- unique.pred[unique.pred > 10^(-8)]
unique.pred <- unique.pred[P[unique.pred] > 10^(-8)]

num.prey <- length(unique.prey)
num.pred <- length(unique.pred)

num.prey.pred <- c(num.prey, num.pred)
num.prey.pred


