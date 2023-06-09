
########################################################################################################
#                                                                                                      #
#  R Script for simulations performed in Adams and Collyer (2016) On the comparison of the strength    #
#  of morphological integration across morphometric datasets. Evolution                                #
#                                                                                                      #
########################################################################################################

# These simulations produce a figure like Fig. 1

# Simulations evaluating PLS for changing N and P
# PLS R ONLY

library(geomorph)  #install latest version from CRAN

######### 1:  Sample Size Dependency
p <- 30 
nspec <- seq(10, 3010, by = 100)
nsets <- 100
mod.gp <- gl(2, p/2)  
pls.iter <- 999
center <- geomorph:::center
sdn <- function(x) sqrt(sum((x-mean(x))^2)/length(x)) # popualtion-level standard deviation

# Simulations
pls.r <- cen.r.100 <- cen.r.1000 <- Z.1000 <- Z.100 <- array(NA, dim = c(nsets, length(nspec)))
for (i in 1:length(nspec)) {
  cat(paste("...#",i ,"#...", sep="")) # print iterations
  for (j in 1:nsets) {
    x <- matrix(rnorm((p*nspec[i])), nrow = nspec[i])
    res <- two.b.pls(x[,which(mod.gp==levels(mod.gp)[1])], 
                     x[,which(mod.gp==levels(mod.gp)[2])],
                     iter = pls.iter,
                     print.progress = FALSE)
    r.rand <- res$random.r
    pls.r[j,i] <- res$r.pls
    cen.r.100[j,i] <- res$r.pls-mean(r.rand[1:100])
    cen.r.1000[j,i] <- res$r.pls-mean(r.rand)
    Z.100[j,i] <- (res$r.pls - mean(r.rand[1:100])) / sdn(r.rand[1:100])
    Z.1000[j,i] <- (res$r.pls - mean(r.rand)) / sdn(r.rand)
    cat(paste(j,".", sep = "")) # print iterations
  }  #end j
}  #end i

# Save results
write.csv(pls.r,"plsN.csv")
write.csv(Z.100,"Z.100N.csv")
write.csv(Z.1000,"Z.1000N.csv")
write.csv(cen.r.100,"cen.r.100N.csv")
write.csv(cen.r.1000,"cen.r.1000N.csv")


######### 2: Variable Dependency
p <- seq(2, 502, by=20)
nspec <- 100
nsets <- 100

#Simulations
pls.r <-cen.r.100 <- cen.r.1000 <- Z.1000 <- Z.100 <- array(NA, dim = c(nsets,length(p)))
for (i in 1:length(p)) {
  cat(paste("...#",i ,"#...", sep="")) # print iterations
  mod.gp <- gl(2, p[i]/2)  
  for (j in 1:nsets) {
    x <- matrix(rnorm((p[i]*nspec)), nrow = nspec)
    var.x <- var(x)  # orig Cov matrix
    res <- two.b.pls(x[,which(mod.gp == levels(mod.gp)[1])],
                     x[,which(mod.gp == levels(mod.gp)[2])],
                     iter = pls.iter,
                     print.progress = FALSE)    
    r.rand <- res$random.r
    pls.r[j,i] <- res$r.pls
    cen.r.100[j,i] <- res$r.pls-mean(r.rand[1:100])
    cen.r.1000[j,i] <- res$r.pls-mean(r.rand)
    Z.100[j,i] <- (res$r.pls-mean(r.rand[1:100])) / sdn(r.rand[1:100])
    Z.1000[j,i] <- (res$r.pls-mean(r.rand)) / sdn(r.rand)
    cat(paste(j,".", sep = "")) # print iterations
  }  #end j
}  #end i

# Save results
write.csv(pls.r,"plsP.csv")
write.csv(Z.100,"Z.100P.csv")
write.csv(Z.1000,"Z.1000P.csv")
write.csv(cen.r.100,"cen.r.100P.csv")
write.csv(cen.r.1000,"cen.r.1000P.csv")

#-------------------------------- Read Previous Results, for plotting -----------------------------------

pls.r<-read.csv("plsN.csv", row.names=1, header=T)
cen.r <- read.csv("cen.r.100N.csv", row.names=1, header=T)
Z.100 <- read.csv("Z.100N.csv", row.names=1, header=T)
Z.1000 <- read.csv("Z.1000N.csv", row.names=1, header=T)

## METHOD ROBUST TO NPERM (but use more for plots)
cor(Z.100[,1], Z.1000[,1])
plot(Z.100[,1], Z.1000[,1])

r.mn.N <- apply(pls.r, 2, mean)
r.sd.N <- apply(pls.r, 2, sdn)
cen.r.mn.N <- apply(cen.r, 2, mean)
cen.r.sd.N <- apply(cen.r, 2, sdn)
z.mn.N <- apply(Z.1000, 2, mean)
z.sd.N <- apply(Z.1000, 2, sdn)

#-------------------------------- Read Previous Results, for plotting -----------------------------------

pls.r<-read.csv("plsP.csv", row.names=1, header=T)
cen.r <- read.csv("cen.r.100P.csv", row.names=1, header=T)
Z.100 <- read.csv("Z.100P.csv", row.names=1, header=T)
Z.1000 <- read.csv("Z.1000P.csv", row.names=1, header=T)

## METHOD ROBUST TO NPERM (but use more for plots)
cor(Z.100[,1], Z.1000[,1])
plot(Z.100[,1], Z.1000[,1])

r.mn.p <- apply(pls.r, 2, mean)
r.sd.p <- apply(pls.r, 2, sdn)
cen.r.mn.p <- apply(cen.r, 2, mean)
cen.r.sd.p <- apply(cen.r, 2, sdn)
z.mn.p <- apply(Z.1000, 2, mean)
z.sd.p <- apply(Z.1000, 2, sdn)

#---------------------------------- Plot the Results (as in Fig. 1) -------------------------------------

par(mfcol=c(3,2))
par(mar = c(1,4.5,1,0))
plot(nspec, r.mn.N, type="l", lwd = 2, ylim = c(0,1), xaxt="n", xlab = "",
     ylab = expression("r"["PLS"]), main = "A")
lines(nspec, r.mn.N - (1.96*r.sd.N), lty = 2)
lines(nspec, r.mn.N + (1.96*r.sd.N), lty = 2)

plot(nspec, cen.r.mn.N, type="l", lwd = 2, ylim = c(-0.2,0.2), xaxt="n", xlab = "",
     ylab = expression("Detrended r"["PLS"]), main = "C")
lines(nspec, cen.r.mn.N - (1.96*cen.r.sd.N), lty = 2)
lines(nspec, cen.r.mn.N + (1.96*cen.r.sd.N), lty = 2)

par(mar=c(4,4.5,1,0))
plot(nspec, z.mn.N, type="l", lwd = 2, ylim = c(-4,4), xlab = "Number of specimens", 
     ylab = "z", main = "E")
lines(nspec, z.mn.N - (1.96*z.sd.N), lty = 2)
lines(nspec, z.mn.N + (1.96*z.sd.N), lty = 2)

par(mar=c(1,2,1,2.5))
plot(p, r.mn.p, type = "l", lwd = 2, ylim = c(0,1), xaxt="n", xlab = "",
     ylab = "", yaxt = "n", main = "B")
lines(p, r.mn.p - (1.96*r.sd.p), lty = 2)
lines(p, r.mn.p + (1.96*r.sd.p), lty = 2)

plot(p, cen.r.mn.p, type = "l", lwd = 2, ylim = c(-0.2,0.2), xaxt = "n", xlab = "",
     ylab = "", yaxt = "n", main = "D")
lines(p, cen.r.mn.p - (1.96*cen.r.sd.p), lty = 2)
lines(p, cen.r.mn.p + (1.96*cen.r.sd.p), lty = 2)

par(mar=c(4,2,1,2.5))
plot(p, z.mn.p, type = "l", lwd = 2, ylim = c(-4,4), xlab = "Number of variables", 
     ylab = "", yaxt = "n", main = "F")
lines(p, z.mn.p - (1.96*z.sd.p), lty = 2)
lines(p, z.mn.p + (1.96*z.sd.p), lty = 2)


