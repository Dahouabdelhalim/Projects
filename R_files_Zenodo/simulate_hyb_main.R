## Want to simulate hybridization through X generations, assuming random mating, no immigration, and no selection

source("simulate_func.R")

## Simulations for 1000 ind, 11 generations, starting proportions of 0.5, 0.75, 0.9

start0.5 <- simulate.hyb(1000,0.5,11)
start0.25 <- simulate.hyb(1000,0.25,11)
start0.05 <- simulate.hyb(1000,0.05,11)
start0.1 <- simulate.hyb(1000,0.1,12)


pdf("sim_5gen_10gen_1000ind.pdf", width=8,height=6)

par(mfrow=c(2,3))

plot(start0.5[[1]][6,], start0.5[[2]][6,], type="n", xlab="", ylab="", main="5 generations, Equal starting ratios", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.5[[1]][6,], start0.5[[2]][6,], col="gray45", cex=1.5)

plot(start0.25[[1]][6,], start0.25[[2]][6,], type="n", xlab="", ylab="", main="5 generations, 25% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.25[[1]][6,], start0.25[[2]][6,], col="gray45", cex=1.5)

plot(start0.1[[1]][6,], start0.1[[2]][6,], type="n", xlab="", ylab="", main="5 generations, 10% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.1[[1]][6,], start0.1[[2]][6,], col="gray45", cex=1.5)

plot(start0.5[[1]][11,], start0.5[[2]][11,], type="n", xlab="", ylab="", main="10 generations, Equal starting ratios", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.5[[1]][11,], start0.5[[2]][11,], col="gray45", cex=1.5)

plot(start0.25[[1]][11,], start0.25[[2]][11,], type="n", xlab="", ylab="", main="10 generations, 25% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.25[[1]][11,], start0.25[[2]][11,], col="gray45", cex=1.5)

plot(start0.1[[1]][11,], start0.1[[2]][11,], type="n", xlab="", ylab="", main="10 generations, 10% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.1[[1]][11,], start0.1[[2]][11,], col="gray45", cex=1.5)

mtext("Proportion of ancestry (q)", side=1, outer=T, line=-1)
mtext("Interspecific ancestry (Q)", side=2, outer=T, line=-1.5)

dev.off()


## Add some migrants

source("simulate_func.R")

## Simulations for 1000 ind, 11 generations, starting proportions of 0.5, 0.75, 0.9

start0.5 <- simulate.hyb(1000,0.5,12,imm.sp1=0.25)
start0.25 <- simulate.hyb(1000,0.25,12,imm.sp1=0.1,imm.sp2=0.1)
start0.05 <- simulate.hyb(1000,0.05,11)
start0.1 <- simulate.hyb(1000,0.1,12)


pdf("sim_5gen_10gen_1000ind.pdf", width=8,height=6)

par(mfrow=c(2,3))

plot(start0.5[[1]][6,], start0.5[[2]][6,], type="n", xlab="", ylab="", main="5 generations, Equal starting ratios", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.5[[1]][6,], start0.5[[2]][6,], col="gray45", cex=1.5)

plot(start0.25[[1]][6,], start0.25[[2]][6,], type="n", xlab="", ylab="", main="5 generations, 25% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.25[[1]][6,], start0.25[[2]][6,], col="gray45", cex=1.5)

plot(start0.1[[1]][6,], start0.1[[2]][6,], type="n", xlab="", ylab="", main="5 generations, 10% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.1[[1]][6,], start0.1[[2]][6,], col="gray45", cex=1.5)

plot(start0.5[[1]][11,], start0.5[[2]][11,], type="n", xlab="", ylab="", main="10 generations, Equal starting ratios", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.5[[1]][11,], start0.5[[2]][11,], col="gray45", cex=1.5)

plot(start0.25[[1]][11,], start0.25[[2]][11,], type="n", xlab="", ylab="", main="10 generations, 25% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.25[[1]][11,], start0.25[[2]][11,], col="gray45", cex=1.5)

plot(start0.1[[1]][11,], start0.1[[2]][11,], type="n", xlab="", ylab="", main="10 generations, 10% Species 1", axes=F, xlim=c(0,1), ylim=c(0,1))
axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
arrows(0,0,0.5,1, length=0, col="gray90")
arrows(0.5,1,1,0, length=0, col="gray90")
points(start0.1[[1]][11,], start0.1[[2]][11,], col="gray45", cex=1.5)

mtext("Proportion of ancestry (q)", side=1, outer=T, line=-1)
mtext("Interspecific ancestry (Q)", side=2, outer=T, line=-1.5)

dev.off()


