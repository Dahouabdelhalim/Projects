#!/usr/bin/env Rscript

source("scripts/plotpredict.R")
source("scripts/summarypop.R")


source("scripts/Gmatrices.R")

axis.type="percent0"

tovar.mat <- get.matrices("Tovar")
tulum.mat <- get.matrices("Tulum")

tovar <- summarypop("data/TovarData.txt")
tulum <- summarypop("data/TulumData.txt")

CI.factor <- 2

pdf("figure4.pdf", width=8, height=8, useDingbats = FALSE)
layout(cbind(1:2,3:4))
par(mar=c(3, 4.2, 3, 1))
plot.data.recenter(
	recenter(tovar, G=tovar.mat$G, Gv=tovar.mat$Gv, P=tovar.mat$P, N=64, Np=12, target="mu.x", normalization="control"),
	main="Tovar", ylab="GA", xlab="", ylim=c(-0.4,0.4), CI.factor=CI.factor, axis.type=axis.type)
par(mar=c(5, 4.2, 1, 1))
plot.data.recenter(
	recenter(tovar, G=tovar.mat$G, Gv=tovar.mat$Gv, P=tovar.mat$P, N=64, Np=12, target="mu.y", normalization="control"),
	main="", ylab="UBA", xlab="Generation", ylim=c(-0.3,0.3), CI.factor=CI.factor, axis.type=axis.type)
par(mar=c(3, 4.2, 3, 1))
plot.data.recenter(
	recenter(tulum, G=tulum.mat$G, Gv=tulum.mat$Gv, P=tulum.mat$P, N=64, Np=12, target="mu.x", normalization="control", G0.boost=TRUE),
	main="Tulum", ylab="", xlab="", ylim=c(-0.4,0.4), CI.factor=CI.factor, axis.type=axis.type)
par(mar=c(5, 4.2, 1, 1))
plot.data.recenter(
	recenter(tulum, G=tulum.mat$G, Gv=tulum.mat$Gv, P=tulum.mat$P, N=64, Np=12, target="mu.y", normalization="control", G0.boost=TRUE),
	main="", ylab="", xlab="Generation", ylim=c(-0.3,0.3), CI.factor=CI.factor, axis.type=axis.type)
dev.off()
