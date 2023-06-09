#!/usr/bin/env Rscript
	
source("scripts/summarypop.R")

source("scripts/plotpredict.R")
source("scripts/Gmatrices.R")

col.sel <- "black"
col.control <- "darkgray"
pch.new <- 19
pch.old <- 1
lty.common <- 1
lty.old <- 3
lty.new <- 2

std.err <- TRUE
std.err.factor <- 2

axis.type="log"

tovar.mat <- get.matrices("Tovar")
tulum.mat <- get.matrices("Tulum")

tovar <- summarypop("data/TovarData.txt")
tulum <- summarypop("data/TulumData.txt")

tovar.old <- read.table("data/TovarSummary_old.txt", header=TRUE, na.string="-")
tulum.old <- read.table("data/TulumSummary_old.txt", header=TRUE, na.string="-")
# There are some differences in the encoding of both new and old data files. 
# The main one is that generations are labelled from 0 to 4 in the new data, and from 1 to 5 in the old data. 

drawline <- function(new.data, old.data, trait="x", rep) {
	my.mu <- paste0("mu.", trait)
	my.se <- paste0("se.", trait)
	
	col <- if (rep == "Control") col.control else col.sel
	gg <- if (rep == "Control") 1:3 else 1:4
	nd  <- new.data$Rep==rep & new.data$Gen %in% c(gg, 5)
	nd1 <- new.data$Rep==rep & new.data$Gen %in% gg
	nd2 <- new.data$Rep==rep & new.data$Gen %in% c(max(gg),5)
	od2 <- old.data$Rep==rep & old.data$Gen %in% c(max(gg)-1, 4)
	od3 <- old.data$Rep==rep & old.data$Gen == 4
	
	lines(gg, new.data[[my.mu]][nd1], type="b", lty=lty.common, col=col, pch=pch.new)
	lines(c(max(gg),5), new.data[[my.mu]][nd2], type="b", lty=lty.new, col=col, pch=c(NA, pch.new))
	lines(c(max(gg),4.9), old.data[[my.mu]][od2], type="b", lty=lty.old, col=col, pch=c(NA, pch.old))

	if (std.err) {
		arrows(x0=c(gg,5), y0=new.data[[my.mu]][nd]-std.err.factor*new.data[[my.se]][nd], y1=new.data[[my.mu]][nd]+std.err.factor*new.data[[my.se]][nd], col=col, lty=1, angle=90, length=0.02, code=3)
		arrows(x0=4.9, y0=old.data[[my.mu]][od3]-std.err.factor*old.data[[my.se]][od3], y1=old.data[[my.mu]][od3]+std.err.factor*old.data[[my.se]][od3], col=col, lty=1, angle=90, length=0.02, code=3)		
	}
}


pdf("figureS3.pdf", width=8, height=8, useDingbats = FALSE)
layout(cbind(1:2,3:4))
par(mar=c(5, 4.2, 3, 1))

G0 <- 0

ylim <- c(2.4,3.2)
plot(NULL, xlim=c(1,5), ylim=ylim, xlab="Generations", ylab=custom.label("GA", axis.type), xaxt="n", yaxt="n")
drawline(tovar, tovar.old, "x", "Control")
drawline(tovar, tovar.old, "x", "Up")
drawline(tovar, tovar.old, "x", "Down")
title("Tovar")
axis(1, at=order(unique(tovar$Gen)), labels=as.character(G0:(G0-1+max(tovar$Gen))))
custom.axis(2, ylim, axis.type)

ylim <- c(5.6, 6.3)
plot(NULL, xlim=c(1,5), ylim=ylim, xlab="Generations", ylab=custom.label("UBA", axis.type), xaxt="n", yaxt="n")
drawline(tovar, tovar.old, "y", "Control")
drawline(tovar, tovar.old, "y", "Up")
drawline(tovar, tovar.old, "y", "Down")
axis(1, at=order(unique(tovar$Gen)), labels=as.character(G0:(G0-1+max(tovar$Gen))))
custom.axis(2, ylim, axis.type)

ylim <- c(2.9,3.5)
plot(NULL, xlim=c(1,5), ylim=ylim, xlab="Generations", ylab=custom.label("GA", axis.type), xaxt="n", yaxt="n")
drawline(tulum, tulum.old, "x", "Control")
drawline(tulum, tulum.old, "x", "Up")
drawline(tulum, tulum.old, "x", "Down")
title("Tulum")
axis(1, at=order(unique(tovar$Gen)), labels=as.character(G0:(G0-1+max(tovar$Gen))))
custom.axis(2, ylim, axis.type)

ylim <- c(5.9,6.5)
plot(NULL, xlim=c(1,5), ylim=ylim, xlab="Generations", ylab=custom.label("UBA", axis.type), xaxt="n", yaxt="n")
drawline(tulum, tulum.old, "y", "Control")
drawline(tulum, tulum.old, "y", "Up")
drawline(tulum, tulum.old, "y", "Down")
axis(1, at=order(unique(tovar$Gen)), labels=as.character(G0:(G0-1+max(tovar$Gen))))
custom.axis(2, ylim, axis.type)

dev.off()
