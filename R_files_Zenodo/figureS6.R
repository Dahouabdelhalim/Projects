#!/usr/bin/env Rscript

source("./scripts/plotpredict.R")
source("scripts/summarypop.R")

library(ellipse)

datadir <- "./data"
cols    <- c(Tovar="blue", Tulum="darkorange")
ltys    <- c(full=1, exself=2)
transp  <- 5

logGA.lim <- c(2.6,3.6)
logUBA.lim <- c(5.5,6.5)

axis.type <- "natural" # "log"

df <- list(
		Tovar = summarypop(paste(datadir, "TovarData.txt", sep="/")),
		Tulum = summarypop(paste(datadir, "TulumData.txt", sep="/"))
)



pdf("figureS6.pdf", width=5, height=5)

par(mar=c(4.5, 4.5 ,1, 1))

plot(NULL, xlim=logGA.lim, ylim=logUBA.lim, xlab=custom.label("GA", axis.type), ylab=custom.label("UBA", axis.type), asp=1, xaxt="n", yaxt="n")
custom.axis(1, logGA.lim,  axis.type)
custom.axis(2, logUBA.lim, axis.type)


for (pop in c("Tovar", "Tulum")) {
	load(paste0("./data/", pop, "_MCMC.RData"))
	full   <- mod$VCV
	load(paste0("./data/", pop, "_MCMC_excludingself.RData"))
	exself <- mod$VCV
	
	genet.VCV <- grep(colnames(full), pattern="animal")
	full   <- full[,genet.VCV]/10000
	exself <- exself[,genet.VCV]/10000
	
	centre <- c(df[[pop]][1, "mu.x"], df[[pop]][1, "mu.y"])
	
	for (i in 1:nrow(full))
		lines(ellipse(matrix(full[i,], ncol=2), centre=centre), col=makeTransparent(cols[pop], transp), lty=ltys["full"])
#~ 	for (i in 1:nrow(exself))
#~ 		lines(ellipse(matrix(exself[i,], ncol=2), centre=centre), col=makeTransparent(cols[pop], transp), lty=ltys["exself"])
	
	lines(ellipse(matrix((colMeans(full)), ncol=2), centre=centre), col=cols[pop], lty=ltys["full"], lwd=3)
	lines(ellipse(matrix((colMeans(exself)), ncol=2), centre=centre), col=cols[pop], lty=ltys["exself"], lwd=3)
	
	points(x=centre[1], y=centre[2], col=cols[pop], cex=2)
}

legend("topleft", col=cols, lty=1, lwd=3, legend=names(cols))
legend("bottomright", col="darkgray", lty=ltys, lwd=3, legend=c("Full data", "W/o self"))


dev.off()
