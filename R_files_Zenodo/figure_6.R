library(e1071)
library(ineq)


rm(list = ls())

setwd(".../Data_and_Code")

load("plant_heights.Rdata")
ls()

print(plant_heights)


rhodo_atl <- plant_heights$rhodo_heights_aboveTL
rhodo_btl <- plant_heights$rhodo_heights_belowTL
abies_atl <- plant_heights$abies_heights_aboveTL
abies_btl <- plant_heights$abies_heights_belowTL



#           Rhododendron above treeline 
# ==============================================

hist(rhodo_atl, breaks=20)
skewness(rhodo_atl)             #g1 with package e1071 on original  data
ineq(rhodo_atl, type="Gini")    # Gini coefficient
Lasym(rhodo_atl)                # Lorenz asymmetry coefficient. it worked in regular R editor

# exponential regression
rcampATL <- hist(rhodo_atl, breaks=20)
rcampATL 
rm(nrow, seq)
nrow <- nrow(as.data.frame(rcampATL$counts));	nrow
seq <- seq(nrow)*10 - 5; seq
rcampATLdf <- data.frame(size = seq, freq = rcampATL$counts);	rcampATLdf 
lines(rcampATLdf$size, rcampATLdf$freq)

# fit non-linear model
mod_rcampATL <- nls(freq ~ exp(a + b * size), data = rcampATLdf, start = list(a = 0, b = 0))
summary(mod_rcampATL)
lines(rcampATLdf$size, predict(mod_rcampATL, list(x = rcampATLdf$size)), col='blue', lwd=2)



#           Rhododendron below treeline 
# ==============================================

hist(rhodo_btl, breaks=20)
skewness(rhodo_btl)
ineq(rhodo_btl, type="Gini")
Lasym(rhodo_btl) 

# exponential regression
rcampBTL <- hist(rhodo_btl, breaks=20)
rcampBTL 
rm(nrow, seq)
nrow <- nrow(as.data.frame(rcampBTL$counts));	nrow
seq <- seq(nrow)*10 - 5; seq
rcampBTLdf <- data.frame(size = seq, freq = rcampBTL$counts);	rcampBTLdf 
lines(rcampBTLdf$size, rcampBTLdf$freq)

# fit non-linear model
mod_rcampBTL <- nls(freq ~ exp(a + b * size), data = rcampBTLdf, start = list(a = 0, b = 0))
summary(mod_rcampBTL)
lines(rcampBTLdf$size, predict(mod_rcampBTL, list(x = rcampBTLdf$size)), col='blue', lwd=2)



#              Abies above treeline 
# ==============================================

hist(abies_atl, breaks=20)
skewness(abies_atl)
ineq(abies_atl, type="Gini")
Lasym(abies_atl) 

# exponential regression
abiesATL <- hist(abies_atl, breaks=20)
abiesATL 
rm(nrow, seq)
nrow <- nrow(as.data.frame(abiesATL$counts));	nrow
seq <- seq(nrow)*10 - 5; seq
abiesATLdf <- data.frame(size = seq, freq = abiesATL$counts);	abiesATLdf 
lines(abiesATLdf$size, abiesATLdf$freq)

# fit non-linear model
mod_abiesATL <- nls(freq ~ exp(a + b * size), data = abiesATLdf, start = list(a = 0, b = 0))
summary(mod_abiesATL)
lines(abiesATLdf$size, predict(mod_abiesATL, list(x = abiesATLdf$size)), col='blue', lwd=2)



#              Abies below treeline 
# ==============================================
hist(abies_btl, breaks=20)
skewness(abies_btl)
ineq(abies_btl, type="Gini")
Lasym(abies_btl) 

# exponential regression
abiesBTL <- hist(abies_btl, breaks=20)
abiesBTL 
rm(nrow, seq)
nrow <- nrow(as.data.frame(abiesBTL$counts));	nrow
seq <- seq(nrow)*10 - 5; seq
abiesBTLdf <- data.frame(size = seq, freq = abiesBTL$counts);	abiesBTLdf 
lines(abiesBTLdf$size, abiesBTLdf$freq)

# fit non-linear model
mod_abiesBTL <- nls(freq ~ exp(a + b * size), data = abiesBTLdf, start = list(a = 0, b = 0))
summary(mod_abiesBTL)
lines(abiesBTLdf$size, predict(mod_abiesBTL, list(x = abiesBTLdf$size)), col='blue', lwd=2)





# MAKE NICE COMPOSITE HISTOGRAM PLOT of , ONE OF UNTRANSFORMED DATA, THE OTHER OF LOG TRANSFORMED DATA
# ------------------------------------------------------------------------------------------------


pdf("Figure_6_reverseJ_10cm_by_histogram.pdf", width=6, height=5)

# dev.new(width=6, height=6)
par(mfrow = c(2, 2))
par(oma = c(2, 2, 1, 0)) # make room  for the overall x and y axis titles
par(mar = c(2,2,1, 1)) # make the plots be closer together
par(ps = 8.5, cex = 1, cex.main = 1)

    hist(rhodo_atl, breaks=20, col="gray50", border="white", xlab="", ylab="", main="")
    text(x = 100, y = 150, cex=0.95, labels = expression(paste(italic("R. campanulatum "), "above treeline")))
    text(x = 100, y = 150*0.875, cex=0.95, labels = expression(paste(italic("G "), "= 0.40,  ", italic("S "), "= 0.74,  ", italic("g"[1]), " = 0.26")))
    text(x = 190, y = 150*1.15, labels = "(a)", font=2)
    lines(rcampATLdf$size, predict(mod_rcampATL, list(x = rcampATLdf$size)), lwd=1)
    
    hist(abies_atl, breaks=20, col="gray50", border="white", xlab="", ylab="", xlim=c(0,200), main="")
    text(x = 100, y = 50, cex=0.95, labels = expression(paste(italic("A. spectabilis "), "above treeline")))
    text(x = 100, y = 50*0.875, cex=0.95, labels = expression(paste(italic("G "), "= 0.48,  ", italic("S "), "= 0.97,  ", italic("g"[1]), " = 2.04")))
    text(x = 190, y = 50*1.15, labels = "(b)", font=2)
    lines(abiesATLdf$size, predict(mod_abiesATL, list(x = abiesATLdf$size)), lwd=1)
    
    hist(rhodo_btl, breaks=20, col="gray50", border="white", xlab="", ylab="", main="")
    text(x = 100, y = 200, cex=0.95, labels = expression(paste(italic("R. campanulatum "), "below treeline")))
    text(x = 100, y = 200*0.875, cex=0.95, labels = expression(paste(italic("G "), "= 0.53,  ", italic("S "), "= 0.67,  ", italic("g"[1]), " = 0.48")))
    text(x = 190, y = 200*1.15, labels = "(c)", font=2)
    lines(rcampBTLdf$size, predict(mod_rcampBTL, list(x = rcampBTLdf$size)), lwd=1)
    
    hist(abies_btl, breaks=20, col="gray50", border="white", xlab="", ylab="", main="")
    text(x = 100, y = 210, cex=0.95, labels = expression(paste(italic("A. spectabilis "), "below treeline")))
    text(x = 100, y = 210*0.875, cex=0.95, labels = expression(paste(italic("G "), "= 0.59,  ", italic("S "), "= 1.08,  ", italic("g"[1]), " = 2.89")))
    text(x = 190, y = 210*1.15, labels = "(d)", font=2)
    lines(abiesBTLdf$size, predict(mod_abiesBTL, list(x = abiesBTLdf$size)), lwd=1)
    
    mtext(side=2, line=0.35, cex=1.15, "Number of plant individuals", outer=TRUE) 
    mtext(side=1, line=0.35, cex=1.15,  "Size class (plant height) at 10 cm increment", outer=TRUE) 

dev.off()
