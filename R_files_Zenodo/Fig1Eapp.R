#Function to convert temperature in Celcius to Kelvin
T.K <- function(x, Tref = 15,kb=8.62E-5) 1/((Tref+273)*kb) - 1/((x + 273)*kb)

#Plot four typical linear growth~temperature relationships and the fitness landscape under a given environmental temperature
#different xbar & Topt, 
dx   <- 7.5
xbar <- seq(dx, 30-dx, length.out = 4)
theta<- xbar + dx/2
xmin <- xbar - dx/2
xbar <- T.K(xbar)
theta<- T.K(theta)
xmin <- T.K(xmin)
Ea   <- 0.65
b0   <- log(2)
b    <- b0 - Ea * theta
ym   <- theta*Ea + b 

pdf('Fig1EappEinter.pdf', height=4, width=8, useDingbats=FALSE)
op <- par(font.lab = 1,
           family  = "serif",
           cex.lab = 1.9,
           oma     = c(2,2,0.1,0.1),
           mgp     = c(1.5,.1,0),
           mar     = c(3,3,2,0.1),
           mfrow   = c(1,3))

temp    <- seq(-2, 35, 0.1)
temp    <- T.K(temp)
yrange  <- c(0.5, 3)

plot(temp, temp**2, 
     xlim = c(T.K(0), T.K(30)), 
     ylim = log(yrange), 
     xlab = expression(paste( italic(x)*' ('*eV^-1*')')),
     ylab = expression(paste( italic(y) ) ),
     xaxs = "i",
     yaxs = "i",
     xaxt = 'n',
     yaxt = 'n',
     type = 'n')

pool <- data.frame(X = c(), Y = c())

#Plot each TPC
for (i in 1:length(xbar)){
  curve(Ea*x + b[i], from = xmin[i], to = theta[i], add = T)
  #add data points
  X <- seq(xmin[i], theta[i], length.out = 4)
  Y <- Ea*X + b[i]
  points(X,Y)
  pool <- rbind(pool, data.frame(X,Y))
}

LMp <- lm(Y ~ X, pool)
abline(LMp, col=2, lwd=1.5)

#Add vertical line of 0
abline(v=0, lty = 3, col = 2)

#add xbar [j]
segments(x0 = xbar[1], y0 = log(0.5), 
         x1 = xbar[1], y1 = Ea*xbar[1] + b[1],
         lty = 2, col = 'blue')

#add b[j]
segments(x0 = theta[1], y0 = Ea*theta[1] + b[1], 
         x1 = 0,   y1 = b[1],
         lty = 2, col = 2)

segments(x0 = 0,      y0 = b[1], 
         x1 = T.K(0), y1 = b[1],
         lty = 2, col = 2)

#add xm[j] and ym[j]
segments(x0 = theta[1], y0 = log(yrange[1]), 
         x1 = theta[1], y1 = Ea*theta[1] + b[1],
         lty = 2, col = 3)
segments(x0 = theta[1], y0 = Ea*theta[1] + b[1], 
         x1 = T.K(0),   y1 = Ea*theta[1] + b[1],
         lty = 2, col = 3)

mtext(expression(paste('A) '*italic(E[app])*' & '*italic(E[intra]))), adj=0, cex=1.4)
mtext(expression(paste(bar(italic( x[j] )))), side = 1, adj = .25, col = 'blue', line = .5)
mtext(expression(paste(italic(x[mj]))),  side = 1, adj = .4,  col = 3, line = .5)
mtext('0', side = 1, adj = .52, col = 2, line=0.1)
mtext(expression(paste(italic( b[j])) ),      side = 2, adj = .97, col = 2)
mtext(expression(paste(italic(y[mj]))), side = 2, adj = .77,  col = 3)

#b) plot Ee (b ~ xbar)
plot(xbar, b, 
     xlim = c(T.K(0), T.K(30)), 
     ylim = log(yrange), 
     xlab = expression(paste(bar(italic(x)) *' ('*eV^-1*')')),
     ylab = expression(paste(italic(b))),
     xaxs = "i",
     yaxs = "i",
     xaxt = 'n',
     yaxt = 'n')
abline(lm(b ~ xbar), lwd = 1.5, col = 2 )

#add xbar [j]
segments(x0 = xbar[1], y0 = log(0.5), 
         x1 = xbar[1], y1 = b[1],
         lty = 2, col = 'blue')

#add b[j]
segments(x0 = T.K(0),  y0 = b[1], 
         x1 = xbar[1], y1 = b[1],
         lty = 2, col = 2)

#Add vertical line of 0
abline(v=0, lty = 3, col = 2)
mtext('0', side = 1, adj = .52, col = 2)
mtext(expression(paste(bar(italic( x[j] )))), side = 1, adj = .25, col = 'blue', line = .5)

mtext(expression(paste(italic( b[j])) ),      side = 2, adj = .97, col = 2)

mtext(expression(paste('B) '* italic(E[inter]))), adj=0, cex=1.4)

#c) plot ym against xm
plot(theta, ym, 
     xlim = c(T.K(0), T.K(30)), 
     ylim = log(yrange), 
     xlab = expression(paste(italic(x[m]) *' ('*eV^-1*')')),
     ylab = expression(paste(italic(y[m]))),
     xaxs = "i",
     yaxs = "i",
     xaxt = 'n',
     yaxt = 'n')
abline(lm(ym ~ theta), lwd = 1.5, col = 3 )

segments(x0 = theta[1], y0 = log(yrange[1]), 
         x1 = theta[1], y1 = Ea*theta[1] + b[1],
        lty = 2, col = 3)        
abline(v=0, lty = 3, col = 2)
mtext(expression(paste(italic(x[mj]))), side = 1, adj = .4,  col = 3, line = .2)
mtext(expression(paste(italic(y[mj]))), side = 2, adj = .77,  col = 3)
mtext('0', side = 1, adj = .52, col = 2, line=0.1)
mtext(expression(paste('C) '* italic(E[L]))), adj=0, cex=1.4)
dev.off()
