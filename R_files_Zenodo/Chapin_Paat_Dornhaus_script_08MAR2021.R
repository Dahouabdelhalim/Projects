### R script for the publication: Brood as booty: the effect of colony size and resource value in social insect contests
### Authors: Kenneth James Chapin, Victor Alexander Paat, Anna Dornhaus
### Contact: chapinkj@gmail.com

##### Packages
require(effects); 
##### LOAD DATA
a <- read.csv('Chapin_Paat_Dornhaus_data.csv'); 
#obs <- read.csv('Chapin_Paat_Dornhaus_observationaldata.csv') 

colnames(a)
for(ii in 1:length(a$workersl)){
        if(a$workersl[ii] > a$workersr[ii]){ # convert "left" and "right" colonies to large (L) and small (S).
                a$workersLG[ii]  <- a$workersl[ii]
                a$workersSM[ii]  <- a$workersr[ii]
                a$broodLG[ii]    <- a$broodl[ii] # note that "larger" colonies may have fewer brood
                a$broodSM[ii]    <- a$broodr[ii]
                a$fightsSM[ii]   <- a$fightr[ii]
                a$fightsLG[ii]   <- a$fightl[ii]
                a$entranceSM[ii] <- a$entrr[ii]
                a$entranceLG[ii] <- a$entrl[ii]
        }else{
                a$workersSM[ii]  <- a$workersl[ii]
                a$workersLG[ii]  <- a$workersr[ii]
                a$broodLG[ii]    <- a$broodr[ii]
                a$broodSM[ii]    <- a$broodl[ii]
                a$fightsSM[ii]   <- a$fightl[ii]
                a$fightsLG[ii]   <- a$fightr[ii]
                a$entranceSM[ii] <- a$entrl[ii]
                a$entranceLG[ii] <- a$entrr[ii]
        }
        a$workerdif[ii] <- (a$workersL[ii] - a$workersS[ii])/a$workersL[ii] # make proportions and ratios
        a$ratiodif[ii] <- ((a$broodL[ii]/a$workersL[ii]) - (a$broodS[ii]/a$workersS[ii]))/(a$broodL[ii]/a$workersl[ii])
        a$brooddif[ii] <-  (a$broodL[ii] - a$broodS[ii])/a$broodL[ii]
}
max(a$brooddif)

fight <- glm(fightsSM ~  workerdif * brooddif * time, family="poisson", data=a, na.action=na.fail)
entr <- glm(entranceSM ~ workerdif * brooddif * time, family="poisson", data=a, na.action=na.fail)
summary(fight)
summary(entr)
ef <- allEffects(fight, xlevels=list(workerdif = seq(0, 0.5, 0.1), brooddif = seq(-1, 1, 0.5), time = seq(0, 180, 15)))
plot(ef, lines=list(multiline=TRUE, se=TRUE), rug=FALSE, colors=grey.colors(8), lwd=4)

#
fight <- glm(fightsSM ~  fightsLG * time, family="poisson", data=a, na.action=na.fail)
entr <- glm(entranceSM ~ entranceLG * time, family="poisson", data=a, na.action=na.fail)
summary(fight)
summary(entr)

ef <- effect("fightsLG:time", fight, xlevels=list(time = seq(0, 180, 15), fightsLG = seq(0, 10, 1)))
plot(ef, lines=list(multiline=TRUE, se=TRUE), rug=FALSE, colors=grey.colors(11), lwd=4)



# Supplemental
x <- unique(data.frame("value" = c(a$broodSM, a$broodLG, a$workersSM, a$workersLG),
                       "measure" = c(rep("brood", length(c(a$broodSM, a$broodLG))), 
                                     rep("workers", length(c(a$workersSM, a$workersLG)))),
                       "ids" = c(paste("SM", a$trial, sep=""), paste("LG", a$trial, sep=""))))

x$measure <- as.factor(x$measure)
x$ids <- as.factor(x$ids)
par(mfrow=c(1,1), mar=c(4,4,1,1))
barplot(x$value ~ x$measure + x$ids, beside=TRUE, las=1, ylim=c(0, 250),
        names.arg = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18),
        ylab="Total Number", xlab="Colonies")
legend("topleft", legend=c("Brood", "Workers"), fill=c("black", "grey"), bty="n");
box()


y <- c(a$broodl, a$broodr)
x <- c(a$workersl, a$workersr)
plot(y ~ x)
mod <- lm(y ~ x)
summary(mod)
abline(mod)


hist(a$smallFighters)


trials <- unique(a$trial)
        par(mfrow=c(1,2), mar=c(4, 4, 1, 1))
        clrs <- rainbow(n=length(trials), 0.6, 0.7)
        plot(0, 0, type="n", ylim=c(0, 15), xlim=c(0, 180), ylab="Individuals Fighting", xlab="Time (minutes)", xaxt="n")
        axis(1, at = seq(0, 180, 15))
        for(ii in 1:length(trials)){
                b <- subset(a, a$trial == trials[ii])
        lines(jitter(b$fightl, 1) ~ jitter(b$time, 3), lty=1, col = clrs[ii], lwd=3)
        lines(jitter(b$fightr, 1) ~ jitter(b$time, 3), lty=2, col = clrs[ii], lwd=3)
}

plot(0, 0, type="n", ylim=c(0, 5), xlim=c(0, 180), ylab="Individuals Fighting", xlab="Time (minutes)", xaxt="n")
axis(1, at = seq(0, 180, 15))
for(ii in 1:length(trials)){
        b <- subset(a, a$trial == trials[ii])
        lines(jitter(b$entrl, 1) ~ jitter(b$time, 3), lty=1, col = clrs[ii], lwd=3)
        lines(jitter(b$entrr, 1) ~ jitter(b$time, 3), lty=2, col = clrs[ii], lwd=3)
}
#

### continous study
obs <- read.csv("Chapin_Paat_Dornhaus_observationaldata.csv") 
obs$start <- as.numeric(obs$start) # data formatting
obs$end <- as.numeric(obs$end) 
obs$id <- as.factor(paste(substr(obs$ant1nest, 0, 1), sep="", obs$trial))  # make unique ids for colonies

cols <- c("grey", "grey","black", "black") # set colors
ltys <- c(1, 2, 1, 2) # set line types
behaviors <- unique(obs$behavior) 
ids <- unique(obs$id)
par(mfrow=c(4, 2), mar=c(0, 0, 3, 1), oma=c(4, 5, 0.5, 0), mgp=c(4, 0.3, 0), tcl=-0.2)
for(ii in 1:length(behaviors)){
        plot(0, 0, ylim=c(0, 0.12), xlim=c(0, 180), type="n", axes=FALSE)
        mtext(behaviors[ii], side = 3, line=0.5, cex=1.3)
        box(col="grey60")
        if(ii %in% c(1, 3, 5, 7)){
                axis(2, las=1, col="grey60")
                mtext("Density", side=2, line=3)
        }
        if(ii %in% c(7, 8)){
                axis(1, at=seq(0, 180, 60), col="grey60")
                mtext("Trial time (minutes)", side=1, line=2)
        }
        
        for(jj in 1:4){
                obse <- na.omit(subset(obs$start/60, obs$behavior==behaviors[ii] & obs$id==ids[jj]))
                if(length(obse) > 1){
                        obser <- density(obse, adjust=0.1)
                        lines(obser$y ~ obser$x, col=cols[jj], lty=ltys[jj], lwd=2)
                        rug(obse)
                }
        }
}


