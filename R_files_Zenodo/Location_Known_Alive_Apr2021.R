######################################################################################################################################
## Script by Arianne F. Messerman
## 13 April, 2021
## This script calculates the effect of relative apparent activity levels on juvenile spotted and marbled salamander (Ambystoma maculatum and A. opacum)
##   known survival at three time points that may be indicative of selective pressure. 
## Please see the associated manuscript for full study details: 
##   Messerman, A.F. and M. Leal. Submitted. The contributions of individual traits to survival among terrestrial juvenile 
##   pond-breeding salamanders. Ecology.
######################################################################################################################################

#setwd()

locs<-read.csv("Location_KnownAlive.csv", header=TRUE)
locs$alive<-as.factor(locs$alive)
locs$tag<-as.factor(locs$tag)
locs$pen<-as.factor(locs$pen)
locs$block<-as.factor(locs$block)
str(locs)

mod.loc<-glm(alive ~ period*prop.new.alive, locs, family = binomial(link = 'logit'))
summary(mod.loc)
library(car)
Anova(mod.loc, type="III")

library(emmeans)
emm.loc<-emmeans(mod.loc, ~period*prop.new.alive)
pairs(emm.loc, simple="period")

mod.loc1<-lm(prop.new.alive ~ period + alive, locs)
summary(mod.loc1)
Anova(mod.loc1)

emm.loc1<-emmeans(mod.loc1, ~period+alive)
pairs(emm.loc1, simple="period")
pairs(emm.loc1, simple="alive")

oct<-subset(locs, locs$period=="October")
nov<-subset(locs, locs$period=="November")
apr<-subset(locs, locs$period=="April")

oct.0<-subset(oct, oct$alive=="0")
oct.1<-subset(oct, oct$alive=="1")
nov.0<-subset(nov, nov$alive=="0")
nov.1<-subset(nov, nov$alive=="1")
apr.0<-subset(apr, apr$alive=="0")
apr.1<-subset(apr, apr$alive=="1")

prop.means<-c(mean(oct.1$prop.new.alive), mean(oct.0$prop.new.alive),
              mean(nov.1$prop.new.alive), mean(nov.0$prop.new.alive),
              mean(apr.1$prop.new.alive), mean(apr.0$prop.new.alive))

prop.sd<-c(sd(oct.1$prop.new.alive), sd(oct.0$prop.new.alive),
              sd(nov.1$prop.new.alive), sd(nov.0$prop.new.alive),
              sd(apr.1$prop.new.alive), sd(apr.0$prop.new.alive))

tiff("PropNewLoc-KnownAlive-300dpi.tiff", width = 10, height = 8, units = 'in', res=300, compression = 'none')
par(mfrow=c(1,1))
par(mar = c(5, 4, 4, 1) + 1.5)
barplot(prop.means, col=c("mediumorchid3","paleturquoise1","mediumorchid3","paleturquoise1","mediumorchid3","paleturquoise1"), 
        ylab="Proportion new locations/occasions alive", bty='l',
        xlab="Time period", ylim=c(0,2),cex.lab=2, cex.axis=2)
axis(1, c(1.3, 3.7, 6.1), cex.axis=2,
     labels=c("October", "November", "April"))
#mtext(side=3,at = c(.7, 1.9, 3.1, 4.3, 5.5, 6.7),c("A", "B", "A", "B", "C", "D"), cex=2)
segments( c(.7, 1.9, 3.1, 4.3, 5.5, 6.7), prop.means+prop.sd, c(.7, 1.9, 3.1, 4.3, 5.5, 6.7),prop.means-prop.sd, col=1, lwd=2)
abline(v=c(2.5, 4.9), lty=2, lwd=2, col="grey")
legend(5.2,1.9, cex=1.5, bty='n', legend=c("Known alive", "Not known alive"), pch=15, col=c("mediumorchid3","paleturquoise1"))
dev.off()