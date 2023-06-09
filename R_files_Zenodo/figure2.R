rm(list=ls())
library(lme4)
library(nlme)
library("emmeans")
library(multcomp)
library(ggplot2)
library(RLRsim)

setwd('/Users/cdesjonq/Documents/RÃ©dactionarticles/1_in-prep/2021_Desjonqu_socialhybrids/submission/211014_PNAS/data/')

xdata <- read.csv('signals-preferences.csv', dec=',')
#full model fitting ----
res <- lmer(peak_pref~sp*sex*treat+year+z.temp+(1|aggregation), data=xdata,REML=FALSE)

#Figure----

coefs <- fixef(res)
randeff <- ranef(res)
xdata$ranef <- NA
for (i in 1:length(xdata$aggregation)){
  xdata$ranef[i] <- randeff$aggregation[which(row.names(randeff$aggregation)==xdata$aggregation[i]), 1]
}

xdata$pred.data <- xdata$peak_pref-(xdata$z.temp*coefs[7]+coefs[5]*as.numeric(xdata$year=='2019')+coefs[6]*as.numeric(xdata$year=='2020')+xdata$ranef)

fems <- which(xdata$sex=='female')
means_f <- tapply(xdata$pred.data[fems], INDEX=paste(xdata$sp[fems],xdata$treat[fems], xdata$loc[fems]), FUN=mean)
sd_f <- tapply(xdata$pred.data[fems], INDEX=paste(xdata$sp[fems],xdata$treat[fems], xdata$loc[fems]), FUN=sd)
n_f <- tapply(xdata$pred.data[fems], INDEX=paste(xdata$sp[fems],xdata$treat[fems], xdata$loc[fems]), FUN=length)

mals <- which(xdata$sex=='male')
means_m <- tapply(xdata$pred.data[mals], INDEX=paste(xdata$sp[mals],xdata$treat[mals], xdata$loc[mals]), FUN=mean)
sd_m <- tapply(xdata$pred.data[mals], INDEX=paste(xdata$sp[mals],xdata$treat[mals], xdata$loc[mals]), FUN=sd)
n_m <- tapply(xdata$pred.data[mals], INDEX=paste(xdata$sp[mals],xdata$treat[mals], xdata$loc[mals]), FUN=length)
factors <- names(means_m)

se_f <- sd_f/sqrt(n_f)
se_m <- sd_m/sqrt(n_m)

ci_inf_f <- means_f-se_f
ci_inf_m <- means_m-se_m
ci_sup_f <- means_f+se_f
ci_sup_m <- means_m+se_m

means_table <- data.frame(means_males=means_m, se_males=se_m, means_females=means_f, se_females=se_f)

bg <- ifelse(substr(factors, start=4, stop=5)=='ho', yes='aquamarine2', no='orange')
bg[substr(factors, start=9, stop=9)!='F'&substr(factors, start=4, stop=5)=='ho'] <- 'blue'
bg[substr(factors, start=9, stop=9)!='F'&substr(factors, start=4, stop=5)=='he'] <- 'brown'

lims <- range(c(means_m, ci_inf_m, ci_sup_m, means_f, ci_sup_f, ci_inf_f), na.rm=TRUE)
bitmap("socialhybrids_zoomedin_controlled.jpg", width = 5000, height = 2000, units = 'px', res=900)
par(mar=c(4.1, 4.1, 1.1, 1.1))
layout(matrix(c(1,2,4,1,3,4), nr = 2, byrow=TRUE), widths=c(2,1,2))
plot(means_f,means_m, ylab='male signal frequency (Hz)', xlab='female peak preference (Hz)', xlim=lims, ylim=lims, col=bg, pch=16, cex=1.5, xaxt='n', yaxt='n',frame=FALSE, font=1.5)
mtext(text = c('(a)', expression(sp[low]), expression(sp[high])), side = 1, line = c(-29,-29,-29), at=c(155,180, 290))
arrows(x0 = ci_sup_f, y0 = means_m, x1 = ci_inf_f, y1 = means_m, length = 0, col=bg, angle=90, code=3)
arrows(x0 = means_f, y0 = ci_sup_m, x1 = means_f, y1 = ci_inf_m, length = 0.0, col=bg, angle=90, code=3)
legend('bottomright', legend = c('mixed allopatric', 'mixed sympatric', 'own allopatric', 'own sympatric'), col=c('brown', 'orange', 'blue', 'aquamarine2'), pch=16, bty='n')

het <- which(substr(factors,start= 4,stop=5)=='he')
hom <- which(substr(factors,start= 4,stop=5)=='ho')

axis(side = 2, at = seq(165, 305, by=20))
axis(side = 1, at = seq(165, 305, by=20))

plot(means_f,means_m, ylab='male signal frequency (Hz)', xlab='', xlim=c(270, 310), ylim=c(270, 310), col=bg, pch=16, cex=1.5, xaxt='n', yaxt='n',frame=FALSE)
mtext(text = c('(b)', expression(sp[high])), side = 1, line = c(-12,-12), at=c(272, 282))
arrows(x0 = ci_sup_f, y0 = means_m, x1 = ci_inf_f, y1 = means_m, length = 0, col=bg, angle=90, code=3)
arrows(x0 = means_f, y0 = ci_sup_m, x1 = means_f, y1 = ci_inf_m, length = 0.0, col=bg, angle=90, code=3)
axis(side = 2, at = seq(270, 310, by=10))
axis(side = 1, at = seq(270, 310, by=10))
arrows(x0 = means_f[c(1,2)], y0 = means_m[c(1,2)], x1 = means_f[c(3,4)], y1 = means_m[c(3,4)], length = 0, col='grey', lty=3, angle=90, code=3)


plot(means_f,means_m, ylab='male signal frequency (Hz)', xlab='female peak preference (Hz)', xlim=c(160, 200), ylim=c(160, 200), col=bg, pch=16, cex=1.5, xaxt='n', yaxt='n',frame=FALSE)
mtext(text = c('(c)', expression(sp[low])), side = 1, line = c(-12,-12), at=c(162, 172))
arrows(x0 = ci_sup_f, y0 = means_m, x1 = ci_inf_f, y1 = means_m, length = 0, col=bg, angle=90, code=3)
arrows(x0 = means_f, y0 = ci_sup_m, x1 = means_f, y1 = ci_inf_m, length = 0.0, col=bg, angle=90, code=3)
axis(side = 2, at = seq(160, 200, by=10))
axis(side = 1, at = seq(160, 200, by=10))
arrows(x0 = means_f[3:14], y0 = means_m[3:14], x1 = means_f[3:14], y1 = means_m[3:14], length = 0, col='grey', lty=3, angle=90, code=3)


females <- xdata[!is.na(xdata$strength),]

head(females)
str(females)

#full model fitting ----
females$log.strength <- log(females$strength)
res <- lmer(log.strength~sp*treat+year+z.temp+(1|aggregation), data=females)

coefs <- fixef(res)
randeff <- ranef(res)
females$ranef <- NA
for (i in 1:length(females$aggregation)){
  females$ranef[i] <- randeff$aggregation[which(row.names(randeff$aggregation)==females$aggregation[i]), 1]
}
females$pred.data <- females$log.strength-(females$z.temp*coefs[6]+coefs[4]*as.numeric(females$year=='2019')+coefs[5]*as.numeric(females$year=='2020')+females$ranef)


means_m <- tapply(females$pred.data, INDEX=paste(females$sp,females$treat, females$loc), FUN=mean)
sd_m <- tapply(females$pred.data, INDEX=paste(females$sp,females$treat, females$loc), FUN=sd)
n_m <- tapply(females$pred.data, INDEX=paste(females$sp,females$treat, females$loc), FUN=length)
factors <- names(means_m)

se_m <- sd_m/sqrt(n_m)

ci_inf_m <- means_m-se_m
ci_sup_m <- means_m+se_m

bg <- ifelse(substr(factors, start=4, stop=5)=='ho', yes='blue', no='orange')

col <- ifelse(females$treat=='homo', yes='aquamarine2', no='orange')
col[females$loc%in%c('BOG', 'OLT', 'PNV')&females$treat=='homo'] <- 'blue'
col[females$loc%in%c('BOG', 'OLT', 'PNV')&females$treat=='hete'] <- 'brown'

bg <- ifelse(substr(factors, start=4, stop=5)=='ho', yes='aquamarine2', no='orange')
bg[substr(factors, start=9, stop=9)!='F'&substr(factors, start=4, stop=5)=='ho'] <- 'blue'
bg[substr(factors, start=9, stop=9)!='F'&substr(factors, start=4, stop=5)=='he'] <- 'brown'

index <- ifelse(substr(factors, start=1, stop=2)=='LF', yes=1, no=3)+ifelse(substr(factors, start=4, stop=5)=='ho', yes=0, no=1)


females$index <- 0
females$index[females$treat=='homo'&females$sp=='LF'] <- 1
females$index[females$treat=='hete'&females$sp=='LF'] <- 2
females$index[females$treat=='homo'&females$sp=='HF'] <- 3
females$index[females$treat=='hete'&females$sp=='HF'] <- 4

#ylim <- range(c(females$pred.data,means_m, ci_inf_m, ci_sup_m))
ylim <- range(c(means_m, ci_inf_m, ci_sup_m))

plot(index,means_m, xaxt='n', ylab='preference strength (log)', xlab='', col=bg, pch=16, ylim=ylim, xlim=c(0.9, 4.3), frame=FALSE)
mtext(text = '(d)', side = 1, line = -29, at=0.3)
#points(jitter(females$index, 0.4), females$pred.data, col=alpha(col, 0.5), pch=16, cex=0.6)
arrows(x0 = index, y0 = ci_sup_m, x1 =index, y1 = ci_inf_m, length = 0, angle=90, code=3, col=bg)

het <- which(substr(factors,start= 4,stop=5)=='he')
hom <- which(substr(factors,start= 4,stop=5)=='ho')

arrows(x0=index[hom], x1=index[het],y0=means_m[hom], y1=means_m[het], length = 0, angle=90, code=3, col='grey', lty=3)
axis(side = 1, at = c(1, 1.5, 3.5, 4), labels = c('',expression(sp[low]), expression(sp[high]), ''))
mtext(side = 1, at = 2.5, text = c('species'), line = 2, font=1)
legend('bottomright', legend = c('mixed allopatric', 'mixed sympatric', 'own allopatric', 'own sympatric'), col=c('brown', 'orange', 'blue', 'aquamarine2'), pch=16, bty='n')

dev.off()

