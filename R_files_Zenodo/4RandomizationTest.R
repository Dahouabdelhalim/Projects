#Compare real and random

options(stringsAsFactors = FALSE)

modeldata<-read.csv("modeldata.csv")
modeldatasexed<-modeldata[-which(is.na(modeldata$sex)),]

proportion_infected_matrix<-read.csv("curveball_proportion_infected.csv", row.names = 1)

realnsim<-merge(modeldatasexed, proportion_infected_matrix, 
                by.x="dolphin_id", by.y=0,all.x=TRUE, all.y=FALSE)

library(logistf)

flog<-logistf(formula = tsd ~ proportion_infected + sex, data = modeldatasexed, pl = TRUE,
              alpha = 0.05, firth = TRUE)

actualp<-flog$prob[2]
actualB<-coef(flog)[2]
actualse<-sqrt(flog$var[2,2])

randp<-list()
randB<-list()
randse<-list()

for (i in 1:1000){
  
  flog<-logistf(formula = tsd ~ realnsim[,i+11] + sex, data = realnsim, pl = TRUE,
                alpha = 0.05, firth = TRUE)
  
  randp[[i]]<-flog$prob[2]
  randB[[i]]<-flog$coefficients[2]
  randse[[i]]<-sqrt(flog$var[2,2])
}


allrandB<-unlist(randB)
allrandse<-unlist(randse)
allrandz<-allrandB/allrandse
actualz<-actualB/actualse

Btile<-round(ecdf(allrandB)(actualB)*100,0)
Ztile<-round(ecdf(allrandz)(actualz)*100,0)



pdf(file="supplemental_histograms.pdf", width=9.5, height=5.9)
par(mfrow=c(1,2), xpd=NA, mar=c(5.1, 4.1, 2.1, 2.1))
hg<-hist(allrandB, probability = TRUE, col="grey85", border=NA,
         main=NA, 
         xlab="Coefficients generated from randomized data", 
         ylab=NA,
         axes=FALSE, 
         cex.lab=1.3)

axis(1, cex.axis=1.3)

segments(actualB,0,actualB, max(hg$density)*1.02, col="red", lty=3, lwd=2)
text(x=actualB, y=max(hg$density)*1.05,
     labels=expression(paste("Observed ", beta, " = 23rd %tile")))

y1<-density(allrandB, adjust=1.25, from=min(hg$breaks), to=max(hg$breaks))
lines(y1$x, y1$y, lwd=2, col="darkgrey", lty=2)

##plot 2
hg<-hist(allrandz, probability = TRUE, col="grey85", border=NA,
         main=NA, 
         xlab="Z-values generated from randomized data", 
         ylab=NA,
         xlim=c(min(allrandz),max(c(allrandz, actualz))),
         axes=FALSE, 
         cex.lab=1.3)

axis(1, cex.axis=1.3)

segments(actualz,0,actualz, max(hg$density)*1.02, col="red", lty=3, lwd=2)
text(x=actualz, y=max(hg$density)*1.05, labels="Observed z-value = 80th %tile")

y1<-density(allrandz, adjust=1.25, from=min(hg$breaks), to=max(hg$breaks))
lines(y1$x, y1$y, lwd=2, col="darkgrey", lty=2)

dev.off()






