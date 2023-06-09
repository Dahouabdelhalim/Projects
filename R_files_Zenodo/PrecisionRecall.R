#This source code reproduces the results 
#of Precision/Recall (PR) curves.

#author: Andres Aldana


#Upload the source list
data <- read.csv("~/RPS_Cont_2013-2020.csv", 
                 stringsAsFactors=FALSE)


#Define the risk factor (descriptor) and class (condition)
descriptor = "RAD"
condition = "PCS"
#The interval for descriptor
ini = 0
fin = 1
tick = 0.01
th = seq(ini,fin,tick)

#A data frame to save PR results for a given threshold
prTab = data.frame(Threshold = th, Precision = 0, Recall = 0)

#Compute PR
for(i in 1: nrow(prTab)){
  
  t = prTab[i,"Threshold"]
  print(t)
  pred = data[data[,descriptor] >= t , ]
  tp = pred[pred$Status == condition,]
  precision =  nrow(tp)/nrow(pred)
  positives = data[data[,"Status"] == condition ,]
  tp = positives[positives[,descriptor] >= t ,]  
  recall = nrow(tp)/nrow(positives)
  prTab[i,"Precision"] = precision
  prTab[i,"Recall"] = recall
  
}

#Plot PR curves
png("PR-PCS_RAD.png", width = 5, height = 5, units = "in", res = 300)
par(mar = c(5.1,4.1,4.1,2.1) + c(0,2,0,2))
plot(NULL, ylim=c(0,1),xlim=c(0,1),font.lab=1,font=1,las=1,main=condition,axes=F,xlab=descriptor,
     ylab="",mar=c(5.1,8.1,4.1,2.1))
ticks = seq(0,1,0.2)
axis(side = 1,at = ticks,labels = ticks,font=1)
inter = 2:100
prec = prTab[inter,"Precision"]
ran = nrow(data[data[,"Status"] == condition,])/nrow(data)
ma = max(c(na.omit(prec), ran))
mi = min(c(na.omit(prec), ran))
lines(prTab$Threshold[inter],(prec-mi)/(ma-mi),col="blue",lwd=2,lty=1)
ticks = seq(mi,ma,(ma-mi)/5)
axis(side=2,at = seq(0,1,0.2),labels = format(ticks,digits = 3),col="blue",col.axis="blue",font=1,las=1)
abline(h=(ran-mi)/(ma-mi), col="blue",lty=2,lwd=1)
rec = prTab[inter,"Recall"]
ma = max(na.omit(rec))
mi = min(c(na.omit(rec),0.45))
lines(prTab$Threshold[inter],(rec-mi)/(ma-mi),col="deeppink",lwd=2,lty=1)
ticks = seq(mi,ma,(ma-mi)/5)
axis(side=4,at = seq(0,1,0.2),labels = format(ticks,digits = 2),col="deeppink",col.axis="deeppink",font=1,las=1)
abline(h=(0.5-mi)/(ma-mi), col = "deeppink",lwd=1,lty=2)
#legend("topright", bty="n", legend=c("Precision","Recall"), lwd=2, col=c("blue","deeppink"))
dev.off()


#CPW - SPW
condition = "PCS"
ini = 1
fin = 10
tick = 1
th = seq(ini,fin,tick)

prTab = data.frame(Threshold = th, Precision = 0, Recall = 0)

for(i in 1: nrow(prTab)){
  
  t = prTab[i,"Threshold"]
  print(t)
  pred = data[data[,"CPW"] >= t & data[,"SPW"] >= 1000*t, ]
  tp = pred[pred$Status == condition,]
  precision =  nrow(tp)/nrow(pred)
  positives = data[data[,"Status"] == condition ,]
  tp = positives[ positives[,"CPW"] >= t & positives[,"SPW"] >= 2000*t,]  
  recall = nrow(tp)/nrow(positives)
  prTab[i,"Precision"] = precision
  prTab[i,"Recall"] = recall
  
}


png("PR-PCS_CSPW.png", width = 5, height = 5, units = "in", res = 300)
par(mar = c(5.1,4.1,4.1,2.1) + c(0,2,0,2))
plot(NULL, ylim=c(0,1),xlim=c(min(th),max(th)),font.lab=1,font=1,las=1,main="",axes=F,xlab="CPW",
     ylab="",mar=c(5.1,8.1,4.1,2.1))
axis(side = 1,at = pretty(th),labels = pretty(th),font=1)
axis(side = 3,at = pretty(th),labels = pretty(th)*2000,font=1,las=1)
mtext("SPW",side=3,outer=T,line=-2,adj=0.525,font=1)
inter = 1:length(th)
prec = prTab[inter,"Precision"]
ran = nrow(data[data[,"Status"] == condition,])/nrow(data)
ma = max(c(na.omit(prec), ran))
mi = min(c(na.omit(prec), ran))
lines(prTab$Threshold[inter],(prec-mi)/(ma-mi),col="blue",lwd=2,lty=1)
ticks = seq(mi,ma,(ma-mi)/5)
axis(side=2,at = seq(0,1,0.2),labels = format(ticks,digits = 3),col="blue",col.axis="blue",font=1,las=1)
abline(h=(ran-mi)/(ma-mi), col="blue",lty=2,lwd=1)
rec = prTab[inter,"Recall"]
ma = max(na.omit(rec),0.55)
mi = min(c(na.omit(rec),0.45))
lines(prTab$Threshold[inter],(rec-mi)/(ma-mi),col="deeppink",lwd=2,lty=1)
ticks = seq(mi,ma,(ma-mi)/5)
axis(side=4,at = seq(0,1,0.2),labels = format(ticks,digits = 2),col="deeppink",col.axis="deeppink",font=1,las=1)
abline(h=(0.5-mi)/(ma-mi), col = "deeppink",lwd=1,lty=2)
dev.off()
