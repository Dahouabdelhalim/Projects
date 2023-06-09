library(ggplot2)
library(gridExtra)

#---------------------------------------------------------- fxns
panel.cor<-function(x,y,digits=2,prefix='',cex.cor,...){
  usr<-par("usr");on.exit(par(usr))  
  par(usr=c(0,1,0,1))
  r<-abs(cor(x,y,method='pearson',use="pairwise.complete.obs"))
  txt<-format(c(r,0.123456789),digits=digits)[1]
  test<-cor.test(x,y,method='pearson',use='pairwise.complete.obs')
  Signif<-ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
  text(0.5,0.75,cex=1.5,paste("r=",txt))
  text(.5,.25,cex=1.5,Signif)
}

panel.lm <-function(x,y,col=par('col'), bg = NA, pch=par('pch'), cex=1, col.smooth='red', ...){
  points(x,y,pch=pch,col=col,bg=bg,cex=cex)
  abline(stats::lm(y~x),col=col.smooth,...)
}

#----------------------------------------------------------------

#setwd('~/Desktop/nectar/amanda_project/')
avgs <- read.delim('F2.avgs.KA.txt')

str(avgs)

ka_sub <- data.frame(avgs$short.stamen, avgs$ln.area, avgs$ln.vol)

pairs(ka_sub,pch=19, cex = 1, upper.panel=panel.cor,lower.panel=panel.lm)

cor.test(avgs$short.stamen,avgs$ln.area) 
cor.test(avgs$ln.area,avgs$ln.vol) 
cor.test(avgs$short.stamen,avgs$ln.vol) 

