#####FOR THIS SECTION, MUST RUN ANALYSIS ON D. SIMULANS FIRST####

library(RColorBrewer)
colors <- brewer.pal(7,"Dark2")
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# par(family="Palatino")

myplot <- function(x, y, shade=0.1, axs=c(1,2), axlbs=c(1,2), boxlab=NA,boxlabsize=NA,boxcol=1,box=T, ...){
  plot(x,y,type="n",axes=F,...)
  if(is.element(1,axs)) abline(v=axis(1,mgp=c(3,0.5,0),tcl=0.3,lwd=0,lwd.ticks=1,labels=ifelse(is.element(1,axlbs),TRUE,FALSE)),col=add.alpha(1,shade))
  if(is.element(2,axs)) abline(h=axis(2,mgp=c(3,0.5,0),tcl=0.3,lwd=0,lwd.ticks=1,las=2,labels=ifelse(is.element(2,axlbs),TRUE,FALSE)),col=add.alpha(1,shade))
  if(is.element(3,axs)) axis(3,mgp=c(3,0.5,0),tcl=0.1,lwd=0,lwd.ticks=1,las=2,labels=NA)
  if(is.element(4,axs)) axis(4,mgp=c(3,0.5,0),tcl=0.1,lwd=0,lwd.ticks=1,las=2,labels=NA)
  if(box==T) box(col=add.alpha(1,0.1))
  if (!is.na(boxlab)){
  rect(par("usr")[1], par("usr")[4], par("usr")[2], 1.08*par("usr")[4], xpd=NA,col = add.alpha(boxcol,0.3),border=add.alpha(boxcol,0.3))
  text((par("usr")[1]+par("usr")[2])/2,1.04*par("usr")[4],boxlab,xpd=NA,font=3)
  }
}

mytimeplot <- function(x, y, shade=0.1, axs=c(1,2), axlbs=c(1,2), boxlab=NA, ...){
  plot(x,y,type="n",axes=F,...)
  if(is.element(1,axs)) abline(v=axis.POSIXct(1,x,mgp=c(3,0.5,0),tcl=0.3,lwd=0,lwd.ticks=1,labels=ifelse(is.element(1,axlbs),TRUE,FALSE)),col=add.alpha(1,shade))
  if(is.element(2,axs)) abline(h=axis(2,mgp=c(3,0.5,0),tcl=0.3,lwd=0,lwd.ticks=1,las=2,labels=ifelse(is.element(2,axlbs),TRUE,FALSE)),col=add.alpha(1,shade))
  if(is.element(3,axs)) axis(3,mgp=c(3,0.5,0),tcl=0.2,lwd=0,lwd.ticks=1,las=2,labels=NA)
  if(is.element(4,axs)) axis(4,mgp=c(3,0.5,0),tcl=0.2,lwd=0,lwd.ticks=1,las=2,labels=NA)
  box(col=add.alpha(1,0.1))
  if (!is.na(boxlab)){
    rect(par("usr")[1], par("usr")[4], par("usr")[2], 1.08*par("usr")[4], xpd=NA,col = add.alpha(1,0.1),border=add.alpha(1,0.1))
    text((par("usr")[1]+par("usr")[2])/2,1.04*par("usr")[4],boxlab,xpd=NA,font=3,cex=boxlabsize)
  }
}

newdat3 <- expand.grid(day=seq(0,15,length=1000),vial=vials$ID)
dens.pred <- predict(dens.m6,newdat3,type="response")
dens.pred.overall <-  predict(dens.m6,newdat3,type="response",re.form=NA)

#tiff("s1.tiff",height=7,width=7, units='in', res=600)
# par(mfrow=c(1,1),mar=c(5,5,3,3))
split.screen(rbind(c(0.05,1,0.05,.99),c(0,.05,0,1),c(0.05,1,0,0.05)))
split.screen(c(6,10))
r.squared <- c()
betas <- c()
for(i in 1:60) {
  screen(i+3)
  par(mar=c(0.1,.1,.4,.1),cex=0.6)
  treat <- vials$treatment[i]
  if(treat=='kin'){
	colour <- rgb(30/252, 30/252, 30/252)
	}else{
		colour <- rgb(175/252, 175/252, 175/252)
	}
  myplot(c(0,12),c(1,1),ylab="Cannibalised (cumulative)",xlab="Experimental Day",ylim=c(0,max(mouthparts.pervial)),axlbs=c(2),boxlab=vials$ID[i],boxlabsize=textsize,boxcol=colour,axs=c())
  points(1:12,mouthparts.pervial[,i],col=colour,pch=16,cex=0.8)
  model <- lm( long$mouthparts[which(long$vial==vials$ID[i])]~0+long$day[which(long$vial==vials$ID[i])]  ) 
  r.squared[i] <- summary(model)$r.squared
  betas[i] <- as.numeric(coef(model)[1])
  abline(model,col=colour)
  # log(seq(1,12,length=100))
  # lines(seq(0,12,length=100),predict.lines[which(newdat$vial==dimnames(mouthparts.pervial)$vial[i])],col=add.alpha(colors[ifelse(treat=="kin",1,2)],0.7))
  # log(seq(0,13,length=100))
  # lines(seq(0,13,length=100),log(seq(1,12,length=100))*exp(fixef(m1)['log(day)'] +  ranef(m1)$vial[i,1]),col=colors[2])
}
screen(3)
mtext("Time (1-12 days)",side=1,line=3.75, cex=1.4)
screen(2)
mtext("Mouthparts (cumulative)",side=2,line=2.75, cex=1.4)
close.screen(all=T)
#dev.off()


#tiff("s2.tiff",height=7,width=7, units='in', res=600)
# par(mfrow=c(1,1),mar=c(5,5,3,3))
split.screen(rbind(c(0.05,1,0.05,.99),c(0,.05,0,1),c(0.05,1,0,0.05)))
split.screen(c(6,10))
for(i in 1:60) {
  screen(i+3)
  par(mar=c(0.1,.1,.4,.1),cex=0.6)
  treat <- vials$treatment[i]
   if(treat=='kin'){
	colour <- rgb(30/252, 30/252, 30/252)
	}else{
		colour <- rgb(175/252, 175/252, 175/252)
	}
  myplot(c(0,12),c(1,1),ylab="Density",xlab="Experimental Day",ylim=c(0,max(density.pervial)),axlbs=c(2),boxlab=vials$ID[i],boxlabsize=textsize,boxcol=colour,axs=c())
  points(1:12,density.pervial[,i],col=colour,pch=16,cex=0.8)
  lines(newdat3$day[which(newdat3$vial==vials$ID[i])],dens.pred[which(newdat3$vial==vials$ID[i])],col=colour)
}
screen(3)
mtext("Time (1-12 days)",side=1,line=3.75, cex=1.4)
screen(2)
mtext("Larval density",side=2,line=2.75, cex=1.4)
close.screen(all=T)
#dev.off()