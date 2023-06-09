#Aspect function, to plot Fourier spectrum of time series

#functions variables and parameters:
  #y = time series to plot
  #log=x : classic Fourier spectrum; log=xy : bilog Spectrum
  #dtrend: if 0, no detrending, if 1, linear detrending, if 2, loess detrending
  #AC : acf applied or not
  #filt : filter to be used : "blackman","hamming","hanning","bartlett","porte","triang","gausswin", or"flattopwin"
  #freqloc =T : to manually determine a frequency on the Fourier graph by clicking twice around the desired area (slope break)
  #recons = to reconstruct a frequency exctracted from the spectrum. Click twice to define bandwidth, then click right to make the reconstructed bandwidth appear
    #on a graph above the plotted time series. Then click 3 times to add parameters values on the graph
  #spect : spect="PSD" or "Amp" or "Norm" for power density spectrum, or power on y-axis, or normalised power.
  #scalcoef : to determine the slope of a part of the spectrum and plot it. Click twice to create to points on the graph, then click again on the two crosses that appear
     #then click again to add the slope coefficient directly on the graph

aspect <- function(y,tstep=1,detrend=1,AC=FALSE,sigtest=FALSE,clevel=0.95,AR1=NULL,filt="none",lfilt=length(y),wingauss=1,freqloc=FALSE,recons=FALSE,units=NULL,labX=NULL,labY=NULL,log="x",col="blue",spect="PSD",scalcoef=F,scalelim=F,nlim=1,wd=getwd())	{
  
  rc <- NA
  sig.in <- y
  y <- ts(y)
  absc <- time(y)
#   if (yr==TRUE) { convan <- 365.25 } else { convan <- 1 }
  setwd(wd)
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  
  if (detrend==0)  {
    y <- y
  }
  
  if (detrend==1)  {
    linfit <- lm(y~time(y))
    y <- linfit$residuals
  }
  
  if (detrend==2)  {
    yloessmod <- loess(y~absc,span=1)
    trend <- predict(yloessmod)
    y <- y-trend
  }
  acfy <- acf(y,lag.max=length(y)/3,plot=F,demean=T)
  
  if (AC)	{
    y <- acf(y,lag.max=length(y)/3,plot=F,demean=T)
    N <- length(y$acf)
    y <- acfy <- ts(y$acf)
    tps <- time(y)*tstep
    if (filt=="blackman") { require(signal); y <- y*blackman(lfilt)}
    if (filt=="hamming") { require(signal); y <- y*hamming(lfilt)}
    if (filt=="hanning") { require(signal); y <- y*hanning(lfilt)}
    if (filt=="bartlett") { require(signal); y <- y*bartlett(lfilt)}
    if (filt=="porte") { require(signal); y <- y*boxcar(lfilt)}
    if (filt=="triang") { require(signal); y <- y*triang(lfilt)}
    if (filt=="gausswin") { require(signal); y <- y*gausswin(length(y),wingauss)}
    if (filt=="flattopwin") { require(signal); y <- y*flattopwin(lfilt)}
    yfreq <- frequency(y)
    freq <- seq.int(from=yfreq/N,by=yfreq/N,length.out=floor(N/2))
    Per <- 1/freq
    ffty <- fft(y)
    if (spect=="PSD") { PSD <- Mod(1/N*ffty[2:floor((N/2)+1)]) ; spectunit <- "Power" }
    if (spect=="Amp") { PSD <- 2*Mod(ffty[2:floor((N/2)+1)])/length(y) ; spectunit <- "Amplitude" }
    if (spect=="Norm") { PSD <- N*(Mod(1/N*ffty[2:floor((N/2)+1)])^2)/var(y) ; spectunit <- "Normalized Power" }
    op <- par(mfrow=c(2,1),mar=c(4,4,1,1))
    plot(tps,acfy,type="l",lwd=2,col=col,xlab=units,ylab="Autocorrelation")
    plot(tim <- Per*tstep,PSD,type="l",lwd=2,col=col,log=log,axes=T,xlab="Period",xlim=c(max(tim),min(tim)),ylab=spectunit)
    if (sigtest==TRUE) {
      if (!is.null(AR1)) ar1 <- AR1
      else ar1 <- acfy[2]
      Per <- 1/freq
      p <- (1-ar1^2)/(1+ar1^2-2*ar1*cos(2*pi*(1/Per)))
      pCL=qchisq(clevel,df=2)*p/2
      lines(Per,p,lty="dashed",lwd=3,col="darkgrey")
      lines(Per,pCL,lty="dashed",lwd=2)
    }
    if (freqloc==TRUE) {
      par(new=T);plot(PSD,type="l",lwd=2,log=log,axes=F,xlab="",ylab="",col=col);points(PSD,type="h",col="darkgrey")
      borne <- locator(2)
      f1 <- borne$x[1]; if (f1<1) {f1 <- 1};legend(f1,max(PSD)/2,legend=round(1/freq[f1],2),box.lty=0);abline(v=f1)
      f2 <- borne$x[2]; if (f2<1) {f2 <- 1};legend(f2,max(PSD)/3,legend=round(1/freq[f2],2),box.lty=0);abline(v=f2)
    }
    if (scalcoef==TRUE) {
      plot(freq,PSD,type="l",lwd=0.5,log="xy",axes=T,xlab="Frequency",ylab=spectunit)
      pos <- identify(freq,PSD,tolerance=1e3,n=2,labels="+",col="grey")
      regfit <- lm(log(PSD[pos[1]:pos[2]])~log(freq[pos[1]:pos[2]]))
      locator(2,type="l",col="red");text(locator(1),coeff <- as.character(round(regfit$coefficients[2],2)))
    }
    par(op)
  }
  
  else	{
    N <- length(y)
    tps <- time(y)
    if (filt=="blackman") { require(signal); y <- y*blackman(lfilt)}
    if (filt=="hamming") { require(signal); y <- y*hamming(lfilt)}
    if (filt=="hanning") { require(signal); y <- y*hanning(lfilt)}
    if (filt=="bartlett") { require(signal); y <- y*bartlett(lfilt)}
    if (filt=="porte") { require(signal); y <- y*boxcar(lfilt)}
    if (filt=="triang") { require(signal); y <- y*triang(lfilt)}
    if (filt=="gausswin") { require(signal); y <- y*gausswin(length(y),wingauss)}
    if (filt=="flattopwin") { require(signal); y <- y*flattopwin(lfilt)}
    yfreq <- frequency(y)
    freq <- seq.int(from=yfreq/N,by=yfreq/N,length.out=floor(N/2))
    Per <- 1/freq
    ffty <- fft(y)
    if (spect=="PSD") { PSD <- Mod(1/N*ffty[2:floor((N/2)+1)])^2 ; spectunit <- "Power" }
    if (spect=="Amp") { PSD <- 2*Mod(ffty[2:floor((N/2)+1)])/length(y) ; spectunit <- "Amplitude" }
    if (spect=="Norm") { PSD <- (Mod(ffty[2:floor((N/2)+1)])^2)/N/var(y) ; spectunit <- "Normalized Power" }
#     if (spect=="Norm") { PSD <- N*(Mod(1/N*ffty[2:floor((N/2)+1)])^2)/(2*var(y)) ; spectunit <- "Normalized Power" }
    plot(tim <- Per*tstep,PSD,type="l",lwd=2,col=col,xlim=c(max(tim),min(tim)),log=log,axes=T,xlab="Period",ylab=spectunit)
    if (sigtest==TRUE) {
      if (!is.null(AR1)) ar1 <- AR1
      else ar1 <- acfy$acf[2]
      p <- (1-ar1^2)/(1+ar1^2-2*ar1*cos(2*pi*(1/Per)))
      #       p = (1 - ar1^2)/abs(1 - ar1 * exp(-2 * (0+1i) * pi * (1/Per)))^2
      pCL=qchisq(clevel,df=2)*p/2
      lines(tim,p,lty="dashed",lwd=3,col="darkgrey")
      lines(tim,pCL,lty="dashed",lwd=2)
    }
    if (freqloc==TRUE) {
      par(new=T);plot(PSD,type="l",lwd=2,log=log,axes=F,xlab="",ylab="",col=col);points(PSD,type="h",col="darkgrey")
      borne <- locator(2)
      f1 <- borne$x[1]; if (f1<1) {f1 <- 1};legend(f1,max(PSD)/2,legend=round(1/freq[f1],2),box.lty=0);abline(v=f1)
      f2 <- borne$x[2]; if (f2<1) {f2 <- 1};legend(f2,max(PSD)/3,legend=round(1/freq[f2],2),box.lty=0);abline(v=f2)
    }
    if (scalcoef==TRUE) {
      plot(freq,PSD,type="l",lwd=0.5,log="xy",axes=T,xlab="Frequency",ylab=spectunit)
      pos <- identify(freq,PSD,tolerance=1e5,n=2,labels="+",col="grey")
      regfit <- lm(log(PSD[pos[1]:pos[2]])~log(freq[pos[1]:pos[2]]))
      locator(2,type="l",col="red");text(locator(1),coeff <- as.character(round(regfit$coefficients[2],2)))
    }
    if (scalelim==TRUE) {
      sclim <- locator(nlim)
      for (i in 1:nlim) {
        abline(v=sclim$x[i],lty="dashed",lwd=1.5)
      }
    }	
  }
  
  
  if (recons)	{
    if (AC==F)   {   
      
      par(new=T);plot(PSD,type="l",lwd=2,log=log,axes=F,xlab="",ylab="",col=col);points(PSD,type="h",col="darkgrey")
      borne <- locator(2)
      f1 <- borne$x[1]; if (f1<1) {f1 <- 1};legend(f1,max(PSD)/2,legend=round(1/freq[f1],2),box.lty=0);abline(v=f1)
      f2 <- borne$x[2]; if (f2<1) {f2 <- 1};legend(f2,max(PSD)/3,legend=round(1/freq[f2],2),box.lty=0);abline(v=f2)
      if(f2<f1) { f1temp <- f1; f1 <- f2; f2 <- f1temp }
      if(f1<1) { f1 <- 1; if(f2<1) f2 <- 1; if(f1>floor(length(PSD))) f1 <- length(PSD); if(f2>floor(length(PSD))) f2 <- length(PSD) }
      ffty[1:f1] <- 0;ffty[(length(ffty)-f1+1):length(ffty)] <- 0
      ffty[f2:(length(ffty)-f2+1)] <- 0
#       var.contrib <- Mod(sum(ffty))^2/Mod(sum(fft(y)))^2
      rc <- Re(fft(ffty, inverse=T)/length(ffty))
      # X11()
      plot(y,type="l",col="darkgrey",lwd=2,ylim=c(min(min(rc),min(y)),max(y)),xlab=labX,ylab=labY) #enlève absc
      rc.contrib  <- 100*var(rc)/var(y)
      rc.sd <- var(rc)
      lines(rc+mean(y),col=col,lwd=2) #on enlève absc
      legend(locator(1),legend=c("Bandwidth:",c(round(1/freq[f1],2),round(1/freq[f2],2))),bty="o",bg="white",box.lty=0,ncol=1,xjust=0)
      legend(locator(1),legend=c("Reconstruction (Ã©cart-type):",round(rc.sd,1),"Ecart-type sÃ©rie initiale:",round(sd(y),1)),bg="white",box.lty=0)
      legend(locator(1),legend=c("Rapport Rc/Total (%):",round(rc.contrib,1)),bg="white",box.lty=0)
#       legend(locator(1),legend=c("Var (%):",round(var.contrib,1)),bg="white",box.lty=0)
    }
    
    if (AC==T)   {   
      if(detrend) {y <- sig.in-trend} else {y <- sig.in}
      N <- length(y)
      y <- ts(y)
      tps <- time(y)
      if (filt=="blackman") { require(signal); y <- y*blackman(lfilt)}
      if (filt=="hamming") { require(signal); y <- y*hamming(lfilt)}
      if (filt=="hanning") { require(signal); y <- y*hanning(lfilt)}
      if (filt=="bartlett") { require(signal); y <- y*bartlett(lfilt)}
      if (filt=="porte") { require(signal); y <- y*boxcar(lfilt)}
      if (filt=="triang") { require(signal); y <- y*triang(lfilt)}
      if (filt=="gausswin") { require(signal); y <- y*gausswin(length(y),wingauss)}
      if (filt=="flattopwin") { require(signal); y <- y*flattopwin(lfilt)}
      yfreq <- frequency(y)
      freq <- seq.int(from=yfreq/N,by=yfreq/N,length.out=floor(N/2))
      ffty <- fft(y)
      if (spect=="PSD") { PSD <- Mod(ffty[2:floor((N/2)+1)])^2 ; spectunit <- "Power" }
      if (spect=="Amp") { PSD <- 2*Mod(ffty[2:floor((N/2)+1)])/length(y) ; spectunit <- "Amplitude" }
      if (spect=="Norm") { PSD <- PSD <- length(y)*abs(ffty[2:floor((N/2)+1)])^2/(2*var(y)) ; spectunit <- "Normalized Power" }
      X11()
      plot(freq,PSD,type="h",col="grey",lwd=3,log=log,axes=T,xlab="Frequency",ylab=spectunit);par(new=T);plot(PSD,type="l",lwd=2,log=log,axes=F,xlab="",ylab="",col=col)
      borne <- locator(2)
      f1 <- borne$x[1]; if (f1<1) {f1 <- 1};legend(f1,max(PSD)/2,legend=round(1/freq[f1],2),box.lty=0);abline(v=f1)
      f2 <- borne$x[2]; if (f2<1) {f2 <- 1};legend(f2,max(PSD)/3,legend=round(1/freq[f2],2),box.lty=0);abline(v=f2)
      if(f2<f1) { f1temp <- f1; f1 <- f2; f2 <- f1temp }
      if(f1<1) { f1 <- 1; if(f2<1) f2 <- 1; if(f1>floor(length(PSD))) f1 <- length(PSD); if(f2>floor(length(PSD))) f2 <- length(PSD) }
      ffty[1:f1] <- 0;ffty[(length(ffty)-f1+1):length(ffty)] <- 0
      ffty[c((f2:(length(ffty)-f2+1)))] <- 0
      rc <- Re(fft(ffty, inverse=T)/length(ffty))
      # X11()
      plot(absc,y,type="l",col="darkgrey",lwd=2,ylim=c(min(min(rc),min(y)),max(y)),xlab=labX,ylab=labY)
      rc.contrib  <- 100*sd(rc)/sd(y)
      rc.sd <- sd(rc)
      lines(absc,rc+mean(y),col=col,lwd=2);legend(locator(1),legend=c("Bandwidth:",c(round(1/freq[f1],2),round(1/freq[f2],2))),bty="o",bg="white",box.lty=0,ncol=1,xjust=0);legend(locator(1),legend=c("Reconstruction (Ã©cart-type):",round(rc.sd,1),"Ecart-type sÃ©rie initiale:",round(sd(y),1)),bg="white",box.lty=0);legend(locator(1),legend=c("Rapport Rc/Total (%):",round(rc.contrib,1)),bg="white",box.lty=0)
    }
    #     if (save==T) { 
    #       if (is.null(name.spec)==TRUE) { name.spec="spectre" }
    #       if (is.null(name.recons)==TRUE) { name.recons="reconstruction" }
    #       write.table(cbind(freq,PSD),file=name.spec,append=FALSE,sep="\\t",dec=",");write.table(cbind(absc,rc),file=name.recons,append=FALSE,sep="\\t",dec=",")
    #     }
    recons <- invisible(rc)
  }
  result <- return(list(freq,PSD,rc))
}