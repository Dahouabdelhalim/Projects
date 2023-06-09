# ================================================================
# R CODE FOR CHARACTERIZING UNCERTAINTY AROUND WASH ESTIMATES
# ================================================================
# For the explanaition of the code, please refer to: 
# Ezbakhe, F. and PÃ©rez-Foguet, A. (2019) Estimating access to
# drinking water and sanitation: The need to account for
# uncertainty in trend analysis. Science of the Total Environment,
# doi: 10.1016/j.scitotenv.2019.133830

rm(list=ls(all=TRUE)); dev.off(); setwd("~/Desktop")

# ================================================================
# LIBRARIES AND FUNCTIONS FOR THE ANALYSIS
# ================================================================

# LIBRARIES
aux<-c("stringr","compositions","nleqslv","dplyr","readxl","PearsonDS","gam","Metrics","hydroGOF","ggplot2","gridExtra","RColorBrewer")
lapply(aux, require, character.only = TRUE)

# FUNCTIONS
#Pre-processing of JMP data as 4-parts compositions
calculate.compositions <- function(df){
  #In Water, y1=Piped, y2=Other improved, y3=Surface, y4=Other unimproved
  #In Sanitation, y1=Sewer, y2=Other improved, y3=Open defecation, y4=Other unimproved
  df<-na.omit(df)
  df<-df[df$x1+df$x3<=100 & df$x1>=df$x2,] 
  y1 <- df$x2/100
  y2 <- df$x1/100-df$x2/100
  y3 <- df$x3/100
  y4 <- round(1-df$x1/100-df$x3/100,3)
  return(data.frame("Source"=df$Source,"Year"=df$Year,"y1"=y1,"y2"=y2,"y3"=y3,"y4"=y4))
}
#Treatment of zeros
treat.zeros <- function(df,e){
  n<-nrow(df)
  y<-c("y1","y2","y3","y4")
  for (i in 1:n){
    aux<-df[i,c(y)]
    if (any(aux==0)){
      df[i,c(y)][df[i,c(y)]==0]<-e
      df[i,c(y)][df[i,c(y)]!=e]<-(1-sum(df[i,c(y)][df[i,c(y)]==e]))*df[i,c(y)][df[i,c(y)]!=e]
    }
  }
  df[,-c(1,2)]<-round(df[,-c(1,2)],4)
  return(df)
}
#Estimation of standard errors
estimate.errors <- function(df,e,emax1,emin1,emax2,emin2){
  n<-ncol(df)
  a1<-emin1^2; b1<-(e/(1-e))*(emax1^2-emin1^2)
  a2<-emin2^2; b2<-(e/(1-e))*(emax2^2-emin2^2)
  df$se1<-NA; df$se2<-NA; df$se3<-NA; df$se4<-NA
  for (i in 1:4){
    aux<-df[,n-4+i]
    df[,n+i]<-sqrt(a2+b2*(1-aux)/aux)*aux/100
    aux<-df[df$Source=="MICS" | df$Source=="DHS" | df$Source=="LSMS" | df$Source=="WHS",][,n-4+i]
    df[df$Source=="MICS" | df$Source=="DHS" | df$Source=="LSMS" | df$Source=="WHS",][,n+i]<-sqrt(a1+b1*(1-aux)/aux)*aux/100
    df[df$Source=="CEN",][,n+i]<-0*df[df$Source=="CEN",][,n-4+i]
  }
  df<-subset(df,select=-c(Source))
  return(df)
}
#Simulation of data
simulate.data <- function(df,e,n){
  t<-df$Year
  aux1 <- list()
  aux2 <- as.data.frame(matrix(NA,nrow=n,ncol=nrow(df))); colnames(aux2)<-t
  for (i in 1:4){
    for (j in 1:nrow(df)){
      mu<-eval(parse(text=paste("df$y",i,sep="")))[j]
      s<-eval(parse(text=paste("df$se",i,sep="")))[j]
      dbeta_param <- function(z,m,s,e) {
        x <- z^2
        y <- numeric(2)
        y[1] <- e + (1-2*e)*x[1]/(x[1]+x[2])-m
        y[2] <- (1-2*e)^2*x[1]*x[2]/((x[1]+x[2])^2*(x[1]+x[2]+1))-s^2
        y
      }
      aux <- nleqslv(c(0.1,0.2), dbeta_param, control=list(trace=0,btol=.0001,delta="newton"),m=mu,s=s,e=e)
      a <- aux$x[1]^2; b <- aux$x[2]^2
      aux2[,j] <- rpearsonI(n,a,b,e,1-2*e)
    }
    aux1[[i]]<-aux2
  }
  aux1<-list("y1"=aux1[[1]],"y2"=aux1[[2]],"y3"=aux1[[3]],"y4"=aux1[[4]])
  return(aux1)
}
#Regressions of data
regression.data <- function(df,dfsim,n,aa){
  tt <- seq(1990,2020,by=0.5)
  if (aa=="Standard"){
    aux1 <- as.data.frame(matrix(NA,nrow=length(tt),ncol=4)); colnames(aux1)<-c("y1","y2","y3","y4")
    aux2 <- as.data.frame(matrix(NA,nrow=n,ncol=length(tt))); colnames(aux2)<-tt
    aux <- list()
    for (i in 1:4){
      #Regression of JMP data
      t<-df$Year
      y<-df[,1+i]
      lmGAM <- gam(y ~ s(t,k), family=gaussian)
      t<-tt
      yGAM <- predict.Gam(lmGAM, data.frame(t))
      aux1[,i]<-yGAM
      #Regression of simulated data
      for (j in 1:n){
        t<-df$Year
        y <- t(dfsim[[i]][j,])
        lmGAM <- gam(y ~ s(t,k), family=gaussian)
        t<-tt
        yGAM<- predict.Gam(lmGAM,data.frame(t))
        aux2[j,]<-yGAM
      }
      aux[[i]]<-aux2
    }
    aux <- list("y1"=aux[[1]],"y2"=aux[[2]],"y3"=aux[[3]],"y4"=aux[[4]])
  }else if (aa=="Compositional"){
    signs <- rbind (c( 1, 1, -1, -1),c(1, -1, 0, 0),c(0, 0, 1, -1)); VV=gsi.buildilrBase(t(signs))
    aux1 <- as.data.frame(matrix(NA,nrow=length(tt),ncol=4)); colnames(aux1)<-c("y1","y2","y3","y4")
    aux2 <- array(numeric(), dim=c(n,length(tt),4)); dimnames(aux2) <- list(1:n,tt,c("y1","y2","y3","y4"))
    aux <- list()
    #Regression of JMP data
    t<-df$Year
    y<-df[,c("y1","y2","y3","y4")]
    cc <- acomp(y); cc <- ilr(cc,VV)
    ccGAM <- as.data.frame(matrix(NA, nrow=length(tt), ncol=3))
    for (i in 1:3){
      t<-df$Year
      lmGAM <- gam(cc[,i] ~ s(t,k), family=gaussian)
      t<-tt
      ccGAM[,i]<- predict.Gam(lmGAM,data.frame(t))
    }
    aux1[,] <- unclass(ilrInv(ccGAM,VV))
    #Regression of simulated data
    for (j in 1:n){
      y<-data.frame(t(dfsim[[1]][j,]),t(dfsim[[2]][j,]),t(dfsim[[3]][j,]),t(dfsim[[4]][j,]))
      cc <- acomp(y); cc <- ilr(cc,VV)
      ccGAM <- as.data.frame(matrix(NA, nrow=length(tt), ncol=3))
      for (i in 1:3){
        t<-df$Year
        lmGAM <- gam(cc[,i] ~ s(t,k), family=gaussian)
        t<-tt
        ccGAM[,i]<- predict.Gam(lmGAM,data.frame(t))
      }
      aux2[j,,] <- unclass(ilrInv(ccGAM,VV))
    }
    aux<-list("y1"=aux2[,,1],"y2"=aux2[,,2],"y3"=aux2[,,3],"y4"=aux2[,,4])
  }
  return(list("JMPdata"=data.frame("Year"=tt,aux1),"SIMdata"=aux))
}
#Plot regressions
plot.results <- function(part,ci,ym){
  col<-c(brewer.pal(5,"OrRd")[5],brewer.pal(5,"Blues")[5],brewer.pal(5,"Greens")[5]); siz<-c(0.009,0.4,0.5); lin<-c("solid","dotted","dotdash")
  xaxis <- "Year"; yaxis <- "Proportion";
  g <- ggplot() + xlab(xaxis) + ylab(yaxis) + ggtitle(paste("x",part,sep="")) + theme_bw()+ scale_y_continuous(limits = c(ym[1], ym[2]), breaks=seq(ym[1], ym[2],ym[3]))
  mm<-c("OLS","GAM")
  for (i in 1:length(mm)){
    #Data
    aux <- eval(parse(text=paste(mm[i],".JMPdata",sep="")))
    t <- aux[,1]; p <- aux[,1+part]
    g <- g + 
      geom_line(data = data.frame(t,p), aes(x=t, y=p),color=col[i],size=siz[2])
    #Simulation
    aux <- eval(parse(text=paste(mm[i],".SIMdata",sep="")))
    t <- as.numeric(colnames(aux[[part]]))
    for (j in 1:min(nrow(aux[[part]]),100)){
      p <- as.numeric(aux[[part]][j,])
      g <- g + 
        geom_line(data = data.frame(t,p), aes(x=t, y=p),color=col[i],size=siz[1],linetype=lin[1])
    }
    ml <- apply(aux[[part]], 2, function(x) quantile(x, (1-ci)/2))
    mu <- apply(aux[[part]], 2, function(x) quantile(x, (1+ci)/2))
    g <- g + 
      geom_line(data = data.frame(t,ml), aes(x=t, y=ml),color=col[i],size=siz[2],linetype=lin[3])+
      geom_line(data = data.frame(t,mu), aes(x=t, y=mu),color=col[i],size=siz[2],linetype=lin[3])
  }
  #Data points
  aux <- eval(parse(text=paste("JMPdata",sep="")))
  t <- aux$Year; p <-aux[,1+part]
  g <- g + geom_point(data = data.frame(t,p), aes(x=t, y=p),size=1.5,shape=21,fill="gray32")
  return(g)
}

# ================================================================
# CASE STUDY AND TYPE OF ANALYSIS
# ================================================================

#Country (Code+Full Name), Service (Water or Sanitation), Setting (Urban or Rural)
country<-c("MAR_Morocco")
service<-c("Sanitation")
setting<-c("Urban")

#Approach (Standard or Compositional)
#ni (number of simulations), re (precision of data), ci (confidence interval)
#rsemax1 and rsemin1 (max and min relative standards for household surveys from MICS, DHS, LSMS and WHS)
#rsemax2 and rsemin2 (max and min relative standards for other household surveys)
approach<-c("Compositional")
ni<-1000; re<-0.5*10^-3; ci<-0.95
rsemax1<-40; rsemin1<-4; rsemax2<-80; rsemin2<-8

# ================================================================
# APPLICATION
# ================================================================

#Load JMP data
JMPdata.allcountries<-read.csv("JMPdata.csv")

#Analysis
JMPdata.original<-subset(JMPdata.allcountries,JMPdata.allcountries$Country==country & JMPdata.allcountries$Service==service & JMPdata.allcountries$Setting==setting)
JMPdata<-calculate.compositions(JMPdata.original)
JMPdata<-treat.zeros(JMPdata,re)
JMPdata<-estimate.errors(JMPdata,re,rsemax1,rsemin1,rsemax2,rsemin2)
SIMdata<-simulate.data(JMPdata,re,ni)
k=1; aux<-regression.data(JMPdata,SIMdata,ni,approach); OLS.JMPdata<-aux[[1]]; OLS.SIMdata<-aux[[2]]
k=4; aux<-regression.data(JMPdata,SIMdata,ni,approach); GAM.JMPdata<-aux[[1]]; GAM.SIMdata<-aux[[2]]

# ================================================================
# SAVE RESULTS
# ================================================================

aux1<-list("OLS.JMPdata"=OLS.JMPdata,"OLS.SIMdata"=OLS.SIMdata)
aux2<-list("GAM.JMPdata"=GAM.JMPdata,"GAM.SIMdata"=GAM.SIMdata)
aux<-list("JMPdata"=JMPdata, "SIMdata"=SIMdata, "ResultsOLS"=aux1, "ResultsGAM"=aux2)
save(aux, file=paste(country,service,setting,approach,"Rdata",sep="."))

# ================================================================
# PLOT RESULTS
# ================================================================

yy<-c(-0.1,1.1,0.1) #Vertical axis from yy1 to yy2, by yy3
aux<-list(); for (i in 1:4){aux[[i]]<-plot.results(i,ci,yy)}
do.call(grid.arrange,c(aux, ncol=2))
