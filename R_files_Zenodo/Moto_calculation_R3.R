### Calculate Moto from otolith d13C ######

All<-read.csv('Otolith_isotope_ratio',header = T,stringsAsFactors=T)

n=10000 #simulation numbers
set.seed(2020) # random choose fit 
#####################Moto estimation for JP bu Monte Carlo Simulations ##################################
otoCjp<-All[6][(All[1]=="JP"),]#otolith values setting for jp pop
dCjp<-runif(n, min = -22, max = -19)#dietary values setting for jp pop
wCjp<-runif(n, min = 0.53, max = 1.05)#water values setting for jp pop 
predMjp<-matrix(NA, ncol=n, nrow=length(otoCjp))
meanjp=NA
sdjp=NA

for(i in 1:length(otoCjp)){
  for (j in 1:n) {
    predMjp[i,j]<-(otoCjp[i]-wCjp[j])/(dCjp[j]-wCjp[j])
  }
  meanjp[i]<-mean(predMjp[i,], na.rm=TRUE)
  sdjp[i]<-sd(predMjp[i,], na.rm=TRUE)
}

#####################Moto estimation for CA by Monte Carlo Simulations ##################################
otoCca<-All[6][(All[1]=="CA"),]#otolith values setting for ca pop
dCca<-runif(n, min = -21.5, max = -18.5)#dietary values setting for ca pop
wCca<-runif(n, min = -0.31, max = 2.2)#water values setting for ca pop
predMca<-matrix(NA, ncol=n, nrow=length(otoCjp))
meanca=NA
sdca=NA

for(i in 1:length(otoCca)){
  for (j in 1:n) {
    predMca[i,j]<-(otoCca[i]-wCca[j])/(dCca[j]-wCca[j])
  }
  meanca[i]<-mean(predMca[i,], na.rm=TRUE)
  sdca[i]<-sd(predMca[i,], na.rm=TRUE)
}

#####################Combining Data for GAM##################################

All$Moto.mean<-c(meanjp,meanca)
All$Moto.sd<-c(sdjp,sdca)
write.csv(All,'Moto.csv')