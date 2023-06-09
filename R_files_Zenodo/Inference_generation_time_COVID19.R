#####################################
# This code uses a Maximum Composite Likelihood approach
# for the inference of COVID-19 generation time
# and a Bayesian approach
# for the probability of pre-symptomatic transmission
# (Ferretti, Wymant et al, Science 2020)
#####################################


# maximum duration of infection in days
M<-30
# growth rate of the epidemic
r<-log(2)/5

# incubation time
# Li et al: mean 5.2, 95% quantile 12.5
log_incubation_sd1<-qnorm(0.95)-sqrt(qnorm(0.95)^2-2*(log(12.5)-log(5.2)))
log_incubation_median1<-log(12.5)-log_incubation_sd1*qnorm(0.95)
# Lauer et al: median 5.2, 2.5-97.5% quantile 2.5-10.5, mean 5.5
log_incubation_sd2<-qnorm(0.975)-sqrt(qnorm(0.975)^2-2*(log(10.5)-log(5.5)))
log_incubation_median2<-log(10.5)-log_incubation_sd2*qnorm(0.975)
# Wallinga et al: Weibull mean 6.4, 2.5-97.5% quantile 2.1-11.1
log_incubation_median3<-NA
log_incubation_sd3<-NA
# Yang et al: median 4.8, 25-75% quantile 3.0-7.2
log_incubation_median4<-log(4.8)
log_incubation_sd4<-(log(7.2)-log_incubation_median4)/qnorm(0.75)

# Use Lauer et al:
log_incubation_median<-log_incubation_median2
log_incubation_sd<-log_incubation_sd2

# incubation time
inc<-function(x){
  dlnorm(x,meanlog = log_incubation_median, sdlog = log_incubation_sd)
}

Inc<-sapply(0:(2*M),function(j){integrate(inc,lower = max(j-0.5,0),upper = j+0.5)$value})

###########
# Composite Likelihood and auxiliary functions
############

getlimits<-function(info){
  v<-list()
  v$T1L<-max(info$T1L,info$s1-2*M,na.rm=T)
  v$T1R<-min(info$T1R,info$s1,info$s2,info$T2R,na.rm=T)
  v$T2L<-max(info$T2L,info$s2-2*M,v$T1L,na.rm=T)
  v$T2R<-min(info$T2R,info$s2,na.rm=T)
  return(v)
}

LtermT<-function(info,W,limits){
  #print(info)
  #print(limits$T1L:limits$T1R)
  #print(limits$T2L:limits$T2R)
  sum(sapply(limits$T1L:limits$T1R,function(t1){
    sum(sapply(max(t1,limits$T2L):limits$T2R,function(t2){
      exp(r*info[["Exponential.growth"]]*t1)*Inc[info$s1-t1+1]*W[t2-t1+1]*Inc[info$s2-t2+1]
    }))
  }))
}

LtermS<-function(info,W,limits){
  sum(sapply(limits$T1L:limits$T1R,function(t1){
    (info$s1<=limits$T2R)*
    sum(sapply(max(t1,limits$T2L,info$s1):limits$T2R,function(t2){
      exp(r*info[["Exponential.growth"]]*t1)*Inc[info$s1-t1+1]*W[t2-t1+1]*Inc[info$s2-t2+1]
    }))
  }))
}

LtermN<-function(info,W,limits){
  sum(sapply(limits$T1L:limits$T1R,function(t1){
    sum(sapply(max(t1,limits$T2L):limits$T2R,function(t2){
      exp(r*info[["Exponential.growth"]]*t1)*Inc[info$s1-t1+1]*W[t2-t1+1]*sum(Inc[c(t2:min(info$Tr,t2+M))-t2+1])
    }))
  }))
}

logL<-function(info,W){
  limits<-getlimits(info)
  #log(LtermT(info,W,limits))-log(LtermN(info,W,limits))
  log(LtermT(info,W,limits))
}

fracSymptomatic<-function(info,W){
  limits<-getlimits(info)
  return(LtermS(info,W,limits)/LtermT(info,W,limits))
}

###########
# load data
###########

info_data<-read.csv("SupplementaryTableTransmissionPairs.csv",sep="\\t",stringsAsFactors = F)
location<-info_data$DiagnosisCountry
info_data[,"Exponential.growth"]<-0+(info_data[,"Exponential.growth"]=="YES")
info_data<-apply(info_data[,c("T1L","T1R","s1","T2L","T2R","s2","Tr","Exponential.growth")],2,function(x){as.numeric(x)})

location[location=="China"]<-"China"
location[location=="Hong Kong"]<-"China (HK)"

barplot(table(location),las=2,main="",cex.names=0.8)
##dev.copy2pdf(file="geographic_distribution_our_cases.pdf")

hist(info_data[,"s2"]-info_data[,"s1"],breaks=c(0:10)+0.5,main="empirical distribution of serial intervals",xlab="serial interval",ylab="",col="darkgrey")
##dev.copy2pdf(file="serial_interval_distribution.pdf")

###########
# optimisation
###########

neglogL2optim<-function(pars){
  #print(pars)
  W<-sapply(0:(4*M),function(j){integrate(function(x){w(x,pars)},lower = max(j-0.5,0),upper = j+0.5)$value})
  return(-sum(apply(info_data,1,function(x){logL(as.list(x),W)})))
}

# shape of generation time distribution
w<-function(x,pars){
  dlnorm(x, meanlog = pars[1], sdlog = exp(pars[2]))
  #dgamma(x, shape = exp(pars[1]), rate = exp(pars[2]))
  #dweibull(x, shape = exp(pars[1]), scale = exp(pars[2]))
}
#run optimiser
pars_optim<-optim(par = c(0,0), fn = neglogL2optim)
pars_optim_lognormal<-pars_optim

# shape of generation time distribution
w<-function(x,pars){
  #dlnorm(x, meanlog = pars[1], sdlog = exp(pars[2]))
  dgamma(x, shape = exp(pars[1]), rate = exp(pars[2]))
  #dweibull(x, shape = exp(pars[1]), scale = exp(pars[2]))
}
#run optimiser
pars_optim<-optim(par = c(0,0), fn = neglogL2optim)
pars_optim_gamma<-pars_optim

# shape of generation time distribution
w<-function(x,pars){
  #dlnorm(x, meanlog = pars[1], sdlog = exp(pars[2]))
  #dgamma(x, shape = exp(pars[1]), rate = exp(pars[2]))
  dweibull(x, shape = exp(pars[1]), scale = exp(pars[2]))
}
#run optimiser
pars_optim<-optim(par = c(0,0), fn = neglogL2optim)
pars_optim_weibull<-pars_optim

nlogL0<-neglogL2optim(pars_optim_weibull$par)
pars<-expand.grid(log(c(30:100)/20),log(c(40:75)/10))
vals<-apply(pars,1,neglogL2optim)
vals_ci<-as.data.frame(cbind(pars,vals)); colnames(vals_ci)<-c("par1","par2","nlogL")
vals_ci$DnlogL<-2*(vals_ci$nlogL-nlogL0)
vals_ci$shape<-exp(vals_ci$par1)
vals_ci$scale<-exp(vals_ci$par2)
library(lattice)
levelplot(DnlogL ~ shape * scale, data = vals_ci, cuts=99, main="parameters Weibull distribution, 95% CI ",col.regions = heat.colors(10000)[10001-(101-c(1:100))^2], panel = function(..., at, contour = FALSE, labels = NULL) {
  panel.levelplot(..., at = at, contour = contour, labels = labels)
  panel.contourplot(..., at = 6, contour = TRUE,labels = FALSE)
})

##pdf(file="generation_time_distribution-WeibullCI.pdf",height=5,width=5)
w<-function(x,pars){
  dweibull(x, shape = exp(pars[1]), scale = exp(pars[2]))
}
W<-sapply(0:(4*M),function(j){integrate(function(x){w(x,pars_optim_weibull$par)},lower = max(j-0.5,0),upper = j+0.5)$value})
plot(c(0:200)/10,sapply(c(0:200)/10,function(x){w(x,pars_optim_weibull$par)}),type="l",ylim=c(0,0.3),lwd=1,main="Distribution of generation times",ylab="Density",xlab="Generation time")
randv<-runif(dim(vals_ci)[1])
vals_w<-c(0:500)/1000
upperlower<-sapply(c(0:200)/10,function(y){
  DnlogL_w<-sapply(vals_w,function(vw){
    min(vals_ci$DnlogL[abs(apply(vals_ci,1,function(x){w(y,x[1:2])})-vw)<1/100])
  })
  return(c(min(vals_w[DnlogL_w<3.84],na.rm=T),max(vals_w[DnlogL_w<3.84],na.rm=T)))
})
lower<-upperlower[1,]
upper<-upperlower[2,]
polygon(c(c(0:200)/10,rev(c(0:200)/10)),c(upper,rev(lower)),col="gray")
lines(c(0:200)/10,sapply(c(0:200)/10,function(x){inc(x)}),col=2,lty=2,lwd=2)
lines(c(0:200)/10,sapply(c(0:200)/10,function(x){w(x,pars_optim_weibull$par)}),lwd=3,col=4)
legend("topright",c("generation time","incubation time"),lwd=c(3,2),col=c(4,2),lty=c(1,2))
##dev.off()


w<-function(x,pars){
  dlnorm(x, meanlog = pars[1], sdlog = exp(pars[2]))
}
W<-sapply(0:(4*M),function(j){integrate(function(x){w(x,pars_optim_lognormal$par)},lower = max(j-0.5,0),upper = j+0.5)$value})
plot(c(0:200)/10, sapply(c(0:200)/10, function(x){w(x,pars_optim_lognormal$par)}),
     type="l", ylim=c(0,0.25), lwd=1, main=NULL, #"Distribution of generation times",
     ylab="probability density",
     xlab=expression("time since infection " * tau * " (days)"), xaxs = "i",
     yaxs = "i")
lines(c(0:200)/10,sapply(c(0:200)/10,function(x){inc(x)}),col=2,lty=2,lwd=2)
lines(c(0:200)/10,sapply(c(0:200)/10,function(x){dlnorm(x,log(4),(sqrt(2*log(4.7/4))))}),col="grey",lwd=2)
lines(c(0:200)/10,sapply(c(0:200)/10,function(x){dlnorm(x,log(7.5)-(log((3.4/7.5)^2+1))/2,(sqrt(log((3.4/7.5)^2+1))))}),col="lightblue",lwd=2)
lines(c(0:200)/10,sapply(c(0:200)/10,function(x){w(x,pars_optim_lognormal$par)}),lwd=1)
w<-function(x,pars){
  dgamma(x, shape = exp(pars[1]), rate = exp(pars[2]))
}
W<-sapply(0:(4*M),function(j){integrate(function(x){w(x,pars_optim_gamma$par)},lower = max(j-0.5,0),upper = j+0.5)$value})
lines(c(0:200)/10,sapply(c(0:200)/10,function(x){w(x,pars_optim_gamma$par)}),lwd=2)
w<-function(x,pars){
  dweibull(x, shape = exp(pars[1]), scale = exp(pars[2]))
}
W<-sapply(0:(4*M),function(j){integrate(function(x){w(x,pars_optim_weibull$par)},lower = max(j-0.5,0),upper = j+0.5)$value})
lines(c(0:200)/10,sapply(c(0:200)/10,function(x){w(x,pars_optim_weibull$par)}),lwd=3)
legend("topright",c("Generation time: Weibull",
                    "Generation time: gamma",
                    "Generation time: lognormal",
                    "Serial interval: Li et al.",
                    "Serial interval: Nishiura et al.",
                    "Incubation period"),
       lwd=c(3,2,1,2,2,2),col=c(1,1,1,"lightblue","grey",2),lty=c(1,1,1,1,1,2))
##dev.copy2pdf(file="generation_time_distribution.pdf", height = 7.5, width = 7.5)


###########
# probability of pre-symptomatic transmission
###########

w<-function(x,pars){
  dlnorm(x, meanlog = pars[1], sdlog = exp(pars[2]))
}
W<-sapply(0:(4*M),function(j){integrate(function(x){w(x,pars_optim_lognormal$par)},lower = max(j-0.5,0),upper = j+0.5)$value})
hist(1-apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)}),breaks=10,col="grey",main="Probability of pre-symptomatic transmission",ylab="Count",xlab="Bayesian posterior probability") 
abline(v=mean(1-apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)})),col=2)
##dev.copy2pdf(file="presyntomatic_transmission_posterior_lognormal.pdf")

w<-function(x,pars){
  dgamma(x, shape = exp(pars[1]), rate = exp(pars[2]))
}
W<-sapply(0:(4*M),function(j){integrate(function(x){w(x,pars_optim_gamma$par)},lower = max(j-0.5,0),upper = j+0.5)$value})
hist(1-apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)}),breaks=10,col="grey",main="Probability of pre-symptomatic transmission",ylab="Count",xlab="Bayesian posterior probability") 
abline(v=mean(1-apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)})),col=2)
##dev.copy2pdf(file="presyntomatic_transmission_posterior_gamma.pdf")


w<-function(x,pars){
  dweibull(x, shape = exp(pars[1]), scale = exp(pars[2]))
}
W<-sapply(0:(4*M),function(j){integrate(function(x){w(x,pars_optim_weibull$par)},lower = max(j-0.5,0),upper = j+0.5)$value})
hist(1-apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)}),breaks=10,col="grey",xlab="Probability that transmission occurred before symptoms",ylab="number of transmission pairs",main="") 
abline(v=mean(1-apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)})),col=2)
##dev.copy2pdf(file="presyntomatic_transmission_posterior_weibull.pdf")

w<-function(x,pars){
  dweibull(x, shape = exp(pars[1]), scale = exp(pars[2]))
}
W<-sapply(0:(4*M),function(j){integrate(function(x){w(x,pars_optim_weibull$par)},lower = max(j-0.5,0),upper = j+0.5)$value})
hist(1-apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)})*0.75/(0.25+apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)})*0.5),breaks=10,col="grey",xlab="Probability that transmission occurred before symptoms",ylab="number of transmission pairs",main="prior 25% pre-symptomatic transmissions") 
abline(v=mean(1-apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)})*0.75/(0.25+apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)})*0.5)),col=2)
##dev.copy2pdf(file="presyntomatic_transmission_posterior_weibull_prior_75sym_25pre.pdf")

w<-function(x,pars){
  dweibull(x, shape = exp(pars[1]), scale = exp(pars[2]))
}
W<-sapply(0:(4*M),function(j){integrate(function(x){w(x,pars_optim_weibull$par)},lower = max(j-0.5,0),upper = j+0.5)$value})
hist(1-apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)})*0.25/(0.75-apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)})*0.5),breaks=10,col="grey",xlab="Probability that transmission occurred before symptoms",ylab="number of transmission pairs",main="prior 75% pre-symptomatic transmissions") 
abline(v=mean(1-apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)})*0.25/(0.75-apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)})*0.5)),col=2)
##dev.copy2pdf(file="presyntomatic_transmission_posterior_weibull_prior_25sym_75pre.pdf")


distr_prob_asymtrans<-c();
k<-0
for(i in round(10*exp(-vals_ci$DnlogL/2))){
  k<-k+1
  print(c(k,i))
  if(i>0){
    for(cc in 1:i){
    prob_asymtrans<-1-apply(info_data,1,function(x){fracSymptomatic(as.list(x),W)})
    W<-sapply(0:(4*M),function(j){integrate(function(x){w(x,unlist(pars[k,]))},lower = max(j-0.5,0),upper = j+0.5)$value})
    distr_prob_asymtrans<-c(distr_prob_asymtrans,mean(prob_asymtrans>runif(40)))
    }
  }
}
quantile(distr_prob_asymtrans,0.025)
quantile(distr_prob_asymtrans,0.975)
hist(distr_prob_asymtrans,breaks=(c((min(distr_prob_asymtrans)*40-1):(max(distr_prob_asymtrans)*40))+0.5)/40,col="darkgreen",main="Pre-symptomatic transmission uncertainty",xlab="Fraction of pre-symptomatic transmissions",freq=F)
distr_prob_asymtrans_df<-data.frame(pa=distr_prob_asymtrans)
library(ggplot2)
p<-ggplot(distr_prob_asymtrans_df) + 
  geom_histogram(aes(pa, y=..density..),binwidth=1/40) + 
  labs(x = expression("R"["P"]*"/(R"["P"]*"+R"["S"]*")"), y = "posterior density") + 
  theme_classic(base_size = 10) + 
  coord_cartesian(expand = F) + xlim(c(0, 0.9))
##ggsave("uncertainty_fraction_presymptomatic.pdf", p, height = 3, width =3.4)

