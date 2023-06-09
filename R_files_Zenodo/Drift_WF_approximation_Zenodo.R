#####Prepare data Eira####
eira=read.table("eira.csv", h=T, sep=";")
str(eira)
names(eira)

{
  #Neutral = 6:72
  sum1=apply(eira[which(eira$year==1),6:72],2,sum, na.rm=T)#Nb of reference alleles period 1
  tot1=apply(!is.na(eira[which(eira$year==1),6:72]),2,sum, na.rm=T)*2#Total number of alleles period 1
  sum2=apply(eira[which(eira$year==2),6:72],2,sum, na.rm=T)#Period 2
  tot2=apply(!is.na(eira[which(eira$year==2),6:72]),2,sum, na.rm=T)*2#Period 2
  sum3=apply(eira[which(eira$year==3),6:72],2,sum, na.rm=T)#Period 3
  tot3=apply(!is.na(eira[which(eira$year==3),6:72]),2,sum, na.rm=T)*2#Period 3

  #VIP
  vip=which(names(eira)=="vgll3"|names(eira)=="six6")
  
  sum4=apply(eira[which(eira$year==1),c(vip)],2,sum, na.rm=T)#Nb of E alleles period 1
  tot4=apply(!is.na(eira[which(eira$year==1),c(vip)]),2,sum, na.rm=T)*2#Total number of alleles period 1
  sum5=apply(eira[which(eira$year==2),c(vip)],2,sum, na.rm=T)#Period 2
  tot5=apply(!is.na(eira[which(eira$year==2),c(vip)]),2,sum, na.rm=T)*2#Period 2
  sum6=apply(eira[which(eira$year==3),c(vip)],2,sum, na.rm=T)#Period 3
  tot6=apply(!is.na(eira[which(eira$year==3),c(vip)]),2,sum, na.rm=T)*2#Period 3

  #Number of neutral SNP
  n=length(sum1)
  #Combine the period in one array
  A_num=array(c(sum1, sum2, sum3), c(n,3))
  tot_num=array(c(tot1, tot2, tot3), c(n,3))
  
  #Combine the period in one array for VIP
  A_num_VIP=array(c(sum4, sum5, sum6), c(2,3))
  tot_num_VIP=array(c(tot4, tot5, tot6), c(2,3))
  
  #List used by runjags
  neutral_data=list(n=n, A_num=A_num, tot_num=tot_num, A_num_VIP=A_num_VIP, tot_num_VIP=tot_num_VIP,
                    lim=c(0,1), lim_y=array(1, c(n,2)), A_plot=A_num[,2:3], tot_plot=tot_num[,2:3],n_y=3, n_VIP=2)
  
}

#####OR prepare data Stryn####
stryn=read.table("stryn.csv", h=T, sep=";")
str(stryn)
names(stryn)

{
  #Neutral = 6:71
  sum1=apply(stryn[which(stryn$year==1),6:71],2,sum, na.rm=T)
  tot1=apply(!is.na(stryn[which(stryn$year==1),6:71]),2,sum, na.rm=T)*2
  sum2=apply(stryn[which(stryn$year==2),6:71],2,sum, na.rm=T)
  tot2=apply(!is.na(stryn[which(stryn$year==2),6:71]),2,sum, na.rm=T)*2
  
  #VIP
  vip=which(names(stryn)=="vgll3"|names(stryn)=="six6")
  
  sum4=apply(stryn[which(stryn$year==1),c(vip)],2,sum, na.rm=T)
  tot4=apply(!is.na(stryn[which(stryn$year==1),c(vip)]),2,sum, na.rm=T)*2
  sum5=apply(stryn[which(stryn$year==2),c(vip)],2,sum, na.rm=T)
  tot5=apply(!is.na(stryn[which(stryn$year==2),c(vip)]),2,sum, na.rm=T)*2
 
  n=length(sum1)
  A_num=array(c(sum1, sum2), c(n,2))
  tot_num=array(c(tot1, tot2), c(n,2))
  
  A_num_VIP=array(c(sum4, sum5), c(2,2))
  tot_num_VIP=array(c(tot4, tot5), c(2,2))

  neutral_data=list(n=n, A_num=A_num, tot_num=tot_num, A_num_VIP=A_num_VIP, tot_num_VIP=tot_num_VIP, 
                    lim=c(0,1), lim_y=array(1, c(n,2)),A_plot=cbind(A_num[,2],NA), tot_plot=cbind(tot_num[,2],NA), n_y=2, n_VIP=2)
}


########Model######
#Take the initial and final allele frequencies of putatively neutral loci to estimate param = t/2N
#Estimate vgll3 and six6 allele frequency

{estimate_t2N="model{
for(i in 1:n){#for each putitavely neutral marker 
  for(y in 1:n_y){#Sampling error for each SNP and year
  #########To calculate drift
    A_num[i,y]~dbinom(P_neutral[i,y], tot_num[i,y])
  }
  
  #Prior
  P_neutral[i,1]~dbeta(1,1)#First year is just drawn from a beta
  
  for(y in 1:(n_y-1)){
    #Next years depend on first years + drift assuming no selection
    P_neutral[i,y+1]~dnorm(P_plot[i,y], Tau_diff[i,y])#Allele frequency changes following a normal distribution
  
    #Dinterval to force allele frequency to be between 0 and 1
    lim_y[i,y] ~ dinterval(P_neutral[i,y+1], lim) 
  
    #Variation in AF follows WF approximation with the normal distribution
    Tau_diff[i,y]=1/var_diff[i,y]

    #Var depends on 1/2N and initial allele frequencies
    var_diff[i,y]=param[y]*(P_plot[i,y]*(1-P_plot[i,y]))#initial allele frequency independently estimated
  
    #Difference in Allele frequency => biased toward mean expected changes under WF
    diff_neutral[i,y]=abs(P_neutral[i,y+1]-P_neutral[i,y])
  
    #Diff for plot
    plot1[i,y]=abs(P_plot[i,y+1]-P_plot[i,y])

  }
  
  ########Allele frequency estimated independently from WF modelling to plot it, from year 2
  for(y in 1:(n_y-1)){
    A_plot[i,y]~dbinom(P_plot[i,y+1], tot_plot[i,y])
    P_plot[i,y+1]~dbeta(1,1)
  }
  P_plot[i,1]=P_neutral[i,1]#Beta(1,1) too
  
}

#VIP SNPs
for(i in 1:n_VIP){
  for(y in 1:n_y){#VIP allele frequency independent over years
    A_num_VIP[i,y]~dbinom(P_neutral_VIP[i,y], tot_num_VIP[i,y])
    P_neutral_VIP[i,y]~dbeta(1,1)#Every year is just drawn from a beta
  }

  #Estimated drift
  for(y in 1:(n_y-1)){
    drift[i,y]~dnorm(0,1/var_drift[i,y])#Expected deviation from 0
    var_drift[i,y]=param[y]*(P_neutral_VIP[i,y]*(1-P_neutral_VIP[i,y]))
  
    drift2[i,y]=abs(drift[i,y])#Absolute values of deviations
  
    diff[i,y]=abs(P_neutral_VIP[i,y+1]-P_neutral_VIP[i,y])#Absolute observed changes
  
    sel[i,y]=diff[i,y]>drift2[i,y]#Proportion of change larger than observed
  }
}

#Estimated drift parameters
for(n in 1:n_y){
  param[n]~dunif(0.0001,1)#prior for the param parameter
}

kappa~dunif(1,5000)

}"}


###Run it####
#with rjags, runjags...
library(runjags)
library(rjags)

res_rev=run.jags(estimate_t2N, monitor=c("param", "P_neutral", "drift", "sel", 
                                         "P_neutral_VIP", "drift2", "diff", "kappa",
                                         "P_plot", "plot1", "diff_neutral"), 
                 data=neutral_data, adapt=30000, 
                 sample=10000, thin=40, burnin = 370000,method = "rjparallel", n.chains=2)

###Model diagnostic...
library(coda);library(runjags)
mcmc.info=as.mcmc.list(res_rev)
#Gelman diagnostic
gel.res=gelman.diag(mcmc.info,multivariate = F)
sum(gel.res$psrf[,1]>1.1, na.rm=T)

