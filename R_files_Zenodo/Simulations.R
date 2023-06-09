
# Function to run diet summary simulations:
# -----------------------------------------
diet.sim <- function(n.taxa, n.samps, n.seqs, avg.taxa, poo.cutoff, p.opt='decline',
  bias=rep(1,n.taxa), seed=1)
{
  # Inputs:
  # - n.taxa: number of taxa in the diet of the population
  # - n.samps: number of scat samples
  # - n.seqs: number of sequence reads obtained per sample
  # - avg.taxa: average number of taxa per sample
  # - poo.cutoff: minimum read threshold required to call a taxon present, as a proportion of n.seqs
  # - p.opt: either "decline", in which case taxa in the diet occur in exponentially declining abundance; or
  #          "equal", in which case taxa occur in equal proportion in the diet
  # - bias: a vector with n.taxa elements specifying the relative sequence recovery bias for each taxa
  # - seed: seed for the random number generated can be set to ensure results are reproducible
  #
  # Outputs:
  # - a matrix

  if(p.opt=='decline') p=(1:n.taxa)^(-1.2);p=p/sum(p)
  if(p.opt=='equal') p=rmultinom(1,100,rep(1/n.taxa,n.taxa))/100

  taxa.per.samp=rbinom(n.samps, n.taxa, avg.taxa/n.taxa)
  taxa.per.samp[taxa.per.samp==0]=1 
 
  which.taxa=lapply(taxa.per.samp, function(x) sample(n.taxa, x, prob=p))
  cnt.dat=matrix(0, n.taxa, n.samps)

  for(i in 1:n.samps) {
   cnt.dat[which.taxa[[i]], i]<- rmultinom(1, n.seqs, (p*bias)[which.taxa[[i]]])
  }

  # Calculate:
  # RRA = relative read abundance
  # FOO = frequency of occurrence (presence/absence)
  # POO = proportion of occurrence
  RRA.per.sample=cnt.dat*NA
  for(i in 1:n.samps) {
    RRA.per.sample[,i]=cnt.dat[,i]/sum(cnt.dat[,i])
  }
  RRA=apply(RRA.per.sample,1,mean)
  FOO=(cnt.dat/n.seqs)>poo.cutoff
  POO=apply(FOO,1,sum)/sum(FOO)

  output=cbind(RRA,POO,p,bias)
  return(output)
}


# R code to produce simulation results in paper:
# -----------------------------------------------
library(fields)
library(ecodist)

n.taxa=40
n.samps=100
n.seqs=1000
n.sims=1000

# Code to calculate the Bray-Curtis dissimilarity measure for "n.sims" simulated count datsets for 24 different cases:
# - 2 levels of average taxa per sample (3 or 20)
# - 3 bias scenarios (no bias, low bias, high bias)
# - 4 methods of diet determination: RRA, POO with a min read threshold of 0.01 (1%), 0.001 (0.1%) and 0.05 (5%)
# The output is a list with 2 elements -- one for each average taxa level -- containing a vector of length 12000
# with the B-C values for the 3 bias scenarios x 4 diet methods x 1000 simulations
bc<- list()
for(avg.taxa in c(3,20)){
  true.p <- diet.sim(n.taxa=n.taxa, n.samps=10000,n.seqs=n.seqs,avg.taxa=avg.taxa,poo.cutoff=0)[,1]
  RRA<- POO1 <- POO.1 <- POO5 <- array(NA,c(n.taxa,3,1000))
  for(i in 1:n.sims) {
     NoB=rep(1,n.taxa)
     LoB=exp(rnorm(n.taxa,0,log(4)/3))
     HiB=exp(rnorm(n.taxa,0,log(20)/3))
     for(b in 1:3){
       b.opt=c('NoB','LoB','HiB')[b]
       # Run simulations with poo.cutoff = 0.01  (note that the RRA results are the same regardless of the poo.cutoff)
       tmp=diet.sim(n.taxa=n.taxa, n.samps=n.samps,n.seqs=n.seqs,avg.taxa=avg.taxa,poo.cutoff=0.01,bias=eval(as.name(b.opt)))
       RRA[,b,i]=tmp[,1]
       POO1[,b,i]=tmp[,2]
       # Run simulations with poo.cutoff = 0.001
       POO.1[,b,i]=diet.sim(n.taxa=n.taxa, n.samps=n.samps,n.seqs=n.seqs,avg.taxa=avg.taxa,poo.cutoff=0.001,bias=eval(as.name(b.opt)))[,2]
       # Run simulations with poos.cutoff = 0.05
       POO5[,b,i]=diet.sim(n.taxa=n.taxa, n.samps=n.samps,n.seqs=n.seqs,avg.taxa=avg.taxa,poo.cutoff=0.05,bias=eval(as.name(b.opt)))[,2]
     }
  }
  out <- list(true.p,RRA,POO1,POO.1,POO5)
  # Use the bcdist function in R package 'ecodist' to calculate the Bray-Curtis dissimilarity for each bias scenario and diet method
  tmp<- NULL
  for(b in 1:3) {
    x=as.matrix(bcdist(t(cbind(out[[1]],out[[2]][,b,]))))[1,]
    tmp=c(tmp,x[-1])
    x=as.matrix(bcdist(t(cbind(out[[1]],out[[3]][,b,]))))[1,]
    tmp=c(tmp,x[-1])
    x=as.matrix(bcdist(t(cbind(out[[1]],out[[4]][,b,]))))[1,]
    tmp=c(tmp,x[-1])
    x=as.matrix(bcdist(t(cbind(out[[1]],out[[5]][,b,]))))[1,]
    tmp=c(tmp,x[-1])
  }
  meth=rep(rep(c("RRA","POO 1%","POO 0.1%","POO 5%"),each=n.sims),times=3)
  bias=rep(c("No","Low","High"),each=n.sims*4)
  names(tmp)<- paste(meth,bias)
  bc[[paste(avg.taxa)]]<- tmp
}


# Code to calculate proportional diet summaries for a dataset simulated using the following 5 bias scenarios:
# - B0: no bias
# - B1: low bias, with bias on taxa 1 fixed at 3SD's above the mean
# - B2: high bias, with bias on taxa 1 fixed at 3SD's above the mean
# - B3: low bias, with bias on taxa 1 fixed at 3SD's below the mean
# - B4: high bias, with bias on taxa 1 fixed at 3SD's below the mean
B0=rep(1,n.taxa)
set.seed(1); B1=c(4,exp(rnorm(n.taxa-1,0,log(4)/3)))
set.seed(2); B2=c(20,exp(rnorm(n.taxa-1,0,log(20)/3)))
set.seed(30); B3=c(.25,exp(rnorm(n.taxa-1,0,log(4)/3)))
set.seed(40); B4=c(1/20,exp(rnorm(n.taxa-1,0,log(20)/3)))
diet.summary <- list()
for(avg.taxa in c(3,20)) {
  true.p <- diet.sim(n.taxa=n.taxa, n.samps=10000,n.seqs=n.seqs,avg.taxa=avg.taxa,poo.cutoff=0)[,1]
  RRA<- POO1 <- POO.1<- POO5 <- matrix(NA,n.taxa, 5)
  i=0
  for(b in c('B0','B1','B2','B3','B4')){
    i=i+1
    tmp=diet.sim(n.taxa=n.taxa, n.samps=n.samps,n.seqs=n.seqs,avg.taxa=avg.taxa,poo.cutoff=0.01,bias=eval(as.name(b)))
    RRA[,i]=tmp[,1]
    POO1[,i]=tmp[,2]
    POO.1[,i]=diet.sim(n.taxa=n.taxa, n.samps=n.samps,n.seqs=n.seqs,avg.taxa=avg.taxa,poo.cutoff=0.001,bias=eval(as.name(b)))[,2]
    POO5[,i]=diet.sim(n.taxa=n.taxa, n.samps=n.samps,n.seqs=n.seqs,avg.taxa=avg.taxa,poo.cutoff=0.05,bias=eval(as.name(b)))[,2]
  }
  diet.summary[[paste(avg.taxa)]]<- cbind(true.p,RRA,POO1,POO.1,POO5)
  colnames(diet.summary[[paste(avg.taxa)]]) <- c('True','RRA B0','RRA B1','RRA B2','RRA B3','RRA B4','POO1 B0','POO1 B1','POO1 B2','POO1 B3','POO1 B4',
       'POO.1 B0','POO.1 B1','POO.1 B2','POO.1 B3','POO.1 B4', 'POO5 B0','POO5 B1','POO5 B2','POO5 B3','POO5 B4')
}



# Code to produce figures for Box 2 of main text:
# -------------------------------------------------

# Plot showing proportion of each taxa in diet, for the exponentially declining option
windows(height=5,width=7)
p=(1:n.taxa)^(-1.2);p=p/sum(p)
plot(1:n.taxa,p,type='b',col=4,pch=16,xlab="Food taxa indicator",ylab="Proportion")


# Plot bias vectors for different scenarios:
# a) random low & high bias
windows(height=5,width=7)
set.seed(1); LoB=exp(rnorm(n.taxa,0,log(4)/3))
set.seed(2); HiB=exp(rnorm(n.taxa,0,log(20)/3))
plot(log(LoB),type='p',pch=1,lwd=2,lty=1,col='deeppink',ylim=c(-3.2,3.2),ylab="Bias factor (log scale)",xlab="Taxa",yaxt='n')
axis(side=2,at=log(c(.05,.25,1,5,20)),labels=c(.05,.25,1,5,20))
points(log(HiB),type='p',pch=2,lwd=2,col='steelblue',lty=2)
legend('topright',c('Low','High'),pch=1:2,col=c('deeppink','steelblue'),pt.lwd=2,cex=1.1,horiz=T)
abline(h=0)
# b) random low & high bias with taxa 1 fixed at 3SD's above/below the mean
windows(height=5,width=7)
bias.names=c('No bias','Low T1+','High T1+','Low T1-','High T1-')
bias.abbr=c('NB','T1+','T1++','T1-','T1--')
plot(log(B1),type='p',pch=1,lwd=c(3,rep(2,length(B2)-1)),col=2,,cex=c(1.3,rep(1,length(B2)-1)),
  ylim=c(-3.2,3.2),ylab="Bias factor (log scale)",xlab="Taxa",yaxt='n')
axis(side=2,at=log(c(.05,.25,1,5,20)),labels=c(.05,.25,1,5,20))
points(log(B2),type='p',pch=2,lwd=c(3,rep(2,length(B2)-1)),col=3,cex=c(1.3,rep(1,length(B2)-1)))
points(log(B3),type='p',pch=3,lwd=c(3,rep(2,length(B2)-1)),col=5,cex=c(1.3,rep(1,length(B2)-1)))
points(log(B4),type='p',pch=4,lwd=c(3,rep(2,length(B2)-1)),col=4,cex=c(1.3,rep(1,length(B2)-1)))
legend('topright',bias.names[-1],pch=1:4,col=c(2,3,5,4),pt.lwd=2,cex=1.1,horiz=T)
abline(h=0)


# Histograms showing number of taxa per sample for one simulation using an average of 3 or 20 taxa per sample
for(avg.taxa in c(3,20)){
  windows(5,5)
  taxa.per.samp=rbinom(n.samps, n.taxa, avg.taxa/n.taxa)
  taxa.per.samp[taxa.per.samp==0]=1 
  hist(taxa.per.samp,xlim=c(0,30),breaks=seq(1,30,1),prob=F,xlab="# taxa per sample",ylab="Frequency",main="")
  if(avg.taxa==3) legend('topright','Mean = 3',bty='n')
  if(avg.taxa==20) legend('topleft','Mean = 20',bty='n')
}


# Bray-Curtis boxplots comparing RRA and POO with min read threshold of 1%, for no, low and high bias scenarios
library(fields)
for(avg.taxa in c(3,20)){
   windows(height=5,width=7)
   tf <- meth%in%c("RRA","POO 1%")
   boxplot(bc[[paste(avg.taxa)]][tf]~factor(bias[tf],levels=c("No","Low","High"))+
     factor(meth[tf],levels=c("RRA","POO 1%")),
     at=c(1:3,4:6+.5),ylim=c(0,0.7),col=rep(c(0,'deeppink','steelblue'),2),ylab="Bray-Curtis dissimilarity",
     names=rep(c("No bias","Low","High"),2), las=2,cex.axis=1.1)
   mtext('RRA',side=3,at=2,line=.5,cex=.8)
   mtext('POO 1%',side=3,at=5.5,line=.5,cex=.8)
}


# Diet summary barplots comparing RRA and POO with a min read threshold of 1%, for bias scenarios with extreme bias on Taxa 1
for(avg.taxa in c(3,20)){
  windows(height=5,width=7)
  tf<- substring(colnames(diet.summary[[paste(avg.taxa)]]),1,4) %in% c('True','RRA ','POO1')
  barplot(diet.summary[[paste(avg.taxa)]][,tf],space=c(0,rep(c(.5,0,0,0,0),2)),col=c(tim.colors(10),rep('grey',n.taxa-10)),
    names.arg=c('',rep(bias.names,2)),cex.names=1.1,las=2,mgp=c(3,.5,0),ylab="Taxa proportion")
  mtext('TRUE',side=3,at=0.4,line=1,cex=.9)
  mtext('RRA',side=3,at=4,line=1,cex=.9)
  mtext('POO 1%',side=3,at=9.5,line=1,cex=.9)
  mtext(round(as.matrix(bcdist(t(diet.summary[[paste(avg.taxa)]])))[1,tf][-1],2),side=3,at=c(2:6,7.5:11.5),line=0,cex=.8)
}



# Code to produce figures for supplementary material
# -----------------------------------------------------

# Bray-Curtis boxlots comparing POO with min read thresholds of 0.1%, 1% and 5%, for no, low and high bias scenarios
for(avg.taxa in c(3,20)){
   windows(height=5,width=7)
   tf <- meth%in%c("POO 0.1%","POO 1%","POO 5%")
   boxplot(bc[[paste(avg.taxa)]][tf]~factor(bias[tf],levels=c("No","Low","High"))+
     factor(meth[tf],levels=c("POO 0.1%","POO 1%","POO 5%")),
     at=c(1:3,4:6+.5,8:10),ylim=c(0,0.7),col=rep(c(0,'deeppink','steelblue'),3),ylab="Bray-Curtis dissimilarity",
     names=rep(c("No bias","Low","High"),3), las=2,cex.axis=1.1)
   mtext('POO 0.1%',side=3,at=2,line=.5,cex=.8)
   mtext('POO 1%',side=3,at=5.5,line=.5,cex=.8)
   mtext('POO 5%',side=3,at=9,line=.5,cex=.8)
}


# Diet summary barplots comparing POO with min read thresholds of 1%, 0.1% and 5%, for bias scenarios with extreme bias on Taxa 1
for(avg.taxa in c(3,20)){
  windows(height=5,width=9)
  tf<- substring(colnames(diet.summary[[paste(avg.taxa)]]),1,4) %in% c('True','POO1','POO.','POO5')
  barplot(diet.summary[[paste(avg.taxa)]][,tf],space=c(0,rep(c(.5,0,0,0,0),3)),col=c(tim.colors(10),rep('grey',n.taxa-10)),
    names.arg=c('',rep(bias.names,3)),cex.names=1.1,las=2,mgp=c(3,.5,0),ylab="Taxa proportion")
  mtext('TRUE',side=3,at=0.4,line=1,cex=.9)
  mtext('POO 1%',side=3,at=4,line=1,cex=.9)
  mtext('POO 0.1%',side=3,at=9.5,line=1,cex=.9)
  mtext('POO 5%',side=3,at=15,line=1,cex=.9)
  mtext(round(as.matrix(bcdist(t(diet.summary[[paste(avg.taxa)]])))[1,tf][-1],2),side=3,at=c(2:6,7.5:11.5,13:17),line=0,cex=.8)
}


