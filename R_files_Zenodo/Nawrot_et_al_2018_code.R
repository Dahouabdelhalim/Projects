# Electronic supplementary material
#
# Nawrot R., Scarponi D., Azzarone M., Dexter T.A., Kusnerik K.M., Wittmer J.M., Amorosi A., Kowalewski M. 2018. 
# Stratigraphic signatures of mass extinctions: ecological and sedimentary determinants. 
# Proceedings of the Royal Society B, 20181191 (doi: 10.1098/rspb.2018.1191)
#
# --------------------------------------------------------------------------------------------
# R code for analysing stratigraphic occurrences of species
# Rafal Nawrot, Florida Museum of Natural History
# email: rnawrot@flmnh.ufl.edu, paleo.nawrot@gmail.com
#
# Last updated: 10 August 2018
#
#
# Facies and sequence stratigraphic interpretations were based on Amorosi et al., 2017
# (doi.org/10.1016/j.marpetgeo.2017.01.020). Please contact me for additional R functions 
# used for plotting these data.
# 
# Note that the algorithm for estimating the number of extinction pulses (Wang & Zhong, 2018; 
# doi:10.1017/pab.2016.30), can take over an hour to compute depending on processing power.
#
# Before running the code dawnload Nawrot_et_al_Dataset_S1.zip and upack it into your working 
# directory.
#
#---------------------------------------------------------------------------------------------

op<-par(no.readonly = TRUE) #save default plotting parameters
library(scales)
library(vegan)

############# Upload data #############################################################
taxa<-read.csv('taxa.csv')
abund<-read.csv('abundances.csv')[,-(1:4)]
cores<-read.csv('abundances.csv')[,1:4]
rownames(abund)<-cores$Sample.ID

# correct the order of levels of Well.ID
wells<-c("204-S7", "205-S5", "205-S9", "205-S14") 
cores$Well.ID<-factor(cores$Well.ID, levels=wells)
nCores<-length(wells)

# species recorded in the Pleistocene strata
inPleis<-taxa$Taxon.ID[taxa$Pleistocene.Record==1]

# species abundant in brackish environments
brackish<-taxa$Taxon.ID[taxa$Brackish==1]

# positions of the Transgressive Surface and Maximum Flooding Surface in the cores 
MFS<-c(12.2, 11.4, 25.3, 26.3)
TS<-c(24, 22.4, 31.4, 31.7)

############# Split data by core ######################################################
# sampling depths per core
depth.pc<-split(cores$Depth.Sampling, cores$Well.ID) 

abund.pc<-split(abund, cores$Well.ID)
# keep only species actually present in a given core
abund.pc<-lapply(abund.pc, function(x) x[,colSums(x)>0])
# add sampling depths as rownames
for (i in 1:nCores){
  rownames(abund.pc[[i]])<-depth.pc[[i]]
}

# convert to presence-absence data
PA.pc<-lapply(abund.pc, function(x) ifelse(x>0,1,0))

# split taxonomic data by core
taxa.pc<-list()
for(i in 1:nCores){
  taxa.pc[[i]]<-droplevels(taxa[taxa$Taxon.ID %in% colnames(abund.pc[[i]]),])
}
names(taxa.pc)<-names(abund.pc)


############# First and last occurrences ##############################################
#----First occurrences (FOs)
calcFO<-function(x){ 
  out<-c()
  for (i in 1:ncol(x)) {
    out[i]<-as.numeric(rownames(x)[max(which(x[,i]>0))])
  }
  names(out)<-colnames(x)
  out
}

depthFO.pc<-lapply(abund.pc, calcFO)


#----Last occurrences (LOs)
calcLO<-function(x){ 
  out<-c()
  for (i in 1:ncol(x)) {
    out[i]<-as.numeric(rownames(x)[min(which(x[,i]>0), na.rm=T)])
  }
  names(out)<-colnames(x)
  out
}

depthLO.pc<-lapply(abund.pc, calcLO)

# number of LOs in each sample
nLOs.pc<-lapply(depthLO.pc, table)

# max number of LO per horizon in ech core
sapply(nLOs.pc,max)


#----Reorder taxa according to the positions of their FO and LO
PA.ord.pc<-PA.pc
depthFO.ord.pc<-depthFO.pc
depthLO.ord.pc<-depthLO.pc

for(i in 1:nCores){
  x<-PA.pc[[i]]
  lo<-depthLO.pc[[i]]
  fo<-depthFO.pc[[i]]
  x<-x[,order(fo, decreasing=T)]
  lo<-lo[order(fo, decreasing=T)]
  fo<-sort(fo, decreasing=T)
  PA.ord.pc[[i]]<-x[,order(lo, decreasing=T)]
  depthFO.ord.pc[[i]]<-fo[order(lo, decreasing=T)]
  depthLO.ord.pc[[i]]<-sort(lo, decreasing=T)
}


#----Plot the observed stratigraphic ranges
# (the plot will be saved as a pdf in your working directory)

cairo_pdf("Observed stratigraphic ranges.pdf")
par(mfrow=c(2,2), tcl=-0.2, mgp=c(1.5,0.3,0), mar=c(5,2.5,3,0.1), cex.axis=0.9)
xlim<-max(sapply(abund.pc, ncol))
for (i in 1:4){
  x<-PA.ord.pc[[i]]
  fo<-depthFO.ord.pc[[i]]
  lo<-depthLO.ord.pc[[i]]
  plot(1:ncol(x), type="n", ylim=c(37,-1), xlim=c(0, xlim), axes=F,
       main="", ylab="Stratigraphic position [m]", xlab="", bty="n")
  axis(2, at=seq(0,35, by=5))
  if(i %in% 1:2) { axis(1, at=seq(0,30, by=10)); m<-15 }
  else { axis(1, at=seq(0,80, by=10)); m<-40 }
  mtext("Species", side=1, line=1.5, at=m, cex=0.8)
  abline(h=0, lty=2)
  abline(h=MFS[i], lty=1, lwd=1.5, col="forestgreen")
  abline(h=TS[i], lty=1, lwd=1.5, col="cornflowerblue")
  mtext(wells[i], side=3, line=0.3, at=5, cex=1)
  
  #add ranges
  segments(1:ncol(x), fo, 1:ncol(x), lo, lwd=1, col=1) 
  
  #extend range to the Pleistocene
  p<-colnames(x) %in% inPleis
  segments((1:ncol(x))[p], fo[p], (1:ncol(x))[p], rep(38, ncol(x))[p], lwd=1, lty=1, col="gray60")
  
  #add occurrences
  for(j in 1:nrow(x)){
    d<-depth.pc[[i]][j]
    points((1:ncol(x))[x[j,]>0], rep(d, rowSums(x)[j]), pch=20, cex=0.3, col=1)
  }
  #add LOs
  points(1:ncol(x), lo, type="p", pch=19, cex=0.4, col=1)

}
par(op)
dev.off()

#----Plot the number of LOs observed at each sampled horizon and the preferred 
# water depth of species versus the stratigraphic position of their LO

cairo_pdf("Positions of LOs and preferred water depth.pdf")
# (A) number of LOs per horizon 
par(mfrow=c(2,4), tcl=-0.2, mgp=c(2,0.3,0), mar=c(2.5, 1.5, 9.5, 1.5), oma=c(0,2.2,0,0.2))
for (i in 1:nCores){
  x<-nLOs.pc[[i]]
  plot(1:5, type="n", xlim=c(0,11), ylim=c(35,0), bty="n", main="", ylab="", xlab="")
  if(i==1) {
    mtext("Stratigraphic position [m]", side=2, line=2, cex=0.7)
    mtext("A", side=2, at=-3, line=2, cex=1.4, las=2)
  }
  abline(h=0, lty=3)
  abline(h=MFS[i], lty=1, lwd=1.5, col="forestgreen")
  abline(h=TS[i], lty=1, lwd=1.5, col="cornflowerblue")
  d<-as.numeric(names(x))
  segments(0, d, x, d, lwd=2, col=1)
  mtext(wells[i], side=2, at=-2, line=-10, cex=0.9, las=2)
  mtext("Number of last occurrences", side=1, line=1.5, cex=0.65)
}

# (B) preferred water depth vs. LO
par(mar=c(9.5, 1.5, 2.5, 1.5))
xlm<-35
for (i in 1:4){
  y<-depthLO.pc[[i]]
  pd<-taxa.pc[[i]]$Pref.Depth
  y<-y[!is.na(pd)]  #species for which with pref. depth estimate is available 
  pd<-pd[!is.na(pd)]
  
  plot(pd, y, type="n", xlim=c(0,xlm), ylim=c(35,0), bty="n",
       main="", ylab="", xlab="")
  abline(h=0, lty=2)
  if(i==1) {
    mtext("Position of last occurrence [m]", side=2, line=2, cex=0.7)
    mtext("B", side=2, at=-3, line=2, cex=1.4, las=2)
  }
  rect(0,36,10,0, col="gray90", border=NA)  #rectangle to mark depths <10 m
  abline(h=MFS[i], lty=1, lwd=1.5, col="forestgreen")
  abline(h=TS[i], lty=1, lwd=1.5, col="cornflowerblue")
  points(pd, y, pch=21, bg="white", cex=0.9)
  
  #add outliers
  if(sum(pd>xlm)>0){
    arrows(x0=33, y0=y[pd>xlm], x1=37, y1=y[pd>xlm], length=0.06, angle=15)
    mtext(paste(round(pd[pd>xlm],1),"m"), side=4, line=0.08, at=y[pd>xlm], cex=0.6)
  }

  # species common in brackish environments 
  points(pd[names(y) %in% brackish], y[names(y) %in% brackish], pch=21, cex=0.9, bg="black")
  
  mtext(wells[i], side=2, at=-2, line=-10, cex=0.9, las=2)
  mtext("Preferred water depth [m]", side=1, line=1.5, cex=0.65)
}
par(op)
dev.off()


############# Null models #############################################################
#----Model 1
# Simulate multiple abundance matrices with the same number of samples and species 
# as in the original data under the assumption of random distribution of species 
# and uniform sampling intensity (constant number of individuals sampled per horizon).
#
# Parameters:
# m     abundance matrix (samples x taxa) with species names stored as column names
# nrep  number of simulations

model1<-function(m, nrep=10000){ 
  N<-nrow(m) #number of samples
  sampleSize<-floor(sum(m)/N)  #average sample size across all samples
  pool.tot<-colSums(m)  #total pool of species in a core
  pool.tot<-pool.tot[order(names(pool.tot))]  #to be sure that taxa are always in the same order
  
  #create 3D matrix to store simulated data
  sim<-array(0, dim=c(nrow(m), ncol(m), nrep))  
  colnames(sim)<-names(pool.tot)
  rownames(sim)<-rownames(m)
  
  for (a in 1:nrep){
    pool<-pool.tot
    #Go from bottom to top of the core and sample without replacement the same number
    #of specimens (= average sample size in the real data) from each horizon
    for (i in nrow(m):1){ 
      p<-rep(names(pool), pool)
      temp<-sample(p, size=sampleSize, replace=F)
      picked<-table(temp)
      sim[i,names(picked),a]<-picked
      #excluding selected specimens from the species pool
      pool[names(pool) %in% names(picked)]<-pool[names(pool) %in% names(picked)]-picked
    }
  }
  sim
}

#----Model 2
# Simulate multiple abundance matrices with the same number of samples and species 
# as in the original data under the assumption of random distribution of species, 
# while allowing the sample size to vary vertically according to the observed changes 
# in fossil abundance.
#
# Parameters:
# m     abundance matrix (samples x taxa) with species names stored as column names
# nrep  number of simulations

model2<-function(m, nrep=10000){ 
  sampleSizes<-rowSums(m)  #sample sizes
  pool.tot<-colSums(m)  #total pool of species in a core
  pool.tot<-pool.tot[order(names(pool.tot))]  #to be sure that taxa are always in the same order
  
  #create 3D matrix to store simulated data
  sim<-array(0, dim=c(nrow(m), ncol(m), nrep))
  colnames(sim)<-names(pool.tot)
  rownames(sim)<-rownames(m)
  
  for (a in 1:nrep){
    pool<-pool.tot
    # Go from bottom to top of the core and sampling without replacement 
    # the same number of specimens as in the original sample
    for (i in nrow(m):1){ 
      p<-rep(names(pool), pool)
      temp<-sample(p, size=sampleSizes[i], replace=F)
      picked<-table(temp)
      sim[i,names(picked),a]<-picked
      #excluding selected specimens from the species pool
      pool[names(pool) %in% names(picked)]<-pool[names(pool) %in% names(picked)]-picked
    }
  }
  sim
}

#---Simulate a bunch of abundance matrices
# can take a while for high number of iterations
res1.abund.pc<-lapply(abund.pc, model1, nrep=1000)
res2.abund.pc<-lapply(abund.pc, model2, nrep=1000)


#----Calculate the positions of LOs for each species based on simulated datasets
# (Note: because in model 1 the average sample size is rounded down to the nearest 
# integer, the total number of specimens in simulated datasets is lower compared to 
# the empirical data. As a result, rare species may not be sampled at all in some 
# iterations, prompting warnings when the minimum sampling depth is calculated 
# in the code below and producing NAs for them in the output.)

res1.depthLO.pc<-list()
for(i in 1:nCores){
  res1.depthLO.pc[[i]]<-apply(res1.abund.pc[[i]], 3, calcLO)
}
names(res1.depthLO.pc)<-wells

res2.depthLO.pc<-list()
for(i in 1:nCores){
  res2.depthLO.pc[[i]]<-apply(res2.abund.pc[[i]], 3, calcLO)
}
names(res2.depthLO.pc)<-wells


#----Distribution of LOs in simulated range charts ("ranked-LO distribution", RLOD)
calcRLOD<-function(x){ 
  temp<-apply(x, 2, sort, decreasing=T, na.last=T)
  out<-list(Median=NA, lCI=NA, uCI=NA)
  #calculate median and percentile-based 95% CIs across all simulations
  out$Median<-apply(temp, 1, median, na.rm=T)
  out$uCI<-apply(temp, 1, quantile, probs=0.975, na.rm=T)
  out$lCI<-apply(temp, 1, quantile, probs=0.025, na.rm=T)
  out
}

res1.RLOD.pc<-lapply(res1.depthLO.pc, calcRLOD)
res2.RLOD.pc<-lapply(res2.depthLO.pc, calcRLOD)


#----Number of LOs observed at each sampling horizon in the simulated data 
calc.nLOs<-function(x, depths){
  #create an array (samples in rows, iterations in columns) to store the results 
  los<-array(0, dim=c(length(depths), ncol(x)))
  rownames(los)<-depths #add sampling depths as rownames
  #tabulate numbers of LOs per sample and store them in correct positions
  temp<-apply(x,2,table)  
  for(j in 1:length(temp)){
    los[rownames(los) %in% names(temp[[j]]), j]<-temp[[j]]
  }
  #calculate median and percentile-based 95% CIs across all simulations
  out<-list(Median=NA, lCI=NA, uCI=NA)
  out$Median<-apply(los, 1, median)
  out$uCI<-apply(los, 1, quantile, probs=0.975)
  out$lCI<-apply(los, 1, quantile, probs=0.025)
  out
}
res1.nLOs.pc<-mapply(calc.nLOs, x=res1.depthLO.pc, depths=depth.pc, SIMPLIFY=F)
res2.nLOs.pc<-mapply(calc.nLOs, x=res2.depthLO.pc, depths=depth.pc, SIMPLIFY=F)


#----Plot the results of the resampling simulations
# fossil abundance (number of specimens) observed in each sample
sampleSize.pc<-lapply(abund.pc, rowSums)

for (i in 1:nCores){
  A<-wells[i]
  
  cairo_pdf(paste0("Resampling simulations for core ", A, ".pdf"), height=3.5)
  par(mfrow=c(1,3), tcl=-0.2, mgp=c(2,0.3,0), mar=c(3.5, 3, 3.5, 0), oma=c(1,4,1,4), cex.lab=1.2)
  
  # (A) stratigraphic trends in fossil abundance
  x<-sampleSize.pc[[A]]
  xlm<-max(sapply(sampleSize.pc, max))
  suppressWarnings(plot(x,names(x), type="o", xlim=c(0.2,xlm), ylim=c(35,0), xaxt="n", log="x", 
       pch=21, cex=0.8, bg="black", main="", ylab="Stratigraphic position [m]", 
       xlab="Fossil abundance", bty="n"))
  axis(1, at=c(1,5,50,500)); abline(h=0, lty=2)
  mtext("A", side=2, at=-2, line=-1, cex=1.3, las=2)
  mtext(wells[A], side=2, line=3.5, cex=1)
  #indicate position of barren samples
  points(x[x==0]+0.5, names(x)[x==0], pch=22, cex=0.8, bg="white")
  abline(h=MFS[i], lty=3); abline(h=TS[i], lty=3)
  legend("bottomright", leg=c("Barren sample"), pch=22, bty="n", cex=1.1)
  
  # (B) number of LOs per sample
  par(mar=c(3.5, 1.5, 3.5, 1.5))
  plot(1:5, type="n", xlim=c(0, 10), ylim=c(35,0), main="", yaxt="n",  
       ylab="", xlab="Number of last occurrences", bty="n")
  axis(2, at=seq(0, 35, by=5), lab=NA); abline(h=0, lty=2)
  mtext("B", side=2, at=-2, line=-1, cex=1.3, las=2)
  #add emipirical data
  x1<-c(nLOs.pc[[A]]) 
  d1<-as.numeric(names(x1))
  segments(0, d1, x1, d1, lwd=2, col="gray60")
  #add simulated data from model 2
  xmax<-res2.nLOs.pc[[A]]$uCI
  xmed<-res2.nLOs.pc[[A]]$Median
  xmin<-res2.nLOs.pc[[A]]$lCI
  ds<-as.numeric(names(xmax))
  points(xmed, ds, type="o", lwd=1, pch=19, cex=0.5, col="red") #median number of expected LO
  polygon(c(rev(xmax), xmin), c(rev(ds), ds), col=alpha("red",0.3), border=NA) #CIs shading
  #add simulated data from model 1
  xmax2<-res1.nLOs.pc[[A]]$uCI
  xmed2<-res1.nLOs.pc[[A]]$Median
  xmin2<-res1.nLOs.pc[[A]]$lCI
  ds<-as.numeric(names(xmax2))
  points(xmed2, ds, type="o", lwd=1, pch=19, cex=0.5, col="blue") #median number of expected LO
  polygon(c(rev(xmax2), xmin2), c(rev(ds), ds), col=alpha("blue",0.3), border=NA) #CIs shading
  #mark when the observed number of LOs is significantly different from the model predictions 
  lci<-mapply(min, xmin, xmin2) #lower CI for both models combined
  uci<-mapply(max, xmax, xmax2) #upper CI for both models combined
  lci<-lci[names(lci) %in% names(x1)]
  uci<-uci[names(uci) %in% names(x1)]
  points(x1[x1>uci | x1<lci], d1[x1>uci | x1<lci], pch=21, cex=0.8, bg=1)
  points(x1[x1<=uci & x1>=lci], d1[x1<=uci & x1>=lci], pch=21, cex=0.8, bg="white")
  abline(h=MFS[i], lty=3); abline(h=TS[i], lty=3)
  
  legend("bottomright", leg=c("Model 1", "Model 2"), lwd=1.5, 
         col=c("blue", "red"), bty="n", cex=1.1)
  
  # (C) ranked distribution of LOs (simulated range charts)
  par(mar=c(3.5, 0, 3.5, 3))
  x<-PA.ord.pc[[A]]
  fo<-depthFO.ord.pc[[A]]
  lo<-depthLO.ord.pc[[A]]
  
  plot(1:ncol(x), type="n", ylim=c(35,0), xlim=c(0, ncol(x)), main="", yaxt="n",
       ylab="", xlab="Species", bty="n")
  axis(2, at=seq(0,35, by=5), lab=NA); abline(h=0, lty=2)
  mtext("C", side=2, at=-2, line=-1, cex=1.3, las=2)
  #add simulated data from model 1
  sim<-res1.RLOD.pc[[A]]$Median
  ord<-1:length(sim)
  points(ord, sim, type="l", lwd=1.5, col="blue") #median ranked-LO
  cis<-c(res1.RLOD.pc[[A]]$lCI, rev(res1.RLOD.pc[[A]]$uCI))
  polygon(c(ord, rev(ord)), cis, col=alpha("blue",0.3), border=NA) #CIs shading
  #add simulated data from model 2
  sim<-res2.RLOD.pc[[A]]$Median
  ord<-1:length(sim)
  points(ord, sim, type="l", lwd=1.5, col=2)
  cis<-c(res2.RLOD.pc[[A]]$lCI, rev(res2.RLOD.pc[[A]]$uCI)) #median ranked-LO
  polygon(c(ord, rev(ord)), cis, col=alpha(2,0.3), border=NA) #CIs shading
  #add observed stratigraphic ranges and LOs
  segments(1:ncol(x), fo, 1:ncol(x), lo, lwd=1, col="gray20")
  p<-colnames(x) %in% inPleis
  segments((1:ncol(x))[p], fo[p], (1:ncol(x))[p], rep(36, ncol(x))[p], 
           lwd=1, lty=1, col="gray60")
  points(1:ncol(x), lo, type="o", pch=19, cex=0.5, lwd=1.5, col=1)
  abline(h=MFS[i], lty=3); abline(h=TS[i], lty=3)
  
  par(op)
  dev.off()
}


############# Stratigraphic abundance ################################################# 
# Stratigraphic abundance of Meldahl, 1990: proportion of stratigraphic intervals 
# (core samples in this case) in which a species was observed
stratab.pc<-lapply(PA.pc, FUN=function(x) colSums(x)/nrow(x))

#----Plot stratigraphic abundance of a species versus the position of its LO
cairo_pdf("Stratigraphic abundance vs. LO.pdf", height=3.5)
par(mfrow=c(1,4), tcl=-0.2, mgp=c(2,0.3,0), mar=c(6, 1.5, 6, 1.5), oma=c(0,2.5,0,1))
xlm<-max(sapply(stratab.pc, max)) #max start. abundance across all cores

for (i in 1:4){
  y<-depthLO.pc[[i]]
  x<-stratab.pc[[i]]
  plot(x,y, type="p", xlim=c(0,xlm), ylim=c(35,0), pch=21, bg="gray", cex=1,
       main="", ylab="", xlab="", bty="n")
  if(i==1) mtext("Position of last occurrence [m]", side=2, line=2, cex=0.8)
  abline(h=0, lty=2)
  abline(h=MFS[i], lty=1, lwd=1.5, col="forestgreen")
  abline(h=TS[i], lty=1, lwd=1.5, col="cornflowerblue")
  mtext(wells[i], side=3, at=0.55, line=0, cex=0.9)
  mtext("Strat. abundance", side=1, line=1.5, cex=0.7)
}
par(op)
dev.off()


############# Likelihood ratio test ################################################### 
# Likelihood ratio test of the null hypothesis of the simultaneous extinction of all species.
# See equations (2) and (3) in Wang & Everson, 2007.
# To make core data consistent with Wang & Everson's framework, the stratigraphic position 
# of all horizons (including the top of the succession) is expressed as a distance from 
# the base of the core

LRtest<-function(bottom, lo, pa, no.sing=T){
  #position of the LO of a species (measured from the bottom of a core)
  y<-bottom-lo 
  #common extinction time of all species if H0 is true (i.e the top of the core)
  t0<-bottom
  #number of fossil finds for each species
  n<-colSums(pa)
  #remove singletons
  if(no.sing){ 
    y<-y[n>=2]
    n<-n[n>=2]
  }
  #the total number of species
  S<-length(n)
  
  out<-rep(NA, 3)
  out[1]<-qchisq(p=0.95, df=2*S)  #critical chi-squared value
  out[2]<- -2*sum(n*log(y/t0))  #test statistic
  out[3]<-pchisq(out[2], df=2*S, lower.tail=F)  #p-value
  out
}

LRres<-data.frame(Core=wells, chisq.95=NA, chisq.stat=NA, p.val=NA)
for (i in 1:nCores){
  LRres[i,2:4]<-LRtest(bottom=max(depth.pc[[i]]), lo=depthLO.pc[[i]], pa=PA.pc[[i]], no.sing=T)
}
LRres


############# Classical confidence intervals on stratigraphic ranges ################## 
# The classical method of Strauss & Sadler, 1989.
# See equations (1) and (2) in Holland, 2003.
classicCI<-function(pa, ranges, cl=0.5){
  H<-colSums(pa) #number of horizons at which each taxon occurs based on P/A matrix
  a<- -1/(H-1)
  out<-(((1-cl)^a)-1)*ranges
  out
}

# length of the observed stratigraphic ranges (with species ordered according to FO and LO)
obsRange.ord.pc<-mapply(function(x,y) (x-y), x=depthFO.ord.pc, y=depthLO.ord.pc)

# length of the classical CIs (50% and 95% CIs)
CI50length.ord.pc<-mapply(classicCI, pa=PA.ord.pc, ranges=obsRange.ord.pc, cl=0.50)
CI95length.ord.pc<-mapply(classicCI, pa=PA.ord.pc, ranges=obsRange.ord.pc, cl=0.95)

# stratigraphic positions of CI endpoints (i.e. the core depth of LO minus the length of CI)
# negative values indicate that CIs extend beyond the modern sedimentary surface (= horizon 0 m)
CI50end.ord.pc<-mapply(function(x,y) x-y, x=depthLO.ord.pc, y=CI50length.ord.pc)
CI95end.ord.pc<-mapply(function(x,y) x-y, x=depthLO.ord.pc, y=CI95length.ord.pc)


############# Generalized confidence intervals on stratigraphic ranges ################ 
# Generalized confidence intervals of Marshall, 1997 with recovery potential
# of species estimated using DCA-based approach of Holland, 2003

#----Dataset preparation
# create occurrence matrix
PA<-ifelse(abund>0,1,0)

sum(rowSums(abund)==0)  #N of barren samples (across all cores)
sum(rowSums(PA)==1)  #N of samples with a single species
sum(colSums(PA)==1)  #N of species found in a single sample

# remove barren samples and samples with only 1 species
abund.dca<-abund[rowSums(abund)>0 & rowSums(PA)>1,]

# remove species occurring in a single sample
abund.dca<-abund.dca[,colSums(ifelse(abund.dca>0,1,0))>1]
dim(abund.dca) #143 samples and 83 species left

# restrict the data to samples and species used in DCA
cores.dca<-droplevels(cores[cores$Sample.ID %in% rownames(abund.dca), ])
taxa.dca<-droplevels(taxa[taxa$Taxon.ID %in% colnames(abund.dca), ])

# transform to presence-absence data
pa.dca<-ifelse(abund.dca>0,1,0) 

# split PA data by core and keep only species present in a given core (needed later)
pa.dca.pc<-split(as.data.frame(pa.dca), cores.dca$Well.ID)
pa.dca.pc<-lapply(pa.dca.pc, function(x) x[,colSums(x)>0])


#----DCA of the total dataset
prop.dca<-abund.dca/rowSums(abund.dca)  #transform data into proportional abundances
logprop.dca<-log(prop.dca+1)  #log-transform

DCA.out<-decorana(logprop.dca) 

cairo_pdf("DCA ordination of the total dataset.pdf")
plot(DCA.out, type="n")
points(DCA.out, pch=19, col="gray80", cex=0.9)
text(DCA.out, "species", lab=gsub("taxon.", "", taxa.dca$Taxon.ID), cex=0.8, col=2)
dev.off()

# extracting scores
sp.scores<-scores(DCA.out, "species", choices=1)
sample.scores<-scores(DCA.out, "sites", choices=1)

# correlation between the DCA axis 1 score of a species and the empirical estimate of
# the preferred water depth based on the ENEA data
wd<-!is.na(taxa.dca$Pref.Depth)  #species with water depth estimates
cor.test(sp.scores[wd], taxa.dca$Pref.Depth[wd])

#----Estimate parameters of species response curves
estResponseCurves<-function(species.scores, sample.scores, data) {
  taxa<-names(species.scores)
  out<-data.frame(Species=taxa, PD=NA, DT=NA, PA=NA)
  
  #Preferred Depth (PD): DCA axis 1 scores for species
  out$PD<-species.scores
  
  for(i in 1:length(taxa)){
    #Depth Tolerance (DT): standard deviation of axis 1 scores of all samples containing the species
    a<-which(data[,taxa[i]]>0)  
    out$DT[i]<-sd(sample.scores[a])
    #Peak Abundance (PA): proportion of samples in which the species was found located within
    #one standard deviation (DT) of the PD, rescaled by a factor of 1.186 
    b<-which(sample.scores<=out$PD[i]+out$DT[i] & sample.scores>=out$PD[i]-out$DT[i]) 
    ok<-intersect(a,b)
    out$PA[i]<-length(ok)/length(b)*1.186 #rescaling proportion of samples
  }
  out
}

sp.curves<-estResponseCurves(sp.scores, sample.scores, data=abund.dca)

#----Estimate sample-level collection probabilities (recovery potential) for each species
coll.prob<-array(NA, dim=dim(abund.dca))
colnames(coll.prob)<-colnames(abund.dca)
rownames(coll.prob)<-rownames(abund.dca)

for(i in 1:ncol(coll.prob)){
  coll.prob[,i]<-dnorm(sample.scores, sp.curves$PD[i], sp.curves$DT[i])*sp.curves$PA[i]
}

# use only species with >=4 occurrences
sp.minOccr<-colnames(pa.dca)[colSums(pa.dca)>=4]
coll.prob<-coll.prob[ ,sp.minOccr]

# number of species for which recovery potential could not be estimated
sum(is.na(coll.prob[1,])) + sum(colSums(coll.prob)==0, na.rm=T)

# split sampling depths by core
depth.dca.pc<-split(cores.dca$Depth.Sampling, cores.dca$Well.ID)

# split collection probabilities by core
coll.prob.pc<-list()
for(i in 1:nCores){
  coll.prob.pc[[i]]<-coll.prob[cores.dca$Well.ID==wells[i], ]
  rownames(coll.prob.pc[[i]])<-depth.dca.pc[[i]] # change names of samples to core depths
}
names(coll.prob.pc)<-wells


#----Plot patterns in DCA sample scores and recovery potential
# split sample scores by core
sample.scores.pc<-split(sample.scores, cores.dca$Well.ID)

#choose for which species recovery potential curves should be shown
sp1<-"Lentidium mediterraneum"
sp2<-"Varicorbula gibba"

cairo_pdf("Strat. trends in DCA scores and recovery potential.pdf")
par(mfrow=c(2,4), tcl=-0.2, mgp=c(2,0.3,0), mar=c(2.5, 1.5, 9.5, 1.5), oma=c(0,2.2,0,0.2))
# (A) stratigraphic trends in DCA sample scores
for (i in 1:4){
  plot(sample.scores.pc[[i]], depth.dca.pc[[i]], type="o", xlim=c(-3,3), ylim=c(35,0), 
       pch=21, bg="gray", cex=0.6, main="", ylab="", xlab="", bty="n")
  if(i==1){
    mtext("Stratigraphic position [m]", side=2, line=2, cex=0.8)
    mtext("A", side=2, at=-3, line=2, cex=1.4, las=2)
  }
  abline(h=0, lty=2)
  abline(h=MFS[i], lty=1, lwd=1.5, col="forestgreen")
  abline(h=TS[i], lty=1, lwd=1.5, col="cornflowerblue")
  mtext(wells[i], side=3, line=0, cex=0.8)
  mtext("DCA axis 1 scores", side=1, line=1.5, cex=0.7)
}

# (B) species recovery potential curves
par(mar=c(9.5, 1.5, 2.5, 1.5))
tax1<-as.character(taxa$Taxon.ID[taxa$Species==sp1])  #taxon.ID for species 1
tax2<-as.character(taxa$Taxon.ID[taxa$Species==sp2])  #taxon.ID for species 2
for (i in 1:4){
  plot(coll.prob.pc[[i]][,tax1], depth.dca.pc[[i]], type="o", xlim=c(-0.20,1), ylim=c(35,0), 
       xaxt="n", pch=21, col="red", cex=0.6, main="", ylab="", xlab="", bty="n")
  if(i==1){
    mtext("Stratigraphic position [m]", side=2, line=2, cex=0.8)
    mtext("B", side=2, at=-3, line=2, cex=1.4, las=2)
    legend("bottomleft", leg=c(sp1, sp2), pch=21, lty=1, col=c("red","blue"), cex=0.8, bty="n")
  }
  axis(1, at=seq(0,1,0.2))
  abline(h=MFS[i], lty=1, lwd=1.5, col="forestgreen")
  abline(h=TS[i], lty=1, lwd=1.5, col="cornflowerblue")
  abline(h=0, lty=2); abline(v=0)

  occ<-abund[cores$Well.ID==wells[i], tax1]>0  #samples in which species 1 was found
  d<-depth.pc[[i]][occ]
  points(rep(-0.15,length(d)), d, pch=21, bg="red", cex=0.6)
  
  points(coll.prob.pc[[i]][,tax2], depth.dca.pc[[i]], type="o", pch=21, col="blue", cex=0.6)
  occ<-abund[cores$Well.ID==wells[i], tax2]>0  #samples in which the species 2 was found
  d<-depth.pc[[i]][occ]
  points(rep(-0.05,length(d)), d, pch=21, bg="blue", cex=0.6)

  mtext(wells[i], side=3, line=0, cex=0.8)
  mtext("Recovery potential", side=1, at=0.5, line=1.5, cex=0.7)
}
par(op)
dev.off()

#----Calculate generalized CIs on stratigraphic ranges
dcaCI<-function(recovery, pa, cl){
  #match collection probabilities to taxa in the occurrence matrix 
  recovery<-recovery[, match(colnames(pa), colnames(recovery)) ] 
  colnames(recovery)<-colnames(pa)
  
  #calculate alpha from equation (2) in Holland, 2003
  H<-colSums(pa)
  a<- -1/(H-1)
  alpha<-((1-cl)^a)-1
  
  #calculate range extensions
  CIends<-rep(NA, ncol(recovery)) 
  for (i in 1:ncol(recovery)){
    #first, check if a species is a singleton or lacks recovery curve estimate
    if(!(sum(recovery[,i], na.rm=T)==0 | sum(is.na(recovery[,i]))>0 | H[i]==1)) {
      #calculate LO and FO
      first<-max(which(pa[,i]>0))
      last<-min(which(pa[,i]>0))
      #when species occurs in the topmost sample assume the CI extends beyond 
      #the top of a core and assign an arbitrary value of -1 (top = horizon 0 m)
      if(last==1){
        CIends[i]<- -1
      } else {
        #total recovery potential within the observed range for each species
        obsRec<-sum(recovery[first:last, i]) 
        #total recovery potential within the confidence limit on the range
        extRec<-obsRec*alpha[i]
        #when the total recovery potential in the rest of the core is lower 
        #than that in the observed range, the CI extends beyond the top of a core
        if(sum(recovery[(last-1):1, i])<extRec){
          CIends[i]<- -1
        } else {
          #cumulatively add recovery potentials for samples located above the observed LO
          r<-cumsum(recovery[(last-1):1, i])
          #position of the last sample for which the summed recovery potential is greater than
          #or equal to that recorded within the observed range marks the end of CI
          the.end<-min(which(r>=extRec))
          CIends[i]<-as.numeric(names(r)[the.end])
        }
      }
    }
  }
  
  names(CIends)<-colnames(recovery)
  CIends
}

# calculate CIs (negative values indicate that CIs extend beyond the top of a core)
dcaCI50end.pc<-mapply(dcaCI, recovery=coll.prob.pc, pa=pa.dca.pc, cl=0.5)
dcaCI95end.pc<-mapply(dcaCI, recovery=coll.prob.pc, pa=pa.dca.pc, cl=0.95)

# reorder CIs dataset according to FOs and LO in the total dataset (for plotting)
dcaCI50end.ord.pc<-list()
for( i in 1:nCores){
  x<-PA.ord.pc[[i]]
  y<-dcaCI50end.pc[[i]]
  a<-y[match(colnames(x), names(y))]
  names(a)<-colnames(x)
  dcaCI50end.ord.pc[[i]]<-a
}
names(dcaCI50end.ord.pc)<-wells

dcaCI95end.ord.pc<-list()
for( i in 1:nCores){
  x<-PA.ord.pc[[i]]
  y<-dcaCI95end.pc[[i]]
  a<-y[match(colnames(x), names(y))]
  names(a)<-colnames(x)
  dcaCI95end.ord.pc[[i]]<-a
}
names(dcaCI95end.ord.pc)<-wells


############# CIs on stratigraphic ranges: summary #################################### 
# proportion of taxa with =>4 occurrences in a given core for which the CIs do not reach
# the topmost sample (i.e. fall below 0.55 m core depth) 
fun<-function(pa, ci){
  ok<-ci[colSums(pa)>=4]  #taxa with =>4 occurrences
  n<-sum(!is.na(ok))  #number of taxa for which CIs could be estimated
  sum(ok>0.55, na.rm=T)/n 
}

# classical CIs
mapply(fun, pa=PA.ord.pc, ci=CI50end.ord.pc)  # 50% CIs
mapply(fun, pa=PA.ord.pc, ci=CI95end.ord.pc)  # 95% CIs

# ordination-based CIs
mapply(fun, pa=PA.ord.pc, ci=dcaCI50end.ord.pc)  # 50% CIs
mapply(fun, pa=PA.ord.pc, ci=dcaCI95end.ord.pc)  # 95% CIs


#----Plot CIs on stratigraphic ranges
# decide if ordination-based CIs should be used; otherwise classical CIs will be plotted
dcCIs<-T
# minimum number of occurrences per core
minOcc<-4 

cairo_pdf("Stratigraphic ranges with CIs.pdf")
par(mfrow=c(2,2), tcl=-0.2, mgp=c(1.5,0.3,0), mar=c(5,2.5,3,0.1), cex.axis=0.9)
for (i in 1:4){
  x<-PA.ord.pc[[i]]
  ok<-colSums(x)>=minOcc
  fo<-depthFO.ord.pc[[i]][ok] 
  lo<-depthLO.ord.pc[[i]][ok]
  if (dcCIs){
    CI50<-dcaCI50end.ord.pc[[i]][ok]
    CI95<-dcaCI95end.ord.pc[[i]][ok]
    col<-4
  } else{
    CI50<-CI50end.ord.pc[[i]][ok]
    CI95<-CI95end.ord.pc[[i]][ok]
    col<-2
  }
  x<-x[,ok] 

  plot(1:ncol(x), type="n", ylim=c(37,-1), xlim=c(0, 35), axes=F,
       main="", ylab="Stratigraphic position [m]", xlab="", bty="n")
  axis(2, at=seq(0,35, by=5))
  if(i %in% 1:2) { axis(1, at=c(0,4,8,12)); m<-6 }
  else { axis(1, at=seq(0,35, by=5)); m<-35/2 }
  mtext("Species", side=1, line=1.5, at=m, cex=0.8)
  abline(h=0, lty=2)
  mtext(wells[i], side=3, line=0.3, at=5, cex=1)
  
  #add the observed ranges
  segments(1:ncol(x), fo, 1:ncol(x), lo, lwd=1, col=1) 
  
  #extend ranges to the Pleistocene
  p<-colnames(x) %in% inPleis
  segments((1:ncol(x))[p], fo[p], (1:ncol(x))[p], rep(38, ncol(x))[p], lwd=1, lty=1, col="gray60")
  
  #add occurrences
  for(j in 1:nrow(x)){
    d<-depth.pc[[i]][j]
    points((1:ncol(x))[x[j,]>0], rep(d, rowSums(x)[j]), pch=20, cex=0.3, col=1)
  }

  #add CIs
  segments(1:ncol(x), lo, 1:ncol(x), CI95, lty=1, lwd=1, col=alpha(col, 0.5))
  if (dcCIs){
    #do not plot the end points of DCA-based CI ends if they extent beyond the section, 
    # as their exact position is not determined
    points((1:ncol(x))[CI95>0], CI95[CI95>0], pch=25, cex=0.4, col=col, bg=col)
    points((1:ncol(x))[CI50>0], CI50[CI50>0], pch=18, cex=0.7, col=col)
  } else {
    points((1:ncol(x)), CI95, pch=25, cex=0.4, col=col, bg=col)
    points((1:ncol(x)), CI50, pch=18, cex=0.7, col=col)
  }
  #adding LOs
  points(1:ncol(x), lo, type="p", pch=19, cex=0.4, col=1)
  
  if(i==2){legend("right", pch=c(25, 23, 21), pt.cex=0.6, cex=0.9, bty="n",
           col=c(alpha(col, 0.5), col, 1), pt.bg=c(col, col, 1),
           leg=c("95% Confidence interval", "50% Confidence interval", "Fossil occurrence"))
  }
}
par(op)
dev.off()


############# Estimating the number of extinction pulses ############################## 
# Two-step algorithm for estimating the number of extinction pulses of Wang & Zhong, 2018

# Source function definitions from the original R code provided in the supplementary
# material to Wang and Zhong’s paper (doi:10.1017/pab.2016.30)
source(textConnection(readLines('https://datadryad.org/bitstream/handle/10255/dryad.121812/code.R?sequence=1')[22:360]))

# Wang and Zhong’s code requires the functions knn.probability, knn.dist, and classprob
# from the package knnflex by Atina Dunlap Brooks 
# Source them from https://github.com/cran/knnflex (knnflex is no longer available on CRAN)
source("https://raw.githubusercontent.com/cran/knnflex/master/R/knn.probability.R")
source("https://raw.githubusercontent.com/cran/knnflex/master/R/knn.dist.R")
source("https://raw.githubusercontent.com/cran/knnflex/master/R/classprob.R")

# reformat the Po dataset and remove singletons
WZdataformat<-function(stratpos, data){
  out<-array(NA, dim=dim(data))
  for (i in 1:ncol(data)){
    occur<-stratpos[data[,i]>0]
    out[1:length(occur),i]<-occur
  }
  #let the strat. position of the samples be measured from the bottom to the top
  out<-40-out
  singles<-which(apply(out, 2, function(x) sum(!is.na(x))==1)) #remove singletons
  out<-out[,-singles]
  out
}

pulse.data<-mapply(WZdataformat, stratpos=depth.pc, data=abund.pc)

# Estimate the number of extinction pulses and corresponding posterior probabilities
# Depending on processing power, this can take over an hour to compute.
npulses.pc<-list()
for (i in 1:nCores){
  data<-pulse.data[[i]]
  
  # the following code was taken from Wang and Zhong’s script
  numtaxa <- dim(data)[2]
  ord <- order(data[1,])
  data <- data[,ord]           
  maxnumpulses <- 5 
  partitionlist <- vector("list", maxnumpulses)
  for(pulse in 1:maxnumpulses)  
    partitionlist[[pulse]] <- getpartitions(numtaxa, pulse) 
  pulseout <- howmanypulses(data=data, maxnumpulses=maxnumpulses)
  # end of Wang and Zhong’s code
  
  npulses.pc[[i]]<-pulseout
}

# plot result for each core
for (i in 1:nCores){
  cairo_pdf(paste0("ML estimates of extinction pulses in ", wells[i], ".pdf"), width=12, height=6)
  plotpulses(pulse.data[[i]], npulses.pc[[i]])
  dev.off()
}

############# END ####################################################################
