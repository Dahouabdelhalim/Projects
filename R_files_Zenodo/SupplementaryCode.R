#1. Function to randomize time series of diversification rates and return
#resulting diversity values and randomized diversification rates.
#Input time series of diversity values (not log-transformed).
#For time series of m intervals, diversity input consists of n=(m+1)
#values: the m values at the bottom interval boundaries, and the
#value at the top boundary of the final interval.
#Calculate diversification rates.
#Randomize diversification rates to construct new diversity trajectory.
#Simulated diversity history constrained to start and end at observed levels
#and to remain at or above one species. Trials falling below this level
#are discarded.
#This function is used to generate results presented in Figures 2, S1, S2, S3D, S4D, and S9.
randomize_1<-function(D)
{
 n<-length(D)
 D<-log(D)
 R<-diff(D) #diversification rates scaled per time interval
 success<-0
 while(success==0)
 {
  RR<-sample(R) #randomize rates
  DR<-numeric() #corresponding random diversity trajectory
  DR[1]<-D[1] #start at observed first value
  DR[2:n]<-DR[1]+cumsum(RR) #calculate
  if (all(DR>=0)) success<-1
 }
 return(list(D=exp(DR),R=RR))
}

#2. Function to randomize time series of speciation and
#extinction rates and return resulting diversity values and
#randomized rates.
#Input time series of diversity;
#speciation rates (p) extinction rates (q) per lineage-million years;
#and interval length (dt).
#For time series of m intervals, diversity input consists of n=(m+1)
#values: the m values at the bottom interval boundaries, and the
#value at the top boundary of the final interval.
#Simulated diversity history constrained to start and end at observed levels
#and to remain at or above one species. Trials falling below this level
#are discarded.
#This function is used to generate results presented in Figure 4.
randomize_2<-function(D,p,q,dt)
{
 n<-length(D)
 D<-log(D)
 m<-length(p) #equal to n-1
 success<-0
 while(success==0)
 {
  ii<-sample(1:m)
  PR<-p[ii]
  QR<-q[ii]
  RR<-(PR-QR) #net diversification rate
  DR<-numeric() #corresponding random diversity trajectory
  DR[1]<-D[1] #start at observed first value
  DR[2:n]<-DR[1]+cumsum(RR*dt)
  if (all(DR>=0)) success<-1
 }
 return(list(D=exp(DR),P=PR,Q=QR))
}

#3. Function to simulate diversification rates and diversity trajectory
#according to empirical AR1 model of diversification rates.
#Input is time series of diversity, as in randomize1(), above.
#This function is used to generate results presented in Figure S5.
simulate_AR1<-function(D)
{
 n<-length(D)
 m<-(n-1)
 D<-log(D)
 R<-diff(D) #diversification rates scaled per time interval
 R1<-R[1:(m-1)]
 R2<-R[2:m]
 g<-glm(R2~R1) #fit AR1 model
 b0<-g$coef[1] #intercept
 b1<-g$coef[2] #slope
 R.pred<-b0+b1*R1 #predicted diversification rates
 R.resid<-R2-R.pred #residuals
 sd.resid<-sd(R.resid) #standard deviation of residuals
 success<-0
 while(success==0)
 {
  RR<-rep(0,m)
  DR<-numeric() #random diversity trajectory
  DR[1]<-D[1] #start at observed first value
  for (i in 1:m)
  {
   if (i==1)
   {
    RR[i]<-sample(R,1) #sample first rate at random. Thereafter, use autoregressive model.
   }
   else
   {
    RR[i]<-b0+b1*RR[i-1]+rnorm(1,sd=sd.resid)
   }
  }
  DR[2:n]<-DR[1]+cumsum(RR)
  if (all(DR>=0)) success<-1
 }
 return(list(D=exp(DR),R=RR))
}

#4. Algorithm for constructing evolutionary lineages from nominal species.
#Requires the following input for each of n species:
#genus: vector of genus names
#name: vector of full taxon names
#FA: vector of first appearance in terms of ordered temporal levels, 1=oldest.
#LA: vector of level of last appearance.
#This function is used to generate data columns on evolutionary lineages and
#results presented in Figure S2J.
make_lineages<-function(genus,name,FA,LA)
{
 n<-length(name)
 tolerances<-c(0,1,-1,2,-2,3,-3) #tolerances: 0=exact match; 1=gap of one level; -1=overlap of one level, etc.
 for (tol in tolerances) #cycle through tolerance levels
 {
  i<-0
  while (i<n) #loop through until no more species join to form lineages
  {
   i<-i+1
   #find candidate descendants: same genus; different nominal species; gap/overlap equal to tolerance
   ii<-which(genus==genus[i]&name!=name[i]&FA-LA[i]==tol)
   if (length(ii)==0) #no putative descendants; go to next species on list
   {}
   else if (length(ii)==1) #exactly one putative descendant
   {
    #redefine range of lineage; LA is now the LA of the descendant
    LA[i]<-LA[ii]
    #rewrite vectors, omitting the descendant; combined lineage keeps name of ancestor
    FA<-FA[-ii]
    LA<-LA[-ii]
    name<-name[-ii]
    genus<-genus[-ii]
    n<-length(name)
    i<-0 #reset counter and restart loop for recursion
   }
   else #multiple candidate descendants
   {
    ix<-sample(ii,1) #pick one at random
    #redefine range of lineage; LA is now the LA of the descendant
    LA[i]<-LA[ix]
    #rewrite vectors, omitting the descendant; combined lineage keeps name of ancestor
    FA<-FA[-ix]
    LA<-LA[-ix]
    name<-name[-ix]
    genus<-genus[-ix]
    n<-length(name)
    i<-0 #reset counter and restart loop for recursion
   }
  }
 }
 return(list(genus=genus,name=name,FA=FA,LA=LA))
}

#5. Diversity-dependent diversification simulation.
#N0: starting diversity
#K: carrying capacity
#D: diversity
#p0: speciation rate at D=0
#q0: extinction rate at D=0
#nint: number of time intervals
#Dt: length of interval, in Myr
#nstep: number of time steps per interval; set to small value to obtain a reasonable approximation to continuous-time process
#a: strength of diversity-dependence of speciation rate
#b: strength of diversity-dependence of extinction rate
#lag: lag in number of intervals
#TO: time of origination of species
#TE: time of extinction of species
#This function is used to generate results of Figures S7 and S8.
dd<-function(N0,K,p0,q0,nint,Dt,nstep,a,b,lag)
{
 tmax<-nint*nstep #ultimate number of fine time steps
 dt<-Dt/nstep #length of fine time step
 success<-0
 while (success==0)
 {
  m<-ntot<-nextant<-N0 #initialize counters
  TO<-rep(-0.1,N0) #species extant at start assigned arbitrary time of origination before t=0
  TE<-rep(NA,N0)
  is.extant<-rep(TRUE,N0)
  t<-0
  D<-rep(0,nint) #diversity at start of interval
  for (iint in 1:nint)
  {
   D[iint]<-nextant<-sum(is.extant)
   if (iint>lag) #if lag intervals have passed, invoke diversity-dependence
   {
    p<-(p0-a*(D[(iint-lag)]/K))*dt #convert from rate per Myr to probability per fine time step
    q<-(q0+b*(D[(iint-lag)]/K))*dt
   }
   else #if lag intervals have not passed, set rates to equilibrial values
   {
    p<-(p0-a)*dt #convert from rate per Myr to probability per fine time step
    q<-(q0+b)*dt
   }
   if (nextant>0)
   {
    for (istep in 1:nstep)
    {
     t<-t+1
     for (i in 1:m)
     {
      if (is.extant[i])
      {
       xr<-runif(1)
       if (xr<p) #branch
       {
        ntot<-ntot+1
        TO[ntot]<-t
        TE[ntot]<-NA
        is.extant[ntot]<-TRUE
       }
       else if (xr>1-q) #extinct
       {
        TE[i]<-t
        is.extant[i]<-FALSE
       }
      }
     }
     m<-ntot
    }
   }
   else break;
  }
  nextant<-sum(is.extant)
  if (nextant>0) success<-1
 }
 ii<-which(is.na(TE)) #find lineages still extant at end
 TE[ii]<-tmax+1 #assign them arbitrary time of extinction after t=tmax
 TO<-TO*dt #convert from fine time steps to Myr
 TE<-TE*dt
 return(list(TO=TO,TE=TE))
}

#6. Function to degrade true duration into preserved stratigraphic range.
#TO: time of origination
#TE: time of extinction
#s: sampling rate per lineage-million-years
#FA: time of first preserved appearance (NA if species never preserved)
#LA: time of last preserved appearance (NA if species never preserved)
#As in function dd(), lower values correspond to older times.
#This function is used to generate results of Figures S7 and S8.
degrade<-function(TO,TE,s)
{
 FA<-LA<-NA
 dur<-TE-TO #duration
 n<-rpois(1,dur*s) #number of sampling levels drawn at random from Poisson distribution
 if (n>0)
 {
  levels<-runif(n,min=TO,max=TE) #sampling levels drawn at random from uniform distribution
  FA<-min(levels)
  LA<-max(levels)
 }
 return(list(FA=FA,LA=LA))
}

#7. Function to perturb level of first and last appearance.
#nlev: number of unique event levels; 1=oldest, nlev=youngest
#FA.lev: observed level of first appearance
#LA.lev: observed level of last appearance
#E: magnitude of perturbation, in number of levels. Perturbations drawn uniformly from integers on [-E,E].
#FA.pert: perturbed level of first appearance
#LA.pert: perturbed level of last appearance
#This function is used to generate results of Figure S9.
perturb<-function(FA.lev,LA.lev,E,nlev)
{
 success<-0
 while (success==0)
 {
  FA.pert<-FA.lev+sample(-E:E,1)
  LA.pert<-LA.lev+sample(-E:E,1)
  if (FA.pert<1) FA.pert<-1 #if FA below oldest observed level, adjust to oldest
  if (LA.pert>nlev) LA.pert<-nlev #if LA above youngest observed level, adjust to youngest
  if (LA.pert>=FA.pert) success<-1 #keep result if FA is not younger than LA
 }
 return(list(FA.pert=FA.pert,LA.pert=LA.pert))
}

#8. Example of implementation of randomization, using baseline analysis (Figure 2).
#Input diversity at base of successive 0.25-Myr intervals.
#Final value is diversity at top of last interval.
D<-c(18,19,19,22,18,16,16,16,15,16,19,22,22,24,24,26,26,24,25,25,30,32,31,32,35,45,49,55,66,72,75,75,78,82,89,81,78,73,71,73,71,70,69,75,75,71,68,74,76,84,91,94,91,88,93,99,94,91,90,93,91,87,85,73,61,59,57,58,61,58,59,62,65,64,67,74,76,75,68,65,62,57,49,44,44,47,46,43,44,45,53,50,51,56,57,54,59,55,56,67,68,79,85,80,80,77,79,69,62,60,61,56,54,55,61,62,65,66,66,66,73,71,58,55,58,58,57,53,54,55,55,67,79,81,76,66,54,48,41,34,25,20,29,35,36,44,46,44,47,54,42,46,48,56,56,61,61,65,67,67,69,67,68,65,58,43,48,51,55,77,71,76,62,60,48,51,47,46,43,43,43,42,52,44,44,43,37,32,32,32,41,41,23,20,20,21,21,21,21,22,25,25,29,30,26,25,18,18,18,19,14,16,24,27,28,31,36,52,53,52,47,46,37,37,34,26,27,21,18,19,18,16,13,15,14,13,14,16,12,13,9,7,6,4,4,4,4,4,5)

nrep<-10000 #number of randomizations
cor.rnd<-rep(NA,nrep) #correlation between diversity and diversification for randomized diversification rates
f<-0.3 #LOWESS smoothing span
lag<-4 #lag, in number of intervals
n<-length(D)-1 #number of intervals
LD<-log(D)[1:n] #log diversity
r<-diff(log(D)) #diversification rates
L<-lowess(1:n,LD,f=f) #LOWESS regression
D.res<-LD-L$y #residuals
cor.obs<-cor(D.res[1:(n-lag)],r[(lag+1):n],method="spearman") #observed correlation between residual diversity and diversification rate
for (i in 1:nrep)
{
 DR<-randomize_1(D)$D #diversity trajectory resulting from randomized diversification rates
 LD<-log(DR)[1:n] #log diversity
 r<-diff(log(DR)) #randomized diversification rates
 L<-lowess(1:n,LD,f=f) #LOWESS regression
 D.res<-LD-L$y #residuals
 cor.rnd[i]<-cor(D.res[1:(n-lag)],r[(lag+1):n],method="spearman") #correlation for randomization
}

d<-density(cor.rnd) #kernel density estimator of null correlations
x<-d$x
y<-d$y
xlim<-range(x,cor.obs)
ylim<-range(0,y)

plot(x,y,type="l",xlim=xlim,ylim=ylim) #distribution of randomized correlations
points(cor.obs,0,pch=15) #observed correlation









