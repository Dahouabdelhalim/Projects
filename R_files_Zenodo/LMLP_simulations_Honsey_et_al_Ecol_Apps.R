## This file contains example code for running simulations to determine the data quality
## requirements and estimate reliability of the LMLP algorithm as described in the 
## Ecological Applications article by Honsey et al. entitled "Accurate estimates of 
## age-at-maturity from the growth trajectories of fishes and other ectotherms". In this 
## file, growth parameters and sample sizes-at-age are loosely based on walleye Sander vitreus
## (age-at-maturity = 5 yrs, immature growth intercept = 100 mm, immature growth slope = 50 mm/yr).
## See Honsey et al. (main text and supplemental file) for details.

## Load required libraries (install packages if needed)
library(compiler)
library(boot)

## Step 1: Create vectors containing values of each data quality factor

# In this example, we will have 5 iterations of the following combination:
# sample size = 250, precision in length-at-age = 14, and g = 0.2

# You can create whatever scenarios you like in R, or you could read in data
# quality factor values from a spreadsheet.

ss<-rep(250,5)
prec.set<-rep(14,5)
g.set<-rep(0.2,5)


## Step 2: Create vectors and matrices for storing values
parms.len<-length(ss) # store length of parameters (for looping)
AAM.est<-AAM.error<-AAM95lb<-AAM95ub<-g.est<-h1.est<-b0.est<-rep(NA,parms.len) # make empty vectors for storing results
Sims.Results<-cbind(ss,g.set,prec.set,AAM.est,AAM.error,AAM95lb,AAM95ub,g.est,h1.est,b0.est) # make matrix for storing results


## Step 3: Store and compile Lester model likelihood function as 'Biphas.Lik.MA' excluding age-at-maturity parameter 
## Optional: include marginal likelihoods for immature growth slope and intercept
###### likelihood function
Biphas.Lik.MA = function(parms) { 
  b0 = parms[1]
  h1 = parms[2]
  g = inv.logit(parms[3])   
  sighat = sqrt(parms[4]) 
  
  age.i = age[age<=mat.age]
  len.i = len[age<=mat.age]
  age.m = age[age>mat.age]
  len.m = len[age>mat.age]
  
  t1 = -b0/h1                                          
  Linf = 3*h1/g
  k = log(1 + g/3)
  t0 = mat.age + 
    suppressWarnings(log(1-(g*(mat.age-t1)/3)))/log(1+g/3)                                  
  mn.i = b0 + h1*age.i
  mn.m = Linf*(1-exp(-k*(age.m-t0)))
  
  b0.lik = dnorm(b0,mean=b0est,sd=25,log=T) #optional (also, distribution can be adjusted if needed)
  h1.lik = dnorm(h1,mean=h1est,sd=5, log=T) #optional (also, distribution can be adjusted if needed)
  L.i = dnorm(len.i,mean=mn.i,sd=sighat,log=T)
  L.m = dnorm(len.m,mean=mn.m,sd=sighat,log=T)
  return(sum(c(L.i,L.m, b0.lik,h1.lik)))
}
cmpfun(Biphas.Lik.MA) #compile function to run faster



## Step 4: Run loop to simulate populations for each iteration across data quality scenarios
for (i in 1:parms.len){
    
    ## Mean lengths at age based on biphasic model
    ## Parameters: h1=50, b0=100, age-at-maturity=5, varying g
    b.0<-100
    h.1<-50
    g.mat<-g.set[i]
    t.1 = -b.0/h.1                                          
    L.inf = 3*h.1/g.mat
    k.mat = log(1 + g.mat/3)
    t=5
    t.0 = t + suppressWarnings(log(1-(g.mat*(t-t.1)/3)))/log(1+g.mat/3)
    prec=prec.set[i]
    sample.size=ss[i]
    
    #Generate mean lengths-at-age for all fish (immature and mature)
    imm.list<-c(1:t) # immature ages
    mat.list<-c((t+1):25) # Mature ages
    imm.means<-sapply(imm.list,function(x) b.0+h.1*x) # generate mean lengths-at-age for immature fish
    mat.means<-sapply(mat.list,function(x) L.inf*(1-exp(-k.mat*(x-t.0)))) # same as above for mature fish
    means<-c(imm.means,mat.means) # concatenate mean lengths-at-age into one vector
    
    ## Generate data
    samp.size<-c(30,70,130,140,150,130,65,45,38,32,28,24,22,18,15,12,10,8,7,6,5,5,4,3,3) #set sample sizes for each age
    mean.samp<-as.data.frame(cbind(means,samp.size)) #make matrix of mean lengths-at-age and sample sizes
    mlen<-rep(mean.samp[,1],mean.samp[,2]) #repeat each mean "sample size" number of times
    a<-(1:25) #vector of ages
    ages<-rep(a,mean.samp[,2]) #repeat each age "sample size" number of times
    lengths<-sapply(mlen,function(x) rnorm(1,mean=x,sd=x/prec)) #generate random normal length data using means & precision
    Data<-as.data.frame(cbind(ages,lengths)) #bind vectors into age and length matrix, covert to data frame
    Data<-Data[sample(nrow(Data),size=sample.size,replace=F),] #draw a random sample from the population
    age = Data$ages #assign 'age' for likelihood function
    len = Data$lengths #assign 'len' for likelihood function
    
    ##Use linear model of first few years to inform marginal likelihoods for early growth slope & intercept
    immData<- Data[ which(age <= (min(age)+3)), ] #choose data within first four ages
    immout<-lm(lengths~ages, data=immData) #linear regression on "immature" data
    b0est<-immout$coefficients[[1]] #name intercept estimate, used for prior likelihood
    h1est<-immout$coefficients[[2]] #name slope estimate, used for prior likelihood
    
    ################ profile likelihood for mat.age
    ##############    Parameter  start values   
    b0 = b0est     #  early growth intercept
    h1 = h1est        # early growth slope 
    g = 0.15        #reproductive investment
    sighat = 25   #standard deviation
    parms=c(b0,h1,logit(g),sighat^2)  #compile parameters
    
    #  range of mat.age values for profile likelihood calculation
    Mat.age = seq(1,16,by=0.05)
    lik<-b0<-h1<-g<-var<-rep(NA,length(Mat.age))
    mat.age.Lik = cbind(Mat.age,lik,b0,h1,g,var) 
    
    ######  optimize Likelihood function for each mat.age value 
    for(j in 1:length(Mat.age)) {
      mat.age = Mat.age[j] 
      L.out = try(optim(par=parms,fn=Biphas.Lik.MA, 
                        control=list(fnscale=-1,reltol=1e-8)), silent=T)
      check<-is.numeric(L.out[[1]])
      
      ## store values only if model worked
      if (check[[1]] == "TRUE"){
        
        #Store parameter values (back-transform g)
        mat.age.Lik[j,2] <- L.out$value
        mat.age.Lik[j,3] <- L.out$par[[1]]
        mat.age.Lik[j,4] <- L.out$par[[2]]
        mat.age.Lik[j,5] <- inv.logit(L.out$par[[3]])
        mat.age.Lik[j,6] <- L.out$par[[4]]
      }}
    mat.age.Lik<-as.data.frame(mat.age.Lik) #convert to data frame for easier referencing
    mat.age.Lik<-mat.age.Lik[which(mat.age.Lik$lik != "NA"),] #remove failed runs
    mle = max(mat.age.Lik$lik) ## Find maximum likelihood
    MLE<-mat.age.Lik[which(mat.age.Lik$lik == mle),] ## maximum likelihood estimates for all parameters
    MLE.mat<-MLE[[1]] ## MLE for age-at-maturity
    H1<-MLE[[4]] ## MLE for immature growth slope
    B0<-MLE[[3]] ## MLE for immature growth intercept
    G<-MLE[[5]] ## MLE for investment in reproduction
    
    ## Confidence intervals in terms of chi-squared
    ndx1 = which(mat.age.Lik$lik>(mle-1.92)) # ~95% CI 
    BP.95CI = c(min(mat.age.Lik$Mat.age[ndx1]),max(mat.age.Lik$Mat.age[ndx1])) # determine 95% CI
    BP.95CI.LB<-BP.95CI[[1]]; BP.95CI.UB<-BP.95CI[[2]] # define 95% CI upper and lower bounds
    
    
    #Store results
    Sims.Results[i,4]<-MLE.mat # store LMLP estimate of AAM
    Sims.Results[i,5]<-abs(MLE.mat-t)/t #store % error in LMLP estimate of AAM
    Sims.Results[i,6]<-BP.95CI.LB # store lower bound of 95% confidence interval
    Sims.Results[i,7]<-BP.95CI.UB # store upper bound of 95% confidence interval
    Sims.Results[i,8]<-G # store estimate of cost to somatic growth of maturity
    Sims.Results[i,9]<-H1 # store estimate of immature growth slope
    Sims.Results[i,10]<-B0 # store estimate of immature growth intercept
    
  }
  
  ## Write results to a file
  Sims.Results<-as.data.frame(Sims.Results)
  write.csv(Sims.Results,file="Sims_Results.csv")
  