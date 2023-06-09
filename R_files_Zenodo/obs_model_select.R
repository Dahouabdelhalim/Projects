source("mle_2004.R")
source("mle_1992.R")
source("shifted_dists.R")


## functions to calc AIC & BIC from optim output

myAIC<-function(out){
    2*length(out$par)+2*out$value
}

myBIC<-function(out, n){
    log(n)*length(out$par)+2*out$value
}

myBICMP<-function(BICs){

    bics<-unlist(BICs)

    ebics<-exp(-0.5*(bics-min(bics)))

    return(round(ebics/sum(ebics), 6))

}

my.interval.approx<-function(fit, neg.LL=TRUE){

    if(dim(fit$hessian)[1]==1){
        fisherinfo <-1/fit$hessian
    }else{
        if(neg.LL) fisherinfo<-solve(fit$hessian)
        if(!neg.LL) fisherinfo<-solve(-fit$hessian)
    }
    prop_sigma<-sqrt(diag(fisherinfo))
    ##if(dim(fit$hessian)[1]!=1)prop_sigma<-diag(prop_sigma)
    upper<-fit$par+1.96*prop_sigma
    lower<-fit$par-1.96*prop_sigma
    interval<-data.frame(value=fit$par, upper=upper, lower=lower)
    return(interval)
}

## read in and reconfigure data
w.dats<-c("BBA2002","walb2004","walb1998", "walb1992")


obs_model_selection_result <- obs.fitted.params<-list()

dat.all<-read.csv("albatross_drysegments.csv")

for (data_set in seq_along(w.dats)){
    w.dat <- w.dats[data_set]
    dat.raw<-subset(dat.all, dataset==w.dat)
    
    if(w.dat=="walb2004"){      
        ## parameters for the observation aggregation scheme for the
        ## particular dataset. set these when reading in the data.
        jmin<-3
        mi<-30
        ss<-10
        
        obs<-dat.raw$duration_n
        cutoff<-dat.raw$cutoff[1]
        l.fun<-nll.2004.num
    }
    
    if(w.dat=="walb1998"){
        ## parameters for the observation aggregation scheme for the
        ## particular dataset. set these when reading in the data.
        jmin<-2
        mi<-30
        ss<-15
        
        obs<-dat.raw$duration_n
        cutoff<-dat.raw$cutoff[1]
        l.fun<-nll.2004.num
        
    }
    
    if(w.dat=="BBA2002"){
        ## parameters for the observation aggregation scheme for the
        ## particular dataset. set these when reading in the data.
        jmin<-0
        mi<-30
        ss<-600
        
        obs<-dat.raw$duration_n
        cutoff<-dat.raw$cutoff[1]
        l.fun<-nll.1992.num
        
    }
    
    if(w.dat=="walb1992"){
        ## parameters for the observation aggregation scheme for the
        ## particular dataset. set these when reading in the data.
        jmin<-0
        mi<-30
        ss<-3600
        
        obs<-dat.raw$duration_n
        cutoff<-dat.raw$cutoff[1]
        l.fun<-nll.1992.num
        
    }
    
    
    obs<-obs[which(obs<cutoff)]
    l<-max(obs)
    dat<-data.frame(matrix(NA, ncol=2, nrow=l-jmin))
    names(dat)<-c("J", "count")
    
    dat$J<-(jmin+1):l
    
    for(i in 1:length(dat$J)){
        j<-dat$J[i]
        dat$count[i]<-length(which(obs==j))
    }
    
    n<-sum(dat$count)
    
    
    ## Fitting the four candidate distributions
    
    params<-1
    dist<-"sexp"
    unitconv<-3600
    w.conv<-1
    params<-1
    out.exp<-optim(params, l.fun, dat=dat, method="Brent",
                   lower=0, upper=2, dist=dist,
                   mi=mi, ss=ss, jmin=jmin,
                   unitconv=unitconv, w.conv=w.conv,
                   control=list(fnscale=1000),hessian=TRUE)
    
    dist<-"pareto"
    unitconv<-NULL
    out.par<-optim(params, l.fun, dat=dat, method="Brent",
                   lower=0, upper=2, dist=dist,
                   mi=mi, ss=ss, jmin=jmin,
                   unitconv=unitconv, w.conv=w.conv,
                   control=list(fnscale=1000),hessian=TRUE)
    
    
    dist<-"sgamma"
    unitconv<-3600
    w.conv<-2
    params<-c(1,1)
    out.gam<-optim(params, l.fun, dat=dat, method="L-BFGS-B",
                   lower=c(0.01, 0.01), upper=c(6, 5), dist=dist,
                   mi=mi, ss=ss, jmin=jmin,
                   unitconv=unitconv, w.conv=w.conv,
                   control=list(fnscale=1000),hessian=TRUE)
    
    
    dist<-"sqexp"
    unitconv<-3600
    w.conv<-2
    params<-c(1,1)
    if (w.dats[data_set]!="walb1992"){
        ##as optim with hessian=TRUE tries to return Hessian for
        ##unconstrained problem, it triggers a stop directive in
        ##psqexp(), so for now we skip interval estimation for the
        ##1992 data
        out.qexp<-optim(params, l.fun, dat=dat, method="L-BFGS-B",
                        lower=c(1, 0.01), upper=c(1.9999, Inf), dist=dist,
                        mi=mi, ss=ss, jmin=jmin,
                        unitconv=unitconv, w.conv=w.conv,
                        control=list(fnscale=1000),hessian=TRUE)
    } else {
        out.qexp<-optim(params, l.fun, dat=dat, method="L-BFGS-B",
                        lower=c(1, 0.01), upper=c(1.9999, Inf), dist=dist,
                        mi=mi, ss=ss, jmin=jmin,
                        unitconv=unitconv, w.conv=w.conv,
                    control=list(fnscale=1000),hessian=FALSE)
    }
    
    
    AICs<-list(exp=myAIC(out.exp),
               pareto=myAIC(out.par),
               gamma=myAIC(out.gam),
               qexp=myAIC(out.qexp))
    
    BICs<-list(exp=myBIC(out.exp, n),
               pareto=myBIC(out.par, n),
               gamma=myBIC(out.gam, n),
               qexp=myBIC(out.qexp, n))
    
    
    
    MPs<-myBICMP(BICs)  
    
    ##pack model selection results
    obs_model_selection_result[[data_set]] <- cbind(AICs, BICs, MPs)

    ##name subtables
    names(obs_model_selection_result)[data_set] <- w.dats[data_set]
    
    ##skip interval approximation for 1992 data (as there is no hessian)
    if (w.dats[data_set]!="walb1992"){
        obs.fitted.params[[data_set]]<-list(exp=my.interval.approx(out.exp),
                                            pareto=my.interval.approx(out.par),
                                            gamma=my.interval.approx(out.gam),
                                            qexp=my.interval.approx(out.qexp))
    }
    
    
    
}

obs_model_selection_result[[1]]
obs_model_selection_result[[2]]
obs_model_selection_result[[3]]
obs_model_selection_result[[4]]

##qq plots for the 3 main data sets:

par(mfrow=c(1,3), bty="n")

for (data_set in 1:3){
    w.dat <- w.dats[data_set]

    mi<-30
    cutoff<-24*3600 ## same as what we put into the data analysis
    dat.raw<-subset(dat.all, dataset==w.dat)

    obs<-dat.raw$duration_sec
    
    obs<-obs[(obs-mi)<cutoff]
    obs<-(obs-mi)/3600
    
    l<-length(obs)
    probs<-(1:l)/(l+1)
    
    shape<-obs.fitted.params[[data_set]]$gamma[1,1]
    rate<-obs.fitted.params[[data_set]]$gamma[2,1]
    
    gam.qs<- qgamma(probs, shape=shape, rate=rate)
    
    plot(sort(gam.qs), sort(obs),
         xlab = 'Theoretical Quantiles',
         ylab = 'Sample Quantiles',
         main =  w.dat)
    abline(0,1, col=2)   
}

