### Here is the implementation of the numerical evaluation of the
### likelihood for the 1992 data presented in the supp mat of the
### Edwards et al 2007 paper

source("shifted_dists.R")

###
### Simulate some data for testing:

make.dat.1992<-function(N.samps=3000, params=0.25, dist="sexp",
                        mi=1/120, ss=1, ret.obs=FALSE){

    rdist<-eval(as.name(paste("r", dist, sep="")))
    formals(rdist)[2] <- params[1]
    if(length(params) == 2) formals(rdist)[3] <- params[2]
    if(dist=="sexp" || dist=="sgamma") formals(rdist)$x0<-mi
    if(dist=="pareto") formals(rdist)$a<-mi

    f.full<-rdist(N.samps)

    starts<-runif(N.samps, 0, ss)
    
    obs<-trunc((f.full-(1-starts))/ss)

    maxJ<-max(obs)
    dat<-data.frame(matrix(NA, nrow=maxJ, ncol=2))
    names(dat)<-c("J", "count")
    for(i in 1:maxJ){
        dat$J[i]<-i
        dat$count[i]<-length(which(obs==i))
    }

    if(ret.obs) return(list(dat=dat, obs=obs))
    else return(dat)
}



### Negative log likelihood
nll.1992.num<-function(dat, params, dist="sexp", mi=1/120, ss=1, jmin=0, unitconv=NULL,
                       w.conv=1, test=FALSE){

    print(paste("params = ", round(params, digits=4), sep=" "))

    ##if(dist=="pareto") stop("pareto not supported in this function")

    if(!is.null(unitconv)) params[w.conv]<-params[w.conv]/unitconv

    pdist<-eval(as.name(paste("p", dist, sep="")))
    ddist<-eval(as.name(paste("d", dist, sep="")))

    formals(pdist)[2] <- params[1]
    formals(ddist)[2] <- params[1]
    if(length(params) == 2){
        formals(pdist)[3] <- params[2]
        formals(ddist)[3] <- params[2]
    }
    if(dist=="sexp" || dist=="sgamma"){
        formals(pdist)$x0<-mi
        formals(ddist)$x0<-mi
    }
    if(dist=="pareto"){
        formals(pdist)$a<-mi
        formals(ddist)$a<-mi
    }

    js <- 1:(max(dat$J)+2)

    CDFs <- sapply(ss*js, FUN=pdist) 
    delCDFs <- diff(CDFs)

    xfx<-function(x)(x*ddist(x))
    xfxintegral<-vector()
    for(i in 1:max(js)) xfxintegral[i]<-integrate(xfx, lower=(ss*i), upper=(ss*(i+1)))$value

    ll.temp<-NULL
    ll.temp<-rep(NA, length(dat$J))
    ll.temp[which(dat$count==0)]<-0
    
    ww<-which(dat$count!=0)
    for(i in 1:length(ww)){
        k<-ww[i]
        j<-dat$J[k]
        tt<-NA
        tt<- -j*delCDFs[j] + (2+j)*delCDFs[j+1]+(xfxintegral[j]-xfxintegral[j+1])/ss
        if(tt<0){
            if(abs(tt)<10^(-10)){
                ll.temp[k] <- dat$count[k] * log( abs(tt) )
                print("tt<=0, <<1")
            }else{
                ll.temp[k]<-dat$count[k] * -1000
                print(paste("tt<=0, j=", j, sep=""))
            }
        }else{
            ##print("correct evaluation")
            ll.temp[k] <- dat$count[k] * log( tt )
        }
    }

    ptemp<-NULL
    ptemp<-ddist(ss)+2*delCDFs[1]-xfxintegral[1]/ss
    if(ptemp>1 || ptemp<0) stop("ptemp out of bounds")
    
    if(sum(ll.temp!=0)){
        ll<- - sum(ll.temp) + sum(dat$count) * log(1 - ptemp) # last bit
    }else{
        ll<- - 100000 + sum(dat$count)*log(1-ptemp)
        print("sum of ll.temp==0")
    }
    
    return(ll)
}

#########################################################################
#########################################################################

### Examples of simulating data from the distributions and fitting the
### distributions to data

if(FALSE){
    ## EXPONENTIAL
    ## some trials
    N.samps<-3000
    mi<-1/120
    ss<-1
    lambda<-0.25
    unitconv<-NULL

    set.seed(1234)
    dat<-make.dat.1992(N.samps, params=lambda, dist="sexp", mi=mi, ss=ss)

    nll.1992.num(dat, params=lambda, dist="sexp", mi=mi, ss=ss, unitconv=3600, w.conv=1)

    l.trials<-seq(0.01, 2, length=100)
    lls<-rep(NA, length(l.trials))
    for(i in 1:length(l.trials)){
        lls[i]<- - nll.1992.num(dat, params=l.trials[i],  dist="sexp",
                                mi=mi, ss=ss, unitconv=unitconv, w.conv=1)
    }


    params<-0.1
    out.num<-optim(params, nll.1992.num, dat=dat, method="Brent",
                   lower=0, upper=1,  dist="sexp", mi=mi, ss=ss,
                   unitconv=unitconv, w.conv=1)

    plot(l.trials, lls, type="l")
    abline(v=out.num$par)
    abline(v=lambda, col=2)


}

#########################################################################
#########################################################################

    
if(FALSE){
    ## GAMMA
    ## some trials
    set.seed(1234)
    N.samps<-3000
    mi<-1/120

    s<-1.2
    r<-0.4

    dat<-make.dat.1992(N.samps, params=c(s,r), dist="sgamma", mi=mi)

    nll.1992.num(dat, params=c(s,r), dist="sgamma", mi=mi)

    s.trials<-seq(0.01, 4, length=50)
    r.trials<-seq(0.01, 3, length=50)
    ll<-matrix(NA, nrow=length(r.trials), ncol=length(s.trials))
    for(i in 1:length(s.trials)){
        for(j in 1:length(r.trials)){
            ll[j,i]<- -nll.1992.num(dat=dat, params=c(s.trials[i],r.trials[j]),
                                    dist="sgamma", mi=mi) 
        }
    }
    filled.contour(r.trials, s.trials, ll, xlab="r", ylab="s")#,
    


    params<-c(1,1)
    out.num<-optim(params, nll.1992.num, dat=dat, method = "L-BFGS-B",
                   lower=c(0.01, 0.01), upper=c(4, 3), dist="sgamma", mi=mi)
                   ##control=list(maxit=1000, factr=1e-12, fnscale=5000))


    filled.contour(r.trials, s.trials, ll, xlab="r", ylab="s",
                   plot.axes = {points(c(out.num$par[2], params[2], r),
                                       c(out.num$par[1], params[1], s),
                                       pch=c(10, 1, 8)); axis(1); axis(2) })
    


}


######################################################################
######################################################################

if(FALSE){
### Pareto 
    
    set.seed(123)
    N.samps<-3000

    lambda<-0.6
    mi<-30
    ss<-3600
    dist<-"pareto"

    out<-make.dat.1992(N.samps, params=lambda, dist=dist,
                       mi=mi, ss=ss, ret.obs=TRUE)
    dat<-out$dat

    nll.1992.num(dat=dat, dist=dist, params=lambda, mi=mi, ss=ss)   

    l.trials<-seq(0.1, 2, length=20)
    nll<-rep(NA, length(l.trials))
    for(i in 1:length(l.trials)){
        nll[i]<- -nll.1992.num(dat=dat, dist=dist, params=l.trials[i],
                               mi=mi, unitconv=NULL, w.conv=1, ss=ss)
    }
    

    params<-1
    out.num<-optim(params, nll.1992.num, dat=dat, method="Brent",
               lower=0, upper=3, dist=dist, mi=mi, unitconv=NULL, w.conv=1, ss=ss)


    plot(l.trials, nll, type="l", xlim=c(0, 2))
    abline(v=out.num$par[1], col="green")
    abline(v=lambda, col=2)


}
