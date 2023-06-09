### Here is the implementation of the numerical evaluation of the
### likelihood for the 2004 data presented in the supp mat of the
### Edwards et al 2007 paper

source("shifted_dists.R")

###
### Simulate some data for testing:

make.dat.2004<-function(N.samps=3000, params=0.25/3600,
                        dist="sexp", mi=30, ss=10, jmin=3, upperlim=NULL){

    rdist<-eval(as.name(paste("r", dist, sep="")))
    formals(rdist)[2] <- params[1]
    if(length(params) == 2) formals(rdist)[3] <- params[2]
    if(dist=="sexp" || dist=="sgamma") formals(rdist)$x0<-mi
    if(dist=="pareto") formals(rdist)$a<-mi
    
    f.full<-rdist(N.samps)

    starts<-runif(N.samps, 0, ss)
    
    obs<-trunc((f.full-(1-starts))/ss)

    if(!is.null(upperlim)) obs<-obs[which(obs<upperlim)]

    maxJ<-max(obs)-jmin
    dat<-data.frame(matrix(NA, nrow=maxJ, ncol=2))
    names(dat)<-c("J", "count")
    for(i in 1:maxJ){
        dat$J[i]<-i+jmin
        dat$count[i]<-length(which(obs==(i+jmin)))
    }

    return(dat)
}


nll.2004.num<-function(dat, params, dist="sexp", ss=10, mi=30, jmin=3,
                       unitconv=NULL, w.conv=1, test=FALSE, ...){

    print(paste("params = ", round(params, digits=4), sep=" "))

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


    js <- 1:(max(dat$J)+1)

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
        ## print(j)
        tt<-NULL
        tt<- -(j-1) * delCDFs[j-1] + (j + 1) * delCDFs[j] +
            ( xfxintegral[j-1] - xfxintegral[j] )/ss
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

    ptemp<- (jmin+1)*delCDFs[jmin] - xfxintegral[jmin]/ss
    if(ptemp>1 || ptemp<0) stop("ptemp out of bounds")

    if(sum(ll.temp!=0)){
        ll<- - sum(ll.temp) + sum(dat$count) * log(1 - ptemp) # last bit
                                        # corresponds
                                        # to
                                        # renormalisation
        if(test){
            print(paste("ll.temp =",  round(sum(ll.temp), 4),
                        "PR3.gam =", round(ptemp, 4),
                        "counts =", sum(dat$count),
                        "normalization = ", 
                        round(sum(dat$count) * log(1 - ptemp), 5),
                        sep=" "))
        }
        
    }else{
        ll<- - 100000 + sum(dat$count)*log(1-ptemp)
        print("sum of ll.temp==0")
    }
    
    return(ll)

}

######################################################################
######################################################################

if(FALSE){
### Exponential example
    set.seed(123)
    N.samps<-10000
    lambda<-0.05
    mi<-30
    ss<-10
    jmin<-3

    dat<-make.dat.2004(N.samps, params=lambda, dist="sexp",
                       mi=mi, ss=ss, jmin=jmin)
    

    l.trials<-seq(0.0001, 0.5, length=100)
    nll<-rep(NA, length(l.trials))
    for(i in 1:length(l.trials)){
        nll[i]<- -nll.2004.num(dat=dat, dist="sexp", params=l.trials[i],
                               mi=mi, ss=ss, jmin=jmin,
                               unitconv=NULL, w.conv=NULL)
    }


    params<-0.1
    out.num<-optim(params, nll.2004.num, dat=dat, method="Brent",
               lower=0, upper=0.3, dist="sexp", mi=mi, ss=ss, jmin=jmin)


    plot(l.trials, nll, type="l", xlim=c(0, 0.3))
    abline(v=out.num$par[1], col="green")
    abline(v=lambda, col=2)


}

######################################################################
######################################################################

if(FALSE){
### Gamma example
    set.seed(123)
    N.samps<-3000

    s<-0.4
    r<-1.2/3600
    mi<-30
    ss<-10
    jmin<-3
    dist<-"sgamma"

    dat<-make.dat.2004(N.samps, params=c(s,r), dist=dist,
                       mi=mi, ss=ss, jmin=jmin)
    

    s.trials<-seq(0.01, 2, length=32)
    r.trials<-seq(0.5, 5, length=30)
    nll<-matrix(NA, nrow=length(r.trials), ncol=length(s.trials))
    for(i in 1:length(s.trials)){
        for(j in 1:length(r.trials)){
            nll[j,i]<- - nll.2004.num(dat=dat, dist=dist,
                                      params=c(s.trials[i],r.trials[j]),
                                      unitconv=3600, w.conv=2,
                                      mi=mi, ss=ss, jmin=jmin)

        }
    }
    filled.contour(r.trials, s.trials, nll, xlab="r", ylab="s")#,
    
    params<-c(1, 2)
    
    out.num<-optim(params, nll.2004.num, dat=dat, method = "L-BFGS-B",
                    lower=c(0.01, 0.01), dist=dist, unitconv=3600, w.conv=2,
                    mi=mi, ss=ss, jmin=jmin,
                    control=list(fnscale=10000) ) ##maxit=1000, factr=1e-10
                    
    filled.contour(r.trials, s.trials, ll, xlab="r", ylab="s",
                   color=terrain.colors,
                   plot.axes = {points(c(out.num$par[2], params[2], r*3600),
                                       c(out.num$par[1], params[1], s),
                                       pch=c(10, 1, 8)); axis(1); axis(2) })


}

######################################################################
######################################################################

if(FALSE){
### Pareto example
    set.seed(123)
    N.samps<-3000

    lambda<-0.9
    mi<-30
    ss<-10
    jmin<-3
    dist<-"pareto"

    dat<-make.dat.2004(N.samps, params=lambda, dist=dist,
                       mi=mi, ss=ss, jmin=jmin)

    nll.2004.num(dat=dat, dist=dist, params=lambda,
                 mi=mi, ss=ss, jmin=jmin, unitconv=NULL, w.conv=NULL)   

    l.trials<-seq(0.1, 3, length=30)
    nll<-rep(NA, length(l.trials))
    for(i in 1:length(l.trials)){
        nll[i]<- -nll.2004.num(dat=dat, dist=dist, params=l.trials[i],
                               mi=mi, ss=ss, jmin=jmin,
                               unitconv=NULL, w.conv=NULL)
    }
    

    params<-0.5
    out.num<-optim(params, nll.2004.num, dat=dat, method="Brent",
               lower=0, upper=3, dist=dist, mi=mi, ss=ss, jmin=jmin)


    plot(l.trials, nll, type="l", xlim=c(0, 3))
    abline(v=out.num$par[1], col="green")
    abline(v=lambda, col=2)


}

if(FALSE){
    ## q-Exp example

    set.seed(1234)
    N.samps<-3000

    q<-1.25 ## shape
    lambda<-1 ## rate

    unitconv<-3600
    w.conv<-2
    
    mi<-30
    ss<-10
    jmin<-3
    dist<-"sqexp"

    dat<-make.dat.2004(N.samps, params=c(q,lambda/unitconv), dist=dist,
                       mi=mi, ss=ss, jmin=jmin, upperlim=24*3600/ss)

    nll.2004.num(dat=dat, dist=dist, params=c(q, lambda),
                 mi=mi, ss=ss, jmin=jmin, unitconv=unitconv, w.conv=w.conv)   

    l.trials<-seq(0.1, 2, length=20)
    q.trials<-seq(1, 1.9, length=20)
    nll<-matrix(NA, nrow=length(l.trials), ncol=length(q.trials))
    for(i in 1:length(l.trials)){
        for(j in 1:length(q.trials)){
            nll[j,i]<- - nll.2004.num(dat=dat, dist=dist,
                                      params=c(q.trials[i],l.trials[j]),
                                      unitconv=unitconv, w.conv=w.conv,
                                      mi=mi, ss=ss, jmin=jmin)

        }
    }
    filled.contour(l.trials, q.trials, nll, xlab="l", ylab="q")#,



    params<-c(1, 1)
    
    out.num<-optim(params, nll.2004.num, dat=dat, method = "L-BFGS-B",
                   lower=c(1, 0.01), upper=c(2, Inf), dist=dist,
                   unitconv=3600, w.conv=2,
                    mi=mi, ss=ss, jmin=jmin,
                    control=list(fnscale=10000) ) ##maxit=1000, factr=1e-10
                    
    filled.contour(l.trials, q.trials, nll, xlab="r", ylab="s",
                   color=terrain.colors,
                   plot.axes = {points(c(out.num$par[2], params[2], lambda),
                                       c(out.num$par[1], params[1], q),
                                       pch=c(10, 1, 8)); axis(1); axis(2) })




}
