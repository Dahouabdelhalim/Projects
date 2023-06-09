################
# This code generates a grid of values
# for the maximum value of R0 for COVID-19
# that can be controlled by intervention (isolation+contact tracing)
# by solving the integral operator version of the Euler-Lotka equation
# (Ferretti, Wymant et al, Science 2020)
################

# incubation from Lauer et al
log_incubation_sd<-qnorm(0.975)-sqrt(qnorm(0.975)^2-2*(log(10.5)-log(5.5)))
log_incubation_median<-log(10.5)-log_incubation_sd*qnorm(0.975)
#S<-function(x){1-(1-fa)*plnorm(x,meanlog = log_incubation_median, sdlog = log_incubation_sd)}
# Michelle's exponential approximation
#S<-function(x){exp(-x/5.2)}
# our Weibull generation time
beta<-function(x){
  dweibull(x, shape = 2.826027, scale = 5.665302)
}

fa<-0.0625
fe<-0.1
delays<-c(0,4,12,24,48,72) # in hours

Yerror<-0.001
miniter<-3
maxiter<-100
nmax<-200
ndiscr<-10 # discretization of the unit for integral
n<-25 #discretization parameters

successmatrix<-list()
for (delay in delays){
# from here on, computes matrix of maxR0
S<-function(x){1-(1-fa)*plnorm(x-delay/24,meanlog = log_incubation_median, sdlog = log_incubation_sd)}
M1<-function(x){beta(x)*(1-ei+ei*S(x))}
M2<-function(x,y){(1-et+et*S(x+y)/S(y))}
v<-c(1:nmax)/ndiscr
Y<-rep(1,nmax)/nmax
m<-matrix(0,nrow=n-1,ncol=n-1)
for(ci in 1:(n-1)){
  for(cj in 1:(n-1)){
    ei<-ci/n ; et<-cj/n*(1-fe) 
    Y<-rep(1,nmax)/nmax;
    Yold<-Y
    for(i in 1:maxiter){
      Y<-M1(v)*sapply(v,function(z){sum(M2(z,v)*Y)/ndiscr})
      Y<-Y/sum(Y)
      if(sum(abs(Y-Yold))<Yerror & i>=miniter){break}
      Yold<-Y
    }
    m[ci,cj]<-(lm(I(M1(v)*sapply(v,function(z){sum(M2(z,v)*Y)/ndiscr})) ~ Y)$coeff[2])
  }
}
colnames(m)<-c(1:(n-1))/n
rownames(m)<-c(1:(n-1))/n
m<-1/m
successmatrix[[as.character(delay)]]<-m

}

for(name in names(successmatrix)){write.table(successmatrix[[name]],file=paste("maxR0controlled_matrix_ae_",name,".txt",sep=""),quote=F)}
