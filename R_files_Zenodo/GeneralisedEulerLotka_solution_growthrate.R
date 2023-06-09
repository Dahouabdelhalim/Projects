################
# This code generates a grid of values
# for the growth rate of COVID-19
# post intervention (isolation+contact tracing)
# by solving the integral operator version of the Euler-Lotka equation
# (Ferretti, Wymant et al, Science 2020)
################


# incubation from Lauer et al
log_incubation_sd<-qnorm(0.975)-sqrt(qnorm(0.975)^2-2*(log(10.5)-log(5.5)))
log_incubation_median<-log(10.5)-log_incubation_sd*qnorm(0.975)
#S<-function(x){1-plnorm(x,meanlog = log_incubation_median, sdlog = log_incubation_sd)}
# Michelle's exponential approximation
#S<-function(x){exp(-x/5.2)}
# our Weibull generation time
beta<-function(x){
  dweibull(x, shape = 2.826027, scale = 5.665302)
}

fa<-0.0625
fe<-0.1
delays<-c(0,4,12,24,48,72) # in hours
listR0<-c(1.7,2,2.5,3,3.5)

Yerror<-0.001
miniter<-3
maxiter<-100
nmax<-240
ndiscr<-10 # discretization of the unit for integral
n<-25 #discretization parameters

# generate matrices
r_matrix<-list()

for(delay in delays){
for(R0 in listR0){
print(c(R0,delay))

S<-function(x){1-(1-fa)*plnorm(x-delay/24,meanlog = log_incubation_median, sdlog = log_incubation_sd)}

# auxiliary functions 
M1<-function(x){R0*beta(x)*(1-ei+ei*S(x))}
M2<-function(x,y){(1-et+et*S(x+y)/S(y))}

# discretization of the integral
v<-c(1:nmax)/ndiscr

# initialization
Y<-rep(1,nmax)/nmax
r<-0

m<-matrix(0,nrow=n-1,ncol=n-1)
for(ci in 1:(n-1)){
  for(cj in 1:(n-1)){
    ei<-ci/n ; et<-cj/n*(1-fe)
    eigen<-function(my_r){
      r<-my_r
      Y<-rep(1,nmax)/nmax;
      Yold<-Y
    for(i in 1:maxiter){
      Y<-M1(v)*exp(-v*r)*sapply(v,function(z){sum(M2(z,v)*Y)/ndiscr})
      Y<-Y/sum(Y)
      if(sum(abs(Y-Yold))<Yerror & i>=miniter){break}
      Yold<-Y
    }
    return(lm(I(M1(v)*exp(-v*r)*sapply(v,function(z){sum(M2(z,v)*Y)/ndiscr})) ~ Y)$coeff[2]-1)
    }
    m[ci,cj]<-tryCatch(uniroot(eigen,interval = c(-2,2))$root, error = function(e) {return(NA)} ) }
}
colnames(m)<-c(1:(n-1))/n
rownames(m)<-c(1:(n-1))/n
r_matrix[[paste(as.character(R0),as.character(delay),sep="_")]]<-m

}
}

for(name in names(r_matrix)){write.table(r_matrix[[name]],file=paste("growthrate_matrix_ae_",name,".txt",sep=""),quote=F)}
