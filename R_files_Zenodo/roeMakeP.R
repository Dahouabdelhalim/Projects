##abbreviations:
# bd=birth date
# bm=body mass
# M=Mother
# Y=Yearling
# Ad=Adult


minsize <- 5 # minimum size of a roe deer
maxsize <- 35 # maximum size of a roe deer
minbd <- 105 # minimum birth date of a roe deer
maxbd <- 165 # maximum birth date of a roe deer


if (TRUE) {  ## debugging
  nb=25 #number of size classes within the model
  nd=25 # number of birth date classes within the model
##  nb=10 #number of size classes within the model
##  nd=10 # number of birth date classes within the model  
} else {
  nb=60 #number of size classes within the model
  nd <- 60# number of birth date classes within the model
}
n.age <- 2 # number of age classes
## 'max.age <- 50 # maximum possible age a roe deer can live until

# boundary points b and mesh points y for size and birth date
Lb<-0.99*minsize; Ub<-1.01*maxsize
bb<-Lb+c(0:nb)*(Ub-Lb)/nb
yb<-0.5*(bb[1:nb]+bb[2:(nb+1)])

bd<-minbd+c(0:nd)*(maxbd-minbd)/nd
yd<-0.5*(bd[1:nd]+bd[2:(nd+1)])

# Part (II) ##############################################################
# Compute the kernel component functions from the fitted models

# growth kernel -- read in intercept (a) and slope (b,c) for the mean (a,b,c) components of the kernel and d for the variance 
gxy<-function(x,y,xbd,a,b,c,d) {
        sigmax2g <-d
        sigmax <- sqrt(sigmax2g)
        mux <- a+b*x+c*xbd
        fac1 <- sqrt(2*pi)*sigmax
        fac2 <- ((y-mux)^2)/(2*sigmax2g)
        return(exp(-fac2)/fac1)
}

# early survival (part of recruitment)
survfaon<-function(Mbm,bd,a,b,c){
      bd[bd<=132]=132
       nkids <- exp(a+b*Mbm+c*bd) 
        nkids <- nkids/(1+nkids)
       nkids[nkids>1]=1
        return(nkids)
}

# yearling survival.
survy<-function(bm,a,b){
  s <- exp(a+b*bm) 
  s <- s/(1+s)
  s[s>1]=1
  return(s)
}

# inheritance Body mass --
ibmxy<-function(x,y,xbd,a,b,c,d) {
        mux <- a+b*x+c*xbd
        sigmax2 <- d
        sigmax <- sqrt(sigmax2)
        fac1<-sqrt(2*pi)*sigmax
        fac2<-((y-mux)^2)/(2*sigmax2)
        f<-exp(-fac2)/fac1
        return(f)
}

# inheritance birth date -- 
ibdxy<-function(x,y,a,c,d) {
  mux <- a+c*x
  sigmax2 <- d
  sigmax <- sqrt(sigmax2)
  fac1<-sqrt(2*pi)*sigmax
  fac2<-((y-mux)^2)/(2*sigmax2)
  f<-exp(-fac2)/fac1
  return(f)
}

# Calculates the mean and variance of a quantity, typically a continuous character
means.and.vars <- function(z,n){
    me <- sum(z*n) / sum(n)
    va <- sum(z*z*n) / sum(n) - me^2
    return(matrix(c(me,va),c(2,1)))
}

# For large matrices it is faster to calculate the dominant eigenvalue and eigenvectors associated with it via iteration rather than using eigen.
get.eigen.stuff <- function(mat){ 
sz <- dim(mat)[1]
t.now <- runif(sz)
t.now <- t.now/sum(t.now)
t.next <- mat%*%t.now
t.next <- t.next/sum(t.next)
i <- 0
while (sum(abs(t.next-t.now))>0.0000001){
i <- i+1
t.now <- t.next
t.next <- mat%*%t.now
lambda <- sum(t.next)/sum(t.now)
t.next <- t.next/sum(t.next)
}
r.now <- runif(sz)
r.now <- r.now/sum(r.now)
r.next <- r.now%*%mat
r.next <- r.next/sum(r.next)
while (sum(abs(r.next-r.now))>0.0000001){
  r.now <- r.next
r.next <- r.now%*%mat
r.next <- r.next/sum(r.next)
}
return(list(lambda,t.next,r.next))
}

# Part (IV) ##############################################################

# parameter list

survf <- c(-0.23685, 0.24834,-0.04256)# intercept,slope_bm,slope_bd
gro<-c(19.485286,0.63290550,-0.04661111,2.319529)# intercept,slope_bm,slope_bd,variance
surv<-c(-4.2426,0.5303,0.824)# intercept,slope_bm, adult survival
inhbd<-c(112.9039, 0.1692, 21.610) # intercept, slope_bd,variance  
inhbm<-c(20.90399-1.15660,0.38507,-0.09642,1.4786+0.1175) # intercept,slope_bm,slope_bd,variance




##inhertiance birth date
Dd <- t(outer(yd,yd,ibdxy,inhbd[1],inhbd[2],inhbd[3]))
##example with 3 classes of birth date and body mass
#   M_bd1 M_bd2 M_bd3
# BD1 x   x     x
# BD2 x   x     x 
# BD3 x   x     x
Dd <- Dd/matrix(as.vector(apply(Dd,2,sum)),nrow=nd,ncol=nd,byrow=TRUE)
Dd <- ifelse(is.nan(Dd),0,Dd)


SFmat=matrix(0,ncol=nb,nrow=nd)# recruitment
#   M_bm1 M_bm2 M_bm3
# BD1 x   x     x
# BD2 x   x     x 
# BD3 x   x     x

Dbmat=matrix(0,ncol=nb,nrow=nb*nd) #inheritance body mass
#          M_bm1  M_bm2 M_bm3
#     Y_bm1  x     x     x
# BD1 Y_bm2  x     x     x
#     Y_bm3  x     x     x
#     Y_bm1  x     x     x
# BD2 Y_bm2  x     x     x
#     Y_bm3  x     x     x
#     Y_bm1  x     x     x
# BD3 Y_bm2  x     x     x
#     Y_bm3  x     x     x


SY=diag(survy(yb,surv[1],surv[2])) ##yearling survival

#   Y_bm1 Y_bm2 Y_bm3
# Ad_bm1 x   0     0
# Ad_bm2 0   x     0 
# Ad_bm3 0   0     x


Gro=matrix(0,ncol=nb,nrow=nb*nd) #growth yea-->ad
#          Y_bm1  Y_bm2 Y_bm3
#     Ad_bm1  x     x     x
# BD1 Ad_bm2  x     x     x
#     Ad_bm3  x     x     x
#     Ad_bm1  x     x     x
# BD2 Ad_bm2  x     x     x
#     Ad_bm3  x     x     x
#     Ad_bm1  x     x     x
# BD3 Ad_bm2  x     x     x
#     Ad_bm3  x     x     x


for (i in 1:nd){#for each birth date
  xbd=yd[i] 
  
  #early survival
  SFmat[i,]<- survfaon(yb,xbd,survf[1],survf[2],survf[3])
  
  #inheritance body mass
  Db <- t(outer(yb,yb,ibmxy,xbd,inhbm[1],inhbm[2],inhbm[3],inhbm[4]))
  Db <- Db/matrix(as.vector(apply(Db,2,sum)),nrow=nb,ncol=nb,byrow=TRUE)
  Db <- ifelse(is.nan(Db),0,Db)
  Dbmat[((i-1)*nb+1):(i*nb),]=Db
  
  #survival and growth yea-->ad
  G <- t(outer(yb,yb,gxy,xbd,gro[1],gro[2],gro[3],gro[4]))
  G <- G/matrix(as.vector(apply(G,2,sum)),nrow=nb,ncol=nb,byrow=TRUE)
  G <- ifelse(is.nan(G),0,G)
  Gro[((i-1)*nb+1):(i*nb),]= G %*% SY
  }



P = F <- matrix(0,nrow=n.age*nb*nd,ncol=n.age*nb*nd)
for (i in 1:nd){

  SFm=matrix(rep(SFmat[i,],nb),ncol=nb,nrow=nb,byrow=T)
  P[(nd*nb+(i-1)*nb+1):(i*nb+nd*nb),((i-1)*nb+1):(i*nb)]=Gro[((i-1)*nb+1):(i*nb),]
  for (c in 1:nd){
    F[((i-1)*nb+1):(i*nb),((c-1)*nb+1+nb*nd):(c*nb+nd*nb)]=Dbmat[((i-1)*nb+1):(i*nb),]*SFm*Dd[i,c]
    
  }}

P[(nb*nd+1):(nb*nd*2),(nb*nd+1):(nb*nd*2)]=diag(rep(surv[3],nb*nd))

