###### Parameter ######
args = commandArgs(trailingOnly=TRUE)
chrmax<-as.numeric(args[1])

###### Library ######

library(diversitree)
library(Matrix)
library(stats4)


###### Function, q_generation_2020 ######
### making Q-matrix from the 4 parameters of karytype evolution
### chrmax: maximum chromosome number

q_generation_2020<-function(k1,k2,k3,k4,chrmax){
  stanum<-(chrmax+1)*(chrmax+2)/2-1
  q<-matrix(rep(0),nrow=stanum,ncol=stanum)
  
  #        k2
  #        ^
  #        |
  # k4 <---+---> k3
  #        |
  #        v
  #        k1
  #
  
  #transition rates from states with limited directions of transition
  #in the beginning are defined one by one.
  
  q[1,2]<- k3
  q[1,1]<- -k3
  q[2,1]<- k4
  q[2,3]<- k2
  q[2,2]<- -k2 -k4
  q[3,2]<- k1 ##2020
  q[3,4]<- 2*k3
  q[3,3]<- -k1-2*k3 ##2020
  q[4,3]<- k4
  q[4,5]<- k3
  q[4,6]<- k2
  q[4,4]<- -k4-k3-k2
  q[5,4]<- 2*k4
  q[5,7]<- 2*k2
  q[5,5]<- -2*k4-2*k2
  
  for(d in 3:(chrmax-1)){
    #d, chromosome number
    
    #d*(d+1)/2 state with no metacentric (no fission, no M-A transition)
    q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
    q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1/2 ##2020
    q[d*(d+1)/2,d*(d+1)/2]<- (-d*k3 -(d*(d-1)*k1/2)) ##2020
    
    for(j in (d+1):(2*d-2)){
      #j, arm number
      #loop for the states with 4 directions
      
      #k,state number
      #up, upper state
      #down, lower state
      #anum, the number of acrocentric
      #mnum, the number of metacentric
      
      k<-d*(d+1)/2+j-d
      up<-(d+1)*(d+2)/2+j-d-1
      down<-d*(d-1)/2+j-d+1
      anum<-2*d-j
      mnum<-j-d
      q[k,k+1]<-anum*k3
      q[k,k-1]<-mnum*k4
      q[k,down]<-anum*(anum-1)*k1/2 ##2020
      q[k,up]<-mnum*k2
      q[k,k]<- -q[k,k+1]-q[k,k-1]-q[k,down]-q[k,up]
    }
    
    #states with 1 acrocentric (no fusion)
    k<-d*(d+1)/2+2*d-1-d
    up<-(d+1)*(d+2)/2+2*d-1-d-1
    q[k,k+1]<-k3
    q[k,k-1]<-(d-1)*k4
    q[k,up]<-(d-1)*k2
    q[k,k]<- -q[k,k+1]-q[k,k-1]-q[k,up]
    
    #states with no acrocentric (no fusion, no A-M transition)
    k<-d*(d+1)/2+2*d-d
    up<-(d+1)*(d+2)/2+2*d-d-1
    q[k,k-1]<-d*k4
    q[k,up]<-d*k2
    q[k,k]<- -q[k,k-1]-q[k,up]
    
  }
  
  d<-chrmax
  
  #state in the chrmax (no fission)
  #no metacentric case (no fission, no M-A transition)
  q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
  q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1/2 ##2020
  q[d*(d+1)/2,d*(d+1)/2]<- -d*k3-d*(d-1)*k1/2 ##2020
  
  for(j in (d+1):(2*d-2)){
    #3 directions
    k<-d*(d+1)/2+j-d
    down<-d*(d-1)/2+j-d+1
    anum<-2*d-j
    mnum<-j-d
    q[k,k+1]<-anum*k3
    q[k,k-1]<-mnum*k4
    q[k,down]<-anum*(anum-1)*k1/2 ##2020
    q[k,k]<- -q[k,k+1]-q[k,k-1]-q[k,down]
  }
  
  #one acrocentric case (no fission, no fusion)
  k<-d*(d+1)/2+2*d-1-d
  q[k,k+1]<-k3
  q[k,k-1]<-(d-1)*k4
  q[k,k]<- -q[k,k+1]-q[k,k-1]
  
  #no acrocentric case (no fission, no fusion, no A-M transition)
  k<-d*(d+1)/2+2*d-d
  q[k,k-1]<-d*k4
  q[k,k]<- -q[k,k-1]
  
  return(q)
}



###### Function for change between karyotype state vs. (arm no., chr no.) ######


#(chr number, arm number) -> state number
#this is (y,x) but not (x,y).
scal<-function(d,j){d*(d+1)/2+j-d}

#state number -> c(chr number, arm number)
chr_arm_cal<-function(s){
  d<- 1
  while(s-d*(d+1)/2 > d){
    d<- d+1
  }
  return(c(d,s-d*(d-1)/2))
}

#state vector -> matrix with two columns of chr number and arm number of each state
chr_arm_vec<-function(svec){
  chr_arm_mat<-c()  
  for(i in 1:length(svec)){
    chr_arm_mat<-rbind(chr_arm_mat,chr_arm_cal(svec[i]))
  }
  return(chr_arm_mat)
}

#state number -> c(arm number,chr number)
arm_chr_cal<-function(s){
  d<- 1
  while(s-d*(d+1)/2 > d){
    d<- d+1
  }
  return(c(s-d*(d-1)/2,d))
}

#state number -> karyotype "(x,y)"
karyotype_cal<-function(s){
  d<- 1
  while(s-d*(d+1)/2 > d){
    d<- d+1
  }
  return(paste("(",s-d*(d-1)/2,",",d,")",sep=""))
}

## Function, chr_arm_list
### vector of states -> list with [[1]] chr no. vector and [[2]] arm no. vector

chr_arm_list<-function(state){
  maxstate<-max(state)
  chr<-rep(0,length=length(state))
  arm<-rep(0,length=length(state))
  s<-1
  d<-1
  while(s <= maxstate){
    for(j in d:(2*d)){
      idents<-which(state==s)
      if(length(idents)>0){
        chr[idents]<-rep(d,length=length(idents))
        arm[idents]<-rep(j,length=length(idents))
      }
      s<-s+1
    }
    d<-d+1
  }
  return(list(chr,arm))
}


###### Functions for Mkn EC ######

##function, tcal
### Aquisition of number of transition rate 
### snum, total number of states
### is, initial state
### ts, transited state
### The no. of transition is
### snum x snum - snum (diagonals)
### left and right entries of diagonals have difference by one

tcal<-function(snum,is,ts){
  if(is<ts){
    (is-1)*(snum-1)+ts-1
  }else{
    (is-1)*(snum-1)+ts
  }
}

## function, qwrite
### Writing transition rate q like q0405
### The order of states are important.
### If you use few states, the (1,1) state are expressed as "q01".
### But if many, q001 or more longer.

qwrite<-function(ord,istate,tstate){
  ord<-as.character(sprintf("%02d",ord))
  qexp<-paste("q%",ord,"d%",ord,"d",sep="")
  sprintf(qexp,istate,tstate)
}

## function, get_cons_and_target
### important function to model karyotype evolution in mkn function
### as constraints of transition rates
### q003002,q002003,q001002 and q002001 are free parameters here,
### which corresponding to k1, k2, k3 and k4, respectively.
### The constraints are made with these 4 parameters.

get_cons_and_target_ec<-function(chrmax){
  
  snum<-(chrmax+1)*(chrmax+2)/2-1
  ord<-as.integer(log10(snum))+1
  ord2<-as.character(sprintf("%02d",ord))
  formula<-paste("q%",ord2,"d%",ord2,"d~%d*%s",sep="")
  kdel<-c(qwrite(ord,3,2),qwrite(ord,2,3),qwrite(ord,1,2),qwrite(ord,2,1))
  con_add<-function(cons,is,ts,kcoef,knum){
    return(c(cons,list(sprintf(formula,is,ts,kcoef,kdel[knum]))))
  }
  
  cons<-list()
  cons<-con_add(cons,3,4,2,3) #q[3,4]<- 2*k3
  cons<-con_add(cons,4,3,1,4) #q[4,3]<- k4
  cons<-con_add(cons,4,5,1,3) #q[4,5]<- k3
  cons<-con_add(cons,4,6,1,2) #q[4,6]<- k2
  cons<-con_add(cons,5,4,2,4) #q[5,4]<- 2*k4
  cons<-con_add(cons,5,7,2,2) #q[5,7]<- 2*k2
  target.i<-c(tcal(snum,3,4),tcal(snum,4,3),tcal(snum,4,5),tcal(snum,4,6),tcal(snum,5,4),tcal(snum,5,7))
  for(d in 3:(chrmax-1)){
    
    s<-scal(d,d)
    down<-scal(d-1,d)
    cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1
    cons<-con_add(cons,s,s+1,d,3) #q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
    target.i<-c(target.i,tcal(snum,s,down),tcal(snum,s,s+1))
    
    for(j in (d+1):(2*d-2)){
      s<-scal(d,j)
      up<-scal(d+1,j)
      down<-scal(d-1,j)
      anum<-2*d-j
      mnum<-j-d
      cons<-con_add(cons,s,down,anum*(anum-1)/2,1)#q[k,down]<-anum*(anum-1)*k1
      cons<-con_add(cons,s,s-1,mnum,4) #q[k,k-1]<-mnum*k4
      cons<-con_add(cons,s,s+1,anum,3) #q[k,k+1]<-anum*k3
      cons<-con_add(cons,s,up,mnum,2)#q[k,up]<-mnum*k2
      target.i<-c(target.i,tcal(snum,s,down),tcal(snum,s,s-1),tcal(snum,s,s+1),tcal(snum,s,up))
    }
    
    s<-scal(d,2*d-1)
    up<-scal(d+1,2*d-1)
    cons<-con_add(cons,s,s-1,d-1,4) #q[k,k-1]<-(d-1)*k4
    cons<-con_add(cons,s,s+1,1,3) #q[k,k+1]<-k3
    cons<-con_add(cons,s,up,d-1,2)#q[k,up]<-(d-1)*k2
    target.i<-c(target.i,tcal(snum,s,s-1),tcal(snum,s,s+1),tcal(snum,s,up))
    
    s<-scal(d,2*d)
    up<-scal(d+1,2*d-1)
    cons<-con_add(cons,s,s-1,d,4) #q[k,k-1]<-d*k4
    cons<-con_add(cons,s,up,d,2)#q[k,up]<-d*k2
    target.i<-c(target.i,tcal(snum,s,s-1),tcal(snum,s,up))
    
  }
  
  d<-chrmax
  
  s<-scal(d,d)
  down<-scal(d-1,d)
  cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1
  cons<-con_add(cons,s,s+1,d,3) #q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
  target.i<-c(target.i,tcal(snum,s,down),tcal(snum,s,s+1))
  
  for(j in (d+1):(2*d-2)){
    s<-scal(d,j)
    down<-scal(d-1,j)
    anum<-2*d-j
    mnum<-j-d
    cons<-con_add(cons,s,down,anum*(anum-1)/2,1)#q[k,down]<-anum*(anum-1)*k1
    cons<-con_add(cons,s,s-1,mnum,4) #q[k,k-1]<-mnum*k4
    cons<-con_add(cons,s,s+1,anum,3) #q[k,k+1]<-anum*k3
    target.i<-c(target.i,tcal(snum,s,down),tcal(snum,s,s-1),tcal(snum,s,s+1))
  }
  
  s<-scal(d,2*d-1)
  cons<-con_add(cons,s,s-1,d-1,4) #q[k,k-1]<-(d-1)*k4
  cons<-con_add(cons,s,s+1,1,3) #q[k,k+1]<-k3
  target.i<-c(target.i,tcal(snum,s,s-1),tcal(snum,s,s+1))
  
  s<-scal(d,2*d)
  cons<-con_add(cons,s,s-1,d,4) #q[k,k-1]<-d*k4
  target.i<-c(target.i,tcal(snum,s,s-1))
  
  return(list(cons,target.i))
}

## Function, get_cons_and_target_null
### For null hypothesis k3 = k4.
### Three parameters here.
### one constraint is added q002001 ~ q001002.

get_cons_and_target_null_ec<-function(chrmax){
  
  snum<-(chrmax+1)*(chrmax+2)/2-1
  ord<-as.integer(log10(snum))+1
  ord2<-as.character(sprintf("%02d",ord))
  formula<-paste("q%",ord2,"d%",ord2,"d~%d*%s",sep="")
  kdel<-c(qwrite(ord,3,2),qwrite(ord,2,3),qwrite(ord,1,2),qwrite(ord,1,2))
  con_add<-function(cons,is,ts,kcoef,knum){
    return(c(cons,list(sprintf(formula,is,ts,kcoef,kdel[knum]))))
  }
  
  cons<-list()
  cons<-con_add(cons,2,1,1,3) #q[2,1]<- q[1,2]
  cons<-con_add(cons,3,4,2,3) #q[3,4]<- 2*k3
  cons<-con_add(cons,4,3,1,4) #q[4,3]<- k4
  cons<-con_add(cons,4,5,1,3) #q[4,5]<- k3
  cons<-con_add(cons,4,6,1,2) #q[4,6]<- k2
  cons<-con_add(cons,5,4,2,4) #q[5,4]<- 2*k4
  cons<-con_add(cons,5,7,2,2) #q[5,7]<- 2*k2
  target.i<-c(tcal(snum,2,1),tcal(snum,3,4),tcal(snum,4,3),tcal(snum,4,5),tcal(snum,4,6),tcal(snum,5,4),tcal(snum,5,7))
  for(d in 3:(chrmax-1)){
    
    #All index numbers should be replaced 
    # d*(d-1) [sum:2d] -> d*(d+1)/2-1 [(sum:d+1)-1]
    
    s<-scal(d,d)
    down<-scal(d-1,d)
    cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1
    cons<-con_add(cons,s,s+1,d,3) #q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
    target.i<-c(target.i,tcal(snum,s,down),tcal(snum,s,s+1))
    
    for(j in (d+1):(2*d-2)){
      s<-scal(d,j)
      up<-scal(d+1,j)
      down<-scal(d-1,j)
      anum<-2*d-j
      mnum<-j-d
      cons<-con_add(cons,s,down,anum*(anum-1)/2,1)#q[k,down]<-anum*(anum-1)*k1
      cons<-con_add(cons,s,s-1,mnum,4) #q[k,k-1]<-mnum*k4
      cons<-con_add(cons,s,s+1,anum,3) #q[k,k+1]<-anum*k3
      cons<-con_add(cons,s,up,mnum,2)#q[k,up]<-mnum*k2
      target.i<-c(target.i,tcal(snum,s,down),tcal(snum,s,s-1),tcal(snum,s,s+1),tcal(snum,s,up))
    }
    
    s<-scal(d,2*d-1)
    up<-scal(d+1,2*d-1)
    cons<-con_add(cons,s,s-1,d-1,4) #q[k,k-1]<-(d-1)*k4
    cons<-con_add(cons,s,s+1,1,3) #q[k,k+1]<-k3
    cons<-con_add(cons,s,up,d-1,2)#q[k,up]<-(d-1)*k2
    target.i<-c(target.i,tcal(snum,s,s-1),tcal(snum,s,s+1),tcal(snum,s,up))
    
    s<-scal(d,2*d)
    up<-scal(d+1,2*d-1)
    cons<-con_add(cons,s,s-1,d,4) #q[k,k-1]<-d*k4
    cons<-con_add(cons,s,up,d,2)#q[k,up]<-d*k2
    target.i<-c(target.i,tcal(snum,s,s-1),tcal(snum,s,up))
    
  }
  
  d<-chrmax
  
  s<-scal(d,d)
  down<-scal(d-1,d)
  cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1
  cons<-con_add(cons,s,s+1,d,3) #q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
  target.i<-c(target.i,tcal(snum,s,down),tcal(snum,s,s+1))
  
  for(j in (d+1):(2*d-2)){
    s<-scal(d,j)
    down<-scal(d-1,j)
    anum<-2*d-j
    mnum<-j-d
    cons<-con_add(cons,s,down,anum*(anum-1)/2,1)#q[k,down]<-anum*(anum-1)*k1
    cons<-con_add(cons,s,s-1,mnum,4) #q[k,k-1]<-mnum*k4
    cons<-con_add(cons,s,s+1,anum,3) #q[k,k+1]<-anum*k3
    target.i<-c(target.i,tcal(snum,s,down),tcal(snum,s,s-1),tcal(snum,s,s+1))
  }
  
  s<-scal(d,2*d-1)
  cons<-con_add(cons,s,s-1,d-1,4) #q[k,k-1]<-(d-1)*k4
  cons<-con_add(cons,s,s+1,1,3) #q[k,k+1]<-k3
  target.i<-c(target.i,tcal(snum,s,s-1),tcal(snum,s,s+1))
  
  s<-scal(d,2*d)
  cons<-con_add(cons,s,s-1,d,4) #q[k,k-1]<-d*k4
  target.i<-c(target.i,tcal(snum,s,s-1))
  
  return(list(cons,target.i))
}

get_exp<-function(form,ind){
  return(as.formula(form)[[ind]])
}

## Function, constrain.kt
### This function replace "constrain" function that produced errors in our analysis.

constrain.kt <-function(f,formulae=NULL,free.i,target.i,names=argnames(f)){
  
  final<-names[free.i]
  rels<-lapply(formulae,get_exp,3)
  names(rels)<-as.character(lapply(formulae,get_exp,2))
  target.names<-names[target.i]
  target.i <- target.i[match(names(rels),target.names)]
  pars.out=rep(0,length=length(names))
  names(pars.out)<-names #really needed?
  
  g<- function(pars, ...,pars.only=FALSE){
    pars.out[free.i]<-pars
    e<-structure(as.list(pars),names=final)
    pars.out[target.i]<-unlist(lapply(rels,eval,e))
    if(pars.only)
      pars.out
    else
      f(pars.out, ...)
  }
  class (g) <- c("constrained",class(f))
  attr(g, "argnames")<-final
  attr(g, "formulae")<-formulae #substitution of formulae attr is not precise 
  attr(g, "func") <- f
  g
}

###### Functions for MuSSE EC ######

## Function, get_cons_and_target_musse
### almost same as get_cons_and_target but lambda and mu (speciation rate and extinction rate, respectively) 
### This function is specific to the test to know one given state have different lambda and mu from the other states.
### sp_state is the state number (calculated with scal)

get_cons_and_target_musse_ec<-function(chrmax,sp_state){
  
  snum<-(chrmax+1)*(chrmax+2)/2-1
  other_state<-1:snum
  other_state<-other_state[-c(1,sp_state)]
  ord<-as.integer(log10(snum))+1
  ord2<-as.character(sprintf("%02d",ord))
  
  formula<-paste("lambda%",ord2,"d~lambda%",ord2,"d",sep="")
  cons<-as.list(sprintf(formula,other_state,1))
  formula<-paste("mu%",ord2,"d~mu%",ord2,"d",sep="")
  cons<-c(cons,as.list(sprintf(formula,other_state,1)))
  
  ## The first snum*2 target.i are for lambda and mu
  
  target.i<-1:(snum*2)
  target.i<-target.i[-c(1,sp_state,(snum+1),(snum+sp_state))]
  
  formula<-paste("q%",ord2,"d%",ord2,"d~%d*%s",sep="")
  kdel<-c(qwrite(ord,3,2),qwrite(ord,2,3),qwrite(ord,1,2),qwrite(ord,2,1))
  con_add<-function(cons,is,ts,kcoef,knum){
    return(c(cons,list(sprintf(formula,is,ts,kcoef,kdel[knum]))))
  }
  
  cons<-con_add(cons,3,4,2,3) #q[3,4]<- 2*k3
  cons<-con_add(cons,4,3,1,4) #q[4,3]<- k4
  cons<-con_add(cons,4,5,1,3) #q[4,5]<- k3
  cons<-con_add(cons,4,6,1,2) #q[4,6]<- k2
  cons<-con_add(cons,5,4,2,4) #q[5,4]<- 2*k4
  cons<-con_add(cons,5,7,2,2) #q[5,7]<- 2*k2
  targets<-c(tcal(snum,3,4),tcal(snum,4,3),tcal(snum,4,5),tcal(snum,4,6),tcal(snum,5,4),tcal(snum,5,7))
  for(d in 3:(chrmax-1)){
    
    s<-scal(d,d)
    down<-scal(d-1,d)
    cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1
    cons<-con_add(cons,s,s+1,d,3) #q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
    targets<-c(targets,tcal(snum,s,down),tcal(snum,s,s+1))
    
    for(j in (d+1):(2*d-2)){
      s<-scal(d,j)
      up<-scal(d+1,j)
      down<-scal(d-1,j)
      anum<-2*d-j
      mnum<-j-d
      cons<-con_add(cons,s,down,anum*(anum-1)/2,1)#q[k,down]<-anum*(anum-1)*k1
      cons<-con_add(cons,s,s-1,mnum,4) #q[k,k-1]<-mnum*k4
      cons<-con_add(cons,s,s+1,anum,3) #q[k,k+1]<-anum*k3
      cons<-con_add(cons,s,up,mnum,2)#q[k,up]<-mnum*k2
      targets<-c(targets,tcal(snum,s,down),tcal(snum,s,s-1),tcal(snum,s,s+1),tcal(snum,s,up))
    }
    
    s<-scal(d,2*d-1)
    up<-scal(d+1,2*d-1)
    cons<-con_add(cons,s,s-1,d-1,4) #q[k,k-1]<-(d-1)*k4
    cons<-con_add(cons,s,s+1,1,3) #q[k,k+1]<-k3
    cons<-con_add(cons,s,up,d-1,2)#q[k,up]<-(d-1)*k2
    targets<-c(targets,tcal(snum,s,s-1),tcal(snum,s,s+1),tcal(snum,s,up))
    
    s<-scal(d,2*d)
    up<-scal(d+1,2*d-1)
    cons<-con_add(cons,s,s-1,d,4) #q[k,k-1]<-d*k4
    cons<-con_add(cons,s,up,d,2)#q[k,up]<-d*k2
    targets<-c(targets,tcal(snum,s,s-1),tcal(snum,s,up))
    
  }
  
  d<-chrmax
  
  s<-scal(d,d)
  down<-scal(d-1,d)
  cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1
  cons<-con_add(cons,s,s+1,d,3) #q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
  targets<-c(targets,tcal(snum,s,down),tcal(snum,s,s+1))
  
  for(j in (d+1):(2*d-2)){
    s<-scal(d,j)
    down<-scal(d-1,j)
    anum<-2*d-j
    mnum<-j-d
    cons<-con_add(cons,s,down,anum*(anum-1)/2,1)#q[k,down]<-anum*(anum-1)*k1
    cons<-con_add(cons,s,s-1,mnum,4) #q[k,k-1]<-mnum*k4
    cons<-con_add(cons,s,s+1,anum,3) #q[k,k+1]<-anum*k3
    targets<-c(targets,tcal(snum,s,down),tcal(snum,s,s-1),tcal(snum,s,s+1))
    
  }
  
  s<-scal(d,2*d-1)
  cons<-con_add(cons,s,s-1,d-1,4) #q[k,k-1]<-(d-1)*k4
  cons<-con_add(cons,s,s+1,1,3) #q[k,k+1]<-k3
  targets<-c(targets,tcal(snum,s,s-1),tcal(snum,s,s+1))
  
  s<-scal(d,2*d)
  cons<-con_add(cons,s,s-1,d,4) #q[k,k-1]<-d*k4
  targets<-c(targets,tcal(snum,s,s-1))
  
  ## The target.i from snum*2+1 are for the transition rate.
  
  target.i<-c(target.i,(targets+((snum)*2)))
  
  return(list(cons,target.i))
}

## Function, get_cons_and_target_musse_null
### almost same as get_cons_and_target but no difference in lambda and mu between states

get_cons_and_target_musse_null_ec<-function(chrmax){
  
  snum<-(chrmax+1)*(chrmax+2)/2-1
  all_state<-2:snum
  ord<-as.integer(log10(snum))+1
  ord2<-as.character(sprintf("%02d",ord))
  
  formula<-paste("lambda%",ord2,"d~lambda%",ord2,"d",sep="")
  cons<-as.list(sprintf(formula,all_state,1))
  formula<-paste("mu%",ord2,"d~mu%",ord2,"d",sep="")
  cons<-c(cons,as.list(sprintf(formula,all_state,1)))
  target.i<-1:(snum*2)
  target.i<-target.i[-c(1,(snum+1))]
  
  formula<-paste("q%",ord2,"d%",ord2,"d~%d*%s",sep="")
  kdel<-c(qwrite(ord,3,2),qwrite(ord,2,3),qwrite(ord,1,2),qwrite(ord,2,1))
  con_add<-function(cons,is,ts,kcoef,knum){
    return(c(cons,list(sprintf(formula,is,ts,kcoef,kdel[knum]))))
  }
  
  cons<-con_add(cons,3,4,2,3) #q[3,4]<- 2*k3
  cons<-con_add(cons,4,3,1,4) #q[4,3]<- k4
  cons<-con_add(cons,4,5,1,3) #q[4,5]<- k3
  cons<-con_add(cons,4,6,1,2) #q[4,6]<- k2
  cons<-con_add(cons,5,4,2,4) #q[5,4]<- 2*k4
  cons<-con_add(cons,5,7,2,2) #q[5,7]<- 2*k2
  targets<-c(tcal(snum,3,4),tcal(snum,4,3),tcal(snum,4,5),tcal(snum,4,6),tcal(snum,5,4),tcal(snum,5,7))
  for(d in 3:(chrmax-1)){
    
    s<-scal(d,d)
    down<-scal(d-1,d)
    cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1
    cons<-con_add(cons,s,s+1,d,3) #q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
    targets<-c(targets,tcal(snum,s,down),tcal(snum,s,s+1))
    
    for(j in (d+1):(2*d-2)){
      s<-scal(d,j)
      up<-scal(d+1,j)
      down<-scal(d-1,j)
      anum<-2*d-j
      mnum<-j-d
      cons<-con_add(cons,s,down,anum*(anum-1)/2,1)#q[k,down]<-anum*(anum-1)*k1
      cons<-con_add(cons,s,s-1,mnum,4) #q[k,k-1]<-mnum*k4
      cons<-con_add(cons,s,s+1,anum,3) #q[k,k+1]<-anum*k3
      cons<-con_add(cons,s,up,mnum,2)#q[k,up]<-mnum*k2
      targets<-c(targets,tcal(snum,s,down),tcal(snum,s,s-1),tcal(snum,s,s+1),tcal(snum,s,up))
    }
    
    s<-scal(d,2*d-1)
    up<-scal(d+1,2*d-1)
    cons<-con_add(cons,s,s-1,d-1,4) #q[k,k-1]<-(d-1)*k4
    cons<-con_add(cons,s,s+1,1,3) #q[k,k+1]<-k3
    cons<-con_add(cons,s,up,d-1,2)#q[k,up]<-(d-1)*k2
    targets<-c(targets,tcal(snum,s,s-1),tcal(snum,s,s+1),tcal(snum,s,up))
    
    s<-scal(d,2*d)
    up<-scal(d+1,2*d-1)
    cons<-con_add(cons,s,s-1,d,4) #q[k,k-1]<-d*k4
    cons<-con_add(cons,s,up,d,2)#q[k,up]<-d*k2
    targets<-c(targets,tcal(snum,s,s-1),tcal(snum,s,up))
    
  }
  
  d<-chrmax
  
  s<-scal(d,d)
  down<-scal(d-1,d)
  cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1
  cons<-con_add(cons,s,s+1,d,3) #q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
  targets<-c(targets,tcal(snum,s,down),tcal(snum,s,s+1))
  
  for(j in (d+1):(2*d-2)){
    s<-scal(d,j)
    down<-scal(d-1,j)
    anum<-2*d-j
    mnum<-j-d
    cons<-con_add(cons,s,down,anum*(anum-1)/2,1)#q[k,down]<-anum*(anum-1)*k1
    cons<-con_add(cons,s,s-1,mnum,4) #q[k,k-1]<-mnum*k4
    cons<-con_add(cons,s,s+1,anum,3) #q[k,k+1]<-anum*k3
    targets<-c(targets,tcal(snum,s,down),tcal(snum,s,s-1),tcal(snum,s,s+1))
  }
  
  s<-scal(d,2*d-1)
  cons<-con_add(cons,s,s-1,d-1,4) #q[k,k-1]<-(d-1)*k4
  cons<-con_add(cons,s,s+1,1,3) #q[k,k+1]<-k3
  targets<-c(targets,tcal(snum,s,s-1),tcal(snum,s,s+1))
  
  s<-scal(d,2*d)
  cons<-con_add(cons,s,s-1,d,4) #q[k,k-1]<-d*k4
  targets<-c(targets,tcal(snum,s,s-1))
  target.i<-c(target.i,(targets+((snum)*2)))
  
  return(list(cons,target.i))
}

###### Preparation of the tree ######

## Load tree and karyotype data
phy.fish<-read.tree("Rabosky_et_al_timetree.tre")

## Rounding the tree
tipnum<-Ntip(phy.fish)
edgeindex<-which(phy.fish$edge[,2]<=tipnum)
excesslen<-node.depth.edgelength(phy.fish)[1:tipnum]-345.0000
phy.fish$edge.length[edgeindex]<-phy.fish$edge.length[edgeindex]-excesslen
names.fish<-phy.fish$tip.label
names.fish<-gsub("_"," ", names.fish)

## Binding names and karyotype states to the tip

load("neot_state_YK2021.Robj")
keep<-which(names.fish %in% names(neot.state))
phy.o.neot<-keep.tip(phy.fish,keep)
phy.o.neot$tip.label<-names(neot.state)
phy.o.neot$tip.state<-neot.state

## Remove the tip with out of range in PCM

snum<-(chrmax+1)*(chrmax+2)/2-1

phy.o.neot<-keep.tip(phy.o.neot,which(neot.state<=snum))
phy.o.neot$tip.state<-phy.o.neot$tip.state[neot.state<=snum]

###### Sampling.f preparation ######
## To compare M2, sampling.f is set as well as M2.

load("val_chrinfo_YK2021.Robj")

neotlist<-c("AULOPIFORMES","MYCTOPHIFORMES","GADIFORMES","OPHIDIIFORMES",
            "MUGILIFORMES","ATHERINIFORMES","BELONIFORMES","CYPRINODONTIFORMES","STEPHANOBERYCIFORMES",
            "BERYCIFORMES","GASTEROSTEIFORMES","BATRACHOIDIFORMES","SYNBRANCHIFORMES","SCORPAENIFORMES",
            "PERCIFORMES","PLEURONECTIFORMES","LOPHIIFORMES","TETRAODONTIFORMES","ZEIFORMES","OPHIDIIFORMES")

## Determine sampling.f

val_chrinfo<-val_chrinfo[val_chrinfo$Order %in% neotlist,]
val_chrinfo$State<-factor(val_chrinfo$State,levels=1:snum)
val_chrinfo<-val_chrinfo[!(is.na(val_chrinfo$State)),]
tipstate<-factor(phy.o.neot$tip.state,levels=1:snum)
sf300<-sum(tipstate==300)/sum(val_chrinfo$State==300)
sfn300<-sum(tipstate!=300)/sum(val_chrinfo$State!=300)
sampling.f<-c(rep(sfn300,299),sf300,rep(sfn300,snum-300))


###### MLE of MuSSE M0 with sampling.f ######

mussefunc<-make.musse(phy.o.neot,phy.o.neot$tip.state,snum,sampling.f=sampling.f,strict=F)
free.i<-c(tcal(snum,1,2),tcal(snum,2,1),tcal(snum,2,3),tcal(snum,3,2))
free.i<-c(1,snum+1,((snum*2)+free.i))
cons_and_tar<-get_cons_and_target_musse_null_ec(chrmax)
mussefunc<-constrain.kt(mussefunc,cons_and_tar[[1]],free.i,cons_and_tar[[2]])
int.p<-starting.point.musse(phy.o.neot,k=snum)
p<-int.p[free.i] #prior
names(p)<-c()

#MLE
print(paste("Start MLE of neot MuSSE M0 prime CM",chrmax, " to compare test SE..."))
timepoint<-c(proc.time()[[3]])
fit.musse.neot<-find.mle(mussefunc,p)
timepoint<-c(timepoint,proc.time()[[3]])

## Save fit.musse.neot
filename<-paste("fit_musse_neot_M0_prime_CM",chrmax,"_samplingf_XXXXXX.Robj",sep="")
save(fit.musse.neot,file=filename)

## Result of the 6 para Estimation

coef_res<-c(coef(fit.musse.neot)[c(6,5,3,4,1,2)]) #2020 definition
names(coef_res)<-c("k1","k2","k3","k4","lambda","mu")
coef_res #coef
fit.musse.neot[2] #lnLik
paste("Chromosome limit,",chrmax, ": Process time,", round(timepoint[2]-timepoint[1],3))