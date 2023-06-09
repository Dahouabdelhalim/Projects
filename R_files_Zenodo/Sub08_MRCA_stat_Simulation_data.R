###### Library ######

library(diversitree)
library(Matrix)
library(stats4)


###### Font preparation ######
subset(fonttable(), FamilyName == "Arial") ##Check not more than one "Arial" in your system
loadfonts(device = "postscript")
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

###### Additional function ######

#For 95% range
chrmax<-35
snum<-(chrmax+1)*(chrmax+2)/2-1
stch.mat<-matrix(rep(0),nrow=snum,ncol=(chrmax))
star.mat<-matrix(rep(0),nrow=snum,ncol=(chrmax*2))
stch.vec<-c()
star.vec<-c()
for(d in 1:chrmax){
  for(j in d:(2*d)){
    stch.vec<-c(stch.vec,d)
    star.vec<-c(star.vec,j)
  }
}

mean95range<-function(st.node.matrix,convert.matrix){
  reschmat<-t(st.node.matrix) %*% convert.matrix
  res<-c()
  for(i in 1:nrow(reschmat)){
    probs<-reschmat[i,]
    cumprob<-0
    range95<-c()
    chrindex<-1:ncol(convert.matrix)
    while(cumprob<0.95){
      index<-which(probs==max(probs))
      range95<-c(range95,chrindex[index])
      cumprob<-cumprob+probs[index]
      chrindex<-chrindex[-index]
      probs<-probs[-index]
    }
    lower<-min(range95)
    upper<-max(range95)
    res<-rbind(res,c(lower,upper))
  }
  rangeoutput<-apply(res,2,mean)
  return(rangeoutput)
}

# for coloring the ASR plot
coloring_prob<-function(ps){
  ordered_asr<-sort(ps,decreasing=T)
  range50<-sum(cumsum(ordered_asr)<0.50)+1
  range75<-sum(cumsum(ordered_asr)<0.75)+1
  range90<-sum(cumsum(ordered_asr)<0.90)+1
  pc<-rep("100%",length(ps))
  pc[ps>=ordered_asr[range90]]<-rep("90%")
  pc[ps>=ordered_asr[range75]]<-rep("75%")
  pc[ps>=ordered_asr[range50]]<-rep("50%")
  return(pc)
}

###### Preparation of the tree ######

load("phy_o_neot_YK2021.Robj")

chrmax<-35
snum<-(chrmax+1)*(chrmax+2)/2-1
sp_state<-300 #state number of (24,24)

tipnum<-Ntip(phy.o.neot)
edgeindex<-which(phy.o.neot$edge[,2]<=tipnum)
meanedge<-mean(node.depth.edgelength(phy.o.neot)[1:tipnum])
excesslen<-node.depth.edgelength(phy.o.neot)[1:tipnum]-meanedge
phy.o.neot$edge.length[edgeindex]<-phy.o.neot$edge.length[edgeindex]-excesslen

###### ASR marginal MuSSE Eurypterygii ######

phy.sim<-phy.o.neot
res.data<-c()
for(trial in 1:100){
  print(paste("Eurypterygii trial",trial))
  filename<-paste("sim_traits_neot_trial",trial,".Robj",sep="")
  load(filename)
  phy.sim$tip.state<-as.vector(sim_traits)
  names(phy.sim$tip.state)<-names(sim_traits)
  mussefunc<-make.musse(phy.sim,phy.sim$tip.state,snum,strict=F)
  free.i<-c(tcal(snum,1,2),tcal(snum,2,1),tcal(snum,2,3),tcal(snum,3,2))
  free.i<-c(1,snum+1,((snum*2)+free.i))
  cons_and_tar<-get_cons_and_target_musse_null_ec(chrmax)
  mussefunc<-constrain.kt(mussefunc,cons_and_tar[[1]],free.i,cons_and_tar[[2]])
  
  load(file="fit_musse_neot_M0_YK2021.Robj")
  asr.mrca<-asr.marginal(mussefunc,coef(fit.musse.neot),nodes=1)
  res.data<-cbind(res.data,asr.mrca)
}

## Save result data
write.table(res.data,file="MRCA_ASR_sim_neot815_XXXXXX.txt")

## 95% range
mean95range(res.data,stch.mat)
mean95range(res.data,star.mat)

neot815asr<-res.data

###### ASR marginal MuSSE Cyprinodontiformes ######
phy.sim<-keep.tip(phy.o.neot,436:515)
res.data<-c()
for(trial in 1:100){
  print(paste("Cyprinodontiformes trial",trial))
  filename<-paste("sim_traits_Cyprinodon_trial",trial,".Robj",sep="")
  load(filename)
  phy.sim$tip.state<-as.vector(sim_traits)
  names(phy.sim$tip.state)<-names(sim_traits)
  mussefunc<-make.musse(phy.sim,phy.sim$tip.state,snum,strict=F)
  free.i<-c(tcal(snum,1,2),tcal(snum,2,1),tcal(snum,2,3),tcal(snum,3,2))
  free.i<-c(1,snum+1,((snum*2)+free.i))
  cons_and_tar<-get_cons_and_target_musse_null_ec(chrmax)
  mussefunc<-constrain.kt(mussefunc,cons_and_tar[[1]],free.i,cons_and_tar[[2]])
  
  load(file="fit_musse_neot_M0_YK2021.Robj")
  asr.mrca<-asr.marginal(mussefunc,coef(fit.musse.neot),nodes=1)
  res.data<-cbind(res.data,asr.mrca)
}

## Save result data
write.table(res.data,file="MRCA_ASR_sim_Cyprinodon_XXXXXX.txt")

## 95% range
mean95range(res.data,stch.mat)
mean95range(res.data,star.mat)

cyprinoasr<-res.data

###### ASR marginal MuSSE Goodeidae ######
phy.sim<-keep.tip(phy.o.neot,452:480)
res.data<-c()
for(trial in 1:100){
  print(paste("Goodeidae trial",trial))
  filename<-paste("sim_traits_Goodeidae_trial",trial,".Robj",sep="")
  load(filename)
  phy.sim$tip.state<-as.vector(sim_traits)
  names(phy.sim$tip.state)<-names(sim_traits)
  mussefunc<-make.musse(phy.sim,phy.sim$tip.state,snum,strict=F)
  free.i<-c(tcal(snum,1,2),tcal(snum,2,1),tcal(snum,2,3),tcal(snum,3,2))
  free.i<-c(1,snum+1,((snum*2)+free.i))
  cons_and_tar<-get_cons_and_target_musse_null_ec(chrmax)
  mussefunc<-constrain.kt(mussefunc,cons_and_tar[[1]],free.i,cons_and_tar[[2]])
  
  load(file="fit_musse_neot_M0_YK2021.Robj")
  asr.mrca<-asr.marginal(mussefunc,coef(fit.musse.neot),nodes=1)
  res.data<-cbind(res.data,asr.mrca)
}

## Save result data
write.table(res.data,file="MRCA_ASR_sim_Goodeidae_XXXXXX.txt")

## 95% range
mean95range(res.data,stch.mat)
mean95range(res.data,star.mat)

godeiasr<-res.data


###### ASR karyograph plot of 3 cases each (S3 Fig) ######

# plot data for Eurypterygii
pdata1<-data.frame(cbind(rep(c("Case 1","Case 2","Case 3"),each=665),
                         matrix(rep(as.vector(t(chr_arm_vec(1:665))),3),byrow=T,ncol=2),
                         c(neot815asr[,1],neot815asr[,2],neot815asr[,3])))
colnames(pdata1)<-c("Case","py","px","ps")
pdata1$Case<-factor(pdata1$Case)
pdata1$px<-as.integer(pdata1$px)
pdata1$py<-as.integer(pdata1$py)
pdata1$ps<-as.numeric(pdata1$ps)
pdata1$pc<-c(coloring_prob(neot815asr[,1]),coloring_prob(neot815asr[,2]),coloring_prob(neot815asr[,3]))
pdata1$pc<-factor(pdata1$pc,level=c("50%","75%","90%","100%"))
pdata1$Taxon<-rep("Eurypterygii")

# plot data for Cyprinodontiformes
pdata2<-data.frame(cbind(rep(c("Case 1","Case 2","Case 3"),each=665),
                         matrix(rep(as.vector(t(chr_arm_vec(1:665))),3),byrow=T,ncol=2),
                         c(cyprinoasr[,1],cyprinoasr[,2],cyprinoasr[,3])))
colnames(pdata2)<-c("Case","py","px","ps")
pdata2$Case<-factor(pdata2$Case)
pdata2$px<-as.integer(pdata2$px)
pdata2$py<-as.integer(pdata2$py)
pdata2$ps<-as.numeric(pdata2$ps)
pdata2$pc<-c(coloring_prob(cyprinoasr[,1]),coloring_prob(cyprinoasr[,2]),coloring_prob(cyprinoasr[,3]))
pdata2$pc<-factor(pdata2$pc,level=c("50%","75%","90%","100%"))
pdata2$Taxon<-rep("Cyprinodontiformes")

# plot data for Goodeidae
pdata3<-data.frame(cbind(rep(c("Case 1","Case 2","Case 3"),each=665),
                         matrix(rep(as.vector(t(chr_arm_vec(1:665))),3),byrow=T,ncol=2),
                         c(godeiasr[,1],godeiasr[,2],godeiasr[,3])))
colnames(pdata3)<-c("Case","py","px","ps")
pdata3$Case<-factor(pdata3$Case)
pdata3$px<-as.integer(pdata3$px)
pdata3$py<-as.integer(pdata3$py)
pdata3$ps<-as.numeric(pdata3$ps)
pdata3$pc<-c(coloring_prob(godeiasr[,1]),coloring_prob(godeiasr[,2]),coloring_prob(meanprob))
pdata3$pc<-factor(pdata3$pc,level=c("50%","75%","90%","100%"))
pdata3$Taxon<-rep("Goodeidae")

# Plot karyograph
pdata<-rbind(pdata1,pdata2,pdata3)
pdata$Taxon<-factor(pdata$Taxon,level=c("Eurypterygii","Cyprinodontiformes","Goodeidae"))
library(ggplot2)
postscript("Fig_simulation_3taxon_300_ASR_XXXXXX.eps", height = 4, width = 7.4,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
xtitle<-expression(paste("Arm number, ",italic("x")))
ytitle<-expression(paste("Chromosome number, ",italic("y")))
ggplot(pdata)+
  geom_point(aes(x=px,y=py,size=ps,col=pc))+
  scale_color_manual(values=c("firebrick3","darkorange2","steelblue1","grey60"))+
  labs(x=xtitle,y=ytitle,size="Probability",color="Top") +
  facet_grid(Taxon~Case)+
  scale_size(range=c(-0.3,0.9),limits=c(0,0.90))+
  guides(size = guide_legend(order = 1),
         color = guide_legend(order = 2))+
  theme(text=element_text(size=8),
        panel.border=element_rect(colour="black",size=0.6,fill=NA),
        strip.background = element_rect(color="black",size=0.5),
        legend.position="right",
        legend.spacing=unit(0,"in"),
        legend.key.size=unit(0.15,"in"),
        legend.box.spacing=unit(0,"in"),
        strip.text=element_text(size=7.5))
dev.off()

###### ASR karyograph plot of mean (S4 Fig) ######

## Plot data
# Eurypterygii
pdata1<-data.frame(cbind(rep("Mean (N = 100)",665),
                         matrix(as.vector(t(chr_arm_vec(1:665))),byrow=T,ncol=2),
                         meanprob))
colnames(pdata1)<-c("Case","py","px","ps")
pdata1$Case<-factor(pdata1$Case)
pdata1$px<-as.integer(pdata1$px)
pdata1$py<-as.integer(pdata1$py)
pdata1$ps<-as.numeric(pdata1$ps)
pdata1$pc<-coloring_prob(meanprob)
pdata1$pc<-factor(pdata1$pc,level=c("50%","75%","90%","100%"))
pdata1$Taxon<-rep("Eurypterygii")

# Cyprinodontiformes
pdata2<-data.frame(cbind(rep("Mean (N = 100)",665),
                         matrix(as.vector(t(chr_arm_vec(1:665))),byrow=T,ncol=2),
                         meanprob))
colnames(pdata2)<-c("Case","py","px","ps")
pdata2$Case<-factor(pdata2$Case)
pdata2$px<-as.integer(pdata2$px)
pdata2$py<-as.integer(pdata2$py)
pdata2$ps<-as.numeric(pdata2$ps)
pdata2$pc<-coloring_prob(meanprob)
pdata2$pc<-factor(pdata2$pc,level=c("50%","75%","90%","100%"))
pdata2$Taxon<-rep("Cyprinodontiformes")

#Goodeidae
pdata3<-data.frame(cbind(rep("Mean (N = 100)",665),
                         matrix(as.vector(t(chr_arm_vec(1:665))),byrow=T,ncol=2),
                         meanprob))
colnames(pdata3)<-c("Case","py","px","ps")
pdata3$Case<-factor(pdata3$Case)
pdata3$px<-as.integer(pdata3$px)
pdata3$py<-as.integer(pdata3$py)
pdata3$ps<-as.numeric(pdata3$ps)
pdata3$pc<-coloring_prob(meanprob)
pdata3$pc<-factor(pdata3$pc,level=c("50%","75%","90%","100%"))
pdata3$Taxon<-rep("Goodeidae")
pdata<-rbind(pdata1,pdata2,pdata3)
pdata$Taxon<-factor(pdata$Taxon,level=c("Eurypterygii","Cyprinodontiformes","Goodeidae"))

# Plot karyograph

levels(pdata$Case)=c("Mean (N = 100)"=expression(paste("Mean (",italic('N')," = 100)")))

library(ggplot2)
postscript("Fig_simulation_3taxon_mean_ASR_XXXXXX.eps", height = 6, width = 5,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)

xtitle<-expression(paste("Arm number, ",italic("x")))
ytitle<-expression(paste("Chromosome number, ",italic("y")))
ggplot(pdata)+
  geom_point(aes(x=px,y=py,size=ps,col=pc))+
  scale_color_manual(values=c("firebrick3","darkorange2","steelblue1","grey60"))+
  labs(x=xtitle,y=ytitle,size="Probability",color="Top") +
  facet_grid(Taxon~Case,labeller=label_parsed)+
  scale_size(range=c(0,2),limits=c(0,0.90))+
  guides(size = guide_legend(order = 1),
         color = guide_legend(order = 2))+
  theme(text=element_text(size=10),
        panel.border=element_rect(colour="black",size=0.6,fill=NA),
        legend.position="right",
        legend.spacing=unit(0,"in"),
        legend.key.size=unit(0.15,"in"),
        legend.box.spacing=unit(0,"in"),
        strip.text.y=element_text(size=6),
        strip.background = element_rect(colour="black")
  )
dev.off()