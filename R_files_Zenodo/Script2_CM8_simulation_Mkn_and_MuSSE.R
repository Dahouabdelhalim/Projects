###### Parameter ######
## Give id number of starting trial in args[1]
## Give id number of end trial in args[2]
chrmax<-8
args = commandArgs(trailingOnly=TRUE)
strial<-as.numeric(args[1])
etrial<-as.numeric(args[2])

###### Library ######

library(diversitree)

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



###### Function for Mkn and MuSSE ######

## function, tcal
### Aquisition of number of transition rate 
### snum, total number of states
### is, initial state
### ts, transited state
### The no. of transition is
### snum x snum - snum (diagonals)
### right entries of diagonals have shifted index

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


## For Mk-n, function get_cons_and_target_ec
### Mk-n, 4 parameter
### Constraints of the model are expressed as formulae by this function.
### This function also output target.i indicating index of target parameters.
### q003002,q002003,q001002 and q002001 are free parameters here,
### which corresponding to k1, k2, k3 and k4, respectively.
### The constraints are made with these 4 parameters.
### chrmax, upper limit of chromosome number (>=4)

get_cons_and_target_ec<-function(chrmax){
  
  #state number
  snum<-(chrmax+1)*(chrmax+2)/2-1
  
  #In diversitree, the transition is written as q001003 to show transition from state1 to state3 for example.
  #The digit of the state number is depending on the total number of the state
  #ord and ord2 shows this digit number.
  ord<-as.integer(log10(snum))+1
  ord2<-as.character(sprintf("%02d",ord))
  
  #formula, that represents "q%03d%03d~%d*%s" for example.
  #this is used in sprintf
  formula<-paste("q%",ord2,"d%",ord2,"d~%d*%s",sep="")
  
  #kdel, that gives the transition rate corresponding to 4 parameters
  kdel<-c(qwrite(ord,3,2),qwrite(ord,2,3),qwrite(ord,1,2),qwrite(ord,2,1))
  
  #con_add, express and add constraint formulae to the given cons vector
  #is, initial state
  #ts, transited state
  #kcoef, coefficient multiplied to kx
  #knum, the index of parameters in the sense, kx
  con_add<-function(cons,is,ts,kcoef,knum){
    return(c(cons,list(sprintf(formula,is,ts,kcoef,kdel[knum]))))
  }
  
  #the states of chromosome number(y) = 1 or 2, which have limited directions
  cons<-list()
  cons<-con_add(cons,3,4,2,3) #q[3,4]<- 2*k3
  cons<-con_add(cons,4,3,1,4) #q[4,3]<- k4
  cons<-con_add(cons,4,5,1,3) #q[4,5]<- k3
  cons<-con_add(cons,4,6,1,2) #q[4,6]<- k2
  cons<-con_add(cons,5,4,2,4) #q[5,4]<- 2*k4
  cons<-con_add(cons,5,7,2,2) #q[5,7]<- 2*k2
  target.i<-c(tcal(snum,3,4),tcal(snum,4,3),tcal(snum,4,5),tcal(snum,4,6),tcal(snum,5,4),tcal(snum,5,7))
  
  #the states of chromosome number(y) = 3 to (limit-1)
  for(d in 3:(chrmax-1)){
    
    # the states of chromosome number(y) = arm number(x), which have two directions
    s<-scal(d,d)
    down<-scal(d-1,d)
    cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1
    cons<-con_add(cons,s,s+1,d,3) #q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
    target.i<-c(target.i,tcal(snum,s,down),tcal(snum,s,s+1))
    
    # the states of y < x < 2*y, which have four directions
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
    
    # the state of x = 2*y-1, which have three directions
    s<-scal(d,2*d-1)
    up<-scal(d+1,2*d-1)
    cons<-con_add(cons,s,s-1,d-1,4) #q[k,k-1]<-(d-1)*k4
    cons<-con_add(cons,s,s+1,1,3) #q[k,k+1]<-k3
    cons<-con_add(cons,s,up,d-1,2)#q[k,up]<-(d-1)*k2
    target.i<-c(target.i,tcal(snum,s,s-1),tcal(snum,s,s+1),tcal(snum,s,up))
    
    # the state of x = 2*y, which have two directions
    s<-scal(d,2*d)
    up<-scal(d+1,2*d-1)
    cons<-con_add(cons,s,s-1,d,4) #q[k,k-1]<-d*k4
    cons<-con_add(cons,s,up,d,2)#q[k,up]<-d*k2
    target.i<-c(target.i,tcal(snum,s,s-1),tcal(snum,s,up))
    
  }
  
  # the states of y = limit
  d<-chrmax
  
  # the state of y = x = limit, which has two directions
  s<-scal(d,d)
  down<-scal(d-1,d)
  cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1
  cons<-con_add(cons,s,s+1,d,3) #q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
  target.i<-c(target.i,tcal(snum,s,down),tcal(snum,s,s+1))
  
  # the states of y = limit < x < 2y, which have three directions
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
  
  # the state of x = 2y-1 = 2*limit-1, which has two directions
  s<-scal(d,2*d-1)
  cons<-con_add(cons,s,s-1,d-1,4) #q[k,k-1]<-(d-1)*k4
  cons<-con_add(cons,s,s+1,1,3) #q[k,k+1]<-k3
  target.i<-c(target.i,tcal(snum,s,s-1),tcal(snum,s,s+1))
  
  # the state of x = 2y = 2*limit, which has one direction
  s<-scal(d,2*d)
  cons<-con_add(cons,s,s-1,d,4) #q[k,k-1]<-d*k4
  target.i<-c(target.i,tcal(snum,s,s-1))
  
  return(list(cons,target.i))
}


## For MuSSE, Function, get_cons_and_target_musse_null_ec
### M0 model, 6 parameters
### Constraints of the model are expressed as formulae by this function.
### This function also output target.i indicating index of target parameters.
### lambda and mu express speciation and extinction rates.
### all labda or all mu are assumed as constant as lambda001 or mu001, respectively.
### q003002,q002003,q001002 and q002001 are free parameters,
### which corresponding to k1, k2, k3 and k4, respectively.
### chrmax, upper limit of chromosome number

get_cons_and_target_musse_null_ec<-function(chrmax){
  
  # state number
  snum<-(chrmax+1)*(chrmax+2)/2-1
  
  #In diversitree, the transition is written as q001003 to show transition from state1 to state3 for example.
  #The speciation and extinction rate are also expressed like lambda001 and mu001 for example.
  #The digit of the state number is depending on the total number of the state
  #ord and ord2 shows this digit number.
  ord<-as.integer(log10(snum))+1
  ord2<-as.character(sprintf("%02d",ord))
  
  #state 1 is given as the representative state for speciation and extinction rate of all states
  all_state<-2:snum
  #This formula represents "lambda%03~lambdad%03" for example.
  formula<-paste("lambda%",ord2,"d~lambda%",ord2,"d",sep="")
  #all constraint formulae for lambda are put in the vector
  cons<-as.list(sprintf(formula,all_state,1))
  
  #This formula represents "mu%03~mu%03" for example.
  formula<-paste("mu%",ord2,"d~mu%",ord2,"d",sep="")
  #all constraint formulae for mu are put in the vector
  cons<-c(cons,as.list(sprintf(formula,all_state,1)))
  
  #target.i for lambda and mu
  target.i<-1:(snum*2)
  target.i<-target.i[-c(1,(snum+1))]
  
  #formula, that represents "q%03d%03d~%d*%s" for example.
  #this is used in sprintf
  formula<-paste("q%",ord2,"d%",ord2,"d~%d*%s",sep="")
  
  #kdel, that gives the transition rate corresponding to 4 parameters
  kdel<-c(qwrite(ord,3,2),qwrite(ord,2,3),qwrite(ord,1,2),qwrite(ord,2,1))
  
  #con_add, express and add constraint formulae to the given cons vector
  #is, initial state
  #ts, transited state
  #kcoef, coefficient multiplied to kx
  #knum, the index of parameters in the sense, kx
  con_add<-function(cons,is,ts,kcoef,knum){
    return(c(cons,list(sprintf(formula,is,ts,kcoef,kdel[knum]))))
  }
  
  #states of y<= 2
  cons<-con_add(cons,3,4,2,3) #q[3,4]<- 2*k3
  cons<-con_add(cons,4,3,1,4) #q[4,3]<- k4
  cons<-con_add(cons,4,5,1,3) #q[4,5]<- k3
  cons<-con_add(cons,4,6,1,2) #q[4,6]<- k2
  cons<-con_add(cons,5,4,2,4) #q[5,4]<- 2*k4
  cons<-con_add(cons,5,7,2,2) #q[5,7]<- 2*k2
  
  #targets has the same numbering as target.i in Mk-n. The number is shifted later.
  targets<-c(tcal(snum,3,4),tcal(snum,4,3),tcal(snum,4,5),tcal(snum,4,6),tcal(snum,5,4),tcal(snum,5,7))
  
  #states of 2<y<limit
  for(d in 3:(chrmax-1)){
    
    #states of y = x, which have two directions
    s<-scal(d,d)
    down<-scal(d-1,d)
    cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1
    cons<-con_add(cons,s,s+1,d,3) #q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
    targets<-c(targets,tcal(snum,s,down),tcal(snum,s,s+1))
    
    #states of y<x<2y, which have four direcction
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
    
    #states of x=2y-1, which have three directions
    s<-scal(d,2*d-1)
    up<-scal(d+1,2*d-1)
    cons<-con_add(cons,s,s-1,d-1,4) #q[k,k-1]<-(d-1)*k4
    cons<-con_add(cons,s,s+1,1,3) #q[k,k+1]<-k3
    cons<-con_add(cons,s,up,d-1,2)#q[k,up]<-(d-1)*k2
    targets<-c(targets,tcal(snum,s,s-1),tcal(snum,s,s+1),tcal(snum,s,up))
    
    #states of x=2y, which have two directions
    s<-scal(d,2*d)
    up<-scal(d+1,2*d-1)
    cons<-con_add(cons,s,s-1,d,4) #q[k,k-1]<-d*k4
    cons<-con_add(cons,s,up,d,2)#q[k,up]<-d*k2
    targets<-c(targets,tcal(snum,s,s-1),tcal(snum,s,up))
    
  }
  
  #states of y=limit
  d<-chrmax
  
  #x=y, two directions
  s<-scal(d,d)
  down<-scal(d-1,d)
  cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1
  cons<-con_add(cons,s,s+1,d,3) #q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
  targets<-c(targets,tcal(snum,s,down),tcal(snum,s,s+1))
  
  #y<x<2y, three directions
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
  
  #x=2y-1, two directions
  s<-scal(d,2*d-1)
  cons<-con_add(cons,s,s-1,d-1,4) #q[k,k-1]<-(d-1)*k4
  cons<-con_add(cons,s,s+1,1,3) #q[k,k+1]<-k3
  targets<-c(targets,tcal(snum,s,s-1),tcal(snum,s,s+1))
  
  #x=2y, one direction
  s<-scal(d,2*d)
  cons<-con_add(cons,s,s-1,d,4) #q[k,k-1]<-d*k4
  targets<-c(targets,tcal(snum,s,s-1))
  
  #target.i is combined with targets (those number are shifted)
  target.i<-c(target.i,(targets+((snum)*2)))
  
  return(list(cons,target.i))
}


## For MuSSE, Function, get_cons_and_target_musse_null_Ki1_ec
### M1 model, 5 parameters
### For null hypothesis with k3 = k4.
### One constraint is added q002001 ~ q001002 to get_cons_and_target_musse_null_ec.
### chrmax, upper limit of chromosome number

get_cons_and_target_musse_null_Ki1_ec<-function(chrmax){
  
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
  
  #for 5par k4 is set as the same transition as k3
  kdel<-c(qwrite(ord,3,2),qwrite(ord,2,3),qwrite(ord,1,2),qwrite(ord,1,2)) #changed for 5 par
  con_add<-function(cons,is,ts,kcoef,knum){
    return(c(cons,list(sprintf(formula,is,ts,kcoef,kdel[knum]))))
  }
  
  cons<-con_add(cons,2,1,1,3) #q[2,1]<- q[1,2]=k3 newly added for 5 par
  cons<-con_add(cons,3,4,2,3) #q[3,4]<- 2*k3
  cons<-con_add(cons,4,3,1,4) #q[4,3]<- k4
  cons<-con_add(cons,4,5,1,3) #q[4,5]<- k3
  cons<-con_add(cons,4,6,1,2) #q[4,6]<- k2
  cons<-con_add(cons,5,4,2,4) #q[5,4]<- 2*k4
  cons<-con_add(cons,5,7,2,2) #q[5,7]<- 2*k2
  targets<-c(tcal(snum,2,1),tcal(snum,3,4),tcal(snum,4,3),tcal(snum,4,5),tcal(snum,4,6),tcal(snum,5,4),tcal(snum,5,7)) #changed for 5 par
  for(d in 3:(chrmax-1)){
    
    s<-scal(d,d)
    down<-scal(d-1,d)
    cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)/2*k1
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


## For MuSSE, Function, get_cons_and_target_musse_testSE_ec
### M2 model, 8 parameters
### Almost same as get_cons_and_target_musse_null_ec.
### But two speciation and two extinction rates were given.
### sp_state specifies states have different rates from the other states.
### This function is specific to the test to know two groups of states have different lambda and mu.
### chrmax, upper limit of chromosome number
### sp_state is the state numbers for different speciation rate (calculated with scal)

get_cons_and_target_musse_testSE_ec<-function(chrmax,sp_state){
  
  snum<-(chrmax+1)*(chrmax+2)/2-1
  
  #other_state is defined as the number of targetted other state from sp_state
  other_state<-1:snum
  other_state<-other_state[-sp_state]
  
  #repstate represents the representative state for the free parameter of the other states.
  repstate<-other_state[1]
  other_state<-other_state[-1]
  
  ord<-as.integer(log10(snum))+1
  ord2<-as.character(sprintf("%02d",ord))
  formula<-paste("lambda%",ord2,"d~lambda%",ord2,"d",sep="")
  cons<-as.list(sprintf(formula,other_state,repstate))
  formula<-paste("mu%",ord2,"d~mu%",ord2,"d",sep="")
  cons<-c(cons,as.list(sprintf(formula,other_state,repstate)))
  
  #when sp_state has more than one state 
  #the constraints within sp_state are made.
  if(length(sp_state)>1){
    f_sp_state<-sp_state[1]
    l_sp_state<-sp_state[-1]
    formula<-paste("lambda%",ord2,"d~lambda%",ord2,"d",sep="")
    cons<-append(cons,as.list(sprintf(formula,l_sp_state,f_sp_state)),l_sp_state-3)
    formula<-paste("mu%",ord2,"d~mu%",ord2,"d",sep="")
    cons<-append(cons,as.list(sprintf(formula,l_sp_state,f_sp_state)),snum+l_sp_state-5)
  }
  
  target.i<-1:(snum*2)
  target.i<-target.i[-c(repstate,sp_state[1],(snum+repstate),(snum+sp_state[1]))]
  
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


## Function, get_cons_and_target_musse_testSE_Ki1_ec
### M3 model, 7 parameters
### The constraint, Ki=1 is given to get_cons_and_target_musse_testSE_ec
### chrmax, upper limit of chromosome number
### sp_state is the state numbers for different speciation rate (calculated with scal)

get_cons_and_target_musse_testSE_Ki1_ec<-function(chrmax,sp_state){
  
  snum<-(chrmax+1)*(chrmax+2)/2-1
  other_state<-1:snum
  other_state<-other_state[-sp_state]
  repstate<-other_state[1]
  ord<-as.integer(log10(snum))+1
  ord2<-as.character(sprintf("%02d",ord))
  
  formula<-paste("lambda%",ord2,"d~lambda%",ord2,"d",sep="")
  cons<-as.list(sprintf(formula,other_state,repstate))
  formula<-paste("mu%",ord2,"d~mu%",ord2,"d",sep="")
  cons<-c(cons,as.list(sprintf(formula,other_state,repstate)))
  
  if(length(sp_state)>1){
    f_sp_state<-sp_state[1]
    l_sp_state<-sp_state[-1]
    formula<-paste("lambda%",ord2,"d~lambda%",ord2,"d",sep="")
    cons<-append(cons,as.list(sprintf(formula,l_sp_state,f_sp_state)),l_sp_state-3)
    formula<-paste("mu%",ord2,"d~mu%",ord2,"d",sep="")
    cons<-append(cons,as.list(sprintf(formula,l_sp_state,f_sp_state)),snum+l_sp_state-5)
  }
  
  target.i<-1:(snum*2)
  target.i<-target.i[-c(repstate,sp_state[1],(snum+repstate),(snum+sp_state[1]))]
  
  formula<-paste("q%",ord2,"d%",ord2,"d~%d*%s",sep="")
  
  #k4 is set as the same transition as k3
  kdel<-c(qwrite(ord,3,2),qwrite(ord,2,3),qwrite(ord,1,2),qwrite(ord,1,2)) #for Ki=1
  con_add<-function(cons,is,ts,kcoef,knum){
    return(c(cons,list(sprintf(formula,is,ts,kcoef,kdel[knum]))))
  }
  
  cons<-con_add(cons,2,1,1,3) #for Ki=1
  cons<-con_add(cons,3,4,2,3) #q[3,4]<- 2*k3
  cons<-con_add(cons,4,3,1,4) #q[4,3]<- k4
  cons<-con_add(cons,4,5,1,3) #q[4,5]<- k3
  cons<-con_add(cons,4,6,1,2) #q[4,6]<- k2
  cons<-con_add(cons,5,4,2,4) #q[5,4]<- 2*k4
  cons<-con_add(cons,5,7,2,2) #q[5,7]<- 2*k2
  targets<-c(tcal(snum,2,1),tcal(snum,3,4),tcal(snum,4,3),tcal(snum,4,5),tcal(snum,4,6),tcal(snum,5,4),tcal(snum,5,7)) #for Ki=1
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


## Function, get_cons_and_target_musse_polyp_ec
### M4 model, 7 parameters
### Almost same as get_cons_and_target_musse_null_ec 
### but a new transition for polyploidization was given.
### Only states with half chromosome number of the limit or less
### have this transition as doubling chromosome and arm number.
### The new transition, k5, is given as q001003.(i.e.(1,1)=>(2,2))
### chrmax, upper limit of chromosome number which should be more than 3.

get_cons_and_target_musse_polyp_ec<-function(chrmax){
  
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
  
  #k5 is given for the new transition for polyploidization as q001003
  kdel<-c(qwrite(ord,3,2),qwrite(ord,2,3),qwrite(ord,1,2),qwrite(ord,2,1),qwrite(ord,1,3))
  con_add<-function(cons,is,ts,kcoef,knum){
    return(c(cons,list(sprintf(formula,is,ts,kcoef,kdel[knum]))))
  }
  
  cons<-con_add(cons,2,5,1,5) #q[2,5]<- 1*k5
  cons<-con_add(cons,3,4,2,3) #q[3,4]<- 2*k3
  cons<-con_add(cons,3,10,1,5) #q[3,10]<- 1*k5
  cons<-con_add(cons,4,3,1,4) #q[4,3]<- k4
  cons<-con_add(cons,4,5,1,3) #q[4,5]<- k3
  cons<-con_add(cons,4,6,1,2) #q[4,6]<- k2
  cons<-con_add(cons,4,12,1,5) #q[4,12]<- 1*k5
  cons<-con_add(cons,5,4,2,4) #q[5,4]<- 2*k4
  cons<-con_add(cons,5,7,2,2) #q[5,7]<- 2*k2
  cons<-con_add(cons,5,14,1,5) #q[5,14]<- 1*k5
  targets<-c(tcal(snum,2,5),tcal(snum,3,4),tcal(snum,3,10),tcal(snum,4,3),
             tcal(snum,4,5),tcal(snum,4,6),tcal(snum,4,12),tcal(snum,5,4),
             tcal(snum,5,7),tcal(snum,5,14))
  half<-trunc(chrmax/2)
  if(half>=3){
    
    #chromosome number with 5 transitions (i.e. polyploidization happens)
    for(d in 3:half){
      
      s<-scal(d,d)
      down<-scal(d-1,d)
      poly<-scal(2*d,2*d)
      cons<-con_add(cons,s,down,d*(d-1)/2,1)#q[d*(d+1)/2,d*(d-1)/2+1]<-d*(d-1)*k1
      cons<-con_add(cons,s,s+1,d,3) #q[d*(d+1)/2,d*(d+1)/2+1]<-d*k3
      cons<-con_add(cons,s,poly,1,5) #q[d*(d+1)/2,2*d*(2*d+1)/2]<-1*k5
      targets<-c(targets,tcal(snum,s,down),tcal(snum,s,s+1),tcal(snum,s,poly))
      
      for(j in (d+1):(2*d-2)){
        s<-scal(d,j)
        up<-scal(d+1,j)
        down<-scal(d-1,j)
        poly<-scal(2*d,2*j)
        anum<-2*d-j
        mnum<-j-d
        cons<-con_add(cons,s,down,anum*(anum-1)/2,1)#q[k,down]<-anum*(anum-1)*k1
        cons<-con_add(cons,s,s-1,mnum,4) #q[k,k-1]<-mnum*k4
        cons<-con_add(cons,s,s+1,anum,3) #q[k,k+1]<-anum*k3
        cons<-con_add(cons,s,up,mnum,2)#q[k,up]<-mnum*k2
        cons<-con_add(cons,s,poly,1,5)
        targets<-c(targets,tcal(snum,s,down),tcal(snum,s,s-1),
                   tcal(snum,s,s+1),tcal(snum,s,up),tcal(snum,s,poly))
      }
      
      s<-scal(d,2*d-1)
      up<-scal(d+1,2*d-1)
      poly<-scal(2*d,4*d-2)
      cons<-con_add(cons,s,s-1,d-1,4) #q[k,k-1]<-(d-1)*k4
      cons<-con_add(cons,s,s+1,1,3) #q[k,k+1]<-k3
      cons<-con_add(cons,s,up,d-1,2)#q[k,up]<-(d-1)*k2
      cons<-con_add(cons,s,poly,1,5)
      targets<-c(targets,tcal(snum,s,s-1),tcal(snum,s,s+1),tcal(snum,s,up),tcal(snum,s,poly))
      
      s<-scal(d,2*d)
      up<-scal(d+1,2*d-1)
      poly<-scal(2*d,4*d)
      cons<-con_add(cons,s,s-1,d,4) #q[k,k-1]<-d*k4
      cons<-con_add(cons,s,up,d,2)#q[k,up]<-d*k2
      cons<-con_add(cons,s,poly,1,5)
      targets<-c(targets,tcal(snum,s,s-1),tcal(snum,s,up),tcal(snum,s,poly))
    }
    nextchr<-half+1
  }else{
    nextchr<-3
  }
  
  #chromosome number with 4 or less transitions
  for(d in nextchr:(chrmax-1)){
    
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
  
  #chromosome number with 3 or less transitions
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


## For constrain.kt (bellow), function get_exp
### Returns selected formulae from the list in lapply

get_exp<-function(form,ind){
  return(as.formula(form)[[ind]])
}

## Function, constrain.kt
### This function replace "constrain{diversitree}" function that produced errors in our analysis.
### You need to specify index of free parameters and target parameters as free.i and targt.i.
### free.i can be determined using tcal (see the example of running scripts)
### target.i is given by function of get_cons_and_target_~
### The other transition than free.i and target.i were set as 0.

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


###### Load phylogeny of 815 species of fishes ######

# This tree is corresponding to one used for Eurypterygii
# with chromoosme limit 35
# It includes 815 species
load("phy_o_neot_YK2021.Robj")

## Because the ultrametric tree has problem of the rounding when it was saved.
## That has slightly different distance of the terminal nodes from the ancestor.
## The loaded file should be rounded first.
## Just the terminal edges were trancated here using mean value.

tipnum<-Ntip(phy.o.neot)
edgeindex<-which(phy.o.neot$edge[,2]<=tipnum)
meanedge<-mean(node.depth.edgelength(phy.o.neot)[1:tipnum])
excesslen<-node.depth.edgelength(phy.o.neot)[1:tipnum]-meanedge
phy.o.neot$edge.length[edgeindex]<-phy.o.neot$edge.length[edgeindex]-excesslen
phy.sim<-phy.o.neot


###### Trial ######

for(trial in strial:etrial){
  
  ###### Trait simulation ######
  
  ## state number
  snum<-(chrmax+1)*(chrmax+2)/2-1
  
  ## determine 4 parameters by runif in log scale
  logcoef<-runif(4,-4,-1)
  coef_res<-10^logcoef
  names(coef_res)<-c("k1","k2","k3","k4")
  print(coef_res)
  
  ## make qmatrix
  qmatrix<-q_generation_2020(coef_res["k1"],coef_res["k2"],coef_res["k3"],coef_res["k4"],chrmax)
  
  ## simulate traits on the phylogeny
  sim_traits<-sim.character(phy.sim,qmatrix,x0=12,model="mkn")
  
  ## save file
  filename<-paste("sim_traits_EC_CM8_trial",trial,"_XXXXXX.Robj",sep="")
  save(sim_traits,file=filename)
  
  ## bind simulated traits on the phylogeny
  phy.sim$tip.state<-as.vector(sim_traits)
  names(phy.sim$tip.state)<-names(sim_traits)
  
  ###### Mk-n MLE #######
  
  ## make mkn object
  
  mknfunc<-make.mkn(phy.sim,phy.sim$tip.state,snum,strict=F)
  
  ## determine free parameters index
  free.i<-c(tcal(snum,1,2),tcal(snum,2,1),tcal(snum,2,3),tcal(snum,3,2))
  
  ## get constraint formulae and target.i for Mk-n
  cons_and_tar<-get_cons_and_target_ec(chrmax)
  
  ## make constraints
  mknfunc<-constrain.kt(mknfunc,cons_and_tar[[1]],free.i,cons_and_tar[[2]])
  
  ## Starting MLE
  p<-c(0.1,0.1,0.1,0.1) #prior
  print(paste("Start MLE of Mk-n of CM8 trial",trial))
  timepoint<-c(proc.time()[[3]])
  fit.mkn<-find.mle(mknfunc,p)
  timepoint<-c(timepoint,proc.time()[[3]])
  mkntime<-round(timepoint[2]-timepoint[1],3)
  
  ## Save the result object
  filename<-paste("fit_mkn_EC_CM8_trial",trial,"_XXXXXX.Robj",sep="")
  save(fit.mkn,file=filename)
  
  ## Output summary result in the text
  coef_res<-c(coef(fit.mkn)[c(4,3,1,2)]) #2020 definition
  names(coef_res)<-c("k1","k2","k3","k4")
  print(paste("Mkn (Process time",mkntime,")"))
  print(coef_res)
  
  ###### MLE MuSSE M0 model #######
  
  ## make musse object
  mussefunc<-make.musse(phy.sim,phy.sim$tip.state,snum,strict=F)
  
  ## determine free parameters index
  free.i<-c(tcal(snum,1,2),tcal(snum,2,1),tcal(snum,2,3),tcal(snum,3,2)) #kx
  free.i<-c(1,snum+1,((snum*2)+free.i)) # speciation and extinction rate and kx (shifted index)
  
  ## get constraint formulae and target.i for MuSSE M0 model
  cons_and_tar<-get_cons_and_target_musse_null_ec(chrmax)
  
  ## make constraint
  mussefunc<-constrain.kt(mussefunc,cons_and_tar[[1]],free.i,cons_and_tar[[2]])
  
  ## get prior
  int.p<-starting.point.musse(phy.sim,k=snum) #prior
  p<-int.p[free.i] #prior
  names(p)<-c()
  
  ## Start MLE
  print(paste("Start MLE of MuSSE of CM8 trial",trial))
  timepoint<-c(proc.time()[[3]])
  fit.musse<-find.mle(mussefunc,p)
  timepoint<-c(timepoint,proc.time()[[3]])
  mussetime<-round(timepoint[2]-timepoint[1],3)
  
  ## Save the result object
  filename<-paste("fit_musse_EC_CM8_trial",trial,"_XXXXXX.Robj",sep="")
  save(fit.musse,file=filename)
  
  ## Output summary result in the text
  coef_res<-c(coef(fit.musse)[c(6,5,3,4,1,2)]) #2020 definition
  names(coef_res)<-c("k1","k2","k3","k4","lambda","mu")
  print(paste("MuSSE (Process time",mussetime,")"))
  print(coef_res)

}
