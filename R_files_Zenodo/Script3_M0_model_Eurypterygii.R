###### Parameter ######
## Give chromosome number limit 
args = commandArgs(trailingOnly=TRUE)
chrmax<-as.numeric(args[1])

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


###### Preparation of the tree ######

## Load tree and karyotype data
phy.fish<-read.tree("Rabosky_et_al_timetree.tre")

## Rounding the tree
## Need to round the ultrametric tree before use
## Round the tree to 345 my, which is assumed in the original paper.
## Just terminal edges are trancated.

tipnum<-Ntip(phy.fish)
edgeindex<-which(phy.fish$edge[,2]<=tipnum)
excesslen<-node.depth.edgelength(phy.fish)[1:tipnum]-345.0000
phy.fish$edge.length[edgeindex]<-phy.fish$edge.length[edgeindex]-excesslen
names.fish<-phy.fish$tip.label
names.fish<-gsub("_"," ", names.fish)

## Binding names and karyotype states to the tip

## Phenotype are already filtered and arranged as neot.state.
## load the neot.state from the file

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
## This part is not important for M0 model.
## This calculates how many species were used in the analysis
## and how many species were listed in the karyoytpe data for each karyotype state.

## Karyotype and taxon information
## val_chrinfo contains all taxon and karyotype information of species used in the analyses.

load("val_chrinfo_YK2021.Robj")

neotlist<-c("AULOPIFORMES","MYCTOPHIFORMES","GADIFORMES","OPHIDIIFORMES",
            "MUGILIFORMES","ATHERINIFORMES","BELONIFORMES","CYPRINODONTIFORMES","STEPHANOBERYCIFORMES",
            "BERYCIFORMES","GASTEROSTEIFORMES","BATRACHOIDIFORMES","SYNBRANCHIFORMES","SCORPAENIFORMES",
            "PERCIFORMES","PLEURONECTIFORMES","LOPHIIFORMES","TETRAODONTIFORMES","ZEIFORMES","OPHIDIIFORMES")

## Determine sampling fruction, sampling.f
### sampling.f is not needed to be given in M0 model
### The calculation for "sfall" here is a coventional way to make sampling.f in our method.
### Because there are too many states,to decide sampling.f value of each state is impossible.
### We averaged sampling.f values for the states having the same speciation and extinction rates.
### Because M0 has constant speciation and extinction rates, all states have the same value as sfall.

val_chrinfo<-val_chrinfo[val_chrinfo$Order %in% neotlist,]
val_chrinfo$State<-factor(val_chrinfo$State,levels=1:snum)
val_chrinfo<-val_chrinfo[!(is.na(val_chrinfo$State)),]
sfall<-Ntip(phy.o.neot)/nrow(val_chrinfo)
sampling.f<-c(rep(sfall,snum))

###### MLE wih MuSSE M0 model ######
## Maximum Likelihood Estimation

## Make musse object with sampling.f
mussefunc<-make.musse(phy.o.neot,phy.o.neot$tip.state,snum,sampling.f=sampling.f,strict=F)

## Determine free.i for 4 parameters of karyotype evolution
free.i<-c(tcal(snum,1,2),tcal(snum,2,1),tcal(snum,2,3),tcal(snum,3,2))

## Determine free.i for 2 parameters of speciation and extinction rates
free.i<-c(1,snum+1,((snum*2)+free.i))

## Get constraint formulae and target.i for M0 model
cons_and_tar<-get_cons_and_target_musse_null_ec(chrmax)

## Make constraints
mussefunc<-constrain.kt(mussefunc,cons_and_tar[[1]],free.i,cons_and_tar[[2]])

## Get priors
int.p<-starting.point.musse(phy.o.neot,k=snum)
p<-int.p[free.i] #prior
names(p)<-c()

## Start MLE
print(paste("Start MLE of M0 neot CM",chrmax, "..."))
timepoint<-c(proc.time()[[3]])
fit.musse.neot<-find.mle(mussefunc,p)
timepoint<-c(timepoint,proc.time()[[3]])

## Save resultant fit.musse.neot R object
filename<-paste("fit_musse_neot_M0_CM",chrmax,"_unisample_XXXXXX.Robj",sep="")
save(fit.musse.neot,file=filename)

## Text output of the result
coef_res<-c(coef(fit.musse.neot)[c(6,5,3,4,1,2)]) #2020 definition
names(coef_res)<-c("k1","k2","k3","k4","lambda","mu")
coef_res #coef
fit.musse.neot[2] #lnLik
paste("Chromosome limit,",chrmax, ": Process time,", round(timepoint[2]-timepoint[1],3))


###### ASR wih MuSSE M0 model ######
## Ansestral State Reconstruction

print(paste("Start ASR of M0 neot CM",chrmax, "..."))
timepoint<-c(proc.time()[[3]])
st.nodes.neot<-asr.marginal(mussefunc,coef(fit.musse.neot))
timepoint<-c(timepoint,proc.time()[[3]])

## Save the resultant R objecct, st.nodes.neot
filename<-paste("st_nodes_neot_M0_CM",chrmax,"_unisample_XXXXXX.Robj",sep="")
save(st.nodes.neot,file=filename)

## Time
paste("Chromosome limit,",chrmax, ": Process time,", round(timepoint[2]-timepoint[1],3))

