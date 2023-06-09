## Use R version >3.5.0, otherwise it takes time in mle
## This script takes time in certain parts annotated by (!).
## It's smart to save R object in each part.

###### Library ######

library(diversitree)
library(expm)
library(Matrix)
library(stats4)

## The parts with (!) need much prcess time.
## Once you get resultant R objects, you should skip these parts.
## You can load the R object used in the manuscript and skip these parts.

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



###### Data preparation for Eurypterygii ######

chrmax<-35
load("phy_o_otop_YK2021.Robj")
load("fit_musse_otop_M2_YK2021.Robj")
load("st_nodes_otop_M2_YK2021.Robj")

phy<-phy.o.otop
fit.musse<-fit.musse.otop
st.nodes<-st.nodes.otop

snum<-(chrmax+1)*(chrmax+2)/2-1

#make q matrix
k1<-coef(fit.musse)[8] #2020 definition
k2<-coef(fit.musse)[7]
k3<-coef(fit.musse)[5]
k4<-coef(fit.musse)[6]
qmat<-q_generation_2020(k1,k2,k3,k4,chrmax)

###### Matrix and vector for unique transitions ######
## tvec and tindseq produced here are important through this script until the end.
## A elements of edge_len2 later is corresponding to the tvec.
## The same type of transition is fitted only once to reduce the process time.
## tvec is a vector for the unique transition
## The transition is numbered in order of q matrix including diagonal. (different from tcal)
## tindseq shows index of tvec corresponding to edges.
## If no transition happens in the edge, input 0 to the entry of tindseq

tvec<-c()
tindseq<-rep(0,nrow(phy$edge))
for(i in 1:nrow(phy$edge)){
  inode<-phy$edge[i,1]
  fnode<-phy$edge[i,2]
  length(st.nodes)
  maxprob<-max(st.nodes[,inode-Ntip(phy)])
  
  # mode of the distribution of acendant node is used as the acendant state, is
  is<-which(st.nodes[,inode-Ntip(phy)]==maxprob)
  if(length(is)>1){
    print(paste("More than one candidate of i state",inode,is))
    is<-sample(is,1)
  }
  if(fnode <= Ntip(phy)){
    # descendant state, fs, is defined simly as tip state
    fs<-phy$tip.state[fnode]
  }else{
    maxprob<-max(st.nodes[,fnode-Ntip(phy)])
    # mode of the distribution of descendant node is used as the descendant state, fs
    fs<-which(st.nodes[,fnode-Ntip(phy)]==maxprob)
  }
  if(length(fs)>1){
    # mostly, this doesn't happen
    print(paste("More than one candidate of f state",fnode,fs))
    fs<-sample(fs,1)
  }
  if(is==fs){
    tindseq[i]<-0
  }else{
    #transition id is given for each tvec element
    t.on.edge<-(fs-1)*snum+is
    tind<-which(tvec==t.on.edge)
    if(length(tind)>0){
      # when tvec has an element of the same transition
      tindseq[i]<-tind
    }else{
      # when no element, newly input the transition to tvec
      tvec<-c(tvec,t.on.edge)
      tindseq[i]<-length(tvec)
    }
  }
}

print("Complete tindmat")

###### (!) Taylor series generation ######
## This is only for generation of edge_len2.
## You can skip and reduce process time once you get the resultant R object, edge_len2.
## This part produces qc_mat, which represent 170 coefficient of Taylor series of each unique transition
## qmatc170 needs much memory space. I remove qmatc170 after the output of qc_mat

## Making list of matrixes for coefficients of all transitions of each term of Talyor series
qmatc170<-list()
for(n in 1:170){
  print(n)
  qmatc170<-c(qmatc170,(qmat %^% n)/(factorial(n)))
}

## Function to get corresponidng 170 coefficients from a transtion (is -> fs)
### qmatind here is a vector of indices of unlisted qmatc170
qc_extract<-function(is,fs,snum){
  tr<-(fs-1)*snum+is
  qmatind<-rep(tr,170)+(0:169)*length(qmat)
  unlist(qmatc170[qmatind])
}

## Generate qc_mat
qc_mat<-matrix(0,nrow=length(tvec),ncol=170)
for(i in 1:length(tvec)){
  is<-((tvec[i]-1) %% snum)+1
  fs<-((tvec[i]-1) %/% snum)+1
  qc_mat[i,]<-qc_extract(is,fs,snum)
}

## Remove qmatc170
remove(qmatc170)

print("Complete qc_mat")

###### (!) Fitting of branch length to each unique transition ######
## This is only for generation of edge_len2.
## You can skip and reduce process time once you get the resultant R object, edge_len2.

edge_len<-c()
errors<-c()
for(tr in 1:length(tvec)){
  print(paste(tr,((tvec[tr]-1) %% snum)+1,((tvec[tr]-1) %/% snum)+1))
  is<-((tvec[tr]-1) %% snum)+1
  fs<-((tvec[tr]-1) %/% snum)+1
  qc<-qc_mat[tr,]
  
  ## Function of log likelihood is defined here using qc for a given transition
  ## As explained in the manuscript, there is conditional control to calculate likelihood
  logl.cal.branch2 <-function(brlen){
    if(brlen>0 && brlen<15000){
      if(brlen >= 10){
        -log(expm(brlen*qmat,method="Higham08.b")[is,fs])
      }else{
        -log(sum(qc*((brlen)^(1:170))))
      }
    }else{
      NA
    }
  }
  
  p<-list(brlen=10)
  timepoint<-c(proc.time()[[3]])
  fit.branch<-try(mle(logl.cal.branch2,start=p))
  timepoint<-c(timepoint,proc.time()[[3]])
  print(paste("Process time", round(timepoint[2]-timepoint[1],3))) #see the process time
  if(inherits(fit.branch, "try-error")){
    errors<-c(errors,tr)
    edge_len[tr]<-NA
  }else{
    edge_len[tr]<-coef(fit.branch)[1]
    print(edge_len[tr])
  }
}

## Better to save edge_len here to use later.
filename<-paste("edge_len_otop_XXXXXX.Robj",sep="")
save(edge_len,file=filename)
print(paste("Error:",paste(errors,collapse=" ")))

###### Check of difference in lnL of ML between different fitting method ######

lnL_dif<-c()
for(i in 1:length(tvec)){
  brlen<-edge_len[i]
  if(brlen<10){
    is<-((tvec[i]-1) %% snum)+1
    fs<-((tvec[i]-1) %/% snum)+1
    lnL_qcmat<- -log(sum(qc_mat[i,]*((brlen)^(1:170))))
    lnL_expm<- -log(expm(brlen*qmat,method="Higham08.b")[is,fs])
    lnL_dif[i]<-(lnL_qcmat-lnL_expm)/lnL_expm
    print(paste(i,brlen,lnL_qcmat,lnL_expm,lnL_dif[i]))
  }else{
    lnL_dif[i]<-0
  }
}
max(lnL_dif)
#6.279472e-15
min(lnL_dif)
#-8.264897e-14

###### (!) Preparation of cumulative prob with actual time ######
## For categorizing the edges
## actual edge length of the phylogeny is used here.
## edge_len is not used here.

#Data:cumseq
#index--edge of phylogeny
#Value--cumulative probability

cumseq<-c()
for(i in 1:nrow(phy$edge)){
  print(i)
  #extract the information of given branch from joint asr
  inode<-phy$edge[i,1]
  fnode<-phy$edge[i,2]
  maxprob<-max(st.nodes[,inode-Ntip(phy)])
  is<-which(st.nodes[,inode-Ntip(phy)]==maxprob)
  
  if(length(is)>1){
    print(paste("More than one candidate of i state",inode,is))
    is<-sample(is,1)
  }
  if(fnode <= Ntip(phy)){
    fs<-phy$tip.state[fnode]
  }else{
    maxprob<-max(st.nodes[,fnode-Ntip(phy)])
    fs<-which(st.nodes[,fnode-Ntip(phy)]==maxprob)
  }
  if(length(fs)>1){
    print(paste("More than one candidate of f state",fnode,fs))
    fs<-sample(fs,1)
  }
  
  # Calculation of the probability series of transition from is.
  realq<-expm(qmat*phy$edge.length[i])
  probvec<-realq[is,]
  cump<-0
  probvec2<-probvec
  
  # Culculate the cumulative probability
  #probvec keeps all elements.
  #probvec2 loses elements summed.
  while(length(probvec2) >= 1){
    minp<-min(probvec2)
    minfs<-which(probvec==minp)
    minfs2<-which(probvec2==minp)
    probvec2<-probvec2[-minfs2]
    if(fs==minfs){
      cumprobvec<-cump
    }
    cump<-cump+sum(probvec[minfs])
  }
  cumseq<-c(cumseq,cumprobvec)
}

## Better to save cummat here to use later.
save(cumseq,file="cumseq_otop_XXXXXX.Robj")

###### (!) Preparation of cumulative prob matrix with fitted time ######
## For categorizing the edges
## edge_len2 is used here but cumseq is not used here.

cum_fitbr_tvec<-c()
for(tr in 1:length(tvec)){
  print(tr)
  is<-((tvec[tr]-1) %% snum)+1
  fs<-((tvec[tr]-1) %/% snum)+1
  realq<-expm(qmat*edge_len[tr])
  # Calculation of the probability.
  probvec<-realq[is,]
  cump<-0
  probvec2<-probvec
  #probvec keeps all elements.
  #probvec2 loses elements summed.
  while(length(probvec2) >= 1){
    minp<-min(probvec2)
    minfs<-which(probvec==minp)
    minfs2<-which(probvec2==minp)
    probvec2<-probvec2[-minfs2]
    if(fs==minfs){
      cumprobvec<-cump
    }
    cump<-cump+sum(probvec[minfs])
  }
  cum_fitbr_tvec<-c(cum_fitbr_tvec,cumprobvec)
}

cumseq_fitbr<-rep(1,length(tindseq))
for(i in 1:length(tindseq)){
  if(tindseq[i]!=0){
    cumseq_fitbr[i]<-cum_fitbr_tvec[tindseq[i]]
  }
}

## Better to save cummat_fitbr here to use later.
save(cumseq_fitbr,file="cumseq_fitbr_otop_XXXXXX.Robj")



###### Applying fitted branch length on phylogeny ######
## From this part, the process time is not problematic.
## Here elenmat and phy_fitted_br are produced.
## elenseq is a vector of fitted branch length of each edge of phylogeny.
## The phylogeny with fitted branch length is produced as phy_fitted_br
## phy_fitted_br should be saved to visualize later.

elenseq<-tindseq
for(i in 1:length(tindseq)){
    if(tindseq[i]!=0){
      elenseq[i]<-edge_len[tindseq[i]]
    }
}

phy_fitted_bl<-phy
phy_fitted_bl$edge.length<-elenseq

## Better to save phy_fitted_bl here to use later.
save(phy_fitted_bl,file="phy_fitted_bl_otop_XXXXXX.Robj")


###### Classification of branches into 4 category (S3 Table) #######
## The fitted branches of each edge is categorized 

catseq<-c()
for(e in 1:length(phy$edge.length)){
  cat<-0
  if(cumseq[e] >= 0.01){
    cat<-1
  }else{
    if(cumseq_fitbr[e] >= 0.01){
      if(elenseq[e]>phy$edge.length[e]){
        cat<-2
      }else{
        cat<-3
      }
    }else{
      cat<-4
    }
    
  }
  catseq<-c(catseq,cat) 
}

## Statistic test for the 4 category

terminal<-which(phy$edge[,2] <= Ntip(phy))
internal<-which(phy$edge[,2] > Ntip(phy))
inor<-sum(catseq[internal]==1)
ilon<-sum(catseq[internal]==2)
isho<-sum(catseq[internal]==3)
iunu<-sum(catseq[internal]==4)
tnor<-sum(catseq[terminal]==1)
tlon<-sum(catseq[terminal]==2)
tsho<-sum(catseq[terminal]==3)
tunu<-sum(catseq[terminal]==4)
c(inor,ilon,isho,iunu,tnor,tlon,tsho,tunu)
c(ilon,iunu+inor+isho,tlon,tunu,tnor+tsho)
fisher.test(matrix(c(inor,ilon+isho+iunu,tnor,tlon+tsho+tunu),nrow=2))
fisher.test(matrix(c(ilon,inor+isho+iunu,tlon,tnor+tsho+tunu),nrow=2))
c(isho,inor+ilon+iunu,tsho,tnor+tlon+tunu)
fisher.test(matrix(c(iunu,inor+ilon+isho,tunu,tnor+tlon+tunu+tsho),nrow=2))
fisher.test(matrix(c(isho,inor+ilon+iunu,tsho,tnor+tlon+tunu+tunu),nrow=2))
binom.test(isho+tsho,nrow(phy$edge),0.01*0.99/2)

## Preparation of edge_col used in the tree drawing

edge_col<-rep("black",length=length(phy$edge))
edge_col[catseq==3]<-rep("dodgerblue3")
edge_col[catseq==2]<-rep("gold3")
edge_col[catseq==4]<-rep("firebrick1")

#Better to save edge_col2 for using to visualize later.
save(edge_col,file="edge_col_otop_XXXXXX.Robj")

###### Stat_of_long_branches (S4 Table) #######

top10fitbrind<-order(elenseq,decreasing=T)[1:10]
tbdata<-data.frame(edge=top10fitbrind)
tbdata$fit<-elenseq[top10fitbrind]
tbdata$time<-phy$edge.length[top10fitbrind]
tbdata$category<-catseq[top10fitbrind]
ts<-tvec[tindseq[top10fitbrind]]
is<-((ts-1) %% snum)+1
fs<-((ts-1) %/% snum)+1
ik<-sapply(is,karyotype_cal)
fk<-sapply(fs,karyotype_cal)
tbdata$karyotype<-paste(ik,fk)

names(phy$tip.state[phy$edge[top10fitbrind,2]])

## Output the table of top 10 species with the longest branch length
write.table(tbdata,file="Table_top10_species_otop_XXXXXX.txt",row.names=F,quote=F,sep="\\t")
