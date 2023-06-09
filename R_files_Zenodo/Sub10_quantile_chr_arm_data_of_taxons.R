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


###### Data preparation for karyotype stat ######

chrinfo<-read.table("Teleostei_karyotype_data_YK2021.csv",header=T,sep=",",na.strings="")
load("neot_state_YK2021.Robj")#817
load("otop_state_YK2021.Robj")#519
load("val_index_YK2021.Robj")

delindex<-which(is.na(match(chrinfo$Species_name,names(neot.state)))==F)
delindex<-c(delindex,which(is.na(match(chrinfo$Species_name,names(otop.state)))==F))
chrinfo2<-chrinfo[-delindex,]
usedindex<-match(names(neot.state),chrinfo$Species_name)
usedindex<-c(usedindex,match(names(otop.state),chrinfo$Species_name))
val_chrinfo<-rbind(chrinfo2[val_index,],chrinfo[usedindex,])#2604
val_chrinfo$Chr_num<-round(val_chrinfo$Chr_num/2)
val_chrinfo$Arm_num<-round(val_chrinfo$Arm_num/2)
val_chrinfo$State<-scal(val_chrinfo$Chr_num,val_chrinfo$Arm_num)
val_chrinfo$State[1269:2604]<-c(neot.state,otop.state)

###### Chrmosome no. and Arm no. stat of each fish group (Fig S1) ######

## Binding Chromosome no. and Arm no. of Eurypterygii and Otophysi

neot_otop_chr_arm<-chr_arm_list(c(neot.state,otop.state))
val_chrinfo$Chr_num[1269:2604]<-neot_otop_chr_arm[[1]]
val_chrinfo$Arm_num[1269:2604]<-neot_otop_chr_arm[[2]]

## Quantile and species no. are checked of each order

order_list<-unique(val_chrinfo$Order)
restable<-c()
for(ord in 1:length(order_list)){
  spnumber<-sum(val_chrinfo$Order==order_list[ord])
  chrq<-quantile(val_chrinfo$Chr_num[val_chrinfo$Order==order_list[ord]])
  armq<-quantile(val_chrinfo$Arm_num[val_chrinfo$Order==order_list[ord]])
  restable<-rbind(restable,c(as.character(order_list[ord]),spnumber,chrq,armq))
}

## Difinition of special taxons

subd_oste<-order_list[5:6]
subd_elop<-order_list[7:8]
subd_otoc<-order_list[9:14]
supo_osta<-order_list[10:14]
seri_otop<-order_list[11:14]
subd_eute<-order_list[15:37]
eurypt<-order_list[19:37]
teleos<-order_list[5:37]
taxcat_names<-c("Subdivision Osteoglossomorpha","Subdivision Elopomorpha",
              "Subdivision Otocephala","Superorder Ostariophysi",
              "Series Otophysi","Subdivision Euteleostei",
              "Eurypterygia","Teleostei")
taxcat_orders<-list(subd_oste,subd_elop,subd_otoc,supo_osta,seri_otop,subd_eute,eurypt,teleos)

## Quantile and species no. are checked of each special taxon

for(i in 1:length(taxcat_orders)){
  cati<-which(is.na(match(val_chrinfo$Order,taxcat_orders[[i]]))==F)
  chrq<-quantile(val_chrinfo$Chr_num[cati])
  armq<-quantile(val_chrinfo$Arm_num[cati])
  restable<-rbind(restable,c(taxcat_names[[i]],length(cati),chrq,armq))
}

write.table(restable,file="Quantile_chr_arm_data_XXXXXX.txt",sep="\\t",quote=F,row.names=F)

