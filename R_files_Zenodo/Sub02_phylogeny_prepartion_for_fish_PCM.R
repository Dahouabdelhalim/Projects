###### Parameter ######
chrmax<-35
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


###### Preparation of trees with valid tips ######

## Load tree and karyotype data 
sp_name_get<-function(x){
  genus<-x[3]
  sec_species<-x[4]
  paste(genus,sec_species)
}
phy.fish<-read.tree("Rabosky_et_al_timetree.tre")
names.fish<-phy.fish$tip.label
names.fish<-gsub("_"," ", names.fish)

chrinfo<-read.table("Teleostei_karyotype_data_YK2021.csv",header=T,sep=",",na.strings="")
sum(!is.na(match(chrinfo$Species_name,names.fish)))

### 1441 species in the tree
### 140 duplications

## Drop tips not analyzed

delindex<-which(is.na(match(names.fish,chrinfo$Species_name))) #tip index of not karyotyped species
chrindex<-match(names.fish,chrinfo$Species_name) #index in chrinfo of each tip
outindex<-c()
outlist<-c("ACIPENSERIFORMES","LEPISOSTEIFORMES","AMIIFORMES","OSTEOGLOSSIFORMES",
           "ELOPIFORMES","ANGUILLIFORMES","CLUPEIFORMES","GONORYNCHIFORMES","OSMERIFORMES",
           "SALMONIFORMES","ESOCIFORMES","ARGENTINIFORMES") #Groups removed
for(i in 1:length(outlist)){
  outindex<-c(outindex,which(chrinfo$Order[chrindex]==outlist[i])) #tip index of species in each outlist
}
delindex<-c(outindex,delindex) #tip index removed
phy.overlapped<-drop.tip(phy.fish,delindex)
#alltree 1336 Species

## Confirmation
test.ind<-match(names.fish[-delindex],as.character(chrinfo$Species_name))
unique(chrinfo$Order[test.ind])

## Generate trees for two groups 
### neoteleostei-817
### otophysi-519

phy.o.neot<-drop.tip(phy.overlapped,818:1336)
phy.o.otop<-drop.tip(phy.overlapped,1:817)
names.over<-names.fish[-delindex]
names.o.neot<-names.over[1:817]
names.o.otop<-names.over[818:1336]

## Confirmation
test.ind<-match(names.o.neot,as.character(chrinfo$Species_name))
unique(chrinfo$Order[test.ind])
plot(phy.o.neot,cex=1)
test.ind<-match(names.o.otop,as.character(chrinfo$Species_name))
unique(chrinfo$Order[test.ind])
plot(phy.o.otop,cex=1)

###### Binding karyotype states on tree ######

## Random sampling when there are more than one report of one species in the karyotype data
index_choice<-function(x,y){
  c_ind<-c()
  match_list<-charmatch(x,y)
  for(i in 1:length(x)){
    if(match_list[i]==0){
      c_ind[i]<-sample(which(y==x[i]),1)
    }else{
      c_ind[i]<-match_list[i]
    }
  }
  c_ind
}
neot.ind<-index_choice(names.o.neot,chrinfo$Species_name) #index in chrinfo of each tip
otop.ind<-index_choice(names.o.otop,chrinfo$Species_name)

## Get haplid chromosome and arm number

neot.chr<-round(chrinfo$Chr_num[neot.ind]/2)
neot.arm<-round(chrinfo$Arm_num[neot.ind]/2)
otop.chr<-round(chrinfo$Chr_num[otop.ind]/2)
otop.arm<-round(chrinfo$Arm_num[otop.ind]/2)

## Convert to karyotype state number

neot.state<-scal(neot.chr,neot.arm)
otop.state<-scal(otop.chr,otop.arm)

## Binding names to the karyotype states

names(neot.state)<-names.o.neot
names(otop.state)<-names.o.otop
save(neot.state,"neot_state_XXXXXX.Robj")#817
save(otop.state,"otop_state_XXXXXX.Robj")#519

## Binding names and karyotype states to the tip

phy.o.neot$tip.label<-names.o.neot
phy.o.neot$tip.state<-neot.state
phy.o.otop$tip.label<-names.o.otop
phy.o.otop$tip.state<-otop.state

## Remove the tip with out of range in PCM

snum<-(chrmax+1)*(chrmax+2)/2-1
phy.o.neot<-drop.tip(phy.o.neot,which(neot.chr>chrmax))
phy.o.neot$tip.state<-phy.o.neot$tip.state[which(neot.chr<=chrmax)]
phy.o.otop<-drop.tip(phy.o.otop,which(otop.chr>chrmax))
phy.o.otop$tip.state<-phy.o.otop$tip.state[which(otop.chr<=chrmax)]

save(phy.o.neot,"phy_o_neot_XXXXXX.Robj")
save(phy.o.otop,"phy_o_otop_XXXXXX.Robj")