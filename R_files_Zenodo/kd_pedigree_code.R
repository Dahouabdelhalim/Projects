library(radiator)
library(vcfR)
library(sequoia)
library(stringi)
sca_info<-read.csv("~/Sistrurus_all_info.csv")
sca_strata<-sca_info[,c(2,10)]
colnames(sca_strata)<-c("INDIVIDUALS","STRATA")

colnames(kd_bad_ids)<-c("INDIVIDUALS")
kd_output<-genomic_converter("C:/Users/smart/Documents/kd_sca_2020.vcf", blacklist.id=kd_bad_ids,
                             strata=sca_strata, output=c('genind','tidy'), filename='kd_sca_2020', parallel.core=1)

names<-c("ID","Sex","Age","County","Locality","Field","Lat","Long","Day","Month","Year","Date","Mass_g","SVL_cm", "Tail_cm","Total_length","Num_rattlesegments","Rattlecomplete")
kd_genind<-kd_output$genind
kd_tidy<-kd_output$tidy.data
kd_ordered_genind<-kd_genind[order(kd_genind@pop),]
kd_info<-subset(sca_info, sca_info$Dna.. %in% row.names(kd_genind@tab))

#this trims down the info objects to only important data, and renames the columns
#then genind objects are reorganized to match order in info object
kd_info_trim<-kd_info[,c(2,4,5,8:16,34:39)]
colnames(kd_info_trim)<-names
kd_genind<-kd_genind[match(kd_info_trim$ID,kd_genind@strata$INDIVIDUALS),]



kd_info_trim$SVL_cm<-as.numeric(as.character(kd_info_trim$SVL_cm))
kd_prob_age<-cut(as.numeric(kd_info_trim$SVL_cm), breaks=c(-Inf,1,20,32,50,500,Inf), labels=0:5)
kd_prob_age<-as.numeric(kd_prob_age)-2
kd_prob_age[is.na(kd_prob_age)]<-0
kd_info_trim$prob_age<-kd_prob_age
kd_life_hist<-cbind(as.character(kd_info_trim$ID),as.character(kd_info_trim$Sex), as.numeric(kd_info_trim$Year)-as.numeric(kd_info_trim$prob_age))
#recodes for sequioa input where 1=female, 2=male, 3=unknown
for (i in 1:length(kd_life_hist[,1])){
  if (as.character(kd_life_hist[i,2])=="F")
    kd_life_hist[i,2]<-1
  else if (as.character(kd_life_hist[i,2])=="M")
    kd_life_hist[i,2]<-2
  else
    kd_life_hist[i,2]<-3
  if (is.na(kd_life_hist[i,3]))
    kd_life_hist[i,3]<- -1
}
kd_life_hist<-as.data.frame(kd_life_hist)



kd_output_thin<-read.vcfR("C:/Users/smart/Documents/kd_sca_2020_filtered.vcf")
loc_test<-vcfR::vcfR2loci(kd_output_thin)
nind<-length(loc_test[,1])
nloci<-length(loc_test[1,])*2
#Create new matrix and format it to calculate pedigrees for all individuals
relate<-matrix(NA, nrow=nind, ncol=nloci)
for (i in 1:nind){
  for (j in 1:((nloci/2))){
    loci<-substring(loc_test[i,j],1,3)
    a1<-substring(loci,1,1)
    a2<-stringi::stri_sub(loci,-1,-1)
    relate[i,((j-1)*2+1)]<-as.numeric(a1)+1
    relate[i,((j-1)*2+2)]<-as.numeric(a2)+1
  }
}
row.names(relate)<-row.names(loc_test)
kd_relate<-relate

kd_seq<-GenoConvert(InFile = kd_relate, InFormat = "col")
rownames(kd_seq)==kd_life_hist$V1
kd_life_hist<-subset(kd_life_hist, kd_life_hist$V1 %in% rownames(kd_seq))
kd_life_hist<-kd_life_hist[match(rownames(kd_seq),as.character(kd_life_hist$V1)),]

kd_pedi<-sequoia(GenoM=kd_seq, LifeHistData = kd_life_hist, MaxMismatch = 2,
                 FindMaybeRel = TRUE, CalcLLR = TRUE, UseAge = "yes", 
                 args.AP = list(Flatten = TRUE, Smooth = FALSE))




kd_info_trim<-read.csv("~/Killdeer_sca_2020_info.csv")
kd_info_trim<-kd_info_trim[,-1]
kd_info_trim<-subset(kd_info_trim, kd_info_trim$ID %in% row.names(kd_seq))
kd_info_trim$Lat<-as.numeric(kd_info_trim$Lat)+rnorm(length(kd_info_trim$Lat),mean=0, sd=0.000002)
kd_info_trim$Long<-as.numeric(kd_info_trim$Long)+rnorm(length(kd_info_trim$Long),mean=0, sd=0.000002)
kd_pedi_pairs<-kd_pedi$MaybeRel
kd_pedi_pairs[,5:6]<-as.character(kd_pedi_pairs[,5:6])
for (i in 1:length(kd_pedi_pairs[,1])){
  for (j in 1:length(kd_info_trim[,1])){
    if (kd_pedi_pairs[i,1]==kd_info_trim[j,1]){
      kd_pedi_pairs[i,5]<-as.character(kd_info_trim$Field[j])
      kd_pedi_pairs[i,7:8]<-c(kd_info_trim$Lat[j],kd_info_trim$Long[j])
    }
    if (kd_pedi_pairs[i,2]==kd_info_trim[j,1]){
      kd_pedi_pairs[i,6]<-as.character(kd_info_trim$Field[j])
      kd_pedi_pairs[i,9:10]<-c(kd_info_trim$Lat[j],kd_info_trim$Long[j])
    }
  }  
}
colnames(kd_pedi_pairs)<-c("ID1","ID2","TopRel","LLR","Loc1","Loc2","Lat1","Long1","Lat2","Long2")
kd_pedi_pairs$is_diff<-NA
for (i in 1:length(kd_pedi_pairs[,1])){
  if (kd_pedi_pairs[i,5] != kd_pedi_pairs[i,6]){
    kd_pedi_pairs[i,11]=TRUE
  }
  if (kd_pedi_pairs[i,5] == kd_pedi_pairs[i,6]){
    kd_pedi_pairs[i,11]=FALSE
  }
}
kd_pedi_pairs$geo_dist<-NA
kd_pedi_pairs$geo_dist<-as.numeric(distHaversine(kd_pedi_pairs[,7:8],kd_pedi_pairs[,9:10]))
kd_pedi_pairs$disp_num<-NA
for (i in 1:length(kd_pedi_pairs[,1])){
  if (kd_pedi_pairs[i,]$TopRel=="PO")
    kd_pedi_pairs[i,]$disp_num<-1
  if (kd_pedi_pairs[i,]$TopRel=="FS")
    kd_pedi_pairs[i,]$disp_num<-2
  if (kd_pedi_pairs[i,]$TopRel=="HS")
    kd_pedi_pairs[i,]$disp_num<-3
  if (kd_pedi_pairs[i,]$TopRel=="GP")
    kd_pedi_pairs[i,]$disp_num<-2
  if (kd_pedi_pairs[i,]$TopRel=="FA")
    kd_pedi_pairs[i,]$disp_num<-4
  if (kd_pedi_pairs[i,]$TopRel=="HA")
    kd_pedi_pairs[i,]$disp_num<-4
}

kd_pedi_pairs$disp_dist<-NA
kd_pedi_pairs$disp_dist<-kd_pedi_pairs$geo_dist/kd_pedi_pairs$disp_num



###########Code for conStruct
# Need to export vcf file as structure file to then convert into construct
kd_str<-as.matrix(read.delim("C:/Users/smart/Desktop/sca_data/kd_str_dec_2018.str", header=F))
kd_const<-structure2conStruct('C:/Users/smart/Desktop/sca_data/kd_str_dec_2018.str',onerowperind = F,start.loci = 3,missing.datum = "-9", outfile="kd_construct6")
kd_Dgeo<-as.matrix(dist(kd_info_const[,4:5], diag=T, upper=T))
kd_const<-kd_const[order(row.names(kd_const)),]
#quick check to ensure gps points match up to correct individuals
row.names(kd_const)==as.character(kd_info_const[,1])
kd_coords<- as.matrix(kd_info_const[,4:5])
kd_construct_1<-conStruct(spatial=T, K=1, freqs=kd_const,geoDist=kd_Dgeo,
                          coords = kd_coords,prefix="kd_construct_k1_", make.figs = T,
                          n.iter=10000)
kd_construct_2<-conStruct(spatial=T, K=2, freqs=kd_const,geoDist=kd_Dgeo,
                          coords = kd_coords,prefix="kd_construct_k2_", make.figs = T,
                          n.iter=10000)
kd_construct_3<-conStruct(spatial=T, K=3, freqs=kd_const,geoDist=kd_Dgeo,
                          coords = kd_coords,prefix="kd_construct_k3_", make.figs = T,
                          n.iter=10000)
kd_construct_4<-conStruct(spatial=T, K=4, freqs=kd_const,geoDist=kd_Dgeo,
                          coords = kd_coords,prefix="kd_construct_k4_", make.figs = T,
                          n.iter=10000)
kd_construct_5<-conStruct(spatial=T, K=5, freqs=kd_const,geoDist=kd_Dgeo,
                          coords = kd_coords,prefix="kd_construct_k5_", make.figs = T)
kd_construct_6<-conStruct(spatial=T, K=6, freqs=kd_const,geoDist=kd_Dgeo,
                          coords = kd_coords,prefix="kd_construct_k6_", make.figs = T)
kd_construct_7<-conStruct(spatial=T, K=7, freqs=kd_const,geoDist=kd_Dgeo,
                          coords = kd_coords,prefix="kd_construct_k7_", make.figs = T)

kd_layer.contributions <- matrix(NA,nrow=7,ncol=7)
load("C:/Users/smart/Documents/R/kd_construct_k1__conStruct.results.Robj")
load("C:/Users/smart/Documents/R/kd_construct_k1__data.block.Robj")
kd_layer.contributions[,1] <- c(calculate.layer.contribution(conStruct.results[[1]],data.block),rep(0,6))
tmp <- conStruct.results[[1]]$MAP$admix.proportions
for(i in 2:7){
  # load the conStruct.results.Robj and data.block.Robj
  #	files saved at the end of a conStruct run
  load(sprintf("C:/Users/smart/Documents/R/kd_construct_k%s__conStruct.results.Robj",i))
  load(sprintf("C:/Users/smart/Documents/R/kd_construct_k%s__data.block.Robj",i))
  
  # match layers up across runs to keep plotting colors consistent
  #	for the same layers in different runs
  tmp.order <- match.layers.x.runs(tmp,conStruct.results[[1]]$MAP$admix.proportions)	
  # calculate layer contributions
  kd_layer.contributions[,i] <- c(calculate.layer.contribution(conStruct.results=conStruct.results[[1]],
                                                               data.block=data.block,
                                                               layer.order=tmp.order),
                                  rep(0,7-i))
  tmp <- conStruct.results[[1]]$MAP$admix.proportions[,tmp.order]
}

#looks like k=2 is most supported, using construct. 2nd make struture orders on longitude
load("kd_construct_k2__conStruct.results.Robj")
load("kd_construct_k2__data.block.Robj")
admix.props <- conStruct.results$chain_1$MAP$admix.proportions
make.structure.plot(admix.proportions = admix.props,
                    sort.by = 1, sample.names=row.names(data.block$coords))
make.structure.plot(admix.proportions = admix.props,
                    sample.order = order(data.block$coords[,2]))
make.structure.plot(admix.proportions = admix.props,
                    sample.names = row.names(data.block$obsCov),
                    mar = c(4.5,4,2,2))
make.admix.pie.plot(admix.proportions = admix.props,
                    coords = data.block$coords,
                    x.lim = c(40.69,40.73),
                    y.lim = c(-83.35,-83.25), radii=3)
random_coords<-data.block$coords+rnorm(214, mean = 0, sd = 0.0025)
make.admix.pie.plot(admix.proportions = admix.props,
                    coords = random_coords,
                    x.lim = c(40.69,40.73),
                    y.lim = c(-83.35,-83.25), radii=3)