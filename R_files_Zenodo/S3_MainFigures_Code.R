####CODE

source("S3_MainFigures_Setup_200826.R")

####HARD SOFT BREAKDOWN----
master_taxon_list<-c()
perc_basal<-c()
for(all in names(allwebs_comb)) {
  
  selected_loc<-allwebs_comb[[all]]  
  
  fw_data<-selected_loc[["FW"]]
  fw_nodes<-selected_loc[["NODE"]]
  
  s.nodes <- traits(fw_nodes,name=ifelse(all=="Caiman","Cayman",all))
  
  s.nodes <- assign.p(s.nodes, p.eql, p.hard, p.med, p.soft, vary=FALSE)
  
  #Calculate percent primary producers
  primprods <- fw_data[,colSums(fw_data) == 0]
  primprods <- (as.character(colnames(primprods)))
  perc_basal<-rbind(perc_basal,cbind(perc_basal=(100*length(primprods)/nrow(fw_data)),web=all))
  
  master_taxon_list<-rbind(master_taxon_list,cbind(s.nodes[,c("Species","HardSoft")],all))
}
perc_basal<-as.data.frame(perc_basal)
perc_basal$perc_basal<-as.numeric(as.character(perc_basal$perc_basal))

library(dplyr)
master_taxon_list_HS<-master_taxon_list %>% 
  group_by(all) %>%
  count(counting=HardSoft) %>% 
  mutate(counting = recode(counting, 
                           "1" = "Soft-bodied",
                           "2" = "Intermediate",
                           "3" = "Hard-bodied"))
master_taxon_list_counts<-as.data.frame(rbind(cbind(as.matrix(master_taxon_list_HS),type="Body type")))
master_taxon_list_counts$n<-as.numeric(as.character(master_taxon_list_counts$n))
#As percentages
master_taxon_list_counts <- as.data.frame(master_taxon_list_counts %>% group_by(all,type) %>% mutate(percent = n/sum(n)))

master_taxon_list_counts$foss<-ifelse(is.element(master_taxon_list_counts$all,names(fossil_comb)),"fossil_web","modern_web")
master_taxon_list_counts$all<-ordered(master_taxon_list_counts$all, levels = web.agesize.order)
foss_m <- ifelse(master_taxon_list_counts$foss == "fossil_web", "red", "blue")

lev_foss<-ifelse(is.element(size.order,names(fossil_comb)),"darkorange2","darkolivegreen4")

#Put legend in order
master_taxon_list_counts$counting<-ordered(master_taxon_list_counts$counting, levels = c("Soft-bodied", "Intermediate", "Hard-bodied"))
master_taxon_list_counts$all2<-recode(master_taxon_list_counts$all,Caiman="Cayman")

####NETWORK WIDE STATS----

OG_full_net_stats<-c()
#Loop different foodwebs
for(all in names(allwebs_comb)) {
  
  selected_loc<-allwebs_comb[[all]]  
  
  fw_data<-selected_loc[["FW"]]
  fw_nodes<-selected_loc[["NODE"]]
  
  s.nodes <- traits(fw_nodes,name=ifelse(all=="Caiman","Cayman",all))
  
  s.nodes <- assign.p(s.nodes, p.eql, p.hard, p.med, p.soft, vary=FALSE)
  
  mean_BodyType<-mean(s.nodes$p.hs)
  
  graphed_list<-list()
  graphed_list[[1]]<-as.data.frame(fw_data)
  
  OG.metrics   <- calc.metrics(graphed_list,1)
  chosen_row<-subset(OG.metrics,type=="Original")
  OG_full_net_stats<-rbind(OG_full_net_stats,cbind(chosen_row,mean_BodyType,web=all))
  print(all)
}



OG_full_net_stats2<-reshape::melt(OG_full_net_stats,id=c("run","type","web"))
OG_full_net_stats2$type<-ifelse(is.element(OG_full_net_stats2$web,names(fossil_comb)),"fossil_web","modern_web")

net.sizes<-subset(OG_full_net_stats2,variable=="S")
net.size<-net.sizes[,c("web","value")]

net.tr.size<-c()
for(i in unique(names(allwebs_comb))) {
  t<-c(web=i,
       trcon=(calc.C(make.trophic.species.web(as.matrix(allwebs_comb[[i]][["FW"]])))),
       trsize=(nrow(make.trophic.species.web(as.matrix(allwebs_comb[[i]][["FW"]]))))
       ) 
  net.tr.size<-rbind(net.tr.size,t)
}

####NODE LEVEL METRICS----

OG_web_stats<-c()

#Loop different foodwebs
for(all in names(allwebs_comb)) {
  
  selected_loc<-allwebs_comb[[all]]  
  
  fw_data<-selected_loc[["FW"]]
  fw_nodes<-selected_loc[["NODE"]]
  
  #Change name to Cayman to match with preservation groups (where necessary)
  s.nodes <- traits(fw_nodes,name=ifelse(all=="Caiman","Cayman",all))
  
  s.nodes <- assign.p(s.nodes, p.eql, p.hard, p.med, p.soft, vary=FALSE)
  
  graphed<-graph_from_adjacency_matrix(as.matrix(fw_data),mode="directed")
  
  avlink<-sum(fw_data)/nrow(fw_data)
  
  btw<-as.data.frame(betweenness(graphed,directed=TRUE,normalized=TRUE))

  deg_dist_in<-as.data.frame(degree(graphed,mode=c("in")))
  deg_dist_out<-as.data.frame(degree(graphed,mode=c("out")))
  
  tl <- NetIndices::TrophInd(as.matrix(fw_data))
  swtl<-calc.swtl(fw_data)
  
  oi<-calc.swoi(fw_data)
  
  OG_size<-as.numeric(as.character(nrow(s.nodes)))
  
  web_list<-cbind(s.nodes,
                  btw=as.numeric(btw[,1]),
                  deg_in=as.numeric(deg_dist_in[,1]),
                  deg_out=as.numeric(deg_dist_out[,1]),
                  deg_in_norm=as.numeric(deg_dist_in[,1]/avlink),
                  deg_out_norm=as.numeric(deg_dist_out[,1]/avlink),
                  TL=tl$TL,
                  OI=oi,
                  SWTL=swtl[,c("SWTL")],
                  web=all,
                  size=OG_size)
  rownames(web_list)<-NULL

  OG_web_stats<-rbind(OG_web_stats,web_list)
  
  print(paste(all))
  print(Sys.time())
}

#Plot plain stats
OG_web_stats_reduced<-OG_web_stats[,c("Species","HardSoft","btw","deg_in","deg_out","deg_in_norm","deg_out_norm",
                                      "TL","OI","SWTL","web","size")]

OG_web_stats_melt<-reshape2::melt(OG_web_stats_reduced,id=c("Species","HardSoft","web"))

OG_HS<-OG_web_stats_melt[,c("Species","HardSoft","web","variable","value")]
OG_HS$SpAtt<-"HardSoft"
colnames(OG_HS)<-c("Species","level","web","variable","value","SpAtt")
OG_test2<-rbind(OG_HS)

##Compare confidence intervals
OG_web_stats_reduced<-OG_web_stats[,c("Species","HardSoft", "btw","deg_in_norm","deg_out_norm",
                                      "TL","OI","SWTL","web")]
#Need to change column names to do.call below
colnames(OG_web_stats_reduced)<-c("Species","HardSoft", "btw","deg_in_norm","deg_out_norm",
                                  "TL","OI","SWTL","web")
OG_web_stats_melt<-reshape2::melt(OG_web_stats_reduced,id=c("Species","HardSoft","web"))
OG_web_stats_melt_HS<-OG_web_stats_melt[,c("Species","HardSoft","web","variable","value")]
colnames(OG_web_stats_melt_HS)<-c("Species","Level","web","variable","value")
OG_web_stats_melt_HS$Type<-"HardSoft"
OG_web_stats_melt_joinedBPHS<-rbind(OG_web_stats_melt_HS)
OG_web_stats_melt_joinedBPHS$Type<-as.factor(OG_web_stats_melt_joinedBPHS$Type)

library(pairwiseCI)
#apcdata<-subset(OG_web_stats_melt_joinedBPHS,Level!=2)
apcdata<-subset(OG_web_stats_melt_joinedBPHS, variable=="btw" | variable=="deg_in_norm" | variable=="deg_out_norm" | variable=="OI" | variable=="SWTL")
apcdata<-as.data.frame(merge(apcdata,subset(SpAtt_names,type=="HardSoft"),by.x=c("Level"),by.y=c("level"),all.x=TRUE))
apcdata$NewLevel <- factor(apcdata$NewLevel, levels = c("Soft-bodied","Mixed","Hard-bodied"))

apc<-pairwiseCI(value~NewLevel,data=apcdata,by=c("web","variable","Type"),alternative = "two.sided", conf.level = 0.95, method = "Param.diff")
apc2<-as.data.frame(print(apc))
#plot(apc)
#Add column describing data
apc3 <- as.data.frame(cbind(apc2,(do.call('rbind', strsplit(as.character(apc2$by),'.',fixed=TRUE)))))
colnames(apc3)<-c("estimate","lower","upper","comparison","by","web","variable","Type")
#Color if overlapping 0
apc3$overlap<-(apc3$lower<=0) & (apc3$upper>=0)
apc3$SigDif<-ifelse(apc3$overlap==TRUE,"No sig. dif.","Sig. dif.")

apc3$comparison<-as.factor(gsub("Mixed","Intermediate",apc3$comparison))
apc3<-merge(apc3,NodeAtt_names,by.x="variable",by.y=1,all.x=TRUE)
apc3$NewNodeAtt<-factor(apc3$NewNodeAtt, levels = as.character(NodeAtt_names[,2]))
apc3$web<-factor(apc3$web,levels=web.agesize.order)


apcdata2 <- merge(apcdata, NodeAtt_names, by.x = "variable", by.y = "OldNodeAtt", all.x = TRUE)
apcdata2$web<-factor(apcdata2$web,levels=web.agesize.order)
apcdata2$NewLevel<-gsub("Mixed","Intermediate",apcdata2$NewLevel)
apcdata2$NewLevel<-as.factor(apcdata2$NewLevel)
apcdata2$NewLevel<-factor(apcdata2$NewLevel,levels=rev(levels(apcdata2$NewLevel)))
apcdata2$NewNodeAtt <- factor(apcdata2$NewNodeAtt, levels = NodeAtt_names$NewNodeAtt)


#Drop outliers for boxplot
apcdata3 <- apcdata2 %>%
  group_by(NewNodeAtt, web, NewLevel) %>%
  summarise(boxplot = list(setNames(
    boxplot.stats(value)$stats,
    c("lower_whisker", "lower_hinge", "median", "upper_hinge", "upper_whisker")
  ))) %>%
  unnest_wider(boxplot)

apcdata2_new<-apcdata2
apcdata2_new$web<-recode(apcdata2_new$web,Caiman="Cayman")
apcdata2_new$NewLevel<-recode(apcdata2_new$NewLevel,"Soft-bodied"="Soft","Intermediate"="Int.","Hard-bodied"="Hard")

apcdata3_new<-apcdata3
apcdata3_new$web<-recode(apcdata3_new$web,Caiman="Cayman")
apcdata3_new$NewLevel<-recode(apcdata3_new$NewLevel,"Soft-bodied"="Soft","Intermediate"="Int.","Hard-bodied"="Hard")



####TROPHIC BREAKDOWN----
trophic_list<-c()
for(all in names(allwebs_comb)) {
  
  selected_loc<-allwebs_comb[[all]]  
  
  fw_data<-selected_loc[["FW"]]
  
  OG_size<-as.numeric(as.character(nrow(fw_data)))
  fw_data<-make.trophic.species.web(fw_data)
  NEW_size<-as.numeric(as.character(nrow(fw_data)))
  
  trophic_list<-rbind(trophic_list,cbind(OG_size,NEW_size,all))
}

trophic_list<-as.data.frame(trophic_list)
trophic_list$OG_size<-as.numeric(as.character(trophic_list$OG_size))
trophic_list$NEW_size<-as.numeric(as.character(trophic_list$NEW_size))

trophic_list$rat<-trophic_list$NEW_size/trophic_list$OG_size

trophic_list$foss<-factor(ifelse(is.element(trophic_list$all,names(fossil_comb)),"Ancient","Modern"))
trophic_list$all<-factor(trophic_list$all, levels = web.agesize.order)
trophic_list$all<-recode(trophic_list$all,Caiman="Cayman")


#First calculate TS and OG size for all webs
tm.list<-c()
for(all in names(allwebs_comb)) {
  
  selected_loc<-allwebs_comb[[all]]  
  
  fw_data<-selected_loc[["FW"]]
  fw_nodes<-selected_loc[["NODE"]]
  
  s.nodes <- traits(fw_nodes,name=ifelse(all=="Caiman","Cayman",all))
  s.nodes2<-s.nodes[,c("Species","HardSoft")]
  s.nodes2$Species<-as.character(rownames(fw_data))
  
  OG_web_size_TrophicSpecies<-calc.troph.species.LIST(fw_data,overlap=0.0001)
  
  troph.merge<-merge(OG_web_size_TrophicSpecies,s.nodes2,by.x="X1",by.y="Species",all.x=TRUE)
  troph.merge2<-merge(troph.merge,s.nodes2,by.x="X2",by.y="Species",all.x=TRUE)
  
  troph.merge2<- troph.merge2 %>% rowwise() %>% mutate(HardSoftMin = min(HardSoft.x, HardSoft.y))
  troph.merge2<- troph.merge2 %>% rowwise() %>% mutate(HardSoftMax = max(HardSoft.x, HardSoft.y))
  
  
  tm.BodyType<-as.data.frame(troph.merge2) %>% 
    group_by(HardSoftMin,HardSoftMax) %>% 
    summarise(x=mean(overlapped_combined,na.rm=TRUE))%>%
    mutate(web=all,type="BodyType")
  colnames(tm.BodyType)<-c("min","max","x","web","type")
  
  tm.list<-rbind(as.data.frame(tm.list),as.data.frame(tm.BodyType))
  
  print(all)
}


tm.list2<-as.data.frame(tm.list)
tm.list2$grp<-paste(tm.list2$min,tm.list2$max,sep="_")

tm.list2$foss<-ifelse(is.element(tm.list2$web,names(fossil_comb)),"fossil_web","modern_web")
tm.list2<-subset(tm.list2, type=="BodyType")

tm.list3<-tm.list2 %>% 
  mutate(grp2 = recode(grp, 
                       "1_1" = "Soft",
                       "1_2" = "S-I",
                       "2_2" = "Int.",
                       "2_3" = "I-H",
                       "3_3" ="Hard",
                       "1_3" = "H-S"))
tm.list3$grp2<-factor(tm.list3$grp2, levels = c("Soft","S-I","Int.","I-H","Hard","H-S"))
tm.list3$web<-factor(tm.list3$web, levels = web.agesize.order)
tm.list3$web<-recode(tm.list3$web,Caiman="Cayman")

####PBDB COMPARISON----
library(data.table)
pb<-fread("S3_pbdb_data_200504.csv", header = T, sep = ',',select=c("phylum","class","order","family","genus","formation")) 
pb2<-as.data.frame(pb)

taxinfo<-fread("S3_AllWebTaxonomy_updated_200903.csv",header=T,sep=",")
#Match taxinfo names with nodelist names
taxinfo$node_name <- gsub("[[:space:]]", ".", taxinfo$node_name)
taxinfo$node_name <- gsub(",", "X.", taxinfo$node_name)
taxinfo$node_name <- str_replace_all(taxinfo$node_name, "[:punct:]", ".")
taxinfo<-subset(taxinfo,node_name!="")

#Number of rank in PBDB - genera first, then other ranks
taxinfomini<-as.data.frame(subset(taxinfo,Food_web!="Messel" & Food_web!="Chengjiang" & Food_web!="Burgess" & kingdom=="Metazoa"))
pb_genera<-sapply(pb2[,c("genus")],tolower)
pb_genera<-unique(pb_genera)
pb_genera<-pb_genera[pb_genera != ""]
us_genera<-sapply(taxinfomini[,c("name")],tolower)
us_genera<-unique(us_genera)
us_genera<-us_genera[us_genera != ""]
us_genera <- gsub("([A-Za-z]+).*", "\\\\1",us_genera)

print(length(intersect(pb_genera,us_genera)))
print(length(us_genera))

ranks<-c("class","order","family")

for(ralev in ranks) {
  pb_genera<-sapply(pb2[,c(ralev)],tolower)
  pb_genera<-unique(pb_genera)
  pb_genera<-pb_genera[pb_genera != ""]
  
  us_genera<-sapply(taxinfomini[,c(ralev)],tolower)
  us_genera<-unique(us_genera)
  us_genera<-us_genera[us_genera != ""]

  print(ralev)
  print(length(intersect(pb_genera,us_genera)))
  print(length(us_genera))
    
}

#Make shelly PBDB data
shelly_dat<-c()  
shelly_comb_PBDB<-list()
ranks<-c("class","order","family")
run_as_trophic_pbdb<-TRUE

for(ralev in ranks) {
  
  #Change level where needed
  pb_genera<-sapply(pb2[,c(ralev)],tolower)
  pb_genera<-unique(pb_genera)
  pb_genera<-pb_genera[pb_genera != ""]
  
  for(foss_web in names(allwebs_comb)) {
    
    selected_loc<-allwebs_comb[[foss_web]]
    
    taxinfosub<-subset(taxinfo,Food_web==foss_web)
    
    fw_data<-selected_loc[["FW"]]
    fw_nodes<-selected_loc[["NODE"]]
    
    OGS<-nrow(fw_nodes)
    
    fw_nodes2<-fw_nodes
    
    fw_nodes2b<-merge(fw_nodes2,taxinfosub,by.x="Species",by.y="node_name",all.x=TRUE)
    fw_nodes2b$genus<-sapply(gsub("([A-Za-z]+).*", "\\\\1",fw_nodes2b[,c(ralev)]),tolower) 
    
    #Change level if needed
    fw_nodes3<-fw_nodes2b[sapply(fw_nodes2b$genus,tolower) %in% pb_genera,]
    
    row_nums<-as.numeric(as.character(rownames(fw_nodes3)))
    
    fw_nodes_shelly<-fw_nodes3[,c("X","Species","Class")]
    fw_data_shelly<-fw_data[row_nums,row_nums]
    
    web_comb<-list(FW=fw_data_shelly,NODE=fw_nodes_shelly)
    
    shelly_comb_PBDB[[foss_web]] <- web_comb
    
    lag_dat<-(as.matrix(shelly_comb_PBDB[[foss_web]][["FW"]]))
    
    if(run_as_trophic_pbdb){
      lag_dat<-make.trophic.species.web(lag_dat)
    }else{}
    
    
    swtl<-calc.mean.swtl(lag_dat)
    SOI<-calc.SOI(lag_dat)    
    diam<-calc.diam(lag_dat)
    C<-calc.C(lag_dat)
    S<-calc.S(lag_dat)
    degree<-calc.norm.mean.degree(lag_dat)
    cc<-calc.cc(lag_dat)
    
    avpath<-calc.mean.path.length(lag_dat)
    
    
    lag_metrics <- as.data.frame(cbind(OGS,
                                       swtl,
                                       
                                       SOI,
                                       diam,
                                       C,
                                       
                                       S,
                                       
                                       degree,
                                       
                                       cc,
                                       
                                       avpath,
                                       
                                       web=foss_web,
                                       rank=ralev
    ))
    
    
    
    shelly_dat<-rbind(shelly_dat,lag_metrics)
    print(foss_web)
    
  }
  print(ralev)
}



shelly_dat2<-reshape::melt(shelly_dat,id=c("web","OGS","S","rank"))
shelly_dat2$OGS<-as.numeric(as.character(shelly_dat2$OGS))
shelly_dat2$S<-as.numeric(as.character(shelly_dat2$S))
shelly_dat2$value<-as.numeric(as.character(shelly_dat2$value))
shelly_dat2$node_loss<-as.numeric(as.character(((shelly_dat2$OGS - shelly_dat2$S) / (shelly_dat2$OGS))))
shelly_dat3<-subset(shelly_dat2, web!="Burgess" & web!="Chengjiang" & web!="Messel")


att_dat<-c()  
att_comb_PBDB<-list()
for(foss_web in names(allwebs_comb)) {
  
  selected_loc<-allwebs_comb[[foss_web]]
  
  taxinfosub<-subset(taxinfo,Food_web==foss_web)
  
  fw_data<-selected_loc[["FW"]]
  fw_nodes<-selected_loc[["NODE"]]
  
  OGS<-nrow(fw_nodes)
  
  fw_nodes2<-traits(fw_nodes,name=ifelse(foss_web=="Caiman","Cayman",foss_web))
  fw_nodes2<-subset(fw_nodes2,HardSoft!=1 & HardSoft!=2)
  
  row_nums<-as.numeric(as.character(rownames(fw_nodes2)))
  
  fw_nodes_shelly<-fw_nodes2[,c("X","Species","Class")]
  fw_data_shelly<-fw_data[row_nums,row_nums]
  
  web_comb<-list(FW=fw_data_shelly,NODE=fw_nodes_shelly)
  
  att_comb_PBDB[[foss_web]] <- web_comb
  
  lag_dat<-(as.matrix(att_comb_PBDB[[foss_web]][["FW"]]))
  
  if(run_as_trophic_pbdb){
    lag_dat<-make.trophic.species.web(lag_dat)
  }else{}
  
  swtl<-calc.mean.swtl(lag_dat)
  
  SOI<-calc.SOI(lag_dat)    
  diam<-calc.diam(lag_dat)
  C<-calc.C(lag_dat)
  
  S<-calc.S(lag_dat)
  degree<-calc.norm.mean.degree(lag_dat)
  
  cc<-calc.cc(lag_dat)
  
  avpath<-calc.mean.path.length(lag_dat)
  
  
  lag_metrics <- as.data.frame(cbind(OGS,
                                     swtl,
                                     
                                     SOI,
                                     diam,
                                     C,
                                     
                                     S,
                                     
                                     degree,
                                     
                                     cc,
                                     
                                     avpath,
                                     
                                     web=foss_web,
                                     rank="Soft-bodied and intermediate removed"
  ))
  
  
  att_dat<-rbind(att_dat,lag_metrics)
  print(foss_web)
  
}
att_dat2<-reshape::melt(att_dat,id=c("web","OGS","S","rank"))
att_dat2$OGS<-as.numeric(as.character(att_dat2$OGS))
att_dat2$S<-as.numeric(as.character(att_dat2$S))
att_dat2$value<-as.numeric(as.character(att_dat2$value))
att_dat2$node_loss<-as.numeric(as.character(((att_dat2$OGS - att_dat2$S) / (att_dat2$OGS))))
att_dat3<-subset(att_dat2, web!="Burgess" & web!="Chengjiang" & web!="Messel")




soft_dat<-c()  
soft_comb_PBDB<-list()
for(foss_web in names(allwebs_comb)) {
  
  selected_loc<-allwebs_comb[[foss_web]]
  
  taxinfosub<-subset(taxinfo,Food_web==foss_web)
  
  fw_data<-selected_loc[["FW"]]
  fw_nodes<-selected_loc[["NODE"]]
  
  OGS<-nrow(fw_nodes)
  
  fw_nodes2<-traits(fw_nodes,name=ifelse(foss_web=="Caiman","Cayman",foss_web))
  fw_nodes2<-subset(fw_nodes2,HardSoft!=1)
  
  row_nums<-as.numeric(as.character(rownames(fw_nodes2)))
  
  fw_nodes_shelly<-fw_nodes2[,c("X","Species","Class")]
  fw_data_shelly<-fw_data[row_nums,row_nums]
  
  web_comb<-list(FW=fw_data_shelly,NODE=fw_nodes_shelly)
  
  soft_comb_PBDB[[foss_web]] <- web_comb
  
  lag_dat<-(as.matrix(soft_comb_PBDB[[foss_web]][["FW"]]))
  
  if(run_as_trophic_pbdb){
    lag_dat<-make.trophic.species.web(lag_dat)
  }else{}
  
  swtl<-calc.mean.swtl(lag_dat)
  
  SOI<-calc.SOI(lag_dat)    
  diam<-calc.diam(lag_dat)
  C<-calc.C(lag_dat)
  
  S<-calc.S(lag_dat)
  degree<-calc.norm.mean.degree(lag_dat)
  
  cc<-calc.cc(lag_dat)
  
  avpath<-calc.mean.path.length(lag_dat)
  
  
  lag_metrics <- as.data.frame(cbind(OGS,
                                     swtl,
                                     
                                     SOI,
                                     diam,
                                     C,
                                     
                                     S,
                                     
                                     degree,
                                     
                                     cc,
                                     
                                     avpath,
                                     
                                     web=foss_web,
                                     rank="Soft-bodied removed"
  ))
  
  
  soft_dat<-rbind(soft_dat,lag_metrics)
  print(foss_web)
  
}
soft_dat2<-reshape::melt(soft_dat,id=c("web","OGS","S","rank"))
soft_dat2$OGS<-as.numeric(as.character(soft_dat2$OGS))
soft_dat2$S<-as.numeric(as.character(soft_dat2$S))
soft_dat2$value<-as.numeric(as.character(soft_dat2$value))
soft_dat2$node_loss<-as.numeric(as.character(((soft_dat2$OGS - soft_dat2$S) / (soft_dat2$OGS))))
soft_dat3<-subset(soft_dat2, web!="Burgess" & web!="Chengjiang" & web!="Messel")


join_att_pbdb<-rbind(cbind(soft_dat3,type="Pres. group method"),cbind(att_dat3,type="Pres. group method"),cbind(shelly_dat3,type="PBDB method"))

join_att_pbdb2<-merge(join_att_pbdb,OG_full_net_stats2,by=c("web","variable"),all.x=TRUE)
join_att_pbdb2$web<-factor(join_att_pbdb2$web, levels = size.order)
join_att_pbdb2<-merge(join_att_pbdb2,NetAtt_names,by.x="variable",by.y=1,all.x=TRUE)

join_att_pbdb2$NewNetAtt<-factor(join_att_pbdb2$NewNetAtt, levels = c("Mean TL","SOI","Diameter","Degree","Clustering","CPL","Connectance"))

join_att_pbdb2$dif<-join_att_pbdb2$value.y-join_att_pbdb2$value.x
join_att_pbdb2$perc<-(join_att_pbdb2$value.y-join_att_pbdb2$value.x)/join_att_pbdb2$value.y

join_att_pbdb_size3<-join_att_pbdb2 %>% group_by(web,NewNetAtt) %>% mutate(meanNL=mean(node_loss,na.rm=T),meanVAL=mean(value.x,na.rm=T))
join_att_pbdb_size3<-join_att_pbdb_size3[,c("web","NewNetAtt","meanNL","meanVAL","value.y")]
join_att_pbdb_size3<-unique(join_att_pbdb_size3[,])

subset(join_att_pbdb_size3, NewNetAtt!="Diameter" & NewNetAtt!="Degree")


####GENERATE BASIC NICHE MODELS-REQUIRES PARAMETERS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----

fw.list<-list(  Burgess_fw,
                Burgess_fw
                ,Caiman_fw
                ,Chengjiang_fw
                ,Chilean_fw
                ,LR_fw
                ,Messel_fw
                ,Sanak_fw
                ,Weddell_fw)

all_niche_stats<-c()
nicheruns<-1000


for(i in c(TRUE,FALSE)) {
  
  fw.list2<-fw.list
  
  #Make trophic species versions
  if(i==TRUE){
    fw.list2<-lapply(fw.list2,make.trophic.species.web)
  }else{
    
  }
  
  t<-calc.metrics(fw.list2,8)
  t2<-subset(t,type=="Fossilized")
  t2$web<-c("Burgess","Caiman","Chengjiang","Chilean","LR","Messel","Sanak","Weddell")
  
  t2$type<-NULL
  t2$run<-0
  
  t3<-reshape::melt(t2,id=c("web","run"))
  t3<-subset(t3,variable!="C" & variable!="S")
  
  t4<-t3
  t4$foss<-ifelse(is.element(t4$web,names(fossil_comb)),"fossil_web","modern_web")
  
  
  library(ggplot2)
  ggplot(data=t4,aes(x=web,y=value,color=foss))+geom_point()+facet_wrap(.~variable,scales="free")
  
  
  library(foreach)
  library(doParallel)
  library(R.utils)
  
  
  print(Sys.time())
  nicmets<-c()
  myCluster <- makeCluster(7)
  registerDoParallel(myCluster)
  nicmets<-foreach(env=1:nicheruns,.combine=rbind,.packages=c("dplyr","data.table","igraph","R.utils","stringr")) %dopar% {
  
    
    nic_web<-list()
    nic_web[[1]]<-matrix(ncol=0,nrow=0)
    for(nicrow in 1:nrow(t2)) {
      
      nic1<-t2[nicrow,]
      
      nam<-nic1[,c("web")]
      S<-as.numeric(as.character(nic1[,c("S")]))
      C<-as.numeric(as.character(nic1[,c("C")]))
      
      t <- withTimeout(niche2(S,C), timeout = 5, onTimeout = "error")
      
      if(class(t)=="matrix"){
        nic_web[[nicrow+1]]<-t
        rownames(nic_web[[nicrow+1]])<-1:S
        colnames(nic_web[[nicrow+1]])<-1:S
      }else{
        nic_web[[nicrow+1]]<-matrix(ncol=0,nrow=0)
        #rownames(nic_web[[nicrow+1]])<-1:S
        #colnames(nic_web[[nicrow+1]])<-1:S
      }
      
      
    }
    
    n<-calc.metrics(nic_web,8)
    n2<-subset(n,type=="Fossilized")
    
    n2$run<-NULL
    
    niche_mets<-(cbind(n2,run=env,web=t2[,c("web")]))
    
    niche_mets
    
    
  }
  
  stopCluster(myCluster)
  closeAllConnections()
  registerDoSEQ()
  print(Sys.time())
  
  #Plot
  nicmets2<-as.data.frame(nicmets)
  
  n3<-as.matrix(nicmets2)
  n3<-as.data.frame(n3)
  n3$web<-as.factor(as.character(n3$web))
  n3$run<-as.factor(as.character(n3$run))
  
  n3[,1]<-as.numeric(as.character(n3[,1]))
  n3[,2]<-as.numeric(as.character(n3[,2]))
  n3[,3]<-as.numeric(as.character(n3[,3]))
  n3[,4]<-as.numeric(as.character(n3[,4]))
  n3[,5]<-as.numeric(as.character(n3[,5]))
  n3[,6]<-as.numeric(as.character(n3[,6]))
  n3[,7]<-as.numeric(as.character(n3[,7]))
  n3[,8]<-as.numeric(as.character(n3[,8]))
  
  n3$type<-NULL
  
  n4<-reshape::melt((n3),id=c("web","run"))
  n4<-subset(n4,variable!="C" & variable!="S")
  
  j<-rbind(cbind(t3,type="real"),cbind(n4,type="niche"))
  j$foss<-ifelse(is.element(j$web,names(fossil_comb)),"fossil_web","modern_web")
  
  ggplot(data=j,aes(x=web,y=value,color=type,shape=type))+geom_point()+facet_wrap(.~variable,scales="free")
  
  n5<-n4 %>%
    group_by(web,variable) %>%
    summarise(avg = mean(value,na.rm=TRUE),
              SD = sd(value, na.rm=TRUE),
              min = quantile(value, probs = 0.025,na.rm=TRUE), max = quantile(value, probs = 0.975,na.rm=TRUE))
  
  
  
  n6<-merge(n5,t4,by=c("web","variable"),all.y=TRUE)
  
  n6$distinct<-ifelse((n6$value<n6$min | n6$value>n6$max),TRUE,FALSE)
  n6$web<-factor(n6$web, levels = size.order)
  
  n6<-merge(n6,NetAtt_names,by.x="variable",by.y="OldNetAtt",all.x=TRUE)
  n6$web<-factor(n6$web, levels = size.order)
  n6$NewNetAtt<-factor(n6$NewNetAtt, levels = as.character(NetAtt_names[,2]))
  
  all_niche_stats<-rbind(all_niche_stats,cbind(n6,trophic_species=i))
  
}

all_niche_stats<-as.data.frame(all_niche_stats)

all_niche_stats$label<-factor(ifelse(all_niche_stats$trophic_species,"Tr. sp. web","Unaltered web"),levels=c("Unaltered web","Tr. sp. web"))
all_niche_stats$web<-factor(all_niche_stats$web, levels = web.agesize.order)
all_niche_stats$NewNetAtt <- factor(all_niche_stats$NewNetAtt, levels = as.character(NetAtt_names[, 2]))

time<-gsub(" ","",gsub('[[:punct:] ]+',' ',Sys.time()))
nam<-paste("Sep7th_",time,"_",".csv",sep="")
write.csv(all_niche_stats,file=nam)

####CUMULATIVE DEGREE DISTRIBUTIONS----

cm_deg_dists<-c()

for(tr in c(TRUE,FALSE)){
  for(i in names(allwebs_comb)) {
    for(j in c("all","in","out")) {
      selected_loc<-allwebs_comb[[i]]  
      fw_data<-selected_loc[["FW"]]
      
      if(tr==TRUE){
        fw_data<-make.trophic.species.web(fw_data)
      }else{
        
      }
      
      
      graph<-graph_from_adjacency_matrix(as.matrix(fw_data))
      
      
      
      avlinks<-(2*ecount(graph))/vcount(graph)
      cm<-as.data.frame(degree_distribution(graph,cumulative = TRUE,mode=j))
      colnames(cm)<-c("Dist")
      
      cm2<-cbind(cm,NumOfLinks=1:nrow(cm),type=j,web=i,trophic_species=tr)
      cm2$NormNOL<-cm2$NumOfLinks/avlinks
      
      cm_deg_dists<-rbind(cm_deg_dists,cm2)
      print(i)
    }
  }
}

cm_deg_dists$type<-gsub("all","All",cm_deg_dists$type)
cm_deg_dists$type<-gsub("in","Generality",cm_deg_dists$type)
cm_deg_dists$type<-gsub("out","Vulnerability",cm_deg_dists$type)

cm_deg_dists$web<-factor(cm_deg_dists$web, levels = web.agesize.order)
cm_deg_dists$label<-factor(ifelse(cm_deg_dists$trophic_species,"Tr. sp. web","Unaltered web"),levels=c("Unaltered web","Tr. sp. web"))
cm_deg_dists$foss<-as.factor(ifelse(is.element(cm_deg_dists$web,names(fossil_comb)),"fossil_web","modern_web"))

cm_deg_dists$web<-recode(cm_deg_dists$web,Caiman="Cayman")
cm_deg_dists_fossdat<-subset(cm_deg_dists,web=="Burgess" | web=="Chengjiang" | web=="Messel")
cm_deg_dists_moddat<-subset(cm_deg_dists,web!="Burgess" | web!="Chengjiang" | web!="Messel")


####CALCULATE ALPHA VERSUS NODE LOSS-REQUIRES PARAMETERS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----
library(foreach)
library(doParallel)
library(R.utils)

#Loop different environment types

fossil.s.nodenumbers.joined<-c()

alpharuns<-5

#For niche time-out


alpha_params<-data.frame(
  stringsAsFactors = FALSE,
  mod = c("H90 I50 S10", "H25 I50 S75", "Std: H75 I50 S25", "H45 I30 S15","H50 I50 S50","H75 I25 S25", "Std: H75 I50 S25, VARIED" ),
  phard = c(0.9, 0.25, 0.75, 0.45,0.5,0.75,0.75),
  pint = c(0.5, 0.5, 0.5, 0.3,0.5,0.25,0.5),
  psoft = c(0.1, 0.75, 0.25, 0.15,0.5,0.25,0.25),
  pvar=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE) #varying param
)

fossil.s.nodenumbers.joined.mod<-c()
for(ap in 1:nrow(alpha_params)) {
  
  apmod<-as.character(alpha_params[ap,"mod"])
  
  print(Sys.time())
  myCluster <- makeCluster(7)
  registerDoParallel(myCluster)
  fossil.s.nodenumbers.joined<-foreach(env=1:99,.combine=rbind,.packages=c("dplyr","data.table","igraph","R.utils","stringr"),.errorhandling = "remove") %dopar% {
    
    #Set environmental attributes to loop
    minienv<-env/10
    e1<-as.numeric(minienv)
    e2<-as.numeric(10-minienv)
    edis<-as.character("beta")
    
    fossil.s.nodenumbers.joined2<-c()
    
    #Loop different foodwebs
    for(all in names(modern_comb)) {
      
      selected_loc<-modern_comb[[all]]  
      
      fw_data<-selected_loc[["FW"]]
      fw_nodes<-selected_loc[["NODE"]]
      
      s.nodes <- traits(fw_nodes,name=ifelse(all=="Caiman","Cayman",all))
      s.nodes<-assign.p(nodes=s.nodes, p.eql = 0.5, p.hard = alpha_params[ap,"phard"], p.med =alpha_params[ap,"pint"], p.soft = alpha_params[ap,"psoft"], vary=alpha_params[ap,"pvar"])
      
      
      
      att<-"p.hs"
      
      #Real atts
      fossil.s.nodes.loop  <- fossilize(s.nodes, alpharuns, trmnt = att, E.shape1 = e1, E.shape2 = e2, E.distr = edis, null=FALSE)
      
      fossil.s.nodenumbers.loop<-calc.att.numbers(fossil.s.nodes.loop,s.nodes,alpharuns)
      
      
      fossil.s.nodenumbers.joined2<-rbind(fossil.s.nodenumbers.joined2,
                                          
                                          cbind(fossil.s.nodenumbers.loop,model=paste(att,".real",sep=""),web=all,edistr=env,model_group=att,shuffled_atts="no"))
    }
    
    
    
    fossil.s.nodenumbers.joined2
    toprint<-as.data.frame(fossil.s.nodenumbers.joined2)
    toprint
    
    #print(env)
    
  }
  
  stopCluster(myCluster)
  closeAllConnections()
  registerDoSEQ()
  print(Sys.time())
  
  fossil.s.nodenumbers.joined.mod<-rbind(fossil.s.nodenumbers.joined.mod,cbind(fossil.s.nodenumbers.joined,mod=apmod))
  
}


fossil.s.nodenumbers.joined.mstr<-as.data.frame(fossil.s.nodenumbers.joined.mod)
fossil.s.nodenumbers.joined.mstr$Freq<-as.numeric(as.character(fossil.s.nodenumbers.joined.mstr$Freq))
fossil.s.nodenumbers.joined.mstr$edistr<-as.numeric(as.character(fossil.s.nodenumbers.joined.mstr$edistr))

fossil.s.nodenumbers.joined.mstr2<-fossil.s.nodenumbers.joined.mstr %>% 
  group_by(att,run,model,web,edistr,mod) %>% 
  summarise(Frequency = sum(Freq))

fossil.s.nodenumbers.joined.mstr2$OG<-ifelse(fossil.s.nodenumbers.joined.mstr2$run==0,"OG","Fossilized")

node_number_analysis<-fossil.s.nodenumbers.joined.mstr2

node_number_analysis_original<-subset(node_number_analysis,OG=="OG")
node_number_analysis_original<-node_number_analysis_original[,c("Frequency","web")]
node_number_analysis_original<-unique(node_number_analysis_original[,])
node_number_analysis2<-merge(node_number_analysis,node_number_analysis_original,by=c("web"))
node_number_analysis2<-subset(node_number_analysis2,OG=="Fossilized")
node_number_analysis2$DifferencePercent<-(as.numeric(as.character(node_number_analysis2$Frequency.y))-as.numeric(as.character(node_number_analysis2$Frequency.x)))/as.numeric(as.character(node_number_analysis2$Frequency.y))
node_number_analysis2$web<-factor(node_number_analysis2$web,levels=web.agesize.order)

####SELECTIVE VERSUS RANDOM FOSSILIZATION-REQUIRES PARAMETERS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----

library(foreach)
library(doParallel)
library(R.utils)

modwebs <- c(
  "Caiman",
  "Chilean",
  "LR",
  "Sanak",
  "Weddell"
)

modelruns <- 200
#envlist <- 1:99
envlist<-seq(from=1,to=99,by=1)
run_as_trophic <- TRUE
pvar<-FALSE


print(Sys.time())
argen <- c()
closeAllConnections()
myCluster <- makeCluster(6)
registerDoParallel(myCluster)
fossil.s.metrics.joined <- c()
argen <- foreach(env = envlist, .combine = rbind, .packages = c("dplyr", "data.table", "igraph", "R.utils", "stringr"), .errorhandling = "remove"
) %dopar% {
  
  
  # Set environmental attributes to loop
  minienv <- env / 10
  e1 <- as.numeric(minienv)
  e2 <- as.numeric(10 - minienv)
  edis <- as.character("beta")
  
  fossil.s.metrics.joined <- c()
  
  # Loop different foodwebs
  for (all in modwebs) {
    selected_loc <- modern_comb[[all]]
    
    fw_data <- selected_loc[["FW"]]
    fw_nodes <- selected_loc[["NODE"]]
    
    s.nodes <- traits(fw_nodes,name=ifelse(all=="Caiman","Cayman",all))
    s.nodes <- assign.p(nodes = s.nodes, p.eql = 0.5, p.hard = 0.75, p.med = 0.5, p.soft = 0.25, vary=pvar)
    
    att <- "p.hs"
    
    # Real atts
    fossil.s.nodes.loop <- fossilize(s.nodes, modelruns, trmnt = att, E.shape1 = e1, E.shape2 = e2, E.distr = edis, null = FALSE)
    fossil.s.webs.loop <- recreate.web(fossil.s.nodes.loop, fw_data, modelruns)
    if (run_as_trophic == TRUE) {
      fossil.s.webs.loop <- lapply(fossil.s.webs.loop, make.trophic.species.web)
    } else {
    }
    fossil.s.metrics.loop <- calc.metrics(fossil.s.webs.loop, modelruns)
    
    # Real atts niche model
    nic_web <- list()
    nic_web[[1]] <- matrix(ncol = 0, nrow = 0)
    nic <- subset(fossil.s.metrics.loop, type == "Fossilized")
    for (nicrow in 1:nrow(nic)) {
      nic1 <- nic[nicrow, ]
      
      S <- as.numeric(as.character(nic1[, c("S")]))
      C <- as.numeric(as.character(nic1[, c("C")]))
      
      t <- withTimeout(niche2(S, C), timeout = 5, onTimeout = "error")
      
      if (class(t) == "matrix") {
        nic_web[[nicrow + 1]] <- t
        rownames(nic_web[[nicrow + 1]]) <- 1:S
        colnames(nic_web[[nicrow + 1]]) <- 1:S
      } else {
        nic_web[[nicrow + 1]] <- matrix(ncol = 0, nrow = 0)

      }
    }
    if (run_as_trophic == TRUE) {
      nic_web <- lapply(nic_web, make.trophic.species.web)
    } else {
    }
    fossil.s.metrics.loop.niche <- calc.metrics(nic_web, nrow(nic))
    
    
    
    # Shuffled atts
    fossil.s.nodes.null.shuffle.loop <- fossilize(s.nodes, modelruns, trmnt = att, E.shape1 = e1, E.shape2 = e2, E.distr = edis, null = "null.shuffle")
    fossil.s.webs.null.shuffle.loop <- recreate.web(fossil.s.nodes.null.shuffle.loop, fw_data, modelruns)
    if (run_as_trophic == TRUE) {
      fossil.s.webs.null.shuffle.loop <- lapply(fossil.s.webs.null.shuffle.loop, make.trophic.species.web)
    } else {
    }
    fossil.s.metrics.null.shuffle.loop <- calc.metrics(fossil.s.webs.null.shuffle.loop, modelruns)
     
    # Shuffled atts niche model
    nic_web <- list()
    nic_web[[1]] <- matrix(ncol = 0, nrow = 0)
    nic <- subset(fossil.s.metrics.null.shuffle.loop, type == "Fossilized")
    for (nicrow in 1:nrow(nic)) {
      nic1 <- nic[nicrow, ]
      
      S <- as.numeric(as.character(nic1[, c("S")]))
      C <- as.numeric(as.character(nic1[, c("C")]))
      
      t <- withTimeout(niche2(S, C), timeout = 5, onTimeout = "error")
      
      if (class(t) == "matrix") {
        nic_web[[nicrow + 1]] <- t
        rownames(nic_web[[nicrow + 1]]) <- 1:S
        colnames(nic_web[[nicrow + 1]]) <- 1:S
      } else {
        nic_web[[nicrow + 1]] <- matrix(ncol = 0, nrow = 0)
      
      }
    }
    if (run_as_trophic == TRUE) {
      nic_web <- lapply(nic_web, make.trophic.species.web)
    } else {
    }
    fossil.s.metrics.null.shuffle.loop.niche <- calc.metrics(nic_web, nrow(nic))
    
    fossil.s.metrics.joined <- rbind(
      fossil.s.metrics.joined,
      cbind(fossil.s.metrics.null.shuffle.loop.niche, model = paste(att, ".null.shuffle.niche", sep = ""), web = all, edistr = env, model_group = att, shuffled_atts = "null.shuffle", niche = "Yes"),
      cbind(fossil.s.metrics.null.shuffle.loop, model = paste(att, ".null.shuffle", sep = ""), web = all, edistr = env, model_group = att, shuffled_atts = "null.shuffle", niche = "No"),
      cbind(fossil.s.metrics.loop.niche, model = paste(att, ".real.niche", sep = ""), web = all, edistr = env, model_group = att, shuffled_atts = "no", niche = "Yes"),
      cbind(fossil.s.metrics.loop, model = paste(att, ".real", sep = ""), web = all, edistr = env, model_group = att, shuffled_atts = "no", niche = "No")
    )
  }
  
  time<-gsub(" ","",gsub('[[:punct:] ]+',' ',Sys.time()))
  nam<-paste("Sep4th_",time,"_",env,".csv",sep="")
  write.csv(fossil.s.metrics.joined,file=nam)
  
  fossil.s.metrics.joined
  
}

stopCluster(myCluster)
closeAllConnections()
registerDoSEQ()
print(Sys.time())

# Plot
dat <- as.data.frame(argen)

dat$type <- as.factor(dat$type)
cols.num <- c(1:8)
dat[cols.num] <- (lapply(dat[cols.num], as.character))
dat[cols.num] <- (lapply(dat[cols.num], as.numeric))

fossil.metrics.withniche <- reshape::melt(dat, id = c("run", "type", "model", "web", "edistr", "model_group", "shuffled_atts", "niche"))


# Drop metrics we aren't using
fossil_model_values <- fossil.metrics.withniche

fossil_model_values$web <- factor(fossil_model_values$web, levels = size.order)
fossil_model_values <- fossil_model_values %>%
  mutate(model = recode(model,
                        "p.hs.null.shuffle.niche" = "Random-based niche models",
                        "p.hs.null.shuffle" = "Random fossilization",
                        "p.hs.real" = "Selective fossilization",
                        "p.hs.real.niche" = "Selective-based niche models"
  ))
fossil_model_values$model <- factor(fossil_model_values$model, levels = c(
  "Random fossilization",
  c("Random-based niche models"),
  c("Selective fossilization"),
  c("Selective-based niche models")
))


fossil_model_values <- merge(fossil_model_values, NetAtt_names, by.x = "variable", by.y = 1, all.x = TRUE)
fossil_model_values$variable <- factor(fossil_model_values$NewNetAtt, levels = as.character(NetAtt_names[, 2]))


# character variables
cols.num <- c("edistr", "value")
fossil_model_values[cols.num] <- (lapply(fossil_model_values[cols.num], as.character))
fossil_model_values[cols.num] <- (lapply(fossil_model_values[cols.num], as.numeric))



# Calculate node loss bins for niche and fossil models
# Filter non-numeric values
fossil_model_values <- fossil_model_values %>% filter_at(vars(value), all_vars(!is.infinite(.)))
OGsize <- subset(fossil_model_values, type == "Original" & variable == "S")
OGsize <- OGsize[, c("web", "value")]
OGsize <- unique(OGsize[, ])
OGsize <- subset(OGsize, value > 0)
colnames(OGsize) <- c("web", "OG_size")
NEWsize <- subset(fossil_model_values, type == "Fossilized" & variable == "S")
NEWsize <- NEWsize[, c("web", "value", "run", "model", "edistr", "model_group", "shuffled_atts", "niche")]
NEWsize <- unique(NEWsize[, ])
colnames(NEWsize) <- c("web", "NEW_size", "run", "model", "edistr", "model_group", "shuffled_atts", "niche")
JOINEDsize <- merge(OGsize, NEWsize, by = "web", all.y = TRUE)
JOINEDsize$node_loss <- as.numeric(as.character(((JOINEDsize$OG_size - JOINEDsize$NEW_size) / (JOINEDsize$OG_size))))
fossil_model_values <- merge(fossil_model_values, JOINEDsize, by = c("web", "run", "model", "edistr", "model_group", "shuffled_atts", "niche"), all.x = TRUE)
fossil_model_values <- subset(fossil_model_values, type == "Fossilized")

#Drop runs with less than 3 nodes
fossil_model_values<-subset(fossil_model_values,NEW_size>=3)

# Drop metrics we aren't using
fossil_model_values2 <- subset(fossil_model_values, variable == "Mean TL" | variable == "SOI" | variable == "Diameter" | variable == "Connectance" | variable == "Degree" | variable == "Clustering" | variable == "CPL")


# Assign to node loss bins
library(Hmisc)
bin.size <- rep(seq(from = 0, to = 1, by = 0.05), 1)
fossil_model_values2$node_loss_bin <- Hmisc::cut2(fossil_model_values2$node_loss, bin.size)
# Drop after comma and brackets
fossil_model_values2$node_loss_bin1 <- gsub("^(.*?),.*", "\\\\1", fossil_model_values2$node_loss_bin)
fossil_model_values2$node_loss_bin1 <- gsub("[()]", "", fossil_model_values2$node_loss_bin1)
fossil_model_values2$node_loss_bin1 <- as.numeric(as.character(gsub("\\\\[|\\\\]", "", fossil_model_values2$node_loss_bin1)))
# Drop before comma and brackets
fossil_model_values2$node_loss_bin2 <- sub(".*,\\\\s*", "\\\\1", fossil_model_values2$node_loss_bin)
fossil_model_values2$node_loss_bin2 <- gsub("[()]", "", fossil_model_values2$node_loss_bin2)
fossil_model_values2$node_loss_bin2 <- as.numeric(as.character(gsub("\\\\[|\\\\]", "", fossil_model_values2$node_loss_bin2)))
fossil_model_values2$node_loss_bin <- (((fossil_model_values2$node_loss_bin1) + (fossil_model_values2$node_loss_bin2)) / 2)
# Drop extra rows
fossil_model_values2$node_loss_bin1 <- NULL
fossil_model_values2$node_loss_bin2 <- NULL

detach(package:Hmisc)
detach(package:plyr)


raw_ribbon_plot <- fossil_model_values2 %>%
  group_by(web, variable, type, shuffled_atts, node_loss_bin, model) %>%
  # Sometimes this function wants summarize or summarise
  summarise(
    ME_avg = mean(value, na.rm = TRUE),
    ME_SD = sd(value, na.rm = TRUE),
    ME_sample_size = n(),
    ME_med = median(value, na.rm = TRUE),
    ME_min = quantile(value, probs = 0.025, na.rm = TRUE), ME_max = quantile(value, probs = 0.975, na.rm = TRUE)
  )


# 1. Generate median and 95% CIs for all webs
library(dplyr)
niche_model_distributions <- subset(fossil_model_values2, niche == "Yes") %>%
  group_by(web, variable, type, shuffled_atts, model_group, node_loss_bin) %>%
  # Sometimes this function wants summarize or summarise
  summarise(
    niche_avg = mean(value, na.rm = TRUE),
    niche_SD = sd(value, na.rm = TRUE),
    niche_sample_size = n(),
    niche_med = median(value, na.rm = TRUE),
    niche_min = quantile(value, probs = 0.025, na.rm = TRUE), niche_max = quantile(value, probs = 0.975, na.rm = TRUE)
  )

### Compare niche models and respective random/selective info loss CIs

# Seperate non-niche webs
non_niche_model_values <- subset(fossil_model_values2, niche == "No")
# Drop NA rows
non_niche_model_values <- non_niche_model_values[complete.cases(non_niche_model_values$value), ]

# Match real+null webs to niche model webs
niche_model_comparisons <- merge(non_niche_model_values, niche_model_distributions,
                                 by = c("web", "variable", "type", "shuffled_atts", "model_group", "node_loss_bin"), all.x = TRUE
)
niche_model_comparisons <- niche_model_comparisons[, c(
  "web", "variable", "type", "shuffled_atts", "run", "model", "edistr",
  "value", "node_loss", "node_loss_bin", "niche_med", "niche_min", "niche_max"
)]


# A. Calculate model error

# Dif between niche med and CIs
niche_model_comparisons$niche_upper_dif <- niche_model_comparisons$niche_max - niche_model_comparisons$niche_med
niche_model_comparisons$niche_lower_dif <- niche_model_comparisons$niche_med - niche_model_comparisons$niche_min

niche_model_comparisons$ME <- ifelse(niche_model_comparisons$value > niche_model_comparisons$niche_med,
                                     ((niche_model_comparisons$value - niche_model_comparisons$niche_med) / niche_model_comparisons$niche_upper_dif),
                                     ((niche_model_comparisons$value - niche_model_comparisons$niche_med) / niche_model_comparisons$niche_lower_dif)
)
niche_model_comparisons <- niche_model_comparisons %>% filter_at(vars(ME), all_vars(!is.infinite(.)))
niche_model_comparisons <- niche_model_comparisons[complete.cases(niche_model_comparisons), ]


# B. Calculate 95% CI of MEs
niche_model_comparisons2 <- niche_model_comparisons %>%
  group_by(web, variable, type, shuffled_atts, node_loss_bin) %>%
  # Sometimes this function wants summarize or summarise
  summarise(
    ME_avg = mean(ME, na.rm = TRUE),
    ME_SD = sd(ME, na.rm = TRUE),
    ME_sample_size = n(),
    ME_med = median(ME, na.rm = TRUE),
    ME_min = quantile(ME, probs = 0.025, na.rm = TRUE), ME_max = quantile(ME, probs = 0.975, na.rm = TRUE)
  )

niche_model_comparisons2 <- niche_model_comparisons2 %>%
  ungroup(shuffled_atts) %>%
  mutate(shuffled_atts = recode(as.character(shuffled_atts),
                                "null.shuffle" = "Random fossilization",
                                "no" = "Selective fossilization"
  ))


#### RANDOM V SELECTIVE NODE LOSS


# Seperate non-niche webs
non_niche_model_values <- subset(fossil_model_values2, niche == "No")
# Drop NA rows
non_niche_model_values <- non_niche_model_values[complete.cases(non_niche_model_values$value), ]
# Seperate real/null
real_model_values <- subset(non_niche_model_values, shuffled_atts == "no")
null_model_values <- subset(non_niche_model_values, shuffled_atts == "null.shuffle")

rs_model_distributions <- null_model_values %>%
  group_by(web, variable, type, shuffled_atts, model_group, node_loss_bin) %>%
  # Sometimes this function wants summarize or summarise
  summarise(
    null_avg = mean(value, na.rm = TRUE),
    null_SD = sd(value, na.rm = TRUE),
    null_sample_size = n(),
    null_med = median(value, na.rm = TRUE),
    null_min = quantile(value, probs = 0.025, na.rm = TRUE), null_max = quantile(value, probs = 0.975, na.rm = TRUE)
  )



# Match real+null webs to niche model webs
rs_model_comparisons <- merge(real_model_values, rs_model_distributions,
                              by = c("web", "variable", "type", "node_loss_bin"), all.x = TRUE
)
rs_model_comparisons <- rs_model_comparisons[, c(
  "web", "variable", "type", "run", "model", "edistr",
  "value", "node_loss", "node_loss_bin", "null_med", "null_min", "null_max"
)]


# A. Calculate model error

# Dif between null med and CIs
rs_model_comparisons$null_upper_dif <- rs_model_comparisons$null_max - rs_model_comparisons$null_med
rs_model_comparisons$null_lower_dif <- rs_model_comparisons$null_med - rs_model_comparisons$null_min

rs_model_comparisons$ME <- ifelse(rs_model_comparisons$value > rs_model_comparisons$null_med,
                                  ((rs_model_comparisons$value - rs_model_comparisons$null_med) / rs_model_comparisons$null_upper_dif),
                                  ((rs_model_comparisons$value - rs_model_comparisons$null_med) / rs_model_comparisons$null_lower_dif)
)
rs_model_comparisons <- rs_model_comparisons %>% filter_at(vars(ME), all_vars(!is.infinite(.)))
rs_model_comparisons <- rs_model_comparisons[complete.cases(rs_model_comparisons), ]

# Look at mean model errors
rs_model_comparisons_means_web <- rs_model_comparisons %>%
  group_by(web, type) %>%
  # Sometimes this function wants summarize or summarise
  summarise(ME_avg2 = mean(abs(ME), na.rm = TRUE))
rs_model_comparisons_means_web$lev<-"web"
rs_model_comparisons_means_variable <- rs_model_comparisons %>%
  group_by(variable, type) %>%
  # Sometimes this function wants summarize or summarise
  summarise(ME_avg2 = mean(abs(ME), na.rm = TRUE))
rs_model_comparisons_means_variable$lev<-"variable"
rs_model_comparisons_means <- rbind(rs_model_comparisons_means_web, rs_model_comparisons_means_variable)

# B. Calculate 95% CI of MEs
rs_model_comparisons2 <- rs_model_comparisons %>%
  group_by(web, variable, type, node_loss_bin) %>%
  # Sometimes this function wants summarize or summarise
  summarise(
    ME_avg = mean(ME, na.rm = TRUE),
    ME_SD = sd(ME, na.rm = TRUE),
    ME_sample_size = n(),
    ME_med = median(ME, na.rm = TRUE),
    ME_min = quantile(ME, probs = 0.025, na.rm = TRUE), ME_max = quantile(ME, probs = 0.975, na.rm = TRUE)
  )



####ANALYSIS OF TAXONOMIC RANKS----
tara<-read.csv("S3_AnalysisOfTaxonomicRanks.csv")

tara2<-tara %>% 
  group_by(Food_web) %>%
  count(counting=Lowest_taxonomic) %>%
  mutate(perc=n/sum(n))%>%
  drop_na(counting)
tara2<-subset(tara2,counting!="")
tara2$counting<-stringr::str_to_title(tara2$counting)
tara2$counting<-paste0(toupper(substr(tara2$counting, 1, 1)), substr(tara2$counting, 2, nchar(tara2$counting)))
tara2$counting<-factor(tara2$counting,levels=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
tara2$Food_web<-factor(tara2$Food_web,levels=web.agesize.order)
