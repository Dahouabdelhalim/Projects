####Caiman web comparison

source("S3_MainFigures_Setup_200826.R")

NetAtt_names<-as.data.frame(rbind(c("trophlev","Mean Tr. Lev."),
                                  c("swtl","Mean TL"),
                                  c("oindex","Mean Omn. Ind."),
                                  c("SOI","SOI"),
                                  c("diam","Diameter"),
                                  c("C","Connectance"),
                                  c("degree","Degree"),
                                  c("cc","Clustering"),
                                  c("avpath","CPL"),
                                  c("meanfp","Mean FP"),
                                  c("S","S")
))
colnames(NetAtt_names)<-c("OldNetAtt","NewNetAtt")

NodeAtt_names<-as.data.frame(rbind(c("SWTL","TL"),
                                   c("OI","OI"),
                                   #c("deg_in","Generality"),
                                   c("deg_in_norm","Generality"),
                                   #c("deg_out","Vulnerability"),
                                   c("deg_out_norm","Vulnerability"),
                                   c("btw","Btw.")
                                   
))
colnames(NodeAtt_names)<-c("OldNodeAtt","NewNodeAtt")


#Import food webs
Cayman_fw<-read.csv("S1_coralreef_caiman_caymanislands.csv", row.names = 1, stringsAsFactors = FALSE)  %>% fixnames.fw
Jamaica_fw<-read.csv("S1_coralreef_caiman_jamaica.csv", row.names = 1, stringsAsFactors = FALSE)   %>% fixnames.fw
Cuba_fw<-read.csv("S1_coralreef_caiman_cuba.csv", row.names = 1, stringsAsFactors = FALSE)   %>% fixnames.fw


#Make node list
Caiman_nodes <- read.csv("S1_coralreef_caiman_meta_updated.csv", stringsAsFactors = FALSE) %>% mutate(Species=Nodename) %>% fixnames.nodes

Cayman_nodes<-Caiman_nodes %>% filter(Species %in% rownames(Cayman_fw))
setdiff(rownames(Cayman_fw),Cayman_nodes$Species)

Jamaica_nodes<-Caiman_nodes %>% filter(Species %in% rownames(Jamaica_fw))
setdiff(rownames(Jamaica_fw),Jamaica_nodes$Species)

Cuba_nodes<-Caiman_nodes %>% filter(Species %in% rownames(Cuba_fw))  
setdiff(rownames(Cuba_fw),Cuba_nodes$Species)


#Combine node and FW
Jamaica_comb<-list(FW=Jamaica_fw,NODE=Jamaica_nodes)
Cayman_comb<-list(FW=Cayman_fw,NODE=Cayman_nodes)
Cuba_comb<-list(FW=Cuba_fw,NODE=Cuba_nodes)

caiman_comb<-list(
  Cayman=Cayman_comb
  ,Jamaica=Jamaica_comb
  ,Cuba=Cuba_comb
)


fplist<-bind_rows(cbind(Cayman_nodes,Food_web="Cayman"),
                  cbind(Jamaica_nodes,Food_web="Jamaica"),
                  cbind(Cuba_nodes,Food_web="Cuba")
                  )
length(unique(fplist$Species))



pgs<-read.csv("S3_PresGrAssignments_updated_200902.csv") %>% filter(Food_web %in% c("Cuba","Jamaica","Cayman"))
pgs2<- pgs %>% filter(node_name %in% fplist$Nodename)
setdiff(pgs$node_name,unique(fplist$Nodename))

fpcount<-pgs2%>%
  group_by(Food_web,BodyType) %>%
  summarise(counting=n()) %>%
  mutate(BodyType = recode(BodyType, 
                           "1" = "Soft-bodied",
                           "2" = "Intermediate",
                           "3" = "Hard-bodied"))

fpcount$Food_web<-factor(fpcount$Food_web, levels = c("Jamaica","Cuba","Cayman"))
fpcount$BodyType<-ordered(fpcount$BodyType, levels = c("Soft-bodied", "Intermediate", "Hard-bodied"))



ggplot(fpcount, aes(x = Food_web, y=counting,fill=factor(BodyType)) ) + 
  scale_fill_brewer(palette="Dark2") +
  xlab("") + ylab("Number of taxa") +
  geom_bar(stat="identity")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill=guide_legend(title="Preservation group")) 


####NETWORK WIDE STATS----

caiman_OG_full_net_stats<-c()
#Loop different foodwebs
for(all in names(caiman_comb)) {
  
  selected_loc<-caiman_comb[[all]]  
  
  fw_data<-selected_loc[["FW"]]
  fw_nodes<-selected_loc[["NODE"]]
  
  s.nodes <- traits(fw_nodes,all)
  s.nodes <- assign.p(s.nodes, p.eql, p.hard, p.med, p.soft, vary=FALSE)
  mean_fp<-as.numeric(as.character(unique(s.nodes$p.mean)))
  
  graphed_list<-list()
  graphed_list[[1]]<-as.data.frame(fw_data)
  
  OG.metrics   <- calc.metrics(graphed_list,1)
  chosen_row<-subset(OG.metrics,type=="Original")
  chosen_row<-cbind(chosen_row,meanfp=mean_fp)
  caiman_OG_full_net_stats<-rbind(caiman_OG_full_net_stats,cbind(chosen_row,web=all))
  print(all)
}


caiman_OG_full_net_stats2<-reshape::melt(caiman_OG_full_net_stats,id=c("run","type","web"))

caiman_OG_full_net_stats2<-caiman_OG_full_net_stats2 %>% left_join(NetAtt_names,by=c("variable"="OldNetAtt"))

caiman_OG_full_net_stats2$web<-factor(caiman_OG_full_net_stats2$web, levels = c("Jamaica","Cuba","Cayman"))


ggplot(caiman_OG_full_net_stats2) + 
  geom_point(aes(x=web,y=value))+
  xlab("Web") +ylab("") + # expand_limits(x = 0, y = 0) +
  facet_wrap(NewNetAtt~.,scales="free")





####NODE LEVEL METRICS----

caiman_OG_web_stats<-c()

#Loop different foodwebs
for(all in names(caiman_comb)) {
  
  selected_loc<-caiman_comb[[all]]  
  
  fw_data<-selected_loc[["FW"]]
  fw_nodes<-selected_loc[["NODE"]]
  
  s.nodes <- traits(fw_nodes,name=all)
  
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
  
  caiman_OG_web_stats<-rbind(caiman_OG_web_stats,web_list)
  
  print(paste(all))
  print(Sys.time())
}

#Plot plain stats
caiman_OG_web_stats_reduced<-caiman_OG_web_stats[,c("Species","HardSoft","btw","deg_in","deg_out","deg_in_norm","deg_out_norm",
                                      "TL","OI","SWTL","web","size")]

caiman_OG_web_stats_melt<-reshape2::melt(caiman_OG_web_stats_reduced,id=c("Species","HardSoft","web"))

caiman_OG_HS<-caiman_OG_web_stats_melt[,c("Species","HardSoft","web","variable","value")]
caiman_OG_HS$SpAtt<-"HardSoft"
colnames(caiman_OG_HS)<-c("Species","level","web","variable","value","SpAtt")
caiman_OG_test2<-rbind(caiman_OG_HS)

##Compare confidence intervals
caiman_OG_web_stats_reduced<-caiman_OG_web_stats[,c("Species","HardSoft", "btw","deg_in_norm","deg_out_norm",
                                      "TL","OI","SWTL","web")]
#Need to change column names to do.call below
colnames(caiman_OG_web_stats_reduced)<-c("Species","HardSoft", "btw","deg_in_norm","deg_out_norm",
                                  "TL","OI","SWTL","web")
caiman_OG_web_stats_melt<-reshape2::melt(caiman_OG_web_stats_reduced,id=c("Species","HardSoft","web"))
caiman_OG_web_stats_melt_HS<-caiman_OG_web_stats_melt[,c("Species","HardSoft","web","variable","value")]
colnames(caiman_OG_web_stats_melt_HS)<-c("Species","Level","web","variable","value")
caiman_OG_web_stats_melt_HS$Type<-"HardSoft"
caiman_OG_web_stats_melt_joinedBPHS<-rbind(caiman_OG_web_stats_melt_HS)
caiman_OG_web_stats_melt_joinedBPHS$Type<-as.factor(caiman_OG_web_stats_melt_joinedBPHS$Type)

caiman_OG_web_stats_melt_joinedBPHS<-caiman_OG_web_stats_melt_joinedBPHS %>% 
  mutate(Level = recode(Level, 
                           "1" = "Soft-bodied",
                           "2" = "Intermediate",
                           "3" = "Hard-bodied")) %>%
  filter(variable!="TL") %>%
  left_join(NodeAtt_names,by=c("variable"="OldNodeAtt"))

caiman_OG_web_stats_melt_joinedBPHS$web<-factor(caiman_OG_web_stats_melt_joinedBPHS$web, levels = c("Jamaica","Cuba","Cayman"))
caiman_OG_web_stats_melt_joinedBPHS$NewNodeAtt<-factor(caiman_OG_web_stats_melt_joinedBPHS$NewNodeAtt, levels = as.character(NodeAtt_names$NewNodeAtt))
caiman_OG_web_stats_melt_joinedBPHS$Level<-ordered(caiman_OG_web_stats_melt_joinedBPHS$Level, levels = c("Soft-bodied", "Intermediate", "Hard-bodied"))


ggplot(caiman_OG_web_stats_melt_joinedBPHS, aes(x=web, y=value,color=as.factor(Level),group=as.factor(web))) + 
  geom_boxplot(position=position_dodge(1),width=0.8,show.legend=F) +
  scale_x_discrete("Web") + ylab("") +
  labs(color = "Preservation group") + scale_color_brewer(palette="Dark2")+
  facet_grid(NewNodeAtt~Level,scales="free")



####CUMULATIVE DEGREE DISTRIBUTIONS----

caiman_cm_deg_dists<-c()

for(tr in c(TRUE,FALSE)){
  for(i in names(caiman_comb)) {
    for(j in c("all","in","out")) {
      selected_loc<-caiman_comb[[i]]  
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
      
      caiman_cm_deg_dists<-rbind(caiman_cm_deg_dists,cm2)
      print(i)
    }
  }
}

caiman_cm_deg_dists$type<-gsub("all","All",caiman_cm_deg_dists$type)
caiman_cm_deg_dists$type<-gsub("in","Generality",caiman_cm_deg_dists$type)
caiman_cm_deg_dists$type<-gsub("out","Vulnerability",caiman_cm_deg_dists$type)
caiman_cm_deg_dists$web<-factor(caiman_cm_deg_dists$web, levels = c("Jamaica","Cuba","Cayman"))


ggplot()+
  scale_x_continuous(trans='log',breaks=c(0,0.01,0.1,1,10)) +
  scale_y_continuous(trans='log',breaks=c(0,0.01,0.1,1)) +
  geom_point(data=caiman_cm_deg_dists,alpha=0.5,aes(x=NormNOL, y=Dist,fill=web),colour="grey",shape=21) +
  xlab("Normalized number of links") + ylab("Cumulative distribution") +
  theme_bw() +
  guides(fill=guide_legend(title="Food web")) +
  facet_grid(type~.,scales="free")

