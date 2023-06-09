####FIGURES

source("S3_MainFigures_Code.R")

####FIG 1: Node level metrics----
ggplot(apcdata2_new, aes(x=NewLevel, y=value,color=as.factor(NewLevel),group=as.factor(NewLevel))) + 
  #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_boxplot() +
  scale_x_discrete("") + ylab("") +
  labs(color = "Preservation group") + scale_color_brewer(palette="Dark2")+
  facet_wrap(NewNodeAtt~web,scales="free",ncol=8)

ggplot(apcdata3_new, aes(x=NewLevel, lower=lower_hinge, upper=upper_hinge, middle=median, ymin=lower_whisker, 
               ymax=upper_whisker,color=as.factor(NewLevel),group=as.factor(NewLevel))) + 
  #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_boxplot(stat="identity",position=position_dodge(1),width=0.8,show.legend=F) +
  scale_x_discrete("") +
  labs(color = "Preservation group") + scale_color_brewer(palette="Dark2")+
  #theme(panel.spacing.y=unit(0, "lines"),panel.spacing.x=unit(0, "lines")) +
  facet_grid(NewNodeAtt~web,scales="free")


####SUP FIG 2A: Node CI comparison----
apc4<-apc3 %>% mutate(comparison=recode(comparison,
                                        "Intermediate-Soft-bodied"="Int-SB",
                                        "Hard-bodied-Soft-bodied"="HB-SB",
                                        "Hard-bodied-Intermediate"="HB-Int"))
apc4$web<-recode(apc4$web,Caiman="Cayman")
ggplot()+
  geom_errorbar(data=subset(apc4,Type=="HardSoft" & variable!="TL"),
                aes(x=comparison,ymin=lower,ymax=upper,color=SigDif)) +
  ylab("Difference in CI") +
  xlab("") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.title=element_blank()) +
  facet_grid(NewNodeAtt~web,scales="free")

####SUP FIG 1A: HARD-SOFT BREAKDOWN----
library(RColorBrewer)
ggplot(subset(master_taxon_list_counts,type=="Body type"), aes(x = all2, y=n,fill=factor(counting)) ) + 
  scale_fill_brewer(palette="Dark2") +
  xlab("") + ylab("Number of taxa") +
  geom_bar(stat="identity")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill=guide_legend(title="Preservation group")) 


####SUP FIG 1C: SP VS TROHIC SP WEB SIZE----
ggplot()+
  geom_abline(slope=1, intercept=0,linetype="dashed",color="grey") +
  geom_point(data=trophic_list,aes(x=OG_size,y=NEW_size,fill=all,color=all),size=3) +
  guides(color=guide_legend(title="Food web"),fill=guide_legend(title="Food web"),shape=guide_legend(title="Web type")) +
  scale_shape_manual(values=c(23,21)) +
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2") +
  coord_cartesian(xlim=c(0,520),ylim=c(0,520)) +
  theme_bw()+
  xlab("Unaltered web size") + ylab("Trophic species web size")

####SUP FIG 1D: TROPHIC OVERLAP----
ggplot(tm.list3) +
  xlab("") + ylab("Mean trophic overlap ") +
  scale_x_discrete() +
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  scale_shape_manual(values=c(23,21)) +
  geom_line(aes(x=grp2,y=x,color=web,group=web),show.legend=FALSE) +
  geom_point(aes(x=grp2,y=x,color=web,shape=foss,group=web,fill=web),size=3,show.legend=FALSE) +
  theme_bw() +
  facet_wrap(vars(web),scales="free",nrow=2)

####SUP FIG X: ANALYSIS OF TAXONOMIC RANKS----
ggplot(tara2, aes(x = Food_web, y=perc*100,fill=counting)) + 
  xlab("") + ylab("Percentage of taxa by lowest taxonomic rank they are assigned") +
  scale_fill_brewer(palette="Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill=guide_legend(title="Rank"))  +
  geom_bar(stat="identity",position="stack")

####SUP FIG X: PBDB VERSUS BODY TYPE----
ggplot(join_att_pbdb2 %>% filter(web!="Caiman"))+
  geom_point(aes(x=0,y=value.y),color="black",fill="black",shape=1,size=2,stroke=2,alpha=0.7)  +
  geom_point(aes(x=node_loss,y=value.x,color=rank,shape=type.x),size=4,stroke=2,alpha=0.7)  +
  geom_point(data=join_att_pbdb_size3 %>% filter(web!="Caiman"),aes(x=meanNL,y=meanVAL),size=4,shape=18,color="black",alpha=0.7) +
  scale_shape_manual(values=c(22,25)) +
  scale_x_continuous("Percentage node loss", labels = function(x) scales::percent(x, accuracy = 1),breaks=c(0,0.25,0.5,0.75)) +
  theme_bw() +
  ylab("") +
  theme(legend.title=element_blank())+
  facet_grid(NewNetAtt~web,scales="free_y")

####SUP FIG X: CUMULATIVE DEGREE DISTRIBUTIONS----
ggplot()+
  scale_x_continuous(trans='log',breaks=c(0,0.01,0.1,1,10)) +
  scale_y_continuous(trans='log',breaks=c(0,0.01,0.1,1)) +
  scale_fill_brewer(palette="Dark2") +
  geom_point(data=cm_deg_dists_moddat,alpha=0.1,aes(x=NormNOL, y=Dist,fill=web),colour="grey",shape=21) +
  geom_point(data=cm_deg_dists_fossdat,alpha=0.8,aes(x=NormNOL, y=Dist,fill=web),colour="white",shape=23,size=2) +
  xlab("Normalized number of links") + ylab("Cumulative distribution") +
  theme_bw() +
  guides(fill=guide_legend(title="Food web")) +
  facet_grid(type~label,scales="free")




####SUP FIG X: NICHE MODELS----
all_niche_stats$web<-recode(all_niche_stats$web,Caiman="Cayman")
ggplot(data=(all_niche_stats))+
  geom_errorbar(aes(x=web,ymin=min,ymax=max),alpha=1,width=0,color="black") +
  geom_point(aes(x=web,y=value,color=distinct),show.legend=F,size=2) +
  ylab("") +xlab("") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(NewNetAtt~label,scales="free")



####SUP FIG X: ENVIRONMENTAL ALPHA VERSUS PERCENT NODE LOSS----
node_number_analysis2$web<-recode(node_number_analysis2$web,Caiman="Cayman")
ggplot(data=node_number_analysis2,aes(x=edistr/10,y=DifferencePercent,color=web))+
  geom_point(alpha=0.05)+
  geom_smooth(alpha=0)+
  xlab("Environmental filter alpha value") +
  scale_y_continuous("Node loss", labels = function(x) scales::percent(x, accuracy = 1),breaks=c(0,0.25,0.5,0.75,1)) +
  guides(color=guide_legend(title="Food web")) +
  theme_bw() +
  facet_grid(mod~.)

####SELECTIVE V RANDOM FOSSILIZATION----
####Raw model values----
raw_ribbon_plot$web<-recode(raw_ribbon_plot$web,Caiman="Cayman")
ggplot(subset(raw_ribbon_plot)) +
  geom_ribbon(aes(x=node_loss_bin,ymin=ME_min,ymax=ME_max,color=model,group=model,fill=model),alpha=0.5) +
  #geom_point(aes(x=node_loss_bin,y=ME_med,color=model)) +
  theme_bw() +
  ylab("") +
  guides(fill=guide_legend(title="Model"),color=guide_legend(title="Model"))  +
  scale_x_continuous("Percentage node loss", labels = function(x) scales::percent(x, accuracy = 1),breaks=c(0,0.25,0.5,0.75)) +
  #ggtitle(modelruns)+
  facet_grid(variable~web,scales="free")

####Relative to respective niche----
niche_model_comparisons2$web<-recode(niche_model_comparisons2$web,Caiman="Cayman")
ggplot(subset(niche_model_comparisons2,variable!="Connectance" & variable!="Diameter" & variable!="Degree" & node_loss_bin<0.975)) +
  geom_hline(yintercept=0,size=2,color="black",alpha=0.5) +
  geom_hline(yintercept=1,size=1,color="black",linetype="dashed",alpha=0.5) +
  geom_hline(yintercept=-1,size=1,color="black",linetype="dashed",alpha=0.5) +
  geom_ribbon(aes(x=node_loss_bin,ymin=ME_min,ymax=ME_max,group=shuffled_atts,fill=shuffled_atts),alpha=0.3) +
  #geom_point(aes(x=node_loss_bin,y=ME_med,color=shuffled_atts)) +
  geom_line(aes(x=node_loss_bin,y=ME_med,color=shuffled_atts,group=shuffled_atts),size=2,alpha=0.8) +
  ylab("Model error") +
  #guides(fill=guide_legend(title="Model"),color=guide_legend(title="Model"))  +
  guides(fill=FALSE,color=FALSE)  +
  scale_x_continuous("Percentage node loss", labels = function(x) scales::percent(x, accuracy = 1),breaks=c(0,0.25,0.5,0.75)) +
  theme_bw() +
  #ggtitle(modelruns)+
  facet_grid(variable~web,scales="free")

####Selective versus random----
rs_model_comparisons2$web<-recode(rs_model_comparisons2$web,Caiman="Cayman")
ggplot(subset(rs_model_comparisons2,node_loss_bin>0.025 & node_loss_bin<0.975)) +
  geom_hline(yintercept=0,size=2,color="black",alpha=0.5) +
  geom_hline(yintercept=1,size=1,color="black",linetype="dashed",alpha=0.5) +
  geom_hline(yintercept=-1,size=1,color="black",linetype="dashed",alpha=0.5) +
  #geom_ribbon(aes(x=node_loss_bin,ymin=ME_min,ymax=ME_max,fill="type"),alpha=0.5) +
  geom_ribbon(aes(x=node_loss_bin,ymin=ME_min,ymax=ME_max),color="darkolivegreen3",fill="darkolivegreen3",alpha=0.3) +
  #geom_point(aes(x=node_loss_bin,y=ME_med)) +
  geom_line(aes(x=node_loss_bin,y=ME_med),color="darkolivegreen4",size=2) +
  ylab("Model error") +
  scale_x_continuous("Percentage node loss", labels = function(x) scales::percent(x, accuracy = 1),breaks=c(0,0.25,0.5,0.75)) +
  #ggtitle(modelruns)+
  facet_grid(variable~web,scales="free")

####Selective versus random model error----
ggplot(subset(rs_model_comparisons2)) +
  geom_hline(yintercept=0,size=2,color="black",alpha=0.5) +
  geom_hline(yintercept=1,size=1,color="black",linetype="dashed",alpha=0.5) +
  geom_hline(yintercept=-1,size=1,color="black",linetype="dashed",alpha=0.5) +
  geom_line(aes(x=node_loss_bin,y=ME_med,color=variable),size=2,alpha=0.8) +
  ylab("Model error") +
  guides(color=guide_legend(title="Metric"))  +
  scale_x_continuous("Percentage node loss", labels = function(x) scales::percent(x, accuracy = 1),breaks=c(0,0.25,0.5,0.75)) +
  theme_bw() +
  #ggtitle(modelruns)+
  facet_grid(.~web,scales="free")
ggplot(subset(rs_model_comparisons2,variable!="Diameter" & variable!="Degree")) +
  geom_hline(yintercept=0,size=2,color="black",alpha=0.5) +
  geom_hline(yintercept=1,size=1,color="black",linetype="dashed",alpha=0.5) +
  geom_hline(yintercept=-1,size=1,color="black",linetype="dashed",alpha=0.5) +
  geom_line(aes(x=node_loss_bin,y=ME_med,color=web),size=2,alpha=0.8) +
  ylab("Model error") +
  guides(color=guide_legend(title="Food web"))  +
  scale_x_continuous(limits=c(0,0.95),"Percentage node loss", labels = function(x) scales::percent(x, accuracy = 1),breaks=c(0,0.25,0.5,0.75)) +
  theme_bw() +
  #ggtitle(modelruns)+
  facet_grid(.~variable,scales="free")



#############################################

####TABLE 1: SUMMARY STATS OF UNALTERED WEBS----
tab_size<-master_taxon_list_counts[,c("all2","counting","percent")]
tab_size$percent<-tab_size$percent*100
colnames(tab_size)<-c("web","BodyType","percent")
tab_size2<-tab_size %>% spread(key=BodyType,value=percent)
tab_size2<-tab_size2 %>% 
  mutate_if(is.numeric, round,digits=0)

#!!!! Changes dependent on species or tr sp web - size and connectance
tab_stats<-OG_full_net_stats[,c("web","S","C")]
tab_stats$web<-recode(tab_stats$web,Caiman="Cayman")
tab_tr_size<-as.data.frame(net.tr.size)
tab_tr_size$trcon<-signif(as.numeric(as.character(tab_tr_size$trcon)),3)
tab_tr_size$web<-recode(tab_tr_size$web,Caiman="Cayman")

#!!!! Changes dependent on species or tr sp web
tab_niche<-subset(all_niche_stats,trophic_species==TRUE)
tab_niche <-  tab_niche[,c("NewNetAtt","web","value","distinct")]
tab_nicheb<-reshape2::dcast(tab_niche, web~NewNetAtt,value.var = "distinct")
tab_nichea<-reshape2::dcast(tab_niche, web~NewNetAtt,value.var = "value")
tab_nichec<-cbind(tab_nichea,tab_nicheb)
tab_nichec<-tab_nichec[ , order(names(tab_nichec))]
tab_nichec[,c("web.1")]<-NULL
tab_niche2<-tab_nichec %>% 
  mutate_if(is.numeric, signif,digits=3)

tab_basal<-perc_basal
tab_basal$perc_basal<-as.numeric(as.character(tab_basal$perc_basal))
tab_basal$web<-recode(tab_basal$web,Caiman="Cayman")
tab_basal<-tab_basal %>% 
  mutate_if(is.numeric, round,digits=0)

#!!!! Changes dependent on species or tr sp web
tab_me<-rs_model_comparisons_means_web
tab_me$web<-recode(tab_me$web,Caiman="Cayman")

tab_env<-data.frame(
  stringsAsFactors = FALSE,
  web = c("Burgess","Chengjiang",
          "Chilean",
          "Cayman","LR",
          "Messel","Sanak","Weddell"),
  Environment = c("Marine - Basinal",
                  "Marine - Offshore","Marine - Intertidal","Marine - Reef",
                  "Freshwater - Lake","Freshwater - Lake","Marine - Nearshore",
                  "Marine - Shelf")
)

tab_join<- left_join(tab_size2,tab_stats,by="web")
tab_join<- left_join(tab_join,tab_tr_size,by="web")
tab_join<- left_join(tab_join,tab_niche2,by="web")
tab_join<- left_join(tab_join,tab_basal,by="web")
tab_join<- left_join(tab_join,tab_env,by="web")
tab_join<- left_join(tab_join,tab_me,by="web")


#write.csv(tab_join,"/Users/jos23/Google Drive/PERSONAL_NETWORK-FOSSILIZATION/TrSpSummaryTable.csv")

####Stats versus size/connectance----
cs_stats<-tab_join[,c("web","trsize","trcon","Clustering","CPL","Degree","Diameter","Mean TL","SOI")]
cs_stats <- cs_stats %>% pivot_longer(-c(web,trsize,trcon),names_to="var",values_to="val")
cs_stats$trsize<-as.numeric(as.character(cs_stats$trsize))
cs_stats$trcon<-as.numeric(as.character(cs_stats$trcon))
cs_stats$val<-as.numeric(as.character(cs_stats$val))
cs_stats$var<-factor(cs_stats$var, levels = c("Mean TL","SOI","Diameter","Degree","Clustering","CPL"))


p1<-ggplot(cs_stats)+geom_point(aes(x=trsize,y=val,color=web))+facet_grid(var~.,scales="free")+labs(color="Food web")+xlab("# Trophic Species") + ylab("")+scale_color_brewer(palette="Dark2")
p2<-ggplot(cs_stats)+geom_point(aes(x=trcon,y=val,color=web))+facet_grid(var~.,scales="free")+labs(color="Food web")+xlab("Connectance") + ylab("")+scale_color_brewer(palette="Dark2")

p1+p2+plot_layout(guides='collect')



