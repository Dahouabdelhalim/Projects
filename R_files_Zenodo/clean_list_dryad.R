##########################################################################
#
#	Script for summarizing data on West Indian extinct mammal dates
#
#	Files required in same folder
#	data: fossil.csv archeo.csv list.csv short_list_edited.csv cuba_silva2007.csv
#	data: hispaniola_sam.csv puerto_rico_sam.csv
#
##########################################################################

library(dplyr)
library(reshape2)
library(ggplot2)
library(viridis)
library(scales)
library(plyr)
rm(list = ls())

#read data
all<-read.csv("list.csv")
all2<-all[,c(2,6:87)]
all.m<-melt(all2, id="Order")
all.m<-all.m[!(is.na(all.m$value) | all.m$value ==""), ]
colnames(all.m)[2:3]<-c("island","species")
write.csv(all.m,"short_list.csv", row.names = F)

#data have been cleaned up
all3<-read.csv("short_list_edited.csv")

#list of spp. names Cuba
cnb<-read.csv("cuba_silva2007.csv")
cnb_ori<-subset(all3, Order != "Chiroptera" & island == "Cuba", select = c(species, extirpated))
cnb_del<-anti_join(cnb_ori, cnb)
cnb_add<-anti_join(cnb, cnb_ori)
cnb_add$Order<-"Rodentia"
cnb_add$island<-"Cuba"

#list of spp. names Hispaniola
hnb<-read.csv("hispaniola_sam.csv")
hnb_ori<-subset(all3, Order != "Chiroptera" & island == "Hispaniola", select = c(species, extirpated))
hnb_del<-anti_join(hnb_ori, hnb)
hnb_add<-anti_join(hnb, hnb_ori)
hnb_add$island<-"Hispaniola"
hnb_add$Order<-"Rodentia"
pnb_ori<-subset(all3, Order != "Chiroptera" & island == "Puerto Rico", select = c(species, extirpated))

#list of spp. names Puerto Rico
pnb<-read.csv("puerto_rico_sam.csv")
pnb_del<-anti_join(pnb_ori, pnb)
pnb_add<-anti_join(pnb, pnb_ori)
pnb_add$island<-"Puerto Rico"
pnb_add$Order<-"Rodentia"
all3<-rbind(all3,cnb_add)
all3<-rbind(all3,hnb_add)
all3<-rbind(all3,pnb_add)
all3<-all3[!(all3$species %in% cnb_del$species),]
all3<-all3[!(all3$species %in% hnb_del$species),]

#write accepted list
write.csv(all3,"island_species_list.csv", row.names = F)

#read fossils
fos<-read.csv("fossil.csv")

#exclude all fossils much older than Holocene
fos1<-subset(fos, synonym=="N" & is.na(age)==F & Extinct =="Y" & age<25000)

#exclude species that went extinct at time of specimens collections (post 1492) 
fos1<-subset(fos1, Material!="specimen")
fos2<-fos1[,c(2, 5:8,10,14)]
fos3<-aggregate(age~island+Order+Species, data=fos2,FUN= length)
colnames(fos3)[3]<-"species"
all4<-merge(all3, fos3, all.x=T, all.y=T)
all4 <- within(all4, extirpated <- ifelse(is.na(extirpated), 1, extirpated))
all4 <- within(all4, age <- ifelse(is.na(age), 0, age))
nspp_ext<-subset(aggregate(species~Order+ island+ extirpated, data=all4,FUN= length), extirpated==1)
nspp<-aggregate(species~Order+ island, data=all4,FUN= length)
colnames(nspp)[3]<-"total"
nspp_ext$extirpated<-NULL
colnames(nspp_ext)[3]<-"extirpated"
nspp_dat<-aggregate(age~Order+ island, data=subset(all4, age>0),FUN= length)
colnames(nspp_dat)[3]<-"dated"
nspp_mt1<-aggregate(age~Order+ island+species, data=subset(all4, age>1 & extirpated ==1),FUN= length)
nspp_mt1<-nspp_mt1[with(nspp_mt1, order(island, Order)), ]

#write out fossils with mor ethan one date
write.csv(nspp_mt1[,1:3], "available_more_than_one.csv", row.names = F)
nspp_gri<-aggregate(age~Order+ island+species, data=subset(all4, age>4 & extirpated ==1),FUN= length)
nspp_gri<-nspp_gri[with(nspp_gri, order(island, Order)), ]

#write out fossils with at least 5 dates
write.csv(nspp_gri[,1:3], "available_GRIWM.csv", row.names = F)
nspp_gri<-aggregate(age~Order+ island, data= nspp_gri,FUN= length)
nspp_mt1<-aggregate(age~Order+ island, data= nspp_mt1,FUN= length)
colnames(nspp_gri)[3]<-"GRIWM"
colnames(nspp_mt1)[3]<-"mt1"
all5<-merge(nspp, nspp_ext, all.x=T, all.y=T)
all5<-merge(all5, nspp_dat, all.x=T, all.y=T)
all5<-merge(all5, nspp_mt1, all.x=T, all.y=T)
all5<-merge(all5, nspp_gri, all.x=T, all.y=T)
all5[is.na(all5)] <- 0
all5$pextinct<-all5$extirpated/all5$total
all5$pdated<-with(all5, dated/total)
all5$pdated1<-with(all5,ifelse(extirpated==0, NA, dated/extirpated))
all5$pmt1 <-with(all5, ifelse(dated==0, 0, mt1*100/dated))
all5$pGRIWM <-with(all5, ifelse(dated==0, 0, GRIWM*100/dated))
all5<-all5[with(all5, order(island, Order)), ]

#write out data available
write.csv(all5, "data_available.csv", row.names = F)

#read archeo data
arch<-read.csv("archeo.csv", fileEncoding="latin1")
arc1<-subset(arch, (!is.na(arch$hage)))
arc1$lon<-arc1$lon*-1
arc2<-aggregate(hage ~ Island, data = arc1, max)
arc3 <- merge(arc1, arc2, by = c("Island", "hage"))
all6<-merge(all5, arc3, by.x="island", by.y="Island")
all6$Order <- ordered(all6$Order, levels = c("Pilosa", "Primates", "Eulipotyphla", "Chiroptera", "Rodentia"))

#plot various summaries
all6$arc<-revalue(all6$arc, c("BA"="Bahamas", "GA"="Greater Antilles", "VI"="Virgin Islands", "LA"="Lesser Antilles"))
perdat1<-ggplot(subset(all6, extirpated>0), aes(x=pdated1, y=extirpated))+geom_bar(stat = "identity", aes(fill = reorder(island, lon)))+scale_fill_viridis(discrete=T, begin=1, end=0, guide = guide_legend(title = "Islands from west to east")) + theme_bw()+facet_grid(Order~., scales="free_y", switch = "y")+scale_y_continuous(breaks=pretty_breaks(n=3), position = "right")+ scale_x_continuous(labels = scales::percent)+xlab("Percentage of extinct species dated on island")+ylab("Number of extinct species")+theme(legend.position="left")
ggsave("per_dated_extinct_by_island1.pdf", h=6, w=7)
perdat2<-ggplot(subset(all6, extirpated>0), aes(x=pdated1, y=extirpated))+geom_bar(stat = "identity", aes(fill = reorder(arc, lon)))+scale_fill_viridis(discrete=T, begin=1, end=0, guide = guide_legend(title = "Archipelago")) + theme_bw()+facet_grid(Order~., scales="free_y", switch = "y")+scale_y_continuous(breaks=pretty_breaks(n=3), position = "right")+ scale_x_continuous(labels = scales::percent)+xlab("Percentage of extinct species dated on island")+ylab("Number of extinct species")+theme(legend.position="left")
ggsave("per_dated_extinct_by_island2.pdf", h=6, w=6)
perdat2r<-ggplot(subset(all6, extirpated>0), aes(x=pdated1, y=extirpated))+geom_bar(stat = "identity", aes(fill = reorder(arc, lon))) +scale_fill_viridis(discrete=T, begin=1, end=0, guide = guide_legend(title = "Archipelago")) + theme_bw()+facet_grid(Order~.)+scale_y_continuous(breaks=pretty_breaks(n=3))+ scale_x_continuous(labels = scales::percent)+xlab("Percentage of extinct species dated on island")+ylab("Number of extinct species")+theme(legend.position="left")
ggsave("per_dated_extinct_by_island2r.pdf", h=6, w=6)
perdat3<-ggplot(all6, aes(x=pdated, y=total))+geom_bar(stat = "identity", aes(fill = reorder(island, lon)))+scale_fill_viridis(discrete=T, begin=1, end=0, guide = guide_legend(title = "Islands from west to east")) + theme_bw()+facet_grid(Order~., scales="free_y", switch = "y")+scale_y_continuous(breaks=pretty_breaks(n=3), position = "right")+ scale_x_continuous(labels = scales::percent)+xlab("Percentage of all species dated on island")+ylab("Number of species")+theme(legend.position="left")
ggsave("per_dated_by_island1.pdf", h=6, w=7)
perdat4<-ggplot(all6, aes(x=pdated, y=total))+geom_bar(stat = "identity", aes(fill = reorder(arc, lon)))+scale_fill_viridis(discrete=T, begin=1, end=0, guide = guide_legend(title = "Archipelago")) + theme_bw()+facet_grid(Order~., scales="free_y", switch = "y")+scale_y_continuous(breaks=pretty_breaks(n=3), position = "right")+ scale_x_continuous(labels = scales::percent)+xlab("Percentage of all species dated on island")+ylab("Number of species")+theme(legend.position="left")
ggsave("per_dated_by_island2.pdf", h=6, w=6)

#clean up all kinds of fossil materials
fos1$mat<-1
fos1$mat[grep("bone", fos1$Material)]<-"bone"
fos1$mat[grep("Bone", fos1$Material)]<-"bone"
fos1$mat[grep("humerus", fos1$Material)]<-"bone"
fos1$mat[grep("femur", fos1$Material)]<-"bone"
fos1$mat[grep("skull", fos1$Material)]<-"bone"
fos1$mat[grep("tibia", fos1$Material)]<-"bone"
fos1$mat[grep("mandible", fos1$Material)]<-"bone"
fos1$mat[grep("rostrum", fos1$Material)]<-"bone"
fos1$mat[grep("vertebra", fos1$Material)]<-"bone"
fos1$mat[grep("scapula", fos1$Material)]<-"bone"
fos1$mat[grep("rib", fos1$Material)]<-"bone"
fos1$mat[grep("Tyto", fos1$Material)]<-"bone"
fos1$mat[grep("guano", fos1$Material)]<-"guano"
fos1$mat[grep("shell", fos1$Material)]<-"shell"
fos1$mat[grep("charcoal", fos1$Material)]<-"charcoal"
fos1$mat[grep("Wood", fos1$Material)]<-"wood"
fos1$mat[grep("wood", fos1$Material)]<-"wood"
fos1$mat[grep("gastropod", fos1$Material)]<-"gastropod"
fos1$mat[grep("tooth", fos1$Material)]<-"tooth"
fos1$mat[grep("Teeth", fos1$Material)]<-"tooth"
fos1$mat[grep("patina", fos1$Material)]<-"patina"
fos1$mat[grep("reported", fos1$Material)]<-"not reported"
fos1$mat[grep("carapace", fos1$Material)]<-"carapace"
fos1$mat[grep("collagen", fos1$Material)]<-"collagen"
fos1$mat[grep("specimen", fos1$Material)]<-"specimen"
sum_spp<-aggregate(Species~Site.Name+island+ mat+Order+ Reference +Lab, data=fos1,FUN= length)
min_age<-aggregate(age~Site.Name+island+ mat+Order+ Reference +Lab, data=fos1,FUN= min)
max_age<-aggregate(age~Site.Name+island+ mat+Order+ Reference +Lab, data=fos1,FUN= max)
min_sig<-aggregate(Sigma~Site.Name+island+ mat+Order+ Reference +Lab, data=fos1,FUN= min)
max_sig<-aggregate(Sigma~Site.Name+island+ mat+Order+ Reference +Lab, data=fos1,FUN= max)
sum<-merge(sum_spp, min_age)
colnames(sum)[8]<-"minage"
sum<-merge(sum, max_age)
colnames(sum)[9]<-"maxage"
sum<-merge(sum, min_sig)
colnames(sum)[10]<-"minsig"
sum<-merge(sum, max_sig)
colnames(sum)[11]<-"maxsig"
colnames(sum)[3]<-"Material"
sum.m<-melt(sum, id=c(1:6,8:11))
sum1<-dcast(sum.m, Site.Name+ island+Material+minage+maxage+minsig+maxsig+Lab+Reference ~ Order)
sum1$Age<-with(sum1, ifelse(minage==maxage, paste(minage, "+/-",minsig), paste(minage, "-",maxage)))
sum1<-sum1[, c(1:3, 15, 10:14, 8:9)]
sum1<-sum1[with(sum1, order(island, Reference)), ]

#print out summary of dates
write.csv(sum1, "summary_dates.csv", row.names = F)
fos4<-aggregate(age ~ island+Species, data = subset(fos1, age<12000), min)
fos5<-aggregate(age ~ island+Species, data = subset(fos1, age<15000), max)
fos6<- merge(fos4, fos5, by = c("island", "Species"))
colnames(fos6)[3]<-"age"
fos6<- merge(fos1, fos6, by = c("island", "Species", "age"), all.x=F, all.y=F)
fos6<-fos6[,c(1:3, 5, 10:11, 13:14,16:18)]
colnames(fos6)[3]<-"age.x"
colnames(fos6)[11]<-"age"
fos6<-merge(fos1, fos6, by = c("island", "Species", "age"), all.x=F, all.y=F)
colnames(fos6)[3]<-"maxage"
colnames(fos6)[18]<-"minage"
fos6$Species<-with(fos6, ifelse(Species=="UNK", "Pilosa sp.", ifelse(Species=="Primate mandible", "Primate sp.", as.character(Species))))
fos6$mintype<-with(fos6, ifelse(Direct.Date.of.mammal.y=="N", ifelse(Archaeo.or.Palaeo.y=="Paleo", "indirect paleo", "indirect archeo" ), ifelse(Direct.Date.of.mammal.y=="Y", ifelse(Archaeo.or.Palaeo.y=="Paleo", "direct paleo", "direct archeo" ),ifelse(Archaeo.or.Palaeo.y=="Paleo", "associated paleo", "associated archeo"))))
fos6$maxtype<-with(fos6, ifelse(Direct.Date.of.mammal.x=="N", ifelse(Archaeo.or.Palaeo.x=="Paleo", "indirect paleo", "indirect archeo" ), ifelse(Direct.Date.of.mammal.x=="Y", ifelse(Archaeo.or.Palaeo.y=="Paleo", "direct paleo", "direct archeo" ),ifelse(Archaeo.or.Palaeo.x=="Paleo", "associated paleo", "associated archeo"))))
fos6$Sigma.y<-with(fos6, ifelse(island=="Marie Galante", ifelse(Order=="Chiroptera","", Sigma.y), Sigma.y))
fos6$Sigma.x<-with(fos6, ifelse(island=="Marie Galante", ifelse(Order=="Chiroptera","", Sigma.x), Sigma.x))
fos6$minAge<-with(fos6, paste(minage, "+/-",Sigma.y))
fos6$maxAge<-with(fos6, paste(maxage, "+/-",Sigma.x))
fos6$Age<-with(fos6, ifelse(minAge==maxAge,maxAge, paste(minAge,";",maxAge)))
fos6$Type<-with(fos6, ifelse(mintype==maxtype,maxtype, paste(mintype,";",maxtype)))
fos6$Reference<-with(fos6, ifelse(Reference.x==Reference.y,as.character(Reference.x), paste(Reference.y,";",Reference.x)))
fos6$Lab<-with(fos6, ifelse(Lab.x==Lab.y,as.character(Lab.x), paste(Lab.y,";",Lab.x)))
fos6<-fos6[,c(2, 1, 30:31, 33, 32)]
fos6<-fos6[with(fos6, order(island, Species)), ]

#print out alternative summary of dates by site
write.csv(fos6, "summary_dates_main.csv", row.names = F)
save.image("data_avail.Rdata")