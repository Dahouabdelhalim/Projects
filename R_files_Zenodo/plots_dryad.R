##########################################################################
#
#	Script to make plots from meta-analyses of West Indian extinct mammal dates
#
#	Files required in same folder
#	prior results: dates.Rdata
#
##########################################################################
library(reshape2)
library(ggplot2)
library(viridis)
library(scales)
library(lattice)
library(plyr)
library(dplyr)
library(grid)
library(R2jags)
rm(list = ls())

#load prior results
load("dates.Rdata")
arc4<-arc3[,1:5]
colnames(arc4)[1]<-"island"
colnames(arc4)[2]<-"age"
fos4<-subset(fos3, age<12000)
fos4$Species<-ifelse(fos4$Species=="UNK", "Pilosa sp.", ifelse(fos4$Species=="Primate mandible", "Primate sp.", as.character(fos4$Species)))
fos4<-fos4[,c(1:3,5:7,11)]
fos4$Species <- factor(fos4$Species, levels = unique(fos4$Species[order(fos4$Species)]))
fos4$island <- ordered(fos4$island, levels = unique(alld$island))
fos4<-subset(fos4, island !="Guana")
arc4$Species<-"Homo sapiens"
arc4$Order<-"Primates"
arc5<-subset(arc3, Island !="Aruba" & Island !="Curacao")
arc5$arc<-revalue(arc5$arc, c("BA"="Bahamas", "GA"="Greater Antilles", "LA"="Lesser Antilles", "VI"="Virgin Islands"))
arc5$arc<-ordered(arc5$arc, levels  = c("Bahamas", "Greater Antilles", "Virgin Islands",  "Lesser Antilles"))
colnames(arc5)[6]<-"Latitude"
age_plot_human1<-ggplot(arc5, aes(x = reorder(Island, lon), fill=Latitude)) +
   geom_boxplot(aes(
       lower = hage - Sigma, 
       upper = hage + Sigma, 
       middle = hage, 
       ymin = hage - 3* Sigma, 
       ymax = hage + 3* Sigma),
     stat = "identity")+theme_bw()+facet_grid(.~arc, drop=TRUE, scales="free_x", switch = "y", space = "free_x")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_fill_viridis(option="C")+ylab("First appearance humans in years before present")+xlab("Islands from west to east")+geom_hline(yintercept = 2000)
ggsave("age_plot_human1.pdf", h=6.5, w=9.5)
age_plot_human2<-ggplot(arc5, aes(x = reorder(Island, lon), fill=Latitude)) +
   geom_boxplot(aes(
       lower = hage - Sigma, 
       upper = hage + Sigma, 
       middle = hage, 
       ymin = hage - 3* Sigma, 
       ymax = hage + 3* Sigma),
     stat = "identity")+theme_bw()+facet_grid(.~arc, drop=TRUE, scales="free_x", switch = "y", space = "free_x")+theme(axis.text.x = element_text(angle = 45, hjust = 1, colour="grey10"))+scale_fill_viridis(option="D")+ylab("First appearance humans in years before present")+xlab("Islands from west to east")+geom_hline(yintercept = 2000)
ggsave("age_plot_human2.pdf", h=6.5, w=9.5)
all2<-rbind(fos4,arc4[arc4$island %in% fos4$island,])
all2$Date<-ifelse(all2$Species=="Homo sapiens", "Human first appearance", "Fauna last appearance")
all2$Type <- ordered(all2$Date, levels = c("Human first appearance", "Fauna last appearance"))
all2<-merge(all2, arc5, by.x="island", by.y="Island", all.x=T)
all2<-all2[,c(1:9,14:16)]
colnames(all2)[7]<-"Sigma"
reverselog_trans <- function(base = exp(1))	{
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
              								}
age_plot<-ggplot(all2, aes(x = Species, fill=Date)) +
   geom_boxplot(aes(
       lower = age - Sigma, 
       upper = age + Sigma, 
       middle = age, 
       ymin = ifelse(age - 3* Sigma>0,age - 3* Sigma,100), 
       ymax = age + 3* Sigma),
     stat = "identity")+coord_flip()+theme_bw()+facet_grid(island~., drop=TRUE, scales="free_y", switch = "y", space = "free_y")+theme(axis.text.y = element_text(face="italic"), legend.position="left") +scale_fill_viridis(discrete = TRUE, option="C")+scale_y_continuous(breaks = c(100, 1000, 10000), trans=reverselog_trans(10))+ylab("Calibrated age in years before present")+scale_x_discrete(position = "top")
ggsave("age_plot_log.pdf", h=15, w=7.5)
age_plot1<-ggplot(all2, aes(x = Species, fill=Date)) +
   geom_boxplot(aes(
       lower = age - Sigma, 
       upper = age + Sigma, 
       middle = age, 
       ymin = ifelse(age - 3* Sigma>0,age - 3* Sigma,100), 
       ymax = age + 3* Sigma),
     stat = "identity")+theme_bw()+facet_wrap(~reorder(island, -Latitude), ncol=10, drop=TRUE, scales="free_x", strip.position = "left", dir="v")+theme(axis.text.x = element_text(face="italic",angle = 45, hjust = 1), legend.position="right") +scale_fill_viridis(discrete = TRUE, option="C")+scale_y_continuous(breaks = c(100, 1000, 10000), trans="log10")+ylab("Calibrated age in years before present")
N1<-all2%>% group_by(reorder(island, -Latitude))%>% summarise(count = length(unique(Species)))
gt <- ggplotGrob(age_plot1)
panelI <- gt$layout$l[grepl("panel", gt$layout$name)]
gt$widths[panelI] <- unit(N1$count, "null")
#gt$widths[panelI][4]<-unit(13,"null")
ggsave("age_plot_log1.pdf", h=8, w=11, gt)
age_plot2<-ggplot(all2, aes(x = Species, fill=Date)) +
   geom_boxplot(aes(
       lower = age - Sigma, 
       upper = age + Sigma, 
       middle = age, 
       ymin = age - 3* Sigma, 
       ymax = age + 3* Sigma),
     stat = "identity")+ coord_flip()+theme_bw()+facet_grid(island~., drop=TRUE, scales="free_y", switch = "y", space = "free_y")+theme(axis.text.y = element_text(face="italic"), legend.position="right") +scale_fill_viridis(discrete = TRUE, option="C")+ylab("Calibrated age in years before present")+scale_y_reverse()+scale_x_discrete(position = "top")
ggsave("age_plot.pdf", h=16, w=8)
age_plot3<-ggplot(all2, aes(x = Species, fill=Date, colour=Date)) +
   geom_boxplot(aes(
       lower = age - Sigma, 
       upper = age + Sigma, 
       middle = age, 
       ymin = age - 3* Sigma, 
       ymax = age + 3* Sigma),
     stat = "identity")+theme_bw()+facet_wrap(~reorder(island, -Latitude), ncol=10, drop=TRUE, scales="free_x", strip.position = "left", dir="v")+theme(axis.text.x = element_text(face="italic",angle = 45, hjust = 1, colour="grey10", size=7), legend.position="right", strip.text.y = element_text(colour="white"), strip.background = element_rect(fill=viridis(4, option="D")[1])) +scale_fill_viridis(discrete = TRUE, option="C")+scale_color_viridis(discrete = TRUE, option="C")+ylab("Calibrated age in years before present")
gt <- ggplotGrob(age_plot3)
panelI <- gt$layout$l[grepl("panel", gt$layout$name)]
gt$widths[panelI] <- unit(N1$count, "null")
#gt$widths[panelI][4]<-unit(13,"null")
ggsave("age_plot3.pdf", h=8, w=11, gt)
age_plot4<-ggplot(all2, aes(x = Species, fill=Date, colour=Date)) +
   geom_boxplot(aes(
       lower = age - Sigma, 
       upper = age + Sigma, 
       middle = age, 
       ymin = age - 3* Sigma, 
       ymax = age + 3* Sigma),
     stat = "identity")+theme_bw()+facet_wrap(~reorder(island, -Latitude), ncol=10, drop=TRUE, scales="free_x", strip.position = "left", dir="v")+theme(axis.text.x = element_text(face="italic",angle = 45, hjust = 1,colour="grey10", size=7), legend.position="right") +scale_fill_viridis(discrete = TRUE, option="C", begin = 0, end = 0.8)+scale_color_viridis(discrete = TRUE, option="C", begin = 0, end = 0.8)+ylab("Calibrated age in years before present")
gt <- ggplotGrob(age_plot4)
panelI <- gt$layout$l[grepl("panel", gt$layout$name)]
gt$widths[panelI] <- unit(N1$count, "null")
gt$widths[panelI][6]<-unit(20,"null")
gt$widths[panelI][9]<-unit(8,"null")
gt$widths[panelI][12]<-unit(12,"null")
gt$widths[panelI][15]<-unit(12,"null")
gt$widths[panelI][18]<-unit(24,"null")
gt$widths[panelI][21]<-unit(12,"null")
ggsave("age_plot4.pdf", h=8, w=11, gt)
attach.jags(car.3)
islands<-as.data.frame(a)
detach.jags()
colnames(islands)<-unique(alld$island)
islands.m<-melt(islands)
colnames(islands.m)<-c("island", "difference")
islands.m<-merge(islands.m, arc5, by.x="island", by.y="Island", all.x=T)
est_plot1<-qplot(difference, data=islands.m, geom="density", fill= arc)+theme_bw()+facet_wrap(~ reorder(island, -Latitude), strip.position = "left", ncol=3)+scale_fill_viridis(discrete = TRUE, option="D")+xlab("Difference between human first appearance and faunal last appearance")+ylab("Posterior frequency distribution")+theme(legend.position ="none", axis.text.y=element_blank(), axis.ticks.y=element_blank())+scale_y_continuous(position = "right", breaks= c(0, 0.00025, 0.0005))+geom_vline(aes(xintercept = 0), colour="grey50")
ggsave("estimate_plot1.pdf", h=10, w=8)
attach.jags(car.5)
islands<-as.data.frame(a)
detach.jags()
colnames(islands)<-unique(alld$island2)
islands.m<-melt(islands)
colnames(islands.m)<-c("island", "difference")
islands.m$island2<-with(islands.m, ifelse(island=="Cuba_b", "Cuba", ifelse(island=="Cuba_a", "Cuba", ifelse(island=="Hispaniola_a", "Hispaniola", ifelse(island=="Hispaniola_b", "Hispaniola", as.character(island))))))
islands.m<-merge(islands.m, arc5, by.x="island2", by.y="Island", all.x=T)
islands.m$island2<-with(islands.m, ifelse(island=="Cuba_b", "Cuba post", ifelse(island=="Cuba_a", "Cuba pre", ifelse(island=="Hispaniola_a", "Hispaniola pre", ifelse(island=="Hispaniola_b", "Hispaniola post", as.character(island))))))
islands.m$island2 <- with(islands.m, ordered(island2, levels=unique(unique(islands.m$island2)[c(1:5, 7,6,8:9, 11, 10, 12:23)]))) 
est_plot2<-qplot(difference, data=islands.m, geom="density", fill= arc)+theme_bw()+facet_wrap(~reorder(island2, lon), strip.position = "right", ncol=3)+scale_fill_viridis(discrete = TRUE, option="D", begin=1, end =0, guide = guide_legend(title = "Archipelago"))+xlab("Difference between human first appearance and faunal last appearance")+ylab("Posterior frequency distribution")+theme(legend.position ="left", axis.text.y=element_blank(), axis.ticks.y=element_blank())+scale_y_continuous(position = "right", breaks= c(0, 0.00025, 0.0005))+geom_vline(aes(xintercept = 0), colour="grey50")
ggsave("estimate_plot2.pdf", h=10, w=8)
save.image("dates.Rdata")
