
setwd("C:/Users/emchenoweth/Desktop/New folder/Step 5 Compare Models Fig")

library(readr)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

library("ggplot2", lib.loc="~/R/win-library/3.5")
is.odd <- function(x) x %% 2 != 0
devtools::install_github("seananderson/ggsidekick", force =TRUE)
1
library(ggsidekick)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#NewFig <- read.csv("NewFig!.csv")
#AllLunges <- read.csv("AllLunges.csv")

New.Fig<-results.Both
New.Fig$null<-New.Fig$Null
##Prey Escape Model
New.Fig<-New.Fig[!New.Fig$Prey.Group == "HerringTI_large_subsurface",]
New.Fig$Tactics<-New.Fig$depth
New.Fig[New.Fig$Prey.Group == "HerringTI_large_shallow",]$Tactics<-"group shallow"
New.Fig<-New.Fig[!is.na(New.Fig$Tactics),]

New.Fig.keep<-New.Fig
New.Fig<-New.Fig.keep
default.lunge.rate<-mean(AllDives[AllDives$p.max>=22.5,]$lunges_min, na.rm = TRUE)
New.Fig$Null<-default.lunge.rate*New.Fig$Null
New.Fig$B<-default.lunge.rate*New.Fig$B
New.Fig$C<-default.lunge.rate*New.Fig$C
New.Fig$D<-default.lunge.rate*New.Fig$D

New.Fig<-New.Fig[order(New.Fig$Null, decreasing= FALSE),]
pA1<-ggplot(New.Fig[!is.na(New.Fig$Tactics),]) +
  theme_sleek()+
  theme(text = element_text(size=10),
        axis.line=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        legend.position = c(0.65,0.33),
        legend.background = element_rect(colour = 'light grey', fill = 'white', linetype='solid'))+labs(color = "Prey", shape = "Tactics" )+
  geom_point(aes(New.Fig[!is.na(New.Fig$Tactics),]$Null,1:nrow(New.Fig[!is.na(New.Fig$Tactics),]), color = Prey.species, shape = New.Fig[!is.na(New.Fig$Tactics),]$Tactics), size =3) +scale_colour_manual(values=cbPalette)+
  #geom_text_repel(aes(null, B , label = New.Fig$X, color = species)) +
  #theme_classic(base_size = 10)+
  #geom_abline (lty = 3, size = 1, intercept = 0, slope = mean(New.Fig$prop.captured..no.whale.behavior)*24.8*0.84)+
  ylab("")+xlab("kJ/min")

pA2<-ggplot(New.Fig[!is.na(New.Fig$Tactics),]) +
  theme_sleek() + theme(text = element_text(size=10),
        axis.line=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")+labs(color = "Prey", shape = "Tactics" )+
  geom_point(aes(New.Fig[!is.na(New.Fig$Tactics),]$B,1:nrow(New.Fig[!is.na(New.Fig$Tactics),]), color = Prey.species, shape = New.Fig[!is.na(New.Fig$Tactics),]$Tactics), size =3) +scale_colour_manual(values=cbPalette)+
  #geom_text_repel(aes(null, B , label = New.Fig$X, color = species)) +
  #theme_classic(base_size = 10)+
  #geom_abline (lty = 3, size = 1, intercept = 0, slope = mean(New.Fig$prop.captured..no.whale.behavior)*24.8*0.84)+
  ylab("")+xlab("kJ/min")+geom_segment(y = 1:nrow(New.Fig[!is.na(New.Fig$Tactics),]), yend = 1:nrow(New.Fig[!is.na(New.Fig$Tactics),]),x = New.Fig[!is.na(New.Fig$Tactics),]$B,xend =New.Fig[!is.na(New.Fig$Tactics),]$Null)+
  scale_x_continuous(lim = c(min(c(New.Fig[!is.na(New.Fig$Tactics),]$Null,New.Fig[!is.na(New.Fig$Tactics),]$B)),max(c(New.Fig[!is.na(New.Fig$Tactics),]$Null,New.Fig[!is.na(New.Fig$Tactics),]$B))) )
  
pA3<-ggplot(New.Fig[!is.na(New.Fig$Tactics),]) +
  theme_sleek() + theme(text = element_text(size=10),
                        axis.line=element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        axis.text.y=element_blank(),
                        axis.ticks=element_blank(),
                        axis.title.y=element_blank(),
                        legend.position = "none")+labs(color = "Prey", shape = "Tactics" )+
  geom_point(aes(New.Fig[!is.na(New.Fig$Tactics),]$C,1:nrow(New.Fig[!is.na(New.Fig$Tactics),]), color = Prey.species, shape = New.Fig[!is.na(New.Fig$Tactics),]$Tactics), size =3) +scale_colour_manual(values=cbPalette)+
  #geom_text_repel(aes(null, B , label = New.Fig$X, color = species)) +
  #theme_classic(base_size = 10)+
  #geom_abline (lty = 3, size = 1, intercept = 0, slope = mean(New.Fig$prop.captured..no.whale.behavior)*24.8*0.84)+
  ylab("")+xlab("kJ/min")+geom_segment(y = 1:nrow(New.Fig[!is.na(New.Fig$Tactics),]), yend = 1:nrow(New.Fig[!is.na(New.Fig$Tactics),]),x = New.Fig[!is.na(New.Fig$Tactics),]$C,xend =New.Fig[!is.na(New.Fig$Tactics),]$B)+
  scale_x_continuous(lim = c(min(c(New.Fig[!is.na(New.Fig$Tactics),]$C,New.Fig[!is.na(New.Fig$Tactics),]$B)),max(c(New.Fig[!is.na(New.Fig$Tactics),]$B,New.Fig[!is.na(New.Fig$Tactics),]$C))) )

pA4<-ggplot(New.Fig[!is.na(New.Fig$Tactics),]) +
  theme_sleek() + theme(text = element_text(size=10),
                        axis.line=element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        axis.text.y=element_blank(),
                        axis.ticks=element_blank(),
                        axis.title.y=element_blank(),
                        legend.position = "none")+labs(color = "Prey", shape = "Tactics" )+
  geom_point(aes(New.Fig[!is.na(New.Fig$Tactics),]$D,1:nrow(New.Fig[!is.na(New.Fig$Tactics),]), color = Prey.species, shape = New.Fig[!is.na(New.Fig$Tactics),]$Tactics), size =3) +scale_colour_manual(values=cbPalette)+
  #geom_text_repel(aes(null, B , label = New.Fig$X, color = species)) +
  #theme_classic(base_size = 10)+
  #geom_abline (lty = 3, size = 1, intercept = 0, slope = mean(New.Fig$prop.captured..no.whale.behavior)*24.8*0.84)+
  ylab("")+xlab("kJ/min")+geom_segment(y = 1:nrow(New.Fig[!is.na(New.Fig$Tactics),]), yend = 1:nrow(New.Fig[!is.na(New.Fig$Tactics),]),x = New.Fig[!is.na(New.Fig$Tactics),]$D,xend =New.Fig[!is.na(New.Fig$Tactics),]$C)+
  scale_x_continuous(lim = c(min(c(New.Fig[!is.na(New.Fig$Tactics),]$D,New.Fig[!is.na(New.Fig$Tactics),]$C)),max(c(New.Fig[!is.na(New.Fig$Tactics),]$D,New.Fig[!is.na(New.Fig$Tactics),]$C))) )


pA5<-ggplot(New.Fig[!is.na(New.Fig$Tactics),]) +
  theme_sleek() + theme(text = element_text(size=10),
                        axis.line=element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        axis.text.y=element_blank(),
                        axis.ticks=element_blank(),
                        axis.title.y=element_blank(),
                        legend.position = "none")+labs(color = "Prey", shape = "Tactics" )+
  geom_vline(xintercept = -3*BMR*60,lty=2)+
  geom_point(aes(New.Fig[!is.na(New.Fig$Tactics),]$E,1:nrow(New.Fig[!is.na(New.Fig$Tactics),]), color = Prey.species, shape = New.Fig[!is.na(New.Fig$Tactics),]$Tactics), size =3) +scale_colour_manual(values=cbPalette)+
  #geom_text_repel(aes(null, B , label = New.Fig$X, color = species)) +
  #theme_classic(base_size = 10)+
  #geom_abline (lty = 3, size = 1, intercept = 0, slope = mean(New.Fig$prop.captured..no.whale.behavior)*24.8*0.84)+
  ylab("")+xlab("kJ/min")+geom_segment(y = 1:nrow(New.Fig[!is.na(New.Fig$Tactics),]), yend = 1:nrow(New.Fig[!is.na(New.Fig$Tactics),]),x = New.Fig[!is.na(New.Fig$Tactics),]$E,xend =New.Fig[!is.na(New.Fig$Tactics),]$D)+
  scale_x_continuous(lim = c(min(c(New.Fig[!is.na(New.Fig$Tactics),]$D,New.Fig[!is.na(New.Fig$Tactics),]$E)),max(c(New.Fig[!is.na(New.Fig$Tactics),]$D,New.Fig[!is.na(New.Fig$Tactics),]$E))) )


pB1<-ggplot(New.Fig)+ ggtitle("Energy Density Mod")+
  geom_point(aes(New.Fig$Null,as.numeric(is.odd(rank(New.Fig$Null, ties.method = "first")))*0.5, color = Prey.species, shape = New.Fig$Tactics), size =3) +scale_colour_manual(values=cbPalette)+
  scale_y_continuous(lim = c(-0.5,1.5))+theme_sleek()+theme(text = element_text(size=10))+
  theme(axis.line=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")+
  ylab("")+xlab("kJ/lunge")

pB2<-ggplot(New.Fig)+ ggtitle("Prey Escape Mod")+
  geom_point(aes(New.Fig$B,as.numeric(is.odd(rank(New.Fig$B, ties.method = "first")))*0.5, color = Prey.species, shape = New.Fig$Tactics), size =3) +scale_colour_manual(values=cbPalette)+
  scale_y_continuous(lim = c(-0.5,1.5))+theme_sleek()+theme(text = element_text(size=10))+
  theme(axis.line=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")+
  ylab("")#+xlab("kJ/lunge")

pB3<-ggplot(New.Fig)+ ggtitle("Whale Tactics Mod")+
  geom_point(aes(New.Fig$C,as.numeric(is.odd(rank(New.Fig$C)))*0.5, color = Prey.species, shape = New.Fig$Tactics), size =3) +scale_colour_manual(values=cbPalette)+
  scale_y_continuous(lim = c(-0.5,1.5))+theme_sleek()+theme(text = element_text(size=10))+
  theme(axis.line=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")+
  ylab("")#+xlab("kJ/lunge")

pB4<-ggplot(New.Fig)+ ggtitle("Meta. Costs Mod")+
  geom_point(aes(New.Fig$D,as.numeric(is.odd(rank(New.Fig$D)))*0.5, color = Prey.species, shape = New.Fig$Tactics), size =3) +scale_colour_manual(values=cbPalette)+
  scale_y_continuous(lim = c(-0.5,1.5))+theme_sleek()+theme(text = element_text(size=10))+
  theme(axis.line=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")+
  ylab("")#+xlab("kJ/lunge")

pB5<-ggplot(New.Fig)+ ggtitle("Lunge Rate Mod")+
  geom_point(aes(New.Fig$E,as.numeric(is.odd(rank(New.Fig$E)))*0.5, color = Prey.species, shape = New.Fig$Tactics), size =3) +scale_colour_manual(values=cbPalette)+
  scale_y_continuous(lim = c(-0.5,1.5))+theme_sleek()+theme(text = element_text(size=10))+
  theme(axis.line=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")+
  ylab("")#+xlab("kJ/lunge")

New.Fig$label<-"NA"
New.Fig$label<-c("m","l","k","j","i","h","f","g","e","c","d","a","b")

pC1<-ggplot(New.Fig)+ ggtitle("Energy Density Model")+
  geom_point(aes(rank(desc(New.Fig$Null), ties.method = "first"),1, color = Prey.species, shape = New.Fig$Tactics), size =7) +scale_colour_manual(values=cbPalette)+
  geom_text(aes(rank(desc(New.Fig$Null), ties.method = "first"),1, label = New.Fig$label), size = 5)+
  scale_y_continuous()+scale_x_continuous(breaks = seq(1, 13, by = 1))+theme_sleek()+#theme(text = element_text(size=10))+
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.position = "none")+
  ylab("")+xlab("rank")

pC2<-ggplot(New.Fig)+ ggtitle("Prey Escape Model")+
  geom_point(aes(rank(desc(New.Fig$B), ties.method = "first"),1, color = Prey.species, shape = New.Fig$Tactics), size =7) +scale_colour_manual(values=cbPalette)+
  geom_text(aes(rank(desc(New.Fig$B), ties.method = "first"),1, label = New.Fig$label), size = 5)+
  scale_y_continuous()+scale_x_continuous(breaks = seq(1, 13, by = 1))+theme_sleek()+theme(text = element_text(size=10))+
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks =element_blank(),
        axis.title=element_blank(),
        legend.position = "none")+
  ylab("")+xlab("rank")

pC3<-ggplot(New.Fig)+ ggtitle("Whale Tactics Model")+
  geom_point(aes(rank(desc(New.Fig$C), ties.method = "first"),1, color = Prey.species, shape = New.Fig$Tactics), size =7) +scale_colour_manual(values=cbPalette)+
  geom_text(aes(rank(desc(New.Fig$C), ties.method = "first"),1, label = New.Fig$label), size = 5)+
  scale_y_continuous()+scale_x_continuous(breaks = seq(1, 13, by = 1))+theme_sleek()+theme(text = element_text(size=10))+
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.position = "none")+
  ylab("")+xlab("rank")

pC4<-ggplot(New.Fig)+ ggtitle("Meta. Costs Model")+
  geom_point(aes(rank(desc(New.Fig$D), ties.method = "first"),1, color = Prey.species, shape = New.Fig$Tactics), size =7) +scale_colour_manual(values=cbPalette)+
  geom_text(aes(rank(desc(New.Fig$D), ties.method = "first"),1, label = New.Fig$label), size = 5)+
  scale_y_continuous()+scale_x_continuous(breaks = seq(1, 13, by = 1))+theme_sleek()+theme(text = element_text(size=10))+
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.position = "none")+
  ylab("")+xlab("rank")
pC5<-ggplot(New.Fig)+ ggtitle("Lunge Rate Model")+
  geom_point(aes(rank(desc(New.Fig$E), ties.method = "first"),1, color = Prey.species, shape = New.Fig$Tactics), size =7) +scale_colour_manual(values=cbPalette)+
  geom_text(aes(rank(desc(New.Fig$E), ties.method = "first"),1, label = New.Fig$label), size = 5)+
  scale_y_continuous()+scale_x_continuous(breaks = seq(1, 13, by = 1))+theme_sleek()+theme(text = element_text(size=10))+
  theme(axis.line=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")+
  ylab("")+xlab("rank")

multiplot(pA1,pA2,pA3,pA4,pA5, cols = 5)
multiplot(pB1,pB2,pB3,pB4,pB5, cols = 5)
multiplot(pC1,pC2,pC3,pC4,pC5, cols = 1)

ggsave("multibottom5.tiff",multiplot(pA1,pA2,pA3,pA4,pA5, cols = 5), width = 8.5, dpi = "retina", height = 4)
#ggsave("multitop2.tiff",multiplot(pB1,pB2,pB3,pB4,pB5, cols = 5), width = 8.5, dpi = "retina")
#ggsave("multitop3.tiff",multiplot(pC1,pC2,pC3,pC4,pC5, cols = 5), width = 8.5, dpi = "retina", height = 1)

#Table

those<-which(names(New.Fig) %in% c("Prey.Group",".energy_m3.mean","CaptureProp_WhaleTactics","WT.energy_gain.mean","lunge.rate","metabolic.costs_lunge","metabolic.costs_min","net.energy.gain_min"))
New.Fig[order(New.Fig$net.energy.gain_min, decreasing = FALSE),those]
New.Fig[order(New.Fig$net.energy.gain_min, decreasing = FALSE),c(1,28,27,45,52,56,54,59)]
Table4<-New.Fig[order(New.Fig$net.energy.gain_min, decreasing = FALSE),c(1,28,27,45,52,56,54,59)]
names(Table4)<-c("Foraging_Scenario","Patch_energy_density","Proportion_captured","Gross_energy_gain","Lunge_rate","Metabolic_costs_per_lunge","Metabolic_costs_per_time", "Net_energy_gain")
Table4$Patch_energy_density<-round(Table4$Patch_energy_density)
Table4$Gross_energy_gain<-round(Table4$Gross_energy_gain)
Table4$Metabolic_costs_per_lunge<-round(Table4$Metabolic_costs_per_lunge)
Table4$Metabolic_costs_per_time<-round(Table4$Metabolic_costs_per_time)
Table4$Net_energy_gain<-round(Table4$Net_energy_gain)
Table4$Proportion_captured<-round(Table4$Proportion_captured,2)
Table4$Lunge_rate<-round(Table4$Lunge_rate,2)


View(Table4)

### Table Supplemental G
TableG<-EnergyCosts
those<-which(names(TableG) %in% c("tag","group","dur","group.size","lunge.depth","dive.depth","max.speed","meta.costs.mean","lunge_min","meta.costs_min"))
TableG[order(TableG$group, decreasing = FALSE),those]
TableG[order(TableG$group, decreasing = FALSE),c(1,2,9,10,12,13,18,19,31,27)]
TableG<-TableG[order(TableG$group, decreasing = FALSE),c(1,2,9,10,12,13,18,19,31,27)]

TableG$dur<-round(TableG$dur/60)
TableG$lunge.depth<-round(TableG$lunge.depth)
TableG$dive.depth<-round(TableG$dive.depth)
TableG$max.speed<-round(TableG$max.speed,1)
TableG$meta.costs.mean<-round(TableG$meta.costs.mean)
TableG$lunge_min<-round(TableG$lunge_min,2)
TableG$ meta.costs_min<-round(TableG$meta.costs_min)
View(TableG)

View(Table4)


