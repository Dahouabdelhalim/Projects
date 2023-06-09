#This code is for analysizing field data 

setwd("~/UV/2016Incubations/Ecology Submission/Resubmission")

library(dplyr)
library(tidyr)
library(ggplot2)
library(nlme)
library(cowplot)
library(lme4)
library(jtools)

comp.light <- read.table("FieldDataLightEpidemics.txt", header=TRUE, sep="\\t", dec=".", strip.white=TRUE, stringsAsFactors = FALSE)

#make one dataset for each parasite      
p.comp.light<-filter(comp.light, parasite=="Pasteuria")
m.comp.light<-filter(comp.light, parasite=="Metschnikowia")

#UV index graph (Figure 3A)
#radbyday.txt not available on dryad (data from Ameriflux site Morgan Monroe State Forest, IN; Novick and Phillips 1999-present)
light<-read.table("radbyday.txt", header=TRUE, sep="\\t", dec=".", strip.white=TRUE, stringsAsFactors =FALSE)
years<-c(2014, 2015, 2016)
light2<-filter(light, year%in%years)
light3<-filter(light2, Julian>195)
light3$year<-as.factor(light3$year)

cbbPalette <- c("#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")

plot3A<-ggplot(light3, aes(x=Julian, y=maxSW, color=year))+
  geom_path(alpha=0.5)+
  theme_classic()+
  xlab("Day of the year")+
  xlim(c(190,330))+
  scale_y_continuous(limits=c(0,1000), breaks=c(0,300,600,900))+
  theme(axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_smooth(method="loess", alpha=0.9, fill=NA)+
  scale_color_manual(values = cbbPalette, name="")+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=14))+
  ylab(bquote('Max incoming SW radiation (W '*m^-2*')'))+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  coord_flip()+theme(legend.position = c(0.9,0.9))+
  theme(axis.text=element_text(color="black", size="14"))+
  annotate("text", label = c("A"), size = 4, x = 325, y = 1)


comp.light$parasite <- factor(comp.light$parasite, levels = c("Pasteuria","Metschnikowia"))

#Compare start dates of the parasites (Figure 3B)
plot3B<-ggplot(comp.light, aes(x=factor(year), y=start.date))+
  geom_boxplot(outlier.color=NULL,outlier.size=1, aes(color=parasite))+
  ylab("Epidemic start date")+
  xlab("")+
  geom_point(size=3, position=position_jitterdodge(), 
             aes(x=factor(year), y=start.date,  fill=parasite,group=parasite, shape=parasite, color=parasite))+
  theme_classic()+
  scale_color_manual(values=c("gray52","gray33"))+
  scale_fill_manual(values=c("gray52", "gray33"))+
  scale_shape_manual(values=c(24, 21))+
  theme(legend.title = element_blank())+
  theme(axis.text = element_text(size=14, color="black"))+
  theme(axis.title.y=element_text(size=14, color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=244, linetype="dashed", colour="black")+
  geom_hline(yintercept=274, linetype="dashed", colour="black")+
  geom_hline(yintercept=305, linetype="dashed", colour="black")+
  geom_hline(yintercept=213, linetype="dashed",colour="black")+
  theme(legend.text = element_text(face="italic"))+ 
  theme(legend.text = element_text(size=14))+
  theme(axis.text = element_text(color="black"))+
  scale_y_continuous(limits=c(190, 330))+
  theme(legend.position = "none")+
  theme(axis.ticks.x = element_blank())+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  annotate("text", label = c("B"), size = 4, x = 0.5, y = 325)


#This model takes years and lake as random effects (make different intercepts)
comp.light1<-filter(comp.light, !is.na(start.date))
comp.light1$fyear<-factor(comp.light1$year)
ModelStarts<-lme(start.date~parasite, random=list(~1|fyear, ~1|lake), data=comp.light1)
summary(ModelStarts)

#Median start dates in each year
S.dates.com<-comp.light1%>%group_by(year, parasite)%>%summarise(med.start=median(start.date))

#Join pasteuria and metsch datasets together and then filter out 
#lakes where epidemics didn't start to compare lakes that had epidemics
#of both parasites
slakes<-inner_join(p.comp.light, m.comp.light, by=c("lake", "year"))
slakes<-filter(slakes, !is.na(start.date.x))
slakes<-filter(slakes, !is.na(start.date.y))
slakes<-select(slakes, year, lake, start.date.x, start.date.y)
colnames(slakes)<-c("year", "lake", "Pasteuria", "Metschnikowia")
slakes.long<-gather(slakes, "parasite", "start.date", 3:4)
t.test(slakes$Pasteuria, slakes$Metschnikowia, paired=TRUE)
slakes.long$parasite <- factor(slakes.long$parasite, levels = c("Pasteuria","Metschnikowia"))

plot3C<-ggplot(slakes.long, aes(x=factor(year), y=start.date, group=parasite, shape=parasite))+
  ylab("")+xlab("")+
  geom_point(size=3,position=position_jitterdodge(jitter.height =0, jitter.width=.2), aes(x=factor(year), y=start.date,fill=parasite,color=parasite, shape=parasite))+
  scale_fill_manual(values=c("gray52","gray33")) + 
  scale_color_manual(values=c("gray52","gray33"))+
  theme_classic()+ 
  scale_shape_manual(values=c(24, 21))+
  theme(legend.title = element_blank())+
  theme(legend.text=element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(strip.background = element_rect(fill="white", colour = "white"))+
  geom_hline(yintercept=244, linetype="dashed", colour="black")+
  geom_hline(yintercept=274, linetype="dashed", colour="black")+
  geom_hline(yintercept=305, linetype="dashed", colour="black")+
  geom_hline(yintercept=213, linetype="dashed",colour="black")+
  theme(legend.text = element_text(face="italic"))+
  annotate("text", x = 3.4, y = 216, label = c("","","1 Aug"), color="black", size=3)+
  annotate("text", x = 3.4, y = 247, label = c("","","1 Sept"), color="black", size=3)+
  annotate("text", x = 3.4, y = 277, label = c("","","1 Oct"), color="black", size=3)+
  annotate("text", x = 3.4, y = 308, label = c("","","1 Nov"), color="black", size=3)+
  theme(axis.text.x = element_text(size=14, color = "black"))+
  scale_y_continuous(limits=c(190,330))+
  theme(panel.spacing.x=unit(0, "lines"))+
  theme(axis.ticks =element_blank() )+
  theme(legend.position = c(0.85,0.94))+
  theme(axis.text.y =element_blank() )+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  annotate("text", label = c("C"), size = 4, x = 0.5, y = 325)

fig3<-plot_grid(plot3A, plot3B, plot3C, ncol = 3, rel_widths = c(1.4, 1.7, 1.7), align = "h", axis="b")
save_plot("LightFigure3.jpg", fig3, base_width = 13, base_height = 4)


###############################################
#Pasteuria lake transparency analysis

#p.comp.light1 only includes lakes that had at least 1 ad320 measurement taken in a given year
p.comp.light1<-filter(p.comp.light, !is.na(m.ad320)) 
p.comp.light1$fyear<-factor(p.comp.light1$year)

#p.comp.light2 only includes lakes where epidemics that started
p.comp.light2<-filter(p.comp.light1, !is.na(start.date))

#Is start date of epidemic important to Pasteuria epidemic size?
past1<-glmer(cbind(p.comp.light2$infcount, p.comp.light2$count.at.max-p.comp.light2$infcount)~start.date+(1|fyear), family="binomial", data=p.comp.light2)
#This model fails to converge.

#Center and scale variables
p.comp.light1$est1perc.sc<-as.numeric(scale(p.comp.light1$est1perc))
p.comp.light1$start.date.sc<-as.numeric(scale(p.comp.light1$start.date))
p.comp.light1$max_Z.sc<-as.numeric(scale(p.comp.light1$max_Z))
p.comp.light1$m.tot.chl.sc<-as.numeric(scale(p.comp.light1$m.tot.chl))
p.comp.light1$m.density.sc<-as.numeric(scale(p.comp.light1$m.density))

#include observation level random effect
p.comp.light1$OLRE <- seq_len(nrow(p.comp.light1))

#Start with a global model and then drop non-significant terms until reaching lowest AIC
#AIC=334.3, only an effect of est1perc
p.glmm.1<-glmer(max.prev~
                  est1perc.sc+fyear+max_Z.sc+m.tot.chl.sc+m.density.sc+(1|OLRE), family="binomial",weight=count.at.max, data=p.comp.light1)
#AIC=333.1
p.glmm.2<-glmer(max.prev~
                  est1perc.sc+m.density.sc+m.tot.chl.sc+fyear+(1|OLRE), family="binomial",weight=count.at.max, data=p.comp.light1)
#AIC=332.8
p.glmm.3<-glmer(max.prev~
                  est1perc.sc+m.tot.chl.sc+fyear+(1|OLRE), family="binomial",weight=count.at.max, data=p.comp.light1)
#AIC=332.0
p.glmm.4<-glmer(max.prev~
                  est1perc.sc+fyear+(1|OLRE), family="binomial",weight=count.at.max, data=p.comp.light1)
#AIC=329.9
p.glmm.5<-glmer(max.prev~
                  est1perc.sc+(1|OLRE), family="binomial",weight=count.at.max, data=p.comp.light1)
drop1(p.glmm.5, test="Chisq")

#For making figure 3A
p.comp.light2<-filter(p.comp.light1, !is.na(start.date.sc))
p.comp.light2$start<-c(rep(1, length(p.comp.light2$parasite)))
both<-full_join(p.comp.light1, p.comp.light2, by=c("year","lake","parasite","start.date", "est1perc", "area", "fyear", "max_Z", "max.prev","m.tot.chl", "est1perc.sc","count.at.max"))
both$start[is.na(both$start)]<-0

A.pasteuria<-effect_plot(p.glmm.5, pred = est1perc.sc, interval=TRUE)+
  xlab("Scaled transparency")+ 
  ylab(expression("Maximum prevalence "~italic(Pasteuria)))+
  theme_classic()+
  scale_y_continuous(limits=c(0,0.13), breaks=c(0,0.03, 0.06, 0.09, 0.12))+
  geom_point(data=both, aes(est1perc.sc, max.prev, size=count.at.max, shape=fyear, color=factor(start)))+
  theme(axis.text = element_text(color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values=c("black","red"), label=c("outbreak","epidemic"), name="Outbreak category")+
  scale_shape_manual(values=c(0,2,8),name="Year")+
  scale_size_binned(name="Hosts counted")+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=14))+
  theme(axis.text = element_text(size=14, color = "black"))+
  theme(axis.title = element_text(size=14, color="black"))+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  annotate("text", label = c("A"), size = 4, x = -1.2, y = 0.125)

#t.test to compare transparency of lakes with and without epidemics
startinglakes.trans<-p.comp.light2$est1perc
p.comp.light3<-filter(p.comp.light1, is.na(start.date))
nostart.trans<-p.comp.light3$est1perc

t.test(startinglakes.trans, nostart.trans)

####################################
#Metschnikowia Analysis

#drop lakes without ad320 measurement
m.comp.light1<-filter(m.comp.light, !is.na(m.ad320))
#drop lakes without start dates
m.comp.light2<-filter(m.comp.light1, !is.na(start.date))
m.comp.light2$fyear<-factor(m.comp.light2$year)

#Center and scale variables
m.comp.light2$est1perc.sc<-as.numeric(scale(m.comp.light2$est1perc))
m.comp.light2$start.date.sc<-as.numeric(scale(m.comp.light2$start.date))
m.comp.light2$max_Z.sc<-as.numeric(scale(m.comp.light2$max_Z))
m.comp.light2$m.tot.chl.sc<-as.numeric(scale(m.comp.light2$m.tot.chl))
m.comp.light2$m.density.sc<-as.numeric(scale(m.comp.light2$m.density))

#include a random effect for each observation
m.comp.light2$OLRE <- seq_len(nrow(m.comp.light2))

#Global binomial model for Metschnikowia epidemics
#AIC=409.4; no effect of est1perc or interaction with year
m.glmm.1<-glmer(max.prev~
                est1perc.sc*fyear+start.date.sc+max_Z.sc+
                m.tot.chl.sc+m.density.sc+(1|OLRE), weight=count.at.max, family="binomial",data=m.comp.light2)

#AIC=405.7; effect of est1perc, start.date, max_z, and density (no chlorophyll or year)
m.glmm.2<-glmer(max.prev~
                est1perc.sc+start.date.sc+max_Z.sc+
                fyear+m.tot.chl.sc+m.density.sc+(1|OLRE), weight=count.at.max, family="binomial",data=m.comp.light2)

#AIC=403.9; effect of est1perc, start.date, max_Z, and density
m.glmm.3<-glmer(max.prev~
                est1perc.sc+start.date.sc+max_Z.sc+
                fyear+m.density.sc+(1|OLRE), family="binomial",weight=count.at.max, data=m.comp.light2)
drop1(m.glmm.3, test="Chisq")

#AIC=405.0; all are significant
m.glmm.4<-glmer(max.prev~
                est1perc.sc+start.date.sc+max_Z.sc+
                  m.density.sc+(1|OLRE), family="binomial",weight=count.at.max, data=m.comp.light2)


A<-effect_plot(m.glmm.3, pred = est1perc.sc, interval=TRUE)+
  xlab("Scaled transparency")+ 
  ylab(expression("Maximum prevalence "~italic(Metschnikowia)~"\\n "))+
  theme_classic()+
  geom_point(data=m.comp.light2, aes(est1perc.sc, max.prev, size=count.at.max, shape=fyear), color="red")+
  theme(axis.text = element_text(color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(limits=c(0,0.64))+
  scale_shape_manual(values=c(0,2,8),name="Year")+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size=14, color = "black"))+
  theme(axis.title = element_text(size=14, color="black"))+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1.1)+
  annotate("text", label = c("B"), size = 4, x = -1.1, y = 0.63)

B<-effect_plot(m.glmm.3, pred = start.date.sc, interval = TRUE)+
  xlab("Scaled epidemic \\n start date")+ylab("")+
  theme_classic()+
  geom_point(data=m.comp.light2, aes(start.date.sc, max.prev, size=count.at.max, shape=fyear), color="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none")+
  scale_y_continuous(limits=c(0,0.64))+
  scale_shape_manual(values=c(0,2,8),name="Year")+
  theme(axis.text = element_text(size=14, color = "black"))+
  theme(axis.title = element_text(size=14, color="black"))+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1.1)+
  annotate("text", label = c("C"), size = 4, x = -2.1, y = 0.63)

C<-effect_plot(m.glmm.3, pred = max_Z.sc, interval = TRUE)+
  xlab("Scaled lake depth")+ylab("")+
  theme_classic()+
  geom_point(data=m.comp.light2, aes(max_Z.sc, max.prev, size=count.at.max, shape=fyear), color="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none")+
  scale_y_continuous(limits=c(0,0.64))+
  scale_shape_manual(values=c(0,2,8),name="Year")+
  theme(axis.text = element_text(size=14, color = "black"))+
  theme(axis.title = element_text(size=14, color="black"))+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1.1)+
  annotate("text", label = c("D"), size = 4, x = -1.2, y = 0.63)

D<-effect_plot(m.glmm.3, pred = m.density.sc, interval = TRUE)+
  xlab("Scaled mean \\n host density")+ylab("")+
  theme_classic()+
  geom_point(data=m.comp.light2, aes(m.density.sc, max.prev, size=count.at.max, shape=fyear), color="red")+
  theme(axis.text = element_text(color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none")+
  scale_y_continuous(limits=c(0,0.64))+
  scale_shape_manual(values=c(0,2,8),name="Year")+
  theme(axis.text = element_text(size=14, color = "black"))+
  theme(axis.title = element_text(size=14, color="black"))+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1.1)+
  annotate("text", label = c("E"), size = 4, x = -1.2, y = 0.63)


top<-plot_grid(A.pasteuria, NULL, nrow=1, rel_widths = c(0.8, 1))
bottom<-plot_grid(NULL,A, B, C, D, nrow = 1, align = "h", axis="b", rel_widths = c(0.04, 1,1,1,1))

fig3<-plot_grid(top, bottom, nrow=2, align="v", axis = "l")
save_plot("LightFigure4.jpg", fig3, base_width = 10, base_height = 8)

#Supplemental Figure S4
FigS4.B<-plot_summs(m.glmm.1, m.glmm.2, m.glmm.3, m.glmm.4)
FigS4.A<-plot_summs(p.glmm.1, p.glmm.2, p.glmm.3, p.glmm.4, p.glmm.5)

FigS4<-plot_grid(FigS4.A, FigS4.B, nrow=1, labels=c("A","B"))
save_plot("LightFigureS4.jpg", FigS4, base_width = 11, base_height = 3)


#Is Metschnikowia start date influenced by transparency?
m.glmm.5<-glm(start.date~
                est1perc.sc+fyear, data=m.comp.light2)
drop1(m.glmm.5, test="F")

