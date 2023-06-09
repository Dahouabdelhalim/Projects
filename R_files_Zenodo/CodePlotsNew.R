library(ggplot2)
library(latex2exp)
library(ggpubr)

#FIGURE 1
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure1/5-2.5/base/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/Figure1/5-2.5/base/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/Figure1/5-2.5/base/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS_A",length(epi.data$FinSize)))

scaling.parm<-400
data.525<-data.frame("Label"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.strategy))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "Strategy"=c(as.character(temp.strategy),as.character(temp.strategy)))

temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure1/5-8/base/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/Figure1/5-8/base/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/Figure1/5-8/base/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS_A",length(epi.data$FinSize)))

scaling.parm2<-400
data.58<-data.frame("Label"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.strategy))), "Value"=c(temp.fs,temp.nsdl/scaling.parm2), "Strategy"=c(as.character(temp.strategy),as.character(temp.strategy)))

data.overall<-data.frame("Scenario"=c(rep("Wuhan",length(data.525$Label)),rep("Delta",length(data.525$Label))), "SummaryMeasure"=c(as.character(data.525$Label), as.character(data.58$Label)), "Value"=c(data.525$Value,data.58$Value), "Strategy"=c(as.character(data.525$Strategy),as.character(data.58$Strategy))     )
data.overall$Scenario<-factor(data.overall$Scenario, levels = c("Wuhan","Delta"))
data.overall$Strategy<-factor(data.overall$Strategy, levels = c("SI","RS","RS_A"))

ggplot(data = data.overall, aes(x=Strategy, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15,17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),6), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+facet_wrap( ~ Scenario, scales = "fixed")+scale_x_discrete(labels=c("SI","ReaS","RepS"))+labs(fill="Model Outcome")
#save 9:5.5

#median and quantile
#SI-Wuhan
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure1/5-2.5/base/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

#SI-Delta
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure1/5-8/base/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

#Rs-Wuhan
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure1/5-2.5/base/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

#RS-Delta
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure1/5-8/base/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

#Rs-Wuhan
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure1/5-2.5/base/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

#RS-Delta
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure1/5-8/base/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)



#FIGURE 2
scaling.parm<-100

temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure2/sa6/20/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("2",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/Figure2/sa6/22/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("4",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/Figure2/sa6/24/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("6",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/Figure2/sa6/26/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("8",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/Figure2/sa6/28/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("1000",length(epi.data$FinSize)))


data.58SA6<-data.frame("SummaryMeasure"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs ,temp.nsdl/scaling.parm), "Threshold"=c(as.character(temp.thresh),as.character(temp.thresh))   )
data.58SA6$Threshold<-factor(data.58SA6$Threshold, levels = c("2","4","6","8","1000"))
ggplot(data = data.58SA6, aes(x=Threshold, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),5), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Class Closure Threshold")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+labs(fill="Model Outcome")+scale_x_discrete(labels = c("2","4","6","8","No Threshold"))


#Median and Quantile

temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure2/sa6/20/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)




temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure2/sa6/22/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure2/sa6/24/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure2/sa6/26/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure2/sa6/28/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)



#FIGURE 3
scaling.parm<-30
temp.fs<-NULL
temp.nsdl<-NULL
temp.test<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure3/5-8/sa7/22/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.test<-c(temp.test,rep("1",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/Figure3/5-8/sa7/46/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.test<-c(temp.test,rep("2",length(epi.data$FinSize)))

data.58SA7<-data.frame("SummaryMeasure"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "NTest"=rep(temp.test,2))

ggplot(data = data.58SA7, aes(x=NTest, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15), size=3, color=rep(c("royalblue","gold"),2), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Number of Tests Per Week")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+labs(fill="Model Outcome")

temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure2/sa6/26/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure3/sa7/10/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

# Figure 3B (Wuhan Vs Delta twice vs single testing)
temp.fs<-NULL
temp.nsdl<-NULL
temp.nt<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure3/5-8/sa7/16/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.nt<-c(temp.nt,rep("1",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/Figure3/5-8/sa7/40/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.nt<-c(temp.nt,rep("2",length(epi.data$FinSize)))
scaling.parm<-300
data.ntD<-data.frame("Label"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "Number of Tests per Week"=c(rep(temp.nt,2)))

scaling.parm2<-300
temp.fs<-NULL
temp.nsdl<-NULL
temp.nt<-NULL
epi.data<-read.csv("~/Desktop/Testing/Figure3/5-2.5/sa7/16/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.nt<-c(temp.nt,rep("1",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/Figure3/5-2.5/sa7/40/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.nt<-c(temp.nt,rep("2",length(epi.data$FinSize)))

data.ntW<-data.frame("Label"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs,temp.nsdl/scaling.parm2), "Strategy"=c(rep(temp.nt,2)))

data.overall<-data.frame("Seeds"=c(rep("Wuhan",length(data.ntW$Label)),rep("Delta",length(data.ntD$Label))), "SummaryMeasure"=c(as.character(data.ntW$Label), as.character(data.ntD$Label)), "Value"=c(data.ntW$Value,data.ntD$Value), "Strategy"=c(as.character(data.ntW$Strategy),as.character(data.ntD$Strategy))     )
data.overall$Scenario<-factor(data.overall$Seeds, levels = c("Wuhan","Delta"))
data.overall$Strategy<-factor(data.overall$Strategy, levels = c("1","2"))



ggplot(data = data.overall, aes(x=Strategy, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),4), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Number of Test Per Week")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+facet_wrap( ~ Scenario, scales = "fixed")+scale_x_discrete(labels=c("1","2"))+labs(fill="Model Outcome")






#FIGURE SA1
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA1/1-8/base/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA1/1-8/base/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA1/1-8/base/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS_A",length(epi.data$FinSize)))

scaling.parm<-300
data.18<-data.frame("Label"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "Strategy"=c(rep(temp.strategy,2)))

scaling.parm2<-300
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA1/10-8/base/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA1/10-8/base/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA1/10-8/base/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS_A",length(epi.data$FinSize)))

data.108<-data.frame("Label"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs,temp.nsdl/scaling.parm2), "Strategy"=c(rep(temp.strategy,2)))

data.overall<-data.frame("Seeds"=c(rep("1 Seed",length(data.18$Label)),rep("10 Seeds",length(data.108$Label))), "SummaryMeasure"=c(as.character(data.18$Label), as.character(data.108$Label)), "Value"=c(data.18$Value,data.108$Value), "Strategy"=c(as.character(data.18$Strategy),as.character(data.108$Strategy))     )
data.overall$Scenario<-factor(data.overall$Seeds, levels = c("1 Seed","10 Seeds"))
data.overall$Strategy<-factor(data.overall$Strategy, levels = c("SI","RS","RS_A"))



ggplot(data = data.overall, aes(x=Strategy, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15,17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),6), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+facet_wrap( ~ Scenario, scales = "fixed")+scale_x_discrete(labels=c("SI","ReaS","RepS"))+labs(fill="Model Outcome")



#1-seed SM
temp.fs<-NULL
temp.nsdl<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA1/1-8/base/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)


temp.fs<-NULL
temp.nsdl<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA1/1-8/base/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA1/1-8/base/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

#10 Seeds
temp.fs<-NULL
temp.nsdl<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA1/10-8/base/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)


temp.fs<-NULL
temp.nsdl<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA1/10-8/base/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA1/10-8/base/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)









#FIGURE SA2
temp.fs<-NULL
temp.nsdl<-NULL
temp.pl<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.pl<-c(temp.pl,rep("0.2",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.pl<-c(temp.pl,rep("0.2",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.pl<-c(temp.pl,rep("0.2",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.pl<-c(temp.pl,rep("0.5",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/8/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.pl<-c(temp.pl,rep("0.5",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/10/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.pl<-c(temp.pl,rep("0.5",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/12/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.pl<-c(temp.pl,rep("0.9",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/14/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.pl<-c(temp.pl,rep("0.9",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/16/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.pl<-c(temp.pl,rep("0.9",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))


finalSizedata.sa1<-data.frame(output=temp.fs,prop_lambda_w=temp.pl,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
finalSizedata.sa1$strategy<-factor(finalSizedata.sa1$strategy, levels = c("SI","RS","RS_A"))
schooldaylost.sa1<-data.frame(output=temp.nsdl,prop_lambda_w=temp.pl,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
schooldaylost.sa1$strategy<-factor(schooldaylost.sa1$strategy, levels = c("SI","RS","RS_A"))

sa1fs<-ggplot(data = finalSizedata.sa1, aes(x=strategy, fill=prop_lambda_w,y=output))+geom_boxplot()+scale_fill_brewer(palette="Accent", name="Between-Classes Contact Rate", labels=c(unname(TeX(c("$0.2\\\\lambda_b$", "$0.5\\\\lambda_b$","$0.9\\\\lambda_b$")))) )+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

sa1nsdl<-ggplot(data = schooldaylost.sa1, aes(x=strategy, fill=prop_lambda_w,y=output))+geom_boxplot()+scale_fill_brewer(palette="Accent", name="Between-Classes Contact Rate", labels=c(unname(TeX(c("$0.2\\\\lambda_b$", "$0.5\\\\lambda_b$","$0.9\\\\lambda_b$")))) )+stat_summary(fun.y=mean, geom="point", shape=15, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("NSDL"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

ggarrange(sa1fs,sa1nsdl,legend = "top",common.legend = TRUE)

#SM
temp.fs<-NULL
temp.nsdl<-NULL
temp.pl<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)


temp.fs<-NULL
temp.nsdl<-NULL
temp.pl<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
temp.pl<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
temp.pl<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
temp.pl<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/8/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
temp.pl<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/10/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
temp.pl<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/12/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
temp.pl<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/14/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
temp.pl<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA2/sa1/16/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)



#FIGURE SA3
#Immune Proportion Children
temp.fs<-NULL
temp.nsdl<-NULL
temp.ic<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA3/sa2/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ic<-c(temp.ic,rep("0.2",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA3/sa2/10/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ic<-c(temp.ic,rep("0.2",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA3/sa2/18/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ic<-c(temp.ic,rep("0.2",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA3/sa2/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ic<-c(temp.ic,rep("0.3",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA3/sa2/12/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ic<-c(temp.ic,rep("0.3",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA3/sa2/20/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ic<-c(temp.ic,rep("0.3",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA3/sa2/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ic<-c(temp.ic,rep("0.4",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA3/sa2/14/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ic<-c(temp.ic,rep("0.4",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA3/sa2/22/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ic<-c(temp.ic,rep("0.4",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))


finalSizedata.sa3<-data.frame(output=temp.fs,prop_immune_child=temp.ic,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
finalSizedata.sa3$strategy<-factor(finalSizedata.sa3$strategy, levels = c("SI","RS","RS_A"))
schooldaylost.sa3<-data.frame(output=temp.nsdl,prop_immune_child=temp.ic,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
schooldaylost.sa3$strategy<-factor(schooldaylost.sa3$strategy, levels = c("SI","RS","RS_A"))

sa3fs.ch<-ggplot(data = finalSizedata.sa3, aes(x=strategy, fill=prop_immune_child,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="Immune Proportion Children", labels=c(unname(TeX(c("$0.2", "$0.3$","$0.4$")))) )+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

sa3nsdl.ch<-ggplot(data = schooldaylost.sa3, aes(x=strategy, fill=prop_immune_child,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="Proportion Immune Children", labels=c(unname(TeX(c("$0.2$", "$0.3$","$0.4$")))) )+stat_summary(fun.y=mean, geom="point", shape=15, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("NSDL"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

ggarrange(sa3fs.ch,sa3nsdl.ch,legend = "top",common.legend = TRUE)

#adult FIG SA4 - ImmuneProportion Adults
temp.fs<-NULL
temp.nsdl<-NULL
temp.ia<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA4/sa3/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ia<-c(temp.ia,rep("0.3",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA4/sa3/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ia<-c(temp.ia,rep("0.3",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA4/sa3/12/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ia<-c(temp.ia,rep("0.3",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA4/sa3/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ia<-c(temp.ia,rep("0.5",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA4/sa3/8/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ia<-c(temp.ia,rep("0.5",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA4/sa3/14/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ia<-c(temp.ia,rep("0.5",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA4/sa3/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ia<-c(temp.ia,rep("0.9",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA4/sa3/10/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ia<-c(temp.ia,rep("0.9",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA4/sa3/16/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ia<-c(temp.ia,rep("0.9",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))


finalSizedata.sa4<-data.frame(output=temp.fs,prop_immune_ad=temp.ic,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
finalSizedata.sa4$strategy<-factor(finalSizedata.sa4$strategy, levels = c("SI","RS","RS_A"))
schooldaylost.sa4<-data.frame(output=temp.nsdl,prop_immune_ad=temp.ic,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
schooldaylost.sa4$strategy<-factor(schooldaylost.sa4$strategy, levels = c("SI","RS","RS_A"))

sa4.ad<-ggplot(data = finalSizedata.sa4, aes(x=strategy, fill=prop_immune_ad,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="Immune Proportion Adults", labels=c(unname(TeX(c("$0.3", "$0.5$","$0.9$")))) )+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

sa4nsdl.ad<-ggplot(data = schooldaylost.sa4, aes(x=strategy, fill=prop_immune_ad,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="Proportion Immune Adults", labels=c(unname(TeX(c("$0.3$", "$0.5$","$0.9$")))) )+stat_summary(fun.y=mean, geom="point", shape=15, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("NSDL"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

ggarrange(sa4.ad,sa4nsdl.ad,legend = "top",common.legend = TRUE)


#FIGURE SA5 - Infectivity Asymptomatic Carriers
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureSA5/sa4/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA5/sa4/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureSA5/sa4/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS_A",length(epi.data$FinSize)))


finalSizedata.sa5<-data.frame(output=temp.fs, strategy=temp.strategy)
finalSizedata.sa5$strategy<-factor(finalSizedata.sa5$strategy, levels = c("SI","RS","RS_A"))
finalSizedata.sa5$gruop1<-c(rep("1",100),rep("2",100),rep("3",100))

schooldaylost.sa5<-data.frame(output=temp.nsdl, strategy=temp.strategy)
schooldaylost.sa5$strategy<-factor(schooldaylost.sa5$strategy, levels = c("SI","RS","RS_A"))
schooldaylost.sa5$gruop1<-c(rep("1",100),rep("2",100),rep("3",100))


sa5fs<-ggplot(data = finalSizedata.sa5, aes(x=strategy, fill=strategy,y=output))+geom_boxplot()+scale_fill_brewer(palette="Accent")+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "none",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

sa5nsdl<-ggplot(data = schooldaylost.sa5, aes(x=strategy, fill=strategy,y=output))+geom_boxplot()+scale_fill_brewer(palette="Accent")+stat_summary(fun.y=mean, geom="point", shape=15, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "none",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("NSDL"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

ggarrange(sa5fs,sa5nsdl,legend = "none")

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################

# NEW PLOTS (REVISION)

#########################
# Numbers of reactive screening

scaling.parm<-30
temp.fs<-NULL
temp.nsdl<-NULL
temp.test<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR1/sa6/16/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.test<-c(temp.test,rep("1",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR1/sa11/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.test<-c(temp.test,rep("2",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR1/sa11/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.test<-c(temp.test,rep("3",length(epi.data$FinSize)))

data.R1<-data.frame("SummaryMeasure"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "NTest"=rep(temp.test,2))
ggplot(data = data.R1, aes(x=NTest, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),3), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Number of Class Screening")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+labs(fill="Model Outcome")


temp.fs<-NULL
temp.nsdl<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR1/sa6/16/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR1/sa11/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

temp.fs<-NULL
temp.nsdl<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR1/sa11/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
median(temp.fs)
quantile(temp.fs, 0.025)
quantile(temp.fs, 0.975)
median(temp.nsdl)
quantile(temp.nsdl, 0.025)
quantile(temp.nsdl, 0.975)

#######################3
# School Threshold - RepS

scaling.parm<-500
temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/28/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("2",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/30/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("5",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/32/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("10",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/34/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("20",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/36/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("50",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/38/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("100",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/40/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("1000",length(epi.data$FinSize)))


data.R2<-data.frame("SummaryMeasure"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "Threshold"=rep(temp.thresh,2))
data.R2$Threshold<-factor(data.R2$Threshold, levels = c("2","5","10" ,"20","50","100","1000"))

ggplot(data = data.R2, aes(x=Threshold, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15,17,15,17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),7), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("School Closure Threshold")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+labs(fill="Model Outcome")+scale_x_discrete(labels = c("2","5","10","20","50","100","No Threshold"))


# ReaS

scaling.parm<-500
temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/14/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("2",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/16/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("5",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/18/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("10",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/20/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("20",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/22/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("50",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/24/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("100",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/26/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("1000",length(epi.data$FinSize)))


data.R2.ReaS<-data.frame("SummaryMeasure"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "Threshold"=rep(temp.thresh,2))
data.R2.ReaS$Threshold<-factor(data.R2.ReaS$Threshold, levels = c("2","5","10" ,"20","50","100","1000"))

ggplot(data = data.R2.ReaS, aes(x=Threshold, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15,17,15,17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),7), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("School Closure Threshold")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+labs(fill="Model Outcome")+scale_x_discrete(labels = c("2","5","10","20","50","100","No Threshold"))

# Sympt Iso

scaling.parm<-500
temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("2",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("5",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("10",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("20",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/8/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("50",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/10/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("100",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR2/sa8/12/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("1000",length(epi.data$FinSize)))


data.R2.SI<-data.frame("SummaryMeasure"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "Threshold"=rep(temp.thresh,2))
data.R2.SI$Threshold<-factor(data.R2.SI$Threshold, levels = c("2","5","10" ,"20","50","100","1000"))

ggplot(data = data.R2.SI, aes(x=Threshold, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15,17,15,17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),7), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("School Closure Threshold")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+labs(fill="Model Outcome")+scale_x_discrete(labels = c("2","5","10","20","50","100","No Threshold"))


#FIGURE R3 - Seeding Time
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR3/sa12/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR3/sa12/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR3/sa12/10/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS_A",length(epi.data$FinSize)))


finalSizedata.R3<-data.frame(output=temp.fs, strategy=temp.strategy)
finalSizedata.R3$strategy<-factor(finalSizedata.R3$strategy, levels = c("SI","RS","RS_A"))
finalSizedata.R3$gruop1<-c(rep("1",100),rep("2",100),rep("3",100))

schooldaylost.R3<-data.frame(output=temp.nsdl, strategy=temp.strategy)
schooldaylost.R3$strategy<-factor(schooldaylost.R3$strategy, levels = c("SI","RS","RS_A"))
schooldaylost.R3$gruop1<-c(rep("1",100),rep("2",100),rep("3",100))

R3fs<-ggplot(data = finalSizedata.R3, aes(x=strategy, fill=strategy,y=output))+geom_boxplot()+scale_fill_brewer(palette="Accent")+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "none",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

R3nsdl<-ggplot(data = schooldaylost.R3, aes(x=strategy, fill=strategy,y=output))+geom_boxplot()+scale_fill_brewer(palette="Accent")+stat_summary(fun.y=mean, geom="point", shape=15, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "none",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("NSDL"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

ggarrange(R3fs,R3nsdl,legend = "none")


################################################
# R4 : probability of showing symptoms
temp.fs<-NULL
temp.nsdl<-NULL
temp.sp<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR4/sa13/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sp<-c(temp.sp,rep("0.2",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR4/sa13/8/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sp<-c(temp.sp,rep("0.2",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR4/sa13/16/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sp<-c(temp.sp,rep("0.2",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR4/sa13/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sp<-c(temp.sp,rep("0.4",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR4/sa13/10/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sp<-c(temp.sp,rep("0.4",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR4/sa13/18/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sp<-c(temp.sp,rep("0.4",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR4/sa13/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sp<-c(temp.sp,rep("0.8",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR4/sa13/14/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sp<-c(temp.sp,rep("0.8",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR4/sa13/22/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sp<-c(temp.sp,rep("0.8",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))


#finalSizedata.R4<-data.frame(output=temp.fs,prob_sympt=temp.sp,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100),rep("10",100),rep("11",100),rep("12",100)))
finalSizedata.R4<-data.frame(output=temp.fs,prob_sympt=temp.sp,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
finalSizedata.R4$strategy<-factor(finalSizedata.R4$strategy, levels = c("SI","RS","RS_A"))
finalSizedata.R4$prob_sympt<-factor(finalSizedata.R4$prob_sympt, levels = c("0.2","0.4","0.8"))
#schooldaylost.R4<-data.frame(output=temp.nsdl,prob_sympt=temp.sp,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100),rep("10",100),rep("11",100),rep("12",100)))
schooldaylost.R4<-data.frame(output=temp.nsdl,prob_sympt=temp.sp,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
schooldaylost.R4$strategy<-factor(schooldaylost.R4$strategy, levels = c("SI","RS","RS_A"))

R4fs.ch<-ggplot(data = finalSizedata.R4, aes(x=strategy, fill=prob_sympt,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="Probability of Symptomatic Infections", labels=c(unname(TeX(c("$0.2", "$0.4$","$0.8$")))) )+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

R4nsdl.ch<-ggplot(data = schooldaylost.R4, aes(x=strategy, fill=prob_sympt,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="Probability of symptomatic infections", labels=c(unname(TeX(c("$0.2$", "$0.4$","$0.8$")))) )+stat_summary(fun.y=mean, geom="point", shape=15, size=3, color="royalblue", fill="black",position=position_dodge(0.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("NSDL"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

ggarrange(R4fs.ch,R4nsdl.ch,legend = "top",common.legend = TRUE)


######################################################
# Reducing Incubation Period

temp.fs<-NULL
temp.nsdl<-NULL
temp.ip<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR5/sa6/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ip<-c(temp.ip,rep("Baseline",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR5/sa6/16/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ip<-c(temp.ip,rep("Baseline",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR5/sa6/26/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ip<-c(temp.ip,rep("Baseline",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR5/sa6/7/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ip<-c(temp.ip,rep("ShorterIP",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR5/sa6/17/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ip<-c(temp.ip,rep("ShorterIP",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR5/sa6/27/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.ip<-c(temp.ip,rep("ShorterIP",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))


finalSizedata.R5<-data.frame(output=temp.fs,inc_per=temp.ip,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100)))
finalSizedata.R5$strategy<-factor(finalSizedata.R5$strategy, levels = c("SI","RS","RS_A"))
finalSizedata.R5$inc_per<-factor(finalSizedata.R5$inc_per, levels = c("Baseline","ShorterIP"))
schooldaylost.R5<-data.frame(output=temp.nsdl,inc_per=temp.ip,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100)))
schooldaylost.R5$strategy<-factor(schooldaylost.R5$strategy, levels = c("SI","RS","RS_A"))
schooldaylost.R5$inc_per<-factor(schooldaylost.R5$inc_per, levels = c("Baseline","ShorterIP"))

R5fs.ch<-ggplot(data = finalSizedata.R5, aes(x=strategy, fill=inc_per,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="Incubation Period", labels=c(unname(TeX(c("Baseline", "ShorterIP")))) )+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

R5nsdl.ch<-ggplot(data = schooldaylost.R5, aes(x=strategy, fill=inc_per,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="Incubation Period", labels=c(unname(TeX(c("Baseline", "ShorterIP")))) )+stat_summary(fun.y=mean, geom="point", shape=15, size=3, color="royalblue", fill="black",position=position_dodge(0.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("NSDL"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

ggarrange(R5fs.ch,R5nsdl.ch,legend = "top",common.legend = TRUE)



#############################################
# R6 - varying t.delay
temp.fs<-NULL
temp.nsdl<-NULL
temp.td<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR6/sa6/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.td<-c(temp.td,rep("1",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR6/sa6/16/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.td<-c(temp.td,rep("1",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR6/sa6/26/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.td<-c(temp.td,rep("1",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR6/sa14/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.td<-c(temp.td,rep("2",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR6/sa14/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.td<-c(temp.td,rep("2",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR6/sa14/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.td<-c(temp.td,rep("2",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR6/sa15/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.td<-c(temp.td,rep("3",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR6/sa15/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.td<-c(temp.td,rep("3",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR6/sa15/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.td<-c(temp.td,rep("3",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))


finalSizedata.R6<-data.frame(output=temp.fs,t_delay=temp.td,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
finalSizedata.R6$strategy<-factor(finalSizedata.R6$strategy, levels = c("SI","RS","RS_A"))
finalSizedata.R6$t_delay<-factor(finalSizedata.R6$t_delay, levels = c("1","2","3"))
schooldaylost.R6<-data.frame(output=temp.nsdl,t_delay=temp.td,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
schooldaylost.R6$strategy<-factor(schooldaylost.R6$strategy, levels = c("SI","RS","RS_A"))
schooldaylost.R6$t_delay<-factor(schooldaylost.R6$t_delay, levels = c("1","2","3"))

R6fs<-ggplot(data = finalSizedata.R6, aes(x=strategy, fill=t_delay,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="Test Result Delay", labels=c(unname(TeX(c("$1$", "$2$","$3$")))) )+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

R6nsdl<-ggplot(data = schooldaylost.R6, aes(x=strategy, fill=t_delay,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="Test Result Delay", labels=c(unname(TeX(c("$1$", "$2$","$3$")))) )+stat_summary(fun.y=mean, geom="point", shape=15, size=3, color="royalblue", fill="black",position=position_dodge(0.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("NSDL"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

ggarrange(R6fs,R6nsdl,legend = "top",common.legend = TRUE)



#################################################################################
# R7 - Compliance
scaling.parm<-15
temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR7/sa18/1/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("0",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR7/sa9/20/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("0.2",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR7/sa9/22/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("0.4",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR7/sa9/24/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("0.6",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR7/sa9/26/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("0.8",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR7/sa9/28/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("1",length(epi.data$FinSize)))


data.R7<-data.frame("SummaryMeasure"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "Threshold"=rep(temp.thresh,2))
data.R7$Threshold<-factor(data.R7$Threshold, levels = c("0","0.2","0.4","0.6" ,"0.8","1"))

ggplot(data = data.R7, aes(x=Threshold, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15,17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),6), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Compliance")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+labs(fill="Model Outcome")+scale_x_discrete(labels = c("0","0.2","0.4","0.6","0.8","1"))

########################################################################3
# R8 - school size

temp.fs<-NULL
temp.nsdl<-NULL
temp.sz<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR8/5-8/sa16/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("200",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/5-8/sa16/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("200",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/5-8/sa16/8/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("200",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/5-8/sa16/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("1000",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/5-8/sa16/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("1000",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/5-8/sa16/10/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("1000",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/5-8/sa16/3/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("2000",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/5-8/sa16/7/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("2000",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/5-8/sa16/11/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("2000",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))


finalSizedata.R8<-data.frame(output=temp.fs,school_size=temp.sz,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
finalSizedata.R8$strategy<-factor(finalSizedata.R8$strategy, levels = c("SI","RS","RS_A"))
finalSizedata.R8$school_size<-factor(finalSizedata.R8$school_size, levels = c("200","1000","2000"))
schooldaylost.R8<-data.frame(output=temp.nsdl,school_size=temp.sz,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
schooldaylost.R8$strategy<-factor(schooldaylost.R8$strategy, levels = c("SI","RS","RS_A"))
schooldaylost.R8$school_size<-factor(schooldaylost.R8$school_size, levels = c("200","1000","2000"))

R8fs<-ggplot(data = finalSizedata.R8, aes(x=strategy, fill=school_size,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="School Size", labels=c(unname(TeX(c("$200$", "$1000$","$2000$")))) )+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

R8nsdl<-ggplot(data = schooldaylost.R8, aes(x=strategy, fill=school_size,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="School Size", labels=c(unname(TeX(c("$1$", "$2$","$3$")))) )+stat_summary(fun.y=mean, geom="point", shape=15, size=3, color="royalblue", fill="black",position=position_dodge(0.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("NSDL"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

ggarrange(R8fs,R8nsdl,legend = "top",common.legend = TRUE)

#school size with weighted seeds
temp.fs<-NULL
temp.nsdl<-NULL
temp.sz<-NULL
temp.str<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR8/1-8/sa16/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("200",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/1-8/sa16/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("200",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/1-8/sa16/8/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("200",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/5-8/sa16/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("1000",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/5-8/sa16/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("1000",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/5-8/sa16/10/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("1000",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/10-8/sa16/3/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("2000",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/10-8/sa16/7/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("2000",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR8/10-8/sa16/11/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.sz<-c(temp.sz,rep("2000",length(epi.data$FinSize)))
temp.str<-c(temp.str,rep("RS_A",length(epi.data$FinSize)))


finalSizedata.R8.2<-data.frame(output=temp.fs,school_size=temp.sz,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
finalSizedata.R8.2$strategy<-factor(finalSizedata.R8.2$strategy, levels = c("SI","RS","RS_A"))
finalSizedata.R8.2$school_size<-factor(finalSizedata.R8.2$school_size, levels = c("200","1000","2000"))
schooldaylost.R8.2<-data.frame(output=temp.nsdl,school_size=temp.sz,strategy=temp.str,gruop1=c(rep("1",100),rep("2",100),rep("3",100),rep("4",100),rep("5",100),rep("6",100),rep("7",100),rep("8",100),rep("9",100)))
schooldaylost.R8.2$strategy<-factor(schooldaylost.R8.2$strategy, levels = c("SI","RS","RS_A"))
schooldaylost.R8.2$school_size<-factor(schooldaylost.R8.2$school_size, levels = c("200","1000","2000"))

R8fs.2<-ggplot(data = finalSizedata.R8.2, aes(x=strategy, fill=school_size,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="School Size", labels=c(unname(TeX(c("$200$", "$1000$","$2000$")))) )+stat_summary(fun.y=mean, geom="point", shape=17, size=3, color="royalblue", fill="black",position=position_dodge(.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

R8nsdl.2<-ggplot(data = schooldaylost.R8.2, aes(x=strategy, fill=school_size,y=output))+geom_boxplot()+scale_fill_brewer(palette="Pastel1", name="School Size", labels=c(unname(TeX(c("$1$", "$2$","$3$")))) )+stat_summary(fun.y=mean, geom="point", shape=15, size=3, color="royalblue", fill="black",position=position_dodge(0.75), aes(group=gruop1)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("NSDL"))+xlab("Testing Strategy")+scale_x_discrete(labels=c("SI","ReaS","RepS"))

ggarrange(R8fs.2,R8nsdl.2,legend = "top",common.legend = TRUE)

#R10 - No class & school closures
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR10/sa6/8/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR10/sa6/18/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR10/sa6/28/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS_A",length(epi.data$FinSize)))

scaling.parm<-10
data.R10<-data.frame("Label"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.strategy))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "Strategy"=c(as.character(temp.strategy),as.character(temp.strategy)))

data.R10$Strategy<-factor(data.R10$Strategy, levels = c("SI","RS","RS_A"))

ggplot(data = data.R10, aes(x=Strategy, fill=Label,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),3), fill="black", aes(group=Label),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+scale_x_discrete(labels=c("SI","ReaS","RepS"))+labs(fill="Model Outcome")
#save 9:5.5



# R9 Omicron scenario: low pre-immunity, shorter incubation period, same R0
#1-test-per-week
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR9/sa17/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR9/sa17/14/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR9/sa17/22/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS_A",length(epi.data$FinSize)))

scaling.parm<-10
data.R9De1<-data.frame("Label"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.strategy))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "Strategy"=c(as.character(temp.strategy),as.character(temp.strategy)))

temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR9/sa17/1/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR9/sa17/9/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR9/sa17/17/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS_A",length(epi.data$FinSize)))

scaling.parm2<-10
data.R9Om1<-data.frame("Label"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.strategy))), "Value"=c(temp.fs,temp.nsdl/scaling.parm2), "Strategy"=c(as.character(temp.strategy),as.character(temp.strategy)))

data.overallR9.1<-data.frame("Scenario"=c(rep("Delta",length(data.R9De1$Label)),rep("Omicron",length(data.R9Om1$Label))), "SummaryMeasure"=c(as.character(data.R9De1$Label), as.character(data.R9Om1$Label)), "Value"=c(data.R9De1$Value,data.R9Om1$Value), "Strategy"=c(as.character(data.R9De1$Strategy),as.character(data.R9Om1$Strategy))     )
data.overallR9.1$Scenario<-factor(data.overallR9.1$Scenario, levels = c("Delta","Omicron"))
data.overallR9.1$Strategy<-factor(data.overallR9.1$Strategy, levels = c("SI","RS","RS_A"))

ggplot(data = data.overallR9.1, aes(x=Strategy, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15,17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),6), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+facet_wrap( ~ Scenario, scales = "fixed")+scale_x_discrete(labels=c("SI","ReaS","RepS"))+labs(fill="Model Outcome")
#save 9:5.5

#1 vs 2 weekly test
scaling.parm<-10
temp.fs<-NULL
temp.nsdl<-NULL
temp.test<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR9/sa17/17/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.test<-c(temp.test,rep("1",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR9/sa17/41/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.test<-c(temp.test,rep("2",length(epi.data$FinSize)))

data.R9NTests<-data.frame("SummaryMeasure"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "NTest"=rep(temp.test,2))

ggplot(data = data.R9NTests, aes(x=NTest, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15), size=3, color=rep(c("royalblue","gold"),2), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Number of Tests Per Week")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+labs(fill="Model Outcome")



# R11 - No testing scenario
# with school closure
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR11/sa18/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("Bs",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR11/base/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR11/base/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR11/base/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS_A",length(epi.data$FinSize)))

scaling.parm<-300
data.R11<-data.frame("Label"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.strategy))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "Strategy"=c(as.character(temp.strategy),as.character(temp.strategy)))

data.R11$Strategy<-factor(data.R11$Strategy, levels = c("Bs","SI","RS","RS_A"))

ggplot(data = data.R11, aes(x=Strategy, fill=Label,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),4), fill="black", aes(group=Label),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+scale_x_discrete(labels=c("No Testing","SI","ReaS","RepS"))+labs(fill="Model Outcome")
#save 9:5.5

#without school closure
temp.fs<-NULL
temp.nsdl<-NULL
temp.strategy<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR11/sa18/1/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("Bs",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR11/sa6/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("SI",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR11/sa6/16/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR11/sa6/26/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.strategy<-c(temp.strategy,rep("RS_A",length(epi.data$FinSize)))

scaling.parm<-10
data.R11<-data.frame("Label"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.strategy))), "Value"=c(temp.fs,temp.nsdl/scaling.parm), "Strategy"=c(as.character(temp.strategy),as.character(temp.strategy)))

data.R11$Strategy<-factor(data.R11$Strategy, levels = c("Bs","SI","RS","RS_A"))

ggplot(data = data.R11, aes(x=Strategy, fill=Label,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),4), fill="black", aes(group=Label),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Testing Strategy")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+scale_x_discrete(labels=c("No Testing","SI","ReaS","RepS"))+labs(fill="Model Outcome")
#save 9:5.5

#Class threshold SI
scaling.parm<-100

temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR12/sa6/0/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("2",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR12/sa6/2/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("4",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR12/sa6/4/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("6",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR12/sa6/6/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("8",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR12/sa6/8/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("1000",length(epi.data$FinSize)))

data.SIcl<-data.frame("SummaryMeasure"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs ,temp.nsdl/scaling.parm), "Threshold"=c(as.character(temp.thresh),as.character(temp.thresh))   )
data.SIcl$Threshold<-factor(data.SIcl$Threshold, levels = c("2","4","6","8","1000"))
ggplot(data = data.SIcl, aes(x=Threshold, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),5), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Class Closure Threshold")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+labs(fill="Model Outcome")+scale_x_discrete(labels = c("2","4","6","8","No Threshold"))


scaling.parm<-100

temp.fs<-NULL
temp.nsdl<-NULL
temp.thresh<-NULL
epi.data<-read.csv("~/Desktop/Testing/FigureR12/sa6/10/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("2",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR12/sa6/12/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("4",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR12/sa6/14/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("6",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR12/sa6/16/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("8",length(epi.data$FinSize)))
epi.data<-read.csv("~/Desktop/Testing/FigureR12/sa6/18/out.csv")
temp.fs<-c(temp.fs,epi.data$FinSizeProp)
temp.nsdl<-c(temp.nsdl,epi.data$n.schooldayslost)
temp.thresh<-c(temp.thresh,rep("1000",length(epi.data$FinSize)))

data.RScl<-data.frame("SummaryMeasure"=c(rep("FinalSize",length(temp.fs)),rep("SchoolDaysLost",length(temp.nsdl))), "Value"=c(temp.fs ,temp.nsdl/scaling.parm), "Threshold"=c(as.character(temp.thresh),as.character(temp.thresh))   )
data.RScl$Threshold<-factor(data.SIcl$Threshold, levels = c("2","4","6","8","1000"))
ggplot(data = data.RScl, aes(x=Threshold, fill=SummaryMeasure,y=Value))+geom_boxplot()+scale_fill_brewer(palette="Dark2", labels=c(" Attack Rate","NSDL") )+stat_summary(fun.y=mean, geom="point", shape=c(17,15,17,15,17,15,17,15,17,15), size=3, color=rep(c("royalblue","gold"),5), fill="black", aes(group=SummaryMeasure),position=position_dodge(.75)) +theme(
  panel.background = element_rect(fill = "white",
                                  colour = "gray",
                                  size = 0.75, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                  colour = "gray"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"),
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1,"cm"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(size=14),
  #legend.position = "top",
  #legend.position = "top",
  legend.title = element_text(size=15, face = "bold"),
  axis.text.x = element_text(size=14, angle = 45, vjust = 0.5, hjust = 0.5),
  axis.title = element_text(size=15),
  axis.text.y = element_text(size=14),
  strip.text = element_text(size = 15)
)+ylab(("Attack Rate"))+xlab("Class Closure Threshold")+scale_y_continuous(sec.axis = sec_axis(~.*scaling.parm, name = "NSDL"))+labs(fill="Model Outcome")+scale_x_discrete(labels = c("2","4","6","8","No Threshold"))


