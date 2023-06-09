#Code to analyze data from incubation experiment

setwd("~/UV/2016Incubations")

#packages needed
library("dplyr")
library("tidyr")
library("ggplot2")
library(lme4)
library(cowplot)
library(viridis)

#DATA320<- read.table("LightExperimentInfections.txt", header=TRUE, sep="\\t", na.strings="", dec=".", strip.white=TRUE, stringsAsFactors = FALSE) 
DATA<- read.table("LightExperimentData.txt", header=TRUE, sep="\\t", na.strings="", dec=".", strip.white=TRUE, stringsAsFactors = FALSE) 

DATA320<-select(DATA, Lake, Month, Vial, Julian, Para, LgtTrt, depth, inf, Prop, Ninf)
DATA320$Month<-as.factor(DATA320$Month)
DATA320$Month <- factor(DATA320$Month,levels(DATA320$Month)[c(2,1,3)])
DATA320$Para<-as.factor(DATA320$Para)
DATA320$Para <- factor(DATA320$Para,levels(DATA320$Para)[c(2,1)])

#for labelling
dat_text <- data.frame(
  label = c("A", "B", "C","D","E","F"),
  Para   = c("Pasteuria", "Metschnikowia", "Pasteuria", "Metschnikowia","Pasteuria", "Metschnikowia"),
  Month = c("July","July","August","August","November","November")
)

#Produce Figure 1A
A<-ggplot(DATA320, aes(x=depth, y=Prop))+
  facet_grid(Month~Para, scales="free_y")+
  geom_boxplot(outlier.color="black",outlier.size=1,aes(fill=LgtTrt))+
  ylab("Proportion infected")+
  scale_fill_manual(values=c("darkgrey", "white"))+ 
  scale_color_manual(values=c("darkgrey", "white"))+
  xlab("Depth of incubation")+theme_bw()+
  theme(panel.border = element_blank())+
  theme(legend.title = element_blank())+
  theme(strip.background =element_rect(fill="white", color="white"))+
  theme(strip.text.x = element_text(face="italic", size=14))+
  theme(strip.text.y=element_text(size=14))+
  theme(legend.text = element_text(size=14))+
  theme(axis.title.x = element_text(size=14))+
  theme(axis.title.y = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(size=12, color="Black"))+
  theme(legend.justification="top")+
  theme(axis.line = element_blank())+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  geom_text(data=dat_text, mapping=aes(x=0.5,y=1,label=label))

#Make parts B and C of figure
LL.Avg<-DATA%>%group_by(Month)%>%summarise(MeanPAR=mean(real.PAR), sdPAR=sd(real.PAR), Mean320=mean(real.320), sd320=sd(real.320))
LL.Avg$Month <- factor(LLA$Month, levels = c("July", "August", "November"))

B<-ggplot(LL.Avg, aes(x=Month, y=MeanPAR))+
  geom_errorbar(aes(ymin=MeanPAR-sdPAR, ymax=MeanPAR+sdPAR), width=0)+
  geom_point(size=2)+
  ylab(expression(paste("Cumulative PAR (W/m"^"2",")")))+
  theme_classic()+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size=12))+
  scale_x_discrete(labels=c("July","Aug","Nov"))+
  scale_y_continuous(breaks=c(20000000,30000000,40000000), labels=c("2","3","4"))+
  theme(axis.text = element_text(size=12, color="Black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate("text", label = c("H"), size = 4, x = 0.6, y = 39000000)+
  theme(axis.line = element_line(size=0.5))

C<-ggplot(LL.Avg, aes(x=Month, y=Mean320))+
  geom_errorbar(aes(ymin=Mean320-sd320, ymax=Mean320+sd320), width=0)+
  geom_point(size=2)+
  ylab(expression(paste("Cumulative UV (320nm, KJ/m"^"2",")")))+
  theme_classic()+
  theme(axis.title.y = element_text(size=14))+
  theme(axis.title.x=element_blank())+
  theme(axis.text = element_text(size=12, color="Black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_blank())+
  annotate("text", label = c("G"), size = 4, x = 0.6, y = 31)+
  theme(axis.line = element_line(size=0.5))

BC<-plot_grid(C, B, nrow = 2, align = "v")

fig1<-ggdraw() +
  draw_plot(A, 0, 0, 0.75, 1) +
  draw_plot(BC, 0.7, 0.03, 0.25, 0.7)

save_plot("LightFigure1.jpg", fig1, base_width = 8, base_height = 7.5)
#############################################

#Analyze at Pasteuria and Metschnikowia separately
Pprops<-filter(DATA320, Para=="Pasteuria")
Mprops<-filter(DATA320, Para=="Metschnikowia")

#Statistics
fullglmm<-glmer(cbind(DATA320$inf, DATA320$Ninf)~Para*LgtTrt+(1|Lake:Month), family=binomial, data=DATA320)
summary(fullglmm)
#This model shows a significant interaction between light interaction (P=0.0048) and parasite 
#Indicating that Metsch is more sensitive to light

#These models are for pasteuria
PastGLMM<-glmer(cbind(Pprops$inf, Pprops$Ninf)~LgtTrt*depth*Month+(1|Lake:Month), 
                family=binomial, data=Pprops)
summary(PastGLMM)
#Fails to converge

PastGLMM<-glmer(cbind(Pprops$inf, Pprops$Ninf)~LgtTrt*depth+LgtTrt*Month+depth*Month+(1|Lake:Month), 
                family=binomial, data=Pprops)
summary(PastGLMM)

#Pasteuria By Month
#July
July<-filter(Pprops, Month=="July")
JulyGLMM<-glmer(cbind(July$inf, July$Ninf)~LgtTrt*depth+(1|Lake), 
                family=binomial, data=July)
summary(JulyGLMM)
#significant effect of light and interaction

#August
August<-filter(Pprops, Month=="August")
AugustGLMM<-glmer(cbind(August$inf, August$Ninf)~LgtTrt*depth+(1|Lake), 
                  family=binomial, data=August)
summary(AugustGLMM)
#significant effect of light and marginal effect of interaction

#November
November<-filter(Pprops, Month=="November")
NovemberGLMM<-glmer(cbind(November$inf, November$Ninf)~LgtTrt*depth+(1|Lake), 
                    family=binomial, data=November)
summary(NovemberGLMM)
#No effect of light or interaction.


#Metschnikowia

#Glmm -- used drop1(test="Chisq") to narrow down to best model
MetschGLMM<-glmer(cbind(Mprops$inf, Mprops$Ninf)~LgtTrt*depth*Month+(1|Lake), 
                  family=binomial, data=Mprops)
#This model fails to converge
MetschGLMM<-glmer(cbind(Mprops$inf, Mprops$Ninf)~LgtTrt*depth+LgtTrt*Month+
                    depth*Month+(1|Lake:Month), 
                  family=binomial, data=Mprops)
MetschGLMM<-glmer(cbind(Mprops$inf, Mprops$Ninf)~LgtTrt*depth+
                    LgtTrt*Month+(1|Lake:Month), 
                  family=binomial, data=Mprops)
MetschGLMM<-glmer(cbind(Mprops$inf, Mprops$Ninf)~LgtTrt*depth+Month+(1|Lake:Month), 
                  family=binomial, data=Mprops)#This is the model used in the paper
summary(MetschGLMM)

#Look at each month separately
MJuly<-filter(Mprops, Month=="July")
MAugust<-filter(Mprops, Month=="August")
MNov<-filter(Mprops, Month=="November")

#JulyGLMM
JulyM<-glmer(cbind(MJuly$inf, MJuly$Ninf)~LgtTrt*depth+(1|Lake), 
             family=binomial, data=MJuly)
JulyM<-glmer(cbind(MJuly$inf, MJuly$Ninf)~LgtTrt+depth+(1|Lake), 
             family=binomial, data=MJuly)
JulyM<-glmer(cbind(MJuly$inf, MJuly$Ninf)~LgtTrt+(1|Lake), 
             family=binomial, data=MJuly) #This is the model in the paper
summary(JulyM)

#AugustGLMM
AugustM<-glmer(cbind(MAugust$inf, MAugust$Ninf)~LgtTrt*depth+(1|Lake), 
               family=binomial, data=MAugust)
AugustM<-glmer(cbind(MAugust$inf, MAugust$Ninf)~LgtTrt+depth+(1|Lake), 
               family=binomial, data=MAugust)
AugustM<-glmer(cbind(MAugust$inf, MAugust$Ninf)~LgtTrt+(1|Lake), 
               family=binomial, data=MAugust)#This is the model in the paper
summary(AugustM)

#GLMM
NovemberM<-glmer(cbind(MNov$inf, MNov$Ninf)~LgtTrt+(1|Lake), 
                 family=binomial, data=MNov)
summary(NovemberM)


###For supplement S3.
small<-select(DATA, Lake, Month, Para, LgtTrt, depth.number, Prop)
small.L<-filter(small, LgtTrt=="Light")
small.D<-filter(small, LgtTrt=="Dark")

LightStuff<-select(DATA, Lake, Month, depth.number, real.PAR, real.320, estL, estP)
LightStuff<-unique(LightStuff[,c("Lake","Month","depth.number", "real.PAR","real.320","estL","estP")])
Month<-c("July","July","August","August","November","November")
Para<-c("Pasteuria", "Metschnikowia","Pasteuria", "Metschnikowia","Pasteuria", "Metschnikowia")
Doses<-c(2000, 100, 2000, 250, 2000, 250)

DoseData<-as.data.frame(cbind(Month, Para, Doses))
L.Doses<-full_join(small.L, DoseData, by=c("Month", "Para"))
D.Doses<-full_join(small.D, DoseData, by=c("Month","Para"))
L.Doses$Doses<-as.numeric(as.character(L.Doses$Doses))
D.Doses$Doses<-as.numeric(as.character(D.Doses$Doses))

L.Doses$DosesL<-L.Doses$Doses*1000
D.Doses$DosesL<-D.Doses$Doses*1000

L.Doses$Prop[L.Doses$Prop==1]<-0.99
D.Doses$Prop[D.Doses$Prop==1]<-0.99
D.Doses<-D.Doses[-153,]
L.Doses$betaL<--log(1-L.Doses$Prop)/L.Doses$DosesL
D.Doses$betaD<--log(1-D.Doses$Prop)/D.Doses$DosesL

AvgDark<-D.Doses%>%group_by(Month, Lake, Para, depth.number)%>%summarize(AvgDarkBeta=mean(betaD))
Both<-full_join(L.Doses, AvgDark, by=c("Month","Lake","Para","depth.number"))
Both$StandardizedBeta<-Both$betaL/Both$AvgDarkBeta
Both.1<-filter(Both, !is.na(Both$StandardizedBeta))
Both.2<-filter(Both.1, StandardizedBeta!=Inf)
Fullgroup<-Both.2%>%group_by(Month, Lake,Para,depth.number)%>%summarize(AvgBeta=mean(StandardizedBeta), sdBeta=sd(StandardizedBeta))
Fullgroup.1<-full_join(Fullgroup, LightStuff, by=c("Lake","Month","depth.number"))

Fullgroup.1$Act.320<-(Fullgroup.1$estL/100)*Fullgroup.1$real.320
Fullgroup.1$Act.par<-(Fullgroup.1$estP/100)*Fullgroup.1$real.PAR

Fullgroup.1$Month <- factor(Fullgroup.1$Month, levels = c("July", "August", "November"))
Fullgroup.1$Para <- factor(Fullgroup.1$Para, levels = c("Pasteuria", "Metschnikowia"))
Fullgroup.1$depth.number<-as.factor(Fullgroup.1$depth.number)
forgraph<-select(Fullgroup.1, Month, Lake, Para, depth.number, AvgBeta, sdBeta, Act.320, Act.par)
forgraph.1<-gather(forgraph, key="Light",value="amount", Act.320, Act.par)
forgraph.1$Month <- factor(forgraph.1$Month, levels = c("July", "August", "November"))
forgraph.1$Para <- factor(forgraph.1$Para, levels = c("Pasteuria", "Metschnikowia"))

forgraph.1$Light[forgraph.1$Light=="Act.320"]<-"320nm UV exposure (KJ/m^2)"
forgraph.1$Light[forgraph.1$Light=="Act.par"]<-"PAR exposure (W/m^2)"


UV<-filter(forgraph.1, Light=="320nm UV exposure (KJ/m^2)")
PAR<-filter(forgraph.1, Light=="PAR exposure (W/m^2)")

dat_text2 <- data.frame(
  label = c("A","A", "C","C", "A","A","C","C", "A","A", "C","C"),
  Para = c("Pasteuria","Pasteuria","Metschnikowia","Metschnikowia", "Pasteuria","Pasteuria","Metschnikowia","Metschnikowia", "Pasteuria","Pasteuria","Metschnikowia","Metschnikowia"),
  Month = c("July","August","November", "July","August","November", "July","August","November", "July","August","November"),
  depth.number = c("0.5","2","0.5","2","0.5","2","0.5","2","0.5","2","0.5","2"))

U<-ggplot(UV, aes(amount, AvgBeta, color=Month, shape=depth.number))+  
  geom_errorbar(aes(ymin=AvgBeta-sdBeta, ymax=AvgBeta+sdBeta), width=0)+
  geom_point(size=4, fill="white", stroke=2)+
  scale_shape_manual(values=c(21,16), name="Incubation depth", labels=c("0.5 m", "2 m"))+
  facet_grid(Para~., scales="free")+
  ylab("Relative infectivity")+
  theme_classic()+
  scale_color_manual(values=c("#D55E00","#0072B2", "#CC79A7"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text.y = element_blank())+
  theme(strip.background=element_rect(color="white"))+
  theme(strip.text.x=element_text(size=14))+
  theme(legend.title = element_blank())+
  theme(legend.key = element_blank())+
  theme(legend.position="none")+
  theme(axis.title.x = element_text(size=14))+
  labs(x=expression("320nm"~UV~exposure~(KJ/m^2)))+
  theme(axis.title.y = element_text(size=14))+
  theme(legend.text = element_blank())+
  theme(axis.text = element_text(size=12, color="Black"))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  annotate("segment", x=-Inf, xend=Inf, y=1, yend=1, lty="dashed")+
  scale_x_continuous(breaks=c(0,5,10), labels=c("0","5","10"))+
  geom_text(data=dat_text2, mapping=aes(x=-0.4,y=3,label=label), size=4, color="black")

dat_text2 <- data.frame(
  label = c("B","B", "D","D", "B","B","D","D", "B","B", "D","D"),
  Para = c("Pasteuria","Pasteuria","Metschnikowia","Metschnikowia", "Pasteuria","Pasteuria","Metschnikowia","Metschnikowia", "Pasteuria","Pasteuria","Metschnikowia","Metschnikowia"),
  Month = c("July","August","November", "July","August","November", "July","August","November", "July","August","November"),
  depth.number = c("0.5","2","0.5","2","0.5","2","0.5","2","0.5","2","0.5","2"))

P<-ggplot(PAR, aes(amount, AvgBeta, color=Month, shape=depth.number))+  
  geom_errorbar(aes(ymin=AvgBeta-sdBeta, ymax=AvgBeta+sdBeta), width=0)+
  geom_point(size=4, fill="white", stroke=2)+
  scale_shape_manual(values=c(21,16), name="Incubation depth", labels=c("0.5 m", "2 m"))+
  facet_grid(Para~., scales="free")+
  theme_classic()+
  scale_color_manual(values=c("#D55E00","#0072B2", "#CC79A7"), name="Month", labels=c("July", "August", "November"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text.y = element_text(size=14))+
  theme(strip.background=element_rect(color="white"))+
  theme(strip.text.x=element_text(size=14))+
  theme(legend.title = element_text(size=14))+
  theme(axis.title.x = element_text(size=14))+
  labs(x=expression(PAR~exposure~(W/m^2~x10^7)))+
  theme(axis.title.y = element_blank())+
  theme(legend.text = element_text(size=12))+
  theme(axis.text = element_text(size=12, color="Black"))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  annotate("segment", x=-Inf, xend=Inf, y=1, yend=1, lty="dashed")+
  scale_x_continuous(breaks=c(10000000,20000000,30000000,40000000), labels=c("1","2","3","4"))+
  geom_text(data=dat_text2, mapping=aes(x=-100,y=3,label=label), size=4, color="black")


legend <- get_legend(P)
k<-plot_grid(U, P+theme(legend.position="none"), nrow = 1, align = "h")
k2<-plot_grid(k,legend,rel_widths=c(2,0.5))

save_plot("LightFigure2.jpg", k2, base_width = 9, base_height = 6)

UV<-filter(forgraph.1, Light=="320nm UV exposure (KJ/m^2)")
UV.P<-filter(UV, Para=="Pasteuria")
UV.M<-filter(UV, Para=="Metschnikowia")

PAR<-filter(forgraph.1, Light=="PAR exposure (W/m^2)")
PAR.P<-filter(PAR, Para=="Pasteuria")
PAR.M<-filter(PAR, Para=="Metschnikowia")

P.UV<-lm(AvgBeta~amount, data=UV.P)
drop1(P.UV, test="F")
M.UV<-lm(AvgBeta~amount, data=UV.M)
drop1(M.UV, test="F")

P.PAR<-lm(AvgBeta~amount, data=PAR.P)
drop1(P.PAR, test="F")
M.PAR<-lm(AvgBeta~amount, data=PAR.M)
drop1(M.PAR, test="F")
