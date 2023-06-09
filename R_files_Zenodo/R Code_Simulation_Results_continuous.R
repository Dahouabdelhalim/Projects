######################################################################################
### Results--single continuous outcome
######################################################################################

##Type I error probabilities:
TypeI_n50=round(readRDS("TypeI_n50.Rds"),digits=4)
TypeI_n75=round(readRDS("TypeI_n75.Rds"),digits=4)
TypeI_n100=round(readRDS("TypeI_n100.Rds"),digits=4)
TypeI_n125=round(readRDS("TypeI_n125.Rds"),digits=4)
TypeI_n150=round(readRDS("TypeI_n150.Rds"),digits=4)

TypeI_single_continuous_outcome=rbind(TypeI_n50,TypeI_n75,TypeI_n100,TypeI_n125,TypeI_n150)
colnames(TypeI_single_continuous_outcome)=c('S=1,T=1','S=2,T=1')
#write.csv(TypeI_single_continuous_outcome,'TypeI_single_continuous_outcome.csv',row.names = TRUE)

##Power:
Power_n50=round(readRDS("Power_n50.Rds"),digits=4)
Power_n75=round(readRDS("Power_n75.Rds"),digits=4)
Power_n100=round(readRDS("Power_n100.Rds"),digits=4)
Power_n125=round(readRDS("Power_n125.Rds"),digits=4)
Power_n150=round(readRDS("Power_n150.Rds"),digits=4)

Power_single_continuous_outcome=rbind(Power_n50,Power_n75,Power_n100,Power_n125,Power_n150)
colnames(Power_single_continuous_outcome)=c('S=1,T=1','S=2,T=1')
#write.csv(Power_single_continuous_outcome,'Power_single_continuous_outcome.csv',row.names = TRUE)

######################################################################################
### Results--multiple continuous outcomes
######################################################################################

#n0=n1=50
ANCOVA_n50_1=readRDS(file="ANCOVA_n50_1.Rds")
ANCOVA_n50_2=readRDS(file="ANCOVA_n50_2.Rds")
ANCOVA_n50_3=readRDS(file="ANCOVA_n50_3.Rds")

TypeI_multiple_continuous_n50=rbind(ANCOVA_n50_1[c(1,3,5),],ANCOVA_n50_2[c(1,3,5),],ANCOVA_n50_3[c(1,3,5),])
TypeI_multiple_continuous_n50=round(TypeI_multiple_continuous_n50,digits=4)
colnames(TypeI_multiple_continuous_n50)=1:9
#write.csv(TypeI_multiple_continuous_n50,'TypeI_multiple_continuous_n50.csv',row.names = TRUE)

#n0=n1=100
ANCOVA_n100_1=readRDS(file="ANCOVA_n100_1.Rds")
ANCOVA_n100_2=readRDS(file="ANCOVA_n100_2.Rds")
ANCOVA_n100_3=readRDS(file="ANCOVA_n100_3.Rds")

TypeI_multiple_continuous_n100=rbind(ANCOVA_n100_1[c(1,3,5),],
                                 ANCOVA_n100_2[c(1,3,5),],
                                 ANCOVA_n100_3[c(1,3,5),])
TypeI_multiple_continuous_n100=round(TypeI_multiple_continuous_n100,digits=4)
colnames(TypeI_multiple_continuous_n100)=1:9
#write.csv(TypeI_multiple_continuous_n100,'TypeI_multiple_continuous_n100.csv',row.names = TRUE)

#n0=n1=150
ANCOVA_n150_1=readRDS(file="ANCOVA_n150_1.Rds")
ANCOVA_n150_2=readRDS(file="ANCOVA_n150_2.Rds")
ANCOVA_n150_3=readRDS(file="ANCOVA_n150_3.Rds")

TypeI_multiple_continuous_n150=rbind(ANCOVA_n150_1[c(1,3,5),],
                                 ANCOVA_n150_2[c(1,3,5),],
                                 ANCOVA_n150_3[c(1,3,5),])
TypeI_multiple_continuous_n150=round(TypeI_multiple_continuous_n150,digits=4)
colnames(TypeI_multiple_continuous_n150)=1:9
#write.csv(TypeI_multiple_continuous_n150,'TypeI_multiple_continuous_n150.csv',row.names = TRUE)

### Power plots--multiple continuous outcomes
library(ggplot2)
library(RColorBrewer)
library(ggpattern)
library(ggpubr)
library(splines)

#power plot: n=50, 100, 150
power_multiple_continuous_ANCOVA_n50=rbind(ANCOVA_n50_1[c(2,4,6),],ANCOVA_n50_2[c(2,4,6),],ANCOVA_n50_3[c(2,4,6),])
colnames(power_multiple_continuous_ANCOVA_n50)=1:9

power_multiple_continuous_ANCOVA_n100=rbind(ANCOVA_n100_1[c(2,4,6),],ANCOVA_n100_2[c(2,4,6),],ANCOVA_n100_3[c(2,4,6),])
colnames(power_multiple_continuous_ANCOVA_n100)=1:9

power_multiple_continuous_ANCOVA_n150=rbind(ANCOVA_n150_1[c(2,4,6),],ANCOVA_n150_2[c(2,4,6),],ANCOVA_n150_3[c(2,4,6),])
colnames(power_multiple_continuous_ANCOVA_n150)=1:9

#maximum power of each row
apply(power_multiple_continuous_ANCOVA_n50, 1, max)
#Sopt_ANCOVA
Sopt_ANCOVA_n50=rep(NA,9)
for (i in 1:9){
  Sopt_ANCOVA_n50[i]=which(power_multiple_continuous_ANCOVA_n50[i,]%in%max(power_multiple_continuous_ANCOVA_n50[i,]))
}
Sopt_ANCOVA_n50 #5 5 3 5 5 5 3 5 5

apply(power_multiple_continuous_ANCOVA_n100, 1, max)
#Sopt_ANCOVA
Sopt_ANCOVA_n100=rep(NA,9)
for (i in 1:9){
  Sopt_ANCOVA_n100[i]=which(power_multiple_continuous_ANCOVA_n100[i,]%in%max(power_multiple_continuous_ANCOVA_n100[i,]))
}
Sopt_ANCOVA_n100 #4 4 4 2 4 4 4 4 4

apply(power_multiple_continuous_ANCOVA_n150, 1, max)
#Sopt_ANCOVA
Sopt_ANCOVA_n150=rep(NA,9)
for (i in 1:9){
  Sopt_ANCOVA_n150[i]=which(power_multiple_continuous_ANCOVA_n150[i,]%in%max(power_multiple_continuous_ANCOVA_n150[i,]))
}
Sopt_ANCOVA_n150 #4 5 4 2 4 5 4 4 4

Sopt_ANCOVA=rbind(Sopt_ANCOVA_n50,Sopt_ANCOVA_n100,Sopt_ANCOVA_n150)

#highlight Sopt_ANCOVA
#display.brewer.all()
cols1 <- c(brewer.pal(9, "Pastel1")[1],brewer.pal(9, "Paired")[1],brewer.pal(9, "Paired")[2])

#9 plots: rhoXY= 0.5,rhoX= 0.6; rhoXY= 0.5,rhoX= 0.7; rhoXY= 0.5,rhoX= 0.8; rhoXY= 0.5,rhoX= 0.9;
#rhoXY= 0.6,rhoX= 0.7; rhoXY= 0.6,rhoX= 0.8; rhoXY= 0.6,rhoX= 0.9;
#rhoXY= 0.7,rhoX= 0.8; rhoXY= 0.7,rhoX= 0.9

#rhoXY= 0.5,rhoX= 0.6
power_multiple_continuous_ANCOVA_1 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_continuous_ANCOVA_n50[1,],
                                              power_multiple_continuous_ANCOVA_n100[1,],
                                              power_multiple_continuous_ANCOVA_n150[1,]))
power_multiple_continuous_ANCOVA_1$n=as.character(power_multiple_continuous_ANCOVA_1$n)
power_multiple_continuous_ANCOVA_1$n=factor(power_multiple_continuous_ANCOVA_1$n, levels=unique(power_multiple_continuous_ANCOVA_1$n))

power_multiple_continuous_ANCOVA_1$S=as.character(power_multiple_continuous_ANCOVA_1$S)
power_multiple_continuous_ANCOVA_1$S=factor(power_multiple_continuous_ANCOVA_1$S, levels=unique(power_multiple_continuous_ANCOVA_1$S))
power_multiple_continuous_ANCOVA_1

power_multiple_continuous_ANCOVA_n50_1=power_multiple_continuous_ANCOVA_1[power_multiple_continuous_ANCOVA_1$n=='n=50',][Sopt_ANCOVA['Sopt_ANCOVA_n50',1],]
power_multiple_continuous_ANCOVA_n100_1=power_multiple_continuous_ANCOVA_1[power_multiple_continuous_ANCOVA_1$n=='n=100',][Sopt_ANCOVA['Sopt_ANCOVA_n100',1],]
power_multiple_continuous_ANCOVA_n150_1=power_multiple_continuous_ANCOVA_1[power_multiple_continuous_ANCOVA_1$n=='n=150',][Sopt_ANCOVA['Sopt_ANCOVA_n150',1],]

Sopt_ANCOVA_multiple_continuous_1=rbind(power_multiple_continuous_ANCOVA_n50_1,
                             power_multiple_continuous_ANCOVA_n100_1,
                             power_multiple_continuous_ANCOVA_n150_1)

Plot_power_multiple_continuous_ANCOVA_1<-ggplot(power_multiple_continuous_ANCOVA_1, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.5~','~rho[X]==rho[Y] ~'=0.6'))+  
  scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_ANCOVA_multiple_continuous_1, 
             aes(x=S,y=power), color=brewer.pal(9, "Paired")[6],size=7)+
  theme(text = element_text(size = 45),
        axis.title=element_text(size=40),
        axis.text=element_text(size=35),
        #axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_rect(fill = "white", linetype = "solid",color = "black"),
        # Add axis line
        axis.line = element_line(color = "black"),
        legend.position = 'top',
        legend.title=element_text(size =38),
        legend.text =element_text(size =38),plot.margin=unit(c(1,1,-0.1,1), "cm"))+
  guides(color = guide_legend(title ="Sample sizes per group",nrow = 1,byrow=TRUE)) #bquote(~n[0]==n[1]~'=')
Plot_power_multiple_continuous_ANCOVA_1

#rhoXY= 0.5,rhoX= 0.7
power_multiple_continuous_ANCOVA_2 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_continuous_ANCOVA_n50[2,],
                                              power_multiple_continuous_ANCOVA_n100[2,],
                                              power_multiple_continuous_ANCOVA_n150[2,]))
power_multiple_continuous_ANCOVA_2$n=as.character(power_multiple_continuous_ANCOVA_2$n)
power_multiple_continuous_ANCOVA_2$n=factor(power_multiple_continuous_ANCOVA_2$n, levels=unique(power_multiple_continuous_ANCOVA_2$n))

power_multiple_continuous_ANCOVA_2$S=as.character(power_multiple_continuous_ANCOVA_2$S)
power_multiple_continuous_ANCOVA_2$S=factor(power_multiple_continuous_ANCOVA_2$S, levels=unique(power_multiple_continuous_ANCOVA_2$S))
power_multiple_continuous_ANCOVA_2

power_multiple_continuous_ANCOVA_n50_2=power_multiple_continuous_ANCOVA_2[power_multiple_continuous_ANCOVA_2$n=='n=50',][Sopt_ANCOVA['Sopt_ANCOVA_n50',2],]
power_multiple_continuous_ANCOVA_n100_2=power_multiple_continuous_ANCOVA_2[power_multiple_continuous_ANCOVA_2$n=='n=100',][Sopt_ANCOVA['Sopt_ANCOVA_n100',2],]
power_multiple_continuous_ANCOVA_n150_2=power_multiple_continuous_ANCOVA_2[power_multiple_continuous_ANCOVA_2$n=='n=150',][Sopt_ANCOVA['Sopt_ANCOVA_n150',2],]

Sopt_ANCOVA_multiple_continuous_2=rbind(power_multiple_continuous_ANCOVA_n50_2,
                             power_multiple_continuous_ANCOVA_n100_2,
                             power_multiple_continuous_ANCOVA_n150_2)

Plot_power_multiple_continuous_ANCOVA_2<-ggplot(power_multiple_continuous_ANCOVA_2, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.5~','~rho[X]==rho[Y] ~'=0.7'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_ANCOVA_multiple_continuous_2, 
             aes(x=S,y=power), color=brewer.pal(9, "Paired")[6],size=7)+
  theme(text = element_text(size = 45),
        axis.title=element_text(size=40),
        axis.text=element_text(size=35),
        #axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_rect(fill = "white", linetype = "solid",color = "black"),
        # Add axis line
        axis.line = element_line(color = "black"),
        legend.position = 'top',
        legend.title=element_text(size =38),
        legend.text =element_text(size =38),plot.margin=unit(c(1,2,-0.1,2), "cm"))+
  guides(color = guide_legend(title ="Sample sizes per group",nrow = 1,byrow=TRUE)) #bquote(~n[0]==n[1]~'=')
Plot_power_multiple_continuous_ANCOVA_2

#rhoXY= 0.5,rhoX= 0.8
power_multiple_continuous_ANCOVA_3 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_continuous_ANCOVA_n50[3,],
                                              power_multiple_continuous_ANCOVA_n100[3,],
                                              power_multiple_continuous_ANCOVA_n150[3,]))
power_multiple_continuous_ANCOVA_3$n=as.character(power_multiple_continuous_ANCOVA_3$n)
power_multiple_continuous_ANCOVA_3$n=factor(power_multiple_continuous_ANCOVA_3$n, levels=unique(power_multiple_continuous_ANCOVA_3$n))

power_multiple_continuous_ANCOVA_3$S=as.character(power_multiple_continuous_ANCOVA_3$S)
power_multiple_continuous_ANCOVA_3$S=factor(power_multiple_continuous_ANCOVA_3$S, levels=unique(power_multiple_continuous_ANCOVA_3$S))
power_multiple_continuous_ANCOVA_3

power_multiple_continuous_ANCOVA_n50_3=power_multiple_continuous_ANCOVA_3[power_multiple_continuous_ANCOVA_3$n=='n=50',][Sopt_ANCOVA['Sopt_ANCOVA_n50',3],]
power_multiple_continuous_ANCOVA_n100_3=power_multiple_continuous_ANCOVA_3[power_multiple_continuous_ANCOVA_3$n=='n=100',][Sopt_ANCOVA['Sopt_ANCOVA_n100',3],]
power_multiple_continuous_ANCOVA_n150_3=power_multiple_continuous_ANCOVA_3[power_multiple_continuous_ANCOVA_3$n=='n=150',][Sopt_ANCOVA['Sopt_ANCOVA_n150',3],]

Sopt_ANCOVA_multiple_continuous_3=rbind(power_multiple_continuous_ANCOVA_n50_3,
                             power_multiple_continuous_ANCOVA_n100_3,
                             power_multiple_continuous_ANCOVA_n150_3)

Plot_power_multiple_continuous_ANCOVA_3<-ggplot(power_multiple_continuous_ANCOVA_3, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.5~','~rho[X]==rho[Y] ~'=0.8'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_ANCOVA_multiple_continuous_3, 
             aes(x=S,y=power), color=brewer.pal(9, "Paired")[6],size=7)+
  theme(text = element_text(size = 45),
        axis.title=element_text(size=40),
        axis.text=element_text(size=35),
        #axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_rect(fill = "white", linetype = "solid",color = "black"),
        # Add axis line
        axis.line = element_line(color = "black"),
        legend.position = 'top',
        legend.title=element_text(size =38),
        legend.text =element_text(size =38),plot.margin=unit(c(1,2,-0.1,2), "cm"))+
  guides(color = guide_legend(title ="Sample sizes per group",nrow = 1,byrow=TRUE)) #bquote(~n[0]==n[1]~'=')
Plot_power_multiple_continuous_ANCOVA_3

#rhoXY= 0.5,rhoX= 0.9;
#rhoXY= 0.6,rhoX= 0.7; rhoXY= 0.6,rhoX= 0.8; rhoXY= 0.6,rhoX= 0.9;
#rhoXY= 0.7,rhoX= 0.8; rhoXY= 0.7,rhoX= 0.9

#rhoXY= 0.5,rhoX= 0.9
power_multiple_continuous_ANCOVA_4 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_continuous_ANCOVA_n50[4,],
                                              power_multiple_continuous_ANCOVA_n100[4,],
                                              power_multiple_continuous_ANCOVA_n150[4,]))
power_multiple_continuous_ANCOVA_4$n=as.character(power_multiple_continuous_ANCOVA_4$n)
power_multiple_continuous_ANCOVA_4$n=factor(power_multiple_continuous_ANCOVA_4$n, levels=unique(power_multiple_continuous_ANCOVA_4$n))

power_multiple_continuous_ANCOVA_4$S=as.character(power_multiple_continuous_ANCOVA_4$S)
power_multiple_continuous_ANCOVA_4$S=factor(power_multiple_continuous_ANCOVA_4$S, levels=unique(power_multiple_continuous_ANCOVA_4$S))
power_multiple_continuous_ANCOVA_4

power_multiple_continuous_ANCOVA_n50_4=power_multiple_continuous_ANCOVA_4[power_multiple_continuous_ANCOVA_4$n=='n=50',][Sopt_ANCOVA['Sopt_ANCOVA_n50',4],]
power_multiple_continuous_ANCOVA_n100_4=power_multiple_continuous_ANCOVA_4[power_multiple_continuous_ANCOVA_4$n=='n=100',][Sopt_ANCOVA['Sopt_ANCOVA_n100',4],]
power_multiple_continuous_ANCOVA_n150_4=power_multiple_continuous_ANCOVA_4[power_multiple_continuous_ANCOVA_4$n=='n=150',][Sopt_ANCOVA['Sopt_ANCOVA_n150',4],]

Sopt_ANCOVA_multiple_continuous_4=rbind(power_multiple_continuous_ANCOVA_n50_4,
                             power_multiple_continuous_ANCOVA_n100_4,
                             power_multiple_continuous_ANCOVA_n150_4)

Plot_power_multiple_continuous_ANCOVA_4<-ggplot(power_multiple_continuous_ANCOVA_4, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.5~','~rho[X]==rho[Y] ~'=0.9'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_ANCOVA_multiple_continuous_4, 
             aes(x=S,y=power), color=brewer.pal(9, "Paired")[6],size=7)+
  theme(text = element_text(size = 45),
        axis.title=element_text(size=40),
        axis.text=element_text(size=35),
        #axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_rect(fill = "white", linetype = "solid",color = "black"),
        # Add axis line
        axis.line = element_line(color = "black"),
        legend.position = 'top',
        legend.title=element_text(size =38),
        legend.text =element_text(size =38),plot.margin=unit(c(1,2,-0.1,2), "cm"))+
  guides(color = guide_legend(title ="Sample sizes per group",nrow = 1,byrow=TRUE)) #bquote(~n[0]==n[1]~'=')
Plot_power_multiple_continuous_ANCOVA_4

#rhoXY= 0.6,rhoX= 0.7
power_multiple_continuous_ANCOVA_5 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_continuous_ANCOVA_n50[5,],
                                              power_multiple_continuous_ANCOVA_n100[5,],
                                              power_multiple_continuous_ANCOVA_n150[5,]))
power_multiple_continuous_ANCOVA_5$n=as.character(power_multiple_continuous_ANCOVA_5$n)
power_multiple_continuous_ANCOVA_5$n=factor(power_multiple_continuous_ANCOVA_5$n, levels=unique(power_multiple_continuous_ANCOVA_5$n))

power_multiple_continuous_ANCOVA_5$S=as.character(power_multiple_continuous_ANCOVA_5$S)
power_multiple_continuous_ANCOVA_5$S=factor(power_multiple_continuous_ANCOVA_5$S, levels=unique(power_multiple_continuous_ANCOVA_5$S))
power_multiple_continuous_ANCOVA_5

power_multiple_continuous_ANCOVA_n50_5=power_multiple_continuous_ANCOVA_5[power_multiple_continuous_ANCOVA_5$n=='n=50',][Sopt_ANCOVA['Sopt_ANCOVA_n50',5],]
power_multiple_continuous_ANCOVA_n100_5=power_multiple_continuous_ANCOVA_5[power_multiple_continuous_ANCOVA_5$n=='n=100',][Sopt_ANCOVA['Sopt_ANCOVA_n100',5],]
power_multiple_continuous_ANCOVA_n150_5=power_multiple_continuous_ANCOVA_5[power_multiple_continuous_ANCOVA_5$n=='n=150',][Sopt_ANCOVA['Sopt_ANCOVA_n150',5],]

Sopt_ANCOVA_multiple_continuous_5=rbind(power_multiple_continuous_ANCOVA_n50_5,
                             power_multiple_continuous_ANCOVA_n100_5,
                             power_multiple_continuous_ANCOVA_n150_5)

Plot_power_multiple_continuous_ANCOVA_5<-ggplot(power_multiple_continuous_ANCOVA_5, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.6~','~rho[X]==rho[Y] ~'=0.7'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_ANCOVA_multiple_continuous_5, 
             aes(x=S,y=power), color=brewer.pal(9, "Paired")[6],size=7)+
  theme(text = element_text(size = 45),
        axis.title=element_text(size=40),
        axis.text=element_text(size=35),
        #axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_rect(fill = "white", linetype = "solid",color = "black"),
        # Add axis line
        axis.line = element_line(color = "black"),
        legend.position = 'top',
        legend.title=element_text(size =38),
        legend.text =element_text(size =38),plot.margin=unit(c(1,2,-0.1,2), "cm"))+
  guides(color = guide_legend(title ="Sample sizes per group",nrow = 1,byrow=TRUE)) #bquote(~n[0]==n[1]~'=')
Plot_power_multiple_continuous_ANCOVA_5

#rhoXY= 0.6,rhoX= 0.8
power_multiple_continuous_ANCOVA_6 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_continuous_ANCOVA_n50[6,],
                                              power_multiple_continuous_ANCOVA_n100[6,],
                                              power_multiple_continuous_ANCOVA_n150[6,]))
power_multiple_continuous_ANCOVA_6$n=as.character(power_multiple_continuous_ANCOVA_6$n)
power_multiple_continuous_ANCOVA_6$n=factor(power_multiple_continuous_ANCOVA_6$n, levels=unique(power_multiple_continuous_ANCOVA_6$n))

power_multiple_continuous_ANCOVA_6$S=as.character(power_multiple_continuous_ANCOVA_6$S)
power_multiple_continuous_ANCOVA_6$S=factor(power_multiple_continuous_ANCOVA_6$S, levels=unique(power_multiple_continuous_ANCOVA_6$S))
power_multiple_continuous_ANCOVA_6

power_multiple_continuous_ANCOVA_n50_6=power_multiple_continuous_ANCOVA_6[power_multiple_continuous_ANCOVA_6$n=='n=50',][Sopt_ANCOVA['Sopt_ANCOVA_n50',6],]
power_multiple_continuous_ANCOVA_n100_6=power_multiple_continuous_ANCOVA_6[power_multiple_continuous_ANCOVA_6$n=='n=100',][Sopt_ANCOVA['Sopt_ANCOVA_n100',6],]
power_multiple_continuous_ANCOVA_n150_6=power_multiple_continuous_ANCOVA_6[power_multiple_continuous_ANCOVA_6$n=='n=150',][Sopt_ANCOVA['Sopt_ANCOVA_n150',6],]

Sopt_ANCOVA_multiple_continuous_6=rbind(power_multiple_continuous_ANCOVA_n50_6,
                             power_multiple_continuous_ANCOVA_n100_6,
                             power_multiple_continuous_ANCOVA_n150_6)

Plot_power_multiple_continuous_ANCOVA_6<-ggplot(power_multiple_continuous_ANCOVA_6, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.6~','~rho[X]==rho[Y] ~'=0.8'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_ANCOVA_multiple_continuous_6, 
             aes(x=S,y=power), color=brewer.pal(9, "Paired")[6],size=7)+
  theme(text = element_text(size = 45),
        axis.title=element_text(size=40),
        axis.text=element_text(size=35),
        #axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_rect(fill = "white", linetype = "solid",color = "black"),
        # Add axis line
        axis.line = element_line(color = "black"),
        legend.position = 'top',
        legend.title=element_text(size =38),
        legend.text =element_text(size =38),plot.margin=unit(c(1,2,-0.1,2), "cm"))+
  guides(color = guide_legend(title ="Sample sizes per group",nrow = 1,byrow=TRUE)) #bquote(~n[0]==n[1]~'=')
Plot_power_multiple_continuous_ANCOVA_6

#rhoXY= 0.6,rhoX= 0.9
power_multiple_continuous_ANCOVA_7 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_continuous_ANCOVA_n50[7,],
                                              power_multiple_continuous_ANCOVA_n100[7,],
                                              power_multiple_continuous_ANCOVA_n150[7,]))
power_multiple_continuous_ANCOVA_7$n=as.character(power_multiple_continuous_ANCOVA_7$n)
power_multiple_continuous_ANCOVA_7$n=factor(power_multiple_continuous_ANCOVA_7$n, levels=unique(power_multiple_continuous_ANCOVA_7$n))

power_multiple_continuous_ANCOVA_7$S=as.character(power_multiple_continuous_ANCOVA_7$S)
power_multiple_continuous_ANCOVA_7$S=factor(power_multiple_continuous_ANCOVA_7$S, levels=unique(power_multiple_continuous_ANCOVA_7$S))
power_multiple_continuous_ANCOVA_7

power_multiple_continuous_ANCOVA_n50_7=power_multiple_continuous_ANCOVA_7[power_multiple_continuous_ANCOVA_7$n=='n=50',][Sopt_ANCOVA['Sopt_ANCOVA_n50',7],]
power_multiple_continuous_ANCOVA_n100_7=power_multiple_continuous_ANCOVA_7[power_multiple_continuous_ANCOVA_7$n=='n=100',][Sopt_ANCOVA['Sopt_ANCOVA_n100',7],]
power_multiple_continuous_ANCOVA_n150_7=power_multiple_continuous_ANCOVA_7[power_multiple_continuous_ANCOVA_7$n=='n=150',][Sopt_ANCOVA['Sopt_ANCOVA_n150',7],]

Sopt_ANCOVA_multiple_continuous_7=rbind(power_multiple_continuous_ANCOVA_n50_7,
                             power_multiple_continuous_ANCOVA_n100_7,
                             power_multiple_continuous_ANCOVA_n150_7)

Plot_power_multiple_continuous_ANCOVA_7<-ggplot(power_multiple_continuous_ANCOVA_7, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.6~','~rho[X]==rho[Y] ~'=0.9'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_ANCOVA_multiple_continuous_7, 
             aes(x=S,y=power), color=brewer.pal(9, "Paired")[6],size=7)+
  theme(text = element_text(size = 45),
        axis.title=element_text(size=40),
        axis.text=element_text(size=35),
        #axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_rect(fill = "white", linetype = "solid",color = "black"),
        # Add axis line
        axis.line = element_line(color = "black"),
        legend.position = 'top',
        legend.title=element_text(size =38),
        legend.text =element_text(size =38),plot.margin=unit(c(1,2,-0.1,2), "cm"))+
  guides(color = guide_legend(title ="Sample sizes per group",nrow = 1,byrow=TRUE)) #bquote(~n[0]==n[1]~'=')
Plot_power_multiple_continuous_ANCOVA_7

#rhoXY= 0.7,rhoX= 0.8
power_multiple_continuous_ANCOVA_8 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_continuous_ANCOVA_n50[8,],
                                              power_multiple_continuous_ANCOVA_n100[8,],
                                              power_multiple_continuous_ANCOVA_n150[8,]))
power_multiple_continuous_ANCOVA_8$n=as.character(power_multiple_continuous_ANCOVA_8$n)
power_multiple_continuous_ANCOVA_8$n=factor(power_multiple_continuous_ANCOVA_8$n, levels=unique(power_multiple_continuous_ANCOVA_8$n))

power_multiple_continuous_ANCOVA_8$S=as.character(power_multiple_continuous_ANCOVA_8$S)
power_multiple_continuous_ANCOVA_8$S=factor(power_multiple_continuous_ANCOVA_8$S, levels=unique(power_multiple_continuous_ANCOVA_8$S))
power_multiple_continuous_ANCOVA_8

power_multiple_continuous_ANCOVA_n50_8=power_multiple_continuous_ANCOVA_8[power_multiple_continuous_ANCOVA_8$n=='n=50',][Sopt_ANCOVA['Sopt_ANCOVA_n50',8],]
power_multiple_continuous_ANCOVA_n100_8=power_multiple_continuous_ANCOVA_8[power_multiple_continuous_ANCOVA_8$n=='n=100',][Sopt_ANCOVA['Sopt_ANCOVA_n100',8],]
power_multiple_continuous_ANCOVA_n150_8=power_multiple_continuous_ANCOVA_8[power_multiple_continuous_ANCOVA_8$n=='n=150',][Sopt_ANCOVA['Sopt_ANCOVA_n150',8],]

Sopt_ANCOVA_multiple_continuous_8=rbind(power_multiple_continuous_ANCOVA_n50_8,
                             power_multiple_continuous_ANCOVA_n100_8,
                             power_multiple_continuous_ANCOVA_n150_8)

Plot_power_multiple_continuous_ANCOVA_8<-ggplot(power_multiple_continuous_ANCOVA_8, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.7~','~rho[X]==rho[Y] ~'=0.8'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_ANCOVA_multiple_continuous_8, 
             aes(x=S,y=power), color=brewer.pal(9, "Paired")[6],size=7)+
  theme(text = element_text(size = 45),
        axis.title=element_text(size=40),
        axis.text=element_text(size=35),
        #axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_rect(fill = "white", linetype = "solid",color = "black"),
        # Add axis line
        axis.line = element_line(color = "black"),
        legend.position = 'top',
        legend.title=element_text(size =38),
        legend.text =element_text(size =38),plot.margin=unit(c(1,2,-0.1,2), "cm"))+
  guides(color = guide_legend(title ="Sample sizes per group",nrow = 1,byrow=TRUE)) #bquote(~n[0]==n[1]~'=')
Plot_power_multiple_continuous_ANCOVA_8

#rhoXY= 0.7,rhoX= 0.9
power_multiple_continuous_ANCOVA_9 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_continuous_ANCOVA_n50[9,],
                                              power_multiple_continuous_ANCOVA_n100[9,],
                                              power_multiple_continuous_ANCOVA_n150[9,]))
power_multiple_continuous_ANCOVA_9$n=as.character(power_multiple_continuous_ANCOVA_9$n)
power_multiple_continuous_ANCOVA_9$n=factor(power_multiple_continuous_ANCOVA_9$n, levels=unique(power_multiple_continuous_ANCOVA_9$n))

power_multiple_continuous_ANCOVA_9$S=as.character(power_multiple_continuous_ANCOVA_9$S)
power_multiple_continuous_ANCOVA_9$S=factor(power_multiple_continuous_ANCOVA_9$S, levels=unique(power_multiple_continuous_ANCOVA_9$S))
power_multiple_continuous_ANCOVA_9

power_multiple_continuous_ANCOVA_n50_9=power_multiple_continuous_ANCOVA_9[power_multiple_continuous_ANCOVA_9$n=='n=50',][Sopt_ANCOVA['Sopt_ANCOVA_n50',9],]
power_multiple_continuous_ANCOVA_n100_9=power_multiple_continuous_ANCOVA_9[power_multiple_continuous_ANCOVA_9$n=='n=100',][Sopt_ANCOVA['Sopt_ANCOVA_n100',9],]
power_multiple_continuous_ANCOVA_n150_9=power_multiple_continuous_ANCOVA_9[power_multiple_continuous_ANCOVA_9$n=='n=150',][Sopt_ANCOVA['Sopt_ANCOVA_n150',9],]

Sopt_ANCOVA_multiple_continuous_9=rbind(power_multiple_continuous_ANCOVA_n50_9,
                             power_multiple_continuous_ANCOVA_n100_9,
                             power_multiple_continuous_ANCOVA_n150_9)

Plot_power_multiple_continuous_ANCOVA_9<-ggplot(power_multiple_continuous_ANCOVA_9, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.7~','~rho[X]==rho[Y] ~'=0.9'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_ANCOVA_multiple_continuous_9, 
             aes(x=S,y=power), color=brewer.pal(9, "Paired")[6],size=7)+
  theme(text = element_text(size = 45),
        axis.title=element_text(size=40),
        axis.text=element_text(size=35),
        #axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_rect(fill = "white", linetype = "solid",color = "black"),
        # Add axis line
        axis.line = element_line(color = "black"),
        legend.position = 'top',
        legend.title=element_text(size =38),
        legend.text =element_text(size =38),plot.margin=unit(c(1,2,-0.1,2), "cm"))+
  guides(color = guide_legend(title ="Sample sizes per group",nrow = 1,byrow=TRUE)) #bquote(~n[0]==n[1]~'=')
Plot_power_multiple_continuous_ANCOVA_9

library(cowplot)
png(file=paste0('Plot_power_multiple_continuous_ANCOVA.png'),width=2750,height=2700,units='px',pointsize =30)
ggdraw() +
  draw_label("ANCOVA",fontface = 'bold',x=0.4,y=1,hjust = 0, vjust = 0, size=80) +
  draw_plot(Plot_power_multiple_continuous_ANCOVA_1, x = 0, y =0.66, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_continuous_ANCOVA_2, x = 0.33, y =0.66, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_continuous_ANCOVA_3, x = 0.66, y =0.66, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_continuous_ANCOVA_4, x = 0, y =0.33, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_continuous_ANCOVA_5, x = 0.33, y =0.33, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_continuous_ANCOVA_6, x = 0.66, y =0.33, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_continuous_ANCOVA_7, x = 0, y =0, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_continuous_ANCOVA_8, x = 0.33, y =0, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_continuous_ANCOVA_9, x = 0.66, y =0, width =0.33, height=0.33)+
  theme(plot.margin = margin(85, 0, 0, 0))
dev.off()