######################################################################################
### Results--single binary outcome
######################################################################################

##Type I error probabilities:
TypeI_u_0.4_n50=cbind(round(readRDS("TypeI_u_0.4_n50.Rds"),digits=4),rep(NA,9))
TypeI_u_0.4_n75=round(readRDS("TypeI_u_0.4_n75.Rds"),digits=4)
TypeI_u_0.4_n100=round(readRDS("TypeI_u_0.4_n100.Rds"),digits=4)
TypeI_u_0.4_n125=round(readRDS("TypeI_u_0.4_n125.Rds"),digits=4)
TypeI_u_0.4_n150=round(readRDS("TypeI_u_0.4_n150.Rds"),digits=4)

TypeI_single_binary_outcome=rbind(TypeI_u_0.4_n50,TypeI_u_0.4_n75,TypeI_u_0.4_n100,TypeI_u_0.4_n125,TypeI_u_0.4_n150)
colnames(TypeI_single_binary_outcome)=c('Model1','Model2','Model3')
#write.csv(TypeI_single_binary_outcome,'TypeI_single_binary_outcome.csv',row.names = TRUE)

##Power:
Power_u_0.4_n50=cbind(round(readRDS("Power_u_0.4_n50.Rds"),digits=4),rep(NA,9))
Power_u_0.4_n75=round(readRDS("Power_u_0.4_n75.Rds"),digits=4)
Power_u_0.4_n100=round(readRDS("Power_u_0.4_n100.Rds"),digits=4)
Power_u_0.4_n125=round(readRDS("Power_u_0.4_n125.Rds"),digits=4)
Power_u_0.4_n150=round(readRDS("Power_u_0.4_n150.Rds"),digits=4)

Power_single_binary_outcome=rbind(Power_u_0.4_n50,Power_u_0.4_n75,Power_u_0.4_n100,Power_u_0.4_n125,Power_u_0.4_n150)
colnames(Power_single_binary_outcome)=c('Model1','Model2','Model3')
#write.csv(Power_single_binary_outcome,'Power_single_binary_outcome.csv',row.names = TRUE)

######################################################################################
### Results--multiple binary outcomes
######################################################################################

#n0=n1=50
GEE_sp_0.4_n50_1=readRDS(file="GEE_sp_0.4_n50_1.Rds")
GEE_sp_0.4_n50_2=readRDS(file="GEE_sp_0.4_n50_2.Rds")
GEE_sp_0.4_n50_3=readRDS(file="GEE_sp_0.4_n50_3.Rds")

TypeI_multiple_binary_n50=rbind(GEE_sp_0.4_n50_1[c(1,5,9),],GEE_sp_0.4_n50_2[c(1,5,9),],GEE_sp_0.4_n50_3[c(1,5,9),])
TypeI_multiple_binary_n50=round(TypeI_multiple_binary_n50,digits=4)
colnames(TypeI_multiple_binary_n50)=1:9
#write.csv(TypeI_multiple_binary_n50,'TypeI_multiple_binary_n50.csv',row.names = TRUE)
#write.csv(GEE_sp_0.4_n50_TypeI>0.053,'TypeI_multiple_binary_n50_ind.csv',row.names = TRUE)

#n0=n1=100
GEE_sp_0.4_n100_1=readRDS(file="GEE_sp_0.4_n100_1.Rds")
GEE_sp_0.4_n100_2=readRDS(file="GEE_sp_0.4_n100_2.Rds")
GEE_sp_0.4_n100_3=readRDS(file="GEE_sp_0.4_n100_3.Rds")

TypeI_multiple_binary_n100=rbind(GEE_sp_0.4_n100_1[c(1,2,9,10,17,18),],
                                 GEE_sp_0.4_n100_2[c(1,2,9,10,17,18),],
                                 GEE_sp_0.4_n100_3[c(1,2,9,10,17,18),])
TypeI_multiple_binary_n100=round(TypeI_multiple_binary_n100,digits=4)
colnames(TypeI_multiple_binary_n100)=1:9
#write.csv(TypeI_multiple_binary_n100,'TypeI_multiple_binary_n100.csv',row.names = TRUE)
#write.csv(TypeI_multiple_binary_n100>0.053,'TypeI_multiple_binary_n100_ind.csv',row.names = TRUE)

#n0=n1=150
GEE_sp_0.4_n150_1=readRDS(file="GEE_sp_0.4_n150_1.Rds")
GEE_sp_0.4_n150_2=readRDS(file="GEE_sp_0.4_n150_2.Rds")
GEE_sp_0.4_n150_3=readRDS(file="GEE_sp_0.4_n150_3.Rds")

TypeI_multiple_binary_n150=rbind(GEE_sp_0.4_n150_1[c(1,2,9,10,17,18),],
                                 GEE_sp_0.4_n150_2[c(1,2,9,10,17,18),],
                                 GEE_sp_0.4_n150_3[c(1,2,9,10,17,18),])
TypeI_multiple_binary_n150=round(TypeI_multiple_binary_n150,digits=4)
colnames(TypeI_multiple_binary_n150)=1:9
#write.csv(TypeI_multiple_binary_n150,'TypeI_multiple_binary_n150.csv',row.names = TRUE)
#write.csv(TypeI_multiple_binary_n150>0.053,'TypeI_multiple_binary_n150_ind.csv',row.names = TRUE)

sum(TypeI_multiple_binary_n50>0.053) #42
sum(TypeI_multiple_binary_n100>0.053) #24/2
sum(TypeI_multiple_binary_n150>0.053) #27/2

### Power plots--multiple binary outcomes
library(ggplot2)
library(RColorBrewer)
library(ggpattern)
library(ggpubr)
library(splines)

#power plot: GEE model 2; n=50, 100, 150
power_multiple_binary_GEE2_n50=rbind(GEE_sp_0.4_n50_1[c(4,8,12),],GEE_sp_0.4_n50_2[c(4,8,12),],GEE_sp_0.4_n50_3[c(4,8,12),])
colnames(power_multiple_binary_GEE2_n50)=1:9

power_multiple_binary_GEE2_n100=rbind(GEE_sp_0.4_n100_1[c(6,14,22),],GEE_sp_0.4_n100_2[c(6,14,22),],GEE_sp_0.4_n100_3[c(6,14,22),])
colnames(power_multiple_binary_GEE2_n100)=1:9

power_multiple_binary_GEE2_n150=rbind(GEE_sp_0.4_n150_1[c(6,14,22),],GEE_sp_0.4_n150_2[c(6,14,22),],GEE_sp_0.4_n150_3[c(6,14,22),])
colnames(power_multiple_binary_GEE2_n150)=1:9

#maximum power of each row
apply(power_multiple_binary_GEE2_n50, 1, max)
#Sopt_GEE2
Sopt_GEE2_n50=rep(NA,9)
for (i in 1:9){
  Sopt_GEE2_n50[i]=which(power_multiple_binary_GEE2_n50[i,]%in%max(power_multiple_binary_GEE2_n50[i,]))
}
Sopt_GEE2_n50 #4 4 3 3 4 4 3 4 4

apply(power_multiple_binary_GEE2_n100, 1, max)
#Sopt_GEE2
Sopt_GEE2_n100=rep(NA,9)
for (i in 1:9){
  Sopt_GEE2_n100[i]=which(power_multiple_binary_GEE2_n100[i,]%in%max(power_multiple_binary_GEE2_n100[i,]))
}
Sopt_GEE2_n100 #4 3 3 3 4 4 4 4 4

apply(power_multiple_binary_GEE2_n150, 1, max)
#Sopt_GEE2
Sopt_GEE2_n150=rep(NA,9)
for (i in 1:9){
  Sopt_GEE2_n150[i]=which(power_multiple_binary_GEE2_n150[i,]%in%max(power_multiple_binary_GEE2_n150[i,]))
}
Sopt_GEE2_n150 #3 3 3 3 5 3 3 5 5

Sopt_GEE2=rbind(Sopt_GEE2_n50,Sopt_GEE2_n100,Sopt_GEE2_n150)

#highlight Sopt_GEE2
#display.brewer.all()
cols1 <- c(brewer.pal(9, "Pastel1")[1],brewer.pal(9, "Paired")[1],brewer.pal(9, "Paired")[2])

#9 plots: rhoXY= 0.5,rhoX= 0.6; rhoXY= 0.5,rhoX= 0.7; rhoXY= 0.5,rhoX= 0.8; rhoXY= 0.5,rhoX= 0.9;
#rhoXY= 0.6,rhoX= 0.7; rhoXY= 0.6,rhoX= 0.8; rhoXY= 0.6,rhoX= 0.9;
#rhoXY= 0.7,rhoX= 0.8; rhoXY= 0.7,rhoX= 0.9

#rhoXY= 0.5,rhoX= 0.6
power_multiple_binary_GEE2_1 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_binary_GEE2_n50[1,],
                                              power_multiple_binary_GEE2_n100[1,],
                                              power_multiple_binary_GEE2_n150[1,]))
power_multiple_binary_GEE2_1$n=as.character(power_multiple_binary_GEE2_1$n)
power_multiple_binary_GEE2_1$n=factor(power_multiple_binary_GEE2_1$n, levels=unique(power_multiple_binary_GEE2_1$n))

power_multiple_binary_GEE2_1$S=as.character(power_multiple_binary_GEE2_1$S)
power_multiple_binary_GEE2_1$S=factor(power_multiple_binary_GEE2_1$S, levels=unique(power_multiple_binary_GEE2_1$S))
power_multiple_binary_GEE2_1

power_multiple_binary_GEE2_n50_1=power_multiple_binary_GEE2_1[power_multiple_binary_GEE2_1$n=='n=50',][Sopt_GEE2['Sopt_GEE2_n50',1],]
power_multiple_binary_GEE2_n100_1=power_multiple_binary_GEE2_1[power_multiple_binary_GEE2_1$n=='n=100',][Sopt_GEE2['Sopt_GEE2_n100',1],]
power_multiple_binary_GEE2_n150_1=power_multiple_binary_GEE2_1[power_multiple_binary_GEE2_1$n=='n=150',][Sopt_GEE2['Sopt_GEE2_n150',1],]

Sopt_GEE2_multiple_binary_1=rbind(power_multiple_binary_GEE2_n50_1,
                             power_multiple_binary_GEE2_n100_1,
                             power_multiple_binary_GEE2_n150_1)

Plot_power_multiple_binary_GEE2_1<-ggplot(power_multiple_binary_GEE2_1, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.5~','~rho[X]==rho[Y] ~'=0.6'))+  
  scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE2_multiple_binary_1, 
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
Plot_power_multiple_binary_GEE2_1

#rhoXY= 0.5,rhoX= 0.7
power_multiple_binary_GEE2_2 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_binary_GEE2_n50[2,],
                                              power_multiple_binary_GEE2_n100[2,],
                                              power_multiple_binary_GEE2_n150[2,]))
power_multiple_binary_GEE2_2$n=as.character(power_multiple_binary_GEE2_2$n)
power_multiple_binary_GEE2_2$n=factor(power_multiple_binary_GEE2_2$n, levels=unique(power_multiple_binary_GEE2_2$n))

power_multiple_binary_GEE2_2$S=as.character(power_multiple_binary_GEE2_2$S)
power_multiple_binary_GEE2_2$S=factor(power_multiple_binary_GEE2_2$S, levels=unique(power_multiple_binary_GEE2_2$S))
power_multiple_binary_GEE2_2

power_multiple_binary_GEE2_n50_2=power_multiple_binary_GEE2_2[power_multiple_binary_GEE2_2$n=='n=50',][Sopt_GEE2['Sopt_GEE2_n50',2],]
power_multiple_binary_GEE2_n100_2=power_multiple_binary_GEE2_2[power_multiple_binary_GEE2_2$n=='n=100',][Sopt_GEE2['Sopt_GEE2_n100',2],]
power_multiple_binary_GEE2_n150_2=power_multiple_binary_GEE2_2[power_multiple_binary_GEE2_2$n=='n=150',][Sopt_GEE2['Sopt_GEE2_n150',2],]

Sopt_GEE2_multiple_binary_2=rbind(power_multiple_binary_GEE2_n50_2,
                             power_multiple_binary_GEE2_n100_2,
                             power_multiple_binary_GEE2_n150_2)

Plot_power_multiple_binary_GEE2_2<-ggplot(power_multiple_binary_GEE2_2, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.5~','~rho[X]==rho[Y] ~'=0.7'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE2_multiple_binary_2, 
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
Plot_power_multiple_binary_GEE2_2

#rhoXY= 0.5,rhoX= 0.8
power_multiple_binary_GEE2_3 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_binary_GEE2_n50[3,],
                                              power_multiple_binary_GEE2_n100[3,],
                                              power_multiple_binary_GEE2_n150[3,]))
power_multiple_binary_GEE2_3$n=as.character(power_multiple_binary_GEE2_3$n)
power_multiple_binary_GEE2_3$n=factor(power_multiple_binary_GEE2_3$n, levels=unique(power_multiple_binary_GEE2_3$n))

power_multiple_binary_GEE2_3$S=as.character(power_multiple_binary_GEE2_3$S)
power_multiple_binary_GEE2_3$S=factor(power_multiple_binary_GEE2_3$S, levels=unique(power_multiple_binary_GEE2_3$S))
power_multiple_binary_GEE2_3

power_multiple_binary_GEE2_n50_3=power_multiple_binary_GEE2_3[power_multiple_binary_GEE2_3$n=='n=50',][Sopt_GEE2['Sopt_GEE2_n50',3],]
power_multiple_binary_GEE2_n100_3=power_multiple_binary_GEE2_3[power_multiple_binary_GEE2_3$n=='n=100',][Sopt_GEE2['Sopt_GEE2_n100',3],]
power_multiple_binary_GEE2_n150_3=power_multiple_binary_GEE2_3[power_multiple_binary_GEE2_3$n=='n=150',][Sopt_GEE2['Sopt_GEE2_n150',3],]

Sopt_GEE2_multiple_binary_3=rbind(power_multiple_binary_GEE2_n50_3,
                             power_multiple_binary_GEE2_n100_3,
                             power_multiple_binary_GEE2_n150_3)

Plot_power_multiple_binary_GEE2_3<-ggplot(power_multiple_binary_GEE2_3, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.5~','~rho[X]==rho[Y] ~'=0.8'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE2_multiple_binary_3, 
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
Plot_power_multiple_binary_GEE2_3

#rhoXY= 0.5,rhoX= 0.9;
#rhoXY= 0.6,rhoX= 0.7; rhoXY= 0.6,rhoX= 0.8; rhoXY= 0.6,rhoX= 0.9;
#rhoXY= 0.7,rhoX= 0.8; rhoXY= 0.7,rhoX= 0.9

#rhoXY= 0.5,rhoX= 0.9
power_multiple_binary_GEE2_4 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_binary_GEE2_n50[4,],
                                              power_multiple_binary_GEE2_n100[4,],
                                              power_multiple_binary_GEE2_n150[4,]))
power_multiple_binary_GEE2_4$n=as.character(power_multiple_binary_GEE2_4$n)
power_multiple_binary_GEE2_4$n=factor(power_multiple_binary_GEE2_4$n, levels=unique(power_multiple_binary_GEE2_4$n))

power_multiple_binary_GEE2_4$S=as.character(power_multiple_binary_GEE2_4$S)
power_multiple_binary_GEE2_4$S=factor(power_multiple_binary_GEE2_4$S, levels=unique(power_multiple_binary_GEE2_4$S))
power_multiple_binary_GEE2_4

power_multiple_binary_GEE2_n50_4=power_multiple_binary_GEE2_4[power_multiple_binary_GEE2_4$n=='n=50',][Sopt_GEE2['Sopt_GEE2_n50',4],]
power_multiple_binary_GEE2_n100_4=power_multiple_binary_GEE2_4[power_multiple_binary_GEE2_4$n=='n=100',][Sopt_GEE2['Sopt_GEE2_n100',4],]
power_multiple_binary_GEE2_n150_4=power_multiple_binary_GEE2_4[power_multiple_binary_GEE2_4$n=='n=150',][Sopt_GEE2['Sopt_GEE2_n150',4],]

Sopt_GEE2_multiple_binary_4=rbind(power_multiple_binary_GEE2_n50_4,
                             power_multiple_binary_GEE2_n100_4,
                             power_multiple_binary_GEE2_n150_4)

Plot_power_multiple_binary_GEE2_4<-ggplot(power_multiple_binary_GEE2_4, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.5~','~rho[X]==rho[Y] ~'=0.9'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE2_multiple_binary_4, 
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
Plot_power_multiple_binary_GEE2_4

#rhoXY= 0.6,rhoX= 0.7
power_multiple_binary_GEE2_5 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_binary_GEE2_n50[5,],
                                              power_multiple_binary_GEE2_n100[5,],
                                              power_multiple_binary_GEE2_n150[5,]))
power_multiple_binary_GEE2_5$n=as.character(power_multiple_binary_GEE2_5$n)
power_multiple_binary_GEE2_5$n=factor(power_multiple_binary_GEE2_5$n, levels=unique(power_multiple_binary_GEE2_5$n))

power_multiple_binary_GEE2_5$S=as.character(power_multiple_binary_GEE2_5$S)
power_multiple_binary_GEE2_5$S=factor(power_multiple_binary_GEE2_5$S, levels=unique(power_multiple_binary_GEE2_5$S))
power_multiple_binary_GEE2_5

power_multiple_binary_GEE2_n50_5=power_multiple_binary_GEE2_5[power_multiple_binary_GEE2_5$n=='n=50',][Sopt_GEE2['Sopt_GEE2_n50',5],]
power_multiple_binary_GEE2_n100_5=power_multiple_binary_GEE2_5[power_multiple_binary_GEE2_5$n=='n=100',][Sopt_GEE2['Sopt_GEE2_n100',5],]
power_multiple_binary_GEE2_n150_5=power_multiple_binary_GEE2_5[power_multiple_binary_GEE2_5$n=='n=150',][Sopt_GEE2['Sopt_GEE2_n150',5],]

Sopt_GEE2_multiple_binary_5=rbind(power_multiple_binary_GEE2_n50_5,
                             power_multiple_binary_GEE2_n100_5,
                             power_multiple_binary_GEE2_n150_5)

Plot_power_multiple_binary_GEE2_5<-ggplot(power_multiple_binary_GEE2_5, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.6~','~rho[X]==rho[Y] ~'=0.7'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE2_multiple_binary_5, 
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
Plot_power_multiple_binary_GEE2_5

#rhoXY= 0.6,rhoX= 0.8
power_multiple_binary_GEE2_6 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_binary_GEE2_n50[6,],
                                              power_multiple_binary_GEE2_n100[6,],
                                              power_multiple_binary_GEE2_n150[6,]))
power_multiple_binary_GEE2_6$n=as.character(power_multiple_binary_GEE2_6$n)
power_multiple_binary_GEE2_6$n=factor(power_multiple_binary_GEE2_6$n, levels=unique(power_multiple_binary_GEE2_6$n))

power_multiple_binary_GEE2_6$S=as.character(power_multiple_binary_GEE2_6$S)
power_multiple_binary_GEE2_6$S=factor(power_multiple_binary_GEE2_6$S, levels=unique(power_multiple_binary_GEE2_6$S))
power_multiple_binary_GEE2_6

power_multiple_binary_GEE2_n50_6=power_multiple_binary_GEE2_6[power_multiple_binary_GEE2_6$n=='n=50',][Sopt_GEE2['Sopt_GEE2_n50',6],]
power_multiple_binary_GEE2_n100_6=power_multiple_binary_GEE2_6[power_multiple_binary_GEE2_6$n=='n=100',][Sopt_GEE2['Sopt_GEE2_n100',6],]
power_multiple_binary_GEE2_n150_6=power_multiple_binary_GEE2_6[power_multiple_binary_GEE2_6$n=='n=150',][Sopt_GEE2['Sopt_GEE2_n150',6],]

Sopt_GEE2_multiple_binary_6=rbind(power_multiple_binary_GEE2_n50_6,
                             power_multiple_binary_GEE2_n100_6,
                             power_multiple_binary_GEE2_n150_6)

Plot_power_multiple_binary_GEE2_6<-ggplot(power_multiple_binary_GEE2_6, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.6~','~rho[X]==rho[Y] ~'=0.8'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE2_multiple_binary_6, 
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
Plot_power_multiple_binary_GEE2_6

#rhoXY= 0.6,rhoX= 0.9
power_multiple_binary_GEE2_7 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_binary_GEE2_n50[7,],
                                              power_multiple_binary_GEE2_n100[7,],
                                              power_multiple_binary_GEE2_n150[7,]))
power_multiple_binary_GEE2_7$n=as.character(power_multiple_binary_GEE2_7$n)
power_multiple_binary_GEE2_7$n=factor(power_multiple_binary_GEE2_7$n, levels=unique(power_multiple_binary_GEE2_7$n))

power_multiple_binary_GEE2_7$S=as.character(power_multiple_binary_GEE2_7$S)
power_multiple_binary_GEE2_7$S=factor(power_multiple_binary_GEE2_7$S, levels=unique(power_multiple_binary_GEE2_7$S))
power_multiple_binary_GEE2_7

power_multiple_binary_GEE2_n50_7=power_multiple_binary_GEE2_7[power_multiple_binary_GEE2_7$n=='n=50',][Sopt_GEE2['Sopt_GEE2_n50',7],]
power_multiple_binary_GEE2_n100_7=power_multiple_binary_GEE2_7[power_multiple_binary_GEE2_7$n=='n=100',][Sopt_GEE2['Sopt_GEE2_n100',7],]
power_multiple_binary_GEE2_n150_7=power_multiple_binary_GEE2_7[power_multiple_binary_GEE2_7$n=='n=150',][Sopt_GEE2['Sopt_GEE2_n150',7],]

Sopt_GEE2_multiple_binary_7=rbind(power_multiple_binary_GEE2_n50_7,
                             power_multiple_binary_GEE2_n100_7,
                             power_multiple_binary_GEE2_n150_7)

Plot_power_multiple_binary_GEE2_7<-ggplot(power_multiple_binary_GEE2_7, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.6~','~rho[X]==rho[Y] ~'=0.9'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE2_multiple_binary_7, 
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
Plot_power_multiple_binary_GEE2_7

#rhoXY= 0.7,rhoX= 0.8
power_multiple_binary_GEE2_8 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_binary_GEE2_n50[8,],
                                              power_multiple_binary_GEE2_n100[8,],
                                              power_multiple_binary_GEE2_n150[8,]))
power_multiple_binary_GEE2_8$n=as.character(power_multiple_binary_GEE2_8$n)
power_multiple_binary_GEE2_8$n=factor(power_multiple_binary_GEE2_8$n, levels=unique(power_multiple_binary_GEE2_8$n))

power_multiple_binary_GEE2_8$S=as.character(power_multiple_binary_GEE2_8$S)
power_multiple_binary_GEE2_8$S=factor(power_multiple_binary_GEE2_8$S, levels=unique(power_multiple_binary_GEE2_8$S))
power_multiple_binary_GEE2_8

power_multiple_binary_GEE2_n50_8=power_multiple_binary_GEE2_8[power_multiple_binary_GEE2_8$n=='n=50',][Sopt_GEE2['Sopt_GEE2_n50',8],]
power_multiple_binary_GEE2_n100_8=power_multiple_binary_GEE2_8[power_multiple_binary_GEE2_8$n=='n=100',][Sopt_GEE2['Sopt_GEE2_n100',8],]
power_multiple_binary_GEE2_n150_8=power_multiple_binary_GEE2_8[power_multiple_binary_GEE2_8$n=='n=150',][Sopt_GEE2['Sopt_GEE2_n150',8],]

Sopt_GEE2_multiple_binary_8=rbind(power_multiple_binary_GEE2_n50_8,
                             power_multiple_binary_GEE2_n100_8,
                             power_multiple_binary_GEE2_n150_8)

Plot_power_multiple_binary_GEE2_8<-ggplot(power_multiple_binary_GEE2_8, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.7~','~rho[X]==rho[Y] ~'=0.8'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE2_multiple_binary_8, 
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
Plot_power_multiple_binary_GEE2_8

#rhoXY= 0.7,rhoX= 0.9
power_multiple_binary_GEE2_9 <- data.frame(n=rep(c('n=50','n=100','n=150'), each=9),
                                      S=rep(1:9,3),
                                      power=c(power_multiple_binary_GEE2_n50[9,],
                                              power_multiple_binary_GEE2_n100[9,],
                                              power_multiple_binary_GEE2_n150[9,]))
power_multiple_binary_GEE2_9$n=as.character(power_multiple_binary_GEE2_9$n)
power_multiple_binary_GEE2_9$n=factor(power_multiple_binary_GEE2_9$n, levels=unique(power_multiple_binary_GEE2_9$n))

power_multiple_binary_GEE2_9$S=as.character(power_multiple_binary_GEE2_9$S)
power_multiple_binary_GEE2_9$S=factor(power_multiple_binary_GEE2_9$S, levels=unique(power_multiple_binary_GEE2_9$S))
power_multiple_binary_GEE2_9

power_multiple_binary_GEE2_n50_9=power_multiple_binary_GEE2_9[power_multiple_binary_GEE2_9$n=='n=50',][Sopt_GEE2['Sopt_GEE2_n50',9],]
power_multiple_binary_GEE2_n100_9=power_multiple_binary_GEE2_9[power_multiple_binary_GEE2_9$n=='n=100',][Sopt_GEE2['Sopt_GEE2_n100',9],]
power_multiple_binary_GEE2_n150_9=power_multiple_binary_GEE2_9[power_multiple_binary_GEE2_9$n=='n=150',][Sopt_GEE2['Sopt_GEE2_n150',9],]

Sopt_GEE2_multiple_binary_9=rbind(power_multiple_binary_GEE2_n50_9,
                             power_multiple_binary_GEE2_n100_9,
                             power_multiple_binary_GEE2_n150_9)

Plot_power_multiple_binary_GEE2_9<-ggplot(power_multiple_binary_GEE2_9, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols1)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.7~','~rho[X]==rho[Y] ~'=0.9'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE2_multiple_binary_9, 
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
Plot_power_multiple_binary_GEE2_9

library(cowplot)
png(file=paste0('Plot_power_multiple_binary_GEE2.png'),width=2750,height=2700,units='px',pointsize =30)
ggdraw() +
  draw_label("GEE Model 2",fontface = 'bold',x=0.4,y=1,hjust = 0, vjust = 0, size=80) +
  draw_plot(Plot_power_multiple_binary_GEE2_1, x = 0, y =0.66, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE2_2, x = 0.33, y =0.66, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE2_3, x = 0.66, y =0.66, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE2_4, x = 0, y =0.33, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE2_5, x = 0.33, y =0.33, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE2_6, x = 0.66, y =0.33, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE2_7, x = 0, y =0, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE2_8, x = 0.33, y =0, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE2_9, x = 0.66, y =0, width =0.33, height=0.33)+
  theme(plot.margin = margin(85, 0, 0, 0))
dev.off()


#power plot: GEE model 1; n=100, 150
power_multiple_binary_GEE1_n100=rbind(GEE_sp_0.4_n100_1[c(5,13,21),],GEE_sp_0.4_n100_2[c(5,13,21),],GEE_sp_0.4_n100_3[c(5,13,21),])
colnames(power_multiple_binary_GEE1_n100)=1:9

power_multiple_binary_GEE1_n150=rbind(GEE_sp_0.4_n150_1[c(5,13,21),],GEE_sp_0.4_n150_2[c(5,13,21),],GEE_sp_0.4_n150_3[c(5,13,21),])
colnames(power_multiple_binary_GEE1_n150)=1:9

#maximum power of each row
apply(power_multiple_binary_GEE1_n100, 1, max)
#Sopt_GEE1
Sopt_GEE1_n100=rep(NA,9)
for (i in 1:9){
  Sopt_GEE1_n100[i]=which(power_multiple_binary_GEE1_n100[i,]%in%max(power_multiple_binary_GEE1_n100[i,]))
}
Sopt_GEE1_n100 #4 3 3 3 4 4 4 4 4

apply(power_multiple_binary_GEE1_n150, 1, max)
#Sopt_GEE1
Sopt_GEE1_n150=rep(NA,9)
for (i in 1:9){
  Sopt_GEE1_n150[i]=which(power_multiple_binary_GEE1_n150[i,]%in%max(power_multiple_binary_GEE1_n150[i,]))
}
Sopt_GEE1_n150 #3 3 3 3 5 3 3 5 5
Sopt_GEE1=rbind(Sopt_GEE1_n100,Sopt_GEE1_n150)

#highlight Sopt_GEE1
#display.brewer.all()
cols2 <- c(brewer.pal(9, "Paired")[1],brewer.pal(9, "Paired")[2])

#rhoXY= 0.5,rhoX= 0.6
power_multiple_binary_GEE1_1 <- data.frame(n=rep(c('n=100','n=150'), each=9),
                                           S=rep(1:9,2),
                                           power=c(power_multiple_binary_GEE1_n100[1,],
                                                   power_multiple_binary_GEE1_n150[1,]))
power_multiple_binary_GEE1_1$n=as.character(power_multiple_binary_GEE1_1$n)
power_multiple_binary_GEE1_1$n=factor(power_multiple_binary_GEE1_1$n, levels=unique(power_multiple_binary_GEE1_1$n))

power_multiple_binary_GEE1_1$S=as.character(power_multiple_binary_GEE1_1$S)
power_multiple_binary_GEE1_1$S=factor(power_multiple_binary_GEE1_1$S, levels=unique(power_multiple_binary_GEE1_1$S))
power_multiple_binary_GEE1_1

power_multiple_binary_GEE1_n100_1=power_multiple_binary_GEE1_1[power_multiple_binary_GEE1_1$n=='n=100',][Sopt_GEE1['Sopt_GEE1_n100',1],]
power_multiple_binary_GEE1_n150_1=power_multiple_binary_GEE1_1[power_multiple_binary_GEE1_1$n=='n=150',][Sopt_GEE1['Sopt_GEE1_n150',1],]

Sopt_GEE1_multiple_binary_1=rbind(power_multiple_binary_GEE1_n100_1,
                                  power_multiple_binary_GEE1_n150_1)

Plot_power_multiple_binary_GEE1_1<-ggplot(power_multiple_binary_GEE1_1, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols2)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.5~','~rho[X]==rho[Y] ~'=0.6'))+  
  scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE1_multiple_binary_1, 
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
Plot_power_multiple_binary_GEE1_1

#rhoXY= 0.5,rhoX= 0.7
power_multiple_binary_GEE1_2 <- data.frame(n=rep(c('n=100','n=150'), each=9),
                                           S=rep(1:9,2),
                                           power=c(power_multiple_binary_GEE1_n100[2,],
                                                   power_multiple_binary_GEE1_n150[2,]))
power_multiple_binary_GEE1_2$n=as.character(power_multiple_binary_GEE1_2$n)
power_multiple_binary_GEE1_2$n=factor(power_multiple_binary_GEE1_2$n, levels=unique(power_multiple_binary_GEE1_2$n))

power_multiple_binary_GEE1_2$S=as.character(power_multiple_binary_GEE1_2$S)
power_multiple_binary_GEE1_2$S=factor(power_multiple_binary_GEE1_2$S, levels=unique(power_multiple_binary_GEE1_2$S))
power_multiple_binary_GEE1_2

power_multiple_binary_GEE1_n100_2=power_multiple_binary_GEE1_2[power_multiple_binary_GEE1_2$n=='n=100',][Sopt_GEE1['Sopt_GEE1_n100',2],]
power_multiple_binary_GEE1_n150_2=power_multiple_binary_GEE1_2[power_multiple_binary_GEE1_2$n=='n=150',][Sopt_GEE1['Sopt_GEE1_n150',2],]

Sopt_GEE1_multiple_binary_2=rbind(power_multiple_binary_GEE1_n100_2,
                                  power_multiple_binary_GEE1_n150_2)

Plot_power_multiple_binary_GEE1_2<-ggplot(power_multiple_binary_GEE1_2, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols2)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.5~','~rho[X]==rho[Y] ~'=0.7'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE1_multiple_binary_2, 
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
Plot_power_multiple_binary_GEE1_2

#rhoXY= 0.5,rhoX= 0.8
power_multiple_binary_GEE1_3 <- data.frame(n=rep(c('n=100','n=150'), each=9),
                                           S=rep(1:9,2),
                                           power=c(power_multiple_binary_GEE1_n100[3,],
                                                   power_multiple_binary_GEE1_n150[3,]))
power_multiple_binary_GEE1_3$n=as.character(power_multiple_binary_GEE1_3$n)
power_multiple_binary_GEE1_3$n=factor(power_multiple_binary_GEE1_3$n, levels=unique(power_multiple_binary_GEE1_3$n))

power_multiple_binary_GEE1_3$S=as.character(power_multiple_binary_GEE1_3$S)
power_multiple_binary_GEE1_3$S=factor(power_multiple_binary_GEE1_3$S, levels=unique(power_multiple_binary_GEE1_3$S))
power_multiple_binary_GEE1_3

power_multiple_binary_GEE1_n100_3=power_multiple_binary_GEE1_3[power_multiple_binary_GEE1_3$n=='n=100',][Sopt_GEE1['Sopt_GEE1_n100',3],]
power_multiple_binary_GEE1_n150_3=power_multiple_binary_GEE1_3[power_multiple_binary_GEE1_3$n=='n=150',][Sopt_GEE1['Sopt_GEE1_n150',3],]

Sopt_GEE1_multiple_binary_3=rbind(power_multiple_binary_GEE1_n100_3,
                                  power_multiple_binary_GEE1_n150_3)

Plot_power_multiple_binary_GEE1_3<-ggplot(power_multiple_binary_GEE1_3, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols2)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.5~','~rho[X]==rho[Y] ~'=0.8'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE1_multiple_binary_3, 
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
Plot_power_multiple_binary_GEE1_3

#rhoXY= 0.5,rhoX= 0.9;
#rhoXY= 0.6,rhoX= 0.7; rhoXY= 0.6,rhoX= 0.8; rhoXY= 0.6,rhoX= 0.9;
#rhoXY= 0.7,rhoX= 0.8; rhoXY= 0.7,rhoX= 0.9

#rhoXY= 0.5,rhoX= 0.9
power_multiple_binary_GEE1_4 <- data.frame(n=rep(c('n=100','n=150'), each=9),
                                           S=rep(1:9,2),
                                           power=c(power_multiple_binary_GEE1_n100[4,],
                                                   power_multiple_binary_GEE1_n150[4,]))
power_multiple_binary_GEE1_4$n=as.character(power_multiple_binary_GEE1_4$n)
power_multiple_binary_GEE1_4$n=factor(power_multiple_binary_GEE1_4$n, levels=unique(power_multiple_binary_GEE1_4$n))

power_multiple_binary_GEE1_4$S=as.character(power_multiple_binary_GEE1_4$S)
power_multiple_binary_GEE1_4$S=factor(power_multiple_binary_GEE1_4$S, levels=unique(power_multiple_binary_GEE1_4$S))
power_multiple_binary_GEE1_4

power_multiple_binary_GEE1_n100_4=power_multiple_binary_GEE1_4[power_multiple_binary_GEE1_4$n=='n=100',][Sopt_GEE1['Sopt_GEE1_n100',4],]
power_multiple_binary_GEE1_n150_4=power_multiple_binary_GEE1_4[power_multiple_binary_GEE1_4$n=='n=150',][Sopt_GEE1['Sopt_GEE1_n150',4],]

Sopt_GEE1_multiple_binary_4=rbind(power_multiple_binary_GEE1_n100_4,
                                  power_multiple_binary_GEE1_n150_4)

Plot_power_multiple_binary_GEE1_4<-ggplot(power_multiple_binary_GEE1_4, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols2)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.5~','~rho[X]==rho[Y] ~'=0.9'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE1_multiple_binary_4, 
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
Plot_power_multiple_binary_GEE1_4

#rhoXY= 0.6,rhoX= 0.7
power_multiple_binary_GEE1_5 <- data.frame(n=rep(c('n=100','n=150'), each=9),
                                           S=rep(1:9,2),
                                           power=c(power_multiple_binary_GEE1_n100[5,],
                                                   power_multiple_binary_GEE1_n150[5,]))
power_multiple_binary_GEE1_5$n=as.character(power_multiple_binary_GEE1_5$n)
power_multiple_binary_GEE1_5$n=factor(power_multiple_binary_GEE1_5$n, levels=unique(power_multiple_binary_GEE1_5$n))

power_multiple_binary_GEE1_5$S=as.character(power_multiple_binary_GEE1_5$S)
power_multiple_binary_GEE1_5$S=factor(power_multiple_binary_GEE1_5$S, levels=unique(power_multiple_binary_GEE1_5$S))
power_multiple_binary_GEE1_5

power_multiple_binary_GEE1_n100_5=power_multiple_binary_GEE1_5[power_multiple_binary_GEE1_5$n=='n=100',][Sopt_GEE1['Sopt_GEE1_n100',5],]
power_multiple_binary_GEE1_n150_5=power_multiple_binary_GEE1_5[power_multiple_binary_GEE1_5$n=='n=150',][Sopt_GEE1['Sopt_GEE1_n150',5],]

Sopt_GEE1_multiple_binary_5=rbind(power_multiple_binary_GEE1_n100_5,
                                  power_multiple_binary_GEE1_n150_5)

Plot_power_multiple_binary_GEE1_5<-ggplot(power_multiple_binary_GEE1_5, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols2)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.6~','~rho[X]==rho[Y] ~'=0.7'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE1_multiple_binary_5, 
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
Plot_power_multiple_binary_GEE1_5

#rhoXY= 0.6,rhoX= 0.8
power_multiple_binary_GEE1_6 <- data.frame(n=rep(c('n=100','n=150'), each=9),
                                           S=rep(1:9,2),
                                           power=c(power_multiple_binary_GEE1_n100[6,],
                                                   power_multiple_binary_GEE1_n150[6,]))
power_multiple_binary_GEE1_6$n=as.character(power_multiple_binary_GEE1_6$n)
power_multiple_binary_GEE1_6$n=factor(power_multiple_binary_GEE1_6$n, levels=unique(power_multiple_binary_GEE1_6$n))

power_multiple_binary_GEE1_6$S=as.character(power_multiple_binary_GEE1_6$S)
power_multiple_binary_GEE1_6$S=factor(power_multiple_binary_GEE1_6$S, levels=unique(power_multiple_binary_GEE1_6$S))
power_multiple_binary_GEE1_6

power_multiple_binary_GEE1_n100_6=power_multiple_binary_GEE1_6[power_multiple_binary_GEE1_6$n=='n=100',][Sopt_GEE1['Sopt_GEE1_n100',6],]
power_multiple_binary_GEE1_n150_6=power_multiple_binary_GEE1_6[power_multiple_binary_GEE1_6$n=='n=150',][Sopt_GEE1['Sopt_GEE1_n150',6],]

Sopt_GEE1_multiple_binary_6=rbind(power_multiple_binary_GEE1_n100_6,
                                  power_multiple_binary_GEE1_n150_6)

Plot_power_multiple_binary_GEE1_6<-ggplot(power_multiple_binary_GEE1_6, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols2)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.6~','~rho[X]==rho[Y] ~'=0.8'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE1_multiple_binary_6, 
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
Plot_power_multiple_binary_GEE1_6

#rhoXY= 0.6,rhoX= 0.9
power_multiple_binary_GEE1_7 <- data.frame(n=rep(c('n=100','n=150'), each=9),
                                           S=rep(1:9,2),
                                           power=c(power_multiple_binary_GEE1_n100[7,],
                                                   power_multiple_binary_GEE1_n150[7,]))
power_multiple_binary_GEE1_7$n=as.character(power_multiple_binary_GEE1_7$n)
power_multiple_binary_GEE1_7$n=factor(power_multiple_binary_GEE1_7$n, levels=unique(power_multiple_binary_GEE1_7$n))

power_multiple_binary_GEE1_7$S=as.character(power_multiple_binary_GEE1_7$S)
power_multiple_binary_GEE1_7$S=factor(power_multiple_binary_GEE1_7$S, levels=unique(power_multiple_binary_GEE1_7$S))
power_multiple_binary_GEE1_7

power_multiple_binary_GEE1_n100_7=power_multiple_binary_GEE1_7[power_multiple_binary_GEE1_7$n=='n=100',][Sopt_GEE1['Sopt_GEE1_n100',7],]
power_multiple_binary_GEE1_n150_7=power_multiple_binary_GEE1_7[power_multiple_binary_GEE1_7$n=='n=150',][Sopt_GEE1['Sopt_GEE1_n150',7],]

Sopt_GEE1_multiple_binary_7=rbind(power_multiple_binary_GEE1_n100_7,
                                  power_multiple_binary_GEE1_n150_7)

Plot_power_multiple_binary_GEE1_7<-ggplot(power_multiple_binary_GEE1_7, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols2)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.6~','~rho[X]==rho[Y] ~'=0.9'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE1_multiple_binary_7, 
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
Plot_power_multiple_binary_GEE1_7

#rhoXY= 0.7,rhoX= 0.8
power_multiple_binary_GEE1_8 <- data.frame(n=rep(c('n=100','n=150'), each=9),
                                           S=rep(1:9,2),
                                           power=c(power_multiple_binary_GEE1_n100[8,],
                                                   power_multiple_binary_GEE1_n150[8,]))
power_multiple_binary_GEE1_8$n=as.character(power_multiple_binary_GEE1_8$n)
power_multiple_binary_GEE1_8$n=factor(power_multiple_binary_GEE1_8$n, levels=unique(power_multiple_binary_GEE1_8$n))

power_multiple_binary_GEE1_8$S=as.character(power_multiple_binary_GEE1_8$S)
power_multiple_binary_GEE1_8$S=factor(power_multiple_binary_GEE1_8$S, levels=unique(power_multiple_binary_GEE1_8$S))
power_multiple_binary_GEE1_8

power_multiple_binary_GEE1_n100_8=power_multiple_binary_GEE1_8[power_multiple_binary_GEE1_8$n=='n=100',][Sopt_GEE1['Sopt_GEE1_n100',8],]
power_multiple_binary_GEE1_n150_8=power_multiple_binary_GEE1_8[power_multiple_binary_GEE1_8$n=='n=150',][Sopt_GEE1['Sopt_GEE1_n150',8],]

Sopt_GEE1_multiple_binary_8=rbind(power_multiple_binary_GEE1_n100_8,
                                  power_multiple_binary_GEE1_n150_8)

Plot_power_multiple_binary_GEE1_8<-ggplot(power_multiple_binary_GEE1_8, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols2)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.7~','~rho[X]==rho[Y] ~'=0.8'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE1_multiple_binary_8, 
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
Plot_power_multiple_binary_GEE1_8

#rhoXY= 0.7,rhoX= 0.9
power_multiple_binary_GEE1_9 <- data.frame(n=rep(c('n=100','n=150'), each=9),
                                           S=rep(1:9,2),
                                           power=c(power_multiple_binary_GEE1_n100[9,],
                                                   power_multiple_binary_GEE1_n150[9,]))
power_multiple_binary_GEE1_9$n=as.character(power_multiple_binary_GEE1_9$n)
power_multiple_binary_GEE1_9$n=factor(power_multiple_binary_GEE1_9$n, levels=unique(power_multiple_binary_GEE1_9$n))

power_multiple_binary_GEE1_9$S=as.character(power_multiple_binary_GEE1_9$S)
power_multiple_binary_GEE1_9$S=factor(power_multiple_binary_GEE1_9$S, levels=unique(power_multiple_binary_GEE1_9$S))
power_multiple_binary_GEE1_9

power_multiple_binary_GEE1_n100_9=power_multiple_binary_GEE1_9[power_multiple_binary_GEE1_9$n=='n=100',][Sopt_GEE1['Sopt_GEE1_n100',9],]
power_multiple_binary_GEE1_n150_9=power_multiple_binary_GEE1_9[power_multiple_binary_GEE1_9$n=='n=150',][Sopt_GEE1['Sopt_GEE1_n150',9],]

Sopt_GEE1_multiple_binary_9=rbind(power_multiple_binary_GEE1_n100_9,
                                  power_multiple_binary_GEE1_n150_9)

Plot_power_multiple_binary_GEE1_9<-ggplot(power_multiple_binary_GEE1_9, aes(x=S, y=power, group=n)) +
  geom_line(aes(color=n),size=2)+
  geom_point(aes(color=n),size=5)+ 
  scale_color_manual(values=cols2)+ 
  xlab('S')+ylab('Power')+
  labs(fill = "",title=bquote(~rho[XY]==0.7~','~rho[X]==rho[Y] ~'=0.9'))+  
  scale_y_continuous(breaks=seq(0,2,by=0.2),limits=c(0,1))+
  geom_point(data=Sopt_GEE1_multiple_binary_9, 
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
Plot_power_multiple_binary_GEE1_9

library(cowplot)
png(file=paste0('Plot_power_multiple_binary_GEE1.png'),width=2750,height=2700,units='px',pointsize =30)
ggdraw() +
  draw_label("GEE Model 1",fontface = 'bold',x=0.4,y=1,hjust = 0, vjust = 0, size=80) +
  draw_plot(Plot_power_multiple_binary_GEE1_1, x = 0, y =0.66, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE1_2, x = 0.33, y =0.66, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE1_3, x = 0.66, y =0.66, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE1_4, x = 0, y =0.33, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE1_5, x = 0.33, y =0.33, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE1_6, x = 0.66, y =0.33, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE1_7, x = 0, y =0, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE1_8, x = 0.33, y =0, width =0.33, height=0.33)+
  draw_plot(Plot_power_multiple_binary_GEE1_9, x = 0.66, y =0, width =0.33, height=0.33)+
  theme(plot.margin = margin(85, 0, 0, 0))
dev.off()