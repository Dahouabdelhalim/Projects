###### Library ######

library(extrafont)
library(dplyr)
library(ggplot2)

###### Font preparation ######

subset(fonttable(), FamilyName == "Arial") ##Check not more than one "Arial" in your system
loadfonts(device = "postscript")

###### Data Preparation ######

rdata<-read.table(file="res_table_simulation_CM8_YK2021.txt",header=F)
colnames(rdata)<-c("Trial","Truek1","Truek2","Truek3","Truek4",
                   "Mknk1","Mknk2","Mknk3","Mknk4","MuSk1","MuSk2","MuSk3","MuSk4",
                   "MuSlamda","MuSmu","Mkntime","MuStime")

rdata<-rdata[order(rdata$Trial,decreasing=F)[1:100],]

pdata0<-data.frame(Comparison=rep("TrvMk"),
                   Parameter=rep(c("k1","k2","k3","k4"),each=100),
                  X=c(rdata$Truek1,rdata$Truek2,rdata$Truek3,rdata$Truek4),
                  Y=c(rdata$Mknk1,rdata$Mknk2,rdata$Mknk3,rdata$Mknk4))

pdata1<-data.frame(Comparison=rep("TrvMu"),
                   Parameter=rep(c("k1","k2","k3","k4"),each=100),
                  X=c(rdata$Truek1,rdata$Truek2,rdata$Truek3,rdata$Truek4),
                  Y=c(rdata$MuSk1,rdata$MuSk2,rdata$MuSk3,rdata$MuSk4))
                  
pdata2<-data.frame(Comparison=rep("MkvMu"),
                   Parameter=rep(c("k1","k2","k3","k4"),each=100),
                  X=c(rdata$Mknk1,rdata$Mknk2,rdata$Mknk3,rdata$Mknk4),
                  Y=c(rdata$MuSk1,rdata$MuSk2,rdata$MuSk3,rdata$MuSk4))


pdata<-rbind(pdata0,pdata1,pdata2)

pdata$Parameter<-factor(pdata$Parameter)
levels(pdata$Parameter)=c("k1"=expression(paste(italic('k'[1]))),
                          "k2"=expression(paste(italic('k'[2]))),
                          "k3"=expression(paste(italic('k'[3]))),
                          "k4"=expression(paste(italic('k'[4]))))
pdata$Comparison<-factor(pdata$Comparison,level=c("TrvMu","TrvMk","MkvMu"))
levels(pdata$Comparison)=c("TrvMu"=expression(paste("True (X) ",italic('vs.')," MuSSE (Y)")),
                           "TrvMk"=expression(paste("True (X) ",italic('vs.')," Mk-n (Y)")),
                           "MkvMu"=expression(paste("Mk-n (X) ",italic('vs.')," MuSSE (Y)")))

##### Correlation coefficient #######

tdata<-pdata %>%
  group_by(Comparison,Parameter) %>%
  summarise(Cor=cor(log10(X),log10(Y),method="pearson"))
tdata$Exp<-paste("italic('r')","'='",sprintf("%.2f",round(tdata$Cor,2)),sep="~")

###### Figure drawing ######

postscript("Fig_CM8_simulation_XXXXXX.eps", height = 6.8, width = 5.1,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)

ggplot()+
  geom_point(data=pdata,aes(x=log10(X),y=log10(Y)),size=0.3)+
  geom_text(data=tdata,aes(x=-1.5,y=-4.5,label=Exp),parse=T)+
  facet_grid(Parameter~Comparison,labeller=label_parsed)+
  coord_cartesian(xlim=c(-5,0),ylim=c(-5,0))+
  labs(x=expression(paste("log"[10],italic(" k"[i])," of X")),y=expression(paste("log"[10],italic(" k"[i])," of Y")))

dev.off()