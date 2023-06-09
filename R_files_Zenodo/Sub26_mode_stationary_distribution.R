###### Library ######
library(diversitree)
library(ggplot2)
library(ggpubr)
library(extrafont)


####### Drawing the 6 conditions for the mode of stationary distribution (Fig S5) ######

chrmax<-30
snum<-(chrmax+1)*(chrmax+2)/2-1
pdata<-data.frame(px=rep(0,length=snum),py=rep(0))
s<-1
for(d in 1:chrmax){
  for(j in d:(2*d)){
    pdata$px[s]<-j
    pdata$py[s]<-d
    s<-s+1
  }
}

mode_drawing<-function(kf,ki){
  #Corresponding (A.24)
  aline<-function(x){ ((1+ki)*x+ki)/(2+ki) }
  bline<-function(x){ ((1+ki)*x-1)/(2+ki) }
  cline<-function(x){  x/2+kf/ki-1/2  }
  dline<-function(x){  x-2*kf/(ki^2)+1 }
  eline<-function(x){  x/2+kf/ki }
  fline<-function(x){  x-2*kf/(ki^2) }
  bdata<-data.frame(bx=1:80)
  ggplot()+
    geom_point(aes(x=pdata$px,y=pdata$py),size=0.2)+
    scale_x_continuous(limit=c(0,80))+
    scale_y_continuous(limit=c(0,40))+
    stat_function(data=bdata,mapping=aes(bx),fun=aline,colour="dodgerblue2")+
    stat_function(aes(bx),fun=bline,colour="dodgerblue2")+
    stat_function(aes(bx),fun=cline,colour="red3")+
    stat_function(aes(bx),fun=eline,colour="red3")+
    stat_function(aes(bx),fun=dline,colour="green4")+
    stat_function(aes(bx),fun=fline,colour="green4")+
    annotate(geom="text",x=60,y=2,parse=T,label=paste0("italic(K[f]) == ",round(kf,3),"~italic(K[i]) == ",round(ki,3)))+
    theme(text=element_text(size=11))+
    labs(x="Arm number, x",y="Chromosome number, y")
}

gA<-mode_drawing(6.25,1)
gB<-mode_drawing(288,24) #Minimum Kf and Ki for (24, 24)
gC<-mode_drawing(8/22,4/22) #Maximum Kf and Ki for (47, 25)
gD<-mode_drawing(1/54,1/27) #Maximum Kf and Ki for (54, 27)

ggarrange(gA,gB,gC,gD,ncol=2,nrow=2,hjust=-0.2,vjust=1.5,
          font.label=list(size=20,familty="Helvetica Neue",face="plain"),labels=c("A","B","C","D"))
ggsave(filename="Fig_V3_Mode_stationary_distributions_XXXXXX.pdf",device="pdf",units="cm",width=20,height=15)


