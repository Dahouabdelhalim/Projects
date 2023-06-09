#load the required library 
library(ggplot2)

##read in data file
#the input file was obtained using the "LD_by_population.sh" script; the output file from this script (LD_decay_20kb_all_populations.txt) was annotated to add the LD average and standard deviation across invasive populations

#plot LD decay (script was modified from https://www.biostars.org/p/300381/)
data_all<-read.csv(file="LD_decay_20kb_all_populations.csv",header=TRUE, sep = ',')
ggplot()+
  geom_ribbon(data=data_all,aes(x=mid,ymax=Invasive_average+Invasive_stdev,ymin=Invasive_average-Invasive_stdev), fill="lightgrey", alpha=0.5)+
  geom_line(data=data_all,aes(x=mid,y=SOR),size=0.7,alpha=1, colour="#e09148ff")+
  geom_line(data=data_all,aes(x=mid,y=MAR),size=0.7,alpha=1, colour="#e09148ff")+
  geom_line(data=data_all,aes(x=mid,y=LHA),size=0.7,alpha=1, colour="#e09148ff")+
  geom_line(data=data_all,aes(x=mid,y=LHB),size=0.7,alpha=1, colour="#e09148ff")+
  geom_line(data=data_all,aes(x=mid,y=GUA),size=0.7,alpha=1, colour="#e09148ff")+
  geom_line(data=data_all,aes(x=mid,y=JIC),size=0.7,alpha=1, colour="#b53e35ff")+
  geom_line(data=data_all,aes(x=mid,y=CAB),size=0.7,alpha=1, colour="#9dc8e9ff")+
  geom_line(data=data_all,aes(x=mid,y=ESM),size=0.7,alpha=1, colour="#415b9eff")+
  geom_line(data=data_all,aes(x=mid,y=POR),size=0.7,alpha=1, colour="#853789ff")+
  geom_line(data=data_all,aes(x=mid,y=Invasive_average),size=0.7,alpha=1, colour="black")+
  labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0,0.5*10^6,1*10^6,1.5*10^6,2*10^6),labels=c("0","0.5","1","1.5","2"))+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
       legend.position="none")
