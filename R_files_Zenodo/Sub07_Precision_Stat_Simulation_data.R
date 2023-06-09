###### Library ######
library(diversitree)
library(Matrix)
library(stats4)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

###### Font preparation ######
subset(fonttable(), FamilyName == "Arial") ##Check not more than one "Arial" in your system
loadfonts(device = "postscript")
###### Input true value ######

load("fit_musse_neot_M0_YK2021.Robj")
true<-coef(fit.musse.neot)
true<-c(true,true[5]/true[6],true[4]/true[3])


###### Precision of MLE on Eurypterygii tree ######

res.coef<-c()
for(trial in 1:250){
  print(trial)
  filename<-paste("fit_musse_sim_neot_trial",trial,".Robj",sep="")
  if(file.exists(filename)){
    load(file=filename)
    res.coef<-rbind(res.coef,coef(fit.musse.sim))
  }
}
write.table(res.coef,file="res_coef_sim_neot815_XXXXXX.txt")

coefdata<-res.coef[1:100,]
colnames(coefdata)<-c("lambda","mu","k3","k4","k2","k1")
coefdata$kf<-coefdata$k2/coefdata$k1
coefdata$ki<-coefdata$k4/coefdata$k3
coefdata$id<-rownames(coefdata)

nordata<-data.frame(id=coefdata$id,
                    k1=coefdata$k1/true[6],k2=coefdata$k2/true[5],
                    k3=coefdata$k3/true[3],k4=coefdata$k4/true[4],
                    kf=coefdata$kf/true[7],ki=coefdata$ki/true[8])

mean_log10_dis<-function(vec)
  mean(abs(log10(vec)))

c(mean_log10_dis(nordata$k1),
  mean_log10_dis(nordata$k2),
  mean_log10_dis(nordata$k3),
  mean_log10_dis(nordata$k4))



mylabels<-c(expression(paste(italic('k'[1]))),
            expression(paste(italic('k'[2]))),
            expression(paste(italic('k'[3]))),
            expression(paste(italic('k'[4]))),
            expression(paste(italic('K'[f]))),
            expression(paste(italic('K'[i]))))
ylab<-expression(paste("log"[10]," fold difference from true value"))
g1<-nordata %>%
  pivot_longer(!id) %>%
  ggplot(aes(x=name,y=log10(value)))+
  geom_hline(yintercept=0)+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=0.02)+
  coord_cartesian(ylim=c(-2,2))+
  scale_x_discrete(labels=mylabels)+
  theme(text=element_text(size=9),
        panel.border=element_rect(colour="black",size=0.6,fill=NA),
        axis.title.y=element_text(size=10),
        plot.title=element_text(size=9))+
   labs(x="Parameters",y=ylab,title="Eurypterygii tree")
#ggsave(g1,file="Precision_of_estimates_neot815.pdf",device="pdf",width=8,height=5)


###### Precision of MLE on Cyprinodontiformes tree ######

res.coef<-c()
for(trial in 1:250){
  print(trial)
  filename<-paste("fit_musse_sim_Cyprinodon_trial",trial,".Robj",sep="")
  if(file.exists(filename)){
    load(file=filename)
    res.coef<-rbind(res.coef,coef(fit.musse.sim))
  }
}
write.table(res.coef,file="res_coef_sim_Cyprinodon_XXXXXX.txt")

coefdata<-res.coef[1:100,]
colnames(coefdata)<-c("lambda","mu","k3","k4","k2","k1")
coefdata$kf<-coefdata$k2/coefdata$k1
coefdata$ki<-coefdata$k4/coefdata$k3
coefdata$id<-rownames(coefdata)

nordata<-data.frame(id=coefdata$id,
                    k1=coefdata$k1/true[6],k2=coefdata$k2/true[5],
                    k3=coefdata$k3/true[3],k4=coefdata$k4/true[4],
                    kf=coefdata$kf/true[7],ki=coefdata$ki/true[8])

mean_log10_dis<-function(vec)
  mean(abs(log10(vec)))

c(mean_log10_dis(nordata$k1),
  mean_log10_dis(nordata$k2),
  mean_log10_dis(nordata$k3),
  mean_log10_dis(nordata$k4))

g2<-nordata %>%
  pivot_longer(!id) %>%
  ggplot(aes(x=name,y=log10(value)))+
  geom_hline(yintercept=0)+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=0.02)+
  coord_cartesian(ylim=c(-2,2))+
  scale_x_discrete(labels=mylabels)+
  theme(text=element_text(size=9),
        panel.border=element_rect(colour="black",size=0.6,fill=NA),
        plot.title=element_text(size=9))+
  labs(x="Parameters",y=NULL,title="Cyprinodontifomes tree")


###### Precision of MLE on Goodeidae tree ######

res.coef<-c()
for(trial in 1:250){
  print(trial)
  filename<-paste("fit_musse_sim_Goodeidae_trial",trial,".Robj",sep="")
  if(file.exists(filename)){
    load(file=filename)
    res.coef<-rbind(res.coef,coef(fit.musse.sim))
  }
}
write.table(res.coef,file="res_coef_sim_Goodeidae_XXXXXX.txt")


coefdata<-res.coef[1:100,]
colnames(coefdata)<-c("lambda","mu","k3","k4","k2","k1")
coefdata$kf<-coefdata$k2/coefdata$k1
coefdata$ki<-coefdata$k4/coefdata$k3
coefdata$id<-rownames(coefdata)

nordata<-data.frame(id=coefdata$id,
                    k1=coefdata$k1/true[6],k2=coefdata$k2/true[5],
                    k3=coefdata$k3/true[3],k4=coefdata$k4/true[4],
                    kf=coefdata$kf/true[7],ki=coefdata$ki/true[8])

mean_log10_dis<-function(vec)
  mean(abs(log10(vec)))

c(mean_log10_dis(nordata$k1),
  mean_log10_dis(nordata$k2),
  mean_log10_dis(nordata$k3),
  mean_log10_dis(nordata$k4))

g3<-nordata %>%
  pivot_longer(!id) %>%
  ggplot(aes(x=name,y=log10(value)))+
  geom_hline(yintercept=0)+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=0.02)+
  coord_cartesian(ylim=c(-2,2))+
  scale_x_discrete(labels=mylabels)+
  theme(text=element_text(size=9),
        panel.border=element_rect(colour="black",size=0.6,fill=NA),
        plot.title=element_text(size=9))+
  labs(x="Parameters",y=NULL,title="Goodeidae tree")

###### Plot of combined figure ######

postscript("Fig_precision_of_estimates_XXXXXX.eps", height = 3, width = 5.1,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
ggarrange(g1,g2,g3,nrow=1,widths=c(1,0.92,0.92))
dev.off()


