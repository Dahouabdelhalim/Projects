#########
##Fig1D##
#########
#glade_NDVI_means is glade_NDVI_data.csv
par(mar=c(4.5,4.5,1,1))
plot(0,0,xlim=c(-25,250),ylim=c(0.475,1.025),type="n",xlab="Distance from Glade Edge (m)",ylab="NDVI",xaxs="i",yaxs="i",axes=FALSE,cex.lab=1.5)
abline(v=0,col="grey",lty=2)
for(i in 1:26){
  lines(NDVI_mean~Glade_Dist,data=Glade_NDVI_means[which(Glade_NDVI_means$Glade_ID==glades[[i]]),],col=Glade_Y[which(Glade_Y$Glade_ID==glades[[i]]),5])
}
axis(1,at=seq(0,250,50),cex.axis=1.25)
axis(2,at=seq(0.5,1,0.1),las=2,cex.axis=1.25)
abline(h=0.475)
abline(v=-25)
lines(NDVI~Dist,data=NDVI_predict,lwd=3)
##############################
##FigS2 and Null Model Tests##
##############################
library(ggplot2)
library(dplyr)
library(plyr)
library(gridExtra)
#
load(FigS2.Rdata) ##Objects: mean_mound_Fij, Mound_dists
muS2A <- ddply(mean_mound_Fij, "Glade_Colony", summarise, grp.mean=mean(mean_Fij))
muS2B <- ddply(Mound_dists[c(which(Mound_dists$Dist_Factor%in%c("Neighbors","<150m"))),], "Dist_Factor", summarise, grp.mean=mean(D))
S2A_2 =mean_mound_Fij %>%
  ggplot( aes(x=mean_Fij, fill=Glade_Colony)) +
  geom_histogram(color="#e9ecef",size=0, alpha=0.75, position = 'dodge',aes(y=..density..),binwidth = 0.05)+
  scale_fill_manual(values=c("dodgerblue", "#404080","#69b3a2"), labels=c("< 60m","Glade","Mound")) +
  labs(fill="", title = "", x = expression(Relatedness (F[ij])), y = "Density")+
  theme(legend.position = "top", legend.text = element_text(size=8),legend.key = element_rect(fill = "lightblue", color = NA),legend.key.size = unit(0.05, "in"),legend.key.width = unit(0.1,"in"),
        panel.background = element_blank(),
        axis.line = element_line(color ="black",size=0.25), axis.ticks = element_line(size=0.25),axis.title = element_text(size=10),axis.text = element_text(size=6))+
  scale_y_continuous(expand = c(0,NA))+ scale_x_continuous(expand = c(-0.125,NA))+
  geom_vline(data=muS2A,aes(xintercept=grp.mean),color=c("dodgerblue", "#404080","#69b3a2"),linetype="dashed",size=0.5)+ guides(color=FALSE)+
  theme(plot.margin = margin(0.,0.,0,0.1,"in"))

S2B_main=ggplot(Mound_dists, aes((geodist/1000),(D)))+geom_point(alpha=0.25)+
  geom_smooth(method = "gam",formula = y~s(x),color="red")+
  theme(axis.title = element_text(size=10),axis.text = element_text(size=6),
        legend.position="none",axis.line = element_line(colour = "black",size=0.25), axis.ticks = element_line(size=0.25),panel.background = element_blank())+
  labs(fill="", title = "", x = "Geographic Distance (km)", y = "Genetic Distance (D)",size=0.5)+
  theme(plot.margin = margin(0.35,0.,0.05,0.1,"in"))

S2B_inset<-ggplot(Mound_dists[which(Mound_dists$geodist<=150),], aes((geodist/1000),(D),color=neighbors))+geom_point(alpha=0.25,size=0.5)+
  scale_color_manual(values=c("red","black"))+
  geom_smooth(method = "gam",formula = y~s(x),color="red")+
  theme(axis.title = element_blank(),axis.text = element_text(size=5),
        legend.position="none",axis.line = element_line(colour = "black",size=0.25), axis.ticks = element_line(size=0.25),panel.background = element_blank(),plot.background = element_blank())+
  labs(fill="", title = "", x = "Geographic Distance (m)", y = "Genetic Distance (D)",size=0.5)

S2C<- Mound_dists[c(which(Mound_dists$Dist_Factor%in%c("Neighbors","<150m"))),] %>%
  ggplot(aes(x=D, fill=Dist_Factor)) +
  theme(panel.background = element_blank(),plot.background = element_blank(),panel.border = element_blank(),
        axis.title = element_text(size=10), axis.text = element_text(size=6),axis.line = element_line(colour = "black",size=0.25), axis.ticks = element_line(size=0.25),
        legend.position="top",legend.text = element_text(size=8),legend.key.size = unit(0.05, "in"),legend.key.width = unit(0.1,"in"),legend.title = element_blank())+
  geom_histogram(alpha=0.75, position = 'dodge',binwidth = 0.0025,size=0,aes(y=..density..)) +
  scale_fill_manual(values=c("dodgerblue", "pink"), labels=c("Non-Neighbors","Neighbors")) +
  labs(fill="", title = "", x = "Genetic Distance (D)", y = "Density")+
  scale_y_continuous(expand = c(0,NA))+ scale_x_continuous(expand = c(-0.125,NA))+
  geom_vline(data=muS2B,aes(xintercept=grp.mean), color=c("blue", "red"),linetype="dashed",size=0.5)+guides(color=FALSE)+
  theme(plot.margin = margin(0,0.05,0.05,0.1,"in"))

#FigS2
grid.arrange(S2A_2,S2B_main + annotation_custom(ggplotGrob(S2B_inset),xmin=2.5,ymin=0.007,xmax=12.6,ymax=0.07),S2C,nrow=1)

###########################################
##Null Model Tests of termite relatedness##
###########################################
##Inter-individual##
n_mound<-length(which(mean_mound_Fij$Glade_Colony=="Mound"))
n_glade<-length(which(mean_mound_Fij$Glade_Colony=="Glade"))
n_other<-length(which(mean_mound_Fij$Glade_Colony=="0 - 60 m"))
##########################################
##Null Model Test of on versus off glade##
on_off_data<-mean_mound_Fij[-c(which(mean_mound_Fij$Glade_Colony=="Mound")),]
on_off_aov<-aov(mean_Fij~Glade_Colony,data=on_off_data)
summary(on_off_aov)

rand_Fs<-NULL
for(i in 1:1000){
  rand_sample_mound<-on_off_data[sample(1:(dim(on_off_data)[1]),n_glade+n_other,replace = FALSE),]
  rand_sample_mound$Rand_Factor<-sample(c(rep(1,n_glade),rep(2,n_other)),n_glade+n_other,replace=FALSE)
  temp_aov<-aov(mean_Fij~as.factor(Rand_Factor),data=rand_sample_mound)
  rand_Fs<-c(rand_Fs,summary(temp_aov)[[1]][[4]][[1]])
  rm(rand_sample_mound,temp_aov)
  if(i%in%seq(0,1000,100))print(i)
}

hist(rand_Fs,xlim=c(0,30),yaxs="i",ylim=c(0,800))
abline(v=summary(on_off_aov)[[1]][[4]][[1]],col="red")
length(which(rand_Fs>summary(on_off_aov)[[1]][[4]][[1]]))/1000

#########################################
##Null Model Test of mound versus glade##
mound_glade_data<-mean_mound_Fij[-c(which(mean_mound_Fij$Glade_Colony=="0 - 60 m")),]
mound_glade_aov<-aov(mean_Fij~Glade_Colony,data=mound_glade_data)
summary(mound_glade_aov)

rand_Fs_2<-NULL
for(i in 1:1000){
  rand_sample_mound<-mound_glade_data[sample(1:(dim(mound_glade_data)[1]),n_glade+n_mound,replace = FALSE),]
  rand_sample_mound$Rand_Factor<-sample(c(rep(1,n_glade),rep(2,n_mound)),n_glade+n_mound,replace=FALSE)
  temp_aov<-aov(mean_Fij~as.factor(Rand_Factor),data=rand_sample_mound)
  rand_Fs_2<-c(rand_Fs_2,summary(temp_aov)[[1]][[4]][[1]])
  if(i%in%seq(0,1000,100))print(i)
  rm(rand_sample_mound,temp_aov)
}

hist(rand_Fs_2)
abline(v=summary(mound_glade_aov)[[1]][[4]][[1]],col="red")
length(which(rand_Fs_2>summary(mound_glade_aov)[[1]][[4]][[1]]))/1000
#############################################
##Null Model Test of mound versus off-glade##
mound_offglade_data<-mean_mound_Fij[-c(which(mean_mound_Fij$Glade_Colony=="Glade")),]
mound_offglade_aov<-aov(mean_Fij~Glade_Colony,data=mound_offglade_data)
summary(mound_offglade_aov)

rand_Fs_3<-NULL
for(i in 1:1000){
  rand_sample_mound<-mound_offglade_data[sample(1:(dim(mound_offglade_data)[1]),n_other+n_mound,replace = FALSE),]
  rand_sample_mound$Rand_Factor<-sample(c(rep(1,n_other),rep(2,n_mound)),n_other+n_mound,replace=FALSE)
  temp_aov<-aov(mean_Fij~as.factor(Rand_Factor),data=rand_sample_mound)
  rand_Fs_3<-c(rand_Fs_3,summary(temp_aov)[[1]][[4]][[1]])
  if(i%in%seq(0,1000,100))print(i)
  rm(rand_sample_mound,temp_aov)
}

hist(rand_Fs_3)
abline(v=summary(mound_offglade_aov)[[1]][[4]][[1]],col="red")
length(which(rand_Fs_3>summary(mound_offglade_aov)[[1]][[4]][[1]]))/1000
###############
##Inter-mound##
D_n_neighbors<-length(which(Mound_dists$Dist_Factor=="Neighbors"))
non_neighbor_dist<-150
D_n_non_neighbors<-length(which(Mound_dists$neighbors=="NO"&Mound_dists$geodist<=non_neighbor_dist))

D_local_data<-Mound_dists[c(which(Mound_dists$geodist<=non_neighbor_dist)),]

D_local_aov<-aov(D~neighbors,data=D_local_data)
summary(D_local_aov)

D_rand_Fs<-NULL
for(i in 1:1000){
  rand_sample_mound<-D_local_data[sample(1:(dim(D_local_data)[1]),D_n_non_neighbors+D_n_neighbors,replace = FALSE),]
  rand_sample_mound$Rand_Factor<-sample(c(rep(1,D_n_non_neighbors),rep(2,D_n_neighbors)),D_n_non_neighbors+D_n_neighbors,replace=FALSE)
  temp_aov<-aov(D~as.factor(Rand_Factor),data=rand_sample_mound)
  D_rand_Fs<-c(D_rand_Fs,summary(temp_aov)[[1]][[4]][[1]])
  rm(rand_sample_mound,temp_aov)
  if(i%in%seq(0,1000,100))print(i)
}

hist(D_rand_Fs,xlim=c(0,30),yaxs="i",xaxs="i")
abline(v=summary(D_local_aov)[[1]][[4]][[1]],col="red",lwd=2)
length(which(D_rand_Fs>summary(D_local_aov)[[1]][[4]][[1]]))/1000

#########
##FigS3##
#########
