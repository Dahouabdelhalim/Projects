#####################
# script for Hasik et al. 2021
#####################

# ENSURE THE FOLLOWING FILES ARE IN THE WORKING DIRECTORY

# Behavioral_Assays__By_Lake.csv
# lake.csv
# Lake_density_experiment_final.csv
# mesocosm.csv
# Mesocosm_assignments.csv
# PO_Assays_raw_data_2019_field_experiment.csv
# PO_Assays_Raw_Data_meso.csv
# Spring_2019_Behavioral_Assay_Distance_Totals.csv
# SummarizedData.csv

library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(grid)
library(gridExtra)
library(car)
library(plyr)
library(dplyr)
library(patchwork)
library(lsmeans)
library(cowplot)
library(lme4)
library(multcomp)
library(drc)
library(MuMIn)

# Load up color blind palette
cbbPalette <- c( "#56B4E9","#009E73", "#D55E00","#CC79A7", "#999999","#000000")

########
# LOAD DATA FILES
#######

#orginal mesocosm data
orig_meso<-read.csv("Mesocosm_assignments.csv")
#mesocosm po data
po_meso<-read.csv("PO_Assays_Raw_Data_meso.csv")
#add activity data
act.data<-read.csv("Spring_2019_Behavioral_Assay_Distance_Totals.csv",stringsAsFactors = T,
                   fileEncoding="UTF-8-BOM")

#original field data
orig_field<-read.csv("Lake_density_experiment_final.csv",
                     fileEncoding="UTF-8-BOM")
#field PO data
po_field<-read.csv("PO_Assays_raw_data_2019_field_experiment.csv")
#field activity data
act.field<-read.csv("Behavioral_Assays__By_Lake.csv")
#load growth rate data
meso.growth<-read.csv("mesocosm.csv",stringsAsFactors = T)
lake.growth<-read.csv("lake.csv")
#change names to match other df's
lake.growth<-lake.growth%>%mutate(
  Lake=case_when(
    Lake=="BobbKidd"~"BobbKiddLake",
    Lake=="Charleston"~"CharlestonLake",
    Lake=="Greenwood"~"GreenwoodLake",
    Lake=="Fayetteville"~"LakeFayetteville",
    Lake=="Wilson"~"LakeWilson",
    Lake=="Lock.and.Dam"~"LockAndDamPond")   
)

#merge mesocosm datasheets
total_meso<-merge(orig_meso,po_meso,by=c("Cage"))
total_meso<-total_meso %>% distinct(ID,.keep_all = TRUE) #remove duplicates from the dataframe
total_meso<-merge(total_meso,act.data,by=c("Individual"))

#merge field datasheets
total_field<-merge(orig_field,po_field,by=c("Lake","Cage"))
total_field<-total_field %>% distinct(ID,.keep_all = TRUE) #remove duplicates from the dataframe

#########################
# MESOCOSM DATA
#########################

#set selection levels and cage as factors
total_meso$Selection.level<-as.factor(total_meso$Selection.level)
total_meso$Cage<-as.factor(total_meso$Cage)

#get cage means of PO
po.cage.mean <- total_meso %>% group_by(Cage,Selection.level,Density) %>% summarise(cage.mean=mean(PO_mean),
                                                                                    cage.sd=sd(PO_mean),
                                                                                    cage.n=n_distinct(PO_mean))%>% 
  mutate(cage.se=cage.sd/sqrt(cage.n))

#get cage means of prot
prot.cage.mean <- total_meso %>% group_by(Cage,Selection.level,Density) %>% summarise(cage.mean.prot=mean(Prot_mean),
                                                                                      cage.sd.prot=sd(Prot_mean),
                                                                                    cage.n.prot=n_distinct(Prot_mean))%>% 
  mutate(cage.se.prot=cage.sd.prot/sqrt(cage.n.prot))

#merge cage po and prot means
cage.df<-merge(po.cage.mean,prot.cage.mean,by=c("Cage","Selection.level","Density")) 

#get cage means of distance
dist.cage.mean <- total_meso %>% group_by(Cage,Selection.level,Density) %>% summarise(cage.mean.dist=mean(Behavioral_Assay_Distance)) #create dataframe for the cage means

#merge cage po,prot,and distance means
cage.df<-merge(cage.df,dist.cage.mean,by=c("Cage","Selection.level","Density")) 

#create dataframe for the treatment means and se's
po.treat.mean <- po.cage.mean %>% 
  group_by(Selection.level,Density) %>% 
  summarise(treat.mean=mean(cage.mean),treat.sd=sd(cage.mean),treat.n=n_distinct(cage.mean))%>%
  mutate(treat.se=treat.sd/sqrt(treat.n))

#merge dataframes
cage.df<-merge(cage.df,po.treat.mean,by=c("Selection.level","Density"))
#set selection level and cage as factors
meso.growth$Selection.level<-as.factor(meso.growth$Selection.level)
meso.growth$Cage<-as.factor(meso.growth$Cage)
#merge dataframes
cage.df.growth<-merge(cage.df,meso.growth,by=c("Cage","Selection.level"))
cage.df.growth

#test for relationship between growth rates and total PO
cor.test(cage.df.growth$cage.mean,cage.df.growth$growth.rate,method = c("pearson", "kendall", "spearman")) #significant correlation of 0.26
meso.scatter<-ggscatter(cage.df.growth, x = "growth.rate", y = "cage.mean",
          add = "reg.line", conf.int = TRUE,     #visualize the relationship
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "growth rate", ylab = "cage mean po")

##mixed model
mod.mixed <- lmer(PO_mean ~ Prot_mean + Selection.level*Density + (1|Cage),
                  data = total_meso)
plot(mod.mixed)
summary(mod.mixed)
Anova(mod.mixed,type = "III")
anova(mod.mixed, type = 3)

#subset data by density
cage1<-subset(total_meso,total_meso$Density==1)
cage2<-subset(total_meso,total_meso$Density==2)
cage4<-subset(total_meso,total_meso$Density==4)
cage10<-subset(total_meso,total_meso$Density==10)
#run mixed models on individual density datasets
mod1<-glm(PO_mean ~ Prot_mean + Selection.level, data=cage1)
mod2<-lmer(PO_mean ~ Prot_mean + Selection.level + (1|Cage),  data=cage2)
mod4<-lmer(PO_mean ~ Prot_mean + Selection.level + (1|Cage),  data=cage4)
mod10<-lmer(PO_mean ~ Prot_mean + Selection.level + (1|Cage),  data=cage10)

lsmeans(mod1,
        ~ Selection.level,
        adjust="tukey")
lsmeans(mod2,
        ~ Selection.level,
        adjust="tukey")
lsmeans(mod4,
        ~ Selection.level,
        adjust="tukey")
lsmeans(mod10,
        ~ Selection.level,
        adjust="tukey")

#####
# no significant results for the mesocosm study, neither for main effects nor the interaction
#####

####
# SCRIPT TO CREATE FIG 2
####

#create dataframes for each level of selection
meso1<-subset(cage.df,cage.df$Selection.level=="1")
meso2<-subset(cage.df,cage.df$Selection.level=="2")
meso3<-subset(cage.df,cage.df$Selection.level=="3")
meso4<-subset(cage.df,cage.df$Selection.level=="4")

total_meso1<-subset(total_meso,total_meso$Selection.level=="1")
total_meso2<-subset(total_meso,total_meso$Selection.level=="2")
total_meso3<-subset(total_meso,total_meso$Selection.level=="3")
total_meso4<-subset(total_meso,total_meso$Selection.level=="4")

brewer.pal(n = 8, name = 'Blues')

meso1.plot<-ggplot(meso1,aes(x=Density,y=treat.mean,ymin=treat.mean-treat.se,
                             ymax=treat.mean+treat.se))+
  geom_point(data=total_meso1,aes(x=Density,y=PO_mean),size=1.5,shape=21,col="black",fill="#6BAED6",
             alpha=2/5,inherit.aes = F)+
  geom_errorbar(aes(ymin=treat.mean-treat.se, ymax=treat.mean+treat.se),col="#6BAED6",
                width=0,cex=1)+
  geom_smooth(data=total_meso1,aes(x=Density,y=PO_mean),
              method="lm",formula = y~x,se=TRUE,fullrange=TRUE,size=.5,color="black",inherit.aes = F)+
  geom_pointrange(shape=21,size=.75,col="#6BAED6",fill="#6BAED6")+
  theme_classic()+xlab("Damselfly density")+ylab("Mesocosm mean PO (\\u0394 od 485 nm / min)")+
  scale_y_continuous(limits=c(0,160),breaks=c(0,40,80,120,160),
                     labels=c(0,40,80,120,160))+
  scale_x_continuous(limits = c(0,11),
                     breaks=c(0,2.5,5.0,7.5,10.0),
                     labels=c(0,2.5,5.0,7.5,10.0))+
  theme(text = element_text(size=15),aspect.ratio = 0.6,legend.position = "none",
        axis.title.y = element_blank(),axis.title.x = element_blank()
        )

meso1.plot

meso2.plot<-ggplot(meso2,aes(x=Density,y=treat.mean,ymin=treat.mean-treat.se,
                             ymax=treat.mean+treat.se))+
  geom_point(data=total_meso2,aes(x=Density,y=PO_mean),size=1.5,shape=21,col="black",fill="#6BAED6",
             alpha=2/5,inherit.aes = F)+
  geom_errorbar(aes(ymin=treat.mean-treat.se, ymax=treat.mean+treat.se),col="#4292C6",
                width=0,cex=1)+
  geom_smooth(data=total_meso2,aes(x=Density,y=PO_mean),
              method="lm",formula = y~x,se=TRUE,fullrange=TRUE,size=.5,color="black",inherit.aes = F)+
  geom_pointrange(size=0.75,col="#4292C6")+
  theme_classic()+xlab("Damselfly density")+ylab("Mesocosm mean PO (\\u0394 od 485 nm / min)")+
  scale_y_continuous(limits=c(0,160),breaks=c(0,40,80,120,160),
                     labels=c(0,40,80,120,160))+
  scale_x_continuous(limits = c(0,11),
                     breaks=c(0,2.5,5.0,7.5,10.0),
                     labels=c(0,2.5,5.0,7.5,10.0))+
  theme(text = element_text(size=15),aspect.ratio = 0.6,legend.position = "none",
        axis.title.y = element_blank(),axis.title.x = element_blank()
  )

meso2.plot

meso3.plot<-ggplot(meso3,aes(x=Density,y=treat.mean,ymin=treat.mean-treat.se,
                             ymax=treat.mean+treat.se))+
  geom_point(data=total_meso3,aes(x=Density,y=PO_mean),size=1.5,shape=21,col="black",fill="#6BAED6",
             alpha=2/5,inherit.aes = F)+
  geom_errorbar(aes(ymin=treat.mean-treat.se, ymax=treat.mean+treat.se),col="#2171B5",
                width=0,cex=1)+
  geom_smooth(data=total_meso3,aes(x=Density,y=PO_mean),
              method="lm",formula = y~x,se=TRUE,fullrange=TRUE,size=.5,color="black",inherit.aes = F)+
  geom_pointrange(size=0.75,col="#2171B5")+
  theme_classic()+xlab("Damselfly density")+ylab("Mesocosm mean PO (\\u0394 od 485 nm / min)")+
  scale_y_continuous(limits=c(0,160),breaks=c(0,40,80,120,160),
                     labels=c(0,40,80,120,160))+
  scale_x_continuous(limits = c(0,11),
                     breaks=c(0,2.5,5.0,7.5,10.0),
                     labels=c(0,2.5,5.0,7.5,10.0))+
  theme(text = element_text(size=15),aspect.ratio = 0.6,legend.position = "none",
        axis.title.y = element_blank(),axis.title.x = element_blank()
  )

meso3.plot

meso4.plot<-ggplot(meso4,aes(x=Density,y=treat.mean,ymin=treat.mean-treat.se,
                             ymax=treat.mean+treat.se))+
  geom_point(data=total_meso4,aes(x=Density,y=PO_mean),size=1.5,shape=21,col="black",fill="#6BAED6",
             alpha=2/5,inherit.aes = F)+
  geom_errorbar(aes(ymin=treat.mean-treat.se, ymax=treat.mean+treat.se),color="#084594",width=0,cex=1)+
  geom_smooth(data=total_meso4,aes(x=Density,y=PO_mean),
              method="lm",formula = y~x,se=TRUE,fullrange=TRUE,size=.5,color="black",inherit.aes = F)+
  geom_pointrange(size=0.75,color="#084594")+
  theme_classic()+xlab("Damselfly density")+ylab("Mesocosm mean PO (\\u0394 od 485 nm / min)")+
  scale_y_continuous(limits=c(0,160),breaks=c(0,40,80,120,160),
                     labels=c(0,40,80,120,160))+
  scale_x_continuous(limits = c(0,11),
                     breaks=c(0,2.5,5.0,7.5,10.0),
                     labels=c(0,2.5,5.0,7.5,10.0))+
  theme(text = element_text(size=15),aspect.ratio = 0.6,legend.position = "none",
        axis.title.y = element_blank(),axis.title.x = element_blank()
  )

meso4.plot

all_meso<-grid.arrange(meso1.plot,meso2.plot,meso3.plot,meso4.plot,ncol=1)

#create common x and y labels

y.grob <- textGrob("Mesocosm mean PO (   od 485 nm / min)", 
                   gp=gpar(fontsize=15), rot=90)

x.grob <- textGrob("Damselfly density", 
                   gp=gpar(fontsize=15))

#add to plot
all_meso<-grid.arrange(arrangeGrob(all_meso, left = y.grob, bottom = x.grob))

##############
# MESO SCATTER-HISTOGRAM PLOT
##############

#check relationship between ativity and po
hist(cage.df$cage.mean) #bit of a left skew
hist(cage.df$cage.mean.dist) #right skew
cor.test(cage.df$cage.mean,cage.df$cage.mean.dist,
         method = c("pearson", "kendall", "spearman")) #non-significant correlation of 0.0006
ggscatter(cage.df, x = "cage.mean.dist", y = "cage.mean",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean Activity", ylab = "Cage PO mean")

#code for the individual parts
hist_top <- ggplot(data=total_meso)+geom_histogram(aes(x=Behavioral_Assay_Distance,color=Selection.level,
                                                       fill=Selection.level),binwidth = 15)+
  xlab("Individual activity rate (mm/3 hrs)")+ylab("Count")+
  theme_classic()+scale_color_manual(values=c("#6BAED6","#4292C6","#2171B5","#084594"),name="Activity level")+
  scale_fill_manual(values=c("#6BAED6","#4292C6","#2171B5","#084594",name="Activity level"))+theme(legend.position = "none",
                                                                                                   text = element_text(size=15))
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())
scatter <- ggplot()+geom_point(data=cage.df,aes(x=cage.mean.dist,y=cage.mean,col=Selection.level,
                                                fill=Selection.level),size=3)+
  ylab("Mesocosm mean PO (   od 485 nm / min)")+xlab("Mesocosm mean activity rate (mm/3 hrs)")+
  theme_classic()+scale_color_manual(values=c("#6BAED6","#4292C6","#2171B5","#084594"),name="Activity level")+
  scale_fill_manual(values=c("#6BAED6","#4292C6","#2171B5","#084594",name="Activity level"))+theme(legend.position = "none",
                                                                                                   text = element_text(size=15))
hist_right <- ggplot(data=total_meso)+geom_histogram(aes(x=PO_mean,color=Selection.level,
                                                         fill=Selection.level),binwidth = 3)+
  xlab("Individual mean PO (   od 485 nm / min)")+ylab("Count")+
  coord_flip()+theme_classic()+scale_color_manual(name="Activity level",
                                                  values=c("#6BAED6","#4292C6","#2171B5","#084594"))+
  scale_fill_manual(values=c("#6BAED6","#4292C6","#2171B5","#084594",name="Activity level"))+
  theme(legend.position = "right",text = element_text(size=15),legend.title = element_text("Activity level"))+guides(color=FALSE)

#put parts together
scatter_hist<-grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 2), heights=c(1, 4))

#new figure 2
fig.2<-grid.arrange(all_meso, scatter_hist, ncol=2, nrow=1, widths=c(1, 2))

####
# SCRIPT TO CREAT FIG. S1
####

#subset data
cage.df.1.only<-subset(cage.df,cage.df$Density==1)
#check relationship between ativity and po
hist(cage.df.1.only$cage.mean) #no skew
hist(cage.df.1.only$cage.mean.dist) #some outliers on the right
cor.test(cage.df.1.only$cage.mean,cage.df.1.only$cage.mean.dist,
         method = c("pearson", "kendall", "spearman")) #non-significant correlation of 0.0006
ggscatter(cage.df.1.only, x = "cage.mean.dist", y = "cage.mean",
          add = "reg.line", conf.int = TRUE) #visualize the relationship

single.scatter <- ggplot()+geom_point(data=cage.df.1.only,aes(x=cage.mean.dist,y=cage.mean,col=Selection.level,
                                                fill=Selection.level),size=3)+
  ylab("Damselfly mean PO (   od 485 nm / min)")+xlab("Damselfly activity rate (mm/3 hrs)")+
  theme_classic()+scale_color_manual(values=c("#6BAED6","#4292C6","#2171B5","#084594"),name="Activity level",
                                     labels=c("1","2","3","4"))+
  scale_fill_manual(values=c("#6BAED6","#4292C6","#2171B5","#084594",name="Activity level",
                             labels=c("1","2","3","4")))+
  theme(legend.position = "right",text = element_text(size=15),aspect.ratio = 1)+guides(fill=FALSE)
single.scatter

#########################
# FIELD DATA
#########################

#set cage as factor
total_field$Cage<-as.factor(total_field$Cage)

#get cage means
field.cage.mean <- total_field %>% group_by(Lake,Cage,Density) %>% summarise(cage.mean=mean(PO_mean),
                                                                             cage.sd=sd(PO_mean),
                                                                             cage.n=n_distinct(PO_mean))%>%
  mutate(cage.se=cage.mean/sqrt(cage.sd))

#get cage means of prot
prot.field.cage.mean <- total_field %>% group_by(Lake,Cage,Density) %>% summarise(cage.mean.prot=mean(Prot_mean),
                                                                                  cage.sd.prot=sd(Prot_mean),
                                                                                  cage.n.prot=n_distinct(Prot_mean))%>%
  mutate(cage.se.prot=cage.sd.prot/sqrt(cage.n.prot))

#merge cage po and prot means
field.df<-merge(field.cage.mean,prot.field.cage.mean,
                by=c("Lake","Density","Cage")) 
field.df.growth<-merge(field.df,lake.growth,by=c("Lake","Cage"))

#test for relationship between growth rates and total po
cor.test(field.df.growth$cage.mean,field.df.growth$growth.rate,method = c("pearson", "kendall", "spearman")) #significant correlation of 0.26
lake.scatter<-ggscatter(field.df.growth, x = "growth.rate", y = "cage.mean",
                        conf.int = FALSE,     #visualize the relationship
                        xlab = "growth rate", ylab = "cage mean po")
lake.scatter

#################
# create Fig. S2
#################
meso.scatter<-ggscatter(cage.df.growth, shape=21,x = "growth.rate", y = "cage.mean",
                        xlab = "Mesocosm mean growth rate / day", 
                        ylab = "Mesocosm mean total PO (\\u0394 od 485 nm / min)",
                        color="black",
                        fill = "grey")+xlim(.008,0.04)+ylim(16,129)+theme(text = element_text(size=10))

meso.scatter

lake.scatter<-ggscatter(field.df.growth, shape=21,x = "growth.rate", y = "cage.mean",
                        xlab = "Cage mean growth rate / day", 
                        ylab = "Cage mean total PO (\\u0394 od 485 nm / min)",
                        color="black",
                        fill = "grey")+xlim(.008,0.04)+ylim(16,129)+theme(text = element_text(size=10))
lake.scatter

fig.s2<-plot_grid(meso.scatter+theme(axis.text.x = element_blank(),
                                     text=element_text(size=10)),
                   NULL,
                   lake.scatter+theme(text=element_text(size=10)),
                   align="hv",
                   hjust=-1,
                   ncol = 1,
                   rel_heights=c(1,0,1))
fig.s2

#####################
# ADD ENVR DATA 
#####################
env.data<-read.csv("SummarizedData.csv",stringsAsFactors = T)
#get means for env data
env.cage.means<-merge(field.df,env.data,by="Lake") #merge data

#create dataframe for the lake means and se's
field.lake.mean <- field.df %>% 
  group_by(Lake) %>% 
  summarise(env.mean=mean(cage.mean),env.sd=sd(cage.mean),env.n=n_distinct(cage.mean))%>%
  mutate(env.se=env.sd/sqrt(env.n))
env.cage.means<-merge(env.cage.means,field.lake.mean,by=c("Lake"))

#create dataframe for the lake activity means
field.act.mean <- act.field %>% 
  group_by(Lake) %>% 
  summarise(activity.mean=mean(Activity.distance.mm))

env.cage.means<-merge(env.cage.means,field.act.mean,by=c("Lake"))

env.cage.means$Density <-as.numeric(as.character(env.cage.means$Density))

total_field<-merge(total_field,env.data,by="Lake") #merge data for mixed model

##mixed modeL with lake to see which lakes differed
mod.lake <- lmer(PO_mean ~ Prot_mean + Lake + (1|Lake:Cage),
                  data = total_field)
plot(mod.lake)
summary(mod.lake)
Anova(mod.lake,type = "III")
anova(mod.lake,type = 3)

## multcomp to see which lakes differed
library(emmeans)
emmeans(mod.lake,pairwise~Lake) #Lake Wilson lower than all others

##mixed model
mod.mixed <- lmer(PO_mean ~ Prot_mean +Prey.DensityL*Density+Fish.Densitym2*Density + (1|Lake:Cage),
                  data = total_field)
plot(mod.mixed)
summary(mod.mixed)
Anova(mod.mixed,type = "III")
anova(mod.mixed, type = 3)

##get R2

r.squaredGLMM(mod.mixed)

##############
# SCRIPT FOR FIG. S2
##############

#check relationship between ativity and po
just.lakes<-env.cage.means%>%
  group_by(Lake)%>%
  summarise(mean.po=mean(env.mean),
            mean.act=mean(activity.mean))
just.lakes
hist(just.lakes$mean.po) 
cor.test(just.lakes$mean.po,just.lakes$mean.act,
         method = c("pearson", "kendall", "spearman")) #non-significant correlation of -0.10
ggscatter(just.lakes, x = "mean.act", y = "mean.po",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mean Activity", ylab = "Lake PO mean")

#code for the individual parts
hist_top <- ggplot(data=act.field)+geom_histogram(aes(x=Activity.distance.mm,color=Lake,
                                                       fill=Lake),binwidth = 15)+
  theme_classic()+scale_color_manual(values=cbbPalette,name="Lake")+ylim(0,15)+
  xlab("Individual Activity Rate (mm/3 hrs)")+
  ylab("Count")+scale_fill_manual(values=cbbPalette,name="Lake",
                                  labels=c("Bobb Kidd", "Charleston", "Greenwood",
                                           "Fayetteville","Wilson","Lock and Dam"))+theme(legend.position = "none")
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())
scatter <- ggplot()+geom_point(data=env.cage.means,aes(x=activity.mean,y=env.mean,col=Lake,
                                                fill=Lake),size=3)+
  ylab("Lake mean PO (   od 485 nm / min)")+xlab("Lake Mean Activity Rate (mm/3 hrs)")+
  theme_classic()+scale_color_manual(values=cbbPalette,name="Lake")+
  scale_fill_manual(values=cbbPalette,name="Lake",
                    labels=c("Bobb Kidd", "Charleston", "Greenwood",
                             "Fayetteville","Wilson","Lock and Dam"))+theme(legend.position = "none")
hist_right <- ggplot(data=po_field)+geom_histogram(aes(x=PO_mean,color=Lake,
                                                         fill=Lake),binwidth = 3)+ylim(0,15)+
  xlab("Individual mean PO (   od 485 nm / min)")+ylab("Count")+coord_flip()+
  theme_classic()+scale_color_manual(values=cbbPalette,name="Lake",
                                     labels=c("Bobb Kidd", "Charleston", "Greenwood",
                                              "Fayetteville","Wilson","Lock and Dam"))+
  scale_fill_manual(values=cbbPalette,name="Lake",labels=c("Bobb Kidd", "Charleston", "Greenwood",
                                                           "Fayetteville","Wilson","Lock and Dam"))

#put parts together
scatter_hist_field<-grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 2), heights=c(1, 4))

######
# PREY
######

#NLS via asymptotic exponential with two terms with individual data 
prey.nls.2<-nls(PO_mean~a*(1-exp(-c*Prey.DensityL)),
                start=list(a=125,c=0.03730772),data=total_field,
                control=nls.control(printEval=FALSE, minFactor=2^-50, warnOnly=TRUE,maxiter=10000))
summary(prey.nls.2)
coef(prey.nls.2)
confint(prey.nls.2)


#plot it
x <- seq(0, max(total_field$Prey.DensityL), length=100)
y <- predict(prey.nls.2, list(Prey.DensityL=x))
plot(total_field$Prey.DensityL,total_field$PO_mean,
     xlim = c(0,130),ylim=c(0,130))
points(x, y, type='l', col='blue')
#calculate R2
sse<-as.vector((summary(prey.nls.2)[[3]])^2*76)
sse
null <- lm(total_field$PO_mean~1)
sst <- as.vector(unlist(summary.aov(null)[[1]][2]))
sst
100*(sst-sse)/sst #R2 is 87.56

#get CI's for model
library(propagate)

#construct dataframe for 95ci's
df_total = data.frame()
for(i in 0:130){
  #run model
  mod<-predictNLS(model = prey.nls.2,newdata = data.frame(Prey.DensityL=i),
                  interval = "confidence")
  p<-mod$prop
  #extract relevant info
  mean<-p[[c(1,13,1)]]
  lwr<-p[[c(1,13,5)]]
  upr<-p[[c(1,13,6)]]
  # add vector to a dataframe
  df <- data.frame(cbind(mean,lwr,upr))
  df_total <- rbind(df_total,df)
}

df_x_total = data.frame()
for(i in 0:130){
  #run
  Prey.DensityL<-i
  # add vector to a dataframe
  df_x <- data.frame(Prey.DensityL)
  df_x_total <- rbind(df_x_total,df_x)
}

df_total<-cbind(df_total,df_x_total)


####
# SCRIPT FOR FIG. 3h
####

field.prey.plot<-ggplot(env.cage.means,aes(x=Prey.DensityL,y=env.mean,color=Lake,
                                           fill=Lake,ymin=env.mean-env.se,
                                           ymax=env.mean+env.se))+
  geom_point(data=total_field,aes(x=Prey.DensityL,y=PO_mean,fill=Lake),col="black",shape=21,size=1.5,
             alpha=3/5,inherit.aes = F)+
  geom_errorbar(aes(ymin=env.mean-env.se, ymax=env.mean+env.se,col=Lake),
                width=0,cex=1)+
  geom_pointrange(aes(color=Lake),size=0.75)+
  geom_ribbon(data=df_total,aes(x=Prey.DensityL,y=mean,ymin=lwr, ymax=upr), alpha=0.2,inherit.aes = F)+
  geom_smooth(data=env.cage.means,aes(x=Prey.DensityL,y=cage.mean),
              method="nls",formula = y~a*(1-exp(-c*x)),
              method.args = list(start=c(a=125,c=0.03730772)),
              se=FALSE,fullrange=TRUE,size=1,color="black",inherit.aes = F)+
  ylab("Cage mean PO (   od 485 nm / min)")+xlab("Prey density / L")+
  theme_classic()+
  scale_color_manual(values=cbbPalette,name="Lake",
                                     labels=c("Bobb Kidd Lake", "Charleston Lake", "Greenwood Lake",
                                              "Lake Fayetteville","Lake Wilson","Lock and Dam Pond"))+
  scale_fill_manual(values=cbbPalette,name="Lake",
                    labels=c("Bobb Kidd Lake", "Charleston Lake", "Greenwood Lake",
                             "Lake Fayetteville","Lake Wilson","Lock and Dam Pond"))+
  scale_y_continuous(breaks=c(30,60,90,120))+
  theme(text = element_text(size=15),aspect.ratio = 1)+
  # labs(tag = "c")+
  guides(fill = guide_legend(override.aes= list(alpha = 0.6)))

field.prey.plot

###########
# FISH
###########

####
# SCRIPT FOR FIG. 3G
####

field.fish.plot<-ggplot(env.cage.means,aes(x=Fish.Densitym2,y=env.mean,color=Lake,
                                           fill=Lake,ymin=env.mean-env.se,
                                           ymax=env.mean+env.se))+
  geom_point(data=total_field,aes(x=Fish.Densitym2,y=PO_mean,fill=Lake),col="black",shape=21,size=1.5,
             alpha=3/5,inherit.aes = F)+
  geom_errorbar(aes(ymin=env.mean-env.se, ymax=env.mean+env.se,col=Lake),
                width=0,cex=1)+
  geom_pointrange(aes(color=Lake),size=0.75)+
  geom_smooth(data=env.cage.means,aes(x=Fish.Densitym2,y=cage.mean),
              method="lm",formula = y~x,se=TRUE,fullrange=TRUE,size=1,color="black",inherit.aes = F)+
  ylab("Cage mean PO (   od 485 nm / min)")+xlab(bquote('Fish density /'~m^2))+
  theme_classic()+scale_color_manual(values=cbbPalette,name="Lake",
                                     labels=c("Bobb Kidd", "Charleston", "Greenwood",
                                              "Fayetteville","Wilson","Lock and Dam"))+
  scale_fill_manual(values=cbbPalette,name="Lake",
                    labels=c("Bobb Kidd", "Charleston", "Greenwood",
                             "Fayetteville","Wilson","Lock and Dam"))+
  scale_y_continuous(breaks=c(30,60,90,120))+
  theme(text = element_text(size=15),aspect.ratio = 1)+
  # labs(tag = "b")+
  guides(fill = guide_legend(override.aes= list(alpha = 0.6)))

field.fish.plot

############
# split dataframe up by lake to make individual density plots
############
envBK<-subset(env.cage.means,env.cage.means$Lake=="BobbKiddLake")
envC<-subset(env.cage.means,env.cage.means$Lake=="CharlestonLake")
envG<-subset(env.cage.means,env.cage.means$Lake=="GreenwoodLake")
envF<-subset(env.cage.means,env.cage.means$Lake=="LakeFayetteville")
envW<-subset(env.cage.means,env.cage.means$Lake=="LakeWilson")
envLD<-subset(env.cage.means,env.cage.means$Lake=="LockAndDamPond")

total_fieldBK<-subset(total_field,total_field$Lake=="BobbKiddLake")
total_fieldC<-subset(total_field,total_field$Lake=="CharlestonLake")
total_fieldG<-subset(total_field,total_field$Lake=="GreenwoodLake")
total_fieldF<-subset(total_field,total_field$Lake=="LakeFayetteville")
total_fieldW<-subset(total_field,total_field$Lake=="LakeWilson")
total_fieldLD<-subset(total_field,total_field$Lake=="LockAndDamPond")

############
# SCRIPT FOR FIGS. 3A-3F
############
bobbkidd.density.plot<-ggplot(envBK,aes(x=Density,y=env.mean,ymin=env.mean-env.se,
                                           ymax=env.mean+env.se))+
  geom_point(data=total_fieldBK,aes(x=Density,y=PO_mean),size=1.5,shape=21,col="black",fill="#56B4E9",
             alpha=3/5,inherit.aes = F)+
  geom_errorbar(col="#56B4E9",aes(ymin=env.mean-env.se, ymax=env.mean+env.se),
                width=0,cex=1)+
  geom_pointrange(col="#56B4E9",size=0.5)+
  geom_smooth(data=total_fieldBK,aes(x=Density,y=PO_mean),
              method="lm",formula = y~x,se=TRUE,fullrange=TRUE,size=1,color="black",inherit.aes = F)+
  ylab("Cage mean PO (\\u0394 od 485 nm / min)")+xlab("Damselfly density")+
  theme_classic()+
  scale_y_continuous(limits=c(0,160),breaks=c(0,40,80,120,160),
                     labels=c(0,40,80,120,160))+
  scale_x_continuous(limits = c(0,11),
                     breaks=c(0,2.5,5.0,7.5,10.0),
                     labels=c(0,2.5,5.0,7.5,10.0))+
  theme(text = element_text(size=15),aspect.ratio = 0.6,legend.position = "none",
         axis.title.y = element_blank(),axis.title.x = element_blank())

bobbkidd.density.plot

charleston.density.plot<-ggplot(envC,aes(x=Density,y=env.mean,ymin=env.mean-env.se,
                                          ymax=env.mean+env.se))+
  geom_point(data=total_fieldC,aes(x=Density,y=PO_mean),size=1.5,shape=21,col="black",fill="#009E73",
             alpha=3/5,inherit.aes = F)+
  geom_errorbar(col="#009E73",aes(ymin=env.mean-env.se, ymax=env.mean+env.se),
                width=0,cex=1)+
  geom_pointrange(col="#009E73",size=0.5)+
  geom_smooth(data=total_fieldC,aes(x=Density,y=PO_mean),
              method="lm",formula = y~x,se=TRUE,fullrange=TRUE,size=1,color="black",inherit.aes = F)+
  ylab("Cage mean PO (\\u0394 od 485 nm / min)")+xlab("Damselfly density")+
  theme_classic()+
  scale_y_continuous(limits=c(0,160),breaks=c(0,40,80,120,160),
                     labels=c(0,40,80,120,160))+
  scale_x_continuous(limits = c(0,11),
                     breaks=c(0,2.5,5.0,7.5,10.0),
                     labels=c(0,2.5,5.0,7.5,10.0))+
  theme(text = element_text(size=15),aspect.ratio = 0.6,legend.position = "none",
        axis.title.y = element_blank(),axis.title.x = element_blank())

charleston.density.plot

greenwood.density.plot<-ggplot(envG,aes(x=Density,y=env.mean,ymin=env.mean-env.se,
                                         ymax=env.mean+env.se))+
  geom_point(data=total_fieldG,aes(x=Density,y=PO_mean),size=1.5,shape=21,col="black",fill="#D55E00",
             alpha=3/5,inherit.aes = F)+
  geom_errorbar(col="#D55E00",aes(ymin=env.mean-env.se, ymax=env.mean+env.se),
                width=0,cex=1)+
  geom_pointrange(col="#D55E00",size=0.5)+
  geom_smooth(data=total_fieldG,aes(x=Density,y=PO_mean),
              method="lm",formula = y~x,se=TRUE,fullrange=TRUE,size=1,color="black",inherit.aes = F)+
  ylab("Cage mean PO (\\u0394 od 485 nm / min)")+xlab("Damselfly density")+
  theme_classic()+
  scale_y_continuous(limits=c(0,160),breaks=c(0,40,80,120,160),
                     labels=c(0,40,80,120,160))+
  scale_x_continuous(limits = c(0,11),
                     breaks=c(0,2.5,5.0,7.5,10.0),
                     labels=c(0,2.5,5.0,7.5,10.0))+
  theme(text = element_text(size=15),aspect.ratio = 0.6,legend.position = "none",
        axis.title.y = element_blank(),axis.title.x = element_blank())

greenwood.density.plot

fayetteville.density.plot<-ggplot(envF,aes(x=Density,y=env.mean,ymin=env.mean-env.se,
                                        ymax=env.mean+env.se))+
  geom_point(data=total_fieldF,aes(x=Density,y=PO_mean),size=1.5,shape=21,col="black",fill="#CC79A7",
             alpha=3/5,inherit.aes = F)+
  geom_errorbar(col="#CC79A7",aes(ymin=env.mean-env.se, ymax=env.mean+env.se),
                width=0,cex=1)+
  geom_pointrange(col="#CC79A7",size=0.5)+
  geom_smooth(data=total_fieldF,aes(x=Density,y=PO_mean),
              method="lm",formula = y~x,se=TRUE,fullrange=TRUE,size=1,color="black",inherit.aes = F)+
  ylab("Cage mean PO (\\u0394 od 485 nm / min)")+xlab("Damselfly density")+
  theme_classic()+
  scale_y_continuous(limits=c(0,160),breaks=c(0,40,80,120,160),
                     labels=c(0,40,80,120,160))+
  scale_x_continuous(limits = c(0,11),
                     breaks=c(0,2.5,5.0,7.5,10.0),
                     labels=c(0,2.5,5.0,7.5,10.0))+
  theme(text = element_text(size=15),aspect.ratio = 0.6,legend.position = "none",
        axis.title.y = element_blank(),axis.title.x = element_blank())

fayetteville.density.plot

wilson.density.plot<-ggplot(envW,aes(x=Density,y=env.mean,ymin=env.mean-env.se,
                                           ymax=env.mean+env.se))+
  geom_point(data=total_fieldW,aes(x=Density,y=PO_mean),size=1.5,shape=21,col="black",fill="#999999",
             alpha=3/5,inherit.aes = F)+
  geom_errorbar(col="#999999",aes(ymin=env.mean-env.se, ymax=env.mean+env.se),
                width=0,cex=1)+
  geom_pointrange(col="#999999",size=0.5)+
  geom_smooth(data=total_fieldW,aes(x=Density,y=PO_mean),
              method="lm",formula = y~x,se=TRUE,fullrange=TRUE,size=1,color="black",inherit.aes = F)+
  ylab("Cage mean PO (\\u0394 od 485 nm / min)")+xlab("Damselfly density")+
  theme_classic()+
  scale_y_continuous(limits=c(0,160),breaks=c(0,40,80,120,160),
                     labels=c(0,40,80,120,160))+
  scale_x_continuous(limits = c(0,11),
                     breaks=c(0,2.5,5.0,7.5,10.0),
                     labels=c(0,2.5,5.0,7.5,10.0))+
  theme(text = element_text(size=15),aspect.ratio = 0.6,legend.position = "none",
        axis.title.y = element_blank(),axis.title.x = element_blank())

wilson.density.plot

ld.density.plot<-ggplot(envLD,aes(x=Density,y=env.mean,ymin=env.mean-env.se,
                                 ymax=env.mean+env.se))+
  geom_point(data=total_fieldLD,aes(x=Density,y=PO_mean),size=1.5,shape=21,col="black",fill="#000000",
             alpha=3/5,inherit.aes = F)+
  geom_errorbar(col="#000000",aes(ymin=env.mean-env.se, ymax=env.mean+env.se),
                width=0,cex=1)+
  geom_pointrange(col="#000000",size=0.5)+
  geom_smooth(data=total_fieldLD,aes(x=Density,y=PO_mean),
              method="lm",formula = y~x,se=TRUE,fullrange=TRUE,size=1,color="black",inherit.aes = F)+
  ylab("Cage mean PO (\\u0394 od 485 nm / min)")+xlab("Damselfly density")+
  theme_classic()+
  scale_y_continuous(limits=c(0,160),breaks=c(0,40,80,120,160),
                     labels=c(0,40,80,120,160))+
  scale_x_continuous(limits = c(0,11),
                     breaks=c(0,2.5,5.0,7.5,10.0),
                     labels=c(0,2.5,5.0,7.5,10.0))+
  theme(text = element_text(size=15),aspect.ratio = 0.6,legend.position = "none",
        axis.title.y = element_blank(),axis.title.x = element_blank())

ld.density.plot

#######
# ARRANGE FIGURE 3

# arrange fish and prey in a single row
fish_prey <- plot_grid(
  field.fish.plot + theme(legend.position="none")+guides(fill = guide_legend(override.aes= list(alpha = 0.6))),
  field.prey.plot + theme(legend.position="none")+guides(fill = guide_legend(override.aes= list(alpha = 0.6))),
  align = 'hv',
  hjust = -1,
  nrow = 2,rel_widths = c(1,1)
)
fish_prey

#extract legend
legend<-get_legend(field.fish.plot+theme(legend.position = "right")+
                     guides(fill = guide_legend(override.aes= list(alpha = 0.6))))

# arrange lakes in a single row
lakerow <- plot_grid(
  bobbkidd.density.plot,
  charleston.density.plot,
  greenwood.density.plot,
  fayetteville.density.plot,
  wilson.density.plot,
  ld.density.plot,
  align = 'v',
  hjust = -.5,
  scale = 1,
  nrow = 6,rel_heights = c(1,1,1,1,1,1)
)
lakerow 
#create common x and y labels

y.grob <- textGrob("Cage mean PO (   od 485 nm / min)", 
                   gp=gpar(fontsize=15), rot=90)

x.grob <- textGrob("Damselfly density", 
                   gp=gpar(fontsize=15))

#add to plot
lakerow<-grid.arrange(arrangeGrob(lakerow, left = y.grob, bottom = x.grob))

# arrange lakes in a single row
lakerow_fishprey <- plot_grid(
  lakerow,fish_prey,legend,
  axis = "l",
  align = 'h',
  hjust = -.8,
  nrow = 1,rel_widths = c(.7,1,.35)
)
lakerow_fishprey