# setwd("C:/Users/Astaroph/Dropbox/Harvard/Research/Gary polarization manuscript work")
eyes<-read.table("Gary_ED_PD_rechecked_v4_final.csv",header=T,sep=",")
#eyes<-read.table("Gary_ED_PD_rechecked_v3_workon.csv",header=T,sep=",")
library(nlme)
library(lattice)
library(lme4)
library(car)
library(lsmeans)
library(pbkrtest)
library(multcomp)
library(car)
library(effects)
eyes=subset(eyes,Notes!="Exclude?")

library(ggplot2)
eyes$Elevation30=eyes$Elevation+30

# ###Graphical exploration of edge detectionand 5E polarization detection by sex and elevation
# graph1=ggplot(
#   data = eyes
#   , mapping = aes(
#     y =ED
#     , x =Elevation
#     , group=Sex, colour=Sex,fill=Sex
#   )
# )+
#   geom_point(
#     mapping = aes(
#       colour =Sex,shape=Sex,fill=Sex
#     )
#     , size=3.5,colour="black"
#   )+
#   geom_smooth(method='glm',
#     mapping = aes(
#       fill =Sex,colour=Sex
#     )
#     , size=1.1
#   )+
#   labs(
#     x = 'Elevation'
#     ,title= 'PD ommatidia over elevation'
#     , y = 'PD', 
#     colour="Sex"
#   )
# graph1
# graph1_2=graph1+scale_colour_manual(values=c("#fdae61","#2c7bb6"))+theme_bw()+
#   scale_fill_manual(values=c("#fdae61","#2c7bb6"))+scale_shape_manual(values=c(21,22))
# graph1_2
# 
# #Changing the linetypes of the treatments
# graph1_3=graph1_2+scale_linetype_manual(values=c("solid","solid"))
# graph1_3
# #Rescaling the axes
# graph1_4=graph1_3+scale_x_continuous(breaks=c(-30,0,30,60))
# graph1_4
# #Resizing and moving the text in the graph
# graph1_5=graph1_4+theme(legend.title = element_text(colour="Black", size=16, face="bold"))+
#   theme(legend.text = element_text(colour="Black", size = 16, face = "bold"))+
#   theme(axis.text=element_text(colour="Black", size=16, face="plain"))+
#   theme(axis.title=element_text(colour="Black", size=16, face="bold"))+
#   theme(plot.title = element_blank())+
#   theme(axis.title.y=element_text(vjust=1))+
#   theme(strip.text=element_text(colour="Black", size = 16, face = "plain"))
# graph1_5
# library(grid)
# graph1_6=graph1_5+theme(panel.spacing = unit(0.5, "lines"),
#                         legend.position="right")
# graph1_6
# 
# #remove gridlines and top/right borders and facets
# graph1_7=graph1_6+ theme(axis.line = element_line(colour = "black"),
#                          panel.grid.major = element_blank(),
#                          panel.grid.minor = element_blank(),
#                          panel.border=element_blank())+
#   theme(axis.line = element_line(color = 'black'))+
#   theme(strip.text=element_blank())+
#   theme(strip.background=element_blank())
# graph1_7

model1<-lmer(sqrt(Only_ED)~Elevation*Sex+
(1|Individual),data=eyes,REML=FALSE)
summary(model1)
Anova(model1)

model2<-lmer(sqrt(Only_ED)~Elevation+Sex+
(1|Individual),data=eyes,REML=FALSE)
summary(model2)
Anova(model2)
anova(model2,model1)


#no significant difference between the models, go with the simpler model2

##Now pbmodcomp, start with the full model and drop factors hierarchically 
model1<-lmer(sqrt(ED)~Elevation*Sex+
(1|Individual),data=eyes,REML=FALSE)
summary(model1)
Anova(model1)
##Dropping the interaction
model2<-lmer(sqrt(ED)~Elevation+Sex+
(1|Individual),data=eyes,REML=FALSE)
summary(model2)
Anova(model2)
anova(model1,model2)

##Now dropping each individual factor
model2.1=update(model2,~.-Elevation)
model2.2=update(model2,~.-Sex)


#enabling parallelization
library(parallel)
(nc <- detectCores())
## Create clusters
cl <- makeCluster(rep("localhost", nc-1))

#testing each model against the next largest model
A=PBmodcomp(model1,model2,nsim=10000,cl=cl)
A
B_1=PBmodcomp(model2,model2.1,nsim=10000,cl=cl)
B_1
B_2=PBmodcomp(model2,model2.2,nsim=10000,cl=cl)
B_2

stopCluster(cl)
options(scipen=999)
mctab=data.frame(rbind(A$test,B_1$test,B_2$test))
modcomptable=data.frame(Model=c("model2","model2","model2.1","model2.1","model2.2",
"model2.2"),Test=c("LRT","PBtest","LRT","PBtest",
"LRT","PBtest"),Statistic=mctab$stat,P_value=mctab$p.value,
Factor=c("Elevation:Sex","Elevation:Sex","Elevation","Elevation","Sex","Sex"),
Significance=c("ns","ns","***","***","**","*"))
modcomptable
write.csv(modcomptable,"PBmodcomp_factor_significance_sqrt_Gary_ONLY_ED_final.csv")



#########Using the effects package to generate fitted values and standard error estimates for
##plotting
library(effects)

model1<-lmer(sqrt(ED)~Elevation*Sex+
               (1|Individual),data=eyes,REML=FALSE)
summary(model1)
Anova(model1)
##Dropping the interaction
model2<-lmer(sqrt(ED)~Elevation+Sex+
               (1|Individual),data=eyes,REML=FALSE)
summary(model2)
Anova(model2)
anova(model1,model2)

A=effect("Sex",model1, KR=TRUE,confidence.level=.95)
plot(A)
Sex_frame=as.data.frame(A)
Sex_frame

##Sex*Elevation
B=effect("Elevation:Sex",model1, KR=TRUE,confidence.level=.95,xlevels=
           list(Elevation=c(-30,0,30,60)),transformation=list(link=sqrt,inverse=function(x) x^2))
plot(B)
Elevation_sex_frame=as.data.frame(B)
Elevation_sex_frame

##################now using ggplot2 to plot these graphs
##
library(ggplot2)

##Figures 5D (ED) and 5E(PD), as well as Figure S3A (PD-only) and B(ED-only)
graph1=ggplot(
  data = Elevation_sex_frame
  , mapping = aes(
    y =fit^2
    , x =Elevation
    , group=Sex, colour=Sex
  )
)+
  geom_ribbon(
    mapping = aes(x=Elevation, ymin = lower^2, 
                  ymax = upper^2,fill=Sex, linetype=Sex),colour=NA
    ,alpha=.2)+
  geom_line(mapping = aes(
      colour =Sex,linetype=Sex
    )
    , size=1.1
  )+
  labs(
    x = 'Elevation'
    ,title= 'ED ommatidia by Sex over elevation'
    , y = 'ED'
  )+
  geom_point(data=eyes,
             mapping = aes( x=Elevation, y=ED,colour =Sex,shape=Sex,fill=Sex),
             size=3.5,colour="black"
  )

graph1
graph1_2=graph1+scale_colour_manual(values=c("#fdae61","#2c7bb6"))+theme_bw()+
  scale_fill_manual(values=c("#fdae61","#2c7bb6"))+scale_shape_manual(values=c(21,22))
graph1_2

#Changing the linetypes of the treatments
graph1_3=graph1_2+scale_linetype_manual(values=c("solid","solid"))
graph1_3
#Rescaling the axes
graph1_4=graph1_3+scale_x_continuous(breaks=c(-30,0,30,60))
graph1_4
#Resizing and moving the text in the graph
graph1_5=graph1_4+theme(legend.title = element_text(colour="Black", size=16, face="bold"))+
  theme(legend.text = element_text(colour="Black", size = 16, face = "bold"))+
  theme(axis.text=element_text(colour="Black", size=16, face="plain"))+
  theme(axis.title=element_text(colour="Black", size=16, face="bold"))+
  theme(plot.title = element_blank())+
  theme(axis.title.y=element_text(vjust=1))+
  theme(strip.text=element_text(colour="Black", size = 16, face = "plain"))
graph1_5
library(grid)
graph1_6=graph1_5+theme(panel.spacing = unit(0.5, "lines"),
                        legend.position="right")
graph1_6

#remove gridlines and top/right borders and facets
graph1_7=graph1_6+ theme(axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border=element_blank())+
  theme(axis.line = element_line(color = 'black'))+
  theme(strip.text=element_blank())+
  theme(strip.background=element_blank())
graph1_7


# ##Now trying with log10
# model1<-lmer(log10(ED)~Elevation*Sex+
# (Elevation|Individual),data=eyes,REML=FALSE)
# summary(model1)
# Anova(model1)
# 
# model2<-lmer(log10(ED)~Elevation+Sex+
# (Elevation|Individual),data=eyes,REML=FALSE)
# summary(model2)
# Anova(model2)
# anova(model1,model2)
# 
# ##Now dropping each individual factor
# model2.1=update(model2,~.-Elevation)
# model2.2=update(model2,~.-Sex)
# 
# #enabling parallelization
# library(parallel)
# (nc <- detectCores())
# ## Create clusters
# cl <- makeCluster(rep("localhost", 3))
# 
# #testing each model against the next largest model
# A=PBmodcomp(model1,model2,nsim=10000,cl=cl)
# A
# B_1=PBmodcomp(model2,model2.1,nsim=10000,cl=cl)
# B_1
# B_2=PBmodcomp(model2,model2.2,nsim=10000,cl=cl)
# B_2

# stopCluster(cl)
# options(scipen=999)
# mctab=data.frame(rbind(A$test,B_1$test,B_2$test))
# modcomptable=data.frame(Model=c("model2","model2","model2.1","model2.1","model2.2",
# "model2.2"),Test=c("LRT","PBtest","LRT","PBtest",
# "LRT","PBtest"),Statistic=mctab$stat,P_value=mctab$p.value,
# Factor=c("Elevation:Sex","Elevation:Sex","Elevation","Elevation","Sex","Sex"),
# Significance=c("ns","ns","**","**",".","ns"))
# modcomptable
# write.csv(modcomptable,"PBmodcomp_factor_significance_Gary_PD_ED.csv")



##Now pbmodcomp for PD and ED ONLY stats, start with the full model and drop factors hierarchically 
model1<-lmer(sqrt(Only_PD)~Elevation*Sex+
(1|Individual),data=eyes,REML=FALSE)
summary(model1)
Anova(model1)
##Dropping the interaction
model2<-lmer(sqrt(Only_PD)~Elevation+Sex+
(1|Individual),data=eyes,REML=FALSE)
summary(model2)
Anova(model2)
anova(model1,model2)

##Now dropping each individual factor
model2.1=update(model2,~.-Elevation)


model2.2=update(model2,~.-Sex)


#enabling parallelization
library(parallel)
(nc <- detectCores())
## Create clusters
cl <- makeCluster(rep("localhost", nc-1))

#testing each model against the next largest model
A=PBmodcomp(model1,model2,nsim=10000,cl=cl)
A
B_1=PBmodcomp(model2,model2.1,nsim=10000,cl=cl)
B_1
B_2=PBmodcomp(model2,model2.2,nsim=10000,cl=cl)
B_2

stopCluster(cl)
options(scipen=999)
mctab=data.frame(rbind(A$test,B_1$test,B_2$test))
modcomptable=data.frame(Model=c("model2","model2","model2.1","model2.1","model2.2",
"model2.2"),Test=c("LRT","PBtest","LRT","PBtest",
"LRT","PBtest"),Statistic=mctab$stat,P_value=mctab$p.value,
Factor=c("Elevation:Sex","Elevation:Sex","Elevation","Elevation","Sex","Sex"),
Significance=c("ns","ns","***","***","**","**"))
modcomptable
write.csv(modcomptable,"PBmodcomp_factor_significance_Gary_ONLY_PD_final.csv")


###Figure 5C
dpsi<-read.table("R_ready_delta_psi_workbook.csv",header=T,sep=",")
dpsi=subset(dpsi,Percent>0)
dpsi_0=subset(dpsi,Elevation==0)
dpsi_30=subset(dpsi,Elevation==30)
str(dpsi)
graph1=ggplot(
  data = dpsi_0
  , mapping = aes(
    y =Percent
    , x =dPsi
    ,  group=Individual,colour=Sex,fill=Sex
  )
)+
  geom_smooth(method='loess',se=FALSE,
              mapping = aes(
                fill =Sex,colour=Sex
              )
              , size=1,linetype='dashed'
  )+
  labs(
    x = 'delta-Psi'
    ,title= '0 Theta Nearest-neighbor Psi differences (delta-Psi)'
    , y = 'Percent', 
    colour="Sex"
  )+
  geom_smooth(method='loess',se=TRUE,
              mapping = aes(group=Sex,
                fill =Sex,colour=Sex
              )
              , size=1.1
  )
graph1
graph1_2=graph1+scale_colour_manual(values=c("#fdae61","#2c7bb6"))+theme_bw()+
  scale_fill_manual(values=c("#fdae61","#2c7bb6"))+scale_shape_manual(values=c(21,22))
graph1_2

#Changing the linetypes of the treatments
graph1_3=graph1_2+scale_linetype_manual(values=c("solid","solid"))
graph1_3
#Rescaling the axes
graph1_4=graph1_3+scale_x_continuous(breaks=c(0,15,30,45))
graph1_4
#Resizing and moving the text in the graph
graph1_5=graph1_4+theme(legend.title = element_text(colour="Black", size=16, face="bold"))+
  theme(legend.text = element_text(colour="Black", size = 16, face = "bold"))+
  theme(axis.text=element_text(colour="Black", size=16, face="plain"))+
  theme(axis.title=element_text(colour="Black", size=16, face="bold"))+
  theme(plot.title = element_blank())+
  theme(axis.title.y=element_text(vjust=1))+
  theme(strip.text=element_text(colour="Black", size = 16, face = "plain"))
graph1_5
library(grid)
graph1_6=graph1_5+theme(panel.spacing = unit(0.5, "lines"),
                        legend.position="right")
graph1_6

#remove gridlines and top/right borders and facets
graph1_7=graph1_6+ theme(axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border=element_blank())+
  theme(axis.line = element_line(color = 'black'))+
  theme(strip.text=element_blank())+
  theme(strip.background=element_blank())
graph1_7




