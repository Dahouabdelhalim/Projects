library(nlme)
library(lme4)
library(car)
library(effects)
library(lsmeans)
library(pbkrtest)
library(effects)
#install.packages("ggplot2")
library("r2glmm")
library(ggplot2)
library(MuMIn)
library(semEff)
library(devtools)

# setwd("C:/Users/richa/Dropbox/Harvard/Research/Gary polarization manuscript work/NY_CCT_specimen_model_polDATA")
# setwd("C:/Users/Astaroph/Dropbox/Harvard/Research/Gary polarization manuscript work/NY_CCT_specimen_model_polDATA")
dop=read.csv('7-27-20_new_spectra_smooth_R_ready_compiled.csv')
robots=read.csv('Robot_Refs_combined_R_ready.csv')
comb=read.csv('7-27-20_new_spectra_smooth_R_ready_compiled_with_robots.csv')
comb_standard=subset(comb,Illumination=='Standard')
dop_standard=subset(dop,Illumination=='Standard')

comb_standard1=subset(comb_standard,Side!='Left')
library(ggplot2)

violin_1=ggplot(
  data = dop_standard
  , mapping = aes(
    y =Reflectance,
    , x=wavelength
    , colour=factor(Sex)
  )
)+
  labs(
    x = 'Wavelength (nm)'
    ,title= 'Reflectance by Wavelength'
    , y = 'Reflectance'
  )+
  # geom_point(mapping = aes(x=wavelength,
  #                           colour =factor(Sex)
  #             )
  #             , size=1.8,alpha=0.05
  # )+
  geom_smooth(method="loess",alpha=.4, se =TRUE)+
  facet_grid(Wing~Illumination)
violin_1
violin_1_2=violin_1+scale_color_manual(values=c("#2c7bb6","#fdae61"))+
  scale_fill_manual(values=c("#2c7bb6","#fdae61"))+theme_bw()
violin_1_2

violin_1_3=violin_1_2

#violin_1_3=violin_1_2+scale_y_log10()
violin_1_3
violin_1_4=violin_1_3
violin_1_4
#Resizing and moving the text in the graph
violin_1_5=violin_1_4+theme(legend.title = element_text(colour="Black", size=12, face="bold"))+
  theme(legend.text = element_text(colour="Black", size = 12, face = "bold"))+
  theme(axis.text=element_text(colour="Black", size=10, face="plain"))+
  theme(axis.title=element_text(colour="Black", size=12, face="bold"))+
  theme(axis.title.y=element_text(vjust=1))+
  theme(strip.text=element_text(colour="Black", size = 10, face = "plain"))
violin_1_5
library(grid)
violin_1_6=violin_1_5+theme(panel.spacing = unit(0.5, "lines"),
                            legend.position="right")
violin_1_6

#remove gridlines and top/right borders and facets
violin_1_7=violin_1_6+ theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border=element_blank())+
  theme(axis.line = element_line(color = 'black'))
#theme(strip.text=element_blank())+
#theme(strip.background=element_blank())
violin_1_7



##this is code to create a dataframe that ggplot inherits in the stat_summary
##function, to make simple standard error bars
mean_se <- function(x, mult = 1) {  
  x <- na.omit(x)
  se <- mult * sqrt(var(x) / length(x))
  mean <- mean(x)
  data.frame(y = mean, ymin = mean - se, ymax = mean + se)
}

library(ggplot2)

##Figure 4B
violin_1=ggplot(
  data = dop_standard
  , mapping = aes(
    y =Reflectance,
    , x=wavelength
    , colour=Sex
  )
)+
  labs(
    x = 'Wavelength (nm)'
    ,title= 'Reflectance by Wavelength'
    , y = 'Reflectance'
  )+
  # geom_point(mapping = aes(x=wavelength,
  #                           colour =Sex
  #             )
  #             , size=1.8,alpha=0.05
  # )+
  # geom_smooth(method="gam",alpha=.4, se =FALSE,mapping=aes(linetype=Wing))+
  stat_summary(fun.data=mean_se, geom="ribbon", alpha=0.15,mapping=aes(linetype=Wing,fill=Sex))+
  stat_summary(fun.data=mean_se, geom="line", alpha=0.95,size=1,mapping=aes(linetype=Wing,color=Sex))+
  facet_grid(.~.)
violin_1
violin_1_2=violin_1+scale_color_manual(values=c("#fdae61","#2c7bb6"))+
  scale_fill_manual(values=c("#fdae61","#2c7bb6"))+theme_bw()
violin_1_2

violin_1_3=violin_1_2

#violin_1_3=violin_1_2+scale_y_log10()
violin_1_3
violin_1_4=violin_1_3
violin_1_4
#Resizing and moving the text in the graph
violin_1_5=violin_1_4+theme(legend.title = element_text(colour="Black", size=12, face="bold"))+
  theme(legend.text = element_text(colour="Black", size = 14, face = "bold"))+
  theme(axis.text=element_text(colour="Black", size=12, face="plain"))+
  theme(axis.title=element_text(colour="Black", size=14, face="bold"))+
  theme(axis.title.y=element_text(vjust=1))+
  theme(strip.text=element_text(colour="Black", size = 12, face = "plain"))
violin_1_5
library(grid)
violin_1_6=violin_1_5+theme(panel.spacing = unit(0.5, "lines"),
                            legend.position="right")
violin_1_6

#remove gridlines and top/right borders and facets
violin_1_7=violin_1_6+ theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border=element_blank())+
  theme(axis.line = element_line(color = 'black'))
#theme(strip.text=element_blank())+
#theme(strip.background=element_blank())
violin_1_7


###############################################
##now code to combine the robot and specimen reflectance plots

library(ggplot2)

comb_standard$Sex <- factor(comb_standard$Sex,levels = c("F", "M", "HNPB depolarized", "HNPB polarized",
                                                         "3M_Blue_R374_depol","3M_Blue_R374_pol",
                                                         "3M_Yellow_R312_depol","3M_Yellow_R312_pol",
                                                         "3M_Red_G280_depol","3M_Red_G280_pol"))
##Figure 3A
violin_1=ggplot(
  data = comb_standard
  , mapping = aes(
    y =Reflectance,
    , x=wavelength
    , colour=Sex
  )
)+
  labs(
    x = 'Wavelength (nm)'
    ,title= 'Reflectance by Wavelength'
    , y = 'Reflectance'
  )+
  # geom_point(mapping = aes(x=wavelength,
  #                           colour =Sex
  #             )
  #             , size=1.8,alpha=0.05
  # )+
  # geom_smooth(method="gam",alpha=.4, se =FALSE,mapping=aes(linetype=Wing))+
  stat_summary(fun.data=mean_se, geom="ribbon", alpha=0.15,mapping=aes(linetype=Wing,fill=Sex))+
  stat_summary(fun.data=mean_se, geom="line", alpha=0.95,size=1,mapping=aes(linetype=Wing,color=Sex))+
  facet_grid(.~.)
violin_1
violin_1_2=violin_1+scale_color_manual(values=c("#fdae61","#2c7bb6","#9E71A8","#632770","#669999","#226666","#D4CB6A","#AAA039","#CA6573","#A23645"))+
  scale_fill_manual(values=c("#fdae61","#2c7bb6","#9E71A8","#632770","#669999","#226666","#D4CB6A","#AAA039","#CA6573","#A23645"))+theme_bw()
violin_1_2

violin_1_3=violin_1_2+scale_linetype_manual(values=c("solid","dashed","dotdash"))

#violin_1_3=violin_1_2+scale_y_log10()
violin_1_3
violin_1_4=violin_1_3
violin_1_4
#Resizing and moving the text in the graph
violin_1_5=violin_1_4+theme(legend.title = element_text(colour="Black", size=12, face="bold"))+
  theme(legend.text = element_text(colour="Black", size = 14, face = "bold"))+
  theme(axis.text=element_text(colour="Black", size=12, face="plain"))+
  theme(axis.title=element_text(colour="Black", size=14, face="bold"))+
  theme(axis.title.y=element_text(vjust=1))+
  theme(strip.text=element_text(colour="Black", size = 12, face = "plain"))
violin_1_5
library(grid)
violin_1_6=violin_1_5+theme(panel.spacing = unit(0.5, "lines"),
                            legend.position="right")
violin_1_6

#remove gridlines and top/right borders and facets
violin_1_7=violin_1_6+ theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border=element_blank())+
  theme(axis.line = element_line(color = 'black'))
#theme(strip.text=element_blank())+
#theme(strip.background=element_blank())
violin_1_7





dop_FW=subset(dop,Wing=="Forewing")
dop_solar=subset(dop_FW,Illumination=="Solar")


###The following code conducts a cross-correlation analysis between robot and specimen 
##reflectance traces to quantify the spectral similarity between the male and female specimens and 
##the robots
###Running cross correlation analysis on robot models and specimens
str(robots)

robot_solar=subset(robots,Illumination=="Solar")
robot_standard=subset(robots,Illumination=='Standard')
Udep=subset(robot_standard,Wing_type=="UV depolarized")
Upol=subset(robot_standard,Wing_type=="UV polarized")
Bdep=subset(robot_standard,Wing_type=="Blue depolarized")
Bpol=subset(robot_standard,Wing_type=="Blue polarized")
Ydep=subset(robot_standard,Wing_type=="Yellow depolarized")

Ypol=subset(robot_standard,Wing_type=="Yellow polarized")
Rdep=subset(robot_standard,Wing_type=="Red depolarized")
Rpol=subset(robot_standard,Wing_type=="Red polarized")

dopFW=subset(dop,Wing=="FW")
dopHW=subset(dop,Wing=="HW")

dop_solar=subset(dopFW,Illumination=="Solar")
dop_standard=subset(dop,Illumination=="Standard")
F1=subset(dop_standard,Individual=="F1")
xcor=ccf(F1$Reflectance,Rdep$Reflectance)
xcor$acf[24]

cross_cor=function(Ind_str,R_string,ind_frame,R_frame){
  indyframe=subset(ind_frame,Individual==Ind_str)
  robotframe=subset(R_frame,Wing_type==R_string)
  xcor=ccf(indyframe$Reflectance,robotframe$Reflectance)
  cor=xcor$acf[24]
  cor
}
cor1=cross_cor("F1","UV depolarized",dop_standard,robot_standard)
cor1
cross_cor_frame=function(color){
  cor1=cross_cor("F1",color,dop_standard,robot_standard)
  cor2=cross_cor("F2",color,dop_standard,robot_standard)
  cor3=cross_cor("F3",color,dop_standard,robot_standard)
  cor4=cross_cor("F4",color,dop_standard,robot_standard)
  cor5=cross_cor("F5",color,dop_standard,robot_standard)
  cor6=cross_cor("F6",color,dop_standard,robot_standard)
  cor7=cross_cor("F7",color,dop_standard,robot_standard)
  cor8=cross_cor("F8",color,dop_standard,robot_standard)
  cor9=cross_cor("F9",color,dop_standard,robot_standard)
  cor10=cross_cor("F10",color,dop_standard,robot_standard)
  cor11=cross_cor("F11",color,dop_standard,robot_standard)
  cor12=cross_cor("F12",color,dop_standard,robot_standard)
  cor13=cross_cor("M1",color,dop_standard,robot_standard)
  cor14=cross_cor("M2",color,dop_standard,robot_standard)
  cor15=cross_cor("M3",color,dop_standard,robot_standard)
  cor16=cross_cor("M4",color,dop_standard,robot_standard)
  cor17=cross_cor("M5",color,dop_standard,robot_standard)
  cor18=cross_cor("M6",color,dop_standard,robot_standard)
  cor19=cross_cor("M7",color,dop_standard,robot_standard)
  cor20=cross_cor("M8",color,dop_standard,robot_standard)
  cor21=cross_cor("M9",color,dop_standard,robot_standard)
  cor22=cross_cor("M10",color,dop_standard,robot_standard)
  cor23=cross_cor("M11",color,dop_standard,robot_standard)
  cor24=cross_cor("M12",color,dop_standard,robot_standard)
  cor25=cross_cor("M13",color,dop_standard,robot_standard)
  cors=c(cor1,cor2,cor3,cor4,cor5,cor6,cor7,cor8,cor9,cor10,cor11,cor12,cor13,cor14,cor15,cor16,cor17,cor18,cor19,
         cor20,cor21,cor22,cor23,cor24,cor25)
  colors=c(color,color,color,color,color,color,color,color,color,color,color,color,color,color,color,color,
           color,color,color,color,color,color,color,color,color)
  individuals=c('F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12'
                ,'M1','M2','M3','M4','M5','M6','M7','M8','M9','M10','M11','M12','M13')
  sex=c("F","F","F","F","F","F","F","F","F","F","F","F",
        "M","M","M","M","M","M","M","M","M","M","M","M","M")
  colframe=data.frame(individuals,sex,colors,cors)
  colframe
}
UV_dep=cross_cor_frame("UV depolarized")
UV_pol=cross_cor_frame("UV polarized")
Blue_dep=cross_cor_frame("Blue depolarized")
Blue_pol=cross_cor_frame("Blue polarized")
Yellow_dep=cross_cor_frame("Yellow depolarized")
Yellow_pol=cross_cor_frame("Yellow polarized")
Red_dep=cross_cor_frame("Red depolarized")
Red_pol=cross_cor_frame("Red polarized")


final_frame=rbind(UV_dep,UV_pol,Blue_dep,Blue_pol,Yellow_dep,Yellow_pol,Red_dep,Red_pol)

#write.csv(final_frame,"Individual_cross_correlation_analysis_with_robots_10_1_20.csv")



###Figure 3b
violin_1=ggplot(
  data = final_frame
  , mapping = aes(
    y =cors
    , x=sex
    , colour=factor(colors)
  )
)+
  labs(
    x = 'Specimen sex'
    ,title= 'Cross-correlation between wing model and specimen reflectance'
    , y = 'Cross-correlation'
  )+
  geom_point(position=position_jitterdodge(dodge.width = .9),
             mapping = aes(x=sex,
                           colour =factor(colors)
             )              , size=1.8,alpha=0.25)+
  geom_violin(alpha=.4, draw_quantiles = c(0.25, 0.5, 0.75),
              mapping=aes(x=sex,fill=factor(colors),colour=factor(colors)))+
  #  geom_smooth(method="loess",alpha=.4, se =TRUE)+
  facet_grid(.~.)
violin_1
violin_1_2=violin_1+scale_color_manual(values=c("#9E71A8","#632770","#669999","#226666","#D4CB6A","#AAA039","#CA6573","#A23645"))+
  scale_fill_manual(values=c("#9E71A8","#632770","#669999","#226666","#D4CB6A","#AAA039","#CA6573","#A23645"))+theme_bw()
violin_1_2

violin_1_3=violin_1_2

#violin_1_3=violin_1_2+scale_y_log10()
violin_1_3
violin_1_4=violin_1_3
violin_1_4
#Resizing and moving the text in the graph
violin_1_5=violin_1_4+theme(legend.title = element_text(colour="Black", size=12, face="bold"))+
  theme(legend.text = element_text(colour="Black", size = 12, face = "bold"))+
  theme(axis.text=element_text(colour="Black", size=10, face="plain"))+
  theme(axis.title=element_text(colour="Black", size=12, face="bold"))+
  theme(axis.title.y=element_text(vjust=1))+
  theme(strip.text=element_text(colour="Black", size = 10, face = "plain"))
violin_1_5
library(grid)
violin_1_6=violin_1_5+theme(panel.spacing = unit(0.5, "lines"),
                            legend.position="right")
violin_1_6

#remove gridlines and top/right borders and facets
violin_1_7=violin_1_6+ theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border=element_blank())+
  theme(axis.line = element_line(color = 'black'))
#theme(strip.text=element_blank())+
#theme(strip.background=element_blank())
violin_1_7


####FWHM calculations
dop_standard=subset(dop,Illumination=='Standard')


mean_se <- function(x, mult = 1) {  
  x <- na.omit(x)
  se <- mult * sqrt(var(x) / length(x))
  mean <- mean(x)
  data.frame(y = mean, ymin = mean - se, ymax = mean + se)
}


violin_1=ggplot(
  data = dop_standard
  , mapping = aes(
    y =Reflectance,
    , x=wavelength
    , colour=Sex
  )
)+
  labs(
    x = 'Wavelength (nm)'
    ,title= 'Reflectance by Wavelength'
    , y = 'Reflectance'
  )+
  # geom_point(mapping = aes(x=wavelength,
  #                           colour =Sex
  #             )
  #             , size=1.8,alpha=0.05
  # )+
  # geom_smooth(method="gam",alpha=.4, se =FALSE,mapping=aes(linetype=Wing))+
  stat_summary(fun.data=mean_se, geom="ribbon", alpha=0.15,mapping=aes(linetype=Wing,fill=Sex))+
  stat_summary(fun.data=mean_se, geom="line", alpha=0.95,size=1,mapping=aes(linetype=Wing,color=Sex))+
  facet_grid(.~.)
violin_1



plot_info<-ggplot_build(violin_1)
summarized_spectra<-plot_info$data[[2]]
summarized_spectra<-subset(summarized_spectra,x<720)
summarized_spectra$Sex<-summarized_spectra$colour
summarized_spectra$Sex[summarized_spectra$Sex=='#F8766D']<-'F'
summarized_spectra$Sex[summarized_spectra$Sex=='#00BFC4']<-'M'
summarized_spectra$Wing<-summarized_spectra$linetype
summarized_spectra$Wing[summarized_spectra$Wing=='solid']<-'Forewing'
summarized_spectra$Wing[summarized_spectra$Wing==22]<-'Hindwing'
sum_FFW<-subset(summarized_spectra,Wing=='Forewing'& Sex=='F')
sum_FHW<-subset(summarized_spectra,Wing=='Hindwing'& Sex=='F')
sum_MFW<-subset(summarized_spectra,Wing=='Forewing'& Sex=='M')
sum_MHW<-subset(summarized_spectra,Wing=='Hindwing'& Sex=='M')

plot(sum_FHW$x,sum_FHW$y)
x=sum_MHW$y


#
data<-sum_FHW
d<-data.frame(x=data$x,y=data$y)
# d <- density(na.omit(x),n=1e4)
xmax <- d$x[d$y==max(d$y)]
x1 <- d$x[d$x < xmax][which.min(abs(d$y[d$x < xmax]-max(d$y)/2))]
x2 <- d$x[d$x > xmax][which.min(abs(d$y[d$x > xmax]-max(d$y)/2))]
plot(c(x1, x2), c(d$y[d$x==x1], d$y[d$x==x2]), col="red")
(FWHM <- x2-x1)
c(xmax,FWHM)
