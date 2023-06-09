library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(ggpubr)

########################################################
#### export from GIS vegetation type classification ####
data <- read.table("Vegetation classification GIS.csv", sep=";", header=T)

## summary table ##
rowSE <- function(x) {sd(x)/sqrt(length(x))} 

summ <- data %>%
  group_by(classes_new) %>%
  summarise(N = n(),
            m.Elev = mean(SRTM_WGS84),
            Q1 = quantile(SRTM_WGS84, probs = 0.25),
            Q2 = quantile(SRTM_WGS84, probs = 0.5),
            Q3 = quantile(SRTM_WGS84, probs = 0.75),
            min = min(SRTM_WGS84),
            max = max(SRTM_WGS84),
            SD = sd(SRTM_WGS84),
            SE = rowSE(SRTM_WGS84))





###############
#### ANOVA ####


## check for normal distribution
ggdensity(data$SRTM_WGS84)
ggqqplot(data$SRTM_WGS84)
## check


aov1<-aov(SRTM_WGS84 ~ classes_new, data=data)
summary(aov1)
#                 Df    Sum Sq  Mean Sq F value Pr(>F)    
#  classes_new     3 142642928 47547643   13347 <2e-16 ***
#  Residuals   37345 133042070     3563                   
#---


## posthoc Tukey test for pairwise differences between the classes
TukeyHSD(aov1, data=data)
#                                           diff        lwr       upr p adj
#Grassland         - Ecotone grassland  12.69899   9.275946  16.12204     0
#Miombo  forest    - Ecotone grassland 139.39950 137.016554 141.78245     0
#Plateau grassland - Ecotone grassland 201.54754 198.180471 204.91461     0
#Miombo  forest    - Grassland         126.70051 123.892278 129.50874     0
#Plateau grassland - Grassland         188.84855 185.168211 192.52889     0
#Plateau grassland - Miombo  forest     62.14804  59.408316  64.88776     0
### all classes differ more or less distinct

####################################################################################


data <- data %>%
  unique()
data$classes_new <- factor(data$classes_new, levels=c("Plateau grassland", "Miombo forest", "Ecotone grassland", "Grassland"))

####scatterplot with smoothing lines
png("20210917 Treelines east-west gradient color2.png",res=600, width=5000, height=2500)
ggplot(data, aes(x=Lon, y=SRTM_WGS84, col=classes_new))+
  facet_wrap(~classes_new, nrow = 1, ncol = 4)+
  theme_classic()+
  theme(aspect.ratio = .8, panel.border = element_rect(colour = "black", fill=NA), legend.position = "none", 
        strip.text = element_text(size=16) ,strip.background = element_rect(colour = "white"), 
        axis.text = element_text(size=11), axis.title = element_text(size=15))+
  scale_x_continuous(name="Longitude [°]", breaks=seq(18,20,0.5), limits=c(min(data$Lon)-.05,max(data$Lon)+.05))+
  scale_y_continuous(name="Elevation [m]", breaks=seq(1200,1650,100), limits=c(1150,1680))+
  
  scale_color_manual(values=c("#42929e", "#10652f","#5bcc13",  "#ffff00"))+
  geom_hline(yintercept=seq(1200,1600,100), col="grey")+
  geom_point(size=1, show.legend = F, alpha=.7)+
  geom_smooth(data=data[data$classes_new=="Ecotone grassland",], formula=y~s(x, k=4, bs="cs"), method = "gam",col="black",  size=1.5)+
  geom_smooth(data=data[data$classes_new=="Grassland",], formula=y~s(x, k=4, bs="cs"), method = "gam", col="black", size=1.5)+
  geom_smooth(data=data[data$classes_new=="Miombo forest",], formula=y~s(x, k=4, bs="cs"), method = "gam", col="black", size=1.5)+
  geom_smooth(data=data[data$classes_new=="Plateau grassland",], formula=y~s(x, k=4, bs="cs"), method = "gam", col="black", size=1.5)
dev.off()
