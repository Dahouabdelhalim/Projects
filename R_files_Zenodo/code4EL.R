#load packages
library(tidyverse)
library(ez)
library(ggsci)
library(gridExtra)
library(scales)
pal=pal_npg("nrc")(9) # set color palette

#repeated measures anova
re1120avers_gw <- read.csv("data4rs.csv")
re1120avers_gw$yr <- factor(re1120avers_gw$yr) 
re1120avers_gw$plot <- factor(re1120avers_gw$plot)
re1120avers_gw$trt <- factor(re1120avers_gw$trt)
re1120avers_gwrm <- ezANOVA(data =re1120avers_gw, wid = plot, dv = rs, within = yr, between = trt)
print(re1120avers_gwrm)

#figure 1
retrt1 <- re1120avers_gw %>% group_by(yr,trt) %>% summarise(r = mean(rs),sd = sd(rs),se = sd/sqrt(3)) 
retrt$trt <- ordered(retrt$trt, levels = c("control","N20","N50","N100"))
retrt$yr <- factor(retrt$yr)
figure1 <-ggplot(data =retrt,aes(x = yr, y = r,group = trt))+
        geom_bar(aes(fill = trt),position = position_dodge(),stat = "identity")+
        geom_errorbar(aes(col = trt, ymin = r-se,ymax = r+se),position = position_dodge(0.9),width =0.3)+
        scale_fill_manual(values = pal)+
        scale_color_manual(values = pal)+
        theme_bw()+
        theme(panel.grid = element_blank(),strip.background = element_blank())+
        theme(legend.position = "top")+
        theme( strip.background = element_blank(),strip.text.x = element_text( size=12))+
        theme(axis.text.x = element_text(size = 12)) +
        theme(axis.text.y = element_text(size=12)) +
        theme(axis.title.x = element_text(size=12)) +
        theme(axis.title.y = element_text( size=12),strip.text = element_text(size = 12))+
        xlab("Year")+
        ylab("Rs")

#figure 2

re1120es <- read.csv("data4es.csv")#calculation of the effect size are descried  in the paper
pal_sub <- c("#4DBBD5FF","#00A087FF","#3C5488FF")#set palette
re1120es$trt <- ordered(re1120es$trt,levels = c("N20","N50","N100"))
es <-ggplot(data = re1120es,aes(x = dur, y = es))+
        geom_point(size =4,aes(col = trt))+
        geom_smooth(method = "lm",fill = "lightgrey",col = "#F39B7FFF")+
        facet_grid(.~trt)+
        scale_x_continuous(breaks = pretty_breaks())+
        stat_smooth_func(method = "lm",geom = "text",parse =T,size =0.4)+
        geom_hline(yintercept = 0)+
        ylab("Effect size of Rs")+
        xlab("Experimental duration (years)")+
        theme_bw()+
        theme(panel.grid = element_blank(),strip.background = element_blank(),panel.spacing = unit(0.9,"cm"))+
        scale_color_manual(values = pal_sub)+
        theme(legend.position = "none")+
        theme( strip.background = element_blank(),strip.text.x = element_text( size=14))+
        theme(axis.text.x = element_text(size = 16)) +
        theme(axis.text.y = element_text(size=16)) +
        theme(axis.title.x = element_text(size=12)) +
        theme(axis.title.y = element_text( size=16),strip.text = element_text(size = 14))

#effect size during different stages
#first 3 years
re1120es3<- subset(re1120es, dur %in% c(1:3))
re1120es3trt <- re1120es3 %>% group_by(trt) %>% summarise(m = mean(es),sd = sd(es),se = sd/sqrt(3))
#years 5-8
re1120es8 <- subset(re1120es, dur %in% c(5:8))
re1120aes8trt <- re1120es8 %>% group_by(trt) %>% summarise(m = mean(es),sd = sd(es),se = sd/sqrt(3))
# years 9-11
re1120es11 <- subset(re1120es, dur %in% c(9:11))
re1120es11trt <- re1120es11 %>% group_by(trt) %>% summarise(m = mean(es),sd = sd(es),se = sd/sqrt(3))
#effect size by durations
esdur <- rbind(re1120es3trt,re1120es8trt,re1120es11trt)
dur <- rep(c("Early stage","Mid-stage","Late stage"),each =3)

stage <-ggplot(data = esdur,aes(x = dur, y =m,fill = trt))+
        geom_bar(stat = "identity",position = position_dodge())+
        scale_fill_manual(values = pal_sub)+
        geom_errorbar(aes(ymin = m-se,ymax = m+se),position = position_dodge(0.9),width =0.2)+
        theme_bw()+
        theme(panel.grid = element_blank(),strip.background = element_blank())+
        scale_color_manual(values = pal_sub)+
        theme(legend.position = "none")+
        ylab("Effect size of Rs")+
        theme( strip.background = element_blank(),strip.text.x = element_text( size=14))+
        theme(axis.text.x = element_text(size = 16)) +
        theme(axis.text.y = element_text(size=16)) +
        theme(axis.title.x = element_text(size=12)) +
        theme(axis.title.y = element_text( size=16),strip.text = element_text(size = 14))+
        scale_y_continuous(labels = scales::number_format(accuracy = 0.1,
                                                          decimal.mark = '.')) 

figure2 <- grid.arrange(arrangeGrob(es+ 
                                           theme(strip.text.x = element_text(size = 5)) + theme(legend.position="none") +
                                           theme( plot.margin=unit(c(0.3,1,0,0), "cm"))
                                   ,                                     stage +
                                           theme(strip.text.x = element_blank(),
                                                 strip.background = element_blank()) +
                                           theme( plot.margin=unit(c(0.3,1,0,0), "cm")) +
                                           theme(legend.position="none") ,nrow = 2,heights = c(6,6)))

#figure 3
rsmic <- read.csv("data4Fun.csv")
rsmic$trt <- ordered(rsmic$trt,levels = c("control","N20","N50","N100"))
figure3 <- ggplot(data = rsmic,aes(x = oturich,y = rs))+
        geom_point(aes(color =trt),size =3)+
        facet_wrap(yr~guild,scales ="free_x")+
        geom_smooth(method = "lm",color= "black")+
        #stat_smooth_func(method = "lm",geom = "text",parse = T,size =0.9)+
        theme_bw()+
        theme(panel.grid = element_blank(),strip.background = element_blank())+
        theme(legend.position = "bottom",legend.text = element_text(size = 16))+
        theme( strip.background = element_blank(),strip.text.x = element_text( size=14))+
        theme(axis.text.x = element_text(size = 16)) +
        theme(axis.text.y = element_text(size=16)) +
        theme(axis.title.x = element_text(size=16)) +
        theme(axis.title.y = element_text( size=16),strip.text = element_text(size = 14))+
        scale_x_continuous(breaks = pretty_breaks())+
        scale_color_npg()+
        theme(panel.spacing = unit(0.9,"cm"))


#figure 4
fineroot <- read.csv("data4Froot.csv")
head(fineroot)
fineroot$trt <- ordered(fineroot$trt, levels = c("control","N20","N50","N100"))

figure4 <- ggplot(data = fineroot,aes(x = fineroot, y = rs))+
        geom_point(size =4,aes(col = factor(n)))+
        geom_smooth(method = "lm",fill ="lightgrey",col ="#F39B7FFF")+
        #stat_smooth_func(method = "lm",geom = "text",parse =T,size =3)+
        scale_color_manual(values = pal)+
        ylab("Soil respiration")+
        xlab("Fine root biomass")+
        theme_bw()+
        theme(panel.grid = element_blank(),strip.background = element_blank())+
        theme(legend.position = "top")+
        theme( strip.background = element_blank(),strip.text.x = element_text( size=14))+
        theme(axis.text.x = element_text(size = 16)) +
        theme(axis.text.y = element_text(size=16)) +
        theme(axis.title.x = element_text(size=14)) +
        theme(axis.title.y = element_text( size=16),strip.text = element_text(size = 14))
