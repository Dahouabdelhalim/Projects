# clear workspace
rm(list = ls(all.names = TRUE)) 

library(ggplot2)
library(ggpubr)
library(cowplot)

# library(showtext)
# showtext_auto()
# font_add(family = "Myriad Pro", regular = "./fonts/Myriad Pro Regular.ttf")
# showtext_auto()

library(extrafont)
font_import(path = "./fonts") 

DataAll <- read.csv("Fischer_et_al_Data_masks_group_average_time_plot.txt")

dataAll.frame <- data.frame(DataAll)


##
customblue <- rgb(20/255, 111/255, 247/255) # "#146ff7" 
customred <- rgb(205/255, 17/255, 23/255)# "#cd1117"

xlimits2use <- c(-7,10) #

xlimitleft2use  <- -640/60 #


### cut the very first data points partially dominated by single subjects
dataUsed <- subset(dataAll.frame,time_sec>-640)

##
###make plots for all variables

sto2plot <- ggplot(data=dataUsed,aes(x=time_sec/60,y=CBF_mean,col=masktype) ) +
  annotate("rect", xmin = -54.45/60, xmax = 0, ymin = -Inf, ymax = Inf, 
           alpha = .5,fill = "grey") +
  geom_vline(xintercept = 3, linetype = "dashed", size = 1.2, col = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2, col = "grey") +
  geom_line(data=dataUsed,aes(x=time_sec/60,y=(StO2_mean),group=masktype,col=masktype),size=1.5)  + 
  geom_ribbon(data=dataUsed,aes(x=time_sec/60,ymin=StO2_mean-(StO2_std/1),ymax=StO2_mean+(StO2_std/1),fill=masktype),alpha=0.3,colour = NA)+
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10),
                     limits = xlimits2use) +
  scale_y_continuous(position = "left",
                     limits = c(-1.6, 3.7))+
  labs(y = expression(atop(paste(Delta*"StO"[2]),"[%]")), x = "time [min]") +
  scale_fill_manual(values=c(customred,customblue)) +
  scale_color_manual(values=c(customred,customblue)) +
  annotate("rect", xmin = 630.49/60, xmax = 600/60, ymin = -Inf, ymax = Inf, 
           alpha = .7, fill = "white") +
  theme_light()   +
  theme(text = element_text(size = 35, family="Myriad pro"), 
        axis.title.y=element_text(size=rel(0.9)), legend.position = "none",
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black")) 

cbfplot <- ggplot(data=dataUsed,aes(x=time_sec/60,y=StO2_mean,col=masktype) ) +
  annotate("rect", xmin = -54.45/60, xmax = 0, ymin = -Inf, ymax = Inf, 
           alpha = .5,fill = "grey") +
  #Xgeom_vline(xintercept = 3, linetype = "dashed", size = 1.2, col = "magenta") +
  geom_segment(aes(x=3,xend=3,y=-Inf,yend=30),linetype = "dashed", size = 1.2, col = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2, col = "grey") +
  geom_line(data=dataUsed,aes(x=time_sec/60,y=(CBF_mean),group=masktype,col=masktype),size=1.5)  + 
  geom_ribbon(data=dataUsed,aes(x=time_sec/60,ymin=CBF_mean-(CBF_std/1),ymax=CBF_mean+(CBF_std/1),fill=masktype),alpha=0.3,colour = NA)+
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10),
                     limits = xlimits2use) +
  scale_y_continuous(position = "left",
                     limits = c(-12, 35))+
  labs(y = expression(atop(paste(Delta*"rCBF"),"[%]")), x = "time [min]") +
  scale_fill_manual(values=c(customred,customblue)) +
  scale_color_manual(values=c(customred,customblue)) +
  annotate("rect", xmin = 630.49/60, xmax = 600/60, ymin = -Inf, ymax = Inf, 
           alpha = .7, fill = "white") +
  theme_light()   +
  theme(text = element_text(size = 35, family="Myriad pro"), 
        axis.title.y=element_text(size=rel(0.9)), legend.position = "none",
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black")) 
  
thcplot <- ggplot(data=dataUsed,aes(x=time_sec/60,y=THC_mean,col=masktype) ) +
  annotate("rect", xmin = -54.45/60, xmax = 0, ymin = -Inf, ymax = Inf, 
           alpha = .5,fill = "grey") +
  geom_vline(xintercept = 3, linetype = "dashed", size = 1.2, col = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2, col = "grey") +
  geom_line(data=dataUsed,aes(x=time_sec/60,y=(THC_mean),group=masktype,col=masktype),size=1.5)  + 
  geom_ribbon(data=dataUsed,aes(x=time_sec/60,ymin=THC_mean-(THC_std/1),ymax=THC_mean+(THC_std/1),fill=masktype),alpha=0.3,colour = NA)+
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10),
                     limits = xlimits2use) +
  scale_y_continuous(position = "left",
                     limits = c(-1.5, 3.2))+
  labs(y = expression(atop(paste(Delta*"tHb"),"["*mu*"M]")), x = "time [min]") +
  scale_fill_manual(values=c(customred,customblue)) +
  scale_color_manual(values=c(customred,customblue)) +
  annotate("rect", xmin = 630.49/60, xmax = 600/60, ymin = -Inf, ymax = Inf, 
           alpha = .7, fill = "white") +
  theme_light()   +
  theme(text = element_text(size = 35, family="Myriad pro"), 
        axis.title.y=element_text(size=rel(0.9)), legend.position = "none",
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black")) 

oefplot <- ggplot(data=dataUsed,aes(x=time_sec/60,y=OEF_mean,col=masktype) ) +
  annotate("rect", xmin = -54.45/60, xmax = 0, ymin = -Inf, ymax = Inf, 
           alpha = .5,fill = "grey") +
  geom_vline(xintercept = 3, linetype = "dashed", size = 1.2, col = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2, col = "grey") +
  geom_line(data=dataUsed,aes(x=time_sec/60,y=(OEF_mean),group=masktype,col=masktype),size=1.5)  + 
  geom_ribbon(data=dataUsed,aes(x=time_sec/60,ymin=OEF_mean-(OEF_std/1),ymax=OEF_mean+(OEF_std/1),fill=masktype),alpha=0.3,colour = NA)+
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10),
                     limits = xlimits2use) +
  scale_y_continuous(position = "left",
                     limits = c(-9, 4.5))+
  labs(y = expression(atop(paste(Delta*"rOEF"),"[%]")), x = "time [min]") +
  scale_fill_manual(values=c(customred,customblue)) +
  scale_color_manual(values=c(customred,customblue)) +
  annotate("rect", xmin = 630.49/60, xmax = 600/60, ymin = -Inf, ymax = Inf, 
           alpha = .7, fill = "white") +
  theme_light()   +
  theme(text = element_text(size = 35, family="Myriad pro"), 
        axis.title.y=element_text(size=rel(0.9)), legend.position = "none",
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black")) 

cmro2plot <- ggplot(data=dataUsed,aes(x=time_sec/60,y=CMRO2_mean,col=masktype) ) +
  annotate("rect", xmin = -54.45/60, xmax = 0, ymin = -Inf, ymax = Inf, 
           alpha = .5,fill = "grey") +
  geom_vline(xintercept = 3, linetype = "dashed", size = 1.2, col = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2, col = "grey") +
  geom_line(data=dataUsed,aes(x=time_sec/60,y=(CMRO2_mean),group=masktype,col=masktype),size=1.5)  + 
  geom_ribbon(data=dataUsed,aes(x=time_sec/60,ymin=CMRO2_mean-(CMRO2_std/1),ymax=CMRO2_mean+(CMRO2_std/1),fill=masktype),alpha=0.3,colour = NA)+
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10),
                     limits = xlimits2use) +
  scale_y_continuous(position = "left",
                     limits = c(-13, 40))+
  labs(y = expression(atop(paste(Delta*"rCMRO"[2]),"[%]")), x = "time [min]") +
  scale_fill_manual(values=c(customred,customblue)) +
  scale_color_manual(values=c(customred,customblue)) +
  annotate("rect", xmin = 630.49/60, xmax = 600/60, ymin = -Inf, ymax = Inf, 
           alpha = .7, fill = "white") +
  theme_light()   +
  theme(text = element_text(size = 35, family="Myriad pro"), 
        axis.title.y=element_text(size=rel(0.9)), legend.position = "none",
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black")) 


hrplot <- ggplot(data=dataUsed,aes(x=time_sec/60,y=HR_mean,col=masktype) ) +
  annotate("rect", xmin = -54.45/60, xmax = 0, ymin = -Inf, ymax = Inf, 
           alpha = .5,fill = "grey") +
  geom_vline(xintercept = 3, linetype = "dashed", size = 1.2, col = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2, col = "grey") +
  geom_line(data=dataUsed,aes(x=time_sec/60,y=(HR_mean),group=masktype,col=masktype),size=1.5)  + 
  geom_ribbon(data=dataUsed,aes(x=time_sec/60,ymin=HR_mean-(HR_std/1),ymax=HR_mean+(HR_std/1),fill=masktype),alpha=0.3,colour = NA)+
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10),
                     limits = xlimits2use) +
  scale_y_continuous(position = "left",
                     limits = c(-4, 13))+
  labs(y = expression(atop(paste(Delta*"HR"),"[1/min]")), x = "time [min]") +
  scale_fill_manual(values=c(customred,customblue)) +
  scale_color_manual(values=c(customred,customblue)) +
  annotate("rect", xmin = 630.49/60, xmax = 600/60, ymin = -Inf, ymax = Inf, 
           alpha = .7, fill = "white") +
  theme_light()   +
  theme(text = element_text(size = 35, family="Myriad pro"), 
        axis.title.y=element_text(size=rel(0.9)), legend.position = "none",
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))  

rrplot <- ggplot(data=dataUsed,aes(x=time_sec/60,y=RR_mean,col=masktype) ) +
  annotate("rect", xmin = -54.45/60, xmax = 0, ymin = -Inf, ymax = Inf, 
           alpha = .5,fill = "grey") +
  geom_vline(xintercept = 3, linetype = "dashed", size = 1.2, col = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2, col = "grey") +
  geom_line(data=dataUsed,aes(x=time_sec/60,y=(RR_mean),group=masktype,col=masktype),size=1.5)  + 
  geom_ribbon(data=dataUsed,aes(x=time_sec/60,ymin=RR_mean-(RR_std/1),ymax=RR_mean+(RR_std/1),fill=masktype),alpha=0.3,colour = NA)+
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10),
                     limits = xlimits2use) +
  scale_y_continuous(position = "left",
                     limits = c(-10, 4))+
  labs(y = expression(atop(paste(Delta*"RR"),"[1/min]")), x = "time [min]") + 
  scale_fill_manual(values=c(customred,customblue)) +
  scale_color_manual(values=c(customred,customblue)) +
  annotate("rect", xmin = 630.49/60, xmax = 600/60, ymin = -Inf, ymax = Inf, 
           alpha = .7, fill = "white") +
  theme_light()   +
  theme(text = element_text(size = 35, family="Myriad pro"), 
        axis.title.y=element_text(size=rel(0.9)), legend.position = "none",
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black")) 

tco2plot <- ggplot(data=dataUsed,aes(x=time_sec/60,y=TCO2_mean,col=masktype) ) +
  annotate("rect", xmin = -54.45/60, xmax = 0, ymin = -Inf, ymax = Inf, 
           alpha = .5,fill = "grey") +
  geom_vline(xintercept = 3, linetype = "dashed", size = 1.2, col = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2, col = "grey") +
  geom_line(data=dataUsed,aes(x=time_sec/60,y=(TCO2_mean),group=masktype,col=masktype),size=1.5)  + 
  geom_ribbon(data=dataUsed,aes(x=time_sec/60,ymin=TCO2_mean-(TCO2_std/1),ymax=TCO2_mean+(TCO2_std/1),fill=masktype),alpha=0.3,colour = NA)+
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10),
                     limits = xlimits2use) +
  scale_y_continuous(position = "left",
                     limits = c(-2, 1.7))+
  labs(y = expression(atop(paste(Delta*"TcCO"[2]),"[mmHg]")), x = "time [min]") + 
  scale_fill_manual(values=c(customred,customblue)) +
  scale_color_manual(values=c(customred,customblue)) +
  annotate("rect", xmin = 630.49/60, xmax = 600/60, ymin = -Inf, ymax = Inf, 
           alpha = .7, fill = "white") +
  theme_light()   +
  theme(text = element_text(size = 35, family="Myriad pro"), 
        axis.title.y=element_text(size=rel(0.9)), legend.position = "none",
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black")) 

spo2plot <- ggplot(data=dataUsed,aes(x=time_sec/60,y=SPO2_mean,col=masktype) ) +
  annotate("rect", xmin = -54.45/60, xmax = 0, ymin = -Inf, ymax = Inf, 
           alpha = .5,fill = "grey") +
  geom_vline(xintercept = 3, linetype = "dashed", size = 1.2, col = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2, col = "grey") +
  geom_line(data=dataUsed,aes(x=time_sec/60,y=(SPO2_mean),group=masktype,col=masktype),size=1.5)  + 
  geom_ribbon(data=dataUsed,aes(x=time_sec/60,ymin=SPO2_mean-(SPO2_std/1),ymax=SPO2_mean+(SPO2_std/1),fill=masktype),alpha=0.3,colour = NA)+
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10),
                     limits = xlimits2use) +
  scale_y_continuous(position = "left",
                     limits = c(-1.35, .85))+
  labs(y = expression(atop(paste(Delta*"SpO"[2]),"[%]")), x = "time [min]") + 
  scale_fill_manual(values=c(customred,customblue)) +
  scale_color_manual(values=c(customred,customblue)) +
  annotate("rect", xmin = 630.49/60, xmax = 600/60, ymin = -Inf, ymax = Inf, 
           alpha = .7, fill = "white") +
  theme_light()   +
  theme(text = element_text(size = 35, family="Myriad pro"), 
        axis.title.y=element_text(size=rel(0.9)), legend.position = "none",
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black")) 

mapplot <- ggplot(data=dataUsed,aes(x=time_sec/60,y=MAP_mean,col=masktype) ) +
  annotate("rect", xmin = -54.45/60, xmax = 0, ymin = -Inf, ymax = Inf, 
           alpha = .5,fill = "grey") +
  geom_vline(xintercept = 3, linetype = "dashed", size = 1.2, col = "magenta") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2, col = "grey") +
  geom_line(data=dataUsed,aes(x=time_sec/60,y=(MAP_mean),group=masktype,col=masktype),size=1.5)  + 
  geom_ribbon(data=dataUsed,aes(x=time_sec/60,ymin=MAP_mean-(MAP_std/1),ymax=MAP_mean+(MAP_std/1),fill=masktype),alpha=0.3,colour = NA)+
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10),
                     limits = xlimits2use) +
  scale_y_continuous(position = "left",
                     limits = c(-10, 13))+
  labs(y = expression(atop(paste(Delta*"MAP"),"[mmHg]")), x = "time [min]") + 
  scale_fill_manual(values=c(customred,customblue)) +
  scale_color_manual(values=c(customred,customblue)) +
  annotate("rect", xmin = 630.49/60, xmax = 600/60, ymin = -Inf, ymax = Inf, 
           alpha = .7, fill = "white") +
  theme_light()   +
  theme(text = element_text(size = 35, family="Myriad pro"), 
        axis.title.y=element_text(size=rel(0.9)), legend.position = "none",
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black")) 



cbfplot2 <- cbfplot +
  geom_segment(aes(x=-7,xend=-54.45/60,y=30,yend=30),size = 1.5, col = "black") +
  geom_segment(aes(x=0,xend=10,y=30,yend=30),size = 1.5, col = "black") 


# Putting the single plots together
fig4pnas <- ggarrange(cbfplot2 + rremove("xlab") + rremove("x.text") + rremove("x.axis") + rremove("x.ticks"), 
          sto2plot  + rremove("xlab") + rremove("x.text") + rremove("x.axis") + rremove("x.ticks"), 
          thcplot  + rremove("xlab") + rremove("x.text") + rremove("x.axis") + rremove("x.ticks"), 
          oefplot  + rremove("xlab") + rremove("x.text") + rremove("x.axis") + rremove("x.ticks"), 
          cmro2plot  + rremove("xlab") + rremove("x.text") + rremove("x.axis") + rremove("x.ticks"),
          hrplot  + rremove("xlab") + rremove("x.text") + rremove("x.axis") + rremove("x.ticks"), 
          rrplot  + rremove("xlab") + rremove("x.text") + rremove("x.axis") + rremove("x.ticks"),
          tco2plot  + rremove("xlab") + rremove("x.text") + rremove("x.axis") + rremove("x.ticks"),
          spo2plot  + rremove("xlab") + rremove("x.text") + rremove("x.axis") + rremove("x.ticks"),
          mapplot   + rremove("xlab") + rremove("x.axis") + rremove("x.ticks"), 
          ncol = 1, nrow = 10,
          align = "v")

fig4pnas2 <- annotate_figure(fig4pnas,
                bottom = text_grob("            time [min]", color = "black", size = 38, family="Myriad pro"))

fig4pnas3 <- annotate_figure(fig4pnas2,
                             top = text_grob(" Mask off                Mask on", color = "black", size = 45,  family="Myriad pro"))

#windows()
fig4pnas3
#Figure was saved with 1303 x 2422 pixels/ 13.57 x 25.23 inch
