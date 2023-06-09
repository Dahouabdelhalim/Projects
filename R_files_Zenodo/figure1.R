###########################################################################################
# load packages
library(rethinking)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

source("Code/aft_calculate_predictions.R")

##############################################################
# This section can be adjusted as needed #####################
##############################################################

# adjust plot parameters for prediction plot
shade_alpha = 0.1 # transparency of HPDI 
y <- 0.6 # y-axis limit

line_width <- 1 # line width
color1 <- "gold"
color2 <- "cadetblue3"
color3 <- "coral1"
shade <- "burlywood"

# axes titles
axes_title_size <- 12  
panel_label_size <- 5
##############################################################

# PLOT
pdf("Outputs/Figure1.pdf",width=8, height=12 )

#GROOMING GIVING
p1gg <- ggplot(NULL, aes(x=soc_nat[[1]], y=agg_mean ) )+
  geom_line(color=color3, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=agg_HPDI[1,], ymax=agg_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("Probability of Dying")+
  scale_x_continuous(name ="", 
                     breaks=c(round( min(soc_nat[[1]]),0): round( max(soc_nat[[1]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))
  
# plot for grooming giving predictions
p2gg <- ggplot(NULL, aes(x=soc_nat[[6]], y=fgg_mean ) )+
  geom_line(color=color1, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=fgg_HPDI[1,], ymax=fgg_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("")+
  scale_x_continuous(name ="Grooming giving posterior mean", 
                     breaks=c(round( min(soc_nat[[6]]),0): round( max(soc_nat[[6]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))


p3gg <- ggplot(NULL, aes(x=soc_nat[[11]], y=mgg_mean ) )+
  geom_line(color=color2, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=mgg_HPDI[1,], ymax=mgg_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("")+
  scale_x_continuous(name ="", 
                     breaks=c(round( min(soc_nat[[11]]),0): round( max(soc_nat[[11]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))

# GROOMING RECEIVING
p1gr <- ggplot(NULL, aes(x=soc_nat[[2]], y=agr_mean ) )+
  geom_line(color=color3, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=agr_HPDI[1,], ymax=agr_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("Probability of Dying")+
  scale_x_continuous(name ="", 
                     breaks=c( round( min(soc_nat[[2]]),0): round( max(soc_nat[[2]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))

p2gr <- ggplot(NULL, aes(x=soc_nat[[7]], y=fgr_mean ) )+
  geom_line(color=color1, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=fgr_HPDI[1,], ymax=fgr_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("")+
  scale_x_continuous(name ="Grooming receiving posterior mean", 
                     breaks=c(round( min(soc_nat[[7]]),0): round( max(soc_nat[[7]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))


p3gr <- ggplot(NULL, aes(x=soc_nat[[12]], y=mgr_mean ) )+
  geom_line(color=color2, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=mgr_HPDI[1,], ymax=mgr_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("")+
  scale_x_continuous(name ="", 
                     breaks=c(round( min(soc_nat[[12]]),0): round( max(soc_nat[[12]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))

# SUPPORT GIVING
p1sg <- ggplot(NULL, aes(x=soc_nat[[3]], y=asg_mean ) )+
  geom_line(color=color3, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=asg_HPDI[1,], ymax=asg_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("Probability of Dying")+
  scale_x_continuous(name ="", 
                     breaks=c(round( min(soc_nat[[3]]),0): round( max(soc_nat[[3]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))

p2sg <- ggplot(NULL, aes(x=soc_nat[[8]], y=fsg_mean ) )+
  geom_line(color=color1, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=fsg_HPDI[1,], ymax=fsg_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("")+
  scale_x_continuous(name ="Support giving posterior mean", 
                     breaks=c(round( min(soc_nat[[8]]),0): round( max(soc_nat[[8]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))


p3sg <- ggplot(NULL, aes(x=soc_nat[[13]], y=msg_mean ) )+
  geom_line(color=color2, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=msg_HPDI[1,], ymax=msg_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("")+
  scale_x_continuous(name ="", 
                     breaks=c(round( min(soc_nat[[13]]),0): round( max(soc_nat[[13]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))


# SUPPORT RECEIVING
p1sr <- ggplot(NULL, aes(x=soc_nat[[4]], y=asr_mean ) )+
  geom_line(color=color3, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=asr_HPDI[1,], ymax=asr_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("Probability of Dying")+
  scale_x_continuous(name ="", 
                     breaks=c(round( min(soc_nat[[4]]),0): round( max(soc_nat[[4]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))

p2sr <- ggplot(NULL, aes(x=soc_nat[[9]], y=fsr_mean ) )+
  geom_line(color=color1, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=fsr_HPDI[1,], ymax=fsr_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("")+
  scale_x_continuous(name ="Support receiving posterior mean", 
                     breaks=c(round( min(soc_nat[[9]]),0): round( max(soc_nat[[9]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))


p3sr <- ggplot(NULL, aes(x=soc_nat[[14]], y=msr_mean ) )+
  geom_line(color=color2, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=msr_HPDI[1,], ymax=msr_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("")+
  scale_x_continuous(name ="", 
                     breaks=c(round( min(soc_nat[[14]]),0): round( max(soc_nat[[14]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))


# FORAGING
p1f <- ggplot(NULL, aes(x=soc_nat[[5]], y=af_mean ) )+
  geom_line(color=color3, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=af_HPDI[1,], ymax=af_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("Probability of Dying")+
  scale_x_continuous(name ="", 
                     breaks=c(round( min(soc_nat[[5]]),0): round( max(soc_nat[[5]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))

# plot for grooming giving predictions
p2f <- ggplot(NULL, aes(x=soc_nat[[10]], y=ff_mean ) )+
  geom_line(color=color1, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=ff_HPDI[1,], ymax=ff_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("")+
  scale_x_continuous(name ="Foraging posterior mean", 
                     breaks=c(round( min(soc_nat[[10]]),0): round( max(soc_nat[[10]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))


p3f <- ggplot(NULL, aes(x=soc_nat[[15]], y=mf_mean ) )+
  geom_line(color=color2, lwd=line_width)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_ribbon( aes(ymin=mf_HPDI[1,], ymax=mf_HPDI[2,]), alpha=shade_alpha)+
  ylim(0, y) +
  ylab("")+
  scale_x_continuous(name ="", 
                     breaks=c(round( min(soc_nat[[15]]),0): round( max(soc_nat[[15]]),0) ) )+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) +
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9) )+
  theme(axis.title.x = element_text(size = axes_title_size),
        axis.title.y = element_text(size = axes_title_size))

# arrange the plot
grid.arrange( arrangeGrob(
  p1gg + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ),
  p2gg + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ),
  p3gg + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ), 
  top=textGrob("Grooming giving",  x = 0, hjust = -0.27),  nrow=1),
  
  arrangeGrob(
    p1gr + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ),
    p2gr + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ),
    p3gr + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ), 
    top=textGrob("Grooming receiving",  x = 0, hjust = -0.23),  nrow=1),
  
  arrangeGrob(
    p1sg + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ),
    p2sg + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ),
    p3sg + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ), 
    top=textGrob("Support giving",  x = 0, hjust = -0.3), nrow=1),
  
  arrangeGrob(
    p1sr + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ),
    p2sr + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ),
    p3sr + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ), 
    top=textGrob("Support receiving", x = 0, hjust = -0.25),  nrow=1),
  
  arrangeGrob(
    p1f + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ), 
    p2f + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ),
    p3f + theme(axis.title.x = element_blank(),axis.title.y = element_blank() ), 
    top=textGrob("Foraging",  x = 0, hjust = -0.5),
    
    # global axes
    bottom = textGrob("Posterior Mean", rot = 0, vjust = 0),  nrow=1),
  left = textGrob("Probability of Dying", rot = 90, vjust = 1),
  nrow=5)


dev.off()

