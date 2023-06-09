require(readxl)
require(ggplot2)
require(gridExtra)

vga <- read.csv('Processed_data/Valve_gape_Absorption.csv')

vga.dec2 <- vga[c(2:9,13,15)]
library(tidyr)
vga.dec2 <- vga.dec2 %>% gather(Time, Valve_gape, -c(Trial, Mussel, Treatment, Interval, Start_tone, Mussel_size, Max_open, Mussel_presence))
vga.dec2 <- vga.dec2 %>% separate(Time, c(NA, "Time"), "Mean_gape_t")
vga.dec2$Time <- as.numeric(vga.dec2$Time)
vga.dec2$Absorption[vga.dec2$Time == 60] <- vga$Mean_abs_t60 - vga$Mean_abs_t0
vga.dec2$Absorption[vga.dec2$Time == 95] <- vga$Mean_abs_t95 - vga$Mean_abs_t60
vga.dec2$Absorption_bsl <- c(vga$Mean_abs_t0, vga$Mean_abs_t0)
vga.dec2$Absorption_10min[vga.dec2$Time == 60] <- (vga.dec2$Absorption[vga.dec2$Time == 60]/60)*10
vga.dec2$Absorption_10min[vga.dec2$Time == 95] <- (vga.dec2$Absorption[vga.dec2$Time == 95]/35)*10
vga.dec2 <- subset(vga.dec2, vga.dec2$Mussel_presence == F | 
                     (vga.dec2$Mussel_presence == T & vga.dec2$Max_open > 0))

# Feeding rate - mussel presence
p.abs2 <- ggplot(data = vga.dec2, aes(x = as.factor(Time), y = Absorption_10min, fill = Mussel_presence))+
  geom_hline(yintercept=0, linetype = 'dashed', color = '#696969')+
  geom_boxplot(outlier.shape = 1)+
  stat_summary(fun.y = "mean", geom="point", color = "red", size = 3, shape = 18,
               aes(group=Mussel_presence), position=position_dodge(.75), show.legend = FALSE)+
  ggtitle('\\u0394 phytoplankton\\nby mussel presence')+
  ylab(bquote('Change in absorption (AU·10 min'^-1*')'))+xlab("Time bin (min)")+
  scale_fill_manual(name = "Mussel\\npresent", values = c("#7CAE00", "#C77CFF"), 
                      labels = c(paste("\\nNo\\n(n = ", length(vga.dec2$Mussel[vga.dec2$Mussel_presence == F])/2, ")\\n", sep = ""), 
                                 paste("\\nYes\\n(n = ", length(vga.dec2$Mussel[vga.dec2$Mussel_presence == T])/2, ")\\n", sep = "")))+
  scale_x_discrete(labels=c("60" = "0-60", "95" = "60-95"))+
  #scale_y_continuous(limits = c(-.01,.006))+
  labs(tag = "a)")+theme_bw()+
  geom_text(aes(label = 'a c', x = .8, y = 0.023), size = 4)+
  geom_text(aes(label = 'b c', x = 1.2, y = 0.023), size = 4)+
  geom_text(aes(label = 'a d', x = 1.8, y = 0.023), size = 4)+
  geom_text(aes(label = 'b d', x = 2.2, y = 0.023), size = 4)
p.abs2

# stats on mussel presence effect on phytoplankton
vga.dec2$Time <- as.factor(vga.dec2$Time)
m1.2 <- glm(Absorption_10min ~ Mussel_presence * Time, data = vga.dec2, na.action = "na.fail")
require(MuMIn)
dredge(m1.2, rank = "AICc")
m1.2.b <- glm(Absorption_10min ~ Mussel_presence + Time, data = vga.dec2, na.action = "na.fail")
summary(m1.2.b)
require(effectsize)
standardize_parameters(m1.2.b)
citation("effectsize")

# Feeding rate - sound
p.abs2.treat <- ggplot(data = subset(vga.dec2, vga.dec2$Mussel_presence == T & vga.dec2$Max_open > .5), aes(x = as.factor(Time), y = Absorption_10min, fill = Treatment))+
  geom_hline(yintercept=0, linetype = 'dashed', color = '#696969')+
  geom_boxplot(outlier.shape = 1)+
  stat_summary(fun.y="mean", geom="point", color = "red", size = 3, shape = 18,
               aes(group=Treatment), position=position_dodge(.75), show.legend = FALSE)+
  ggtitle('\\u0394 phytoplankton\\nby sound exposure')+
  ylab(bquote('Change in absorption (AU·10 min'^-1*')'))+xlab("Time (min)")+
  scale_fill_discrete(name = "Treatment", 
                      labels = c(paste("\\nControl\\n(n = ", length(vga.dec2$Mussel[vga.dec2$Mussel_presence == T & vga.dec2$Treatment == 'controle' & vga.dec2$Max_open > 0.5])/2, ")\\n", sep = ""), 
                                 paste("\\nExposure\\n(n = ", length(vga.dec2$Mussel[vga.dec2$Mussel_presence == T & vga.dec2$Treatment == 'exposure' & vga.dec2$Max_open > 0.5])/2, ")\\n", sep = "")))+
  scale_x_discrete(labels=c("60" = "0-60", "95" = "60-95"))+
  #scale_y_continuous(limits = c(-.01,.006))+
  labs(tag = "b)")+theme_bw()+
  geom_text(aes(label = ' ', x = .8, y = 0.023), size = 4)
p.abs2.treat

# Stats on filtration effects by sound treatment
vga.dec2.p <- subset(vga.dec2, vga.dec2$Mussel_presence == T & !is.na(vga.dec2$Valve_gape) & vga.dec2$Valve_gape > .5)
## stats on mussel presence
hist(vga.dec2.p$Absorption, breaks = 15)
#vga.dec2.p$Time <- as.factor(vga.dec2.p$Time)
m2.2 <- glm(Absorption_10min ~ Time * Treatment + Mussel_size + Time + Absorption_bsl + Valve_gape, data = vga.dec2.p, family = "gaussian", na.action = "na.fail")
require(MuMIn)
dredge(m2.2, rank = "AICc")
m2.2.b <- glm(Absorption_10min ~ Mussel_size + Absorption_bsl + Treatment, data = vga.dec2.p, family = "gaussian", na.action = "na.fail")
#plot(m2.2.b)
#shapiro.test(residuals(m2.2.b))
summary(m2.2.b)
standardize_parameters(m2.2.b)

# Valve gape
p.gape.treat <- ggplot(data = subset(vga.dec2, vga.dec2$Mussel_presence == T & vga.dec2$Max_open > .5), aes(x = as.factor(Time), y = Valve_gape, fill = Treatment))+
  geom_boxplot(outlier.shape = 1)+
  stat_summary(fun.y="mean", geom="point", color = "red", size = 3, shape = 18,
               aes(group=Treatment), position=position_dodge(.75), show.legend = FALSE)+
  ggtitle('Mussels\\' valve gape\\nby sound exposure')+ylab('Mean valve gape (mm)')+xlab("Time (min)")+
  scale_fill_discrete(name = "Treatment", 
                      labels = c(paste("\\nControl\\n(n = ", length(vga.dec2$Mussel[vga.dec2$Mussel_presence == T & vga.dec2$Treatment == 'controle' & vga.dec2$Max_open > .5])/2, ")\\n", sep = ""), 
                                 paste("\\nExposure\\n(n = ", length(vga.dec2$Mussel[vga.dec2$Mussel_presence == T & vga.dec2$Treatment == 'exposure' & vga.dec2$Max_open > .5])/2, ")\\n", sep = "")))+
  scale_x_discrete(labels=c("60" = "0-60", "95" = "60-95"))+
  labs(tag = "c)")+theme_bw()+
  geom_text(aes(label = 'a c', x = .8, y = 5.3), size = 4)+
  geom_text(aes(label = 'b c', x = 1.2, y = 5.3), size = 4)+
  geom_text(aes(label = 'a d', x = 1.8, y = 5.3), size = 4)+
  geom_text(aes(label = 'b d', x = 2.2, y = 5.3), size = 4)
p.gape.treat

# Stats on valve gape effects by sound treatment
vga.dec2.p <- subset(vga.dec2, vga.dec2$Mussel_presence == T & !is.na(vga.dec2$Valve_gape))
## stats on mussel presence
hist(vga.dec2.p$Valve_gape, breaks = 15)
vga.dec2.p$Time <- as.factor(vga.dec2.p$Time)
m3 <- glm(Valve_gape ~ Time * Treatment + Mussel_size + Time, data = vga.dec2.p, na.action = "na.fail")
require(MuMIn)
dredge(m3, rank = "AICc")
m3.b <- glm(Valve_gape ~ Time + Treatment + Mussel_size, data = vga.dec2.p, na.action = "na.fail")
#plot(m3.b)
#shapiro.test(residuals(m3))
summary(m3.b)
standardize_parameters(m3.b)

## Make combined plot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend_mussel_p <-g_legend(p.abs2+theme(legend.position = "bottom"))
legend_sound <-g_legend(p.gape.treat+theme(legend.position = "bottom"))

overview.p <- grid.arrange(p.abs2+theme(legend.position = "none"), 
                         p.abs2.treat+theme(legend.position = "none"), 
                         p.gape.treat+theme(legend.position = "none"), 
                         nrow = 1, widths = c(1/3, 1/3, 1/3))
overview.l <- grid.arrange(legend_mussel_p, 
                           legend_sound,
                           nrow = 1, widths = c(.4, .6))
overview <- grid.arrange(overview.p, overview.l, heights = c(.85, .15))
ggsave('Figures/Fig_5abc_HQ.jpg', overview, png(), units = "in", width = 8.3, height = 5, dpi = 600)

