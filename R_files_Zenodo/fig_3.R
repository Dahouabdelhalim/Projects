### Figures to display the exposure schedules
require(ggplot2)

# Create DF
df <- data.frame('Time' = c((0*60):(100*60))/60)
df$Exposure <- NA
df$Frequency <- NA
exposure_no <- 1
interval <- 40
pulse <- 12
for(i in seq((25),(5420/60),((interval+5)/60))){ # 5 = actual pulse duration
  df$Exposure[df$Time >= (i-(pulse*.5/60)) & df$Time < (i+((pulse*.5)/60))] <- exposure_no
  exposure_no <- exposure_no + 1
}
df$Frequency[df$Exposure < max(df$Exposure, na.rm = T)] <- '200 Hz'
df$Frequency[df$Exposure == max(df$Exposure, na.rm = T)] <- '350 Hz'
df$Frequency[is.na(df$Exposure)] <- 'Silence'

p1 <- ggplot(data = df, aes(x = Time, y = 1, color = as.factor(Frequency)))+
  geom_point(shape = "|")+
  scale_color_manual(breaks = c("200 Hz", "350 Hz", "Silence"),
                     values=c("#F8766D", "#00BA38", "#d7e1e0"), name = "Exposure")+
  scale_x_continuous(limits = c(0,95), breaks = seq(0,100,10), expand = c(0.01, 0.01))+
  labs(subtitle = "a) Pulse train - Interval: 40 s; Start frequency: 200 Hz")+
  geom_vline(xintercept=0, linetype="solid")+
  geom_vline(xintercept=60, linetype="solid")+
  geom_vline(xintercept=95, linetype="solid")+
  theme_bw()+theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.text.x=element_blank())+
  guides(color = guide_legend(override.aes = list(shape = 15)))

df$Frequency[df$Exposure < max(df$Exposure, na.rm = T)] <- '350 Hz'
df$Frequency[df$Exposure == max(df$Exposure, na.rm = T)] <- '200 Hz'

p2 <- ggplot(data = df, aes(x = Time, y = 1, color = as.factor(Frequency)))+
  geom_point(shape = "|")+
  scale_color_manual(breaks = c("200 Hz", "350 Hz", "Silence"),
                     values=c("#F8766D", "#00BA38", "#d7e1e0"), name = "Exposure")+
  scale_x_continuous(limits = c(0,95), breaks = seq(0,100,10), expand = c(0.01, 0.01))+
  labs(subtitle = "b) Pulse train - Interval: 40 s; Start frequency: 350 Hz")+
  geom_vline(xintercept=0, linetype="solid")+
  geom_vline(xintercept=60, linetype="solid")+
  geom_vline(xintercept=95, linetype="solid")+
  theme_bw()+theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.text.x=element_blank())+
  guides(color = guide_legend(override.aes = list(shape = 15)))

df$Frequency[df$Exposure < max(df$Exposure, na.rm = T)] <- 'Silence'
df$Frequency[df$Exposure == max(df$Exposure, na.rm = T)] <- '350 Hz'
df$Time <- df$Time - 0.25

p5 <- ggplot(data = df, aes(x = Time, y = 1, color = as.factor(Frequency)))+
  geom_point(shape = "|")+
  scale_color_manual(breaks = c("350 Hz", "Silence"),
                     values=c("#00BA38", "#d7e1e0"), name = "Exposure")+
  scale_x_continuous(limits = c(0,95), breaks = seq(0,100,10), expand = c(0.01, 0.01))+
  labs(subtitle = "e) Control - Last frequency: 350 Hz")+
  geom_vline(xintercept=0, linetype="solid")+
  geom_vline(xintercept=60, linetype="solid")+
  geom_vline(xintercept=95, linetype="solid")+
  theme_bw()+theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.text.x=element_blank())+
  guides(color = guide_legend(override.aes = list(shape = 15)))

df$Frequency[df$Exposure == max(df$Exposure, na.rm = T)] <- '200 Hz'

p6 <- ggplot(data = df, aes(x = Time, y = 1, color = as.factor(Frequency)))+
  geom_point(shape = "|")+
  scale_color_manual(breaks = c("200 Hz", "Silence"),
                     values=c("#F8766D", "#d7e1e0"), name = "Exposure")+
  scale_x_continuous(limits = c(0,95), breaks = seq(0,100,10), expand = c(0.01, 0.01),
                     name = "Time (min)")+
  labs(subtitle = "f) Control - Last frequency: 200 Hz")+
  geom_vline(xintercept=0, linetype="solid")+
  geom_vline(xintercept=60, linetype="solid")+
  geom_vline(xintercept=95, linetype="solid")+
  theme_bw()+theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   #axis.title.x=element_blank(),
                   #axis.text.x=element_blank()
                   )+
  guides(color = guide_legend(override.aes = list(shape = 15)))

# Create DF
df <- data.frame('Time' = c((0*60):(100*60))/60)
df$Exposure <- NA
df$Frequency <- NA
exposure_no <- 1
interval <- 160
for(i in seq((25),(5461/60),((interval+5)/60))){
  df$Exposure[df$Time >= i-(pulse*.5/60) & df$Time < (i+(pulse*.5/60))] <- exposure_no
  exposure_no <- exposure_no + 1
}
df$Frequency[df$Exposure < max(df$Exposure, na.rm = T)] <- '200 Hz'
df$Frequency[df$Exposure == max(df$Exposure, na.rm = T)] <- '350 Hz'
df$Frequency[is.na(df$Exposure)] <- 'Silence'

p3 <- ggplot(data = df, aes(x = Time, y = 1, color = as.factor(Frequency)))+
  geom_point(shape = "|")+
  scale_color_manual(breaks = c("200 Hz", "350 Hz", "Silence"),
                     values=c("#F8766D", "#00BA38", "#d7e1e0"), name = "Exposure")+
  scale_x_continuous(limits = c(0,95), breaks = seq(0,100,10), expand = c(0.01, 0.01))+
  labs(subtitle = "c) Pulse train - Interval: 160 s; Start frequency: 200 Hz")+
  geom_vline(xintercept=0, linetype="solid")+
  geom_vline(xintercept=60, linetype="solid")+
  geom_vline(xintercept=95, linetype="solid")+
  theme_bw()+theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.text.x=element_blank())+
  guides(color = guide_legend(override.aes = list(shape = 15)))

df$Frequency[df$Exposure < max(df$Exposure, na.rm = T)] <- '350 Hz'
df$Frequency[df$Exposure == max(df$Exposure, na.rm = T)] <- '200 Hz'

p4 <- ggplot(data = df, aes(x = Time, y = 1, color = as.factor(Frequency)))+
  geom_point(shape = "|")+
  scale_color_manual(breaks = c("200 Hz", "350 Hz", "Silence"),
                     values=c("#F8766D", "#00BA38", "#d7e1e0"), name = "Exposure")+
  scale_x_continuous(limits = c(0,95), breaks = seq(0,100,10), expand = c(0.01, 0.01),
                     name = 'Time (min)')+
  labs(subtitle = "d) Pulse train - Interval: 160 s; Start frequency: 350 Hz")+
  geom_vline(xintercept=0, linetype="solid")+
  geom_vline(xintercept=60, linetype="solid")+
  geom_vline(xintercept=95, linetype="solid")+
  theme_bw()+theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.text.x=element_blank(),
                   legend.position="bottom")+
  guides(color = guide_legend(override.aes = list(shape = 15)))
p4

require(gridExtra)
# Get legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(p4)

overview <- grid.arrange(p1+theme(legend.position="none"),
             p2+theme(legend.position="none"),
             p3+theme(legend.position="none"),
             p4+theme(legend.position="none"),
             p5+theme(legend.position="none"),
             p6+theme(legend.position="none"), nrow = 6, heights=c(rep(0.1515,5),.2425))
overview2 <- grid.arrange(overview, mylegend, nrow = 2, heights=c(.9,0.1), top = 'Sound treatments')
ggsave('Figures/Fig_3_HQ_no_symbol.jpg', overview2, png(), units = "in", width = 8.3, height = 4.2, dpi = 600)

