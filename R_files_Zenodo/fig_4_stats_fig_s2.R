# Calibration and delta absorption

require(readxl)
require(ggplot2)

df <- read.csv('Meta_data/Absorption_calibration.csv')
df$concentration <- df$concentration/10000000 # to Liters

p.calib <- ggplot(data = df, aes(x = concentration, y = absorption))+
  geom_segment(aes(x = 0, xend = 24, y = 0.0148601, yend = (0.0148601 + 0.0089459*24)), color="blue", size = 1) +
  geom_point(shape = 1)+
  ggtitle('Calibration absorption measures')+xlab(expression(e^7~Cells/L))+ylab('Absorption (AU)')+
  ylim(0,0.23)+
  labs(tag = "a)")+theme_bw()
p.calib

hist(df$absorption)

m1 <- glm(absorption ~ concentration, data = df, family = "gaussian")
#plot(m1)
shapiro.test(residuals(m1))
summary(m1)
require(effectsize)
standardize_parameters(m1)

# Absorption
vga <- read.csv('Processed_data/Valve_gape_Absorption.csv')

# Link absorption to valve gape
vga$d.abs_tot <- vga$Mean_abs_t95 - vga$Mean_abs_t0
vga$Mean_gape_tot <- (vga$Mean_gape_t60*(60/95)) + (vga$Mean_gape_t95*(35/90))

p.vga3 <- ggplot(data = vga, aes(x = Mean_gape_tot, y = d.abs_tot))+
  #trendline for mussel-size = 4 cm
  geom_segment(aes(x = 0, xend = 4.4, y = 0.036212+(4 * -0.012178), yend = 0.036212+(4 * -0.012178) + -0.004013*4.4), color="blue", size = 1) +
  geom_point(shape = 1)+
  ggtitle("Absorption and valve gape (0-95 min)")+ylab('\\u0394 absorption (AU)')+xlab('Mean valve gape (mm)')+
  scale_y_continuous(#breaks = seq(0,0.25,0.05), 
    limits = c(min(c(vga$d.abs_tot,vga$d.abs_tot), na.rm = T),max(c(vga$d.abs_tot,vga$d.abs_tot), na.rm = T)))+
  labs(tag = "b)")+theme_bw()
p.vga3

hist(vga$d.abs_tot)

# Stats
vga.s <- subset(vga, !is.na(Mean_gape_tot))
m2.f <- glm(d.abs_tot ~ Mean_gape_tot * Mussel_size, data = vga.s, family = "gaussian", na.action = "na.fail")
require(MuMIn)
dredge(m2.f, rank = "AICc")
m2 <- glm(d.abs_tot ~ Mean_gape_tot + Mussel_size, data = vga.s, family = "gaussian", na.action = "na.fail")
#plot(m2)
shapiro.test(residuals(m2))
hist(residuals(m2))
summary(m2)
standardize_parameters(m2)

## Also determine r^2
require(rsq)
rsq(m2, adj = T, type = 'v')

require(gridExtra)
overview <- grid.arrange(p.calib, p.vga3, nrow = 1, widths = c(0.5, 0.5))
overview
ggsave('Figures/fig_4ab_HQ.jpg', overview, png(), units = "in", width = 8.3, height = 4, dpi = 600)


# Extra plots for supplement
p.s1 <- ggplot(data = subset(vga, vga$Mussel_size >= 3.5 & vga$Mussel_size < 4.0), aes(x = Mean_gape_tot, y = d.abs_tot))+
  geom_segment(aes(x = 0, xend = 4.23, y = -0.008953, yend = -0.008953 + (-0.002992 * 4.23)), color="blue", size = 1)+
  geom_point(shape = 1)+
  ylab('\\u0394 absorption (AU)')+
  xlab('')+
  ggtitle('Absorption and valve gape (0 - 95 min)\\nSize: 3.5 - 3.9 cm')+
  scale_y_continuous(#breaks = seq(0,0.25,0.05), 
    limits = c(min(c(vga$d.abs_tot,vga$d.abs_tot), na.rm = T),max(c(vga$d.abs_tot,vga$d.abs_tot), na.rm = T)))+
  scale_x_continuous(
    limits = c(0,max(vga$Mean_gape_tot, na.rm = T)))+
  labs(tag = "a)")+theme_bw()+theme(axis.title.x = element_blank())
p.s1

p.s2 <- ggplot(data = subset(vga, vga$Mussel_size >= 4 & vga$Mussel_size < 4.5), aes(x = Mean_gape_tot, y = d.abs_tot))+
  geom_segment(aes(x = 0, xend = 4.372, y = -0.016005, yend = -0.016005 + (-0.003525 * 4.372)), color="blue", size = 1)+
  geom_point(shape = 1)+
  ylab('')+
  xlab('')+
  ggtitle('\\nSize: 4.0 - 4.4 cm')+
  scale_y_continuous(#breaks = seq(0,0.25,0.05), 
    limits = c(min(c(vga$d.abs_tot,vga$d.abs_tot), na.rm = T),max(c(vga$d.abs_tot,vga$d.abs_tot), na.rm = T)))+
  scale_x_continuous(
    limits = c(0,max(vga$Mean_gape_tot, na.rm = T)))+
  labs(tag = "b)")+theme_bw()+theme(axis.title.x = element_blank())

p.s3 <- ggplot(data = subset(vga, vga$Mussel_size >= 4.5 & vga$Mussel_size < 5), aes(x = Mean_gape_tot, y = d.abs_tot))+
  geom_segment(aes(x = 0, xend = 3.63, y = -0.014994, yend = -0.014994 + (-0.007634 * 3.63)), color="blue", size = 1)+
  geom_point(shape = 1)+
  ylab('\\u0394 absorption (AU)')+
  xlab('Mean valve gape (mm)')+
  ggtitle('Size: 4.5 - 4.9 cm')+
  scale_y_continuous(#breaks = seq(0,0.25,0.05), 
    limits = c(min(c(vga$d.abs_tot,vga$d.abs_tot), na.rm = T),max(c(vga$d.abs_tot,vga$d.abs_tot), na.rm = T)))+
  scale_x_continuous(
    limits = c(0,max(vga$Mean_gape_tot, na.rm = T)))+
  labs(tag = "c)")+theme_bw()

p.s4 <- ggplot(data = subset(vga, vga$Mussel_size >= 5), aes(x = Mean_gape_tot, y = d.abs_tot))+
  geom_segment(aes(x = 0.05, xend = 3.4, y = -0.022328 + (-0.003801 * 0.05), yend = -0.022328 + (-0.003801 * 3.4)), color="blue", size = 1)+
  geom_point(shape = 1)+
  ylab('')+
  xlab('Mean valve gape (mm)')+
  ggtitle('Size: 5.0 - 6.0 cm')+
  scale_y_continuous(#breaks = seq(0,0.25,0.05), 
    limits = c(min(c(vga$d.abs_tot,vga$d.abs_tot), na.rm = T),max(c(vga$d.abs_tot,vga$d.abs_tot), na.rm = T)))+
  scale_x_continuous(
    limits = c(0,max(vga$Mean_gape_tot, na.rm = T)))+
  labs(tag = "d)")+theme_bw()

psizes <- grid.arrange(p.s1, p.s2, p.s3, p.s4, nrow = 2, widths = c(0.5, 0.5), heights = c(0.5, 0.5))
ggsave('Figures/Fig_s2a-d_HQ.jpg', psizes, png(), units = "in", width = 8.3, height = 6.5, dpi = 600)
