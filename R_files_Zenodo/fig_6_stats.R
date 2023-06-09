######
## Valve gape reaction to sound
######

require(tidyr)
require(ggplot2)
require(gridExtra)

vo <- read.csv('Processed_data/Valve_opening_sound.csv')[2:15]
vo <- subset(vo, vo$Max_open > 0.5)
vo.l <- vo %>% gather(When, Valve_gape, -c(Trial, Mussel, Treatment, Interval, Start_tone, Mussel_size, Max_open))
vo.l$Exposure <- substr(vo.l$When, 15, 17)
vo.l$When <- substr(vo.l$When, 11, 13)
vo.l$Exposure[vo.l$Exposure == 'fst'] <- 'First'
vo.l$Exposure[vo.l$Exposure == 'sec'] <- 'Second last'
vo.l$Exposure[vo.l$Exposure == 'lst'] <- 'Last'
vo.l$When[vo.l$When == 'bef'] <- 'Before'
vo.l$When[vo.l$When == 'aft'] <- 'After'
vo.l$When <- as.factor(vo.l$When)
vo.l$When <- factor(vo.l$When, levels = c("Before", "After"))

p.gape.p1 <- ggplot(data = subset(vo.l, vo.l$Exposure == 'First' & vo.l$Treatment == 'exposure' & vo.l$Interval >= 20), 
                    aes(x = When, y = Valve_gape))+
  geom_boxplot(outlier.shape = 1, fill = '#00BFC4')+
  stat_summary(fun.y = "mean", geom="point", color = "red", size = 3, shape = 18, show.legend = FALSE)+
  ggtitle('Valve gape\\nFirst exposure')+
  ylab('Valve gape (mm)')+xlab('Period')+
  scale_y_continuous(limits = c(0,6))+
  labs(tag = "a)")+theme_bw()+
  geom_text(aes(label = '*', x = 1.5, y = 5.6), size = 9)
p.gape.p1

p.gape.pm <- ggplot(data = subset(vo.l, vo.l$Exposure == 'Second last' & vo.l$Treatment == 'exposure' & vo.l$Interval >= 20), 
                    aes(x = When, y = Valve_gape))+
  geom_boxplot(outlier.shape = 1, fill = '#00BFC4')+
  stat_summary(fun.y="mean", geom="point", color = "red", size = 3, shape = 18, show.legend = FALSE)+
  ggtitle('\\nSecond last exposure')+
  ylab(' ')+xlab('Period')+
  scale_y_continuous(limits = c(0,6))+
  labs(tag = "b)")+guides(fill = F)+theme_bw()
#p.gape.pm

p.gape.p2 <- ggplot(data = subset(vo.l, vo.l$Exposure == 'Last' & vo.l$Treatment == 'exposure' & vo.l$Interval >= 20), 
                    aes(x = When, y = Valve_gape))+
  geom_boxplot(outlier.shape = 1, fill = '#00BFC4')+ # 
  stat_summary(fun.y="mean", geom="point", color = "red", size = 3, shape = 18, show.legend = FALSE)+
  ggtitle('\\nLast exposure')+
  ylab(' ')+xlab('Period')+
  scale_y_continuous(limits = c(0,6))+
  labs(tag = "c)")+
  theme_bw()+
  geom_text(aes(label = "+", x = 1.5, y = 5.8), size = 6)
p.gape.p2


overview <- grid.arrange(p.gape.p1, p.gape.pm, p.gape.p2, nrow = 1, widths = c((1/3), (1/3), (1/3)))
ggsave('Figures/Fig_6abc_HQ.jpg', overview, png(), units = "in", width = 8.3, height = 4, dpi = 600)

# Stats
vo.l1 <- subset(vo.l, vo.l$Exposure == 'First' & vo.l$Treatment == 'exposure' & vo.l$Interval >= 20)
hist(vo.l1$Valve_gape, breaks = 23)
m1 <- glm(Valve_gape ~ When, data = vo.l1, family = 'gaussian')
#plot(m1)
summary(m1)
require(effectsize)
standardize_parameters(m1)

vo.l2 <- subset(vo.l, vo.l$Exposure == 'Second last' & vo.l$Treatment == 'exposure' & vo.l$Interval >= 20)
hist(vo.l1$Valve_gape, breaks = 23)
m2 <- glm(Valve_gape ~ When, data = vo.l2, family = 'gaussian')
#plot(m2)
summary(m2)
standardize_parameters(m2)

vo.l3 <- subset(vo.l, vo.l$Exposure == 'Last' & vo.l$Treatment == 'exposure' & vo.l$Interval >= 20)
hist(vo.l3$Valve_gape, breaks = 23)
m3 <- glm(Valve_gape ~ When, data = vo.l3, family = 'gaussian')
#plot(m3)
summary(m3)
standardize_parameters(m3)

# sample size
nrow(subset(vo.l, vo.l$Exposure == 'Last' & vo.l$Treatment == 'exposure' & vo.l$Interval >= 20))/2