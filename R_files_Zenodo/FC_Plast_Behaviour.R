# The genetic basis of variation in sexual aggression: evolution versus social plasticity

# FC Plasticity Behaviour Analysis
# Andrew M. Scott


library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(car)
library(grid)
library(gridExtra)
library(emmeans)

plast_behav <- read_csv("FC_plast_behav.csv")

plast_behav$Treatment <- as.factor(plast_behav$Treatment)
plast_behav$Observer <- as.factor(plast_behav$Observer)
plast_behav$`Mounting Duration` <- as.numeric(plast_behav$`Mounting Duration`)
plast_behav$`Pursuit Duration` <- as.numeric(plast_behav$`Pursuit Duration`)

## Mating analysis

plast_behav_mate <- glmmTMB(Mated ~ 1 + Treatment + (1 | Day) + Block + (1 | Dish), # Block changed to fixed effect to solve convergence issue
                            data = plast_behav, family = "binomial")

plot(simulateResiduals(plast_behav_mate)) 

summary(plast_behav_mate) 
car::Anova(plast_behav_mate) 

mate_emm <- emmeans(plast_behav_mate, "Treatment")
eff_size(mate_emm, sigma = sigma(plast_behav_mate), df.residual(plast_behav_mate))

## Behaviour analysis (Mounting and pursuit duration)

plast_behav_durations <- plast_behav[complete.cases(plast_behav), ]

plast_behav_mount <- glmmTMB(`Mounting Duration` ~ 1 + Observer + Treatment + 
                               (1 | Day) + Block + (1 | Dish), 
                             data = plast_behav_durations, family = tweedie())

plot(simulateResiduals(plast_behav_mount)) 

summary(plast_behav_mount) 
car::Anova(plast_behav_mount) 

mount_emm <- emmeans(plast_behav_mount, "Treatment")
eff_size(mount_emm, sigma = sigma(plast_behav_mount), df.residual(plast_behav_mount))

plast_behav_pursuit <- glmmTMB(`Pursuit Duration` ~ 1 + Observer + Treatment + 
                                 Day + Block + (1 | Dish), 
                               data = plast_behav_durations, family = tweedie())

plot(simulateResiduals(plast_behav_pursuit)) 

summary(plast_behav_pursuit) 
car::Anova(plast_behav_pursuit) 

pursuit_emm <- emmeans(plast_behav_pursuit, "Treatment")
eff_size(pursuit_emm, sigma = sigma(plast_behav_pursuit), df.residual(plast_behav_pursuit))

########
# Plots

mating_summarized <- data.frame(Treatment = c("EXP", "ISO"), Prop_mated = c(1/90, 6/90))
mating_summarized$StdErr <- sqrt((mating_summarized$Prop_mated*(1 - mating_summarized$Prop_mated))/90)

mating_plot <- ggplot(data = mating_summarized, 
                      aes(x = Treatment, 
                          y = Prop_mated, 
                          color = Treatment)) +
  geom_point(shape = 95, size = 20) +
  geom_errorbar(aes(ymin = Prop_mated - StdErr, ymax = Prop_mated + StdErr), width = 0.075) +
  scale_color_manual(limits = c("EXP", "ISO"), values = c("blue", "red")) +
  scale_x_discrete(limits = c("EXP", "ISO"), labels = c("Experienced", "Isolated")) +
  scale_y_continuous(limits = c(0, .4), breaks = c(0, 0.1, 0.2, 0.3)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16)) +
  labs(x = "Treatment",
       y = "Proportion males forcibly copulated")


mount_plot <- ggplot(data = plast_behav_durations, 
                     aes(x = Treatment,
                         y = (`Mounting Duration`/30),  # seconds per minute
                         color = Treatment)) +
  geom_violin() +
  geom_point(stat = "summary", fun = "mean", shape = 95, size = 10) +
  scale_color_manual(limits = c("EXP", "ISO"), values = c("blue", "red")) +
  scale_x_discrete(limits = c("EXP", "ISO"), labels = c("Experienced", "Isolated")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16)) +
  labs(x = "Treatment",
       y = "Mounting (s/min)")


pursuit_plot <- ggplot(data = plast_behav_durations, 
                       aes(x = Treatment,
                           y = (`Pursuit Duration`/30),  # seconds per minute
                           color = Treatment)) +
  geom_violin() +
  geom_point(stat = "summary", fun = "mean", shape = 95, size = 10) +
  scale_color_manual(limits = c("EXP", "ISO"), values = c("blue", "red")) +
  scale_x_discrete(limits = c("EXP", "ISO"), labels = c("Experienced", "Isolated")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16)) +
  labs(x = "Treatment",
       y = "Pursuit (s/min)")

blank <- grid.rect(gp=gpar(col="white"))
grid.arrange(arrangeGrob(blank,
                         mating_plot,
                         ncol = 2, 
                         widths = c(2/3, 1/3)),
             arrangeGrob(pursuit_plot + theme(axis.title.x = element_blank()), 
                         mount_plot + theme(axis.title.x = element_blank()),
                         ncol = 2, 
                         widths = c(1/2, 1/2),
                         bottom = textGrob("Treatment", gp = gpar(fontsize = 16))))


grid.arrange(arrangeGrob(blank,
                         mating_plot,
                         mount_plot + theme(axis.title.x = element_blank()), 
                         pursuit_plot + theme(axis.title.x = element_blank()),
                         ncol = 2, 
                         widths = c(2/3, 1/3, 1/2, 1/2),
                         bottom = textGrob("Treatment", gp = gpar(fontsize = 16))))


