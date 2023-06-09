# The genetic basis of variation in sexual aggression: evolution versus social plasticity

# FC candidate validation
# Andrew M. Scott


library(tidyverse)
library(lme4)
library(car)
library(glmmTMB)
library(DHARMa)
library(gridExtra)
library(pbkrtest)
library(grid)
library(ggpubr)
library(MASS)

options(contrasts = c("contr.sum", "contr.poly"))

#######################################
#### Mating data

fc.data <- read_csv("FC_Candidate_Val_Mating.csv")

### Figures

fc.summarized <- fc.data %>%
  group_by(Gene_name, Treatment) %>%
  summarise(Trials = length(Mated),
            Num.mated = sum(Mated))

fc.summarized$Gene_name <- as.factor(fc.summarized$Gene_name)

fc.summarized$Prop.mated <- fc.summarized$Num.mated/fc.summarized$Trials

fc.summarized$Std.error <- sqrt((fc.summarized$Prop.mated*(1 - fc.summarized$Prop.mated))/fc.summarized$Trials)

# Overlap genes
fc_val_plot <- ggplot(data = subset(fc.summarized, Gene_name == "CG14153"|
                                      Gene_name == "Drsl4"|
                                      Gene_name == "GstZ1"|
                                      Gene_name == "Nepl18"|
                                      Gene_name == "verm"), 
                      aes(x = Treatment, 
                          y = Prop.mated,
                          color = Treatment,
                          fill = Treatment)) +
  geom_bar(stat = "identity", width = 0.75, fill = "white", lwd = 0.8) + 
  geom_errorbar(aes(ymin = Prop.mated - Std.error, ymax = Prop.mated + Std.error), width = 0.075) +
  facet_grid(~ Gene_name) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(limits = c("RNAi_elav", "RNAi_CyO", "TripCon_elav"),
                   labels = c("RNAi/elav-Gal4", "RNAi/+", "+/elav-Gal4")) +
  labs(x = "Genotype", 
       y = "Proportion mated") +
  scale_colour_manual(limits = c("RNAi_elav", "RNAi_CyO", "TripCon_elav"), values = c("black", "grey65", "grey65")) +
  scale_fill_manual(limits = c("RNAi_elav", "RNAi_CyO", "TripCon_elav"), values = c("black", "grey65", "grey65")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.6))

fc_val_plot

# Bonus genes (non-overlapped)
fc_val_plot.ex <- ggplot(data = subset(fc.summarized, Gene_name == "Lsp2"|
                                         Gene_name == "Nazo"), 
                         aes(x = Treatment, 
                             y = Prop.mated,
                             color = Treatment,
                             fill = Treatment)) +
  geom_bar(stat = "identity", width = 0.75, fill = "white", lwd = 0.8) + 
  geom_errorbar(aes(ymin = Prop.mated - Std.error, ymax = Prop.mated + Std.error), width = 0.075) +
  facet_grid(~ Gene_name) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(limits = c("RNAi_elav", "RNAi_CyO", "TripCon_elav"),
                   labels = c("RNAi/elav-Gal4", "RNAi/+", "+/elav-Gal4")) +
  labs(x = "Genotype", 
       y = "Proportion mated") +
  scale_colour_manual(limits = c("RNAi_elav", "RNAi_CyO", "TripCon_elav"), values = c("black", "grey65", "grey65")) +
  scale_fill_manual(limits = c("RNAi_elav", "RNAi_CyO", "TripCon_elav"), values = c("black", "grey65", "grey65")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.6))

fc_val_plot.ex

grid.arrange(arrangeGrob(fc_val_plot + theme(axis.title.x = element_blank()), 
                         fc_val_plot.ex + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
                         ncol = 2, 
                         bottom = textGrob("Genotype", gp = gpar(fontsize = 16)), 
                         widths = c(5, 2)))


### Analysis

fc.data$Day <- as.factor(fc.data$Day)
fc.data$Rack <- as.factor(fc.data$Rack)
fc.data$Treatment <- as.factor(fc.data$Treatment)

contrast1 <- c(-0.5, 1, -0.5) # compare trt to controls
contrast2 <- c(-1, 0, 1) # compare controls
contrasts <- rbind(contrast1, contrast2)

# compute generalized inverse to get contrast matrix
con.mat <- ginv(contrasts)
colnames(con.mat) <- c("Exp_vs_Controls", "Con_vs_con")

# apply to trt var
contrasts(fc.data$Treatment) <- con.mat


# Gene 1 - CG14153

g1.glm <- glmmTMB(Mated ~ 1 + Treatment + Day + (1 | Rack),
                  family = binomial(),
                  data = subset(fc.data, Gene == 1))
summary(g1.glm)

g1.glm.2 <- glmmTMB(Mated ~ 1 + Treatment + Day + (1 | Treatment:Day:Rack), 
                    family = binomial(),
                    data = subset(fc.data, Gene == 1))

#g1.glm.2 <- glmmTMB(Mated ~ 1 + Day + Rack*Treatment,
#                    family = binomial(),
#                    data = subset(fc.data, Gene == 1))
summary(g1.glm.2)

ggplot(data = fc.data, aes(y = Mated, x = as.numeric(Rack), 
                           color = Treatment)) +
  geom_smooth() +
  geom_point() +
  geom_jitter(height = 0.1, width = 0.05)


# Gene 2 - verm

g2.glm <- glmmTMB(Mated ~ 1 + Treatment + Day + (1 | Rack),
                  family = binomial(),
                  data = subset(fc.data, Gene == 2))
summary(g2.glm)

# Gene 5 - Drsl4

g5.glm <- glmmTMB(Mated ~ 1 + Treatment + Day + (1 | Rack),
                  family = binomial(),
                  data = subset(fc.data, Gene == 5))
summary(g5.glm)

# Gene 6 - GstZ1

g6.glm <- glmmTMB(Mated ~ 1 + Treatment + Day + (1 | Rack),
                  family = binomial(),
                  data = subset(fc.data, Gene == 6))
summary(g6.glm)

# Gene 7 - Nepl18

g7.glm <- glmmTMB(Mated ~ 1 + Treatment + Day + (1 | Rack),
                  family = binomial(),
                  data = subset(fc.data, Gene == 7))
summary(g7.glm)


## Extra genes

# Lsp2

g3.glm <- glmmTMB(Mated ~ 1 + Treatment + Day + (1 | Rack),
                  family = binomial(),
                  data = subset(fc.data, Gene == 3))
summary(g3.glm)

# Nazo

g4.glm <- glmmTMB(Mated ~ 1 + Treatment + (1 | Rack),
                  family = binomial(),
                  data = subset(fc.data, Gene == 4))
summary(g4.glm) # not sig - too low sample size


#############################################
#### Pursuit Data

pursuit.data <- read_csv("FC_Candidate_Val_Pursuit.csv")

pursuit.data$Treatment <- as.factor(pursuit.data$Treatment)
pursuit.data$Day <- as.factor(pursuit.data$Day)
pursuit.data$Rack <- as.factor(pursuit.data$Rack)
pursuit.data$Prop_trials_pursuing <- as.numeric(pursuit.data$Prop_trials_pursuing)

dim(pursuit.data)

# remove cases where there was a mating right at the start of the trial (i.e. no pursuit obs)
pursuit.data <- pursuit.data[complete.cases(pursuit.data), ]
dim(pursuit.data)

# Figure
pursuit.fig <- ggplot(data = subset(pursuit.data, Gene_name == "CG14153"|
                                      Gene_name == "Drsl4"|
                                      Gene_name == "GstZ1"|
                                      Gene_name == "Nepl18"|
                                      Gene_name == "verm"), 
                      aes(x = Treatment, 
                          y = Prop_trials_pursuing,
                          color = Treatment)) +
  geom_boxplot(width = 0.6, lwd = 0.8, outlier.shape = NA, fatten = 1) + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.9, alpha = 0.3, stroke = 0) +
  facet_grid(~ Gene_name) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(limits = c("RNAi_elav", "RNAi_CyO", "TripCon_elav"),
                   labels = c("RNAi/elav-Gal4", "RNAi/+", "+/elav-Gal4")) +
  labs(x = "Genotype", 
       y = "Proportion of observations with pursuit") +
  scale_colour_manual(limits = c("RNAi_elav", "RNAi_CyO", "TripCon_elav"), values = c("black", "grey65", "grey65")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.30), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))

pursuit.fig

# Extra genes
pursuit.fig.ex <- ggplot(data = subset(pursuit.data, Gene_name == "Lsp2"|
                                         Gene_name == "Nazo"), 
                         aes(x = Treatment, 
                             y = Prop_trials_pursuing,
                             color = Treatment)) +
  geom_boxplot(width = 0.6, lwd = 0.8, outlier.shape = NA, fatten = 1) + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.9, alpha = 0.3, stroke = 0) +
  facet_grid(~ Gene_name) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(limits = c("RNAi_elav", "RNAi_CyO", "TripCon_elav"),
                   labels = c("RNAi/elav-Gal4", "RNAi/+", "+/elav-Gal4")) +
  labs(x = "Genotype", 
       y = "Proportion of observations with pursuit") +
  scale_colour_manual(limits = c("RNAi_elav", "RNAi_CyO", "TripCon_elav"), values = c("black", "grey65", "grey65")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.30), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))

grid.arrange(arrangeGrob(pursuit.fig + theme(axis.title.x = element_blank()), 
                         pursuit.fig.ex + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
                         ncol = 2, 
                         bottom = textGrob("Genotype", gp = gpar(fontsize = 16)), 
                         widths = c(5, 2)))

grid.arrange(arrangeGrob(fc_val_plot + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                         fc_val_plot.ex + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank()),
                         ncol = 2, 
                         widths = c(5, 2)),
             arrangeGrob(pursuit.fig + theme(axis.title.x = element_blank(), strip.background = element_blank(),
                                             strip.text.x = element_blank()), 
                         pursuit.fig.ex + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background = element_blank(),
                                                strip.text.x = element_blank()),
                         ncol = 2, 
                         bottom = textGrob("Genotype", gp = gpar(fontsize = 16)), 
                         widths = c(5, 2)),
             ncol = 1,
             heights = c(1, 1.3))

pursuit.data.N <- pursuit.data %>%
  group_by(Gene_name) %>%
  summarise(n = length(Trial_id))

# Models

contrast1 <- c(-0.5, 1, -0.5) # compare trt to controls
contrast2 <- c(-1, 0, 1) # compare controls

contrasts <- rbind(contrast1, contrast2)

# compute generalized inverse to get contrast matrix
con.mat <- ginv(contrasts)
colnames(con.mat) <- c("Exp_vs_Controls", "Con_vs_con")

contrasts(pursuit.data$Treatment) <- con.mat

# add a trial-level variable
pursuit.data$Trial_id <- as.factor(1:length(pursuit.data$Gene))

pursuit.long <- pursuit.data %>%
  gather(Observation, pursuit, `1`:`12`)

pursuit.long$Observation <- as.numeric(pursuit.long$Observation)
pursuit.long$pursuit <- as.numeric(pursuit.long$pursuit)
pursuit.long <- pursuit.long[complete.cases(pursuit.long), ] # remove NAs (i.e. obs after a mating)


# Gene 1 - CG14153

pursuit.binom.g1.glm <- glmmTMB(pursuit ~ 1 + Mated + Treatment + Day + 
                                  Rack + Observation + (1 + Observation | Trial_id),
                                data = subset(pursuit.long, Gene_name == "CG14153"),
                                family = binomial())
summary(pursuit.binom.g1.glm)

plot(simulateResiduals(pursuit.binom.g1.glm))

# Gene 2 - verm

pursuit.binom.g2.glm <- glmmTMB(pursuit ~ 1 + Mated + Treatment + # no day effect as this was the first test - had to throw out first day 
                                  Observation + (1|Rack) + (1| Trial_id), 
                                data = subset(pursuit.long, Gene_name == "verm"),
                                family = binomial())
summary(pursuit.binom.g2.glm)

plot(simulateResiduals(pursuit.binom.g2.glm))

# Gene 5 - Drsl4

pursuit.binom.g5.glm <- glmmTMB(pursuit ~ 1 + Mated + Treatment + Day + 
                                  Rack + Observation + (1| Trial_id),
                                data = subset(pursuit.long, Gene_name == "Drsl4"),
                                family = binomial())
summary(pursuit.binom.g5.glm)

plot(simulateResiduals(pursuit.binom.g5.glm))

# Gene 6 - GstZ1

pursuit.binom.g6.glm <- glmmTMB(pursuit ~ 1 + Mated + Treatment + Day + 
                                  Rack + Observation + (1| Trial_id),
                                data = subset(pursuit.long, Gene_name == "GstZ1"),
                                family = binomial())
summary(pursuit.binom.g6.glm)

plot(simulateResiduals(pursuit.binom.g6.glm))

# Gene 7 - Nepl18

pursuit.binom.g7.glm <- glmmTMB(pursuit ~ 1 + Mated + Treatment + Day + 
                                  Rack + Observation + (1| Trial_id),
                                data = subset(pursuit.long, Gene_name == "Nepl18"),
                                family = binomial())
summary(pursuit.binom.g7.glm)

plot(simulateResiduals(pursuit.binom.g7.glm))

## Extra genes

# Gene 3 - Lsp2

pursuit.binom.g3.glm <- glmmTMB(pursuit ~ 1 + Mated + Treatment + Day +  ### Have to remove Rack, Obs to get model to fit
                                  (1| Trial_id),
                                data = subset(pursuit.long, Gene_name == "Lsp2"),
                                family = binomial())
summary(pursuit.binom.g3.glm) # Only con-con sig

plot(simulateResiduals(pursuit.binom.g3.glm))

# Gene 4 - Nazo

pursuit.binom.g4.glm <- glmmTMB(pursuit ~ 1 + Mated + Treatment + 
                                  (1| Trial_id),
                                data = subset(pursuit.long, Gene_name == "Nazo"),
                                family = binomial())
summary(pursuit.binom.g4.glm) # Both contrasts sig

plot(simulateResiduals(pursuit.binom.g4.glm)) 

