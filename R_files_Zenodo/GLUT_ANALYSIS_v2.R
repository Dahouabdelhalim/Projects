# GLUT Fed/Fasted and Tissue Analysis # Ver 2
# 20200308

# 1. Model (selection using AICc)
# 2. Normality assumptions of model residuals
# 3. Post-hoc tukey contrast/multiple comparisons
# 4. Plot

library('plyr')
library('dplyr')
library('MuMIn') # AIC, AICc
library('moments') # Skewness
library('emmeans')
library('lme4')
library('lmerTest') # Makes lme4 object 
library('ggplot2') # le sigh
library('LMERConvenienceFunctions')
library('scales') # displaying data as exponents
library('ggpubr')
library('extrafont')

# Data ######################################
setwd("C:\\\\Users\\\\Raafay\\\\Google Drive\\\\1 - Research\\\\Welch Lab\\\\Paper Data\\\\5_Paper_Data-PM,Cell")

# UTF8 encoding required for read.csv to avoid weird header letters
ffGLUTcell <- read.csv("ffGLUTdata_cell_livG1cell.csv", header = T, fileEncoding="UTF-8-BOM") # Cell fraction data
ffGLUTpm <- read.csv("ffGLUTdata_pm.csv", header = T, fileEncoding="UTF-8-BOM") # Plasma membrane data

# GLUT Intensity across all Tissue tested #
ffGLUT1cell <- filter(ffGLUTcell, ffGLUTcell$GLUT == 'GLUT1') # GLUT1 cell: mus, hrt
ffGLUT2cell <- filter(ffGLUTcell, ffGLUTcell$GLUT == 'GLUT2') # GLUT2 cell: mus, hrt, liv
ffGLUT3cell <- filter(ffGLUTcell, ffGLUTcell$GLUT == 'GLUT3') # GLUT3 cell: mus, hrt, liv
ffGLUT5cell <- filter(ffGLUTcell, ffGLUTcell$GLUT == 'GLUT5') # GLUT5 cell: mus, hrt, liv

ffGLUT1pm <- filter(ffGLUTpm, ffGLUTpm$GLUT == 'GLUT1') # GLUT1 PM: mus, hrt
ffGLUT2pm <- filter(ffGLUTpm, ffGLUTpm$GLUT == 'GLUT2') # GLUT2 PM: mus, hrt, liv
ffGLUT3pm <- filter(ffGLUTpm, ffGLUTpm$GLUT == 'GLUT3') # GLUT3 PM: mus, hrt, liv
ffGLUT5pm <- filter(ffGLUTpm, ffGLUTpm$GLUT == 'GLUT5') # GLUT5 PM: mus, hrt, liv

# AICc ######################################
# Function runs the following possible model combinations to determine AIC/AICc scores
mix.mods <- function(x, y){
  lm1 <- lmer(y ~ Condition + (1|Blot), data = x)
  lm2 <- lmer(y ~ Tissue + (1|Blot), data = x)
  lm3 <- lmer(y ~ Condition + Tissue + (1|Blot), data = x)
  lm4 <- lmer(y ~ Condition * Tissue + (1|Blot), data = x)
  lm5 <- lm(y ~ Condition * Tissue, data = x)
  aic.mods <- AIC(lm1, lm2, lm3, lm4, lm5) # runs AIC on all models with selected dataset
  aicc.mods <- AICc(lm1, lm2, lm3, lm4, lm5) %>% # runs AICc on all models with selected dataset
    mutate(AIC = aic.mods$AIC)
  return(aicc.mods)
}

G1pm.aicc <- mix.mods(ffGLUT1pm, ffGLUT1pm$Intensity) %>% print() # 4
G1cell.aicc <- mix.mods(ffGLUT1cell, ffGLUT1cell$Intensity) %>% print() # 4

G2pm.aicc <- mix.mods(ffGLUT2pm, ffGLUT2pm$Intensity) %>% print() # 4
G2cell.aicc <- mix.mods(ffGLUT2cell, ffGLUT2cell$Intensity) %>% print() # 4

G3pm.aicc <- mix.mods(ffGLUT3pm, ffGLUT3pm$Intensity) %>% print() # 4
G3cell.aicc <- mix.mods(ffGLUT3cell, ffGLUT3cell$Intensity) %>% print() # 4

G5pm.aicc <- mix.mods(ffGLUT5pm, ffGLUT5pm$Intensity) %>% print() # 4 
G5cell.aicc <- mix.mods(ffGLUT5cell, sqrt(ffGLUT5cell$Intensity)) %>% print() # 4


###################################### GLUT1 ######################################
#################################### WholeCell ####################################

######## Model ########
lm.G1cell  <- lmer(Intensity ~ Condition * Tissue + (1|Blot), data=ffGLUT1cell, REML = T) # no transformation
loglm.G1cell  <- lmer(log(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT1cell, REML = T) # log Intensity
sqrtlm.G1cell <- lmer(sqrt(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT1cell, REML =T) # sqrt Intensity

##### Assumptions #####
  #1 Shapiro.Wilkes test
shapiro.test(resid(lm.G1cell)) # p = 0.02183 (non-normal)
shapiro.test(resid(loglm.G1cell)) # p = 0.06889 (normal)
shapiro.test(resid(sqrtlm.G1cell)) # p = 0.1600 (normal)

  #2 Visual analysis of residual qqplot
par(mfrow = c(1, 3))
qqnorm(resid(lm.G1cell), main = "No Transformation")
qqline(resid(lm.G1cell)) 
qqnorm(resid(loglm.G1cell), main = "Log") # closest to line
qqline(resid(loglm.G1cell))
qqnorm(resid(sqrtlm.G1cell), main = "Sqrt") # pretty good too
qqline(resid(sqrtlm.G1cell))

  #3 fitted v. residuals
par(mfrow = c(1, 3))
plot(fitted(lm.G1cell), resid(lm.G1cell)) %>% abline(0,0) # most skewed
plot(fitted(loglm.G1cell), resid(loglm.G1cell)) %>% abline(0,0)
plot(fitted(sqrtlm.G1cell), resid(sqrtlm.G1cell)) %>% abline(0,0)
# all similar

  #4 histogram
par(mfrow = c(1, 3))
hist(resid(lm.G1cell))
hist(resid(loglm.G1cell)) # look similarly normal - left skewed
hist(resid(sqrtlm.G1cell)) # look similarly normal - left skewed

######## ANOVA ########
an.G1cell <- anova(sqrtlm.G1cell, ddf = "Ken") %>% print()
# Condition significant

###### Post Hoc #######
em.G1cell  <- emmeans(lm.G1cell, ~ Condition * Tissue, type = "response") # all estimated means

# em.G1cell  <- emmeans(loglm.G1cell, ~ Condition * Tissue, type = "response") # log lm creates ratios

mc.G1cell <- contrast(em.G1cell, "tukey", reverse = T) %>% 
  print()
tis.G1cell <- contrast(em.G1cell, "tukey", reverse = F, by = "Tissue") %>% 
  print()
con.G1cell <- contrast(em.G1cell, "tukey", reverse = F, by = "Condition") %>% 
  print()
# Fed/Fast flight muscle significant

pairs(em.G1cell)

######## Plot #########
max.G1cell <- max(data.frame(em.G1cell)$emmean+data.frame(em.G1cell)$SE)
min.G1cell <- min(data.frame(em.G1cell)$emmean-data.frame(em.G1cell)$SE)
gg.G1cell.bar <- ggplot(data.frame(em.G1cell), 
                        aes(x=interaction(Condition,Tissue),
                            y=emmean, 
                            fill=interaction(Condition,Tissue))) + 
  coord_cartesian(ylim = c(0, 1.16*max.G1cell), clip = "off") + 
  geom_bar(stat="identity", position="dodge", color='black') + 
  geom_errorbar(aes(ymin=pmax(emmean-SE, 0),
                    ymax=emmean+SE), width=.2, position='dodge') + 
  annotate("text", label = c("*"), x = c(1.5), y = 1.15*max.G1cell) + # fed/fast significance
  annotate("text", label = c("a", "b", "c", "c"), x = c(1, 2, 5, 6), y = 1.07*max.G1cell) + # tissue significance
  annotate("text", label = c("6", "6", "2", "2", "3", "3"), fontface = "italic", 
           x = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25), y = 0.5e07) + # sample sizes
  # annotate("text", x=c(2, 4), y=0.5e07, color = c("black"), family='sans', # ratio
  #          label = paste(round((1/ratio.G1cell$ratio), 2))) + #, round(ratio.G1cell$SE, 2), sep=" ± ")) +
  annotate("text", label = c("Flight\\nMuscle", "Heart", "Liver"), x = c(1.5, 3.5, 5.5), y = 0, vjust = c(2, 4, 4), hjust = 0.4) + # tissue names
  geom_segment(aes(x = 1, y = 1.12*max.G1cell, xend = 2, yend = 1.12*max.G1cell)) +
  # geom_segment(aes(x = 2.5, xend = 2.5, y = 0,  yend = 1.2*max.G1cell), linetype=2) + # line diving tissue
  scale_y_continuous("Protein Abundance\\n(arb. intensity units)", 
                     position = "left",
                     labels = scales::comma_format(accuracy=2, scale=0.0000001)) +
  scale_x_discrete("Whole Tissue Homogenate",
                   labels = c("Fed", "Fasted",
                              "Fed", "Fasted",
                              "Fed", "Fasted"),
                   limits = c("Fed.flightmuscle", "Fasted.flightmuscle",
                              "Fed.heart", "Fasted.heart", 
                              "Fed.liver", "Fasted.liver")) + 
  scale_fill_manual(values = c("lightpink", "indianred2",
                               "slategray1", "slategray3",
                               "cornsilk", "burlywood1")) +   
  theme(legend.position = "none", #plot.margin = unit(c(2, 1, 1, 1), "cm"), 
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(colour = 0),
        panel.grid.minor = element_line(colour = 0),
        # panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        # panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "gray90"), #grid
        plot.margin=unit(c(0.1, 0.1, 1, 0.1),"cm"),
        axis.title.x = element_text(vjust = 95))
print(gg.G1cell.bar)



###################################### GLUT1 ######################################
################################# Plasma Membrane #################################

######## Model ########
lm.G1pm  <- lmer(Intensity ~ Condition * Tissue + (1|Blot), data=ffGLUT1pm, REML = T) # no transformation
loglm.G1pm  <- lmer(log(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT1pm, REML = T) # log Intensity
sqrtlm.G1pm <- lmer(sqrt(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT1pm, REML =T) # sqrt Intensity

##### Assumptions #####
#1 Shapiro.Wilkes test
shapiro.test(resid(lm.G1pm)) # p = 0.4987 (normal)
shapiro.test(resid(loglm.G1pm)) # p = 0.6085 (normal)
shapiro.test(resid(sqrtlm.G1pm)) # p = 0.2749 (normal)

#2 Visual analysis of residual qqplot
par(mfrow = c(1, 3))
qqnorm(resid(lm.G1pm), main = "No Transformation") # 2nd closest
qqline(resid(lm.G1pm)) 
qqnorm(resid(loglm.G1pm), main = "Log") # closest to line
qqline(resid(loglm.G1pm))
qqnorm(resid(sqrtlm.G1pm), main = "Sqrt") # furthest from line
qqline(resid(sqrtlm.G1pm))

#3 fitted v. residuals
par(mfrow = c(1, 3))
plot(fitted(lm.G1pm), resid(lm.G1pm)) %>% abline(0,0) # most skewed
plot(fitted(loglm.G1pm), resid(loglm.G1pm)) %>% abline(0,0)
plot(fitted(sqrtlm.G1pm), resid(sqrtlm.G1pm)) %>% abline(0,0)
# all similar

#4 histogram
par(mfrow = c(1, 3))
hist(resid(lm.G1pm))# most normal looking; but hard to tell
hist(resid(loglm.G1pm)) 
hist(resid(sqrtlm.G1pm)) 

######## ANOVA ########
an.G1pm <- anova(loglm.G1pm, ddf = "Ken") %>% print()
# Condition almost significant
# Tissue significant
# Interaction significant

###### Post Hoc #######
em.G1pm  <- emmeans(lm.G1pm, ~ Condition * Tissue, type = "response") # all estimated means

#em.G1pm  <- emmeans(loglm.G1pm, ~ Condition * Tissue, type = "response") # log lm for ratios
mc.G1pm <- contrast(em.G1pm, "tukey", reverse = F) %>% 
  print()
tis.G1pm <- contrast(em.G1pm, "tukey", reverse = F, by = "Tissue") %>% 
  print()
con.G1pm <- contrast(em.G1pm, "tukey", reverse = F, by = "Condition") %>% 
  print()
  
# Fed/Fast flight muscle significant


####### Plot ##########
# GLUT 1 # 
# Plasma Membrane # 
max.G1pm <- max(data.frame(em.G1pm)$emmean+data.frame(em.G1pm)$SE)
gg.G1pm.bar <- ggplot(data.frame(em.G1pm), 
                      aes(x=interaction(Condition,Tissue), 
                          y=emmean, 
                          fill=interaction(Condition, Tissue))) + 
  geom_bar(stat="identity", position="dodge", color='black') + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2, position='dodge') + 
  coord_cartesian(ylim = c(0, 1.16*max.G1pm), clip = "off") +
  annotate("text", label = c("a", "b"), x = c(1, 3), y = 1.09*max.G1pm) + # significance
  annotate("text", label = c("*"), x = c(1.5), y = 1.15*max.G1pm) + # significance
  annotate("text", label = c("6", "6", "4", "5", "3", "3"), fontface = "italic", 
           x = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25), y = 0.4e07) + # sample sizes
  # annotate("text", label = paste(round((1/ratio.G1pm$ratio), 2)),
  #          x=c(2, 4), y=0.4e07, color = c("black")) +  # ratio
  #, round(ratio.G1pm$SE, 2), sep=" ± ")) +
  annotate("text", label = c("Flight\\nMuscle", "Heart", "Liver"), 
           x = c(1.5, 3.5, 5.5), y = 0, vjust = c(2, 4, 4), hjust = 0.4) + # tissue names
  annotate("text", label = c("N/A", "N/A"), x = c(5, 6), y = 0, fontface = "italic", size = 3) + # n/a's
  geom_segment(aes(x = 1, y = 1.13*max.G1pm, xend = 2, yend = 1.13*max.G1pm)) +
#  geom_segment(aes(x = 3, y = 1.13*max.G1pm, xend = 4, yend = 1.13*max.G1pm)) +
#  geom_segment(aes(x = 2.5, xend = 2.5, y = 0,  yend = 1.2*max.G1pm), linetype=2) +
#  labs(title = "Plasma Membrane\\n") +
  scale_y_continuous("Protein Abundance\\n(arb. intensity units)", 
                     position = "right",
                     labels = scales::comma_format(accuracy=2, scale=0.0000001)) +
  scale_x_discrete(name = "Plasma Membrane Fraction",
                   labels = c("Fed", "Fasted",
                              "Fed", "Fasted", 
                              "Fed", "Fasted"),
                   limits = c("Fed.flightmuscle", "Fasted.flightmuscle",
                              "Fed.heart", "Fasted.heart", 
                              "Fed.liver", "Fast.liver")) + 
  scale_fill_manual(values = c("lightpink", "indianred2",
                               "slategray1", "slategray3",
                               "cornsilk", "burlywood1")) +   
  theme(legend.position = "none", #plot.margin = unit(c(2, 1, 1, 1), "cm"), 
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(colour = 0),
        panel.grid.minor = element_line(colour = 0),
        plot.margin=unit(c(0.1, 0.1, 1, 0.1),"cm"),
        # panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        # panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "gray90"), # grid
        
        axis.title.x = element_text(vjust = 95))
print(gg.G1pm.bar)

###################################### GLUT2 ######################################
#################################### WholeCell ####################################

######## Model ########
lm.G2cell  <- lmer(Intensity ~ Condition * Tissue + (1|Blot), data=ffGLUT2cell, REML = T) # no transformation
noreml.lm.G2cell  <- lmer(Intensity ~ Condition * Tissue + (1|Blot), data=ffGLUT2cell, REML = F) # no transformation
# turning REML off b/c it doesn't work with contrast here
# turning it back on for Anova with Ken method
loglm.G2cell  <- lmer(log(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT2cell, REML = T) # log Intensity
sqrtlm.G2cell <- lmer(sqrt(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT2cell, REML = T) # sqrt Intensity


##### Assumptions #####
#1 Shapiro.Wilkes test
shapiro.test(resid(lm.G2cell)) # p = 0.6558 (normal)
shapiro.test(resid(loglm.G2cell)) # p = 0.7434 (normal)
shapiro.test(resid(sqrtlm.G2cell)) # p = 0.8388 (normal)

#2 Visual analysis of residual qqplot
par(mfrow = c(1, 3))
qqnorm(resid(lm.G2cell), main = "No Transformation")
qqline(resid(lm.G2cell)) 
qqnorm(resid(loglm.G2cell), main = "Log") # closest to line
qqline(resid(loglm.G2cell))
qqnorm(resid(sqrtlm.G2cell), main = "Sqrt")
qqline(resid(sqrtlm.G2cell))
# all v. similar

#3 fitted v. residuals
par(mfrow = c(1, 3))
plot(fitted(lm.G2cell), resid(lm.G2cell)) %>% abline(0,0) 
plot(fitted(loglm.G2cell), resid(loglm.G2cell)) %>% abline(0,0)
plot(fitted(sqrtlm.G2cell), resid(sqrtlm.G2cell)) %>% abline(0,0)
# all similar

#4 histogram
par(mfrow = c(1, 3))
hist(resid(lm.G2cell))
hist(resid(loglm.G2cell)) 
hist(resid(sqrtlm.G2cell)) # most normal looking

######## ANOVA ########
an.G2cell <- anova(lm.G2cell, ddf = "Ken") %>% print()
# Condition significant

###### Post Hoc #######
em.G2cell  <- emmeans(lm.G2cell, ~ Condition * Tissue, type = "response") # all estimated means

#em.G2cell  <- emmeans(loglm.G2cell, ~ Condition * Tissue, type = "response") # log lm for ratios
mc.G2cell <- contrast(em.G2cell, "tukey", reverse = T) %>% 
  print()
tis.G2cell <- contrast(em.G2cell, "tukey", reverse = F, by = "Tissue") %>% 
  print()
con.G2cell <- contrast(em.G2cell, "tukey", reverse = F, by = "Condition") %>% 
  print()
# No significances

######## Plot #########
max.G2cell <- max(data.frame(em.G2cell)$emmean+data.frame(em.G2cell)$SE)
gg.G2cell.bar <- ggplot(data.frame(em.G2cell), 
                        aes(x=interaction(Condition,Tissue), 
                            y=emmean, 
                            fill=interaction(Condition, Tissue))) + 
  coord_cartesian(ylim = c(0, 1.21*max.G2cell), clip = "off") + 
  geom_bar(stat="identity", position="dodge", color='black') + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2, position='dodge') + 
  #  labs(title = "GLUT2", subtitle = "Whole Cell") +
  annotate("text", label = c("*"), x = c(1.5), y = 1.15*max.G2cell) +
  # annotate("text", label = paste(round((1/ratio.G2cell$ratio), 2)), #, round(ratio.G2cell$SE, 2), sep="\\n±"), # ratio
  #          x=c(2, 4, 6), y=0.9e+07, color = c("black"), family="sans", 
  #          hjust="middle", lineheight=0.85) +
  annotate("text", label = c("Flight\\nMuscle", "Heart", "Liver"), 
           x = c(1.5, 3.5, 5.5), y = 0, vjust = c(2, 4, 4), hjust = 0.4) + # tissue names
  annotate("text", label = c("4", "4", "2", "2", "3", "3"), fontface = "italic", 
           x = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25), y = 0.9e07) +
  geom_segment(aes(x = 1, xend = 2, y = 1.1*max.G2cell, yend = 1.1*max.G2cell)) +
#  geom_segment(aes(x = 2.5, xend = 2.5, y = 0,  yend = 1.25*max.G2cell), linetype=2) +
#  geom_segment(aes(x = 4.5, xend = 4.5, y = 0,  yend = 1.25*max.G2cell), linetype=2) +
  scale_y_continuous("Protein Abundance\\n(arb. intensity units)",
                     position = "left", 
                     labels = scales::comma_format(accuracy=2, scale=0.0000001)) +
  scale_x_discrete(name = "Whole Tissue Homogenate",
                   labels = c("Fed", "Fasted",
                              "Fed", "Fasted", 
                              "Fed", "Fasted"),
                   limits = c("Fed.flightmuscle", "Fasted.flightmuscle", 
                              "Fed.heart", "Fasted.heart", 
                              "Fed.liver", "Fasted.liver")) + 
  scale_fill_manual(values = c("lightpink", "indianred2",
                               "slategray1", "slategray3",
                               "cornsilk", "burlywood1")) +   
  theme(legend.position = "none", #plot.margin = unit(c(2, 1, 1, 1), "cm"), 
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(colour = 0),
        panel.grid.minor = element_line(colour = 0),
        plot.margin=unit(c(0.1, 0.1, 1, 0.1),"cm"),
        # panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        # panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        axis.title.x = element_text(vjust = 93))
print(gg.G2cell.bar)


###################################### GLUT2 ######################################
################################# Plasma Membrane #################################

######## Model ########
lm.G2pm  <- lmer(Intensity ~ Condition * Tissue + (1|Blot), data=ffGLUT2pm, REML = T) # no transformation
loglm.G2pm  <- lmer(log(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT2pm, REML = T) # log Intensity
sqrtlm.G2pm <- lmer(sqrt(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT2pm, REML =T) # sqrt Intensity

##### Assumptions #####
#1 Shapiro.Wilkes test
shapiro.test(resid(lm.G2pm)) # p = 0.4442 (normal)
shapiro.test(resid(loglm.G2pm)) # p = 0.5870 (normal)
shapiro.test(resid(sqrtlm.G2pm)) # p = 0.8418 (normal)

#2 Visual analysis of residual qqplot
par(mfrow = c(1, 3))
qqnorm(resid(lm.G2pm), main = "No Transformation") # closest to line
qqline(resid(lm.G2pm)) 
qqnorm(resid(loglm.G2pm), main = "Log") 
qqline(resid(loglm.G2pm))
qqnorm(resid(sqrtlm.G2pm), main = "Sqrt") # furthest from line
qqline(resid(sqrtlm.G2pm))

#3 fitted v. residuals
par(mfrow = c(1, 3))
plot(fitted(lm.G2pm), resid(lm.G2pm)) %>% abline(0,0) # most even above/below
plot(fitted(loglm.G2pm), resid(loglm.G2pm)) %>% abline(0,0)
plot(fitted(sqrtlm.G2pm), resid(sqrtlm.G2pm)) %>% abline(0,0) # similarly even
# all similar

#4 histogram
par(mfrow = c(1, 3))
hist(resid(lm.G2pm))# normal looking 
hist(resid(loglm.G2pm)) 
hist(resid(sqrtlm.G2pm)) # also normal looking

######## ANOVA ########
an.G2pm <- anova(lm.G2pm, ddf = "Ken") %>% print()
# Nothing sig.

###### Post Hoc #######
em.G2pm  <- emmeans(lm.G2pm, ~ Condition * Tissue, type = "response") # all estimated means

#em.G2pm  <- emmeans(loglm.G2pm, ~ Condition * Tissue, type = "response") # log lm for ratios
mc.G2pm <- contrast(em.G2pm, "tukey", reverse = F) %>% 
  print()
tis.G2pm <- contrast(em.G2pm, "tukey", reverse = F, by = "Tissue") %>% 
  print()
con.G2pm <- contrast(em.G2pm, "tukey", reverse = F, by = "Condition") %>% 
  print()
# No significances

####### Plot ##########
max.G2pm <- max(data.frame(em.G2pm)$emmean+data.frame(em.G2pm)$SE)
gg.G2pm.bar <- ggplot(data.frame(em.G2pm), 
                      aes(x=interaction(Condition,Tissue), 
                          y=emmean, 
                          fill=interaction(Condition,Tissue)) )+ 
  coord_cartesian(ylim = c(0, 1.11*max.G2pm), clip = "off") + 
  geom_bar(stat="identity", position="dodge", color='black') + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2, position='dodge') + 
  annotate("text", label = c(""), x = c(1.5), y = 1.15*max.G2pm) +
  # annotate("text", label = paste(round((1/ratio.G2pm$ratio), 2)), # round(ratio.G2pm$SE, 2), sep="\\n±"),   # ratio
  #          x=c(2, 4, 6), y=0.8e+07, color = c("black"), family="sans",
  #          hjust="middle", lineheight=0.85) +
  annotate("text", label = c("Flight\\nMuscle", "Heart", "Liver"), 
           x = c(1.5, 3.5, 5.5), y = 0, vjust = c(2, 4, 4), hjust = 0.4) + # tissue names
  annotate("text", label = c("4", "4", "5", "5", "3", "3"), fontface = "italic", 
           x = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25), y = 0.7e07) + # sample sizes
  # geom_segment(aes(x = 2.5, xend = 2.5, y = 0,  yend = 1.15*max.G2pm), linetype=2) +
  # geom_segment(aes(x = 4.5, xend = 4.5, y = 0,  yend = 1.15*max.G2pm), linetype=2) +
  #  labs(title = "GLUT2", subtitle = "Plasma Membrane") +
  scale_y_continuous("Protein Abundance\\n(arb. intensity units)",
                     position = "right", 
                     labels = scales::comma_format(accuracy=2, scale=0.0000001)) +
  scale_x_discrete(name = "Plasma Membrane Fraction",
                   labels = c("Fed", "Fasted",
                              "Fed", "Fasted", 
                              "Fed", "Fasted"),
                   limits = c("Fed.flightmuscle", "Fasted.flightmuscle", 
                              "Fed.heart", "Fasted.heart", 
                              "Fed.liver", "Fasted.liver")) + 
  scale_fill_manual(values = c("lightpink", "indianred2",
                               "slategray1", "slategray3",
                               "cornsilk", "burlywood1")) + 
  theme(legend.position = "none", #plot.margin = unit(c(2, 1, 1, 1), "cm"), 
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(colour = 0),
        panel.grid.minor = element_line(colour = 0),
        plot.margin=unit(c(0.1, 0.1, 1, 0.1),"cm"),
        # panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        # panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        axis.title.x = element_text(vjust = 93))
print(gg.G2pm.bar)


###################################### GLUT3 ######################################
#################################### WholeCell ####################################

######## Model ########
lm.G3cell  <- lmer(Intensity ~ Condition * Tissue + (1|Blot), data=ffGLUT3cell, REML = T) # no transformation
noremllm.G3cell  <- lmer(Intensity ~ Condition * Tissue + (1|Blot), data=ffGLUT3cell, REML = F) # no transformation

loglm.G3cell  <- lmer(log(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT3cell, REML = T) # log Intensity
sqrtlm.G3cell <- lmer(sqrt(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT3cell, REML =T) # sqrt Intensity

##### Assumptions #####
#1 Shapiro.Wilkes test
shapiro.test(resid(lm.G3cell)) # p = 0.9073 (normal)
shapiro.test(resid(loglm.G3cell)) # p = 0.9627 (normal)
shapiro.test(resid(sqrtlm.G3cell)) # p = 0.7772 (normal)

#2 Visual analysis of residual qqplot
par(mfrow = c(1, 3))
qqnorm(resid(lm.G3cell), main = "No Transformation") # closest to line
qqline(resid(lm.G3cell)) 
qqnorm(resid(loglm.G3cell), main = "Log") 
qqline(resid(loglm.G3cell))
qqnorm(resid(sqrtlm.G3cell), main = "Sqrt")
qqline(resid(sqrtlm.G3cell))

#3 fitted v. residuals
par(mfrow = c(1, 3))
plot(fitted(lm.G3cell), resid(lm.G3cell)) %>% abline(0,0) 
plot(fitted(loglm.G3cell), resid(loglm.G3cell)) %>% abline(0,0)
plot(fitted(sqrtlm.G3cell), resid(sqrtlm.G3cell)) %>% abline(0,0)# most skewed
# all similar

#4 histogram
par(mfrow = c(1, 3))
hist(resid(lm.G3cell))# most normal looking
hist(resid(loglm.G3cell)) 
hist(resid(sqrtlm.G3cell)) 

######## ANOVA ########
an.G3cell <- anova(lm.G3cell, ddf = "Ken") %>% print()
an.G3cell <- anova(noremllm.G3cell, ddf = "Sat") %>% print()
# Condition significant
# Interaction almost significant

###### Post Hoc #######
em.G3cell  <- emmeans(lm.G3cell, ~ Condition * Tissue, type = "response") # all estimated means

#em.G3cell  <- emmeans(loglm.G3cell, ~ Condition * Tissue, type = "response") # log lm for ratios
mc.G3cell <- contrast(em.G3cell, "tukey", reverse = T) %>% 
  print()
tis.G3cell <- contrast(em.G3cell, "tukey", reverse = F, by = "Tissue") %>% 
  print()
con.G3cell <- contrast(em.G3cell, "tukey", reverse = F, by = "Condition") %>% 
  print()
# Fed/Fast liver significant


######### Plot ############
max.G3cell <- max(data.frame(em.G3cell)$emmean+data.frame(em.G3cell)$SE)
gg.G3cell.bar <- ggplot(data.frame(em.G3cell), 
                        aes(x=interaction(Condition,Tissue), 
                            y=emmean, fill=interaction(Condition, Tissue))) + 
  coord_cartesian(ylim = c(0, 1.2*max.G3cell), clip = "off") + 
  geom_bar(stat="identity", position="dodge", color='black') + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2, position='dodge') + 
  #  labs(title = "GLUT3", subtitle = "Whole Cell") +
  annotate("text", label = c("b", "a"), x = c(3, 5), y = 1.10*max.G3cell) +
  annotate("text", label = c("*"), x = c(1.5, 5.5), y = 1.2*max.G3cell) +
  # annotate("text", label = paste(round((1/ratio.G3cell$ratio), 2)), #, round(ratio.G3cell$SE, 2), sep="\\n±"), # ratio
  #          x=c(2, 4, 6), y=0.5e+07, color = c("black"),
  #          hjust="middle", lineheight=0.85) +
  #  annotate("text", label = c(1, 1, 1), x=c(1, 3, 5), y=1e+03, color = c("black")) + 
  annotate("text", label = c("Flight\\nMuscle", "Heart", "Liver"), 
           x = c(1.5, 3.5, 5.5), y = 0, vjust = c(2, 4, 4), hjust = 0.4) + # tissue names
  annotate("text", label = c("4", "4", "2", "2", "3", "3"), fontface = "italic", 
           x = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25), y = 0.5e07) + # sample sizes
  # geom_segment(aes(x = 1.1, xend = 1.9, y = 1.1*max.G3cell, yend = 1.1*max.G3cell)) +
  geom_segment(aes(x = 5, xend = 6, y = 1.17*max.G3cell, yend = 1.17*max.G3cell)) +
  geom_segment(aes(x = 1, xend = 2, y = 1.17*max.G3cell, yend = 1.17*max.G3cell)) +
  # geom_segment(aes(x = 2.5, xend = 2.5, y = 0,  yend = 1.25*max.G3cell), linetype=2) +
  # geom_segment(aes(x = 4.5, xend = 4.5, y = 0,  yend = 1.25*max.G3cell), linetype=2) +
  scale_y_continuous("Protein Abundance\\n(arb. intensity units)",
                     position = "left", 
                     labels = scales::comma_format(accuracy=2, scale=0.0000001)) +
  scale_x_discrete(name = "Whole Tissue Homogenate",
                   labels = c("Fed", "Fasted",
                              "Fed", "Fasted", 
                              "Fed", "Fasted"),
                   limits = c("Fed.flightmuscle", "Fasted.flightmuscle", 
                              "Fed.heart", "Fasted.heart", 
                              "Fed.liver", "Fasted.liver")) + 
  scale_fill_manual(values = c("lightpink", "indianred2",
                               "slategray1", "slategray3",
                               "cornsilk", "burlywood1")) +   
  theme(legend.position = "none", #plot.margin = unit(c(2, 1, 1, 1), "cm"), 
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(colour = 0),
        panel.grid.minor = element_line(colour = 0),
        plot.margin=unit(c(0.1, 0.1, 1, 0.1),"cm"),
        # panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        # panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        axis.title.x = element_text(vjust = 93))
print(gg.G3cell.bar)


###################################### GLUT3 ######################################
################################# Plasma Membrane #################################

######## Model ########
lm.G3pm  <- lmer(Intensity ~ Condition * Tissue + (1|Blot), data=ffGLUT3pm, REML = T) # no transformation
loglm.G3pm  <- lmer(log(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT3pm, REML = T) # log Intensity
sqrtlm.G3pm <- lmer(sqrt(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT3pm, REML =T) # sqrt Intensity

##### Assumptions #####
#1 Shapiro.Wilkes test
shapiro.test(resid(lm.G3pm)) # p = 0.8613 (normal)
shapiro.test(resid(loglm.G3pm)) # p = 0.2845 (normal)
shapiro.test(resid(sqrtlm.G3pm)) # p = 0.2398 (normal)

#2 Visual analysis of residual qqplot
par(mfrow = c(1, 3))
qqnorm(resid(lm.G3pm), main = "No Transformation") # closest to line
qqline(resid(lm.G3pm)) 
qqnorm(resid(loglm.G3pm), main = "Log") 
qqline(resid(loglm.G3pm))
qqnorm(resid(sqrtlm.G3pm), main = "Sqrt")
qqline(resid(sqrtlm.G3pm))
# all similar

#3 fitted v. residuals
par(mfrow = c(1, 3))
plot(fitted(lm.G3pm), resid(lm.G3pm)) %>% abline(0,0) 
plot(fitted(loglm.G3pm), resid(loglm.G3pm)) %>% abline(0,0)
plot(fitted(sqrtlm.G3pm), resid(sqrtlm.G3pm)) %>% abline(0,0)
# all similar

#4 histogram
par(mfrow = c(1, 3))
hist(resid(lm.G3pm))# most normal looking
hist(resid(loglm.G3pm)) 
hist(resid(sqrtlm.G3pm)) 

######## ANOVA ########
an.G3pm <- anova(lm.G3pm, ddf = "Ken") %>% print()
# Condition  significant
# Interaction significant

###### Post Hoc #######
em.G3pm  <- emmeans(lm.G3pm, ~ Condition * Tissue, type = "response") # all estimated means

#em.G3pm  <- emmeans(loglm.G3pm, ~ Condition * Tissue, type = "response") # log lm for ratios
mc.G3pm <- contrast(em.G3pm, "tukey", reverse = F) %>% 
  print()
tis.G3pm <- contrast(em.G3pm, "tukey", reverse = F, by = "Tissue") %>% 
  print()
con.G3pm <- contrast(em.G3pm, "tukey", reverse = F, by = "Condition") %>% 
  print()# Fed/Fast liver significant

######## Plots ########
# GLUT 3 #
# Plasma Membrane # 
max.G3pm <- max(data.frame(em.G3pm)$emmean+data.frame(em.G3pm)$SE)
gg.G3pm.bar <- ggplot(data.frame(em.G3pm), 
                      aes(x=interaction(Condition,Tissue), 
                          y=emmean, 
                          fill=interaction(Condition,Tissue))) + 
  coord_cartesian(ylim = c(0, 1.2*max.G3pm), clip = "off") + 
  geom_bar(stat="identity", position="dodge", color='black') + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2, position='dodge') + 
  annotate("text", label = c("*"), x = c(5.5), y = 1.15*max.G3pm) +
  # annotate("text", label = paste(round((1/ratio.G3pm$ratio), 2)), # round(ratio.G3pm$SE, 2), sep="\\n±"), # ratio
  #          x=c(2, 4, 6), y=0.5e+07, color = c("black"),
  #          hjust="middle", lineheight=0.85) +
  annotate("text", label = c("Flight\\nMuscle", "Heart", "Liver"), 
           x = c(1.5, 3.5, 5.5), y = 0, vjust = c(2, 4, 4), hjust = 0.4) + # tissue names
  annotate("text", label = c("4", "4", "5", "5", "3", "3"), fontface = "italic", 
           x = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25), y = 0.5e07) + # sample sizes
  geom_segment(aes(x = 5,  xend = 6, y = 1.1*max.G3pm, yend = 1.1*max.G3pm)) +
  # geom_segment(aes(x = 2.5, xend = 2.5, y = 0,  yend = 1.25*max.G3pm), linetype=2) +
  # geom_segment(aes(x = 4.5, xend = 4.5, y = 0,  yend = 1.25*max.G3pm), linetype=2) +
  #  labs(title = "GLUT2", subtitle = "Plasma Membrane") +
  scale_y_continuous("Protein Abundance\\n(arb. intensity units)",
                     position = "right", 
                     labels = scales::comma_format(accuracy=2, scale=0.0000001)) +
  scale_x_discrete(name = "Plasma Membrane Fraction",
                   labels = c("Fed", "Fasted",
                              "Fed", "Fasted", 
                              "Fed", "Fasted"),
                   limits = c("Fed.flightmuscle", "Fasted.flightmuscle", 
                              "Fed.heart", "Fasted.heart", 
                              "Fed.liver", "Fasted.liver")) + 
  scale_fill_manual(values = c("lightpink", "indianred2",
                               "slategray1", "slategray3",
                               "cornsilk", "burlywood1")) + 
  theme(legend.position = "none", #plot.margin = unit(c(2, 1, 1, 1), "cm"), 
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(colour = 0),
        panel.grid.minor = element_line(colour = 0),
        plot.margin=unit(c(0.1, 0.1, 1, 0.1),"cm"),
        # panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        # panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        axis.title.x = element_text(vjust = 93))
print(gg.G3pm.bar)


###################################### GLUT5 ######################################
#################################### WholeCell ####################################

######## Model ########
lm.G5cell  <- lmer(Intensity ~ Condition * Tissue + (1|Blot), data=ffGLUT5cell, REML = F) # no transformation
loglm.G5cell  <- lmer(log(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT5cell, REML = T) # log Intensity
sqrtlm.G5cell <- lmer(sqrt(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT5cell, REML =T) # sqrt Intensity

##### Assumptions #####
#1 Shapiro.Wilkes test
shapiro.test(resid(lm.G5cell)) # p = 1.98e-05 (non-normal)
shapiro.test(resid(loglm.G5cell)) # p = 0.0001718 (non-normal)
shapiro.test(resid(sqrtlm.G5cell)) # p = 5.504e-05 (non-normal)


#2 Visual analysis of residual qqplot
par(mfrow = c(1, 3))
qqnorm(resid(lm.G5cell), main = "No Transformation") # closest to line
qqline(resid(lm.G5cell)) 
qqnorm(resid(loglm.G5cell), main = "Log")
qqline(resid(loglm.G5cell))
qqnorm(resid(sqrtlm.G5cell), main = "Sqrt")
qqline(resid(sqrtlm.G5cell))

# outliers present

#3 fitted v. residuals
par(mfrow = c(1, 3))
plot(fitted(lm.G5cell), resid(lm.G5cell)) %>% abline(0,0) 
plot(fitted(loglm.G5cell), resid(loglm.G5cell)) %>% abline(0,0)
plot(fitted(sqrtlm.G5cell), resid(sqrtlm.G5cell)) %>% abline(0,0)
# outliers present; hard to judge (log/sqrt best)

#4 histogram
par(mfrow = c(1, 3))
hist(resid(lm.G5cell))
hist(resid(loglm.G5cell)) # most normal looking
hist(resid(sqrtlm.G5cell)) 
# outliers present; hard to tell (norm/sqrt best)



######## ANOVA ########
an.G5cell <- anova(lm.G5cell, ddf = "Ken") %>% print()
# no significances

###### Post Hoc #######
em.G5cell  <- emmeans(loglm.G5cell, ~ Condition * Tissue, type = "response") # all estimated means (log lm for ratios)
mc.G5cell <- contrast(em.G5cell, "tukey", reverse = F) %>% 
  print()
tis.G5cell <- contrast(em.G5cell, "tukey", reverse = F, by = "Tissue") %>% 
  print()
con.G5cell <- contrast(em.G5cell, "tukey", reverse = F, by = "Condition") %>% 
  print()
# no significances

# no difference in final outcome (transformations are v. similar). Choosing log transform as it gives best error bars


# Removing Outliers - significances don't look right with 'outliers' removed
# noout.G5cell <- romr.fnc(lm.G5cell, data=ffGLUT5cell, trim=2.5) # above 2.5 stdev removed; flightmuscle H12 Fasted
# 
# noout.lmG5cell <- lmer((Intensity) ~ Condition * Tissue + (1|Blot), data=noout.G5cell$data, REML = T) 
# noout.logG5cell <- lmer(log(Intensity) ~ Condition * Tissue + (1|Blot), data=noout.G5cell$data, REML = F)
# noout.sqrtG5cell <- lmer(sqrt(Intensity) ~ Condition * Tissue + (1|Blot), data=noout.G5cell$data, REML = T)
# 
# shapiro.test(resid(noout.lmG5cell)) # p-value = 0.03778 (non-normal)
# shapiro.test(resid(noout.logG5cell)) # p-value = 0.02584 (non-normal)
# shapiro.test(resid(noout.sqrtG5cell)) # p-value = 0.008846 (non-normal)
# par(mfrow = c(1, 3))
# qqnorm(resid(noout.lmG5cell), main = "no Transformation")
# qqline(resid(noout.lmG5cell))
# qqnorm(resid(noout.logG5cell), main = "log") # closest to line
# qqline(resid(noout.logG5cell))
# qqnorm(resid(noout.sqrtG5cell), main = "sqrt")
# qqline(resid(noout.sqrtG5cell))
# par(mfrow = c(1, 3))
# plot(fitted(noout.lmG5cell), resid(noout.lm.G5cell)) %>% abline(0,0) # very skewed
# plot(fitted(noout.logG5cell), resid(noout.logG5cell)) %>% abline(0,0) # least skewed
# plot(fitted(noout.sqrtG5cell), resid(noout.sqrtG5cell)) %>% abline(0,0)
# par(mfrow = c(1, 3))
# hist(resid(noout.lmG5cell))# most normal looking
# hist(resid(noout.logG5cell)) 
# hist(resid(noout.sqrtG5cell)) 
# 
# noout.anG5cell <- anova(noout.logG5cell, ddf = "Ken") %>% print()
# 
# noout.emG5cell <- emmeans(noout.logG5cell, ~ Condition * Tissue, type = "response") # all estimated means
# noout.mcG5cell <- contrast(noout.emG5cell, "tukey", reverse = T) %>% print()
# noout.tisG5cell <- contrast(noout.emG5cell, "tukey", reverse = T, by = "Tissue") %>% print()
# noout.conG5cell <- contrast(noout.emG5cell, "tukey", reverse = F, by = "Condition") %>% print()


######## Plot #########
max.G5cell <- max(data.frame(em.G5cell)$response+data.frame(em.G5cell)$SE)
gg.G5cell.bar <- ggplot(data.frame(em.G5cell), 
                        aes(x=interaction(Condition,Tissue), 
                            y=response, 
                            fill=interaction(Condition,Tissue))) + 
  coord_cartesian(ylim = c(0, 1.0*max.G5cell), clip = "off") + 
  geom_bar(stat="identity", position="dodge", color='black') + 
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=.2, position='dodge') + 
  #annotate("text", label = c("*"), x = c(1.5), y = 0.92*max.G5cell) + # significance
  # annotate("text", 
  #          label = paste(round((1/ratio.G5cell$ratio), 2)), # round(ratio.G5cell$SE, 2), sep="\\n±"), # ratio
  #          x=c(2, 4, 6), y=0.055e+08, color = c("black"),
  #          hjust="middle", lineheight=0.85) +
  annotate("text", label = c("Flight\\nMuscle", "Heart", "Liver"), 
           x = c(1.5, 3.5, 5.5), y = 0, vjust = c(2, 4, 4), hjust = 0.4) + # tissue names
  annotate("text", label = c("4", "4", "2", "2", "3", "3"), fontface = "italic", 
           x = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25), y = 0.055e+08) + # sample sizes
  #geom_segment(aes(x = 1.2, y = 0.9*max.G5cell, xend = 1.8, yend = 0.9*max.G5cell)) +
  # geom_segment(aes(x = 2.5, xend = 2.5, y = 0,  yend = 1.05*max.G5cell), linetype=2) +
  # geom_segment(aes(x = 4.5, xend = 4.5, y = 0,  yend = 1.05*max.G5cell), linetype=2) +
  scale_y_continuous("Protein Abundance\\n(arb. intensity units)",
                     position = "left", 
                     labels = scales::comma_format(accuracy=2, scale=0.0000001)) +
  scale_x_discrete(name = "Whole Tissue Homogenate",
                   labels = c("Fed", "Fasted",
                              "Fed", "Fasted", 
                              "Fed", "Fasted"),
                   limits = c("Fed.flightmuscle", "Fasted.flightmuscle", 
                              "Fed.heart", "Fasted.heart", 
                              "Fed.liver", "Fasted.liver")) + 
  scale_fill_manual(values = c("lightpink", "indianred2",
                               "slategray1", "slategray3",
                               "cornsilk", "burlywood1")) +   
  theme(legend.position = "none", #plot.margin = unit(c(2, 1, 1, 1), "cm"), 
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(colour = 0),
        panel.grid.minor = element_line(colour = 0),
        plot.margin=unit(c(0.1, 0.1, 1, 0.1),"cm"),
        # panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        # panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        axis.title.x = element_text(vjust = 93))
print(gg.G5cell.bar)


###################################### GLUT5 ######################################
################################# Plasma Membrane #################################

######## Model ########
lm.G5pm  <- lmer(Intensity ~ Condition * Tissue + (1|Blot), data=ffGLUT5pm, REML = T) # no transformation
loglm.G5pm  <- lmer(log(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT5pm, REML = T) # log Intensity
sqrtlm.G5pm <- lmer(sqrt(Intensity) ~ Condition * Tissue + (1|Blot), data=ffGLUT5pm, REML =T) # sqrt Intensity

##### Assumptions #####
#1 Shapiro.Wilkes test
shapiro.test(resid(lm.G5pm)) # p = 0.5957 (normal)
shapiro.test(resid(loglm.G5pm)) # p = 0.4303 (normal)
shapiro.test(resid(sqrtlm.G5pm)) # p = 0.6866 (normal)

#2 Visual analysis of residual qqplot
par(mfrow = c(1, 3))
qqnorm(resid(lm.G5pm), main = "No Transformation") # closest to line
qqline(resid(lm.G5pm)) 
qqnorm(resid(loglm.G5pm), main = "Log") 
qqline(resid(loglm.G5pm))
qqnorm(resid(sqrtlm.G5pm), main = "Sqrt") 
qqline(resid(sqrtlm.G5pm))

#3 fitted v. residuals
par(mfrow = c(1, 3))
plot(fitted(lm.G5pm), resid(lm.G5pm)) %>% abline(0,0)
plot(fitted(loglm.G5pm), resid(loglm.G5pm)) %>% abline(0,0)
plot(fitted(sqrtlm.G5pm), resid(sqrtlm.G5pm)) %>% abline(0,0)
# all similar

#4 histogram
par(mfrow = c(1, 3))
hist(resid(lm.G5pm))# most normal looking; but hard to tell
hist(resid(loglm.G5pm)) 
hist(resid(sqrtlm.G5pm)) # may also be normal

######## ANOVA ########
an.G5pm <- anova(lm.G5pm, ddf = "Ken") %>% print()
# No significant data

###### Post Hoc #######
em.G5pm  <- emmeans(lm.G5pm, ~ Condition * Tissue, type = "response") # all estimated means

#em.G5pm  <- emmeans(loglm.G5pm, ~ Condition * Tissue, type = "response") # log lm for ratios
mc.G5pm <- contrast(em.G5pm, "tukey", reverse = T) %>% 
  print()
tis.G5pm <- contrast(em.G5pm, "tukey", reverse = F, by = "Tissue") %>% 
  print()
con.G5pm <- contrast(em.G5pm, "tukey", reverse = F, by = "Condition") %>% 
  print()
# No significant data

####### Plot ##########
max.G5pm <- max(data.frame(em.G5pm)$emmean+data.frame(em.G5pm)$SE)
gg.G5pm.bar <- ggplot(data.frame(em.G5pm), 
                      aes(x=interaction(Condition,Tissue), 
                          y=emmean, fill=interaction(Condition, Tissue))) + 
  coord_cartesian(ylim = c(0, 1.05*max.G5pm), clip = "off") + 
  geom_bar(stat="identity", position="dodge", color='black') + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2, position ='dodge') + 
  # annotate("text", label = paste(round((1/ratio.G5pm$ratio), 2)), # round(ratio.G5pm$SE, 2), sep="\\n±"), # ratio
  #          x=c(2, 4, 6), y=0.05e+08, color = c("black"),
  #          hjust="middle", lineheight=0.85) +
  annotate("text", label = c("Flight\\nMuscle", "Heart", "Liver"), 
           x = c(1.5, 3.5, 5.5), y = 0, vjust = c(2, 4, 4), hjust = 0.4) + # tissue names
  annotate("text", label = c("4", "4", "5", "5", "3", "3"), fontface = "italic", 
           x = c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25), y = 0.055e+08) + # sample sizes
  # geom_segment(aes(x = 2.5, xend = 2.5, y = 0,  yend = 1.1*max.G5pm), linetype=2) +
  # geom_segment(aes(x = 4.5, xend = 4.5, y = 0,  yend = 1.1*max.G5pm), linetype=2) +
  scale_y_continuous("Protein Abundance\\n(arb. intensity units)",
                     position = "right", 
                     labels = scales::comma_format(accuracy=2, scale=0.0000001)) +
  scale_x_discrete(name = "Plasma Membrane Fraction",
                   labels = c("Fed", "Fasted",
                              "Fed", "Fasted", 
                              "Fed", "Fasted"),
                   limits = c("Fed.flightmuscle", "Fasted.flightmuscle", 
                              "Fed.heart", "Fasted.heart", 
                              "Fed.liver", "Fasted.liver")) + 
  scale_fill_manual(values = c("lightpink", "indianred2",
                               "slategray1", "slategray3",
                               "cornsilk", "burlywood1")) +   
  theme(legend.position = "none", #plot.margin = unit(c(2, 1, 1, 1), "cm"), 
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(colour = 0),
        panel.grid.minor = element_line(colour = 0),
        plot.margin=unit(c(0.1, 0.1, 1, 0.1),"cm"),
        # panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        # panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
        axis.title.x = element_text(vjust = 93))
print(gg.G5pm.bar)



############################## Combined Plots ###################################

glut1.fig <- annotate_figure(ggarrange(gg.G1cell.bar, gg.G1pm.bar,
                                       labels = c("A)", "B)"), 
#                                       label.y = 1.075, label.x = c(0.18, 0), 3 for copying
                                        label.y = 1.073, label.x = c(0.175, -0.005), # for PDF
                                       font.label = list(size = 13)),
                             top = text_grob("GLUT1\\n", size=14)) %>%
  print() # Export to clipboard at 600x388 px (6 in x 4.04 in)
glut2.fig <- annotate_figure(ggarrange(gg.G2cell.bar, gg.G2pm.bar, 
                                       labels = c("A)", "B)"), 
                                       label.y = 1.073, label.x = c(0.18, -0.01), 
                                       font.label = list(size = 13)),
                             top = text_grob("GLUT2\\n", size=16)) %>%
  print() # export to clipboard at 600x388
glut3.fig <- annotate_figure(ggarrange(gg.G3cell.bar, gg.G3pm.bar, 
                                       labels = c("A)", "B)"), 
                                       label.y = 1.073, label.x = c(0.16, -0.01), 
                                       font.label = list(size = 13)),
                             top = text_grob("GLUT3\\n", size=16)) %>%
  print()
glut5.fig <- annotate_figure(ggarrange(gg.G5cell.bar, gg.G5pm.bar, 
                                       labels = c("A)", "B)"), 
                                       label.y = 1.073, label.x = c(0.16, -0.01), 
                                       font.label = list(size = 13)),
                             top = text_grob("GLUT5\\n", size=16)) %>%
  print()



GLUTs.fig <- ggarrange(glut1.fig, glut2.fig, glut3.fig, glut5.fig) %>%
  print()
