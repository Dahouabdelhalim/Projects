##########################################################################################
# Vibrations emitted in presence of reproductives affect social organization in termites #                                                      #
##########################################################################################

# Louis PAILLER 

###########################
# PACKAGES LOADING #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###########################
#library(MASS)
library(multcomp)
library(multcompView)
library(lme4)
library(car) 
library(emmeans)
library(ggplot2)
library(DHARMa)
library(reshape2)
library(gridExtra)

###########################
# DATA LOADING #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###########################
data_distibution <- read.csv("Data_distribution.csv", sep = ";")


d <- read.csv("Data_body_shaking.csv", sep = ";")


data_identity <- read.csv("Data_identity.csv", sep = ";")


dc <- read.csv("Data_comparison_non_bs_vs_bs.csv", sep = ";")
d2_during <- subset(dc, dc$timing=="during")
d2_during$treatment <- factor(d2_during$treatment, levels=c("Control_W", "Control_WR", "Original_W", "Original_WR", "Randomised_W", "Randomised_WR"))




#---------------------------------------------------------------------------------------------------------------------------------------------------------------
##############################
### DISTRIBUTION BODY-shaker
### FIGURE SUPPLEMENTAL 1
##############################
data_distribution_during <- subset(data_distibution, data_distibution$timing=="during")
data_distribution_after <- subset(data_distibution, data_distibution$timing=="after")


ggplot(data_distribution_during, aes(x=nb_bs)) +
  geom_histogram(bins = 40) + 
  scale_x_continuous(breaks=seq(0, 40, 1)) +
  scale_y_continuous(breaks=seq(0, 1375, 25))

ggplot(data_distribution_after, aes(x=nb_bs)) +
  geom_histogram(bins = 40) + 
  scale_x_continuous(breaks=seq(0, 40, 1)) +
  scale_y_continuous(breaks=seq(0, 1375, 25))




#---------------------------------------------------------------------------------------------------------------------------------------------------------------




#---------------------------------------------------------------------------------------------------------------------------------------------------------------
##############################
### OCCURRENCE OF BODY SHAKING
### FIGURE 2
##############################
d$log_nb_total_bs <- log(d$nb_total_bs+1)
d1_during <- subset(d, d$timing=="during")
d1_during$treatment <- factor(d1_during$treatment, levels=c("Control_W", "Control_WR", "Original_W", "Original_WR", "Randomised_W", "Randomised_WR"))



# TREATMENT

model1 <- lmer(log_nb_total_bs ~ treatment + (1 | colony_id),  data = d1_during)
hist(resid(model1))
testResiduals(model1)
anova <- Anova(model1)
anova
emmeans_model1 <- emmeans(model1, pairwise ~ treatment, adjust="tukey")
emmeans_model1

model1_cld <- cld(object = emmeans_model1, adjust = "Tukey", Letters = letters, alpha = 0.05)
model1_cld


ggplot(d1_during, aes(x=treatment, log_nb_total_bs)) +
  geom_boxplot() +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0.2), size=1.5) +
  ylim(0, 6) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# DURING vs. AFTER
d1_Original_W <- subset(d, (treatment=="Original_W"))
d1_Original_WR <- subset(d, (treatment=="Original_WR"))
d1_Randomised_W <- subset(d, (treatment=="Randomised_W"))
d1_Randomised_WR <- subset(d, (treatment=="Randomised_WR"))


ggplot(d1_Original_W, aes(x=timing, log_nb_total_bs)) +
  geom_boxplot() +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0.2), size=1.5) +
  ylim(0, 6) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


ggplot(d1_Original_WR, aes(x=timing, log_nb_total_bs)) +
  geom_boxplot() +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0.2), size=1.5) +
  ylim(0, 6) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


ggplot(d1_Randomised_W, aes(x=timing, log_nb_total_bs)) +
  geom_boxplot() +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0.2), size=1.5) +
  ylim(0, 6) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


ggplot(d1_Randomised_WR, aes(x=timing, log_nb_total_bs)) +
  geom_boxplot() +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0.2), size=1.5) +
  ylim(0, 6) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


t.test(log_nb_total_bs~timing, data=d1_Original_W, paired=TRUE, conf.level=0.95)

t.test(log_nb_total_bs~timing, data=d1_Original_WR, paired=TRUE, conf.level=0.95)

t.test(log_nb_total_bs~timing, data=d1_Randomised_W, paired=TRUE, conf.level=0.95)

t.test(log_nb_total_bs~timing, data=d1_Randomised_WR, paired=TRUE, conf.level=0.95)




##########################
### NUMBER OF BODY-shaker
# FIGURE 3
##########################

# TREATMENT
model2 <- lmer(nb_total_body_shaker ~ treatment + (1 | colony_id),  data = d1_during)
hist(resid(model2))
testResiduals(model2)
anova <- Anova(model2)
anova
emmeans_model2 <- emmeans(model2, pairwise ~ treatment, adjust="tukey")
emmeans_model2

model2_cld <- cld(object = emmeans_model2, adjust = "Tukey", Letters = letters, alpha = 0.05)
model2_cld


ggplot(d1_during, aes(x=treatment, nb_total_body_shaker)) +
  geom_boxplot() +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0.2), size=1.5) +
  ylim(0, 20) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



# DURING vs. AFTER

ggplot(d1_Original_W, aes(x=timing, nb_total_body_shaker)) +
  geom_boxplot() +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0.2), size=1.5) +
  ylim(0, 20) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


ggplot(d1_Original_WR, aes(x=timing, nb_total_body_shaker)) +
  geom_boxplot() +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0.2), size=1.5) +
  ylim(0, 20) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


ggplot(d1_Randomised_W, aes(x=timing, nb_total_body_shaker)) +
  geom_boxplot() +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0.2), size=1.5) +
  ylim(0, 20) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


ggplot(d1_Randomised_WR, aes(x=timing, nb_total_body_shaker)) +
  geom_boxplot() +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0.2), size=1.5) +
  ylim(0, 20) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


t.test(nb_total_body_shaker~timing, data=d1_Original_W, paired=TRUE, conf.level=0.95)

t.test(nb_total_body_shaker~timing, data=d1_Original_WR, paired=TRUE, conf.level=0.95)

t.test(nb_total_body_shaker~timing, data=d1_Randomised_W, paired=TRUE, conf.level=0.95)

t.test(nb_total_body_shaker~timing, data=d1_Randomised_WR, paired=TRUE, conf.level=0.95)




#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#######################
### CORRELATION PLOT X6
### FIGURE SUPPLEMENTAL 3
#######################

# Control R-

d_Control_W <- subset(d, d$treatment=="Control_W")

d_Control_W_during <- subset(d_Control_W, d_Control_W$timing=="during")
d_Control_W_after <- subset(d_Control_W, d_Control_W$timing=="after")

res_during <- cor.test(d_Control_W_during$nb_total_body_shaker, d_Control_W_during$log_nb_total_bs, method = "pearson")
res_during


ggplot(d_Control_W_during, aes(nb_total_body_shaker, log_nb_total_bs)) +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0), size=2) +
  geom_smooth(method = "lm", se = FALSE) +
  ylim(0, 6) +
  xlim(0, 11) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())




# Control R+
d_Control_WR <- subset(d, d$treatment=="Control_WR")
d_Control_WR_during <- subset(d_Control_WR, d_Control_WR$timing=="during")
res_during <- cor.test(d_Control_WR_during$nb_total_body_shaker, d_Control_WR_during$log_nb_total_bs, method = "pearson")
res_during


ggplot(d_Control_WR_during, aes(nb_total_body_shaker, log_nb_total_bs)) +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0), size=2) +
  geom_smooth(method = "lm", se = FALSE) +
  ylim(0, 6) +
  xlim(0, 11) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



# Original R-
## During playback
d_Original_W <- subset(d, d$treatment=="Original_W")
d_Original_W_during <- subset(d_Original_W, d_Original_W$timing=="during")
d_Original_W_after <- subset(d_Original_W, d_Original_W$timing=="after")
res_during <- cor.test(d_Original_W_during$nb_total_body_shaker, d_Original_W_during$log_nb_total_bs, method = "pearson")
res_during

res_after <- cor.test(d_Original_W_after$nb_total_body_shaker, d_Original_W_after$log_nb_total_bs, method = "pearson")
res_after


ggplot(d_Original_W_during, aes(nb_total_body_shaker, log_nb_total_bs)) +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0), size=2) +
  geom_smooth(method = "lm", se = FALSE) +
  ylim(0, 6) +
  xlim(0, 11) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



## After playback
ggplot(d_Original_W_after, aes(nb_total_body_shaker, log_nb_total_bs)) +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0), size=2) +
  geom_smooth(method = "lm", se = FALSE) +
  ylim(0, 6) +
  xlim(0, 11) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())




# Original R+
## During playback

d_Original_WR <- subset(d, d$treatment=="Original_WR")
d_Original_WR_during <- subset(d_Original_WR, d_Original_WR$timing=="during")
d_Original_WR_after <- subset(d_Original_WR, d_Original_WR$timing=="after")
res_during <- cor.test(d_Original_WR_during$nb_total_body_shaker, d_Original_WR_during$log_nb_total_bs, method = "pearson")
res_during

res_after <- cor.test(d_Original_WR_after$nb_total_body_shaker, d_Original_WR_after$log_nb_total_bs, method = "pearson")
res_after


ggplot(d_Original_WR_during, aes(nb_total_body_shaker, log_nb_total_bs)) +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0), size=2) +
  geom_smooth(method = "lm", se = FALSE) +
  ylim(0, 6) +
  xlim(0, 11) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



## After playback
ggplot(d_Original_WR_after, aes(nb_total_body_shaker, log_nb_total_bs)) +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0), size=2) +
  geom_smooth(method = "lm", se = FALSE) +
  ylim(0, 6) +
  xlim(0, 11) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



# Randomised R-
## During playback

d_Randomised_W <- subset(d, d$treatment=="Randomised_W")
d_Randomised_W_during <- subset(d_Randomised_W, d_Randomised_W$timing=="during")
d_Randomised_W_after <- subset(d_Randomised_W, d_Randomised_W$timing=="after")
res_during <- cor.test(d_Randomised_W_during$nb_total_body_shaker, d_Randomised_W_during$log_nb_total_bs, method = "pearson")
res_during

res_after <- cor.test(d_Randomised_W_after$nb_total_body_shaker, d_Randomised_W_after$log_nb_total_bs, method = "pearson")
res_after


ggplot(d_Randomised_W_during, aes(nb_total_body_shaker, log_nb_total_bs)) +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0), size=2) +
  geom_smooth(method = "lm", se = FALSE) +
  ylim(0, 6) +
  xlim(0, 11) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



## After playback
ggplot(d_Randomised_W_after, aes(nb_total_body_shaker, log_nb_total_bs)) +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0), size=2) +
  geom_smooth(method = "lm", se = FALSE) +
  ylim(0, 6) +
  xlim(0, 11) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())




# Randomised R+
## During playback

d_Randomised_WR <- subset(d, d$treatment=="Randomised_WR")
d_Randomised_WR_during <- subset(d_Randomised_WR, d_Randomised_WR$timing=="during")
d_Randomised_WR_after <- subset(d_Randomised_WR, d_Randomised_WR$timing=="after")
res_during <- cor.test(d_Randomised_WR_during$nb_total_body_shaker, d_Randomised_WR_during$log_nb_total_bs, method = "pearson")
res_during

res_after <- cor.test(d_Randomised_WR_after$nb_total_body_shaker, d_Randomised_WR_after$log_nb_total_bs, method = "pearson")
res_after


ggplot(d_Randomised_WR_during, aes(nb_total_body_shaker, log_nb_total_bs)) +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0), size=2) +
  geom_smooth(method = "lm", se = FALSE) +
  ylim(0, 6) +
  xlim(0, 11) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



## After playback
ggplot(d_Randomised_WR_after, aes(nb_total_body_shaker, log_nb_total_bs)) +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0), size=2) +
  geom_smooth(method = "lm", se = FALSE) +
  ylim(0, 6) +
  xlim(0, 11) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


#---------------------------------------------------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------------------------------------------------------------
######################################
### BODY SHAKER IDENTITY
# FIGURE 4
######################################
data_identity <- read.csv("Data_identity.csv", sep = ";")


data_identity1 <- subset(data_identity, !(treatment=="Control_W"))
data_identity1 <- subset(data_identity1, !(treatment=="Control_WR"))

data_identity1$treatment <- factor(data_identity1$treatment, levels=c("Original_W", "Original_WR", "Randomised_W", "Randomised_WR"))
data_identity1$identity <- factor(data_identity1$identity, levels=c("never", "changing", "always"))



# Never
data_identity2 <- subset(data_identity1, data_identity1$identity=="never")

ggplot(data_identity2, aes(treatment, count)) +
  geom_boxplot() +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0.1), size=1.5) +
  ylim(0, 20) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


model3 <- lmer(count ~ treatment + (1 | colony_id),  data = data_identity2)
hist(resid(model3))
testResiduals(model3)
anova <- Anova(model3)
anova
emmeans_model3 <- emmeans(model3, pairwise ~ treatment, adjust="tukey")
emmeans_model3

model3_cld <- cld(object = emmeans_model3, adjust = "Tukey", Letters = letters, alpha = 0.05)
model3_cld



# Changing
data_identity2 <- subset(data_identity1, data_identity1$identity=="changing")

ggplot(data_identity2, aes(treatment, count)) +
  geom_boxplot() +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0.1), size=1.5) +
  ylim(0, 20) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

model4 <- lmer(count ~ treatment + (1 | colony_id),  data = data_identity2)
hist(resid(model4))
testResiduals(model4)
anova <- Anova(model4)
anova
emmeans_model4 <- emmeans(model4, pairwise ~ treatment, adjust="tukey")
emmeans_model4

model4_cld <- cld(object = emmeans_model4, adjust = "Tukey", Letters = letters, alpha = 0.05)
model4_cld


# Always
data_identity2 <- subset(data_identity1, data_identity1$identity=="always")

ggplot(data_identity2, aes(treatment, count)) +
  geom_boxplot() +
  geom_jitter(shape = 1, color = "black", position = position_jitter(width = 0.1), size=1.5) +
  ylim(0, 20) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

model5 <- lmer(count ~ treatment + (1 | colony_id),  data = data_identity2)
hist(resid(model5))
testResiduals(model5)
anova <- Anova(model5)
anova





#---------------------------------------------------------------------------------------------------------------------------------------------------------------
######################################
### NON BS vs. BS DISTANCE ENTRE INDIV
# Figure 5
######################################


### DURING
ggplot(d2_during, aes(x=treatment, dist_between_indiv_non_bs_vs_bs, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 0.04) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")



model6 <- lmer(dist_between_indiv_non_bs_vs_bs ~ treatment+non_bs_vs_bs + (1 | colony_id), data = d2_during)
testResiduals(model6)
anova <- Anova(model6)
anova
emmeans_model6 <- emmeans(model6, pairwise ~ treatment, adjust="tukey")
emmeans_model6

model6_cld <- cld(object = emmeans_model6, adjust = "Tukey", Letters = letters, alpha = 0.05)
model6_cld





# DURING vs. AFTER

d_Original_W <- subset(dc, (treatment=="Original_W"))
d_Original_W$timing <- factor(d_Original_W$timing, levels=c("during", "after"))
d_Original_WR <- subset(dc, (treatment=="Original_WR"))
d_Original_WR$timing <- factor(d_Original_WR$timing, levels=c("during", "after"))
d_Randomised_W <- subset(dc, (treatment=="Randomised_W"))
d_Randomised_W$timing <- factor(d_Randomised_W$timing, levels=c("during", "after"))
d_Randomised_WR <- subset(dc, (treatment=="Randomised_WR"))
d_Randomised_WR$timing <- factor(d_Randomised_WR$timing, levels=c("during", "after"))


ggplot(d_Original_W, aes(x=timing, dist_between_indiv_non_bs_vs_bs, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 0.04) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

ggplot(d_Original_WR, aes(x=timing, dist_between_indiv_non_bs_vs_bs, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 0.04) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

ggplot(d_Randomised_W, aes(x=timing, dist_between_indiv_non_bs_vs_bs, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 0.04) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

ggplot(d_Randomised_WR, aes(x=timing, dist_between_indiv_non_bs_vs_bs, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 0.04) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")



model7 <- lmer(dist_between_indiv_non_bs_vs_bs ~ non_bs_vs_bs+timing + (1 | colony_id), data = d_Original_W)
testResiduals(model7)
anova <- Anova(model7)
anova

model8 <- lmer(dist_between_indiv_non_bs_vs_bs ~ non_bs_vs_bs+timing + (1 | colony_id), data = d_Original_WR)
testResiduals(model8)
anova <- Anova(model8)
anova

model9 <- lmer(dist_between_indiv_non_bs_vs_bs ~ non_bs_vs_bs+timing + (1 | colony_id), data = d_Randomised_W)
testResiduals(model9)
anova <- Anova(model9)
anova

model10 <- lmer(dist_between_indiv_non_bs_vs_bs ~ non_bs_vs_bs+timing + (1 | colony_id), data = d_Randomised_WR)
testResiduals(model10)
anova <- Anova(model10)
anova




#---------------------------------------------------------------------------------------------------------------------------------------------------------------
##################################
### NON BS vs. BS DISTANCE COVERED
# Figure 6
##################################

dc$log_dist_covered <- log(dc$distance_covered_indiv_non_bs_vs_bs+1)
d1_during <- subset(dc, dc$timing=="during")
d1_during$treatment <- factor(d1_during$treatment, levels=c("Control_W", "Control_WR", "Original_W", "Original_WR", "Randomised_W", "Randomised_WR"))


### DURING
ggplot(d1_during, aes(x=treatment, log_dist_covered, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 1.5) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


model11 <- lmer(log_dist_covered ~ treatment+non_bs_vs_bs + (1 | colony_id), data = d1_during)
testResiduals(model11)
anova <- Anova(model11)
anova


d_Original_W <- subset(dc, (treatment=="Original_W"))
d_Original_W$timing <- factor(d_Original_W$timing, levels=c("during", "after"))
d_Original_WR <- subset(dc, (treatment=="Original_WR"))
d_Original_WR$timing <- factor(d_Original_WR$timing, levels=c("during", "after"))
d_Randomised_W <- subset(dc, (treatment=="Randomised_W"))
d_Randomised_W$timing <- factor(d_Randomised_W$timing, levels=c("during", "after"))
d_Randomised_WR <- subset(dc, (treatment=="Randomised_WR"))
d_Randomised_WR$timing <- factor(d_Randomised_WR$timing, levels=c("during", "after"))


ggplot(d_Original_W, aes(x=timing, log_dist_covered, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 1.5) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

ggplot(d_Original_WR, aes(x=timing, log_dist_covered, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 1.5) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

ggplot(d_Randomised_W, aes(x=timing, log_dist_covered, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 1.5) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

ggplot(d_Randomised_WR, aes(x=timing, log_dist_covered, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 1.5) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")




model12 <- lmer(log_dist_covered ~ non_bs_vs_bs+timing + (1 | colony_id), data = d_Original_W)
testResiduals(model12)
anova <- Anova(model12)
anova

model13 <- lmer(log_dist_covered ~ non_bs_vs_bs+timing + (1 | colony_id), data = d_Original_WR)
testResiduals(model13)
anova <- Anova(model13)
anova

model14 <- lmer(log_dist_covered ~ non_bs_vs_bs+timing + (1 | colony_id), data = d_Randomised_W)
testResiduals(model14)
anova <- Anova(model14)
anova

model15 <- lmer(log_dist_covered ~ non_bs_vs_bs+timing + (1 | colony_id), data = d_Randomised_WR)
testResiduals(model15)
anova <- Anova(model15)
anova



#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#########################################
### DISTANCE PIEZOELECTRIC BS vs. NON BS
# FIGURE SUPPLEMENTAL 6
#########################################

d1_during <- subset(dc, dc$timing=="during")
d1_during$treatment <- factor(d1_during$treatment, levels=c("Control_W", "Control_WR", "Original_W", "Original_WR", "Randomised_W", "Randomised_WR"))


### DURING
ggplot(d1_during, aes(x=treatment, dist_between_indiv_piezo_non_bs_vs_bs, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 0.03) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


model16 <- lmer(dist_between_indiv_piezo_non_bs_vs_bs ~ treatment+non_bs_vs_bs + (1 | colony_id), data = d1_during)
testResiduals(model16)
anova <- Anova(model16)
anova


d_Original_W <- subset(dc, (treatment=="Original_W"))
d_Original_W$timing <- factor(d_Original_W$timing, levels=c("during", "after"))
d_Original_WR <- subset(dc, (treatment=="Original_WR"))
d_Original_WR$timing <- factor(d_Original_WR$timing, levels=c("during", "after"))
d_Randomised_W <- subset(dc, (treatment=="Randomised_W"))
d_Randomised_W$timing <- factor(d_Randomised_W$timing, levels=c("during", "after"))
d_Randomised_WR <- subset(dc, (treatment=="Randomised_WR"))
d_Randomised_WR$timing <- factor(d_Randomised_WR$timing, levels=c("during", "after"))


ggplot(d_Original_W, aes(x=timing,  dist_between_indiv_piezo_non_bs_vs_bs, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 0.03) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

ggplot(d_Original_WR, aes(x=timing, dist_between_indiv_piezo_non_bs_vs_bs, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 0.03) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

ggplot(d_Randomised_W, aes(x=timing, dist_between_indiv_piezo_non_bs_vs_bs, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 0.03) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

ggplot(d_Randomised_WR, aes(x=timing, dist_between_indiv_piezo_non_bs_vs_bs, fill=non_bs_vs_bs)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7) , shape = 1, size = 1.6) +
  ylim(0, 0.03) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold")) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme( axis.line = element_line(colour = "black", size = 0.6,linetype = "solid")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic")) +
  theme(strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) + 
  theme(panel.background = element_rect(fill = "white",colour = "white", size = 0, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=10)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")



model17 <- lmer(dist_between_indiv_piezo_non_bs_vs_bs ~ non_bs_vs_bs+timing + (1 | colony_id), data = d_Original_W)
testResiduals(model17)
anova <- Anova(model17)
anova

model18 <- lmer(dist_between_indiv_piezo_non_bs_vs_bs ~ non_bs_vs_bs+timing + (1 | colony_id), data = d_Original_WR)
testResiduals(model18)
anova <- Anova(model18)
anova

model19 <- lmer(dist_between_indiv_piezo_non_bs_vs_bs ~ non_bs_vs_bs+timing + (1 | colony_id), data = d_Randomised_W)
testResiduals(model19)
anova <- Anova(model19)
anova

model20 <- lmer(dist_between_indiv_piezo_non_bs_vs_bs ~ non_bs_vs_bs+timing + (1 | colony_id), data = d_Randomised_WR)
testResiduals(model20)
anova <- Anova(model20)
anova


