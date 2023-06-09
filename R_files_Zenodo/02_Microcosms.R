library(ggplot2)
library(doBy) #for summaryBy() function to get mean/sd/se/n
library(car)
library(emmeans)

d <- micro_redo
d <- as.data.frame(d)
d <- subset(d, AFTER_p1 == "Focal")

##### Hurdle Model ######

d <- subset(d, Treatment != "Clip")

# Creating a binomial variable for lesion presence
d$Lesion_P <- ifelse(d$AFTER_pclean_redo> 0, "1","0")
d$Lesion_P <- as.integer(d$Lesion_P)

# Running binomial model on lesion presence
m <- glm(Lesion_P ~ Treatment, data=d, family = binomial)
summary(m)

#Running a linear model on the non-zero lesion values
m <- glm(AFTER_pclean_redo ~ Treatment, data = subset(d, AFTER_pclean_redo > 0), family = Gamma)
summary(m)


#### Boxplot #####

ggplot(d, aes(x = Treatment, y= AFTER_pclean_redo, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1.5, alpha=0.8) +
  theme_classic() +
  geom_hline(yintercept = 0) + 
  theme(legend.position="none")+
  ylab("New Lesion Area (cm2)")+ xlab(NULL) +
  scale_x_discrete(labels=c("Isopod","Clipped", "Control"))+
  scale_fill_manual(values=c("lightsteelblue1", "lightsteelblue3", "lightsteelblue4")) +
  theme(text=element_text(size=10))
ggsave("microcosm_newlesion_boxplot_redo.png", dpi=300, units="in", width=3.5, height=2.5)
