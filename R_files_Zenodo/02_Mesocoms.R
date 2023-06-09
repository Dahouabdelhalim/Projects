library(ggplot2)
library(emmeans)

d <- mesocosm
d <- as.data.frame(d)
d <- subset(d, AFTER_Lesions != 44)

##### Lesion Area Change #######

#Stats
m <- lm(Lesion_Area_Change ~ Treatment, data=d)
summary(m)
emm <- emmeans(m, ~Treatment)
contrast(emm, method = "tukey")
shapiro.test(residuals(m))

hist(residuals(m))

#There is an outlier
cooksd <-- cooks.distance(m)

ggplot(d, aes(x=Treatment, y=Lesion_Area_Change)) +
  geom_boxplot(fill = "olivedrab3", color = "black") 

#It is Clip8, which I will go ahead and remove 

d <- subset(d, AFTER_Lesions != 44)
#Residuals now normal. 
#Stats re-run with Clip8 removed from the dataset

#want a box plot to see the negative change values for Bug 

ggplot(d, aes(x = Treatment, y= Lesion_Area_Change, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1.5, alpha=0.8) +
  theme_classic() +
  geom_hline(yintercept = 0) + 
  theme(legend.position="none")+
  ylab("Change in Total Lesion Area (cm2)")+ xlab(NULL) +
  scale_x_discrete(labels=c("Isopod","Clipped", "Control"))+
  scale_fill_manual(values=c("lightsteelblue1", "lightsteelblue3", "lightsteelblue4")) +
  theme(text=element_text(size=10))
ggsave("Lesion_Area_Mesocosms_BoxPlot.png", dpi=300, units="in", width=3.5, height=2.5)

