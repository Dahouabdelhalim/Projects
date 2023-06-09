library(ggplot2)
library(doBy) #for summaryBy() function to get mean/sd/se/n
library(lme4) # mixed effects models
library(Hotelling)


##### P. resecata feeding trial #####
d <- subset(feeding_trial_FINAL)
d <- as.data.frame(d)

ggplot(d, aes(x = Leaf_Health, y= Actual_Area_Consumed, fill = Leaf_Health)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1.5, alpha=0.8) +
  theme_classic() +
  theme(legend.position="none")+
  ylab("Area Consumed (cm2)")+ xlab(NULL) +
  scale_fill_manual(values=c("darkolivegreen3","darkgoldenrod4")) +
  scale_x_discrete(labels=c("Green","Lesion"))+
  theme(text=element_text(size=12))
ggsave("Field_Area_Consumed_BOXPLOT.png", dpi=300, units="in", width=3.5, height=2.5)


#Consumption Stats for Actual_Area_Consumed
m <- lm(Actual_Area_Consumed ~ Leaf_Health, data=d)
summary(m)


#checking residuals
#residuals should be greater than 0.05 
shapiro.test(residuals(m1))
hist(residuals(m1))

##### Penetrometer Figure ########

d <- penetrometer_2020
d <- as.data.frame(d)

ggplot(d, aes(x = Leaf_Health, y = Penetrometer_Mass, fill = Leaf_Health)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1.5, alpha=0.8) +
  theme_classic() +
  theme(legend.position="none")+
  ylab("Penetrometer Mass (mg)")+ xlab(NULL) +
  scale_fill_manual(values=c("darkolivegreen3","darkgoldenrod4")) +
  scale_x_discrete(labels=c("Green","Lesion"))+
  theme(text=element_text(size=12))
ggsave("Penetrometer_Field_Trial_BOXPLOT.png", dpi=300, units="in", width= 3.5, height= 2.5)
       

#Penetrometer Stats 
m <- lm(Penetrometer_Mass ~ Leaf_Health, data=d)
summary(m)
shapiro.test(residuals(m))


##### A. valida feeding trial #####

d <- valida_data
d <- as.data.frame(d)

##### Percent of Total Consumption #######

ggplot(d, aes(x = Treatment, y= Percent_Consumed, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1.5, alpha=0.8) +
  theme_classic() +
  theme(legend.position="none")+
  ylab("% of Total COnsumption")+ xlab(NULL) +
  scale_fill_manual(values=c("darkolivegreen3","darkgoldenrod4")) +
  scale_x_discrete(labels=c("Green","Lesion"))+
  theme(text=element_text(size=12))
ggsave("Valida_Percent_Total_Consumed_BOXPLOT.png", dpi=300, units="in", width=3.5, height=2.5)


m <- lm(Percent_Consumed ~ Treatment, data=d)
summary(m)

#checking residuals
#residuals should be greater than 0.05 
shapiro.test(residuals(m1))
hist(residuals(m1))



