library(ggplot2)
library(lme4)
library(dplyr)
library(tidyr)
library(car)
library(AICcmodavg)
library(effects)
library(lattice)  
library(agricolae)

##########################################

#provisioning
fish=read.csv("Byerlyetal2020_provisioningdata1.csv") 
fullfish=read.csv("Byerlyetal2020_provisioningdata2.csv")

#shealer index
mean(fullfish$shealer_mass, na.rm=TRUE)
sd(fullfish$shealer_mass, na.rm=TRUE)

length=as.factor(fullfish$shealer_length)
mean(fullfish$shealer_length, na.rm=TRUE)
sd(fullfish$shealer_length, na.rm=TRUE)


#split-plot ANOVA
#model effect of covariates on hourly feeding rate (by day)
hourly_mod=aov(hourly ~ day * island + nest, data=fish)
summary(hourly_mod)
day=as.factor(fish$day)
plot(day, fish$hourly)

#model effect of covariates on time of day of feeding (full feeding history)
full_mod=aov(hour ~ day + island * island + nest, data=fullfish)
summary(full_mod)

#model effect of covariates on prey volume (full feeding history)
#vol=as.integer(fullfish$shealer_mass, na.remove=TRUE)
volume_mod=aov(fullfish$shealer_mass ~ day + island + hour * island + nest, data=fullfish)
summary(volume_mod)
plot(fullfish$island, fullfish$shealer_mass)

#summaries by island
summary=fullfish%>%
  group_by(island) %>%
  summarise(size=mean(shealer_mass, na.rm=TRUE),
            SD=sd(shealer_mass, na.rm=TRUE))

#plot 1
day=as.factor(fish1$day)
#figures for chick provisioning rate by day
rate=fish1 %>% 
  ggplot(aes(x=day, y=fish1$hourly)) +
  geom_boxplot (aes(group=day))+
  labs(title ="a.", 
       x = "Days post-hatch", y = "Fish per hour\\n") + 
  geom_jitter(alpha = 0.3, color = "black") +
  pubtheme
plot(rate)

#plot 2
mass=as.factor(fullfish$shealer_mass)
#percent prey category by island: stacked bar chart
data=data.frame(fullfish$shealer_mass, fullfish$island)
island=ggplot(data, aes(x=fullfish$island, y=fullfish$shealer_mass, fill=mass)) + 
  geom_bar(position="fill", stat="identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title ="b.", x="Island", y = "% Fish Mass\\n") + pubtheme + 
  scale_fill_grey (name="Index Wet Mass", labels=c("tiny: <0.10 g", "small: 0.1 - 0.4 g", "medium: 0.4 - 1.2 g", "large: >1.2 g"))

grid.arrange(rate, island, nrow=2, ncol=1) 