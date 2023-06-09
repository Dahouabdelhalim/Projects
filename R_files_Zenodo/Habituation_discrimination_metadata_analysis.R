
#Analysis of metadata associated with habituation-discrimination playback experiments on golden rocket frogs and Kai rocket frogs

#Associated with paper:
#Tumulty, Lange, Bee. Identity signals, identity reception, and the evolution of social recognition in a Neotropical frog



### Load packages and data ###

library(ggplot2)

data <- read.csv("Habituation_discrimination_metadata.csv")

beebeimeta <- data[data$Species=="beebei",]
kaieimeta <- data[data$Species=="kaiei",]


### Descriptive stats ###

#Time to reach habituation critera
mean(beebeimeta$H_Time)
range(beebeimeta$H_Time)

mean(kaieimeta$H_Time)
range(kaieimeta$H_Time)


#SPL Difference between habituation and discrimination stimuli
print(beebeiSPLdiff <- beebeimeta$H_Stim_SPL - beebeimeta$D_Stim_SPL)
print(kaiSPLdiff <- kaieimeta$H_Stim_SPL - kaieimeta$D_Stim_SPL)

range(abs(beebeiSPLdiff), na.rm = T)
range(abs(kaiSPLdiff), na.rm = T)



### Relationship between time to reach habituation critera and nearest neighbor distance ###

#Plot relationship
ggplot(data, aes(x=N_Dist, y=H_Time, fill=Species))+
  geom_point(pch=21, size=3)+
  geom_smooth(method=lm, se=FALSE, aes(col=Species))+
  xlab("Nearest neighbor distance (m)")+
  ylab("Habituation time (min)")+
  scale_x_continuous(expand=c(0,0), limits=c(0,30))+
  scale_fill_manual(values = c("#EBCC2A", "#3B9AB2"))+
  scale_color_manual(values = c("#EBCC2A", "#3B9AB2"))+
  theme_bw()+
  theme(panel.grid = element_blank())


#Spearman's rank correlation between habituation time and nearest neighbor distance
cor.test(beebeimeta$H_Time, beebeimeta$N_Dist, method="spearman", continuity=FALSE)
cor.test(kaieimeta$H_Time, kaieimeta$N_Dist, method="spearman",  continuity=FALSE)


