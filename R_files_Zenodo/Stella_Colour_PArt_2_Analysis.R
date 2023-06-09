



dfB<-read.csv("WBTest.csv")

hist(dfB$Diff)
qqnorm(dfB$Diff)
qqline(dfB$Diff)

dfB$Treatment<-as.factor(dfB$Treatment)

dfB$Treatment<-ifelse(dfB$Treatment == "1", "Matched Pair", dfB$Treatment)
dfB$Treatment<-ifelse(dfB$Treatment == "2", "Unmatched Pair", dfB$Treatment)
dfB$Treatment<-ifelse(dfB$Treatment == "3", "Solitary", dfB$Treatment)

GlobalBlack <- lmer(Diff~Treatment * I(Time^2) + (1|Trial), data = dfB)
Anova(GlobalBlack)

scatter <- ggplot(dfB, aes(x=Time, y=Diff, colour=Treatment))
scatter + 
  coord_cartesian(ylim=c(-50, 150))+ 
  coord_cartesian(xlim=c(0, 300))+
  geom_point() + stat_smooth(method= "lm",formula = y ~ x + I(x^2),se=TRUE) + 
  labs(x = "Time (s)", y = "Colour Difference")




dfW<-read.csv("BWTest.csv")

hist(dfW$Diff)
qqnorm(dfW$Diff)
qqline(dfW$Diff)

dfW$Treatment<-as.factor(dfW$Treatment)

dfW$Treatment<-ifelse(dfW$Treatment == "4", "Matched Pair", dfW$Treatment)
dfW$Treatment<-ifelse(dfW$Treatment == "2", "Unmatched Pair", dfW$Treatment)
dfW$Treatment<-ifelse(dfW$Treatment == "3", "Solitary", dfW$Treatment)

GlobalWhite <- lmer(Diff~Treatment * I(Time^2) + (1|Trial), data = dfW)
Anova(GlobalWhite)

scatter <- ggplot(dfW, aes(x=Time, y=Diff, colour=Treatment))
scatter + 
  coord_cartesian(ylim=c(-50, 150))+ 
  coord_cartesian(xlim=c(0, 300))+
  geom_point() + stat_smooth(method= "lm",formula = y ~ x + I(x^2),se=TRUE) + 
  labs(x = "Time (s)", y = "Colour Difference")




dfWW<-read.csv("WWTest.csv")

hist(dfWW$Diff)

dfWW$Treatment<-as.factor(dfWW$Treatment)

GlobalWW <- lmer(Diff~Treatment * Time + (1|Trial), data = dfWW)
Anova(GlobalWW)

scatter <- ggplot(dfWW, aes(x=Time, y=Diff, colour=Group))
scatter + 
  coord_cartesian(ylim=c(-50, 150))+ 
  coord_cartesian(xlim=c(0, 300))+
  geom_point() + stat_smooth(method= "lm",se=TRUE) + 
  labs(x = "Time (s)", y = "Colour Difference")


dfBB<-read.csv("BBTest.csv")

hist(dfBB$Diff)

dfBB$Treatment<-as.factor(dfBB$Treatment)

GlobalBB <- lmer(Diff~Treatment * Time + (1|Trial), data = dfBB)
Anova(GlobalBB)

scatter <- ggplot(dfBB, aes(x=Time, y=Diff, colour=Group))
scatter + 
  coord_cartesian(ylim=c(-50, 150))+ 
  coord_cartesian(xlim=c(0, 300))+
  geom_point() + stat_smooth(method= "lm",se=TRUE) + 
  labs(x = "Time (s)", y = "Colour Difference")
