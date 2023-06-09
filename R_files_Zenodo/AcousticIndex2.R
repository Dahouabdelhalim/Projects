#acoustic index results
library(data.table)
library(MuMIn)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(sjPlot)


setwd("~/Dropbox/Papers in progress/Kinabalu Recorder/Data/ManuscriptData")
x101 <- read.csv("SortedDateAcoustciIndexPub.csv")
class(x101$postime)

levels(x101$Site)
length(x101$V1) # n = 126528

#x101[!is.na(x101$Site),]

x5 <- x101[x101$Habitat == "Cloud",]
x6 <- x101[x101$Habitat == "Fagacea",]
length(x5$X)
length(x6$X)

ggplot(x5, aes(x=postime, y=Temperature, color=Site)) +
  geom_point(alpha = 0.2) +
  theme_classic() +
  scale_color_manual(values = safe_cols) +
  scale_x_discrete(breaks=c("0000-01-01 00:00:00","0000-01-01 06:00:00","0000-01-01 12:00:00","0000-01-01 18:00:00","0000-01-01 23:50:00"), labels=c("00:00", "06:00", "12:00","18:00","24:00"))+
  ylim(c(5,30))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title.x = element_blank(),legend.title = element_blank(),legend.position = "top")


ggplot(x6, aes(x=postime, y=Temperature, color=Site)) +
  geom_point(alpha = 0.2) +
  theme_classic() +
  scale_color_manual(values = safe_cols) +
  scale_x_discrete(breaks=c("0000-01-01 00:00:00","0000-01-01 06:00:00","0000-01-01 12:00:00","0000-01-01 18:00:00","0000-01-01 23:50:00"), labels=c("00:00", "06:00", "12:00","18:00","24:00"))+
  ylim(c(5,30))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title.x = element_blank(),legend.title = element_blank(),legend.position = "top")

#mixed model#############################################
#mixed model#############################################
#mixed model#############################################
library(lme4)
#correct pressure values that are incorrect
x101[x101$Pressure > 120000, ] <- NA
tax <- x101[,c("Temperature","Pressure")]
tax <- scale(tax)
colnames(tax) <- c("STemperature","SPressure")
x101 <- cbind(x101, tax)

x4 <- x101

set.seed(17778)
#NormalisedDifferenceSoundscape#############################################
model1 <- lmer(NormalisedDifferenceSoundscape ~ Habitat + STemperature + SPressure + Month +(1|Site) ,data = x4)

model2 <- lmer(NormalisedDifferenceSoundscape ~ Habitat + STemperature + SPressure + (1|Site) ,data = x4)

model3 <- lmer(NormalisedDifferenceSoundscape ~ Habitat + STemperature  +(1|Site) ,data = x4)

model4 <- lmer(NormalisedDifferenceSoundscape ~ Habitat  +(1|Site) ,data = x4)

model5 <- lmer(NormalisedDifferenceSoundscape ~ STemperature +(1|Site) ,data = x4)

model6 <- lmer(NormalisedDifferenceSoundscape ~  SPressure  +(1|Site) ,data = x4)

model7 <- lmer(NormalisedDifferenceSoundscape ~  Month +(1|Site) ,data = x4)

model8 <- lmer(NormalisedDifferenceSoundscape ~ 1 +(1|Site) ,data = x4)

AIC(model1, model2, model3, model4, model5, model6, model7,model8)

sink("NormalisedDifferenceSoundscape.txt")
sjPlot::plot_model(model1)
r.squaredGLMM(model1)
model1
summary(model1)
sink()
sjPlot::plot_model(model1, pred.type = "re")
plot(NormalisedDifferenceSoundscape ~ Temperature, data = x3)
abline(model4)

normalised <- model1

#RightBioacousticIndex#############################################
model1 <- lmer(BioacousticIndex ~ Habitat + STemperature + SPressure + Month +(1|Site) ,data = x4)

model2 <- lmer(BioacousticIndex ~ Habitat + STemperature + SPressure + (1|Site) ,data = x4)

model3 <- lmer(BioacousticIndex ~ Habitat + STemperature  +(1|Site) ,data = x4)

model4 <- lmer(BioacousticIndex ~ Habitat  +(1|Site) ,data = x4)

model5 <- lmer(BioacousticIndex ~ STemperature +(1|Site) ,data = x4)

model6 <- lmer(BioacousticIndex ~  SPressure  +(1|Site) ,data = x4)

model7 <- lmer(BioacousticIndex ~  Month +(1|Site) ,data = x4)

model8 <- lmer(BioacousticIndex ~ 1 +(1|Site) ,data = x4)

AIC(model1,  model2, model3, model4, model5, model6, model7,model8)

BioA <- model1

sink("BioacousticIndex.txt")
sjPlot::plot_model(model1)
r.squaredGLMM(model1)
summary(model1)
model1
sink()
sjPlot::plot_model(model2, pred.type = "re")
plot(RightBioacousticIndex ~ Temperature, data = x3)
abline(model4)

ggplot(data = x3, aes(x = Temperature , y =RightBioacousticIndex)) +
  geom_smooth(method = "lm")


normalised
BioA


p <- list()

p1 <- plot_model(normalised,show.values = TRUE,vline.color = "black", value.size = 4, col = "black",title = "Normalised Difference Index") +
  theme_classic()  
  #theme(axis.text.y=element_blank())

p2 <- plot_model(BioA,vline.color = "black", show.values = TRUE,value.size = 4, col = "black",title = "Bioacoustic Index") + 
  theme_classic() +
  theme(axis.text.y=element_blank())


p4 <-ggplot(data = x4, aes(x = Pressure , y =BioacousticIndex)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess", col = "black")  +
  xlab("Barometric Pressure (Pa)") +
  ylab("Bioacoustic Index") +
  theme_classic()

title  = expression("Temperature ("*~degree*C*")")
p3 <- ggplot(data = x4, aes(x = Temperature , y =NormalisedDifferenceSoundscape)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess", col = "white")  +
  xlab(title) +
  ylab("Normalised Difference Index") +
  theme_classic()

library(gridExtra)
library(ggpubr)

p[[1]] <- p1
p[[2]] <- p2
p[[3]] <- p3
p[[4]] <- p4

setwd("~/Dropbox/Papers in progress/Kinabalu Recorder/Figures")
ggarrange(plotlist = p[c(1:4)],nrow = 2, ncol = 2)
ggsave("MixedModelresult.png",width = 20, height =15, units= "cm",  dpi = "print")

sink("LMERresults.txt")
print("#########Normalised#########")
print("Summary")
summary(normalised)
print("model results")
normalised
print("r^2")
r.squaredGLMM(normalised)
print("Coefficients")
coef(summary(normalised))
print("Confidence")
confint(normalised)
print("#########Bioacoustic#########")
print("Summary")
summary(BioA)
print("model results")
BioA
print("r^2")
r.squaredGLMM(BioA)
print("Coefficients")
coef(summary(BioA))
print("Confidence")
confint(BioA)
sink()
