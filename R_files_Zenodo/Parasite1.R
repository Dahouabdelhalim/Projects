#Code by Caleb Nielebeck to analyze survivorship and questing data in Amblyomma americanum, Dermacentor variabilis, and Ixodes scapularis at 20-30 degrees C

# Import survivorship data
library(readxl)
Survivorship <- read_excel("Survivorship.xlsx")

#AIC
library(AICcmodavg)
#individual models~survival
Species.mod <- lm(Time ~ Species, data = Survivorship)
RH.mod <- lm(Time ~ RH, data = Survivorship)
#Combination and interaction models
combination.mod <- lm(Time ~ Species + RH, data = Survivorship)
interaction.mod <- lm(Time ~ Species*RH, data = Survivorship)
#prepare the models for AIC
models <- list(Species.mod, RH.mod, combination.mod, interaction.mod)
model.names <- c('Species.mod', 'RH.mod', 'combination.mod', 'interaction.mod')
#AICc Analysis and table
aictab(cand.set = models, modnames = model.names)
#individual model summary with R2
summary(Species.mod)
summary(RH.mod)
summary(combination.mod)
summary(interaction.mod)

#Kaplan-Meier analysis
library(survival)
library(survminer)
#First setup survival objects
surv <- Surv(time = Survivorship[['Time']], event = Survivorship[['Status']])
#Create Object
surv.rh <- survfit(surv ~ Species+RH, data = Survivorship, type = 'kaplan-meier')
print(surv.rh,print.rmean=TRUE)
#Pairwise log-rank test
data <- pairwise_survdiff(Surv(Time, Status) ~ Species+RH, data = Survivorship)
data <- data$p.value
data <- data.frame(data)
library(writexl)
write_xlsx(data,"C:\\\\Users\\\\cmnieleb\\\\Documents\\\\Project Survivorship\\\\1.0\\\\Pairwise Survivorship Curves.xlsx")
#Plot survivorship curves with CIs
Survival <- ggsurvplot(surv.rh,data=Survivorship,conf.int = TRUE, legend.title = element_blank(),legend.labs = c("RH32","RH58","RH84","RH32","RH58",
                                                                                                                  "RH84","RH32","RH58","RH84"),
                         ggtheme = theme(
                          axis.text = element_text(size = 15),
                          axis.title = element_text(size = 20, face = "bold"),
                          strip.text.y = element_text(size = 20)),
                            palette = c("#F8766D","#00BA38","#619CFF","#F8766D","#00BA38","#619CFF","#F8766D","#00BA38","#619CFF")) +
  xlab("Day") +
  ylab("Survivorship")
Survival$plot + facet_grid(Species~.)


# Import questing height data
library(readxl)
QH <- read_excel("QuestingHeight.xlsx")

#Subset data by species
QHamblyomma <- subset(QH, Species == "Amblyomma")
QHdermacentor <- subset(QH, Species == "Dermacentor")
QHixodes <- subset(QH, Species == "Ixodes")

#Kruskal-Wallis and Mann-Whitney U-tests or Dunn's posthoc tests
library(sjstats)
library(coin)
library(dunn.test)
kruskal.test(QH~RH,data=QHamblyomma)
mwu(QHamblyomma,QH,RH)
dunn.test(QHamblyomma$QH,QHamblyomma$RH)
kruskal.test(QH~RH,data=QHdermacentor)
mwu(QHdermacentor,QH,RH)
dunn.test(QHdermacentor$QH,QHdermacentor$RH)
kruskal.test(QH~RH,data=QHixodes)
mwu(QHixodes,QH,RH)
dunn.test(QHixodes$QH,QHixodes$RH)

#Plot
library(ggplot2)
QHplot <- ggplot(QH, 
                 aes(y = QH, x = RH, fill = RH)) + 
  #geom_boxplot(size=1) +
  stat_summary(fun=mean,geom="errorbar",aes(ymax=..y..,ymin=..y..)) +
  geom_jitter(size = 1.5, aes(color = RH)) +
  ylab("Questing Height (cm)") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        legend.position = "none")
QHplot + facet_grid(~Species)


# Import questing frequency data
library(readxl)
QF <- read_excel("QuestingFrequency.xlsx")

#Subset data by species
QFamblyomma <- subset(QF, Species == "Amblyomma")
QFdermacentor <- subset(QF, Species == "Dermacentor")
QFixodes <- subset(QF, Species == "Ixodes")

#Kruskal-Wallis and Mann-Whitney U-tests or Dunn's posthoc tests
library(sjstats)
library(coin)
library(dunn.test)
kruskal.test(QF~RH,data=QFamblyomma)
mwu(QFamblyomma,QF,RH)
dunn.test(QFamblyomma$QF,QFamblyomma$RH)
kruskal.test(QF~RH,data=QFdermacentor)
mwu(QFdermacentor,QF,RH)
dunn.test(QFdermacentor$QF,QFdermacentor$RH)
kruskal.test(QF~RH,data=QFixodes)
mwu(QFixodes,QF,RH)
dunn.test(QFixodes$QF,QFixodes$RH)

#Plot
library(ggplot2)
QFplot <- ggplot(QF, 
                 aes(y = QF, x = RH, fill = RH)) + 
  #geom_boxplot(size=1) +
  stat_summary(fun=mean,geom="errorbar",aes(ymax=..y..,ymin=..y..)) +
  geom_jitter(size = 2, aes(color = RH)) +
  ylab("Questing Frequency") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        legend.position = "none")
QFplot + facet_grid(~Species)

#Linear survivorship curve
library(ggplot2)
field_surv_lineplot <- ggplot(QF, aes(x = Day, y = Survivorship * 100)) + 
  geom_line(size = 2, stat = "identity", aes(color = RH)) +
  geom_point(size = 4, aes(color = RH)) +
  #scale_x_continuous(breaks = pretty(0,30)) +
  xlab("Day") +
  ylab("Survivorship (%)") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18, face = "bold"),
        legend.title = element_blank(), 
        legend.text = element_text(size=18),
        legend.position = "top",
        strip.text.x = element_text(size = 20)) +
  scale_x_continuous(minor_breaks = seq(0, 30, 1), breaks = seq(0, 30, 5)) +
  scale_y_continuous(breaks = seq(0, 100, 20))
field_surv_lineplot + facet_grid(~Species,scales="free")
