### Article  ACNH, Coroller & Flinois


# Packages
library(tidyverse)
library(reshape2)
library(ggpubr)
library(wesanderson)
library(ggthemes)
library(ggpmisc)
library(devtools)
library(ade4)
library(corrplot)
library(ggrepel)
library(factoextra)
library(arules)
library(effects)
library(rstatix)

source_gist("524eade46135f6348140") #ajout à devtools pour formule et r2 sur ggplot2
getwd()
setwd("D:/Documents/ACNH article")

# Import data

dataACNH<- read_delim("Reponses.csv", ";", 
                      escape_double = FALSE, trim_ws = TRUE)
colnames(dataACNH)[colnames(dataACNH) == "Age (années) :"] <- "Age"
colnames(dataACNH[,4]) <- "Age"
resACNH <- dataACNH %>% select(3,13:54,61:74,55)
datACPasAC <- read_delim("QuestionACPasAC.csv", "\\t", escape_double = FALSE, trim_ws = TRUE)
datParticipant <- dataACNH %>% select(3:4, 7:12)
datjeux <- read_delim("jeux.tsv", "\\t", escape_double = FALSE, 
                      trim_ws = TRUE)
Joueurs_recent<- read_delim("Joueurs recent.csv", 
                          ";", escape_double = FALSE, trim_ws = TRUE)
ResACNHrecent <- Joueurs_recent %>% select(3,13:54,62:75,55:56)

# Explore data

summary(dataACNH$Age)
table(datjeux)

datjeux2 <- datjeux %>%
  mutate(game = strsplit(gsub("[][\\"]", "", game), ",")) %>%
  unnest(game)

table(datjeux2)

# Test Age

aggregate(x= dataACNH$Age,
          by= list(dataACNH$Joueur),
          FUN=mean)

levene_test(data = dataACNH, Age ~ Joueur)
wilcox.test(data = dataACNH, Age ~ Joueur )
boxplotage = 
  ggplot(dataACNH, aes(x=Joueur,y=Age, color = Joueur),
         position = "dodge")+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=1.4, notch=TRUE)+
  scale_color_manual(values=wes_palette(n=3, name="Cavalcanti1"))+
  theme_base()
boxplotage

# Test Naturalistic score

wilcox.test(data = dataACNH, Fibre_natura ~ Joueur )

boxplotfibre = 
  ggplot(dataACNH, aes(x=Joueur,y=Fibre_natura, color = Joueur),
         position = "dodge")+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=1.4, notch=TRUE)+
  scale_color_manual(values=wes_palette(n=3, name="Cavalcanti1"))+
  theme_base()+
  geom_hline( linetype="dashed",aes(yintercept = 10),color= "green", size=0.2)+
  geom_hline( linetype="dashed",aes(yintercept = 0),color= "red", size=0.2)+ 
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
boxplotfibre

dataACNH$Fibre_natura <- as.factor(dataACNH$Fibre_natura)
genact<- dataACNH %>% select(3,7:54,61:74,55,59)
meltgact <- melt(genact)
colnames(meltgact) <- c("ID","Peche","Ichtyo", "Poissonnerie", "Entomo", "Paleo","Fibre_natura","Joueur","Botanique","Question", "value")
merggact <- left_join(meltgact, datACPasAC, by="Question")
scoregact <- aggregate(merggact$value, by=list(Type=merggact$Type, Joueur=merggact$Joueur,ID=merggact$ID,AC_PasAC = merggact$AC_PasAC, Peche = merggact$Peche, Ichtyo = merggact$Ichtyo, Poissonnerie = merggact$Poissonnerie, Entomo = merggact$Entomo, Paleo = merggact$Paleo, Fibre_natura = merggact$Fibre_natura, Botanique = merggact$Botanique), FUN=mean)

scoregact$Fibre_natura <- as.numeric(factor(scoregact$Fibre_natura))
summary(lm(data = scoregact, x ~ Fibre_natura * Joueur))

regfibre = 
  ggplot(scoregact, aes (Fibre_natura, x, color = Joueur))+
  geom_jitter(height = 0.01)+
  geom_smooth(method = lm, se = FALSE)+
  scale_color_manual(labels = c("No", "Yes"),values=c("Non" = wes_palette("Zissou1")[1],"Oui"=wes_palette("Zissou1")[5]))+
  geom_hline( linetype="dashed",aes(yintercept = 1),color= "green", size=0.2)+
  geom_hline( linetype="dashed",aes(yintercept = 0),color= "red", size=0.2)+
  stat_smooth_func(geom="text",method="lm",hjust=0, vjust = c(2,0),parse=TRUE)+
  labs(color = "Animal Crossing \\n Player", title = "General scores", subtitle = "Function of Naturalist index", x = "Naturalist index", y = "Score")

regfibre

# boxplot N score ~ Player

scoregact$Fibre_natura <- as.character(factor(scoregact$Fibre_natura))
boxplotfibr = 
  ggplot(scoregact, aes(x = Fibre_natura ,y = x, color = Joueur), position = "dodge")+
  geom_boxplot()+
  scale_color_manual(values=wes_palette(n=2, name="Cavalcanti1"))+
  theme_base()+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
boxplotfibr

scoregact$Fibre_natura <- as.numeric(factor(scoregact$Fibre_natura))

# Hobbies factors

Peche <- as.factor(dataACNH$Peche)
levels(Peche)
Ichtyo <- as.factor(dataACNH$Ichtyo)
levels(Ichtyo)
Entomo <- as.factor(dataACNH$Entomo)
levels(Entomo)
Paleo <- as.factor(dataACNH$Paleo)
levels(Paleo)
Poissonnerie <- as.factor(dataACNH$Poissonnerie)
levels(Poissonnerie) 
Botanique <- as.factor(dataACNH$Botanique)
levels(Botanique)

act <- dataACNH %>% select(3,7:11,13:54,61:74,55,59)
meltact <- melt(act)
colnames(meltact) <- c("ID","Peche","Ichtyo", "Poissonnerie", "Entomo", "Paleo","Joueur","Botanique","Question", "value")
mergact <- left_join(meltact, datACPasAC, by="Question")
scoreact <- aggregate(mergact$value, by=list(Type=mergact$Type, Joueur=mergact$Joueur,ID=mergact$ID,AC_PasAC = mergact$AC_PasAC, Peche = mergact$Peche, Ichtyo = mergact$Ichtyo, Poissonnerie = mergact$Poissonnerie, Entomo = mergact$Entomo, Paleo = mergact$Paleo, Botanique = mergact$Botanique), FUN=mean)

levels(scoreact$Peche) <- c("Regular", "Occasion", "Non")
levels(scoreact$Ichtyo) <- c("Passion", "Oui aqua + poissons", "Oui poissons", "Occasion poissons", "Non+aqua", "Non")
levels(scoreact$Entomo) <- c("Passion", "Occasion", "Non")
levels(scoreact$Paleo) <- c("Passion", "Occasion", "Non")
levels(scoreact$Poissonnerie) <- c("Regular", "Occasion", "Arret", "Non")
levels(scoreact$Botanique) <- c("Passion", "Occasion+jardin", "Occasion-jardin", "Non")


# Fish score

scoreactpoi <- filter(scoreact, Type %in% "Poisson")

## Ichtiology score

scoreactpoi$Ichtyo <- as.numeric(factor(scoreactpoi$Ichtyo))

summary(lm(data = scoreactpoi, x ~ Ichtyo * Joueur))

regactpoi = 
  ggplot(scoreactpoi, aes (Ichtyo, x, color = Joueur))+
  geom_jitter(height = 0.01)+
  geom_smooth(method = lm, se = FALSE)+
  scale_color_manual(labels = c("No", "Yes"),values=c("Non" = wes_palette("Zissou1")[1],"Oui"=wes_palette("Zissou1")[5]))+
  geom_hline( linetype="dashed",aes(yintercept = 1),color= "green", size=0.2)+
  geom_hline( linetype="dashed",aes(yintercept = 0),color= "red", size=0.2)+
  stat_smooth(geom="text",method="lm",hjust=0,parse=TRUE)+
  labs(color = "Animal Crossing \\n Player", title = "Scores for Fishes", subtitle = "Function of Ichtyology/Aquariophilia interest", x = "Ichtyology/Aquariophilia Interest", y = "Score")

regactpoi

## Fishing score

scoreactpoi$Peche <- as.numeric(factor(scoreactpoi$Peche))

summary(lm(data = scoreactpoi, x ~ Peche * Joueur))

regactpoi2 = 
  ggplot(scoreactpoi, aes (Peche, x, color = Joueur))+
  geom_jitter(height = 0.01)+
  geom_smooth(method = lm, se = FALSE)+
  scale_color_manual(labels = c("No", "Yes"),values=c("Non" = wes_palette("Zissou1")[1],"Oui"=wes_palette("Zissou1")[5]))+
  geom_hline( linetype="dashed",aes(yintercept = 1),color= "green", size=0.2)+
  geom_hline( linetype="dashed",aes(yintercept = 0),color= "red", size=0.2)+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE)+
  labs(color = "Animal Crossing \\n Player", title = "Scores for Fishes", subtitle = "Function of Fishing practice frequency", x = "Fishing practice frequency", y = "Score")

regactpoi2

## Fish consumption

scoreactpoi$Poissonnerie <- as.numeric(factor(scoreactpoi$Poissonnerie))

summary(lm(data = scoreactpoi, x ~ Poissonnerie * Joueur))

regactpoi3 = 
  ggplot(scoreactpoi, aes (Poissonnerie, x, color = Joueur))+
  geom_jitter(height = 0.01)+
  geom_smooth(method = lm, se = FALSE)+
  scale_color_manual(labels = c("No", "Yes"),values=c("Non" = wes_palette("Zissou1")[1],"Oui"=wes_palette("Zissou1")[5]))+
  geom_hline( linetype="dashed",aes(yintercept = 1),color= "green", size=0.2)+
  geom_hline( linetype="dashed",aes(yintercept = 0),color= "red", size=0.2)+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE)+
  labs(color = "Animal Crossing \\n Player", title = "Scores for Fishes", subtitle = "Function of Fish consumption frequency", x = "Fish consumption frequency", y = "Score")

regactpoi3

# Insects score

scoreactins <- filter(scoreact, Type %in% "Insecte")

# Fossils score

scoreactpal <- filter(scoreact, Type %in% "Fossile")

## Entomology factor

scoreactins$Entomo <- as.numeric(factor(scoreactins$Entomo))

summary(lm(data = scoreactins, x ~ Entomo * Joueur))

regactins = 
  ggplot(scoreactins, aes (Entomo, x, color = Joueur))+
  geom_jitter(height = 0.01)+
  geom_smooth(method = lm, se = FALSE)+
  scale_color_manual(labels = c("No", "Yes"),values=c("Non" = wes_palette("Zissou1")[1],"Oui"=wes_palette("Zissou1")[5]))+
  geom_hline( linetype="dashed",aes(yintercept = 1),color= "green", size=0.2)+
  geom_hline( linetype="dashed",aes(yintercept = 0),color= "red", size=0.2)+
  stat_smooth_func(geom="text",method="lm",hjust=0, vjust = c(2,0),parse=TRUE)+
  labs(color = "Animal Crossing \\n Player", title = "Scores for Insects", subtitle = "Function of Entomology interest", x = "Entomology Interest", y = "Score")

regactins

#ANCOVA
## Prepare data

meltres <- melt(resACNH)
colnames(meltres) <- c("ID","Joueur","Question", "value")
options(max.print = 30000)

mergeres <- left_join(meltres, datACPasAC, by="Question")
scorecat <- aggregate(mergeres$value, by=list(Type=mergeres$Type,ID=mergeres$ID,AC_PasAC = mergeres$AC_PasAC, Joueur=mergeres$Joueur), FUN=sum)
scorecat2 <- aggregate(meltres$value, by=list(ID=meltres$ID, Joueur=meltres$Joueur), FUN=sum)
restot <- left_join(scorecat, datParticipant, by="ID")
restotAge <- left_join(scorecat2, datParticipant, by="ID")
colnames(restot)[colnames(restot) == "Age (années) :"] <- "Age"
attach(restot)

res.an1 <- aov(x ~ Fibre_natura * AC_PasAC * Joueur * Age, data = restot)
summary(res.an1)
summary.aov(res.an1)
summary.lm(res.an1)

glm1 <- glm(x ~ Fibre_natura + AC_PasAC + Joueur + Age + Joueur*Age + Fibre_natura*Joueur + AC_PasAC*Joueur + Age*Fibre_natura, data = restot)
plot(allEffects(glm1)) 

shapiro.test(res.an1$residuals) #Not Normal
library(car)
leveneTest(x ~ AC_PasAC * Joueur, data = restot
) # Heteroscedasticity

# 1. Homogeneity of variances
plot(res.an1, which = 1)

# 2. Normality
plot(res.an1, which = 2)

shapiro.test(x)
shapiro.test(Fibre_natura)

# Discretize Age
restot2 <- restot 
restot2$Age <- discretize(restot2$Age, method = "frequency")
res.an2 <- aov(x ~ Type * Age * AC_PasAC * Joueur + Fibre_natura, data = restot2)
summary(res.an2)
summary.lm(res.an2)



ggplot(restot, aes(Age , x ))+
  geom_point()

length(unique(restot[["Age"]]))

ggplot(restotAge, aes(Age , x/56 ))+
  geom_bin2d(bins = 180, binwidth = c(2, 0.04))+
  geom_hline( linetype="dashed",aes(yintercept = 1),color= "green", size=0.1)+
  geom_hline( linetype="dashed",aes(yintercept =0),color= "red", size=0.1)+
  ylab("Score")+
  theme_bw()
       
agenat <- table(restot2$Age, restot2$Fibre_natura)
prop.table(agenat,1)
aovagenat <- aov(Fibre_natura ~ Age, restot2)
summary(aovagenat)
agenat <- data.frame(Age = factor(restot2$Age), Nat = restot2$Fibre_natura)
tapply(agenat$Nat,agenat$Age,mean)
agenat <- data.frame(Age = factor(restot2$Age), score = restot2$x)
tapply(agenat$score,agenat$Age,mean)

shapiro.test(res.an2$residuals)

# 1. Homogeneity of variances
plot(res.an2, which = 1)

# 2. Normality
plot(res.an1, which = 2)

# ANOVA on Non Player for AC/Not AC organisms

restotNJ <- filter(restot2, Joueur %in% "Non")
res.an3 <- aov(x ~ AC_PasAC, data = restotNJ)
summary(res.an3)
summary.lm(res.an3)

shapiro.test(res.an3$residuals)
leveneTest(x ~ AC_PasAC, data = restotNJ)
shapiro.test(res.an2$residuals)

# 1. Homogeneity of variances
plot(res.an3, which = 1)

# 2. Normality
plot(res.an3, which = 2)

scorecat <- aggregate(mergeres$value, by=list(Type=mergeres$Type,Joueur=mergeres$Joueur,ID=mergeres$ID,AC_PasAC = mergeres$AC_PasAC), FUN=mean)


#ANOVA on Age for AC organisms

restotAC <- filter(restot2, AC_PasAC %in% "AC")

res.an4 <- aov(x ~ Age, data = restotAC)
summary(res.an4)
summary.lm(res.an4)

TukeyHSD(res.an4)
shapiro.test(res.an3$residuals)
leveneTest(x ~ AC_PasAC, data = restotNJ)
shapiro.test(res.an2$residuals)

# 1. Homogeneity of variances
plot(res.an3, which = 1)

# 2. Normality
plot(res.an3, which = 2)

scorecat <- aggregate(mergeres$value, by=list(Type=mergeres$Type,Joueur=mergeres$Joueur,ID=mergeres$ID,AC_PasAC = mergeres$AC_PasAC), FUN=mean)

?aov

### ANOVA on Recent vs Ex Players

meltresrecent <- melt(ResACNHrecent)
colnames(meltresrecent) <- c("ID","Joueur","Recent", "Question", "value")

mergeresrecent <- left_join(meltresrecent, datACPasAC, by="Question")
scorecatrecent <- aggregate(mergeresrecent$value, by=list(Type=mergeresrecent$Type,ID=mergeresrecent$ID,AC_PasAC = mergeresrecent$AC_PasAC, Joueur=mergeresrecent$Joueur, Recent=mergeresrecent$Recent), FUN=sum)
restotrecent <- left_join(scorecatrecent, datParticipant, by="ID")

res.an5 <- aov(x ~ Recent, data = restotrecent)
summary(res.an5)
summary.lm(res.an5)

ggplot(restotrecent, aes(y=x, x=Recent))+
  geom_boxplot()+
  geom_smooth(colour="red", method="lm", fill="red") +
  ylab("Score")+
  xlab("Recent Player or Ex Player")

# Boxplots

boxplot = 
  ggplot(scorecat, aes(Type,x, color = Joueur),
         position = "dodge")+
  geom_boxplot(aes(fill=scorecat$AC_PasAC),
               outlier.colour="black", outlier.shape=16,
               outlier.size=1.4, notch=TRUE)+
  scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+
  scale_color_manual(values=wes_palette(n=3, name="Cavalcanti1"))+
  theme_base()+
  geom_hline( linetype="dashed",aes(yintercept = 1),color= "green", size=0.2)+
  geom_hline( linetype="dashed",aes(yintercept = 0),color= "red", size=0.2)+ 
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
boxplot

boxplot2 = 
  ggplot(restot2, aes(Age,x),
         position = "dodge")+
  geom_boxplot()+
  xlab("Age (years)")+
  ylab("Score")+
  geom_text(x=1,y=9.2, label = "a" , family = "serif", color = "red")+
  geom_text(x=2,y=9.2, label = "b" , family = "serif", color = "darkgreen")+
  geom_text(x=3,y=9.2, label = "b" , family = "serif", color = "darkgreen")+
  theme_base()+
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
boxplot2

boxplot3 = 
  ggplot(scorecatrecent, aes(Recent,x, fill = Recent),
         position = "dodge")+
  geom_boxplot(
               outlier.colour="black", outlier.shape=16,
               outlier.size=1.4, notch=TRUE)+
  scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+
  scale_color_manual(values=wes_palette(n=3, name="Cavalcanti1"))+
  theme_base()+
  geom_hline( linetype="dashed",aes(yintercept = 1),color= "green", size=0.2)+
  geom_hline( linetype="dashed",aes(yintercept = 0),color= "red", size=0.2)+ 
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
boxplot3

### Iterations

n <- 5000
scorecat <- as.data.frame(scorecat)
sampl500 <- scorecat %>% group_by(AC_PasAC) %>% sample_n(n, replace = TRUE)

ggplot(sampl500, aes(Type,x, color = AC_PasAC),
       position = "dodge")+
  geom_point(aes(shape=Joueur))+
  geom_jitter()+
  scale_color_manual(values=wes_palette(n=3, name="Cavalcanti1"))+
  theme_base()+
  geom_hline( linetype="dashed",aes(yintercept = 1),color= "green", size=0.2)+
  geom_hline( linetype="dashed",aes(yintercept = 0),color= "red", size=0.2)+ 
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

boxplot2 = 
  ggplot(sampl500, aes(Type,x, color = AC_PasAC),
         position = "dodge")+
  geom_boxplot(aes(fill=Joueur),
               outlier.colour="black", outlier.shape=16,
               outlier.size=1.4, notch=TRUE)+
  scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+
  scale_color_manual(values=wes_palette(n=3, name="Cavalcanti1"))+
  theme_base()+
  geom_hline( linetype="dashed",aes(yintercept = 1),color= "green", size=0.2)+
  geom_hline( linetype="dashed",aes(yintercept = 0),color= "red", size=0.2)+ 
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
boxplot2

boxplotgen = 
  ggplot(sampl500, aes(Joueur,x),
         position = "dodge")+
  geom_boxplot(aes(fill=AC_PasAC),
               outlier.colour="black", outlier.shape=16,
               outlier.size=1.4, notch=TRUE)+
  scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1"))+
  theme_base()+
  geom_hline( linetype="dashed",aes(yintercept = 1),color= "green", size=0.2)+
  geom_hline( linetype="dashed",aes(yintercept = 0),color= "red", size=0.2)+ 
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
boxplotgen

wilcox.test(data = sampl500, x ~ Joueur )

wilcox.test(data = sampl500, x ~ AC_PasAC )

Simon <- aov(x ~ Joueur * AC_PasAC, data = sampl500)
summary(Simon)

count(scorecat, Joueur)
count(sampl500, Joueur)

taborganAC= filter(sampl500, AC_PasAC=="AC")
taborganNAC= filter(sampl500, AC_PasAC=="PasAC")

taborganAC # With AC organisms
taborganNAC # With Not AC organisms

JAC = filter(taborganAC, Joueur=="Oui") # Players and results AC
NJAC = filter(taborganAC, Joueur=="Non") # Non players and results AC

JNAC = filter(taborganNAC, Joueur=="Oui") # Players and results NAC
NJNAC = filter(taborganNAC, Joueur=="Non") # Non Players and results NAC

JAC
NJAC
JNAC
NJNAC

# Fish
JPoisAC= filter(JAC, Type=="Poisson") # J  AC
NJPoisAC= filter(NJAC, Type=="Poisson")# NJ   AC
JPoisNAC= filter(JNAC, Type=="Poisson")# J NAC
NJPoisNAC= filter(NJNAC, Type=="Poisson")# NJ  NAC

XJPoisAC = JPoisAC$x
XNJPoisAC = NJPoisAC$x
XJPoisNAC = JPoisNAC$x
XNJPoisNAC = NJPoisNAC$x

wilcox.test(XJPoisAC,XNJPoisAC)
wilcox.test(XJPoisNAC,XNJPoisNAC)

# Plants
JFleurAC= filter(JAC, Type=="Fleur") # J  AC
NJFleurAC= filter(NJAC, Type=="Fleur")# NJ   AC
JFleurNAC= filter(JNAC, Type=="Fleur")# J NAC
NJFleurNAC= filter(NJNAC, Type=="Fleur")# NJ  NAC

XJFleurAC = JFleurAC$x
XNJFleurAC = NJFleurAC$x
XJFleurNAC = JFleurNAC$x
XNJFleurNAC = NJFleurNAC$x

wilcox.test(XJFleurAC,XNJFleurAC)
wilcox.test(XJFleurNAC,XNJFleurNAC)



# Insects
JInsecteAC= filter(JAC, Type=="Insecte") # J  AC
NJInsecteAC= filter(NJAC, Type=="Insecte")# NJ   AC
JInsecteNAC= filter(JNAC, Type=="Insecte")# J NAC
NJInsecteNAC= filter(NJNAC, Type=="Insecte")# NJ  NAC

XJInsecteAC = JInsecteAC$x
XNJInsecteAC = NJInsecteAC$x
XJInsecteNAC = JInsecteNAC$x
XNJInsecteNAC = NJInsecteNAC$x

wilcox.test(XJInsecteAC,XNJInsecteAC)
wilcox.test(XJInsecteNAC,XNJInsecteNAC)


# Fossils
JFossileAC= filter(JAC, Type=="Fossile") # J  AC
NJFossileAC= filter(NJAC, Type=="Fossile")# NJ   AC
JFossileNAC= filter(JNAC, Type=="Fossile")# J NAC
NJFossileNAC= filter(NJNAC, Type=="Fossile")# NJ  NAC

XJFossileAC = JFossileAC$x
XNJFossileAC = NJFossileAC$x
XJFossileNAC = JFossileNAC$x
XNJFossileNAC = NJFossileNAC$x

wilcox.test(XJFossileAC,XNJFossileAC)
wilcox.test(XJFossileNAC,XNJFossileNAC)

boxplot2 = 
  
  ggplot(sampl500, aes(x=AC_PasAC, y=x, color = AC_PasAC),
         position = "dodge")+
  geom_violin(aes(fill=Joueur), position="dodge",
              outlier.colour="black", outlier.shape=16,
              outlier.size=1, notch=TRUE)+
  scale_color_manual(values=c("black","red"))+
  geom_jitter(color="black", size=0.01, alpha=0.05) +
  theme_classic()+
  scale_fill_brewer(palette = "Accent")+
  facet_wrap(~factor(Type), ncol = 2)+
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill=Joueur,alpha=0.1), position="dodge")+
  
  stat_compare_means(aes(group=Joueur),method = "wilcox", label.y = 1.1)+
  stat_compare_means(aes(group= AC_PasAC),method = "wilcox", label.y = 1.2, label.x = 1.3)+ #compare ac pas ac joueur ? ou joueur tout court ?
  
  geom_hline( linetype="dashed",aes(yintercept = 1),color= "green", size=0.1)+
  geom_hline( linetype="dashed",aes(yintercept =0),color= "darkred", size=0.1)+ 
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
boxplot2



boxplottotale = 
  
  ggplot(sampl500, aes(x=AC_PasAC, y=x, color = AC_PasAC),
         position = "dodge")+
  geom_violin(aes(fill=Joueur), position= "dodge",
              outlier.colour="black", outlier.shape=16,
              outlier.size=1, notch=TRUE)+
  scale_color_manual(values=c("black","red"))+
  geom_jitter(color="black", size=0.01, alpha=0.05) +
  theme_classic()+
  scale_fill_brewer(palette = "Accent")+
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill=Joueur,alpha=0.1), position="dodge")+
  
  stat_compare_means(aes(group=Joueur),method = "wilcox", label.y = 1.1)+
  stat_compare_means(aes(group= AC_PasAC),method = "wilcox", label.y = 1.2, label.x = 1.3)+ #compare ac pas ac joueur ? ou joueur tout court ?
  
  geom_hline( linetype="dashed",aes(yintercept = 1),color= "green", size=0.1)+
  geom_hline( linetype="dashed",aes(yintercept =0),color= "darkred", size=0.1)+ 
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
boxplottotale

summary(sampl500)
count(sampl500$AC_PasAC)


#??? Tests iterations

n <- 5000
n_sim <- 10000
coef <- matrix(ncol = 2, nrow = n_sim)
i = 0
#for (i in 1:n_sim) {
  repeat {
    i = i+1
  sampl500 <- scorecat %>% group_by(AC_PasAC) %>% 
    sample_n(n, replace = TRUE)
  taborganAC <- filter(sampl500, AC_PasAC=="AC")
  taborganNAC <- filter(sampl500, AC_PasAC=="PasAC")
  
  taborganAC #contient les organismes d'AC
  taborganNAC #contient les organismes NAC
  
  JAC <- filter(taborganAC, Joueur=="Oui") # Joueurs et résultats AC
  NJAC <- filter(taborganAC, Joueur=="Non") # Non Joueurs et résultats AC
  
  JNAC <- filter(taborganNAC, Joueur=="Oui") # Joueurs et résultats NAC
  NJNAC <- filter(taborganNAC, Joueur=="Non") # Non Joueurs et résultats NAC
  
  JFleurAC <- filter(JAC, Type=="Fleur") # J  AC
  NJFleurAC <- filter(NJAC, Type=="Fleur")# NJ   AC
  JFleurNAC <- filter(JNAC, Type=="Fleur")# J NAC
  NJFleurNAC <- filter(NJNAC, Type=="Fleur")# NJ  NAC
  
  XJFleurAC <- JFleurAC$x
  XNJFleurAC <- NJFleurAC$x
  XJFleurNAC <- JFleurNAC$x
  XNJFleurNAC <- NJFleurNAC$x
  coef[i,] <- wilcox.test(XJFleurAC,XNJFleurAC, exact=FALSE )$p.value
  
  if (i==n_sim) {break}
  }
#}

for(i in 1:nrow(coef)){
  if(coef[,1][i]<0.05){
    coef[,2][i]=1
  }
}
mean(coef[,2])

### environ 55% des tests sont significatifs: on ne peut pas conclure sur les fleurs AC

### 93% de tests significatifs pour les Poissons non AC

###??? 91,4% de tests significatifs pour les fossiles Non AC

### 99,9% de tests significatifs pour les insectes Non AC

#??? Iterations tests insects

n2 <- 500
n_sim2 <- 1000
coef <- matrix(ncol = 2, nrow = n_sim2)
i = 0
#for (i in 1:n_sim) {
repeat {
  i = i+1
  sampl500 <- scorecat %>% group_by(AC_PasAC) %>% 
    sample_n(n2, replace = TRUE)
  taborganAC <- filter(sampl500, AC_PasAC=="AC")
  taborganNAC <- filter(sampl500, AC_PasAC=="PasAC")
  
  taborganAC #contient les organismes d'AC
  taborganNAC #contient les organismes NAC
  
  tabinsAC <- filter(taborganAC, Type=="Insecte") 
  tabinsNAC <- filter(taborganNAC, Type=="Insecte")

  
  XInsAC <- tabinsAC$x
  XInsNAC <- tabinsNAC$x
 
  coef[i,] <- wilcox.test(XInsAC,XInsNAC, exact=FALSE )$p.value
  
  if (i==n_sim2) {break}
}
#}

for(i in 1:nrow(coef)){
  if(coef[,1][i]<0.05){
    coef[,2][i]=1
  }
}
mean(coef[,2])

# ACP

### Prepare data

ACPtab <- dcast(scorecat, Joueur + AC_PasAC + ID ~ Type , value.var = 'x')

Entomo <- select(scoreactins, ID, Entomo, Joueur, AC_PasAC)
ACPtab <- left_join(x = ACPtab , y = scoreactpoi , by = c("ID", "Joueur", "AC_PasAC") ) %>% 
  select(one_of(c("Joueur" , "AC_PasAC", "ID" ,  "Peche", "Ichtyo", "Poissonnerie", "Poisson", "Insecte", "Fleur", "Fossile")))

ACPtab <- left_join(x = ACPtab, y = Entomo, by = c("ID", "Joueur", "AC_PasAC"),  all.y = FALSE)

scoreactpal$Paleo <- as.numeric(factor(scoreactpal$Paleo))
Paleo <- select(scoreactpal, ID, Paleo, Joueur, AC_PasAC)

ACPtab <- left_join(x = ACPtab, y = Paleo, by = c("ID", "Joueur", "AC_PasAC"),  all.y = FALSE)

scoreactfle <- filter(scoreact, Type %in% "Fleur")
scoreactfle$Botanique <- as.numeric(factor(scoreactfle$Botanique))
Botanique <- select(scoreactfle, ID, Botanique, Joueur, AC_PasAC)

ACPtab <- left_join(x = ACPtab, y = Botanique, by = c("ID", "Joueur", "AC_PasAC"),  all.y = FALSE)

Agedat <- select(dataACNH, Age, ID)

ACPtab <- left_join(x = ACPtab, y = Agedat, by = "ID",  all.y = FALSE)

which(is.na(ACPtab), arr.ind=TRUE) # NA sur ID 042

ACPtab <- na.omit(ACPtab)

### ACP

pca <-ACPtab %>%
  select(Peche, Ichtyo, Poissonnerie, Insecte, Fleur, Fossile, Poisson, Entomo, Paleo, Botanique, Age ) %>% 
  dudi.pca(scannf = FALSE, nf= ncol(.)) 
pca 

### Inertias

res <- inertia.dudi(pca, col.inertia = TRUE)
res 

### Proper values

inertia <- res$tot.inertia 
round(inertia,3)

### Absolute contributions 

cont.abs <- res$col.abs
round(cont.abs,2)

### Correlation plot

correl <- ACPtab %>%
  select(Peche, Ichtyo, Poissonnerie, Insecte, Fleur, Fossile, Poisson, Entomo, Paleo, Botanique, Age ) %>% 
  cor()
corrplot(correl, type="upper", order="hclust", tl.col="black", tl.srt=45)

### Individuals plot

indiv <- ACPtab %>%
  select(Joueur, AC_PasAC) %>%
  bind_cols(pca$li)

indiv %>%
  ggplot(aes(x = Axis1, y = Axis2, color=Joueur,shape=AC_PasAC)) +
  geom_point(size=3) +
  labs( title = "Graphique des individus de l'ACP")+
  theme_bw()

fviz_pca_ind(pca, label="none", habillage= ACPtab$Joueur,
             addEllipses=TRUE, ellipse.type = "confidence", palette = "Dark2")

#Correlation plot

s.corcircle(pca$co, xax = 1, yax = 2, sub = "Cercle des corrélations des axes 1 et 2 de l'ACP")


fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

## ACP 2

ACPtab2 <- ACPtab
ACPtab2$Joueur <- as.numeric(factor(ACPtab2$Joueur))
pca2 <-ACPtab2 %>%
  select(Joueur, Ichtyo, Insecte, Fleur, Fossile, Poisson, Entomo, Paleo, Botanique, Age ) %>% 
  dudi.pca(scannf = FALSE, nf= ncol(.)) 
pca2 

### Inertias 2

res2 <- inertia.dudi(pca2, col.inertia = TRUE)
res2 

### Proper values 2

inertia2 <- res2$tot.inertia 
round(inertia2,3)

### Absolute contributions 2

cont.abs2 <- res2$col.abs
round(cont.abs2,2)

### Correlation plot 2

correl2 <- ACPtab2 %>%
  select(Joueur, Ichtyo, Insecte, Fleur, Fossile, Poisson, Entomo, Paleo, Botanique, Age ) %>% 
  cor()
corrplot(correl2, type="upper", order="hclust", tl.col="black", tl.srt=45)

### Individual plot 2

fviz_pca_ind(pca2, label="none", habillage= ACPtab$AC_PasAC,
             addEllipses=TRUE, ellipse.type = "confidence", palette = "Dark2")

#Coreelqtion circle 2

fviz_pca_var(pca2,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

## ACP 3

ACPtab3 <- ACPtab2
Fibnat <- select(dataACNH, Fibre_natura, ID)
Fibnat$Fibre_natura <- as.numeric(as.character(Fibnat$Fibre_natura))

ACPtab3 <- left_join(x = ACPtab3, y = Fibnat, by = "ID",  all.y = FALSE)

pca3 <-ACPtab3 %>%
  select(Joueur, Fibre_natura, Ichtyo, Insecte, Fleur, Fossile, Poisson, Entomo, Paleo, Botanique, Age ) %>% 
  dudi.pca(scannf = FALSE, nf= ncol(.)) 
pca3

summary(ACPtab3)

### Inertias 3

res3 <- inertia.dudi(pca3, col.inertia = TRUE)
res3

### Proper values 3

inertia3 <- res3$tot.inertia 
round(inertia3,3)

### Absolute contributions 3

cont.abs3 <- res3$col.abs
round(cont.abs3,2)

### Correlation plot 3

correl3 <- ACPtab3 %>%
  select(Joueur, Fibre_natura, Ichtyo, Insecte, Fleur, Fossile, Poisson, Entomo, Paleo, Botanique, Age ) %>% 
  cor()
corrplot(correl3, type="upper", order="hclust", tl.col="black", tl.srt=45)

### Individual plot 3

fviz_pca_ind(pca3, label="none", habillage= ACPtab$AC_PasAC,
             addEllipses=TRUE, ellipse.type = "confidence", palette = "Dark2")

#Correlation plot 3

fviz_pca_var(pca3,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             labels = c("Joueur" = "Player", "Fossile" = "Fossil"),
             labelsize = 3,
             repel = TRUE  )   # Avoid text overlapping 


## ACP 4

ACPtab4 <- ACPtab3


ACPtab4 <- filter(ACPtab4, AC_PasAC %in% "AC")

pca4 <-ACPtab4 %>%
  select(Joueur, Fibre_natura, Ichtyo, Insecte, Fleur, Fossile, Poisson, Entomo, Paleo, Botanique, Age ) %>% 
  dudi.pca(scannf = FALSE, nf= ncol(.)) 
pca4

summary(ACPtab4)

### Inertias 4

res4 <- inertia.dudi(pca4, col.inertia = TRUE)
res4

### Proper values 4

inertia4 <- res4$tot.inertia 
round(inertia4,3)

### Absolute contributions 4

cont.abs4 <- res4$col.abs
round(cont.abs4,2)

### Correlation plot 4

correl4 <- ACPtab4 %>%
  select(Joueur, Fibre_natura, Ichtyo, Insecte, Fleur, Fossile, Poisson, Entomo, Paleo, Botanique, Age ) %>% 
  cor()
corrplot(correl4, type="upper", order="hclust", tl.col="black", tl.srt=45)

### Individuals plot 4

fviz_pca_ind(pca4, label="none", habillage= ACPtab$Joueur,
             addEllipses=TRUE, ellipse.type = "confidence", palette = "Dark2")

#Correlation circle 4

fviz_pca_var(pca4,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE  )   # Avoid text overlapping 

#Correlation circle 4 (pour Sim)

fviz_pca_var(pca4,
             col.var = "contrib",geom= c("arrow","point","text"), 
             gradient.cols = c("#F7230C","#1C100B" ), pointsize= "contrib",
             repel = TRUE  )+
  theme_pubr()


fviz_pca_var(pca4,
             col.var = "contrib",geom= c("arrow","point"), 
             gradient.cols = c("#73C2FB","#F7230C" ), pointsize= "contrib",
             repel = TRUE  )+
  theme_pubr()

