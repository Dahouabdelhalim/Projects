library(dplyr)
library(ggplot2)

setwd("~/Documents/")
dat <- read.csv(file = "dat_with_guilds.csv", header = T, stringsAsFactors = F)

# Data set split by Diet.5Cat

#deleting records where distance = 50 and = 100
dat<-dat %>% filter(Distance<50)

unique(dat$Diet.5Cat)


library(tidyr)
library(gridExtra)
dat2 <- tidyr::complete(dat, Replicate, Stratum, Diet.5Cat)
dat2$Number[is.na(dat2$Number)] <- 0

# Average abundance *per point* in *each stratum*, for birds of each *diet group*

step1 <- dat2 %>% group_by(Replicate, Stratum, Diet.5Cat) %>% summarise(
  tot.abd = sum(Number),
  tot.spr = n_distinct(Species, na.rm=T))

step2 <- step1 %>% group_by(Stratum, Diet.5Cat) %>%  summarise(
  av.abd =  mean(tot.abd),
  se.abd = sd(tot.abd)/sqrt(length(tot.abd)),
  av.spr =  mean(tot.spr),
  se.spr = sd(tot.spr)/sqrt(length(tot.spr)),
  ss = length(tot.abd))

# Abundance

stratnames<-c("Con. Tea", "Org. Tea", "Shade Tea", "Fragment", "Rainforest")

guilds <- c("Frugivore","Insectivore", "Omnivore", "Granivore","Raptors")
names(guilds) <- unique(step2$Diet.5Cat)


p1 <- ggplot(step2, aes(x=factor(Stratum, labels = stratnames), y = av.abd))+
  geom_bar(stat = "identity", fill = "gray", colour = "black", width = 0.6)+
  geom_errorbar(data = step2, aes( ymin = av.abd - 2*se.abd,
                                   ymax = av.abd + 2*se.abd,
                                   width = 0.1))+
  facet_wrap(.~Diet.5Cat, ncol = 1, labeller= labeller( Diet.5Cat = guilds), scales = "free")+
  theme_classic() +
  theme(strip.text = element_text(size=11),
        axis.text.x = element_text(size = 10, angle = 0, colour = "black"),
        axis.text.y = element_text(size = 10, angle = 0, colour = "black"),  
        axis.title.x = element_text(size = 10, angle = 0, colour = "black"),
        axis.title.y = element_text(size = 10, angle = 90, colour = "black"))+
  ylab("Abundance (birds/point)")+xlab("")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.5)))


# Richness
p2 <- ggplot(step2, aes(x=factor(Stratum, labels = stratnames), y = av.spr))+
  geom_bar(stat = "identity", fill = "gray", colour = "black", width = 0.6)+
  geom_errorbar(data = step2, aes( ymin = av.spr - 2*se.spr,ymax = av.spr + 2*se.spr,
                                   width = 0.1))+
  facet_wrap(.~Diet.5Cat, ncol = 1, labeller= labeller( Diet.5Cat = guilds), scales = "free")+
  theme_classic() +
  theme(strip.text = element_text(size=11),
        axis.text.x = element_text(size = 10, angle = 0, colour = "black"),
        axis.text.y = element_text(size = 10, angle = 0, colour = "black"),  
        axis.title.x = element_text(size = 10, angle = 0, colour = "black"),
        axis.title.y = element_text(size = 10, angle = 90, colour = "black"))+
  ylab("Species richess (species/point)")+xlab("")+
  scale_y_continuous(expand = expansion(mult = c(0.000002, 0.5)))

px<-grid.arrange(p2,p1, nrow = 1)
ggsave("Fig5.jpg", px, width=8, height=11)

#### Fruit Nectar birds ####

library(multcomp) ####

fruinect.spr <- step1 %>% group_by(Stratum) %>% filter(Diet.5Cat=="FruiNect") 
fruinect.spr$Stratum <- as.factor(fruinect.spr$Stratum)

# Richness

glmod.fruinect <-glm(data=fruinect.spr, tot.spr~Stratum, family = poisson)
summary(glmod.fruinect)
Tukeyres1<-glht(glmod.fruinect, mcp(Stratum="Tukey"))
Tukeyres1
plot(Tukeyres1)
cld(Tukeyres1) #Compact letter display

# Abundance
glmod.fruinect <-glm(data=fruinect.spr, tot.abd~Stratum, family = poisson)
summary(glmod.fruinect)
Tukeyres1<-glht(glmod.fruinect, mcp(Stratum="Tukey"))
Tukeyres1
plot(Tukeyres1)
cld(Tukeyres1) #Compact letter display

#### Invertebrate birds ####

invt <- step1 %>% group_by(Stratum) %>% filter(Diet.5Cat=="Invertebrate") 
invt$Stratum <- as.factor(invt$Stratum)

# Richness

glmod.invt <-glm(data=invt, tot.spr~Stratum, family = poisson)
summary(glmod.invt)
Tukeyres1<-glht(glmod.invt, mcp(Stratum="Tukey"))
Tukeyres1
plot(Tukeyres1)
cld(Tukeyres1) #Compact letter display

# Abundance
glmod.invt<-glm(data=invt, tot.abd~Stratum, family = poisson)
summary(glmod.invt)
Tukeyres1<-glht(glmod.invt, mcp(Stratum="Tukey"))
Tukeyres1
plot(Tukeyres1)
cld(Tukeyres1) #Compact letter display


#### Omnivorous birds ####

omni <- step1 %>% group_by(Stratum) %>% filter(Diet.5Cat=="Omnivore") 
omni$Stratum <- as.factor(omni$Stratum)

# Richness

glmod.omni <-glm(data=omni, tot.spr~Stratum, family = poisson)
summary(glmod.omni)
Tukeyres1<-glht(glmod.omni, mcp(Stratum="Tukey"))
Tukeyres1
plot(Tukeyres1)
cld(Tukeyres1) #Compact letter display

# Abundance
glmod.omni<-glm(data=omni, tot.abd~Stratum, family = poisson)
summary(glmod.omni)
Tukeyres1<-glht(glmod.omni, mcp(Stratum="Tukey"))
Tukeyres1
plot(Tukeyres1)
cld(Tukeyres1) #Compact letter display


#### Granivorous birds ####

gran <- step1 %>% group_by(Stratum) %>% filter(Diet.5Cat=="PlantSeed") 
gran$Stratum <- as.factor(gran$Stratum)

# Richness

glmod.gran <-glm(data=gran, tot.spr~Stratum, family = poisson)
summary(glmod.gran)
Tukeyres1<-glht(glmod.gran, mcp(Stratum="Tukey"))
Tukeyres1
plot(Tukeyres1)
cld(Tukeyres1) #Compact letter display

# Abundance
glmod.gran<-glm(data=gran, tot.abd~Stratum, family = poisson)
summary(glmod.gran)
Tukeyres1<-glht(glmod.gran, mcp(Stratum="Tukey"))
Tukeyres1
plot(Tukeyres1)
cld(Tukeyres1) #Compact letter display

#### Raptors ####

rapt <- step1 %>% group_by(Stratum) %>% filter(Diet.5Cat=="VertFishScav") 
rapt$Stratum <- as.factor(rapt$Stratum)

# Richness

glmod.rapt <-glm(data=rapt, tot.spr~Stratum, family = poisson)
summary(glmod.rapt)
Tukeyres1<-glht(glmod.rapt, mcp(Stratum="Tukey"))
Tukeyres1
plot(Tukeyres1)
cld(Tukeyres1) #Compact letter display

# Abundance
glmod.rapt<-glm(data=rapt, tot.abd~Stratum, family = poisson)
summary(glmod.rapt)
Tukeyres1<-glht(glmod.rapt, mcp(Stratum="Tukey"))
Tukeyres1
plot(Tukeyres1)
cld(Tukeyres1) #Compact letter display
