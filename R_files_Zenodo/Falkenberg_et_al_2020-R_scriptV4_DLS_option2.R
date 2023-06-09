#LOADING PACKAGES

library(ggplot2)
library(nlme)
library(gridExtra)
library(grid)
library(egg)
library(gtable)
library(lsmeans)
library(vegan)


#IMPORTING DATA

O2<- read.csv("Falkenberg_et_al_2020-O2_Data.csv")
snails<- read.csv("Falkenberg_et_al_2020-Snail&Feeding_Data.csv")
ibutton.data<-read.csv("Falkenberg_et_al_2020-IButtonData.csv")



####Temperature Data####

#linear model to establish treatment differences 
temp.model<-lme(Value~treatment, random = ~ 1|Tank,data=ibutton.data)
summary(temp.model)
anova(temp.model)

plot(temp.model)
plot(residuals(temp.model))
qqnorm(residuals(temp.model))

#calculating mean temperature and SE
mean<-tapply(ibutton.data$Value,ibutton.data$treatment, mean )
n<-tapply(ibutton.data$Value,ibutton.data$treatment, length)
SD<-tapply(ibutton.data$Value,ibutton.data$treatment, sd)
SE<-SD/sqrt(n)

mean
SE

#######Oxygen Consumption######

#recording which round each snail is in so that it is paired with the correct control later
O2$Round<-as.numeric(O2$Round)

for(i in 1:length(snails$basket_no)) {
  a<-O2$Round[which(O2$SnailID==snails$basket_no[i])]
  snails$round[i]<-a[1]
  }

#Calculating  rate of O2 consumption using a linear model

str(O2)

O2$SnailID<-as.factor((O2$SnailID))

Names<-levels(O2$SnailID)
new<-as.data.frame(Names)

for(i in 1:length(new$Names)) {
  model<-lm(Value~delta_t,   
            data=O2[which(O2$SnailID==new$Names[i]),])
  new$rate[i]<- model$coefficients[2]
  new$intercept[i]<- model$coefficients[1]
}
  
coef<-new$rate

length(coef)

Rate<-as.data.frame(coef)
Rate$SnailID<-new$Names

snails$O2rate<-NA



#adding a raw O2 Rate column to our information on the snails

for(i in 1:length(Rate$SnailID)){
  snails$O2rate[which(snails$basket_no==Rate$SnailID[i])]<- -Rate$coef[i]
}

str(snails)

snails$tank.n<-as.factor(snails$tank.n)

#removing any snails that have died before Aug21 and reducing the number of controls, this way there are the same number of controls as Lunella or chlorostoma in each tank.

snails<-snails[which(snails$O2rate!="NA"),]

#setting aside data for feeding analysis. will return to this dataframe later
feeding<-snails

#back to oxygen

#Correcting for water column O2 consumption (blank controls) and water volume
control<-snails[which(snails$species=="Control"),]
snails<-snails[which(snails$species!="Control"),]

#this formula comes from Archimedes Principle and 1.0236 is the density of sea water
snails$SnailV<-(snails$wet_weight-snails$sub_weight)/1.0236 
snails$WaterV<-(55-snails$SnailV)/1000

str(snails)

snails$O2Rates2<-snails$O2rate*snails$WaterV

str(control)
control$O2Rates2<-control$O2rate*(55/1000)

for(i in 1:length(snails$basket_no)){
  a<-snails$tank.n[i]
  b<-snails$round[i]
  c<-which(control$tank.n == a & control$round==b)
  d<-control$O2Rates2[c]
  snails$O2Rate[i]<- snails$O2Rates2[i]-d
}



snails$final<-100*snails$O2Rate/snails$wet_weight
snails$species<-as.character(snails$species)
snails$species<-as.factor(snails$species)

snails$O2Consumption<-snails$final

#Prepping O2 Data for Primer####

hist(snails$O2Consumption)

str(snails)

ReadyForPrimer<-snails[,c(16,4,3,1,2)]

str(ReadyForPrimer)

write.csv(ReadyForPrimer, "ReadyForPrimer_O2Consumption.csv")


#FIGURE O2 consumption w/Mean and SE#####

treatment<-c("Ambient","Ambient","Elevated","Elevated")
O2.Summary<-as.data.frame(treatment)
O2.Summary$species<-c("Chlorostoma","Lunella","Chlorostoma","Lunella")

O2.Summary$mean<-tapply(snails$final,snails$treatment:snails$species, mean )
O2.Summary$n<-tapply(snails$final,snails$treatment:snails$species, length)
O2.Summary$SD<-tapply(snails$final,snails$treatment:snails$species, sd)
O2.Summary$SE<-O2.Summary$SD/sqrt(O2.Summary$n)
O2.Summary$median<-tapply(snails$final,snails$treatment:snails$species, median)

plotC.mean.SE<-
  ggplot(data = O2.Summary, aes(y = mean,x =treatment, color=treatment))+
  scale_color_manual(values=c( "#80cdc1","#dfc27d" ), 
                    name="Treatment",
                    breaks=c("Ambient", "Elevated"),
                    labels=c("Ambient", "Elevated"))+
  geom_pointrange(ymin=O2.Summary$mean-O2.Summary$SE, 
                  ymax=O2.Summary$mean+O2.Summary$SE,
                  size=2.5)+  
  theme_classic(base_size = 20)+
  ylab(bquote('Oxygen Consumption( '*mu~g~O[2]~ g^-1~min^-1*')'))+
  xlab("Snail species")+
  theme(legend.position="none",
        strip.text.x =element_text(face = "italic"),
        axis.title.x = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank())+
  facet_grid(.~species, switch = "x")+
  ylim(.015,.055)
  #annotate("text", x = .6, y = 0.055, label = c("C.", ""))

dat_text <- data.frame(
  label = c("C.", " "),
  species   = c("Chlorostoma", "Lunella"))

plotC.mean.SE<- plotC.mean.SE + geom_text(
        data    = dat_text,
        mapping = aes(x = 0.45, y = 0.054, label = label),
        colour = c( "Black","Black"),
        size = 6,
        hjust   = -0.1,
        vjust   = -1)

plotC.mean.SE




#####Survival#####

#Survival August###
#extracting the number of snails of each species alive in each tank on Aug 21

alive<-table(snails$tank.n, 
             snails$species, useNA = "ifany")
alive<-as.data.frame(alive)
levels(alive$Var2)

#removing control baskets as there are no snails in them
alive<-alive[which(alive$Var2!="Control"),]

#re-associating each tank with the appropriate treatment
for(i in 1:length(alive$Var1)){
  
  if (alive$Var1[i]  == "1") {
    alive$treatment[i]= "Elevated"
  } else if (alive$Var1[i]  == "4") {
    alive$treatment[i]= "Elevated"
  } else if (alive$Var1[i]  == "5") {
    alive$treatment[i]= "Elevated"  
  } else if (alive$Var1[i]  == "8") {
    alive$treatment[i]= "Elevated" 
  } else if (alive$Var1[i]  == "9") {
    snails$treatment[i]= "Elevated"  
  } else if (alive$Var1[i]  == "12") {
    alive$treatment[i]= "Elevated"
  } else {
    alive$treatment[i]= "Ambient"
  }}

alive$Proportion<-alive$Freq/3


#FIGURE Alive w/Mean and SE#####

# recall: Var2 == Species; Var1 == Tank



alive$final<-alive$Proportion*100
alive$treatment<-as.factor(alive$treatment)
alive$species<-as.factor(alive$Var2)
alive$tank<-as.factor(alive$Var1)

treatment<-c("Ambient","Ambient","Elevated","Elevated")
Alive.Summary<-as.data.frame(treatment)
Alive.Summary$species<-c("Chlorostoma","Lunella","Chlorostoma","Lunella")

mean<-tapply(alive$final,alive$treatment:alive$Var2, mean )

Alive.Summary$mean<-tapply(alive$final,alive$treatment:alive$Var2, mean )
Alive.Summary$n<-tapply(alive$final,alive$treatment:alive$Var2, length)
Alive.Summary$SD<-tapply(alive$final,alive$treatment:alive$Var2,sd)
Alive.Summary$SE<-Alive.Summary$SD/sqrt(Alive.Summary$n)
Alive.Summary$median<-tapply(alive$final,alive$treatment:alive$Var2,median)

plotA.mean.SE<-
  ggplot(data = Alive.Summary, aes(y = mean,x =treatment, color=treatment))+
  scale_color_manual(guide = 'legend',
                     values=c( "#80cdc1","#dfc27d" ), 
                     name="",
                     breaks=c("Ambient", "Elevated"),
                     labels=c("Current", "Future"))+
  geom_pointrange(ymin=Alive.Summary$mean-Alive.Summary$SE, 
                  ymax=Alive.Summary$mean+Alive.Summary$SE,
                  size= 2.5)+  
  theme_classic(base_size = 20)+
  ylab("Survival (%)")+
  xlab("Treatment")+
  theme(strip.text.x =element_blank(),
        axis.title.x = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank(),
        legend.justification=c(0,0), 
        legend.position=c(0.03,0.06),
        legend.key.size = unit(1.5, "cm"))+
  facet_grid(.~species)+
  ylim(60,100)
  #annotate("text", x = .6, y = 100, label = c("A.",""))

dat_text_A <- data.frame(
  label = c("A.", " "),
  species   = c("Chlorostoma", "Lunella"))

plotA.mean.SE<- plotA.mean.SE + geom_text(
  data    = dat_text_A,
  mapping = aes(x = 0.45, y = 99.4, label = label),
  colour = c( "Black","Black"),
  size=6,
  hjust   = -0.1,
  vjust   = -1)

plotA.mean.SE


Survival<-alive

colnames(Survival)[1] <- 'Tank'
colnames(Survival)[2] <- 'Species'
colnames(Survival)[3] <- 'No.Alive'

write.csv(Survival, "ReadyForPrimer_Survival.csv")



######Feeding Rates#######

#Correcting for algal growth
#creating dataframes for just the control algae and the herbivore exposed algae
snails2<-feeding[which(feeding$species!="Control"),]
snailsControl<-feeding[which(feeding$species=="Control"),]
snails2$species<-as.character(snails2$species)
snails2$species<-as.factor(snails2$species)
str(snails2)
snailsControl$species<-as.character(snailsControl$species)
snailsControl$species<-as.factor(snailsControl$species)
str(snailsControl)


#calculating correction factors based on controls
for(i in 1:length(snailsControl$basket_no)){
  snailsControl$CF[i]<- snailsControl$algae_weight_end[i]/snailsControl$algae_weight_start[i]
}


#using mean correction factors to figure out how much algae was consumed per day
for(i in 1:length(snails2$basket_no)){
  a<-snails2$tank.n[i]
  b<-which(snailsControl$tank.n == a )
  snails2$FeedingRate[i]<- ((snails2$algae_weight_start[i]*(mean(snailsControl$CF[b])))-snails2$algae_weight_end[i])/7
}


#cleaning up dataframe
str(snails2)
Feeding<- snails2
str(Feeding)

#FIGURE Feed w/Mean and SE#####

Feeding$final<-Feeding$FeedingRate/Feeding$wet_weight*100


hist(Feeding$final)

treatment<-c("Ambient","Ambient","Elevated","Elevated")
Feeding.Summary<-as.data.frame(treatment)
Feeding.Summary$species<-c("Chlorostoma","Lunella","Chlorostoma","Lunella")

mean<-tapply(Feeding$final,Feeding$treatment:Feeding$species, mean )

Feeding.Summary$mean<-tapply(Feeding$final,Feeding$treatment:Feeding$species, mean )
Feeding.Summary$n<-tapply(Feeding$final,Feeding$treatment:Feeding$species, length)
Feeding.Summary$SD<-tapply(Feeding$final,Feeding$treatment:Feeding$species,sd)
Feeding.Summary$SE<-Feeding.Summary$SD/sqrt(Feeding.Summary$n)
Feeding.Summary$median<-tapply(Feeding$final,Feeding$treatment:Feeding$species,median)

plotB.mean.SE<-
  ggplot(data = Feeding.Summary, aes(y = mean,x =treatment, color=treatment))+
  scale_color_manual(values=c( "#80cdc1","#dfc27d" ), 
                     name="Treatment",
                     breaks=c("Ambient", "Elevated"),
                     labels=c("Ambient", "Elevated"))+
  geom_pointrange(ymin=Feeding.Summary$mean-Feeding.Summary$SE, 
                  ymax=Feeding.Summary$mean+Feeding.Summary$SE,
                  size=2.5)+  
  theme_classic(base_size = 20)+
  ylab(bquote('Feeding (mg algae ' ~g^-1~day^-1*')'))+
    xlab("Snail species")+
  theme(legend.position="none",
        strip.text.x =element_blank(),
        axis.title.x = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank())+
  facet_grid(.~species)+
  ylim(-.03,.12)#+
  #annotate("text", x = .6, y = 0.12, label = c("B.",""))

dat_text_B <- data.frame(
  label = c("B.", " "),
  species   = c("Chlorostoma", "Lunella"))

plotB.mean.SE<- plotB.mean.SE + geom_text(
  data    = dat_text_B,
  mapping = aes(x = 0.45, y = 0.115, label = label),
  colour = c( "Black","Black"),
  size=6,
  hjust   = -0.1,
  vjust   = -1)

plotB.mean.SE


#######Assembling Compound Figure#######


g1 <- ggplotGrob(plotA.mean.SE)
g2 <- ggplotGrob(plotB.mean.SE)
g3 <- ggplotGrob(plotC.mean.SE)

g <- rbind(g1, g2, g3, size = "first")
g$widths <- unit.pmax(g2$widths, g3$widths)
grid.draw(g)

ggsave("final_graphs_combined.png", g, width = 8.47, height = 20.65)


