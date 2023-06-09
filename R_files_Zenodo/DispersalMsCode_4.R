# Code for "Fungi are more dispersal-limited than bacteria among flowers
# April 29, 2021
# Authors: R. Vannette and M. McMunn

# Setup --read in data from compiled datasheet
#########################
library(ggplot2)
library(forcats)
library(lme4)
library(nlme)
library(visreg)
library(gdata)
library(plyr)
library(reshape2)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dat <- read.csv("FieldData5.19.20.csv")
dim(dat)
dat$Nectary_location<-reorder.factor(dat$Nectary_location, new.order=c("Exposed", "Within_corolla_short","Within_corolla_mid", "Within_corolla_long"))


###########
#Marshall's analyses

########
#####co-occurence of fungi and bacteria
######
# Correlation between bacteria and fungi in flowers
cor.test(dat$BactPres, dat$YeastPres, method = "pearson")
#bacteria and fungi occur in the same flowers more often than expected by chance
#this could Chi-square test, report conditional probabilites

cor.test(dat$LogTotBact, dat$LogTotFungi, method="pearson")

#given bacteria present, are there more or less bacteria when there is fungus?
#more bacteria when fungi present
dat_fungiPres <- dat[dat$LogTotFungi>0,]
dat_bactPres <- dat[dat$LogTotBact>0,]
with(dat_bactPres, tapply(LogTotBact ,YeastPres ,mean, na.rm=TRUE))
anova(with(dat_bactPres, lm(LogTotBact~as.factor(YeastPres))))

#the opposite, given yeast present, are there more or less yeast when there are bacteria?
#more yeast when bacteria present
with(dat_fungiPres, tapply(LogTotFungi ,BactPres ,mean, na.rm=TRUE))
anova(with(dat_fungiPres, lm(LogTotFungi~as.factor(BactPres))))




#calculate probabilities
#clean vectors without NA's in either column
presAbs<-data.frame(BactPres = dat$BactPres, YeastPres = dat$YeastPres)
presAbs<-na.omit(presAbs)

#how often yeast
sum(presAbs$YeastPres) / length(presAbs$YeastPres)
#19.9%

#how often bacteria
sum(presAbs$BactPres) / length(presAbs$BactPres)
#49.3%

#joint probability if independent
(sum(presAbs$BactPres) / length(presAbs$BactPres))*
  (sum(presAbs$YeastPres) / length(presAbs$YeastPres))
#9.8%
Xsq<-chisq.test(presAbs$YeastPres, presAbs$BactPres)
Xsq$statistic
Xsq$p.value

#how often both
sum(presAbs$BactPres&presAbs$YeastPres)/nrow(presAbs)
#14.4%

#how much more frequent is co-occerence than expected?
(sum(presAbs$BactPres&presAbs$YeastPres)/nrow(presAbs)) /

((sum(presAbs$BactPres) / length(presAbs$BactPres))*
  (sum(presAbs$YeastPres) / length(presAbs$YeastPres)))
#46.7% elevated co-occurence

#given yeast, how often bacteria
sum(presAbs[presAbs$YeastPres,"BactPres" ]) / length(dat[dat$YeastPres ,"BactPres" ])
#69.0%


#given bacteria, how often yeast
sum(presAbs[presAbs$BactPres,"YeastPres" ]) / length(dat[dat$BactPres ,"YeastPres" ])
#28.7%

#given one or other, how often both?
sum(presAbs$YeastPres & presAbs$BactPres) /
  nrow(presAbs[presAbs$YeastPres | presAbs$BactPres, ])
#26.2%


#how often one or other
nrow(presAbs[presAbs$YeastPres | presAbs$BactPres, ]) / nrow(presAbs)
#54.8%

#how often just bacteria, but not yeast
sum(presAbs$BactPres&!presAbs$YeastPres) / nrow(presAbs)
#34.8%

#how often just yeast, but not bacteria
sum(presAbs$YeastPres&!presAbs$BactPres) / nrow(presAbs)
#5.5%

#how often neither
sum(!presAbs$YeastPres&!presAbs$BactPres) / nrow(presAbs)
#45.2%




####################
# Figure 3a: plot of microbial incidence among species
####################

fungip <- ddply(dat, .(Flower.Abbreviation),summarize,
                meanf = mean(YeastPres, na.rm=TRUE),
                meanb = mean(BactPres, na.rm=TRUE))
ggplot(dat, aes(x=fct_reorder(Flower.Abbreviation, BactPres, .fun='mean', na.rm=TRUE) , y=BactPres))+
  #geom_point()+
  #stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=0.15, color="black")+
  stat_summary(fun.data = mean_cl_normal, geom = "point", width=0.15, color="blue")+
  geom_point(data=fungip, aes(x=fct_reorder(Flower.Abbreviation, meanb, .fun='mean', na.rm=TRUE) , y=meanf), color="orange")+
  # stat_summary(data=comb2,aes(x=fct_reorder(Flower.Abbreviation, BactPres, .fun='mean', na.rm=TRUE), y=YeastPres), fun.data = mean_cl_normal, geom = "point", width=0.15, color="orange")+
  # facet_wrap(~Site)+
  ylim(0,1)+
  theme_bw()+  
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Prop. nectar w. microbes")+
  xlab("Plant species")+
  ggsave("Plantspecies w proportion microbes.pdf", width=6, height=4)

myord <- order(meanb$Flower.Abbreviation, .fun='mean')

dat_both_melt <- melt(dat[,c("Flower.Abbreviation","Site", "Year", "BactPres", "YeastPres")], id.vars=c("Flower.Abbreviation","Site", "Year"))

out <- glm(value ~ variable *(Site+ Flower.Abbreviation)+Year, data=dat_both_melt, family="binomial")
drop1(out, .~., test="Chi")
summary(out)

# examine just bacteria
summary(with(dat , glm(BactPres~ Site+Year+ Flower.Abbreviation,family = "binomial" )))
# bacteria differ between sites (less abundant at STebbins) and between years but interaction is not significant. 
#Significant variation among plant species

# examine just fungi
summary(with(dat , glm(YeastPres~ Site+Year+ Flower.Abbreviation,family = "binomial" )))
# fungi do not differ in incidence between sites or years but do vary among plant species



####################
# Figure 3b: plot of microbial abundance among species
####################

dat_noze <- dat[dat$LogTotBact>0|dat$LogTotFungi>0,]
microbeab <- ddply(dat_noze, .(Flower.Abbreviation),summarize,
                   meanb = mean(LogTotBact, na.rm=TRUE),
                   meanf = mean(LogTotFungi, na.rm=TRUE),
                   seb = sd(LogTotBact)/length(LogTotBact), 
                   sef = sd(LogTotFungi)/length(LogTotFungi)
)
microbeab<-microbeab[-43,]
microbe_all <-merge(fungip, microbeab, by="Flower.Abbreviation")

ggplot(microbe_all, aes(x=fct_reorder(Flower.Abbreviation, meanb.x) , y=meanb.y))+
  geom_point( color="blue")+
  geom_point(aes(y=meanf.y), color="orange")+
  geom_errorbar(aes(ymin=meanb.y-seb, ymax=meanb.y+seb), color="blue")+
  geom_errorbar(aes(ymin=meanf.y-sef, ymax=meanf.y+sef), color="orange")+
  theme_bw()+  
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Average Log10(CFU+1) in flowers w. microbes")+
  xlab("Plant species")+
  ggsave("Plantspecies w abundance microbes.pdf", width=6, height=4)



dat_both_melt <- melt(dat_noze[,c("Flower.Abbreviation","Site", "Year", "LogTotBact", "LogTotFungi")], id.vars=c("Flower.Abbreviation","Site", "Year"))

out <- lm(value ~ variable *(Flower.Abbreviation)+Site+Year, data=dat_both_melt)
drop1(out, .~., test="F")
summary(out)


#### Temporal models (supplementary materials) ########

# presence
drop1(with(dat , glm(BactPres~ Site+Julian *Year,family = "binomial" )), .~., test="Chi")
summary(with(dat , glm(BactPres~ Site+Julian *Year,family = "binomial" )))

drop1(with(dat , glm(YeastPres~ Site*Julian +Year,family = "binomial" )), .~., test="Chi")
summary(with(dat , glm(YeastPres~ Site*Julian+Year+Julian,family = "binomial" )))

summary(with(dat, lm(YeastPres~Site)))       
with(dat, tapply(BactPres ,Site ,mean, na.rm=TRUE))
with(dat, tapply(YeastPres ,Site ,mean, na.rm=TRUE))

#### use AIC for model selection for bacteria presence ####
AllInt <- glm(BactPres~ Julian*Site+Site*Year+Julian *Year,data=dat, family = "binomial" )
NoInt <- glm(BactPres~ Julian+Site+Year,data=dat, family = "binomial" )
RedInt <- glm(BactPres~ Julian+Site+Site*Year+Julian *Year,data=dat, family = "binomial" )

AIC(AllInt, NoInt, RedInt)
drop1(RedInt, .~., test="Chi")
summary(RedInt)

#best-fit model for bacterial presence 

#### AIC model selection for fungi presence ####
AllInt <- glm(YeastPres~ Julian*Site+Site*Year+Julian *Year,data=dat, family = "binomial" )
NoInt <- glm(YeastPres~ Julian+Site+Year,data=dat, family = "binomial" )
RedInt <- glm(YeastPres~ Julian*Site+Site+Year*Julian *Year,data=dat, family = "binomial" )

AIC(AllInt, NoInt, RedInt)
drop1(RedInt, .~., test="Chi")
summary(RedInt)

#### use AIC for model selection for bacteria abundance ####
AllInt <- lm(LogTotBact~ Julian*Site+Site*Year+Julian *Year,data=dat )
NoInt <- lm(LogTotBact~ Julian+Site+Year,data=dat )
RedInt <- lm(LogTotBact~ Julian+Site+Site*Year+Julian *Year,data=dat )

AIC(AllInt, NoInt, RedInt)
drop1(RedInt, .~., test="F")
summary(RedInt)

#### use AIC for model selection for fungi abundance ####
AllInt <- lm(LogTotFungi~ Julian*Site+Site*Year+Julian *Year,data=dat )
NoInt <- lm(LogTotFungi~ Julian+Site+Year,data=dat )
RedInt <- lm(LogTotFungi~ Julian*Site+Site+Year*Julian *Year,data=dat )

AIC(AllInt, NoInt, RedInt)
drop1(RedInt, .~., test="F")
summary(RedInt)


# add plant species; does this remove effect of site or date? 
drop1(with(dat , glm(BactPres~ Flower.Abbreviation+Site+Julian+Julian *Year,family = "binomial" )), .~., test="Chi")
summary(with(dat , glm(BactPres~ Flower.Abbreviation+Site+Julian +Year,family = "binomial" )))

fit <- glm(BactPres~ Flower.Abbreviation+Site*Year+Julian *Year,family = "binomial", data=dat)
drop1(fit, .~., test="Chi")
vif(fit)
visreg(fit, "Julian", by="Year")
visreg(fit, "Site", by="Year")

drop1(with(dat , glm(YeastPres~ Lognect+Flower.Abbreviation+Site+Julian +Year,family = "binomial" )), .~., test="Chi")
summary(with(dat , glm(YeastPres~ Flower.Abbreviation+Site+Julian +Year,family = "binomial" )))

fit <- glm(YeastPres~ Site+Julian*Year,family = "binomial", data=dat)
drop1(fit, .~., test="Chi")

vif(fit)
visreg(fit)
visreg(fit, "Julian", by="Year")
visreg(fit, "Julian", by="Site")

summary(with(dat , glm(BactPres~ Site+Year+Julian,family = "binomial" )))
summary(with(dat , glm(YeastPres~ Site+Year+Julian,family = "binomial" )))
# significant positive effect of Julian date on the abundance of both yeasts and bacteria

# abundance
drop1(with(dat , lm(LogTotBact~ Site*Year+Julian+Julian *Year )), .~., test="F")
summary(with(dat , lm(LogTotBact~ Site*Year+Julian+Julian *Year )))

drop1(with(dat , lm(LogTotFungi~ Site+Year+Julian *Year )), .~., test="F")
summary(with(dat , lm(LogTotFungi~ Site+Year+Julian *Year)))

# add plant species
drop1(with(dat , lm(LogTotBact~ Site*Year+Julian+Julian *Year+Flower.Abbreviation )), .~., test="F")
summary(with(dat , lm(LogTotBact~ Site*Year+Julian+Julian *Year+Flower.Abbreviation )))


drop1(with(dat , lm(LogTotFungi~ Site+Year+Julian *Year +Flower.Abbreviation)), .~., test="F")
summary(with(dat , lm(LogTotFungi~ Site+Year+Julian *Year +Flower.Abbreviation)))



 # estimate if there is a temporal trend
m1 <- lm(LogTotBact~ Julian, dat)
par(mfrow=c(2,2))
plot(m1)
# examine residuals
par(mfrow=c(1,1))
plot(residuals(m1))
acf(residuals(m1))

#look for the same in fungi
m1 <- lm(LogTotFungi~ Julian, dat)
par(mfrow=c(2,2))
plot(m1)
# examine residuals
par(mfrow=c(1,1))
plot(residuals(m1))
acf(residuals(m1))

library(mgcv)
mdl.ac <- bam(LogTotBact ~Julian+Site+Year+Flower.Abbreviation, data=dat, 
              correlation = corAR1(form=~Julian),
              na.action=na.omit)
summary(mdl.ac)

mdl.ac <- bam(LogTotFungi ~Julian+Site+Year+Flower.Abbreviation, data=dat, 
              correlation = corAR1(form=~Julian),
              na.action=na.omit)
summary(mdl.ac)

#  accounting for temporal autocorrelation (Julian date) does not significantly affect model results 



#### Figure 4a: plot of microbial incidence by flower traits ####

dat_both_melt <- melt(dat[,c("Flower.Abbreviation","Site", "Year", "Nectary_location", "BactPres", "YeastPres")], 
                      id.vars=c("Flower.Abbreviation", "Site", "Year", "Nectary_location"))


out <- glm(value ~ variable *(Nectary_location+Site)+Year, data=dat_both_melt, family="binomial")
drop1(out, .~., test="Chi")

summary(out)

out <- glm(value ~ variable *(Site)+Year, data=dat_both_melt, family="binomial")
drop1(out, .~., test="Chi")

summary(out)

# redo this analysis with just bacteria or just fungi
out <- glm(BactPres ~  Nectary_location, data=dat, family="binomial")
drop1(out, .~., test="Chi")

summary(out)
# bacteria found more commonly in exposed or long corollas
out <- glm(YeastPres ~  Nectary_location, data=dat, family="binomial")
drop1(out, .~., test="Chi")

summary(out)

out <- glm(LogTotBact ~  Nectary_location, data=dat)
drop1(out, .~., test="F")
summary(out)

out <- glm(LogTotFungi ~  Nectary_location, data=dat)
drop1(out, .~., test="F")
summary(out)


# species-level analyses: 
# redo this analysis with just bacteria or just fungi
out <- glm(BactPres ~  Flower.Abbreviation, data=dat, family="binomial")
drop1(out, .~., test="Chi")

summary(out)
# bacteria found more commonly in exposed or long corollas
out <- glm(YeastPres ~  Flower.Abbreviation, data=dat, family="binomial")
drop1(out, .~., test="Chi")

summary(out)

out <- glm(LogTotBact ~  Flower.Abbreviation, data=dat)
drop1(out, .~., test="F")
summary(out)

out <- glm(LogTotFungi ~  Flower.Abbreviation, data=dat)
drop1(out, .~., test="F")
summary(out)


# Figure


ggplot(dat_both_melt, aes(x=Nectary_location, y=value, color=variable))+
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=0.15)+
  stat_summary(fun.data = mean_cl_normal, geom = "point", width=0.15)+
  theme_bw()+
  scale_color_manual(values=c("blue", "orange"))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab("P(microbial incidence)")+
  ylim(0,1)+
  facet_wrap(~variable, scales="free")+
  ggsave("MicrobeIncidence_location.pdf", width=5.5, height=3.5)


####################
# Figure 2a: plot of microbial abundance by flower traits
####################
dat_noze <- dat[dat$LogTotBact>0|dat$LogTotFungi>0,]

dat_both_melt <- melt(dat_noze[,c("Flower.Abbreviation","Site", "Nectary_location","Year", "Color", "LogTotBact", "LogTotFungi")], 
                      id.vars=c("Flower.Abbreviation", "Site", "Year", "Nectary_location", "Color"))

dat_both_melt <-dat_both_melt[is.na(dat_both_melt$Color)==F,]

ggplot(dat_both_melt, aes(x=Color, y=value, color=variable))+
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=0.15)+
  stat_summary(fun.data = mean_cl_normal, geom = "point", width=0.15)+
  theme_bw()+
  scale_color_manual(values=c("blue", "orange"))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab("Log10(CFU+1)")+
  ylim(0,4)+
  ggsave("MicrobeAbundance_color.pdf", width=5.5, height=3.5)

out <- lm(value ~ variable *(Color+Year+ Site), data=dat_both_melt)
drop1(out, .~., test="F")



ggplot(dat_both_melt, aes(x=Nectary_location, y=value, color=variable))+
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=0.15)+
  stat_summary(fun.data = mean_cl_normal, geom = "point", width=0.15)+
  theme_bw()+
  scale_color_manual(values=c("blue", "orange"))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab("Log10(CFU+1)")+
  #ylim(0,3.5)+
  facet_wrap(~variable, scales="free")+
  ggsave("MicrobeAbundance_location.pdf", width=5.5, height=3.5)

out <- lm(value ~ variable *(Nectary_location)+Year+Site, data=dat_both_melt)
drop1(out, .~., test="F")




#####################
# Supplemental figures and tables
#####################

# Supplementary Figure 1a Seasonal patterns in microbial incidence
# NOTE: only one species was sampled at late dates, so seasonal trends in yeast/bacterial abundance at the end of the season 
# are skewed by this sampling artifact. 

dat_both_melt2 <- melt(dat[,c("Flower.Abbreviation","Site","Julian", "Year", "BactPres", "YeastPres")], id.vars=c("Flower.Abbreviation", "Site", "Julian", "Year"))

p1a <- ggplot(data = dat_both_melt2, aes(x=Julian, y=value, color=variable))+
  # geom_line()+
  theme_bw()+
  #geom_point(data=comb2, aes(x=Julian, y=value, color=variable), alpha=0.2)+
  geom_point(alpha=0.2)+
  geom_smooth(method = "loess")+
  facet_grid(Year~Site)+
  ylab("P(microbes)")+
  scale_color_manual(values=c("blue", "orange"))+
  ggsave("Prop(CFUs) by site and year loess.pdf", height=3, width=5)

p1a

# Supplementary Figure 1b Seasonal patterns in microbial abundance
# NOTE: only one species was sampled at late dates, so seasonal trends in yeast/bacterial abundance at the end of the season 
# are skewed by this sampling artifact. 
dat_both_melt <- melt(dat[,c("Flower.Abbreviation","Site","Julian", "Year", "LogTotBact", "LogTotFungi")], id.vars=c("Flower.Abbreviation", "Site", "Julian", "Year"))
#dat_both_melt_nolate <- dat_both_melt[dat_both_melt$Julian<200,]
dat_both_melt_noz <- dat_both_melt[dat_both_melt$value>0,]
dat_both_melt_noz <- dat_both_melt_noz[is.na(dat_both_melt_noz$value)==FALSE,]

p2 <- ggplot(data = dat_both_melt_noz, aes(x=Julian, y=value, color=variable))+
  # geom_line()+
  theme_bw()+
  #geom_point(data=comb2, aes(x=Julian, y=value, color=variable), alpha=0.2)+
  geom_point(alpha=0.2)+
  geom_smooth(method = "loess")+
  facet_grid(Year~Site, drop=TRUE)+
  ylab("Log10(CFU+1)")+
  scale_color_manual(values=c("blue", "orange"))+
  ggsave("CFUs by site and year loess.pdf", height=3, width=5)

p2

# Pearson correlations for each species between bacterial and fungal CFUS in all flowers
corrResSp = list()
counter = 0 

datCor<-dat
datCor$Flower.Abbreviation <-drop.levels(dat$Flower.Abbreviation)

for (i in unique(datCor$Flower.Abbreviation)){
  counter = counter + 1
  x = as.numeric(datCor[datCor$Flower.Abbreviation == i,]$LogTotFungi)
  y = as.numeric(datCor[datCor$Flower.Abbreviation == i,]$LogTotBact)
  corrResSp[[counter]] <-cor.test(x,y,method="pearson", na.action=na.omit)
}

corResDF <- data.frame(NULL)
for (i in 1:length(corrResSp)){
  corResDF[i,1] <- unique(datCor$Flower.Abbreviation)[i]
  corResDF[i,2] <-corrResSp[[i]]$p.value
  corResDF[i,3] <- corrResSp[[i]]$estimate
}

corResDF$adj <- p.adjust(corResDF$V2, method='fdr')
corResDF$plotSig<-ifelse(corResDF$adj<0.05 , "black", "grey70")

pdf("Abundance correlation Coef within spp.pdf")
plot(corResDF$V3, col = corResDF$plotSig, pch = 16, xlab = "index", ylab = "R", 
    xlim=c(-5,40))
abline(h=0)
text(x=seq_along(corResDF$V3), y=corResDF$V3, labels=corResDF$V1, 
     cex=0.8, font=2, col=corResDF$plotSig, pos=2)
dev.off()
# Pearson correlations for each species between bacterial and fungal presence in all flowers
corrResSp = list()
counter = 0 

for (i in unique(dat$Flower.Abbreviation)){
  counter = counter + 1
  x = as.numeric(datCor[datCor$Flower.Abbreviation == i,]$YeastPres)
  y = as.numeric(datCor[datCor$Flower.Abbreviation == i,]$BactPres)
  corrResSp[[counter]] <-cor.test(x,y,method="pearson", na.action=na.omit)
}

corResDFpresence <- data.frame(NULL)
for (i in 1:length(corrResSp)){
  corResDFpresence[i,1] <- unique(datCor$Flower.Abbreviation)[i]
  corResDFpresence[i,2] <-corrResSp[[i]]$p.value
  corResDFpresence[i,3] <- corrResSp[[i]]$estimate
}

corResDFpresence$adj <- p.adjust(corResDFpresence$V2, method='fdr')
corResDFpresence$plotSig<-ifelse(corResDFpresence$adj<0.05 , "black", "grey70")

pdf("PresAbs Coef within spp - PA.pdf")
plot(corResDFpresence$V3, col = corResDFpresence$plotSig, pch = 16, 
     xlab = "index", ylab = "R", xlim=c(-5,40))
abline(h=0)
text(x=seq_along(corResDFpresence$V3), y=corResDFpresence$V3, labels=corResDFpresence$V1, 
     cex=0.8, font=2, col=corResDFpresence$plotSig, pos=2)
dev.off()





