###################################################################################################################
###### Title: Tropical-temperate comparisons in insect seed predation vary between study levels and years #########
###### Author: Wenlan Wu et al.
###### Updated Date: 14 August 2022
###################################################################################################################

rm(list = ls()) # remove all the data in the Global Environment

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Rcode/Wenlan_granivory_Dryad") # CHANGE to your working directory


##########################################  1. Community-wide levels  analyses ####################################

Community_wide <- read.csv('Seed_per_batch_2019_2020.csv', header = T, row.names = 1) 

head(Community_wide)
# Seedmass_insect = The mass of seeds detected with insect predators, which related to Incidence_of_seed_predator
# Seedmass_attack = The mass of seeds showing signs of insect attack, which related to Seed_predation_rate


str(Community_wide)

Community_wide$Site <- as.factor(Community_wide$Site)

##### 1.1 calculate Annual and Community-wide measurements (part of Table 1) #####

##### 1.1.1 Incidence of seed predator in 2019  ####

# CBS

Community_wide2019 <- subset(Community_wide, Year=="2019")

Community_wide2019_CBS <- subset(Community_wide2019, Site=="CBS")


head(Community_wide2019_CBS)


sum(Community_wide2019_CBS$Seednumber.total) # 15614

sum(Community_wide2019_CBS$Seedmass.total) # 1093.5g


sum(Community_wide2019_CBS$Seedmass_insect)/sum(Community_wide2019_CBS$Seedmass.total) # 0.1792985


mean(Community_wide2019_CBS$Incidence_of_seed_predator) # 0.1755801

SD <- sd(Community_wide2019_CBS$Incidence_of_seed_predator) # 0.02841021

SE <- SD/sqrt(7)  # 0.01073805



# XSBN

Community_wide2019_XSBN <- subset(Community_wide2019, Site=="XSBN")


sum(Community_wide2019_XSBN$Seednumber.total) # 15128

sum(Community_wide2019_XSBN$Seedmass.total) # 3958.592g


sum(Community_wide2019_XSBN$Seedmass_insect)/sum(Community_wide2019_XSBN$Seedmass.total) # 0.08401151


mean(Community_wide2019_XSBN$Incidence_of_seed_predator) # 0.09729547

SD <- sd(Community_wide2019_XSBN$Incidence_of_seed_predator) # 0.08767949

SE <- SD/sqrt(26)  # 0.01719536




##### 1.1.2 Incidence of seed predator in 2020  ####

# CBS

Community_wide2020 <- subset(Community_wide, Year=="2020")

Community_wide2020_CBS <- subset(Community_wide2020, Site=="CBS")


sum(Community_wide2020_CBS$Seednumber.total) # 148586

sum(Community_wide2020_CBS$Seedmass.total) # 3262.18g


sum(Community_wide2020_CBS$Seedmass_insect)/sum(Community_wide2020_CBS$Seedmass.total) # 0.02001744


mean(Community_wide2020_CBS$Incidence_of_seed_predator) # 0.01843841

SD <- sd(Community_wide2020_CBS$Incidence_of_seed_predator) # 0.01533273

SE <- SD/sqrt(13)  # 0.004252534


# XSBN

Community_wide2020_XSBN <- subset(Community_wide2020, Site=="XSBN")


sum(Community_wide2020_XSBN$Seednumber.total) # 23589

sum(Community_wide2020_XSBN$Seedmass.total) # 4346.749g


sum(Community_wide2020_XSBN$Seedmass_insect)/sum(Community_wide2020_XSBN$Seedmass.total) # 0.06717901


mean(Community_wide2020_XSBN$Incidence_of_seed_predator) # 0.06784181

SD <- sd(Community_wide2020_XSBN$Incidence_of_seed_predator) # 0.04494625

SE <- SD/sqrt(27)  # 0.008649911



##### 1.1.3 Seed predation rates in 2020 ####

# CBS

sum(Community_wide2020_CBS$Seedmass_attack)/sum(Community_wide2020_CBS$Seedmass.total) # 0.05231729


mean(Community_wide2020_CBS$Seed_predation_rate) # 0.05727119

SD <- sd(Community_wide2020_CBS$Seed_predation_rate) # 0.03673546

SE <- SD/sqrt(13)  # 0.01018858


# XSBN

sum(Community_wide2020_XSBN$Seedmass_attack)/sum(Community_wide2020_XSBN$Seedmass.total) # 0.247312


mean(Community_wide2020_XSBN$Seed_predation_rate) # 0.2465907

SD <- sd(Community_wide2020_XSBN$Seed_predation_rate) # 0.08440964

SE <- SD/sqrt(27)  # 0.01624464


####  *Simulated data of the Korean Pine Pinus koraiensis ####

# change the seed predation rates of late September (2020/9/28) and early October (2020/10/14) batches in CBS in 2020

Community_wide2020.pine <- Community_wide2020

Community_wide2020.pine[Community_wide2020.pine$Site == 'CBS' & Community_wide2020.pine$Date == '2020/9/28', ]$Seed_predation_rate <- 0.197

Community_wide2020.pine[Community_wide2020.pine$Site == 'CBS' & Community_wide2020.pine$Date == '2020/10/14', ]$Seed_predation_rate <- 0.197


Community_wide2020_CBS_pine <- subset(Community_wide2020.pine, Site=="CBS")

mean(Community_wide2020_CBS_pine$Seed_predation_rate) # 0.08104514

SD <- sd(Community_wide2020_CBS_pine$Seed_predation_rate) # 0.06253433

SE <- SD/sqrt(13)  # 0.0173439



#####  1.2 Generalized Least Squares (GLS) analyse  #####

#### 1.2.1 Incidence of seed predator in 2019 ####

# Visualization using the violin Plot (Figure 2a)

library(vioplot)

# par(mfrow=c(1,2), mar=c(2,5,2,0), oma=c(1,0,1,2), bty="l") # mar and oma settings correspond to lower, left, upper and right

vioplot(Community_wide2019$Incidence_of_seed_predator ~ Community_wide2019$Site,
        
        col=c("#00BFC4", "#F8766D"), ylim=c(0,0.35),
        
        xlab="", ylab="Incidence of seed predator") # XSBN < CBS


head(Community_wide2019)

library(lubridate)

# Converse the Date to the number, i.e. the days since 1 January of the year.

x <- as.Date(Community_wide2019$Date)

Community_wide2019$Date <- yday(x)

str(Community_wide2019)


# run gls

library(nlme)  

Community.mod2019 <- gls(log(sqrt(Incidence_of_seed_predator)+1) ~ Site + Date, weights=varIdent(form= ~1|Site),
                         
                         correlation=corAR1(form = ~ Date | Site), data=Community_wide2019)

# model diagnostic plot

plot(Community.mod2019) # ok

qqnorm(Community.mod2019$residuals) # basically ok

qqline(Community.mod2019$residuals)

# model results

summary(Community.mod2019)

anova(Community.mod2019)  # Site <.001 (part of Table 2)



#### 1.2.2 Incidence of seed predator in 2020 ####

# Visualization using the violin Plot (Figure 2b)

vioplot(Community_wide2020$Incidence_of_seed_predator ~ Community_wide2020$Site,
        
        col=c("#00BFC4", "#F8766D"), ylim=c(0,0.35), yaxt="n",
        
        xlab = "", ylab = "")

axis(side=1, labels=c('CBS', 'XSBN'), at=c(1, 2)) # XSBN > CBS



# Converse the Date to the number, i.e. the days since 1 January of the year.

x <- as.Date(Community_wide2020$Date)

Community_wide2020$Date <- yday(x)

# run gls 

Community.mod2020 <- gls(log(sqrt(Incidence_of_seed_predator)+1) ~ Site + Date, 
                         
                         correlation=corAR1(form = ~ Date | Site), data=Community_wide2020)

# model diagnostic plot

plot(Community.mod2020) # ok 

qqnorm(Community.mod2020$residuals) # ok

qqline(Community.mod2020$residuals)

# model results

summary(Community.mod2020)

anova(Community.mod2020)       # Site= 0.0001 (part of Table 2)



##### 1.2.3 Seed predation rates in 2020 ####

# Visualization using the violin Plot (Figure S1a)

# par(mfrow=c(1,2),mar=c(2,5,2,0),oma=c(1,0,1,2),bty="l")

vioplot(Community_wide2020$Seed_predation_rate ~ Community_wide2020$Site,
        
        col=c("#00BFC4", "#F8766D"), ylim=c(0,0.5),
        
        xlab="", ylab="Seed predation rate")  # XSBN > CBS 



# run gls

Community.mod2020.Spr <- gls(log(sqrt(Seed_predation_rate)+1) ~ Site + Date, weights=varIdent(form= ~1|Site),
                             
                             correlation = corAR1(form = ~ Date | Site), data = Community_wide2020)


# model diagnostic plot

plot(Community.mod2020.Spr) # ok 

qqnorm(Community.mod2020.Spr$residuals) # basically ok

qqline(Community.mod2020.Spr$residuals)

# model results

summary(Community.mod2020.Spr)

anova(Community.mod2020.Spr)       # Site= <.0001 (part of Table 2)



#### *Simulated data of the Korean Pine Pinus koraiensis ####

# Visualization using the violin Plot (Figure S1b)

vioplot(Community_wide2020.pine$Seed_predation_rate ~ Community_wide2020.pine$Site,
        
        col=c("#00BFC4", "#F8766D"), ylim=c(0,0.5), 
        
        yaxt = "n", xlab="", ylab="")

axis(side = 1, labels = c('CBS', 'XSBN'), at = c(1, 2)) # XSBN > CBS 



# Converse the Date to the number, i.e. the days since 1 January of the year.

x <- as.Date(Community_wide2020.pine$Date)

Community_wide2020.pine$Date <- yday(x)

str(Community_wide2020.pine)


# run gls

Community.mod2020.pine <- gls(log(sqrt(Seed_predation_rate)+1) ~ Site + Date,
                              
                              correlation = corAR1(form = ~ Date | Site), data=Community_wide2020.pine)

# model diagnostic plot

plot(Community.mod2020.pine$residuals) # basically ok

qqnorm(Community.mod2020.pine$residuals)

qqline(Community.mod2020.pine$residuals) # basically ok

# model results

summary(Community.mod2020.pine)

anova(Community.mod2020.pine)       # Site= <.0001




#########################################  2. Cross-species levels analyses  ######################################

Cross_species <- read.csv('Seed_per_species_2019_2020.csv', header = T, row.names = 1)

head(Cross_species)
# Seedmass_insect = The mass of seeds detected with insect predators, which related to Incidence_of_seed_predator
# Seedmass_attack = The mass of seeds showing signs of insect attack, which related to Seed_predation_rate
# Seednumber_insect = The number of seeds detected with insect predators, which related to Incidence_of_seed_predator.unweighted
# Seednumber_attack = The number of seeds showing signs of insect attack, which related to Seed_predation_rate.unweighted


str(Cross_species)

Cross_species$Species <- as.factor(Cross_species$Species)

Cross_species$Site <- as.factor(Cross_species$Site)

Cross_species$Fruit_type <- as.factor(Cross_species$Fruit_type)

Cross_species$Year <- as.factor(Cross_species$Year)

##### 2.1 calculate Cross-species measurements (part of Tables 1 and S1)  #####

##### 2.1.1 Incidence of seed predator in 2019  ####

# Screen for species with >= 50 seeds in 2019

Cross_species2019 <- subset(Cross_species, Year=="2019")

library(dplyr)

Cross_species2019_50 <- filter(.data = Cross_species2019, Seednumber >= 50) # 54 species, CBS=7  XSBN=47

summary(Cross_species2019_50)


# XSBN

Cross_species2019_50_XSBN <- subset(Cross_species2019_50, Site=="XSBN")

# weighted by seed mass

mean(Cross_species2019_50_XSBN$Incidence_of_seed_predator) # 0.07641621

SD <- sd(Cross_species2019_50_XSBN$Incidence_of_seed_predator) # sd = 0.1000772

SE <- SD/sqrt(47)  # se = 0.01459775

# unweighted by seed mass

mean(Cross_species2019_50_XSBN$Incidence_of_seed_predator.unweighted) # 0.06397508

SD <- sd(Cross_species2019_50_XSBN$Incidence_of_seed_predator.unweighted) # sd = 0.08347236

SE <- SD/sqrt(47)  # se = 0.0121757



# CBS

Cross_species2019_50_CBS <- subset(Cross_species2019_50, Site=="CBS")

# weighted by seed mass

mean(Cross_species2019_50_CBS$Incidence_of_seed_predator) # 0.09417326

SD <- sd(Cross_species2019_50_CBS$Incidence_of_seed_predator) # sd = 0.08006068

SE <- SD/sqrt(7)  # se = 0.03026009

# unweighted by seed mass

mean(Cross_species2019_50_CBS$Incidence_of_seed_predator.unweighted) # 0.08679483

SD <- sd(Cross_species2019_50_CBS$Incidence_of_seed_predator.unweighted) # sd = 0.0771702

SE <- SD/sqrt(7)  # se = 0.0291676



##### 2.1.2 Incidence of seed predator in 2020  ####

# Screen for species with >= 50 seeds in 2020

Cross_species2020 <- subset(Cross_species, Year=="2020")

Cross_species2020_50 <- filter(.data = Cross_species2020, Seednumber >= 50) # 65 species, CBS=9  XSBN=56

summary(Cross_species2020_50)


# XSBN

Cross_species2020_50_XSBN <- subset(Cross_species2020_50, Site=="XSBN")

# weighted by seed mass

mean(Cross_species2020_50_XSBN$Incidence_of_seed_predator) # 0.05105389

SD <- sd(Cross_species2020_50_XSBN$Incidence_of_seed_predator) # sd =  0.05929449

SE <- SD/sqrt(56)  # se = 0.007923559

# unweighted by seed mass

mean(Cross_species2020_50_XSBN$Incidence_of_seed_predator.unweighted) # 0.04263666

SD <- sd(Cross_species2020_50_XSBN$Incidence_of_seed_predator.unweighted) # sd = 0.0570705

SE <- SD/sqrt(56)  # se = 0.007626367



# CBS

Cross_species2020_50_CBS <- subset(Cross_species2020_50, Site=="CBS")

# weighted by seed mass

mean(Cross_species2020_50_CBS$Incidence_of_seed_predator) # 0.02417197

SD <- sd(Cross_species2020_50_CBS$Incidence_of_seed_predator) # sd = 0.02358349

SE <- SD/sqrt(9)  # se = 0.007861162

# unweighted by seed mass

mean(Cross_species2020_50_CBS$Incidence_of_seed_predator.unweighted) # 0.02426702

SD <- sd(Cross_species2020_50_CBS$Incidence_of_seed_predator.unweighted) # sd = 0.02306726

SE <- SD/sqrt(9)  # se = 0.007689087



##### 2.1.3 Seed predation rate in 2020  ####

# XSBN

# weighted by seed mass

mean(Cross_species2020_50_XSBN$Seed_predation_rate) # 0.1937077

SD <- sd(Cross_species2020_50_XSBN$Seed_predation_rate) # sd = 0.1667389

SE <- SD/sqrt(56)  # se = 0.02228142


# unweighted by seed mass

mean(Cross_species2020_50_XSBN$Seed_predation_rate.unweighted) # 0.1522043

SD <- sd(Cross_species2020_50_XSBN$Seed_predation_rate.unweighted) # sd = 0.1465309

SE <- SD/sqrt(56)  # se = 0.01958102



# CBS

# weighted by seed mass

mean(Cross_species2020_50_CBS$Seed_predation_rate) # 0.07011097

SD <- sd(Cross_species2020_50_CBS$Seed_predation_rate) # sd = 0.05531559

SE <- SD/sqrt(9)  # se = 0.01843853


# unweighted by seed mass

mean(Cross_species2020_50_CBS$Seed_predation_rate.unweighted) # 0.07833145

SD <- sd(Cross_species2020_50_CBS$Seed_predation_rate.unweighted) # sd = 0.05206354

SE <- SD/sqrt(9)  # se = 0.01735451



####  *Simulated data of the Korean Pine Pinus koraiensis ####

# Add simulated data of the Korean Pine Pinus koraiensis in 2020
# seed predation rate = 0.197, mean seedmass = 0.49 g, fruit type = Dry_fruit

Cross_species2020_50.pine <- rbind(Cross_species2020_50, 
                                   data.frame("Species" = "Pinus_koraiensis", "Site" = "CBS", "Fruit_type"="Dry_fruit",
                                              "Year"="2020", "Seednumber_insect"=NA, "Seednumber_attack"=NA, "Seedmass.total"=NA,
                                              "Seednumber"=NA, "Seedmass_insect"=NA, "Seedmass_attack"=NA, "Seedmass.mean"=0.49,
                                              "Incidence_of_seed_predator"=NA, "Seed_predation_rate"=0.197, 
                                              "Incidence_of_seed_predator.unweighted"=NA, "Seed_predation_rate.unweighted"=0.197)
                                   )

str(Cross_species2020_50.pine)

Cross_species2020_50_CBS.pine<- subset(Cross_species2020_50.pine, Site == "CBS")


# weighted by seed mass

mean(Cross_species2020_50_CBS.pine$Seed_predation_rate) # 0.08279988

SD <- sd(Cross_species2020_50_CBS.pine$Seed_predation_rate) # sd = 0.06580211

SE <- SD/sqrt(10)  # se = 0.02080845

# unweighted by seed mass

mean(Cross_species2020_50_CBS.pine$Seed_predation_rate.unweighted) # 0.09019831

SD <- sd(Cross_species2020_50_CBS.pine$Seed_predation_rate.unweighted) # sd = 0.06178718

SE <- SD/sqrt(10)  # se = 0.01953882



#####  2.2 Phylogenetic Generalized Least Squares (PGLS) analyses  #####

##### 2.2.1 Incidence of seed predator in 2019 #####

#### Visualization using the violin Plot


# weighted by seed mass (Figure 3a)

# par(mfrow=c(1,2), mar=c(2,5,2,0), oma=c(1,0,1,2), bty="l")

vioplot(Cross_species2019_50$Incidence_of_seed_predator ~ Cross_species2019_50$Site,
        
        col=c("#00BFC4", "#F8766D"), ylim=c(0,0.5),
        
        xlab="", ylab="Incidence of seed predator")


# unweighted by seed mass (Figure S3a)

# par(mfrow=c(1,2), mar=c(2,5,2,0), oma=c(1,0,1,2), bty="l")

vioplot(Cross_species2019_50$Incidence_of_seed_predator.unweighted ~ Cross_species2019_50$Site,
        
        col=c("#00BFC4", "#F8766D"), ylim=c(0,0.4),
        
        xlab="", ylab="Incidence of seed predator")



#### Construct the plant phylogenetic tree 

library (V.PhyloMaker)

Taxonomy2019 <- read.csv("Taxonomy2019_50.csv")

Tree2019 <- phylo.maker(sp.list=Taxonomy2019, tree=GBOTB.extended, nodes=nodes.info.1, scenarios="S3")


#### run pgls

library(ape)

row.names(Cross_species2019_50) <- Cross_species2019_50$Species

Tree2019 <- Tree2019$scenario.3

str(Cross_species2019_50)

plot(Tree2019)

Cross_species2019_50 <- Cross_species2019_50[ , colSums(is.na(Cross_species2019_50)) < nrow(Cross_species2019_50)]  # Remove rows with NA only

library(geiger)

name.check(Cross_species2019_50, Tree2019) # Compares taxa in data and tree

Tree2019$node.label <- NULL

library(caper)

Comparative2019 <- comparative.data(Tree2019, Cross_species2019_50, Species, vcv=TRUE, vcv.dim=2) # combine phylogenies with datasets


# weighted by seed mass

Cross.mod2019 <- pgls(log(sqrt(Incidence_of_seed_predator)+1) ~ Fruit_type + Site + Seedmass.mean, 
                      
                      data=Comparative2019, lambda="ML")

# model diagnostic plot

par(mfrow=c(2,2))

plot(Cross.mod2019) # basically ok, although fitted value plot was not well

par(mfrow=c(1,1))

# model results

summary(Cross.mod2019)

anova(Cross.mod2019)  # lambda = 0.87 (part of Table 3)


# unweighted by seed mass

Cross.mod2019.unweighted <- pgls(log(sqrt(Incidence_of_seed_predator.unweighted)+1) ~ Fruit_type + Site + Seedmass.mean, 
                                 
                                 data=Comparative2019, lambda="ML")

# model diagnostic plot

par(mfrow=c(2,2))

plot(Cross.mod2019.unweighted) # basically ok, although fitted value plot was not well

par(mfrow=c(1,1))

# model results

summary(Cross.mod2019.unweighted)

anova(Cross.mod2019.unweighted) # lambda = 0.83 (part of Table S2)



##### 2.2.2 Incidence of seed predator in 2020 #####

#### Visualization using the violin Plot 


# weighted by seed mass (Figure 3b)

vioplot(Cross_species2020_50$Incidence_of_seed_predator ~ Cross_species2020_50$Site,
        
        col=c("#00BFC4", "#F8766D"), ylim=c(0,0.5), xlab="", ylab="")

axis(side=1, labels=c('CBS', 'XSBN'), at=c(1, 2))


# unweighted by seed mass (Figure S3b)

vioplot(Cross_species2020_50$Incidence_of_seed_predator.unweighted ~ Cross_species2020_50$Site,
        
        col=c("#00BFC4", "#F8766D"), ylim=c(0,0.4), xlab="", ylab="")

axis(side=1, labels=c('CBS', 'XSBN'), at=c(1, 2))



#### Construct the plant phylogenetic tree

Taxonomy2020 <- read.csv("Taxonomy2020_50.csv")

Tree2020 <- phylo.maker(sp.list=Taxonomy2020, tree=GBOTB.extended, nodes=nodes.info.1, scenarios="S3")


#### run pgls

row.names(Cross_species2020_50) <- Cross_species2020_50$Species

Tree2020 <- Tree2020$scenario.3

str(Cross_species2020_50)

plot(Tree2020)

name.check(Cross_species2020_50, Tree2020)

Tree2020$node.label<- NULL

Comparative2020 <- comparative.data(Tree2020, Cross_species2020_50, Species, vcv=TRUE, vcv.dim=2)


# weighted by seed mass

Cross.mod2020 <- pgls(log(sqrt(Incidence_of_seed_predator)+1) ~ Fruit_type + Site + Seedmass.mean, 
                      
                      data=Comparative2020, lambda="ML")

# model diagnostic plot

par(mfrow=c(2,2))

plot(Cross.mod2020) # bacially ok

par(mfrow=c(1,1))

# model results

summary(Cross.mod2020)

anova(Cross.mod2020)  # lambda = 0.49 (part of Table 3)


# unweighted by seed mass

Cross.mod2020.unweighted <- pgls(log(sqrt(Incidence_of_seed_predator.unweighted)+1) ~ Fruit_type + Site + Seedmass.mean , 
                                 
                                 data=Comparative2020, lambda ="ML")

# model diagnostic plot

par(mfrow=c(2,2))

plot(Cross.mod2020.unweighted) # ok

par(mfrow=c(1,1))

# model results

summary(Cross.mod2020.unweighted)

anova(Cross.mod2020.unweighted)  # lambda = 0.55 (part of Table S2)



##### 2.2.3 Seed predation rate in 2020 #####

#### Visualization using the violin Plot (Figure S2)

# weighted by seed mass

# par(mfrow=c(1,2), mar=c(2,5,2,0), oma=c(1,0,1,2), bty="l")

vioplot(Cross_species2020_50$Seed_predation_rate ~ Cross_species2020_50$Site,
        
        col=c("#00BFC4", "#F8766D"), ylim=c(0,0.7), xlab="", ylab="Seed predation rate")

axis(side=1, labels=c('CBS', 'XSBN'), at=c(1, 2))


# unweighted by seedmass

vioplot(Cross_species2020_50$Seed_predation_rate.unweighted ~ Cross_species2020_50$Site,
        
        col=c("#00BFC4", "#F8766D"), ylim=c(0,0.7), xlab="", ylab="")





#### run pgls

# weighted by seed mass

Cross.mod2020.Spr <- pgls(log(sqrt(Seed_predation_rate)+1) ~ Fruit_type + Site + Seedmass.mean, 
                          
                          data=Comparative2020, lambda="ML")

# model diagnostic plot

par(mfrow=c(2,2))

plot(Cross.mod2020.Spr) # basically ok

par(mfrow=c(1,1))

# model results

summary(Cross.mod2020.Spr)

anova(Cross.mod2020.Spr) # lambda = 0.64 (part of Table 3)


# unweighted by seed mass

Cross.mod2020.Spr.unweighted <- pgls(log(sqrt(Seed_predation_rate.unweighted)+1) ~ Fruit_type + Site + Seedmass.mean, 
                                     
                                     data=Comparative2020, lambda="ML")

# model diagnostic plot

par(mfrow=c(2,2))

plot(Cross.mod2020.Spr.unweighted) # basically ok

par(mfrow=c(1,1))

# model results

summary(Cross.mod2020.Spr.unweighted)

anova(Cross.mod2020.Spr.unweighted)  # lambda = 0.66 (part of Table S2)



####  *Simulated data of the Korean Pine Pinus koraiensis ####

#### Visualization using the violin Plot (Figure S4)

# weighted by seed mass

# par(mfrow=c(1,2), mar=c(2,5,2,0), oma=c(1,0,1,2), bty="l")

vioplot(Cross_species2020_50.pine$Seed_predation_rate ~ Cross_species2020_50.pine$Site,
        
        col=c("#00BFC4", "#F8766D"), ylim=c(0,0.7), xlab="", ylab="Seed predation rate")


# unweighted by seed mass

vioplot(Cross_species2020_50.pine$Seed_predation_rate.unweighted ~ Cross_species2020_50.pine$Site,
        
        col=c("#00BFC4", "#F8766D"), ylim=c(0,0.7), xlab="", ylab="")

axis(side=1, labels=c('CBS', 'XSBN'), at=c(1, 2))



#### Construct  plant phylogenetic tree

# Add the taxonomy of the Korean Pine Pinus koraiensis in 'Taxonomy2020'

Taxonomy2020.pine <- rbind(Taxonomy2020, 
                           data.frame("Species" = "Pinus_koraiensis", "Genus" = "Pinus", "Family"="Pinaceae")
                          )

str(Taxonomy2020.pine)

Tree2020.pine <- phylo.maker(sp.list=Taxonomy2020.pine, tree=GBOTB.extended, nodes=nodes.info.1, scenarios="S3")

####  run pgls 

Tree2020.pine <- Tree2020.pine$scenario.3

plot(Tree2020.pine)

row.names(Cross_species2020_50.pine) <- Cross_species2020_50.pine$Species

name.check(Tree2020.pine, Cross_species2020_50.pine)

Tree2020.pine$node.label <- NULL

Comparative2020.pine <- comparative.data(Tree2020.pine, Cross_species2020_50.pine, Species, vcv=TRUE, vcv.dim=2)


# weighted by seed mass

Cross.mod2020.pine <- pgls(log(sqrt(Seed_predation_rate)+1) ~ Fruit_type + Site + Seedmass.mean,
                           
                           data=Comparative2020.pine, lambda = "ML")

# model diagnostic plot

par(mfrow=c(2,2))

plot(Cross.mod2020.pine) # basically ok

par(mfrow=c(1,1))

# model results

summary(Cross.mod2020.pine)

anova(Cross.mod2020.pine) # lambda = 0.80 (part of Table 3)



# unweighted by seed mass

Cross.mod2020.pine.unweighted <- pgls(log(sqrt(Seed_predation_rate.unweighted)+1) ~ Fruit_type + Site + Seedmass.mean, 
                                      
                                      data=Comparative2020.pine, lambda="ML")

# model diagnostic plot

par(mfrow=c(2,2))

plot(Cross.mod2020.pine.unweighted) # basically ok

par(mfrow=c(1,1))

# model results

summary(Cross.mod2020.pine.unweighted)

anova(Cross.mod2020.pine.unweighted) # lambda = 0.82 (part of Table S2)

###################################################### END #######################################################