# Add authors, Long-term trends in gastropod abundance and biodiversity: disentangling effects of press versus pulse disturbances
# Global Ecology and Biogeography, contact details 
# Last edited: 16 August 2021

#Supporting Information - R Code used to conduct analyses including:

#(1) temporal trends in temperature,
#(2) GLMMS to evaluate temporal trends in gastropod abundance, point biodiversity, and gamma biodiversity,
#(3) GLMMS to evaluate responses of gastropods to disturbance and global warming,
#(4) GLMMs to evaluate resistance of gastropod abundance and biodiversity to hurricanes,
#(5) variation partitioning to determine unique and shared variation in
#    gastropod abundance, biodiversity, and composition explained by hurricane 
#    and temperature effects, and
#(6) code to determine predicted values based on negative binomial models
#    as well as the confidence intervals associated with those predicted #values.



# Origianl Gastropod data are archived publicly in accordance with NSF guidelines at
# https://luq.lter.network/data/luqmetadata107 and
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-luq.107.9996737)




# Versions of Packages used for analysis
# MASS Version 7.3.54
# car Version 3.0.11
# lme4 Version 1.1.26
# nlme Version 3.1.152
# vegan Version 2.5.7
# picante Version 1.8.2
# foreign Version 0.8.71
# ggplot2 Version 3.3.2



#(1) temporal trends in temperature

# Variable definitions for "Survey_average.txt"
#Year = calander year
# Hurricane = Hurricane Identity after which data were collected (Hugo, Georges, or Maria)
#Airport_max_temp Daily average ambient temperature for each year
#Max_temp = Daily average understory temperature for each year
# TAH = number of years after last major hurricane (years since Hugo, Georges, or Maria)
# Aulalt = Mean survey abundance for Austroselenites alticola for all 40 points combined each year
# Alcalt = Mean survey abundance for Alcadia alta for all 40 points combined each year
# Alcstr = Mean survey abundance for Alcadia striata for all 40 points combined each year
# Carcar = Mean survey abundance for Caracolus caracolla for all 40 points combined each year
# Carmar = Mean survey abundance for Caracolus marginella for all 40 points combined each year
# Cepsqu = Mean survey abundance for Cepolis squamosa for all 40 points combined each year
# Gaenig = Mean survey abundance for Gaeotis nigrolineata for all 40 points combined each year
# Lamgra = Mean survey abundance for Lamellaxis gracilis for all 40 points combined each year
# Megcro = Mean survey abundance for Megalomastoma croceum for all 40 points combined each year
# Nentri = Mean survey abundance for Nenia tridens for all 40 points combined each year
# Olegla = Mean survey abundance for Oleacina glabra for all 40 points combined each year
# Olepla = Mean survey abundance for Oleacina playa for all 40 points combined each year
# Obeter = Mean survey abundance for Obeliscus terabraster for all 40 points combined each year
# Polacu = Mean survey abundance for Polydontes acutangula for all 40 points combined each year
# Plapor = Mean survey abundance for Platysuccinea portoricensis for all 40 points combined each year
# Suboct = Mean survey abundance for Subulina octona for all 40 points combined each year
# Vagocc = Mean survey abundance for Vaginulus occidentalis for all 40 points combined each year
# Totabu = Mean survey abundance for Total gastropod abundance for all 40 points combined each year
# Rich = Mean survey species richness for all 40 points combined each year
# Diversity = Mean survey Hill-transformed Shannon diversity for all 40 points combined each year
# Even = Mean survey Hill-transformed Camargo evenness for all 40 points combined each year
# Dom = Mean survey Hill-transformed Berger-Parker dominance for all 40 points combined each year

###############################################
#
#   Temporal Trends in Temperature
#
#
#

#Data set with annual gastropod and environmental data
Gastropods <- read.table("Survey_average.txt", sep = '\\t', header = T) 


#Detection of linear temporal trends in average daily maximum understory temperature
model901 =  glm(Max_temp ~ Year, data = Gastropods)
summary(model901)


#Detection of linear temporal trends in average daily maximum ambient temperature
model902 =  glm(Airport_max_temp ~ Year, data = Gastropods)
summary(model902)





#(2) Generalized linear mixed effects models and general linear mixed-effect models
#    to evaluate temporal trends in gastropod abundance, point biodiversity, gamma biodiversity

# Variable definitions for "Long-term Point data.txt"
# Year = calander year
# Point = Point number at which data were collected, 40 unique points surveyed each year
# Hurricane = Hurricane Identity after which data were collected (Hugo, Georges, or Maria)
# Airport_max_temp Daily average ambient temperature for each year
# Max_temp = Daily average understory temperature for each year
# TAH = number of years after last major hurricane (years since Hugo, Georges, or Maria)
# Aulalt = Mean survey abundance for Austroselenites alticola at each point each year
# Alcalt = Mean survey abundance for Alcadia alta at each point each year
# Alcstr = Mean survey abundance for Alcadia striata at each point each year
# Carcar = Mean survey abundance for Caracolus caracolla at each point each year
# Carmar = Mean survey abundance for Caracolus marginella at each point each year
# Cepsqu = Mean survey abundance for Cepolis squamosa at each point each year
# Gaenig = Mean survey abundance for Gaeotis nigrolineata at each point each year
# Lamgra = Mean survey abundance for Lamellaxis gracilis at each point each year
# Megcro = Mean survey abundance for Megalomastoma croceum at each point each year
# Nentri = Mean survey abundance for Nenia tridens at each point each year
# Olegla = Mean survey abundance for Oleacina glabra at each point each year
# Olepla = Mean survey abundance for Oleacina playa at each point each year
# Obeter = Mean survey abundance for Obeliscus terabraster at each point each year
# Polacu = Mean survey abundance for Polydontes acutangula at each point each year
# Plapor = Mean survey abundance for Platysuccinea portoricensis at each point each year
# Suboct = Mean survey abundance for Subulina octona at each point each year
# Vagocc = Mean survey abundance for Vaginulus occidentalis at each point each year
# Totabu = Mean survey abundance for Total gastropod abundance at each point each year
# Rich = Mean survey species richness at each point each year
# Diversity = Mean survey Hill-transformed Shannon diversity at each point each year
# Even = Mean survey Hill-transformed Camargo evenness at each point each year
# Dom = Mean survey Hill-transformed Berger-Parker dominance at each point each year

###############################################
#
#      Temporal Trends in Gastropod Abundance
#    Point Biodiversity, and Gamma Biodiversity
#
#
#

library(MASS)
library(car)
library(lme4)
library(nlme)

Gastropods_point <- read.table("Long-term Point data.txt", sep = '\\t', header = T)
Gastropods_point$Hurricane = factor(Gastropods_point$Hurricane)
Gastropods_point$Point = factor(Gastropods_point$Point)

Gastropods <- read.table("Survey_average.txt", sep = '\\t', header = T) 
Gastropods$Hurricane = factor(Gastropods$Hurricane)

#Generalized linear mixed-effects model evaluating temporal trends in Austroselenites alticola abundance
model401 =  glmer.nb(Aulalt ~ Year + (Year|Point), data = Gastropods_point)
summary(model401)

#Generalized linear mixed-effects model evaluating temporal trends in Alcadia altla abundance
model402 =  glmer.nb(Alcalt ~ Year + (Year|Point), data = Gastropods_point)
summary(model402)

#Generalized linear mixed-effects model evaluating temporal trends in Alcadia striata abundance
model403 =  glmer.nb(Alcstr ~ Year + (Year|Point), data = Gastropods_point)
summary(model403)

#Generalized linear mixed-effects model evaluating temporal trends in Caracolus caracolla abundance
model404 =  glmer.nb(Carcar ~ Year + (Year|Point), data = Gastropods_point)
summary(model404)

#Generalized linear mixed-effects model evaluating temporal trends in Caracollus marginella abundance
model405 =  glmer.nb(Carmar ~ Year + (Year|Point), data = Gastropods_point)
summary(model405)


#Generalized linear mixed-effects model evaluating temporal trends in Cepolis squamosa abundance
model406 =  glmer.nb(Cepsqu ~ Year + (Year|Point), data = Gastropods_point)
summary(model406)


#Generalized linear mixed-effects model evaluating temporal trends in Gaeotis nigrolineata abundance
model407 =  glmer.nb(Gaenig ~ Year + (Year|Point), data = Gastropods_point)
summary(model407)

#Generalized linear mixed-effects model evaluating temporal trends in Lamellaxis gracilis abundance
model408 =  glmer.nb(Lamgra ~ Year + (Year|Point), data = Gastropods_point)
summary(model408)

#Generalized linear mixed-effects model evaluating temporal trends in Megalomastoma croceum abundance
model409 =  glmer.nb(Megcro ~ Year + (Year|Point), data = Gastropods_point)
summary(model409)


#Generalized linear mixed-effects model evaluating temporal trends in Nenia tridens abundance
model410 =  glmer.nb(Nentri ~ Year + (Year|Point), data = Gastropods_point)
summary(model410)

#Generalized linear mixed-effects model evaluating temporal trends in Oleacina glabra abundance
model411 =  glmer.nb(Olegla ~ Year + (Year|Point), data = Gastropods_point)
summary(model411)

#Generalized linear mixed-effects model evaluating temporal trends in Oleacina playa abundance
model412 =  glmer.nb(Olepla ~ Year + (Year|Point), data = Gastropods_point)
summary(model412)

#Generalized linear mixed-effects model evaluating temporal trends in Obeliscus terabraster abundance
model413 =  glmer.nb(Obeter ~ Year + (Year|Point), data = Gastropods_point)
summary(model413)

#Generalized linear mixed-effects model evaluating temporal trends in Polydontes acutangula abundance
model414 =  glmer.nb(Polacu ~ Year + (Year|Point), data = Gastropods_point)
summary(model414)

#Generalized linear mixed-effects model evaluating temporal trends in Platysuccinea portoricensis abundance
model415 =  glmer.nb(Plapor ~ Year + (Year|Point), data = Gastropods_point)
summary(model415)

#Generalized linear mixed-effects model evaluating temporal trends in Subulina octona abundance
model416 =  glmer.nb(Suboct ~ Year + (Year|Point), data = Gastropods_point)
summary(model416)

#Generalized linear mixed-effects model evaluating temporal trends in Vaginulls(Diplosolenodes) occidentalis abundance
model417 =  glmer.nb(Vagocc ~ Year + (Year|Point), data = Gastropods_point)
summary(model417)


#Generalized linear mixed-effects model evaluating temporal trends in total gastropod abundance
model418 =  glmer.nb(Totabu ~ Year + (Year|Point), data = Gastropods_point)
summary(model418)


#General linear mixed-effects model evaluating temporal trends in gastropod Point species ricness 
model421 =  lme(Rich ~ Year, random = ~Year|Point, data = Gastropods_point)
summary(model421)

#General linear mixed-effects model evaluating temporal trends in gastropod Point Shannon diversity 
model422 =  lme(Diversity ~ Year, random = ~Year|Point, data = Gastropods_point)
summary(model422)

#General linear mixed-effects model evaluating temporal trends in gastropod Point Camargo evenness
model423 =  lme(Even ~ Year, random = ~Year|Point, data = Gastropods_point)
summary(model423)

#General linear mixed-effects model evaluating temporal trends in gastropod Point Berger-Parker dominance
model424 =  lme(Dom ~ Year, random = ~Year|Point, data = Gastropods_point)
summary(model424)


#General linear mixed-effects model evaluating temporal trends in gastropod gamma species ricness 
model225 =  lme(Rich ~ Year, random = ~1|Year, data = Gastropods)
summary(model225)

#General linear mixed-effects model evaluating temporal trends in gastropod gamma Shannon diversity 
model226 =  lme(Diversity ~ Year, random = ~1|Year, data = Gastropods)
summary(model226)

#General linear mixed-effects model evaluating temporal trends in gastropod gamma Camargo evenness
model227 =  lme(Even ~ Year, random = ~1|Year, data = Gastropods)
summary(model227)

#General linear mixed-effects model evaluating temporal trends in gastropod gamma Berger-Parker dominance
model228 =  lme(Dom ~ Year, random = ~1|Year, data = Gastropods)
summary(model228)





#(3) GLMMS to evaluate responses of gastropods to disturbance and global warming,

#Set of analyses evaluating responses of gastropod abundance and biodiversity to
#interannual variation in hurricane-induced disturbances and temperature

# Variable definitions for "Long-term Point data.txt"
# Year = calander year
# Point = Point number at which data were collected, 40 unique points surveyed each year
# Hurricane = Hurricane Identity after which data were collected (Hugo, Georges, or Maria)
# Airport_max_temp Daily average ambient temperature for each year
# Max_temp = Daily average understory temperature for each year
# TAH = number of years after last major hurricane (years since Hugo, Georges, or Maria)
# Aulalt = Mean survey abundance for Austroselenites alticola at each point each year
# Alcalt = Mean survey abundance for Alcadia alta at each point each year
# Alcstr = Mean survey abundance for Alcadia striata at each point each year
# Carcar = Mean survey abundance for Caracolus caracolla at each point each year
# Carmar = Mean survey abundance for Caracolus marginella at each point each year
# Cepsqu = Mean survey abundance for Cepolis squamosa at each point each year
# Gaenig = Mean survey abundance for Gaeotis nigrolineata at each point each year
# Lamgra = Mean survey abundance for Lamellaxis gracilis at each point each year
# Megcro = Mean survey abundance for Megalomastoma croceum at each point each year
# Nentri = Mean survey abundance for Nenia tridens at each point each year
# Olegla = Mean survey abundance for Oleacina glabra at each point each year
# Olepla = Mean survey abundance for Oleacina playa at each point each year
# Obeter = Mean survey abundance for Obeliscus terabraster at each point each year
# Polacu = Mean survey abundance for Polydontes acutangula at each point each year
# Plapor = Mean survey abundance for Platysuccinea portoricensis at each point each year
# Suboct = Mean survey abundance for Subulina octona at each point each year
# Vagocc = Mean survey abundance for Vaginulus occidentalis at each point each year
# Totabu = Mean survey abundance for Total gastropod abundance at each point each year
# Rich = Mean survey species richness at each point each year
# Diversity = Mean survey Hill-transformed Shannon diversity at each point each year
# Even = Mean survey Hill-transformed Camargo evenness at each point each year
# Dom = Mean survey Hill-transformed Berger-Parker dominance at each point each year


###############################################
#
#   Effects of disturbance and global
#   warming on gastropod abundance and
#            biodiversity
#
#

library(MASS)
library(car)
library(lme4)
library(nlme)


Gastropods_point <- read.table("Long-term Point data.txt", sep = '\\t', header = T)
Gastropods_point$Hurricane = factor(Gastropods_point$Hurricane)
Gastropods_point$Point = factor(Gastropods_point$Point)

Gastropods <- read.table("Survey_average.txt", sep = '\\t', header = T) 
Gastropods$Hurricane = factor(Gastropods$Hurricane)


#Analysis for Austroselenites alticola abundance
model101 =  glmer.nb(Aulalt ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model101, type="II")


#Analysis for Alcadia striata abundance
model102 =  glmer.nb(Alcalt ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model102, type="II")


#Analysis for  Alcadia alta abundance
model103 =  glmer.nb(Alcstr ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model103, type="II")


#Analysis for Caracolus caracollus abundance
model104 =  glmer.nb(Carcar ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model104, type="II")


#Analysis for Caracolus marginella abundance
model105 =  glmer.nb(Carmar ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model105, type="II")


#Analysis for Cepolis squamosa abundance
model106 =  glmer.nb(Cepsqu ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model106, type="II")


#Analysis for Gaeotis nigrolineata abundance
model107 =  glmer.nb(Gaenig ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model107, type="II")


#Analysis for Lamellaxis gracilis abundance
model108 =  glmer.nb(Lamgra ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model108, type="II")


#Analysis for Megalomastoma croceum abundance
model109 =  glmer.nb(Megcro ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model109, type="II")


#Analysis for Nenia tridens abundance
model110 =  glmer.nb(Nentri ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model110, type="II")


#Analysis for Oleacina glabra abundance
model111 =  glmer.nb(Olegla ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model111, type="II")


#Analysis for Oleacina playa abundance
model112 =  glmer.nb(Olepla ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model112, type="II")


#Analysis for Obeter terebraster abundance
model113 =  glmer.nb(Obeter ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model113, type="II")


#Analysis for Polydontes acutangula abundance
model114 =  glmer.nb(Polacu ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model114, type="II")


#Analysis for Platysuccinea portoricensis abundance
model115 =  glmer.nb(Plapor ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model115, type="II")


#Analysis for Subulina octona abundance
model116 =  glmer.nb(Suboct ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model116, type="II")


#Analysis for Diplosolenodes occidentalis abundance
model117 =  glmer.nb(Vagocc ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model117, type="II")


#Analysis for Total gastropod abundance
model118 =  glmer.nb(Totabu ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane + (Year|Point), data = Gastropods_point)
Anova(model118, type="II")


#Analysis for Point Species richness
model121 =  lme(Rich ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane, random = ~Year|Point, data = Gastropods_point)
Anova(model121, type="II")


#Analysis for Point Shannon diversity
model122 =  lme(Diversity ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane, random = ~Year|Point, data = Gastropods_point)
Anova(model122, type="II")


#Analysis for Point Camargo evenness
model123 =  lme(Even ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane, random = ~Year|Point, data = Gastropods_point)
Anova(model123, type="II")


#Analysis for Point Berger-Parker dominance
model124 =  lme(Dom ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane, random = ~Year|Point, data = Gastropods_point)
Anova(model124, type="II")


#Analysis for Gamma Species richness
model125 =  lme(Rich ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane, random = ~Year|Point, data = Gastropods)
Anova(model125, type="II")


#Analysis for Gamma Shannon diversity
model126 =  lme(Diversity ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane, random = ~Year|Point, data = Gastropods)
Anova(model126, type="II")


#Analysis for Gamma Camargo evenness
model127 =  lme(Even ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane, random = ~Year|Point, data = Gastropods)
Anova(model127, type="II")


#Analysis for Gamma Berger-Parker dominance
model128 =  lme(Dom ~ Hurricane + Airport_max_temp + Max_temp + TAH + TAH*Hurricane, random = ~Year|Point, data = Gastropods)
Anova(model128, type="II")





#(4) GLMMs to evaluate resistance of gastropod abundance and biodiversity to hurricanes,


# Variable definitions for "Hurricane_Resistance_data.txt"
# Point = Point number at which data were collected, 40 unique points surveyed each year
# Disturbance = year before or year after a hurricane
# Hurricane = Hurricane Identity after which data were collected (Hugo, Georges, or Maria)
# Aulalt = Mean survey abundance for Austroselenites alticola at each point the year before or year after a hurricane
# Alcalt = Mean survey abundance for Alcadia alta at each point the year before or year after a hurricane
# Alcstr = Mean survey abundance for Alcadia striata at each point the year before or year after a hurricane
# Carcar = Mean survey abundance for Caracolus caracolla at each point the year before or year after a hurricane
# Carmar = Mean survey abundance for Caracolus marginella at each point the year before or year after a hurricane
# Cepsqu = Mean survey abundance for Cepolis squamosa at each point the year before or year after a hurricane
# Gaenig = Mean survey abundance for Gaeotis nigrolineata at each point the year before or year after a hurricane
# Lamgra = Mean survey abundance for Lamellaxis gracilis at each point the year before or year after a hurricane
# Megcro = Mean survey abundance for Megalomastoma croceum at each point the year before or year after a hurricane
# Nentri = Mean survey abundance for Nenia tridens at each point the year before or year after a hurricane
# Olegla = Mean survey abundance for Oleacina glabra at each point the year before or year after a hurricane
# Olepla = Mean survey abundance for Oleacina playa at each point the year before or year after a hurricane
# Obeter = Mean survey abundance for Obeliscus terabraster at each point the year before or year after a hurricane
# Polacu = Mean survey abundance for Polydontes acutangula at each point the year before or year after a hurricane
# Plapor = Mean survey abundance for Platysuccinea portoricensis at each point the year before or year after a hurricane
# Suboct = Mean survey abundance for Subulina octona at each point the year before or year after a hurricane
# Vagocc = Mean survey abundance for Vaginulus occidentalis at each point the year before or year after a hurricane
# Totabu = Mean survey abundance for Total gastropod abundance at each point the year before or year after a hurricane
# Rich = Mean survey species richness at each point the year before or year after a hurricane
# Diversity = Mean survey Hill-transformed Shannon diversity at each point the year before or year after a hurricane
# Even = Mean survey Hill-transformed Camargo evenness at each point the year before or year after a hurricane
# Dom = Mean survey Hill-transformed Berger-Parker dominance at each point the year before or year after a hurricane


###############################################
#
#   Evaluation of Gastropod resistance
#            to hurricanes   
#
#

library(MASS)
library(car)
library(lme4)
library(nlme)

Resist <- read.table("Hurricane_Resistance_data.txt", sep = '\\t', header = T)
Resist$Hurricane = factor(Resist$Hurricane)
Resist$Disturbance = factor(Resist$Disturbance)


#Analysis for Austroselenites alticola abundance
model001 =  glmer.nb(Aulalt ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model001, type="II")


#Analysis for Alcadia striata abundance
model002 =  glmer.nb(Alcalt ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model002, type="II")


#Analysis for  Alcadia alta abundance
model003 =  glmer.nb(Alcstr ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model003, type="II")


#Analysis for Caracolus caracollus abundance
model004 =  glmer.nb(Carcar ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model004, type="II")


#Analysis for Caracolus marginella abundance
model005 =  glmer.nb(Carmar ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model005, type="II")


#Analysis for Cepolis squamosa abundance
model006 =  glmer.nb(Cepsqu ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model006, type="II")

#Analysis for Gaeotis nigrolineata abundance
model007 =  glmer.nb(Gaenig ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model007, type="II")


#Analysis for Lamellaxis gracilis abundance
model008 =  glmer.nb(Lamgra ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model008, type="II")


#Analysis for Megalomastoma croceum abundance
model009 =  glmer.nb(Megcro ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model009, type="II")


#Analysis for Nenia tridens abundance
model010 =  glmer.nb(Nentri ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model010, type="II")


#Analysis for Oleacina glabra abundance
model011 =  glmer.nb(Olegla ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model011, type="II")


#Analysis for Oleacina playa abundance
model012 =  glmer.nb(Olepla ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model012, type="II")

#Analysis for Obeter terebraster abundance
model013 =  glmer.nb(Obeter ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model013, type="II")


#Analysis for Polydontes acutangula abundance
model014 =  glmer.nb(Polacu ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model014, type="II")


#Analysis for Platysuccinea portoricensis abundance
model015 =  glmer.nb(Plapor ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model015, type="II")


#Analysis for Subulina octona abundance
model016 =  glmer.nb(Suboct ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model016, type="II")


#Analysis for Diplosolenodes occidentalis abundance
model017 =  glmer.nb(Vagocc ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model017, type="II")


#Analysis for Total gastropod abundance
model018 =  glmer.nb(Totabu ~ Hurricane*Disturbance + (1|Point), data = Resist)
Anova(model018, type="II")


#Analysis for Point Species richness
model019 =  lme(Rich ~ Hurricane*Disturbance, random = ~1|Point, data = Resist)
Anova(model019, type="II")


#Analysis for Point Shannon diversity
model020 =  lme(Diversity ~ Hurricane*Disturbance, random = ~1|Point, data = Resist)
Anova(model020, type="II")


#Analysis for Point Camargo evenness
model021 =  lme(Even ~ Hurricane*Disturbance, random = ~1|Point, data = Resist)
Anova(model021, type="II")


#Analysis for Point Berger-Parker dominance
model022 =  lme(Dom ~ Hurricane*Disturbance, random = ~1|Point, data = Resist)
Anova(model022, type="II")




#(5) variation partitioning to determine unique and shared variation in
#    gastropod abundance, biodiversity, and composition explained by hurricane 
#    and temperature effects, and

# Each file of abundance or biodiversity for variation partitioning only includes
# a single column of abundance or biodiversity values indicated by the file name.
# Value are for each year from 1993 through 2019 and are based on the combined
# data for all 40 points from the LFDP

# "Composition_props.txt" are the proportional abundances for each of 17 species
# for each year from 1993 through 2019, calculated for data combined from
# all 40 points

# Variable definitions for "Gastropod_explanatory.txt"
# Year = calander year
# Hurricane = Hurricane Identity after which data were collected (Hugo, Georges, or Maria)
# Airport_max_temp Daily average ambient temperature for each year
# Max_temp = Daily average understory temperature for each year
# TAH = number of years after last major hurricane (years since Hugo, Georges, or Maria)


###############################################
#
#   Gastropod - Variation Partioning
#
#
#


library(vegan)
library(picante)

# Data sets with annual abundance values for each species and for all gastropods as a group
Aulalt_abund <- read.table("Aulalt_abund.txt", sep = '\\t', header = T)
Alcalt_abund <- read.table("Alcalt_abund.txt", sep = '\\t', header = T)
Alcstr_abund <- read.table("Alcstr_abund.txt", sep = '\\t', header = T)
Carcar_abund <- read.table("Carcar_abund.txt", sep = '\\t', header = T)
Carmar_abund <- read.table("Carmar_abund.txt", sep = '\\t', header = T)
Cepsqu_abund <- read.table("Cepsqu_abund.txt", sep = '\\t', header = T)
Gaenig_abund <- read.table("Gaenig_abund.txt", sep = '\\t', header = T)
Lamgra_abund <- read.table("Lamgra_abund.txt", sep = '\\t', header = T)
Megcro_abund <- read.table("Megcro_abund.txt", sep = '\\t', header = T)
Nentri_abund <- read.table("Nentri_abund.txt", sep = '\\t', header = T)
Obeter_abund <- read.table("Obeter_abund.txt", sep = '\\t', header = T)
Olegla_abund <- read.table("Olegla_abund.txt", sep = '\\t', header = T)
Olepla_abund <- read.table("Olepla_abund.txt", sep = '\\t', header = T)
Plapor_abund <- read.table("Plapor_abund.txt", sep = '\\t', header = T)
Polacu_abund <- read.table("Polacu_abund.txt", sep = '\\t', header = T)
Suboct_abund <- read.table("Suboct_abund.txt", sep = '\\t', header = T)
Vagocc_abund <- read.table("Vagocc_abund.txt", sep = '\\t', header = T)
Total_abund <- read.table("Total_abund.txt", sep = '\\t', header = T)


# Data sets with annual values for biodiversity metrics
Richness <- read.table("Richness.txt", sep = '\\t', header = T)
Diversity <- read.table("Diversity.txt", sep = '\\t', header = T)
Evenness <- read.table("Evenness.txt", sep = '\\t', header = T)
Dominance <- read.table("Dominance.txt", sep = '\\t', header = T)


# Data set with proportional abundances each year for gastropod species
Composition_props <- read.table("Composition_props.txt", sep = '\\t', header = T)

# Data for 4 biodiversity metrics each year
Biodiversity <- read.table("Biodiversity.txt", sep = '\\t', header = T)

# Data set with explanatory variables
Gastropod_explanatory <- read.table("Gastropod_explanatory.txt", sep = '\\t', header = T)
Gastropod_explanatory$Hurricane = factor(Gastropod_explanatory$Hurricane)

#creation of data for the hurricane partition
Gastropod_hurricane <- model.matrix(~ Hurricane + TAH, Gastropod_explanatory)[,-1]

#creation of data for the ambient temperature partition
Gastropod_ambient <- model.matrix(~ Airport_max_temp, Gastropod_explanatory)[, -1]

#creation of data for the understory temperature partition
Gastropod_understory <- model.matrix(~ Max_temp, Gastropod_explanatory)[,-1]


# Variation partitioning for abundance of Austroselenites alticola

mod1 <- varpart(Aulalt_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory) 
mod1

showvarparts(3, bg=2:4)
plot(mod1, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.aulalt.hurr <- rda(Aulalt_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.aulalt.hurr, step=200, perm.max=200)
RsquareAdj(rda.aulalt.hurr)

rda.aulalt.ambient <- rda(Aulalt_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.aulalt.ambient, step=200, perm.max=200)
RsquareAdj(rda.aulalt.ambient)

rda.aulalt.understory <- rda(Aulalt_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.aulalt.understory, step=200, perm.max=200)
RsquareAdj(rda.aulalt.understory)

# Redundancy analyses to determine significane of whole partitions

rda.aulalt.hurr2 <- rda(Aulalt_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.aulalt.hurr2, step=200, perm.max=200)
RsquareAdj(rda.aulalt.hurr2)

rda.aulalt.ambient2 <- rda(Aulalt_abund ~ as.matrix(Gastropod_ambient))
anova(rda.aulalt.ambient2, step=200, perm.max=200)
RsquareAdj(rda.aulalt.ambient2)

rda.aulalt.under2 <- rda(Aulalt_abund ~ as.matrix(Gastropod_understory))
anova(rda.aulalt.under2, step=200, perm.max=200)
RsquareAdj(rda.aulalt.under2)

# Redundancy analysis to determine significane of full model

rda.aulalt.total <- rda(Aulalt_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.aulalt.total, step=200, perm.max=200)
RsquareAdj(rda.aulalt.total)


# Variation partitioning for abundance of Alcadia alticola

mod2 <- varpart(Alcalt_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod2

showvarparts(3, bg=2:4)
plot(mod2, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.alcalt.hurr <- rda(Alcalt_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.alcalt.hurr, step=200, perm.max=200)
RsquareAdj(rda.alcalt.hurr)

rda.alcalt.ambient <- rda(Alcalt_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.alcalt.ambient, step=200, perm.max=200)
RsquareAdj(rda.alcalt.ambient)

rda.alcalt.understory <- rda(Alcalt_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.alcalt.understory, step=200, perm.max=200)
RsquareAdj(rda.alcalt.understory)

# Redundancy analyses to determine significane of whole partitions

rda.alcalt.hurr2 <- rda(Alcalt_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.alcalt.hurr2, step=200, perm.max=200)
RsquareAdj(rda.alcalt.hurr2)

rda.alcalt.ambient2 <- rda(Alcalt_abund ~ as.matrix(Gastropod_ambient))
anova(rda.alcalt.ambient2, step=200, perm.max=200)
RsquareAdj(rda.alcalt.ambient2)

rda.alcalt.under2 <- rda(Alcalt_abund ~ as.matrix(Gastropod_understory))
anova(rda.alcalt.under2, step=200, perm.max=200)
RsquareAdj(rda.alcalt.under2)

# Redundancy analysis to determine significane of full model

rda.alcalt.total <- rda(Alcalt_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.alcalt.total, step=200, perm.max=200)
RsquareAdj(rda.alcalt.total)



# Variation partitioning for abundance of Alcadia striata

mod3 <- varpart(Alcstr_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod3

showvarparts(3, bg=2:4)
plot(mod3, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.alcstr.hurr <- rda(Alcstr_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.alcstr.hurr, step=200, perm.max=200)
RsquareAdj(rda.alcstr.hurr)

rda.alcstr.ambient <- rda(Alcstr_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.alcstr.ambient, step=200, perm.max=200)
RsquareAdj(rda.alcstr.ambient)

rda.alcstr.understory <- rda(Alcstr_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.alcstr.understory, step=200, perm.max=200)
RsquareAdj(rda.alcstr.understory)

# Redundancy analyses to determine significane of whole partitions

rda.alcstr.hurr2 <- rda(Alcstr_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.alcstr.hurr2, step=200, perm.max=200)
RsquareAdj(rda.alcstr.hurr2)

rda.alcstr.ambient2 <- rda(Alcstr_abund ~ as.matrix(Gastropod_ambient))
anova(rda.alcstr.ambient2, step=200, perm.max=200)
RsquareAdj(rda.alcstr.ambient2)

rda.alcstr.under2 <- rda(Alcstr_abund ~ as.matrix(Gastropod_understory))
anova(rda.alcstr.under2, step=200, perm.max=200)
RsquareAdj(rda.alcstr.under2)

# Redundancy analysis to determine significane of full model

rda.alcstr.total <- rda(Alcstr_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.alcstr.total, step=200, perm.max=200)
RsquareAdj(rda.alcstr.total)



# Variation partitioning for abundance of Caracolus caracolla

mod4 <- varpart(Carcar_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod4

showvarparts(3, bg=2:4)
plot(mod4, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.carcar.hurr <- rda(Carcar_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.carcar.hurr, step=200, perm.max=200)
RsquareAdj(rda.carcar.hurr)

rda.carcar.ambient <- rda(Carcar_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.carcar.ambient, step=200, perm.max=200)
RsquareAdj(rda.carcar.ambient)

rda.carcar.understory <- rda(Carcar_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.carcar.understory, step=200, perm.max=200)
RsquareAdj(rda.carcar.understory)


# Redundancy analyses to determine significane of whole partitions

rda.carcar.hurr2 <- rda(Carcar_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.carcar.hurr2, step=200, perm.max=200)
RsquareAdj(rda.carcar.hurr2)

rda.carcar.ambient2 <- rda(Carcar_abund ~ as.matrix(Gastropod_ambient))
anova(rda.carcar.ambient2, step=200, perm.max=200)
RsquareAdj(rda.carcar.ambient2)

rda.carcar.under2 <- rda(Carcar_abund ~ as.matrix(Gastropod_understory))
anova(rda.carcar.under2, step=200, perm.max=200)
RsquareAdj(rda.carcar.under2)

# Redundancy analysis to determine significane of full model

rda.carcar.total <- rda(Carcar_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.carcar.total, step=200, perm.max=200)
RsquareAdj(rda.carcar.total)



# Variation partitioning for abundance of Caracolus marginella

mod5 <- varpart(Carmar_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod5

showvarparts(3, bg=2:4)
plot(mod5, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.carmar.hurr <- rda(Carmar_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.carmar.hurr, step=200, perm.max=200)
RsquareAdj(rda.carmar.hurr)

rda.carmar.ambient <- rda(Carmar_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.carmar.ambient, step=200, perm.max=200)
RsquareAdj(rda.carmar.ambient)

rda.carmar.understory <- rda(Carmar_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.carmar.understory, step=200, perm.max=200)
RsquareAdj(rda.carmar.understory)


# Redundancy analyses to determine significane of whole partitions

rda.carmar.hurr2 <- rda(Carmar_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.carmar.hurr2, step=200, perm.max=200)
RsquareAdj(rda.carmar.hurr2)

rda.carmar.ambient2 <- rda(Carmar_abund ~ as.matrix(Gastropod_ambient))
anova(rda.carmar.ambient2, step=200, perm.max=200)
RsquareAdj(rda.carmar.ambient2)

rda.carmar.under2 <- rda(Carmar_abund ~ as.matrix(Gastropod_understory))
anova(rda.carmar.under2, step=200, perm.max=200)
RsquareAdj(rda.carmar.under2)

# Redundancy analysis to determine significane of full model

rda.carmar.total <- rda(Carmar_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.carmar.total, step=200, perm.max=200)
RsquareAdj(rda.carmar.total)



# Variation partitioning for abundance of Cepolus squamosa

mod6 <- varpart(Cepsqu_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod6

showvarparts(3, bg=2:4)
plot(mod6, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.cepsqu.hurr <- rda(Cepsqu_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.cepsqu.hurr, step=200, perm.max=200)
RsquareAdj(rda.cepsqu.hurr)

rda.cepsqu.ambient <- rda(Cepsqu_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.cepsqu.ambient, step=200, perm.max=200)
RsquareAdj(rda.cepsqu.ambient)

rda.cepsqu.understory <- rda(Cepsqu_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.cepsqu.understory, step=200, perm.max=200)
RsquareAdj(rda.cepsqu.understory)


# Redundancy analyses to determine significane of whole partitions

rda.cepsqu.hurr2 <- rda(Cepsqu_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.cepsqu.hurr2, step=200, perm.max=200)
RsquareAdj(rda.cepsqu.hurr2)

rda.cepsqu.ambient2 <- rda(Cepsqu_abund ~ as.matrix(Gastropod_ambient))
anova(rda.cepsqu.ambient2, step=200, perm.max=200)
RsquareAdj(rda.cepsqu.ambient2)

rda.cepsqu.under2 <- rda(Cepsqu_abund ~ as.matrix(Gastropod_understory))
anova(rda.cepsqu.under2, step=200, perm.max=200)
RsquareAdj(rda.cepsqu.under2)

# Redundancy analysis to determine significane of full model

rda.cepsqu.total <- rda(Cepsqu_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.cepsqu.total, step=200, perm.max=200)
RsquareAdj(rda.cepsqu.total)



# Variation partitioning for abundance of Gaeotis nigrolineata

mod7 <- varpart(Gaenig_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod7

showvarparts(3, bg=2:4)
plot(mod7, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.gaenig.hurr <- rda(Gaenig_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.gaenig.hurr, step=200, perm.max=200)
RsquareAdj(rda.gaenig.hurr)

rda.gaenig.ambient <- rda(Gaenig_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.gaenig.ambient, step=200, perm.max=200)
RsquareAdj(rda.gaenig.ambient)

rda.gaenig.understory <- rda(Gaenig_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.gaenig.understory, step=200, perm.max=200)
RsquareAdj(rda.gaenig.understory)


# Redundancy analyses to determine significane of whole partitions

rda.gaenig.hurr2 <- rda(Gaenig_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.gaenig.hurr2, step=200, perm.max=200)
RsquareAdj(rda.gaenig.hurr2)

rda.gaenig.ambient2 <- rda(Gaenig_abund ~ as.matrix(Gastropod_ambient))
anova(rda.gaenig.ambient2, step=200, perm.max=200)
RsquareAdj(rda.gaenig.ambient2)

rda.gaenig.under2 <- rda(Gaenig_abund ~ as.matrix(Gastropod_understory))
anova(rda.gaenig.under2, step=200, perm.max=200)
RsquareAdj(rda.gaenig.under2)

# Redundancy analysis to determine significane of full model

rda.gaenig.total <- rda(Gaenig_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.gaenig.total, step=200, perm.max=200)
RsquareAdj(rda.gaenig.total)



# Variation partitioning for abundance of Lamellaxis gracilis

mod8 <- varpart(Lamgra_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod8

showvarparts(3, bg=2:4)
plot(mod8, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.lamgra.hurr <- rda(Lamgra_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.lamgra.hurr, step=200, perm.max=200)
RsquareAdj(rda.lamgra.hurr)

rda.lamgra.ambient <- rda(Lamgra_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.lamgra.ambient, step=200, perm.max=200)
RsquareAdj(rda.lamgra.ambient)

rda.lamgra.understory <- rda(Lamgra_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.lamgra.understory, step=200, perm.max=200)
RsquareAdj(rda.lamgra.understory)


# Redundancy analyses to determine significane of whole partitions

rda.lamgra.hurr2 <- rda(Lamgra_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.lamgra.hurr2, step=200, perm.max=200)
RsquareAdj(rda.lamgra.hurr2)

rda.lamgra.ambient2 <- rda(Lamgra_abund ~ as.matrix(Gastropod_ambient))
anova(rda.lamgra.ambient2, step=200, perm.max=200)
RsquareAdj(rda.lamgra.ambient2)

rda.lamgra.under2 <- rda(Lamgra_abund ~ as.matrix(Gastropod_understory))
anova(rda.lamgra.under2, step=200, perm.max=200)
RsquareAdj(rda.lamgra.under2)

# Redundancy analysis to determine significane of full model

rda.lamgra.total <- rda(Lamgra_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.lamgra.total, step=200, perm.max=200)
RsquareAdj(rda.lamgra.total)




# Variation partitioning for abundance of Megalomastoma croceum

mod9 <- varpart(Megcro_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod9

showvarparts(3, bg=2:4)
plot(mod9, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.megcro.hurr <- rda(Megcro_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.megcro.hurr, step=200, perm.max=200)
RsquareAdj(rda.megcro.hurr)

rda.megcro.ambient <- rda(Megcro_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.megcro.ambient, step=200, perm.max=200)
RsquareAdj(rda.megcro.ambient)

rda.megcro.understory <- rda(Megcro_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.megcro.understory, step=200, perm.max=200)
RsquareAdj(rda.megcro.understory)


# Redundancy analyses to determine significane of whole partitions

rda.megcro.hurr2 <- rda(Megcro_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.megcro.hurr2, step=200, perm.max=200)
RsquareAdj(rda.megcro.hurr2)

rda.megcro.ambient2 <- rda(Megcro_abund ~ as.matrix(Gastropod_ambient))
anova(rda.megcro.ambient2, step=200, perm.max=200)
RsquareAdj(rda.megcro.ambient2)

rda.megcro.under2 <- rda(Megcro_abund ~ as.matrix(Gastropod_understory))
anova(rda.megcro.under2, step=200, perm.max=200)
RsquareAdj(rda.megcro.under2)

# Redundancy analysis to determine significane of full model

rda.megcro.total <- rda(Megcro_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.megcro.total, step=200, perm.max=200)
RsquareAdj(rda.megcro.total)



# Variation partitioning for abundance of Nenia tridens

mod10 <- varpart(Nentri_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod10

showvarparts(3, bg=2:4)
plot(mod10, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.nentri.hurr <- rda(Nentri_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.nentri.hurr, step=200, perm.max=200)
RsquareAdj(rda.nentri.hurr)

rda.nentri.ambient <- rda(Nentri_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.nentri.ambient, step=200, perm.max=200)
RsquareAdj(rda.nentri.ambient)

rda.nentri.understory <- rda(Nentri_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.nentri.understory, step=200, perm.max=200)
RsquareAdj(rda.nentri.understory)


# Redundancy analyses to determine significane of whole partitions

rda.nentri.hurr2 <- rda(Nentri_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.nentri.hurr2, step=200, perm.max=200)
RsquareAdj(rda.nentri.hurr2)

rda.nentri.ambient2 <- rda(Nentri_abund ~ as.matrix(Gastropod_ambient))
anova(rda.nentri.ambient2, step=200, perm.max=200)
RsquareAdj(rda.nentri.ambient2)

rda.nentri.under2 <- rda(Nentri_abund ~ as.matrix(Gastropod_understory))
anova(rda.nentri.under2, step=200, perm.max=200)
RsquareAdj(rda.nentri.under2)

# Redundancy analysis to determine significane of full model

rda.nentri.total <- rda(Nentri_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.nentri.total, step=200, perm.max=200)
RsquareAdj(rda.nentri.total)



# Variation partitioning for abundance of Obeliscus terebraster

mod11 <- varpart(Obeter_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod11

showvarparts(3, bg=2:4)
plot(mod11, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.obeter.hurr <- rda(Obeter_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.obeter.hurr, step=200, perm.max=200)
RsquareAdj(rda.obeter.hurr)

rda.obeter.ambient <- rda(Obeter_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.obeter.ambient, step=200, perm.max=200)
RsquareAdj(rda.obeter.ambient)

rda.obeter.understory <- rda(Obeter_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.obeter.understory, step=200, perm.max=200)
RsquareAdj(rda.obeter.understory)


# Redundancy analyses to determine significane of whole partitions

rda.obeter.hurr2 <- rda(Obeter_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.obeter.hurr2, step=200, perm.max=200)
RsquareAdj(rda.obeter.hurr2)

rda.obeter.ambient2 <- rda(Obeter_abund ~ as.matrix(Gastropod_ambient))
anova(rda.obeter.ambient2, step=200, perm.max=200)
RsquareAdj(rda.obeter.ambient2)

rda.obeter.under2 <- rda(Obeter_abund ~ as.matrix(Gastropod_understory))
anova(rda.obeter.under2, step=200, perm.max=200)
RsquareAdj(rda.obeter.under2)

# Redundancy analysis to determine significane of full model

rda.obeter.total <- rda(Obeter_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.obeter.total, step=200, perm.max=200)
RsquareAdj(rda.obeter.total)



# Variation partitioning for abundance of Oleacina glabra

mod12 <- varpart(Olegla_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod12

showvarparts(3, bg=2:4)
plot(mod12, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.olegla.hurr <- rda(Olegla_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.olegla.hurr, step=200, perm.max=200)
RsquareAdj(rda.olegla.hurr)

rda.olegla.ambient <- rda(Olegla_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.olegla.ambient, step=200, perm.max=200)
RsquareAdj(rda.olegla.ambient)

rda.olegla.understory <- rda(Olegla_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.olegla.understory, step=200, perm.max=200)
RsquareAdj(rda.olegla.understory)


# Redundancy analyses to determine significane of whole partitions

rda.olegla.hurr2 <- rda(Olegla_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.olegla.hurr2, step=200, perm.max=200)
RsquareAdj(rda.olegla.hurr2)

rda.olegla.ambient2 <- rda(Olegla_abund ~ as.matrix(Gastropod_ambient))
anova(rda.olegla.ambient2, step=200, perm.max=200)
RsquareAdj(rda.olegla.ambient2)

rda.olegla.under2 <- rda(Olegla_abund ~ as.matrix(Gastropod_understory))
anova(rda.olegla.under2, step=200, perm.max=200)
RsquareAdj(rda.olegla.under2)

# Redundancy analysis to determine significane of full model

rda.olegla.total <- rda(Olegla_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.olegla.total, step=200, perm.max=200)
RsquareAdj(rda.olegla.total)



# Variation partitioning for abundance of Oleacina playa

mod13 <- varpart(Olepla_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod13

showvarparts(3, bg=2:4)
plot(mod13, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.olepla.hurr <- rda(Olepla_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.olepla.hurr, step=200, perm.max=200)
RsquareAdj(rda.olepla.hurr)

rda.olepla.ambient <- rda(Olepla_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.olepla.ambient, step=200, perm.max=200)
RsquareAdj(rda.olepla.ambient)

rda.olepla.understory <- rda(Olepla_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.olepla.understory, step=200, perm.max=200)
RsquareAdj(rda.olepla.understory)


# Redundancy analyses to determine significane of whole partitions

rda.olepla.hurr2 <- rda(Olepla_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.olepla.hurr2, step=200, perm.max=200)
RsquareAdj(rda.olepla.hurr2)

rda.olepla.ambient2 <- rda(Olepla_abund ~ as.matrix(Gastropod_ambient))
anova(rda.olepla.ambient2, step=200, perm.max=200)
RsquareAdj(rda.olepla.ambient2)

rda.olepla.under2 <- rda(Olepla_abund ~ as.matrix(Gastropod_understory))
anova(rda.olepla.under2, step=200, perm.max=200)
RsquareAdj(rda.olepla.under2)

# Redundancy analysis to determine significane of full model

rda.olepla.total <- rda(Olepla_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.olepla.total, step=200, perm.max=200)
RsquareAdj(rda.olepla.total)




# Variation partitioning for abundance of Polydontes acutangula

mod14 <- varpart(Polacu_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod14

showvarparts(3, bg=2:4)
plot(mod14, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.polacu.hurr <- rda(Polacu_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.polacu.hurr, step=200, perm.max=200)
RsquareAdj(rda.polacu.hurr)

rda.polacu.ambient <- rda(Polacu_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.polacu.ambient, step=200, perm.max=200)
RsquareAdj(rda.polacu.ambient)

rda.polacu.understory <- rda(Polacu_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.polacu.understory, step=200, perm.max=200)
RsquareAdj(rda.polacu.understory)


# Redundancy analyses to determine significane of whole partitions

rda.polacu.hurr2 <- rda(Polacu_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.polacu.hurr2, step=200, perm.max=200)
RsquareAdj(rda.polacu.hurr2)

rda.polacu.ambient2 <- rda(Polacu_abund ~ as.matrix(Gastropod_ambient))
anova(rda.polacu.ambient2, step=200, perm.max=200)
RsquareAdj(rda.polacu.ambient2)

rda.polacu.under2 <- rda(Polacu_abund ~ as.matrix(Gastropod_understory))
anova(rda.polacu.under2, step=200, perm.max=200)
RsquareAdj(rda.polacu.under2)

# Redundancy analysis to determine significane of full model

rda.polacu.total <- rda(Polacu_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.polacu.total, step=200, perm.max=200)
RsquareAdj(rda.polacu.total)




# Variation partitioning for abundance of Platysuccinea portoricensis

mod15 <- varpart(Plapor_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod15

showvarparts(3, bg=2:4)
plot(mod15, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.plapor.hurr <- rda(Plapor_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.plapor.hurr, step=200, perm.max=200)
RsquareAdj(rda.plapor.hurr)

rda.plapor.ambient <- rda(Plapor_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.plapor.ambient, step=200, perm.max=200)
RsquareAdj(rda.plapor.ambient)

rda.plapor.understory <- rda(Plapor_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.plapor.understory, step=200, perm.max=200)
RsquareAdj(rda.plapor.understory)


# Redundancy analyses to determine significane of whole partitions

rda.plapor.hurr2 <- rda(Plapor_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.plapor.hurr2, step=200, perm.max=200)
RsquareAdj(rda.plapor.hurr2)

rda.plapor.ambient2 <- rda(Plapor_abund ~ as.matrix(Gastropod_ambient))
anova(rda.plapor.ambient2, step=200, perm.max=200)
RsquareAdj(rda.plapor.ambient2)

rda.plapor.under2 <- rda(Plapor_abund ~ as.matrix(Gastropod_understory))
anova(rda.plapor.under2, step=200, perm.max=200)
RsquareAdj(rda.plapor.under2)

# Redundancy analysis to determine significane of full model

rda.plapor.total <- rda(Plapor_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.plapor.total, step=200, perm.max=200)
RsquareAdj(rda.plapor.total)




# Variation partitioning for abundance of Subulina octona

mod16 <- varpart(Suboct_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod16

showvarparts(3, bg=2:4)
plot(mod16, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.suboct.hurr <- rda(Suboct_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.suboct.hurr, step=200, perm.max=200)
RsquareAdj(rda.suboct.hurr)

rda.suboct.ambient <- rda(Suboct_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.suboct.ambient, step=200, perm.max=200)
RsquareAdj(rda.suboct.ambient)

rda.suboct.understory <- rda(Suboct_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.suboct.understory, step=200, perm.max=200)
RsquareAdj(rda.suboct.understory)


# Redundancy analyses to determine significane of whole partitions

rda.suboct.hurr2 <- rda(Suboct_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.suboct.hurr2, step=200, perm.max=200)
RsquareAdj(rda.suboct.hurr2)

rda.suboct.ambient2 <- rda(Suboct_abund ~ as.matrix(Gastropod_ambient))
anova(rda.suboct.ambient2, step=200, perm.max=200)
RsquareAdj(rda.suboct.ambient2)

rda.suboct.under2 <- rda(Suboct_abund ~ as.matrix(Gastropod_understory))
anova(rda.suboct.under2, step=200, perm.max=200)
RsquareAdj(rda.suboct.under2)

# Redundancy analysis to determine significane of full model

rda.suboct.total <- rda(Suboct_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.suboct.total, step=200, perm.max=200)
RsquareAdj(rda.suboct.total)




# Variation partitioning for abundance of Diplosolenodes occidentalis

mod17 <- varpart(Vagocc_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod17

showvarparts(3, bg=2:4)
plot(mod17, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.vagocc.hurr <- rda(Vagocc_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.vagocc.hurr, step=200, perm.max=200)
RsquareAdj(rda.vagocc.hurr)

rda.vagocc.ambient <- rda(Vagocc_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.vagocc.ambient, step=200, perm.max=200)
RsquareAdj(rda.vagocc.ambient)

rda.vagocc.understory <- rda(Vagocc_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.vagocc.understory, step=200, perm.max=200)
RsquareAdj(rda.vagocc.understory)


# Redundancy analyses to determine significane of whole partitions

rda.vagocc.hurr2 <- rda(Vagocc_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.vagocc.hurr2, step=200, perm.max=200)
RsquareAdj(rda.vagocc.hurr2)

rda.vagocc.ambient2 <- rda(Vagocc_abund ~ as.matrix(Gastropod_ambient))
anova(rda.vagocc.ambient2, step=200, perm.max=200)
RsquareAdj(rda.vagocc.ambient2)

rda.vagocc.under2 <- rda(Vagocc_abund ~ as.matrix(Gastropod_understory))
anova(rda.vagocc.under2, step=200, perm.max=200)
RsquareAdj(rda.vagocc.under2)

# Redundancy analysis to determine significane of full model

rda.vagocc.total <- rda(Vagocc_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.vagocc.total, step=200, perm.max=200)
RsquareAdj(rda.vagocc.total)




# Variation partitioning for abundance of total gastropod abundance

mod18 <- varpart(Total_abund, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod18

showvarparts(3, bg=2:4)
plot(mod18, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.total.hurr <- rda(Total_abund ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.total.hurr, step=200, perm.max=200)
RsquareAdj(rda.total.hurr)

rda.total.ambient <- rda(Total_abund ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.total.ambient, step=200, perm.max=200)
RsquareAdj(rda.total.ambient)

rda.total.understory <- rda(Total_abund ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.total.understory, step=200, perm.max=200)
RsquareAdj(rda.total.understory)


# Redundancy analyses to determine significane of whole partitions

rda.total.hurr2 <- rda(Total_abund ~ as.matrix(Gastropod_hurricane))
anova(rda.total.hurr2, step=200, perm.max=200)
RsquareAdj(rda.total.hurr2)

rda.total.ambient2 <- rda(Total_abund ~ as.matrix(Gastropod_ambient))
anova(rda.total.ambient2, step=200, perm.max=200)
RsquareAdj(rda.total.ambient2)

rda.total.under2 <- rda(Total_abund ~ as.matrix(Gastropod_understory))
anova(rda.total.under2, step=200, perm.max=200)
RsquareAdj(rda.total.under2)

# Redundancy analysis to determine significane of full model

rda.total.total <- rda(Total_abund ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.total.total, step=200, perm.max=200)
RsquareAdj(rda.total.total)




# Variation partitioning for species richness

mod19 <- varpart(Richness, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod19

showvarparts(3, bg=2:4)
plot(mod19, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.rich.hurr <- rda(Richness ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.rich.hurr, step=200, perm.max=200)
RsquareAdj(rda.rich.hurr)

rda.rich.ambient <- rda(Richness ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.rich.ambient, step=200, perm.max=200)
RsquareAdj(rda.rich.ambient)

rda.rich.understory <- rda(Richness ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.rich.understory, step=200, perm.max=200)
RsquareAdj(rda.rich.understory)


# Redundancy analyses to determine significane of whole partitions

rda.rich.hurr2 <- rda(Richness ~ as.matrix(Gastropod_hurricane))
anova(rda.rich.hurr2, step=200, perm.max=200)
RsquareAdj(rda.rich.hurr2)

rda.rich.ambient2 <- rda(Richness ~ as.matrix(Gastropod_ambient))
anova(rda.rich.ambient2, step=200, perm.max=200)
RsquareAdj(rda.rich.ambient2)

rda.rich.under2 <- rda(Richness ~ as.matrix(Gastropod_understory))
anova(rda.rich.under2, step=200, perm.max=200)
RsquareAdj(rda.rich.under2)

# Redundancy analysis to determine significane of full model

rda.rich.total <- rda(Richness ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.rich.total, step=200, perm.max=200)
RsquareAdj(rda.rich.total)




# Variation partitioning for Shannon diversity

mod20 <- varpart(Diversity, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod20

showvarparts(3, bg=2:4)
plot(mod20, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.diver.hurr <- rda(Diversity ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.diver.hurr, step=200, perm.max=200)
RsquareAdj(rda.diver.hurr)

rda.diver.ambient <- rda(Diversity ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.diver.ambient, step=200, perm.max=200)
RsquareAdj(rda.diver.ambient)

rda.diver.understory <- rda(Diversity ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.diver.understory, step=200, perm.max=200)
RsquareAdj(rda.diver.understory)


# Redundancy analyses to determine significane of whole partitions

rda.diver.hurr2 <- rda(Diversity ~ as.matrix(Gastropod_hurricane))
anova(rda.diver.hurr2, step=200, perm.max=200)
RsquareAdj(rda.diver.hurr2)

rda.diver.ambient2 <- rda(Diversity ~ as.matrix(Gastropod_ambient))
anova(rda.diver.ambient2, step=200, perm.max=200)
RsquareAdj(rda.diver.ambient2)

rda.diver.under2 <- rda(Diversity ~ as.matrix(Gastropod_understory))
anova(rda.diver.under2, step=200, perm.max=200)
RsquareAdj(rda.diver.under2)

# Redundancy analysis to determine significane of full model

rda.diver.total <- rda(Diversity ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.diver.total, step=200, perm.max=200)
RsquareAdj(rda.diver.total)




# Variation partitioning for Camargo evenness

mod21 <- varpart(Evenness, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod21

showvarparts(3, bg=2:4)
plot(mod21, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.even.hurr <- rda(Evenness ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.even.hurr, step=200, perm.max=200)
RsquareAdj(rda.even.hurr)

rda.even.ambient <- rda(Evenness ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.even.ambient, step=200, perm.max=200)
RsquareAdj(rda.even.ambient)

rda.even.understory <- rda(Evenness ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.even.understory, step=200, perm.max=200)
RsquareAdj(rda.even.understory)


# Redundancy analyses to determine significane of whole partitions

rda.even.hurr2 <- rda(Evenness ~ as.matrix(Gastropod_hurricane))
anova(rda.even.hurr2, step=200, perm.max=200)
RsquareAdj(rda.even.hurr2)

rda.even.ambient2 <- rda(Evenness ~ as.matrix(Gastropod_ambient))
anova(rda.even.ambient2, step=200, perm.max=200)
RsquareAdj(rda.even.ambient2)

rda.even.under2 <- rda(Evenness ~ as.matrix(Gastropod_understory))
anova(rda.even.under2, step=200, perm.max=200)
RsquareAdj(rda.even.under2)

# Redundancy analysis to determine significane of full model

rda.even.total <- rda(Evenness ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.even.total, step=200, perm.max=200)
RsquareAdj(rda.even.total)




# Variation partitioning for Berger-Parker dominance

mod22 <- varpart(Dominance, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod22

showvarparts(3, bg=2:4)
plot(mod22, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.dom.hurr <- rda(Dominance ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.dom.hurr, step=200, perm.max=200)
RsquareAdj(rda.dom.hurr)

rda.dom.ambient <- rda(Dominance ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.dom.ambient, step=200, perm.max=200)
RsquareAdj(rda.dom.ambient)

rda.dom.understory <- rda(Dominance ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.dom.understory, step=200, perm.max=200)
RsquareAdj(rda.dom.understory)


# Redundancy analyses to determine significane of whole partitions

rda.dom.hurr2 <- rda(Dominance ~ as.matrix(Gastropod_hurricane))
anova(rda.dom.hurr2, step=200, perm.max=200)
RsquareAdj(rda.dom.hurr2)

rda.dom.ambient2 <- rda(Dominance ~ as.matrix(Gastropod_ambient))
anova(rda.dom.ambient2, step=200, perm.max=200)
RsquareAdj(rda.dom.ambient2)

rda.dom.under2 <- rda(Dominance ~ as.matrix(Gastropod_understory))
anova(rda.dom.under2, step=200, perm.max=200)
RsquareAdj(rda.dom.under2)

# Redundancy analysis to determine significane of full model

rda.dom.total <- rda(Dominance ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.dom.total, step=200, perm.max=200)
RsquareAdj(rda.dom.total)



# Variation partitioning for species composition based on proportional abundances

mod23 <- varpart(Composition_props, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod23

showvarparts(3, bg=2:4)
plot(mod23, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.comp.hurr <- rda(Composition_props ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.comp.hurr, step=200, perm.max=200)
RsquareAdj(rda.comp.hurr)

rda.comp.ambient <- rda(Composition_props ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.comp.ambient, step=200, perm.max=200)
RsquareAdj(rda.comp.ambient)

rda.comp.understory <- rda(Composition_props ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.comp.understory, step=200, perm.max=200)
RsquareAdj(rda.comp.understory)


# Redundancy analyses to determine significane of whole partitions

rda.comp.hurr2 <- rda(Composition_props ~ as.matrix(Gastropod_hurricane))
anova(rda.comp.hurr2, step=200, perm.max=200)
RsquareAdj(rda.comp.hurr2)

rda.comp.ambient2 <- rda(Composition_props ~ as.matrix(Gastropod_ambient))
anova(rda.comp.ambient2, step=200, perm.max=200)
RsquareAdj(rda.comp.ambient2)

rda.comp.under2 <- rda(Composition_props ~ as.matrix(Gastropod_understory))
anova(rda.comp.under2, step=200, perm.max=200)
RsquareAdj(rda.comp.under2)

# Redundancy analysis to determine significane of full model

rda.comp.total <- rda(Composition_props ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.comp.total, step=200, perm.max=200)
RsquareAdj(rda.comp.total)



# Variation partitioning for all 4 biodiversity metrics

mod24 <- varpart(Biodiversity, ~ Hurricane + TAH, ~ Airport_max_temp, ~ Max_temp, data = Gastropod_explanatory)
mod24

showvarparts(3, bg=2:4)
plot(mod24, bg=2:4)

# Redundancy analyses to determine significane of unique partitions

rda.biodiv.hurr <- rda(Biodiversity ~ as.matrix(Gastropod_hurricane) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.biodiv.hurr, step=200, perm.max=200)
RsquareAdj(rda.biodiv.hurr)

rda.biodiv.ambient <- rda(Biodiversity ~ as.matrix(Gastropod_ambient) + Condition(as.matrix(Gastropod_hurricane)) + Condition(as.matrix(Gastropod_understory)))
anova(rda.biodiv.ambient, step=200, perm.max=200)
RsquareAdj(rda.biodiv.ambient)

rda.biodiv.understory <- rda(Biodiversity ~ as.matrix(Gastropod_understory) + Condition(as.matrix(Gastropod_ambient)) + Condition(as.matrix(Gastropod_hurricane)))
anova(rda.biodiv.understory, step=200, perm.max=200)
RsquareAdj(rda.biodiv.understory)


# Redundancy analyses to determine significane of whole partitions

rda.biodiv.hurr2 <- rda(Biodiversity ~ as.matrix(Gastropod_hurricane))
anova(rda.biodiv.hurr2, step=200, perm.max=200)
RsquareAdj(rda.biodiv.hurr2)

rda.biodiv.ambient2 <- rda(Biodiversity ~ as.matrix(Gastropod_ambient))
anova(rda.biodiv.ambient2, step=200, perm.max=200)
RsquareAdj(rda.biodiv.ambient2)

rda.biodiv.under2 <- rda(Biodiversity ~ as.matrix(Gastropod_understory))
anova(rda.biodiv.under2, step=200, perm.max=200)
RsquareAdj(rda.biodiv.under2)

# Redundancy analysis to determine significane of full model

rda.biodiv.total <- rda(Biodiversity ~ as.matrix(Gastropod_hurricane) + as.matrix(Gastropod_ambient) + as.matrix(Gastropod_understory))
anova(rda.biodiv.total, step=200, perm.max=200)
RsquareAdj(rda.biodiv.total)



#(6) code to determine predicted values based on negative binomial models
#    as well as the confidence intervals associated with those predicted #values.

# This script was used to calcuate expected values for negative binomial models to plot regression lines on Figure 4
# so that the regression lines and confidence intervals are based on the same model as that used to determine 
# significance (i.e. a negative binomial error distribution)

###############################################################
#############################################################
#####
#####        Caclulating predicted values for
#####      significant negative binomial models
##
##
##


library(MASS)
library(car)
library(lme4)
library(nlme)
library(foreign)
library(ggplot2)

Gastropods <- read.table("Survey_average.txt", sep = '\\t', header = T)
Gastropods$Year = factor(Gastropods$Year)


ggplot(Gastropods, aes(Aulalt)) + geom_histogram(binwidth = 1)
ggplot(Gastropods, aes(Lamgra)) + geom_histogram(binwidth = 1)
ggplot(Gastropods, aes(Nentri)) + geom_histogram(binwidth = 1)
ggplot(Gastropods, aes(Polacu)) + geom_histogram(binwidth = 1)
ggplot(Gastropods, aes(Totabu)) + geom_histogram(binwidth = 1)


with(Gastropods, tapply(Aulalt, Year, function(x) {sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x)) }))
with(Gastropods, tapply(Lamgra, Year, function(x) {sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x)) }))
with(Gastropods, tapply(Nentri, Year, function(x) {sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x)) }))
with(Gastropods, tapply(Polacu, Year, function(x) {sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x)) }))
with(Gastropods, tapply(Totabu, Year, function(x) {sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x)) }))


# Negative binomial models for relationships of abundance through time (same models as those above)

model201 =  glm.nb(Aulalt ~ Year, data = Gastropods)
model208 =  glm.nb(Lamgra ~ Year, data = Gastropods)
model210 =  glm.nb(Nentri ~ Year, data = Gastropods)
model214 =  glm.nb(Polacu ~ Year, data = Gastropods)
model218 =  glm.nb(Totabu ~ Year, data = Gastropods)

# Creation of variables that include predicted values and confidence intervals for each nb model

est201 <- cbind(Estimate = coef(model201), confint(model201))
est208 <- cbind(Estimate = coef(model208), confint(model208))
est210 <- cbind(Estimate = coef(model210), confint(model210))
est214 <- cbind(Estimate = coef(model214), confint(model214))
est218 <- cbind(Estimate = coef(model218), confint(model218))


# Creation of graph for Austroselenites alticola

newdata1_201 <- data.frame(Year = mean(Gastropods$Year))
newdata1_201$phat <- predict(model201, newdata1_201, type = "response")

newdata2_201 <- data.frame(Year = rep(seq(from = min(Gastropods$Year), to = max(Gastropods$Year), length.out = 100)))
newdata2_201 <- cbind(newdata2_201, predict(model201, newdata2_201, type = "link", se.fit=TRUE))

newdata2_201 <- within(newdata2_201, {
  Aulalt <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


ggplot(newdata2_201, aes(Year, Aulalt)) + geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) + 
  geom_line(size = 1) + labs(x = "Year", y = "Predicted Austrocelenites alticola abundance")


# Creation of graph for Lamellaxis gracilis

newdata1_208 <- data.frame(Year = mean(Gastropods$Year))
newdata1_208$phat <- predict(model208, newdata1_208, type = "response")

newdata2_208 <- data.frame(Year = rep(seq(from = min(Gastropods$Year), to = max(Gastropods$Year), length.out = 100)))
newdata2_208 <- cbind(newdata2_208, predict(model208, newdata2_208, type = "link", se.fit=TRUE))

newdata2_208 <- within(newdata2_208, {
  Lamgra <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


ggplot(newdata2_208, aes(Year, Lamgra)) + geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) + 
  geom_line(size = 1) + labs(x = "Year", y = "Predicted Lamellaxis gracilis abundance")




# Creation of graph for Nenia tridens

newdata1_210 <- data.frame(Year = mean(Gastropods$Year))
newdata1_210$phat <- predict(model210, newdata1_210, type = "response")

newdata2_210 <- data.frame(Year = rep(seq(from = min(Gastropods$Year), to = max(Gastropods$Year), length.out = 100)))
newdata2_210 <- cbind(newdata2_210, predict(model210, newdata2_210, type = "link", se.fit=TRUE))

newdata2_210 <- within(newdata2_210, {
  Nentri <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


ggplot(newdata2_210, aes(Year, Nentri)) + geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) + 
  geom_line(size = 1) + labs(x = "Year", y = "Predicted Nenia tridens abundance")



# Creation of graph for Polydontes acutangula

newdata1_214 <- data.frame(Year = mean(Gastropods$Year))
newdata1_214$phat <- predict(model214, newdata1_214, type = "response")

newdata2_214 <- data.frame(Year = rep(seq(from = min(Gastropods$Year), to = max(Gastropods$Year), length.out = 100)))
newdata2_214 <- cbind(newdata2_214, predict(model214, newdata2_214, type = "link", se.fit=TRUE))

newdata2_214 <- within(newdata2_214, {
  Polacu <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


ggplot(newdata2_214, aes(Year, Polacu)) + geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) + 
  geom_line(size = 1) + labs(x = "Year", y = "Predicted Polydontes acutangula abundance")




# Creation of graph for total gastropod abundances

newdata1_218 <- data.frame(Year = mean(Gastropods$Year))
newdata1_218$phat <- predict(model218, newdata1_218, type = "response")

newdata2_218 <- data.frame(Year = rep(seq(from = min(Gastropods$Year), to = max(Gastropods$Year), length.out = 100)))
newdata2_218 <- cbind(newdata2_218, predict(model218, newdata2_218, type = "link", se.fit=TRUE))

newdata2_218 <- within(newdata2_218, {
  Totabu <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


ggplot(newdata2_218, aes(Year, Totabu)) + geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) + 
  geom_line(size = 1) + labs(x = "Year", y = "Predicted Total abundance abundance")


