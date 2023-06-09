#Density and Incidence Models and Figures
#Pycnopodia range wide assessment
#created by Dr. Sarah Gravem (2021)

require(pscl) # for zero inflated models
library(tidyverse)
library(ggplot2)
library(ggpattern)
library(plotrix) # for plot textures
library(car) # for Anaova function
library(parameters) # for model parameters
library(plotrix)#calculate std.error in summary tables
library(car) #Anova() function for models
library(readxl)

setwd("C:/Users/hamiltsa/Desktop/Pycnopodia/Manuscript/FinalDataAndScripts")
Master<- read.csv("MasterPycno_ToShare.csv")

#change all instances of central and northern coastal British Columbia to British Columbia* and all washinton outer coast to Washinton*
Master$RegionFull <- str_replace_all(Master$RegionFull, "central and northern coastal British Columbia", "British Columbia*")
Master$RegionFull <- str_replace_all(Master$RegionFull, "Washington outer coast", "Washington*")

#change factor orders because it's nice later
RegionFull_order=c("Aleutians", "west Gulf of Alaska" ,"east Gulf of Alaska","southeast Alaska" , 
         "British Columbia*", "Salish Sea", 
         "Washington*" , "Oregon",
         "northern California" , "central California" , "southern California", "Baja California" )
Master$RegionFull <- factor (Master$RegionFull, levels = RegionFull_order)

#change PopPhase order because it's nice later
PopPhase_order=c("Historic", "Decline" ,"Current")
Master$PopulationPhase <- factor (Master$PopulationPhase, levels = PopPhase_order)

#Create Shallow and Deep Datasets
Shallow <- Master %>%
  filter(DepthBin == 'shallow')
Deep <- Master %>%
  filter(DepthBin == 'deep')

#FULL SUMMARY TABLES FOR ALL DATA
#will need later for bar charts too
ShalSum_ByReg<-Shallow %>%
  #filter(PopulationPhase != "Decline") %>% #Can drop  DECLINE phase for easier table building. don't drop for figure (lots of gaps here that will prevent models, plus uninteresting)
  filter(source != "CCIRA_Dive_BC") %>% #CCIRA data were inflating current data, with no analog in historical 
  complete(RegionFull ,  nesting(PopulationPhase), fill = list(density_m2 = NA, pres_abs = NA)) %>% #this fills in missing combinations with NA
  group_by(RegionFull, PopulationPhase) %>% # this tells which columns to group by
  summarize (NDens = sum(!is.na(density_m2)), # this is the list of statistics to calculate. na.rm = TRUE means ignore NAs
             Density = mean(density_m2, na.rm=TRUE), 
             SDDens = sd(density_m2, na.rm=TRUE),
             SEDens = std.error(density_m2, na.rm=TRUE),
             NInc = sum(!is.na(pres_abs)),
             Incidence = mean(pres_abs),
             SDInc = sd(pres_abs),
             SEInc = std.error(pres_abs))

ShalSum<-Shallow %>%
  filter(PopulationPhase != "Decline") %>% #Dropped DECLINE phase (lots of gaps here that will prevent models, plus uninteresting)
  filter(source != "CCIRA_Dive_BC") %>% #CCIRA data were inflating current data, with no analog in historical
  complete(PopulationPhase, fill = list(density_m2 = NA, pres_abs = NA)) %>% #this fills in missing combinations with NA
  group_by(PopulationPhase) %>% # this tells which columns to group by
  summarize (NDens = sum(!is.na(density_m2)),
             Density = mean(density_m2, na.rm=TRUE), # this is the list of statistics to calculate. na.rm = TRUE means ignore NAs
             SDDens = sd(density_m2, na.rm=TRUE),
             SEDens = std.error(density_m2, na.rm=TRUE),
             NInc = sum(!is.na(pres_abs)),
             Incidence = mean(pres_abs),
             SDInc = sd(pres_abs),
             SEInc = std.error(pres_abs))

DeepSum_ByReg<-Deep %>%
  filter(PopulationPhase != "Decline") %>% #Dropped DECLINE phase (lots of gaps here that will prevent models, plus uninteresting)
  complete(RegionFull ,  nesting(PopulationPhase), fill = list(density_m2 = NA, pres_abs = NA)) %>% #this fills in missing combinations with NA
  group_by(RegionFull, PopulationPhase) %>% # this tells which columns to group by
  summarize (NDens = sum(!is.na(density_m2)),
             Density = mean(density_m2, na.rm=TRUE), # this is the list of statistics to calculate. na.rm = TRUE means ignore NAs
             SDDens = sd(density_m2, na.rm=TRUE),
             SEDens = std.error(density_m2, na.rm=TRUE),
             NInc = sum(!is.na(pres_abs)),
             Incidence = mean(pres_abs),
             SDInc = sd(pres_abs),
             SEInc = std.error(pres_abs))

DeepSum<-Deep %>%
  filter(PopulationPhase != "Decline") %>% #Dropped DECLINE phase (lots of gaps here that will prevent models, plus uninteresting)
  complete(PopulationPhase, fill = list(density_m2 = NA, pres_abs = NA)) %>% #this fills in missing combinations with NA
  group_by(PopulationPhase) %>% # this tells which columns to group by
  summarize (NDens = sum(!is.na(density_m2)),
             Density = mean(density_m2, na.rm=TRUE), # this is the list of statistics to calculate. na.rm = TRUE means ignore NAs
             SDDens = sd(density_m2, na.rm=TRUE),
             SEDens = std.error(density_m2, na.rm=TRUE),
             NInc = sum(!is.na(pres_abs)),
             Incidence = mean(pres_abs),
             SDInc = sd(pres_abs),
             SEInc = std.error(pres_abs))

GlobalSum<-Master %>%
  filter(PopulationPhase != "Decline") %>% #Dropped DECLINE phase (lots of gaps here that will prevent models, plus uninteresting)
  filter(source != "CCIRA_Dive_BC") %>% #CCIRA data were inflating current data, with no analog in historical
  complete(PopulationPhase, fill = list(density_m2 = NA, pres_abs = NA)) %>% #this fills in missing combinations with NA
  group_by(PopulationPhase) %>% # this tells which columns to group by
  summarize (NDens = sum(!is.na(density_m2)),
             Density = mean(density_m2, na.rm=TRUE), # this is the list of statistics to calculate. na.rm = TRUE means ignore NAs
             SDDens = sd(density_m2, na.rm=TRUE),
             SEDens = std.error(density_m2, na.rm=TRUE),
             NInc = sum(!is.na(pres_abs)),
             Incidence = mean(pres_abs),
             SDInc = sd(pres_abs),
             SEInc = std.error(pres_abs))

# #Save files
# library("writexl")
# write_xlsx(ShalSum_ByReg, 'Summary_Shall_ByPhaseRegion.xlsx')
# write_xlsx(ShalSum, 'Summary_Shall_ByPhase.xlsx')
# write_xlsx(DeepSum_ByReg, 'Summary_Deep_ByPhaseRegion.xlsx')
# write_xlsx(DeepSum, 'Summary_Deep_ByPhase.xlsx')
# write_xlsx(GlobalSum, 'Summary_Global_ByPhase.xlsx')
# 

#__________________
# FIGURES
 
#Label columns from the summary data
Dens_SS <- ShalSum_ByReg$Density
Inc_SS <- ShalSum_ByReg$Incidence
Reg_SS <-as.factor(ShalSum_ByReg$RegionFull)
Phase_SS <-as.factor(ShalSum_ByReg$PopulationPhase)
#Label columns from the summary data
Dens_DS <- DeepSum_ByReg$Density
Inc_DS <- DeepSum_ByReg$Incidence
Reg_DS <-as.factor(DeepSum_ByReg$RegionFull)
Phase_DS <-as.factor(DeepSum_ByReg$PopulationPhase)

#change order of levels for graphing
Phase_SS<-factor(Phase_SS, levels=c("Historic", "Decline", "Current"))

Reg_SS<-factor(Reg_SS, levels=c("Aleutians", "west Gulf of Alaska" ,"east Gulf of Alaska","southeast Alaska" , 
                                "British Columbia*", "Salish Sea", 
                                "Washington*" , "Oregon",
                                "northern California" , "central California" , "southern California", "Baja California" ))
Phase_DS<-factor(Phase_DS, levels=c("Historic", "Decline", "Current"))

Reg_DS<-factor(Reg_DS, levels=c("Aleutians", "west Gulf of Alaska" ,"east Gulf of Alaska","southeast Alaska" , 
                                "British Columbia*", "Salish Sea", 
                                "Washington*" , "Oregon",
                                "northern California" , "central California" , "southern California", "Baja California" ))

#colors for regions
cols12 <-   c("#9E9E9E","#F06292", "#9C27B0", 
              "#3F51B5", "#2196F3", "#00BCD4", 
              "#4CAF50", "#CDDC39", "#FDD835", 
              "#F99800", "#F44336", "#795548")


#Shallow Density Figure, by Region and Phase

Bar_ShalDens_PhaseReg <- ggplot(data = ShalSum_ByReg, aes(x=Reg_SS, y=Dens_SS, fill = Reg_SS, pattern = Phase_SS)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(preserve = 'single'),
    color = "black", 
    pattern_fill = "darkgray",
    pattern_angle = 45,
    pattern_density = 0.05,
    pattern_spacing = 0.02,
    pattern_key_scale_factor = 0.9) + 
  geom_errorbar( aes(x=Reg_SS, ymin=Density-SEDens, ymax=Density+SEDens), position = position_dodge(preserve = 'single')) +
  scale_pattern_manual(values = c(Historic = "none", Decline = "stripe", Current = "crosshatch")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = FALSE) +
  theme_classic()+
  labs (title = "a", pattern = "Phase")+
  scale_fill_manual(values= cols12, name="Region")+
  xlab("Region")+ 
  ylab(expression(paste("Density of" , italic(" Pycnopodia"), "(m" ^ "-2", ")"))) + 
  theme(aspect.ratio = 3/5)+
  theme(axis.title = element_text(size=12)) + 
  theme(axis.text = element_text(size=10)) + 
  theme(legend.text = element_text(size=10)) + 
  theme(plot.title = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 45 , hjust = 1))
plot(Bar_ShalDens_PhaseReg)
ggsave("Bar_ShalDens_PhaseReg.png")
ggsave("Bar_ShalDens_PhaseReg.eps")

#Deep Density Figure, by Region and Phase
Bar_DeepDens_PhaseReg <- ggplot(data = DeepSum_ByReg, aes(x=Reg_DS, y=Dens_DS, fill = Reg_DS, pattern = Phase_DS)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(preserve = 'single'),
    color = "black", 
    pattern_fill = "darkgray",
    pattern_angle = 45,
    pattern_density = 0.05,
    pattern_spacing = 0.02,
    pattern_key_scale_factor = 0.9) + 
  geom_errorbar( aes(x=Reg_DS, ymin=Density-SEDens, ymax=Density+SEDens), position = position_dodge(preserve = 'single')) +
  scale_pattern_manual(values = c(Historic = "none", Decline = "stripe", Current = "crosshatch")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = FALSE) +
  theme_classic()+
  labs (title = "c", pattern = "Phase")+
  scale_fill_manual(values= cols12, name="Region")+
  xlab("Region")+ 
  ylab(expression(paste("Density of" , italic(" Pycnopodia"), "(m" ^ "-2", ")"))) + 
  theme(aspect.ratio = 3/5)+
  theme(axis.title = element_text(size=12)) + 
  theme(axis.text = element_text(size=10)) + 
  theme(legend.text = element_text(size=10)) + 
  theme(plot.title = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 45 , hjust = 1))
plot(Bar_DeepDens_PhaseReg)
ggsave("Bar_DeepDens_PhaseReg.png")
ggsave("Bar_DeepDens_PhaseReg.eps")

#Shallow Incidence Figure, by Region and Phase
Bar_ShalInc_PhaseReg <- ggplot(data = ShalSum_ByReg, aes(x=Reg_SS, y=Inc_SS, fill = Reg_SS, pattern = Phase_SS)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(preserve = 'single'),
    color = "black", 
    pattern_fill = "darkgray",
    pattern_angle = 45,
    pattern_density = 0.05,
    pattern_spacing = 0.02,
    pattern_key_scale_factor = 0.9) + 
  geom_errorbar( aes(x=Reg_SS, ymin=Incidence-SEInc, ymax=Incidence+SEInc), position = position_dodge(preserve = 'single')) +
  scale_pattern_manual(values = c(Historic = "none", Decline = "stripe", Current = "crosshatch")) +
  guides(pattern = FALSE,
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme_classic()+
  labs (title = "b", pattern = "Phase")+
  scale_fill_manual(values= cols12, name="Region")+
  xlab("Region")+ 
  ylab(expression(paste("Occurrence of" , italic(" Pycnopodia"), " (survey" ^ "-1", ")"))) + 
  theme(aspect.ratio = 3/5)+
  theme(axis.title = element_text(size=12)) + 
  theme(axis.text = element_text(size=10)) + 
  theme(legend.text = element_text(size=10)) + 
  theme(plot.title = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 45 , hjust = 1))
plot(Bar_ShalInc_PhaseReg)
ggsave("Bar_ShalInc_PhaseReg.png")
ggsave("Bar_ShalInc_PhaseReg.eps")

#Deep Incidence Figure, by Region and Phase
Bar_DeepInc_PhaseReg <- ggplot(data = DeepSum_ByReg, aes(x=Reg_DS, y=Inc_DS, fill = Reg_DS, pattern = Phase_DS)) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(preserve = 'single'),
    color = "black", 
    pattern_fill = "darkgray",
    pattern_angle = 45,
    pattern_density = 0.05,
    pattern_spacing = 0.02,
    pattern_key_scale_factor = 0.9) + 
  geom_errorbar( aes(x=Reg_DS, ymin=Incidence-SEInc, ymax=Incidence+SEInc), position = position_dodge(preserve = 'single')) +
  scale_pattern_manual(values = c(Historic = "none", Decline = "stripe", Current = "crosshatch")) +
  guides(pattern = FALSE,
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme_classic()+
  labs (title = "d", pattern = "Phase")+
  scale_fill_manual(values= cols12, name="Region")+
  xlab("Region")+ 
  ylab(expression(paste("Incidence of" , italic(" Pycnopodia"), " (survey" ^ "-1", ")"))) + 
  theme(aspect.ratio = 3/5)+
  theme(axis.title = element_text(size=12)) + 
  theme(axis.text = element_text(size=10)) + 
  theme(legend.text = element_text(size=10)) + 
  theme(plot.title = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 45 , hjust = 1))
plot(Bar_DeepInc_PhaseReg)
ggsave("Bar_DeepInc_PhaseReg.png")
ggsave("Bar_DeepInc_PhaseReg.eps")


#__________________
#Create Density  Datasets for stats
Dens <- Master %>% 
  filter(!is.na(density_m2)) %>%
  filter(source != "CCIRA_Dive_BC") %>% #CCIRA data were inflating current data, with no analog in historical 
  filter(density_m2 < 13) %>% #Dropped one crazy density outlier (13.5m-2 in Salish Sea)
  filter(event_year != "2020") %>% #Dropped 2020 -data  thin and misleading bc of outliers
  filter(PopulationPhase != "Decline") %>% #Dropped DECLINE phase (lots of gaps here that will prevent models, plus uninteresting)
  mutate(totalpycnos = round(totalpycnos, 0)) %>%  # this rounds total pycnos to an integer (0 sig digits)
  mutate (LogArea = log10(area + 1)) # add column for log area plus one
#View(Dens)

#Create Shallow and Deep Density Datasets
Dens_Shall<- Dens %>%
  filter(DepthBin == "shallow")  %>%
  filter (Region_Shall != "SouthernCalifornia", Region_Shall != "Baja") # had to exclude from models bc can't model all zeroes. 
  
#View(Dens_Shall)

Dens_Deep<- Dens %>%
  filter(DepthBin == "deep")
#View(Dens_Deep)


#SHALLOW DENSITY MODEL
#This model is zero-inflated. Had to use R. (Resids vs. fits in JMP were terrible, and no phase diffs within any region)
#Poisson dist. Use pYcno counts with area offset. 
#Model (zero-inflated glm with poisson dist and log link
#Lumped western alaska and Aleutians using Reg_Shall column
#see older density code for models that didn't work

# Figure out the number of cells in each combination, and df's for difference lumping schemes
DensShallNs<-Dens_Shall %>%
  complete(Region_Shall ,  nesting(PopulationPhase), fill = list(density_m2 = NA)) %>% #this fills in missing combinations with NA
  group_by(Region_Shall, PopulationPhase) %>% # this tells which columns to group by
  summarize (NDens = sum(!is.na(density_m2)), # this is the list of statistics to calculate. na.rm = TRUE means ignore NAs
             Density = mean(density_m2, na.rm=TRUE), 
             SDDens = sd(density_m2, na.rm=TRUE),
             SEDens = std.error(density_m2, na.rm=TRUE))

#Name shallow density columns (_SD)
LogArea_SD<-Dens_Shall$LogArea
PopPhase_SD<-Dens_Shall$PopulationPhase
Reg_SD<-Dens_Shall$Region_Shall # simplifies regions now that I need lots more df's per cell with zeroinfle model
Dens_SD <- Dens_Shall$density_m2
Pycnos_SD <-Dens_Shall$totalpycnos

# how many zeroes? 
100*sum(Pycnos_SD == 0)/nrow(Dens_Shall)
[1] 34.19573
#34% zeroes! need zero inflation model

library(pscl)
M_Dens_Shall <- zeroinfl(Pycnos_SD ~ Reg_SD * PopPhase_SD | ## Predictor for the Poisson process (regular part)
                 Reg_SD * PopPhase_SD, ## Predictor for the Bernoulli process; (zeroes)
               dist = 'poisson',
               offset = LogArea_SD, 
               data = Dens_Shall)

summary(M_Dens_Shall)

N <- Dens_Shall %>%
  summarize (NDens = sum(!is.na(density_m2)))
#Run ANOVA on Shallow density data
Anova(M_Dens_Shall)
#post-hoc tests for the above ANOVA
library(lsmeans)
lsmeans(M_Dens_Shall, pairwise~PopPhase_SD|Reg_SD, adjust="tukey")
lsmeans(M_Dens_Shall, pairwise~Reg_SD|PopPhase_SD, adjust="tukey")

#check for overdispersion
E2 <- resid(M_Dens_Shall, type = "pearson")
N  <- nrow(Dens_Shall)
p  <- length(coef(M_Dens_Shall))   
sum(E2^2) / (N - p)
[1] 4.940174
#overdispersed but a ton better than other models

#__________________
#DeepDensity

#label factors Density Deep model _DD
LogArea_DD<-Dens_Deep$LogArea
PopPhase_DD<-Dens_Deep$PopulationPhase
Dens_DD <- Dens_Deep$density_m2
Pycnos_DD <-Dens_Deep$totalpycnos

#Summary table for modeled data
DensDeepNs<-Dens_Deep %>%
  complete(RegionFull , nesting(PopulationPhase), fill = list(density_m2 = NA)) %>% #this fills in missing combinations with NA
  group_by(RegionFull, PopulationPhase) %>% # this tells which columns to group by
  summarize (NDens = sum(!is.na(density_m2)), # this is the list of statistics to calculate. na.rm = TRUE means ignore NAs
             Density = mean(density_m2, na.rm=TRUE), 
             SDDens = sd(density_m2, na.rm=TRUE),
             SEDens = std.error(density_m2, na.rm=TRUE))
View(DensDeepNs)
DensDeepNs

#zerores
100*sum(Pycnos_DD == 0)/nrow(Dens_Deep)
#[1] 64.60752
#64% zeroes! 
library(pscl)
M_Dens_Deep <- zeroinfl(Pycnos_DD ~ PopPhase_DD | ## Predictor for the Poisson process
                 PopPhase_DD, ## Predictor for the Bernoulli 0's process;
               dist = 'poisson',
               offset = LogArea_DD, 
               data = Dens_Deep)

#overdispersion test
E2 <- resid(M_Dens_Deep, type = "pearson")
N  <- nrow(Dens_Deep)
p  <- length(coef(M_Dens_Deep))   
sum(E2^2) / (N - p)
#[1] 8.903779
#Overdispersed but better than before
Anova(M_Dens_Deep)
summary(M_Dens_Deep)
N <- Dens_Deep %>%
  summarize (NDens = sum(!is.na(density_m2)))

#__________________

#Create Incidence Datasets for stats
Inc <- Master %>% 
  filter(PopulationPhase != "Decline") %>% #Dropped DECLINE phase (lots of gaps here that will prevent models, plus uninteresting)
  mutate(totalpycnos = round(totalpycnos, 0)) %>%  # this rounds total pycnos to an integer (0 sig digits)
  mutate (LogArea = log10(area + 1)) # add column for log area plus one
#View(Inc)

#Create Shallow  Incidence Datasets
Inc_Shall<- Inc %>%
  filter(DepthBin == "shallow") %>%
  filter (RegionFull != "southern California", RegionFull != "Baja California") # had to exclude from models bc can't model all zeroes. 

#__________________
#Shallow Incidence MODEL 
#Label variables _SI
Area_SI<-Inc_Shall$area
LogArea_SI<-Inc_Shall$LogArea
PopPhase_SI<-Inc_Shall$PopulationPhase
Reg_SI<-Inc_Shall$Region_Shall # grouped PNW and W Alaska
Inc_SI <- Inc_Shall$pres_abs
Pycnos_SI <-Inc_Shall$totalpycnos

#Summary table for modeled data
IncShallNs<-Inc_Shall %>%
  complete(Region_Shall , nesting(PopulationPhase), fill = list(pres_abs = NA)) %>% #this fills in missing combinations with NA
  group_by(Region_Shall, PopulationPhase) %>% # this tells which columns to group by
  summarize (NInc = sum(!is.na(pres_abs)),
             Incidence = mean(pres_abs),
             SDInc = sd(pres_abs),
             SEInc = std.error(pres_abs))
View(IncShallNs)


#Model
M_Inc_Shall <- glm(Inc_SI ~ Area_SI + PopPhase_SI * Reg_SI , family = binomial(link = "logit"),
    data = Inc_Shall)
summary(M_Inc_Shall)
Anova(M_Inc_Shall)
N <- Inc_Shall %>%
  summarize (NDens = sum(!is.na(density_m2)))

#overdispersion test
E2 <- resid(M_Inc_Shall, type = "pearson")
N  <- nrow(Inc_Shall)
p  <- length(coef(M_Inc_Shall))   
sum(E2^2) / (N - p)
[1] 0.1662725
#underdispersed?
#Post-hoc tests for ANOVA
lsmeans(M_Inc_Shall, pairwise~PopPhase_SI|Reg_SI, adjust="tukey")
$lsmeans
lsmeans(M_Inc_Shall, pairwise~Reg_SI|PopPhase_SI, adjust="tukey")

#__________________
#Deep  Incidence MODEL
#unbalanced 1s and 0s : 2,500 vs. 12,500

Inc_Deep<- Inc %>%
  filter(DepthBin == "deep") 
#View(Inc_Deep)

#Label variables _DI
LogArea_DI<-Inc_Deep$LogArea
Area_DI<-Inc_Deep$area
PopPhase_DI<-Inc_Deep$PopulationPhase
Reg_DI<-Inc_Deep$Region_Deep
Inc_DI <- Inc_Deep$pres_abs
Pycnos_DI <-Inc_Deep$totalpycnos

#Summary table for modeled data
IncDeepNs<-Inc_Deep %>%
  complete(Region_Deep , nesting(PopulationPhase), fill = list(pres_abs = NA)) %>% #this fills in missing combinations with NA
  group_by(Region_Deep, PopulationPhase) %>% # this tells which columns to group by
  summarize (NInc = sum(!is.na(pres_abs)),
             Incidence = mean(pres_abs),
             SDInc = sd(pres_abs),
             SEInc = std.error(pres_abs))
M_Inc_Deep <- glm(Inc_DI ~ Area_DI + PopPhase_DI * Reg_DI, family = binomial(link = "logit"),
                   data = Inc_Deep)
summary(M_Inc_Deep)
Anova(M_Inc_Deep)
plot(M_Inc_Deep)
#Looks ok

#overdispersion test DHARMa
library(DHARMa)
sim_resids_DI <- simulateResiduals(M_Inc_Deep, refit=T) 
testDispersion(sim_resids_DI, type = 'bootstrap')
plot(sim_resids_DI)
#Post hoc tests
lsmeans(M_Inc_Deep, pairwise~PopPhase_DI|Reg_DI, adjust="tukey")
lsmeans(M_Inc_Deep, pairwise~Reg_DI|PopPhase_DI, adjust="tukey")

