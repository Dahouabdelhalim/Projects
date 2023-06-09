#Continued Decline Models and Figures
#Pycnopodia rangewide assessment
#Dr. Sarah Gravem, 2021


library(pscl) # for zero inflated models
library(lmtest) # for tukeys
library(tidyverse) # for data manip
library(ggplot2) # for plotting
library(car) # for Anova function
library(parameters) # for model parameters
library(plotrix)#calculate std.error in summary tables


setwd("C:/Users/hamiltsa/Desktop/Pycnopodia/Manuscript/FinalDataAndScripts")
Master <- read.csv("MasterPycno_ToShare.csv")
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

#select  data for continued decline figure 

ContDecline <- Master %>% 
  filter(!is.na(density_m2)) %>% # remove nas in density
  filter(density_m2 < 13) %>% #Dropped one crazy density outlier (13.5m-2 in Salish Sea)
  filter(DepthBin == "shallow")  %>%
  filter(PopulationPhase == "Current") %>% # select current data only (2017-2020)
  filter(event_year != "2020") %>%  #Dropped 2020 -data  thin and misleading bc of outliers
  filter(source != "CCIRA_Dive_BC") %>% #CCIRA data were inflating current data, with no analog in historical 
  mutate(totalpycnos = round(totalpycnos, 0)) %>%  # this rounds total pycnos to an integer (0 sig digits). need for model
  mutate (LogArea = log10(area + 1)) # add column for log area plus one.  need for model
View(ContDecline)


#__________________
# Continued decline FIGURE
 
#Label columns from the summary data
Dens_Cur <- ContDecline$density_m2
Inc_Cur <- ContDecline$pres_abs
Reg_Cur <-as.factor(ContDecline$RegionFull)
Phase_Cur <-as.factor(ContDecline$PopulationPhase)
Date_Cur <-  as.Date(ContDecline$date, format = "%m/%d/%Y")
Pycnos_Cur <- ContDecline$totalpycnos
LogArea_Cur <- ContDecline$LogArea

#change order of levels for graphing
Phase_Cur<-factor(Phase_Cur, levels=c("Historic", "Decline", "Current"))

Reg_Cur<-factor(Reg_Cur, levels=c("Aleutians", "west Gulf of Alaska" ,"east Gulf of Alaska","southeast Alaska" , 
                                "British Columbia*", "Salish Sea", 
                                "Washington*" , "Oregon",
                                "northern California" , "central California" , "southern California", "Baja California" ))

#colors for regions
cols12 <-   c("#9E9E9E","#F06292", "#9C27B0", 
              "#3F51B5", "#2196F3", "#00BCD4", 
              "#4CAF50", "#CDDC39", "#FDD835", 
              "#F99800", "#F44336", "#795548")


#Continued Decline Density

Line_ShalDens_Current <- ggplot(data = ContDecline, aes(x=Date_Cur, y=Dens_Cur, color = Reg_Cur)) +
  geom_point(shape =1, size = 3)+
  geom_smooth(method='lm', se = FALSE)+
  theme_classic()+
  scale_color_manual(values= cols12, name="Region")+
  xlab("Date")+ 
  ylab(expression(paste("Density of" , italic(" Pycnopodia"), "(m" ^ "-2", ")"))) + 
  theme(aspect.ratio = 3/3)+
  theme(axis.title = element_text(size=12)) + 
  theme(axis.text = element_text(size=10)) + 
  theme(legend.text = element_text(size=10)) + 
  theme(plot.title = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 45 , hjust = 1))
plot(Line_ShalDens_Current)
ggsave("ContinuedDecline.png")
ggsave("ContinuedDecline.eps")


# model
#keep only regions with decent data (not flat-lined because those models are too zero-inflated)
ContDecline_ModData <- ContDecline %>%
  filter (RegionFull != "Aleutians", RegionFull != "west Gulf of Alaska", 
          RegionFull != "Washington*", RegionFull != "Oregon", 
          RegionFull != "northern California", RegionFull != "central California", 
          RegionFull != "southern California", RegionFull != "Baja California")

#Re-run vector labels with trimmed data 

Dens_CurM <- ContDecline_ModData$density_m2
Inc_CurM <- ContDecline_ModData$pres_abs
Reg_CurM <-as.factor(ContDecline_ModData$RegionFull)
Phase_CurM <-as.factor(ContDecline_ModData$PopulationPhase)
Date_CurM <-  as.Date(ContDecline_ModData$date, format = "%m/%d/%Y")
Pycnos_CurM <- ContDecline_ModData$totalpycnos
LogArea_CurM <- ContDecline_ModData$LogArea

#how many zeroes? 
100*sum(Pycnos_CurM == 0)/nrow(ContDecline_ModData)
#[1] 73.61111

#THIS ZEROINFL DIDN"T WORK
M_Dens_Current <- zeroinfl(Pycnos_CurM ~ Date_CurM * Reg_CurM | ## Predictor for the Poisson process (regular part)
                        Date_CurM * Reg_CurM , ## Predictor for the Bernoulli process; (zeroes)
                         dist = 'poisson',
                         offset = LogArea_CurM, 
                         data = ContDecline_ModData)
#threw errors. I think this is because there are SO many zeroes. I have droppd everything I can think of and it still doesnt work/. 
# Warning message:
#   In value[[3L]](cond) :
#   system is computationally singular: reciprocal condition number = 1.62342e-17FALSE

#THIS GLM WORKED> NOT PERFECT BC OF ZEROINFL BUT BEST I CAN DO>
M2_Dens_Current <- glm(Pycnos_CurM ~ Date_CurM * Reg_CurM , 
                      family = poisson(),
                           offset = LogArea_CurM, 
                           data = ContDecline_ModData)
summary(M2_Dens_Current)
plot(M2_Dens_Current)
#doesn't look great but kinda the best we can do in this situation
anova(M2_Dens_Current)
Anova(M2_Dens_Current)

#Summary table for model
ContDeclineSum<-ContDecline_ModData %>%
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
View(ContDeclineSum)

