#R Code that accompanies the manuscript titled: 
#The relationship between water quality and the proportion of watershed protection: a case study in Québec, Canada 
#Published in Freshwater Science in 2021
#First author: Dalal E.L. Hanna, contact at: dalal.e.hanna@gmail.com 


#Clear code####
rm(list=ls())

#Load relevant libraries####
library(ggplot2)
library(gridExtra)
library(car)
library(lme4)
library (multcomp)
library(MuMIn)
library(dplyr)
library(vegan)
library(corrplot)
library(effsize)


#Set working directory####
#Make sure to adjust to insert your own working directory for data to load properly 
setwd("/Users/dhann028/Desktop/PhD/Chapter2/Chp2Statistics/")


#IQBP DATA ANALYSIS####
#IQBP - Loading and exploring Data####
#load water quality data from government website 
IQBP<-read.csv("BQMA/eaux_canada-master/data/IQPB.csv")
head(IQBP)
colnames(IQBP)[3] <- "StationID"
str(IQBP)


#Explanation of column headers in IQBP
#CF (UFC/100 ml); total coliforms, a unit above 200 UFC/100ml is considered to be too much for recreational activities
#CHL-AA (µg/l); chlorophyll-a, a concentration of above 8.6ug/l is considered to be too high
#COD (mg/l); dissolved organic carbon
#COND (µS/cm), conductivity	
#NH3 (mg/l); ammonia
#NOX (mg/l); nitrites and nitrates
#NTOT (mg/l); total nitrogen, anything above 1mg/L
#PH (pH); pH
#PHEO (µg/l)	
#P-T-PER (mg/l): total phosphorous, anytime above 0.03mg/L is considered to have chronic effects on aquatic organisms
#SS (mg/l): suspended solids, anytime above 13mgl/L is considered too high
#tEMP (°C): temperature
#TURB (UTN): turbidity

IQBP_Summary <-IQBP %>% group_by(StationID) %>% summarize(IQBP_mean = mean(IQBP, na.rm=TRUE),
                                                          # IQBP_SampleSize = n(IQBP),
                                                          IQBP_sd = sd(IQBP,na.rm=TRUE),
                                                          CF_mean = mean(CF, na.rm=TRUE),
                                                          CF_sd = sd(CF,na.rm=TRUE),
                                                          CHl_AA_mean = mean(CHl_AA, na.rm=TRUE),
                                                          CHl_AA_sd = sd(CHl_AA,na.rm=TRUE),
                                                          DOC_mean = mean(COD, na.rm=TRUE),
                                                          DOC_sd = sd(COND,na.rm=TRUE),
                                                          Cond_mean = mean(COD, na.rm=TRUE),
                                                          Cond_sd = sd(COND,na.rm=TRUE),
                                                          NH3_mean = mean(NH3, na.rm=TRUE),
                                                          NH3_sd = sd(NH3,na.rm=TRUE),
                                                          NOX_mean = mean(NOX, na.rm=TRUE),
                                                          NOX_sd = sd(NOX,na.rm=TRUE),
                                                          NTOT_mean = mean(NTOT, na.rm=TRUE),
                                                          NTOT_sd = sd(NTOT,na.rm=TRUE),
                                                          PH_mean = mean(PH, na.rm=TRUE),
                                                          PH_sd = sd(PH,na.rm=TRUE),
                                                          PHEO_mean = mean(PHEO, na.rm=TRUE),
                                                          PHEO_sd = sd(PHEO,na.rm=TRUE),
                                                          PTPER_mean = mean(PTPER, na.rm=TRUE),
                                                          PTPER_sd = sd(PTPER,na.rm=TRUE),
                                                          SS_mean = mean(SS, na.rm=TRUE),
                                                          SS_sd = sd(SS,na.rm=TRUE),
                                                          TEMP_mean = mean(TEMP, na.rm=TRUE),
                                                          TEMP_sd = sd(TEMP,na.rm=TRUE),
                                                          TURB_mean = mean(TURB, na.rm=TRUE),
                                                          TURB_sd = sd(TURB,na.rm=TRUE))

write.csv(IQBP_Summary, file = "IQBP_Summary.csv")

#Dataset used in submission 1
#IQBP_GISvalues<-read.csv("HannaChp2_IQBP_GISData.csv")

#Dataset used for submission 2
IQBP_GISvalues<-read.csv("IQBP_AllData_Oct2020_Edited.csv")


IQBP_FullSummary <- left_join(IQBP_Summary, IQBP_GISvalues, by = "StationID")


#Remove points that did not line up properly with hydrosheds river network, or were located at a border of quebec,
#beyond which we didn't have land use data
table(IQBP_FullSummary$Notes)
IQBP_FullSummary<-subset(IQBP_FullSummary, Notes == "D")
IQBP_FullSummary<-subset(IQBP_FullSummary,  StationsToRemove != "X")

#Drives the dataset down to 183 observations 

#Add columns of different categories of data 
IQBP_FullSummary$ProtectionCategory <- ifelse(IQBP_FullSummary$PercentPro == 0, "0%",
                                              ifelse(IQBP_FullSummary$PercentPro > 0 & IQBP_FullSummary$PercentPro < 13, "1-12%", 
                                                     ifelse(IQBP_FullSummary$PercentPro > 13 & IQBP_FullSummary$PercentPro < 25, "13-24%", 
                                                            ifelse(IQBP_FullSummary$PercentPro > 25 & IQBP_FullSummary$PercentPro < 50, "25-49%", 
                                                                   ifelse(IQBP_FullSummary$PercentPro > 50 & IQBP_FullSummary$PercentPro < 75, "50-74%", "75-100%")))))

IQBP_FullSummary$ProtectionCategorical <- ifelse(IQBP_FullSummary$PercentPro == 0, "No IUCN Category\\nII protection", "Some level of\\nIUCN Category II protection")


table(IQBP_FullSummary$ProtectionCategory)


IQBP_FullSummary$ProtectionCategoryWeighted <- ifelse(IQBP_FullSummary$RelPro == 0, "0%",
                                              ifelse(IQBP_FullSummary$RelPro > 0 & IQBP_FullSummary$RelPro < 13, "1-12%", 
                                                     ifelse(IQBP_FullSummary$RelPro > 13 & IQBP_FullSummary$RelPro < 25, "13-24%", 
                                                            ifelse(IQBP_FullSummary$RelPro > 25 & IQBP_FullSummary$RelPro < 50, "25-49%", 
                                                                   ifelse(IQBP_FullSummary$RelPro > 50 & IQBP_FullSummary$RelPro < 75, "50-74%", "75-100%")))))

IQBP_FullSummary$ForestCategory <- ifelse(IQBP_FullSummary$RelFor == 0, "0%",
                                   ifelse(IQBP_FullSummary$RelFor > 0 & IQBP_FullSummary$RelFor < 25, "1-24%", 
                                   ifelse(IQBP_FullSummary$RelFor >= 25 & IQBP_FullSummary$RelFor < 50, "25-49%", 
                                   ifelse(IQBP_FullSummary$RelFor >= 50 & IQBP_FullSummary$RelFor < 75, "50-74%","75-100%"))))


IQBP_FullSummary$"Inverse distance weighted % watershed forest" <- ifelse(IQBP_FullSummary$RelFor == 0, "0%",
                                          ifelse(IQBP_FullSummary$RelFor > 0 & IQBP_FullSummary$RelFor < 25, "1-24%", 
                                                 ifelse(IQBP_FullSummary$RelFor >= 25 & IQBP_FullSummary$RelFor < 50, "25-49%", 
                                                        ifelse(IQBP_FullSummary$RelFor >= 50 & IQBP_FullSummary$RelFor < 75, "50-74%","75-100%"))))

IQBP_FullSummary$AgriculturalCategory <- ifelse(IQBP_FullSummary$RelAgr == 0, "0%",
                                                ifelse(IQBP_FullSummary$RelAgr > 0 & IQBP_FullSummary$RelAgr < 25, "1-24%", 
                                                       ifelse(IQBP_FullSummary$RelAgr >= 25 & IQBP_FullSummary$RelAgr < 50, "25-49%", 
                                                              ifelse(IQBP_FullSummary$RelAgr >= 50 & IQBP_FullSummary$RelAgr < 75, "50-74%","75-100%"))))


IQBP_FullSummary$AntropogenicCategory <- ifelse(IQBP_FullSummary$PercentAnt == 0, "0%",
                                                ifelse(IQBP_FullSummary$PercentAnt > 0 & IQBP_FullSummary$PercentAnt < 11, "1-10%", 
                                                       ifelse(IQBP_FullSummary$PercentAnt >= 11 & IQBP_FullSummary$PercentAnt < 50, "11-49%", 
                                                              ifelse(IQBP_FullSummary$PercentAnt >= 50 & IQBP_FullSummary$PercentAnt < 75, "50-74%","75-100%"))))
table(IQBP_FullSummary$AntropogenicCategory)


IQBP_FullSummary$WaterQualityCategory <- ifelse(IQBP_FullSummary$IQBP_mean > 0 & IQBP_FullSummary$IQBP_mean < 20, "Very Bad",
                                                ifelse(IQBP_FullSummary$IQBP_mean >= 20 & IQBP_FullSummary$IQBP_mean < 39, "Bad",
                                                       ifelse(IQBP_FullSummary$IQBP_mean >= 40 & IQBP_FullSummary$IQBP_mean < 59, "Doubtful",
                                                              ifelse(IQBP_FullSummary$IQBP_mean >= 60 & IQBP_FullSummary$IQBP_mean < 79, "Satisfactory", "Good"))))

#Create a "natural land cover" variable that combines areas with Forest, wetlands, aquatic, and barren land
IQBP_FullSummary$RelPercentNatLand <- rowSums(cbind(IQBP_FullSummary$RelFor, IQBP_FullSummary$RelWet, IQBP_FullSummary$RelAqu, IQBP_FullSummary$RelBar))


#Watershed sizes
#Calculate the surface area of the watershed in square km
#Given the latitude correction, each pixel is about 0.0058 km 2 
IQBP_FullSummary$WatersheadAreakm2 <- IQBP_FullSummary$FlowAccumu*0.0058
range(IQBP_FullSummary$WatersheadAreakm2)
hist(IQBP_FullSummary$WatersheadAreakm2)


#Calculate the surface area of protection in each watershed
#Percent protection represents the percentage of the full watershed with protection
IQBP_FullSummary$ProtectionArea <- IQBP_FullSummary$WatersheadAreakm2 * IQBP_FullSummary$PercentPro/100

#Assess the number of assessed watersheds, and nested sub-watersheds
IQBPBasinsTable <- table(IQBP_FullSummary$IQBPBasins)
IQBPBasinsTable<- IQBP_FullSummary %>%
  group_by(IQBPBasins) %>%
  summarise (n = n())

IQBPMultipleBasinsTable <- subset(IQBPBasinsTable, n > 1)

range(IQBPMultipleBasinsTable$n)
mean(IQBPMultipleBasinsTable$n)
#Number of basins assessed:
sum(IQBPBasinsTable$n)
#Number of watersheds with nested sub-watersheds:
sum(IQBPBasinsTable$n > 1)
#Number of watersheds without nested sub-watersheds
sum(IQBPBasinsTable[IQBPBasinsTable$n < 1.1,]$n)

#Number of watersheds nested within nested sub-waterheds:
sum(IQBPBasinsTable[IQBPBasinsTable$n > 1,]$n)



#IQBP - Range of protection, water quality and land use in protected watersheds# 
IQBP_ProSubset <- IQBP_FullSummary %>% filter(PercentPro > 0)
range(IQBP_ProSubset$IQBP_mean)
mean(IQBP_ProSubset$IQBP_mean)

range(IQBP_ProSubset$PercentFor)
mean(IQBP_ProSubset$PercentFor)

range(IQBP_ProSubset$PercentAgg)
mean(IQBP_ProSubset$PercentAgg)

range(IQBP_ProSubset$PercentAnt)
mean(IQBP_ProSubset$PercentAnt)

#Values in unprotected watersheds
IQBP_UnProSubset <- IQBP_FullSummary %>% filter(PercentPro < 0.00001)
range(IQBP_UnProSubset$IQBP_mean)
mean(IQBP_UnProSubset$IQBP_mean)

#Extract copy of data
#write.csv(IQBP_FullSummary, file = "IQBP_FullSummary.csv")

#Create a subset of the data with the same amount of points with 0% protection as those with some level of protection
IQBP_EqualProAndNot<- IQBP_FullSummary %>% 
  group_by(ProtectionCategorical) %>% 
  do(sample_n(.,37))


#IQBP - Visuzaling raw data####
#Figure 3A - IQBP Relative Protection####
#Relative Protection
#color 
IQBP_RelPro <- ggplot(IQBP_FullSummary, aes(x=RelPro, y=IQBP_mean)) +
  geom_point(size = 3, aes(colour = ForestCategory, shape = ForestCategory)) +
  scale_color_manual(name="Inverse distance\\nweighted %\\nwatershed forest", 
                     breaks =c("0%", "1-24%", "25-49%","50-74%", "75-100%"),
                     labels = c("0%", "1-24%", "25-49%","50-74%", "75-100%"),
                      values=c("#a6611a", "#dfc27d",  "grey","#80dcd1","#018571")) +
  scale_shape_discrete(name="Inverse distance\\nweighted %\\nwatershed forest", 
                        breaks =c("0%", "1-24%", "25-49%","50-74%", "75-100%"),
                        labels = c("0%", "1-24%", "25-49%","50-74%", "75-100%"))+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=16,face="bold"), 
        panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.title=element_text(size=13, face = "bold"), legend.text=element_text(size=12), 
        #legend.position = c(0.6, 0.3), 
        legend.background = element_blank(),legend.box.background = element_rect(colour = "white"))+
  geom_hline(yintercept = c(19, 39, 59, 79), colour = c("black", "black", "black", "black"), linetype = "dashed") +
  annotate("text", x = 50, y = 8, angle = 360, size = 4.5, label = "Very bad        ")+
  annotate("text", x = 50, y = 29, angle = 360,size = 4.5, label = "Bad")+
  annotate("text", x = 50, y = 49, angle = 360,size = 4.5, label = "Questionable               ")+
  annotate("text", x = 50, y = 69, angle = 360, size = 4.5,label = "Satisfactory              ")+
  annotate("text", x = 50, y = 90, angle = 360,size = 4.5, label = "Good    ")+
  annotate("text", x = 50, y = 99, angle = 360,size = 4.5, label = "")+
  xlab("Inverse distance weighted %\\nwatershed protection") +
  ylab("Water quality\\n(IQBP indicator)")
  
IQBP_RelPro

ggsave(IQBP_RelPro, file="Figure3A_IQBPRelPro.png", scale=1, dpi = 600)



#IQBP - Exploring Data Distribution####
hist(IQBP_FullSummary$IQBP_mean)
IQBP_FullSummary$logIQBPmean<- log(IQBP_FullSummary$IQBP_mean)
hist(IQBP_FullSummary$logIQBPmean)

IQBP_FullSummary$sqrtIQBPmean<- sqrt(IQBP_FullSummary$IQBP_mean)
hist(IQBP_FullSummary$sqrtIQBPmean)

IQBP_FullSummary$logitIQBPmean<- logit(IQBP_FullSummary$IQBP_mean, percent=TRUE, adjust=0)
hist(IQBP_FullSummary$logitIQBPmean)


#Z Correct the data so that all variables are on the same scale 
IQBP_FullSummary$ZlogitIQBP<-scale(IQBP_FullSummary$logitIQBPmean, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZPercentPro<-scale(IQBP_FullSummary$PercentPro, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZPercentAgg<-scale(IQBP_FullSummary$PercentAgg, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZPercentWet<-scale(IQBP_FullSummary$PercentWet, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZPercentBar<-scale(IQBP_FullSummary$PercentBar, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZPercentAqu<-scale(IQBP_FullSummary$PercentAqu, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZPercentAnt<-scale(IQBP_FullSummary$PercentAnt, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZPercentCut<-scale(IQBP_FullSummary$PercentCut, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZPercentFor<-scale(IQBP_FullSummary$PercentFor, center = TRUE, scale =TRUE)

IQBP_FullSummary$ZRelPro<-scale(IQBP_FullSummary$RelPro, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZFlowAccumu<-scale(IQBP_FullSummary$FlowAccumu, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZRelAgr<-scale(IQBP_FullSummary$RelAgr, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZRelWet<-scale(IQBP_FullSummary$RelWet, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZRelBar<-scale(IQBP_FullSummary$RelBar, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZRelAqu<-scale(IQBP_FullSummary$RelAqu, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZRelAnt<-scale(IQBP_FullSummary$RelAnt, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZRelCut<-scale(IQBP_FullSummary$RelCut, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZRelFor<-scale(IQBP_FullSummary$RelFor, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZLatitude<-scale(IQBP_FullSummary$Latitude, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZLongitude<-scale(IQBP_FullSummary$Longitude, center = TRUE, scale =TRUE)
IQBP_FullSummary$ZRelPercentNatLand <-scale(IQBP_FullSummary$RelPercentNatLand, center = TRUE, scale =TRUE)


#IQBP - Mixed Model Analysis####
#with weighted Protection##
#Look at correlations between variables
IQBP_Corr <- select(IQBP_FullSummary, "ZlogitIQBP", "ZRelPro", "ZFlowAccumu", "ZRelWet", "ZRelBar", "ZRelAqu", "ZRelCut", "ZRelFor", "ZRelAnt", "ZRelAgr", "ZLatitude", "ZLongitude", "ZRelPercentNatLand")
plot(IQBP_Corr)
cor(IQBP_Corr, use="complete.obs", method = "kendall")
#Agriculture is 61% correlated to IQBP, and forests are 58% correlated to IQBP. Everything else is lower. As long as all is below
#80% this is considered to be OK

#Create a model that contains all variables 
M0_IQBP <- lmer(ZlogitIQBP ~ ZRelPro + ZRelAgr +ZRelFor+  ZRelWet + ZRelCut + ZRelAnt + ZRelAqu + ZRelBar + ZFlowAccumu + ZLatitude + ZLongitude + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)

#Use the MuMIn package for model selection and averaging
options(na.action = "na.fail") 
resultsM0_IQBP <- dredge(M0_IQBP)
options(na.action = "na.omit") 
top.modelsIQBP <- get.models(resultsM0_IQBP, subset = delta < 5)
top.modelsIQBP[1]
model.avg(top.modelsIQBP)
aic_impIQBP <- importance(top.modelsIQBP)
aic_impIQBP <- as.data.frame((aic_impIQBP))
aic_impIQBP$variables <- rownames(aic_impIQBP)
aic_impIQBP$VariableNumbers <- c("K","J","I","H","G","F","E","D","C","B","A")

#Figure 3B - IQBP Weighted Rel Importance####
# Plotting models Delta AIC < 5
IQBPWeightedRelImportance <-   ggplot(data=aic_impIQBP, aes(y=aic_impIQBP[,1],x=VariableNumbers)) +
  geom_bar(position="dodge", stat="identity") + 
  coord_flip() + xlab("Predictor variable") + ylab("Relative importance (proportion)")  +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15,face="bold"),
        panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #Add information about model estimates 
  annotate("text", x = 12.5, y = 1.12, angle = 0, label = "")+
  annotate("text", x = 12, y = 1.08, angle = 0, size = 4.5, label = "Average\\nSlope")+
  annotate("text", x = 11.15, y = 1.08, angle = 0,size = 4.5, label = "      0.64 + * ")+ #Forest
  annotate("text", x = 10.15, y = 1.08, angle = 0,size = 4.5, label = "      0.30 + * ")+ #Aqu
  annotate("text", x = 9.15, y = 1.08, angle = 0,size = 4.5, label = "      0.22 + * ")+#Cut
  annotate("text", x = 8.15, y = 1.08, angle = 0, size = 4.5,label = "    - 0.12 + * ") + #Bar
  annotate("text", x = 7.15, y = 1.08, angle = 0, size = 4.5,label = "      0.14 + * ")+ #Ant
  annotate("text", x = 6.15, y = 1.08, angle = 0, size = 4.5,label = "  0.09  ")+ # Longitude
  annotate("text", x = 5.15, y = 1.08, angle = 0, size = 4.5,label = " - 0.11  ")+  #Agri
  annotate("text", x = 4.15, y = 1.08, angle = 0, size = 4.5,label = "- 0.04  ")+ #Watershed Size
  annotate("text", x = 3.15, y = 1.08, angle = 0, size = 4.5,label = "  0.02  ")+  #Protection
  annotate("text", x = 2.15, y = 1.08, angle = 0, size = 4.5,label = "  0.00  ")+ #wetland
  annotate("text", x = 1.15, y = 1.08, angle = 0, size = 4.5, label = "  0.02  ")+ #Latitude
  
  scale_x_discrete(labels= c( "Latitude",
                            "IDW% wetlands",
                           "IDW% protection",
                             "Watershed size",
                             "IDW% agricultural",
                           "Longitude",
                             "IDW% anthropogenic",
                           "IDW% rock & lichens",
                             "IDW% planted, cut,\\n& burned forests", 
                             "IDW% aquatic", 
                             "IDW% forest")) 
IQBPWeightedRelImportance
ggsave(IQBPWeightedRelImportance, file="Figure3B_IQBPWeightedRelImportance.tiff", scale=1, dpi = 600)

grid.arrange(IQBP_RelPro, IQBPWeightedRelImportance, labels = c("A", "B"), ncol = 2)


#Validate best fit model 
M1_IQBP <- lmer(ZlogitIQBP ~ZRelAnt +ZRelAqu + ZRelBar + ZRelCut + ZRelFor + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)

vif(M1_IQBP)
#VIFs below 5 or 10 are considered acceptable 

par(mfrow=c(3,1)) 
#A) Look at homogeneity
ResidM1_IQBP<-resid(M1_IQBP)
FM1_IQBP<-fitted(M1_IQBP)
plot(x=FM1_IQBP, y= ResidM1_IQBP, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model

#Random factor
boxplot(ResidM1_IQBP~IQBPBasins, ylab= "Normalized residuals", data = IQBP_FullSummary, xlab = "Nested Basins")
abline(0,0, lty=2)

#In the model
plot(x=IQBP_FullSummary$ZRelAqu, y= ResidM1_IQBP, xlab = "Percent Upland Aquatic (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelAnt, y= ResidM1_IQBP, xlab = "Percent Upland Anthropogenic (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelBar, y= ResidM1_IQBP, xlab = "Percent Upland Barren Land and Rocks (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelCut, y= ResidM1_IQBP, xlab = "Percent Upland Deforestation and Burns (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelFor, y= ResidM1_IQBP, xlab = "Percent Upland Forest (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)

#Not in the model
plot(x=IQBP_FullSummary$ZRelPro, y= ResidM1_IQBP, xlab = "Percent Upland Protected (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZLatitude, y= ResidM1_IQBP, xlab = "Latitude", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelWet, y= ResidM1_IQBP, xlab = "Percent Upland Wetland (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZLongitude, y= ResidM1_IQBP, xlab = "Longitude", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelAgr, y= ResidM1_IQBP, xlab = "Percent Upland Agriculture (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZFlowAccumu, y= ResidM1_IQBP, xlab = "Watershed Size (Flow Accumulation)", ylab = "Normalized residuals")
abline(0,0, lty=2)

#C) Look at normality of residuals
hist(ResidM1_IQBP)

#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(M1_IQBP)

#Assess significance of fixed factors 
summary(M1_IQBP)
#p value for each variable
M1_IQBP <- lmer(ZlogitIQBP ~ ZRelAnt +ZRelAqu + ZRelBar + ZRelCut + ZRelFor + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)

#Aquat
MAqu_IQBP <- lmer(ZlogitIQBP ~ ZRelAnt+ ZRelBar + ZRelCut + ZRelFor +  (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)
anova(M1_IQBP, MAqu_IQBP)
#Ant
MAnt_IQBP <- lmer(ZlogitIQBP ~ ZRelAqu + ZRelBar + ZRelCut + ZRelFor +  (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)
anova(M1_IQBP, MAnt_IQBP)
#Barren land 
MBar_IQBP<- lmer(ZlogitIQBP ~ ZRelAnt +ZRelAqu + ZRelCut + ZRelFor + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)
anova(M1_IQBP, MBar_IQBP)
#Cut
MCut_IQBP <- lmer(ZlogitIQBP ~ ZRelAnt +ZRelAqu + ZRelBar + ZRelFor + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)
anova(M1_IQBP, MCut_IQBP)
#Forest
MFor_IQBP<- lmer(ZlogitIQBP ~ ZRelAnt +ZRelAqu + ZRelBar + ZRelCut  + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)
anova(M1_IQBP, MFor_IQBP)

#IQBP - T test comparing watersheds with and without protection####
#Visualize Water Quality at sites without and without protection  
#Boxplot with 2 levels showing differences
#Figure 5A####
IQBP_ProCatViolin <- ggplot(IQBP_FullSummary, aes(x=ProtectionCategorical, y =IQBP_mean)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point", shape=19, size=3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=16,face="bold"),
       panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
       panel.background = element_blank(), axis.line = element_line(colour = "black"), 
       legend.title=element_text(size=15), legend.text=element_text(size=14))+
  geom_hline(yintercept = c(19, 39, 59, 79), colour = c("black", "black", "black", "black"), linetype = "dashed") +
  annotate("text", x = 2.7, y = 8, angle = 360, size= 4.5, label = "Very bad")+
  annotate("text", x = 2.7, y = 29, angle = 360, size= 4.5,label = "Bad")+
  annotate("text", x = 2.7, y = 50, angle = 360, size= 4.5,label = "Questinonable    ")+
  annotate("text", x = 2.7, y = 69, angle = 360,size= 4.5, label = "Satisfactory  ")+
  annotate("text", x = 2.7, y = 90, angle = 360,size= 4.5, label = "Good")+
  annotate("text", x = 2.7, y = 8, angle = 360,size= 4.5, label = " ") +
  annotate("text", x = 3, y = 101, angle = 360, label = " ") +
  xlab("Protection status") +
  ylab("Water quality\\n(IQBP indicator)")
IQBP_ProCatViolin

ggsave("Figure5A_IQBPProCatViolin.tiff")


#A) compare water quality between protected and unprotected sites 
#Welch's T-test (allows for unequal sample sizes)
t.test(IQBP_FullSummary$logitIQBPmean~ IQBP_FullSummary$ProtectionCategorical, var.equal=FALSE)

#Verify model assumptions
#is the data from both groups normally distributed 
hist(with(IQBP_FullSummary,logitIQBPmean[ProtectionCategorical == "No Protection"]))
with(IQBP_FullSummary, shapiro.test(logitIQBPmean[ProtectionCategorical == "No Protection"]));
hist(with(IQBP_FullSummary, logitIQBPmean[ProtectionCategorical == "Some level of protection"]))
with(IQBP_FullSummary, shapiro.test(logitIQBPmean[ProtectionCategorical == "Some level of protection"]))
#Not "some level of protection", so move on to a non-parametric test 


#Mann_Whitney U Test
#unequal sample sizes also OK: https://stats.stackexchange.com/questions/40342/mann-whitney-u-test-with-unequal-sample-sizes
IQBP_FullSummary$ProtectionCategorical <- as.factor(IQBP_FullSummary$ProtectionCategorical)
wilcox.test(IQBP_mean ~ ProtectionCategorical, data = IQBP_FullSummary)


#Compute effect size
#Compute Vargha and Delaney's A to assess effect size 
VD.A(d = IQBP_FullSummary$IQBP_mean,
     f = IQBP_FullSummary$ProtectionCategorical)
#Effect size is "small" 

#B) determine if  water quality related to % natural land cover within the watershed?
#Run mixed model 
IQBPlmer <- lmer(IQBP_mean ~ RelPercentNatLand+ (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)


#Validate  model 

#A) Look at homogeneity
ResidIQBPlmer<-resid(IQBPlmer)
FitIQBPlmer<-fitted(IQBPlmer)
plot(x=FitIQBPlmer, y= ResidIQBPlmer, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model 
#Random factor
boxplot(ResidIQBPlmer~IQBPBasins, ylab= "Normalized residuals", data = IQBP_FullSummary, xlab = "Nested Basins")
abline(0,0, lty=2)

#Fixed factors not in the model
plot(x=IQBP_FullSummary$ZRelAnt, y= ResidIQBPlmer, xlab = "Percent Upland Agricultural (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZFlowAccumu, y= ResidIQBPlmer, xlab = "Percent Upland Agricultural (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelAgr, y= ResidIQBPlmer, xlab = "Percent Upland Agricultural (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)

#C) Look at normality of residuals
hist(ResidIQBPlmer)

#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(IQBPlmer)
#Assess significance of fixed factors 
summary(IQBPlmer)


#p value for Natural Land Cover
IQBPlmerRandomEffOnly <- lmer(IQBP_mean ~ 1 + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)
anova(IQBPlmer, IQBPlmerRandomEffOnly)


#Visualize findings
#visualize regression
#Figure 5C####
IQBP_NatLand<-ggplot(IQBP_FullSummary, aes(x=RelPercentNatLand, y=IQBP_mean)) +
  geom_point() +
 #geom_abline(intercept = -0.004193, slope = 0.854515) +
  #geom_segment(aes(x =0, xend = 100, y = -0.004193 + 0.854515, yend =-0.004193 + 0.854515*100)) +
  #geom_smooth(method ='lm', colour = "black") +
  geom_segment(aes(x =0, xend = 100, y = 11.69776 +  0.83670, yend = 11.69776 +  0.83670*100)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=13,face="bold"), 
        panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title=element_text(size=12), legend.text=element_text(size=11))+
  geom_hline(yintercept = c(19, 39, 59, 79), colour = c("black", "black", "black", "black"), linetype = "dashed") +
  annotate("text", x = 105, y = 8, angle = 360, label = "Very bad")+
  annotate("text", x = 105, y = 29, angle = 360, label = "Bad")+
  annotate("text", x = 105, y = 49, angle = 360, label = "Questinonable    ")+
  annotate("text", x = 105, y = 69, angle = 360, label = "Satisfactory")+
  annotate("text", x = 105, y = 90, angle = 360, label = "Good")+
  annotate("text", x = 108, y = 102, angle = 360, label = "")+
  xlab("Inverse distance weighted % upland natural land cover") +
  ylab("Water quality \\n(IQBP indicator)")
IQBP_NatLand
#Figure5C_IQBPNatLandCover

#C) compare natural land cover between protected and unprotected watersheds 

#compare statistically 
#T-test
t.test(IQBP_FullSummary$RelPercentNatLand~ IQBP_FullSummary$ProtectionCategorical, var.equal=FALSE)

#Verify model assumptions
#is the data from both groups normally distributed 
hist(with(IQBP_FullSummary,RelPercentNatLand[ProtectionCategorical == "No IUCN Category\\nII protection"]))
with(IQBP_FullSummary, shapiro.test(RelPercentNatLand[ProtectionCategorical == "No IUCN Category\\nII protection"]));
#No, so move on to a non-parametric test 


#Mann_Whitney U Test
wilcox.test(RelPercentNatLand ~ ProtectionCategorical, data = IQBP_FullSummary)
#There is a significant difference 

#Compute effect size 
VD.A(d = IQBP_FullSummary$RelPercentNatLand,
     f = IQBP_FullSummary$ProtectionCategorical)


#Visualize 
#Figure 5B####
IQBP_NatLandViolin<- ggplot(IQBP_FullSummary, aes(x=ProtectionCategorical, y = RelPercentNatLand)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point", shape=19, size=3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=16,face="bold"), 
        panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title=element_text(size=15), legend.text=element_text(size=14))+
  annotate("text", x = 2.5, y = 101, angle = 90, label = " ") +
  xlab("Protection status") +
  ylab("Inverse distance weighted %\\n upland natural land cover")
IQBP_NatLandViolin
#Figure5B_IQBPProNatViolin

#IQBP - Chi Squared tests####
#http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r

#on sites with protection
#Reorder factors 
IQBP_ProSubset$WaterQualityCategory <- as.factor(IQBP_ProSubset$WaterQualityCategory)
levels(IQBP_ProSubset$WaterQualityCategory)
IQBP_ProSubset$WaterQualityCategory <- factor(IQBP_ProSubset$WaterQualityCategory, 
                                              levels(IQBP_ProSubset$WaterQualityCategory)[c(5,1,2,4,3)])
#Chi-Square Test
IQBPProchisq<- chisq.test(IQBP_ProSubset$ProtectionCategory,IQBP_ProSubset$WaterQualityCategory, correct= FALSE)
IQBPProchisq

IQBP_ProSubsetTable <- IQBP_ProSubset %>%
  group_by(ProtectionCategory, WaterQualityCategory) %>%
  summarize(n = n()) %>%
  mutate(freq = n/sum(n))

#Contingency Table 
IQBP_ProSubsetTable<-as.data.frame(IQBP_ProSubsetTable)
IQBP_ProSubsetTable4Spread<- IQBP_ProSubsetTable[,-3]
IQBP_ProSubsetContingencyTable<-spread(IQBP_ProSubsetTable4Spread,ProtectionCategory, freq)
IQBP_ProSubsetContingencyTable[is.na(IQBP_ProSubsetContingencyTable)] <- 0
IQBP_ProSubsetContingencyTable

write.csv(IQBP_ProSubsetContingencyTable, "IQBP_ProSubsetContingencyTable.csv")

#Using all data 
IQBPchisq<- chisq.test(IQBP_FullSummary$ProtectionCategory,IQBP_FullSummary$WaterQualityCategory, correct= FALSE)
IQBPchisq

IQBP_Table <- IQBP_FullSummary %>%
  group_by(ProtectionCategory, WaterQualityCategory) %>%
  summarize(n = n()) %>%
  mutate(freq = n/sum(n))

IQBP_MeanWQ <- IQBP_FullSummary %>%
  group_by(ProtectionCategory) %>%
  summarize(meanWQ = mean(IQBP_mean), minWQ = min(IQBP_mean), maxWQ = max(IQBP_mean), seWQ = se(IQBP_mean))

#Contingency Table 
IQBP_Table<-as.data.frame(IQBP_Table)
IQBP_Table4Spread<- IQBP_Table[,-3]
IQBP_ContingencyTable<-spread(IQBP_Table4Spread,ProtectionCategory, freq)
IQBP_ContingencyTable[is.na(IQBP_ContingencyTable)] <- 0
IQBP_ContingencyTable

write.csv(IQBP_ContingencyTable, "IQBP_ProSubsetContingencyTable.csv")

#Comparision between protected and unprotected sites
IQBPchisqProNot<- chisq.test(IQBP_FullSummary$ProtectionCategorical,IQBP_FullSummary$WaterQualityCategory, correct= FALSE)
IQBPchisqProNot 
#Small Contingency Table 
IQBP_SmallTable <- IQBP_FullSummary %>%
  group_by(ProtectionCategorical, WaterQualityCategory) %>%
  summarize(n = n()) %>%
  mutate(freq = n/sum(n))

#Contingency Table 
IQBP_SmallTable<-as.data.frame(IQBP_SmallTable)
IQBP_SmallTable4Spread<- IQBP_SmallTable[,-3]
IQBP_SmallContingencyTable<-spread(IQBP_SmallTable4Spread,ProtectionCategorical, freq)
IQBP_SmallContingencyTable[is.na(IQBP_SmallContingencyTable)] <- 0
IQBP_SmallContingencyTable

write.csv(IQBP_SmallContingencyTable, "IQBP_SmallContingencyTable.csv")



#ISBs DATA ANLAYSIS####
#ISBs -Load and clean up data####
#SB_All <- read.csv("SB_Tout.csv")
SB_All<-read.csv("Hanna_Chp2SBGISData2.csv", fileEncoding= 'utf8' )
IQBR<-read.csv("IQBR.csv", fileEncoding= 'utf8' )


SB_All<-left_join(SB_All, IQBR, by = "Identifian")
colnames(SB_All)[colnames(SB_All) == 'IQBR.y'] <- 'IQBR'

#Add a column of categories of protection
SB_All$ProtectionCategory <- ifelse(SB_All$PercentPro == 0, "0%",
                                    ifelse(SB_All$PercentPro > 0 & SB_All$PercentPro < 12, "1-12%", 
                                           ifelse(SB_All$PercentPro > 12 & SB_All$PercentPro < 25, "13-24%",
                                           ifelse(SB_All$PercentPro > 25 & SB_All$PercentPro < 50, "25-49%", 
                                                  ifelse(SB_All$PercentPro > 50 & SB_All$PercentPro < 75, "50-74%","75-100%")))))
table(SB_All$ProtectionCategory)

#SB_All$ProtectionCategory <- ifelse(SB_All$RelPro == 0, "0%",
#                                    ifelse(SB_All$RelPro > 0 & SB_All$RelPro < 25, "1-24%", 
#                                           ifelse(SB_All$RelPro > 25 & SB_All$RelPro < 50, "25-49%", 
#                                                  ifelse(SB_All$RelPro > 50 & SB_All$RelPro < 75, "50-74%","75-100%"))))


SB_All$ProtectionCategorical <- ifelse(SB_All$PercentPro == 0, "No IUCN Category\\nII protection", "Some level of IUCN\\nCategory II protection")
SB_ProCategorical<- SB_All %>%
  group_by(ProtectionCategorical) %>%
  summarise (n = n())

table(SB_All$ProtectionCategorical)


SB_All$AgriCategory <- ifelse(SB_All$RelAgr == 0, "0%",
                              ifelse(SB_All$RelAgr > 0 & SB_All$RelAgr < 25, "1-24%", 
                                     ifelse(SB_All$RelAgr > 25 & SB_All$RelAgr < 50, "25-49%", 
                                            ifelse(SB_All$RelAgr > 50 & SB_All$RelAgr < 75, "50-74%","75-100%"))))

SB_All$ForestCategory <- ifelse(SB_All$RelFor == 0, "0%",
                                ifelse(SB_All$RelFor > 0 & SB_All$RelFor < 25, "1-24%", 
                                       ifelse(SB_All$RelFor > 25 & SB_All$RelFor < 50, "25-49%", 
                                              ifelse(SB_All$RelFor > 50 & SB_All$RelFor < 75, "50-74%","75-100%"))))

SB_All$AnthropogenicCategory <- ifelse(SB_All$RelAnt == 0, "0%",
                                       ifelse(SB_All$RelAnt > 0 & SB_All$RelAnt < 25, "1-24%", 
                                              ifelse(SB_All$RelAnt > 25 & SB_All$RelAnt < 50, "25-49%", 
                                                     ifelse(SB_All$RelAnt > 50 & SB_All$RelAnt < 75, "50-74%","75-100%"))))


SB_All$WaterQualityCategory <- ifelse(SB_All$ISVB > 0 & SB_All$ISVB < 45, "Bad",
                                                ifelse(SB_All$ISVB >= 45 & SB_All$ISVB < 75, "Precarious", "Good"))
                
#Create a "natural land cover" variable that combines areas with Forest, wetlands, aquatic, and barren land
SB_All$RelPercentNatLand <- rowSums(cbind(SB_All$RelFor, SB_All$RelWet, SB_All$RelAqu, SB_All$RelBar))
                                      
#Calculate the surface area of the watershed in square km
#Each pixel is about 0.0058 km 2 
SB_All$WatersheadAreakm2 <- SB_All$FlowAccumu*0.0058
range(SB_All$WatersheadAreakm2)
hist(SB_All$WatersheadAreakm2)

#Calculate the surface area of protection in each watershed
#Percent protection represents the percentage of the full watershed with protection

#Assessing the number of assessed watersheds, and nested sub-watersheds
ISBsBasinsTable <- table(SB_All$SBWBasins)
ISBsBasinsTable<- SB_All %>%
  group_by(SBWBasins) %>%
  summarise (n = n())

#Number of basins assessed:
sum(ISBsBasinsTable$n)
#Number of watersheds with nested sub-watersheds:
sum(ISBsBasinsTable$n > 1)
#Number of watersheds without nested sub-watersheds
sum(ISBsBasinsTable[ISBsBasinsTable$n < 1.1,]$n)
#Number of watersheds nested within nested sub-waterheds:
sum(ISBsBasinsTable[ISBsBasinsTable$n > 1,]$n)

#ISBs - Calculate range of protection, water quality, and land use in protected watersheds
SB_AllSubset <- SB_All %>% select(ISVB, Identifian, PercentPro, FlowAccumu, WatersheadAreakm2, PercentAgg, PercentFor, PercentAnt)
SB_ProSubset <- SB_All%>% filter(PercentPro > 0)
range(SB_ProSubset$ISVB)
range(SB_ProSubset$WatersheadAreakm2)

SB_UnProSubset <- SB_All %>% filter(PercentPro < 0.00001)


#ISBs -Visuzaling raw data####
#Distance weighted 
#Figure 4A####
ISVB_RelPro<-ggplot(SB_All, aes(x=RelPro, y=ISVB, color = AgriCategory, shape = AgriCategory)) +  
  geom_point(size = 3) +
  scale_color_manual(name="Inverse distance\\nweighted %\\nwatershed\\nagriculture", 
                   breaks =c("0%", "1-24%", "25-49%","50-74%"),
                     labels = c("0%", "1-24%", "25-49%","50-74%"),
                     values=c("#fed98e",  "#fe9929","#cc4c02", "#993404")) +
  scale_shape_discrete(name="Inverse distance\\nweighted %\\nwatershed\\nagriculture", 
                       breaks =c("0%", "1-24%", "25-49%","50-74%", "75-100%"),
                     labels = c("0%", "1-24%", "25-49%","50-74%", "75-100%"))+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=16,face="bold"), 
        panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.title=element_text(size=15, face="bold"), legend.text=element_text(size=14),
        #legend.position = c(0.7, 0.2), 
        legend.background = element_blank(),legend.box.background = element_rect(colour = "white"))+
  geom_hline(yintercept = c(45, 75), colour = c("black", "black"), linetype = "dashed") +
  
   annotate("text", x = 110, y = 35, angle = 360, size = 5, label = "Bad  ")+
  annotate("text", x = 110, y = 60, angle = 360, size = 5,label = "Precarious              ")+
  annotate("text", x = 110, y = 90, angle = 360, size = 5,label = "Good ")+
  xlab("Inverse distance weighted\\n % watershed protection") +
 # xlab("Distance Weighted Percent Watershed Protection\\n (Indicator of the Amount of Protection\\nin the Sampled Points' Watersheds)") +
  #ylab("Indicator of Benthic Invertebrate Health \\n(ISBs)")
  ylab("Water quality\\n(ISBs indicator)")
ISVB_RelPro

ggsave(ISVB_RelProUpdated, file="Figure4A_ISBsRelPro.pdf", scale=1, dpi = 600)



#ISBs - Data Distribution Exploration####
hist(SB_All$ISVB)
hist(SB_All$RelPro, main= "", xlab= "Percent protection of watershed (weighted)")
#hist(SB_All$IQBR_Point, main= "", xlab= "River bank quality index")
hist(SB_All$IQBR, main= "", xlab= "River bank quality index (weighted)")
#hist(SB_All$FlowAccumu, main= "", xlab= "Relative watershed size")
hist(SB_All$FlowAccumu, main= "", xlab= "Relative watershed size (weighted)")
hist(SB_All$RelAgr, main= "", xlab= "Percent agricultural of watershed (weighted)")
hist(SB_All$RelWet, main= "", xlab= "Percent wetland of watershed (weighted)")
hist(SB_All$RelBar, main= "", xlab= "Percent rock and shrubs of watershed (weighted)")
hist(SB_All$RelAqu, main= "", xlab= "Percent aquatic of watershed (weighted)")
hist(SB_All$RelAnt, main= "", xlab= "Percent anthropogenic of watershed (weighted)")
hist(SB_All$RelCut, main= "", xlab= "Percent cut and regenerating forests of watershed (weighted)")
hist(SB_All$RelFor, main= "", xlab= "Percent forest of watershed (weighted)")


#Try to adjust data distribution of dependant variable to normalize it
SB_All$logitISVB<- logit(SB_All$ISVB, percent=TRUE, adjust=0)
hist(SB_All$logitISVB)

#Scale and centre data###
SB_All$ZlogitISVB<-scale(SB_All$logitISVB, center = TRUE, scale =TRUE)
SB_All$ZISVB<-scale(SB_All$logitISVB, center = TRUE, scale =TRUE)
SB_All$ZPercentPro<-scale(SB_All$PercentPro, center = TRUE, scale =TRUE)
SB_All$ZFlowAccumu<-scale(SB_All$FlowAccumu, center = TRUE, scale =TRUE)
SB_All$ZPercentAgg<-scale(SB_All$PercentAgg, center = TRUE, scale =TRUE)
SB_All$ZPercentWet<-scale(SB_All$PercentWet, center = TRUE, scale =TRUE)
SB_All$ZPercentBar<-scale(SB_All$PercentBar, center = TRUE, scale =TRUE)
SB_All$ZPercentAqu<-scale(SB_All$PercentAqu, center = TRUE, scale =TRUE)
SB_All$ZPercentAnthro<-scale(SB_All$PercentAnt, center = TRUE, scale =TRUE)
SB_All$ZPercentCut<-scale(SB_All$PercentCut, center = TRUE, scale =TRUE)
SB_All$ZPercentFor<-scale(SB_All$PercentFor, center = TRUE, scale =TRUE)

#Weighted data
SB_All$ZRelPro<-scale(SB_All$RelPro, center = TRUE, scale =TRUE)
SB_All$ZRelAgr<-scale(SB_All$RelAgr, center = TRUE, scale =TRUE)
SB_All$ZRelWet<-scale(SB_All$RelWet, center = TRUE, scale =TRUE)
SB_All$ZRelBar<-scale(SB_All$RelBar, center = TRUE, scale =TRUE)
SB_All$ZRelAqu<-scale(SB_All$RelAqu, center = TRUE, scale =TRUE)
SB_All$ZRelAnthro<-scale(SB_All$RelAnt, center = TRUE, scale =TRUE)
SB_All$ZRelCut<-scale(SB_All$RelCut, center = TRUE, scale =TRUE)
SB_All$ZRelFor<-scale(SB_All$RelFor, center = TRUE, scale =TRUE)

#Other variables
SB_All$ZCond<-scale(SB_All$StreamConductivity, center = TRUE, scale =TRUE)
SB_All$ZDO<-scale(SB_All$StreamDO, center = TRUE, scale =TRUE)
SB_All$ZpH<-scale(SB_All$StreampH, center = TRUE, scale =TRUE)
SB_All$ZTemp<-scale(SB_All$StreamTemperature, center = TRUE, scale =TRUE)
SB_All$ZSpeed<-scale(SB_All$AverageSpeed, center = TRUE, scale =TRUE)
SB_All$ZAltitude<-scale(SB_All$Altitude, center = TRUE, scale =TRUE)
SB_All$ZLatitude<-scale(SB_All$Latitude, center = TRUE, scale =TRUE)
SB_All$ZLongitude<-scale(SB_All$Longitude, center = TRUE, scale =TRUE)
SB_All$ZFlowAccumu <-scale(SB_All$FlowAccumu, center = TRUE, scale =TRUE)
SB_All$ZIQBR<-scale(SB_All$IQBR, center = TRUE, scale =TRUE)
SB_All$ZRelPercentNatLand<-scale(SB_All$RelPercentNatLand, center = TRUE, scale =TRUE)



#ISBs - Mixed Model Analysis####
#Create a model that contains all variables 
#SB_Model<-subset(SB_All[,c("Identifian","SBWBasins","ZlogitISVB", "ZIQBR", "ZRelPro", "ZRelAgr", "ZRelFor", "ZRelWet", "ZRelCut", "ZRelAnthro", "ZRelAqu", "ZFlowAccumu", "ZTemp", "ZAltitude", "ZLatitude", "ZLongitude" )])
M1W <- lmer(ZlogitISVB ~ ZRelPro +ZRelAgr +ZRelFor + ZRelWet + ZRelCut + ZRelAnthro+ ZRelAqu + ZRelBar+ ZFlowAccumu +ZLatitude +ZLongitude + (1|SBWBasins), data=SB_All,REML= FALSE)

SB_All4Corr <- select(SB_All, "ZlogitISVB", "ZRelPro", "ZFlowAccumu", "ZRelWet", "ZRelBar", "ZRelAqu", "ZRelCut", "ZRelAnthro",
                      "ZRelAgr", "ZRelFor", "ZLatitude", "ZLongitude", "ZRelPercentNatLand")
plot(SB_All4Corr)
cor(SB_All4Corr, use="complete.obs", method = "kendall")
#No variables about 43% so we're OK to include them all

#Use the MuMIn package for model selection and averaging
#OK to model average with mixed model - https://stats.stackexchange.com/questions/351917/model-averaging-in-mixed-models
options(na.action = "na.fail") 
resultsM1W <- dredge(M1W)
options(na.action = "na.omit") 

top.models <- get.models(resultsM1W, subset = delta < 5)
top.models[1]
model.avg(top.models)
aic_imp <- importance(top.models)
aic_imp <- as.data.frame((aic_imp))
aic_imp$variables <- rownames(aic_imp)
aic_imp$VariableNumbers <- c("K","I","J","H","G","F","E","D","C","B","A")

#Figure 4B####
# Plotting models Delta AIC < 5
ISBsWeightedRelImportance <-   ggplot(data=aic_imp, aes(y=aic_imp[,1],x=VariableNumbers)) +
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() + xlab("Predictor Variable") + ylab("Relative importance")  +
  theme(axis.text=element_text(size=12.5), axis.title=element_text(size=15,face="bold"),
        panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #Add information about model estimates 
  annotate("text", x = 12.5, y = 1.15, angle = 0, label = "")+
  annotate("text", x = 12, y = 1.08, angle = 0, label = "Average\\nSlope")+
  annotate("text", x = 11.14, y = 1.08, angle = 0, label = "     - 0.37 + *")+ #Agricultural 
  annotate("text", x = 10.15, y = 1.08, angle = 0, label = "      0.24 + *")+ #wetland
  annotate("text", x = 9.15, y = 1.08, angle = 0, label = "     - 0.25 + *")+ #anthro
  annotate("text", x = 8.15, y = 1.08, angle = 0, label = "   0.15  ")+ #Longitude
  annotate("text", x = 7.15, y = 1.08, angle = 0, label = "     0.12 +") + #Forest
  annotate("text", x = 6.15, y = 1.08, angle = 0, label = "- 0.11") + # Latitude
  annotate("text", x = 5.15, y = 1.08, angle = 0, label = "- 0.08")+ #Watershed size
  annotate("text", x = 4.15, y = 1.08, angle = 0, label = "- 0.05")+ # Aquatic
  annotate("text", x = 3.15, y = 1.08, angle = 0, label = "  0.03") + #Cut
  annotate("text", x = 2.15, y = 1.08, angle = 0, label = "  0.03")+ #Pro
  annotate("text", x = 1.15, y = 1.08, angle = 0, label = "  0.00")+ #Bar
  scale_x_discrete(labels= c("IDW% rock\\n& lichens",
                             "IDW% protection",
                             "IDW% planted, cut,\\n& burned forests", 
                             "IDW% aquatic", 
                             "Watershed size",
                             "Latitude",
                             "IDW% forest",
                             "Longitude",
                             "IDW%\\nanthropogenic",
                             "IDW% wetlands", 
                             "IDW% agricultural"))

ISBsWeightedRelImportance
ggsave(ISBsWeightedRelImportance, file="Figure4B_ISBsWeightedRelImportance.png", scale=1, dpi = 600)


#Build and Validate best fit model 
M1W1 <- lmer(ZlogitISVB ~ ZRelAgr +ZRelFor+ ZRelAnthro +ZRelWet + (1|SBWBasins), data=SB_All,REML= FALSE)

M1W2 <- lmer(ZlogitISVB ~ ZRelAgr +ZRelFor+ ZRelAnthro +ZRelWet + ZLongitude + (1|SBWBasins), data=SB_All,REML= FALSE)
anova(M1W1, M1W2)
#No statistical difference between models so pick the one with the least variables given principle of parsimony. 

#Check multicollinearity 
vif(M1W1)

#A) Look at homogeneity
ResidM1W1<-resid(M1W1)
FM1W1<-fitted(M1W1)
plot(x=FM1W1, y= ResidM1W1, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model
#Random factor
boxplot(ResidM1W1~SBWBasins, ylab= "Normalized residuals", data = SB_All, xlab = "Nested Basins")
abline(0,0, lty=2)

#In the model
plot(x=SB_All$ZRelAgr, y= ResidM1W1, xlab = "Percent Agricultural", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelWet, y= ResidM1W1, xlab = "Percent Wetland", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelAnthro, y= ResidM1W1, xlab = "Percent Anthropogenic", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelFor, y= ResidM1W, xlab = "Percent Forest", ylab = "Normalized resitduals")
abline(0,0, lty=2)

#Not in the model 
plot(x=SB_All$ZRelPro, y= ResidM1W1, xlab = "Percent Protection", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZFlowAccumu, y= ResidM1W1, xlab = "Flow Accumulation", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelCut, y= ResidM1W1, xlab = "Percent Cut", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelBar, y= ResidM1W1, xlab = "Percent Barren Land", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelAqu, y= ResidM1W1, xlab = "Percent Aquatic", ylab = "Normalized residuals")
abline(0,0, lty=2)

#C) Look at normality of residuals
hist(ResidM1W1)
#this model meets its assumptions

#Extract R2 values and derive significance of variables 
summary(M1W1)
r.squaredGLMM(M1W1)

#significance of each variable

#Agri
MAgrW1 <- lmer(ZlogitISVB ~ ZRelFor+ ZRelAnthro +ZRelWet  + (1|SBWBasins), data=SB_All,REML= FALSE)
anova(M1W1, MAgrW1)
#For
MForW1 <- lmer(ZlogitISVB ~ ZRelAgr + ZRelAnthro +ZRelWet + (1|SBWBasins), data=SB_All,REML= FALSE)
anova(M1W1, MForW1)
#Ant
MAntW1 <- lmer(ZlogitISVB ~ ZRelAgr +ZRelFor +ZRelWet+ (1|SBWBasins), data=SB_All,REML= FALSE)
anova(M1W1, MAntW1)
#Wet
MWetW1 <- lmer(ZlogitISVB ~ ZRelAgr +ZRelFor+ ZRelAnthro+ (1|SBWBasins), data=SB_All,REML= FALSE)
anova(M1W1, MWetW1)

#ISBs - T test comparing watersheds with and without protection####
#Visualize Water Quality at sites without and without protection  
#Boxplot with 2 levels showing differences
#Figure6A####
#Figure6A_ISBsProCatViolin
ISBs_ProCatViolin <- ggplot(SB_All, aes(x=ProtectionCategorical, y =ISVB)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point", shape=19, size=3.5) +
  geom_hline(yintercept = c(45, 75), colour = c("black", "black"), linetype = "dashed") +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=16,face="bold"), 
        panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title=element_text(size=15), legend.text=element_text(size=14))+
  annotate("text", x = 2.7, y = 35, angle = 360, size = 6, label = "Bad")+
  annotate("text", x = 2.7, y = 60, angle = 360, size = 6,label = "Precarious   ")+
  annotate("text", x = 2.7, y = 90, angle = 360, size = 6,label = "Good")+
  annotate("text", x = 3, y = 40, angle = 360, size = 6,label = "")+
  
  xlab("Protection status") +
  ylab("Water quality\\n(ISBs indicator)")
ISBs_ProCatViolin

#A) compare water quality between protected and unprotected sites 
#Welch's T-test (allows for unequal sample sizes)
t.test(SB_All$logitISVB~ SB_All$ProtectionCategorical, var.equal=FALSE)

#Verify model assumptions
levels(SB_All$ProtectionCategorical)
SB_All$ProtectionCategorical <- as.factor(SB_All$ProtectionCategorical)
with(SB_All, shapiro.test(logitISVB[ProtectionCategorical == "No IUCN Category\\nII protection"]))
with(SB_All, shapiro.test(logitISVB[ProtectionCategorical == "Some level of IUCN\\nCategory II protection"]))
#yes

#Do the two populations have the same variance (although this isn't technically an assumption of the Welch's test)
var.test(logitISVB ~ ProtectionCategorical, data = SB_All)
# p =0.06 (so just passes the test)


#calculate effect size
#https://cran.r-project.org/web/packages/effsize/effsize.pdf
cohen.d(SB_All$logitISVB,SB_All$ProtectionCategorical)

#B) determine if  water quality related to % natural land cover within the watershed?

cor(SB_All4Corr, use="complete.obs", method = "kendall")
#No variables (aside from forest) too highly correlated with percent nat land to keep them in the model 

SB_NatLand <- lmer(ZlogitISVB ~ ZRelPro +ZRelAgr  +ZRelPercentNatLand + 
                     ZRelCut + ZRelAnthro+  ZFlowAccumu +ZLatitude +ZLongitude + (1|SBWBasins), data=SB_All,REML= FALSE)

#Select best fit model
options(na.action = "na.fail") 
resultsSB_NatLand  <- dredge(SB_NatLand )
options(na.action = "na.omit") 
top.modelsSB_NatLand  <- get.models(resultsSB_NatLand , subset = delta < 5)
top.modelsSB_NatLand [1]
model.avg(top.modelsSB_NatLand )
aic_impSB_NatLand  <- importance(top.modelsSB_NatLand)
aic_impSB_NatLand  <- as.data.frame((aic_impSB_NatLand))
aic_impSB_NatLand $variables <- rownames(aic_impSB_NatLand)
aic_impSB_NatLand $VariableNumbers <- c("K","J","I","H","G","F","E","D","C","B","A")

#Stick with full model to estimate parameters
#Check multicollinearity 
vif(SB_NatLand)

#A) Look at homogeneity
ResidSB_NatLand<-resid(SB_NatLand)
FitSB_NatLand<-fitted(SB_NatLand)
plot(x=FitSB_NatLand, y= ResidSB_NatLand, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model
#Random factor
boxplot(ResidSB_NatLand~SBWBasins, ylab= "Normalized residuals", data = SB_All, xlab = "Nested Basins")
abline(0,0, lty=2)

#In the model
plot(x=SB_All$ZRelAgr, y= ResidSB_NatLand, xlab = "Percent Agricultural", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelAnthro, y= ResidSB_NatLand, xlab = "Percent Anthropogenic", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelPro, y= ResidSB_NatLand, xlab = "Percent Protection", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZFlowAccumu, y= ResidSB_NatLand, xlab = "Flow Accumulation", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelCut, y= ResidSB_NatLand, xlab = "Percent Cut", ylab = "Normalized residuals")
abline(0,0, lty=2)

#C) Look at normality of residuals
hist(ResidSB_NatLand)
#this model meets its assumptions

#Extract R2 values and derive significance of variables 
summary(SB_NatLand)
r.squaredGLMM(SB_NatLand)
coef(SB_NatLand)

#significance of natural land cover 

#Agri
SB_NatLandNoNat <-  lmer(ZlogitISVB ~ ZRelPro +ZRelAgr  +
                           #ZRelPercentNatLand + 
                              ZRelCut + ZRelAnthro+  ZFlowAccumu +ZLatitude +ZLongitude + (1|SBWBasins), data=SB_All,REML= FALSE)


anova(SB_NatLand,SB_NatLandNoNat)

#Looking at natural land cover only 
SB_NatLandOnly<-  lmer(ISVB ~ 
                           RelPercentNatLand + 
                            (1|SBWBasins), data=SB_All,REML= FALSE)

#A) Look at homogeneity
ResidSB_NatLandOnly<-resid(SB_NatLandOnly)
FitSB_NatLandOnly<-fitted(SB_NatLandOnly)
plot(x=FitSB_NatLandOnly, y= ResidSB_NatLandOnly, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model
#Random factor
boxplot(ResidSB_NatLandOnly~SBWBasins, ylab= "Normalized residuals", data = SB_All, xlab = "Nested Basins")
abline(0,0, lty=2)


#C) Look at normality of residuals
hist(ResidSB_NatLandOnly)
#this model meets its assumptions

#Model results
summary(SB_NatLandOnly)
r.squaredGLMM(SB_NatLandOnly)

#Significance of NatLandCover
SB_NatLandRandomEffect<-  lmer(ISVB ~ 
                        1 + 
                         (1|SBWBasins), data=SB_All,REML= FALSE)

anova(SB_NatLandOnly, SB_NatLandRandomEffect)

#visualize regression
#Figure 6C####
ISBs_NatLand<-ggplot(SB_All, aes(x=RelPercentNatLand, y=ISVB)) +
  geom_point() +
  #geom_segment(aes(x =0, xend = 100, y = -0.001995 + 0.077324, yend = -0.001995 + 0.077324*100)) +
  #geom_segment(aes(x =0, xend = 100, y = -0.005835 + 0.383704, yend = -0.005835 + 0.383704*100)) +
  #geom_segment(aes(x =0, xend = 100, y = 0.480853 + 0.013357, yend = 0.480853 + 0.013357*100)) +
  #Intercepts from model with untransformed data
  geom_segment(aes(x =0, xend = 100, y = 60.17434 + 0.24260, yend = 60.17434 + 0.24260*100)) +
  geom_hline(yintercept = c(45, 75), colour = c("black", "black"), linetype = "dashed") +
  theme(axis.text=element_text(size=13),axis.title=element_text(size=14,face="bold"), 
        panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title=element_text(size=15), legend.text=element_text(size=14))+
  annotate("text", x = 105, y = 35, angle = 360, size = 4, label = "Bad")+
  annotate("text", x = 105, y = 60, angle = 360, size = 4,label = "Precarious  ")+
  annotate("text", x = 105, y = 90, angle = 360, size = 4,label = "Good")+  
  annotate("text", x = 108, y = 40, angle = 360, size = 4,label = "")+  
  
  xlab("Inverse distance weighted % upland natural land cover") +
  ylab("Water quality\\n(ISBs indicator)")
ISBs_NatLand

#Figure6C_SBNatLandCover
#C) compare natural land cover between protected and unprotected watersheds 
#Visualize 
#Figure6B####
ISBs_NatLandViolin<- ggplot(SB_All, aes(x=ProtectionCategorical, y = RelPercentNatLand)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point", shape=19, size=3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=16,face="bold"),
        panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title=element_text(size=15), legend.text=element_text(size=14))+
  annotate("text", x = 2.5, y = 101, angle = 90, label = " ") +
  xlab("Protection status") +
  ylab("Inverse distance weighted\\n % upland natural land cover")
ISBs_NatLandViolin

#compare statistically 
#T-test (Welch's Test due to different group sizes)
t.test(SB_All$RelPercentNatLand~ SB_All$ProtectionCategorical, var.equal=FALSE)

#Verify model assumptions
#is the data from both groups normally distributed 
hist(with(SB_All,RelPercentNatLand[ProtectionCategorical == "No IUCN Category\\nII Protection"]))
with(SB_All, shapiro.test(RelPercentNatLand[ProtectionCategorical == "No IUCN Category\\nII Protection"]));
#No, so move on to a non-parametric test 


#Mann_Whitney U Test
wilcox.test(RelPercentNatLand ~ ProtectionCategorical, data = SB_All)
#There is a significant difference 

#Compute effect size
VD.A(d = SB_All$RelPercentNatLand,
     f = SB_All$ProtectionCategorical)

#Boxplot to illustrate differences in means 
ISVB_ProCat2 <- ggplot(SB_All, aes(x=ProtectionCategorical, y =ISVB)) +
  geom_boxplot() 
ISVB_ProCat2

t.test(SB_All$logitISVB ~ SB_All$ProtectionCategorical)

#Verify model assumptions
#is the data from both groups normally distributed 


#Chi Squared test - ISBs####
#http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r

#on sites with protection
#Reorder factors if necessary
SB_ProSubset$WaterQualityCategory <- as.factor(SB_ProSubset$WaterQualityCategory)
levels(SB_ProSubset$WaterQualityCategory)

#Chi-Square Test
SBProchisq<- chisq.test(SB_ProSubset$ProtectionCategory,SB_ProSubset$WaterQualityCategory, correct= FALSE)
SBProchisq

SB_ProSubsetTable <- SB_ProSubset %>%
  group_by(ProtectionCategory, WaterQualityCategory) %>%
  summarize(n = n()) %>%
  mutate(freq = n/sum(n))

#Contingency Table 
SB_ProSubsetTable<-as.data.frame(SB_ProSubsetTable)
SB_ProSubsetTable4Spread<- SB_ProSubsetTable[,-3]
SB_ProSubsetContingencyTable<-spread(SB_ProSubsetTable4Spread,ProtectionCategory, freq)
SB_ProSubsetContingencyTable[is.na(SB_ProSubsetContingencyTable)] <- 0
SB_ProSubsetContingencyTable

write.csv(SB_ProSubsetContingencyTable, "SB_ProSubsetContingencyTable.csv")

#Comparision between protected and unprotected sites
SBchisqProNot<- chisq.test(SB_All$ProtectionCategorical,SB_All$WaterQualityCategory, correct= FALSE)
SBchisqProNot 

#Small Contingency Table 
SB_SmallTable <- SB_All %>%
  group_by(ProtectionCategorical, WaterQualityCategory) %>%
  summarize(n = n()) %>%
  mutate(freq = n/sum(n))

#Contingency Table 
SB_SmallTable<-as.data.frame(SB_SmallTable)
SB_SmallTable4Spread<- SB_SmallTable[,-3]
SB_SmallContingencyTable<-spread(SB_SmallTable4Spread,ProtectionCategorical, freq)
SB_SmallContingencyTable[is.na(SB_SmallContingencyTable)] <- 0
SB_SmallContingencyTable

write.csv(SB_SmallContingencyTable, "SB_SmallContingencyTable.csv")


##All data#
SB_Table <- SB_All %>%
  group_by(ProtectionCategory, WaterQualityCategory) %>%
  summarize(n = n()) %>%
  mutate(freq = n/sum(n))

SBchisq<- chisq.test(SB_All$WaterQualityCategory,SB_All$ProtectionCategory, correct= FALSE)
SBchisq

SB_MeanWQ <- SB_All %>%
  group_by(ProtectionCategory) %>%
  summarize(meanWQ = mean(ISVB), minWQ = min(ISVB), maxWQ = max(ISVB), seWQ = se(ISVB))



#Supporting Information####
#SI - IQBP Analyses#### 
#Section 3.2 - Analysis of subsets of data####
#Subsets with same amount of protected and unprotected watersheds####
#Look at data distribution
hist(IQBP_EqualProAndNot$IQBP_mean)
hist(IQBP_EqualProAndNot$logitIQBPmean)

#Assess overall model 
M1IQBPEqualProAndNot <- lmer(ZlogitIQBP ~ ZRelPro + ZRelAgr +ZRelFor+  ZRelWet + 
                               ZRelCut + ZRelAnt + ZRelAqu + ZRelBar + ZFlowAccumu + 
                               ZLatitude + ZLongitude + (1|IQBPBasins), data=IQBP_EqualProAndNot,REML= FALSE)

#Singular fit so we move onto a model that only has protected sites 
M1IQBPEqualProAndNot <- lmer(ZlogitIQBP ~ ZRelPro  + (1|IQBPBasins), data=IQBP_EqualProAndNot,REML= FALSE)

#Validate the model 
#A) Look at homogeneity
ResidM1IQBPEqualProAndNot<-resid(M1IQBPEqualProAndNot)
FitM1IQBPEqualProAndNot<-fitted(M1IQBPEqualProAndNot)
plot(x=FitM1IQBPEqualProAndNot, y= ResidM1IQBPEqualProAndNot, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)
#Normal that we see a trend here given that there are missing covariates that do explain a good proportion of variation 

#B) Look at independence
#Plot residuals vs each covariate 

#Random factor
boxplot(ResidM1IQBPEqualProAndNot~IQBPBasins, ylab= "Normalized residuals", data = IQBP_EqualProAndNot, xlab = "Nested Basins")
abline(0,0, lty=2)


#C) Look at normality of residuals
hist(ResidM1IQBPEqualProAndNot)

#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(M1IQBPEqualProAndNot)

summary(M1IQBPEqualProAndNot)

#Compare to a model without protection
M1IQBPEqualProAndNotNoPro <- lmer(ZlogitIQBP ~ 1  + (1|IQBPBasins), data=IQBP_EqualProAndNot,REML= FALSE)
anova(M1IQBPEqualProAndNotNoPro, M1IQBPEqualProAndNot)

#Try the same model including agriculture 
M1IQBPEqualProAndNot2 <- lmer(ZlogitIQBP ~ ZRelPro + ZRelAgr + (1|IQBPBasins), data=IQBP_EqualProAndNot,REML= FALSE)

#Validate the model 
#A) Look at homogeneity
ResidM1IQBPEqualProAndNot2<-resid(M1IQBPEqualProAndNot2)
FitM1IQBPEqualProAndNot2<-fitted(M1IQBPEqualProAndNot2)
plot(x=FitM1IQBPEqualProAndNot2, y= ResidM1IQBPEqualProAndNot2, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)
#THis model does meet this assumption, unlike the one that only included protection 

#B) Look at independence
#Plot residuals vs each covariate 

#Random factor
boxplot(ResidM1IQBPEqualProAndNot2~IQBPBasins, ylab= "Normalized residuals", data = IQBP_EqualProAndNot, xlab = "Nested Basins")
abline(0,0, lty=2)


#C) Look at normality of residuals
hist(ResidM1IQBPEqualProAndNot2)

#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(M1IQBPEqualProAndNot2)

summary(M1IQBPEqualProAndNot2)
#Protection is no longer considered significant, which demonstrates its lack of importance as a predictor variable 



#Look at natural land cover model 
IQBPNatEqualProAndNot <- lmer(ZlogitIQBP ~ ZRelPercentNatLand  + ZRelCut +ZRelAnt + 
                                ZFlowAccumu + ZLatitude + ZLongitude + (1|IQBPBasins), data=IQBP_EqualProAndNot,REML= FALSE)

#Validate  model 
vif(IQBPNatEqualProAndNot)
#VIFs below 5 are considered acceptable 

#A) Look at homogeneity
ResidIQBPNatEqualProAndNot<-resid(IQBPNatEqualProAndNot)
FitIQBPNatEqualProAndNot<-fitted(IQBPNatEqualProAndNot)
plot(x=FitIQBPNatEqualProAndNot, y= ResidIQBPNatEqualProAndNot, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model 
#Random factor
boxplot(ResidIQBPNatEqualProAndNot~IQBPBasins, ylab= "Normalized residuals", data = IQBP_EqualProAndNot, xlab = "Nested Basins")
abline(0,0, lty=2)

#Fixed factors in the model
plot(x=IQBP_EqualProAndNot$ZRelAnt, y= ResidIQBPNatEqualProAndNot, xlab = "Percent Upland Agricultural (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_EqualProAndNot$ZFlowAccumu, y= ResidIQBPNatEqualProAndNot, xlab = "Percent Upland Agricultural (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)

#Fixed factors not in the model
plot(x=IQBP_EqualProAndNot$ZRelAgr, y= ResidIQBPNatEqualProAndNot, xlab = "Percent Upland Agricultural (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)

#C) Look at normality of residuals
hist(ResidIQBPNatEqualProAndNot)


#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(IQBPNatEqualProAndNot)

#Assess significance of fixed factors 
summary(IQBPNatEqualProAndNot)


#p value for Natural Land Cover
IQBPNatEqualProAndNotNoNat <- lmer(ZlogitIQBP ~ ZRelCut +ZRelAnt + 
                                     ZFlowAccumu + ZLatitude + ZLongitude + (1|IQBPBasins), data=IQBP_EqualProAndNot,REML= FALSE)

anova(IQBPNatEqualProAndNot, IQBPNatEqualProAndNotNoNat)


#Chi-square Comparison between protected sites and a subset of unprotected sites of equal size (i.e 37 locations per category)
#Comparision between protected and unprotected sites 
IQBPchisqEqualProNot<- chisq.test(IQBP_EqualProAndNot$ProtectionCategorical,IQBP_EqualProAndNot$WaterQualityCategory, correct= FALSE)
IQBPchisqEqualProNot 
#Just significant 

#Comparision between protected and unprotected sites using different protection category 
IQBPchisqEqualProNot2<- chisq.test(IQBP_EqualProAndNot$ProtectionCategory,IQBP_EqualProAndNot$WaterQualityCategory, correct= FALSE)
IQBPchisqEqualProNot2 
#Not significant 


#Subset only including Watersheds with protection ####
#Mixed model 
hist(IQBP_ProSubset$IQBP_mean)
hist(IQBP_ProSubset$logIQBPmean)
hist(IQBP_ProSubset$logitIQBPmean)

#Check for correlation among variables 
IQBP_ProSubsetCorr <- select(IQBP_ProSubset, "ZlogitIQBP", "ZRelPro", "ZFlowAccumu", "ZRelWet", "ZRelBar", "ZRelAqu", 
                             "ZRelCut", "ZRelFor", "ZRelAnt", "ZRelAgr", "ZLatitude", "ZLongitude", "ZRelPercentNatLand")
cor(IQBP_ProSubsetCorr)
#Agriculture and forest are 89% correlated so both cannot be included in the model 

#Select the best-fit model among the two 
MPro_IQBPAgr <- lmer(ZlogitIQBP ~ ZRelPro + 
                       ZRelAgr +
                       # ZRelFor+  
                       ZRelWet + 
                       ZRelCut + ZRelAnt + ZRelAqu + ZRelBar + ZFlowAccumu + 
                       ZLatitude + ZLongitude + (1|IQBPBasins), data=IQBP_ProSubset,REML= FALSE)
summary(MPro_IQBPAgr)
#AIC 51.8
r.squaredGLMM(MPro_IQBPAgr)
#R2: 0.816, 0.883

MPro_IQBPFor <- lmer(ZlogitIQBP ~ ZRelPro + 
                       #ZRelAgr +
                       ZRelFor+  
                       ZRelWet + 
                       ZRelCut + ZRelAnt + ZRelAqu + ZRelBar + ZFlowAccumu + 
                       ZLatitude + ZLongitude + (1|IQBPBasins), data=IQBP_ProSubset,REML= FALSE)
summary(MPro_IQBPFor)
#AIC 51.3
r.squaredGLMM(MPro_IQBPFor)
#R2: 0.823, 0.881


#Model averaging/selection
options(na.action = "na.fail") 
resultsMPro_IQBPFor <- dredge(MPro_IQBPFor)
options(na.action = "na.omit") 
#Many models are showing singular fits, meaning that they are overfitted and should ideally not be considered. 
top.modelsMPro_IQBPFor<- get.models(resultsMPro_IQBPFor, subset = delta < 5)
top.modelsMPro_IQBPFor[1]
model.avg(top.modelsMPro_IQBPFor)

aic_impMPro_IQBPFor <- importance(top.modelsMPro_IQBPFor)
aic_impMPro_IQBPFor <- as.data.frame((aic_impMPro_IQBPFor))
aic_impMPro_IQBPFor$variables <- rownames(aic_impMPro_IQBPFor)
aic_impMPro_IQBPFor$VariableNumbers <- c("J","I","H","G","F","E","D","C","B","A")
aic_impMPro_IQBPFor


#Build and Validate best fit model 
MPro_IQBP_Best<- lmer(ZlogitIQBP ~  ZRelAqu + ZRelBar + ZRelFor +
                        (1|IQBPBasins), 
                      data=IQBP_ProSubset,REML= FALSE)


#Validate the model 
vif(MPro_IQBPFor_Best)


#A) Look at homogeneity
ResidMPro_IQBP_Best<-resid(MPro_IQBP_Best)
FitMPro_IQBP_Best<-fitted(MPro_IQBP_Best)
plot(x=FitMPro_IQBP_Best, y= ResidMPro_IQBP_Best, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model

#Random factor
boxplot(ResidMPro_IQBP_Best~IQBPBasins, ylab= "Normalized residuals", data = IQBP_ProSubset, xlab = "Nested Basins")
abline(0,0, lty=2)


#C) Look at normality of residuals
hist(ResidMPro_IQBP_Best)

#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(MPro_IQBP_Best)

summary(MPro_IQBP_Best)


#Subsets of watershed sizes ####
#Explore data distribution by watershed sizes
IQBP_Watersheds10000 <- IQBP_FullSummary %>% filter(WatersheadAreakm2 <= 10000)
hist(IQBP_Watersheds10000$WatersheadAreakm2)
#176/183 (96%)

IQBP_Watersheds1000<- IQBP_FullSummary %>% filter(WatersheadAreakm2 <= 1000)
hist(IQBP_Watersheds1000$WatersheadAreakm2)
#122/183 (67%)

IQBP_WatershedsBelow100<- IQBP_FullSummary %>% filter(WatersheadAreakm2 <= 100)
hist(IQBP_WatershedsBelow100$WatersheadAreakm2)
#30/183 (16%)

IQBP_Watersheds100to1000<- IQBP_FullSummary %>% filter(WatersheadAreakm2 >= 100 & WatersheadAreakm2 <= 1000)
hist(IQBP_Watersheds100to1000$WatersheadAreakm2)
#92/183 (50%)

IQBP_Watersheds1000to10000<- IQBP_FullSummary %>% filter(WatersheadAreakm2 >= 1000 & WatersheadAreakm2 <= 10000)
hist(IQBP_Watersheds1000to10000$WatersheadAreakm2)
#54/183 (30%)

IQBP_Watersheds100to500<- IQBP_FullSummary %>% filter(WatersheadAreakm2 >= 100 & WatersheadAreakm2 <= 500)
hist(IQBP_Watersheds100to500$WatersheadAreakm2)
#60/183 (33%)

#Assess relationship in watersheds between 100 and 1000km2 (50% of watersheds)
hist(IQBP_Watersheds100to1000$IQBP_mean)
hist(IQBP_Watersheds100to1000$logitIQBPmean)

IQBP_Watersheds100to1000Cor <- IQBP_Watersheds100to1000 %>% select(ZRelPro, ZRelAgr, ZRelFor,  ZRelWet, 
                                                                   ZRelCut, ZRelAnt, ZRelAqu, ZRelBar, ZFlowAccumu, 
                                                                   ZLatitude, Longitude)
cor(IQBP_Watersheds100to1000Cor)
#No correlations above  49% so we're good to proceed 

MW100to1000_IQBP <- lmer(ZlogitIQBP ~ ZRelPro + ZRelAgr +ZRelFor+  ZRelWet + 
                           ZRelCut + ZRelAnt + ZRelAqu + ZRelBar + 
                           ZFlowAccumu + 
                           ZLatitude + ZLongitude + (1|IQBPBasins), data=IQBP_Watersheds100to1000,REML= FALSE)


#Assess the relative importance of the variables included in the model 
options(na.action = "na.fail") 
resultsMW100to1000_IQBP <- dredge(MW100to1000_IQBP)
#Some of the models have a singular fit 

options(na.action = "na.omit") 
top.modelsIQBPMW100to1000<- get.models(resultsMW100to1000_IQBP, subset = delta < 5)
top.modelsIQBPMW100to1000[1]
model.avg(top.modelsIQBPMW100to1000)
#Several variables have a coef above 1 but this isn't a problem if the Vifs of the best fit model we end up using is OK
aic_impIQBPMW100to1000 <- importance(top.modelsIQBPMW100to1000)
aic_impIQBPMW100to1000 <- as.data.frame((aic_impIQBPMW100to1000))
aic_impIQBPMW100to1000$variables <- rownames(aic_impIQBPMW100to1000)
aic_impIQBPMW100to1000$VariableNumbers <- c("K","J","I","H","G","F","E","D","C","B","A")


#Build and Validate best fit model 
M1W100to1000_IQBP<- lmer(ZlogitIQBP ~  ZLatitude + 
                           ZRelAnt +
                           ZRelAqu +
                           ZRelBar +
                           ZRelCut +
                           ZRelFor +
                           (1|IQBPBasins), 
                         data=IQBP_Watersheds100to1000,REML= FALSE)


#Validate the model 
vif(M1W100to1000_IQBP)

#A) Look at homogeneity
ResidM1W100to1000_IQBP<-resid(M1W100to1000_IQBP)
FM1W100to1000_IQBP<-fitted(M1W100to1000_IQBP)
plot(x=FM1W100to1000_IQBP, y= ResidM1W100to1000_IQBP, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model

#Random factor
boxplot(ResidM1W100to1000_IQBP~IQBPBasins, ylab= "Normalized residuals", data = IQBP_Watersheds100to1000, xlab = "Nested Basins")
abline(0,0, lty=2)


#C) Look at normality of residuals
hist(ResidM1W100to1000_IQBP)

#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(M1W100to1000_IQBP)

summary(M1W100to1000_IQBP)

#Assess relationship in watersheds between 1000km2 and 10000km2 (30% of data)
hist(IQBP_Watersheds1000to10000$IQBP_mean)
hist(IQBP_Watersheds1000to10000$logitIQBPmean)

IQBP_Watersheds1000to10000Cor <- IQBP_Watersheds1000to10000 %>% select(ZRelPro, ZRelAgr, ZRelFor,  ZRelWet, 
                                                                       ZRelCut, ZRelAnt, ZRelAqu, ZRelBar, ZFlowAccumu, 
                                                                       ZLatitude, Longitude)
cor(IQBP_Watersheds1000to10000Cor)
#Forest and agriculture are over 80% correlated so one needs to be removed from the model 

#Try model with Forest
MW1000to10000_IQBP_For <- lmer(ZlogitIQBP ~ ZRelPro + 
                                 # ZRelAgr +
                                 ZRelFor+  
                                 ZRelWet + 
                                 ZRelCut + ZRelAnt + ZRelAqu + ZRelBar + 
                                 ZFlowAccumu + 
                                 ZLatitude + ZLongitude + (1|IQBPBasins), data=IQBP_Watersheds1000to10000,REML= FALSE)

#Try model with Agriculture
MW1000to10000_IQBP_Agr <- lmer(ZlogitIQBP ~ ZRelPro + 
                                 ZRelAgr +
                                 #ZRelFor+  
                                 ZRelWet + 
                                 ZRelCut + ZRelAnt + ZRelAqu + ZRelBar + 
                                 ZFlowAccumu + 
                                 ZLatitude + ZLongitude + (1|IQBPBasins), data=IQBP_Watersheds1000to10000,REML= FALSE)



#Compare their AIC and R2 values to decide which one to retain 
summary(MW1000to10000_IQBP_For) 
r.squaredGLMM(MW1000to10000_IQBP_For)
#AIC = 66.1,
#R2 - 0.80, 0.85

summary(MW1000to10000_IQBP_Agr)
r.squaredGLMM(MW1000to10000_IQBP_Agr)
#AIC = 67.5
#R2 - 0.78, 0.85

#Because the model with forest has a slightly lower AIC and higher R2 we move forward with this model 

#Assess the relative importance of the variables included in the model 
options(na.action = "na.fail") 
resultsMW1000to10000_IQBP_For <- dredge(MW1000to10000_IQBP_For)
options(na.action = "na.omit") 
top.modelsIQBPMW1000to10000For<- get.models(resultsMW1000to10000_IQBP_For, subset = delta < 5)
top.modelsIQBPMW1000to10000For[1]
model.avg(top.modelsIQBPMW1000to10000For)

aic_impIQBPMW1000to10000For <- importance(top.modelsIQBPMW1000to10000For)
aic_impIQBPMW1000to10000For <- as.data.frame((aic_impIQBPMW1000to10000For))
aic_impIQBPMW1000to10000For$variables <- rownames(aic_impIQBPMW1000to10000For)
aic_impIQBPMW1000to10000For$VariableNumbers <- c("K","J","I","H","G","F","E","D","C","B","A")


#Build and Validate best fit model 
MW1000to10000_IQBP_ForBest<- lmer(ZlogitIQBP ~  
                                    ZFlowAccumu +
                                    ZRelAqu +
                                    ZRelBar +
                                    ZRelCut +
                                    ZRelFor +
                                    (1|IQBPBasins), 
                                  data=IQBP_Watersheds1000to10000,REML= FALSE)


#Validate the model 
vif(MW1000to10000_IQBP_ForBest)


#A) Look at homogeneity
ResidMW1000to10000_IQBP_ForBest<-resid(MW1000to10000_IQBP_ForBest)
FitMW1000to10000_IQBP_ForBest<-fitted(MW1000to10000_IQBP_ForBest)
plot(x=FitMW1000to10000_IQBP_ForBest, y= ResidMW1000to10000_IQBP_ForBest, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model

#Random factor
boxplot(ResidMW1000to10000_IQBP_ForBest~IQBPBasins, ylab= "Normalized residuals", data = IQBP_Watersheds1000to10000, xlab = "Nested Basins")
abline(0,0, lty=2)


#C) Look at normality of residuals
hist(ResidMW1000to10000_IQBP_ForBest)

#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(MW1000to10000_IQBP_ForBest)

summary(MW1000to10000_IQBP)


#Model with all variables
MW1000to10000_IQBP <- lmer(ZlogitIQBP ~ ZRelPro + 
                             ZRelAgr +
                             ZRelFor+  
                             ZRelWet + 
                             ZRelCut + ZRelAnt + ZRelAqu + ZRelBar + 
                             ZFlowAccumu + 
                             ZLatitude + ZLongitude + (1|IQBPBasins), data=IQBP_Watersheds1000to10000,REML= FALSE)

#Assess the relative importance of the variables included in the model 
options(na.action = "na.fail") 
resultsMW1000to10000_IQBP <- dredge(MW1000to10000_IQBP)
#Two models have singular fits, but these may not be problematic if they are not in the top models 

options(na.action = "na.omit") 
top.modelsIQBPMW1000to10000<- get.models(resultsMW1000to10000_IQBP, subset = delta < 5)
top.modelsIQBPMW1000to10000[1]
model.avg(top.modelsIQBPMW1000to10000)

aic_impIQBPMW1000to10000 <- importance(top.modelsIQBPMW1000to10000)
aic_impIQBPMW1000to10000 <- as.data.frame((aic_impIQBPMW1000to10000))
aic_impIQBPMW1000to10000$variables <- rownames(aic_impIQBPMW1000to10000)
aic_impIQBPMW1000to10000$VariableNumbers <- c("K","J","I","H","G","F","E","D","C","B","A")


#Build and Validate best fit model 
MW1000to10000_IQBP_Best<- lmer(ZlogitIQBP ~  
                                 ZFlowAccumu +
                                 ZRelAqu +
                                 ZRelBar +
                                 ZRelCut +
                                 ZRelFor +
                                 (1|IQBPBasins), 
                               data=IQBP_Watersheds1000to10000,REML= FALSE)


#Validate the model 
vif(MW1000to10000_IQBP_Best)


#A) Look at homogeneity
ResidMW1000to10000_IQBP_Best<-resid(MW1000to10000_IQBP_Best)
FitMW1000to10000_IQBP_Best<-fitted(MW1000to10000_IQBP_Best)
plot(x=FitMW1000to10000_IQBP_Best, y= ResidMW1000to10000_IQBP_Best, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model

#Random factor
boxplot(ResidMW1000to10000_IQBP_Best~IQBPBasins, ylab= "Normalized residuals", data = IQBP_Watersheds1000to10000, xlab = "Nested Basins")
abline(0,0, lty=2)


#C) Look at normality of residuals
hist(ResidMW1000to10000_IQBP_Best)

#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(MW1000to10000_IQBP_Best)

summary(MW1000to10000_IQBP)



#Watersheds between 100 and 10 000km2 (all watersheds of 1000km2 plus or minus 1 order of magnitude - 
#reduces the dataset from 183 observations to 146 observations
IQBP_Watersheds1000 <- IQBP_FullSummary %>% filter(WatersheadAreakm2 >= 100 & WatersheadAreakm2 <= 10000)
hist(IQBP_Watersheds1000$WatersheadAreakm2)

hist(IQBP_Watersheds1000$ZlogitIQBP)

MW1000_IQBP <- lmer(ZlogitIQBP ~ ZRelPro + ZRelAgr +ZRelFor+  ZRelWet + 
                      ZRelCut + ZRelAnt + ZRelAqu + ZRelBar + ZFlowAccumu + 
                      ZLatitude + ZLongitude + (1|IQBPBasins), data=IQBP_Watersheds1000,REML= FALSE)
summary(MW1000_IQBP)

#Assess the relative importance of the variables included in the model 
options(na.action = "na.fail") 
resultsMW1000_IQBP <- dredge(MW1000_IQBP)
options(na.action = "na.omit") 
top.modelsIQBPMW1000<- get.models(resultsMW1000_IQBP, subset = delta < 5)
top.modelsIQBPMW1000[1]
model.avg(top.modelsIQBPMW1000)
aic_impIQBPMW100 <- importance(top.modelsIQBPMW1000)
aic_impIQBPMW100 <- as.data.frame((aic_impIQBPMW100))
aic_impIQBPMW100$variables <- rownames(aic_impIQBPMW100)
aic_impIQBPMW100$VariableNumbers <- c("K","J","I","H","G","F","E","D","C","B","A")

# Plotting models Delta AIC < 5
IQBPWeightedRelImportance_W1000 <-   ggplot(data=aic_impIQBPMW100, aes(y=aic_imp[,1],x=VariableNumbers)) +
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() + xlab("Predictor Variable") + ylab("Relative importance")  +
  theme(axis.text=element_text(size=10.5), axis.title=element_text(size=12,face="bold")) +
  #Add information about model estimates 
  annotate("text", x = 12.5, y = 1.15, angle = 0, label = "")+
  annotate("text", x = 12, y = 1.08, angle = 0, label = "Average\\nSlope")+
  annotate("text", x = 11.14, y = 1.08, angle = 0, label = "      0.30 + *")+ #Aqua
  annotate("text", x = 10.15, y = 1.08, angle = 0, label = "     - 0.13 + *")+ #Bar
  annotate("text", x = 9.15, y = 1.08, angle = 0, label = "       0.23 + *")+ #Cut
  annotate("text", x = 8.15, y = 1.08, angle = 0, label = "   0.73 + * ")+ #For
  annotate("text", x = 7.15, y = 1.08, angle = 0, label = "     0.21 + *") + #Ant
  annotate("text", x = 6.15, y = 1.08, angle = 0, label = " 0.07") + # Longitude
  annotate("text", x = 5.15, y = 1.08, angle = 0, label = " 0.09")+ #Latitute
  annotate("text", x = 4.15, y = 1.08, angle = 0, label = "- 0.14")+ #Agri
  annotate("text", x = 3.15, y = 1.08, angle = 0, label = " - 0.03") + #Wet
  annotate("text", x = 2.15, y = 1.08, angle = 0, label = " - 0.06")+ #Watershed Size
  annotate("text", x = 1.15, y = 1.08, angle = 0, label = "  0.01")+ #Protection 
  scale_x_discrete(labels= c("Percent Watershed Protection\\n(Distance Weighted)",
                             "Watershed size",
                             "Percent Watershed Wetlands\\n(Distance Weighted)", 
                             "Percent Watershed Agricultural\\n(Distance Weighted)",
                             "Latitude",
                             "Longitude",
                             "Percent Watershed Anthropogenic\\n(Distance Weighted)",
                             "Percent Watershed Forest\\n(Distance Weighted)",
                             "Percent Watershed\\nPlanted Forests, Cuts & Burns\\n(Distance Weighted)", 
                             "Percent Watershed Barren Land\\n(Distance Weighted)",
                             "Percent Watershed Aquatic\\n(Distance Weighted)"))

IQBPWeightedRelImportance_W1000
ggsave(ISBsWeightedRelImportance, file="FigureS5_IQBPW1000_WeightedRelImportance.png", scale=1, dpi = 600)


#Build and Validate best fit model 
M1W1000_IQBP<- lmer(ZlogitIQBP ~ ZRelAnt + ZRelAqu + ZRelBar + ZRelCut +ZRelFor +  (1|IQBPBasins), data=IQBP_Watersheds1000,REML= FALSE)


#Validate the model 
vif(M1W1000_IQBP)
#VIFs are below 5 this is valid 

#A) Look at homogeneity
ResidM1W1000_IQBP<-resid(M1W1000_IQBP)
FM1W1000_IQBP<-fitted(M1W1000_IQBP)
plot(x=FM1W1000_IQBP, y= ResidM1W1000_IQBP, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model

#Random factor
boxplot(ResidM1W1000_IQBP~IQBPBasins, ylab= "Normalized residuals", data = IQBP_Watersheds1000, xlab = "Nested Basins")
abline(0,0, lty=2)


#C) Look at normality of residuals
hist(ResidM1W1000_IQBP)

#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(M1W1000_IQBP)

summary(M1W1000_IQBP)


#Section 3.3 - Percent Watershed Protection Data Visualization####
#Figure S1####
IQBP_Pro <- ggplot(IQBP_FullSummary, aes(x=PercentPro, y=IQBP_mean)) +
  #geom_point(aes(shape = AgriculturalCategory), size= 3.5) +
  geom_point(size = 3, aes(color = ForestCategory, shape = ForestCategory)) +
  scale_color_manual(name="%\\nwatershed forest", 
                     breaks =c("0%", "1-24%", "25-49%","50-74%", "75-100%"),
                     labels = c("0%", "1-24%", "25-49%","50-74%", "75-100%"),
                     values=c("#a6611a", "#dfc27d",  "grey","#80dcd1","#018571")) +
  scale_shape_discrete(name="%\\nwatershed forest", 
                       breaks =c("0%", "1-24%", "25-49%","50-74%", "75-100%"),
                       labels = c("0%", "1-24%", "25-49%","50-74%", "75-100%"))+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=16,face="bold"), 
        panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.title=element_text(size=13, face = "bold"), legend.text=element_text(size=12)) +
  geom_hline(yintercept = c(19, 39, 59, 79), colour = c("black", "black", "black", "black"), linetype = "dashed") +
  #colour = c("black", "darkorange3", "darkorange", "green4"),
  annotate("text", x = 75, y = 8, angle = 360, label = "Very bad    ")+
  annotate("text", x = 75, y = 29, angle = 360, label = "Bad")+
  annotate("text", x = 75, y = 50, angle = 360, label = "Questionable             ")+
  annotate("text", x = 75, y = 70, angle = 360, label = "Satisfactory        ")+
  annotate("text", x = 75, y = 90, angle = 360, label = "Good")+
  xlab("% watershed protection") +
  ylab("Water quality\\n(IQBP indicator)")
IQBP_Pro

ggsave(IQBP_Pro, file="FigureS1_IQBPPro.tiff", dpi = 600)


#Section 3.4 - Percent Watershed Protection Mixed Model Analysis####
#Look at correlations between variables
IQBP_Corr2 <- select(IQBP_FullSummary, "ZlogitIQBP", "ZPercentPro", "ZFlowAccumu", "ZPercentWet", "ZPercentBar", 
                     "ZPercentAqu", "ZPercentCut", "ZPercentFor", "ZPercentAnt", "ZPercentAgg", "ZLatitude", "ZLongitude")
plot(IQBP_Corr2)
cor(IQBP_Corr2, use="complete.obs", method = "kendall")
#Agriculture is 63% correlated to IQBP, and forests are 48% correlated to IQBP. Everything else is lower

#Create a model that contains all variables 
M0_IQBP2 <- lmer(ZlogitIQBP ~ ZPercentPro + ZPercentAgg +ZPercentFor+  ZPercentWet + ZPercentCut + ZPercentAnt + 
                   ZPercentAqu + ZPercentBar + ZFlowAccumu + ZLatitude + ZLongitude + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)
summary(M0_IQBP2)


#Use the MuMIn package for model selection and averaging
options(na.action = "na.fail") 
resultsM0_IQBP2 <- dredge(M0_IQBP2)
options(na.action = "na.omit") 
top.modelsIQBP2 <- get.models(resultsM0_IQBP2, subset = delta < 5)
top.modelsIQBP2[1]
model.avg(top.modelsIQBP2)
aic_impIQBP2 <- importance(top.modelsIQBP2)
aic_impIQBP2 <- as.data.frame((aic_impIQBP2))
aic_impIQBP2$variables <- rownames(aic_impIQBP2)
aic_impIQBP2$VariableNumbers <- c("K","J","I","H","G","F","E","D","C","B","A")

# Plotting models Delta AIC < 5
#Fig. S3####
IQBPProRelImportance <-   ggplot(data=aic_impIQBP2, aes(y=aic_impIQBP[,1],x=VariableNumbers)) +
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() + xlab("Predictor variable") + ylab("Relative importance")  +
  theme(axis.text=element_text(size=12.5), axis.title=element_text(size=15,face="bold"),
        panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #Add information about model estimates 
  annotate("text", x = 12.5, y = 1.12, angle = 0, label = "")+
  annotate("text", x = 12, y = 1.08, angle = 0, label = "Average\\nSlope")+
  annotate("text", x = 11.15, y = 1.08, angle = 0, label = "      0.39 + *")+ #Forest
  annotate("text", x = 10.15, y = 1.08, angle = 0, label = "      0.36 + *")+ #Aquatic
  annotate("text", x = 9.15, y = 1.08, angle = 0, label = "       0.26 + *")+ #Cut
  annotate("text", x = 8.15, y = 1.08, angle = 0, label = "     0.15 + *") + #Longitude
  annotate("text", x = 7.15, y = 1.08, angle = 0, label = "    - 0.23 + *")+ #Agricultural
  annotate("text", x = 6.15, y = 1.08, angle = 0, label = "    - 0.08 + * ")+ #Barren land 
  annotate("text", x = 5.15, y = 1.08, angle = 0, label = " 0.14 ") + #Anthro
  annotate("text", x = 4.15, y = 1.08, angle = 0, label = "- 0.05  ") + #Watershed Size
  annotate("text", x = 3.15, y = 1.08, angle = 0, label = "- 0.05  ")+ #Latitute
  annotate("text", x = 2.15, y = 1.08, angle = 0, label = "  0.04  ")+ #Wetland
  annotate("text", x = 1.15, y = 1.08, angle = 0, label = "  0.04  ")+ #Protection
  scale_x_discrete(labels= c("% protection",
                             "% wetlands",
                             "Latitude",
                             "Watershed size",
                             "% anthropogenic",
                             "% rock and lichen", 
                             "% agricultural",
                             "Longitude", 
                             "% planted, cut,\\n& burned forests", 
                             "% aquatic",
                             "% forest"))

IQBPProRelImportance

ggsave(IQBPProRelImportance, file="FigureS3_IQBPRelImportance.pdf", scale=1, dpi = 600)

#Validate best fit model 
M1_IQBP2 <- lmer(ZlogitIQBP ~  ZPercentAgg +ZPercentFor + ZPercentCut + 
                   ZPercentAqu + ZLongitude + ZPercentBar + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)

vif(M1_IQBP2)
#All VIFs below 5 so this is OK 

par(mfrow=c(1,1)) 
#A) Look at homogeneity
ResidM1_IQBP2<-resid(M1_IQBP2)
FM1_IQBP2<-fitted(M1_IQBP2)
plot(x=FM1_IQBP2, y= ResidM1_IQBP2, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model

#Random factor
boxplot(ResidM1_IQBP2~IQBPBasins, ylab= "Normalized residuals", data = IQBP_FullSummary, xlab = "Nested Basins")
abline(0,0, lty=2)

#In the model
plot(x=IQBP_FullSummary$ZLongitude, y= ResidM1_IQBP2, xlab = "Longitude", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelAgr, y= ResidM1_IQBP2, xlab = "Percent Upland Agriculture (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelAqu, y= ResidM1_IQBP2, xlab = "Percent Upland Aquatic (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelCut, y= ResidM1_IQBP2, xlab = "Percent Upland Deforestation and Burns (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelFor, y= ResidM1_IQBP2, xlab = "Percent Upland Forest (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)

#Not in the model
plot(x=IQBP_FullSummary$ZRelPro, y= ResidM1_IQBP2, xlab = "Percent Upland Protected (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZLatitude, y= ResidM1_IQBP2, xlab = "Latitude", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelWet, y= ResidM1_IQBP2, xlab = "Percent Upland Wetland (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelAnt, y= ResidM1_IQBP2, xlab = "Percent Upland Anthropogenic (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZRelBar, y= ResidM1_IQBP2, xlab = "Percent Upland Barren Land and Rocks (weighted)", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=IQBP_FullSummary$ZFlowAccumu, y= ResidM1_IQBP2, xlab = "Watershed Size (Flow Accumulation)", ylab = "Normalized residuals")
abline(0,0, lty=2)

#C) Look at normality of residuals
hist(ResidM1_IQBP2)


#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(M1_IQBP2)

#Assess significance of fixed factors 
summary(M1_IQBP2)

#p value for each variable
#agriculture
MAgg_IQBP2<- lmer(ZlogitIQBP ~  ZPercentFor + ZPercentCut + 
                    ZPercentAqu + ZLongitude + ZPercentBar + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)
anova(M1_IQBP2, MAgg_IQBP2)

#Forest
MFor_IQBP2<- lmer(ZlogitIQBP ~  ZPercentAgg  + ZPercentCut + 
                    ZPercentAqu + ZLongitude + ZPercentBar + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)
anova(M1_IQBP2, MFor_IQBP2)
#Cut
MCut_IQBP2 <- lmer(ZlogitIQBP ~  ZPercentAgg +ZPercentFor +
                     ZPercentAqu + ZLongitude + ZPercentBar + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)
anova(M1_IQBP2, MCut_IQBP2)

#Aquat
MAqu_IQBP2 <- lmer(ZlogitIQBP ~  ZPercentAgg +ZPercentFor + ZPercentCut + 
                     ZLongitude + ZPercentBar + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)

anova(M1_IQBP2, MAqu_IQBP2)

#Longitude
MLon_IQBP2 <- lmer(ZlogitIQBP ~  ZPercentAgg +ZPercentFor + ZPercentCut + 
                     ZPercentAqu  + ZPercentBar + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)
anova(M1_IQBP2, MLon_IQBP2)

#Bar
MBar_IQBP2 <- lmer(ZlogitIQBP ~  ZPercentAgg +ZPercentFor + ZPercentCut + 
                     ZPercentAqu  + ZLongitude + (1|IQBPBasins), data=IQBP_FullSummary,REML= FALSE)
anova(M1_IQBP2, MBar_IQBP2)


#SI - ISBs Analyses####
#Section 3.2 - Analyses of subsets of data#### 
#Subset of data with same amount of protected and unprotected watersheds####
#Create a subset of the data with the same amount of points with 0% protection as those with some level of protection
SB_EqualProAndNot <- SB_All %>% 
  group_by(ProtectionCategorical) %>% 
  do(sample_n(.,30))

SB_EqualProAndNot <- as.data.frame(SB_EqualProAndNot)

#Look at data distribution
hist(SB_EqualProAndNot$ISVB)
hist(SB_EqualProAndNot$logitISVB)

#Check correlations
SB_EqualProAndNotCorr<- SB_EqualProAndNot %>% select(ZRelPro, ZRelAgr, ZRelFor, ZRelWet, 
                                                     ZRelCut, ZRelAnthro, ZRelAqu, ZRelBar, 
                                                     ZFlowAccumu,
                                                     ZLatitude,ZLongitude)

cor(SB_EqualProAndNotCorr)
#Nothing above 58% so OK to keep everything in the model

M1SBEqualProAndNot <- lmer(ZlogitISVB ~ ZRelPro + ZRelAgr +ZRelFor+  ZRelWet + 
                             ZRelCut + ZRelAnthro + ZRelAqu + ZRelBar + ZFlowAccumu + 
                             ZLatitude + ZLongitude + (1|SBWBasins), data=SB_EqualProAndNot,REML= FALSE)

#Singular fit, so we proceed with a model that only includes protection
M1SBEqualProAndNot <- lmer(ZlogitISVB ~ ZRelPro + (1|SBWBasins), data=SB_EqualProAndNot,REML= FALSE)


#Validate the model 
#VIF irrelevant 
vif(M1SBEqualProAndNot)


#A) Look at homogeneity
ResidM1SBEqualProAndNot<-resid(M1SBEqualProAndNot)
FitM1SBEqualProAndNot<-fitted(M1SBEqualProAndNot)
plot(x=FitM1SBEqualProAndNot, y= ResidM1SBEqualProAndNot, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)
#It is no surprise that there is a trend in the data given that we have excluded other important variables

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model

#Random factor
boxplot(ResidM1SBEqualProAndNot~SBWBasins, ylab= "Normalized residuals", data = SB_EqualProAndNot, xlab = "Nested Basins")
abline(0,0, lty=2)


#C) Look at normality of residuals
hist(ResidM1SBEqualProAndNot)

#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(M1SBEqualProAndNot)

summary(M1SBEqualProAndNot)



##Chi-Square on Subset of data with equal pro and unpro#
SB_EqualProAndNotTable <- SB_EqualProAndNot %>%
  group_by(ProtectionCategory, WaterQualityCategory) %>%
  summarize(n = n()) %>%
  mutate(freq = n/sum(n))

SBchisqEqualProAndNot<- chisq.test(SB_EqualProAndNot$WaterQualityCategory,SB_EqualProAndNot$ProtectionCategory, correct= FALSE)
SBchisqEqualProAndNot
#Does not make this test significant... 

SB_EqualProAndNotTable2 <- SB_EqualProAndNot %>%
  group_by(ProtectionCategorical, WaterQualityCategory) %>%
  summarize(n = n()) %>%
  mutate(freq = n/sum(n))

SBchisqEqualProAndNot2<- chisq.test(SB_EqualProAndNot$WaterQualityCategory,SB_EqualProAndNot$ProtectionCategorical, correct= FALSE)
SBchisqEqualProAndNot2
#Does not make this test significant either... 



#Subset of data that only includes Protected watesrheds####
hist(SB_ProSubset$ISVB)
hist(SB_ProSubset$logitISVB)
#Check correlations among variables included in the full model 
SB_ProSubetCorr <- SB_ProSubset %>% select(ZRelPro, ZRelAgr, ZRelFor, ZRelWet, 
                                           ZRelCut, ZRelAnthro, ZRelAqu, ZRelBar, 
                                           ZFlowAccumu,
                                           ZLatitude,ZLongitude)

cor(SB_ProSubetCorr)
#% Forest and %Protected are 86% correlated so both cannot be included in the model
#We remove forest from the model 

MPro_SB <- lmer(ZlogitISVB ~ ZRelPro + 
                  ZRelAgr +
                  #ZRelFor+  
                  ZRelWet + 
                  ZRelCut + 
                  ZRelAnthro + 
                  ZRelAqu + 
                  ZRelBar + 
                  ZFlowAccumu + 
                  ZLatitude + 
                  ZLongitude + (1|SBWBasins), data=SB_ProSubset,REML= FALSE)

#The model has a singular fit, which means it is overfitted, so we re-run the model only including protection

MPro_SB_Best <- lmer(ZlogitISVB ~ ZRelPro + 
                       (1|SBWBasins), data=SB_ProSubset,REML= FALSE)


#Validate the model 
#VIF not relevant given that there aren't covariates 

#A) Look at homogeneity
ResidMPro_SB_Best<-resid(MPro_SB_Best)
FitMPro_SB_Best<-fitted(MPro_SB_Best)
plot(x=FitMPro_SB_Best, y= ResidMPro_SB_Best, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model

#Random factor
boxplot(ResidMPro_SB_Best~SBWBasins, ylab= "Normalized residuals", data = SB_ProSubset, xlab = "Nested Basins")
abline(0,0, lty=2)


#C) Look at normality of residuals
hist(ResidMPro_SB_Best)

#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(MPro_SB_Best)

summary(MPro_SB_Best)






#Subsets of watershed sizes####
#Look at data distribution by watershed size 
#Watersheds between 100 and 1000km2 
ISBs_Watersheds100to1000 <- SB_All %>% filter(WatersheadAreakm2 >= 100 & WatersheadAreakm2 <= 1000)
#Reduces sample to 33/145 observations (23%)

#Watersheds below 100km2 
ISBs_WatershedsBelow100<- SB_All %>% filter(WatersheadAreakm2 <= 100)
#Reduces sample to 111/145 observations (77%)

#Create mixed model with the data for watersheds below 100km2#
#Look at data distribution
hist(ISBs_WatershedsBelow100$ZlogitISVB)

#Check for correlations between variables
ISBs_WatershedsBelow100Cor <- ISBs_WatershedsBelow100 %>% select(ZRelPro, ZRelAgr, ZRelFor,  ZRelWet, 
                                                                 ZRelCut, ZRelAnthro, ZRelAqu, ZRelBar, ZFlowAccumu, 
                                                                 ZLatitude, ZLongitude)
cor(ISBs_WatershedsBelow100Cor)
#Nothing higher than 69 so we're OK to proceed

MB100_ISBs <- lmer(ZlogitISVB ~ ZRelPro + ZRelAgr +ZRelFor+  ZRelWet + 
                     ZRelCut + ZRelAnthro + ZRelAqu + ZRelBar + ZFlowAccumu + 
                     ZLatitude + ZLongitude + (1|SBWBasins), data=ISBs_WatershedsBelow100,REML= FALSE)


#Assess the relative importance of the variables included in the model 
options(na.action = "na.fail") 
resultsMB100_ISBs <- dredge(MB100_ISBs)
options(na.action = "na.omit") 
top.modelsMB100_ISBs <- get.models(resultsMB100_ISBs, subset = delta < 5)
top.modelsMB100_ISBs[1]
model.avg(top.modelsMB100_ISBs)
aic_impMB100_ISBs  <- importance(top.modelsMB100_ISBs)
aic_impMB100_ISBs  <- as.data.frame((aic_impMB100_ISBs))
aic_impMB100_ISBs$variables <- rownames(aic_impMB100_ISBs)
aic_impMB100_ISBs$VariableNumbers <- c("K","J","I","H","G","F","E","D","C","B","A")


#Build and Validate best fit model 
MB100_ISBs_Best<- lmer(ZlogitISVB ~ ZRelAgr + 
                         ZRelAnthro + 
                         ZRelWet + 
                         (1|SBWBasins), data=ISBs_WatershedsBelow100,REML= FALSE)


#Validate the model 
vif(MB100_ISBs_Best)
#VIFs are below 5 this is valid 

#A) Look at homogeneity
ResidMB100_ISBs_Best<-resid(MB100_ISBs_Best)
FitMB100_ISBs_Best<-fitted(MB100_ISBs_Best)
plot(x=FitMB100_ISBs_Best, y= ResidMB100_ISBs_Best, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model

#Random factor
boxplot(ResidMB100_ISBs_Best~SBWBasins, ylab= "Normalized residuals", data = ISBs_WatershedsBelow100, xlab = "Nested Basins")
abline(0,0, lty=2)


#C) Look at normality of residuals
hist(ResidMB100_ISBs_Best)

#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(MB100_ISBs_Best)

summary(MB100_ISBs_Best)


#Create mixed model with the data for watersheds below between 100km2 and 1000km2##
#Look at data distribution
hist(ISBs_Watersheds100to1000$ZISVB)

#Check for correlations between variables
ISBs_Watersheds100to1000Cor <- ISBs_Watersheds100to1000 %>% select(ZRelPro, ZRelAgr, ZRelFor,  ZRelWet, 
                                                                   ZRelCut, ZRelAnthro, ZRelAqu, ZRelBar, ZFlowAccumu, 
                                                                   ZLatitude, ZLongitude)
cor(ISBs_Watersheds100to1000Cor)
#Nothing higher than 58% correlated to so its OK to proceed

M100to1000_ISBs <- lmer(ZlogitISVB ~ ZRelPro + ZRelAgr +ZRelFor+  ZRelWet + 
                          ZRelCut + ZRelAnthro + ZRelAqu + ZRelBar + ZFlowAccumu + 
                          ZLatitude + ZLongitude + (1|SBWBasins), data=ISBs_Watersheds100to1000,REML= FALSE)
#The model is singularly fit so we proceed with the simplified model 

MB100_ISBs_Best <- lmer(ZlogitISVB ~ ZRelPro + 
                          (1|SBWBasins), data=ISBs_Watersheds100to1000,REML= FALSE)

#Still has problems with singular fit... 

#validate model 

#A) Look at homogeneity
ResidMB100_ISBs_Best<-resid(MB100_ISBs_Best)
FitMB100_ISBs_Best<-fitted(MB100_ISBs_Best)
plot(x=FitMB100_ISBs_Best, y= ResidMB100_ISBs_Best, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model

#Random factor
boxplot(ResidMB100_ISBs_Best~SBWBasins, ylab= "Normalized residuals", data = ISBs_WatershedsBelow100, xlab = "Nested Basins")
abline(0,0, lty=2)


#C) Look at normality of residuals
hist(ResidMB100_ISBs_Best)

#Determine Conditional and Marginal R2 of the model
r.squaredGLMM(MB100_ISBs_Best)
summary(MB100_ISBs_Best)


#Section 3.3 - Percent Watershed Protection Data Visualization####
#Figure S2####
#Not Weighted 
ISVB_Pro<-ggplot(SB_All, aes(x=PercentPro, y=ISVB, colour = AgriCategory, shape = AgriCategory)) +  
  geom_point(size = 3) +
  scale_color_manual(values=c("#fed98e",  "#fe9929","#cc4c02", "#993404")) +
  geom_hline(yintercept = c(45, 75), colour = c("black", "black"), linetype = "dashed") +
  scale_color_manual(name="% watershed\\nagriculture", 
                     breaks =c("0%", "1-24%", "25-49%","50-74%"),
                     labels = c("0%", "1-24%", "25-49%","50-74%"),
                     values=c("#fed98e",  "#fe9929","#cc4c02", "#993404")) +
  scale_shape_discrete(name="% watershed\\nagriculture", 
                       breaks =c("0%", "1-24%", "25-49%","50-74%", "75-100%"),
                       labels = c("0%", "1-24%", "25-49%","50-74%", "75-100%"))+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=16,face="bold"), 
        panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.title=element_text(size=15, face="bold"), legend.text=element_text(size=14),
        #legend.position = c(0.7, 0.2), 
        legend.background = element_blank(),legend.box.background = element_rect(colour = "white"))+
  annotate("text", x = 110, y = 35, angle = 360, label = "Bad", size = 5)+
  annotate("text", x = 110, y = 60, angle = 360, label = "Precarious         ", size = 5)+
  annotate("text", x = 110, y = 90, angle = 360, label = "Good", size = 5)+
  xlab("% watershed protection") +
  ylab("Water quality \\n(ISBs indicator)")
ISVB_Pro

ggsave(ISVB_Pro, file="FigureS2_ISBsPro.tiff", scale=1, dpi = 600)



#Section 3.4. - Percent Watershed Protection Mixed Model Analysis####
#Create a model that contains all variables 
M1W2 <- lmer(ZlogitISVB ~ ZPercentPro +ZPercentAgg  + ZPercentWet +ZPercentFor + ZPercentCut + ZPercentAnthro+ ZPercentAqu + ZPercentBar+ ZFlowAccumu +ZLatitude +ZLongitude + (1|SBWBasins), data=SB_All,REML= FALSE)
summary(M1W2)
#Use the MuMIn package for model selection and averaging
#OK to model average with mixed model - https://stats.stackexchange.com/questions/351917/model-averaging-in-mixed-models
options(na.action = "na.fail") 
resultsM1W2 <- dredge(M1W2)
options(na.action = "na.omit") 

top.models2 <- get.models(resultsM1W2, subset = delta < 5)
top.models2[1]
model.avg(top.models2)

aic_imp2 <- importance(top.models2)
aic_imp2 <- as.data.frame((aic_imp2))
aic_imp2$variables <- rownames(aic_imp2)
aic_imp2$VariableNumbers <- c("K","J","I","H","G","F","E","D","C","B","A")

# Visualize findings Plotting models Delta AIC < 5
#Figure S4####
ISBsProImportance <-   ggplot(data=aic_imp2, aes(y=aic_imp2[,1],x=VariableNumbers)) +
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() + xlab("Predictor variable") + ylab("Relative importance")  +
  theme(axis.text=element_text(size=12.5), axis.title=element_text(size=15,face="bold"),
        panel.grid.major = (element_blank()), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #Add information about model estimates 
  annotate("text", x = 12.5, y = 1.12, angle = 0, label = "")+
  annotate("text", x = 12, y = 1.08, angle = 0, label = "Average\\nSlope")+
  annotate("text", x = 11.15, y = 1.08, angle = 0, label = "       0.30 + *")+ #wetland
  annotate("text", x = 10.15, y = 1.08, angle = 0, label = "       0.43 + *")+ #Forest
  annotate("text", x = 9.15, y = 1.08, angle = 0, label = "     0.15 +")+ #Longitude
  annotate("text", x = 8.15, y = 1.08, angle = 0, label = "- 0.10")+ #Watershed size
  annotate("text", x = 7.15, y = 1.08, angle = 0, label = "- 0.17") + #Ant
  annotate("text", x = 6.15, y = 1.08, angle = 0, label = "- 0.20") + #Agr
  annotate("text", x = 5.15, y = 1.08, angle = 0, label = "- 0.07")+ #Lat
  annotate("text", x = 4.15, y = 1.08, angle = 0, label = "  0.04")+ #Cut
  annotate("text", x = 3.15, y = 1.08, angle = 0, label = "  0.04") + #Bar
  annotate("text", x = 2.15, y = 1.08, angle = 0, label = "  0.05")+ # Pro
  annotate("text", x = 1.15, y = 1.08, angle = 0, label = "  0.03")+ #Aqu
  scale_x_discrete(labels= c("% aquatic","% protection", "% rocks & lichen",
                             "% planted, cut\\n& burned forests","Latitude",
                             "% agricultural","% anthropogenic",
                             "Watershed size","Longitude", 
                             "% wetlands",
                             "% forest"))
ISBsProImportance
ggsave(ISBsProImportance, file="FigureS4_ISBsWeightedRelImportance.pdf", scale=1, dpi = 900)


#Validate Best fit model
M1W2.2 <- lmer(ZlogitISVB ~ ZPercentFor  + ZPercentWet  +ZLongitude + (1|SBWBasins), data=SB_All,REML= FALSE)

#Assess variance inflation 
vif(M1W2.2)
#Variance inflation (Multicoliniarity)  is too high, this model does not meet its assumptions

par(mfrow=c(1,1)) 
#A) Look at homogeneity
ResidM1W2.2<-resid(M1W2.2)
FM1W2.2<-fitted(M1W2.2)
plot(x=FM1W2.2, y= ResidM1W2.2, xlab = "Fitted values", ylab = "Normalized residuals")
abline(0,0, lty=2)

#B) Look at independence
#Plot residuals vs each covariate in the model and not in the model
#Random factor
boxplot(ResidM1W2.2~SBWBasins, ylab= "Normalized residuals", data = SB_All, xlab = "Nested Basins")
abline(0,0, lty=2)

#In the model
plot(x=SB_All$ZRelWet, y= ResidM1W2.2, xlab = "Percent Wetland", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelFor, y= ResidM1W2.2, xlab = "Percent Forest", ylab = "Normalized resitduals")
abline(0,0, lty=2)
plot(x=SB_All$ZLongitude, y= ResidM1W2.2, xlab = "Percent Forest", ylab = "Normalized resitduals")
abline(0,0, lty=2)

#Not in the model
plot(x=SB_All$ZRelAgr, y= ResidM1W2.2, xlab = "Percent Agricultural", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelAnthro, y= ResidM1W2.2, xlab = "Percent Anthropogenic", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelPro, y= ResidM1W2.2, xlab = "Percent Protection", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZFlowAccumu, y= ResidM1W2.2, xlab = "Flow Accumulation", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelCut, y= ResidM1W2.2, xlab = "Percent Cut", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelBar, y= ResidM1W2.2, xlab = "Percent Barren Land", ylab = "Normalized residuals")
abline(0,0, lty=2)
plot(x=SB_All$ZRelAqu, y= ResidM1W2.2, xlab = "Percent Aquatic", ylab = "Normalized residuals")
abline(0,0, lty=2)

#C) Look at normality of residuals
hist(ResidM1W2.2)

#Model meets assumptions 

#Determine Conditional and Marginal R2 of the model, as well as significance of fixed factors 
r.squaredGLMM(M1W2.2)
summary(M1W2.2)  

#calculate significance of each variable in the best fit model
#significance of each variable
#For
M1W2.2 <- lmer(ZlogitISVB ~ ZPercentFor  + ZPercentWet  +ZLongitude + (1|SBWBasins), data=SB_All,REML= FALSE)
MForW2 <-lmer(ZlogitISVB ~ ZPercentWet  +ZLongitude + (1|SBWBasins), data=SB_All,REML= FALSE)
anova(M1W2.2, MForW2)

#Wetland
MWetW2 <- lmer(ZlogitISVB ~ ZPercentFor  +ZLongitude + (1|SBWBasins), data=SB_All,REML= FALSE)
anova(M1W2.2, MWetW2)

#Longi
MLongW2 <- lmer(ZlogitISVB ~ ZPercentFor  + ZPercentWet  + (1|SBWBasins), data=SB_All,REML= FALSE)
anova(M1W2.2, MLongW2)

#National and provincial parks land use####
ParkLandUse<-read.csv("NatProvParksLandUse.csv")
library(tidyr)

ParkLandUse_Long <- gather(ParkLandUse,Park, ProportionofLandUse, Parc_de_la_Gatineau:Rserve_de_parc_national_du_Canada_de_l_Archipel_de_Mingan,
                           factor_key = TRUE)

TotalLandUse<-ParkLandUse_Long %>% group_by(Park) %>% summarize(TotalLanduse = sum(ProportionofLandUse))

ParkLandUse_Long_Test <- left_join(ParkLandUse_Long, TotalLandUse, by = "Park")

ParkLandUse_Long_Test$PercentLanduse <- (ParkLandUse_Long_Test$ProportionofLandUse/ParkLandUse_Long_Test$TotalLanduse)*100

ParkLandUse_Long_Test$PercentLanduse<-format(round(ParkLandUse_Long_Test$PercentLanduse, 2), nsmall = 2)
ParkLandUse_Long_Test

ParkLandUse_Long_Test$PercentLanduse <- as.numeric(ParkLandUse_Long_Test$PercentLanduse)
str(ParkLandUse_Long_Test)
ParkLandUseSummary <- ParkLandUse_Long_Test %>% group_by(LABEL) %>% summarise(max = max(PercentLanduse), min = min(PercentLanduse), mean = mean(PercentLanduse))


