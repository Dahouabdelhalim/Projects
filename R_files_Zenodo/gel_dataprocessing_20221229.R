## Marcus Griffiths, Alex Liu (2022)
## Data processing for 3D gel imager

#################
## User setup  ## ------------------------------------------------------------------------------------------------------
#################
## 1) Set working directory --------------------------------------------------------------------------------------------
getwd()
setwd("T:/FOR/for_AlexL/3DGel/TaGAA/10-26-22 sandbox/TaGAA_N_Accessions") #PC
setwd("C:\\\\Users\\\\USERNAME\\\\FILEPATH") #PC
setwd("/Users/USERNAME/FILEPATH") #macOS

# 2) Plantdata .csv files share the same PlantID column and identifier barcode eg. 011119_NAM1_B73_p1_t100_r1
# 3) Install or load following packages
library(tidyverse)      #includes packages ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats
library(ggcorrplot)     #for correlation matrix plots
library(lmerTest)       #lme4, Matrix, lmer, step # for heritability & stat analysis
library(broom)          #for summarizing statistical model objects in tidy tibbles
library(agricolae)      #for tukeyhsd plots
library(ggpubr)         #for regression equations and r squared values
library(factoextra)     #for pca plots
library(FactoMineR)     #for pca plots

##enter experiment parameter
exp <- "GAA"
species <- "Ta"
expname <- paste0(species,exp)       #experiment name for filenames

## Import temporal plant measure .csv files
RootGIA3Ddata <- read_csv(file=paste(expname,"_dataprocessing/",expname,"_GIA3D_data.csv", sep=""), na = c("NA", "na", "n.a.", ""))
RootDynamicdata <- read_csv(file=paste(expname,"_dataprocessing/",expname,"_DynamicRoots_data.csv", sep=""), na = c("NA", "na", "n.a.", ""))
## Import plant harvest measures .csv files
Biomassdata <- read_csv(file=paste(expname,"_dataprocessing/",expname,"_Biomass_data.csv", sep=""), na = c("NA", "na", "n.a.", ""))
ShootImagedata <- read_csv(file=paste(expname,"_dataprocessing/",expname,"_ShootImage_data.csv", sep=""), na = c("NA", "na", "n.a.", ""))
ExtraRootdata <- read_csv(file=paste(expname,"_dataprocessing/",expname,"_ExtraRoot_data.csv", sep=""), na = c("NA", "na", "n.a.", ""))

  
## merge temporal measures
Plantdata <- RootGIA3Ddata %>%
  full_join(RootDynamicdata, by = "PlantID") %>% full_join(ExtraRootdata, by = "PlantID")


Plantdata <-  Plantdata %>%
  separate(PlantID,
           into = c("PlantID", "Day"),
           sep = '_d')

Plantdata <- Plantdata %>%
  full_join(Biomassdata, by = "PlantID") %>%
  full_join(ShootImagedata, by = "PlantID")

rm(RootGIA3Ddata
   ,RootDynamicdata
   ,Biomassdata
   ,ShootImagedata
   ,ExtraRootdata
   )

##########################
## Extract PlantID info ## ---------------------------------------------------------------------------------------------
##########################
name <- c("Plant1"
          ,"Treatment1"
          ,"Rep1")

Plantdata <-  Plantdata %>%
  separate(PlantID,
           into = c("Date", "Experiment", "Geno", "Plant", "Treatment", "Rep"),
           # into = c("Date", "Experiment", "Geno", "Plant", "Rep"),
           sep = '_') %>%
  separate(`Plant`,
           into = c("Plant1","Plant"),
           sep = 'p') %>%
 separate(`Treatment`,
          into = c("Treatment1","Treatment"),
          sep = 't') %>%
  separate(`Rep`,
           into = c("Rep1","Rep"),
           sep = 'r') %>%
  select(-one_of(name))

rm(name)

###################
## offset ## 
## if sample day goes to 21, those plants should be timeshifted back 4 days back to d17
###################

# concept: 
# if Day == 21:
# Day of every entry with plant # x = Day-4
# get list of rows to target
# change the day to numeric
Plantdata <- Plantdata %>% group_by(Plant) %>% mutate(Day = as.numeric(Day))

shifted <- which(Plantdata$Day == 21)
late_plants <- Plantdata$Plant[shifted]
#Plantdata$Plant %in% late_plants
Plantdata <- Plantdata %>% 
  mutate(Day = replace(Day, Plant %in% late_plants, Day-4)) %>% 
  ungroup()

###################
## Sample counts ## ----------------------------------------------------------------------------------------------------
###################
SampleCount <- Plantdata %>% 
  group_by(Experiment
           ,Geno
           ,Treatment
           ,Day
           ) %>% 
  count(Geno)

dir.create(paste0(expname,"_dataprocessing/output"),showWarnings = TRUE)
write_csv(SampleCount, paste0(expname,"_dataprocessing/output/",expname,"_samplecount.csv"))

####################
## Trait renaming ## ---------------------------------------------------------------------------------------------------
####################
Plantdata <- Plantdata %>% 
  rename(RootConvexHullVolume3D_mm3 = "ConvexHullVolume3D(mm^3)") %>%
  rename(RootSolidity3D = Solidity3D) %>%
  rename(RootBushiness3D = Bushiness3D) %>%
  rename(RootDepth3D_mm = "Depth3D(mm)") %>%
  rename(RootMaximumNetworkWidth3D_mm = "MaximumNetworkWidth3D(mm)") %>%
  rename(RootLengthDistribution3D_mm = "LengthDistribution3D(mm)") %>%
  rename(RootWidthDepth3D_Ratio = "WidthDepthRatio3D") %>%
  rename(RootTotalCount = "TotalRootNumber") %>%
  rename(RootTotalVolume_mm3 = "TotalRootVolume(mm^3)") %>%
  rename(RootTotalLength_mm = "TotalRootLength(mm)") %>%
  rename(RootTotalDepth_mm = "TotalRootDepth(mm)") %>%
  rename(RootMeanVolume_mm3 = "MeanRootVolume(mm^3)") %>%
  rename(RootMeanLength_mm = "MeanRootLength(mm)") %>%
  rename(RootMeanDepth_mm = "MeanRootDepth(mm)") %>%
  rename(RootMeanTortuosity = MeanRootTortuosity) %>%
  rename(RootMeanRadius_mm = "MeanRootRadius(mm)") %>%
  rename(RootMeanSoilAngle_degrees = MeanRootSoilAngle) %>%
  rename(RootMeanBranchingAngle_degrees = MeanRootBranchingAngle) %>%
  rename(RootMedianVolume_mm3 = "MedianRootVolume(mm^3)") %>%
  rename(RootMedianLength_mm = "MedianRootLength(mm)") %>%
  rename(RootMedianDepth_mm = "MedianRootDepth(mm)") %>%
  rename(RootMedianTortuosity = MedianRootTortuosity) %>%
  rename(RootMedianRadius_mm = "MedianRootRadius(mm)") %>%
  rename(RootMedianSoilAngle_degrees = MedianRootSoilAngle) %>%
  rename(RootMedianBranchingAngle_degrees = MedianRootBranchingAngle) %>%
  rename(RootCountLateral = TotalLateralRootNumber) %>%
  rename(RootVolumeLateral_mm3 = "TotalLateralRootVolume(mm^3)") %>%
  rename(RootLengthLateral_mm = "TotalLateralRootLength(mm)") %>%
  rename(RootTotalDepthLateral_mm = "TotalLateralRootDepth(mm)") %>%
  rename(RootMeanVolumeLateral_mm3 = "MeanLateralRootVolume(mm^3)") %>%
  rename(RootMeanLengthLateral_mm = "MeanLateralRootLength(mm)") %>%
  rename(RootMeanDepthLateral_mm = "MeanLateralRootDepth(mm)") %>%
  rename(RootMeanTortuosityLateral = MeanLateralRootTortuosity) %>%
  rename(RootMeanRadiusLateral_mm = "MeanLateralRootRadius(mm)") %>%
  rename(RootMeanSoilAngleLateral_degrees = MeanLateralRootSoilAngle) %>%
  rename(RootMeanBranchingAngleLateral_degrees = MeanLateralRootBranchingAngle) %>%
  rename(RootMedianVolumeLateral_mm3 = "MedianLateralRootVolume(mm^3)") %>%
  rename(RootMedianLengthLateral_mm = "MedianLateralRootLength(mm)") %>%
  rename(RootMedianDepthLateral_mm = "MedianLateralRootDepth(mm)") %>%
  rename(RootMedianTortuosityLateral = MedianLateralRootTortuosity) %>%
  rename(RootMedianRadiusLateral_mm = "MedianLateralRootRadius(mm)") %>%
  rename(RootMedianSoilAngleLateral_degrees = MedianLateralRootSoilAngle) %>%
  rename(RootMedianBranchingAngleLateral_degrees = MedianLateralRootBranchingAngle) %>%
  rename(RootCountFirstOrderLateral = TotalFirstOrderLateralRootNumber) %>%
  rename(RootVolumeFirstOrderLateral_mm3 = "TotalFirstOrderLateralRootVolume(mm^3)") %>%
  rename(RootLengthFirstOrderLateral_mm = "TotalFirstOrderLateralRootLength(mm)") %>%
  rename(RootTotalDepthFirstOrderLateral_mm = "TotalFirstOrderLateralRootDepth(mm)") %>%
  rename(RootMeanVolumeFirstOrderLateral_mm3 = "MeanFirstOrderLateralRootVolume(mm^3)") %>%
  rename(RootMeanLengthFirstOrderLateral_mm = "MeanFirstOrderLateralRootLength(mm)") %>%
  rename(RootMeanDepthFirstOrderLateral_mm = "MeanFirstOrderLateralRootDepth(mm)") %>%
  rename(RootMeanTortuosityFirstOrderLateral = MeanFirstOrderLateralRootTortuosity) %>%
  rename(RootMeanRadiusFirstOrderLateral_mm = "MeanFirstOrderLateralRootRadius(mm)") %>%
  rename(RootMeanSoilAngleFirstOrderLateral_degrees = MeanFirstOrderLateralRootSoilAngle) %>%
  rename(RootMeanBranchingAngleFirstOrderLateral_degrees = MeanFirstOrderLateralRootBranchingAngle) %>%
  rename(RootMedianVolumeFirstOrderLateral_mm = "MedianFirstOrderLateralRootVolume(mm^3)") %>%
  rename(RootMedianLengthFirstOrderLateral_mm = "MedianFirstOrderLateralRootLength(mm)") %>%
  rename(RootMedianDepthFirstOrderLateral_mm = "MedianFirstOrderLateralRootDepth(mm)") %>%
  rename(RootMedianTortuosityFirstOrderLateral = "MedianFirstOrderLateralRootTortuosity") %>%
  rename(RootMedianRadiusFirstOrderLateral_mm = "MedianFirstOrderLateralRootRadius(mm)") %>%
  rename(RootMedianSoilAngleFirstOrderLateral_degrees = MedianFirstOrderLateralRootSoilAngle) %>%
  rename(RootMedianBranchingAngleFirstOrderLateral_degrees = MedianFirstOrderLateralRootBranchingAngle) %>%
  rename(RootDensityFirstOrderLateral_TL = DensityFirstOrderLateralRoot_TL) %>%
  rename(RootDensityFirstOrderLateral_BRTL = DensityFirstOrderLateralRoot_BRTL) %>%
  rename(InterbranchMeanDistance_mm = "MeanInterbranchDistance(mm)") %>%
  rename(InterbranchMedianDistance_mm = "MedianInterbranchDistance(mm)") %>%
  rename(RootVolumePrimary_mm3 = "PrimaryRootVolume(mm^3)") %>%
  rename(RootLengthPrimary_mm = "PrimaryRootLength(mm)") %>%
  rename(RootDepthPrimary_mm = "PrimaryRootDepth(mm)") %>%
  rename(RootTortuosityPrimary = "PrimaryRootTortuosity") %>%
  rename(RootRadiusPrimary = "PrimaryRootRadius(mm)") %>%
  rename(RootSoilAnglePrimary_degrees = "PrimaryRootSoilAngle")

########################
## Trait calculations ## -----------------------------------------------------------------------------------------------
########################
Plantdata$Day <- as.numeric(Plantdata$Day)

## arrange & calculate timeperiod between scans
Plantdata <- Plantdata  %>%
  arrange(Date
          ,Experiment
          ,Geno
          ,Plant
          ,Treatment
          ,Rep
          ,Day) %>%
  mutate(Timeperiod = (Day-lag(Day))) ##order of rows is important for lag function

Plantdata$Timeperiod[Plantdata$Timeperiod < 0 | is.na(Plantdata$Timeperiod)] <- 5


# replace NAs in Extended_Root with 0s
Plantdata$Extended_Root[is.na(Plantdata$Extended_Root)] <- 0


## reorder columns
Plantdata <- Plantdata %>% select(Experiment
                                  ,Date
                                  ,Geno
                                  ,Plant
                                  ,Treatment
                                  ,Rep
                                  ,Day
                                  ,Timeperiod
                                  ,RootDW_g
                                  ,ShootDW_g
                                  ,everything()
                                  )
# adding manually measured primary root lengths beyond camera FOV

# root traits
Plantdata <- Plantdata %>% ##add comments for traits
  mutate(RootLengthSecondOrderLateral_mm = RootLengthLateral_mm-RootLengthFirstOrderLateral_mm
         ,RootTotalGrowthRate_mm.h = (RootTotalLength_mm-lag(RootTotalLength_mm))/(Timeperiod*24) 
         ,RootLateralGrowthRate_mm.h = (RootLengthLateral_mm-lag(RootLengthLateral_mm))/(Timeperiod*24)
         ,SpecRootLength_m.g = ((RootTotalLength_mm+Extended_Root)/1000)/(RootDW_g)
         ## get correct units
         ,SpecLeafArea_g.m2 = (ShootDW_g)/(PlantCVShootArea_unit2/10000)
         ,TotalMass_g = (RootDW_g+ShootDW_g)
         ,RootMassFract_g.g = (RootDW_g/TotalMass_g)
         )

# get the negatives of lag
# there has to be a better way of doing this
Plantdata$RootTotalGrowthRate_mm.h[Plantdata$RootTotalGrowthRate_mm.h < 0 |is.na(Plantdata$RootTotalGrowthRate_mm.h)] <- 
  Plantdata$RootTotalLength_mm[Plantdata$RootTotalGrowthRate_mm.h < 0 |is.na(Plantdata$RootTotalGrowthRate_mm.h)]/
  Plantdata$Timeperiod[Plantdata$RootTotalGrowthRate_mm.h < 0 |is.na(Plantdata$RootTotalGrowthRate_mm.h)]/
  24
Plantdata$RootLateralGrowthRate_mm.h[Plantdata$RootLateralGrowthRate_mm.h < 0 | is.na(Plantdata$RootLateralGrowthRate_mm.h)] <- 
  Plantdata$RootLengthLateral_mm[Plantdata$RootLateralGrowthRate_mm.h < 0 | is.na(Plantdata$RootLateralGrowthRate_mm.h)]/
  Plantdata$Timeperiod[Plantdata$RootLateralGrowthRate_mm.h < 0 | is.na(Plantdata$RootLateralGrowthRate_mm.h)]/
  24

## save it
write_csv(Plantdata, paste0(expname,"_dataprocessing/output/",expname,"_plantdata_unfiltered.csv"))

###################
## Data curation ## ----------------------------------------------------------------------------------------------------
###################

any(is.na(Plantdata))

  ### Alex only throw out 1 reps per day rather than whole genotype?

Plantdata <- Plantdata %>% replace(is.na(.), 0)
Plantdata <- Plantdata %>% mutate_all(function(x) ifelse(is.infinite(x), 0, x))

temporaldata <- Plantdata %>% select(-Experiment
                                     ,-Date
                                     ,-Plant
                                     ,-Timeperiod
                                     ,-PlantCVShootArea_unit2
                                     ,-PlantCVLeafCount
                                     ,-SpecLeafArea_g.m2
                                     ,-RootDW_g
                                     ,-ShootDW_g
                                     ,-SpecRootLength_m.g
                                     ,-TotalMass_g
                                     ,-RootMassFract_g.g
                                     #,-RootDensityFirstOrderLateral_BRTL #inf value
                                     )

## remove plants with only one rep otherwise you'll break the ANOVA
# one_reps <- SampleCount[which(SampleCount$n ==1), 2:3]

# remove plants with any of the timepoints with one rep otherwise you'll break the ANOVA
plants_rm <- SampleCount$Geno[which(SampleCount$n ==1)]
one_reps <- SampleCount[SampleCount$Geno %in% plants_rm, 2:3]

# remove plants that are missing a day 
temporaldata <- anti_join(temporaldata, one_reps)
Plantdata <- anti_join(Plantdata, one_reps)

write_csv(Plantdata, paste0(expname,"_dataprocessing/output/",expname,"_plantdata.csv"))

####################
## Means & Errors ## ---------------------------------------------------------------------------------------------------
####################
dat <- Plantdata
dat$Timeperiod <- dat$Timeperiod %>% as_factor()
dat <- dat %>% mutate(across(where(is.numeric), ~na_if(., Inf)), across(where(is.numeric), ~na_if(., -Inf)))

## define standard error function --------------------------------------------------------------------------------------
stderror <- function(x) sd(na.omit(x))/sqrt(length(na.omit(x)))

## calculate means & errors --------------------------------------------------------------------------------------------
dat_mean <- dat %>% 
  group_by(Geno
           ,Treatment
           ,Day
           ) %>% 
  summarise(
    across(
      .cols = where(is.numeric)
      ,mean
      ,.names = "mean_{.col}"
      ,na.rm = TRUE
    )
  ) %>% 
  ungroup()

dat_se <- dat %>% 
  group_by(Geno
           ,Treatment
           ,Day
           ) %>% 
  summarise(
    across(
      .cols = where(is.numeric)
      ,stderror
      ,.names = "se_{.col}"
    )
  ) %>% 
  ungroup()

dat_meanErrors <- merge(dat_mean, dat_se, by=c("Geno"
                                               ,"Treatment"
                                               ,"Day"
                                               ))

write_csv(dat_meanErrors, paste0(expname,"_dataprocessing/output/",expname,"_plantdata_MeansErrors.csv"))

rm(dat
   ,dat_mean
   ,dat_se
   )

########################
## correlation matrix ## -----------------------------------------------------------------------------------------------
########################
dir.create(paste0(expname,"_dataprocessing/output/corrplot"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/corrplot/png"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/corrplot/pdf"),showWarnings = TRUE)

######################
## correlation matrix all data
## 

## keep variable columns only
## Remove Inf value traits if present

# dat <- Plantdata %>% 
#   filter(Day==17) %>%
#   # filter(Treatment=="HighN") %>%
#   select(-Experiment
#         ,-Date
#         ,-Geno
#         ,-Plant
#         ,-Treatment
#         ,-Day
#         ,-Timeperiod
#         ,-Rep
#         )

dat <- Plantdata %>%
  filter(Day==17) %>%
  filter(Treatment=="HighN") %>%
  #keep variable columns only
  select(RootTotalLength_mm
         ,RootTotalVolume_mm3
         ,RootTotalCount
         ,RootDepth3D_mm
         ,RootBushiness3D
         ,RootMaximumNetworkWidth3D_mm
         ,RootConvexHullVolume3D_mm3
         # ,RootWidthDepth3D_Ratio
         ,RootLengthDistribution3D_mm
         ,RootSolidity3D
         # ,RootMeanRadius_mm
         ,SpecRootLength_m.g
         # ,RootTotalGrowthRate_mm.h
         ,InterbranchMeanDistance_mm
         ,RootMedianBranchingAngleLateral_degrees
  )

dat <- dat %>% mutate_all(function(x) ifelse(is.infinite(x), 0, x))

corr <- cor(dat, use = "pairwise.complete.obs")
pmat <- cor_pmat(dat) #compute matrix of correlation p-values

ggcorrplot(corr
           ,method = "circle"
           ,hc.order = TRUE
           #,type = "lower"
           #,lab = TRUE
           ,p.mat = pmat
           ,outline.col = "white"
           ) +
  theme(legend.position="top")
ggsave(paste0(expname,"_dataprocessing/output/corrplot/png/",expname,"_corrplot_plantdata_order_subset.png"), width=8, height=8, dpi=300)
ggsave(paste0(expname,"_dataprocessing/output/corrplot/pdf/",expname,"_corrplot_plantdata_order_subset.pdf"), width=8, height=8, useDingbats=FALSE)
# ggsave(paste0(expname,"_dataprocessing/output/corrplot/png/",expname,"_corrplot_plantdata_order.png"), width=28, height=28, dpi=300)
# ggsave(paste0(expname,"_dataprocessing/output/corrplot/pdf/",expname,"_corrplot_plantdata_order.pdf"), width=28, height=28, useDingbats=FALSE)

corr %>% as_tibble() %>% write_csv(paste0(expname,"_dataprocessing/output/corrplot/",expname,"_corr_plantdata_subset.csv")) 
pmat %>% as_tibble() %>% write_csv(paste0(expname,"_dataprocessing/output/corrplot/",expname,"_pmat_plantdata_subset.csv"))

rm(dat
   ,corr
   ,pmat
   )

##################
## PCA analysis ## -----------------------------------------------------------------------------------------------------
##################
dir.create(paste0(expname,"_dataprocessing/output/pca"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/pca/png"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/pca/pdf"),showWarnings = TRUE)

## keep variable columns only

# dat <- Plantdata %>%
#   filter(Day==17) %>%
#   # filter(Treatment=="HighN") %>%
#   select(-Experiment
#          ,-Date
#          ,-Geno
#          ,-Plant
#          ,-Treatment
#          ,-Day
#          ,-Timeperiod
#          ,-Rep
#   )

dat <- Plantdata %>%
  filter(Day==17) %>%
  filter(Treatment=="HighN") %>%
  #keep variable columns only
  select(RootTotalLength_mm
         ,RootTotalVolume_mm3
         ,RootTotalCount
         ,RootDepth3D_mm
         ,RootBushiness3D
         ,RootMaximumNetworkWidth3D_mm
         ,RootConvexHullVolume3D_mm3
         ,RootWidthDepth3D_Ratio
         ,RootLengthDistribution3D_mm
         ,RootSolidity3D
         ,RootMeanRadius_mm
         ,SpecRootLength_m.g
         # ,RootTotalGrowthRate_mm.h
         ,InterbranchMeanDistance_mm
         ,RootMedianBranchingAngleLateral_degrees
  )

## Remove NA and Inf value traits if present
dat <- dat %>% replace(is.na(.), 0)
dat <- dat %>% mutate_all(function(x) ifelse(is.infinite(x), 0, x))

res.pca <- PCA(dat, graph = FALSE)
get_eig(res.pca)
eig_output <- get_eig(res.pca)
write.csv(eig_output, paste0(expname,"_dataprocessing/output/pca/PCA_eig_plantdata_subset.csv"))

fviz_screeplot(res.pca
               ,addlabels = TRUE
               ,ggtheme = theme_bw()
               ,title = ""
               )
ggsave(paste0(expname,"_dataprocessing/output/pca/png/",expname,"_elbow_plantdata_subset.png"), width=7, height=5, dpi=300)
ggsave(paste0(expname,"_dataprocessing/output/pca/pdf/",expname,"_elbow_plantdata_subset.pdf"), width=7, height=5)

var <- get_pca_var(res.pca)
var_contrib <- as.data.frame(var$contrib)
write.csv(var_contrib, paste0(expname,"_dataprocessing/output/pca/PCA_contribution_plantdata_subset.csv"))

fviz_pca_var(res.pca
             ,col.var="contrib"
             ,gradient.cols = c("blue", "red")
             ,ggtheme = theme_bw()
             ,repel = FALSE # Avoid text overlapping
             ,title = ""
             )
ggsave(paste0(expname,"_dataprocessing/output/pca/png/",expname,"_PCA_plantdata_subset.png"), width=7, height=5, dpi=300)
ggsave(paste0(expname,"_dataprocessing/output/pca/pdf/",expname,"_PCA_plantdata_subset.pdf"), width=7, height=5)

fviz_contrib(res.pca
             ,choice = "var"
             ,axes = 1
             ,top = 10
             ,title = ""
             )
ggsave(paste0(expname,"_dataprocessing/output/pca/png/",expname,"_contrib_plantdata_subset_PC1.png"), width=7, height=5, dpi=300)
ggsave(paste0(expname,"_dataprocessing/output/pca/pdf/",expname,"_contrib_plantdata_subset_PC1.pdf"), width=7, height=5)

rm(dat
   ,res.pca
   ,eig_output
   ,var
   ,var_contrib
   )

###################
## Stat analysis ## ----------------------------------------------------------------------------------------------------
###################
dir.create(paste0(expname,"_dataprocessing/output/stat_analysis"),showWarnings = TRUE)

## define significance stars function
make_significance_stars <- function(pval) {    #make star function
  stars = ""
  if(pval <= 0.001)
    stars = "***"
  if(pval > 0.001 & pval <= 0.01)
    stars = "**"
  if(pval > 0.01 & pval <= 0.05)
    stars = "*"
  if(pval > 0.05)
    stars = "ns"
  stars
}

## define Day * Geno + (1|Rep) anova function
compute_aov_DayGeno <- function(x)
{
  dat <- x
  dat[[1]] <- factor(dat[[1]])
  dat[[2]] <- factor(dat[[2]])
  dat[[3]] <- factor(dat[[3]])
  
  resultOutput <- data.frame()
  
  for (i in 4:length(colnames(dat)))
  {
    print(colnames(dat)[i])
    # if(i != 31 && i != 32 && i != 36 && i != 37) #drop variables which are not present in both treatments
    {
      anova_output <- anova(lmer(dat[[i]] ~ Day * Geno + (1|Rep), data=dat))
      anova_output <- broom::tidy(anova_output)
      anova_output <- anova_output %>%
        mutate(trait = colnames(dat)[i]) %>%
        mutate(signif = sapply(p.value, function(x) make_significance_stars(x))) %>%
        mutate(statistic_2dp = formatC(statistic, digits=2,format="f")) %>%
        mutate(fvalue_significance = paste(statistic_2dp, signif))
      resultOutput <- bind_rows(resultOutput, anova_output)
    }
  }
  
  return(resultOutput)
}

## Day * Genotype
dir.create(paste0(expname,"_dataprocessing/output/stat_analysis/DayxGeno"),showWarnings = TRUE)

dat <- temporaldata %>% select(everything()
                               ,-Treatment
                               )

dat <- dat %>% replace(is.na(.), 0)
dat <- dat %>% mutate_all(function(x) ifelse(is.infinite(x), 0, x))

resultOutput <- compute_aov_DayGeno(dat)
write_csv(resultOutput, paste0(expname,"_dataprocessing/output/stat_analysis/DayxGeno/",expname,"_aov_dat_DayxGeno_raw.csv"))
resultOutputSummaryTable <- resultOutput %>% select(term, fvalue_significance, trait) %>%
  spread(term, fvalue_significance)
write_csv(resultOutputSummaryTable, paste0(expname,"_dataprocessing/output/stat_analysis/DayxGeno/",expname,"_aov_dat_DayxGeno_summary.csv"))

rm(resultOutput
   ,resultOutputSummaryTable
   ,dat
  )

########################
## TukeyHSD (posthoc) ## -----------------------------------------------------------------------------------------------
########################
dir.create(paste0(expname,"_dataprocessing/output/stat_analysis_posthoc/"),showWarnings = TRUE)

## define Day * Geno + Rep tukeyHSD function
compute_aov_tukey_DayGeno <- function(x)
{
  dat <- x
  dat[[1]] <- factor(dat[[1]])
  dat[[2]] <- factor(dat[[2]])
  # rep needs to be factored too
  dat[[3]] <- factor(dat[[3]])
  
  for (i in 4:length(colnames(dat)))
  {
    print(colnames(dat)[i])
    # if(i != 31 && i != 32 && i != 44 && i != 45) #drop variables which are not present in both treatments
    {
      print(i)
      summarytable <- tidy(aov(dat[[i]] ~ Day * Geno + Rep, data = dat, na.action=na.omit))
      write_csv(summarytable, paste0(expname,"_dataprocessing/output/stat_analysis_posthoc/DayxGeno/",colnames(dat)[i],"_aov.csv"))
      paste('made it')
      tempdf <- dat[c(1,2,3,i)];
      thsd <- tidy(TukeyHSD(aov(tempdf[[4]] ~ Day * Geno + Rep, data=tempdf, na.action=na.omit)))
      write_csv(thsd, paste0(expname,"_dataprocessing/output/stat_analysis_posthoc/DayxGeno/TukeyHSD/",colnames(dat)[i],"_TukeyHSD.csv"))
    }
  }
}

## Day * Genotype
dat <- temporaldata %>%
  select(everything()
         ,-Treatment
  )

dir.create(paste0(expname,"_dataprocessing/output/stat_analysis_posthoc/DayxGeno/"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/stat_analysis_posthoc/DayxGeno/TukeyHSD/"),showWarnings = TRUE)

compute_aov_tukey_DayGeno(dat)

rm(dat)

#########################
## Interaction boxplots## ----------------------------------------------------------------------------------------------
#########################
# read in Plantdata if just starting here
Plantdata <- read_csv(paste0(expname,"_dataprocessing/output/",expname,"_Plantdata.csv"))

## Day * Geno
boxplot_fun_DayGeno = function(y) {
  ggplot(dat, aes(x = interaction(dat[factors]), y=.data[[y]])) +
    # ggplot(dat, aes(fct_reorder(Group, .data[[y]]), y=.data[[y]])) +
    geom_boxplot(aes(
      color=Geno
    )) +
    ylab(y) +
    theme_bw() +
    theme(
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,axis.text.x = element_text(angle=45,hjust=1)
    )
  ggsave(paste0(expname,"_dataprocessing/output/plant_boxplots_all/DayxGeno/png/",expname,"_",y,"_HighN.png"), width=6, height=3.1968, dpi=300)
  ggsave(paste0(expname,"_dataprocessing/output/plant_boxplots_all/DayxGeno/pdf/",expname,"_",y,"_HighN.pdf"), width=6, height=3.1968, useDingbats=FALSE)
}

dat <- Plantdata %>% select(Geno,Day,everything())
dat <- Plantdata %>% filter(Treatment=="HighN") %>% select(Geno,Day,everything())
# dat <- Plantdata %>% filter(Treatment=="TraceN") %>% select(Day,Geno,everything())

dir.create(paste0(expname,"_dataprocessing/output/plant_boxplots_all/"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/plant_boxplots_all/DayxGeno/"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/plant_boxplots_all/DayxGeno/pdf"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/plant_boxplots_all/DayxGeno/png"),showWarnings = TRUE)

colnames(dat)

factors = names(dat[1:2])
factors = set_names(factors)
factors

traits = names(dat[9:length(dat)])
traits = set_names(traits)
traits

dat <- unite(dat
             ,all_of(factors)
             ,col = "Group"
             ,sep = "_"
             ,remove = FALSE
)

all_crossed <- expand.grid(traits = traits,
                           stringsAsFactors = FALSE
)

map(all_crossed$traits, boxplot_fun_DayGeno)

rm(dat
   ,factors
   ,traits
   ,all_crossed
)

## Geno * Treatment
d6data <- Plantdata %>% filter(Day=="6")
d9data <- Plantdata %>% filter(Day=="9")
d13data <- Plantdata %>% filter(Day=="13")
d17data <- Plantdata %>% filter(Day=="17")

## define Geno * Treatment boxplot function
boxplot_fun_GenoTreatment = function(y) {
  # ggplot(dat, aes(x = interaction(dat[factors]), y=.data[[y]])) +
    ggplot(dat, aes(fct_reorder(Group, .data[[y]]), y=.data[[y]])) +
    geom_boxplot(aes(color=Treatment)) +
    ylab(y) +
    theme_bw() +
    theme(
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,axis.text.x = element_text(angle=45,hjust=1)
    )
  ggsave(paste0(expname,"_dataprocessing/output/plant_boxplots_all/GenoxTreatment/d17/png/",expname,"_",y,".png"), width=6, height=3.1968, dpi=300)
  ggsave(paste0(expname,"_dataprocessing/output/plant_boxplots_all/GenoxTreatment/d17/pdf/",expname,"_",y,".pdf"), width=6, height=3.1968, useDingbats=FALSE)
}

dir.create(paste0(expname,"_dataprocessing/output/plant_boxplots_all/"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/plant_boxplots_all/GenoxTreatment/"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/plant_boxplots_all/GenoxTreatment/d17/"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/plant_boxplots_all/GenoxTreatment/d17/pdf"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/plant_boxplots_all/GenoxTreatment/d17/png"),showWarnings = TRUE)

dat <- d17data %>% select(Geno
                          ,Treatment
                          ,everything())
colnames(dat)

factors = names(dat[1:2])
factors = set_names(factors)
factors

traits = names(dat[9:length(dat)])
traits = set_names(traits)
traits

dat <- unite(dat
             ,all_of(factors)
             ,col = "Group"
             ,sep = "_"
             ,remove = FALSE
             )

all_crossed <- expand.grid(traits = traits,
                           stringsAsFactors = FALSE
                           )

map(all_crossed$traits, boxplot_fun_GenoTreatment)

rm(dat
   ,d6data
   ,d9data
   ,d13data
   ,d17data
   ,traits
   ,factors
   ,all_crossed
   )

## Day * Treatment
dir.create(paste0(expname,"_dataprocessing/output/plant_boxplots_all/"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/plant_boxplots_all/DayxTreatment/"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/plant_boxplots_all/DayxTreatment/pdf"),showWarnings = TRUE)
dir.create(paste0(expname,"_dataprocessing/output/plant_boxplots_all/DayxTreatment/png"),showWarnings = TRUE)

boxplot_fun_DayTreatment = function(y) {
  ggplot(dat, aes(x = interaction(dat[factors]), y=.data[[y]])) +
  # ggplot(dat, aes(fct_reorder(Group, .data[[y]]), y=.data[[y]])) +
    geom_boxplot(aes(
      color=Treatment
    )) +
    ylab(y) +
    theme_bw() +
    theme(
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,axis.text.x = element_text(angle=45,hjust=1)
    )
  ggsave(paste0(expname,"_dataprocessing/output/plant_boxplots_all/DayxTreatment/png/",expname,"_",y,".png"), width=6, height=3.1968, dpi=300)
  ggsave(paste0(expname,"_dataprocessing/output/plant_boxplots_all/DayxTreatment/pdf/",expname,"_",y,"_.pdf"), width=6, height=3.1968, useDingbats=FALSE)
}

dat <- Plantdata %>% select(Day
                            ,Treatment
                            ,everything()
                            )
colnames(dat)

factors = names(dat[1:2])
factors = set_names(factors)
factors

traits = names(dat[9:length(dat)])
traits = set_names(traits)
traits

dat <- unite(dat
             ,all_of(factors)
             ,col = "Group"
             ,sep = "_"
             ,remove = FALSE
             )

all_crossed <- expand.grid(traits = traits,
                           stringsAsFactors = FALSE
                           )

map(all_crossed$traits, boxplot_fun_DayTreatment)

rm(dat
   ,factors
   ,traits
   ,all_crossed
   )

















