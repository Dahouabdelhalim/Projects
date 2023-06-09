##############################################################################################################
### Supplementary Information for ###

# Title: Species out of sight: elucidating the determinants of research effort in global reptiles
# Authors: Jhonny J M Guedes1, Mario R Moura2, José A. F. Diniz-Filho3
# 1 Programa de Pós-Graduação em Ecologia e Evolução, Departamento de Ecologia, 
# Universidade Federal de Goiás (UFG) – Campus Samambaia, 74690-900, Goiânia, GO, Brazil
# 2 Instituto de Ciências Biológicas, Dept. de Ecologia, UFG – Campus Samambaia, 74690-900, Goiânia, GO, Brazil
# 3 Departamento de Ciências Biológicas, Universidade Federal da Paraíba, Areia, Paraiba, Brazil
# Corresponding author: Jhonny J M Guedes; E-mai: jhonnyguds@gmail.com

##############################################################################################################

# Load and install needed package
needed_packages <- c("glmmTMB", # for fitting GLMMMs
                     'DHARMa', # for residual diagnostics
                     "bbmle", # for AICtab (model comparison)
                     "dplyr", # for data manipulation
                     "data.table", # enhanced data.frame
                     "ggplot2", # for plotting
                     'phylobase','phylosignal') # for phylogenetic analysis
new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)

# Clean global environment
rm(list=ls()); gc()

# Set working directory
setwd("WRITE HERE YOUR WORKING DIRECTORY")

# STEPS IN THIS SCRIPT 

#  1. Load and understand the dataset.
#  2. Get basic summary data for research effort.
#  3. Prepare variables that will be used, check for multicollinearity and standardise.
#  4. Analyse the data through a Generalized Linear Mixed Effect (GLMM) model.
#  5. Run a GLMM for each biogeographic realm separately.
#  6. Check whether GLMM residuals show phylogenetic autocorrelation.
#  7. Map research effort.
#  8. Perform an additional analysis for threat status as assessed vs. non-assessed.
#  9. Perform an additional analysis using polygon-based within-range predictors.
# 10. Create a bubble plot for grid- and polygon-based model coefficients

#####

# STEP 1 - Load and understand the dataset.
##############################################################################################################
# STEP 1 - Load and understand the dataset.

mydata <- data.table::fread("ResearchEffortData.csv", stringsAsFactors = T, na.strings = "", encoding = 'UTF-8')
names(mydata)

## We have 27 columns in the dataset, each of which is explained below:
# Species: accepeted scientific name based on May 2021 version of the Reptile Database (RDB).
# RDB_sciname: edited scientific name - nomenclature follows the RDB
# MoL_sciname: edited scientific name - nomenclature follows the Map of Life project
# Comb_sciname: edited scientific name combined from RDB and MoL (for matching operations)
# Year: year of species description
# Genus: genus each species belong to 
# Family: family each species belong to 
# Order: order each species belong to           
# Suborder: suborder each species belong to
# Synonyms: all synonyms attributed to a currently valid species; data extracted from the Reptile Database.
# AllUniqueNames: all unique names attributed to each species (including current species name)
# AllAmbigNames: all ambiguos names attributed to each species (including current species name)
# Npapers1stApproach: number of papers extracted from Scopus using 1st approach (see main text)
# Npapers2ndApproach: number of papers extracted from Scopus using 2nd approach (see main text)
# Npapers3rdApproach: number of papers extracted from Scopus using 3rd approach (see main text)
# RangeSize: range size as the number of 110x110km grid cells 
# RangeRarity: a metric of endemism richness (see main text)
# Elevation: avraged elevation across the range of each species
# HumanDensity: averaged human density across the range of each species
# TotBiodInst: total number of biodiversity institutions within the range of each species
# Avg_StdPPP: weighted purchasing power based on species-country intersactions (see main text)
# MainBiogRealm: main biogeographic realm for each species based on the number of grid cells in each region
# ThreatStatus: a combination of IUCN and predicted (IUNC) threat statuses (see main text)
# BodyMass: body mass in grams for each species
# SourceBodyMass: references from where body mass information was retrieved
# Habit: species' habit (fossorial, terrestrial, aquatic, arboreal, or a combination of such levels)
# SourceHabitat: references from where habit information was retrieved

# Quickly inspect missing data across columns
NA_Prop<-round((colSums(is.na(mydata))/nrow(mydata)), digits = 3) # compute the % of NA in each column
NA_Prop

# Plot the relationshipt between the three different search approaches
(f1 <- ggplot(mydata, (aes(y=log(Npapers1stApproach+1), x=log(Npapers2ndApproach+1)))) + 
    geom_point(size = .7) + 
    ylab('log N papers 1st approach') + xlab('log N papers 2nd approach') +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=8, face="bold"),
          axis.line = element_line(colour="black"),
          axis.text = element_text(size=6, colour = "black")) +
    # show the corresponding linear equations and the R2 by groups
    ggpubr::stat_cor(method = "spearman"))

(f2 <- ggplot(mydata, (aes(y=log(Npapers1stApproach+1), x=log(Npapers3rdApproach+1)))) + 
    geom_point(size = .7) + 
    ylab('log N papers 1st approach') + xlab('log N papers 3rd approach') +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=8, face="bold"),
          axis.line = element_line(colour="black"),
          axis.text = element_text(size=6, colour = "black")) +
    # show the corresponding linear equations and the R2 by groups
    ggpubr::stat_cor(method = "spearman"))

(f3 <- ggplot(mydata, (aes(y=log(Npapers3rdApproach+1), x=log(Npapers2ndApproach+1)))) + 
    geom_point(size = .7) + 
    ylab('log N papers 3rd approach') + xlab('log N papers 2nd approach') +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=8, face="bold"),
          axis.line = element_line(colour="black"),
          axis.text = element_text(size=6, colour = "black")) +
    # show the corresponding linear equations and the R2 by groups
    ggpubr::stat_cor(method = "spearman"))

# Make a multipanel plot
(figS3 <- cowplot::plot_grid(f1, f2, f3, align = "hv", nrow = 3, ncol = 1))
# Export the figure
ggsave(paste0(getwd(), "/FigureS3.NpapersTypesCorr.pdf"), plot=figS3, width=6, height=5, units="in", dpi = "print")

# The different search approaches are highly correlated. Subsequent analysis will be based only on the 1st approach. 
n_col<-which(names(mydata)=="Npapers1stApproach")[1] # get the number of the column holding the nº of papers from the 1st approch
colnames(mydata)[n_col] <- 'Npapers'

rm(f1,f2,f3,figS1)

#####

# STEP 2 - Get basic summary data for research effort
##############################################################################################################
# STEP 2 - Get basic summary data for research effort

names(mydata)
summary(mydata$Npapers)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000    0.000    1.000    8.478    3.000 2130.000

sum(mydata$Npapers) # a total of 89,280 papers were retrived from the Scopus database

top10 <- dplyr::top_n(mydata, 10, Npapers); sum(top10$Npapers); rm(top10)
# the top-10 species in terms of N papers had a total of 13,449 (15% of all) publications

# Check frequency of species with zero papers
length(which(mydata$Npapers==0)) # 4,027 species with no publication at all.
length(which(mydata$Npapers==0))/length(mydata$Npapers) # 38.2% of all species

# Check spcies with 10 papers or less
length(which(mydata$Npapers<=10)) # 9,529 species with up to 10 publications


# Plot the number of publications across years and label the top-10 most researched species.
(fig1 <- ggplot(mydata, aes(x = Year, y = Npapers, label = Species))+
  geom_point(alpha=.7, size=1)+
  geom_text(aes(label=ifelse(Npapers>800, as.character(Species),'')), size = 1.5, hjust = -0.05, vjust = .05)+
  labs(x = NULL, y = "Number of papers per species")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=7, face="bold"),
        axis.line = element_line(colour="black"),
        axis.text = element_text(size=6, colour='black')))

ggsave(paste0(getwd(), "/Figure1.PapersPerSpp.pdf"), plot=fig1, width=4, height=3, units="in", dpi = "print")
ggsave(paste0(getwd(), "/Figure1.PapersPerSpp.png"), plot=fig1, width=4, height=3, units="in", dpi = "print")
# This image was latter edited in InkScape for minor aesthetics adjustments (see differences in main text).

# Plot the distribution of research effort across realms
(figS4 <- ggplot(mydata[!is.na(mydata$MainBiogRealm), ], 
             aes(x = log10(Npapers1stApproach+1)))+
    geom_histogram(fill='grey50', color='grey50')+
    scale_fill_brewer(type = 'qual', palette = 2)+
    labs(x=expression(Log[10]~ "number of papers"), y="Count")+
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.position = 'none')+
    facet_wrap(~MainBiogRealm, nrow = 4, ncol = 2, scales = "free_y"))

# save figure - add 'device = cairo_pdf' to save special character in pdf, otherwise there is some bug 
ggsave(paste0(getwd(), "/FigureS4.ResEffortByRealm.pdf"), plot=figS4, width=5, height=6, units="in", bg = 'transparent', dpi = "print", device = cairo_pdf)


# Summarize and plot the number of papers per family 
PapersByFamily <- mydata%>%
  group_by(Suborder, Family)%>%
  summarise(nSpp = sum(!is.na(Species)),
            total = sum(Npapers, na.rm = T),
            mean = mean(Npapers, na.rm = T),
            median = median(Npapers, na.rm = T),
            sd = sd(Npapers))%>%
  mutate(se = sd / sqrt(nSpp))
#fwrite(PapersByFamily, 'PapersByFamily.csv') # Table S2 in suplementary material

(figS5a <- ggplot(PapersByFamily, aes(x = log(nSpp), y = log(total+1)))+
    geom_point(alpha=.7, size = 1.5)+
    geom_smooth(method = 'lm', fullrange = FALSE, col = 'black')+
    labs(x = NULL, y = "Total number of papers \\nper family")+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=8, face="bold"),
          axis.line = element_line(colour="black"),
          axis.text = element_text(size=6, colour='black'))+
    ggpubr::stat_cor(method = "spearman", cor.coef.name = 'rho', size = 3,
                     label.x.npc = 0.7, label.y.npc = 0.2, p.accuracy = 0.0001))

(figS5b <- ggplot(PapersByFamily, aes(x = log(nSpp), y = log(mean+1)))+
    geom_pointrange(aes(ymin = log(mean-se+1), ymax = log(mean+se+1)), size = .25, alpha = 0.7)+
    geom_smooth(method = 'lm', fullrange = FALSE, col = 'black')+
    labs(x = "Number of species per family", y = "Avg. number of papers \\nper species per family")+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=8, face="bold"),
          axis.line = element_line(colour="black"),
          axis.text = element_text(size=6, colour='black'))+
    ggpubr::stat_cor(method = "spearman", cor.coef.name = 'rho', size = 3,
                     label.x.npc = 0.7, label.y.npc = 0.8, p.accuracy = 0.0001))


(figS4 <- cowplot::plot_grid(figS5a, figS5b, align = "hv", nrow = 2, ncol = 1))
# Export the figure
ggsave(paste0(getwd(), "/FigureS5.PapersByFamily.pdf"), plot=figS5, width=4, height=4, units="in", dpi = "print")


# Make a circular barplot showing the total number of publications per family,
# median research effort per species within each family, and
# a density plot of overall research effort.

library(tidyverse)

data <- as.data.frame(PapersByFamily)

# Reange suborder order to match reptile phylogeny
levels(data$Suborder)
data$Suborder <- factor(data$Suborder,
                        levels = c("Crocodylia","Cryptodira","Pleurodira","Rhynchocephalia",
                                   "Dibamidae","Gekkota","Scincoidea","Lacertiformes",
                                   "Anguimorpha","Iguania","Serpentes"))

# Remove unecessary columns
data <- data[ , c('Suborder','Family','nSpp','total','median')]
data$median <- round(log10(data$median+1), 2) # log-10 transform and round avg. number of papers
data$total <- round(log10(data$total+1), 2) # log-10 transform and round total number of papers

# Order data in drecreasing total number of publications
data <- data %>% arrange(Suborder, desc(total))

# Set a number of 'empty bar' (= graphic space) to add at the end of each group
empty_bar <- 1
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$Suborder), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$Suborder <- rep(levels(data$Suborder), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(Suborder)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
# Id substractd 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar 
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Prepare a data frame for base lines
base_data <- data %>% 
  group_by(Suborder) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=median(c(start, end)))

# Prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]
summary(grid_data)

# Set automatic breaks for the axes and colorbars (for the fill component):
myBreaks<-c(ceiling(min(data$median, na.rm=T)*100)/100, 
            floor(max(data$median, na.rm=T)*100)/100)
myBreaks<-c(myBreaks[1], ((myBreaks[2]+myBreaks[1])/2), myBreaks[2])

# Check maximun and minum values
range(data$total, na.rm = T) # 0.00; 4.21
range(data$median, na.rm = T) # 0.00 2.86

# Make the circular plot
(fig2 <- ggplot() + 
    
    # Get the graphic space properly: 
    geom_bar(data=data, aes(x=as.factor(id), y=total, fill=median), color="black", lwd=0.1, stat="identity") +
    
    # Add thin dotted line segments:
    geom_segment(data=grid_data, aes(x=0, y=4, xend=nrow(data), yend=4), colour="#BEBEBE", alpha=.5, size=0.3, linetype ="dotted", inherit.aes=FALSE) +
    geom_segment(data=grid_data, aes(x=0, y=3, xend=nrow(data), yend=3), colour="#BEBEBE", alpha=.5, size=0.3, linetype ="dotted", inherit.aes=FALSE) +
    geom_segment(data=grid_data, aes(x=0, y=2, xend=nrow(data), yend=2), colour="#BEBEBE", alpha=.5, size=0.3, linetype ="dotted", inherit.aes=FALSE) +
    geom_segment(data=grid_data, aes(x=0, y=1, xend=nrow(data), yend=1), colour="#BEBEBE", alpha=.5, size=0.3, linetype ="dotted", inherit.aes=FALSE) +
    geom_segment(data=grid_data, aes(x=0, y=0, xend=nrow(data), yend=0), colour="#BEBEBE", alpha=.5, size=0.3, linetype ="dotted", inherit.aes=FALSE) +
    
    # Replot bars on the top of the segments plotted before:
    geom_bar(data=data, aes(x=as.factor(id), y=total, fill=median), color="black", lwd=0.1, stat="identity") +
    #scale_fill_gradientn(colours=colorspace::diverge_hcl(7),na.value='transparent') + 
    scale_fill_gradientn(colours=colorspace::sequential_hcl(n=7, h1=0, c1=0, l1=10, l2=98, p1=1.3, rev = TRUE)) + 
    
    geom_text(data=base_data, aes(x = title, y = -1.5, label=Suborder), hjust=0.5, 
              colour = "black", alpha=0.8, size=2.5, fontface="bold", inherit.aes=F) +
    #geom_text(data=label_data, aes(x=id, y=max(data$total,na.rm = T)+0.05, label=paste(Family, '(', nSpp, ')'), hjust=hjust), color="black",
     #         fontface="plain", alpha=0.9, size=2.7, angle=label_data$angle, inherit.aes=T) +
    geom_text(data=label_data, aes(x=id, y=data$total+0.1, label=paste0(Family, ' (',nSpp,')'), hjust=hjust), color="black",
              fontface="plain", alpha=0.9, size=2.5, angle=label_data$angle, inherit.aes=F) +
  
    geom_segment(data=base_data, aes(x = start, y = -.125, xend = end, yend = -.125), colour="black", alpha=0.8, size=0.6, inherit.aes=FALSE) +
    
    # Add the tick mark values in the y-axis (add raw values instead of log10-transformed ones to ease interpretation):
    annotate("text", x = rep((max(data$id)+1),5), y = c(0.0, 1, 2, 3, 4), label = c("1", "10", "100", "1,000", "10,000"), color="black", size=2 , angle=0, fontface="bold", hjust=1) +
    annotate("text", label=bquote('Median number of\\n papers per species'), x=0, y=-5, color="black", fontface="bold", size=2.3, angle=0, hjust=1) +
    ylim(-6, 6) +
    
    # Customize the theme:
    theme_minimal() +
    theme(
      legend.position = 'right',
      legend.justification = c(0.5, 0.5),
      legend.text=element_text(size=rel(0.5), hjust=0.5),
      legend.key.width=unit(0.3,"cm"), 
      legend.key.height=unit(0.3,"cm"),
      legend.background=element_blank(),
      legend.title = element_blank(),
      legend.spacing.x = unit(0.1, 'cm'),
      legend.spacing.y = unit(0.1, 'cm'),
      plot.background=element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(0,5), "cm")) +
    coord_polar()+
    # Adjust colourbar for the legend:
    guides(fill=guide_colourbar(nbin=7, barwidth=.6, barheight=3.5,frame.colour="black",
                                frame.linewidth=1, ticks.linewidth=1,draw.ulim=T, draw.llim=T, 
                                ticks=T, ticks.colour="black")))

ggsave(paste0(getwd(), "/Figure2.TotalAndAvgPapersByFamily.pdf"), plot=fig2, width=8, height=8, units="in", bg = "transparent")

# Get a density plot of overall reaserch effort.
# It will be included in figure 2.
(d.plot <- ggplot(mydata, aes(log10(Npapers+1)))+
    #geom_density(alpha=0.3, color = 'grey40', fill='grey40', outline.type = "full")+
    geom_histogram(alpha=0.5, bins = 15, color = 'grey50', fill='grey50')+
    labs(x=expression(Log[10]~ "number of papers per species"), y = 'Count')+
    theme_classic()+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10)))
ggsave(paste0(getwd(), "/Figure2b.ResEffortDistrib.pdf"), plot=d.plot, width=4, height=3, units="in", bg = "transparent")

# Figure latter edited in the software Inkscape for minor aestethics adjustments.

# Clean work environment
rm(list=setdiff(ls(), c("mydata")))

#####

# STEP 3 - Prepare variables that will be used, check for multicollinearity and standardise.
##############################################################################################################
# STEP 3 - Prepare variables that will be used, check for multicollinearity and standardise.

# Check for skewed distributions and kurtosis of predictors (transform data if necessary).
# Values for skewness and kurtosis between -2 and +2 are considered acceptable
# to prove normal univariate distribution (George & Mallery, 2010) 
e1071::skewness(mydata$Year, na.rm = T); e1071::kurtosis(mydata$Year, na.rm = T) # okay
e1071::skewness(mydata$BodyMass, na.rm = T); e1071::kurtosis(mydata$BodyMass, na.rm = T) 
e1071::skewness(mydata$RangeSize, na.rm = T); e1071::kurtosis(mydata$RangeSize, na.rm = T) 
e1071::skewness(mydata$RangeRarity, na.rm = T); e1071::kurtosis(mydata$RangeRarity, na.rm = T) 
e1071::skewness(mydata$Elevation, na.rm = T); e1071::kurtosis(mydata$Elevation, na.rm = T)
e1071::skewness(mydata$HumanDensity, na.rm = T); e1071::kurtosis(mydata$HumanDensity, na.rm = T)
e1071::skewness(mydata$TotBiodInst, na.rm = T); e1071::kurtosis(mydata$TotBiodInst, na.rm = T)
e1071::skewness(mydata$Avg_StdPPP, na.rm = T); e1071::kurtosis(mydata$Avg_StdPPP, na.rm = T)
# Conclusion: log10 transform all predictors but year of description

# Log-10 transform and standardize (z-tranformation) predictors (mean 0 and SD = 1):
for (i in 1) {
  mydata$StdYear <- scale(mydata$Year, center = T, scale = T)
  mydata$LogBodyMass <- scale(log10(mydata$BodyMass), center = T, scale = T)
  mydata$LogRangeSize <- scale(log10(mydata$RangeSize), center = T, scale = T)
  mydata$LogRangeRarity <- scale(log10(mydata$RangeRarity), center = T, scale = T)
  mydata$LogElevation <- scale(log10(mydata$Elevation+1), center = T, scale = T)
  mydata$LogHumD <- scale(log10(mydata$HumanDensity+1), center = T, scale = T)
  mydata$LogNumBioInst <- scale(log10(mydata$TotBiodInst+1), center = T, scale = T)
  mydata$LogStdPPP <- scale(log10(mydata$Avg_StdPPP), center = T, scale = T)
}

# Check categorical predictors that will be tested:
names(mydata)

### Habit ###
levels(mydata$Habit)

# Convert habitat into a metric of verticality, where:
# fos = 0
# fos & aqu OR fos & ter = 0.25
# aqu OR ter OR ter & aqu OR fos & arb OR fos & ter & arb = 0.5 
# ter & aqu & arb OR aqu & arb OR ter & arb = 0.75
# arb = 1
mydata$Verticality <- factor(mydata$Habit,labels = c(0.5,0.75,1,0,0.25,0.5,0.25,0.25,0.5,0.5,0.5,0.75,0.75))
mydata$Verticality <- as.numeric(as.character(mydata$Verticality)) # convert to numeric
mydata$Verticality <- scale(mydata$Verticality, center = T, scale = T) # z-transform variable to make it comparable

### IUCN status ###
summary(mydata$ThreatStatus)
levels(mydata$ThreatStatus) 

# IUCN categories are an ordinal index:
# LC (1), NT (2), VU (3), EN (4), CR (5), EW(6). 
# Note: extinct and data deficient species are removed (see main text for explanation).
mydata$IUCN <- as.integer(as.character(factor(mydata$ThreatStatus,labels = c(5,NA,4,6,NA,1,2,3))))
mydata$IUCN <- scale(mydata$IUCN, center = T, scale = T) # z-transform variable to make it comparable

# Check multicollinearity among the predictor variables:
vars <- mydata[ , StdYear:IUCN]
vars$Verticality2 <- vars$Verticality^2
usdm::vif(vars); rm(vars) # Table S1 in supplementary material
# Conclusion: continuous predictors show low (VIF < 5) multicollinearity (keep them all)

### Biog. Realm ###
levels(mydata$MainBiogRealm) # Australasia is missing an 'a'; fixing it next 
mydata$MainBiogRealm <- factor(mydata$MainBiogRealm,
                               labels = c("Afrotropic","Australasia","IndoMalay","Neartic","Neotropic","Oceania","Paleartic"))

# Select only response, fixed, random, and grouping variables as a new dataset
names(mydata)
dataset <- mydata%>%
  select(Comb_sciname,Npapers,Family,MainBiogRealm,StdYear:IUCN)
summary(dataset)

# Remove NAs in the new dataset
dataset <- droplevels(dataset[complete.cases(dataset) , ]) # n = 7,102 species

# save the new dataset
save(dataset, file='Data.Rdata')

#####

# STEP 4 - Analyse the data through a Generalized Linear Mixed Effect (GLMM) model.
##############################################################################################################
# STEP 4 - Analyse the data through a global Generalized Linear Mixed Effect Model.

names(dataset)

# model formula
formula <- Npapers ~ StdYear+LogBodyMass+LogRangeSize+LogRangeRarity+LogElevation+LogHumD+
  LogStdPPP+LogNumBioInst+Verticality+I(Verticality^2)+IUCN+(1|Family)


### 1 - fit a poisson model
poismod <- glmmTMB(formula = formula,
                   data = dataset,
                   family = poisson)

# Check how well the model fits the data
simulateResiduals(fittedModel = poismod, plot = T, n = 1000) 
testDispersion(poismod, plot = F, type = 'PearsonChisq', alternative = 'greater') # overdispersion detected

### 2 - fit a negative binomial model (to account for overdispersion)
nbmod <- glmmTMB(formula = formula,
                 data = dataset,
                 family = nbinom2) # assumes the variance increases quadratically with the mean.

# Model validation 
simulateResiduals(fittedModel = nbmod, plot = T, n = 1000) # Much better fit than the Poisson model
performance::r2(nbmod) # R2 (cond) = 0.76, R2 (marg) = 0.59 

summary(nbmod) # model results
save(nbmod, file = 'nbGLMM_global.Rdata')

rm(list=setdiff(ls(), c("mydata", "dataset", "nbmod")))

#####

# STEP 5 - Run a GLMM for each biogeographic realm separately.
##############################################################################################################
# STEP 5 - Run a GLMM for each biogeographic realm separately.

names(dataset)

# Extract coefficients from the global model - will be added to the realm-specific models
# to create a tree plot with model coefficients.
#load("nbGLMM_global.Rdata")
res <- summary(nbmod)

# store coefficients
out.global <- as.data.frame(res$coefficients$cond)
out.global$Predictors <- rownames(out.global) # add predictor names into a new column
out.global <- mutate(out.global, 
                     Realm = 'Global', # add realm names into a new column
                     lower95 = Estimate - 1.96 * `Std. Error`, # get lower CI95%
                     upper95 = Estimate + 1.96 * `Std. Error`) # get upper CI95%
out.global <- out.global[ , c("Realm","Predictors","Estimate","Std. Error","z value","Pr(>|z|)","lower95","upper95")]


### FIT NEGATIVE BINOMIAL GLMMs FOR EACH BIOGEOGRAPHIC REALM ###

# Model formula
formula <- Npapers ~ StdYear+LogBodyMass+LogRangeSize+LogRangeRarity+LogElevation+LogHumD+
  LogStdPPP+LogNumBioInst+Verticality+I(Verticality^2)+IUCN+(1|Family)

# Create a vector holdingn biog. realm names to iterate through
levels(dataset$MainBiogRealm)
biorealm <- levels(dataset$MainBiogRealm)
biorealm <- biorealm[-6] # remove Oceania; too small sample size

# Create an empty list to store the outputs:
outs <- list()

for (i in seq_along(biorealm)) {
  # create subdataset containing only observations from a single biog. realm
  subdataset = droplevels(dataset[dataset$MainBiogRealm==biorealm[i] , ])
  
  # fit the NB mixed model
  nbglmm = glmmTMB(formula = formula, data = subdataset, family = nbinom2)
  
  print(paste0(biorealm[j], " ", performance::r2(nbglmm)))
  
  # residual diagnostics
  pdf(paste("ResidPlot", biorealm[i], '.pdf', sep = '_'))
  simulationOutput = simulateResiduals(fittedModel = nbglmm, plot = T)
  dev.off()
  
  # save model 
  save(nbglmm, file = paste("nbGLMM", biorealm[i], '.Rdata', sep = '_'))
  
  # get result
  res = summary(nbglmm)
  
  # store coefficients
  coefs = as.data.frame(res$coefficients$cond)
  coefs$Predictors <- rownames(coefs) # add predictor names into a new column
  coefs = mutate(coefs, 
                 Realm = paste(biorealm[i]), # add realm names into a new column
                 lower95 = Estimate - 1.96 * `Std. Error`, # get lower CI95%
                 upper95 = Estimate + 1.96 * `Std. Error`) # get upper CI95%
  coefs = coefs[ , c("Realm","Predictors","Estimate","Std. Error","z value","Pr(>|z|)","lower95","upper95")]
  
  outs[[i]] <- coefs
  
  rm(subdataset, nbglmm, simulationOutput, coefs)
  
}

# Pseudo R2:
# Afrotropic: Conditional R2 = 0.69458750696114; Marginal R2 = 0.503410783704316 
# Australasia: Conditional R2 = 0.777046613403484; Marginal R2 = 0.60761540174943   
# IndoMalay: Conditional R2 = 0.78650810101583; Marginal R2 = 0.610865411100905 
# Neartic: Conditional R2 = 0.832044878983144; Marginal R2 = 0.681750793506492  
# Neotropic: Conditional R2 = 0.704783779296954; Marginal R2 = 0.544241381264331    
# Paleartic: Conditional R2 = 0.819897416399419; Marginal R2 = 0.691803415472978  

results <- bind_rows(outs) # combine results of each biorealm
results <- bind_rows(outs, out.global) # join with results of the global model
results[,c(3:5,7:8)] <- round(results[,c(3:5,7:8)], digits = 3)
results[,6] <- round(results[,6], digits = 4)

fwrite(results, 'gridBasedModelResults.csv') # save model coefficients

rm(list=setdiff(ls(), c("mydata", "dataset", "results")))

#####

# STEP 6 - Check whether GLMM residuals show phylogenetic autocorrelation.
##############################################################################################################
# STEP 6 - Check whether GLMM residuals show phylogenetic autocorrelation.

# NOTE: The following analyses were performed in a cluster using 15 cores.
# Running it on a personal computer may take much pc memory and time.

# Load and install needed package
needed_packages <- c('glmmTMB', # for model fitting
                     'dplyr', # for data manipulation
                     'data.table', # enhanced data.frame
                     'foreach', # for looping construct
                     'doParallel', # for parallel computing
                     'snow', # for setting cl
                     'geiger','phytools',"phylobase","phylosignal") # for phylogenetic analysis
new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)

# Clean global environment
rm(list=ls()); gc()

# Set working directory
setwd("DEFINE HERE YOUR WORKING DIRECTORY")

# Load the 100 fully-sampled phylogeny trimmed to include only species in the global dataset
load("Subset100Supertrees_7102taxaReptiles.Rdata")

#--------------------------------------------------------------------------------#
# i - Phylogenetic correlation in residuals of the global model
#--------------------------------------------------------------------------------#

cl <- snow::makeMPIcluster(count = 15)

# Load needed files
load("Data.Rdata") # dataset used in to fit the negative binomial GLMM
load('nbGLMM_global.Rdata') # GLMM result for the global model

# Extract residuals of the final model:
mod_resid<-resid(nbmod)

# Add model resids into dataset
dataset <- as.data.frame(cbind(dataset, mod_resid)) 
rownames(dataset) <- dataset[ , 'Comb_sciname'] # spp name as row name

# Check if we have the same spp in our data as in the tree
obj <- name.check(supertrees[[1]], dataset)

# Prepare workspace for parallel computing:
registerDoParallel(cl)
getDoParWorkers()

# NOTE: running time can be too long (few days) for the global model. Therefore, we suggest
# partitioning the analyses - e.g., using i = 1:2, 3:4, ..., 9-10, then 
# merging correlograms into a single dataset.

{
  a<-Sys.time(); a
  
  PhyCorrGlobal<-foreach(i = 1:10, 
                         .export = 'rbind',
                         .packages = c("data.table", "phylobase", "phylosignal"))  %dopar% {
                           
                           # Select one trimmed fully-sampled tree:
                           my_tree<-supertrees[[i]]
                           
                           # Create a phylo4 object including GLMM model residuals:
                           phylo4d_filter<-phylobase::phylo4d(x=my_tree, data.frame(GLMM_resid=dataset$mod_resid))
                           
                           # Compute the phylogenetic correlogram:
                           phy.cor<-phylosignal::phyloCorrelogram(p4d=phylo4d_filter, trait=names(tdata(phylo4d_filter)),
                                                                  dist.phylo="patristic", n.points=14, ci.bs=100)
                           
                           correlogram_data<-as.data.frame(phy.cor[[1]])
                           names(correlogram_data)<-c("dist.class", "lower_ci", "upper_ci", "coef")
                           correlogram_data$Iter<-i
                           correlogram_data$N_class<-1:14
                           correlogram_data
                           
                         }
  snow::stopCluster(cl) # terminate cluster
  
  b<-Sys.time(); b-a
  
}

# Extract the average correlogram output across iterations:
PhyCorrGlobal<-as.data.table(rbindlist(PhyCorrGlobal))
AvgPhyCorrGlobal<-PhyCorrGlobal[, .(Distance=mean(dist.class, na.rm=T),
                                    Lower_CI=mean(lower_ci, na.rm=T),
                                    Upper_CI=mean(upper_ci, na.rm=T),
                                    MoranI_coef=mean(coef, na.rm=T)),
                                by = .(N_class)]

# Export the results:
save(PhyCorrGlobal, AvgPhyCorrGlobal, file="Correlogram_global.Rdata")

rm(list=setdiff(ls(), c('tree_set')))

#--------------------------------------------------------------------------------#
# ii - Phylogenetic correlation in residuals of the Afrotropical model 
#--------------------------------------------------------------------------------#

cl <- snow::makeMPIcluster(count = 15)

# Load needed files
load("Data.Rdata") # dataset used in to fit the GLMM
load("nbGLMM_Afrotropic_.Rdata") # GLMM result for the Afrotropic realm

# Add model resids into dataset
subdataset <- droplevels(as.data.frame(dataset[dataset$MainBiogRealm=='Afrotropic' , ]))

# Extract residuals of the final model:
mod_resid<-resid(nbglmm)
dataset <- as.data.frame(cbind(subdataset, mod_resid)) 

# Update taxonomic names to match those used in the phylogeny:
rownames(dataset) <- dataset[ , 'Comb_sciname'] # spp name as row name

# Check if we have the same spp in our data as in the tree
obj <- name.check(supertrees[[1]], dataset)

# Crop tip labels based on species on the dataset
phylo_tree_set<-list()

# Remove spp that are in the phylogeny but not in the dataset
for (i in 1:100) {
  phylo_tree_set[[i]] <- drop.tip(supertrees[[i]], obj$tree_not_data)
} 
name.check(phylo_tree_set[[1]], dataset)
rm(obj)

# Prepare workspace for parallel computing:
registerDoParallel(cl)
getDoParWorkers()

{
  
  a<-Sys.time(); a
  
  PhyCorrAfro<-foreach(i = 1:100, 
                       .export = 'rbind',
                       .packages = c("data.table", "phylobase", "phylosignal"))  %dopar% {
                         
                         # Select one trimmed fully-sampled tree:
                         my_tree<-phylo_tree_set[[i]]
                         
                         # Create a phylo4 object including GLMM model residuals:
                         phylo4d_filter<-phylobase::phylo4d(x=my_tree, data.frame(GLMM_resid=dataset$mod_resid))
                         
                         # Compute the phylogenetic correlogram:
                         phy.cor<-phylosignal::phyloCorrelogram(p4d=phylo4d_filter, trait=names(tdata(phylo4d_filter)),
                                                                dist.phylo="patristic", n.points=14, ci.bs=100)
                         
                         correlogram_data<-as.data.frame(phy.cor[[1]])
                         names(correlogram_data)<-c("dist.class", "lower_ci", "upper_ci", "coef")
                         correlogram_data$Iter<-i
                         correlogram_data$N_class<-1:14
                         correlogram_data
                         
                       }
  
  snow::stopCluster(cl) # terminate cluster
  
  b<-Sys.time(); b-a # Time difference of 18.9 hours for 100 trees using 15 cores
  
} 

# Extract the average correlogram output across iterations:
PhyCorrAfro<-as.data.table(rbindlist(PhyCorrAfro))
AvgPhyCorrAfro<-PhyCorrAfro[, .(Distance=mean(dist.class, na.rm=T),
                                Lower_CI=mean(lower_ci, na.rm=T),
                                Upper_CI=mean(upper_ci, na.rm=T),
                                MoranI_coef=mean(coef, na.rm=T)),
                            by = .(N_class)]

# Export the results:
save(PhyCorrAfro, AvgPhyCorrAfro, file="Correlogram_Afro.Rdata")

rm(list=setdiff(ls(), c('supertrees')))

#--------------------------------------------------------------------------------#
# iii - Phylogenetic correlation in residuals of the Australasia model 
#--------------------------------------------------------------------------------#

cl <- snow::makeMPIcluster(count = 15)

# Load needed files
load("Data.Rdata") # dataset used in to fit the GLMM
load("nbGLMM_Australasia_.Rdata") # GLMM result for the Australasia realm

# Add model resids into dataset
subdataset <- droplevels(as.data.frame(dataset[dataset$MainBiogRealm=='Australasia' , ]))

# Extract residuals of the final model:
mod_resid<-resid(nbglmm)
dataset <- as.data.frame(cbind(subdataset, mod_resid)) 
rownames(dataset) <- dataset[ , 'Comb_sciname'] # spp name as row name

# Check if we have the same spp in our data as in the tree
obj <- name.check(supertrees[[1]], dataset)
#dataset <- dataset[!dataset$Comb_sciname %in% obj$data_not_tree, ] # remove spp not in the tree

# Crop tip labels based on species on the dataset
phylo_tree_set<-list()

# Remove spp that are in the phylogeny but not in the dataset
for (i in 1:100) {
  phylo_tree_set[[i]] <- drop.tip(supertrees[[i]], obj$tree_not_data)
} 
name.check(phylo_tree_set[[1]], dataset)
rm(obj)

# Prepare workspace for parallel computing:
registerDoParallel(cl)
getDoParWorkers()

{
  a<-Sys.time(); a
  
  PhyCorrAust<-foreach(i = 1:100, 
                       .export = 'rbind',
                       .packages = c("data.table", "phylobase", "phylosignal"))  %dopar% {
                         
                         # Select one trimmed fully-sampled tree:
                         my_tree<-phylo_tree_set[[i]]
                         
                         # Create a phylo4 object including GLMM model residuals:
                         phylo4d_filter<-phylobase::phylo4d(x=my_tree, data.frame(GLMM_resid=dataset$mod_resid))
                         
                         # Compute the phylogenetic correlogram:
                         phy.cor<-phylosignal::phyloCorrelogram(p4d=phylo4d_filter, trait=names(tdata(phylo4d_filter)),
                                                                dist.phylo="patristic", n.points=14, ci.bs=100)
                         
                         correlogram_data<-as.data.frame(phy.cor[[1]])
                         names(correlogram_data)<-c("dist.class", "lower_ci", "upper_ci", "coef")
                         correlogram_data$Iter<-i
                         correlogram_data$N_class<-1:14
                         correlogram_data
                         
                       }
  stopCluster(cl) # terminate cluster
  
  b<-Sys.time(); b-a # Time difference of 7.058551 hours using 15 cores
}

# Extract the average correlogram output across iterations:
PhyCorrAust<-as.data.table(rbindlist(PhyCorrAust))
AvgPhyCorrAust<-PhyCorrAust[, .(Distance=mean(dist.class, na.rm=T),
                                Lower_CI=mean(lower_ci, na.rm=T),
                                Upper_CI=mean(upper_ci, na.rm=T),
                                MoranI_coef=mean(coef, na.rm=T)),
                            by = .(N_class)]

# Export the results:
save(PhyCorrAust, AvgPhyCorrAust, file="Correlogram_Aust.Rdata")

rm(list=setdiff(ls(), c('supertrees')))

#--------------------------------------------------------------------------------#
# iv - Phylogenetic correlation in residuals of the IndoMalay model 
#--------------------------------------------------------------------------------#

cl <- snow::makeMPIcluster(count = 15)

# Load needed files
load("Data.Rdata") # dataset used in to fit the GLMM
load("nbGLMM_IndoMalay_.Rdata") # GLMM result for the IndoMalay realm

# Add model resids into dataset
subdataset <- droplevels(as.data.frame(dataset[dataset$MainBiogRealm=='IndoMalay' , ]))

# Extract residuals of the final model:
mod_resid<-resid(nbglmm)
dataset <- as.data.frame(cbind(subdataset, mod_resid)) 
rownames(dataset) <- dataset[ , 'Comb_sciname'] # spp name as row name

# Check if we have the same spp in our data as in the tree
obj <- name.check(supertrees[[1]], dataset)

# Crop tip labels based on species on the dataset
phylo_tree_set<-list()

# Remove spp that are in the phylogeny but not in the dataset
for (i in 1:100) {
  phylo_tree_set[[i]] <- drop.tip(supertrees[[i]], obj$tree_not_data)
}
name.check(phylo_tree_set[[1]], dataset)
rm(obj, nbgmlm)

# Prepare workspace for parallel computing:
registerDoParallel(cl)
getDoParWorkers()

{
  a<-Sys.time(); a
  
  PhyCorrIndo<-foreach(i = 1:100, 
                       .export = 'rbind',
                       .packages = c("data.table", "phylobase", "phylosignal"))  %dopar% {
                         
                         # Select one trimmed fully-sampled tree:
                         my_tree<-phylo_tree_set[[i]]
                         
                         # Create a phylo4 object including GLMM model residuals:
                         phylo4d_filter<-phylobase::phylo4d(x=my_tree, data.frame(GLMM_resid=dataset$mod_resid))
                         
                         # Compute the phylogenetic correlogram:
                         phy.cor<-phylosignal::phyloCorrelogram(p4d=phylo4d_filter, trait=names(tdata(phylo4d_filter)),
                                                                dist.phylo="patristic", n.points=14, ci.bs=100)
                         
                         correlogram_data<-as.data.frame(phy.cor[[1]])
                         names(correlogram_data)<-c("dist.class", "lower_ci", "upper_ci", "coef")
                         correlogram_data$Iter<-i
                         correlogram_data$N_class<-1:14
                         correlogram_data
                         
                       }
  stopCluster(cl) # terminate cluster
  
  b<-Sys.time(); b-a # Time difference of 1.12 days
  
}

# Extract the average correlogram output across iterations:
PhyCorrIndo<-as.data.table(rbindlist(PhyCorrIndo))
AvgPhyCorrIndo<-PhyCorrIndo[, .(Distance=mean(dist.class, na.rm=T),
                                Lower_CI=mean(lower_ci, na.rm=T),
                                Upper_CI=mean(upper_ci, na.rm=T),
                                MoranI_coef=mean(coef, na.rm=T)),
                            by = .(N_class)]

# Export the results:
save(PhyCorrIndo, AvgPhyCorrIndo, file="Correlogram_Indo.Rdata")

rm(list=setdiff(ls(), c('supertrees')))

#--------------------------------------------------------------------------------#
# v - Phylogenetic correlation in residuals of the Neartic model 
#--------------------------------------------------------------------------------#

cl <- snow::makeMPIcluster(count = 15)

# Load needed files
load("Data.Rdata") # dataset used in to fit the GLMM
load("nbGLMM_Neartic_.Rdata") # GLMM result for the Neartic realm

# Add model resids into dataset
subdataset <- droplevels(as.data.frame(dataset[dataset$MainBiogRealm=='Neartic' , ]))

# Extract residuals of the final model:
mod_resid<-resid(nbglmm)
dataset <- as.data.frame(cbind(subdataset, mod_resid)) 
rownames(dataset) <- dataset[ , 'Comb_sciname'] # spp name as row name

# Check if we have the same spp in our data as in the tree
obj <- name.check(supertrees[[1]], dataset)

# Crop tip labels based on species on the dataset
phylo_tree_set<-list()

# Remove spp that are in the phylogeny but not in the dataset
for (i in 1:100) {
  phylo_tree_set[[i]] <- drop.tip(supertrees[[i]], obj$tree_not_data)
} 
name.check(phylo_tree_set[[1]], dataset)
rm(obj, nbglmm)

# Prepare workspace for parallel computing:
registerDoParallel(cl)
getDoParWorkers()

{
  a<-Sys.time(); a
  
  PhyCorrNea<-foreach(i = 1:100, 
                      .export = 'rbind',
                      .packages = c("data.table", "phylobase", "phylosignal"))  %dopar% {
                        
                        # Select one trimmed fully-sampled tree:
                        my_tree<-phylo_tree_set[[i]]
                        
                        # Create a phylo4 object including GLMM model residuals:
                        phylo4d_filter<-phylobase::phylo4d(x=my_tree, data.frame(GLMM_resid=dataset$mod_resid))
                        
                        # Compute the phylogenetic correlogram:
                        phy.cor<-phylosignal::phyloCorrelogram(p4d=phylo4d_filter, trait=names(tdata(phylo4d_filter)),
                                                               dist.phylo="patristic", n.points=14, ci.bs=100)
                        
                        correlogram_data<-as.data.frame(phy.cor[[1]])
                        names(correlogram_data)<-c("dist.class", "lower_ci", "upper_ci", "coef")
                        correlogram_data$Iter<-i
                        correlogram_data$N_class<-1:14
                        correlogram_data
                        
                      }
  stopCluster(cl) # terminate cluster
  
  b<-Sys.time(); b-a # Time difference of 24.0 mins
}

# Extract the average correlogram output across iterations:
PhyCorrNea<-as.data.table(rbindlist(PhyCorrNea))
AvgPhyCorrNea<-PhyCorrNea[, .(Distance=mean(dist.class, na.rm=T),
                              Lower_CI=mean(lower_ci, na.rm=T),
                              Upper_CI=mean(upper_ci, na.rm=T),
                              MoranI_coef=mean(coef, na.rm=T)),
                          by = .(N_class)]

# Export the results:
save(PhyCorrNea, AvgPhyCorrNea, file="Correlogram_Near.Rdata")

rm(list=setdiff(ls(), c('supertrees')))

#--------------------------------------------------------------------------------#
# vi - Phylogenetic correlation in residuals of the Neotropic model 
#--------------------------------------------------------------------------------#

cl <- snow::makeMPIcluster(count = 15)

# Load needed files
load("Data.Rdata") # dataset used in to fit the GLMM
load("nbGLMM_Neotropic_.Rdata") # GLMM result for the Neotropic realm

# Add model resids into dataset
subdataset <- droplevels(as.data.frame(dataset[dataset$MainBiogRealm=='Neotropic' , ]))

# Extract residuals of the final model:
mod_resid<-resid(nbglmm)
dataset <- as.data.frame(cbind(subdataset, mod_resid)) 
rownames(dataset) <- dataset[ , 'Comb_sciname'] # spp name as row name

# Check if we have the same spp in our data as in the tree
obj <- name.check(supertrees[[1]], dataset)

# Crop tip labels based on species on the dataset
phylo_tree_set<-list()

# Remove spp that are in the phylogeny but not in the dataset
for (i in 1:100) {
  phylo_tree_set[[i]] <- drop.tip(supertrees[[i]], obj$tree_not_data)
} 
name.check(phylo_tree_set[[1]], dataset)
rm(obj, nbglmm)

# Prepare workspace for parallel computing:
registerDoParallel(cl)
getDoParWorkers()

{
  a<-Sys.time(); a
  
  PhyCorrNeot50<-foreach(i = 1:100, 
                         .export = 'rbind',
                         .packages = c("data.table", "phylobase", "phylosignal"))  %dopar% {
                           
                           # Select one trimmed fully-sampled tree:
                           my_tree<-phylo_tree_set[[i]]
                           
                           # Create a phylo4 object including GLMM model residuals:
                           phylo4d_filter<-phylobase::phylo4d(x=my_tree, data.frame(GLMM_resid=dataset$mod_resid))
                           
                           # Compute the phylogenetic correlogram:
                           phy.cor<-phylosignal::phyloCorrelogram(p4d=phylo4d_filter, trait=names(tdata(phylo4d_filter)),
                                                                  dist.phylo="patristic", n.points=14, ci.bs=100)
                           
                           correlogram_data<-as.data.frame(phy.cor[[1]])
                           names(correlogram_data)<-c("dist.class", "lower_ci", "upper_ci", "coef")
                           correlogram_data$Iter<-i
                           correlogram_data$N_class<-1:14
                           correlogram_data
                           
                         }
  stopCluster(cl) # terminate cluster
  
  b<-Sys.time(); b-a # i:1-50
}

# Extract the average correlogram output across iterations:
PhyCorrNeot<-as.data.table(rbindlist(PhyCorrNeot))
AvgPhyCorrNeot<-PhyCorrNeot[, .(Distance=mean(dist.class, na.rm=T),
                                Lower_CI=mean(lower_ci, na.rm=T),
                                Upper_CI=mean(upper_ci, na.rm=T),
                                MoranI_coef=mean(coef, na.rm=T)),
                            by = .(N_class)]

# Export the results:
save(PhyCorrNeot, AvgPhyCorrNeot, file="Correlogram_Neot.Rdata")

rm(list=setdiff(ls(), c('supertrees')))

#--------------------------------------------------------------------------------#
# vii - Phylogenetic correlation in residuals of the Paleartic model 
#--------------------------------------------------------------------------------#

cl <- snow::makeMPIcluster(count = 15)

# Load needed files
load("Data.Rdata") # dataset used in to fit the GLMM
load("nbGLMM_Paleartic_.Rdata") # GLMM result for the Paleartic realm

# Add model resids into dataset
subdataset <- droplevels(as.data.frame(dataset[dataset$MainBiogRealm=='Paleartic' , ]))

# Extract residuals of the final model:
mod_resid<-resid(nbglmm)
dataset <- as.data.frame(cbind(subdataset, mod_resid)) 
rownames(dataset) <- dataset[ , 'Comb_sciname'] # spp name as row name

# Check if we have the same spp in our data as in the tree
obj <- name.check(supertrees[[1]], dataset)

# Crop tip labels based on species on the dataset
phylo_tree_set<-list()

# Remove spp that are in the phylogeny but not in the dataset
for (i in 1:100) {
  phylo_tree_set[[i]] <- drop.tip(supertrees[[i]], obj$tree_not_data)
} 
name.check(phylo_tree_set[[1]], dataset)
rm(obj, nbglmm)

# Prepare workspace for parallel computing:
registerDoParallel(cl)
getDoParWorkers()

{
  a<-Sys.time(); a
  
  PhyCorrPal<-foreach(i = 1:100, 
                      .export = 'rbind',
                      .packages = c("data.table", "phylobase", "phylosignal"))  %dopar% {
                        
                        # Select one trimmed fully-sampled tree:
                        my_tree<-phylo_tree_set[[i]]
                        
                        # Create a phylo4 object including GLMM model residuals:
                        phylo4d_filter<-phylobase::phylo4d(x=my_tree, data.frame(GLMM_resid=dataset$mod_resid))
                        
                        # Compute the phylogenetic correlogram:
                        phy.cor<-phylosignal::phyloCorrelogram(p4d=phylo4d_filter, trait=names(tdata(phylo4d_filter)),
                                                               dist.phylo="patristic", n.points=14, ci.bs=100)
                        
                        correlogram_data<-as.data.frame(phy.cor[[1]])
                        names(correlogram_data)<-c("dist.class", "lower_ci", "upper_ci", "coef")
                        correlogram_data$Iter<-i
                        correlogram_data$N_class<-1:14
                        correlogram_data
                        
                      }
  stopCluster(cl) # terminate cluster
  
  b<-Sys.time(); b-a # Time difference of 1.15 hours
}

# Extract the average correlogram output across iterations:
PhyCorrPal<-as.data.table(rbindlist(PhyCorrPal))
AvgPhyCorrPal<-PhyCorrPal[, .(Distance=mean(dist.class, na.rm=T),
                              Lower_CI=mean(lower_ci, na.rm=T),
                              Upper_CI=mean(upper_ci, na.rm=T),
                              MoranI_coef=mean(coef, na.rm=T)),
                          by = .(N_class)]

# Export the results:
save(PhyCorrPal, AvgPhyCorrPal, file="Correlogram_Pal.Rdata")


### Plot the correlograms (Figure S6) ### 

setwd("~/Documents/Projetos-Manuscritos/Em andamento/DOUTORADO/Correlates of taxonomic biases/Figures and Tables/Correlograms_Realms_May2022")
load('Correlogram_Afro.Rdata')
load('Correlogram_Aust.Rdata')
load('Correlogram_Indo.Rdata')
load('Correlogram_Near.Rdata')
load('Correlogram_Neot.Rdata')
load('Correlogram_Pal.Rdata')
load('Correlogram_global.Rdata')

Corr_list <- list(AvgPhyCorrAfro,AvgPhyCorrAust,AvgPhyCorrIndo,AvgPhyCorrNea,AvgPhyCorrNeot,AvgPhyCorrPal,AvgPhyCorrGlobal)
MyCorrelograms<-list() # to save plots into a list

for (i in seq_along(Corr_list)) {
  MyCorrelograms[[i]] <- ggplot(Corr_list[[i]], aes(x = Distance, y = MoranI_coef)) +
    geom_point(size=1.5)+
    geom_linerange(aes(ymin = Lower_CI, ymax = Upper_CI))+
    geom_line()+
    geom_hline(yintercept=0, linetype="dashed", color="black") +
    ylim(c(-1, 1)) +
    ylab("Moran's I - GLMM residuals") +
    xlab("Phylogenetic distance (mya)") +
    theme(panel.grid.minor = element_blank(), # remove minor gridlines
          panel.grid.major = element_blank(), # remove major gridlines
          panel.background = element_blank(), # white background
          axis.line = element_line(colour="black"), # axis lines aesthetitcs
          axis.text.y = element_text(hjust=0.5, vjust=0.5, angle=0, size=6),
          axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=6),
          axis.ticks.y=element_blank(),
          axis.title.y=element_text(size=8, colour="black", face="bold", margin=margin(t=0, r=5, b=0, l=0)), # margin between axis.title and axis.values
          axis.title.x=element_text(size=8, colour="black", face="bold", margin=margin(t=5, r=0, b=0, l=0)), # margin between axis.title and axis.values
          legend.position="none",
          plot.background=element_blank(),
          panel.spacing=unit(0,"null"))
}

# Create the multipanel plot with both correlograms:
Multipanel_plot<-ggpubr::ggarrange(MyCorrelograms[[1]],MyCorrelograms[[2]],MyCorrelograms[[3]],
                                   MyCorrelograms[[4]],MyCorrelograms[[5]],MyCorrelograms[[6]],MyCorrelograms[[7]],
                                   labels='auto', align="hv",label.x = .5,
                                   font.label=list(size=10, colour = "black"), ncol=2, nrow=4); Multipanel_plot

# Save to disk:
ggsave(paste0(getwd(), "/FigureS6.PhyloCorrelogram.pdf"), plot=Multipanel_plot, width=6, height=8, units="in", bg = 'transparent', dpi = "print")
# Note: this figure was slightly edited in InkScape for minor aesthetic adjustments  

#####

# STEP 7 - Map research effort.
##############################################################################################################
# STEP 7 - Map research effort.

rm(list=ls())

# Load needed packages
library(rgdal)
library(raster)
library(sp)

# Load assemblage dataset (available in the Dryad repository)
rept_110km <- data.table::fread("ReptSpatialData.csv", stringsAsFactors = T, encoding = 'UTF-8')
rept_110km$Cell_Id110 <- as.factor(as.character(rept_110km$Cell_Id110))
summary(rept_110km)

# Add range size data (from 'mydata') to rept_110km dataset
mydata <- data.table::fread("ResearchEffortData.csv", stringsAsFactors = T, encoding = 'UTF-8', na.strings = c("", NA))
names(mydata)

mydata <- mydata[ , c('Comb_sciname','RangeSize')]
mydata$RangeSize <- log10(mydata$RangeSize+1) # log10 tranform range data; needs +1 to avoid weights of zero
colnames(mydata)[2] <- 'LogRangeSize' 
summary(mydata)

# Merge data with assemblage data using both reptile name sources
AssembData <- left_join(rept_110km, mydata, by = c('GARD_Std_sciname'='Comb_sciname'))
AssembData2 <- left_join(rept_110km, mydata, by = c('RDB_Std_sciname'='Comb_sciname'))
AssembData$LogRangeSize <- ifelse(is.na(AssembData$LogRangeSize),
                     yes=paste(AssembData2$LogRangeSize),
                     no=paste(AssembData$LogRangeSize))
AssembData$LogRangeSize <- as.numeric(AssembData$LogRangeSize)
summary(AssembData$N_papers) # 9,988 NAs
summary(AssembData$LogRangeSize) # 30,949 NAs

# Npapers and range must have the same length, so remove NAs in both columns
AssembData <- AssembData[!is.na(AssembData$LogRangeSize) , ] # this way we also remove all NAs in Npapers

# Multiply the number of papers by the inverse of species'range to obtain a weighted average metric
AssembData$w.papers <- AssembData$N_papers * (1/AssembData$LogRangeSize)

dat <- AssembData %>%
  group_by(Cell_Id110, Long, Lat) %>% # group by cell ID
  summarise(SppRich = NROW(na.omit(N_papers)), # This is needed instead of, e.g. n_distinct(RDB_Std_sciname), becuase the latter dischard spp with NAs in Npapers
            W.AvgPapers = round(mean(w.papers), digits = 2))

# Prepare for plotting research effort in space.
# Define projection to be used
cea<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # cylindrical equal area projection

# Load grid cells (unzip 'gridcells_110km.zip' into a folder called 'Shapefiles' in your working directory):
grid_cells <- rgdal::readOGR(dsn="Shapefiles", layer='gridcells_110km')
raster::crs(grid_cells) # cylindrical equal area grid cells
# Transform coordinates of simple feature
grid_cells_sf <- sf::st_transform(sf::st_as_sf(grid_cells), crs=cea) # certify that data is in the same projection
grid_cells_sf$Cell_Id110 <- as.factor(as.character(grid_cells_sf$Cell_Id110)) # convert grid IDs to factors

# Merge grid cells and data to plot weighted average number of papers
dat <- dat[ , c("Cell_Id110",'W.AvgPapers')]
grid_cells_sf <- left_join(grid_cells_sf, dat, by="Cell_Id110" )

# Load the shapefile of the Ecoregions 2017 (available at: https://ecoregions.appspot.com/):
study_area <- rgdal::readOGR(dsn="Shapefiles", layer='Ecoregions2017') # change directory as needed
study_area <- rmapshaper::ms_simplify(input = study_area) # simplify polygon for faster plotting
raster::crs(study_area)
# Transform native projection into cylindrical equal area projection
study_area <- sp::spTransform(study_area, raster::crs(grid_cells))

# compare projections (make sure both have the same projection)
raster::compareCRS(study_area, grid_cells)

# Convert SpatialPolygonsDataframe into sf objects with CEA projection
study_area_sf <- sf::st_transform(sf::st_as_sf(study_area), crs=cea)

# Map reseach effort
(fig4a <- ggplot2::ggplot() +
    
    # Add the raster layer in the background:
    geom_sf(data=grid_cells_sf, aes(fill=log(W.AvgPapers+1)), colour=NA) +
    
    # All options below are useful to control coloramp:
    scale_fill_gradientn(colours=grDevices::hcl.colors(n = 50, palette="viridis", rev=T), 
                         na.value='transparent') + 
    
    # Add polygon boundaries for the study area:
    geom_sf(data=study_area_sf, fill=NA, colour="black", size=0.3)+
    
    # Specify other aesthetics:
    theme(axis.line=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_blank(),
          plot.background=element_blank(),
          plot.margin=unit(c(0, 0, 0, 0), "cm"),  # top, right, bottom, left
          panel.spacing=unit(c(0, 0, 0, 0), "cm"),  # top, right, bottom, left
          panel.border=element_blank(),
          legend.justification=c(1, 1),
          legend.position=c(0.52, 0.5),
          legend.direction="vertical",
          legend.title=element_text(size=4.5, hjust=0, face="bold"),
          legend.text=element_text(size=4, hjust=0), #size=rel(0.8)
          legend.background=element_blank()
    ) +
    
    # Adjust colourbar for the legend:
    guides(fill=guide_colourbar(title="Weighted avg.\\nnumber of papers",label=T, nbin=5, 
                                barwidth=.5, barheight=3.5,frame.colour="black",
                                frame.linewidth=1, ticks.linewidth=1,draw.ulim=T, draw.llim=T, 
                                ticks=T, ticks.colour="black")))
  
# Export the figure:
ggsave(paste0(getwd(), "/Figure4A.MapNpapers.pdf"), plot=fig4a, width=7, height=5, units="in", dpi = "print")


# Create a line plot with number of papers x latitude 
rept_110km <- data.table::fread("ReptSpatialData.csv", stringsAsFactors = T, encoding = 'UTF-8')
rept_110km$Cell_Id110 <- as.factor(as.character(rept_110km$Cell_Id110))

# Summarize data by latitude based on each species
PapersByLatAndTaxa <- rept_110km%>%
  group_by(Lat, GARD_sciname)%>%
  summarise(nGrids = n(),
            CumPapers = sum(N_papers, na.rm = T))%>%
  mutate(PapersPerGrid = CumPapers/nGrids) # get the # of papers each spp contribute by latitude

PapersByLat <- PapersByLatAndTaxa%>%
  group_by(Lat)%>%
  summarise(nPapers = sum(PapersPerGrid))
# nPapers will still be inflated because the number of papers for given species 
# can be 'shared' across multiple latitude values.
range(PapersByLat$nPapers) # 9 to 34273

(fig4b <- ggplot(PapersByLat, aes(x=Lat, y=(nPapers))) +
    geom_smooth(method = "loess", se = F, col = 'black', span = 0.5)+
    scale_y_continuous(breaks = c(0,10000,20000,30000,40000), 
                       labels = c('0', '10k', '20k', '30k', '40k'),
                       limits = c(0,35000))+
    scale_x_continuous(breaks = c(-40,0,40), labels = c('40ºS', '0º', '40ºN'))+
    labs(x='Latitude',y='Number of publications')+
    theme_classic()+
    coord_flip())
ggsave(paste0(getwd(), "/Figure4B.nPapersByLat.pdf"), plot=fig4b, width=4, height=3, units="in", dpi = "print")

# Both figures were combined in the program Inkscape for minor aesthetic adjustments


## Create box plots showing median research effort per biodiversity hotspots
## and color hotspots based on geographic realm associations

# Load biodiversity hotspots shape file
library(sf)
hotspots <- st_read(paste0(getwd(), '/Shapefiles/hotspots_2016_1.shp'))
# keep only hotspot area (36 in total)
hotspots <- hotspots[hotspots$Type=='hotspot area', ]; hotspots$Type <- NULL
crs(hotspots) # +proj=longlat +datum=WGS84 +no_defs

# Load grid shapefile
grids <- st_read(paste0(getwd(), '/Shapefiles/gridcells_110km.shp'))
crs(grids) # +proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs 
names(grids)
grids <- grids[ , c(2,8,11)] # keep cell id 110, wwwf realm, and geometry information

# convert hotspots projection to the same as grids
hotspots_cea <- st_transform(hotspots, crs = crs(grids))
compareCRS(grids, hotspots_cea)

# keep only grids with species on it
grids_cropped <- grids[grids$Cell_Id110 %in% rept_110km$Cell_Id110 , ]

plot(grids_cropped$geometry)
plot(hotspots_cea$geometry, add = T, col = 'red')

# get the number of papers per grid
PapersPerSppByGrid <- rept_110km %>% 
  group_by(Cell_Id110) %>% 
  summarise(TotPapers = sum(N_papers, na.rm = T),
            nSpp = n()) %>%
  mutate(PapersPerSpp = round(TotPapers/nSpp, 2))

# add info 
grids_cropped <- left_join(grids_cropped, PapersPerSppByGrid[,c(1,4)], by = 'Cell_Id110'); rm(PapersPerSppByGrid)

# join both sf_objects based on their intersection
gridHotspotInt <- st_join(grids_cropped, hotspots_cea, join = st_intersects)

# exclude geometry column, then all lines without a hotspot match
summary(gridHotspotInt)
gridHotspotInt <- st_drop_geometry(gridHotspotInt)
gridHotspotInt <- gridHotspotInt[!is.na(gridHotspotInt$NAME) , ]

levels(as.factor(gridHotspotInt$WWF_Realm)) # ok
levels(as.factor(gridHotspotInt$NAME)) # ok

myColors <- c("#6E7EF0","#75FAC1","#BF4140","#DFE03A","#FAAD66","#F065D5","#16FA1C")
names(myColors)<-levels(as.factor(gridHotspotInt$WWF_Realm))

# make a boxplot showing median number of papers per hotspot
(fig5 <- ggplot(gridHotspotInt[!is.na(gridHotspotInt$WWF_Realm), ], 
                 aes(x=reorder(NAME, PapersPerSpp), y=PapersPerSpp, fill=WWF_Realm))+
    geom_boxplot(outlier.size=0.5,  outlier.alpha = 0.7, lwd = 0.4, fatten = 1, colour="black", alpha=0.8)+
    scale_fill_manual(name = "WWF_Realm", values = myColors)+
    labs(x = 'Hotspots', y = 'Number of papers per species')+
    theme_bw()+
    theme(axis.title = element_text(size=10, face="bold"),
          axis.line = element_line(colour="black"),
          axis.text = element_text(size=6, colour = "black"),
          legend.position = c(0.7, 0.2),
          legend.title = element_blank(),
          legend.text = element_text(size = 5),
          legend.background = element_blank())+
    coord_flip())
ggsave(paste0(getwd(), "/Figure5.PapersPerSppByHotspot.pdf"), plot=fig5, width=6, height = 7, units="in", dpi = "print")
ggsave(paste0(getwd(), "/Figure5.PapersPerSppByHotspot.png"), plot=fig5, width=6, height = 7, units="in", dpi = "print")

# Clean global environment
rm(list=ls()); gc()

#####

# STEP 8 - Perform an additional analysis for threat status as assessed vs. non-assessed.
##############################################################################################################
# STEP 8 - Perform an additional analysis for threat status as assessed vs. non-assessed.

# During the modelling of research effort we removed data deficient species (see main text for explanation).
# Now, let's create a new variable distinguishing between species with and
# without a threat status, and then perform basic and applied analysis with this data.

# (Re)load raw data
mydata <- data.table::fread("ResearchEffortData.csv", stringsAsFactors = T, na.strings = "", encoding = 'UTF-8')
n_col<-which(names(mydata)=="Npapers1stApproach")[1] # get the number of the column holding the nº of papers from the 1st approch
colnames(mydata)[n_col] <- 'Npapers'

levels(mydata$ThreatStatus)
mydata$ThreatStatus <- factor(mydata$ThreatStatus,
                              levels = c("LC","NT","VU","EN","CR","EW","EX","DD"))

# Plot the number of publications by threat status
# get summary data by threat status
PapersByThreat <- mydata %>%
  group_by(ThreatStatus) %>%
  summarise(n = n(), 
            mean = round(mean(Npapers, na.rm = T),digits = 2),
            sd = round(sd(Npapers, na.rm = T),digits = 2)) %>%
  mutate(se = sd / sqrt(n))

# Plot 
(figS1 <- ggplot(PapersByThreat, aes(x=ThreatStatus, y=mean))+
  geom_bar(stat = 'identity', alpha=0.8)+
  geom_errorbar( aes(x=ThreatStatus, ymin=mean-se, ymax=mean+se), width=0.2, colour="black", alpha=0.3, size=0.7)+
  labs(x = 'Threat status', y = 'Avg. number of papers per species')+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=10, face="bold"),
        axis.line = element_line(colour="black"),
        axis.text = element_text(size=6, colour = "black"),
        legend.position = 'none'))

ggsave(paste0(getwd(), "/FigureS1.AvgPapersByThreatStatus.pdf"), plot=figS1, width=4, height = 3, units="in", dpi = "print")


### Model research effort using assessed status instead of threat status
mydata <- data.table::fread("ResearchEffortData.csv", stringsAsFactors = T, na.strings = "", encoding = 'UTF-8')
n_col<-which(names(mydata)=="Npapers1stApproach")[1] # get the number of the column holding the nº of papers from the 1st approch
colnames(mydata)[n_col] <- 'Npapers'


# Log-10 transform and standardize (z-tranformation) predictors (mean 0 and SD = 1):
for (i in 1) {
  mydata$StdYear <- scale(mydata$Year, center = T, scale = T)
  mydata$LogBodyMass <- scale(log10(mydata$BodyMass), center = T, scale = T)
  mydata$LogRangeSize <- scale(log10(mydata$RangeSize), center = T, scale = T)
  mydata$LogRangeRarity <- scale(log10(mydata$RangeRarity), center = T, scale = T)
  mydata$LogElevation <- scale(log10(mydata$Elevation+1), center = T, scale = T)
  mydata$LogHumD <- scale(log10(mydata$HumanDensity+1), center = T, scale = T)
  mydata$LogNumBioInst <- scale(log10(mydata$TotBiodInst+1), center = T, scale = T)
  mydata$LogStdPPP <- scale(log10(mydata$Avg_StdPPP), center = T, scale = T)
}

# Check categorical predictors that will be tested:
names(mydata)

### Habit ###
levels(mydata$Habit)
mydata$Verticality <- factor(mydata$Habit,labels = c(0.5,0.75,1,0,0.25,0.5,0.25,0.25,0.5,0.5,0.5,0.75,0.75))
mydata$Verticality <- as.numeric(as.character(mydata$Verticality)) # convert to numeric
mydata$Verticality <- scale(mydata$Verticality, center = T, scale = T) # z-transform variable to make it comparable

### Assessed vs. non-assessed status ###
levels(mydata$ThreatStatus)
mydata$Assessed <- as.factor(ifelse(mydata$ThreatStatus=="DD" | is.na(mydata$ThreatStatus),
                                    yes = 'No', no = 'Yes'))
summary(mydata$Assessed) # okay

### Biog. Realm ###
levels(mydata$MainBiogRealm) # Australasia is missing an 'a'; fixing it next 
mydata$MainBiogRealm <- factor(mydata$MainBiogRealm,
                               labels = c("Afrotropic","Australasia","IndoMalay","Neartic","Neotropic","Oceania","Paleartic"))

# Select response, fixed, random, and grouping variables in a new dataset
names(mydata)
dataset <- mydata%>%
  select(Comb_sciname,Npapers,Family,MainBiogRealm,StdYear:Assessed)
summary(dataset)

# Remove NAs in the new dataset
dataset <- droplevels(dataset[complete.cases(dataset) , ]) # n = 8,269 species
summary(dataset)

# Fit the negative binomial model
# model formula
formula <- Npapers ~ StdYear+LogBodyMass+LogRangeSize+LogRangeRarity+LogElevation+LogHumD+
  LogStdPPP+LogNumBioInst+Verticality+I(Verticality^2)+Assessed+(1|Family)

nbmod <- glmmTMB(formula = formula,
                 data = dataset,
                 family = nbinom2) # assumes the variance increases quadratically with the mean.

# Model validation 
simulateResiduals(fittedModel = nbmod, plot = T, n = 1000) # good fit
performance::r2(nbmod) # R2 (cond) = 0.75, R2 (marg) = 0.59 

summary(nbmod) # model results

rm(list=setdiff(ls(), c("mydata", "dataset", "nbmod")))


### Run a GLMM for each biogeographic realm separately.
names(dataset)

# Extract coefficients from the global model - will be added to the realm-specific models
# to create a tree plot with model coefficients.
res <- summary(nbmod)

# store coefficients
out.global <- as.data.frame(res$coefficients$cond)
out.global$Predictors <- rownames(out.global) # add predictor names into a new column
out.global <- mutate(out.global, 
                     Realm = 'Global', # add realm names into a new column
                     lower95 = Estimate - 1.96 * `Std. Error`, # get lower CI95%
                     upper95 = Estimate + 1.96 * `Std. Error`) # get upper CI95%
out.global <- out.global[ , c("Realm","Predictors","Estimate","Std. Error","z value","Pr(>|z|)","lower95","upper95")]


### FIT NEGATIVE BINOMIAL GLMMs FOR EACH BIOGEOGRAPHIC REALM ###

# Model formula
formula <- Npapers ~ StdYear+LogBodyMass+LogRangeSize+LogRangeRarity+LogElevation+LogHumD+
  LogStdPPP+LogNumBioInst+Verticality+I(Verticality^2)+Assessed+(1|Family)

# Create a vector holdingn biog. realm names to iterate through
levels(dataset$MainBiogRealm)
biorealm <- levels(dataset$MainBiogRealm)
biorealm <- biorealm[-6] # remove Oceania; too small sample size

# Create an empty list to store the outputs:
outs <- list()

for (i in seq_along(biorealm)) {
  # create subdataset containing only observations from a single biog. realm
  subdataset = droplevels(dataset[dataset$MainBiogRealm==biorealm[i] , ])
  
  # fit the NB mixed model
  nbglmm = glmmTMB(formula = formula, data = subdataset, family = nbinom2)
  
  # get global model pseudo R2
  print(performance::r2(nbglmm))
  
  # get result
  res = summary(nbglmm)
  
  # store coefficients
  coefs = as.data.frame(res$coefficients$cond)
  coefs$Predictors <- rownames(coefs) # add predictor names into a new column
  coefs = mutate(coefs, 
                 Realm = paste(biorealm[i]), # add realm names into a new column
                 lower95 = Estimate - 1.96 * `Std. Error`, # get lower CI95%
                 upper95 = Estimate + 1.96 * `Std. Error`) # get upper CI95%
  coefs = coefs[ , c("Realm","Predictors","Estimate","Std. Error","z value","Pr(>|z|)","lower95","upper95")]
  
  outs[[i]] <- coefs
  
  rm(subdataset, nbglmm, coefs)
  
}

results <- bind_rows(outs) # combine results of each biorealm
results <- bind_rows(outs, out.global) # join with results of the global model
results[,c(3:5,7:8)] <- round(results[,c(3:5,7:8)], digits = 3)
results[,6] <- round(results[,6], digits = 4)

rm(list=setdiff(ls(), c("mydata", "dataset", "results")))


### Create a supplementary Tree Plot based on the new threat status
library(ggplot2)
library(viridis)

results <- droplevels(results[results$Predictors!="(Intercept)", ]) # remove intercept
results <- droplevels(results[results$Predictors!="Verticality", ]) # # remove linear verticality (we are interested only in the unimodal relationship; see main text)
results$Predictors <- as.factor(results$Predictors)
levels(results$Predictors) # rename variable names
results$Predictors <- factor(results$Predictors,
                             levels = c("LogBodyMass","I(Verticality^2)","LogRangeSize",
                                        "LogElevation","AssessedYes",
                                        "LogRangeRarity","LogHumD","LogStdPPP",
                                        "LogNumBioInst","StdYear"),  
                             labels = c('Body mass',"Verticality^2","Range size",
                                        "Elevation","Assessed by IUCN","Range rarity",
                                        "Human density","Purchase Power\\nParity",
                                        "Number of biod.\\ninstitutions","Year of description"))

# Reorder levels to plot and add labels:
results$Realm <- as.factor(results$Realm)
levels(results$Realm)
results$Realm <- factor(results$Realm,
                        levels=c("Global","Afrotropic","Australasia","IndoMalay","Neartic",
                                 "Neotropic","Paleartic"))

# Define colors to be used in the plot
myColors <- viridis_pal(option="plasma")(7) 
names(myColors)<-levels(results$Realm)

(fig <-  
    ggplot(results, aes(x = Realm, y = Estimate, ymin = lower95, ymax = upper95))+
    geom_pointrange(aes(col = Realm, shape = Realm), size = 0.1, na.rm = T)+
    scale_colour_manual(name = "Realm", values = myColors)+
    scale_shape_manual(values=c(0,1,2,3,4,5,6))+
    geom_errorbar(aes(ymin=lower95, ymax=upper95, col = Realm), width=0.05, na.rm = T)+
    geom_hline(yintercept =0, linetype=2)+
    labs(x=NULL, y='Average Std. Coefficient (CI95%)')+
    scale_x_discrete(limits = rev(levels(results$Realm)))+
    facet_wrap(~Predictors, strip.position = "left", scales = "free_y", nrow=22)+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=10),
          axis.text = element_text(size=8),
          axis.text.y = element_blank(),
          axis.line = element_line(colour="black"),
          axis.ticks.y.left = element_blank(),
          plot.background=element_rect(fill = "white"),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle=0, hjust = 1, size = 8),
          legend.key = element_blank(),
          legend.key.size = unit(.7,"line"),
          legend.title = element_blank(),
          legend.text=element_text(size=8),
          legend.justification = "center",
          legend.background = element_rect(fill = NA),
          legend.direction = "horizontal",
          legend.position=c(.7,.45))+
    coord_flip())

# Export the figure:
ggsave(paste0(getwd(), "/FigureS7.TreePlot_AssessedStatus.pdf"), plot=fig, width=7, height=5, units="in", dpi = "print")
# This figure was latter edited in the software InkScape for aesthetic adjustments.

rm(list=setdiff(ls(), c()))

#####

# STEP 9 - Perform an additional analysis using polygon-based within-range predictors.
##############################################################################################################
# STEP 9 - Perform an additional analysis using polygon-based within-range predictors.

### Load raw data holding polygon-based predictors
mydata <- fread('polyBasedData.csv', na.strings = "")

### Load grid-based predictors and create a correlogram plot among spatial predictors
gridBased <- data.table::fread("ResearchEffortData.csv", stringsAsFactors = T, na.strings = "", encoding = 'UTF-8')
names(gridBased)
colnames(gridBased)[13] <- 'Npapers'

# Select only spatial-dependent variables that will vary between 
# different approaches (grid- or polygon-based)
gridBased <- gridBased[ , c(4,16:21)]
gridBased$Comb_sciname <- gsub("_", " ", gridBased$Comb_sciname)
gridBased$Comb_sciname <- stringr::str_to_sentence(gridBased$Comb_sciname)
gridBased$Comb_sciname <- stringr::str_trim(gridBased$Comb_sciname)

names(mydata)
polyBased <- mydata[ , c(1,9:14)]

setdiff(gridBased$Comb_sciname, polyBased$Comb_sciname) # rename species
which(gridBased$Comb_sciname=="Gerrhosaurus intermedius") # 4517
gridBased[4517, 'Comb_sciname'] <- 'Gerrhosaurus nigrolineatus'

# join dataframes
corr <- full_join(gridBased, polyBased, by = 'Comb_sciname'); rm(gridBased, polyBased)
# rename variables to facilitate understanding of correlogram plot
names(corr)
colnames(corr)[2] <- "RangeSizeGB"; colnames(corr)[3] <- "RangeRarityGB"
colnames(corr)[4] <- "ElevationGB"; colnames(corr)[5] <- "HumanDensityGB"
colnames(corr)[6] <- "TotBiodInstGB"; colnames(corr)[7] <- "PurchasePowerGB"
colnames(corr)[8] <- "RangeSizePB"; colnames(corr)[9] <- "TotBiodInstPB"
colnames(corr)[10] <- "ElevationPB"; colnames(corr)[11] <- "HumanDensityPB"
colnames(corr)[12] <- "RangeRarityPB"; colnames(corr)[13] <- "PurchasePowerPB"

# reorder columns paring the same variable (but computed through different approaches)
corr <- corr[ , c(1,2,8,3,12,4,10,5,11,6,9,7,13)]

# complete cases only - remove NAs
corr <- corr[complete.cases(corr), ] # 8,448 species

library(corrplot)

# Computing correlation matrix
M <- cor(corr[ , 2:13])

pdf(file=paste0(getwd(), "/FigureS2.CorrelogramPlot.pdf"), width = 7, height = 7)
corrplot(M, method="color", type="upper", addCoef.col = "black", diag = F, 
         tl.col = "black", tl.srt = 45, tl.cex = .8, number.cex = .5)
dev.off()
rm(M,corr)

# Get verticality and ordinal threat status
levels(as.factor(mydata$Habit))
mydata$Verticality <- factor(mydata$Habit,labels = c(0.5,0.75,1,0,0.25,0.5,0.25,0.25,0.5,0.5,0.5,0.75,0.75))
mydata$Verticality <- as.numeric(as.character(mydata$Verticality)) # convert to numeric

levels(as.factor(mydata$ThreatStatus))
mydata$IUCN <- as.integer(as.character(factor(mydata$ThreatStatus,labels = c(5,NA,4,6,NA,1,2,3))))

### Prepare variables that will be used
# Check for skewed distributions and kurtosis of predictors (transform data if necessary).
# Values for skewness and kurtosis between -2 and +2 are considered acceptable
# to prove normal univariate distribution (George & Mallery, 2010) 
e1071::skewness(mydata$Year, na.rm = T); e1071::kurtosis(mydata$Year, na.rm = T) # okay
e1071::skewness(mydata$BodyMass, na.rm = T); e1071::kurtosis(mydata$BodyMass, na.rm = T) 
e1071::skewness(mydata$polykm2, na.rm = T); e1071::kurtosis(mydata$polykm2, na.rm = T) 
e1071::skewness(mydata$range.rarity, na.rm = T); e1071::kurtosis(mydata$range.rarity, na.rm = T) 
e1071::skewness(mydata$elev.avg, na.rm = T); e1071::kurtosis(mydata$elev.avg, na.rm = T)
e1071::skewness(mydata$humd.avg, na.rm = T); e1071::kurtosis(mydata$humd.avg, na.rm = T)
e1071::skewness(mydata$bint_counts, na.rm = T); e1071::kurtosis(mydata$bint_counts, na.rm = T)
e1071::skewness(mydata$W_avgPPP, na.rm = T); e1071::kurtosis(mydata$W_avgPPP, na.rm = T)
# Conclusion: log10 transform all predictors but year of description

# Log-10 transform and standardize (z-tranformation) predictors (mean 0 and SD = 1):
summary(mydata)
# for elevation (with negative values), we'll translate, then transform - [log10(y+1 - min(y))]
for (i in 1) {
  mydata$StdYear <- scale(mydata$Year, center = T, scale = T)
  mydata$LogBodyMass <- scale(log10(mydata$BodyMass), center = T, scale = T)
  mydata$LogRangeSize <- scale(log10(mydata$polykm2), center = T, scale = T)
  mydata$LogRangeRarity <- scale(log10(mydata$range.rarity), center = T, scale = T)
  mydata$LogElevation <- scale(log10(mydata$elev.avg+1 - min(mydata$elev.avg, na.rm = T)), center = T, scale = T)
  mydata$LogHumD <- scale(log10(mydata$humd.avg+1), center = T, scale = T)
  mydata$LogNumBioInst <- scale(log10(mydata$bint_counts+1), center = T, scale = T)
  mydata$LogStdPPP <- scale(log10(mydata$W_avgPPP), center = T, scale = T)
  mydata$Verticality <- scale(mydata$Verticality, center = T, scale = T)
  mydata$IUCN <- scale(mydata$IUCN, center = T, scale = T)
}

# Check multicollinearity among the predictor variables, including squared verticality:
names(mydata)
vars <- mydata[ , Verticality:LogStdPPP]
vars$Verticality2 <- vars$Verticality^2 # add verticality with quadratic term (see main text for explanation)
usdm::vif(vars); rm(vars) # Table S1 in supplementary material
# Conclusion: continuous predictors show low (VIF < 4) multicollinearity (keep them all)

# Biog. Realm - fix small typo in Australasia
levels(as.factor(mydata$MainBiogRealm))
mydata$MainBiogRealm <- factor(mydata$MainBiogRealm,
                               labels = c("Afrotropic","Australasia","IndoMalay","Neartic","Neotropic","Oceania","Paleartic"))

# Select only response, fixed, random, and grouping variables as a new dataset
names(mydata)
dataset <- mydata%>%
  select(Comb_sciname,Npapers,Family,MainBiogRealm,Verticality:LogStdPPP)
summary(dataset)

# Remove NAs in the new dataset
dataset <- droplevels(dataset[complete.cases(dataset) , ]) # n = 7,293 species

dataset %>% group_by(MainBiogRealm) %>% summarise(n = n())
# MainBiogRealm     n
# Afrotropic     1551
# Australasia    1123
# IndoMalay      1312
# Neartic         397
# Neotropic      2318
# Oceania          17
# Paleartic       575

### Run a GLMM for each biogeographic realm separately as well as a global model

# To allow more direct comparisons between grid- and polygon-based models,
# let's also run a negative binomial model for the polygon-based analysis.

# Run the models iterativelly for each biogeographic realm and also for a global model.
# model formula
formula <- Npapers ~ StdYear+LogBodyMass+LogRangeSize+LogRangeRarity+LogElevation+LogHumD+
  LogStdPPP+LogNumBioInst+Verticality+I(Verticality^2)+IUCN+(1|Family)

# Create a vector holdingn biog. realm names to iterate through
levels(dataset$MainBiogRealm)
biorealm <- levels(dataset$MainBiogRealm)
biorealm <- biorealm[-6] # remove Oceania; too small sample size

# Create an empty list to store the outputs:
outs <- list()

library(glmmTMB) # for running glmm

for (i in 1) { # make global model
  # fit the NB mixed model
  nbglmm = glmmTMB(formula = formula, data = dataset, family = nbinom2)
  
  # get global model pseudo R2
  print(paste0('Global', " ", performance::r2(nbglmm)))# cond. R2 = 0.722, marg. R2 = 0.523
  
  # get result
  res = summary(nbglmm)
  
  # store coefficients
  coefs = as.data.frame(res$coefficients$cond)
  coefs$Predictors <- rownames(coefs) # add predictor names into a new column
  coefs = mutate(coefs, 
                 Realm = 'Global', # add realm names into a new column
                 lower95 = Estimate - 1.96 * `Std. Error`, # get lower CI95%
                 upper95 = Estimate + 1.96 * `Std. Error`) # get upper CI95%
  coefs = coefs[ , c("Realm","Predictors","Estimate","Std. Error","z value","Pr(>|z|)","lower95","upper95")]
  
  outs[[i]] <- coefs
  
  rm(nbglmm, res, coefs)
  
  for (j in seq_along(biorealm)) {
    # create subdataset containing only observations from a single biog. realm
    subdataset = droplevels(dataset[dataset$MainBiogRealm==biorealm[j] , ])
    
    # fit the NB mixed model
    nbglmm = glmmTMB(formula = formula, data = subdataset, family = nbinom2)
    
    print(paste0(biorealm[j], " ", performance::r2(nbglmm)))
    
    # get result
    res = summary(nbglmm)
    
    # store coefficients
    coefs = as.data.frame(res$coefficients$cond)
    coefs$Predictors <- rownames(coefs) # add predictor names into a new column
    coefs = mutate(coefs, 
                   Realm = paste(biorealm[j]), # add realm names into a new column
                   lower95 = Estimate - 1.96 * `Std. Error`, # get lower CI95%
                   upper95 = Estimate + 1.96 * `Std. Error`) # get upper CI95%
    coefs = coefs[ , c("Realm","Predictors","Estimate","Std. Error","z value","Pr(>|z|)","lower95","upper95")]
    
    outs[[j+1]] <- coefs
    
    rm(subdataset, nbglmm, res, coefs)
    
  }
}

# Pseudo R2:
# Global: Conditional R2 = 0.736530447895947; Marginal R2 = 0.569220836547044   
# Afrotropic: Conditional R2 = 0.693531580685252; Marginal R2 = 0.502381733316956
# Australasia: Conditional R2 = 0.749618618818637; Marginal R2 = 0.602305492088446 
# IndoMalay: Conditional R2 = 0.761507324903332; Marginal R2 = 0.590326695312412   
# Neartic: Conditional R2 = 0.804407583822889; Marginal R2 = 0.65776972242307  
# Neotropic: Conditional R2 = 0.688354832205891; Marginal R2 = 0.519926499170894   
# Paleartic: Conditional R2 = 0.782425433324134; Marginal R2 = 0.631646007113889 


results <- bind_rows(outs) # combine results
names(results)
results[,c(3:5,7:8)] <- round(results[,c(3:5,7:8)], digits = 3)
results[,6] <- round(results[,6], digits = 4)

# save results
fwrite(results, file = 'polyBasedModelResults.csv')

rm(list=setdiff(ls(), c()))

#####

# Step 10. Create a bubble plot for grid- and polygon-based model coefficients
##############################################################################################################
# Step 10. Create a bubble plot for grid- and polygon-based model coefficients

# load grid-model results
gridBasedOuts <- fread('gridBasedModelResults.csv')
# create new column informing approach used for computing spatial variables
gridBasedOuts$Approach <- 'grid-based'

# load polygon-model results
polyBasedOuts <- fread('polyBasedModelResults.csv')
polyBasedOuts$Approach <- 'polygon-based'

# merge datasets
outs <- as.data.frame(rbind(gridBasedOuts, polyBasedOuts)); rm(gridBasedOuts, polyBasedOuts)

# make a bubble plot based on coefficient values 
names(outs)

# create column informing whether a given predictor was significant or not;
# this can be use to fill or not the bubbles
outs$sig <- NA
for(i in 1:nrow(outs)){
  lower <- as.numeric(outs[i, 'lower95'])
  upper <- as.numeric(outs[i, 'upper95'])
  
  if(lower < 0 & upper < 0 || lower > 0 & upper > 0){
    outs[i, 'sig'] <- TRUE # for significant predictors
  } else{
    outs[i, 'sig'] <- FALSE # for non-significant predictors
  }
  rm(lower, upper)
}

# create column informing the direction of the coefficient (negative or positive)
outs$dir <- NA
for(i in 1:nrow(outs)){
  if(outs[i, 'Estimate'] < 0 ){
    outs[i, 'dir'] <- FALSE # for negative effects
  } else{
    outs[i, 'dir'] <- TRUE # for positive effects
  }
}

# Reorder realms
levels(as.factor(outs$Realm))
outs$Realm <- factor(outs$Realm, 
                     levels = c("Global","Afrotropic","Australasia","IndoMalay", 
                                "Neartic","Neotropic","Paleartic"))

# Remove intercept and linear verticality
outs <- outs[!outs$Predictors=="(Intercept)", ] # remove intercept
outs <- outs[!outs$Predictors=="Verticality", ] # remove intercept

# Reorder predictors  
levels(as.factor(outs$Predictors))
outs$Predictors <- factor(outs$Predictors, 
                          levels = c("LogBodyMass","I(Verticality^2)","LogRangeSize",
                                     "LogElevation","IUCN","LogRangeRarity","LogHumD",
                                     "LogStdPPP","LogNumBioInst","StdYear"),
                          labels = c("Body mass","Verticality","Range size","Elevation", 
                                     "Threat status","Range rarity","Human density",
                                     "Purchase Power\\nParity","Number of biod.\\ninstitutions",
                                     "Year of description"))

# Define colors to be used (https://www.canva.com/colors/color-wheel/)
myColors <- c("#2478DB", "#E4461B") # blue and red
names(myColors)<-levels(as.factor(outs$dir))

summary(abs(outs$Estimate))

(p <- ggplot(outs, aes(x = Realm, y = Estimate,
                       size = abs(Estimate), color = dir))+
    geom_point(aes(shape=as.factor(sig)))+ 
    scale_color_manual(name = "dir", values = myColors, guide = 'none')+
    scale_shape_manual(values = c(1, 16), guide='none')+ # open vs. closed circle
    scale_size(name = '\\u03B2 GLMM')+
    labs(x=NULL, y=NULL)+
    theme_classic()+
    facet_grid(Predictors~Approach, scales = 'free_y')+
    theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
          strip.text.x = element_text(size = 7),
          strip.text.y = element_text(angle = 0, size = 7),
          panel.spacing = unit(.2, "lines")))

# save figure as pdf
ggsave(paste0(getwd(), "/Figure5.bubble_plot.pdf"), plot=p, width=7, height=15, units="in", dpi = "print", device = cairo_pdf)
# this figure was edited in InkScape to improve aesthetics.

rm(list=setdiff(ls(), c()))


##### END OF THE SCRIPT #####