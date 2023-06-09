### Supporting Information to ###

# Title: Unwrapping broken tails: biological and environmental correlates of predation pressure in limbless reptiles
# Authors: Moura, MR., HC Costa, AD Abegg, E Alaminos, T Angarita-Sierra, W Azevedo, H Cabral, P Carvalho, S Cechin, N Citeli,
#                      ACM Dourado, AFV Duarte, EMX Freire, PCA Garcia, R Montero, ACM Moraes-da-Silva, R Mol, JM Pleguezuelos,
#                      RFD Sales, DJ Santana, L Santos, VTC Silva, V Sudré, DC Passos, P Passos, R Perez, P Prado, A Prudente,
#                      O Torres-Carvajal, JJ Torres-Ramirez, V Wallach, GR Winck, and JJM Guedes
# Journal: Journal of Animal Ecology
# * Corresponding author: mariormoura@gmail.com
# https://mariormoura.wordpress.com/

# STEPS IN THIS SCRIPT 
#  1. Load the dataset and standardise covariates.
#  2. Build histograms of body size across other species-level predictors.
#  3. Plot the examined species in the geographical space.
#  4. Build the barplots of frequency of tail loss across levels of biological attributes.
#  5. Plot autotomy frequency per species across the phylogeny.
#  6. Analyse determinants of autotomy probability using Generalised Linear Mixed Effect Model.
#  7. Get predicted probabilities of tail loss across values of continuous covariates.
#  8. Check if model residuals show phylogenetic or spatial autocorrelation.
#  9. Inspect sex- and age-dependent effects on per species autotomy frequency.
#  10. Make boxplots to check sexual size dimorphism among species.

# First, clean workspace and install recquired packages:
rm(list=ls()); gc()
setwd("DefineYourDirectory")

# Install and load needed packages:
needed_packages<-c("ape", "boot", "car", "compiler", "CoordinateCleaner", "data.table", "devtools", "dplyr", "geiger",
                   "GGally", "ggnewscale", "ggtree", "ggplot2",  "ggpubr", "lattice", "lme4", "lmtest", "maptools", 
                   "parallel", "phangorn", "phylobase", "phylosignal", "phytools", "pgirmess", "plyr", "readr",
                   "rgdal", "raster", "reshape2", "sf", "SoDA", "sp", "SpatialPack", "stringr", "tidyr", "tools")
new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)
install_github("kassambara/easyGgplot2", force = T)

#####

# STEP 1 - Load the dataset and standardise covariates
##############################################################################################################
# STEP 1 - Load the dataset and standardise covariates
rm(list=ls())# Load the dataset and understand it:

Data<-fread("CompiledDataCleanned.txt", stringsAsFactors=TRUE, encoding="UTF-8")
names(Data)
# SpeciesName_ReptileDatabase: valid species name, as informed in the ReptileDatabase.
# Genus: taxonomic genus.
# Family: taxonomic family.
# Suborder: taxonomic suborder.
# SpecimenID: collection label for the respective specimen.
# CollectionAcronym: collection acronym.
# Collaborator: person who examined the specimen.
# Locality: local name for the collection site.
# Municipality: city in which the specimen was collected.
# AdmUnit: state, province, or department in which the specimen was collected.
# Country: country in which the specimen was collected.
# Latitude: latitude in decimal degrees.
# Longitude: longitude in decimal degrees.
# AutotomyPresence: binary, 0 = tail intact, 1 = tail autotomised.
# SVL_mm: snout-vent length (mm).
# Sex: factor (male, female, or unknown).
# LifeStage: factor (adult, juvenile)
# Fossorial: species habitat use, binary, 0 = fossorial habitat not used, 1 = fossorial habitat used.
# Aquatic: species habitat use, binary, 0 = aquatic habitat not used, 1 = aquatic habitat used.
# Terrestrial: species habitat use, binary, 0 = terrestrial habitat not used, 1 = terrestrial habitat used.
# Arboreal: species habitat use, binary, 0 = arboreal habitat not used, 1 = arboreal habitat used.
# Activity: activity time pattern, factor (diurnal, nocturnal, cathemeral).
# AMT: annual mean temperature at the collection site, in Celsius (extracted from WorldClim v.2, https://www.worldclim.org)
# AMP: average monthly precipitation at the collection site, in mm (extracted from WorldClim v.2, https://www.worldclim.org).
# Biome: biome where the specimen was collected (extracted from https://ecoregions.appspot.com).
# Tropicality: major climatic zone of the collection site, binary, 0 = temperate, 1 = tropical.

# Compute the number of specimens and autotomy events per species: 
Data$AutotomyPresence<-as.integer(as.character(Data$AutotomyPresence))
Autotomy_N<-Data %>% 
  dplyr::group_by(SpeciesName_ReptileDatabase) %>%
  dplyr::summarise(TotAutotomy=sum(AutotomyPresence, na.rm = T),
                   TotOccurrences=length(SpeciesName_ReptileDatabase))

# Merge species-level information on number of specimens examined and autotomy events:
Data<-merge(Data, Autotomy_N, by="SpeciesName_ReptileDatabase", all.X=T)

# Check the total number of species and specimens in the dataset:
Data %>% 
  dplyr::group_by(Suborder) %>% 
  dplyr::summarise(Records=length(SpeciesName_ReptileDatabase),
                   Richness=length(unique(SpeciesName_ReptileDatabase)))

# Set the minimum sampling effort and filter the dataset:
Data<-Data[Data$TotAutotomy>=5,] # at least five autotomy events per species
Data<-Data[Data$TotOccurrences>=30,] # at least 30 examined specimens per species
Data<-droplevels(Data)
Data %>% 
  dplyr::group_by(Suborder) %>% 
  dplyr::summarise(Records=length(SpeciesName_ReptileDatabase),
                   Richness=length(unique(SpeciesName_ReptileDatabase)))

# Convert LifeStage to binary:
summary(Data$LifeStage)
Data$LifeStageBinary<-factor(Data$LifeStage,
                             levels=c("Juvenile", "Adult"),
                             labels=c(0,1))
Data$LifeStageBinary<-as.numeric(as.character(Data$LifeStageBinary))

# Convert Sex to ordinal:
summary(Data$Sex)
Data[Data$SpecimenID=="MNHN 1992.4699",]$Sex<-"Unknown" # single female juvenile Cynisca leucura converted to unknown sex to avoid discarding the specimen
Data$SexBinary<-factor(Data$Sex,
                       levels=c("Male", "Unknown", "Female"),
                       labels=c(0,0.5,1))
Data$SexBinary<-as.numeric(as.character(Data$SexBinary))
Data %>% 
  dplyr::group_by(Suborder, Sex) %>% 
  dplyr::summarise(Records=length(SpeciesName_ReptileDatabase),
                   Richness=length(unique(SpeciesName_ReptileDatabase)))

# Convert Tropicality to binary:
summary(Data$Tropicality)
Data$Tropicality<-factor(Data$Tropicality,
                         levels = c("Temperate", "Tropical"),
                         labels = c(0, 1))
Data$Tropicality<-as.numeric(as.character(Data$Tropicality))
summary(Data$Tropicality)

# Convert Activity pattern to ordinal, but first set Cynisca leucura as cathemeral to avoid the removal of this species:
Data[Data$Activity=="" & Data$SpeciesName_ReptileDatabase=="Cynisca leucura",]$Activity<-"cathemeral" 
Data<-droplevels(Data)
Data$Diurnality<-factor(Data$Activity,
                        levels = c("nocturnal", "cathemeral", "diurnal"),
                        labels = c(0, 0.5, 1))
Data$Diurnality<-as.numeric(as.character(Data$Diurnality))

# Convert HabitatUse to ordinal:
Data<-as.data.frame(Data)
Data$Verticality<-as.character(NA)
Data[which(Data$Fossorial==1 & Data$Terrestrial==0 & Data$Aquatic==0 & Data$Arboreal==0),]$Verticality<-0 # Fossorial
#Data[which(Data$Fossorial==1 & Data$Terrestrial==1 & Data$Aquatic==0 & Data$Arboreal==0),]$Verticality<-0.25 # Semifossorial, condition not present in this data
Data[which(Data$Fossorial==0 & Data$Terrestrial==1 & Data$Aquatic==0 & Data$Arboreal==0),]$Verticality<-0.5 # Terrestrial (Terrestrial)
Data[which(Data$Fossorial==0 & Data$Terrestrial==0 & Data$Aquatic==1 & Data$Arboreal==0),]$Verticality<-0.5 # Aquatic (Terrestrial)
Data[which(Data$Fossorial==0 & Data$Terrestrial==1 & Data$Aquatic==1 & Data$Arboreal==0),]$Verticality<-0.5 # Semiaquatic (Terrestrial)
Data[which(Data$Fossorial==0 & Data$Terrestrial==1 & Data$Aquatic==0 & Data$Arboreal==1),]$Verticality<-0.75 # Semiarboreal
#Data[which(Data$Fossorial==0 & Data$Terrestrial==1 & Data$Aquatic==1 & Data$Arboreal==1),]$Verticality<-0.75 # Semiarboreal, condition not present in this data
#Data[which(Data$Fossorial==0 & Data$Terrestrial==0 & Data$Aquatic==1 & Data$Arboreal==1),]$Verticality<-0.75 # Semiarboreal, condition not present in this data
Data[which(Data$Fossorial==0 & Data$Terrestrial==0 & Data$Aquatic==0 & Data$Arboreal==1),]$Verticality<-1 # Arboreal
summary(as.factor(Data$Verticality))

# Create a function to standardise predictors in the range of 0-1:
range01 <- function(x){(x-min(x,na.rm=T))/(max(x, na.rm=T)-min(x,na.rm=T))}

# Standardise body size separately for each species, sex, and life-stage:
Data$SpeciesLifeStage<-paste0(Data$SpeciesName_ReptileDatabase, "_", Data$LifeStage, "_", Data$Sex)

# How many samples in each 'SpeciesLifeStage' category?
SummaryTable <- Data %>% dplyr::group_by(SpeciesLifeStage) %>% dplyr::summarise(N = length(SpeciesName_ReptileDatabase))
summary(SummaryTable$N)

# Select 'SpeciesLifeStage' combinations with at least five samples (exclude species-sex-life-stage categories with a single specimen):
Data<-droplevels(Data[Data$SpeciesLifeStage %in% SummaryTable[SummaryTable$N>=5,]$SpeciesLifeStage, ])
Data<-as.data.table(Data)
Data<-Data[ , group_StdSVL := range01(SVL_mm), "SpeciesLifeStage"]
Data$StdSVL<-range01(log10(Data$SVL_mm)) # Temp is based on collection month whereas AMT is the annual average

# Standardise environmental variables in the range of 0-1:
Data$StdTemp<-range01(log10(Data$AMT)) 
Data$StdPrec<-range01(log10(Data$AMP+1))

# Re-extract the number of records per species:
Data$AutotomyPresence<-as.integer(as.character(Data$AutotomyPresence))
Autotomy_N<-Data %>% 
  dplyr::group_by(SpeciesName_ReptileDatabase) %>%
  dplyr::summarise(TotAutotomy=sum(AutotomyPresence, na.rm = T),
                   TotOccurrences=length(SpeciesName_ReptileDatabase))

# Use only species with at least five autotomy records and 30 examined specimens:
Data<-merge(Data[,-c("TotOccurrences", "TotAutotomy")], Autotomy_N, by="SpeciesName_ReptileDatabase", all.X=T)
Data<-Data[Data$TotAutotomy>=5,] # at least five autotomy events per species
Data<-Data[Data$TotOccurrences>=30,] # at least 30 examined specimens per species
Data<-droplevels(Data)
Data %>% 
  dplyr::group_by(Suborder) %>% 
  dplyr::summarise(Records=length(SpeciesName_ReptileDatabase),
                   Richness=length(unique(SpeciesName_ReptileDatabase)))

# Inspect the multicollinearity among predictors:
PredictorData<-Data[, c("LifeStageBinary", "group_StdSVL", "SexBinary", "Verticality", "Diurnality", "StdTemp", "StdPrec", "Tropicality")]
usdm::vif(PredictorData) # using body size rescaled within each species-sex-lifestage

PredictorData<-Data[, c("LifeStageBinary", "StdSVL", "SexBinary", "Verticality", "Diurnality", "StdTemp", "StdPrec", "Tropicality")]
usdm::vif(PredictorData) # using body size rescaled across all species

#####

# STEP 2 - Build histograms of body size across other species-level predictors
##############################################################################################################
# STEP 2 - Build histograms of body size across other species-level predictors
rm(list=setdiff(ls(), c("Data")))

# Histogram of body size across other species-level predictors:
Data$fVerticality<-factor(Data$Verticality, levels=c("0", "0.5", "0.75", "1"), labels=c("Low", "Medium", "High", "Very high"))
myColours<-grDevices::hcl.colors(n = nlevels(Data$fVerticality), palette="viridis", rev=F)
names(myColours)<-levels(Data$fVerticality)
MyHist1<-easyGgplot2::ggplot2.histogram(data=Data, xName='SVL_mm', groupName='fVerticality',
                                        addMeanLine=TRUE,
                                        meanLineSize=1,
                                        meanLineColor=myColours,
                                        groupColours=myColours,
                                        alpha=0.1,
                                        addDensityCurve = TRUE,
                                        removePanelGrid=TRUE,
                                        removePanelBorder=TRUE,
                                        axisLine=c(0.5, "solid", "black"),
                                        showLegend=TRUE,
                                        binwidth=0.01,
                                        scale="density",
                                        backgroundColor="white") +
  labs(x="Snout-vent length", y="Density") +
  theme(legend.position=c(0.75, 0.75),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.background=element_blank(),
        panel.spacing = unit(0, "lines"), 
        axis.line = element_line(colour="black", size=0.4),
        axis.ticks = element_line(colour="black", size=0.4),
        axis.title = element_text(size=14, colour="black"),
        axis.text.y.left = element_text(size=12, colour="black", face="plain"),
        axis.text.x.bottom =  element_text(size=12, colour="black", face="plain"))

Data$fDiurnality<-factor(Data$Diurnality, levels=c("0", "0.5", "1"), labels=c("Low", "Medium", "High"))
myColours<-grDevices::hcl.colours(n = nlevels(Data$fDiurnality), palette="viridis", rev=F)
names(myColours)<-levels(Data$fDiurnality)
MyHist2<-easyGgplot2::ggplot2.histogram(data=Data, xName='SVL_mm', groupName='fDiurnality',
                                        addMeanLine=TRUE,
                                        meanLineSize=1,
                                        meanLineColor=myColours,
                                        groupColours=myColours,
                                        alpha=0.1,
                                        addDensityCurve = TRUE,
                                        removePanelGrid=TRUE,
                                        removePanelBorder=TRUE,
                                        axisLine=c(0.5, "solid", "black"),
                                        showLegend=TRUE,
                                        binwidth=0.01,
                                        scale="density",
                                        backgroundColor="white") +
  labs(x="Snout-vent length", y="Density") +
  theme(legend.position=c(0.75, 0.75),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.background=element_blank(),
        panel.spacing = unit(0, "lines"), 
        axis.line = element_line(colour="black", size=0.4),
        axis.ticks = element_line(colour="black", size=0.4),
        axis.title = element_text(size=14, colour="black"),
        axis.text.y.left = element_text(size=12, colour="black", face="plain"),
        axis.text.x.bottom =  element_text(size=12, colour="black", face="plain"))

Multipanel_plot<-ggpubr::ggarrange(MyHist1, MyHist2, labels=c("a)", "b)"), align="hv",
                                   font.label=list(size=12, colour = "black"), ncol=2, nrow=1)

# Save to disk:
ggsave("FigureS1.png", plot=Multipanel_plot, width=10, height=4, units="in", bg = "transparent")
ggsave("FigureS1.pdf", plot=Multipanel_plot, width=10, height=4, units="in", bg = "transparent")

#####

# STEP 3 - Plot the examined species in the geographical space
##############################################################################################################
# STEP 3 - Plot the examined species in the geographical space
rm(list=setdiff(ls(), c("Data")))

# Duplicate data selecting only variables need for plotting the map:
Data<-Data[order(Data$AutotomyPresence, decreasing=F), ] # To plot autotomy events on top
Data<-Data[order(Data$Suborder, decreasing=T), ] # To plot amphisbaenian occurrences on top

# Load shapefile on WWF Ecoregions, available at: https://storage.googleapis.com/teow2016/Ecoregions2017.zip
wgs84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # set the coordinate reference system parameters
wwf_biomes<-readOGR(dsn="DefineYourDirectory", layer='Ecoregions2017') # change directory as needed
wwf_biomes<-spTransform(wwf_biomes, wgs84) # double-check if wwf_biomes are in the same crs as the specimen database
wwf_biomes<-sf::st_as_sf(wwf_biomes) # convert to sf object (faster to plot)

# Load the spatial data on country borders:
gadm_countries<-rnaturalearth::ne_countries(scale = 10, type = "countries")
gadm_countries<-sf::st_as_sf(gadm_countries) # convert to sf object

# Relabel levels of 'BIOME_NAME' in Ecoregion2017 shapefile:
unique(levels(as.factor(wwf_biomes$BIOME_NAME)))
wwf_biomes$BiomeFill<-factor(wwf_biomes$BIOME_NAME,
                             levels=c("Montane Grasslands & Shrublands",
                                      "Tundra",
                                      "Boreal Forests/Taiga",
                                      "Tropical & Subtropical Moist Broadleaf Forests",   
                                      "Mediterranean Forests, Woodlands & Scrub",
                                      "Mangroves",
                                      "Temperate Broadleaf & Mixed Forests",
                                      "Temperate Conifer Forests",
                                      "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                      "Tropical & Subtropical Coniferous Forests",
                                      "Temperate Grasslands, Savannas & Shrublands",
                                      "Tropical & Subtropical Dry Broadleaf Forests",
                                      "Deserts & Xeric Shrublands",
                                      "Flooded Grasslands & Savannas",
                                      "N/A"),
                             labels=c("Montane grass/shrub",
                                      "Tundra",
                                      "Boreal Forests/Taiga",
                                      "Trop/subtrop moist broadleaf",          
                                      "Mediterranean",
                                      "Mangroves",
                                      "Temperate broadleaf/mixed",
                                      "Temperate coniferous",
                                      "Trop/subtrop grass/savan/shrub",
                                      "Trop/subtrop coniferous",
                                      "Temperate grass/savan/shrub",
                                      "Trop/subtrop dry broadleaf",
                                      "Deserts/xeric shrub",
                                      "Flooded grass/savan",
                                      ""))

# Get a colour ramp for biome types:
colourCount<-nlevels(as.factor(wwf_biomes$BiomeFill)) # double the number of colours
colour_set<-RColorBrewer::brewer.pal(2*colourCount, "Greys")[1:7]  # use only the lighter-half of the white-to-black colour ramp
getPalette = colorRampPalette(colour_set) # interpolate hexadecimal colours 
biomes_fill<-scale_fill_manual(values=rev(getPalette(colourCount)))

# Plot the map:
MyMap<-ggplot(data=Data) + 
  
  # Add polygon boundaries for biomes:
  geom_sf(data=wwf_biomes, aes(fill=BiomeFill), colour=NA) +  biomes_fill +
  
  # Add polygon boundaries for the country:
  geom_sf(data=gadm_countries, fill=NA, colour="black", size=0.15) +
  
  # Restrict the map to study area extent and define projection:
  coord_sf(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
           xlim = c((min(Data$Longitude))-1, (max(Data$Longitude))+1),
           ylim = c((min(Data$Latitude))-1, (max(Data$Latitude))+1), expand = FALSE) +
  
  # Other aesthetics
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.position="none",
        panel.background=element_rect(fill="#d1e5f0", colour=NA),
        plot.background=element_rect(fill="transparent", colour=NA),
        panel.border=element_rect(fill=NA, colour="black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.placement = "none")

# Add a new colour scale for the fill aesthetics:
MyMap2 <- MyMap + ggnewscale::new_scale_fill()

# Prepare variables to guide aesthetics of occurrence records:
Data$TailCondition<-factor(Data$AutotomyPresence, levels=c(0,1), labels=c("Without autotomy", "With autotomy"))
myColours <- c("#67001f", "#053061", "#d6604d", "#4393c3")
Data$fill_var<-paste0(Data$TailCondition, " ", Data$Suborder)
names(myColours)<-levels(as.factor(Data$fill_var))

# Built the final map (add the occurrence records):
MyMap2 <- MyMap2 +
  geom_point(data=Data, aes(x=Longitude, y=Latitude,
                            shape=TailCondition, size=TailCondition,
                            colour=fill_var, fill=fill_var), alpha=0.5) +
  scale_fill_manual(values = myColours) +
  scale_colour_manual(values = myColours) +
  scale_size_manual(values=c(1.5, 1)) +
  scale_shape_manual(values=c(4, 21))

# Save to disk:
ggsave("Figure1.png", plot=MyMap2, width=11, height=6, units="in", bg = "transparent")
ggsave("Figure1.pdf", plot=MyMap2, width=11, height=6, units="in", bg = "transparent")
# OBS: minor editions performed outside R, in InkScape.

#####

# STEP 4 - Build the barplots of frequency of tail loss across levels of biological attributes 
##############################################################################################################
# STEP 4 - Build the barplots of examined species across levels of biological attributes 
rm(list=setdiff(ls(), c("Data")))

### BARPLOT 1 - Frequency of autotomy across examined species:
summary_table<-as.data.frame(Data %>% dplyr::group_by(Suborder, SpeciesName_ReptileDatabase) %>% 
                               dplyr::summarise(N_Records=count(AutotomyPresence),
                                                N_Total=length(SpeciesName_ReptileDatabase))) %>% ungroup()
summary_table<-data.frame("Suborder"=summary_table$Suborder,
                          "Category"=summary_table$SpeciesName_ReptileDatabase,
                          "Autotomy"=summary_table$N_Records$x,
                          "N_Records"=summary_table$N_Records$freq,
                          "N_Total"=summary_table$N_Total)
summary_table$PropRecords<-summary_table$N_Records/summary_table$N_Total

# Create a column to inform the number of examined species per species:
summary_table$N_Label<-NA
for(i in 1:nlevels(summary_table$Category)){
  
  selected_row<-which(summary_table$Category==levels(summary_table$Category)[i])
  
  if(length(selected_row)<=2){
    
    summary_table$N_Label[selected_row[1]]<-paste(summary_table$N_Total[selected_row[1]])
    
  } else {
    
    midrow<-ceiling(mean(selected_row))
    summary_table$N_Label[midrow]<-paste(summary_table$N_Total[midrow])
    
  }
}

# Relabel tail condition:
summary_table$TailCondition<-factor(summary_table$Autotomy,
                                    levels=c(0,1),
                                    labels=c("Scar-free tail", "Scarred tail"))

# Reorder species in increasing order of autotomy frequency:
summary_table<-summary_table[order(summary_table$PropRecords, decreasing=F),]
summary_table<-summary_table[order(summary_table$Autotomy, decreasing=T),]
summary_table$CategoryLabel<-factor(summary_table$Category, levels=rev(unique(summary_table$Category)))

# Get the average, max, and min autotomy frequency across all species and within each Suborder:
summary_table %>% dplyr::group_by(Autotomy) %>% dplyr::summarize(Count=sum(N_Records),
                                                                 Avg=mean(PropRecords, na.rm=T),
                                                                 Min=min(PropRecords, na.rm=T),
                                                                 Max=max(PropRecords, na.rm=T))
summary_table %>% dplyr::group_by(Suborder, Autotomy) %>% dplyr::summarize(Avg=mean(PropRecords, na.rm=T),
                                                                           Min=min(PropRecords, na.rm=T),
                                                                           Max=max(PropRecords, na.rm=T))

# Prepare colours:
summary_table$fill_var<-paste0(summary_table$TailCondition, " ", summary_table$Suborder)
myColours <- c("#fddbc7", "#d1e5f0", "#f4a582", "#92c5de")
names(myColours)<-levels(summary_table$fill_var)
fillScale<-scale_fill_manual(values = myColours)

# Build the plot:
MyBarplot<-list()
MyBarplot[[1]]<-ggplot(data=summary_table, aes(x=CategoryLabel, y=PropRecords, fill=fill_var)) +
  geom_bar(position="stack", stat="identity") +
  geom_text(y=1.025, aes(label=N_Label), size=2, colour="gray25", position=position_stack(vjust=0.5)) +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Proportion of specimens with tail autotomy") +
  fillScale +
  theme(panel.grid.minor = element_blank(), # remove minor grid lines
        panel.grid.major = element_blank(), # remove major grid lines
        panel.background = element_blank(), # white background
        axis.line = element_line(colour="black"), # axis lines aesthetics
        axis.text.y = element_text(hjust=1, vjust=0.5, angle=0, size=8, face="italic"),
        axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=8, margin=margin(t=0, r=0, b=5, l=0)),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=10, colour="black", face="bold", margin=margin(t=0, r=5, b=0, l=0)), # margin between axis.title and axis.values
        legend.position="none",
        plot.background=element_rect(fill="transparent", colour=NA),
        plot.margin=unit(c(0, 0.65, 0.05, 0.3), "cm"), # top, right, bottom, left
        panel.spacing = unit(0.1, "lines")) +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1), 
                     expand=expansion(mult=c(0, .05), add=c(0, 0)),
                     labels=scales::number_format(accuracy = 0.01, decimal.mark='.')) +
  scale_x_discrete(position="bottom")


### BARPLOT 2 - Frequency of autotomy across levels of life-stage and sex:
# Get the marginal totals per life-stage category:
summary_table<-as.data.frame(Data %>% dplyr::group_by(Suborder, LifeStage) %>% 
                               dplyr::summarise(N_Records=count(AutotomyPresence),
                                                N_Total=length(SpeciesName_ReptileDatabase))) %>% ungroup()
summary_table1<-data.frame("Suborder"=summary_table$Suborder,
                           "Category"=summary_table$LifeStage,
                           "Autotomy"=summary_table$N_Records$x,
                           "N_Records"=summary_table$N_Records$freq,
                           "N_Total"=summary_table$N_Total)

# Get the marginal totals per sex category:
summary(Data$Sex)
summary_table<-as.data.frame(Data %>% dplyr::filter(Sex!="Unknown") %>%  dplyr::group_by(Suborder, Sex) %>% 
                               dplyr::summarise(N_Records=count(AutotomyPresence),
                                                N_Total=length(SpeciesName_ReptileDatabase))) %>% ungroup()
summary_table2<-data.frame("Suborder"=summary_table$Suborder,
                           "Category"=summary_table$Sex,
                           "Autotomy"=summary_table$N_Records$x,
                           "N_Records"=summary_table$N_Records$freq,
                           "N_Total"=summary_table$N_Total)

# Get the marginal totals per life-stage and sex category combined:
Data$Category<-paste0(Data$LifeStage, "_", Data$Sex)
summary_table<-as.data.frame(Data %>% dplyr::filter(Sex!="Unknown") %>% dplyr::group_by(Suborder, Category) %>% 
                               dplyr::summarise(N_Records=count(AutotomyPresence),
                                                N_Total=length(SpeciesName_ReptileDatabase)))
summary_table3<-data.frame("Suborder"=summary_table$Suborder,
                           "Category"=summary_table$Category ,
                           "Autotomy"=summary_table$N_Records$x,
                           "N_Records"=summary_table$N_Records$freq,
                           "N_Total"=summary_table$N_Total)

# Bind all summary tables and get the proportion of records in each level of autotomy:
summary_table<-rbind(summary_table1, summary_table2, summary_table3)
rm(summary_table1, summary_table2, summary_table3)
summary_table$PropRecords<-summary_table$N_Records/summary_table$N_Total

# Relabel tail condition:
summary_table$TailCondition<-factor(summary_table$Autotomy, levels=c(0,1), labels=c("Scar-free tail", "Scarred tail"))

# Specify the order of appearance of each life-stage and sex category in the barplot:
summary_table$CategoryLabel<-paste0(summary_table$Suborder, "_", summary_table$Category)
levels(as.factor(summary_table$CategoryLabel))
summary_table$CategoryLabel<-factor(summary_table$CategoryLabel,
                                    levels=rev(c("Amphisbaenia_Adult", "Amphisbaenia_Juvenile", "Serpentes_Adult", "Serpentes_Juvenile", 
                                                 "Amphisbaenia_Female", "Amphisbaenia_Male", "Serpentes_Female", "Serpentes_Male",
                                                 
                                                 "Amphisbaenia_Adult_Female", "Amphisbaenia_Adult_Male", 
                                                 "Amphisbaenia_Juvenile_Female", "Amphisbaenia_Juvenile_Male",
                                                 
                                                 "Serpentes_Adult_Female", "Serpentes_Adult_Male", 
                                                 "Serpentes_Juvenile_Female", "Serpentes_Juvenile_Male")),
                                    
                                    labels=rev(c("Adult amphisb.", "Juvenile amphisb.", "Adult snakes", "Juvenile snakes", 
                                                 "Female amphisb.", "Male amphisb.", "Female snakes", "Male snakes",
                                                 
                                                 "Ad. fem. amphisb.", "Ad. male amphisb.",
                                                 "Juv. fem. amphisb.", "Juv. male amphisb.",
                                                 
                                                 "Ad. female snakes", "Ad. male snakes",
                                                 "Juv. female snakes", "Juv. male snakes")))

# Create a column to inform the number of examined species per category of life-stage and sex:
summary_table$N_Label<-NA
for(i in 1:nlevels(summary_table$CategoryLabel)){
  
  selected_row<-which(summary_table$CategoryLabel==levels(summary_table$CategoryLabel)[i])
  
  if(length(selected_row)<=2){
    
    summary_table$N_Label[selected_row[1]]<-paste(summary_table$N_Total[selected_row[1]])
    
  } else {
    
    midrow<-ceiling(mean(selected_row))
    summary_table$N_Label[midrow]<-paste(summary_table$N_Total[midrow])
    
  }
}

# Specify the colours for the barplot:
summary_table$fill_var<-paste0(summary_table$TailCondition, "_", summary_table$Suborder)
myColours <- c("#fddbc7", "#d1e5f0", "#f4a582", "#92c5de") # amph_not_uro #snake_not_uro #amph_uro # snake_uro
names(myColours)<-levels(summary_table$fill_var)
fillScale<-scale_fill_manual(values = myColours)

# Build the plot:
MyBarplot[[2]]<-ggplot(data=summary_table, aes(x=CategoryLabel, y=PropRecords, fill=fill_var)) +
  geom_bar(position="stack", stat="identity") +
  geom_text(y=1.025, aes(label=N_Label), size=2, colour="gray25", position=position_stack(vjust=0.5)) +
  coord_flip() +
  xlab("") + ylab("") +
  fillScale +
  theme(panel.grid.minor = element_blank(), # remove minor gridlines
        panel.grid.major = element_blank(), # remove major gridlines
        panel.background = element_blank(), # white background
        axis.line = element_line(colour="black"), # axis lines aesthetics
        axis.text.y = element_text(hjust=1, vjust=0.5, angle=0, size=8),
        axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=8),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position="none",
        plot.background=element_rect(fill="transparent", colour=NA),
        plot.margin=unit(c(0, 0.65, 0.05, 0), "cm"), # top right bottom left
        panel.spacing = unit(0.1, "lines")) +
  geom_vline(xintercept=8.5, color= "gray50", size=0.5) +  # add a dotted line at x=1 after flip
  geom_vline(xintercept=12.5, color= "gray50", size=0.5) +  # add a dotted line at x=1 after flip
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1), 
                     expand=expansion(mult=c(0, .05), add=c(0, 0)),
                     labels=scales::number_format(accuracy = 0.01, decimal.mark='.'))


### BARPLOT 3 - Frequency of autotomy across levels of activity pattern:
Data$Activity<-factor(Data$Diurnality, levels=c(1, 0.5, 0), labels=c("Diurnal", "Cathemeral", "Nocturnal"))
summary_table<-as.data.frame(Data %>% dplyr::group_by(Suborder, Activity) %>% 
                               dplyr::summarise(N_Records=count(AutotomyPresence), N_Total=length(SpeciesName_ReptileDatabase))) %>% ungroup()
summary_table<-data.frame("Suborder"=summary_table$Suborder,
                          "Category"=summary_table$Activity,
                          "Autotomy"=summary_table$N_Records$x,
                          "N_Records"=summary_table$N_Records$freq,
                          "N_Total"=summary_table$N_Total)
summary_table$PropRecords<-summary_table$N_Records/summary_table$N_Total

# Relabel tail condition:
summary_table$TailCondition<-factor(summary_table$Autotomy,
                                    levels=c(0,1),
                                    labels=c("Scar-free tail", "Scarred tail"))

# Specify the order of appearance of each activity pattern category in the barplot:
summary_table$CategoryLabel<-paste0(summary_table$Suborder, "_", summary_table$Category)
summary_table$CategoryLabel<-factor(summary_table$CategoryLabel,
                                    levels=rev(c("Amphisbaenia_Diurnal", "Amphisbaenia_Cathemeral", "Amphisbaenia_Nocturnal",
                                                 "Serpentes_Diurnal", "Serpentes_Cathemeral", "Serpentes_Nocturnal")),
                                    labels=rev(c("Diurnal amphisb.", "  Cathem. amphisb.", "Nocturnal amphisb.",
                                                 "Diurnal snakes", "Cathem. snakes", "  Nocturnal snakes")))

# Create a column to inform the number of examined species per category of activity pattern:
summary_table$N_Label<-NA
for(i in 1:nlevels(summary_table$CategoryLabel)){
  
  selected_row<-which(summary_table$CategoryLabel==levels(summary_table$CategoryLabel)[i])
  
  if(length(selected_row)<=2){
    
    summary_table$N_Label[selected_row[1]]<-paste(summary_table$N_Total[selected_row[1]])
    
  } else {
    
    midrow<-ceiling(mean(selected_row))
    summary_table$N_Label[midrow]<-paste(summary_table$N_Total[midrow])
    
  }
}

# Get summary values irrespective of taxonomic suborder:
summary_table %>% dplyr::group_by(Category,Autotomy) %>% dplyr::summarise(mean(PropRecords))

# Specify the colours for the barplot:
summary_table$fill_var<-paste0(summary_table$TailCondition, "_", summary_table$Suborder)
myColours <- c("#fddbc7", "#d1e5f0", "#f4a582", "#92c5de") # amph_not_uro #snake_not_uro #amph_uro # snake_uro
names(myColours)<-levels(summary_table$fill_var)
fillScale<-scale_fill_manual(values = myColours)

# Build the plot:
MyBarplot[[3]]<-ggplot(data=summary_table, aes(x=CategoryLabel, y=PropRecords, fill=fill_var)) +
  geom_bar(position="stack", stat="identity") +
  geom_text(y=1.025, aes(label=N_Label), size=2, colour="gray25", position=position_stack(vjust=0.5)) +
  coord_flip() +  
  xlab("") + ylab("") +
  fillScale +
  theme(panel.grid.minor = element_blank(), # remove minor grid lines
        panel.grid.major = element_blank(), # remove major grid lines
        panel.background = element_blank(), # white background
        axis.line = element_line(colour="black"), # axis lines aesthetics
        axis.text.y = element_text(hjust=1, vjust=0.5, angle=0, size=8),
        axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=8),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        #axis.title.y=element_text(size=10, colour="black", face="bold", margin=margin(t=0, r=2, b=0, l=0)), # margin between axis.title and axis.values
        #axis.title.x=element_text(size=10, colour="black", face="bold", margin=margin(t=5, r=0, b=0, l=0)), # margin between axis.title and axis.values
        legend.position="none",
        plot.background=element_rect(fill="transparent", colour=NA),
        plot.margin=unit(c(0, 0.65, 0.05, 0), "cm"), # top right bottom left
        panel.spacing = unit(0.1, "lines")) +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1), 
                     expand=expansion(mult=c(0, .05), add=c(0, 0)),
                     labels=scales::number_format(accuracy = 0.01, decimal.mark='.')) +
  scale_x_discrete(position="bottom")


### BARPLOT 4 - Frequency of autotomy across levels of habitat use:
# Get the autotomy frequency among fossorial species:
summary_table<-as.data.frame(Data %>% dplyr::group_by(Suborder, Fossorial) %>% 
                               dplyr::summarise(N_Records=count(AutotomyPresence),
                                                N_Total=length(SpeciesName_ReptileDatabase))) %>% ungroup()
summary_table1<-data.frame("Suborder"=summary_table$Suborder,
                           "Category"=summary_table$Fossorial,
                           "Autotomy"=summary_table$N_Records$x,
                           "N_Records"=summary_table$N_Records$freq,
                           "N_Total"=summary_table$N_Total)
summary_table1$Habitat<-"Fossorial"
summary_table1$PropRecords<-summary_table1$N_Records/summary_table1$N_Total
summary_table1<-summary_table1[summary_table1$Category==1,]

# Get the autotomy frequency among terrestrial species:
summary_table<-as.data.frame(Data %>% dplyr::group_by(Suborder, Terrestrial) %>% 
                               dplyr::summarise(N_Records=count(AutotomyPresence),
                                                N_Total=length(SpeciesName_ReptileDatabase))) %>% ungroup()
summary_table2<-data.frame("Suborder"=summary_table$Suborder,
                           "Category"=summary_table$Terrestrial,
                           "Autotomy"=summary_table$N_Records$x,
                           "N_Records"=summary_table$N_Records$freq,
                           "N_Total"=summary_table$N_Total)
summary_table2$Habitat<-"Terrestrial"
summary_table2$PropRecords<-summary_table2$N_Records/summary_table2$N_Total
summary_table2<-summary_table2[summary_table2$Category==1,]

# Get the autotomy frequency among aquatic species:
summary_table<-as.data.frame(Data %>% dplyr::group_by(Suborder, Aquatic) %>% 
                               dplyr::summarise(N_Records=count(AutotomyPresence),
                                                N_Total=length(SpeciesName_ReptileDatabase))) %>% ungroup()
summary_table3<-data.frame("Suborder"=summary_table$Suborder,
                           "Category"=summary_table$Aquatic,
                           "Autotomy"=summary_table$N_Records$x,
                           "N_Records"=summary_table$N_Records$freq,
                           "N_Total"=summary_table$N_Total)
summary_table3$Habitat<-"Aquatic"
summary_table3$PropRecords<-summary_table3$N_Records/summary_table3$N_Total
summary_table3<-summary_table3[summary_table3$Category==1,]

# Get the autotomy frequency among arboreal species:
summary_table<-as.data.frame(Data %>% dplyr::group_by(Suborder, Arboreal) %>% 
                               dplyr::summarise(N_Records=count(AutotomyPresence),
                                                N_Total=length(SpeciesName_ReptileDatabase))) %>% ungroup()
summary_table4<-data.frame("Suborder"=summary_table$Suborder,
                           "Category"=summary_table$Arboreal,
                           "Autotomy"=summary_table$N_Records$x,
                           "N_Records"=summary_table$N_Records$freq,
                           "N_Total"=summary_table$N_Total)
summary_table4$Habitat<-"Arboreal"
summary_table4$PropRecords<-summary_table4$N_Records/summary_table4$N_Total
summary_table4<-summary_table4[summary_table4$Category==1,]

# Bind all summary tables:
summary_table<-rbind(summary_table1, summary_table3, summary_table2, summary_table4)
rm(summary_table1, summary_table3, summary_table2, summary_table4)

# Specify the order of appearance of each habitat-use category in the barplot:
summary_table$CategoryLabel<-paste0(summary_table$Suborder, "_", summary_table$Habitat)
summary_table$CategoryLabel<-factor(summary_table$CategoryLabel,
                                    levels=rev(c("Amphisbaenia_Fossorial", "Serpentes_Fossorial",
                                                 "Serpentes_Aquatic", "Serpentes_Terrestrial", "Serpentes_Arboreal")),
                                    
                                    labels=rev(c("Fossorial amphisb.", "Fossorial snakes",
                                                 "Aquatic snakes", "Terrestrial snakes", "Arboreal snakes")))

# Create a column to inform the number of examined species per category of activity pattern:
summary_table$N_Label<-NA
for(i in 1:nlevels(summary_table$CategoryLabel)){
  
  selected_row<-which(summary_table$CategoryLabel==levels(summary_table$CategoryLabel)[i])
  
  if(length(selected_row)<=2){
    
    summary_table$N_Label[selected_row[1]]<-paste(summary_table$N_Total[selected_row[1]])
    
  } else {
    
    midrow<-ceiling(mean(selected_row))
    summary_table$N_Label[midrow]<-paste(summary_table$N_Total[midrow])
    
  }
}

# Relabel tail condition:
summary_table$TailCondition<-factor(summary_table$Autotomy,
                                    levels=c(0,1),
                                    labels=c("Scar-free tail", "Scarred tail"))

# Get summary values irrespective of taxonomic suborder:
summary_table %>% dplyr::group_by(Habitat,Autotomy) %>% dplyr::summarise(mean(PropRecords))

# Specify the colours for the barplot:
summary_table$fill_var<-paste0(summary_table$TailCondition, "_", summary_table$Suborder)
myColours <- c("#fddbc7", "#d1e5f0", "#f4a582", "#92c5de") # amph_not_uro #snake_not_uro #amph_uro # snake_uro
names(myColours)<-levels(summary_table$fill_var)
fillScale<-scale_fill_manual(values = myColours)

# Build the plot:
MyBarplot[[4]]<-ggplot(data=summary_table, aes(x=CategoryLabel, y=PropRecords, fill=fill_var)) +
  geom_bar(position="stack", stat="identity") +
  geom_text(y=1.025, aes(label=N_Label), size=2, colour="gray25", position=position_stack(vjust=0.5)) +
  coord_flip() +  
  xlab("") + ylab("") +
  fillScale +
  theme(panel.grid.minor = element_blank(), # remove minor grid lines
        panel.grid.major = element_blank(), # remove major grid lines
        panel.background = element_blank(), # white background
        axis.line = element_line(colour="black"), # axis lines aesthetics
        axis.text.y = element_text(hjust=1, vjust=0.5, angle=0, size=8),
        axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=8),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position="none",
        plot.background=element_rect(fill="transparent", colour=NA),
        plot.margin=unit(c(0, 0.65, 0.05, 0), "cm"), # top right bottom left
        panel.spacing = unit(0.1, "lines")) +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1), 
                     expand=expansion(mult=c(0, .05), add=c(0, 0)),
                     labels=scales::number_format(accuracy = 0.01, decimal.mark='.')) +
  scale_x_discrete(position="bottom")


### BARPLOT 5 - Frequency of autotomy across levels of biome:
summary_table<-as.data.frame(Data %>% dplyr::group_by(Suborder, Biome) %>% 
                               dplyr::summarise(N_Records=count(AutotomyPresence),
                                                N_Total=length(SpeciesName_ReptileDatabase))) %>% ungroup()
summary_table<-data.frame("Suborder"=summary_table$Suborder,
                          "Category"=summary_table$Biome,
                          "Autotomy"=summary_table$N_Records$x,
                          "N_Records"=summary_table$N_Records$freq,
                          "N_Total"=summary_table$N_Total)
summary_table$PropRecords<-summary_table$N_Records/summary_table$N_Total

# Split the autotomy frequency per biome between snakes and amphisbaenians:
summary_table_snakes<-summary_table[summary_table$Suborder=="Serpentes",]
summary_table_amphisb<-summary_table[summary_table$Suborder=="Amphisbaenia",]

# Specify the order of appearance of each biome in the barplot:
summary_table_amphisb[summary_table_amphisb$Autotomy==0,]
summary_table_amphisb$CategoryLabel<-factor(summary_table_amphisb$Category,
                                            levels=c("Flooded Grasslands & Savannas",
                                                     "Tropical & Subtropical Dry Broadleaf Forests",
                                                     "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                                     "Mangroves",
                                                     "Tropical & Subtropical Moist Broadleaf Forests",
                                                     "Temperate Grasslands, Savannas & Shrublands",
                                                     "Montane Grasslands & Shrublands",
                                                     "Deserts & Xeric Shrublands"),
                                            
                                            labels=c("Flooded grass/savan",
                                                     "Trop/subtrop dry broadleaf",
                                                     "Trop/subtrop grass/savan/shrub",
                                                     "Mangroves",
                                                     "Trop/subtrop moist broadleaf",
                                                     "Temperate grass/savan/shrub",
                                                     "Montane grass/shrub",
                                                     "Deserts/xeric shrub"))

summary_table_snakes[summary_table_snakes$Autotomy==0,]
summary_table_snakes$CategoryLabel<-factor(summary_table_snakes$Category,
                                           levels=rev(c("Mediterranean Forests, Woodlands & Scrub",
                                                        "Tropical & Subtropical Dry Broadleaf Forests",
                                                        "Tropical & Subtropical Moist Broadleaf Forests",
                                                        "Mangroves",
                                                        "Flooded Grasslands & Savannas",
                                                        "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                                        "Deserts & Xeric Shrublands",
                                                        "Temperate Grasslands, Savannas & Shrublands")),
                                           
                                           labels=rev(c("Mediterranean",
                                                        "Trop/subtrop dry broadleaf",
                                                        "Trop/subtrop moist broadleaf",
                                                        "Mangroves",
                                                        "Flooded grass/savan",
                                                        "Temperate grass/savan/shrub",
                                                        "Deserts/xeric shrub",
                                                        "Trop/subtrop grass/savan/shrub")))

# Create a column to inform the number of examined species per biome:
summary_table_amphisb$N_Label<-NA
for(i in 1:nlevels(summary_table_amphisb$Category)){
  
  selected_row<-which(summary_table_amphisb$Category==levels(summary_table_amphisb$Category)[i])
  
  if(length(selected_row)<=2){
    
    summary_table_amphisb$N_Label[selected_row[1]]<-paste(summary_table_amphisb$N_Total[selected_row[1]])
    
  } else {
    
    midrow<-ceiling(mean(selected_row))
    summary_table_amphisb$N_Label[midrow]<-paste(summary_table_amphisb$N_Total[midrow])
    
  }
}
summary_table_snakes$N_Label<-NA
for(i in 1:nlevels(summary_table_snakes$Category)){
  
  selected_row<-which(summary_table_snakes$Category==levels(summary_table_snakes$Category)[i])
  
  if(length(selected_row)<=2){
    
    summary_table_snakes$N_Label[selected_row[1]]<-paste(summary_table_snakes$N_Total[selected_row[1]])
    
  } else {
    
    midrow<-ceiling(mean(selected_row))
    summary_table_snakes$N_Label[midrow]<-paste(summary_table_snakes$N_Total[midrow])
    
  }
}

# Build the plots:
MyBarplot[[5]]<-ggplot(data=summary_table_amphisb, aes(x=CategoryLabel, y=PropRecords, fill=as.factor(Autotomy))) +
  geom_bar(position="stack", stat="identity") +
  geom_text(y=1.025, aes(label=N_Label), size=2, colour="gray25", position=position_stack(vjust=0.5)) +
  coord_flip() +  
  xlab("") + ylab("Proportion of specimens with scarred tail") +
  scale_fill_manual(values = c("#fddbc7","#f4a582")) +
  theme(panel.grid.minor = element_blank(), # remove minor grid lines
        panel.grid.major = element_blank(), # remove major grid lines
        panel.background = element_blank(), # white background
        axis.line = element_line(colour="black"), # axis lines aesthetics
        axis.text.y = element_text(hjust=1, vjust=0.5, angle=0, size=7.5),
        axis.text.x = element_blank(), #element_text(hjust=0.5, vjust=0.5, angle=0, size=8),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        #axis.title.y=element_text(size=10, colour="black", face="bold", margin=margin(t=0, r=2, b=0, l=0)), # margin between axis.title and axis.values
        #axis.title.x=element_text(size=10, colour="black", face="bold", margin=margin(t=5, r=0, b=0, l=0)), # margin between axis.title and axis.values
        legend.position="none",
        plot.background=element_rect(fill="transparent", colour=NA),
        plot.margin=unit(c(0, 0.65, 0.05, 0), "cm"), # top right bottom left
        panel.spacing = unit(0.1, "lines")) +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1), 
                     expand=expansion(mult=c(0, .05), add=c(0, 0)),
                     labels=scales::number_format(accuracy = 0.01, decimal.mark='.')) +
  scale_x_discrete(position="bottom") +
  annotate("text", label=bquote('Amphisbaenians'), x=1, y=0.99, color="black", size=2.5, hjust=1)

MyBarplot[[6]]<-ggplot(data=summary_table_snakes, aes(x=CategoryLabel, y=PropRecords, fill=as.factor(Autotomy))) +
  geom_bar(position="stack", stat="identity") +
  geom_text(y=1.025, aes(label=N_Label), size=2, colour="gray25", position=position_stack(vjust=0.5)) +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Proportion of specimens with tail autotomy") +
  scale_fill_manual(values = c("#d1e5f0","#92c5de")) +
  theme(panel.grid.minor = element_blank(), # remove minor grid lines
        panel.grid.major = element_blank(), # remove major grid lines
        panel.background = element_blank(), # white background
        axis.line = element_line(colour="black"), # axis lines aesthetics
        axis.text.y = element_text(hjust=1, vjust=0.5, angle=0, size=7.5),
        axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=8, margin=margin(t=0, r=0, b=5, l=0)),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=10, colour="black", face="bold", margin=margin(t=0, r=5, b=0, l=0)), # margin between axis.title and axis.values
        legend.position="none",
        plot.background=element_rect(fill="transparent", colour=NA),
        plot.margin=unit(c(0, 0.65, 0.05, 0), "cm"), # top right bottom left
        panel.spacing = unit(0.1, "lines")) +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1), 
                     expand=expansion(mult=c(0, .05), add=c(0, 0)),
                     labels=scales::number_format(accuracy = 0.01, decimal.mark='.')) +
  scale_x_discrete(position="bottom") +
  annotate("text", label=bquote('Snakes'), x=1, y=0.99, color="black", size=2.5, hjust=1)


# Build the multipanel plot:
plot_layout <- rbind(c(1,1,2,2), #1
                     c(1,1,2,2), #2
                     c(1,1,2,2), #3
                     c(1,1,2,2), #4
                     c(1,1,2,2), #5
                     c(1,1,2,2), #6
                     c(1,1,2,2), #7
                     c(1,1,2,2), #8
                     c(1,1,3,3), #9
                     c(1,1,3,3), #10
                     c(1,1,3,3), #11
                     c(1,1,3,3), #12
                     c(1,1,4,4), #13
                     c(1,1,4,4), #14
                     c(1,1,4,4), #15
                     c(1,1,5,5), #16
                     c(1,1,5,5), #17
                     c(1,1,5,5), #18
                     c(1,1,5,5), #18
                     c(1,1,6,6), #19
                     c(1,1,6,6), #20
                     c(1,1,6,6), #21
                     c(1,1,6,6), #22
                     c(1,1,6,6), #23
                     c(1,1,6,6)) #24
Multipanel_plot<-gridExtra::grid.arrange(grobs = MyBarplot, layout_matrix = plot_layout,
                                         padding = unit(-0.05, "line"), align="hv")

# Add small caption letters for each panel:
Multipanel_plot2 <- as_ggplot(Multipanel_plot) + 
  cowplot::draw_plot_label(label = c("a)", "b)", "c)", "d)", "e)"), size = 10,
                           x = c(-0.005, 0.48, 0.48, 0.48, 0.48),
                           y = c(1.01, 1.01, 0.685, 0.525, 0.39)); Multipanel_plot2 # Add labels

# Export to disk:
ggsave("Figure2.png", plot=Multipanel_plot2, width=11, height=6, units="in", bg = "transparent")
ggsave("Figure2.pdf", plot=Multipanel_plot2, width=11, height=6, units="in", bg = "transparent")

#####

# STEP 5 - Plot autotomy frequency per species along the phylogeny
##############################################################################################################
# STEP 5 - Plot autotomy frequency per species along the phylogeny
rm(list=setdiff(ls(), c("Data")))

# The function will remove special characters, lowercase the names, and replace white space by underline symbol
get_std_sciname<-function(x){
  edited_name<-stringr::str_to_lower(x) # lowercase every name
  edited_name<-gsub('\\\\(', '', edited_name) # remove parentheses
  edited_name<-gsub('\\\\)', '', edited_name) # remove parentheses
  edited_name<-gsub("'", '"', edited_name) # replace single straight quotes by double straight quotes
  edited_name<-gsub('"', '', edited_name) # remove double straight quotes
  edited_name<-as.factor(stringr::str_squish(edited_name)) # remove whitespace from start and end of string
  edited_name<-as.factor(gsub('([[:punct:]])|\\\\s+', '_', edited_name)) # replace white space by underline symbol
}

# Load a set of 100 fully-sampled phylogenies for reptiles, available at: https://vertlife.org/phylosubsets
phylo_tree_set<-ape::read.tree('Tonini_etal_100trees.trees')

# First, check for consistency of binomials between phylogenetic tree and the autotomy database:
phylo_names<-as.data.frame(phylo_tree_set[[1]]$tip.label); names(phylo_names)<-"tip.label"
phylo_names$edited_name<-get_std_sciname(phylo_names$tip.label)
Data$edited_name<-get_std_sciname(Data$SpeciesName_ReptileDatabase)

# Identify species included in the autotomy dataset but absent in the phylogeny:
NotInPhylo<-merge(Data, phylo_names, by="edited_name", all.x=T)
NotInPhylo<-unique(NotInPhylo[is.na(NotInPhylo$tip.label),1:2])
NotInPhylo$edited_name

# Update some taxonomic names to match those used in the phylogeny:
Data$PhyloName<-Data$edited_name
Data$PhyloName<-gsub('chironius_gouveai', 'chironius_bicarinatus', Data$PhyloName)  # Entiauspe-Neto, Omar Machado; Mariana Lyra, Claudia Koch, Fernando  Quintela, Arthur Diesel Abegg, Daniel Loebmann 2020. Taxonomic Revision of Chironius bicarinatus (Wied 1820) (Serpentes: Colubridae), with Description of a New Species.\\ Herpetological Monographs 34 (1): 98-115
Data$PhyloName<-gsub('leptodeira_ornata', 'leptodeira_septentrionalis', Data$PhyloName) # Reptile Database
Data$PhyloName<-gsub('palusophis_bifossatus', 'mastigodryas_bifossatus', Data$PhyloName) # Reptile Database
Data$PhyloName<-gsub('zamenis_scalaris', 'rhinechis_scalaris', Data$PhyloName) # Reptile Database

# Get the standardised species name for matching operations:
SpeciesSubset<-unique(Data[,.(PhyloName)])
row.names(SpeciesSubset)<-SpeciesSubset$PhyloName

# For each one of the 100 fully-sampled trees, extract the 43 species included in the autotomy database:
tree_set<-list()
for(i in 1:100){
  
  # Get one tree and standardised tip labels for matching operations:
  phylo_tree<-phylo_tree_set[[i]]
  phylo_tree$tip.label<-get_std_sciname(phylo_tree$tip.label)
  
  # Get the phylogeny for snakes only:
  my_tree<-(geiger::treedata(phy=phylo_tree, data=SpeciesSubset, sort=T))$phy
  my_tree$tip.label<-as.character(my_tree$tip.label)
  
  # Store the subsampled tree i in a list:
  tree_set[[i]]<-my_tree
}

# Trim the consensus tree, available at: https://datadryad.org/stash/dataset/doi:10.5061/dryad.db005
consensus_tree<-ape::read.tree('squam_shl_new_Consensus_9755.tre')
consensus_tree$tip.label<-get_std_sciname(consensus_tree$tip.label)
my_tree<-(geiger::treedata(phy=consensus_tree, data=SpeciesSubset, sort=T))$phy
my_tree$tip.label<-as.character(my_tree$tip.label)

# Export the trimmed trees as an RData object (for latter use):
save(tree_set, file = "Subset100Trees_43taxa.Rdata")
save(my_tree, file = "ConsensusTree_43taxa.Rdata")
Data<-Data[,-c("edited_name")]
rm(list=setdiff(ls(), c("Data")))

# Load the trimmed consensus tree created in the previous step:
load("ConsensusTree_43taxa.Rdata")

# Get new tip labels for better aesthetics:
renamed_labels<-as.data.frame(my_tree$tip.label); names(renamed_labels)<-"tip.label"
renamed_labels$Genus<-gsub("([A-Za-z]+).*", "\\\\1", renamed_labels$tip.label) # extract the genus
renamed_labels$Genus<-str_to_title(renamed_labels$Genus) # upper case the genus initial
Epithet<-(str_split(renamed_labels$tip.label, "_")) # split species names in multiple columns
renamed_labels$Epithet<-NA # create a empty column to hold the epithet value
for(i in 1:nrow(renamed_labels)){
  
  if(length(Epithet[[i]])>=3){renamed_labels$Epithet[i]<-Epithet[[i]][[3]]
  
  } else {
    
    renamed_labels$Epithet[i]<-Epithet[[i]][[2]]
  }
} # fill the epithet value
renamed_labels$Binomial<-paste0(renamed_labels$Genus, " " ,renamed_labels$Epithet)
my_tree$tip.label<-renamed_labels$Binomial

# Get the autotomy frequency per species:
summary_table<-as.data.frame(Data %>% dplyr::group_by(SpeciesName_ReptileDatabase, Suborder) %>% 
                               dplyr::summarise(N_Records=count(AutotomyPresence),
                                                N_Total=length(SpeciesName_ReptileDatabase))) %>% ungroup()
summary_table<-data.frame("Species"=summary_table$SpeciesName_ReptileDatabase,
                          "Group"=summary_table$Suborder,
                          "Autotomy"=summary_table$N_Records$x,
                          "N_Records"=summary_table$N_Records$freq,
                          "N_Total"=summary_table$N_Total)
summary_table$PropRecords<-summary_table$N_Records/summary_table$N_Total
summary_table<-summary_table[summary_table$Autotomy==1,]

# Adjust binomials that differ between the phylogenetic tree and the dataset:
my_tree$tip.label<-gsub('Leptodeira septentrionalis', 'Leptodeira ornata', my_tree$tip.label) # Reptile Database
my_tree$tip.label<-gsub('Mastigodryas bifossatus', 'Palusophis bifossatus', my_tree$tip.label) # Reptile Database
my_tree$tip.label<-gsub('Rhinechis scalaris', 'Zamenis scalaris', my_tree$tip.label) # Reptile Database
row.names(summary_table)<-summary_table$Species

# Get species present in both phylogenetic tree and autotomy database:
consensus_tree<-(geiger::treedata(phy=my_tree, data=summary_table, sort=T))$phy
summary_table<-as.data.frame((geiger::treedata(phy=my_tree, data=summary_table, sort=T))$data)
consensus_tree$tip.label<-as.character(consensus_tree$tip.label)

# Compute the autotomy frequency per species:
summary_table<-summary_table[,c("Species", "Group", "N_Records", "N_Total", "PropRecords")]
summary_table$N_Records<-as.numeric(summary_table$N_Records)
summary_table$N_Total<-as.numeric(summary_table$N_Total)
summary_table$PropRecords<-as.numeric(summary_table$PropRecords)

# Specify colours for tree branches representing amphisbaenians or snakes:
groupInfo<-list(Amphisbaenia=as.character(summary_table[(summary_table$Group=="Amphisbaenia"),1]),
                Serpentes=as.character(summary_table[(summary_table$Group=="Serpentes"),1]))
consensus_tree <- groupOTU(consensus_tree, groupInfo) # group tip labels

# Plot the tree with branches coloured according to each group:
b<-ggtree(consensus_tree, layout="circular", branch.length=T, aes(color=group), size=.75) +
  scale_color_manual(values=c("gray", "#b2182b", "#2166ac")); b

# Create a new variable containing values of autotomy frequency for amphisbaenians only:
summary_table$NewVar<-NA
summary_table[summary_table$Group=="Amphisbaenia",]$NewVar<-summary_table[summary_table$Group=="Amphisbaenia",]$PropRecords
mybreaks<-c(# Min
  round(ceiling(min(summary_table$NewVar*100, na.rm=T)),3),
  
  # Midrange
  (((round(floor(max(summary_table$NewVar*100, na.rm=T)),3) - round(floor(min(summary_table$NewVar*100, na.rm=T)),3))/2)+
     (round(floor(min(summary_table$NewVar*100, na.rm=T)),3))),
  
  # Max
  round(floor(max(summary_table$NewVar*100, na.rm=T)),3))
mybreaks<-mybreaks/100 # (min, midrange, max)

# Plot the external heatmap showing autotomy frequency for amphisbaenian species:
b1 <- gheatmap(b, summary_table[, "NewVar", drop=F], offset=0, width=.1, colnames_angle=90, colnames=F, colnames_offset_y=.25) +
  scale_fill_distiller(type = "seq", 
                       palette = "Reds",
                       direction = 1,
                       aesthetics = "fill",
                       breaks=mybreaks,
                       na.value=NA) +
  theme(legend.direction = "horizontal") +
  guides(group="none", color="none",
         fill=guide_colorbar(title="Autotomy Freq.", title.position="top", title.hjust=.5,
                             label=T, nbin=5, direction="horizontal", label.position = "bottom",
                             draw.ulim=T, draw.llim=T, frame.colour="black", ticks=F)); b1
legend_amphisb<-cowplot::get_legend(b1) 

# Same as above, but for snakes:
summary_table$NewVar<-NA
summary_table[summary_table$Group=="Serpentes",]$NewVar<-summary_table[summary_table$Group=="Serpentes",]$PropRecords
mybreaks<-c(# Min
  round(ceiling(min(summary_table$NewVar*100, na.rm=T)),3),
  (((round(floor(max(summary_table$NewVar*100, na.rm=T)),3) - round(floor(min(summary_table$NewVar*100, na.rm=T)),3))/2)+
     (round(floor(min(summary_table$NewVar*100, na.rm=T)),3))),
  round(floor(max(summary_table$NewVar*100, na.rm=T)),3))
mybreaks<-mybreaks/100 # (min, midrange, max)
b1 <- b1 + ggnewscale::new_scale_fill()
b2 <- gheatmap(b1, summary_table[, "NewVar", drop=F], offset=0, width=.1, colnames_angle=90, colnames=F, colnames_offset_y=.25) +
  scale_fill_distiller(type = "seq", 
                       palette = "Blues",
                       direction = 1,
                       aesthetics = "fill",
                       breaks=mybreaks,
                       #limits=c(0.05,.75),
                       na.value=NA) +
  theme(legend.direction = "horizontal") +
  guides(group="none", color="none",
         fill=guide_colorbar(title="", title.position="top", title.hjust=.5,
                             label=T, nbin=5, direction="horizontal", label.position = "bottom",
                             draw.ulim=T, draw.llim=T, frame.colour="black", ticks=F)); b2
legend_snakes<-cowplot::get_legend(b2) 

# Rename grobs of the second legend, and extract only the second grob (red colorbar):
legend_snakes$layout[,7]<-c("A", "B", "legend.box.background")
legend_snakes<-gtable::gtable_filter(legend_snakes, "B")

# Add tip labels:
b3 <- b2 + 
  geom_tiplab2(size=4, align=T, offset=25,  fontface="italic") +
  theme(plot.margin=unit(c(3, 0, 3, 0), "cm"), # top, right, bottom, left
        legend.position="none"); b3

# Build the phylogeny plot:
MyPhylogeny <- cowplot::ggdraw() +
  cowplot::draw_plot(b3) +
  cowplot::draw_plot(legend_snakes, x=0.35, y=0.2, width=0.35, height=0.35, halign=0.5) +
  cowplot::draw_plot(legend_amphisb, x=0.3283, y=0.255, width=0.35, height=0.35, halign=0.5)

# Export to disk:
ggsave("Figures/Figure3.png", plot=MyPhylogeny, width=10, height=10, units="in", bg = "transparent", limitsize=F)
ggsave("Figures/Figure3.pdf", plot=MyPhylogeny, width=10, height=10, units="in", bg = "transparent", limitsize=F)
# Other aesthetic components were included outside the R environment, including photos of species.

#####

# STEP 6 - Analyse determinants of autotomy probability using Generalised Linear Mixed Effect Model
##############################################################################################################
# STEP 6 - Analyse determinants of autotomy probability using Generalised Linear Mixed Effect Model
rm(list=setdiff(ls(), c("Data")))
set.seed(123)

# Separate data only on sexed specimens to run analysis without unsexed specimens:
# Data<-droplevels(Data[Data$SexBinary!=0.5,])

# Check how each covariate is being read by R:
summary(Data)
Data$Verticality<-as.numeric(as.character(Data$Verticality))

# Run the GLMM for determinants of autotomy probability:
m1full_group_StdSVL <- glmer(AutotomyPresence ~ 
                               
                               # Fixed effects
                               Diurnality + Verticality + SexBinary + LifeStageBinary + StdTemp + StdPrec + Biome + Tropicality + group_StdSVL + Suborder + 
                               
                               # Using only species as Random effect
                               (1|SpeciesName_ReptileDatabase),
                             
                             data=Data, family=binomial, control=glmerControl(optimizer = "Nelder_Mead"), nAGQ=100)

# Perform variable selection using likelihood-ratio tests to find the most parsimonious model:
var_selection_output<-list()
updated_model_group_StdSVL<-m1full_group_StdSVL
for(i in 1:10){ # 10 predictors
  
  # Backward selection based on Likelihood Ratio Test:
  variable_importance<-drop1(updated_model_group_StdSVL, test="Chisq")
  variable_importance<-variable_importance[order(variable_importance$`Pr(Chi)`, decreasing=F),]
  
  # Store the results of each iteration:
  var_selection_output[[i]]<-variable_importance
  
  # Proceed with the backward selection procedure until all remaining predictors are significant:
  if(variable_importance$`Pr(Chi)`[(nrow(variable_importance)-1)]>=0.05){ # detect non-significant variables
    
    # Specify the model formula:
    formula_drop1<-as.formula(paste("AutotomyPresence ~", paste(row.names(variable_importance)[1:(nrow(variable_importance)-2)], collapse= "+"), " + (1 | SpeciesName_ReptileDatabase)"))
    
    # Update the (GLMM) model in test:
    updated_model_group_StdSVL<-glmer(formula_drop1,  data=Data, family=binomial, control=glmerControl(optimizer = "Nelder_Mead"), nAGQ=100)
    
    # Remove unused terms:
    rm(variable_importance, formula_drop1)
    
  } # end of if condition
  
} # end of for loop

# Check the outputs for the full and the most parsimonious models: 
var_selection_output[[1]] # full model
var_selection_output[[length(var_selection_output)]] # most parsimonious model

# Define the final model selected and remove unnecessary objects:
final_model<-updated_model_group_StdSVL
summary(final_model)
rsq::rsq.glmm(final_model)

# Get the confidence intervals (CI) around each predictor:
se_1 <- sqrt(diag(vcov(final_model)))
tab_1 <- cbind(Est = fixef(final_model), LL = fixef(final_model) - 1.96 * se_1, UL = fixef(final_model) + 1.96 *se_1); tab_1
round(exp(tab_1), 3) # exponentiate the estimates and CIs for odds ratios instead of coefficients on the logit scale

# Export the results:
save(final_model, Data, variable_importance, var_selection_output, tab_1, file="GLMM_output_group_StdSVL.Rdata")

# Same as above, but for body size rescaled across all species:
# In using body size rescaled across all species, it is necessary to add nested random effect (Species + Species:Sex + Species:LifeStage):
m1full_StdSVL <- glmer(AutotomyPresence ~ 
                         
                         # Fixed effects
                         Diurnality + Verticality + SexBinary + LifeStageBinary + StdTemp + StdPrec + Biome + Tropicality + StdSVL +  Suborder +  
                         
                         # Random effects
                         ((1|SpeciesName_ReptileDatabase) + (1|SpeciesName_ReptileDatabase:Sex) + (1|SpeciesName_ReptileDatabase:LifeStage)),
                       
                       data=Data, family=binomial, control=glmerControl(optimizer = "Nelder_Mead"), nAGQ=1)

var_selection_output_StdSVL<-list()
updated_model_StdSVL<-m1full_StdSVL
for(i in 1:10){ # 10 predictors
  
  variable_importance_StdSVL<-drop1(updated_model_StdSVL, test="Chisq")
  variable_importance_StdSVL<-variable_importance_StdSVL[order(variable_importance_StdSVL$`Pr(Chi)`, decreasing=F),]
  var_selection_output_StdSVL[[i]]<-variable_importance_StdSVL
  
  if(variable_importance_StdSVL$`Pr(Chi)`[(nrow(variable_importance_StdSVL)-1)]>=0.05){
    formula_drop1<-as.formula(paste("AutotomyPresence ~", paste(row.names(variable_importance_StdSVL)[1:(nrow(variable_importance_StdSVL)-2)], collapse= "+"),
                                    " + ((1|SpeciesName_ReptileDatabase) + (1|SpeciesName_ReptileDatabase:Sex) + (1|SpeciesName_ReptileDatabase:LifeStage))"))
    updated_model_StdSVL<-glmer(formula_drop1,  data=Data, family=binomial, control=glmerControl(optimizer = "Nelder_Mead"), nAGQ=1)
    rm(variable_importance_StdSVL, formula_drop1)
  } # end of if condition
  
} # end of for loop
var_selection_output_StdSVL[[length(var_selection_output_StdSVL)]]
final_model<-updated_model_StdSVL
summary(final_model)
rsq::rsq.glmm(final_model)
se_1 <- sqrt(diag(vcov(final_model)))
tab_1 <- cbind(Est = fixef(final_model), LL = fixef(final_model) - 1.96 * se_1, UL = fixef(final_model) + 1.96 *se_1); tab_1
round(exp(tab_1), 3)
save(final_model, Data, variable_importance_StdSVL, var_selection_output_StdSVL, tab_1, file="GLMM_output_StdSVL.Rdata")

#####

# STEP 7 - Get predicted probabilities of tail loss across values of continuous covariates
##############################################################################################################
# STEP 7 - Get predicted probabilities of tail loss across values of continuous covariates
rm(list=ls())

# Load the Data and the output of GLMM:
load("GLMM_output_group_StdSVL.Rdata")

# Get the predicted probabilities from graphing. This is done separately for each predictor:
tmpdat <- Data

# Instead of holding the predictor constant (e.g. using the average), use a set of values along the predictor range.
# Probabilities are then estimated for each value to build a confidence interval of probabilities.
jvalues <- with(Data, seq(from = min(group_StdSVL), to = max(group_StdSVL), length.out = 100))

# Calculate the predicted probabilities and store in a list
pp_group_StdSVL <- lapply(jvalues, function(j) {
  tmpdat$group_StdSVL <- j
  predict(final_model, newdata = tmpdat, type = "response")
})

# Get the average marginal predicted probability across a few different predictor values:
sapply(pp_group_StdSVL[c(1, 20, 40, 60, 80, 100)], mean)

# Get the means with lower and upper quartiles:
plotdat <- t(sapply(pp_group_StdSVL, function(x) {
  c(M = mean(x), quantile(x, c(0.25, 0.75)))
}))

# Add in the predictor values and convert to data frame:
plotdat <- as.data.frame(cbind(plotdat, jvalues))
colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "StdSVL")
head(plotdat)

# Plot the autotomy probability across SVL:
Myplot<-list()
Myplot[[1]]<-ggplot(plotdat, aes(x = StdSVL, y = PredictedProbability)) +
  geom_linerange(aes(ymin = Lower, ymax = Upper), colour="gray50") + 
  geom_line(size = 2, colour="black") + 
  ylim(c(0, 1)) +
  ylab("Probability of tail loss") +
  xlab("Within-species std. body size") +
  theme(panel.grid.minor = element_blank(), # remove minor gridlines
        panel.grid.major = element_blank(), # remove major gridlines
        panel.background = element_blank(), # white background
        axis.line = element_line(colour="black"), # axis lines aesthetics
        axis.text.y = element_text(hjust=0.5, vjust=0.5, angle=0, size=10),
        axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=10),
        axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=12, colour="black", face="bold", margin=margin(t=0, r=5, b=0, l=0)), # margin between axis.title and axis.values
        axis.title.x=element_text(size=12, colour="black", face="bold", margin=margin(t=5, r=0, b=0, l=0)), # margin between axis.title and axis.values
        legend.position="none",
        plot.background=element_blank(),
        panel.spacing=unit(0,"null"))

# Same as above, but for temperature:
jvalues <- with(Data, seq(from = min(StdTemp), to = max(StdTemp), length.out = 100))
pp_StdTemp <- lapply(jvalues, function(j) {
  tmpdat$StdTemp <- j
  predict(final_model, newdata = tmpdat, type = "response")
})
sapply(pp_StdTemp[c(1, 20, 40, 60, 80, 100)], mean)
plotdat2 <- t(sapply(pp_StdTemp, function(x) {c(M = mean(x), quantile(x, c(0.25, 0.75)))}))
plotdat2 <- as.data.frame(cbind(plotdat2, jvalues))
colnames(plotdat2) <- c("PredictedProbability", "Lower", "Upper", "StdTemp")

# Convert temperature back to raw values:
summary(Data$AMT)
plotdat2$AMT<-(plotdat2$StdTemp*(max(Data$AMT, na.rm=T)-min(Data$AMT, na.rm=T)))+min(Data$AMT, na.rm=T)

# Build the plot:
Myplot[[2]]<-ggplot(plotdat2, aes(x = AMT, y = PredictedProbability)) +
  geom_linerange(aes(ymin = Lower, ymax = Upper), color="gray50") + 
  geom_line(size = 2, color="black") + 
  ylim(c(0, 1)) +
  ylab("Probability of tail loss") +
  xlab("Annual average temperature (°C)") +
  theme(panel.grid.minor = element_blank(), # remove minor gridlines
        panel.grid.major = element_blank(), # remove major gridlines
        panel.background = element_blank(), # white background
        axis.line = element_line(colour="black"), # axis lines aesthetics
        axis.text.y = element_text(hjust=0.5, vjust=0.5, angle=0, size=10),
        axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=10),
        axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=12, colour="black", face="bold", margin=margin(t=0, r=5, b=0, l=0)), # margin between axis.title and axis.values
        axis.title.x=element_text(size=12, colour="black", face="bold", margin=margin(t=5, r=0, b=0, l=0)), # margin between axis.title and axis.values
        legend.position="none",
        plot.background=element_blank(),
        panel.spacing=unit(0,"null"))

# Specify a function to add number of samples to a boxplot:
stat_box_data <- function(x, lower_limit=0.025) {
  return(data.frame(y = 0.025, label = paste(format(length(x), big.mark=",", decimal.mark=".", scientific=FALSE))))
}

# Create a boxplot function to get median predictor values between autotomised and non-autotomised specimens:       
Autotomy_Boxplot<-function(x_var, y_var, yAxis_lab){
  
  x_var<-as.factor(x_var)
  
  # Get the x and y axis position of the post hoc small caption letters:
  value_max = as.data.frame(cbind(x_var, y_var)) %>% 
    dplyr::group_by(x_var) %>%
    dplyr::summarize(max_value = max(as.numeric(as.character(y_var))))
  
  # Relabel the levels of x_var:
  value_max$x_var<-factor(value_max$x_var, levels = c(1, 2), labels = c(0,1))
  
  # Get the multiple comparisons among groups:
  posthoctest<-agricolae::kruskal(y=y_var, trt=x_var, alpha = 0.05, p.adj="bonferroni", group=TRUE) 
  sig.letters <- posthoctest$groups[rev(order(row.names(posthoctest$groups))), ]
  
  # Build the plot:
  ggplot(data = data.frame(x_var=x_var, y_var=y_var), aes(x=x_var, y=y_var)) +
    
    # Add geom layers:
    geom_boxplot(aes(fill=x_var), colour="black", outlier.alpha=0.05, notch=F, na.rm=T) + 
    geom_text(data=value_max, aes(x=value_max$x_var, y=(value_max$max_value)*1.05, label=sig.letters$groups), 
              position=position_dodge(width = 2), size=3, vjust=0) +
    stat_summary(fun.data=stat_box_data, size=3, geom="text", hjust=0.5, vjust=0.1) +
    
    # Set the scales:
    scale_fill_manual(values=c("gray90", "gray40")) + 
    scale_colour_manual(values=c("gray90", "gray40")) + 
    scale_x_discrete(breaks=c(0,1), labels=c("Non-autotomised\\ntail", "Autotomised\\ntail")) +
    
    # Set other aesthetics:
    labs(x="", y=yAxis_lab) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          panel.spacing = unit(0, "lines"), 
          axis.title = element_text(size=8, colour="black"),
          axis.line = element_line(colour="black"),
          axis.text = element_text(size=7, colour="black"),
          axis.text.x = element_text(size=8, colour="black"),
          axis.title.x = element_text(margin=margin(t=2, r=0, b=0, l=0)),
          axis.title.y = element_text(margin=margin(t=0, r=2, b=0, l=0)),
          legend.position="none",
          plot.background=element_rect(fill="transparent", colour=NA),
          strip.background = element_blank(),
          strip.placement = "none",
          strip.text.x = element_text(size=8, color="black"))
}

# Build boxplots to use as inset figures:
MyBoxplot<-list()
MyBoxplot[[1]]<-Autotomy_Boxplot(Data$AutotomyPresence, Data$group_StdSVL, "Body size") + 
  scale_fill_manual(values=c("gray80","gray30")) +
  scale_color_manual(values=c("gray80","gray30")); MyBoxplot[[1]]

MyBoxplot[[2]]<-Autotomy_Boxplot(Data$AutotomyPresence, Data$AMT, "Temperature") + 
  scale_fill_manual(values=c("gray80","gray30")) +
  scale_color_manual(values=c("gray80","gray30")); MyBoxplot[[2]]

# Build the multipanel plot with geographical patterns of beta-diversity and biotic homogenization:
PanelA<-ggpubr::ggarrange(Myplot[[1]] +
                            annotation_custom(ggplotGrob(MyBoxplot[[1]]), xmin=-0.05, xmax=0.7, ymin=0.25, ymax=1.05),
                          labels="", align="hv",
                          font.label=list(size=14, colour = "black"), ncol=1, nrow=1)
PanelB<-ggpubr::ggarrange(Myplot[[2]] +
                            annotation_custom(ggplotGrob(MyBoxplot[[2]]), 
                                              xmin=7.50, xmax=22, ymin=0.25, ymax=1.05),
                          labels="", align="hv",
                          font.label=list(size=14, colour = "black"), ncol=1, nrow=1)

Multipanel_plot<-ggpubr::ggarrange(PanelA, PanelB,
                                   labels=c("a)", "b)"), align="hv",
                                   font.label=list(size=12, colour = "black"), ncol=2, nrow=1); Multipanel_plot

# Save to disk:
ggsave("Figure4.png", plot=Multipanel_plot, width=10, height=4, units="in", bg = "transparent")
ggsave("Figure4.pdf", plot=Multipanel_plot, width=10, height=4, units="in", bg = "transparent")

#####

# STEP 8 - Check if model residuals show phylogenetic or spatial autocorrelation
##############################################################################################################
# STEP 8 - Check if model residuals show phylogenetic or spatial autocorrelation
rm(list=ls())

# Load the autotomy dataset, the output of GLMM, and the trimmed phylogenies:
load("Subset100Trees_43taxa.Rdata") # tree_set object built in step 5
load("GLMM_output_group_StdSVL.Rdata") # glmm outputs computed in step 6

# Extract the GLMM residuals for each specimen:
Data$autotomy_resid<-resid(final_model)

# The function will remove special characters, lowercase the names, and replace white space by underline symbol
get_std_sciname<-function(x){
  edited_name<-stringr::str_to_lower(x) # lowercase every name
  edited_name<-gsub('\\\\(', '', edited_name) # remove parentheses
  edited_name<-gsub('\\\\)', '', edited_name) # remove parentheses
  edited_name<-gsub("'", '"', edited_name) # replace single straight quotes by double straight quotes
  edited_name<-gsub('"', '', edited_name) # remove double straight quotes
  edited_name<-as.factor(stringr::str_squish(edited_name)) # remove whitespace from start and end of string
  edited_name<-as.factor(gsub('([[:punct:]])|\\\\s+', '_', edited_name)) # replace white space by underline symbol
}

# Update taxonomic names to match those used in the phylogeny:
Data$PhyloName<-get_std_sciname(Data$SpeciesName_ReptileDatabase)
Data$PhyloName<-gsub('chironius_gouveai', 'chironius_bicarinatus', Data$PhyloName)  # Entiauspe-Neto, Omar Machado; Mariana Lyra, Claudia Koch, Fernando  Quintela, Arthur Diesel Abegg, Daniel Loebmann 2020. Taxonomic Revision of Chironius bicarinatus (Wied 1820) (Serpentes: Colubridae), with Description of a New Species.\\ Herpetological Monographs 34 (1): 98-115
Data$PhyloName<-gsub('leptodeira_ornata', 'leptodeira_septentrionalis', Data$PhyloName) # Reptile Database
Data$PhyloName<-gsub('palusophis_bifossatus', 'mastigodryas_bifossatus', Data$PhyloName) # Reptile Database
Data$PhyloName<-gsub('zamenis_scalaris', 'rhinechis_scalaris', Data$PhyloName) # Reptile Database

# For all specimens, create a unique tip.label nested within the respective species:
Data$PhyloName2<-NA
Data<-as.data.frame(Data)
for(i in 1:nlevels(as.factor(Data$PhyloName))){
  
  # Identify the number of specimens in the species i:
  SelectedRows<-Data[Data$PhyloName==levels(as.factor(Data$PhyloName))[i],]$PhyloName
  
  # Add the new label to represent phylogeny tips:
  Data[Data$PhyloName==levels(as.factor(Data$PhyloName))[i],]$PhyloName2<-paste0(SelectedRows, "_", c(1:length(SelectedRows)))
  
}

# Prepare workspace for parallel computing:
library(foreach)
library(doParallel)
cl<-makePSOCKcluster(detectCores()-1, type="SOCK")
registerDoParallel(cl)
getDoParWorkers()

# Merge the node number where each species need to added in the phylogeny (do not try to run it locally):
dir.create("Correlograms", showWarnings = FALSE)
foreach(i = 1:100,
        .export = 'rbind',
        .packages = c("data.table", "ggtree", "phangorn", "geiger", "phylobase", "phylosignal"))  %dopar% {
          
          # Select one trimmed fully-sampled tree:
          my_tree<-tree_set[[i]]
          
          # Identify the node holding each species: 
          TreeNodes<-ggtree::fortify(my_tree)
          TreeNodes<-subset(TreeNodes, isTip)
          TreeNodes<-TreeNodes[,c("label", "node")]
          Data1<-merge(Data, TreeNodes, by.x="PhyloName", by.y="label", all.x=T)
          
          # Add specimens to their respective species-level node in the tree:
          my_tree <- phangorn::add.tips(my_tree, Data1$PhyloName2, Data1$node)
          
          # Drop the non-used tips (e.g., the binomials without specimen number):
          Datasubset<-as.data.table(Data1)
          SpeciesSubset<-unique(Datasubset[,.(PhyloName2, autotomy_resid)])
          row.names(SpeciesSubset)<-SpeciesSubset$PhyloName
          my_tree_subset<-(geiger::treedata(phy=my_tree, data=SpeciesSubset, sort=T))
          
          # Create a phylo4 object including GLMM model residuals:
          phylo4d_filter<-phylobase::phylo4d(x=my_tree_subset$phy, 
                                             data.frame(GLMM_resid=as.numeric(as.character(my_tree_subset$data[,2]))))
          
          # Compute the phylogenetic correlogram:
          phy.cor<-phylosignal::phyloCorrelogram(p4d=phylo4d_filter, trait=names(tdata(phylo4d_filter)), dist.phylo="patristic", n.points=14, ci.bs=100)
          correlogram_data<-as.data.frame(phy.cor[[1]])
          names(correlogram_data)<-c("dist.class", "lower_ci", "upper_ci", "coef")
          correlogram_data$Iter<-i
          correlogram_data$N_class<-1:14
          
          # Export each iteration:
          return(
            save(correlogram_data, file=paste0("Correlograms/Tree_", i,".RData"))
          )
        }

# Load RData files with the phylogenetic correlogram of each subsampled tree:
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
correlogram_output<-lapply(as.list(list.files(path="Correlograms/", pattern='Tree_', full.names=T)), function(x) loadRData(x))
correlogram_output<-as.data.table(rbindlist(correlogram_output))

# Compute the average output for the phylogenetic correlogram across the 100 trees:
avg_phylo_correlogram<-correlogram_output[, .(Distance=mean(dist.class, na.rm=T),
                                            Lower_CI=mean(lower_ci, na.rm=T),
                                            Upper_CI=mean(upper_ci, na.rm=T),
                                            MoranI_coef=mean(coef, na.rm=T)),
                                        by = .(N_class)]

# Load the output of the 'avg_phylo_correlogram':
avg_phylo_correlogram<-loadRData("PhyloCorrelogram.Rdata")

# Plot the phylogenetic correlogram:
MyCorrelograms<-list()
MyCorrelograms[[1]]<-ggplot(avg_phylo_correlogram, aes(x = Distance, y = MoranI_coef)) +
  geom_point(size=3)+
  geom_line()+
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  ylim(c(-1, 1)) +
  ylab("Moran's I \\u2013 GLMM residuals") +
  xlab("Phylogenetic distance (mya)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        axis.text.y = element_text(hjust=0.5, vjust=0.5, angle=0, size=10),
        axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=10),
        axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=12, colour="black", face="bold", margin=margin(t=0, r=5, b=0, l=0)),
        axis.title.x=element_text(size=12, colour="black", face="bold", margin=margin(t=5, r=0, b=0, l=0)),
        legend.position="none",
        plot.background=element_blank(),
        panel.spacing=unit(0,"null"))

# Get the geographical coordinates converted for the Euclidean space:
coords<-Data[,c("Latitude", "Longitude")]
coords<-SoDA::geoXY(latitude=coords$Latitude, longitude=coords$Longitude, lat0=0, lon0=0, unit = 1)
coords<-coords/1000 # Convert the distance to km 
spatial_correlogram<-as.data.frame(pgirmess::correlog(coords=coords, z=Data$autotomy_resid, method="Moran", nbclass=14))
names(spatial_correlogram)<-c("Distance", "MoranI_coef", "pvalue", "N_pairs")

# Plot the Moran's I spatial correlogram for the GLMM residuals:
MyCorrelograms[[2]]<-ggplot(spatial_correlogram, aes(x = Distance, y = MoranI_coef)) +
  geom_point(size=3)+
  geom_line()+
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  ylim(c(-1, 1)) +
  ylab("") +
  xlab("Geographical distance (km)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        axis.text.y = element_text(hjust=0.5, vjust=0.5, angle=0, size=10),
        axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=10),
        axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=12, colour="black", face="bold", margin=margin(t=0, r=5, b=0, l=0)),
        axis.title.x=element_text(size=12, colour="black", face="bold", margin=margin(t=5, r=0, b=0, l=0)),
        legend.position="none",
        plot.background=element_blank(),
        panel.spacing=unit(0,"null"))

# Create the multipanel plot with both correlograms:
Multipanel_plot<-ggpubr::ggarrange(MyCorrelograms[[1]], MyCorrelograms[[2]],
                                   labels=c("a)", "b)"), align="hv",
                                   font.label=list(size=12, colour = "black"), ncol=2, nrow=1); Multipanel_plot

# Save to disk:
ggsave("FigureS2.png", plot=Multipanel_plot, width=10, height=4, units="in", bg = "transparent")
ggsave("FigureS2.pdf", plot=Multipanel_plot, width=10, height=4, units="in", bg = "transparent")

#####

# STEP 9 - Inspect sex- and age-dependent effects on per species autotomy frequency.
##############################################################################################################
# STEP 9 - Inspect sex- and age-dependent effects on per species autotomy frequency.
rm(list=ls())

# Load the Data and the output of GLMM:
load("GLMM_output_group_StdSVL.Rdata")
rm(final_model, tab_1, var_selection_output, variable_importance)

# Check if sexual size dimorphism is related to autotomy frequency:
summary_table<-Data %>% 
  dplyr::group_by(Suborder, SpeciesName_ReptileDatabase, Sex) %>% 
  dplyr::summarise(MedianSVL=median(SVL_mm),
                   MeanSVL=mean(SVL_mm),
                   NAutotomy=sum(AutotomyPresence),
                   NTotal=length(SpeciesName_ReptileDatabase))

# Get the median size and the number of autotomy records per species and sex:
summary_table$AutotomyFreq<-summary_table$NAutotomy/summary_table$NTotal
summary_table<-droplevels(summary_table[summary_table$Sex!="Unknown",])
male_df<-droplevels(summary_table[summary_table$Sex=="Male", c("SpeciesName_ReptileDatabase", "MedianSVL", "MeanSVL", "AutotomyFreq", "NTotal", "Suborder")])
female_df<-droplevels(summary_table[summary_table$Sex=="Female", c("SpeciesName_ReptileDatabase", "MedianSVL", "MeanSVL", "AutotomyFreq", "NTotal")])
names(male_df)<-c("Species", "MaleMedianSVL", "MaleMeanSVL", "MaleAutotomyFreq", "NTotalMales", "Suborder")
names(female_df)<-c("Species", "FemaleMedianSVL", "FemaleMeanSVL", "FemaleAutotomyFreq", "NTotalFemales")

# Combine the datasets:
SSD_Table<-merge(male_df, female_df, by="Species", all=T)
SSD_Table$NTotal<-SSD_Table$NTotalFemales+SSD_Table$NTotalMales
SSD_Table$PropFemales<-SSD_Table$NTotalFemales/SSD_Table$NTotal
rm(male_df, female_df, summary_table)

# Compute the autotomy frequency for each species and sex:
SSD_Table$MedianSSD<-SSD_Table$FemaleMedianSVL/SSD_Table$MaleMedianSVL
SSD_Table$MeanSSD<-SSD_Table$FemaleMeanSVL/SSD_Table$MaleMeanSVL
SSD_Table$AutotomyFreqSD<-SSD_Table$FemaleAutotomyFreq/SSD_Table$MaleAutotomyFreq

# Remove undefined and NA values (for species with only one sampled sex or one sex with autotomy):
SSD_Table<-SSD_Table[!is.infinite(SSD_Table$AutotomyFreqSD),]
SSD_Table<-SSD_Table[SSD_Table$AutotomyFreqSD!=0,]
SSD_Table<-SSD_Table[complete.cases(SSD_Table),]

# Perform tests of sexual dimorphism the 
SSD_is_present<-list()
for(i in 1:nlevels(Data$SpeciesName_ReptileDatabase)){
  
  # Create a data subset without unsexed specimens:
  Data_subset<-droplevels(Data[Data$SpeciesName_ReptileDatabase==levels(Data$SpeciesName_ReptileDatabase)[i] & Data$Sex!="Unknown",])
  
  if(nrow(Data_subset)>=29){ # at least 29 specimens
    SSD_test<-agricolae::kruskal(y=Data_subset$SVL_mm, trt=Data_subset$Sex, alpha = 0.05, p.adj="bonferroni", group=TRUE) 
    SSD_is_present[[i]]<-data.frame(Species=levels(Data$SpeciesName_ReptileDatabase)[i], pvalue=SSD_test$statistics[1,3])
  }
  
}
SSD_is_present<-rbindlist(SSD_is_present)

# Include information on significance of SSD:
SSD_Table<-merge(SSD_Table, SSD_is_present, by="Species", all.x=T)
SSD_Table<-SSD_Table %>% mutate(SSD_Present = cut(pvalue, breaks=c(0, 0.05, 1), include.lowest=TRUE))
SSD_Table$SSD_Present<-factor(SSD_Table$SSD_Present, 
                              levels = c(levels(SSD_Table$SSD_Present)[1], levels(SSD_Table$SSD_Present)[2]),
                              labels = c("with SSD", "without SSD"))
SSD_Table$GroupSSD<-paste0(SSD_Table$Suborder, " ", SSD_Table$SSD_Present)
SSD_Table<-SSD_Table[order(SSD_Table$Suborder, decreasing=F),]

# Plot difference in autotomy frequency between sexes against female proportion in the sample:
place_label<-function(x, perc.dist, origin){
  if(origin=="max"){ return(max(x, na.rm=T) - (perc.dist*(max(x, na.rm=T)-min(x, na.rm=T))))}
  if(origin=="min"){ return(min(x, na.rm=T) + (perc.dist*(max(x, na.rm=T)-min(x, na.rm=T))))}
}

# Create the plot of sexual differences in size and autotomy frequency: 
SSD_Table$Suborder<-factor(SSD_Table$Suborder, levels=c("Amphisbaenia", "Serpentes"), labels=c("Amphisbaenians", "Snakes"))
SSD_plot<-list()
SSD_plot[[1]]<-ggplot(data=SSD_Table, aes(y=log10(AutotomyFreqSD), x=MedianSSD, group=Suborder, fill=Suborder)) + 
  geom_hline(yintercept=0, linetype="dashed", color="gray70") +
  geom_vline(xintercept=1, linetype="dashed", color="gray70") +
  stat_smooth(method="lm", alpha=0.2, aes(fill=Suborder, color=Suborder), show.legend=FALSE) + 
  geom_point(aes(fill=Suborder, shape=Suborder, size=Suborder), color="Black") +
  scale_shape_manual(values=c(22,21)) +
  scale_size_manual(values=c(4,4)) +
  scale_fill_manual(values=c("#d6604d", "#4393c3")) +
  scale_colour_manual(values=c("#d6604d", "#4393c3")) +
  scale_y_continuous(limits=c(
    place_label(log10(SSD_Table$AutotomyFreqSD), -0.1, "min"),
    place_label(log10(SSD_Table$AutotomyFreqSD), -0.1, "max"))) +
  labs(x="Sexual Size Dimorphism", y="") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=12, face="bold"),
        axis.line = element_line(colour="black"),
        axis.text = element_text(size=10),
        legend.justification = c(1, 0),
        legend.position=c(.9, .15),
        legend.title=element_blank(),
        legend.spacing.y = unit(1, 'cm'),
        legend.text=element_text(size=10),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.key = element_blank(),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  
  annotate("text", label=bquote('All sexed specimens'), x=place_label(SSD_Table$MedianSSD, 0.5, "min"), y=1.05, color="gray30", size=4, hjust=.5)+
  annotate("text", label=bquote('More autotomy in females'), x=1.8, y=0.05, color="gray50", size=3, hjust=1)+
  annotate("text", label=bquote('More autotomy in males'), x=1.8, y=-0.05, color="gray50", size=3, hjust=1) +
  annotate("text", label=bquote('Larger females'), x=1.015, y=place_label(log10(SSD_Table$AutotomyFreqSD), -0.1, "min"), color="gray50", size=3, hjust=0)+
  annotate("text", label=bquote('Larger males'), x=0.985, y=place_label(log10(SSD_Table$AutotomyFreqSD), -0.1, "min"), color="gray50", size=3, hjust=1) +
  ggpubr::stat_cor(aes(color=Suborder), size=3, method = "pearson", digits=3, hjust=1,
                   label.y=c(
                     place_label(log10(SSD_Table$AutotomyFreqSD), -0.04, "min"),
                     place_label(log10(SSD_Table$AutotomyFreqSD), -0.1, "min")
                   ),
                   label.x=c(
                     place_label(SSD_Table$MedianSSD, 0.0, "max"),
                     place_label(SSD_Table$MedianSSD, 0.0, "max")
                   )
  ); SSD_plot[[1]]

SSD_plot[[2]]<-ggplot(data=SSD_Table, aes(y=log10(AutotomyFreqSD), x=PropFemales, group=Suborder, fill=Suborder)) + 
  geom_hline(yintercept=0, linetype="dashed", color="gray70") +
  geom_vline(xintercept=0.5, linetype="dashed", color="gray70") +
  stat_smooth(method="lm", alpha=0.2, aes(fill=Suborder, color=Suborder), show.legend=FALSE) + 
  geom_point(aes(fill=Suborder, shape=Suborder, size=Suborder), color="Black") +
  scale_shape_manual(values=c(22,21)) +
  scale_size_manual(values=c(4,4)) +
  scale_fill_manual(values=c("#d6604d", "#4393c3")) +
  scale_colour_manual(values=c("#d6604d", "#4393c3")) +
  scale_y_continuous(limits=c(
    place_label(log10(SSD_Table$AutotomyFreqSD), -0.1, "min"),
    place_label(log10(SSD_Table$AutotomyFreqSD), -0.1, "max"))) +
  labs(x="Proportion of female in the sample", y="") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=12, face="bold"),
        axis.line = element_line(colour="black"),
        axis.text = element_text(size=10),
        legend.justification = c(0, 0),
        legend.position="none") +
  
  annotate("text", label=bquote('All sexed specimens'), x=place_label(SSD_Table$PropFemales, 0.5, "min"), y=1.05, color="gray30", size=4, hjust=.5)+
  annotate("text", label=bquote('More females'), x=0.51, y=place_label(log10(SSD_Table$AutotomyFreqSD), -0.1, "min"), color="gray50", size=3, hjust=0)+
  annotate("text", label=bquote('More males'), x=0.49, y=place_label(log10(SSD_Table$AutotomyFreqSD), -0.1, "min"), color="gray50", size=3, hjust=1) +
  ggpubr::stat_cor(aes(color=SSD_Table$Suborder), size=3, method = "pearson", digits=3, hjust=0,
                   label.y=c(
                     place_label(log10(SSD_Table$AutotomyFreqSD), -0.04, "min"),
                     place_label(log10(SSD_Table$AutotomyFreqSD), -0.1, "min")
                   ),
                   label.x=c(
                     place_label(SSD_Table$PropFemales, 0.0, "min"),
                     place_label(SSD_Table$PropFemales, 0.0, "min")
                   )
  )

# Same as above, but using only adults:
summary_table<-Data %>% dplyr::filter(LifeStageBinary==1) %>%
  dplyr::group_by(Suborder, SpeciesName_ReptileDatabase, Sex) %>% 
  dplyr::summarise(MedianSVL=median(SVL_mm),
                   MeanSVL=mean(SVL_mm),
                   NAutotomy=sum(AutotomyPresence),
                   NTotal=length(SpeciesName_ReptileDatabase))
summary_table$AutotomyFreq<-summary_table$NAutotomy/summary_table$NTotal
summary_table<-droplevels(summary_table[summary_table$Sex!="Unknown",])
male_df<-droplevels(summary_table[summary_table$Sex=="Male", c("SpeciesName_ReptileDatabase", "MedianSVL", "MeanSVL", "AutotomyFreq", "NTotal", "Suborder")])
female_df<-droplevels(summary_table[summary_table$Sex=="Female", c("SpeciesName_ReptileDatabase", "MedianSVL", "MeanSVL", "AutotomyFreq", "NTotal")])
names(male_df)<-c("Species", "MaleMedianSVL", "MaleMeanSVL", "MaleAutotomyFreq", "NTotalMales", "Suborder")
names(female_df)<-c("Species", "FemaleMedianSVL", "FemaleMeanSVL", "FemaleAutotomyFreq", "NTotalFemales")
SSD_Table2<-merge(male_df, female_df, by="Species", all=T)
SSD_Table2$NTotal<-SSD_Table2$NTotalFemales+SSD_Table2$NTotalMales
SSD_Table2$PropFemales<-SSD_Table2$NTotalFemales/SSD_Table2$NTotal
rm(male_df, female_df, summary_table)
SSD_Table2$MedianSSD<-SSD_Table2$FemaleMedianSVL/SSD_Table2$MaleMedianSVL
SSD_Table2$MeanSSD<-SSD_Table2$FemaleMeanSVL/SSD_Table2$MaleMeanSVL
SSD_Table2$AutotomyFreqSD<-SSD_Table2$FemaleAutotomyFreq/SSD_Table2$MaleAutotomyFreq
SSD_Table2<-SSD_Table2[!is.infinite(SSD_Table2$AutotomyFreqSD),]
SSD_Table2<-SSD_Table2[SSD_Table2$AutotomyFreqSD!=0,]
SSD_Table2<-SSD_Table2[complete.cases(SSD_Table2),]
SSD_Table2<-merge(SSD_Table2, SSD_is_present, by="Species", all.x=T)
SSD_Table2<-SSD_Table2 %>% mutate(SSD_Present = cut(pvalue, breaks=c(0, 0.05, 1), include.lowest=TRUE))
SSD_Table2$SSD_Present<-factor(SSD_Table2$SSD_Present, 
                               levels = c(levels(SSD_Table2$SSD_Present)[1], levels(SSD_Table2$SSD_Present)[2]),
                               labels = c("with SSD", "without SSD"))
SSD_Table2$GroupSSD<-paste0(SSD_Table2$Suborder, " ", SSD_Table2$SSD_Present)
SSD_Table2<-SSD_Table2[order(SSD_Table2$Suborder, decreasing=F),]
SSD_Table2$Suborder<-factor(SSD_Table2$Suborder, levels=c("Amphisbaenia", "Serpentes"), labels=c("Amphisbaenians", "Snakes"))

SSD_plot[[3]]<-ggplot(data=SSD_Table2, aes(y=log10(AutotomyFreqSD), x=MedianSSD, group=Suborder, fill=Suborder)) + 
  geom_hline(yintercept=0, linetype="dashed", color="gray70") +
  geom_vline(xintercept=1, linetype="dashed", color="gray70") +
  stat_smooth(method="lm", alpha=0.2, aes(fill=Suborder, color=Suborder), show.legend=FALSE) + 
  geom_point(aes(fill=Suborder, shape=Suborder, size=Suborder), color="Black") +
  scale_shape_manual(values=c(22,21)) +
  scale_size_manual(values=c(4,4)) +
  scale_fill_manual(values=c("#d6604d", "#4393c3")) +
  scale_colour_manual(values=c("#d6604d", "#4393c3")) +
  scale_y_continuous(limits=c(
    place_label(log10(SSD_Table2$AutotomyFreqSD), -0.1, "min"),
    place_label(log10(SSD_Table2$AutotomyFreqSD), -0.1, "max"))) +
  labs(x="Female to male median size ratio", y="") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=12, face="bold"),
        axis.line = element_line(colour="black"),
        axis.text = element_text(size=10),
        legend.justification = c(1, 0),
        legend.position=c(.9, .2),
        legend.title=element_blank(),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.text=element_text(size=10),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.key = element_blank(),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  
  annotate("text", label=bquote('Only sexed adults'), x=place_label(SSD_Table2$MedianSSD, 0.5, "min"), y=1.03, color="gray30", size=4, hjust=.5)+
  annotate("text", label=bquote('More autotomy in females'), x=1.55, y=0.05, color="gray50", size=3, hjust=1)+
  annotate("text", label=bquote('More autotomy in males'), x=1.55, y=-0.05, color="gray50", size=3, hjust=1) +
  annotate("text", label=bquote('Larger females'), x=1.015, y=place_label(log10(SSD_Table2$AutotomyFreqSD), -0.1, "min"), color="gray50", size=3, hjust=0)+
  annotate("text", label=bquote('Larger males'), x=0.985, y=place_label(log10(SSD_Table2$AutotomyFreqSD), -0.1, "min"), color="gray50", size=3, hjust=1) +
  ggpubr::stat_cor(aes(color=SSD_Table2$Suborder), size=3, method = "pearson", digits=3, hjust=1,
                   label.y=c(
                     place_label(log10(SSD_Table2$AutotomyFreqSD), -0.04, "min"),
                     place_label(log10(SSD_Table2$AutotomyFreqSD), -0.1, "min")
                   ),
                   label.x=c(
                     place_label(SSD_Table2$MedianSSD, 0.0, "max"),
                     place_label(SSD_Table2$MedianSSD, 0.0, "max")
                   )
  )

SSD_plot[[4]]<-ggplot(data=SSD_Table2, aes(y=log10(AutotomyFreqSD), x=PropFemales, group=Suborder, fill=Suborder)) + 
  geom_hline(yintercept=0, linetype="dashed", color="gray70") +
  geom_vline(xintercept=0.5, linetype="dashed", color="gray70") +
  stat_smooth(method="lm", alpha=0.2, aes(fill=Suborder, color=Suborder), show.legend=FALSE) + 
  geom_point(aes(fill=Suborder, shape=Suborder, size=Suborder), color="Black") +
  scale_shape_manual(values=c(22,21)) +
  scale_size_manual(values=c(4,4)) +
  scale_fill_manual(values=c("#d6604d", "#4393c3")) +
  scale_colour_manual(values=c("#d6604d", "#4393c3")) +
  scale_y_continuous(limits=c(
    place_label(log10(SSD_Table2$AutotomyFreqSD), -0.1, "min"),
    place_label(log10(SSD_Table2$AutotomyFreqSD), -0.1, "max"))) +
  labs(x="Proportion of females in the sample", y="") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=12, face="bold"),
        axis.line = element_line(colour="black"),
        axis.text = element_text(size=10),
        legend.justification = c(0, 0),
        legend.position="none") +
  
  annotate("text", label=bquote('Only sexed adults'), x=place_label(SSD_Table2$PropFemales, 0.5, "min"), y=1.03, color="gray30", size=4, hjust=.5)+
  annotate("text", label=bquote('More females'), x=0.515, y=place_label(log10(SSD_Table2$AutotomyFreqSD), -0.1, "min"), color="gray50", size=3, hjust=0)+
  annotate("text", label=bquote('More males'), x=0.485, y=place_label(log10(SSD_Table2$AutotomyFreqSD), -0.1, "min"), color="gray50", size=3, hjust=1) +
  ggpubr::stat_cor(aes(color=SSD_Table2$Suborder), size=3, method = "pearson", digits=3, hjust=1,
                   label.y=c(
                     place_label(log10(SSD_Table2$AutotomyFreqSD), -0.04, "min"),
                     place_label(log10(SSD_Table2$AutotomyFreqSD), -0.10, "min")
                   ),
                   label.x=c(
                     place_label(SSD_Table2$PropFemales, 0.0, "max"),
                     place_label(SSD_Table2$PropFemales, 0.0, "max")
                   )
  )

# Build the multipanel plot:
Multipanel_plot<-ggpubr::ggarrange(SSD_plot[[1]] + ylab("LR (Autotomy Freq. Females to Males)") + xlab(""),
                                   SSD_plot[[2]] + xlab(""), 
                                   SSD_plot[[3]] + theme(legend.position="none") + ylab("LR (Autotomy Freq. Females to Males)"), 
                                   SSD_plot[[4]], 
                                   labels=c("a)", "b)", "c)", "d)"), align="hv",
                                   font.label=list(size=12, colour = "black"), ncol=2, nrow=2); Multipanel_plot

# Export to disk:
ggsave("Figure5.png", plot=Multipanel_plot, width=10, height=8, units="in", bg = "transparent")
ggsave("Figure5.pdf", plot=Multipanel_plot, width=10, height=8, units="in", bg = "transparent")


# Same as above, but for adults versus juveniles:
summary_table<-Data %>% dplyr::group_by(Suborder, SpeciesName_ReptileDatabase, LifeStage) %>% dplyr::summarise(NAutotomy=sum(AutotomyPresence),
                                                                                                               NTotal=length(SpeciesName_ReptileDatabase))
summary_table$AutotomyFreq<-summary_table$NAutotomy/summary_table$NTotal
adult_df<-droplevels(summary_table[summary_table$LifeStage=="Adult", c("SpeciesName_ReptileDatabase", "AutotomyFreq", "NTotal", "Suborder")])
juvenile_df<-droplevels(summary_table[summary_table$LifeStage=="Juvenile", c("SpeciesName_ReptileDatabase", "AutotomyFreq", "NTotal")])
names(adult_df)<-c("Species", "AdultAutotomyFreq", "NTotalAdults", "Suborder")
names(juvenile_df)<-c("Species", "JuvenileAutotomyFreq", "NTotalJuveniles")
spp_autotomy<-Data %>% dplyr::group_by(SpeciesName_ReptileDatabase) %>% dplyr::summarise(NAutotomy=sum(AutotomyPresence),NTotal=length(SpeciesName_ReptileDatabase))
spp_autotomy$AutotomyFreq<-spp_autotomy$NAutotomy/spp_autotomy$NTotal
LS_Table<-merge(adult_df, juvenile_df, by="Species", all=T)
LS_Table<-merge(LS_Table, spp_autotomy, by.x="Species", by.y="SpeciesName_ReptileDatabase", all.x=T)
LS_Table$NTotal<-LS_Table$NTotalJuveniles+LS_Table$NTotalAdults
LS_Table$PropAdults<-LS_Table$NTotalAdults/LS_Table$NTotal
LS_Table$AutotomyFreqLS<-LS_Table$AdultAutotomyFreq/LS_Table$JuvenileAutotomyFreq
rm(adult_df, juvenile_df, summary_table, spp_autotomy)
LS_Table<-LS_Table[!is.infinite(LS_Table$AutotomyFreqLS),]
LS_Table<-LS_Table[LS_Table$AutotomyFreqLS!=0,]
LS_Table<-LS_Table[complete.cases(LS_Table),]
LS_Table$PropJuveniles<-LS_Table$NTotalJuveniles/LS_Table$NTotal
LS_Table$Suborder<-factor(LS_Table$Suborder, levels=c("Amphisbaenia", "Serpentes"), labels=c("Amphisbaenians", "Snakes"))

MyPlot<-list()
MyPlot[[1]]<-ggplot(data=LS_Table, aes(y=PropJuveniles, x=AutotomyFreq, group=Suborder, fill=Suborder)) + 
  stat_smooth(method="lm", alpha=0.2, aes(fill=Suborder, color=Suborder), show.legend=FALSE) + 
  geom_point(aes(fill=Suborder, shape=Suborder, size=Suborder), color="Black") +
  scale_shape_manual(values=c(22,21)) +
  scale_size_manual(values=c(4,4)) +
  scale_fill_manual(values=c("#d6604d", "#4393c3")) +
  scale_colour_manual(values=c("#d6604d", "#4393c3")) +
  scale_y_continuous(limits=c(
    place_label(LS_Table$JuvenileAutotomyFreq, -0.6, "min"),
    place_label(LS_Table$JuvenileAutotomyFreq, -0.1, "max"))) +
  labs(x="Proportion of juveniles in the sample", y="Juvenile Autotomy Frequency") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=12, face="bold"),
        axis.line = element_line(colour="black"),
        axis.text = element_text(size=10),
        legend.justification = c(0, 0),
        legend.position=c(.05, .15),
        legend.title=element_blank(),
        legend.spacing.y = unit(1, 'cm'),
        legend.text=element_text(size=10),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  
  ggpubr::stat_cor(aes(color=Suborder), size=3, method = "pearson", digits=3, hjust=0,
                   label.y=c(
                     place_label(LS_Table$JuvenileAutotomyFreq, -0.44, "min"),
                     place_label(LS_Table$JuvenileAutotomyFreq, -0.51, "min")
                   ),
                   label.x=c(
                     place_label(LS_Table$PropJuveniles, 0.0, "min"),
                     place_label(LS_Table$PropJuveniles, 0.0, "min")
                   )
  ); MyPlot[[1]]


MyPlot[[1]]<-ggplot(data=LS_Table, aes(y=JuvenileAutotomyFreq, x=PropJuveniles, group=Suborder, fill=Suborder)) + 
  stat_smooth(method="lm", alpha=0.2, aes(fill=Suborder, color=Suborder), show.legend=FALSE) + 
  geom_point(aes(fill=Suborder, shape=Suborder, size=Suborder), color="Black") +
  scale_shape_manual(values=c(22,21)) +
  scale_size_manual(values=c(4,4)) +
  scale_fill_manual(values=c("#d6604d", "#4393c3")) +
  scale_colour_manual(values=c("#d6604d", "#4393c3")) +
  scale_y_continuous(limits=c(
    place_label(LS_Table$JuvenileAutotomyFreq, -0.6, "min"),
    place_label(LS_Table$JuvenileAutotomyFreq, -0.1, "max"))) +
  labs(x="Proportion of juveniles in the sample", y="Juvenile Autotomy Frequency") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=12, face="bold"),
        axis.line = element_line(colour="black"),
        axis.text = element_text(size=10),
        legend.justification = c(0, 0),
        legend.position=c(.05, .15),
        legend.title=element_blank(),
        legend.spacing.y = unit(1, 'cm'),
        legend.text=element_text(size=10),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  
  ggpubr::stat_cor(aes(color=Suborder), size=3, method = "pearson", digits=3, hjust=0,
                   label.y=c(
                     place_label(LS_Table$JuvenileAutotomyFreq, -0.44, "min"),
                     place_label(LS_Table$JuvenileAutotomyFreq, -0.51, "min")
                   ),
                   label.x=c(
                     place_label(LS_Table$PropJuveniles, 0.0, "min"),
                     place_label(LS_Table$PropJuveniles, 0.0, "min")
                   )
  ); MyPlot[[1]]


MyPlot[[2]]<-ggplot(data=LS_Table, aes(y=log10(AutotomyFreqLS), x=PropAdults, group=Suborder, fill=Suborder)) + 
  stat_smooth(method="lm", alpha=0.2, aes(fill=Suborder, color=Suborder), show.legend=FALSE) + 
  geom_hline(yintercept=0, linetype="dashed", color="gray70") +
  geom_vline(xintercept=0.5, linetype="dashed", color="gray70") +
  geom_point(aes(fill=Suborder, shape=Suborder, size=Suborder), color="Black") +
  scale_shape_manual(values=c(22,21)) +
  scale_size_manual(values=c(4,4)) +
  scale_fill_manual(values=c("#d6604d", "#4393c3")) +
  scale_colour_manual(values=c("#d6604d", "#4393c3")) +
  scale_y_continuous(limits=c(
    place_label(log10(LS_Table$AutotomyFreqLS), -0.1, "min"),
    place_label(log10(LS_Table$AutotomyFreqLS), -0.1, "max"))) +
  labs(x="Proportion of adults in the sample", y="LR (Autotomy Freq. Adults to Juv.)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=12, face="bold"),
        axis.line = element_line(colour="black"),
        axis.text = element_text(size=10),
        legend.justification = c(1, 0),
        legend.position="none",
        legend.title=element_blank(),
        legend.spacing.y = unit(1, 'cm'),
        legend.background = element_blank(),
        legend.text=element_text(size=10),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.key = element_blank(),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  
  annotate("text", label=bquote('More autotomy in adults'), x=0.9, y=0.05, color="gray50", size=3, hjust=1)+
  annotate("text", label=bquote('More autotomy in juveniles'), x=0.9, y=-0.05, color="gray50", size=3, hjust=1) +
  annotate("text", label=bquote('More adults'), x=0.515, y=place_label(log10(LS_Table$AutotomyFreqLS), -0.1, "min"), color="gray50", size=3, hjust=0)+
  annotate("text", label=bquote('More juveniles'), x=0.485, y=place_label(log10(LS_Table$AutotomyFreqLS), -0.1, "min"), color="gray50", size=3, hjust=1) +
  ggpubr::stat_cor(aes(color=LS_Table$Suborder), size=3, method = "pearson", digits=3, hjust=1,
                   label.y=c(
                     place_label(log10(LS_Table$AutotomyFreqLS), -0.04, "min"),
                     place_label(log10(LS_Table$AutotomyFreqLS), -0.10, "min")
                   ),
                   label.x=c(
                     place_label(LS_Table$PropAdults, 0.0, "max"),
                     place_label(LS_Table$PropAdults, 0.0, "max")
                   )
  ); MyPlot[[2]]

# Build the multipanel plot:
Multipanel_plot<-ggpubr::ggarrange(MyPlot[[1]],
                                   MyPlot[[2]], 
                                   labels=c("a)", "b)"), align="hv",
                                   font.label=list(size=12, colour = "black"), ncol=2, nrow=1); Multipanel_plot

# Export to disk:
ggsave("FigureS3.png", plot=Multipanel_plot, width=10, height=4, units="in", bg = "transparent")
ggsave("FigureS3.pdf", plot=Multipanel_plot, width=10, height=4, units="in", bg = "transparent")

#####

# STEP 10 - Check sexual size dimorphism across species
##############################################################################################################
# STEP 10 - Check sexual size dimorphism across species
rm(list=ls())

# Load the Data and the output of GLMM:
load("GLMM_output_group_StdSVL.Rdata")
rm(final_model, tab_1, var_selection_output, variable_importance)

# Create a boxplot function to plot the median SVL between males and females of each species:       
SVL_Boxplot<-function(x_var, y_var, xAxis_lab){
  
  # Just for intra-function tests:
  x_var<-as.factor(x_var)
  
  # Get the x and y axis position of the posthoc small caption letters
  value_max = as.data.frame(cbind(x_var, y_var)) %>% 
    dplyr::group_by(x_var) %>%
    dplyr::summarize(max_value = max(as.numeric(as.character(y_var))*1.01),
                     label = length(y_var))
  
  # Relabel the levels of x_var:
  value_max$x_var<-factor(value_max$x_var, levels = c(1, 2), labels = c(0,1))
  
  # Get the multiple comparisons among groups:
  posthoctest<-agricolae::kruskal(y=y_var, trt=x_var, alpha = 0.05, p.adj="bonferroni", group=TRUE) 
  sig.letters <- posthoctest$groups[order(row.names(posthoctest$groups)), ]
  
  # Build the plot:
  ggplot(data = data.frame(x_var=x_var, y_var=y_var), aes(x=x_var, y=y_var)) +
    
    # Add geom objects, the number of observations and error bars:
    geom_boxplot(aes(fill=x_var), colour="black", outlier.alpha=0.25, notch=F, na.rm=T) + 
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.3) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width=0.3) +
    geom_text(data=value_max, aes(x=x_var, 
                                  y=(min(y_var, na.rm=T)-(max(y_var, na.rm=T)-min(y_var, na.rm=T))*.1),
                                  label=label), size=3, vjust=0) +
    
    # Add small caption letters for the Kruskal-Wallis test:
    geom_text(data=value_max, aes(x=x_var, y=max_value, label=sig.letters$groups), 
              position=position_dodge(width = 2), size=3, vjust=0) +
    
    # Set the scales:
    scale_x_discrete(breaks=c(0,1), labels=c("Male", "Female")) +
    scale_y_continuous(limits=c(
      (min(y_var, na.rm=T)-(max(y_var, na.rm=T)-min(y_var, na.rm=T))*.15), # min value minus 10% less than the y_var range
      (max(y_var, na.rm=T)+(max(y_var, na.rm=T)-min(y_var, na.rm=T))*.1) # max value plus 10% less than the y_var range
    )) +
    
    # Other aesthetics:
    labs(x=xAxis_lab, y="") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          panel.spacing = unit(.2, "lines"), 
          axis.title = element_text(size=8, colour="black"),
          axis.line = element_line(colour="black"),
          axis.text = element_text(size=7, colour="black"),
          axis.text.x = element_text(size=8, colour="black"),
          axis.title.x = element_text(margin=margin(t=2, r=0, b=0, l=0)),
          axis.title.y = element_text(margin=margin(t=0, r=2, b=0, l=0)),
          legend.position="none",
          plot.margin=unit(c(0.5, 0, 0.2, 0), "cm"), # top, right, bottom, left
          plot.background=element_rect(fill="transparent", colour=NA),
          strip.background = element_blank(),
          strip.placement = "none",
          strip.text.x = element_text(size=8, color="black"))
}

# Build the boxplots of sexual size dimorphism for each species:
MyBoxplot<-list()
SSD_is_present<-list()
for(i in 1:nlevels(Data$SpeciesName_ReptileDatabase)){
  
  # Create a data subset without unsexed specimens:
  Data_subset<-droplevels(Data[Data$SpeciesName_ReptileDatabase==levels(Data$SpeciesName_ReptileDatabase)[i] & Data$Sex!="Unknown",])
  
  # Perform Kruskal–Wallis tests of differences in median SVL between male and female specimens:
  if(nrow(Data_subset)>=25){ # at least 25 specimens
    SSD_test<-agricolae::kruskal(y=Data_subset$SVL_mm, trt=Data_subset$Sex, alpha = 0.05, p.adj="bonferroni", group=TRUE) 
    SSD_is_present[[i]]<-data.frame(Species=levels(Data$SpeciesName_ReptileDatabase)[i], pvalue=SSD_test$statistics[1,3])
  }
  
  # Build the boxplots for amphisbaenians:
  if(nrow(Data_subset)>=25 & Data[Data$SpeciesName_ReptileDatabase==levels(Data$SpeciesName_ReptileDatabase)[i],]$Suborder[1]=="Amphisbaenia"){
    
    xAxis_lab<-levels(Data$SpeciesName_ReptileDatabase)[i]
    nlevels(Data_subset$Sex)
    
    MyBoxplot[[i]]<-SVL_Boxplot(Data_subset$SexBinary, Data_subset$SVL_mm, xAxis_lab) +
      scale_fill_manual(values=c("#f4a582", "#b2182b"))
  }
  
  # Build the boxplots for snakes:
  if(nrow(Data_subset)>=25 & Data[Data$SpeciesName_ReptileDatabase==levels(Data$SpeciesName_ReptileDatabase)[i],]$Suborder[1]=="Serpentes"){
    
    xAxis_lab<-levels(Data$SpeciesName_ReptileDatabase)[i]
    MyBoxplot[[i]]<-SVL_Boxplot(Data_subset$SexBinary, Data_subset$SVL_mm, xAxis_lab) +
      scale_fill_manual(values=c("#92c5de", "#2166ac"))
  }
}
names(SSD_is_present)<-levels(Data$SpeciesName_ReptileDatabase)
SSD_is_present<-rbindlist(SSD_is_present)

# Remove any empty element within the list and build multipanel plot:
names(MyBoxplot)<-levels(Data$SpeciesName_ReptileDatabase)
MyBoxplot<-MyBoxplot[lengths(MyBoxplot) > 0L]
Multipanel_plot<-ggpubr::ggarrange(
  
  # Amphisbaenian species:
  MyBoxplot[[1]] + ylab("Snout Vent Length (mm)"), 
  MyBoxplot[[2]], MyBoxplot[[3]], MyBoxplot[[4]], MyBoxplot[[5]],  MyBoxplot[[6]],
  
  MyBoxplot[[7]] + ylab("Snout Vent Length (mm)"), MyBoxplot[[14]], 
  
  # Snake species:
  MyBoxplot[[8]], MyBoxplot[[9]], MyBoxplot[[10]], MyBoxplot[[11]],
  MyBoxplot[[12]] + ylab("Snout Vent Length (mm)"),
  MyBoxplot[[13]], MyBoxplot[[15]], MyBoxplot[[16]], MyBoxplot[[17]], MyBoxplot[[18]],
  
  MyBoxplot[[19]] + ylab("Snout Vent Length (mm)"), 
  MyBoxplot[[20]], MyBoxplot[[21]], MyBoxplot[[22]], MyBoxplot[[23]], MyBoxplot[[24]],
  
  MyBoxplot[[25]] + ylab("Snout Vent Length (mm)"), 
  MyBoxplot[[26]], MyBoxplot[[27]], MyBoxplot[[28]], MyBoxplot[[29]], MyBoxplot[[30]],
  
  MyBoxplot[[31]] + ylab("Snout Vent Length (mm)"), 
  MyBoxplot[[32]], MyBoxplot[[33]], MyBoxplot[[34]], MyBoxplot[[35]], MyBoxplot[[36]], 
  
  MyBoxplot[[37]] + ylab("Snout Vent Length (mm)"), 
  MyBoxplot[[38]], MyBoxplot[[39]], MyBoxplot[[40]], 
  
  align="hv", font.label=list(size=12, colour = "black"), ncol=6, nrow=7)

# Export to disk:
ggsave("FigureS4.png", plot=Multipanel_plot, width=12, height=15, units="in", bg = "transparent", limitsize=F)
ggsave("FigureS4.pdf", plot=Multipanel_plot, width=12, height=15, units="in", bg = "transparent", limitsize=F)

####