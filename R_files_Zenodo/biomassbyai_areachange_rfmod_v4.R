# Created by MRU on July 22, 2021
# To estimate changes in the area of in climatic zones and average biomass by climatic zone

# Libraries: --------------------------------------------------------------
library(raster)
library(spatialEco)
library(tidyr)
library(SpaDES)
library(ncdf4)
library(dplyr)
library(stringr)
library(tidyverse)
library(Hmisc)
library(rasterVis)
library(rgdal)
library(rgeos)
library(mapview)
library(boot)
library(rethinking)
library(bayestestR)
library(sf)
library(janitor)

# Directories: ------------------------------------------------------------

biomassdir = "/Volumes/GoogleDrive-107504600749023033793/My Drive/GitHub/tropicalbiomass/final_codedata/"
biomassdir = "H:/My Drive/GitHub/tropicalbiomass/final_codedata/"

# Tropics shape: ----------------------------------------------------------

tropsbiomassext <- extent(-180,180,-23.4, 23.4)
tropsbiomasspoly <- as(tropsbiomassext, 'SpatialPolygons')
crs(tropsbiomasspoly) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# Climate data: ----------------------------------------------

## GCMs:
### Precipitation data:
i.gcmP = stack(paste0(biomassdir, "climate_data/Global_Mgcm_P.30yr.nc"),
               varname = "rcp85mean")
### MCWD data:
i.gcmMCWD.pm = stack(paste0(biomassdir, "climate_data/Global_Mgcm_MCWDpm.30yr.nc"),
                     varname = "rcp85mean")

### CRU lower and upper Q:
i.gcmP.low = raster::stack(paste0(biomassdir, "climate_data/Global_Mgcm_P.30yr_rcp85cru_lowQ.nc"))
i.gcmMCWD.low = raster::stack(paste0(biomassdir, "climate_data/Global_Mgcm_MCWDpm.30yr_rcp85cru_lowQ.nc"))

i.gcmP.up = raster::stack(paste0(biomassdir, "climate_data/Global_Mgcm_P.30yr_rcp45cru_upQ.nc"))
i.gcmMCWD.up = raster::stack(paste0(biomassdir, "climate_data/Global_Mgcm_MCWDpm.30yr_rcp45cru_upQ.nc"))

# Biomass data: -----------------------------------------------------------

# Baccini:
biomasstrops2_nohumanwater = raster(paste0(biomassdir, "biomass_data/whrc_biomass_tropicsv2_alltrops_noESAhumanwater.tif"))
biomasstrops2_nohumanwater = crop(biomasstrops2_nohumanwater, tropsbiomasspoly)

# Xu and Saatchi:
xutrops_nohumanwater = raster(paste0(biomassdir, "biomass_data/Xu_biomass_tropics_noESAhumanwater.tif"))
xutrops_nohumanwater = crop(xutrops_nohumanwater, biomasstrops2_nohumanwater)

# ESA 05 deg:
esatrops500_nohumanwater = raster(paste0(biomassdir, "biomass_data/ESA_biomass_tropics_noESAhumanwater500m.tif"))
esatrops500_nohumanwater = crop(esatrops500_nohumanwater, biomasstrops2_nohumanwater)

# ESA water mask: ---------------------------------------------------------

esawatermask = raster(paste0(biomassdir, "ESACCI-LC-L4-LCCS-Map-300m-P1Y-2010-v2.0.7_watermask05MRU.tif"))
esawatermask = crop(esawatermask, biomasstrops2_nohumanwater)

# Climate zone from P and MCWD: ------------------------------------------------------

# Crop to the tropics:
i.gcmP.trops = crop(i.gcmP.up, biomasstrops2_nohumanwater)
i.gcmMCWD.trops = crop(i.gcmMCWD.up, biomasstrops2_nohumanwater)

AI19501979 <- i.gcmP.trops[[1]]/i.gcmMCWD.trops[[1]]
AI19501979[AI19501979 < -3.8 | AI19501979 == Inf] <- 6
AI19501979[AI19501979 >= -3.8  & AI19501979 < -1.8] <- 5
AI19501979[AI19501979 >= -1.8 & AI19501979 < -1] <- 4
AI19501979[AI19501979 >= -1 & AI19501979 < -0.25] <- 3
AI19501979[AI19501979 >= -0.25 & AI19501979 < -1/19] <- 2
AI19501979[AI19501979 >= -1/19 & AI19501979 <= 0 ] <- 1
AI19501979[AI19501979 == 6 & i.gcmP.trops[[1]] >= 1700] <- 7   # separate 6 by precipitation
plot(AI19501979)

AI19802009 <- i.gcmP.trops[[2]]/i.gcmMCWD.trops[[2]]
AI19802009[AI19802009 < -3.8 | AI19802009 == Inf] <- 6
AI19802009[AI19802009 >= -3.8  & AI19802009 < -1.8] <- 5
AI19802009[AI19802009 >= -1.8 & AI19802009 < -1] <- 4
AI19802009[AI19802009 >= -1 & AI19802009 < -0.25] <- 3
AI19802009[AI19802009 >= -0.25 & AI19802009 < -1/19] <- 2
AI19802009[AI19802009 >= -1/19 & AI19802009 <= 0 ] <- 1
AI19802009[AI19802009 == 6 & i.gcmP.trops[[2]] >= 1700] <- 7   # separate 6 by precipitation
plot(AI19802009)

AI2040 <- i.gcmP.trops[[3]]/i.gcmMCWD.trops[[3]]
AI2040[AI2040 < -3.8 | AI2040 == Inf] <- 6
AI2040[AI2040 >= -3.8  & AI2040 < -1.8] <- 5
AI2040[AI2040 >= -1.8 & AI2040 < -1] <- 4
AI2040[AI2040 >= -1 & AI2040 < -0.25] <- 3
AI2040[AI2040 >= -0.25 & AI2040 < -1/19] <- 2
AI2040[AI2040 >= -1/19 & AI2040 <= 0 ] <- 1
AI2040[AI2040 == 6 & i.gcmP.trops[[3]] >= 1700] <- 7   # separate 6 by precipitation
plot(AI2040)

AI2070 <- i.gcmP.trops[[4]]/i.gcmMCWD.trops[[4]]
AI2070[AI2070 < -3.8 | AI2070 == Inf] <- 6
AI2070[AI2070 >= -3.8  & AI2070 < -1.8] <- 5
AI2070[AI2070 >= -1.8 & AI2070 < -1] <- 4
AI2070[AI2070 >= -1 & AI2070 < -0.25] <- 3
AI2070[AI2070 >= -0.25 & AI2070 < -1/19] <- 2
AI2070[AI2070 >= -1/19 & AI2070 <= 0 ] <- 1
AI2070[AI2070 == 6 & i.gcmP.trops[[4]] >= 1700] <- 7   # separate 6 by precipitation
plot(AI2070)

AI2100 <- i.gcmP.trops[[5]]/i.gcmMCWD.trops[[5]]
AI2100[AI2100 < -3.8 | AI2100 == Inf] <- 6
AI2100[AI2100 >= -3.8  & AI2100 < -1.8] <- 5
AI2100[AI2100 >= -1.8 & AI2100 < -1] <- 4
AI2100[AI2100 >= -1 & AI2100 < -0.25] <- 3
AI2100[AI2100 >= -0.25 & AI2100 < -1/19] <- 2
AI2100[AI2100 >= -1/19 & AI2100 <= 0 ] <- 1
AI2100[AI2100 == 6 & i.gcmP.trops[[5]] >= 1700] <- 7   # separate 6 by precipitation
plot(AI2100)

# Biomass by climate zone: ----------------------------------------------------------

# Resample climate zone data to biomass resolution:
AI19802009 = crop(AI19802009, biomasstrops2_nohumanwater)
AI19802009.500m = raster::resample(AI19802009, biomasstrops2_nohumanwater, method = "ngb")
AI19802009.500m = raster::mask(AI19802009.500m, biomasstrops2_nohumanwater)
plot(AI19802009.500m)

# Create a list of rasters where each raster correspond to each climatic zone:
AI19802009.byAI = lapply(c(1:7), function(x){   
  AI19802009.thisai = AI19802009.500m
  AI19802009.thisai[AI19802009.thisai != x] = NA
  return(AI19802009.thisai)
})

# Create raster layers for each climatic zone:
biomass19802009.byai = function(whichai, biomassras, biomassdataset){
  thefilename = paste0(biomassdir, "biomass_byai_partialresults/biomassbyai", biomassdataset,"_19802009_nohumanwater/biomass", biomassdataset,"_19802009_ai_", whichai, "_nohumanwater.tif")
  print(thefilename)
  biomassthisai = raster::mask(biomassras, AI19802009.byAI[[whichai]], filename = thefilename)  #, filename = thefilename
  plot(biomassthisai)
  return(biomassthisai)
}
lapply(c(1:7), biomass19802009.byai, biomassras = biomasstrops2_nohumanwater, biomassdataset = "Bacc")
lapply(c(1:7), biomass19802009.byai, biomassras = xutrops_nohumanwater, biomassdataset = "Xu")
lapply(c(1:7), biomass19802009.byai, biomassras = esatrops500_nohumanwater, biomassdataset = "ESA")

# Load the created rasters with climatic zones for each dataset:
biomass.allaiBacc.list = as.list(stack(list.files(paste0(biomassdir, "biomass_byai_partialresults/biomassbyaiBacc_19802009_nohumanwater/"), full.names = TRUE)))
biomass.allaiXu.list = as.list(stack(list.files(paste0(biomassdir, "biomass_byai_partialresults/biomassbyaiXu_19802009_nohumanwater/"), full.names = TRUE)))
biomass.allaiESA.list = as.list(stack(list.files(paste0(biomassdir, "biomass_byai_partialresults/biomassbyaiESA_19802009_nohumanwater/"), full.names = TRUE)))

# Estimate the average, lower and upper biomass by climate zone for each dataset:
biomassavg.byai = function(biomass.allai.list, biomassdataset){
  lapply(biomass.allai.list, function(x){
    biomass.thisaidf = as.data.frame(x)
    colnames(biomass.thisaidf) = c("biomass")
    biomass.quant.thisai = quantile(biomass.thisaidf[,1], na.rm = TRUE)
    biomass.clboot.thisai = smean.cl.boot(biomass.thisaidf[,1])
    biomass.sd.thisai = sd(biomass.thisaidf[,1], na.rm = TRUE)
    biomass.thisai = data.frame(filename = names(x),
                                biomass.avg = biomass.clboot.thisai[1],
                                biomass.low = biomass.clboot.thisai[2],
                                biomass.upp = biomass.clboot.thisai[3],
                                biomass.med = biomass.quant.thisai[3],
                                biomass.lowq = biomass.quant.thisai[2],
                                biomass.upq = biomass.quant.thisai[4],
                                biomass.sd = biomass.sd.thisai,
                                pixelcount = sum(!is.na(biomass.thisaidf))) ##!! this is not gonna be for x, but for the new AI cropped data
    write.table(biomass.thisai, file= paste0(biomassdir, "biomass_byai_partialresults/biomassavg", biomassdataset,"_allaimcwd_19802009_nohumanwater.csv"), 
                append = T, row.names=F, col.names=F, sep = ",")
    return(biomass.thisai)
  })
}
## run function:
biomassdatasets <- c("Bacc", "Xu", "ESA")
biomass.allai.lists <- list(biomass.allaiBacc.list, biomass.allaiXu.list, biomass.allaiESA.list)
lapply(c(1:3), function(y){
  biomassavg.byai(biomass.allai.list = biomass.allai.lists[[y]], 
                  biomassdataset = biomassdatasets[y])
})

# Fit RF model by dataset: ------------------------------------------------

biomass_rfmod <- function(biomassdataset, biomassdataset_name){
  # Resample climate data to biomass resolution:
  mcwd1980_500m = raster::resample(i.gcmMCWD.trops[[2]], biomassdataset, method = "ngb")
  p1980_500m = raster::resample(i.gcmP.trops[[2]], biomassdataset, method = "ngb")
  mcwd500_nohumanwater = raster::mask(mcwd1980_500m, biomassdataset)
  p500_nohumanwater = raster::mask(p1980_500m, biomassdataset)
  
  # Stack all data:
  biomassclim500m_nohumanwater = stack(biomassdataset, mcwd500_nohumanwater, p500_nohumanwater)
  # Turn all data into a data frame:
  biomasstrops_df = sampleRandom(biomassclim500m_nohumanwater, 
                                 size = 50000, 
                                 na.rm = TRUE,
                                 xy = TRUE)
  colnames(biomasstrops_df) = c("lon", "lat", "biomass", "MCWD", "P")
  
  # Data partition into train and test:
  set.seed(222)
  ind <- sample(2, nrow(biomasstrops_df), replace = TRUE, prob = c(0.7, 0.3))
  biomasstrops_df = as_tibble(biomasstrops_df) %>%
    bind_cols(ind = ind)
  TrainSet <- biomasstrops_df %>% filter(ind == 1) %>% select(-ind)
  ValidSet<- biomasstrops_df %>% filter(ind == 2) %>% select(-ind)
  
  rfrg <- ranger(
    formula = biomass ~ MCWD + P,
    data = TrainSet,
    num.trees = 300,
    min.node.size = 10,
    keep.inbag = TRUE,
    quantreg = TRUE)
  predictions(rfrg)
  print(rfrg)
  saveRDS(rfrg, paste0(biomassdir, "rf_models/rfrg_", biomassdataset_name,"2.rds"))
  
  return(rfrg)
}

# Estimate changes in area by climate zone and region: -----------------------------------------

# The regions shapes:
territories = readOGR(paste0(biomassdir, "shapes/world-administrative-boundaries/world-administrative-boundaries.shp"))
unique(territories$continent)
africa <- subset(territories, continent == "Africa")
america <- subset(territories, continent == "Americas")
asia <- subset(territories, continent %in% c("Asia", "Oceania"))
theamazon = readOGR(paste0(biomassdir, "shapes/amapoly_ivb/amapoly_ivb.shp"))
mapview(theamazon)

# Function to estimate changes in area by climate zone by region:
aiarea_changes_byregion = function(allaiscomb, fileout, region, 
                                   rastert1, rastert2, rastert3, rastert4, rastert5){
  
  # Cut the rasters to the region's extent:
  raster_stack = raster::stack(rastert1, rastert2, rastert3, rastert4, rastert5)
  raster_region = raster::mask(raster_stack, region)
  
  ai_t1 = raster_region[[1]]
  ai_t1[ai_t1 != allaiscomb$ai1950] = NA
  ai_t1[!is.na(ai_t1)] = 1
  ai_t2 = raster_region[[2]]
  ai_t2[ai_t2 != allaiscomb$ai2010] = NA
  ai_t2[!is.na(ai_t2)] = 1
  ai_t3 = raster_region[[3]]
  ai_t3[ai_t3 != allaiscomb$ai2040] = NA
  ai_t3[!is.na(ai_t3)] = 1
  ai_t4 = raster_region[[4]]
  ai_t4[ai_t4 != allaiscomb$ai2070] = NA
  ai_t4[!is.na(ai_t4)] = 1
  ai_t5 = raster_region[[5]]
  ai_t5[ai_t5 != allaiscomb$ai2100] = NA
  ai_t5[!is.na(ai_t5)] = 1
  
  thisaisall = stack(ai_t1, ai_t2, ai_t3, ai_t4, ai_t5)
  thisaisall_sum = sum(thisaisall)
  print(thisaisall_sum)
  thisaisall_freq = freq(thisaisall_sum)
  print(thisaisall_freq)
  print(length(which(thisaisall_freq[,1] == 5)))
  
  if(is.na(minValue(thisaisall_sum)) & is.na(maxValue(thisaisall_sum))){
    print("all nas")
    thisais_all_df = data.frame(allaiscomb$ai1950, allaiscomb$ai2010, allaiscomb$ai2040, allaiscomb$ai2070, allaiscomb$ai2100, 0, 0)
    colnames(thisais_all_df) = c("ai1950", "ai2010", "ai2040", "ai2070", "ai2100", "npixs", "area")
  }
  else if(length(which(thisaisall_freq[,1] == 5)) > 0){
    thisaisall_common = c(thisaisall_freq[which(thisaisall_freq[,1] == 5), 2])
    thisaisall_common_ras = thisaisall_sum
    thisaisall_common_ras[thisaisall_common_ras != 5] = NA
    #plot(thisaisall_common_ras)
    thisaisall_common_area = cellStats(area(thisaisall_common_ras, na.rm=TRUE, weights=FALSE), stat = sum)
    thisais_all_df = data.frame(allaiscomb$ai1950, allaiscomb$ai2010, allaiscomb$ai2040, allaiscomb$ai2070, allaiscomb$ai2100, thisaisall_common, thisaisall_common_area)
    colnames(thisais_all_df) = c("ai1950", "ai2010", "ai2040", "ai2070", "ai2100", "npixs", "area")
  }
  else{
    thisais_all_df = data.frame(allaiscomb$ai1950, allaiscomb$ai2010, allaiscomb$ai2040, allaiscomb$ai2070, allaiscomb$ai2100, 0, 0)
    colnames(thisais_all_df) = c("ai1950", "ai2010", "ai2040", "ai2070", "ai2100", "npixs", "area")
  }
  write.table(thisais_all_df, file=paste0(biomassdir, fileout), 
              append = TRUE, row.names=FALSE, col.names=FALSE, sep = ",")
  return(thisais_all_df)
}

# Create dataframe with all the climatic zones combinations:
ais2010 = data.frame(ai1950 = rep(c(1:7), each = 7), ai2010 = rep(c(1:7), 7))
ais2040 = data.frame(ai1950 = rep(c(1:7), each = 7), ai2040 = rep(c(1:7), 7))
ais2070 = data.frame(ai1950 = rep(c(1:7), each = 7), ai2070 = rep(c(1:7), 7))
ais2100 = data.frame(ai1950 = rep(c(1:7), each = 7), ai2100 = rep(c(1:7), 7))
aisall = left_join(ais2010, ais2040, by = "ai1950") %>%
  left_join(ais2070, by = "ai1950") %>%
  left_join(ais2100, by = "ai1950")
aisall.list <- split(aisall, seq(nrow(aisall)))

# Apply water mask to the climate zones rasters:
AI19501979 = raster::mask(AI19501979, esawatermask)
AI19802009 = raster::mask(AI19802009, esawatermask)
AI2040 = raster::mask(AI2040, esawatermask)
AI2070 = raster::mask(AI2070, esawatermask)
AI2100 = raster::mask(AI2100, esawatermask)

# Run the function:
lapply(aisall.list[1:length(aisall.list)], aiarea_changes_byregion, 
       region = america, fileout = "areachange_partialresults/aiareachange_rcp45_cruup_america_noESAwater.csv", 
       rastert1 = AI19501979, rastert2 = AI19802009, rastert3 = AI2040, 
       rastert4 = AI2070, rastert5 = AI2100)
lapply(aisall.list[1:length(aisall.list)], aiarea_changes_byregion, 
       region = africa, fileout = "areachange_partialresults/aiareachange_rcp45_cruup_africa_noESAwater.csv", 
       rastert1 = AI19501979, rastert2 = AI19802009, rastert3 = AI2040, 
       rastert4 = AI2070, rastert5 = AI2100)
lapply(aisall.list[1:length(aisall.list)], aiarea_changes_byregion, 
       region = asia, fileout = "areachange_partialresults/aiareachange_rcp45_cruup_asia_noESAwater.csv", 
       rastert1 = AI19501979, rastert2 = AI19802009, rastert3 = AI2040, 
       rastert4 = AI2070, rastert5 = AI2100)


# Summarize and organize changes in area by climate zone and region: --------------------------------------

# function to organize the data biomass and climate zone area changes:
areachanges = function(aiareadata){
  aiareachange_df = read_csv(aiareadata, col_names = FALSE) %>%
    dplyr::rename("ai1950" = X1, "ai2010" = X2, "ai2040" = X3, "ai2070" = X4, "ai2100" = X5, "n_pixs" = X6, "area" = X7) %>%
    filter(n_pixs > 0)
  aiareachange_dflong = aiareachange_df %>%
    pivot_longer( cols = starts_with("ai"),
                  names_to = "timeframe",
                  names_prefix = "ai",
                  values_to = "ai",
                  values_drop_na = FALSE) %>%
    mutate(ai = as.factor(ai)) %>%
    group_by(timeframe, ai) %>%
    summarise(n_pixs = sum(n_pixs, na.rm = TRUE),
              area = sum(area, na.rm = TRUE))
  
  aiareachange_dflong_total = aiareachange_dflong %>% 
    filter(timeframe == 1950) %>%
    ungroup() %>%
    dplyr::select(n_pixs, area)
  
  aiareachange_dflong2 = aiareachange_dflong %>%
    mutate(area_perc = (area*100)/sum(aiareachange_dflong_total$area)) %>%
    mutate(label_y = cumsum(area))
  
  #write_csv(biomassarea_df, "/Volumes/GoogleDrive/Shared drives/Data/GeoDataBase/Biomass_project/biomassbycs/", thisfilename, ".csv")
  return(aiareachange_dflong2)
}

# load data:
aiarea_neotrops_file = paste0(biomassdir, "areachange_partialresults/aiareachange_rcp45_cruup_america_noESAwater.csv")
aiarea_afr_file = paste0(biomassdir, "areachange_partialresults/aiareachange_rcp45_cruup_africa_noESAwater.csv")
aiarea_asia_file = paste0(biomassdir, "areachange_partialresults/aiareachange_rcp45_cruup_asia_noESAwater.csv")

# run function (skip if the data was already cleaned and organized):
## by region:
areachange_allreg = lapply(c(aiarea_neotrops_file, aiarea_afr_file, aiarea_asia_file),
                           areachanges)
reg_names = c("Neotropics", "Africa", "Asia")
areachange_allreg2 = lapply(c(1:3), function(x){
  temp = areachange_allreg[[x]]
  temp$region = reg_names[x]
  return(temp)
})
areachange_allreg3 = bind_rows(areachange_allreg2, .id = "column_label")

write_csv(areachange_allreg3, paste0(biomassdir,  "areachange_partialresults/aiareachange_rcp45_cruup_noESAwater_summ_allregions.csv"))

## load data and reorganize:
area_avg = read_csv(paste0(biomassdir,  "areachange_partialresults/aiareachange_rcp45_noESAwater_summ_allregions.csv")) %>%
  mutate(ID = "avg2")
area_low = read_csv(paste0(biomassdir,  "areachange_partialresults/aiareachange_rcp45_crulow_noESAwater_summ_allregions.csv")) %>%
  mutate(ID = "low")
area_upp = read_csv(paste0(biomassdir,  "areachange_partialresults/aiareachange_rcp45_cruup_noESAwater_summ_allregions.csv")) %>%
  mutate(ID = "upp")

areachange = bind_rows(area_avg, area_low) %>%
  bind_rows(area_upp)

# estimate total area and carbon by region to use later:
regionarea = areachange %>% 
  filter(timeframe == 1950) %>%
  group_by(region) %>%
  summarise(totalarea = sum(area, na.rm = TRUE)/1000000) %>%
  adorn_totals("row") %>%
  mutate(region = ifelse(region == "Total", "Pantropical", region))

createtable = function(avgdf, lowdf, uppdf, whichvar){
  print(head(avgdf))
    ## calculate for the pantropics:
    avg_pant = avgdf %>%
      dplyr::select("ai", "timeframe", "area") %>%
      group_by(ai, timeframe) %>%
      summarise(area = sum(area, na.rm = TRUE)) %>%
      mutate(area = round(area/1000000, digits = 2)) %>%
      pivot_wider(names_from = timeframe, values_from = area) %>%
      mutate(region = "Pantropical") %>%
      relocate(region, .before = "ai") %>%
      dplyr::rename("X1950" = `1950`, "X1980" = `2010`, "X2010" = `2040`, "X2040" = `2070`, "X2070" = `2100`)
    low_pant = lowdf %>%
      dplyr::select("ai", "timeframe", "area") %>%
      group_by(ai, timeframe) %>%
      summarise(area = sum(area, na.rm = TRUE)) %>%
      mutate(area = round(area/1000000, digits = 2)) %>%
      pivot_wider(names_from = timeframe, values_from = area)%>%
      mutate(region = "Pantropical") %>%
      relocate(region, .before = "ai") %>%
      dplyr::rename("X1950" = `1950`, "X1980" = `2010`, "X2010" = `2040`, "X2040" = `2070`, "X2070" = `2100`)
    upp_pant = uppdf %>%
      dplyr::select("ai", "timeframe", "area") %>%
      group_by(ai, timeframe) %>%
      summarise(area = sum(area, na.rm = TRUE)) %>%
      mutate(area = round(area/1000000, digits = 2)) %>%
      pivot_wider(names_from = timeframe, values_from = area) %>%
      mutate(region = "Pantropical") %>%
      relocate(region, .before = "ai") %>%
      dplyr::rename("X1950" = `1950`, "X1980" = `2010`, "X2010" = `2040`, "X2040" = `2070`, "X2070" = `2100`)
    
    ## total by region:
    avg_totreg = avgdf %>%
      group_by(region, timeframe) %>%
      summarise(area = sum(area, na.rm = TRUE)) %>%
      mutate(area = round(area/1000000, digits = 2)) %>%
      pivot_wider(names_from = timeframe, values_from = area) %>%
      dplyr::rename("X1950" = `1950`, "X1980" = `2010`, "X2010" = `2040`, "X2040" = `2070`, "X2070" = `2100`) %>%
      left_join(regionarea, by = "region") %>%
      adorn_totals("row") %>%
      mutate(Change = round(X2070 - X1950, 1),
             Change_perc1950 = round((Change*100)/X1950, 1),
             Change_percregion = round((Change*100)/totalarea, 1)) %>%
      dplyr::select(-totalarea)
    low_totreg = lowdf %>%
      group_by(region, timeframe) %>%
      summarise(area = sum(area, na.rm = TRUE)) %>%
      mutate(area = round(area/1000000, digits = 2)) %>%
      pivot_wider(names_from = timeframe, values_from = area) %>%
      dplyr::rename("X1950" = `1950`, "X1980" = `2010`, "X2010" = `2040`, "X2040" = `2070`, "X2070" = `2100`) %>%
      left_join(regionarea, by = "region") %>%
      adorn_totals("row") %>%
      mutate(Change = round(X2070 - X1950, 1),
             Change_perc1950 = round((Change*100)/X1950, 1),
             Change_percregion = round((Change*100)/totalarea, 1)) %>%
      dplyr::select(-totalarea)
    upp_totreg = uppdf %>%
      group_by(region, timeframe) %>%
      summarise(area = sum(area, na.rm = TRUE)) %>%
      mutate(area = round(area/1000000, digits = 2)) %>%
      pivot_wider(names_from = timeframe, values_from = area) %>%
      dplyr::rename("X1950" = `1950`, "X1980" = `2010`, "X2010" = `2040`, "X2040" = `2070`, "X2070" = `2100`) %>%
      left_join(regionarea, by = "region") %>%
      adorn_totals("row") %>%
      mutate(Change = round(X2070 - X1950, 1),
             Change_perc1950 = round((Change*100)/X1950, 1),
             Change_percregion = round((Change*100)/totalarea, 1)) %>%
      dplyr::select(-totalarea)
    
    ## now organize the data by regions and add the pantropics rows:
    avg2 = avgdf %>%
      dplyr::select("region", "ai", "timeframe", "area") %>%
      mutate(area = round(area/1000000, digits = 2)) %>%
      pivot_wider(names_from = timeframe, values_from = area) %>%
      dplyr::rename("X1950" = `1950`, "X1980" = `2010`, "X2010" = `2040`, "X2040" = `2070`, "X2070" = `2100`) %>%
      bind_rows(avg_pant) %>%
      left_join(regionarea, by = "region") %>%
      mutate(Change = round(X2070 - X1950, 2),
             Change_perc1950 = round((Change*100)/X1950, 1),
             Change_percregion = round((Change*100)/totalarea, 1)) %>%
      dplyr::select(-totalarea) %>%
      bind_rows(avg_totreg)
    low2 = lowdf %>%
      dplyr::select("region", "ai", "timeframe", "area") %>%
      mutate(area = round(area/1000000, digits = 2)) %>%
      pivot_wider(names_from = timeframe, values_from = area) %>%
      dplyr::rename("X1950" = `1950`, "X1980" = `2010`, "X2010" = `2040`, "X2040" = `2070`, "X2070" = `2100`) %>%
      bind_rows(low_pant) %>%
      left_join(regionarea, by = "region") %>%
      mutate(Change = round(X2070 - X1950, 2),
             Change_perc1950 = round((Change*100)/X1950, 1),
             Change_percregion = round((Change*100)/totalarea, 1)) %>%
      dplyr::select(-totalarea) %>%
      bind_rows(low_totreg)
    upp2 = uppdf %>%
      dplyr::select("region", "ai", "timeframe", "area") %>%
      mutate(area = round(area/1000000, digits = 2)) %>%
      pivot_wider(names_from = timeframe, values_from = area) %>%
      dplyr::rename("X1950" = `1950`, "X1980" = `2010`, "X2010" = `2040`, "X2040" = `2070`, "X2070" = `2100`) %>%
      bind_rows(upp_pant) %>%
      left_join(regionarea, by = "region") %>%
      mutate(Change = round(X2070 - X1950, 2),
             Change_perc1950 = round((Change*100)/X1950, 1),
             Change_percregion = round((Change*100)/totalarea, 1)) %>%
      dplyr::select(-totalarea) %>%
      bind_rows(upp_totreg)
  
  ## now organize the table with the uncertainty stuff in parenthesis:
  final = data.frame(Region = avg2$region,
                     AI = as.factor(avg2$ai),
                     X1950_1979 = paste0(avg2$X1950, " (", low2$X1950, "-", upp2$X1950, ")"),
                     X1980_2009 = paste0(avg2$X1980, " (", low2$X1980, "-", upp2$X1980, ")"),
                     X2010_2039 = paste0(avg2$X2010, " (", low2$X2010, "-", upp2$X2010, ")"),
                     X2040_2069 = paste0(avg2$X2040, " (", low2$X2040, "-", upp2$X2040, ")"),
                     X2070_2099 = paste0(avg2$X2070, " (", low2$X2070, "-", upp2$X2070, ")"),
                     Change19502099 = paste0(avg2$Change, " (", low2$Change, "-", upp2$Change, ")"),
                     Change_perc1950 = paste0(avg2$Change_perc1950, " (", low2$Change_perc1950, "-", upp2$Change_perc1950, ")"),
                     Change_percregion = paste0(avg2$Change_percregion, " (", low2$Change_percregion, "-", upp2$Change_percregion, ")")) %>%
    mutate(AI = recode(as.factor(AI), '1' = "Hyperarid", '2' = "Arid", '3' = "Semiarid", '4' = "Dry Subhumid",
                       '5' = "Humid Savanna", '6' = "Humid Seasonal", '7' = "Humid"))
  return(final)
}

area_results = createtable(avgdf = area_avg, lowdf = area_low, uppdf = area_upp)

write_csv(area_results, file = paste0(biomassdir,  "areachange_partialresults/areachanges_rcp45_nohumanwater_noESAwater.csv"))

