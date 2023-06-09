### load all necessary csv files
#-------------------------------
# This project was run under R version 4.1.2 Bird Hippie

#### incub_resp ####
# respiration_over_50d_incubation.csv contains drift-corrected GC data of each sampling point
incub_resp <- read.csv2("respiration_over_50d_incubation.csv")
incub_resp[c(1:12,26,32)] <- lapply(incub_resp[c(1:12,26,32)], as.factor)

incub_resp$date <- as.POSIXct(incub_resp$date)  
incub_resp$date_time <- as.POSIXct(incub_resp$date_time)

# specify levels
levels(incub_resp$landuse)[levels(incub_resp$landuse)=="F"] <- "Forest"
levels(incub_resp$landuse)[levels(incub_resp$landuse)=="G"] <- "Grassland"
levels(incub_resp$landuse)[levels(incub_resp$landuse)=="CM"] <- "Cropland" 
levels(incub_resp$site)[levels(incub_resp$site)=="NO"] <- "Site MB"
levels(incub_resp$site)[levels(incub_resp$site)=="SI"] <- "Site SI"
levels(incub_resp$site)[levels(incub_resp$site)=="CD"] <- "Site CD"

# unordered version needed for lmem on RRTN over time 
incub_resp_unordered <- incub_resp

# ordered version needed for plotting
incub_resp$landuse <- ordered(incub_resp$landuse, levels = c("Forest", "Grassland", "Cropland"))


#### sample_data ####
# sample_data.csv contains cumulative respiration, CN data of each individual sample, CUE, Cmic data
sample_data <- read.csv2("sample_data.csv")

sample_data[1:12] <- lapply(sample_data[1:12], as.factor)

levels(sample_data$landuse)[levels(sample_data$landuse)=="F"] <- "Forest"
levels(sample_data$landuse)[levels(sample_data$landuse)=="G"] <- "Grassland"
levels(sample_data$landuse)[levels(sample_data$landuse)=="CM"] <- "Cropland"
levels(sample_data$site)[levels(sample_data$site)=="NO"] <- "Site MB"
levels(sample_data$site)[levels(sample_data$site)=="SI"] <- "Site SI"
levels(sample_data$site)[levels(sample_data$site)=="CD"] <- "Site CD"

sample_data$qCO2_ugC_mgCmic_h <- (sample_data$cumul_resp_incub_ngC_g_total_time/(49*24))/((sample_data$Cmic_ugC_gDW+sample_data$Cmic_V1)/2)
sample_data$cumul_resp_Corg <- (sample_data$cumul_resp_incub_ngC_g_total_time/1000000)/(sample_data$Corg_perc/100)
sample_data$Cmic_Corg_perc <- ((sample_data$Cmic_V1/1000000)/(sample_data$Corg_perc/100))*100
sample_data$Cmic.Nmic_V1 <- ((sample_data$Cmic_V1)/(sample_data$Nmic_V1))

# unordered version needed for lmem on RRTN over time 
sample_data_unordered <- sample_data

# ordered version needed for plotting
sample_data$landuse <- ordered(sample_data$landuse, levels = c("Forest", "Grassland", "Cropland"))


#### gensoil ####
# general_soil_parameters_per_sample.csv contains general soil parameter data for each individual mb sample
# 3 field replicates per plot without treatments, n = 27
gensoil <- read.csv2("general_soil_parameters_per_sample.csv")

gensoil[1:7] <- lapply(gensoil[1:7], as.factor)

gensoil$landuse <- ordered(gensoil$landuse, levels = c("F", "G", "CM"))
levels(gensoil$landuse)[levels(gensoil$landuse)=="F"] <- "Forest"
levels(gensoil$landuse)[levels(gensoil$landuse)=="G"] <- "Grassland"
levels(gensoil$landuse)[levels(gensoil$landuse)=="CM"] <- "Cropland"
levels(gensoil$site)[levels(gensoil$site)=="NO"] <- "Site MB"
levels(gensoil$site)[levels(gensoil$site)=="SI"] <- "Site SI"
levels(gensoil$site)[levels(gensoil$site)=="CD"] <- "Site CD"

gensoil$Cmic_Corg_perc <- ((gensoil$Cmic_V1_ugC_gsoil/1000000)/(gensoil$Corg_perc/100))*100
gensoil$Cmic.Nmic <- (gensoil$Cmic_V1_ugC_gsoil/gensoil$Nmic_V1_ugN_gsoil)


#### compsoil ####
# general_soil_parameters_per_sample.csv contains general soil parameter data for each plot 
# measured on one pooled composite sample of the three replicates per plot n = 9
compsoil <- read.csv2("general_soil_parameters_per_plot.csv")

compsoil[1:6] <- lapply(compsoil[1:6], as.factor)

compsoil$landuse <- ordered(compsoil$landuse, levels = c("F", "G", "CM"))
levels(compsoil$landuse)[levels(compsoil$landuse)=="F"] <- "Forest"
levels(compsoil$landuse)[levels(compsoil$landuse)=="G"] <- "Grassland"
levels(compsoil$landuse)[levels(compsoil$landuse)=="CM"] <- "Cropland"
levels(compsoil$site)[levels(compsoil$site)=="NO"] <- "Site MB"
levels(compsoil$site)[levels(compsoil$site)=="SI"] <- "Site SI"
levels(compsoil$site)[levels(compsoil$site)=="CD"] <- "Site CD"
