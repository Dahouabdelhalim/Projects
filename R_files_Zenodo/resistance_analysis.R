# resistance_analysis.R
# Test relative role of geographic vs. climatic isolation on beta diversity patterns
# Author: Jairo PatiÃ±o
# Modified by: Jun Ying Lim
# Reference: Lim et al. Climatic niche conservatism shapes the ecological assembly of Hawaiian arthropod communities
# NOTE: remember to set the working directory to the folder containing all the scripts and data

## PACKAGES  =========
# Packages necessary to install resistanceGA
# library("devtools")
# library("tinytex")
# devtools::install_github("wpeterman/ResistanceGA", build_vignettes = FALSE)
library(ResistanceGA)
library(doParallel)
library(raster) 
library(rgdal)
library(vegan)

## INPUT MODEL RUN PARAMETERS  =========
rep_no <- 80 # 14, 51, 80, 90, 92
sites <- c("Laupahoehoe", "Stainback")
taxa <- c("Araneae", "Coleoptera", "Hemiptera", "Lepidoptera", "Orthoptera", "Psocoptera")

## ITERATE ANALYSES OVER SITE / TAXA =========
for(i in 1:2){
   site <- sites[i]
   for(j in 1:6){
     taxon <- taxa[j]
     print(paste0("BEGIN ANALYSIS FOR TAXON=", taxon, " SITE=", site, " REPLICATE=", rep_no))

## IMPORT CLIAMTIC DATA =============
if(site == "Stainback"){
  annPrecip <- readRDS("annPrecip_steinb.rds") 
  annTemp <- readRDS("annTemp_steinb.rds")   
}
if(site == "Laupahoehoe"){
  annPrecip <- readRDS("annPrecip_laup.rds") 
  annTemp <- readRDS("annTemp_laup.rds")   
}

# create raster stack object
clim.stack <- stack(annPrecip, annTemp) 

## IMPORT OTU TABLES =============
otu_tab <- readRDS(paste0("data/resistance_analysis/", tolower(site), "_OTU_", taxon, "_r", rep_no, ".rds"))
plot_bray <- vegdist(t(otu_tab), method = "bray")
plot_bray_vector <- as.vector(as.dist(plot_bray, diag = FALSE, upper = FALSE))

# get locations by each sampling plot
site_data <- read.csv("siteData.csv")
site_data_subset <- site_data[site_data$site == site,]
site_coords <- site_data_subset[,c("longitude", "latitude")]
rownames(site_coords) <- site_data_subset$site.id
site_data_locales <- SpatialPoints(site_coords)

## RUN ANALYSIS =============

# set up parallel cluster
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

res.dir <- paste0("r",rep_no, "_", tolower(taxon), "_", tolower(site))
dir.create(res.dir)
GA.inputs <- GA.prep(ASCII.dir = clim.stack,
                     Results.dir = "all.comb",
                     method = "LL",
                     max.cat = 500,
                     max.cont = 500,
                     seed = 1,
                     select.trans = list("A","A"),
                     parallel = cl,
                     maxiter = 1000)


gdist.inputs <- gdist.prep(nrow(site_data_subset),
                           samples = site_data_locales,
                           response = plot_bray_vector,
                           method = 'commuteDistance') ## optimize using commute distance


# export stuff to cluster
clusterExport(cl = cl, varlist = c("GA.inputs", "gdist.inputs", "site_data_subset", "clim.stack", "plot_bray_vector", "site_data_locales"))


# run the final optimization

allComb_RESULTS.gdist <- all_comb(gdist.inputs = gdist.inputs,
                                  GA.inputs = GA.inputs,
                                  results.dir = res.dir,
                                  max.combination = 2,
                                  sample.prop = 0.75)

stopCluster(cl)

}
}