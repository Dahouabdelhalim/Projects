library(dplyr)
library(purrr)
library(tidyverse)
library(stringr)
library(SYNCSA)
require(readr)

source('01_extracting_matrices_env.R')
source('04_Single_site_convergence.R')
source('05_All_region_convergence.R')


'%nin%' <- function(x,y)!('%in%'(x,y))

B_axis<- read.table('RES_pca_individuals_0.7.7_ranks_cleaned.txt', sep="\\t", header = TRUE)
visits<-read_csv("visits_owners.csv")

#=========input the biomass of species in bromeliads ====
biomass<- read_csv("biomass_syncsa_ready.csv",
                   col_types = cols(
                     dataset_id = col_character(), 
                     species_id = col_character(), 
                     bwg_name = col_character(), 
                     bromeliad_id = col_character(),
                     biomass = col_number()
                   ))%>% 
  as.data.frame()

Biom <-biomass %>% 
  pivot_wider(id_cols = bromeliad_id, 
              names_from=c(dataset_id,species_id), 
              names_sep = "_", values_from=biomass, 
              values_fill = 0) %>% 
  mutate(bromeliad_id = as.character(bromeliad_id))


#======input raw traits, select those with enough coverage

traits<- read_csv("traits_pca_ready.csv") 

Tachet_traits <- traits %>% 
  select(species_id, matches("^[A-Z]{2}\\\\d")) 

Tra <- traits[-which(rowSums(is.na(Tachet_traits))>16),]%>% 
  mutate(species_id = as.double(species_id))

#====input the bromeliads=====

Broms.all.full<-read_csv("bromeliads_syncsa_ready.csv", 
                         col_types = cols(
                           visit_id = col_character(),
                           bromeliad_id = col_character(),
                           dataset_name = col_character(),
                           total_detritus = col_number(),
                           open.canopy = col_number(),
                           actual_water = col_number(),
                           sites = col_character(),
                           total_detritus_center= col_number(),
                           logtotal_detritus= col_number(),
                           logtotal_detritus_center= col_number(),
                           include.detritus= col_number(),
                           open.canopy_center= col_number(),
                           include.canopy= col_number(),
                           actual_water_center = col_number(),
                           logactual_water = col_number(),
                           logactual_water_center = col_number(),
                           include.water= col_number()
                         ) 
) %>% 
  as.data.frame()

#### format the trait axes #####

B_axis2 <- B_axis%>%mutate(Axis.1=as.character(Axis.1))%>%mutate(Axis.2=as.character(Axis.2))%>%
  mutate(Axis.3=as.character(Axis.3))%>%mutate(Axis.4=as.character(Axis.4))%>%
  select(species_id, Axis.1, Axis.2, Axis.3, Axis.4)%>%mutate(species_id = as.character(species_id))%>%
  mutate(Axis.1=as.numeric(Axis.1))%>%mutate(Axis.2=as.numeric(Axis.2))%>%
  mutate(Axis.3=as.numeric(Axis.3))%>%mutate(Axis.4=as.numeric(Axis.4))



###############################################
#                                             #
#       Total detritus                        #
#                                             #
###############################################

Broms.all <-Broms.all.full %>% 
  filter(include.detritus==1)

#Then we make a vector of all the regions
field.sites<-Broms.all$sites%>%unique()%>%as.vector()

##  ACROSS SITES ###

n<-1
output.site<-data.frame(TCAP.rho=numeric(n),TCDAP.rho=numeric(n), TDAP.rho=numeric(n), RE.rho=numeric(n),
                        TCAP.p=numeric(n), TCDAP.p=numeric(n),TDAP.p=numeric(n), RE.p=numeric(n),
                        best.TCAP.traits=numeric(n), best.TCAPtrait.rho=numeric(n), 
                        best.TDAP.traits=numeric(n), best.TDAPtrait.rho=numeric(n),
                        best.RAO.traits=numeric(n), best.RAOtrait.rho=numeric(n),
                        cv.env=numeric(n), min.env=numeric(n), max.env=numeric(n),
                        bestTCAP.trait.rho2=numeric(n),bestTCAP.trait.p2=numeric(n),
                        bestTDAP.trait.rho2=numeric(n),bestTDAP.trait.p2=numeric(n),
                        bestRAO.trait.rho2=numeric(n),bestRAo.trait.p2=numeric(n))

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "total_detritus_center", stock = "biomass"))
total_det.biomass.linear.axis<-output.site %>% unlist() %>% as.data.frame()
write_csv(total_det.biomass.linear.axis, "total_det_biomass_linear_axis.csv")

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "logtotal_detritus_center", ro.method = "procrustes", stock = "biomass"))
total_det.biomass.log.axis<-output.site %>% unlist() %>% as.data.frame()
write_csv(total_det.biomass.log.axis, "total_det_biomass_log_axis.csv")

##all sites but axes singly

n<-1
output.site<-data.frame(TCAP.rho=numeric(n),TCDAP.rho=numeric(n), TDAP.rho=numeric(n), RE.rho=numeric(n), TCAP.p=numeric(n), TCDAP.p=numeric(n),TDAP.p=numeric(n), RE.p=numeric(n))

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "logtotal_detritus_center", trait.name = "Axis.1",  ro.method = "procrustes", stock = "biomass"))

all_total_detritus.biomass.log.axis1<-output.site %>% 
  as.data.frame() %>% 
  mutate(sites = "All sites",
         detA1TCAP.rho =as.character(V1), 
         detA1TCAP.p =as.character(V5), 
         detA1TDAP.rho =as.character(V3), 
         detA1TDAP.p =as.character(V7), 
         detA1TCAP.rho =as.numeric(detA1TCAP.rho), 
         detA1TCAP.p =as.numeric(detA1TCAP.p), 
         detA1TDAP.rho =as.numeric(detA1TDAP.rho), 
         detA1TDAP.p =as.numeric(detA1TDAP.p)) %>% 
  select(sites, detA1TCAP.rho, detA1TCAP.p, detA1TDAP.rho, detA1TDAP.p)

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "logtotal_detritus_center", trait.name = "Axis.2", ro.method = "procrustes", stock = "biomass"))

all_total_detritus.biomass.log.axis2<-output.site %>% 
  as.data.frame() %>% 
  mutate(sites = "All sites",
         detA2TCAP.rho =as.character(V1), 
         detA2TCAP.p =as.character(V5), 
         detA2TDAP.rho =as.character(V3), 
         detA2TDAP.p =as.character(V7), 
         detA2TCAP.rho =as.numeric(detA2TCAP.rho), 
         detA2TCAP.p =as.numeric(detA2TCAP.p), 
         detA2TDAP.rho =as.numeric(detA2TDAP.rho), 
         detA2TDAP.p =as.numeric(detA2TDAP.p)) %>% 
  select(sites, detA2TCAP.rho, detA2TCAP.p, detA2TDAP.rho, detA2TDAP.p)

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "logtotal_detritus_center", trait.name = "Axis.3", ro.method = "procrustes", stock = "biomass"))

all_total_detritus.biomass.log.axis3<-output.site %>% 
  as.data.frame() %>% 
  mutate(sites = "All sites",
         detA3TCAP.rho =as.character(V1), 
         detA3TCAP.p =as.character(V5), 
         detA3TDAP.rho =as.character(V3), 
         detA3TDAP.p =as.character(V7), 
         detA3TCAP.rho =as.numeric(detA3TCAP.rho), 
         detA3TCAP.p =as.numeric(detA3TCAP.p), 
         detA3TDAP.rho =as.numeric(detA3TDAP.rho), 
         detA3TDAP.p =as.numeric(detA3TDAP.p)) %>% 
  select(sites, detA3TCAP.rho, detA3TCAP.p, detA3TDAP.rho, detA3TDAP.p)

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "logtotal_detritus_center", trait.name = "Axis.4", ro.method = "procrustes", stock = "biomass"))

all_total_detritus.biomass.log.axis4<-output.site %>% 
  as.data.frame() %>% 
  mutate(sites = "All sites",
         detA4TCAP.rho =as.character(V1), 
         detA4TCAP.p =as.character(V5), 
         detA4TDAP.rho =as.character(V3), 
         detA4TDAP.p =as.character(V7), 
         detA4TCAP.rho =as.numeric(detA4TCAP.rho), 
         detA4TCAP.p =as.numeric(detA4TCAP.p), 
         detA4TDAP.rho =as.numeric(detA4TDAP.rho), 
         detA4TDAP.p =as.numeric(detA4TDAP.p)) %>% 
  select(sites, detA4TCAP.rho, detA4TCAP.p, detA4TDAP.rho, detA4TDAP.p)

alldet1234<-all_total_detritus.biomass.log.axis1 %>% 
  full_join(all_total_detritus.biomass.log.axis2) %>% 
  full_join(all_total_detritus.biomass.log.axis3) %>% 
  full_join(all_total_detritus.biomass.log.axis4)

##  WITHIN SITES ###


###Repeated detrital analysis, using Trait Syndromes in lieu of traits###

n<-length(field.sites)

output.site<-data.frame(TCAP.rho=numeric(n),TCDAP.rho=numeric(n), TDAP.rho=numeric(n), RE.rho=numeric(n),
                        TCAP.p=numeric(n), TCDAP.p=numeric(n),TDAP.p=numeric(n), RE.p=numeric(n),
                        best.TCAP.traits=numeric(n), best.TCAPtrait.rho=numeric(n), 
                        best.TDAP.traits=numeric(n), best.TDAPtrait.rho=numeric(n),
                        best.RAO.traits=numeric(n), best.RAOtrait.rho=numeric(n),
                        cv.env=numeric(n), min.env=numeric(n), max.env=numeric(n),
                        bestTCAP.trait.rho2=numeric(n),bestTCAP.trait.p2=numeric(n),
                        bestTDAP.trait.rho2=numeric(n),bestTDAP.trait.p2=numeric(n),
                        bestRAO.trait.rho2=numeric(n),bestRAo.trait.p2=numeric(n))


for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "linear", "axis", "total_detritus", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_total_det.biomass.linear.axis<-output.site
write_csv(single_total_det.biomass.linear.axis, "single_total_det_biomass_linear_axis.csv")

output.site<-data.frame(TCAP.rho=numeric(n),TCDAP.rho=numeric(n), TDAP.rho=numeric(n), RE.rho=numeric(n),
                        TCAP.p=numeric(n), TCDAP.p=numeric(n),TDAP.p=numeric(n), RE.p=numeric(n),
                        best.TCAP.traits=numeric(n), best.TCAPtrait.rho=numeric(n), 
                        best.TDAP.traits=numeric(n), best.TDAPtrait.rho=numeric(n),
                        best.RAO.traits=numeric(n), best.RAOtrait.rho=numeric(n),
                        cv.env=numeric(n), min.env=numeric(n), max.env=numeric(n),
                        bestTCAP.trait.rho2=numeric(n),bestTCAP.trait.p2=numeric(n),
                        bestTDAP.trait.rho2=numeric(n),bestTDAP.trait.p2=numeric(n),
                        bestRAO.trait.rho2=numeric(n),bestRAo.trait.p2=numeric(n))

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "log", "axis", "total_detritus", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_total_det.biomass.log.axis<-output.site
write_csv(single_total_det.biomass.log.axis, "single_total_det_biomass_log_axis.csv")


#within site single trait axes
n<-length(field.sites)
output.site<-data.frame(TCAP.rho=numeric(n),TCDAP.rho=numeric(n), TDAP.rho=numeric(n), RE.rho=numeric(n),
                        TCAP.p=numeric(n), TCDAP.p=numeric(n),TDAP.p=numeric(n), RE.p=numeric(n))

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "log", "axis", "total_detritus", trait.name = "Axis.1", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_total_det.biomass.log.axis1<-output.site %>% 
  rownames_to_column(var="sites") %>% 
  mutate(detA1TCAP.rho = as.numeric(TCAP.rho),
         detA1TCAP.p = as.numeric(TCAP.p),
         detA1TDAP.rho = as.numeric(TDAP.rho),
         detA1TDAP.p = as.numeric(TDAP.p)) %>% 
  select(sites, detA1TCAP.rho, detA1TCAP.p, detA1TDAP.rho, detA1TDAP.p)

write_csv(single_total_det.biomass.log.axis1, "single_total_det_biomass_log_axis1.csv")

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "log", "axis", "total_detritus", trait.name = "Axis.2", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_total_det.biomass.log.axis2<-output.site %>% 
  rownames_to_column(var="sites") %>% 
  mutate(detA2TCAP.rho = as.numeric(TCAP.rho),
         detA2TCAP.p = as.numeric(TCAP.p),
         detA2TDAP.rho = as.numeric(TDAP.rho),
         detA2TDAP.p = as.numeric(TDAP.p)) %>% 
  select(sites, detA2TCAP.rho, detA2TCAP.p, detA2TDAP.rho, detA2TDAP.p)

write_csv(single_total_det.biomass.log.axis2, "single_total_det_biomass_log_axis2.csv")

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "log", "axis", "total_detritus", trait.name = "Axis.3", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_total_det.biomass.log.axis3<-output.site %>% 
  rownames_to_column(var="sites") %>% 
  mutate(detA3TCAP.rho = as.numeric(TCAP.rho),
         detA3TCAP.p = as.numeric(TCAP.p),
         detA3TDAP.rho = as.numeric(TDAP.rho),
         detA3TDAP.p = as.numeric(TDAP.p)) %>% 
  select(sites, detA3TCAP.rho, detA3TCAP.p, detA3TDAP.rho, detA3TDAP.p)

write_csv(single_total_det.biomass.log.axis3, "single_total_det_biomass_log_axis3.csv")

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "log", "axis", "total_detritus", trait.name = "Axis.4", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_total_det.biomass.log.axis4<-output.site %>% 
  rownames_to_column(var="sites") %>% 
  mutate(detA4TCAP.rho = as.numeric(TCAP.rho),
         detA4TCAP.p = as.numeric(TCAP.p),
         detA4TDAP.rho = as.numeric(TDAP.rho),
         detA4TDAP.p = as.numeric(TDAP.p)) %>% 
  select(sites, detA4TCAP.rho, detA4TCAP.p, detA4TDAP.rho, detA4TDAP.p)

write_csv(single_total_det.biomass.log.axis4, "single_total_det_biomass_log_axis4.csv")

single_total_det.biomass.log.axis1234<-single_total_det.biomass.log.axis1 %>% 
  left_join(single_total_det.biomass.log.axis2) %>% 
  left_join(single_total_det.biomass.log.axis3) %>% 
  left_join(single_total_det.biomass.log.axis4) %>% 
  full_join(alldet1234)


###############################################
#                                             #
#       Open=1, Canopy = 0                    #
#                                             #
###############################################

Broms.all <-Broms.all.full %>% 
  filter(include.canopy==1)

#Then we make a vector of all the regions
field.sites<-Broms.all$sites%>%unique()%>%as.vector()

## ACROSS SITES

n<-1
output.site<-data.frame(TCAP.rho=numeric(n),TCDAP.rho=numeric(n), TDAP.rho=numeric(n), RE.rho=numeric(n),
                        TCAP.p=numeric(n), TCDAP.p=numeric(n),TDAP.p=numeric(n), RE.p=numeric(n),
                        best.TCAP.traits=numeric(n), best.TCAPtrait.rho=numeric(n), 
                        best.TDAP.traits=numeric(n), best.TDAPtrait.rho=numeric(n),
                        best.RAO.traits=numeric(n), best.RAOtrait.rho=numeric(n),
                        cv.env=numeric(n), min.env=numeric(n), max.env=numeric(n),
                        bestTCAP.trait.rho2=numeric(n),bestTCAP.trait.p2=numeric(n),
                        bestTDAP.trait.rho2=numeric(n),bestTDAP.trait.p2=numeric(n),
                        bestRAO.trait.rho2=numeric(n),bestRAo.trait.p2=numeric(n))

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "open.canopy_center",ro.method = "procrustes", stock = "biomass"))
canopy.biomass.linear.axis<-output.site %>% unlist() %>% as.data.frame()
write_csv(canopy.biomass.linear.axis, "canopy_biomass_linear.axis.csv")

##all sites but axes singly

n<-1
output.site<-data.frame(TCAP.rho=numeric(n),TCDAP.rho=numeric(n), TDAP.rho=numeric(n), RE.rho=numeric(n),
                        TCAP.p=numeric(n), TCDAP.p=numeric(n),TDAP.p=numeric(n), RE.p=numeric(n))

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "open.canopy", trait.name = "Axis.1", ro.method = "procrustes", stock = "biomass"))

all_canopy.biomass.linear.axis1<-output.site %>% 
  as.data.frame() %>% 
  mutate(sites = "All sites",
         canA1TCAP.rho =as.character(V1), 
         canA1TCAP.p =as.character(V5), 
         canA1TDAP.rho =as.character(V3), 
         canA1TDAP.p =as.character(V7), 
         canA1TCAP.rho =as.numeric(canA1TCAP.rho), 
         canA1TCAP.p =as.numeric(canA1TCAP.p), 
         canA1TDAP.rho =as.numeric(canA1TDAP.rho), 
         canA1TDAP.p =as.numeric(canA1TDAP.p)) %>% 
  select(sites, canA1TCAP.rho, canA1TCAP.p, canA1TDAP.rho, canA1TDAP.p)

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "open.canopy", trait.name = "Axis.2", ro.method = "procrustes", stock = "biomass"))

all_canopy.biomass.linear.axis2<-output.site %>% 
  as.data.frame() %>% 
  mutate(sites = "All sites",
         canA2TCAP.rho =as.character(V1), 
         canA2TCAP.p =as.character(V5), 
         canA2TDAP.rho =as.character(V3), 
         canA2TDAP.p =as.character(V7), 
         canA2TCAP.rho =as.numeric(canA2TCAP.rho), 
         canA2TCAP.p =as.numeric(canA2TCAP.p), 
         canA2TDAP.rho =as.numeric(canA2TDAP.rho), 
         canA2TDAP.p =as.numeric(canA2TDAP.p)) %>% 
  select(sites, canA2TCAP.rho, canA2TCAP.p, canA2TDAP.rho, canA2TDAP.p)

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "open.canopy", trait.name = "Axis.3", ro.method = "procrustes", stock = "biomass"))

all_canopy.biomass.linear.axis3<-output.site %>% 
  as.data.frame() %>% 
  mutate(sites = "All sites",
         canA3TCAP.rho =as.character(V1), 
         canA3TCAP.p =as.character(V5), 
         canA3TDAP.rho =as.character(V3), 
         canA3TDAP.p =as.character(V7), 
         canA3TCAP.rho =as.numeric(canA3TCAP.rho), 
         canA3TCAP.p =as.numeric(canA3TCAP.p), 
         canA3TDAP.rho =as.numeric(canA3TDAP.rho), 
         canA3TDAP.p =as.numeric(canA3TDAP.p)) %>% 
  select(sites, canA3TCAP.rho, canA3TCAP.p, canA3TDAP.rho, canA3TDAP.p)

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "open.canopy", trait.name = "Axis.4", ro.method = "procrustes", stock = "biomass"))

all_canopy.biomass.linear.axis4<-output.site %>% 
  as.data.frame() %>% 
  mutate(sites = "All sites",
         canA4TCAP.rho =as.character(V1), 
         canA4TCAP.p =as.character(V5), 
         canA4TDAP.rho =as.character(V3), 
         canA4TDAP.p =as.character(V7), 
         canA4TCAP.rho =as.numeric(canA4TCAP.rho), 
         canA4TCAP.p =as.numeric(canA4TCAP.p), 
         canA4TDAP.rho =as.numeric(canA4TDAP.rho), 
         canA4TDAP.p =as.numeric(canA4TDAP.p)) %>% 
  select(sites, canA4TCAP.rho, canA4TCAP.p, canA4TDAP.rho, canA4TDAP.p)

allcan1234<-all_canopy.biomass.linear.axis1 %>% 
  full_join(all_canopy.biomass.linear.axis2) %>% 
  full_join(all_canopy.biomass.linear.axis3) %>% 
  full_join(all_canopy.biomass.linear.axis4)

## WITHIN SITES

n<-length(field.sites)
output.site<-data.frame(TCAP.rho=numeric(n),TCDAP.rho=numeric(n), TDAP.rho=numeric(n), RE.rho=numeric(n),
                        TCAP.p=numeric(n), TCDAP.p=numeric(n),TDAP.p=numeric(n), RE.p=numeric(n),
                        best.TCAP.traits=numeric(n), best.TCAPtrait.rho=numeric(n), 
                        best.TDAP.traits=numeric(n), best.TDAPtrait.rho=numeric(n),
                        best.RAO.traits=numeric(n), best.RAOtrait.rho=numeric(n),
                        cv.env=numeric(n), min.env=numeric(n), max.env=numeric(n),
                        bestTCAP.trait.rho2=numeric(n),bestTCAP.trait.p2=numeric(n),
                        bestTDAP.trait.rho2=numeric(n),bestTDAP.trait.p2=numeric(n),
                        bestRAO.trait.rho2=numeric(n),bestRAo.trait.p2=numeric(n))

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "linear", "axis", "open.canopy", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_canopy.biomass.linear.axis<-output.site
write_csv(single_canopy.biomass.linear.axis, "single_canopy_biomass_linear_axis.csv")

#within site single trait axes
n<-length(field.sites)
output.site<-data.frame(TCAP.rho=numeric(n),TCDAP.rho=numeric(n), TDAP.rho=numeric(n), RE.rho=numeric(n),
                        TCAP.p=numeric(n), TCDAP.p=numeric(n),TDAP.p=numeric(n), RE.p=numeric(n))

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "linear", "axis", "open.canopy", trait.name = "Axis.1", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_canopy.biomass.linear.axis1<-output.site %>% 
  rownames_to_column(var="sites") %>% 
  mutate(canA1TCAP.rho = as.numeric(TCAP.rho),
         canA1TCAP.p = as.numeric(TCAP.p),
         canA1TDAP.rho = as.numeric(TDAP.rho),
         canA1TDAP.p = as.numeric(TDAP.p)) %>% 
  select(sites, canA1TCAP.rho, canA1TCAP.p, canA1TDAP.rho, canA1TDAP.p)

write_csv(single_canopy.biomass.linear.axis1, "single_canopy_biomass_linear_axis1.csv")

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "linear", "axis", "open.canopy", trait.name = "Axis.2", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_canopy.biomass.linear.axis2<-output.site %>% 
  rownames_to_column(var="sites") %>% 
  mutate(canA2TCAP.rho = as.numeric(TCAP.rho),
         canA2TCAP.p = as.numeric(TCAP.p),
         canA2TDAP.rho = as.numeric(TDAP.rho),
         canA2TDAP.p = as.numeric(TDAP.p)) %>% 
  select(sites, canA2TCAP.rho, canA2TCAP.p, canA2TDAP.rho, canA2TDAP.p)

write_csv(single_canopy.biomass.linear.axis2, "single_canopy_biomass_linear_axis2.csv")

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "linear", "axis", "open.canopy", trait.name = "Axis.3", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_canopy.biomass.linear.axis3<-output.site %>% 
  rownames_to_column(var="sites") %>% 
  mutate(canA3TCAP.rho = as.numeric(TCAP.rho),
         canA3TCAP.p = as.numeric(TCAP.p),
         canA3TDAP.rho = as.numeric(TDAP.rho),
         canA3TDAP.p = as.numeric(TDAP.p)) %>% 
  select(sites, canA3TCAP.rho, canA3TCAP.p, canA3TDAP.rho, canA3TDAP.p)

write_csv(single_canopy.biomass.linear.axis3, "single_canopy_biomass_linear_axis3.csv")

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "linear", "axis", "open.canopy", trait.name = "Axis.4", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_canopy.biomass.linear.axis4<-output.site %>% 
  rownames_to_column(var="sites") %>% 
  mutate(canA4TCAP.rho = as.numeric(TCAP.rho),
         canA4TCAP.p = as.numeric(TCAP.p),
         canA4TDAP.rho = as.numeric(TDAP.rho),
         canA4TDAP.p = as.numeric(TDAP.p)) %>% 
  select(sites, canA4TCAP.rho, canA4TCAP.p, canA4TDAP.rho, canA4TDAP.p)

write_csv(single_canopy.biomass.linear.axis4, "single_canopy_biomass_linear_axis4.csv")

single_canopy.biomass.linear.axis1234<-single_canopy.biomass.linear.axis1 %>% 
  left_join(single_canopy.biomass.linear.axis2) %>% 
  left_join(single_canopy.biomass.linear.axis3) %>% 
  left_join(single_canopy.biomass.linear.axis4) %>% 
  full_join(allcan1234)




##############################################
#                                             #
#       actual_water                          #
#                                             #
###############################################

Broms.all <-Broms.all.full %>% 
  filter(include.water==1)

#Then we make a vector of all the regions
field.sites<-Broms.all$sites%>%unique()%>%as.vector()

###remove this:
#single.site.syncsa("Pitilla", "linear", "axis", "actual_water", ro.method = "procrustes", stock = "biomass")


##  ACROSS SITES


n<-1
output.site<-data.frame(TCAP.rho=numeric(n),TCDAP.rho=numeric(n), TDAP.rho=numeric(n), RE.rho=numeric(n),
                        TCAP.p=numeric(n), TCDAP.p=numeric(n),TDAP.p=numeric(n), RE.p=numeric(n),
                        best.TCAP.traits=numeric(n), best.TCAPtrait.rho=numeric(n), 
                        best.TDAP.traits=numeric(n), best.TDAPtrait.rho=numeric(n),
                        best.RAO.traits=numeric(n), best.RAOtrait.rho=numeric(n),
                        cv.env=numeric(n), min.env=numeric(n), max.env=numeric(n),
                        bestTCAP.trait.rho2=numeric(n),bestTCAP.trait.p2=numeric(n),
                        bestTDAP.trait.rho2=numeric(n),bestTDAP.trait.p2=numeric(n),
                        bestRAO.trait.rho2=numeric(n),bestRAo.trait.p2=numeric(n))

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "actual_water_center", ro.method = "procrustes", stock = "biomass"))
actualw.biomass.linear.axis<-output.site %>% unlist() %>% as.data.frame()
write_csv(actualw.biomass.linear.axis, "actualw_biomass_linear_axis.csv")

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "logactual_water_center", ro.method = "procrustes",stock = "biomass"))
actualw.biomass.log.axis<-output.site %>% unlist() %>% as.data.frame()
write_csv(actualw.biomass.log.axis, "actualw_biomass_log_axis.csv")

##all sites but axes singly

n<-1
output.site<-data.frame(TCAP.rho=numeric(n),TCDAP.rho=numeric(n), TDAP.rho=numeric(n), RE.rho=numeric(n), TCAP.p=numeric(n), TCDAP.p=numeric(n),TDAP.p=numeric(n), RE.p=numeric(n))

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "logactual_water_center", trait.name = "Axis.1", ro.method = "procrustes", stock = "biomass"))

all_actual_water.biomass.log.axis1<-output.site %>% 
  as.data.frame() %>% 
  mutate(sites = "All sites",
         watA1TCAP.rho =as.character(V1), 
         watA1TCAP.p =as.character(V5), 
         watA1TDAP.rho =as.character(V3), 
         watA1TDAP.p =as.character(V7), 
         watA1TCAP.rho =as.numeric(watA1TCAP.rho), 
         watA1TCAP.p =as.numeric(watA1TCAP.p), 
         watA1TDAP.rho =as.numeric(watA1TDAP.rho), 
         watA1TDAP.p =as.numeric(watA1TDAP.p)) %>% 
  select(sites, watA1TCAP.rho, watA1TCAP.p, watA1TDAP.rho, watA1TDAP.p)

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "logactual_water_center", trait.name = "Axis.2", ro.method = "procrustes", stock = "biomass"))

all_actual_water.biomass.log.axis2<-output.site %>% 
  as.data.frame() %>% 
  mutate(sites = "All sites",
         watA2TCAP.rho =as.character(V1), 
         watA2TCAP.p =as.character(V5), 
         watA2TDAP.rho =as.character(V3), 
         watA2TDAP.p =as.character(V7), 
         watA2TCAP.rho =as.numeric(watA2TCAP.rho), 
         watA2TCAP.p =as.numeric(watA2TCAP.p), 
         watA2TDAP.rho =as.numeric(watA2TDAP.rho), 
         watA2TDAP.p =as.numeric(watA2TDAP.p)) %>% 
  select(sites, watA2TCAP.rho, watA2TCAP.p, watA2TDAP.rho, watA2TDAP.p)

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "logactual_water_center", trait.name = "Axis.3", ro.method = "procrustes", stock = "biomass"))

all_actual_water.biomass.log.axis3<-output.site %>% 
  as.data.frame() %>% 
  mutate(sites = "All sites",
         watA3TCAP.rho =as.character(V1), 
         watA3TCAP.p =as.character(V5), 
         watA3TDAP.rho =as.character(V3), 
         watA3TDAP.p =as.character(V7), 
         watA3TCAP.rho =as.numeric(watA3TCAP.rho), 
         watA3TCAP.p =as.numeric(watA3TCAP.p), 
         watA3TDAP.rho =as.numeric(watA3TDAP.rho), 
         watA3TDAP.p =as.numeric(watA3TDAP.p)) %>% 
  select(sites, watA3TCAP.rho, watA3TCAP.p, watA3TDAP.rho, watA3TDAP.p)

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "logactual_water_center", trait.name = "Axis.4", ro.method = "procrustes", stock = "biomass"))

all_actual_water.biomass.log.axis4<-output.site %>% 
  as.data.frame() %>% 
  mutate(sites = "All sites",
         watA4TCAP.rho =as.character(V1), 
         watA4TCAP.p =as.character(V5), 
         watA4TDAP.rho =as.character(V3), 
         watA4TDAP.p =as.character(V7), 
         watA4TCAP.rho =as.numeric(watA4TCAP.rho), 
         watA4TCAP.p =as.numeric(watA4TCAP.p), 
         watA4TDAP.rho =as.numeric(watA4TDAP.rho), 
         watA4TDAP.p =as.numeric(watA4TDAP.p)) %>% 
  select(sites, watA4TCAP.rho, watA4TCAP.p, watA4TDAP.rho, watA4TDAP.p)

allwat1234<-all_actual_water.biomass.log.axis1 %>% 
  full_join(all_actual_water.biomass.log.axis2) %>% 
  full_join(all_actual_water.biomass.log.axis3) %>% 
  full_join(all_actual_water.biomass.log.axis4)

## WITHIN SITES

n<-length(field.sites)
output.site<-data.frame(TCAP.rho=numeric(n),TCDAP.rho=numeric(n), TDAP.rho=numeric(n), RE.rho=numeric(n),
                        TCAP.p=numeric(n), TCDAP.p=numeric(n),TDAP.p=numeric(n), RE.p=numeric(n),
                        best.TCAP.traits=numeric(n), best.TCAPtrait.rho=numeric(n), 
                        best.TDAP.traits=numeric(n), best.TDAPtrait.rho=numeric(n),
                        best.RAO.traits=numeric(n), best.RAOtrait.rho=numeric(n),
                        cv.env=numeric(n), min.env=numeric(n), max.env=numeric(n),
                        bestTCAP.trait.rho2=numeric(n),bestTCAP.trait.p2=numeric(n),
                        bestTDAP.trait.rho2=numeric(n),bestTDAP.trait.p2=numeric(n),
                        bestRAO.trait.rho2=numeric(n),bestRAo.trait.p2=numeric(n))


for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "linear", "axis", "actual_water", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}
single_actualw.biomass.linear.axis<-output.site
write_csv(single_actualw.biomass.linear.axis, "single_actualw_biomass_linear_axis.csv")


for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "log", "axis", "actual_water", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_actualw.biomass.log.axis<-output.site
write_csv(single_actualw.biomass.log.axis, "single_actualw_biomass_log_axis.csv")


#within site single trait axes
n<-length(field.sites)
output.site<-data.frame(TCAP.rho=numeric(n),TCDAP.rho=numeric(n), TDAP.rho=numeric(n), RE.rho=numeric(n),
                        TCAP.p=numeric(n), TCDAP.p=numeric(n),TDAP.p=numeric(n), RE.p=numeric(n))

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "log", "axis", "actual_water", trait.name = "Axis.1", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_total_wat.biomass.log.axis1<-output.site %>% 
  rownames_to_column(var="sites") %>% 
  mutate(watA1TCAP.rho = as.numeric(TCAP.rho),
         watA1TCAP.p = as.numeric(TCAP.p),
         watA1TDAP.rho = as.numeric(TDAP.rho),
         watA1TDAP.p = as.numeric(TDAP.p)) %>% 
  select(sites, watA1TCAP.rho, watA1TCAP.p, watA1TDAP.rho, watA1TDAP.p)

write_csv(single_total_wat.biomass.log.axis1, "single_actual_water_biomass_log_axis1.csv")

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "log", "axis", "actual_water", trait.name = "Axis.2", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_total_wat.biomass.log.axis2<-output.site %>% 
  rownames_to_column(var="sites") %>% 
  mutate(watA2TCAP.rho = as.numeric(TCAP.rho),
         watA2TCAP.p = as.numeric(TCAP.p),
         watA2TDAP.rho = as.numeric(TDAP.rho),
         watA2TDAP.p = as.numeric(TDAP.p)) %>% 
  select(sites, watA2TCAP.rho, watA2TCAP.p, watA2TDAP.rho, watA2TDAP.p)

write_csv(single_total_wat.biomass.log.axis2, "single_actual_water_biomass_log_axis2.csv")

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "log", "axis", "actual_water", trait.name = "Axis.3", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_total_wat.biomass.log.axis3<-output.site %>% 
  rownames_to_column(var="sites") %>% 
  mutate(watA3TCAP.rho = as.numeric(TCAP.rho),
         watA3TCAP.p = as.numeric(TCAP.p),
         watA3TDAP.rho = as.numeric(TDAP.rho),
         watA3TDAP.p = as.numeric(TDAP.p)) %>% 
  select(sites, watA3TCAP.rho, watA3TCAP.p, watA3TDAP.rho, watA3TDAP.p)

write_csv(single_total_wat.biomass.log.axis3, "single_actual_water_biomass_log_axis3.csv")

for (i in 1:n)
{
  output.site[i,]<-single.site.syncsa(field.sites[i], "log", "axis", "actual_water", trait.name = "Axis.4", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_total_wat.biomass.log.axis4<-output.site %>% 
  rownames_to_column(var="sites") %>% 
  mutate(watA4TCAP.rho = as.numeric(TCAP.rho),
         watA4TCAP.p = as.numeric(TCAP.p),
         watA4TDAP.rho = as.numeric(TDAP.rho),
         watA4TDAP.p = as.numeric(TDAP.p)) %>% 
  select(sites, watA4TCAP.rho, watA4TCAP.p, watA4TDAP.rho, watA4TDAP.p)

write_csv(single_total_wat.biomass.log.axis4, "single_actual_water_biomass_log_axis4.csv")

single_total_wat.biomass.log.axis1234<-single_total_wat.biomass.log.axis1 %>% 
  left_join(single_total_wat.biomass.log.axis2) %>% 
  left_join(single_total_wat.biomass.log.axis3) %>% 
  left_join(single_total_wat.biomass.log.axis4) %>% 
  full_join(allwat1234)





###############################################
#                                             #
#       save as RDS                           #
#                                             #
###############################################



syncsa_summary_list<- list(
  total_det.biomass.linear.axis = total_det.biomass.linear.axis,
  total_det.biomass.log.axis  = total_det.biomass.log.axis,
  single_total_det.biomass.linear.axis = single_total_det.biomass.linear.axis,
  single_total_det.biomass.log.axis = single_total_det.biomass.log.axis,
  canopy.biomass.linear.axis = canopy.biomass.linear.axis,
  single_canopy.biomass.linear.axis = single_canopy.biomass.linear.axis,
  actualw.biomass.linear.axis = actualw.biomass.linear.axis,
  actualw.biomass.log.axis = actualw.biomass.log.axis,
  single_actualw.biomass.linear.axis = single_actualw.biomass.linear.axis,
  single_actualw.biomass.log.axis = single_actualw.biomass.log.axis
  )



#----allout---

detout <- syncsa_summary_list$single_total_det.biomass.log.axis %>% 
  rownames_to_column(var = "sites") %>% 
  select(sites, TCAP.rho, TCAP.p, TDAP.rho, TDAP.p, best.TCAP.traits, best.TDAP.traits) %>% 
  rbind(c("line",0,0,0,0,NA,NA)) %>% 
  mutate(detbest.TCAP.traits =as.character(best.TCAP.traits), 
         detbest.TDAP.traits =as.character(best.TDAP.traits),
         detTCAP.rho =as.numeric(TCAP.rho), 
         detTCAP.p =as.numeric(TCAP.p), 
         detTDAP.rho =as.numeric(TDAP.rho), 
         detTDAP.p =as.numeric(TDAP.p)) %>% 
  select(sites, detTCAP.rho, detTCAP.p, detTDAP.rho, detTDAP.p, detbest.TCAP.traits, detbest.TDAP.traits)

watout <- syncsa_summary_list$single_actualw.biomass.log.axis %>% 
  rownames_to_column(var = "sites") %>% 
  select(sites, TCAP.rho, TCAP.p, TDAP.rho, TDAP.p, best.TCAP.traits, best.TDAP.traits) %>% 
  rbind(c("line",0,0,0,0,NA,NA)) %>% 
  mutate(watbest.TCAP.traits =as.character(best.TCAP.traits), 
         watbest.TDAP.traits =as.character(best.TDAP.traits),
         watTCAP.rho =as.numeric(TCAP.rho), 
         watTCAP.p =as.numeric(TCAP.p), 
         watTDAP.rho =as.numeric(TDAP.rho), 
         watTDAP.p =as.numeric(TDAP.p)) %>% 
  select(sites, watTCAP.rho, watTCAP.p, watTDAP.rho, watTDAP.p, watbest.TCAP.traits, watbest.TDAP.traits)


canout <- syncsa_summary_list$single_canopy.biomass.linear.axis %>% 
  rownames_to_column(var = "sites") %>% 
  select(sites, TCAP.rho, TCAP.p, TDAP.rho, TDAP.p, best.TCAP.traits, best.TDAP.traits) %>% 
  rbind(c("line",0,0.051,0,0.051,NA,NA)) %>% #a hack to create a phantom level of significance
  mutate(canbest.TCAP.traits =as.character(best.TCAP.traits), 
         canbest.TDAP.traits =as.character(best.TDAP.traits),
         canTCAP.rho =as.numeric(TCAP.rho), 
         canTCAP.p =as.numeric(TCAP.p), 
         canTDAP.rho =as.numeric(TDAP.rho), 
         canTDAP.p =as.numeric(TDAP.p)) %>% 
  select(sites, canTCAP.rho, canTCAP.p, canTDAP.rho, canTDAP.p, canbest.TCAP.traits, canbest.TDAP.traits)

alldet<-syncsa_summary_list$total_det.biomass.log.axis %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(detTCAP.rho =as.character(V1), 
         detTCAP.p =as.character(V5), 
         detTDAP.rho =as.character(V3), 
         detTDAP.p =as.character(V7), 
         detbest.TCAP.traits =as.character(V9), 
         detbest.TDAP.traits =as.character(V11),
         detTCAP.rho =as.numeric(detTCAP.rho), 
         detTCAP.p =as.numeric(detTCAP.p), 
         detTDAP.rho =as.numeric(detTDAP.rho), 
         detTDAP.p =as.numeric(detTDAP.p)
  ) %>%  
  mutate(sites = "All sites") %>% 
  select(sites, detTCAP.rho, detTCAP.p, detTDAP.rho, detTDAP.p, detbest.TCAP.traits, detbest.TDAP.traits)


allwat<-syncsa_summary_list$actualw.biomass.log.axis %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(watTCAP.rho =as.character(V1), 
         watTCAP.p =as.character(V5), 
         watTDAP.rho =as.character(V3), 
         watTDAP.p =as.character(V7), 
         watbest.TCAP.traits =as.character(V9), 
         watbest.TDAP.traits =as.character(V11),
         watTCAP.rho =as.numeric(watTCAP.rho), 
         watTCAP.p =as.numeric(watTCAP.p), 
         watTDAP.rho =as.numeric(watTDAP.rho), 
         watTDAP.p =as.numeric(watTDAP.p)
  ) %>% 
  mutate(sites = "All sites") %>% 
  select(sites, watTCAP.rho, watTCAP.p, watTDAP.rho, watTDAP.p, watbest.TCAP.traits, watbest.TDAP.traits)

allcan<-syncsa_summary_list$canopy.biomass.linear.axis %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(canTCAP.rho =as.character(V1), 
         canTCAP.p =as.character(V5), 
         canTDAP.rho =as.character(V3), 
         canTDAP.p =as.character(V7), 
         canbest.TCAP.traits =as.character(V9), 
         canbest.TDAP.traits =as.character(V11),
         canTCAP.rho =as.numeric(canTCAP.rho), 
         canTCAP.p =as.numeric(canTCAP.p), 
         canTDAP.rho =as.numeric(canTDAP.rho), 
         canTDAP.p =as.numeric(canTDAP.p)
  ) %>% 
  mutate(sites = "All sites") %>% 
  select(sites, canTCAP.rho, canTCAP.p, canTDAP.rho, canTDAP.p, canbest.TCAP.traits, canbest.TDAP.traits)

wat<-full_join(watout,allwat) 
det<-full_join(detout,alldet) 
can<-full_join(canout,allcan)

allout<-wat %>% 
  full_join(det, by = c("sites"="sites")) %>% 
  left_join(can) %>% 
  mutate(detTCAP.sig = ifelse(detTCAP.p<0.0505,"sig", 
                              (ifelse(detTCAP.p<0.0705,"marg", "ns")))) %>% 
  mutate(watTCAP.sig = ifelse(watTCAP.p<0.0505,"sig", 
                              (ifelse(watTCAP.p<0.0705,"marg", "ns")))) %>% 
 mutate(canTCAP.sig = ifelse(canTCAP.p<0.0505,"sig",
                              (ifelse(canTCAP.p<0.0705,"marg", "ns")))) %>%
  mutate(detTDAP.sig = ifelse(detTDAP.p<0.0505,"sig", 
                              (ifelse(detTDAP.p<0.0705,"marg", "ns")))) %>% 
  mutate(watTDAP.sig = ifelse(watTDAP.p<0.0505,"sig", 
                              (ifelse(watTDAP.p<0.0705,"marg", "ns")))) %>% 
  mutate(canTDAP.sig = ifelse(canTDAP.p<0.0505,"sig", 
                              (ifelse(canTDAP.p<0.0705,"marg", "ns")))) 

#-------allaxes1234----------------

allaxes1234 <- single_total_wat.biomass.log.axis1234 %>% 
  full_join(single_total_det.biomass.log.axis1234, by = c("sites"="sites")) %>% 
  left_join(single_canopy.biomass.linear.axis1234) %>% 
  rbind(c("line", rep(0,64)))

##-----write files

saveRDS(syncsa_summary_list, "syncsa_summary_sept2021.rds")
write_csv(allout, "allout.csv")
write_csv(allaxes1234, "allaxes1234.csv")


