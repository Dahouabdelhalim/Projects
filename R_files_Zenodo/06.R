library(dplyr)
library(purrr)
library(tidyverse)
library(stringr)
library(SYNCSA)

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
  filter(include.detritus==1) %>% 
  filter(sites %nin% c("Cusuco", "DOM_LOW", "EVDW", "EVPC", "EVTAB", "SONLOW", "Pitilla", "SABA_MID", "Quintana_Roo", "SONHI"))

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

output.site<-as.data.frame(multi.site.syncsa("linear", "axis", "logtotal_detritus_center", ro.method = "procrustes", stock = "biomass"))
total_det.biomass.log.axis<-output.site %>% unlist() %>% as.data.frame()
write_csv(total_det.biomass.log.axis, "total_det_biomass_log_axis_south.csv")

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
  output.site[i,]<-single.site.syncsa(field.sites[i], "log", "axis", "total_detritus", ro.method = "procrustes", stock = "biomass")
  rownames(output.site)[i]<-field.sites[i]
}

single_total_det.biomass.log.axis<-output.site
write_csv(single_total_det.biomass.log.axis, "single_total_det_biomass_log_axis_south.csv")


###############################################
#                                             #
#       save                                  #
#                                             #
###############################################


#----allout---

detout <- single_total_det.biomass.log.axis %>% 
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

alldet<-total_det.biomass.log.axis %>% 
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

det<-full_join(detout,alldet) %>% 
  mutate(detTCAP.sig = ifelse(detTCAP.p<0.0505,"sig", 
                              (ifelse(detTCAP.p<0.10,"marg", "ns"))))%>%
  mutate(detTDAP.sig = ifelse(detTDAP.p<0.0505,"sig", 
                              (ifelse(detTCAP.p<0.10,"marg", "ns")))) 

#-------allaxes1234----------------

allaxes1234_det_south <- alldet1234 %>% 
  rbind(c("line", rep(0,64)))

##-----write files

write_csv(det, "all_det_south.csv")
write_csv(allaxes1234_det_south, "allaxes1234_det_south.csv")


