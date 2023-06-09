# read data from csv, filter, and organize

library(tidyverse)
Clarkia_Field_2010 <- read.csv("ClarkiaFieldData2010.csv")
Clarkia_Field_2010$Taxon <- as.factor(Clarkia_Field_2010$Taxon)
Clarkia_Field_2010$Population <- as.factor(Clarkia_Field_2010$Population)
Clarkia_Field_2010$Reproductive_status <- as.factor(Clarkia_Field_2010$Reproductive_status)
Clarkia_Field_filtered <- Clarkia_Field_2010 %>% filter(Photo_DayGE > 2 & Photo_DayGE < 50)
Clarkia_Field_filtered <- Clarkia_Field_filtered %>% filter(Tleaf_DayGE > 0 & Tleaf_DayGE < 43)
Clarkia_Field_filtered <- Clarkia_Field_filtered %>% filter(Cond_DayGE > 0 & Cond_DayGE < 2)
Clarkia_Field_filtered <- Clarkia_Field_filtered %>% filter(Trans_DayGE > 0 & Trans_DayGE < 0.04)
Clarkia_Field_filtered <- Clarkia_Field_filtered %>% filter(WUE_DayGE > 0 & WUE_DayGE < 1.5)
unguiculata_Rep <- Clarkia_Field_filtered %>% filter(Taxon == "unguiculata", Reproductive_status == "Reproductive")
unguiculata_Veg <- Clarkia_Field_filtered %>% filter(Taxon == "unguiculata", Reproductive_status == "Vegetative")
exilis_Rep <- Clarkia_Field_filtered %>% filter(Taxon == "exilis", Reproductive_status == "Reproductive")
exilis_Veg <- Clarkia_Field_filtered %>% filter(Taxon == "exilis", Reproductive_status == "Vegetative")
xantiana_Rep <- Clarkia_Field_filtered %>% filter(Taxon == "xantiana", Reproductive_status == "Reproductive")
xantiana_Veg <- Clarkia_Field_filtered %>% filter(Taxon == "xantiana", Reproductive_status == "Vegetative")
parviflora_Rep <- Clarkia_Field_filtered %>% filter(Taxon == "parviflora", Reproductive_status == "Reproductive")
parviflora_Veg <- Clarkia_Field_filtered %>% filter(Taxon == "parviflora", Reproductive_status == "Vegetative")

Clarkia_Field_filtered$Reproductive_status <- factor(Clarkia_Field_filtered$Reproductive_status, levels = c("Vegetative", "Reproductive"))
Clarkia_Field_filtered$Taxon <- factor(Clarkia_Field_filtered$Taxon, levels = c("exilis", "unguiculata", "parviflora", "xantiana"))
CeCu <- Clarkia_Field_filtered %>% filter(Taxon == "unguiculata" | Taxon == "exilis")
CeCu$Population <- as.character(CeCu$Population)
CeCu$Population[CeCu$Taxon == "unguiculata" & CeCu$Population == "StarkCreek"] <- "StarkCreekCU"
CeCu$Population <- as.character(CeCu$Population)
CeCu$Population[CeCu$Taxon == "exilis" & CeCu$Population == "StarkCreek"] <- "StarkCreekCE"
CeCu$Population <- factor(CeCu$Population)
CxxCxp <- Clarkia_Field_filtered %>% filter(Taxon == "xantiana" | Taxon == "parviflora")
CxxCxp$Population <- as.character(CxxCxp$Population)
CxxCxp$Population[CxxCxp$Taxon == "xantiana" & CxxCxp$Population == "SawmillRoad"] <- "SawmillRoadCXX"
CxxCxp$Population <- as.character(CxxCxp$Population)
CxxCxp$Population[CxxCxp$Taxon == "parviflora" & CxxCxp$Population == "SawmillRoad"] <- "SawmillRoadCXP"
CxxCxp$Population <- factor(CxxCxp$Population)
unguiculata_Rep_MC_Std <- unguiculata_Rep %>% mutate(
        MCstd_Biomass = (Biomass - mean(Biomass, na.rm = TRUE)) / sd(Biomass, na.rm = TRUE),
        MCstd_Photo_DayGE = (Photo_DayGE - mean(Photo_DayGE, na.rm = TRUE)) / sd(Photo_DayGE, na.rm = TRUE),
        MCstd_Trans_DayGE = (Trans_DayGE - mean(Trans_DayGE, na.rm = TRUE)) / sd(Trans_DayGE, na.rm = TRUE),
        MCstd_WUE_DayGE = (WUE_DayGE - mean(WUE_DayGE, na.rm = TRUE)) / sd(WUE_DayGE, na.rm = TRUE),
        MC_Tleaf_DayGE = (Tleaf_DayGE - mean(Tleaf_DayGE, na.rm = TRUE)),
        MC_DayGE_Node = (DayGE_Node - mean(DayGE_Node, na.rm = TRUE)),
        MCstd_Photo_DayGE_squared = MCstd_Photo_DayGE^2,
        MCstd_Trans_DayGE_squared = MCstd_Trans_DayGE^2,
        MCstd_WUE_DayGE_squared = MCstd_WUE_DayGE^2
)
exilis_Rep_MC_Std <- exilis_Rep %>% mutate(
        MCstd_Biomass = (Biomass - mean(Biomass, na.rm = TRUE)) / sd(Biomass, na.rm = TRUE),
        MCstd_Photo_DayGE = (Photo_DayGE - mean(Photo_DayGE, na.rm = TRUE)) / sd(Photo_DayGE, na.rm = TRUE),
        MCstd_Trans_DayGE = (Trans_DayGE - mean(Trans_DayGE, na.rm = TRUE)) / sd(Trans_DayGE, na.rm = TRUE),
        MCstd_WUE_DayGE = (WUE_DayGE - mean(WUE_DayGE, na.rm = TRUE)) / sd(WUE_DayGE, na.rm = TRUE),
        MC_Tleaf_DayGE = (Tleaf_DayGE - mean(Tleaf_DayGE, na.rm = TRUE)),
        MC_DayGE_Node = (DayGE_Node - mean(DayGE_Node, na.rm = TRUE)),
        MCstd_Photo_DayGE_squared = MCstd_Photo_DayGE^2,
        MCstd_Trans_DayGE_squared = MCstd_Trans_DayGE^2,
        MCstd_WUE_DayGE_squared = MCstd_WUE_DayGE^2
)
xantiana_Rep_MC_Std <- xantiana_Rep %>% mutate(
        MCstd_Biomass = (Biomass - mean(Biomass, na.rm = TRUE)) / sd(Biomass, na.rm = TRUE),
        MCstd_Photo_DayGE = (Photo_DayGE - mean(Photo_DayGE, na.rm = TRUE)) / sd(Photo_DayGE, na.rm = TRUE),
        MCstd_Trans_DayGE = (Trans_DayGE - mean(Trans_DayGE, na.rm = TRUE)) / sd(Trans_DayGE, na.rm = TRUE),
        MCstd_WUE_DayGE = (WUE_DayGE - mean(WUE_DayGE, na.rm = TRUE)) / sd(WUE_DayGE, na.rm = TRUE),
        MC_Tleaf_DayGE = (Tleaf_DayGE - mean(Tleaf_DayGE, na.rm = TRUE)),
        MC_DayGE_Node = (DayGE_Node - mean(DayGE_Node, na.rm = TRUE)),
        MCstd_Photo_DayGE_squared = MCstd_Photo_DayGE^2,
        MCstd_Trans_DayGE_squared = MCstd_Trans_DayGE^2,
        MCstd_WUE_DayGE_squared = MCstd_WUE_DayGE^2
)
parviflora_Rep_MC_Std <- parviflora_Rep %>% mutate(
        MCstd_Biomass = (Biomass - mean(Biomass, na.rm = TRUE)) / sd(Biomass, na.rm = TRUE),
        MCstd_Photo_DayGE = (Photo_DayGE - mean(Photo_DayGE, na.rm = TRUE)) / sd(Photo_DayGE, na.rm = TRUE),
        MCstd_Trans_DayGE = (Trans_DayGE - mean(Trans_DayGE, na.rm = TRUE)) / sd(Trans_DayGE, na.rm = TRUE),
        MCstd_WUE_DayGE = (WUE_DayGE - mean(WUE_DayGE, na.rm = TRUE)) / sd(WUE_DayGE, na.rm = TRUE),
        MC_Tleaf_DayGE = (Tleaf_DayGE - mean(Tleaf_DayGE, na.rm = TRUE)),
        MC_DayGE_Node = (DayGE_Node - mean(DayGE_Node, na.rm = TRUE)),
        MCstd_Photo_DayGE_squared = MCstd_Photo_DayGE^2,
        MCstd_Trans_DayGE_squared = MCstd_Trans_DayGE^2,
        MCstd_WUE_DayGE_squared = MCstd_WUE_DayGE^2
)
unguiculata_Veg_MC_Std <- unguiculata_Veg %>% mutate(
        MCstd_Biomass = (Biomass - mean(Biomass, na.rm = TRUE)) / sd(Biomass, na.rm = TRUE),
        MCstd_Photo_DayGE = (Photo_DayGE - mean(Photo_DayGE, na.rm = TRUE)) / sd(Photo_DayGE, na.rm = TRUE),
        MCstd_Trans_DayGE = (Trans_DayGE - mean(Trans_DayGE, na.rm = TRUE)) / sd(Trans_DayGE, na.rm = TRUE),
        MCstd_WUE_DayGE = (WUE_DayGE - mean(WUE_DayGE, na.rm = TRUE)) / sd(WUE_DayGE, na.rm = TRUE),
        MC_Tleaf_DayGE = (Tleaf_DayGE - mean(Tleaf_DayGE, na.rm = TRUE)),
        MC_DayGE_Node = (DayGE_Node - mean(DayGE_Node, na.rm = TRUE)),
        MCstd_Photo_DayGE_squared = MCstd_Photo_DayGE^2,
        MCstd_Trans_DayGE_squared = MCstd_Trans_DayGE^2,
        MCstd_WUE_DayGE_squared = MCstd_WUE_DayGE^2
)
exilis_Veg_MC_Std <- exilis_Veg %>% mutate(
        MCstd_Biomass = (Biomass - mean(Biomass, na.rm = TRUE)) / sd(Biomass, na.rm = TRUE),
        MCstd_Photo_DayGE = (Photo_DayGE - mean(Photo_DayGE, na.rm = TRUE)) / sd(Photo_DayGE, na.rm = TRUE),
        MCstd_Trans_DayGE = (Trans_DayGE - mean(Trans_DayGE, na.rm = TRUE)) / sd(Trans_DayGE, na.rm = TRUE),
        MCstd_WUE_DayGE = (WUE_DayGE - mean(WUE_DayGE, na.rm = TRUE)) / sd(WUE_DayGE, na.rm = TRUE),
        MC_Tleaf_DayGE = (Tleaf_DayGE - mean(Tleaf_DayGE, na.rm = TRUE)),
        MC_DayGE_Node = (DayGE_Node - mean(DayGE_Node, na.rm = TRUE)),
        MCstd_Photo_DayGE_squared = MCstd_Photo_DayGE^2,
        MCstd_Trans_DayGE_squared = MCstd_Trans_DayGE^2,
        MCstd_WUE_DayGE_squared = MCstd_WUE_DayGE^2
)
xantiana_Veg_MC_Std <- xantiana_Veg %>% mutate(
        MCstd_Biomass = (Biomass - mean(Biomass, na.rm = TRUE)) / sd(Biomass, na.rm = TRUE),
        MCstd_Photo_DayGE = (Photo_DayGE - mean(Photo_DayGE, na.rm = TRUE)) / sd(Photo_DayGE, na.rm = TRUE),
        MCstd_Trans_DayGE = (Trans_DayGE - mean(Trans_DayGE, na.rm = TRUE)) / sd(Trans_DayGE, na.rm = TRUE),
        MCstd_WUE_DayGE = (WUE_DayGE - mean(WUE_DayGE, na.rm = TRUE)) / sd(WUE_DayGE, na.rm = TRUE),
        MC_Tleaf_DayGE = (Tleaf_DayGE - mean(Tleaf_DayGE, na.rm = TRUE)),
        MC_DayGE_Node = (DayGE_Node - mean(DayGE_Node, na.rm = TRUE)),
        MCstd_Photo_DayGE_squared = MCstd_Photo_DayGE^2,
        MCstd_Trans_DayGE_squared = MCstd_Trans_DayGE^2,
        MCstd_WUE_DayGE_squared = MCstd_WUE_DayGE^2
)
parviflora_Veg_MC_Std <- parviflora_Veg %>% mutate(
        MCstd_Biomass = (Biomass - mean(Biomass, na.rm = TRUE)) / sd(Biomass, na.rm = TRUE),
        MCstd_Photo_DayGE = (Photo_DayGE - mean(Photo_DayGE, na.rm = TRUE)) / sd(Photo_DayGE, na.rm = TRUE),
        MCstd_Trans_DayGE = (Trans_DayGE - mean(Trans_DayGE, na.rm = TRUE)) / sd(Trans_DayGE, na.rm = TRUE),
        MCstd_WUE_DayGE = (WUE_DayGE - mean(WUE_DayGE, na.rm = TRUE)) / sd(WUE_DayGE, na.rm = TRUE),
        MC_Tleaf_DayGE = (Tleaf_DayGE - mean(Tleaf_DayGE, na.rm = TRUE)),
        MC_DayGE_Node = (DayGE_Node - mean(DayGE_Node, na.rm = TRUE)),
        MCstd_Photo_DayGE_squared = MCstd_Photo_DayGE^2,
        MCstd_Trans_DayGE_squared = MCstd_Trans_DayGE^2,
        MCstd_WUE_DayGE_squared = MCstd_WUE_DayGE^2
)

remove(unguiculata_Veg, unguiculata_Rep, exilis_Veg, exilis_Rep,
       parviflora_Veg, parviflora_Rep, xantiana_Veg, xantiana_Rep,
       Clarkia_Field_2010, Clarkia_Field_filtered)
