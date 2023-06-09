#!/usr/bin/env RScript

library("tidyverse")
library("readxl")
library("ggpubr")

##### Data import ####
#Products recommended in one recommendation

recom_resellers <- read_csv2("data/100_all_country_clean_reseller/retailer_surveys.csv")

recom_product_id <- read_csv2("data/03_product_names_cleaned/20_recom_with_id/all_countries_recom_commercial_names_ids.csv",
                              col_types = cols(.default = "c")) %>%
    mutate(chemical_quantity = sub(",", ".", chemical_quantity, fixed = TRUE))

sheet_products_ec <- read_xlsx("data/30_extract_pest_dose_infos_from_TS/01_EC_products_pest_to_be_enriched_WIP_v2.xlsx") %>% 
    select(product_id, 
           pest_id, 
           crop,
           crop_found,
           pest_found,
           ts_product_qty_min, 
           ts_product_qty_max,
           ts_product_unit,
           ts_denomiator_qty,
           ts_denominator_unit) %>% mutate(pest_found = as.numeric(pest_found), country = "ecuador")

sheet_products_pe <- read_xlsx("data/30_extract_pest_dose_infos_from_TS/01_PE_products_pest_to_be_enriched_WIP_4.xlsx") %>% 
    select(product_id, 
           pest_id, 
           crop,
           crop_found,
           pest_found,
           ts_product_qty_min, 
           ts_product_qty_max,
           ts_product_unit,
           ts_denomiator_qty,
           ts_denominator_unit) %>% 
    mutate(country = "peru")

sheet_products_bo <- read_xlsx("data/30_extract_pest_dose_infos_from_TS/01_BO_products_pest_to_be_enriched_WIP.xlsx") %>% 
    select(product_id, 
           pest_id, 
           crop,
           crop_found,
           pest_found,
           ts_product_qty_min, 
           ts_product_qty_max,
           ts_product_unit,
           ts_denomiator_qty,
           ts_denominator_unit) %>% 
    mutate(country = "bolivia")

sheet_products_all <- bind_rows(sheet_products_ec, sheet_products_pe, sheet_products_bo)


##### Filter data ####

sheet_products_all_nona <- sheet_products_all %>% filter(!is.na(crop_found), !is.na(pest_found))

# Data at the reseller's recommendation level 
recom_dosis <- recom_product_id %>% 
    select(`_submission__uuid`, product_id, chemical_quantity, chem_qty_units, 
           water_quantity, water_unit, frec_num_veces, frecuencia_aplicacion)

recom_dosis_complete <- recom_dosis %>% 
    filter(!is.na(chemical_quantity) & !is.na(chem_qty_units) &
           !is.na(water_quantity) & !is.na(water_unit) &
           !is.na(frec_num_veces) & !is.na(frecuencia_aplicacion))

recom_dosis_all_ids <- recom_resellers %>% select(`_uuid`, crop, pest_cplx_id) %>% 
    inner_join(recom_dosis_complete, by = c("_uuid" = "_submission__uuid"))


# Filter products with dosis corresponding to pests in their technical sheets
sheet_prods_with_dosis <- sheet_products_all_nona %>%
    filter(!is.na(ts_product_qty_min) & 
             !is.na(ts_product_qty_max) & 
             !is.na(ts_product_unit)&
             !is.na(ts_denomiator_qty)&
             !is.na(ts_denominator_unit)
             )

## Combine data from resellers with data from technical sheets

dosis_recom_ts <- recom_dosis_all_ids %>% inner_join(sheet_prods_with_dosis, by = c("product_id", "crop", "pest_cplx_id" = "pest_id"))

##### Conversion from local units to international units ##### 

unit_conversion <- read_xlsx("data/02_equivalence_synonyms/all_unit_conversion.xlsx")

## Should be done at the country level because variable across countries

## Ecuador
unit_conv_ec <- unit_conversion %>% filter(country == "ecuador") %>% 
    select(unit, coef_to_ml)

dosis_ml_ec <- dosis_recom_ts %>% filter(country == "ecuador") %>% 
    left_join(unit_conv_ec, by = c("chem_qty_units" = "unit")) %>% 
    mutate(chem_quantity_ml = as.numeric(chemical_quantity) * as.numeric(coef_to_ml))

dosis_ml_ec_2 <- dosis_ml_ec %>% rename("coef_to_ml_resell" = "coef_to_ml") %>% 
    left_join(unit_conv_ec, by = c("ts_product_unit" = "unit")) %>% 
    mutate(ts_product_qty_min_ml = ts_product_qty_min * as.numeric(coef_to_ml),
           ts_product_qty_max_ml = ts_product_qty_max * as.numeric(coef_to_ml))
## Peru 
unit_conv_pe <- unit_conversion %>% filter(country == "peru") %>% 
    select(unit, coef_to_ml)

dosis_ml_pe <- dosis_recom_ts %>% filter(country == "peru") %>% 
    left_join(unit_conv_pe, by = c("chem_qty_units" = "unit")) %>% 
    mutate(chem_quantity_ml = as.numeric(chemical_quantity) * as.numeric(coef_to_ml))

# Error in one data sheet --> mg instead of g
dosis_ml_pe <- dosis_ml_pe %>% mutate(ts_product_unit = ifelse(product_id == "pe-7", "g", ts_product_unit))

dosis_ml_pe_2 <- dosis_ml_pe %>% rename("coef_to_ml_resell" = "coef_to_ml") %>% 
    left_join(unit_conv_pe, by = c("ts_product_unit" = "unit")) %>% 
    mutate(ts_product_qty_min_ml = ts_product_qty_min * as.numeric(coef_to_ml),
           ts_product_qty_max_ml = ts_product_qty_max * as.numeric(coef_to_ml))

## Bolivia 
unit_conv_bo <- unit_conversion %>% filter(country == "bolivia") %>% 
    select(unit, coef_to_ml)

dosis_ml_bo <- dosis_recom_ts %>% filter(country == "bolivia") %>% 
    left_join(unit_conv_bo, by = c("chem_qty_units" = "unit")) %>% 
    mutate(chem_quantity_ml = as.numeric(chemical_quantity) * as.numeric(coef_to_ml))

dosis_ml_bo_2 <- dosis_ml_bo %>% rename("coef_to_ml_resell" = "coef_to_ml") %>% 
    left_join(unit_conv_bo, by = c("ts_product_unit" = "unit")) %>% 
    mutate(ts_product_qty_min_ml = ts_product_qty_min * as.numeric(coef_to_ml),
           ts_product_qty_max_ml = ts_product_qty_max * as.numeric(coef_to_ml))

## Combine data converted across countries
dosis_ml_all <- bind_rows(dosis_ml_ec_2, dosis_ml_pe_2, dosis_ml_bo_2)

# Remove TS info with ha
dosis_ml_all <- dosis_ml_all %>% filter(ts_denominator_unit != "ha")

# Compute dilution ratios
dilution_ratio_all <- dosis_ml_all %>% 
    mutate(chem_qty_l = chem_quantity_ml / 1000, # to transform in liters
           dilution_ratio_reseller = chem_qty_l/as.numeric(water_quantity),
           chem_qty_l_min = ts_product_qty_min_ml / 1000, # to transform in liters
           dilution_ratio_ts_min = chem_qty_l_min/as.numeric(ts_denomiator_qty),
           chem_qty_l_max = ts_product_qty_max_ml / 1000, # to transform in liters
           dilution_ratio_ts_max = chem_qty_l_max/as.numeric(ts_denomiator_qty) )

dilution_categories <- dilution_ratio_all %>% mutate(correct_dose = ifelse(
    dilution_ratio_reseller <= dilution_ratio_ts_max & dilution_ratio_reseller >= dilution_ratio_ts_min,
    "correct",
    ifelse(
        dilution_ratio_reseller > dilution_ratio_ts_max,
        "overuse",
         ifelse(dilution_ratio_reseller < dilution_ratio_ts_min,
        "underuse", NA)
    )))

dilution_categories_normalized <- dilution_categories %>% 
    mutate(mean_ts_dilution = (dilution_ratio_ts_min + dilution_ratio_ts_max)/2,
           normalized_dilution_reseller = dilution_ratio_reseller/mean_ts_dilution,
           normalized_dilution_min = dilution_ratio_ts_min/mean_ts_dilution,
           normalized_dilution_max = dilution_ratio_ts_max/mean_ts_dilution)

prop_cat <- dilution_categories_normalized %>% group_by(correct_dose) %>% 
    mutate(total = nrow(.)) %>% summarize(nb = n(), total = mean(total)) %>% 
    mutate(prop = (nb/total)*100, x = "1")

percent_cat <- ggplot(prop_cat, mapping = aes(x = x, y =prop, fill = correct_dose, group = reorder(correct_dose, c(2,1,3))))+
    geom_bar(stat = "identity", position="stack", width = 0.8, alpha = 0.5)+
    geom_text(aes(label=paste0(nb, " \\n (",sprintf("%1.1f", prop),"%)")), 
              position=position_stack(vjust=0.5), colour="black") + 
    scale_fill_manual(values = c("#00a261","#d71e36" , "#ff7f27"))+
    theme_void()+
    theme(legend.position="bottom")

barplot<-ggplot(dilution_categories_normalized, aes(x = normalized_dilution_reseller, fill = correct_dose)) +
    geom_histogram(alpha = 0.5, position = "identity", binwidth = 0.05) +
    scale_x_log10(breaks = c(1/8, 1/4, 1/3, 1/2, 1/1.5, 1, 1.5, 2,3,4,8), labels = c("÷8", "÷4","÷3","÷2","÷1.5", "1", "x1.5", "x2", "x3", "x4", "x8" ))+
    scale_fill_manual(values = c("#00a261","#d71e36" , "#ff7f27"))+
    theme_minimal()+
    labs(x = "Normalized dilution ratio (log)")  +  theme(legend.position="none")

total_plot <- ggarrange(barplot,percent_cat, ncol = 2, widths = c(7,1))

ggsave(total_plot, filename = "output/Fig_2_underuse_overuse.pdf", device = "pdf", width = 10, height = 6)

## Export plot data
dilution_categories_normalized_export <- dilution_categories_normalized %>% 
    select(`_uuid`, product_id, normalized_dilution_reseller, correct_dose) 

proportion_dose_type_export <- prop_cat %>% select(-x) %>% 
    rename(dose_type = correct_dose, nb_recommendations = nb,
           total_recommendations = total, proportion = prop)

