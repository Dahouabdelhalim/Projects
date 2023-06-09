#!/usr/bin/env RScript

library("tidyverse")
library("readxl")

#Recommandations by reseller
recom_resellers <- read_csv2("data/100_all_country_clean_reseller/retailer_surveys.csv")

#Products recmmended in one recommendation
recom_resellers_products <- read_csv2("data/03_product_names_cleaned/20_recom_with_id/all_countries_recom_commercial_names_ids.csv",
                                      col_types = cols(.default = "c")) 

# Homogenize toxic categories across countries for recommendations
recom_prod_toxcat <- recom_resellers_products %>% mutate(toxic_cat_homogenized = substr(toxic_cat, start = 0, stop = 3),
                                                         toxic_cat_homogenized = trimws(toxic_cat_homogenized))
recom_prod_toxcat_nona <- recom_prod_toxcat %>% filter(!is.na(toxic_cat_homogenized)) %>% select(-country)

# Homogenize toxic categories across countries for list of unique products
products_unq <- read_xls("data/04_pesticide_list_ingredients/04_ingr_homogenized_enriched/homogenized-active-ingredients-all-countries.xls")%>% 
    select(product_id, toxic_cat) %>% 
    mutate(toxic_cat_homogenized = substr(toxic_cat, start = 0, stop = 3),
           toxic_cat_homogenized = trimws(toxic_cat_homogenized))

products_unq_toxcat_nona <- products_unq %>% filter(!is.na(toxic_cat_homogenized))

##### Model 1 - Expert model ####

## Create list of crops and pests to be used as the base for models
recom_sample_for_model <- recom_resellers %>% inner_join(recom_prod_toxcat_nona, by = c("_uuid" = "_submission__uuid"))

sheet_products_ec <- read_xlsx("data/30_extract_pest_dose_infos_from_TS/01_EC_products_pest_to_be_enriched_WIP_v2.xlsx") %>% 
    select(product_id, pest_id, crop, crop_found, pest_found) %>% 
    mutate(pest_found = as.numeric(pest_found), country = "ecuador")
sheet_products_pe <- read_xlsx("data/30_extract_pest_dose_infos_from_TS/01_PE_products_pest_to_be_enriched_WIP_4.xlsx") %>% 
    select(product_id, pest_id, crop, crop_found, pest_found) %>% 
    mutate(country = "peru")
sheet_products_bo <- read_xlsx("data/30_extract_pest_dose_infos_from_TS/01_BO_products_pest_to_be_enriched_WIP.xlsx") %>% 
    select(product_id, pest_id, crop, crop_found, pest_found)%>% 
    mutate(country = "bolivia")

sheet_products_all <- bind_rows(sheet_products_ec, sheet_products_pe, sheet_products_bo)
sheet_products_all_nona <- sheet_products_all %>% filter(!is.na(crop_found), !is.na(pest_found))
right_products_all <- sheet_products_all_nona %>% filter(crop_found == 1 & pest_found == 1) %>% 
    inner_join(products_unq_toxcat_nona, by = "product_id")

right_products_all_num <- right_products_all %>% mutate(toxic_cat_numeric = ifelse(
    toxic_cat_homogenized == "ia", 1, ifelse(
        toxic_cat_homogenized == "ib", 1, ifelse(
            toxic_cat_homogenized == "ii", 2, ifelse(
                toxic_cat_homogenized == "iii", 3, ifelse(
                    toxic_cat_homogenized == "iv", 4, NA
                )
            )
        )
    )
))


expert_model_by_country <- function(recom_prod = recom_sample_for_model,
                                    country_sel = "ecuador", 
                                    right_products = right_products_all_num)
{
    recom_prod_country <- recom_prod %>% filter(country == country_sel) %>% 
        select(`_uuid`, crop, pest_cplx_id) 
    right_products_country <- right_products %>% filter(country == country_sel) 
    
    # Some crop/pest combination do not have right products because no data sheet mentioned them
    # so we need to get rid of those for comparison
    recom_right_prod <- recom_prod_country %>% 
        inner_join(right_products_country, by = c("crop", "pest_cplx_id" = "pest_id")) 
    
    expert_model_toxcat <- recom_right_prod %>% group_by(`_uuid`, crop, pest_cplx_id) %>% 
        filter(toxic_cat_numeric == max(toxic_cat_numeric)) %>% sample_n(1) %>%  #The higher the toxic category the less harmful; the better the recommendation
       ungroup %>%  select(`_uuid`,toxic_cat_homogenized)
    
    return(expert_model_toxcat)
}

expert_mod_output_bo <- expert_model_by_country(recom_prod = recom_sample_for_model,
                                                country_sel = "bolivia", 
                                                right_products = right_products_all_num)
expert_mod_output_ec <- expert_model_by_country(recom_prod = recom_sample_for_model,
                                                country_sel = "ecuador", 
                                                right_products = right_products_all_num)
expert_mod_output_pe <- expert_model_by_country(recom_prod = recom_sample_for_model,
                                                country_sel = "peru", 
                                                right_products = right_products_all_num)

expert_model_toxcat_prop_all <- bind_rows(expert_mod_output_bo, expert_mod_output_ec, expert_mod_output_pe) %>% 
    mutate(toxic_cat_homogenized = ifelse(toxic_cat_homogenized == "ia", "i", 
                                          ifelse(toxic_cat_homogenized == "ib", "i", toxic_cat_homogenized))) %>% 
    rename(toxcat_expert = toxic_cat_homogenized)


# To filter the observed dataset with the same samples from the expert model
recom_used_expert_uuid <-  recom_sample_for_model %>%  
    inner_join(right_products_all_num, by = c("crop", "pest_cplx_id" = "pest_id")) %>% 
    select(`_uuid`) %>% distinct()

##### Observed recommended toxic categories ####

observed_toxcat_prop <- recom_resellers %>% filter(`_uuid` %in% recom_used_expert_uuid$`_uuid`) %>% 
    inner_join(recom_prod_toxcat_nona, by = c("_uuid" = "_submission__uuid")) %>%
    select(`_uuid`, index, toxic_cat_homogenized) %>% 
    filter(index == 1) %>% ## Keep only the first product to match the expert model
    mutate(toxic_cat_homogenized = ifelse(toxic_cat_homogenized == "ia", "i", 
                                          ifelse(toxic_cat_homogenized == "ib", "i", toxic_cat_homogenized))) %>% 
    select(-index) %>% 
    rename(toxcat_obs = toxic_cat_homogenized)


##### Plot results alluvial####

all_models <- expert_model_toxcat_prop_all %>% inner_join(observed_toxcat_prop, by = "_uuid") %>% 
    group_by(toxcat_obs, toxcat_expert) %>% summarize(freq = n())

all_models_values <- all_models %>% ungroup %>% mutate(total = sum(freq), prop_percent = freq/total*100)

all_models_values_export <- all_models_values %>%
    rename(nb_recommendation = freq, 
           total_recommendations = total,
           proportion = prop_percent)

write.table(all_models_values, "output/table/Fig_4_expert_model_data.csv", sep =";", dec = ".", row.names = FALSE)

# Alluvial plot
library(alluvial)

colors_toxicity <- c(  "#d71e36", "#ff7f27", "#015ba9", "#00a261")

pdf(paste0("output/Fig-3_alluvial_chart.pdf"), width = 7, height = 12)
alluvial::alluvial(all_models[,1:2], freq = all_models$freq,
                   xw=0.2, gap.width=0.15,cw = 0.05,
                   col=ifelse(
                       all_models$toxcat_obs == "i", colors_toxicity[1], ifelse(
                           all_models$toxcat_obs == "ii", colors_toxicity[2], ifelse(
                               all_models$toxcat_obs == "iii", colors_toxicity[3],  ifelse(
                                   all_models$toxcat_obs == "iv", colors_toxicity[4], "black"
                           )))),border = NA,
                   blocks=TRUE, cex = 0.7, cex.axis = 1.5
)
dev.off()


all_models %>% group_by(toxcat_obs) %>% summarise(sum = sum(freq)) %>% 
    ungroup() %>% mutate(tot = sum (sum), perc = round((sum/tot)*100,1))


all_models %>% group_by(toxcat_expert) %>% summarise(sum = sum(freq)) %>% 
    ungroup() %>% mutate(tot = sum (sum), perc = round((sum/tot)*100, 2))

