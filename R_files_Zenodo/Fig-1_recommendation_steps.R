#!/usr/bin/env RScript

library("tidyverse")
library("readxl")
library("xlsx")

# If several products are recommended by the reseller for the same pest, what is the proportion that should be appropriate ?
TRESHOLD_ADEQUATE_PRODUCT_TYPE <- 1
PROPORTION_PRODUCTS_PER_RECOM_TO_BE_RIGHT <- 1
PROPORTION_DOSE_PER_RECOM_TO_BE_RIGHT <- 1

#### STEP 1 : evaluating pest identification by resellers ####

## Load pest synonyms list
pest_synonyms <- read_xlsx("data/40_pest_synonyms/all_countries_pest_synonyms_v3.xlsx")


pest_kobo_id <- pest_synonyms  %>% distinct(pest_id, pest_name_kobo)

## Load pest names harmonized
pest_names_harmonized_bo <- read_xls("data/05_pest_names_clean/BO_raw_pest_names.xls") 
pest_names_harmonized_ec <- read_xls("data/05_pest_names_clean/EC_raw_pest_names.xls") 
pest_names_harmonized_pe <- read_xls("data/05_pest_names_clean/PE_raw_pest_names.xls") 

## Load all recomendations by resellers with pest ids
recoms <- read_csv2("data/100_all_country_clean_reseller/retailer_surveys.csv")


recoms_names_harm <- recoms %>% left_join(pest_names_harmonized_bo, by = "_uuid") %>%
    mutate(id_pest = ifelse(!is.na(id_plaga), id_plaga, id_pest)) %>% 
    select(-id_plaga) %>% 
    left_join(pest_names_harmonized_ec, by = "_uuid") %>%
    mutate(id_pest = ifelse(!is.na(id_plaga), id_plaga, id_pest)) %>% 
    select(-id_plaga) %>% 
    left_join(pest_names_harmonized_pe, by = "_uuid") %>%
    mutate(id_pest = ifelse(!is.na(id_plaga), id_plaga, id_pest)) %>% 
    select(-id_plaga)

## Add pest ids (per country) to the recommendations based on the kobo 
## Here the same pests species have different IDs depending on the country because
## pests synonyms are different depending on the region
recoms_with_pest_id <- recoms_names_harm %>% left_join(pest_kobo_id, by = c("pest_label_kobo" = "pest_name_kobo")) %>% 
    mutate(country = ifelse(
        country == "ecuador", "EC", 
        ifelse(country == "peru", "PE",
               ifelse(country == "bolivia", "BO", ""))
    ))

dl <- recoms_with_pest_id %>%
    count(country, pest_label_kobo, pest_id)

## Join the synonyms with the pests accounting for differences between countries
recom_with_synonymy <- recoms_with_pest_id %>% 
    left_join(pest_synonyms, by = c("id_pest" = "pest_name_given_by_reseller","pest_label_kobo" = "pest_name_kobo", "country"))

## Identifiy the pest synonyms for which local expert are unsure about the synonymy
## To compute the "unknown" barplot section

unknown_synonymy <- recom_with_synonymy %>% filter(synonymy_harmonized_qs == "unknown" | synonymy_harmonized_qs == "error") %>% group_by(`_uuid`) %>% 
    summarize() %>% mutate(no_synonymy = 1)

recom_with_synonymy_bin <- recom_with_synonymy %>% 
    mutate(synonymy_harmonized_qs = ifelse(synonymy_harmonized_qs == "yes", 1, 0))

recom_synonymy_checked <- recom_with_synonymy_bin %>% group_by(`_uuid`) %>%
    summarize(synonymy_mean = mean(synonymy_harmonized_qs)) %>% left_join(unknown_synonymy,by = "_uuid") %>% 
    mutate(synonymy_mean = ifelse(is.na(no_synonymy), synonymy_mean, -1)) %>% select(-no_synonymy) %>% 
    mutate(pest_identification = ifelse(synonymy_mean > 0, "right", ifelse(
        synonymy_mean == 0, "wrong", ifelse(
            synonymy_mean == -1, "unknown", NA
        )
    )))

pest_identification <- recom_synonymy_checked %>% group_by(pest_identification) %>% summarize(n = n())

# Create list of uuid with right answers to continue the decision tree
recom_correctly_id_pest <- recom_synonymy_checked %>%  mutate(
    synonymy_mean = ifelse(is.na(synonymy_mean), 0, synonymy_mean) # NA corresponds to no answer given by the reseller --> false answer (0)
) %>% filter(synonymy_mean > 0) %>% select(`_uuid`)

# write.xlsx2(recom_correctly_id_pest, "data/40_pest_synonyms/50_ALL_recom_id_corrrectly_identified_pest.xlsx")
# write.xlsx2(recom_synonymy_checked, "data/40_pest_synonyms/50_ALL_recom_id_pest.xlsx")



#### STEP 2 : evaluating adequacy between product type and pest type ####

#Products recommended in one recommendation
recom_product_id <- read_csv2("data/03_product_names_cleaned/20_recom_with_id/all_countries_recom_commercial_names_ids.csv",
                              col_types = cols(.default = "c"))



pest_phylum <- read_xlsx("data/40_pest_synonyms/30_ALL_pest_complex_id_common_names.xlsx") %>% 
    select(pest_complex_id, phylum)

recom_pest_id <- read_csv2("data/100_all_country_clean_reseller/retailer_surveys.csv") %>% 
    select(`_uuid`, pest_cplx_id, country)

recom_with_unknown_products <- recom_pest_id %>% anti_join(recom_product_id, by = c("_uuid" = "_submission__uuid"))
recom_with_unknown_products_formatted <- recom_with_unknown_products %>% select(`_uuid`) %>% 
    mutate(product_adequate = "unknown")

# Compare product_type with pest phyum for each product
correct_product_for_pest <- recom_product_id %>% inner_join(recom_pest_id, by = c("_submission__uuid" = "_uuid")) %>% 
    left_join(pest_phylum, by = c("pest_cplx_id" = "pest_complex_id")) %>% 
    select(`_submission__uuid`, type, pest_cplx_id, phylum) %>% 
    mutate(phylum = ifelse(phylum == "fungus",
                           "fungicide",
                           ifelse(phylum == "insect", "insecticide", NA))) %>% 
    mutate(correct_product = ifelse(type == phylum, 1, 0))

# Group responses by recommendation (each recommendation may have several products)
correct_prod_pest_per_recom <-correct_product_for_pest %>% group_by(`_submission__uuid`) %>% 
    mutate(total_products = n()) %>% 
    summarize(correct_prod_per_recom = sum(correct_product),
              total_products_per_recom = mean(total_products),
              proportion_correct_prod_per_recom = correct_prod_per_recom/total_products_per_recom)


# Qualify recommendations as being right if all products correspond to the phylum

adequate_prod_recom <- correct_prod_pest_per_recom %>% 
    mutate(product_adequate = ifelse(proportion_correct_prod_per_recom >= TRESHOLD_ADEQUATE_PRODUCT_TYPE, "right", "wrong")) %>% 
    select(`_submission__uuid`, product_adequate) %>% rename("_uuid" = "_submission__uuid") %>% 
    bind_rows(recom_with_unknown_products_formatted) %>% 
    # Keep only the reommendations with the correct pest identification by reseller
    inner_join(recom_correctly_id_pest, by = c("_uuid")) 

product_type_adequacy <- adequate_prod_recom %>% group_by(product_adequate) %>% summarize(n = n())

# Create list of uuid with right answers to continue with the decision tree
recom_correct_prod_type <- adequate_prod_recom %>%  filter(product_adequate == "right") %>% 
    select(`_uuid`)


#### STEP 3 : evaluating adequacy between pest species and product data sheet ####


pool_products_with_same_ingr <- function(product_list){
    pooled_products <- product_list %>% select(product_id, ingr_CAS, ingr_id) %>% 
        pivot_wider(id_cols = product_id, names_from = ingr_id, values_from = ingr_CAS, names_prefix = "ingr_") %>% 
        group_by(ingr_1, ingr_2) %>%
        mutate(pooled_product_id = paste0("pooled-", group_indices()))
    return(pooled_products)
}


product_unq <- read_xls("data/04_pesticide_list_ingredients/04_ingr_homogenized_enriched/homogenized-active-ingredients-all-countries.xls")

recom_pest_id_ok_prev_step <- read_csv2("data/100_all_country_clean_reseller/retailer_surveys.csv",
                                        col_types = cols(.default = "c")) %>%  
    inner_join(recom_correct_prod_type, by = "_uuid")

## PERU

#Products
product_pe <- product_unq %>% filter(country == "PE") 
pooled_prod_pe <- pool_products_with_same_ingr(product_list = product_pe)

#Technical sheet info from product
pests_pe_product <- read_xlsx("data/30_extract_pest_dose_infos_from_TS/01_PE_products_pest_to_be_enriched_WIP_4.xlsx")
products_pooled_pest_pe <- pests_pe_product %>% inner_join(pooled_prod_pe, by = "product_id")

# Recommendations by resellers
recom_prod_type_ok_pe <- recom_pest_id_ok_prev_step %>% filter(country == "peru")
recom_to_prod <- recom_product_id %>% select(`_submission__uuid`, product_id)

## Combine the different levels of data and identify the products-crop-pest without data sheet
recom_prod_ts <- recom_prod_type_ok_pe %>% left_join(recom_to_prod, by = c("_uuid" = "_submission__uuid")) %>% 
    left_join(products_pooled_pest_pe, by = c("product_id", "crop", "pest_cplx_id" = "pest_id")) %>% 
    mutate(ts_not_found = ifelse(is.na(crop_found) | is.na(pest_found), 1, 0))

## Determine if correct product based on pooled products with same ingredients
recom_prod_ts_pooled <- recom_prod_ts %>% group_by(product_id, crop, pest_cplx_id) %>% 
    mutate(crop_found_pooled = sum(crop_found),
           pest_found_pooled = sum(pest_found),
           right_product = ifelse(crop_found_pooled > 0 & pest_found_pooled > 0, 1, 0))

# Determine the neglected pests for Peru
recom_prod_ts_pooled %>% group_by(pest_cplx_id) %>% summarize(total_right_products = sum(right_product, na.rm = TRUE)) ## No neglected pests for Peru

ts_not_found <- recom_prod_ts_pooled %>% group_by(`_uuid`) %>% summarize(ts_not_found_recom = mean(ts_not_found))

recom_without_ts <- recom_prod_ts_pooled %>%  group_by(`_uuid`) %>% 
    summarize(ts_not_found_recom = mean(ts_not_found)) %>% 
    filter(ts_not_found_recom > 0)

recom_prod_ts_pooled_summ <- recom_prod_ts_pooled %>% ungroup %>%
    filter(!`_uuid` %in% recom_without_ts$`_uuid`) %>%  group_by(`_uuid`) %>% 
    summarise(nb_prod = n(),
              right_products_prop = sum(right_product)/nb_prod)

recom_prod_ts_ok <- recom_prod_ts_pooled_summ %>% 
    mutate(type = ifelse(right_products_prop >= PROPORTION_PRODUCTS_PER_RECOM_TO_BE_RIGHT, "right", "wrong")) %>% 
    select(`_uuid`, type)

prod_datasheet_adequacy_pe <- recom_without_ts %>% mutate(type = "unknown") %>%
    bind_rows(recom_prod_ts_ok) %>% select(-ts_not_found_recom)

# write.xlsx2(recom_prod_ts_pooled_summ, "data/50_fig_1_steps/30_PE_correct_prod_TS_proportion.xlsx")


## ECUADOR

#Products
product_ec <- product_unq %>% filter(country == "EC") 
pooled_prod_ec <- pool_products_with_same_ingr(product_list = product_ec)

#Technical sheet info from product
pests_ec_product <- read_xlsx("data/30_extract_pest_dose_infos_from_TS/01_EC_products_pest_to_be_enriched_WIP_v2.xlsx")
products_pooled_pest_ec <- pests_ec_product %>% inner_join(pooled_prod_ec, by = "product_id")

# Recommendations by resellers
recom_prod_type_ok_ec <- recom_pest_id_ok_prev_step %>% filter(country == "ecuador")
recom_to_prod <- recom_product_id %>% select(`_submission__uuid`, product_id)

## Combine the different levels of data and identify the products-crop-pest without data sheet
recom_prod_ts <- recom_prod_type_ok_ec %>% left_join(recom_to_prod, by = c("_uuid" = "_submission__uuid")) %>% 
    left_join(products_pooled_pest_ec, by = c("product_id", "crop", "pest_cplx_id" = "pest_id")) %>% 
    mutate(ts_not_found = ifelse(is.na(crop_found) | is.na(pest_found), 1, 0))

## Determine if correct product based on pooled products with same ingredients
recom_prod_ts_pooled <- recom_prod_ts %>%
    mutate(pest_found = ifelse(pest_found == "spodoptera", 1, pest_found)) %>% 
    group_by(product_id, crop, pest_cplx_id) %>% 
    mutate(crop_found_pooled = sum(as.numeric(crop_found), na.rm = T),
           pest_found_pooled = sum(as.numeric(pest_found),  na.rm = T),
           right_product = ifelse(crop_found_pooled > 0 & pest_found_pooled > 0, 1, 0))

# Determine the neglected pests for Ecuador
neglected_pests_id_ec <- recom_prod_ts_pooled %>% group_by(pest_cplx_id) %>% mutate(total_right_products = sum(right_product, na.rm = TRUE)) %>% 
    filter(total_right_products == 0) %>% distinct(pest_cplx_id) %>% 
    mutate(neglected_pests = 1)

## Change neglected pests as unknown instead of wrong

ts_not_found <- recom_prod_ts_pooled %>% group_by(`_uuid`) %>% summarize(ts_not_found_recom = mean(ts_not_found))

recom_with_neglected_pests <- recom_prod_ts_pooled %>% inner_join(neglected_pests_id_ec, by = "pest_cplx_id") %>% 
    group_by(`_uuid`) %>% 
    summarize(neglected_pests_recom = mean(neglected_pests)) %>% distinct(`_uuid`)

recom_without_ts <- recom_prod_ts_pooled %>%  group_by(`_uuid`) %>% 
    summarize(ts_not_found_recom = mean(ts_not_found)) %>% 
    filter(ts_not_found_recom > 0) %>% distinct(`_uuid`)

recom_without_ts_with_neglected <- bind_rows(recom_without_ts, recom_with_neglected_pests) %>% distinct()

## Right and wrong classification
recom_prod_ts_pooled_summ <- recom_prod_ts_pooled %>% ungroup %>%
    filter(!`_uuid` %in% recom_without_ts_with_neglected$`_uuid`) %>%  group_by(`_uuid`) %>% 
    summarise(nb_prod = n(),
              right_products_prop = sum(right_product)/nb_prod)

recom_prod_ts_ok <- recom_prod_ts_pooled_summ %>% 
    mutate(type = ifelse(right_products_prop >= PROPORTION_PRODUCTS_PER_RECOM_TO_BE_RIGHT, "right", "wrong")) %>% 
    select(`_uuid`, type)

## Add unknown cat
prod_datasheet_adequacy_ec <- recom_without_ts_with_neglected %>% mutate(type = "unknown") %>%
    bind_rows(recom_prod_ts_ok) 



# write.xlsx2(recom_prod_ts_pooled_summ, "data/50_fig_1_steps/30_EC_correct_prod_TS_proportion.xlsx")

## BOLIVIA

#Products
product_bo <- product_unq %>% filter(country == "BO") 
pooled_prod_bo <- pool_products_with_same_ingr(product_list = product_bo)

#Technical sheet info from product
pests_bo_product <- read_xlsx("data/30_extract_pest_dose_infos_from_TS/01_BO_products_pest_to_be_enriched_WIP.xlsx")
products_pooled_pest_bo <- pests_bo_product %>% inner_join(pooled_prod_bo, by = "product_id")

# Recommendations by resellers
recom_prod_type_ok_bo <- recom_pest_id_ok_prev_step %>% filter(country == "bolivia")
recom_to_prod <- recom_product_id %>% select(`_submission__uuid`, product_id)

## Combine the different levels of data and identify the products-crop-pest without data sheet
recom_prod_ts <- recom_prod_type_ok_bo %>% left_join(recom_to_prod, by = c("_uuid" = "_submission__uuid")) %>% 
    left_join(products_pooled_pest_bo, by = c("product_id", "crop", "pest_cplx_id" = "pest_id")) %>% 
    mutate(ts_not_found = ifelse(is.na(crop_found) | is.na(pest_found), 1, 0))

## Determine if correct product based on pooled products with same ingredients
recom_prod_ts_pooled <- recom_prod_ts %>% group_by(product_id, crop, pest_cplx_id) %>% 
    mutate(crop_found_pooled = sum(crop_found),
           pest_found_pooled = sum(pest_found),
           right_product = ifelse(crop_found_pooled > 0 & pest_found_pooled > 0, 1, 0))

# Determine the neglected pests for Bolivia
neglected_pests_id_bo <- recom_prod_ts_pooled %>% group_by(pest_cplx_id) %>% summarize(total_right_products = sum(right_product, na.rm = TRUE)) %>% 
    filter(total_right_products == 0) ## No neglected pests for Bolivia

ts_not_found <- recom_prod_ts_pooled %>% group_by(`_uuid`) %>% summarize(ts_not_found_recom = mean(ts_not_found))

recom_without_ts <- recom_prod_ts_pooled %>%  group_by(`_uuid`) %>% 
    summarize(ts_not_found_recom = mean(ts_not_found)) %>% 
    filter(ts_not_found_recom > 0)

recom_prod_ts_pooled_summ <- recom_prod_ts_pooled %>% ungroup %>%
    filter(!`_uuid` %in% recom_without_ts$`_uuid`) %>%  group_by(`_uuid`) %>% 
    summarise(nb_prod = n(),
              right_products_prop = sum(right_product)/nb_prod)

recom_prod_ts_ok <- recom_prod_ts_pooled_summ %>% 
    mutate(type = ifelse(right_products_prop >= PROPORTION_PRODUCTS_PER_RECOM_TO_BE_RIGHT, "right", "wrong")) %>% 
    select(`_uuid`, type)

prod_datasheet_adequacy_bo <- recom_without_ts %>% mutate(type = "unknown") %>%
    bind_rows(recom_prod_ts_ok) %>% select(-ts_not_found_recom)

# write.xlsx2(recom_prod_ts_pooled_summ, "data/50_fig_1_steps/30_BO_correct_prod_TS_proportion.xlsx")


## Merge results from all countries

prod_datasheet_adequacy_all <- bind_rows(prod_datasheet_adequacy_pe, prod_datasheet_adequacy_ec, prod_datasheet_adequacy_bo) 

prod_datasheet_ok <- prod_datasheet_adequacy_all %>% filter(type == "right") %>% select(`_uuid`)

prod_datasheet_adequacy <- prod_datasheet_adequacy_all %>% group_by(type) %>% summarise(n = n())


#### STEP 4 : evaluating adequacy between dose proposed by resellers and provided by data sheets ####
source("Fig-2_dosage.R")

# Keep dilution results only for recommendations tha have correctly passed all the previous steps
dose_categories <- dilution_categories %>% filter(`_uuid` %in% prod_datasheet_ok$`_uuid`)

dose_categories_binary <- dose_categories %>% mutate(correct_dose = ifelse(correct_dose != "correct", 0, 1))

dose_categories_binary_summarized <- dose_categories_binary %>%  group_by(`_uuid`) %>% summarize(correct_dose = mean(correct_dose)) %>% 
    mutate(dose_type = ifelse(correct_dose >= PROPORTION_DOSE_PER_RECOM_TO_BE_RIGHT, "right", "wrong"))

right_dose_recom <- dose_categories_binary_summarized %>%  
    distinct(`_uuid`)

unknown_doses <- prod_datasheet_ok %>% anti_join(right_dose_recom, by = "_uuid") %>% 
    mutate(dose_type = "unknown")

all_doses <- dose_categories_binary_summarized %>% select(-correct_dose) %>% bind_rows(unknown_doses) %>% 
    group_by(dose_type) %>% summarise(n = n())


## COMBINE all recommendations and their stauts at each recom step
step_1_recom <- recom_synonymy_checked %>% mutate(step = "pest_id", answer_type = pest_identification) %>% 
    select(`_uuid`, step, answer_type)

step_2_recom <- adequate_prod_recom %>% mutate(step = "product_type", answer_type = product_adequate) %>% 
    select(`_uuid`, step, answer_type)

step_3_recom <- prod_datasheet_adequacy_all %>% mutate(step = "product_datasheet", answer_type = type) %>% 
    select(`_uuid`, step, answer_type)

step_4_recom <- dose_categories_binary_summarized %>% bind_rows(unknown_doses)%>% mutate(step = "dose_datasheet", answer_type = dose_type) %>% 
    select(`_uuid`, step, answer_type)

allrecom_allsteps <- bind_rows(step_1_recom, step_2_recom, step_3_recom, step_4_recom)
# write_csv(allrecom_allsteps, path = "data/fig_1_all_steps_all_types.csv")

#### FINAL STEP : Plot all barplots ####

step_1 <- pest_identification %>% mutate(step = "pest_id") %>% rename(answer_type = pest_identification)
step_2 <- product_type_adequacy %>% mutate(step = "product_type")%>% rename(answer_type = product_adequate)
step_3 <- prod_datasheet_adequacy %>% mutate(step = "product_datasheet")%>% rename(answer_type = type)
step_4 <- all_doses %>% mutate(step = "dose_datasheet")%>% rename(answer_type = dose_type)


## Because we do not account for unknown responses in the proportion computation
data_plot_known <- bind_rows(step_1, step_2, step_3, step_4) %>% 
    filter(answer_type != "unknown") %>% 
    group_by(step) %>% mutate(total = sum(n), prop = n/total)
    
    
data_plot_unknown <- bind_rows(step_1, step_2, step_3, step_4) %>% 
    filter(answer_type == "unknown") %>% 
    group_by(step) %>% mutate(prop = 0)

data_plot <- bind_rows(data_plot_known, data_plot_unknown)


plot_a <- ggplot(data = data_plot, 
                 mapping = aes(x = factor(step, levels=c("pest_id", "product_type", "product_datasheet", "dose_datasheet")), 
                               y = n , 
                               fill = factor(answer_type, levels=rev(c("right", "wrong", "unknown")))))+
    geom_bar(stat="identity", position = "stack", width = 0.4)+
    geom_text(aes(label = paste0(n, " \\n ", round(prop*100,1),"%")),position=position_stack(vjust=0.5), colour="black") +
    scale_fill_manual(values = rev(c("#7fd0b0ff",  "#eb8e9aff","grey80")))+
    scale_x_discrete(position = "top") +    theme_minimal()+
    xlab("")+ ylab("Number of recommendation")+
    scale_y_reverse()+
    theme(legend.position = "none")


data_plot_merged <- data_plot %>% group_by(answer_type) %>% summarize(n = sum(n)) %>% 
    mutate(n = ifelse(answer_type == "right",pull(data_plot[data_plot$answer_type == "right" & data_plot$step == "dose_datasheet", "n"]), n))  # because we only take the last step to compute the right as they are cumulated

data_plot_merged_uk <- data_plot_merged %>% filter(answer_type == "unknown") %>% 
    mutate(total = sum(n), prop = 0, step = "Overall")
data_plot_merged_knw <- data_plot_merged %>% filter(answer_type != "unknown") %>% 
    mutate(total = sum(n), prop = n/total, step = "Overall")

data_plot_merged_corr <- bind_rows(data_plot_merged_uk, data_plot_merged_knw)

plot_b <- ggplot(data = data_plot_merged_corr, 
                 mapping = aes(x = step, y = n , 
                               fill = factor(answer_type, levels=rev(c("right", "wrong", "unknown")))))+
    geom_bar(stat="identity", position = "stack", width = 0.4)+
    geom_text(aes(label = paste0(n, " \\n ", round(prop*100,1),"%")),position=position_stack(vjust=0.5), colour="black") +
    scale_fill_manual(values = rev(c("#7fd0b0ff",  "#eb8e9aff","grey80")))+
    scale_x_discrete(position = "top") +    theme_minimal()+
    xlab("")+ ylab("Number of recommendation")+
    scale_y_reverse()+
    theme(legend.position = "none")

plot_complete <- ggarrange(plot_a, plot_b, ncol = 2, widths = c(6,2))

ggsave(plot_complete, filename = "output/fig_1_all_stepsv2.pdf", device = "pdf", width = 10, height = 7)


