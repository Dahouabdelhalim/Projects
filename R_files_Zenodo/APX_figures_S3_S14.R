#!/usr/bin/env RScript

library(tidyverse)
colors_org <- c( "black", "#d71e36", "#ff7f27", "#015ba9", "#00a261")

#Recommandations by reseller
recom_resellers <- read_csv2("data/100_all_country_clean_reseller/retailer_surveys.csv")

## Biopesticide proportion in recommendations

org_prop <- recom_resellers %>%
    count(organic_option_available) %>%
    filter(!is.na(.)) %>%
    mutate(n_tot = sum(n),
           perc = round(n/n_tot*100, 1))

org_prop$organic_option_available[org_prop$organic_option_available == "si"] <- "Yes"
org_prop$organic_option_available[org_prop$organic_option_available == "no"] <- "No"

orgplot <- ggplot(data = org_prop, mapping = aes(x = " ", y = n, fill = organic_option_available, group = organic_option_available))+
    geom_bar(stat = "identity", position = "stack", alpha = .7)+
    coord_flip()+
    labs(x = "", y = "Number of surveys", fill = "Organic option available")+
    theme_minimal()+
    scale_fill_manual(values = colors_org[c(2,4)])+
    geom_text(aes(label = paste0(perc, "%")), position = position_stack(vjust=0.5))+
    theme(legend.position = "top", plot.background = element_rect(fill = "white",colour = "white"))

ggsave(filename = "output/appendix_fig/Fig-S11_organic_option_available.png",
       plot = orgplot, device = "png",
       width = 8, height = 2)

### Type of retailers

# Small mistake in "distribuidor" name from the survey
recom_resellers$reseller_type[recom_resellers$reseller_type == "distribuido"] <- "distribuidor"

retail_type <- recom_resellers %>% count(country, reseller_type) %>% 
    mutate(country = str_to_title(country))

plot_retailers <- ggplot(data = retail_type, mapping = aes(x = reseller_type, y = n))+
    geom_bar(stat = "identity", fill = "grey50")+
    coord_flip()+
    facet_wrap(~country)+
    labs(y = "Number of surveys", x = "Type of retailer")

ggsave(plot = plot_retailers, filename = "output/appendix_fig/retailer_types_per_country.pdf", width = 6, height = 3, device = "pdf")
