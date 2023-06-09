#!/usr/bin/env RScript

library("tidyverse")
library("readxl")
library("xlsx")

ingr_raw <- read_xls("data/04_pesticide_list_ingredients/04_ingr_homogenized_enriched/homogenized-active-ingredients-all-countries.xls")


## Join prouct id with number of times recommended
all_recom <- read_csv2("data/03_product_names_cleaned/20_recom_with_id/all_countries_recom_commercial_names_ids.csv",
                       col_types = cols(.default = "c")) 


all_recom_nb <- all_recom %>% group_by(product_id) %>% summarize(nb_recom = n())

ingr <- ingr_raw %>% left_join(all_recom_nb, by = "product_id") %>% filter(type == "insecticide" | type == "fungicide") %>%
    filter(!is.na(nb_recom))

## Rename the factors that are identical across pesticide types to nest them 
ingr %>% group_by(type, pesticide_family) %>% summarise(n = n()) %>% ungroup  %>% filter(duplicated(pesticide_family))

ingr$pesticide_family[ingr$pesticide_family == "other" & ingr$type == "fungicide"] <- "other_fungicide"
ingr$pesticide_family[ingr$pesticide_family == "other" & ingr$type == "insecticide"] <- "other_insecticide"
ingr$pesticide_family[ingr$pesticide_family == "dicarboximide" & ingr$type == "insecticide"] <- "dicarboximide_insect"
ingr$pesticide_family[ingr$pesticide_family == "dicarboximide" & ingr$type == "fungicide"] <- "dicarboximide_fung"
ingr$pesticide_family[ingr$pesticide_family == "carbamate" & ingr$type == "fungicide"] <- "carbamate_fung"
ingr$pesticide_family[ingr$pesticide_family == "carbamate" & ingr$type == "insecticide"] <- "carbamate_insect"
ingr$pesticide_family[ingr$pesticide_family == "antibiotic" & ingr$type == "fungicide"] <- "antibiotic_fung"

## Same for duplicated names across active ingredients
ingr %>% group_by(pesticide_family, active_ingr_only_reconc) %>% summarise(n = n()) %>% ungroup  %>% filter(duplicated(active_ingr_only_reconc))

ingr$active_ingr_only_reconc[ingr$pesticide_family == "antibiotic_fung" & ingr$active_ingr_only_reconc == "kasugamycin"] <- "kasugamycin_fung"
ingr$active_ingr_only_reconc[ingr$pesticide_family == "dicarboximide_insect" & ingr$active_ingr_only_reconc == "captan"] <- "captan_ins"

# Same for duplicates among commercial names

duplicated_names <-  ingr %>% group_by(country, product_id, ingr_uniq_id) %>% summarise(nb_ingr = n()) %>% 
    ungroup %>% filter(duplicated(product_id))

for(i in 1:nrow(duplicated_names)){
    rownb <- which(ingr$ingr_uniq_id == duplicated_names$ingr_uniq_id[i])
    ingr$product_id[rownb] <- paste0(duplicated_names$product_id[i], "-", duplicated_names$country[i])
}

# Last exception is copper which is both an active ingredient and a family
ingr$active_ingr_only_reconc[ingr$active_ingr_only_reconc == "copper"] <- "copper_ing"

# Get rid of products without pesticide family
ingr_all <- ingr %>% filter(!is.na(pesticide_family))

# Keep first ingredient only if products has two or more ingredients; this prevent to overrepresent products
ingr <- ingr_all %>% filter(ingr_id == 1)

## Uppercase first letter of future labels
ingr <- ingr %>% mutate(active_ingr_only_reconc = paste(toupper(substring(active_ingr_only_reconc, 1, 1)), substring(active_ingr_only_reconc, 2), sep = ""))
ingr <- ingr %>% mutate(pesticide_family = paste(toupper(substring(pesticide_family, 1, 1)), substring(pesticide_family, 2), sep = ""))

## Export unique families to match mode of action
pesticide_families <- ingr %>% distinct(type, pesticide_family, active_ingr_only_reconc) %>% 
    arrange(type, pesticide_family, active_ingr_only_reconc)
# write_excel_csv2(pesticide_families, "data/pesticide_families.csv")

## Add mode of action
moa <- read_xlsx("data/90_mode_of_action/pesticides_mode_of_action_filled.xlsx") 

ingr_moa <- ingr %>% inner_join(moa, by = c("type", "pesticide_family", "active_ingr_only_reconc"))


ingr_moa %>% group_by(type) %>% distinct(moa_group_pesticide) %>% summarise(nb = n())


### Check Rotterdam annex 3 for banned products
library(janitor)

eu_list <- read_xlsx("data/80_pop_list_stockhlom/EU_not_approved_chemicals.xlsx")%>% clean_names()

eu_list_no_cas <- eu_list %>% filter(is.na(cas_registry_number))
eu_list_cas <- eu_list %>% filter(!is.na(cas_registry_number)) %>% 
    select(cas_registry_number, status_under_reg_ec_no_1107_2009)

# eu_list_cas %>% inner_join(ingr_moa, by = c("cas_registry_number" = "ingr_CAS"))
ingr_moa <- ingr_moa %>% left_join(eu_list_cas, by = c("ingr_CAS" = "cas_registry_number")) %>% 
    mutate(legal_status = ifelse(is.na(status_under_reg_ec_no_1107_2009), "ok", status_under_reg_ec_no_1107_2009))

## Some chemicals have several CAS numbers because they have several isomers --> Manual check hereafter
ingr_moa[ingr_moa$active_ingr_only_reconc == "Λ-cyhalothrin", ]$legal_status <- "Not Approved"
ingr_moa[ingr_moa$active_ingr_only_reconc == "Methomyl", ]$legal_status <- "Not Approved"
ingr_moa[ingr_moa$active_ingr_only_reconc == "Beta−cyfluthrin", ]$legal_status <- "Not Approved"


## RESULTS NUMBERS ###
ingr_moa %>% group_by(type) %>% summarise(n = n(), nb_recom_sum = sum(nb_recom)) 
ingr_moa %>% group_by(moa_group_pesticide) %>% summarise(n = n(), nb_recom_sum = sum(nb_recom)) 

ingr_moa %>% group_by(type, moa_group_pesticide) %>% count %>%  group_by(type) %>% count 
ingr_moa %>% group_by(country, moa_group_pesticide, type) %>% count %>% group_by(country, type) %>% count

ingr_moa %>% group_by(type,moa_group_pesticide) %>% summarise(n = n(), nb_recom_sum = sum(nb_recom)) %>% 
    arrange(desc(nb_recom_sum))

ingr_moa %>% group_by(type,active_ingr_only_reconc) %>% summarise(n = n(), nb_recom_sum = sum(nb_recom)) %>% 
    arrange(desc(nb_recom_sum))

## MoA and legal status of products - product number/percent

moa_hazard_nb <- ingr_moa %>% group_by(legal_status, moa_group_pesticide) %>% 
    summarise(nb_prod = n()) %>%
    group_by(moa_group_pesticide) %>% 
    mutate( nb_prod_moa = sum(nb_prod),
            prod_perc = nb_prod/sum(nb_prod)*100) %>%
    group_by(moa_group_pesticide) %>%
    mutate(order = sum(nb_prod)) %>%
    arrange(moa_group_pesticide, legal_status)

moa_hazard_nb$legal_status[moa_hazard_nb$legal_status == "ok"] <- "Approved"
moa_hazard_nb$legal_status[moa_hazard_nb$legal_status == "Not Approved"] <- "Not approved"

colors_org <- c( "black", "#d71e36", "#ff7f27", "#015ba9", "#00a261")

plot_approved_moa <- ggplot(data = moa_hazard_nb, mapping = aes(x = reorder(moa_group_pesticide, order),
                                                                y = nb_prod,
                                                                fill = legal_status))+
    geom_bar(stat = "identity", position = "stack")+
    geom_text(aes(label = paste0(round(prod_perc,1), "%")),
              position = position_stack(vjust=0.5),
              size = 2.5)+
    labs(x = "Pesticide mode of action", y = "Number of products", fill = "Approved in the EU")+
    coord_flip()+
    scale_fill_manual(values = c("#00a261","#ff7f27" ))+
    theme_minimal()+
    theme(legend.position = "top", plot.background = element_rect(fill = "white",colour = "white"))

ggsave(filename = "Fig-S13_MoA-approved_EU_nbprod.png",
       plot = plot_approved_moa,
       device = "png",
       width = 10,
       height = 6)

## MoA and type - recommendation number/percent

moa_hazard_nb_recom <- ingr_moa %>% group_by(type, moa_group_pesticide) %>% 
    summarise(nb_recom_sum = sum(nb_recom)) %>%
    group_by(type) %>% 
    mutate( nb_recom_moa = sum(nb_recom_sum),
            recom_perc = nb_recom_sum/sum(nb_recom_sum)*100,
            type = str_to_title(type))

plot_approved_moa_recom <- ggplot(
    data = moa_hazard_nb_recom, mapping = aes(x = reorder(moa_group_pesticide, nb_recom_sum),
                                                                y = nb_recom_sum))+
    geom_bar(stat = "identity", position = "stack", fill = "grey65")+
    geom_text(aes(label = paste0(round(recom_perc,1), "%")),
              position = position_stack(vjust=0.5),
              size = 2.5)+
    labs(x = "Pesticide mode of action", y = "Number of recommendations", fill = "Approved in the EU")+
    coord_flip()+
    facet_wrap(~type, scales = "free_y")+
    theme_minimal()+
    theme(legend.position = "top", plot.background = element_rect(fill = "white",colour = "white"))

ggsave(filename = "Fig-S14_MoA-type_nbrecom.png",
       plot = plot_approved_moa_recom,
       device = "png",
       width = 8,
       height = 5)


ingr_moa %>% group_by(legal_status) %>% summarise(nb_prod = n(), nb_recom = sum(nb_recom)) %>% 
    mutate(nb_prop = nb_prod/sum(nb_prod)*100,
           recom_prop = nb_recom/sum(nb_recom)*100)

table_recom_banned <- ingr_moa %>% group_by(legal_status, type) %>% summarise(nb_prod = n(), nb_recom = sum(nb_recom)) %>% group_by(type) %>% 
    mutate(nb_prop = nb_prod/sum(nb_prod)*100,
           recom_prop = nb_recom/sum(nb_recom)*100)

write.csv2(x = table_recom_banned, file = "number_proportion_products_recom_bannedEU.csv", row.names = F)



export_table_pesticide_active_ingr <- ingr_moa %>% group_by(type, active_ingr_only_reconc) %>% mutate(nb_recom_across_countries = sum(nb_recom)) %>% 
    group_by(country, legal_status, type, active_ingr_only_reconc, nb_recom_across_countries) %>% 
    summarise(nb_prod = n(), nb_recom_sum = sum(nb_recom), nb_recom_across = mean (nb_recom_across_countries)) %>%
    arrange(type, desc(nb_recom_across_countries)) 

plot_insecticides <- ggplot(data = filter(export_table_pesticide_active_ingr, type == "insecticide"), mapping = 
                                aes(x = reorder(active_ingr_only_reconc, (nb_recom_across_countries)), y = nb_recom_sum, fill = country))+
    geom_bar(stat = "identity",  position ="stack")+
    coord_flip()+
    xlab("Active substance")+
    ylab("Number of times recommended")+
    ylim(c(0,280))+
    ggtitle("Insecticides")+
    theme(text = element_text(size=10)) 



plot_fungicides <- ggplot(data = filter(export_table_pesticide_active_ingr, type == "fungicide"), mapping =
                              aes(x = reorder(active_ingr_only_reconc, (nb_recom_across_countries)), y = nb_recom_sum, fill = country))+
    geom_bar(stat = "identity",  position ="stack")+
    coord_flip()+
    xlab("Active substance")+
    ylab("Number of times recommended")+
    ylim(c(0,280))+
    ggtitle("Fungicides")+
    theme(text = element_text(size=10)) 



library(ggpubr)
plot_pesticides <- ggarrange(plotlist = list(plot_insecticides, plot_fungicides), common.legend = TRUE)


ggsave(filename = "all_products_active_ingr_unapproved_EU_stacked.pdf",plot = plot_pesticides, device = "pdf", width = 8, height = 9)


## Create edge list with useful values to plot
prod <- ingr_moa %>%group_by(type, moa_group_pesticide, active_ingr_only_reconc, product_id, toxic_cat, legal_status) %>%
    summarize(count = n(), nb_recom = sum(nb_recom)) %>% mutate(root = "root") 

edges_root <- prod %>% ungroup() %>% group_by(root, type) %>% 
    summarize(count = n(), nb_recom = sum(nb_recom)) %>% 
    rename(from = 1, to = 2)%>% 
    mutate(group = to, legal_status = "ok", children_group = '')

edges_type <- prod %>% ungroup() %>% group_by(type, moa_group_pesticide) %>% 
    summarize(count = n(), nb_recom = sum(nb_recom)) %>%
    rename(from = 1, to = 2) %>% 
    mutate(group = from, legal_status = "ok", children_group = to)

edges_family <- prod %>% ungroup() %>% group_by(moa_group_pesticide, active_ingr_only_reconc) %>% 
    summarize(count = n(), nb_recom = sum(nb_recom)) %>% 
    rename(from = 1, to = 2)%>% 
    left_join(edges_type, by = c("from" = "to")) %>% 
    select(-count.y, -from.y, -nb_recom.y) %>% 
    rename(count = count.x, nb_recom = nb_recom.x) %>% 
    mutate(legal_status = "ok", children_group = to)

edges_name <- prod %>% ungroup() %>% group_by(active_ingr_only_reconc, product_id, legal_status) %>% 
    summarize(count = n(), nb_recom = sum(nb_recom)) %>% 
    rename(from = 1, to = 2)%>% 
    left_join(edges_family, by = c("from" = "to")) %>% 
    select(-count.y, -from.y, -legal_status.y, -nb_recom.y) %>% 
    rename(count = count.x, legal_status = legal_status.x, nb_recom = nb_recom.x) %>% 
    mutate(children_group = '')

edge_list <- rbind(edges_root, edges_type, edges_family, edges_name) %>% filter(!is.na(from), !is.na(to)) %>% 
    mutate(children_group = ifelse(nb_recom < 10, '', children_group))  ## Show labels only for products that have been recommended more than 10 times


### Plot with tidygraph

library(tidygraph)
library(ggraph)
library(igraph)
library(RColorBrewer) 
library(gridExtra) 

## Create a graph in order to determine what are the leaves and roots
graph_node_detec <- graph_from_data_frame( edge_list, vertices = unique(c(edge_list$from, edge_list$to)))
graph_node_detec_tbl <- as_tbl_graph(graph_node_detec)

## Determine leaf, cut and root
mygraph <- graph_node_detec_tbl %>% mutate(leaf = node_is_leaf(), root = node_is_source(), middle = node_is_cut())

# #Customize the nodes
mygraph_angles <- mygraph %>%as_tibble() %>%  mutate(id = NA, angle = NA, hjust = NA)

mygraph_angles$group <- edge_list$from[match(mygraph_angles$name, edge_list$to)]

# Arrange the labels for leaves
leaves <- which(mygraph_angles$leaf == TRUE)
nleaves <- length(leaves)
mygraph_angles$id[leaves] <- seq(1:nleaves)
mygraph_angles$angle <- 360 * mygraph_angles$id / nleaves
mygraph_angles$hjust<-ifelse( mygraph_angles$angle < -90, 0, 1)
mygraph_angles$angle<-ifelse(mygraph_angles$angle < -90, mygraph_angles$angle+180, mygraph_angles$angle)
mygraph_angles$nb_recom <- edge_list$nb_recom[ match( mygraph_angles$name, edge_list$to ) ]

#Recreate graph
mygraph <- graph_from_data_frame(edge_list, vertices = mygraph_angles)

mygraph <- mygraph %>% as_tbl_graph %>%  activate(nodes)%>% arrange(group, desc(nb_recom))



# Make the plot
hyperbolic_tree <- ggraph(mygraph, layout = 'dendrogram', circular = F) + 
    geom_edge_diagonal(aes(width = nb_recom, col = legal_status)) +
    scale_edge_width_continuous(range = c(0.7,40))+
    geom_node_text(aes(filter = middle, x = x, y=y,  label=name,
                       angle = 0,
                       hjust=0),
                   size=6, alpha=1,  check_overlap = TRUE) +  scale_edge_colour_manual(values = c("#ea8e9a","grey85"))  +
    scale_size_continuous( range = c(1,15) ) +
    theme_void() + coord_flip()+ scale_y_reverse()+scale_x_reverse()+
    theme(
        plot.margin=unit(c(0,0,0,0),"cm"),
        
        legend.position = "none"
    ) +
    expand_limits(x = c(-1.03, 1.03), y = c(-1.03, 1.03))


##########################
## CHeck nb of multiple products per recom

all_recom <- read_csv2("data/03_product_names_cleaned/20_recom_with_id/all_countries_recom_commercial_names_ids.csv",
                       col_types = cols(.default = "c")) %>%
                       select(product_id, `_submission__uuid`)


all_recom_nb <- all_recom %>% group_by(product_id) %>% summarize(nb_recom = n())



products_recom_firstingr <- all_recom %>% left_join(ingr_raw, by = "product_id") %>% 
    filter(type == "insecticide" | type == "fungicide") %>% filter(ingr_id == 1) %>% 
    group_by(`_submission__uuid`) %>% mutate(nb_products_per_recom = n()) 

ingr_moa2 <- ingr_moa %>% select(product_id, moa_group_pesticide, moa_group_name)

products_recom_moa <- products_recom_firstingr %>% left_join(ingr_moa2, by = "product_id")


product_type <- products_recom_moa %>% ungroup %>% select(product_id, type) %>%
    distinct(product_id, .keep_all = TRUE)

product_moa <- products_recom_moa %>% ungroup %>% select(product_id, moa_group_pesticide) %>%
    distinct(product_id, .keep_all = TRUE)

recom_more_prod <- products_recom_moa %>%  filter(nb_products_per_recom > 1) %>% 
    ungroup %>% select(`_submission__uuid`, product_id) %>% 
    group_by(`_submission__uuid`) %>% mutate(product = paste0("prod_",row_number())) %>% mutate(nb = 1) %>% 
    pivot_wider( names_from = product, values_from = product_id)

prod_cocktail_inter <- recom_more_prod %>% group_by(prod_1, prod_2) %>% 
    summarise(width = n()) %>% 
    rename(from = prod_1, 
           to = prod_2) %>% 
    left_join(product_type, by = c("from" = "product_id")) %>% rename(type_from = type) %>% 
    left_join(product_type, by = c("to" = "product_id")) %>% rename(type_to = type) %>% 
    mutate(combined_type = ifelse(type_from == type_to, "0", "1")) %>% 
    left_join(product_moa, by = c("from" = "product_id")) %>% rename(moa_from = moa_group_pesticide) %>% 
    left_join(product_moa, by = c("to" = "product_id")) %>% rename(moa_to = moa_group_pesticide) %>% 
    mutate(combined_moa = ifelse(moa_from == moa_to, "same_moa", "diff_moa")) %>% 
    mutate(type_cocktail = ifelse(combined_type == 1, "diff_type", combined_moa)) %>% 
    na.omit()
    

graph_cocktail <- graph_from_data_frame(prod_cocktail_inter)
graph_cocktail <- graph_cocktail %>% as_tbl_graph %>%  activate(nodes) 


## Proportion of recom with 1, 2 and 3+ products
recom_nb_prod <- products_recom_moa %>% 
    group_by(`_submission__uuid`, nb_products_per_recom) %>% 
    summarize() %>% 
    mutate(nb_products_per_recom = 
               ifelse(nb_products_per_recom > 3, 3, nb_products_per_recom)) %>% 
    group_by(nb_products_per_recom) %>% 
    summarize(qty = n()) %>% 
    mutate(total = sum(qty), prop = qty/total*100)

recom_nb_prod$nb_products_per_recom[recom_nb_prod$nb_products_per_recom == 3] <- "3+"

plot_recom_nb_prod <- ggplot(data = recom_nb_prod, mapping = aes(x = nb_products_per_recom, y = qty))+
    geom_bar(stat = "identity", position = "stack", alpha = .7)+
    labs(y = "Number of surveys", x = "Number of products recommended")+
    theme_minimal()+
    geom_text(aes(label = paste0(round(prop,1), "%")), position = position_stack(vjust=0.5))+
    theme(legend.position = "top", plot.background = element_rect(fill = "white",colour = "white"))

ggsave(filename = "Fig-S12_number_of_prod_per_recom.png",
       plot = plot_recom_nb_prod, device = "png",
       width = 3, height = 4)


## Proportion of products in mixture
n_prod_mixture <- products_recom_moa %>% group_by(`_submission__uuid`, nb_products_per_recom, product_id) %>% summarize() %>% 
    filter(nb_products_per_recom > 1) %>% ungroup() %>% distinct(product_id) %>% nrow() 
n_prod_total <- products_recom_moa %>% ungroup %>% distinct(product_id) %>% nrow()
n_prod_mixture/n_prod_total

polar_barchart_data <- mygraph  %>% activate(nodes) %>% arrange(node_topo_order()) %>% 
    graph_join(graph_cocktail, by = c("name")) %>% mutate(name2 = fct_relevel(name, name)) %>%  filter(leaf == TRUE)

cockail_interactions <- ggraph(polar_barchart_data, layout = 'linear', circular = F) + 
    geom_edge_arc(aes(width = width, col = type_cocktail), fold =TRUE, alpha = 0.4) + 
    theme_void() +
    scale_edge_width_continuous(range = c(0.5,5))+
    scale_edge_colour_manual(values = c("#00a261","#015ba9" , "#ff7f27")) +
    # scale_edge_alpha_manual(values = c(0.3, 0.8)) +
    coord_flip()+ scale_y_reverse()+scale_x_reverse()+
    theme(
        plot.margin=unit(c(0,0,0,0),"cm"),
        legend.position = "none"
    ) + 
    expand_limits(x = c(-1.03, 1.03), y = c(-1.03, 1.03))


## Add a plot with proportion of cocktail scenarios for each pesticide type 

prop_cocktail_scenarios <- prod_cocktail_inter %>% group_by(type_cocktail) %>% summarise(n_recom = sum(width)) %>% 
    ungroup %>% mutate(tot = sum(n_recom), prop = n_recom/tot*100, x = "x") %>% arrange(desc(prop))


porp_moa <- ggplot(data = prop_cocktail_scenarios, mapping = aes(x = x, y = prop, fill = factor(type_cocktail, levels=c("diff_moa", "same_moa", "diff_type"))))+
    geom_bar(stat = "identity", position="stack", width = 0.8, alpha = 0.5)+
    geom_text(aes(label=paste0(n_recom, " \\n (", sprintf("%1.1f", prop),"%)")), 
              position=position_stack(vjust=0.5), colour="black") + 
    scale_fill_manual(values = c("#00a261", "#ff7f27", "#015ba9"))+
    theme_void()+
    theme(legend.position="bottom")

library(ggpubr)
## Combine plots
pdf(file = "output/FIG_3.pdf", width = 21, height = 25)
    ggarrange( hyperbolic_tree,cockail_interactions,porp_moa, nrow = 1, widths = c(6,1, 1))
dev.off()

## Export data from graph 
graph_data_export <- edge_list %>% 
    rename(node_from = from, node_to = to, nb_children = count, nb_recommendations = nb_recom, approved_in_EU = legal_status) %>% 
    select(-children_group)
proportion_cocktail_type_export <- prop_cocktail_scenarios %>% select(-x)
cocktail_data_arc_export <- prod_cocktail_inter %>% 
    rename(product_from = from, product_to = to, nb_recommendations = width) %>% 
    select(-combined_moa, -combined_type)

write.table(graph_data_export, "output/table/Fig_3_dataset_graph_tree.csv", sep =";", dec = ".", row.names = FALSE)
write.table(cocktail_data_arc_export, "output/table/Fig_3_dataset_graph_cocktail_arcs.csv", sep =";", dec = ".", row.names = FALSE)
write.table(proportion_cocktail_type_export, "output/table/Fig_3_dataset_cocktail_proportions.csv", sep =";", dec = ".", row.names = FALSE)
