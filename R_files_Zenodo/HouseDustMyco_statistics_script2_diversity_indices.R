
#########################Filename: HouseDustMyco_statistics_script2_diversity_indices.R

#Publication: Analyzing indoor mycobiomes through a large-scale citizen science study in Norway
#Authors: Pedro M. Martin-Sanchez, Eva Lena F. Estensmo, Luis N. Morgado, Sundy Maurice, Ingeborg B. Engh, Inger Skrede and Håvard Kauserud

##Related to the sections: 
        #Materials and Methods- 2.5. Statistical analyses
        #Results- 3.1. Data features and overall fungal diversity
        #Supplemental Information

## Goals: 
#(1) To check the variation of diversity indices by house compartments and geographical regions
#(2) To construct a GLM model to check the effect of variables on the species richness at sample level   

#Inputs: 
        #The rarefied fungal OTU table ("OTU_table_rar2000x10.xlsx")
        #The metadata file for dust samples ("metadata_file_dust_samples.xlsx")
#Outputs: 
        #Plots for the Figure 2 (Box plots showing Richness, Shannon index and betadisper index by house compartments)
        #Plots for the Figure S6 (Box plots showing Evenness, Inverse Simpson and Chao1 by house compartments)
        #Plots for the Figure S7 (Box plots showing Richness, Evenness and Shannon index in outdoor dust samples by regions in Norway)


#########################
#PREPARATIONS
#########################

#Load packages

library(openxlsx)
library(tidyverse)
library(vegan)
library(ggpubr)
library(rstatix) 

#Set your path
setwd("C:/WP1_data_analyses/script_preparation/")

#Load the rarefied OTU table (6,632 OTUs from 807 dust samples)
OTU_tab_rar=read.xlsx("OTU_table_rar2000x10.xlsx")

#Different format after moving the column "sample" to the row names
OTU_tab_rar_for_analyses=OTU_tab_rar %>% remove_rownames() %>% column_to_rownames("sample")

#Hellinger transformed OTU table
otutable_hel = OTU_tab_rar_for_analyses %>% decostand("hellinger")

#Load the metadata
metadata=read.xlsx("metadata_file_dust_samples.xlsx")
#Remark: Some geographic variables(latitude, longitude, municipality, county and region) as well as the construction year of the study houses
# were excluded from this file to protect personal data of volunteers 


#########################
#CALCULATION OF DIVERSITY INDICES
#########################

#Alpha diversity indices (Species Richness, Shannon index and Inverse Simpson index)
richness=OTU_tab_rar %>% remove_rownames() %>% column_to_rownames("sample")%>% specnumber() %>% as.data.frame() %>% rename(richness_div=".")
shannon=OTU_tab_rar %>% remove_rownames() %>% column_to_rownames("sample")%>% diversity() %>% as.data.frame() %>% rename(shannon_div=".")
invSimpson=OTU_tab_rar %>% remove_rownames() %>% column_to_rownames("sample")%>% diversity("inv") %>% as.data.frame() %>% rename(invsimpson_div=".")

#Beta diversity per house compartments ("location") using the betadisper function
beta_diversity_per_location<-betadisper(dist(otutable_hel), metadata$location, sqrt.dist=T)

#Extract the betadisper distances
betadisp_dist_per_sample<-beta_diversity_per_location$distances %>% as.data.frame()%>% rename(betadisp_div=".")

#Combine the calculated diversity indices as a data frame. Calculation of Evenness.
div_index=cbind(richness, shannon, invSimpson, betadisp_dist_per_sample) %>%rownames_to_column("sample") %>%  mutate(evenness_div=shannon_div/log(richness_div ))

#Merge them to the metadata 
metadata_enhanced=metadata  %>% left_join(div_index) 


#########################
#VARIATION OF DIVERSITY INDICES by INDOOR VS. OUTDOOR
#########################

#Data preparation
div.indoor_vs_outdoor=metadata_enhanced %>%
        select(contains("div"), indoor_outdoor) %>% 
        gather(div.Ind, div.ind.val,-indoor_outdoor)

#Box plots including statistics (ANOVA)
div.indoor_vs_outdoor %>% 
        ggplot(aes(x=indoor_outdoor,y=div.ind.val,fill=indoor_outdoor))+
        geom_boxplot()+
        facet_wrap(~div.Ind, scales = "free_y")+
        stat_compare_means(method = "anova")
#Remark: There is a consistent trend with higher alpha diversity indoors,
# while the betadisper showed the opposite trend 


#########################
#VARIATION OF DIVERSITY INDICES by HOUSE COMPARTMENTS ("location")
#########################

#Data preparation
div.location=metadata_enhanced %>%
        select(contains("div"), location) %>% 
        gather(div.Ind, div.ind.val,-location)

#Reorder the categories in the factor "location"  
div.location$location <- factor(div.location$location,
                                         levels = c('outside','central', 'bathroom'),ordered = TRUE)

#Box plots including statistics (ANOVA - global comparison of the means, parametric test)
div.location %>% 
        ggplot(aes(x=location,y=div.ind.val,fill=location))+
        geom_boxplot()+
        facet_wrap(~div.Ind, scales = "free_y")+
        stat_compare_means(method = "anova")# by default: Kruskal-Wallis (non-parametric) 
####Remark: These plots were used in the Figure 2 and the Supplemental Figure S6

#To add pairwise comparisons of the means -> t-test (parametric) or Wilcoxon test (non-parametric)
#Specify the comparisons we want
location_comparisons <- list( c("outside", "central"), c("outside", "bathroom"), c("central", "bathroom") )

#Box plots including t-test p-values for pairwise comparisons
div.location %>%
        ggplot(aes(x=location,y=div.ind.val,fill=location))+
        geom_boxplot()+
        facet_wrap(~div.Ind, scales = "free_y")+
        stat_compare_means(comparisons = location_comparisons, method = "t.test")#change by wilcox.test (default)

#Example for checking these results separately (e.g. for Richeness)
div.location_richness = div.location %>% filter (div.Ind== "richness_div")

compare_means(div.ind.val ~ location, method= "t.test", data = div.location_richness)
# A tibble: 3 x 8
#.y.         group1  group2             p      p.adj p.format p.signif method
#<chr>       <chr>   <chr>          <dbl>      <dbl> <chr>    <chr>    <chr> 
#1 div.ind.val outside central  0.000000639 0.0000013  6.4e-07  ****     T-test
#2 div.ind.val outside bathroom 0.000000103 0.00000031 1.0e-07  ****     T-test
#3 div.ind.val central bathroom 0.452       0.45       0.45     ns       T-test 

## Tukey HSD tests for the Figures 2 and S6

# Tukey HSD test for different house compartments (location) - Richness data
tukey_location_richness <- aov(div.ind.val ~ location, data = div.location_richness) %>%
        tukey_hsd()
tukey_location_richness
# A tibble: 3 x 9
#term     group1  group2   null.value estimate conf.low conf.high       p.adj p.adj.signif
#* <chr>    <chr>   <chr>         <dbl>    <dbl>    <dbl>     <dbl>       <dbl> <chr>       
#1 location outside central           0    27.1     14.1       40.0 0.00000315  ****        
#2 location outside bathroom          0    31.1     18.2       44.0 0.000000065 ****        
#3 location central bathroom          0     4.02    -8.84      16.9 0.743       ns  

# Tukey HSD test for different house compartments (location) - Shannon data
div.location_shannon = div.location %>% filter (div.Ind== "shannon_div")

tukey_location_shannon <- aov(div.ind.val ~ location, data = div.location_shannon) %>%
        tukey_hsd()
tukey_location_shannon
# A tibble: 3 x 9
#term     group1  group2   null.value estimate conf.low conf.high    p.adj p.adj.signif
#* <chr>    <chr>   <chr>         <dbl>    <dbl>    <dbl>     <dbl>    <dbl> <chr>       
#1 location outside central           0    0.478  0.317       0.640 1.92e-11 ****        
#2 location outside bathroom          0    0.629  0.468       0.790 0.       ****        
#3 location central bathroom          0    0.150 -0.00986     0.311 7.11e- 2 ns

# Tukey HSD test for different house compartments (location) -  Betadisp data
div.location_betadisp = div.location %>% filter (div.Ind== "betadisp_div")

tukey_location_betadisp <- aov(div.ind.val ~ location, data = div.location_betadisp) %>%
        tukey_hsd()
tukey_location_betadisp
# A tibble: 3 x 9
#term     group1  group2   null.value estimate conf.low conf.high p.adj p.adj.signif
#* <chr>    <chr>   <chr>         <dbl>    <dbl>    <dbl>     <dbl> <dbl> <chr>       
#1 location outside central           0 -0.0647  -0.0747    -0.0547 0     ****        
#2 location outside bathroom          0 -0.0571  -0.0671    -0.0472 0     ****        
#3 location central bathroom          0  0.00757 -0.00238    0.0175 0.175 ns  

# Tukey HSD test for different house compartments (location) -  Evenness data
div.location_evenness = div.location %>% filter (div.Ind== "evenness_div")

tukey_location_evenness <- aov(div.ind.val ~ location, data = div.location_evenness) %>%
        tukey_hsd()
tukey_location_evenness
# A tibble: 3 x 9
#term     group1  group2   null.value estimate conf.low conf.high    p.adj p.adj.signif
#* <chr>    <chr>   <chr>         <dbl>    <dbl>    <dbl>     <dbl>    <dbl> <chr>       
#1 location outside central           0   0.0737  0.0472     0.100  3.65e-10 ****        
#2 location outside bathroom          0   0.105   0.0785     0.132  0.       ****        
#3 location central bathroom          0   0.0312  0.00482    0.0577 1.55e- 2 *

# Tukey HSD test for different house compartments (location) -  Inverse Simpson data
div.location_invSimpson = div.location %>% filter (div.Ind== "invsimpson_div")

tukey_location_invSimpson <- aov(div.ind.val ~ location, data = div.location_invSimpson) %>%
        tukey_hsd()
tukey_location_invSimpson
# A tibble: 3 x 9
#term     group1  group2   null.value estimate conf.low conf.high    p.adj p.adj.signif
#* <chr>    <chr>   <chr>         <dbl>    <dbl>    <dbl>     <dbl>    <dbl> <chr>       
#1 location outside central           0     2.61    1.20       4.01 4.72e- 5 ****        
#2 location outside bathroom          0     4.66    3.25       6.07 2.03e-14 ****        
#3 location central bathroom          0     2.06    0.655      3.46 1.74e- 3 **


#########################
#VARIATION OF DIVERSITY INDICES in OUTDOOR SAMPLES by GEOGRAPHICAL REGIONS
#########################

#Remark: NOTE THAT THE VARIABLE "REGION" HAS BEEN EXCLUDED OF THE METADATA FILE
# In order to protect the personal data of the volunteers
# Therefore, THIS ANALYSIS CANNOT BE RUN!!
  

#Data preparation
div.region=metadata_enhanced %>%
        filter(indoor_outdoor== "outdoor")%>% #for outdoor samples
        select(richness_div,evenness_div,shannon_div, region) %>% 
        gather(div.Ind, div.ind.val,-region)

#reorder the categories in the factor "region" 
div.region$region <- factor(div.region$region,
                            levels = c('West','East','South', 'Mid', 'North','Svalbard'), ordered = TRUE)

#reorder the categories in the variables div.Ind.val for ordering the panels in the plot 
div.region$div.Ind <- factor(div.region$div.Ind,
                             levels = c('richness_div','evenness_div', 'shannon_div'),ordered = TRUE)

#Set a color blind friendly palette:
cbf_1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Box plots including statistics comparing mean values: ANOVA (global comparison)
div.region %>% 
        ggplot(aes(x=region,y=div.ind.val,fill=region))+
        geom_boxplot()+
        facet_wrap(~div.Ind, scales = "free_y")+
        stat_compare_means(method = "anova")+#by default Kruskal-Wallis
        scale_fill_manual(values =cbf_1)+
        labs(title= "Alpha diversity_OUTDOOR", subtitle = "Region")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
        theme_bw()+
        theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank())
####Remark: These plots were used in the Supplemental Figure S7


#To add pairwise comparisons of the means -> t-test (parametric) or Wilcoxon test (non-parametric)
#Specify the comparisons we want
regions_comparisons <- list( c("West", "East"), c("West", "South"), c("West", "Mid"), c("West", "North"), c("West", "Svalbard"),
                             c("East", "South"),c("East", "Mid"),c("East", "North"),c("East", "Svalbard"),
                             c("South", "Mid"),c("South", "North"),c("South", "Svalbard"),
                             c("Mid", "North"),c("Mid", "Svalbard"),
                             c("North", "Svalbard"))

#Box plots including t-test p-values for pairwise comparisons
div.region %>% 
        ggplot(aes(x=region,y=div.ind.val,fill=region))+
        geom_boxplot()+
        facet_wrap(~div.Ind, scales = "free_y")+
        stat_compare_means(comparisons = regions_comparisons, method = "t.test")+#wilcox.test (default)
        scale_fill_manual(values =cbf_1)+
        labs(title= "Alpha diversity_OUTDOOR", subtitle = "Region")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
        theme_bw()+
        theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Example for checking these results separately (e.g. for Richeness)
div.region_richness = div.region %>% filter (div.Ind== "richness_div")

compare_means(div.ind.val ~ region, method= "t.test", data = div.region_richness)
# A tibble: 15 x 8
#.y.         group1 group2        p p.adj p.format p.signif method
# <chr>       <chr>  <chr>     <dbl> <dbl> <chr>    <chr>    <chr> 
# 1 div.ind.val West   East     0.453      1 0.453    ns       T-test
# 2 div.ind.val West   South    0.435      1 0.435    ns       T-test
# 3 div.ind.val West   Mid      0.333      1 0.333    ns       T-test
# 4 div.ind.val West   North    0.435      1 0.435    ns       T-test
# 5 div.ind.val West   Svalbard 0.187      1 0.187    ns       T-test
# 6 div.ind.val East   South    0.720      1 0.720    ns       T-test
# 7 div.ind.val East   Mid      0.0967     1 0.097    ns       T-test
# 8 div.ind.val East   North    0.130      1 0.130    ns       T-test
# 9 div.ind.val East   Svalbard 0.240      1 0.240    ns       T-test
#10 div.ind.val South  Mid      0.141      1 0.141    ns       T-test
#11 div.ind.val South  North    0.182      1 0.182    ns       T-test
#12 div.ind.val South  Svalbard 0.219      1 0.219    ns       T-test
#13 div.ind.val Mid    North    0.824      1 0.824    ns       T-test
#14 div.ind.val Mid    Svalbard 0.123      1 0.123    ns       T-test
#15 div.ind.val North  Svalbard 0.141      1 0.141    ns       T-test 

# Tukey HSD test for different regions (Figure S7) - Richness data
tukey_region_richness <- aov(div.ind.val ~ region, data = div.region_richness) %>%
        tukey_hsd()
tukey_region_richness
# A tibble: 15 x 9
#term   group1 group2   null.value estimate conf.low conf.high p.adj p.adj.signif
#* <chr>  <chr>  <chr>         <dbl>    <dbl>    <dbl>     <dbl> <dbl> <chr>       
#1 region West   East              0    -8.42    -40.6      23.7 0.975 ns          
#2 region West   South             0   -14.2     -84.0      55.6 0.992 ns          
#3 region West   Mid               0    16.3     -29.2      61.8 0.908 ns          
#4 region West   North             0    12.3     -42.1      66.6 0.987 ns          
#5 region West   Svalbard          0   -69.5    -208.       68.6 0.699 ns          
#6 region East   South             0    -5.78    -71.2      59.6 1     ns          
#7 region East   Mid               0    24.7     -13.8      63.2 0.439 ns          
#8 region East   North             0    20.7     -28.0      69.4 0.827 ns          
#9 region East   Svalbard          0   -61.1    -197.       74.8 0.79  ns          
#10 region South  Mid               0    30.5     -42.4     103.  0.836 ns          
#11 region South  North             0    26.5     -52.3     105.  0.929 ns          
#12 region South  Svalbard          0   -55.3    -205.       94.1 0.896 ns          
#13 region Mid    North             0    -4.04    -62.4      54.3 1     ns          
#14 region Mid    Svalbard          0   -85.8    -225.       53.9 0.491 ns          
#15 region North  Svalbard          0   -81.7    -225.       61.1 0.571 ns
##Remark: All differences tested were not significant

# Tukey HSD test for different regions - Evenness data
div.region_evenness = div.region %>% filter (div.Ind== "evenness_div")

tukey_region_evenness <- aov(div.ind.val ~ region, data = div.region_evenness) %>%
        tukey_hsd()
tukey_region_evenness
#Remark: the only significant difference was for East-North p.adj=0.00217 **

# Tukey HSD test for different regions - Shannon data
div.region_shannon = div.region %>% filter (div.Ind== "shannon_div")

tukey_region_shannon <- aov(div.ind.val ~ region, data = div.region_shannon) %>%
        tukey_hsd()
tukey_region_shannon
#Remark: the only significant difference was for East-North p.adj=0.00751 **


#########################
# GLM model - Richness
#########################
##To check the effect of the selected variables on the Species Richness at sample level

# Selection of variables. We included the sequencing depth (seq_depth) as offset because this variable is correlated to richness
# Seq-depth was calculated on the quality-filtered fungal OTU table before rarefaction, and later included in the metadata file
data_for_glm_analyses<-metadata_enhanced %>%
        select(indoor_outdoor, location,# region, county, municipality,latitude, longitude,
               BIO1, BIO4, BIO9, BIO10, BIO11, BIO12, growing_season_length, swe_4, sca_2,
               ar50,geonorge_bedrock_nutrient, solar_radiation, area,
               building_type, building_material, ventilation_type, water_damage, moisture_problem, odor_problem, pest_type,
               people, females, children, asthma, allergy_type, pet_type, dust_coverage, seq_depth)

#Remark: NOTE THAT THE VARIABLES "LATITUDE", "LONGITUDE", "MUNICIPALITY", COUNTY" AND "REGION" HAVE BEEN EXCLUDED OF THE METADATA FILE
# In order to protect the personal data (addresses) of the volunteers
# Therefore, this example shows slightly different results compared to those detailed in the article considering al these variables  

# Preparation of Richness data
richness_per_sample<-metadata_enhanced %>% select(richness_div)

# Construct the null model used as reference in the first screening (below) 
richness_null_model=glm(richness_per_sample$richness_div~1, offset = seq_depth,
                        data =data_for_glm_analyses ) 
AIC(richness_null_model) #### AIC  18974.12
BIC(richness_null_model) #### BIC  18983.51

# First screening of selected variable in a loop
mylist=colnames(data_for_glm_analyses)[1:30]
results=data.frame(matrix(nrow=length(mylist),ncol=4))
colnames(results)=c("Variable","p", "AIC", "BIC")

for (i in (1:length(mylist))) {
        sub=data.frame(data_for_glm_analyses[,colnames(data_for_glm_analyses)==mylist[i]], data_for_glm_analyses$seq_depth)
        colnames(sub)=c("variable", "seq_depth")
        model_data<-glm(richness_per_sample$richness_div~variable, offset = seq_depth, data=sub)
        print(mylist[i])
        print(model_data)
        plot(resid(model_data))
        results$Variable[i]=mylist[i]
        results$p[i]=summary(model_data)$coefficients[8]
        results$AIC[i] =AIC(model_data)
        results$BIC[i] =BIC(model_data)
}

Results_richness_GLM<-results %>% arrange(AIC) # Get and see the results
Results_richness_GLM

#Variable                                 p       AIC       BIC
#1                  seq_depth  0.000000e+00  8966.459  8980.539
#2                   location  3.156173e+00 18952.452 18971.225
#3             indoor_outdoor  7.481675e-05 18960.389 18974.469
#4                      BIO10  4.552681e-04 18963.794 18977.874
#5                       BIO4  3.015283e-03 18967.297 18981.377
#6                       BIO9  8.862550e-03 18969.251 18983.331
#7                      sca_2  9.499919e-03 18969.375 18983.455
#8                      BIO12  1.658832e-02 18970.365 18984.445
#9                   children  1.673343e-02 18970.381 18984.461
#10                    people  5.644844e-02 18972.472 18986.552
#11                  pet_type  2.996007e+03 18972.884 18996.351
#12              allergy_type  2.006499e+04 18972.918 19015.158
#13           solar_radiation  8.042550e-02 18973.057 18987.137
#14                    asthma  1.063331e-01 18973.506 18987.586
#15          ventilation_type  3.842284e+03 18973.601 18997.068
#16              water_damage  1.177635e-01 18973.668 18987.748
#17          moisture_problem  1.190588e-01 18973.685 18987.765
#18     growing_season_length  2.026914e-01 18974.495 18988.574
#19                     BIO11  2.290633e-01 18974.671 18988.751
#20                   females  2.384661e-01 18974.728 18988.808
#21                     swe_4  2.725896e-01 18974.915 18988.995
#22 geonorge_bedrock_nutrient -1.553496e+00 18975.129 18993.903
#23                      BIO1  3.219805e-01 18975.139 18989.218
#24              odor_problem  3.648747e-01 18975.299 18989.379
#25                      area  7.241008e-01 18975.997 18990.077
#26             dust_coverage  9.109700e-01 18976.110 18990.190
#27         building_material  1.045075e+00 18976.470 18995.243
#28             building_type -4.016841e-01 18977.745 18996.518
#29                      ar50  8.710586e+03 18979.091 19011.944
#30                 pest_type  6.278603e+03 18981.053 19013.907


# The variables wih AIC < 18974.12 (value for null model) and p-value<0.05 are shown below
#(those variables that can improve the explanatory power of the GLM model) 

#Variable                                 p       AIC       BIC
#3             indoor_outdoor  7.481675e-05 18960.389 18974.469
#4                      BIO10  4.552681e-04 18963.794 18977.874
#5                       BIO4  3.015283e-03 18967.297 18981.377
#6                       BIO9  8.862550e-03 18969.251 18983.331
#7                      sca_2  9.499919e-03 18969.375 18983.455
#8                      BIO12  1.658832e-02 18970.365 18984.445
#9                   children  1.673343e-02 18970.381 18984.461

#Fix the best explanatory variable (indoor_outdoor) in the model 
richness_indoor_outdoor_model=glm(richness_per_sample$richness_div~indoor_outdoor,offset = seq_depth, data=data_for_glm_analyses)
summary(richness_indoor_outdoor_model) 

#Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
#(Intercept)             -30883       1312 -23.538  < 2e-16 ***
#indoor_outdooroutdoor     9098       2285   3.981 7.48e-05 ***
#        ---
#        Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#(Dispersion parameter for gaussian family taken to be 931284059)
#Null deviance: 7.6444e+11  on 806  degrees of freedom
#Residual deviance: 7.4968e+11  on 805  degrees of freedom
#AIC: 18960
#Number of Fisher Scoring iterations: 2

AIC(richness_indoor_outdoor_model)### AIC 18960.39
BIC(richness_indoor_outdoor_model)### BIC 18974.47
plot(resid(richness_indoor_outdoor_model))

## Second loop trying to improve our model
for (i in (1:length(mylist))) {
        sub=data.frame(data_for_glm_analyses[,colnames(data_for_glm_analyses)==mylist[i]], data_for_glm_analyses$seq_depth, data_for_glm_analyses$indoor_outdoor)
        colnames(sub)=c("variable","seq_depth","indoor_outdoor")
        model_data<-glm(richness_per_sample$richness_div~variable+indoor_outdoor, offset = seq_depth, data=sub)
        print(mylist[i])
        print(model_data)
        plot(resid(model_data))
        results$Variable[i]=mylist[i]
        results$p[i]=summary(model_data)$coefficients[8]
        results$AIC[i] =AIC(model_data)
        results$BIC[i] =BIC(model_data)
}

Results_richness_GLM_2nd_round<-results %>% arrange(BIC)
Results_richness_GLM_2nd_round

#Variable             p       AIC       BIC
#1                  seq_depth -1.436811e+04  8910.447  8929.221
#2                      BIO10  3.561774e+00 18949.755 18968.528
#3                   location  3.156173e+00 18952.452 18971.225
#4                       BIO4  2.988382e+00 18953.474 18972.248
#5                       BIO9 -2.633072e+00 18955.460 18974.233
#6                      sca_2  2.596555e+00 18955.650 18974.423
#7             indoor_outdoor  7.481675e-05 18960.389 18974.469
#8                      BIO12 -2.403735e+00 18956.610 18975.383
#9                   children  2.396327e+00 18956.645 18975.419
#10                    people  1.910348e+00 18958.734 18977.507
#11           solar_radiation  1.792264e+00 18959.171 18977.944
#12                    asthma -1.632111e+00 18959.719 18978.493
#13              water_damage  1.601969e+00 18959.817 18978.590
#14          moisture_problem  1.584247e+00 18959.873 18978.647
#15     growing_season_length -1.257819e+00 18960.802 18979.576
#16                     BIO11 -1.191774e+00 18960.964 18979.738
#17                   females  1.171911e+00 18961.011 18979.785
#18                     swe_4 -1.117885e+00 18961.135 18979.909
#19                      BIO1  1.022112e+00 18961.341 18980.114
#20              odor_problem  9.384121e-01 18961.505 18980.279
#21                      area  3.781096e-01 18962.245 18981.019
#22             dust_coverage -1.286508e-01 18962.372 18981.145
#23 geonorge_bedrock_nutrient  2.283877e+03 18961.317 18984.784
#24         building_material  2.285737e+03 18962.673 18986.140
#25                  pet_type  6.428510e+03 18959.063 18987.223
#26             building_type  2.287640e+03 18964.028 18987.495
#27          ventilation_type  3.351891e+03 18959.785 18987.945
#28                      ar50  3.604913e+03 18965.236 19002.782
#29                 pest_type  5.630095e+03 18967.269 19004.816
#30              allergy_type  2.006499e+04 18958.916 19005.849

#Remark: There is no further variable that improves the model (AIC<18960.39) significantly (p<0.05) 


#################### END OF THE SCRIPT ##########################