### 00. visualize effect sizes from BEDASSLE runs ------------------------------

# import geographic distance between populations
setwd("~/Desktop/DRYAD")
geo_dist <- read.csv(file = "05_geographic_distance_matrix.csv", header = T)

# add row names to the data.frame, remove row names from the first column
rownames(geo_dist) <- geo_dist[,1]
geo_dist <- geo_dist[,-1]


# load BEDASSLE
# install.packages("/Users/Kyle/Documents/Post_Doc/Collaborations/Mike_Rotter/bedassle.tar.gz", repos = NULL, type = "source")
library(matrixcalc)
library(emdbook)
library(BEDASSLE)


# import BEDASSLE runs
setwd("~/Desktop/DRYAD/BEDASSLE_runs/")

### list all files in "bedassle" folder
files <- list.files(pattern = ".Robj")

# subset MCMC outputs
my_runs <- files[grepl("MCMC", files)]


### 01. build function to combine data from all 10 MCMC runs -------------------

Find.Effect.Sizes <- function(MCMC_run_name){
    
    # import MCMC output
    show(load(MCMC_run_name))
    
    # extract estimated effect sizes after 20% burn-in
    aE1_pbi <- as.vector(aE[1, 201:1000]) # herbivores
    aE2_pbi <- as.vector(aE[2, 201:1000]) # vegetation
    aE3_pbi <- as.vector(aE[3, 201:1000]) # climate (temp)
    aE4_pbi <- as.vector(aE[4, 201:1000]) # climate (precipitation)
    aE5_pbi <- as.vector(aE[5, 201:1000]) # climate (seasonality)
    
    ### combine data to plot parameter estimates
    my_df <- data.frame(
        estimate = c(aE1_pbi, aE4_pbi,aE5_pbi, aE2_pbi, aE3_pbi ),
        factor = rep(c("Herbivores", "Precipitation", "Seasonality", "Vegetation", "Temperature" ), each = 800))
    my_df$run <- MCMC_run_name
    
    return(my_df)
}


# test function
Find.Effect.Sizes(MCMC_run_name = my_runs[1])
Find.Effect.Sizes(MCMC_run_name = my_runs[2])

# compile aE:aD ratios across all runs
output <- lapply(my_runs, Find.Effect.Sizes)

# combine into single df
output2 <- do.call(rbind, output)

# convert "run" to factor
output2$run <- as.factor(output2$run)
levels(output2$run) <-  1:length(my_runs)

str(output2)


### 02. summarize effect sizes for all predictors ------------------------------

### herbivores
median(output2$estimate[output2$factor == "Herbivores"]) 
    # 0.117
quantile(output2$estimate[output2$factor == "Herbivores"], probs = c(0.25, 0.75)) 
    # 0.04796607 0.24472826 
quantile(output2$estimate[output2$factor == "Herbivores"], probs = c(0.025, 0.975)) 
    # 0.004363976 0.928455958 


### precipitation
median(output2$estimate[output2$factor == "Precipitation"]) 
    # 0.3100998
quantile(output2$estimate[output2$factor == "Precipitation"], probs = c(0.25, 0.75)) 
    # 0.1230880 0.6734638 
quantile(output2$estimate[output2$factor == "Precipitation"], probs = c(0.025, 0.975)) 
    # 0.01160812 2.12264233 


### seasonality
median(output2$estimate[output2$factor == "Seasonality"]) 
    # 0.6944969
quantile(output2$estimate[output2$factor == "Seasonality"], probs = c(0.25, 0.75)) 
    # 0.2982597 1.3690181 
quantile(output2$estimate[output2$factor == "Seasonality"], probs = c(0.025, 0.975)) 
    # 0.02790929 3.52836969 


### vegetation
median(output2$estimate[output2$factor == "Vegetation"]) 
    # 1.068133
quantile(output2$estimate[output2$factor == "Vegetation"], probs = c(0.25, 0.75)) 
    # 0.556957 1.879267 
quantile(output2$estimate[output2$factor == "Vegetation"], probs = c(0.025, 0.975)) 
    # 0.1028116 4.7286598 


### temperature
median(output2$estimate[output2$factor == "Temperature"]) 
    # 1.125604
quantile(output2$estimate[output2$factor == "Temperature"], probs = c(0.25, 0.75)) 
    # 0.4880137 2.1242250 
quantile(output2$estimate[output2$factor == "Temperature"], probs = c(0.025, 0.975)) 
    # 0.04503373 5.04414234 


### distance
# median(aD_pbi) / max(geo_dist) * 1000000 # median effect per 1000km

Find.aD <- function(MCMC_run_name){
    
    # import MCMC output
    show(load(MCMC_run_name))
    
    # extract estimated effect sizes after 20% burn-in
    aD_pbi <- as.vector(aD[201:1000]) # distance
    
    return(aD_pbi)
}

# compile aD parameters across all runs
distance_output <- sapply(my_runs, Find.aD)
aD_pbi <- as.vector(distance_output)

# summarize
mean(aD_pbi) # 1.766839
    # average effect of max geodistance = 8,700km
median(aD_pbi) # 1.428426

median(aD_pbi) / max(geo_dist) * 1000000 
    # median effect per 1000km = 0.165152

aD1000 <- median(aD_pbi) / max(geo_dist) * 1000000 

### comparisons
median(output2$estimate[output2$factor == "Herbivores"]) / aD1000 

median(output2$estimate[output2$factor == "Precipitation"]) / aD1000 

median(output2$estimate[output2$factor == "Seasonality"]) / aD1000 

median(output2$estimate[output2$factor == "Vegetation"]) / aD1000 

median(output2$estimate[output2$factor == "Temperature"]) / aD1000 


### 03. plot Figures -----------------------------------------------------------

# add plotting order
output2$factor[output2$factor == "Herbivores"] <- "1_Herbivores"
output2$factor[output2$factor == "Precipitation"] <- "2_Precipitation"
output2$factor[output2$factor == "Seasonality"] <- "3_Seasonality"
output2$factor[output2$factor == "Vegetation"] <- "4_Vegetation"
output2$factor[output2$factor == "Temperature"] <- "5_Temperature"

# add plotting colors
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]}

my_colors <-  rainbow(n = 5)

# order plotting colors
my_colors <- c(my_colors[2], my_colors[4], my_colors[5], my_colors[3], my_colors[1])


### FIGURE 5 - plot effect sizes -----------------------------------------------
library(ggplot2)
library(ggbeeswarm)

# function to plot 95% CI from bootstrapping distribution
Plot.IQR.CI <- function(x) {
    y = median(x, na.rm = T)
    ymin <- as.numeric(quantile(x, na.rm = T, probs=c(0.25, 0.75))[1])
    ymax <- as.numeric(quantile(x, na.rm = T, probs=c(0.25, 0.75))[2])
    return(c(y = y, ymin=ymin,ymax=ymax))}

Plot.95.CI <- function(x) {
    y = median(x, na.rm = T)
    ymin <- as.numeric(quantile(x, na.rm = T, probs=c(0.025, 0.975))[1])
    ymax <- as.numeric(quantile(x, na.rm = T, probs=c(0.025, 0.975))[2])
    return(c(y = y, ymin=ymin,ymax=ymax))}


dev.new(width = 6, height = 5, noRStudioGD = TRUE)

ggplot(output2, aes(x = factor, y = estimate)) + 
    geom_quasirandom(aes(color = factor), alpha = 0.25, size = 1) +
    scale_colour_manual(values = my_colors) +
    stat_summary(fun.data = Plot.95.CI, geom = "errorbar", color = "darkgray", size = 1, width = 0.25) + 
    stat_summary(fun.data = Plot.IQR.CI, geom = "pointrange", color = "black", size = 1) +
    ylim(c(0,6)) + 
    xlab("") +
    ylab("Relative effect size (aE) on \\n phytochemical arsenal differentiation") +
    theme(axis.title.x = element_text(size = 12, face = "bold")) +
    theme(axis.title.y = element_text(size = 12, face = "bold")) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(size = 12)) +
    theme(panel.background = element_blank()) +
    theme(axis.text.y = element_text(size = 12)) +
    scale_x_discrete(labels = c("Herbivores", "Precipitation", "Seasonality", "Vegetation", "Temperature"))

# quartz.save(file = "Figure5_Bedassle_effect_sizes_09_06_2022.jpg", type = "jpg", dpi = 300)


