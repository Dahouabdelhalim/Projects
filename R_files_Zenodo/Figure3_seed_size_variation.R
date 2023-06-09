### 00. import all seed size data from ImageJ particle analysis ----------------

setwd("~/Desktop/Dryad/")
dat <- read.csv(file = "2A_seed_size_data.csv", header = T, stringsAsFactors = F)

nrow(dat) # 33093
length(unique(dat$maternal_line)) # 291


### 02. import long term climate data from PRISM -------------------------------
COV_d1 <- read.csv(file = "PRISM_precip_1967_2017.csv", header = T, stringsAsFactors = F)


### 03. calculate COV for precipitation at each of the sites -------------------
Find.COV <- function(x){
    COV = sd(x)/mean(x)
    return(COV)}

# find COV for Mean Annual Precipitation MAP for each population
tapply(COV_d1$MAP_inches, COV_d1$Population, Find.COV)

precip_COV_df <- data.frame(
    precip_COV_1967 = tapply(COV_d1$MAP_inches, COV_d1$Population, Find.COV),
    Population = names(tapply(COV_d1$MAP_inches, COV_d1$Population, Find.COV))
    )

row.names(precip_COV_df) <- NULL

# remove underscores from Population names
precip_COV_df$Population <- gsub(pattern = "_", replacement = "", x = precip_COV_df$Population)

# change Huckaby naming convention
precip_COV_df$Population[precip_COV_df$Population == "HuckabyTrailhead"] <- "Huckaby"


### 04. calculate COV for seed size for each maternal line ---------------------
seed_COV_df <- data.frame(
    seed_COV = tapply(dat$seed_area, dat$maternal_line, Find.COV),
    full_line = names(tapply(dat$seed_area, dat$maternal_line, Find.COV)))

seed_COV_df$Population <- sapply(strsplit(seed_COV_df$full_line, "_"), function(x) x[1])
seed_COV_df$maternal_line <- sapply(strsplit(seed_COV_df$full_line, "_"), function(x) x[2])

row.names(seed_COV_df) <- NULL


### 05. merge seed size COV df with precipitation COV df -----------------------
sort(unique(seed_COV_df$Population)) == sort(unique(precip_COV_df$Population))

dat2 <- merge(seed_COV_df, precip_COV_df, by.x = "Population", by.y = "Population")


# import additional climate data
clim_dat <- read.csv(file = "Plantago_climate_data_simple.csv", header = T, stringsAsFactors = F)


# remove underscores from Population names
clim_dat$Population <- gsub(pattern = "_", replacement = "", x = clim_dat$Population)

# change Huckaby naming convention
clim_dat$Population[clim_dat$Population == "HuckabyTrailhead"] <- "Huckaby"

dat3 <- merge(dat2, clim_dat, by.x = "Population", by.y = "Population")


### 06. explore seed size variation --------------------------------------------
library(emmeans)
library(scales)
library(lmerTest)
library(ggplot2)
library(ggeffects)

# build linear model 
m1 <- lmer(seed_COV ~ precip_COV_1967 + AI + (1|Population), data = dat3)
summary(m1)

# build prediction data.frame
pred_df <- ggpredict(m1, terms = c("precip_COV_1967"))

dev.new()

ggplot(pred_df, aes(x, predicted)) + 
    geom_jitter(data = dat3, aes(x = precip_COV_1967, y = seed_COV), 
        width = 0.00025, size = 3, pch = 16, col = "darkgray", alpha = 0.5) +
    geom_line(aes(linetype=group, color=group), size = 1.5) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha = 0.25) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_x_continuous(breaks = c(0.18, 0.21, 0.24, 0.27, 0.3, 0.33), labels = c(0.18, 0.21, 0.24, 0.27, 0.3, 0.33), limits = c(0.18, 0.33)) +
    theme_light() +
    xlab("Precipitation COV (1967-2017)") +
    theme(axis.title.x = element_text(size = 12, face = "bold")) +
    ylab("Seed size COV") +
    theme(axis.title.y = element_text(size = 12, face = "bold")) +
    theme(legend.title = element_blank()) +
    theme(legend.position = "none") 

# quartz.save(file = "Figure3_seedCOV_and_precip_COV.jpg", type = "jpg", dpi = 300)



