
# --------------------------------------------------------------------------------
# Install and load packages

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

packages <- c("nlme", "emmeans", "tidyverse", "car")
ipak(packages)

# suppress correlation matrix from summary output for gls
assignInNamespace("print.correlation", function(x, title) return(), ns="nlme")

# global options
options(scipen = 20)


# --------------------------------------------------------------------------------
# Load data

setwd("path-to-directory") # edit file path for working directory

dat <- read.csv("fox_linear_volume_data.csv", stringsAsFactors=TRUE)
glimpse(dat)


# --------------------------------------------------------------------------------
# Clean data

dat <- dat %>%
    mutate(endocranial_volume_cr = endocranial_volume^(1/3),
                zygomatic_width_norm = zygomatic_width / centroid_size,
                cranial_vault_width_norm = cranial_vault_width / centroid_size,
                upper_jaw_width_norm = upper_jaw_width / centroid_size,
                total_skull_length_norm = total_skull_length / centroid_size,
                snout_length_norm = snout_length / centroid_size,
                cranial_vault_height_norm = cranial_vault_height / centroid_size,
                endocranial_volume_cr_norm = endocranial_volume_cr / centroid_size
                )
glimpse(dat)

# --------------------------------------------------------------------------------
#Assess heteroskedascity of variances

leveneTest(centroid_size ~ population, data = dat)

# GLS model testing the hypothesis of different centroid size between populations

model_c <- gls(centroid_size ~ population + sex,
               weights = varIdent(form = ~ 1 | population),
               na.action = "na.exclude",
               control=list(max.iter = 100),
               data = dat)

summary(model_c)
plot(model_c)

c_means <- emmeans(model_c, specs = ~ population)
contr_c <- summary(rbind(contrast(c_means, method = "pairwise", adjust = "none")), infer = TRUE)
contr_c$p.value.adj <- p.adjust(contr_c$p.value, method = "holm")
contr_c


# --------------------------------------------------------------------------------
# GLS model comparing size corrected values (divided by centroid size) among populations

norm_dat_long <- dat %>%
  pivot_longer(
    cols = ends_with("norm"),
    names_to = "response",
    values_to = "Y"
  )

glimpse(norm_dat_long)

model1 <- gls(Y ~ response * (population + sex),
              weights = varIdent(form = ~ 1 | response),  
              correlation = corSymm(form = ~ 1 | id),
              na.action = "na.exclude",
              data = norm_dat_long)

summary(model1)
plot(model1)

main_population_means <- emmeans(model1, specs = ~ population | response)
contr1 <- summary(rbind(contrast(main_population_means, method = "pairwise", adjust = "none")), infer = TRUE)
contr1$p.value.adj <- p.adjust(contr1$p.value, method = "holm")
contr1


# --------------------------------------------------------------------------------
# GLS model with population * sex interaction term (looking at the effects of sexual dimorphism)

raw_dat_long <- dat %>%
  pivot_longer(
    cols = c("zygomatic_width", "cranial_vault_width", "upper_jaw_width", 
             "total_skull_length", "snout_length", "cranial_vault_height", "endocranial_volume_cr"),
    names_to = "response",
    values_to = "Y"
  )

model2 <- gls(Y ~ response * (population * sex),
              weights = varIdent(form = ~ 1 | response), 
              correlation = corSymm(form = ~ 1 | id),
              na.action = "na.exclude",
              data = raw_dat_long)

summary(model2)
plot(model2) 

interaction_means <- emmeans(model2, specs = ~ sex * population | response, mode = "df.error")
# get contrasts between sexes for each level of population and response
tab <- summary(rbind(contrast(interaction_means, method = "revpairwise", adjust = "none"), adjust = "none"), infer = TRUE)

# get row indices of just the within population contrasts
row_index <- c(which(str_count(tab$contrast, "Domesticated") == 2), 
               which(str_count(tab$contrast, "Unselected") == 2), 
               which(str_count(tab$contrast, "Wild") == 2))

# subset rows
tab_sub <- tab[row_index, ]
# adjust p-values
tab_sub$p.value.adj <- p.adjust(tab_sub$p.value, method = "holm") 
tab_sub

##---------------------------------------------------------------------------------------
#Get percent differences by population for model 1
##numerators
DU.diff <- subset(contr1, contrast == "Domesticated - Unselected", select = c(estimate, lower.CL, upper.CL))
DW.diff <- subset(contr1, contrast == "Domesticated - Wild", select = c(estimate, lower.CL, upper.CL))
UW.diff <- subset(contr1, contrast == "Unselected - Wild", select = c(estimate, lower.CL, upper.CL))

#denominators
model.means <- summary(main_population_means)

Domest.means <- subset(model.means, population == "Domesticated", select = c(emmean, lower.CL, upper.CL))
Unselected.means <- subset(model.means, population == "Unselected", select = c(emmean, lower.CL, upper.CL))
Wild.means <- subset(model.means, population == "Wild", select = c(emmean, lower.CL, upper.CL))

Domest.means <- as.matrix(Domest.means, nrow = 6, ncol = 3)
Unselected.means <- as.matrix (Unselected.means, nrow = 6, ncol = 3)
Wild.means <- as.matrix(Wild.means, norw = 6, ncol = 3)

DU.list <- list(Domest.means,Unselected.means)
denom.DU <- apply(simplify2array(DU.list), 1:2, mean)
DW.list <- list(Domest.means,Wild.means)
denom.DW <- apply(simplify2array(DW.list), 1:2, mean)
UW.list <- list(Unselected.means,Wild.means)
denom.UW <- apply(simplify2array(DU.list), 1:2, mean)

#%differences

DUp.diff <- (DU.diff/denom.DU)*100
DWp.diff <- (DW.diff/denom.DW)*100
UWp.diff <- (UW.diff/denom.UW)*100

##----------------------------------------------------------------------------------------------------
# save data and fitted models
save(dat, model_c, model1, model2, tab_sub, file = "gls_models_and_data.Rdata", compress = "gzip")


