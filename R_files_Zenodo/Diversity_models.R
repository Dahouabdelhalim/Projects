## Diversity equations
library(glmmTMB)
library(vegan)

full_dataset <- read.csv("OTU_matrix_by_species.csv")
terpene_dataset <- read.csv("Terpene_OTU.csv")

###### Full dataset models ######

## Generate diversity metrics for maternal family models
full_dataset <- cbind(full_dataset[c(1:13)],
                              Richness = specnumber(full_dataset[,c(14:ncol(full_dataset))]),
                              Shannon = diversity(full_dataset[,c(14:ncol(full_dataset))]))

## Create column to estimate residual effects
full_dataset$Residual <- factor(1:nrow(full_dataset))

## Calculate True Diversity
full_dataset$ExpSpec <- exp(full_dataset$Shannon)
full_dataset[full_dataset$Richness == 0,]$ExpSpec <- 0

## Total tips as covariate
full_dataset$Tot_tips <- rowSums(full_dataset[,13:41])

## Models to test effects on species richness
full_data_rich_mod <- glmmTMB(Richness ~ Site + Tot_tips + (1 | Family) + (1 | Grid) + (1 | Residual), data = full_dataset, family = poisson(link = log))
full_data_rich_mod2 <- glmmTMB(Richness ~ Site + Tot_tips + (1 | Grid) + (1 | Residual), data = full_dataset, family = poisson(link = log))
full_data_rich_mod3 <- glmmTMB(Richness ~ Tot_tips + (1 | Family) + (1 | Grid) + (1 | Residual), data = full_dataset, family = poisson(link = log))

## Test for effect of Family
anova(full_data_rich_mod, full_data_rich_mod2, test = "LRT")

## Test for effect of Site
anova(full_data_rich_mod, full_data_rich_mod3, test = "LRT")

## Models to test effects on species diversity
full_data_diversity_mod <- glmmTMB(ExpSpec ~ Site + Tot_tips + (1 | Family) + (1 | Grid), data = full_dataset)
full_data_diversity_mod2 <- glmmTMB(ExpSpec ~ Site + Tot_tips + (1 | Grid), data = full_dataset)
full_data_diversity_mod3 <- glmmTMB(ExpSpec ~ Tot_tips + (1 | Family) + (1 | Grid), data = full_dataset)

## Test for effect of Family
anova(full_data_diversity_mod, full_data_diversity_mod2)

## Test for effect of Site
anova(full_data_diversity_mod, full_data_diversity_mod3)

###### Terpene models ######

## Generate diversity metrics for terpene models
terpene_dataset$Richness <- vegan::specnumber(terpene_dataset[,c(14:42)])
terpene_dataset$SpecDiv <- vegan::diversity(terpene_dataset[,c(14:42)])
terpene_dataset$ExpSpec <- exp(terpene_dataset$SpecDiv)
terpene_dataset[terpene_dataset$Richness == 0,]$ExpSpec <- 0
terpene_dataset$Residual <- 1:nrow(terpene_dataset)

## Models to test effects on species richness
terpene_data_rich_mod <- glmmTMB(Richness ~ Tot_tips + Chemodiversity + Site + (1 | Grid) + (1 | Residual), data = terpene_dataset, family = "poisson")
terpene_data_rich_mod2 <- glmmTMB(Richness ~ Tot_tips +  Site + (1 | Grid) + (1 | Residual), data = terpene_dataset, family = "poisson")
terpene_data_rich_mod3 <- glmmTMB(Richness ~ Tot_tips + Chemodiversity + (1 | Grid) + (1 | Residual), data = terpene_dataset, family = "poisson")

## Test for effect of Family
anova(terpene_data_rich_mod, terpene_data_rich_mod2)

## Test for effect of Site
anova(terpene_data_rich_mod, terpene_data_rich_mod3)

## Models to test effects on species diversity
terpene_data_diversity_mod <- glmmTMB(ExpSpec ~ Tot_tips + Chemodiversity + Site + (1 | Grid), data = terpene_dataset)
terpene_data_diversity_mod2 <- glmmTMB(ExpSpec ~ Tot_tips + Site + (1 | Grid), data = terpene_dataset)
terpene_data_diversity_mod3 <- glmmTMB(ExpSpec ~ Tot_tips + Chemodiversity + (1 | Grid), data = terpene_dataset)

## Test for effect of Family
anova(terpene_data_diversity_mod, terpene_data_diversity_mod2)

## Test for effect of Site
anova(terpene_data_diversity_mod, terpene_data_diversity_mod3)


