# ---------------------
# Statistical analyses.
# ---------------------

# 0) Preparation.
#
# a) Save the individual worksheets in "raw_supplement.xlsx" as the following CSV files:
#   - "Flower_screens.csv",
#   - "Bee_screens.csv"
#   - "Site_surveys.csv"
#
# b) Place these files in the same folder as this R script.
#
# c) Set working directory to the directory this R script can be found in.





# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------





# 1) Initialisation.

# a) Load data.
Flower_screens <- read.csv("Flower_screens.csv")
Bee_screens <- read.csv("Bee_screens.csv")
Site_surveys <- read.csv("Site_surveys.csv")

# b) Standardise column names for site and week.
names(Flower_screens)[names(Flower_screens)=="Site"]      <- "site"
names(Bee_screens)[names(Bee_screens)=="Collection.site"] <- "site"
names(Site_surveys)[names(Site_surveys)=="Site"]          <- "site"
names(Flower_screens)[names(Flower_screens)=="Week"]      <- "week"
names(Bee_screens)[names(Bee_screens)=="Week"]            <- "week"
names(Site_surveys)[names(Site_surveys)=="week_number"]   <- "week"

# c) Load packages used.
library(lme4)
library(glmmTMB)
library(DHARMa)
library(XNomial)
library(multcomp)
library(vegan)

# d) Define some objects for convenience.
# Vector of site names.
site.vec <- c("Lansing", "McDaniels", "Whipple")

# List of the four most common bee genera, and all other bee genera combined.
group.list <- list("Apis"="Apis",
                   "Bombus"="Bombus",
                   "Ceratina"="Ceratina",
                   "Lasioglossum"="Lasioglossum")
group.list[["Others"]] <- setdiff(unique(Bee_screens$PollinatorGenus), unlist(group.list))

# Vector of parasite group names for flower and bee samples.
# Note the following index for each parasite group.
#   1. Microsporidia
#   2. Trypanosoma
#   3. Nosema bombi
#   4. Nosema ceranae
#   5. Crithidia bombi
#   6. Crithidia expoeki
#   7. Neogregarine
#   8. 3-7 combined.
parasite.name.f.vec <- c("Microsporidia_f", "Trypanosoma_f", "Nosema_bombi_f", "Nosema_ceranae_f",
                         "Crithidia_bombi_f", "Crithidia_expoeki_f", "Apicystis_f", "any_bee_parasite_f")
parasite.name.b.vec <- c("Microsporidia", "Trypanosoma", "Nosema.bombi", "Nosema.ceranae",
                         "Crithidia.bombi", "Crithidia.expoeki", "Apicystis", "any_bee_parasite_b")





# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------





# 3) Temporal trends in parasite prevalence (Figure 2, Supplementary Tables 4 and 5).

# a) Bees.
# Specify parasite of interest. Use the parasite index for convenience.
parasite <- parasite.name.b.vec[8]

# Fit GLMM model, print model summary, and test significance of week number as predictor.
model <- glmer(as.formula(paste(parasite, "~ week + (1 | site)", sep="")), family=binomial, data=Bee_screens,
               control=glmerControl(nAGQ0initStep=F, tolPwrss=1e-12))
summary(model)
drop1(model, test="Chisq")

# Generate scaled residuals, and perform Durbin-Watson test for temporal autocorrelation.
# There may be some slight difference from the reported values due to randomness when generating scaled residuals.
resid <- simulateResiduals(model)
resid <- recalculateResiduals(resid, group=Bee_screens$week)
testTemporalAutocorrelation(resid, time=sort(unique(Bee_screens$week)), plot=F)



# b) Flowers.
# Specify parasite of interest. Use the parasite index for convenience.
parasite <- parasite.name.f.vec[8]

# Fit GLMM model, print model summary, and test significance of week number as predictor.
model <- glmer(as.formula(paste(parasite, "~ week + (1 | site)", sep="")), family=binomial, data=Flower_screens,
               control=glmerControl(nAGQ0initStep=F, tolPwrss=1e-12))
summary(model)
drop1(model, test="Chisq")

# Generate scaled residuals, and perform Durbin-Watson test for temporal autocorrelation.
# There may be some slight difference from the reported values due to randomness when generating scaled residuals.
resid <- simulateResiduals(model)
resid <- recalculateResiduals(resid, group=Flower_screens$week)
testTemporalAutocorrelation(resid, time=sort(unique(Flower_screens$week)), plot=F)





# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------





# 4) Turnover in pollinator community (Supplementary Figure 1).

# a) Initialisation and data preparation.

# Specify site and week partiton.
# Site can be "Lansing", "McDaniels" or "Whipple".
# Week partitions can be 1:8, 9:16 or 17:24
site <- "Lansing"
partition <- 1:8

# Prepare data for abundance of pollinator groups (4 most common genera + all others combined)
# First, restrict the bee screens to those from the chosen site and week partition.
group.abundance.df <- Bee_screens[Bee_screens$site==site
                                  & Bee_screens$week %in% partition,
                                  "PollinatorGenus", drop=F]
# Next, prepare the five groups.
group.abundance.df$group <- NA
for(group in names(group.list)){
  group.abundance.df[group.abundance.df$PollinatorGenus %in% group.list[[group]], "group"] <- group
}
group.abundance.df$group <- factor(group.abundance.df$group, levels=names(group.list))



# c) Turnover in pollinator community.

# Perform the exact multinomial tests.
xmonte(obs=as.vector(table(group.abundance.df$group)), expr=rep(1, length(group.list)), ntrials=1e+7, statName="LLR", detail=2)

# Perform the posthoc exact binomial tests against the hypothesis of equal multinomial proportions.
for(group in names(group.list)){
  # Multiply p-value by number of groups for Bonferroni correction.
  print( binom.test( sum(group.abundance.df$group==group), nrow(group.abundance.df), p=1/length(group.list) )$p.value * length(group.list) )
}





# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------



# 5) Parasite prevalence in each bee group, and posthoc tests (Supplementary Figure 2, Supplementary Tables 6, 7 and 8).

# a) Initialisation and data preparation.

# Specify parasite of interest. Use the parasite index for convenience.
# For Crithidia.expoeki (parasite = 6), exclude Ceratina since there were no positives.
# The fastest way to do so is to remove Ceratina from Bee_screens and group.list.
parasite <- parasite.name.b.vec[8]

# Prepare data on parasite prevalence in bees.
prevalence.df <- Bee_screens[, c("site", "week", "PollinatorGenus", parasite)]
prevalence.df$group <- NA
for(group in names(group.list)){
  prevalence.df[prevalence.df$PollinatorGenus %in% group.list[[group]], "group"] <- group
}
prevalence.df$group <- factor(prevalence.df$group, levels=names(group.list))

# Shift week for each group such that the median week number among all samples of that group is week 0.
week.median.vec <- rep(NA, length(group.list)); names(week.median.vec) <- names(group.list)
for(group in names(group.list)){
  week.median.vec[group] <- median( prevalence.df[prevalence.df$PollinatorGenus %in% group.list[[group]], "week"] )
}
prevalence.df$week.shifted <- prevalence.df$week - week.median.vec[prevalence.df$group]



# b) Effects of bee groups and week number (shifted) on parasite prevalence.

# Fit GLMM model, and print model summary.
model <- glmmTMB(as.formula(paste(parasite, "~ week.shifted*group + (1|site)", sep="")), family=binomial, data=prevalence.df)
summary(model)
# Test significance of the interaction.
drop1(model, test="Chisq")
# Test significance of the main effects of bee group.
drop1(model, ~group, test="Chisq")



# c) Pairwise contrast between bee groups. Uses glht() for multiple testing corrections.

# This is how to make glht() work for glmmTMB models.
glht_glmmTMB <- function (model, ..., component="cond") {
  glht(model, ...,
       coef. = function(x) fixef(x)[[component]],
       vcov. = function(x) vcov(x)[[component]],
       df = NULL)
}

# May need to re-run this until the "Completion with error > abseps" warning disappears, since there is some randomness involved.
summary(glht_glmmTMB(model, mcp(group="Tukey")))



# d) Temporal trend of each group. Uses glht() for multiple testing corrections.

# Define linear contrast matrix.
linfct <- matrix(0, nrow=5, ncol=length(fixef(model)$cond))
rownames(linfct) <- names(group.list)
colnames(linfct) <- names(fixef(model)$cond)
linfct[, "week.shifted"] <- 1
for(group in rownames(linfct)[-1]){
  linfct[group, paste("week.shifted:group", group, sep="")] <- 1
}

# Use linear contrast matrix to test significance of temporal trends for each group.
summary(glht_glmmTMB(model, linfct=linfct))





# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------





# 6) Tests involving bee Shannon diversity.
# - Variation in bee Shannon diversity across sites (Supplementary Figure 5a).
# - Effects on bee Shannon diversity on parasite prevalence in bees (Figure 3d, Supplementary Table 9).

# a) Initialisation and data preparation.

# Specify parasite of interest. Use the parasite index for convenience.
parasite <- parasite.name.b.vec[8]

# Decide on subsample size when calculating Shannon diversity. If we choose not to subsample, set it to NA.
subsample.size <- 36

# Exclude Apis and Bombus when calculating prevalence?
exclude.common <- F


# Prepare data on parasite prevalence in bees and bee Shannon diversity.
combined.df <- c()
# Loop over sites.
for(site in site.vec){
  # Unique weeks for that site.
  week.vec <- sort(unique(Bee_screens[Bee_screens$site==site, "week"]))
  
  # Loop over weeks.
  for(week in week.vec){
    # Restrict the bee screens to those from that site and week.
    Bee_screens.subset <- Bee_screens[Bee_screens$site==site & Bee_screens$week==week, c("PollinatorName", "PollinatorGenus", parasite)]  
    
    # Calculate number of infected and infected bees.
    # Do we want to include Apis and Bombus?
    if(!exclude.common){
      # Yes.
      infected=sum(Bee_screens.subset[, parasite])
      uninfected=sum(!Bee_screens.subset[, parasite])
    } else {
      # No, exclude them.
      infected=sum(Bee_screens.subset[!(Bee_screens.subset$PollinatorGenus %in% c("Apis", "Bombus")), parasite])
      uninfected=sum(!Bee_screens.subset[!(Bee_screens.subset$PollinatorGenus %in% c("Apis", "Bombus")), parasite])
    }
    
    # Should we subsample?
    if(!is.na(subsample.size)){
      
      # If yes, is there enough samples given the subsample size?
      if(nrow(Bee_screens.subset)>=subsample.size){
        
        # If yes, subsample and calculate the Shannon diversity. Repeat 100 times, and take the median.
        diversity.vec <- c()
        for(i in 1:100){
          diversity.vec <- c( diversity.vec,
                              diversity(rrarefy(table(Bee_screens.subset$PollinatorName),
                                                subsample.size)) )
        }
        
        # Add new row to the dummy data frame.
        combined.df <- rbind(combined.df,
                             data.frame(shannon=median(diversity.vec),
                                        infected=infected,
                                        uninfected=uninfected,
                                        site=site,
                                        week=week))
      }
      
    } else {
      # If no, don't subsample.
      
      # Add new row to the dummy data frame.
      combined.df <- rbind(overall.prevalence.df,
                           data.frame(shannon=diversity(table(Bee_screens.subset$PollinatorName)),
                                      infected=sum(Bee_screens.subset[, parasite]),
                                      uninfected=sum(!Bee_screens.subset[, parasite]),
                                      site=site,
                                      week=week))
    }
  }
}



# c) Variation in Shannon diversity across sites.
# There may be some slight difference from the reported values due to randomness when subsampling.
oneway.test(shannon ~ site, data=combined.df)



# d) Effects of Shannon diversity on parasite prevalence.
# Fit GLMM, print model summary, and test significance of Shannon diversity as predictor.
# There may be some slight difference from the reported values due to randomness when subsampling.
model <- glmer(cbind(infected, uninfected) ~ shannon + (1|site), family=binomial, data=combined.df,
               control=glmerControl(nAGQ0initStep=F, tolPwrss=1e-12))
summary(model)
drop1(model, test="Chisq")





# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------





# 7) Temporal trends in floral abundance (Figure 4a).

# a) Data preparation.

# Prepare floral abundance data.
total.floral.abundance.df <- c() 
# Loop over sites.
for(site in site.vec){
  # Unique weeks for that site.
  week.vec <- sort(unique(Site_surveys[Site_surveys$site==site, "week"]))
  
  # Loop over weeks.
  for(week in week.vec){
    # Unique quadrants for that site and week.
    zone.vec <- sort(unique(Site_surveys[Site_surveys$site==site & Site_surveys$week==week, "Zone"]))
    
    # Loop over quadrants.
    for(zone in zone.vec){
      
      # Sum the floral abundance of all species in the quadrant.
      # Add new row to the dummy data frame.
      total.floral.abundance.df <- rbind(total.floral.abundance.df,
                                         data.frame(site=site,
                                                    week=week,
                                                    zone.count=round(sum( Site_surveys[Site_surveys$site==site
                                                                                       & Site_surveys$week==week
                                                                                       & Site_surveys$Zone==zone,
                                                                                       "floral.abundance"] ))))
    }
  }
}



# b) Temporal trends in floral abundance.

# Fit negative binomial GLMM, print model summary, and test significance of week number as predictor.
model <- glmmTMB(zone.count ~ week + (1|site), data=total.floral.abundance.df, family=nbinom2)
summary(model)
drop1(model, test="Chisq")

# Generate scaled residuals, and perform Durbin-Watson test for temporal autocorrelation.
# There may be some slight difference from the reported values due to randomness when generating scaled residuals.
resid <- simulateResiduals(model)
resid <- recalculateResiduals(resid, group=total.floral.abundance.df$week)
testTemporalAutocorrelation(resid, time=sort(unique(total.floral.abundance.df$week)), plot=F)





# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------





# 8) Effects on floral unit abundance on parasite prevalence (Figure 4b, Supplementary Table 11).

# a) Initialisation and data preparation.

# Specify parasite of interest. Use the parasite index for convenience.
parasite <- parasite.name.f.vec[8]

# Prepare data for parasite prevalence on flowers
overall.prevalence.df <- c()
# Loop over sites.
for(site in site.vec){
  # Unique weeks for that site.
  week.vec <- sort(unique(Flower_screens[Flower_screens$site==site, "week"]))
  
  total.vec <- c()
  positive.vec <- c()
  # Loop over weeks.
  for(week in week.vec){
    total.vec    <- c(total.vec,
                      sum(Flower_screens$site==site & Flower_screens$week==week))
    positive.vec <- c(positive.vec,
                      sum(Flower_screens[Flower_screens$site==site & Flower_screens$week==week, parasite]))
  }
  # Add new row to the dummy data frame.
  overall.prevalence.df <- rbind(overall.prevalence.df,
                                 data.frame(site=site, week=week.vec, total=total.vec, prevalence=positive.vec/total.vec))
}

# Prepare floral abundance data.
# Unlike the previous section, we now only have one value per site and week, i.e. we average over all 3 quadrants.
total.floral.abundance.df <- c()
# Loop over sites.
for(site in site.vec){
  # Unique weeks for that site.
  week.vec <- sort(unique(Site_surveys[Site_surveys$site==site, "week"]))
  
  zone.count.vec <- c()
  # Loop over weeks.
  for(week in week.vec){
    zone.count.vec <- c(zone.count.vec,
                        sum(Site_surveys[Site_surveys$site==site & Site_surveys$week==week, "floral.abundance"])/3)
  }
  # Add new row to the dummy data frame.
  total.floral.abundance.df <- rbind(total.floral.abundance.df,
                                     data.frame(site=site, week=week.vec, zone.count=zone.count.vec))
}



# c) Effects of log10(floral abundance) on parasite prevalence.

# Merge the prevalence and abundance data frame.
combined.df <- merge(overall.prevalence.df, total.floral.abundance.df,
                     by=c("site", "week"))

# Fit GLMM, print model summary, and test significance of log10(floral abundance) as predictor.
model <- glmer(prevalence ~ log10(zone.count) + (1|site), family=binomial, data=combined.df, weights=total,
               control=glmerControl(optimizer="Nelder_Mead", nAGQ0initStep=F, tolPwrss=1e-14))
summary(model)
drop1(model, test="Chisq")





# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------





# 9) Tests involving flower Shannon diversity.
# - Variation in flower Shannon diversity across sites (Supplementary Figure 5b).
# - Temporal trends in flower Shannon diversity.
# - Effects of flower Shannon diversity on parasite prevalence on flowers (Supplementary Table 10).



# a) Initialisation and data preparation.

# Specify parasite of interest. Use the parasite index for convenience.
parasite <- parasite.name.f.vec[8]

# Prepare data for parasite prevalence on flowers.
overall.prevalence.df <- c()
# Loop over sites.
for(site in site.vec){
  # Unique weeks for that site.
  week.vec <- sort(unique(Flower_screens[Flower_screens$site==site, "week"]))
  
  total.vec <- c()
  positive.vec <- c()
  # Loop over weeks.
  for(week in week.vec){
    # Total number of samples.
    total.vec    <- c(total.vec,
                      sum(Flower_screens$site==site & Flower_screens$week==week))
    # Number of positive samples.
    positive.vec <- c(positive.vec,
                      sum(Flower_screens[Flower_screens$site==site & Flower_screens$week==week, parasite]))
  }
  
  # Add new row to the dummy data frame.
  overall.prevalence.df <- rbind(overall.prevalence.df,
                                 data.frame(site=site, week=week.vec, total=total.vec, prevalence=positive.vec/total.vec))
}

# Prepare data for flower Shannon diversity.
shannon.df <- c()
# Loop over sites.
for(site in site.vec){
  # Unique weeks for that site.
  week.vec <- sort(unique(Site_surveys[Site_surveys$site==site, "week"]))
  
  shannon.vec <- c()
  # Loop over weeks.
  for(week in week.vec){
    # Restrict the site surveys to that site and week.
    Site_surveys.subset <- Site_surveys[Site_surveys$site==site & Site_surveys$week==week,]
    # We convert FlowerTaxonomicName from factor to character, so as not to generate NA for flower species that do not show up for that site and week.
    shannon.vec <- c(shannon.vec,
                     diversity(tapply(Site_surveys.subset$floral.abundance, as.character(Site_surveys.subset$FlowerTaxonomicName), sum)))
  }
  
  # Add new row to the dummy data frame.
  shannon.df <- rbind(shannon.df,
                      data.frame(site=site, week=week.vec, shannon=shannon.vec))
}



# b) Variation in Shannon diversity across sites.

oneway.test(shannon ~ site, shannon.df)



# c) Temporal trends in Shannon diversity.

# Fit LMM, print model summary, and test significance of week number as predictor.
model <- lmer(shannon ~ week + (1|site), data=shannon.df)
summary(model)
drop1(model, test="Chisq")

# Generate scaled residuals, and perform Durbin-Watson test for temporal autocorrelation.
# There may be some slight difference from the reported values due to randomness when generating scaled residuals.
resid <- simulateResiduals(model)
resid <- recalculateResiduals(resid, group=shannon.df$week)
testTemporalAutocorrelation(resid, time=week.vec, plot=F)



# d) Effects of Shannon diversity on parasite prevalence.

# Merge prevalence and Shannon diversity data.
combined.df <- merge(overall.prevalence.df, shannon.df,
                     by=c("site", "week"))

# Fit GLMM, print model summary, and test significance of Shannon diversity as predictor.
model <- glmer(prevalence ~ shannon + (1|site), family=binomial, data=combined.df, weights=total,
               control=glmerControl(nAGQ0initStep=F, tolPwrss=1e-12))
summary(model)
drop1(model, test="Chisq")