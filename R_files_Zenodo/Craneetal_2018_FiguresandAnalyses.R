# Generate Figures, Tables, and Analyses associated with
# Crane, Cox, Kisare, and Patek, 2018


# load required packages
library(ggplot2)
library(tidyverse)
library(RVAideMemoire)
library(lme4)
library(stringr)


# Import data. All data are available for download on DataDryad
feedingtrials_all <- read_csv("Craneetal_2018_FeedingTrials.csv")
strikelocations_all <- read_csv("Craneetal_2018_StrikeLocations.csv")
ninjabottrials <- read_csv("Craneetal_2018_NinjabotTrials.csv")


# filter feedingtrials as described in methods
feedingtrials <- feedingtrials_all %>%
  filter(success == 1,   # remove unsuccessful strike sequences
         is.na(notes),   # remove one lost video
         unknown_count / total_strikes < 0.5)   # remove strike sequences where the location 
                                                # of fewer than half of strikes can be identified
strikelocations <- strikelocations_all %>%
  filter(trial_id %in% feedingtrials$trial_id,   # filter strike data according to the same parameters
         rewatch == 0)   # exclude rewatches of videos






#################################################
#################### Methods ####################
#################################################


#### Mantis shrimp carapace length and body length correlation ####
# Consider one measurement for each mantis shrimp
bodysizetest <- feedingtrials %>%   
  group_by(mantis_id) %>%
  summarize(body = mantis_length_mm[1],
            carapace = mantis_carapace_mm[1])

# Generate and report linear model of carapace length in terms of body length
summary(lm(carapace ~ body, data = bodysizetest))
nrow(bodysizetest)   # n





#### test size-matching of mantis shrimp and snails ####
# Kruskal Wallis comparing the total number of strikes mantis shrimp delivered to each shell shape in 2015
kruskal.test(data = filter(feedingtrials, year == 2015), total_strikes ~ as.factor(snail_genus))
sum(feedingtrials$year == 2015)   # n








############# Figure S1 ###################
# Size matching of snails and mantis shrimp
# Figure S1A
feedingtrials %>%
  mutate(snail_genus = factor(snail_genus, ordered = TRUE, levels = c("Nerita", "Cenchritis", "Cerithium"))) %>%  # convert snail genus to an ordered factor
  ggplot(aes(x = mantis_carapace_mm, y = snail_length_mm, col = as.factor(year))) +
    geom_smooth(method = lm, se = FALSE) +
    geom_point() + 
    facet_grid(. ~ snail_genus) +
    labs(x = "carapace length (mm)",
         y = "snail length (mm)",
         col = "year") +
    theme_classic() +
    scale_color_manual(values = c("darkgray", "black"))



# Figure S1B
feedingtrials %>%
  mutate(snail_genus = factor(snail_genus, ordered = TRUE, levels = c("Nerita", "Cenchritis", "Cerithium"))) %>%  # convert snail genus to an ordered factor
  ggplot(aes(x = snail_genus, y = total_strikes, col = as.factor(year))) +
    geom_boxplot() +
   labs(x = "snail genus",
        y = "total number of strikes",
        col = "year")  +
    scale_y_continuous(expand = c(0,0), limits = c(0, 475)) +
   theme_classic() +
   scale_color_manual(values = c("darkgray", "black"))











#### Preliminary test comparing two groups of Cerithium spp. snails ####
# Test whether the proportion of strikes to each region of a shell differed significantly between
# mantis shrimp tested in 2014 and 2015 to test validity of collapsing across all Cerithium 
# strike sequences
sum(feedingtrials$snail_genus == "Cerithium")   # n

# compare proportion of strikes to the aperture between 2014 and 2015
wilcox.test(aperture_proportion ~ year, data = filter(feedingtrials, snail_genus == "Cerithium"))
# compare proportion of strikes to the whorls between 2014 and 2015
wilcox.test(whorls_proportion ~ year, data = filter(feedingtrials, snail_genus == "Cerithium"))
# compare proportion of strikes to the apex between 2014 and 2015
wilcox.test(apex_proportion ~ year, data = filter(feedingtrials, snail_genus == "Cerithium"))







#### Preliminary test comparing effect of mantis shrimp sex ####
# Test whether the proportion of strikes to each region of a shell differed significantly between
# male and female mantis shrimp fed Cerithium spp. snails
sum(feedingtrials$snail_genus == "Cerithium")   # n
sum(feedingtrials$snail_genus == "Cerithium" & feedingtrials$mantis_sex == "M")   # male mantis shrimp
sum(feedingtrials$snail_genus == "Cerithium" & feedingtrials$mantis_sex == "F")   # female mantis shrimp


# compare proportion of strikes to the aperture between male and female mantis shrimp
wilcox.test(aperture_proportion ~ mantis_sex, data = filter(feedingtrials, snail_genus == "Cerithium"))
# compare proportion of strikes to the whorls between male and female mantis shrimp
wilcox.test(whorls_proportion ~ mantis_sex, data = filter(feedingtrials, snail_genus == "Cerithium"))
# compare proportion of strikes to the apex between male and female mantis shrimp
wilcox.test(apex_proportion ~ mantis_sex, data = filter(feedingtrials, snail_genus == "Cerithium"))










#### Figure 3 ####
# The amount of damage to the aperture and apex associated with where mantis shrimp ate
# Figure 3A
feedingtrials %>%
  filter(feeding_location == "aperture",   # filter to only mantis shrimp that ate from the aperture
         !is.na(aperture_damage)) %>%      # exclude snail whose aperture damage was not measured
  ggplot(aes(x = aperture_damage/4)) +
    geom_bar() +
    scale_x_continuous(limits = c(-0.2,1.4),
                       expand = c(0,0), breaks = c(0,.25,.5,.75,1, 1.25)) +
    labs(x = "aperture damage rotations",
         y = "number of snails",
         title = "shells eaten from the aperture") +
    scale_y_continuous(limits = c(0,4.5), expand = c(0,0)) +
    theme_classic() +
    geom_vline(xintercept = 0.625, color = "red", linetype = "dashed")


# Figure 3B
feedingtrials$aperture_damage[feedingtrials$aperture_damage>4] <- 5  # collapse all damage measurements above one rotation
feedingtrials %>%
  filter(feeding_location != "aperture",  # filter out mantis shrimp that ate from the aperture
         !is.na(aperture_damage)) %>%     # exclude snails whose aperture damage was not measured
  ggplot(aes(x = aperture_damage/4)) +
  geom_bar() +
  scale_x_continuous(limits = c(-0.2,1.4), 
                     expand = c(0,0), breaks = c(0,.25,.5,.75,1, 1.25)) +
  labs(x = "damage rotations around axis", 
       y = "number of snails",
       title = "shells not eaten from the aperture") +
  scale_y_continuous(limits = c(0,28), expand = c(0,0)) +
  theme_classic() +
  geom_vline(xintercept = 0.625, color = "red", linetype = "dashed") +
  theme(text = element_text(size = 10))


# Figure 3C
feedingtrials %>%
  filter(feeding_location == "apex",   # filter to only shells consumed from the apex
         !is.na(apex_loss_proportion)) %>%    # exclude shells that were not measured
  ggplot(aes(x = sample(1:length(apex_loss_proportion), replace = FALSE), # generate a random order for the x-axis
             y = apex_loss_proportion)) +
    geom_point() +
    geom_hline(yintercept = .2, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0) +
    labs(x = "individual snails", y = "damage to apex / shell length", 
       title = "shells eaten from the apex") +
    scale_y_continuous(limits = c(-0.1,0.57), expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank())

# Figure 3D
feedingtrials %>%
  filter(feeding_location != "apex",   # exclude shells consumed from the apex
         !is.na(apex_loss_proportion)) %>%   # exclude shells that were not measured
  ggplot(aes(x = sample(1:length(apex_loss_proportion), replace = FALSE), # generate a random order for the x-axis
             y = apex_loss_proportion)) +
  geom_point() +
  geom_hline(yintercept = .2, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0) +
  labs(x = "individual snails", y = "damage to apex / shell length", 
       title = "shells not eaten from the apex") +
  scale_y_continuous(limits = c(-0.1,0.57), expand = c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())








##############################################################################
#### Do mantis shrimp preferentially strike shells in specific locations? ####
##############################################################################


### Figure 4 ####
# Generate each panel of Figure 4, which shows every strike location
# for every feeding trial, separated by snail genus

figure4 <- function(genus){
  genusDF <- filter(feedingtrials, snail_genus == genus)   # subset feedingtrials by the specified genus
  
  genusstrikesDF <- strikelocations %>%
    filter(trial_id %in% genusDF$trial_id,   # subset strikelocations by specified genus
           !is.na(strike_location)) %>%      # exclude strikes that did not impact the shell
    mutate(minutes = time_sec / 60,          # convert time in seconds to minutes
           strike_location = factor(strike_location,
                                    levels = c("aperture", "whorls","apex", "unknown"),
                                    ordered = T))     # convert strike location to an ordered categorical variable
           
  
  # reorder data by the total time of the strike sequence
  # the numbers of this ordering will form the basis of the y-axis
  genusDF <- genusDF[order(genusDF$total_time_sec),]   # generate order
  
  genusstrikesDF$id_num <- rep(0, nrow(genusstrikesDF))  # initialize blank column for order
  
  # label every strike according to the order generated above
  for(i in 1:nrow(genusDF)){    
    genusstrikesDF$id_num <- ifelse(
      genusstrikesDF$trial_id == as.character(genusDF$trial_id[i]), 
      i, genusstrikesDF$id_num)
  }
  rm(i)
  

  # Because the triangles that identify apex strikes center slightly differently, we need to shift them to sit in line with the other points
  # Depending on the scaling of the figure, how much you need to shift the triangles by will vary
  genusstrikesDF <- mutate(genusstrikesDF, 
                           id_num = ifelse(strike_location == "apex",
                                           id_num - 0.1, id_num))              
  
  # generate figure
  ggplot(genusstrikesDF, aes(x = minutes, y = id_num, 
                             shape = strike_location, col = strike_location)) +
    geom_point() +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 140), 
                       breaks = seq(0, 150, 30)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, nrow(genusDF) + 1)) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_colour_manual(values=c("#FF5553", "#E0A100", "#4881B2", "#000000")) +
    scale_shape_manual(values=c(19,15,17,4)) +
    labs(x = "time (min)",
         y = "individual mantis shrimp",
         title = str_c(as.character(genus), " spp."))
}

# Figure 4A
figure4(genus = "Nerita")
# Figure 4B
figure4(genus = "Cenchritis")
# Figure 4C
figure4(genus = "Cerithium")








#### Figure 5 ####
# Illustrate how the frequency of strikes to each locations differs from expected (5A)
# and how it varies across shell shapes (5B)

# Figure 5A
# calculate ratio of observed to expected strikes for every region.
# Add one strike to each region to account for strike sequences where 
# specific regions were never struck.
g.feedingtrials <- feedingtrials %>%
  mutate(aperture = aperture_count + 1,  # Add one strike to the tallies of strikes.
         apex = apex_count + 1,
         whorls = whorls_count +1,
         totalknown = total_known_strikes +3) %>%  # Add three strikes to the total tally, because we have added one strike to each region.
  mutate(g.aperture_proportion = aperture / totalknown,  # Calculate proportion of strikes to each region.
         g.apex_proportion = apex / totalknown,
         g.whorls_proportion = whorls / totalknown) %>%
  mutate(aperture_o2e = g.aperture_proportion / aperture_expected,   # Compare observed proportion of strikes to what would be predicted by surface area.
         apex_o2e = g.apex_proportion / apex_expected,
         whorls_o2e = g.whorls_proportion / whorls_expected)


# reogranize data to prepare for figure 5A
byregion <- g.feedingtrials %>%
  select(mantis_id, snail_genus,  apex_o2e,  aperture_o2e, whorls_o2e) %>%
  gather(apex_o2e,  aperture_o2e, whorls_o2e, key = location, value = obs.exp) %>%   # Give strikes to each region their own row
  mutate(
    location = factor(str_replace(location, "_o2e", ""), ordered = T,
                    levels = c("aperture", "whorls", "apex")),   # Convert the strike regions to an ordered factor.
    snail_genus = factor(snail_genus, ordered = T,
                         levels = c("Nerita", "Cenchritis", "Cerithium")),  #C Convert snail genus to an ordered factor.
    logo2e = log(base = 10, obs.exp)    # log transform the observed:expected ratio
  )


# generate figure 5A
ggplot(byregion, aes(x = location, y = logo2e, fill = location)) +
  geom_boxplot(outlier.size = 0.9) +
  labs(x = "shell region",
       y = "observed / expected strikes") +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line.x = element_blank()) +
  scale_fill_manual(values = c("#FF5553", "#FFB700", "#4881B2")) +
  facet_grid(. ~ snail_genus)
rm(byregion)




# Figure 5B
# filter to just mantis shrimp that successfully handled all three shell shapes
feedingtrials_paired <- feedingtrials %>%
  filter(year == "2015",  # filter to only 2015 experiments
         mantis_id %in% feedingtrials$mantis_id[feedingtrials$snail_genus == "Cerithium"],  # filter to animals that successfully handled Cerithium
         mantis_id %in% feedingtrials$mantis_id[feedingtrials$snail_genus == "Cenchritis"], # filter to animals that successfully handled Cenchritis
         mantis_id %in% feedingtrials$mantis_id[feedingtrials$snail_genus == "Nerita"])     # filter to animals that successfully handled Nerita


# Pair data and generate figure 5B
feedingtrials_paired %>%
  # reorganize data so every region for every strike sequence has its own row
  gather(apex_proportion, whorls_proportion, aperture_proportion, key = location, value = proportion) %>%
  
  # rename strike locations to be ordered factors
  mutate(location = str_replace(location, "_proportion", "")) %>%
  mutate(location = factor(location, levels = c("aperture", "whorls", "apex"),
                         ordered = TRUE),
         snail_genus = factor(snail_genus, levels = c("Nerita", "Cenchritis", "Cerithium"), ordered = TRUE)) %>% # make snail genus and ordered factor
 
   # Generate Figure 5B
  ggplot(aes(x = snail_genus, y = proportion, fill = location)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
    labs(x = "snail genus", y = "proportion of strikes") +
    facet_grid(. ~ location) +
    scale_fill_manual(values = c("#FF5553", "#FFB700", "#4881B2")) +
    theme(legend.position = "none")








#### Table 1 ####
# G Test of Goodness of fit to assess whether mantis shrimp struck regions at the same frequency as would be
# predicted by surface area measurements

GTest <- function(genus){
  genus.g.feedingtrials <- g.feedingtrials %>%
    filter(snail_genus == genus) %>%   # filter by genus
    mutate(G = rep(NA),    # generate columns to hold G statistic and df
           df = rep(NA))
  
  # compile expected proportions of strikes for each region
  predicted <- c(genus.g.feedingtrials$aperture_expected[1], 
                 genus.g.feedingtrials$whorls_expected[1], 
                 genus.g.feedingtrials$apex_expected[1])
  
  # generate G test statistic for every strike sequence
  for(i in 1:nrow(genus.g.feedingtrials)){
    genus.g.feedingtrials$G[i] = G.test(x = c(genus.g.feedingtrials$aperture[i], 
                                              genus.g.feedingtrials$whorls[i], 
                                              genus.g.feedingtrials$apex[i]),
                                        p = predicted)$statistic
    genus.g.feedingtrials$df[i] = G.test(x = c(genus.g.feedingtrials$aperture[i], 
                                              genus.g.feedingtrials$whorls[i], 
                                              genus.g.feedingtrials$apex[i]),
                                        p = predicted)$parameter
  }
  
  # Sum G statistics and df for all strike sequences
  total.g = sum(genus.g.feedingtrials$G)   # sum g statistics
  total.df = sum(genus.g.feedingtrials$df) # sum df
  p =  pchisq(total.g, df = total.df, lower.tail = FALSE)   # calculate p from these totals
  n = nrow(genus.g.feedingtrials)   # calculate number of strike sequences
  print(c(str_c("n = ", n),
          str_c("G = ", total.g), 
          str_c("df = ", total.df), 
          str_c("p = ", p)))
}

# Calculate G Statistic for Nerita spp.
GTest(genus = "Nerita")
# Calculate G statistics for Cenchritis muricatus
GTest(genus = "Cenchritis")
# Calculate G statistic for Cerithium spp.
GTest(genus = "Cerithium")


# Observed proportions of strikes
feedingtrials %>%
  mutate(snail_genus = factor(snail_genus, c("Nerita", "Cenchritis", "Cerithium"), ordered = TRUE)) %>%  # convert genera to ordered factor
  gather(apex_proportion, aperture_proportion, whorls_proportion, key = region, value = proportion) %>%  # gather so each region has its own row
  group_by(snail_genus, region) %>%    # group by genus and strike region
  summarize(             # calculate median, interquartile range (iqr), min, and max
    median = round(median(proportion), digits = 2), 
    low.iqr = round(quantile(proportion)[2], digits = 2),
    hi.iqr = round(quantile(proportion)[4], digits = 2),
    min = round(min(proportion), digits = 2),
    max = round(max(proportion), digits = 2)
  ) %>%
  mutate(region = str_replace(region, "_proportion", ""))

# expected proportions of strikes
feedingtrials %>%
  mutate(snail_genus = factor(snail_genus, c("Nerita", "Cenchritis", "Cerithium"), ordered = TRUE)) %>%   # convert genera to ordered factor
  gather(apex_expected, aperture_expected, whorls_expected, key = region, value = expected) %>%    # gather so each region has its own row
  group_by(snail_genus, region) %>%    # group by genus and strike region
  summarize(    # report expected values
    expected = expected[1]
  ) %>%
  mutate(region = str_replace(region, "_expected", ""))
  











#########################################################################################################
#### Do mantis shrimp target different regions across shell shapes and throughout a strike sequence? ####
#########################################################################################################


#### Table 2 ####
# Compares the proportions of strikes to each region between mantis shrimp handling each snail genus
# Includes only animals that handled all three shell shapes, so that tests can be paired
# A paired Friedman Rank Sum Test is conducted, and if significant, post-hoc paired Wilcoxon Signed
# Rank Tests are conducted to follow

# sample size
length(unique(feedingtrials_paired$mantis_id))

# Calculate tests for aperture
# Paired Friedman Rank Sum Test
friedman.test(aperture_proportion ~ snail_genus | mantis_id, data = feedingtrials_paired)

# spread paired feeding trials data so every mantis shrimp is a row, to facilitiate paired comparisons
posthocDF <- feedingtrials_paired %>%
  select(mantis_id, snail_genus, aperture_proportion) %>%
  spread(key = snail_genus, value = aperture_proportion)
# Post hoc pairwise Wilcoxon Signed Rank Test
pairwise.wilcox.test(x = c(posthocDF$Nerita, posthocDF$Cerithium,
                           posthocDF$Cenchritis),
                     g = rep(c("Nerita", "Cerithium", "Cenchritis"), each = 19),
                     paired = T,
                     p.adj = "none")
# Bonferroni correction:
0.05/3




# Calculate tests for whorls
friedman.test(whorls_proportion ~ snail_genus | mantis_id, data = feedingtrials_paired)



# Calculate tests for apex
friedman.test(apex_proportion ~ snail_genus | mantis_id, data = feedingtrials_paired)

# spread paired feeding trials data so every mantis shrimp is a row, to facilitiate paired comparisons across
# each feeding
posthocDF <- feedingtrials_paired %>%
  select(mantis_id, snail_genus, apex_proportion) %>%
  spread(key = snail_genus, value = apex_proportion)
# Post hoc pairwise Wilcoxon Signed Rank Test
pairwise.wilcox.test(x = c(posthocDF$Nerita, posthocDF$Cerithium,
                           posthocDF$Cenchritis),
                     g = rep(c("Nerita", "Cerithium", "Cenchritis"), each = 19),
                     paired = T,
                     p.adj = "none")
# Bonferroni correction:
0.05/3


rm(posthocDF)








#### Table 3 ####
# Test whether mantis shrimp change where they tend to strike a shell throughout a strike sequence.
# Fit two generalized mixed models: with or without strike region as a fixed effect.
# Include strike sequence as a random effect.
# Compare models with a likelihood ratio test.

# organize data to fit model
strikelocations_change <- strikelocations %>%
  filter(!is.na(strike_location),   # exclude strike that did not hit the shell
         strike_location != "unknown") %>%    # exclude strikes whose location could not be identified
  mutate(strike_location = factor(strike_location, ordered = F, levels = c("apex", "whorls", "aperture")),   # convert strike_location to unordered factor
         snail_genus = factor(snail_genus, ordered = F),   # convert snail genus to unordered factor
         trial_id = factor(trial_id),    # convert trial_id to factor
         time_min = time_sec/60)     # Convert seconds to minutes

# create function to generate and compare 2 models
strategychange <- function(genus){
  # create model with strike location as a fixed effect
  model1 <- glmer(time_min ~ strike_location + (1|trial_id), data = filter(strikelocations_change, snail_genus == genus),
                    family = Gamma(link = "log"))
  # create model without strike location as a fixed effect
  model2 <- glmer(time_min ~ 1 + (1|trial_id), data = filter(strikelocations_change, snail_genus == genus),
                      family = Gamma(link = "log"))
  

  # calculate number of strikes
  print(str_c("# strikes = ",
              nrow(filter(strikelocations_change, snail_genus == genus))))
  
  # calculate number of strike sequences
  print(str_c("# strikes sequences = ",
              length(unique(strikelocations_change$trial_id[strikelocations_change$snail_genus == genus]))))
  
  # Use likelihood ratio test to compare models
  print(anova(model1, model2))
  if(anova(model1, model2)[2,8] < 0.05){
    print(summary(model1))
  }
}

strategychange(genus = "Nerita")
strategychange(genus = "Cenchritis")
strategychange(genus = "Cerithium")















####################################################################################################
#### Do strike locations correspond to regions of the shell that break more easily than others? ####
####################################################################################################


#### Table 4 ####
# Summarize damage to the aperture and apex for snails from each treatment group struck by Ninjabot.
# Use results of mantis shrimp behavioral experiments to calculate how close to being able to be eaten shells are.

# Looking at aperture damage:
ninjabottrials <- ninjabottrials %>%
  mutate(aperture_eat = recode(as.character(aperture_damage),
                               "0" = "none",
                               "1" = "1/3",
                               "2" = "2/3",
                               "3" = "eatable")) %>%
  mutate(aperture_eat = factor(aperture_eat, ordered = TRUE,
                               levels = c("none", "1/3", "2/3", "eatable")))

# Looking at apex damage:
ninjabottrials <- ninjabottrials %>%
  mutate(apex_div_eat = apex_loss_proportion / 0.2,  # divide amount of apex damage by threshold (0.2)
         apex_div_eat = as.character(cut(apex_div_eat, c(-1, 1/6, 3/6, 5/6, 10), labels = FALSE)),
         apex_eat = recode(apex_div_eat,
                           "1" = "none",
                           "2" = "1/3",
                           "3" = "2/3",
                           "4" = "eatable")) %>%
  mutate(apex_eat = factor(apex_eat, ordered = TRUE,
                           levels = c("none", "1/3", "2/3", "eatable")))



# Generate summary table for aperture damage
ninjabottrials %>%
  mutate(treatment = factor(treatment, ordered = TRUE, levels = c("aperture", "whorls", "apex", "random"))) %>%  # convert Ninjabot treatments to ordered factor
  group_by(treatment) %>%
  summarise(
    no_damage = sum(aperture_eat == "none"),
    one_third = sum(aperture_eat == "1/3"),
    two_thirds = sum(aperture_eat == "2/3"),
    eatable = sum(aperture_eat == "eatable")
  )

# Generate summary table for apex damage
ninjabottrials %>%
  mutate(treatment = factor(treatment, ordered = TRUE, levels = c("aperture", "whorls", "apex", "random"))) %>%
  group_by(treatment) %>%
  summarise(
    no_damage = sum(apex_eat == "none"),
    one_third = sum(apex_eat == "1/3"),
    two_thirds = sum(apex_eat == "2/3"),
    eatable = sum(apex_eat == "eatable")
  )









#### Damage caused by Ninjabot in each of the treatment groups ####
# Do the four Ninjabot treatment groups differ in the quantity of damage that they caused?
# Kruskal Wallis Rank Sum Test with a post-hoc Wilcoxon Rank Sum Test with a Bonferroni adjustment

# compile how close to eatable a shell is
ninjabottrials$eat <- rep(NA, nrow(ninjabottrials))
for(i in 1:nrow(ninjabottrials)){
  if (ninjabottrials$treatment[i] == "apex") {  # if the shell was struck at the apex, focus on apex damage
    ninjabottrials$eat[i] <- ninjabottrials$apex_eat[i]
  } else if(ninjabottrials$treatment[i] == "aperture") {  # if the shell was struck on the aperture, focus on aperture damage
    ninjabottrials$eat[i] <- ninjabottrials$aperture_eat[i]
  } else {
    ninjabottrials$eat[i] <- max(c(ninjabottrials$aperture_eat[i], ninjabottrials$apex_eat[i])) # otherwise (whorls and random groups), focus on the most damaged region
  }
}
rm(i)

ninjabottrials$treatment <- factor(ninjabottrials$treatment, ordered = FALSE)   # convert treatment to unordered factor

# Kruskal Wallis Rank Sum Test
kruskal.test(ninjabottrials$eat ~ ninjabottrials$treatment)

# post-hoc Wilcoxon Rank Sum Test
pairwise.wilcox.test(x = as.numeric(ninjabottrials$eat),
                     g = ninjabottrials$treatment,
                     p.adj = "none",
                     paired = FALSE)
# Bonferroni correction
0.05/6










#### Figure 6 ####
# reorganize Ninjabot data so every row is a strike
ninjabotdamage <- ninjabottrials %>%
  filter(!(treatment == "random")) %>%    # exclude shells struck at randome regions
  
  # convert every column which identifies whether damage occured into a single column
  gather(ends_with("spall_location"), key = "strike", value = "location") %>%
  
  # cut unnecessary words in the key column so every strike is just identified by number
  mutate(strike = str_replace(strike, "_spall_location", "")) %>%
  mutate(strike = str_replace(strike, "strike", "")) %>%
  mutate(strike = as.integer(strike),
         treatment = factor(treatment, levels = c("aperture", "whorls", "apex"), ordered = TRUE)) %>% # convert treatments to ordered factor
  
  # tally damage for each strike in each treatment group
  group_by(treatment, strike) %>%    # group data by treatment and strike number
  summarise(total_count = n(),       # calculate the total number of shells in that treatment group
            damage_count = sum(treatment == location, na.rm = TRUE),   # calculate the number of shells damaged on each strike for each treatment
            damage_proportion = damage_count / total_count)            # calculate the proportion of shells damaged on each strike for each treatment


# make figure 6
ggplot(ninjabotdamage, aes(x = strike, y = damage_proportion, color = treatment, shape = treatment)) +
  geom_point() +
  geom_line() +
  labs(x = "strike number",
       y = "proportion of strikes causing damage") +
  theme_classic() +
  scale_color_manual(values = c("#ff5451", "#ffb500", "#477fb0")) +
  scale_shape_manual(values = c(19, 15, 17)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,20)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))
