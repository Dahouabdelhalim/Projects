# Load libraries ------
library(tidyverse)
library(patchwork)
# Theme for plots
theme_gsk <- function() {
  theme_minimal()+
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.tag = element_text(face = "bold")
    ) 
}

# Toggle to true if you want script to generate files/figures
# leave it on false if you don't want to override existing files.
write_objects <- FALSE

# Read in the data and clean out empty rows. ------
biomass <- read_csv("data/phase2-harvest.csv")

# Get rid of "empty" rows
biomass <- biomass %>% arrange(focal_species, source_soil, replicate) %>%  
  filter(!(grepl("Empty", biomass$source_soil))) 
# Get a sense for the structure of the dataframe and the data itself

glimpse(biomass)
# 6 columns:
# 1. "number" - this just has a unique number for each pot, and was used to randomize the map
# 1. cont. Not really meaningful for this analysis, except as a way to keep track of each pot.
# 2. pot_id - the unique ID for each pot. soil_focalspecies_replicate
# 3. Source soil - identifies which soil was used to inoculate the pot
# 3. cont. Eight unique soils: AC, FE, HO, SA, PL, UR, and two controls.
# 3. cont. the two controls are "aastr" and "abfld" and refer to sterile and field soil.
# 4. focal_species: identity of the species growing in the pot.
# 5. replicate: which block the pot came from.
# 6. abg_dry_g: dry biomass, in grams. 

# There are 48 combinations of soil x focal species, grown in ten replicate racks
# for a total of 480 rows:
dim(biomass)
# Unfortunately 6 pots did not have live plants at the two-week mark,
# so these plots will be excluded from further analysis: 
biomass %>% filter(is.na(abg_dry_g)) %>% select(pot_id)
num_NA_pots_main <- biomass %>% filter(is.na(abg_dry_g)) %>% nrow
# We can eliminate these pots from the datafame:

biomass <- na.omit(biomass)

# Let's explore the distribution of the data. ------
# Exploring the distribution of biomasses 
biom_hist <- ggplot(biomass) +
  geom_histogram(aes(x = abg_dry_g)) +
  facet_wrap(.~focal_species)
biom_hist_log <- biom_hist + scale_x_log10() + ggtitle("logged x-axis")


biom_hist/biom_hist_log

# Rename columns with more intuitive names - helps down the line
biomass <- biomass %>%
  # Rename soil names; add a column that adds consp or hetsp source
  mutate(pointcol = ifelse(source_soil == "aastr", "Sterile", "Phase 1\\nHeterospecific\\ncultivated"),
         pointcol = ifelse(source_soil == "abfld", "Field", pointcol),
         pointcol = ifelse(str_extract(focal_species, "..")== source_soil, "Phase 1\\nConspecific\\ncultivated", pointcol),
         source_soil = ifelse(source_soil == "abfld", "Field", source_soil),
         source_soil = ifelse(source_soil == "aastr", "Sterile", source_soil))


# Make a "Cleveland dotplot" type plot of the biomass values (Figure 1) -------
# Point will be plotted at the Median, with error bars extending
# from the LQR to the UQR.
# When maximum or minimum values in a group are more then 1.5*IQR
# away from L/UQR, they will be plotted. 
biomass_for_cleveland <- biomass %>%
  group_by(focal_species, source_soil, pointcol) %>%
  summarize(
    # Get the mean and SEM for each (focal*soil source) group:
    mean_bm = mean(abg_dry_g),
    se_bm = sd(abg_dry_g)/sqrt(n()),
    # Also get the median and IQRs
    median_bm = median(abg_dry_g),
    lqr = quantile(abg_dry_g, 0.25),
    uqr = quantile(abg_dry_g, 0.75),
    # Get outliers
    min_val = min(abg_dry_g, na.rm = T),
    max_val = max(abg_dry_g, na.rm = T),
    out_low = ifelse(min_val < lqr - 1.5*(uqr-lqr), min_val, NA),
    out_upp = ifelse(max_val > uqr + 1.5*(uqr-lqr), max_val, NA),
    n = n()) 



# add a column that sets the y-value of the points -- 
# needed, because I want to add the outlier points
biomass_for_cleveland$plot_y <- rev(c(seq(0.6, 1.3, by = .1),
                                      seq(1.6, 2.3, by = .1)+1,
                                      seq(2.6, 3.3, by = .1)+2,
                                      seq(3.6, 4.3, by = .1)+3,
                                      seq(4.6, 5.3, by = .1)+4,
                                      seq(5.6, 6.3, by = .1)+5))
# Make a new data frame that includes the outlier points
# Outliers here are defined as those points that are either
# lower than (LQR - 1.5*IQR), or higher than (UQR + 1.5*IQR)
outliers <- left_join(biomass, biomass_for_cleveland) %>%
  filter(abg_dry_g < lqr - (1.5*(uqr-lqr)) | 
           abg_dry_g > uqr + (1.5*(uqr-lqr))) %>%
  select(focal_species, source_soil, abg_dry_g, plot_y)

(biomass_dotplot <- ggplot(biomass_for_cleveland, 
                           aes(y = plot_y, x = median_bm, fill = pointcol, 
                               label = source_soil, shape = pointcol)) +
    scale_x_log10() +
    # add median value
    geom_point(size = 4, stroke = .9) +
    # add dashed lines connecting outlier to median (if outlier exists)
    geom_errorbarh(aes(xmin = out_low, xmax = lqr),
                   height = 0, linetype = 3, size = .25) +
    geom_errorbarh(aes(xmin = uqr, xmax = out_upp),
                   height = 0, linetype = 3, size = .25) +
    # add outlying points (if they exist)
    geom_point(aes(y = plot_y, x = out_low), shape = 21, fill = "grey50", size = .5) +
    geom_point(aes(y = plot_y, x = out_upp), shape = 21, fill = "grey50", size = .5) +
    geom_point(data = outliers, aes(x = abg_dry_g, y = plot_y),
               shape = 21, fill = "grey50", size = .5, inherit.aes = F) +
    
    # add solid line connecting median to lower and upper quantile
    geom_errorbarh(aes(xmin = lqr, xmax = uqr), height = .1) +
    
    # Set color and shape for all the points
    scale_fill_manual(values = c("#009E73", "grey10", "grey90", "#D55E00"),
                      name = "Inoculum\\nsource") +
    scale_shape_manual(values = c(22,21,21,23),
                       name = "Inoculum\\nsource") +
    # add thin line between each focal species
    geom_hline(yintercept = c(2, 4, 6, 8, 10), linetype = 1, size = 0.25) +
    # add species names along y-axis
    scale_y_continuous(breaks = c(1,3,5,7,9,11),
                       labels = rev(unique(biomass_for_cleveland$focal_species))) +
    # misc. theme settings.
    ylab("Focal species") +
    xlab("Aboveground biomass (g)") +
    theme_gsk() +
    theme(#legend.justification=c(1,0), legend.position=c(1,.15),
      legend.position = "top",
      # legend.background = element_rect(colour = "grey50"),
      legend.text = element_text(size = 8)) + 
    NULL)  

if(write_objects) {
  saveRDS(biomass_dotplot, file = "manuscript/figures/biomass_dotplot.Rds")
}

# Now, the goal is to use the biomass values to estimate
# the degree of microbially mediated stabilization and
# fitness differences.
# Reshaping data for calculating Stabilization and FDs ----

# So first, let's make a column that is simply the log AGB
# of each pot, and then go from there. 
biomass$log_agb <- log(biomass$abg_dry_g)

# Now, we need to convert the "long format" biomass data frame  
# into a "wide format". Here's a function for that:
make_wide_biomass <- function(df) {
  # This takes the biomass dataframe, and converts into a wide
  # dataframe, such that each row is one replicate rack,
  # and each column represents the log AGB of one species in 
  # one soil background in that rack.
  df %>% mutate(source_soil = ifelse(source_soil == "aastr", "str", source_soil),
                source_soil = ifelse(source_soil == "abfld", "fld", source_soil),
                pair = paste0(source_soil, "_", focal_species)) %>%
    select(replicate, pair, log_agb) %>%
    spread(pair, log_agb)
}
biomass_wide <- make_wide_biomass(biomass)

# Calculate the Stabilization between each species pair ----
# Recall that following the definition of the m terms in Bever 1997,
# and the analysis of this model in Kandlikar 2019,
# stablization = -0.5*(log(m1A) + log(m2B) - log(m1B) - log(m2A))
# Here is a function that does this calculation 
# (recall that after the data reshaping above, 
# each column represents a given value of log(m1A))
calculate_stabilization <- function(df) {
  df %>% 
    # In df, each row is one repilcate/rack, and each column
    # represents the growth of one species in one soil type.
    # e.g. "FE_ACWR" is the growth of FE in ACWR-cultivated soil.
    mutate(AC_FE = -0.5*(AC_ACWR - AC_FEMI - FE_ACWR + FE_FEMI),
           AC_HO = -0.5*(AC_ACWR - AC_HOMU - HO_ACWR + HO_HOMU),
           AC_SA = -0.5*(AC_ACWR - AC_SACO - SA_ACWR + SA_SACO),
           AC_PL = -0.5*(AC_ACWR - AC_PLER - PL_ACWR + PL_PLER),
           AC_UR = -0.5*(AC_ACWR - AC_URLI - UR_ACWR + UR_URLI),
           FE_HO = -0.5*(FE_FEMI - FE_HOMU - HO_FEMI + HO_HOMU),
           FE_SA = -0.5*(FE_FEMI - FE_SACO - SA_FEMI + SA_SACO),
           FE_PL = -0.5*(FE_FEMI - FE_PLER - PL_FEMI + PL_PLER),
           FE_UR = -0.5*(FE_FEMI - FE_URLI - UR_FEMI + UR_URLI),
           HO_PL = -0.5*(HO_HOMU - HO_PLER - PL_HOMU + PL_PLER),
           HO_SA = -0.5*(HO_HOMU - HO_SACO - SA_HOMU + SA_SACO),
           HO_UR = -0.5*(HO_HOMU - HO_URLI - UR_HOMU + UR_URLI),
           SA_PL = -0.5*(SA_SACO - SA_PLER - PL_SACO + PL_PLER),
           SA_UR = -0.5*(SA_SACO - SA_URLI - UR_SACO + UR_URLI),
           PL_UR = -0.5*(PL_PLER - PL_URLI - UR_PLER + UR_URLI)) %>%
    select(replicate, AC_FE:PL_UR) %>% 
    gather(pair, stabilization, AC_FE:PL_UR)
}

stabilization_values <- calculate_stabilization(biomass_wide)
# Seven values are NA; we can omit these 
stabilization_values <- stabilization_values %>% 
  filter(!(is.na(stabilization)))

# Calculating the Fitness difference between each species pair ----
# Similarly, we can now calculate the fitness difference between
# each pair. Recall that FD = 0.5*(log(m1A)+log(m1B)-log(m2A)-log(m2B))
# But recall that here, IT IS IMPORTANT THAT
# m1A = (m1_soilA - m1_fieldSoil)!

# The following function does this calculation:
calculate_fitdiffs <- function(df) {
  df %>% 
    mutate(AC_FE = 0.5*((AC_ACWR-Field_ACWR) + (FE_ACWR-Field_ACWR) - (AC_FEMI-Field_FEMI) - (FE_FEMI-Field_FEMI)),
           AC_HO = 0.5*((AC_ACWR-Field_ACWR) + (HO_ACWR-Field_ACWR) - (AC_HOMU-Field_HOMU) - (HO_HOMU-Field_HOMU)),
           AC_SA = 0.5*((AC_ACWR-Field_ACWR) + (SA_ACWR-Field_ACWR) - (AC_SACO-Field_SACO) - (SA_SACO-Field_SACO)),
           AC_PL = 0.5*((AC_ACWR-Field_ACWR) + (PL_ACWR-Field_ACWR) - (AC_PLER-Field_PLER) - (PL_PLER-Field_PLER)),
           AC_UR = 0.5*((AC_ACWR-Field_ACWR) + (UR_ACWR-Field_ACWR) - (AC_URLI-Field_URLI) - (UR_URLI-Field_URLI)),
           FE_HO = 0.5*((FE_FEMI-Field_FEMI) + (HO_FEMI-Field_FEMI) - (FE_HOMU-Field_HOMU) - (HO_HOMU-Field_HOMU)),
           FE_SA = 0.5*((FE_FEMI-Field_FEMI) + (SA_FEMI-Field_FEMI) - (FE_SACO-Field_SACO) - (SA_SACO-Field_SACO)),
           FE_PL = 0.5*((FE_FEMI-Field_FEMI) + (PL_FEMI-Field_FEMI) - (FE_PLER-Field_PLER) - (PL_PLER-Field_PLER)),
           FE_UR = 0.5*((FE_FEMI-Field_FEMI) + (UR_FEMI-Field_FEMI) - (FE_URLI-Field_URLI) - (UR_URLI-Field_URLI)),
           HO_PL = 0.5*((HO_HOMU-Field_HOMU) + (PL_HOMU-Field_HOMU) - (HO_PLER-Field_PLER) - (PL_PLER-Field_PLER)),
           HO_SA = 0.5*((HO_HOMU-Field_HOMU) + (SA_HOMU-Field_HOMU) - (HO_SACO-Field_SACO) - (SA_SACO-Field_SACO)),
           HO_UR = 0.5*((HO_HOMU-Field_HOMU) + (UR_HOMU-Field_HOMU) - (HO_URLI-Field_URLI) - (UR_URLI-Field_URLI)),
           SA_PL = 0.5*((SA_SACO-Field_SACO) + (PL_SACO-Field_SACO) - (SA_PLER-Field_PLER) - (PL_PLER-Field_PLER)),
           SA_UR = 0.5*((SA_SACO-Field_SACO) + (UR_SACO-Field_SACO) - (SA_URLI-Field_URLI) - (UR_URLI-Field_URLI)),
           PL_UR = 0.5*((PL_PLER-Field_PLER) + (UR_PLER-Field_PLER) - (PL_URLI-Field_URLI) - (UR_URLI-Field_URLI))) %>%
    select(replicate, AC_FE:PL_UR) %>% 
    gather(pair, fitdiff_fld, AC_FE:PL_UR)
}
fd_values <- calculate_fitdiffs(biomass_wide)
# Twenty-two values are NA; let's omit these.
fd_values <- fd_values %>% filter(!(is.na(fitdiff_fld)))

# Now, generate statistical summaries of SD and FD

stabiliation_summary <- stabilization_values %>% group_by(pair) %>%
  summarize(mean_sd = mean(stabilization),
            sem_sd = sd(stabilization)/sqrt(n()),
            n_sd = n())

fitdiff_summary <- fd_values %>% group_by(pair) %>%
  summarize(mean_fd = mean(fitdiff_fld),
            sem_fd = sd(fitdiff_fld)/sqrt(n()),
            n_fd = n())

# Combine the two separate data frames.
sd_fd_summary <- left_join(stabiliation_summary, fitdiff_summary)

# Some of the FDs are negative, let's flip these to be positive
# and also flip the label so that the first species in the name
# is always the fitness superior.
sd_fd_summary <- sd_fd_summary %>% 
  # if mean_fd is < 0, the following command gets the absolute
  # value and also flips around the species code so that
  # the fitness superior is always the first species in the code
  mutate(pair = ifelse(mean_fd < 0,
                       paste0(str_extract(pair, "..$"),
                              "_",
                              str_extract(pair, "^..")), 
                       pair),
         mean_fd = abs(mean_fd))

# Make a column that indicates the net outcome (coex. or exclusion)
sd_fd_summary <- sd_fd_summary %>% mutate(
  stabilize = ifelse(mean_sd - 2*sem_sd > 0, "yes", "no"), 
  fitness = ifelse(mean_fd - 2*sem_fd > 0, "yes", "no"), 
  outcome = ifelse(mean_fd - 2*sem_fd >
                     mean_sd + 2*sem_sd, "exclusion", "neutral"),
  outcome2 = ifelse(mean_fd > mean_sd, "exclusion", "coexistence"),
  outcome2 = ifelse(mean_sd < 0, "exclusion or priority effect", outcome2)
)

# Plot SD vs. FD -----------
(sd_fd_xpyplot <- ggplot(sd_fd_summary, aes(x = mean_sd, y = mean_fd, label = pair,
                                            fill = outcome)) +
   geom_ribbon(aes(x = seq(0, 1.55, length.out = 15),
                   ymin = rep(0, 15),
                   ymax = seq(0, 1.55, length.out = 15)),
               fill = "#FFF2CC")  +
   geom_errorbar(aes(ymin = (mean_fd)-2*sem_fd,
                     ymax = (mean_fd)+2*sem_fd),
                 size = .25) +
   geom_errorbarh(aes(xmin = mean_sd-2*sem_sd,
                      xmax = mean_sd+2*sem_sd),
                  size = .25) +
   ggrepel::geom_text_repel(segment.size = 0.025, color = "#610B0B", size = 3,
                            nudge_x = c(.15, -.2, .12, .2, .3,
                                        -0.2, -.25, -.05, .28, -.1,
                                        -.2, 0.15, .20, -.15, -.15),
                            nudge_y = c(.1, .08, .12, -.05, .05,
                                        -.075, -.45, .1, .1, .2,
                                        -.05, -.05, -.25, .08, 0.1)
   ) +
   ylim(c(-.6, 1.5)) +
   xlim(c(-.5, 1.5)) +
   geom_abline(linetype = 2) +
   # geom_label() +
   geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
   geom_point(shape = 21, size = 2, stroke = 1) +
   scale_fill_manual(values = c("grey50", "white"),
                     labels = c("Lower bound of fitness difference\\nestimate is greater than upper\\nbound of stabilization estimate",
                                "Confidence intervals of fitness\\ndifference and stabilization\\nestimates overlap"),
                     name = "") + 
   xlab(bquote(atop("Microbially mediated stabilization",
                    -frac(1,2)~(m["1A"]-m["1B"]-m["2A"]+m["2B"]) ))) +
   ylab(bquote(atop("Microbially mediated fitness difference",
                    frac(1,2)~(m["1A"]+m["1B"]-m["2A"]-m["2B"])))) + 
   annotate("text", x = 1.27, y = .2, label = "coexistence", 
            vjust = 0, hjust = 1, size = 4, fontface = "bold.italic") +
   # annotate("text", x = 1.27, y = .15, label = "(assuming otherwise", 
   #          hjust = 1, size = 2.75, color = "grey45", fontface = "italic") + 
   # annotate("text", x = 1.27, y = .07, label = "equal competitors)", 
   #          hjust = 1, size = 2.75, color = "grey45", fontface = "italic") + 
   
   annotate("text", x = 1.27, y = 1.5, label = "exclusion", 
            vjust = 0, hjust = 1, size = 4, fontface = "bold.italic") + 
   # annotate("text", x = 1.27, y = 1.45, label = "(assuming otherwise", 
   #          hjust = 1, size = 2.75, color = "grey45", fontface = "italic") + 
   # annotate("text", x = 1.27, y = 1.37, label = "equal competitors)", 
   #          hjust = 1, size = 2.75, color = "grey45", fontface = "italic") + 
   
   theme_gsk() +
   theme(axis.line = element_line(size = 0),
         legend.justification=c(1,0), legend.position=c(1, 0.025),
         legend.background = element_rect(colour = "grey50"),
         legend.key.height = unit(1, 'cm'),
         legend.text = element_text(size = 7.5),
         legend.title = element_text(size = 0)) +
   NULL)
if(write_objects) {
  saveRDS(sd_fd_xpyplot, file = "manuscript/figures/sd_fd_xyplot_revised.Rds")
}

# Now, answering some basic questions: ------
# 1. How much lower was the biomass of ACWR in sterile soil than
# average growth with any live soil?
acwr_sterile_mean <- biomass %>% filter(focal_species == "ACWR" & 
                                          source_soil == "Sterile") %>%
  summarize(mean = mean(abg_dry_g)) %>% unlist

acwr_live_mean <- biomass %>% filter(focal_species == "ACWR" & 
                                       source_soil != "Sterile") %>%
  summarize(mean = mean(abg_dry_g)) %>% unlist

acwr_field_mean <- biomass %>% filter(focal_species == "ACWR" & 
                                        source_soil == "Field") %>%
  summarize(mean = mean(abg_dry_g)) %>% unlist
acwr_nonFieldLive_mean <- biomass %>% filter(focal_species == "ACWR" & 
                                               source_soil != "Field" &
                                               source_soil != "Sterile") %>%
  summarize(mean = mean(abg_dry_g)) %>% unlist
acwr_field_nonField_ratio <- acwr_nonFieldLive_mean/acwr_field_mean


nonacwr_field_mean <- biomass %>% filter(focal_species != "ACWR" & 
                                           source_soil == "Field") %>%
  summarize(mean = mean(abg_dry_g)) %>% unlist
nonacwr_nonFieldLive_mean <- biomass %>% filter(focal_species != "ACWR" & 
                                                  source_soil != "Field" &
                                                  source_soil != "Sterile") %>%
  summarize(mean = mean(abg_dry_g)) %>% unlist
nonacwr_field_nonField_ratio <- nonacwr_nonFieldLive_mean/nonacwr_field_mean

# 2. How different is growth in conspecific soils
# vs. growth in heterospecific soils, across all species?
mean_growth_in_conspecific <- biomass %>% filter(pointcol == "Phase 1\\nConspecific\\ncultivated") %>% 
  group_by(focal_species) %>% summarize(mean_self = mean(abg_dry_g)) %>% 
  select(mean_self) %>% unlist
mean_growth_in_heterospecific <- biomass %>% filter(pointcol == "Phase 1\\nHeterospecific\\ncultivated") %>% 
  group_by(focal_species) %>% summarize(mean_nonself = mean(abg_dry_g)) %>%
  select(mean_nonself) %>% unlist

# 3. What is the outcome of running an ANOVA
# for the main experiment?
biomass_lm <- lm(log(abg_dry_g) ~ source_soil*focal_species, biomass)
biomass_lm_anova <- car::Anova(biomass_lm)

# Save these one-off objects as a list to be called in the
# RMarkdown file. 
saveRDS(object = list(biomass_lm_anova = biomass_lm_anova,
                      nonacwr_field_nonField_ratio = nonacwr_field_nonField_ratio,
                      acwr_field_nonField_ratio = acwr_field_nonField_ratio,
                      mean_growth_in_conspecific = mean_growth_in_conspecific,
                      mean_growth_in_heterospecific = mean_growth_in_heterospecific),
        file = "manuscript/figures/relevant_objects.Rds")


# Make a table of outcomes (table 1 of paper) --------
tbl_to_print <- sd_fd_summary %>%
  mutate(`Species pair` = c("\\\\textit{A. wrangelianus/\\nF.  microstachys}",
                            "\\\\textit{A. wrangelianus/\\nH.  murinum}",
                            "\\\\textit{A. wrangelianus/\\nP.  erecta}",
                            "\\\\textit{A. wrangelianus/\\nS. columbariae}",
                            "\\\\textit{A. wrangelianus/\\nU.  lindleyi}",
                            "\\\\textit{H.  murinum/\\nF.  microstachys}",
                            "\\\\textit{P.  erecta/\\nF.  microstachys}",
                            "\\\\textit{F.  microstachys/\\nS. columbariae}",
                            "\\\\textit{U.  lindleyi/\\nF.  microstachys}",
                            "\\\\textit{H.  murinum/\\nP.  erecta}",
                            "\\\\textit{H.  murinum/\\nS. columbariae}",
                            "\\\\textit{H.  murinum/\\nU.  lindleyi}",
                            "\\\\textit{P.  erecta/\\nU.  lindleyi}",
                            "\\\\textit{S. columbariae/\\nP.  erecta}",
                            "\\\\textit{S. columbariae/\\nU.  lindleyi}")) %>%
  mutate_if(is.numeric, round, 3) %>% 
  # mutate(arrangement = c(2, 1, 2, 3, 3,1, 4,4,4,2,1, 3, 4,4,4)) %>%
  mutate(arrangement = c(3,3,3, 2,2,3, 1,2,1, 3,3,2, 2,4,1)) %>%
  arrange(arrangement) %>%
  # mutate(annotation = c(rep("Strong fitness differences drive exclusion", 6),
  #                       rep("Interplay of strong fitness differences and stabilization results in net neutral effects", 3),
  #                       rep("Weak stabilization and fitness differences", 6))) %>%
  mutate(pair = str_replace(pair, "_", "\\\\\\\\_")) %>%
  mutate(outcome2 = ifelse(outcome == "exclusion", paste0("\\\\textbf{", outcome2, "}"), outcome2)) %>%
  mutate(Stabilization = ifelse(mean_sd - 2*sem_sd > 0, paste0("\\\\textbf{",mean_sd,"}"), mean_sd)) %>%
  mutate(`Fitness Difference` = ifelse(mean_fd - 2*sem_fd > 0, paste0("\\\\textbf{",mean_fd,"}"), mean_fd)) %>%
  mutate(Stabilization = paste0(Stabilization,"\\\\quad{\\\\footnotesize(",mean_sd-2*sem_sd, "--", mean_sd+2*sem_sd,")}" ),
         `Fitness Difference` = paste0(`Fitness Difference`," {(\\\\footnotesize",mean_fd-2*sem_fd, "--", mean_sd+2*sem_fd,")}")) %>%
  select(`Species pair`, 
         Code = pair, Stabilization, 
         `Fitness Difference`, 
         `Net effect of PSF` = outcome2) %>%
  # arrange(`Species pair`) %>%
  mutate_all(kableExtra::linebreak, align = "c")  

tbl_to_print <- kableExtra::kable(tbl_to_print,  booktabs = T, escape = F, align = "lcccc", format = "latex",
                                  label = "tab:psf_greenhouse_tab1", 
                                  caption.short = "Microbially mediated stabilization and fitness differences among the fifteen species pairs in our study.",
                                  caption = "Microbially mediated stabilization and fitness differences among the fifteen species pairs in our study. Bold terms in the Stabilization and Fitness Difference columns indicate those values whose confidence intervals do not overlap zero. The net effect of plant-soil feedbacks reflects the relative magnitude of stabilization vs. fitness differences in Bever et al. (1997)'s plant-soil feedback model. ") %>%
  kableExtra::kable_styling(font_size = 10.5)  %>%
  kableExtra::column_spec(1, width = "1.35in") %>%
  kableExtra::column_spec(3, width = ".85in") %>%
  kableExtra::column_spec(4, width = ".85in") %>%
  kableExtra::column_spec(5, width = "1.1in") %>%
  kableExtra::pack_rows("Plant-soil feedbacks tend to promote coexistence (stabilization > fitness difference)\\n", 
                        1,3, indent = F, 
                        latex_gap_space = "1em") %>%
  kableExtra::pack_rows("Plant-soil feedbacks tend to promote exclusion (fitness difference > stabilization)\\n", 
                        4,8,
                        indent = F, latex_gap_space = "1em") %>%
  kableExtra::pack_rows("Strong evidence that plant-soil feedbacks promote exclusion\\n(lower bound fitness difference estimate > upper bound of stabilization estimate)\\n", 
                        9,14,
                        indent = F, latex_gap_space = "1em") %>%
  kableExtra::pack_rows("Plant-soil feedbacks tend to destabilize plant interactions (stabilization < 0)\\n", 
                        15,16,
                        indent = F, latex_gap_space = "1em") %>%
  
  kableExtra::row_spec(c(1,4,6,9,11,13), extra_latex_after = "\\\\rowcolor{gray!10}")

if(write_objects){
  saveRDS(tbl_to_print, "manuscript/figures/outcomes_table.Rds")
  saveRDS(tbl_to_print, "~/grad/format-thesis/manuscript-rmds/psf-greenhouse/figures/outcomes_table.Rds")
  
}


# Comparing ND/FD decomp to feasibility analysis ---------
biomass <- read_csv("data/phase2-harvest.csv")

# Get rid of "empty" rows
biomass <- biomass %>% arrange(focal_species, source_soil, replicate) %>%  
  filter(!(grepl("Empty", biomass$source_soil))) 
biomass$focal_species_2let <- str_extract(biomass$focal_species, "..")
biomass$log_agb <- log(biomass$abg_dry_g)
combinations <- combn(unique(biomass$focal_species_2let),2)
biomass_arranged <- tibble(pair = NA,
                           rep  = NA,
                           B1A = NA, B1B = NA, B1F = NA, B2A = NA, B2B = NA, B2F = NA)


for(current_combination in 1:ncol(combinations)) {
  pair <- combinations[,current_combination]
  sp1 <- pair[1]; sp2 <- pair[2]
  
  # subset the biomass data file
  biomass_sub <- 
    biomass %>% 
    filter(source_soil %in% c(pair, "abfld") & focal_species_2let %in% pair)
  
  bm1_1 <- biomass_sub %>% filter(source_soil == pair[1], focal_species_2let == pair[1]) %>%
    select(log_agb) %>% unlist
  bm1_2 <- biomass_sub %>% filter(source_soil == pair[2], focal_species_2let == pair[1]) %>%
    select(log_agb) %>% unlist
  bm1_f <- biomass_sub %>% filter(source_soil == "abfld", focal_species_2let == pair[1]) %>%
    select(log_agb) %>% unlist
  
  bm2_1 <- biomass_sub %>% filter(source_soil == pair[1], focal_species_2let == pair[2]) %>%
    select(log_agb) %>% unlist
  bm2_2 <- biomass_sub %>% filter(source_soil == pair[2], focal_species_2let == pair[2]) %>%
    select(log_agb) %>% unlist
  bm2_f <- biomass_sub %>% filter(source_soil == "abfld", focal_species_2let == pair[2]) %>%
    select(log_agb) %>% unlist
  
  new_df <- tibble(pair = paste0(pair[1], "_", pair[2]),
                   rep  = unique(biomass_sub$replicate),
                   B1A = bm1_1, B1B = bm1_2, B1F = bm1_f, 
                   B2A = bm2_1, B2B = bm2_2, B2F = bm2_f)
  biomass_arranged <- bind_rows(biomass_arranged, new_df)
}
biomass_arranged <- biomass_arranged %>% filter(!is.na(pair))
predict_each_rep <- biomass_arranged %>% 
  mutate(IS = (B1A-B1F) - (B1B-B1F) - (B2A-B2F) + (B2B-B2F)) %>% 
  mutate(stabilization = -0.5*IS) %>% 
  mutate(fitdiff = 0.5*((B1A-B1F) + (B1B-B1F) - (B2A-B2F) - (B2B-B2F))) %>% 
  mutate(p1_num = ((B2B-B2F) - (B1B-B1F))) %>% 
  mutate(p2_num = ((B1A-B1F) - (B2A-B2F))) %>%
  mutate(check = (p1_num+p2_num - IS < .001)) %>%
  mutate(p1 = p1_num/IS, p2 = p2_num/IS) %>%
  mutate(coex_stfd = ifelse(abs(fitdiff) < stabilization, "coex", "exclude")) %>%
  mutate(coex_IS = ifelse(IS < 0 & p1*p2 > 0, "coex", "exclude"))
# predict_each_rep %>% View
which(predict_each_rep$coex_stfd != predict_each_rep$coex_IS)

table_s1.1 <- predict_each_rep %>% 
  mutate(m1A = B1A-B1F, m1B = B1B-B1F,
         m2A = B2A-B2F, m2B = B2B-B2F,
         pair = str_replace(pair, "_", "\\\\\\\\_")) %>%
  mutate_if(is.numeric, round, 3) %>% 
  select( rep, pair, 
          `$m_{1A}$` = m1A, `$m_{1B}$` = m1B, 
          `$m_{2A}$` = m2A, `$m_{2B}$` = m2B,
          `$I_S$` = IS, 
          `$\\\\hat{p_1}$` = p1, `$\\\\hat{p_2}$` =p2, 
          `outcome \\n (feasibility)` = coex_IS,
          stabilization, 
          `fitness difference` = fitdiff, `outcome (IGR)` = coex_stfd) %>% 
  filter(!is.na(`outcome (IGR)`)) %>%
  mutate(rep = str_remove(rep, "R"),
         rep = as.numeric(rep)) %>%
  arrange(pair, rep) %>%
  mutate_all(kableExtra::linebreak)  %>%
  kableExtra::kable(booktabs = T, longtable = T, label = "tab:psfgreenhouse-tabs31",
                    caption.short = "Coexistence consequences of soil microbes via microbially mediated stabilization and fitness differences vs. via comparing $I_S$ to $p^*$", 
                    caption = "Evaluating the coexistence consequences of soil microbes via microbially mediated stabilization and fitness differences vs. via comparing $I_S$ to $p^*$ yield identical results",
                    linesep = c(rep("",9),  "\\\\addlinespace\\\\addlinespace",
                                rep("",9),  "\\\\addlinespace\\\\addlinespace",
                                rep("",9),  "\\\\addlinespace\\\\addlinespace",
                                rep("",7),  "\\\\addlinespace\\\\addlinespace",
                                rep("",7),  "\\\\addlinespace\\\\addlinespace",
                                rep("",9),  "\\\\addlinespace\\\\addlinespace",
                                rep("",9),  "\\\\addlinespace\\\\addlinespace",
                                rep("",7),  "\\\\addlinespace\\\\addlinespace",
                                rep("",7),  "\\\\addlinespace\\\\addlinespace",
                                rep("",7),  "\\\\addlinespace\\\\addlinespace",
                                rep("",7),  "\\\\addlinespace\\\\addlinespace",
                                rep("",7),  "\\\\addlinespace\\\\addlinespace",
                                rep("",7),  "\\\\addlinespace\\\\addlinespace",
                                rep("",7),  "\\\\addlinespace\\\\addlinespace",
                                rep("",5)),
                    escape = F,format = "latex", align = "cccccc|cccc|ccc") %>%
  kableExtra::column_spec(c(10,12,13), width = ".1in") %>% 
  kableExtra::kable_styling(font_size = 10)  %>%
  kableExtra::row_spec(c(10:19, 30:37, 46:55, 66:73, 82:89, 98:105, 114:121), extra_latex_after = "\\\\rowcolor{gray!10}")





if(write_objects){
  saveRDS(table_s1.1, "manuscript/figures/table_s1.1.Rds")
  saveRDS(table_s1.1, "~/grad/format-thesis/manuscript-rmds/psf-greenhouse/figures/table_s1.1.Rds")
  
}


# Conceptual figure of when coexistence is inferred -----------
# Commented out because this figure was not included in the final ms
# dummy_df <- data.frame(mean_sd = c(.55, .65, 1, .35),
#                        sem_sd = c(.1, .05, .04, .08),
#                        mean_fd = c(.25, .8, .35, 1.3),
#                        sem_fd = c(.1, .05, .075, .04),
#                        outcome = c("neutral", "neutral", "coexist", "exclusion"))
# nsp <- nrow(dummy_df)
# 
# (sd_fd_xy_conceptual <- ggplot(dummy_df, aes(x = mean_sd, y = mean_fd, fill = outcome)) + 
#     geom_ribbon(aes(x = seq(0, 1.55, length.out = 4),
#                     ymin = rep(0, 4),
#                     ymax = seq(0, 1.55, length.out = 4)),
#                 fill = "#FFF2CC")  +
#     geom_errorbar(aes(ymin = (mean_fd)-2*sem_fd,
#                       ymax = (mean_fd)+2*sem_fd),
#                   size = .25, width = .01) +
#     geom_errorbar(aes(ymin = (mean_fd)-2*sem_fd,
#                       ymax = (mean_fd)+2*sem_fd,
#                       x = mean_sd+2*sem_sd),
#                   size = .35, width = 0, linetype = 3) +
#     geom_errorbar(aes(ymin = (mean_fd)-2*sem_fd,
#                       ymax = (mean_fd)+2*sem_fd,
#                       x = mean_sd-2*sem_sd),
#                   size = .35, width = 0, linetype = 3) +
#     
#     geom_errorbarh(aes(xmin = mean_sd-2*sem_sd,
#                        xmax = mean_sd+2*sem_sd),
#                    size = .25, height = .01) +
#     geom_errorbarh(aes(xmin = mean_sd-2*sem_sd,
#                        xmax = mean_sd+2*sem_sd,
#                        y = c(mean_fd+2*sem_fd)),
#                    size = .35, height = 0, linetype = 3) +
#     geom_errorbarh(aes(xmin = mean_sd-2*sem_sd,
#                        xmax = mean_sd+2*sem_sd,
#                        y = c(mean_fd-2*sem_fd)),
#                    size = .35, height = 0, linetype = 3) +
#     
#     xlim(c(-.1, 1.55)) +
#     ylim(c(-.1, 1.55)) +
#     geom_abline(linetype = 2) +
#     # geom_label() +
#     geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
#     geom_point(shape = 21, size = 4, stroke = 1.5) +
#     scale_fill_manual(values = c( "#CC79A7", "grey50","white"),
#                       labels = c("PSFs drive coexistence", 
#                                  "PSFs drive exclusion",
#                                  "No strong evidence that PSFs
# drive exclusion or coexistence"),
#                       name = "Coexistence outcome") + 
#     xlab(bquote(atop("Microbially mediated stabilization",
#                      -frac(1,2)~(m["1A"]-m["1B"]-m["2A"]+m["2B"]) ))) +
#     ylab(bquote(atop("Microbially mediated fitness difference",
#                      frac(1,2)~(m["1A"]+m["1B"]-m["2A"]-m["2B"])))) + 
#     annotate("text", x = 1.1, y = .05, label = "coexistence", 
#              vjust = 0, hjust = 1, size = 4, fontface = "bold.italic") +
#     annotate("text", x = 1.1, y = 1.25, label = "exclusion", 
#              vjust = 0, hjust = 1, size = 4, fontface = "bold.italic") + 
#     
#     theme_gsk() +
#     theme(axis.line = element_line(size = 0),
#           legend.justification=c(1,0), legend.position=c(1,.4),
#           legend.background = element_rect(colour = "grey50"),
#           legend.text = element_text(size = 7)) + 
#     NULL)
# 
# if(write_objects){
#   saveRDS(sd_fd_xy_conceptual, "manuscript/figures/sd_fd_xy_conceptual.Rds")
# }

# Variation experiment -----
var_expt <- read_csv("data/variation-harvest.csv")
# Two plants died and have NA in abg_dry_g, let's eliminate these
num_NA_pots_var <- var_expt %>% filter((is.na(abg_dry_g))) %>% nrow

var_expt <- var_expt %>% filter(!(is.na(abg_dry_g)))

ggplot(var_expt) + 
  geom_boxplot(aes(y = abg_dry_g, x = focal_species))

# Make a dataset of median, lqr, uqr, and outlier
# as we had done for the main experiment
variation_expt_data_sum <- var_expt %>%
  group_by(source_soil, focal_species) %>%
  summarize(median_bm = median(abg_dry_g),
            lqr = quantile(abg_dry_g, 0.25),
            uqr = quantile(abg_dry_g, 0.75),
            min_val = min(abg_dry_g, na.rm = T),
            max_val = max(abg_dry_g, na.rm = T),
            out_low = ifelse(min_val < lqr - 1.5*(uqr-lqr), min_val, NA),
            out_upp = ifelse(max_val > uqr + 1.5*(uqr-lqr), max_val, NA),
            n = n())

# I want to also extract the rows for PL_PLER and PL_FEMI from the 
# dataframe used for the earlier dotplot. 
# Additionally I want to show the averages from the variation experiment
variation_expt_data_splevel <- var_expt %>%
  group_by(focal_species) %>%
  summarize(median_bm = median(abg_dry_g),
            lqr = quantile(abg_dry_g, 0.25),
            uqr = quantile(abg_dry_g, 0.75),
            min_val = min(abg_dry_g, na.rm = T),
            max_val = max(abg_dry_g, na.rm = T),
            out_low = ifelse(min_val < lqr - 1.5*(uqr-lqr), min_val, NA),
            out_upp = ifelse(max_val > uqr + 1.5*(uqr-lqr), max_val, NA),
            n = n()) %>%
  mutate(source_soil = "PL_separate_averaged") %>%
  select(source_soil, everything())


var_expt_for_cleveland <- bind_rows(variation_expt_data_sum,
                                    # The rows for PLER and FEMI in PL soil in the main experiment
                                    # variation_expt_data_splevel,
                                    biomass_for_cleveland %>% filter(focal_species %in% c("PLER", "FEMI"),
                                                                     source_soil == "PL"),
                                    .id = "which_df"
) %>% arrange(focal_species) %>% ungroup %>%
  mutate(plot_y = c(seq(0.8,1.2, length.out = 6),
                    seq(1.8, 2.2, length.out = 6)))

# Make the "Cleveland plot" for the variation experiment -----
(var_expt_cleveland  <- ggplot(var_expt_for_cleveland, aes(x = median_bm, y = plot_y, shape = which_df,
                                                           fill = which_df)) + 
   scale_x_log10() +
   # add median value
   geom_point(size = 4, stroke = .9) +
   scale_shape_manual(name = "Inoculum source",
                      values = c(21, 24), #, 22),
                      labels = c("Single PLER\\nmonoculture", 
                                 #"Mean across single\\nPLER monoculture\\n", )) + 
                                 "Homogenized across\\nPLER monocultures")) +
   scale_fill_manual(name = "Inoculum source",
                     values = c("#56B4E9", "#F0E442"),#, "#CC79A7"),
                     labels = c("Single PLER\\nmonoculture", 
                                # "Mean across single\\nPLER monoculture\\n" )) +
                                "Homogenized across\\nPLER monocultures")) + 
   # add dashed lines connecting outlier to median (if outlier exists)
   geom_errorbarh(aes(xmin = out_low, xmax = lqr),
                  height = 0, linetype = 3, size = .25) +
   geom_errorbarh(aes(xmin = uqr, xmax = out_upp),
                  height = 0, linetype = 3, size = .25) +
   # add outlying points (if they exist)
   geom_point(aes(y = plot_y, x = out_low), shape = 21, fill = "grey50", size = .5) +
   geom_point(aes(y = plot_y, x = out_upp), shape = 21, fill = "grey50", size = .5) +
   # add solid line connecting median to lower and upper quantile
   geom_errorbarh(aes(xmin = lqr, xmax = uqr), height = .01) +
   theme_gsk() +
   theme(legend.justification=c(0,0), legend.position=c(.65,.40),
         # legend.position = "right",
         legend.background = element_rect(colour = "grey50"),
         legend.text = element_text(size = 7),
         legend.title = element_text(size = 8)
   ) +
   xlab("Aboveground biomass (g)") + 
   ylab("Focal species") + 
   scale_y_continuous(breaks = c(1,2),
                      labels = c("FEMI", "PLER")) +
   
   NULL)
if(write_objects) {
  saveRDS(var_expt_cleveland, file = "manuscript/figures/var_expt_cleveland.Rds")
}
# Compare variances across groups in Var Expt ------
# The goal here is to ask whether variance is higher in
# the 50 plants from the variation experiment vs. the 10 PL_PLER
# or PL_FEMI plants. 

var_and_main <- bind_rows(var_expt %>% rename(number = pot_id) %>% select(experiment_id, everything()),
                          biomass %>% rename(pot_tag = pot_id) %>% 
                            filter(focal_species %in% c("PLER", "FEMI"),
                                   source_soil == "PL") %>%
                            select(-pointcol) %>% mutate(experiment_id = "main") %>% 
                            select(experiment_id, everything())) %>%
  mutate(experiment_id = as.factor(experiment_id))
var_and_main_pl <- var_and_main %>% filter(focal_species == "PLER")

summary(aov(log(abg_dry_g)~source_soil, data = var_and_main %>% filter(focal_species == "PLER")))
pl_var_leveneP <- car::leveneTest(abg_dry_g~source_soil, data = var_and_main %>% filter(focal_species == "PLER"))
pl_var_leveneP

summary(aov(log(abg_dry_g)~source_soil, data = var_and_main %>% filter(focal_species == "FEMI")))
fe_var_leveneP <- car::leveneTest(abg_dry_g~source_soil, data = var_and_main %>% filter(focal_species == "FEMI"))
fe_var_leveneP

