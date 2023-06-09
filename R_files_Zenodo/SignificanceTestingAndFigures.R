## SENM Significance Testing and Figures
## Ali Galezo
## Created: 14 Dec 2018
## Last modified: 25 Feb 2020

library(dplyr)
library(ggplot2)
library(lazyeval)
library(gridExtra)
library(kableExtra)

## Prepare environment.
options(stringsAsFactors = FALSE)
num_sim <- 1000

## Make folder to save output into.
dir.create("Results")

## Load files.
observed <- read.csv("real_network_metrics_20190930.csv")
expected <- read.csv("all_random_metrics_20190930.csv")
dolphins <- unique(observed$ego)


########################################################################
## Summary stats
########################################################################

observed %>%
  group_by(sex) %>%
  summarize(mean_degree = mean(mixdegree),
            sd_degree = sd(mixdegree),
            n = n(),
            se_degree = sd_degree/sqrt(n),
            min_degree = min(mixdegree),
            max_degree = max(mixdegree))

########################################################################
## Observed vs. expected social network metrics
########################################################################

metrics <- c("mixdegree","ss_degree","os_degree","mixstrength","ss_strength","os_strength","os_strength_kin","os_strength_nonkin","mixcc","percent_close_kin")

combinations <- data.frame(metric = rep(metrics, each = 2),
                           sex = rep(c("FEMALE","MALE"), length(metrics)))

sigTester <- function(metric, sex){
  # This function calculates the p-value according to the rule P = (b + 1)/(m + 1),
  # Where b = number times a permuted value is more extreme than the observed value, and m = number of permutations.
  # It returns the p-value, observed mean, and permuted mean.
  sex_ <- sex
  n <- nrow(expected %>% filter(sex == sex_) %>% distinct(ego))
  expected_means <-
    expected %>%
    group_by(sex, iteration) %>%
    summarize_(mean = interp(~mean(var, na.rm = TRUE), var = as.name(metric))) %>%
    filter(sex == sex_) %>%
    pull(mean)
  observed_mean <- 
    observed %>%
    filter(sex == sex_) %>%
    summarize_(observed_mean = interp(~mean(var, na.rm = TRUE), var = as.name(metric))) %>%
    pull(observed_mean)
  extreme <- ifelse(observed_mean >= median(expected_means),
                    sum(expected_means >= observed_mean),
                    sum(observed_mean >= expected_means))
  p <- (extreme + 1) / (num_sim + 1)
  return(data.frame(metric = metric,
                    sex = sex,
                    n = n,
                    p = p,
                    obs_mean = observed_mean,
                    exp_mean = mean(expected_means),
                    exp_min = min(expected_means),
                    exp_max = max(expected_means)
                    ))
}

results <- do.call("rbind",
                   mapply(sigTester,
                          metric = combinations$metric,
                          sex = combinations$sex,
                          SIMPLIFY = FALSE))
results$adjusted_p <- p.adjust(results$p, "bonferroni")

results %>%
  mutate(obs_mean = round(obs_mean, 2),
         exp_min = round(exp_min, 2),
         exp_max = round(exp_max, 2),
         adjusted_p = round(adjusted_p, 5)) %>%
  mutate(expected_range = paste(exp_min, exp_max, sep = " - ")) %>%
  select(metric, sex, obs_mean, expected_range, adjusted_p) %>%
  mutate(sex = recode(sex, FEMALE = "Female", MALE = "Male")) %>%
  mutate(metric = recode(metric,
                         mixdegree = "Degree",
                         ss_degree = "Same-sex degree",
                         os_degree = "Opposite-sex degree",
                         mixstrength = "Strength",
                         ss_strength = "Same-sex strength",
                         os_strength= "Opposite-sex strength",
                         os_strength_kin = "Opposite-sex strength (kin only)",
                         os_strength_nonkin = "Opposite-sex strength (non-kin only)",
                         mixcc = "Clustering coefficient",
                         percent_close_kin = "Proportion of close kin associated with")) %>%
  rename('Social network metric' = metric,
         'Sex' = sex,
         'Observed mean' = obs_mean,
         'Expected range' = expected_range,
         'p' = adjusted_p) %>%
  kable(format = "html", escape = FALSE) %>%
  kable_styling(font_size = 10) %>%
  save_kable(., file = paste(getwd(), "/Results/SocialMetricTable", sep = ""), self_contained = FALSE)

########################################################################
## Figure 3 plots (observed vs. expected social network metrics)
########################################################################

## Titles
titles <- data.frame(short = metrics,
                     long = c("Number of Associates", "Number of Same-Sex Associates", "Number of Opposite-Sex Associates",
                              "Association Strength", "Same-Sex Strength", "Opposite-Sex Strength", "Opposite-Sex Strength (kin only)", "Opposite-Sex Strength (nonkin only)",
                              "Clustering Coefficient","Proportion of Available Close Kin\\nAssociated With"))

## Axes
axes <- list(ss_degree = list(breaks_x = c(9,12,15,18,21), limits_x = c(8, 22)),
             os_degree = list(breaks_x = c(7,10,13,16,19), limits_x = c(6.9, 19.5)),
             mixstrength = list(breaks_x = seq(0.5, 1.5, 0.25), limits_x = c(0.5,1.5)),
             mixcc = list(breaks_x = seq(0.25,0.45,0.05), limits_x = c(0.23,0.45)),
             percent_close_kin = list(breaks_x = seq(0,0.7,0.2), limits_x = c(-0.05,0.7))
)

## Function to generate plots
SENMplot <- function(metric, sex){
  sex_ <- sex
  expected_vals <-
    expected %>%
    group_by(sex, iteration) %>%
    summarize_(mean = interp(~mean(var, na.rm = TRUE), var = as.name(metric))) %>%
    filter(sex == sex_) %>%
    pull(mean)
  observed_val <-
    observed %>%
    group_by(sex) %>%
    summarize_(observed_mean = interp(~mean(var, na.rm = TRUE), var = as.name(metric))) %>%
    filter(sex == sex_) %>%
    pull(observed_mean)
  observed_se <- 
    observed %>%
    filter(sex == sex_) %>%
    summarize_(n = ~n(),
               observed_sd = interp(~sd(var, na.rm = TRUE), var = as.name(metric))) %>%
    mutate(observed_se = observed_sd / sqrt(n)) %>%
    pull(observed_se)
  plot_data <- data.frame(value = c(expected_vals, observed_val),
                          identity = c(rep(" Expected Mean ", num_sim),
                                       " Observed Mean \\u00B1 SE "))
  ggplot(plot_data) +
    geom_density(aes(x = value, fill = identity)) +
    geom_vline(xintercept = observed_val, color = "black", linetype = "solid", size = 0.8) +
    geom_vline(xintercept = observed_val + observed_se, color = "black", linetype = "dashed", size = 0.8) +
    geom_vline(xintercept = observed_val - observed_se, color = "black", linetype = "dashed", size = 0.8) +
    xlab(titles[titles$short == metric,]$long) +
    ylab("Density") +
    ggtitle(tools::toTitleCase(tolower(sex))) +
    scale_x_continuous(breaks = axes[[metric]]$breaks_x, limits = axes[[metric]]$limits_x) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c("gray", "black")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title=element_blank())
}

## Plot all together
ssdegree <- ggpubr::ggarrange(SENMplot("ss_degree", "FEMALE") + theme(legend.position = "none"),
                              SENMplot("ss_degree", "MALE") + theme(legend.position = "none"),
                              ncol = 2, nrow = 1)
osdegree <- ggpubr::ggarrange(SENMplot("os_degree", "FEMALE") + theme(legend.position = "none"),
                              SENMplot("os_degree", "MALE") + theme(legend.position = "none"),
                              ncol = 2, nrow = 1)
strength <- ggpubr::ggarrange(SENMplot("mixstrength", "FEMALE") + theme(legend.position = "none"),
                              SENMplot("mixstrength", "MALE") + theme(legend.position = "none"),
                              ncol = 2, nrow = 1)
clustering <- ggpubr::ggarrange(SENMplot("mixcc", "FEMALE") + theme(legend.position = "none"),
                                SENMplot("mixcc", "MALE") + theme(legend.position = "none"),
                                ncol = 2, nrow = 1)
propkin <- ggpubr::ggarrange(SENMplot("percent_close_kin", "FEMALE") + theme(legend.position = "none"),
                             SENMplot("percent_close_kin", "MALE") + theme(legend.position = "none"),
                             ncol = 2, nrow = 1)
legend <- ggpubr::get_legend(SENMplot("percent_close_kin", "FEMALE"))

all <- ggpubr::ggarrange(ssdegree,
                         osdegree,
                         strength,
                         clustering,
                         propkin,
                         ncol = 1,
                         labels = c("A","B","C","D","E"),
                         legend.grob = legend,
                         legend = "bottom")
ggsave("SENM_Results.tif", plot = all, device = "tiff", path = "Results", width = 177, height = 238, units = "mm", dpi = 300)

########################################################################
## Compare 1) strength of same-sex bonds
## 2) clustering coefficients in males vs. females
## 3) degree in males vs. females
########################################################################

## Calculate differential between observed and expected (positive value: observed is higher than expected)
combined <- merge(x = expected, y = observed, by = "ego")
combined$ss_strength_differential <- combined$ss_strength.y - combined$ss_strength.x
combined$cc_differential <- combined$mixcc.y - combined$mixcc.x
combined$degree_differential <- combined$mixdegree.y - combined$mixdegree.x
combined <-
  combined %>%
  group_by(ego) %>%
  summarize(strength_differential = mean(ss_strength_differential),
            cc_differential = mean(cc_differential),
            degree_differential = mean(degree_differential),
            sex = max(sex.x))

## Strength of same-sex bonds:
combined %>%
  group_by(sex) %>% 
  summarize(mean = mean(strength_differential, na.rm = TRUE),
            sd = sd(strength_differential, na.rm = TRUE),
            n = n(),
            se = sd/sqrt(n))
observed %>%
  select(ego, sex, ss_strength) %>%
  merge(x = ., y = combined[c("ego","strength_differential")], by = "ego", all.x = TRUE) %>%
  mutate(percent_of_observed = strength_differential/ss_strength) %>%
  group_by(sex) %>%
  summarize(mean = mean(percent_of_observed)*100,
            sd = sd(percent_of_observed)*100,
            n = n(),
            se = sd/sqrt(n))
perm::permTS(combined %>% filter(sex == "FEMALE") %>% pull(strength_differential),
             combined %>% filter(sex == "MALE") %>% pull(strength_differential),
             alternative = "two.sided", exact = TRUE)

## Clustering coefficients:
combined %>%
  group_by(sex) %>% 
  summarize(mean = mean(cc_differential, na.rm = TRUE),
            sd = sd(cc_differential, na.rm = TRUE),
            n = n(),
            se = sd/sqrt(n))
perm::permTS(combined %>% filter(sex == "FEMALE") %>% filter(!is.na(cc_differential)) %>% pull(cc_differential),
             combined %>% filter(sex == "MALE") %>% filter(!is.na(cc_differential)) %>% pull(cc_differential),
             alternative = "two.sided", exact = TRUE)

## Degree:
combined %>%
  group_by(sex) %>% 
  summarize(mean = mean(degree_differential, na.rm = TRUE),
            sd = sd(degree_differential, na.rm = TRUE),
            n = n(),
            se = sd/sqrt(n))
perm::permTS(combined %>% filter(sex == "FEMALE") %>% filter(!is.na(degree_differential)) %>% pull(degree_differential),
             combined %>% filter(sex == "MALE") %>% filter(!is.na(degree_differential)) %>% pull(degree_differential),
             alternative = "two.sided", exact = TRUE)


########################################################################
## Compare strength/abundance of mixed sex vs. opposite sex bonds
########################################################################

## Compare strengh of same-sex vs opposite-sex bonds
strengths <- data.frame(id = as.factor(rep(observed$ego, 2)),
                        sex = as.factor(rep(observed$sex, 2)),
                        class = as.factor(c(rep("opposite", nrow(observed)), rep("same", nrow(observed)))),
                        strength = c(observed$os_strength, observed$ss_strength))
strengths %>%
  group_by(sex, class) %>%
  summarize(mean = mean(strength),
            sd = sd(strength),
            n = n(),
            se = sd/sqrt(n))
coin::symmetry_test(strength ~ class | id, data = strengths, alternative = "two.sided")

## Compare abundance of same-sex vs. opposite-sex bonds
degrees <- data.frame(id = as.factor(rep(observed$ego, 2)),
                      sex = as.factor(rep(observed$sex, 2)),
                      class = as.factor(c(rep("opposite", nrow(observed)), rep("same", nrow(observed)))),
                      degree = c(observed$os_degree, observed$ss_degree))
degrees %>%
  group_by(sex, class) %>%
  summarize(mean = mean(degree),
            sd = sd(degree),
            n = n(),
            se = sd/sqrt(n))
coin::symmetry_test(degree ~ class | id, data = degrees, alternative = "two.sided")
