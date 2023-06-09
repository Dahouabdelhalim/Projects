## Activity budgets of surveyed all-juvenile groups
## Ali Galezo

library(dplyr)
library(ggplot2)

## Set seed for reproducible Fisher's exact tests
set.seed(123)

## Load group size and activity data for all-juvenile surveyed groups
juviesurveys <- read.csv("juvenile_surveys_activity.csv")

## Restrict to groups of 2+ individuals
juviesurveys <- juviesurveys %>% filter(group_size > 1)

## Chi square test of independence
## Null: no relationship between group comp (all male/all female/mixed) and activity budget
chi <- data.frame(juviesurveys$group_sex, juviesurveys$first.five.activity)
chi <- table(droplevels(chi))
fisher.test(chi, simulate.p.value = TRUE, B = 1000)

## Significant result --> do post-hoc pairwise tests with Bonferonni corrections
fem.mix <- chi[c("All Female", "Mixed"), ]
mal.fem <- chi[c("All Male", "All Female"), ]
mal.mix <- chi[c("All Male", "Mixed"), ]
mal.fem.p <- fisher.test(mal.fem, B = 1000, simulate.p.value = TRUE)$p.value
fem.mix.p <- fisher.test(fem.mix, B = 1000, simulate.p.value = TRUE)$p.value
mal.mix.p <-fisher.test(mal.mix, B = 1000, simulate.p.value = TRUE)$p.value
mal.fem.p.bonferroni <- mal.fem.p * 3
fem.mix.p.bonferroni <- fem.mix.p * 3
mal.mix.p.bonferroni <- mal.mix.p * 3
data.frame(value = c("Male vs. Female", "Female vs. Mixed", "Male vs. Mixed"),
           adjusted_p = c(round(mal.fem.p.bonferroni, 5),
                          round(fem.mix.p.bonferroni, 5),
                          round(mal.mix.p.bonferroni, 5)))

## Reorder factors for figure
juviesurveys$group_sex <- factor(juviesurveys$group_sex, levels = c("All Female","Mixed","All Male"))

## Figure
juviesurveys %>%
  group_by(group_sex) %>%
  count(first.five.activity) %>%
  mutate(num_samps = sum(n)) %>%
  mutate(prop = n/sum(n),
         se = sqrt(prop * (1-prop) / num_samps)) %>%
  filter(first.five.activity != "Other") %>%
  ggplot(., aes(x = first.five.activity, y = prop, fill = group_sex)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = prop-se, ymax = prop+se), width=.2, position=position_dodge(.9)) +
  theme_classic() +
  xlab("Predominant Activity") + ylab("Proportioon of Groups \\u00B1 SE") +
  guides(fill = guide_legend(title = "Group Type")) +
  scale_y_continuous(limits = c(0,0.55), expand = c(0,0)) +
  scale_fill_grey(start = 0, end = 1) +
  theme(legend.position = c(0.9, 0.85))
ggsave("SurveyActivityBudget.tif", plot = last_plot(), device = "tiff", width = 129, height = 100, units = "mm", dpi = 300)




