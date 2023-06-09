################################################################################ 

# SCRIPT 7: ASSESSING INFLUENCE OF PREDICTORS ON TEST STATISTIC AICC (MARKOV
# CHARACTER SET ONLY).

# AUTHOR: Jeremy D. Wilson

# ARTICLE: Wilson et al. 2022 - Chronogram or phylogram for ancestral state
# estimation? Model-fit statistics indicate the branch lengths underlying a
# binary characterâ€™s evolution

################################################################################

# LOAD PACKAGES-----------------------------------------------------------------

require(tidyverse)
require(mgsub)
require(geiger)
require(FossilSim)
require(phangorn)
require(PerformanceAnalytics)
require(psych)
require(adephylo)
require(gridExtra)


# CALCULATING TREE AND CHARACTER PROPERTIES-------------------------------------

#Tree branch score distance.

tree.distance.after.rescaling <- function(Chrono, Phylo){
  rescaled_chrono <- geiger::rescale(Chrono, "depth", 1)
  phylo_depth <- tree.max(Phylo)
  phylo_mean_patristic <- mean(distRoot(Phylo))
  rescale_value <- phylo_depth/phylo_mean_patristic
  rescaled_phylo <- geiger::rescale(Phylo, "depth", rescale_value)
  kf_distance <- KF.dist(rescaled_phylo, rescaled_chrono, rooted = T)
  return(kf_distance)
}

data$tree_distance <- pmap_dbl(list(data$chronogram,
                                    data$phylogram),
                               tree.distance.after.rescaling)

# Level of assymetry between character transition rates.

data <- mutate(data, rate_difference = as.numeric(
  abs(Q_rate_1/(Q_rate_1 + Q_rate_2) - Q_rate_2/(Q_rate_1 + Q_rate_2))))


# Proportion of tips in rare state. 

data$min_state_prop <- map_dbl(data$tip_states,
                               function(x) {min(prop.table(table(x)))})


# Average branch length (on chrono)

chrono.meanBL.after.rescaling <- function(Chrono){
  rescaled_chrono <- geiger::rescale(Chrono, "depth", 1)
  mean_BL <- mean(rescaled_chrono$edge.length)
  return(mean_BL)
}

data$chrono_mean_BL <- map_dbl(data$chronogram,
                               chrono.meanBL.after.rescaling)

# Branch length homogeneit (on chrono)

sd.p=function(x){sd(x)*sqrt((length(x)-1)/length(x))}

BL.SD.after.rescaling <- function(Chrono){
  rescaled_chrono <- geiger::rescale(Chrono, "depth", 1)
  BL_SD_chrono <- sd.p(rescaled_chrono$edge.length)
  return(BL_SD_chrono)
}

data$chrono_SD_BL <- map_dbl(data$chronogram,
                             BL.SD.after.rescaling)



# CORRELATION TEST--------------------------------------------------------------

correlation_test <-
  data %>%
  transmute("X1" = as.numeric(tree_size),
            "X2" = as.numeric(chronogram_depth),
            "X3" = as.numeric(chrono_mean_BL),
            "X4" = as.numeric(chrono_SD_BL),
            "X5" = as.numeric(rate_difference),
            "X6" = as.numeric(min_state_prop),
            "X7" = as.numeric(tree_distance)
  ) %>% {
    list(
      "corr.test" = corr.test(.,
                              use = "pairwise",
                              method = "spearman",
                              adjust = "none",
                              alpha = .05),
      "chart" = chart.Correlation(.,
                                  method="spearman",
                                  histogram=TRUE,
                                  pch=16)
    )
  }

predictor_list <- c("tree_size",
                    "chronogram_depth",
                    "tree_distance",
                    "chrono_mean_BL",
                    "rate_difference")

data %>%
  transmute(
    chosen_tree_AICc_assessment = as.numeric(mgsub(data$chosen_tree_AICc_assessment,
                                                   c("Incorrect", "Correct"),
                                                   c(0, 1))),
    tree_size = as.numeric(tree_size),
    chronogram_depth = as.numeric(chronogram_depth),
    tree_distance = as.numeric(tree_distance),
    chrono_mean_BL = as.numeric(chrono_mean_BL),
    rate_difference = as.numeric(rate_difference),
  ) %>%
  
  My.stepwise::My.stepwise.glm(data = .,
                               Y = "chosen_tree_AICc_assessment",
                               variable.list = predictor_list,
                               myfamily = "binomial")

# FIGURES-----------------------------------------------------------------------

#tree_size

p_tree_size <- data %>%
  mutate(chosen_tree_AICc_assessment = factor(chosen_tree_AICc_assessment,
                                              levels = c("Incorrect",
                                                         "Correct"))) %>%
  ggplot +
  geom_histogram(mapping = aes(x = tree_size,
                               fill = chosen_tree_AICc_assessment,
                               colour = chosen_tree_AICc_assessment),
                 position = "fill",
                 binwidth = 100,
                 boundary = 0) +
  geom_hline(aes(yintercept=0.5),
             size=0.5,
             linetype="dashed",
             color = "red") +
  scale_fill_manual(values = scales::alpha(c("#fc8d59", "#99d594"), 0.5)) +
  scale_colour_manual(values = c("#fc8d59","#99d594")) +
  scale_x_continuous(breaks = seq(0, 1000, 200)) +
  ggtitle("X1: Tree size") +
  theme(panel.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text = element_text(face = "bold", size = 11),
        axis.title = element_blank(),
        legend.position = "none"
  )

#chronogram depth

p_tree_depth <- data %>%
  mutate(chosen_tree_AICc_assessment = factor(chosen_tree_AICc_assessment,
                                              levels = c("Incorrect",
                                                         "Correct"))) %>%
  ggplot +
  geom_histogram(mapping = aes(x = chronogram_depth,
                               fill = chosen_tree_AICc_assessment,
                               colour = chosen_tree_AICc_assessment),
                 position = "fill",
                 binwidth = 10,
                 boundary = 0) +
  geom_hline(aes(yintercept=0.5),
             size=0.5,
             linetype="dashed",
             color = "red") +
  scale_fill_manual(values = scales::alpha(c("#fc8d59", "#99d594"), 0.5)) +
  scale_colour_manual(values = c("#fc8d59","#99d594")) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  ggtitle("X2: Tree depth") +
  theme(panel.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text = element_text(face = "bold", size = 11),
        axis.title = element_blank(),
        legend.position = "none"
  )


# chronogram mean branch length

p_mean_BL <- data %>%
  mutate(chosen_tree_AICc_assessment = factor(chosen_tree_AICc_assessment,
                                              levels = c("Incorrect",
                                                         "Correct"))) %>%
  ggplot +
  geom_histogram(mapping = aes(x = chrono_mean_BL,
                               fill = chosen_tree_AICc_assessment,
                               colour = chosen_tree_AICc_assessment),
                 position = "fill",
                 breaks=c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12),
                 boundary = 0) +
  geom_hline(aes(yintercept=0.5),
             size=0.5,
             linetype="dashed",
             color = "red") +
  scale_fill_manual(values = scales::alpha(c("#fc8d59", "#99d594"), 0.5)) +
  scale_colour_manual(values = c("#fc8d59","#99d594")) +
  scale_x_continuous(breaks = seq(0, 0.12, 0.04)) +
  ggtitle("X3: Average branch length") +
  theme(panel.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text = element_text(face = "bold", size = 11),
        axis.title = element_blank(),
        legend.position = "none"
  )

#tree distance.

p_tree_distance <- data %>%
  mutate(chosen_tree_AICc_assessment = factor(chosen_tree_AICc_assessment,
                                              levels = c("Incorrect",
                                                         "Correct")),
         tree_distance2 = ifelse(tree_distance > 3.5, 3.6, tree_distance)) %>%
  ggplot +
  geom_histogram(mapping = aes(x = tree_distance2,
                               fill = chosen_tree_AICc_assessment,
                               colour = chosen_tree_AICc_assessment),
                 position = "fill",
                 breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
                 boundary = 0) +
  geom_hline(aes(yintercept=0.5),
             size=0.5,
             linetype="dashed",
             color = "red") +
  scale_fill_manual(values = scales::alpha(c("#fc8d59", "#99d594"), 0.5)) +
  scale_colour_manual(values = c("#fc8d59","#99d594")) +
  scale_x_continuous(breaks = seq(0, 4, 1)) +
  ggtitle("X7: Difference in tree shape") +
  theme(panel.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text = element_text(face = "bold", size = 11),
        axis.title = element_blank(),
        legend.position = "none"
  )

# Final figure.

final_plot <- grid.arrange(p_tree_size, p_tree_depth, p_mean_BL, p_tree_distance, ncol = 4)

ggsave("Fig_4_predictors.pdf",
       plot = final_plot,
       device = "pdf",
       dpi = 1000,
       width = 12,
       height = 4)
