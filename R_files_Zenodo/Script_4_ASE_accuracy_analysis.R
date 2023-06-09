################################################################################ 

# SCRIPT 4: ANALYSING WHETHER ASRs ON CORRECT BLs ARE MORE ACCURATE.

# AUTHOR: Jeremy D. Wilson

# ARTICLE: Wilson et al. 2022 - Chronogram or phylogram for ancestral state
# estimation? Model-fit statistics indicate the branch lengths underlying a
# binary characterâ€™s evolution

################################################################################

# LOAD PACKAGES-----------------------------------------------------------------


require(tidyverse)
require(scales)
require(gridExtra)
load("trees&characters&metrics.Rdata")


# IDENTIFY RECONSTRUCTION WITH BEST AICc FOR CHRONO & PHYLO---------------------


best.marginal.recon <- function(AICcs, MarginalReconstruction){
  min <- as.numeric(which.min(AICcs))
  return(MarginalReconstruction[[min]])
}


# Markov characters.

data$chrono_best_marginal_recon <- map2(transpose(data[, c("chrono_ER_AICc",
                                                           "chrono_ARD_AICc")]),
                                        transpose(data[, c("chrono_ER_marginals",
                                                           "chrono_ARD_marginals")])
                                        , best.marginal.recon)


data$phylo_best_marginal_recon <- map2(transpose(data[, c("phylo_ER_AICc",
                                                          "phylo_ARD_AICc")]),
                                       transpose(data[, c("phylo_ER_marginals",
                                                          "phylo_ARD_marginals")])
                                       , best.marginal.recon)


# Hidden Rates characters.

data$HRM_chrono_best_marginal_recon <- map2(transpose(data[, c("HRM_chrono_ER_AICc",
                                                               "HRM_chrono_ARD_AICc")]),
                                            transpose(data[, c("HRM_chrono_ER_marginals",
                                                               "HRM_chrono_ARD_marginals")])
                                            , best.marginal.recon)


data$HRM_phylo_best_marginal_recon <- map2(transpose(data[, c("HRM_phylo_ER_AICc",
                                                              "HRM_phylo_ARD_AICc")]),
                                           transpose(data[, c("HRM_phylo_ER_marginals",
                                                              "HRM_phylo_ARD_marginals")])
                                           , best.marginal.recon)


# Amplified hidden rates characters

data$HRMs_chrono_best_marginal_recon <- map2(transpose(data[, c("HRMs_chrono_ER_AICc",
                                                                "HRMs_chrono_ARD_AICc")]),
                                             transpose(data[, c("HRMs_chrono_ER_marginals",
                                                                "HRMs_chrono_ARD_marginals")])
                                             , best.marginal.recon)


data$HRMs_phylo_best_marginal_recon <- map2(transpose(data[, c("HRMs_phylo_ER_AICc",
                                                               "HRMs_phylo_ARD_AICc")]),
                                            transpose(data[, c("HRMs_phylo_ER_marginals",
                                                               "HRMs_phylo_ARD_marginals")])
                                            , best.marginal.recon)


# CALCULATE AVERAGE NODE ERROR FOR CHRONO & PHYLO-------------------------------


av.ASR.error <- function(TrueNodeStates, EstimatedStates){
  raw_error <- 0
  for(i in seq_along(TrueNodeStates)){
    raw_error <- raw_error + (1 - EstimatedStates[i, TrueNodeStates[i]])
  }
  average_error <- raw_error/length(TrueNodeStates)
  return(average_error)
}


# Markov characters.

data$av_ASR_error_chrono <- map2_dbl(data$node_states,
                                     data$chrono_best_marginal_recon,
                                     av.ASR.error)


data$av_ASR_error_phylo <- map2_dbl(data$node_states,
                                    data$phylo_best_marginal_recon,
                                    av.ASR.error)


# Hidden rates characters.

data$HRM_av_ASR_error_chrono <- map2_dbl(data$HRM_node_states,
                                         data$HRM_chrono_best_marginal_recon,
                                         av.ASR.error)


data$HRM_av_ASR_error_phylo <- map2_dbl(data$HRM_node_states,
                                        data$HRM_phylo_best_marginal_recon,
                                        av.ASR.error)

# Amplified hidden rates characters.

data$HRMs_av_ASR_error_chrono <- map2_dbl(data$HRMs_node_states,
                                          data$HRMs_chrono_best_marginal_recon,
                                          av.ASR.error)


data$HRMs_av_ASR_error_phylo <- map2_dbl(data$HRMs_node_states,
                                         data$HRMs_phylo_best_marginal_recon,
                                         av.ASR.error)


# ASSESS WHETHER AV. NODE ERROR LESS ON CORRECT BL'S---------------------------


# Columns for tests to compare error on correct vs's incorrect BLs, and 
# to look at the decrease in error associated with the choice of the correct BLs.

data <-
  data %>%
  mutate(av_correct_ASR_error = if_else(character_branch_lengths == "chronogram",
                                        av_ASR_error_chrono,
                                        av_ASR_error_phylo),
         av_incorrect_ASR_error = if_else(character_branch_lengths == "chronogram",
                                          av_ASR_error_phylo,
                                          av_ASR_error_chrono),
         error_change = av_correct_ASR_error - av_incorrect_ASR_error,
         HRM_av_correct_ASR_error = if_else(character_branch_lengths == "chronogram",
                                            HRM_av_ASR_error_chrono,
                                            HRM_av_ASR_error_phylo),
         HRM_av_incorrect_ASR_error = if_else(character_branch_lengths == "chronogram",
                                              HRM_av_ASR_error_phylo,
                                              HRM_av_ASR_error_chrono),
         HRM_error_change = HRM_av_correct_ASR_error - HRM_av_incorrect_ASR_error,
         HRMs_av_correct_ASR_error = if_else(character_branch_lengths == "chronogram",
                                             HRMs_av_ASR_error_chrono,
                                             HRMs_av_ASR_error_phylo),
         HRMs_av_incorrect_ASR_error = if_else(character_branch_lengths == "chronogram",
                                               HRMs_av_ASR_error_phylo,
                                               HRMs_av_ASR_error_chrono),
         HRMs_error_change = HRMs_av_correct_ASR_error - HRMs_av_incorrect_ASR_error,
  )


# Paired Wilcoxon: is average node error lower on the correct BLs?

av_error_wilcoxon_test <-
  data %>% 
  {
    list("Markov" = wilcox.test(.$av_correct_ASR_error,
                                .$av_incorrect_ASR_error,
                                alternative = "less",
                                paired = T),
         "HRM" = wilcox.test(.$HRM_av_correct_ASR_error,
                             .$HRM_av_incorrect_ASR_error,
                             alternative = "less",
                             paired = T),
         "HRMs" = wilcox.test(.$HRMs_av_correct_ASR_error,
                              .$HRMs_av_incorrect_ASR_error,
                              alternative = "less",
                              paired = T)
    )
  }

# Generating columns to find counts and proportion of replicates where the
# correct BLs led to more accurate ASRs.

data <-
  data %>%
  mutate(
    best_av_ASR_error = as.factor(if_else(
      av_ASR_error_chrono < av_ASR_error_phylo,
      "chronogram", "phylogram"
    )),
    HRM_best_av_ASR_error = as.factor(if_else(
      HRM_av_ASR_error_chrono < HRM_av_ASR_error_phylo,
      "chronogram", "phylogram"
    )),
    HRMs_best_av_ASR_error = as.factor(if_else(
      HRMs_av_ASR_error_chrono < HRMs_av_ASR_error_phylo,
      "chronogram", "phylogram"
    )),
  ) 

# Proportions.

best_av_ASR_error_props <- purrr::map(list("Markov" = table(data$character_branch_lengths,
                                                            data$best_av_ASR_error),
                                           "Hidden_rates" = table(data$character_branch_lengths,
                                                                  data$HRM_best_av_ASR_error),
                                           "Scaled_hidden_rates" = table(data$character_branch_lengths,
                                                                         data$HRMs_best_av_ASR_error)
),
prop.table)

# Find mean improvement for each character set.

mean_error_change <- data %>%
  dplyr::select(error_change,
                HRM_error_change,
                HRMs_error_change) %>%
  summarise(M_mean = mean(error_change),
            HRM_mean = mean(HRM_error_change),
            HRMs_mean = mean(HRMs_error_change),
  )


# CALCULATE AVERAGE DYNAMIC NODE ERROR FOR CHRONO & PHYLO----------------------


av.dnode.error <- function(TrueNodeStates, EstimatedStates1, EstimatedStates2){
  node_prob_change <- abs(EstimatedStates1[,1] - EstimatedStates2[,1])
  indices <- which(node_prob_change %in% sort(node_prob_change, decreasing = T)[1:5])
  error_per_node <- vector("numeric", length = length(TrueNodeStates))
  for(i in seq_along(TrueNodeStates)){
    error_per_node[i] <- 1 - EstimatedStates1[i, TrueNodeStates[i]]
  }
  x <- mean(error_per_node[indices])
  return(x)
}

# Markov characters.

data$av_dnode_error_chrono <- pmap_dbl(list(data$node_states,
                                            data$chrono_best_marginal_recon,
                                            data$phylo_best_marginal_recon),
                                       av.dnode.error)


data$av_dnode_error_phylo <- pmap_dbl(list(data$node_states,
                                           data$phylo_best_marginal_recon,
                                           data$chrono_best_marginal_recon),
                                      av.dnode.error)


# Hidden rates characters.

data$HRM_av_dnode_error_chrono <- pmap_dbl(list(data$HRM_node_states,
                                                data$HRM_chrono_best_marginal_recon,
                                                data$HRM_phylo_best_marginal_recon),
                                           av.dnode.error)


data$HRM_av_dnode_error_phylo <- pmap_dbl(list(data$HRM_node_states,
                                               data$HRM_phylo_best_marginal_recon,
                                               data$HRM_chrono_best_marginal_recon),
                                          av.dnode.error)


# Amplified hidden rates characters.

data$HRMs_av_dnode_error_chrono <- pmap_dbl(list(data$HRMs_node_states,
                                                 data$HRMs_chrono_best_marginal_recon,
                                                 data$HRMs_phylo_best_marginal_recon),
                                            av.dnode.error)


data$HRMs_av_dnode_error_phylo <- pmap_dbl(list(data$HRMs_node_states,
                                                data$HRMs_phylo_best_marginal_recon,
                                                data$HRMs_chrono_best_marginal_recon),
                                           av.dnode.error)


# ASSESS WHETHER AV. VOLATILE NODE ERROR IS LESS ON CORRECT BL'S----------------


# Columns for tests to compare error on correct vs's incorrect BLs, and 
# to look at the decrease in error associated with the choice of the correct BLs.

data <-
  data %>%
  mutate(av_correct_dnode_error = if_else(character_branch_lengths == "chronogram",
                                          av_dnode_error_chrono,
                                          av_dnode_error_phylo),
         av_incorrect_dnode_error = if_else(character_branch_lengths == "chronogram",
                                            av_dnode_error_phylo,
                                            av_dnode_error_chrono),
         dnode_error_change = av_correct_dnode_error - av_incorrect_dnode_error,
         HRM_av_correct_dnode_error = if_else(character_branch_lengths == "chronogram",
                                              HRM_av_dnode_error_chrono,
                                              HRM_av_dnode_error_phylo),
         HRM_av_incorrect_dnode_error = if_else(character_branch_lengths == "chronogram",
                                                HRM_av_dnode_error_phylo,
                                                HRM_av_dnode_error_chrono),
         HRM_dnode_error_change = HRM_av_correct_dnode_error - HRM_av_incorrect_dnode_error,
         HRMs_av_correct_dnode_error = if_else(character_branch_lengths == "chronogram",
                                               HRMs_av_dnode_error_chrono,
                                               HRMs_av_dnode_error_phylo),
         HRMs_av_incorrect_dnode_error = if_else(character_branch_lengths == "chronogram",
                                                 HRMs_av_dnode_error_phylo,
                                                 HRMs_av_dnode_error_chrono),
         HRMs_dnode_error_change = HRMs_av_correct_dnode_error - HRMs_av_incorrect_dnode_error,
  )


# Paired Wilcoxon: is average dnode error lower on the correct BLs?

true_tree_dnode_av_error_wilcoxon_test <-
  data %>% 
  {
    list("Markov" = wilcox.test(.$av_correct_dnode_error,
                                .$av_incorrect_dnode_error,
                                alternative = "less",
                                paired = T),
         "HRM" = wilcox.test(.$HRM_av_correct_dnode_error,
                             .$HRM_av_incorrect_dnode_error,
                             alternative = "less",
                             paired = T),
         "HRMs" = wilcox.test(.$HRMs_av_correct_dnode_error,
                              .$HRMs_av_incorrect_dnode_error,
                              alternative = "less",
                              paired = T)
    )
  }


# Generating columns to find counts and proportion of replicates where the
# correct BLs led to more accurate dnodes.

data <-
  data %>%
  mutate(
    best_av_dnode_error = as.factor(if_else(
      av_dnode_error_chrono < av_dnode_error_phylo,
      "chronogram", "phylogram"
    )),
    HRM_best_av_dnode_error = as.factor(if_else(
      HRM_av_dnode_error_chrono < HRM_av_dnode_error_phylo,
      "chronogram", "phylogram"
    )),
    HRMs_best_av_dnode_error = as.factor(if_else(
      HRMs_av_dnode_error_chrono < HRMs_av_dnode_error_phylo,
      "chronogram", "phylogram"
    )),
  ) 

# Proportions.

best_av_dnode_error_props <- purrr::map(list("Markov" = table(data$character_branch_lengths,
                                                              data$best_av_dnode_error),
                                             "Hidden_rates" = table(data$character_branch_lengths,
                                                                    data$HRM_best_av_dnode_error),
                                             "Scaled_hidden_rates" = table(data$character_branch_lengths,
                                                                           data$HRMs_best_av_dnode_error)
),
prop.table)

# Find mean improvement for each character set.

mean_dnode_error_change <- data %>%
  dplyr::select(dnode_error_change,
                HRM_dnode_error_change,
                HRMs_dnode_error_change) %>%
  summarise(M_mean = mean(dnode_error_change),
            HRM_mean = mean(HRM_dnode_error_change),
            HRMs_mean = mean(HRMs_dnode_error_change),
  )



# FIGURES-----------------------------------------------------------------------


# Average node error plot.

plot_1 <- data %>%
  dplyr::select(phylogram_type, error_change, HRM_error_change, HRMs_error_change) %>%
  gather(variable, value, -phylogram_type) %>%
  ggplot(aes(x = variable, y = value)) +
  geom_hline(aes(yintercept=0), size=0.5,linetype="dashed", color = "red") +
  geom_jitter(alpha = 0.1, size = 1.5) +
  geom_boxplot(colour = "#FFFF66",lwd = 1, fill = NA, outlier.colour = NA) +
  facet_grid(cols = vars(phylogram_type),
             labeller = labeller(
               phylogram_type = set_names(
                 c("Phylograms: Autocorrelated rates",
                   "Phylograms: Uncorrelated rates"),
                 c("autocorrelated", "uncorrelated")))) +
  scale_x_discrete(labels = c("error_change" = "M",
                              "HRM_error_change" = "HR",
                              "HRMs_error_change" = "AHR")) +
  scale_y_continuous(breaks = c(-0.05, -0.025, 0.0, 0.025, 0.05),
                     limits = c(-0.06, 0.06)) +
  ggtitle("Average node error") +
  theme(panel.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 16),
        axis.title = element_blank(),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_blank())


# Swing node error plot.

plot_2 <- data %>%
  dplyr::select(phylogram_type, dnode_error_change,
                HRM_dnode_error_change, HRMs_dnode_error_change) %>%
  gather(variable, value, -phylogram_type) %>%
  ggplot(aes(x = variable, y = value)) +
  geom_hline(aes(yintercept=0),size=0.5,linetype="dashed",color = "red") +
  geom_jitter(alpha = 0.1, size = 1.5) +
  geom_boxplot(colour = "#FFFF66", lwd = 1,fill = NA, outlier.colour = NA) +
  facet_grid(cols = vars(phylogram_type),
             labeller = labeller(
               phylogram_type = set_names(
                 c("Phylograms: Autocorrelated rates",
                   "Phylograms: Uncorrelated rates"),
                 c("autocorrelated", "uncorrelated")))) +
  scale_x_discrete(labels = c("dnode_error_change" = "M",
                              "HRM_dnode_error_change" = "HR",
                              "HRMs_dnode_error_change" = "AHR")) +
  scale_y_continuous(limits = c(-1, 1)) +
  ggtitle("Swing node error") +
  theme(panel.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 16),
        axis.title = element_blank(),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_blank())


# Final figure.

final_plot <- grid.arrange(plot_1, plot_2, ncol = 2)

ggsave("Figure_2.pdf",
       plot = final_plot,
       device = "pdf",
       dpi = 1000,
       width = 12,
       height = 4.51)