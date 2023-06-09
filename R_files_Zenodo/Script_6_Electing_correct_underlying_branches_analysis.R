################################################################################ 

# SCRIPT 6: ANALYSING WHETHER PHYLOGENETIC SIGNAL OR MODEL-FIT STATISTICS CAN
# ACCURATELY IDENTIFY THE BRANCH LENGTHS UNDERLYING A BINARY CHARACTER,

# AUTHOR: Jeremy D. Wilson

# ARTICLE: Wilson et al. 2022 - Chronogram or phylogram for ancestral state
# estimation? Model-fit statistics indicate the branch lengths underlying a
# binary character’s evolution

################################################################################

# LOAD PACKAGES-----------------------------------------------------------------

source("Script 4 - Analysis of average error on different branch length sets.R") # Or run this script first...
require(tidyverse)
require(reshape2)
require(gridExtra)


# GENERATE ASSESSMENT COLUMNS, COUNTS AND SIGNIFICANCE TESTS--------------------


data <- mutate(data,
               chosen_tree_lambda = if_else(
                 lambda_chronogram > lambda_phylogram, "chronogram", if_else(
                   lambda_chronogram == lambda_phylogram, "Equal", "phylogram")
               ),
               chosen_tree_lambda_assessment = if_else(
                 character_branch_lengths == "chronogram" & chosen_tree_lambda == "chronogram" |
                   character_branch_lengths == "phylogram" & chosen_tree_lambda == "phylogram",
                 "Correct", if_else(
                   character_branch_lengths == "chronogram" & chosen_tree_lambda == "phylogram" |
                     character_branch_lengths == "phylogram" & chosen_tree_lambda == "chronogram",
                   "Incorrect", "Equal"
                 )
               ),
               chosen_tree_Fritz_D = if_else(
                 Fritz_D_chronogram < Fritz_D_phylogram,
                 "chronogram", "phylogram"
               ),
               chosen_tree_Fritz_D_assessment = if_else(
                 character_branch_lengths == chosen_tree_Fritz_D,
                 "Correct", "Incorrect"
               ),
               chosen_tree_Borges_D = if_else(
                 Borges_D_chrono > Borges_D_phylo,
                 "chronogram", "phylogram"
               ),
               chosen_tree_Borges_D_assessment = if_else(
                 character_branch_lengths == chosen_tree_Borges_D,
                 "Correct", "Incorrect"
               ),
               chosen_tree_AICc = if_else(
                 pmin(chrono_ER_AICc,
                      chrono_ARD_AICc) < pmin(phylo_ER_AICc,
                                              phylo_ARD_AICc),
                 "chronogram", "phylogram"
               ),
               chosen_tree_AICc_assessment = if_else(
                 character_branch_lengths == chosen_tree_AICc,
                 "Correct", "Incorrect"
               ),
               chosen_tree_BIC = if_else(
                 pmin(chrono_ER_BIC,
                      chrono_ARD_BIC) < pmin(phylo_ER_BIC,
                                             phylo_ARD_BIC),
                 "chronogram", "phylogram"
               ),
               chosen_tree_BIC_assessment = if_else(
                 character_branch_lengths == chosen_tree_BIC,
                 "Correct", "Incorrect"
               ),
               HRM_chosen_tree_lambda = if_else(
                 HRM_lambda_chronogram > HRM_lambda_phylogram, "chronogram", if_else(
                   HRM_lambda_chronogram == HRM_lambda_phylogram, "Equal", "phylogram")
               ),
               HRM_chosen_tree_lambda_assessment = if_else(
                 character_branch_lengths == "chronogram" & HRM_chosen_tree_lambda == "chronogram" |
                   character_branch_lengths == "phylogram" & HRM_chosen_tree_lambda == "phylogram",
                 "Correct", if_else(
                   character_branch_lengths == "chronogram" & HRM_chosen_tree_lambda == "phylogram" |
                     character_branch_lengths == "phylogram" & HRM_chosen_tree_lambda == "chronogram",
                   "Incorrect", "Equal"
                 )
               ),
               HRM_chosen_tree_Fritz_D = if_else(
                 HRM_Fritz_D_chronogram < HRM_Fritz_D_phylogram,
                 "chronogram", "phylogram"
               ),
               HRM_chosen_tree_Fritz_D_assessment = if_else(
                 character_branch_lengths == HRM_chosen_tree_Fritz_D,
                 "Correct", "Incorrect"
               ),
               HRM_chosen_tree_Borges_D = if_else(
                 HRM_Borges_D_chrono > HRM_Borges_D_phylo,
                 "chronogram", "phylogram"
               ),
               HRM_chosen_tree_Borges_D_assessment = if_else(
                 character_branch_lengths == HRM_chosen_tree_Borges_D,
                 "Correct", "Incorrect"
               ),
               HRM_chosen_tree_AICc = if_else(
                 pmin(HRM_chrono_ER_AICc,
                      HRM_chrono_ARD_AICc) < pmin(HRM_phylo_ER_AICc,
                                                  HRM_phylo_ARD_AICc),
                 "chronogram", "phylogram"
               ),
               HRM_chosen_tree_AICc_assessment = if_else(
                 character_branch_lengths == HRM_chosen_tree_AICc,
                 "Correct", "Incorrect"
               ),
               HRM_chosen_tree_BIC = if_else(
                 pmin(HRM_chrono_ER_BIC,
                      HRM_chrono_ARD_BIC) < pmin(HRM_phylo_ER_BIC,
                                                 HRM_phylo_ARD_BIC),
                 "chronogram", "phylogram"
               ),
               HRM_chosen_tree_BIC_assessment = if_else(
                 character_branch_lengths == HRM_chosen_tree_BIC,
                 "Correct", "Incorrect"
               ),
               HRMs_chosen_tree_lambda = if_else(
                 HRMs_lambda_chronogram > HRMs_lambda_phylogram, "chronogram", if_else(
                   HRMs_lambda_chronogram == HRMs_lambda_phylogram, "Equal", "phylogram")
               ),
               HRMs_chosen_tree_lambda_assessment = if_else(
                 character_branch_lengths == "chronogram" & HRMs_chosen_tree_lambda == "chronogram" |
                   character_branch_lengths == "phylogram" & HRMs_chosen_tree_lambda == "phylogram",
                 "Correct", if_else(
                   character_branch_lengths == "chronogram" & HRMs_chosen_tree_lambda == "phylogram" |
                     character_branch_lengths == "phylogram" & HRMs_chosen_tree_lambda == "chronogram",
                   "Incorrect", "Equal"
                 )
               ),
               HRMs_chosen_tree_Fritz_D = if_else(
                 HRMs_Fritz_D_chronogram < HRMs_Fritz_D_phylogram,
                 "chronogram", "phylogram"
               ),
               HRMs_chosen_tree_Fritz_D_assessment = if_else(
                 character_branch_lengths == HRMs_chosen_tree_Fritz_D,
                 "Correct", "Incorrect"
               ),
               HRMs_chosen_tree_Borges_D = if_else(
                 HRMs_Borges_D_chrono > HRMs_Borges_D_phylo,
                 "chronogram", "phylogram"
               ),
               HRMs_chosen_tree_Borges_D_assessment = if_else(
                 character_branch_lengths == HRMs_chosen_tree_Borges_D,
                 "Correct", "Incorrect"
               ),
               HRMs_chosen_tree_AICc = if_else(
                 pmin(HRMs_chrono_ER_AICc,
                      HRMs_chrono_ARD_AICc) < pmin(HRMs_phylo_ER_AICc,
                                                   HRMs_phylo_ARD_AICc),
                 "chronogram", "phylogram"
               ),
               HRMs_chosen_tree_AICc_assessment = if_else(
                 character_branch_lengths == HRMs_chosen_tree_AICc,
                 "Correct", "Incorrect"
               ),
               HRMs_chosen_tree_BIC = if_else(
                 pmin(HRMs_chrono_ER_BIC,
                      HRMs_chrono_ARD_BIC) < pmin(HRMs_phylo_ER_BIC,
                                                  HRMs_phylo_ARD_BIC),
                 "chronogram", "phylogram"
               ),
               HRMs_chosen_tree_BIC_assessment = if_else(
                 character_branch_lengths == HRMs_chosen_tree_BIC,
                 "Correct", "Incorrect"
               )
)


metric_accuracy_chi_squared_tests <- 
  map(
    list(
      Markov = list(lambda = data$chosen_tree_lambda,
                    Fritz_D = data$chosen_tree_Fritz_D,
                    Borges_D = data$chosen_tree_Borges_D,
                    AICc = data$chosen_tree_AICc,
                    BIC = data$chosen_tree_BIC),
      HRM = list(lambda = data$HRM_chosen_tree_lambda,
                 Fritz_D = data$HRM_chosen_tree_Fritz_D,
                 Borges_D = data$HRM_chosen_tree_Borges_D,
                 AICc = data$HRM_chosen_tree_AICc,
                 BIC = data$HRM_chosen_tree_BIC),
      HRMs = list(lambda = data$HRMs_chosen_tree_lambda,
                  Fritz_D = data$HRMs_chosen_tree_Fritz_D,
                  Borges_D = data$HRMs_chosen_tree_Borges_D,
                  AICc = data$HRMs_chosen_tree_AICc,
                  BIC = data$HRMs_chosen_tree_BIC)
    ),
    map,   
    function(x){
      chisq.test(table(data$character_branch_lengths, x, exclude = "Equal"))
    }
  )

metric_bias_chi_squared_tests <- 
  map(
    list(
      Markov = list(lambda = "chosen_tree_lambda_assessment",
                    Fritz_D = "chosen_tree_Fritz_D_assessment",
                    Borges_D = "chosen_tree_Borges_D_assessment",
                    AICc = "chosen_tree_AICc_assessment",
                    BIC = "chosen_tree_BIC_assessment"),
      HRM = list(lambda = "HRM_chosen_tree_lambda_assessment",
                 Fritz_D = "HRM_chosen_tree_Fritz_D_assessment",
                 Borges_D = "HRM_chosen_tree_Borges_D_assessment",
                 AICc = "HRM_chosen_tree_AICc_assessment",
                 BIC = "HRM_chosen_tree_BIC_assessment"),
      HRMs = list(lambda = "HRMs_chosen_tree_lambda_assessment",
                  Fritz_D = "HRMs_chosen_tree_Fritz_D_assessment",
                  Borges_D = "HRMs_chosen_tree_Borges_D_assessment",
                  AICc = "HRMs_chosen_tree_AICc_assessment",
                  BIC = "HRMs_chosen_tree_BIC_assessment")
    ),
    map,   
    function(y){
      modified_data <- 
        data %>%
        group_by(character_branch_lengths, .data[[y]]) %>%
        summarise(n = n()) %>%
        spread(character_branch_lengths, n)
      modified_data <- as.data.frame(modified_data)
      rownames(modified_data) <- modified_data[[y]]
      Xsq <-chisq.test(modified_data[-1])
      return(Xsq)
    }
  )


# CHECKING HOW SIMILAR RESULTS FROM AICc AND BIC ARE----------------------------

mean(data$chosen_tree_AICc_assessment == data$chosen_tree_BIC_assessment)
mean(data$HRM_chosen_tree_AICc_assessment == data$HRM_chosen_tree_BIC_assessment)
mean(data$HRMs_chosen_tree_AICc_assessment == data$HRMs_chosen_tree_BIC_assessment)

# COLUMNS TO ASSESS CHANGE IN AICc BETWEEN CORRECT AND INCORRECT BLs------------

data <- 
  data %>%
  mutate(
    AICc_correct = if_else(character_branch_lengths == "chronogram",
                           pmin(chrono_ER_AICc, chrono_ARD_AICc),
                           pmin(phylo_ER_AICc, phylo_ARD_AICc)
    ),
    AICc_incorrect = if_else(character_branch_lengths == "phylogram",
                             pmin(chrono_ER_AICc, chrono_ARD_AICc),
                             pmin(phylo_ER_AICc, phylo_ARD_AICc)
    ),
    AICc_change = AICc_correct - AICc_incorrect,
    HRM_AICc_correct = if_else(character_branch_lengths == "chronogram",
                               pmin(HRM_chrono_ER_AICc, HRM_chrono_ARD_AICc),
                               pmin(HRM_phylo_ER_AICc, HRM_phylo_ARD_AICc)
    ),
    HRM_AICc_incorrect = if_else(character_branch_lengths == "phylogram",
                                 pmin(HRM_chrono_ER_AICc, HRM_chrono_ARD_AICc),
                                 pmin(HRM_phylo_ER_AICc, HRM_phylo_ARD_AICc)
    ),
    HRM_AICc_change = HRM_AICc_correct - HRM_AICc_incorrect,
    HRMs_AICc_correct = if_else(character_branch_lengths == "chronogram",
                                pmin(HRMs_chrono_ER_AICc, HRMs_chrono_ARD_AICc),
                                pmin(HRMs_phylo_ER_AICc, HRMs_phylo_ARD_AICc)
    ),
    HRMs_AICc_incorrect = if_else(character_branch_lengths == "phylogram",
                                  pmin(HRMs_chrono_ER_AICc, HRMs_chrono_ARD_AICc),
                                  pmin(HRMs_phylo_ER_AICc, HRMs_phylo_ARD_AICc)
    ),
    HRMs_AICc_change = HRMs_AICc_correct - HRMs_AICc_incorrect
  )


# TESTING FOR CORRELATION BETWEEN CHANGE IN AICc AND CHANGE IN ERROR------------


correlation_tests_average_error <- list("Markov" = cor.test(data$AICc_change,
                                                            data$error_change,
                                                            method = "spearman"),
                                        "Hidden rates" = cor.test(data$HRM_AICc_change,
                                                                  data$HRM_error_change,
                                                                  method = "spearman"),
                                        "Scaled hidden rates" = cor.test(data$HRMs_AICc_change,
                                                                         data$HRMs_error_change,
                                                                         method = "spearman")
)


# FIGURE COMPARING UTILITY OF DIFFERENT TEST STATISTICS-------------------------

# Lambda.

p_lambda <- data %>%
  melt(measure.vars = c(Markov = "chosen_tree_lambda_assessment",
                        Hidden_Rates = "HRM_chosen_tree_lambda_assessment",
                        Scaled_Hidden_Rates = "HRMs_chosen_tree_lambda_assessment"
  ),
  variable.name = "Character_set"
  ) %>%
  mutate(value = factor(value, levels = c("Equal", "Incorrect", "Correct"))) %>%
  ggplot(aes(x = Character_set)) +
  geom_hline(aes(yintercept=0.5), size=0.5, linetype="dashed", color = "red") +
  geom_bar(mapping = aes(fill = value, colour = value), position = "fill") +
  scale_fill_manual(values = scales::alpha(c("grey50", "#fc8d59", "#99d594"), 0.5)) +
  scale_colour_manual(values = c("grey50", "#fc8d59","#99d594")) +
  scale_x_discrete(labels = c("chosen_tree_lambda_assessment" = "M",
                              "HRM_chosen_tree_lambda_assessment" = "HR",
                              "HRMs_chosen_tree_lambda_assessment" = "SHR")
  ) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtitle("λ") +
  theme(panel.background = element_rect(colour = "black"),
        #panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 16),
        axis.title = element_blank(),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_blank(),
        legend.position = "none"
        #legend.text = element_text(face = "bold", size = 12)
  )

# Fritz D.

p_fritz <- data %>%
  melt(measure.vars = c(Markov = "chosen_tree_Fritz_D_assessment",
                        Hidden_Rates = "HRM_chosen_tree_Fritz_D_assessment",
                        Scaled_Hidden_Rates = "HRMs_chosen_tree_Fritz_D_assessment"
  ),
  variable.name = "Character_set"
  ) %>%
  mutate(value = factor(value, levels = c("Incorrect", "Correct"))) %>%
  ggplot(aes(x = Character_set)) +
  geom_hline(aes(yintercept=0.5), size=0.5, linetype="dashed", color = "red") +
  geom_bar(mapping = aes(fill = value, colour = value), position = "fill") +
  scale_fill_manual(values = scales::alpha(c("#fc8d59", "#99d594"), 0.5)) +
  scale_colour_manual(values = c("#fc8d59","#99d594")) +
  scale_x_discrete(labels = c("chosen_tree_Fritz_D_assessment" = "M",
                              "HRM_chosen_tree_Fritz_D_assessment" = "HR",
                              "HRMs_chosen_tree_Fritz_D_assessment" = "SHR")
  ) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtitle("D") +
  theme(panel.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 16),
        axis.title = element_blank(),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_blank(),
        legend.position = "none"
        #legend.text = element_text(face = "bold", size = 12)
  )

# Borges d.

p_delta <- data %>%
  melt(measure.vars = c(Markov = "chosen_tree_Borges_D_assessment",
                        Hidden_Rates = "HRM_chosen_tree_Borges_D_assessment",
                        Scaled_Hidden_Rates = "HRMs_chosen_tree_Borges_D_assessment"
  ),
  variable.name = "Character_set"
  ) %>%
  mutate(value = factor(value, levels = c("Incorrect", "Correct"))) %>%
  ggplot(aes(x = Character_set)) +
  geom_hline(aes(yintercept=0.5), size=0.5, linetype="dashed", color = "red") +
  geom_bar(mapping = aes(fill = value, colour = value), position = "fill") +
  scale_fill_manual(values = scales::alpha(c("#fc8d59", "#99d594"), 0.5)) +
  scale_colour_manual(values = c("#fc8d59","#99d594")) +
  scale_x_discrete(labels = c("chosen_tree_Borges_D_assessment" = "M",
                              "HRM_chosen_tree_Borges_D_assessment" = "HR",
                              "HRMs_chosen_tree_Borges_D_assessment" = "SHR")
  ) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtitle("δ") +
  theme(panel.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 16),
        axis.title = element_blank(),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(face = "bold", size = 12)
  )

# AICc.

p_AICc <- data %>%
  melt(measure.vars = c(Markov = "chosen_tree_AICc_assessment",
                        Hidden_Rates = "HRM_chosen_tree_AICc_assessment",
                        Scaled_Hidden_Rates = "HRMs_chosen_tree_AICc_assessment"
  ),
  variable.name = "Character_set"
  ) %>%
  mutate(value = factor(value, levels = c("Incorrect", "Correct"))) %>%
  ggplot(aes(x = Character_set)) +
  geom_hline(aes(yintercept=0.5), size=0.5, linetype="dashed", color = "red") +
  geom_bar(mapping = aes(fill = value, colour = value), position = "fill") +
  scale_fill_manual(values = scales::alpha(c("#fc8d59", "#99d594"), 0.5)) +
  scale_colour_manual(values = c("#fc8d59","#99d594")) +
  scale_x_discrete(labels = c("chosen_tree_AICc_assessment" = "M",
                              "HRM_chosen_tree_AICc_assessment" = "HR",
                              "HRMs_chosen_tree_AICc_assessment" = "SHR")
  ) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtitle("AICc") +
  xlab("Character set") +
  labs(fill = "BL selection", colour = "BL selection") +
  theme(panel.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 16),
        axis.title = element_blank(),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_blank(),
        legend.position = "none"
        #legend.text = element_text(face = "bold", size = 12)
  )

# BIC.

p_BIC <- data %>%
  melt(measure.vars = c(Markov = "chosen_tree_BIC_assessment",
                        Hidden_Rates = "HRM_chosen_tree_BIC_assessment",
                        Scaled_Hidden_Rates = "HRMs_chosen_tree_BIC_assessment"
  ),
  variable.name = "Character_set"
  ) %>%
  mutate(value = factor(value, levels = c("Incorrect", "Correct"))) %>%
  ggplot(aes(x = Character_set)) +
  geom_hline(aes(yintercept=0.5), size=0.5, linetype="dashed", color = "red") +
  geom_bar(mapping = aes(fill = value, colour = value), position = "fill") +
  scale_fill_manual(values = scales::alpha(c("#fc8d59", "#99d594"), 0.5)) +
  scale_colour_manual(values = c("#fc8d59","#99d594")) +
  scale_x_discrete(labels = c("chosen_tree_BIC_assessment" = "M",
                              "HRM_chosen_tree_BIC_assessment" = "HR",
                              "HRMs_chosen_tree_BIC_assessment" = "SHR")
  ) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtitle("BIC") +
  xlab("Character set") +
  labs(fill = "BL selection", colour = "BL selection") +
  theme(panel.background = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 16),
        axis.title = element_blank(),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_blank(),
        legend.position = "none"
        #legend.text = element_text(face = "bold", size = 12)
  )


final_plot_2 <- grid.arrange(p_lambda, p_fritz, p_delta, p_AICc, p_BIC, ncol = 5)

ggsave("Figure_4.pdf",
       plot = final_plot_2,
       device = "pdf",
       dpi = 1000,
       width = 12,
       height = 4.51)


# FIGURE LOOKING AT RELATIONSHIP BETWEEN CHANGE IN AICc AND CHANGE IN ERROR-----
# BETWEEN UNDERLYING BRANCH LENGTHS---------------------------------------------

# Markov

p_Markov <- data %>%
  ggplot(., aes(AICc_change, error_change)) +
  geom_point(size = 2, colour = alpha("black", 0.05)) +
  geom_smooth(method=glm, lwd = 0.8, fullrange = F) +
  geom_hline(aes(yintercept=0), size=0.5, linetype="dashed", color = "red") +
  geom_vline(aes(xintercept=0), size=0.5, linetype="dashed", color = "red") +
  scale_y_continuous(breaks = c(-0.05, 0.0, 0.05),
                     limits = c(-0.06, 0.06)) +
  #  scale_x_continuous(breaks = c(-80, -60, -40, -20, 0, 10),
  #                     limits = c(-88, 11)) +
  scale_x_continuous(breaks = c(-80, -60, -40, -20, 0, 20, 40),
                     limits = c(-88, 41)) +
  theme(panel.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(face = "bold", size = 12),
  )

# Hidden Rates.

p_HRM <- data %>%
  ggplot(., aes(HRM_AICc_change, HRM_error_change)) +
  geom_point(size = 2, colour = alpha("black", 0.05))  +
  geom_smooth(method=glm, lwd = 0.8, fullrange = F) +
  geom_hline(aes(yintercept=0), size=0.5, linetype="dashed", color = "red") +
  geom_vline(aes(xintercept=0), size=0.5, linetype="dashed", color = "red") +
  scale_y_continuous(breaks = c(-0.05, 0.0, 0.05),
                     limits = c(-0.06, 0.06)) +
  #scale_x_continuous(breaks = c(-80, -60, -40, -20, 0, 10),
  #                   limits = c(-88, 11)) +
  scale_x_continuous(breaks = c(-80, -60, -40, -20, 0, 20, 40),
                     limits = c(-88, 41)) +
  theme(panel.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(face = "bold", size = 12)
  )


# Scaled Hidden Rates.

p_HRMs <- data %>%
  ggplot(., aes(HRMs_AICc_change, HRMs_error_change)) +
  geom_point(size = 2, colour = alpha("black", 0.05)) +
  geom_smooth(method=glm, lwd = 0.8, fullrange = F) +
  geom_hline(aes(yintercept=0), size=0.5, linetype="dashed", color = "red") +
  geom_vline(aes(xintercept=0), size=0.5, linetype="dashed", color = "red") +
  scale_y_continuous(breaks = c(-0.05, 0.0,  0.05),
                     limits = c(-0.06, 0.06)) +
  scale_x_continuous(breaks = c(-80, -60, -40, -20, 0, 20, 40),
                     limits = c(-88, 41)) +
  theme(panel.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(face = "bold", size = 12),
  )

final_plot_3 <- grid.arrange(p_Markov, p_HRM, p_HRMs, nrow = 3)

ggsave("Figure_S5.pdf",
       plot = final_plot_3,
       device = "pdf",
       dpi = 1000,
       width = 12,
       height = 4.51)