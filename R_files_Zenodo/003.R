# Author: Bicheng Dong
# Date: Mon Dec 05 2022
# Last Modified by:
# Last Modified time: Mon Dec 05 2022
# Email: bcdong@bjfu.edu.cn
# Description: R codes for figure s4-s8

# -------------------------------- preparation ------------------------------- #
# cleaning memory
cat("\\014")
rm(list = ls())
gc()

# loading packages
library(AER)
library(broom)
library(ciTools)
library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(purrr)
library(readr)
library(readxl)
library(skimr)
library(stringr)
library(tidyr)
library(VGAM)
library(scales)
library(dplyr)

# set work directory as your path
# warning: setup your own path before running codes !!!
file_path <- "your path"

setwd(file_path)
getwd()

# loading data
invasion.list <- fread("./Data/002.alien.plants_online.csv")
str(invasion.list)

#
invasion.list$invasion_status <- as.factor(invasion.list$invasion_status)
levels(invasion.list$invasion_status)

# ---------------------------------------------------------------------------- #
#                                   table s1                                   #
# ---------------------------------------------------------------------------- #
# -------------------------- testing overdispersion -------------------------- #

fit_lm.lf_form <- lm(Life_form_new ~ invasion_status, data = invasion.list)

fit_lm.res <- lm(scale(log(residence_time_to_2020)) ~ invasion_status, data = invasion.list)

fit_lm.suit <- lm(scale(log(hab_suit_prop_max + 0.001)) ~ invasion_status, data = invasion.list)

summary(fit_lm.lf_form)
summary(fit_lm.res)
summary(fit_lm.suit)

# https://fukamilab.github.io/BIO202/04-B-binary-data.html
# Overdispersion (variance is larger than mean):
# Needs correction when Phi (= D/(n-P)) > 1.5,
# and quick fix is to use family = quasipoisson,
# but cost is that se of parameters will be multiplied by sqrt(Phi).
# Negative binomial GLM for count data, with overdispersion

# ------------------------------ contrast matrix ----------------------------- #
# user-specific contrast matrix
# help form Dr. Duo Chen at Uni. of Konstanz
mat_comp <- cbind(
      c(1 / 3, 1 / 3, 1 / 3),
      c(-1, 1 / 2, 1 / 2),
      c(0, -1, 1)
)

# solution of contrast_matrix
mymat_comp <- solve(t(mat_comp))

# merge contrast matrix with raw data
contrast_comp <- as.data.frame(mymat_comp[, 2:3])

# merge contrast matrix with raw data
contrast_comp <- contrast_comp %>%
      mutate(invasion_status = c("1", "2", "3"))

names(contrast_comp)[1:2] <- c("C1_nonnat_nat", "C2_noninv_inv")
contrast_comp

# merge invasion.list with contrast_comp
invasion.list01 <- invasion.list %>%
      left_join(contrast_comp, by = "invasion_status")

View(invasion.list01)

# ---------------------------------------------------------------------------- #
#                                likelihood test                               #
# ---------------------------------------------------------------------------- #
# use vgam package to cope with zero-truncated data.

# ---------------------------- min. residence time --------------------------- #
# full model
fit_lm.res <- lm(scale(log(residence_time_to_2020)) ~ C1_nonnat_nat + C2_noninv_inv, data = invasion.list01)

# model 1
# drop C2_noninv_inv from full model
fit_lm.res01 <- lm(scale(log(residence_time_to_2020)) ~ C1_nonnat_nat, data = invasion.list01)

# mode 2
# drop C1_nonnat_nat from model 1
fit_lm.res02 <- lm(scale(log(residence_time_to_2020)) ~ 1, data = invasion.list01)

# Error in anova.vglm(fit_vglm.n.province01, fit_vglm.n.province, type = "III",:
# argument 'type' must 'I' or 1 for multiple fits

# Naturalized non-invasive vs. invasive
# anova(model 1 vs. full model)
anova01_lm.res <- anova(fit_lm.res01, fit_lm.res, test = "LRT")
anova01_lm.res

# Non-naturalized vs. naturalized
# anova(model 2 vs. model 1)
anova02_lm.res <- anova(fit_lm.res02, fit_lm.res01, test = "LRT")
anova02_lm.res

# --------------------------- climatic suitability --------------------------- #
# full model
fit_lm.suit <- lm(scale(log(hab_suit_prop_max + 0.001)) ~ C1_nonnat_nat + C2_noninv_inv, data = invasion.list01)

# model 1
# drop C2_noninv_inv from full model
fit_lm.suit01 <- lm(scale(log(hab_suit_prop_max + 0.001)) ~ C1_nonnat_nat, data = invasion.list01)

# mode 2
# drop C1_nonnat_nat from model 1
fit_lm.suit02 <- lm(scale(log(hab_suit_prop_max + 0.001)) ~ 1, data = invasion.list01)

# Naturalized non-invasive vs. invasive
# anova(model 1 vs. full model)
anova01_lm.suit <- anova(fit_lm.suit01, fit_lm.suit, test = "LRT")
anova01_lm.suit

# Non-naturalized vs. naturalized
# anova(model 2 vs. model 1)
anova02_lm.suit <- anova(fit_lm.suit02, fit_lm.suit01, test = "LRT")
anova02_lm.suit

# ------------------------------- anova summary ------------------------------ #
lm_models <- list(
      anova02_lm.res = anova02_lm.res,
      anova01_lm.res = anova01_lm.res,
      anova02_lm.suit = anova02_lm.suit,
      anova01_lm.suit = anova01_lm.suit
)

res_lm_tidy <- map_dfr(lm_models, tidy, .id = "variable")

res_lm_tidy <- res_lm_tidy %>%
      mutate(contrast.name = case_when(
            grepl("anova01", variable) ~ "C1_noninv_inv",
            grepl("anova02", variable) ~ "C1_nonnat_nat"
      ))


View(res_lm_tidy)

# ---------------------------------- output ---------------------------------- #
write.csv(res_lm_tidy, "./Output/Table_s1.csv", row.names = FALSE)

# ---------------------------------------------------------------------------- #
#                                   figure s4                                  #
# ---------------------------------------------------------------------------- #
# ------------------------------- species rank ------------------------------- #
# species summary
summary_sp <- invasion.list01 %>%
      mutate_at(vars(invasion_status), ~ case_when(
            . == 1 ~ "non_nat",
            . == 2 ~ "non_invas",
            . == 3 ~ "invas"
      ))

# the top 20 species
summary_sp01 <- summary_sp %>%
      select(TPL_names, invasion_status, n.company) %>%
      arrange(desc(n.company)) %>%
      head(20)

summary_sp02 <- summary_sp01 %>%
      mutate(TPL_names = factor(TPL_names, levels = rev(TPL_names)))

str(summary_sp02)

# ----------------------------- species rank plot ---------------------------- #
# https://blog.csdn.net/Cccrush/article/details/119671619

# display colors
map_color <- c("#ff8884", "#6f80be", "#a9b8c6")
show_col(map_color)

# ggplot syntax
plot.sp <- ggplot(
      summary_sp02,
      aes(
            fill = invasion_status,
            y = TPL_names,
            x = n.company
      )
) +
      geom_bar(position = "stack", stat = "identity") +
      theme_classic() +
      theme(
            panel.background = element_rect(fill = NA),
            legend.position = "none",
            axis.line = element_line(
                  size = 1,
                  linetype = "solid"
            ),
            axis.ticks = element_line(
                  colour = "black",
                  linetype = "solid",
                  size = 1
            ),
            axis.text = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            axis.text.y = element_text(
                  family = "serif",
                  colour = "black",
                  face = "italic",
                  size = 14
            ),
            axis.title = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            legend.title = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            legend.text = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            )
      ) +
      scale_x_continuous(
            limits = c(0, 200),
            breaks = seq(0, 200, by = 50)
      ) +
      labs(
            x = "Number of online nurseries",
            y = "Taxon"
      ) +
      scale_fill_manual(
            name = "",
            labels = c(
                  "Invasive",
                  "Naturalized non-invasive",
                  "Non-naturalized"
            ),
            values = map_color
      )

plot.sp

# -------------------------------- family rank ------------------------------- #
# family summary
summary_family <- invasion.list01 %>%
      group_by(TPL_Family, invasion_status) %>%
      dplyr::summarize(n.species = length(unique(TPL_names)))

# rename levels of invasion status
summary_family <- summary_family %>%
      mutate_at(vars(invasion_status), ~ case_when(
            . == 1 ~ "non_nat",
            . == 2 ~ "non_invas",
            . == 3 ~ "invas"
      ))

# pivot data from long to wide
summary_family <- summary_family %>%
      tidyr::pivot_wider(
            names_from = invasion_status,
            values_from = n.species
      )

# total number of species across invasion status
summary_family <- summary_family %>%
      rowwise() %>%
      mutate(n.species = sum(c_across(non_nat:invas), na.rm = TRUE)) %>%
      data.frame()

# summary_family[is.na(summary_family)] <- 0
View(summary_family)

# china_data_province$company_number[is.na(china_data_province$company_number)]<-0
summary_family01 <- summary_family %>%
      mutate_at(vars(non_nat:n.species), ~ str_replace_na(., 0)) %>%
      mutate_at(vars(non_nat:n.species), ~ as.numeric(.))

View(summary_family01)

# the top 20 families
summary_family02 <- summary_family01 %>%
      arrange(desc(n.species)) %>%
      head(20)

summary_family02 <- summary_family02 %>%
      mutate(TPL_Family = factor(TPL_Family, levels = rev(TPL_Family)))

View(summary_family02)

# pivot data from wide to long
summary_family03 <- summary_family02 %>%
      select(-n.species) %>%
      tidyr::pivot_longer(
            cols = non_nat:invas,
            names_to = "invasion_status",
            values_to = "n.species"
      )

View(summary_family03)

# ----------------------------- family rank plot ----------------------------- #
# dislay colour
map_color <- c("#ff8884", "#6f80be", "#a9b8c6")
show_col(map_color)

# ggplot syntax
plot.family <- ggplot(
      summary_family03,
      aes(
            fill = invasion_status,
            y = TPL_Family,
            x = n.species
      )
) +
      geom_bar(position = "stack", stat = "identity") +
      theme_classic() +
      theme(
            panel.background = element_rect(fill = NA),
            axis.line = element_line(
                  size = 1,
                  linetype = "solid"
            ),
            axis.ticks = element_line(
                  colour = "black",
                  linetype = "solid", size = 1
            ),
            axis.text = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            axis.title = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            legend.title = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            legend.text = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            legend.position = c(0.6, 0.2)
      ) +
      scale_x_continuous(
            limits = c(0, 140),
            breaks = seq(0, 140, by = 20)
      ) +
      labs(x = "Number of species", y = "Family") +
      scale_fill_manual(
            name = "Invasion status",
            labels = c(
                  "Invasive",
                  "Naturalized non-invasive",
                  "Non-naturalized"
            ),
            values = map_color
      )

(plot.family)

# ---------------------------------- output ---------------------------------- #
plot.rank <- plot.family +
      plot.sp +
      plot_layout(ncol = 2) +
      plot_annotation(tag_levels = "A")

ggexport(plot.rank,
      filename = "./Output/figure_s4.png",
      width = 3200,
      height = 2000,
      pointsize = 10,
      res = 300
)

# ---------------------------------------------------------------------------- #
#                                   figure s8                                  #
# ---------------------------------------------------------------------------- #
# ---------------------------- mean residence time --------------------------- #
# data summary
summary_sp <- invasion.list01 %>%
      mutate_at(vars(invasion_status), ~ case_when(
            . == 1 ~ "non_nat",
            . == 2 ~ "non_invas",
            . == 3 ~ "invas"
      ))

View(summary_sp)

# subset variables
summary_mrt01 <- summary_sp %>%
      select(TPL_names, invasion_status, residence_time_to_2020) %>%
      filter(!is.na(residence_time_to_2020)) %>%
      data.frame()

View(summary_mrt01)

# encoded as a order factor
summary_mrt02 <- summary_mrt01 %>%
      mutate_at(vars(invasion_status), ~ factor(., levels = c("non_nat", "non_invas", "invas")))

levels(summary_mrt02$invasion_status)

str(summary_mrt02)

# ----------------------------------- plot ----------------------------------- #
# user-defined display color
map_color <- c("#a9b8c6", "#6f80be", "#ff8884")

# user-defined contrast segments
plot.scales <- data.frame(
      plot.name = c("res_time", "suit"),
      ymax = c(160, 0.8)
)

plot.scales

# r$> table(summary_mrt02$invasion_status)
#   non_nat non_invas     invas
#       341       154        45

# user-defined x labels
x_axis_mrt <- c(
      "non_nat" = "Non-naturalized\\n(341)",
      "non_invas" = "Naturalized\\nnon-invasiv\\n(154)",
      "invas" = "Invasive\\n(45)"
)

# dodge overlapping objects side-to-side
pd <- position_dodge(width = 0.5)

# ggplot syntax
plot_mrt <- ggplot(
      summary_mrt02,
      aes(
            fill = invasion_status,
            y = residence_time_to_2020,
            x = invasion_status
      )
) +
      geom_boxplot(position = pd, outlier.shape = 1) +
      stat_summary(
            fun = mean, geom = "point",
            shape = 20, size = 4,
            color = "black", fill = "black"
      ) +
      theme_classic() +
      theme(
            panel.background = element_rect(fill = NA),
            legend.position = "none",
            axis.line = element_line(
                  size = 1,
                  linetype = "solid"
            ),
            axis.ticks = element_line(
                  colour = "black",
                  linetype = "solid",
                  size = 1
            ),
            axis.text = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            axis.title = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            legend.title = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            legend.text = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            )
      ) +
      scale_y_continuous(
            limits = c(0, 160),
            breaks = seq(0, 160, by = 20)
      ) +
      labs(
            x = "",
            y = "Min. residence time (year)"
      ) +
      scale_fill_manual(
            name = "",
            labels = c(
                  "Invasive",
                  "Naturalized non-invasive",
                  "Non-naturalized"
            ),
            values = map_color
      ) +
      scale_x_discrete(labels = x_axis_mrt)

(plot_mrt)

# adding significance symbols
plot_mrt01 <- plot_mrt +
      labs(size = 10) +
      geom_segment(aes(
            x = 0.7,
            y = plot.scales$ymax[1],
            xend = 0.7 + 0.6,
            yend = plot.scales$ymax[1] + 0
      )) +
      geom_segment(aes(
            x = 1.7,
            y = plot.scales$ymax[1],
            xend = 2.7 + 0.6,
            yend = plot.scales$ymax[1] + 0
      )) +
      geom_segment(aes(
            x = 1.7,
            y = plot.scales$ymax[1] * 0.9,
            xend = 1.7 + 0.6,
            yend = plot.scales$ymax[1] * 0.9 + 0
      )) +
      geom_segment(aes(
            x = 2.7,
            y = plot.scales$ymax[1] * 0.9,
            xend = 2.7 + 0.6,
            yend = plot.scales$ymax[1] * 0.9 + 0
      )) +
      annotate("text",
            x = 3.45,
            y = plot.scales$ymax[1] * 0.985,
            size = 6,
            label = "***"
      ) +
      annotate("text",
            x = 3.45,
            y = plot.scales$ymax[1] * 0.885,
            size = 6,
            label = "ns"
      )

(plot_mrt01)


# --------------------------- climatic suitability --------------------------- #
# data summary
summary_sp <- invasion.list01 %>%
      mutate_at(vars(invasion_status), ~ case_when(
            . == 1 ~ "non_nat",
            . == 2 ~ "non_invas",
            . == 3 ~ "invas"
      ))

View(summary_sp)

# subset variables
summary_cs01 <- summary_sp %>%
      select(TPL_names, invasion_status, hab_suit_prop_max) %>%
      filter(!is.na(hab_suit_prop_max)) %>%
      data.frame()

# encoded as a order factor
summary_cs02 <- summary_cs01 %>%
      mutate_at(vars(invasion_status), ~ factor(., levels = c("non_nat", "non_invas", "invas")))

levels(summary_cs02$invasion_status)
str(summary_cs02)

# ----------------------------------- plot ----------------------------------- #
# user-defined display color
map_color <- c("#a9b8c6", "#6f80be", "#ff8884")

# user-defined contrast segments
plot.scales <- data.frame(
      plot.name = c("res_time", "suit"),
      ymax = c(160, 0.8)
)

plot.scales

# r$> table(summary_cs02$invasion_status)
#  non_nat non_invas     invas
#       421       151        45

# user-defined x labels
x_axis_cs <- c(
      "non_nat" = "Non-naturalized\\n(421)",
      "non_invas" = "Naturalized\\nnon-invasive\\n(151)",
      "invas" = "Invasive\\n(45)"
)

# dodge overlapping objects side-to-side
pd <- position_dodge(width = 0.5)

# ggplot syntax
plot_cs <- ggplot(
      summary_cs02,
      aes(
            fill = invasion_status,
            y = hab_suit_prop_max,
            x = invasion_status
      )
) +
      geom_boxplot(position = pd, outlier.shape = 1) +
      stat_summary(
            fun = mean, geom = "point",
            shape = 20, size = 4,
            color = "black", fill = "black"
      ) +
      theme_classic() +
      theme(
            panel.background = element_rect(fill = NA),
            legend.position = "none",
            axis.line = element_line(size = 1, linetype = "solid"),
            axis.ticks = element_line(
                  colour = "black",
                  linetype = "solid",
                  size = 1
            ),
            axis.text = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            axis.title = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            legend.title = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            legend.text = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            )
      ) +
      scale_y_continuous(
            limits = c(0, 0.8),
            breaks = seq(0, 0.8, by = 0.2)
      ) +
      labs(x = "", y = "Climatic suitability") +
      scale_fill_manual(
            name = "",
            values = map_color,
            labels = c(
                  "Invasive",
                  "Naturalized non-invasive",
                  "Non-naturalized"
            )
      ) +
      scale_x_discrete(labels = x_axis_cs)

(plot_cs)

plot_cs01 <- plot_cs +
      labs(size = 10) +
      geom_segment(aes(
            x = 0.7,
            y = plot.scales$ymax[2],
            xend = 0.7 + 0.6,
            yend = plot.scales$ymax[2] + 0
      )) +
      geom_segment(aes(
            x = 1.7,
            y = plot.scales$ymax[2],
            xend = 2.7 + 0.6,
            yend = plot.scales$ymax[2] + 0
      )) +
      geom_segment(aes(
            x = 1.7,
            y = plot.scales$ymax[2] * 0.9,
            xend = 1.7 + 0.6,
            yend = plot.scales$ymax[2] * 0.9 + 0
      )) +
      geom_segment(aes(
            x = 2.7,
            y = plot.scales$ymax[2] * 0.9,
            xend = 2.7 + 0.6,
            yend = plot.scales$ymax[2] * 0.9 + 0
      )) +
      annotate("text",
            x = 3.45,
            y = plot.scales$ymax[2] * 0.985,
            size = 6,
            label = "***"
      ) +
      annotate("text",
            x = 3.45,
            y = plot.scales$ymax[2] * 0.885,
            size = 6,
            label = "ns"
      )

(plot_cs01)

# ---------------------------------- output ---------------------------------- #
figure.s8 <- plot_mrt01 +
      plot_cs01 +
      plot_layout(ncol = 1) +
      plot_annotation(tag_levels = "A")

figure.s8

ggexport(figure.s8,
      filename = "./Output/figure_s8.png",
      width = 2000,
      height = 3000,
      pointsize = 10,
      res = 300
)
