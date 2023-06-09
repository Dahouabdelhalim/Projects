# Author: Bicheng Dong
# Date: Mon Dec 05 2022
# Last Modified by:
# Last Modified time: Mon Dec 05 2022
# Email: bcdong@bjfu.edu.cn
# Description: R codes for figure s5-s7


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
#                                 figures s7                                   #
# ---------------------------------------------------------------------------- #
# ------------------------------- data summary ------------------------------- #
# seeds
# data grouped by life form and invasion status
summary.seed_life.form <- invasion.list %>%
      filter(seed_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new, invasion_status) %>%
      dplyr::summarize(n.species = length(unique(TPL_names))) %>%
      mutate(products = "seeds")

View(summary.seed_life.form)

# data grouped by life form
summary.seed_only.life.form <- invasion.list %>%
      filter(seed_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new) %>%
      dplyr::summarize(n.species = length(unique(TPL_names))) %>%
      mutate(products = "seed")

View(summary.seed_only.life.form)

# data grouped by life form
summary.seed_across.product <- invasion.list %>%
      filter(seed_max == 1 | adult_max == 1 | vegetative_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new) %>%
      dplyr::summarize(n.species = length(unique(TPL_names)))

View(summary.seed_across.product)

# live plants
summary.adult_life.form <- invasion.list %>%
      filter(adult_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new, invasion_status) %>%
      dplyr::summarize(n.species = length(unique(TPL_names))) %>%
      mutate(products = "adult")

View(summary.adult_life.form)

# vegetative organs
summary.vegetative_life.form <- invasion.list %>%
      filter(vegetative_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new, invasion_status) %>%
      dplyr::summarize(n.species = length(unique(TPL_names))) %>%
      mutate(products = "vegetative")

View(summary.vegetative_life.form)

# merge all summarized data
summary.life.form <- rbind(
      summary.seed_life.form,
      summary.adult_life.form,
      summary.vegetative_life.form
)

View(summary.life.form)

# rename levels of invasion status and life form
summary.life.form01 <- summary.life.form %>%
      mutate_at(vars(invasion_status), ~ case_when(
            . == 1 ~ "non_nat",
            . == 2 ~ "non_invas",
            . == 3 ~ "invas"
      )) %>%
      mutate_at(vars(Life_form_new), ~ case_when(
            . == 1 ~ "Short-lived",
            . == 2 ~ "Long-lived",
            . == 3 ~ "Woody"
      ))

View(summary.life.form01)

# life form has missing values, we need to supplement them with zero values
newly.added_data01 <- data.frame(
      Life_form_new = c("Short-lived", "Short-lived", "Short-lived", "Long-lived", "Woody"),
      invasion_status = c("non_nat", "non_invas", "invas", "invas", "invas"),
      n.species = rep(0, 5),
      products = rep("vegetative", 5)
)

newly.added_data01

# 合并两组数据
summary.life.form02 <- rbind(summary.life.form01, newly.added_data01)
View(summary.life.form02)


# data with a clear category of life form grouped by life form
summary.life.form_s1 <- invasion.list %>%
      filter(seed_max == 1 | adult_max == 1 | vegetative_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new, invasion_status) %>%
      dplyr::summarize(total.species = length(unique(TPL_names)))

# rename levels of invasion status and life form
summary.life.form_s2 <- summary.life.form_s1 %>%
      mutate_at(vars(invasion_status), ~ case_when(
            . == 1 ~ "non_nat",
            . == 2 ~ "non_invas",
            . == 3 ~ "invas"
      )) %>%
      mutate_at(vars(Life_form_new), ~ case_when(
            . == 1 ~ "Short-lived",
            . == 2 ~ "Long-lived",
            . == 3 ~ "Woody"
      ))

# merge summary.life.form02 with summary.life.form_s2
summary.life.form03 <- summary.life.form02 %>%
      left_join(summary.life.form_s2,
            by = c("Life_form_new", "invasion_status")
      )

View(summary.life.form03)

# calculating proportion of no. of species at one particular status (one life form x invasion status) within one type of product to no. of species at one particular status (one life form x invasion status)
summary.life.form04 <- summary.life.form03 %>%
      mutate(proportion = n.species / total.species)

View(summary.life.form04)

# encoded as a factor variable
summary.life.form04 <- summary.life.form04 %>%
      mutate_at(
            vars(products),
            ~ factor(., levels = c("seeds", "adult", "vegetative"))
      ) %>%
      mutate_at(
            vars(Life_form_new),
            ~ factor(., levels = c("Short-lived", "Long-lived", "Woody"))
      ) %>%
      mutate_at(
            vars(invasion_status),
            ~ factor(., levels = c("non_nat", "non_invas", "invas"))
      )

levels(summary.life.form04$products)

levels(summary.life.form04$Life_form_new)

levels(summary.life.form04$invasion_status)

# ---------------------------------------------------------------------------- #
#                             figure_s7: facet-plot                            #
# ---------------------------------------------------------------------------- #
# New facet label names for lf variable
lf.labs <- c("Short-lived\\nherbs", "Long-lived\\nherbs", "Woody")
names(lf.labs) <- c("Short-lived", "Long-lived", "Woody")

# New facet label names for invasion_status variable
status.labs <- c("Non-naturalized", "Naturalized\\nnon-invasive", "Invasive")
names(status.labs) <- c("non_nat", "non_invas", "invas")

# user-defined colors
map_color02 <- c("#a9b8c6", "#6f80be", "#ff8884")

# ggplot syntax
plot.life.form <- ggplot(
      summary.life.form04,
      aes(y = proportion, x = products, fill = invasion_status)
) +
      geom_bar(width = .7, position = "dodge", stat = "identity") +
      facet_grid(Life_form_new ~ invasion_status,
            labeller = labeller(
                  Life_form_new = lf.labs,
                  invasion_status = status.labs
            )
      ) +
      theme_classic() +
      theme(
            panel.background = element_rect(fill = NA),
            legend.position = "none",
            strip.text.x = element_text(
                  family = "serif",
                  size = 14
            ),
            strip.text.y = element_text(
                  family = "serif",
                  size = 14
            ),
            strip.background = element_rect(
                  size = 1,
                  linetype = "solid"
            ),
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
            plot.title = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            ),
            axis.text.x = element_text(
                  angle = 90,
                  hjust = 1,
                  vjust = 0.5
            )
      ) +
      scale_fill_manual(
            name = "",
            labels = c(
                  "Invasive",
                  "Naturalized non-invasive",
                  "Non-naturalized"
            ),
            values = map_color02
      ) +
      scale_x_discrete(
            name = "",
            labels = c(
                  "Seeds",
                  "Live plants",
                  "Veg. organs"
            )
      ) +
      scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, by = 0.2)
      ) +
      labs(y = "Proportion of taxa")


(plot.life.form)

# ---------------------------------- output ---------------------------------- #
ggexport(plot.life.form,
      filename = "./Output/figure_s7.jpeg",
      width = 3000,
      height = 3000,
      pointsize = 10,
      res = 300
)

# ---------------------------------------------------------------------------- #
#                                   figure s5a                                 #
# ---------------------------------------------------------------------------- #
# All taxa that are available for online trade are considered here.
# data grouped by life form and invasion status
summary.invasion_x_lf <- invasion.list %>%
      filter(seed_max == 1 | adult_max == 1 | vegetative_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new, invasion_status) %>%
      dplyr::summarize(n.species = length(unique(TPL_names)))

# rename levels of invasion status and life form
summary.invasion_x_lf01 <- summary.invasion_x_lf %>%
      mutate_at(vars(invasion_status), ~ case_when(
            . == 1 ~ "non_nat",
            . == 2 ~ "non_invas",
            . == 3 ~ "invas"
      )) %>%
      mutate_at(vars(Life_form_new), ~ case_when(
            . == 1 ~ "Short-lived",
            . == 2 ~ "Long-lived",
            . == 3 ~ "Woody"
      ))

(summary.invasion_x_lf01)

# data grouped by invasion status
summary.invasion_x_lf_s1 <- invasion.list %>%
      filter(seed_max == 1 | adult_max == 1 | vegetative_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(invasion_status) %>%
      dplyr::summarize(total.species = length(unique(TPL_names)))

# rename levels of invasion status
summary.invasion_x_lf_s2 <- summary.invasion_x_lf_s1 %>%
      mutate_at(vars(invasion_status), ~ case_when(
            . == 1 ~ "non_nat",
            . == 2 ~ "non_invas",
            . == 3 ~ "invas"
      ))

# merge summary.invasion_x_lf01 with summary.invasion_x_lf_s2
summary.invasion_x_lf02 <- summary.invasion_x_lf01 %>%
      left_join(summary.invasion_x_lf_s2, by = "invasion_status")

summary.invasion_x_lf02

# calculating proportion of no. of species at one particular status (one life form x invasion status) within one type of product to no. of species at one particular status (one life form x invasion status)
summary.invasion_x_lf03 <- summary.invasion_x_lf02 %>%
      mutate(proportion = n.species / total.species)

summary.invasion_x_lf03

# encoded as a factor variable
summary.invasion_x_lf03 <- summary.invasion_x_lf03 %>%
      mutate_at(
            vars(Life_form_new),
            ~ factor(., levels = c("Short-lived", "Long-lived", "Woody"))
      ) %>%
      mutate_at(
            vars(invasion_status),
            ~ factor(., levels = c("non_nat", "non_invas", "invas"))
      )

View(summary.invasion_x_lf03)

# ----------------------------------- plot ----------------------------------- #
# label names of lf variable
lf.labs <- c("Short-lived\\nherbs", "Long-lived\\nherbs", "Woody")
names(lf.labs) <- c("Short-lived", "Long-lived", "Woody")

# label names of status variable
status.labs <- c("Non-naturalized", "Naturalized\\nnon-invasive", "Invasive")
names(status.labs) <- c("non_nat", "non_invas", "invas")

# user-defined colour
map_color03 <- c("#7bc0f7", "#ffdb69", "#f18226")

# ggplot syntax
plot.invasion_x_lf <- ggplot(
      summary.invasion_x_lf03,
      aes(
            x = invasion_status,
            y = proportion,
            fill = Life_form_new
      )
) +
      geom_bar(stat = "identity") +
      theme_classic() +
      theme(
            panel.background = element_rect(fill = NA),
            legend.position = "none",
            strip.text.x = element_text(family = "serif", size = 14),
            strip.text.y = element_text(family = "serif", size = 14),
            strip.background = element_rect(size = 1, linetype = "solid"),
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
            plot.title = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            )
      ) +
      scale_fill_manual(
            name = "",
            labels = c(
                  "Short-lived herbs",
                  "Long-lived herbs",
                  "Woody species"
            ),
            values = map_color03
      ) +
      scale_x_discrete(
            name = "",
            labels = c(
                  "Non-naturalized\\n(673)",
                  "Naturalized\\nnon-invasive (152)",
                  "Invasive\\n(44)"
            )
      ) +
      scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, by = 0.2)
      ) +
      labs(y = "Proportion of taxa")


plot.invasion_x_lf

# ---------------------------------------------------------------------------- #
#                                   figure s5b                                 #
# ---------------------------------------------------------------------------- #
# All taxa that are available for online trade are considered here.
# data grouped by life form
summary.seed_x_life.form <- invasion.list %>%
      filter(seed_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new) %>%
      dplyr::summarize(n.species = length(unique(TPL_names))) %>%
      mutate(products = "seeds")

# live plants
summary.adult_x_life.form <- invasion.list %>%
      filter(adult_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new) %>%
      dplyr::summarize(n.species = length(unique(TPL_names))) %>%
      mutate(products = "adult")

# veg. organs
summary.vegetative_x_life.form <- invasion.list %>%
      filter(vegetative_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new) %>%
      dplyr::summarize(n.species = length(unique(TPL_names))) %>%
      mutate(products = "vegetative")

# merge three subsets
summary.product_x_life.form <- rbind(
      summary.seed_x_life.form,
      summary.adult_x_life.form,
      summary.vegetative_x_life.form
)

(summary.product_x_life.form)

# renames levels of life form
summary.product_x_life.form01 <- summary.product_x_life.form %>%
      mutate_at(vars(Life_form_new), ~ case_when(
            . == 1 ~ "Short-lived",
            . == 2 ~ "Long-lived",
            . == 3 ~ "Woody"
      ))

# supplement zero values of short-lived species
new.row01 <- data.frame(
      Life_form_new = c("Short-lived"),
      n.species = rep(0, 1),
      products = rep("vegetative", 1)
)

(new.row01)

# merge summary.product_x_life.form01 and new.row01
summary.product_x_life.form02 <- rbind(summary.product_x_life.form01, new.row01)
View(summary.product_x_life.form02)

# data grouped by products
summary.product_x_life.form_s1 <- summary.product_x_life.form02 %>%
      group_by(products) %>%
      dplyr::summarize(total.species = sum(n.species))

# merge summary.product_x_life.form01 with summary.product_x_life.form_s1
summary.product_x_life.form02 <- summary.product_x_life.form01 %>% left_join(summary.product_x_life.form_s1, by = "products")

summary.product_x_life.form02

# calculating proportion of no. of species at one particular status (one life form x product type) within one type of product to no. of species at one particular status (one life form x product type)
summary.product_x_life.form03 <- summary.product_x_life.form02 %>% mutate(proportion = n.species / total.species)

summary.product_x_life.form03

# encoded as a factor variable
summary.product_x_life.form03 <- summary.product_x_life.form03 %>%
      mutate_at(
            vars(products),
            ~ factor(., levels = c("seeds", "adult", "vegetative"))
      ) %>%
      mutate_at(
            vars(Life_form_new),
            ~ factor(., levels = c("Short-lived", "Long-lived", "Woody"))
      )

View(summary.product_x_life.form03)

levels(summary.product_x_life.form03$products)

# ----------------------------------- plot ----------------------------------- #
#  label names of lf variable
lf.labs <- c("Short-lived\\nherbs", "Long-lived\\nherbs", "Woody")
names(lf.labs) <- c("Short-lived", "Long-lived", "Woody")

# label names of status variable
status.labs <- c("Non-naturalized", "Naturalized\\nnon-invasive", "Invasive")
names(status.labs) <- c("non_nat", "non_invas", "invas")

# user-defined colour
map_color03 <- c("#7bc0f7", "#ffdb69", "#f18226")

# ggplot syntax
plot.product.type_x_lf <- ggplot(
      summary.product_x_life.form03,
      aes(
            y = proportion,
            x = products,
            fill = Life_form_new
      )
) +
      geom_bar(stat = "identity") +
      theme_classic() +
      theme(
            panel.background = element_rect(fill = NA),
            legend.position = "right",
            legend.justification = "left",
            legend.key.height = unit(1, "cm"),
            legend.key.width = unit(1, "cm"),
            strip.text.x = element_text(family = "serif", size = 14),
            strip.text.y = element_text(family = "serif", size = 14),
            strip.background = element_rect(size = 1, linetype = "solid"),
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
            plot.title = element_text(
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
      scale_fill_manual(
            name = "Life form",
            labels = c(
                  "Short-lived herbs",
                  "Long-lived herbs",
                  "Woody species"
            ),
            values = map_color03
      ) +
      scale_x_discrete(name = "", labels = c(
            "Seeds\\n(347)",
            "Live plants\\n(720)",
            "Veg. organs\\n(67)"
      )) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      labs(y = "Proportion of taxa")


plot.product.type_x_lf

# ---------------------------------- output ---------------------------------- #
figure_s5 <- plot.invasion_x_lf +
      plot.product.type_x_lf +
      plot_layout(ncol = 2) +
      plot_annotation(tag_levels = "A")

figure_s5

ggexport(figure_s5,
      filename = "./Output/figure_s5.png",
      width = 4000,
      height = 2000,
      pointsize = 10,
      res = 300
)

# ---------------------------------------------------------------------------- #
#                                   figure s6                                  #
# ---------------------------------------------------------------------------- #
# -------------------------------- figure s6a -------------------------------- #
# the subset of taxa sold only as one type of products are considered here.
# merge species sold only as seed, as live plants, or as vegetative organs
invasion.list.only <- rbind(
      invasion.list %>%
            filter(seed_max == 1 & adult_max == 0 & vegetative_max == 0),
      invasion.list %>%
            filter(seed_max == 0 & adult_max == 1 & vegetative_max == 0),
      invasion.list %>%
            filter(seed_max == 0 & adult_max == 0 & vegetative_max == 1)
)

# data grouped by life form and invasions status
summary.invasion_x_lf.1term <- invasion.list.only %>%
      filter(seed_max == 1 | adult_max == 1 | vegetative_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new, invasion_status) %>%
      dplyr::summarize(n.species = length(unique(TPL_names)))


# renames levels of life form and invasion status
summary.invasion_x_lf.1term01 <- summary.invasion_x_lf.1term %>%
      mutate_at(vars(invasion_status), ~ case_when(
            . == 1 ~ "non_nat",
            . == 2 ~ "non_invas",
            . == 3 ~ "invas"
      )) %>%
      mutate_at(vars(Life_form_new), ~ case_when(
            . == 1 ~ "Short-lived",
            . == 2 ~ "Long-lived",
            . == 3 ~ "Woody"
      ))

View(summary.invasion_x_lf.1term01)

# data grouped by invasion status
summary.invasion_x_lf.1term_s1 <- invasion.list.only %>%
      filter(seed_max == 1 | adult_max == 1 | vegetative_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(invasion_status) %>%
      dplyr::summarize(total.species = length(unique(TPL_names)))

# renames levels of invasion status
summary.invasion_x_lf.1term_s2 <- summary.invasion_x_lf.1term_s1 %>%
      mutate_at(vars(invasion_status), ~ case_when(
            . == 1 ~ "non_nat",
            . == 2 ~ "non_invas",
            . == 3 ~ "invas"
      ))

View(summary.invasion_x_lf.1term_s2)

# summary.invasion_x_lf.1term01 with summary.invasion_x_lf.1term_s2
summary.invasion_x_lf.1term02 <- summary.invasion_x_lf.1term01 %>%
      left_join(summary.invasion_x_lf.1term_s2, by = "invasion_status")
summary.invasion_x_lf.1term02

# calculating proportion of no. of species at one particular status (one life form x invasion status) within one type of product to no. of species at one particular status (one life form x invasion status)
summary.invasion_x_lf.1term03 <- summary.invasion_x_lf.1term02 %>%
      mutate(proportion = n.species / total.species)
summary.invasion_x_lf.1term03

# encoded as a factor variable
summary.invasion_x_lf.1term03 <- summary.invasion_x_lf.1term03 %>%
      mutate_at(
            vars(Life_form_new),
            ~ factor(., levels = c("Short-lived", "Long-lived", "Woody"))
      ) %>%
      mutate_at(
            vars(invasion_status),
            ~ factor(., levels = c("non_nat", "non_invas", "invas"))
      )

View(summary.invasion_x_lf.1term03)

# -------------------------------- figure s6a -------------------------------- #
# label names of life form variable
lf.labs <- c("Short-lived\\nherbs", "Long-lived\\nherbs", "Woody")
names(lf.labs) <- c("Short-lived", "Long-lived", "Woody")

# label names of invasion status
status.labs <- c("Non-naturalized", "Naturalized\\nnon-invasive", "Invasive")
names(status.labs) <- c("non_nat", "non_invas", "invas")

# user-defined colour
map_color03 <- c("#7bc0f7", "#ffdb69", "#f18226")

# ggplot syntax
plot.invasion_x_lf.1term <- ggplot(
      summary.invasion_x_lf.1term03,
      aes(
            y = proportion,
            x = invasion_status,
            fill = Life_form_new
      )
) +
      geom_bar(stat = "identity") +
      theme_classic() +
      theme(
            panel.background = element_rect(fill = NA),
            legend.position = "none",
            strip.text.x = element_text(family = "serif", size = 14),
            strip.text.y = element_text(family = "serif", size = 14),
            strip.background = element_rect(size = 1, linetype = "solid"),
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
            plot.title = element_text(
                  family = "serif",
                  colour = "black",
                  size = 14
            )
      ) +
      scale_fill_manual(
            name = "",
            labels = c(
                  "Short-lived herbs",
                  "Long-lived herbs",
                  "Woody species"
            ),
            values = map_color03
      ) +
      scale_x_discrete(
            name = "",
            labels = c(
                  "Non-naturalized\\n(530)",
                  "Naturalized\\nnon-invasive (81)",
                  "Invasive\\n(28)"
            )
      ) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      labs(y = "Proportion of taxa")


(plot.invasion_x_lf.1term)

# -------------------------------- figure s6b -------------------------------- #
# the subset of taxa sold only as one type of products are considered here.
# seeds
summary.seed_x_lf.1term <- invasion.list.only %>%
      filter(seed_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new) %>%
      dplyr::summarize(n.species = length(unique(TPL_names))) %>%
      mutate(products = "seeds")

# live plants
summary.adult_x_lf.1term <- invasion.list.only %>%
      filter(adult_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new) %>%
      dplyr::summarize(n.species = length(unique(TPL_names))) %>%
      mutate(products = "adult")

# veg. organs
summary.vegetative_x_lf.1term <- invasion.list.only %>%
      filter(vegetative_max == 1) %>%
      filter(!is.na(Life_form_new)) %>%
      group_by(Life_form_new) %>%
      dplyr::summarize(n.species = length(unique(TPL_names))) %>%
      mutate(products = "vegetative")


# merge species sold only as seed, as live plants, or as vegetative organs
summary.product_x_lf.1term <- rbind(
      summary.seed_x_lf.1term,
      summary.adult_x_lf.1term,
      summary.vegetative_x_lf.1term
)

summary.product_x_lf.1term

# renames levels of invasion status
summary.product_x_lf.1term01 <- summary.product_x_lf.1term %>%
      mutate_at(vars(Life_form_new), ~ case_when(
            . == 1 ~ "Short-lived",
            . == 2 ~ "Long-lived",
            . == 3 ~ "Woody"
      ))

# supplement zero values of short-lived and long-lived species
new.row01.1term <- data.frame(
      Life_form_new = c("Short-lived", "Long-lived"),
      n.species = rep(0, 2),
      products = rep("vegetative", 2)
)

new.row01.1term

summary.product_x_lf.1term02 <- rbind(summary.product_x_lf.1term01, new.row01.1term)
View(summary.product_x_lf.1term02)

# data grouped by product type
summary.product_x_lf.1term_s1 <- summary.product_x_lf.1term02 %>%
      group_by(products) %>%
      dplyr::summarize(total.species = sum(n.species))


# summary.product_x_lf.1term02 with summary.product_x_lf.1term_s1
summary.product_x_lf.1term03 <- summary.product_x_lf.1term02 %>%
      left_join(summary.product_x_lf.1term_s1, by = "products")
View(summary.product_x_lf.1term03)

# calculating proportion of no. of species at one particular status (one life form x product type) within one type of product to no. of species at one particular status (one life form x product type)
summary.product_x_lf.1term04 <- summary.product_x_lf.1term03 %>%
      mutate(proportion = n.species / total.species)
View(summary.product_x_lf.1term04)


# encoded as a factor variable
summary.product_x_lf.1term04 <- summary.product_x_lf.1term04 %>%
      mutate_at(
            vars(products),
            ~ factor(., levels = c("seeds", "adult", "vegetative"))
      ) %>%
      mutate_at(
            vars(Life_form_new),
            ~ factor(., levels = c("Short-lived", "Long-lived", "Woody"))
      )

levels(summary.product_x_lf.1term04$products)

levels(summary.product_x_lf.1term04$Life_form_new)

# ----------------------------------- plot ---------------------------------plot
# label names of lf variable
lf.labs <- c("Short-lived\\nherbs", "Long-lived\\nherbs", "Woody")
names(lf.labs) <- c("Short-lived", "Long-lived", "Woody")

# label names of status variable
status.labs <- c("Non-naturalized", "Naturalized\\nnon-invasive", "Invasive")
names(status.labs) <- c("non_nat", "non_invas", "invas")

# user-defined colour
map_color03 <- c("#7bc0f7", "#ffdb69", "#f18226")

# ggplot synstax
plot.product_x_lf.1term04 <- ggplot(
      summary.product_x_lf.1term04,
      aes(
            y = proportion,
            x = products,
            fill = Life_form_new
      )
) +
      geom_bar(stat = "identity") +
      theme_classic() +
      theme(
            panel.background = element_rect(fill = NA),
            legend.position = "right",
            legend.justification = "left",
            legend.key.height = unit(1, "cm"),
            legend.key.width = unit(1, "cm"),
            strip.text.x = element_text(family = "serif", size = 14),
            strip.text.y = element_text(family = "serif", size = 14),
            strip.background = element_rect(size = 1, linetype = "solid"),
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
            plot.title = element_text(
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
      scale_fill_manual(
            name = "Life form",
            labels = c(
                  "Short-lived herbs",
                  "Long-lived herbs",
                  "Woody species"
            ),
            values = map_color03
      ) +
      scale_x_discrete(
      name = "",
      labels = c(
            "Seeds\\n(146)",
            "Live plants\\n(490)",
            "Veg. organs\\n(3)"
      )) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      labs(y = "Proportion of taxa")


(plot.product_x_lf.1term04)

# ---------------------------------- output ---------------------------------- #
figure_s6 <- plot.invasion_x_lf.1term +
      plot.product_x_lf.1term04 +
      plot_layout(ncol = 2) +
      plot_annotation(tag_levels = "A")

figure_s6

ggexport(figure_s6,
      filename = "./Output/figure_s6.png",
      width = 4000,
      height = 2000,
      pointsize = 10,
      res = 300
)
