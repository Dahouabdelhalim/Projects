# Author: Bicheng Dong
# Date: Sat Dec 03 2022
# Last Modified by:
# Last Modified time: Sat Dec 03 2022
# Email: bcdong@bjfu.edu.cn
# Description: R codes for table 1 and figure 1
# Results of generalized linear model (GLM) testing effects of invasion status
# (non-naturalized, naturalized non-invasive and invasive) on the
# online-nursery availability of cultivated alien plant taxa.
# warning: setup your own path before running codes !!!

# -------------------------------- preparation ------------------------------- #
# cleaning memory
cat("\\014")
rm(list = ls())
gc()

# loading packages
library(data.table)
library(ggplot2)
library(ggpubr)
library(scales)
library(tidyr)
library(dplyr)

# set work directory as your path
# warning: setup your own path before running codes !!!
file_path <- "your path"

setwd(file_path)
getwd()

# loading data
invasion.list <- fread("./Data/001.alien.plants_in_China.csv")
str(invasion.list)

#
invasion.list <- data.frame(invasion.list)
invasion.list$invasion_status <- as.factor(invasion.list$invasion_status)

# data table of invasion status
# 1: non_naturalized
# 2: naturalized but non_invasive
# 3: invasive
table(invasion.list$invasion_status)
str(invasion.list)

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

contrast_comp <- contrast_comp %>%
      mutate(invasion_status = c("1", "2", "3"))

names(contrast_comp)[1:2] <- c("C1_nonnat_nat", "C2_noninv_inv")
contrast_comp

# merge invasion.list with contrast_comp
invasion.list01 <- invasion.list %>%
      left_join(contrast_comp, by = "invasion_status")

View(invasion.list01)

# ------------------------------ likelihood test ----------------------------- #
# main code of table 1
# full model
fit_glm.ali <- glm(alibaba_status ~ C1_nonnat_nat + C2_noninv_inv,
      family = binomial("cloglog"), data = invasion.list01
)

# model 1
# drop C2_noninv_inv from full model
fit_glm.ali01 <- glm(alibaba_status ~ C1_nonnat_nat,
      family = binomial("cloglog"), data = invasion.list01
)

# model 2
# drop C1_nonnat_nat from model 1
fit_glm.ali02 <- glm(alibaba_status ~ 1,
      family = binomial("cloglog"), data = invasion.list01
)

# Naturalized non-invasive vs. invasive
# anova(model 1 vs. full model)
anova01_glm.ali <- anova(fit_glm.ali01, fit_glm.ali, test = "Chisq")
anova01_glm.ali

# Non-naturalized vs. naturalized
# anova(model 2 vs. model 1)
anova02_glm.ali <- anova(fit_glm.ali02, fit_glm.ali01, test = "Chisq")
anova02_glm.ali

# ------------------------------- data summary ------------------------------- #
#
invasion.list01$invasion_status <- as.factor(invasion.list01$invasion_status)
levels(invasion.list01$invasion_status)

# data summary
# alibaba_status: 1: online availability; 0: not online availability
invasion.list.summary <- plyr::ddply(invasion.list01,
      plyr::.(invasion_status),
      summarise,
      n_not_sold = length(alibaba_status[alibaba_status == "0"]),
      n_sold = length(alibaba_status[alibaba_status == "1"]),
      n_all_species = n_not_sold + n_sold,
      proportion = n_sold / n_all_species
)

# add a new column named pt_proportion, which is proportion of number of alien
# plants within each invasion status to total number of alien plants in China
invasion.list.summary <- invasion.list.summary %>%
      mutate(pt_proportion = n_all_species / sum(n_all_species))

(invasion.list.summary)

# ----------------------------------- plot ----------------------------------- #
# display colors
map_color <- c("#a9b8c6", "#6f80be", "#ff8884")
show_col(map_color)

# main code of figure 1
alibaba.status_plot <- ggplot(
      invasion.list.summary,
      aes(
            x = invasion_status,
            y = proportion,
            width = pt_proportion,
            fill = as.factor(invasion_status)
      )
) +
      geom_bar(stat = "identity") +
            scale_fill_manual(values = map_color) +
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
                  colour = "black", size = 10
            ),
            axis.title = element_text(
                  family = "serif",
                  colour = "black", size = 10
            )
      ) +
      labs(
            y = "Proportion of taxa\\navailable for online nurseries",
            x = "",
            size = 10
      ) +
      scale_y_continuous(
            limits = c(0, 0.4),
            breaks = seq(0, 0.4, by = 0.1)
      ) +
      geom_text(aes(
            x = invasion_status,
            y = proportion + 0.015,
            label = round(proportion, 2),
            size = 8
      )) +
      scale_x_discrete(labels = c(
            "Non-naturalized",
            "Naturalized\\nnon-invasive",
            "Invasive"
      )) +
      geom_segment(aes(
            x = 0.7, y = 0.4,
            xend = 0.7 + 0.6, yend = 0.4 + 0
      )) +
      geom_segment(aes(
            x = 1.7, y = 0.4,
            xend = 2.7 + 0.6, yend = 0.4 + 0
      )) +
      geom_segment(aes(
            x = 1.7, y = 0.4 * 0.9,
            xend = 1.7 + 0.6, yend = 0.4 * 0.9 + 0
      )) +
      geom_segment(aes(
            x = 2.7, y = 0.4 * 0.9,
            xend = 2.7 + 0.6, yend = 0.4 * 0.9 + 0
      )) +
      annotate("text",
            x = 3.45, y = 0.4 * 0.985,
            size = 4, label = "***"
      ) +
      annotate("text",
            x = 3.45, y = 0.4 * 0.885,
            size = 4, label = " * "
      )


(alibaba.status_plot)

# ---------------------------------- output ---------------------------------- #
ggexport(alibaba.status_plot,
      filename = "./Output/figure01.png",
      width = 2200,
      height = 2000,
      pointsize = 10,
      res = 600
)