# Author: Bicheng Dong
# Date: Sun Dec 04 2022
# Last Modified by:
# Last Modified time: Sun Dec 04 2022
# Email: bcdong@bjfu.edu.cn
# Description: R codes for tables 2-3 and figures 2-3
# Results of generalized linear model (GLM) testing effects of invasion status
# (non-naturalized, naturalized non-invasive and invasive) on the
# the number of online nurseries that offered a taxon for sale, and the number
# of provinces and the number of cities where these online nurseries were located,
# and effects on the proportions of online nurseries that offered a taxon for sale
# as seeds, live plants or vegetative organs.

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
library(ggsci)
library(patchwork)
library(plyr)
library(purrr)
library(RColorBrewer)
library(readr)
library(readxl)
library(scales)
library(skimr)
library(tidyr)
library(VGAM)
library(dplyr)

# set work directory as your path
# warning: setup your own path before running codes !!!
file_path <- "your path"

setwd(file_path)
getwd()

# loading data
invasion.list <- fread("./Data/002.alien.plants_online.csv")
str(invasion.list)

# proportion of seed, adult (living plants), and vegetative organs
invasion.list <- invasion.list %>%
    mutate(
        proportion_seed = seed_1 / (seed_1 + seed_0),
        proportion_adult = adult_1 / (adult_1 + adult_0),
        proportion_vegetative = vegetative_1 / (vegetative_1 + vegetative_0)
    )

str(invasion.list)

#
invasion.list$n.company <- as.numeric(invasion.list$n.company)
invasion.list$n.province <- as.numeric(invasion.list$n.province)
invasion.list$n.city <- as.numeric(invasion.list$n.city)

invasion.list$proportion_seed <- as.numeric(invasion.list$proportion_seed)
invasion.list$proportion_adult <- as.numeric(invasion.list$proportion_adult)
invasion.list$proportion_vegetative <- as.numeric(invasion.list$proportion_vegetative)

#
invasion.list$invasion_status <- as.factor(invasion.list$invasion_status)
levels(invasion.list$invasion_status)

# data view
str(invasion.list)
View(invasion.list)

range(invasion.list$n.company)
range(invasion.list$n.province)
range(invasion.list$n.city)
range(invasion.list$proportion_seed)
range(invasion.list$proportion_adult)
range(invasion.list$proportion_vegetative)

#
# refs:
# https://stackoverflow.com/questions/12953045/warning-non-integer-successes-in-a-binomial-glm-survey-packages

# over-dispersion problem
# http://biometry.github.io/APES/LectureNotes/2016-JAGS/Overdispersion/OverdispersionJAGS.html
# Faraway, “Extending the Linear Model with R: Generalized Linear, Mixed Effects and Nonparametric Regression Models.” p88

# https://fukamilab.github.io/BIO202/04-B-binary-data.html
# Overdispersion (variance is larger than mean):
# Needs correction when Phi (= D/(n-P)) > 1.5,
# and quick fix is to use family = quasipoisson,
# but cost is that se of parameters will be multiplied by sqrt(Phi).
# Negative binomial GLM for count data, with overdispersion

fit_glm.n.company <- glm(n.company ~ invasion_status,
    family = poisson(),
    data = invasion.list
)

fit_glm.n.province <- glm(n.province ~ invasion_status,
    family = poisson(),
    data = invasion.list
)

fit_glm.n.city <- glm(n.city ~ invasion_status,
    family = poisson(),
    data = invasion.list
)

dispersiontest(fit_glm.n.company)
dispersiontest(fit_glm.n.province)
dispersiontest(fit_glm.n.city)

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

# -------------------- no. of companies (online nurseries) ------------------- #

# full model
fit_vglm.n.company <- vglm(n.company ~ C1_nonnat_nat + C2_noninv_inv,
    family = pospoisson(), data = invasion.list01
)

# model 1
# drop C2_noninv_inv from full model
fit_vglm.n.company01 <- vglm(n.company ~ C1_nonnat_nat,
    family = pospoisson(), data = invasion.list01
)

# mode 2
# drop C1_nonnat_nat from model 1
fit_vglm.n.company02 <- vglm(n.company ~ 1,
    family = pospoisson(), data = invasion.list01
)

# Error in anova.vglm(fit_vglm.n.province01, fit_vglm.n.province, type = "III",:
# argument 'type' must 'I' or 1 for multiple fits

# Naturalized non-invasive vs. invasive
# anova(model 1 vs. full model)
anova01_vglm.n.company <- anova(fit_vglm.n.company01, fit_vglm.n.company,
    type = "I", test = "LRT"
)

anova01_vglm.n.company

# Non-naturalized vs. naturalized
# anova(model 2 vs. model 1)
anova02_vglm.n.company <- anova(fit_vglm.n.company02, fit_vglm.n.company01,
    type = "I", test = "LRT"
)

anova02_vglm.n.company

# ----------------------------- no. of provinces ----------------------------- #
# full model
fit_vglm.n.province <- vglm(n.province ~ C1_nonnat_nat + C2_noninv_inv,
    family = pospoisson(), data = invasion.list01
)

# model 1
# drop C2_noninv_inv from full model
fit_vglm.n.province01 <- vglm(n.province ~ C1_nonnat_nat,
    family = pospoisson(), data = invasion.list01
)

# mode 2
# drop C1_nonnat_nat from model 1
fit_vglm.n.province02 <- vglm(n.province ~ 1,
    family = pospoisson(), data = invasion.list01
)

# Naturalized non-invasive vs. invasive
# anova(model 1 vs. full model)
anova01_vglm.n.province <- anova(fit_vglm.n.province01, fit_vglm.n.province,
    type = "I", test = "LRT"
)

anova01_vglm.n.province

# Non-naturalized vs. naturalized
# anova(model 2 vs. model 1)
anova02_vglm.n.province <- anova(fit_vglm.n.province02, fit_vglm.n.province01,
    type = "I", test = "LRT"
)

anova02_vglm.n.province

# ------------------------------- no. of cities ------------------------------ #
# full model
fit_vglm.n.city <- vglm(n.city ~ C1_nonnat_nat + C2_noninv_inv,
    family = pospoisson(), data = invasion.list01
)

# model 1
# drop C2_noninv_inv from full model
fit_vglm.n.city01 <- vglm(n.city ~ C1_nonnat_nat,
    family = pospoisson(), data = invasion.list01
)

# model 2
# drop C1_nonnat_nat from model 1
fit_vglm.n.city02 <- vglm(n.city ~ 1,
    family = pospoisson(), data = invasion.list01
)

# Naturalized non-invasive vs. invasive
# anova(model 1 vs. full model)
anova01_vglm.n.city <- anova(fit_vglm.n.city01, fit_vglm.n.city,
    type = "I", test = "LRT"
)

anova01_vglm.n.city

# Non-naturalized vs. naturalized
# anova(model 2 vs. model 1)
anova02_vglm.n.city <- anova(fit_vglm.n.city02, fit_vglm.n.city01,
    type = "I", test = "LRT"
)

anova02_vglm.n.city


# ------------------ proportion of nurseries that sold seeds ----------------- #
# full model
fit_glm.proportion_seed <- glm(proportion_seed ~ C1_nonnat_nat + C2_noninv_inv,
    family = quasibinomial(), data = invasion.list01
)

# model 1
# drop C2_noninv_inv from full model
fit_glm.proportion_seed01 <- glm(proportion_seed ~ C1_nonnat_nat,
    family = quasibinomial(), data = invasion.list01
)

# mode 2
# drop C1_nonnat_nat from model 1
fit_glm.proportion_seed02 <- glm(proportion_seed ~ 1,
    family = quasibinomial(), data = invasion.list01
)

# Naturalized non-invasive vs. invasive
# anova(model 1 vs. full model)
anova01_glm.proportion_seed <- anova(fit_glm.proportion_seed01, fit_glm.proportion_seed, test = "Chisq")

anova01_glm.proportion_seed

# Non-naturalized vs. naturalized
# anova(model 2 vs. model 1)
anova02_glm.proportion_seed <- anova(fit_glm.proportion_seed02, fit_glm.proportion_seed01, test = "Chisq")

anova02_glm.proportion_seed

# -------------- proportion of nurseries that sold living plants ------------- #
# adult = living plants
# full model
fit_glm.proportion_adult <- glm(proportion_adult ~ C1_nonnat_nat + C2_noninv_inv, family = quasibinomial(), data = invasion.list01)

# model 1
# drop C2_noninv_inv from full model
fit_glm.proportion_adult01 <- glm(proportion_adult ~ C1_nonnat_nat,
    family = quasibinomial(), data = invasion.list01
)

# model 2
# drop C1_nonnat_nat from model 1
fit_glm.proportion_adult02 <- glm(proportion_adult ~ 1,
    family = quasibinomial(), data = invasion.list01
)

# Naturalized non-invasive vs. invasive
# anova(model 1 vs. full model)
anova01_glm.proportion_adult <- anova(fit_glm.proportion_adult01, fit_glm.proportion_adult, test = "Chisq")

anova01_glm.proportion_adult

# Non-naturalized vs. naturalized
# anova(model 2 vs. model 1)
anova02_glm.proportion_adult <- anova(fit_glm.proportion_adult02, fit_glm.proportion_adult01, test = "Chisq")

anova02_glm.proportion_adult

# ------------ proportion of nurseries that sold vegetative organs ----------- #
# full model
fit_glm.proportion_vegetative <- glm(proportion_vegetative ~ C1_nonnat_nat + C2_noninv_inv, family = quasibinomial(), data = invasion.list01)

# model 1
# drop C2_noninv_inv from full model
fit_glm.proportion_vegetative01 <- glm(proportion_vegetative ~ C1_nonnat_nat, family = quasibinomial(), data = invasion.list01)

# mode 2
# drop C1_nonnat_nat from model 1
fit_glm.proportion_vegetative02 <- glm(proportion_vegetative ~ 1, family = quasibinomial(), data = invasion.list01)

# Naturalized non-invasive vs. invasive
# anova(model 1 vs. full model)
anova01_glm.proportion_vegetative <- anova(fit_glm.proportion_vegetative01, fit_glm.proportion_vegetative, type = "I", test = "Chisq")

anova01_glm.proportion_vegetative

# Non-naturalized vs. naturalized
# anova(model 2 vs. model 1)
anova02_glm.proportion_vegetative <- anova(fit_glm.proportion_vegetative02, fit_glm.proportion_vegetative01, type = "I", test = "Chisq")

anova02_glm.proportion_vegetative

# ------------------------------- anova summary ------------------------------ #
#
glm_models <- list(
    anova02_vglm.n.company = anova02_vglm.n.company,
    anova01_vglm.n.company = anova01_vglm.n.company,
    anova02_vglm.n.province = anova02_vglm.n.province,
    anova01_vglm.n.province = anova01_vglm.n.province,
    anova02_vglm.n.city = anova02_vglm.n.city,
    anova01_vglm.n.city = anova01_vglm.n.city,
    anova02_glm.proportion_seed = anova02_glm.proportion_seed,
    anova01_glm.proportion_seed = anova01_glm.proportion_seed,
    anova02_glm.proportion_adult = anova02_glm.proportion_adult,
    anova01_glm.proportion_adult = anova01_glm.proportion_adult,
    anova02_glm.proportion_vegetative = anova02_glm.proportion_vegetative,
    anova01_glm.proportion_vegetative = anova01_glm.proportion_vegetative
)

#
res_glm_tidy <- map_dfr(glm_models, tidy, .id = "variable")

#
res_glm_tidy <- res_glm_tidy %>%
    mutate(contrast.name = case_when(
        grepl("anova01", variable) ~ "C1_noninv_inv",
        grepl("anova02", variable) ~ "C1_nonnat_nat"
    ))

# data view
View(res_glm_tidy)

# save data
write.csv(res_glm_tidy, "./Output/anova.summary_table2_3.csv", row.names = FALSE)

# ---------------------------------------------------------------------------- #
#                                     plot                                     #
# ---------------------------------------------------------------------------- #
# display colors
map_color <- c("#a9b8c6", "#6f80be", "#ff8884")
show_col(map_color)

# -------------------- user-specific ggplot batch function ------------------- #
# refs
# https://blog.csdn.net/weixin_43151909/article/details/107185648

ggplot.batch <- function(variable, status, map_color, data_source) {

    # subset columns
    data.used <- data_source %>% select(variable, status)
    read_name <- names(data.used)

    # rename columns
    colnames(data.used) <- c("variable", "status")

    # remove na data
    data.used <- data.used %>% filter(!is.na(variable))

    # main code for numeric variable
    if (is.numeric(data.used$variable)) {

        # data summary
        data.used_01 <- plyr::ddply(data.used, plyr::.(status), summarise,
            n = length(variable),
            mean = mean(variable),
            sd = sd(variable),
            se = sd / sqrt(n)
        )

        # encoded as a factor
        data.used_01$status <- as.factor(data.used_01$status)

        # ggplot syntax
        isNum_plot.01 <- ggplot(
            data = data.used_01, aes(x = status, y = mean, fill = status)
        ) +
            geom_bar(colour = "black", stat = "identity") +
            geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2
            ) +
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
                    linetype = "solid", size = 1
                ),
                axis.text.x = element_text(
                    family = "serif",
                    colour = "black", size = 10
                ),
                axis.text.y = element_text(
                    family = "serif",
                    colour = "black", size = 10
                ),
                axis.title = element_text(
                    family = "serif",
                    colour = "black", size = 10
                )
            ) +
            labs(x = read_name[2], y = read_name[1])

        return(isNum_plot.01)
    }

    # main code for factor variable
    if (is.factor(data.used$variable)) {

        # encoded as a factor
        data.used$variable <- as.factor(data.used$variable)

        # data summary
        data.used_02 <- plyr::ddply(data.used, plyr::.(status), summarise,
            n_not_sold = length(variable[variable == "0"]),
            n_sold = length(variable[variable == "1"]),
            n_all_species = n_not_sold + n_sold,
            proportion = n_sold / n_all_species
        )
        (data.used_02)

        # ggplot syntax
        isFact_plot.01 <- ggplot(
            data.used_02,
            aes(x = status, y = proportion, fill = as.factor(status))
        ) +
            geom_bar(colour = "black", stat = "identity") +
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
                axis.text.x = element_text(
                    family = "serif",
                    colour = "black",
                    size = 10
                ),
                axis.text.y = element_text(
                    family = "serif",
                    colour = "black",
                    size = 10
                ),
                axis.title = element_text(
                    family = "serif",
                    colour = "black",
                    size = 10
                )
            ) +
            labs(x = read_name[2], y = "proportion")

        return(isFact_plot.01)
    }
}

# ------------------------------ batching plots ------------------------------ #
# variable list
glm.vars <- c(
    "n.company",
    "n.province",
    "n.city",
    "proportion_seed",
    "proportion_adult",
    "proportion_vegetative"
)

# variable names
glm.vars <- set_names(glm.vars)
(glm.vars)

# using "map" to transform all variables to ggplot-batch function
status <- "invasion_status"
ggplot_groups <- glm.vars %>% map(~ ggplot.batch(
    variable = .x,
    status = status,
    map_color = map_color,
    data_source = invasion.list01
))

(ggplot_groups)

# ggplot list
names(ggplot_groups)

# ---------------------- adding plots with more elements --------------------- #
# ranges of variables
range(invasion.list01$n.company, na.rm = TRUE)
range(invasion.list01$n.province, na.rm = TRUE)
range(invasion.list01$n.city, na.rm = TRUE)
range(invasion.list01$proportion_seed, na.rm = TRUE)
range(invasion.list01$proportion_adult, na.rm = TRUE)
range(invasion.list01$proportion_vegetative, na.rm = TRUE)

# line segments
# geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df)

# data of line segments
plot.scales <- data.frame(
    plot.name = c(
        "n.company",
        "n.province",
        "n.city",
        "proportion_seed",
        "proportion_adult",
        "proportion_vegetative"
    ),
    ymax = c(30, 6, 10, 1, 1, 0.05)
)

# no. of online nurseries
bar_plot_n.company <- ggplot_groups$n.company +
    labs(
        y = "No. online nurseries",
        x = "",
        size = 10
    ) +
    theme(axis.text.x = element_blank()) +
    scale_y_continuous(
        limits = c(0, plot.scales$ymax[1]),
        breaks = seq(0, plot.scales$ymax[1],
            by = 5
        )
    ) +
    scale_x_discrete(labels = c(
        "Non-naturalized",
        "Naturalized\\nnon-invasive",
        "Invasive"
    )) +
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
        size = 4,
        label = "***"
    ) +
    annotate("text",
        x = 3.45,
        y = plot.scales$ymax[1] * 0.885,
        size = 4,
        label = "***"
    )

# no. of provinces
bar_plot_n.province <- ggplot_groups$n.province +
    labs(
        y = "No. provinces",
        x = "",
        size = 10
    ) +
    theme(axis.text.x = element_blank()) +
    scale_y_continuous(
        limits = c(0, 6),
        breaks = seq(0, 6, by = 1)
    ) +
    scale_x_discrete(labels = c(
        "Non-naturalized",
        "Naturalized\\nnon-invasive",
        "Invasive"
    )) +
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
        size = 4,
        label = "***"
    ) +
    annotate("text",
        x = 3.45,
        y = plot.scales$ymax[2] * 0.885,
        size = 4,
        label = "ns"
    )

# no. of cities
bar_plot_n.city <- ggplot_groups$n.city +
    labs(
        y = "No. cities",
        x = "",
        size = 10
    ) +
    scale_y_continuous(
        limits = c(0, 10),
        breaks = seq(0, 10, by = 2)
    ) +
    scale_x_discrete(labels = c(
        "Non-naturalized",
        "Naturalized\\nnon-invasive",
        "Invasive"
    )) +
    geom_segment(aes(
        x = 0.7,
        y = plot.scales$ymax[3],
        xend = 0.7 + 0.6,
        yend = plot.scales$ymax[3] + 0
    )) +
    geom_segment(aes(
        x = 1.7,
        y = plot.scales$ymax[3],
        xend = 2.7 + 0.6,
        yend = plot.scales$ymax[3] + 0
    )) +
    geom_segment(aes(
        x = 1.7,
        y = plot.scales$ymax[3] * 0.9,
        xend = 1.7 + 0.6,
        yend = plot.scales$ymax[3] * 0.9 + 0
    )) +
    geom_segment(aes(
        x = 2.7,
        y = plot.scales$ymax[3] * 0.9,
        xend = 2.7 + 0.6,
        yend = plot.scales$ymax[3] * 0.9 + 0
    )) +
    annotate("text",
        x = 3.45,
        y = plot.scales$ymax[3] * 0.985,
        size = 4,
        label = "***"
    ) +
    annotate("text",
        x = 3.45,
        y = plot.scales$ymax[3] * 0.885,
        size = 4,
        label = "ns"
    )

# proportion of online nurseries that sold seeds
bar_plot_proportion_seed <- ggplot_groups$proportion_seed +
    labs(
        y = "Proportion of nurseries\\nthat sold seeds",
        x = "",
        size = 10
    ) +
    theme(axis.text.x = element_blank()) +
    scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2)
    ) +
    scale_x_discrete(labels = c(
        "Non-naturalized",
        "Naturalized\\nnon-invasive",
        "Invasive"
    )) +
    geom_segment(aes(
        x = 0.7,
        y = plot.scales$ymax[4],
        xend = 0.7 + 0.6,
        yend = plot.scales$ymax[4] + 0
    )) +
    geom_segment(aes(
        x = 1.7,
        y = plot.scales$ymax[4],
        xend = 2.7 + 0.6,
        yend = plot.scales$ymax[4] + 0
    )) +
    geom_segment(aes(
        x = 1.7,
        y = plot.scales$ymax[4] * 0.9,
        xend = 1.7 + 0.6,
        yend = plot.scales$ymax[4] * 0.9 + 0
    )) +
    geom_segment(aes(
        x = 2.7,
        y = plot.scales$ymax[4] * 0.9,
        xend = 2.7 + 0.6,
        yend = plot.scales$ymax[4] * 0.9 + 0
    )) +
    annotate("text",
        x = 3.45,
        y = plot.scales$ymax[4] * 0.985,
        size = 4,
        label = "***"
    ) +
    annotate("text",
        x = 3.45,
        y = plot.scales$ymax[4] * 0.885,
        size = 4,
        label = "*"
    )

# proportion of online nurseries that sold living plants
bar_plot_proportion_adult <- ggplot_groups$proportion_adult +
    labs(
        y = "Proportion of nurseries\\nthat sold live plants",
        x = "",
        size = 10
    ) +
    theme(axis.text.x = element_blank()) +
    scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2)
    ) +
    scale_x_discrete(labels = c(
        "Non-naturalized",
        "Naturalized\\nnon-invasive",
        "Invasive"
    )) +
    geom_segment(aes(
        x = 0.7,
        y = plot.scales$ymax[5],
        xend = 0.7 + 0.6,
        yend = plot.scales$ymax[5] + 0
    )) +
    geom_segment(aes(
        x = 1.7,
        y = plot.scales$ymax[5],
        xend = 2.7 + 0.6,
        yend = plot.scales$ymax[5] + 0
    )) +
    geom_segment(aes(
        x = 1.7,
        y = plot.scales$ymax[5] * 0.9,
        xend = 1.7 + 0.6,
        yend = plot.scales$ymax[5] * 0.9 + 0
    )) +
    geom_segment(aes(
        x = 2.7,
        y = plot.scales$ymax[5] * 0.9,
        xend = 2.7 + 0.6,
        yend = plot.scales$ymax[5] * 0.9 + 0
    )) +
    annotate("text",
        x = 3.45,
        y = plot.scales$ymax[5] * 0.985,
        size = 4,
        label = "***"
    ) +
    annotate("text",
        x = 3.45,
        y = plot.scales$ymax[5] * 0.885,
        size = 4,
        label = "*"
    )

# proportion of online nurseries that sold vegetative organs
bar_plot_proportion_vegetative <- ggplot_groups$proportion_vegetative +
    labs(
        y = "Proportion of nurseries\\nthat sold vegetative organs",
        x = "",
        size = 10
    ) +
    scale_y_continuous(
        limits = c(0, 0.05),
        breaks = seq(0, 0.05, by = 0.01)
    ) +
    scale_x_discrete(labels = c(
        "Non-naturalized",
        "Naturalized\\nnon-invasive",
        "Invasive"
    )) +
    geom_segment(aes(
        x = 0.7,
        y = plot.scales$ymax[6],
        xend = 0.7 + 0.6,
        yend = plot.scales$ymax[6] + 0
    )) +
    geom_segment(aes(
        x = 1.7,
        y = plot.scales$ymax[6],
        xend = 2.7 + 0.6,
        yend = plot.scales$ymax[6] + 0
    )) +
    geom_segment(aes(
        x = 1.7,
        y = plot.scales$ymax[6] * 0.9,
        xend = 1.7 + 0.6,
        yend = plot.scales$ymax[6] * 0.9 + 0
    )) +
    geom_segment(aes(
        x = 2.7,
        y = plot.scales$ymax[6] * 0.9,
        xend = 2.7 + 0.6,
        yend = plot.scales$ymax[6] * 0.9 + 0
    )) +
    annotate("text",
        x = 3.45,
        y = plot.scales$ymax[6] * 0.985,
        size = 4,
        label = "**"
    ) +
    annotate("text",
        x = 3.45,
        y = plot.scales$ymax[6] * 0.885,
        size = 4,
        label = "ns"
    )

# plot lists
bar_plot_n.company

bar_plot_n.province

bar_plot_n.city

bar_plot_proportion_seed

bar_plot_proportion_adult

bar_plot_proportion_vegetative

# ---------------------------------- output ---------------------------------- #
# figure 02
figure02 <- bar_plot_n.company +
    bar_plot_n.province +
    bar_plot_n.city +
    plot_layout(ncol = 1) +
    plot_annotation(tag_levels = "A")

(figure02)

# figure 03
figure03 <- bar_plot_proportion_seed +
    bar_plot_proportion_adult +
    bar_plot_proportion_vegetative +
    plot_layout(ncol = 1) +
    plot_annotation(tag_levels = "A")

(figure03)

# export figure02
ggexport(figure02,
    filename = "./Output/figure02.png",
    width = 2400,
    height = 4800,
    pointsize = 10,
    res = 600
)

# export figure03
ggexport(figure03,
    filename = "./Output/figure03.png",
    width = 2400,
    height = 4800,
    pointsize = 10,
    res = 600
)
