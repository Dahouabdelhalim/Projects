
# This script replicates Tables S-12 and S-13 in the supplementary materials, which 
# analyze match-level data on cards and referee identity 

# Change the working directory as necessary
setwd("~/Dropbox/iraq-experiment/science/replication/data")

library(tidyverse)
library(lfe)
library(broom)
library(stargazer)
library(glue)
library(ggpubr)
library(viridis)
library(pbapply)

# Read in data ------------------------------------------------------------

df <- read.csv("match-data.csv")


# Table S-12 ------------------------------------------------------------

t.test(df$total_cards[df$type == "mixed"],df$total_cards[df$type != "mixed"] )

t.test(df$yellow_cards[df$type == "mixed"],df$yellow_cards[df$type != "mixed"] )

t.test(df$red_cards[df$type == "mixed"],df$red_cards[df$type != "mixed"] )


# Table S-13 ------------------------------------------------------------

# Regress total, yellow, and red cards individually on match type
# with control for referee identity 


total <- lm(total_cards ~ type + total_goals + mus_ref, data = df)

yellow <- lm(yellow_cards ~ type + total_goals + mus_ref, data = df)

red <- lm(red_cards ~ type + total_goals + mus_ref, data = df)

covlabs <- c("Mixed Status", "Both Treated", "Intercept")

stargazer(total,yellow, red,
          font.size = "small",
          add.lines = list(c("Team Clustered S.E.", rep("\\\\checkmark", 4,
                                                        "Individual Controls", rep("\\\\checkmark", 4)))),
          no.space = F , 
          dep.var.labels = c("Total Cards", "Yellow", "Red"), 
          #omit = c("factor", "goals"),
          omit.stat = c("f", "ser"),
          digits = 2)


