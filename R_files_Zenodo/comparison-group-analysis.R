# This script replicates the  analyses involving the comparison group, producing 
# Figures S-7, S-8 and S-9.
# The script makes use of the datasets waves-1-and-2-comparison.csv, wave-2-comparison.csv,
# and mosul-outcome-wave-2-comparison.csv

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

# Data from waves one and two (training and event attendance outcomes)

merged_pilot <- read.csv("waves-1-and-2-comparison.csv",
                         header=T, stringsAsFactors=FALSE, fileEncoding="latin1")


# Data from wave two (survey-based outcomes)

df <- read.csv("wave-2-comparison.csv",
                         header=T, stringsAsFactors=FALSE, fileEncoding="latin1")

# Data from wave two (Mosul coupon outcome)

mosul_df <- read.csv("mosul-outcome-wave-2-comparison.csv",
                     header=T, stringsAsFactors=FALSE, fileEncoding="latin1")
  

# Creating new "treatment" variables ------------------------------------------------------------
# binary indicators for 1) comparison group, and 2) control groups (comparison  + control) 

merged_pilot <- merged_pilot %>% mutate(comparison = case_when(type == "comparison" ~ "1",
                                                                TRUE ~ "0"),
                                        control = case_when(type == "control" ~ "1",
                                                                TRUE ~ "0") )


df <- df %>% mutate(comparison = case_when(type == "comparison" ~ "1",
                                                                TRUE ~ "0"),
                    control = case_when(type == "control" ~ "1",
                                                                 TRUE ~ "0") )


mosul_df <- mosul_df %>% mutate(comparison = case_when(type == "comparison" ~ "1",
                                            TRUE ~ "0"),
                                control = case_when(type == "control" ~ "1",
                                             TRUE ~ "0"))

mosul_df$block <- ifelse(mosul_df$comparison == 1, NA, mosul_df$block)


# Function for Estimating ATE------------------------------------------------------------


ate <- function(dta, dv, treat, covars, cluster = "team", extravars_extract = NULL) {
  
  # Create OLS formula 
  rhs <- paste(c(treat, covars), collapse = " + ")
  f <- paste(dv, "~", rhs, "| 0 | 0 |", cluster)
  cat("Formula used:", f, "\\n")
  
  # Estimate regression
  m <- felm(formula(f), data = dta)
  
  # Extract results of interest 
  tidy(m, conf.int = T) %>% 
    filter(term %in% c("(Intercept)", treat, extravars_extract))
  
}


# Covariates (note: block exlcuded as there was no 
# randomization block in the comparison leaague (similar to a pure control group)

covariates <- c( "edu" , "church" , "income" , "status1" ,
                 "marital" , "isis_abuse", "player_type",
                 "birth.year")

# Analysis #1: Comparing comparison group to treated group----------------------


# Attend mixed social event

m1 <- ate(merged_pilot[merged_pilot$type != "control",], dv = "attend", treat = "treated", 
          covars = c(covariates),
          extravars_extract = "treated1")


m1

# Training with Muslims at the 6 month mark

m2 <- ate(merged_pilot[merged_pilot$type != "control",], dv = "train_t1", treat = "treated",
          covars = c(covariates, "train_t0"),
          extravars_extract = "treated1")


m2


# Visiting Mosul

m3 <- summary(felm(went_mosul ~ treated  | 0 | 0 | team,
                   data = mosul_df[mosul_df$type != "control",]),
              extravars_extract = "treated1")


m3

# Donating to mixed NGO

m4 <- ate(df[df$type != "control",], dv = "donate_t1_new1", treat = "treated", 
          covars = c(covariates, "donate_t0_new1"),
          extravars_extract = "treated1")

m4



## Registering to a mixed team next season 

m5 <- ate(df[df$type != "control",], dv = "own_group_preference", treat = "treated",
          covars = c(covariates),
          extravars_extract = "treated1")

m5



## Attitudinal indices

# Coexistence

m6 <- ate(df[df$type != "control",], dv = "dv_secular_t1", treat = "treated", 
          covars = c(covariates, "dv_secular_t0"),
          extravars_extract = "treated1")


m6


# Blaming Muslims 

m7 <- ate(df[df$type != "control",], dv = "dv_blame_t1", treat = "treated", 
          covars = c(covariates, "dv_blame_t0"),
          extravars_extract = "treated1")


m7

# Acceptance of Muslims as neighbors

m8 <- ate(df[df$type != "control",], dv = "dv_neigh_t1", treat = "treated", 
          covars = c(covariates, "dv_neigh_t0"),
          extravars_extract = "treated1")


m8

## Visualizing treatment effects

means <- c((m1)$estimate[2],
           (m2)$estimate[2],
           (m3)$coefficients[2],
           (m4)$estimate[2],
           (m5)$estimate[2],
           (m6)$estimate[2],
           (m7)$estimate[2],
           (m8)$estimate[2])

ses <- c((m1)$std.error[2],
         (m2)$std.error[2],
         0.07007981204625308991 ,
         (m4)$std.error[2],
         (m5)$std.error[2],
         (m6)$std.error[2],
         (m7)$std.error[2],
         (m8)$std.error[2])



names <- c("attend event", "train w/ muslims", "visit mosul", "donate mixed NGO", 
           "register mixed team", "natl_unity_index", "blame_index", "neighbor_index")

merged <- data.frame(names, means, ses)

merged <- merged %>%  mutate(lb = means - ses * 1.96, 
                             ub = means + ses * 1.96) 


# Plotting ATE for Treated vs. Comparison (Fig. S-8) 

ggplot(merged, aes(x = reorder(names, - means), y = means)) + 
  geom_line(aes(group = names), color = "dark grey" ) + 
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0, color = "grey") + 
  # facet_wrap(~key) +
  theme_minimal() +
  geom_point(color = "black") + 
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(legend.position = "none", text = element_text(size=15)) +
  xlab("outcome") +
  ylab("treatment effect")




# Analysis #2: Comparing comparison group to control group----------------------

# Attend mixed social event

m1 <- ate(merged_pilot[merged_pilot$type != "treated",], dv = "attend", treat = "control", 
          covars = c(covariates),
          extravars_extract = "control1")


m1

# Training with Muslims at the 6 month mark

m2 <- ate(merged_pilot[merged_pilot$type != "treated",], dv = "train_t1", treat = "control",
          covars = c(covariates, "train_t0"),
          extravars_extract = "control1")


m2

# Visiting Mosul

m3 <- summary(felm(went_mosul ~ control  | 0 | 0 | team,, 
                   data = mosul_df[mosul_df$type != "treated",]),
              extravars_extract = "control1")


m3



# Donating to mixed NGO

m4 <- ate(df[df$type != "treated",], dv = "donate_t1_new1", treat = "control", 
          covars = c(covariates, "donate_t0_new1"),
          extravars_extract = "control1")

m4



## Registering to a mixed team next season 

m5 <- ate(df[df$type != "treated",], dv = "own_group_preference", treat = "control",
          covars = c(covariates),
          extravars_extract = "control1")

m5



## Attitudinal indices

# National Unity

m6 <- ate(df[df$type != "treated",], dv = "dv_secular_t1", treat = "control", 
          covars = c(covariates, "dv_secular_t0"),
          extravars_extract = "control1")


m6


# Blaming Muslims 

m7 <- ate(df[df$type != "treated",], dv = "dv_blame_t1", treat = "control", 
          covars = c(covariates, "dv_blame_t0"),
          extravars_extract = "control1")


m7

# Acceptance of Muslims as neighbors

m8 <- ate(df[df$type != "treated",], dv = "dv_neigh_t1", treat = "control", 
          covars = c(covariates, "dv_neigh_t0"),
          extravars_extract = "control1")


m8


## Visualizing treatment effects

means <- c((m1)$estimate[2],
           (m2)$estimate[2],
           (m3)$coefficients[2],
           (m4)$estimate[2],
           (m5)$estimate[2],
           (m6)$estimate[2],
           (m7)$estimate[2],
           (m8)$estimate[2])

ses <- c((m1)$std.error[2],
         (m2)$std.error[2],
         0.05557 ,
         (m4)$std.error[2],
         (m5)$std.error[2],
         (m6)$std.error[2],
         (m7)$std.error[2],
         (m8)$std.error[2])



names <- c("attend event", "train w/ muslims", "visit mosul", "donate mixed NGO", 
           "register mixed team", "natl_unity_index", "blame_index", "neighbor_index")

merged <- data.frame(names, means, ses)

merged <- merged %>%  mutate(lb = means - ses * 1.96, 
                             ub = means + ses * 1.96) 


# Plotting ATE for Control vs. Comparison (Fig. S-7) 
ggplot(merged, aes(x = reorder(names, - means), y = means)) + 
  geom_line(aes(group = names), color = "dark grey" ) + 
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0, color = "grey") + 
  # facet_wrap(~key) +
  theme_minimal() +
  geom_point(color = "black") + 
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(legend.position = "none", text = element_text(size=15)) +
  xlab("outcome") +
  ylab("treatment effect")


# Directly testing ATEs from Analysis 1 and 2 ------------------------------------------------------------

# Ccontrol vs. comparison, and treated vs. comparison. 

###### Attending mixed event

m1 <- felm(attend ~ treated + comparison + edu + church + income + status1 + marital + isis_abuse
           + player_type + birth.year | 0 | 0 | team, data = merged_pilot)


diff1 <- coef(m1)[2] - coef(m1)[3]
se1 <- sqrt( vcov(m1)[2, 2] + vcov(m1)[3, 3] - 2 * vcov(m1)[2, 3])
t.stat <- abs(diff1 / se1)

# pvalue
p1 <- 2 * (1 - pnorm(t.stat))

###### Training with Muslims
m2 <- felm(train_t1 ~ treated + comparison + edu + church + income + status1 + marital + isis_abuse
           + player_type + birth.year + train_t0  | 0 | 0 | team, data = merged_pilot )


diff2 <- coef(m2)[2] - coef(m2)[3]
se2 <- sqrt( vcov(m2)[2, 2] + vcov(m2)[3, 3] - 2 * vcov(m2)[2, 3])
t.stat <- abs(diff2 / se2)

# pvalue
p2 <- 2 * (1 - pnorm(t.stat))

###### Visiting Mosul 

m3 <- felm(went_mosul ~ treated + comparison, data = mosul_df )

summary(m3)

diff3 <- coef(m3)[2] - coef(m3)[3]
se3 <- sqrt( vcov(m3)[2, 2] + vcov(m3)[3, 3] - 2 * vcov(m3)[2, 3])
t.stat <- abs(diff3 / se3)

# pvalue
p3 <- 2 * (1 - pnorm(t.stat))

###### Donating to mixed NGO
m4 <- felm(donate_t1_new1 ~ treated + comparison + edu + church + income + status1 + marital + isis_abuse
           + player_type + birth.year + donate_t0_new1  | 0 | 0 | team, data = df )


diff4 <- coef(m4)[2] - coef(m4)[3]
se4 <- sqrt( vcov(m4)[2, 2] + vcov(m4)[3, 3] - 2 * vcov(m4)[2, 3])
t.stat <- abs(diff4 / se4)

# pvalue
p4 <- 2 * (1 - pnorm(t.stat))


###### Register mixed team
m5 <- felm(own_group_preference ~ treated + comparison + edu + church + income + status1 + marital + isis_abuse
           + player_type + birth.year  | 0 | 0 | team, data = df )


diff5 <- coef(m5)[2] - coef(m5)[3]
se5 <- sqrt( vcov(m5)[2, 2] + vcov(m5)[3, 3] - 2 * vcov(m5)[2, 3])
t.stat <- abs(diff5 / se5)

# pvalue
p5 <- 2 * (1 - pnorm(t.stat))

###### National unity index
m6 <- felm(dv_secular_t1 ~ treated + comparison + edu + church + income + status1 + marital + isis_abuse
           + player_type + birth.year + dv_secular_t0 | 0 | 0 | team, data = df )


diff6 <- coef(m6)[2] - coef(m6)[3]
se6 <- sqrt( vcov(m6)[2, 2] + vcov(m6)[3, 3] - 2 * vcov(m6)[2, 3])
t.stat <- abs(diff6 / se6)

# pvalue
p6 <- 2 * (1 - pnorm(t.stat))

###### Muslim Blame index
m7 <- felm(dv_blame_t1 ~ treated + comparison + edu + church + income + status1 + marital + isis_abuse
           + player_type + birth.year + dv_blame_t0 | 0 | 0 | team, data = df )


diff7 <- coef(m7)[2] - coef(m7)[3]
se7 <- sqrt( vcov(m7)[2, 2] + vcov(m7)[3, 3] - 2 * vcov(m7)[2, 3])
t.stat <- abs(diff7 / se7)

# pvalue
p7 <- 2 * (1 - pnorm(t.stat))


###### Muslim Neighbor index
m8 <- felm(dv_neigh_t1 ~ + comparison + edu + church + income + status1 + marital + isis_abuse
           + player_type + birth.year + dv_neigh_t0 | 0 | 0 | team, data = df )


diff8 <- coef(m8)[2] - coef(m8)[3]
se8 <- sqrt( vcov(m8)[2, 2] + vcov(m8)[3, 3] - 2 * vcov(m8)[2, 3])
t.stat <- abs(diff8 / se8)

# pvalue
p8 <- 2 * (1 - pnorm(t.stat))

### Visualizing treatment effects


means <- c(diff1, diff2, diff3,  diff4, diff5, diff6, diff7, diff8)

ses <- c(se1, se2, se3, se4, se5, se6, se7, se8)


names <- c("attend event", "train w/ muslims", "visit mosul", "donate mixed NGO", 
           "register mixed team", "natl_unity_index", "blame_index", "neighbor_index")

merged <- data.frame(names, means, ses)

merged <- merged %>%  mutate(lb = means - ses * 1.96, 
                             ub = means + ses * 1.96) 



# Direct comparing ATEs (Fig. S-9) 
ggplot(merged, aes(x = reorder(names, - means), y = means)) + 
  geom_line(aes(group = names), color = "dark grey" ) + 
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0, color = "grey") + 
  theme_minimal() +
  geom_point(color = "black") + 
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(legend.position = "none", text = element_text(size=15)) +
  xlab("outcome") +
  ylab("treatment effect") 



