## Allison Galezo
## Sexual Segregation and Aggregation Statistic for Juveniles

library(dplyr)
library(vegan)
library(ggplot2)

## Set seed for reproducibility.
set.seed(123)

## Load group composition data for surveyed groups with no adults
juvie_surveys <- read.csv("juvenile_surveys_SSAS.csv")

## Create contingency table
ct <- juvie_surveys %>% select(juvenile_females, juvenile_males)

## Function for SSAS (Sexual Segregation and Aggregation Statistic)
## Where ct is a contingency table, 1st column is # females and 2nd column is # males
## See Bonenfant et al. 2007 ("Testing sexual segregation and aggregation: old ways are best") for derivation
SSAS <- function(ct){
  ss.X <- as.numeric(sum(ct[,2]))
  ss.Y <- as.numeric(sum(ct[,1]))
  ss.N <- ss.X + ss.Y
  ss.xi <- ct[,2]
  ss.yi <- ct[,1]
  ss.ni <- rowSums(ct)
  1 - ((ss.N/(ss.X*ss.Y)) * sum((ss.yi*ss.xi)/ss.ni))
}

## Observed SSAS
obs <- SSAS(ct)

## Generate permuted distribution of SSAS values by permuting contingency table
num_sims <- 1000
perms <- permatfull(ct,
                    fixedmar = "both", # preserve both row totals (group size) and column totals (sex ratio)
                    mtype = "count",
                    times = num_sims)
permutedSSAS <- sapply(perms[[3]], SSAS)

## P-value
b <- sum(permutedSSAS > obs) + 1
m <- num_sims + 1
p <- b / m
p

## Sample size (number of surveys)
nrow(juvie_surveys)

## Plot
plotme <- data.frame(value = c(permutedSSAS, obs),
                     identity = c(rep("Expected", num_sims), "Observed"))
ggplot(plotme) +
  geom_density(aes(x = value, fill = identity)) +
  geom_vline(xintercept = obs, color = "black", linetype = "twodash", size = 0.8, show_guide=FALSE) +
  xlab("Sexual Segregation and\\nAggregation Statistic") +
  ylab("Density") +
  xlim(0.7, 0.86) +
  scale_fill_manual(values = c("gray", "black")) +
  guides(fill = guide_legend(title=NULL)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = c(0.18,0.9),
        legend.background = element_rect(fill = "transparent")) +
  scale_y_continuous(expand = c(0,0))
ggsave("JuvieSSAS.tif", plot = last_plot(), device = "tiff", width = 85, height = 85, units = "mm", dpi = 300)


