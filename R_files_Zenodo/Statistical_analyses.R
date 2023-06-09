

# Data and script to accompany the manuscript: 

# "Analysis of within-individual variation in extra-pair paternity in blue tits (Cyanistes caeruleus) 
#  shows low repeatability and little effect of changes in neighborhood "



#  Load required packages
sapply(c('lme4','multcomp','data.table','rptR','reshape2','dplyr','broom','car','ggplot2','patchwork','ggpubr'),
       function(x) suppressMessages(require (x ,character.only =TRUE) ) )

#  Set working directoty
setwd("C:/Users/kbeck/Desktop/data to submit")

#   1. Repeatability analyses
####--------------------------------------------------------------------------------------------------

# Load data
load("data_repeatability.RData")

# Examples are all shown for the data on females
# For the data on males:
# rep_males <- data_repeatability[[2]] 
rep_females <- data_repeatability[[1]]   # includes all females
# If one wants to select only females with completely genotyped nests: 
# rep_females <- rep_females[complete==1,] 


# Extra-pair offspring (EPO)
print(rpt(EPO ~ (1 | mother) , grname = "mother", 
          data = rep_females,  datatype = "Poisson", nboot = 1000, npermut = 0, parallel = TRUE))
# EPP yes/no
print(rpt(epp_bin ~ (1 | mother) , grname = "mother", 
          data = rep_females,  datatype = "Binary", nboot = 1000, npermut = 0, parallel = TRUE))

# Adjusted repeatability for clutch size
# Example for the number of extra-pair offspring
print(rpt(EPO ~ clutch + (1 | mother), grname = "mother", 
                  data = rep_females,  datatype = "Poisson", nboot = 1000, npermut = 0, parallel = TRUE))  

# Including nestbox identity
# Example for the number of extra-pair offspring
print(rpt(EPO ~ (1 | mother) + (1|box) , grname = c("mother","box"), 
          data = rep_females,  datatype = "Poisson", nboot = 1000, npermut = 0, parallel = TRUE))  

# Including pair identity
# Example for the number of extra-pair offspring
print(rpt(EPO ~ (1 | mother) + (1|pair) , grname = c("mother","pair"), 
          data = rep_females,  datatype = "Poisson", nboot = 1000, npermut = 0, parallel = TRUE))

####--------------------------------------------------------------------------------------------------




#   2. Between-year changes in EPP in relation to changes in the breeding environment
####--------------------------------------------------------------------------------------------------

# Load data
load("data_environmental_changes.RData")

# Examples are all shown for the data on females 
# For the data on males:
# data_yearling_adult <- data_environmental_changes[[2]]   # Data including males turning from yearling to adult
# data_only_adult <- data_environmental_changes[[3]]       # Data including only adult males 

data_female <- data_environmental_changes[[1]]

# Change in EPO (for change in number of EP partners, simply change dependent variable in d_EP_males)
delta_f1 <- lmer(d_EPO ~ d_sum_neig + d_tar_soc_male + scale(d_area) +
                   + pair_consistent + d_prop_fam_males + d_prop_fam_females + d_prop_juv +
                   d_mean_neig_tar + fam_soc_male + fam_ep_male + d_prop_juv_fem + d_female_neig_tar +
                   (1|year_) + (1|focal_female), data=data_female) 
K <- diag(length(fixef(delta_f1)))[-1,]; rownames(K) <- names(fixef(delta_f1))[- 1]
glht(delta_f1,  linfct = K) %>% summary %>% tidy %>% data.table
vif(delta_f1)
qqnorm(residuals(delta_f1))
plot(delta_f1)
acf(resid(delta_f1))


# delta EPP yes/no
delta_f2 <- glmer(abs(d_ep_bin) ~ abs(d_sum_neig) + abs(d_tar_soc_male) + abs(scale(d_area)) +
                    + pair_consistent + abs(d_prop_fam_males) + abs(d_prop_fam_females) + abs(d_prop_juv) +
                    abs(d_mean_neig_tar) + fam_soc_male + fam_ep_male + abs(d_prop_juv_fem) + abs(d_female_neig_tar) +
                    (1|year_) + (1|focal_female), data=data_female, family="binomial",glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
K <- diag(length(fixef(delta_f2)))[-1,]; rownames(K) <- names(fixef(delta_f2))[-1]
glht(delta_f2,  linfct = K) %>% summary %>% tidy %>% data.table
vif(delta_f2)

####--------------------------------------------------------------------------------------------------



#  Plot
####--------------------------------------------------------------------------------------------------

# Example for female data

# Change in number of EPO
####-----------------------------------------------------------------------

# Example for "Territory size"
# Plot raw variables
p11 <- ggplot(data_female, aes(d_area, d_EPO)) + 
  geom_jitter() + 
  scale_y_continuous(breaks=seq(-5,5,5)) + theme_classic() +
  theme(axis.text=element_text(size=15, color = "black"), axis.title=element_text(size=15),legend.position="none")
p11 <- p11 + labs(x = "\\nΔ Territory size\\n", y= "\\nΔ EPY")
# Explanatory variable
p111 <- ggplot(data_female, aes(d_area)) + 
  geom_histogram(col="black", fill="grey50") + 
  theme_void() + theme(legend.position="none")   
# Dependent variable
p1 <- ggplot(data_female, aes(d_EPO)) + 
  geom_histogram(col="black", fill="grey50") + 
  theme_void() + theme(legend.position="none") + rotate()
pl1 <- p111 + p11 + plot_layout(ncol=1)


# Example for "Former social partner present"
# Plot raw variables
data_female$fam_soc_male2 <- ifelse(data_female$fam_soc_male==1, "No", "Yes")
p11 <- ggplot(data_female, aes(fam_soc_male2, d_EPO), fill="grey50") + 
  geom_boxplot(fill="grey50") + 
  scale_y_continuous(breaks=seq(-5,5,5)) + theme_classic() +
  theme(axis.text=element_text(size=15, color = "black"), axis.title=element_text(size=15),legend.position="none")
p11 <- p11 + labs(x = "\\n Former social partner present\\n", y= "\\nΔ EPY")
# Explanatory variable
data_female$fam_soc_male <- as.numeric(data_female$fam_soc_male)
p111 <- ggplot(data_female, aes(fam_soc_male)) + 
  geom_histogram(col="black", fill="grey50") + 
  theme_void() + theme(legend.position="none")   # theme_void -> remove all axis
# Dependent variable
p1 <- ggplot(data_female, aes(d_EPO)) + 
  geom_histogram(col="black", fill="grey50") + 
  theme_void() + theme(legend.position="none") + rotate()
pl2 <- p111 + p11 + plot_layout(ncol=1)


####--------------------------------------------------------------------------------------------------




