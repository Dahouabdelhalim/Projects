# How female-female competition affects male-male competition: 
# insights on post-copulatory sexual selection from socially polyandrous species.    

# 11.23.2021

# Sara E. Lipshutz, Samuel J. Torneo, Kimberly A. Rosvall

# Set path
setwd()

# Read in datasheet - you'll have to specify the correct path
sperm <-read.csv("jacana_sperm_measurements.csv") 
sperm.avg <-read.csv("sperm_measurements_averages.csv")

# subset dataset for each species
spinosa.sperm.avg = subset(sperm.avg, Species == "Northern")
jacana.sperm.avg = subset(sperm.avg, Species == "Wattled")

#Examine structure of the data
str(sperm)
str(sperm.avg)

# Export figures to pptx using officer and rvg 
require (officer)
require(rvg)
require(ggpubr)
require(ggplot2)

# Compare species - rough overlook
library(ggplot2)
ggplot(sperm, aes(x=Tail, color = Species)) + geom_histogram()
# spionsa have longer tail
ggplot(sperm, aes(x=Midpiece, color = Species)) + geom_histogram()
# spionsa have longer midpiece
ggplot(sperm, aes(x=Head, color = Species)) + geom_histogram()
# Head length seems like its the same between species


# Check for outliers - decide whether to exclude any samples
ggplot(sperm, aes(x=Individual, y=Tail, fill=Breeding)) + geom_boxplot() # WAM7 looks like an outlier
ggplot(sperm, aes(x=Individual, y=Midpiece, fill=Breeding)) + geom_boxplot()
ggplot(sperm, aes(x=Individual, y=Head, fill=Breeding)) + geom_boxplot()

library(outliers)
grubbs.test(jacana.sperm$Tail, type = 10)
#G = 2.94924, U = 0.87752, p-value = 0.08774
#alternative hypothesis: lowest value 56.328 is an outlier

grubbs.test(spinosa.sperm$Tail, type = 10) 
#G = 2.89249, U = 0.92395, p-value = 0.182
#alternative hypothesis: lowest value 65.37 is an outlier

#### Decided not to exclude WAM7 - outlier test not significant


# Validation of 10 sperm per individual

# Dataset: 50-60 sperm samples from 1 individual of each species, 
# with measurements on total length, tail, midpiece, and head, along with 
# 10 sperm samples from those same individuals

method <- read.csv("sperm_sample_comparison.csv")

# Total length - methods are comparable
ggplot(method, aes(x=Method, y=Total.Length, fill=Species)) + geom_boxplot()
ggboxplot(method, x = "Method", y="Total.Length", color="Species", add = "jitter")
method.total.length <- aov(Total.Length ~ Method * Species, data = method)
anova(method.total.length)

# Assessing technical repeatability of sperm within an individual

# Repeatability using rptr
#install.packages("rptR")
library(rptR)

repeatability <- read.csv("Repeatability samples.csv")

rpt(Tail ~ (1 | Sample), grname = "Sample", data = repeatability, datatype = "Gaussian", nboot = 1000, npermut = 0)
# R  = 0.851, SE = 0.07, CI = [0.669, 0.938], P  = 2.65e-07 [LRT], NA [Permutation]
rpt(Midpiece ~ (1 | Sample), grname = "Sample", data = repeatability, datatype = "Gaussian", nboot = 1000, npermut = 0)
# R  = 0.803, SE = 0.089, CI = [0.574, 0.918]. P  = 7.68e-07 [LRT]. NA [Permutation]
rpt(Head ~ (1 | Sample), grname = "Sample", data = repeatability, datatype = "Gaussian", nboot = 1000, npermut = 0)
# R  = 0.924, SE = 0.04, CI = [0.82, 0.968], P  = 1.12e-11 [LRT]. NA [Permutation]


# MANOVA
# Followed this R tutorial: https://www.datanovia.com/en/lessons/one-way-manova-in-r/
library(tidyverse)
library(ggpubr)
library(rstatix)
library(car)
library(broom)

ggboxplot(sperm.avg, x = "Species", y = c("Tail", "Midpiece","Head"),  merge = TRUE, palette = "jco")

# Confirm that we do not violate sample size assumption (the n in each cell > the number of outcome variables)
sperm.avg %>%
  group_by(Species) %>%
  get_summary_stats(Tail, Midpiece, Head, type = "mean_sd")

sperm.avg %>%
  group_by(Breeding) %>%
  get_summary_stats(Tail, Midpiece, Head, type = "mean_sd")

#Check for outliers
sperm.avg %>%
  group_by(Species) %>%
  identify_outliers(Tail)

# Univariate normality assumption
sperm.avg %>%
  group_by(Species) %>%
  shapiro_test(Tail, Midpiece, Head) %>%
  arrange(variable)

# QQ plot
ggqqplot(sperm.avg, "Tail", facet.by = "Species",
         ylab = " Tail", ggtheme = theme_bw())

ggqqplot(sperm.avg, "Midpiece", facet.by = "Species",
         ylab = " Midpiece", ggtheme = theme_bw())

ggqqplot(sperm.avg, "Head", facet.by = "Species",
         ylab = " Head", ggtheme = theme_bw())

# Linearity assumption
# Create a scatterplot matrix by group
library(GGally)
results <- sperm.avg %>%
  select(Tail, Midpiece, Head, Species) %>%
  group_by(Species) %>%
  doo(~ggpairs(.) + theme_bw(), result = "plots")
results
results$plots
# Looks like there are several non linear relationships - we will lose power

# Homogeneity of covariances
box_m(sperm.avg[, c("Tail", "Midpiece","Head")], sperm.avg$Species)
# statistic p.value parameter method                                             
# <dbl>   <dbl>     <dbl> <chr>                                              
#   1      4.18   0.652         6 Box's M-test for Homogeneity of Covariance Matrices
# Covariation is homogeneic!

# Homogeneity of variance
sperm.avg %>% 
  gather(key = "variable", value = "value", Tail, Midpiece, Head) %>%
  group_by(variable) %>%
  levene_test(value ~ Species)

# variable   df1   df2 statistic     p
# <chr>    <int> <int>     <dbl> <dbl>
#   1 Head         1    16     0.147 0.707
# 2 Midpiece     1    16     0.676 0.423
# 3 Tail         1    16     0.118 0.736


# MANOVA for species comparison

model <- lm(cbind(Tail, Midpiece, Head) ~ Species, sperm.avg)
Manova(model, test.statistic = "Pillai") # recommended for unbalanced design

# Group the data by variable
grouped.data <- sperm.avg %>%
  gather(key = "variable", value = "value", Tail, Midpiece, Head) %>%
  group_by(variable)
# Type II MANOVA Tests: Pillai test statistic
#         Df    test stat approx F num Df den Df   Pr(>F)    
# Species  1    0.8373   24.017      3     14   8.78e-06 ***

grouped.data %>% anova_test(value ~ Species)
# variable    Effect    DFn   DFd    F      p     `p<.05`   ges
# 1 Head     Species     1    16  0.195   0.665      ""     0.012
# 2 Midpiece Species     1    16  24.0   0.000162   "*"     0.6  
# 3 Tail     Species     1    16  48.1   0.00000336 "*"     0.75 

# Between species, there was a significant difference in Tail length (F(1, 16) = 48.1, p < 0.0001 ) 
# and midpiece length (F(1, 16) = 24.0, p = 0.00016 ), but not head length (F(1, 16) = 0.195, p = 0.665 )

# Bonferroni multiple testing correction: divide 0.05 by # tests (3), so significance criteria is p < 0.01666
p <- c(0.665,0.000162,0.00000336)
p.adjust (p, method = "bonferroni")
# 1.000e+00 4.860e-04 1.008e-05


# MANOVA for breeding stage comparison
model <- lm(cbind(Tail, Midpiece, Head) ~ Breeding, sperm.avg)
Manova(model, test.statistic = "Pillai") # recommended for unbalanced design
#           Df test stat approx F num Df den Df Pr(>F)
# Breeding  1   0.12001   0.6364      3     14 0.6039

# Group the data by variable
grouped.data <- sperm.avg %>%
  gather(key = "variable", value = "value", Tail, Midpiece, Head) %>%
  group_by(variable)

grouped.data %>% anova_test(value ~ Breeding)
# variable Effect     DFn   DFd     F     p `p<.05`   ges
# 1 Head     Breeding     1    16 1.37  0.26  ""      0.079
# 2 Midpiece Breeding     1    16 0.424 0.524 ""      0.026
# 3 Tail     Breeding     1    16 0.033 0.859 ""      0.002

p <- c(0.26,0.524,0.859)
p.adjust (p, method = "bonferroni")
#  0.78 1.00 1.00


# Avg Tail length  - species and breeding stage comparisons
ggboxplot(sperm, x = "Species", y = "Tail", add = "dotplot", color = "Individual")

tail <- ggboxplot(sperm.avg, x = "Species", y = "Tail", add = "dotplot")
tail + theme(text=element_text(size=rel(2.2))) + theme(axis.title.x = element_text(size=16))

shapiro.test(jacana.sperm.avg$Tail) # W = 0.89545, p-value = 0.3043
t.test(sperm.avg$Tail ~ sperm.avg$Species) #t = 6.3913, df = 9.7457, p-value = 8.89e-05

ggboxplot(sperm.avg, x = "Species", y = "Tail", color = "Breeding", add = "dotplot")

species.tail.length <- lm(Tail ~ Species + Breeding, data = sperm.avg)
anova(species.tail.length)
#           Df  Sum Sq Mean Sq F value    Pr(>F)    
#Species    1 225.046 225.046 49.2224 4.168e-06 ***
#Breeding   1   6.281   6.281  1.3739    0.2594    
#Residuals 15  68.580   4.572 


# Plots that go into Figure 2 - ppt version - Tail Length
tail.length <- ggboxplot(sperm.avg, x = "Species", y = "Tail", add = "dotplot")
editable_graph <- dml(ggobj = tail.length)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "tail.length.pptx")



# Avg Midpiece Length
ggboxplot(sperm, x = "Species", y = "Midpiece", add = "dotplot", color = "Individual")

midpiece <- ggboxplot(sperm.avg, x = "Species", y = "Midpiece", add = "dotplot")
midpiece + theme(text=element_text(size=rel(2.2))) + theme(axis.title.x = element_text(size=16))

shapiro.test(jacana.sperm.avg$Midpiece) # W = 0.91635, p-value = 0.4416
t.test(sperm.avg$Midpiece ~ sperm.avg$Species) #t = 4.7361, df = 11.546, p-value = 0.0005371

ggboxplot(sperm.avg, x = "Species", y = "Midpiece", color = "Breeding", add = "dotplot")

species.midpiece.length <- lm(Midpiece ~ Species + Breeding, data = sperm.avg)
anova(species.midpiece.length)
#         Df Sum Sq Mean Sq F value    Pr(>F)    
#Species    1 5.2092  5.2092 26.6509 0.0001158 ***
#Breeding   1 0.5456  0.5456  2.7914 0.1154960    
#Residuals 15 2.9319  0.1955     

# Plots that go into Figure 2 - ppt version - Midpiece Length
midpiece.length <- ggboxplot(sperm.avg, x = "Species", y = "Midpiece", add = "dotplot")
midpiece.length
editable_graph <- dml(ggobj = midpiece.length)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "midpiece.length.pptx")



# Head
ggboxplot(sperm, x = "Species", y = "Head", add = "dotplot", color = "Individual")

head <- ggboxplot(sperm.avg, x = "Species", y = "Head", add = "dotplot")
head + theme(text=element_text(size=rel(2.2))) + theme(axis.title.x = element_text(size=16))

shapiro.test(jacana.sperm.avg$Head) # W = 0.96137, p-value = 0.8304
t.test(sperm.avg$Head ~ sperm.avg$Species) # t = -0.42472, df = 11.303, p-value = 0.679

ggboxplot(sperm.avg, x = "Species", y = "Head", color = "Breeding", add = "dotplot")

species.head.length <- lm(Head ~ Species + Breeding, data = sperm.avg)
anova(species.head.length)
#Df Sum Sq Mean Sq F value Pr(>F)
#Species    1 0.0880 0.08804  0.2006 0.6607
#Breeding   1 0.6353 0.63534  1.4474 0.2476
#Residuals 15 6.5843 0.43895        

# Plots that go into Figure 2 - ppt version - Head Length
head.length <- ggboxplot(sperm.avg, x = "Species", y = "Head", add = "dotplot")
head.length
editable_graph <- dml(ggobj = head.length)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "head.length.pptx")


# Do testes differ between the species?
ggboxplot(sperm.avg, x = "Species", y = "Testes.Mass", add = "dotplot", color = "Breeding")

# ANCOVA
# Adequate sample size. Rule of thumb: the n in each cell > the number of outcome variables.
sperm.avg %>%
  group_by(Species) %>%
  get_summary_stats(Testes.Mass, type = "mean_sd")

sperm.avg %>%
  group_by(Breeding) %>%
  get_summary_stats(Testes.Mass, type = "mean_sd")


# Also need to test for homogeneity of variance
library(car)
leveneTest(sperm.avg$Testes.Mass ~ sperm.avg$Species)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
# group  1   0.7793 0.3904
#       16    
# Variance is equal. 

ancova_model_species <- aov(sperm.avg$Testes.Mass ~ sperm.avg$Species + sperm.avg$Somatic.Mass)
Anova(ancova_model_species, type = "III")
# Anova Table (Type III tests)
# 
# Response: sperm.avg$Testes.Mass
#                         Sum Sq Df F value Pr(>F)
# (Intercept)               460  1  0.0146 0.9055
# sperm.avg$Species        6270  1  0.1988 0.6621
# sperm.avg$Somatic.Mass  10206  1  0.3236 0.5779
# Residuals              473105 15               

ancova_model_breeding <- aov(sperm.avg$Testes.Mass ~ sperm.avg$Somatic.Mass + sperm.avg$Breeding)
Anova(ancova_model_breeding, type = "III")
# Anova Table (Type III tests)
# 
# Response: sperm.avg$Testes.Mass
# Sum Sq Df F value Pr(>F)
# (Intercept)             15937  1  0.5410 0.4734
# sperm.avg$Somatic.Mass    344  1  0.0117 0.9154
# sperm.avg$Breeding      37524  1  1.2739 0.2768
# Residuals              441850 15   



# Testes volume
ggboxplot(sperm.avg, x = "Species", y = "Avg.Testes.Vol", add = "dotplot", color = "Breeding")

species.testes.vol <- aov(Avg.Testes.Vol ~ Species * Breeding, data = sperm.avg,na.action= na.exclude)
summary(species.testes.vol)
# Df   Sum Sq Mean Sq F value Pr(>F)
# Species           1   464494  464494   0.401  0.541
# Breeding          1   789882  789882   0.682  0.428
# Species:Breeding  1  1733443 1733443   1.497  0.249
# Residuals        10 11579282 1157928               
# 4 observations deleted due to missingness

species.testes.vol <- aov(Avg.Testes.Vol ~ Species + Breeding, data = sperm.avg,na.action= na.exclude)
summary(species.testes.vol)
# Df   Sum Sq Mean Sq F value Pr(>F)
# Species      1   464494  464494   0.384  0.548
# Breeding     1   789882  789882   0.653  0.436
# Residuals   11 13312725 1210248               
# 4 observations deleted due to missingness 




# Coefficient of variation of each trait - variation across breeding season?
# Intermale CV for both breeding stages and combined

spinosa.sperm = subset(sperm, Species == "Northern")
jacana.sperm = subset(sperm, Species == "Wattled")
spinosa.sperm.court = subset(sperm, Species == "Northern" & Breeding == "Copulation")
spinosa.sperm.inc = subset(sperm, Species == "Northern" & Breeding == "Incubation")
jacana.sperm.court = subset(sperm, Species == "Wattled" & Breeding == "Copulation")
jacana.sperm.inc = subset(sperm, Species == "Wattled" & Breeding == "Incubation")

# Intermale CV tail-
sd(spinosa.sperm$Tail, na.rm=TRUE)/
  mean(spinosa.sperm$Tail, na.rm=TRUE)*100 # spinosa CV =  4.459243
sd(spinosa.sperm.court$Tail, na.rm=TRUE)/
  mean(spinosa.sperm.court$Tail, na.rm=TRUE)*100 # spinosa court CV = 3.677851
sd(spinosa.sperm.inc$Tail, na.rm=TRUE)/
  mean(spinosa.sperm.inc$Tail, na.rm=TRUE)*100 # spinosa inc CV = 5.030138

sd(jacana.sperm$Tail, na.rm=TRUE)/
  mean(jacana.sperm$Tail, na.rm=TRUE)*100 # spinosa CV =  5.746834
sd(jacana.sperm.court$Tail, na.rm=TRUE)/
  mean(jacana.sperm.court$Tail, na.rm=TRUE)*100 # spinosa court CV =  4.439126
sd(jacana.sperm.inc$Tail, na.rm=TRUE)/
  mean(jacana.sperm.inc$Tail, na.rm=TRUE)*100 # spinosa inc CV = 6.43182

# Intermale CV Midpiece
sd(spinosa.sperm$Midpiece, na.rm=TRUE)/
  mean(spinosa.sperm$Midpiece, na.rm=TRUE)*100 # spinosa CV =  8.875077
sd(spinosa.sperm.court$Midpiece, na.rm=TRUE)/
  mean(spinosa.sperm.court$Midpiece, na.rm=TRUE)*100 # spinosa court CV = 7.498641
sd(spinosa.sperm.inc$Midpiece, na.rm=TRUE)/
  mean(spinosa.sperm.inc$Midpiece, na.rm=TRUE)*100 # spinosa inc CV = 9.705709

sd(jacana.sperm$Midpiece, na.rm=TRUE)/
  mean(jacana.sperm$Midpiece, na.rm=TRUE)*100 # spinosa CV =   8.065452
sd(jacana.sperm.court$Midpiece, na.rm=TRUE)/
  mean(jacana.sperm.court$Midpiece, na.rm=TRUE)*100 # spinosa court CV =  6.628148
sd(jacana.sperm.inc$Midpiece, na.rm=TRUE)/
  mean(jacana.sperm.inc$Midpiece, na.rm=TRUE)*100 # spinosa inc CV = 9.632959

# Intermale CV Head
sd(spinosa.sperm$Head, na.rm=TRUE)/
  mean(spinosa.sperm$Head, na.rm=TRUE)*100 # spinosa CV =  8.43109
sd(spinosa.sperm.court$Head, na.rm=TRUE)/
  mean(spinosa.sperm.court$Head, na.rm=TRUE)*100 # spinosa court CV = 8.785489
sd(spinosa.sperm.inc$Head, na.rm=TRUE)/
  mean(spinosa.sperm.inc$Head, na.rm=TRUE)*100 # spinosa inc CV = 7.77009

sd(jacana.sperm$Head, na.rm=TRUE)/
  mean(jacana.sperm$Head, na.rm=TRUE)*100 # spinosa CV =   9.335285
sd(jacana.sperm.court$Head, na.rm=TRUE)/
  mean(jacana.sperm.court$Head, na.rm=TRUE)*100 # spinosa court CV =  9.366753
sd(jacana.sperm.inc$Head, na.rm=TRUE)/
  mean(jacana.sperm.inc$Head, na.rm=TRUE)*100 # spinosa inc CV = 9.343497



# CV equality - test for difference in intraspecific + CVs between northern and wattled

#install.packages("cvequality")
library(cvequality)
library(ggbeeswarm)
library(ggplot2)
library(knitr)

kable(head(sperm.avg), caption = "Preview of first few rows of the sperm data")
spinosa.sperm = subset(sperm, Species == "Northern")
jacana.sperm = subset(sperm, Species == "Wattled")

sperm.court = subset(sperm.avg, sperm.avg$Breeding == "Copulating")
sperm.inc = subset(sperm.avg, sperm.avg$Breeding == "Incubating")

# Supplementary Table 1
# Copulation - species diffs
ggplot(sperm.court, aes(Species, Tail)) + geom_boxplot() + geom_quasirandom(alpha = 0.05) + theme_bw()
species_court_tail_length_cv_test_MSLR <- with(sperm.court, mslr_test(nr = 1e4,Tail, Species))
species_court_tail_length_cv_test_MSLR # MSLRT = 0.4430512, p =  0.5056534

ggplot(sperm.court, aes(Species, Midpiece)) + geom_boxplot() + geom_quasirandom(alpha = 0.05) + theme_bw()
species_court_midpiece_length_cv_test_MSLR <- with(sperm.court, mslr_test(nr = 1e4,Midpiece, Species))
species_court_midpiece_length_cv_test_MSLR # MSLRT = 0.02178721, p =  0.8826546

ggplot(sperm.court, aes(Species, Head)) + geom_boxplot() + geom_quasirandom(alpha = 0.05) + theme_bw()
species_court_head_length_cv_test_MSLR <- with(sperm.court, mslr_test(nr = 1e4,Head, Species))
species_court_head_length_cv_test_MSLR # MSLRT = 0.2724925, p =  0.6016646

# Incubation- species diffs
ggplot(sperm.inc, aes(Species, Tail)) + geom_boxplot() + geom_quasirandom(alpha = 0.05) + theme_bw()
species_inc_tail_length_cv_test_MSLR <- with(sperm.inc, mslr_test(nr = 1e4,Tail, Species))
species_inc_tail_length_cv_test_MSLR # MSLRT = 0.7983578, p =  0.3715848

ggplot(sperm.inc, aes(Species, Midpiece)) + geom_boxplot() + geom_quasirandom(alpha = 0.05) + theme_bw()
species_inc_midpiece_length_cv_test_MSLR <- with(sperm.inc, mslr_test(nr = 1e4,Midpiece, Species))
species_inc_midpiece_length_cv_test_MSLR # MSLRT = 0.4125012, p =  0.5207027

ggplot(sperm.inc, aes(Species, Head)) + geom_boxplot() + geom_quasirandom(alpha = 0.05) + theme_bw()
species_inc_head_length_cv_test_MSLR <- with(sperm.inc, mslr_test(nr = 1e4,Head, Species))
species_inc_head_length_cv_test_MSLR # MSLRT = 1.44006, p =  0.2301296

# Bonferroni multiple testing correction: 
p <- c(0.5056534,0.8826546,0.6016646,0.3715848,0.5207027,0.2301296)
p.adjust (p, method = "bonferroni")


# Supplementary Table 2
# test of difference in CV in breeding stage, for each species
# Northern Jacana
ggplot(spinosa.sperm, aes(Breeding, Tail)) + geom_boxplot() + geom_quasirandom(alpha = 0.05) + theme_bw()
#spinosa_breeding_tail_length_cv_test <-  with(spinosa.sperm, asymptotic_test(Tail,Breeding))
spinosa_breeding_tail_length_cv_test_MSLRT <- with(spinosa.sperm, mslr_test(nr = 1e4,Tail, Breeding))
spinosa_breeding_tail_length_cv_test_MSLRT # MSLRT = 5.273314, p =  0.02463542

ggplot(spinosa.sperm, aes(Breeding, Midpiece)) + geom_boxplot() + geom_quasirandom(alpha = 0.05) + theme_bw()
#breeding_midpiece_length_cv_test <-  with(spinosa.sperm, asymptotic_test(Midpiece,Breeding))
breeding_midpiece_length_cv_test_MSLRT <- with(spinosa.sperm, mslr_test(nr = 1e4,Midpiece, Breeding))
breeding_midpiece_length_cv_test_MSLRT # MSLRT = 3.444121, p =  0.06347718

ggplot(spinosa.sperm, aes(Breeding, Head)) + geom_boxplot() + geom_quasirandom(alpha = 0.05) + theme_bw()
#breeding_midpiece_length_cv_test <-  with(spinosa.sperm, asymptotic_test(Head,Breeding))
breeding_head_length_cv_test_MSLRT <- with(spinosa.sperm, mslr_test(nr = 1e4,Head, Breeding))
breeding_head_length_cv_test_MSLRT # MSLRT = 0.7968582, p = 0.3720344

# Wattled Jacana
ggplot(jacana.sperm, aes(Breeding, Tail)) + geom_boxplot() + geom_quasirandom(alpha = 0.05) + theme_bw()
#breeding_tail_length_cv_test <-  with(jacana.sperm, asymptotic_test(Tail,Breeding))
breeding_tail_length_cv_test_MSLRT <- with(jacana.sperm, mslr_test(nr = 1e4,Tail, Breeding))
breeding_tail_length_cv_test_MSLRT # MSLRT = 4.519644, p = 0.03350781

ggplot(jacana.sperm, aes(Breeding, Midpiece)) + geom_boxplot() + geom_quasirandom(alpha = 0.05) + theme_bw()
#breeding_midpiece_length_cv_test <-  with(jacana.sperm, asymptotic_test(Midpiece,Breeding))
breeding_midpiece_length_cv_test_MSLRT <- with(jacana.sperm, mslr_test(nr = 1e4,Midpiece, Breeding))
breeding_midpiece_length_cv_test_MSLRT # MSLRT = 4.572663, p = 0.03248604

ggplot(jacana.sperm, aes(Breeding, Head)) + geom_boxplot() + geom_quasirandom(alpha = 0.05) + theme_bw()
#breeding_midpiece_length_cv_test <-  with(jacana.sperm, asymptotic_test(Midpiece,Breeding))
breeding_head_length_cv_test_MSLRT <- with(jacana.sperm, mslr_test(nr = 1e4,Head, Breeding))
breeding_head_length_cv_test_MSLRT # MSLRT = 0.0009138081, p =  0.9758842

# Bonferroni multiple testing correction: 
p <- c(0.02463542,0.06347718,0.3720344,0.03350781,0.03248604,0.9758842)
p.adjust (p, method = "bonferroni")
# 0.1478125 0.3808631 1.0000000 0.2010469 0.1949162 1.0000000

# Intra-ejaculate Coefficient of Variation

# Tail Intra-ejaculate CV 
ggboxplot(sperm.avg, x = "Species", y = "CV.Tail", add = "dotplot", color = "Breeding")

CV.Tail <- ggboxplot(sperm.avg, x = "Species", y = "CV.Tail", add = "dotplot", color = "Breeding", ylim = c(0,15))
CV.Tail
editable_graph <- dml(ggobj = CV.Tail)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "CV.Tail.pptx")
# CV for tail length changes with breeding stage, marginally different between species.
# Tail length CV is lower during copulation than incubation

species.CV.tail <- aov(CV.Tail ~ Species * Breeding, data = sperm.avg)
summary(species.CV.tail)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# Species           1  4.338   4.338   4.294 0.05719 . 
# Breeding          1 12.131  12.131  12.010 0.00379 **
# Species:Breeding  1  0.565   0.565   0.559 0.46693   
# Residuals        14 14.141   1.010 

# Interaction not significant so we deleted

species.CV.tail <- aov(CV.Tail ~ Species + Breeding, data = sperm.avg)
summary(species.CV.tail)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# Species      1  4.338   4.338   4.424 0.05272 . 
# Breeding     1 12.131  12.131  12.373 0.00311 **
# Residuals   15 14.706   0.980    


# Midpiece Intra-ejaculate CV 
ggboxplot(sperm.avg, x = "Species", y = "CV.Midpiece", add = "dotplot", color = "Breeding")

species.CV.midpiece <- aov(CV.Midpiece ~ Species * Breeding, data = sperm.avg)
summary(species.CV.midpiece)
# Df Sum Sq Mean Sq F value Pr(>F)
# Species           1  10.83  10.826   1.701  0.213
# Breeding          1  11.96  11.957   1.879  0.192
# Species:Breeding  1   0.04   0.039   0.006  0.939
# Residuals        14  89.08   6.363 

species.CV.midpiece <- aov(CV.Midpiece ~ Species + Breeding, data = sperm.avg)
summary(species.CV.midpiece)
#          Df Sum Sq Mean Sq F value Pr(>F)
# Species      1  10.83  10.826   1.822  0.197
# Breeding     1  11.96  11.957   2.012  0.176
# Residuals   15  89.12   5.941  

CV.midpiece <- ggboxplot(sperm.avg, x = "Species", y = "CV.Midpiece", add = "dotplot", color = "Breeding", ylim = c(0,15))
CV.midpiece
editable_graph <- dml(ggobj = CV.midpiece)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "CV.midpiece.pptx")
# CV for midpiece length does not significantly differ w/ breeding stage or species


# Head Intra-ejaculate CV 
species.CV.head <- aov(CV.Head ~ Species * Breeding, data = sperm.avg)
summary(species.CV.head)
# Df Sum Sq Mean Sq F value Pr(>F)
# Species           1   2.03   2.027   0.502  0.490
# Breeding          1   0.91   0.913   0.226  0.642
# Species:Breeding  1   9.12   9.125   2.261  0.155
# Residuals        14  56.50   4.036

species.CV.head <- aov(CV.Head ~ Species + Breeding, data = sperm.avg)
summary(species.CV.head)
#           Df Sum Sq Mean Sq F value Pr(>F)
# Species      1   2.03   2.027   0.463  0.506
# Breeding     1   0.91   0.913   0.209  0.654
# Residuals   15  65.62   4.375  

ggboxplot(sperm.avg, x = "Species", y = "CV.Head", add = "dotplot", color = "Breeding",ylim=c(0,14))
CV.Head <- ggboxplot(sperm.avg, x = "Species", y = "CV.Head", add = "dotplot", color = "Breeding",ylim=c(0,15))
editable_graph <- dml(ggobj = CV.Head)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "CV.Head.pptx")
# CV for head length does not significantly differ w/ breeding stage or species
CV.Head


