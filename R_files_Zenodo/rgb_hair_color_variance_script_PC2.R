RGB_traits_hind_scores <- read.csv("RGB_traits_hind_scores.csv")
#OR
RGB_traits_fore_scores <- read.csv("RGB_traits_fore_scores.csv")
#OR
RGB_traits_head_scores <- read.csv("RGB_traits_head_scores.csv")
#OR
RGB_traits_tail_scores <- read.csv("RGB_traits_tail_scores.csv")
#OR
RGB_traits_torso_scores <- read.csv("RGB_traits_torso_scores.csv")

RGB_traits_hind_scores_cats = subset(RGB_traits_hind_scores, Clade =="Catarrhine", select=PC2)
RGB_traits_hind_scores_strep = subset(RGB_traits_hind_scores, Clade =="Strepsirhine", select=PC2)
RGB_traits_hind_scores_plat = subset(RGB_traits_hind_scores, Clade =="Platyrrhine", select=PC2)
#OR
RGB_traits_fore_scores_cats = subset(RGB_traits_fore_scores, Clade =="Catarrhine", select=PC2)
RGB_traits_fore_scores_strep = subset(RGB_traits_fore_scores, Clade =="Strepsirhine", select=PC2)
RGB_traits_fore_scores_plat = subset(RGB_traits_fore_scores, Clade =="Platyrrhine", select=PC2)
#OR
RGB_traits_head_scores_cats = subset(RGB_traits_head_scores, Clade =="Catarrhine", select=PC2)
RGB_traits_head_scores_strep = subset(RGB_traits_head_scores, Clade =="Strepsirhine", select=PC2)
RGB_traits_head_scores_plat = subset(RGB_traits_head_scores, Clade =="Platyrrhine", select=PC2)
#OR
RGB_traits_tail_scores_cats = subset(RGB_traits_tail_scores, Clade =="Catarrhine", select=PC2)
RGB_traits_tail_scores_strep = subset(RGB_traits_tail_scores, Clade =="Strepsirhine", select=PC2)
RGB_traits_tail_scores_plat = subset(RGB_traits_tail_scores, Clade =="Platyrrhine", select=PC2)
#OR
RGB_traits_torso_scores_cats = subset(RGB_traits_torso_scores, Clade =="Catarrhine", select=PC2)
RGB_traits_torso_scores_strep = subset(RGB_traits_torso_scores, Clade =="Strepsirhine", select=PC2)
RGB_traits_torso_scores_plat = subset(RGB_traits_torso_scores, Clade =="Platyrrhine", select=PC2)


View(RGB_traits_hind_scores_strep)
View(RGB_traits_hind_scores_cats)
View(RGB_traits_hind_scores_plat)
#OR
View(RGB_traits_fore_scores_strep)
View(RGB_traits_fore_scores_cats)
View(RGB_traits_fore_scores_plat)
#OR
View(RGB_traits_head_scores_strep)
View(RGB_traits_head_scores_cats)
View(RGB_traits_head_scores_plat)
#OR
View(RGB_traits_tail_scores_strep)
View(RGB_traits_tail_scores_cats)
View(RGB_traits_tail_scores_plat)
#OR
View(RGB_traits_torso_scores_strep)
View(RGB_traits_torso_scores_cats)
View(RGB_traits_torso_scores_plat)

# For catarrhines and platyrrhines, sample n times (when n=number of strep species in body region pc scores file) without replacement and calculate the variance for each subset, then create a histogram of these simulated variances
set.seed(4324)
hind_color_cats = replicate (1000, var(sample(RGB_traits_hind_scores_cats$PC2, size = 23, replace =F)))
hind_color_plat = replicate (1000, var(sample(RGB_traits_hind_scores_plat$PC2, size = 23, replace =F)))
#OR
fore_color_cats = replicate (1000, var(sample(RGB_traits_fore_scores_cats$PC2, size = 23, replace =F)))
fore_color_plat = replicate (1000, var(sample(RGB_traits_fore_scores_plat$PC2, size = 23, replace =F)))
#OR
head_color_cats = replicate (1000, var(sample(RGB_traits_head_scores_cats$PC2, size = 24, replace =F)))
head_color_plat = replicate (1000, var(sample(RGB_traits_head_scores_plat$PC2, size = 24, replace =F)))
#OR
tail_color_cats = replicate (1000, var(sample(RGB_traits_tail_scores_cats$PC2, size = 23, replace =F)))
tail_color_plat = replicate (1000, var(sample(RGB_traits_tail_scores_plat$PC2, size = 23, replace =F)))
#OR
torso_color_cats = replicate (1000, var(sample(RGB_traits_torso_scores_cats$PC2, size = 24, replace =F)))
torso_color_plat = replicate (1000, var(sample(RGB_traits_torso_scores_plat$PC2, size = 24, replace =F)))

#make a histogram
library(ggplot2)

hist_hind_color_cats = ggplot() + geom_histogram(aes(hind_color_cats), binwidth = 0.05)
hist_hind_color_plat = ggplot() + geom_histogram(aes(hind_color_plat), binwidth = 0.05)
#OR
hist_fore_color_cats = ggplot() + geom_histogram(aes(fore_color_cats), binwidth = 0.05)
hist_fore_color_plat = ggplot() + geom_histogram(aes(fore_color_plat), binwidth = 0.05)
#OR
hist_head_color_cats = ggplot() + geom_histogram(aes(head_color_cats), binwidth = 0.05)
hist_head_color_plat = ggplot() + geom_histogram(aes(head_color_plat), binwidth = 0.05)
#OR
hist_tail_color_cats = ggplot() + geom_histogram(aes(tail_color_cats), binwidth = 0.05)
hist_tail_color_plat = ggplot() + geom_histogram(aes(tail_color_plat), binwidth = 0.05)
#OR
hist_torso_color_cats = ggplot() + geom_histogram(aes(torso_color_cats), binwidth = 0.05)
hist_torso_color_plat = ggplot() + geom_histogram(aes(torso_color_plat), binwidth = 0.05)

#Calculate the real strep variance
RGB_traits_hind_scores_strep = subset(RGB_traits_hind_scores, Clade =="Strepsirhine", select=PC2)
strep_var1 = var(RGB_traits_hind_scores_strep)
#OR
strep_var2 = var(RGB_traits_fore_scores_strep)
#OR
strep_var3 = var(RGB_traits_head_scores_strep)
#OR
strep_var4 = var(RGB_traits_tail_scores_strep)
#OR
strep_var5 = var(RGB_traits_torso_scores_strep)

#plot the real strep variance as a vertical line; if it's outside of 95% distribution of random cata values then they are significantly different
hist_hind_color_cats+geom_vline(xintercept = strep_var1, color ="blue")
hist_hind_color_plat+geom_vline(xintercept = strep_var1, color ="blue")  #This one narrowly falls within distribution (the upper CI = 1.209, strep var = 1.107)
#OR
hist_fore_color_cats+geom_vline(xintercept = strep_var2, color ="blue")
hist_fore_color_plat+geom_vline(xintercept = strep_var2, color ="blue")
#OR
hist_head_color_cats+geom_vline(xintercept = strep_var3, color ="blue")
hist_head_color_plat+geom_vline(xintercept = strep_var3, color ="blue")
#OR
hist_tail_color_cats+geom_vline(xintercept = strep_var4, color ="blue")  #This one falls well within range (lower CI=0.402, strep var=0.469)
hist_tail_color_plat+geom_vline(xintercept = strep_var4, color ="blue")  #This one falls well within range (lower CI=0.342, strep var=0.469)
#OR
hist_torso_color_cats+geom_vline(xintercept = strep_var5, color ="blue") #This one falls within range (lower CI=1.414, strep var=1.56)
hist_torso_color_plat+geom_vline(xintercept = strep_var5, color ="blue") #This one falls well within range (upper CI=2.064, strep var=1.56)

###Calculate 95% CI
#example:
quantile(torso_color_plat,.95)
