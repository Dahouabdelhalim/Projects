#loading packages for data analysis and plotting (Figure 1 and 2)
library(tidyverse)
library(ggbeeswarm) 
library(ggthemes)
library(ggpubr) 
library(dunn.test)

# reading Source_data
GSF<-read.table("Source_data.csv", header = T, sep=",") 
names(GSF)
# -------------------------------- Figure 2 ---------------------------------------

GSF_climate <- GSF %>% 
  filter(!(KG_climate == "E")) %>% 
  drop_na(Stemflow, KG_climate)

# Test of normality for global stemflow
# P<0.05 indicates the data is not in normal distribution; P>0.05 indicates that the data obey a normal distribution
shapiro.test(GSF$Stemflow)

# Spread data according to KG_climate (koppen-Geier Climate)
SF_Climate<- spread(GSF, key = KG_climate, value = Stemflow)

# test of normality of stemflow in each climate group
shapiro.test(SF_Climate$A)
shapiro.test(SF_Climate$B)
shapiro.test(SF_Climate$C)
shapiro.test(SF_Climate$D) 

# Set comparisons which will be used in multiple comparisons in Figure 2
my_comparisons <- list(c("A", "B"),
                       c("A", "C"), 
                       c("A", "D"),
                      c("B", "C"),
                      c("B", "D"),
                      c("C", "D"))

# Kruskal-Wallis Test
kruskal.test(Stemflow ~ KG_climate, data = GSF)

# dunn.test for multiple comparisons
dunn.test(x = GSF$Stemflow, g = GSF$KG_climate, method = "bonferroni", kw = TRUE, label = TRUE, altp = TRUE)

# Plotting Figure 2; differences in Semflow between climate types was compared using Kruskal-Wallis rank sum test 
SF_KG <- ggplot(GSF_climate, aes(x = reorder(KG_climate, Stemflow, FUN = mean), y = Stemflow, na.rm = TRUE, fill = KG_climate))+ #reorder(x,y,FUN) is used to reorder x according to mean of y
  geom_violin(alpha = 0.3, 
              scale = "area")+
  geom_quasirandom(alpha = 0.5, size = 3,
                   varwidth = TRUE,
                   aes(color = KG_climate)) + 
  geom_boxplot(width = .1,
               notch = TRUE,
               fill = "orange",
               outlier.shape = NA,
               alpha = 0.5) +
  stat_compare_means(comparisons = my_comparisons, step.increase = 0.05, hide.ns = FALSE, vjust = 0.5, size = 6, tip.length = 0.02, bracket.size = 0.5, label = "p.signif", label.y = c(47, 43, 45, 45, 49, 41))+
  labs(x = "Climate type", y = "Stemflow percentage (%)") + scale_y_continuous(limits = c(0, 50))+ 
  theme_wsj() + 
  theme(axis.title.x = element_text(vjust = -0.5), axis.title.y = element_text(vjust = 2))+
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 30), legend.position = "none")

# print Figure 2
print(SF_KG)

# - ------------------------------ Figure 3 -----------------------------------------#
# Spread data according to community type (tree and shrub)
SF_community_type<- spread(GSF, key = Community_type, value = Stemflow)

# test of normality of stemflow in Tree and Shrub
shapiro.test(SF_community_type$Tree)
shapiro.test(SF_community_type$Shrub)

# Kruskal-Wallis Test
kruskal.test(Stemflow ~ Community_type, data = GSF)

# dunn.test for multiple comparisons
dunn.test(x = GSF$Stemflow, g = GSF$Community_type, method = "bonferroni", kw = TRUE, label = TRUE, altp = TRUE)

# Plotting Figure 3
SF_Tree_Shrub <- ggplot(data = GSF, aes(x = Community_type, y = Stemflow, fill = Community_type)) +
  geom_violin(alpha = 0.3) +
  geom_quasirandom(alpha = 0.6, size = 3,
                   varwidth = TRUE,
                   aes(color = Community_type)) +
  geom_boxplot(width = .15,
               notch = TRUE,
               fill = "orange",
               outlier.shape = NA,
               alpha = 0.5)+
  stat_compare_means(comparisons = list(c("Shrub", "Tree")), step.increase = 0.05, hide.ns = FALSE, vjust = 0.5, size = 6, tip.length = 0.02, bracket.size = 0.5, label = "p.signif", label.y = c(51))+ 
  labs(x = "Community type", y = "Stemflow percentage (%)") +
  theme_wsj() +
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 30), legend.position = "none")

# print Figure 3
print(SF_Tree_Shrub)
