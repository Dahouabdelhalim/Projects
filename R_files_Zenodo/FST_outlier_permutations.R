##############################################################
####script modified from Rennison et al. (2020) "Ecological factors and morphological traits are associated with repeated genomic differentiation between lake and stream stickleback"

library(plyr)
library(dplyr)
library(ggplot2)

#permutations of the Cuba data
FST_Cuba_5percent <- read.csv("FST_outlier_permutations/Cuba_set_all_markers.with_outlier_info.txt", sep = ' ')

rowstat <- list()
perm <- 10000

for(i in 1:perm){
  y <- as.list(FST_Cuba_5percent[,(3:6)]) 
  y1 <- lapply(y, function(y){
    y1 <- y
    y1 <- sample(y1, replace=FALSE)
    return(y1)
  })
  y2 <- data.frame(y1)
  y2$overalltrue <- rowSums(y2 == "TRUE")
  rowstat[[i]] <- y2$overalltrue
}

rowstat_Cuba_5percent <- as.data.frame(rowstat)

perm_2_repeated_counts_Cuba_5percent<-as.data.frame(sapply(rowstat_Cuba_5percent[names(rowstat_Cuba_5percent)], function(x) sum(x>1)))
perm_3_repeated_counts_Cuba_5percent<-as.data.frame(sapply(rowstat_Cuba_5percent[names(rowstat_Cuba_5percent)], function(x) sum(x>2)))

##rename columns 
names(perm_2_repeated_counts_Cuba_5percent)[names(perm_2_repeated_counts_Cuba_5percent) == "sapply(rowstat_Cuba_5percent[names(rowstat_Cuba_5percent)], function(x) sum(x > 1))"] <- "outlier_counts"
names(perm_3_repeated_counts_Cuba_5percent)[names(perm_3_repeated_counts_Cuba_5percent) == "sapply(rowstat_Cuba_5percent[names(rowstat_Cuba_5percent)], function(x) sum(x > 2))"] <- "outlier_counts"

#empirical estimates for outliers shared between more than 1 or more than 2 population pairs
empirical_2_repeated_counts_Cuba_5percent<-sum(FST_Cuba_5percent$outlier_count > 1)
empirical_3_repeated_counts_Cuba_5percent<-sum(FST_Cuba_5percent$outlier_count > 2)

#calculate significance
"see https://genomicsclass.github.io/book/pages/permutation_tests.html"
pval_Cuba_2_5percent<-(sum(perm_2_repeated_counts_Cuba_5percent$outlier_counts > empirical_2_repeated_counts_Cuba_5percent) + 1) / (length(perm_2_repeated_counts_Cuba_5percent$outlier_counts) + 1)
pval_Cuba_3_5percent<-(sum(perm_3_repeated_counts_Cuba_5percent$outlier_counts > empirical_3_repeated_counts_Cuba_5percent) + 1) / (length(perm_3_repeated_counts_Cuba_5percent$outlier_counts) + 1)

#plot results (permutations and empirical)
ggplot(data = perm_2_repeated_counts_Cuba_5percent, aes(x=outlier_counts)) +
  geom_density(adjust = 3,colour="black", fill="grey97")+
  geom_vline(xintercept=empirical_2_repeated_counts_Cuba_5percent,color="red", size=0.5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlab("\\nNumber of shared FST outliers")

ggplot(data = perm_3_repeated_counts_Cuba_5percent, aes(x=outlier_counts)) +
  geom_density(adjust = 3,colour="black", fill="grey97")+
  geom_vline(xintercept=empirical_3_repeated_counts_Cuba_5percent,color="red", size=0.5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlab("\\nNumber of shared FST outliers")

#repeated outlier counts based on the data and expected counts based on permutations (outliers repeated in 2 or more populations, and 3 or more populations)
#these counts are used for barplots below

empirical_2_repeated_counts_Cuba_5percent
empirical_3_repeated_counts_Cuba_5percent

median(perm_2_repeated_counts_Cuba_5percent$outlier_counts)
median(perm_3_repeated_counts_Cuba_5percent$outlier_counts)

#P values
pval_Cuba_2_5percent
pval_Cuba_3_5percent


#permutations of the Invasive set1 (first 4 random population pairs)
FST_set1_5percent <- read.csv("FST_outlier_permutations/Florida_set1_all_markers.with_outlier_info.txt", sep = ' ')

rowstat <- list()
perm <- 10000

for(i in 1:perm){
  y <- as.list(FST_set1_5percent[,(3:6)]) 
  y1 <- lapply(y, function(y){
    y1 <- y
    y1 <- sample(y1, replace=FALSE)
    return(y1)
  })
  y2 <- data.frame(y1)
  y2$overalltrue <- rowSums(y2 == "TRUE")
  rowstat[[i]] <- y2$overalltrue
}

rowstat_set1_5percent <- as.data.frame(rowstat)

perm_2_repeated_counts_set1_5percent<-as.data.frame(sapply(rowstat_set1_5percent[names(rowstat_set1_5percent)], function(x) sum(x>1)))
perm_3_repeated_counts_set1_5percent<-as.data.frame(sapply(rowstat_set1_5percent[names(rowstat_set1_5percent)], function(x) sum(x>2)))

##rename columns 
names(perm_2_repeated_counts_set1_5percent)[names(perm_2_repeated_counts_set1_5percent) == "sapply(rowstat_set1_5percent[names(rowstat_set1_5percent)], function(x) sum(x > 1))"] <- "outlier_counts"
names(perm_3_repeated_counts_set1_5percent)[names(perm_3_repeated_counts_set1_5percent) == "sapply(rowstat_set1_5percent[names(rowstat_set1_5percent)], function(x) sum(x > 2))"] <- "outlier_counts"

#empirical estimates for outliers shared between more than 1 or more than 2 population pairs
empirical_2_repeated_counts_set1_5percent<-sum(FST_set1_5percent$outlier_count > 1)
empirical_3_repeated_counts_set1_5percent<-sum(FST_set1_5percent$outlier_count > 2)

#calculate significance
"see https://genomicsclass.github.io/book/pages/permutation_tests.html"
pval_set1_2_5percent<-(sum(perm_2_repeated_counts_set1_5percent$outlier_counts > empirical_2_repeated_counts_set1_5percent) + 1) / (length(perm_2_repeated_counts_set1_5percent$outlier_counts) + 1)
pval_set1_3_5percent<-(sum(perm_3_repeated_counts_set1_5percent$outlier_counts > empirical_3_repeated_counts_set1_5percent) + 1) / (length(perm_3_repeated_counts_set1_5percent$outlier_counts) + 1)


#plot results (permutations and empirical)
ggplot(data = perm_2_repeated_counts_set1_5percent, aes(x=outlier_counts)) +
  geom_density(adjust = 3,colour="black", fill="grey97")+
  geom_vline(xintercept=empirical_2_repeated_counts_set1_5percent,color="red", size=0.5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlab("\\nNumber of shared FST outliers")

ggplot(data = perm_3_repeated_counts_set1_5percent, aes(x=outlier_counts)) +
  geom_density(adjust = 3,colour="black", fill="grey97")+
  geom_vline(xintercept=empirical_3_repeated_counts_set1_5percent,color="red", size=0.5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlab("\\nNumber of shared FST outliers")

#repeated outlier counts based on the data and expected counts based on permutations (outliers repeated in 2 or more populations, and 3 or more populations)
#these counts are used for barplots below

empirical_2_repeated_counts_set1_5percent
empirical_3_repeated_counts_set1_5percent

median(perm_2_repeated_counts_set1_5percent$outlier_counts)
median(perm_3_repeated_counts_set1_5percent$outlier_counts)

pval_set1_2_5percent
pval_set1_3_5percent


#permutations of the Invasive set2 (second 4 random population pairs)
FST_set2_5percent <- read.csv("FST_outlier_permutations/Florida_set2_all_markers.with_outlier_info.txt", sep = ' ')

rowstat <- list()
perm <- 10000

for(i in 1:perm){
  y <- as.list(FST_set2_5percent[,(3:6)]) 
  y1 <- lapply(y, function(y){
    y1 <- y
    y1 <- sample(y1, replace=FALSE)
    return(y1)
  })
  y2 <- data.frame(y1)
  y2$overalltrue <- rowSums(y2 == "TRUE")
  rowstat[[i]] <- y2$overalltrue
}

rowstat_set2_5percent <- as.data.frame(rowstat)

perm_2_repeated_counts_set2_5percent<-as.data.frame(sapply(rowstat_set2_5percent[names(rowstat_set2_5percent)], function(x) sum(x>1)))
perm_3_repeated_counts_set2_5percent<-as.data.frame(sapply(rowstat_set2_5percent[names(rowstat_set2_5percent)], function(x) sum(x>2)))

##rename columns 
names(perm_2_repeated_counts_set2_5percent)[names(perm_2_repeated_counts_set2_5percent) == "sapply(rowstat_set2_5percent[names(rowstat_set2_5percent)], function(x) sum(x > 1))"] <- "outlier_counts"
names(perm_3_repeated_counts_set2_5percent)[names(perm_3_repeated_counts_set2_5percent) == "sapply(rowstat_set2_5percent[names(rowstat_set2_5percent)], function(x) sum(x > 2))"] <- "outlier_counts"


#empirical estimates for outliers shared between more than 1 or more than 2 population pairs
empirical_2_repeated_counts_set2_5percent<-sum(FST_set2_5percent$outlier_count > 1)
empirical_3_repeated_counts_set2_5percent<-sum(FST_set2_5percent$outlier_count > 2)

#calculate significance
"see https://genomicsclass.github.io/book/pages/permutation_tests.html"
pval_set2_2_5percent<-(sum(perm_2_repeated_counts_set2_5percent$outlier_counts > empirical_2_repeated_counts_set2_5percent) + 1) / (length(perm_2_repeated_counts_set2_5percent$outlier_counts) + 1)
pval_set2_3_5percent<-(sum(perm_3_repeated_counts_set2_5percent$outlier_counts > empirical_3_repeated_counts_set2_5percent) + 1) / (length(perm_3_repeated_counts_set2_5percent$outlier_counts) + 1)


#plot results (permutations and empirical)
ggplot(data = perm_2_repeated_counts_set2_5percent, aes(x=outlier_counts)) +
  geom_density(adjust = 3,colour="black", fill="grey97")+
  geom_vline(xintercept=empirical_2_repeated_counts_set2_5percent,color="red", size=0.5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlab("\\nNumber of shared FST outliers")

ggplot(data = perm_3_repeated_counts_set2_5percent, aes(x=outlier_counts)) +
  geom_density(adjust = 3,colour="black", fill="grey97")+
  geom_vline(xintercept=empirical_3_repeated_counts_set2_5percent,color="red", size=0.5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlab("\\nNumber of shared FST outliers")

#repeated outlier counts based on the data and expected counts based on permutations (outliers repeated in 2 or more populations, and 3 or more populations)
#these counts are used for barplots below

empirical_2_repeated_counts_set2_5percent
empirical_3_repeated_counts_set2_5percent

median(perm_2_repeated_counts_set2_5percent$outlier_counts)
median(perm_3_repeated_counts_set2_5percent$outlier_counts)

pval_set2_2_5percent
pval_set2_3_5percent


#permutations of the Invasive set3 (third 4 random population pairs)

FST_set3_5percent <- read.csv("FST_outlier_permutations/Florida_set3_all_markers.with_outlier_info.txt", sep = ' ')

rowstat <- list()
perm <- 10000

for(i in 1:perm){
  y <- as.list(FST_set3_5percent[,(3:6)]) 
  y1 <- lapply(y, function(y){
    y1 <- y
    y1 <- sample(y1, replace=FALSE)
    return(y1)
  })
  y2 <- data.frame(y1)
  y2$overalltrue <- rowSums(y2 == "TRUE")
  rowstat[[i]] <- y2$overalltrue
}

rowstat_set3_5percent <- as.data.frame(rowstat)

perm_2_repeated_counts_set3_5percent<-as.data.frame(sapply(rowstat_set3_5percent[names(rowstat_set3_5percent)], function(x) sum(x>1)))
perm_3_repeated_counts_set3_5percent<-as.data.frame(sapply(rowstat_set3_5percent[names(rowstat_set3_5percent)], function(x) sum(x>2)))

##rename columns 
names(perm_2_repeated_counts_set3_5percent)[names(perm_2_repeated_counts_set3_5percent) == "sapply(rowstat_set3_5percent[names(rowstat_set3_5percent)], function(x) sum(x > 1))"] <- "outlier_counts"
names(perm_3_repeated_counts_set3_5percent)[names(perm_3_repeated_counts_set3_5percent) == "sapply(rowstat_set3_5percent[names(rowstat_set3_5percent)], function(x) sum(x > 2))"] <- "outlier_counts"


#empirical estimates for outliers shared between more than 1 or more than 2 population pairs
empirical_2_repeated_counts_set3_5percent<-sum(FST_set3_5percent$outlier_count > 1)
empirical_3_repeated_counts_set3_5percent<-sum(FST_set3_5percent$outlier_count > 2)

#calculate significance
"see https://genomicsclass.github.io/book/pages/permutation_tests.html"
pval_set3_2_5percent<-(sum(perm_2_repeated_counts_set3_5percent$outlier_counts > empirical_2_repeated_counts_set3_5percent) + 1) / (length(perm_2_repeated_counts_set3_5percent$outlier_counts) + 1)
pval_set3_3_5percent<-(sum(perm_3_repeated_counts_set3_5percent$outlier_counts > empirical_3_repeated_counts_set3_5percent) + 1) / (length(perm_3_repeated_counts_set3_5percent$outlier_counts) + 1)


#plot results (permutations and empirical)
ggplot(data = perm_2_repeated_counts_set3_5percent, aes(x=outlier_counts)) +
  geom_density(adjust = 3,colour="black", fill="grey97")+
  geom_vline(xintercept=empirical_2_repeated_counts_set3_5percent,color="red", size=0.5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlab("\\nNumber of shared FST outliers")

ggplot(data = perm_3_repeated_counts_set3_5percent, aes(x=outlier_counts)) +
  geom_density(adjust = 3,colour="black", fill="grey97")+
  geom_vline(xintercept=empirical_3_repeated_counts_set3_5percent,color="red", size=0.5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlab("\\nNumber of shared FST outliers")

#repeated outlier counts based on the data and expected counts based on permutations (outliers repeated in 2 or more populations, and 3 or more populations)
#these counts are used for barplots below
empirical_2_repeated_counts_set3_5percent
empirical_3_repeated_counts_set3_5percent

median(perm_2_repeated_counts_set3_5percent$outlier_counts)
median(perm_3_repeated_counts_set3_5percent$outlier_counts)

pval_set3_2_5percent
pval_set3_3_5percent


###barplots for observed and expected values (invasive sets 1-3 are averaged)
results <- read.csv("FST_outlier_permutations/results_native_invasive.csv", sep = ',')
results$Category = factor(results$Category, levels=c("Observed", "Expected"))
results_native<-results[results$Dataset == 'Native',]
results_invasive<-results[results$Dataset == 'Invasive_average',]


ggplot(data=results_native, aes(Category,X2_or_more_outliers))+
  geom_bar(stat="identity",aes(fill=Category), colour="black", width=0.8)+
  theme_bw()+
  theme(axis.title=element_text(size=14),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")+
  ylim(0,350)+
  labs(x = NULL, y = "Number of windows") +
  scale_fill_manual(values = c("#9498a7ff","#a7d0a8ff"))

ggplot(data=results_native, aes(Category,X3_or_more_outliers))+
  geom_bar(stat="identity",aes(fill=Category), colour="black", width=0.8)+
  theme_bw()+
  theme(axis.title=element_text(size=14),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")+
  ylim(0,50)+
  labs(x = NULL, y = "Number of windows") +
  scale_fill_manual(values = c("#9498a7ff","#a7d0a8ff"))


ggplot(data=results_invasive, aes(Category,X2_or_more_outliers))+
  geom_bar(stat="identity",aes(fill=Category), colour="black", width=0.8)+
  geom_errorbar(aes(ymin=X2_or_more_outliers, ymax=X2_or_more_outliers+X2_or_more_outliers_SE,colour=Category), width=0, size=1.5)+
  theme_bw()+
  theme(axis.title=element_text(size=14),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")+
  ylim(0,350)+
  labs(x = NULL, y = "Number of windows") +
  scale_fill_manual(values = c("#9498a7ff","#a7d0a8ff"))+
  scale_colour_manual(values = c("#9498a7ff","#a7d0a8ff"))

ggplot(data=results_invasive, aes(Category,X3_or_more_outliers))+
  geom_bar(stat="identity",aes(fill=Category), colour="black", width=0.8)+
  geom_errorbar(aes(ymin=X3_or_more_outliers, ymax=X3_or_more_outliers+X3_or_more_outliers_SE,colour=Category), width=0, size=1.5)+
  theme_bw()+
  theme(axis.title=element_text(size=14),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")+
  ylim(0,50)+
  labs(x = NULL, y = "Number of windows") +
  scale_fill_manual(values = c("#9498a7ff","#a7d0a8ff"))+
  scale_colour_manual(values = c("#9498a7ff","#a7d0a8ff"))

#Bonferroni corrections for permutation test results, Native populations set and 3 Invasive population sets
pvalues<-c("9.999e-05","9.999e-05","9.999e-05","9.999e-05","9.999e-05","9.999e-05","9.999e-05","9.999e-05")
p.adjust(pvalues, method="bonferroni")
