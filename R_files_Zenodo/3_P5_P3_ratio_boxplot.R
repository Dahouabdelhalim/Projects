library(ggplot2)
library(plyr)

# Load data
pedipalp_measurements <- read.csv(paste(getwd(), "/all_measurements_and_ratios.csv", sep=""))

# Convert Species to a factor variable
pedipalp_measurements$Species=as.factor(pedipalp_measurements$Species)
head(pedipalp_measurements)

# Convert Sex to a factor variable
pedipalp_measurements$Sex=as.factor(pedipalp_measurements$Sex)
head(pedipalp_measurements)

# Subset by Sex
test_subset_group_males <- pedipalp_measurements[pedipalp_measurements$Sex== "male",]
test_subset_group_females <- pedipalp_measurements[pedipalp_measurements$Sex== "female",]

# Plot data

pdf(file = paste(getwd(), '/P5_P3_ratio_boxplot.pdf', sep = ''))
ggplot(pedipalp_measurements, aes(x=Species, y=P5_P3_ratio)) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', stackratio=1, dotsize=0.5) +
  stat_summary(fun=mean, geom="point", shape=3, size=5, color="red")
dev.off()

pdf(file = paste(getwd(), '/P5_P3_ratio_males_boxplot.pdf', sep = ''))
ggplot(test_subset_group_males, aes(x=Species, y=P5_P3_ratio)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, size=2.5) + 
  scale_y_continuous(name="P5_P3_ratio",breaks=seq(0.7, 1.2, 0.1), limits=c(0.7,1.2)) +
  stat_summary(fun=mean, geom="point", shape=3, size=5, color="red") +
  theme(aspect.ratio=1)
dev.off()

pdf(file = paste(getwd(), '/P5_P3_ratio_females_boxplot.pdf', sep = ''))
ggplot(test_subset_group_females, aes(x=Species, y=P5_P3_ratio)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, size=2.5) + 
  scale_y_continuous(name="P5_P3_ratio",breaks=seq(0.7, 1.2, 0.1), limits=c(0.7,1.2)) +
  stat_summary(fun=mean, geom="point", shape=3, size=5, color="red") +
  theme(aspect.ratio=1)
dev.off()

pdf(file = paste(getwd(), '/P5_P3_ratio_boxplots_both_Sexes.pdf', sep = ''))
ggplot(pedipalp_measurements, aes(x=Species, y=P5_P3_ratio, fill=Sex)) + 
  geom_boxplot()+
  scale_y_continuous(name="P5_P3_ratio",breaks=seq(0.7, 1.2, 0.1), limits=c(0.7,1.2)) +
  stat_summary(fun=mean, geom="point", shape=3, size=5, color="red") +
  theme(aspect.ratio=1)
dev.off()
