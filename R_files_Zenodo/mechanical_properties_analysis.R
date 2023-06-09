#set and confirm working directory
getwd()
setwd("D:/oa_ma/open data files for submission")
setwd("/Volumes/NO NAME/oa_ma/Material analysis")
getwd()

#load relevant packages for analysis and visualization
library(dplyr)
library(ggplot2)
library(metafor)

#import data from CSV file and check all lines were added to variable
data <- read.csv("mechanical_properties_data.csv", header = TRUE); data$ï..ReferenceShort

#calculate variances for each study. s.pool uses the equation for sigma (usually an intermediate in other standardized mean difference effect sizes). variance uses the equations normally found with lnRR in OA analyses
data$variance <- ((data$SDReduc^2)/((data$nReduc*(data$AverageReduc^2))))+((data$SDAmb^2)/((data$nAmb*(data$AverageAmb^2))))

#filter out empty datasets
data <- filter(data, data$ReferenceShort != "Chadwick et al. 2019")
data <- filter(data, data$ReferenceShort != "Landes and Zimmer 2012")
unique(data$Order)

#make sure relevant parameters are recognized as numeric variables by R
data$lnRR <- as.numeric(as.character(data$lnRR))
data$delta.ph <- as.numeric(as.character(data$delta.ph))
data$lnLength <- as.numeric(as.character(data$lnLength))
data$pco2.reduc <- as.numeric(as.character(data$pco2.reduc))

data$pco2.bin <- as.factor(data$pco2.bin)

#separate into metrics, hardness and thickness
hardness <- filter(data, data$Metric == "Hardness"); hardness$lnRR
thickness <- filter(data, data$Metric == "Thickness")

#separate hardness into pco2 ranges
hist(hardness$pco2.reduc, main = "High pCO2 level for hardness dataset", breaks = 5) # looks like n = 12 for 500-100, n = 6 for 1500-2000, and n = 1 for 2500-3000
hardness.500 <- filter(hardness, pco2.reduc < 500) # n = 0
hardness.1000 <- filter(hardness, pco2.reduc >= 0 & pco2.reduc < 1000) #n = 12
hardness.1500 <- filter(hardness, pco2.reduc >= 1000 & pco2.reduc < 1500) #n = 0
hardness.2000 <- filter(hardness, pco2.reduc >= 1500 & pco2.reduc < 2000) # n = 6
hardness.high.pco2 <- filter(hardness, pco2.reduc >= 2000) # n = 1

#similarly, separate thickness datatset by high pCO2 levels
hist(thickness$pco2.reduc, main = "High pCO2 level for thickness dataset", breaks =10)
thickness.500 <- filter(thickness, pco2.reduc < 500)
thickness.1000 <- filter(thickness, pco2.reduc >= 500 & pco2.reduc < 1000) # n = 1000
thickness.1500 <- filter(thickness, pco2.reduc >= 1000 & pco2.reduc < 1500) # n = 1
thickness.2000 <- filter(thickness, pco2.reduc >= 1500 & pco2.reduc < 2000) # n = 6
thickness.high.pco2 <- filter(thickness, pco2.reduc >= 2000) # n = 1

hardness.overall <- rma(lnRR, variance, method = "DL", data = hardness); hardness.overall
rma(lnRR, variance, method = "DL", data = hardness.500) #no data
rma(lnRR, variance, method = "DL", data = hardness.1000) #negative and statistically significant
rma(lnRR, variance, method = "DL", data = hardness.1500) #no data
rma(lnRR, variance, method = "DL", data = hardness.2000) #negative and statistically significant
rma(lnRR, variance, method = "DL", data = hardness.high.pco2) #sample size of 1

hardness.pco2.bins <- data.frame(name = c("500-999", "1000-1499", "1500-1999", "2000+"),
                               mean = c(-0.0796, 0,  -0.2450, -0.0496),
                               lci = c(-0.1548, 0, -0.4255, -0.2482),
                               uci = c(-0.0044, 0, -0.0645, 0.1490),
                               counts = c(12, 0,  6, 1),
                               stat = c("p = .0380", "No data",  "p = .0078", "p = .6246"))
hardness.pco2.bins$name <- factor(hardness.pco2.bins$name,
                                  levels=c("500-999", "1000-1499", "1500-1999","2000+"))

#thickness
thickness.overall <- rma(lnRR, variance, method = "DL", data = thickness); thickness.overall
rma(lnRR, variance, method = "DL", mods = ~factor(Order)-1, data = thickness)
rma(lnRR, variance, method = "DL", mods = cbind(delta.ph), data = thickness)
rma(lnRR, variance, method = "DL", mods = cbind(lnLength), data = thickness)

rma(lnRR, variance, method = "DL", data = thickness.500)
rma(lnRR, variance, method = "DL", data = thickness.1000)
rma(lnRR, variance, method = "DL", data = thickness.1500)
rma(lnRR, variance, method = "DL", data = thickness.2000)
rma(lnRR, variance, method = "DL", data = thickness.high.pco2)


thickness.pco2.bins <- data.frame(name = c("500-999", "1000-1499", "1500-1999", "2000+"),
                                 mean = c(-0.0476, -0.1645,  -0.0891, -0.0645),
                                 lci = c(-0.1039, -0.6950, -0.1926, -0.1541),
                                 uci = c(0.0087, 0.3661, 0.0144, 0.0250),
                                 counts = c(10, 1,  6, 1),
                                 stat = c("p = .0974",  "p = .5435", "p = .0917", "p = .1578"))

thickness.pco2.bins$name <- factor(thickness.pco2.bins$name,
                                   levels=c("500-999", "1000-1499", "1500-1999", "2000+"))

#plot both together
all_data_mat <- rbind(hardness.pco2.bins, thickness.pco2.bins)
all_data_mat

all_data_plot_mat <- ggplot(all_data_mat, aes(x = name, y = mean, group =1))+
  geom_point(shape=21, size=3, fill = "white")+
  geom_errorbar(width=.1, aes(ymin = lci, ymax = uci))+
  geom_text(aes(x = name, y = uci, label = counts, vjust = -0.5))+
  geom_text(aes(x = name, y = uci, label = stat, vjust = -5.0))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  ylim(-0.7, 0.7)+
  theme_classic()+
  xlab("pCO2 level for acidified treatment (\\u03BCatm)")+
  ylab("Mean effect size (lnRR)")+
  theme(axis.title=element_text(size = 12), axis.text = element_text(size = 12))+
  #facet_grid(identity~.)
  facet_wrap(~identity, ncol = 1)

all_data_plot_mat
all_data_plot_mat_labels <- data.frame(label = c("A", "B"),
                                       identity = c("Biomechanics", "Total cuticle thickness"),
                                       x = c(0.5, 0.5),
                                       y = c(0.6, 0.6))
all_data_plot_mat + geom_text(data = all_data_plot_mat_labels,
                              mapping = aes(x = x, y = y, label = label), size = 6)