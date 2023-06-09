#set and confirm working directory
getwd()
setwd("D:/oa_ma/open data files for submission")
#setwd("/Volumes/NO NAME/oa_ma")
getwd()

#load relevant packages for analysis and visualization
library(dplyr)
library(ggplot2)
library(metafor)

#import data from CSV file
data <- read.csv("element_analysis_dataset_complete.csv", header = TRUE); data

#calculate variances for each study. variance uses the equations normally found with lnRR in OA analyses
data$variance <- ((data$SDReduc^2)/((data$n.reduc*(data$AverageReduc^2))))+((data$SDAmb^2)/((data$n.amb*(data$AverageAmb^2))))


names(data) #check all data is present
levels(data$Order) #should only be decapods and sessilian barnacles
levels(data$Anatomy)


#filter data so only Ca + Mg, for decapods and sessilian barnacles, are present
data <- filter(data, data$Order == "Sessilia" | data$Order == "Decapoda")
data <- filter(data, data$Element == "Calcium" | data$Element == "Magnesium")
unique(data$ReferenceShort)


#makes sure the relevant things are numeric variables
data$lnRR <- as.numeric(as.character(data$lnRR))
data$pco2.reduc <- as.numeric(as.character(data$pco2.reduc))

#subgroup by acidified pCO2 level
co21 <- filter(data, pco2.reduc >= 0 & pco2.reduc < 500)
c022 <- filter(data, pco2.reduc >= 500 & pco2.reduc < 1000)
co23 <- filter(data, pco2.reduc >= 1000 & pco2.reduc < 1500)
c024 <- filter(data, pco2.reduc >= 1500 & pco2.reduc < 2000)
co25 <- filter(data, pco2.reduc > 2000)

#use the rma function from metafor to run meta analyses on each CO2 bin. variance is the within-study variance metric calculated above. "DL" is the DerSimion-Laird method to approximate the between-study variance
rma(lnRR, variance, method = "DL", data = co21)
rma(lnRR, variance, method = "DL", data = c022)
rma(lnRR, variance, method = "DL", data = co23)
rma(lnRR, variance, method = "DL", data = c024)
rma(lnRR, variance, method = "DL", data = co25)

#note: REML and DL moethods produce similar outputs, either can be used I'd say
#save RMA results in a data frame to plot
ion.pco2.bins <- data.frame(name = c("500-999", "1000-1499", "1500-1999", "2000+"),
                                  mean = c(0.0002, 0.0173,  -0.0793, 0.0135),
                                  lci = c(-0.0263, -0.0639, -0.1196, -0.0230),
                                  uci = c(0.0267, 0.0986, -0.0390, 0.0500),
                                  counts = c(30, 21,  29, 16),
                                  stat = c("p = .9865",  "p = .6761", "p = .0001", "p = .4694"))

#rename and reorder bin levels for visualization
ion.pco2.bins$name <- factor(ion.pco2.bins$name,
                                   levels=c("500-999", "1000-1499", "1500-1999", "2000+"))

#now, filter so we have calcium only pCO2 bins
calcium <- filter(data, data$Element == "Calcium")
calcium$adj.lnrr <- calcium$lnRR/(calcium$variance+0.0042) #this effect size accounts for the between study and within study variance when plotting for the study duration analyses

max(abs(calcium$lnRR))
min(calcium$lnRR)
rma(lnRR, variance, data = calcium, method = "DL")

#filter calcium data points into subgroups based on high CO2 levels
calc.500 <- filter(calcium, pco2.reduc >= 0 & pco2.reduc < 500)
calc.1000 <- filter(calcium, pco2.reduc >= 500 & pco2.reduc < 1000)
calc.1500 <- filter(calcium, pco2.reduc >= 1000 & pco2.reduc < 1500)
calc.2000 <- filter(calcium, pco2.reduc >= 1500 & pco2.reduc < 2000)
calc.high <- filter(calcium, pco2.reduc > 2000)

max(abs(calc.2000$lnRR))
max(abs(calc.high$lnRR))

#calculate summary effect sizes for each CO2 bin based on random effects model using DL approximation of between study variance
rma(lnRR, variance, method = "DL", data = calc.500)
rma(lnRR, variance, method = "DL", data = calc.1000)
rma(lnRR, variance, method = "DL", data = calc.1500)
rma(lnRR, variance, method = "DL", data = calc.2000)
rma(lnRR, variance, method = "DL", data = calc.high)

#sensitivity test: remove Page et al. (2017) and determine the result
calc.pageless <- filter(calc.2000, calc.2000$ï..ReferenceShort != "Page et al 2017")
rma(lnRR, variance, method = "DL", data = calc.pageless)

#back to main analysis: save Ca summary effect sizes and details for each pCO2 bin in a data frame
calcium.pco2.bins <- data.frame(name = c("500-999", "1000-1499", "1500-1999", "2000+"),
                            mean = c(0.0065, 0.0193,  -0.0796, -0.0170),
                            lci = c(-0.0193, -0.1141, -0.1324, -0.0483),
                            uci = c(0.0323, 0.1528, -0.0268, 0.0143),
                            counts = c(16, 12,  15, 8),
                            stat = c("p = .6210",  "p = .7663", "p = .0031", "p = .2879"))

#reorder pCO2 bin names for plotting
calcium.pco2.bins$name <- factor(ion.pco2.bins$name,
                             levels=c("500-999", "1000-1499", "1500-1999", "2000+"))

#subgroup analyses by biological predictors, using mixed effect models (factors as moderators, "DL" to approximate between study variance)
#subgroup analysis: by order
rma(lnRR, variance, data = calcium, method = 'DL', mods = ~factor(Order)-1) #subgroup by order; stat sig decrease for Decapods
rma(lnRR, variance, data = calc.500, method = 'DL', mods = ~factor(Order)-1)
rma(lnRR, variance, data = calc.1000, method = 'DL', mods = ~factor(Order)-1) 
rma(lnRR, variance, data = calc.1500, method = 'DL', mods = ~factor(Order)-1)
rma(lnRR, variance, data = calc.2000, method = 'DL', mods = ~factor(Order)-1)
rma(lnRR, variance, data = calc.high, method = 'DL', mods = ~factor(Order)-1)

#subgroup analysis: by biogeography
rma(lnRR, variance, data = calcium, method = 'DL', mods = ~factor(Region)-1) #subgroup by biogeographic region
rma(lnRR, variance, data = calc.500, method = 'DL', mods = ~factor(Region)-1)
rma(lnRR, variance, data = calc.1000, method = 'DL', mods = ~factor(Region)-1)
rma(lnRR, variance, data = calc.1500, method = 'DL', mods = ~factor(Region)-1)
rma(lnRR, variance, data = calc.2000, method = 'DL', mods = ~factor(Region)-1)
rma(lnRR, variance, data = calc.high, method = 'DL', mods = ~factor(Region)-1)

#subgroup analysis: by developmental stage
rma(lnRR, variance, data = calcium, method = 'DL', mods = ~factor(DevStage)-1) #subgroup by developmental stage
rma(lnRR, variance, data = calc.500, method = 'DL', mods = ~factor(DevStage)-1)
rma(lnRR, variance, data = calc.1000, method = 'DL', mods = ~factor(DevStage)-1)
rma(lnRR, variance, data = calc.1500, method = 'DL', mods = ~factor(DevStage)-1)
rma(lnRR, variance, data = calc.2000, method = 'DL', mods = ~factor(DevStage)-1)
rma(lnRR, variance, data = calc.high, method = 'DL', mods = ~factor(DevStage)-1)

#subgroup analysis: by anatomical region
rma(lnRR, variance, data = calcium, method = 'DL', mods = ~factor(Anatomy)-1) #subgroup by anatomy
rma(lnRR, variance, data = calc.500, method = 'DL', mods = ~factor(Anatomy)-1)
rma(lnRR, variance, data = calc.1000, method = 'DL', mods = ~factor(Anatomy)-1)
rma(lnRR, variance, data = calc.1500, method = 'DL', mods = ~factor(Anatomy)-1)
rma(lnRR, variance, data = calc.2000, method = 'DL', mods = ~factor(Anatomy)-1)
rma(lnRR, variance, data = calc.high, method = 'DL', mods = ~factor(Anatomy)-1)

#subgroup analysis: by study duration
hist(calcium$Length)
rma(lnRR, variance, data = calcium, method = 'DL', mods = ~cbind(lnLength)-1) #subgroup by study duration
rma(lnRR, variance, data = calc.500, method = 'DL', mods = ~cbind(lnLength)-1)
rma(lnRR, variance, data = calc.1000, method = 'DL', mods = ~cbind(lnLength)-1)
rma(lnRR, variance, data = calc.1500, method = 'DL', mods = ~cbind(lnLength)-1)
rma(lnRR, variance, data = calc.2000, method = 'DL', mods = ~cbind(lnLength)-1)
rma(lnRR, variance, data = calc.high, method = 'DL', mods = ~cbind(lnLength)-1)

#magnesium only data points, similarly subgrouped by high pCO2 bins
magnesium <- filter(data, data$Element == "Magnesium")
rma(lnRR, variance, data = magnesium, method = "DL")

#Mg data points subsetted by high CO2 levels
magn.500 <- filter(magnesium, pco2.reduc >= 0 & pco2.reduc < 500)
magn.1000 <- filter(magnesium, pco2.reduc >= 500 & pco2.reduc < 1000)
magn.1500 <- filter(magnesium, pco2.reduc >= 1000 & pco2.reduc < 1500)
magn.2000 <- filter(magnesium, pco2.reduc >= 1500 & pco2.reduc < 2000)
magn.high <- filter(magnesium, pco2.reduc > 2000)

#random effects models run on CO2 bins for summary effect sizes
rma(lnRR, variance, method = "DL", data = magn.500)
rma(lnRR, variance, method = "DL", data = magn.1000)
rma(lnRR, variance, method = "DL", data = magn.1500)
rma(lnRR, variance, method = "DL", data = magn.2000)
rma(lnRR, variance, method = "DL", data = magn.high)

#similar to Ca dataset, removal of Page et al. (2017) values and reanalysis of effect size
magn.2000.pageless <- filter(magn.2000, magn.2000$ï..ReferenceShort != "Page et al 2017")
rma(lnRR, variance, method = "DL", data = magn.2000.pageless)

#save summary effect sizes and details into a data frame
magnesium.pco2.bins <- data.frame(name = c("500-999", "1000-1499", "1500-1999", "2000+"),
                                mean = c(-0.0196, 0.0152,  -0.0799, 0.0483),
                                lci = c(-0.0844, -0.0678, -0.1460, -0.0212),
                                uci = c(0.0451, 0.0981, -0.0138, 0.1177),
                                counts = c(14, 9,  14, 8),
                                stat = c("p = .5521",  "p = .7204", "p = .0178", "p = .1731"))

#rename and reorder bin names for plotting
magnesium.pco2.bins$name <- factor(ion.pco2.bins$name,
                                 levels=c("500-999", "1000-1499", "1500-1999", "2000+"))

#subgroup analyses (biological moderators) for each CO2 bin
#subgroup analysis: by order
rma(lnRR, variance, data = magnesium, method = 'DL', mods = ~factor(Order)-1) #subgroup by order
rma(lnRR, variance, data = magn.500, method = 'DL', mods = ~factor(Order)-1)
rma(lnRR, variance, data = magn.1000, method = 'DL', mods = ~factor(Order)-1)
rma(lnRR, variance, data = magn.1500, method = 'DL', mods = ~factor(Order)-1)
rma(lnRR, variance, data = magn.2000, method = 'DL', mods = ~factor(Order)-1)
rma(lnRR, variance, data = magn.high, method = 'DL', mods = ~factor(Order)-1)

#subgroup analysis: by biogeography
rma(lnRR, variance, data = magnesium, method = 'DL', mods = ~factor(Region)-1) #subgroup by biogeographic region
rma(lnRR, variance, data = magn.500, method = 'DL', mods = ~factor(Region)-1)
rma(lnRR, variance, data = magn.1000, method = 'DL', mods = ~factor(Region)-1)
rma(lnRR, variance, data = magn.1500, method = 'DL', mods = ~factor(Region)-1)
rma(lnRR, variance, data = magn.2000, method = 'DL', mods = ~factor(Region)-1)
rma(lnRR, variance, data = magn.high, method = 'DL', mods = ~factor(Region)-1)

#subgroup analysis: by developmental stage
rma(lnRR, variance, data = magnesium, method = 'DL', mods = ~factor(DevStage)-1) #subgroup by developmental stage
rma(lnRR, variance, data = magn.500, method = 'DL', mods = ~factor(DevStage)-1)
rma(lnRR, variance, data = magn.1000, method = 'DL', mods = ~factor(DevStage)-1)
rma(lnRR, variance, data = magn.1500, method = 'DL', mods = ~factor(DevStage)-1)
rma(lnRR, variance, data = magn.2000, method = 'DL', mods = ~factor(DevStage)-1)
rma(lnRR, variance, data = magn.high, method = 'DL', mods = ~factor(DevStage)-1)

#subgroup analysis: by anatomical region
rma(lnRR, variance, data = magnesium, method = 'DL', mods = ~factor(Anatomy)-1) #subgroup by anatomy
rma(lnRR, variance, data = magn.500, method = 'DL', mods = ~factor(Anatomy)-1)
rma(lnRR, variance, data = magn.1000, method = 'DL', mods = ~factor(Anatomy)-1)
rma(lnRR, variance, data = magn.1500, method = 'DL', mods = ~factor(Anatomy)-1)
rma(lnRR, variance, data = magn.2000, method = 'DL', mods = ~factor(Anatomy)-1)
rma(lnRR, variance, data = magn.high, method = 'DL', mods = ~factor(Anatomy)-1)

#finally, combined plot as displayed in Figure 2
#try ggplot faceting instead
all_data <- rbind(ion.pco2.bins, calcium.pco2.bins, magnesium.pco2.bins)
all_data

all_data_plot <- ggplot(all_data, aes(x = name, y = mean, group =1))+
  geom_point(shape=21, size=3, fill = "white")+
  geom_errorbar(width=.1, aes(ymin = lci, ymax = uci))+
  geom_text(aes(x = name, y = uci, label = counts, vjust = -0.5))+
  geom_text(aes(x = name, y = uci, label = stat, vjust = -5.0))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  ylim(-0.3, 0.3)+
  theme_classic()+
  xlab("pCO2 level for acidified treatment (\\u03BCatm)")+
  ylab("Mean effect size (lnRR)")+
  theme(axis.title=element_text(size = 12), axis.text = element_text(size = 12))+
  #facet_grid(identity~.)
  facet_wrap(~identity, ncol = 1)
all_data_plot
all_data_plot_labels <- data.frame(label = c("A", "B", "C"),
                                   identity = c("Combined", "Calcium", "Magnesium"),
                                   x = c(0.5, 0.5, 0.5),
                                   y = c(0.25, 0.25, 0.25))
all_data_plot + geom_text(data = all_data_plot_labels,
                          mapping = aes(x = x, y = y, label = label), size = 6)

all_data_plot + annotate("text", x = 0., y = 0.25, label = "A", size = 5)