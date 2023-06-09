################# Get this party started #################
###Set Working Directory###
setwd("C:/Users/sasha/OneDrive/Graduate Studies/Thesis/Think Before They Squeak/Submissions/Raw Data/")

###Load packages###
library(ape)
library(caper) #PGLS Modeling
library(dplyr) #Y
library(geiger) #Phylogenetic tree and pgls
library(ggmap) #Mapping with ggplot
library(ggplot2) #Plotting
library(ggpubr) #ggplot extension for combining plots
library(glmm)
library(phytools)
library(plyr)
library(reshape2) #Melting, reshaping, and merging dataframes
library(tidyr) #Y
library(tidyverse)
library(viridis) #Beautiful colour schemes for ggplot - colour blind approved?

################# Load and Merge Datasets #################
#Call Data
Squ_Calls <- read.csv("Squirrel_Calls.csv", header=T, na.strings=c("","NA"))
Squ_Calls$ScName <- Squ_Calls$Scientific_Name
Squ_Calls$Scientific_Name <- as.factor(gsub(" ", "_", Squ_Calls$Scientific_Name))
str(Squ_Calls)

#Ecological Data
Squ_Eco <- read.csv("Squirrel_Ecological_Traits.csv", header = T, na.strings = c("","NA"))
Squ_Eco$Scientific_Name <- gsub(" ", "_", Squ_Eco$Scientific_Name)
str(Squ_Eco)

#Merge Datasets
Squirrels <- merge(Squ_Calls, Squ_Eco, by = "Scientific_Name", all.x = T)
Squirrels <- Squirrels[, -which(names(Squirrels) %in% c("Common_Name.y"))]

#Cleaning up categorical variables
Squirrels$Openness <- factor(Squirrels$H_Openness)
Squirrels$Sociality <- factor(Squirrels$Sociality)
Squirrels$Soc <- as.factor(ifelse(Squirrels$Sociality == "CFTM", "Social", 
                                  ifelse(Squirrels$Sociality == "CFWM", "Social", 
                                         ifelse(Squirrels$Sociality == "Colonial", "Social", 
                                                ifelse(Squirrels$Sociality == "HFCM", "Social", 
                                                       ifelse(Squirrels$Sociality == "HFTM", "Social", 
                                                              ifelse(Squirrels$Sociality == "Alternate", "Social",
                                                                     ifelse(Squirrels$Sociality == "Monogamous", "Social", 
                                                                            "Solitary"))))))))

Squirrels <- Squirrels[!(Squirrels$Sex == "Male"),]
Squirrels <- Squirrels[!(Squirrels$Age == "Pup"),]

################# Frequencies #################

#Fundamental Frequency
Freq_fun_temp <- Squirrels %>%
  group_by(Scientific_Name) %>%
  slice(which.max(Fun_Freq_kHz))
Freq_fun <- as.data.frame(Freq_fun_temp[,c("Scientific_Name", "ScName", "Mass_g_Females", "Time_Partitioning", "Openness", "Soc", "Fun_Freq_kHz", "Year", "Ana_Hz_Max")])
Freq_fun$Ana_Max <- log(Freq_fun$Ana_Hz_Max + 1)
Freq_fun$Time_Partitioning[Freq_fun$Time_Partitioning == "Crepuscular"] <- NA
Freq_fun <- Freq_fun %>%
  drop_na(Time_Partitioning, Soc, Ana_Max)
row.names(Freq_fun) <- Freq_fun$Scientific_Name
Freq_fun$Time_Partitioning <- factor(Freq_fun$Time_Partitioning)

###Dominant Frequency
Freq_dom_temp <- Squirrels %>%
  group_by(Scientific_Name) %>%
  slice(which.max(Dom_Freq_kHz))
Freq_dom <- as.data.frame(Freq_dom_temp[,c("Scientific_Name", "ScName", "Mass_g_Females", "Time_Partitioning", "Openness", "Soc", "Dom_Freq_kHz", "Year", "Ana_Hz_Max")])
Freq_dom$Ana_Max <- log(Freq_dom$Ana_Hz_Max + 1)
Freq_dom$Time_Partitioning[Freq_dom$Time_Partitioning == "Crepuscular"] <- NA
Freq_dom <- Freq_dom %>%
  drop_na(Time_Partitioning, Soc, Ana_Max)
row.names(Freq_dom) <- Freq_dom$Scientific_Name
Freq_dom$Time_Partitioning <- factor(Freq_dom$Time_Partitioning)

#Minimum Dominant Frequency
Freq_min <- Squirrels %>%
  group_by(Scientific_Name) %>%
  slice(which.min(Min_Freq_kHz))
Freq_min <- as.data.frame(Freq_min[,c("Scientific_Name", "ScName", "Mass_g_Females", "Time_Partitioning", "Openness", "Soc", "Min_Freq_kHz", "Year", "Ana_Hz_Max")])
Freq_min$Ana_Max <- log(Freq_min$Ana_Hz_Max + 1)
Freq_min$Time_Partitioning[Freq_min$Time_Partitioning == "Crepuscular"] <- NA
Freq_min <- Freq_min %>%
  drop_na(Time_Partitioning, Soc, Ana_Max)
row.names(Freq_min) <- Freq_min$Scientific_Name
Freq_min$Time_Partitioning <- factor(Freq_min$Time_Partitioning)

#Maximum Frequency
Freq_max <- Squirrels %>%
  group_by(Scientific_Name) %>%
  top_n(1, Max_Freq_kHz)
Freq_max <- as.data.frame(Freq_max[,c("Scientific_Name", "ScName", "Mass_g_Females", "Time_Partitioning", "Openness", "Soc", "Max_Freq_kHz", "Year", "Ana_Hz_Max")])
Freq_max$Ana_Max <- log(Freq_max$Ana_Hz_Max + 1)
Freq_max$Time_Partitioning[Freq_max$Time_Partitioning == "Crepuscular"] <- NA
Freq_max <- Freq_max %>%
  drop_na(Time_Partitioning, Soc, Ana_Max)
row.names(Freq_max) <- Freq_max$Scientific_Name
Freq_max$Time_Partitioning <- factor(Freq_max$Time_Partitioning)

#Maximum Visible Harmonic
Harm_max <- Squirrels %>%
  group_by(Scientific_Name) %>%
  top_n(1, Highest_Harmonic_kHz)
Harm_max <- Harm_max[!duplicated(Harm_max$Scientific_Name),]
Harm_max <- as.data.frame(Harm_max[,c("Scientific_Name", "ScName", "Mass_g_Females", "Time_Partitioning", "Openness", "Soc", "Highest_Harmonic_kHz", "Year", "Ana_Hz_Max")])
Harm_max$Ana_Max <- log(Harm_max$Ana_Hz_Max + 1)
Harm_max$Time_Partitioning[Harm_max$Time_Partitioning == "Crepuscular"] <- NA
Harm_max <- Harm_max %>%
  drop_na(Time_Partitioning, Soc, Ana_Max)
row.names(Harm_max) <- Harm_max$Scientific_Name
Harm_max$Time_Partitioning <- factor(Harm_max$Time_Partitioning)

###Combining Frequencies
Frequencies <- merge(Freq_fun, Freq_dom, by = "ScName")
Frequencies <- merge(Frequencies, Freq_max, by = "ScName")
Frequencies <- merge(Frequencies, Freq_min, by = "ScName")
Frequencies <- merge(Frequencies, Harm_max, by = "ScName")
colnames(Frequencies) <- c("ScName", "Scientific.Name.ff", "Body.Size...Females..g..ff", "Temporality.ff", "Openness.ff", "Soc.ff", 
                           "Fundamental.Freq..kHz.", "Year.ff", "Ana_Max.Hz.ff", "Scientific.Name.df", "Body.Size...Females..g..df",
                           "Temporality.df", "Openness.df", "Soc.df", "Dom.Frequency..mean", "Year.df", "Ana_Max.Hz.df", 
                           "Scientific.Name.mxf", "Body.Size...Females..g..mxf", "Temporality.mxf", "Openness.mxf", "Soc.mxf", 
                           "Max.Freq..kHz.", "Year.mxf", "Ana_Max.Hz.mxf", "Scientific.Name.mnf", "Body.Size...Females..g..mnf", 
                           "Temporality.mnf", "Openness.mnf", "Soc.mnf", "Min.Freq..kHz.", "Year.mnf", "Ana_Max.Hz.mnf", 
                           "Scientific.Name.hmf", "Body.Size...Females..g..hmf", "Temporality.hmf", "Openness.hmf", "Soc.hmf", 
                           "Mean_of_Highest_Harmonic_Visible", "Year.hmf", "Ana_Max.Hz.hmf")

################# Phylogenetic Tree #################
#Load Tree
VertLife_Trees <- read.nexus("output_200402.nex")
Squirrel_Tree <- VertLife_Trees$tree_2033 
#sample(1:100, 1) #Randomly selecting tree #321: tree_2033

#Add tips to tree
Squirrels_sci_name <- unique(Squirrels$Scientific.Name)
name.check(Squirrel_Tree, Squirrels_sci_name)
Squirrel_Tree <- bind.tip(Squirrel_Tree, "Glaucomys_sabrinus_coloratus", edge.length = 1, where = 65, position = 1)
Squirrel_Tree <- bind.tip(Squirrel_Tree, "Xerospermophilus_spilosoma_annectens", edge.length = 1, where = 40, position = 1)
Squirrel_Tree <- bind.tip(Squirrel_Tree, "Xerospermophilus_spilosoma_marginatus", edge.length = 1, where = 40, position = 1)
Squirrel_Tree <- bind.tip(Squirrel_Tree, "Spermophilus_musicus", edge.length = 1, where = 28, position = 1)


#Rename any tips
Squirrel_Tree$tip.label[Squirrel_Tree$tip.label == "Sciurus_aberti"] <- "Sciurus_aberti_kaibabensis"
Squirrel_Tree$tip.label[Squirrel_Tree$tip.label == "Sciurus_niger"] <- "Sciurus_niger_rufiventer"
Squirrel_Tree$tip.label[Squirrel_Tree$tip.label == "Callosciurus_erythraeus"] <- "Callosciurus_erythraeus_thaiwanensis"

#Drop tips
Squirrel_Tree <- drop.tip(Squirrel_Tree, 41)


par(mar = c(0,0,0,0))
plot(Squirrel_Tree, type = "p", no.margin = TRUE, cex = 0.5, direction = "right", 
     use.edge.length = TRUE)
#nodelabels()
add.scale.bar(cex = 0.7, font = 2, col = "red")

library(png)
library(grid)

Ptero_image <- readPNG("Glaucomys_volans.png")
Sciu_image <- readPNG("Tamiasciurus_hudsonicus.png")
Call_image <- readPNG("Callosciurus.png")
Marm_image <- readPNG("Marmota_marmota.png")
Sper_image <- readPNG("Urocitellus_richardsonii.png")
Tam_image <- readPNG("Tamias_striatus.png")
Xer_image <- readPNG("Xerus_inauris.png")

Ptero <- Squirrel_Tree$tip.label[c(66:70)] 
Sciu <- Squirrel_Tree$tip.label[c(61:65)]
Callo <- Squirrel_Tree$tip.label[c(57:60)]
Marm <- Squirrel_Tree$tip.label[c(43:52)]
Sper <- Squirrel_Tree$tip.label[c(20:42,53:58)]
Tam <- Squirrel_Tree$tip.label[c(2:20)]
Xer <- Squirrel_Tree$tip.label[c(1)]

#Long form for publication
#dev.new(width = 900, height = 1200, unit = "px")
tiff("Fig2.tiff", units="in", width=6,height=8, res=300)
plot(Squirrel_Tree,tip.color = ifelse(Squirrel_Tree$tip.label %in% Ptero, "#D55E00", 
                                      ifelse(Squirrel_Tree$tip.label %in% Sciu, "#E69F00",
                                             ifelse(Squirrel_Tree$tip.label %in% Callo, "#999999",
                                                    ifelse(Squirrel_Tree$tip.label %in% Marm,"#56B4E9",
                                                           ifelse(Squirrel_Tree$tip.label %in% Sper, "#0072B2",
                                                                  ifelse(Squirrel_Tree$tip.label %in% Tam, "#009E73", "#CC79A7")))))),
     type = "phylogram", no.margin = TRUE, cex = 0.65, direction = "right") 

add.scale.bar(x = 55, y = 1.3, cex = 0.65)

add.arrow(Squirrel_Tree, "Glaucomys_sabrinus_coloratus", col = "#D55E00", lwd = 1, offset = 1, arrl = 5, hedl = 2)
add.arrow(Squirrel_Tree, "Tamiasciurus_hudsonicus", col = "#E69F00", lwd = 1, offset = 2.5, arrl = 12, hedl = 2)
add.arrow(Squirrel_Tree, "Callosciurus_nigrovittatus", col = "#999999", lwd = 1, offset = 2.5, arrl = 5, hedl = 2)
add.arrow(Squirrel_Tree, "Marmota_vancouverensis", col = "#56B4E9", lwd = 1, offset = 2.5, arrl = 8, hedl = 2)
add.arrow(Squirrel_Tree, "Urocitellus_richardsonii", col = "#0072B2", lwd = 1, offset = 2.5, arrl = 8, hedl = 2)
add.arrow(Squirrel_Tree, "Tamias_sonomae", col = "#009E73", lwd = 1, offset = 2.5, arrl = 8, hedl = 2)
add.arrow(Squirrel_Tree, "Xerus_inauris", col = "#CC79A7", lwd = 1, offset = 3, arrl = 8, hedl = 2)

grid.raster(Ptero_image, 0.85, 0.97, width = 0.15)
grid.raster(Call_image, 0.83, 0.78, width = 0.15)
grid.raster(Marm_image, 0.84, 0.65, width = 0.17)
grid.raster(Sper_image, 0.85, 0.45, width = 0.08)
grid.raster(Tam_image, 0.82, 0.17, width = 0.1)
grid.raster(Xer_image, 0.75, 0.07, width = 0.1)
grid.raster(Sciu_image, 0.92, 0.87, width = 0.1)

dev.off()


################# Minimum Frequency Models #################
combined_LFQ <- comparative.data(
  Squirrel_Tree, Freq_min, names.col = "Scientific.Name", 
  vcv = TRUE, na.omit = TRUE)

#Null
model0_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ Ana_Min.Hz, data = combined_LFQ, lambda = "ML")
#summary(model0_LFQ)
AIC(model0_LFQ)

#Body mass
model1_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Ana_Min.Hz, data = combined_LFQ, lambda = "ML")
summary(model1_LFQ)
AIC(model1_LFQ)

#Temporality
model2_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ Temporality + Ana_Min.Hz, data = combined_LFQ, lambda = "ML")
#summary(model2_LFQ)
AIC(model2_LFQ)

#Sociality
model3_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ Soc + Ana_Min.Hz, data = combined_LFQ, lambda = "ML")
#summary(model3_LFQ)
AIC(model3_LFQ)

#Life History
model4_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ Openness + Ana_Min.Hz, data = combined_LFQ, lambda = "ML")
#summary(model4_LFQ)
AIC(model4_LFQ)

#Body mass, temporality
model5_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Temporality + Ana_Min.Hz,
                   data = combined_LFQ, lambda = "ML")
summary(model5_LFQ)
AIC(model5_LFQ)

#Body mass, sociality
model6_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Soc + Ana_Min.Hz, 
                   data = combined_LFQ, lambda = "ML")
summary(model6_LFQ)
AIC(model6_LFQ)

#Body mass, life.history
model7_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Openness + Ana_Min.Hz, 
                   data = combined_LFQ, lambda = "ML")
summary(model7_LFQ)
AIC(model7_LFQ)

#Body mass, temporality, sociality
model8_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Temporality +Soc + Ana_Min.Hz, 
                   data = combined_LFQ, lambda = "ML")
#summary(model8_LFQ)
AIC(model8_LFQ)

#Body mass, temporality, life history
model9_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Temporality + Openness + Ana_Min.Hz, 
                   data = combined_LFQ, lambda = "ML")
#summary(model9_LFQ)
AIC(model9_LFQ)

#Body mass, sociality, life history
model10_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Soc + Openness + Ana_Min.Hz, 
                    data = combined_LFQ, lambda = "ML")
#summary(model10_LFQ)
AIC(model10_LFQ)

#Body mass, temporality, sociality, life-history FULL MODEL
model11_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Temporality + Soc + Openness + Ana_Min.Hz, 
                    data = combined_LFQ, lambda = "ML")
summary(model11_LFQ)
AIC(model11_LFQ)
anova(model11_LFQ)

#Temporality, sociality
model12_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ Temporality + Soc + Ana_Min.Hz, data = combined_LFQ, lambda = "ML")
#summary(model12_LFQ)
AIC(model12_LFQ)

#Temporality, life history
model13_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ Temporality + Openness + Ana_Min.Hz, data = combined_LFQ, lambda = "ML")
#summary(model13_LFQ)
AIC(model13_LFQ)

#Sociality, life history
model14_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ Soc + Openness + Ana_Min.Hz, data = combined_LFQ, lambda = "ML")
#summary(model14_LFQ)
AIC(model14_LFQ)

#Temporality, sociality, life-history
model15_LFQ <- pgls(log(Min.Freq..kHz. + 1) ~ Temporality + Soc + Openness + Ana_Min.Hz, data = combined_LFQ, lambda = "ML")
#summary(model15_LFQ)
AIC(model15_LFQ)


################# Maximum Frequency Models #################
combined_HFQ <- comparative.data(
  Squirrel_Tree, Freq_max, names.col = "Scientific.Name", 
  vcv = TRUE, na.omit = FALSE)

#Null
model0_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ Ana_Max.Hz, 
                   data = combined_HFQ, lambda = "ML")
summary(model0_HFQ)
AIC(model0_HFQ)

#Body mass
model1_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Ana_Max.Hz, 
                   data = combined_HFQ, lambda = "ML")
summary(model1_HFQ)
AIC(model1_HFQ)

#Temporality
model2_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ Temporality + Ana_Max.Hz, data = combined_HFQ, 
                   lambda = "ML")
#summary(model2_HFQ)
AIC(model2_HFQ)

#Sociality
model3_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ Soc + Ana_Max.Hz, data = combined_HFQ, lambda = "ML")
#summary(model3_HFQ)
AIC(model3_HFQ)

#Life History
model4_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ Openness + Ana_Max.Hz, data = combined_HFQ, 
                   lambda = "ML")
#summary(model4_HFQ)
AIC(model4_HFQ)

#Body mass, temporality
model5_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Temporality + Ana_Max.Hz,
                   data = combined_HFQ, lambda = "ML")
summary(model5_HFQ)
AIC(model5_HFQ)
#model5_HFQ_t <- pgls(log(Max.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Temporality,
#                   data = combined_HFQ, lambda = "ML")
#summary(model5_HFQ_t)
#AIC(model5_HFQ_t)

#Body mass, sociality
model6_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Soc + Ana_Max.Hz, 
                   data = combined_HFQ, lambda = "ML")
#summary(model6_HFQ)
AIC(model6_HFQ)

#Body mass, life.history
model7_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Openness + Ana_Max.Hz, 
                   data = combined_HFQ, lambda = "ML")
#summary(model7_HFQ)
AIC(model7_HFQ)

#Body mass, temporality, sociality
model8_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Temporality +
                     Soc + Ana_Max.Hz, data = combined_HFQ, lambda = "ML")
#summary(model8_HFQ)
AIC(model8_HFQ)

#Body mass, temporality, life history
model9_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Temporality + 
                     Openness + Ana_Max.Hz, data = combined_HFQ, lambda = "ML")
#summary(model9_HFQ)
AIC(model9_HFQ)

#Body mass, sociality, life history
model10_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Soc + 
                      Openness + Ana_Max.Hz, data = combined_HFQ, lambda = "ML")
#summary(model10_HFQ)
AIC(model10_HFQ)

#Body mass, temporality, sociality, life-history
model11_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Temporality + Soc + Openness* + Ana_Max.Hz, 
                    data = combined_HFQ, lambda = "ML")
summary(model11_HFQ)
AIC(model11_HFQ)
anova(model11_HFQ)

#Temporality, sociality
model12_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ Temporality + Soc + 
                      Ana_Max.Hz, 
                    data = combined_HFQ, lambda = "ML")
#summary(model12_HFQ)
AIC(model12_HFQ)

#Temporality, life history
model13_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ Temporality + Openness+ 
                      Ana_Max.Hz, 
                    data = combined_HFQ, lambda = "ML")
#summary(model13_HFQ)
AIC(model13_HFQ)

#Sociality, life history
model14_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ Soc + Openness + Ana_Max.Hz, 
                    data = combined_HFQ, lambda = "ML")
#summary(model14_HFQ)
AIC(model14_HFQ)

#Temporality, sociality, life-history
model15_HFQ <- pgls(log(Max.Freq..kHz. + 1) ~ Temporality + Soc + Openness+ 
                      Ana_Max.Hz, 
                    data = combined_HFQ, lambda = "ML")
#summary(model15_HFQ)
AIC(model15_HFQ)


################# Dominant Frequency Models #################

combined_DFQ <- comparative.data(
  Squirrel_Tree, Freq_dom, names.col = "Scientific.Name", 
  vcv = TRUE, na.omit = FALSE)


#Null
model0_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ Ana_Max.Hz, data = combined_DFQ, lambda = "ML")
#summary(model0_DFQ)
AIC(model0_DFQ)

#Body mass
model1_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ log(Body.Size...Females..g.) + Ana_Max.Hz, data = combined_DFQ, lambda = "ML")
summary(model1_DFQ)
AIC(model1_DFQ)

#Temporality
model2_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ Temporality + Ana_Max.Hz, data = combined_DFQ, lambda = "ML")
#summary(model2_DFQ)
AIC(model2_DFQ)

#Sociality
model3_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ Soc + Ana_Max.Hz, data = combined_DFQ, lambda = "ML")
#summary(model3_DFQ)
AIC(model3_DFQ)

#Life History
model4_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ Openness + Ana_Max.Hz, data = combined_DFQ, lambda = "ML")
#summary(model4_DFQ)
AIC(model4_DFQ)

#Body mass, temporality
model5_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ log(Body.Size...Females..g.) + Temporality + Ana_Max.Hz, 
                   data = combined_DFQ, lambda = "ML")
summary(model5_DFQ)
AIC(model5_DFQ)

#Body mass, sociality
model6_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ log(Body.Size...Females..g.) + Soc + Ana_Max.Hz, 
                   data = combined_DFQ, lambda = "ML")
#summary(model6_DFQ)
AIC(model6_DFQ)

#Body mass, life.history
model7_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ log(Body.Size...Females..g.) + Openness + Ana_Max.Hz, 
                   data = combined_DFQ, lambda = "ML")
#summary(model7_DFQ)
AIC(model7_DFQ)

#Body mass, temporality, sociality
model8_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ log(Body.Size...Females..g.) +Temporality + Soc + Ana_Max.Hz, 
                   data = combined_DFQ, lambda = "ML")
#summary(model8_DFQ)
AIC(model8_DFQ)

#Body mass, temporality, life history
model9_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ log(Body.Size...Females..g.) + Temporality + Openness + Ana_Max.Hz, 
                   data = combined_DFQ, lambda = "ML")
summary(model9_DFQ)
AIC(model9_DFQ)

#Body mass, sociality, life history
model10_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ log(Body.Size...Females..g.) + Soc + Openness + Ana_Max.Hz, 
                    data = combined_DFQ, lambda = "ML")
#summary(model10_DFQ)
AIC(model10_DFQ)

#Body mass, temporality, sociality, life-history
model11_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ log(Body.Size...Females..g.) + Temporality + Soc + Openness + Ana_Max.Hz, 
                    data = combined_DFQ, lambda = "ML")
summary(model11_DFQ)
AIC(model11_DFQ)
anova(model11_DFQ)

#Temporality, sociality
model12_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ Temporality + Soc + Ana_Max.Hz, 
                    data = combined_DFQ, lambda = "ML")
#summary(model12_DFQ)
AIC(model12_DFQ)

#Temporality, life history
model13_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ Temporality + Openness + Ana_Max.Hz, 
                    data = combined_DFQ, lambda = "ML")
#summary(model13_DFQ)
AIC(model13_DFQ)

#Sociality, life history
model14_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ Soc + Openness + Ana_Max.Hz, data = combined_DFQ, lambda = "ML")
#summary(model14_DFQ)
AIC(model14_DFQ)

#Temporality, sociality, life-history
model15_DFQ <- pgls(log(Dom.Frequency..mean. + 1) ~ Temporality + Soc + Openness + Ana_Max.Hz, 
                    data = combined_DFQ, lambda = "ML")
#summary(model15_DFQ)
AIC(model15_DFQ)

################# Fundamental Frequency Models #################

combined_FFQ <- comparative.data(
  Squirrel_Tree, Freq_fun, names.col = "Scientific.Name", 
  vcv = TRUE, na.omit = FALSE)

#Null
model0_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ Ana_Max.Hz, data = combined_FFQ, lambda = "ML")
#summary(model0_FFQ)
AIC(model0_FFQ)

#Body mass
model1_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Ana_Max.Hz, data = combined_FFQ, lambda = "ML")
summary(model1_FFQ)
AIC(model1_FFQ)

#Temporality
model2_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ Temporality + Ana_Max.Hz, data = combined_FFQ, lambda = "ML")
#summary(model2_FFQ)
AIC(model2_FFQ)

#Sociality
model3_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ Soc + Ana_Max.Hz, data = combined_FFQ, lambda = "ML")
#summary(model3_FFQ)
AIC(model3_FFQ)

#Life History
model4_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ Openness + Ana_Max.Hz, data = combined_FFQ, lambda = "ML")
#summary(model4_FFQ)
AIC(model4_FFQ)

#Body mass, temporality
model5_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Temporality + Ana_Max.Hz, 
                   data = combined_FFQ, lambda = "ML")
summary(model5_FFQ)
AIC(model5_FFQ)

#Body mass, sociality
model6_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Soc + Ana_Max.Hz, 
                   data = combined_FFQ, lambda = "ML")
#summary(model6_FFQ)
AIC(model6_FFQ)

#Body mass, life.history
model7_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Openness + Ana_Max.Hz, 
                   data = combined_FFQ, lambda = "ML")
#summary(model7_FFQ)
AIC(model7_FFQ)

#Body mass, temporality, sociality
model8_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) +Temporality + Soc + Ana_Max.Hz, 
                   data = combined_FFQ, lambda = "ML")
#summary(model8_FFQ)
AIC(model8_FFQ)

#Body mass, temporality, life history
model9_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Temporality + Openness + Ana_Max.Hz, 
                   data = combined_FFQ, lambda = "ML")
summary(model9_FFQ)
AIC(model9_FFQ)

#Body mass, sociality, life history
model10_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Soc + Openness + Ana_Max.Hz, 
                    data = combined_FFQ, lambda = "ML")
#summary(model10_FFQ)
AIC(model10_FFQ)

#Body mass, temporality, sociality, life-history
model11_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ log(Body.Size...Females..g.) + Temporality + Soc + Openness + Ana_Max.Hz, 
                    data = combined_FFQ, lambda = "ML")
summary(model11_FFQ)
AIC(model11_FFQ)
anova(model11_FFQ)

#Temporality, sociality
model12_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ Temporality + Soc + Ana_Max.Hz, 
                    data = combined_FFQ, lambda = "ML")
#summary(model12_FFQ)
AIC(model12_FFQ)

#Temporality, life history
model13_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ Temporality + Openness + Ana_Max.Hz, 
                    data = combined_FFQ, lambda = "ML")
#summary(model13_FFQ)
AIC(model13_FFQ)

#Sociality, life history
model14_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ Soc + Openness + Ana_Max.Hz, data = combined_FFQ, lambda = "ML")
#summary(model14_FFQ)
AIC(model14_FFQ)

#Temporality, sociality, life-history
model15_FFQ <- pgls(log(Fundamental.Freq..kHz. + 1) ~ Temporality + Soc + Openness + Ana_Max.Hz, 
                    data = combined_FFQ, lambda = "ML")
#summary(model15_FFQ)
AIC(model15_FFQ)

################# Harmonic Frequency Models #################

combined_HMF <- comparative.data(
  Squirrel_Tree, Harm_max, names.col = "Scientific.Name", 
  vcv = TRUE, na.omit = FALSE)

#Null
model0_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ Ana_Max.Hz, data = combined_HMF, lambda = "ML")
#summary(model0_HMF)
AIC(model0_HMF)

#Body mass
model1_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ log(Body.Size...Females..g.) + Ana_Max.Hz, data = combined_HMF, lambda = "ML")
summary(model1_HMF)
AIC(model1_HMF)

#Temporality
model2_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ Temporality + Ana_Max.Hz, data = combined_HMF, lambda = "ML")
#summary(model2_HMF)
AIC(model2_HMF)

#Sociality
model3_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ Soc + Ana_Max.Hz, data = combined_HMF, lambda = "ML")
#summary(model3_HMF)
AIC(model3_HMF)

#Life History
model4_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ Openness + Ana_Max.Hz, data = combined_HMF, lambda = "ML")
#summary(model4_HMF)
AIC(model4_HMF)

#Body mass, temporality
model5_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ log(Body.Size...Females..g.) + Temporality + Ana_Max.Hz, 
                   data = combined_HMF, lambda = "ML")
summary(model5_HMF)
AIC(model5_HMF)

#Body mass, sociality
model6_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ log(Body.Size...Females..g.) + Soc + Ana_Max.Hz, 
                   data = combined_HMF, lambda = "ML")
summary(model6_HMF)
AIC(model6_HMF)

#Body mass, life.history
model7_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ log(Body.Size...Females..g.) + Openness + Ana_Max.Hz, 
                   data = combined_HMF, lambda = "ML")
summary(model7_HMF)
AIC(model7_HMF)

#Body mass, temporality, sociality
model8_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ log(Body.Size...Females..g.) +Temporality + Soc + Ana_Max.Hz, 
                   data = combined_HMF, lambda = "ML")
#summary(model8_HMF)
AIC(model8_HMF)

#Body mass, temporality, life history
model9_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ log(Body.Size...Females..g.) + Temporality + Openness + Ana_Max.Hz, 
                   data = combined_HMF, lambda = "ML")
summary(model9_HMF)
AIC(model9_HMF)

#Body mass, sociality, life history
model10_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ log(Body.Size...Females..g.) + Soc + Openness + Ana_Max.Hz, 
                    data = combined_HMF, lambda = "ML")
summary(model10_HMF)
AIC(model10_HMF)

#Body mass, temporality, sociality, life-history
model11_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ log(Body.Size...Females..g.) + Openness + Temporality + Soc + Ana_Max.Hz, 
                    data = combined_HMF, lambda = "ML")
summary(model11_HMF)
AIC(model11_HMF)
anova(model11_HMF)

#Temporality, sociality
model12_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ Temporality + Soc + Ana_Max.Hz, 
                    data = combined_HMF, lambda = "ML")
#summary(model12_HMF)
AIC(model12_HMF)

#Temporality, life history
model13_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ Temporality + Openness + Ana_Max.Hz, 
                    data = combined_HMF, lambda = "ML")
#summary(model13_HMF)
AIC(model13_HMF)

#Sociality, life history
model14_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ Soc + Openness + Ana_Max.Hz, data = combined_HMF, lambda = "ML")
#summary(model14_HMF)
AIC(model14_HMF)

#Temporality, sociality, life-history
model15_HMF <- pgls(log(Mean_of_Highest_Harmonic_Visible) ~ Temporality + Soc + Openness, 
                    data = combined_HMF, lambda = "ML")
#summary(model15_HMF)
AIC(model15_HMF)
################# Bandwidth Models #################

combined_BW <- comparative.data(
  Squirrel_Tree, Freq_BW, names.col = "Scientific.Name")

#Null
model0_BW <- pgls(log(BW) ~ 1, 
                  data = combined_BW, lambda = "ML")
summary(model0_BW)
AIC(model0_BW)
#Body mass
model1_BW <- pgls(log(BW) ~ log(Body.Size...Females..g.), 
                  data = combined_BW, lambda = "ML")
summary(model1_BW)
AIC(model1_BW)
#Temporality
model2_BW <- pgls(log(BW) ~ Temporality, 
                  data = combined_BW, lambda = "ML")
summary(model2_BW)
AIC(model2_BW)
#Sociality
model3_BW <- pgls(log(BW) ~ Soc, 
                  data = combined_BW, lambda = "ML")
summary(model3_BW)
AIC(model3_BW)
#Life History
model4_BW <- pgls(log(BW) ~ Life.History, 
                  data = combined_BW, lambda = "ML")
summary(model4_BW)
AIC(model4_BW)
#Body mass, temporality
model5_BW <- pgls(log(BW) ~ log(Body.Size...Females..g.) + 
                    Temporality, data = combined_BW, lambda = "ML")
summary(model5_BW)
AIC(model5_BW)
#Body mass, sociality
model6_BW <- pgls(log(BW) ~ log(Body.Size...Females..g.) + Soc,
                  data = combined_BW, lambda = "ML")
summary(model6_BW)
AIC(model6_BW)
#Body mass, life.history
model7_BW <- pgls(log(BW) ~ log(Body.Size...Females..g.) + 
                    Life.History, data = combined_BW, lambda = "ML")
summary(model7_BW)
AIC(model7_BW)
#Body mass, temporality, sociality
model8_BW <- pgls(log(BW) ~ log(Body.Size...Females..g.) + 
                    Temporality + Soc, data = combined_BW, lambda = "ML")
summary(model8_BW)
AIC(model8_BW)
#Body mass, temporality, life history
model9_BW <- pgls(log(BW) ~ log(Body.Size...Females..g.) + 
                    Temporality + Life.History, data = combined_BW, lambda = "ML")
summary(model9_BW)
AIC(model9_BW)
#Body mass, sociality, life history
model10_BW <- pgls(log(BW) ~ log(Body.Size...Females..g.) + Soc + 
                     Life.History, data = combined_BW, lambda = "ML")
summary(model10_BW)
AIC(model10_BW)
#Body mass, temporality, sociality, life-history
model11_BW <- pgls(log(BW) ~ log(Body.Size...Females..g.) + 
                     Temporality + Soc + Life.History + Year, data = combined_BW, lambda = "ML")
summary(model11_BW)
AIC(model11_BW)
#Temporality, sociality
model12_BW <- pgls(log(BW) ~ Temporality + Soc, 
                   data = combined_BW, lambda = "ML")
summary(model12_BW)
AIC(model12_BW)
#Temporality, life history
model13_BW <- pgls(log(BW) ~ Temporality + Life.History, 
                   data = combined_BW, lambda = "ML")
summary(model13_BW)
AIC(model13_BW)
#Sociality, life history
model14_BW <- pgls(log(BW) ~ Soc + Life.History, 
                   data = combined_BW, lambda = "ML")
summary(model14_BW)
AIC(model14_BW)
#Temporality, sociality, life-history
model15_BW <- pgls(log(BW) ~ Temporality + Soc + Life.History, 
                   data = combined_BW, lambda = "ML")
summary(model15_BW)
AIC(model15_BW)


################# Let's Plot These Models! #################
par(mfrow = c(2,2))
#Minimum
plot(model1_LFQ)
#Maximum
plot(model8_HFQ)
plot(model2_HFQ)
plot(model5_HFQ)
plot(model5_HFQ_t)
#Dominant
plot(model8_DFQ)
plot(model2_DFQ)
#Bandwidth
plot(model0_BW)
plot(model2_BW)
#Full models
plot(model11_LFQ)
plot(model11_HFQ)
plot(model11_DFQ)
plot(model11_BW)
par(mfrow = c(1,1))

################# Alarm Calls Only #################
Alarm_Types <- c("Aerial Predator", "Alarm", "Alarm; Mating; Locomotion", "Snake Mobbing",
                 "Terrestrial Predator")
Alarm_Types <- c("Aerial Predator", "Alarm", "Terrestrial Predator")
Alarm <- filter(Squirrels, Call.Type %in% Alarm_Types)

Alarm$Scientific.Name <- factor(Alarm$Scientific.Name)
Alarm_names <- unique(Alarm$Scientific.Name)
name.check(Squirrel_Tree, Alarm_names)

combined_LFQ <- comparative.data(
  Squirrel_Tree, Freq_min, names.col = "Scientific.Name")

Alarm$Sex[Alarm$Sex == "Unknown"] <- NA
Alarm$Sex[Alarm$Sex == "Both"] <- NA
Alarm$Sex <- factor(Alarm$Sex)

Alarm$Age[Alarm$Age == "All"] <- NA
Alarm$Age[Alarm$Age == "Pup"] <- NA
Alarm$Age[Alarm$Age == "Unknown"] <- NA
Alarm$Age[Alarm$Age == "Both"] <- NA
Alarm$Age <- factor(Alarm$Age)

library(lme4)
alarm1 <- lmer(log(Max.Freq..kHz.) ~ Sex + Age + Year + 
                 log(Body.Size...Females..g. + 1) + (1|Scientific.Name), data = Alarm)
alarm2 <- lmer(log(Dom.Frequency..mean.) ~ Sex + Age + Year + 
                 log(Body.Size...Females..g. + 1) + (1|Scientific.Name), data = Alarm)
library(car)
Anova(alarm1)

###### YOU ARE RIGHT HERE#####

plot(log(Max.Freq..kHz.) ~ log(Body.Size...Females..g.+ 1), data = Alarm)
plot(log(Dom.Frequency..mean.) ~ Year, data = Alarm)
plot(log(Max.Freq..kHz.) ~ Sex, data = Alarm)

ggplot(Alarm, aes(x = Year, y = value, color = variable)) +
  geom_point(aes(y = log(Max.Freq..kHz.), col = "Maximum")) +
  geom_smooth(aes(y = log(Max.Freq..kHz.), col = "Maximum"),
              method = "lm", se = F) +
  geom_point(aes(y = log(Dom.Frequency..mean.), col = "Dominant")) +
  geom_smooth(aes(y = log(Dom.Frequency..mean.), col = "Dominant"),
              method = "lm", se = F) +
  geom_point(aes(y = log(Min.Freq..kHz. + 1), col = "Minimum")) +
  geom_smooth(aes(y = log(Min.Freq..kHz. + 1), col = "Minimum"),
              method = "lm", se = F) +
  labs(color = "Frequency Type", x = "Year", y = "log(Frequency (kHz))") +
  theme_classic(base_size = 9) +
  theme(legend.position = c(0.90, 0.85),
        legend.background = element_rect(linetype = "solid", colour = "lightblue")) +
  scale_color_brewer(palette = "Dark2")
