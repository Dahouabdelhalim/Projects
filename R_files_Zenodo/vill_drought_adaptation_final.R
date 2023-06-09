# Analyses for, "Avoidance of seasonal drought in a neotropical understory herb"
# By: Julia Harenčár 
# February 23, 2022

setwd("/Users/Julia/Library/CloudStorage/GoogleDrive-jharenca@ucsc.edu/My Drive/GitHub/vill_drought_adaptation")
library(dplyr)
library(tidyverse)
library(ggplot2)
library(Hmisc)

#### Field Observations ####
# importing data and removing empty rows
trait <- read.csv("vill_vs_alle_210505.csv", header = T) %>% filter(!is.na(Date))
vill.trait <-trait %>% filter(Spp=='vill')
alle.trait <-trait %>% filter(Spp=='alle')

### stomatal density, pore length, and SPI
shapiro.test(vill.trait$stomatal.density) # not normal; W = 0.72726, p-value = 0.01183
shapiro.test(alle.trait$stomatal.density) # W = 0.85093, p-value = 0.1975
qqnorm(vill.trait$stomatal.density)
qqnorm(alle.trait$stomatal.density)
# wilcox.test(trait$stomatal.density~trait$Spp) # W = 20, p-value = 0.4102
t.test(trait$stomatal.density~trait$Spp) # t = -0.057385, df = 7.2793, p-value = 0.9558

#shapiro.test(trait$stomatal.pore.length- mean(trait$stomatal.pore.length, na.rm=T)) # normal; W = 0.97349, p-value = 0.9194
shapiro.test(vill.trait$stomatal.pore.length) # W = 0.90601, p-value = 0.4107
shapiro.test(alle.trait$stomatal.pore.length) # W = 0.84841, p-value = 0.1896
qqnorm(vill.trait$stomatal.pore.length)
qqnorm(alle.trait$stomatal.pore.length)
t.test(trait$stomatal.pore.length~trait$Spp) # t = -2.1398, df = 8.9529, p-value = 0.06119; means: A-0.0216 V-0.0262

shapiro.test(vill.trait$SPI) # W = 0.8632, p-value = 0.2004
shapiro.test(alle.trait$SPI) # W = 0.87428, p-value = 0.2842
qqnorm(vill.trait$SPI)
qqnorm(alle.trait$SPI)
t.test(trait$SPI~trait$Spp) # t = -2.2678, df = 7.0147, p-value = 0.05759 ; means A-0.039  V-0.056

shapiro.test(vill.trait$gmax) # W = 0.84967, p-value = 0.1564
shapiro.test(alle.trait$gmax) # W = 0.83326, p-value = 0.1471
qqnorm(vill.trait$gmax)
qqnorm(alle.trait$gmax)
# wilcox.test(trait$gmax~trait$Spp) # W = 14, p-value = 0.9307
t.test(trait$gmax~trait$Spp) # t = -0.68085, df = 7.2512, p-value = 0.5171

# #scatterplot of gmax against spi
# ggplot(trait, aes(x=gmax, y=SPI, color=Spp)) +
#   geom_point() +
#   theme_bw() +
#   scale_colour_manual('Species', labels= c(expression(italic("C. allenii")), 
#                                 expression(italic("C. villosissimus"))), 
#                       values = c("#009E73", "#E69F00")) +
#   theme(legend.text.align = 0) +
#   labs(x= bquote(~g[max]))

### delta 13 C
shapiro.test(vill.trait$d13C) # W = W = 0.8475, p-value = 0.1502
shapiro.test(alle.trait$d13C) # W = 0.89814, p-value = 0.363
qqnorm(vill.trait$d13C)
qqnorm(alle.trait$d13C)
# wilcox.test(trait$d13C~ trait$Spp) # W = 6, p-value = 0.06494
t.test(trait$d13C ~ trait$Spp) # t = -1.8714, df = 8.2019, p-value = 0.09729

### Leaf N
shapiro.test(vill.trait$leaf.N) # W = 0.95847, p-value = 0.8079
shapiro.test(alle.trait$leaf.N) # W = 0.92478, p-value = 0.5404
qqnorm(vill.trait$leaf.N)
qqnorm(alle.trait$leaf.N)
t.test(trait$leaf.N ~ trait$Spp) # t = -2.5547, df = 8.7165, p-value = 0.03172

### Stem density (sig)
shapiro.test(vill.trait$Stem.density) # W = 0.91364, p-value = 0.4608
shapiro.test(alle.trait$Stem.density) # W = 0.83586, p-value = 0.1204
qqnorm(vill.trait$Stem.density)
qqnorm(alle.trait$Stem.density)
# wilcox.test(trait$Stem.density ~ trait$Spp) # W = 31, p-value = 0.04113
t.test(trait$Stem.density ~ trait$Spp) # t = 2.7607, df = 9.921, p-value = 0.02025; means: alle = 0.13048169 vill = 0.09252183 

### Hydraulic Conductivity
# removing empty rows and dividing by area to get conductivity (rather than conductance)
hyd.cond <- trait %>% filter(!is.na(theo.kh)) %>% mutate(theo.kh = theo.kh/0.000002414)
vill.hyd.cond <-hyd.cond %>% filter(Spp=='vill')
alle.hyd.cond <-hyd.cond %>% filter(Spp=='alle')

shapiro.test(vill.hyd.cond$theo.kh) # W = 0.87205, p-value = 0.2345
shapiro.test(alle.hyd.cond$theo.kh) # W = 0.90849, p-value = 0.4744
qqnorm(vill.hyd.cond$theo.kh)
qqnorm(alle.hyd.cond$theo.kh)
# wilcox.test(hyd.cond$theo.kh ~ hyd.cond$Spp) # W = 9, p-value = 0.6095
t.test(hyd.cond$theo.kh ~ hyd.cond$Spp) # t = -1.1073, df = 7.3146, p-value = 0.3032; means: alle = 10.45360 vill = 16.94283 

### Rhizome water content
shapiro.test(vill.trait$Rhizome.water.content) # W = 0.9746, p-value = 0.9218
shapiro.test(alle.trait$Rhizome.water.content) # W = 0.91917, p-value = 0.4994
qqnorm(vill.trait$Rhizome.water.content) 
qqnorm(alle.trait$Rhizome.water.content)  
t.test(trait$Rhizome.water.content ~ trait$Spp) # t = 0.26028, df = 9.2048, p-value = 0.8004; means: alle = 86.68057 vill = 87.56916  

### visually checking for correlation
pairs(~SPI+d13C+leaf.N+Stem.density+theo.kh+Rhizome.water.content, data = trait)

# checking for correlations between SPI and other traits

cor(na.omit(trait$SPI),trait[c(1:5,7:12),c( 'd13C','leaf.N','Stem.density','Rhizome.water.content')], method = "pearson")
# highest Pearson CC = 0.608 (leaf.N and SPI)
# highest Spearman CC = 0.573 (leaf.N and SPI)

# (((gmax and d13C sig correlated with pearson only)))

# testing significance of correlation
rcorr(na.omit(trait$SPI),trait$leaf.N[c(1:5,7:12)], type = "spearman")
# pearson - significant; Pearson = 0.61; p=0.0473
# spearman - not significant; Spearman = 0.57; p = 0.0655

# checking for correlation between SPI and theo.kh
cor(na.omit(trait$SPI[2:12]), na.omit(trait$theo.kh), method = "pearson") 
# Pearson CC = 0.784
# Spearman CC = 0.673
# testing significance of correlation
rcorr(na.omit(trait$SPI[2:12]), na.omit(trait$theo.kh), type = "pearson")
# pearson significant\\; p = 0.0073
# spearman significant; p = 0.033

# checking for correlations between hydraulic conductivity and other traits 
cor(na.omit(trait$theo.kh),trait[c(1,3,4,6:12),c( 'd13C','leaf.N','Stem.density','Rhizome.water.content')], method = "spearman")
# highest Pearson CC = 0.75 (theo.kh & d13C)
# highest Spearman CC = 0.442 

# testing significance of correlation
rcorr(na.omit(trait$theo.kh), trait$d13C[c(1,3,4,6:12)], type = "spearman")
# pearson significant; CC = 0.75, p = 0.0119
# spearman significant; CC = 0.44, p = 0.20

# double checking for correlations in traits assessed in Avila-Lovera et al. 2022
cor(trait[,c( 'd13C','leaf.N','Stem.density','Rhizome.water.content')],trait[,c( 'd13C','leaf.N','Stem.density','Rhizome.water.content')], method = "pearson")
# highest pearson CC = 0.43
# highest Spearman CC = 0.4 

#### Greenhouse Observations ####
### Delta 13
dc <- read.csv("delta13.csv", header = T)
# calculate average of 2 values
dc <- dc %>% filter(spp %in% c("VV", "AA")) %>% mutate(mean_delta_PDB = (Delta_PDB_1 + Delta_PDB_2)/2)
vill.dc <- dc %>% filter(spp=="VV")
alle.dc <- dc %>% filter(spp=="AA")

shapiro.test(vill.dc$mean_delta_PDB)# W = 0.87753, p-value = 0.1222
shapiro.test(alle.dc$mean_delta_PDB)# W = 0.89086, p-value = 0.1734
qqnorm(vill.dc$mean_delta_PDB)
qqnorm(alle.dc$mean_delta_PDB)
# wilcox.test(dc$mean_delta_PDB ~ dc$spp)# W = 48, p-value = 0.9097 
t.test(dc$mean_delta_PDB ~ dc$spp) # t = -0.18873, df = 17.909, p-value = 0.8524

### LMA
lma <- read.csv("raw_chen_lma.csv", header = T)

# remove hybrid data, unnecessary columns, and one allenii value for equal sample sizes
AVlma <- lma[2:411,c('Category','mean_g_m2')] %>% 
  dplyr::filter(!Category == 'F2')
vill.LMA <- AVlma %>% filter(Category=="VV")
alle.LMA <- AVlma %>% filter(Category=="AA")

# AVlma %>% dplyr::filter(Category == "VV") %>% dplyr::select(mean_g_m2) #checking number of individuals
shapiro.test(vill.LMA$mean_g_m2) # W = 0.88778, p-value = 0.001886, NOT normal
shapiro.test(alle.LMA$mean_g_m2) # W = 0.98019, p-value = 0.7654
qqnorm(vill.LMA$mean_g_m2) 
qqnorm(alle.LMA$mean_g_m2) 
# wilcox.test(AVlma$mean_g_m2) # V = 2485, p-value = 3.637e-13
t.test(AVlma$mean_g_m2~AVlma$Category) # t = 5.4763, df = 67.205, p-value = 6.988e-07 means allenii = 51.875  vill = 37.915  

#### Growth Rate 
gr <- read.csv("AV_growth_rate.csv", header = T)
# removing hybrids and calculating growth rate as cm/day
# also removing rows with less than 14 days total growth)
gr <- gr %>% 
  filter(spp %in% c("alle", "vill")) %>% 
  filter(days > 14) %>% 
  mutate(GR_cm.day = total_growth/days) 
# rounding growth rate to 10 place (cm growth measured to 10 place)
gr$GR_cm.day <- round(gr$GR_cm.day, 1)
vill.gr <- gr %>% filter(spp=="vill")
alle.gr <- gr %>% filter(spp=="alle")

shapiro.test(vill.gr$GR_cm.day) # W = 0.9607, p-value = 0.2548
shapiro.test(alle.gr$GR_cm.day) # W = 0.86201, p-value = 0.00853
qqnorm(vill.gr$GR_cm.day)
qqnorm(alle.gr$GR_cm.day)
# wilcox.test(gr$GR_cm.day~gr$spp, exact=F) # W = 104.5, p-value = 2.046e-05 
t.test(gr$GR_cm.day~gr$spp) # t = -5.5588, df = 51.981, p-value = 9.495e-07

### Asat
AQ_results <- read.csv("A_Q_results2_10-500.csv", header=T)
vill.AQ <- AQ_results %>% filter(spp=="V")
alle.AQ <- AQ_results %>% filter(spp=="A")

shapiro.test(vill.AQ$A_sat) # W = 0.88269, p-value = 0.3217
shapiro.test(alle.AQ$A_sat) # W = 0.91366, p-value = 0.4899
qqnorm(alle.AQ$A_sat)
qqnorm(vill.AQ$A_sat)
t.test(AQ_results$A_sat ~ AQ_results$spp) # t = -2.2536, df = 7.3549, p-value = 0.05709

### checking for correlations not assessed in Avila-Lovera et. al. 2022
# use random data subsets when data lengths differ
# trait correlations with growth rate
# d13C
cor(gr$GR_cm.day[c(1:10,45:54)], dc$mean_delta_PDB, method = "pearson")
# not correlated; Spearmans CC = -0.002
# not correlated; Pearsons CC = -0.2245

# LMA
cor(gr$GR_cm.day[c(1:20,35:54)], AVlma$mean_g_m2[c(1:20, 51:70)], method = "pearson")
# not correlated; Spearmans CC = -0.137
# not correlated; Pearsons CC = -0.07841181

# Asat
cor(gr$GR_cm.day[c(1:5,50:54)], AQ_results$A_sat, method = "pearson")
# not correlated; Spearmans CC = -0.122
# not correlated; Pearsons CC = -0.1314389

# trait correlations with Asat
# LMA
cor(AQ_results$A_sat, AVlma$mean_g_m2[c(1:5, 66:70)], method = "pearson")
# not correlated; Spearmans CC = -.0.067
# not correlated; Pearson CC =-0.1374052

# d13C
cor(AQ_results$A_sat, dc$mean_delta_PDB[c(1:5, 11:15)], method = "pearson")
# not correlated; Spearmans CC = 0.103
# not correlated; Pearson CC = -0.0126

#### Boxplot plot panels ####
#load datasets
lma <- read.csv("raw_chen_lma.csv", header = T)
AVlma <- lma[,c('Category','mean_g_m2')] %>% filter(!Category == 'F2')
#trait <- read.csv("vill_vs_alle_210407.csv", header = T) %>% filter(!is.na(Date))
trait <- read.csv("vill_vs_alle_210505.csv", header = T) %>% filter(!is.na(Date))
# removing empty rows and dividing by area to get conductivity (rather than conductance)
hyd.cond <- trait %>% filter(!is.na(theo.kh)) %>% mutate(theo.kh = theo.kh/0.000002414)
gr <- read.csv("AV_growth_rate.csv", header = T)
gr <- gr %>% 
  filter(spp %in% c("alle", "vill")) %>% 
  mutate(GR_cm.day = total_growth/days) 
dc <- read.csv("delta13.csv", header = T)
dc <- dc %>% filter(spp %in% c("VV", "AA")) %>% mutate(mean_delta_PDB = (Delta_PDB_1 + Delta_PDB_2)/2)
AQ_results <- read.csv("A_Q_results2_10-500.csv", header=T)

## Field traits
# preparing for plotting
dev.off()

pdf("Field_trait_boxplots.pdf", 
    width = 8, height = 2.5, 
    bg = "white") 

par(mfrow = c(1, 4))
# setting inner margins to allow space for different axes and inner margins to uniform 1
par(mar = c(3, 5, 0, 0), oma = c(1, 1, 1, 1)) 

# boxplot of leaf N
boxplot(trait$leaf.N ~ trait$Spp, 
        ylab = bquote("Leaf N (%) "),
        cex.lab = 1.3,
        las = 1,
        col = c("#009E73", "#E69F00"),
        names = c("", ""))

# boxplot of stem density
boxplot(trait$Stem.density ~ trait$Spp, 
        ylab = bquote("Stem density" ~(g~""~cm^ -3 )), 
        cex.lab = 1.3,
        las = 1,
        col = c("#009E73", "#E69F00"),
        names = c("", ""))

# boxplot of stomata pore size index
boxplot(trait$SPI~trait$Spp, 
        ylab = "SPI ", 
        cex.lab = 1.3,
        las = 1,
        col = c("#009E73", "#E69F00"),
        names = c("", ""))

# boxplot of hydraulic conductivity  
boxplot(hyd.cond$theo.kh ~ hyd.cond$Spp, 
        ylab = bquote(~k[h*","*theo]~ ~(kg~""~m^ -1~MPa^ -1~s^ -1)), 
        cex.lab = 1.3,
        las = 1,
        col = c("#009E73", "#E69F00"),
        names = c("", ""),
        ylim = c(0, 35))

dev.off()

# boxplot of stomata gmax for the supplement
pdf("gmax_boxplot.pdf", 
    width = 4, height = 3.5, 
    bg = "white") 
# setting inner margins to allow space for different axes and inner margins to uniform 1
par(mar = c(3, 5, 0, 0), oma = c(1, 1, 1, 1)) 

boxplot(trait$gmax~trait$Spp, 
        ylab = bquote(~g[max]~""~(mmol~""~m^ -2~s^ -1)),
        cex.lab = 1.3,
        las = 1,
        col = c("#009E73", "#E69F00"),
        names = c("", ""))
dev.off()

## Greenhouse traits
dev.off()
pdf("GH_trait_boxplots.pdf", 
    width = 8, height = 2.5, 
    bg = "white") 

par(mfrow = c(1, 4))

# setting inner margins to allow space for different axes and inner margins to uniform 1
par(mar = c(3, 5, 0, 0), oma = c(1, 1, 1, 1)) 

# boxplot of growth rate
#quartz(width=4.4, height=5) #from before multipanel 
boxplot(gr$GR_cm.day ~ gr$spp, 
        ylab = bquote("Growth rate" ~(cm~""~day^ -1)), 
        cex.lab = 1.3,
        las = 1,
        col = c("#009E73", "#E69F00"),
        names = c("", ""),
        ylim = c(0,1.43))

# boxplot of LMA 
boxplot(AVlma$mean_g_m2~AVlma$Category, 
        ylab = bquote("LMA"~(g~""~cm^ -2)), 
        cex.lab = 1.3,
        las = 1,
        col = c("#009E73", "#E69F00"),
        names = c("", ""))

# boxplot of light saturated PSN 
boxplot(AQ_results$A_sat ~ AQ_results$spp, 
        ylab = bquote(~A[sat]~""~(µmol~""~m^ -2~s^ -1)), 
        cex.lab = 1.3,
        las = 1,
        col = c("#009E73", "#E69F00"),
        names = c("", ""),
        ylim = c(0,8.5))

# boxplot of delta13
boxplot(dc$mean_delta_PDB ~ dc$spp, 
        cex.lab = 1.3,
        ylab = expression(paste(delta^{13}, "C (\\u2030)")), 
        las = 1,
        col = c("#009E73", "#E69F00"),
        names = c("", ""))

dev.off()

#### Dry Down ####
## plotting leaves over time
# importing data 
library(readxl) 
library(tidyverse)
# get sheet names for iteration in loop
sheet_names <- excel_sheets("dry_down_full_dataset.xls")
# initiate empty dataframe for loop to fill
dry_down <- data.frame()

# for loop to generate single data frame from data across all sheets
for (i in 1:length(sheet_names)){
  temp <- read_excel("dry_down_full_dataset.xls", sheet = i) %>% # read in sheet i
    select(number:comments) %>% # select only cols of interest
    mutate(date=sheet_names[i]) %>%  # add a col with the name of the data sheet as the date for all data from the sheet
    mutate_all(as.character) # convert all data to character for joining purposes
  dry_down <- bind_rows(dry_down, temp) # append sheet i to the previous sheets in a single df
}

#adjusting df structure
dry_down$`total leaf` <- as.numeric(dry_down$`total leaf`)
dry_down$nosh <- as.numeric(dry_down$nosh)
dry_down$green <- as.numeric(dry_down$green)
dry_down$date <- as.Date(dry_down$date,"%b %d %Y") 

# converting date to number of days since last hydration
# last hydration was on July 30
dry_down$days_LH <- dry_down$date - as.Date("2007-07-30")

# adding rehydration treatment info column for plotting purposes
library("stringr")

for(i in 1:nrow(dry_down)){
  if(is.na(dry_down$comments[i]==TRUE)){
    dry_down$rehydration[i] <- "none"
  }else{
    dry_down$rehydration[i] <- str_sub(dry_down$comments[i], 1, 2)
  }
}

# Remove leaked on plants 
dry_down <- dry_down %>% filter(!number %in% c('19', '42', '46', '157', '154', '191', '194', '170', '192', '195'))

# remove all rehydration data after 2 months (incomplete data after that duration)
# calculate day number that is 60 days after rehydration
W1_60 <- (as.Date("2008-01-14") - as.Date("2007-07-30")) + 60
W2_60 <- (as.Date("2008-04-21") - as.Date("2007-07-30")) + 60
# removing incomplete data 60+ days past rehydration
dry_down <- dry_down %>% filter(!(rehydration == "W2" & days_LH > W1_60 ))
dry_down <- dry_down %>% filter(!(rehydration == "W3" & days_LH > W2_60 ))

# get sample sizes
temp <- dry_down %>% 
  select(c(cross, date, green, comments, nosh)) %>% 
  filter(cross %in% c('A','V'), date=='2007-06-02', !is.na(green), is.na(comments), !is.na(nosh))
table(temp$cross) # 45 of each

# rehydration sample sizes
temp <- dry_down %>% 
  select(c(cross, date, green, rehydration, nosh)) %>% 
  filter(cross %in% c('V'), date=="2008-06-08", !is.na(green), !is.na(nosh))
table(temp$rehydration)

# summarizing the green leaf data 
GreenSummary <- dry_down %>% 
  select(c(cross, nosh, green, date, comments, rehydration, days_LH)) %>% 
  filter(cross %in% c('A','V'), date!=23, !is.na(green)) %>% 
  group_by(days_LH, cross, rehydration) %>% 
  summarise(mean_green = mean(green), sd_green = sd(green))

#plotting green leaf data
library(ggplot2)
leaves <- ggplot(GreenSummary, aes(x=days_LH, y=mean_green, shape=cross, colour=cross)) +
  geom_line(aes(linetype=rehydration)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean_green-sd_green, ymax=mean_green+sd_green),
                width = 5, position=position_dodge(0.05)) +
  theme_bw() +
  theme(legend.position = c(0.87, 0.87), legend.margin=margin(t=-0.5, r=0, b=0.4, l=0, unit= "cm")) +
  labs(x = "Days since last watering", y = "Number of green leaves", color = "", shape ="", linetype="") +
  scale_colour_manual(labels= c(expression(italic("C. allenii")), 
                                expression(italic("C. villosissimus"))), 
                      values = c("#009E73", "#E69F00")) +
  scale_shape_discrete(labels= c(expression(italic("C. allenii")), 
                                 expression(italic("C. villosissimus")))) +
  scale_linetype_manual(values = c("solid", "twodash", "twodash")) +
  guides(linetype = FALSE) +
  theme(legend.text.align = 0) 


#### Comparison of differences to broader Costus variance ####
setwd("/Users/juliaharencar/Google Drive (jharenca@ucsc.edu)/GitHub/vill_drought_adaptation")
library(dplyr)
library(tidyverse)

# importing data, selecting relevant coloumns/traits, and removing alle and vill
full.dat <- read.csv('dataset_climate_soil_traits.csv', header = T)
dat <- full.dat %>% 
  select(Site, Spp, SLA, Stem.density, Rhizome.water.content, d13C, leaf.N) # %>%
filter(!Spp %in% c("vill", "alle"))

dat$Spp<-factor(dat$Spp)
dat$Site<-factor(dat$Site)

# calculate species means for all spp within a site
spp.avs <- dat %>% 
  group_by(Site, Spp) %>% 
  summarise(mean_LMA = mean(1/SLA*10000),
            mean_Stem.density = mean(Stem.density),
            mean_Rhizome.water.content = mean(Rhizome.water.content),
            mean_d13C = mean(d13C),
            mean_leaf.N = mean(leaf.N))

# calculate std dev per site of spp means
site_stdv <- spp.avs %>% 
  group_by(Site) %>% 
  summarise(stdv_LMA = sd(mean_LMA, na.rm = T),
            stdv_Stem.density = sd(mean_Stem.density, na.rm = T),
            stdv_Rhizome.water.content = sd(mean_Rhizome.water.content, na.rm = T),
            stdv_d13C = sd(mean_d13C, na.rm = T),
            stdv_leaf.N = sd(mean_leaf.N, na.rm = T))

write.table(site_stdv, file = "per_site_stdev.csv", sep=",")

plot(dat$Stem.density)

