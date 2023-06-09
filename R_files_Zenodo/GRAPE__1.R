

##############################
##### Processing Raw Data #####
##############################

d2 <- read.csv("raw_data_07.04.2021.csv") 

##### Creating Across-year Averages USING YRS 1, 2 AND 3 "tave.dat" - skipping total acidity
tave.dat <- with(d2, data.frame(vineyard, pinot, mngt_int))
tave.dat$Nave <- with(d2, (Grape_N_15+Grape_N_16+Grape_N_17)/3)
tave.dat$Tartave <- with(d2, (Tart_15+Tart_16+Tart_17)/3)
tave.dat$Malicave <- with(d2, (Malic_15+Malic_16+Malic_17)/3)
tave.dat$sugave <- with(d2, (sug_15+sug_16+sug_17)/3) 

##### Adjust data scales for lavaan
tave.dat$sclNave = with(tave.dat, Nave/100)
tave.dat$sclsugave = with(tave.dat, sugave/10)

## Averaging pH values involves additional steps
# step 1: convert pH values to hydrogen ion values
# concentration of H+ = 10^(-X) , where X is the pH
d2$H1 = with(d2, 10^(-pH_15))
d2$H2 = with(d2, 10^(-pH_16))
d2$H3 = with(d2, 10^(-pH_17))
d2$Have = with(d2, (H1+H2+H3)/3)
# step 2: convert Have to pHave and place in tave.dat
# pH = -log10(concentration of H+)
tave.dat$pHave = with(d2, -log10(Have))


### Temporal correlations of environmental variables ####
library("PerformanceAnalytics")
chart.Correlation(d2[, 24:38], use = "pairwise.complete.obs", histogram=TRUE, pch=19, method = "pearson")

# Soil Nitrogen
N.dat <- subset(d2, select = soil_N_15:soil_N_17)
chart.Correlation(N.dat, histogram=TRUE, pch=19, method = "pearson")

#plantsdiv
plantsdiv.dat <- subset(d2, select = plantsdiv_15:plantsdiv_17)
chart.Correlation(plantsdiv.dat, histogram=TRUE, pch=19, method = "pearson")  

#plantscov
plantscov.dat <- subset(d2, select = plantscov_15:plantscov_17)
chart.Correlation(plantscov.dat, histogram=TRUE, pch=19, method = "pearson")  

#plantscov
plantscov.dat <- subset(d2, select = plantscov_15:plantscov_17)
chart.Correlation(plantscov.dat, histogram=TRUE, pch=19, method = "pearson")  

#covNfix
covNfix.dat <- subset(d2, select = covNfix_15:covNfix_17)
chart.Correlation(covNfix.dat, histogram=TRUE, pch=19, method = "pearson")  


##### Creating Across-year Averages "tave.dat" - skip total acidity


# create plant cover variable
tave.dat$pltcovave = with(d2, (plantscov_15+plantscov_16+plantscov_17)/3)
str(tave.dat)

# create plant diversity variables
tave.dat$pltdivave = with(d2, (plantsdiv_15+plantsdiv_16+plantsdiv_17)/3)

# create Nfix cover variable
tave.dat$covNfixave = with(d2, (covNfix_15+covNfix_16+covNfix_17)/3)
str(tave.dat)

# create relNfix cover variable
tave.dat$soilNave = with(d2, (soil_N_15+soil_N_16+soil_N_17)/3)
str(tave.dat)


## Recode vars to roughly same scale
tave.dat$sclpltcovave <- tave.dat$pltcovave/100
tave.dat$sclcovNfixave <- tave.dat$covNfixave/10


# Check variances
varTable(tave.dat)

out.dat <- with(tave.dat, data.frame( vineyard, pinot, mngt_int, sclNave, Tartave, Malicave, sclsugave, pltdivave, soilNave, sclcovNfixave))
labels <- c("Site", "Pinot", "MgtInten", "N", "Tart", "Malic",  "Sugars", "PltRich", "SoilN", "NfixCov") 
colnames(out.dat) <- labels

write.csv(out.dat, "avg_data_archive_fix.csv", row.names=F)

pub.dat <- read.csv("avg_data_archive_fix.csv")
names(pub.dat)
head(pub.dat)

