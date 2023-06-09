#dataprep for plots

install.packages("ggpubr")

setwd("D:/Toegepaste Biologie Projecten/Stage/Revision/Modelling/Final Models/Seed1/Final_models/M_0.5_F_lqp_Set_147_E/plots")

library(tidyverse)
library(ggplot2)
library(ggpubr)
Coniferous1 <- read.csv('Mustela_putorius_0_coniferous.dat')
Coniferous2 <- read.csv('Mustela_putorius_1_coniferous.dat')
Coniferous3 <- read.csv('Mustela_putorius_2_coniferous.dat')
Coniferous4 <- read.csv('Mustela_putorius_3_coniferous.dat')
Coniferous5 <- read.csv('Mustela_putorius_4_coniferous.dat')
Coniferous6 <- read.csv('Mustela_putorius_5_coniferous.dat')
Coniferous7 <- read.csv('Mustela_putorius_6_coniferous.dat')
Coniferous8 <- read.csv('Mustela_putorius_7_coniferous.dat')
Coniferous9 <- read.csv('Mustela_putorius_8_coniferous.dat')
Coniferous10 <- read.csv('Mustela_putorius_9_coniferous.dat')



Humanfootprint1 <- read.csv('Mustela_putorius_0_humanfootprint.dat')
Humanfootprint2 <- read.csv('Mustela_putorius_1_humanfootprint.dat')
Humanfootprint3 <- read.csv('Mustela_putorius_2_humanfootprint.dat')
Humanfootprint4 <- read.csv('Mustela_putorius_3_humanfootprint.dat')
Humanfootprint5 <- read.csv('Mustela_putorius_4_humanfootprint.dat')
Humanfootprint6 <- read.csv('Mustela_putorius_5_humanfootprint.dat')
Humanfootprint7 <- read.csv('Mustela_putorius_6_humanfootprint.dat')
Humanfootprint8 <- read.csv('Mustela_putorius_7_humanfootprint.dat')
Humanfootprint9 <- read.csv('Mustela_putorius_8_humanfootprint.dat')
Humanfootprint10 <- read.csv('Mustela_putorius_9_humanfootprint.dat')

Snowdays1 <- read.csv('Mustela_putorius_0_snowdays.dat')
Snowdays2 <- read.csv('Mustela_putorius_1_snowdays.dat')
Snowdays3 <- read.csv('Mustela_putorius_2_snowdays.dat')
Snowdays4 <- read.csv('Mustela_putorius_3_snowdays.dat')
Snowdays5 <- read.csv('Mustela_putorius_4_snowdays.dat')
Snowdays6 <- read.csv('Mustela_putorius_5_snowdays.dat')
Snowdays7 <- read.csv('Mustela_putorius_6_snowdays.dat')
Snowdays8 <- read.csv('Mustela_putorius_7_snowdays.dat')
Snowdays9 <- read.csv('Mustela_putorius_8_snowdays.dat')
Snowdays10 <- read.csv('Mustela_putorius_9_snowdays.dat')

Water1 <- read.csv('Mustela_putorius_0_water.dat')
Water2 <- read.csv('Mustela_putorius_1_water.dat')
Water3 <- read.csv('Mustela_putorius_2_water.dat')
Water4 <- read.csv('Mustela_putorius_3_water.dat')
Water5 <- read.csv('Mustela_putorius_4_water.dat')
Water6 <- read.csv('Mustela_putorius_5_water.dat')
Water7 <- read.csv('Mustela_putorius_6_water.dat')
Water8 <- read.csv('Mustela_putorius_7_water.dat')
Water9 <- read.csv('Mustela_putorius_8_water.dat')
Water10 <- read.csv('Mustela_putorius_9_water.dat')

Openlandscape1 <- read.csv('Mustela_putorius_0_openlandscape.dat')
Openlandscape2 <- read.csv('Mustela_putorius_1_openlandscape.dat')
Openlandscape3 <- read.csv('Mustela_putorius_2_openlandscape.dat')
Openlandscape4 <- read.csv('Mustela_putorius_3_openlandscape.dat')
Openlandscape5 <- read.csv('Mustela_putorius_4_openlandscape.dat')
Openlandscape6 <- read.csv('Mustela_putorius_5_openlandscape.dat')
Openlandscape7 <- read.csv('Mustela_putorius_6_openlandscape.dat')
Openlandscape8 <- read.csv('Mustela_putorius_7_openlandscape.dat')
Openlandscape9 <- read.csv('Mustela_putorius_8_openlandscape.dat')
Openlandscape10 <- read.csv('Mustela_putorius_9_openlandscape.dat')

Soilmoisture1 <- read.csv('Mustela_putorius_0_soilmoisture.dat')
Soilmoisture2 <- read.csv('Mustela_putorius_1_soilmoisture.dat')
Soilmoisture3 <- read.csv('Mustela_putorius_2_soilmoisture.dat')
Soilmoisture4 <- read.csv('Mustela_putorius_3_soilmoisture.dat')
Soilmoisture5 <- read.csv('Mustela_putorius_4_soilmoisture.dat')
Soilmoisture6 <- read.csv('Mustela_putorius_5_soilmoisture.dat')
Soilmoisture7 <- read.csv('Mustela_putorius_6_soilmoisture.dat')
Soilmoisture8 <- read.csv('Mustela_putorius_7_soilmoisture.dat')
Soilmoisture9 <- read.csv('Mustela_putorius_8_soilmoisture.dat')
Soilmoisture10 <- read.csv('Mustela_putorius_9_soilmoisture.dat')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Merging bootstrapped runs

Coniferoustotal <- cbind(Coniferous1, Coniferous2$y, Coniferous3$y, Coniferous4$y, Coniferous5$y, Coniferous6$y, Coniferous7$y, Coniferous8$y, Coniferous9$y, Coniferous10$y)
Humanfootprinttotal <- cbind(Humanfootprint1, Humanfootprint2$y, Humanfootprint3$y, Humanfootprint4$y, Humanfootprint5$y, Humanfootprint6$y, Humanfootprint7$y, Humanfootprint8$y, Humanfootprint9$y, Humanfootprint10$y)
Openlandscapetotal <- cbind(Openlandscape1, Openlandscape2$y, Openlandscape3$y, Openlandscape4$y, Openlandscape5$y, Openlandscape6$y, Openlandscape7$y, Openlandscape8$y, Openlandscape9$y, Openlandscape10$y)
Snowdaystotal <- cbind(Snowdays1, Snowdays2$y, Snowdays3$y, Snowdays4$y, Snowdays5$y, Snowdays6$y, Snowdays7$y, Snowdays8$y, Snowdays9$y, Snowdays10$y)
Watertotal <- cbind(Water1, Water2$y, Water3$y, Water4$y, Water5$y, Water6$y, Water7$y, Water8$y, Water9$y, Water10$y)
Soilmoisturetotal <- cbind(Soilmoisture1, Soilmoisture2$y, Soilmoisture3$y, Soilmoisture4$y, Soilmoisture5$y, Soilmoisture6$y, Soilmoisture7$y, Soilmoisture8$y, Soilmoisture9$y, Soilmoisture10$y)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Naming columns

colnames(Coniferoustotal) <- c('variable', 'x', 'y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8', 'y9','y10')
colnames(Humanfootprinttotal) <- c('variable', 'x', 'y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8', 'y9','y10')
colnames(Openlandscapetotal) <- c('variable', 'x', 'y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8', 'y9','y10')
colnames(Snowdaystotal) <- c('variable', 'x', 'y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8', 'y9','y10')
colnames(Watertotal) <- c('variable', 'x', 'y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8', 'y9','y10')
colnames(Soilmoisturetotal) <- c('variable', 'x', 'y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8', 'y9','y10')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Adding average collumn

Coniferoustotal$Mean <- rowMeans(Coniferoustotal[,3:12])
Humanfootprinttotal$Mean <- rowMeans(Humanfootprinttotal[,3:12])
Openlandscapetotal$Mean <- rowMeans(Openlandscapetotal[,3:12])
Snowdaystotal$Mean <- rowMeans(Snowdaystotal[,3:12])
Watertotal$Mean <- rowMeans(Watertotal[,3:12])
Soilmoisturetotal$Mean <- rowMeans(Soilmoisturetotal[,3:12])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

