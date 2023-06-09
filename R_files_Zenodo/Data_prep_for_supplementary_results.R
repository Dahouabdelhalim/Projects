#Combining all plots for the analysis on south Sweden Supplementary files
Reponsecurves_Suplementary_files_South <- ggarrange(Snowdaysplot, Coniferousplot, Humanfootprintplot, Openlandscapeplot, Soilmoistureplot, Waterplot,
                                                    labels = c("a)", "b)", "c)", "d)", "f)", "g"),
                                                    ncol = 3, nrow = 2)
Reponsecurves_Suplementary_files_South 
setwd("D:/Toegepaste Biologie Projecten/Stage/Modelling/R_MAXENT/South_maxent/plots")

ggsave(filename = "Responsecurves Maintext.jpeg",
       path = "D:/Toegepaste Biologie Projecten/Stage/Revision/Results/plots",
       units =c("cm"),
       width = 30,
       height = 14.6,
       dpi = 450)

#dataprep for plots
library(tidyverse)
library(ggplot2)
Coniferous1 <- read.csv('Mustela_putorius_0_coniferous_south.dat')
Coniferous2 <- read.csv('Mustela_putorius_1_coniferous_south.dat')
Coniferous3 <- read.csv('Mustela_putorius_2_coniferous_south.dat')
Coniferous4 <- read.csv('Mustela_putorius_3_coniferous_south.dat')
Coniferous5 <- read.csv('Mustela_putorius_4_coniferous_south.dat')
Coniferous6 <- read.csv('Mustela_putorius_5_coniferous_south.dat')
Coniferous7 <- read.csv('Mustela_putorius_6_coniferous_south.dat')
Coniferous8 <- read.csv('Mustela_putorius_7_coniferous_south.dat')
Coniferous9 <- read.csv('Mustela_putorius_8_coniferous_south.dat')
Coniferous10 <- read.csv('Mustela_putorius_9_coniferous_south.dat')



Humanfootprint1 <- read.csv('Mustela_putorius_0_humanfootprint_south.dat')
Humanfootprint2 <- read.csv('Mustela_putorius_1_humanfootprint_south.dat')
Humanfootprint3 <- read.csv('Mustela_putorius_2_humanfootprint_south.dat')
Humanfootprint4 <- read.csv('Mustela_putorius_3_humanfootprint_south.dat')
Humanfootprint5 <- read.csv('Mustela_putorius_4_humanfootprint_south.dat')
Humanfootprint6 <- read.csv('Mustela_putorius_5_humanfootprint_south.dat')
Humanfootprint7 <- read.csv('Mustela_putorius_6_humanfootprint_south.dat')
Humanfootprint8 <- read.csv('Mustela_putorius_7_humanfootprint_south.dat')
Humanfootprint9 <- read.csv('Mustela_putorius_8_humanfootprint_south.dat')
Humanfootprint10 <- read.csv('Mustela_putorius_9_humanfootprint_south.dat')

Snowdays1 <- read.csv('Mustela_putorius_0_snowdays_south.dat')
Snowdays2 <- read.csv('Mustela_putorius_1_snowdays_south.dat')
Snowdays3 <- read.csv('Mustela_putorius_2_snowdays_south.dat')
Snowdays4 <- read.csv('Mustela_putorius_3_snowdays_south.dat')
Snowdays5 <- read.csv('Mustela_putorius_4_snowdays_south.dat')
Snowdays6 <- read.csv('Mustela_putorius_5_snowdays_south.dat')
Snowdays7 <- read.csv('Mustela_putorius_6_snowdays_south.dat')
Snowdays8 <- read.csv('Mustela_putorius_7_snowdays_south.dat')
Snowdays9 <- read.csv('Mustela_putorius_8_snowdays_south.dat')
Snowdays10 <- read.csv('Mustela_putorius_9_snowdays_south.dat')

Water1 <- read.csv('Mustela_putorius_0_water_south.dat')
Water2 <- read.csv('Mustela_putorius_1_water_south.dat')
Water3 <- read.csv('Mustela_putorius_2_water_south.dat')
Water4 <- read.csv('Mustela_putorius_3_water_south.dat')
Water5 <- read.csv('Mustela_putorius_4_water_south.dat')
Water6 <- read.csv('Mustela_putorius_5_water_south.dat')
Water7 <- read.csv('Mustela_putorius_6_water_south.dat')
Water8 <- read.csv('Mustela_putorius_7_water_south.dat')
Water9 <- read.csv('Mustela_putorius_8_water_south.dat')
Water10 <- read.csv('Mustela_putorius_9_water_south.dat')

Openlandscape1 <- read.csv('Mustela_putorius_0_openlandscape_south.dat')
Openlandscape2 <- read.csv('Mustela_putorius_1_openlandscape_south.dat')
Openlandscape3 <- read.csv('Mustela_putorius_2_openlandscape_south.dat')
Openlandscape4 <- read.csv('Mustela_putorius_3_openlandscape_south.dat')
Openlandscape5 <- read.csv('Mustela_putorius_4_openlandscape_south.dat')
Openlandscape6 <- read.csv('Mustela_putorius_5_openlandscape_south.dat')
Openlandscape7 <- read.csv('Mustela_putorius_6_openlandscape_south.dat')
Openlandscape8 <- read.csv('Mustela_putorius_7_openlandscape_south.dat')
Openlandscape9 <- read.csv('Mustela_putorius_8_openlandscape_south.dat')
Openlandscape10 <- read.csv('Mustela_putorius_9_openlandscape_south.dat')

Soilmoisture1 <- read.csv('Mustela_putorius_0_soilmoisture_south.dat')
Soilmoisture2 <- read.csv('Mustela_putorius_1_soilmoisture_south.dat')
Soilmoisture3 <- read.csv('Mustela_putorius_2_soilmoisture_south.dat')
Soilmoisture4 <- read.csv('Mustela_putorius_3_soilmoisture_south.dat')
Soilmoisture5 <- read.csv('Mustela_putorius_4_soilmoisture_south.dat')
Soilmoisture6 <- read.csv('Mustela_putorius_5_soilmoisture_south.dat')
Soilmoisture7 <- read.csv('Mustela_putorius_6_soilmoisture_south.dat')
Soilmoisture8 <- read.csv('Mustela_putorius_7_soilmoisture_south.dat')
Soilmoisture9 <- read.csv('Mustela_putorius_8_soilmoisture_south.dat')
Soilmoisture10 <- read.csv('Mustela_putorius_9_soilmoisture_south.dat')

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

