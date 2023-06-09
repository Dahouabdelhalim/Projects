### The arrival and spread of the European firebug Pyrrhocoris apterus 
### in Australia as documented by citizen scientists
### Luis Mata1, Blythe Vogel, Estibaliz Palma and Mallik Malapatil
### Code by Estibaliz Palma | v1 6 August 2021

PA_obs = read.csv("firebug_Australia.csv", header=TRUE, sep=",", na.strings=TRUE)

PA_obs <- PA_obs[PA_obs$quality_grade == "research",] # 98 records around Melbourne
PA_obs$year = as.numeric(substr(PA_obs$observed_on, 1, 4))
PA_obs$month = as.numeric(substr(PA_obs$observed_on, 6, 7))
PA_obs <- PA_obs[-c(19:20),] # remove records in the water

PA_obs_2018 <- PA_obs[PA_obs$year == 2018,] # 3 observations
PA_obs_2019 <- PA_obs[PA_obs$year == 2019,] # 16 observations
PA_obs_2020 <- PA_obs[PA_obs$year == 2020,] # 65 observations
PA_obs_2021 <- PA_obs[PA_obs$year == 2021,] # 12 observations


PA_obs_accumulated <- data.frame(year = c(rep(2018,1),rep(2019,12),rep(2020,12),rep(2021,7)), 
                                 month = c(12, rep(1:12, 2), 1:7),
                                 Nobs = c(nrow(PA_obs_2018),
                                          nrow(PA_obs_2019[PA_obs_2019$month==1,]),
                                          nrow(PA_obs_2019[PA_obs_2019$month==2,]),
                                          nrow(PA_obs_2019[PA_obs_2019$month==3,]),
                                          nrow(PA_obs_2019[PA_obs_2019$month==4,]),
                                          nrow(PA_obs_2019[PA_obs_2019$month==5,]),
                                          nrow(PA_obs_2019[PA_obs_2019$month==6,]),
                                          nrow(PA_obs_2019[PA_obs_2019$month==7,]),
                                          nrow(PA_obs_2019[PA_obs_2019$month==8,]),
                                          nrow(PA_obs_2019[PA_obs_2019$month==9,]),
                                          nrow(PA_obs_2019[PA_obs_2019$month==10,]),
                                          nrow(PA_obs_2019[PA_obs_2019$month==11,]),
                                          nrow(PA_obs_2019[PA_obs_2019$month==12,]),
                                          nrow(PA_obs_2020[PA_obs_2020$month==1,]),
                                          nrow(PA_obs_2020[PA_obs_2020$month==2,]),
                                          nrow(PA_obs_2020[PA_obs_2020$month==3,]),
                                          nrow(PA_obs_2020[PA_obs_2020$month==4,]),
                                          nrow(PA_obs_2020[PA_obs_2020$month==5,]),
                                          nrow(PA_obs_2020[PA_obs_2020$month==6,]),
                                          nrow(PA_obs_2020[PA_obs_2020$month==7,]),
                                          nrow(PA_obs_2020[PA_obs_2020$month==8,]),
                                          nrow(PA_obs_2020[PA_obs_2020$month==9,]),
                                          nrow(PA_obs_2020[PA_obs_2020$month==10,]),
                                          nrow(PA_obs_2020[PA_obs_2020$month==11,]),
                                          nrow(PA_obs_2020[PA_obs_2020$month==12,]),
                                          nrow(PA_obs_2021[PA_obs_2021$month==1,]),
                                          nrow(PA_obs_2021[PA_obs_2021$month==2,]),
                                          nrow(PA_obs_2021[PA_obs_2021$month==3,]),
                                          nrow(PA_obs_2021[PA_obs_2021$month==4,]),
                                          nrow(PA_obs_2021[PA_obs_2021$month==5,]),
                                          nrow(PA_obs_2021[PA_obs_2021$month==6,]),
                                          nrow(PA_obs_2021[PA_obs_2021$month==7,])) )
PA_obs_accumulated$Nobs_accumulated = NA
PA_obs_accumulated$Nobs_accumulated[1] <- PA_obs_accumulated$Nobs[1]
for (i in 2:nrow(PA_obs_accumulated)){
  PA_obs_accumulated$Nobs_accumulated[i] <- sum(PA_obs_accumulated$Nobs[1:i])
}

## Plot accumulated records from Dec2018 to July2021
jpeg(filename = "Accumulated_observations.jpg", width=1100, height=400)

  par(mar=c(4,4,4,4))
  my_bar <- barplot(PA_obs_accumulated$Nobs_accumulated, border=F, yaxt='n', las=2, ylim=c(0,100), col="dodgerblue")
  text(cex=0.75, x=my_bar+0.01, y=-2, paste(PA_obs_accumulated$month,PA_obs_accumulated$year), xpd=TRUE, srt=45, adj=1)
  points(my_bar, PA_obs_accumulated$Nobs_accumulated, type="l", lwd=5, col="darkorchid")
  axis(2, at=c(0,20,40,60,80,100), labels=c(0,20,40,60,80,100), las=1, cex.axis=2, padj=0.5)
  
dev.off()

## The first few records are curated to improve the accuracy of their location [they all should belong to Brimbank]
PA_obs$latitude[1:6] <- "-37.76662"
PA_obs$longitude[1:6] <- "144.8222"
PA_obs <- PA_obs[,c(3,22,23,32,36,37)]

PA_obs_2018 <- PA_obs[PA_obs$year == 2018,] # 3 observations
PA_obs_2019 <- PA_obs[PA_obs$year <= 2019,] # 19 observations
PA_obs_2020 <- PA_obs[PA_obs$year <= 2020,] # 84 observations
PA_obs_2021 <- PA_obs[PA_obs$year <= 2021,] # 96 observations

write.csv(PA_obs_2018, file="PA_obs_2018.csv")
write.csv(PA_obs_2019, file="PA_obs_2019.csv")
write.csv(PA_obs_2020, file="PA_obs_2020.csv")
write.csv(PA_obs_2021, file="PA_obs_2021.csv")

# I upload files 'PA_obs_2018.csv', 'PA_obs_2019.csv', 'PA_obs_2020.csv' and 'PA_obs_2021.csv' to QGIS
# Intersect them with layer 'LGA_POLYGON.shp' 
# Processing>ToolBox>Vector general>Join by location
# This generate objects 'PA_obs_2018_LGA.csv', 'PA_obs_2019_LGA.csv', 'PA_obs_2020_LGA.csv' and 'PA_obs_2021_LGA.csv'
# which contains all the accumulated individual records by year 2018/2019/2020/2021 and at which LGA each record was collected 

## Summarize the records of each year by LGA:
PA_obs_2018_LGA = read.csv("PA_obs_2018_LGA.csv", header=TRUE, sep=",", na.strings=TRUE)
PA_obs_2019_LGA = read.csv("PA_obs_2019_LGA.csv", header=TRUE, sep=",", na.strings=TRUE)
PA_obs_2020_LGA = read.csv("PA_obs_2020_LGA.csv", header=TRUE, sep=",", na.strings=TRUE)
PA_obs_2021_LGA = read.csv("PA_obs_2021_LGA.csv", header=TRUE, sep=",", na.strings=TRUE)

PA_obs_2018_LGA$occ <- 1
Heatmap_Nobs_2018_LGA <- aggregate(PA_obs_2018_LGA$occ ~ PA_obs_2018_LGA$LGA_NAME, FUN=sum)
names(Heatmap_Nobs_2018_LGA) <- c("LGA_NAME", "Nobs")
write.csv(Heatmap_Nobs_2018_LGA, file="Heatmap_Nobs_2018_LGA.csv")

PA_obs_2019_LGA$occ <- 1
Heatmap_Nobs_2019_LGA <- aggregate(PA_obs_2019_LGA$occ ~ PA_obs_2019_LGA$LGA_NAME, FUN=sum)
names(Heatmap_Nobs_2019_LGA) <- c("LGA_NAME", "Nobs")
write.csv(Heatmap_Nobs_2019_LGA, file="Heatmap_Nobs_2019_LGA.csv")

PA_obs_2020_LGA$occ <- 1
Heatmap_Nobs_2020_LGA <- aggregate(PA_obs_2020_LGA$occ ~ PA_obs_2020_LGA$LGA_NAME, FUN=sum)
names(Heatmap_Nobs_2020_LGA) <- c("LGA_NAME", "Nobs")
write.csv(Heatmap_Nobs_2020_LGA, file="Heatmap_Nobs_2020_LGA.csv")

PA_obs_2021_LGA$occ <- 1
Heatmap_Nobs_2021_LGA <- aggregate(PA_obs_2021_LGA$occ ~ PA_obs_2021_LGA$LGA_NAME, FUN=sum)
names(Heatmap_Nobs_2021_LGA) <- c("LGA_NAME", "Nobs")
write.csv(Heatmap_Nobs_2021_LGA, file="Heatmap_Nobs_2021_LGA.csv")

# I upload files 'Heatmap_Nobs_2018_LGA.csv', 'Heatmap_Nobs_2019_LGA.csv', 'Heatmap_Nobs_2020_LGA.csv' and 'Heatmap_Nobs_2021_LGA.csv' 
# files to QGIS, and use the column 'Nobs' to generate the color palette for the heatmaps

