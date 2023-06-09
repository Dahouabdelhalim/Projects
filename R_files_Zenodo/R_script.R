# Script by Camilla Fløjgaard
# Supporting Information for manuscript: 
# Title: Exploring a natural baseline for large-herbivore biomass in ecological restoration

# Load the data file in the Supporting information and recreate analyses, tables and figures
# from the manuscript

# Large herbivore biomass and NPP. 

rm(list=ls(all=TRUE))
#insert path
setwd()  
library(readxl)
library(raster)
library(data.table)
library(ggplot2)
library(rworldmap)

# NPP:
# Run this bit first time to get the raster with NAs
npp1 <- raster("MOD17A3_Science_NPP_mean_00_15.tif")
NAvalue(npp1) <- 65535
# From readme.txt file here: http://files.ntsg.umt.edu/data/NTSG_Products/MOD17/GeoTIFF/MOD17A3/
# Multiply with 0.1 to get g Carbon/m2/yr. 
npp <- npp1*0.1
plot(npp)
writeRaster(npp,'MOD17A3_Science_NPP_mean_00_15_NA.tif', overwrite=TRUE)
# ...or load raster: 
npp <- raster("MOD17A3_Science_NPP_mean_00_15_NA.tif")

# Load biomass data:
MyData = data.frame(read_xlsx('data file.xlsx', sheet="Sheet1")) 
dim(MyData) # 302 15
MyData <- MyData[-which(gsub("([A-Za-z]+).*", "\\\\1", MyData$Note)=="Excluded"),]
dim(MyData) # 288 15
length(unique(MyData$Name_Ecosystem_new)) #145

# Make geodesic buffers:
# Run this bit first time, then save the table. 
pts <- MyData[,c(5,4)]
make_GeodesicBuffer <- function(pts, width) {
  
  # A) Construct buffers as points at given distance and bearing ---------------
  
  dg <- seq(from = 0, to = 360, by = 5)
  
  # Construct equidistant points defining circle shapes (the "buffer points")
  buff.XY <- geosphere::destPoint(p = pts, 
                                  b = rep(dg, each = length(pts)), 
                                  d = width)
  
  # B) Make SpatialPolygons -------------------------------------------------
  
  # Group (split) "buffer points" by id
  buff.XY <- as.data.frame(buff.XY)
  id  <- rep(1:dim(pts)[1], times = length(dg))
  lst <- split(buff.XY, id)
  
  # Make SpatialPolygons out of the list of coordinates
  poly   <- lapply(lst, sp::Polygon, hole = FALSE)
  polys  <- lapply(list(poly), sp::Polygons, ID = NA)
  spolys <- sp::SpatialPolygons(Srl = polys, 
                                proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  # Disaggregate (split in unique polygons)
  spolys <- sp::disaggregate(spolys)
  return(spolys)
}

pts_buf_1km <- make_GeodesicBuffer(as.matrix(pts), width = 1*10^3)
pts_buf_5km <- make_GeodesicBuffer(as.matrix(pts), width = 5*10^3)
pts_buf_10km <- make_GeodesicBuffer(as.matrix(pts), width = 10*10^3)
pts_buf_50km <- make_GeodesicBuffer(as.matrix(pts), width = 50*10^3)
pts_buf_100km <- make_GeodesicBuffer(as.matrix(pts), width = 100*10^3)

npp.table <- data.frame(pts)
npp.table$mean_1km <- extract(npp, pts_buf_1km, fun=mean, na.rm=TRUE, df=TRUE)[,2]
npp.table$mean_5km <- extract(npp, pts_buf_5km, fun=mean, na.rm=TRUE, df=TRUE)[,2]
npp.table$mean_10km <- extract(npp, pts_buf_10km, fun=mean, na.rm=TRUE, df=TRUE)[,2]
npp.table$mean_50km <- extract(npp, pts_buf_50km, fun=mean, na.rm=TRUE, df=TRUE)[,2]
npp.table$mean_100km <- extract(npp, pts_buf_100km, fun=mean, na.rm=TRUE, df=TRUE)[,2]

npp.table$median_1km <- extract(npp, pts_buf_1km, fun=median, na.rm=TRUE, df=TRUE)[,2]
npp.table$median_5km <- extract(npp, pts_buf_5km, fun=median, na.rm=TRUE, df=TRUE)[,2]
npp.table$median_10km <- extract(npp, pts_buf_10km, fun=median, na.rm=TRUE, df=TRUE)[,2]
npp.table$median_50km <- extract(npp, pts_buf_50km, fun=median, na.rm=TRUE, df=TRUE)[,2]
npp.table$median_100km <- extract(npp, pts_buf_100km, fun=median, na.rm=TRUE, df=TRUE)[,2]

#write the table and skip the above next time: 
#write.table(npp.table, file="npp_table.txt",sep="\\t")
npp.table <- read.csv("npp_table.txt",sep="\\t", header=TRUE)

MyData2 <- cbind(MyData, npp.table[,c(3:12)])
MyData2$X20kg_TotalHerbivoreBiomass_kgkm2 <- as.numeric(MyData2$X20kg_TotalHerbivoreBiomass_kgkm2)
names(MyData2)[7] <- "Min20kgTotalHerbivoreBiomass_kgkm2"
MyData2[,16:25]=MyData2[,16:25]*1000#converting NPP from g/m2 to kg/km2


#### Aggregate by Name_ecosystem_new: ####
MyDatamean <- aggregate(.~Name_Ecosystem_new+Continent, MyData2[,c(2:5,6,16:25)], FUN = mean, na.action=na.pass)
dim(MyDatamean) #145 15

########################

Result <- data.frame(matrix(ncol = 6, nrow = 90))
names(Result) <- c("Continent","r2","p", "intercept_log","intercept", "slope")
for (i in c(1:10)){
  lm1 <- lm(log(MyDatamean$TotalHerbivoreBiomass_kgkm2)~log(MyDatamean[,5+i]))
  x <- summary(lm1)
  Result$r2[i] <- x$r.squared
  Result$p[i] <- x$coefficients[2,4]
  Result$intercept_log[i] <- x$coefficients[1,1]
  Result$intercept[i] <- exp(Result$intercept_log[i])
  Result$slope[i] <- x$coefficients[2,1]
  Result$Continent[i] <- "Global"
}

continent <- c("Africa", "Asia", "Europe", "NorthAmerica", "SouthAmerica")
for (j in c(1:5)){
  MyDatamean2 <- MyDatamean[which(MyDatamean$Continent==continent[j]),]
  for (i in c(1:10)){
    lm1 <- lm(log(MyDatamean2$TotalHerbivoreBiomass_kgkm2)~log(MyDatamean2[,5+i]))
    x <- summary(lm1)
    Result$r2[i+(j*10)] <- x$r.squared
    Result$p[i+(j*10)] <- x$coefficients[2,4]
    Result$intercept_log[i+(j*10)] <- x$coefficients[1,1]
    Result$intercept[i+(j*10)] <- exp(Result$intercept_log[i+(j*10)])
    Result$slope[i+(j*10)] <- x$coefficients[2,1]
    Result$Continent[i+(j*10)] <- continent[j]
  }
}

# lm for Africa using only data from >=2000:
MyDataAf <- MyData2[MyData2$Continent=="Africa",]
MyDataAf$Year <- as.numeric(MyDataAf$Year)
MyDataAf2000 <- MyDataAf[!is.na(MyDataAf$Year),]
MyDataAf2000 <- MyDataAf2000[MyDataAf2000$Year>=2000,]
MyDataAf2000mean <- aggregate(.~Name_Ecosystem_new+Continent, MyDataAf2000[,c(2:5,6,16:25)], FUN = mean, na.action=na.pass)
for (i in c(1:10)){
  lm1 <- lm(log(MyDataAf2000mean$TotalHerbivoreBiomass_kgkm2)~log(MyDataAf2000mean[,5+i]))
  x <- summary(lm1)
  Result$r2[i+60] <- x$r.squared
  Result$p[i+60] <- x$coefficients[2,4]
  Result$intercept_log[i+60] <- x$coefficients[1,1]
  Result$intercept[i+60] <- exp(Result$intercept_log[i+60])
  Result$slope[i+60] <- x$coefficients[2,1]
  Result$Continent[i+60] <- "Africa >=2000"
}


# lm for Africa using only ecosystems with megaherbivores:
MyDataAfmega <- MyDataAf[which(MyDataAf$Megafauna>=1),]
str(MyDataAfmega)
MyDataAfmega$Megafauna <- as.numeric(MyDataAfmega$Megafauna)
MyDataAfmega <- MyDataAfmega[-which(is.na(MyDataAfmega$Megafauna)),]
MyDataAfmegamean <- aggregate(.~Name_Ecosystem_new+Continent, MyDataAfmega[,c(2:5,6,16:25)], FUN = mean, na.action=na.pass)
for (i in c(1:10)){
  lm1 <- lm(log(MyDataAfmegamean$TotalHerbivoreBiomass_kgkm2)~log(MyDataAfmegamean[,5+i]))
  x <- summary(lm1)
  Result$r2[i+70] <- x$r.squared
  Result$p[i+70] <- x$coefficients[2,4]
  Result$intercept_log[i+70] <- x$coefficients[1,1]
  Result$intercept[i+70] <- exp(Result$intercept_log[i+70])
  Result$slope[i+70] <- x$coefficients[2,1]
  Result$Continent[i+70] <- "Africa megaherbivores"
}

# lm for Africa using herbivores > 20 kg: 
str(MyData2)
MyDataAf20 <- MyDataAf[!is.na(MyDataAf$Min20kgTotalHerbivoreBiomass_kgkm2),]
MyDataAf20mean <- aggregate(.~Name_Ecosystem_new+Continent, MyDataAf20[,c(2:5,7,16:25)], FUN = mean, na.action=na.pass)
for (i in c(1:10)){
  lm1 <- lm(log(MyDataAf20mean$Min20kgTotalHerbivoreBiomass_kgkm2)~log(MyDataAf20mean[,5+i]))
  x <- summary(lm1)
  Result$r2[i+80] <- x$r.squared
  Result$p[i+80] <- x$coefficients[2,4]
  Result$intercept_log[i+80] <- x$coefficients[1,1]
  Result$intercept[i+80] <- exp(Result$intercept_log[i+80])
  Result$slope[i+80] <- x$coefficients[2,1]
  Result$Continent[i+80] <- "Africa min 20 kg"
}


# find the best model for each continent:
Result$NPP <- cbind(as.vector(rep(c("meanNPP1km", "meanNPP5km","meanNPP10km", "meanNPP50km", "meanNPP100km", 
  "medianNPP1km", "medianNPP5km","medianNPP10km", "medianNPP50km", "medianNPP100km"), times=9)))
Result <- Result[,c(1,7,2,3,4,5,6)]
Result_summary1 <- Result[which(Result$NPP=="meanNPP1km"),]
Result_summary2 <- setDT(Result)[,.SD[which.max(r2)],keyby=Continent]
Result_summary2 <- Result_summary2[which(Result_summary2$p<=0.05),]
Result_summary <- rbind(Result_summary1,Result_summary2)
table(MyDatamean$Continent)
Result_summary$n <- NA
Result_summary$n[Result_summary$Continent=="Global"]<- 146
Result_summary$n[Result_summary$Continent=="Africa"]<- 48
Result_summary$n[Result_summary$Continent=="Asia"]<- 37
Result_summary$n[Result_summary$Continent=="Europe"]<- 17
Result_summary$n[Result_summary$Continent=="NorthAmerica"]<- 32
Result_summary$n[Result_summary$Continent=="SouthAmerica"]<- 11
Result_summary$n[Result_summary$Continent=="Africa >=2000"]<- dim(MyDataAf2000mean)[1]
Result_summary$n[Result_summary$Continent=="Africa megaherbivores"]<- dim(MyDataAfmegamean)[1]
Result_summary$n[Result_summary$Continent=="Africa min 20 kg"]<- dim(MyDataAf20mean)[1]

Result_summary$CurrentBiomassmean <- NA
Result_summary$CurrentBiomasssd <- NA
Result_summary$CurrentBiomassmean[which(Result_summary$Continent== "Africa")] <- 
  mean(MyDatamean$TotalHerbivoreBiomass_kgkm2[which(MyDatamean$Continent== "Africa")])
Result_summary$CurrentBiomasssd[which(Result_summary$Continent== "Africa")] <-
  sd(MyDatamean$TotalHerbivoreBiomass_kgkm2[which(MyDatamean$Continent== "Africa")])
Result_summary$CurrentBiomassmean[which(Result_summary$Continent== "Global")]<- 
  mean(MyDatamean$TotalHerbivoreBiomass_kgkm2)
Result_summary$CurrentBiomasssd[which(Result_summary$Continent== "Global")]<- 
  sd(MyDatamean$TotalHerbivoreBiomass_kgkm2)
Result_summary$CurrentBiomassmean[which(Result_summary$Continent== "Europe")] <- 
  mean(MyDatamean$TotalHerbivoreBiomass_kgkm2[which(MyDatamean$Continent== "Europe")])
Result_summary$CurrentBiomasssd[which(Result_summary$Continent== "Europe")]<- 
  sd(MyDatamean$TotalHerbivoreBiomass_kgkm2[which(MyDatamean$Continent== "Europe")])
Result_summary$CurrentBiomassmean[which(Result_summary$Continent== "NorthAmerica")] <- 
  mean(MyDatamean$TotalHerbivoreBiomass_kgkm2[which(MyDatamean$Continent== "NorthAmerica")])
Result_summary$CurrentBiomasssd[which(Result_summary$Continent== "NorthAmerica")]<- 
  sd(MyDatamean$TotalHerbivoreBiomass_kgkm2[which(MyDatamean$Continent== "NorthAmerica")])
Result_summary$CurrentBiomassmean[which(Result_summary$Continent== "SouthAmerica")] <- 
  mean(MyDatamean$TotalHerbivoreBiomass_kgkm2[which(MyDatamean$Continent== "SouthAmerica")])
Result_summary$CurrentBiomasssd[which(Result_summary$Continent== "SouthAmerica")]<- 
  sd(MyDatamean$TotalHerbivoreBiomass_kgkm2[which(MyDatamean$Continent== "SouthAmerica")])
Result_summary$CurrentBiomassmean[which(Result_summary$Continent== "Asia")] <- 
  mean(MyDatamean$TotalHerbivoreBiomass_kgkm2[which(MyDatamean$Continent== "Asia")])
Result_summary$CurrentBiomasssd[which(Result_summary$Continent== "Asia")]<- 
  sd(MyDatamean$TotalHerbivoreBiomass_kgkm2[which(MyDatamean$Continent== "Asia")])
Result_summary$CurrentBiomassmean[which(Result_summary$Continent== "Africa >=2000")] <- mean(MyDataAf2000mean$TotalHerbivoreBiomass_kgkm2)
Result_summary$CurrentBiomasssd[which(Result_summary$Continent== "Africa >=2000")]<- 
  sd(MyDataAf2000mean$TotalHerbivoreBiomass_kgkm2)
Result_summary$CurrentBiomassmean[which(Result_summary$Continent== "Africa megaherbivores")] <- 
  mean(MyDataAfmegamean$TotalHerbivoreBiomass_kgkm2)
Result_summary$CurrentBiomasssd[which(Result_summary$Continent== "Africa megaherbivores")]<-
  sd(MyDataAfmegamean$TotalHerbivoreBiomass_kgkm2)
Result_summary$CurrentBiomassmean[which(Result_summary$Continent== "Africa min 20 kg")] <- 
  mean(MyDataAf20mean$Min20kgTotalHerbivoreBiomass_kgkm2)
Result_summary$CurrentBiomasssd[which(Result_summary$Continent== "Africa min 20 kg")]<- 
  sd(MyDataAf20mean$Min20kgTotalHerbivoreBiomass_kgkm2)


write.table(Result, file="modelcomparison.txt",sep="\\t") 
write.table(Result_summary, file="modelcomparison_summary.txt",sep="\\t") 


### make some plots... 

MyDatamean$Continent[MyDatamean$Continent=="SouthAmerica"] <- "S. America"
MyDatamean$Continent[MyDatamean$Continent=="NorthAmerica"] <- "N. America"
H = 5.1
W = 6

### log-log with regression lines for continents:
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73","grey", "#0072B2", "#D55E00")
MyDatamean$signif <- "ns"
MyDatamean$signif[which(MyDatamean$Continent %in% c("Africa", "N. America"))] <- "sig"
MyDatamean$signif <- factor(MyDatamean$signif, levels=c("sig", "ns"))

pdf("BiomassNPP_CI_1km.pdf", height = H, width=W)
ggplot(NULL) +
  geom_point(data=MyDatamean, aes(x=mean_1km, y=TotalHerbivoreBiomass_kgkm2, color=Continent)) + 
  scale_colour_manual(values=cbbPalette)+
  geom_smooth(data=MyDatamean, aes(x=mean_1km, y=TotalHerbivoreBiomass_kgkm2,color= Continent, linetype=as.factor(signif)),
              method="lm", se=TRUE, alpha=0.2,fullrange=FALSE)+
  geom_smooth(data=MyDatamean, aes(x=mean_1km, y=TotalHerbivoreBiomass_kgkm2, color="Global"),
              method="lm", se=TRUE,alpha=0.2, fullrange=FALSE)+
  theme(legend.position = "bottom",panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.title=element_blank(),legend.key = element_rect(fill = "white"))+
  scale_x_continuous(trans='log', breaks = c(50000,500000, 1000000, 2000000), name = expression(NPP ~ (kg ~ C ~ km^{-2} ~ yr^{-1}))) +
  scale_y_continuous(trans='log', breaks = c(10, 100, 250,500,1000,2500,5000,10000), name = expression(Large ~ herbivore ~ biomass ~ (kg ~ km^{-2})))+
  coord_cartesian(xlim=c(50000,2200000))
dev.off()


#### overview of biomass: 
MyDatamean$new <- MyDatamean$Continent
MyDatamean$new[which(MyDatamean$Name_Ecosystem_new %in% 
                       c("Faia Brava","Kostilkovo","Groenlanden","Millingerwaard","Lika Plains","Oostvaardersplassen","Mols laboratory"))] <- "RewEu"
MyDatamean$NPP <- "low NPP"
MyDatamean$NPP[which(MyDatamean$mean_1km>500000)] <- "medium NPP"
MyDatamean$NPP[which(MyDatamean$mean_1km>1000000)] <- "high NPP"

MyDatamean$NPP <- factor(MyDatamean$NPP, c("low NPP", "medium NPP", "high NPP"))
MyDatamean$new <- factor(MyDatamean$new, c("Africa", "Asia", "Europe", "RewEu", "N. America", "S. America"))
library(cowplot)
theme_set(theme_cowplot())
pdf("biomass by NPP.pdf", height = 4, width=10)
ggplot(data = MyDatamean, aes(x=new, y=TotalHerbivoreBiomass_kgkm2)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.3))+
  geom_point(aes(y=TotalHerbivoreBiomass_kgkm2, group=new), position = position_dodge(width=0.75))+
  background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
  panel_border("black")+ # and a border around each panel
  ylab(expression(Large ~ herbivore ~ biomass ~ (kg ~ km^{-2})))+facet_wrap(~NPP,ncol = 3)+
  xlab("")
dev.off()

# Additional tables and figures - Figure S1 - Map fo NPP and study sites: 
nppkgkm <- npp*1000
#writeRaster(nppkgkm,'NPPkgkm', overwrite=TRUE)
nppkgkm <- raster("NPPkgkm.tif")

H = 4
W = 7
pdf("S5_NPPandSites.pdf", height = H, width=W)
par(bg=NA,mar=c(6,6,2,6),oma=c(0,0,0,1))
# add axes titles: xlab="Longitude", ylab="Latitude"?
plot(nppkgkm, axes=F, xlab="Longitude", ylab="Latitude",
     legend.args = list(text = expression(NPP ~ (kg ~ C ~ km^{-2} ~ yr^{-1}))))
axis(1, at=c(-150,-100, -50, 0, 50, 100, 150), 
     labels=c(expression(paste(-150*degree)),expression(paste(-100*degree)),
              expression(paste(-50*degree)),expression(paste(0*degree)),
              expression(paste(50*degree)),expression(paste(100*degree)),
              expression(paste(150*degree))))
axis(2, at=c(-100, -50, 0, 50, 100), 
     labels=c(expression(paste(-100*degree)),expression(paste(-50*degree)),
              expression(paste(0*degree)),expression(paste(50*degree)),
              expression(paste(100*degree))))

points(MyData[,c(5,4)], pch=19, cex=0.3)
dev.off()


