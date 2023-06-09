# Delineate ecotone based on synchrony analyses and extract statistics about it

rm(list=ls())

library(raster)
library(wsyn)
library(rgdal)
library(maptools)
library(geosphere)
library(nlme)
library(dplyr)
library(fields)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

certaintyfiles <- list.files("../Outputs/Intermediate", pattern="certainty")
certaintyfiles <- certaintyfiles[!grepl(".aux", certaintyfiles)]
certaintyfiles <- certaintyfiles[!grepl("allYrs", certaintyfiles)]
dem <- raster("../Data/10n090w_20101117_gmted_med075.tif")
certaintyAllYrs <- raster("../Outputs/Intermediate/certaintyRasterallYrs.tif")
groundtruth <- read.csv("../Data/forestPlots_ACG.csv")

groundtruth$pch <- NA
groundtruth$pch[groundtruth$forestType=="dry"] <- 1
groundtruth$pch[groundtruth$forestType=="transition"] <- 2
groundtruth$pch[groundtruth$forestType=="wet"] <- 3
groundtruth$pch[groundtruth$forestType=="cloud"] <- 3

thresh <- 0.2 #below this threshold = ecotone
rcl <- c(-Inf, thresh, 1,
         thresh, Inf, NA)
ecotoneAllYrs <- reclassify(certaintyAllYrs, rcl)


guan <- readOGR("../Data/guanacaste.shp")
prj <- "+proj=lcc +datum=NAD27 +ellps=clrk66 +lat_0=10.4667 +lat_1=9.9333 +lat_2=11.0000 +lon_0=-84.3333 +x_0=500000.000 + y_0=271820.522 +units=m"
crs(guan) <- CRS(prj) 
guan <-spTransform(guan, CRS(proj4string(ecotoneAllYrs)), use_ob_tran=TRUE)


lsib <- readOGR("../Data/Department_of_State_Large-Scale_International_Boundary_28LSIB29/Department of State Large-Scale International Boundary (LSIB).shp")
str(lsib@data)
lsib <- lsib[lsib$Name == "COSTA RICA",]
lsib <- spTransform(lsib, CRS(proj4string(ecotoneAllYrs)), use_ob_tran=TRUE)

#loop over years ----------------------------------------------------------------------------------
smooth <- FALSE

#initialize outputs
area_km2 <- NULL
perim_km <- NULL
elev_med <- NULL
elev_range <- NULL


for(ii in 1:length(certaintyfiles)){
  
  certainty.ii <- raster(paste0("../Outputs/Intermediate/", certaintyfiles[ii]))
  
  if(smooth){
    rcl <- c(-Inf, thresh, 1,
             thresh, Inf, 0)
    ecotone.rast.ii <- reclassify(certainty.ii, rcl)
    plot(ecotone.rast.ii, main=certaintyfiles[ii])
    ecotone.rast.sm <- focal(ecotone.rast.ii, w=matrix(1, 3, 3), fun=mean, na.rm=T)
    ecotone.rast.ii <- reclassify(ecotone.rast.sm, rcl=matrix(c(-Inf,1/2,NA,1/2,Inf,1)))
    plot(ecotone.rast.ii, main=certaintyfiles[ii])
    basename <- paste0("ecotone", substr(certaintyfiles[ii], 16, 19), "t", thresh*100,"sm")
  }
  else{
    rcl <- c(-Inf, thresh, 1,
             thresh, Inf, NA)
    ecotone.rast.ii <- reclassify(certainty.ii, rcl)
    basename <- paste0("ecotone", substr(certaintyfiles[ii], 16, 19), "t", thresh*100)
  }
  
  if(ii==1){
    ecotone.nyears <- ecotone.rast.ii
    ecotone.nyears[is.na(ecotone.rast.ii)] <- 0
  }
  else{
    tmp <- ecotone.rast.ii
    tmp[is.na(ecotone.rast.ii)] <- 0
    ecotone.nyears <- ecotone.nyears + tmp
  }
  
  ecotone.poly.ii <- rasterToPolygons(ecotone.rast.ii, dissolve=TRUE)
  #ecotone.poly.ii <- disaggregate(ecotone.poly.ii)
  
  ecotone.poly.ii$area_km2 <- raster::area(ecotone.poly.ii)/1e6
  ecotone.poly.ii$perim_km <- perimeter(ecotone.poly.ii)/1000

  
  dem.ii <- crop(dem, ecotone.poly.ii)
  dem.ii <- mask(dem.ii, ecotone.poly.ii)
  
  
  # writeOGR(ecotone.poly.ii, dsn="../Outputs/Ecotones", layer=basename, driver="ESRI Shapefile"
  #          , overwrite_layer = TRUE)
  # writeRaster(ecotone.rast.ii, paste0("../Outputs/Ecotones/",basename,".tif"), format="GTiff")
  
  area_km2 <- c(area_km2, sum(ecotone.poly.ii$area_km2))
  perim_km <- c(perim_km, sum(ecotone.poly.ii$perim_km))
  elev_med <- c(elev_med, median(values(dem.ii), na.rm=T))
  elev_range <- c(elev_range, quantile(values(dem.ii), 0.75, na.rm=TRUE) - 
                    quantile(values(dem.ii), 0.25, na.rm=TRUE))
  
}


ecotone.nyears <- mask(ecotone.nyears, guan)
dem2 <- crop(dem, guan)
dem2 <- mask(dem2, guan)

png("../Outputs/fig_ecotone_nyears.png", width=6.5, height=6.5, units="in", res=300)

plot(ecotone.nyears)
contour(dem2, add=T, levels=c(200,300,400,500,600,800,1000,1200,1400))

dev.off()


png("../Outputs/FigSx_validation.png", width=6.5, height=5.7, units="in", res=300)
par(mar=c(2.1,2.1,1.1,1.1))
plot(ecotone.nyears)
points(groundtruth$longitude, groundtruth$latitude, pch=groundtruth$pch)
legend("topright", pch=1:3, legend=c("dry","ecotone","rain"), bty="n", ncol=1)

#inset map
par(fig=c(0.65, 0.84, 0.075, 0.24), mar=c(0.1,0.1,0.1,0.1), new=TRUE)
plot(lsib, col="grey", border=NA, bty="o")
plot(guan, col="black", border=NA, add=TRUE)


dev.off()



## look at ecotone in particularly wet/dry years

ecotone.dry <- raster(paste0("../Outputs/Intermediate/", "certaintyRaster","2001",".tif"))
ecotone.wet <- raster(paste0("../Outputs/Intermediate/", "certaintyRaster","2008",".tif"))

rcl <- c(-Inf, thresh, 1,
         thresh, Inf, NA)

ecotone.dry <- reclassify(ecotone.dry, rcl)
ecotone.wet <- reclassify(ecotone.wet, rcl)


png("../Outputs/fig_Sx_dryYear_vs_wetYear.png", units="in", width=6.5, height=3.5, res=300)

par(mar=c(2.5,2.5,1.9,1), mfrow=c(1,2))

plot(ecotone.dry, col="darkgoldenrod", legend=FALSE)
plot(guan, add=TRUE)
mtext("Dry year (2001)", line=0.2)

plot(ecotone.wet, col="springgreen3", legend=FALSE)
plot(guan, add=TRUE)
mtext("Wet year (2008)", line=0.2)

dev.off()





## Test for trends in ecotone characteristics


ecotone.stats <- data.frame(year = 2000:2021
                            ,area_km2 = area_km2
                            ,perim_area = perim_km/area_km2
                            ,elev_med = elev_med
                            ,elev_range = elev_range
)


plot(ecotone.stats$year, ecotone.stats$elev_med)
summary(gls(elev_med ~ year, data=ecotone.stats, correlation=corAR1(form=~1)))

plot(ecotone.stats$year, ecotone.stats$elev_range)
summary(gls(elev_range ~ year, data=ecotone.stats, correlation=corAR1(form=~1)))

plot(ecotone.stats$year, ecotone.stats$area_km2)
summary(gls(area_km2 ~ year, data=ecotone.stats, correlation=corAR1(form=~1)))

plot(ecotone.stats$year, ecotone.stats$perim_area)
summary(gls(perim_area ~ year, data=ecotone.stats, correlation=corAR1(form=~1)))


## test for relationships with climate

climate <- read.csv("../Data/guanClimate.csv")

df <- left_join(ecotone.stats, climate, by=c("year" = "years"))
df$totalPr_m1 <- c(NA, df$totalPr[1:(nrow(df)-1)])
df$meanPET_m1 <- c(NA, df$meanPET[1:(nrow(df)-1)])
df$drySeasonPr_m1 <- c(NA, df$drySeasonPr[1:(nrow(df)-1)])
df$mei_m1 <- c(NA, df$mei[1:(nrow(df)-1)])

cormat <- cor(df, use="pairwise.complete.obs") #use correlation to filter what tests we want to do
cormat.sub <- cormat[c(6,7,9:14), 2:5]

pal <- colorRampPalette(colors=c("red","white","blue"))

coords <- expand.grid(seq(0,1,length.out=8),seq(0,1,length.out=4))


png("../Outputs/fig2_correlations_rough.png", width=6.5, height=4, units="in", res=300)
par(mar=c(9.1,6.1,1.1,1.1))
image.plot(cormat.sub, zlim=c(-1,1), col=pal(50), xaxt="n", yaxt="n")
axis(2, at=seq(0,1,length.out=4), labels=c("Area","Perim:Area","Median Elev.","Elev. range"), las=2)
axis(1, at=seq(0,1,length.out=8), labels=c("Total PPT","Mean PET","Dry Season PPT","MEI","Total PPT(t-1)",
                                           "Mean PET(t-1)","Dry Season PPT(t-1)","MEI(t-1)"), las=2)
text(x=coords[,1], y=coords[,2], labels=as.character(round(cormat.sub,2)))
dev.off()

# summary(lm(perim_area ~ meanPET, data=df)) #not significant
# summary(lm(perim_area ~ mei, data=df)) #possibly significant
# summary(lm(elev_med ~ meanPET, data=df)) #possibly significant
# summary(lm(elev_med ~ totalPr, data=df)) #possibly significant
# summary(lm(elev_med ~ mei, data=df)) #possibly significant
# summary(lm(elev_range ~ totalPr, data=df)) #possibly significant
# summary(lm(elev_range ~ meanPET, data=df)) #likely significant

# plot(lm(perim_area ~ mei, data=df)) #looks reasonable
# plot(lm(elev_med ~ meanPET, data=df))
# plot(df$meanPET, df$elev_med)


# acf(residuals(lm(perim_area ~ mei, data=df)))
# acf(residuals(lm(elev_med ~ meanPET, data=df)))
# acf(residuals(lm(elev_med ~ totalPr, data=df)))
# acf(residuals(lm(elev_med ~ mei, data=df)))
# acf(residuals(lm(elev_range ~ totalPr, data=df)))    
# acf(residuals(lm(elev_range ~ meanPET, data=df)))    
# 
# 
# hist(residuals(lm(perim_area ~ mei, data=df)))
# hist(residuals(lm(elev_med ~ meanPET, data=df)))
# hist(residuals(lm(elev_med ~ totalPr, data=df)))
# hist(residuals(lm(elev_med ~ mei, data=df)))
# hist(residuals(lm(elev_range ~ totalPr, data=df)))    
# hist(residuals(lm(elev_range ~ meanPET, data=df)))    


library(nlme)

df2 <- df[complete.cases(df),]

#do model selection to compare lagged vs non-lagged

#area
summary(gls(area_km2 ~ mei, data=df2, correlation=corAR1(form = ~1))) #significant, AIC = 129.41 <-----------
summary(gls(area_km2 ~ totalPr_m1, data=df2, correlation=corAR1(form = ~1))) #significant, AIC = 143.00 <-----------
summary(gls(area_km2 ~ meanPET_m1, data=df2, correlation=corAR1(form = ~1))) #significant, AIC = 126.62 <-----------
summary(gls(area_km2 ~ mei_m1, data=df2, correlation=corAR1(form = ~1))) #not quite significant, AIC = 130.82


#perim:area -- none with correlation > 0.3
summary(gls(perim_area ~ drySeasonPr_m1, data=df2, correlation=corAR1(form = ~1))) #not significant

#median elevation
summary(gls(elev_med ~ meanPET, data=df2, correlation=corAR1(form = ~1))) #almost significant, AIC = 169.90
summary(gls(elev_med ~ mei, data=df2, correlation=corAR1(form = ~1))) #not significant, AIC = 166.50
summary(gls(elev_med ~ totalPr_m1, data=df2, correlation=corAR1(form = ~1))) #not significant, AIC = 179.67
summary(gls(elev_med ~ meanPET_m1, data=df2, correlation=corAR1(form = ~1))) #significant, AIC = 169.14 <-----------


#elevation range

summary(gls(elev_range ~ totalPr, data=df2, correlation=corAR1(form = ~1))) #not significant, AIC = 195.42
summary(gls(elev_range ~ meanPET, data=df2, correlation=corAR1(form = ~1))) #significant, AIC=181.84 <-----------
summary(gls(elev_range ~ drySeasonPr, data=df2, correlation=corAR1(form = ~1))) #almost significant, AIC = 191.91
summary(gls(elev_range ~ meanPET_m1, data=df2, correlation=corAR1(form = ~1))) #not significant, AIC=187.09




png("../Outputs/fig3_sigEffects_rough.png", width=6.5, height=4, units="in", res=300)

par(mfrow=c(2,3), mar=c(4.1,4.1,1.1,1.1))

plot(df2$mei, df2$area_km2, pch=19, xlab="MEI", ylab="Area")
abline(gls(area_km2 ~ mei, data=df2, correlation=corAR1(form = ~1)), col="grey")
mtext("a)", at=par("usr")[1] - 0.03*diff(par("usr")[1:2]), cex=2/3, line=0.1)

plot(df2$totalPr_m1, df2$area_km2, pch=19, xlab="Total PPT(t-1)", ylab="Area")
abline(gls(area_km2 ~ totalPr_m1, data=df2, correlation=corAR1(form = ~1)), col="grey")
mtext("b)", at=par("usr")[1] - 0.03*diff(par("usr")[1:2]), cex=2/3, line=0.1)

plot(df2$meanPET_m1, df2$area_km2, pch=19, xlab="Mean PET(t-1)", ylab="Area")
abline(gls(area_km2 ~ meanPET_m1, data=df2, correlation=corAR1(form = ~1)), col="grey")
mtext("c)", at=par("usr")[1] - 0.03*diff(par("usr")[1:2]), cex=2/3, line=0.1)

plot(df2$meanPET_m1, df2$elev_med, pch=19, xlab="Mean PET(t-1)", ylab="Median elevation")
abline(gls(elev_med ~ meanPET_m1, data=df2, correlation=corAR1(form = ~1)), col="grey")
mtext("d)", at=par("usr")[1] - 0.03*diff(par("usr")[1:2]), cex=2/3, line=0.1)

plot(df2$meanPET, df2$elev_range, pch=19, xlab="Mean PET", ylab="Elevation range")
abline(gls(elev_range ~ meanPET, data=df2, correlation=corAR1(form = ~1)), col="grey")
mtext("e)", at=par("usr")[1] - 0.03*diff(par("usr")[1:2]), cex=2/3, line=0.1)

dev.off()




png("../Outputs/figSx_climateTimeSeries.png", width=3.25, height=5.5, units="in", res=300)

par(mfrow=c(4,1), mar=c(2.1,4.1,1.1,1.1), mgp=c(2.1,0.9,0), xpd=T)

plot(df2$year, df2$totalPr, type="b", ylab=expression(Total~PPT~"("*kg~m^-2*")"), xlab="")
plot(df2$year, df2$meanPET, type="b", ylab=expression(Mean~PET~"("*kg~m^-2*")"), xlab="")
plot(df2$year, df2$drySeasonPr, type="b", ylab=expression(Dry~Season~PPT~"("*kg~m^-2*")"), xlab="")
plot(df2$year, df2$mei, type="b", ylab="MEI", xlab="")

dev.off()




pdf("../Outputs/Fig3_scatterplots_sigEnvEffects.pdf", onefile=TRUE)

plot(df2$meanPET, df2$elev_med, pch=16, xlab="Annual mean PET", ylab="Median elevation of ecotone")
plot(df2$totalPr, df2$elev_med, pch=16, xlab="Annual total ppt", ylab="Median elevation of ecotone")
plot(df2$meanPET, df2$elev_range, pch=16, xlab="Annual mean PET", ylab="Elevation range of ecotone")
plot(df2$totalPr, df2$elev_range, pch=16, xlab="Annual total ppt", ylab="Elevation range of ecotone")

dev.off()




## Is topography related to ecotone --------------------------

library(boot)
twi.raw <- raster("../Data/guan_TWI.tif")

twi <- resample(twi.raw, ecotoneAllYrs)


dem2 <- mask(dem, guan)
dem2 <- resample(dem2, ecotoneAllYrs)


png("../Outputs/figSx_topography.png", res=300, units="in", width=6.5, height=2.8)

par(mfrow=c(1,2), mar=c(0.6,1.1,1.1,1.1))

plot(dem2, xaxt="n", yaxt="n")
mtext("Elevation", line=0.1)
pp <- par("usr")
text(pp[1]+0.05*diff(pp[1:2]), pp[4]-0.05*diff(pp[3:4]), "a)")
contour(dem2, add=T, levels=c(200,300,400,500,600,800,1000,1200,1400))


plot(twi, col=viridis(64), xaxt="n", yaxt="n")
mtext("Topographic Wetness Index", line=0.1)
pp <- par("usr")
text(pp[1]+0.05*diff(pp[1:2]), pp[4]-0.05*diff(pp[3:4]), "b)")

dev.off()


df.twi <- data.frame(ecotone_nyears = values(ecotone.nyears),
                     ecotoneAllYrs = values(ecotoneAllYrs),
                     twi = values(twi),
                     elev = values(dem2))

df.twi$ecotoneAllYrs[is.na(df.twi$ecotoneAllYrs)] <- 0
df.twi <- df.twi[complete.cases(df.twi),]


## logistic regression with whole study area

logreg <- glm(ecotoneAllYrs ~ twi, data=df.twi, family = "binomial")
summary(logreg)

xseq <- seq(min(df.twi$twi), max(df.twi$twi), by=0.1)
ypred <- inv.logit(predict(logreg, newdata=data.frame(twi=xseq)))

plot(df.twi$twi, df.twi$ecotoneAllYrs)
lines(xseq, ypred)


#logistic regression with restricted elevation range

etrange <- range(df.twi$elev[df.twi$ecotoneAllYrs==1])

df.twi2 <- df.twi[df.twi$elev >= min(etrange) & df.twi$elev <= max(etrange),]

logreg2 <- glm(ecotoneAllYrs ~ twi, data=df.twi2, family = "binomial")
summary(logreg2)

ypred2 <- inv.logit(predict(logreg2, newdata=data.frame(twi=xseq)))

plot(df.twi2$twi, df.twi2$ecotoneAllYrs)
lines(xseq, ypred) #very similar results to with full elevation range


#logistic GAM
library(mgcv)

loggam <- gam(ecotoneAllYrs ~ s(twi), data=df.twi, family="binomial")
summary(loggam)
gam.check(loggam)

plot(loggam)

loggam2 <- gam(ecotoneAllYrs ~ s(twi) + s(elev), data=df.twi, family="binomial")
summary(loggam2)
gam.check(loggam2)

concurvity(loggam2)


png("../Outputs/fig_topographyEffects.png", width=6.5, height=3.25, units="in", res=300)
par(mfrow=c(1,2), mar=c(3.1,3.1,1.1,1.1), mgp=c(2, 0.8, 0))
plot(loggam2, select=1, ylim=c(-5,5), xlab="Topographic Wetness Index", ylab="Effect on ecotone probability", rug=FALSE)
qq <- par("usr")
text(qq[1]+0.05*(diff(qq[1:2])), qq[4]-0.05*diff(qq[3:4]),"a)")
plot(loggam2, select=2, ylim=c(-60,20), xlim=c(100,1300), xlab="Elevation", ylab="Effect on ecotone probability", rug=FALSE)
qq <- par("usr")
text(qq[1]+0.05*(diff(qq[1:2])), qq[4]-0.05*diff(qq[3:4]),"b)")
dev.off()
