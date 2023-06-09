rm(list = ls(all=TRUE))
require(lmodel2)


file = read.csv (file = "ArmorMorphometrics.csv", header = TRUE, sep = ",") 

attach(file)
head(file)
tail(file)

logSL <- log(length)
logGeoL <- log(geoL)
logGeoW <- log(geoW)
logGeoH <- log(geoH)

regGeoL <- lmodel2(formula = logGeoL ~ logSL, data = file, nperm = 99)
regGeoL
regGeoW <- lmodel2(formula = logGeoW ~ logSL, data = file, nperm = 99)
regGeoW
regGeoH <- lmodel2(formula = logGeoH ~ logSL, data = file, nperm = 99)
regGeoH

### BMD (bone mineral density) Code ###
file = read.csv (file = "PoacherBMD_R.csv", header = TRUE, sep = ",") 
#note that 160 mm SL individual may be an outlier, remove this if needed - does not effect overall pattern

attach(file)
head(file)
tail(file)

regBMD <- lmodel2(formula = logBMD ~ logSL, data = file, nperm = 99)
regBMD
plot(regBMD)

mod <- aov(data = file, logBMD ~ logSL * scaleNum)
summary(mod)
mod <- aov(data = file, logBMD ~ logSL + scaleNum)
summary(mod)
mod <- aov(data = file, logBMD ~ logSL * bin)
summary(mod)

TukeyHSD(mod, "bin", ordered = TRUE, conf.level = 0.95)
plot(TukeyHSD(mod, "bin"))

