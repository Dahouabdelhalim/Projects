# Analyses for Van Buskirk & Smith. 2021. Ecological causes of fluctuating
# natural selection on habitat choice in an amphibian. Evolution.
#
# To install package "brms" you must first install Rtools from
# this site: https://cran.r-project.org/bin/windows/Rtools/.
# Probably best to put it in the c:\\ directory.
# After that, it should work to install.packages("brms").
# Brms uses Rcpp. If you're having trouble, it may be because
# you have an older version of Rcpp. Just remove Rcpp manually
# and reinstall it fresh. Then brms installs properly.
#
library(glmmTMB)
library(lme4)
library(brms)
library(data.table)
library(jpeg)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))    # this tells brms where to look for Rtools.
#
# Set your working directory.
setwd("C:/a_josh/Isleroyale/Manuscripts/Fluctuating selection MS/dryad package")
data.wd   <- getwd()
figs.wd   <- paste0(getwd(), "/figures")
models.wd <- paste0(getwd(), "/models")
#
# Load useful functions.
source(file="VanBuskirk_Smith_2021 sourced functions.R")
#
# Set colors. These can be explored at https://html-color-codes.info/old/colorpicker.html
#
main.color      <- "#FEA100"                                  # an orange color
main.color.pale <- adjustcolor(main.color, alpha.f = 0.08)
#
# This color scheme runs from purple up through orange
#
lower.color     <- "#1E0047"                                  # a dark purple color
lowmid.color    <- "#EFE4FF"
highmid.color   <- "#FFD386"
upper.color     <- main.color
#
# This color scheme runs from blue through orange
#
lower.color     <- "#04006E"                                  # a dark bluish color
lowmid.color    <- "#C8C7FF"
highmid.color   <- "#FFD386"
upper.color     <- main.color
#






#############################################
#
#
# Import the list of pools.
pools <- read.table("Pool list.txt", header = TRUE)
pools$max.area <- ifelse(pools$pool == 58.1, 0.16, pools$max.area)
pools          <- pools[, c("pool", "location", "max.area")]
#
#
##############################################
#
# Import the dragonfly data.
aeshna        <- read.table("Aeshna numbers.txt", header = FALSE, na.strings = ".")
names(aeshna) <- c("island", "pool", "year", "juncea23", "juncea3", "umbrosa")  # juncea23 is nr age 2 and 3 years old, juncea3 is number age 3, umbrosa is nr of 2-yr-old Aeshna umbrosa.
aeshna        <- aeshna[aeshna$island == "ng", ]
aeshna$aeshna <- ifelse(!is.na(aeshna$juncea3), aeshna$juncea23 + aeshna$umbrosa, NA)
aeshna        <- aeshna[ , c("pool", "year", "aeshna") ]
#
# Total Aeshna pop size each year.
total.aeshna          <- aggregate(aeshna[ , "aeshna"], by = list(aeshna$year), sum, na.rm = TRUE)
names(total.aeshna)   <- c("year", "total.aeshna")
aeshna                <- merge(pools, aeshna, by="pool")
aeshna$aeshna.density <- aeshna$aeshna / aeshna$max.area
#
# Mean Aeshna density in each pool over all years.
pool.list               <- data.frame(year = rep(1986:1998, each = nrow(pools)), pools[rep(1:nrow(pools), times = 13), ] )
new.aeshna              <- merge(pool.list, aeshna[,c("pool", "year", "aeshna.density")], by = c("pool", "year"), all=TRUE)
new.aeshna$aeshna.density <- ifelse(is.na(new.aeshna$aeshna.density) & new.aeshna$year>1985, 0, new.aeshna$aeshna.density)
pool.mean.aeshna        <- aggregate(new.aeshna[ , "aeshna.density"], by = list(new.aeshna$pool, new.aeshna$location, new.aeshna$max.area), mean, na.rm=TRUE)
names(pool.mean.aeshna) <- c("pool", "location", "area", "aeshna.density")
#
#
######################################
#
# Import and prepare the wave height data.
#   Keep dates during the main tadpole period (26 May through 15 July).
#   Calculate the mean wave height during that period.
#   Values are expressed as a fraction of the total height of the shore.
waves      <- read.table("Wave heights 1983-1998.txt", header = TRUE)
date.range <- range(waves$day)
dates      <- seq(date.range[1], date.range[2], 2)
new.w      <- waves[1:length(dates),]
new.w$day  <- dates
index      <- seq(1, date.range[2]-11, length.out = length(dates))
for (i in 1:length(dates)) {
  for (j in 2:dim(waves)[2]) {
    new.w[i,j] <- mean(c(waves[index[i],j], waves[index[i]+1,j]), na.rm = TRUE)
    }
  }
new.w$mean.ht <- 36 - apply(new.w[,-1], 1, mean, na.rm = TRUE)
new.w$N       <- apply(new.w[,-1], 1, samp.size)[1,] - 1
new.w         <- new.w[-1,-c(2:17)]
#
waves <- subset(waves, day > 25 & day < 75)
waves <- data.frame(year = 1983:1998, mean.wave.ht = (36 - apply(waves[,-1], 2, mean, na.rm = TRUE))/ 36)
rownames(waves) <- NULL
#
#
######################################
#
# Import and prepare the pool area/volume data.
#   Keep dates during the main tadpole period (26 May through 15 July).
#
drying <- read.table("Pool sizes 1983-1998.txt", header = TRUE)
drying <- subset(drying, date > 25 & date < 75)
dry    <- data.frame(year = 1983:1998, prop.dry = apply(drying[,34:49], 2, mean, na.rm = TRUE) )
rownames(dry) <- NULL
#
#
#########################################
#
# Import pool temperature data from Edwards pool 12.
#
t1        <- read.table("Pool temperatures 1983-1998.txt", header = FALSE, stringsAsFactors = FALSE)
names(t1) <- c("date", paste(rep(c("min", "max"), 8), rep(1983:1998, each = 2), sep = ""))
t2        <- data.frame(date = t1$date, t.1983 = NA, t.1984 = NA, t.1985 = NA, t.1986 = NA, t.1987 = NA, t.1988 = NA, t.1989 = NA, t.1990 = NA, t.1991 = NA, t.1992 = NA, t.1993 = NA, t.1994 = NA, t.1995 = NA, t.1996 = NA, t.1997 = NA, t.1998 = NA)
count     <- 1
for(i in 1983:1998) {
  count <- count+1
  t2[1,count] <- mean(c(t1[1,((count-1)*2)], t1[1,1+((count-1)*2)]), na.rm = TRUE)
  t2[2,count] <- mean(c(t1[1:2,((count-1)*2)], t1[1:2,1+((count-1)*2)]), na.rm = TRUE)
  for(j in 3:length(t1$date)) {
    t2[j, count] <- mean(c(t1[(j-2):j,((count-1)*2)], t1[(j-2):j,1+((count-1)*2)]), na.rm = TRUE)
    }
  }
t3    <- t2    # data.frame t3 will have three-day moving average temperature centered on the stated date.
count <- 1
for(i in 1983:1998) {
  count <- count+1
  t3[1, count] <- mean(t2[1:2,count], na.rm = TRUE)
  t3[length(t1$date), count] <- mean(t2[(length(t1$date)-1):length(t1$date),count], na.rm = TRUE)
  for(j in 2:(length(t1$date)-1)) {
    t3[j,count] <- mean(t2[(j-1):(j+1), count], na.rm = TRUE)
    }
  }
t4         <- data.frame(date = t3$date, N = apply(t3[,-1], 1, samp.size)[1,], temp = apply(t3[,-1], 1, mean, na.rm = TRUE), SD.temp = apply(t3[,-1], 1, sd, na.rm = TRUE))     # data.frame t4 has 3-day moving averages, averaged across years, with SD across years.
t4$SE.temp <- t4$SD.temp / sqrt(t4$N)
t4         <- NA.to.zero(t4)
#
#
# Import the tadpole brood data
#
broods        <- read.table("All tadpole broods 1983-1998.txt", header = FALSE, na.strings = ".")
colnames(broods) <- c("island", "species", "year", "pool", "brood", "catches", "origin", "date1", "inter1", "fate1", "nbrd11", "nbrd12", "npool1", "area1", "den1", "geoden1", "temp1", "size11", "size12", "date2", "inter2", "fate2", "nbrd21", "nbrd22", "npool2", "area2", "den2", "geoden2", "temp2", "size21", "size22", "date3", "inter3", "fate3", "nbrd31", "nbrd32", "npool3", "area3", "den3", "geoden3", "temp3", "size31", "size32", "date4", "inter4", "fate4", "nbrd41", "nbrd42", "npool4", "area4", "den4", "geoden4", "temp4", "size41", "size42")
broods$area1  <- ifelse(broods$area1 == -1, NA, broods$area1)
broods$area2  <- ifelse(broods$area2 == -1, NA, broods$area2)
broods$area3  <- ifelse(broods$area3 == -1, NA, broods$area3)
broods$area4  <- ifelse(broods$area4 == -1, NA, broods$area4)
broods$npool1 <- ifelse(broods$npool1 == -1, NA, broods$npool1)
broods$npool2 <- ifelse(broods$npool2 == -1, NA, broods$npool2)
broods$npool3 <- ifelse(broods$npool3 == -1, NA, broods$npool3)
broods$npool4 <- ifelse(broods$npool4 == -1, NA, broods$npool4)
broods$size11 <- ifelse(broods$size11 == -1, NA, broods$size11)
broods$size21 <- ifelse(broods$size21 == -1, NA, broods$size21)
broods$size31 <- ifelse(broods$size31 == -1, NA, broods$size31)
broods$size41 <- ifelse(broods$size41 == -1, NA, broods$size41)
broods$size12 <- ifelse(broods$size12 == -1, NA, broods$size12)
broods$size22 <- ifelse(broods$size22 == -1, NA, broods$size22)
broods$size32 <- ifelse(broods$size32 == -1, NA, broods$size32)
broods$size42 <- ifelse(broods$size42 == -1, NA, broods$size42)
broods$temp1  <- ifelse(broods$temp1 == -1, NA, broods$temp1)
broods$temp2  <- ifelse(broods$temp2 == -1, NA, broods$temp2)
broods$temp3  <- ifelse(broods$temp3 == -1, NA, broods$temp3)
broods$temp4  <- ifelse(broods$temp4 == -1, NA, broods$temp4)
broods        <- broods[broods$island=="ng" & broods$species=="pt", ]   # get rid of Pseudacris crucifer and Edwards data.
#
# Isolate first capture intervals.
#
interval.1 <- broods[ , c("year", "pool", "brood", "origin", "date1", "inter1", "fate1", "nbrd11", "nbrd12", "npool1", "area1", "den1", "temp1", "size11", "size12", "catches") ]
colnames(interval.1) <- c("year", "pool", "brood", "origin", "first.day", "elapsed", "fate", "nbrood1", "nbrood2", "npool1", "init.area", "ave.dens", "mean.temp", "size1", "size2", "catches")
interval.1$interval <- 1
interval.1$interval <- ifelse(interval.1$catches == 1, 0, interval.1$interval)    # Those that were caught only once are listed here as interval=0
interval.1$catches  <- NULL
#
# Second intervals.
interval.2 <- broods[ , c("year", "pool", "brood", "origin", "date2", "inter2", "fate2", "nbrd21", "nbrd22", "npool2", "area2", "den2", "temp2", "size21", "size22") ]
colnames(interval.2) <- c("year", "pool", "brood", "origin", "first.day", "elapsed", "fate", "nbrood1", "nbrood2", "npool1", "init.area", "ave.dens", "mean.temp", "size1", "size2")
interval.2$interval <- 2
keep <- ifelse(is.na(interval.2$first.day), FALSE, TRUE)
interval.2 <- interval.2[ keep, ]
#
# Third intervals.
interval.3 <- broods[ , c("year", "pool", "brood", "origin", "date3", "inter3", "fate3", "nbrd31", "nbrd32", "npool3", "area3", "den3", "temp3", "size31", "size32") ]
colnames(interval.3) <- c("year", "pool", "brood", "origin", "first.day", "elapsed", "fate", "nbrood1", "nbrood2", "npool1", "init.area", "ave.dens", "mean.temp", "size1", "size2")
interval.3$interval <- 3
keep <- ifelse(is.na(interval.3$first.day), FALSE, TRUE)
interval.3 <- interval.3[ keep, ]
#
# Fourth intervals.
interval.4 <- broods[ , c("year", "pool", "brood", "origin", "date4", "inter4", "fate4", "nbrd41", "nbrd42", "npool4", "area4", "den4", "temp4", "size41", "size42") ]
colnames(interval.4) <- c("year", "pool", "brood", "origin", "first.day", "elapsed", "fate", "nbrood1", "nbrood2", "npool1", "init.area", "ave.dens", "mean.temp", "size1", "size2")
interval.4$interval <- 4
keep <- ifelse(is.na(interval.4$first.day), FALSE, TRUE)
interval.4 <- interval.4[ keep, ]
#
# Bring all intervals together, calculate survival and growth rates, discard broods that washed in.
#
all.intervals <- rbind(interval.1, interval.2, interval.3, interval.4)
all.intervals$rainwash <- ifelse(all.intervals$fate == 3, "y", "n")
all.intervals$wavewash <- ifelse(all.intervals$fate == 4, "y", "n")
all.intervals$dryup    <- ifelse(all.intervals$fate == 5, "y", "n")
all.intervals$brood.dens <- all.intervals$nbrood1 / all.intervals$init.area
all.intervals$survival <- ifelse( all.intervals$elapsed > 0 & all.intervals$nbrood2 > 0, exp(log(all.intervals$nbrood2/all.intervals$nbrood1)/all.intervals$elapsed), 0 )
# If none survive at end of interval, then assume that on day (elapsed-1) there were (0.5*N/elapsed) tadpoles remaining
all.intervals$survival <- ifelse( all.intervals$elapsed > 0 & all.intervals$nbrood2 == 0, exp(log((all.intervals$nbrood1/(2*all.intervals$elapsed))/all.intervals$nbrood1)/(all.intervals$elapsed-1)), all.intervals$survival )
all.intervals$survival <- ifelse( all.intervals$interval == 0, NA, all.intervals$survival )
all.intervals$growth   <- exp(log(all.intervals$size2/all.intervals$size1)/all.intervals$elapsed)
all.intervals <- all.intervals[all.intervals$origin != 2, ]
# Estimate the day on which tadpoles were 3.5 mm (taken as their birthdate).
# The following growth model (growth = initdens meantemp) using data from 452 broods less than 4 mm size, gives a
# model with R2=0.0663: EstimatedGrowthRate = 0.0359 - 0.00001467 * initdens + 0.00085304 * meantemp.
all.intervals$est.growth  <- ifelse(all.intervals$interval == 1 & !is.na(all.intervals$growth), all.intervals$growth - 1, 0.0359-0.00001467*(all.intervals$npool1)/(all.intervals$init.area)+0.00085304*all.intervals$mean.temp)
all.intervals$est.day.3.5 <- ifelse(all.intervals$interval < 2, all.intervals$first.day - (log(all.intervals$size1/3.5)/all.intervals$est.growth), NA)
# Some estimates of date-at-3.5mm are badly off for larger tads. Therefore, we keep estimates only small tadpoles.
# This plot >> plot(all.intervals$size1, all.intervals$est.day.3.5) << showed that they should be good up to size1 = 5 or 6.
all.intervals$est.day.3.5 <- ifelse(all.intervals$size1 < 5, all.intervals$est.day.3.5, NA)
all.intervals$est.day.3.5 <- ifelse(all.intervals$est.day.3.5 < 5 & all.intervals$first.day == 47, 34.5, all.intervals$est.day.3.5)   # two odd estimates to be corrected.
all.intervals$est.day.3.5 <- ifelse(all.intervals$est.day.3.5 > 110 & all.intervals$first.day == 48, 61.5, all.intervals$est.day.3.5)
#
# Pull out the first captures (intervals 0 and 1), connect first.day and growth to metamorphic emergence and size, and calculate fitness.
first.captures <- subset(all.intervals, all.intervals$interval < 2)
first.captures$StandFirstDay <- 62.9 + 0.75145*first.captures$est.day.3.5 - 186.49*first.captures$est.growth   # These two relationships come from analysis of 1983 and 1989 broods
first.captures$StandGrowth   <- 10.338 + 28.84*first.captures$est.growth
first.captures$StandFirstDay <- (first.captures$StandFirstDay - mean(first.captures$StandFirstDay, na.rm=TRUE)) / sd(first.captures$StandFirstDay, na.rm=TRUE)
first.captures$StandGrowth   <- (first.captures$StandGrowth - mean(first.captures$StandGrowth, na.rm=TRUE)) / sd(first.captures$StandGrowth, na.rm=TRUE)
logitP <- -1.969 + 0.6013*first.captures$StandGrowth - 0.1978*first.captures$StandFirstDay    # using Dave Smith's (1987, Ecology) data to get new fitness
first.captures$DcsFitness <- exp(logitP) / (1 + exp(logitP))                                  # this is estimated survival to reproductive maturity at 2 years of age, assuming metamorphic date is proportional to hatching date and metamorphic size is proportional to growth.
first.captures$surv       <- first.captures$nbrood2 / first.captures$nbrood1
first.captures$surv       <- ifelse(first.captures$surv > 1, 1, first.captures$surv)
first.captures$NewFitness <- first.captures$surv * first.captures$DcsFitness
first.captures$surv       <- first.captures$DcsFitness <- NULL
#
# Tally up the number of broods followed for each growth interval.
interval.means <- cbind(aggregate(all.intervals[ , c("survival") ], by=list(all.intervals$interval), length), aggregate(all.intervals[ , c("survival", "growth", "nbrood1", "size1") ], by=list(all.intervals$interval), mean, na.rm=TRUE)[,2:5])
names(interval.means)[1:2] <- c("interval", "N.broods")
interval.means
#
# Count up proportion of broods that dry and wave-wash for each pool. Then combine with the Aeshna density data.
all.fates           <- expand.grid(unique(pool.mean.aeshna$pool), unique(all.intervals$year), stringsAsFactors = FALSE)
names(all.fates)    <- c("pool", "year")
nr.of.broods        <- as.data.frame(table(all.intervals$pool, all.intervals$year), stringsAsFactors = FALSE)
names(nr.of.broods) <- c("pool", "year", "broods")
all.fates           <- merge(nr.of.broods, all.fates, by = c("pool", "year"), all = TRUE)
nr.dried            <- as.data.frame(table(all.intervals[all.intervals$fate == 5, "pool"], all.intervals[all.intervals$fate == 5, "year"]))
nr.wave.washed      <- as.data.frame(table(all.intervals[all.intervals$fate == 4, "pool"], all.intervals[all.intervals$fate == 4, "year"]))
names(nr.dried)       <- c("pool", "year", "dried")
names(nr.wave.washed) <- c("pool", "year", "waves")
all.fates           <- merge(all.fates, nr.dried, by = c("pool", "year"), all = TRUE)
all.fates           <- merge(all.fates, nr.wave.washed, by = c("pool", "year"), all = TRUE)
brood.present       <- subset(all.fates, all.fates$broods > 0)
brood.present$dried <- ifelse(is.na(brood.present$dried), 0, brood.present$dried)
brood.present$waves <- ifelse(is.na(brood.present$waves), 0, brood.present$waves)
brood.present$dried <- brood.present$dried / brood.present$broods
brood.present$waves <- brood.present$waves / brood.present$broods
brood.present       <- aggregate(brood.present[,c("broods", "dried","waves")], by = list(brood.present$pool), mean, na.rm = TRUE)
all.pools           <- aggregate(all.fates[,c("broods", "dried","waves")], by = list(all.fates$pool), mean, na.rm = TRUE)
names(all.pools)[1] <- names(brood.present)[1] <- "pool"
all.pools           <- merge(pools, all.pools, by = "pool", all = TRUE)
all.pools           <- merge(all.pools, pool.mean.aeshna[,c("pool", "aeshna.density")], by = "pool", all = TRUE)
all.pools           <- subset(all.pools, pool < 400)[ , c(1,2,3,7,4,5)]
brood.present       <- merge(pools, brood.present, by = "pool", all = TRUE)
brood.present       <- subset(brood.present, broods > 0)
#
#
# Set aside the first growth interval. Center and standardize trait values.
#
fc <- merge(pools, first.captures, by="pool")
    # prepare covariate location
fc$loc          <- fc$location
mean.loc.kept   <- mean(fc$loc, na.rm = TRUE)
sd.loc.kept     <- sd(fc$loc, na.rm = TRUE)
fc$loc.s        <- (fc$loc - mean.loc.kept) / sd.loc.kept
fc$loc2         <- fc$loc.s^2
mean.loc2.kept  <- mean(fc$loc2, na.rm = TRUE)
sd.loc2.kept    <- sd(fc$loc2, na.rm = TRUE)
fc$loc2.s       <- (fc$loc2 - mean.loc2.kept) / sd.loc2.kept
    # prepare covariate area
fc$area         <- log(fc$max.area)
mean.area.kept  <- mean(fc$area, na.rm = TRUE)
sd.area.kept    <- sd(fc$area, na.rm = TRUE)
fc$area.s       <- (fc$area - mean.area.kept) / sd.area.kept
fc$area2        <- fc$area.s^2
mean.area2.kept <- mean(fc$area2, na.rm = TRUE)
sd.area2.kept   <- sd(fc$area2, na.rm = TRUE)
fc$area2.s      <- (fc$area2 - mean.area2.kept) / sd.area2.kept
    # prepare covariate date
fc$date         <- fc$est.day.3.5
mean.date.kept  <- mean(fc$date, na.rm = TRUE)
sd.date.kept    <- sd(fc$date, na.rm = TRUE)
fc$date.s       <- (fc$date - mean.date.kept) / sd.date.kept
fc$date2        <- fc$date.s^2
mean.date2.kept <- mean(fc$date2, na.rm = TRUE)
sd.date2.kept   <- sd(fc$date2, na.rm = TRUE)
fc$date2.s      <- (fc$date2 - mean.date2.kept) / sd.date2.kept
#
fc$rel.fitness <- fc$NewFitness / mean(fc$NewFitness, na.rm=TRUE)
fc$round.fitness <- round( fc$rel.fitness * 10 )
#
# Calculate average tadpole density by year.
#
fc$pooldens  <- fc$npool1 / fc$init.area
temp  <- fc[,c("year", "pooldens", "nbrood1")]
temp  <- temp[complete.cases(temp),]
temp1 <- fc[,c("pooldens", "year", "nbrood1")]
temp1$wtdensity <- temp1$pooldens * temp1$nbrood1
temp1 <- aggregate(wtdensity ~ year, data = temp1, sum)
tempw <- aggregate(nbrood1 ~ year, data = fc, sum)
dens4 <- cbind(aggregate(fc$pooldens, by = list(year = fc$year), median, na.rm = TRUE), aggregate(fc$pooldens, by = list(year = fc$year), mean, na.rm = TRUE)[,2], temp1$wtdensity / tempw$nbrood1)
names(dens4)       <- c("year", "meddens", "meandens", "wtmeandens")
dens4$l.meandens   <- log(dens4$meandens)
dens4$l.wtmeandens <- log(dens4$wtmeandens)
dens4              <- merge(dens4, aggregate(fc$nbrood1, by = list(year = fc$year), sum), by = "year")
dens4              <- rename(dens4, "x", "Nr.tads")
dens4              <- cbind(dens4, Nr.broods = apply(table(fc$year, fc$pool), 1, sum))
temp1 <- fc[,c("rel.fitness", "year", "nbrood1")]
temp1 <- temp1[complete.cases(temp1),]
temp1$wtrelfit <- temp1$rel.fitness * temp1$nbrood1
tempw <- aggregate(nbrood1 ~ year, data = temp1, sum)
temp1 <- aggregate(wtrelfit ~ year, data = temp1, sum)
dens4 <- merge(dens4, data.frame(year = temp1$year, aggregate(rel.fitness ~ year, data = fc, mean), wtmeanfit = temp1$wtrelfit / tempw$nbrood1), by = "year", all = TRUE)
dens4$year.1 <- NULL
dens4 <- rename(dens4, "rel.fitness", "meanfit")
cor(dens4, use = "pairwise.complete.obs")    # This shows that mean density and mean weighted density are tightly correlated. Median density also positive, and total nr of tadpoles.
c.dens4 <- matrix(as.numeric(format.p.val(cor(dens4, use = "pairwise.complete.obs"), 3)), ncol = 10)
rownames(c.dens4) <- colnames(c.dens4) <- names(dens4)
#
#
# Figure S1 -- photographs of study area, frog, and tadpole.
#
setwd(data.wd)
panelA <- readJPEG(source = "Fig S1 panel A.jpg")
panelB <- readJPEG(source = "Fig S1 panel B.jpg")
panelC <- readJPEG(source = "Fig S1 panel C.jpg")
upper.marg <- 0.05
lower.marg <- 0.05
margin     <- 0
middle     <- (1-(upper.marg+lower.marg))
ht1        <- 812   # heights of the three photos
ht2        <- 791
ht3        <- 740
wd1        <- 1200
middle.pixels <- (margin*2 + ht1 + ht2 + ht3)
vert.marg  <- c(1-upper.marg, 1 - upper.marg - (middle * (ht1/middle.pixels)), 1 - upper.marg - (middle * (ht1/middle.pixels)) - (middle * (margin/middle.pixels)), 1 - upper.marg - (middle * (ht1/middle.pixels)) - (middle * (ht2/middle.pixels)) - (middle * (margin/middle.pixels)), 1 - upper.marg - (middle * (ht1/middle.pixels)) - (middle * (ht2/middle.pixels)) - (middle * (2*margin/middle.pixels)), 1 - upper.marg - (middle * (ht1/middle.pixels)) - (middle * (ht2/middle.pixels)) - (middle * (ht3/middle.pixels)) - (middle * (2*margin/middle.pixels)) )
PDF.size.ratio <- (middle.pixels/(1-(upper.marg+lower.marg))) / (wd1/0.8)
setwd(figs.wd)
pdf("Fig S1.pdf", width = 7/PDF.size.ratio, height = 7, useDingbats = FALSE)
plot.new()
  par(new = "TRUE", xpd = TRUE, plt = c(0.1, 0.9, vert.marg[2], vert.marg[1]), cex.axis = 1, las = 1 )
  plot(NULL, bg = "transparent", type = "n", bty = "n", xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", xlim = c(1,2), ylim = c(1,2) )
  rasterImage(panelA, xleft = 1, ybottom = 1, xright = 2, ytop = 2, interpolate = FALSE)
  mtext(paste("Van Buskirk & Smith, Figure S1,", Sys.Date()), side = 3, line = 1, cex = 0.5, at = 1.5)
  text("A", x = 1.96, y = 1.94, cex = 1.2, col = "white")
  box(lwd = 1.3)
#
  par(new = "TRUE", xpd = TRUE, plt = c(0.1, 0.9, vert.marg[4], vert.marg[3]), cex.axis = 1, las = 1 )
  plot(NULL, bg = "transparent", type = "n", bty = "n", xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", xlim = c(1,2), ylim = c(1,2) )
  rasterImage(panelB, xleft = 1, ybottom = 1, xright = 2, ytop = 2, interpolate = FALSE)
  text("B", x = 1.96, y = 1.94, cex = 1.2, col = "white")
  box(lwd = 1.3)
#
  par(new = "TRUE", xpd = TRUE, plt = c(0.1, 0.9, vert.marg[6], vert.marg[5]), cex.axis = 1, las = 1 )
  plot(NULL, bg = "transparent", type = "n", bty = "n", xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", xlim = c(1,2), ylim = c(1,2) )
  rasterImage(panelC, xleft = 1, ybottom = 1, xright = 2, ytop = 2, interpolate = FALSE)
  text("C", x = 1.96, y = 1.94, cex = 1.2, col = "white")
  box(lwd = 1.3)
dev.off()
#
#
# Figure 1 -- maps and a study area photograph.
#
setwd(data.wd)
panelA <- readJPEG(source = "Fig 1 panel A.jpg")
panelB <- readJPEG(source = "Fig 1 panel B.jpg")
panelC <- readJPEG(source = "Fig 1 panel C.jpg")
upper.marg    <- 0.05
lower.marg    <- 0.05
side.marg     <- 0.1
margin        <- 0
middle        <- (1-(upper.marg+lower.marg))
dimA          <- dim(panelA)[1:2]  # heights of the three photos
dimB          <- dim(panelB)[1:2]
dimC          <- dim(panelC)[1:2]
scaled.width  <- dimA[2] + dimB[2]
scaled.height <- dimA[1] + ((scaled.width / dimC[2]) * 813)
hori.marg     <- c(0.1, 0.1 + 0.8 * dimA[2]/(dimA[2]+dimB[2]), 0.9)
vert.marg     <- c(1-upper.marg, (1 - upper.marg) - (dimA[1]/scaled.height)*middle, lower.marg)
line.width    <- 2
setwd(figs.wd)
pdf("Fig 1.pdf", width = 10 * (scaled.width/(1-2*side.marg))/(scaled.height/(1-upper.marg+lower.marg)), height = 10, useDingbats = FALSE)
plot.new()
  par(new = "TRUE", xpd = TRUE, plt = c(hori.marg[1], hori.marg[2], vert.marg[2], vert.marg[1]), cex.axis = 1, las = 1 )
  plot(NULL, bg = "transparent", type = "n", bty = "n", xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", xlim = c(1,2), ylim = c(1,2) )
  rasterImage(panelA, xleft = 0.97, ybottom = 0.97, xright = 2.03, ytop = 2.03, interpolate = FALSE)
  mtext(paste("Van Buskirk & Smith, Figure 1,", Sys.Date()), side = 3, line = 1, cex = 0.8, at = 1.5)
  text("A", x = 1.97, y = 1.97, cex = 1.6, col = "black")
  box(lwd = line.width)
#
  par(new = "TRUE", xpd = TRUE, plt = c(hori.marg[2], hori.marg[3], vert.marg[2], vert.marg[1]), cex.axis = 1, las = 1 )
  plot(NULL, bg = "transparent", type = "n", bty = "n", xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", xlim = c(1,2), ylim = c(1,2) )
  rasterImage(panelB, xleft = 0.967, ybottom = 0.967, xright = 2.033, ytop = 2.033, interpolate = FALSE)
  text("B", x = 1.97, y = 1.97, cex = 1.6, col = "black")
  box(lwd = line.width)
#
  par(new = "TRUE", xpd = TRUE, plt = c(hori.marg[1], hori.marg[3], vert.marg[3], vert.marg[2]), cex.axis = 1, las = 1 )
  plot(NULL, bg = "transparent", type = "n", bty = "n", xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", xlim = c(1,2), ylim = c(1,2) )
  rasterImage(panelC, xleft = 0.96, ybottom = 0.96, xright = 2.04, ytop = 2.04, interpolate = FALSE)
  text("C", x = 1.997, y = 1.988, cex = 1.6, col = "black")
  box(lwd = line.width)
dev.off()
#
#
# Figure S2  -- seasonal time course of water temperature and tadpole density.
#
# Organize data on density of tadpoles through the season.
all.intervals$pool.dens      <- all.intervals$npool1 / all.intervals$init.area
first.intervals              <- subset(all.intervals, interval == 1)
first.intervals              <- first.intervals[,c("year", "pool", "brood", "first.day", "brood.dens", "pool.dens")]
first.intervals$older.dens   <- first.intervals$pool.dens - first.intervals$brood.dens
first.intervals$l.older.dens <- log(1 + first.intervals$older.dens)
first.intervals$l.pool.dens  <- log(1 + first.intervals$pool.dens)
first.intervals$l.brood.dens <- log(1 + first.intervals$brood.dens)
first.intervals$week <- 4                                                                 # week centered on 4 May
first.intervals$week <- ifelse(first.intervals$first.day > 7, 11, first.intervals$week)   # week centered on 11 May
first.intervals$week <- ifelse(first.intervals$first.day > 14, 18, first.intervals$week)   # week centered on 18 May
first.intervals$week <- ifelse(first.intervals$first.day > 21, 25, first.intervals$week)   # week centered on 25 May
first.intervals$week <- ifelse(first.intervals$first.day > 28, 32, first.intervals$week)   # week centered on 1 June
first.intervals$week <- ifelse(first.intervals$first.day > 35, 39, first.intervals$week)   # week centered on 8 June
first.intervals$week <- ifelse(first.intervals$first.day > 42, 46, first.intervals$week)   # week centered on 15 June
first.intervals$week <- ifelse(first.intervals$first.day > 49, 53, first.intervals$week)   # week centered on 22 June
first.intervals$week <- ifelse(first.intervals$first.day > 56, 60, first.intervals$week)   # week centered on 29 June
first.intervals$week <- ifelse(first.intervals$first.day > 63, 67, first.intervals$week)   # week centered on 6 July
first.intervals$week <- ifelse(first.intervals$first.day > 70, 74, first.intervals$week)   # week centered on 13 July
first.intervals$week <- ifelse(first.intervals$first.day > 77, 81, first.intervals$week)   # week centered on 20 July
first.intervals$week <- ifelse(first.intervals$first.day > 84, 88, first.intervals$week)   # all the rest: week centered on 27 July
f1    <- aggregate(first.intervals[, c("pool.dens", "brood.dens", "older.dens", "l.pool.dens", "l.brood.dens", "l.older.dens")], by = list(first.intervals$week, first.intervals$year), mean, na.rm = TRUE)
f2    <- aggregate(f1[, c("pool.dens", "brood.dens", "older.dens", "l.pool.dens", "l.brood.dens", "l.older.dens")], by = list(f1$Group.1), mean, na.rm = TRUE)
f2.sd <- aggregate(f1[, c("pool.dens", "brood.dens", "older.dens", "l.pool.dens", "l.brood.dens", "l.older.dens")], by = list(f1$Group.1), sd, na.rm = TRUE)
names(f2.sd)[2:7] <- c("pool.dens.SD", "brood.dens.SD", "older.dens.SD", "1.pool.dens.SD", "1.brood.dens.SD", "1.older.dens.SD")
f2.N  <- aggregate(f1[, c("pool.dens")], by = list(f1$Group.1), samp.size)
f3    <- merge(f2, f2.sd, by = "Group.1")
f3    <- cbind(f3, f2.N[,2][,1])
names(f3)[c(1,14)] <- c("week", "N")
f3$pool.dens.SE  <- f3$pool.dens.SD / sqrt(f3$N)
f3$brood.dens.SE <- f3$brood.dens.SD / sqrt(f3$N)
f3$older.dens.SE <- f3$older.dens.SD / sqrt(f3$N)
f4               <- f3[ , c(1,14,2,8,15,3,9,16,4,10,17)]
f4               <- NA.to.zero(f4)
sample.sizes     <- data.frame(week = f4$week, N.broods = aggregate(first.intervals[, c("older.dens")], by = list(first.intervals$week), samp.size)[,2][,1])
#
all.intervals$pool.dens  <- all.intervals$npool1 / all.intervals$init.area
all.intervals.dens       <- merge(all.intervals, pools, by = "pool", all = TRUE)
all.intervals.dens       <- all.intervals.dens[,c(1,2,24,25,26,27)]
all.intervals.dens$larea <- log10(all.intervals.dens$max.area)
all.intervals.dens$ldens <- log10(all.intervals.dens$pool.dens)
all.intervals.dens       <- merge(all.intervals.dens, dens4, by = "year", all = TRUE)
names(all.intervals.dens) <- c("year", "pool", "date", "pool.dens", "loc", "area", "larea", "ldens", "junk1", "year.dens", "junk2", "lyear.dens", "junk3", "Nr.tads", "junk4", "junk5", "junk6")
all.intervals.dens$date  <- all.intervals.dens$junk1 <- all.intervals.dens$junk2 <- all.intervals.dens$junk3 <- all.intervals.dens$junk4 <- all.intervals.dens$junk5<- all.intervals.dens$junk6 <- NULL
all.intervals.dens       <- all.intervals.dens[complete.cases(all.intervals.dens), ]
#
line.color.pale <- adjustcolor(main.color, alpha.f = 0.2)
x.limits        <- c(18, 90)
tick.label.size <- 1
axis.label.size <- 1.2
#
setwd(figs.wd)
pdf("Fig S2.pdf", useDingbats = FALSE, width = 5.5, height = 7.5)
#
# Top panel is temperature through time
plot.new()
par(new = "TRUE", plt = c(0.25, 0.85, 0.56, 0.86))
  y.limits <- c(11.5, 23)
  plot(c(15,20), c(10,20), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = x.limits, ylim = y.limits)
  polygon(x = c(t4$date[-1], rev(t4$date)[-1]), y = c(t4$temp[-1]-t4$SE.temp[-1], rev(t4$temp[-1]+t4$SE.temp[-1])), border = NA, col = line.color.pale)
  lines(t4$date[-1], t4$temp[-1], lwd = 2.5, col = main.color)
  box(lwd = 1.3)
  axis(side = 1, at = c(21, 41, 61, 81), labels = FALSE, tck = -0.025)
  x.text <- c(paste(make.dates(21), collapse = " "), paste(make.dates(41), collapse = " "), paste(make.dates(61), collapse = " "), paste(make.dates(81), collapse = " "))
  mtext(text = x.text, side = 1, line = 0.4, las = 1, at = seq(20,80,20), cex = tick.label.size)
  mtext(text = "Date", side = 1, line = 1.55, cex = axis.label.size)
  axis(side = 2, labels = FALSE, tck = -0.025)
  mtext(text = c("14","18","22"), side = 2, line = 0.4, las = 1, at = c(14,18,22), cex = tick.label.size)
  mtext( text = "Water temperature (°C)", side = 2, line = 2.1, cex = axis.label.size)
  mtext( text = "A", at = 22, side = 4, line = 0.3, las = 1, cex = 1.6)
  mtext( text = paste("Van Buskirk & Smith, Figure S2,", Sys.Date()), side = 3, line = 4, cex = 0.65)
#
# Lower panel is tadpole density through time
par(new = "TRUE", plt = c(0.25, 0.85, 0.12, 0.47))
  plot(c(15,20), c(10,20), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = x.limits, ylim = log(c(1, 65)))
  polygon(x = c(f4$week, rev(f4$week)), y = c(log(1+f4$older.dens)-log(1+f4$older.dens.SE), rev(log(1+f4$older.dens)+log(1+f4$older.dens.SE))), border = NA, col = line.color.pale)
  lines(f4$week, log(1+f4$older.dens), lwd = 2.5, col = main.color)
  axis(side = 1, at = c(21, 41, 61, 81), labels = FALSE, tck = -0.025)
  x.text <- c(paste(make.dates(21), collapse = " "), paste(make.dates(41), collapse = " "), paste(make.dates(61), collapse = " "), paste(make.dates(81), collapse = " "))
  mtext(text = x.text, side = 1, line = 0.4, las = 1, at = seq(20,80,20), cex = tick.label.size)
  mtext(text = "First capture date of brood", side = 1, line = 1.55, cex = axis.label.size)
  axis(side = 2, at = log((1+c(0, 10, 50))), labels = FALSE, tck = -0.025)
  axis(side = 2, at = log((1+c(1, 2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 60, 70, 80))), labels = FALSE, tck = -0.015)
  mtext(text = c("0","10","50"), side = 2, line = 0.4, las = 1, at = log(1+c(0,10,50)), cex = tick.label.size)
  mtext( text = "Density of older tadpoles (nr./m  )", side = 2, line = 2.1, cex = axis.label.size)
  mtext( text = "2", side = 2, at = log(74), line = 2.4, cex = 0.9)
  # sample sizes
  text(f4$N, x = f4$week, y = log(1.1), cex = 0.8)
  text(c("nr. of", "years"), x = 19.7, y = c(log(c(1.75, 1.3))), adj = 0.5, cex = 0.75)
  text(sample.sizes$N.broods, x = sample.sizes$week, y = log(57), cex = 0.8)
  text(c("nr. of", "broods"), x = 20.3, y = c(log(c(45, 35))), adj = 0.5, cex = 0.75)
  mtext( text = "B", at = log(50), side = 4, line = 0.3, las = 1, cex = 1.6)
dev.off()
#
#
# Fig. 3 -- distribution of Aeshna, wave wash, and drying.
#
no.aeshna  <- subset(all.pools, aeshna.density == 0)
yes.aeshna <- subset(all.pools, aeshna.density > 0)
no.dry     <- subset(brood.present, dried == 0)
yes.dry    <- subset(brood.present, dried > 0)
no.waves   <- subset(brood.present, waves == 0)
yes.waves  <- subset(brood.present, waves > 0)
#
x.limits        <- c(0.02,1)
y.limits        <- c(0.01, 17.12)
axis.label.size <- 1.15
tick.label.size <- 1
small.point.size<- 0.65
magnif          <- c(2.5,20,20)      # size of bubbles (for Aeshna, washed, and dried in that order).
#
setwd(figs.wd)
pdf("Fig 3.pdf", useDingbats = FALSE, width = 5.0, height = 9)
plot.new()
# Bottom panel is proportion of broods wave-washed.
par(new = "TRUE", plt = c(0.25, 0.85, 0.10, 0.35))
  plot(c(1,2), c(1,2), log =  "y", xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = x.limits, ylim = y.limits)
  points(yes.dry$location, yes.dry$max.area, cex = sqrt(1.5+magnif[3]*yes.dry$dried), col = "black", xaxt = "n", yaxt = "n", bg = main.color, pch = 21 )
  points(no.dry$location, no.dry$max.area, cex = small.point.size, pch = 19 )
  box(lwd = 1.3)
  axis(side = 1, labels = FALSE, tck = -0.025)
  mtext(text = c("0.0","0.2","0.4","0.6","0.8","1.0"), side = 1, line = 0.4, las = 1, at = c(0,0.2,0.4,0.6,0.8,1), cex = tick.label.size)
  axis(side = 2, at = c(0.01,0.1,1,10), labels = FALSE, tck = -0.025)
  axis(side = 2, at = c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20), labels = FALSE, tck = -0.011)
  mtext(text = c("0.01","0.1","1.0","10"), side = 2, line = 0.4, las = 1, at = c(0.01,0.1,1,10), cex = tick.label.size)
  mtext( text = "Pool location", side = 1, line = 1.7, at = 0.5, cex = axis.label.size)
  mtext( text = "Pool surface area (m  )", side = 2, line = 2.1, cex = axis.label.size)
  mtext( text = "2", at = 9.4, side = 2, line = 2.5, cex = 0.9)
  text(labels = "C. Dried", x = -0.02, y = 12.5, pos = 4, cex = axis.label.size)
  # Now add the legend lower left
  x.values   <- c(0.03, 0.23)
  ticks.at   <- c(0,0.5,1)
  tick.marks <- x.values[1] + (x.values[2]-x.values[1])*ticks.at
  diameters  <- sqrt(1.5+magnif[3]*ticks.at)
  lines(x = x.values, y = c(0.038, 0.038), lwd = 1.3)
  lines(x = c(tick.marks[1], tick.marks[1]), y = c(0.043, 0.036), lwd = 1.2)
  lines(x = c(tick.marks[2], tick.marks[2]), y = c(0.043, 0.036), lwd = 1.2)
  lines(x = c(tick.marks[3], tick.marks[3]), y = c(0.043, 0.036), lwd = 1.2)
  text(labels = ticks.at, x = tick.marks, y = 0.056, cex = 0.8)
  text(labels = c("Prop. broods", "dried"), x = 0.28, y = c(0.018,0.012), pos = 4, cex = 0.8)
  points(x = tick.marks, y = c(0.017, 0.017, 0.017), cex = diameters, col = "black", xaxt = "n", yaxt = "n", bg = main.color, pch = 21 )
  # Arrows at the bottom
  par(xpd = TRUE)
  mtext(text = c("Lake", "Forest"), side = 1, line = 1.4, at = c(0.1, 0.91), adj = c(0,1), cex = tick.label.size, las = 0)
  arrows(x0 = 0.09, y0 = 0.0019, x1 = -0.01, y1 = 0.0019, angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  arrows(x0 = 0.92, y0 = 0.0019, x1 = 1.02, y1 = 0.0019, angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  par(xpd = FALSE)
# Middle panel is proportion of broods dried.
par(new = "TRUE", plt = c(0.25, 0.85, 0.40, 0.65))
  plot(c(1,2), c(1,2), log =  "y", xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = x.limits, ylim = y.limits)
  points(yes.waves$location, yes.waves$max.area, cex = sqrt(1.5+magnif[2]*yes.dry$dried), col = "black", xaxt = "n", yaxt = "n", bg = main.color, pch = 21 )
  points(no.waves$location, no.waves$max.area, cex = small.point.size, pch = 19 )
  box(lwd = 1.3)
  axis(side = 1, labels = FALSE, tck = -0.025)
  mtext(text = c("0.0","0.2","0.4","0.6","0.8","1.0"), side = 1, line = 0.4, las = 1, at = c(0,0.2,0.4,0.6,0.8,1), cex = tick.label.size)
  axis(side = 2, at = c(0.01,0.1,1,10), labels = FALSE, tck = -0.025)
  axis(side = 2, at = c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20), labels = FALSE, tck = -0.011)
  mtext(text = c("0.01","0.1","1.0","10"), side = 2, line = 0.4, las = 1, at = c(0.01,0.1,1,10), cex = tick.label.size)
  mtext( text = "Pool surface area (m  )", side = 2, line = 2.1, cex = axis.label.size)
  mtext( text = "2", at = 9.4, side = 2, line = 2.5, cex = 0.9)
  text(labels = "B. Washed by waves", x = -0.02, y = 12.5, pos = 4, cex = axis.label.size)
  # Now add the legend lower left
  x.values   <- c(0.03, 0.23)
  ticks.at   <- c(0,0.5,1)
  tick.marks <- x.values[1] + (x.values[2]-x.values[1])*ticks.at
  diameters  <- sqrt(1.5+magnif[2]*ticks.at)
  lines(x = x.values, y = c(0.038, 0.038), lwd = 1.3)
  lines(x = c(tick.marks[1], tick.marks[1]), y = c(0.043, 0.036), lwd = 1.2)
  lines(x = c(tick.marks[2], tick.marks[2]), y = c(0.043, 0.036), lwd = 1.2)
  lines(x = c(tick.marks[3], tick.marks[3]), y = c(0.043, 0.036), lwd = 1.2)
  text(labels = ticks.at, x = tick.marks, y = 0.056, cex = 0.8)
  text(labels = c("Prop. broods", "washed out"), x = 0.28, y = c(0.018,0.012), pos = 4, cex = 0.8)
  points(x = tick.marks, y = c(0.017, 0.017, 0.017), cex = diameters, col = "black", xaxt = "n", yaxt = "n", bg = main.color, pch = 21 )
# Top panel is Aeshna mean density 1987-2014.
par(new = "TRUE", plt = c(0.25, 0.85, 0.70, 0.95))
  plot(no.aeshna$location, no.aeshna$area, log =  "y", xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = x.limits, ylim = y.limits)
  points(yes.aeshna$location, yes.aeshna$max.area, cex = sqrt(1.5+magnif[1]*yes.aeshna$aeshna.density), col = "black", xaxt = "n", yaxt = "n", bg = main.color, pch = 21 )
  points(no.aeshna$location, no.aeshna$max.area, cex = small.point.size, pch = 19 )
  box(lwd = 1.3)
  axis(side = 1, labels = FALSE, tck = -0.025)
  mtext(text = c("0.0","0.2","0.4","0.6","0.8","1.0"), side = 1, line = 0.4, las = 1, at = c(0,0.2,0.4,0.6,0.8,1), cex = tick.label.size)
  axis(side = 2, at = c(0.01,0.1,1,10), labels = FALSE, tck = -0.025)
  axis(side = 2, at = c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20), labels = FALSE, tck = -0.011)
  mtext(text = c("0.01","0.1","1.0","10"), side = 2, line = 0.4, las = 1, at = c(0.01,0.1,1,10), cex = tick.label.size)
  mtext( text = "Pool surface area (m  )", side = 2, line = 2.1, cex = axis.label.size)
  mtext( text = "2", at = 9.4, side = 2, line = 2.5, cex = 0.9)
  text(labels = substitute(paste("A. ", italic("Aeshna"), " density")), x = -0.02, y = 12.5, pos = 4, cex = axis.label.size)
  # Now add the legend lower left
  x.values   <- c(0.03, 0.23)
  ticks.at   <- c(0,5,10)
  tick.marks <- x.values[1] + (x.values[2]-x.values[1])*(ticks.at/ticks.at[3])
  diameters  <- sqrt(1.5+magnif[1]*ticks.at)
  lines(x = x.values, y = c(0.037, 0.037), lwd = 1.3)
  lines(x = c(tick.marks[1], tick.marks[1]), y = c(0.042, 0.035), lwd = 1.2)
  lines(x = c(tick.marks[2], tick.marks[2]), y = c(0.042, 0.035), lwd = 1.2)
  lines(x = c(tick.marks[3], tick.marks[3]), y = c(0.042, 0.035), lwd = 1.2)
  text(labels = ticks.at, x = tick.marks, y = 0.054, cex = 0.8)
  text(labels = c("Density of", "dragonflies"), x = 0.28, y = c(0.018,0.011), pos = 4, cex = 0.8)
  points(x = tick.marks, y = c(0.017, 0.017, 0.017), cex = diameters, col = "black", xaxt = "n", yaxt = "n", bg = main.color, pch = 21 )
  mtext(text = paste("Van Buskirk & Smith, Figure 3,", Sys.Date()), side = 3, line = 1.3, las = 1, at = 0.1, cex = 0.7)
dev.off()
#
#
# Sorting out sample sizes
#
nr.total     <- dim(interval.1)[1]
nr.washed.in <- dim(interval.1[interval.1$origin == 2,])[1]
int.1.loc    <- interval.1[interval.1$origin == 1,]
nr.vanished  <- dim(int.1.loc[int.1.loc$fate == 9,])[1]
int.1.no.2nd <- int.1.loc[int.1.loc$fate == 6,]
nr.no.2nd    <- dim(int.1.no.2nd)[1]
int.1.okay   <- int.1.loc[int.1.loc$fate < 6,]
nr.okay      <- dim(int.1.okay)[1]
nr.complete  <- samp.size(fc$round.fitness)[1]
nr.no.date   <- nr.okay - nr.complete
what.happened <- data.frame(fate = c("Total", "DSQ: washed in", "Vanished no reason", "No second measure", "No hatch date", "Complete"), number = c(nr.total, nr.washed.in, nr.vanished, nr.no.2nd, nr.no.date, nr.complete))
what.happened
what.happened[1,2] - sum(what.happened[-1,2])    # all broods accounted for.
#
#
# Fig. S4 -- How do conditions change from year to year?
#
conditions        <- merge(merge(merge(total.aeshna, waves, by = "year", all = TRUE), dry, by = "year", all = TRUE), dens4, by = "year", all = TRUE)
conditions        <- subset(conditions, year<1999)
conditions$l.Nr.tads <- log(conditions$Nr.tads)
d                 <- subset(fc, interval < 2)
d$pool.dens.1     <- d$npool1 / d$init.area
year.means        <- aggregate(d[,c("loc", "init.area", "est.day.3.5", "pool.dens.1", "mean.temp")], by = list(d$year), mean, na.rm = TRUE)
year.N            <- aggregate(d[,c("loc")], by = list(d$year), length)
year.means        <- merge(year.N, year.means, by = "Group.1")
names(year.means) <- c("year", "N", "loc", "area", "date", "density", "temp")
conditions        <- merge(conditions, year.means)
conditions        <- conditions[,c(1,15,8,9,14,3,4,2,16:18)]
names(conditions) <- c("year", "N", "l.meandens.year", "l.wtmeandens.year", "l.Ntads.year", "mean.waves.year", "prop.dry.year", "aeshna.year", "loc.chosen", "area.chosen", "date.chosen")
conditions$aeshna.year <- ifelse(conditions$aeshna.year == 0, NA, conditions$aeshna.year)
#
# Names ending with "year" are the average values of those causative agents during that year.
# Names ending in "chosen" are the average values of the pools chosen by frogs in that year.
# "N" is the number of broods detected in that year. This is the sample size for the "chosen" variables.
#
# Correlations among the four ecological risks.
#
conditions$l.aeshna.year <- log(conditions$aeshna.year)
cor(conditions[ , c("prop.dry.year", "mean.waves.year", "l.aeshna.year", "l.Ntads.year")], use = "pairwise.complete.obs")
critical.r(16)     # critical value of r
critical.r(13)     # critical value for correlations involving Aeshna (13 years).
cor.test(x = conditions$prop.dry.year, y = conditions$l.Ntads.year, use = "pairwise.complete.obs")
cor.test(x = conditions$prop.dry.year, y = conditions$mean.waves.year, use = "pairwise.complete.obs")
#
# Fig. S4.
#
x.to.y <- function(x, rx, ry) {
  # function to get predicted y from one or more x values.
  # x is one or more numeric values.
  # rx and ry are the ranges of the x variables and y variables.
  slope <- (ry[2]-ry[1]) / (rx[2]-rx[1])
  ry[1] + ((x - rx[1]) * slope)
  }
color.list      <- data.frame(agent = c("waves", "dry", "aeshna", "dens", "loc", "area", "date", "temp"), col = c("#240B95", "#E4AB00", "#7E1001", "black", "blue", "#006D0F", "#BE008B", "#00F235"), stringsAsFactors = FALSE)
color.list$col  <- rep("black", 8)      # just make everything black.
tick.label.size <- 0.9
point.size      <- 1.4
line.size       <- 1.2
axis.label.size <- 1.05
xlimits         <- c(1982.5, 1998.6)
left.edge       <- 0.08
right.edge      <- 0.05
buffer          <- 0.06
panel.width     <- (1 - ((3 * buffer) + (left.edge + right.edge))) / 4
x.margins       <- c(left.edge, left.edge + panel.width, left.edge + panel.width + buffer, left.edge + 2 * panel.width + buffer, left.edge + 2 * panel.width + 2 * buffer, left.edge + 3 * panel.width + 2 * buffer, left.edge + 3 * panel.width + 3 * buffer, (1 - right.edge))
#
setwd(figs.wd)
pdf("Fig S4.pdf", useDingbats = FALSE, width = 13.0, height = 6.2)
plot.new()
#
# Panel A: Location
par(new = "TRUE", plt = c(x.margins[1], x.margins[2], 0.55, 0.90))
  ylimits <- c(0.965, 1.01) * range(conditions$loc.chosen)
  ylimits[1] <- 0.46
  plot(c(1,2), c(1,2), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = xlimits, ylim = ylimits, xaxs="i", yaxs="i")
  points(x = conditions$year, y = conditions$loc.chosen, pch = 20, cex = point.size, col = color.list$col[5])
  lines(x = conditions$year, y = conditions$loc.chosen, lwd = line.size, col = color.list$col[5])
  box(lwd = 1.2)
  axis(side = 1, labels = FALSE, at = seq(1983,1998,5), tck = -0.025)
  axis(side = 1, labels = FALSE, at = c(1984,1985,1986,1987,1989,1990,1991,1992,1994,1995,1996,1997), tck = -0.018)
  mtext(text = c("1983","1988","1993","1998"), side = 1, line = 0.4, las = 1, at = seq(1983,1998,5), cex = tick.label.size)
  axis(side = 2, at = c(0.46,0.5,0.54,0.58), labels = FALSE, tck = -0.022)
  axis(side = 2, at = seq(0.47,0.61,0.01), labels = FALSE, tck = -0.018)
  mtext(text = c("0.46", "0.50", "0.54", "0.58"), side = 2, line = 0.35, las = 1, at = c(0.46,0.5,0.54,0.58), cex = tick.label.size)
  mtext( text = "Pool location", side = 2, line = 2.0, cex = axis.label.size)
  text(labels = "A", y = ylimits[1] + 0.94*diff(ylimits), x = 1983.7, cex = 1.4)
  mtext( text = paste("Van Buskirk & Smith, Figure S4,", Sys.Date()), side = 3, line = 2.1, cex = 0.6)
#
# Panel B: Pool area
par(new = "TRUE", plt = c(x.margins[3], x.margins[4], 0.55, 0.90))
  ylimits <- c(0.96, 1.03) * range(conditions$area.chosen)
  plot(c(1,2), c(1,2), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = xlimits, ylim = ylimits, xaxs="i", yaxs="i")
  points(x = conditions$year, y = conditions$area.chosen, pch = 20, cex = point.size, col = color.list$col[6])
  lines(x = conditions$year, y = conditions$area.chosen, lwd = line.size, col = color.list$col[6])
  box(lwd = 1.2)
  axis(side = 1, labels = FALSE, at = seq(1983,1998,5), tck = -0.025)
  axis(side = 1, labels = FALSE, at = c(1984,1985,1986,1987,1989,1990,1991,1992,1994,1995,1996,1997), tck = -0.018)
  mtext(text = c("1983","1988","1993","1998"), side = 1, line = 0.4, las = 1, at = seq(1983,1998,5), cex = tick.label.size)
  axis(side = 2, labels = FALSE, tck = -0.022)
  mtext(text = c("1.0", "1.2", "1.4"), side = 2, line = 0.35, las = 1, at = c(1,1.2,1.4), cex = tick.label.size)
  mtext(text = "Pool area (m  )", side = 2, line = 1.7, cex = axis.label.size)
  mtext(text = "2", side = 2, line = 1.9, at = 1.385, cex = axis.label.size-0.25)
  text(labels = "B", y = ylimits[1] + 0.94*diff(ylimits), x = 1983.7, cex = 1.4)
#
# Panel C: Date chosen
par(new = "TRUE", plt = c(x.margins[5], x.margins[6], 0.55, 0.90))
  ylimits <- c(0.97, 1.03) * range(conditions$date.chosen)
  plot(c(1,2), c(1,2), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = xlimits, ylim = ylimits, xaxs="i", yaxs="i")
  points(x = conditions$year, y = conditions$date.chosen, pch = 20, cex = point.size, col = color.list$col[7])
  lines(x = conditions$year, y = conditions$date.chosen, lwd = line.size, col = color.list$col[7])
  box(lwd = 1.2)
  axis(side = 1, labels = FALSE, at = seq(1983,1998,5), tck = -0.025)
  axis(side = 1, labels = FALSE, at = c(1984,1985,1986,1987,1989,1990,1991,1992,1994,1995,1996,1997), tck = -0.018)
  mtext(text = c("1983","1988","1993","1998"), side = 1, line = 0.4, las = 1, at = seq(1983,1998,5), cex = tick.label.size)
  axis(side = 2, at = c(41,46,51,56), labels = FALSE, tck = -0.018)
  mtext(text = c(paste(make.dates(41)[[1]], make.dates(41)[[2]], sep = " "), paste(make.dates(46)[[1]], make.dates(46)[[2]], sep = " "), paste(make.dates(51)[[1]], make.dates(51)[[2]], sep = " "), paste(make.dates(56)[[1]], make.dates(56)[[2]], sep = " ")), side = 2, line = 0.3, las = 1, at = c(41,46,51,56), cex = tick.label.size-0.1)
  mtext( text = "Hatching date", side = 2, line = 2.3, cex = axis.label.size)
  text(labels = "C", y = ylimits[1] + 0.94*diff(ylimits), x = 1983.7, cex = 1.4)
#
# Panel D: Sample size
par(new = "TRUE", plt = c(x.margins[7], x.margins[8], 0.55, 0.90))
  ylimits <- c(0.93, 1.03) * range(conditions$N)
  ylimits[1] <- 29.5
  plot(c(1,2), c(1,2), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = xlimits, ylim = ylimits, xaxs="i", yaxs="i")
  points(x = conditions$year, y = conditions$N, pch = 20, cex = point.size, col = color.list$col[7])
  lines(x = conditions$year, y = conditions$N, lwd = line.size, col = color.list$col[7])
  box(lwd = 1.2)
  axis(side = 1, labels = FALSE, at = seq(1983,1998,5), tck = -0.025)
  axis(side = 1, labels = FALSE, at = c(1984,1985,1986,1987,1989,1990,1991,1992,1994,1995,1996,1997), tck = -0.018)
  mtext(text = c("1983","1988","1993","1998"), side = 1, line = 0.4, las = 1, at = seq(1983,1998,5), cex = tick.label.size)
  axis(side = 2, at = c(30,50,70,90), labels = FALSE, tck = -0.022)
  axis(side = 2, at = seq(40,105,5), labels = FALSE, tck = -0.018)
  mtext(text = c("30","50","70","90"), side = 2, line = 0.35, las = 1, at = c(30,50,70,90), cex = tick.label.size)
  mtext(text = c("Number of broods"), side = 2, line = 1.5, cex = axis.label.size)
  text(labels = "D", y = ylimits[1] + 0.94*diff(ylimits), x = 1983.7, cex = 1.4)
#
# Panel E: Drying
par(new = "TRUE", plt = c(x.margins[1], x.margins[2], 0.13, 0.48))
  ylimits <- c(0, 1.04 * max(conditions$prop.dry.year))
  plot(c(1,2), c(1,2), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = xlimits, ylim = ylimits, xaxs="i", yaxs="i")
  points(x = conditions$year, y = conditions$prop.dry.year, pch = 20, cex = point.size, col = color.list$col[2])
  lines(x = conditions$year, y = conditions$prop.dry.year, lwd = line.size, col = color.list$col[2])
  box(lwd = 1.2)
  axis(side = 1, labels = FALSE, at = seq(1983,1998,5), tck = -0.025)
  axis(side = 1, labels = FALSE, at = c(1984,1985,1986,1987,1989,1990,1991,1992,1994,1995,1996,1997), tck = -0.018)
  mtext(text = c("1983","1988","1993","1998"), side = 1, line = 0.4, las = 1, at = seq(1983,1998,5), cex = tick.label.size)
  mtext(text = "Year", side = 1, line = 1.5, cex = axis.label.size)
  axis(side = 2, labels = FALSE, at = c(0,0.05,0.1,0.15,0.2), tck = -0.025)
  axis(side = 2, labels = FALSE, at = seq(0,0.20,0.01), tck = -0.018)
  mtext(text = c("0.00","0.05","0.10","0.15","0.20"), side = 2, line = 0.35, las = 1, at = c(0,0.05,0.1,0.15,0.2), cex = tick.label.size)
  mtext( text = "Proportion of pools drying", side = 2, line = 2.0, cex = axis.label.size)
  text(labels = "E", y = ylimits[1] + 0.94*diff(ylimits), x = 1983.7, cex = 1.4)
#
# Panel F: Wave wash
par(new = "TRUE", plt = c(x.margins[3], x.margins[4], 0.13, 0.48))
  ylimits <- c(0, 1.06 * max(conditions$mean.waves.year))
  plot(c(1,2), c(1,2), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = xlimits, ylim = ylimits, xaxs="i", yaxs="i")
  points(x = conditions$year, y = conditions$mean.waves.year, pch = 20, cex = point.size, col = color.list$col[1])
  lines(x = conditions$year, y = conditions$mean.waves.year, lwd = line.size, col = color.list$col[1])
  box(lwd = 1.2)
  axis(side = 1, labels = FALSE, at = seq(1983,1998,5), tck = -0.025)
  axis(side = 1, labels = FALSE, at = c(1984,1985,1986,1987,1989,1990,1991,1992,1994,1995,1996,1997), tck = -0.018)
  mtext(text = c("1983","1988","1993","1998"), side = 1, line = 0.4, las = 1, at = seq(1983,1998,5), cex = tick.label.size)
  mtext(text = "Year", side = 1, line = 1.5, cex = axis.label.size)
  axis(side = 2, labels = FALSE, at = c(0,0.05,0.1,0.15), tck = -0.025)
  axis(side = 2, labels = FALSE, at = seq(0,0.21,0.01), tck = -0.018)
  mtext(text = c("0.00","0.05","0.10","0.15","0.20"), side = 2, line = 0.35, las = 1, at = c(0,0.05,0.1,0.15,0.2), cex = tick.label.size)
  mtext( text = "Mean wave height", side = 2, line = 2.0, cex = axis.label.size)
  text(labels = "F", y = ylimits[1] + 0.94*diff(ylimits), x = 1983.7, cex = 1.4)
#
# Panel G: Aeshna numbers
par(new = "TRUE", plt = c(x.margins[5], x.margins[6], 0.13, 0.48))
  ylimits <- c(0.97, 1.02) * range(log(conditions$aeshna.year), na.rm = TRUE)
  plot(c(1,2), c(1,2), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = xlimits, ylim = ylimits, xaxs="i", yaxs="i")
  points(x = conditions$year, y = log(conditions$aeshna.year), pch = 20, cex = point.size, col = color.list$col[3])
  lines(x = conditions$year, y = log(conditions$aeshna.year), lwd = line.size, col = color.list$col[3])
  box(lwd = 1.2)
  axis(side = 1, labels = FALSE, at = seq(1983,1998,5), tck = -0.025)
  axis(side = 1, labels = FALSE, at = c(1984,1985,1986,1987,1989,1990,1991,1992,1994,1995,1996,1997), tck = -0.018)
  mtext(text = c("1983","1988","1993","1998"), side = 1, line = 0.4, las = 1, at = seq(1983,1998,5), cex = tick.label.size)
  mtext(text = "Year", side = 1, line = 1.5, cex = axis.label.size)
  axis(side = 2, labels = FALSE, at = log(c(100,200,300,400,500)), tck = -0.025)
  mtext(text = seq(100,500,100), side = 2, line = 0.35, las = 1, at = log(seq(100,500,100)), cex = tick.label.size)
  mtext( text = "Nr. of dragonflies", side = 2, line = 1.85, cex = axis.label.size)
  text(labels = "G", y = ylimits[1] + 0.94*diff(ylimits), x = 1983.7, cex = 1.4)
#
# Panel H: Tadpole numbers
par(new = "TRUE", plt = c(x.margins[7], x.margins[8], 0.13, 0.48))
  ylimits <- c(0.96, 1.02) * range(conditions$l.Ntads.year, na.rm = TRUE)
  plot(c(1,2), c(1,2), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = xlimits, ylim = ylimits, xaxs="i", yaxs="i")
  points(x = conditions$year, y = conditions$l.Ntads.year, pch = 20, cex = point.size, col = color.list$col[4])
  lines(x = conditions$year, y = conditions$l.Ntads.year, lwd = line.size, col = color.list$col[4])
  box(lwd = 1.2)
  axis(side = 1, labels = FALSE, at = seq(1983,1998,5), tck = -0.025)
  axis(side = 1, labels = FALSE, at = c(1984,1985,1986,1987,1989,1990,1991,1992,1994,1995,1996,1997), tck = -0.018)
  mtext(text = c("1983","1988","1993","1998"), side = 1, line = 0.4, las = 1, at = seq(1983,1998,5), cex = tick.label.size)
  mtext(text = "Year", side = 1, line = 1.5, cex = axis.label.size)
  axis(side = 2, labels = FALSE, at = log(c(500,1000,5000)), tck = -0.025)
  axis(side = 2, labels = FALSE, at = log(c(400,600,700,800,900,2000,3000,4000,6000,7000,8000)), tck = -0.018)
  mtext(text = c(500,1000,5000), side = 2, line = 0.35, las = 1, at = log(c(500,1000,5000)), cex = tick.label.size)
  mtext(text = "Nr. of tadpoles", side = 2, line = 2, cex = axis.label.size)
  text(labels = "H", y = ylimits[1] + 0.94*diff(ylimits), x = 1983.7, cex = 1.4)
dev.off()
#
#
# For each brood, determine whether another brood appeared nearby on the same night.
#
setwd(data.wd)
dist.obj.to.df <- function(d){
  # Make a data.frame out of a distance object.
  size <- attr(d, "Size")
  return(data.frame(subset(expand.grid(row = 2:size, col = 1:(size-1)), row > col), distance = as.numeric(d), row.names = NULL))
  }
pool.coords   <- read.table("Pool XYZ coordinates.txt", header = TRUE)
pool.coords   <- subset(pool.coords, pool > 1000)
pool.labels   <- pool.coords$pool - 1000
dist.obj      <- dist(pool.coords[,c("X", "Y")])
pool.dist     <- dist.obj.to.df(dist.obj)
pool.dist$row <- pool.labels[pool.dist$row]
pool.dist$col <- pool.labels[pool.dist$col]
day.window    <- 1            # hatching date must be within this number of days to potentially belong to the same brood.
nearest.pools <- data.frame(pool = pool.labels, nearest.pool.1 = 0, nearest.pool.2 = 0, nearest.pool.3 = 0, dist.nearest.1 = 0, dist.nearest.2 = 0, dist.nearest.3 = 0)
count         <- 0
for (i in pool.labels) {
  count     <- count + 1
  distances <- pool.dist[pool.dist$row == i | pool.dist$col == i, ]
  distances <- distances[order(distances$distance), ][1:3,]
  nearest.pools[count, 2:4] <- apply(distances[ , 1:2], 1, sum) - i
  nearest.pools[count, 5:7] <- distances[ , 3]
  nearest.pools[count, 1]   <- i
  }
res   <- data.frame(year= 1983:1998, nr.broods = NA, nr.1.m = 0, nr.2.m = 0, nr.3.m = 0, nr.4.m = 0, nr.nearest.pool.1 = 0, nr.nearest.pool.2 = 0, nr.nearest.pool.3 = 0)
count <- 0
for (j in 1983:1998) {
  count <- count + 1
  d.sub <- fc[fc$year == j,]
  res$nr.broods[count] <- dim(d.sub)[1]
  for (i in 1:dim(d.sub)[1]) {              # start looping through each brood in this year
    d.sub$day.diff <- d.sub$est.day.3.5[i] - d.sub$est.day.3.5
    for (k in 1:dim(d.sub)[1]) {              # start looping through all other broods in this year
      if (i != k) {                                  # is this a different brood
        if (!is.na(d.sub$day.diff[k]) & abs(d.sub$day.diff[k]) < day.window) {   # did this brood hatch within 1 day?
          which.this.pool    <- which.min(abs(d.sub$pool[i] - nearest.pools$pool))   # for some reason, match did not work.
          if (d.sub$pool[k] == nearest.pools[which.this.pool, 2]) res$nr.nearest.pool.1[count] <- res$nr.nearest.pool.1[count] + 1    # if this brood is in the closest pool
          if (d.sub$pool[k] == nearest.pools[which.this.pool, 3]) res$nr.nearest.pool.2[count] <- res$nr.nearest.pool.2[count] + 1    # if this brood is in the second closest pool
          if (d.sub$pool[k] == nearest.pools[which.this.pool, 4]) res$nr.nearest.pool.3[count] <- res$nr.nearest.pool.3[count] + 1    # if this brood is in the third closest pool
          distances <- c(pool.dist$distance[d.sub$pool[i] == pool.dist$row & d.sub$pool[k] == pool.dist$col], pool.dist$distance[d.sub$pool[i] == pool.dist$col & d.sub$pool[k] == pool.dist$row])
          if (length(distances) > 0) {
            if (distances <= 1 ) res$nr.1.m[count] <- res$nr.1.m[count] + 1
            if (distances <= 2 ) res$nr.2.m[count] <- res$nr.2.m[count] + 1
            if (distances <= 3 ) res$nr.3.m[count] <- res$nr.3.m[count] + 1
            if (distances <= 4 ) res$nr.4.m[count] <- res$nr.4.m[count] + 1
            } 
          } # end if this brood hatched within 1 day.
        }   # end if this is not the same brood.
      }     # end k looping through other broods.
    }       # end i looping through all broods.
  }         # end j looping through years.
res[,3:9] <- res[,3:9] / res[,2]
res[,8]   <- res[,7] + res[,8]
res[,9]   <- res[,8] + res[,9]
colnames(res) <- c("year", "broods", "p.1.m", "p.2.m", "p.3.m", "p.4.m", "p.nearest.pool", "p.2.nearest.pools", "p.3.nearest.pools")
res
#
#
###################
#
# Now prepare the individual-level dataset.
# Each line is an individual tadpole caught at the beginning of the first capture interval.
#
# Read in the individual growth rates for all tadpoles that survive the growth interval.
# Assumes that rank-order remains the same during the growth interval.
# This was produced directly from data frame "fc" (the broods dataset) using "2014-10-13 individual tadpole growth rates.R"
#
ind.growth        <- read.table("Growth assuming rank order.txt", header = FALSE)
names(ind.growth) <- c("year", "pool", "brood", "interval", "growth")
fc$pooldens       <- fc$npool1 / fc$init.area
indivs            <- merge(fc[,c("year", "pool", "brood", "interval", "loc", "max.area", "origin", "first.day", "fate", "est.day.3.5", "nbrood1", "nbrood2", "survival", "pooldens")], ind.growth, by = c("year", "pool", "brood", "interval"), all = TRUE)
indivs$ind.surv   <- 1
fc1          <- fc
fc1$nr.diff  <- fc1$nbrood1 - fc1$nbrood2
zeros        <- data.frame()       # Tadpoles that died should have survival = zero
for (i in 1:dim(fc1)[1]) {
  if (!is.na(fc1$nr.diff[i])) {    # If we visited the pool at the begining and end of the interval.
    if (fc1$nr.diff[i] > 0) {      # If at least one tadpole died during the interval.
      temp <- data.frame(ind.surv = rep(0, fc1$nr.diff[i]), growth = NA, year = fc1$year[i], pool = fc1$pool[i], brood = fc1$brood[i], interval = fc1$interval[i], loc = fc1$loc[i], max.area = fc1$max.area[i], origin = fc1$origin[i], first.day = fc1$first.day[i], fate = fc1$fate[i], est.day.3.5 = fc1$est.day.3.5[i], nbrood1 = fc1$nbrood1[i], nbrood2 = fc1$nbrood2[i], pooldens = fc1$pooldens[i], survival = fc1$survival[i])
      zeros <- rbind(zeros, temp)
      }
    }
  }
fc2      <- rbind(zeros, indivs[, c("year", "pool", "brood", "interval", "loc", "max.area", "origin", "first.day", "fate", "est.day.3.5", "nbrood1", "nbrood2", "pooldens", "survival", "ind.surv", "growth")])
fc2      <- fc2[order(fc2$year, fc2$pool, fc2$brood, fc2$interval),]
fc2$date <- fc2$est.day.3.5
fc2$area <- log(fc2$max.area)
#
# Prepare fitness measure
#
fc2$StandGrowth   <- (fc2$growth - mean(fc2$growth, na.rm = TRUE)) / sd(fc2$growth, na.rm = TRUE)
fc2$StandFirstDay <- (fc2$date - mean(fc2$date, na.rm = TRUE)) / sd(fc2$date, na.rm = TRUE)
logitP            <- -1.969 + 0.6013*fc2$StandGrowth - 0.1978*fc2$StandFirstDay       # using Dave Smith's (1987, Ecology) data to get "fitness".
fc2$fitness       <- exp(logitP) / (1 + exp(logitP))                                  # this is estimated survival to reproductive maturity at 2 years of age, assuming metamorphic date is proportional to hatching date and metamorphic size is proportional to growth.
fc2$fitness       <- ifelse(fc2$ind.surv == 0, 0, fc2$fitness)
fc2$rel.fitness   <- fc2$fitness / mean(fc2$fitness, na.rm=TRUE)
fc2$round.fitness <- round( fc2$rel.fitness * 10 )
#
# Fig S3 -- zero-inflated distribution of fitness.
#
setwd(figs.wd)
pdf("Fig S3.pdf", useDingbats = FALSE, width = 5.0, height = 5)
plot.new()
par(new = "TRUE", plt = c(0.25, 0.85, 0.25, 0.85))
  hist(fc2$round.fitness, breaks=c(50), col = main.color, lwd = 1, xlab = NULL, ylab = NULL, xaxt="n", yaxt="n", main = NA, xlim = c(0,80) )
  box(lwd = 1.3)
  axis(side = 1, labels = FALSE, tck = -0.025)
  mtext(text = c("0","20","40","60","80"), side = 1, line = 0.4, las = 1, at = c(0,20,40,60,80), cex = 1)
  axis(side = 2, labels = FALSE, tck = -0.025)
  mtext(text = c("0","5,000","10,000","15,000"), side = 2, line = 0.5, las = 1, at = c(0,5000,10000,15000), cex = 1)
  mtext( text = "Expected fitness", side = 1, line = 1.6, cex = 1.3)
  mtext( text = "Number of tadpoles", side = 2, line = 3.3, cex = 1.3)
  text(c("N = 42,172 tadpoles"), x = 25, y = 15000, pos = 4, cex = 1)
  mtext(text = paste("Van Buskirk & Smith, Figure S3,", Sys.Date()), side = 3, line = 2.8, cex = 0.7)
dev.off()
#
#
#
################################
################################
#
# Selection analyses
#
#
years   <- unique(fc2$year)
years   <- years[order(years)]
#
# Prepare data for selection models.
#
fc3 <- fc2[ , c("year", "pool", "brood", "interval", "loc", "area", "date", "fitness", "rel.fitness", "round.fitness", "pooldens")]
fc3 <- rename(fc3, "area", "larea")
fc3$lfitness <- log(1 + fc3$fitness)
fc3 <- fc3[complete.cases(fc3), ]
#
# standardize the covariates
mean.loc.indivs.kept    <- mean(fc3$loc, na.rm = TRUE)
sd.loc.indivs.kept      <- sd(fc3$loc, na.rm = TRUE)
fc3$loc.s               <- (fc3$loc - mean.loc.indivs.kept) / sd.loc.indivs.kept
sd.loc.s2.indivs.kept   <- sd(fc3$loc.s^2, na.rm = TRUE)
fc3$loc.s2              <- fc3$loc.s^2 / sd.loc.s2.indivs.kept
mean.larea.indivs.kept  <- mean(fc3$larea, na.rm = TRUE)
sd.larea.indivs.kept    <- sd(fc3$larea, na.rm = TRUE)
fc3$larea.s             <- (fc3$larea- mean.larea.indivs.kept) / sd.larea.indivs.kept
sd.larea.s2.indivs.kept <- sd(fc3$larea.s^2, na.rm = TRUE)
fc3$larea.s2            <- fc3$larea.s^2 / sd.larea.s2.indivs.kept
mean.date.indivs.kept   <- mean(fc3$date, na.rm = TRUE)
sd.date.indivs.kept     <- sd(fc3$date, na.rm = TRUE)
fc3$date.s              <- (fc3$date- mean.date.indivs.kept) / sd.date.indivs.kept
sd.date.s2.indivs.kept  <- sd(fc3$date.s^2, na.rm = TRUE)
fc3$date.s2             <- fc3$date.s^2 / sd.date.s2.indivs.kept
fc3 <- fc3[ , c("year", "pool", "brood", "interval", "loc", "larea", "date", "loc.s", "loc.s2", "larea.s", "larea.s2", "date.s", "date.s2", "lfitness", "rel.fitness", "round.fitness", "pooldens")]
setwd(models.wd)
#
#
# Begin with the method of Chevin et al. (2015, Evolution).
# Estimate whether fitness optimum moves around from one year to the next (and also
# the directional selection coefficient and maximum fitness).
#
# This is done separately for each of the three traits.
#
# Each model takes 7-8 hours.
#
# Pool location.
start.time <- Sys.time()
set.seed(501)
t.loc <- brm(bf(lfitness ~ loc.s + loc.s2
          + ( 1 | pool ) + ( 1 | year )
          + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year ) ),
          family = gaussian(),
          control = list(max_treedepth = 20),
          iter = 1600,
          warmup = 600,
          cores = 2,
          data = fc3,
          file = "t_loc" )
end.time <- Sys.time()
end.time - start.time
allz.t  <- fc1[ , c("year", "nbrood1", "loc")]
allz.t  <- allz.t[complete.cases(allz.t),]
allz.t$loc.s <- (allz.t$loc - mean.loc.indivs.kept) / sd.loc.indivs.kept
meanz.t <- aggregate(allz.t$loc.s, by = list(year = allz.t$year), mean)
C1      <- Chevin.2015(ps = posterior_samples(t.loc), meanz.t, allz.t = allz.t[,c("year", "loc.s")], outliers = TRUE, threshold = 1)
#
#
# Pool size.
start.time <- Sys.time()
set.seed(502)
t.area <- brm(bf(lfitness ~ larea.s + larea.s2
          + ( 1 | pool ) + ( 1 | year )
          + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year ) ),
          family = gaussian(),
          control = list(max_treedepth = 20),
          iter = 1600,
          warmup = 600,
          cores = 2,
          data = fc3,
          file = "t_area" )
end.time <- Sys.time()
end.time - start.time
allz.t  <- fc1[ , c("year", "nbrood1", "area")]
allz.t  <- allz.t[complete.cases(allz.t),]
allz.t$larea.s <- (allz.t$area - mean.larea.indivs.kept) / sd.larea.indivs.kept
meanz.t <- aggregate(allz.t$larea.s, by = list(year = allz.t$year), mean)
C2      <- Chevin.2015(ps = posterior_samples(t.area), meanz.t, allz.t = allz.t[,c("year", "larea.s")], outliers = TRUE, threshold = 1)
#
#
# Hatching date.
start.time <- Sys.time()
set.seed(503)
t.date <- brm(bf(lfitness ~ date.s + date.s2
          + ( 1 | pool ) + ( 1 | year )
          + ( 0 + date.s | year ) + ( 0 + date.s2 | year ) ),
          family = gaussian(),
          control = list(max_treedepth = 20),
          iter = 1600,
          warmup = 600,
          cores = 2,
          data = fc3,
          file = "t_date" )
end.time <- Sys.time()
end.time - start.time
allz.t        <- fc1[ , c("year", "nbrood1", "date")]
allz.t        <- allz.t[complete.cases(allz.t),]
allz.t$date.s <- (allz.t$date - mean.date.indivs.kept) / sd.date.indivs.kept
meanz.t <- aggregate(allz.t$date.s, by = list(year = allz.t$year), mean)
C3      <- Chevin.2015(ps = posterior_samples(t.date), meanz.t, allz.t = allz.t[,c("year", "date.s")], outliers = TRUE, threshold = 1)
#
#
# Figure S9 -- Evolutionary loads for the three traits.
# Calculated as per Chevin et et al. (2015) and de Villemereuil et al. (2020, PNAS).
#
# There are two kinds of load. One is the proportional decline in mean population
# fitness due to the fact that the mean trait value is not at the optimum (theta).
# This could be due to lag (if theta is moving directionally), fluctuations in theta,
# or other causes of maladaptation. This load is what is calculated here as yearly
# load. The total load.
# The load due to fluctuating theta is given on page 31975 of de Villemereuil et al.
# (2020) as var(theta) / 2*(omega^2 + 1). Var(theta) is the variance among years.
# This quantity is stated to be "proportional to the decline in log mean fitness."
# Both types of load are calculated here as medians +/- HPDI.
#
get.loads <- function(C, quantiles) {
  yearly.load           <- apply(C, 2, quantile, quantiles, na.rm = TRUE)
  yearly.load           <- apply(yearly.load, 2, format.p.val, 6)
  yearly.load           <- as.data.frame(apply(yearly.load, 2, as.numeric))
  rownames(yearly.load) <- c("lci", "mean.load", "uci")
  yearly.load           <- cbind(year = c(1983,1985:1998), t(yearly.load))
  return(list(year.load = yearly.load, mean.loads = as.numeric(format.p.val(apply(yearly.load[,c(2,3,4)], 2, mean), 6))))
  }
fun1 <- function(x) { format.p.val(x, 6) }
#
q1 <- quantile(C1$fluctuation.load, c(0.025, 0.5, 0.975), na.rm = TRUE)
q2 <- quantile(C2$fluctuation.load, c(0.025, 0.5, 0.975), na.rm = TRUE)
q3 <- quantile(C3$fluctuation.load, c(0.025, 0.5, 0.975), na.rm = TRUE)
quantiles <- c(0.1, 0.5, 0.9)
tl.loc  <- get.loads(C = C1$total.load, quantiles)                                    # total load
ll.loc  <- get.loads(C = C1$lag.load, quantiles)                                      # lag load
fl.loc  <- as.numeric(fun1(quantile(C1$fluctuation.load, quantiles, na.rm = TRUE)))   # lag load due to fluctuation theta
tl.area <- get.loads(C = C2$total.load, quantiles)
ll.area <- get.loads(C = C2$lag.load, quantiles)
fl.area <- as.numeric(fun1(quantile(C2$fluctuation.load, quantiles, na.rm = TRUE)))
tl.date <- get.loads(C = C3$total.load, quantiles)
ll.date <- get.loads(C = C3$lag.load, quantiles)
fl.date <- as.numeric(fun1(quantile(C3$fluctuation.load, quantiles, na.rm = TRUE)))
#   To get off the log scale to absolute fitness decline:    1 - 1/exp(yl.date$mean.loads)
#
# What fraction of the lag load is due to the fluctuating optimum?
median(C1$fluctuation.load, na.rm = TRUE) / ll.loc$mean.loads[2]     # location
median(C2$fluctuation.load, na.rm = TRUE) / ll.area$mean.loads[2]    # area
median(C3$fluctuation.load, na.rm = TRUE) / ll.date$mean.loads[2]    # date
#
# What fraction of the total load is due to the fluctuating optimum?
median(C1$fluctuation.load, na.rm = TRUE) / tl.loc$mean.loads[2]     # location
median(C2$fluctuation.load, na.rm = TRUE) / tl.area$mean.loads[2]    # area
median(C3$fluctuation.load, na.rm = TRUE) / tl.date$mean.loads[2]    # date
#
# Fig S9.
#
draw.plot <- function(tl, ll, fl, xx, error.line.width, colors) {
  points(xx[1], tl[2], cex = 1.7, pch = 19, col = colors[1])
  lines(x = c(xx[1], xx[1]), y = c(tl[1], tl[3]), lwd = error.line.width, col = colors[1])
  points(xx[2], ll[2], cex = 1.7, pch = 19, col = colors[2])
  lines(x = c(xx[2], xx[2]), y = c(ll[1], ll[3]), lwd = error.line.width, col = colors[2])
  points(xx[3], fl[2], cex = 1.7, pch = 19, col = colors[3])
  lines(x = c(xx[3], xx[3]), y = c(fl[1], fl[3]), lwd = error.line.width, col = colors[3])
  }
interval        <- 0.31
tick.width      <- 1.2
tick.label.size <- 0.9
axis.label.size <- 1
plotxlimits     <- c(0, 12)
x.seq           <- c(1,2,3,5,6,7,9,10,11)
plotylimits     <- c(0,0.046)
tl.color <- "#FFBB45"
ll.color <- "#C57D00"
fl.color <- "mediumpurple4"
setwd(figs.wd)
pdf(paste("Fig S9.pdf", sep = ""), useDingbats = FALSE, width = 6.5, height = 3)
plot.new()
par(new = "TRUE", plt = c(0.15, 0.9, 0.12, 0.87))
  plot(c(1,2), c(1,2), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = plotxlimits, ylim = plotylimits, xaxs="i", yaxs="i")
  mtext(text = c("Pool location", "Pool size", "Hatching date"), at = c(2, 6, 10), side = 1, line = 0.5, cex = axis.label.size)
  axis(side = 2, tck = -0.018, lwd = tick.width, labels = FALSE)
  mtext(text = c("0.00","0.01","0.02","0.03","0.04"), side = 2, line = 0.38, at = c(0, 0.01, 0.02, 0.03, 0.04), las = 1, cex = tick.label.size)
  mtext(text = "Load (prop. decline in log-fitness)", side = 2, line = 2.4, cex = axis.label.size)
  draw.plot(tl = tl.loc$mean.loads, ll = ll.loc$mean.loads, fl = fl.loc, xx = x.seq[1:3], error.line.width = 1.1, colors = c(tl.color, ll.color,fl.color))
  draw.plot(tl = tl.area$mean.loads, ll = ll.area$mean.loads, fl = fl.area, xx = x.seq[4:6], error.line.width = 1.1, colors = c(tl.color, ll.color,fl.color))
  draw.plot(tl = 100*tl.date$mean.loads, ll = 100*ll.date$mean.loads, fl = 100*fl.date, xx = x.seq[7:9], error.line.width = 1.1, colors = c(tl.color, ll.color,fl.color))
  interval <- 0.003
  tops     <- c(0.028, 0.007, 0.008)
  text(c("total", "load", "lag", "load", "fluctuating", "optimum"), x = c(1.7,1.7,1.3,1.3,4.1,4.1), y = c(tops[1], tops[1]-interval, tops[2], tops[2]-interval, tops[3], tops[3]-interval), cex = 0.8, col = c(tl.color, tl.color, ll.color, ll.color, fl.color, fl.color))
  box(lwd = tick.width)
  mtext(text = paste("Van Buskirk & Smith, Figure S9,", Sys.Date()), side = 3, line = 1.2, cex = 0.6)
dev.off()
#
#
# Figure 5 -- estimated optimum (theta) and the observed trait distribution.
#
quantiles        <- c(0.1, 0.9)
wt.means         <- aggregate.weighted.mean(x = c("loc", "area", "date"), w = "nbrood1", by = "year", dat = fc1)
wt.means$loc.s   <- (wt.means$wtMean.loc - mean.loc.indivs.kept) / sd.loc.indivs.kept
wt.means$area.s  <- (wt.means$wtMean.area - mean.larea.indivs.kept) / sd.larea.indivs.kept
wt.means$date.s  <- (wt.means$wtMean.date - mean.date.indivs.kept) / sd.date.indivs.kept
# prepare quantiles of the three traits.
for (i in 1:dim(fc1)[1]) {
  d <- data.frame(count = 1:fc1$nbrood1[i], year = fc1$year[i], loc = fc1$loc[i], larea = fc1$area[i], date = fc1$date[i])
  if (i == 1) {
    res <- d
    } else {
    res <- rbind(res, d)
    }
  }
res$loc   <- (res$loc - mean.loc.indivs.kept) / sd.loc.indivs.kept
res$larea <- (res$larea - mean.larea.indivs.kept) / sd.larea.indivs.kept
res$date  <- (res$date - mean.date.indivs.kept) / sd.date.indivs.kept
trait.quants <- aggregate(res[,c("loc", "larea", "date")], by = list(year = res$year), mean, na.rm = TRUE)
names(trait.quants) <- c("year", "mean.loc", "mean.larea", "mean.date")
trait.quants <- merge(trait.quants, aggregate(res[,c("loc", "larea")], by = list(year = res$year), quantile, quantiles), by = "year")
res1         <- res[complete.cases(res),]
trait.quants <- merge(trait.quants, aggregate(res1[,"date"], by = list(year = res1$year), quantile, quantiles), by = "year")
trait.quants <- cbind(trait.quants[,c(1:4)], cbind(cbind(as.data.frame(trait.quants$loc), as.data.frame(trait.quants$larea)), as.data.frame(trait.quants$x)))
names(trait.quants)[5:10] <- c("loc.lci", "loc.uci", "larea.lci", "larea.uci", "date.lci", "date.uci")
#
drawing.function <- function(C, trait.dist, trait.quants, quantiles, years, main.color, main.color.lite, trait.color) {
  # Get the distribution of estimated optima from the fitted model.
  theta.quants <- data.frame(i = 1:length(years), year = years, median.theta = NA, lci = NA, uci = NA)
  for (i in 1:length(years)) {
    theta.quants[i,3]   <- median(C$theta[,i], na.rm = TRUE)
    theta.quants[i,4:5] <- quantile(C$theta[,i], quantiles, na.rm = TRUE)
    }
  # Draw the posterior distribution of model-estimated theta-t values.
  polygon(x = c(theta.quants$year, rev(theta.quants$year)), y = c(theta.quants$lci, rev(theta.quants$uci)), col = main.color.lite, border = NA)
  lines(x = theta.quants$year, y = theta.quants$median.theta, lwd = 2, col = main.color)
  # Draw the mean and quantiles of the trait values in the population
  points(x = trait.quants$year, y = trait.quants[,2], cex = 0.85, pch = 19, col = trait.color)
  for (i in 1:length(trait.quants$year)) {
    lines(x = c(trait.quants$year[i], trait.quants$year[i]), y = c(trait.quants[i,3], trait.quants[i,4]), lwd = 1.3, col = trait.color)
    }
  }   # end function
#
xlimits          <- c(-10,10)
plotxlimits      <- c(1982.5, 1998.5)
trait.color      <- "purple"  # lower.color
main.color.lite  <- adjustcolor(main.color, alpha.f = 0.2)
tick.width       <- 1.3
tick.label.size  <- 1
axis.label.size  <- 1.2
#
setwd(figs.wd)
pdf(paste("Fig 5.pdf", sep = ""), useDingbats = FALSE, width = 10, height = 4.5)
plot.new()
  # 
  # Panel A, Pool location.
  trait.dist <- fc3[, c("year", "loc.s")]
  wtMeans    <- wt.means[, c("year", "loc.s")]
  C          <- C1
  C$theta[C$theta > xlimits[2]] <- xlimits[2]
  C$theta[C$theta < xlimits[1]] <- xlimits[1]
  plotylimits     <- c(-2.5, 2.7)
  par(new = "TRUE", plt = c(0.08, 0.33, 0.15, 0.84))
  plot(c(1,2), c(1,2), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = plotxlimits, ylim = plotylimits, xaxs="i", yaxs="i")
  axis(side = 1, tck = -0.023, lwd = tick.width, labels = FALSE)
  axis(side = 1, tck = -0.015, lwd = tick.width, labels = FALSE, at = c(1983,1984,1986:1989,1991:1994, 1996:1998))
  mtext(text = c(1985,1990,1995), side = 1, line = 0.35, at = c(1985,1990,1995), cex = tick.label.size)
  mtext(text = c("Year"), side = 1, line = 1.4, cex = axis.label.size)
  axis(side = 2, tck = -0.023, lwd = tick.width, at = (c(0.2, 0.4, 0.6, 0.8) - mean.loc.indivs.kept) / sd.loc.indivs.kept, labels = FALSE)
  axis(side = 2, tck = -0.015, lwd = tick.width, at = (c(0.3, 0.5, 0.7, 0.9) - mean.loc.indivs.kept) / sd.loc.indivs.kept, labels = FALSE)
  mtext(text = c(0.2,0.4,0.6,0.8), side = 2, line = 0.35, at = (c(0.2, 0.4, 0.6, 0.8) - mean.loc.indivs.kept) / sd.loc.indivs.kept, las = 1, cex = tick.label.size)
  mtext(text = expression(paste("Value of trait or ", theta)), side = 2, line = 2.35, cex = axis.label.size)
  # arrows at the left side.
  par(xpd = TRUE)
  mtext(text = c("Lake", "Forest"), side = 2, line = 1.5, at = c(plotylimits[1]+0.85, plotylimits[2]-0.85), adj = c(0,1), cex = tick.label.size-0.1, las = 0)
  arrows(x0 = plotxlimits[1]-2.5, y0 = plotylimits[1]+0.72, x1 = plotxlimits[1]-2.5, y1 = plotylimits[1]+0.08, angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  arrows(x0 = plotxlimits[1]-2.5, y0 = plotylimits[2]-0.72, x1 = plotxlimits[1]-2.5, y1 = plotylimits[2]-0.08, angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  par(xpd = FALSE)
  drawing.function(C, trait.dist, trait.quants = trait.quants[,c(1,2,5,6)], quantiles, years = c(1983,1985:1998), main.color, main.color.lite, trait.color)
  box(lwd = tick.width)
  mtext(text = "A. Pool location", side = 3, line = 0.3, at = 1982.8, adj = 0, cex = tick.label.size+0.1)
  mtext(text = paste("Van Buskirk & Smith, Figure 5,", Sys.Date()), side = 3, line = 2.2, cex = 0.6)
  # 
  # Panel B, Pool size.
  trait.dist <- fc3[, c("year", "larea.s")]
  C          <- C2
  C$theta[C$theta > xlimits[2]] <- xlimits[2]
  C$theta[C$theta < xlimits[1]] <- xlimits[1]
  plotylimits     <- c(-2.84569, 2.3)
  par(new = "TRUE", plt = c(0.37, 0.62, 0.15, 0.84))
  plot(c(1,2), c(1,2), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = plotxlimits, ylim = plotylimits, xaxs="i", yaxs="i")
  axis(side = 1, tck = -0.023, lwd = tick.width, labels = FALSE)
  axis(side = 1, tck = -0.015, lwd = tick.width, labels = FALSE, at = c(1983,1984,1986:1989,1991:1994, 1996:1998))
  mtext(text = c(1985,1990,1995), side = 1, line = 0.35, at = c(1985,1990,1995), cex = tick.label.size)
  mtext(text = c("Year"), side = 1, line = 1.4, cex = axis.label.size)
  axis(side = 2, tck = -0.023, lwd = tick.width, at = (log(c(0.1,1,10)) - mean.larea.indivs.kept) / sd.larea.indivs.kept, labels = FALSE)
  axis(side = 2, tck = -0.015, lwd = tick.width, at = (log(c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2:9)) - mean.larea.indivs.kept) / sd.larea.indivs.kept, labels = FALSE)
  mtext(text = c("0.1", "1.0", "10"), side = 2, line = 0.35, at = (log(c(0.1,1,10)) - mean.larea.indivs.kept) / sd.larea.indivs.kept, las = 1, cex = tick.label.size)
  drawing.function(C, trait.dist, trait.quants[,c(1,3,7:8)], quantiles, years = c(1983,1985:1998), main.color, main.color.lite, trait.color)
  box(lwd = tick.width)
  mtext(text = "B. Pool size (m  )", side = 3, line = 0.3, at = 1982.8, adj = 0, cex = tick.label.size+0.1)
  mtext(text = "2", side = 3, line = 0.55, at = 1990.6, adj = 0, cex = tick.label.size-0.15)
  # 
  # Panel C, Hatching date.
  trait.dist <- fc3[, c("year", "date.s")]
  C          <- C3
  C$theta[C$theta > xlimits[2]] <- xlimits[2]
  C$theta[C$theta < xlimits[1]] <- xlimits[1]
  plotylimits     <- c(-3.5, 2.5)
  par(new = "TRUE", plt = c(0.69, 0.94, 0.15, 0.84))
  plot(c(1,2), c(1,2), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = plotxlimits, ylim = plotylimits, xaxs="i", yaxs="i")
  axis(side = 1, tck = -0.023, lwd = tick.width, labels = FALSE)
  axis(side = 1, tck = -0.015, lwd = tick.width, labels = FALSE, at = c(1983,1984,1986:1989,1991:1994, 1996:1998))
  mtext(text = c(1985,1990,1995), side = 1, line = 0.35, at = c(1985,1990,1995), cex = tick.label.size)
  mtext(text = c("Year"), side = 1, line = 1.4, cex = axis.label.size)
  axis(side = 2, tck = -0.023, lwd = tick.width, at = (c(20,40,60) - mean.date.indivs.kept) / sd.date.indivs.kept, labels = FALSE)
  axis(side = 2, tck = -0.015, lwd = tick.width, at = (c(30,50,70) - mean.date.indivs.kept) / sd.date.indivs.kept, labels = FALSE)
  mtext(text = c(paste(make.dates(20)[[1]], make.dates(20)[[2]], sep = " "), paste(make.dates(40)[[1]], make.dates(40)[[2]], sep = " "), paste(make.dates(60)[[1]], make.dates(60)[[2]], sep = " ")), side = 2, line = 0.35, at = (c(20,40,60) - mean.date.indivs.kept) / sd.date.indivs.kept, las = 1, cex = tick.label.size)
  drawing.function(C, trait.dist, trait.quants[,c(1,4,9:10)], quantiles, years = c(1983,1985:1998), main.color, main.color.lite, trait.color)
  box(lwd = tick.width)
  mtext(text = "C. Hatching date", side = 3, line = 0.3, at = 1982.8, adj = 0, cex = tick.label.size+0.1)
  text(c("observed", "trait", "value"), x = 1987.7, y = c(2.25, 1.91, 1.57), cex = 0.85, col = trait.color)
  text(c("estimated", "optimum"), x = 1993, y = c(-2.4, -2.745), cex = 0.85, col = main.color)
dev.off()
#
#
# Lande and Arnold approach.
#
# Begin with reduced model that has only linear terms. This is used to estimate directional
# selection gradients and their heterogeneity among years. Non-linear selection will be
# estimated later using the full regression model (Eqn 1 in the MS).
# This approach is recommended by Lande and Arnold (1983, p. 1218, lower left).
#
# Each model takes about 8 min to run. Once you have run the following code once, the 
# fitted glmmTMB model objects will be loaded from disk thereafter.
#
setwd(models.wd)
tC <- tryCatch.W.E( load("m1_DIRGRAD.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  start.time <- Sys.time()
  m1.DIRGRAD <- glmmTMB(round.fitness ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ) ,
            zi = ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.DIRGRAD, file = "m1_DIRGRAD.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_DIRGRAD.txt") }
#
tC <- tryCatch.W.E( load("m1_DIRGRAD.COND.no.year.by.date.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.DIRGRAD.COND.no.year.by.date <- glmmTMB(round.fitness ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) ,
            zi = ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.DIRGRAD.COND.no.year.by.date, file = "m1_DIRGRAD.COND.no.year.by.date.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_DIRGRAD.COND.no.year.by.date.txt") }
#
tC <- tryCatch.W.E( load("m1_DIRGRAD.COND.no.year.by.area.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.DIRGRAD.COND.no.year.by.area <- glmmTMB(round.fitness ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + date.s | year ) ,
            zi = ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.DIRGRAD.COND.no.year.by.area, file = "m1_DIRGRAD.COND.no.year.by.area.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_DIRGRAD.COND.no.year.by.area.txt") }
#
tC <- tryCatch.W.E( load("m1_DIRGRAD.COND.no.year.by.loc.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.DIRGRAD.COND.no.year.by.loc <- glmmTMB(round.fitness ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ) ,
            zi = ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.DIRGRAD.COND.no.year.by.loc, file = "m1_DIRGRAD.COND.no.year.by.loc.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_DIRGRAD.COND.no.year.by.loc.txt") }
#
tC <- tryCatch.W.E( load("m1_DIRGRAD.COND.only.year.pool.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  start.time <- Sys.time()
  m1.DIRGRAD.COND.only.year.pool <- glmmTMB(round.fitness ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) ,
            zi = ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.DIRGRAD.COND.only.year.pool, file = "m1_DIRGRAD.COND.only.year.pool.txt" )
  end.time <- Sys.time()
  end.time - start.time
  }
if(is.null(attr(tC$value, "class"))) { load("m1_DIRGRAD.COND.only.year.pool.txt") }
#
tC <- tryCatch.W.E( load("m1_DIRGRAD.COND.no.year.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.DIRGRAD.COND.no.year <- glmmTMB(round.fitness ~ loc.s + larea.s + date.s
            + ( 1 | pool )  ,
            zi = ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.DIRGRAD.COND.no.year, file = "m1_DIRGRAD.COND.no.year.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_DIRGRAD.COND.no.year.txt") }
#
tC <- tryCatch.W.E( load("m1_DIRGRAD.COND.no.pool.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.DIRGRAD.COND.no.pool <- glmmTMB(round.fitness ~ loc.s + larea.s + date.s
            + ( 1 | year ),
            zi = ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.DIRGRAD.COND.no.pool, file = "m1_DIRGRAD.COND.no.pool.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_DIRGRAD.COND.no.pool.txt") }
#
tC <- tryCatch.W.E( load("m1_DIRGRAD.ZI.no.year.by.date.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.DIRGRAD.ZI.no.year.by.date <- glmmTMB(round.fitness ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ) ,
            zi = ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) ,
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.DIRGRAD.ZI.no.year.by.date, file = "m1_DIRGRAD.ZI.no.year.by.date.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_DIRGRAD.ZI.no.year.by.date.txt") }
#
tC <- tryCatch.W.E( load("m1_DIRGRAD.ZI.no.year.by.area.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.DIRGRAD.ZI.no.year.by.area <- glmmTMB(round.fitness ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ) ,
            zi = ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.DIRGRAD.ZI.no.year.by.area, file = "m1_DIRGRAD.ZI.no.year.by.area.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_DIRGRAD.ZI.no.year.by.area.txt") }
#
tC <- tryCatch.W.E( load("m1_DIRGRAD.ZI.no.year.by.loc.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.DIRGRAD.ZI.no.year.by.loc <- glmmTMB(round.fitness ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ) ,
            zi = ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.DIRGRAD.ZI.no.year.by.loc, file = "m1_DIRGRAD.ZI.no.year.by.loc.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_DIRGRAD.ZI.no.year.by.loc.txt") }
#
tC <- tryCatch.W.E( load("m1_DIRGRAD.ZI.only.year.pool.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  start.time <- Sys.time()
  m1.DIRGRAD.ZI.only.year.pool <- glmmTMB(round.fitness ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ) ,
            zi = ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) ,
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.DIRGRAD.ZI.only.year.pool, file = "m1_DIRGRAD.ZI.only.year.pool.txt" )
  end.time <- Sys.time()
  end.time - start.time
  }
if(is.null(attr(tC$value, "class"))) { load("m1_DIRGRAD.ZI.only.year.pool.txt") }
#
tC <- tryCatch.W.E( load("m1_DIRGRAD.ZI.no.year.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.DIRGRAD.ZI.no.year <- glmmTMB(round.fitness ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ) ,
            zi = ~ loc.s + larea.s + date.s
            + ( 1 | pool ) ,
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.DIRGRAD.ZI.no.year, file = "m1_DIRGRAD.ZI.no.year.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_DIRGRAD.ZI.no.year.txt") }
#
tC <- tryCatch.W.E( load("m1_DIRGRAD.ZI.no.pool.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  start.time <- Sys.time()
  m1.DIRGRAD.ZI.no.pool <- glmmTMB(round.fitness ~ loc.s + larea.s + date.s
            + ( 1 | pool ) + ( 1 | year ) + ( 0 + loc.s | year ) + ( 0 + larea.s | year ) + ( 0 + date.s | year ) ,
            zi = ~ loc.s + larea.s + date.s
            + ( 1 | year ) ,
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.DIRGRAD.ZI.no.pool, file = "m1_DIRGRAD.ZI.no.pool.txt" )
  end.time <- Sys.time()
  end.time - start.time
  }
if(is.null(attr(tC$value, "class"))) { load("m1_DIRGRAD.ZI.no.pool.txt") }
#
#
# Table 1 (part): directional selection gradients.
#
nr.digits <- 3
a1 <- anova(m1.DIRGRAD, m1.DIRGRAD.COND.no.year.by.date)             # test year * date
a2 <- anova(m1.DIRGRAD, m1.DIRGRAD.COND.no.year.by.area)             # test year * area
a3 <- anova(m1.DIRGRAD, m1.DIRGRAD.COND.no.year.by.loc)              # test year * loc
a4 <- anova(m1.DIRGRAD.COND.only.year.pool, m1.DIRGRAD.COND.no.year) # test year
a5 <- anova(m1.DIRGRAD.COND.only.year.pool, m1.DIRGRAD.COND.no.pool) # test pool
b1 <- anova(m1.DIRGRAD, m1.DIRGRAD.ZI.no.year.by.date)               # test year * date
b2 <- anova(m1.DIRGRAD, m1.DIRGRAD.ZI.no.year.by.area)               # test year * area
b3 <- anova(m1.DIRGRAD, m1.DIRGRAD.ZI.no.year.by.loc)                # test year * loc
b4 <- anova(m1.DIRGRAD.ZI.only.year.pool, m1.DIRGRAD.ZI.no.year)     # test year
b5 <- anova(m1.DIRGRAD.ZI.only.year.pool, m1.DIRGRAD.ZI.no.pool)     # test pool
#
m1.COND.RE <- rbind(c(a5$Chisq[2], a5[[8]][2]), c(a4$Chisq[2], a4[[8]][2]), c(a3$Chisq[2], a3[[8]][2]), c(a2$Chisq[2], a2[[8]][2]), c(a1$Chisq[2], a1[[8]][2]))
m1.ZI.RE   <- rbind(c(b5$Chisq[2], b5[[8]][2]), c(b4$Chisq[2], b4[[8]][2]), c(b3$Chisq[2], b3[[8]][2]), c(b2$Chisq[2], b2[[8]][2]), c(b1$Chisq[2], b1[[8]][2]))
m1.COND.RE[,c(1)] <- format.p.val(m1.COND.RE[,c(1)], 1)
m1.COND.RE[,c(2)] <- format.p.val(m1.COND.RE[,c(2)], 4)
m1.ZI.RE[,c(1)]   <- format.p.val(m1.ZI.RE[,c(1)], 1)
m1.ZI.RE[,c(2)]   <- format.p.val(m1.ZI.RE[,c(2)], 4)
m1.CI             <- confint(m1.DIRGRAD)
m1.COND.RE        <- cbind(format.p.val(m1.CI[5:9,3]^2, ndigits = (nr.digits + 1)), m1.COND.RE)
m1.ZI.RE          <- cbind(format.p.val(m1.CI[14:18,3]^2, ndigits = (nr.digits + 1)), m1.ZI.RE)
m1.CI[,c(1,2,3)]  <- apply(m1.CI[,c(1,2,3)], 2, format.p.val, nr.digits)
glmmTMB.summary.table <- data.frame(Model.part = c("Logistic", rep("", 5), "Poisson", rep("", 5)),
 Source = c("", "Intercept", "Loc", "Size", "Date", "", "", "Intercept", "Loc", "Size", "Date", ""),
 Estimate = c("", paste(m1.CI[10:13,3], " (", m1.CI[10:13,1], ",", m1.CI[10:13,2], ")", sep = ""), "", "", paste(m1.CI[1:4,3], " (", m1.CI[1:4,1], ",", m1.CI[1:4,2], ")", sep = ""), ""),
 RandEffect = c("", "Pool", "Year", "Location(Year)", "Area(Year)", "Date(Year)", "", "Pool", "Year", "Location(Year)", "Area(Year)", "Date(Year)"),
 VarComp = c("", m1.ZI.RE[,1], "", m1.COND.RE[,1]),
 LR.stat = c("", m1.ZI.RE[,2], "", m1.COND.RE[,2]), stringsAsFactors = FALSE)
#
glmmTMB.summary.table   # Coefficients for the directional selection gradients.
#
#
# Full model for getting the quadratic and correlational selection gradients.
# Each model takes about 18 minutes.
#
tC <- tryCatch.W.E( load("m1_cor_full_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1 <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ) ,
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1, file = "m1_cor_full_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_full_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_COND_no.yr.by.area.by.date_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.COND.no.yr.by.area.by.date <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) ,
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE),
            data = fc3)
  save(m1.COND.no.yr.by.area.by.date, file = "m1_cor_COND_no.yr.by.area.by.date_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_COND_no.yr.by.area.by.date_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_COND_no.yr.by.loc.by.date_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.COND.no.yr.by.loc.by.date <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + larea.s:date.s | year ) ,
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE),
            data = fc3)
  save(m1.COND.no.yr.by.loc.by.date, file = "m1_cor_COND_no.yr.by.loc.by.date_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_COND_no.yr.by.loc.by.date_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_COND_no.yr.by.loc.by.area_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.COND.no.yr.by.loc.by.area <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ) ,
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.COND.no.yr.by.loc.by.area, file = "m1_cor_COND_no.yr.by.loc.by.area_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_COND_no.yr.by.loc.by.area_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_COND_no.yr.by.higher.order_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.COND.no.yr.by.higher.order <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE),
            data = fc3)
  save(m1.COND.no.yr.by.higher.order, file = "m1_cor_COND_no.yr.by.higher.order_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_COND_no.yr.by.higher.order_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_COND_no.yr.by.date2_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.COND.no.yr.by.date2 <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.COND.no.yr.by.date2, file = "m1_cor_COND_no.yr.by.date2_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_COND_no.yr.by.date2_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_COND_no.yr.by.date_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.COND.no.yr.by.date <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.COND.no.yr.by.date, file = "m1_cor_COND_no.yr.by.date_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_COND_no.yr.by.date_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_COND_no.yr.by.area2_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.COND.no.yr.by.area2 <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.COND.no.yr.by.area2, file = "m1_cor_COND_no.yr.by.area2_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_COND_no.yr.by.area2_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_COND_no.yr.by.area_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.COND.no.yr.by.area <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE),
            data = fc3)
  save(m1.COND.no.yr.by.area, file = "m1_cor_COND_no.yr.by.area_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_COND_no.yr.by.area_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_COND_no.yr.by.loc2_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.COND.no.yr.by.loc2 <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE),
            data = fc3)
    save(m1.COND.no.yr.by.loc2, file = "m1_cor_COND_no.yr.by.loc2_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_COND_no.yr.by.loc2_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_COND_no.yr.by.loc_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.COND.no.yr.by.loc <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE),
            data = fc3)
  save(m1.COND.no.yr.by.loc, file = "m1_cor_COND_no.yr.by.loc_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_COND_no.yr.by.loc_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_ZI_no.yr.by.area.by.date_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.ZI.no.yr.by.area.by.date <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) ,
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE),
            data = fc3)
  save(m1.ZI.no.yr.by.area.by.date, file = "m1_cor_ZI_no.yr.by.area.by.date_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_ZI_no.yr.by.area.by.date_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_ZI_no.yr.by.loc.by.date_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.ZI.no.yr.by.loc.by.date <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + larea.s:date.s | year ) ,
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE),
            data = fc3)
  save(m1.ZI.no.yr.by.loc.by.date, file = "m1_cor_ZI_no.yr.by.loc.by.date_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_ZI_no.yr.by.loc.by.date_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_ZI_no.yr.by.loc.by.area_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.ZI.no.yr.by.loc.by.area <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ) ,
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.ZI.no.yr.by.loc.by.area, file = "m1_cor_ZI_no.yr.by.loc.by.area_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_ZI_no.yr.by.loc.by.area_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_ZI_no.yr.by.higher.order_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.ZI.no.yr.by.higher.order <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE),
            data = fc3)
  save(m1.ZI.no.yr.by.higher.order, file = "m1_cor_ZI_no.yr.by.higher.order_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_ZI_no.yr.by.higher.order_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_ZI_no.yr.by.date2_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.ZI.no.yr.by.date2 <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.ZI.no.yr.by.date2, file = "m1_cor_ZI_no.yr.by.date2_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_ZI_no.yr.by.date2_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_ZI_no.yr.by.date_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.ZI.no.yr.by.date <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  # This one returned a false convergence warning. Selecting profile=TRUE solved the problem.
  save(m1.ZI.no.yr.by.date, file = "m1_cor_ZI_no.yr.by.date_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_ZI_no.yr.by.date_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_ZI_no.yr.by.area2_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.ZI.no.yr.by.area2 <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.ZI.no.yr.by.area2, file = "m1_cor_ZI_no.yr.by.area2_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_ZI_no.yr.by.area2_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_ZI_no.yr.by.area_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.ZI.no.yr.by.area <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE),
            data = fc3)
  save(m1.ZI.no.yr.by.area, file = "m1_cor_ZI_no.yr.by.area_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_ZI_no.yr.by.area_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_ZI_no.yr.by.loc2_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.ZI.no.yr.by.loc2 <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE),
            data = fc3)
    save(m1.ZI.no.yr.by.loc2, file = "m1_cor_ZI_no.yr.by.loc2_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_ZI_no.yr.by.loc2_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_ZI_no.yr.by.loc_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.ZI.no.yr.by.loc <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year )
            + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year )
            + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
            + ( 0 + date.s | year )  + ( 0 + date.s2 | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE),
            data = fc3)
  save(m1.ZI.no.yr.by.loc, file = "m1_cor_ZI_no.yr.by.loc_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_ZI_no.yr.by.loc_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_only.yr.and.pool_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.only.yr.and.pool <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year ),
            family = poisson,
            control=glmmTMBControl(profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  #  this one returned two warnings: Choosing profile=TRUE fixed the problem:
  save(m1.only.yr.and.pool, file = "m1_cor_only.yr.and.pool_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_only.yr.and.pool_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_COND_only.yr_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.COND.only.yr <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year ),
            family = poisson,
            data = fc3)
  save(m1.COND.only.yr, file = "m1_cor_COND_only.yr_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_COND_only.yr_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_COND_only.pool_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.COND.only.pool <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  save(m1.COND.only.pool, file = "m1_cor_COND_only.pool_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_COND_only.pool_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_ZI_only.yr_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.ZI.only.yr <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | year ),
            family = poisson,
            data = fc3)
  save(m1.ZI.only.yr, file = "m1_cor_ZI_only.yr_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_ZI_only.yr_ind.txt") }
#
tC <- tryCatch.W.E( load("m1_cor_ZI_only.pool_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.ZI.only.pool <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ) + ( 1 | year ),
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
            + ( 1 | pool ),
            family = poisson,
            control=glmmTMBControl(parallel=2, profile=TRUE, optimizer=optim, optArgs=list(method="BFGS")),
            data = fc3)
  # This one returned several warnings: Choosing profile=TRUE made them go away.
  save(m1.ZI.only.pool, file = "m1_cor_ZI_only.pool_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_cor_ZI_only.pool_ind.txt") }
#
# Prepare summary table. Table 1, in part.
#
nr.digits <- 3
#
a1  <- anova(m1, m1.COND.no.yr.by.area.by.date)                        # test year * area * date
a2  <- anova(m1, m1.COND.no.yr.by.loc.by.date)                         # test year * loc * date
a3  <- anova(m1, m1.COND.no.yr.by.loc.by.area)                         # test year * loc * area
a4  <- anova(m1.COND.no.yr.by.higher.order, m1.COND.no.yr.by.loc2)     # test year * loc2
a5  <- anova(m1.COND.no.yr.by.loc2, m1.COND.no.yr.by.loc)              # test year * loc
a6  <- anova(m1.COND.no.yr.by.higher.order, m1.COND.no.yr.by.area2)    # test year * area2
a7  <- anova(m1.COND.no.yr.by.area2, m1.COND.no.yr.by.area)            # test year * area
a8  <- anova(m1.COND.no.yr.by.higher.order, m1.COND.no.yr.by.date2)    # test year * date2
a9  <- anova(m1.COND.no.yr.by.date2, m1.COND.no.yr.by.date)            # test year * date
a10 <- anova(m1.only.yr.and.pool, m1.COND.only.pool)                   # test year
a11 <- anova(m1.only.yr.and.pool, m1.COND.only.yr)                     # test pool
#
b1  <- anova(m1, m1.ZI.no.yr.by.area.by.date)                          # test year * area * date
b2  <- anova(m1, m1.ZI.no.yr.by.loc.by.date)                           # test year * loc * date
b3  <- anova(m1, m1.ZI.no.yr.by.loc.by.area)                           # test year * loc * area
b4  <- anova(m1.ZI.no.yr.by.higher.order, m1.ZI.no.yr.by.loc2)         # test year * loc2
b5  <- anova(m1.ZI.no.yr.by.loc2, m1.ZI.no.yr.by.loc)                  # test year * loc
b6  <- anova(m1.ZI.no.yr.by.higher.order, m1.ZI.no.yr.by.area2)        # test year * area2
b7  <- anova(m1.ZI.no.yr.by.area2, m1.ZI.no.yr.by.area)                # test year * area
b8  <- anova(m1.ZI.no.yr.by.higher.order, m1.ZI.no.yr.by.date2)        # test year * date2
b9  <- anova(m1.ZI.no.yr.by.date2, m1.ZI.no.yr.by.date)                # test year * date
b10 <- anova(m1.only.yr.and.pool, m1.ZI.only.pool)                     # test year
b11 <- anova(m1.only.yr.and.pool, m1.ZI.only.yr)                       # test pool
#
m1.COND.RE <- rbind(c(a11$Chisq[2], a11[[8]][2]), c(a10$Chisq[2], a10[[8]][2]), c(a5$Chisq[2], a5[[8]][2]), c(a4$Chisq[2], a4[[8]][2]), c(a7$Chisq[2], a7[[8]][2]), c(a6$Chisq[2], a6[[8]][2]), c(a9$Chisq[2], a9[[8]][2]), c(a8$Chisq[2], a8[[8]][2]), c(a3$Chisq[2], a3[[8]][2]), c(a2$Chisq[2], a2[[8]][2]), c(a1$Chisq[2], a1[[8]][2]))
m1.ZI.RE   <- rbind(c(b11$Chisq[2], b11[[8]][2]), c(b10$Chisq[2], b10[[8]][2]), c(b5$Chisq[2], b5[[8]][2]), c(b4$Chisq[2], b4[[8]][2]), c(b7$Chisq[2], b7[[8]][2]), c(b6$Chisq[2], b6[[8]][2]), c(b9$Chisq[2], b9[[8]][2]), c(b8$Chisq[2], b8[[8]][2]), c(b3$Chisq[2], b3[[8]][2]), c(b2$Chisq[2], b2[[8]][2]), c(b1$Chisq[2], b1[[8]][2]))
m1.COND.RE[,c(1)] <- format.p.val(m1.COND.RE[,c(1)], 1)
m1.COND.RE[,c(2)] <- format.p.val(m1.COND.RE[,c(2)], 2)
m1.ZI.RE[,c(1)]   <- format.p.val(m1.ZI.RE[,c(1)], 1)
m1.ZI.RE[,c(2)]   <- format.p.val(m1.ZI.RE[,c(2)], 2)
m1.CI             <- confint(m1)
m1.COND.RE        <- cbind(format.p.val(m1.CI[11:21,3]^2, ndigits = (nr.digits + 1)), m1.COND.RE)
m1.ZI.RE          <- cbind(format.p.val(m1.CI[32:42,3]^2, ndigits = (nr.digits + 1)), m1.ZI.RE)
m1.CI[,c(1,2,3)]  <- apply(m1.CI[,c(1,2,3)], 2, format.p.val, nr.digits)
glmmTMB.summary.table <- data.frame(Model.part = c("Logistic", rep("", 11), "Poisson", rep("", 11)), Source = c("", "Intercept", "Loc", "Loc^2", "Size", "Size^2", "Date", "Date^2", "Loc:Size", "Loc:Date", "Size:Date", "", "", "Intercept", "Loc", "Loc^2", "Size", "Size^2", "Date", "Date^2", "Loc:Size", "Loc:Date", "Size:Date", ""), Estimate = c("", paste(m1.CI[22:31,3], " (", m1.CI[22:31,1], ",", m1.CI[22:31,2], ")", sep = ""), "", "", paste(m1.CI[1:10,3], " (", m1.CI[1:10,1], ",", m1.CI[1:10,2], ")", sep = ""), ""),  RandEffect = c("", "Pool", "Year", "Location(Year)", "Location^2(Year)", "Sire(Year)", "Size^2(Year)", "Date(Year)", "Date^2(Year)", "Loc:Size(Year)", "Loc:Date(Year)", "Size:Date(Year)", "", "Pool", "Year", "Location(Year)", "Location^2(Year)", "Sire(Year)", "Size^2(Year)", "Date(Year)", "Date^2(Year)", "Loc:Size(Year)", "Loc:Date(Year)", "Size:Date(Year)"), VarComp = c("", m1.ZI.RE[,1], "", m1.COND.RE[,1]), LR.stat = c("", m1.ZI.RE[,2], "", m1.COND.RE[,2]), stringsAsFactors = FALSE)
glmmTMB.summary.table   # Summary table for the non-linear terms.
#
#
# Figure 2 -- where and when the eggs appeared.
#
setwd(data.wd)
fieldwork.dates <- read.table("Dates of fieldwork 1983-1998.txt", header = TRUE)  # generally first and last dates on which we measured pool depths; not necessarily dates when we were on the study area.
# Fix situations in which broods were discovered before or after we have pool depth data
fieldwork.dates$last.date  <- ifelse(fieldwork.dates$year == 1985, 70, fieldwork.dates$last.date)
fieldwork.dates$last.date  <- ifelse(fieldwork.dates$year == 1990, 89, fieldwork.dates$last.date)
fieldwork.dates$last.date  <- ifelse(fieldwork.dates$year == 1992, 74, fieldwork.dates$last.date)
fieldwork.dates$last.date  <- ifelse(fieldwork.dates$year == 1996, 90, fieldwork.dates$last.date)
fieldwork.dates$first.date <- ifelse(fieldwork.dates$year == 1995, 27, fieldwork.dates$first.date)
fieldwork.dates$first.date <- ifelse(fieldwork.dates$year == 1997, 22, fieldwork.dates$first.date)
fieldwork.dates$first.date <- ifelse(fieldwork.dates$year == 1998, 21, fieldwork.dates$first.date)
nr.broods        <- aggregate(fc[,c("location", "max.area", "est.day.3.5")], by = list(fc$year, fc$pool), length)
nr.broods        <- aggregate(nr.broods$location, by = list(nr.broods$Group.2), mean)
names(nr.broods) <- c("pool", "mean.nr.broods")
nr.broods        <- merge(nr.broods, pools, by = "pool", all = TRUE)
unused.pools     <- nr.broods[ is.na(nr.broods$mean.nr.broods), ]
unused.pools$mean.nr.broods <- NULL
unused.pools     <- unused.pools[ complete.cases(unused.pools), ]
used.pools       <- nr.broods[ complete.cases(nr.broods), ]
# Fit a quadratic model to the nr-of-broods data.
nr.broods$mean.nr.broods <- NA.to.zero(nr.broods$mean.nr.broods)
nr.broods$location2 <- nr.broods$location^2
nr.broods$logarea   <- log(nr.broods$max.area)
nr.broods$logarea2  <- nr.broods$logarea^2
m.surf           <- lm(mean.nr.broods ~ location + location2 + logarea + logarea2 + location:logarea, data = nr.broods)
xlimits          <- c(0, 1.02)
ylimits          <- c(-0.15, 0.1) + range(log(pools$max.area), na.rm = TRUE)
newdat           <- expand.grid(logarea = seq(ylimits[1], ylimits[2], length.out = 30), location = seq(xlimits[1], xlimits[2], length.out = 30))
newdat$location2 <- newdat$location^2
newdat$logarea2  <- newdat$logarea^2
newdat$pred      <- predict(m.surf, newdat)
years            <- unique(fc$year)
years            <- years[order(years)]
#
# Logistic regression to predict appearance of eggs in a pool as a function of loc,
# area, and date.
#
fc1          <- fc
pool.by.year <- table(fc$year, fc$pool, round(fc$date))        # Table counting up all the broods in each pool (columns) by year (row) by date
pool.choice  <- data.frame(year=integer(), pool=numeric(), loc=numeric(), area=numeric(), date=numeric(), nr.broods=integer())
pools        <- pools[ complete.cases(pools), ]
for (i in 1:length(fieldwork.dates$year)) {
  temp <- expand.grid(year = fieldwork.dates$year[i], pool = pools$pool, loc = NA, area = NA, date = seq(fieldwork.dates$first.date[i], fieldwork.dates$last.date[i], 1), nr.broods = 0)
  pool.choice <- rbind(pool.choice, temp)
  }
which.pools      <- match(pool.choice$pool, pools$pool)
pool.choice$loc  <- pools$location[which.pools]
pool.choice$area <- pools$max.area[which.pools]
fc1      <- fc1[, c("est.day.3.5", "pool", "year")]
fc1$date <- round(fc$est.day.3.5)
fc1      <- fc1[complete.cases(fc1), ]
for (i in 1:dim(fc1)[1]) { pool.choice$nr.broods[fc1$date[i] == pool.choice$date & fc1$pool[i] == pool.choice$pool & fc1$year[i] == pool.choice$year] <- 1  }
# standardize exactly as for the selection analysis.
pool.choice$loc.s    <- (pool.choice$loc - mean.loc.indivs.kept) / sd.loc.indivs.kept
pool.choice$loc.s2   <- pool.choice$loc.s^2 / sd.loc.s2.indivs.kept
pool.choice$larea    <- log(pool.choice$area)
pool.choice$larea.s  <- (pool.choice$larea - mean.larea.indivs.kept) / sd.larea.indivs.kept
pool.choice$larea.s2 <- pool.choice$larea.s^2 / sd.larea.s2.indivs.kept
pool.choice$date.s   <- (pool.choice$date - mean.date.indivs.kept) / sd.date.indivs.kept
pool.choice$date.s2  <- pool.choice$date.s^2
#
# Predict pool use as a function of loc, area, and date: takes 21 min.
setwd(models.wd)
tC <- tryCatch.W.E( load("m_choice.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m.choice <- glmmTMB(nr.broods ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s
          + ( 1 | pool ) + ( 1 | year )
          + ( 0 + loc.s | year ) + ( 0 + loc.s2 | year )
          + ( 0 + larea.s | year ) + ( 0 + larea.s2 | year )
          + ( 0 + date.s | year ) + ( 0 + date.s2 | year )
          + ( 0 + loc.s:larea.s | year ) + ( 0 + loc.s:date.s | year ) + ( 0 + larea.s:date.s | year ) ,
          family = binomial, data = pool.choice)
  save(m.choice, file = "m_choice.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m_choice.txt") }
#
# Predict from "m.choice" across a grid.
area.buffer  <- c(0.955, 1.06)
grid.density <- 25
df.predict          <- expand.grid(year = NA, pool = NA, loc = seq(xlimits[1], xlimits[2], length.out = grid.density), larea = seq(ylimits[1], ylimits[2], length.out = grid.density), date = mean(fc$est.day.3.5, na.rm = TRUE))
df.predict$loc.s    <- (df.predict$loc - mean.loc.indivs.kept) / sd.loc.indivs.kept
df.predict$loc.s2   <- df.predict$loc.s^2 / sd.loc.s2.indivs.kept
df.predict$larea.s  <- (df.predict$larea - mean.larea.indivs.kept) / sd.larea.indivs.kept
df.predict$larea.s2 <- df.predict$larea.s^2 / sd.larea.s2.indivs.kept
df.predict$date.s   <- (df.predict$date - mean.date.indivs.kept) / sd.date.indivs.kept
df.predict$date.s2  <- df.predict$date.s^2 # / sd.larea.s2.indivs.kept
df.predict          <- cbind(df.predict, Predict.choice = predict(m.choice, type = "response", re.form = NA, newdata = df.predict))
#
bubble.color.pale  <- adjustcolor(main.color, alpha.f = 0.5)
axis.label.size    <- 1.15
contour.label.size <- 1
tick.label.size    <- 0.9
magnif             <- 1.9      # size of bubbles
color.gradient     <- function(x, colors = c(lower.color, lowmid.color, highmid.color, upper.color), colsteps = 100) { return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )  }
plot.colors        <- c(color.gradient(25:100), rep(plot.colors[100], 24))
setwd(figs.wd)
pdf("Fig 2.pdf", useDingbats = FALSE, width = 5.5, height = 8.6)
plot.new()
# Top panel (A) is mean number of broods per year in every pool.
par(new = "TRUE", plt = c(0.25, 0.85, 0.575, 0.92))
  filled.contour3(x = unique(df.predict$loc), y = unique(df.predict$larea), z = matrix(df.predict$Predict.choice, nrow = length(unique(df.predict$loc)), ncol = length(unique(df.predict$larea)), byrow = TRUE), levels = seq(min(df.predict$Predict.choice), max(df.predict$Predict.choice), length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  used.pools <- used.pools[order(used.pools$mean.nr.broods, decreasing = TRUE), ]
  # overlay contours and white-out the areas that should not have contours
  contour(x = unique(df.predict$loc), y = unique(df.predict$larea), z = matrix(df.predict$Predict.choice, nrow = length(unique(df.predict$loc)), ncol = length(unique(df.predict$larea)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, nlevels = 6)
  points(used.pools$location, log(used.pools$max.area), cex = 1, col = "black", xaxt = "n", yaxt = "n", bg = "red", pch = 21 )
  points(unused.pools$location, log(unused.pools$max.area), col = "black", bg = "white",  cex = 0.8, pch = 21 )
  polygon(x = c(0.75, xlimits[2], xlimits[2]), y = log(c(0.02, 0.2, 0.02)), border = NA, col = "white")
  polygon(x = c(xlimits[1], xlimits[1], 0.4), y = log(c(0.018, 0.08, 0.018)), border = NA, col = "white")
  polygon(x = c(xlimits[1], xlimits[1], 0.33), y = c(log(1.6), ylimits[2], ylimits[2]), border = NA, col = "white")
  polygon(x = c(0.86, xlimits[2], xlimits[2]), y = c(ylimits[2], ylimits[2], log(2)), border = NA, col = "white")
  box(lwd = 1.3)
  draw.contour.labels(xlimits = xlimits, ylimits, xfraction = 0.105, yfraction = 0.06, x.corners = c(0.23, 0.47, 0.8), y.corners = log(c(2, 2.8, 1.8)), contour.labels = c("0.01", "0.02", "0.03"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  axis(side = 1, labels = FALSE, tck = -0.022)
  mtext(text = c("0.0","0.2","0.4","0.6","0.8","1.0"), side = 1, line = 0.4, las = 1, at = c(0,0.2,0.4,0.6,0.8,1), cex = tick.label.size)
  mtext(text = "Pool location", side = 1, line = 1.4, cex = axis.label.size)
  axis(side = 2, at = log(c(0.01,0.1,1,10)), labels = FALSE, tck = -0.022)
  axis(side = 2, at = log(c(0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20)), labels = FALSE, tck = -0.011)
  mtext(text = c("0.1","1.0","10"), side = 2, line = 0.4, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext( text = "Pool surface area (m  )", side = 2, line = 2.0, cex = axis.label.size)
  mtext( text = "2", at = log(4.8), side = 2, line = 2.4, cex = 0.9)
  text(labels = "A", x = 0.925, y = log(13), pos = 4, cex = 1.4)
  # arrows at the bottom
  par(xpd = TRUE)
  mtext(text = c("Lake", "Forest"), side = 1, line = 1.3, at = c(0.11, 0.9), adj = c(0,1), cex = tick.label.size, las = 0)
  arrows(x0 = 0.1, y0 = log(0.0115), x1 = 0.0, y1 = log(0.0115), angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  arrows(x0 = 0.91, y0 = log(0.0115), x1 = 1.01, y1 = log(0.0115), angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  par(xpd = FALSE)
  mtext( text = paste("Van Buskirk & Smith, Figure 2,", Sys.Date()), side = 3, line = 1.5, cex = 0.7)
#
# Middle panel (B) is histogram of estimated hatching dates.
par(new = "TRUE", plt = c(0.25, 0.85, 0.42, 0.50))
  hist(fc$est.day.3.5, breaks=c(80), col = main.color, lwd = 1, xlab = NULL, ylab = NULL, xaxt="n", yaxt="n", main = NA, xlim = c(15,93) )
  box(lwd = 1.3)
  axis(side = 1, at = seq(20,80,20), labels = FALSE, tck = -0.045)
  axis(side = 2, labels = FALSE, tck = -0.045)
  mtext(text = c("0","20","40"), side = 2, line = 0.5, las = 1, at = c(0,20,40), cex = tick.label.size)
  mtext( text = "Freq.", side = 2, line = 2.0, cex = axis.label.size)
  text(labels = "B", x = 87.5, y = 30, pos = 4, cex = 1.4)
#
# Bottom panel (C) is coverage in time with broods marked as ticks.
par(new = "TRUE", plt = c(0.25, 0.85, 0.09, 0.40))
  plot(c(1,2), c(1983,1998), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = c(15,93))
  for(i in 1983:1998) {
    pres <- subset(fieldwork.dates, year == i)
    rect(xleft = pres$first.date-1, ybottom = i-0.3, xright = pres$last.date+1, ytop = i+0.3, col = bubble.color.pale, border = NA)
    }
 # points(fc$est.day.3.5, fc$year, cex = 0.6, pch = 19)
  for(i in 1:length(years)) {
    d.sub <- subset(fc, year == years[i])[ , c("year", "est.day.3.5")]
    d.sub <- d.sub[complete.cases(d.sub),]
    for(j in 1:dim(d.sub)[1]) {
      lines(x = c(d.sub$est.day.3.5[j], d.sub$est.day.3.5[j]), y = c(years[i]-0.3, years[i]+0.3))
      }
    }
  brood.sample <- aggregate(fc$est.day.3.5, by = list(fc$year), FUN = function(x) samp.size(x))[,2][,1]
  text(labels = brood.sample, x = 16, y = 1983:1998, cex = 0.75)   # sample size along the left
  text(labels = "C", x = 87.5, y = 1997.6, pos = 4, cex = 1.4)
  box(lwd = 1.3)
  axis(side = 1, at = seq(20,80,20), labels = FALSE, tck = -0.022)
  x.text <- c(paste(make.dates(20), collapse = " "), paste(make.dates(40), collapse = " "), paste(make.dates(60), collapse = " "), paste(make.dates(80), collapse = " "))
  mtext(text = x.text, side = 1, line = 0.4, las = 1, at = seq(20,80,20), cex = tick.label.size)
  mtext(text = "Hatching date", side = 1, line = 1.55, cex = axis.label.size)
  axis(side = 2, at = seq(1983,1998, 1), labels = FALSE, tck = -0.022)
  mtext(text = seq(1983,1998, 3), side = 2, line = 0.4, las = 1, at = seq(1983,1998, 3), cex = tick.label.size)
  mtext( text = "Year", side = 2, line = 2.2, cex = axis.label.size)
dev.off()
#
#
# What causes annual variation in selection? Drying, wave wash, competition, and dragonflies.
#
d <- merge(merge(merge(merge(fc3, total.aeshna, by = "year", all = TRUE), waves, by = "year", all = TRUE), dry, by = "year", all = TRUE), dens4, by = "year", all = TRUE)
d <- d[d$year < 1999, ]
d$total.aeshna     <- ifelse(d$year < 1986, NA, d$total.aeshna)
d$l.pooldens       <- log(d$pooldens)
d$l.total.aeshna   <- log(d$total.aeshna)
d$l.Nr.tads        <- log(d$Nr.tads)
mean.lpooldens.kept <- mean(d$l.pooldens, na.rm = TRUE)
mean.aeshna.kept   <- mean(d$l.total.aeshna, na.rm = TRUE)
mean.meanwave.kept <- mean(d$mean.wave.ht, na.rm = TRUE)
mean.propdry.kept  <- mean(d$prop.dry, na.rm = TRUE)
mean.meandens.kept <- mean(d$l.meandens, na.rm = TRUE)
mean.Nrtads.kept   <- mean(d$l.Nr.tads, na.rm = TRUE)
sd.lpooldens.kept  <- sd(d$l.pooldens, na.rm = TRUE)
sd.aeshna.kept     <- sd(d$l.total.aeshna, na.rm = TRUE)
sd.meanwave.kept   <- sd(d$mean.wave.ht, na.rm = TRUE)
sd.propdry.kept    <- sd(d$prop.dry, na.rm = TRUE)
sd.meandens.kept   <- sd(d$l.meandens, na.rm = TRUE)
sd.Nrtads.kept     <- sd(d$l.Nr.tads, na.rm = TRUE)
d$lpooldens.c <- (d$l.pooldens - mean.lpooldens.kept) / sd.lpooldens.kept  # standardize covariates
d$aeshna.c    <- (d$l.total.aeshna - mean.aeshna.kept) / sd.aeshna.kept
d$meanwave.c  <- (d$mean.wave.ht - mean.meanwave.kept) / sd.meanwave.kept
d$propdry.c   <- (d$prop.dry - mean.propdry.kept) / sd.propdry.kept
d$meandens.c  <- (d$l.meandens - mean.meandens.kept) / sd.meandens.kept
d$Nr.tads.c   <- (d$l.Nr.tads - mean.Nrtads.kept) / sd.Nrtads.kept
rbind(apply(d, 2, mean, na.rm = TRUE), apply(d, 2, sd, na.rm = TRUE))
#
# glmmTMB models testing how ecological factors modify selection.
# Each model takes about 45 minutes to run the first time; thereafter, model object will be loaded from disk.
#
setwd(models.wd)
tC <- tryCatch.W.E( load("m1_drying_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.drying <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s + propdry.c + propdry.c:loc.s + propdry.c:loc.s2 + propdry.c:larea.s + propdry.c:larea.s2 + propdry.c:date.s + propdry.c:date.s2 + propdry.c:loc.s:larea.s + propdry.c:loc.s:date.s + propdry.c:larea.s:date.s
            + ( 1 | pool ) ,
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s + propdry.c + propdry.c:loc.s + propdry.c:loc.s2 + propdry.c:larea.s + propdry.c:larea.s2 + propdry.c:date.s + propdry.c:date.s2 + propdry.c:loc.s:larea.s + propdry.c:loc.s:date.s + propdry.c:larea.s:date.s
            + ( 1 | pool ) ,
            family = poisson,
            control = glmmTMBControl(parallel = 2, profile = TRUE, optimizer = optim, optArgs = list(method = "BFGS")),
            data = d)
  save(m1.drying, file = "m1_drying_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_drying_ind.txt") }
#
#
tC <- tryCatch.W.E( load("m1_mean.wave_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.mean.wave <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s + meanwave.c + meanwave.c:loc.s + meanwave.c:loc.s2 + meanwave.c:larea.s + meanwave.c:larea.s2 + meanwave.c:date.s + meanwave.c:date.s2 + meanwave.c:loc.s:larea.s + meanwave.c:loc.s:date.s + meanwave.c:larea.s:date.s
            + ( 1 | pool ) ,
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s + meanwave.c + meanwave.c:loc.s + meanwave.c:loc.s2 + meanwave.c:larea.s + meanwave.c:larea.s2 + meanwave.c:date.s + meanwave.c:date.s2 + meanwave.c:loc.s:larea.s + meanwave.c:loc.s:date.s + meanwave.c:larea.s:date.s
            + ( 1 | pool ) ,
            family = poisson,
            control = glmmTMBControl(parallel = 4, profile = TRUE, optimizer = optim, optArgs = list(method = "BFGS")),
            data = d)
  save(m1.mean.wave, file = "m1_mean.wave_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_mean.wave_ind.txt") }
#
#
d.sub <- d[ ,c("pool", "loc.s", "loc.s2", "larea.s", "larea.s2", "date.s", "date.s2", "round.fitness", "aeshna.c")]
d.sub <- d.sub[complete.cases(d.sub), ]
tC <- tryCatch.W.E( load("m1_Aeshna_ind.txt") ) 
if(!is.null(attr(tC$value, "class"))) {
  m1.Aeshna <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s + aeshna.c + aeshna.c:loc.s + aeshna.c:loc.s2 + aeshna.c:larea.s + aeshna.c:larea.s2 + aeshna.c:date.s + aeshna.c:date.s2 + aeshna.c:loc.s:larea.s + aeshna.c:loc.s:date.s + aeshna.c:larea.s:date.s
            + ( 1 | pool ) ,
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s + aeshna.c + aeshna.c:loc.s + aeshna.c:loc.s2 + aeshna.c:larea.s + aeshna.c:larea.s2 + aeshna.c:date.s + aeshna.c:date.s2 + aeshna.c:loc.s:larea.s + aeshna.c:loc.s:date.s + aeshna.c:larea.s:date.s
            + ( 1 | pool ) ,
            family = poisson,
            control = glmmTMBControl(parallel = 4, profile = TRUE, optimizer = optim, optArgs = list(method = "BFGS")),
            data = d.sub)
  save(m1.Aeshna, file = "m1_Aeshna_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_Aeshna_ind.txt") }
#
#
tC <- tryCatch.W.E( load("m1_Nrtads_ind.txt") )
if(!is.null(attr(tC$value, "class"))) {
  m1.Nrtads <- glmmTMB(round.fitness ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s + Nr.tads.c + Nr.tads.c:loc.s + Nr.tads.c:loc.s2 + Nr.tads.c:larea.s + Nr.tads.c:larea.s2 + Nr.tads.c:date.s + Nr.tads.c:date.s2 + Nr.tads.c:loc.s:larea.s + Nr.tads.c:loc.s:date.s + Nr.tads.c:larea.s:date.s
            + ( 1 | pool ) ,
            zi = ~ loc.s + loc.s2 + larea.s + larea.s2 + date.s + date.s2 + loc.s:larea.s + loc.s:date.s + larea.s:date.s + Nr.tads.c + Nr.tads.c:loc.s + Nr.tads.c:loc.s2 + Nr.tads.c:larea.s + Nr.tads.c:larea.s2 + Nr.tads.c:date.s + Nr.tads.c:date.s2 + Nr.tads.c:loc.s:larea.s + Nr.tads.c:loc.s:date.s + Nr.tads.c:larea.s:date.s
            + ( 1 | pool ) ,
            family = poisson,
            control = glmmTMBControl(parallel = 4, profile = TRUE, optimizer = optim, optArgs = list(method = "BFGS")),
            data = d)
  save(m1.Nrtads, file = "m1_Nrtads_ind.txt" )
  }
if(is.null(attr(tC$value, "class"))) { load("m1_Nrtads_ind.txt") }
#
#
# Table 2.
# Importance of random effects not tested for the causal models.
#
n.int           <- 3  # number of digits
#
s1              <- summary(m1.drying)
ss.drying       <- paste(s1$ngrps$cond, " pools and ", s1$nobs, " observations.", sep = "")
m.CI            <- confint(m1.drying)         # Mean annual proportion pools drying
m.CI[,c(1,2,3)] <- apply(m.CI[,c(1,2,3)], 2, format.p.val, ndigits = n.int)
tab.drying      <- data.frame(source = c("Intercept", "Loc", "Loc^2", "Size", "Size^2", "Date", "Date^2", "Loc:Size", "Loc:Date", "Size:Date", "Drying", "Dry:Loc", "Dry:Loc^2", "Dry:Size", "Dry:Size^2", "Dry:Date", "Dry:Date^2", "Dry:Loc:Size", "Dry:Loc:Date", "Dry:Size:Date"), logistic = paste(m.CI[c(22:28,30:32,29,33:41),3], " (", m.CI[c(22:28,30:32,29,33:41),1], ", ", m.CI[c(22:28,30:32,29,33:41),2], ")", sep = ""), Poisson = paste(m.CI[c(1:7,9:11,8,12:20),3], " (", m.CI[c(1:7,9:11,8,12:20),1], ", ", m.CI[c(1:7,9:11,8,12:20),2], ")", sep = ""), stringsAsFactors = FALSE)
#
s1              <- summary(m1.mean.wave)
ss.mean.wave    <- paste(s1$ngrps$cond, " pools and ", s1$nobs, " observations.", sep = "")
m.CI            <- confint(m1.mean.wave)      # Mean wave height for the year
m.CI[,c(1,2,3)] <- apply(m.CI[,c(1,2,3)], 2, format.p.val, ndigits = n.int)
tab.mean.wave   <- data.frame(source = c("Intercept", "Loc", "Loc^2", "Size", "Size^2", "Date", "Date^2", "Loc:Size", "Loc:Date", "Size:Date", "Meanwave", "Meanwave:Loc", "Meanwave:Loc^2", "Meanwave:Size", "Meanwave:Size^2", "Meanwave:Date", "Meanwave:Date^2", "Meanwave:Loc:Size", "Meanwave:Loc:Date", "Meanwave:Size:Date"), logistic = paste(m.CI[c(22:28,30:32,29,33:41),3], " (", m.CI[c(22:28,30:32,29,33:41),1], ", ", m.CI[c(22:28,30:32,29,33:41),2], ")", sep = ""), Poisson = paste(m.CI[c(1:7,9:11,8,12:20),3], " (", m.CI[c(1:7,9:11,8,12:20),1], ", ", m.CI[c(1:7,9:11,8,12:20),2], ")", sep = ""), stringsAsFactors = FALSE)
#
s1              <- summary(m1.Aeshna)
ss.Aeshna       <- paste(s1$ngrps$cond, " pools and ", s1$nobs, " observations.", sep = "")
m.CI            <- confint(m1.Aeshna)         # Total annual Aeshna numbers
m.CI[,c(1,2,3)] <- apply(m.CI[,c(1,2,3)], 2, format.p.val, ndigits = n.int)
tab.aeshna      <- data.frame(source = c("Intercept", "Loc", "Loc^2", "Size", "Size^2", "Date", "Date^2", "Loc:Size", "Loc:Date", "Size:Date", "Aeshna", "Aeshna:Loc", "Aeshna:Loc^2", "Aeshna:Size", "Aeshna:Size^2", "Aeshna:Date", "Aeshna:Date^2", "Aeshna:Loc:Size", "Aeshna:Loc:Date", "Aeshna:Size:Date"), logistic = paste(m.CI[c(22:28,30:32,29,33:41),3], " (", m.CI[c(22:28,30:32,29,33:41),1], ", ", m.CI[c(22:28,30:32,29,33:41),2], ")", sep = ""), Poisson = paste(m.CI[c(1:7,9:11,8,12:20),3], " (", m.CI[c(1:7,9:11,8,12:20),1], ", ", m.CI[c(1:7,9:11,8,12:20),2], ")", sep = ""), stringsAsFactors = FALSE)
#
s1              <- summary(m1.Nrtads)
ss.density      <- paste(s1$ngrps$cond, " pools and ", s1$nobs, " observations.", sep = "")
m.CI            <- confint(m1.Nrtads)        # Total number of tadpoles on the study area
m.CI[,c(1,2,3)] <- apply(m.CI[,c(1,2,3)], 2, format.p.val, ndigits = n.int)
tab.dens        <- data.frame(source = c("Intercept", "Loc", "Loc^2", "Size", "Size^2", "Date", "Date^2", "Loc:Size", "Loc:Date", "Size:Date", "nr.tads", "nr.tads:Loc", "nr.tads:Loc^2", "nr.tads:Size", "nr.tads:Size^2", "nr.tads:Date", "nr.tads:Date^2", "nr.tads:Loc:Size", "nr.tads:Loc:Date", "nr.tads:Size:Date"), logistic = paste(m.CI[c(22:28,30:32,29,33:41),3], " (", m.CI[c(22:28,30:32,29,33:41),1], ", ", m.CI[c(22:28,30:32,29,33:41),2], ")", sep = ""), Poisson = paste(m.CI[c(1:7,9:11,8,12:20),3], " (", m.CI[c(1:7,9:11,8,12:20),1], ", ", m.CI[c(1:7,9:11,8,12:20),2], ")", sep = ""), stringsAsFactors = FALSE)
#
ss.drying; tab.drying
ss.mean.wave; tab.mean.wave
ss.Aeshna; tab.aeshna
ss.density; tab.dens
#
# Rearranging for Table 2 in the MS.
data.frame(Model.part = c("Logistic", tab.drying$source, "Poisson", tab.drying$source), Drying = c("", tab.drying$logistic, "", tab.drying$Poisson), Waves = c("", tab.mean.wave$logistic, "", tab.mean.wave$Poisson), Aeshna = c("", tab.aeshna$logistic, "", tab.aeshna$Poisson), PopSize = c("", tab.dens$logistic, "", tab.dens$Poisson) )
#
#
# Figure 4 -- univariate selection surfaces estimated by the full model.
#
years         <- unique(fc3$year)
years         <- years[order(years)]
range.loc     <- range(fc3$loc, na.rm = TRUE)
range.larea   <- range(fc3$larea, na.rm = TRUE)
range.area    <- exp(range.larea)
range.date    <- range(fc3$date, na.rm = TRUE)
range.loc.s   <- range(fc3$loc.s)
range.larea.s <- range(fc3$larea.s)
range.date.s  <- range(fc3$date.s)
indiv.points  <- aggregate(fc3[ , c("loc", "larea", "date", "loc.s", "larea.s", "date.s")], by = list(year = fc3$year, pool = fc3$pool, brood = fc3$brood), mean)
indiv.points  <- indiv.points[complete.cases(indiv.points),]
#
# Predict the overall shape of selection from the glmmTMB model.
# Population-level predictions. Takes 70 sec.
fine.grid          <- 40
pred.loc           <- data.frame(loc.s = set.up.interval(range.loc.s, fine.grid), larea.s = 0, larea.s2 = 0, date.s = 0, date.s2 = 0, year = NA, pool = NA)
pred.loc$loc.s2    <- pred.loc$loc.s^2 / sd.loc.s2.indivs.kept
pred.loc           <- cbind(pred.loc, Estimate = predict(m1, type = "response", re.form = NA, newdata = pred.loc) / 10)
pred.area          <- data.frame(larea.s = set.up.interval(range.larea, fine.grid), loc.s = 0, loc.s2 = 0, date.s = 0, date.s2 = 0, year = NA, pool = NA)
pred.area$larea.s2 <- pred.area$larea.s^2 / sd.larea.s2.indivs.kept
pred.area          <- cbind(pred.area, Estimate = predict(m1, type = "response", re.form = NA, newdata = pred.area) / 10)
pred.date          <- data.frame(date.s = set.up.interval(range.date.s, fine.grid), larea.s = 0, larea.s2 = 0, loc.s = 0, loc.s2 = 0, year = NA, pool = NA)
pred.date$date.s2  <- pred.date$date.s^2 / sd.date.s2.indivs.kept
pred.date$date     <- (pred.date$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
pred.date          <- cbind(pred.date, Estimate = predict(m1, type = "response", re.form = NA, newdata = pred.date) / 10)
#
#
# Predict separate selection surfaces for each year. Takes 40 sec.
start.time <- Sys.time()
fine.grid  <- 17
newdat     <- data.frame(year = integer(), loc.s = numeric(), larea.s = numeric(), date.s = numeric())
for(i in years) {
  d.sub1       <- subset(fc3, year == i)
  r.loc.s      <- range(d.sub1$loc.s, na.rm = TRUE)
  r.larea.s    <- range(d.sub1$larea.s, na.rm = TRUE)
  r.date.s     <- range(d.sub1$date.s, na.rm = TRUE)
  if (r.date.s[1] > 0) r.date.s[1] <- 0
  intval.loc   <- set.up.interval(r.loc.s, fine.grid)
  intval.larea <- set.up.interval(r.larea.s, fine.grid)
  intval.date  <- set.up.interval(r.date.s, fine.grid)
  newdat       <- rbind(newdat, expand.grid(year = i, loc.s = intval.loc, larea.s = intval.larea, date.s = intval.date) )
  }
newdat$loc.s2   <- newdat$loc.s^2 / sd.loc.s2.indivs.kept
newdat$larea.s2 <- newdat$larea.s^2 / sd.larea.s2.indivs.kept
newdat$date.s2  <- newdat$date.s^2 / sd.date.s2.indivs.kept
newdat$pool     <- NA
newdat          <- cbind(newdat, Estimate = predict(m1, type = "response", re.form = NULL, newdata = newdat) / 10)
tick.label.size <- 0.9
axis.label.size <- 1.1
year.list       <- unique(newdat$year)
setwd(figs.wd)
pdf("Fig 4.pdf", useDingbats = FALSE, width = 6, height = 8)
plot.new()
#
# Panel A: Location.
par(new = "TRUE", plt = c(0.23, 0.85, 0.7, 0.92))
  predict.loc     <- subset(newdat, newdat$larea.s == 0 & newdat$date.s == 0 )
  predict.loc$loc <- predict.loc$loc.s * sd.loc.indivs.kept + mean.loc.indivs.kept
  y.limits        <- c(0, 3.1)
  y.int           <- y.limits[2] / 6
  y.limits[1]     <- - y.int
  count           <- 0
  for (i in dens4$year) {
    count <- count + 1
    d.sub <- subset(predict.loc, year == i)
    if (count == 1) {
      plot(d.sub$loc, d.sub$Estimate, xlim = range.loc, pch = "", ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt = "n", yaxt = "n", yaxs = "i" )
      box(lwd = 1.2)
      axis(side = 1, tck = -0.025, labels = FALSE)
      axis(side = 2, tck = -0.025, labels = FALSE)
      mtext(text = c("0.3", "0.5", "0.7", "0.9"), side = 1, line = 0.3, las = 1, at = c(0.3, 0.5, 0.7, 0.9), cex = tick.label.size)
      mtext(text = "Pool location", side = 1, line = 1.2, las = 1, cex = axis.label.size)
      mtext(text = c("0", "1", "2", "3"), side = 2, line = 0.5, las = 1, at = c(0,1,2,3), cex = tick.label.size)
      mtext(text = "Relative fitness", side = 2, line = 1.3, cex = axis.label.size)
      }
    lines(d.sub$loc, d.sub$Estimate, col = main.color, lwd = 2)
    }
  # Ticks along the bottom illustrating observations.
  for(j in 1:dim(indiv.points)[1]) { lines(x = c(indiv.points$loc[j], indiv.points$loc[j]), y = c(y.limits[1], y.limits[1]/2))  }
  # Overlay the overall estimated fitness surface.
  lines(pred.loc$loc.s * sd.loc.indivs.kept + mean.loc.indivs.kept, pred.loc$Estimate, col = "black", lwd = 3.5)
  box(lwd=1.2)
  text("A", x = range.loc[1] + (0.985*(range.loc[2]-range.loc[1])), y = y.limits[1] + 0.91*(y.limits[2]-y.limits[1]), cex = 1.5)
  # Arrows at the bottom
  par(xpd = TRUE)
  mtext(text = c("Lake", "Forest"), side = 1, line = 1.15, at = c(range.loc[1]+0.09, range.loc[2]-0.08), adj = c(0,1), cex = tick.label.size, las = 0)
  arrows(x0 = range.loc[1]+0.08, y0 = -1.2, x1 = range.loc[1], y1 = -1.2, angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  arrows(x0 = range.loc[2]-0.07, y0 = -1.2, x1 = range.loc[2]+0.01, y1 = -1.2, angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  par(xpd = FALSE)
  mtext(paste("Van Buskirk & Smith, Figure 4,", Sys.Date()), side = 3, line = 2, cex = 0.7)
#
# Panel B: Pool size.
par(new = "TRUE", plt = c(0.23, 0.85, 0.4, 0.62))
  predict.area       <- subset(newdat, loc.s == 0 & date.s == 0)
  predict.area$larea <- predict.area$larea.s * sd.larea.indivs.kept + mean.larea.indivs.kept
  y.limits      <- range(predict.area$Estimate)
  y.limits[2]   <- 1.05 * 3.1 # y.limits[2]
  y.int         <- y.limits[2] / 6
  y.limits[1]   <- - y.int
  count         <- 0
  for (i in dens4$year) {
    count <- count + 1
    d.sub <- subset(predict.area, year == i)
    if (count == 1) {
      plot(d.sub$larea, d.sub$Estimate, xlim = range.larea, pch = "", ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt = "n", yaxt = "n", yaxs = "i" )
      box(lwd = 1.2)
      axis(side = 1, tck = -0.025, at = log(c(0.1, 1, 10, 20)), labels = FALSE)
      axis(side = 1, tck = -0.015, at = log(c(0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 2, 3, 4, 5, 6, 7, 8, 9)), labels = FALSE)
      axis(side = 2, tck = -0.025, labels = FALSE)
      mtext(text = c("0.1", "1", "10"), side = 1, line = 0.3, las = 1, at = log(c(0.1, 1, 10)), cex = tick.label.size)
      mtext(text = "Pool size (m  )", side = 1, line = 1.2, las = 1, cex = axis.label.size)
      mtext(text = "2", side = 1, line = 0.9, at = 0.613, las = 1, cex = axis.label.size-0.25)
      mtext(text = c("0", "1", "2", "3"), side = 2, line = 0.5, las = 1, at = c(0,1,2,3), cex = tick.label.size)
      mtext(text = "Relative fitness", side = 2, line = 1.3, cex = axis.label.size)
      }
    lines(d.sub$larea, d.sub$Estimate, col = main.color, lwd = 2)
    }
  for(j in 1:dim(indiv.points)[1]) { lines(x = c(indiv.points$larea[j], indiv.points$larea[j]), y = c(y.limits[1], y.limits[1]/2))  }
  lines(pred.area$larea.s, pred.area$Estimate, col = "black", lwd = 3.5)
  box(lwd=1.2)
  text("B", x = range.larea[1] + (0.985*(range.larea[2]-range.larea[1])), y = y.limits[1] + 0.91*(y.limits[2]-y.limits[1]), cex = 1.5)
#
# Panel C: Date of hatching.
par(new = "TRUE", plt = c(0.23, 0.85, 0.1, 0.32))
  predict.date      <- subset(newdat, loc.s == 0 & larea.s == 0)
  predict.date$date <- predict.date$date.s * sd.date.indivs.kept + mean.date.indivs.kept
  y.limits            <- range(predict.date$Estimate)
  y.limits[2]         <- 3.1
  y.int               <- y.limits[2] / 6
  y.limits[1]         <- - y.int
  count               <- 0
  for (i in dens4$year) {
    count <- count + 1
    d.sub <- subset(predict.date, year == i)
    if (count == 1) {
      plot(d.sub$date, d.sub$Estimate, xlim = range.date, pch = "", ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt = "n", yaxt = "n", yaxs = "i" )
      box(lwd = 1.2)
      axis(side = 1, tck = -0.025, labels = FALSE)
      axis(side = 2, tck = -0.025, labels = FALSE)
      mtext(text = c(paste(make.dates(30)[[1]], make.dates(30)[[2]], sep = " "), paste(make.dates(50)[[1]], make.dates(50)[[2]], sep = " "), paste(make.dates(70)[[1]], make.dates(70)[[2]], sep = " ")), side = 1, line = 0.3, las = 1, at = c(30,50,70), cex = tick.label.size)
      mtext(text = "Hatching date", side = 1, line = 1.2, las = 1, cex = axis.label.size)
      mtext(text = c("0", "1", "2", "3"), side = 2, line = 0.5, las = 1, at = c(0,1,2,3), cex = tick.label.size)
      mtext(text = "Relative fitness", side = 2, line = 1.3, cex = axis.label.size)
      for(j in 1:dim(indiv.points)[1]) { lines(x = c(indiv.points$date[j], indiv.points$date[j]), y = c(y.limits[1], y.limits[1]/2)) }
      }
    lines(d.sub$date, d.sub$Estimate, col = main.color, lwd = 2)
    }
  lines(pred.date$date, pred.date$Estimate, col = "black", lwd = 3.5)
  box(lwd=1.2)
  text("C", x = range.date[1] + (0.985*(range.date[2]-range.date[1])), y = y.limits[1] + 0.91*(y.limits[2]-y.limits[1]), cex = 1.5)
dev.off()
#
#
# Figure 6 -- How does selection on pool location change with ecological agents of selection?
#
fine.grid  <- 18
xlimits    <- c(0.2, 1)
xlim.range <- (xlimits - mean.loc.indivs.kept)/sd.loc.indivs.kept
#
# Drying
d.dry              <- d[ , c("propdry.c", "prop.dry", "loc", "loc.s", "loc.s2", "larea.s", "larea.s2", "date.s", "date.s2", "pool")]
d.dry              <- d.dry[complete.cases(d.dry), ]
d.dry              <- aggregate(d.dry[,c("prop.dry", "loc")], by = list(loc.s = d.dry$loc.s, propdry.c = d.dry$propdry.c), mean, na.rm = TRUE)
newdat             <- expand.grid(pool = NA, propdry.c = set.up.interval(c(-0.05, 0.05) + range(d.dry$propdry.c), fine.grid), loc.s = set.up.interval(xlim.range, fine.grid), larea.s = 0, larea.s2 = 0, date.s = 0, date.s2 = 0)
newdat$loc.s2      <- newdat$loc.s^2 / sd.loc.s2.indivs.kept
newdat.dry         <- cbind(newdat, Estimate = predict(m1.drying, type = "response", re.form = NA, newdata = newdat)/10 )
newdat.dry$loc     <- (newdat.dry$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.dry$propdry <- sd.propdry.kept * newdat.dry$propdry.c + mean.propdry.kept
newdat.dry         <- newdat.dry[,c("loc", "propdry", "Estimate")]
#
# Mean wave height
d.wv1              <- d[ , c("meanwave.c", "mean.wave.ht", "round.fitness", "loc", "loc.s", "loc.s2", "larea.s", "larea.s2", "date.s", "date.s2", "pool")]
d.wv1              <- d.wv1[ complete.cases(d.wv1), ]
d.wv1              <- aggregate(d.wv1[,c("mean.wave.ht", "loc")], by = list(loc.s = d.wv1$loc.s, meanwave.c = d.wv1$meanwave.c), mean, na.rm = TRUE)
d.wv1              <- rename(d.wv1, "mean.wave.ht", "meanwave")
newdat             <- expand.grid(pool = NA, meanwave.c = set.up.interval(c(-0.08,0.08) + range(d.wv1$meanwave.c), fine.grid), loc.s = set.up.interval(xlim.range, fine.grid), larea.s = 0, larea.s2 = 0, date.s = 0, date.s2 = 0)
newdat$loc.s2      <- newdat$loc.s^2 / sd.loc.s2.indivs.kept
newdat.wv1         <- cbind(newdat, Estimate = predict(m1.mean.wave, type = "response", re.form = NA, newdata = newdat)/10)
newdat.wv1$loc     <- (newdat.wv1$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.wv1$meanwave<- sd.meanwave.kept * newdat.wv1$meanwave.c + mean.meanwave.kept
newdat.wv1         <- newdat.wv1[,c("loc", "meanwave", "Estimate")]
#
# Aeshna
d.dfly             <- d[ , c("aeshna.c", "total.aeshna", "round.fitness", "loc", "loc.s", "loc.s2", "larea.s", "larea.s2", "date.s", "date.s2", "pool")]
d.dfly             <- d.dfly[ complete.cases(d.dfly), ]
d.dfly             <- aggregate(d.dfly[,c("total.aeshna", "loc")], by = list(loc.s = d.dfly$loc.s, aeshna.c = d.dfly$aeshna.c), mean, na.rm = TRUE)
d.dfly             <- rename(d.dfly, "total.aeshna", "aeshna")
newdat             <- expand.grid(pool = NA, aeshna.c = set.up.interval(c(-0.10,0.09) + range(d.dfly$aeshna.c), fine.grid), loc.s = set.up.interval(xlim.range, fine.grid), larea.s = 0, larea.s2 = 0, date.s = 0, date.s2 = 0)
newdat$loc.s2      <- newdat$loc.s^2 / sd.loc.s2.indivs.kept
newdat.dfly        <- cbind(newdat, Estimate = predict(m1.Aeshna, type = "response", re.form = NA, newdata = newdat)/10)
newdat.dfly$loc    <- (newdat.dfly$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.dfly$aeshna <- exp(newdat.dfly$aeshna.c * sd.aeshna.kept + mean.aeshna.kept)
newdat.dfly        <- newdat.dfly[,c("loc", "aeshna", "Estimate")]
#
# Competition (number of tadpoles)
d.dens             <- d[ , c("Nr.tads.c", "Nr.tads", "round.fitness", "loc", "loc.s", "loc.s2", "larea.s", "larea.s2", "date.s", "date.s2", "pool")]
d.dens             <- d.dens[ complete.cases(d.dens), ]
d.dens             <- aggregate(d.dens[,c("Nr.tads", "loc")], by = list(loc.s = d.dens$loc.s, Nr.tads.c = d.dens$Nr.tads.c), mean, na.rm = TRUE)
d.dens             <- rename(d.dens, "Nr.tads", "dens")
newdat             <- expand.grid(pool = NA, Nr.tads.c = set.up.interval(c(-0.1,0.1) + range(d.dens$Nr.tads.c), fine.grid), loc.s = set.up.interval(xlim.range, fine.grid), larea.s = 0, larea.s2 = 0, date.s = 0, date.s2 = 0)
newdat$loc.s2      <- newdat$loc.s^2 / sd.loc.s2.indivs.kept
newdat.dens        <- cbind(newdat, Estimate = predict(m1.Nrtads, type = "response", re.form = NA, newdata = newdat)/10)
newdat.dens$loc    <- (newdat.dens$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.dens$dens   <- exp(newdat.dens$Nr.tads.c * sd.Nrtads.kept + mean.Nrtads.kept)
newdat.dens        <- newdat.dens[,c("loc", "dens", "Estimate")]
#
newdat.dry[newdat.dry$Estimate > 1.9, "Estimate"] <- 1.9    # limit extreme predictionsd at the edges.
newdat.wv1[newdat.wv1$Estimate > 1.9, "Estimate"] <- 1.9
newdat.dfly[newdat.dfly$Estimate > 1.9, "Estimate"] <- 1.9
newdat.dens[newdat.dens$Estimate > 1.9, "Estimate"] <- 1.9
#
contour.label.size <- 0.85
tick.label.size    <- 0.9
axis.label.size    <- 1.1
small.point.size   <- 0.6
color.gradient     <- function(x, colors = c(lower.color, lowmid.color, highmid.color, upper.color), colsteps = 100) { return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )  }
color.range        <- range(newdat.dry$Estimate, newdat.wv1$Estimate, newdat.dfly$Estimate, newdat.dens$Estimate)
plot.colors        <- color.gradient(1:100)
#
pdf("Fig 6.pdf", useDingbats = FALSE, width = 8.5, height = 8)
plot.new()
  # Panel A: Drying
par(new = "TRUE", plt = c(0.10, 0.43, 0.54, 0.89))
  ylimits <- range(newdat.dry$propdry)
  filled.contour3(x = unique(newdat.dry$loc), y = unique(newdat.dry$propdry), z = matrix(newdat.dry$Estimate, nrow = length(unique(newdat.dry$propdry)), ncol = length(unique(newdat.dry$loc)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(newdat.dry$loc), y = unique(newdat.dry$propdry), z = matrix(newdat.dry$Estimate, nrow = length(unique(newdat.dry$propdry)), ncol = length(unique(newdat.dry$loc)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  points(d.dry$loc, d.dry$prop.dry, cex = small.point.size, col = "black", bg = "red", lwd = 0.8, xaxt = "n", yaxt = "n", pch = 21)
  box(lwd = 1.2)
  draw.contour.labels(xlimits = xlimits, ylimits, xfraction = 0.09, yfraction = 0.06, x.corners = c(0.35, 0.59, 0.525), y.corners = c(0.047, 0.128, 0.172), contour.labels = c("0.8", "1.0", "1.2"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  axis(side = 1, tck = -0.025, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 1, line = 0.3, at = c(0.2,0.4,0.6,0.8,1), las = 1, cex = tick.label.size)
  axis(side = 2, tck = -0.015, at = c(0.02,0.06,0.1,0.14,0.18), labels = FALSE)
  mtext(text = c("0.02", "0.06", "0.10", "0.14", "0.18"), side = 2, line = 0.5, las = 1, at = c(0.02,0.06,0.1,0.14,0.18), cex = tick.label.size)
  mtext(text = "Proportion of pools dried", side = 2, line = 2.1, cex = axis.label.size)
  mtext(text = "A. Pool drying", side = 3, line = 0.25, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size+0.1)
  mtext(text = paste("Van Buskirk & Smith, Figure 6,", Sys.Date()), side = 3, line = 3, las = 1, cex = 0.8)
  # Panel B: Wave wash
par(new = "TRUE", plt = c(0.54, 0.87, 0.54, 0.89))
  ylimits <- range(newdat.wv1$meanwave)
  filled.contour3(x = unique(newdat.wv1$loc), y = unique(newdat.wv1$meanwave), z = matrix(newdat.wv1$Estimate, nrow = length(unique(newdat.wv1$meanwave)), ncol = length(unique(newdat.wv1$loc)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(newdat.wv1$loc), y = unique(newdat.wv1$meanwave), z = matrix(newdat.wv1$Estimate, nrow = length(unique(newdat.wv1$meanwave)), ncol = length(unique(newdat.wv1$loc)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  polygon(x = c(0.77, xlimits[2], xlimits[2]), y = c(ylimits[2], ylimits[2], 0.157), border = NA, col = "white")
  points(d.wv1$loc, d.wv1$meanwave, cex = small.point.size, col = "black", bg = "red", lwd = 0.8, xaxt = "n", yaxt = "n", pch = 21)
  box(lwd = 1.2)
  draw.contour.labels(xlimits = xlimits, ylimits, xfraction = 0.09, yfraction = 0.06, x.corners = c(0.447, 0.735, 0.573), y.corners = c(0.175, 0.116, 0.077), contour.labels = c("0.6", "0.8", "1.2"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  axis(side = 1, tck = -0.025, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 1, line = 0.3, at = c(0.2,0.4,0.6,0.8,1), las = 1, cex = tick.label.size)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.08", "0.12", "0.16", "0.20"), side = 2, line = 0.5, las = 1, at = c(0.08,0.12,0.16,0.2), cex = tick.label.size)
  mtext(text = "Mean wave height", side = 2, line = 2.1, cex = axis.label.size)
  mtext(text = "B. Wave wash", side = 3, line = 0.25, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size+0.1)
  # Panel C: Dragonfly density
par(new = "TRUE", plt = c(0.10, 0.43, 0.105, 0.455))
  ylimits    <- log10(range(newdat.dfly$aeshna))
  ylimits[2] <- log10(500)
  filled.contour3(x = unique(newdat.dfly$loc), y = unique(log10(newdat.dfly$aeshna)), z = matrix(newdat.dfly$Estimate, nrow = length(unique(newdat.dfly$aeshna)), ncol = length(unique(newdat.dfly$loc)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(newdat.dfly$loc), y = unique(log10(newdat.dfly$aeshna)), z = matrix(newdat.dfly$Estimate, nrow = length(unique(newdat.dfly$aeshna)), ncol = length(unique(newdat.dfly$loc)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  points(d.dfly$loc, log10(d.dfly$aeshna), cex = small.point.size, col = "black", bg = "red", lwd = 0.8, xaxt = "n", yaxt = "n", pch = 21)
  box(lwd = 1.2)
  draw.contour.labels(xlimits = xlimits, ylimits, xfraction = 0.09, yfraction = 0.06, x.corners = c(0.795, 0.613, 0.516), y.corners = log10(c(130, 130, 340)), contour.labels = c("0.6", "1.0", "1.2"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  axis(side = 1, tck = -0.025, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 1, line = 0.3, at = c(0.2,0.4,0.6,0.8,1), las = 1, cex = tick.label.size)
  mtext(text = "Pool location", side = 1, line = 1.3, at = 0.593, las = 1, cex = axis.label.size)
  axis(side = 2, tck = -0.022, at = log10(c(150,300,500)), labels = FALSE)
  axis(side = 2, tck = -0.012, at = log10(c(200,400)), labels = FALSE)
  mtext(text = c("150", "300", "500"), side = 2, line = 0.5, las = 1, at = log10(c(150,300,500)), cex = tick.label.size)
  mtext(text = "Number of dragonflies", side = 2, line = 2, cex = axis.label.size)
  par(xpd = TRUE)
  mtext(text = c("Lake", "Forest"), side = 1, line = 1.3, las = 1, at = c(xlimits[1]+0.09, xlimits[2]-0.09), adj = c(0,1), cex = tick.label.size-0.05, las = 0)
  arrows(x0 = xlimits[1]+0.08, y0 = log10(87.5), x1 = xlimits[1], y1 = log10(87.5), angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  arrows(x0 = xlimits[2]-0.08, y0 = log10(87.5), x1 = xlimits[2], y1 = log10(87.5), angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  par(xpd = FALSE)
  mtext(text = "C. Predation", side = 3, line = 0.25, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size+0.1)
  # Panel D: Mean tadpole density, lower right
par(new = "TRUE", plt = c(0.54, 0.87, 0.105, 0.455))
  ylimits <- log10(range(newdat.dens$dens))
  filled.contour3(x = unique(newdat.dens$loc), y = unique(log10(newdat.dens$dens)), z = matrix(newdat.dens$Estimate, nrow = length(unique(newdat.dens$dens)), ncol = length(unique(newdat.dens$loc)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(newdat.dens$loc), y = unique(log10(newdat.dens$dens)), z = matrix(newdat.dens$Estimate, nrow = length(unique(newdat.dens$dens)), ncol = length(unique(newdat.dens$loc)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  polygon(x = c(0.76, xlimits[2], xlimits[2]), y = c(ylimits[1], log10(1100), ylimits[1]), border = NA, col = "white")
  polygon(x = c(xlimits[1], xlimits[1], 0.36), y = c(log10(1500), ylimits[1], ylimits[1]), border = NA, col = "white")
  points(d.dens$loc, log10(d.dens$dens), cex = small.point.size, col = "black", bg = "red", lwd = 0.8, xaxt = "n", yaxt = "n", pch = 21)
  box(lwd = 1.2)
  draw.contour.labels(xlimits = xlimits, ylimits, xfraction = 0.09, yfraction = 0.06, x.corners = c(0.7525, 0.8, 0.705), y.corners = log10(c(4150, 1700, 930)), contour.labels = c("0.6", "1.0", "1.4"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  axis(side = 1, tck = -0.025, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 1, line = 0.3, at = c(0.2,0.4,0.6,0.8,1), las = 1, cex = tick.label.size)
  mtext(text = "Pool location", side = 1, line = 1.3, at = 0.593, las = 1, cex = axis.label.size)
  axis(side = 2, tck = -0.022, at = log10(c(500,1000,5000)), labels = FALSE)
  axis(side = 2, tck = -0.012, at = log10(c(400,600,700,800,900,2000,3000,4000,6000,7000,8000)), labels = FALSE)
  mtext(text = c("500", "1000", "5000"), side = 2, line = 0.5, las = 1, at = log10(c(500,1000,5000)), cex = tick.label.size)
  mtext(text = "Number of tadpoles", side = 2, line = 2.2, cex = axis.label.size)
  par(xpd = TRUE)
  mtext(text = c("Lake", "Forest"), side = 1, line = 1.3, las = 1, at = c(xlimits[1]+0.09, xlimits[2]-0.09), adj = c(0,1), cex = tick.label.size-0.05, las = 0)
  arrows(x0 = xlimits[1]+0.08, y0 = log10(307), x1 = xlimits[1], y1 = log10(307), angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  arrows(x0 = xlimits[2]-0.08, y0 = log10(307), x1 = xlimits[2], y1 = log10(307), angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  par(xpd = FALSE)
  mtext(text = "D. Competition", side = 3, line = 0.25, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size+0.1)
  #
  # Right-side scale.
  par(new = "TRUE", xpd = TRUE, plt = c(0.90, 0.925, 0.2975, 0.6975), las = 1, cex.axis = 1, tck = -0.25 )
  lower.limit <- color.range[1]
  upper.limit <- color.range[2]
  plot.breaks <- round( seq(lower.limit, upper.limit, length.out = 101), 2)
  image.scale(z = matrix(newdat.dens$Estimate, nrow = length(unique(newdat.dens$dens)), ncol = length(unique(newdat.dens$loc)), byrow = TRUE),
     zlim = c(lower.limit, upper.limit), col = plot.colors, xlab = "", ylab = "", cex.lab = 1.5, breaks = plot.breaks, horiz = FALSE, yaxt = "n", xaxt = "n")
    # Label the vertical axis.
  scale.labels <- c(0.5, 1.0, 1.5)
  axis(4, at = scale.labels, labels = FALSE)
  mtext(text = c("0.5", "1.0", "1.5"), side = 4, line = 0.45, las = 1, at = scale.labels, cex = 0.85)
  box(lwd = 1.2)
    # Label at the bottom.
  y.first.line <- lower.limit - 0.02
  text(labels = c("Relative", "fitness"), pos = 1, cex = 0.85, x = 0.6, y = c(y.first.line, y.first.line-0.1))
dev.off()
#
#
# Figure S5 -- How does selection on pool size change with the four ecological agents?
#
fine.grid  <- 18
xlimits    <- log(c(0.035, 18.5))
xlim.range <- (xlimits - mean.larea.indivs.kept)/sd.larea.indivs.kept
#
# Drying
d.dry              <- d[ , c("propdry.c", "prop.dry", "round.fitness", "loc.s", "loc.s2", "larea", "larea.s", "larea.s2", "date.s", "date.s2", "pool")]
d.dry              <- d.dry[ complete.cases(d.dry), ]
d.dry              <- aggregate(d.dry[,c("prop.dry", "larea")], by = list(larea.s = d.dry$larea.s, propdry.c = d.dry$propdry.c), mean, na.rm = TRUE)
d.dry              <- rename(d.dry, "prop.dry", "propdry")
d.dry$area         <- exp(d.dry$larea)
newdat             <- expand.grid(pool = NA, propdry.c = set.up.interval(c(-0.05, 0.05) + range(d.dry$propdry.c), fine.grid), larea.s = set.up.interval(xlim.range, fine.grid), loc.s = 0, loc.s2 = 0, date.s = 0, date.s2 = 0)
newdat$larea.s2    <- newdat$larea.s^2 / sd.larea.s2.indivs.kept
newdat.dry         <- cbind(newdat, Estimate = predict(m1.drying, type = "response", re.form = NA, newdata = newdat)/10 )
newdat.dry$area    <- exp((newdat.dry$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
newdat.dry$propdry <- sd.propdry.kept * newdat.dry$propdry.c + mean.propdry.kept
newdat.dry         <- newdat.dry[,c("area", "propdry", "Estimate")]
#
# Mean wave height
d.wv1              <- d[ , c("meanwave.c", "mean.wave.ht", "round.fitness", "larea", "loc.s", "loc.s2", "larea.s", "larea.s2", "date.s", "date.s2", "pool")]
d.wv1              <- d.wv1[ complete.cases(d.wv1), ]
d.wv1              <- aggregate(d.wv1[,c("mean.wave.ht", "larea")], by = list(larea.s = d.wv1$larea.s, meanwave.c = d.wv1$meanwave.c), mean, na.rm = TRUE)
d.wv1              <- rename(d.wv1, "mean.wave.ht", "meanwave")
d.wv1$area         <- exp(d.wv1$larea)
newdat             <- expand.grid(pool = NA, meanwave.c = set.up.interval(c(-0.08,0.08) + range(d.wv1$meanwave.c), fine.grid), larea.s = set.up.interval(xlim.range, fine.grid), loc.s = 0, loc.s2 = 0, date.s = 0, date.s2 = 0)
newdat$larea.s2    <- newdat$larea.s^2 / sd.larea.s2.indivs.kept
newdat.wv1         <- cbind(newdat, Estimate = predict(m1.mean.wave, type = "response", re.form = NA, newdata = newdat)/10)
newdat.wv1$area    <- exp((newdat.wv1$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
newdat.wv1$meanwave<- sd.meanwave.kept * newdat.wv1$meanwave.c + mean.meanwave.kept
newdat.wv1         <- newdat.wv1[,c("area", "meanwave", "Estimate")]
#
# Aeshna
d.dfly             <- d[ , c("aeshna.c", "total.aeshna", "round.fitness", "larea", "loc.s", "loc.s2", "larea.s", "larea.s2", "date.s", "date.s2", "pool")]
d.dfly             <- d.dfly[ complete.cases(d.dfly), ]
d.dfly             <- aggregate(d.dfly[,c("total.aeshna", "larea")], by = list(larea.s = d.dfly$larea.s, aeshna.c = d.dfly$aeshna.c), mean, na.rm = TRUE)
d.dfly             <- rename(d.dfly, "total.aeshna", "aeshna")
d.dfly$area        <- exp(d.dfly$larea)
newdat             <- expand.grid(pool = NA, aeshna.c = set.up.interval(c(-0.10,0.09) + range(d.dfly$aeshna.c), fine.grid), larea.s = set.up.interval(xlim.range, fine.grid), loc.s = 0, loc.s2 = 0, date.s = 0, date.s2 = 0)
newdat$larea.s2    <- newdat$larea.s^2 / sd.larea.s2.indivs.kept
newdat.dfly        <- cbind(newdat, Estimate = predict(m1.Aeshna, type = "response", re.form = NA, newdata = newdat)/10)
newdat.dfly$area   <- exp((newdat.dfly$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
newdat.dfly$aeshna <- exp(newdat.dfly$aeshna.c * sd.aeshna.kept + mean.aeshna.kept)
newdat.dfly        <- newdat.dfly[,c("area", "aeshna", "Estimate")]
newdat.dfly[newdat.dfly$Estimate > 1.6, "Estimate"] <- 1.6
#
# Competition (number of tadpoles)
d.dens             <- d[ , c("Nr.tads.c", "Nr.tads", "round.fitness", "loc.s", "loc.s2", "larea", "larea.s", "larea.s2", "date.s", "date.s2", "pool")]
d.dens             <- d.dens[ complete.cases(d.dens), ]
d.dens             <- aggregate(d.dens[,c("Nr.tads", "larea")], by = list(larea.s = d.dens$larea.s, Nr.tads.c.c = d.dens$Nr.tads.c), mean, na.rm = TRUE)
d.dens             <- rename(d.dens, "Nr.tads", "dens")
d.dens$area        <- exp(d.dens$larea)
newdat             <- expand.grid(pool = NA, Nr.tads.c = set.up.interval(c(-0.1,0.1) + range(d.dens$Nr.tads.c), fine.grid), larea.s = set.up.interval(xlim.range, fine.grid), loc.s = 0, loc.s2 = 0, date.s = 0, date.s2 = 0)
newdat$larea.s2    <- newdat$larea.s^2 / sd.larea.s2.indivs.kept
newdat.dens        <- cbind(newdat, Estimate = predict(m1.Nrtads, type = "response", re.form = NA, newdata = newdat)/10)
newdat.dens$area   <- exp((newdat.dens$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
newdat.dens$dens   <- exp(newdat.dens$Nr.tads.c * sd.Nrtads.kept + mean.Nrtads.kept)
newdat.dens        <- newdat.dens[,c("area", "dens", "Estimate")]
newdat.dens[newdat.dens$Estimate > 20, "Estimate"] <- 20
#
#
xlimits            <- c(0.035, 18.5)
contour.label.size <- 0.85
point.size         <- 0.5
tick.label.size    <- 0.9
axis.label.size    <- 1.1
color.gradient     <- function(x, colors = c(lower.color, lowmid.color, highmid.color, upper.color), colsteps = 100) { return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )  }
plot.colors        <- color.gradient(1:100)
color.range        <- range(newdat.dry$Estimate, newdat.wv1$Estimate, newdat.dfly$Estimate, newdat.dens$Estimate)
#
setwd(figs.wd)
pdf("Fig S5.pdf", useDingbats = FALSE, width = 8.5, height = 8)
plot.new()
  # Panel A: Drying
par(new = "TRUE", plt = c(0.10, 0.43, 0.54, 0.89))
  ylimits <- range(newdat.dry$propdry)
  filled.contour3(x = unique(log10(newdat.dry$area)), y = unique(newdat.dry$propdry), z = matrix(newdat.dry$Estimate, nrow = length(unique(newdat.dry$propdry)), ncol = length(unique(newdat.dry$area)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = log10(xlimits), xaxs="i")
  contour(x = unique(log10(newdat.dry$area)), y = unique(newdat.dry$propdry), z = matrix(newdat.dry$Estimate, nrow = length(unique(newdat.dry$propdry)), ncol = length(unique(newdat.dry$area)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  polygon(x = c(log10(8), log10(xlimits[2]), log10(xlimits[2])), y = c(ylimits[1], 0.14, ylimits[1]), border = NA, col = "white")
  draw.contour.labels(xlimits = log10(xlimits), ylimits, xfraction = 0.09, yfraction = 0.06, x.corners = log10(c(0.12, 1.87, 4.2)), y.corners = c(0.047, 0.047, 0.13), contour.labels = c("0.6", "0.8", "1.0"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(log10(d.dry$area), d.dry$propdry, cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.025, at = log10(c(0.1,1,10)), labels = FALSE)
  axis(side = 1, tck = -0.01, at = log10(c(0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9)), labels = FALSE)
  mtext(text = c("0.1", "1.0", "10"), side = 1, line = 0.3, at = log10(c(0.1,1,10)), las = 1, cex = tick.label.size)
  axis(side = 2, tck = -0.015, at = c(0.02,0.06,0.1,0.14,0.18), labels = FALSE)
  mtext(text = c("0.02", "0.06", "0.10", "0.14", "0.18"), side = 2, line = 0.5, las = 1, at = c(0.02,0.06,0.1,0.14,0.18), cex = tick.label.size)
  mtext(text = "Proportion of pools dried", side = 2, line = 2.1, cex = axis.label.size)
  mtext(text = "A. Pool drying", side = 3, line = 0.25, las = 1, at = log10(xlimits[1]), adj = 0, cex = axis.label.size+0.1)
  mtext(text = paste("Van Buskirk & Smith, Figure S5,", Sys.Date()), side = 3, line = 3, las = 1, cex = 0.9)
  #
  # Panel B: Wave wash
par(new = "TRUE", plt = c(0.54, 0.87, 0.54, 0.89))
  ylimits <- range(newdat.wv1$meanwave)
  filled.contour3(x = unique(log10(newdat.wv1$area)), y = unique(newdat.wv1$meanwave), z = matrix(newdat.wv1$Estimate, nrow = length(unique(newdat.wv1$meanwave)), ncol = length(unique(newdat.wv1$area)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = log10(xlimits), xaxs="i")
  contour(x = unique(log10(newdat.wv1$area)), y = unique(newdat.wv1$meanwave), z = matrix(newdat.wv1$Estimate, nrow = length(unique(newdat.wv1$meanwave)), ncol = length(unique(newdat.wv1$area)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  polygon(x = log10(c(xlimits[1], xlimits[1], 0.043)), y = c(0.14, ylimits[2], ylimits[2]), border = NA, col = "white")
  polygon(x = c(log10(7), log10(xlimits[2]), log10(xlimits[2])), y = c(ylimits[2], ylimits[2], 0.12), border = NA, col = "white")
  draw.contour.labels(xlimits = log10(xlimits), ylimits, xfraction = 0.09, yfraction = 0.06, x.corners = log10(c(0.61, 0.083, 6.5)), y.corners = c(0.17, 0.078, 0.081), contour.labels = c("0.6", "0.8", "1.0"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(log10(d.wv1$area), d.wv1$meanwave, cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.025, at = log10(c(0.1,1,10)), labels = FALSE)
  axis(side = 1, tck = -0.01, at = log10(c(0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9)), labels = FALSE)
  mtext(text = c("0.1", "1.0", "10"), side = 1, line = 0.3, at = log10(c(0.1,1,10)), las = 1, cex = tick.label.size)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.08", "0.12", "0.16", "0.20"), side = 2, line = 0.5, las = 1, at = c(0.08,0.12,0.16,0.2), cex = tick.label.size)
  mtext(text = "Mean wave height", side = 2, line = 2.1, cex = axis.label.size)
  mtext(text = "B. Wave wash", side = 3, line = 0.25, las = 1, at = log10(xlimits[1]), adj = 0, cex = axis.label.size+0.1)
  #
  # Panel C: Dragonfly density
par(new = "TRUE", plt = c(0.10, 0.43, 0.105, 0.455))
  ylimits    <- log10(range(newdat.dfly$aeshna))
  ylimits[2] <- log10(500)
  filled.contour3(x = unique(log10(newdat.dfly$area)), y = unique(log10(newdat.dfly$aeshna)), z = matrix(newdat.dfly$Estimate, nrow = length(unique(newdat.dfly$aeshna)), ncol = length(unique(newdat.dfly$area)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = log10(xlimits), xaxs="i")
  contour(x = unique(log10(newdat.dfly$area)), y = unique(log10(newdat.dfly$aeshna)), z = matrix(newdat.dfly$Estimate, nrow = length(unique(newdat.dfly$aeshna)), ncol = length(unique(newdat.dfly$area)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  polygon(x = log10(c(xlimits[1], xlimits[1], 0.05)), y = c(log10(370), ylimits[2], ylimits[2]), border = NA, col = "white")
  polygon(x = c(log10(8), log10(xlimits[2]), log10(xlimits[2])), y = c(ylimits[2], ylimits[2], log10(320)), border = NA, col = "white")
  polygon(x = c(log10(8), log10(xlimits[2]), log10(xlimits[2])), y = c(ylimits[1], log10(160), ylimits[1]), border = NA, col = "white")
  draw.contour.labels(xlimits = log10(xlimits), ylimits, xfraction = 0.09, yfraction = 0.06, x.corners = log10(c(0.17, 1.45, 3.4)), y.corners = log10(c(130, 340, 130)), contour.labels = c("0.8", "1.2", "1.4"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(log10(d.dfly$area), log10(d.dfly$aeshna), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.025, at = log10(c(0.1,1,10)), labels = FALSE)
  axis(side = 1, tck = -0.01, at = log10(c(0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9)), labels = FALSE)
  mtext(text = c("0.1", "1.0", "10"), side = 1, line = 0.3, at = log10(c(0.1,1,10)), las = 1, cex = tick.label.size)
  mtext(text = "Pool size (m  )", side = 1, line = 1.3, las = 1, cex = axis.label.size)
  mtext(text = "2", side = 1, line = 1.0, at = log10(2.3), las = 1, cex = 0.85)
  axis(side = 2, tck = -0.025, at = log10(c(150,500)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log10(c(200,300,400)), labels = FALSE)
  mtext(text = c("150", "300", "500"), side = 2, line = 0.5, las = 1, at = log10(c(150,300,500)), cex = tick.label.size)
  mtext(text = "Number of dragonflies", side = 2, line = 2, cex = axis.label.size)
  mtext(text = "C. Predation", side = 3, line = 0.25, las = 1, at = log10(xlimits[1]), adj = 0, cex = axis.label.size+0.1)
  #
  # Panel D: Mean tadpole density, lower right
par(new = "TRUE", plt = c(0.54, 0.87, 0.105, 0.455))
  ylimits <- log10(range(newdat.dens$dens))
  filled.contour3(x = unique(log10(newdat.dens$area)), y = unique(log10(newdat.dens$dens)), z = matrix(newdat.dens$Estimate, nrow = length(unique(newdat.dens$dens)), ncol = length(unique(newdat.dens$area)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = log10(xlimits), xaxs="i")
  contour(x = unique(log10(newdat.dens$area)), y = unique(log10(newdat.dens$dens)), z = matrix(newdat.dens$Estimate, nrow = length(unique(newdat.dens$dens)), ncol = length(unique(newdat.dens$area)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  polygon(x = c(log10(xlimits[1]), log10(xlimits[1]), log10(0.07)), y = c(ylimits[2], log10(6100), ylimits[2]), border = NA, col = "white")
  polygon(x = c(log10(6), log10(xlimits[2]), log10(xlimits[2])), y = c(ylimits[1], log10(3000), ylimits[1]), border = NA, col = "white")
  polygon(x = c(log10(xlimits[1]), log10(xlimits[1]), log10(0.14)), y = c(ylimits[1], log10(1500), ylimits[1]), border = NA, col = "white")
  draw.contour.labels(xlimits = log10(xlimits), ylimits, xfraction = 0.09, yfraction = 0.06, x.corners = log10(c(0.10, 4.2, 2.8)), y.corners = log10(c(1600, 1600, 800)), contour.labels = c("0.5", "1.0", "1.5"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(log10(d.dens$area), log10(d.dens$dens), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.025, at = log10(c(0.1,1,10)), labels = FALSE)
  axis(side = 1, tck = -0.01, at = log10(c(0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9)), labels = FALSE)
  mtext(text = c("0.1", "1.0", "10"), side = 1, line = 0.3, at = log10(c(0.1,1,10)), las = 1, cex = tick.label.size)
  mtext(text = "Pool size (m  )", side = 1, line = 1.3, las = 1, cex = axis.label.size)
  mtext(text = "2", side = 1, line = 1.0, at = log10(2.3), las = 1, cex = 0.85)
  axis(side = 2, tck = -0.022, at = log10(c(500,1000,5000)), labels = FALSE)
  axis(side = 2, tck = -0.012, at = log10(c(400,600,700,800,900,2000,3000,4000,6000,7000,8000)), labels = FALSE)
  mtext(text = c("500", "1000", "5000"), side = 2, line = 0.5, las = 1, at = log10(c(500,1000,5000)), cex = tick.label.size)
  mtext(text = "Number of tadpoles", side = 2, line = 2.2, cex = axis.label.size)
  mtext(text = "D. Competition", side = 3, line = 0.25, las = 1, at = log10(xlimits[1]), adj = 0, cex = axis.label.size+0.1)
  #
  # Right-side scale.
  par(new = "TRUE", xpd = TRUE, plt = c(0.90, 0.925, 0.2975, 0.6975), las = 1, cex.axis = 1, tck = -0.25 )
  lower.limit <- color.range[1]
  upper.limit <- color.range[2]
  plot.breaks <- round( seq(lower.limit, upper.limit, length.out = 101), 2)
  image.scale(z = matrix(newdat.dens$Estimate, nrow = length(unique(newdat.dens$dens)), ncol = length(unique(newdat.dens$loc)), byrow = TRUE),
     zlim = c(lower.limit, upper.limit), col = plot.colors, xlab = "", ylab = "", cex.lab = 1.5, breaks = plot.breaks, horiz = FALSE, yaxt = "n", xaxt = "n")
    # Label the vertical axis.
  scale.labels <- c(0.5, 1, 1.5, 2)
  axis(4, at = scale.labels, labels = FALSE)
  mtext(text = c("0.5", "1.0", "1.5", "2.0"), side = 4, line = 0.45, las = 1, at = scale.labels, cex = 0.85)
  box(lwd = 1.2)
    # Label at the bottom.
  y.first.line <- lower.limit - 0.02
  text(labels = c("Relative", "fitness"), pos = 1, cex = 0.85, x = 0.6, y = c(y.first.line, y.first.line-0.12))
dev.off()
#
#
# Figure S6 -- How does selection on date-of-hatching change with the four ecological agents?
#
fine.grid  <- 18
xlimits    <- c(21.5, 84)
xlim.range <- (xlimits - mean.date.indivs.kept)/sd.date.indivs.kept
#
# Drying
d.dry              <- d[ , c("propdry.c", "prop.dry", "loc.s", "loc.s2", "larea.s", "larea.s2", "date", "date.s", "date.s2", "pool")]
d.dry              <- d.dry[ complete.cases(d.dry), ]
d.dry              <- aggregate(d.dry[,c("prop.dry", "date")], by = list(date.s = d.dry$date.s, propdry.c = d.dry$propdry.c), mean, na.rm = TRUE)
d.dry              <- rename(d.dry, "prop.dry", "propdry")
newdat             <- expand.grid(pool = NA, propdry.c = set.up.interval(c(-0.05, 0.05) + range(d.dry$propdry.c), fine.grid), date.s = set.up.interval(xlim.range, fine.grid), larea.s = 0, larea.s2 = 0, loc.s = 0, loc.s2 = 0)
newdat$date.s2     <- newdat$date.s^2 / sd.date.s2.indivs.kept
newdat.dry         <- cbind(newdat, Estimate = predict(m1.drying, type = "response", re.form = NA, newdata = newdat)/10 )
newdat.dry$date    <- (newdat.dry$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.dry$propdry <- sd.propdry.kept * newdat.dry$propdry.c + mean.propdry.kept
newdat.dry         <- newdat.dry[,c("date", "propdry", "Estimate")]
#
# Mean wave height
d.wv1              <- d[ , c("meanwave.c", "mean.wave.ht", "loc.s", "loc.s2", "larea.s", "larea.s2", "date", "date.s", "date.s2", "pool")]
d.wv1              <- d.wv1[ complete.cases(d.wv1), ]
d.wv1              <- aggregate(d.wv1[,c("mean.wave.ht", "date")], by = list(date.s = d.wv1$date.s, meanwave.c = d.wv1$meanwave.c), mean, na.rm = TRUE)
d.wv1              <- rename(d.wv1, "mean.wave.ht", "meanwave")
newdat             <- expand.grid(pool = NA, meanwave.c = set.up.interval(c(-0.08,0.08)+range(d.wv1$meanwave.c), fine.grid), loc.s = 0, loc.s2 = 0, larea.s = 0, larea.s2 = 0, date.s = set.up.interval(xlim.range, fine.grid) )
newdat$date.s2     <- newdat$date.s^2 / sd.date.s2.indivs.kept
newdat.wv1         <- cbind(newdat, Estimate = predict(m1.mean.wave, type = "response", re.form = NA, newdata = newdat)/10)
newdat.wv1$date    <- (newdat.wv1$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.wv1$meanwave<- newdat.wv1$meanwave.c * sd.meanwave.kept + mean.meanwave.kept
#
# Aeshna
d.dfly             <- d[ , c("aeshna.c", "total.aeshna", "loc.s", "loc.s2", "larea.s", "larea.s2", "date", "date.s", "date.s2", "pool")]
d.dfly             <- d.dfly[ complete.cases(d.dfly), ]
d.dfly             <- aggregate(d.dfly[,c("total.aeshna", "date")], by = list(date.s = d.dfly$date.s, aeshna.c = d.dfly$aeshna.c), mean, na.rm = TRUE)
d.dfly             <- rename(d.dfly, "total.aeshna", "aeshna")
newdat             <- expand.grid(pool = NA, aeshna.c = set.up.interval(c(-0.1,0.09)+range(d.dfly$aeshna.c), fine.grid), loc.s = 0, loc.s2 = 0, larea.s = 0, larea.s2 = 0, date.s = set.up.interval(xlim.range, fine.grid) )
newdat$date.s2     <- newdat$date.s^2 / sd.date.s2.indivs.kept
newdat.dfly        <- cbind(newdat, Estimate = predict(m1.Aeshna, type = "response", re.form = NA, newdata = newdat)/10)
newdat.dfly$date   <- (newdat.dfly$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.dfly$aeshna <- exp(newdat.dfly$aeshna.c * sd.aeshna.kept + mean.aeshna.kept)
#
# Competition (density of tadpoles)
d.dens             <- d[ , c("Nr.tads.c", "Nr.tads", "round.fitness", "loc.s", "loc.s2", "larea", "larea.s", "larea.s2", "date.s", "date.s2", "date", "pool")]
d.dens             <- d.dens[ complete.cases(d.dens), ]
d.dens             <- aggregate(d.dens[,c("Nr.tads", "date")], by = list(date.s = d.dens$date.s, Nr.tads.c.c = d.dens$Nr.tads.c), mean, na.rm = TRUE)
d.dens             <- rename(d.dens, "Nr.tads", "dens")
newdat             <- expand.grid(pool = NA, Nr.tads.c = set.up.interval(c(-0.1,0.1)+range(d.dens$Nr.tads.c), fine.grid), loc.s = 0, loc.s2 = 0, larea.s = 0, larea.s2 = 0, date.s = set.up.interval(xlim.range, fine.grid) )
newdat$date.s2     <- newdat$date.s^2 / sd.date.s2.indivs.kept
newdat.dens        <- cbind(newdat, Estimate = predict(m1.Nrtads, type = "response", re.form = NA , newdata = newdat)/10)
newdat.dens$date   <- (newdat.dens$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.dens$dens   <- exp(newdat.dens$Nr.tads.c * sd.Nrtads.kept + mean.Nrtads.kept)
newdat.dry[newdat.dry$Estimate > 2,   "Estimate"] <- 2
newdat.wv1[newdat.wv1$Estimate > 2,   "Estimate"] <- 2
newdat.dfly[newdat.dfly$Estimate > 2, "Estimate"] <- 2
newdat.dens[newdat.dens$Estimate > 2, "Estimate"] <- 2
#
xlimits            <- c(23, 76)
contour.label.size <- 0.85
point.size         <- 0.5
tick.label.size    <- 0.9
axis.label.size    <- 1.1
color.gradient     <- function(x, colors = c(lower.color, lowmid.color, highmid.color, upper.color), colsteps = 100) { return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )  }
plot.colors        <- color.gradient(1:100)
color.range        <- range(newdat.dry$Estimate, newdat.wv1$Estimate, newdat.dfly$Estimate, newdat.dens$Estimate)
#
setwd(figs.wd)
pdf("Fig S6.pdf", useDingbats = FALSE, width = 8.5, height = 8)
plot.new()
  # A. Drying
par(new = "TRUE", plt = c(0.10, 0.43, 0.54, 0.89))
  ylimits <- range(newdat.dry$propdry)
  filled.contour3(x = unique(newdat.dry$date), y = unique(newdat.dry$propdry), z = matrix(newdat.dry$Estimate, nrow = length(unique(newdat.dry$propdry)), ncol = length(unique(newdat.dry$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(newdat.dry$date), y = unique(newdat.dry$propdry), z = matrix(newdat.dry$Estimate, nrow = length(unique(newdat.dry$propdry)), ncol = length(unique(newdat.dry$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 7 )
  polygon(x = c(xlimits[1], xlimits[1], 30), y = c(0.06, ylimits[1], ylimits[1]), border = NA, col = "white")
  draw.contour.labels(xlimits = xlimits, ylimits, xfraction = 0.09, yfraction = 0.06, x.corners = c(52.8, 43.6, 33.5), y.corners = c(0.13, 0.13, 0.13), contour.labels = c("0.5", "1.0", "1.5"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(d.dry$date, d.dry$propdry, cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c(paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ")), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.015, at = c(0.02,0.06,0.1,0.14,0.18), labels = FALSE)
  mtext(text = c("0.02", "0.06", "0.10", "0.14", "0.18"), side = 2, line = 0.5, las = 1, at = c(0.02,0.06,0.1,0.14,0.18), cex = tick.label.size)
  mtext(text = "Proportion of pools dried", side = 2, line = 2.1, cex = axis.label.size)
  mtext(text = "A. Pool drying", side = 3, line = 0.25, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size+0.1)
  mtext(text = paste("Van Buskirk & Smith, Figure S6,", Sys.Date()), side = 3, line = 3, las = 1, cex = 0.8)
  #
  # B. Wave wash
par(new = "TRUE", plt = c(0.54, 0.87, 0.54, 0.89))
  ylimits <- range(newdat.wv1$meanwave)
  filled.contour3(x = unique(newdat.wv1$date), y = unique(newdat.wv1$meanwave), z = matrix(newdat.wv1$Estimate, nrow = length(unique(newdat.wv1$meanwave)), ncol = length(unique(newdat.wv1$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(newdat.wv1$date), y = unique(newdat.wv1$meanwave), z = matrix(newdat.wv1$Estimate, nrow = length(unique(newdat.wv1$meanwave)), ncol = length(unique(newdat.wv1$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  polygon(x = c(xlimits[1], xlimits[1], 32), y = c(ylimits[2], 0.142, ylimits[2]), border = NA, col = "white")
  polygon(x = c(70, xlimits[2], xlimits[2]), y = c(ylimits[2], ylimits[2], 0.15), border = NA, col = "white")
  draw.contour.labels(xlimits = xlimits, ylimits, xfraction = 0.09, yfraction = 0.06, x.corners = c(38.25, 54.75, 34.6), y.corners = c(0.173, 0.181, 0.078), contour.labels = c("0.5", "1.0", "1.5"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(d.wv1$date, d.wv1$meanwave, cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c(paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ")), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.08", "0.12", "0.16", "0.20"), side = 2, line = 0.5, las = 1, at = c(0.08,0.12,0.16,0.2), cex = tick.label.size)
  mtext(text = "Mean wave height", side = 2, line = 2.1, cex = axis.label.size)
  mtext(text = "B. Wave wash", side = 3, line = 0.25, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size+0.1)
  #
  # C. Dragonfly numbers
par(new = "TRUE", plt = c(0.10, 0.43, 0.105, 0.455))
  ylimits    <- log10(range(newdat.dfly$aeshna))
  ylimits[2] <- log10(500)
  filled.contour3(x = unique(newdat.dfly$date), y = unique(log10(newdat.dfly$aeshna)), z = matrix(newdat.dfly$Estimate, nrow = length(unique(newdat.dfly$aeshna)), ncol = length(unique(newdat.dfly$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(newdat.dfly$date), y = unique(log10(newdat.dfly$aeshna)), z = matrix(newdat.dfly$Estimate, nrow = length(unique(newdat.dfly$aeshna)), ncol = length(unique(newdat.dfly$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  polygon(x = c(xlimits[1], xlimits[1], 26), y = c(ylimits[2], log10(350), ylimits[2]), border = NA, col = "white")
  polygon(x = c(71, xlimits[2], xlimits[2]), y = c(ylimits[2], ylimits[2], log10(310)), border = NA, col = "white")
  polygon(x = c(73, xlimits[2], xlimits[2]), y = c(ylimits[1], log10(200), ylimits[1]), border = NA, col = "white")
  polygon(x = c(xlimits[1], xlimits[1], 29), y = c(log10(220), ylimits[1], ylimits[1]), border = NA, col = "white")
  draw.contour.labels(xlimits, ylimits, xfraction = 0.09, yfraction = 0.06, x.corners = c(60.5, 37.5, 35.3), y.corners = log10(c(130, 130, 340)), contour.labels = c("0.5", "1.0", "1.5"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(d.dfly$date, log10(d.dfly$aeshna), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c(paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ")), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  mtext(text = "Hatching date", side = 1, line = 1.3, las = 1, cex = axis.label.size)
  axis(side = 2, tck = -0.025, at = log10(c(150,500)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log10(c(200,300,400)), labels = FALSE)
  mtext(text = c("150", "300", "500"), side = 2, line = 0.5, las = 1, at = log10(c(150,300,500)), cex = tick.label.size)
  mtext(text = "Number of dragonflies", side = 2, line = 2, cex = axis.label.size)
  mtext(text = "C. Predation", side = 3, line = 0.25, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size+0.1)
  #
  # D. Number of tadpoles, lower right
par(new = "TRUE", plt = c(0.54, 0.87, 0.105, 0.455))
  ylimits <- log10(range(newdat.dens$dens))
  filled.contour3(x = unique(newdat.dens$date), y = unique(log10(newdat.dens$dens)), z = matrix(newdat.dens$Estimate, nrow = length(unique(newdat.dens$dens)), ncol = length(unique(newdat.dens$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(newdat.dens$date), y = unique(log10(newdat.dens$dens)), z = matrix(newdat.dens$Estimate, nrow = length(unique(newdat.dens$dens)), ncol = length(unique(newdat.dens$date)), byrow = TRUE), ylab = "", drawlabels = T, add = TRUE, col = "black", nlevels = 6 )
  polygon(x = c(64, xlimits[2], xlimits[2]), y = c(ylimits[1], log10(2500), ylimits[1]), border = NA, col = "white")
  polygon(x = c(xlimits[1], xlimits[1], 39), y = c(log10(2200), ylimits[1], ylimits[1]), border = NA, col = "white")
  draw.contour.labels(xlimits, ylimits, xfraction = 0.09, yfraction = 0.06, x.corners = c(55, 62, 51), y.corners = log10(c(4200, 1820, 910)), contour.labels = c("0.5", "1.0", "1.5"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(d.dens$date, log10(d.dens$dens), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c(paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ")), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  mtext(text = "Hatching date", side = 1, line = 1.3, las = 1, cex = axis.label.size)
  axis(side = 2, tck = -0.022, at = log10(c(500,1000,5000)), labels = FALSE)
  axis(side = 2, tck = -0.012, at = log10(c(400,600,700,800,900,2000,3000,4000,6000,7000,8000)), labels = FALSE)
  mtext(text = c("500", "1000", "5000"), side = 2, line = 0.5, las = 1, at = log10(c(500,1000,5000)), cex = tick.label.size)
  mtext(text = "Number of tadpoles", side = 2, line = 2.2, cex = axis.label.size)
  mtext(text = "D. Competition", side = 3, line = 0.25, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size+0.1)
  #
  # Right-side scale.
  par(new = "TRUE", xpd = TRUE, plt = c(0.90, 0.925, 0.2975, 0.6975), las = 1, cex.axis = 1, tck = -0.25 )
  lower.limit <- color.range[1]
  upper.limit <- color.range[2]
  plot.breaks <- round( seq(lower.limit, upper.limit, length.out = 101), 2)
  image.scale(z = matrix(newdat.dens$Estimate, nrow = length(unique(newdat.dens$dens)), ncol = length(unique(newdat.dens$loc)), byrow = TRUE),
     zlim = c(lower.limit, upper.limit), col = plot.colors, xlab = "", ylab = "", cex.lab = 1.5, breaks = plot.breaks, horiz = FALSE, yaxt = "n", xaxt = "n")
    # Label the vertical axis.
  scale.labels <- c(0.5, 1, 1.5, 2)
  axis(4, at = scale.labels, labels = FALSE)
  mtext(text = c("0.5", "1.0", "1.5", "2.0"), side = 4, line = 0.45, las = 1, at = scale.labels, cex = 0.85)
  box(lwd = 1.2)
    # Label at the bottom.
  y.first.line <- lower.limit - 0.02
  text(labels = c("Relative", "fitness"), pos = 1, cex = 0.85, x = 0.6, y = c(y.first.line, y.first.line-0.1))
dev.off()
#
#
# Figure S7 -- How correlational selection on location/pool size is affected by ecological factors.
#
nr.extreme.years <- 3
fine.grid        <- 16
xlimits          <- c(0.2, 1.0)
ylimits          <- c(0.028, 20)
xlim.range       <- (xlimits - mean.loc.indivs.kept) / sd.loc.indivs.kept
ylim.range       <- (log(ylimits) - mean.larea.indivs.kept) / sd.larea.indivs.kept
d1               <- aggregate(d[,c("loc", "larea", "prop.dry", "mean.wave.ht", "total.aeshna", "Nr.tads")], by = list(year = d$year, pool = d$pool), mean, na.rm = TRUE)
d1$area          <- exp(d1$larea)
#
# Drying
dry                 <- dry[order(dry$prop.dry),]
low.dry.pools       <- d1[d1$year %in% dry$year[1:6], ]                 # Set aside the pools that were used in drier and wetter years.
low.dry.pools       <- aggregate(low.dry.pools[,c("loc", "area")], by = list(low.dry.pools$pool), mean)
hi.dry.pools        <- d1[d1$year %in% dry$year[11:16], ]
hi.dry.pools        <- aggregate(hi.dry.pools[,c("loc", "area")], by = list(hi.dry.pools$pool), mean)
low.dry.value       <- mean(dry$prop.dry[1:nr.extreme.years])           # low-drying years
low.dry.c           <- (low.dry.value - mean.propdry.kept) / sd.propdry.kept
hi.dry.value        <- mean(dry$prop.dry[(17-nr.extreme.years):16])     # high-drying years
hi.dry.c            <- (hi.dry.value - mean.propdry.kept) / sd.propdry.kept
newdat.low          <- expand.grid(pool = NA, loc.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), date.s = 0, prop.dry = low.dry.value, propdry.c = low.dry.c)
newdat.low$loc.s2   <- newdat.low$loc.s^2 / sd.loc.s2.indivs.kept
newdat.low$larea.s2 <- newdat.low$larea.s^2 / sd.larea.s2.indivs.kept
newdat.low$date.s2  <- 0
newdat.low          <- cbind(newdat.low, Estimate = predict(m1.drying, type = "response", re.form = NA, newdata = newdat.low)/10)
newdat.low$loc      <- (newdat.low$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.low$area     <- exp((newdat.low$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
newdat.hi           <- expand.grid(pool = NA, loc.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), date.s = 0, prop.dry = hi.dry.value, propdry.c = hi.dry.c)
newdat.hi$loc.s2    <- newdat.hi$loc.s^2 / sd.loc.s2.indivs.kept
newdat.hi$larea.s2  <- newdat.hi$larea.s^2 / sd.larea.s2.indivs.kept
newdat.hi$date.s2   <- 0
newdat.hi           <- cbind(newdat.hi, Estimate = predict(m1.drying, type = "response", re.form = NA, newdata = newdat.hi)/10)
newdat.hi$loc       <- (newdat.hi$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.hi$area      <- exp((newdat.hi$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
dat.dry.low         <- newdat.low
dat.dry.hi          <- newdat.hi
#
# Waves
waves               <- waves[order(waves$mean.wave.ht),]
low.wave.pools      <- d1[d1$year %in% waves$year[1:6], ]                          # Set aside the pools that were used in least and most wavy years.
low.wave.pools      <- aggregate(low.wave.pools[,c("loc", "area")], by = list(low.wave.pools$pool), mean)
hi.wave.pools       <- d1[d1$year %in% waves$year[11:16], ]
hi.wave.pools       <- aggregate(hi.wave.pools[,c("loc", "area")], by = list(hi.wave.pools$pool), mean)
low.wave.value      <- mean(waves$mean.wave.ht[1:nr.extreme.years])                # calmest years
low.wave.c          <- (low.wave.value - mean.meanwave.kept) / sd.meanwave.kept
hi.wave.value       <- mean(waves$mean.wave.ht[(17-nr.extreme.years):16])          # most wavy years
hi.wave.c           <- (hi.wave.value - mean.meanwave.kept) / sd.meanwave.kept
newdat.low          <- expand.grid(pool = NA, loc.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), date.s = 0, meanwave = low.wave.value, meanwave.c = low.wave.c)
newdat.low$loc.s2   <- newdat.low$loc.s^2 / sd.loc.s2.indivs.kept
newdat.low$larea.s2 <- newdat.low$larea.s^2 / sd.larea.s2.indivs.kept
newdat.low$date.s2  <- 0
newdat.low          <- cbind(newdat.low, Estimate = predict(m1.mean.wave, type = "response", re.form = NA, newdata = newdat.low)/10)
newdat.low$loc      <- (newdat.low$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.low$area     <- exp((newdat.low$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
newdat.hi           <- expand.grid(pool = NA, loc.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), date.s = 0, meanwave = hi.wave.value, meanwave.c = hi.wave.c)
newdat.hi$loc.s2    <- newdat.hi$loc.s^2 / sd.loc.s2.indivs.kept
newdat.hi$larea.s2  <- newdat.hi$larea.s^2 / sd.larea.s2.indivs.kept
newdat.hi$date.s2   <- 0
newdat.hi           <- cbind(newdat.hi, Estimate = predict(m1.mean.wave, type = "response", re.form = NA, newdata = newdat.hi)/10)
newdat.hi$loc       <- (newdat.hi$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.hi$area      <- exp((newdat.hi$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
dat.waves.low       <- newdat.low
dat.waves.hi        <- newdat.hi
#
# Aeshna
dflies              <- subset(total.aeshna, year < 1999 & year > 1985)
dflies              <- dflies[order(dflies$total.aeshna),]
low.dfly.pools      <- d1[d1$year %in% dflies$year[1:6], ]                      # Set aside the pools that were used in least and most Aeshna years.
low.dfly.pools      <- aggregate(low.dfly.pools[,c("loc", "area")], by = list(low.dfly.pools$pool), mean)
hi.dfly.pools       <- d1[d1$year %in% dflies$year[8:13], ]
hi.dfly.pools       <- aggregate(hi.dfly.pools[,c("loc", "area")], by = list(hi.dfly.pools$pool), mean)
low.dflies.value    <- mean(dflies$total.aeshna[1:nr.extreme.years])            # years with fewest dragonflies
low.dflies.c        <- (log(low.dflies.value) - mean.aeshna.kept) / sd.aeshna.kept
hi.dflies.value     <- mean(dflies$total.aeshna[(14-nr.extreme.years):13])      # years with most dragonflies
hi.dflies.c         <- (log(hi.dflies.value) - mean.aeshna.kept) / sd.aeshna.kept
newdat.low          <- expand.grid(pool = NA, loc.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), date.s = 0, aeshna = low.dflies.value, aeshna.c = low.dflies.c)
newdat.low$loc.s2   <- newdat.low$loc.s^2 / sd.loc.s2.indivs.kept
newdat.low$larea.s2 <- newdat.low$larea.s^2 / sd.larea.s2.indivs.kept
newdat.low$date.s2  <- 0
newdat.low          <- cbind(newdat.low, Estimate = predict(m1.Aeshna, type = "response", re.form = NA, newdata = newdat.low)/10)
newdat.low$loc      <- (newdat.low$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.low$area     <- exp((newdat.low$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
newdat.hi           <- expand.grid(pool = NA, loc.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), date.s = 0, aeshna = hi.dflies.value, aeshna.c = hi.dflies.c)
newdat.hi$loc.s2    <- newdat.hi$loc.s^2 / sd.loc.s2.indivs.kept
newdat.hi$larea.s2  <- newdat.hi$larea.s^2 / sd.larea.s2.indivs.kept
newdat.hi$date.s2   <- 0
newdat.hi           <- cbind(newdat.hi, Estimate = predict(m1.Aeshna, type = "response", re.form = NA, newdata = newdat.hi)/10)
newdat.hi$loc       <- (newdat.hi$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.hi$area      <- exp((newdat.hi$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
dat.dfly.low        <- newdat.low
dat.dfly.hi         <- newdat.hi
#
# Tadpole numbers
dens4               <- dens4[order(dens4$Nr.tads),]
low.dens.pools      <- d1[d1$year %in% dens4$year[1:6], ]                  # Set aside the pools that were used in lowest and highest density years.
low.dens.pools      <- aggregate(low.dens.pools[,c("loc", "area")], by = list(low.dens.pools$pool), mean)
hi.dens.pools       <- d1[d1$year %in% dens4$year[11:16], ]
hi.dens.pools       <- aggregate(hi.dens.pools[,c("loc", "area")], by = list(hi.dens.pools$pool), mean)
low.density.value   <- mean(dens4$Nr.tads[1:nr.extreme.years])             # least crowded years
low.ldensity.c      <- (log(low.density.value) - mean.Nrtads.kept) / sd.Nrtads.kept
hi.density.value    <- mean(dens4$Nr.tads[(17-nr.extreme.years):16])       # most crowded years
hi.ldensity.c       <- (log(hi.density.value) - mean.Nrtads.kept) / sd.Nrtads.kept
newdat.low          <- expand.grid(pool = NA, loc.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), date.s = 0, Nr.tads = low.density.value, Nr.tads.c = low.ldensity.c)
newdat.low$loc.s2   <- newdat.low$loc.s^2 / sd.loc.s2.indivs.kept
newdat.low$larea.s2 <- newdat.low$larea.s^2 / sd.larea.s2.indivs.kept
newdat.low$date.s2  <- 0
newdat.low          <- cbind(newdat.low, Estimate = predict(m1.Nrtads, type = "response", re.form = NA, newdata = newdat.low)/10)
newdat.low$loc      <- (newdat.low$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.low$area     <- exp((newdat.low$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
newdat.hi           <- expand.grid(pool = NA, loc.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), date.s = 0, Nr.tads = hi.density.value, Nr.tads.c = hi.ldensity.c)
newdat.hi$loc.s2    <- newdat.hi$loc.s^2 / sd.loc.s2.indivs.kept
newdat.hi$larea.s2  <- newdat.hi$larea.s^2 / sd.larea.s2.indivs.kept
newdat.hi$date.s2   <- 0
newdat.hi           <- cbind(newdat.hi, Estimate = predict(m1.Nrtads, type = "response", re.form = NA, newdata = newdat.hi)/10)
newdat.hi$loc       <- (newdat.hi$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.hi$area      <- exp((newdat.hi$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
dat.dens.low        <- newdat.low
dat.dens.hi         <- newdat.hi
#
lower.cutoff <- 0.1
upper.cutoff <- 1.6
dat.dry.hi$Estimate[dat.dry.hi$Estimate < lower.cutoff] <- lower.cutoff
dat.waves.hi$Estimate[dat.waves.hi$Estimate < lower.cutoff] <- lower.cutoff
dat.dfly.hi$Estimate[dat.dfly.hi$Estimate < lower.cutoff] <- lower.cutoff
dat.dens.hi$Estimate[dat.dens.hi$Estimate < lower.cutoff] <- lower.cutoff
dat.dry.low$Estimate[dat.dry.low$Estimate < lower.cutoff] <- lower.cutoff
dat.waves.low$Estimate[dat.waves.low$Estimate < lower.cutoff] <- lower.cutoff
dat.dfly.low$Estimate[dat.dfly.low$Estimate < lower.cutoff] <- lower.cutoff
dat.dens.low$Estimate[dat.dens.low$Estimate < lower.cutoff] <- lower.cutoff
dat.dry.low$Estimate[dat.dry.low$Estimate > upper.cutoff] <- upper.cutoff
dat.waves.low$Estimate[dat.waves.low$Estimate > upper.cutoff] <- upper.cutoff
dat.dfly.low$Estimate[dat.dfly.low$Estimate > upper.cutoff] <- upper.cutoff
dat.dens.low$Estimate[dat.dens.low$Estimate > upper.cutoff] <- upper.cutoff
color.range        <- range(dat.dry.low$Estimate, dat.waves.low$Estimate, dat.dfly.low$Estimate, dat.dens.low$Estimate, dat.dry.hi$Estimate, dat.waves.hi$Estimate, dat.dfly.hi$Estimate, dat.dens.hi$Estimate)
point.size         <- 0.6
contour.label.size <- 0.85
tick.label.size    <- 0.9
axis.label.size    <- 1.2
draw.the.polygons  <- function(xlimits, ylimits) {
     polygon(x = c(xlimits[1], xlimits[1], 0.6), y = c(log(7), log(ylimits[2]), log(ylimits[2])), border = NA, col = "white")
     polygon(x = c(0.88, xlimits[2], xlimits[2]), y = c(log(ylimits[2]), log(ylimits[2]), log(4)), border = NA, col = "white")
     polygon(x = c(0.78, xlimits[2], xlimits[2]), y = c(log(ylimits[1]), log(0.17), log(ylimits[1])), border = NA, col = "white")
     polygon(x = c(xlimits[1], xlimits[1], 0.3), y = c(log(0.15), log(ylimits[1]), log(ylimits[1])), border = NA, col = "white")
     }
color.gradient     <- function(x, colors = c(lower.color, lowmid.color, highmid.color, upper.color), colsteps = 100) { return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )  }
plot.colors        <- c(rep(color.gradient(1:2)[1], 10), color.gradient(1:90))
plot.colors        <- color.gradient(1:100)
#
setwd(figs.wd)
pdf("Fig S7.pdf", useDingbats = FALSE, width = 8, height = 14)
plot.new()
  # Drying left side, panel A
par(new = "TRUE", plt = c(0.13, 0.48, 0.77, 0.94))
  filled.contour3(x = unique(dat.dry.low$loc), y = log(unique(dat.dry.low$area)), z = matrix(dat.dry.low$Estimate, nrow = length(unique(dat.dry.low$area)), ncol = length(unique(dat.dry.low$loc)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dry.low$loc), y = log(unique(dat.dry.low$area)), z = matrix(dat.dry.low$Estimate, nrow = length(unique(dat.dry.low$area)), ncol = length(unique(dat.dry.low$loc)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  draw.the.polygons(xlimits, ylimits)
  points(low.dry.pools$loc, log(low.dry.pools$area), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.093, yfraction = 0.06465, x.corners = c(0.755, 0.675, 0.717), y.corners = c(log(0.057), log(1), log(7.5)), contour.labels = c("0.4", "0.8", "1.0"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 1, line = 0.3, las = 1, at = c(0.2,0.4,0.6,0.8,1), cex = tick.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = "Pool size (m  )", side = 2, line = 1.7, cex = axis.label.size)
  mtext(text = "2", side = 2, at = log(3.0), line = 2.1, cex = axis.label.size-0.3)
  mtext(text = paste("A. Low pool drying (", format.p.val(low.dry.value, 2), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  mtext(text = paste("Figure S7, glmmTMB,", Sys.Date()), side = 3, line = 3, las = 1, cex = 0.85)
  # Drying right side, panel B
par(new = "TRUE", plt = c(0.57, 0.92, 0.77, 0.94))
  filled.contour3(x = unique(dat.dry.hi$loc), y = log(unique(dat.dry.hi$area)), z = matrix(dat.dry.hi$Estimate, nrow = length(unique(dat.dry.hi$area)), ncol = length(unique(dat.dry.hi$loc)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dry.hi$loc), y = log(unique(dat.dry.hi$area)), z = matrix(dat.dry.hi$Estimate, nrow = length(unique(dat.dry.hi$area)), ncol = length(unique(dat.dry.hi$loc)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  draw.the.polygons(xlimits, ylimits)
  box(lwd = 1.2)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.093, yfraction = 0.06465, x.corners = c(0.47, 0.78, 0.81), y.corners = c(log(3.7), log(1.5), log(0.28)), contour.labels = c("0.4", "0.8", "1.2"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.dry.pools$loc, log(hi.dry.pools$area), cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 1, line = 0.3, las = 1, at = c(0.2,0.4,0.6,0.8,1), cex = tick.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = paste("B. High pool drying (", format.p.val(hi.dry.value, 2), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  #
  # Waves left side, panel C
par(new = "TRUE", plt = c(0.13, 0.48, 0.55, 0.72))
  filled.contour3(x = unique(dat.waves.low$loc), y = log(unique(dat.waves.low$area)), z = matrix(dat.waves.low$Estimate, nrow = length(unique(dat.waves.low$area)), ncol = length(unique(dat.waves.low$loc)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.waves.low$loc), y = log(unique(dat.waves.low$area)), z = matrix(dat.waves.low$Estimate, nrow = length(unique(dat.waves.low$area)), ncol = length(unique(dat.waves.low$loc)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  draw.the.polygons(xlimits, ylimits)
  box(lwd = 1.2)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.093, yfraction = 0.06465, x.corners = c(0.465, 0.76, 0.665), y.corners = c(log(5.8), log(1.1), log(0.07)), contour.labels = c("0.4", "1.0", "1.2"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(low.wave.pools$loc, log(low.wave.pools$area), cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 1, line = 0.3, las = 1, at = c(0.2,0.4,0.6,0.8,1), cex = tick.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = "Pool size (m  )", side = 2, line = 1.7, cex = axis.label.size)
  mtext(text = "2", side = 2, at = log(3), line = 2.1, cex = axis.label.size-0.3)
  mtext(text = paste("C. Low wave height (", format.p.val(low.wave.value, 2), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  # Waves right side, panel D
par(new = "TRUE", plt = c(0.57, 0.92, 0.55, 0.72))
  filled.contour3(x = unique(dat.waves.hi$loc), y = log(unique(dat.waves.hi$area)), z = matrix(dat.waves.hi$Estimate, nrow = length(unique(dat.waves.hi$area)), ncol = length(unique(dat.waves.hi$loc)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.waves.hi$loc), y = log(unique(dat.waves.hi$area)), z = matrix(dat.waves.hi$Estimate, nrow = length(unique(dat.waves.hi$area)), ncol = length(unique(dat.waves.hi$loc)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  draw.the.polygons(xlimits, ylimits)
  box(lwd = 1.2)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.093, yfraction = 0.06465, x.corners = c(0.39, 0.57, 0.65), y.corners = c(log(0.23), log(1.0), log(9.5)), contour.labels = c("0.4", "0.7", "0.8"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.wave.pools$loc, log(hi.wave.pools$area), cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 1, line = 0.3, las = 1, at = c(0.2,0.4,0.6,0.8,1), cex = tick.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = paste("D. High wave height (", format.p.val(hi.wave.value, 2), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  #
  # Dragonflies left side, panel E
par(new = "TRUE", plt = c(0.13, 0.48, 0.33, 0.50))
  filled.contour3(x = unique(dat.dfly.low$loc), y = log(unique(dat.dfly.low$area)), z = matrix(dat.dfly.low$Estimate, nrow = length(unique(dat.dfly.low$area)), ncol = length(unique(dat.dfly.low$loc)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dfly.low$loc), y = log(unique(dat.dfly.low$area)), z = matrix(dat.dfly.low$Estimate, nrow = length(unique(dat.dfly.low$area)), ncol = length(unique(dat.dfly.low$loc)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", levels = c(0.5,1,1.5) )
  draw.the.polygons(xlimits, ylimits)
  points(low.dfly.pools$loc, log(low.dfly.pools$area), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.093, yfraction = 0.06465, x.corners = c(0.69, 0.61, 0.842), y.corners = c(log(6.7), log(1.0), log(0.31)), contour.labels = c("0.5", "1.0", "1.5"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.wave.pools$loc, log(hi.wave.pools$area), cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 1, line = 0.3, las = 1, at = c(0.2,0.4,0.6,0.8,1), cex = tick.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = "Pool size (m  )", side = 2, line = 1.7, cex = axis.label.size)
  mtext(text = "2", side = 2, at = log(3), line = 2.1, cex = axis.label.size-0.3)
  mtext(paste("E. Low dragonfly abundance (", format.p.val(low.dflies.value, 0), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  # Dragonflies right side, panel F
par(new = "TRUE", plt = c(0.57, 0.92, 0.33, 0.50))
  filled.contour3(x = unique(dat.dfly.hi$loc), y = log(unique(dat.dfly.hi$area)), z = matrix(dat.dfly.hi$Estimate, nrow = length(unique(dat.dfly.hi$area)), ncol = length(unique(dat.dfly.hi$loc)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dfly.hi$loc), y = log(unique(dat.dfly.hi$area)), z = matrix(dat.dfly.hi$Estimate, nrow = length(unique(dat.dfly.hi$area)), ncol = length(unique(dat.dfly.hi$loc)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 4 )
  draw.the.polygons(xlimits, ylimits)
  points(hi.dfly.pools$loc, log(hi.dfly.pools$area), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.093, yfraction = 0.06465, x.corners = c(0.60, 0.847), y.corners = c(log(1.2), log(0.31)), contour.labels = c("0.5", "1.0"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.wave.pools$loc, log(hi.wave.pools$area), cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 1, line = 0.3, las = 1, at = c(0.2,0.4,0.6,0.8,1), cex = tick.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = paste("F. High dragonfly abundance (", format.p.val(hi.dflies.value, 0), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  #
  # Tadpole density left side, panel G
par(new = "TRUE", plt = c(0.13, 0.48, 0.11, 0.28))
  filled.contour3(x = unique(dat.dens.low$loc), y = log(unique(dat.dens.low$area)), z = matrix(dat.dens.low$Estimate, nrow = length(unique(dat.dens.low$area)), ncol = length(unique(dat.dens.low$loc)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dens.low$loc), y = log(unique(dat.dens.low$area)), z = matrix(dat.dens.low$Estimate, nrow = length(unique(dat.dens.low$area)), ncol = length(unique(dat.dens.low$loc)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  draw.the.polygons(xlimits, ylimits)
  points(low.dens.pools$loc, log(low.dens.pools$area), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.093, yfraction = 0.06465, x.corners = c(0.285, 0.455, 0.58), y.corners = c(log(2.3), log(4), log(0.8)), contour.labels = c("0.4", "1.2", "1.4"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.wave.pools$loc, log(hi.wave.pools$area), cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 1, line = 0.3, las = 1, at = c(0.2,0.4,0.6,0.8,1), cex = tick.label.size)
  mtext(text = "Pool location", side = 1, line = 1.4, cex = axis.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = "Pool size (m  )", side = 2, line = 1.7, cex = axis.label.size)
  mtext(text = "2", side = 2, at = log(3), line = 2.1, cex = axis.label.size-0.3)
  mtext(text = paste("G. Low tadpole numbers (", format.p.val(low.density.value, 0), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  par(xpd = TRUE)
  mtext(text = c("Lake", "Forest"), side = 1, line = 1.3, las = 1, at = c(xlimits[1]+0.09, xlimits[2]-0.09), adj = c(0,1), cex = tick.label.size, las = 0)
  arrows(x0 = xlimits[1]+0.08, y0 = log(0.01), x1 = xlimits[1], y1 = log(0.01), angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  arrows(x0 = xlimits[2]-0.08, y0 = log(0.01), x1 = xlimits[2], y1 = log(0.01), angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  par(xpd = FALSE)
  # Tadpole density right side, panel H
par(new = "TRUE", plt = c(0.57, 0.92, 0.11, 0.28))
  filled.contour3(x = unique(dat.dens.hi$loc), y = log(unique(dat.dens.hi$area)), z = matrix(dat.dens.hi$Estimate, nrow = length(unique(dat.dens.hi$area)), ncol = length(unique(dat.dens.hi$loc)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dens.hi$loc), y = log(unique(dat.dens.hi$area)), z = matrix(dat.dens.hi$Estimate, nrow = length(unique(dat.dens.hi$area)), ncol = length(unique(dat.dens.hi$loc)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  draw.the.polygons(xlimits, ylimits)
  points(hi.dens.pools$loc, log(hi.dens.pools$area), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.093, yfraction = 0.06465, x.corners = c(0.48, 0.75, 0.38), y.corners = c(log(3.45), log(1.17), log(0.265)), contour.labels = c("0.4", "0.8", "1.0"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.wave.pools$loc, log(hi.wave.pools$area), cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 1, line = 0.3, las = 1, at = c(0.2,0.4,0.6,0.8,1), cex = tick.label.size)
  mtext(text = "Pool location", side = 1, line = 1.4, cex = axis.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = paste("H. High tadpole numbers (", hi.density.value, ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  par(xpd = TRUE)
  mtext(text = c("Lake", "Forest"), side = 1, line = 1.3, las = 1, at = c(xlimits[1]+0.09, xlimits[2]-0.09), adj = c(0,1), cex = tick.label.size, las = 0)
  arrows(x0 = xlimits[1]+0.08, y0 = log(0.01), x1 = xlimits[1], y1 = log(0.01), angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  arrows(x0 = xlimits[2]-0.08, y0 = log(0.01), x1 = xlimits[2], y1 = log(0.01), angle = 30, code = 2, lwd = 1.1, length = 0.1 )
  par(xpd = FALSE)
dev.off()
#
#
# Figure 7 -- How correlational selection on location/date is affected by ecological factors.
#
nr.extreme.years <- 3
fine.grid        <- 19
xlimits          <- c(23, 73.5)
ylimits          <- c(0.18, 1.0)
xlim.range       <- (xlimits - mean.date.indivs.kept) / sd.date.indivs.kept
ylim.range       <- (ylimits - mean.loc.indivs.kept) / sd.loc.indivs.kept
d1               <- aggregate(d[,c("loc", "date", "prop.dry", "mean.wave.ht", "total.aeshna", "Nr.tads")], by = list(year = d$year, pool = d$pool), mean, na.rm = TRUE)
#
# Drying
dry                 <- dry[order(dry$prop.dry),]
low.dry.pools       <- d1[d1$year %in% dry$year[1:6], ]                 # Set aside the pools that were used in drier and wetter years.
low.dry.pools       <- aggregate(low.dry.pools[,c("pool")], by = list(loc = low.dry.pools$loc, date = low.dry.pools$date), mean)
hi.dry.pools        <- d1[d1$year %in% dry$year[11:16], ]
hi.dry.pools        <- aggregate(hi.dry.pools[,c("pool")], by = list(loc = hi.dry.pools$loc, date = hi.dry.pools$date), mean)
low.dry.value       <- mean(dry$prop.dry[1:nr.extreme.years])           # low-drying years
low.dry.c           <- (low.dry.value - mean.propdry.kept) / sd.propdry.kept
hi.dry.value        <- mean(dry$prop.dry[(17-nr.extreme.years):16])     # high-drying years
hi.dry.c            <- (hi.dry.value - mean.propdry.kept) / sd.propdry.kept
newdat.low          <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), loc.s = set.up.interval(ylim.range, fine.grid), larea.s = 0, prop.dry = low.dry.value, propdry.c = low.dry.c)
newdat.low$loc.s2   <- newdat.low$loc.s^2 / sd.loc.s2.indivs.kept
newdat.low$larea.s2 <- 0
newdat.low$date.s2  <- newdat.low$date.s^2 / sd.date.s2.indivs.kept
newdat.low          <- cbind(newdat.low, Estimate = predict(m1.drying, type = "response", re.form = NA, newdata = newdat.low)/10)
newdat.low$loc      <- (newdat.low$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.low$date     <- (newdat.low$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.hi           <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), loc.s = set.up.interval(ylim.range, fine.grid), larea.s = 0, prop.dry = hi.dry.value, propdry.c = hi.dry.c)
newdat.hi$loc.s2    <- newdat.hi$loc.s^2 / sd.loc.s2.indivs.kept
newdat.hi$larea.s2  <- 0
newdat.hi$date.s2   <- newdat.hi$date.s^2 / sd.date.s2.indivs.kept
newdat.hi           <- cbind(newdat.hi, Estimate = predict(m1.drying, type = "response", re.form = NA, newdata = newdat.hi)/10)
newdat.hi$loc       <- (newdat.hi$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.hi$date      <- (newdat.hi$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
dat.dry.low         <- newdat.low
dat.dry.hi          <- newdat.hi
#
# Waves
waves               <- waves[order(waves$mean.wave.ht),]
low.wave.pools      <- d1[d1$year %in% waves$year[1:6], ]                          # Set aside the pools that were used in least and most wavy years.
low.wave.pools      <- aggregate(low.wave.pools[,c("pool")], by = list(loc = low.wave.pools$loc, date = low.wave.pools$date), mean)
hi.wave.pools       <- d1[d1$year %in% waves$year[11:16], ]
hi.wave.pools       <- aggregate(hi.wave.pools[,c("pool")], by = list(loc = hi.wave.pools$loc, date = hi.wave.pools$date), mean)
low.wave.value      <- mean(waves$mean.wave.ht[1:nr.extreme.years])                # calmest years
low.wave.c          <- (low.wave.value - mean.meanwave.kept) / sd.meanwave.kept
hi.wave.value       <- mean(waves$mean.wave.ht[(17-nr.extreme.years):16])          # most wavy years
hi.wave.c           <- (hi.wave.value - mean.meanwave.kept) / sd.meanwave.kept
newdat.low          <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), loc.s = set.up.interval(ylim.range, fine.grid), larea.s = 0, meanwave = low.wave.value, meanwave.c = low.wave.c)
newdat.low$loc.s2   <- newdat.low$loc.s^2 / sd.loc.s2.indivs.kept
newdat.low$larea.s2 <- 0
newdat.low$date.s2  <- newdat.low$date.s^2 / sd.date.s2.indivs.kept
newdat.low          <- cbind(newdat.low, Estimate = predict(m1.mean.wave, type = "response", re.form = NA, newdata = newdat.low)/10)
newdat.low$loc      <- (newdat.low$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.low$date     <- (newdat.low$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.hi           <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), loc.s = set.up.interval(ylim.range, fine.grid), larea.s = 0, meanwave = hi.wave.value, meanwave.c = hi.wave.c)
newdat.hi$loc.s2    <- newdat.hi$loc.s^2 / sd.loc.s2.indivs.kept
newdat.hi$larea.s2  <- 0
newdat.hi$date.s2   <- newdat.hi$date.s^2 / sd.date.s2.indivs.kept
newdat.hi           <- cbind(newdat.hi, Estimate = predict(m1.mean.wave, type = "response", re.form = NA, newdata = newdat.hi)/10)
newdat.hi$loc       <- (newdat.hi$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.hi$date      <- (newdat.hi$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
dat.waves.low       <- newdat.low
dat.waves.hi        <- newdat.hi
#
# Aeshna
dflies              <- subset(total.aeshna, year < 1999 & year > 1985)
dflies              <- dflies[order(dflies$total.aeshna),]
low.dfly.pools      <- d1[d1$year %in% dflies$year[1:6], ]                      # Set aside the pools that were used in least and most Aeshna years.
low.dfly.pools      <- aggregate(low.dfly.pools[,c("pool")], by = list(loc = low.dfly.pools$loc, date = low.dfly.pools$date), mean)
hi.dfly.pools       <- d1[d1$year %in% dflies$year[8:13], ]
hi.dfly.pools       <- aggregate(hi.dfly.pools[,c("pool")], by = list(loc = hi.dfly.pools$loc, date = hi.dfly.pools$date), mean)
low.dflies.value    <- mean(dflies$total.aeshna[1:nr.extreme.years])            # years with fewest dragonflies
low.dflies.c        <- (log(low.dflies.value) - mean.aeshna.kept) / sd.aeshna.kept
hi.dflies.value     <- mean(dflies$total.aeshna[(14-nr.extreme.years):13])      # years with most dragonflies
hi.dflies.c         <- (log(hi.dflies.value) - mean.aeshna.kept) / sd.aeshna.kept
newdat.low          <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), loc.s = set.up.interval(ylim.range, fine.grid), larea.s = 0, aeshna = low.dflies.value, aeshna.c = low.dflies.c)
newdat.low$loc.s2   <- newdat.low$loc.s^2 / sd.loc.s2.indivs.kept
newdat.low$larea.s2 <- 0
newdat.low$date.s2  <- newdat.low$date.s^2 / sd.date.s2.indivs.kept
newdat.low          <- cbind(newdat.low, Estimate = predict(m1.Aeshna, type = "response", re.form = NA, newdata = newdat.low)/10)
newdat.low$loc      <- (newdat.low$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.low$date     <- (newdat.low$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.hi           <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), loc.s = set.up.interval(ylim.range, fine.grid), larea.s = 0, aeshna = hi.dflies.value, aeshna.c = hi.dflies.c)
newdat.hi$loc.s2    <- newdat.hi$loc.s^2 / sd.loc.s2.indivs.kept
newdat.hi$larea.s2  <- 0
newdat.hi$date.s2   <- newdat.hi$date.s^2 / sd.date.s2.indivs.kept
newdat.hi           <- cbind(newdat.hi, Estimate = predict(m1.Aeshna, type = "response", re.form = NA, newdata = newdat.hi)/10)
newdat.hi$loc       <- (newdat.hi$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.hi$date      <- (newdat.hi$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
dat.dfly.low        <- newdat.low
dat.dfly.hi         <- newdat.hi
#
# Tadpole density
dens4               <- dens4[order(dens4$Nr.tads),]
low.dens.pools      <- d1[d1$year %in% dens4$year[1:6], ]                  # Set aside the pools that were used in lowest and highest density years.
low.dens.pools      <- aggregate(low.dens.pools[,c("pool")], by = list(loc = low.dens.pools$loc, date = low.dens.pools$date), mean)
hi.dens.pools       <- d1[d1$year %in% dens4$year[11:16], ]
hi.dens.pools       <- aggregate(hi.dens.pools[,c("pool")], by = list(loc = hi.dens.pools$loc, date = hi.dens.pools$date), mean)
low.density.value   <- mean(dens4$Nr.tads[1:nr.extreme.years])            # least crowded years
low.ldensity.c      <- (log(low.density.value) - mean.Nrtads.kept) / sd.Nrtads.kept
hi.density.value    <- mean(dens4$Nr.tads[(17-nr.extreme.years):16])      # most crowded years
hi.ldensity.c       <- (log(hi.density.value) - mean.Nrtads.kept) / sd.Nrtads.kept
newdat.low          <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), loc.s = set.up.interval(ylim.range, fine.grid), larea.s = 0, Nr.tads = low.density.value, Nr.tads.c = low.ldensity.c)
newdat.low$loc.s2   <- newdat.low$loc.s^2 / sd.loc.s2.indivs.kept
newdat.low$larea.s2 <- 0
newdat.low$date.s2  <- newdat.low$date.s^2 / sd.date.s2.indivs.kept
newdat.low          <- cbind(newdat.low, Estimate = predict(m1.Nrtads, type = "response", re.form = NA, newdata = newdat.low)/10)

newdat.low$loc      <- (newdat.low$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.low$date     <- (newdat.low$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.hi           <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), loc.s = set.up.interval(ylim.range, fine.grid), larea.s = 0, Nr.tads = hi.density.value, Nr.tads.c = hi.ldensity.c)
newdat.hi$loc.s2    <- newdat.hi$loc.s^2 / sd.loc.s2.indivs.kept
newdat.hi$larea.s2  <- 0
newdat.hi$date.s2  <- newdat.hi$date.s^2 / sd.date.s2.indivs.kept
newdat.hi          <- cbind(newdat.hi, Estimate = predict(m1.Nrtads, type = "response", re.form = NA, newdata = newdat.hi)/10)
newdat.hi$loc      <- (newdat.hi$loc.s * sd.loc.indivs.kept) + mean.loc.indivs.kept
newdat.hi$date     <- (newdat.hi$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
dat.dens.low       <- newdat.low
dat.dens.hi        <- newdat.hi
#
lower.cutoff <- 0.0001
upper.cutoff <- 1.75
dat.dry.hi$Estimate[dat.dry.hi$Estimate < lower.cutoff]       <- lower.cutoff
dat.waves.hi$Estimate[dat.waves.hi$Estimate < lower.cutoff]   <- lower.cutoff
dat.dfly.hi$Estimate[dat.dfly.hi$Estimate < lower.cutoff]     <- lower.cutoff
dat.dens.hi$Estimate[dat.dens.hi$Estimate < lower.cutoff]     <- lower.cutoff
dat.dry.low$Estimate[dat.dry.low$Estimate < lower.cutoff]     <- lower.cutoff
dat.waves.low$Estimate[dat.waves.low$Estimate < lower.cutoff] <- lower.cutoff
dat.dfly.low$Estimate[dat.dfly.low$Estimate < lower.cutoff]   <- lower.cutoff
dat.dens.low$Estimate[dat.dens.low$Estimate < lower.cutoff]   <- lower.cutoff
dat.dry.low$Estimate[dat.dry.low$Estimate > upper.cutoff]     <- upper.cutoff
dat.waves.low$Estimate[dat.waves.low$Estimate > upper.cutoff] <- upper.cutoff
dat.dfly.low$Estimate[dat.dfly.low$Estimate > upper.cutoff]   <- upper.cutoff
dat.dens.low$Estimate[dat.dens.low$Estimate > upper.cutoff]   <- upper.cutoff
dat.dry.hi$Estimate[dat.dry.hi$Estimate > upper.cutoff]      <- upper.cutoff
dat.waves.hi$Estimate[dat.waves.hi$Estimate > upper.cutoff]  <- upper.cutoff
dat.dfly.hi$Estimate[dat.dfly.hi$Estimate > upper.cutoff]    <- upper.cutoff
dat.dens.hi$Estimate[dat.dens.hi$Estimate > upper.cutoff]    <- upper.cutoff
color.range        <- range(dat.dry.low$Estimate, dat.waves.low$Estimate, dat.dfly.low$Estimate, dat.dens.low$Estimate, dat.dry.hi$Estimate, dat.waves.hi$Estimate, dat.dfly.hi$Estimate, dat.dens.hi$Estimate)
point.size         <- 0.6
contour.label.size <- 0.85
tick.label.size    <- 0.9
axis.label.size    <- 1.2
draw.the.polygons  <- function(xlimits, ylimits) {
     polygon(x = c(xlimits[1], xlimits[1], 32), y = c(0.92, ylimits[2], ylimits[2]), border = NA, col = "white")
     polygon(x = c(63, xlimits[2], xlimits[2]), y = c(ylimits[2], ylimits[2], 0.81), border = NA, col = "white")
     polygon(x = c(67.5, xlimits[2], xlimits[2]), y = c(ylimits[1], 0.38, ylimits[1]), border = NA, col = "white")
     polygon(x = c(xlimits[1], xlimits[1], 29), y = c(0.38, ylimits[1], ylimits[1]), border = NA, col = "white")
     }
color.gradient <- function(x, colors = c(lower.color, lowmid.color, highmid.color, upper.color), colsteps = 100) { return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )  }
plot.colors    <- color.gradient(1:100)
#
setwd(figs.wd)
pdf("Fig 7.pdf", useDingbats = FALSE, width = 8, height = 14)
plot.new()
  # Drying left side, panel A
par(new = "TRUE", plt = c(0.13, 0.48, 0.77, 0.94))
  filled.contour3(x = unique(dat.dry.low$date), y = unique(dat.dry.low$loc), z = matrix(dat.dry.low$Estimate, nrow = length(unique(dat.dry.low$loc)), ncol = length(unique(dat.dry.low$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dry.low$date), y = unique(dat.dry.low$loc), z = matrix(dat.dry.low$Estimate, nrow = length(unique(dat.dry.low$loc)), ncol = length(unique(dat.dry.low$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  draw.the.polygons(xlimits, ylimits)
  draw.contour.labels(xlimits, ylimits, xfraction = 0.09, yfraction = 0.06465, x.corners = c(25, 27.8, 57), y.corners = c(0.495, 0.6, 0.78), contour.labels = c("0.6", "0.8", "1.2"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(low.dry.pools$date, low.dry.pools$loc, cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 2, line = 0.5, at = c(0.2,0.4,0.6,0.8,1), las = 1, cex = tick.label.size)
  mtext(text = "Pool location", side = 2, line = 1.7, at = 0.585, cex = axis.label.size)
  par(xpd = TRUE)
  mtext(text = c("Lake", "Forest"), side = 2, line = 1.65, las = 0, at = c(ylimits[1]+0.075, ylimits[2]-0.065), adj = c(0,1), cex = tick.label.size-0.05)
  arrows(y0 = ylimits[1]+0.07, x0 = 15.5, y1 = ylimits[1]-0.01, x1 = 15.5, angle = 30, code = 2, lwd = 1.0, length = 0.1 )
  arrows(y0 = ylimits[2]-0.06, x0 = 15.5, y1 = ylimits[2]+0.02, x1 = 15.5, angle = 30, code = 2, lwd = 1.0, length = 0.1 )
  par(xpd = FALSE)
  mtext(text = paste("A. Low pool drying (", format.p.val(low.dry.value, 2), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  mtext(text = paste("Van Buskirk & Smith, Figure 7,", Sys.Date()), side = 3, line = 3, las = 1, cex = 0.8)
  # Drying right side, panel B
par(new = "TRUE", plt = c(0.57, 0.92, 0.77, 0.94))
  filled.contour3(x = unique(dat.dry.hi$date), y = unique(dat.dry.hi$loc), z = matrix(dat.dry.hi$Estimate, nrow = length(unique(dat.dry.hi$loc)), ncol = length(unique(dat.dry.hi$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dry.hi$date), y = unique(dat.dry.hi$loc), z = matrix(dat.dry.hi$Estimate, nrow = length(unique(dat.dry.hi$loc)), ncol = length(unique(dat.dry.hi$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  draw.the.polygons(xlimits, ylimits)
  box(lwd = 1.2)
  draw.contour.labels(xlimits, ylimits, xfraction = 0.09, yfraction = 0.06465, x.corners = c(64, 27.8, 23.5), y.corners = c(0.345, 0.545, 0.47), contour.labels = c("0.5", "1.0", "1.5"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.dry.pools$date, hi.dry.pools$loc, cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 2, line = 0.5, at = c(0.2,0.4,0.6,0.8,1), las = 1, cex = tick.label.size)
  mtext(text = paste("B. High pool drying (", format.p.val(hi.dry.value, 2), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  #
  # Waves left side, panel C
par(new = "TRUE", plt = c(0.13, 0.48, 0.55, 0.72))
  filled.contour3(x = unique(dat.waves.low$date), y = unique(dat.waves.low$loc), z = matrix(dat.waves.low$Estimate, nrow = length(unique(dat.waves.low$loc)), ncol = length(unique(dat.waves.low$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.waves.low$date), y = unique(dat.waves.low$loc), z = matrix(dat.waves.low$Estimate, nrow = length(unique(dat.waves.low$loc)), ncol = length(unique(dat.waves.low$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  draw.the.polygons(xlimits, ylimits)
  box(lwd = 1.2)
  draw.contour.labels(xlimits, ylimits, xfraction = 0.09, yfraction = 0.06465, x.corners = c(68, 63.7, 27.5), y.corners = c(0.4, 0.25, 0.3), contour.labels = c("0.4", "1.0", "1.4"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(low.wave.pools$date, low.wave.pools$loc, cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 2, line = 0.5, at = c(0.2,0.4,0.6,0.8,1), las = 1, cex = tick.label.size)
  mtext(text = "Pool location", side = 2, line = 1.7, at = 0.585, cex = axis.label.size)
  par(xpd = TRUE)
  mtext(text = c("Lake", "Forest"), side = 2, line = 1.65, las = 0, at = c(ylimits[1]+0.075, ylimits[2]-0.065), adj = c(0,1), cex = tick.label.size-0.05)
  arrows(y0 = ylimits[1]+0.07, x0 = 15.5, y1 = ylimits[1]-0.01, x1 = 15.5, angle = 30, code = 2, lwd = 1.0, length = 0.1 )
  arrows(y0 = ylimits[2]-0.06, x0 = 15.5, y1 = ylimits[2]+0.02, x1 = 15.5, angle = 30, code = 2, lwd = 1.0, length = 0.1 )
  par(xpd = FALSE)
  mtext(text = paste("C. Low wave height (", format.p.val(low.wave.value, 2), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  # Waves right side, panel D
par(new = "TRUE", plt = c(0.57, 0.92, 0.55, 0.72))
  filled.contour3(x = unique(dat.waves.hi$date), y = unique(dat.waves.hi$loc), z = matrix(dat.waves.hi$Estimate, nrow = length(unique(dat.waves.hi$loc)), ncol = length(unique(dat.waves.hi$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.waves.hi$date), y = unique(dat.waves.hi$loc), z = matrix(dat.waves.hi$Estimate, nrow = length(unique(dat.waves.hi$loc)), ncol = length(unique(dat.waves.hi$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", levels = c(0.5,1,1.5) )
  draw.the.polygons(xlimits, ylimits)
  box(lwd = 1.2)
  draw.contour.labels(xlimits, ylimits, xfraction = 0.09, yfraction = 0.06465, x.corners = c(34, 58, 56.3), y.corners = c(0.72, 0.90, 0.62), contour.labels = c("0.5", "1.0", "1.5"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.wave.pools$date, hi.wave.pools$loc, cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 2, line = 0.5, at = c(0.2,0.4,0.6,0.8,1), las = 1, cex = tick.label.size)
  mtext(text = paste("D. High wave height (", format.p.val(hi.wave.value, 2), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  #
  # Dragonflies left side, panel E
par(new = "TRUE", plt = c(0.13, 0.48, 0.33, 0.50))
  filled.contour3(x = unique(dat.dfly.low$date), y = unique(dat.dfly.low$loc), z = matrix(dat.dfly.low$Estimate, nrow = length(unique(dat.dfly.low$loc)), ncol = length(unique(dat.dfly.low$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dfly.low$date), y = unique(dat.dfly.low$loc), z = matrix(dat.dfly.low$Estimate, nrow = length(unique(dat.dfly.low$loc)), ncol = length(unique(dat.dfly.low$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  draw.the.polygons(xlimits, ylimits)
  box(lwd = 1.2)
  draw.contour.labels(xlimits, ylimits, xfraction = 0.09, yfraction = 0.06465, x.corners = c(52, 63, 57), y.corners = c(0.86, 0.6, 0.78), contour.labels = c("0.8", "1.0", "1.2"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(low.dfly.pools$date, low.dfly.pools$loc, cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 2, line = 0.5, at = c(0.2,0.4,0.6,0.8,1), las = 1, cex = tick.label.size)
  mtext(text = "Pool location", side = 2, line = 1.7, at = 0.585, cex = axis.label.size)
  par(xpd = TRUE)
  mtext(text = c("Lake", "Forest"), side = 2, line = 1.65, las = 0, at = c(ylimits[1]+0.075, ylimits[2]-0.065), adj = c(0,1), cex = tick.label.size-0.05)
  arrows(y0 = ylimits[1]+0.07, x0 = 15.5, y1 = ylimits[1]-0.01, x1 = 15.5, angle = 30, code = 2, lwd = 1.0, length = 0.1 )
  arrows(y0 = ylimits[2]-0.06, x0 = 15.5, y1 = ylimits[2]+0.02, x1 = 15.5, angle = 30, code = 2, lwd = 1.0, length = 0.1 )
  par(xpd = FALSE)
  mtext(paste("E. Low dragonfly abundance (", format.p.val(low.dflies.value, 0), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  # Dragonflies right side, panel F
par(new = "TRUE", plt = c(0.57, 0.92, 0.33, 0.50))
  filled.contour3(x = unique(dat.dfly.hi$date), y = unique(dat.dfly.hi$loc), z = matrix(dat.dfly.hi$Estimate, nrow = length(unique(dat.dfly.hi$loc)), ncol = length(unique(dat.dfly.hi$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dfly.hi$date), y = unique(dat.dfly.hi$loc), z = matrix(dat.dfly.hi$Estimate, nrow = length(unique(dat.dfly.hi$loc)), ncol = length(unique(dat.dfly.hi$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  draw.the.polygons(xlimits, ylimits)
  box(lwd = 1.2)
  draw.contour.labels(xlimits, ylimits, xfraction = 0.09, yfraction = 0.06465, x.corners = c(61, 25.3, 57), y.corners = c(0.47, 0.388, 0.35), contour.labels = c("0.5", "1.0", "1.5"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.dfly.pools$date, hi.dfly.pools$loc, cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 2, line = 0.5, at = c(0.2,0.4,0.6,0.8,1), las = 1, cex = tick.label.size)
  mtext(text = paste("F. High dragonfly abundance (", format.p.val(hi.dflies.value, 0), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  #
  # Tadpole density left side, panel G
par(new = "TRUE", plt = c(0.13, 0.48, 0.11, 0.28))
  filled.contour3(x = unique(dat.dens.low$date), y = unique(dat.dens.low$loc), z = matrix(dat.dens.low$Estimate, nrow = length(unique(dat.dens.low$loc)), ncol = length(unique(dat.dens.low$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dens.low$date), y = unique(dat.dens.low$loc), z = matrix(dat.dens.low$Estimate, nrow = length(unique(dat.dens.low$loc)), ncol = length(unique(dat.dens.low$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  draw.the.polygons(xlimits, ylimits)
  box(lwd = 1.2)
  draw.contour.labels(xlimits, ylimits, xfraction = 0.09, yfraction = 0.06465, x.corners = c(25, 59.5, 64), y.corners = c(0.78, 0.33, 0.565), contour.labels = c("1.0", "1.2", "1.4"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(low.dens.pools$date, low.dens.pools$loc, cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  mtext(text = "Hatching date", side = 1, line = 1.4, cex = axis.label.size)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 2, line = 0.5, at = c(0.2,0.4,0.6,0.8,1), las = 1, cex = tick.label.size)
  mtext(text = "Pool location", side = 2, line = 1.7, at = 0.585, cex = axis.label.size)
  mtext(text = paste("G. Low tadpole numbers (", format.p.val(low.density.value, 0), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  par(xpd = TRUE)
  mtext(text = c("Lake", "Forest"), side = 2, line = 1.65, las = 0, at = c(ylimits[1]+0.075, ylimits[2]-0.065), adj = c(0,1), cex = tick.label.size-0.05)
  arrows(y0 = ylimits[1]+0.07, x0 = 15.5, y1 = ylimits[1]-0.01, x1 = 15.5, angle = 30, code = 2, lwd = 1.0, length = 0.1 )
  arrows(y0 = ylimits[2]-0.06, x0 = 15.5, y1 = ylimits[2]+0.02, x1 = 15.5, angle = 30, code = 2, lwd = 1.0, length = 0.1 )
  par(xpd = FALSE)
  # Tadpole density right side, panel H
par(new = "TRUE", plt = c(0.57, 0.92, 0.11, 0.28))
  filled.contour3(x = unique(dat.dens.hi$date), y = unique(dat.dens.hi$loc), z = matrix(dat.dens.hi$Estimate, nrow = length(unique(dat.dens.hi$loc)), ncol = length(unique(dat.dens.hi$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = ylimits, yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dens.hi$date), y = unique(dat.dens.hi$loc), z = matrix(dat.dens.hi$Estimate, nrow = length(unique(dat.dens.hi$loc)), ncol = length(unique(dat.dens.hi$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  draw.the.polygons(xlimits, ylimits)
  box(lwd = 1.2)
  draw.contour.labels(xlimits, ylimits, xfraction = 0.09, yfraction = 0.06465, x.corners = c(27, 31.6, 43), y.corners = c(0.66, 0.3, 0.19), contour.labels = c("0.2", "0.8", "1.0"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.dens.pools$date, hi.dens.pools$loc, cex = point.size, pch = 21, bg = "red")
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  mtext(text = "Hatching date", side = 1, line = 1.4, cex = axis.label.size)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.2", "0.4", "0.6", "0.8", "1.0"), side = 2, line = 0.5, at = c(0.2,0.4,0.6,0.8,1), las = 1, cex = tick.label.size)
  mtext(text = paste("H. High tadpole numbers (", format.p.val(hi.density.value, 0), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
dev.off()
#
#
# Figure S8 -- How correlational selection on location/date is affected by ecological factors.
#
nr.extreme.years <- 3
fine.grid        <- 19
xlimits          <- c(23, 73.5)
ylimits          <- c(0.028, 20)
xlim.range       <- (xlimits - mean.date.indivs.kept) / sd.date.indivs.kept
ylim.range       <- (log(ylimits) - mean.larea.indivs.kept) / sd.larea.indivs.kept
d1               <- aggregate(d[,c("date", "larea", "prop.dry", "mean.wave.ht", "total.aeshna", "Nr.tads")], by = list(year = d$year, pool = d$pool), mean, na.rm = TRUE)
d1$area          <- exp(d1$larea)
#
# Drying
dry                 <- dry[order(dry$prop.dry),]
low.dry.pools       <- d1[d1$year %in% dry$year[1:6], ]                 # Set aside the pools that were used in drier and wetter years.
low.dry.pools       <- aggregate(low.dry.pools[,c("pool")], by = list(date = low.dry.pools$date, area = low.dry.pools$area), mean)
hi.dry.pools        <- d1[d1$year %in% dry$year[11:16], ]
hi.dry.pools        <- aggregate(hi.dry.pools[,c("pool")], by = list(date = hi.dry.pools$date, area = hi.dry.pools$area), mean)
low.dry.value       <- mean(dry$prop.dry[1:nr.extreme.years])           # low-drying years
low.dry.c           <- (low.dry.value - mean.propdry.kept) / sd.propdry.kept
hi.dry.value        <- mean(dry$prop.dry[(17-nr.extreme.years):16])     # high-drying years
hi.dry.c            <- (hi.dry.value - mean.propdry.kept) / sd.propdry.kept
newdat.low          <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), loc.s = 0, prop.dry = low.dry.value, propdry.c = low.dry.c)
newdat.low$loc.s2   <- 0
newdat.low$larea.s2 <- newdat.low$larea.s^2 / sd.larea.s2.indivs.kept
newdat.low$date.s2  <- newdat.low$date.s^2 / sd.date.s2.indivs.kept
newdat.low          <- cbind(newdat.low, Estimate = predict(m1.drying, type = "response", re.form = NA, newdata = newdat.low)/10)
newdat.low$area     <- exp((newdat.low$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
newdat.low$date     <- (newdat.low$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.hi           <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), loc.s = 0, prop.dry = hi.dry.value, propdry.c = hi.dry.c)
newdat.hi$loc.s2    <- 0
newdat.hi$larea.s2  <- newdat.hi$larea.s^2 / sd.larea.s2.indivs.kept
newdat.hi$date.s2   <- newdat.hi$date.s^2 / sd.date.s2.indivs.kept
newdat.hi           <- cbind(newdat.hi, Estimate = predict(m1.drying, type = "response", re.form = NA, newdata = newdat.hi)/10)
newdat.hi$area      <- exp((newdat.hi$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
newdat.hi$date      <- (newdat.hi$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
dat.dry.low         <- newdat.low
dat.dry.hi          <- newdat.hi
#
# Waves
waves               <- waves[order(waves$mean.wave.ht),]
low.wave.pools      <- d1[d1$year %in% waves$year[1:6], ]                     # Set aside the pools that were used in less wavey and wavier years.
low.wave.pools      <- aggregate(low.wave.pools[,c("pool")], by = list(date = low.wave.pools$date, area = low.wave.pools$area), mean)
hi.wave.pools       <- d1[d1$year %in% waves$year[11:16], ]
hi.wave.pools       <- aggregate(hi.wave.pools[,c("pool")], by = list(date = hi.wave.pools$date, area = hi.wave.pools$area), mean)
low.wave.value      <- mean(waves$mean.wave.ht[1:nr.extreme.years])           # low-wave years
low.wave.c          <- (low.wave.value - mean.meanwave.kept) / sd.meanwave.kept
hi.wave.value       <- mean(waves$mean.wave.ht[(17-nr.extreme.years):16])     # high-wave years
hi.wave.c           <- (hi.wave.value - mean.meanwave.kept) / sd.meanwave.kept
newdat.low          <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), loc.s = 0, meanwave = low.wave.value, meanwave.c = low.wave.c)
newdat.low$loc.s2   <- 0
newdat.low$larea.s2 <- newdat.low$larea.s^2 / sd.larea.s2.indivs.kept
newdat.low$date.s2  <- newdat.low$date.s^2 / sd.date.s2.indivs.kept
newdat.low          <- cbind(newdat.low, Estimate = predict(m1.mean.wave, type = "response", re.form = NA, newdata = newdat.low)/10)
newdat.low$date     <- (newdat.low$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.low$area     <- exp((newdat.low$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
newdat.hi           <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), loc.s = 0, meanwave = hi.wave.value, meanwave.c = hi.wave.c)
newdat.hi$loc.s2    <- 0
newdat.hi$larea.s2  <- newdat.hi$larea.s^2 / sd.larea.s2.indivs.kept
newdat.hi$date.s2   <- newdat.hi$date.s^2 / sd.date.s2.indivs.kept
newdat.hi           <- cbind(newdat.hi, Estimate = predict(m1.mean.wave, type = "response", re.form = NA, newdata = newdat.hi)/10)
newdat.hi$date      <- (newdat.hi$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.hi$area      <- exp((newdat.hi$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
dat.waves.low       <- newdat.low
dat.waves.hi        <- newdat.hi
#
# Dragonflies
aeshna             <- total.aeshna[total.aeshna$year < 1999 & total.aeshna$year > 1985, ]
aeshna             <- aeshna[order(aeshna$total.aeshna),]
low.dfly.pools     <- d1[d1$year %in% aeshna$year[1:6], ]                    # Set aside the pools that were used in low- and high-aeshna years.
low.dfly.pools     <- aggregate(low.dfly.pools[,c("pool")], by = list(date = low.dfly.pools$date, area = low.dfly.pools$area), mean)
hi.dfly.pools      <- d1[d1$year %in% aeshna$year[8:13], ]
hi.dfly.pools      <- aggregate(hi.dfly.pools[,c("pool")], by = list(date = hi.dfly.pools$date, area = hi.dfly.pools$area), mean)
low.dfly.value     <- mean(aeshna$total.aeshna[1:nr.extreme.years])           # low-dfly years
low.dfly.c         <- (log(low.dfly.value) - mean.aeshna.kept) / sd.aeshna.kept
hi.dfly.value      <- mean(aeshna$total.aeshna[(14-nr.extreme.years):13])     # high-dfly years
hi.dfly.c          <- (log(hi.dfly.value) - mean.aeshna.kept) / sd.aeshna.kept
newdat.low          <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), loc.s = 0, aeshna = low.dfly.value, aeshna.c = low.dfly.c)
newdat.low$loc.s2   <- 0
newdat.low$larea.s2 <- newdat.low$larea.s^2 / sd.larea.s2.indivs.kept
newdat.low$date.s2  <- newdat.low$date.s^2 / sd.date.s2.indivs.kept
newdat.low          <- cbind(newdat.low, Estimate = predict(m1.Aeshna, type = "response", re.form = NA, newdata = newdat.low)/10)
newdat.low$date     <- (newdat.low$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.low$area     <- exp((newdat.low$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
newdat.hi           <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), loc.s = 0, aeshna = hi.dfly.value, aeshna.c = hi.dfly.c)
newdat.hi$loc.s2    <- 0
newdat.hi$larea.s2  <- newdat.hi$larea.s^2 / sd.larea.s2.indivs.kept
newdat.hi$date.s2   <- newdat.hi$date.s^2 / sd.date.s2.indivs.kept
newdat.hi           <- cbind(newdat.hi, Estimate = predict(m1.Aeshna, type = "response", re.form = NA, newdata = newdat.hi)/10)
newdat.hi$date      <- (newdat.hi$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.hi$area      <- exp((newdat.hi$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
dat.dfly.low        <- newdat.low
dat.dfly.hi         <- newdat.hi
#
# Tadpole numbers
dens4               <- dens4[order(dens4$Nr.tads),]
low.dens.pools      <- d1[d1$year %in% dens4$year[1:6], ]                     # Set aside the pools that were used in low-density and high-density years.
low.dens.pools      <- aggregate(low.dens.pools[,c("pool")], by = list(date = low.dens.pools$date, area = low.dens.pools$area), mean)
hi.dens.pools       <- d1[d1$year %in% dens4$year[11:16], ]
hi.dens.pools       <- aggregate(hi.dens.pools[,c("pool")], by = list(date = hi.dens.pools$date, area = hi.dens.pools$area), mean)
low.dens.value      <- mean(dens4$Nr.tads[1:nr.extreme.years])               # low-density years
low.dens.c          <- (log(low.dens.value) - mean.Nrtads.kept) / sd.Nrtads.kept
hi.dens.value       <- mean(dens4$Nr.tads[(17-nr.extreme.years):16])         # high-density years
hi.dens.c           <- (log(hi.dens.value) - mean.Nrtads.kept) / sd.Nrtads.kept
newdat.low          <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), loc.s = 0, Nr.tads = low.dens.value, Nr.tads.c = low.dens.c)
newdat.low$loc.s2   <- 0
newdat.low$larea.s2 <- newdat.low$larea.s^2 / sd.larea.s2.indivs.kept
newdat.low$date.s2  <- newdat.low$date.s^2 / sd.date.s2.indivs.kept
newdat.low          <- cbind(newdat.low, Estimate = predict(m1.Nrtads, type = "response", re.form = NA, newdata = newdat.low)/10)
newdat.low$date     <- (newdat.low$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.low$area     <- exp((newdat.low$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
newdat.hi           <- expand.grid(pool = NA, date.s = set.up.interval(xlim.range, fine.grid), larea.s = set.up.interval(ylim.range, fine.grid), loc.s = 0, Nr.tads = hi.dens.value, Nr.tads.c = hi.dens.c)
newdat.hi$loc.s2    <- 0
newdat.hi$larea.s2  <- newdat.hi$larea.s^2 / sd.larea.s2.indivs.kept
newdat.hi$date.s2   <- newdat.hi$date.s^2 / sd.date.s2.indivs.kept
newdat.hi           <- cbind(newdat.hi, Estimate = predict(m1.Nrtads, type = "response", re.form = NA, newdata = newdat.hi)/10)
newdat.hi$date      <- (newdat.hi$date.s * sd.date.indivs.kept) + mean.date.indivs.kept
newdat.hi$area      <- exp((newdat.hi$larea.s * sd.larea.indivs.kept) + mean.larea.indivs.kept)
dat.dens.low        <- newdat.low
dat.dens.hi         <- newdat.hi
#
upper.cutoff <- 1.9
dat.dry.low$Estimate[dat.dry.low$Estimate > upper.cutoff]     <- upper.cutoff
dat.dry.hi$Estimate[dat.dry.hi$Estimate > upper.cutoff]       <- upper.cutoff
dat.waves.low$Estimate[dat.waves.low$Estimate > upper.cutoff] <- upper.cutoff
dat.waves.hi$Estimate[dat.waves.hi$Estimate > upper.cutoff]   <- upper.cutoff
dat.dfly.low$Estimate[dat.dfly.low$Estimate > upper.cutoff]   <- upper.cutoff
dat.dfly.hi$Estimate[dat.dfly.hi$Estimate > upper.cutoff]     <- upper.cutoff
dat.dens.low$Estimate[dat.dens.low$Estimate > upper.cutoff]   <- upper.cutoff
dat.dens.hi$Estimate[dat.dens.hi$Estimate > upper.cutoff]     <- upper.cutoff
#
point.size         <- 0.6
contour.label.size <- 0.85
tick.label.size    <- 0.9
axis.label.size    <- 1.2
draw.the.polygons  <- function(xlimits, ylimits) {
     polygon(x = c(xlimits[1], xlimits[1], 30.5), y = log(c(5.0, ylimits[2], ylimits[2])), border = NA, col = "white")
     polygon(x = c(55, xlimits[2], xlimits[2]), y = log(c(ylimits[2], ylimits[2], 2)), border = NA, col = "white")
     polygon(x = c(65.5, xlimits[2], xlimits[2]), y = log(c(ylimits[1], 0.08, ylimits[1])), border = NA, col = "white")
     polygon(x = c(xlimits[1], xlimits[1], 41), y = log(c(0.13, ylimits[1], ylimits[1])), border = NA, col = "white")
     }
color.gradient <- function(x, colors = c(lower.color, lowmid.color, highmid.color, upper.color), colsteps = 100) { return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )  }
plot.colors    <- color.gradient(1:100)
color.range    <- c(min(dat.dry.low$Estimate, dat.waves.low$Estimate, dat.dfly.low$Estimate, dat.dens.low$Estimate, dat.dry.hi$Estimate, dat.waves.hi$Estimate, dat.dfly.hi$Estimate, dat.dens.hi$Estimate), max(dat.dry.low$Estimate, dat.waves.low$Estimate, dat.dfly.low$Estimate, dat.dens.low$Estimate, dat.dry.hi$Estimate, dat.waves.hi$Estimate, dat.dfly.hi$Estimate, dat.dens.hi$Estimate))
#
setwd(figs.wd)
pdf("Fig S8.pdf", useDingbats = FALSE, width = 8, height = 14)
plot.new()
  # Panel A: Drying left side
par(new = "TRUE", plt = c(0.13, 0.48, 0.77, 0.94))
  filled.contour3(x = unique(dat.dry.low$date), y = unique(log(dat.dry.low$area)), z = matrix(dat.dry.low$Estimate, nrow = length(unique(dat.dry.low$area)), ncol = length(unique(dat.dry.low$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dry.low$date), y = unique(log(dat.dry.low$area)), z = matrix(dat.dry.low$Estimate, nrow = length(unique(dat.dry.low$area)), ncol = length(unique(dat.dry.low$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  draw.the.polygons(xlimits, ylimits)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.09, yfraction = 0.06465, x.corners = c(28, 31.3, 47), y.corners = log(c(0.8, 8, 10)), contour.labels = c("0.4", "0.6", "1.0"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(low.dry.pools$date, log(low.dry.pools$area), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = "Pool size (m  )", side = 2, line = 1.7, cex = axis.label.size)
  mtext(text = "2", side = 2, at = log(2.95), line = 2.1, cex = axis.label.size-0.3)
  mtext(text = paste("A. Low pool drying (", format.p.val(low.dry.value, 2), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  mtext(text = paste("Van Buskirk & Smith, Figure S8,", Sys.Date()), side = 3, line = 3, las = 1, cex = 0.8)
  # Panel B: Drying right side
par(new = "TRUE", plt = c(0.57, 0.92, 0.77, 0.94))
  filled.contour3(x = unique(dat.dry.hi$date), y = unique(log(dat.dry.hi$area)), z = matrix(dat.dry.hi$Estimate, nrow = length(unique(dat.dry.hi$area)), ncol = length(unique(dat.dry.hi$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dry.hi$date), y = unique(log(dat.dry.hi$area)), z = matrix(dat.dry.hi$Estimate, nrow = length(unique(dat.dry.hi$area)), ncol = length(unique(dat.dry.hi$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  draw.the.polygons(xlimits, ylimits)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.09, yfraction = 0.06465, x.corners = c(59, 31.2, 67), y.corners = log(c(0.9, 0.20, 0.175)), contour.labels = c("0.5", "1.0", "1.5"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.dry.pools$date, log(hi.dry.pools$area), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.015, labels = FALSE)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = paste("B. High pool drying (", format.p.val(hi.dry.value, 2), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  #
  # Panel C: Waves left side
par(new = "TRUE", plt = c(0.13, 0.48, 0.55, 0.72))
  filled.contour3(x = unique(dat.waves.low$date), y = unique(log(dat.waves.low$area)), z = matrix(dat.waves.low$Estimate, nrow = length(unique(dat.waves.low$area)), ncol = length(unique(dat.waves.low$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.waves.low$date), y = unique(log(dat.waves.low$area)), z = matrix(dat.waves.low$Estimate, nrow = length(unique(dat.waves.low$area)), ncol = length(unique(dat.waves.low$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  draw.the.polygons(xlimits, ylimits)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.09, yfraction = 0.06465, x.corners = c(65.5, 68, 63.5), y.corners = log(c(1.5, 0.38, 0.05)), contour.labels = c("0.4", "0.8", "1.4"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(low.wave.pools$date, log(low.wave.pools$area), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = "Pool size (m  )", side = 2, line = 1.7, cex = axis.label.size)
  mtext(text = "2", side = 2, at = log(2.95), line = 2.1, cex = axis.label.size-0.3)
  mtext(text = paste("C. Low wave height (", format.p.val(low.wave.value, 2), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  #
  # Panel D: Waves right side
par(new = "TRUE", plt = c(0.57, 0.92, 0.55, 0.72))
  filled.contour3(x = unique(dat.waves.hi$date), y = unique(log(dat.waves.hi$area)), z = matrix(dat.waves.hi$Estimate, nrow = length(unique(dat.waves.hi$area)), ncol = length(unique(dat.waves.hi$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.waves.hi$date), y = unique(log(dat.waves.hi$area)), z = matrix(dat.waves.hi$Estimate, nrow = length(unique(dat.waves.hi$area)), ncol = length(unique(dat.waves.hi$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  draw.the.polygons(xlimits, ylimits)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.09, yfraction = 0.06465, x.corners = c(30, 45, 63), y.corners = log(c(7, 9, 1.3)), contour.labels = c("0.4", "0.6", "0.8"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.wave.pools$date, log(hi.wave.pools$area), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = paste("D. High wave height (", format.p.val(hi.wave.value, 2), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  #
  # Panel E: Dragonflies left side
par(new = "TRUE", plt = c(0.13, 0.48, 0.33, 0.50))
  filled.contour3(x = unique(dat.dfly.low$date), y = unique(log(dat.dfly.low$area)), z = matrix(dat.dfly.low$Estimate, nrow = length(unique(dat.dfly.low$area)), ncol = length(unique(dat.dfly.low$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dfly.low$date), y = unique(log(dat.dfly.low$area)), z = matrix(dat.dfly.low$Estimate, nrow = length(unique(dat.dfly.low$area)), ncol = length(unique(dat.dfly.low$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 8 )
  draw.the.polygons(xlimits, ylimits)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.09, yfraction = 0.06465, x.corners = c(24.5, 55.0, 66), y.corners = log(c(0.7, 0.076, 1.1)), contour.labels = c("0.6", "0.8", "1.2"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(low.dfly.pools$date, log(low.dfly.pools$area), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = "Pool size (m  )", side = 2, line = 1.7, cex = axis.label.size)
  mtext(text = "2", side = 2, at = log(2.95), line = 2.1, cex = axis.label.size-0.3)
  mtext(paste("E. Low dragonfly abundance (", format.p.val(low.dfly.value, 0), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  #
  # Panel F: Dragonflies right side
par(new = "TRUE", plt = c(0.57, 0.92, 0.33, 0.50))
  filled.contour3(x = unique(dat.dfly.hi$date), y = unique(log(dat.dfly.hi$area)), z = matrix(dat.dfly.hi$Estimate, nrow = length(unique(dat.dfly.hi$area)), ncol = length(unique(dat.dfly.hi$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dfly.hi$date), y = unique(log(dat.dfly.hi$area)), z = matrix(dat.dfly.hi$Estimate, nrow = length(unique(dat.dfly.hi$area)), ncol = length(unique(dat.dfly.hi$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  draw.the.polygons(xlimits, ylimits)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.09, yfraction = 0.06465, x.corners = c(59, 60.7, 63.5), y.corners = log(c(0.95, 0.24, 0.09)), contour.labels = c("0.5", "1.0", "1.5"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.dfly.pools$date, log(hi.dfly.pools$area), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = paste("F. High dragonfly abundance (", format.p.val(hi.dfly.value, 0), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  #
  # Panel G: Tadpole density left side
par(new = "TRUE", plt = c(0.13, 0.48, 0.11, 0.28))
  filled.contour3(x = unique(dat.dens.low$date), y = unique(log(dat.dens.low$area)), z = matrix(dat.dens.low$Estimate, nrow = length(unique(dat.dens.low$area)), ncol = length(unique(dat.dens.low$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dens.low$date), y = unique(log(dat.dens.low$area)), z = matrix(dat.dens.low$Estimate, nrow = length(unique(dat.dens.low$area)), ncol = length(unique(dat.dens.low$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  draw.the.polygons(xlimits, ylimits)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.09, yfraction = 0.06465, x.corners = c(36.5, 45, 60), y.corners = log(c(8, 8, 0.12)), contour.labels = c("0.5", "1.0", "1.5"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(low.dens.pools$date, log(low.dens.pools$area), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  mtext(text = "Hatching date", side = 1, line = 1.4, cex = axis.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = "Pool size (m  )", side = 2, line = 1.7, cex = axis.label.size)
  mtext(text = "2", side = 2, at = log(2.95), line = 2.1, cex = axis.label.size-0.3)
  mtext(text = paste("G. Low tadpole numbers (", format.p.val(low.dens.value, 0), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
  # Panel H: Tadpole density right side
par(new = "TRUE", plt = c(0.57, 0.92, 0.11, 0.28))
  filled.contour3(x = unique(dat.dens.hi$date), y = unique(log(dat.dens.hi$area)), z = matrix(dat.dens.hi$Estimate, nrow = length(unique(dat.dens.hi$area)), ncol = length(unique(dat.dens.hi$date)), byrow = TRUE), levels = seq(color.range[1], color.range[2], length = 101), col = plot.colors, axes = FALSE, frame.plot = TRUE, ylab = "", ylim = log(ylimits), yaxs="i", xlim = xlimits, xaxs="i")
  contour(x = unique(dat.dens.hi$date), y = unique(log(dat.dens.hi$area)), z = matrix(dat.dens.hi$Estimate, nrow = length(unique(dat.dens.hi$area)), ncol = length(unique(dat.dens.hi$date)), byrow = TRUE), ylab = "", drawlabels = FALSE, add = TRUE, col = "black", nlevels = 6 )
  draw.the.polygons(xlimits, ylimits)
  draw.contour.labels(xlimits, log(ylimits), xfraction = 0.09, yfraction = 0.06465, x.corners = c(63, 51, 61.25), y.corners = log(c(0.95, 0.035, 0.035)), contour.labels = c("0.4", "1.0", "1.2"), label.bg = "white", text.size = contour.label.size, draw.box = TRUE, label.line.width = 0.8)
  points(hi.dens.pools$date, log(hi.dens.pools$area), cex = point.size, pch = 21, bg = "red")
  box(lwd = 1.2)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(text = c( paste(unlist(make.dates(30)), collapse = " "), paste(unlist(make.dates(50)), collapse = " "), paste(unlist(make.dates(70)), collapse = " ") ), side = 1, line = 0.3, las = 1, at = c(30, 50, 70), cex = tick.label.size)
  mtext(text = "Hatching date", side = 1, line = 1.4, cex = axis.label.size)
  axis(side = 2, tck = -0.025, at = log(c(0.1,1,10)), labels = FALSE)
  axis(side = 2, tck = -0.01, at = log(c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,20)), labels = FALSE)
  mtext(text = c("0.1", "1", "10"), side = 2, line = 0.5, las = 1, at = log(c(0.1,1,10)), cex = tick.label.size)
  mtext(text = paste("H. High tadpole numbers (", format.p.val(hi.dens.value, 0), ")", sep = ""), side = 3, line = 0.2, las = 1, at = xlimits[1], adj = 0, cex = axis.label.size-0.1)
dev.off()
#
