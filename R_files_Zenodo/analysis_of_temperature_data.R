# Analysis of Judith Bachmann's 2015 transplant experiment
#
# HOBO data logger temperature data
#
# Start by loading the functions at the bottom of this script.
#
# Set working directory to the location where all the files are.
#
setwd(" ... your working directory here ...")




#############################
#
# Start with some other environmental variables.
# These data are just manually entered here from Table S1 in the MS.
# These tests were reported in "Results" under "Environmental conditions".
#
env <- as.data.frame(matrix(c(491, 25, 7.9, 8.1, 0.13,
         558, 3300, 7.25, 1.32, 0.14,
         648, 488, 7.87, 8.77, 0.17,
         519, 1006, 7.23, 6, 0.15,
         2342, 185, 9.27, 10.01, 0.04,
         2004, 64, 6.87, 7.5, 0,
         2388, 262, 7.1, 7.98, 0.1,
         2557, 152, 7.37, 7.84, 0.06), nrow = 8, ncol = 5, byrow = TRUE))
names(env)    <- c("elev", "clutches", "pH", "O2", "sal")
env$pond      <- c("ellw", "feer", "munt", "siec", "bern", "bide", "flue", "muot")
env$elev.fact <- ifelse(env$elev > 2000, "high", "low")
m1 <- lm(env$clutches ~ env$elev)
m2 <- lm(env$pH ~ env$elev)
m3 <- lm(env$O2 ~ env$elev)
m4 <- lm(env$sal ~ env$elev)
m5 <- lm(env$clutches ~ env$elev.fact)
m6 <- lm(env$pH ~ env$elev.fact)
m7 <- lm(env$O2 ~ env$elev.fact)
m8 <- lm(env$sal ~ env$elev.fact)
summary(m1)
summary(m2)
summary(m3)
summary(m4)   # salinity is the only variable that changes significantly with elevation.
summary(m5)
summary(m6)
summary(m7)
summary(m8)   # salinity is the only variable that varies significantly between elevations.
par(mfrow = c(2, 2))
plot(env$elev, env$clutches)
plot(env$elev, env$pH)
plot(env$elev, env$O2)
plot(env$elev, env$sal)
dev.off()
#
#
#
# Import the dates of the experiment
#
dates             <- read.table("inventory_transplant_loggers.txt", stringsAsFactors = FALSE, header = TRUE)
dates$start.Jdate <- month.day.to.Julian(dates$start.month, dates$start.day)
dates$end.Jdate   <- month.day.to.Julian(dates$end.month, dates$end.day)
dates             <- dates[ order(dates$name), c(1,2,3,8,9)]
names(dates)[3]   <- "pond.block"
dates
#
# Get a list of all the temperature files, figure out which ones actually contain
# temperature data, and recover ID information for each file.
#
filelist   <- list.files()
data.files <- filelist[grep("csv", filelist)]
#
# Get the pond and block for each file. 
#
ID.1       <- read.table(textConnection(data.files), sep="_", col.names=c("pond.block", "B"))
n.files    <- dim(ID.1)[1]
ID.1$pond  <- substr(ID.1$pond.block, 1, 4)
ID.1$block <- substr(ID.1$pond.block, 5, 5)
ID         <- ID.1[ , c(3,4,1)]
#
# Import all the files, put them together into a data frame, calculate Julian date and
# time, clean up the data frame.
#
collatedfiles <- c()
numberrows    <- c()
file.count    <- 0
for (i in data.files) {
  file.count    <- file.count + 1
  f1            <- read.csv(i, colClasses="character", header=TRUE)
  f1            <- f1[, -c(1,2)]
  if(dim(f1)[2] > 2) f1 <- f1[ , c(1,2,3)]
  ID.1          <- matrix(unlist(strsplit(f1$date, split=c("-"))), ncol = 3, nrow = dim(f1)[1], byrow = TRUE)
  f1$year       <- as.numeric(ID.1[,1])
  f1$month      <- as.numeric(ID.1[,2])
  if(i == "flue7_8009.csv") {
    ID.2          <- matrix(unlist(strsplit(ID.1[,3], split=c(":"))), ncol = 2, nrow = dim(f1)[1], byrow = TRUE)
    } else {
    ID.2          <- matrix(unlist(strsplit(ID.1[,3], split=c(":"))), ncol = 3, nrow = dim(f1)[1], byrow = TRUE)
    }
  f1$min        <- as.numeric(ID.2[,2])
  ID.2[,1]      <- sub("  ", " ", ID.2[,1])   # in one file, there were two spaces instead of 1 between two entries
  ID.3          <- matrix(unlist(strsplit(ID.2[,1], split=c(" "))), ncol = 2, nrow = dim(f1)[1], byrow = TRUE)
  f1$hour       <- as.numeric(ID.3[,2])
  f1$day        <- as.numeric(ID.3[,1])
  f1$time       <- f1$hour + f1$min/60
  if(i == "flue7_8009.csv") {
    f1$time <- ifelse(f1$am.pm == "PM", f1$time + 12, f1$time)
    }
  f1$Jdate      <- month.day.to.Julian(f1$month, f1$day) + f1$time/24
  f1$pond       <- ID$pond[file.count]
  f1$block      <- as.numeric(ID$block[file.count])
  f1$pond.block <- ID$pond.block[file.count]
  f1$temp       <- as.numeric(f1$temp)
  collatedfiles <- rbind(collatedfiles, f1[,c("pond", "block", "pond.block", "Jdate", "temp")])
  }
d <- collatedfiles
d$day = floor(d$Jdate)
head(d, n = 15)
dim(d)
#
#
# Edit the data block-by-block to include only temperature data that
# were recorded while the datalogger was in the water.
#
d$temp <- ifelse(d$pond.block == "bern1" & d$Jdate < 177.4458, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern1" & d$Jdate > 239.3625, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern2" & d$Jdate < 177.4792, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern2" & d$Jdate > 239.3542, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern3" & d$Jdate < 177.5243, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern3" & d$Jdate > 239.3576, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern4" & d$Jdate < 177.4868, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern4" & d$Jdate > 239.3618, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern5" & d$Jdate < 177.4993, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern5" & d$Jdate > 218.8326, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern6" & d$Jdate < 177.5000, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern6" & d$Jdate > 239.3750, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern7" & d$Jdate < 177.4924, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern7" & d$Jdate > 239.3674, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern8" & d$Jdate < 177.4951, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bern8" & d$Jdate > 224.5368, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bide1" & d$Jdate < 176.7215, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bide1" & d$Jdate > 211.5340, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bide2" & d$Jdate < 176.6500, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bide2" & d$Jdate > 211.7125, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bide3" & d$Jdate < 176.6562, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bide3" & d$Jdate > 211.6979, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bide4" & d$Jdate < 176.6576, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bide4" & d$Jdate > 211.6368, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bide5" & d$Jdate < 176.6778, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bide6" & d$Jdate < 176.7028, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bide6" & d$Jdate > 217.4528, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bide7" & d$Jdate < 176.7083, NA, d$temp)
d$temp <- ifelse(d$pond.block == "bide7" & d$Jdate > 211.7083, NA, d$temp)
d$temp <- ifelse(d$pond.block == "ellw1" & d$Jdate < 100.6042, NA, d$temp)
d$temp <- ifelse(d$pond.block == "ellw1" & d$Jdate > 155.6042, NA, d$temp)
d$temp <- ifelse(d$pond.block == "ellw2" & d$Jdate < 100.62014, NA, d$temp)
d$temp <- ifelse(d$pond.block == "ellw2" & d$Jdate > 155.57847, NA, d$temp)
d$temp <- ifelse(d$pond.block == "ellw3" & d$Jdate < 100.62014, NA, d$temp)
d$temp <- ifelse(d$pond.block == "ellw3" & d$Jdate > 155.57847, NA, d$temp)
d$temp <- ifelse(d$pond.block == "ellw4" & d$Jdate < 100.52153, NA, d$temp)
d$temp <- ifelse(d$pond.block == "ellw4" & d$Jdate > 145.58403, NA, d$temp)
d$temp <- ifelse(d$pond.block == "ellw7" & d$Jdate < 100.64097, NA, d$temp)
d$temp <- ifelse(d$pond.block == "ellw7" & d$Jdate > 148.41181, NA, d$temp)
d$temp <- ifelse(d$pond.block == "ellw8" & d$Jdate < 100.65764, NA, d$temp)
d$temp <- ifelse(d$pond.block == "ellw8" & d$Jdate > 155.69931, NA, d$temp)
d$temp <- ifelse(d$pond.block == "feer1" & d$Jdate < 100.54444, NA, d$temp)
d$temp <- ifelse(d$pond.block == "feer1" & d$Jdate > 145.54444, NA, d$temp)
d$temp <- ifelse(d$pond.block == "feer2" & d$Jdate < 100.48681, NA, d$temp)
d$temp <- ifelse(d$pond.block == "feer3" & d$Jdate < 100.51042, NA, d$temp)
d$temp <- ifelse(d$pond.block == "feer3" & d$Jdate > 145.55208, NA, d$temp)
d$temp <- ifelse(d$pond.block == "feer4" & d$Jdate < 115.0194, NA, d$temp)
d$temp <- ifelse(d$pond.block == "feer4" & d$Jdate > 152.68611, NA, d$temp)
d$temp <- ifelse(d$pond.block == "feer5" & d$Jdate < 100.62222, NA, d$temp)
d$temp <- ifelse(d$pond.block == "feer5" & d$Jdate > 145.70556, NA, d$temp)
d$temp <- ifelse(d$pond.block == "feer6" & d$Jdate < 100.64653, NA, d$temp)
d$temp <- ifelse(d$pond.block == "feer6" & d$Jdate > 145.58403, NA, d$temp)
d$temp <- ifelse(d$pond.block == "feer7" & d$Jdate < 100.54375, NA, d$temp)
d$temp <- ifelse(d$pond.block == "feer7" & d$Jdate > 145.60625, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue2" & d$Jdate < 176.8806, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue2" & d$Jdate > 239.5681, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue3" & d$Jdate < 176.8736, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue3" & d$Jdate > 239.5403, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue4" & d$Jdate < 176.8528, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue4" & d$Jdate > 239.5194, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue5" & d$Jdate < 176.9007, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue5" & d$Jdate > 239.5049, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue6" & d$Jdate < 176.8611, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue6" & d$Jdate > 239.5486, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue7" & d$Jdate < 176.8972, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue7" & d$Jdate > 239.5431, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue8" & d$Jdate < 176.8979, NA, d$temp)
d$temp <- ifelse(d$pond.block == "flue8" & d$Jdate > 239.5229, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt1" & d$Jdate < 100.76250, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt1" & d$Jdate > 159.63750, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt2" & d$Jdate < 100.78333, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt2" & d$Jdate > 159.65833, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt3" & d$Jdate < 100.76389, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt3" & d$Jdate > 159.63889, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt4" & d$Jdate < 100.74792, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt4" & d$Jdate > 155.45625, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt4" & d$Jdate > 108.8521 & d$Jdate < 113.6854, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt5" & d$Jdate < 100.74583, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt5" & d$Jdate > 159.60000, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt6" & d$Jdate < 100.7500, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt6" & d$Jdate > 159.5833, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt7" & d$Jdate < 100.74306, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt7" & d$Jdate > 159.55556, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt8" & d$Jdate < 100.75556, NA, d$temp)
d$temp <- ifelse(d$pond.block == "munt8" & d$Jdate > 155.44306, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot1" & d$Jdate < 177.7222, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot1" & d$Jdate > 226.7222, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot2" & d$Jdate < 177.7250, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot2" & d$Jdate > 218.4542, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot3" & d$Jdate < 177.7403, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot3" & d$Jdate > 218.4486, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot4" & d$Jdate < 177.7354, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot4" & d$Jdate > 226.7354, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot5" & d$Jdate < 177.7424, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot5" & d$Jdate > 218.4299, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot6" & d$Jdate < 177.7424, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot6" & d$Jdate > 218.4507, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot7" & d$Jdate < 177.7181, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot7" & d$Jdate > 226.6764, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot8" & d$Jdate < 177.7535, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot8" & d$Jdate > 207, NA, d$temp)
d$temp <- ifelse(d$pond.block == "muot8" & d$Jdate > 218.4201, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec1" & d$Jdate < 100.8125, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec1" & d$Jdate > 168.4375, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec2" & d$Jdate < 100.81319, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec2" & d$Jdate > 168.43819, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec3" & d$Jdate < 100.80556, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec3" & d$Jdate > 168.40972, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec4" & d$Jdate < 100.82083, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec4" & d$Jdate > 168.42500, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec5" & d$Jdate < 100.84236, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec5" & d$Jdate > 168.38403, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec6" & d$Jdate < 100.84514, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec6" & d$Jdate > 168.47014, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec7" & d$Jdate < 100.85000, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec7" & d$Jdate > 168.47500, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec8" & d$Jdate < 100.85694, NA, d$temp)
d$temp <- ifelse(d$pond.block == "siec8" & d$Jdate > 168.48194, NA, d$temp)
# 
#
##################
#
# If you want, DELETE all data after 5 weeks from the start of the block (35 days).
#
list.of.blocks     <- unique(d$pond.block)
dates$highest.date <- dates$lowest.date <- NA
for(i in 1:length(list.of.blocks)) {
  d.sub     <- subset(d, d$pond.block == list.of.blocks[i])
  which.row <- match(list.of.blocks[i], dates$pond.block)
  dates[which.row, c("lowest.date","highest.date")] <- range(d.sub$Jdate)
  }
dates$date.5.weeks <- dates$lowest.date + 35
d <- merge(d, dates[, c("pond.block", "date.5.weeks")], by = "pond.block")
d <- subset(d, d$Jdate < d$date.5.weeks)
#
#
# Round all readings from all ponds to 48 half-hour time intervals.
#
d$Jdate1 <- round((d$Jdate * 48)*1)
d$Jdate1 <- round(d$Jdate1/(48 * 1), 4)
ss <- strsplit(as.character(d$Jdate1), "[.]")
right <- numeric()
for(i in 1:length(ss)) {
  sss <- unlist(ss[i])
  if(length(sss) == 1) right[i] <- 0
  if(length(sss) == 2) right[i] <- as.numeric(paste(".", sss[2], sep = ""))
  }
#
#
# One part of feer4 is wildly wrong for unknown reasons.
#
d$temp <- ifelse(d$pond.block == "feer4" & d$Jdate > 107.6444 & d$Jdate < 115.0194, NA, d$temp)
#
# Interpolate data for feer4 based on other nearby blocks in Feerbach.
#
d.sub <- subset(d, d$pond == "feer")
feer  <- data.frame(Jdate1 = unique(d.sub$Jdate1)[order(unique(d.sub$Jdate1))])
feer  <- merge(feer, subset(d.sub, d.sub$pond.block == "feer4")[ , c("Jdate", "Jdate1", "temp")], by = "Jdate1")
feer  <- merge(feer, subset(d.sub, d.sub$pond.block == "feer1")[ , c("Jdate1", "temp")], by = "Jdate1")
feer  <- merge(feer, subset(d.sub, d.sub$pond.block == "feer2")[ , c("Jdate1", "temp")], by = "Jdate1")
names(feer) <- c("Jdate1", "Jdate", "feer4", "feer1", "feer2")
# The rounding did not work properly for block feer3. Fix that.
d.sub$Jdate1 <- ifelse(d.sub$pond.block == "feer3", round((d.sub$Jdate+0.0001) * 48), d.sub$Jdate1)
d.sub$Jdate1 <- ifelse(d.sub$pond.block == "feer3", round(d.sub$Jdate1/(48 * 1), 4), d.sub$Jdate1)
feer <- merge(feer, subset(d.sub, d.sub$pond.block == "feer3")[ , c("Jdate1", "temp")], by = "Jdate1")
feer <- merge(feer, subset(d.sub, d.sub$pond.block == "feer5")[ , c("Jdate1", "temp")], by = "Jdate1")
feer <- merge(feer, subset(d.sub, d.sub$pond.block == "feer6")[ , c("Jdate1", "temp")], by = "Jdate1")
names(feer)[6:8] <- c("feer3", "feer5", "feer6")
feer <- merge(feer, subset(d.sub, d.sub$pond.block == "feer7")[ , c("Jdate1", "temp")], by = "Jdate1")
names(feer)[9] <- c("feer7")
feer <- feer[complete.cases(feer[,-3]), ]
MS.dev <- function(x, y) {
  SS <- sum((x-y)^2, na.rm = TRUE)
  N  <- samp.size(x)[1]
  SS/N
  }
res <- data.frame(other.block = c("feer1", "feer2", "feer3", "feer5", "feer6", "feer7"), MS.dev = NA)
for(i in 1:6) {
  res$MS.dev[i] <- MS.dev(feer[,"feer4"], feer[,(3+i)])
  }
res     # This result indicates that blocks 1-5 are the closest to feer4, overall.
        # Could use the average of those four blocks to fill in feer4.
# feer$feer4 <- ifelse(is.na(feer$feer4), mean(c(feer$feer1,feer$feer2,feer$feer3,feer$feer5)), feer$feer4)
#         Better is to predict feer4 based on other blocks.
pred.feer4 <- lm(feer4 ~ feer1 + feer5 + feer6 + feer7, data=feer)
cc         <- coefficients(pred.feer4)
feer$feer4 <- ifelse(is.na(feer$feer4), cc[1] + cc[2]*feer$feer1 + cc[3]*feer$feer5 + cc[4]*feer$feer6 + cc[5]*feer$feer7, feer$feer4)
# Transfer the predicted values for feer4 over to dataset d.
for(i in 1:dim(d)[1]) {
  if(d$pond.block[i] == "feer4" & is.na(d$temp[i])) {
    for(j in 1:dim(feer)[1]) {
      if(feer$Jdate[j] == d$Jdate[i]) {
        d$temp[i] <- feer[j, "feer4"]
        } } } }
d <- d[complete.cases(d), ]
#
#
# Calculate the mean, min, and max for each day. Infer whether the skies were sunny or cloudy.
#
daily.mean <- aggregate(d$temp, by = list(d$pond, d$block, d$pond.block, d$day), mean, na.rm = TRUE)
daily.min  <- aggregate(d$temp, by = list(d$pond, d$block, d$pond.block, d$day), min, na.rm = TRUE)
daily.max  <- aggregate(d$temp, by = list(d$pond, d$block, d$pond.block, d$day), max, na.rm = TRUE)
daily.N    <- aggregate(d$temp, by = list(d$pond, d$block, d$pond.block, d$day), length)
t.data     <- merge(daily.mean, daily.min, by = c("Group.1", "Group.2", "Group.3", "Group.4"))
t.data     <- merge(t.data, daily.max, by = c("Group.1", "Group.2", "Group.3", "Group.4"))
t.data     <- merge(t.data, daily.N, by = c("Group.1", "Group.2", "Group.3", "Group.4"))
colnames(t.data) <- c("pond", "block", "pond.block", "day", "mean", "min", "max", "N")
t.data           <- subset(t.data, t.data$N > 40)    # delete the date if too many readings are missing
t.data$weather   <- ifelse((t.data$max - t.data$min) > 8, 'sun', 'unk')
t.data$weather   <- ifelse((t.data$max - t.data$min) < 5, 'clo', t.data$weather)
head(t.data)
#
#
# Prepare data for Table S1 in the MS.
#
block.means <- aggregate(t.data[,c("mean", "min", "max")], by = list(t.data$pond, t.data$pond.block, t.data$weather), mean)
block.means <- block.means[ block.means$Group.3 != "unk", ]
names(block.means)[1:3] <- c("pond", "pond.block", "weather")
pond.means  <- aggregate(block.means[,c("mean", "min", "max")], by = list(block.means$pond, block.means$weather), mean)
pond.sd     <- aggregate(block.means[,c("mean", "min", "max")], by = list(block.means$pond, block.means$weather), sd)
names(pond.sd)[3:5] <- c("mean.sd", "min.sd", "max.sd")
pond.N      <- aggregate(block.means[,c("mean", "min", "max")], by = list(block.means$pond, block.means$weather), length)
names(pond.N)[3] <- c("N")
pond.means  <- merge(pond.means, pond.N[,1:3], by = c("Group.1", "Group.2"))
pond.means  <- merge(pond.means, pond.sd, by = c("Group.1", "Group.2"))
pond.means$mean.se <- pond.means$mean.sd / sqrt(pond.means$N)
pond.means$min.se  <- pond.means$min.sd / sqrt(pond.means$N)
pond.means$max.se  <- pond.means$max.sd / sqrt(pond.means$N)
pond.means  <- pond.means[ , c(1,2,6,3,10,4,11,5,12)]
names(pond.means)[1:2] <- c("pond", "weather")
Sunny       <- subset(pond.means, pond.means$weather == "sun")
Sunny       <- Sunny[,c(1,3,6,7,8,9)]
names(Sunny)[2:6] <- c("N.SUN", "min.SUN", "min.se.SUN", "max.SUN", "max.se.SUN")
Cloudy      <- subset(pond.means, pond.means$weather == "clo")
Cloudy      <- Cloudy[,c(1,3,6,7,8,9)]
names(Cloudy)[2:6] <- c("N.CLO", "min.CLO", "min.se.CLO", "max.CLO", "max.se.CLO")
pond.means  <- merge(Sunny, Cloudy, by = c("pond"), all = TRUE)
pH.data <- as.data.frame(matrix(c(7.90, 0.25, 8.1, 0.21, 0.13, 0,
             7.23, 0.07, 6.0, 0.67, 0.15, 0,
             7.25, 0.05, 1.32,0.24, 0.14, 0,
             7.87, 0.15, 8.77, 1.31, 0.17, 0,
             6.87, 0.38, 7.50, 0.49,0.00, 0,
             9.27, 0.23, 10.01, 0.04, 0.04, 0,
             7.10, 0.00, 7.98, 0.08, 0.10, 0,
             7.37, 0.22, 7.84, 0.09, 0.06, 0), nrow = 8, ncol = 6, byrow = TRUE))
pH.data$pond <- c("ellw","siec","feer","munt","bide","bern","flue","muot")
names(pH.data)[1:6] <- c("pH", "pH.se", "O2", "O2.se", "salinity", "salinity.se")
pond.means <- merge(pond.means, pH.data, by = "pond")[c(3,4,6,8,1,2,5,7), ]
pond.means    # This is something like Table S1 in the MS.
#
#
#  Summarize mean, min, max on sunny and cloudy days
#
list.of.ponds <- unique(d$pond)
w             <- data.frame()
for(i in 1:length(list.of.ponds)) {
  d.sub     <- subset(t.data, t.data$pond == list.of.ponds[i])
  weat      <- weather.each.day(d.sub)[ , c("day", "weather")]
  weat$pond <- list.of.ponds[i]
  w         <- rbind(w, weat)
  }
temp.data     <- merge(t.data[ , c("pond", "day", "mean", "min", "max")], w, by = c("pond", "day"))
pond.means    <- aggregate(temp.data[,c("mean","min","max")], by = list(temp.data$pond, temp.data$weather), mean, na.rm = TRUE)
pond.means    <- merge(pond.means, pond.N[,1:3], by = c("Group.1", "Group.2"))
names(pond.means)[1:2] <- c("pond", "weather")
pond.means <- pond.means[pond.means$weather != "unk",]
pond.means$elev <- c("high", "low", "low", "high", "low", "low", "high", "high", "low", "low", "high", "low", "high", "low")
elev.means <- aggregate(pond.means[,c("mean","min","max")], by = list(pond.means$elev, pond.means$weather), mean)
elev.sds   <- aggregate(pond.means[,c("mean","min","max")], by = list(pond.means$elev, pond.means$weather), sd)
elev.N     <- aggregate(pond.means[,c("mean","min","max")], by = list(pond.means$elev, pond.means$weather), length)
names(elev.N)[3]     <- c("N")
names(elev.sds)[3:5] <- c("mean.sd", "min.sd", "max.sd")
elev.means <- merge(elev.means, elev.N[,1:3], by = c("Group.1", "Group.2"))
elev.means <- merge(elev.means, elev.sds, by = c("Group.1", "Group.2"))
elev.means$mean.se <- elev.means$mean.sd/sqrt(elev.means$N)
elev.means$min.se  <- elev.means$min.sd/sqrt(elev.means$N)
elev.means$max.se  <- elev.means$max.sd/sqrt(elev.means$N)
elev.means <- elev.means[ c(3,4,1,2), c(1,2,6,3,10,7,4,11,8,5,12,9)]
names(elev.means)[1:2] <- c("elev", "weather")
elev.means
#
#
#
#################################################################
#
# Now draw Fig. S7 -- the daily time course of temperature.
#
ss <- strsplit(as.character(d$Jdate1), "[.]")
right <- numeric()
for(i in 1:length(ss)) {
  sss <- unlist(ss[i])
  if(length(sss) == 1) right[i] <- 0
  if(length(sss) == 2) right[i] <- as.numeric(paste(".", sss[2], sep = ""))
  }
d$right.side <- right
d1        <- merge(d[,c("pond.block", "pond", "day", "right.side", "temp")], t.data[,c("pond.block", "pond", "day", "weather")], by = c("pond.block", "pond", "day"))
d2        <- aggregate(d1$temp, by = list(d1$pond.block, d1$pond, d1$right.side, d1$weather), mean)
names(d2) <- c("pond.block", "pond", "right.side", "weather", "temp")
d3        <- aggregate(d2$temp, by = list(d2$pond, d2$right.side, d2$weather), mean)
names(d3) <- c("pond", "right.side", "weather", "temp")
d3$elev   <- ifelse(d3$pond == "bern", "high", "low")
d3$elev   <- ifelse(d3$pond == "bide", "high", d3$elev)
d3$elev   <- ifelse(d3$pond == "flue", "high", d3$elev)
d3$elev   <- ifelse(d3$pond == "muot", "high", d3$elev)
d4        <- aggregate(d3$temp, by = list(d3$elev, d3$right.side, d3$weather), mean)
names(d4) <- c("elev", "right.side", "weather", "temp")
#
low.color    <- "#FEA100"       # an orange color
high.color   <- "#1A0778"       # a dark blue
#
pdf(" .... your path .... /Fig S7.pdf", width = 6.5, height = 8, useDingbats = FALSE)
plot.new()
  # Upper panel :: time course through the day
  par(new = "TRUE", plt = c(0.15, 0.85, 0.67, 0.92))
  x.ticks     <- c(0.3, 1.4)
  x.offset    <- 0.015
  point.size  <- 2
  line.width  <- 2.5
  error.width <- 2
  elev.means$weat <- ifelse(elev.means$weather == "clo", x.ticks[1], x.ticks[2])
  y.limits   <- c(10,26)
  x.limits   <- c(0,2)
  plot(x = x.ticks, y = y.limits, pch = "", xlim = x.limits, ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n")
  box(lwd = 1.3)
  axis(side = 1, labels = FALSE, at = x.ticks, tck = -0.015)
  mtext(text = c("Cloudy","Clear"), side = 1, line = 0.4, at = x.ticks, cex = 0.9)
  mtext(text = c("Weather"), side = 1, line = 1.6, las = 1, cex = 1.2)
  axis(side = 2, labels = FALSE, tck = -0.015)
  mtext(text = c("10","15","20","25"), side = 2, line = 0.6, las = 1, at = c(10,15,20,25), cex = 0.9)
  mtext(text = c("Temperature (°C)"), side = 2, line = 2.4, las = 3, cex = 1.2)
  # Lay down the points, error bars, and lines
  points(x = x.ticks-x.offset, y = elev.means[elev.means$elev == "low", "min"], pch = 21, lwd = 2.5, col = low.color, cex = point.size)
  lines(x = c(x.ticks[1], x.ticks[1])-x.offset, y = c(elev.means[elev.means$elev == "low" & elev.means$weather == "clo", "min"] - elev.means[elev.means$elev == "low" & elev.means$weather == "clo", "min.se"], elev.means[elev.means$elev == "low" & elev.means$weather == "clo", "min"] + elev.means[elev.means$elev == "low" & elev.means$weather == "clo", "min.se"]), lwd = error.width, col = low.color)
  lines(x = c(x.ticks[2], x.ticks[2])-x.offset, y = c(elev.means[elev.means$elev == "low" & elev.means$weather == "sun", "min"] - elev.means[elev.means$elev == "low" & elev.means$weather == "sun", "min.se"], elev.means[elev.means$elev == "low" & elev.means$weather == "sun", "min"] + elev.means[elev.means$elev == "low" & elev.means$weather == "sun", "min.se"]), lwd = error.width, col = low.color)
  lines(x = x.ticks-x.offset, y = elev.means[elev.means$elev == "low", "min"], lty = 2, lwd = line.width, col = low.color)
  points(x = x.ticks-x.offset, y = elev.means[elev.means$elev == "low", "max"], pch = 19, col = low.color, cex = point.size)
  lines(x = c(x.ticks[1], x.ticks[1])-x.offset, y = c(elev.means[elev.means$elev == "low" & elev.means$weather == "clo", "max"] - elev.means[elev.means$elev == "low" & elev.means$weather == "clo", "max.se"], elev.means[elev.means$elev == "low" & elev.means$weather == "clo", "max"] + elev.means[elev.means$elev == "low" & elev.means$weather == "clo", "max.se"]), lwd = error.width, col = low.color)
  lines(x = c(x.ticks[2], x.ticks[2])-x.offset, y = c(elev.means[elev.means$elev == "low" & elev.means$weather == "sun", "max"] - elev.means[elev.means$elev == "low" & elev.means$weather == "sun", "max.se"], elev.means[elev.means$elev == "low" & elev.means$weather == "sun", "max"] + elev.means[elev.means$elev == "low" & elev.means$weather == "sun", "max.se"]), lwd = error.width, col = low.color)
  lines(x = x.ticks-x.offset, y = elev.means[elev.means$elev == "low", "max"], lty = 1, lwd = line.width, col = low.color)
  points(x = x.ticks+x.offset, y = elev.means[elev.means$elev == "high", "min"], pch = 21, lwd = 2.5, col = high.color, cex = point.size)
  lines(x = c(x.ticks[1], x.ticks[1])+x.offset, y = c(elev.means[elev.means$elev == "high" & elev.means$weather == "clo", "min"] - elev.means[elev.means$elev == "high" & elev.means$weather == "clo", "min.se"], elev.means[elev.means$elev == "high" & elev.means$weather == "clo", "min"] + elev.means[elev.means$elev == "high" & elev.means$weather == "clo", "min.se"]), lwd = error.width, col = high.color)
  lines(x = c(x.ticks[2], x.ticks[2])+x.offset, y = c(elev.means[elev.means$elev == "high" & elev.means$weather == "sun", "min"] - elev.means[elev.means$elev == "high" & elev.means$weather == "sun", "min.se"], elev.means[elev.means$elev == "high" & elev.means$weather == "sun", "min"] + elev.means[elev.means$elev == "high" & elev.means$weather == "sun", "min.se"]), lwd = error.width, col = high.color)
  lines(x = x.ticks+x.offset, y = elev.means[elev.means$elev == "high", "min"], lty = 2, lwd = line.width, col = high.color)
  points(x = x.ticks+x.offset, y = elev.means[elev.means$elev == "high", "max"], pch = 19, col = high.color, cex = point.size)
  lines(x = c(x.ticks[1], x.ticks[1])+x.offset, y = c(elev.means[elev.means$elev == "high" & elev.means$weather == "clo", "max"] - elev.means[elev.means$elev == "high" & elev.means$weather == "clo", "max.se"], elev.means[elev.means$elev == "high" & elev.means$weather == "clo", "max"] + elev.means[elev.means$elev == "high" & elev.means$weather == "clo", "max.se"]), lwd = error.width, col = high.color)
  lines(x = c(x.ticks[2], x.ticks[2])+x.offset, y = c(elev.means[elev.means$elev == "high" & elev.means$weather == "sun", "max"] - elev.means[elev.means$elev == "high" & elev.means$weather == "sun", "max.se"], elev.means[elev.means$elev == "high" & elev.means$weather == "sun", "max"] + elev.means[elev.means$elev == "high" & elev.means$weather == "sun", "max.se"]), lwd = error.width, col = high.color)
  lines(x = x.ticks+x.offset, y = elev.means[elev.means$elev == "high", "max"], lty = 1, lwd = line.width, col = high.color)
  text(c("Maximum", "temperatures", "Minimum", "temperatures"), x = x.ticks[2] + 0.05, y = c(24.7,23,12.5,10.8), pos = 4, cex = 0.95)
  text("A", x = 0.04, y = 24.7, cex = 1.8)
  mtext(text = c("Fig. S7", "made by 'analysis of Judith temperature data.R'"), at = c(0.1, x.ticks[2]), side = 3, line = 2.1, cex = c(1,0.6))
  #
  # Lower panel :: time course through the day
  par(new = "TRUE", plt = c(0.15, 0.85, 0.11, 0.57))
  d4        <- subset(d4, d4$right.side > 0)   # cloudy and sunnny mis-classified for records made at midnight
  y.limits  <- range(d4$temp)
  x.limits  <- range(d4$right.side)
  plot(d4$right.side, d4$temp, pch = "", xlim = x.limits, ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n")
  box(lwd = 1.3)
  axis(side = 1, labels = FALSE, at = c(0,0.25,0.5,0.75,1), tck = -0.015)
  mtext(text = c("00:00","06:00","12:00","18:00","24:00"), side = 1, line = 0.4, at = c(0,0.25,0.5,0.75,1), cex = 0.9)
  mtext(text = c("Time of day"), side = 1, line = 1.65, las = 1, cex = 1.2)
  axis(side = 2, labels = FALSE, tck = -0.015)
  mtext(text = c("12","16","20","24"), side = 2, line = 0.6, las = 1, at = c(12,16,20,24), cex = 0.9)
  mtext(text = c("Temperature (°C)"), side = 2, line = 2.4, las = 3, cex = 1.2)
  d.sub <- subset(d4, d4$elev == "low" & d4$weather == "sun")
  lines(d.sub$right.side, d.sub$temp, lwd = 4, col = low.color)
  d.sub <- subset(d4, d4$elev == "low" & d4$weather == "clo")
  lines(d.sub$right.side, d.sub$temp, lwd = 4, lty = 2, col = low.color)
  d.sub <- subset(d4, d4$elev == "high" & d4$weather == "sun")
  lines(d.sub$right.side, d.sub$temp, lwd = 4, col = high.color)
  d.sub <- subset(d4, d4$elev == "high" & d4$weather == "clo")
  lines(d.sub$right.side, d.sub$temp, lwd = 4, lty = 2, col = high.color)
  text(c("High", "elevation", "Low", "elevation"), x = c(0.43, 0.43, 0.7, 0.7), y = c(22.7,22,15.5,14.8), cex = 1.0, col = c(high.color, high.color, low.color, low.color))
  text("B", x = 0.04, y = 23.5, cex = 1.8)
  #
dev.off()







##########################################################
#
# Load these functions first.
# 
#
weather.each.day <- function(d.sub) {
  ss        <- table(d.sub$day, d.sub$weather)
  ss.df     <- as.data.frame(table(d.sub$day, d.sub$weather))
  sun.days  <- ss.df[ ss.df$Var2 == "sun", ]$Freq
  clo.days  <- ss.df[ ss.df$Var2 == "clo", ]$Freq
  if(length(clo.days) > 1) sun.cloud <- data.frame(day = as.numeric(attr(ss, "dimnames")[[1]]), nr.loggers = apply(ss, 1, sum), nr.cloudy = clo.days)
  if(length(sun.days) > 1) sun.cloud <- data.frame(day = as.numeric(attr(ss, "dimnames")[[1]]), nr.loggers = apply(ss, 1, sum), nr.sunny = sun.days)
  if(length(clo.days) > 1 & length(sun.days) > 1) sun.cloud <- data.frame(day = as.numeric(attr(ss, "dimnames")[[1]]), nr.loggers = apply(ss, 1, sum), nr.sunny = sun.days, nr.cloudy = clo.days)
  if(length(sun.days) > 1) sun.cloud$weather <- ifelse((sun.cloud$nr.sunny / sun.cloud$nr.loggers) > 0.5, "sun", "unk")
  if(length(clo.days) > 1) sun.cloud$weather <- ifelse((sun.cloud$nr.cloudy / sun.cloud$nr.loggers) > 0.5, "clo", sun.cloud$weather)
  sun.cloud
  }
#
##########################################################
samp.size <- function(x) {
  # Calculate sample sizes, NAs and non-NAs seperately.
  n = length(x) - sum(is.na(x))
  nas = sum(is.na(x))
  out = c(n, nas)
  names(out) = c("", "NAs")
  out
  }
##########################################################
month.day.to.Julian <- function(month, day, year = 1981) {
  # Function to calulate the Julian date of the year, given a month and day.
  # If no year is given we assume this is not a leap year.
  # If the year is given, the function accounts for whether it is a leap year.
  # Returns a vector of integers.
  #
  # This function is located in "2013-05-27 functions for Julian dates.R"
  #
  leap.years <- c(1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020, 2024, 2028)
  julian <- day
  julian <- ifelse (month>1,  (day+31), day)
  julian <- ifelse (month>2,  (day+59), julian)
  julian <- ifelse (month>3,  (day+90), julian)
  julian <- ifelse (month>4,  (day+120), julian)
  julian <- ifelse (month>5,  (day+151), julian)
  julian <- ifelse (month>6,  (day+181), julian)
  julian <- ifelse (month>7,  (day+212), julian)
  julian <- ifelse (month>8,  (day+243), julian)
  julian <- ifelse (month>9,  (day+273), julian)
  julian <- ifelse (month>10, (day+304), julian)
  julian <- ifelse (month>11, (day+334), julian)
  julian <- ifelse(year %in% leap.years & month > 2, julian + 1, julian)
  return(julian)
  }


