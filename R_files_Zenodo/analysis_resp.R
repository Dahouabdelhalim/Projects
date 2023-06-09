#### INITIAL DATA PREP ------------------------------------------------------------------------

library(lme4)
library(multcomp)
library(car)
library(lubridate)
library(magrittr)
library(viridis)

# Import main and hi-res respirometric datasets
varclasses <- c(rep("factor",4),rep("numeric",4))
resp <- read.delim("data_resp_full_190109.txt", colClasses=varclasses)
hiresp <- read.delim("data_hiresp_full_181003.txt", colClasses=varclasses)

# Fix up variables for both datasets
resp$weight <- resp$weight/10
resp$fill.time <- paste(substr(resp$fill.time,1,2),substr(resp$fill.time,3,4),"00",sep=":")
resp$push.time <- paste(substr(resp$push.time,1,2),substr(resp$push.time,3,4),"00",sep=":")
hiresp$weight <- hiresp$weight/10
hiresp$fill.time <- paste(substr(hiresp$fill.time,1,2),substr(hiresp$fill.time,3,4),"00",sep=":")
hiresp$push.time <- paste(substr(hiresp$push.time,1,2),substr(hiresp$push.time,3,4),"00",sep=":")

# Fill and push times in lubridate time format
resp$push <- ymd_hms(paste(resp$push.date, resp$push.time, sep=" "))
resp$fill <- ymd_hms(paste(date(resp$push.date)-1, resp$fill.time, sep=" "))
hiresp$push <- ymd_hms(paste(hiresp$push.date, hiresp$push.time, sep=" "))
hiresp$fill <- ymd_hms(paste(date(hiresp$push.date)-1, hiresp$fill.time, sep=" "))

# Time spent respiring in syringe
resp$syr.time <- as.numeric(resp$push-resp$fill)
hiresp$syr.time <- as.numeric(hiresp$push-hiresp$fill)

# Days since pupation:
resp$day <- rep(NA, nrow(resp))
for (i in levels(resp$ID)) {resp[resp$ID==i,]$day <- date(resp[resp$ID==i,]$push.date)-min(date(resp[resp$ID==i,]$push.date))+1 }
hiresp$day <- rep(NA, nrow(hiresp))
for (i in levels(hiresp$ID)) {hiresp[hiresp$ID==i,]$day <- date(hiresp[hiresp$ID==i,]$push.date)-min(date(hiresp[hiresp$ID==i,]$push.date))+1 }

# Treatment: D(irect), 0, 7 or r(ef)
resp$tmt <- as.factor(substr(resp$ID,1,1))
hiresp$tmt <- as.factor(substr(hiresp$ID,1,1))

# For main-series inds., the first measurement in the hi-res treatment is actually their second measurement:
hiresp[hiresp$tmt=="0",]$day <- hiresp[hiresp$tmt=="0",]$day + 2




#### CORRECTING FOR REFS; MERGING DATASETS ----------------------------------------------------------------------------------

shiftby1 <- function(x) {
  xnew <- rep(0, length(x))
  xnew[(2):length(x)] <- x[1:(length(x)-1)]
  xnew }

# Combine measurements between refs into ref-batches
resp$refbatch <- shiftby1(cumsum(resp$ref))
hiresp$refbatch <- shiftby1(cumsum(hiresp$ref))

# Check that this worked
par(mfrow=c(2,1))
plot(resp$refbatch, pch=19, col=as.numeric(resp$ref)+1)
plot(hiresp$refbatch, pch=19, col=as.numeric(hiresp$ref)+1)
par(mfrow=c(1,1))

# Subtract last CO2 and O2 value of each refbatch (i.e. the ref) from each measurement in the batch
for (i in unique(resp$refbatch)) {
  resp[resp$refbatch==i,]$co2.raw <- resp[resp$refbatch==i,]$co2.raw - resp[resp$refbatch==i,][nrow(resp[resp$refbatch==i,]),]$co2.raw
  resp[resp$refbatch==i,]$o2.raw <- resp[resp$refbatch==i,]$o2.raw - resp[resp$refbatch==i,][nrow(resp[resp$refbatch==i,]),]$o2.raw}
for (i in unique(hiresp$refbatch)) {
  hiresp[hiresp$refbatch==i,]$co2.raw <- hiresp[hiresp$refbatch==i,]$co2.raw - hiresp[hiresp$refbatch==i,][nrow(hiresp[hiresp$refbatch==i,]),]$co2.raw
  hiresp[hiresp$refbatch==i,]$o2.raw <- hiresp[hiresp$refbatch==i,]$o2.raw - hiresp[hiresp$refbatch==i,][nrow(hiresp[hiresp$refbatch==i,]),]$o2.raw}

# Remove refs from dataset; also remove the few cold-treatment inds. that were measured from the hi-res dataset
resp <- droplevels(subset(resp, tmt!="r"))
hiresp <- droplevels(subset(hiresp, tmt!="7" & tmt!="r"))

# Copy 15-day measurements from hi-res dataset to main dataset, and vice versa for 1-day measurements
missingno <- c("0006","0007","0008","0020","0021","0022","0031","0047","0048","0061","0068","0082")
resp <- droplevels(rbind( resp, subset(hiresp, ID %in% missingno & day==15) ))
hiresp <- droplevels(rbind( hiresp, subset(resp, ID %in% missingno & day==1) ))
# Both datasets are now complete.

# Individuals 0070 and 7087 were obviously dead after day 120 and 90, respectively, so remove those measurements
resp <- subset(resp, ID!="0070" | day<150 )
resp <- subset(resp, ID!="7087" | day<100 )

# Trim unnecessary columns
resp <- resp[,c(2,5,7:13)]
hiresp <- hiresp[,c(2,5,7:13)]

# Sort datasets by first treatment, then ID, then measurement time (facilitates figure drawing and maintainted sanity)
resp <- resp[with(resp, order(tmt, ID, day)),]
hiresp <- hiresp[with(hiresp, order(tmt, ID, day)),]
rownames(resp) <- NULL
rownames(hiresp) <- NULL

# Importing main protocol for extracting individual population, sex, etc.
reftable <- read.delim("temp_resp_ind_list.txt", colClasses=c("factor","factor","factor","factor","character") )
termin <- read.delim("data_termin_main_190301.txt",
                     colClasses = c("character","factor","factor","factor","factor","factor","numeric","character",
                                    "numeric","character","numeric","character","numeric","factor","factor","character"))

resp$pop <- rep(NA, nrow(resp))
resp$sex <- rep(NA, nrow(resp))
resp$family <- rep(NA, nrow(resp))
resp$status <- rep(NA, nrow(resp))
for (i in levels(resp$ID)) { resp[resp$ID==i,]$pop <- paste(termin[termin$id==i,]$pop) }
for (i in levels(resp$ID)) { resp[resp$ID==i,]$sex <- paste(termin[termin$id==i,]$sex) }
for (i in levels(resp$ID)) { resp[resp$ID==i,]$family <- paste(termin[termin$id==i,]$family) }
for (i in levels(resp$ID)) { resp[resp$ID==i,]$status <- paste(termin[termin$id==i,]$status) }

# Subtract push date from eclosion date to get the time from first measurement to eclosion
resp$ecl.date <- rep(NA, nrow(resp))
for (i in levels(resp$ID)) { resp[resp$ID==i,]$ecl.date <- paste(termin[termin$id==i,]$ecl.date) }
resp[resp$ecl.date=="",]$ecl.date <- NA
for (i in levels(resp$ID)) { resp[resp$ID==i,]$ecl.date <- date(resp[resp$ID==i,]$ecl.date) - min(date(resp[resp$ID==i,]$push))+1 }
resp$ecl.date <- as.numeric(resp$ecl.date)

# Do the same for the nondiapausing cohort
dirdates <- read.delim("nondiapause_times_190314.txt")
hiresp$ecl.date <- rep(NA, nrow(hiresp))
for (i in levels(hiresp$ID)) { hiresp[hiresp$ID==i & hiresp$tmt=="D",]$ecl.date <- paste(dirdates[dirdates$ID==i,]$ad.date) }
for (i in levels(hiresp$ID)) { hiresp[hiresp$ID==i,]$ecl.date <- date(hiresp[hiresp$ID==i,]$ecl.date) - min(date(hiresp[hiresp$ID==i,]$push))+1 }
hiresp$ecl.date <- as.numeric(hiresp$ecl.date)
hiresp[hiresp$ID=="0082",]$ecl.date <- resp[resp$ID=="0082",]$ecl.date[1]

# Convenience variable for approx. days in cold: round "day" to nearest multiple of 15, except for the 220-day batch
resp$roughday <- round(resp$day/15)*15
resp[resp$roughday==225 | resp$roughday==210,]$roughday <- 220




#### METABOLIC CALCULATIONS ---------------------------------------------------------------------------------------------

# Volume and time correction (filled 10 ml; pushed 8 ml; 10/8=1.25)
resp$co2.corr <- (resp$co2.raw*1.25)/resp$syr.time
resp$o2.corr <- (resp$o2.raw*1.25)/resp$syr.time
hiresp$co2.corr <- (hiresp$co2.raw*1.25)/hiresp$syr.time
hiresp$o2.corr <- (hiresp$o2.raw*1.25)/hiresp$syr.time

# Respiratory quotient
resp$RQ <- resp$co2.raw/resp$o2.raw
hiresp$RQ <- hiresp$co2.raw/hiresp$o2.raw

# Gross metabolic rate (µW)
resp$MR <- (((16+(5.164*resp$RQ))*resp$o2.corr)/3600)*1000000
hiresp$MR <- (((16+(5.164*hiresp$RQ))*hiresp$o2.corr)/3600)*1000000

# Net metabolic rate (µW/mg)
resp$MR.mass <- resp$MR/resp$weight
hiresp$MR.mass <- hiresp$MR/hiresp$weight

# Roughly integrate under metabolic curve for all surviving cold-treatment individuals
respsums <- data.frame(ID=NA, tmt=NA, pop=NA, family=NA, sex=NA, weight=NA, joules=NA)

for (i in levels(resp$ID)) {
  thisind <- subset(resp, tmt=="7" & ID==i & ID!="7087") # Exclude 7087, which died partway
  
  if (nrow(thisind)>1) {
    temp <- data.frame( thisind$MR, thisind$day )
    temp$time <- 0
    temp$integ <- 0
    
    for (t in 2:nrow(temp)) { 
      temp[t,]$time <- temp[t,2] - temp[t-1,2]
      temp[t,]$integ <- mean(c(temp[t,1], temp[t-1,1])) * temp[t,3] * 60*60*24*(1/1000000)
    }
    
    respsums <- rbind(respsums, c(i, as.character(thisind$tmt[1]), thisind$pop[1], thisind$family[1],
                                  thisind$sex[1], thisind$weight[1], sum(temp$integ)))  } }

# Fix up results table
respsums <- respsums[2:nrow(respsums),]
respsums$joules <- as.numeric(respsums$joules)
respsums$weight <- as.numeric(respsums$weight)



#### TOTAL ENERGY USE & STATISTICAL ANALYSIS ---------------------------------------------------------------------------------------------

# Analyse population effect
joulemod <- lmer(joules ~ weight+sex+pop + (1|family), data=respsums)
Anova(joulemod)
# Re-run without sex effect (non-significant)
joulemod1 <- lmer(joules ~ weight+pop + (1|family), data=respsums)
Anova(joulemod1)
comparisonsJoule<-glht(joulemod1, mcp(pop="Tukey"))
summary(comparisonsJoule) # Pairwise differences

coef <- fixef(joulemod1)




#### GENERATING FIGURES ---------------------------------------------------------------------------------------------

# FIGURE 4 PREP
nice.diapause <- subset(hiresp, tmt=="0" & ID!="0082") # Remove number 0082, which failed to enter diapause properly
nice.direct <- subset(hiresp, tmt=="D")
dircols <- magma(2, begin=0.2, end=0.5)

resp.cold <- list(R=subset(resp, tmt=="7" & pop=="R"),
                  V=subset(resp, tmt=="7" & pop=="V"),
                  B=subset(resp, tmt=="7" & pop=="B"),
                  L=subset(resp, tmt=="7" & pop=="L"))
resp.warm <- list(R=subset(resp, tmt=="0" & pop=="R"),
                  V=subset(resp, tmt=="0" & pop=="V"),
                  B=subset(resp, tmt=="0" & pop=="B"),
                  L=subset(resp, tmt=="0" & pop=="L"))


# Exporting figure 4
pdf(file="Fig4.pdf", pointsize=20)

# 4A: Metabolic rate per group over time
means.diap <- tapply(nice.diapause$MR.mass, nice.diapause$day, mean)
means.dir <- tapply(nice.direct$MR.mass, nice.direct$day, mean)
sds.diap <- tapply(nice.diapause$MR.mass, nice.diapause$day, sd)
sds.dir <- tapply(nice.direct$MR.mass, nice.direct$day, sd)
plot(hiresp$day, hiresp$MR.mass, col="white", ylab="Metabolic rate (µW/mg)", xlab="Days since pupation", xlim=c(0,17), bty="n")
# Bars show standard deviation:
arrows(x0=c(1,3,5,7,9,11,13,15,17),y0=means.diap+sds.diap,
       x1=c(1,3,5,7,9,11,13,15,17),y1=means.diap-sds.diap,
       code=3, angle=90, length=0.05, col=dircols[1])
arrows(x0=c(1,3,5,7,9,11,13,15),y0=means.dir+sds.dir,
       x1=c(1,3,5,7,9,11,13,15),y1=means.dir-sds.dir,
       code=3, angle=90, length=0.05, col=dircols[2])
points(means.diap~names(means.diap), col=dircols[1], pch=19, type="b")
points(means.dir~names(means.dir), col=dircols[2], pch=17, type="b")

legend("topleft", inset=0.01, legend=c("Diapausing","Nondiapausing"), lty="solid", col=dircols, pch=c(19,17))

# 4B: RQ per group over time
means.diap <- tapply(nice.diapause$RQ, nice.diapause$day, mean)
means.dir <- tapply(nice.direct$RQ, nice.direct$day, mean)
sds.diap <- tapply(nice.diapause$RQ, nice.diapause$day, sd)
sds.dir <- tapply(nice.direct$RQ, nice.direct$day, sd)
plot(hiresp$day, hiresp$RQ, col="white", ylab="Respiratory Quotient", xlab="Days since pupation", xlim=c(0,17), bty="n")
abline(h=1, lty="dotted", col="grey")
# Bars show standard deviation:
arrows(x0=c(1,3,5,7,9,11,13,15,17),y0=means.diap+sds.diap,
       x1=c(1,3,5,7,9,11,13,15,17),y1=means.diap-sds.diap,
       code=3, angle=90, length=0.05, col=dircols[1])
arrows(x0=c(1,3,5,7,9,11,13,15),y0=means.dir+sds.dir,
       x1=c(1,3,5,7,9,11,13,15),y1=means.dir-sds.dir,
       code=3, angle=90, length=0.05, col=dircols[2])
points(means.diap~names(means.diap), col=dircols[1], pch=19, type="b")
points(means.dir~names(means.dir), col=dircols[2], pch=17, type="b")

# 4C: MR per population over time (cold treatment)
se <- function(x) sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)]))
means.cold <- lapply(resp.cold, function(pop) { tapply(pop$MR.mass, pop$roughday, mean) } )
ses.cold <- lapply(resp.cold, function(pop) { tapply(pop$MR.mass, pop$roughday, se) } )
plot(MR.mass~roughday, data=subset(resp,tmt=="7"), col="white", ylab="Metabolic rate (µW/mg)", xlab="Days since pupation", bty="n")
points(means.cold$R~names(means.cold$R), type="b", pch=21, col="blue", bg="blue")
points(means.cold$V~names(means.cold$V), type="b", pch=22, col="cornflowerblue", bg="cornflowerblue")
points(means.cold$B~names(means.cold$B), type="b", pch=25, col="orange", bg="orange")
points(means.cold$L~names(means.cold$L), type="b", pch=24, col="red", bg="red")
# Bars show standard error:
arrows(x0=c(0,15,30,60,90,120,150,180,220),y0=means.cold$R+ses.cold$R,
       x1=c(0,15,30,60,90,120,150,180,220),y1=means.cold$R-ses.cold$R,
       code=3, angle=90, length=0.05, col="blue")
arrows(x0=c(0,15,30,60,90,120,150,180,220),y0=means.cold$V+ses.cold$V,
       x1=c(0,15,30,60,90,120,150,180,220),y1=means.cold$V-ses.cold$V,
       code=3, angle=90, length=0.05, col="cornflowerblue")
arrows(x0=c(0,15,30,60,90,120,150,180,220),y0=means.cold$B+ses.cold$B,
       x1=c(0,15,30,60,90,120,150,180,220),y1=means.cold$B-ses.cold$B,
       code=3, angle=90, length=0.05, col="orange")
arrows(x0=c(0,15,30,60,90,120,150,180,220),y0=means.cold$L+ses.cold$L,
       x1=c(0,15,30,60,90,120,150,180,220),y1=means.cold$L-ses.cold$L,
       code=3, angle=90, length=0.05, col="red")
legend("topright", legend=c("Stockholm","Highlands","Gotland","Skåne"), col=c("blue","cornflowerblue","orange","red"),
       lty=rep("solid",4), pch=15:18, inset=0.02)

# 4D: RQ per population over time (cold treatment)
se <- function(x) sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)]))
means.cold <- lapply(resp.cold, function(pop) { tapply(pop$RQ, pop$roughday, mean) } )
ses.cold <- lapply(resp.cold, function(pop) { tapply(pop$RQ, pop$roughday, se) } )
plot(RQ~roughday, data=subset(resp,tmt=="7"), col="white", ylab="Respiratory Quotient", xlab="Days since pupation", bty="n")
points(means.cold$R~names(means.cold$R), type="b", pch=21, col="blue", bg="blue")
points(means.cold$V~names(means.cold$V), type="b", pch=22, col="cornflowerblue", bg="cornflowerblue")
points(means.cold$B~names(means.cold$B), type="b", pch=25, col="orange", bg="orange")
points(means.cold$L~names(means.cold$L), type="b", pch=24, col="red", bg="red")
abline(h=1, lty="dotted", col="grey")
# Bars show standard error:
arrows(x0=c(0,15,30,60,90,120,150,180,220),y0=means.cold$R+ses.cold$R,
       x1=c(0,15,30,60,90,120,150,180,220),y1=means.cold$R-ses.cold$R,
       code=3, angle=90, length=0.05, col="blue")
arrows(x0=c(0,15,30,60,90,120,150,180,220),y0=means.cold$V+ses.cold$V,
       x1=c(0,15,30,60,90,120,150,180,220),y1=means.cold$V-ses.cold$V,
       code=3, angle=90, length=0.05, col="cornflowerblue")
arrows(x0=c(0,15,30,60,90,120,150,180,220),y0=means.cold$B+ses.cold$B,
       x1=c(0,15,30,60,90,120,150,180,220),y1=means.cold$B-ses.cold$B,
       code=3, angle=90, length=0.05, col="orange")
arrows(x0=c(0,15,30,60,90,120,150,180,220),y0=means.cold$L+ses.cold$L,
       x1=c(0,15,30,60,90,120,150,180,220),y1=means.cold$L-ses.cold$L,
       code=3, angle=90, length=0.05, col="red")
legend("topright", legend=c("Stockholm","Highlands","Gotland","Skåne"), col=c("blue","cornflowerblue","orange","red"),
       pt.bg=c("blue","cornflowerblue","orange","red"), lty=rep("solid",4), pch=c(21,22,25,24), inset=0.02)

# 4E: Individual trajectories of warm-treatment inds.
plot(MR.mass~day, data=subset(resp, tmt=="0"), col="white", ylab="Metabolic rate (µW/mg)", xlab="Days since pupation", bty="n")
warminds <- levels(droplevels(subset(resp, tmt=="0"))$ID)
for (i in warminds) {
  thisind <- subset(resp, ID==i)
  if("R" %in% thisind$pop) {
    points(MR.mass~day, data=thisind, type="b", pch=21, col="blue", bg="blue") }
  if("V" %in% thisind$pop) {
    points(MR.mass~day, data=thisind, type="b", pch=22, col="cornflowerblue", bg="cornflowerblue") }
  if("B" %in% thisind$pop) {
    points(MR.mass~day, data=thisind, type="b", pch=25, col="orange", bg="orange") }
  if("L" %in% thisind$pop) {
    points(MR.mass~day, data=thisind, type="b", pch=24, col="red", bg="red") }
  lastind <- subset(thisind, day==max(thisind$day))
  if (lastind$status=="ok") points(lastind$ecl.date, lastind$MR.mass, pch=8)
  if (lastind$status=="fe") points(lastind$ecl.date, lastind$MR.mass, pch=4)
  segments(lastind$day, lastind$MR.mass, lastind$ecl.date, lastind$MR.mass, lty="dotted") }

# 4F: Total winter energy use (cold treatment)
plotmod <- lmer(joules ~ weight+pop + (1|family), data=respsums)
coef <- fixef(plotmod)
plot(joules ~ weight, data=respsums, pch="", ylab="Winter energy expenditure (J)", xlab="Pupal weight", ylim=c(74,140), bty="n")
abline(coef[1], coef[2], col="orange", lty="dotted", lwd=1)
abline(coef[1]+coef[3], coef[2], col="red", lty="dotdash", lwd=1)
abline(coef[1]+coef[4], coef[2], col="blue", lty="solid", lwd=1)
abline(coef[1]+coef[5], coef[2], col="cornflowerblue", lty="dashed", lwd=1)
points(joules~weight,subset(respsums, pop=="B"), col="white", bg="orange", pch=25)
points(joules~weight,subset(respsums, pop=="L"), col="white", bg="red", pch=24)
points(joules~weight,subset(respsums, pop=="V"), col="white", bg="cornflowerblue", pch=22)
points(joules~weight,subset(respsums, pop=="R"), col="white", bg="blue", pch=21)
legend("topleft", inset=0.02, legend=c("Stockholm","Highlands","Gotland","Skåne"), col=c("blue","cornflowerblue","orange","red"),
       pt.bg=c("blue","cornflowerblue","orange","red"),  lty=c("solid","dashed","dotted","dotdash"), pch=c(21,22,25,24))


dev.off()



