#### SETTING UP -----------------------------------------------------------------------------

library(viridis)
library(lubridate)
library(magrittr)
library(lme4)
library(dplyr)
library(scales)
library(car)
library(boot)
library(multcomp)

varclasses <- rep("numeric",16)
varclasses[c(1,8,10,12,16)] <- "character"
varclasses[c(2:6,14,15)] <- "factor"
termin <- read.delim("data_termin_main_190301.txt", colClasses = varclasses)

# Fix missing values
termin[termin$end.date=="" ,]$end.date <- NA
termin[termin$ecl.date=="" ,]$ecl.date <- NA

# Fix time variables
termin$coldmonths <- as.numeric(paste(termin$tmt))
termin$st.date <- date(termin$st.date)
termin$end.date <- date(termin$end.date)
termin$ecl.date <- date(termin$ecl.date)

# Precise time spent in cold
termin$colddays <- as.numeric( termin$end.date - termin$st.date )
termin[termin$tmt=="0",]$colddays <- 0

# Time to eclosion (calculated from start date for tmt 0; from end date for all other treatments)
termin$ecl.time <- as.numeric( termin$ecl.date - termin$end.date )
termin[termin$tmt==0,]$ecl.time <- termin[termin$tmt==0,]$ecl.date - termin[termin$tmt==0,]$st.date

# Weight loss, as proportion of starting weight
termin$weightloss <- (termin$st.weight-termin$end.weight)/termin$st.weight*100

# Dataset summary
summary(termin)



#### STATISTICAL ANALYSES ----------------------------------------------------------------

## Mortality

termin$success <- 0
termin[termin$status=="ok",]$success <- 1

succmod1 <- glm(success ~ pop*tmt, family="binomial", data=termin)
Anova(succmod1)
# Remove interaction to analyze interpopulation differences:
succmod2 <- glm(success ~ pop+tmt, family="binomial", data=termin)
Anova(succmod2)
comparisonsSuc <- glht(succmod2, mcp(pop="Tukey"))
summary(comparisonsSuc) # Pairwise differences


## Time to eclose

# Three months or less (exclude nondiapausing individuals in 0-cold treatment):
termin_diap <- subset(termin, coldmonths < 4) 
termin_diap[termin_diap$tmt=="0" & !is.na(termin_diap$ecl.time) & termin_diap$ecl.time<=20,]$ecl.time <- NA
devmodA <- lmer(ecl.time ~ sex + colddays * pop + (1|family), data=termin_diap ) 
Anova(devmodA) # Interaction nonsignificant; remove
devmodA1 <- lmer(ecl.time ~ sex + colddays + pop + (1|family), data=termin_diap )
Anova(devmodA1)
comparisonsA1<-glht(devmodA1, mcp(pop="Tukey"))
summary(comparisonsA1) # Pairwise contrasts before three months

# After three months:
devmodB <- lmer(ecl.time ~ sex + colddays * pop + (1|family), data=subset(termin, coldmonths>3) ) # Interaction nonsignificant; remove
devmodB1 <- lmer(ecl.time ~ sex + colddays + pop + (1|family), data=subset(termin, coldmonths>3) )
Anova(devmodB1)
comparisonsB1<-glht(devmodB1, mcp(pop="Tukey"))
summary(comparisonsB1) # Pairwise contrasts after three months



#### TIME TO ECLOSE WITH BOOTSTRAPPED CONFIDENCE INTERVALS ----------------------------------


# Function for bootstrapping development times for a given treatment
bootconfint <- function(input) {
  # Resample dataset and fit by-population mixed models
  boots <- lapply(1:1000, function(i) dplyr::sample_n(input, nrow(input), replace = TRUE))
  bootmods <- lapply(boots, function(d) purrr::possibly(lmer,otherwise = NA)(ecl.time ~ pop-1 + (1|family), data=d)  )
  # Extract coefficients ( = mean dev. time per population) from simulated model set
  bootset <- list(B=unlist(lapply(bootmods, function(m) if(length(fixef(m)==4)) {return(fixef(m)[1])} else {return(NA)} )),
                  L=unlist(lapply(bootmods, function(m) if(length(fixef(m)==4)) {return(fixef(m)[2])} else {return(NA)} )),
                  R=unlist(lapply(bootmods, function(m) if(length(fixef(m)==4)) {return(fixef(m)[3])} else {return(NA)} )),
                  V=unlist(lapply(bootmods, function(m) if(length(fixef(m)==4)) {return(fixef(m)[4])} else {return(NA)} )) )
  # Save bootstrapped distribution of means as confidence interval
  bootquants <- lapply(bootset, function(s) quantile(s, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) )
  return(bootquants)  }

# Function for extracting median (i=2) or confidence limits (i= 1 and 3) for a given population
quantsperpop <- function(pop, i) {
  c( bootquants.0[[pop]][i], bootquants.1[[pop]][i], bootquants.2[[pop]][i], bootquants.3[[pop]][i],
     bootquants.4[[pop]][i], bootquants.5[[pop]][i], bootquants.6[[pop]][i], bootquants.8[[pop]][i] )   }

# Estimate confidence intervals for each treatment separately:
#set.seed(413)
#bootquants.0 <- bootconfint( subset(termin, tmt==0 & status!="dead" & ecl.time>20) )
#bootquants.1 <- bootconfint( subset(termin, tmt==1 & status!="dead") )
#bootquants.2 <- bootconfint( subset(termin, tmt==2 & status!="dead") )
#bootquants.3 <- bootconfint( subset(termin, tmt==3 & status!="dead") )
#bootquants.4 <- bootconfint( subset(termin, tmt==4 & status!="dead") )
#bootquants.5 <- bootconfint( subset(termin, tmt==5 & status!="dead") )
#bootquants.6 <- bootconfint( subset(termin, tmt==6 & status!="dead") )
#bootquants.8 <- bootconfint( subset(termin, tmt==8 & status!="dead") )

# Export bootstrapped CIs to file
#devboot.lower <- data.frame(R=quantsperpop("R",1), V=quantsperpop("V",1), B=quantsperpop("B",1), L=quantsperpop("L",1))
#devboot.means <- data.frame(R=quantsperpop("R",2), V=quantsperpop("V",2), B=quantsperpop("B",2), L=quantsperpop("L",2))
#devboot.upper <- data.frame(R=quantsperpop("R",3), V=quantsperpop("V",3), B=quantsperpop("B",3), L=quantsperpop("L",3))
#write.table(devboot.lower, file="Devtime confint lower.txt", row.names=FALSE, quote=FALSE, sep="\\t")
#write.table(devboot.means, file="Devtime confint means.txt", row.names=FALSE, quote=FALSE, sep="\\t")
#write.table(devboot.upper, file="Devtime confint upper.txt", row.names=FALSE, quote=FALSE, sep="\\t")




#### GENERATING FIGURES -----------------------------------------------------

# Figure 3 prep

# Import bootstrapped CIs from file
devboot.lower <- read.delim("Devtime confint lower.txt")
devboot.means <- read.delim("Devtime confint means.txt")
devboot.upper <- read.delim("Devtime confint upper.txt")

# Import eclosion times of nondiapausing individuals and bootstrap a confidence interval
dirtimes <- read.delim("nondiapause_times_190314.txt")$ecl.time
set.seed(413)
boot.dirtimes <- apply(array(1:1000), 1, function(i) mean(sample(dirtimes, length(dirtimes), replace=TRUE)))
bootq.dirtimes <- quantile(boot.dirtimes, c(0.025,0.5,0.975))


# Exporting figure 3
pdf(file="Fig3.pdf", width=8, height=6, pointsize=20)

# Figure 3a: development success rates
succ.cols <- magma(3, begin=0.75, end=0.3)
spineplot(status~tmt, data=subset(termin, pop=="R"), ylevels=c(3,2,1), xlab="Months in cold", ylab="", main="Stockholm",
          yaxlabels=c("OK","FE","DP"), col=succ.cols, border=NA)
spineplot(status~tmt, data=subset(termin, pop=="V"), ylevels=c(3,2,1), xlab="Months in cold", ylab="", main="Highlands",
          yaxlabels=c("OK","FE","DP"), col=succ.cols, border=NA)
spineplot(status~tmt, data=subset(termin, pop=="B"), ylevels=c(3,2,1), xlab="Months in cold", ylab="", main="Gotland",
          yaxlabels=c("OK","FE","DP"), col=succ.cols, border=NA)
spineplot(status~tmt, data=subset(termin, pop=="L"), ylevels=c(3,2,1), xlab="Months in cold", ylab="", main="Skåne",
          yaxlabels=c("OK","FE","DP"), col=succ.cols, border=NA)

# Figure 3b: pupal development times
xat <- c(0:6,8)
plot(1, pch="", xlim=c(-0.3,8.3), ylim=c(12,65), xlab="Months in cold", ylab="Pupal development (days)", xaxt="n")
axis(side=1, at=xat)
rect(-5, bootq.dirtimes[1], 9, bootq.dirtimes[3], col="grey90", border=NA) # Nondiapause confint
abline(h=bootq.dirtimes[2], lty="dotted") # Nondiapause mean
arrows(0:6-0.09, devboot.lower$R[1:7], 0:6-0.09, devboot.upper$R[1:7], angle=90, code=3, length=0.05, col="blue")
arrows(xat-0.03, devboot.lower$V, xat-0.03, devboot.upper$V, angle=90, code=3, length=0.05, col="cornflowerblue")
arrows(xat+0.03, devboot.lower$B, xat+0.03, devboot.upper$B, angle=90, code=3, length=0.05, col="orange")
arrows(xat+0.09, devboot.lower$L, xat+0.09, devboot.upper$L, angle=90, code=3, length=0.05, col="red")
points(xat-0.09, devboot.means$R, type="p", pch=21, col="blue", bg="blue")
points(xat-0.03, devboot.means$V, type="p", pch=22, col="cornflowerblue", bg="cornflowerblue")
points(xat+0.03, devboot.means$B, type="p", pch=25, col="orange", bg="orange")
points(xat+0.09, devboot.means$L, type="p", pch=24, col="red", bg="red")
legend("topright", legend=c("Stockholm","Highlands","Gotland","Skåne"), col=c("blue","cornflowerblue","orange","red"),
       pch=c(21,22,25,24), pt.bg=c("blue","cornflowerblue","orange","red"), inset=0.02)

dev.off()