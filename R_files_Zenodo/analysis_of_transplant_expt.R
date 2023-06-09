# Analysis of Judith Bachmann's 2015 tranplant experiment
# This script accompanies: Bachmann & Van Buskirk, 2021. Adaptation to elevation but limited local adaptation in an amphibian. Evolution 75:in press.


# Begin by loading the functions at the bottom.
#
# Define the working directory in which all the data are kept.
#
setwd(" .... your path here .... ")
#
library(lme4)
library(pscl)
library(MCMCglmm)
library(effects)
library(lmerTest)
library(raster)
library(fields)
#
# Define colors that will be used throughout. I figured out these colors at http://html-color-codes.info/
#
low.color       <- "#FEA100"       # an orange color
low.color.pale  <- adjustcolor(low.color, alpha.f = 0.08)
high.color      <- "#1A0778"       # a dark blue
high.color.pale <- adjustcolor(high.color, alpha.f = 0.08)
away.color      <- "#C00F0F"       # a red
away.color.dark <- "#820404"
home.color      <- "#29AB5B"
home.color.dark <- "#0C5F2C"
low.away.color  <- "#FFCD00"
high.away.color <- "#378EAE"
#
#
# Import and analyze data from the lab experiment
#
lab               <- read.table("lab expt data.txt", header = TRUE, stringsAsFactors = FALSE)
lab$dev.42        <- (42-25) / lab$days42
lab$dev.45        <- (45-25) / lab$days45
lab$day.42        <- ifelse(lab$time.42 == "PM", (lab$day.42 + 0.5), lab$day.42)
lab$date42.julian <- month.day.to.Julian(lab$month.42, lab$day.42)
lab$date45.julian <- month.day.to.Julian(lab$month.45, lab$day.45)
lab$date.startexp.julian <- month.day.to.Julian(lab$month.start, lab$day.start)
lab               <- lab[ , c(2, 3, 16, 4, 26, 19, 22, 23, 12, 13)]
fam.means         <- aggregate(lab[ , c(6:10)], by = list(lab$pop, lab$family, lab$elevation, lab$round, lab$date.startexp.julian), mean, na.rm = TRUE)
source.means      <- aggregate(fam.means[ , c(6:10)], by = list(fam.means$Group.1, fam.means$Group.3, fam.means$Group.4, fam.means$Group.5), mean, na.rm = TRUE)
source.means      <- rename(source.means, c("Group.1", "Group.2", "Group.3", "Group.4"), c("pop", "elev", "round", "start.date"))
source.sd         <- aggregate(fam.means[ , c(6:10)], by = list(fam.means$Group.1, fam.means$Group.3, fam.means$Group.4, fam.means$Group.5), sd, na.rm = TRUE)
source.sd         <- rename(source.sd, names(source.sd), c("pop", "elev", "round", "start.date", "surv.SD", "dev.42.SD", "dev.45.SD", "mass.42.SD", "mass.45.SD"))
source.N          <- aggregate(fam.means[ , 6], by = list(fam.means$Group.1, fam.means$Group.3, fam.means$Group.4, fam.means$Group.5), length)
names(source.N)   <- c("pop", "elev", "round", "start.date", "N")
source.means      <- cbind(source.N, source.means, source.sd)
source.means$surv.SE   <- source.means$surv.SD / sqrt(source.means$N)
source.means$mass42.SE <- source.means$mass.42.SD / sqrt(source.means$N)
source.means$dev42.SE  <- source.means$dev.42.SD / sqrt(source.means$N)
source.means <- source.means[ , c(1:3, 5, 10, 24, 11, 26, 13, 25)]
lab <- rename(lab, "family", "clutch")
source.means
#
# Univariate tests with "source.elev" included in fixed effects. This is in Table S3.
#
lab$stand.dev  <- (lab$dev.42 - mean(lab$dev.42, na.rm = TRUE)) / sd(lab$dev.42, na.rm = TRUE)
lab$stand.mass <- (lab$mass.42 - mean(lab$mass.42, na.rm = TRUE)) / sd(lab$mass.42, na.rm = TRUE)
#
lab$source.elev <- ifelse(lab$elev > 2000, "high", "low")
m.surv.full     <- glmer(survived ~ factor(round) + source.elev + factor(round):source.elev + (1|pop/clutch), family = "binomial", data = lab)
m.surv.noclutch <- glmer(survived ~ factor(round) + source.elev + factor(round):source.elev + (1|pop), family = "binomial", data = lab)
m.surv.nopop    <- glm(survived   ~ factor(round) + source.elev + factor(round):source.elev, family = "binomial", data = lab)
summary(m.surv.full)
as.data.frame(coef(summary(m.surv.full)))
anova(m.surv.noclutch, m.surv.nopop)   # test of source population
anova(m.surv.full, m.surv.noclutch)    # test of clutch within population
#
m.dev.full     <- lmer(stand.dev ~ factor(round) + source.elev + factor(round):source.elev + (1|pop/clutch), data = lab)
m.dev.noclutch <- lmer(stand.dev ~ factor(round) + source.elev + factor(round):source.elev + (1|pop), data = lab)
m.dev.nopop    <- glm(stand.dev  ~ factor(round) + source.elev + factor(round):source.elev, data = lab)
summary(m.dev.full)
anova(m.dev.full, type = 3, ddf = c("Satterthwaite"))
anova(m.dev.noclutch, m.dev.nopop)   # test of source population
anova(m.dev.full, m.dev.noclutch)    # test of clutch within population
#
m.mass.full     <- lmer(stand.mass ~ factor(round) + source.elev + factor(round):source.elev + (1|pop/clutch), data = lab)
m.mass.noclutch <- lmer(stand.mass ~ factor(round) + source.elev + factor(round):source.elev + (1|pop), data = lab)
m.mass.nopop    <- glm(stand.mass  ~ factor(round) + source.elev + factor(round):source.elev, data = lab)
summary(m.mass.full)
anova(m.mass.full, type = 3, ddf = c("Satterthwaite"))
anova(m.mass.noclutch, m.mass.nopop)   # test of source population
anova(m.mass.full, m.mass.noclutch)    # test of clutch within population
#
#
####################
# Draw figure showing results of the lab experiment: Fig. S8 in the MS.
#
first.round.color  <- "#A00C11"   # a dark red
second.round.color <- "#03A8D6"   # a blue
source.means$x <- c(1.8, 0.8, 3.0, 4.0, 2.0, 1.0, 3.2, 4.2)    # run these two lines just once!!!
source.means   <- source.means[ order(source.means$x), ]
first.round    <- c(1,3,5,7)
second.round   <- c(2,4,6,8)
source.means$surv.lwr <- source.means$survived - source.means$surv.SE
source.means$surv.upr <- source.means$survived + source.means$surv.SE
source.means$mass.lwr <- source.means$mass.42 - source.means$mass42.SE
source.means$mass.upr <- source.means$mass.42 + source.means$mass42.SE
source.means$dev.lwr <- source.means$dev.42 - source.means$dev42.SE
source.means$dev.upr <- source.means$dev.42 + source.means$dev42.SE
ylimits.dev  <- c( min(source.means$dev.lwr), max(source.means$dev.upr))
ylimits.mass <- c( min(source.means$mass.lwr), max(source.means$mass.upr))
ylimits.surv <- c( min(source.means$surv.lwr), max(source.means$surv.upr))
xlimits <- c(0.2, 4.7)
#
pdf(" .... your path here ..../Fig S8.pdf", useDingbats = FALSE, width = 5, height = 9)
plot.new()
# Bottom panel is mass at emergence.
par(new = "TRUE", plt = c(0.25, 0.85, 0.10, 0.35))
  plot(source.means$x, source.means$ylimits.mass, pch = '', xlim = xlimits, ylim = ylimits.mass, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n" )
  rect(xleft = 0, ybottom = 0.8*ylimits.mass[1], xright = 2.5, ytop = 1.1*ylimits.mass[2], col = low.color.pale, border = NA)
  rect(xleft = 2.5, ybottom = 0.8*ylimits.mass[1], xright = 1.1*xlimits[2], ytop = 1.1*ylimits.mass[2], col = high.color.pale, border = NA)
  box(lwd = 1.4)
  abline(v = 2.5, lwd = 1.4)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = "440", side = 2, line = 0.6, las = 1, at = 440, cex=1); mtext(text = "460", side = 2, line = 0.5, las = 1, at = 460, cex=1); mtext(text = "480", side = 2, line = 0.5, las = 1, at = 480, cex=1); mtext(text = "500", side = 2, line = 0.5, las = 1, at = 500, cex=1); mtext(text = "520", side = 2, line = 0.5, las = 1, at = 520, cex=1); mtext(text = "540", side = 2, line = 0.5, las = 1, at = 540, cex=1); mtext(text = "560", side = 2, line = 0.5, las = 1, at = 560, cex=1)
  mtext(text = "Source population", side = 1, line = 1.8, las = 1, cex = 1.2)
  mtext(text = c("feer","siec","bern","flue"), side = 1, line = 0.4, at = c(0.9, 1.9, 3.1, 4.1), col = c(low.color, low.color, high.color, high.color), las = 1, cex = 1.0)
  mtext(text = "Mass at metamorphosis (mg)", side = 2, line = 2.7, cex = 1.1)
  for(i in 1:dim(source.means)[1]) {
    lines(x = c(source.means$x[i], source.means$x[i]), y = c(source.means$mass.lwr[i], source.means$mass.upr[i]), col = "black", lwd = 1)
    }
  lines(x = source.means$x[first.round], y = source.means$mass.42[first.round], col = first.round.color, lwd = 1.2)
  lines(x = source.means$x[second.round], y = source.means$mass.42[second.round], col = second.round.color, lwd = 1.2)
  points(x = source.means$x[first.round], y = source.means$mass.42[first.round], pch = 21, col = "white", bg = first.round.color, cex = 2.3, lwd = 2)
  points(x = source.means$x[second.round], y = source.means$mass.42[second.round], pch = 21, col = "white", bg = second.round.color, cex = 2.3, lwd = 2)
  text(labels = "C", x = 0.3, y = 552, cex = 1.6)
  text(labels = c("low elevation", "source pops", "high elevation", "source pops"), x = c(1.3, 1.3, 3.7, 3.7), y = c(455, 444, 545, 534), col = c(low.color, low.color, high.color, high.color), cex = 0.8)
#
# Middle panel is development rate.
par(new = "TRUE", plt = c(0.25, 0.85, 0.40, 0.65))
  plot(source.means$x, source.means$ylimits.dev, pch = '', xlim = xlimits, ylim = ylimits.dev, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n" )
  rect(xleft = 0, ybottom = 0.8*ylimits.dev[1], xright = 2.5, ytop = 1.1*ylimits.dev[2], col = low.color.pale, border = NA)
  rect(xleft = 2.5, ybottom = 0.8*ylimits.dev[1], xright = 1.1*xlimits[2], ytop = 1.1*ylimits.dev[2], col = high.color.pale, border = NA)
  box(lwd = 1.4)
  abline(v = 2.5, lwd = 1.4)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = "0.48", side = 2, line = 0.6, las = 1, at = 0.48, cex=1); mtext(text = "0.50", side = 2, line = 0.5, las = 1, at = 0.5, cex=1); mtext(text = "0.52", side = 2, line = 0.5, las = 1, at = 0.52, cex=1); mtext(text = "0.54", side = 2, line = 0.5, las = 1, at = 0.54, cex=1); mtext(text = "0.56", side = 2, line = 0.5, las = 1, at = 0.56, cex=1); mtext(text = "0.58", side = 2, line = 0.5, las = 1, at = 0.58, cex=1)
  mtext( text = c("feer","siec","bern","flue"), side = 1, line = 0.4, at = c(0.9, 1.9, 3.1, 4.1), col = c(low.color, low.color, high.color, high.color), las = 1, cex = 1.0)
  mtext(text = "Development rate (stages/day)", side = 2, line = 2.7, cex = 1.1)
  for(i in 1:dim(source.means)[1]) {
    lines(x = c(source.means$x[i], source.means$x[i]), y = c(source.means$dev.lwr[i], source.means$dev.upr[i]), col = "black", lwd = 1)
    }
  lines(x = source.means$x[first.round], y = source.means$dev.42[first.round], col = first.round.color, lwd = 1.2)
  lines(x = source.means$x[second.round], y = source.means$dev.42[second.round], col = second.round.color, lwd = 1.2)
  points(x = source.means$x[first.round], y = source.means$dev.42[first.round], pch = 21, col = "white", bg = first.round.color, cex = 2.3, lwd = 2)
  points(x = source.means$x[second.round], y = source.means$dev.42[second.round], pch = 21, col = "white", bg = second.round.color, cex = 2.3, lwd = 2)
  text(labels = "B", x = 0.3, y = 0.583, cex = 1.6)
  text(labels = c("short duration", "(March)", "long duration", "(June)"), x = c(1.5, 1.5, 3.8, 3.8), y = c(.525, .515, .504, .494), col = c(first.round.color, first.round.color, second.round.color, second.round.color), cex = 0.8)
#
# Top panel is survival.
par(new = "TRUE", plt = c(0.25, 0.85, 0.70, 0.95))
  plot(source.means$x, source.means$survived, pch = '', xlim = xlimits, ylim = ylimits.surv, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n" )
  rect(xleft = 0, ybottom = 0.8*ylimits.surv[1], xright = 2.5, ytop = 1.1*ylimits.surv[2], col = low.color.pale, border = NA)
  rect(xleft = 2.5, ybottom = 0.8*ylimits.surv[1], xright = 1.1*xlimits[2], ytop = 1.1*ylimits.surv[2], col = high.color.pale, border = NA)
  box(lwd = 1.4)
  abline(v = 2.5, lwd = 1.4)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = "0.7", side = 2, line = 0.5, las = 1, at = 0.7, cex=1); mtext(text = "0.8", side = 2, line = 0.5, las = 1, at = 0.8, cex=1); mtext(text = "0.9", side = 2, line = 0.5, las = 1, at = 0.9, cex=1); mtext(text = "1.0", side = 2, line = 0.5, las = 1, at = 1, cex=1)
  mtext( text = c("feer","siec","bern","flue"), side = 1, line = 0.4, at = c(0.9, 1.9, 3.1, 4.1), col = c(low.color, low.color, high.color, high.color), las = 1, cex = 1.0)
  mtext(text = "Survival to metamorphosis", side = 2, line = 2.7, cex = 1.1)
  for(i in 1:dim(source.means)[1]) {
    lines(x = c(source.means$x[i], source.means$x[i]), y = c(source.means$surv.lwr[i], source.means$surv.upr[i]), col = "black", lwd = 1)
    }
  lines(x = source.means$x[first.round], y = source.means$survived[first.round], col = first.round.color, lwd = 1.2)
  lines(x = source.means$x[second.round], y = source.means$survived[second.round], col = second.round.color, lwd = 1.2)
  points(x = source.means$x[first.round], y = source.means$survived[first.round], pch = 21, col = "white", bg = first.round.color, cex = 2.3, lwd = 2)
  points(x = source.means$x[second.round], y = source.means$survived[second.round], pch = 21, col = "white", bg = second.round.color, cex = 2.3, lwd = 2)
  text(labels = "A", x = 0.3, y = 0.98, cex = 1.6)
  # Label the panel at the top
  mtext(text = c("Fig. S8", "drawn by 'analysis of Judith transplant.R'"), cex = c(0.7, 0.45), at = c(0.4, 3.5), line = 1.3)
dev.off()
#
###########################################################
#
############################################################
#
#
#
# Analysis of the main transplant experiment.
#
d <- read.table(file = "field expt data.txt", header = TRUE, stringsAsFactors = FALSE)
d$date.startexp <- d$date.cage.removed <- d$date.stage42 <- d$date.stage45 <- NULL
# Correct some errors in the dataset
#    home.elev mislabelled for some cages in "ellw".
d$home.elev <- ifelse(d$test.pond == "ellw" & d$source.pop == "flue", "away_elev", d$home.elev)
#    home.pond is mislabelled for siec cage 32.
d$home.pond <- ifelse(d$test.pond == "siec" & d$cage == 32 & d$source.pop == "feer", "away_pond", d$home.pond)
#    home.elev and home.pond are mislabelled for siec cage 2.
d$home.elev <- ifelse(d$test.pond == "siec" & d$cage == 2 & d$source.pop == "flue", "away_elev", d$home.elev)
d$home.pond <- ifelse(d$test.pond == "siec" & d$cage == 2 & d$source.pop == "flue", "away_pond", d$home.pond)
#    home.elev and home.pond are mislabelled for feer cage 20.
d$home.elev <- ifelse(d$test.pond == "feer" & d$cage == 20, "home_elev", d$home.elev)
d$home.pond <- ifelse(d$test.pond == "feer" & d$cage == 20, "home_pond", d$home.pond)
#    One tadpole in cage 30 has the wrong home.elev.
d$home.elev <- ifelse(d$test.pond == "feer" & d$cage == 30, "home_elev", d$home.elev)
#
# Errors in the date of stage 42 or 45:
d$days42 <- ifelse(d$ID == "bide1301", 38, d$days42)        # age at stage 42 should be 38, not 35
d$days42 <- ifelse(d$ID == "bern1509", 46, d$days42)        # age at stage 42 should be 46, not 45
d$days42 <- ifelse(d$ID == "ellw0501", 52, d$days42)        # age at stage 42 should be 52, not 54
d$days45 <- ifelse(d$ID == "muot2905", 44, d$days45)        # age at stage 45 should be 44, not 41
d$days45 <- ifelse(d$ID == "flue3108", 65, d$days45)        # age at stage 45 should be 65, not 95
d$days45 <- ifelse(d$ID == "siec2701", 73, d$days45)        # age at stage 45 should be 73, not 162
d$days45 <- ifelse(d$ID == "siec0605", 100, d$days45)       # age at stage 45 should be 100, not 162
d$days45 <- ifelse(d$ID == "siec0606", 102, d$days45)       # age at stage 45 should be 102, not 164
d$days45 <- ifelse(d$ID == "siec0607", 112, d$days45)       # age at stage 45 should be 112, not 174
d$days45 <- ifelse(d$ID == "siec0608", 102, d$days45)       # age at stage 45 should be 102, not 162
d$days45 <- ifelse(d$ID == "siec0505", 105, d$days45)       # age at stage 45 should be 105, not 167
d$days45 <- ifelse(d$ID == "siec0601", 97.5, d$days45)      # age at stage 45 should be 97.5
d$days45 <- ifelse(d$ID == "siec0702", 106, d$days45)       # age at stage 45 should be 106, not 168
d$days45 <- ifelse(d$ID == "siec0709", 100, d$days45)       # age at stage 45 should be 100, not 162
d$days45 <- ifelse(d$ID == "siec0711", 100, d$days45)       # age at stage 45 should be 100, not 162
d$days45 <- ifelse(d$ID == "siec0805", 98, d$days45)        # age at stage 45 should be 98, not 160
d$days45 <- ifelse(d$ID == "siec0908", 102, d$days45)       # age at stage 45 should be 102, not 164
d$days45 <- ifelse(d$ID == "siec0904", 100, d$days45)       # age at stage 45 should be 100, not 162
d$days45 <- ifelse(d$ID == "siec0905", 101.5, d$days45)     # age at stage 45 should be 101.5, not 163
d$days45 <- ifelse(d$ID == "siec0906", 100, d$days45)       # age at stage 45 should be 100, not 162
d$days45 <- ifelse(d$ID == "siec1704", 103, d$days45)       # age at stage 45 should be 103, not 165
d$days45 <- ifelse(d$ID == "siec2206", 100, d$days45)       # age at stage 45 should be 100, not 164
d$days45 <- ifelse(d$ID == "siec2211", 100.5, d$days45)     # age at stage 45 should be 100.5, not 162
d$days45 <- ifelse(d$ID == "siec2710", 106, d$days45)       # age at stage 45 should be 106, not 168
d$days45 <- ifelse(d$ID == "siec3008", 102, d$days45)       # age at stage 45 should be 102, not 164
d$days45 <- ifelse(d$ID == "siec3009", 100, d$days45)       # age at stage 45 should be 100, not 162
d$days45 <- ifelse(d$ID == "siec0401", 98, d$days45)        # age at stage 45 should be 98, not 160
d$days45 <- ifelse(d$ID == "siec0102", 101.5, d$days45)     # age at stage 45 should be 101.5, not 163
d$days45 <- ifelse(d$ID == "ellw1209", 100.5, d$days45)     # age at stage 45 should be 100.5, not 162
d$days45 <- ifelse(d$ID == "ellw0101", 60, d$days45)        # age at stage 45 should be 60, not 57
d$days45 <- ifelse(d$ID == "ellw0102", 60, d$days45)        # age at stage 45 should be 60, not 57
d$days45 <- ifelse(d$ID == "ellw0103", 60, d$days45)        # age at stage 45 should be 60, not 57
d$days45 <- ifelse(d$ID == "ellw0104", 60, d$days45)        # age at stage 45 should be 60, not 57
d$days45 <- ifelse(d$ID == "ellw0601", 60, d$days45)        # age at stage 45 should be 60, not 57
d$days45 <- ifelse(d$ID == "ellw0602", 60, d$days45)        # age at stage 45 should be 60, not 57
d$days45 <- ifelse(d$ID == "ellw0603", 60, d$days45)        # age at stage 45 should be 60, not 57
d$days45 <- ifelse(d$ID == "ellw0604", 60, d$days45)        # age at stage 45 should be 60, not 57
d$days45 <- ifelse(d$ID == "ellw0605", 60, d$days45)        # age at stage 45 should be 60, not 57
d$days42 <- ifelse(d$ID == "ellw0701", 52, d$days42)        # age at stage 42 should be 52, not 54 (mass42 probably also too small)
d$days42 <- ifelse(d$ID == "ellw0702", 52, d$days42)        # age at stage 42 should be 52, not 54 (mass42 probably also too small)
d$days42 <- ifelse(d$ID == "ellw0901", 54, d$days42)        # age at stage 42 should be 54, not 57 (mass42 probably also too small)
d$days42 <- ifelse(d$ID == "ellw1601", 52, d$days42)        # age at stage 42 should be 52, not 54 (mass42 probably also too small)
d$days42 <- ifelse(d$ID == "siec0301", 61, d$days42)        # age at stage 42 should be 61, not 64
d$days42 <- ifelse(d$ID == "siec1101", 67, d$days42)        # age at stage 42 should be 67, not 70
d$days45 <- ifelse(d$ID == "siec1301", 60, d$days45)        # age at stage 45 should be 60, not 57
d$days42 <- ifelse(d$ID == "siec2006", 61, d$days42)        # age at stage 42 should be 61, not 64
d$days42 <- ifelse(d$ID == "siec2101", 54, d$days42)        # age at stage 42 should be 54, not 57
d$days42 <- ifelse(d$ID == "siec2102", 54, d$days42)        # age at stage 42 should be 54, not 57
d$days42 <- ifelse(d$ID == "siec2303", 61, d$days42)        # age at stage 42 should be 61, not 64
d$days42 <- ifelse(d$ID == "siec2501", 61, d$days42)        # age at stage 45 should be 61, not 64
d$days42 <- ifelse(d$ID == "siec3001", 67, d$days42)        # age at stage 45 should be 67, not 70
# Next are errors in mass.
# Many of these were caught and measured at or near stage 45, and the stage 45 mass was
# substituted for stage 42 mass. I estimated stage 42 mass or stage 45 mass for all those with a
# ratio (42/45) of < 1.25.
d$mass.45<- ifelse(d$ID == "feer2401", 480, d$mass.45)      # lost much too much weight
d$mass.42<- ifelse(d$ID == "ellw2005", 387, d$mass.42)      # gained weight
d$mass.45<- ifelse(d$ID == "feer3106", 354, d$mass.45)      # gained weight
d$mass.45<- ifelse(d$ID == "munt0408", 443, d$mass.45)      # did not lose enough weight
d$mass.42<- ifelse(d$ID == "siec0403", 526, d$mass.42)      # weight at 42 too high
d$mass.42<- ifelse(d$ID == "ellw0501", 130, d$mass.42)      # weight at 42 too low
d$mass.42<- ifelse(d$ID == "siec2101", 500, d$mass.42)      # same weight
d$mass.45<- ifelse(d$ID == "ellw1801", 65, d$mass.45)       # lost too much weight
d$mass.42<- ifelse(d$ID == "ellw1904", 250, d$mass.42)      # no loss of weight
d$mass.42<- ifelse(d$ID == "ellw1903", 281, d$mass.42)      # little loss of weight
d$mass.45<- ifelse(d$ID == "ellw2003", 123, d$mass.45)      # little loss of weight
d$mass.42<- ifelse(d$ID == "ellw1601", 190, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "ellw0901", 352, d$mass.42)      # no loss of weight
d$mass.45<- ifelse(d$ID == "ellw1701", 190, d$mass.45)      # little loss of weight
d$mass.45<- ifelse(d$ID == "ellw1702", 200, d$mass.45)      # little loss of weight
d$mass.45<- ifelse(d$ID == "ellw1902", 150, d$mass.45)      # little loss of weight
d$mass.45<- ifelse(d$ID == "ellw1401", 150, d$mass.45)      # little loss of weight
d$mass.42<- ifelse(d$ID == "ellw0701", 200, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "ellw0702", 160, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "muot1109", 180, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "muot1410", 270, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "siec2502", 400, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "siec2501", 400, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "siec2601", 300, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "muot3001", 100, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "feer2702", 549, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "bide3203", 470, d$mass.42)      # little loss of weight
d$mass.45<- ifelse(d$ID == "feer1607", 320, d$mass.45)      # little loss of weight
d$mass.42<- ifelse(d$ID == "bern1509", 340, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "ellw0101", 300, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "ellw0102", 300, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "ellw0103", 330, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "ellw0104", 290, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "ellw0602", 280, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "ellw1201", 220, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "ellw1202", 230, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "siec0301", 320, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "siec1101", 490, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "siec2006", 330, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "siec2102", 500, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "siec3001", 410, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "siec2701", 330, d$mass.42)      # same weight
d$mass.42<- ifelse(d$ID == "siec2303", 340, d$mass.42)      # same weight
d$mass.45<- ifelse(d$ID == "siec2602", 170, d$mass.45)      # little loss of weight
d$mass.42<- ifelse(d$ID == "siec2001", 330, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "siec2002", 280, d$mass.42)      # little loss of weight
d$mass.45<- ifelse(d$ID == "feer2701", 427, d$mass.45)      # little loss of weight
d$mass.45<- ifelse(d$ID == "flue1806", 407, d$mass.45)      # little loss of weight
d$mass.45<- ifelse(d$ID == "flue3009", 370, d$mass.45)      # little loss of weight
d$mass.45<- ifelse(d$ID == "flue1901", 350, d$mass.45)      # little loss of weight
d$mass.45<- ifelse(d$ID == "ellw1602", 110, d$mass.45)      # little loss of weight
d$mass.42<- ifelse(d$ID == "siec1801", 250, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "bern2906", 210, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "ellw0401", 210, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "ellw2002", 290, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "feer0403", 756, d$mass.42)      # little loss of weight
d$mass.45<- ifelse(d$ID == "flue1503", 370, d$mass.45)      # little loss of weight
d$mass.45<- ifelse(d$ID == "flue1704", 400, d$mass.45)      # little loss of weight
d$mass.45<- ifelse(d$ID == "flue1706", 390, d$mass.45)      # little loss of weight
d$mass.45<- ifelse(d$ID == "munt3010", 270, d$mass.45)      # little loss of weight
d$mass.42<- ifelse(d$ID == "muot0601", 220, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "muot1907", 230, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "muot3210", 230, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "siec2301", 420, d$mass.42)      # little loss of weight
d$mass.42<- ifelse(d$ID == "muot1001", 250, d$mass.42)      # little loss of weight
#
# Other errors that I corrected in the original Excel spreadsheet:
# The last 5 flue tadpoles in ellw cage 7 were listed as source.elev = "low"
# All flue tadpoles in ellw cage 11 were listed as source.elev = "low"
# All flue tadpoles in siec cage 2 were listed as source.elev = "low"
# All feer tadpoles in feer cage 20 were listed as source.elev = "high"
#
d$ratio <- d$mass.42 / d$mass.45
d$DevRate42 <- (42-25) / d$days42
d$DevRate45 <- (45-25) / d$days45
#
# "Fitness" according to Altwegg & Reyer (2003), using age and mass at stage 42
d$normMetaAge  <- (d$days42 - mean(d$days42, na.rm=TRUE)) / sd(d$days42, na.rm=TRUE)
d$normMetaMass <- (d$mass.42 - mean(d$mass.42, na.rm=TRUE)) / sd(d$mass.42, na.rm=TRUE)
logitP         <- -0.693147 + 0.87*d$normMetaMass - 0.021*d$normMetaAge   # using Res Altwegg's Evolution result, mean fitness is 0.3520
d$fitness      <- exp(logitP) / (1 + exp(logitP))
d$fitness      <- ifelse(is.na(d$fitness), 0, d$fitness)                  # if tadpole died before stage 42, fitness = 0; now the overall mean fitness if 0.2269
d$rel.fitness  <- d$fitness / mean(d$fitness, na.rm = TRUE)
d$survival     <- ifelse(d$survival == 1 & d$fitness == 0, 0, d$survival) # those that died in the lab before stage 42 get survival = 0
d$round.fitness <- round(10*d$rel.fitness)
length(d$round.fitness); mean(d$round.fitness); sd(d$round.fitness)
cor(d$round.fitness, d$fitness)
#
#
# Draw Fig. S3 -- distribution of fitness.
#
pdf(" .... your path here ..../Fig S3.pdf", width = 6.5, height = 6)
  plot.new()
  par(new = "TRUE", xpd = FALSE, plt = c(0.2, 0.85, 0.2, 0.85), cex.axis = 1, tck = -0.02, las = 1 )
  hist(d$round.fitness, breaks=c(0,3,6,9,12,15,18,21,24,27,30,33,36,39,42), col="grey40", lwd = 1.5, las = 1, cex.axis = 1.2, ylab = NULL, xlab = NULL, cex.lab = 1.8, mar = c(15, 10, 4, 2), main = "" )
  text(expression(paste(italic("N"), " = 2760 tadpoles")), x = 13, y = 700, cex = 1)
  mtext(text = "Rounded relative fitness", side = 1, line = 2.3, cex = 1.5)
  mtext(text = "Number of tadpoles", side = 2, line = 3.1, cex = 1.5, las = 3)
  mtext(text = c("Fig. S3", "drawn by 'analysis of Judith transplant.R'"), side = 3, at = c(1,20), line = 3.1, cex = c(0.7, 0.5))
dev.off()
#
#
# Create standardized response variables: fit.stand is standardized overall, fit.s.bytestpond is standardized by test pond.
#
d$fit.stand         <- (d$fitness - mean(d$fitness, na.rm = TRUE)) / sd(d$fitness, na.rm = TRUE)
d$surv.stand        <- (d$survival - mean(d$survival, na.rm = TRUE)) / sd(d$survival, na.rm = TRUE)
d$dev.stand         <- (d$DevRate42 - mean(d$DevRate42, na.rm = TRUE)) / sd(d$DevRate42, na.rm = TRUE)
d$mass.stand        <- (d$mass.42 - mean(d$mass.42, na.rm = TRUE)) / sd(d$mass.42, na.rm = TRUE)
d$fit.s.bytestpond  <- d$fitness
d$surv.s.bytestpond <- d$survival
d$dev.s.bytestpond  <- d$DevRate42
d$mass.s.bytestpond <- d$mass.42
d <- rescale.within.groups(d, group.list = c("test.pond"), response.list = c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond"), scale = TRUE)
#
#
#########################################################
#
# Cage means.
#
cage.means <- aggregate(d[,c("survival", "mass.42", "DevRate42")], by = list(d$test.pond, d$source.pop, d$test.elev, d$source.elev, d$cage), mean, na.rm = TRUE)
names(cage.means)[1:5] <- c("test.pond", "source.pop", "test.elev", "source.elev", "cage")
head(cage.means)
#
#
#########################################################
#
# Summarize time spent in the cages and in the lab.
# Then draw a figure illustrating time in the field and the lab: Figure S2 in the manuscript.
#
quantile(d$days_pond, c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(d$days_lab, c(0.05, 0.5, 0.95), na.rm = TRUE)
mean(d$days_pond, na.rm = TRUE)
mean(d$days_lab, na.rm = TRUE)
sd(d$days_lab, na.rm = TRUE)
prop.pond <- d$days_pond / d$days42
mean(prop.pond, na.rm = TRUE)
quantile(prop.pond, c(0.05, 0.5, 0.95), na.rm = TRUE)
#
time.course   <- d[,c("test.pond", "elev", "test.elev", "ID", "days_pond", "days_lab")]
time.course   <- time.course[ complete.cases(time.course), ]
time.course   <- time.course[ order(time.course$elev, time.course$days_pond, time.course$days_lab), ]
rownames(time.course) <- NULL
dim(time.course)
time.course$x   <- 1:dim(time.course)[1]
pond.list       <- unique(time.course$test.pond)
pond.limits     <- data.frame(pond = pond.list, stringsAsFactors = FALSE)
temp            <- merge(aggregate(time.course$x, by = list(time.course$test.pond), min), aggregate(time.course$x, by = list(time.course$test.pond), max), by = "Group.1")
n.size          <- aggregate(time.course$x, by = list(time.course$test.pond), length)
names(temp)     <- c("pond", "low", "high")
names(n.size)   <- c("pond", "N")
pond.limits     <- merge(merge(pond.limits, temp, by = "pond"), n.size, by = "pond")
pond.limits     <- pond.limits[ order(pond.limits$low), ]
line.width      <- 0.3
tick.label.size <- 1
y.limits <- c(0, (1 + dim(time.course)[1]))
x.limits <- c(0, 100)
pdf(" .... your path here .... /Fig S2.pdf", width = 8, height = 6) 
plot.new()
  par(new = "TRUE", plt = c(0.16, 0.89, 0.15, 0.87) )
  plot(x = c(1,2), y = c(0,0), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = x.limits, ylim = y.limits, xaxs = "i", yaxs = "i")
  for (i in 1:dim(time.course)[1]) {
    if (time.course$test.elev[i] == "high") {
      lines(x = c(0, time.course[i, "days_pond"]), y = c(i, i), cex = line.width, col = high.color)
      lines(x = c(time.course[i, "days_pond"], (time.course[i, "days_pond"] + time.course[i, "days_lab"])), y = c(i, i), cex = line.width, col = high.color.pale)
      } else {
      lines(x = c(0, time.course[i, "days_pond"]), y = c(i, i), cex = line.width, col = low.color)
      lines(x = c(time.course[i, "days_pond"], (time.course[i, "days_pond"] + time.course[i, "days_lab"])), y = c(i, i), cex = line.width, col = low.color.pale)
      }
    }
  axis(side = 1, labels = FALSE, tck = -0.015)
  mtext(text = c("0", "20", "40", "60", "80", "100"), side = 1, line = 0.4, at = c(0,20,40,60,80,100), cex = tick.label.size)
  mtext( text = "Days after the experiment began", side = 1, line = 1.6, cex = 1.2)
  # label the ponds on the left side, and sample size on the right.
  for (i in 1:length(pond.list)) {
    mtext(text = pond.limits$pond[i], side = 2, line = 0.4, las = 1, at = (pond.limits$low[i] + pond.limits$high[i])/2, cex = tick.label.size)
    text(paste("N = ", pond.limits$N[i], sep = ""), x = 88, y = (pond.limits$low[i] + pond.limits$high[i])/2, pos = 4, cex = tick.label.size - 0.2)
    }
  mtext( text = "Test site", side = 2, las = 0, line = 2.7, cex = 1.2)
  # horizontal lines between the ponds.
  for (i in 1:7) {
    abline(h = mean(pond.limits$low[i+1], pond.limits$high[i]), lwd = 1)
    }
  box(lwd = 1.3)
  mtext( text = c("Fig. S2", "drawn by 'analysis of Judith transplant.R'"), side = 3, line = 2.6, cex = c(0.9, 0.6), at = c(0,70))
dev.off()
#
#
#########################################################
#
# MCMCglmm test for unique and parallel local adaptation.
#
# There are several steps here:
# 1) First, pull aside the data. These come from the four test ponds that also contributed source tadpoles.
#    In each test pond, we want to compare home-pond versus home-elevation versus away-elevation.
# 2) Create a new categorical variable called "treatment" with three levels: home, home.elev, away.elev.
# 3) MCMCglmm with the following model: fitness ~ treatment + test.elev/test.pond + treatment:test.pond
# 4) Next step is to perform tests of the extent of both kinds of local adaptation. This will be done by
#    calculating predicted values for the three treatments from each MCMC iteration. Then differences
#    among treatments are used to evaluate local adaptation.
#
d.sub               <- subset(d, d$test.pond %in% c("feer", "siec", "bern", "flue"))
d.sub$treatment     <- "A_home"
d.sub$treatment     <- ifelse(d.sub$home.pond == "away_pond", "B_away.pond", d.sub$treatment)
d.sub$treatment     <- ifelse(d.sub$home.elev == "away_elev", "C_away.elev", d.sub$treatment)
# recalculate round.fitness for this subset of the dataset.
d.sub$round.fitness <- round(10 * d.sub$fitness/mean(d.sub$fitness, na.rm = TRUE))
# standardize data by test pond. Not necessary because analysis is on the Poisson fitness measure.
d.sub               <- rescale.within.groups(d.sub, group.list = c("test.pond"), response.list = c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond"), scale = TRUE)
priors              <- list(R=list(V=diag(2), n=2, fix=2), G=list(G1=list(V=diag(c(1, 0.000001)), n=2), G2=list(V=diag(c(1, 0.000001)), n=2, fix=1))  )
#
# The MCMCglmm models take 8-9 minutes.
#
set.seed(1234)
fit.full <- MCMCglmm(round.fitness ~ trait + trait:test.elev + trait:treatment + trait:test.elev:treatment,
            random = ~idh(trait):block + idh(trait):cage:block, data = d.sub, family = "zipoisson", prior = priors,
            rcov = ~idh(trait):units, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
summary(fit.full)
save(list="fit.full", file = "unique.vs.parallel.fit.full.txt", ascii = TRUE)
#
# Next, devise specific tests for the two types of local adaptation. This requires predicting from the MCMCglmm model.
#
# Load the model.
load("unique.vs.parallel.fit.full.txt")
#
# Now fit the same model with zeroinfl() from package pscl, including no random effects.
#
m.zi <- zeroinfl(round.fitness ~ test.elev + treatment + test.elev:treatment, data = d.sub)
summary(m.zi)
#
# Next, predict from the ZIP model. This will be done in three ways:
# 1) Using predict.MCMCglmm, from the package.
# 2) Using the pscl object. These predicted values are very different from the predict.MCMCglmm values.
# 3) "By hand," using a function I wrote. These results are similar (but not the same as) the
#     predict.MCMCglmm results. In other words, the proportion zero is greatly underestimated, so
#     that estimated fitness is much too high.
#     When the "by-hand" function is given the binomial-intercept from the pscl analysis, the
#     results are believable and similar to those of pscl.
#
# Predict using predict.MCMCglmm():
dat1   <- as.data.frame(matrix(unlist(rep(expand.grid(block = unique(d.sub$block), cage = 1:8, stringsAsFactors = FALSE), each = 6)), ncol = 2))
dim1   <- dim(dat1)[1]
newdat <- data.frame(test.elev = rep(c("low", "high"), each = 3), treatment = rep(c("A_home", "B_away.pond", "C_away.elev"), 2), round.fitness = 10, block = "bern5", cage = 1)
pred.from.mcmcglmm <- predict(fit.full, newdata = newdat, interval = "confidence")
pred.from.mcmcglmm <- data.frame(treat = c("low.home","low.away","low.away.elev","high.home","high.away","high.away.elev"), pred = pred.from.mcmcglmm)
#
# Predict using predict.pscl():
pred.from.pscl <- predict(m.zi, newdata = newdat, interval = "confidence")
pred.from.pscl <- data.frame(treat = c("low.home","low.away","low.away.elev","high.home","high.away","high.away.elev"), pred = pred.from.pscl)
#
# Predict using my custom-built function:
fixed.pois <- fit.full$Sol[,c(1,3,5,7,9,11)]
fixed.zi   <- fit.full$Sol[,c(2,4,6,8,10,12)]
     # v.treatment is "A_home", "B_away.pond", "C_away.elev".
     # v.test.elev is "low", "high".
     # v.test.pond is "feer", "siec", "bern", or "flue".
     # Run through the six treatment combinations. The zero-inflated intercept is from pscl, not MCMCglmm.
pred1 <- zi.predict(fixed.pois, fixed.zi, v.treatment="A_home", v.test.elev="low", binom.intercept = coef(m.zi)[7])
pred2 <- zi.predict(fixed.pois, fixed.zi, v.treatment="B_away.pond", v.test.elev="low", binom.intercept = coef(m.zi)[7])
pred3 <- zi.predict(fixed.pois, fixed.zi, v.treatment="C_away.elev", v.test.elev="low", binom.intercept = coef(m.zi)[7])
pred4 <- zi.predict(fixed.pois, fixed.zi, v.treatment="A_home", v.test.elev="high", binom.intercept = coef(m.zi)[7])
pred5 <- zi.predict(fixed.pois, fixed.zi, v.treatment="B_away.pond", v.test.elev="high", binom.intercept = coef(m.zi)[7])
pred6 <- zi.predict(fixed.pois, fixed.zi, v.treatment="C_away.elev", v.test.elev="high", binom.intercept = coef(m.zi)[7])
pred.from.custom <- rbind(pred1[[1]], pred2[[1]], pred3[[1]], pred4[[1]], pred5[[1]], pred6[[1]])
pred.from.mcmcglmm
pred.from.custom
pred.from.pscl
#
# From here on, I will use the MCMC results from pred.from.custom.
# I do this because the values of fitness predicted by the MCMCglmm results are too
# high, and they're high because the predicted number of zeros is much too low.
# The results show that parameter estimates in pscl are all quite similar to MCMCglmm,
# except for the intercept of the zero-inflated component. This explains why the
# predicted incidence of zeros is so low in the MCMCglmm model.
#
#
#### Testing unique local adaptation.
# What follows are the fraction of MCMC iterations in which home-pond tadpoles
# have higher fitness than away-pond tadpoles within the same elevation.
#
t1.vs.t2.low <- pred1[[2]] - pred2[[2]]
length(t1.vs.t2.low[t1.vs.t2.low>0]) / length(t1.vs.t2.low)      # LOW elevation
t1.vs.t2.high <- pred4[[2]] - pred5[[2]]
length(t1.vs.t2.high[t1.vs.t2.high>0]) / length(t1.vs.t2.high)   # HIGH elevation
#
#
#### Test #1 for total local adaptation
# What follows are the fraction of MCMC iterations in which home-pond tadpoles
# have higher fitness than away-elevation.
#
t1.vs.t3.low  <- pred1[[2]] - pred3[[2]]
length(t1.vs.t3.low[t1.vs.t3.low>0]) / length(t1.vs.t3.low)      # LOW elevation
t1.vs.t3.high  <- pred4[[2]] - pred6[[2]]
length(t1.vs.t3.high[t1.vs.t3.high>0]) / length(t1.vs.t3.high)   # HIGH elevation
#
#
#### Test #2 for parallel local adaptation
# What follows are the fraction of MCMC iterations in which home-elevation tadpoles
# have higher fitness than away-elevation tadpoles.
#
t1t2.vs.t3.low  <- (pred1[[2]] + pred2[[2]])/2 - pred3[[2]]
length(t1t2.vs.t3.low[t1t2.vs.t3.low>0]) / length(t1t2.vs.t3.low)      # LOW elevation
t1t2.vs.t3.high  <- (pred4[[2]]+pred5[[2]])/2 - pred6[[2]]
length(t1t2.vs.t3.high[t1t2.vs.t3.high>0]) / length(t1t2.vs.t3.high)   # HIGH elevation
#
#
# What follows are the proportional decline in estimated fitness (starting from the
# fitness of the population estimated at home) due to unique local adaptation and
# parallel adaptation to elevation.
#
prop.unique.low   <- (pred1[[2]] - pred2[[2]]) / pred1[[2]]                             # LOW elevation
prop.unique.low   <- remove.outliers(prop.unique.low, units = 2, remove.na = "na")
quantile(prop.unique.low, c(0.9, 0.5, 0.1), na.rm = TRUE)
prop.parallel.low <- ((pred1[[2]] + pred2[[2]])/2 - pred3[[2]]) / pred1[[2]]
prop.parallel.low <- remove.outliers(prop.parallel.low, units = 2, remove.na = "na")
quantile(prop.parallel.low, c(0.9, 0.5, 0.1), na.rm = TRUE)
prop.unique.high  <- (pred4[[2]] - pred5[[2]]) / pred4[[2]]                             # HIGH elevation
prop.unique.high  <- remove.outliers(prop.unique.high, units = 2, remove.na = "na")
quantile(prop.unique.high, c(0.9, 0.5, 0.1), na.rm = TRUE)
prop.parallel.high <- ((pred4[[2]] + pred5[[2]])/2 - pred6[[2]]) / pred4[[2]]
prop.parallel.high <- remove.outliers(prop.parallel.high, units = 2, remove.na = "na")
quantile(prop.parallel.high, c(0.9, 0.5, 0.1), na.rm = TRUE)#
#
# Absolute magnitude of the difference (in units of expected fitness) between
# parallel adaptation and unique local adaptation, by elevation
#
quantile( (t1t2.vs.t3.low  - t1.vs.t2.low), c(0.9, 0.5, 0.1), na.rm = TRUE)   # low
quantile( (t1t2.vs.t3.high - t1.vs.t2.high), c(0.9, 0.5, 0.1), na.rm = TRUE)  # high
#
#
# Draw Figure S6 -- showing the posterior density distributions of the three
# kinds of local adaptation (total, unique, and parallel), each at high and low elevation.
# These data come from predicting specific treatment values for each MCMCglmm iteration.
#
col.tot.text    <- "red" # "#FF00D5"         # reddish
col.tot         <- adjustcolor(col.tot.text, alpha.f = 0.1)
col.unique.text <- "#038E42"         # greenish
col.unique      <- adjustcolor(col.unique.text, alpha.f = 0.1)
col.paral.text  <- "#1E00FF"         # blue
col.paral       <- adjustcolor(col.paral.text, alpha.f = 0.1)
x.limits        <- c(-12, 19.5)
line.width      <- 1.3
axis.label.size <- 1.15
#
pdf(" .... your path here .... /Fig S6.pdf", width = 6, height = 8)
plot.new()
par(new = "TRUE", plt = c(0.2, 0.8, 0.58, 0.9))
  y.limits     <- c(0,0.26)
  plot(x = c(1,2), y = c(0,0), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = x.limits, ylim = y.limits)
  abline(v = 0)
  h1 <- hist(t1.vs.t3.high, breaks = 20, plot = FALSE)
  h3 <- hist(t1t2.vs.t3.high, breaks = 20, plot = FALSE)
  h2 <- hist(t1.vs.t2.high, breaks = 20, plot = FALSE)
  polygon(x = c((min(h1$mids)-1), h1$mids, (max(h1$mids)+1)), y = c(0, h1$density, 0), col = col.tot, border = NA)
  polygon(x = c((min(h2$mids)-1), h2$mids, (max(h2$mids)+1)), y = c(0, h2$density, 0), col = col.unique, border = NA)
  polygon(x = c((min(h3$mids)-1), h3$mids, (max(h3$mids)+1)), y = c(0, h3$density, 0), col = col.paral, border = NA)
  lines(x = c((min(h1$mids)-1), h1$mids, (max(h1$mids)+1)), y = c(0, h1$density, 0), col = col.tot.text, lwd = line.width)
  lines(x = c((min(h2$mids)-1), h2$mids, (max(h2$mids)+1)), y = c(0, h2$density, 0), col = col.unique.text, lwd = line.width)
  lines(x = c((min(h3$mids)-1), h3$mids, (max(h3$mids)+1)), y = c(0, h3$density, 0), col = col.paral.text, lwd = line.width)
    # Draw three colored arrows
  arrow.lines(vect1 = t1.vs.t3.high, which.arrow = "total", y.lims = y.limits, arr.col = col.tot.text, arr.len = 0.1, text.size = 0.8, digits = 3)
  arrow.lines(vect1 = t1.vs.t2.high, which.arrow = "unique", y.lims = y.limits, arr.col = col.unique.text, arr.len = 0.07, text.size = 0.8, digits = 3)
  arrow.lines(vect1 = t1t2.vs.t3.high, which.arrow = "parallel", y.lims = y.limits, arr.col = col.paral.text, arr.len = 0.1, text.size = 0.8, digits = 3)
  abline(h = 0)
  box(lwd = 1.1)
  axis(side = 1, labels = FALSE, tck = -0.015)
  mtext(text = c("-10","0","10","20"), side = 1, line = 0.4, las = 1, at = c(-10,0,10,20), cex=1)
  axis(side = 2, labels = FALSE, tck = -0.015)
  mtext(text = c("0.0","0.1","0.2"), side = 2, line = 0.4, las = 1, at = c(0,0.1,0.2), cex=1)
  mtext( text = "Posterior density", side = 2, line = 2.1, cex = axis.label.size)
  text(labels = c("High elevation"), x = x.limits[1]-0.035*(x.limits[2]-x.limits[1]), y = y.limits[2]-0.05*(y.limits[2]-y.limits[1]), cex = 1.1, pos = 4)
  mtext(text = "Fig. S6", side = 3, at = -8, line = 2.8, cex = 0.7)
  mtext(text = 'drawn by "analysis of Judith lab expt.R"', side = 3, at = 10.1, line = 2.8, cex = 0.5)
par(new = "TRUE", plt = c(0.2, 0.8, 0.19, 0.51))
  y.limits     <- c(0,0.205)
  plot(x = c(1,2), y = c(0,0), xlab = "", ylab = "", xaxt="n", yaxt="n", pch = NA, xlim = x.limits, ylim = y.limits)
  abline(v = 0)
  h1 <- hist(t1.vs.t3.low, breaks = 20, plot = FALSE)
  h3 <- hist(t1t2.vs.t3.low, breaks = 20, plot = FALSE)
  h2 <- hist(t1.vs.t2.low, breaks = 20, plot = FALSE)
  polygon(x = c((min(h1$mids)-1), h1$mids, (max(h1$mids)+1)), y = c(0, h1$density, 0), col = col.tot, border = NA)
  polygon(x = c((min(h2$mids)-1), h2$mids, (max(h2$mids)+1)), y = c(0, h2$density, 0), col = col.unique, border = NA)
  polygon(x = c((min(h3$mids)-1), h3$mids, (max(h3$mids)+1)), y = c(0, h3$density, 0), col = col.paral, border = NA)
  lines(x = c((min(h1$mids)-1), h1$mids, (max(h1$mids)+1)), y = c(0, h1$density, 0), col = col.tot.text, lwd = line.width)
  lines(x = c((min(h2$mids)-1), h2$mids, (max(h2$mids)+1)), y = c(0, h2$density, 0), col = col.unique.text, lwd = line.width)
  lines(x = c((min(h3$mids)-1), h3$mids, (max(h3$mids)+1)), y = c(0, h3$density, 0), col = col.paral.text, lwd = line.width)
    # Draw three colored arrows
  arrow.lines(vect1 = t1.vs.t3.low, which.arrow = "total", y.lims = y.limits, arr.col = col.tot.text, arr.len = 0.1, text.size = 0.8, digits = 3)
  arrow.lines(vect1 = t1.vs.t2.low, which.arrow = "unique", y.lims = y.limits, arr.col = col.unique.text, arr.len = 0.1, text.size = 0.8, digits = 3)
  arrow.lines(vect1 = t1t2.vs.t3.low, which.arrow = "parallel", y.lims = y.limits, arr.col = col.paral.text, arr.len = 0.1, text.size = 0.8, digits = 3)
  abline(h = 0)
  box(lwd = 1.1)
  axis(side = 1, labels = FALSE, tck = -0.015)
  mtext(text = c("-10","0","10","20"), side = 1, line = 0.4, las = 1, at = c(-10,0,10,20), cex=1)
  axis(side = 2, labels = FALSE, tck = -0.015)
  mtext(text = c("0.0","0.1","0.2"), side = 2, line = 0.4, las = 1, at = c(0,0.1,0.2), cex=1)
  mtext( text = "Difference in expected fitness", side = 1, line = 1.7, cex = axis.label.size)
  mtext( text = "Posterior density", side = 2, line = 2.1, cex = axis.label.size)
  text(labels = c("Low elevation"), x = x.limits[1]-0.035*(x.limits[2]-x.limits[1]), y = y.limits[2]-0.05*(y.limits[2]-y.limits[1]), cex = 1.1, pos = 4)
dev.off()


#
#
###########################################
#
# Tests for local adaptation.
# (1) Within test ponds, do local tads do better than foreign tads?
# (2) Do tadpoles from the four sources do better when tested at home than away?
#
# Test (1) use the four test sites that were sources, but not the other four test sites. Discard those tested at a foreign elevation.
# Test (2) should be done using all test ponds, but discard cases in which a source pop is tested at a foreign elevation.
#
# Test (1) is done on data standardized by test pond. The analysis asks about performance relative to other individuals in the pond.
# Test (2) is done on data standardized by source pond. The analysis asks whether tadpoles of a source do relatively poorly when they are in away sites.
#
dat.for.1   <- subset(d, d$test.pond %in% c("feer", "flue", "bern", "siec"))
dat.for.1   <- subset(dat.for.1, dat.for.1$home.elev == "home_elev")
# Data must be re-standardized within test ponds because the tadpoles at foreign elevation are no longer in the dataset.
dat.for.1   <- rescale.within.groups(dat.for.1, group.list = c("test.pond"), response.list = c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond"), scale = TRUE)
dat.for.1$home.pond <- factor(dat.for.1$home.pond)
dat.for.1$test.elev <- factor(dat.for.1$test.elev)
#
#################################################
#
# Test 1: local versus foreign
#
# MCMCglmm analyses
#
# MCMCglmm analyses of rounded fitness, test 1.
# This is not on standardized fitness because that distribution is shifted left.
# MCMCglmm calls are currently disabled.
#
# recalculate round.fitness for this subset of the dataset.
dat.for.1$round.fitness   <- round(10 * dat.for.1$fitness/mean(dat.for.1$fitness, na.rm = TRUE))
# set the priors.
priors                    <- list(  R=list(V=diag(2), n=2, fix=2), G=list(G1=list(V=diag(c(1, 0.000001)), n=2), G2=list(V=diag(c(1, 0.000001)), n=2, fix=1), G3=list(V=diag(c(1, 0.000001)), n=2, fix=1)  )  )
priors.no.cage            <- list(R=list(V=diag(2), n=2, fix=2), G=list(G1=list(V=diag(c(1, 0.000001)), n=2), G2=list(V=diag(c(1, 0.000001)), n=2, fix=1))  )
priors.no.cage.block      <- list(R=list(V=diag(2), n=2, fix=2), G=list(G1=list(V=diag(c(1, 0.000001)), n=2))  )
priors.no.cage.block.pond <- list(R=list(V=diag(2), n=2, fix=2)  )

set.seed(1)
fit.full.1 <- MCMCglmm(round.fitness ~ trait + trait:test.elev + trait:home.pond + trait:test.elev:home.pond, random = ~idh(trait):test.pond + idh(trait):block:test.pond + idh(trait):cage:block:test.pond, data = dat.for.1, family = "zipoisson", prior = priors, rcov = ~idh(trait):units, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
          # Fit same model with the zeroinfl function from the pscl package, to get the correct Poisson intercept.
fit.full.1.zi <- zeroinfl(round.fitness ~ test.elev + home.pond + test.elev:home.pond, data = dat.for.1)
summary(fit.full.1.zi)
#
set.seed(11)
fit.no.cage.1 <- MCMCglmm(round.fitness ~ trait + trait:test.elev + trait:home.pond + trait:test.elev:home.pond, random = ~idh(trait):test.pond + idh(trait):block:test.pond, data = dat.for.1, family = "zipoisson", prior = priors.no.cage, rcov = ~idh(trait):units, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(111)
fit.no.cage.block.1 <- MCMCglmm(round.fitness ~ trait + trait:test.elev + trait:home.pond + trait:test.elev:home.pond, random = ~idh(trait):test.pond, data = dat.for.1, family = "zipoisson", prior = priors.no.cage.block, rcov = ~idh(trait):units, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(1111)
fit.no.cage.block.pond.1 <- MCMCglmm(round.fitness ~ trait + trait:test.elev + trait:home.pond + trait:test.elev:home.pond, data = dat.for.1, family = "zipoisson", prior = priors.no.cage.block.pond, rcov = ~idh(trait):units, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
# Save the MCMCglmm models for test 1
save(list="fit.full.1", file = "loc.adapt1.fit.full.txt", ascii = TRUE)
save(list="fit.no.cage.1", file = "loc.adapt1.fit.no.cage.txt", ascii = TRUE)
save(list="fit.no.cage.block.1", file = "loc.adapt1.fit.no.cage.block.txt", ascii = TRUE)
save(list="fit.no.cage.block.pond.1", file = "loc.adapt1.fit.no.cage.block.pond.txt", ascii = TRUE)
#
#
# MCMCglmm analyses of survival, test 1 (each model takes 2-3 minutes).
#
priors <- list(R=list(V=1, fix=1), G=list(G1=list(V = 1, nu = 0.002), G2=list(V = 1, nu = 0.002), G3=list(V = 1, nu = 0.002)  )  )
priors.no.cage <- list(R=list(V=1, fix=1), G=list(G1=list(V = 1, nu = 0.002), G2=list(V = 1, nu = 0.002) )  )
priors.no.cage.block <- list(R=list(V=1, fix=1), G=list(G1=list(V = 1, nu = 0.002)  )  )
priors.no.cage.block.pond <- list(R=list(V=1, fix=1))
start.time <- Sys.time()
set.seed(1)
surv.full.1 <- MCMCglmm(survival ~ test.elev*home.pond, random = ~ test.pond + block:test.pond + cage:block:test.pond, data = dat.for.1, family = "categorical", prior = priors, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(11)
surv.no.cage.1 <- MCMCglmm(survival ~ test.elev*home.pond, random = ~ test.pond + block:test.pond, data = dat.for.1, family = "categorical", prior = priors.no.cage, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(111)
surv.no.cage.block.1 <- MCMCglmm(survival ~ test.elev*home.pond, random = ~ test.pond, data = dat.for.1, family = "categorical", prior = priors.no.cage.block, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(1111)
surv.no.cage.block.pond.1 <- MCMCglmm(survival ~ test.elev*home.pond, data = dat.for.1, family = "categorical", prior = priors.no.cage.block.pond, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
# Save the four models for survival, test 1.
save(list="surv.full.1", file = "loc.adapt1.surv.full.txt")
save(list="surv.no.cage.1", file = "loc.adapt1.surv.no.cage.txt")
save(list="surv.no.cage.block.1", file = "loc.adapt1.surv.no.cage.block.txt")
save(list="surv.no.cage.block.pond.1", file = "loc.adapt1.surv.no.cage.block.pond.txt")
#
#
# MCMCglmm analyses of development rate, test 1 (each model takes 30 seconds).
#
priors               <- list(R=list(V = 1e-16, nu = -2), G=list(G1=list(V = 1, nu = 0.002), G2=list(V = 1, nu = 0.002), G3=list(V = 1, nu = 0.002)  )  )
priors.no.cage       <- list(R=list(V = 1e-16, nu = -2), G=list(G1=list(V = 1, nu = 0.002), G2=list(V = 1, nu = 0.002) )  )
priors.no.cage.block <- list(R=list(V = 1e-16, nu = -2), G=list(G1=list(V = 1, nu = 0.002)  )  )
priors.no.cage.block.pond <- list(R=list(V = 1e-16, nu = -2))
dat.dev              <- dat.for.1[ , c("dev.stand.bypond", "test.elev", "home.pond", "block", "cage", "test.pond")]
dat.dev              <- dat.dev[complete.cases(dat.dev), ]
set.seed(1)
dev.full.1 <- MCMCglmm(dev.stand.bypond ~ test.elev*home.pond, random = ~ test.pond + block:test.pond + cage:block:test.pond , data = dat.dev, family = "gaussian", prior = priors, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(11)
dev.no.cage.1 <- MCMCglmm(dev.stand.bypond ~ test.elev*home.pond, random = ~ test.pond + block:test.pond, data = dat.dev, family = "gaussian", prior = priors.no.cage, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(111)
dev.no.cage.block.1 <- MCMCglmm(dev.stand.bypond ~ test.elev*home.pond, random = ~ test.pond , data = dat.dev, family = "gaussian", prior = priors.no.cage.block, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(1111)
dev.no.cage.block.pond.1 <- MCMCglmm(dev.stand.bypond ~ test.elev*home.pond, data = dat.dev, family = "gaussian", prior = priors.no.cage.block.pond, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
# Save four models for development, test 1.
save(list="dev.full.1", file = "loc.adapt1.dev.full.txt")
save(list="dev.no.cage.1", file = "loc.adapt1.dev.no.cage.txt")
save(list="dev.no.cage.block.1", file = "loc.adapt1.dev.no.cage.block.txt")
save(list="dev.no.cage.block.pond.1", file = "loc.adapt1.dev.no.cage.block.pond.txt")
#
#
# MCMCglmm analyses of mass at metamorphosis, test 1 (each model takes 30 seconds).
#
dat.mass <- dat.for.1[ , c("mass.stand.bypond", "test.elev", "home.pond", "block", "cage", "test.pond")]
dat.mass <- dat.mass[complete.cases(dat.mass), ]
set.seed(1)
mass.full.1 <- MCMCglmm(mass.stand.bypond ~ test.elev*home.pond, random = ~ test.pond + block:test.pond + cage:block:test.pond , data = dat.mass, family = "gaussian", prior = priors, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(11)
mass.no.cage.1 <- MCMCglmm(mass.stand.bypond ~ test.elev*home.pond, random = ~ test.pond + block:test.pond, data = dat.mass, family = "gaussian", prior = priors.no.cage, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(111)
mass.no.cage.block.1 <- MCMCglmm(mass.stand.bypond ~ test.elev*home.pond, random = ~ test.pond , data = dat.mass, family = "gaussian", prior = priors.no.cage.block, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(1111)
mass.no.cage.block.pond.1 <- MCMCglmm(mass.stand.bypond ~ test.elev*home.pond, data = dat.mass, family = "gaussian", prior = priors.no.cage.block.pond, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
# Save four models for mass, test 2
save(list="mass.full.1", file = "loc.adapt1.mass.full.txt")
save(list="mass.no.cage.1", file = "loc.adapt1.mass.no.cage.txt")
save(list="mass.no.cage.block.1", file = "loc.adapt1.mass.no.cage.block.txt")
save(list="mass.no.cage.block.pond.1", file = "loc.adapt1.mass.no.cage.block.pond.txt")
#
#
###################
#
# Now produce a summary of the models, for Table 2 (left side) and Table S2 (left side) in the manuscript.
#
make.MCMCglmm.summary(m.full = fit.full.1,  m2 = fit.no.cage.1,  m3 = fit.no.cage.block.1,  m4 = fit.no.cage.block.pond.1)
make.MCMCglmm.summary(m.full = surv.full.1, m2 = surv.no.cage.1, m3 = surv.no.cage.block.1, m4 = surv.no.cage.block.pond.1)
make.MCMCglmm.summary(m.full = dev.full.1,  m2 = dev.no.cage.1,  m3 = dev.no.cage.block.1,  m4 = dev.no.cage.block.pond.1)
make.MCMCglmm.summary(m.full = mass.full.1, m2 = mass.no.cage.1, m3 = mass.no.cage.block.1, m4 = mass.no.cage.block.pond.1)
#
#
###################
#
#  Now test (2): Home versus away criterion: Checking for unique local adaptation within elevation.
#
dat.for.2 <- subset(d, d$home.elev == "home_elev")
dat.for.2$fit.stand.bysource  <- dat.for.2$fitness
dat.for.2$surv.stand.bysource <- dat.for.2$survival
dat.for.2$dev.stand.bysource  <- dat.for.2$DevRate42
dat.for.2$mass.stand.bysource <- dat.for.2$mass.42
dat.for.2                     <- rescale.within.groups(dat.for.2, group.list = c("source.pop"), response.list = c("fit.stand.bysource", "surv.stand.bysource", "dev.stand.bysource", "mass.stand.bysource"), scale = TRUE)
dat.for.2$home.pond           <- factor(dat.for.2$home.pond)
dat.for.2$test.elev           <- factor(dat.for.2$test.elev)
dat.for.2$round.fitness       <- round(10 * dat.for.2$fitness/mean(dat.for.2$fitness, na.rm = TRUE))
#
# MCMCglmm analysis, test2 for local adaptation (home versus away).
# Fitness, test 2 (each model takes 6.5 - 8 minutes).
#
priors                    <- list(R=list(V=diag(2), n=2, fix=2), G=list(G1=list(V=diag(c(1, 0.000001)), n=2), G2=list(V=diag(c(1, 0.000001)), n=2, fix=1), G3=list(V=diag(c(1, 0.000001)), n=2, fix=1)  )  )
priors.no.cage            <- list(R=list(V=diag(2), n=2, fix=2), G=list(G1=list(V=diag(c(1, 0.000001)), n=2), G2=list(V=diag(c(1, 0.000001)), n=2, fix=1))  )
priors.no.cage.block      <- list(R=list(V=diag(2), n=2, fix=2), G=list(G1=list(V=diag(c(1, 0.000001)), n=2))  )
priors.no.cage.block.pond <- list(R=list(V=diag(2), n=2, fix=2)  )
set.seed(1)
fit.full.2 <- MCMCglmm(round.fitness ~ trait + trait:test.elev + trait:home.pond + trait:test.elev:home.pond, random = ~idh(trait):test.pond + idh(trait):block:test.pond + idh(trait):cage:block:test.pond, data = dat.for.2, family = "zipoisson", prior = priors, rcov = ~idh(trait):units, burnin = 100000, thin = 100, nitt   = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
# Fit the same model from pscl package to get the logistic intercept
fit.full.2.zi <- zeroinfl(round.fitness ~ test.elev + home.pond + test.elev:home.pond, data = dat.for.2)
summary(fit.full.2.zi)
#
set.seed(11)
fit.no.cage.2 <- MCMCglmm(round.fitness ~ trait + trait:test.elev + trait:home.pond + trait:test.elev:home.pond, random = ~idh(trait):test.pond + idh(trait):block:test.pond, data = dat.for.2, family = "zipoisson", prior = priors.no.cage, rcov = ~idh(trait):units, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(111)
fit.no.cage.block.2 <- MCMCglmm(round.fitness ~ trait + trait:test.elev + trait:home.pond + trait:test.elev:home.pond, random = ~idh(trait):test.pond, data = dat.for.2, family = "zipoisson", prior = priors.no.cage.block, rcov = ~idh(trait):units, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(1111)
fit.no.cage.block.pond.2 <- MCMCglmm(round.fitness ~ trait + trait:test.elev + trait:home.pond + trait:test.elev:home.pond, data = dat.for.2, family = "zipoisson", prior = priors.no.cage.block.pond, rcov = ~idh(trait):units, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
# Save the four MCMCglmm model objects for fitness, test 2
save(list="fit.full.2", file = "loc.adapt2.fit.full.txt", ascii = TRUE)
save(list="fit.no.cage.2", file = "loc.adapt2.fit.no.cage.txt", ascii = TRUE)
save(list="fit.no.cage.block.2", file = "loc.adapt2.fit.no.cage.block.txt", ascii = TRUE)
save(list="fit.no.cage.block.pond.2", file = "loc.adapt2.fit.no.cage.block.pond.txt", ascii = TRUE)
#
#
# MCMCglmm analyses of survival, test 2 (each model takes 2-3 minutes).
#
priors <- list(R=list(V=1, fix=1), G=list(G1=list(V = 1, nu = 0.002), G2=list(V = 1, nu = 0.002), G3=list(V = 1, nu = 0.002)  )  )
priors.no.cage <- list(R=list(V=1, fix=1), G=list(G1=list(V = 1, nu = 0.002), G2=list(V = 1, nu = 0.002) )  )
priors.no.cage.block <- list(R=list(V=1, fix=1), G=list(G1=list(V = 1, nu = 0.002)  )  )
priors.no.cage.block.pond <- list(R=list(V=1, fix=1))
set.seed(1)
surv.full.2 <- MCMCglmm(survival ~ test.elev*home.pond, random = ~ test.pond + block:test.pond + cage:block:test.pond, data = dat.for.2, family = "categorical", prior = priors, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(11)
surv.no.cage.2 <- MCMCglmm(survival ~ test.elev*home.pond, random = ~ test.pond + block:test.pond, data = dat.for.2, family = "categorical", prior = priors.no.cage, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(111)
surv.no.cage.block.2 <- MCMCglmm(survival ~ test.elev*home.pond, random = ~ test.pond, data = dat.for.2, family = "categorical", prior = priors.no.cage.block, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(1111)
surv.no.cage.block.pond.2 <- MCMCglmm(survival ~ test.elev*home.pond, data = dat.for.2, family = "categorical", prior = priors.no.cage.block.pond, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
# Save the four models for survival, test 2.
save(list="surv.full.2", file = "loc.adapt2.surv.full.txt")
save(list="surv.no.cage.2", file = "loc.adapt2.surv.no.cage.txt")
save(list="surv.no.cage.block.2", file = "loc.adapt2.surv.no.cage.block.txt")
save(list="surv.no.cage.block.pond.2", file = "loc.adapt2.surv.no.cage.block.pond.txt")
#
#
# MCMCglmm analyses of development rate, test 2 (each model takes 30-60 seconds).
#
priors               <- list(R=list(V = 1e-16, nu = -2), G=list(G1=list(V = 1, nu = 0.002), G2=list(V = 1, nu = 0.002), G3=list(V = 1, nu = 0.002)  )  )
priors.no.cage       <- list(R=list(V = 1e-16, nu = -2), G=list(G1=list(V = 1, nu = 0.002), G2=list(V = 1, nu = 0.002) )  )
priors.no.cage.block <- list(R=list(V = 1e-16, nu = -2), G=list(G1=list(V = 1, nu = 0.002)  )  )
priors.no.cage.block.pond <- list(R=list(V = 1e-16, nu = -2))
dat.dev    <- dat.for.2[ , c("dev.stand.bysource", "test.elev", "home.pond", "test.pond", "block", "cage")]
dat.dev    <- dat.dev[complete.cases(dat.dev), ]
#
set.seed(1)
dev.full.2 <- MCMCglmm(dev.stand.bysource ~ test.elev*home.pond, random = ~ test.pond + block:test.pond + cage:block:test.pond, data = dat.dev, family = "gaussian", prior = priors, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(11)
dev.no.cage.2 <- MCMCglmm(dev.stand.bysource ~ test.elev*home.pond, random = ~ test.pond + block:test.pond, data = dat.dev, family = "gaussian", prior = priors.no.cage, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(111)
dev.no.cage.block.2 <- MCMCglmm(dev.stand.bysource ~ test.elev*home.pond, random = ~ test.pond, data = dat.dev, family = "gaussian", prior = priors.no.cage.block, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(1111)
dev.no.cage.block.pond.2 <- MCMCglmm(dev.stand.bysource ~ test.elev*home.pond, data = dat.dev, family = "gaussian", prior = priors.no.cage.block.pond, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
# Save four models for development, test 2.
save(list="dev.full.2", file = "loc.adapt2.dev.full.txt")
save(list="dev.no.cage.2", file = "loc.adapt2.dev.no.cage.txt")
save(list="dev.no.cage.block.2", file = "loc.adapt2.dev.no.cage.block.txt")
save(list="dev.no.cage.block.pond.2", file = "loc.adapt2.dev.no.cage.block.pond.txt")
#
#
# MCMCglmm analyses of mass at metamorphosis, test 2. (each model takes 30-60 seconds)
#
dat.mass <- dat.for.2[ , c("mass.stand.bysource", "test.elev", "home.pond", "test.pond", "block", "cage")]
dat.mass <- dat.mass[complete.cases(dat.mass), ]
#
set.seed(1)
mass.full.2 <- MCMCglmm(mass.stand.bysource ~ test.elev*home.pond, random = ~ test.pond + block:test.pond + cage:block:test.pond , data = dat.mass, family = "gaussian", prior = priors, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(11)
mass.no.cage.2 <- MCMCglmm(mass.stand.bysource ~ test.elev*home.pond, random = ~ test.pond + block:test.pond, data = dat.mass, family = "gaussian", prior = priors.no.cage, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(111)
mass.no.cage.block.2 <- MCMCglmm(mass.stand.bysource ~ test.elev*home.pond, random = ~ test.pond , data = dat.mass, family = "gaussian", prior = priors.no.cage.block, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
set.seed(1111)
mass.no.cage.block.pond.2 <- MCMCglmm(mass.stand.bysource ~ test.elev*home.pond, data = dat.mass, family = "gaussian", prior = priors.no.cage.block.pond, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE)
#
# Save four models for mass, test 2
save(list="mass.full.2", file = "loc.adapt2.mass.full.txt")
save(list="mass.no.cage.2", file = "MCMCglmm results/loc.adapt2.mass.no.cage.txt")
save(list="mass.no.cage.block.2", file = "loc.adapt2.mass.no.cage.block.txt")
save(list="mass.no.cage.block.pond.2", file = "loc.adapt2.mass.no.cage.block.pond.txt")
#
fit.full.2$DIC; fit.no.cage.2$DIC; fit.no.cage.block.2$DIC; fit.no.cage.block.pond.2$DIC       #  DIC values: full model is always best by far
surv.full.2$DIC; surv.no.cage.2$DIC; surv.no.cage.block.2$DIC; surv.no.cage.block.pond.2$DIC
dev.full.2$DIC; dev.no.cage.2$DIC; dev.no.cage.block.2$DIC; dev.no.cage.block.pond.2$DIC
mass.full.2$DIC; mass.no.cage.2$DIC; mass.no.cage.block.2$DIC; mass.no.cage.block.pond.2$DIC
#
make.MCMCglmm.summary(m.full = fit.full.2, m2 = fit.no.cage.2, m3 = fit.no.cage.block.2, m4 = fit.no.cage.block.pond.2)
make.MCMCglmm.summary(surv.full.2, m2 = surv.no.cage.2, m3 = surv.no.cage.block.2, m4 = surv.no.cage.block.pond.2)
make.MCMCglmm.summary(dev.full.2, m2 = dev.no.cage.2, m3 = dev.no.cage.block.2, m4 = dev.no.cage.block.pond.2)
make.MCMCglmm.summary(mass.full.2, m2 = mass.no.cage.2, m3 = mass.no.cage.block.2, m4 = mass.no.cage.block.pond.2)
#
#
#
# End of tests for unique local adaptation
############################







############################
#
# Draw Fig. 3: adaptation to elevation.
# Two panels: one for each of the two kinds of tests.
#
# First, carry out the calculations (means and SEs)
# For Test 1.
dat.for.1           <- subset(d, d$test.pond %in% c("feer", "flue", "bern", "siec"))
dat.for.1           <- subset(dat.for.1, dat.for.1$home.elev == "home_elev")
dat.for.1           <- rescale.within.groups(dat.for.1, group.list = c("test.pond"), response.list = c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond"), scale = TRUE)
cage.means.w.block  <- aggregate(dat.for.1[ , c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond")], by = list(dat.for.1$test.pond, dat.for.1$cage, dat.for.1$source.pop, dat.for.1$home.pond, dat.for.1$block), mean, na.rm = TRUE)
cage.means.w.block  <- rename(cage.means.w.block, c("Group.1", "Group.2", "Group.3", "Group.4", "Group.5"), c("test.pond", "cage", "source.pop", "home.pond", "block"))
home.half           <- subset(cage.means.w.block, cage.means.w.block$home.pond == "home_pond")
home.half           <- rename(home.half, c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond"), c("home.stand.fit", "home.stand.surv", "home.stand.dev", "home.stand.mass"))
away.half           <- subset(cage.means.w.block, cage.means.w.block$home.pond == "away_pond")
away.half           <- rename(away.half, c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond"), c("foreign.stand.fit", "foreign.stand.surv", "foreign.stand.dev", "foreign.stand.mass"))
cage.means.by.block <- merge(home.half, away.half, by = c("test.pond", "block"))
cage.means.by.block$diff.stand.fit  <- cage.means.by.block$home.stand.fit - cage.means.by.block$foreign.stand.fit
cage.means.by.block$diff.stand.surv <- cage.means.by.block$home.stand.surv - cage.means.by.block$foreign.stand.surv
cage.means.by.block$diff.stand.dev  <- cage.means.by.block$home.stand.dev - cage.means.by.block$foreign.stand.dev
cage.means.by.block$diff.stand.mass <- cage.means.by.block$home.stand.mass - cage.means.by.block$foreign.stand.mass
cage.means.by.block          <- cage.means.by.block[ , c(1,2,17:20)]
cage.means.by.pond           <- aggregate(cage.means.by.block[ , -c(1,2)], by = list(cage.means.by.block$test.pond), mean, na.rm = TRUE)
cage.se.by.pond              <- aggregate(cage.means.by.block[ , -c(1,2)], by = list(cage.means.by.block$test.pond), sd, na.rm = TRUE)
cage.se.by.pond[,-1]         <- cage.se.by.pond[,-1] / sqrt(8)
names(cage.se.by.pond)       <- c("pond", "se.stand.fit", "se.stand.surv", "se.stand.dev", "se.stand.mass")
cage.means.by.pond$elev      <- c("high", "low", "high", "low")
names(cage.means.by.pond)[1] <- "pond"
cage.means.by.pond           <- merge(cage.means.by.pond, cage.se.by.pond, by = "pond")
diff.means.for.1             <- cage.means.by.pond[ , c(1,6,2,7,3,8,4,9,5,10)]
cage.means.for.1 <- aggregate(dat.for.1[ , c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond")], by = list(dat.for.1$test.pond, dat.for.1$cage, dat.for.1$source.pop, dat.for.1$home.pond), mean, na.rm = TRUE)
cage.means.for.1 <- rename(cage.means.for.1, c("Group.1", "Group.2", "Group.3", "Group.4"), c("test.pond", "cage", "source.pop", "home.pond"))
pop.means.for.1  <- aggregate(cage.means.for.1[ , c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond")], by = list(cage.means.for.1$test.pond, cage.means.for.1$source.pop, cage.means.for.1$home.pond), mean, na.rm = TRUE)
pop.sd.for.1     <- aggregate(cage.means.for.1[ , c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond")], by = list(cage.means.for.1$test.pond, cage.means.for.1$source.pop, cage.means.for.1$home.pond), sd, na.rm = TRUE)
pop.N.for.1      <- aggregate(cage.means.for.1[ , c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond")], by = list(cage.means.for.1$test.pond, cage.means.for.1$source.pop, cage.means.for.1$home.pond), length)
pop.means.for.1  <- rename(pop.means.for.1, c("Group.1", "Group.2", "Group.3", "fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond"), c("test.pond", "source.pop", "home.pond", "fit.mean", "surv.mean", "dev.mean", "mass.mean"))
pop.sd.for.1     <- rename(pop.sd.for.1, c("Group.1", "Group.2", "Group.3", "fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond"), c("test.pond", "source.pop", "home.pond", "fit.SD", "surv.SD", "dev.SD", "mass.SD"))
pop.N.for.1      <- rename(pop.N.for.1, c("Group.1", "Group.2", "Group.3", "fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond"), c("test.pond", "source.pop", "home.pond", "fit.N", "surv.N", "dev.N", "mass.N"))
pop.means.for.1  <- cbind(pop.means.for.1, pop.sd.for.1[ , -c(1,2,3)], pop.N.for.1[ , -c(1,2,3)])
pop.means.for.1$fit.SE  <- pop.means.for.1$fit.SD / sqrt(pop.means.for.1$fit.N)
pop.means.for.1$surv.SE <- pop.means.for.1$surv.SD / sqrt(pop.means.for.1$surv.N)
pop.means.for.1$dev.SE  <- pop.means.for.1$dev.SD / sqrt(pop.means.for.1$dev.N)
pop.means.for.1$mass.SE <- pop.means.for.1$mass.SD / sqrt(pop.means.for.1$mass.N)
pop.means.for.1         <- pop.means.for.1[ order(pop.means.for.1$test.pond, pop.means.for.1$source.pop), ]
pop.means.for.1$x.seq   <- c(5, 6, 1, 2, 8, 7, 4, 3)
pop.means.for.1         <- pop.means.for.1[ order(pop.means.for.1$x.seq), ]
#
#
# For Test 2.
dat.for.2            <- subset(d, d$home.elev == "home_elev")
dat.for.2$fit.stand.bysource  <- dat.for.2$fitness
dat.for.2$surv.stand.bysource <- dat.for.2$survival
dat.for.2$dev.stand.bysource  <- dat.for.2$DevRate42
dat.for.2$mass.stand.bysource <- dat.for.2$mass.42
dat.for.2 <- rescale.within.groups(dat.for.2, group.list = c("source.pop"), response.list = c("fit.stand.bysource", "surv.stand.bysource", "dev.stand.bysource", "mass.stand.bysource"), scale = TRUE)
cage.means.for.2 <- aggregate(dat.for.2[ , c("fit.stand.bysource", "surv.stand.bysource", "dev.stand.bysource", "mass.stand.bysource")], by = list(dat.for.2$block, dat.for.2$test.pond, dat.for.2$cage, dat.for.2$source.pop, dat.for.2$home.pond), mean, na.rm = TRUE)
cage.means.for.2 <- rename(cage.means.for.2, c("Group.1", "Group.2", "Group.3", "Group.4", "Group.5"), c("block", "test.pond", "cage", "source.pop", "home.pond"))
diff.means.for.2 <- subset(cage.means.for.2, cage.means.for.2$home.pond == "home_pond")
diff.means.for.2$diff.mass <- diff.means.for.2$diff.dev <- diff.means.for.2$diff.surv <- diff.means.for.2$diff.fit <- NA
block.nr <- as.integer(matrix(unlist(strsplit(unique(diff.means.for.2$block), split = "")), ncol = 5, byrow = TRUE)[,5])
for (i in 1:length(block.nr)) {
  away.data  <- subset(cage.means.for.2, cage.means.for.2$source.pop == diff.means.for.2$source.pop[i] & cage.means.for.2$home.pond == "away_pond")
  away.block <- as.integer(matrix(unlist(strsplit(unique(away.data$block), split = "")), ncol = 5, byrow = TRUE)[,5])
  away.data  <- subset(away.data, away.block == block.nr[i])
  for (j in 6:9) {
    diff.means.for.2[i, (j+4)] <- mean(diff.means.for.2[i, j] - away.data[, 6], na.rm = TRUE)
    }
  }
aggregate(diff.means.for.2[ , c("diff.fit", "diff.surv", "diff.dev", "diff.mass")], by = list(diff.means.for.2$source.pop), mean, na.rm = TRUE)
pop.means.for.2  <- aggregate(cage.means.for.2[ , c("fit.stand.bysource", "surv.stand.bysource", "dev.stand.bysource", "mass.stand.bysource")], by = list(cage.means.for.2$test.pond, cage.means.for.2$source.pop, cage.means.for.2$home.pond), mean, na.rm = TRUE)
pop.means.for.2  <- rename(pop.means.for.2, c("Group.1", "Group.2", "Group.3", "fit.stand.bysource", "surv.stand.bysource", "dev.stand.bysource", "mass.stand.bysource"), c("test.pond", "source.pop", "home.pond", "fit.mean", "surv.mean", "dev.mean", "mass.mean"))
diff.means.for.2 <- subset(pop.means.for.2, pop.means.for.2$home.pond == "home_pond")
diff.means.for.2$se.mass <- diff.means.for.2$se.dev <- diff.means.for.2$se.surv <- diff.means.for.2$se.fit <- diff.means.for.2$diff.mass <- diff.means.for.2$diff.dev <- diff.means.for.2$diff.surv <- diff.means.for.2$diff.fit <- NA
pops <- diff.means.for.2$source.pop
for (i in 1:length(pops)) {
  away.data <- subset(pop.means.for.2, pop.means.for.2$source.pop == pops[i] & pop.means.for.2$home.pond == "away_pond")
  for (j in 4:7) {
    diffs <- diff.means.for.2[diff.means.for.2$source.pop == pops[i], j] - away.data[,j]
    diff.means.for.2[diff.means.for.2$source.pop == pops[i], (j+4)] <- mean(diffs)
    diff.means.for.2[diff.means.for.2$source.pop == pops[i], (j+8)] <- sd(diffs) / sqrt(3)
    }
  }
diff.means.for.2 <- diff.means.for.2[ , c(2,8,12,9,13,10,14,11,15)]
pop.sd.for.2     <- aggregate(cage.means.for.2[ , c("fit.stand.bysource", "surv.stand.bysource", "dev.stand.bysource", "mass.stand.bysource")], by = list(cage.means.for.2$test.pond, cage.means.for.2$source.pop, cage.means.for.2$home.pond), sd, na.rm = TRUE)
pop.N.for.2      <- aggregate(cage.means.for.2[ , c("fit.stand.bysource", "surv.stand.bysource", "dev.stand.bysource", "mass.stand.bysource")], by = list(cage.means.for.2$test.pond, cage.means.for.2$source.pop, cage.means.for.2$home.pond), length)
pop.sd.for.2     <- rename(pop.sd.for.2, c("Group.1", "Group.2", "Group.3", "fit.stand.bysource", "surv.stand.bysource", "dev.stand.bysource", "mass.stand.bysource"), c("test.pond", "source.pop", "home.pond", "fit.SD", "surv.SD", "dev.SD", "mass.SD"))
pop.N.for.2      <- rename(pop.N.for.2, c("Group.1", "Group.2", "Group.3", "fit.stand.bysource", "surv.stand.bysource", "dev.stand.bysource", "mass.stand.bysource"), c("test.pond", "source.pop", "home.pond", "fit.N", "surv.N", "dev.N", "mass.N"))
pop.means.for.2  <- cbind(pop.means.for.2, pop.sd.for.2[ , -c(1,2,3)], pop.N.for.2[ , -c(1,2,3)])
pop.means.for.2$fit.SE  <- pop.means.for.2$fit.SD / sqrt(pop.means.for.2$fit.N)
pop.means.for.2$surv.SE <- pop.means.for.2$surv.SD / sqrt(pop.means.for.2$surv.N)
pop.means.for.2$dev.SE  <- pop.means.for.2$dev.SD / sqrt(pop.means.for.2$dev.N)
pop.means.for.2$mass.SE <- pop.means.for.2$mass.SD / sqrt(pop.means.for.2$mass.N)
pop.means.for.2         <- pop.means.for.2[ order(pop.means.for.2$test.pond, pop.means.for.2$source.pop), ]
pop.means.for.2$x.seq   <- c(9, 14, 10, 15, 2, 6, 1, 7, 11, 13, 3, 8, 12, 16, 4, 5)
pop.means.for.2         <- pop.means.for.2[ order(pop.means.for.2$x.seq), ]
#
#
# Draw Fig. S5.
#
margin.limits    <- c(0.1, 0.26, 0.3, 0.46, 0.5, 0.66, 0.7, 0.86)
home.ponds       <- c(1,4,7,10)
away.ponds       <- c(2,5,8,11)
home.ponds.lower <- c(1,6,11,16)
away.ponds.lower <- (1:19)[-c(5,10,15,home.ponds.lower)]
x.axis.label     <- c(1.7, 1.9)  # lines below bottom at which the axis label should be placed, upper and lower panels
top.pond.margin  <- 0.1           # lines above the top of figure at which pond labels should be placed
top.line.margin  <- 0.95          # height above the top of the figure at which line should be drawn
top.label.loc    <- 1.3           # lines above the top of the figure at which upper label should be placed
home.point.size  <- 1.8
away.point.size  <- 1.25
pdf(" .... your path here .... /Fig S5.pdf", useDingbats = FALSE, width = 18, height = 10)
plot.new()
cols <- c("fit.mean", "fit.SE")
  par(new = "TRUE", plt = c(margin.limits[1], margin.limits[2], 0.53, 0.83))   # Top panel local/foreign
  e          <- data.frame(test.pond = pop.means.for.1$test.pond, source = pop.means.for.1$source.pop, x.seq = pop.means.for.1$x.seq, y = pop.means.for.1[ , cols[1]], SE = pop.means.for.1[ , cols[2]], stringsAsFactors = FALSE)
  e$x.seq    <- ifelse(e$x.seq > 6, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 4, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 2, (e$x.seq + 1), e$x.seq)
  y.limits   <- c(min(e$y - e$SE), max(e$y + e$SE))
  plot(x = e$x.seq, y = e$y, xlim = c(0,12), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -1, ybottom = -1, xright = 6, ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = 6, ybottom = -1, xright = 13, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  box(lwd = 1.3)
  abline(v = 3);  abline(v = 6, lwd = 1.5);  abline(v = 9)
  mtext( text = c("Response in SD units"), side = 2, line = 2.3, cex = 1.2)
  # Text indicating test ponds at the top
  mtext( text = c("feer", "siec", "bern", "flue"), side = 3, line = top.pond.margin, cex = 1, at = c(1.05 ,4, 6.8, 9.9), col = c(low.color, low.color, high.color, high.color), adj = 0)
  mtext( text = c("Test site"), side = 3, line = top.label.loc, cex = 1)
  par(xpd = TRUE); lines(x = c(1.1, 11.2), y = c(0.585, 0.585)*top.line.margin, lwd = 1.1); par(xpd = FALSE)
  axis(side = 2, labels = FALSE, tck = -0.015)
  mtext( text = "Expected fitness", side = 3, line = 3, cex = 1.3)
  mtext(text = c("-0.2", "0.0", "0.2", "0.4"), side = 2, line = 0.6, las = 1, at = c(-0.20, 0, 0.2, 0.4), cex = 1)
  par(xpd = TRUE)
  text(x = c(e$x.seq - 0.05), y = -0.41, labels = c(e$source), cex = 0.85, srt = 55, adj = c(1, NA), col = c(rep(low.color,4), rep(high.color,4)))
  par(xpd = FALSE)
  for (i in e$x.seq[1:4]) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = low.color)
  for (i in e$x.seq[5:8]) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = high.color)
  points(x = e$x.seq[e$x.seq %in% home.ponds], y = e$y[e$x.seq %in% home.ponds], bg = c(low.color, low.color, high.color, high.color), col = NA, pch = 21, cex = home.point.size)
  points(x = e$x.seq[e$x.seq %in% away.ponds], y = e$y[e$x.seq %in% away.ponds], bg = c(low.away.color, low.away.color, high.away.color, high.away.color), pch = 21, cex = away.point.size)
  mtext( text = "Source population", side = 1, line = x.axis.label[1], las = 1, cex = 1.2)
  text(labels = "A", x = 11.7, y = 0.43, cex = 1.2)
  # Label at the top
  mtext( text = c("made by 'analysis of Judith transplant.R'"), side = 3, line = 7.4, cex = 0.6)
#
  par(new = "TRUE", plt = c(margin.limits[1], margin.limits[2], 0.11, 0.41))    # Bottom panel home/away.
  e          <- data.frame(test.pond = pop.means.for.2$test.pond, source = pop.means.for.2$source.pop, x.seq = pop.means.for.2$x.seq, y = pop.means.for.2[ , cols[1]], SE = pop.means.for.2[ , cols[2]], stringsAsFactors = FALSE)
  e$x.seq    <- ifelse(e$x.seq > 12, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 8, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 4, (e$x.seq + 1), e$x.seq)
  y.limits   <- c( min(e$y - e$SE), max(e$y + e$SE) )
  plot(x = e$x.seq, y = e$y, xlim = c(0,20), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  # Shade the background by elevation
  rect(xleft = -1, ybottom = (y.limits[1] - 1), xright = 10, ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = 10, ybottom = (y.limits[1] - 1), xright = 21, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  box(lwd = 1.3)
  abline(v = 5);  abline(v = 10, lwd = 1.5);  abline(v = 15)
  mtext( text = c("Response in SD units"), side = 2, line = 2.3, cex = 1.2)
  # Text indicating source ponds at the top
  mtext( text = c("feer", "siec", "bern", "flue"), side = 3, line = top.pond.margin, cex = 1, at = c(1.3 ,6.6, 11.4, 16.4), col = c(low.color, low.color, high.color, high.color), adj = 0)
  mtext( text = c("Source population"), side = 3, line = top.label.loc, cex = 1)
  par(xpd = TRUE); lines(x = c(1.3, 18.4), y = c(1.33, 1.33)*top.line.margin, lwd = 1.1); par(xpd = FALSE)
  axis(side = 2, labels = FALSE, tck = -0.015)
  mtext(text = c("-0.5", "0.0", "0.5", "1.0"), side = 2, line = 0.6, las = 1, at = c(-0.5,0,0.5,1), cex=1)
  par(xpd = TRUE)
  text(x = c(e$x.seq - 0.05), y = -0.665, labels = c(e$test.pond), cex = 0.85, srt = 55, adj = c(1, NA), col = c(rep(low.color,4), rep(high.color,4)))
  par(xpd = FALSE)
  # Error bars first so that they are behind points:
  for (i in e$x.seq[1:8])  lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = low.color)
  for (i in e$x.seq[9:16]) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = high.color)
  # Population means
  points(x = e$x.seq[e$x.seq %in% home.ponds.lower], y = e$y[e$x.seq %in% home.ponds.lower], bg = c(low.color, low.color, high.color, high.color), col = NA, pch = 21, cex = home.point.size)
  points(x = e$x.seq[e$x.seq %in% away.ponds.lower], y = e$y[e$x.seq %in% away.ponds.lower], bg = c(rep(low.away.color, 6), rep(high.away.color, 6)), pch = 21, cex = away.point.size)
  # Label the x-axis.
  mtext( text = "Test site", side = 1, line = x.axis.label[2], las = 1, cex = 1.2)
  text(labels = "E", x = 19.5, y = 1.02, cex = 1.2)
#
# Second block of panels: survival
#
cols <- c("surv.mean", "surv.SE")
  par(new = "TRUE", plt = c(margin.limits[3], margin.limits[4], 0.53, 0.83))
  e          <- data.frame(test.pond = pop.means.for.1$test.pond, source = pop.means.for.1$source.pop, x.seq = pop.means.for.1$x.seq, y = pop.means.for.1[ , cols[1]], SE = pop.means.for.1[ , cols[2]], stringsAsFactors = FALSE)
  e$x.seq    <- ifelse(e$x.seq > 6, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 4, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 2, (e$x.seq + 1), e$x.seq)
  y.limits   <- c( min(e$y - e$SE), max(e$y + e$SE) )
  plot(x = e$x.seq, y = e$y, xlim = c(0,12), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -1, ybottom = -1, xright = 6, ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = 6, ybottom = -1, xright = 13, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  box(lwd = 1.3)
  abline(v = 3);  abline(v = 6, lwd = 1.5);  abline(v = 9)
  # Text indicating test ponds at the top
  mtext( text = c("feer", "siec", "bern", "flue"), side = 3, line = top.pond.margin, cex = 1, at = c(1.05 ,4, 6.8, 9.9), col = c(low.color, low.color, high.color, high.color), adj = 0)
  mtext( text = c("Test site"), side = 3, line = top.label.loc, cex = 1)
  par(xpd = TRUE); lines(x = c(1.1, 11.2), y = c(0.47, 0.47)*top.line.margin, lwd = 1.1); par(xpd = FALSE)
  axis(side = 2, labels = FALSE, tck = -0.015)
  mtext( text = "Survival to metamorphosis", side = 3, line = 3, cex = 1.3)
  mtext(text = c("-0.2", "0.0", "0.2"), side = 2, line = 0.6, las = 1, at = c(-0.20, 0, 0.2), cex = 1)
  par(xpd = TRUE)
  text(x = c(e$x.seq - 0.05), y = -0.375, labels = c(e$source), cex = 0.85, srt = 55, adj = c(1, NA), col = c(rep(low.color,4), rep(high.color,4)))
  par(xpd = FALSE)
  for (i in e$x.seq[1:4]) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = low.color)
  for (i in e$x.seq[5:8]) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = high.color)
  # Population means
  points(x = e$x.seq[e$x.seq %in% home.ponds], y = e$y[e$x.seq %in% home.ponds], bg = c(low.color, low.color, high.color, high.color), col = NA, pch = 21, cex = home.point.size)
  points(x = e$x.seq[e$x.seq %in% away.ponds], y = e$y[e$x.seq %in% away.ponds], bg = c(low.away.color, low.away.color, high.away.color, high.away.color), pch = 21, cex = away.point.size)
  mtext( text = "Source population", side = 1, line = x.axis.label[1], las = 1, cex = 1.2)
  text(labels = "B", x = 11.7, y = 0.351, cex = 1.2)
  mtext( text = c("Figure S5"), side = 3, line = 7.2, cex = 1)
#
  par(new = "TRUE", plt = c(margin.limits[3], margin.limits[4], 0.11, 0.41))
  e          <- data.frame(test.pond = pop.means.for.2$test.pond, source = pop.means.for.2$source.pop, x.seq = pop.means.for.2$x.seq, y = pop.means.for.2[ , cols[1]], SE = pop.means.for.2[ , cols[2]], stringsAsFactors = FALSE)
  e$x.seq    <- ifelse(e$x.seq > 12, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 8, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 4, (e$x.seq + 1), e$x.seq)
  y.limits   <- c( min(e$y - e$SE), max(e$y + e$SE) )
  plot(x = e$x.seq, y = e$y, xlim = c(0,20), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -1, ybottom = (y.limits[1] - 1), xright = 10, ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = 10, ybottom = (y.limits[1] - 1), xright = 21, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  box(lwd = 1.3)
  abline(v = 5);  abline(v = 10, lwd = 1.5);  abline(v = 15)
  # Text indicating source ponds at the top
  mtext( text = c("feer", "siec", "bern", "flue"), side = 3, line = top.pond.margin, cex = 1, at = c(1.3 ,6.6, 11.4, 16.4), col = c(low.color, low.color, high.color, high.color), adj = 0)
  mtext( text = c("Source population"), side = 3, line = top.label.loc, cex = 1)
  par(xpd = TRUE); lines(x = c(1.3, 18.4), y = c(0.553, 0.553)*top.line.margin, lwd = 1.1); par(xpd = FALSE)
  axis(side = 2, labels = FALSE, tck = -0.015, at = c(-0.6, -0.3, 0.0, 0.3))
  mtext(text = c("-0.6", "-0.3", "0.0", "0.3"), side = 2, line = 0.6, las = 1, at = c(-0.6,-0.3,0,0.3), cex=1)
  par(xpd = TRUE); text(x = c(e$x.seq - 0.05), y = -0.78, labels = c(e$test.pond), cex = 0.85, srt = 55, adj = c(1, NA), col = c(rep(low.color,4), rep(high.color,4))); par(xpd = FALSE)
  for (i in e$x.seq[1:8])  lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = low.color)
  for (i in e$x.seq[9:16]) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = high.color)
  points(x = e$x.seq[e$x.seq %in% home.ponds.lower], y = e$y[e$x.seq %in% home.ponds.lower], bg = c(low.color, low.color, high.color, high.color), col = NA, pch = 21, cex = home.point.size)
  points(x = e$x.seq[e$x.seq %in% away.ponds.lower], y = e$y[e$x.seq %in% away.ponds.lower], bg = c(rep(low.away.color, 6), rep(high.away.color, 6)), pch = 21, cex = away.point.size)
  mtext( text = "Test site", side = 1, line = x.axis.label[2], las = 1, cex = 1.2)
  text(labels = "F", x = 19.5, y = 0.353, cex = 1.2)
#
# Third block of panels: development rate
#
cols <- c("dev.mean", "dev.SE")
  par(new = "TRUE", plt = c(margin.limits[5], margin.limits[6], 0.53, 0.83))
  e          <- data.frame(test.pond = pop.means.for.1$test.pond, source = pop.means.for.1$source.pop, x.seq = pop.means.for.1$x.seq, y = pop.means.for.1[ , cols[1]], SE = pop.means.for.1[ , cols[2]], stringsAsFactors = FALSE)
  e$x.seq    <- ifelse(e$x.seq > 6, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 4, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 2, (e$x.seq + 1), e$x.seq)
  y.limits   <- c( min(e$y - e$SE), max(e$y + e$SE) )
  plot(x = e$x.seq, y = e$y, xlim = c(0,12), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -1, ybottom = -1, xright = 6, ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = 6, ybottom = -1, xright = 13, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  box(lwd = 1.3)
  abline(v = 3);  abline(v = 6, lwd = 1.5);  abline(v = 9)
  mtext( text = c("feer", "siec", "bern", "flue"), side = 3, line = top.pond.margin, cex = 1, at = c(1.05 ,4, 6.8, 9.9), col = c(low.color, low.color, high.color, high.color), adj = 0)
  mtext( text = c("Test site"), side = 3, line = top.label.loc, cex = 1)
  par(xpd = TRUE); lines(x = c(1.1, 11.2), y = c(0.784, 0.784)*top.line.margin, lwd = 1.1); par(xpd = FALSE)
  axis(side = 2, labels = FALSE, tck = -0.015)
  mtext(text = "Development rate", side = 3, line = 3, cex = 1.3)
  mtext(text = c("-0.4", "0.0", "0.4"), side = 2, line = 0.6, las = 1, at = c(-0.40, 0, 0.4), cex = 1)
  text(x = c(1.35, 1.35, 1.35, 1.35), y = c(0.36, 0.28, -0.38, -0.46), labels = c("local", "population", "foreign", "population"), cex = c(0.8, 0.8, 0.8, 0.8))
  par(xpd = TRUE); text(x = c(e$x.seq - 0.05), y = -0.81, labels = c(e$source), cex = 0.85, srt = 55, adj = c(1, NA), col = c(rep(low.color,4), rep(high.color,4))); par(xpd = FALSE)
  for (i in e$x.seq[1:4]) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = low.color)
  for (i in e$x.seq[5:8]) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = high.color)
  points(x = e$x.seq[e$x.seq %in% home.ponds], y = e$y[e$x.seq %in% home.ponds], bg = c(low.color, low.color, high.color, high.color), col = NA, pch = 21, cex = home.point.size)
  points(x = e$x.seq[e$x.seq %in% away.ponds], y = e$y[e$x.seq %in% away.ponds], bg = c(low.away.color, low.away.color, high.away.color, high.away.color), pch = 21, cex = away.point.size)
  mtext( text = "Source population", side = 1, line = x.axis.label[1], las = 1, cex = 1.2)
  text(labels = "C", x = 11.7, y = 0.559, cex = 1.2)
  #
  par(new = "TRUE", plt = c(margin.limits[5], margin.limits[6], 0.11, 0.41))
  e          <- data.frame(test.pond = pop.means.for.2$test.pond, source = pop.means.for.2$source.pop, x.seq = pop.means.for.2$x.seq, y = pop.means.for.2[ , cols[1]], SE = pop.means.for.2[ , cols[2]], stringsAsFactors = FALSE)
  e$x.seq    <- ifelse(e$x.seq > 12, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 8, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 4, (e$x.seq + 1), e$x.seq)
  y.limits   <- c( min(e$y - e$SE), max(e$y + e$SE) )
  plot(x = e$x.seq, y = e$y, xlim = c(0,20), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -1, ybottom = (y.limits[1] - 1), xright = 10, ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = 10, ybottom = (y.limits[1] - 1), xright = 21, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  box(lwd = 1.3)
  abline(v = 5);  abline(v = 10, lwd = 1.5);  abline(v = 15)
  mtext(text = c("feer", "siec", "bern", "flue"), side = 3, line = top.pond.margin, cex = 1, at = c(1.3 ,6.6, 11.4, 16.4), col = c(low.color, low.color, high.color, high.color), adj = 0)
  mtext(text = c("Source population"), side = 3, line = top.label.loc, cex = 1)
  par(xpd = TRUE); lines(x = c(1.3, 18.4), y = c(1.74, 1.74)*top.line.margin, lwd = 1.1); par(xpd = FALSE)
  axis(side = 2, labels = FALSE, tck = -0.015, at = c(-1.5,-1,-0.5,0,0.5,1))
  mtext(text = c("-1.5", "-1.0", "-0.5", "0.0", "0.5", "1.0"), side = 2, line = 0.6, las = 1, at = c(-1.5,-1,-0.5,0,0.5,1), cex=1)
  text(x = c(6.7, 6.7, 8, 8), y = c(-0.7, -0.85, 0.90, 0.76), labels = c("home", "site", "away", "sites"), cex = 0.8)
  par(xpd = TRUE); text(x = c(e$x.seq - 0.05), y = -1.75, labels = c(e$test.pond), cex = 0.85, srt = 55, adj = c(1, NA), col = c(rep(low.color,4), rep(high.color,4))); par(xpd = FALSE)
  for (i in e$x.seq[1:8])  lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = low.color)
  for (i in e$x.seq[9:16]) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = high.color)
  points(x = e$x.seq[e$x.seq %in% home.ponds.lower], y = e$y[e$x.seq %in% home.ponds.lower], bg = c(low.color, low.color, high.color, high.color), col = NA, pch = 21, cex = home.point.size)
  points(x = e$x.seq[e$x.seq %in% away.ponds.lower], y = e$y[e$x.seq %in% away.ponds.lower], bg = c(rep(low.away.color, 6), rep(high.away.color, 6)), pch = 21, cex = away.point.size)
  mtext( text = "Test site", side = 1, line = x.axis.label[2], las = 1, cex = 1.2)
  text(labels = "G", x = 19.5, y = 1.225, cex = 1.2)
#
# Fourth block of panels: mass at metamorphosis
#
cols <- c("mass.mean", "mass.SE")
  par(new = "TRUE", plt = c(margin.limits[7], margin.limits[8], 0.53, 0.83))
  e          <- data.frame(test.pond = pop.means.for.1$test.pond, source = pop.means.for.1$source.pop, x.seq = pop.means.for.1$x.seq, y = pop.means.for.1[ , cols[1]], SE = pop.means.for.1[ , cols[2]], stringsAsFactors = FALSE)
  e$x.seq    <- ifelse(e$x.seq > 6, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 4, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 2, (e$x.seq + 1), e$x.seq)
  y.limits   <- c( min(e$y - e$SE), max(e$y + e$SE) )
  plot(x = e$x.seq, y = e$y, xlim = c(0,12), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -1, ybottom = -1, xright = 6, ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = 6, ybottom = -1, xright = 13, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  box(lwd = 1.3)
  abline(v = 3);  abline(v = 6, lwd = 1.5);  abline(v = 9)
  mtext(text = c("feer", "siec", "bern", "flue"), side = 3, line = top.pond.margin, cex = 1, at = c(1.05 ,4, 6.8, 9.9), col = c(low.color, low.color, high.color, high.color), adj = 0)
  mtext(text = c("Test site"), side = 3, line = top.label.loc, cex = 1)
  par(xpd = TRUE); lines(x = c(1.1, 11.2), y = c(0.544, 0.544)*top.line.margin, lwd = 1.1); par(xpd = FALSE)
  axis(side = 2, labels = FALSE, tck = -0.015)
  mtext( text = "Mass at metamorphosis", side = 3, line = 3, cex = 1.3)
  mtext(text = c("-0.4", "0.0", "0.4"), side = 2, line = 0.6, las = 1, at = c(-0.40, 0, 0.4), cex = 1)
  par(xpd = TRUE); text(x = c(e$x.seq - 0.05), y = -0.62, labels = c(e$source), cex = 0.85, srt = 55, adj = c(1, NA), col = c(rep(low.color,4), rep(high.color,4))); par(xpd = FALSE)
  for (i in e$x.seq[1:4]) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = low.color)
  for (i in e$x.seq[5:8]) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = high.color)
  points(x = e$x.seq[e$x.seq %in% home.ponds], y = e$y[e$x.seq %in% home.ponds], bg = c(low.color, low.color, high.color, high.color), col = NA, pch = 21, cex = home.point.size)
  points(x = e$x.seq[e$x.seq %in% away.ponds], y = e$y[e$x.seq %in% away.ponds], bg = c(low.away.color, low.away.color, high.away.color, high.away.color), pch = 21, cex = away.point.size)
  mtext(text = "Source population", side = 1, line = x.axis.label[1], las = 1, cex = 1.2)
  par(xpd = TRUE); text(labels = c("Local-vs-foreign", "comparison"), x = c(15, 13.8), y = c(-0.05, -0.05), srt = -90, cex = 1.3); par(xpd = FALSE)
  text(labels = "D", x = 11.7, y = 0.372, cex = 1.2)
  #
  par(new = "TRUE", plt = c(margin.limits[7], margin.limits[8], 0.11, 0.41))
  e          <- data.frame(test.pond = pop.means.for.2$test.pond, source = pop.means.for.2$source.pop, x.seq = pop.means.for.2$x.seq, y = pop.means.for.2[ , cols[1]], SE = pop.means.for.2[ , cols[2]], stringsAsFactors = FALSE)
  e$x.seq    <- ifelse(e$x.seq > 12, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 8, (e$x.seq + 1), e$x.seq)
  e$x.seq    <- ifelse(e$x.seq > 4, (e$x.seq + 1), e$x.seq)
  y.limits   <- c( min(e$y - e$SE), max(e$y + e$SE) )
  plot(x = e$x.seq, y = e$y, xlim = c(0,20), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -1, ybottom = (y.limits[1] - 1), xright = 10, ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = 10, ybottom = (y.limits[1] - 1), xright = 21, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  box(lwd = 1.3)
  abline(v = 5);  abline(v = 10, lwd = 1.5);  abline(v = 15)
  mtext(text = c("feer", "siec", "bern", "flue"), side = 3, line = top.pond.margin, cex = 1, at = c(1.3 ,6.6, 11.4, 16.4), col = c(low.color, low.color, high.color, high.color), adj = 0)
  mtext(text = c("Source population"), side = 3, line = top.label.loc, cex = 1)
  par(xpd = TRUE); lines(x = c(1.3, 18.4), y = c(1.605, 1.605)*top.line.margin, lwd = 1.1); par(xpd = FALSE)
  axis(side = 2, labels = FALSE, tck = -0.015)
  mtext(text = c("-1.0", "-0.5", "0.0", "0.5", "1.0"), side = 2, line = 0.6, las = 1, at = c(-1,-0.5,0,0.5,1), cex=1)
  par(xpd = TRUE); text(x = c(e$x.seq - 0.05), y = -1.23, labels = c(e$test.pond), cex = 0.85, srt = 55, adj = c(1, NA), col = c(rep(low.color,4), rep(high.color,4))); par(xpd = FALSE)
  for (i in e$x.seq[1:8])  lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = low.color)
  for (i in e$x.seq[9:16]) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = high.color)
  points(x = e$x.seq[e$x.seq %in% home.ponds.lower], y = e$y[e$x.seq %in% home.ponds.lower], bg = c(low.color, low.color, high.color, high.color), col = NA, pch = 21, cex = home.point.size)
  points(x = e$x.seq[e$x.seq %in% away.ponds.lower], y = e$y[e$x.seq %in% away.ponds.lower], bg = c(rep(low.away.color, 6), rep(high.away.color, 6)), pch = 21, cex = away.point.size)
  mtext(text = "Test site", side = 1, line = x.axis.label[2], las = 1, cex = 1.2)
  par(xpd = TRUE); text(labels = c("Home-vs-away", "comparison"), x = c(25, 23), y = c(0.05, 0.05), srt = -90, cex = 1.3); par(xpd = FALSE)
  text(labels = "H", x = 19.5, y = 1.2, cex = 1.2)
dev.off()
#
#
# Draw Fig. 3.
# Currently, this prints only fitness. Would be easy to modify for other responses.
#
# First, recalculate treatment means, centered by test.pond and then returned to the original measurement units.
dat.for.1   <- subset(d, d$test.pond %in% c("feer", "flue", "bern", "siec"))
dat.for.1   <- subset(dat.for.1, dat.for.1$home.elev == "home_elev")
dat.for.1[ , c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond")] <- dat.for.1[,c("fitness", "survival", "DevRate42", "mass.42")]
dat.for.1   <- rescale.within.groups(dat.for.1, group.list = c("test.pond"), response.list = c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond"), scale = FALSE)
dat.for.1   <- rename(dat.for.1, c("fit.s.bytestpond", "surv.s.bytestpond", "dev.s.bytestpond", "mass.s.bytestpond"), c("fit", "surv", "dev", "mass"))
grand.means <- apply(dat.for.1[,c("fitness", "survival", "DevRate42", "mass.42")], 2, mean, na.rm = TRUE)
grand.sds   <- apply(dat.for.1[,c("fitness", "survival", "DevRate42", "mass.42")], 2, sd, na.rm = TRUE)
dat.for.1[ , c("fit")]  <- dat.for.1[ , c("fit")] + grand.means[1]
dat.for.1[ , c("surv")] <- dat.for.1[ , c("surv")] + grand.means[2]
dat.for.1[ , c("dev")]  <- dat.for.1[ , c("dev")] + grand.means[3]
dat.for.1[ , c("mass")] <- dat.for.1[ , c("mass")] + grand.means[4]
d2 <- aggregate(dat.for.1[,c("fit", "surv", "dev", "mass")], by = list(dat.for.1$source.elev, dat.for.1$home.pond, dat.for.1$source.pop, dat.for.1$test.pond, dat.for.1$block, dat.for.1$cage), mean, na.rm = TRUE)
d3 <- aggregate(d2[,c("fit", "surv", "dev", "mass")], by = list(d2$Group.1, d2$Group.2, d2$Group.3, d2$Group.4), mean, na.rm = TRUE)
d4.means <- aggregate(d3[,c("fit", "surv", "dev", "mass")], by = list(d3$Group.1, d3$Group.2), mean, na.rm = TRUE); names(d4.means)[1:2] <- c("test.elev", "local.foreign")
d4.sd    <- aggregate(d3[,c("fit", "surv", "dev", "mass")], by = list(d3$Group.1, d3$Group.2), sd, na.rm = TRUE)
d4.sd[,-c(1:2)] <- d4.sd[,-c(1:2)] / sqrt(2); names(d4.sd) <- c("test.elev", "local.foreign", "fit.SE", "surv.SE", "dev.SE", "mass.SE")
d5 <- merge(d4.means, d4.sd, by = c("test.elev", "local.foreign"))
d5$local.foreign <- c("foreign", "local", "foreign", "local"); d5$x.seq <- c(2.1, 2.0, 1.1, 1.0); d.test.1 <- d5[,c(2,1,3,7,4,8,5,9,6,10,11)]
d.test.1 <- d.test.1[order(d.test.1$x.seq),]
d.test.1   # treatment means and SEs (n = 2), centered by test pond.
#
dat.for.2   <- subset(d, d$home.elev == "home_elev")
dat.for.2$fit.stand.bysource <- dat.for.2$fitness; dat.for.2$surv.stand.bysource <- dat.for.2$survival; dat.for.2$dev.stand.bysource  <- dat.for.2$DevRate42; dat.for.2$mass.stand.bysource <- dat.for.2$mass.42
dat.for.2   <- rescale.within.groups(dat.for.2, group.list = c("source.pop"), response.list = c("fit.stand.bysource", "surv.stand.bysource", "dev.stand.bysource", "mass.stand.bysource"), scale = FALSE)
dat.for.2   <- rename(dat.for.2, c("fit.stand.bysource", "surv.stand.bysource", "dev.stand.bysource", "mass.stand.bysource"), c("fit", "surv", "dev", "mass"))
grand.means <- apply(dat.for.2[,c("fitness", "survival", "DevRate42", "mass.42")], 2, mean, na.rm = TRUE)
grand.sds   <- apply(dat.for.2[,c("fitness", "survival", "DevRate42", "mass.42")], 2, sd, na.rm = TRUE)
dat.for.2[ , c("fit")]  <- dat.for.2[ , c("fit")] + grand.means[1]
dat.for.2[ , c("surv")] <- dat.for.2[ , c("surv")] + grand.means[2]
dat.for.2[ , c("dev")]  <- dat.for.2[ , c("dev")] + grand.means[3]
dat.for.2[ , c("mass")] <- dat.for.2[ , c("mass")] + grand.means[4]
d2 <- aggregate(dat.for.2[,c("fit", "surv", "dev", "mass")], by = list(dat.for.2$source.elev, dat.for.2$home.pond, dat.for.2$source.pop, dat.for.2$test.pond, dat.for.2$block, dat.for.2$cage), mean, na.rm = TRUE)
d3 <- aggregate(d2[,c("fit", "surv", "dev", "mass")], by = list(d2$Group.1, d2$Group.2, d2$Group.3, d2$Group.4), mean, na.rm = TRUE)
d4 <- aggregate(d3[,c("fit", "surv", "dev", "mass")], by = list(d3$Group.1, d3$Group.2, d3$Group.3), mean, na.rm = TRUE)
d5.means <- aggregate(d4[,c("fit", "surv", "dev", "mass")], by = list(d4$Group.1, d4$Group.2), mean, na.rm = TRUE); names(d5.means)[1:2] <- c("source.elev", "home.pond")
d5.sd    <- aggregate(d4[,c("fit", "surv", "dev", "mass")], by = list(d4$Group.1, d4$Group.2), sd, na.rm = TRUE)
d5.sd[,-c(1:2)] <- d5.sd[,-c(1:2)]/sqrt(2); names(d5.sd) <- c("source.elev", "home.pond", "fit.SE", "surv.SE", "dev.SE", "mass.SE")
d.test.2 <- merge(d5.means, d5.sd, by = c("source.elev", "home.pond"))
d.test.2$home.away <- c("away", "home", "away", "home"); d.test.2 <- rename(d.test.2, "source.elev", "test.elev"); d.test.2$x.seq <- c(2.1, 2.0, 1.1, 1.0); d.test.2 <- d.test.2[, c(11,1,3,7,4,8,5,9,6,10,12)]
d.test.2 <- d.test.2[order(d.test.2$x.seq),]
d.test.2   # treatment means and SEs (n = 2 source populations), centered by test pond.
#
axis.label.size <- 1.1
tick.label.size <- 1
home.pt.size    <- 2
away.pt.size    <- 1.35
#
#
pdf(" .... your path here .... /Fig 3.pdf", width = 8, height = 7)
plot.new()
#
# Left panel is test (1) from above: fitness within ponds, local versus foreign criterion.
#
  par(new = "TRUE", plt = c(0.15, 0.45, 0.25, 0.75))
  e <- d.test.1[ , c("local.foreign", "test.elev", "x.seq", "fit", "fit.SE")]
  names(e)[4:5] <- c("y", "SE")
  y.limits <- c( min(e$y - e$SE), max(e$y + e$SE) )
  plot(x = e$x.seq, y = e$y, xlim = c(0.6, 2.5), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = 0, ybottom = -1, xright = mean(e$x.seq), ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = mean(e$x.seq), ybottom = -1, xright = 4, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  abline(v = mean(e$x.seq), lwd = 1.3);  abline(h = 0, lty = 2)
  box(lwd = 1.3)
  # Label the panel at the top
  mtext(text = c("A. Local vs. foreign criterion"), side = 3, line = 0.4, cex = 1.10, at = 0.55, adj = 0)
  # Horizontal axis on the bottom
  mtext( text = "Elevation of test site", side = 1, line = 1.5, cex = 1.05)
  mtext( text = c("Low", "High"), side = 1, line = 0.3, cex = axis.label.size, at = c(mean(e$x.seq[c(1,2)]) , mean(e$x.seq[c(3,4)])), col = c(low.color, high.color))
  # Vertical axis on the left
  axis(side = 2, labels = FALSE, tck = -0.015)
  mtext(text = "Expected fitness", side = 2, line = 2.4, las = 3, cex = axis.label.size)
  mtext(text = c("0.16", "0.20", "0.24", "0.28"), side = 2, line = 0.4, las = 1, at = c(0.16, 0.2, 0.24, 0.28), cex = tick.label.size)
  text(labels = c("local", "foreign"), x = c(0.99, 1.11), y = c(e$y[c(1,2)]), cex = 0.8, pos = c(2,2))
  # Error bars behind points, and then the points:
  for (i in 1:dim(e)[1]) {
    color <- ifelse(e$test.elev[i] == "low", low.color, high.color)
    lines(x = rep(e$x.seq[i], 2), y = c(e$y[i] - e$SE[i], e$y[i] + e$SE[i]), col = color, lwd = 1.3)  }
  points(x = e$x.seq, y = e$y, bg = c(low.color, low.away.color, high.color, high.away.color), col = c(NA, "black", NA, "black"), pch = 21, cex = c(home.pt.size, away.pt.size, home.pt.size, away.pt.size))
  # Labels at the very top
  mtext(text = c("Figure 3"), side = 3, line = 3.8, cex = 0.9)
#
# Right panel is test (2) from above: fitness within ponds, home versus away criterion.
#
  par(new = "TRUE", plt = c(0.58, 0.88, 0.25, 0.75))
  e <- d.test.2[ , c("home.away", "test.elev", "x.seq", "fit", "fit.SE")]
  names(e)[4:5] <- c("y", "SE")
  y.limits <- c( min(e$y - e$SE), max(e$y + e$SE) )
  plot(x = e$x.seq, y = e$y, xlim = c(0.6, 2.5), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = 0, ybottom = -1, xright = mean(e$x.seq), ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = mean(e$x.seq), ybottom = -1, xright = 4, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  abline(v = mean(e$x.seq), lwd = 1.3);  abline(h = 0, lty = 2)
  box(lwd = 1.3)
  # Label the panel at the top
  mtext(text = c("B. Home vs. away criterion"), side = 3, line = 0.4, cex = 1.10, at = 0.55, adj = 0)
  # Horizontal axis on the bottom
  mtext( text = "Elevation of source population", side = 1, line = 1.5, cex = 1.05)
  mtext( text = c("Low", "High"), side = 1, line = 0.3, cex = axis.label.size, at = c(mean(e$x.seq[c(1,2)]) , mean(e$x.seq[c(3,4)])), col = c(low.color, high.color))
  # Vertical axis on the left
  axis(side = 2, labels = FALSE, tck = -0.015)
  mtext( text = "Expected fitness", side = 2, line = 2.4, las = 3, cex = axis.label.size)
  mtext(text = c("0.22", "0.24", "0.26", "0.28"), side = 2, line = 0.4, las = 1, at = c(0.22, 0.24, 0.26, 0.28), cex = tick.label.size)
  text(labels = c("home", "away"), x = c(0.99, 1.09), y = c(e$y[c(1,2)]) + 0.002, cex = 0.8, pos = c(2,4))
  text(labels = c("pond", "pond"), x = c(0.99, 1.09), y = c(e$y[c(1,2)]) - 0.002, cex = 0.8, pos = c(2,4))
  # Error bars first so that they are behind points:
  for (i in 1:dim(e)[1]) {
    color <- ifelse(e$test.elev[i] == "low", low.color, high.color)
    lines(x = rep(e$x.seq[i], 2), y = c(e$y[i] - e$SE[i], e$y[i] + e$SE[i]), col = color, lwd = 1.2)
    }
  # Population mean values of home - foreign difference.
  points(x = e$x.seq, y = e$y, bg = c(low.color, low.away.color, high.color, high.away.color), col = c(NA, "black", NA, "black"), pch = 21, cex = c(home.pt.size, away.pt.size, home.pt.size, away.pt.size))
  mtext(text = c("made by 'analysis of Judith transplant.R'"), side = 3, line = 4, cex = 0.6)
dev.off()



#
#
#
########################################
# Adaptation to elevation.
#
d <- rename(d, c("fit.stand", "surv.stand", "dev.stand", "mass.stand"), c("fit.s.overall", "surv.s.overall", "dev.s.overall", "mass.s.overall"))
#
# MCMCglmm analysis on expected fitness.
# Takes about 14-18 minutes for each model on fitness.
#
priors <- list(
  R=list(V=diag(2), n=2, fix=2),
  G=list(G1=list(V=diag(c(1, 0.000001)), n=2),
         G2=list(V=diag(c(1, 0.000001)), n=2, fix=1),
         G3=list(V=diag(c(1, 0.000001)), n=2, fix=1)  )  )
priors.no.cage            <- list(R=list(V=diag(2), n=2, fix=2), G=list(G1=list(V=diag(c(1, 0.000001)), n=2), G2=list(V=diag(c(1, 0.000001)), n=2, fix=1))  )
priors.no.cage.block      <- list(R=list(V=diag(2), n=2, fix=2), G=list(G1=list(V=diag(c(1, 0.000001)), n=2))  )
priors.no.cage.block.pond <- list(R=list(V=diag(2), n=2, fix=2)  )
set.seed(1)
fit.full <- MCMCglmm(round.fitness ~ trait + trait:test.elev + trait:source.elev + trait:test.elev:source.elev, random = ~idh(trait):test.pond + idh(trait):block:test.pond + idh(trait):cage:block:test.pond, data = d, family = "zipoisson", prior = priors, rcov = ~idh(trait):units, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
# Fit same model from pscl package to get the logistic intercept
fit.full.zi <- zeroinfl(round.fitness ~ test.elev + source.elev + test.elev:source.elev, data = d)
summary(fit.full.zi)
#
set.seed(11)
fit.no.cage <- MCMCglmm(round.fitness ~ trait + trait:test.elev + trait:source.elev + trait:test.elev:source.elev, random = ~idh(trait):test.pond + idh(trait):block:test.pond, data = d, family = "zipoisson", prior = priors.no.cage, rcov = ~idh(trait):units, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
set.seed(111)
fit.no.cage.block <- MCMCglmm(round.fitness ~ trait + trait:test.elev + trait:source.elev + trait:test.elev:source.elev, random = ~idh(trait):test.pond, data = d, family = "zipoisson", prior = priors.no.cage.block, rcov = ~idh(trait):units, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
set.seed(1111)
fit.no.cage.block.pond <- MCMCglmm(round.fitness ~ trait + trait:test.elev + trait:source.elev + trait:test.elev:source.elev, data = d, family = "zipoisson", prior = priors.no.cage.block.pond, rcov = ~idh(trait):units, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
# Save all four models
save(list="fit.full", file = "elev.adapt.fit.full.txt", ascii = TRUE)
save(list="fit.no.cage", file = "elev.adapt.fit.no.cage.txt", ascii = TRUE)
save(list="fit.no.cage.block", file = "elev.adapt.fit.no.cage.block.txt", ascii = TRUE)
save(list="fit.no.cage.block.pond", file = "elev.adapt.fit.no.cage.block.pond.txt", ascii = TRUE)
#
#
# MCMCglmm analyses of survival (about 5-6 minutes for each model)
#
priors <- list(R=list(V=1, fix=1), G=list(G1=list(V = 1, nu = 0.002), G2=list(V = 1, nu = 0.002), G3=list(V = 1, nu = 0.002)  )  )
priors.no.cage <- list(R=list(V=1, fix=1), G=list(G1=list(V = 1, nu = 0.002), G2=list(V = 1, nu = 0.002) )  )
priors.no.cage.block <- list(R=list(V=1, fix=1), G=list(G1=list(V = 1, nu = 0.002)  )  )
priors.no.cage.block.pond <- list(R=list(V=1, fix=1))
#
set.seed(1)
surv.full <- MCMCglmm(survival ~ test.elev*source.elev, random = ~ test.pond + block:test.pond + cage:block:test.pond , data = d, family = "categorical", prior = priors, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
set.seed(11)
surv.no.cage <- MCMCglmm(survival ~ test.elev*source.elev, random = ~ test.pond + block:test.pond, data = d, family = "categorical", prior = priors.no.cage, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
set.seed(111)
surv.no.cage.block <- MCMCglmm(survival ~ test.elev*source.elev, random = ~ test.pond , data = d, family = "categorical", prior = priors.no.cage.block, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
set.seed(1111)
surv.no.cage.block.pond <- MCMCglmm(survival ~ test.elev*source.elev, data = d, family = "categorical", prior = priors.no.cage.block.pond, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
# Save all four models
save(list="surv.full", file = "elev.adapt.surv.full.txt", ascii = TRUE)
save(list="surv.no.cage", file = "elev.adapt.surv.no.cage.txt", ascii = TRUE)
save(list="surv.no.cage.block", file = "elev.adapt.surv.no.cage.block.txt", ascii = TRUE)
save(list="surv.no.cage.block.pond", file = "elev.adapt.surv.no.cage.block.pond.txt", ascii = TRUE)
#
#
# MCMCglmm analyses of development rate (about 1-2 minutes each model)
#
priors               <- list(R=list(V = 1e-16, nu = -2), G=list(G1=list(V = 1, nu = 0.002), G2=list(V = 1, nu = 0.002), G3=list(V = 1, nu = 0.002)  )  )
priors.no.cage       <- list(R=list(V = 1e-16, nu = -2), G=list(G1=list(V = 1, nu = 0.002), G2=list(V = 1, nu = 0.002) )  )
priors.no.cage.block <- list(R=list(V = 1e-16, nu = -2), G=list(G1=list(V = 1, nu = 0.002)  )  )
priors.no.cage.block.pond <- list(R=list(V = 1e-16, nu = -2))
dat.dev <- d[ , c("DevRate42", "test.elev", "source.elev", "test.pond", "block", "cage")]
dat.dev <- dat.dev[ complete.cases(dat.dev), ]
#
set.seed(1)
dev.full <- MCMCglmm(DevRate42 ~ test.elev*source.elev, random = ~ test.pond + block:test.pond + cage:block:test.pond , data = dat.dev, family = "gaussian", prior = priors, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
set.seed(11)
dev.no.cage <- MCMCglmm(DevRate42 ~ test.elev*source.elev, random = ~ test.pond + block:test.pond, data = dat.dev, family = "gaussian", prior = priors.no.cage, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
set.seed(111)
dev.no.cage.block <- MCMCglmm(DevRate42 ~ test.elev*source.elev, random = ~ test.pond , data = dat.dev, family = "gaussian", prior = priors.no.cage.block, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
set.seed(1111)
dev.no.cage.block.pond <- MCMCglmm(DevRate42 ~ test.elev*source.elev, data = dat.dev, family = "gaussian", prior = priors.no.cage.block.pond, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
# Save all four development models
save(list="dev.full", file = "elev.adapt.dev.full.txt", ascii = TRUE)
save(list="dev.no.cage", file = "elev.adapt.dev.no.cage.txt", ascii = TRUE)
save(list="dev.no.cage.block", file = "elev.adapt.dev.no.cage.block.txt", ascii = TRUE)
save(list="dev.no.cage.block.pond", file = "elev.adapt.dev.no.cage.block.pond.txt", ascii = TRUE)
#
#
# MCMCglmm analyses of mass at metamorphosis (1-2 minutes each)
#
dat.mass <- d[ , c("mass.42", "test.elev", "source.elev", "test.pond", "block", "cage")]
dat.mass <- dat.mass[ complete.cases(dat.mass), ]
#
set.seed(1)
mass.full <- MCMCglmm(mass.42 ~ test.elev*source.elev, random = ~ test.pond + block:test.pond + cage:block:test.pond , data = dat.mass, family = "gaussian", prior = priors, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
set.seed(11)
mass.no.cage <- MCMCglmm(mass.42 ~ test.elev*source.elev, random = ~ test.pond + block:test.pond, data = dat.mass, family = "gaussian", prior = priors.no.cage, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
set.seed(111)
mass.no.cage.block <- MCMCglmm(mass.42 ~ test.elev*source.elev, random = ~ test.pond , data = dat.mass, family = "gaussian", prior = priors.no.cage.block, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
set.seed(1111)
mass.no.cage.block.pond <- MCMCglmm(mass.42 ~ test.elev*source.elev, data = dat.mass, family = "gaussian", prior = priors.no.cage.block.pond, burnin = 100000, thin = 100, nitt = 300000, pl = FALSE, pr = TRUE, verbose = FALSE )
#
# Save all four mass models
save(list="mass.full", file = "elev.adapt.mass.full.txt", ascii = TRUE)
save(list="mass.no.cage", file = "elev.adapt.mass.no.cage.txt", ascii = TRUE)
save(list="mass.no.cage.block", file = "elev.adapt.mass.no.cage.block.txt", ascii = TRUE)
save(list="mass.no.cage.block.pond", file = "elev.adapt.mass.no.cage.block.pond.txt", ascii = TRUE)
#
# Full model is always best by far.
fit.full$DIC;  fit.no.cage$DIC;  fit.no.cage.block$DIC;  fit.no.cage.block.pond$DIC
surv.full$DIC; surv.no.cage$DIC; surv.no.cage.block$DIC; surv.no.cage.block.pond$DIC
dev.full$DIC;  dev.no.cage$DIC;  dev.no.cage.block$DIC;  dev.no.cage.block.pond$DIC
mass.full$DIC; mass.no.cage$DIC; mass.no.cage.block$DIC; mass.no.cage.block.pond$DIC
#
# Model summaries for Table 1 in the manuscript.
make.MCMCglmm.summary(m.full = fit.full, m2 = fit.no.cage, m3 = fit.no.cage.block, m4 = fit.no.cage.block.pond)
make.MCMCglmm.summary(surv.full, m2 = surv.no.cage, m3 = surv.no.cage.block, m4 = surv.no.cage.block.pond)
make.MCMCglmm.summary(dev.full, m2 = dev.no.cage, m3 = dev.no.cage.block, m4 = dev.no.cage.block.pond)
make.MCMCglmm.summary(mass.full, m2 = mass.no.cage, m3 = mass.no.cage.block, m4 = mass.no.cage.block.pond)
#
#
######################
#
# Cage means, pond means, and treatment means -- needed for the figures that follow.
#
# First center data by test.pond (within test elevations) and then add back the grand mean for each elevation
d.low        <- subset(d, d$test.elev == "low")
fit.mean.low <- mean(d.low$fitness, na.rm=TRUE); surv.mean.low <- mean(d.low$survival, na.rm=TRUE); dev.mean.low <- mean(d.low$DevRate42, na.rm=TRUE); mass.mean.low <- mean(d.low$mass.42, na.rm=TRUE)
d.low$c.fit  <- d.low$fitness; d.low$c.dev <- d.low$DevRate42; d.low$c.mass <- d.low$mass.42; d.low$c.surv <- d.low$survival
d.low        <- rescale.within.groups(d.low, group.list = c("test.pond"), response.list = c("c.fit", "c.dev", "c.mass", "c.surv"), scale = FALSE)
d.low$c.fit  <- d.low$c.fit + fit.mean.low; d.low$c.dev <- d.low$c.dev + dev.mean.low; d.low$c.mass <- d.low$c.mass + mass.mean.low; d.low$c.surv <- d.low$c.surv + surv.mean.low
d.hi        <- subset(d, d$test.elev == "high")
fit.mean.hi <- mean(d.hi$fitness, na.rm=TRUE); surv.mean.hi <- mean(d.hi$survival, na.rm=TRUE); dev.mean.hi <- mean(d.hi$DevRate42, na.rm=TRUE); mass.mean.hi <- mean(d.hi$mass.42, na.rm=TRUE)
d.hi$c.fit  <- d.hi$fitness; d.hi$c.dev <- d.hi$DevRate42; d.hi$c.mass <- d.hi$mass.42; d.hi$c.surv <- d.hi$survival
d.hi        <- rescale.within.groups(d.hi, group.list = c("test.pond"), response.list = c("c.fit", "c.dev", "c.mass", "c.surv"), scale = FALSE)
d.hi$c.fit  <- d.hi$c.fit + fit.mean.hi; d.hi$c.dev <- d.hi$c.dev + dev.mean.hi; d.hi$c.mass <- d.hi$c.mass + mass.mean.hi; d.hi$c.surv <- d.hi$c.surv + surv.mean.hi
d.cent      <- rbind(d.low, d.hi)
#
cage.means <- aggregate(d.cent[ , c("c.fit", "c.dev", "c.mass", "c.surv")], by = list(d.cent$test.pond, d.cent$source.pop, d.cent$home.pond, d.cent$test.elev, d.cent$source.elev, d.cent$cage), mean, na.rm = TRUE)
names(cage.means)[1:6] <- c("test.pond", "source.pop", "home.pond", "test.elev", "source.elev", "cage")
cage.means <- cage.means[order(cage.means$test.pond, cage.means$cage), ]
head(cage.means, n=30)
#
# Create test.pond.means. This is used for making an appendix figure showing all the individual sources in the eight ponds.
testpond.means     <- aggregate(cage.means[ , c("c.fit", "c.surv", "c.dev", "c.mass")], by = list(cage.means$test.pond, cage.means$source.pop, cage.means$home.pond, cage.means$test.elev, cage.means$source.elev), mean, na.rm = TRUE)
testpond.means     <- rename(testpond.means, c("Group.1", "Group.2", "Group.3", "Group.4", "Group.5"), c("test.pond", "source", "home.pond", "test.elev", "source.elev") )
testpond.N         <- aggregate(cage.means[ , c("c.fit", "c.surv", "c.dev", "c.mass")], by = list(cage.means$test.pond, cage.means$source.pop, cage.means$home.pond, cage.means$test.elev, cage.means$source.elev), length)
names(testpond.N)  <- c("test.pond", "source", "home.pond", "test.elev", "source.elev", "fit.N", "surv.N", "dev.N", "N")
testpond.sd        <- aggregate(cage.means[ , c("c.fit", "c.surv", "c.dev", "c.mass")], by = list(cage.means$test.pond, cage.means$source.pop, cage.means$home.pond, cage.means$test.elev, cage.means$source.elev), sd, na.rm = TRUE)
names(testpond.sd) <- c("test.pond", "source", "home.pond", "test.elev", "source.elev", "fit.SD", "surv.SD", "dev.SD", "mass.SD")
testpond.means     <- cbind(testpond.means, testpond.N[ , -c(1:7)], testpond.sd[ , -c(1:5)])
testpond.means$fit.SE  <- testpond.means$fit.SD / sqrt(testpond.means$N)
testpond.means$surv.SE <- testpond.means$surv.SD / sqrt(testpond.means$N)
testpond.means$dev.SE  <- testpond.means$dev.SD / sqrt(testpond.means$N)
testpond.means$mass.SE <- testpond.means$mass.SD / sqrt(testpond.means$N)
testpond.means <- testpond.means[ , c(1:5, 11, 6, 16, 7, 17, 8, 18, 9, 19)]
testpond.means <- testpond.means[order(testpond.means$test.pond, testpond.means$home.pond, testpond.means$source.elev), ]
#
# Next create treat.means. This is used for the interaction diagram.
source.means <- aggregate(cage.means[ , c("c.fit", "c.surv", "c.dev", "c.mass")], by = list(cage.means$test.pond, cage.means$source.pop, cage.means$test.elev, cage.means$source.elev), mean, na.rm = TRUE)
names(source.means)[1:4] <- c("test.pond", "source", "test.elev", "source.elev")
pond.means <- aggregate(source.means[ , c("c.fit", "c.surv", "c.dev", "c.mass")], by = list(source.means$test.pond, source.means$test.elev, source.means$source.elev), mean, na.rm = TRUE)
names(pond.means)[1:3] <- c("test.pond", "test.elev", "source.elev")
treat.means <- aggregate(pond.means[ , c("c.fit", "c.surv", "c.dev", "c.mass")], by = list(pond.means$test.elev, pond.means$source.elev), mean, na.rm = TRUE)
treat.sds   <- aggregate(pond.means[ , c("c.fit", "c.surv", "c.dev", "c.mass")], by = list(pond.means$test.elev, pond.means$source.elev), sd, na.rm = TRUE)
names(treat.sds)[3:6] <- c("fit.SD", "surv.SD", "dev.SD", "mass.SD")
treat.means <- merge(treat.means, treat.sds, by = c("Group.1", "Group.2"))
treat.N     <- aggregate(pond.means[ , c("c.fit", "c.surv", "c.dev", "c.mass")], by = list(pond.means$test.elev, pond.means$source.elev), length)
names(treat.N)[6] <- c("N")
treat.means <- merge(treat.means, treat.N[c(1,2,6)], by = c("Group.1", "Group.2"))
names(treat.means)[1:2] <- c("test.elev", "source.elev")
treat.means$fit.SE  <- treat.means$fit.SD / sqrt(treat.means$N)
treat.means$surv.SE <- treat.means$surv.SD / sqrt(treat.means$N)
treat.means$dev.SE  <- treat.means$dev.SD / sqrt(treat.means$N)
treat.means$mass.SE <- treat.means$mass.SD / sqrt(treat.means$N)
treat.means <- treat.means[ , c(1,2,11,3,7,12,4,8,13,5,9,14,6,10,15)]
treat.means$source.elev <- ordered(treat.means$source.elev, levels = c("low","high"))    # order factors as you want them in the figures
treat.means$test.elev   <- ordered(treat.means$test.elev, levels = c("low","high"))
#
#
##########################
#
# Fig. S4, adaptation to elevation.
# All the populations in each test site.
# 
#
testpond.means$x.seq  <- c(22, 23, 24, 21, 26, 27, 28, 29, 3, 4, 1, 2, 8, 9, 7, 6, 32, 33, 34, 31, 13, 14, 11, 12, 36, 37, 38, 39, 18, 19, 17, 16)
testpond.means        <- testpond.means[order(testpond.means$x.seq), ]
vector.of.pond.colors <- c(low.color, low.color, high.color, high.color, low.color, low.color, high.color, high.color, low.color, low.color, high.color, high.color, low.color, low.color, high.color, high.color, high.color, high.color, low.color, low.color, high.color, high.color, low.color, low.color, high.color, high.color, low.color, low.color, high.color, high.color, low.color, low.color)
home.ponds            <- c(6, 16, 21, 31)
home.elev             <- c(1, 2, 7, 11, 12, 17, 22, 26, 27, 32, 36, 37)
away.ponds            <- c(3, 4, 8, 9, 13, 14, 18, 19, 23, 24, 28, 29, 33, 34, 38, 39)
plot.margins          <- c(0.9, 0.72, 0.69, 0.51, 0.48, 0.3, 0.27, 0.09)
#
pdf(" .... your path here .... /Fig S4.pdf", width = 8, height = 13)
plot.new()
# Top panel: fitness
par(new = "TRUE", plt = c(0.22, 0.88, plot.margins[2], plot.margins[1]))
  columns  <- c("c.fit", "fit.SE")
  e        <- data.frame(test.pond = testpond.means$test.pond, source = testpond.means$source, x.seq = testpond.means$x.seq, y = testpond.means[ , columns[1]], SE = testpond.means[ , columns[2]], stringsAsFactors = FALSE)
  y.limits <- c( min(e$y - e$SE), max(e$y + e$SE) )
  plot(x = e$x.seq, y = e$y, xlim = c(0,40), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -2, ybottom = y.limits[1] - 50, xright = mean(e$x.seq), ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = mean(e$x.seq), y.limits[1] - 50, xright = 42, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  abline(v = 5); abline(v = 10); abline(v = 15); abline(v = 20); abline(v = 25); abline(v = 30); abline(v = 35)
  box(lwd = 1.3)
    # Text at the top indicating test ponds
  mtext( text = c("ellw", "feer", "munt", "siec", "bern", "bide", "flue", "muot"), side = 3, line = 0.26, cex = 0.95, at = c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5), adj = 0.5, col = c(rep(low.color, 4), rep(high.color, 4)))
  mtext( text = c("Test site"), side = 3, line = 1.6, cex = 1, adj = 0.5)
  par(xpd = TRUE); lines(x = c(1.1, 20), y = c(0.49, 0.49), col = low.color, lwd = 1.5); lines(x = c(20, 39), y = c(0.49, 0.49), col = high.color, lwd = 1.5); par(xpd = FALSE)
    # Vertical axis on the left
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext( text = "Expected fitness", side = 2, line = 2.2, las = 3, cex = 1.2)
  mtext(text = c("0.1", "0.2", "0.3", "0.4"), side = 2, line = 0.6, las = 1, at = c(0.1, 0.2, 0.3, 0.4), cex = 1)
    # Source ponds listed at the bottom
    # par(xpd = TRUE); text(x = c(e$x.seq), y = 0.02, labels = c(e$source), cex = 0.75, srt = 65, adj = c(1, NA), col = vector.of.pond.colors); par(xpd = FALSE)
    # Error bars first so that they are behind points:
  for (i in home.ponds) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = home.color.dark)
  for (i in home.elev) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = home.color)
  for (i in away.ponds) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = away.color.dark)
    # Points: population means
  points(x = e$x.seq[e$x.seq %in% home.ponds], y = e$y[e$x.seq %in% home.ponds], bg = home.color.dark, pch = 21, cex = 1.4)
  points(x = e$x.seq[e$x.seq %in% home.elev], y = e$y[e$x.seq %in% home.elev], bg = home.color, pch = 21, cex = 1.4)
  points(x = e$x.seq[e$x.seq %in% away.ponds], y = e$y[e$x.seq %in% away.ponds], bg = away.color, pch = 21, cex = 1.4)
  text(labels = c("away", "elevation"), x = 27.65, y = c(0.142, 0.119), col = away.color, cex = 0.8)
  text(labels = c("home", "elevation"), x = 27.5, y = c(0.295, 0.27), col = home.color, cex = 0.8)
  text(labels = c("home", "pond"), x = 21.7, y = c(0.256, 0.231), col = home.color.dark, cex = 0.8)
    # Label the upper margin:
  mtext(text = c("Fid. S4", "drawn by 'analysis of Judith transplant.R'"), side = 3, line = 4, cex = c(0.8, 0.6), at = c(3,25))
  text(labels = "A", x = 40, y = 0.97 * y.limits[2], cex = 2)
# Second panel: survival
par(new = "TRUE", plt = c(0.22, 0.88, plot.margins[4], plot.margins[3]))
  columns  <- c("c.surv", "surv.SE")
  e        <- data.frame(test.pond = testpond.means$test.pond, source = testpond.means$source, x.seq = testpond.means$x.seq, y = testpond.means[ , columns[1]], SE = testpond.means[ , columns[2]], stringsAsFactors = FALSE)
  y.limits <- c( min(e$y - e$SE), max(e$y + e$SE) )
  plot(x = e$x.seq, y = e$y, xlim = c(0,40), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -2, ybottom = y.limits[1] - 50, xright = mean(e$x.seq), ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = mean(e$x.seq), y.limits[1] - 50, xright = 42, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  abline(v = 5); abline(v = 10); abline(v = 15); abline(v = 20); abline(v = 25); abline(v = 30); abline(v = 35)
  box(lwd = 1.3)
    # Text at the top indicating test ponds
  mtext( text = c("ellw", "feer", "munt", "siec", "bern", "bide", "flue", "muot"), side = 3, line = 0.26, cex = 0.95, at = c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5), adj = 0.5, col = c(rep(low.color, 4), rep(high.color, 4)))
    # Vertical axis
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = "Proportion surviving", side = 2, line = 2.15, las = 3, cex = 1.2)
  mtext(text = c("0.3", "0.5", "0.7", "0.9"), side = 2, line = 0.6, las = 1, at = c(0.3, 0.5, 0.7, 0.9), cex = 1)
    # par(xpd = TRUE); text(x = c(e$x.seq), y = 0.205, labels = c(e$source), cex = 0.75, srt = 65, adj = c(1, NA), col = vector.of.pond.colors); par(xpd = FALSE)
    # Error bars first so that they are behind points:
  for (i in home.ponds) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = home.color.dark)
  for (i in home.elev) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = home.color)
  for (i in away.ponds) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = away.color.dark)
    # Population means
  points(x = e$x.seq[e$x.seq %in% home.ponds], y = e$y[e$x.seq %in% home.ponds], bg = home.color.dark, pch = 21, cex = 1.4)
  points(x = e$x.seq[e$x.seq %in% home.elev], y = e$y[e$x.seq %in% home.elev], bg = home.color, pch = 21, cex = 1.4)
  points(x = e$x.seq[e$x.seq %in% away.ponds], y = e$y[e$x.seq %in% away.ponds], bg = away.color, pch = 21, cex = 1.4)
  text(labels = "B", x = 40, y = 0.97 * y.limits[2], cex = 2)
# Third panel: development rate
par(new = "TRUE", plt = c(0.22, 0.88, plot.margins[6], plot.margins[5]))
  columns  <- c("c.dev", "dev.SE")
  e        <- data.frame(test.pond = testpond.means$test.pond, source = testpond.means$source, x.seq = testpond.means$x.seq, y = testpond.means[ , columns[1]], SE = testpond.means[ , columns[2]], stringsAsFactors = FALSE)
  y.limits <- c( min(e$y - e$SE), max(e$y + e$SE) )
  plot(x = e$x.seq, y = e$y, xlim = c(0,40), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -2, ybottom = y.limits[1] - 50, xright = mean(e$x.seq), ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = mean(e$x.seq), y.limits[1] - 50, xright = 42, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  abline(v = 5); abline(v = 10); abline(v = 15); abline(v = 20); abline(v = 25); abline(v = 30); abline(v = 35)
  box(lwd = 1.3)
    # Text at the top indicating test ponds
  mtext( text = c("ellw", "feer", "munt", "siec", "bern", "bide", "flue", "muot"), side = 3, line = 0.26, cex = 0.95, at = c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5), adj = 0.5, col = c(rep(low.color, 4), rep(high.color, 4)))
    # Vertical axis on the left
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext( text = "Development rate", side = 2, line = 3.4, las = 3, cex = 1.2)
  mtext( text = "(stages per day)", side = 2, line = 2.36, las = 3, cex = 1.2)
  mtext(text = c("0.26", "0.30", "0.34"), side = 2, line = 0.6, las = 1, at = c(0.26, 0.30, 0.34), cex = 1)
    # par(xpd = TRUE); text(x = c(e$x.seq), y = 0.24, labels = c(e$source), cex = 0.75, srt = 65, adj = c(1, NA), col = vector.of.pond.colors); par(xpd = FALSE)
    # Error bars first so that they are behind points:
  for (i in home.ponds) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = home.color.dark)
  for (i in home.elev) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = home.color)
  for (i in away.ponds) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = away.color.dark)
    # Population means
  points(x = e$x.seq[e$x.seq %in% home.ponds], y = e$y[e$x.seq %in% home.ponds], bg = home.color.dark, pch = 21, cex = 1.4)
  points(x = e$x.seq[e$x.seq %in% home.elev], y = e$y[e$x.seq %in% home.elev], bg = home.color, pch = 21, cex = 1.4)
  points(x = e$x.seq[e$x.seq %in% away.ponds], y = e$y[e$x.seq %in% away.ponds], bg = away.color, pch = 21, cex = 1.4)
  text(labels = "C", x = 40, y = 0.983 * y.limits[2], cex = 2)
# Bottom panel: mass
par(new = "TRUE", plt = c(0.22, 0.88, plot.margins[8], plot.margins[7]))
  columns  <- c("c.mass", "mass.SE")
  e        <- data.frame(test.pond = testpond.means$test.pond, source = testpond.means$source, x.seq = testpond.means$x.seq, y = testpond.means[ , columns[1]], SE = testpond.means[ , columns[2]], stringsAsFactors = FALSE)
  y.limits <- c( min(e$y - e$SE), max(e$y + e$SE) )
  plot(x = e$x.seq, y = e$y, xlim = c(0,40), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -2, ybottom = y.limits[1] - 50, xright = mean(e$x.seq), ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = mean(e$x.seq), y.limits[1] - 50, xright = 42, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  abline(v = 5); abline(v = 10); abline(v = 15); abline(v = 20); abline(v = 25); abline(v = 30); abline(v = 35)
  box(lwd = 1.3)
    # Text at the top indicating test ponds
  mtext( text = c("ellw", "feer", "munt", "siec", "bern", "bide", "flue", "muot"), side = 3, line = 0.26, cex = 0.95, at = c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5), adj = 0.5, col = c(rep(low.color, 4), rep(high.color, 4)))
    # Vertical axis on the left
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext( text = "Mass at", side = 2, line = 3.4, las = 3, cex = 1.2)
  mtext( text = "metamorphosis (mg)", side = 2, line = 2.36, las = 3, cex = 1.2)
  mtext(text = c("350", "400", "450", "500", "550"), side = 2, line = 0.6, las = 1, at = c(350,400,450,500,550), cex = 1)
  par(xpd = TRUE)
  text(x = c(e$x.seq), y = 333, labels = c(e$source), cex = 0.75, srt = 65, adj = c(1, NA), col = vector.of.pond.colors)
  par(xpd = FALSE)
    # Error bars first so that they are behind points:
  for (i in home.ponds) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = home.color.dark)
  for (i in home.elev) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = home.color)
  for (i in away.ponds) lines(x = rep(i, 2), y = c(e$y[e$x.seq == i] - e$SE[e$x.seq == i], e$y[e$x.seq == i] + e$SE[e$x.seq == i]), col = away.color.dark)
    # Population means
  points(x = e$x.seq[e$x.seq %in% home.ponds], y = e$y[e$x.seq %in% home.ponds], bg = home.color.dark, pch = 21, cex = 1.4)
  points(x = e$x.seq[e$x.seq %in% home.elev], y = e$y[e$x.seq %in% home.elev], bg = home.color, pch = 21, cex = 1.4)
  points(x = e$x.seq[e$x.seq %in% away.ponds], y = e$y[e$x.seq %in% away.ponds], bg = away.color, pch = 21, cex = 1.4)
    # Label x axis:
  mtext( text = "Source population", side = 1, line = 1.9, las = 1, cex = 1.2)
  text(labels = "D", x = 40, y = 0.98 * y.limits[2], cex = 2)
dev.off()
#
#
###########################################################
#
# Fig. 2 in the MS.
# Four simple interaction diagrams.
#
# Combine the treatment means and fitted estimates into the same data frame.
cage.means              <- aggregate(d[,c("fitness", "survival", "DevRate42", "mass.42")], by = list(test.pond = d$test.pond, cage = d$cage, source.pop = d$source.pop, test.elev = d$test.elev, home.elev = d$home.elev), mean, na.rm = TRUE)
cage.means              <- rename(cage.means, c("fitness", "survival", "DevRate42", "mass.42"), c("fit", "surv", "dev", "mass"))
test.source.pond.means  <- aggregate(cage.means[,c("fit", "surv", "dev", "mass")], by = list(test.pond = cage.means$test.pond, source.pop = cage.means$source.pop, test.elev = cage.means$test.elev, home.elev = cage.means$home.elev), mean, na.rm = TRUE)
test.pond.means         <- aggregate(test.source.pond.means[,c("fit", "surv", "dev", "mass")], by = list(test.pond = test.source.pond.means$test.pond, test.elev = test.source.pond.means$test.elev, home.elev = test.source.pond.means$home.elev), mean, na.rm = TRUE)
treat.means             <- aggregate(test.pond.means[,c("fit", "surv", "dev", "mass")], by = list(test.elev = test.pond.means$test.elev, home.elev = test.pond.means$home.elev), mean, na.rm = TRUE)
treat.sd                <- aggregate(test.pond.means[,c("fit", "surv", "dev", "mass")], by = list(test.elev = test.pond.means$test.elev, home.elev = test.pond.means$home.elev), sd, na.rm = TRUE)
treat.sd                <- rename(treat.sd, c("fit", "surv", "dev", "mass"), c("fit.sd", "surv.sd", "dev.sd", "mass.sd"))
treat.sd$fit.se         <- treat.sd$fit.sd / 2; treat.sd$surv.se <- treat.sd$surv.sd / 2; treat.sd$dev.se <- treat.sd$dev.sd / 2; treat.sd$mass.se <- treat.sd$mass.sd / 2
treat.sd$fit.sd         <- treat.sd$surv.sd <- treat.sd$dev.sd <- treat.sd$mass.sd <- NULL
treat.means             <- merge(treat.means, treat.sd, by = c("test.elev", "home.elev"))
treat.means$source.elev <- c("low", "high", "high", "low")
point.offset            <- 0.12 / 2
treat.means$x           <- c(7.5-point.offset, 7.5+point.offset, 2.5+point.offset, 2.5-point.offset)
treat.means$color       <- c(low.color, high.color, high.color, low.color)
treat.means             <- treat.means[order(treat.means$x), ]
plot.margins      <- c(0.92, 0.74, 0.705, 0.525, 0.49, 0.31, 0.275, 0.095)
error.bar.width   <- 1.5
symbol.size       <- 1.9
symbol.line.width <- 1.2
axis.label.size   <- 1.1
#
pdf(" .... your path here .... /Fig 2.pdf", useDingbats = FALSE, width = 5, height = 9)
plot.new()
# panel A: fitness
par(new = "TRUE", plt = c(0.25, 0.85, plot.margins[2], plot.margins[1]))
  dat       <- treat.means[ , c("x", "color", "fit", "fit.se")]
  dat       <- rename(dat, c("fit", "fit.se"), c("y", "y.se"))
  dat$lower <- dat[ , "y"] - dat[ , "y.se"]
  dat$upper <- dat[ , "y"] + dat[ , "y.se"]
  y.limits  <- c(min(dat$lower), max(dat$upper))
  plot(x = c(2.5,7.5), y = c(0,1), xlim = c(1,9), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -1, ybottom = 0-10, xright = 5, ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = 5, ybottom = -10, xright = 11, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  abline(v = 5, lwd = 1.1)
  box(lwd = 1.3)
    # lines connecting points, error bars, and points
  lines(x = dat$x[c(2,4)], y = dat$y[c(2,4)], col = high.color, lwd = 2.5, lty = 2)
  lines(x = c(dat$x[2], dat$x[2]), y = c(dat$lower[2], dat$upper[2]), lwd = error.bar.width, col = high.color)
  lines(x = c(dat$x[4], dat$x[4]), y = c(dat$lower[4], dat$upper[4]), lwd = error.bar.width, col = high.color)
  points(x = dat$x[c(2,4)], y = dat$y[c(2,4)], lwd = symbol.line.width, pch = 21, cex = symbol.size, col = high.color, bg = high.color)
  lines(x = dat$x[c(1,3)], y = dat$y[c(1,3)], col = low.color, lwd = 2.5, lty = 1)
  lines(x = c(dat$x[1], dat$x[1]), y = c(dat$lower[1], dat$upper[1]), lwd = error.bar.width, col = low.color)
  lines(x = c(dat$x[3], dat$x[3]), y = c(dat$lower[3], dat$upper[3]), lwd = error.bar.width, col = low.color)
  points(x = dat$x[c(1,3)], y = dat$y[c(1,3)], lwd = symbol.line.width, pch = 21, cex = symbol.size, col = low.color, bg = low.color)
    # Label the axes
  axis(side = 1, tck = -0.015, at = c(2.5, 7.5), labels = FALSE)
  mtext(c("Low", "High"), side = 1, line = 0.2, outer = FALSE, cex = 1.0, at = c(2.5,7.5), col = c(low.color, high.color))
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.1", "0.2", "0.3"), side = 2, line = 0.4, las = 1, at = c(0.1,0.2,0.3), cex = 0.95)
  mtext("Expected fitness", side = 2, line = 2.1, outer = FALSE, cex = axis.label.size)
  text(x = c(7.9, 6.5), y = c(0.15, 0.25), label = c("low source", "high source"), cex = 0.9, col = c(low.color,high.color))
  mtext(c("Fig. 2", "drawn by 'analysis of Judith transplant.R'"), side = 3, line = 2.3, outer = FALSE, cex = c(0.8, 0.6), at = c(2,7))
  text(labels = "A", x = 8.9, y = 0.97 * y.limits[2], cex = 1.3)
# panel B: survival
par(new = "TRUE", plt = c(0.25, 0.85, plot.margins[4], plot.margins[3]))
  dat       <- treat.means[ , c("x", "color", "surv", "surv.se")]
  dat       <- rename(dat, c("surv", "surv.se"), c("y", "y.se"))
  dat$lower <- dat[ , "y"] - dat[ , "y.se"]
  dat$upper <- dat[ , "y"] + dat[ , "y.se"]
  y.limits  <- c(min(dat$lower), max(dat$upper))
  plot(x = c(2.5,7.5), y = c(0,1), xlim = c(1,9), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -1, ybottom = 0-10, xright = 5, ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = 5, ybottom = -10, xright = 11, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  abline(v = 5, lwd = 1.1)
  box(lwd = 1.3)
    # lines connecting points, error bars, and points
  lines(x = dat$x[c(2,4)], y = dat$y[c(2,4)], col = high.color, lwd = 2.5, lty = 2)
  lines(x = c(dat$x[2], dat$x[2]), y = c(dat$lower[2], dat$upper[2]), lwd = error.bar.width, col = high.color)
  lines(x = c(dat$x[4], dat$x[4]), y = c(dat$lower[4], dat$upper[4]), lwd = error.bar.width, col = high.color)
  points(x = dat$x[c(2,4)], y = dat$y[c(2,4)], lwd = symbol.line.width, pch = 21, cex = symbol.size, col = high.color, bg = high.color)
  lines(x = dat$x[c(1,3)], y = dat$y[c(1,3)], col = low.color, lwd = 2.5, lty = 1)
  lines(x = c(dat$x[1], dat$x[1]), y = c(dat$lower[1], dat$upper[1]), lwd = error.bar.width, col = low.color)
  lines(x = c(dat$x[3], dat$x[3]), y = c(dat$lower[3], dat$upper[3]), lwd = error.bar.width, col = low.color)
  points(x = dat$x[c(1,3)], y = dat$y[c(1,3)], lwd = symbol.line.width, pch = 21, cex = symbol.size, col = low.color, bg = low.color)
    # Label the axes
  axis(side = 1, tck = -0.015, at = c(2.5, 7.5), labels = FALSE)
  mtext(c("Low", "High"), side = 1, line = 0.2, outer = FALSE, cex = 1.0, at = c(2.5,7.5), col = c(low.color, high.color))
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.4", "0.5", "0.6", "0.7", "0.8"), side = 2, line = 0.4, las = 1, at = c(0.4,0.5,0.6,0.7,0.8), cex = 0.95)
  mtext("Survival to", side = 2, line = 3.0, outer = FALSE, cex = axis.label.size)
  mtext("metamorphosis", side = 2, line = 2.0, outer = FALSE, cex = axis.label.size)
  text(labels = "B", x = 8.9, y = 0.975 * y.limits[2], cex = 1.3)
    #
# panel C: development rate
par(new = "TRUE", plt = c(0.25, 0.85, plot.margins[6], plot.margins[5]))
  dat       <- treat.means[ , c("x", "color", "dev", "dev.se")]
  dat       <- rename(dat, c("dev", "dev.se"), c("y", "y.se"))
  dat$lower <- dat[ , "y"] - dat[ , "y.se"]
  dat$upper <- dat[ , "y"] + dat[ , "y.se"]
  y.limits <- c(min(dat$lower), max(dat$upper))
  plot(x = c(2.5,7.5), y = c(0,1), xlim = c(1,9), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -1, ybottom = 0-10, xright = 5, ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = 5, ybottom = -10, xright = 11, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  abline(v = 5, lwd = 1.1)
  box(lwd = 1.3)
    # lines connecting points, error bars, and points
  lines(x = dat$x[c(2,4)], y = dat$y[c(2,4)], col = high.color, lwd = 2.5, lty = 2)
  lines(x = c(dat$x[2], dat$x[2]), y = c(dat$lower[2], dat$upper[2]), lwd = error.bar.width, col = high.color)
  lines(x = c(dat$x[4], dat$x[4]), y = c(dat$lower[4], dat$upper[4]), lwd = error.bar.width, col = high.color)
  points(x = dat$x[c(2,4)], y = dat$y[c(2,4)], lwd = symbol.line.width, pch = 21, cex = symbol.size, col = high.color, bg = high.color)
  lines(x = dat$x[c(1,3)], y = dat$y[c(1,3)], col = low.color, lwd = 2.5, lty = 1)
  lines(x = c(dat$x[1], dat$x[1]), y = c(dat$lower[1], dat$upper[1]), lwd = error.bar.width, col = low.color)
  lines(x = c(dat$x[3], dat$x[3]), y = c(dat$lower[3], dat$upper[3]), lwd = error.bar.width, col = low.color)
  points(x = dat$x[c(1,3)], y = dat$y[c(1,3)], lwd = symbol.line.width, pch = 21, cex = symbol.size, col = low.color, bg = low.color)
    # Label the axes
  axis(side = 1, tck = -0.015, at = c(2.5, 7.5), labels = FALSE)
  mtext(c("Low", "High"), side = 1, line = 0.2, outer = FALSE, cex = 1.0, at = c(2.5,7.5), col = c(low.color, high.color))
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("0.26", "0.30", "0.34"), side = 2, line = 0.4, las = 1, at = c(0.26, 0.3, 0.34), cex = 0.95)
  mtext("Development rate", side = 2, line = 3.0, outer = FALSE, cex = axis.label.size)
  mtext("(stages per day)", side = 2, line = 2.0, outer = FALSE, cex = axis.label.size)
  text(labels = "C", x = 8.9, y = 0.988 * y.limits[2], cex = 1.3)
    #
# panel D: mass at metamorphosis
par(new = "TRUE", plt = c(0.25, 0.85, plot.margins[8], plot.margins[7]))
  dat       <- treat.means[ , c("x", "color", "mass", "mass.se")]
  dat       <- rename(dat, c("mass", "mass.se"), c("y", "y.se"))
  dat$lower <- dat[ , "y"] - dat[ , "y.se"]
  dat$upper <- dat[ , "y"] + dat[ , "y.se"]
  y.limits  <- c(min(dat$lower), max(dat$upper))
  plot(x = c(2.5,7.5), y = c(0,1), xlim = c(1,9), ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = -1, ybottom = 0-10, xright = 5, ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = 5, ybottom = -10, xright = 11, ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  abline(v = 5, lwd = 1.1)
  box(lwd = 1.3)
    # lines connecting points, error bars, and points
  lines(x = dat$x[c(2,4)], y = dat$y[c(2,4)], col = high.color, lwd = 2.5, lty = 2)
  lines(x = c(dat$x[2], dat$x[2]), y = c(dat$lower[2], dat$upper[2]), lwd = error.bar.width, col = high.color)
  lines(x = c(dat$x[4], dat$x[4]), y = c(dat$lower[4], dat$upper[4]), lwd = error.bar.width, col = high.color)
  points(x = dat$x[c(2,4)], y = dat$y[c(2,4)], lwd = symbol.line.width, pch = 21, cex = symbol.size, col = high.color, bg = high.color)
  lines(x = dat$x[c(1,3)], y = dat$y[c(1,3)], col = low.color, lwd = 2.5, lty = 1)
  lines(x = c(dat$x[1], dat$x[1]), y = c(dat$lower[1], dat$upper[1]), lwd = error.bar.width, col = low.color)
  lines(x = c(dat$x[3], dat$x[3]), y = c(dat$lower[3], dat$upper[3]), lwd = error.bar.width, col = low.color)
  points(x = dat$x[c(1,3)], y = dat$y[c(1,3)], lwd = symbol.line.width, pch = 21, cex = symbol.size, col = low.color, bg = low.color)
    # Label the axes
  axis(side = 1, tck = -0.015, at = c(2.5, 7.5), labels = FALSE)
  mtext(c("Low", "High"), side = 1, line = 0.2, outer = FALSE, cex = 1.0, at = c(2.5,7.5), col = c(low.color, high.color))
  mtext("Elevation of test site", side = 1, line = 1.4, outer = FALSE, cex = axis.label.size)
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(text = c("400", "500", "600"), side = 2, line = 0.4, las = 1, at = c(400,500,600), cex = 0.95)
  mtext("Mass at", side = 2, line = 3.0, outer = FALSE, cex = axis.label.size)
  mtext("metamorphosis (mg)", side = 2, line = 2.0, outer = FALSE, cex = axis.label.size)
  text(labels = "D", x = 8.9, y = 0.984 * y.limits[2], cex = 1.3)
dev.off()
#
#
############################
#
# Relative importance of unique and parallel local adaptation.
# Draw Fig. 4.
# This section looks at the gradient from home pond --> home habitat --> away habitat
#
cage.means        <- aggregate(d[,c("fitness", "survival", "DevRate42", "mass.42")], by = list(d$test.pond, d$cage, d$source.pop, d$test.elev, d$home.elev), mean, na.rm = TRUE)
names(cage.means) <- c("test.pond", "cage", "source.pop", "test.elev", "home.elev", "fit", "surv", "dev", "mass")
#
# Set up three datasets for the three kinds of tests (home pond, home elev, away elev).
#
home.pond.cages <- subset(cage.means, test.pond == source.pop)
home.pond.means <- aggregate(home.pond.cages[, c("fit")], by = list(home.pond.cages$test.pond, home.pond.cages$test.elev), mean, na.rm = TRUE)
home.pond.sds   <- aggregate(home.pond.cages[, c("fit")], by = list(home.pond.cages$test.pond, home.pond.cages$test.elev), sd, na.rm = TRUE)
home.pond.sds$x <- home.pond.sds$x / sqrt(8)
names(home.pond.means) <- c("test.pond", "test.elev", "fit")
names(home.pond.sds) <- c("test.pond", "test.elev", "fit.se")
home.pond       <- merge(home.pond.means, home.pond.sds)
home.pond$lower <- home.pond$fit - home.pond$fit.se
home.pond$upper <- home.pond$fit + home.pond$fit.se
home.pond$x.val <- c(3,1,4,2)
home.pond       <- home.pond[order(home.pond$x.val),]
home.pond
#
home.elev.cages <- subset(cage.means, home.elev == "home_elev" & test.pond != source.pop)
home.elev.means <- aggregate(home.elev.cages[, c("fit")], by = list(home.elev.cages$test.pond, home.elev.cages$source.pop, home.elev.cages$test.elev), mean, na.rm = TRUE)
home.elev.sds   <- aggregate(home.elev.cages[, c("fit")], by = list(home.elev.cages$test.pond, home.elev.cages$source.pop, home.elev.cages$test.elev), sd, na.rm = TRUE)
home.elev.sds$x <- home.elev.sds$x / sqrt(8)
names(home.elev.means) <- c("test.pond", "source.pop", "test.elev", "fit")
names(home.elev.sds)   <- c("test.pond", "source.pop", "test.elev", "fit.se")
home.elev       <- merge(home.elev.means, home.elev.sds)
home.elev$lower <- home.elev$fit - home.elev$fit.se
home.elev$upper <- home.elev$fit + home.elev$fit.se
home.elev$x.val <- c(4,3,4.1,1,2,2.1,3.1,1.1,2.2,3.2,4.2,1.2)
home.elev       <- home.elev[order(home.elev$x.val),]
home.elev
#
away.elev.cages <- subset(cage.means, home.elev == "away_elev")
away.elev.means <- aggregate(away.elev.cages[, c("fit")], by = list(away.elev.cages$test.pond, away.elev.cages$source.pop, away.elev.cages$test.elev), mean, na.rm = TRUE)
away.elev.sds   <- aggregate(away.elev.cages[, c("fit")], by = list(away.elev.cages$test.pond, away.elev.cages$source.pop, away.elev.cages$test.elev), sd, na.rm = TRUE)
away.elev.sds$x <- away.elev.sds$x / sqrt(8)
names(away.elev.means) <- c("test.pond", "source.pop", "test.elev", "fit")
names(away.elev.sds)   <- c("test.pond", "source.pop", "test.elev", "fit.se")
away.elev       <- merge(away.elev.means, away.elev.sds)
away.elev$x.val <- c(1,2,1.1,2.1,3,4,3.1,4.1,1.2,2.2,3.2,4.2,1.3,2.3,3.3,4.3)
away.elev$lower <- away.elev$fit - away.elev$fit.se
away.elev$upper <- away.elev$fit + away.elev$fit.se
away.elev       <- away.elev[order(away.elev$x.val),]
away.elev
#
home.pond$x.val   <- c(2,12,23,33)
home.el         <- home.elev[c(5,3,10,8),]
home.el$x.val   <- c(4,14,25,35)
away.el         <- away.elev[c(10,14,12,16,1,5,3,7),]
away.el$x.val   <- c(6,7,16,17,27,28,37,38)
#
x.limits          <- c(1, 39)
y.limits          <- c(0, 0.38)
error.bar.width   <- 1.1
symbol.line.width <- 1
symbol.size       <- 2
boundaries        <- c(0,9.5,20,30.5,40)
low.away.color    <- "#FFCD00"
high.away.color   <- "#378EAE"
pdf(" .... your path here .... /Fig 4.pdf", width = 9, height = 5.5, useDingbats = FALSE)
plot.new()
par(new = "TRUE", plt = c(0.25, 0.85, 0.2, 0.8))
  plot(x = c(0,1), y = c(0,1), xlim = x.limits, ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "")
  rect(xleft = 0-10, ybottom = 0-10, xright = boundaries[3], ytop = 1.1*y.limits[2], col = low.color.pale, border = NA)
  rect(xleft = boundaries[3], ybottom = -10, xright = 10+x.limits[2], ytop = 1.1*y.limits[2], col = high.color.pale, border = NA)
  abline(v = boundaries[2:4], lwd = c(1.1, 2, 1.1))
  abline(h = 0, lwd = 1.1, lty = 2)
  box(lwd = 1.3)
    # Text at the top indicating test ponds
  x.locations <- c(5, 14.75, 25.25, 35)
  mtext( text = c("feer", "siec", "bern", "flue"), side = 3, line = 0.26, cex = 0.95, at = x.locations, adj = 0.5, col = c(rep(low.color, 2), rep(high.color, 2)))
  mtext( text = c("Test site"), side = 3, line = 1.6, cex = 1, adj = 0.5)
  par(xpd = TRUE); lines(x = c(4, boundaries[3]), y = c(0.425, 0.425), col = low.color, lwd = 1.5); lines(x = c(boundaries[3], 36), y = c(0.425, 0.425), col = high.color, lwd = 1.5); par(xpd = FALSE)
    # horizontal axis
  mtext(c("Source population"), side = 1, line = 1.8, cex = 1.3 )
  par(xpd = TRUE)
  text(x = home.pond$x.val, y = -0.02, labels = home.pond$test.pond, cex = 0.75, srt = 65, adj = c(1, NA), col = c(low.color,low.color,high.color,high.color))
  text(x = home.el$x.val, y = -0.02, labels = home.el$source.pop, cex = 0.75, srt = 65, adj = c(1, NA), col = c(rep(low.color,2), rep(high.color,2)))
  text(x = away.el$x.val, y = -0.02, labels = away.el$source.pop, cex = 0.75, srt = 65, adj = c(1, NA), col = c(rep(high.color,4), rep(low.color,4)))
  par(xpd = FALSE)
    # vertical axis
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(c("0.0", "0.1", "0.2", "0.3"), side = 2, las = 1, line = 0.4, cex = 1.0, at = c(0,0.1,0.2,0.3) )
  mtext("Expected fitness", side = 2, line = 2, cex = 1.3 )
    # error bars
  for(i in 1:2)  { lines(x = c(home.pond$x.val[i], home.pond$x.val[i]), y = c(home.pond$lower[i], home.pond$upper[i]), lwd = error.bar.width, col = low.color) }
  for(i in 3:4)  { lines(x = c(home.pond$x.val[i], home.pond$x.val[i]), y = c(home.pond$lower[i], home.pond$upper[i]), lwd = error.bar.width, col = high.color) }
  for(i in 1:2)  { lines(x = c(home.el$x.val[i], home.el$x.val[i]), y = c(home.el$lower[i], home.el$upper[i]), lwd = error.bar.width, col = low.color) }
  for(i in 3:4)  { lines(x = c(home.el$x.val[i], home.el$x.val[i]), y = c(home.el$lower[i], home.el$upper[i]), lwd = error.bar.width, col = high.color) }
  for(i in 1:4)  { lines(x = c(away.el$x.val[i], away.el$x.val[i]), y = c(away.el$lower[i], away.el$upper[i]), lwd = error.bar.width, col = high.color) }
  for(i in 5:8) { lines(x = c(away.el$x.val[i], away.el$x.val[i]), y = c(away.el$lower[i], away.el$upper[i]), lwd = error.bar.width, col = low.color) }
    # points
  points(x = home.pond$x.val[1:2], y = home.pond$fit[1:2], lwd = symbol.line.width, pch = 21, cex = symbol.size, col = low.color, bg = low.color)
  points(x = home.el$x.val[1:2], y = home.el$fit[1:2], lwd = symbol.line.width, pch = 21, cex = 0.7*(symbol.size), col = "black", bg = low.away.color)
  points(x = away.el$x.val[1:4], y = away.el$fit[1:4], lwd = symbol.line.width, pch = 21, cex = 0.7*(symbol.size), col = "black", bg = high.away.color)
  points(x = home.pond$x.val[3:4], y = home.pond$fit[3:4], lwd = symbol.line.width, pch = 21, cex = symbol.size, col = high.color, bg = high.color)
  points(x = home.el$x.val[3:4], y = home.el$fit[3:4], lwd = symbol.line.width, pch = 21, cex = 0.7*(symbol.size), col = "black", bg = high.away.color)
  points(x = away.el$x.val[5:8], y = away.el$fit[5:8], lwd = symbol.line.width, pch = 21, cex = 0.7*(symbol.size), col = "black", bg = low.away.color)
    # labels in the middle of the figure
  text(x = c(4, 1.9, 1.9, 7.4, 7.4), y = c(0.335, 0.21, 0.195, 0.13, 0.115), labels = c("local", "same", "elevation", "away", "elevation"), cex = 0.75)
  text(x = c(32, 37.5, 37.5, 36, 36), y = c(0.265, 0.225, 0.21, 0.09, 0.075), labels = c("local", "same", "elevation", "away", "elevation"), cex = 0.75)
    # label the figure at the top
  mtext(text = c("Fig. 4", "drawn by 'analysis of Judith transplant.R'"), side = 3, line = 4, cex = c(0.8, 0.6), at = c(3,35))
dev.off()
#
#
################################
#
# Analyze the review of transplant experiments in plants and animals.
#
te <- read.table("transplant experiments across elevation.txt", header = TRUE, stringsAsFactors = FALSE)
te <- te[order(te$taxon, te$species),]
#
# Figure out how this sample compares with Halbritter et al's (2018, JEB)
proc.means(te, cats = c("taxon", "halbritter"), vars = "species", stats = "N")
length(unique(te$author[te$taxon == "anim"]))                             # Nr. distinct papers in the animal dataset
length(unique(te$species[te$taxon == "anim"]))                            # Nr. species in the animal dataset
length(unique(te$author[te$taxon == "plant" & te$halbritter == "yes"]))   # Nr. distinct papers in the plant dataset, in and out of the Halbritter sa,mple.
length(unique(te$author[te$taxon == "plant" & te$halbritter == "no"]))
length(unique(te$species[te$taxon == "plant" & te$halbritter == "yes"]))  # Nr. species in the plant dataset, in and out of the Halbritter sample.
length(unique(te$species[te$taxon == "plant" & te$halbritter == "no"]))
length(unique(te$author[te$taxon == "plant"]))  # Nr. papers in the plant dataset
length(unique(te$species[te$taxon == "plant"]))  # Nr. species in the plant dataset
combo <- paste(te$halbritter[te$taxon == "plant"], te$author[te$taxon == "plant"], sep = "_")
combo <- unique(combo)
combo[order(combo)]
combo <- paste(te$halbritter[te$taxon == "plant"], te$species[te$taxon == "plant"], sep = "_")
combo <- unique(combo)
combo[order(combo)]
#    Halbritter et al. (2018): 14 publications on 18 species
#    Not in Halbritter: 30 publications on 47 plant species (one of which is the same as in Halbritter et al.'s sample)
#
te$G.pos.fit <- te$G.pos.fit / te$G.N.fit      # these are the garden-specific tests (Local vs Foreign)
te$G.neg.fit <- te$G.neg.fit / te$G.N.fit
te$G.pos.surv <- te$G.pos.surv / te$G.N.surv
te$G.neg.surv <- te$G.neg.surv / te$G.N.surv
te$G.pos.size <- te$G.pos.size/ te$G.N.size
te$G.neg.size <- te$G.neg.size / te$G.N.size
te$G.pos.repr <- te$G.pos.repr / te$G.N.repr
te$G.neg.repr <- te$G.neg.repr / te$G.N.repr
te$P.pos.fit <- te$P.pos.fit / te$P.N.fit      # these are the population-specific tests (Home versus Away)
te$P.neg.fit <- te$P.neg.fit / te$P.N.fit
te$P.pos.surv <- te$P.pos.surv / te$P.N.surv
te$P.neg.surv <- te$P.neg.surv / te$P.N.surv
te$P.pos.size <- te$P.pos.size/ te$P.N.size
te$P.neg.size <- te$P.neg.size / te$P.N.size
te$P.pos.repr <- te$P.pos.repr / te$P.N.repr
te$P.neg.repr <- te$P.neg.repr / te$P.N.repr
te1           <- te
te1[ ,c(6,7,9,10,12,13,15,16,18,19,21,22,24,25,27,28)] <- apply(te[,c(6,7,9,10,12,13,15,16,18,19,21,22,24,25,27,28)], 2, format.p.val, 2)
# Print out Appendix Table S4.
te1[,c("species", "G.N.fit", "G.pos.fit", "G.neg.fit", "G.N.surv", "G.pos.surv", "G.neg.surv", "G.N.size", "G.pos.size", "G.neg.size", "G.N.repr", "G.pos.repr", "G.neg.repr", "P.N.fit", "P.pos.fit", "P.neg.fit", "P.N.surv", "P.pos.surv", "P.neg.surv", "P.N.size", "P.pos.size", "P.neg.size", "P.N.repr", "P.pos.repr", "P.neg.repr", "author")]
te1[,c("species", "G.pos.fit", "G.neg.fit", "G.pos.surv", "G.neg.surv", "G.pos.size", "G.neg.size", "G.pos.repr", "G.neg.repr", "P.pos.fit", "P.neg.fit", "P.pos.surv", "P.neg.surv", "P.pos.size", "P.neg.size", "P.pos.repr", "P.neg.repr", "author")]
pm   <- proc.means(te, cats = "taxon", vars = names(te)[c(6,7,9,10,12,13,15,16,18,19,21,22,24,25,27,28)], stats = c("N", "mean", "SD"))
pm.a <- subset(pm, taxon == "anim")
pm.a <- matrix(pm.a[,-1], ncol = 3, byrow = TRUE)
pm.p <- subset(pm, taxon == "plant")
pm.p <- matrix(pm.p[,-1], ncol = 3, byrow = TRUE)
pm   <- as.data.frame(cbind(pm.p, pm.a))
names(pm) <- c("plant.N", "plant.mean", "plant.SD", "animal.N", "animal.mean", "animal.SD")
pm[,c(2,3,5,6)] <- apply(pm[,c(2,3,5,6)], 2, format.p.val, 3)
pm$variable <- c("G.pos.fit", "G.neg.fit", "G.pos.surv", "G.neg.surv", "G.pos.size", "G.neg.size", "G.pos.repr", "G.neg.repr", "P.pos.fit", "P.neg.fit", "P.pos.surv", "P.neg.surv", "P.pos.size", "P.neg.size", "P.pos.repr", "P.neg.repr")
summry <- pm[,c(7,1:6)]
summry
#
#
# Draw Appendix Figure S9
# Illustrates results of the numerous plant reciprocal transplant studies across elevation.
#
which.cols       <- c(6,7,9,10,12,13,15,16,18,19,21,22,24,25,27,28)
plant.pos.col    <- low.color
plant.neg.col    <- high.color
anim.pos.col     <- "#FF4900"
anim.neg.col     <- "#AB02C5"
lwidth           <- 0
border.col       <- NULL
tick.label.size  <- 0.9
axis.label.size  <- 1.2
panel.label.size <- 1.1
buffers          <- 0.0
top.margin       <- 0.1
bottom.margin    <- 0.1
vert.open.space  <- 1.09
arrow.width      <- 1.6
arrow.length     <- 0.23
arrow.height     <- 5.5
panel            <- (1 - (top.margin+bottom.margin) - (3*buffers)) / 4
fig.marg         <- c(bottom.margin, bottom.margin+panel, bottom.margin+panel+buffers, bottom.margin+2*panel+buffers, bottom.margin+2*panel+2*buffers, bottom.margin+3*panel+2*buffers, bottom.margin+3*panel+3*buffers, bottom.margin+4*panel+3*buffers)
horiz.marg       <- c(0.15, 0.5, 0.55, 0.9)
xlimits          <- c(-1,1)
pdf(" .... your path here .... /Fig S9.pdf", width = 8, height = 11.0, useDingbats = FALSE)
plot.new()
# Local/foreign: "fitness"
par(new = "TRUE", plt = c(horiz.marg[1], horiz.marg[2], fig.marg[7], fig.marg[8]))
  h1 <- hist(te[te$taxon == "plant",which.cols[1]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2 <- hist(te[te$taxon == "plant",which.cols[2]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h3 <- hist(te[te$taxon == "anim",which.cols[1]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h4 <- hist(te[te$taxon == "anim",which.cols[2]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2$breaks <- -1 * h2$breaks
  ylimits <- c(0, vert.open.space * max(c(h1$counts+h3$counts, h2$counts+h4$counts)))
  plot(x = c(0,1), y = c(0,1), xlim = xlimits, ylim = ylimits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "", xaxs = "i", yaxs = "i")
  for (i in 1:length(h1$counts)) {
    rect(xleft = h1$breaks[i], ybottom = 0, xright = h1$breaks[i+1], ytop = h1$counts[i], lwd = lwidth, border = border.col, col = plant.pos.col)
    rect(xleft = h1$breaks[i], ybottom = h1$counts[i], xright = h1$breaks[i+1], ytop = h1$counts[i]+h3$counts[i], lwd = lwidth, border = border.col, col = anim.pos.col)
    rect(xleft = h2$breaks[i+1], ybottom = 0, xright = h2$breaks[i], ytop = h2$counts[i], lwd = lwidth, border = border.col, col = plant.neg.col)
    rect(xleft = h2$breaks[i+1], ybottom = h2$counts[i], xright = h2$breaks[i], ytop = h2$counts[i]+h4$counts[i], lwd = lwidth, border = border.col, col = anim.neg.col)
    }
  box(lwd = 1.3)
  # Add arrows indicating adaptation and maladaptation.
  arrows(-0.13, arrow.height, x1 = -0.13 - arrow.length, y1 = arrow.height, lwd = arrow.width, col = high.color, length = 0.1)
  text(c("maladapta-", "tion to elev."), x = -0.09 - arrow.length, y = c(arrow.height+0.32,arrow.height-0.32), pos = 2, cex = 0.8, col = high.color)
  arrows(0.10, arrow.height, x1 = 0.10 + arrow.length, y1 = arrow.height, lwd = arrow.width, col = low.color, length = 0.1)
  text(c("adaptation", "to elevation"), x = 0.1 + arrow.length, y = c(arrow.height+0.36, arrow.height-0.36), pos = 4, cex = 0.8, col = low.color)
  axis(side = 1, tck = -0.015, labels = FALSE, at = c(-1,-0.5,0,0.5,1))
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(c("2", "4", "6", "8", "10"), side = 2, las = 1, line = 0.4, cex = tick.label.size, at = c(2,4,6,8,10) )
  text("A. Fitness", x = -1.03, y = 0.93 * ylimits[2], cex = panel.label.size, pos = 4)
  mtext(c("Local-foreign criterion", "within test sites"), side = 3, line = c(1.6,0.45), cex = axis.label.size+0.1 )
  mtext("Fig. S9", side = 3, line = 4, cex = 0.6)
# Local/foreign: survival
par(new = "TRUE", plt = c(horiz.marg[1], horiz.marg[2], fig.marg[5], fig.marg[6]))
  h1 <- hist(te[te$taxon == "plant",which.cols[3]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2 <- hist(te[te$taxon == "plant",which.cols[4]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h3 <- hist(te[te$taxon == "anim",which.cols[3]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h4 <- hist(te[te$taxon == "anim",which.cols[4]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2$breaks <- -1 * h2$breaks
  ylimits <- c(0, vert.open.space * max(c(h1$counts+h3$counts, h2$counts+h4$counts)))
  plot(x = c(0,1), y = c(0,1), xlim = xlimits, ylim = ylimits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "", xaxs = "i", yaxs = "i")
  for (i in 1:length(h1$counts)) {
    rect(xleft = h1$breaks[i], ybottom = 0, xright = h1$breaks[i+1], ytop = h1$counts[i], lwd = lwidth, border = border.col, col = plant.pos.col)
    rect(xleft = h1$breaks[i], ybottom = h1$counts[i], xright = h1$breaks[i+1], ytop = h1$counts[i]+h3$counts[i], lwd = lwidth, border = border.col, col = anim.pos.col)
    rect(xleft = h2$breaks[i+1], ybottom = 0, xright = h2$breaks[i], ytop = h2$counts[i], lwd = lwidth, border = border.col, col = plant.neg.col)
    rect(xleft = h2$breaks[i+1], ybottom = h2$counts[i], xright = h2$breaks[i], ytop = h2$counts[i]+h4$counts[i], lwd = lwidth, border = border.col, col = anim.neg.col)
    }
  box(lwd = 1.3)
  axis(side = 1, tck = -0.015, labels = FALSE, at = c(-1,-0.5,0,0.5,1))
  mtext("Number of studies", side = 2, line = 1.7, at = 0, cex = axis.label.size )
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(c("10", "20", "30", "40"), side = 2, las = 1, line = 0.4, cex = tick.label.size, at = c(10,20,30,40) )
  # Add and label rectangles identifying animal and plant studies
  rect(-0.9, 30, -0.8, 33, lwd = lwidth, border = border.col, col = anim.neg.col)
  rect(-0.9, 25, -0.8, 28, lwd = lwidth, border = border.col, col = plant.neg.col)
  text("animals", x = -0.81, y = 31.5, pos = 4, cex = 0.8, col = anim.neg.col)
  text("plants", x = -0.81, y = 26.5, pos = 4, cex = 0.8, col = plant.neg.col)
  rect(0.8, 30, 0.9, 33, lwd = lwidth, border = border.col, col = anim.pos.col)
  rect(0.8, 25, 0.9, 28, lwd = lwidth, border = border.col, col = plant.pos.col)
  text("animals", x = 0.81, y = 31.5, pos = 2, cex = 0.8, col = anim.pos.col)
  text("plants", x = 0.81, y = 26.5, pos = 2, cex = 0.8, col = plant.pos.col)
  text("C. Survival", x = -1.03, y = 0.93 * ylimits[2], cex = panel.label.size, pos = 4)
# Local/foreign: plant size
par(new = "TRUE", plt = c(horiz.marg[1], horiz.marg[2], fig.marg[3], fig.marg[4]))
  h1 <- hist(te[te$taxon == "plant",which.cols[5]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2 <- hist(te[te$taxon == "plant",which.cols[6]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h3 <- hist(te[te$taxon == "anim",which.cols[5]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h4 <- hist(te[te$taxon == "anim",which.cols[6]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2$breaks <- -1 * h2$breaks
  ylimits <- c(0, vert.open.space * max(c(h1$counts+h3$counts, h2$counts+h4$counts)))
  plot(x = c(0,1), y = c(0,1), xlim = xlimits, ylim = ylimits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "", xaxs = "i", yaxs = "i")
  for (i in 1:length(h1$counts)) {
    rect(xleft = h1$breaks[i], ybottom = 0, xright = h1$breaks[i+1], ytop = h1$counts[i], lwd = lwidth, border = border.col, col = plant.pos.col)
    rect(xleft = h1$breaks[i], ybottom = h1$counts[i], xright = h1$breaks[i+1], ytop = h1$counts[i]+h3$counts[i], lwd = lwidth, border = border.col, col = anim.pos.col)
    rect(xleft = h2$breaks[i+1], ybottom = 0, xright = h2$breaks[i], ytop = h2$counts[i], lwd = lwidth, border = border.col, col = plant.neg.col)
    rect(xleft = h2$breaks[i+1], ybottom = h2$counts[i], xright = h2$breaks[i], ytop = h2$counts[i]+h4$counts[i], lwd = lwidth, border = border.col, col = anim.neg.col)
    }
  box(lwd = 1.3)
  axis(side = 1, tck = -0.015, labels = FALSE, at = c(-1,-0.5,0,0.5,1))
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(c("10", "20", "30", "40"), side = 2, las = 1, line = 0.4, cex = tick.label.size, at = c(10,20,30,40) )
  text("E. Size", x = -1.03, y = 0.93 * ylimits[2], cex = panel.label.size, pos = 4)
# Local/foreign: reproduction
par(new = "TRUE", plt = c(horiz.marg[1], horiz.marg[2], fig.marg[1], fig.marg[2]))
  h1 <- hist(te[te$taxon == "plant", which.cols[7]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2 <- hist(te[te$taxon == "plant", which.cols[8]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h3 <- hist(te[te$taxon == "anim", which.cols[7]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h4 <- hist(te[te$taxon == "anim", which.cols[8]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2$breaks <- -1 * h2$breaks
  ylimits <- c(0, vert.open.space * max(c(h1$counts+h3$counts, h2$counts+h4$counts)))
  plot(x = c(0,1), y = c(0,1), xlim = xlimits, ylim = ylimits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "", xaxs = "i", yaxs = "i")
  for (i in 1:length(h1$counts)) {
    rect(xleft = h1$breaks[i], ybottom = 0, xright = h1$breaks[i+1], ytop = h1$counts[i], lwd = lwidth, border = border.col, col = plant.pos.col)
    rect(xleft = h1$breaks[i], ybottom = h1$counts[i], xright = h1$breaks[i+1], ytop = h1$counts[i]+h3$counts[i], lwd = lwidth, border = border.col, col = anim.pos.col)
    rect(xleft = h2$breaks[i+1], ybottom = 0, xright = h2$breaks[i], ytop = h2$counts[i], lwd = lwidth, border = border.col, col = plant.neg.col)
    rect(xleft = h2$breaks[i+1], ybottom = h2$counts[i], xright = h2$breaks[i], ytop = h2$counts[i]+h4$counts[i], lwd = lwidth, border = border.col, col = anim.neg.col)
    }
  box(lwd = 1.3)
  axis(side = 1, tck = -0.015, labels = FALSE, at = c(-1,-0.5,0,0.5,1))
  mtext(c("-1.0", "-0.5", "0.0", "0.5", "1.0"), side = 1, line = 0.2, cex = tick.label.size, at = c(-1,-0.5,0,0.5,1) )
  mtext("Proportion of test sites", side = 1, line = 1.4, cex = axis.label.size )
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(c("0", "5", "10", "15", "20"), side = 2, las = 1, line = 0.4, cex = tick.label.size, at = c(0,5,10,15,20) )
  text("G. Fecundity", x = -1.03, y = 0.93 * ylimits[2], cex = panel.label.size, pos = 4)
#
# Right side of figure.
# Home/away: "fitness"
par(new = "TRUE", plt = c(horiz.marg[3], horiz.marg[4], fig.marg[7], fig.marg[8]))
  h1 <- hist(te[te$taxon == "plant", which.cols[9]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2 <- hist(te[te$taxon == "plant", which.cols[10]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h3 <- hist(te[te$taxon == "anim", which.cols[9]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h4 <- hist(te[te$taxon == "anim", which.cols[10]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2$breaks <- -1 * h2$breaks
  ylimits <- c(0, vert.open.space * max(c(h1$counts+h3$counts, h2$counts+h4$counts)))
  plot(x = c(0,1), y = c(0,1), xlim = xlimits, ylim = ylimits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "", xaxs = "i", yaxs = "i")
  for (i in 1:length(h1$counts)) {
    rect(xleft = h1$breaks[i], ybottom = 0, xright = h1$breaks[i+1], ytop = h1$counts[i], lwd = lwidth, border = border.col, col = plant.pos.col)
    rect(xleft = h1$breaks[i], ybottom = h1$counts[i], xright = h1$breaks[i+1], ytop = h1$counts[i]+h3$counts[i], lwd = lwidth, border = border.col, col = anim.pos.col)
    rect(xleft = h2$breaks[i+1], ybottom = 0, xright = h2$breaks[i], ytop = h2$counts[i], lwd = lwidth, border = border.col, col = plant.neg.col)
    rect(xleft = h2$breaks[i+1], ybottom = h2$counts[i], xright = h2$breaks[i], ytop = h2$counts[i]+h4$counts[i], lwd = lwidth, border = border.col, col = anim.neg.col)
    }
  box(lwd = 1.3)
  axis(side = 1, tck = -0.015, labels = FALSE, at = c(-1,-0.5,0,0.5,1))
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(c("2", "4", "6"), side = 2, las = 1, line = 0.4, cex = tick.label.size, at = c(2,4,6) )
  text("B. Fitness", x = -1.03, y = 0.93 * ylimits[2], cex = panel.label.size, pos = 4)
  mtext(c("Home-away criterion","across test sites"), side = 3, line = c(1.6,0.45), cex = axis.label.size+0.1 )
# Home/away: survival
par(new = "TRUE", plt = c(horiz.marg[3], horiz.marg[4], fig.marg[5], fig.marg[6]))
  h1 <- hist(te[te$taxon == "plant",which.cols[11]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2 <- hist(te[te$taxon == "plant",which.cols[12]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h3 <- hist(te[te$taxon == "anim",which.cols[11]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h4 <- hist(te[te$taxon == "anim",which.cols[12]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2$breaks <- -1 * h2$breaks
  ylimits <- c(0, vert.open.space * max(c(h1$counts+h3$counts, h2$counts+h4$counts)))
  plot(x = c(0,1), y = c(0,1), xlim = xlimits, ylim = ylimits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "", xaxs = "i", yaxs = "i")
  for (i in 1:length(h1$counts)) {
    rect(xleft = h1$breaks[i], ybottom = 0, xright = h1$breaks[i+1], ytop = h1$counts[i], lwd = lwidth, border = border.col, col = plant.pos.col)
    rect(xleft = h1$breaks[i], ybottom = h1$counts[i], xright = h1$breaks[i+1], ytop = h1$counts[i]+h3$counts[i], lwd = lwidth, border = border.col, col = anim.pos.col)
    rect(xleft = h2$breaks[i+1], ybottom = 0, xright = h2$breaks[i], ytop = h2$counts[i], lwd = lwidth, border = border.col, col = plant.neg.col)
    rect(xleft = h2$breaks[i+1], ybottom = h2$counts[i], xright = h2$breaks[i], ytop = h2$counts[i]+h4$counts[i], lwd = lwidth, border = border.col, col = anim.neg.col)
    }
  box(lwd = 1.3)
  axis(side = 1, tck = -0.015, labels = FALSE, at = c(-1,-0.5,0,0.5,1))
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(c("5", "10", "15", "20", "25"), side = 2, las = 1, line = 0.4, cex = tick.label.size, at = c(5,10,15,20,25) )
  text("D. Survival", x = -1.03, y = 0.93 * ylimits[2], cex = panel.label.size, pos = 4)
# Home/away: plant size
par(new = "TRUE", plt = c(horiz.marg[3], horiz.marg[4], fig.marg[3], fig.marg[4]))
  h1 <- hist(te[te$taxon == "plant",which.cols[13]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2 <- hist(te[te$taxon == "plant",which.cols[14]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h3 <- hist(te[te$taxon == "anim",which.cols[13]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h4 <- hist(te[te$taxon == "anim",which.cols[14]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2$breaks <- -1 * h2$breaks
  ylimits <- c(0, vert.open.space * max(c(h1$counts+h3$counts, h2$counts+h4$counts)))
  plot(x = c(0,1), y = c(0,1), xlim = xlimits, ylim = ylimits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "", xaxs = "i", yaxs = "i")
  for (i in 1:length(h1$counts)) {
    rect(xleft = h1$breaks[i], ybottom = 0, xright = h1$breaks[i+1], ytop = h1$counts[i], lwd = lwidth, border = border.col, col = plant.pos.col)
    rect(xleft = h1$breaks[i], ybottom = h1$counts[i], xright = h1$breaks[i+1], ytop = h1$counts[i]+h3$counts[i], lwd = lwidth, border = border.col, col = anim.pos.col)
    rect(xleft = h2$breaks[i+1], ybottom = 0, xright = h2$breaks[i], ytop = h2$counts[i], lwd = lwidth, border = border.col, col = plant.neg.col)
    rect(xleft = h2$breaks[i+1], ybottom = h2$counts[i], xright = h2$breaks[i], ytop = h2$counts[i]+h4$counts[i], lwd = lwidth, border = border.col, col = anim.neg.col)
    }
  box(lwd = 1.3)
  axis(side = 1, tck = -0.015, labels = FALSE, at = c(-1,-0.5,0,0.5,1))
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(c("10", "20", "30"), side = 2, las = 1, line = 0.4, cex = tick.label.size, at = c(10,20,30) )
  text("F. Size", x = -1.03, y = 0.93 * ylimits[2], cex = panel.label.size, pos = 4)
# Home/away: reproduction
par(new = "TRUE", plt = c(horiz.marg[3], horiz.marg[4], fig.marg[1], fig.marg[2]))
  h1 <- hist(te[te$taxon == "plant",which.cols[15]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2 <- hist(te[te$taxon == "plant",which.cols[16]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h3 <- hist(te[te$taxon == "anim",which.cols[15]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h4 <- hist(te[te$taxon == "anim",which.cols[16]], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), plot = F)
  h2$breaks <- -1 * h2$breaks
  ylimits <- c(0, vert.open.space * max(c(h1$counts+h3$counts, h2$counts+h4$counts)))
  plot(x = c(0,1), y = c(0,1), xlim = xlimits, ylim = ylimits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "", xaxs = "i", yaxs = "i")
  for (i in 1:length(h1$counts)) {
    rect(xleft = h1$breaks[i], ybottom = 0, xright = h1$breaks[i+1], ytop = h1$counts[i], lwd = lwidth, border = border.col, col = plant.pos.col)
    rect(xleft = h1$breaks[i], ybottom = h1$counts[i], xright = h1$breaks[i+1], ytop = h1$counts[i]+h3$counts[i], lwd = lwidth, border = border.col, col = anim.pos.col)
    rect(xleft = h2$breaks[i+1], ybottom = 0, xright = h2$breaks[i], ytop = h2$counts[i], lwd = lwidth, border = border.col, col = plant.neg.col)
    rect(xleft = h2$breaks[i+1], ybottom = h2$counts[i], xright = h2$breaks[i], ytop = h4$counts[i]+h2$counts[i], lwd = lwidth, border = border.col, col = anim.neg.col)
    }
  box(lwd = 1.3)
  axis(side = 1, tck = -0.015, labels = FALSE, at = c(-1,-0.5,0,0.5,1))
  mtext(c("-1.0", "-0.5", "0.0", "0.5", "1.0"), side = 1, line = 0.2, cex = tick.label.size, at = c(-1,-0.5,0,0.5,1) )
  mtext("Proportion of populations", side = 1, line = 1.4, cex = axis.label.size )
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(c("0", "4", "8", "12"), side = 2, las = 1, line = 0.4, cex = tick.label.size, at = c(0,4,8,12) )
  text("H. Fecundity", x = -1.03, y = 0.93 * ylimits[2], cex = panel.label.size, pos = 4)
dev.off()
#
#
##################################
# Characterize the environment along elevational gradients in previous studies and our study.
# Import MAT from Bioclim. This dataset comes from this page: https://worldclim.org/data/worldclim21.html
r.mat <- raster(" .... your path here .... /bio_1.bil")     # make a raster; units: °C * 10
#
# Import population locations.
#
pops         <- read.table("locations of source populations.txt", header = TRUE, stringsAsFactors = FALSE)
pops         <- pops[!is.na(pops$lat), ]     #   delete the single study for which population locations cannot be determined
pops$studyID <- paste(pops$authors, pops$species, sep="_")
study.list   <- unique(pops$studyID)
study.list   <- study.list[order(study.list)]
#
# Calculate the growing season for each site. The number of days over 5 degrees C.
# First import monthly WorldClim temperature data. The data come from: https://worldclim.org/data/worldclim21.html
# Reference: Fick SE, Hijmans RJ. 2017. WorldClim 2: new 1km spatial resolution climate surfaces for global land areas. International Journal of Climatology 37:4302-4315
# 1970-2000, resolution 30 seconds.
#
p1     <- data.frame(studyID = pops$studyID, pop = pops$pop, stringsAsFactors = FALSE)
rast   <- raster(" .... your path here .... /wc2.1_30s_tavg_01.tif")
p1$t01 <- extract(rast, SpatialPoints(pops[,c("long","lat")]) )
rast   <- raster(" .... your path here .... /wc2.1_30s_tavg_02.tif")
p1$t02 <- extract(rast, SpatialPoints(pops[,c("long","lat")]) )
rast   <- raster(" .... your path here .... /wc2.1_30s_tavg_03.tif")
p1$t03 <- extract(rast, SpatialPoints(pops[,c("long","lat")]) )
rast   <- raster(" .... your path here .... /wc2.1_30s_tavg_04.tif")
p1$t04 <- extract(rast, SpatialPoints(pops[,c("long","lat")]) )
rast   <- raster(" .... your path here .... /wc2.1_30s_tavg_05.tif")
p1$t05 <- extract(rast, SpatialPoints(pops[,c("long","lat")]) )
rast   <- raster(" .... your path here .... /wc2.1_30s_tavg_06.tif")
p1$t06 <- extract(rast, SpatialPoints(pops[,c("long","lat")]) )
rast   <- raster(" .... your path here .... /wc2.1_30s_tavg_07.tif")
p1$t07 <- extract(rast, SpatialPoints(pops[,c("long","lat")]) )
rast   <- raster(" .... your path here .... /wc2.1_30s_tavg_08.tif")
p1$t08 <- extract(rast, SpatialPoints(pops[,c("long","lat")]) )
rast   <- raster(" .... your path here .... /wc2.1_30s_tavg_09.tif")
p1$t09 <- extract(rast, SpatialPoints(pops[,c("long","lat")]) )
rast   <- raster(" .... your path here .... /wc2.1_30s_tavg_10.tif")
p1$t10 <- extract(rast, SpatialPoints(pops[,c("long","lat")]) )
rast   <- raster(" .... your path here .... /wc2.1_30s_tavg_11.tif")
p1$t11 <- extract(rast, SpatialPoints(pops[,c("long","lat")]) )
rast   <- raster(" .... your path here .... /wc2.1_30s_tavg_12.tif")
p1$t12 <- extract(rast, SpatialPoints(pops[,c("long","lat")]) )
p1$GS.5 <- growing.season(lay.date = NA, sites = p1$studyID, dat = p1[, 3:14], threshold = 5)$season.length
# Plot worldwide temperature map. Put the sampled localities on top of it,.
# plot(rast)
# map("world", add=TRUE)
# points(SpatialPoints(pops[,c("long","lat")]), col = "red" )
#
# Get the MAT from bioclim
#
p1$MAT <- extract(r.mat, SpatialPoints(pops[,c("long","lat")]) )/10
#
# This isolates those populations for which the coordinates are over water. There should be none.
missing.mat <- pops[is.na(pops$MAT),]
missing.mat[ , c("study.nr", "pop", "species", "lat", "long", "elev")]
#
p1     <- p1[, -c(3:14)]
p2     <- merge(pops[,c("studyID","pop","lat","long","elev")], p1, by = c("studyID","pop"))
head(p2)
#
# Now calculate gradient length and steepness for MAT and Growing Season
#
res <- data.frame(studyID = study.list, n.pops = NA, elev.range = NA, km.dist.max = NA, range.MAT = NA, km.dist.MAT = NA, range.GS = NA, km.dist.GS = NA, stringsAsFactors = FALSE)
for (i in 1:length(study.list)) {
  d.sub <- subset(p2, p2$studyID == study.list[i])
  res$studyID[i]     <- study.list[i]
  res$n.pops[i]      <- dim(d.sub)[1]
  all.distances      <- rdist.earth(cbind(d.sub$long, d.sub$lat), cbind(d.sub$long, d.sub$lat), miles = FALSE, R = NULL)
  res$elev.range[i]  <- max(dist(d.sub$elev))
  res$km.dist.max[i] <- max(all.distances)
  # Now Mean Annual Temperature -- range and steepness
  WhichMin           <- which.min(d.sub$MAT)
  WhichMax           <- which.max(d.sub$MAT)
  res$range.MAT[i]   <- d.sub$MAT[WhichMax] - d.sub$MAT[WhichMin]
  res$km.dist.MAT[i] <- rdist.earth(cbind(d.sub$long[WhichMax], d.sub$lat[WhichMax]), cbind(d.sub$long[WhichMin], d.sub$lat[WhichMin]), miles = FALSE, R = NULL)
  # Now Growing Season (5 degrees cutoff) -- range and steepness
  WhichMin           <- which.min(d.sub$GS.5)
  WhichMax           <- which.max(d.sub$GS.5)
  res$range.GS[i]    <- d.sub$GS.5[WhichMax] - d.sub$GS.5[WhichMin]
  res$km.dist.GS[i]  <- rdist.earth(cbind(d.sub$long[WhichMax], d.sub$lat[WhichMax]), cbind(d.sub$long[WhichMin], d.sub$lat[WhichMin]), miles = FALSE, R = NULL)
  }
res$steep.MAT  <- res$range.MAT / res$km.dist.MAT
res$steep.GS   <- res$range.GS / res$km.dist.GS
res$steep.GS   <- ifelse(res$km.dist.GS == 0, 0, res$steep.GS)  # situation where Growing Season is 365 days. 
res$elev.range <- format.p.val(res$elev.range, 0)
res[,c(4,6,8,9,10)] <- apply(res[,c(4,6,8,9,10)], 2, format.p.val, 3)
te$studyID     <- paste(te$author, te$species, sep="_")
summry         <- merge(te, res, by = "studyID", all = TRUE)
summry[,c(7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29)] <- apply(summry[,c(7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29)], 2, format.p.val, 2)
summry$G.fit   <- paste(summry$G.pos.fit, summry$G.neg.fit, sep = "/")
summry$G.surv  <- paste(summry$G.pos.surv, summry$G.neg.surv, sep = "/")
summry$G.size  <- paste(summry$G.pos.size, summry$G.neg.size, sep = "/")
summry$G.repr  <- paste(summry$G.pos.repr, summry$G.neg.repr, sep = "/")
summry$P.fit   <- paste(summry$P.pos.fit, summry$P.neg.fit, sep = "/")
summry$P.surv  <- paste(summry$P.pos.surv, summry$P.neg.surv, sep = "/")
summry$P.size  <- paste(summry$P.pos.size, summry$P.neg.size, sep = "/")
summry$P.repr  <- paste(summry$P.pos.repr, summry$P.neg.repr, sep = "/")
s1             <- summry
my.fun         <- function(x) { ifelse(x == "NA/NA", "-/-", x) }
summry[,c(41:48)] <- apply(summry[,c(41:48)], 2, my.fun)
summry$km.dist.max <- format.p.val(summry$km.dist.max, 1)
summry         <- summry[order(summry$taxon, summry$species), ]
summry         <- summry[,c("species", "G.fit", "G.surv", "G.size", "G.repr", "P.fit", "P.surv", "P.size", "P.repr", "km.dist.max", "elev.range", "range.MAT", "range.GS", "steep.MAT", "steep.GS", "author")]
summry
#
#
# Now plot histograms of environmental steepness (MAT and Growing Season > 5 degrees),
# showing how our study compares.
# Fig. 5 in the Evolution MS.
#
hist.breaks     <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4)
vert.open.space <- 1.04
border.col      <- NULL
axis.label.size <- 1.1
pdf(" .... your path here .... /Fig 5.pdf", width = 6, height = 8, useDingbats = FALSE)
plot.new()
# Panel A: range in GS
  par(new = "TRUE", plt = c(0.15,0.9, 0.55,0.9))
  h1 <- hist(s1$range.GS[s1$taxon == "plant"], plot = FALSE, breaks = 193 * hist.breaks )
  h2 <- hist(s1$range.GS[s1$taxon == "anim"], plot = FALSE, breaks = 193 * hist.breaks )
  xlimits <- c(-2, 272.5)
  ylimits <- c(0, vert.open.space * max(c(h1$counts+h2$counts)))
  plot(x = c(0,1), y = c(0,1), xlim = xlimits, ylim = ylimits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "", xaxs = "i", yaxs = "i")
  for (i in 1:length(h1$counts)) {
    rect(xleft = h1$breaks[i], ybottom = 0, xright = h1$breaks[i+1], ytop = h1$counts[i], lwd = lwidth, border = border.col, col = plant.pos.col)
    rect(xleft = h1$breaks[i], ybottom = h1$counts[i], xright = h1$breaks[i+1], ytop = h1$counts[i]+h2$counts[i], lwd = lwidth, border = border.col, col = plant.neg.col)
    }
  box(lwd = 1.3)
  axis(side = 1, tck = -0.02, labels = FALSE)
  mtext(c("0", "50", "100", "150", "200", "250"), side = 1, line = 0.2, cex = tick.label.size, at = c(0,50,100,150,200,250) )
  mtext("Range in growing season (days)", side = 1, line = 1.4, cex = axis.label.size )
  axis(side = 2, tck = -0.01, labels = FALSE)
  mtext(c("0", "5", "10", "15"), side = 2, las = 1, line = 0.4, cex = tick.label.size, at = c(0,5,10,15) )
  arrows(s1$range.GS[s1$author == "Bachmann.VanBuskirk.2021"], 0.7 * ylimits[2], s1$range.GS[s1$author == "Bachmann.VanBuskirk.2021"], 0.5 * ylimits[2], lwd = 1.3, length = 0.12)
  text(c("this", "study"), x = s1$range.GS[s1$author == "Bachmann.VanBuskirk.2021"], y = c(0.77 * ylimits[2], 0.7 * ylimits[2]), pos = 3, cex = 1)
  text("A", x = 0.955 * xlimits[2], y = 0.93 * ylimits[2], cex = 1.5)
  mtext("Figure 5             drawn by 'analysis of Judith transplant.R'", side = 3, line = 2.4, cex = 0.6 )
#
# Panel B: growing season steepness
  par(new = "TRUE", plt = c(0.15,0.9, 0.1,0.45))
  h1 <- hist(log10(1+s1$steep.GS[s1$taxon == "plant"]), plot = FALSE, breaks = hist.breaks)
  h2 <- hist(log10(1+s1$steep.GS[s1$taxon == "anim"]), plot = FALSE, breaks = hist.breaks)
  xlimits <- c(-0.011, 1.412)
  ylimits <- c(0, vert.open.space * max(c(h1$counts+h2$counts)))
  plot(x = c(0,1), y = c(0,1), xlim = xlimits, ylim = ylimits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "", xaxs = "i", yaxs = "i")
  for (i in 1:length(h1$counts)) {
    rect(xleft = h1$breaks[i], ybottom = 0, xright = h1$breaks[i+1], ytop = h1$counts[i], lwd = lwidth, border = border.col, col = plant.pos.col)
    rect(xleft = h1$breaks[i], ybottom = h1$counts[i], xright = h1$breaks[i+1], ytop = h1$counts[i]+h2$counts[i], lwd = lwidth, border = border.col, col = plant.neg.col)
    }
  box(lwd = 1.3)
  axis(side = 1, tck = -0.02, labels = FALSE, at = log10(1 + c(0,1,10, 20)))
  axis(side = 1, tck = -0.01, labels = FALSE, at = log10(1 + c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,2,3,4,5,6,7,8,9)))
  mtext(c("0", "1", "10", "20"), side = 1, line = 0.2, cex = tick.label.size, at = log10(1+c(0,1,10,20)) )
  mtext("Growing season steepness (1 + days/km)", side = 1, line = 1.4, cex = axis.label.size )
  axis(side = 2, tck = -0.01, labels = FALSE)
  mtext(c("0", "4", "8", "12"), side = 2, las = 1, line = 0.4, cex = tick.label.size, at = c(0,4,8,12) )
  mtext("Number of studies", side = 2, line = 1.6, at = 16.2, cex = axis.label.size )
  arrows(log10(1 + s1$steep.GS[s1$author == "Bachmann.VanBuskirk.2021"]), 0.7 * ylimits[2], log10(1 + s1$steep.GS[s1$author == "Bachmann.VanBuskirk.2021"]), 0.5 * ylimits[2], lwd = 1.3, length = 0.12)
  text(c("this", "study"), x = log10(1 + s1$steep.GS[s1$author == "Bachmann.VanBuskirk.2021"]), y = c(0.77 * ylimits[2], 0.7 * ylimits[2]), pos = 3, cex = 1)
  x.position <- c(xlimits[1] + 0.9*diff(xlimits), xlimits[1] + 0.93*diff(xlimits))
  y.position <- c(0.48*ylimits[2], 0.64*ylimits[2])
  y.interval <- 0.025 * ylimits[2]
  rect(x.position[1], y.position[2]-y.interval, x.position[2], y.position[2]+y.interval, lwd = lwidth, border = border.col, col = plant.neg.col)
  rect(x.position[1], y.position[1]-y.interval, x.position[2], y.position[1]+y.interval, lwd = lwidth, border = border.col, col = plant.pos.col)
  text(c("animals", expression(paste(italic("N"), " = 8"))), x = x.position[1], y = c(y.position[2]+1.2*y.interval, y.position[2]-1.2*y.interval), pos = 2, cex = 1, col = plant.neg.col)
  text(c("plants", expression(paste(italic("N"), " = 71"))), x = x.position[1], y = c(y.position[1]+1.2*y.interval, y.position[1]-1.2*y.interval), pos = 2, cex = 1, col = plant.pos.col)
  text("B", x = 0.955 * xlimits[2], y = 0.93 * ylimits[2], cex = 1.5)
dev.off()
#
#
# Figure out the range in elevation and GS in our study. Where do these values fall relative to other studies?
#
steep.adapt <- data.frame(species = s1$species, studyID = s1$studyID, taxon = s1$taxon, elev.range = s1$elev.range, GS.range = s1$range.GS, MAT.range = s1$range.MAT, steep.MAT = s1$steep.MAT, steep.GS = s1$steep.GS, G.pos.mean = apply(s1[ , names(s1)[grep("G[.]pos", names(s1))]], 1, mean, na.rm = TRUE), G.neg.mean = apply(s1[ , names(s1)[grep("G[.]neg", names(s1))]], 1, mean, na.rm = TRUE), P.pos.mean = apply(s1[ , names(s1)[grep("P[.]pos", names(s1))]], 1, mean, na.rm = TRUE), P.neg.mean = apply(s1[ , names(s1)[grep("P[.]neg", names(s1))]], 1, mean, na.rm = TRUE), stringsAsFactors = FALSE)
proc.means(steep.adapt, cats = c("taxon"), c("G.pos.mean", "G.neg.mean", "P.pos.mean", "P.neg.mean"), stats = c("N", "mean", "SD"))
(our.stdy <- s1$elev.range[s1$author == "Bachmann.VanBuskirk.2021"])
1 - sum(s1$elev.range > our.stdy, na.rm = TRUE) / length(s1$elev.range > our.stdy)
(our.stdy <- s1$range.GS[s1$author == "Bachmann.VanBuskirk.2021"])
1 - sum(s1$range.GS > our.stdy, na.rm = TRUE) / length(s1$range.GS > our.stdy)
(our.stdy <- s1$steep.GS[s1$author == "Bachmann.VanBuskirk.2021"])
1 - sum(s1$steep.GS > our.stdy, na.rm = TRUE) / length(s1$steep.GS > our.stdy)
#
cor.test(steep.adapt$P.pos.mean, steep.adapt$G.pos.mean)
#
#######################
#
# Draw Fig. S10.
#
mII             <- modelII(steep.adapt$P.pos.mean, steep.adapt$G.pos.mean)   # major axis regression
x.limits        <- y.limits <- c(-0.07, 1.03)
edge.width      <- 1
tick.label.size <- 1
axis.label.size <- 1.2
jitter.amount   <- 2.5
sa              <- steep.adapt[order(steep.adapt$GS.range, decreasing = TRUE), c("taxon", "GS.range", "G.pos.mean", "P.pos.mean")]
sa              <- sa[complete.cases(sa), ]
symbol.size     <- 1.1 + (sa$GS.range / 70)
sa$j.G.pos.mean <- sa$G.pos.mean
sa$j.P.pos.mean <- sa$P.pos.mean
sa$x            <- paste(sa$G.pos.mean, sa$P.pos.mean, sep="_")
uniq.x          <- unique(sa$x)
set.seed(125)
for (i in 1:length(uniq.x)) {
  d.sub <- subset(sa, x == uniq.x[i])
  if (dim(d.sub)[1] > 1) {
    sa$j.P.pos.mean[which(sa$x == uniq.x[i])] <- jitter(d.sub$P.pos.mean, factor = jitter.amount)
    sa$j.G.pos.mean[which(sa$x == uniq.x[i])] <- jitter(d.sub$G.pos.mean, factor = jitter.amount)
    }
  }
sa$color <- ifelse(sa$taxon == "plant", plant.pos.col, anim.pos.col)
pdf(" .... your path here .... /Fig S10.pdf", width = 6, height = 6, useDingbats = FALSE)
plot.new()
par(new = "TRUE", plt = c(0.17, 0.87, 0.12, 0.82))
  plot(x = c(0,1), y = c(0,1), xlim = x.limits, ylim = y.limits, xlab = "", ylab = "", tck = 0, xaxt="n", yaxt="n", pch = "", xaxs = "i", yaxs = "i")
  add.short.line(a = mII$Intercept, b = mII$betacoeff, x.vals = sa$P.pos.mean, lwidth = 1.3, ltype = 2)
  points(sa$j.P.pos.mean, sa$j.G.pos.mean, cex = symbol.size, col = "black", bg = sa$color, lwd = edge.width, pch = 21)
  box(lwd = 1.3)
  axis(side = 1, tck = -0.015, labels = FALSE)
  mtext(c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), side = 1, line = 0.3, cex = tick.label.size, at = c(0,0.2,0.4,0.6,0.8,1) )
  mtext(c("Proportion of populations", "performing better at home elevation"), side = 1, line = c(1.4,2.4), cex = axis.label.size )
  axis(side = 2, tck = -0.015, labels = FALSE)
  mtext(c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), side = 2, line = 0.45, las = 1, cex = tick.label.size, at = c(0,0.2,0.4,0.6,0.8,1) )
  mtext(c("Proportion of test sites", "showing a home elevation advantage"), side = 2, line = c(3,2), cex = axis.label.size )
  # arrow indicating our study.
  arrows(0.59, 0.75, 0.533, 0.693, length = 0.1, lwd = 1.3)
  text(c("this", "study"), x = 0.615, y = c(0.81, 0.775), cex = 0.9)
  par(xpd = TRUE)
  # animals versus plants at the top
  bubble.sizes <- 1.1 + (c(0, 125, 250) / 70)
  points(x = c(x.limits[1] + 0.07, x.limits[1] + 0.07), y = c(1.08, 1.145), cex = bubble.sizes[2]-0.2, lwd = edge.width, bg = c(plant.pos.col, anim.pos.col), pch = 21)
  text(c("plants", "animals"), cex = 1, x = c(x.limits[1] + 0.09, x.limits[1] + 0.09), y = c(1.08, 1.145), pos = 4)
  # scale at the top
  points(x = c(0.35, 0.475, 0.6), y = c(1.09, 1.09, 1.09), cex = bubble.sizes, lwd = edge.width, bg = plant.pos.col, pch = 21)
  lines(x = c(0.35, 0.6), y = c(1.155, 1.155), lwd = 1.3)      # caption line at top
  lines(x = c(0.35, 0.35), y = c(1.155, 1.165), lwd = 1.3)       # tick marks
  lines(x = c(0.475, 0.475), y = c(1.155, 1.165), lwd = 1.3)
  lines(x = c(0.6, 0.6), y = c(1.155, 1.165), lwd = 1.3)
  mtext(c("0", "125", "250"), side = 3, at = c(0.35, 0.475, 0.6), line = 2.6, cex = 0.9)
  mtext(c("range in growing", "season (days)"), side = 3, at = 0.67, line = c(1.8, 0.8), adj = 0, cex = 1)
  par(xpd = FALSE)
  mtext("Fig. S10                   drawn by 'analysis of Judith transplant.R'", side = 3, line = 4.3, cex = 0.6)
dev.off()












##########################################
# Load these functions first
##########################################
#
#
zi.predict <- function(fixed.pois, fixed.zi, v.treatment, v.test.elev, binom.intercept = NULL) {
  # Custom function for returning the predicted value from a ZIP model in MCMCglmm.
  # The current model is: round.fitness ~ trait + trait:test.elev + trait:treatment + trait:treatment:test.elev.
  # Args:
  #  fixed.pois and fixed.zi are MCMC estimates of fixed effects for the Poisson and
  #          zero-inflated parts of the model. Each is a column vector as long as niter.
  #  v.treatment, v.test.elev, binom.intercept
  #  v.ttreatment is the desired value of the fixed effect "treatment". Possible values are "A_home", "B_away.pond", "C_away.elev".
  #  v.test.elev is the desired value of the fixed effect "test.elev". Possible values are "low" and "high".
  #  binom.intercept is an alternative intercept that can be included in the binomical part of the
  #          model. Default is NULL. I do this because the estimate from MCMCglmm is way too small.
  #
  n.reps <- dim(fixed.pois)[1]
  #
  # Set dummy variables
  # test.elev dummy variable:
  test.elev.dm <- rep(0, n.reps)
  if(v.test.elev == "low") test.elev.dm <- rep(1, n.reps)
  #
  # treatment dummy variables:
  treatmentB.dm <- treatmentC.dm <- rep(0, n.reps)
  if(v.treatment == "B_away.pond") treatmentB.dm <- rep(1, n.reps)
  if(v.treatment == "C_away.elev") treatmentC.dm <- rep(1, n.reps)
  #
  # calculate the two parts
  pois.part <- fixed.pois[,1] + fixed.pois[,2]*test.elev.dm + fixed.pois[,3]*treatmentB.dm + fixed.pois[,4]*treatmentC.dm + fixed.pois[,5]*test.elev.dm*treatmentB.dm + fixed.pois[,6]*test.elev.dm*treatmentC.dm
  if(!is.null(binom.intercept)) zi.part <- binom.intercept + fixed.zi[,2]*test.elev.dm   + fixed.zi[,3]*treatmentB.dm   + fixed.zi[,4]*treatmentC.dm   + fixed.zi[,5]*test.elev.dm*treatmentB.dm   + fixed.zi[,6]*test.elev.dm*treatmentC.dm
  if(is.null(binom.intercept))  zi.part <- fixed.zi[,1]   + fixed.zi[,2]*test.elev.dm   + fixed.zi[,3]*treatmentB.dm   + fixed.zi[,4]*treatmentC.dm   + fixed.zi[,5]*test.elev.dm*treatmentB.dm   + fixed.zi[,6]*test.elev.dm*treatmentC.dm
  #
  res <- data.frame(poisson.part = exp(pois.part), poisson.zero = dpois(0, exp(pois.part)), plogis.l = plogis(zi.part), plogis.minusl = plogis(-zi.part))
  colnames(res) <- c("poisson.part", "poisson.zero", "plogis.l", "plogis.minusl")
  res$prob.zero <- res$plogis.l + res$plogis.minusl*res$poisson.zero
  res$predict.fitness <- res$poisson.part * (1 - res$prob.zero)
  df.sum <- c(mean(res$predict.fitness), quantile(res$predict.fitness, c(0.5, 0.025, 0.975)))
  names(df.sum)[1:2] <- c("mean", "median")
  list(summary = df.sum, predict.fit = res$predict.fitness)
  }
#
#########################################
#
#
month.day.to.Julian <- function(month, day, year = 1981) {
  # Function to calculate the Julian date of the year, given a month and day.
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
#
#########################################
#
#
rename <- function(dat, oldnames, newnames) {
     # This function takes a data frame, renames one or more columns, and returns the data frame.
     # Examples of calling this function:
     #    rename(data, "old1", "new1")
     #    rename(data, c("old1", "old2"), c("new1", "new2"))
     #
  if( length(oldnames) != length(newnames) ) {
    message("Error: oldnames and newnames must have the same length")
    flush.console()
    break
    }
  for(i in 1:length(oldnames)) {
    colnames(dat)[which(names(dat) == oldnames[i])] <- newnames[i]
    }
  return(dat)
  }
#
#########################################
#
#
rescale.within.groups <- function(data, group.list, response.list, scale = TRUE) {
  # Function to rescale/standardize/normalize some of the variables in a data frame.
  # Arguments:
  #  data          = original data frame.
  #  group.list    = names of variables by which rescaling is done.
  #  response.list = names of variables that should be rescaled.
  #  scale         = rescaled SD=1 (TRUE) or centered only (FALSE)
  #
  new.data <- data
  n.group.vars   <- length(group.list)
  group.var.cols <- match(group.list, names(new.data))
  n.responses    <- length(response.list)
  response.cols  <- match(response.list, colnames(new.data))
      # set up a new variable that defines all the groups
  data$new.group <- factor(data[,group.var.cols[1]])
  if (n.group.vars > 1) {
    for (i in 2:n.group.vars) {
      data$new.group <- paste(data$new.group, data[,group.var.cols[i]], sep=".")
      }
    }
      # call function "group.normalize" to normalize each variable in turn by the new grouping variable
  for (i in response.cols) {
    if (scale==TRUE)  new.data[,i] <- as.numeric(group.normalize(data[,i], data$new.group ))
    if (scale==FALSE) new.data[,i] <- as.numeric(group.center(data[,i], data$new.group ))
    }      # end looping through i response variables
  return(new.data)
  }        # end function "normalize.within.groups"
#
#########################################
#
#
group.normalize <- function(var,grp) { return((var - tapply(var, grp, mean, na.rm = TRUE)[grp])/tapply(var, grp, sd, na.rm = TRUE)[grp]) }
#
#########################################
#
#
arrow.lines <- function(vect1, which.arrow, y.lims, arr.col, arr.len = 0.1, text.size, digits = 3) {
  # This function draws an arrow to a line, and text to the right.
  # Arguments:
  #  vect1       = a vector of numbers from which the median should be calculated (this will be the end of the arrow).
  #  which.arrow = "parallel", "unique", or "total".
  #  y.limits    = two numbers specifying limits of the y axis.
  #  arr.col     = color.
  #  arr.len     = arrow length.
  #  text.size   = cex for the text.
  #  digits      = how many digits to show when printing the median value.
  #  
  med         <- median(vect1, na.rm = TRUE)
  y.range     <- y.limits[2] - y.limits[1]
  y.line.ends <- y.range/25
  if(which.arrow == "parallel") {
    y1 <- 0.95 * y.range
    lines(x = c(med, med), y = c(y1 - y.line.ends, y1 + y.line.ends), col = arr.col)
    }
  if(which.arrow == "unique") {
    y1 <- 0.85 * y.range
    lines(x = c(med, med), y = c(y1 - y.line.ends, y1 + y.line.ends), col = arr.col)
    }
  if(which.arrow == "total") {
    y1 <- 0.75 * y.range
    lines(x = c(med, med), y = c(y1 - y.line.ends, y1 + y.line.ends), col = arr.col)
    }
  arrows(0, y1, med-0.2, y1, col = arr.col, arr.len)
  text(x = med - 0.5, y = y1, labels = paste(format(med, digits = digits), ", ", which.arrow, sep = ""), col = arr.col, cex = text.size, pos = 4)
  }
#
#########################################
#
#
remove.outliers <- function(x, med.mean = "med", units = 2, remove.na = "remove") {
   # This function is in "various useful functions.R"
   # This function removes outliers that are greater than a certain number of SD units away from the mean or median.
   #
   # Arguments:
   #  x         = a vector of observations from which outliers should be removed.
   #  med.mean  = "mean" or "med", indicating whether outliers should be measured from the median or mean.
   #  units     = threshold distance from the mean/med used to gauge outliers.
   #  remove.na = "remove" or "na", indicating whether outliers should be removed from x or replaced with NA.
   #
   std.dev     <- sd(x, na.rm = TRUE)
   threshold.x <- median(x, na.rm = TRUE)
   if(med.mean == "mean") threshold.x <- mean(x, na.rm = TRUE)
   scaled.dev <- abs(x - threshold.x) / std.dev
   if(remove.na == "remove") x <- x[ !(scaled.dev > units) ]    # units SD away from the med/mean is considered an outlier and removed.
   if(remove.na == "na") x[ (scaled.dev > units) ] <- NA        # units SD away from the med/mean is considered an outlier and replaced with NA.
   x
   }
#
#########################################
#
#
make.MCMCglmm.summary <- function(m.full, m2, m3, m4) {
  # One-off function to make a table summarizing the MCMCglmm models in the transplant experiment
  # Arguments are four models: full, no cage, no cage or block, and no cage, block, or test pond.
  LR.cage  <- mean(m2$Deviance) - mean(m.full$Deviance)
  LR.block <- mean(m3$Deviance) - mean(m2$Deviance)
  LR.pond  <- mean(m4$Deviance) - mean(m3$Deviance)
  summary.ran.eff          <- data.frame(random.effect = c("pond", "block", "cage"), Deviance = c(LR.pond, LR.block, LR.cage), P.value = c(1 - pchisq(LR.pond, 1), 1 - pchisq(LR.block, 1), 1 - pchisq(LR.cage, 1)), stringsAsFactors = FALSE)
  summary.ran.eff$Deviance <- format.p.val(summary.ran.eff$Deviance, 2)
  summary.ran.eff$P.value  <- format.p.val(x = summary.ran.eff$P.value, 5)
  sum1                     <- summary(m.full)
  if(m.full$Residual$original.family == "categorical") {
    s1 <- sum1$solutions
    r1 <- sum1$Gcovariances
    CI <- paste("(", format.p.val(s1[,2], 3), ", ", format.p.val(s1[,3], 3), ")", sep = "")
    df <- data.frame(process = c(rep("binomial",7), ""),
                       source = c("intercept", "test.elev", "local.foreign", "test.by.local", "pond  ", "block  ", "cage  ", paste("N =", length(m.full$family))),
                       level = c("", "low", "local", "low.local", rep("", 4)),
                       estimate = c(s1[1:4,1], r1[c(1:3),1], " "),
                       HPDinterval = c(CI[1:4],"", "", "", ""),
                       test.statistic = c(rep("", 4), summary.ran.eff$Deviance, ""),
                       MCMCp = c(s1[1:4,5], summary.ran.eff$P.value, " "))
    }
  if(m.full$Residual$original.family == "gaussian") {
    s1 <- sum1$solutions
    r1 <- sum1$Gcovariances
    CI <- paste("(", format.p.val(s1[,2], 3), ", ", format.p.val(s1[,3], 3), ")", sep = "")
    df <- data.frame(process = c(rep("gaussian",7), ""),
                       source = c("intercept", "test.elev", "local.foreign", "test.by.local", "pond  ", "block  ", "cage  ", paste("N =", length(m.full$family))),
                       level = c("", "low", "local", "low.local", "", "", "", ""),
                       estimate = c(s1[1:4,1], r1[c(1:3),1], " "),
                       HPDinterval = c(CI[1:4],"", "", "", ""),
                       test.statistic = c(rep("", 4), summary.ran.eff$Deviance, ""),
                       MCMCp = c(s1[1:4,5], summary.ran.eff$P.value, " ") )
    }
  if(m.full$Residual$original.family == "zipoisson") {
    s1 <- sum1$solutions[c(2,4,6,8,1,3,5,7),c(1,2,3,5)]
    r1 <- sum1$Gcovariances
    CI <- paste("(", format.p.val(s1[,2], 3), ", ", format.p.val(s1[,3], 3), ")", sep = "")
    df <- data.frame(process = c(rep("ZI",7), rep("Pois",7), ""),
                     source = c("intercept", "test.elev", "local.foreign", "test.by.local", "pond  ", "block  ", "cage  ", "intercept", "test.elev", "local.foreign", "test.by.local", "pond  ", "block  ", "cage  ", paste("N =", length(m.full$family)/2)),
                     level = c("", "low", "local", "low.local", "", "", "", "", "low", "local", "low.local", "", "", "", ""),
                     estimate = c(s1[1:4,1], r1[c(2,4,6),1], s1[5:8,1], r1[c(1,3,5),1], " "),
                     HPDinterval = c(CI[1:4],"", "", "", CI[5:8], "", "", "", ""),
                     test.statistic = c(rep("", 11), summary.ran.eff$Deviance, ""),
                     MCMCp = c(s1[1:4,4], "", "", "", s1[5:8,4], summary.ran.eff$P.value, " ") )
    }
  rownames(df) <- NULL
  df$estimate  <- format.p.val(df$estimate, 4); df$MCMCp <- format.p.val(df$MCMCp, 4)
  return(df)
  }
#
#########################################
#
#
format.p.val <- function(x, ndigits = 5) {
  # Located in "functions for meta-analysis" and "various useful functions".
  # Format p-values to a specified number of digits.
  # x is a number or vector.
  ff <- format(x, scientific = FALSE)
  suppressWarnings(round((as.numeric(ff) * (10^ndigits)))/(10^ndigits))
  }
#
#########################################
#
#
group.center <- function(var,grp) { return(var - tapply(var, grp, mean, na.rm = TRUE)[grp]) }
#
#########################################
#
#
proc.means <- function(data, cats, vars, id.vars = NULL, stats) {
  # Function to perform a SAS-like "proc means" calculation of summary stats. Columns should be properly named.
  # Args:
  #  data    = a data frame containing all the data.
  #  cats    = the variables defining categories/classes within which summary stats should be calculated (can be variable names or column numbers).
  #            Currently, cats must be specified. Can be assigned to an invariant column in the data set, like this: data = cbind(data, data.frame(const = 1))
  #  vars    = the variables for which summary stats should be calculated (can be variable names or column numbers).
  #  id.vars = extra variables to be carried over (can be variable names or column numbers). Should have constant values within cats.
  #  stats   = names of the statistics to calculate. Currently can include the keywords: "N" "n" "mean" "median" "med" "SD" "sd" "VAR" or "var".
  #
  # Use it like this:
  # data <- proc.means(data, vars = "response", cats = c("year", "species"), stats = c("N", "mean", "SD"))
  #
  # Make sure the statistics we're requesting are on the list of possibilities.
  for(i in 1:length(stats)) {
    if(!stats[i] %in% c("N", "n", "mean", "median", "med", "SD", "sd", "VAR", "var")) stop('One or more requested stats is not on the list of possibilities \\n\\n') 
    }
  #
  # Prepare column lists corresponding to variables, categories (classes), and id variables.
  if(is.character(cats[1]))    cat.cols <- match(cats, names(data))
  if(is.character(vars[1]))    var.cols <- match(vars, names(data))
  if(is.character(id.vars[1])) id.cols  <- match(id.vars, names(data))
  if(is.numeric(cats[1]))      cat.cols <- cats
  if(is.numeric(vars[1]))      var.cols <- vars
  if(is.numeric(id.vars[1]))   id.cols  <- id.vars
  if(exists("id.cols"))        cat.cols <- c(cat.cols, id.cols)
  cat.names     <- colnames(data)[cat.cols]
  category.list <- list()
  for (i in 1:length(cat.cols)) category.list[[i]] <- data[,cat.cols[i]]
  #
  # Remove these data frames if they already exist.
  if(exists("d.agg.mean")) rm(d.agg.mean)
  if(exists("d.agg.med")) rm(d.agg.med)
  if(exists("d.agg.SD")) rm(d.agg.SD)
  if(exists("d.agg.var")) rm(d.agg.var)
  if(exists("d.agg.N")) rm(d.agg.N)
  #
  # Perform the aggregation:
  my.funct <- function(x) { length(x[!is.na(x)]) }    # function to count only non-missing observations
  datalist <- list()
  if("N" %in% stats | "n" %in% stats) {
    d.agg.N <- aggregate(data[,var.cols], by = category.list, my.funct)
    colnames(d.agg.N) <- c(cat.names, paste(vars, "N", sep = "_") )
    datalist$N <- d.agg.N
    }
  if("mean" %in% stats) {
    d.agg.mean <- aggregate(data[,var.cols], by = category.list, mean, na.rm = TRUE)
    colnames(d.agg.mean) <- c(cat.names,  paste(vars, "mean", sep = "_"))
    datalist$mean <- d.agg.mean
    }
  if("med" %in% stats | "median" %in% stats) {
    d.agg.med <- aggregate(data[,var.cols], by = category.list, median, na.rm = TRUE)
    colnames(d.agg.med) <- c(cat.names, paste(vars, "med", sep = "_") )
    datalist$med <- d.agg.med
    }
  if("SD" %in% stats | "sd" %in% stats) {
    d.agg.SD <- aggregate(data[,var.cols], by = category.list, sd, na.rm = TRUE)
    colnames(d.agg.SD) <- c(cat.names, paste(vars, "SD", sep = "_") )
    datalist$SD <- d.agg.SD
    }
  if("VAR" %in% stats | "var" %in% stats) {
    d.agg.var <- aggregate(data[,var.cols], by = category.list, var, na.rm = TRUE)
    colnames(d.agg.var) <- c(cat.names, paste(vars, "var", sep = "_") )
    datalist$var <- d.agg.var
    }
  #
  new.dat <- Reduce(function(x,y) {merge(x,y)}, datalist)
  #
  # rearrange column sequence:
  n.cats  <- length(category.list)
  n.stats <- length(stats)
  n.vars  <- length(var.cols)
  col.seq <- 1:dim(new.dat)[2]
  for(i in 1:n.vars) {
    index <- (n.cats + 1) + ((i-1)*n.stats)
    col.seq[index:(index+(n.stats-1))] <- n.cats + ((i-1) + (1:n.stats * (n.vars) + (1 - n.vars)))
    }
  new.dat <- new.dat[ , col.seq]
  return(new.dat)
  }
#
#########################################
#
#
growing.season <- function(lay.date, sites, dat, threshold = 5) {
  # This function estimates the length of the growing season and growing degree days, given a set of monthly temperatures.
  # Season length and growing degree days are currently estimated by two methods: linear interpolation between monthly
  # means, and a spline interpolation fit to the 12 monthly means. In its current form, the function returns only the
  # spline interpolation estimates.
  # Arguments:
  #   lay.date is a vector with the Julian date of the start (egg laying, germination). This can be NA.
  #   sites is a character vector with site names.
  #   dat is a data.frame containing the estimated mean temperatures each month, jan-dec (12 columns, rows must be the same as the number of sites).
  #   threshold is the threshold temperature above which the organism can grow.
  #
  # This function is located in "2017-05 functions for climate analyses.R".
  #
  res <- data.frame(site = sites, lay.date = lay.date, first.date = NA, last.date = NA, temp.sum.season = 0, temp.sum.60d = NA, temp.60d = NA, season.length.spline = NA, temp.sum.season.spline = NA, strings.AsFactors = FALSE)
  for(i in 1:dim(res)[1] ) {
    prof <- as.data.frame(cbind( c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335), t(dat[i, ]) ))
    colnames(prof) <- c("date", "temp")
    for(j in 1:11) {   # Here begins the linear interpolation method
      if(prof$temp[j] < threshold & prof$temp[j+1] >= threshold) {
        res$first.date[i] <- prof$date[j] + ((prof$date[j+1] - prof$date[j]) * ((threshold - prof$temp[j])/(prof$temp[j+1] - prof$temp[j])))
        res$temp.sum.season[i] <- res$temp.sum.season[i] + ((prof$date[j+1] - res$first.date[i]) * (prof$temp[j+1]-threshold)) / 2
        }
      if(prof$temp[j] > threshold & prof$temp[j+1] > threshold) {
        res$temp.sum.season[i] <- res$temp.sum.season[i] + (prof$date[j+1] - prof$date[j]) * ((prof$temp[j+1]+prof$temp[j]) / 2)
        if(j == 1)  res$first.date[i] <- 0
        if(j == 11) res$last.date[i]  <- 365
        }
      if(prof$temp[j] >= threshold & prof$temp[j+1] < threshold) {
        res$last.date[i] <- prof$date[j] + ((prof$date[j+1] - prof$date[j]) * ((prof$temp[j] - threshold )/(prof$temp[j] - prof$temp[j+1])))
        res$temp.sum.season[i] <- res$temp.sum.season[i] + ((res$last.date[i] - prof$date[j]) * (prof$temp[j]-threshold)) / 2
        }
      }  # end j looping through months
    # Here begins the spline interpolation method
    spl <- smooth.spline(x = prof$date, prof$temp)
    if(!is.na(lay.date[i])) {
      day.prof <- data.frame(date = seq(lay.date[i]+14, lay.date[i]+73), temp = predict(object = spl, x = seq(lay.date[i]+14, lay.date[i]+73))$y)
      res$temp.sum.60d[i] <- sum(day.prof$temp - threshold)
      res$temp.60d[i] <- mean(day.prof$temp, na.rm = TRUE)
      }
    day.prof <- data.frame(date = seq(1,365), temp = predict(object = spl, x = seq(1,365))$y)
    day.prof <- day.prof[ day.prof$temp > threshold , ]
    res$season.length.spline[i]   <- dim(day.prof)[1]
    res$temp.sum.season.spline[i] <- sum(day.prof$temp - threshold)
    }    # end i looping through sites
  res$season.length <- res$last.date - res$first.date
  res <- res[ , c(1,2,8,9,7,6)]
  names(res)[c(3,4)] <- c("season.length", "temp.sum.season")
  res$site <- as.character(res$site)
  return(res)
  }
#
#########################################
#
#
modelII <- function(XjArray,YjArray){
# Calculate MODEL II Regression paramaters. Also called "reduced major axis
# regression" (Prentice 1987) or "geometric mean regression" (Webb et al. 1981).
# This gives same results as "SMA" (standard major axis) in package "lmodel2".
#
# XjArray and YjArray are 2wo one dimensional arrays XjArray and YjArray containing the X and Y vectors.
#
  sumXjYj <- 0
  sumXj <- 0
  sumYj <- 0
  n <- 0
  n <- length(XjArray)
  sumXjSquared <- 0
  sumYjSquared <- 0
  covariancexy <- 0
  for(i in 1:n){
    sumXjYj <- sumXjYj + XjArray[i] * YjArray[i]
    sumXj <- sumXj + XjArray[i]
    sumYj <- sumYj + YjArray[i]
    sumXjSquared <- sumXjSquared + XjArray[i]^2
    sumYjSquared <- sumYjSquared + YjArray[i]^2   }   
# Mean of X and Y vectors
  meanyj <- sumYj / n
  meanxj <- sumXj / n   
# Create covariance
  for(i in 1:n){ covariancexy <- covariancexy + ((XjArray[i] - meanxj) * (YjArray[i] - meanyj))   }
  covariancexy <- covariancexy / n      
# get variance of X and Y (SD)
  varXj <- (n * sumXjSquared-sumXj^2)/(n*(n - 1))
  varYj <- (n * sumYjSquared-sumYj^2)/(n*(n - 1))
  sdxij <- (sumXjSquared)-(sumXj^2/n)
  sdxik <- (sumYjSquared)-(sumYj^2/n)
# make beta 'sgn function to return sign with magnitude of 1
  betacoeff <- sign(covariancexy) * ((varYj^0.5) / (varXj^0.5))
# 'make intercept
  Intercept <- meanyj - meanxj * betacoeff   
# Make R the pearson produce moment correlation coefficient
  if (varYj==0 | varXj==0){
    corrCoeff <- 0
    } else {
    corrCoeff <- (sumXjYj - ((sumXj * sumYj) / n)) / ((sdxij * sdxik)^0.5)   }   
# Make sample variances of betacoefficient and intercept
  variancebeta <- (varYj / varXj) * ((1 - (corrCoeff ^ 2)) / n)
  varianceintercept <- (varYj / n) * (1 - corrCoeff) * (2 + ((meanxj ^ 2) * ((1 + corrCoeff) / varXj)))
  sdbeta <- variancebeta^0.5
  sdintercept <- varianceintercept^0.5
  list(betacoeff=betacoeff,
       Intercept=Intercept,
       sdbeta=sdbeta,
       sdintercept=sdintercept,
       meanxj=meanxj,
       meanyj=meanyj,
       corrCoeff=corrCoeff)
  }
#
#########################################
#
#
add.short.line <- function(a, b, x.vals, color = "black", lwidth = 1, ltype = 1, x.log = FALSE) {
  #
  # This function adds regression lines that extend only as far as the data
  # The plot should already have been created.
  # Arguments:
  #  a and b are the intercept and slope of the regression.
  #  x.vals is a vector containing the x values of the data. Line will go from min(x.vals) to max(x.vals). NAs are allowed.
  #  Other arguments are as for normal plotting.
  #
  range <- c(min(x.vals, na.rm=TRUE), max(x.vals, na.rm=TRUE))
  newdat <- data.frame(predict.y=NA, x=seq(range[1], range[2], length.out=length(x.vals) ) )
  if(x.log == FALSE) newdat$predict.y <- a + newdat$x * b
  if(x.log == TRUE)  newdat$predict.y <- a + log(newdat$x) * b
  lines(newdat$x, newdat$predict.y, col=color, lwd = lwidth, lty = ltype)
  }
#
#########################################
