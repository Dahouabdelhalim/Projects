## Marta Shocket, Stanford University / UCLA, marta.shocket@gmail.com
## Started Jan 2018, Updated August 2020
##
## Purpose: Plot trait thermal responses for many Culex and Aedes species
##


# Set working directory
setwd("~/Fitting Traits")


##########
###### Make list of trait averages
##########

TraitAvgList <- list(Cpip.pLA.avg, Cqui.pLA.avg, Ctar.pLA.avg,
                     Cpip.MDR.avg, Cqui.MDR.avg, Ctar.MDR.avg,
                     Cpip.a.avg, Cqui.a.avg, data.a.Ctar,
                     data.EFGC.Cpip, Cqui.ER.avg,
                     Cpip.lf.avg, Cqui.lf.avg, Ctar.lf.avg,
                     CpipWNV.PDR.avg, CtarWEEV.PDR.avg,
                     CpipWNV.c.avg, CpipWNV.bc.avg, CtarWNV.b.avg, CuniWNV.bc.avg, 
                     CtarWEEV.c.avg, CtarWEEV.b.avg, CtarWEEV.bc.avg, CtarSLEV.c.avg, CtarSLEV.b.avg, 
                     CpipSINV.c.avg, AtaeSINV.c.avg, AtaeRVFV.bc.avg, AtriEEEV.bc.avg,
                     Atri.pLA.avg, Atae.lf.avg,
                     Cpip.pO.avg, Cqui.pO.avg,
                     Cthe.EV.avg, Avex.EV.avg)
save(TraitAvgList, file = "TraitAvgList.RSave")


##########
###### Figures of mosquito traits
##########

# Load/subset data, saved fits, and parameters for plotting

Temp.xs <- seq(1, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

# Load averaged trait data
load("TraitAvgList.RSave")
Cpip.pLA.avg <- TraitAvgList[[1]]
Cqui.pLA.avg <- TraitAvgList[[2]]
Ctar.pLA.avg <- TraitAvgList[[3]]
Cpip.MDR.avg <- TraitAvgList[[4]]
Cqui.MDR.avg <- TraitAvgList[[5]]
Ctar.MDR.avg <- TraitAvgList[[6]]
Cpip.a.avg <- TraitAvgList[[7]]
Cqui.a.avg <- TraitAvgList[[8]]
data.a.Ctar <- TraitAvgList[[9]]
data.EFOC.Cpip <- TraitAvgList[[10]] 
Cqui.EPR.avg <- TraitAvgList[[11]]
Cpip.lf.avg <- TraitAvgList[[12]]
Cqui.lf.avg <- TraitAvgList[[13]]
Ctar.lf.avg <- TraitAvgList[[14]]

CpipWNV.PDR.avg <- TraitAvgList[[15]]
CtarWEEV.PDR.avg <- TraitAvgList[[16]]
CpipWNV.c.avg <- TraitAvgList[[17]]
CpipWNV.bc.avg <- TraitAvgList[[18]]
CtarWNV.b.avg <- TraitAvgList[[19]]
CuniWNV.bc.avg <- TraitAvgList[[20]]
CtarWEEV.c.avg <- TraitAvgList[[21]]
CtarWEEV.b.avg <- TraitAvgList[[22]]
CtarWEEV.bc.avg <- TraitAvgList[[23]]
CtarSLEV.c.avg <- TraitAvgList[[24]]
CtarSLEV.b.avg <- TraitAvgList[[25]]
CpipSINV.c.avg <- TraitAvgList[[26]]
AtaeSINV.c.avg <- TraitAvgList[[27]]
AtaeRVFV.bc.avg <- TraitAvgList[[28]]
AtriEEEV.bc.avg <- TraitAvgList[[29]]

Atri.pLA.avg <- TraitAvgList[[30]]
Atae.lf.avg <- TraitAvgList[[31]]

Cpip.pO.avg <- TraitAvgList[[32]]
Cqui.pO.avg <- TraitAvgList[[33]]

Cthe.EV.avg <- TraitAvgList[[34]]
Avex.EV.avg <- TraitAvgList[[35]]

# jitter overlapping pO points
Cpip.pO.avg$Temp[2] <- 20.25 # is actually 20
Cqui.pO.avg$Temp[4] <- 19.75 # is actually 20

Cpip.pO.avg$Temp[1] <- 24.75 # is actually 25
Cqui.pO.avg$Temp[2] <- 25.25 # is actually 25

# Load raw data for traits that without averages (because no replicate data at same temps)
data.MDR <- read.csv("TraitData_MDR.csv")
data.MDR.Atri <- subset(data.MDR, host.code == "Atri")
data.MDR.Avex <- subset(data.MDR, host.code == "Avex")

data.EV <- read.csv("TraitData_EV.csv")
data.EV.Cpip <- subset(data.EV, host.code == "Cpip")
data.EV.Cqui <- subset(data.EV, host.code == "Cqui") 
data.EV.Cqui <- data.EV.Cqui[4:22,] # remove Oda data

data.pLA <- read.csv("TraitData_pLA.csv") 
data.pLA.Atri <- subset(data.pLA, host.code == "Atri")
data.pLA.Avex <- subset(data.pLA, host.code == "Avex")

data.Cmel <- read.csv("Cmeldata.csv")
data.pLA.Cmel <- subset(data.Cmel, trait.name == "pLA")
data.MDR.Cmel <- subset(data.Cmel, trait.name == "1/MDR")
data.pO.Cmel <- subset(data.Cmel, trait.name == "pO")
data.a.Cmel <- subset(data.Cmel, trait.name == "GCD")

# Load saved thermal responses
load("jagsout_a_Cpip_inf.Rdata")
load("jagsout_a_Cqui_inf.Rdata")
load("jagsout_a_Ctar_inf.Rdata")

load("jagsout_pLA_Cpip_inf.Rdata")
load("jagsout_pLA_Cqui_inf.Rdata")
load("jagsout_pLA_Ctar_inf.Rdata")

load("jagsout_MDR_Cpip_inf.Rdata")
load("jagsout_MDR_Cqui_inf.Rdata")
load("jagsout_MDR_Ctar_inf.Rdata")

load("jagsout_EFOC_Cpip_inf.Rdata")
load("jagsout_EPR_Cqui_inf.Rdata")
load("jagsout_pO_Cpip_inf.Rdata")
load("jagsout_pO_Cqui_inf.Rdata")

load("jagsout_EV_Cpip_inf.Rdata")
load("jagsout_EV_Cqui_inf.Rdata")

load("jagsout_lf_Cpip_inf.Rdata")
load("jagsout_lf_Cqui_inf.Rdata")
load("jagsout_lf_Ctar_inf.Rdata")

load("jagsout_a_Cmel_inf.Rdata")
load("jagsout_pLA_Atri_inf.Rdata")
load("jagsout_pLA_Avex_inf.Rdata")
load("jagsout_pLA_Cmel_inf.Rdata")
load("jagsout_MDR_Atri_inf.Rdata")
load("jagsout_MDR_Avex_inf.Rdata")
load("jagsout_MDR_Cmel_inf.Rdata")
load("jagsout_pO_Cmel_inf.Rdata")
load("jagsout_EV_Cthe_inf.Rdata")
load("jagsout_EV_Avex_inf.Rdata")
load("jagsout_lf_Atae_inf.Rdata")


######### Make polygon coordinates
cord.x <- c(Temp.xs, rev(Temp.xs))

Cpip.a.cord.y <- c(a.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(a.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Cqui.a.cord.y <- c(a.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(a.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Ctar.a.cord.y <- c(a.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(a.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Cmel.a.cord.y <- c(a.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(a.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

Cpip.pLA.cord.y <- c(pLA.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(pLA.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Cqui.pLA.cord.y <- c(pLA.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(pLA.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Ctar.pLA.cord.y <- c(pLA.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(pLA.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

Cpip.MDR.cord.y <- c(MDR.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(MDR.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Cqui.MDR.cord.y <- c(MDR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(MDR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Ctar.MDR.cord.y <- c(MDR.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(MDR.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

Cpip.EFOC.cord.y <- c(EFOC.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(EFOC.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Cqui.EPR.cord.y <- c(EPR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(EPR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Cpip.pO.cord.y <- c(pO.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(pO.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Cqui.pO.cord.y <- c(pO.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(pO.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

Cpip.EV.cord.y <- c(EV.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(EV.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Cqui.EV.cord.y <- c(EV.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(EV.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

Atri.pLA.cord.y <- c(pLA.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(pLA.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Avex.pLA.cord.y <- c(pLA.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(pLA.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Cmel.pLA.cord.y <- c(pLA.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(pLA.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Atri.MDR.cord.y <- c(MDR.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(MDR.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Avex.MDR.cord.y <- c(MDR.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(MDR.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Cmel.MDR.cord.y <- c(MDR.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(MDR.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Cmel.pO.cord.y <- c(pO.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(pO.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Cthe.EV.cord.y <- c(EV.Cthe.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(EV.Cthe.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Avex.EV.cord.y <- c(EV.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(EV.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

Cpip.lf.top <- as.numeric(lf.Cpip.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"])
Cpip.lf.top[1:140] <- as.numeric(Cpip.lf.top[141])
Cpip.lf.bot <- rev(as.numeric(lf.Cpip.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"]))
Cpip.lf.bot[441:301] <- as.numeric(Cpip.lf.bot[300])

Cqui.lf.top <- as.numeric(lf.Cqui.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"])
Cqui.lf.top[1:150] <- as.numeric(Cqui.lf.top[151])
Cqui.lf.bot <- rev(as.numeric(lf.Cqui.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"]))
Cqui.lf.bot[441:291] <- as.numeric(Cqui.lf.bot[290])

Ctar.lf.top <- as.numeric(lf.Ctar.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"])
Ctar.lf.top[1:130] <- as.numeric(Ctar.lf.top[131])
Ctar.lf.bot <- rev(as.numeric(lf.Ctar.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"]))
Ctar.lf.bot[441:311] <- as.numeric(Ctar.lf.bot[310])

Atae.lf.top <- as.numeric(lf.Atae.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"])
Atae.lf.top[1:230] <- as.numeric(Atae.lf.top[231])
Atae.lf.bot <- rev(as.numeric(lf.Atae.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"]))
Atae.lf.bot[441:231] <- as.numeric(Atae.lf.bot[230])

Cpip.lf.cord.y <- c(Cpip.lf.top, Cpip.lf.bot)
Cqui.lf.cord.y <- c(Cqui.lf.top, Cqui.lf.bot)
Ctar.lf.cord.y <- c(Ctar.lf.top, Ctar.lf.bot)
Atae.lf.cord.y <- c(Atae.lf.top, Atae.lf.bot)

############################## Figure 3: MDR, pLA, a, and lf in 3 main Cx. spp.
# Set plot settings
par(mfrow = c(2,2), mar = c(1.5, 4.5, 2, 1), oma = c(4, 0, 0, 0))

######### A: MDR
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,0.17), data = Cpip.MDR.avg, xaxt = "n", pch = 19, main = expression(paste("Mosquito Development Rate (",italic(MDR),")")),
     ylab = expression(paste("Rate (day"^-1,")")), xlab = "", cex.lab = 1.15, col = "white")
axis(side = 1, labels = FALSE)
polygon(cord.x, Cpip.MDR.cord.y, col = rgb(0.25, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, Cqui.MDR.cord.y, col = rgb(0.9, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, Ctar.MDR.cord.y, col = rgb(0.25, 0.25, 0.9, 0.2), border = NA)
lines(MDR.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "grey25", lwd = 1.5)
lines(MDR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "red2", lwd = 1.5)
lines(MDR.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "blue", lwd = 1.5)
arrows(Cpip.MDR.avg$Temp, Cpip.MDR.avg$upper, Cpip.MDR.avg$Temp, Cpip.MDR.avg$lower, length = 0, lwd = 1.5, col = "grey25")
arrows(Cqui.MDR.avg$Temp, Cqui.MDR.avg$upper, Cqui.MDR.avg$Temp, Cqui.MDR.avg$lower, length = 0, lwd = 1.5, col = "red2")
arrows(Ctar.MDR.avg$Temp, Ctar.MDR.avg$upper, Ctar.MDR.avg$Temp, Ctar.MDR.avg$lower, length = 0, lwd = 1.5, col = "blue")
points(Cpip.MDR.avg$mean ~ Cpip.MDR.avg$Temp, pch = 19, col = "grey25", cex = 1.25)
points(Cqui.MDR.avg$mean ~ Cqui.MDR.avg$Temp, pch = 19, col = "red2", cex = 1.25)
points(Ctar.MDR.avg$mean ~ Ctar.MDR.avg$Temp, pch = 19, col = "blue", cex = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.1)

a <- expression(paste(italic(Cx.)," ",italic(pipiens)))
b <- expression(paste(italic(Cx.)," ",italic(quinque.)))
c <- expression(paste(italic(Cx.)," ",italic(tarsalis)))
legend(x = 6, y = 0.1675, bty = "n", legend = c(a, b, c), pch = 19, col = c("grey25", "red2", "blue"), cex = 1.1, lwd = 1.5)
legend(x = 3, y = 0.1675, bty = "n", legend = c("", "", ""), border = c("white", "white", "white"), cex = 1.1, 
       fill = c(rgb(0.25, 0.25, 0.25, 0.2), rgb(0.9, 0.25, 0.25, 0.3), rgb(0.25, 0.25, 0.9, 0.3)))

######### B: PLA
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1.05), data = Cpip.pLA.avg, xaxt = "n", pch = 19, main = expression(paste("Larval-to-Adult Survival (",italic(pLA),")")),
     ylab = "Survival probability", xlab = "", cex.lab = 1.15, col = "white")
axis(side = 1, labels = FALSE)
polygon(cord.x, Cpip.pLA.cord.y, col = rgb(0.25, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, Cqui.pLA.cord.y, col = rgb(0.9, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, Ctar.pLA.cord.y, col = rgb(0.25, 0.25, 0.9, 0.2), border = NA)
lines(pLA.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "grey25", lwd = 1.5)
lines(pLA.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "red2", lwd = 1.5)
lines(pLA.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "blue", lwd = 1.5)
arrows(Cpip.pLA.avg$Temp, Cpip.pLA.avg$upper, Cpip.pLA.avg$Temp, Cpip.pLA.avg$lower, length = 0, lwd = 1.5, col = "grey25")
arrows(Cqui.pLA.avg$Temp, Cqui.pLA.avg$upper, Cqui.pLA.avg$Temp, Cqui.pLA.avg$lower, length = 0, lwd = 1.5, col = "red2")
arrows(Ctar.pLA.avg$Temp, Ctar.pLA.avg$upper, Ctar.pLA.avg$Temp, Ctar.pLA.avg$lower, length = 0, lwd = 1.5, col = "blue")
points(Cpip.pLA.avg$mean ~ Cpip.pLA.avg$Temp, pch = 19, col = "grey25", cex = 1.25)
points(Cqui.pLA.avg$mean ~ Cqui.pLA.avg$Temp, pch = 19, col = "red2", cex = 1.25)
points(Ctar.pLA.avg$mean ~ Ctar.pLA.avg$Temp, pch = 19, col = "blue", cex = 1.25)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.1)

######### C: a
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,0.5), data = Cpip.a.avg, pch = 19, main = expression(paste("Biting Rate (",italic(a),")")),
     ylab = expression(paste("Rate (day"^-1,")")), xlab = "", cex.lab = 1.15, col = "white")
polygon(cord.x, Cpip.a.cord.y, col = rgb(0.25, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, Cqui.a.cord.y, col = rgb(0.9, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, Ctar.a.cord.y, col = rgb(0.25, 0.25, 0.9, 0.2), border = NA)
lines(a.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "grey25", lwd = 1.5)
lines(a.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "red2", lwd = 1.5)
lines(a.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "blue", lwd = 1.5)
arrows(Cpip.a.avg$Temp, Cpip.a.avg$upper, Cpip.a.avg$Temp, Cpip.a.avg$lower, length = 0, lwd = 1.5, col = "grey25") # only pip has error bars
points(Cpip.a.avg$mean ~ Cpip.a.avg$Temp, pch = 19, col = "grey25", cex = 1.25)
points(Cqui.a.avg$mean ~ Cqui.a.avg$Temp, pch = 19, col = "red2", cex = 1.25)
points(1/data.a.Ctar$trait ~ data.a.Ctar$T, pch = 19, col = "blue", cex = 1.25)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.1)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

######### D: lf
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,115), data = Cpip.lf.avg, pch = 19, main = expression(paste("Adult Lifespan (",italic(lf),")")),
     ylab = expression(paste("Time (days)")), xlab = "", cex.lab = 1.15, col = "white")
polygon(cord.x, Cpip.lf.cord.y, col = rgb(0.25, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, Cqui.lf.cord.y, col = rgb(0.9, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, Ctar.lf.cord.y, col = rgb(0.25, 0.25, 0.9, 0.2), border = NA)
lines(lf.Cpip.out.inf$BUGSoutput$summary[141:441, "mean"] ~ Temp.xs[141:441], lty = 1, col = "grey25", lwd = 1.5)
lines(lf.Cqui.out.inf$BUGSoutput$summary[151:441, "mean"] ~ Temp.xs[151:441], lty = 1, col = "red2", lwd = 1.5)
lines(lf.Ctar.out.inf$BUGSoutput$summary[131:441, "mean"] ~ Temp.xs[131:441], lty = 1, col = "blue", lwd = 1.5)
lines(rep(lf.Cpip.out.inf$BUGSoutput$summary[141, "mean"],141) ~ Temp.xs[1:141], lty = 1, col = "grey25", lwd = 1.5)
lines(rep(lf.Cqui.out.inf$BUGSoutput$summary[151, "mean"],151) ~ Temp.xs[1:151], lty = 1, col = "red2", lwd = 1.5)
lines(rep(lf.Ctar.out.inf$BUGSoutput$summary[131, "mean"],131) ~ Temp.xs[1:131], lty = 1, col = "blue", lwd = 1.5)
arrows(Cpip.lf.avg$Temp, Cpip.lf.avg$upper, Cpip.lf.avg$Temp, Cpip.lf.avg$lower, length = 0, lwd = 1.5, col = "grey25")
arrows(Cqui.lf.avg$Temp, Cqui.lf.avg$upper, Cqui.lf.avg$Temp, Cqui.lf.avg$lower, length = 0, lwd = 1.5, col = "red2")
arrows(Ctar.lf.avg$Temp, Ctar.lf.avg$upper, Ctar.lf.avg$Temp, Ctar.lf.avg$lower, length = 0, lwd = 1.5, col = "blue")
points(Cpip.lf.avg$mean ~ Cpip.lf.avg$Temp, pch = 19, col = "grey25", cex = 1.25)
points(Cqui.lf.avg$mean ~ Cqui.lf.avg$Temp, pch = 19, col = "red2", cex = 1.25)
points(Ctar.lf.avg$mean ~ Ctar.lf.avg$Temp, pch = 19, col = "blue", cex = 1.25)
legend("topleft", legend = "D", bty = "n", adj = 1, cex = 1.1)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)


############################## Figure 4: Fecundity Traits in 3 main Cx. spp.

# Set plot settings
par(mfrow = c(3,1), mar = c(1.5, 4.5, 2, 1), oma = c(4, 0, 0, 0))

######### A: pO
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = Cpip.pO.avg, xaxt = "n", pch = 19, main = expression(paste("Proportion Ovipositing (",italic(pO),")")),
     ylab = "Proportion ovipositing", xlab = "", cex.lab = 1.2, col = "white", cex.axis = 1.1)
axis(side = 1, labels = FALSE)
polygon(cord.x, Cpip.pO.cord.y, col = rgb(0.25, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, Cqui.pO.cord.y, col = rgb(0.9, 0.25, 0.25, 0.2), border = NA)
lines(pO.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "grey25", lwd = 1.5)
lines(pO.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "red2", lwd = 1.5)
arrows(Cpip.pO.avg$Temp, Cpip.pO.avg$upper, Cpip.pO.avg$Temp, Cpip.pO.avg$lower, length = 0, lwd = 1.5, col = "grey25")
arrows(Cqui.pO.avg$Temp, Cqui.pO.avg$upper, Cqui.pO.avg$Temp, Cqui.pO.avg$lower, length = 0, lwd = 1.5, col = "red2")
points(Cpip.pO.avg$mean ~ Cpip.pO.avg$Temp, pch = 19, col = "grey25", cex = 1.25)
points(Cqui.pO.avg$mean ~ Cqui.pO.avg$Temp, pch = 19, col = "red2", cex = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

a <- expression(paste(italic(Cx.)," ",italic(pip.)))
b <- expression(paste(italic(Cx.)," ",italic(qui.)))
legend(x = 33, y = 0.9, bty = "n", legend = c(a, b), pch = 19, col = c("grey25", "red2"), cex = 1.1, lwd = 1.5)
legend(x = 30, y = 0.9, bty = "n", legend = c("", ""), border = c("white", "white"), cex = 1.1, 
       fill = c(rgb(0.25, 0.25, 0.25, 0.2), rgb(0.9, 0.25, 0.25, 0.3)))

######### B: EFD
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,320), data = Cqui.EPR.avg, xaxt = "n", pch = 19, main = expression(paste("Fecundity (",italic(EFGC),", ",italic(ER),")")),
     ylab = expression(paste("Eggs (see caption for units)")), xlab = "", cex.lab = 1.2, col = "white" , cex.axis = 1.1)
axis(side = 1, labels = FALSE)
polygon(cord.x, Cpip.EFOC.cord.y, col = rgb(0.25, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, Cqui.EPR.cord.y, col = rgb(0.9, 0.25, 0.25, 0.2), border = NA)
lines(EFOC.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "grey25", lwd = 1.5)
lines(EPR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "red2", lwd = 1.5)
arrows(Cqui.EPR.avg$Temp, Cqui.EPR.avg$upper, Cqui.EPR.avg$Temp, Cqui.EPR.avg$lower, length = 0, lwd = 1.5, col = "red2")
points(Cqui.EPR.avg$mean ~ Cqui.EPR.avg$Temp, pch = 19, col = "red2", cex = 1.25)
points(data.EFOC.Cpip$trait ~ data.EFOC.Cpip$T, pch = 19, col = "grey25", cex = 1.25)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.2)

######### C: EV
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.EV.Cpip, pch = 19, main = expression(paste("Egg Viability (",italic(EV),")")),
     ylab = "Proportion hatching", xlab = "", cex.lab = 1.2, col = "white", cex.axis = 1.1)
polygon(cord.x, Cpip.EV.cord.y, col = rgb(0.25, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, Cqui.EV.cord.y, col = rgb(0.9, 0.25, 0.25, 0.2), border = NA)
lines(EV.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "grey25", lwd = 1.5)
lines(EV.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "red2", lwd = 1.5)
arrows(Cqui.EV.avg$Temp, Cqui.EV.avg$upper, Cqui.EV.avg$Temp, Cqui.EV.avg$lower, length = 0, lwd = 1.5, col = "red2")
points(data.EV.Cqui$trait ~ data.EV.Cqui$T, pch = 19, col = "red2", cex = 1.25)
points(data.EV.Cpip$trait ~ data.EV.Cpip$T, pch = 19, col = "grey25", cex = 1.25)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)


######## Supplemental Figure for mosquito traits in other vector spp
par(mfrow = c(3,2), mar = c(1, 4.5, 2, 1), oma = c(4, 0, 0, 0), las = 1)

######### A: MDR
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,0.18), data = data.MDR.Atri, xaxt = "n", pch = 19, main = expression(paste("Mosquito Development Rate (",italic(MDR),")")),
     ylab = expression(paste("Rate (day"^-1,")")), xlab = "", cex.lab = 1.2, col = "white")
axis(side = 1, labels = FALSE)
polygon(cord.x, Atri.MDR.cord.y, col = rgb(0.5, 0.1, 0.8, 0.3), border = NA)
polygon(cord.x, Avex.MDR.cord.y, col = rgb(0.2, 0.5, 0.5, 0.3), border = NA)
polygon(cord.x, Cmel.MDR.cord.y, col = rgb(0.7, 0.3, 0.1, 0.25), border = NA)
lines(MDR.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "darkorchid", lwd = 1.5)
lines(MDR.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "cyan4", lwd = 1.5)
lines(MDR.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "saddlebrown", lwd = 1.5)
points(1/data.MDR.Atri$trait ~ data.MDR.Atri$T, pch = 19, col = "darkorchid", cex = 1.25)
points(1/data.MDR.Avex$trait ~ data.MDR.Avex$T, pch = 19, col = "cyan4", cex = 1.25)
points(1/data.MDR.Cmel$trait ~ data.MDR.Cmel$T, pch = 19, col = "saddlebrown", cex = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

a <- expression(paste(italic(Ae.)," ",italic(tri.)))
b <- expression(paste(italic(Ae.)," ",italic(vex.)))
c <- expression(paste(italic(Cs.)," ",italic(mel.)))
legend(x = 9, y = 0.19, bty = "n", legend = c(a, b, c), pch = 19, col = c("darkorchid", "cyan4", "saddlebrown"), cex = 1.3, lwd = 1.5)
legend(x = 6, y = 0.19, bty = "n", legend = c("", "", ""), fill = c(col = rgb(0.5, 0.1, 0.8, 0.2), rgb(0.2, 0.5, 0.5, 0.3), rgb(0.7, 0.3, 0.1, 0.3)), border = c("white", "white", "white"), cex = 1.3)

######### B: PLA
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1.05), data = Atri.pLA.avg, xaxt = "n", pch = 19, main = expression(paste("Larval-to-Adult Survival (",italic(pLA),")")),
     ylab = "Survival probability", xlab = "", cex.lab = 1.2, col = "white")
axis(side = 1, labels = FALSE)
polygon(cord.x, Atri.pLA.cord.y, col = rgb(0.5, 0.1, 0.8, 0.3), border = NA)
polygon(cord.x, Avex.pLA.cord.y, col = rgb(0.2, 0.5, 0.5, 0.4), border = NA)
polygon(cord.x, Cmel.pLA.cord.y, col = rgb(0.7, 0.3, 0.1, 0.25), border = NA)
lines(pLA.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "darkorchid", lwd = 1.5)
lines(pLA.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "cyan4", lwd = 1.5)
lines(pLA.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "saddlebrown", lwd = 1.5)
arrows(Atri.pLA.avg$Temp, Atri.pLA.avg$upper, Atri.pLA.avg$Temp, Atri.pLA.avg$lower, length = 0, lwd = 1.5, col = "darkorchid")
points(Atri.pLA.avg$mean ~ Atri.pLA.avg$Temp, pch = 19, col = "darkorchid", cex = 1.25)
points(data.pLA.Avex$trait ~ data.pLA.Avex$T, pch = 19, col = "cyan4", cex = 1.25)
points(data.pLA.Cmel$trait ~ data.pLA.Cmel$T, pch = 19, col = "saddlebrown", cex = 1.25)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.2)

#legend(x = 30.5, y = 1.12, bty = "n", legend = c(a, b, c), pch = 19, col = c("darkorchid", "cyan4", "saddlebrown"), cex = 1.3, lwd = 1.5)
#legend(x = 27.5, y = 1.12, bty = "n", legend = c("", "", ""), fill = c(col = rgb(0.5, 0.1, 0.8, 0.2), rgb(0.2, 0.5, 0.5, 0.3), rgb(0.7, 0.3, 0.1, 0.3)), border = c("white", "white", "white"), cex = 1.3)


######### C: Biting Rate
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,0.25), data = data.a.Cmel, xaxt = "n", pch = 19, main = expression(paste("Biting Rate (",italic(a),")")),
     ylab = expression(paste("Rate (day"^-1,")")), xlab = "", cex.lab = 1.2, col = "white", cex.axis = 1.1)
axis(side = 1, labels = FALSE)
polygon(cord.x, Cmel.a.cord.y, col = rgb(0.7, 0.3, 0.1, 0.25), border = NA)
lines(a.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "saddlebrown", lwd = 1.5)
points(1/data.a.Cmel$trait ~ data.a.Cmel$T, pch = 19, col = "saddlebrown", cex = 1.25)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.2)
d <- expression(paste(italic(Cs.)," ",italic(melanura)))
text(x = 38, y = 0.2, labels = d, cex = 1.3)

######### D: lf
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,100), data = Atae.lf.avg, xaxt = "n", pch = 19, main = expression(paste("Adult Lifespan (",italic(lf),")")),
     ylab = expression(paste("Time (days)")), xlab = "", cex.lab = 1.2, col = "white")
axis(side = 1, labels = FALSE)
polygon(cord.x, Atae.lf.cord.y, col = rgb(0.1, 0.9, 0.25, 0.2), border = NA)
lines(lf.Atae.out.inf$BUGSoutput$summary[211:441, "mean"] ~ Temp.xs[211:441], lty = 1, col = "mediumseagreen", lwd = 1.5)
lines(rep(lf.Atae.out.inf$BUGSoutput$summary[211, "mean"],210) ~ Temp.xs[1:210], lty = 1, col = "mediumseagreen", lwd = 1.5)
arrows(Atae.lf.avg$Temp, Atae.lf.avg$upper, Atae.lf.avg$Temp, Atae.lf.avg$lower, length = 0, lwd = 1.5, col = "mediumseagreen")
points(Atae.lf.avg$mean ~ Atae.lf.avg$Temp, pch = 19, col = "mediumseagreen", cex = 1.25)
legend("topleft", legend = "D", bty = "n", adj = 1, cex = 1.2)
e <- expression(paste(italic(Ae.)," ",italic(taeniorhynchus)))
text(x = 17, y = 55, labels = e, cex = 1.3)

######### E: pO
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pO.Cmel, pch = 19, main = expression(paste("Proportion Ovipositing (",italic(pO),")")),
     ylab = "Proportion ovipositing", xlab = "", cex.lab = 1.2, col = "white", cex.axis = 1.1)
polygon(cord.x, Cmel.pO.cord.y, col = rgb(0.7, 0.3, 0.1, 0.25), border = NA)
lines(pO.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "saddlebrown", lwd = 1.5)
points(data.pO.Cmel$trait ~ data.pO.Cmel$T, pch = 19, col = "saddlebrown", cex = 1.25)
legend("topleft", legend = "E", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)
text(x = 37, y = 0.8, labels = d, cex = 1.3)

######### F: EV
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.EV.Cthe, pch = 19, main = expression(paste("Egg Viability (",italic(EV),")")),
     ylab = "Proportion hatching", xlab = "", cex.lab = 1.2, col = "white", cex.axis = 1.1)
polygon(cord.x, Cthe.EV.cord.y, col = rgb(0.9, 0.4, 0.5, 0.3), border = NA)
polygon(cord.x, Avex.EV.cord.y, col = rgb(0.2, 0.5, 0.5, 0.3), border = NA)
lines(EV.Cthe.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "hotpink", lwd = 1.5)
lines(EV.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "cyan4", lwd = 1.5)
arrows(Cthe.EV.avg$Temp, Cthe.EV.avg$upper, Cthe.EV.avg$Temp, Cthe.EV.avg$lower, length = 0, lwd = 1.5, col = "hotpink")
arrows(Avex.EV.avg$Temp, Avex.EV.avg$upper, Avex.EV.avg$Temp, Avex.EV.avg$lower, length = 0, lwd = 1.5, col = "cyan4")
points(Cthe.EV.avg$mean ~ Cthe.EV.avg$Temp, pch = 19, col = "hotpink", cex = 1.25)
points(Avex.EV.avg$mean ~ Avex.EV.avg$Temp, pch = 19, col = "cyan4", cex = 1.25)
legend("topleft", legend = "F", bty = "n", adj = 1, cex = 1.3)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

b <- expression(paste(italic(Ae.)," ",italic(vexans)))
d <- expression(paste(italic(Cx.)," ",italic(theileri)))
legend(x = 15, y = 0.4, bty = "n", legend = c(b, d), pch = 19, col = c("cyan4", "hotpink"), cex = 1.3, lwd = 1.5)
legend(x = 12, y = 0.4, bty = "n", legend = c("", ""), fill = c(col = rgb(0.2, 0.5, 0.5, 0.3), rgb(0.7, 0.3, 0.1, 0.3)), border = c("white", "white"), cex = 1.3)

#polygon(c(25, 30, 30, 25), c(0, 0, .2, .2), col = rgb(0.9, 0.4, 0.5, 0.3), border = NA)


##########
###### Figure 5: PDR
##########

data.PDR <- read.csv("TraitData_PDR.csv")

# Subset Data
data.PDR.CquiWNV <- subset(data.PDR, joint.code == "CquiWNV")
data.PDR.CtarWNV <- subset(data.PDR, joint.code == "CtarWNV")
data.PDR.CuniWNV <- subset(data.PDR, joint.code == "CuniWNV")
data.PDR.CtarSLEV <- subset(data.PDR, joint.code == "CtarSLEV")
data.PDR.AtaeSINV <- subset(data.PDR, joint.code == "AtaeSINV")
data.PDR.AtaeRVFV <- subset(data.PDR, joint.code == "AtaeRVFV")
data.PDR.AtriEEEV <- subset(data.PDR, joint.code == "AtriEEEV")

# Load saved thermal responses
load("jagsout_PDR_CpipWNV_inf.Rdata")
load("jagsout_PDR_CquiWNV_inf.Rdata")
load("jagsout_PDR_CtarWNV_inf.Rdata")
load("jagsout_PDR_CuniWNV_inf.Rdata")
load("jagsout_PDR_CtarWEEV_inf.Rdata")
load("jagsout_PDR_CtarSLEV_inf.Rdata")
load("jagsout_PDR_AtriEEEV_inf.Rdata")
load("jagsout_PDR_AtaeRVFV_inf.Rdata")

######### Make polygon coordinates
cord.x <- c(Temp.xs, rev(Temp.xs))

CpipWNV.PDR.cord.y <- c(PDR.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(PDR.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
CquiWNV.PDR.cord.y <- c(PDR.CquiWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(PDR.CquiWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
CtarWNV.PDR.cord.y <- c(PDR.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(PDR.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
CuniWNV.PDR.cord.y <- c(PDR.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(PDR.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

CtarWEEV.PDR.cord.y <- c(PDR.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(PDR.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
CtarSLEV.PDR.cord.y <- c(PDR.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(PDR.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

AtriEEEV.PDR.cord.y <- c(PDR.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(PDR.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
AtaeRVFV.PDR.cord.y <- c(PDR.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(PDR.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))


par(mfrow = c(3,1), mar = c(1, 4.5, 2, 1), oma = c(4, 0, 0, 0), las = 1)

######### A: PDR for WNV
# "Pathogen Development Rate (",italic(PDR),") - ",
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,0.305), data = CpipWNV.PDR.avg, xaxt = "n", pch = 19, main = expression(paste("WNV")),
     ylab = expression(paste("Rate (day"^-1,")")), xlab = "", cex.lab = 1.15, col = "white", cex.lab = 1.2)
axis(side = 1, labels = FALSE)
polygon(cord.x, CuniWNV.PDR.cord.y, col = rgb(1, 0.8, 0, 0.3), border = NA)
polygon(cord.x, CpipWNV.PDR.cord.y, col = rgb(0.25, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, CquiWNV.PDR.cord.y, col = rgb(0.9, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, CtarWNV.PDR.cord.y, col = rgb(0.25, 0.25, 0.9, 0.2), border = NA)
lines(PDR.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "grey25", lwd = 1.5)
lines(PDR.CquiWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "red2", lwd = 1.5)
lines(PDR.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "blue", lwd = 1.5)
lines(PDR.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "darkorange", lwd = 1.5)
arrows(CpipWNV.PDR.avg$Temp,CpipWNV.PDR.avg$upper, CpipWNV.PDR.avg$Temp, CpipWNV.PDR.avg$lower, length = 0, lwd = 1.5, col = "grey25")
points(CpipWNV.PDR.avg$mean ~ CpipWNV.PDR.avg$Temp, pch = 19, col = "grey25", cex = 1.25)
points(1/trait ~ T, data.PDR.CquiWNV, pch = 19, col = "red2", cex = 1.25)
points(1/trait ~ T, data.PDR.CtarWNV, pch = 19, col = "blue", cex = 1.25)
points(1/trait ~ T, data.PDR.CuniWNV, pch = 19, col = "darkorange", cex = 1.25)
a <- expression(paste(italic(Cx.)," ",italic(pipiens)))
b <- expression(paste(italic(Cx.)," ",italic(quinquefasciatus)))
c <- expression(paste(italic(Cx.)," ",italic(tarsalis)))
d <- expression(paste(italic(Cx.)," ",italic(univittatus)))
legend(x = 9, y = 0.31, bty = "n", legend = c(a, b, c, d), pch = 19, col = c("grey25", "red2", "blue", "darkorange"), cex = 1.2, lwd = 1.5)
legend(x = 6, y = 0.31, bty = "n", legend = c("", "", "", ""), fill = c(rgb(0.25, 0.25, 0.25, 0.4), col = rgb(0.9, 0.25, 0.25, 0.2), rgb(0.25, 0.25, 0.9, 0.4), rgb(1, 0.8, 0, 0.4)), border = c("white", "white", "white"), cex = 1.2)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.1)

######### B: PDR for Ctar
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,0.305), data = CpipWNV.PDR.avg, xaxt = "n", pch = 19, main = expression(paste(italic(Cx.)," ",italic(tarsalis))),
     ylab = expression(paste("Rate (day"^-1,")")), xlab = "", cex.lab = 1.15, col = "white", cex.lab = 1.2)
axis(side = 1, labels = FALSE)
polygon(cord.x, CtarWNV.PDR.cord.y, col = rgb(0.2, 0.5, 1, 0.2), border = NA)
polygon(cord.x, CtarSLEV.PDR.cord.y, col = rgb(0.1, 0.2, 0.9, 0.3), border = NA)
polygon(cord.x, CtarWEEV.PDR.cord.y, col = rgb(0.2, 0.7, 1, 0.2), border = NA)
lines(PDR.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "blue", lwd = 1.5)
lines(PDR.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "dodgerblue", lwd = 1.5)
lines(PDR.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "blue4", lwd = 1.5)
arrows(CtarWEEV.PDR.avg$Temp, CtarWEEV.PDR.avg$upper, CtarWEEV.PDR.avg$Temp, CtarWEEV.PDR.avg$lower, length = 0, lwd = 1.5, col = "dodgerblue")
points(CtarWEEV.PDR.avg$mean ~ CtarWEEV.PDR.avg$Temp, pch = 19, col = "dodgerblue", cex = 1.25)
points(1/trait ~ T, data.PDR.CtarWNV, pch = 19, col = "blue", cex = 1.25)
points(1/trait ~ T, data.PDR.CtarSLEV, pch = 19, col = "blue4", cex = 1.25)
legend(x = 8, y = 0.3, bty = "n", legend = c("WNV", "WEEV", "SLEV"), pch = 19, col = c("blue", "dodgerblue", "blue4"), cex = 1.2, lwd = 1.5)
legend(x = 5, y = 0.3, bty = "n", legend = c("", "", ""), fill = c(rgb(0.2, 0.5, 1, 0.3), rgb(0.2, 0.7, 1, 0.3),  rgb(0.1, 0.2, 0.9, 0.4)), border = c("white", "white", "white"), cex = 1.2)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.1)

######### C: Everything else in Aedes
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,0.47), data = CpipWNV.PDR.avg, pch = 19, main = expression(paste("Other vector-virus pairs")),
     ylab = expression(paste("Rate (day"^-1,")")), xlab = "", cex.lab = 1.15, col = "white", cex.lab = 1.2)
polygon(cord.x, AtaeRVFV.PDR.cord.y, col = rgb(0.1, 0.9, 0.25, 0.2), border = NA)
polygon(cord.x, AtriEEEV.PDR.cord.y, col = rgb(0.5, 0.1, 0.8, 0.2), border = NA)
lines(PDR.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "darkorchid", lwd = 1.5)
lines(PDR.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "mediumseagreen", lwd = 1.5)
points(1/trait ~ T, data.PDR.AtriEEEV, pch = 19, col = "darkorchid", cex = 1.25)
points(1/trait ~ T, data.PDR.AtaeSINV, pch = 19, col = "darkgreen", cex = 1.25)
points(1/trait ~ T, data.PDR.AtaeRVFV, pch = 19, col = "mediumseagreen", cex = 1.25)
a <-expression(paste("EEEV | ",italic(Ae.)," ",italic(tri.)))
b <-expression(paste("RVFV | ",italic(Ae.)," ",italic(tae.)))
c <-expression(paste("SINV | ",italic(Ae.)," ",italic(tae.)))
legend(x = 10, y = 0.475, bty = "n", legend = c(a, b, c), pch = 19, col = c("darkorchid", "mediumseagreen", "darkgreen"), pt.cex = c(1.25, 1.25, 1.4), cex = 1.2, lwd = c(1.5, 1.5, 0))
legend(x = 7, y = 0.475, bty = "n", legend = c("", "", ""), fill = c(rgb(0.5, 0.1, 0.8, 0.2), rgb(0.1, 0.9, 0.25, 0.2), "white"), border = c("white", "white", "white"), cex = 1.2)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.1)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)


##########
###### Figures 6: bc
##########

load("jagsout_bc_CpipWNV_inf.Rdata")
load("jagsout_c_CpipWNV_inf.Rdata")
load("jagsout_b_CtarWNV_inf.Rdata")
load("jagsout_bc_CuniWNV_inf.Rdata")
load("jagsout_bc_CtarWEEV_inf.Rdata")
load("jagsout_c_CtarWEEV_inf.Rdata")
load("jagsout_b_CtarWEEV_inf.Rdata")
load("jagsout_c_CtarSLEV_inf.Rdata")
load("jagsout_b_CtarSLEV_inf.Rdata")
load("jagsout_bc_AtriEEEV_inf.Rdata")
load("jagsout_bc_AtaeRVFV_inf.Rdata")
load("jagsout_c_CpipSINV_inf.Rdata")
load("jagsout_c_AtaeSINV_inf.Rdata")

######### Make polygon coordinates
cord.x <- c(Temp.xs, rev(Temp.xs))

CpipWNV.bc.cord.y <- c(bc.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(bc.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
CpipWNV.c.cord.y <- c(c.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(c.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

CtarWNV.b.cord.y <- c(b.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(b.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

CuniWNV.bc.cord.y <- c(bc.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(bc.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

CtarWEEV.bc.cord.y <- c(bc.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(bc.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
CtarWEEV.c.cord.y <- c(c.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(c.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
CtarWEEV.b.cord.y <- c(b.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(b.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
CtarSLEV.c.cord.y <- c(c.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(c.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
CtarSLEV.b.cord.y <- c(b.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(b.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

AtriEEEV.bc.cord.y <- c(bc.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(bc.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
AtaeRVFV.bc.cord.y <- c(bc.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(bc.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

CpipSINV.c.cord.y <- c(c.CpipSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(c.CpipSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
AtaeSINV.c.cord.y <- c(c.AtaeSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(c.AtaeSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))

a <- expression(paste(italic(Cx.)," ",italic(pipiens)))
b <- expression(paste(italic(Cx.)," ",italic(tarsalis)))
c <- expression(paste(italic(Cx.)," ",italic(univittatus)))

par(mfrow = c(3,3), mar = c(1, 4.5, 1, 1), oma = c(4, 3, 3, 0), las = 1)

######### A: c for WNV
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CpipWNV.bc.avg, xaxt = "n", pch = 19, main = "",
     ylab = "Proportion", xlab = "", cex.lab = 1.15, col = "white", cex.lab = 1.2)
axis(side = 1, labels = FALSE)
polygon(cord.x, CpipWNV.c.cord.y, col = rgb(0.25, 0.25, 0.25, 0.2), border = NA)
lines(c.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "grey25", lwd = 1.5)
arrows(CpipWNV.c.avg$Temp,CpipWNV.c.avg$upper, CpipWNV.c.avg$Temp, CpipWNV.c.avg$lower, length = 0, lwd = 1.5, col = "grey25")
points(CpipWNV.c.avg$mean ~ CpipWNV.c.avg$Temp, pch = 19, col = "grey25", cex = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.1)
text(x = 30, y = 0.2, labels = a, cex = 1.2)
mtext(text = "WNV", side = 2, line = 5, cex = 1.2, las = 0)
mtext(text = expression(paste("Infection efficiency (",italic(c),")")), side = 3, line = 2, cex = 1.2)
mtext(text = expression(paste("# infected / # exposed")), side = 3, line = 0.5, cex = 0.9)
#mtext(text = expression(paste("Infection rate (",italic(c),")")), side = 3, line = 1, cex = 1.2)

######### B: b for WNV
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CpipWNV.bc.avg, xaxt = "n", pch = 19, main = "",
     ylab = "", xlab = "", cex.lab = 1.15, col = "white", cex.lab = 1.2)
axis(side = 1, labels = FALSE)
polygon(cord.x, CtarWNV.b.cord.y, col = rgb(0.25, 0.25, 0.9, 0.2), border = NA)
lines(b.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "blue", lwd = 1.5)
arrows(CtarWNV.b.avg$Temp,CtarWNV.b.avg$upper, CtarWNV.b.avg$Temp, CtarWNV.b.avg$lower, length = 0, lwd = 1.5, col = "blue")
points(CtarWNV.b.avg$mean ~ CtarWNV.b.avg$Temp, pch = 19, col = "blue", cex = 1.25)
text(x = 24, y = 0.1, labels = b, cex = 1.2)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.1)
mtext(text = expression(paste("Transmission efficiency (",italic(b),")")), side = 3, line = 2, cex = 1.2)
mtext(text = expression(paste("# transmitting / # infected")), side = 3, line = 0.5, cex = 0.9)

######### C: bc for WNV
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CpipWNV.bc.avg, xaxt = "n", pch = 19, main = "",
     ylab = "", xlab = "", cex.lab = 1.15, col = "white", cex.lab = 1.2)
axis(side = 1, labels = FALSE)
polygon(cord.x, CpipWNV.bc.cord.y, col = rgb(0.25, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, CuniWNV.bc.cord.y, col = rgb(1, 0.8, 0, 0.3), border = NA)
lines(bc.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "grey25", lwd = 1.5)
lines(bc.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "darkorange", lwd = 1.5)
arrows(CpipWNV.bc.avg$Temp,CpipWNV.bc.avg$upper, CpipWNV.bc.avg$Temp, CpipWNV.bc.avg$lower, length = 0, lwd = 1.5, col = "grey25")
arrows(CuniWNV.bc.avg$Temp,CuniWNV.bc.avg$upper, CuniWNV.bc.avg$Temp, CuniWNV.bc.avg$lower, length = 0, lwd = 1.5, col = "darkorange")
points(CpipWNV.bc.avg$mean ~ CpipWNV.bc.avg$Temp, pch = 19, col = "grey25", cex = 1.25)
points(CuniWNV.bc.avg$mean ~ CuniWNV.bc.avg$Temp, pch = 19, col = "darkorange", cex = 1.25)
text(x = 22.5, y = 0.5, labels = c, cex = 1.2)
text(x = 27, y = 0.07, labels = a, cex = 1.2)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.1)
mtext(text = expression(paste("Vector competence (",italic(bc),")")), side = 3, line = 2, cex = 1.2)
mtext(text = expression(paste("# transmitting / # exposed")), side = 3, line = 0.5, cex = 0.9)

# rgb(1, 0.8, 0, 0.3) yellow
# rgb(0.2, 0.5, 0.5, 0.3) darkcyan
# rgb(0.25, 0.25, 0.9, 0.2), rgb(0.6, 0.4, 0.9, 0.3) purple slate

######### D: c for WEEV
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CpipWNV.bc.avg, xaxt = "n", pch = 19, main = "",
     ylab = "Proportion", xlab = "", cex.lab = 1.15, col = "white", cex.lab = 1.2)
axis(side = 1, labels = FALSE)
polygon(cord.x, CtarWEEV.c.cord.y, col = rgb(0.25, 0.5, 0.9, 0.2), border = NA)
lines(c.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "dodgerblue", lwd = 1.5)
arrows(CtarWEEV.c.avg$Temp,CtarWEEV.c.avg$upper, CtarWEEV.c.avg$Temp, CtarWEEV.c.avg$lower, length = 0, lwd = 1.5, col = "dodgerblue")
points(CtarWEEV.c.avg$mean ~ CtarWEEV.c.avg$Temp, pch = 19, col = "dodgerblue", cex = 1.25)
legend("topleft", legend = "D", bty = "n", adj = 1, cex = 1.1)
mtext(text = "WEEV", side = 2, line = 5, cex = 1.2, las = 0)
text(x = 20, y = 0.6, labels = b, cex = 1.2)

######### E: b for WEEV
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CpipWNV.bc.avg, xaxt = "n", pch = 19, main = "",
     ylab = "", xlab = "", cex.lab = 1.15, col = "white", cex.lab = 1.2)
axis(side = 1, labels = FALSE)
polygon(cord.x, CtarWEEV.b.cord.y, col = rgb(0.25, 0.5, 0.9, 0.2), border = NA)
lines(b.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "dodgerblue", lwd = 1.5)
arrows(CtarWEEV.b.avg$Temp,CtarWEEV.b.avg$upper, CtarWEEV.b.avg$Temp, CtarWEEV.b.avg$lower, length = 0, lwd = 1.5, col = "dodgerblue")
points(CtarWEEV.b.avg$mean ~ CtarWEEV.b.avg$Temp, pch = 19, col = "dodgerblue", cex = 1.25)
legend("topleft", legend = "E", bty = "n", adj = 1, cex = 1.1)
text(x = 22, y = 0.75, labels = b, cex = 1.2)

######### F: bc for WEEV
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CpipWNV.bc.avg, xaxt = "n", pch = 19, main = "",
     ylab = "", xlab = "", cex.lab = 1.15, col = "white", cex.lab = 1.2)
axis(side = 1, labels = FALSE)
polygon(cord.x, CtarWEEV.bc.cord.y, col = rgb(0.25, 0.5, 0.9, 0.2), border = NA)
lines(bc.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "dodgerblue", lwd = 1.5)
arrows(CtarWEEV.bc.avg$Temp,CtarWEEV.bc.avg$upper, CtarWEEV.bc.avg$Temp, CtarWEEV.bc.avg$lower, length = 0, lwd = 1.5, col = "dodgerblue")
points(CtarWEEV.bc.avg$mean ~ CtarWEEV.bc.avg$Temp, pch = 19, col = "dodgerblue", cex = 1.25)
legend("topleft", legend = "F", bty = "n", adj = 1, cex = 1.1)
text(x = 22, y = 0.5, labels = b, cex = 1.2)

######### G: c for SINV + SLEV
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CpipWNV.bc.avg, pch = 19, main = "",
     ylab = "Proportion", xlab = "", col = "white", cex.lab = 1.2)
polygon(cord.x, CpipSINV.c.cord.y, col = rgb(0.25, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, AtaeSINV.c.cord.y, col = rgb(0.2, 0.5, 0.3, 0.3), border = NA)
polygon(cord.x, CtarSLEV.c.cord.y, col = rgb(0.1, 0.2, 0.9, 0.3), border = NA)
lines(c.CpipSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "grey40", lwd = 1.5)
lines(c.AtaeSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "darkgreen", lwd = 1.5)
lines(c.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "blue4", lwd = 1.5)
arrows(CpipSINV.c.avg$Temp,CpipSINV.c.avg$upper, CpipSINV.c.avg$Temp, CpipSINV.c.avg$lower, length = 0, lwd = 1.5, col = "grey40")
arrows(AtaeSINV.c.avg$Temp,AtaeSINV.c.avg$upper, AtaeSINV.c.avg$Temp, AtaeSINV.c.avg$lower, length = 0, lwd = 1.5, col = "darkgreen")
arrows(CtarSLEV.c.avg$Temp,CtarSLEV.c.avg$upper, CtarSLEV.c.avg$Temp, CtarSLEV.c.avg$lower, length = 0, lwd = 1.5, col = "slateblue4")
points(CpipSINV.c.avg$mean ~ CpipSINV.c.avg$Temp, pch = 19, col = "grey40", cex = 1.25)
points(AtaeSINV.c.avg$mean ~ AtaeSINV.c.avg$Temp, pch = 19, col = "darkgreen", cex = 1.25)
points(CtarSLEV.c.avg$mean ~ CtarSLEV.c.avg$Temp, pch = 19, col = "blue4", cex = 1.25)
d <-expression(paste("SINV | ",italic(Cx.)," ",italic(pip.)))
e <-expression(paste("SINV | ",italic(Ae.)," ",italic(tae.)))
f <-expression(paste("SLEV | ",italic(Cx.)," ",italic(tar.)))
legend(x = 17.3, y = 1.05, bty = "n", legend = c(d, e, f), pch = 19, col = c("grey40", "darkgreen", "blue4"), cex = 1.2, lwd = 1.5)
legend(x = 15.3, y = 1.05, bty = "n", legend = c("", "", ""), fill = c(rgb(0.25, 0.25, 0.25, 0.4), rgb(0.2, 0.5, 0.3, 0.35), rgb(0.1, 0.2, 0.9, 0.3)), border = c("white", "white", "white"), cex = 1.2)
legend("topleft", legend = "G", bty = "n", adj = 1, cex = 1.1)
mtext(text = "Other viruses", side = 2, line = 5, cex = 1.2, las = 0)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

######### H: b for SLEV
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CpipWNV.bc.avg, pch = 19, main = "",
     ylab = "", xlab = "", col = "white", cex.lab = 1.2)
polygon(cord.x, CtarSLEV.b.cord.y, col = rgb(0.1, 0.2, 0.9, 0.3), border = NA)
lines(b.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "blue4", lwd = 1.5)
arrows(CtarSLEV.b.avg$Temp, CtarSLEV.b.avg$upper, CtarSLEV.b.avg$Temp, CtarSLEV.b.avg$lower, length = 0, lwd = 1.5, col = "blue4")
points(CtarSLEV.b.avg$mean ~ CtarSLEV.b.avg$Temp, pch = 19, col = "blue4", cex = 1.25)
legend("topleft", legend = "H", bty = "n", adj = 1, cex = 1.1)
g <- expression(paste("SLEV in ",italic(Cx.)," ",italic(tarsalis)))
text(x = 25, y = 0.9, labels = g, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

######### I: bc for RVFV + EEEV
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CpipWNV.bc.avg, pch = 19, main = "",
     ylab = "", xlab = "", col = "white", cex.lab = 1.2)
polygon(cord.x, AtriEEEV.bc.cord.y, col = rgb(0.5, 0.1, 0.8, 0.2), border = NA)
polygon(cord.x, AtaeRVFV.bc.cord.y, col = rgb(0.1, 0.9, 0.25, 0.2), border = NA)
lines(bc.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "darkorchid", lwd = 1.5)
lines(bc.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lty = 1, col = "mediumseagreen", lwd = 1.5)
arrows(AtriEEEV.bc.avg$Temp, AtriEEEV.bc.avg$upper, AtriEEEV.bc.avg$Temp, AtriEEEV.bc.avg$lower, length = 0, lwd = 1.5, col = "darkorchid")
arrows(AtaeRVFV.bc.avg$Temp, AtaeRVFV.bc.avg$upper, AtaeRVFV.bc.avg$Temp, AtaeRVFV.bc.avg$lower, length = 0, lwd = 1.5, col = "dodgerblue")
points(AtriEEEV.bc.avg$mean ~ AtriEEEV.bc.avg$Temp, pch = 19, col = "darkorchid", cex = 1.25)
points(AtaeRVFV.bc.avg$mean ~ AtaeRVFV.bc.avg$Temp, pch = 19, col = "mediumseagreen", cex = 1.25)
g <-expression(paste("EEEV | ",italic(Ae.)," ",italic(tri.)))
h <-expression(paste("RVFV | ",italic(Ae.)," ",italic(tae.)))
legend(x = 11, y = 1.05, bty = "n", legend = c(g, h), pch = 19, col = c("darkorchid", "mediumseagreen"), cex = 1.2, lwd = 1.5)
legend(x = 7, y = 1.05, bty = "n", legend = c("", ""), fill = c(rgb(0.5, 0.1, 0.8, 0.3), rgb(0.1, 0.9, 0.25, 0.3)), border = c("white", "white"), cex = 1.2)
legend("topleft", legend = "I", bty = "n", adj = 1, cex = 1.1)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)