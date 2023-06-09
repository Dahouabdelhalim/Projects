################################################################################
# From Balisi and Van Valkenburgh (2020) Communications Biology
# CODE FOR FIGURE 2:
## Ecomorphology of extinct and survived canid species over 17 time slices
################################################################################


# read in data 
dogdata = read.csv("datasetS1_m1L-BL-DL-m1BS-mass-dur_MB_BVV.csv", stringsAsFactors=FALSE, row.names=1)

# pull m1BS and log10mass from dogdata 
data = dogdata[, 4:5]
# weed out taxa missing one variable or the other 
data.comp = data[complete.cases(data), ]

# convert dataframe columns into their own vectors 
logmass = data.comp$log10mass
names(logmass) = rownames(data.comp)
m1bs = data.comp$m1BS
names(m1bs) = rownames(data.comp)

# Time slices 

# Boundaries (see Supplementary Table 4)
# OREL / WHIT = 32.2 Ma 
# WHIT / EEAK = 30 
# EEAK / LEAK = 27.9 
# LEAK / ELAK = 23.8 
# ELAK / LLAK = 19.5 
# LLAK / EHMF = 18.8 
# EHMF / LHMF = 17.5 
# LHMF / EBAR = 15.9 
# EBAR / LBAR = 14.8 
# LBAR / ECLA = 12.5 
# ECLA / LCLA = 10.75  # boundary set by myself 
# LCLA / EEHP = 9 
# EEHP / LEHP = 7.5
# LEHP / ELHP = 6.7
# ELHP / LLHP = 5.9
# LLHP / BLAN = 4.7 
# BLAN / IRVI = 1.7 
# IRVI / RANC = 0.45 
# RANC / HOLO = 0.01 

# read in species lists per time slice
OREL.csv = read.csv("01_Orellan.csv")
OREL = paste(OREL.csv$GENUS, sep="_", OREL.csv$SPECIES)
WHIT.csv = read.csv("02_Whitneyan.csv")
WHIT = paste(WHIT.csv$GENUS, sep="_", WHIT.csv$SPECIES)
EEAK.csv = read.csv("03_EEAK.csv")
EEAK = paste(EEAK.csv$GENUS, sep="_", EEAK.csv$SPECIES)
LEAK.csv = read.csv("04_LEAK.csv")
LEAK = paste(LEAK.csv$GENUS, sep="_", LEAK.csv$SPECIES)
ELAK.csv = read.csv("05_ELAK.csv")
ELAK = paste(ELAK.csv$GENUS, sep="_", ELAK.csv$SPECIES)
LLAK.csv = read.csv("06_LLAK.csv")
LLAK = paste(LLAK.csv$GENUS, sep="_", LLAK.csv$SPECIES)
EHMF.csv = read.csv("07_EHMF.csv")
EHMF = paste(EHMF.csv$GENUS, sep="_", EHMF.csv$SPECIES)
LHMF.csv = read.csv("08_LHMF.csv")
LHMF = paste(LHMF.csv$GENUS, sep="_", LHMF.csv$SPECIES)
EBAR.csv = read.csv("09_EBAR.csv")
EBAR = paste(EBAR.csv$GENUS, sep="_", EBAR.csv$SPECIES)
LBAR.csv = read.csv("10_LBAR.csv")
LBAR = paste(LBAR.csv$GENUS, sep="_", LBAR.csv$SPECIES)
ECLA.csv = read.csv("11_ECLA.csv")
ECLA = paste(ECLA.csv$GENUS, sep="_", ECLA.csv$SPECIES)
LCLA.csv = read.csv("12_LCLA.csv")
LCLA = paste(LCLA.csv$GENUS, sep="_", LCLA.csv$SPECIES)
EEHP.csv = read.csv("13_EEHP.csv")
EEHP = paste(EEHP.csv$GENUS, sep="_", EEHP.csv$SPECIES)
LEHP.csv = read.csv("14_LEHP.csv")
LEHP = paste(LEHP.csv$GENUS, sep="_", LEHP.csv$SPECIES)
EHMP.csv = merge(EEHP.csv, LEHP.csv, all=TRUE)  # merging of EEHP and LEHP 
EHMP = paste(EHMP.csv$GENUS, sep="_", EHMP.csv$SPECIES) 
ELHP.csv = read.csv("15_ELHP.csv")
ELHP = paste(ELHP.csv$GENUS, sep="_", ELHP.csv$SPECIES)
LLHP.csv = read.csv("16_LLHP.csv")
LLHP = paste(LLHP.csv$GENUS, sep="_", LLHP.csv$SPECIES)
LHMP.csv = merge(ELHP.csv, LLHP.csv, all=TRUE)  # merging of ELHP and LLHP 
LHMP = paste(LHMP.csv$GENUS, sep="_", LHMP.csv$SPECIES)
BLAN.csv = read.csv("17_BLAN.csv")
BLAN = paste(BLAN.csv$GENUS, sep="_", BLAN.csv$SPECIES)
IRVI.csv = read.csv("18_IRVI.csv")
IRVI = paste(IRVI.csv$GENUS, sep="_", IRVI.csv$SPECIES)
RANC.csv = read.csv("19_RANC.csv")
RANC = paste(RANC.csv$GENUS, sep="_", RANC.csv$SPECIES)
HOLO.csv = read.csv("20_HOLO.csv")
HOLO = paste(HOLO.csv$GENUS, sep="_", HOLO.csv$SPECIES)


# COLOR-CODED BY SUBFAMILY
# also colors changed for red-blue color-blind folks
# 05/28/2019: changed font sizes for publication
dogcolors = c("green", "red", "blue")

par(mfrow=c(5,4), oma=c(4,4,0,0)+0.1, mar=c(1,2,0,0)+0.1)

plot(10^logmass[c(OREL, WHIT)], m1bs[c(OREL, WHIT)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002))
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[OREL][OREL %in% WHIT], m1bs[OREL][OREL %in% WHIT], 
       pch=c(19, 15, 17)[OREL.csv$X.SUB.FAMILY][OREL %in% WHIT],
       col=dogcolors[OREL.csv$X.SUB.FAMILY][OREL %in% WHIT], cex=3)
points(10^logmass[OREL][! OREL %in% WHIT], m1bs[OREL][! OREL %in% WHIT], 
       pch=c(0, 0, 2)[OREL.csv$X.SUB.FAMILY][! OREL %in% WHIT], 
       col=dogcolors[OREL.csv$X.SUB.FAMILY][! OREL %in% WHIT], cex=3)
axis(side=1, labels=FALSE)
axis(side=2, labels=TRUE, cex.axis=2)
box(which="plot", bty="l")
text("32.2 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(WHIT, EEAK)], m1bs[c(WHIT, EEAK)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[WHIT][WHIT %in% EEAK], m1bs[WHIT][WHIT %in% EEAK], 
       pch=c(19, 15, 17)[WHIT.csv$X.SUB.FAMILY][WHIT %in% EEAK], 
       col=dogcolors[WHIT.csv$X.SUB.FAMILY][WHIT %in% EEAK], cex=3)
points(10^logmass[WHIT][! WHIT %in% EEAK], m1bs[WHIT][! WHIT %in% EEAK], 
       pch=c(0, 0, 2)[WHIT.csv$X.SUB.FAMILY][! WHIT %in% EEAK], 
       col=dogcolors[WHIT.csv$X.SUB.FAMILY][! WHIT %in% EEAK], cex=3)
axis(side=1, labels=FALSE)
axis(side=2, labels=FALSE)
box(which="plot", bty="l")
text("30 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(EEAK, LEAK)], m1bs[c(EEAK, LEAK)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[EEAK][EEAK %in% LEAK], m1bs[EEAK][EEAK %in% LEAK], 
       pch=c(15, 19, 17)[EEAK.csv$X.SUB.FAMILY][EEAK %in% LEAK], 
       col=rainbow(3)[EEAK.csv$X.SUB.FAMILY][EEAK %in% LEAK], cex=3)
points(10^logmass[EEAK][! EEAK %in% LEAK], m1bs[EEAK][! EEAK %in% LEAK], 
       pch=c(0, 1, 2)[EEAK.csv$X.SUB.FAMILY][! EEAK %in% LEAK], 
       col=rainbow(3)[EEAK.csv$X.SUB.FAMILY][! EEAK %in% LEAK], cex=3) 
axis(side=1, labels=FALSE)
axis(side=2, labels=FALSE)
box(which="plot", bty="l")
text("27.9 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(LEAK, ELAK)], m1bs[c(LEAK, ELAK)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[LEAK][LEAK %in% ELAK], m1bs[LEAK][LEAK %in% ELAK], 
       pch=c(15, 19, 17)[LEAK.csv$X.SUB.FAMILY][LEAK %in% ELAK], 
       col=rainbow(3)[LEAK.csv$X.SUB.FAMILY][LEAK %in% ELAK], cex=3)
points(10^logmass[LEAK][! LEAK %in% ELAK], m1bs[LEAK][! LEAK %in% ELAK], 
       pch=c(0, 1, 2)[LEAK.csv$X.SUB.FAMILY][! LEAK %in% ELAK], 
       col=rainbow(3)[LEAK.csv$X.SUB.FAMILY][! LEAK %in% ELAK], cex=3) 
axis(side=1, labels=FALSE)
axis(side=2, labels=FALSE)
box(which="plot", bty="l")
text("23.8 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(ELAK, LLAK)], m1bs[c(ELAK, LLAK)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[ELAK][ELAK %in% LLAK], m1bs[ELAK][ELAK %in% LLAK], 
       pch=c(15, 19, 17)[ELAK.csv$X.SUB.FAMILY][ELAK %in% LLAK], 
       col=rainbow(3)[ELAK.csv$X.SUB.FAMILY][ELAK %in% LLAK], cex=3)
points(10^logmass[ELAK][! ELAK %in% LLAK], m1bs[ELAK][! ELAK %in% LLAK], 
       pch=c(0, 1, 2)[ELAK.csv$X.SUB.FAMILY][! ELAK %in% LLAK], 
       col=rainbow(3)[ELAK.csv$X.SUB.FAMILY][! ELAK %in% LLAK], cex=3) 
axis(side=1, labels=FALSE)
axis(side=2, labels=TRUE, cex.axis=2)
box(which="plot", bty="l")
text("19.5 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(LLAK, EHMF)], m1bs[c(LLAK, EHMF)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[LLAK][LLAK %in% EHMF], m1bs[LLAK][LLAK %in% EHMF], 
       pch=c(15, 19, 17)[LLAK.csv$X.SUB.FAMILY][LLAK %in% EHMF], 
       col=rainbow(3)[LLAK.csv$X.SUB.FAMILY][LLAK %in% EHMF], cex=3)
points(10^logmass[LLAK][! LLAK %in% EHMF], m1bs[LLAK][! LLAK %in% EHMF], 
       pch=c(0, 1, 2)[LLAK.csv$X.SUB.FAMILY][! LLAK %in% EHMF], 
       col=rainbow(3)[LLAK.csv$X.SUB.FAMILY][! LLAK %in% EHMF], cex=3) 
axis(side=1, labels=FALSE)
axis(side=2, labels=FALSE)
box(which="plot", bty="l")
text("18.8 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(EHMF, LHMF)], m1bs[c(EHMF, LHMF)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[EHMF][EHMF %in% LHMF], m1bs[EHMF][EHMF %in% LHMF], 
       pch=c(15, 19, 17)[EHMF.csv$X.SUB.FAMILY][EHMF %in% LHMF], 
       col=rainbow(3)[EHMF.csv$X.SUB.FAMILY][EHMF %in% LHMF], cex=3)
points(10^logmass[EHMF][! EHMF %in% LHMF], m1bs[EHMF][! EHMF %in% LHMF], 
       pch=c(0, 1, 2)[EHMF.csv$X.SUB.FAMILY][! EHMF %in% LHMF], 
       col=rainbow(3)[EHMF.csv$X.SUB.FAMILY][! EHMF %in% LHMF], cex=3) 
axis(side=1, labels=FALSE)
axis(side=2, labels=FALSE)
box(which="plot", bty="l")
text("17.5 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(LHMF, EBAR)], m1bs[c(LHMF, EBAR)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[LHMF][LHMF %in% EBAR], m1bs[LHMF][LHMF %in% EBAR], 
       pch=c(15, 19, 17)[LHMF.csv$X.SUB.FAMILY][LHMF %in% EBAR], 
       col=rainbow(3)[LHMF.csv$X.SUB.FAMILY][LHMF %in% EBAR], cex=3)
points(10^logmass[LHMF][! LHMF %in% EBAR], m1bs[LHMF][! LHMF %in% EBAR], 
       pch=c(0, 1, 2)[LHMF.csv$X.SUB.FAMILY][! LHMF %in% EBAR], 
       col=rainbow(3)[LHMF.csv$X.SUB.FAMILY][! LHMF %in% EBAR], cex=3) 
axis(side=1, labels=FALSE)
axis(side=2, labels=FALSE)
box(which="plot", bty="l")
text("15.9 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(EBAR, LBAR)], m1bs[c(EBAR, LBAR)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[EBAR][EBAR %in% LBAR], m1bs[EBAR][EBAR %in% LBAR], 
       pch=c(15, 19, 17)[EBAR.csv$X.SUB.FAMILY][EBAR %in% LBAR], 
       col=rainbow(3)[EBAR.csv$X.SUB.FAMILY][EBAR %in% LBAR], cex=3)
points(10^logmass[EBAR][! EBAR %in% LBAR], m1bs[EBAR][! EBAR %in% LBAR], 
       pch=c(0, 1, 2)[EBAR.csv$X.SUB.FAMILY][! EBAR %in% LBAR], 
       col=rainbow(3)[EBAR.csv$X.SUB.FAMILY][! EBAR %in% LBAR], cex=3) 
axis(side=1, labels=FALSE)
axis(side=2, labels=TRUE, cex.axis=2)
box(which="plot", bty="l")
text("14.8 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(LBAR, ECLA)], m1bs[c(LBAR, ECLA)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[LBAR][LBAR %in% ECLA], m1bs[LBAR][LBAR %in% ECLA], 
       pch=c(15, 19, 17)[LBAR.csv$X.SUB.FAMILY][LBAR %in% ECLA], 
       col=rainbow(3)[LBAR.csv$X.SUB.FAMILY][LBAR %in% ECLA], cex=3)
points(10^logmass[LBAR][! LBAR %in% ECLA], m1bs[LBAR][! LBAR %in% ECLA], 
       pch=c(0, 1, 2)[LBAR.csv$X.SUB.FAMILY][! LBAR %in% ECLA], 
       col=rainbow(3)[LBAR.csv$X.SUB.FAMILY][! LBAR %in% ECLA], cex=3) 
axis(side=1, labels=FALSE)
axis(side=2, labels=FALSE)
box(which="plot", bty="l")
text("12.5 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(ECLA, LCLA)], m1bs[c(ECLA, LCLA)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[ECLA][ECLA %in% LCLA], m1bs[ECLA][ECLA %in% LCLA], 
       pch=c(15, 19, 17)[ECLA.csv$X.SUB.FAMILY][ECLA %in% LCLA], 
       col=rainbow(3)[ECLA.csv$X.SUB.FAMILY][ECLA %in% LCLA], cex=3)
points(10^logmass[ECLA][! ECLA %in% LCLA], m1bs[ECLA][! ECLA %in% LCLA], 
       pch=c(0, 1, 2)[ECLA.csv$X.SUB.FAMILY][! ECLA %in% LCLA], 
       col=rainbow(3)[ECLA.csv$X.SUB.FAMILY][! ECLA %in% LCLA], cex=3) 
axis(side=1, labels=FALSE)
axis(side=2, labels=FALSE)
box(which="plot", bty="l")
text("10.75 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(LCLA, EHMP)], m1bs[c(LCLA, EHMP)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[LCLA][LCLA %in% EHMP], m1bs[LCLA][LCLA %in% EHMP], 
       pch=c(15, 19, 17)[LCLA.csv$X.SUB.FAMILY][LCLA %in% EHMP], 
       col=rainbow(3)[LCLA.csv$X.SUB.FAMILY][LCLA %in% EHMP], cex=3)
points(10^logmass[LCLA][! LCLA %in% EHMP], m1bs[LCLA][! LCLA %in% EHMP], 
       pch=c(0, 1, 2)[LCLA.csv$X.SUB.FAMILY][! LCLA %in% EHMP], 
       col=rainbow(3)[LCLA.csv$X.SUB.FAMILY][! LCLA %in% EHMP], cex=3) 
axis(side=1, labels=FALSE)
axis(side=2, labels=FALSE)
box(which="plot", bty="l")
text("9 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(EHMP, LHMP)], m1bs[c(EHMP, LHMP)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[EHMP][EHMP %in% LHMP], m1bs[EHMP][EHMP %in% LHMP], 
       pch=c(15, 19, 17)[EHMP.csv$X.SUB.FAMILY][EHMP %in% LHMP], 
       col=rainbow(3)[EHMP.csv$X.SUB.FAMILY][EHMP %in% LHMP], cex=3)
points(10^logmass[EHMP][! EHMP %in% LHMP], m1bs[EHMP][! EHMP %in% LHMP], 
       pch=c(0, 1, 2)[EHMP.csv$X.SUB.FAMILY][! EHMP %in% LHMP], 
       col=rainbow(3)[EHMP.csv$X.SUB.FAMILY][! EHMP %in% LHMP], cex=3) 
axis(side=1, labels=FALSE)
axis(side=2, labels=TRUE, cex.axis=2)
box(which="plot", bty="l")
text("6.7 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(LHMP, BLAN)], m1bs[c(LHMP, BLAN)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[LHMP][LHMP %in% BLAN], m1bs[LHMP][LHMP %in% BLAN], 
       pch=c(15, 19, 17)[LHMP.csv$X.SUB.FAMILY][LHMP %in% BLAN], 
       col=rainbow(3)[LHMP.csv$X.SUB.FAMILY][LHMP %in% BLAN], cex=3)
points(10^logmass[LHMP][! LHMP %in% BLAN], m1bs[LHMP][! LHMP %in% BLAN], 
       pch=c(0, 1, 2)[LHMP.csv$X.SUB.FAMILY][! LHMP %in% BLAN], 
       col=rainbow(3)[LHMP.csv$X.SUB.FAMILY][! LHMP %in% BLAN], cex=3) 
axis(side=1, labels=TRUE, cex.axis=2)
axis(side=2, labels=FALSE)
box(which="plot", bty="l")
text("4.7 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(BLAN, IRVI)], m1bs[c(BLAN, IRVI)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[BLAN][BLAN %in% IRVI], m1bs[BLAN][BLAN %in% IRVI], 
       pch=19, col=rainbow(3)[2], cex=3)
# switched color code to be all green for all Caninae
# and shapes to be all filled circles for living Caninae
points(10^logmass[BLAN][! BLAN %in% IRVI], m1bs[BLAN][! BLAN %in% IRVI], 
       pch=c(0, 1, 2)[BLAN.csv$X.SUB.FAMILY][! BLAN %in% IRVI], 
       col=rainbow(3)[BLAN.csv$X.SUB.FAMILY][! BLAN %in% IRVI], cex=3)
axis(side=1, labels=TRUE, cex.axis=2)
axis(side=2, labels=FALSE)
box(which="plot", bty="l")
text("1.7 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(IRVI, RANC)], m1bs[c(IRVI, RANC)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[IRVI][IRVI %in% RANC], m1bs[IRVI][IRVI %in% RANC], pch=19, 
       col=rainbow(3)[2], cex=3)
points(10^logmass[IRVI][! IRVI %in% RANC], m1bs[IRVI][! IRVI %in% RANC], pch=1, 
       col=rainbow(3)[2], cex=3) 
axis(side=1, labels=TRUE, cex.axis=2)
axis(side=2, labels=FALSE)
box(which="plot", bty="l")
text("0.45 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(10^logmass[c(RANC, HOLO)], m1bs[c(RANC, HOLO)], type="n", log="x",
     xlim=range(10^logmass, na.rm=T), axes=FALSE,
     ylim=c(range(m1bs, na.rm=T)[1], range(m1bs, na.rm=T)[2]+0.002)) 
rect(14.488, 0.111, 50, 0.15, col="whitesmoke", border=NA)
points(10^logmass[RANC][RANC %in% HOLO], m1bs[RANC][RANC %in% HOLO], pch=19, 
       col=rainbow(3)[2], cex=3)
points(10^logmass[RANC][! RANC %in% HOLO], m1bs[RANC][! RANC %in% HOLO], pch=1, 
       col=rainbow(3)[2], cex=3)
axis(side=1, labels=TRUE, cex.axis=2)
axis(side=2, labels=TRUE, cex.axis=2)
box(which="plot", bty="l")
text("0.01 Ma", x=10^0.4, y=0.14, cex=2.5)

plot(logmass[c(LLHP, BLAN)], m1bs[c(LLHP, BLAN)], type="n", bty="n", axes=FALSE)

plot(logmass[c(LLHP, BLAN)], m1bs[c(LLHP, BLAN)], type="n", bty="n", axes=FALSE,
     xlim=range(logmass, na.rm=T), ylim=range(m1bs, na.rm=T))
legend(0.8, 0.14, legend=c("extinct", "survived"), pch=c(1, 19), col="black", 
       bty="n", cex=2.5, x.intersp=0.5, y.intersp=0.5, pt.cex=3)

plot(logmass[c(LLHP, BLAN)], m1bs[c(LLHP, BLAN)], type="n", axes=FALSE,
     xlim=range(logmass, na.rm=T), ylim=range(m1bs, na.rm=T))
legend(0.2, 0.14, legend=c("Hesperocyoninae", "Borophaginae", "Caninae"), 
       pch=c(17, 15, 19), cex=2.5, col=c("blue", "red", "green"), bty="n",
       x.intersp=0.5, y.intersp=0.5, pt.cex=3)

title(xlab="body mass (kg)", ylab="carnivory (m1BS)", outer=TRUE, line=2, 
      cex.lab=3)



# the actual logistic regression 
# check that this is bootstrapped properly - May 5, 2018


# survived = EEAK in LEAK; extinct = EEAK not in LEAK 
# 27.9 Ma

EEAK0 = data.comp[EEAK, ][! EEAK %in% LEAK, ]
EEAK0$surv = 0
EEAK1 = data.comp[EEAK, ][EEAK %in% LEAK, ]
EEAK1$surv = 1
EEAKbind = rbind(EEAK0, EEAK1)
EEAKbind = EEAKbind[complete.cases(EEAKbind), ]

cor(EEAKbind)
with(EEAKbind, tapply(m1BS, surv, median))  # median of EEAKbind$m1BS by EEAKbind$surv
  # 0 = 0.1040; 1 = 0.1015
with(EEAKbind, tapply(log10mass, surv, median))  # median of EEAKbind$mass by EEAKbind$surv
  # 0 = 0.9410; 1 = 0.6345

glm.out.inter = glm(surv ~ m1BS * log10mass, family=binomial(logit), data=EEAKbind)
summary(glm.out.inter)

# added March 6, 2018
summary(glm.out.inter, se="boot", R=10000, bsmethod="xy")
# this doesn't work--just outputs the same results as summary(glm.out.inter)

# added May 5, 2018
glm.out.add = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=EEAKbind)
summary(glm.out.add, se="boot", R=10000, bsmethod="xy")


################################################################################
# added March 6, 2018
library(boot)

EEAKbind$surv = as.factor(EEAKbind$surv)

lmfit = function(data, indices) {
  fit = glm(surv ~ m1BS * log10mass, family=binomial(logit), data=data[indices, ])
  return(coef(fit))
}

results = boot(data=EEAKbind, statistic=lmfit, R=10000, strata=EEAKbind$surv)

R = 10000  # number of bootstraps

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.4166583      0.5855414      0.5946405      0.4055594 
################################################################################

1 - pchisq(37.096-33.611, df=26-23)  # p = 0.3227128 for interactive model
1 - pchisq(37.096-35.664, df=26-24)  # p = 0.4887032 for additive model
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # no terms producing a sig reduction in deviance from null 

plot(glm.out$fitted)


# survived = LEAK in ELAK; extinct = LEAK not in ELAK 
# 23.8 Ma

LEAK0 = data.comp[LEAK, ][! LEAK %in% ELAK, ]
LEAK0$surv = 0
LEAK1 = data.comp[LEAK, ][LEAK %in% ELAK, ]
LEAK1$surv = 1
LEAKbind = rbind(LEAK0, LEAK1)
LEAKbind = LEAKbind[complete.cases(LEAKbind), ]

cor(LEAKbind)
with(LEAKbind, tapply(m1BS, surv, median))
with(LEAKbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=LEAKbind)
summary(glm.out)

################################################################################
# added March 6, 2018
LEAKbind$surv = as.factor(LEAKbind$surv)

results = boot(data=LEAKbind, statistic=lmfit, R=10000, strata=LEAKbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.5683432      0.4292571      0.3988601      0.5946405 
################################################################################

1 - pchisq(24.435-22.968, df=19-16)  # p = 0.6899074 for interactive model
1 - pchisq(24.435-24.276, df=19-17)  # p = 0.923578 for additive model
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # no terms producing a sig reduction in deviance from null 

plot(glm.out$fitted)


# survived = ELAK in LLAK; extinct = ELAK not in LLAK 
# 19.5 Ma

ELAK0 = data.comp[ELAK, ][! ELAK %in% LLAK, ]
ELAK0$surv = 0
ELAK1 = data.comp[ELAK, ][ELAK %in% LLAK, ]
ELAK1$surv = 1
ELAKbind = rbind(ELAK0, ELAK1)
ELAKbind = ELAKbind[complete.cases(ELAKbind), ]

cor(ELAKbind)
with(ELAKbind, tapply(m1BS, surv, median))
with(ELAKbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=ELAKbind)
summary(glm.out)

################################################################################
# added March 6, 2018
ELAKbind$surv = as.factor(ELAKbind$surv)

results = boot(data=ELAKbind, statistic=lmfit, R=10000, strata=ELAKbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.3782622      0.6195380      0.6374363      0.3668633 
################################################################################

1 - pchisq(25.864-21.547, df=18-15)  # p = 0.2292059 for interactive model
1 - pchisq(25.864-21.547, df=18-16)  # p = 0.1154982 for additive model
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # m1bs p = 0.01824 *******************************
# m1bs produces a SIGNIFICANT reduction in deviance from null 

plot(glm.out$fitted)


# survived = LLAK in EHMF; extinct = LLAK not in EHMF 
# 18.8 Ma

LLAK0 = data.comp[LLAK, ][! LLAK %in% EHMF, ]
LLAK0$surv = 0
LLAK1 = data.comp[LLAK, ][LLAK %in% EHMF, ]
LLAK1$surv = 1
LLAKbind = rbind(LLAK0, LLAK1)
LLAKbind = LLAKbind[complete.cases(LLAKbind), ]

cor(LLAKbind)
with(LLAKbind, tapply(m1BS, surv, median))
with(LLAKbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=LLAKbind)
summary(glm.out)

################################################################################
# added March 6, 2018
LLAKbind$surv = as.factor(LLAKbind$surv)

results = boot(data=LLAKbind, statistic=lmfit, R=10000, strata=LLAKbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.5342466      0.4583542      0.4919508      0.5221478 
################################################################################

1 - pchisq(16.636-16.292, df=11-8)  # p = 0.9515515 for interactive model
1 - pchisq(16.636-16.293, df=11-9)  # p = 0.8424003 for additive model
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # m1bs:log10mass p = 0.06896 #####################
# interaction between m1bs and logmass produces NEARLY sig reduction in deviance from null 

plot(glm.out$fitted)


# survived = EHMF in LHMF; extinct = EHMF not in LHMF 
# 17.5 Ma

EHMF0 = data.comp[EHMF, ][! EHMF %in% LHMF, ]
EHMF0$surv = 0
EHMF1 = data.comp[EHMF, ][EHMF %in% LHMF, ]
EHMF1$surv = 1
EHMFbind = rbind(EHMF0, EHMF1)
EHMFbind = EHMFbind[complete.cases(EHMFbind), ]

cor(EHMFbind)
with(EHMFbind, tapply(m1BS, surv, median))
with(EHMFbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=EHMFbind)
summary(glm.out)

################################################################################
# added March 6, 2018
EHMFbind$surv = as.factor(EHMFbind$surv)

results = boot(data=EHMFbind, statistic=lmfit, R=10000, strata=EHMFbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.4137586      0.5756424      0.5708429      0.4158584 
################################################################################

1 - pchisq(20.597-19.497, df=16-13)  # p = 0.7770741 for interactive model
1 - pchisq(20.597-19.640, df=16-14)  # p = 0.6197123 for additive model
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # no terms producing a sig reduction in deviance from null 

plot(glm.out$fitted)


# survived = LHMF in EBAR; extinct = LHMF not in EBAR 
# 15.9 Ma

LHMF0 = data.comp[LHMF, ][! LHMF %in% EBAR, ]
LHMF0$surv = 0
LHMF1 = data.comp[LHMF, ][LHMF %in% EBAR, ]
LHMF1$surv = 1
LHMFbind = rbind(LHMF0, LHMF1)
LHMFbind = LHMFbind[complete.cases(LHMFbind), ]

cor(LHMFbind)
with(LHMFbind, tapply(m1BS, surv, median))
with(LHMFbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=LHMFbind)
summary(glm.out)

################################################################################
# added March 6, 2018
LHMFbind$surv = as.factor(LHMFbind$surv)

results = boot(data=LHMFbind, statistic=lmfit, R=10000, strata=LHMFbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.7167283      0.2823718             NA             NA 
# What do these NAs mean??
# They persist even without stratification.
# Probably sample size too small...
################################################################################

1 - pchisq(8.9974-4.9819, df=7-4)  # p = 0.2597953 for interactive model
1 - pchisq(8.9974-8.7474, df=7-5)  # p = 0.8824969 for additive model
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # m1bs:log10mass p = 0.05232  ####################
# interaction between m1bs and logmass produces a NEARLY sig reduction in deviance from null 

plot(glm.out$fitted)


# survived = EBAR in LBAR; extinct = EBAR not in LBAR 
# 14.8 Ma

EBAR0 = data.comp[EBAR, ][! EBAR %in% LBAR, ]
EBAR0$surv = 0
EBAR1 = data.comp[EBAR, ][EBAR %in% LBAR, ]
EBAR1$surv = 1
EBARbind = rbind(EBAR0, EBAR1)
EBARbind = EBARbind[complete.cases(EBARbind), ]

cor(EBARbind)
with(EBARbind, tapply(m1BS, surv, median))
with(EBARbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=EBARbind)
summary(glm.out)

################################################################################
# added March 6, 2018
EBARbind$surv = as.factor(EBARbind$surv)

results = boot(data=EBARbind, statistic=lmfit, R=10000, strata=EBARbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.5497450      0.4576542      0.4530547      0.5503450 
################################################################################

1 - pchisq(22.915-21.800, df=17-14)  # p = 0.7734543 for interactive model
1 - pchisq(22.915-21.812, df=17-15)  # p = 0.576085 for additive model
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # no terms producing a sig reduction in deviance from null 

plot(glm.out$fitted)


# survived = LBAR in ECLA; extinct = LBAR not in ECLA
# 12.5 Ma

LBAR0 = data.comp[LBAR, ][! LBAR %in% ECLA, ]
LBAR0$surv = 0
LBAR1 = data.comp[LBAR, ][LBAR %in% ECLA, ]
LBAR1$surv = 1
LBARbind = rbind(LBAR0, LBAR1)
LBARbind = LBARbind[complete.cases(LBARbind), ]

cor(LBARbind)
with(LBARbind, tapply(m1BS, surv, median))
with(LBARbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=LBARbind)
summary(glm.out)

################################################################################
# added March 6, 2018
LBARbind$surv = as.factor(LBARbind$surv)

results = boot(data=LBARbind, statistic=lmfit, R=10000, strata=LBARbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.4757524      0.5144486      0.5518448      0.4639536 
################################################################################

1 - pchisq(25.008-22.132, df=18-15)  # p = 0.4111413 for interactive model
1 - pchisq(25.008-22.680, df=18-16)  # p = 0.3122347 for additive model
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # no terms producing a sig reduction in deviance from null 

plot(glm.out$fitted)


# survived = ECLA in LCLA; extinct = ECLA not in LCLA
# 10.75 Ma

ECLA0 = data.comp[ECLA, ][! ECLA %in% LCLA, ]
ECLA0$surv = 0
ECLA1 = data.comp[ECLA, ][ECLA %in% LCLA, ]
ECLA1$surv = 1
ECLAbind = rbind(ECLA0, ECLA1)
ECLAbind = ECLAbind[complete.cases(ECLAbind), ]

cor(ECLAbind)
with(ECLAbind, tapply(m1BS, surv, median))
with(ECLAbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=ECLAbind)
summary(glm.out)

################################################################################
# added March 6, 2018
ECLAbind$surv = as.factor(ECLAbind$surv)

results = boot(data=ECLAbind, statistic=lmfit, R=10000, strata=ECLAbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.5323468      0.4722528      0.4765523      0.5521448 
################################################################################

1 - pchisq(20.190-20.063, df=14-11)  # p = 0.9884112 for interactive model
1 - pchisq(20.190-20.063, df=14-12)  # p = 0.9384741 for additive model
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # no terms producing a sig reduction in deviance from null 

plot(glm.out$fitted)


# survived = LCLA in EHMP; extinct = LCLA not in EHMP
# 9 Ma

LCLA0 = data.comp[LCLA, ][! LCLA %in% EHMP, ]
LCLA0$surv = 0
LCLA1 = data.comp[LCLA, ][LCLA %in% EHMP, ]
LCLA1$surv = 1
LCLAbind = rbind(LCLA0, LCLA1)
LCLAbind = LCLAbind[complete.cases(LCLAbind), ]

cor(LCLAbind)
with(LCLAbind, tapply(m1BS, surv, median))
with(LCLAbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=LCLAbind)
summary(glm.out)

################################################################################
# added March 6, 2018
LCLAbind$surv = as.factor(LCLAbind$surv)

results = boot(data=LCLAbind, statistic=lmfit, R=10000, strata=LCLAbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.5041496      0.4758524      0.2689731      0.6855314 
################################################################################

1 - pchisq(12.8910-5.9655, df=10-7)  # p = 0.07431063  ################ marginal
1 - pchisq(12.891-8.041, df=10-8)  # p = 0.08847812  ################## marginal
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # no terms producing a sig reduction in deviance from null 

plot(glm.out$fitted)


# survived = EHMP in LHMP; extinct = EHMP not in LHMP
# 6.7 Ma

EHMP0 = data.comp[EHMP, ][! EHMP %in% LHMP, ]
EHMP0$surv = 0
EHMP1 = data.comp[EHMP, ][EHMP %in% LHMP, ]
EHMP1$surv = 1
EHMPbind = rbind(EHMP0, EHMP1)
EHMPbind = EHMPbind[complete.cases(EHMPbind), ]

cor(EHMPbind)
with(EHMPbind, tapply(m1BS, surv, median))
with(EHMPbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=EHMPbind)
summary(glm.out)

################################################################################
# added March 6, 2018
EHMPbind$surv = as.factor(EHMPbind$surv)

results = boot(data=EHMPbind, statistic=lmfit, R=10000, strata=EHMPbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.4242576      0.5740426      0.6207379      0.4226577 
################################################################################

1 - pchisq(17.945-12.640, df=12-9)  # p = 0.1507783 for interactive model
1 - pchisq(17.945-14.199, df=12-10)  # p = 0.153662 for additive model
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # m1bs p = 0.06391  ##############################
# m1BS produces a NEARLY sig reduction in deviance from null 

plot(glm.out$fitted)


# survived = LHMP in BLAN; extinct = LHMP not in BLAN
# 4.7 Ma

LHMP0 = data.comp[LHMP, ][! LHMP %in% BLAN, ]
LHMP0$surv = 0
LHMP1 = data.comp[LHMP, ][LHMP %in% BLAN, ]
LHMP1$surv = 1
LHMPbind = rbind(LHMP0, LHMP1)
LHMPbind = LHMPbind[complete.cases(LHMPbind), ]

cor(LHMPbind)
with(LHMPbind, tapply(m1BS, surv, median))
with(LHMPbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=LHMPbind)
summary(glm.out)

################################################################################
# added March 6, 2018
LHMPbind$surv = as.factor(LHMPbind$surv)

results = boot(data=LHMPbind, statistic=lmfit, R=10000, strata=LHMPbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.6178382      0.3804620      0.4034597             NA 
################################################################################

1 - pchisq(9.5347-6.6741, df=8-5)  # p = 0.413621 for interactive model
1 - pchisq(9.5347-6.9068, df=8-6)  # p = 0.2687564 for additive model
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # no terms producing a sig reduction in deviance from null 

plot(glm.out$fitted)


# survived = BLAN in IRVI; extinct = BLAN not in IRVI
# 1.7 Ma

BLAN0 = data.comp[BLAN, ][! BLAN %in% IRVI, ]
BLAN0$surv = 0
BLAN1 = data.comp[BLAN, ][BLAN %in% IRVI, ]
BLAN1$surv = 1
BLANbind = rbind(BLAN0, BLAN1)
BLANbind = BLANbind[complete.cases(BLANbind), ]

cor(BLANbind)
with(BLANbind, tapply(m1BS, surv, median))
with(BLANbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=BLANbind)
summary(glm.out)

################################################################################
# added March 6, 2018
BLANbind$surv = as.factor(BLANbind$surv)

results = boot(data=BLANbind, statistic=lmfit, R=10000, strata=BLANbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.4765523      0.5265473      0.5636436             NA 
################################################################################

1 - pchisq(11.457-10.273, df=8-5)  # p = 0.7568443 for interactive model
1 - pchisq(11.457-10.273, df=8-6)  # p = 0.5532197 for additive model
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # no terms producing a sig reduction in deviance from null 

plot(glm.out$fitted)


# survived = IRVI in RANC; extinct = IRVI not in RANC
# 0.45 Ma

IRVI0 = data.comp[IRVI, ][! IRVI %in% RANC, ]
IRVI0$surv = 0
IRVI1 = data.comp[IRVI, ][IRVI %in% RANC, ]
IRVI1$surv = 1
IRVIbind = rbind(IRVI0, IRVI1)
IRVIbind = IRVIbind[complete.cases(IRVIbind), ]

cor(IRVIbind)
with(IRVIbind, tapply(m1BS, surv, median))
with(IRVIbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=IRVIbind)
summary(glm.out)

################################################################################
# added March 6, 2018
IRVIbind$surv = as.factor(IRVIbind$surv)

results = boot(data=IRVIbind, statistic=lmfit, R=10000, strata=IRVIbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.6090391      0.4691531      0.4269573             NA 
################################################################################

1 - pchisq(15.158-11.200, df=10-7)  # p = 0.2660352 for interactive model
1 - pchisq(15.158-12.399, df=10-8)  # p = 0.2517044 for additive model
# model appears to have performed poorly: 
# no significant reduction in deviance - no significant difference from null model 

anova(glm.out, test="Chisq")  # no terms producing a sig reduction in deviance from null 

plot(glm.out$fitted)


# survived = RANC in HOLO; extinct = RANC not in HOLO
# RANC / HOLO

RANC0 = data.comp[RANC, ][! RANC %in% HOLO, ]
RANC0$surv = 0
RANC1 = data.comp[RANC, ][RANC %in% HOLO, ]
RANC1$surv = 1
RANCbind = rbind(RANC0, RANC1)
RANCbind = RANCbind[complete.cases(RANCbind), ]

cor(RANCbind)
with(RANCbind, tapply(m1BS, surv, median))
with(RANCbind, tapply(log10mass, surv, median))

glm.out = glm(surv ~ m1BS + log10mass, family=binomial(logit), data=RANCbind)
# note: swapping m1BS and logmass produces the same summary(glm.out) 
##  BUT different anova(glm.out)
summary(glm.out)

################################################################################
# added March 6, 2018
RANCbind$surv = as.factor(RANCbind$surv)

results = boot(data=RANCbind, statistic=lmfit, R=10000)#, strata=RANCbind$surv)

extrema = apply(X=results$t, MARGIN=1, FUN=">", results$t0)

pvals = rowSums(extrema)/(R+1)

# with strata
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.4876512      0.4442556      0.4612539             NA

# without strata
# (Intercept)           m1BS      log10mass m1BS:log10mass 
#   0.4940506      0.4608539      0.4412559      0.5823418 
################################################################################

1 - pchisq(1.2891e+01 - 5.3791e-10, df=10-7)  # p = 0.004878353 ****************
1 - pchisq(1.2891e+01 - 6.1904e-10, df=10-8)  # p = 0.001587651 ****************
# model appears to have performed well: 
# significant reduction in deviance - significant difference from the null model 

anova(glm.out, test="Chisq")  # p = 0.0003302 **********************************
# m1bs ON ITS OWN produced a significant reduction in deviance 
# of 12.891 on 1 degree of freedom 
# but if I swap m1bs and log10mass... both factors contribute

plot(glm.out$fitted)