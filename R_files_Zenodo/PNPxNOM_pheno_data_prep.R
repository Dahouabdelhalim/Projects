# getting phenotype data ready for PNPxNOM QTL mapping 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(here)

# read in data
colour<-read.table(file = "PNPxNOM_colourscores_done.csv", header = T, sep = ",") # 132 F2 males (+ 2 added parental type scores)
teeth<-read.table(file = "PNPxNOM_teeth_done.csv", header = T, sep = ",") # 161 F2s (132 males, 29 females), 20 parentals
guts<-read.table(file = "PNPxNOM_guts_done.csv", header = T, sep = ",") # 161 F2s (132 males, 29 females), 12 parentals
morphology<-read.table(file = "PNPxNOM_morphology_done_adj.csv", header = T, sep = ",") # 158 F2s (129 males (194026, 194206, 194219 missing), 29 females, 1 NA (194205)), 25 parentals 


# colour
#@@@@@@@
# remove unneccessary columns (and parental type rows)
colourF2s<-colour[-c(133:135),-c(20,22:25)]


# teeth 
#@@@@@@@
# calculate scores per jaw (cumulative score of tooth shape divided by number of teeth (excluding broken/missing teeth))
# L1/U1 are unicuspid (=0), L5/U5 are bicuspid (=1), L2/U2-L4/U4 in 0.25 steps inbetween. 
# L6/U6 are tricuspid (transgressive). not included in score but separately as presence/absence.
teeth[25]<-(teeth$L1*0+teeth$L2*0.25+teeth$L3*0.5+teeth$L4*0.75+teeth$L5*1)/(10-(teeth$Lmiss+teeth$Lbrok+teeth$L6))
colnames(teeth)[25]<-"Lscore"
teeth[26]<-(teeth$U1*0+teeth$U2*0.25+teeth$U3*0.5+teeth$U4*0.75+teeth$U5*1)/(10-(teeth$Umiss+teeth$Ubrok+teeth$U6))
colnames(teeth)[26]<-"Uscore"
teeth[27]<-ifelse(teeth$L6>0 | teeth$U6>0, 'yes', 'no')
colnames(teeth)[27]<-"tricuspid(T)"
# now keep only scores and F2s for mapping
teethF2s<-teeth[c(1:161),c(1,18:21,23:27)]

# for comparisons with a few parental species individuals
# subset to F2 males
teethF2males<-teethF2s[teethF2s$sex==1,]
# subset to F2 females
teethF2females<-teethF2s[teethF2s$sex==0,]
# subset to parental species (only males)
teethPNPlab<-teeth[c(174:177),c(1,18:21,23:27)]
teethPNPwild<-teeth[c(166:169),c(1,18:21,23:27)]
teethNOMlab<-teeth[c(170:173),c(1,18:21,23:27)]
teethNOMwild<-teeth[c(162:165),c(1,18:21,23:27)]

# plot count histograms of all types
par(mfrow=c(7,2))
  hist(x = teethF2s$Lscore, breaks = 20, xlab="shape score lower jaw F2 all", main="", col="grey50", freq = F, xlim = c(0,1))
    lines(density(teethF2s$Lscore, na.rm=T))
  hist(x = teethF2s$Uscore, breaks = 20, xlab="shape score upper jaw F2 all", main="", col="grey50", freq = F, xlim = c(0,1))
    lines(density(teethF2s$Lscore, na.rm=T))
  hist(x = teethF2males$Lscore, breaks = 20, xlab="shape score lower jaw F2 males", main="", col="grey50", freq = F, xlim = c(0,1))
   lines(density(teethF2males$Lscore, na.rm=T))
  hist(x = teethF2males$Uscore, breaks = 20, xlab="shape score upper jaw F2 males", main="", col="grey50", freq = F, xlim = c(0,1))
    lines(density(teethF2males$Lscore, na.rm=T))
  hist(x = teethF2females$Lscore, breaks = 20, xlab="shape score lower jaw F2 females", main="", col="grey50", freq = F, xlim = c(0,1))
   lines(density(teethF2females$Lscore, na.rm=T))
    hist(x = teethF2females$Uscore, breaks = 20, xlab="shape score upper jaw F2 females", main="", col="grey50", freq = F, xlim = c(0,1))
  lines(density(teethF2females$Lscore, na.rm=T))
  hist(x = teethPNPlab$Lscore, breaks = 5, xlab="shape score lower jaw PNP males lab", main="", col="grey50", freq = F, xlim = c(0,1))
  hist(x = teethPNPlab$Uscore, breaks = 5, xlab="shape score upper jaw PNP males lab", main="", col="grey50", freq = F, xlim = c(0,1))
  hist(x = teethNOMlab$Lscore, breaks = 5, xlab="shape score lower jaw NOM males lab", main="", col="grey50", freq = F, xlim = c(0,1))
  hist(x = teethNOMlab$Uscore, breaks = 5, xlab="shape score upper jaw NOM males lab", main="", col="grey50", freq = F, xlim=c(0,1))
  hist(x = teethPNPwild$Lscore, breaks = 5, xlab="shape score lower jawPNP males wild", main="", col="grey50", freq = F, xlim = c(0,1))
  hist(x = teethPNPwild$Uscore, breaks = 5, xlab="shape score upper jaw PNP males wild", main="", col="grey50", freq = F, xlim=c(0,1))
  hist(x = teethNOMwild$Lscore, breaks = 5, xlab="shape score lower jaw NOM males wild", main="", col="grey50", freq = F, xlim = c(0,1))
  hist(x = teethNOMwild$Uscore, breaks = 5, xlab="shape score upper jaw NOM males wild", main="", col="grey50", freq = F, xlim=c(0,1))


# intestine length
#@@@@@@@@@@@@@@@@@

# plot against SL
# read in standard length
stdlength<-morphology[,1:4]
# calculate mean SL (from the two measurements)
stdlength[5]<-(stdlength[2]+stdlength[3])/2
colnames(stdlength)[5]<-"meanSL"
# combine the data
gutsSL<-merge(x = stdlength[c(1,5)], y = guts, by="id", all=F)

# plot it
  cols=c("red", "blue")[as.factor(gutsSL$sex)]
  pchs=c(19,3,4,2,6)[as.factor(gutsSL$type)]
  plot(gutsSL$meanSL, gutsSL$gutlength, xlab="SL in cm", ylab="gut length in cm", main="", cex.main=1, pch=pchs, col=cols, xlim=c(40,80), ylim=c(0,40))
  legend(x = "topleft", legend = c("F2 males", "F2females", "NOMlabmales", "NOMwildmales", "PNPlabmales", "PNPwildmales"), pch=c(19,19,3,4,2,6), cex=0.75, col=c("blue", "red", "blue", "blue", "blue", "blue"))

# test if SL significantly correlated with gutlength
cor.test(x = gutsSL$meanSL, y = gutsSL$gutlength)
cor.test(x = log10(gutsSL$meanSL), y = log10(gutsSL$gutlength))
# check if slopes different
summary(aov(log10(gutsSL$gutlength)~log10(gutsSL$meanSL)*gutsSL$type))
# significant effect of type, but no interaction, so do size-correction for all together
# do size correction
gutsSL$gutlengthcorr<-residuals(lm(log10(gutsSL$gutlength)~log10(gutsSL$meanSL), na.action=na.exclude))

# plot again
  cols=c("red", "blue")[as.factor(gutsSL$sex)]
  pchs=c(19,3,4,2,6)[as.factor(gutsSL$type)]
  plot(log10(gutsSL$meanSL), gutsSL$gutlengthcorr, xlab="log(SL)", ylab="size-corrected gut length", main="", cex.main=1, pch=pchs, col=cols)
  legend(x = "topleft", legend = c("F2 males", "F2females", "NOMlabmales", "NOMwildmales", "PNPlabmales", "PNPwildmales"), pch=c(19,19,3,4,2,6), cex=0.75, col=c("blue", "red", "blue", "blue", "blue", "blue"))

# the same only within F2s
gutsSLF2s<-gutsSL[gutsSL$type=="F2",]
# test if SL significantly correlated with gutlength
cor.test(x = gutsSLF2s$meanSL, y = gutsSLF2s$gutlength)
cor.test(x = log10(gutsSLF2s$meanSL), y = log10(gutsSLF2s$gutlength))
# check if slopes different
summary(aov(log10(gutsSLF2s$gutlength)~log10(gutsSLF2s$meanSL)*gutsSLF2s$sex))
# no significant effect of or interaction with sex, so do size-correction for males and females together
# do size correction
gutsSLF2s$gutlengthcorr<-residuals(lm(log10(gutsSLF2s$gutlength)~log10(gutsSLF2s$meanSL), na.action=na.exclude))


# external morphology 
#@@@@@@@@@@@@@@@@@@@@

# get rid of last two columns
morphologyall<-morphology[,-c(68:69)]

# function to calculate means of the two nearest measurements (see also PNPxNOM_morphology_data_prep_old.R script)
findClosestANDcalculateMean <- function(x) {
  x <- sort(x)
  xnew<-x[seq.int(which.min(diff(as.numeric(x), lag = 1)), length.out = 2)]
  (xnew[,1]+xnew[,2])/2
}

# apply it to every trait

#SL
morphologyallSL<-morphologyall[,c(1,2:4)]
morphologyallSL<-morphologyallSL[!(is.na(morphologyallSL[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallSL$meanSL <- 0
for (i in 1:nrow(morphologyallSL)) {
  morphologyallSL1 <- sort(morphologyallSL[i,2:4])
  morphologyallSL2<-morphologyallSL1[seq.int(which.min(diff(as.numeric(morphologyallSL1), lag = 1)), length.out = 2)]
  morphologyallSL3<-(morphologyallSL2[,1]+morphologyallSL2[,2])/2
  morphologyallSL$meanSL[i]<-morphologyallSL3
}

#BD
morphologyallBD<-morphologyall[,c(1,6:8)]
morphologyallBD<-morphologyallBD[!(is.na(morphologyallBD[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallBD$meanBD <- 0
for (i in 1:nrow(morphologyallBD)) {
  morphologyallBD1 <- sort(morphologyallBD[i,2:4])
  morphologyallBD2<-morphologyallBD1[seq.int(which.min(diff(as.numeric(morphologyallBD1), lag = 1)), length.out = 2)]
  morphologyallBD3<-(morphologyallBD2[,1]+morphologyallBD2[,2])/2
  morphologyallBD$meanBD[i]<-morphologyallBD3
}
#HL
morphologyallHL<-morphologyall[,c(1,10:12)]
morphologyallHL<-morphologyallHL[!(is.na(morphologyallHL[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallHL$meanHL <- 0
for (i in 1:nrow(morphologyallHL)) {
  morphologyallHL1 <- sort(morphologyallHL[i,2:4])
  morphologyallHL2<-morphologyallHL1[seq.int(which.min(diff(as.numeric(morphologyallHL1), lag = 1)), length.out = 2)]
  morphologyallHL3<-(morphologyallHL2[,1]+morphologyallHL2[,2])/2
  morphologyallHL$meanHL[i]<-morphologyallHL3
}
#HW
morphologyallHW<-morphologyall[,c(1,14:16)]
morphologyallHW<-morphologyallHW[!(is.na(morphologyallHW[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallHW$meanHW <- 0
for (i in 1:nrow(morphologyallHW)) {
  morphologyallHW1 <- sort(morphologyallHW[i,2:4])
  morphologyallHW2<-morphologyallHW1[seq.int(which.min(diff(as.numeric(morphologyallHW1), lag = 1)), length.out = 2)]
  morphologyallHW3<-(morphologyallHW2[,1]+morphologyallHW2[,2])/2
  morphologyallHW$meanHW[i]<-morphologyallHW3
}
#LJL
morphologyallLJL<-morphologyall[,c(1,18:20)]
morphologyallLJL<-morphologyallLJL[!(is.na(morphologyallLJL[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallLJL$meanLJL <- 0
for (i in 1:nrow(morphologyallLJL)) {
  morphologyallLJL1 <- sort(morphologyallLJL[i,2:4])
  morphologyallLJL2<-morphologyallLJL1[seq.int(which.min(diff(as.numeric(morphologyallLJL1), lag = 1)), length.out = 2)]
  morphologyallLJL3<-(morphologyallLJL2[,1]+morphologyallLJL2[,2])/2
  morphologyallLJL$meanLJL[i]<-morphologyallLJL3
}
#LJW
morphologyallLJW<-morphologyall[,c(1,22:24)]
morphologyallLJW<-morphologyallLJW[!(is.na(morphologyallLJW[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallLJW$meanLJW <- 0
for (i in 1:nrow(morphologyallLJW)) {
  morphologyallLJW1 <- sort(morphologyallLJW[i,2:4])
  morphologyallLJW2<-morphologyallLJW1[seq.int(which.min(diff(as.numeric(morphologyallLJW1), lag = 1)), length.out = 2)]
  morphologyallLJW3<-(morphologyallLJW2[,1]+morphologyallLJW2[,2])/2
  morphologyallLJW$meanLJW[i]<-morphologyallLJW3
}
#SnL
morphologyallSnL<-morphologyall[,c(1,26:28)]
morphologyallSnL<-morphologyallSnL[!(is.na(morphologyallSnL[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallSnL$meanSnL <- 0
for (i in 1:nrow(morphologyallSnL)) {
  morphologyallSnL1 <- sort(morphologyallSnL[i,2:4])
  morphologyallSnL2<-morphologyallSnL1[seq.int(which.min(diff(as.numeric(morphologyallSnL1), lag = 1)), length.out = 2)]
  morphologyallSnL3<-(morphologyallSnL2[,1]+morphologyallSnL2[,2])/2
  morphologyallSnL$meanSnL[i]<-morphologyallSnL3
}
#SnW
morphologyallSnW<-morphologyall[,c(1,30:32)]
morphologyallSnW<-morphologyallSnW[!(is.na(morphologyallSnW[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallSnW$meanSnW <- 0
for (i in 1:nrow(morphologyallSnW)) {
  morphologyallSnW1 <- sort(morphologyallSnW[i,2:4])
  morphologyallSnW2<-morphologyallSnW1[seq.int(which.min(diff(as.numeric(morphologyallSnW1), lag = 1)), length.out = 2)]
  morphologyallSnW3<-(morphologyallSnW2[,1]+morphologyallSnW2[,2])/2
  morphologyallSnW$meanSnW[i]<-morphologyallSnW3
}
#EyL
morphologyallEyL<-morphologyall[,c(1,34:36)]
morphologyallEyL<-morphologyallEyL[!(is.na(morphologyallEyL[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallEyL$meanEyL <- 0
for (i in 1:nrow(morphologyallEyL)) {
  morphologyallEyL1 <- sort(morphologyallEyL[i,2:4])
  morphologyallEyL2<-morphologyallEyL1[seq.int(which.min(diff(as.numeric(morphologyallEyL1), lag = 1)), length.out = 2)]
  morphologyallEyL3<-(morphologyallEyL2[,1]+morphologyallEyL2[,2])/2
  morphologyallEyL$meanEyL[i]<-morphologyallEyL3
}
#EyD
morphologyallEyD<-morphologyall[,c(1,38:40)]
morphologyallEyD<-morphologyallEyD[!(is.na(morphologyallEyD[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallEyD$meanEyD <- 0
for (i in 1:nrow(morphologyallEyD)) {
  morphologyallEyD1 <- sort(morphologyallEyD[i,2:4])
  morphologyallEyD2<-morphologyallEyD1[seq.int(which.min(diff(as.numeric(morphologyallEyD1), lag = 1)), length.out = 2)]
  morphologyallEyD3<-(morphologyallEyD2[,1]+morphologyallEyD2[,2])/2
  morphologyallEyD$meanEyD[i]<-morphologyallEyD3
}
#ChD
morphologyallChD<-morphologyall[,c(1,42:44)]
morphologyallChD<-morphologyallChD[!(is.na(morphologyallChD[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallChD$meanChD <- 0
for (i in 1:nrow(morphologyallChD)) {
  morphologyallChD1 <- sort(morphologyallChD[i,2:4])
  morphologyallChD2<-morphologyallChD1[seq.int(which.min(diff(as.numeric(morphologyallChD1), lag = 1)), length.out = 2)]
  morphologyallChD3<-(morphologyallChD2[,1]+morphologyallChD2[,2])/2
  morphologyallChD$meanChD[i]<-morphologyallChD3
}
#POD
morphologyallPOD<-morphologyall[,c(1,46:48)]
morphologyallPOD<-morphologyallPOD[!(is.na(morphologyallPOD[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallPOD$meanPOD <- 0
for (i in 1:nrow(morphologyallPOD)) {
  morphologyallPOD1 <- sort(morphologyallPOD[i,2:4])
  morphologyallPOD2<-morphologyallPOD1[seq.int(which.min(diff(as.numeric(morphologyallPOD1), lag = 1)), length.out = 2)]
  morphologyallPOD3<-(morphologyallPOD2[,1]+morphologyallPOD2[,2])/2
  morphologyallPOD$meanPOD[i]<-morphologyallPOD3
}
#IOW
morphologyallIOW<-morphologyall[,c(1,50:52)]
morphologyallIOW<-morphologyallIOW[!(is.na(morphologyallIOW[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallIOW$meanIOW <- 0
for (i in 1:nrow(morphologyallIOW)) {
  morphologyallIOW1 <- sort(morphologyallIOW[i,2:4])
  morphologyallIOW2<-morphologyallIOW1[seq.int(which.min(diff(as.numeric(morphologyallIOW1), lag = 1)), length.out = 2)]
  morphologyallIOW3<-(morphologyallIOW2[,1]+morphologyallIOW2[,2])/2
  morphologyallIOW$meanIOW[i]<-morphologyallIOW3
}
#POW
morphologyallPOW<-morphologyall[,c(1,54:56)]
morphologyallPOW<-morphologyallPOW[!(is.na(morphologyallPOW[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallPOW$meanPOW <- 0
for (i in 1:nrow(morphologyallPOW)) {
  morphologyallPOW1 <- sort(morphologyallPOW[i,2:4])
  morphologyallPOW2<-morphologyallPOW1[seq.int(which.min(diff(as.numeric(morphologyallPOW1), lag = 1)), length.out = 2)]
  morphologyallPOW3<-(morphologyallPOW2[,1]+morphologyallPOW2[,2])/2
  morphologyallPOW$meanPOW[i]<-morphologyallPOW3
}
#CPL
morphologyallCPL<-morphologyall[,c(1,58:60)]
morphologyallCPL<-morphologyallCPL[!(is.na(morphologyallCPL[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallCPL$meanCPL <- 0
for (i in 1:nrow(morphologyallCPL)) {
  morphologyallCPL1 <- sort(morphologyallCPL[i,2:4])
  morphologyallCPL2<-morphologyallCPL1[seq.int(which.min(diff(as.numeric(morphologyallCPL1), lag = 1)), length.out = 2)]
  morphologyallCPL3<-(morphologyallCPL2[,1]+morphologyallCPL2[,2])/2
  morphologyallCPL$meanCPL[i]<-morphologyallCPL3
}
#CPD
morphologyallCPD<-morphologyall[,c(1,62:64)]
morphologyallCPD<-morphologyallCPD[!(is.na(morphologyallCPD[,2])),] # anything that has NA in first field has all NAs, so exlude these
morphologyallCPD$meanCPD <- 0
for (i in 1:nrow(morphologyallCPD)) {
  morphologyallCPD1 <- sort(morphologyallCPD[i,2:4])
  morphologyallCPD2<-morphologyallCPD1[seq.int(which.min(diff(as.numeric(morphologyallCPD1), lag = 1)), length.out = 2)]
  morphologyallCPD3<-(morphologyallCPD2[,1]+morphologyallCPD2[,2])/2
  morphologyallCPD$meanCPD[i]<-morphologyallCPD3
}

# now merge all means by id
allmeans<-Reduce(function(x,y) merge(x, y, by = "id", all.x = TRUE, all.y = TRUE), 
                 list(morphologyallSL[,c(1,5)], morphologyallBD[,c(1,5)], morphologyallHL[,c(1,5)], morphologyallHW[,c(1,5)], morphologyallLJL[,c(1,5)], morphologyallLJW[,c(1,5)],
                      morphologyallSnL[,c(1,5)], morphologyallSnW[,c(1,5)], morphologyallEyL[,c(1,5)], morphologyallEyD[,c(1,5)], morphologyallChD[,c(1,5)], morphologyallPOD[,c(1,5)],
                      morphologyallIOW[,c(1,5)], morphologyallPOW[,c(1,5)], morphologyallCPL[,c(1,5)], morphologyallCPD[,c(1,5)], morphologyall[,c(1,66,67)]))

# write out
write.csv(x = allmeans, file = "allPhenoMeans.csv")


# check if size correction is needed 
# also checking for outliers and weird distributions...

# subset to F2s
allmeansF2s<-allmeans[allmeans$Type=="F2",]
allmeansF2s<-droplevels(allmeansF2s)

# define colours 
factor(allmeansF2s$Sex)
cols2<-c("red","blue")[factor(allmeansF2s$Sex)]

# plot (log) traits against (log) SL
par(mfrow=c(4,4))
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanBD), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(BD)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanHL), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(HL)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanHW), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(HW)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanLJL), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(LJL)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanLJW), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(LJW)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanSnL), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(SnL)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanSnW), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(SnW)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanEyL), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(EyL)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanEyD), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(EyD)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanChD), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(ChD)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanPOD), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(POD)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanIOW), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(IOW)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanPOW), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(POW)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanCPL), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(CPL)")
  plot(log10(allmeansF2s$meanSL), log10(allmeansF2s$meanCPD), pch=19, col=cols2, xlab="log10(SL)", ylab="log10(CPD)")
  plot.new()
  legend("center", legend=c("female", "male"), col=c("red", "blue"), pch=19, bty = "o")

# perform correlation tests and homogeneity of slopes tests

cor.test(log10(allmeansF2s$meanBD), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanBD)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex
cor.test(log10(allmeansF2s$meanHL), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanHL)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex
cor.test(log10(allmeansF2s$meanHW), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanHW)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex
cor.test(log10(allmeansF2s$meanLJL), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanLJL)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex
cor.test(log10(allmeansF2s$meanLJW), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanLJW)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex
cor.test(log10(allmeansF2s$meanSnL), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanSnL)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex (only marginally significant)
cor.test(log10(allmeansF2s$meanSnW), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanSnW)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex
cor.test(log10(allmeansF2s$meanEyL), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanEyL)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex
cor.test(log10(allmeansF2s$meanEyD), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanEyD)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex
cor.test(log10(allmeansF2s$meanChD), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanChD)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex
cor.test(log10(allmeansF2s$meanPOD), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanPOD)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex
cor.test(log10(allmeansF2s$meanIOW), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanIOW)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex
cor.test(log10(allmeansF2s$meanPOW), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanPOW)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex
cor.test(log10(allmeansF2s$meanCPL), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanCPL)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated, no significant interaction of SL and Sex
cor.test(log10(allmeansF2s$meanCPD), log10(allmeansF2s$meanSL))
summary(aov(log10(allmeansF2s$meanCPD)~log10(allmeansF2s$meanSL)*allmeansF2s$Sex))
# highly correlated,  significant interaction of SL and Sex!!!

# only CPD has a significant interaction of Type with sex (and SnL marginally significant)


# size correction for external morphology
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

allmeansF2s$BDcorrF2<-residuals(lm(log10(allmeansF2s$meanBD)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$HLcorrF2<-residuals(lm(log10(allmeansF2s$meanHL)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$HWcorrF2<-residuals(lm(log10(allmeansF2s$meanHW)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$LJLcorrF2<-residuals(lm(log10(allmeansF2s$meanLJL)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$LJWcorrF2<-residuals(lm(log10(allmeansF2s$meanLJW)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$SnLcorrF2<-residuals(lm(log10(allmeansF2s$meanSnL)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$SnWcorrF2<-residuals(lm(log10(allmeansF2s$meanSnW)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$EyLcorrF2<-residuals(lm(log10(allmeansF2s$meanEyL)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$EyDcorrF2<-residuals(lm(log10(allmeansF2s$meanEyD)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$ChDcorrF2<-residuals(lm(log10(allmeansF2s$meanChD)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$PODcorrF2<-residuals(lm(log10(allmeansF2s$meanPOD)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$IOWcorrF2<-residuals(lm(log10(allmeansF2s$meanIOW)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$POWcorrF2<-residuals(lm(log10(allmeansF2s$meanPOW)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$CPLcorrF2<-residuals(lm(log10(allmeansF2s$meanCPL)~log10(allmeansF2s$meanSL), na.action=na.exclude))
allmeansF2s$CPDcorrF2<-residuals(lm(log10(allmeansF2s$meanCPD)~log10(allmeansF2s$meanSL), na.action=na.exclude))

# plot them
par(mfrow=c(4,4))
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$BDcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="BDcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$HLcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="HLcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$HWcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="HWcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$LJLcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="LJLcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$LJWcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="LJWcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$SnLcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="SnLcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$SnWcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="SnWcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$EyLcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="EyLcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$EyDcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="EyDcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$ChDcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="ChDcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$PODcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="PODcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$IOWcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="IOWcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$POWcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="POWcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$CPLcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="CPLcorr")
  plot(log10(allmeansF2s$meanSL), (allmeansF2s$CPDcorrF2), pch=19, col=cols2, xlab="log10(SL)", ylab="CPDcorr")
  plot.new()
  legend("center", legend=c("female", "male"), col=c("red", "blue"), pch=19, bty = "o")

# plot histograms for every trait within F2s (as in r/qtl but nicer) 
# all F2s
  fcol=names(allmeansF2s[,20:34])
  par(mfrow=c(4,4))
  for (i in 20:34){
    hist(x = allmeansF2s[,i], breaks = 20, xlab=fcol[i-19], main="", col="grey50", freq = F)
    lines(density(allmeansF2s[,i], na.rm=T))
  }

# only males
allmeansF2sMALES<-allmeansF2s[allmeansF2s$Sex=="male",]
  fcol=names(allmeansF2sMALES[,20:34])
  par(mfrow=c(4,4))
  for (i in 20:34){
    hist(x = allmeansF2sMALES[,i], breaks = 20, xlab=fcol[i-19], main="", col="grey50", freq = F)
    lines(density(allmeansF2sMALES[,i], na.rm=T))
  }



# combine all data for QTL mapping analyses 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# combine F2 data: colourF2s, teethF2s, gutSLF2s, and allmeansF2s
# by id, all=T
mappingdata1<-merge(teethF2s, gutsSLF2s, by="id", all=T)
mappingdata2<-merge(mappingdata1, allmeansF2s, by="id", all=T)
mappingdata3<-merge(mappingdata2, colourF2s, by="id", all=T)
# some checks 
#mappingdata3$sex.x==mappingdata3$sex.y
#mappingdata3$meanSL.x==mappingdata3$meanSL.y
# looks ok
# remove doubles
mappingdata4<-mappingdata3[,-c(7,11,12,14,15,17:34,50,51)]
# re-arrange columns
mappingdata5<-mappingdata4[,c(1,6,27,12:26,11,10,7,8,2:5,9,28:43)]
# females are currently missing family info 
# read in ids
infos<-read.table("PNPxNOM_7more_ids.txt", header = T)
mappingdata6<-merge(mappingdata5, infos, by="id", all=F)
mappingdataOUT<-mappingdata6[,c(1,2,45,4:43)]
# rename colnames
names(mappingdataOUT)
colnames(mappingdataOUT)<-c("id", "sex", "family", "BD", "HL", "HW", "LJL", "LJW", "SnL", "SnW", "EyL", "EyD", "ChD", "POD", "IOW",  "POW", "CPL",  "CPD", 
                            "gutlength", "fattytissue", "toothshapeL",  "toothShapeU", "toothDensityL" , "toothDensityU",   
                            "Gap", "InnerRows", "tricuspid(tg)", "Rdorsalfin1", "Rdorsalfin2", "Rhead", "Rdorsum1",  "Rdorsum2",      
                             "Rflanks1(tg)", "Rcheek(tg)", "Rgillcover(tg)","Rnose","Yflanks", "Yupperlip(tg)", "Ycheek(tg)","Ygillcover",   
                             "Ypelvicfin(tg)","Ynose","No_stripes")

# write out
write.csv(x = mappingdataOUT, file = "allPhenos.csv")
