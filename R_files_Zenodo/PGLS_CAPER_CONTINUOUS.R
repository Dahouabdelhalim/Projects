##First, we need to load the packages to read the phylogeny & run the analyses.
library(ape)
library(nlme)
library(geiger)
library(caper)

##load the data & tree
phy <- read.tree("C_tree.nwk")
data <-read.csv("Continuous.csv",header=TRUE)

##Now the data & the tree should coincide in both number of species and names. We can verify this using name.check.
name.check(phy,data)


##defining the quantitative characters
#Lablength
Lablength <- data$Lablength
names(Lablength) <- data$Species
class(Lablength)
Lablength <- as.factor(Lablength)
class(Lablength)
Lablength

#Labwidth
Labwidth <- data$Labwidth
names(Labwidth) <- data$Species
class(Labwidth)
Labwidth <- as.factor(Labwidth)
class(Labwidth)
Labwidth

#Labnotch
Labnotch <- data$Labnotch
names(Labnotch) <- data$Species
class(Labnotch)
Labnotch <- as.factor(Labnotch)
class(Labnotch)
Labnotch

#Laterallength
Laterallength <- data$Laterallength
names(Laterallength) <- data$Species
class(Laterallength)
Laterallength <- as.factor(Laterallength)
class(Laterallength)
Laterallength

#Lateralwidth
Lateralwidth <- data$Lateralwidth
names(Lateralwidth) <- data$Species
class(Lateralwidth)
Lateralwidth <- as.factor(Lateralwidth)
class(Lateralwidth)
Lateralwidth

#FTLength
FTLength <- data$FTLength
names(FTLength) <- data$Species
class(FTLength)
FTLength <- as.factor(FTLength)
class(FTLength)
FTLength

#Fillength
Fillength <- data$Fillength
names(Fillength) <- data$Species
class(Fillength)
Fillength <- as.factor(Fillength)
class(Fillength)
Fillength

#Antherlength
Antherlength <- data$Antherlength
names(Antherlength) <- data$Species
class(Antherlength)
Antherlength <- as.factor(Antherlength)
class(Antherlength)
Antherlength

#For combining phylogenies with datasets and ensure consistent structure and ordering for use in functions.
comp.data<-comparative.data(phy, data, names.col="Species", vcv.dim=2, warn.dropped=TRUE)

#1. CLLvsCLW
CLLvsCLW <-pgls(Lablength~Labwidth, data=comp.data, lambda="ML")
summary(CLLvsCLW)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLLvsCLW.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLLvsCLW<-pgls.profile(CLLvsCLW, which="lambda")
plot(lm.lk_CLLvsCLW)
# 3. Close the file
dev.off()


#2. CLLvsCLND
CLLvsCLND <-pgls(Lablength~Labnotch, data=comp.data, lambda="ML")
summary(CLLvsCLND)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLLvsCLND.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLLvsCLND<-pgls.profile(CLLvsCLND, which="lambda")
plot(lm.lk_CLLvsCLND)
# 3. Close the file
dev.off()


#3. CLLvsLSL
CLLvsLSL <-pgls(Lablength~Laterallength, data=comp.data, lambda="ML")
summary(CLLvsLSL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLLvsLSL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLLvsLSL<-pgls.profile(CLLvsLSL, which="lambda")
plot(lm.lk_CLLvsLSL)
# 3. Close the file
dev.off()


#4. CLLvsLSW
CLLvsLSW <-pgls(Lablength~Lateralwidth, data=comp.data, lambda="ML")
summary(CLLvsLSW)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLLvsLSW.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLLvsLSW<-pgls.profile(CLLvsLSW, which="lambda")
plot(lm.lk_CLLvsLSW)
# 3. Close the file
dev.off()


#5. CLLvsFTL
CLLvsFTL <-pgls(Lablength~FTLength, data=comp.data, lambda="ML")
summary(CLLvsFTL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLLvsFTL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLLvsFTL<-pgls.profile(CLLvsFTL, which="lambda")
plot(lm.lk_CLLvsFTL)
# 3. Close the file
dev.off()


#6. CLLvsFL
CLLvsFL <-pgls(Lablength~Fillength, data=comp.data, lambda="ML")
summary(CLLvsFL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLLvsFL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLLvsFL<-pgls.profile(CLLvsFL, which="lambda")
plot(lm.lk_CLLvsFL)
# 3. Close the file
dev.off()
 
       
#7. CLLvsAL
CLLvsAL <-pgls(Lablength~Antherlength, data=comp.data, lambda="ML")
summary(CLLvsAL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLLvsAL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLLvsAL<-pgls.profile(CLLvsAL, which="lambda")
plot(lm.lk_CLLvsAL)
# 3. Close the file
dev.off()


#8. CLWvsCLND
CLWvsCLND <-pgls(Labwidth~Labnotch, data=comp.data, lambda="ML")
summary(CLWvsCLND)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLWvsCLND.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLWvsCLND<-pgls.profile(CLWvsCLND, which="lambda")
plot(lm.lk_CLWvsCLND)
# 3. Close the file
dev.off()


#9. CLWvsLSL
CLWvsLSL <-pgls(Labwidth~Laterallength, data=comp.data, lambda="ML")
summary(CLWvsLSL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLWvsLSL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLWvsLSL<-pgls.profile(CLWvsLSL, which="lambda")
plot(lm.lk_CLWvsLSL)
# 3. Close the file
dev.off()


#10. CLWvsLSW
CLWvsLSW <-pgls(Labwidth~Lateralwidth, data=comp.data, lambda="ML")
summary(CLWvsLSW)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLWvsLSW.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLWvsLSW<-pgls.profile(CLWvsLSW, which="lambda")
plot(lm.lk_CLWvsLSW)
# 3. Close the file
dev.off()


#11. CLWvsFTL
CLWvsFTL <-pgls(Labwidth~FTLength, data=comp.data, lambda="ML")
summary(CLWvsFTL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLWvsFTL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLWvsFTL<-pgls.profile(CLWvsFTL, which="lambda")
plot(lm.lk_CLWvsFTL)
# 3. Close the file
dev.off()


#12. CLWvsFL
CLWvsFL <-pgls(Labwidth~Fillength, data=comp.data, lambda="ML")
summary(CLWvsFL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLWvsFL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLWvsFL<-pgls.profile(CLWvsFL, which="lambda")
plot(lm.lk_CLWvsFL)
# 3. Close the file
dev.off()


#13. CLWvsAL
CLWvsAL <-pgls(Labwidth~Antherlength, data=comp.data, lambda="ML")
summary(CLWvsAL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLWvsAL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLWvsAL<-pgls.profile(CLWvsAL, which="lambda")
plot(lm.lk_CLWvsAL)
# 3. Close the file
dev.off()


#14. CLNDvsLSL
CLNDvsLSL <-pgls(Labnotch~Laterallength, data=comp.data, lambda="ML")
summary(CLNDvsLSL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLNDvsLSL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLNDvsLSL<-pgls.profile(CLNDvsLSL, which="lambda")
plot(lm.lk_CLNDvsLSL)
# 3. Close the file
dev.off()


#15. CLNDvsLSW
CLNDvsLSW <-pgls(Labnotch~Lateralwidth, data=comp.data, lambda="ML")
summary(CLNDvsLSW)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLNDvsLSW.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLNDvsLSW<-pgls.profile(CLNDvsLSW, which="lambda")
plot(lm.lk_CLNDvsLSW)
# 3. Close the file
dev.off()


#16. CLNDvsFTL
CLNDvsFTL <-pgls(Labnotch~FTLength, data=comp.data, lambda="ML")
summary(CLNDvsFTL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLNDvsFTL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLNDvsFTL<-pgls.profile(CLNDvsFTL, which="lambda")
plot(lm.lk_CLNDvsFTL)
# 3. Close the file
dev.off()


#17. CLNDvsFL
CLNDvsFL <-pgls(Labnotch~Fillength, data=comp.data, lambda="ML")
summary(CLNDvsFTL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLNDvsFL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLNDvsFL<-pgls.profile(CLNDvsFL, which="lambda")
plot(lm.lk_CLNDvsFL)
# 3. Close the file
dev.off()


#18. CLNDvsAL
CLNDvsAL <-pgls(Labnotch~Antherlength, data=comp.data, lambda="ML")
summary(CLNDvsFTL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLNDvsAL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLNDvsAL<-pgls.profile(CLNDvsAL, which="lambda")
plot(lm.lk_CLNDvsAL)
# 3. Close the file
dev.off()


#19. LSLvsLSW
LSLvsLSW <-pgls(Laterallength~Lateralwidth, data=comp.data, lambda="ML")
summary(LSLvsLSW)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_LSLvsLSW.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_LSLvsLSW<-pgls.profile(LSLvsLSW, which="lambda")
plot(lm.lk_LSLvsLSW)
# 3. Close the file
dev.off()


#20. LSLvsFTL
LSLvsFTL <-pgls(Laterallength~FTLength, data=comp.data, lambda="ML")
summary(LSLvsFTL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_LSLvsFTL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_LSLvsFTL<-pgls.profile(LSLvsFTL, which="lambda")
plot(lm.lk_LSLvsFTL)
# 3. Close the file
dev.off()


#21. LSLvsFL
LSLvsFL <-pgls(Laterallength~Fillength, data=comp.data, lambda="ML")
summary(LSLvsFL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_ LSLvsFL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_LSLvsFL<-pgls.profile(LSLvsFL, which="lambda")
plot(lm.lk_LSLvsFL)
# 3. Close the file
dev.off()


#22. LSLvsAL
LSLvsAL <-pgls(Laterallength~Antherlength, data=comp.data, lambda="ML")
summary(LSLvsAL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_LSLvsAL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_LSLvsAL<-pgls.profile(LSLvsAL, which="lambda")
plot(lm.lk_LSLvsAL)
# 3. Close the file
dev.off()


#23. LSWvsFTL
LSWvsFTL <-pgls(Lateralwidth~FTLength, data=comp.data, lambda="ML")
summary(LSWvsFTL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_LSWvsFTL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_LSWvsFTL<-pgls.profile(LSWvsFTL, which="lambda")
plot(lm.lk_LSWvsFTL)
# 3. Close the file
dev.off()


#24. LSWvsFL
LSWvsFL <-pgls(Lateralwidth~Fillength, data=comp.data, lambda="ML")
summary(LSWvsFL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_LSWvsFL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_LSWvsFL<-pgls.profile(LSWvsFL, which="lambda")
plot(lm.lk_LSWvsFL)
# 3. Close the file
dev.off()


#25. LSWvsAL
LSWvsAL <-pgls(Lateralwidth~Antherlength, data=comp.data, lambda="ML")
summary(LSWvsAL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_LSWvsAL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_LSWvsAL<-pgls.profile(LSWvsAL, which="lambda")
plot(lm.lk_LSWvsAL)
# 3. Close the file
dev.off()


#26. FTLvsFL
FTLvsFL <-pgls(FTLength~Fillength, data=comp.data, lambda="ML")
summary(FTLvsFL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_FTLvsFL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_FTLvsFL<-pgls.profile(FTLvsFL, which="lambda")
plot(lm.lk_FTLvsFL)
# 3. Close the file
dev.off()


#27. FTLvsAL
FTLvsAL <-pgls(FTLength~Antherlength, data=comp.data, lambda="ML")
summary(FTLvsAL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_FTLvsAL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_FTLvsAL<-pgls.profile(FTLvsAL, which="lambda")
plot(lm.lk_FTLvsAL)
# 3. Close the file
dev.off()


#28. FLvsAL
FLvsAL <-pgls(Fillength~Antherlength, data=comp.data, lambda="ML")
summary(FLvsAL)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_FLvsAL.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_FLvsAL<-pgls.profile(FLvsAL, which="lambda")
plot(lm.lk_FLvsAL)
# 3. Close the file
dev.off()


#29. CLLvsCLW+CLND+LSL+LSW+FTL+FL+AL
CLLvsRest <-pgls(Lablength~Labwidth+Labnotch+Laterallength+Lateralwidth+FTLength+Fillength+Antherlength, data=comp.data, lambda="ML")
summary(CLLvsRest)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLLvsRest.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLLvsRest<-pgls.profile(CLLvsRest, which="lambda")
plot(lm.lk_CLLvsRest)
# 3. Close the file
dev.off()


#30. CLWvsCLL+CLND+LSL+LSW+FTL+FL+AL
CLWvsRest <-pgls(Labwidth~Lablength+Labnotch+Laterallength+Lateralwidth+FTLength+Fillength+Antherlength, data=comp.data, lambda="ML")
summary(CLWvsRest)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLWvsRest.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLWvsRest<-pgls.profile(CLWvsRest, which="lambda")
plot(lm.lk_CLWvsRest)
# 3. Close the file
dev.off()


#31. CLNDvsCLL+CLW+LSL+LSW+FTL+FL+AL
CLNDvsRest <-pgls(Labnotch~Lablength+Labwidth+Laterallength+Lateralwidth+FTLength+Fillength+Antherlength, data=comp.data, lambda="ML")
summary(CLNDvsRest)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_CLNDvsRest.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_CLNDvsRest<-pgls.profile(CLNDvsRest, which="lambda")
plot(lm.lk_CLNDvsRest)
# 3. Close the file
dev.off()


#32. LSLvsCLL+CLW+CLND+LSW+FTL+FL+AL
LSLvsRest <-pgls(Laterallength~Lablength+Labwidth+Labnotch+Lateralwidth+FTLength+Fillength+Antherlength, data=comp.data, lambda="ML")
summary(LSLvsRest)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_LSLvsRest.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_LSLvsRest<-pgls.profile(LSLvsRest, which="lambda")
plot(lm.lk_LSLvsRest)
# 3. Close the file
dev.off()


#33. LSWvsCLL+CLW+CLND+LSL+LSW+FTL+FL+AL
LSWvsRest <-pgls(Lateralwidth~Lablength+Labwidth+Labnotch+Laterallength+FTLength+Fillength+Antherlength, data=comp.data, lambda="ML")
summary(LSWvsRest)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_LSWvsRest.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_LSWvsRest<-pgls.profile(LSWvsRest, which="lambda")
plot(lm.lk_LSWvsRest)
# 3. Close the file
dev.off()


#34. FTLvsCLL+CLW+CLND+LSL+LSW+FL+AL
FTLvsRest <-pgls(FTLength~Lablength+Labwidth+Labnotch+Laterallength+Labwidth+Fillength+Antherlength, data=comp.data, lambda="ML")
summary(FTLvsRest)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_FTLvsRest.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_FTLvsRest<-pgls.profile(FTLvsRest, which="lambda")
plot(lm.lk_FTLvsRest)
# 3. Close the file
dev.off()


#35. FLvsCLL+CLW+CLND+LSL+LSW+FTL+AL
FLvsRest <-pgls(Fillength~Lablength+Labwidth+Labnotch+Laterallength+Labwidth+FTLength+Antherlength, data=comp.data, lambda="ML")
summary(FLvsRest)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_FLvsRest.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_FLvsRest<-pgls.profile(FLvsRest, which="lambda")
plot(lm.lk_FLvsRest)
# 3. Close the file
dev.off()


#36. ALvsCLL+CLW+CLND+LSL+LSW+FTL+FL
ALvsRest <-pgls(Antherlength~Lablength+Labwidth+Labnotch+Laterallength+Labwidth+FTLength+Fillength, data=comp.data, lambda="ML")
summary(ALvsRest)
##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("lm.lk_ALvsRest .pdf", width = 8.27, height = 11.69)
# 2. Create the plot
lm.lk_ALvsRest <-pgls.profile(ALvsRest, which="lambda")
plot(lm.lk_ALvsRest)
# 3. Close the file
dev.off()