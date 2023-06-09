library(phytools)
library(ape)
phy <- read.tree("C_tree.nwk")
dat <- read.csv("Continuous.csv", header = T, row.names = 1)

#checking the names of tree tips and data file if they are matching.
name.check(phy,dat)

#Labellum length
# 1. Open pdf file
pdf("Lablength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
Lablength <- dat$Lablength
names(Lablength) <- rownames(dat)
obj<-contMap(phy,Lablength,plot=FALSE)
plot(obj,type="phylogram",legend=0.7*max(nodeHeights(phy)),
     fsize=c(0.7,0.9), outline = F)
# 3. Close the file
dev.off()

#Labellum width
# 1. Open pdf file
pdf("Labwidth.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
Labwidth <- dat$Labwidth
names(Labwidth) <- rownames(dat)
obj<-contMap(phy,Labwidth,plot=FALSE)
plot(obj,type="phylogram",legend=0.7*max(nodeHeights(phy)),
     fsize=c(0.7,0.9), outline = F)
# 3. Close the file
dev.off()

#Labellum notch
# 1. Open pdf file
pdf("Labnotch.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
Labnotch <- dat$Labnotch
names(Labnotch) <- rownames(dat)
obj<-contMap(phy,Labnotch,plot=FALSE)
plot(obj,type="phylogram",legend=0.7*max(nodeHeights(phy)),
     fsize=c(0.7,0.9), outline = F)
# 3. Close the file
dev.off()

#Lateral staminode length
# 1. Open pdf file
pdf("Labterallength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
Laterallength <- dat$Laterallength
names(Laterallength) <- rownames(dat)
obj<-contMap(phy,Laterallength,plot=FALSE)
plot(obj,type="phylogram",legend=0.7*max(nodeHeights(phy)),
     fsize=c(0.7,0.9), outline = F)
# 3. Close the file
dev.off()

#Lateral staminode width
# 1. Open pdf file
pdf("Labteralwidth.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
Lateralwidth <- dat$Lateralwidth
names(Lateralwidth) <- rownames(dat)
obj<-contMap(phy,Lateralwidth,plot=FALSE)
plot(obj,type="phylogram",legend=0.7*max(nodeHeights(phy)),
     fsize=c(0.7,0.9), outline = F)
# 3. Close the file
dev.off()

#Floral tube length
# 1. Open pdf file
pdf("FTLength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
FTLength <- dat$FTLength
names(FTLength) <- rownames(dat)
obj<-contMap(phy,FTLength,plot=FALSE)
plot(obj,type="phylogram",legend=0.7*max(nodeHeights(phy)),
     fsize=c(0.7,0.9), outline = F)
# 3. Close the file
dev.off()

#Filament length
# 1. Open pdf file
pdf("Fillength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
Fillength <- dat$Fillength
names(Fillength) <- rownames(dat)
obj<-contMap(phy,Fillength,plot=FALSE)
plot(obj,type="phylogram",legend=0.7*max(nodeHeights(phy)),
     fsize=c(0.7,0.9), outline = F)
# 3. Close the file
dev.off()

#Anther length
# 1. Open pdf file
pdf("Antherlength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
Antherlength <- dat$Antherlength
names(Antherlength) <- rownames(dat)
obj<-contMap(phy,Antherlength,plot=FALSE)
plot(obj,type="phylogram",legend=0.7*max(nodeHeights(phy)),
     fsize=c(0.7,0.9), outline = F)
# 3. Close the file
dev.off()