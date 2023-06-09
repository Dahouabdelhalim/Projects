#Supplementary figures for Arctic Science manuscript (simple descriptive bar plots)
#Ariane Castagner, finalized June 2022

#Import package(s)
library(ggplot2)
library(plyr)

#Set wd (please change the path to your local directory for the project)

setwd("...\\\\ExcessIceBetaRegression\\\\")

#Source functions

source("src//matthickness.R")
source("src//VIthickness.R")
source("src//surficialmat.R")
source("src//bhdepths.R")

#Import data
data <- read.csv("data//CryostratigraphyBHdataMI.csv")

#Part 1 - Bar plot showing frequency of materials in all boreholes
  #Make frequency table based on the thickness of materials

matthickness <- matthickness(data)
matthickness$Perc.thickness <- as.numeric(format(round(matthickness$Perc.thickness, 1), nsmall=1))
matthickness$Material <- factor(matthickness$Material, levels = matthickness$Material[order(matthickness$Perc.thickness, decreasing = TRUE)])

  #Plot
p1 <- ggplot(data=matthickness, aes(x=Material, y=Perc.thickness)) +
  geom_bar(stat="identity", color="black", fill = "cadetblue") +
  geom_text(aes(label=Perc.thickness), vjust=-0.5, color="black", size=4) +
  theme(text=element_text(size=15)) +
  ylab("Percentage of total thickness (%)") +
  scale_y_continuous(breaks=c(0, 5, 10, 15, 20, 25, 30, 35)) +
  theme_minimal()

ggsave(filename = "output//figs//BarplotsP1.jpg", device="jpg", plot=p1, scale = 1, dpi = 300)

#Part 2 - Bar plot showing frequency of visible ice classes in all boreholes

VIthickness <- VIthickness(data)
VIthickness$Perc.thickness <- as.numeric(format(round(VIthickness$Perc.thickness, 1), nsmall=1))
VIthickness$Visible.ice <- factor(VIthickness$Visible.ice, levels = VIthickness$Visible.ice[order(VIthickness$Perc.thickness, decreasing = TRUE)])
VIthickness <- na.omit(VIthickness)

#Plot
p2 <- ggplot(data=VIthickness, aes(x=Visible.ice, y=Perc.thickness)) +
  geom_bar(stat="identity", color="black", fill = "cadetblue") +
  geom_text(aes(label=Perc.thickness), vjust=-0.5, color="black", size=4) +
  theme(text=element_text(size=15)) +
  ylab("Percentage of total thickness (%)") +
  xlab("Visible ice") +
  theme_minimal()

ggsave(filename = "output//figs//BarplotsP2.jpg", device="jpg", plot=p2, scale = 1, dpi = 300)

#Part 3 - Bar plot showing frequency of surficial materials

surficialmat <- surficialmat(data)
surficialmat$Perc <- as.numeric(format(round(surficialmat$Perc, 1), nsmall=1))
surficialmat <- surficialmat[which(surficialmat$Perc > 0),]
surficialmat$Var1 <- factor(surficialmat$Var1, levels = surficialmat$Var1[order(surficialmat$Perc, decreasing = TRUE)])

#Plot
p3 <- ggplot(data=surficialmat, aes(x=Var1, y=Perc)) +
  geom_bar(stat="identity", color="black", fill = "cadetblue") +
  geom_text(aes(label=Perc), vjust=-0.5, color="black", size=4) +
  theme(text=element_text(size=15)) +
  ylab("Percentage of boreholes (%)") +
  xlab("Material") +
  theme_minimal()

ggsave(filename = "output//figs//BarplotsP3.jpg", device="jpg", plot=p3, scale = 1, dpi = 300)

#Part 4 - Bar plot showing the frequency of borehole depths

bhdepths <- bhdepths(data)
bhdepths$Perc <- as.numeric(format(round(bhdepths$Perc, 1), nsmall=1))
bhdepths <- bhdepths[which(bhdepths$Perc > 0),]
bhdepths$Depths <- factor(bhdepths$Depths, levels = bhdepths$Depths[order(bhdepths$Perc, decreasing = TRUE)])

#Plot
p4 <- ggplot(data=bhdepths, aes(x=Depths, y=Perc)) +
  geom_bar(stat="identity", color="black", fill = "cadetblue") +
  geom_text(aes(label=Perc), vjust=-0.5, color="black", size=4) +
  theme(text=element_text(size=15)) +
  ylab("Percentage of boreholes (%)") +
  xlab("Depth (m)") +
  theme_minimal()

p4_2 <- ggplot(data=bhdepths, aes(x=Depths, y=Perc, group=1)) +
  geom_line(size=1) +
  geom_point(size=3, shape=18) +
  theme(text=element_text(size=16)) +
  xlab("Depth (m)") +
  ylab("Boreholes (%) deeper than") +
  scale_y_continuous(breaks=seq(30, 100, 10)) +
  scale_x_discrete(limits=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) +
  theme_minimal()

ggsave(filename = "output//figs//BarplotsP4.jpg", device="jpg", plot=p4, scale = 1, dpi = 300)
ggsave(filename = "output//figs//BarplotsP4_2.jpg", device="jpg", plot=p4_2, scale = 1, dpi = 300)