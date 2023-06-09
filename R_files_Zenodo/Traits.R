library(ggplot2)
library(dplyr)
library(reshape)
library(RColorBrewer)
library(Hmisc)
library(PerformanceAnalytics)
library(corrplot)


setwd("working directory")

traits<- read.table("phenotypes_200703.txt", header = TRUE)

#reconfigure matrix

mtraits<- melt(traits)

#plot trait histograms with vertical lines for parentals and F1 means
  str(mtraits)
  head(mtraits)
  parent.f1.avg <- mtraits %>%
    filter(Type == "P" | Type == "S" | Type == "F1") %>%
    group_by(Type, variable) %>%
    summarise( mean.value = mean(value))

  pdf("histograms.pdf", width=7.08,height=5.5)  
  ggplot(mtraits, aes(value)) + theme_bw(base_size = 10) +
    geom_histogram(data = filter(mtraits, Type == "F2")) +
    geom_vline(data = filter(parent.f1.avg,Type == "F1"), aes(xintercept=mean.value), color="#fd8d3c", linetype="dashed", size = 1) +
    geom_vline(data = filter(parent.f1.avg,Type == "P"), aes(xintercept=mean.value), color="#e31a1c", linetype="longdash", size = 1) +
    geom_vline(data = filter(parent.f1.avg,Type == "S"), aes(xintercept=mean.value), color="#fecc5c", linetype="dotted", size = 1) +
    facet_wrap(~ variable, scales="free")
  dev.off()
  
#analyze trait correlations

#subset F2 data and remove ID
F2traits<-subset(traits,subset=traits$Type=="F2", select = -c(Type,ID))  

correlations <- rcorr(as.matrix(F2traits), type = "spearman")

chart.Correlation(as.matrix(F2traits), histogram=FALSE, pch=1)

pdf("correlations.pdf", width = 3.14, height = 3.14)

corrplot(correlations$r, type="lower", order="original", 
         p.mat = correlations$P, sig.level = 0.000326797, insig = "blank", tl.col="black", cl.cex = 0.7, tl.cex = 0.8)

dev.off()
