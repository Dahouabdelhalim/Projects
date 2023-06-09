rm(list=ls())
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library("survminer")
library(ggpubr)

##  Daten einlesen
setwd("D:/ConsumerDemandPaper/MaxPaidPrice")
MPP <- read.table ('MaxPaidPrice.txt', header= TRUE, sep= '\\t', dec= ',', as.is= FALSE)


str(MPP)
attach(MPP)
names(MPP)
MPP$Run<-as.factor(MPP$Run)


MaxPaidPrice <- survfit(Surv(MaxNumberOfNP) ~ Run, cluster=animal, data=MPP)
print(MaxPaidPrice)
summary(MaxPaidPrice,time=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,4041,42,43,444,45,46,47,48,49,5051,52,53,54,55,56,57,58,59,60,61,62,63,64))

autoplot(MaxPaidPrice)
MaxPaidPrice

PaidPrice<-ggsurvplot(MaxPaidPrice,
          pval = F, conf.int = T, legend = c("right"),
          xlim=c(1,65),
          break.time.by = 4, # break time axis by 1
          palette = c("#999900","#33FF00", "#339900","#CCCC00","#3399CC","#00CCFF","#0000CC","#0099FF"),
          ylab="Animals [%]", xlab="Maximum Price Paid",
          risk.table = TRUE, # Add risk table
          risk.table.col = "Run", risk.table.title="Animals", 
           legend.title="", tables.height = 0.4,
          surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw(), # Change ggplot2 theme
          ) 
PaidPrice
ggpar(PaidPrice, 
      font.main = c(18, "bold"),
      font.x = c(18, "plain"),
      font.y = c(18, "bold"),
      font.caption = c(18, "bold"), 
      font.legend = c(18, "bold"), 
      font.tickslab = c(15, "bold"))