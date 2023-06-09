setwd("~/Dropbox/Things we generated during Shuxi's visit/R code/survival_plots")
library(survival)
library(ggplot2)
library(ggfortify)
##Survival Curve Drawn
rm(list=ls())
# setwd("C:/Users/Shuxi ZENG/Dropbox/Third Year/FPCA_New/FPCA_1005")
load("ELADSIGC data for Shuxi 9.22.2019.RData")

f = 1
age.truncate.female = 18
source("Collapse_Data.R")
max_age=aggregate(gender.data$age_sample,by=list(gender.data$sname),max)$x
treat_ind=aggregate(gender.data$adv_cumulative,by=list(gender.data$sname),max)$x
coarse_summary=treat_ind
coarse_summary[treat_ind>=2]=2
coarse_summary=as.factor(coarse_summary)

survival.data=data.frame(age=max_age,group=coarse_summary,state=1)
km_fit <- survfit(Surv(age, state) ~ 1, data=survival.data)

#plot(km_fit, xlab="Age in year",xlim=c(3.9,25),main = 'Kaplan Meyer Plot, Pooled Sample') #base graphics is always ready

pdf("../../Figures/Survivalplots/KMplot_pool.pdf")
autoplot(km_fit)+xlim(c(3.9,25))+ggtitle("Kaplan Meyer Plot, Pooled Sample")+
  ylab("Survival rate")+xlab("Age in year")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

km_fit_adv <- survfit(Surv(age, state) ~ coarse_summary, data=survival.data)

pdf("../../Figures/Survivalplots/KMplot_group.pdf")
autoplot(km_fit_adv)+xlim(c(3.9,25))+ggtitle("Kaplan Meyer Plot, by Adversity")+
  ylab("Survival rate")+xlab("Age in year")+ guides(fill = FALSE)+
  labs(colour = "Number of Adv") +
  scale_color_manual(labels = c("0", "1","2+"),values=c(2,3,4))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()