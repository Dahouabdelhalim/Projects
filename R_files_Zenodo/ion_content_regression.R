setwd("D:/oa_ma") #set working directory
setwd("/Volumes/NO NAME/oa_ma")
getwd() #check working directory

#load relevant packages for analysis and visualization
library(dplyr)
library(ggplot2)

data <- read.csv("ca_mg_regression.csv", header = TRUE) #load data into 'data' variable
data$ï..ReferenceShort #check all data points are accounted for - 41 data points from 17 unique articles

#filter ratio data into variables depending on high CO2 range
atm.500 <- filter(data, data$pco2.bin == "500-999")
atm.1000 <- filter(data, pco2.bin == "1000-1499")
atm.1500 <- filter(data, pco2.bin == "1500-1999")
atm.2000 <- filter(data, pco2.bin == "2000+")

#incrmenetally plot points from each bins as well as regression lines
ggplot(atm.500, aes(x = calc.lnrr, y = magn.lnrr))+
  geom_point(aes(color = "#000000"), size = 2)+
  geom_abline(color = "#000000", slope = coef(fit.500)[[2]], intercept = coef(fit.500)[[1]], size = 1.5)+
  geom_point(data = atm.1000, aes(x = calc.lnrr, y = magn.lnrr, color = "#E69F00"), size = 2)+
  geom_abline(color = "#E69F00", slope = coef(fit.1000)[[2]], intercept = coef(fit.1000)[[1]], size = 1.5)+
  geom_point(data = atm.1500, aes(x = calc.lnrr, y = magn.lnrr, color = "#56B4E9"), size = 2)+
  geom_abline(color = "#56B4E9", slope = coef(fit.1500)[[2]], intercept = coef(fit.1500)[[1]], size = 1.5)+
  geom_point(data = atm.2000, aes(x = calc.lnrr, y = magn.lnrr, color = "#009E73"), size = 2)+
  geom_abline(color = "#009E73", slope = coef(fit.2000)[[2]], intercept = coef(fit.2000)[[1]], size = 1.5)+
  xlab("Ca effect size (lnRR)")+
  ylab("Mg effect size (lnRR)")+
  scale_fill_discrete(labels = c("500-999", "1000-1499", "1500-1999", "2000+"))+
  scale_color_manual(name = "pCO2 bin (\\u03BCatm)",labels = c("500-999", "1000-1499", "1500-1999", "2000+"),
                     values = c("#000000", "#E69F00", "#56B4E9", "#009E73"))+
  theme(axis.title=element_text(size = 15), axis.text=element_text(size=15), legend.text=element_text(size=15))+
  theme_classic()

#summary numbers for the linear regressions for each pCO2 bins
summary(lm(atm.500$magn.lnrr ~ atm.500$calc.lnrr))
summary(lm(atm.1000$magn.lnrr ~ atm.1000$calc.lnrr))
summary(lm(atm.1500$magn.lnrr ~ atm.1500$calc.lnrr))
summary(lm(atm.2000$magn.lnrr ~ atm.2000$calc.lnrr))
