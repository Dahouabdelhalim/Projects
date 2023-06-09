# #----------------------------------------------------
# This code produces graphs for the paper:
# CieÅ›lik, A., Gurshev, O. Friends with or without benefits? 
# An empirical evaluation of bilateral trade and economic integration
# between some of the post-Soviet economies. Eurasian Econ Rev (2022). 
# https://doi.org/10.1007/s40822-022-00213-9
#
# Date revised: 10.07.2022
#----------------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(forcats)
library(wesanderson)
#----------------------------------------------------
# directory 
setwd()
getwd()
#----------------------------------------------------
# Loading data

data <-
  read.table('cu.csv',
             sep = ",",
             dec = ".",
             header = TRUE)


exp <- 
  read.table('exp.csv',
             sep = ",",
             dec = ".",
             header = TRUE)

imp <- 
  read.table('imp.csv',
             sep = ",",
             dec = ".",
             header = TRUE)

data$logcu <- log(data$cu*1000000000)

#---------------------------------------------------
# This code uses cu.csv and produces Figure A1: Russia and CU

ggplot(data, aes(x=year)) + 
  geom_line(aes(y = logcu, color ="CU/EUEA"), size = 2) +
  geom_line(aes(y = Russia, color ="Russia"), size = 2) +
  labs (x ="year", y = "Log of total trade", color = 'Entity') +
  geom_vline(xintercept = 2010 , color = "red", linetype = "dashed") + 
  geom_vline(xintercept = 2012 , color = "red", linetype = "dashed") + 
  geom_vline(xintercept = 2015 , color = "red", linetype = "dashed") +
  geom_vline(xintercept = 2016 , color = "red", linetype = "dashed") +
  theme_classic() + 
  scale_x_continuous(breaks = scales::pretty_breaks(10)) + scale_color_brewer(palette="Paired")
#---------------------------------------------------
# This code uses cu.csv and produces Figure A2: total trade

ggplot(data, aes(x=year)) + 
  geom_line(aes(y = logcu, color = "CU/EUEA"), size = 2, linetype = 4) +
  geom_line(aes(y = Belarus, color = "Belarus"), size = 1) +
  geom_line(aes(y = Kazakhstan, color = "Kazakhstan"), size = 1) +
  geom_line(aes(y = Russia, color = "Russia"), size = 1) +
  geom_line(aes(y = Kyrgyzstan, color = "Kyrgyzstan"), size = 1) +
  geom_line(aes(y = Armenia, color = "Armenia"), size = 1) +
  labs (x ="year", y = "Log of total trade", color = 'Entity') +
  geom_vline(xintercept = 2010 , color = "red", linetype = "dashed") + 
  geom_vline(xintercept = 2012 , color = "red", linetype = "dashed") + 
  geom_vline(xintercept = 2015 , color = "red", linetype = "dashed") +
  geom_vline(xintercept = 2016 , color = "red", linetype = "dashed") +
  theme_classic() + 
  scale_x_continuous(breaks = scales::pretty_breaks(10)) + scale_color_brewer(palette="Paired")

#---------------------------------------------------
# Graph without Russia, total trade

ggplot(data, aes(x=year)) + 
  geom_line(aes(y = logcu, color = "CU/EUEA"), size = 2, linetype = 4) +
  geom_line(aes(y = Belarus, color = "Belarus"), size = 1) +
  geom_line(aes(y = Kazakhstan, color = "Kazakhstan"), size = 1) +
  geom_line(aes(y = Kyrgyzstan, color = "Kyrgyzstan"), size = 1) +
  geom_line(aes(y = Armenia, color = "Armenia"), size = 1) +
  labs (x ="year", y = "Log of total trade", color = 'Entity') +
  geom_vline(xintercept = 2010 , color = "red", linetype = "dashed") + 
  geom_vline(xintercept = 2012 , color = "red", linetype = "dashed") + 
  geom_vline(xintercept = 2015 , color = "red", linetype = "dashed") +
  geom_vline(xintercept = 2016 , color = "red", linetype = "dashed") +
  theme_classic() + 
  scale_x_continuous(breaks = scales::pretty_breaks(10)) + scale_color_brewer(palette="Paired")

#---------------------------------------------------
# This code uses exp.csv and produces Figure A3: exports

ggplot(exp, aes(x=year)) + 
  geom_line(aes(y = cu, color = "CU/EUEA"), size = 2, linetype = 4) +
  geom_line(aes(y = Belarus, color = "Belarus"), size = 1) +
  geom_line(aes(y = Kazakhstan, color = "Kazakhstan"), size = 1) +
  geom_line(aes(y = Russia, color = "Russia"), size = 1) +
  geom_line(aes(y = Kyrgyzstan, color = "Kyrgyzstan"), size = 1) +
  geom_line(aes(y = Armenia, color = "Armenia"), size = 1) +
  labs (x ="year", y = "Log of total exports", color = 'Entity') +
  geom_vline(xintercept = 2010 , color = "red", linetype = "dashed") + 
  geom_vline(xintercept = 2012 , color = "red", linetype = "dashed") + 
  geom_vline(xintercept = 2015 , color = "red", linetype = "dashed") +
  geom_vline(xintercept = 2016 , color = "red", linetype = "dashed") +
  theme_classic() + 
  scale_x_continuous(breaks = scales::pretty_breaks(10)) + scale_color_brewer(palette="Paired")

#---------------------------------------------------
# This code uses imp.csv and produces Figure A4: imports

ggplot(imp, aes(x=year)) + 
  geom_line(aes(y = cu, color = "CU/EUEA"), size = 2, linetype = 4) +
  geom_line(aes(y = Belarus, color = "Belarus"), size = 1) +
  geom_line(aes(y = Kazakhstan, color = "Kazakhstan"), size = 1) +
  geom_line(aes(y = Russia, color = "Russia"), size = 1) +
  geom_line(aes(y = Kyrgyzstan, color = "Kyrgyzstan"), size = 1) +
  geom_line(aes(y = Armenia, color = "Armenia"), size = 1) +
  labs (x ="year", y = "Log of total imports", color = 'Entity') +
  geom_vline(xintercept = 2010 , color = "red", linetype = "dashed") + 
  geom_vline(xintercept = 2012 , color = "red", linetype = "dashed") + 
  geom_vline(xintercept = 2015 , color = "red", linetype = "dashed") +
  geom_vline(xintercept = 2016 , color = "red", linetype = "dashed") +
  theme_classic() + 
  scale_x_continuous(breaks = scales::pretty_breaks(10)) + scale_color_brewer(palette="Paired")
