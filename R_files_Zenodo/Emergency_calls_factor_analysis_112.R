library(MuMIn)
library(AER)
library(glmmTMB)
library(Hmisc)
library(pscl)
library(dplyr)
library(AICcmodavg)
library(lme4)
library(nlme)
library(MASS)

setwd(' ') # Mihai's Dropbox


data.112 <- read.csv('dataset112_30Nov2021.csv')
str(data.112)

data.112$forMN <- as.numeric(data.112$forMN)
data.112$urbMN <- as.numeric(data.112$urbMN)
data.112$agMN <- as.numeric(data.112$agMN)
summary(data.112)

# remove datapoints with NA values (from land cover extractions at 10 x 10 km)
data.112 <- na.omit(data.112)
summary(data.112)

c <- rcorr(as.matrix(data.112[,5:22]), type="pearson") 
str(c)
write.csv(c$r, "corelatii_toate.csv")

# subset data for bears based on bear distribution
urs.112 <- subset(data.112, range_urs == 1 & relief != "C")
str(urs.112)
summary(urs.112)
hist(urs.112$urs)

sum(urs.112$urs < 1)
sum(urs.112$urs >= 1)

urs.112.withdata <- subset(urs.112, urs >= 1)
hist(urs.112.withdata$urs, breaks = 100)

### investigate correlations between variables
urs.num <- as.data.frame(urs.112[,5:22])
str(urs.num)
summary(urs.num)
curs <- rcorr(as.matrix(urs.num), type="pearson") 
str(curs)
write.csv(curs$r, "corelatii_urs.csv")


#### MODELS (HURDLE) FOR URS - NULL MODEL + ALL MODEL + 30MODELS 

#scale(Populatie) + scale(pop2010to2019) + scale(Suprafata)+
#  + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
#  + scale(perc_spverzi) + scale(perc_parcurb) + 
#  + scale(Proc_pad) + scale(Proc_agri) + scale(forMN) + scale(agMN) +
#  + scale(urbMN) + scale(Conectivitate_natur.) + relief -1

mod.hurs.0 <-hurdle(urs ~ 1, data = urs.112, 
                    dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.0)

mod.hurs.100 <- hurdle(urs ~ scale(Populatie) + scale(pop2010to2019) +
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(forMN) +
                         + scale(urbMN) + scale(Conectivitate_natur.) + scale(StuSup) +
                         + scale(PopAgri) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.100)

mod.hurs.101 <- hurdle(urs ~ scale(Populatie) + scale(pop2010to2019) +
                         + scale(StuSup) +
                         + scale(PopAgri) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.101)

mod.hurs.102 <- hurdle(urs ~ scale(Suprafata) + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(forMN) +
                         + scale(urbMN) + scale(Conectivitate_natur.) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.102)

mod.hurs.1 <- hurdle(urs ~ scale(Populatie) + 
                       + scale(D_strazi_urb_2018) + 
                       + scale(perc_spverzi) + scale(StuSup)
                       -1, data = urs.112, 
                     dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.1)

mod.hurs.2 <- hurdle(urs ~ scale(Suprafata)+
                         + scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(StuSup) + scale(PopAgri) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.2)

mod.hurs.3 <- hurdle(urs ~ scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(PopAgri) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.3)

mod.hurs.4 <- hurdle(urs ~ scale(Populatie) +
                         + scale(rata_crestere_supr_construita) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(PopAgri) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.4)

mod.hurs.5 <- hurdle(urs ~ scale(pop2010to2019) +
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                         + scale(StuSup) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.5)

mod.hurs.6 <- hurdle(urs ~ scale(Populatie) + scale(pop2010to2019) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(StuSup) + scale(PopAgri)  -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.6)

mod.hurs.7 <- hurdle(urs ~ scale(pop2010to2019) +
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.7)

mod.hurs.8 <- hurdle(urs ~ scale(pop2010to2019) + scale(Suprafata)+
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                         + scale(StuSup) + scale(PopAgri) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.8)

mod.hurs.9 <- hurdle(urs ~ scale(Suprafata)+
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.9)

mod.hurs.10 <- hurdle(urs ~ scale(Populatie) 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(StuSup) + scale(PopAgri) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.10)

mod.hurs.11 <- hurdle(urs ~ + scale(forMN) +
                         + scale(Conectivitate_natur.) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.11)

mod.hurs.12 <- hurdle(urs ~ + scale(agMN) +
                         + scale(Conectivitate_natur.) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.12)

mod.hurs.13 <- hurdle(urs ~ scale(agMN) +
                         + scale(urbMN) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.13)

mod.hurs.14 <- hurdle(urs ~ + scale(forMN) +
                         + scale(urbMN) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.14)

mod.hurs.15 <- hurdle(urs ~ scale(forMN) +
                         + scale(urbMN) + scale(Conectivitate_natur.) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.15)

mod.hurs.16 <- hurdle(urs ~ scale(Proc_pad) +
                         + scale(Conectivitate_natur.) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.16)

mod.hurs.17 <- hurdle(urs ~ scale(Proc_agri) +
                         + scale(Conectivitate_natur.) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.17)

mod.hurs.18 <- hurdle(urs ~ scale(Proc_agri) +
                         + scale(urbMN) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.18)

mod.hurs.19 <- hurdle(urs ~ scale(Proc_pad) +
                         + scale(urbMN) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.19)

mod.hurs.20 <- hurdle(urs ~ scale(Proc_pad) +
                         + scale(urbMN) + scale(Conectivitate_natur.) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.20)

mod.hurs.21 <- hurdle(urs ~ scale(Populatie) + 
                         + scale(D_strazi_urb_2018) + 
                         + scale(urbMN) + scale(Conectivitate_natur.) 
                         + scale(StuSup) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.21)

mod.hurs.22 <- hurdle(urs ~ scale(Suprafata)+
                         + scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(forMN) +
                         + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.22)

mod.hurs.23 <- hurdle(urs ~ scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(agMN) +
                         + scale(urbMN) + scale(PopAgri) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.23)

mod.hurs.24 <- hurdle(urs ~ scale(Populatie) + 
                         + scale(rata_crestere_supr_construita) +
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(Proc_pad) +
                         + scale(urbMN) + scale(PopAgri) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.24)

mod.hurs.25 <- hurdle(urs ~ scale(pop2010to2019) +
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                         + scale(Proc_agri) +
                         + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.25)

mod.hurs.26 <- hurdle(urs ~ scale(Populatie) + scale(pop2010to2019) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(Proc_agri) +
                         + scale(urbMN) + scale(StuSup) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.26)

mod.hurs.27 <- hurdle(urs ~ scale(pop2010to2019) +
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(agMN) +
                         + scale(Conectivitate_natur.) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.27)

mod.hurs.28 <- hurdle(urs ~ scale(pop2010to2019) + scale(Suprafata)+
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                         + scale(forMN) +
                         + scale(urbMN) + scale(StuSup) + scale(PopAgri) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.28)

mod.hurs.29 <- hurdle(urs ~ scale(Suprafata)+
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(Proc_pad) +
                         + scale(Conectivitate_natur.) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.29)

mod.hurs.30 <- hurdle(urs ~ scale(Populatie) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(Proc_pad) +
                         + scale(urbMN) + scale(Conectivitate_natur.) 
                         + scale(StuSup) + scale(PopAgri) -1, data = urs.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hurs.30)

#MODEL SELECTION
model.sel(mod.hurs.0, mod.hurs.100, mod.hurs.101, mod.hurs.102, 
          mod.hurs.1, mod.hurs.2, mod.hurs.3, mod.hurs.4, mod.hurs.5, 
          mod.hurs.6, mod.hurs.7, mod.hurs.8, mod.hurs.9, mod.hurs.10, 
          mod.hurs.11, mod.hurs.12, mod.hurs.13, mod.hurs.14, mod.hurs.15, 
          mod.hurs.16, mod.hurs.17, mod.hurs.18, mod.hurs.19, mod.hurs.20, 
          mod.hurs.21, mod.hurs.22, mod.hurs.23, mod.hurs.24, mod.hurs.25, 
          mod.hurs.26, mod.hurs.27, mod.hurs.28, mod.hurs.29, mod.hurs.30)
mshurs <- model.sel(mod.hurs.0, mod.hurs.100, mod.hurs.101, mod.hurs.102,
                    mod.hurs.1, mod.hurs.2, mod.hurs.3, mod.hurs.4, mod.hurs.5, 
                    mod.hurs.6, mod.hurs.7, mod.hurs.8, mod.hurs.9, mod.hurs.10, 
                    mod.hurs.11, mod.hurs.12, mod.hurs.13, mod.hurs.14, mod.hurs.15, 
                    mod.hurs.16, mod.hurs.17, mod.hurs.18, mod.hurs.19, mod.hurs.20, 
                    mod.hurs.21, mod.hurs.22, mod.hurs.23, mod.hurs.24, mod.hurs.25, 
                    mod.hurs.26, mod.hurs.27, mod.hurs.28, mod.hurs.29, mod.hurs.30)

write.csv(mshurs, "model_selection_all_URS.csv")

# model averaging
avgmod <- model.avg(mshurs,cumsum(weight)<= 0.95)
summary(avgmod)



################################################################
#### MODELS (HURDLE) FOR CERB - NULL MODEL + ALL MODEL + 30MODELS  

# subset data for cerb based on !!!bear distribution
cerb.112 <- subset(data.112, range_urs == 1 & relief != "C")
str(cerb.112)
summary(cerb.112)
hist(cerb.112$cerb)

sum(cerb.112$cerb < 1)
sum(cerb.112$cerb >= 1)

cerb.112.withdata <- subset(cerb.112, cerb >= 1)
hist(cerb.112.withdata$cerb, breaks = 100)



mod.hcerb.0 <-hurdle(cerb ~ 1, data = cerb.112, 
                    dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.0)

mod.hcerb.100 <- hurdle(cerb ~ scale(Populatie) + scale(pop2010to2019) +
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(forMN) +
                         + scale(urbMN) + scale(Conectivitate_natur.) + scale(StuSup) +
                         + scale(PopAgri) -1, data = cerb.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.100)

mod.hcerb.101 <- hurdle(cerb ~ scale(Populatie) + scale(pop2010to2019) +
                         + scale(StuSup) +
                         + scale(PopAgri) -1, data = cerb.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.101)

mod.hcerb.102 <- hurdle(cerb ~ scale(Suprafata) + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(forMN) +
                         + scale(urbMN) + scale(Conectivitate_natur.) -1, data = cerb.112, 
                       dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.102)

mod.hcerb.1 <- hurdle(cerb ~ scale(Populatie) + 
                       + scale(D_strazi_urb_2018) + 
                       + scale(perc_spverzi) + scale(StuSup)
                     -1, data = cerb.112, 
                     dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.1)

mod.hcerb.2 <- hurdle(cerb ~ scale(Suprafata)+
                       + scale(D_strazi_urb_2018) + 
                       + scale(perc_spverzi) + scale(perc_parcurb) + 
                       + scale(StuSup) + scale(PopAgri) -1, data = cerb.112, 
                     dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.2)

mod.hcerb.3 <- hurdle(cerb ~ scale(D_strazi_urb_2018) + 
                       + scale(perc_spverzi) + scale(perc_parcurb) + 
                       + scale(PopAgri) -1, data = cerb.112, 
                     dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.3)

mod.hcerb.4 <- hurdle(cerb ~ scale(Populatie) +
                       + scale(rata_crestere_supr_construita) + 
                       + scale(perc_spverzi) + scale(perc_parcurb) + 
                       + scale(PopAgri) -1, data = cerb.112, 
                     dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.4)

mod.hcerb.5 <- hurdle(cerb ~ scale(pop2010to2019) +
                       + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                       + scale(StuSup) -1, data = cerb.112, 
                     dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.5)

mod.hcerb.6 <- hurdle(cerb ~ scale(Populatie) + scale(pop2010to2019) + 
                       + scale(perc_spverzi) + scale(perc_parcurb) + 
                       + scale(StuSup) + scale(PopAgri)  -1, data = cerb.112, 
                     dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.6)

mod.hcerb.7 <- hurdle(cerb ~ scale(pop2010to2019) +
                       + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                       + scale(perc_spverzi) + scale(perc_parcurb) + 
                       -1, data = cerb.112, 
                     dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.7)

mod.hcerb.8 <- hurdle(cerb ~ scale(pop2010to2019) + scale(Suprafata)+
                       + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                       + scale(StuSup) + scale(PopAgri) -1, data = cerb.112, 
                     dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.8)

mod.hcerb.9 <- hurdle(cerb ~ scale(Suprafata)+
                       + scale(perc_spverzi) + scale(perc_parcurb) + 
                       -1, data = cerb.112, 
                     dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.9)

mod.hcerb.10 <- hurdle(cerb ~ scale(Populatie) 
                      + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(StuSup) + scale(PopAgri) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.10)

mod.hcerb.11 <- hurdle(cerb ~ + scale(forMN) +
                        + scale(Conectivitate_natur.) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.11)

mod.hcerb.12 <- hurdle(cerb ~ + scale(agMN) +
                        + scale(Conectivitate_natur.) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.12)

mod.hcerb.13 <- hurdle(cerb ~ scale(agMN) +
                        + scale(urbMN) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.13)

mod.hcerb.14 <- hurdle(cerb ~ + scale(forMN) +
                        + scale(urbMN) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.14)

mod.hcerb.15 <- hurdle(cerb ~ scale(forMN) +
                        + scale(urbMN) + scale(Conectivitate_natur.) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.15)

mod.hcerb.16 <- hurdle(cerb ~ scale(Proc_pad) +
                        + scale(Conectivitate_natur.) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.16)

mod.hcerb.17 <- hurdle(cerb ~ scale(Proc_agri) +
                        + scale(Conectivitate_natur.) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.17)

mod.hcerb.18 <- hurdle(cerb ~ scale(Proc_agri) +
                        + scale(urbMN) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.18)

mod.hcerb.19 <- hurdle(cerb ~ scale(Proc_pad) +
                        + scale(urbMN) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.19)

mod.hcerb.20 <- hurdle(cerb ~ scale(Proc_pad) +
                        + scale(urbMN) + scale(Conectivitate_natur.) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.20)

mod.hcerb.21 <- hurdle(cerb ~ scale(Populatie) + 
                        + scale(D_strazi_urb_2018) + 
                        + scale(urbMN) + scale(Conectivitate_natur.) 
                      + scale(StuSup) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.21)

mod.hcerb.22 <- hurdle(cerb ~ scale(Suprafata)+
                        + scale(D_strazi_urb_2018) + 
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(forMN) +
                        + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.22)

mod.hcerb.23 <- hurdle(cerb ~ scale(D_strazi_urb_2018) + 
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(agMN) +
                        + scale(urbMN) + scale(PopAgri) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.23)

mod.hcerb.24 <- hurdle(cerb ~ scale(Populatie) + 
                        + scale(rata_crestere_supr_construita) +
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(Proc_pad) +
                        + scale(urbMN) + scale(PopAgri) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.24)

mod.hcerb.25 <- hurdle(cerb ~ scale(pop2010to2019) +
                        + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                        + scale(Proc_agri) +
                        + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.25)

mod.hcerb.26 <- hurdle(cerb ~ scale(Populatie) + scale(pop2010to2019) + 
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(Proc_agri) +
                        + scale(urbMN) + scale(StuSup) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.26)

mod.hcerb.27 <- hurdle(cerb ~ scale(pop2010to2019) +
                        + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(agMN) +
                        + scale(Conectivitate_natur.) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.27)

mod.hcerb.28 <- hurdle(cerb ~ scale(pop2010to2019) + scale(Suprafata)+
                        + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                        + scale(forMN) +
                        + scale(urbMN) + scale(StuSup) + scale(PopAgri) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.28)

mod.hcerb.29 <- hurdle(cerb ~ scale(Suprafata)+
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(Proc_pad) +
                        + scale(Conectivitate_natur.) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.29)

mod.hcerb.30 <- hurdle(cerb ~ scale(Populatie) + 
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(Proc_pad) +
                        + scale(urbMN) + scale(Conectivitate_natur.) 
                      + scale(StuSup) + scale(PopAgri) -1, data = cerb.112, 
                      dist = "negbin", zero.dist = "binomial")
summary(mod.hcerb.30)


#MODEL SELECTION
model.sel(mod.hcerb.0, mod.hcerb.100, mod.hcerb.101, mod.hcerb.102,
          mod.hcerb.1, mod.hcerb.2, mod.hcerb.3, mod.hcerb.4, mod.hcerb.5, 
          mod.hcerb.6, mod.hcerb.7, mod.hcerb.8, mod.hcerb.9, mod.hcerb.10, 
          mod.hcerb.11, mod.hcerb.12, mod.hcerb.13, mod.hcerb.14, mod.hcerb.15, 
          mod.hcerb.16, mod.hcerb.17, mod.hcerb.18, mod.hcerb.19, mod.hcerb.20, 
          mod.hcerb.21, mod.hcerb.22, mod.hcerb.23, mod.hcerb.24, mod.hcerb.25, 
          mod.hcerb.26, mod.hcerb.27, mod.hcerb.28, mod.hcerb.29, mod.hcerb.30)
mshcerb <- model.sel(mod.hcerb.0, mod.hcerb.100, mod.hcerb.101, mod.hcerb.102,
                    mod.hcerb.1, mod.hcerb.2, mod.hcerb.3, mod.hcerb.4, mod.hcerb.5, 
                    mod.hcerb.6, mod.hcerb.7, mod.hcerb.8, mod.hcerb.9, mod.hcerb.10, 
                    mod.hcerb.11, mod.hcerb.12, mod.hcerb.13, mod.hcerb.14, mod.hcerb.15, 
                    mod.hcerb.16, mod.hcerb.17, mod.hcerb.18, mod.hcerb.19, mod.hcerb.20, 
                    mod.hcerb.21, mod.hcerb.22, mod.hcerb.23, mod.hcerb.24, mod.hcerb.25, 
                    mod.hcerb.26, mod.hcerb.27, mod.hcerb.28, mod.hcerb.29, mod.hcerb.30)
#str(ms)
write.csv(mshcerb, "model_selection_all_CERB.csv")

# model averaging
# model averaging
avgmod <- model.avg(mshcerb,cumsum(weight)<= 0.95)
summary(avgmod)


#################################################################
#### MODELS (GLM) FOR CAPRIOR - NULL MODEL + ALL MODEL + 30MODELS           

mod.glmcap.0 <-glm.nb(caprioara ~ 1, data = data.112)
summary(mod.glmcap.0)

mod.glmcap.100 <- glm.nb(caprioara ~ scale(Populatie) + scale(pop2010to2019) +
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + 
                         + scale(forMN) +
                         + scale(urbMN) + scale(Conectivitate_natur.) + relief, data = data.112)
summary(mod.glmcap.100)

mod.glmcap.101 <- glm.nb(caprioara ~ scale(Populatie) + scale(pop2010to2019) +
                          + scale(StuSup) +
                          + scale(PopAgri), data = data.112)
summary(mod.glmcap.101)

mod.glmcap.102 <- glm.nb(caprioara ~ scale(Suprafata) + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(forMN) +
                          + scale(urbMN) + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmcap.102)

mod.glmcap.1 <- glm.nb(caprioara ~ scale(Populatie) + 
                       + scale(D_strazi_urb_2018) + 
                       +scale(perc_spverzi) + scale(StuSup) + relief,
                       data = data.112)
summary(mod.glmcap.1)

mod.glmcap.2 <- glm.nb(caprioara ~ scale(Suprafata)+
                       + scale(D_strazi_urb_2018) + 
                       + scale(perc_spverzi) + scale(perc_parcurb) + scale(StuSup) + scale(PopAgri), data = data.112)
summary(mod.glmcap.2)

mod.glmcap.3 <- glm.nb(caprioara ~ scale(D_strazi_urb_2018) + 
                       + scale(perc_spverzi) + scale(perc_parcurb) + scale(PopAgri) +
                       + relief, data = data.112)
summary(mod.glmcap.3)

mod.glmcap.4 <- glm.nb(caprioara ~ scale(Populatie) +
                       + scale(rata_crestere_supr_construita) + 
                       + scale(perc_spverzi) + scale(perc_parcurb) + scale(PopAgri) +
                       + relief, data = data.112)
summary(mod.glmcap.4)

mod.glmcap.5 <- glm.nb(caprioara ~ scale(pop2010to2019) +
                       + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + scale(StuSup) +  
                       + relief, data = data.112)
summary(mod.glmcap.5)

mod.glmcap.6 <- glm.nb(caprioara ~ scale(Populatie) + scale(pop2010to2019) + 
                       + scale(perc_spverzi) + scale(perc_parcurb) + scale(StuSup) + scale(PopAgri)   
                       , data = data.112)
summary(mod.glmcap.6)

mod.glmcap.7 <- glm.nb(caprioara ~ scale(pop2010to2019) +
                       + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018)  
                       + scale(perc_spverzi) + scale(perc_parcurb) 
                       , data = data.112)
summary(mod.glmcap.7)

mod.glmcap.8 <- glm.nb(caprioara ~ scale(pop2010to2019) + scale(Suprafata)+
                       + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + scale(StuSup) + scale(PopAgri) + 
                       + relief, data = data.112)
summary(mod.glmcap.8)

mod.glmcap.9 <- glm.nb(caprioara ~ scale(Suprafata)+
                       + scale(perc_spverzi) + scale(perc_parcurb) 
                       , data = data.112)
summary(mod.glmcap.9)

mod.glmcap.10 <- glm.nb(caprioara ~ scale(Populatie) 
                      + scale(perc_spverzi) + scale(perc_parcurb) + scale(StuSup) + scale(PopAgri) +
                      + relief, data = data.112)
summary(mod.glmcap.10)

mod.glmcap.11 <- glm.nb(caprioara ~ + scale(forMN) +
                        + scale(Conectivitate_natur.) + relief, data = data.112)
summary(mod.glmcap.11)

mod.glmcap.12 <- glm.nb(caprioara ~ + scale(agMN) +
                        + scale(Conectivitate_natur.) + relief, data = data.112)
summary(mod.glmcap.12)

mod.glmcap.13 <- glm.nb(caprioara ~ scale(agMN) +
                        + scale(urbMN) +relief, data = data.112)
summary(mod.glmcap.13)

mod.glmcap.14 <- glm.nb(caprioara ~ + scale(forMN) +
                        + scale(urbMN) +relief, data = data.112)
summary(mod.glmcap.14)

mod.glmcap.15 <- glm.nb(caprioara ~ scale(forMN) +
                        + scale(urbMN) + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmcap.15)

mod.glmcap.16 <- glm.nb(caprioara ~ scale(Proc_pad) +
                        + scale(Conectivitate_natur.) +relief, data = data.112)
summary(mod.glmcap.16)

mod.glmcap.17 <- glm.nb(caprioara ~ scale(Proc_agri) +
                        + scale(Conectivitate_natur.) +relief, data = data.112)
summary(mod.glmcap.17)

mod.glmcap.18 <- glm.nb(caprioara ~ scale(Proc_agri) +
                        + scale(urbMN) +relief, data = data.112)
summary(mod.glmcap.18)

mod.glmcap.19 <- glm.nb(caprioara ~ scale(Proc_pad) +
                        + scale(urbMN) +relief, data = data.112)
summary(mod.glmcap.19)

mod.glmcap.20 <- glm.nb(caprioara ~ scale(Proc_pad) +
                        + scale(urbMN) + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmcap.20)

mod.glmcap.21 <- glm.nb(caprioara ~ scale(Populatie) + 
                        + scale(D_strazi_urb_2018) + 
                        + scale(urbMN) + scale(Conectivitate_natur.) + scale(StuSup) + relief,
                        data = data.112)
summary(mod.glmcap.21)

mod.glmcap.22 <- glm.nb(caprioara ~ scale(Suprafata)+
                        + scale(D_strazi_urb_2018) + 
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(forMN) +
                        + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri), data = data.112)
summary(mod.glmcap.22)

mod.glmcap.23 <- glm.nb(caprioara ~ scale(D_strazi_urb_2018) + 
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(agMN) +
                        + scale(urbMN) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmcap.23)

mod.glmcap.24 <- glm.nb(caprioara ~ scale(Populatie) + 
                        + scale(rata_crestere_supr_construita) +
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(Proc_pad) +
                        + scale(urbMN) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmcap.24)

mod.glmcap.25 <- glm.nb(caprioara ~ scale(pop2010to2019) +
                        + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                        + scale(Proc_agri) +
                        + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmcap.25)

mod.glmcap.26 <- glm.nb(caprioara ~ scale(Populatie) + scale(pop2010to2019) + 
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(Proc_agri) +
                        + scale(urbMN) + scale(StuSup), data = data.112)
summary(mod.glmcap.26)

mod.glmcap.27 <- glm.nb(caprioara ~ scale(pop2010to2019) +
                        + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(agMN) +
                        + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmcap.27)

mod.glmcap.28 <- glm.nb(caprioara ~ scale(pop2010to2019) + scale(Suprafata)+
                        + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                        + scale(forMN) +
                        + scale(urbMN) + scale(StuSup) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmcap.28)

mod.glmcap.29 <- glm.nb(caprioara ~ scale(Suprafata)+
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(Proc_pad) +
                        + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmcap.29)

mod.glmcap.30 <- glm.nb(caprioara ~ scale(Populatie) + 
                        + scale(perc_spverzi) + scale(perc_parcurb) + 
                        + scale(Proc_pad) +
                        + scale(urbMN) + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri) + relief, 
                        data = data.112)
summary(mod.glmcap.30)

#MODEL SELECTION
model.sel(mod.glmcap.0, mod.glmcap.100, mod.glmcap.101, mod.glmcap.102,
          mod.glmcap.1, mod.glmcap.2, mod.glmcap.3, mod.glmcap.4, mod.glmcap.5, 
          mod.glmcap.6, mod.glmcap.7, mod.glmcap.8, mod.glmcap.9, mod.glmcap.10, 
          mod.glmcap.11, mod.glmcap.12, mod.glmcap.13, mod.glmcap.14, mod.glmcap.15, 
          mod.glmcap.16, mod.glmcap.17, mod.glmcap.18, mod.glmcap.19, mod.glmcap.20, 
          mod.glmcap.21, mod.glmcap.22, mod.glmcap.23, mod.glmcap.24, mod.glmcap.25, 
          mod.glmcap.26, mod.glmcap.27, mod.glmcap.28, mod.glmcap.29, mod.glmcap.30)
msglmcap <- model.sel(mod.glmcap.0, mod.glmcap.100, mod.glmcap.101, mod.glmcap.102,
                    mod.glmcap.1, mod.glmcap.2, mod.glmcap.3, mod.glmcap.4, mod.glmcap.5, 
                    mod.glmcap.6, mod.glmcap.7, mod.glmcap.8, mod.glmcap.9, mod.glmcap.10, 
                    mod.glmcap.11, mod.glmcap.12, mod.glmcap.13, mod.glmcap.14, mod.glmcap.15, 
                    mod.glmcap.16, mod.glmcap.17, mod.glmcap.18, mod.glmcap.19, mod.glmcap.20, 
                    mod.glmcap.21, mod.glmcap.22, mod.glmcap.23, mod.glmcap.24, mod.glmcap.25, 
                    mod.glmcap.26, mod.glmcap.27, mod.glmcap.28, mod.glmcap.29, mod.glmcap.30)
#str(ms)
write.csv(msglmcap, "model_selection_all_CAP.csv")

# model averaging
avgmod <- model.avg(msglmcap,cumsum(weight)<= 0.95)
summary(avgmod)




#################################################################
#### MODELS (GLM) FOR MISTRET - NULL MODEL + ALL MODEL + 30MODELS           

mod.glmmis.0 <-glm.nb(mistret ~ 1, data = data.112)
summary(mod.glmmis.0)

mod.glmmis.100 <- glm.nb(mistret ~ scale(Populatie) + scale(pop2010to2019) +
                           + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                           + scale(perc_spverzi) + scale(perc_parcurb) + 
                           + scale(forMN) +
                           + scale(urbMN) + scale(Conectivitate_natur.) + relief, data = data.112)
summary(mod.glmmis.100)

mod.glmmis.101 <- glm.nb(mistret ~ scale(Populatie) + scale(pop2010to2019) +
                           + scale(StuSup) +
                           + scale(PopAgri), data = data.112)
summary(mod.glmmis.101)

mod.glmmis.102 <- glm.nb(mistret ~ scale(Suprafata) + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                           + scale(perc_spverzi) + scale(perc_parcurb) + 
                           + scale(forMN) +
                           + scale(urbMN) + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmmis.102)

mod.glmmis.1 <- glm.nb(mistret ~ scale(Populatie) + 
                         + scale(D_strazi_urb_2018) + 
                         +scale(perc_spverzi) + scale(StuSup) + relief,
                       data = data.112)
summary(mod.glmmis.1)

mod.glmmis.2 <- glm.nb(mistret ~ scale(Suprafata)+
                         + scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + scale(StuSup) + scale(PopAgri), data = data.112)
summary(mod.glmmis.2)

mod.glmmis.3 <- glm.nb(mistret ~ scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + scale(PopAgri) +
                         + relief, data = data.112)
summary(mod.glmmis.3)

mod.glmmis.4 <- glm.nb(mistret ~ scale(Populatie) +
                         + scale(rata_crestere_supr_construita) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + scale(PopAgri) +
                         + relief, data = data.112)
summary(mod.glmmis.4)

mod.glmmis.5 <- glm.nb(mistret ~ scale(pop2010to2019) +
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + scale(StuSup) +  
                         + relief, data = data.112)
summary(mod.glmmis.5)

mod.glmmis.6 <- glm.nb(mistret ~ scale(Populatie) + scale(pop2010to2019) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + scale(StuSup) + scale(PopAgri)   
                       , data = data.112)
summary(mod.glmmis.6)

mod.glmmis.7 <- glm.nb(mistret ~ scale(pop2010to2019) +
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018)  
                       + scale(perc_spverzi) + scale(perc_parcurb) 
                       , data = data.112)
summary(mod.glmmis.7)

mod.glmmis.8 <- glm.nb(mistret ~ scale(pop2010to2019) + scale(Suprafata)+
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + scale(StuSup) + scale(PopAgri) + 
                         + relief, data = data.112)
summary(mod.glmmis.8)

mod.glmmis.9 <- glm.nb(mistret ~ scale(Suprafata)+
                         + scale(perc_spverzi) + scale(perc_parcurb) 
                       , data = data.112)
summary(mod.glmmis.9)

mod.glmmis.10 <- glm.nb(mistret ~ scale(Populatie) 
                        + scale(perc_spverzi) + scale(perc_parcurb) + scale(StuSup) + scale(PopAgri) +
                          + relief, data = data.112)
summary(mod.glmmis.10)

mod.glmmis.11 <- glm.nb(mistret ~ + scale(forMN) +
                          + scale(Conectivitate_natur.) + relief, data = data.112)
summary(mod.glmmis.11)

mod.glmmis.12 <- glm.nb(mistret ~ + scale(agMN) +
                          + scale(Conectivitate_natur.) + relief, data = data.112)
summary(mod.glmmis.12)

mod.glmmis.13 <- glm.nb(mistret ~ scale(agMN) +
                          + scale(urbMN) +relief, data = data.112)
summary(mod.glmmis.13)

mod.glmmis.14 <- glm.nb(mistret ~ + scale(forMN) +
                          + scale(urbMN) +relief, data = data.112)
summary(mod.glmmis.14)

mod.glmmis.15 <- glm.nb(mistret ~ scale(forMN) +
                          + scale(urbMN) + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmmis.15)

mod.glmmis.16 <- glm.nb(mistret ~ scale(Proc_pad) +
                          + scale(Conectivitate_natur.) +relief, data = data.112)
summary(mod.glmmis.16)

mod.glmmis.17 <- glm.nb(mistret ~ scale(Proc_agri) +
                          + scale(Conectivitate_natur.) +relief, data = data.112)
summary(mod.glmmis.17)

mod.glmmis.18 <- glm.nb(mistret ~ scale(Proc_agri) +
                          + scale(urbMN) +relief, data = data.112)
summary(mod.glmmis.18)

mod.glmmis.19 <- glm.nb(mistret ~ scale(Proc_pad) +
                          + scale(urbMN) +relief, data = data.112)
summary(mod.glmmis.19)

mod.glmmis.20 <- glm.nb(mistret ~ scale(Proc_pad) +
                          + scale(urbMN) + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmmis.20)

mod.glmmis.21 <- glm.nb(mistret ~ scale(Populatie) + 
                          + scale(D_strazi_urb_2018) + 
                          + scale(urbMN) + scale(Conectivitate_natur.) + scale(StuSup) + relief,
                        data = data.112)
summary(mod.glmmis.21)

mod.glmmis.22 <- glm.nb(mistret ~ scale(Suprafata)+
                          + scale(D_strazi_urb_2018) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(forMN) +
                          + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri), data = data.112)
summary(mod.glmmis.22)

mod.glmmis.23 <- glm.nb(mistret ~ scale(D_strazi_urb_2018) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(agMN) +
                          + scale(urbMN) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmmis.23)

mod.glmmis.24 <- glm.nb(mistret ~ scale(Populatie) + 
                          + scale(rata_crestere_supr_construita) +
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(Proc_pad) +
                          + scale(urbMN) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmmis.24)

mod.glmmis.25 <- glm.nb(mistret ~ scale(pop2010to2019) +
                          + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                          + scale(Proc_agri) +
                          + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmmis.25)

mod.glmmis.26 <- glm.nb(mistret ~ scale(Populatie) + scale(pop2010to2019) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(Proc_agri) +
                          + scale(urbMN) + scale(StuSup), data = data.112)
summary(mod.glmmis.26)

mod.glmmis.27 <- glm.nb(mistret ~ scale(pop2010to2019) +
                          + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(agMN) +
                          + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmmis.27)

mod.glmmis.28 <- glm.nb(mistret ~ scale(pop2010to2019) + scale(Suprafata)+
                          + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                          + scale(forMN) +
                          + scale(urbMN) + scale(StuSup) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmmis.28)

mod.glmmis.29 <- glm.nb(mistret ~ scale(Suprafata)+
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(Proc_pad) +
                          + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmmis.29)

mod.glmmis.30 <- glm.nb(mistret ~ scale(Populatie) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(Proc_pad) +
                          + scale(urbMN) + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri) + relief, 
                        data = data.112)
summary(mod.glmmis.30)

#MODEL SELECTION
model.sel(mod.glmmis.0, mod.glmmis.100, mod.glmmis.101, mod.glmmis.102, 
          mod.glmmis.1, mod.glmmis.2, mod.glmmis.3, mod.glmmis.4, mod.glmmis.5, 
          mod.glmmis.6, mod.glmmis.7, mod.glmmis.8, mod.glmmis.9, mod.glmmis.10, 
          mod.glmmis.11, mod.glmmis.12, mod.glmmis.13, mod.glmmis.14, mod.glmmis.15, 
          mod.glmmis.16, mod.glmmis.17, mod.glmmis.18, mod.glmmis.19, mod.glmmis.20, 
          mod.glmmis.21, mod.glmmis.22, mod.glmmis.23, mod.glmmis.24, mod.glmmis.25, 
          mod.glmmis.26, mod.glmmis.27, mod.glmmis.28, mod.glmmis.29, mod.glmmis.30)
msglmmis <- model.sel(mod.glmmis.0, mod.glmmis.100, mod.glmmis.101, mod.glmmis.102,
                      mod.glmmis.1, mod.glmmis.2, mod.glmmis.3, mod.glmmis.4, mod.glmmis.5, 
                      mod.glmmis.6, mod.glmmis.7, mod.glmmis.8, mod.glmmis.9, mod.glmmis.10, 
                      mod.glmmis.11, mod.glmmis.12, mod.glmmis.13, mod.glmmis.14, mod.glmmis.15, 
                      mod.glmmis.16, mod.glmmis.17, mod.glmmis.18, mod.glmmis.19, mod.glmmis.20, 
                      mod.glmmis.21, mod.glmmis.22, mod.glmmis.23, mod.glmmis.24, mod.glmmis.25, 
                      mod.glmmis.26, mod.glmmis.27, mod.glmmis.28, mod.glmmis.29, mod.glmmis.30)
#str(ms)
write.csv(msglmmis, "model_selection_all_MIS.csv")

# model averaging
#avgmod <- model.avg(msglmmis)
#summary(avgmod)




#################################################################
#### MODELS (GLM) FOR VULPE - NULL MODEL + ALL MODEL + 30MODELS           

mod.glmvul.0 <-glm.nb(vulpe ~ 1, data = data.112)
summary(mod.glmvul.0)

mod.glmvul.100 <- glm.nb(vulpe ~ scale(Populatie) + scale(pop2010to2019) +
                           + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                           + scale(perc_spverzi) + scale(perc_parcurb) + 
                           + scale(forMN) +
                           + scale(urbMN) + scale(Conectivitate_natur.) + relief, data = data.112)
summary(mod.glmvul.100)

mod.glmvul.101 <- glm.nb(vulpe ~ scale(Populatie) + scale(pop2010to2019) +
                           + scale(StuSup) +
                           + scale(PopAgri), data = data.112)
summary(mod.glmvul.101)

mod.glmvul.102 <- glm.nb(vulpe ~ scale(Suprafata) + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                           + scale(perc_spverzi) + scale(perc_parcurb) + 
                           + scale(forMN) +
                           + scale(urbMN) + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmvul.102)

mod.glmvul.1 <- glm.nb(vulpe ~ scale(Populatie) + 
                         + scale(D_strazi_urb_2018) + 
                         +scale(perc_spverzi) + scale(StuSup) + relief,
                       data = data.112)
summary(mod.glmvul.1)

mod.glmvul.2 <- glm.nb(vulpe ~ scale(Suprafata)+
                         + scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + scale(StuSup) + scale(PopAgri), data = data.112)
summary(mod.glmvul.2)

mod.glmvul.3 <- glm.nb(vulpe ~ scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + scale(PopAgri) +
                         + relief, data = data.112)
summary(mod.glmvul.3)

mod.glmvul.4 <- glm.nb(vulpe ~ scale(Populatie) +
                         + scale(rata_crestere_supr_construita) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + scale(PopAgri) +
                         + relief, data = data.112)
summary(mod.glmvul.4)

mod.glmvul.5 <- glm.nb(vulpe ~ scale(pop2010to2019) +
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + scale(StuSup) +  
                         + relief, data = data.112)
summary(mod.glmvul.5)

mod.glmvul.6 <- glm.nb(vulpe ~ scale(Populatie) + scale(pop2010to2019) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + scale(StuSup) + scale(PopAgri)   
                       , data = data.112)
summary(mod.glmvul.6)

mod.glmvul.7 <- glm.nb(vulpe ~ scale(pop2010to2019) +
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018)  
                       + scale(perc_spverzi) + scale(perc_parcurb) 
                       , data = data.112)
summary(mod.glmvul.7)

mod.glmvul.8 <- glm.nb(vulpe ~ scale(pop2010to2019) + scale(Suprafata)+
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + scale(StuSup) + scale(PopAgri) + 
                         + relief, data = data.112)
summary(mod.glmvul.8)

mod.glmvul.9 <- glm.nb(vulpe ~ scale(Suprafata)+
                         + scale(perc_spverzi) + scale(perc_parcurb) 
                       , data = data.112)
summary(mod.glmvul.9)

mod.glmvul.10 <- glm.nb(vulpe ~ scale(Populatie) 
                        + scale(perc_spverzi) + scale(perc_parcurb) + scale(StuSup) + scale(PopAgri) +
                          + relief, data = data.112)
summary(mod.glmvul.10)

mod.glmvul.11 <- glm.nb(vulpe ~ + scale(forMN) +
                          + scale(Conectivitate_natur.) + relief, data = data.112)
summary(mod.glmvul.11)

mod.glmvul.12 <- glm.nb(vulpe ~ + scale(agMN) +
                          + scale(Conectivitate_natur.) + relief, data = data.112)
summary(mod.glmvul.12)

mod.glmvul.13 <- glm.nb(vulpe ~ scale(agMN) +
                          + scale(urbMN) +relief, data = data.112)
summary(mod.glmvul.13)

mod.glmvul.14 <- glm.nb(vulpe ~ + scale(forMN) +
                          + scale(urbMN) +relief, data = data.112)
summary(mod.glmvul.14)

mod.glmvul.15 <- glm.nb(vulpe ~ scale(forMN) +
                          + scale(urbMN) + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmvul.15)

mod.glmvul.16 <- glm.nb(vulpe ~ scale(Proc_pad) +
                          + scale(Conectivitate_natur.) +relief, data = data.112)
summary(mod.glmvul.16)

mod.glmvul.17 <- glm.nb(vulpe ~ scale(Proc_agri) +
                          + scale(Conectivitate_natur.) +relief, data = data.112)
summary(mod.glmvul.17)

mod.glmvul.18 <- glm.nb(vulpe ~ scale(Proc_agri) +
                          + scale(urbMN) +relief, data = data.112)
summary(mod.glmvul.18)

mod.glmvul.19 <- glm.nb(vulpe ~ scale(Proc_pad) +
                          + scale(urbMN) +relief, data = data.112)
summary(mod.glmvul.19)

mod.glmvul.20 <- glm.nb(vulpe ~ scale(Proc_pad) +
                          + scale(urbMN) + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmvul.20)

mod.glmvul.21 <- glm.nb(vulpe ~ scale(Populatie) + 
                          + scale(D_strazi_urb_2018) + 
                          + scale(urbMN) + scale(Conectivitate_natur.) + scale(StuSup) + relief,
                        data = data.112)
summary(mod.glmvul.21)

mod.glmvul.22 <- glm.nb(vulpe ~ scale(Suprafata)+
                          + scale(D_strazi_urb_2018) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(forMN) +
                          + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri), data = data.112)
summary(mod.glmvul.22)

mod.glmvul.23 <- glm.nb(vulpe ~ scale(D_strazi_urb_2018) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(agMN) +
                          + scale(urbMN) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmvul.23)

mod.glmvul.24 <- glm.nb(vulpe ~ scale(Populatie) + 
                          + scale(rata_crestere_supr_construita) +
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(Proc_pad) +
                          + scale(urbMN) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmvul.24)

mod.glmvul.25 <- glm.nb(vulpe ~ scale(pop2010to2019) +
                          + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                          + scale(Proc_agri) +
                          + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmvul.25)

mod.glmvul.26 <- glm.nb(vulpe ~ scale(Populatie) + scale(pop2010to2019) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(Proc_agri) +
                          + scale(urbMN) + scale(StuSup), data = data.112)
summary(mod.glmvul.26)

mod.glmvul.27 <- glm.nb(vulpe ~ scale(pop2010to2019) +
                          + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(agMN) +
                          + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmvul.27)

mod.glmvul.28 <- glm.nb(vulpe ~ scale(pop2010to2019) + scale(Suprafata)+
                          + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                          + scale(forMN) +
                          + scale(urbMN) + scale(StuSup) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmvul.28)

mod.glmvul.29 <- glm.nb(vulpe ~ scale(Suprafata)+
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(Proc_pad) +
                          + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmvul.29)

mod.glmvul.30 <- glm.nb(vulpe ~ scale(Populatie) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(Proc_pad) +
                          + scale(urbMN) + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri) + relief, 
                        data = data.112)
summary(mod.glmvul.30)


#modele cu erori 727


#MODEL SELECTION
model.sel(mod.glmvul.0, mod.glmvul.100, mod.glmvul.101, mod.glmvul.102,
          mod.glmvul.1, mod.glmvul.2, mod.glmvul.3, mod.glmvul.4, mod.glmvul.5, 
          mod.glmvul.6,  mod.glmvul.8, mod.glmvul.9, mod.glmvul.10, 
          mod.glmvul.11, mod.glmvul.12, mod.glmvul.13, mod.glmvul.14, mod.glmvul.15, 
          mod.glmvul.16, mod.glmvul.17, mod.glmvul.18, mod.glmvul.19, mod.glmvul.20, 
          mod.glmvul.21, mod.glmvul.22, mod.glmvul.23, mod.glmvul.24, mod.glmvul.25, 
          mod.glmvul.26, mod.glmvul.28, mod.glmvul.29, mod.glmvul.30)
msglmvul <- model.sel(mod.glmvul.0, mod.glmvul.100, mod.glmvul.101, mod.glmvul.102,
                      mod.glmvul.1, mod.glmvul.2, mod.glmvul.3, mod.glmvul.4, mod.glmvul.5, 
                      mod.glmvul.6,  mod.glmvul.8, mod.glmvul.9, mod.glmvul.10, 
                      mod.glmvul.11, mod.glmvul.12, mod.glmvul.13, mod.glmvul.14, mod.glmvul.15, 
                      mod.glmvul.16, mod.glmvul.17, mod.glmvul.18, mod.glmvul.19, mod.glmvul.20, 
                      mod.glmvul.21, mod.glmvul.22, mod.glmvul.23, mod.glmvul.24, mod.glmvul.25, 
                      mod.glmvul.26, mod.glmvul.28, mod.glmvul.29, mod.glmvul.30)
#str(ms)
write.csv(msglmvul, "model_selection_all_VUL.csv")

# model averaging
avgmod <- model.avg(msglmvul, cumsum(weight)<= 0.95)
summary(avgmod)




#################################################################
#### MODELS (GLM) FOR SARPE - NULL MODEL + ALL MODEL + 30MODELS           

mod.glmsar.0 <-glm.nb(sarpe ~ 1, data = data.112)
summary(mod.glmsar.0)

mod.glmsar.100 <- glm.nb(sarpe ~ scale(Populatie) + scale(pop2010to2019) +
                           + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                           + scale(perc_spverzi) + scale(perc_parcurb) + 
                           + scale(forMN) +
                           + scale(urbMN) + scale(Conectivitate_natur.) + relief, data = data.112)
summary(mod.glmsar.100)

mod.glmsar.101 <- glm.nb(sarpe ~ scale(Populatie) + scale(pop2010to2019) +
                           + scale(StuSup) +
                           + scale(PopAgri), data = data.112)
summary(mod.glmsar.101)

mod.glmsar.102 <- glm.nb(sarpe ~ scale(Suprafata) + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                           + scale(perc_spverzi) + scale(perc_parcurb) + 
                           + scale(forMN) +
                           + scale(urbMN) + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmsar.102)

mod.glmsar.1 <- glm.nb(sarpe ~ scale(Populatie) + 
                         + scale(D_strazi_urb_2018) + 
                         +scale(perc_spverzi) + scale(StuSup) + relief,
                       data = data.112)
summary(mod.glmsar.1)

mod.glmsar.2 <- glm.nb(sarpe ~ scale(Suprafata)+
                         + scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + scale(StuSup) + scale(PopAgri), data = data.112)
summary(mod.glmsar.2)

mod.glmsar.3 <- glm.nb(sarpe ~ scale(D_strazi_urb_2018) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + scale(PopAgri) +
                         + relief, data = data.112)
summary(mod.glmsar.3)

mod.glmsar.4 <- glm.nb(sarpe ~ scale(Populatie) +
                         + scale(rata_crestere_supr_construita) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + scale(PopAgri) +
                         + relief, data = data.112)
summary(mod.glmsar.4)

mod.glmsar.5 <- glm.nb(sarpe ~ scale(pop2010to2019) +
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + scale(StuSup) +  
                         + relief, data = data.112)
summary(mod.glmsar.5)

mod.glmsar.6 <- glm.nb(sarpe ~ scale(Populatie) + scale(pop2010to2019) + 
                         + scale(perc_spverzi) + scale(perc_parcurb) + scale(StuSup) + scale(PopAgri)   
                       , data = data.112)
summary(mod.glmsar.6)

mod.glmsar.7 <- glm.nb(sarpe ~ scale(pop2010to2019) +
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018)  
                       + scale(perc_spverzi) + scale(perc_parcurb) 
                       , data = data.112)
summary(mod.glmsar.7)

mod.glmsar.8 <- glm.nb(sarpe ~ scale(pop2010to2019) + scale(Suprafata)+
                         + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + scale(StuSup) + scale(PopAgri) + 
                         + relief, data = data.112)
summary(mod.glmsar.8)

mod.glmsar.9 <- glm.nb(sarpe ~ scale(Suprafata)+
                         + scale(perc_spverzi) + scale(perc_parcurb) 
                       , data = data.112)
summary(mod.glmsar.9)

mod.glmsar.10 <- glm.nb(sarpe ~ scale(Populatie) 
                        + scale(perc_spverzi) + scale(perc_parcurb) + scale(StuSup) + scale(PopAgri) +
                          + relief, data = data.112)
summary(mod.glmsar.10)

mod.glmsar.11 <- glm.nb(sarpe ~ + scale(forMN) +
                          + scale(Conectivitate_natur.) + relief, data = data.112)
summary(mod.glmsar.11)

mod.glmsar.12 <- glm.nb(sarpe ~ + scale(agMN) +
                          + scale(Conectivitate_natur.) + relief, data = data.112)
summary(mod.glmsar.12)

mod.glmsar.13 <- glm.nb(sarpe ~ scale(agMN) +
                          + scale(urbMN) +relief, data = data.112)
summary(mod.glmsar.13)

mod.glmsar.14 <- glm.nb(sarpe ~ + scale(forMN) +
                          + scale(urbMN) +relief, data = data.112)
summary(mod.glmsar.14)

mod.glmsar.15 <- glm.nb(sarpe ~ scale(forMN) +
                          + scale(urbMN) + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmsar.15)

mod.glmsar.16 <- glm.nb(sarpe ~ scale(Proc_pad) +
                          + scale(Conectivitate_natur.) +relief, data = data.112)
summary(mod.glmsar.16)

mod.glmsar.17 <- glm.nb(sarpe ~ scale(Proc_agri) +
                          + scale(Conectivitate_natur.) +relief, data = data.112)
summary(mod.glmsar.17)

mod.glmsar.18 <- glm.nb(sarpe ~ scale(Proc_agri) +
                          + scale(urbMN) +relief, data = data.112)
summary(mod.glmsar.18)

mod.glmsar.19 <- glm.nb(sarpe ~ scale(Proc_pad) +
                          + scale(urbMN) +relief, data = data.112)
summary(mod.glmsar.19)

mod.glmsar.20 <- glm.nb(sarpe ~ scale(Proc_pad) +
                          + scale(urbMN) + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmsar.20)

mod.glmsar.21 <- glm.nb(sarpe ~ scale(Populatie) + 
                          + scale(D_strazi_urb_2018) + 
                          + scale(urbMN) + scale(Conectivitate_natur.) + scale(StuSup) + relief,
                        data = data.112)
summary(mod.glmsar.21)

mod.glmsar.22 <- glm.nb(sarpe ~ scale(Suprafata)+
                          + scale(D_strazi_urb_2018) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(forMN) +
                          + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri), data = data.112)
summary(mod.glmsar.22)

mod.glmsar.23 <- glm.nb(sarpe ~ scale(D_strazi_urb_2018) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(agMN) +
                          + scale(urbMN) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmsar.23)

mod.glmsar.24 <- glm.nb(sarpe ~ scale(Populatie) + 
                          + scale(rata_crestere_supr_construita) +
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(Proc_pad) +
                          + scale(urbMN) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmsar.24)

mod.glmsar.25 <- glm.nb(sarpe ~ scale(pop2010to2019) +
                          + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                          + scale(Proc_agri) +
                          + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmsar.25)

mod.glmsar.26 <- glm.nb(sarpe ~ scale(Populatie) + scale(pop2010to2019) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(Proc_agri) +
                          + scale(urbMN) + scale(StuSup), data = data.112)
summary(mod.glmsar.26)

mod.glmsar.27 <- glm.nb(sarpe ~ scale(pop2010to2019) +
                          + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(agMN) +
                          + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmsar.27)

mod.glmsar.28 <- glm.nb(sarpe ~ scale(pop2010to2019) + scale(Suprafata)+
                          + scale(rata_crestere_supr_construita) + scale(D_strazi_urb_2018) + 
                          + scale(forMN) +
                          + scale(urbMN) + scale(StuSup) + scale(PopAgri) + relief, data = data.112)
summary(mod.glmsar.28)

mod.glmsar.29 <- glm.nb(sarpe ~ scale(Suprafata)+
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(Proc_pad) +
                          + scale(Conectivitate_natur.), data = data.112)
summary(mod.glmsar.29)

mod.glmsar.30 <- glm.nb(sarpe ~ scale(Populatie) + 
                          + scale(perc_spverzi) + scale(perc_parcurb) + 
                          + scale(Proc_pad) +
                          + scale(urbMN) + scale(Conectivitate_natur.) + scale(StuSup) + scale(PopAgri) + relief, 
                        data = data.112)
summary(mod.glmsar.30)

#MODEL SELECTION
model.sel(mod.glmsar.0, mod.glmsar.100, mod.glmsar.101, mod.glmsar.102,
          mod.glmsar.1, mod.glmsar.2, mod.glmsar.3, mod.glmsar.4, mod.glmsar.5, 
          mod.glmsar.6, mod.glmsar.7, mod.glmsar.8, mod.glmsar.9, mod.glmsar.10, 
          mod.glmsar.11, mod.glmsar.12, mod.glmsar.13, mod.glmsar.14, mod.glmsar.15, 
          mod.glmsar.16, mod.glmsar.17, mod.glmsar.18, mod.glmsar.19, mod.glmsar.20, 
          mod.glmsar.21, mod.glmsar.22, mod.glmsar.23, mod.glmsar.24, mod.glmsar.25, 
          mod.glmsar.26, mod.glmsar.27, mod.glmsar.28, mod.glmsar.29, mod.glmsar.30)
msglmsar <- model.sel(mod.glmsar.0, mod.glmsar.100, mod.glmsar.101, mod.glmsar.102,
                      mod.glmsar.1, mod.glmsar.2, mod.glmsar.3, mod.glmsar.4, mod.glmsar.5, 
                      mod.glmsar.6, mod.glmsar.7, mod.glmsar.8, mod.glmsar.9, mod.glmsar.10, 
                      mod.glmsar.11, mod.glmsar.12, mod.glmsar.13, mod.glmsar.14, mod.glmsar.15, 
                      mod.glmsar.16, mod.glmsar.17, mod.glmsar.18, mod.glmsar.19, mod.glmsar.20, 
                      mod.glmsar.21, mod.glmsar.22, mod.glmsar.23, mod.glmsar.24, mod.glmsar.25, 
                      mod.glmsar.26, mod.glmsar.27, mod.glmsar.28, mod.glmsar.29, mod.glmsar.30)
#str(ms)
write.csv(msglmsar, "model_selection_all_SAR.csv")

# model averaging
avgmod <- model.avg(msglmsar, cumsum(weight)<= 0.95)
summary(avgmod)


