
## Working directory must be set to where downloaded supplementary material it stored
setwd("/Users/gdunshea/andypreds/Bauman et. al. 2019 - Supplementary Material")

library(lattice)
library(reshape)
library(nlme)
library(MCMCglmm)
library(MuMIn)

## Loading biomass removal models
# standardized mass loss in grams
load(file="tba.rda")
load(file="tba1.rda")
load(file="tba2.rda")
load(file="tba3.rda")
model.sel( tba,tba1,tba2,tba3, rank="DIC")
summary(tba)

# standardized mass loss in terms of proportion of the thallus removed
load(file="tbap.rda")
load(file="tbap1.rda")
load(file="tbap2.rda")
load(file="tbap3.rda")
model.sel( tbap,tbap1,tbap2,tbap3, rank="DIC")
summary(tbap)


## Loading all bites per event models
load(file = "firstb.rda")
load(file = "firstb1.rda")
load(file = "firstb2.rda")
load(file = "firstb3.rda")
load(file = "firstb4.rda")
load(file = "firstb5.rda")
load(file = "firstb6.rda")
load(file = "firstb7.rda")
load(file = "firstb8.rda")
load(file = "firstb9.rda")
load(file = "firstb10.rda")
load(file = "firstb11.rda")
load(file = "firstb12.rda")
load(file = "firstb13.rda")
load(file = "firstb14.rda")
load(file = "firstb141.rda")
load(file = "firstb142.rda")
load(file = "firstb143.rda")
load(file = "firstb1411.rda")
load(file = "firstb1412.rda")
load(file = "firstb1413.rda")
load(file = "firstb1414.rda")
model.sel(firstb, firstb1, firstb2,  firstb3, firstb4, firstb5, firstb6,firstb7,firstb8,firstb9,firstb10
          ,firstb11,firstb12,firstb13,firstb14,firstb141,firstb142,firstb143, firstb1411, firstb1412, firstb1413,firstb1414, rank="DIC")

summary(firstb143)


### Loading all feeding rates models
load(file = "sv.rda")
load(file = "sv1.rda")
load(file = "sv2.rda")
load(file = "sv3.rda")
model.sel(sv, sv1, sv2,sv3,rank="DIC")
summary(sv)



## loading all group size models:

load(file = "pgrpsv.rda")
load(file = "pgrpsv1.rda")
load(file = "pgrpsv2.rda")
load(file = "pgrpsv3.rda")
load(file = "pgrpsv4.rda")
load(file = "pgrpsv5.rda")
load(file = "pgrpsv11.rda")
load(file = "pgrpsv12.rda")

model.sel(pgrpsv, pgrpsv1, pgrpsv2,pgrpsv3, pgrpsv4, pgrpsv5, pgrpsv11, pgrpsv12, rank="DIC")
summary(pgrpsv)
     