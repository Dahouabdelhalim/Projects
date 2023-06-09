
########## Scripts for SEM analysis performed ##########

## NOTE: 
# Analyses were performed several years ago, using piecewiseSEM package v. 1.0.1. 
# Main commands in new versions of this package have changed substantially (please see: 
# https://cran.r-project.org/web/packages/piecewiseSEM/piecewiseSEM.pdf). 
# Here we provide our original commands in order to ensure reproducibility. 
# We encourage interested readers to check: 
# https://mran.microsoft.com/snapshot/2015-12-10/web/packages/piecewiseSEM/piecewiseSEM.pdf

## NOTE 2: 
# Variables in the Database file (DATA_SEM&MNKA.xlsx) are provided with long names, following the 
# ‘Best practices for computer code and statistics’ detailed in Author guidelines:
# https://esajournals.onlinelibrary.wiley.com/hub/journal/15577015/resources/author-guidelines-ecm 
# Equivalences of such variable names with those in the scripts below are:
# IDg:      PLOT
# PR:       log.TimePredExclus
# OD:       log.TimeOdegusExclus
# SM:       log.TimeSmallMammExclus
# PPt0:     Precipit.t0
# PPt1:     Precipit.t1
# ln.At0:   log.Cov_Shrubs.t0
# ln.At1:   log.Cov_Shrubs.t1
# ln.at0:   log.Cov_Ephemerals.t0
# ln.at1:   log.Cov_Ephemerals.t1
# ln.Ht1:   log.MNKA_Herbivores.t1
# ln.Ht0:   log.MNKA_Herbivores.t2
# In addition to the excel file with the database with long names, we also included a txt file (DATA_SEM.txt)
# with shorter variable names as used in the scripts below for uploading and directly executing the analyses.


# Calling the necessary libraries
library("nlme")
library("piecewiseSEM")

# Loading the Data Base
DatFJ2<-read.table("DATA_SEM.txt",sep="\\t",dec=",",header=T)


## M0: Initial theoretical SEM with all a priori expected effects

# Create list of models implied by the SEM:
M0<-list(
  lme(ln.At0~ln.At1+log(PPt0)+log(PPt1)+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~log(PPt1)+ln.at1+ln.Ht2, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~ln.at1+log(PPt0)+log(PPt1)+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~log(PPt1)+ln.Ht2, random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+ln.at1+SM+OD+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M0, DatFJ2) # Run model fitting
sem.coefs(M0, DatFJ2, standardize="scale") # Extract path coefficients 
sem.model.fits(M0) #Get R2 for individual models

########################################################################################

## M0b: As M0, but without log-transforming PPt and PPt1 for consistency with M16 and M16b (see below)

M0b<-list(
  lme(ln.At0~ln.At1+PPt0+PPt1+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~PPt1+ln.at1+ln.Ht2, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~ln.at1+PPt0+PPt1+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~PPt1+ln.Ht2, random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+ln.at1+SM+OD+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M0b, DatFJ2)
sem.coefs(M0b, DatFJ2, standardize="scale")
sem.model.fits(M0b)

########################################################################################

## M1: As M0, but removing non-significant paths 
M1<-list(
  lme(ln.At0~ln.At1+log(PPt0)+log(PPt1)+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~log(PPt1)+ln.at1, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~ln.at1+log(PPt0)+log(PPt1), random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~log(PPt1), random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+SM, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M1, DatFJ2, corr.errors = c("SM~~ln.Ht2","OD~~ln.Ht2","PR~~ln.Ht2"),add.vars = c("PR","OD"))
sem.coefs(M1, DatFJ2, standardize="scale")
sem.model.fits(M1)

########################################################################################

## M2: As M1, but adding the effects of PPt1 and PR on Ht1, of PR on At1, and a non-directed path (correlation) between at1 and Ht2, 
## as suggested by the sem.missing.paths() function 

M2<-list(
  lme(ln.At0~ln.At1+log(PPt0)+log(PPt1)+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~log(PPt1)+ln.at1+PR, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~ln.at1+log(PPt0)+log(PPt1), random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~log(PPt1), random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+log(PPt1)+SM+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M2, DatFJ2, corr.errors = c("SM~~ln.Ht2","OD~~ln.Ht2","PR~~ln.Ht2","at1~~Ht2"),add.vars = c("OD"))
sem.coefs(M2, DatFJ2, standardize="scale")
sem.model.fits(M2)

########################################################################################

## M3: As M2, bur adding the effect of at1 on Ht1, as suggested by the sem.missing.paths() function

M3<-list(
  lme(ln.At0~ln.At1+log(PPt0)+log(PPt1)+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~log(PPt1)+ln.at1+PR, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~ln.at1+log(PPt0)+log(PPt1), random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~log(PPt1), random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+ln.at1+log(PPt1)+SM+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M3, DatFJ2, corr.errors = c("SM~~ln.Ht2","OD~~ln.Ht2","PR~~ln.Ht2","at1~~Ht2"),add.vars = c("OD"))
sem.coefs(M3, DatFJ2, standardize="scale")
sem.model.fits(M3)

########################################################################################

## M4: as M3, but adding non-directed paths (correlations) between PPt0 and At1, and between at1 and Ht1, 
## as suggested by the sem.missing.paths() function 

M4<-list(
  lme(ln.At0~ln.At1+log(PPt0)+log(PPt1)+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~log(PPt1)+ln.at1+PR, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~ln.at1+log(PPt0)+log(PPt1), random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~log(PPt1), random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+ln.at1+log(PPt1)+SM+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M4, DatFJ2, corr.errors = c("SM~~ln.Ht2","OD~~ln.Ht2","PR~~ln.Ht2","at1~~Ht2","PPt0~~At1","PPt0~~at1","PPt0~~Ht1"),add.vars = c("OD"))
sem.coefs(M4, DatFJ2, standardize="scale")
sem.model.fits(M4)

########################################################################################

## M5: as M4, but removing the non-significant effect of PPt1 on At0 

M5<-list(
  lme(ln.At0~ln.At1+log(PPt0)+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~log(PPt1)+ln.at1+PR, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~ln.at1+log(PPt0)+log(PPt1), random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~log(PPt1), random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+ln.at1+log(PPt1)+SM+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M5, DatFJ2, corr.errors = c("SM~~ln.Ht2","OD~~ln.Ht2","PR~~ln.Ht2","at1~~Ht2","PPt0~~At1","PPt0~~at1","PPt0~~Ht1"),add.vars = c("OD"))
sem.coefs(M5, DatFJ2, standardize="scale")
sem.model.fits(M5)

########################################################################################

### M6: as M4, but not log-transforming PPt nor PPt1

M6<-list(
  lme(ln.At0~ln.At1+PPt0+PPt1+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~PPt1+ln.at1+PR, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~ln.at1+PPt0+PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+ln.at1+log(PPt1)+SM+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M6, DatFJ2, corr.errors = c("SM~~ln.Ht2","OD~~ln.Ht2","PR~~ln.Ht2","ln.at1~~Ht2","PPt0~~ln.At1","PPt0~~ln.at1","PPt0~~ln.Ht1"),add.vars = c("OD"))
sem.coefs(M6, DatFJ2, standardize="scale")
sem.model.fits(M6)

########################################################################################

## M7: as M6, but removing the non-significant effect of PPt1 on At0 

M7<-list(
  lme(ln.At0~ln.At1+PPt0+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~PPt1+ln.at1+PR, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~ln.at1+PPt0+PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+ln.at1+log(PPt1)+SM+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M7, DatFJ2, corr.errors = c("SM~~ln.Ht2","OD~~ln.Ht2","PR~~ln.Ht2","at1~~Ht2","PPt0~~At1","PPt0~~at1","PPt0~~Ht1"),add.vars = c("OD"))
sem.coefs(M7, DatFJ2, standardize="scale")
sem.model.fits(M7)
# In spite of removing an non-significant effect, this model performs worst than M6

########################################################################################

## M8: as M6 but removing nin-significant correlation between PPt0 and ln.Ht2

M8<-list(
  lme(ln.At0~ln.At1+PPt0+PPt1+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~PPt1+ln.at1+PR, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~ln.at1+PPt0+PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+ln.at1+log(PPt1)+SM+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M8, DatFJ2, corr.errors = c("SM~~ln.Ht2","OD~~ln.Ht2","PR~~ln.Ht2","at1~~Ht2","PPt0~~At1","PPt0~~at1"),add.vars = c("OD"))
sem.coefs(M8, DatFJ2, standardize="scale")
sem.model.fits(M8)
# This model represent a simplified version of M6, but with a slightly higher AICc value

########################################################################################

## M9: as M8, but adding a marginally significant effect of OD on ln.Ht1, as suggested by the sem.missing.paths() function 

M9<-list(
  lme(ln.At0~ln.At1+PPt0+PPt1+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~PPt1+ln.at1+PR, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~ln.at1+PPt0+PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+ln.at1+log(PPt1)+SM+PR+OD, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M9, DatFJ2, corr.errors = c("SM~~ln.Ht2","OD~~ln.Ht2","PR~~ln.Ht2","at1~~Ht2","PPt0~~At1","PPt0~~at1","PPt0~~ln.Ht2"))
sem.coefs(M9, DatFJ2, standardize="scale")
sem.model.fits(M9)

########################################################################################

## M10: as M6, but adding an effect of PR on ln.At0, as suggested by the sem.missing.paths() function 

M10<-list(
  lme(ln.At0~ln.At1+PPt0+PPt1+ln.at0+ln.Ht1+PR, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~PPt1+ln.at1+PR, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~ln.at1+PPt0+PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+ln.at1+log(PPt1)+SM+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M10, DatFJ2, corr.errors = c("SM~~ln.Ht2","OD~~ln.Ht2","PR~~ln.Ht2","at1~~Ht2","PPt0~~At1","PPt0~~at1","PPt0~~ln.Ht2"),add.vars = c("OD"))
sem.coefs(M10, DatFJ2, standardize="scale")
sem.model.fits(M10)

########################################################################################

## M11: as M6, but removing the effect of ln.at1 on ln.at0, adding an effect of OD on ln.Ht1, and 
## turning the directed path between ln.at1 and ln.Ht1 into an undirected path (correlation)
## as suggested by the sem.missing.paths() function

M11<-list(
  lme(ln.At0~ln.At1+PPt0+PPt1+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~PPt1+ln.at1+PR, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~PPt0+PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+PPt1+SM+OD+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M11, DatFJ2, corr.errors = c("SM~~ln.Ht2","OD~~ln.Ht2","PR~~ln.Ht2","ln.at1~~Ht2","ln.at1~~Ht1","PPt0~~ln.At1","PPt0~~ln.at1","PPt0~~ln.Ht2"))
sem.coefs(M11, DatFJ2, standardize="scale")
sem.model.fits(M11)

########################################################################################

## M12: as M11, but removing the effect of PR on ln.Ht1

M12<-list(
  lme(ln.At0~ln.At1+PPt0+PPt1+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~PPt1+ln.at1+PR, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~PPt0+PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+PPt1+SM+OD, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M12, DatFJ2, corr.errors = c("SM~~ln.Ht2","OD~~ln.Ht2","PR~~ln.Ht2","ln.at1~~Ht2","ln.at1~~Ht1","PPt0~~ln.At1","PPt0~~ln.at1","PPt0~~ln.Ht2"))
sem.coefs(M12, DatFJ2, standardize="scale")
sem.model.fits(M12)

########################################################################################

## M13: as M11, but the directed path (effect) of PR on ln.At1 is changed by an undirected path (correlation) to maintain temporal coherence

M13<-list(
  lme(ln.At0~ln.At1+PPt0+PPt1+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~PPt1+ln.at1+PR, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~PPt0+PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+PPt1+SM+OD+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M13, DatFJ2, corr.errors = c("SM~~ln.Ht2","OD~~ln.Ht2","PR~~ln.Ht2","ln.at1~~Ht2","ln.at1~~Ht1","PPt0~~ln.At1","PPt0~~ln.at1","PPt0~~ln.Ht2","ln.At1~~PR"))
sem.coefs(M13, DatFJ2, standardize="scale")
sem.model.fits(M13)

########################################################################################

## M14: as M13, but removing non-directed paths between exogenous variables

M14<-list(
  lme(ln.At0~ln.At1+PPt0+PPt1+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~PPt1+ln.at1, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~PPt0+PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+PPt1+SM+OD+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M14, DatFJ2, corr.errors = c("ln.at1~~ln.Ht2","ln.at1~~ln.Ht1","PPt0~~ln.At1","PPt0~~ln.at1","ln.At1~~PR"))
sem.coefs(M14, DatFJ2, standardize="scale")
sem.model.fits(M14)

########################################################################################

## M15: as M14, but replacing the non-directed path (correlation) between ln.Ht2 and ln.at1 by an effect of the former on the latter

M15<-list(
  lme(ln.At0~ln.At1+PPt0+PPt1+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~PPt1+ln.at1, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~PPt0+PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~PPt1+ln.Ht2, random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+PPt1+SM+OD+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M15, DatFJ2, corr.errors = c("ln.at1~~ln.Ht1","PPt0~~ln.At1","PPt0~~ln.at1","ln.At1~~PR"))
sem.coefs(M15, DatFJ2, standardize="scale")
sem.model.fits(M15)
# No further improvement: the effect added is just marginally significant and the corresponding AIC-value higher

########################################################################################

## M16: as M14, but replacing the non-directed path (correlation) between ln.at1 and ln.Ht1 by an effect of the formen on the latter

M16<-list(
  lme(ln.At0~ln.At1+PPt0+PPt1+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~PPt1+ln.at1, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~PPt0+PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+ln.at1+PPt1+SM+OD+PR, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M16, DatFJ2, corr.errors = c("ln.at1~~ln.Ht2","PPt0~~ln.At1","PPt0~~ln.at1","ln.At1~~PR"))
sem.coefs(M16, DatFJ2, standardize="scale")
sem.model.fits(M16)

########################################################################################

## M17: as M16, but adding an effect of ln.at1 on ln.at0, as suggested by the sem.missing.paths() function, 
# and log-transforming PPt0 in that equation to improve residuals' behavior 

M17<-list(
  lme(ln.At0~ln.At1+PPt0+PPt1+ln.at0+ln.Ht1, random=~1 | IDg, data=DatFJ2),
  lme(ln.At1~PPt1+ln.at1, random=~1 | IDg, data=DatFJ2), 
  lme(ln.at0~ln.at1+log(PPt0)+PPt1, random=~1 | IDg, data=DatFJ2),
  lme(ln.at1~log(PPt1), random=~1 | IDg, data=DatFJ2),
  lme(ln.Ht1~ln.Ht2+ln.at1+PPt1+SM+OD, random=~1 | IDg, data=DatFJ2) 
)
sem.fit(M16b, DatFJ2, corr.errors = c("ln.at1~~ln.Ht2","ln.At1~~PR","PPt0~~ln.At1","PPt0~~ln.at1"),add.vars = c("PR"))
sem.coefs(M16b, DatFJ2, standardize="scale")
sem.model.fits(M16b)













