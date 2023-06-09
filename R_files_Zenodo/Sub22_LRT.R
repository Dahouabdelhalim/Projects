###### Library ######

library(diversitree)
library(Matrix)
library(stats4)

###### MO vs M1 ######

print("M0 vs M1, Eurypterygii")
load("fit_musse_neot_M0_YK2021.Robj")
test.fit<-fit.musse.neot
load("fit_musse_neot_M1_YK2021.Robj")
null.fit<-fit.musse.neot
rm(fit.musse.neot)
## Likelihood ratio test
anova(test.fit,null.fit)
rm(test.fit)
rm(null.fit)

print("M0 vs M1, Otophysi")
load("fit_musse_otop_M0_YK2021.Robj")
test.fit<-fit.musse.otop
load("fit_musse_otop_M1_YK2021.Robj")
null.fit<-fit.musse.otop
rm(fit.musse.otop)
## Likelihood ratio test
anova(test.fit,null.fit)
rm(test.fit)
rm(null.fit)

###### M2 vs M0' ######

print("M2 vs M0', Eurypterygii")
load("fit_musse_neot_M2_YK2021.Robj")
test.fit<-fit.musse.neot
load("fit_musse_neot_M0_prime_YK2021.Robj")
null.fit<-fit.musse.neot
rm(fit.musse.neot)
## Likelihood ratio test
anova(test.fit,null.fit)
rm(test.fit)
rm(null.fit)

print("M2 vs M0', Eupercaria")
load("fit_musse_neot_groupA_M2_YK2021.Robj")
test.fit<-fit.musse.neot.groupA
load("fit_musse_neot_null_groupA_M0_prime_YK2021.Robj")
null.fit<-fit.musse.neot.null.groupA
rm(fit.musse.neot)
## Likelihood ratio test
anova(test.fit,null.fit)
rm(test.fit)
rm(null.fit)

print("M2 vs M0', Otophysi")
load("fit_musse_otop_M2_YK2021.Robj")
test.fit<-fit.musse.otop
load("fit_musse_otop_M0_prime_YK2021")
null.fit<-fit.musse.otop
rm(fit.musse.otop)
## Likelihood ratio test
anova(test.fit,null.fit)
rm(test.fit)
rm(null.fit)


###### M2 vs M3 ######

print("M2 vs M3, Eurypterygii")
load("fit_musse_neot_M2_YK2021.Robj")
test.fit<-fit.musse.neot
load("fit_musse_neot_M3_YK2021.Robj")
null.fit<-fit.musse.neot
rm(fit.musse.neot)
## Likelihood ratio test
anova(test.fit,null.fit)
rm(test.fit)
rm(null.fit)

print("M2 vs M3, Otophysi")
load("fit_musse_otop_M2_YK2021.Robj")
test.fit<-fit.musse.otop
load("fit_musse_otop_M3_YK2021.Robj")
null.fit<-fit.musse.otop
rm(fit.musse.otop)
## Likelihood ratio test
anova(test.fit,null.fit)
rm(test.fit)
rm(null.fit)