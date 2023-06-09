######################################################################################################################################
## Script by Arianne F. Messerman
## 13 April, 2021
## This script examines relative change in body size metrics between spotted and marbled salamanders (Ambystoma maculatum and A. opacum)
##   during a semi-natural outdoor experiment.
## Please see the associated manuscript for full study details: 
##   Messerman, A.F. and M. Leal. Submitted. The contributions of individual traits to survival among terrestrial juvenile 
##   pond-breeding salamanders. Ecology.
######################################################################################################################################

#setwd()

size<-read.csv("Phys_Surv_Change_Size.csv", header=TRUE)
str(size)

library(car)

mod.mass<-lm(log(delta.mass)~spp, data=size)
summary(mod.mass)
anova(mod.mass)
Anova(mod.mass, type="II")#0.617

mod.svl<-lm(log(delta.svl)~spp, data=size)
summary(mod.svl)
anova(mod.svl)
Anova(mod.svl, type="II")#0.152

mod.tl<-lm(log(delta.tl+0.01)~spp, data=size)
summary(mod.tl)
anova(mod.tl)
Anova(mod.tl, type="II")#0.146

AMMA<-subset(size, size$spp=="AMMA")
AMOP<-subset(size, size$spp=="AMOP")

#Species Mass Change
mean(AMMA$delta.mass)#4.629
sd(AMMA$delta.mass)#4.607
mean(AMOP$delta.mass)#3.476
sd(AMOP$delta.mass)#2.976
mean(size$delta.mass)#4.008
sd(size$delta.mass)#3.780

#Species SVL Change
mean(AMMA$delta.svl)#17.000
sd(AMMA$delta.svl)#11.939
mean(AMOP$delta.svl)#10.964
sd(AMOP$delta.svl)#3.865
mean(size$delta.svl)#13.750
sd(size$delta.svl)#8.939

#Species TL Change
mean(AMMA$delta.tl)#27.818
sd(AMMA$delta.tl)#16.606
mean(AMOP$delta.tl)#17.071
sd(AMOP$delta.tl)#8.931
mean(size$delta.tl)#22.031
sd(size$delta.tl)#13.879
