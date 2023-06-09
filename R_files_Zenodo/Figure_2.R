
## FIGURE 2A ##

# G-test goodness-of-fit
library(DescTools)

# DC vs FLO
observed=c(8,22)    
expected=c(0.5, 0.5)
GTest(x=observed,p=expected,correct="none")

# DC vs DC+Lin
observed=c(24,6)    
expected=c(0.5, 0.5)
GTest(x=observed,p=expected,correct="none")

# DC vs DC+(Lin:JL)
observed=c(6,16)    
expected=c(0.5, 0.5)
GTest(x=observed,p=expected,correct="none")

# Left vs right
observed=c(41,41)    
expected=c(0.5, 0.5)
GTest(x=observed,p=expected,correct="none") 



## FIGURE 2B ##
# Proportion of eggs laid
# Wilcoxon Signed-Rank Test
library(MASS) 

# DC vs FLO, file: dc_flo.csv
wilcox.test(dc_flo$DC, dc_flo$FLO, paired=TRUE) 

# DC vs DC+Lin, file: dc_lin.csv
wilcox.test(dc_lin$DC, dc_lin$DCLin, paired=TRUE) 

# DC vs DC+(Lin:JL), file: dc_mix.csv
wilcox.test(dc_mix$DC, dc_mix$DCMix, paired=TRUE)



## FIGURE 2C ## 
# Proportion of eggs laid, file: oviposition.csv

# MIRABILIS
# Female fecundity: flowering and non-flowering treatments
mirabilisovip<-ovipositiondata %>%
  filter(plant_species=="mirabilis")
mirab_totaleggs<-aov(total_eggs ~ flowering_status, data=mirabilisovip)
summary(mirab_totaleggs)

# Oviposition: flowering vs non-flowering treatments
mirbovip<-ovipositiondata %>%
  filter(plant_species=="mirabilis")

mirbtreatments<-paste(mirbovip$plant_population, mirbovip$flowering_status)
mirbKW<-cbind(mirbovip, mirbtreatments)

kruskal.test(proportion_laid ~ mirbtreatments, data=mirbKW)


# O. HARRINGTONII
#Female fecundity: flowering and non-flowering treatments
harryovip<-ovipositiondata %>%
  filter(plant_species=="harringtonii")
oh_totaleggs<-aov(total_eggs ~ flowering_status, data=harryovip)
summary(oh_totaleggs)

# Oviposition: flowering vs non-flowering treatments

#Sorting treatments
bloomingFlo<-ovipositiondata %>%
  filter(plant_species=="harringtonii" & plant_population=="FLO" & flowering_status=="Y" & total_eggs > 0)

nbloomingFlo<-ovipositiondata %>%
  filter(plant_species=="harringtonii" & plant_population=="FLO" & flowering_status=="N" & total_eggs > 0)

bloomingB<-ovipositiondata %>%
  filter(plant_species=="harringtonii" & plant_population=="BLOOM" & flowering_status=="Y" & total_eggs > 0)

nbloomingB<-ovipositiondata %>%
  filter(plant_species=="harringtonii" & plant_population=="BLOOM" & flowering_status=="N" & total_eggs > 0)

harringtoniiovip<-ovipositiondata %>%
  filter(plant_species=="harringtonii")

o.h.treatments<-paste(harringtoniiovip$plant_population, harringtoniiovip$flowering_status)
forKW<-cbind(harringtoniiovip, o.h.treatments)

#Test
kruskal.test(proportion_laid ~ o.h.treatments, data=forKW)
