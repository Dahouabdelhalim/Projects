## FIGURE S2A ##

# Proportion survival, file: survival_plants.csv
xtabs(~dead+diet,data=survival_plants)
n<-table(survival_plants$diet) 
n
xtabs(~diet[dead=="1"],data=survival_plants)/n 

conting<-xtabs(~diet[dead=="1"],data=survival_plants)
chi<-chisq.test(conting)
conting/(chi$expected)
(chi$expected)/n



## FIGURE S2B ##
# Larval masses, file: larval masses_plants.csv
kruskal.test(larvalmasses_plants$weight~larvalmasses_plants$diet,data=larvalmasses_plants)

# pairwise comparisons
library("dunn.test")
library(FSA)
dunnTest(larvalmasses_plants$weight~larvalmasses_plants$diet,data=larvalmasses_plants,
         method="bonferroni")