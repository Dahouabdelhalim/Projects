## FIGURE 5 ##
# file: phenolicsDC.csv

library(vegan)

phenolics.dc.matrix<-as.matrix(phenolics.dc[,3:7]) ##response variables in a sample x species matrix
phenolics.dc.matrix

phenolics.dc.permanova<-adonis2(phenolics.dc.matrix~tissue, data=phenolics.dc, 
                                permutations = 999, method="bray")
phenolics.dc.permanova



# Comparisons

# OENOTHEIN B, file: oenotheinB.csv
kruskal.test(oenotheinB$abundance~oenotheinB$tissue,data=oenotheinB)
# pairwise comparisons
library("dunn.test")
library(FSA)
dunnTest(oenotheinB$abundance~oenotheinB$tissue,data=oenotheinB,method="holm")

# QUERCETIN, file: quercetin.csv
kruskal.test(quercetin$abundance~quercetin$tissue,data=quercetin)
# pairwise comparisons
dunnTest(quercetin$abundance~quercetin$tissue,data=quercetin,method="holm")

# OENOTHEIN A, file: oenotheinA.csv
kruskal.test(oenotheinA$abundance~oenotheinA$tissue,data=oenotheinA)
# pairwise comparisons
dunnTest(oenotheinA$abundance~oenotheinA$tissue,data=oenotheinA,method="holm")

# KAEMPFEROL, file: kaempferol.csv
kruskal.test(kaempferol$abundance~kaempferol$tissue,data=kaempferol)
# pairwise comparisons
dunnTest(kaempferol$abundance~kaempferol$tissue,data=kaempferol,method="holm")

# MYRICETIN, file: myricetin.csv
kruskal.test(myricetin$abundance~myricetin$tissue,data=myricetin)
# pairwise comparisons
dunnTest(myricetin$abundance~myricetin$tissue,data=myricetin,method="holm")

# NARINGENIN, file: naringenin.csv
kruskal.test(naringenin$abundance~naringenin$tissue,data=naringenin)
# pairwise comparisons
dunnTest(naringenin$abundance~naringenin$tissue,data=naringenin,method="holm")

# ALL COMPOUNDS COMBINED, file: total.csv
kruskal.test(total$abundance~total$tissue,data=total)
# pairwise comparisons
dunnTest(total$abundance~total$tissue,data=total,method="holm")
