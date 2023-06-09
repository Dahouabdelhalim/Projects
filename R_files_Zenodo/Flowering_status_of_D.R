## Tests for flower stage effects on resistance for Dianthus pavonius

## 1. Fisher exact test of effect of stage on infeciton in A1B1 plants
## only plants with full tray data and flower status are used.

#Matrix of all A1B1 plants that were scored (2019 and 2020 flowering data combined)
mat<-matrix(c(220,163,30,8),2,2)
colnames(mat)<-c("H","D")
rownames(mat)<-c("V","F")
mat
fisher.test(mat)

## Matrix of all A1B1 plants scored in 2020
mat20<-matrix(c(139,73,9,1),2,2)
colnames(mat20)<-c("H","D")
rownames(mat20)<-c("V","F")
mat20
fisher.test(mat20)


## Matrix of all A1B1 plants scored in 2019
mat19<-matrix(c(81,90,21,7),2,2)
colnames(mat19)<-c("H","D")
rownames(mat19)<-c("V","F")
mat19
fisher.test(mat19)
