###Food-Egg-fecundity experiment Matteo
##test whether there are differences in the occurences of specific maintentance terms in control versus treatment

#Food
mat1<-matrix(c(0,1,29-0,43-1),2)
mat2<-matrix(c(0,1,29-0,43-1),2)
mat3<-matrix(c(3,5,29-3,43-5),2)
mat6<-matrix(c(0,3,29-0,43-3),2)
mat7<-matrix(c(0,0,29-0,43-0),2)#na
mat8<-matrix(c(0,1,29-0,43-1),2)
mat9<-matrix(c(0,0,29-0,43-0),2)#na
mat10<-matrix(c(1,1,29-1,43-1),2)
mat11<-matrix(c(0,1,29-0,43-1),2)
mat12<-matrix(c(0,1,29-0,43-1),2)
mat15<-matrix(c(2,3,29-2,43-3),2)
mat16<-matrix(c(3,2,29-3,43-2),2)
###if not stated otherwise p=1 
chisq.test(mat1)
chisq.test(mat2)
chisq.test(mat3)
chisq.test(mat6) #X-squared = 0.7255, df = 1, p-value = 0.3943
chisq.test(mat7)#na
chisq.test(mat8)
chisq.test(mat9)#na
chisq.test(mat10)
chisq.test(mat11)
chisq.test(mat12)
chisq.test(mat15)
chisq.test(mat16) #X-squared = 0.21113, df = 1, p-value = 0.6459


#Egg
mat17<-matrix(c(8,5,194-8,178-5),2)	
chisq.test(mat17)#X-squared = 0.16579, df = 1, p-value = 0.6839
mat18<-matrix(c(9,6,194-9,178-6),2)	
chisq.test(mat18)#X-squared = 0.12775, df = 1, p-value = 0.7208
mat19<-matrix(c(19,11,194-19,178-11),2)	
chisq.test(mat19)#X-squared = 1.1842, df = 1, p-value = 0.2765
mat22<-matrix(c(2,4,194-2,178-4),2)	
chisq.test(mat22)#X-squared = 0.26861, df = 1, p-value = 0.6043
mat23<-matrix(c(0,0,194-0,178-0),2)	
chisq.test(mat23)#na
mat24<-matrix(c(1,0,194-1,178-0),2)	
chisq.test(mat24)#X-squared = 3.1613e-28, df = 1, p-value = 1
mat25<-matrix(c(0,0,194-0,178-0),2)	
chisq.test(mat25)#na
mat26<-matrix(c(17,2,194-17,178-2),2)	
chisq.test(mat26)#################X-squared = 9.6568, df = 1, p-value = 0.001887
mat27<-matrix(c(2,0,194-2,178-0),2)	
chisq.test(mat27)#X-squared = 0.42071, df = 1, p-value = 0.5166
mat28<-matrix(c(11,4,194-11,178-4),2)	
chisq.test(mat28)#X-squared = 1.9956, df = 1, p-value = 0.1578
mat31<-matrix(c(29,24,194-29,178-24),2)	
chisq.test(mat31)#X-squared = 0.065246, df = 1, p-value = 0.7984
mat32<-matrix(c(9,3,194-9,178-3),2)	
chisq.test(mat32)#X-squared = 1.7345, df = 1, p-value = 0.1878
mat33<-matrix(c(9,12,194-9,178-12),2)	
chisq.test(mat33)#X-squared = 0.42617, df = 1, p-value = 0.5139
##FDR correction for Egg results
setwd("D:/OneDrive/studies/matteo/egg-food")
eggp<-read.table("p-values_Egg.txt", header=FALSE)
eggp
p.adjust(eggp$V2,method="fdr")
#0.8871111 0.8871111 0.6912500 0.8871111 1.0000000 0.0188700 0.8871111 0.6260000 0.8871111 0.6260000


#Zusammengefasst 
#Food
mat40<-matrix(c(6,14,29-6,43-14),2)
chisq.test(mat40)#X-squared = 0.69642, df = 1, p-value = 0.404
mat41<-matrix(c(2,3,29-2,43-3),2)
chisq.test(mat41)#X-squared = 8.5023e-30, df = 1, p-value = 1
mat42<-matrix(c(3,2,29-3,43-2),2)
chisq.test(mat42)#X-squared = 0.21113, df = 1, p-value = 0.6459
mat43<-matrix(c(1,2,29-1,43-2),2)
chisq.test(mat43)#X-squared = 2.0143e-30, df = 1, p-value = 1

#Egg
mat34<-matrix(c(57,34,194-57,178-34),2)	
chisq.test(mat34)
mat35<-matrix(c(29,24,194-29,178-24),2)	
chisq.test(mat35)#X-squared = 0.065246, df = 1, p-value = 0.7984
mat36<-matrix(c(9,3,194-9,178-3),2)	
chisq.test(mat36)#X-squared = 1.7345, df = 1, p-value = 0.1878
mat37<-matrix(c(9,12,194-9,178-12),2)	
chisq.test(mat37)#X-squared = 0.42617, df = 1, p-value = 0.5139
