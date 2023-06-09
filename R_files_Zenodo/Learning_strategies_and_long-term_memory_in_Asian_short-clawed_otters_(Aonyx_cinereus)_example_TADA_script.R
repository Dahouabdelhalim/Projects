#To load the NBDA package, you first need to install the package "devtools" in the usual way
#Next load it up as follows:
library(devtools)
#Then download and install my NBDA package from GitHub:
devtools::install_github("whoppitt/NBDA")
library(NBDA)

#Define ILVs of each individual in the network
sex1<-c(0,0,0,1,1)
age1<-c(13,11,11,12,12)

#Put age and sex into a columns
sex1<-cbind(sex1)
age1<-cbind(age1)

#And add it to the asoc vector
asoc1<-c("sex1","age1")

#Create association matrix of social network
amg1<-matrix(data=c(0.00, 0.60, 0.43, 0.50, 0.52,	
			  0.60, 0.00, 0.47, 0.51, 0.49,	
			  0.43, 0.47, 0.00, 0.42, 0.39,	
			  0.50, 0.51, 0.42, 0.00, 0.41,	
			  0.52, 0.49, 0.39, 0.41, 0.00), nrow=5)
#Create order of solve data objects
oag1t1<-c(2,1,4,5,3)
oag1t2<-c(1,2,4,5,3)
oag1t3<-c(3,2,4,5,1)
oag1t4<-c(1,2,3,4,5)
oag1t5<-c(1,2,3,4,5)

#Check the dimensions of the social network
dim(amg1)

#Now we want to specify asoc_ilv = asoc- the ILVs for the additive model and multi_ilv=asoc, the ILVs for the multiplicative model, though we do not have these in the same model together

#Make array with room for 6 networks full of 0s
assMatrix<-array(0,dim=c(5,5,10))
#Enter network into slot 1 for task 1
assMatrix[,,1]<-amg1
#And 1s in slot 6 for the group network
assMatrix[,,6]<-1
ng1t1 <-nbdaData(label="ng1t1",assMatrix=assMatrix, asoc_ilv=asoc1, multi_ilv=asoc1,orderAcq=oag1t1, idname=c("Flint","Magma","Fossil","Earth","Moon"),timeAcq=c(56,57,62,62,65),endTime=176)

#Make array with room for 6 networks full of 0s
assMatrix<-array(0,dim=c(5,5,10))
#Enter network into slot 2 for task 2
assMatrix[,,2]<-amg1
#And 1s in slot 6 for the group network
assMatrix[,,7]<-1
ng1t2 <-nbdaData(label="ng1t2",assMatrix=assMatrix, asoc_ilv=asoc1, multi_ilv=asoc1, orderAcq=oag1t2, idname=c("Flint","Magma","Fossil","Earth","Moon"),timeAcq=c(75,76,76,81,85),endTime=261)

#Make array with room for 6 networks full of 0s
assMatrix<-array(0,dim=c(5,5,10))
#Enter network into slot 3 for task 3
assMatrix[,,3]<-amg1
#And 1s in slot 6 for the group network
assMatrix[,,8]<-1
ng1t3 <-nbdaData(label="ng1t3",assMatrix=assMatrix, asoc_ilv=asoc1, multi_ilv=asoc1, orderAcq=oag1t3, idname=c("Flint","Magma","Fossil","Earth","Moon"),timeAcq=c(55,55,58,62,155),endTime=643)

#Make array with room for 6 networks full of 0s
assMatrix<-array(0,dim=c(5,5,10))
#Enter network into slot 4 for task 4
assMatrix[,,4]<-amg1
#And 1s in slot 6 for the group network
assMatrix[,,9]<-1
ng1t4 <-nbdaData(label="ng1t4",assMatrix=assMatrix, asoc_ilv=asoc1, multi_ilv=asoc1, orderAcq=oag1t4, idname=c("Flint","Magma","Fossil","Earth","Moon"),timeAcq=c(30,32,36,41,43),endTime=545)

#Make array with room for 6 networks full of 0s
assMatrix<-array(0,dim=c(5,5,10))
#Enter network into slot 5 for task 5
assMatrix[,,5]<-amg1
#And 1s in slot 6 for the group network
assMatrix[,,10]<-1
ng1t5 <-nbdaData(label="ng1t5",assMatrix=assMatrix, asoc_ilv=asoc1, multi_ilv=asoc1, orderAcq=oag1t5, idname=c("Flint","Magma","Fossil","Earth","Moon"),timeAcq=c(65,67,70,72,74),endTime=735)


#Repeat for group 2 (Newquay)
sex2<-c(0,1,1,1,0,0,1,1)
age2<-c(9,7,4,1,0.75,0.75,0.75,0.75)
sex2<-cbind(sex2)
age2<-cbind(age2)
asoc2<-c("sex2","age2")
amg2<-matrix(data=c(0.00, 0.67, 0.56, 0.52, 0.52, 0.52, 0.50, 0.45,
			  0.67, 0.00, 0.59, 0.57, 0.61, 0.54, 0.52, 0.53,	
			  0.56, 0.59, 0.00, 0.52, 0.43, 0.44, 0.42, 0.45,
     			  0.52, 0.57, 0.52, 0.00, 0.55, 0.55, 0.56, 0.51,	
			  0.52, 0.61, 0.43, 0.55, 0.00, 0.67, 0.67, 0.67,	
			  0.52, 0.54, 0.44, 0.55, 0.67, 0.00, 0.63, 0.68,
			  0.50, 0.52, 0.42, 0.56, 0.67, 0.63, 0.00, 0.59,
  			  0.45, 0.53, 0.45, 0.51, 0.67, 0.68, 0.59, 0.00), nrow=8)
oag2t1<-c(1,3,2,6,4)
oag2t2<-c(2,5,7,8,4,1,3)
oag2t3<-c(3,1,4,2,5,6,7,8)
oag2t4<-c(4,5,2,6,8,3,7,1)
oag2t5<-c(1,4,6,8,3,5,2,7)

dim(amg2)

assMatrix<-array(0,dim=c(8,8,10))
assMatrix[,,1]<-amg2
assMatrix[,,6]<-1
ng2t1 <-nbdaData(label="ng2t1",assMatrix=assMatrix, asoc_ilv=asoc2, multi_ilv=asoc2,orderAcq=oag2t1, idname=c("Tope","Jam","Pod","Meg","Rosie","Dot","Charlie","Biscuit"),timeAcq=c(20,77,81,155,165),endTime=406)

assMatrix<-array(0,dim=c(8,8,10))
assMatrix[,,2]<-amg2
assMatrix[,,7]<-1
ng2t2 <-nbdaData(label="ng2t2",assMatrix=assMatrix, asoc_ilv=asoc2, multi_ilv=asoc2, orderAcq=oag2t2, idname=c("Tope","Jam","Pod","Meg","Rosie","Dot","Charlie","Biscuit"),timeAcq=c(21,25,25,38,40,48,82,93),endTime=251)

assMatrix<-array(0,dim=c(8,8,10))
assMatrix[,,3]<-amg2
assMatrix[,,8]<-1
ng2t3 <-nbdaData(label="ng2t3",assMatrix=assMatrix, asoc_ilv=asoc2, multi_ilv=asoc2, orderAcq=oag2t3, idname=c("Tope","Jam","Pod","Meg","Rosie","Dot","Charlie","Biscuit"),timeAcq=c(11,13,16,62,124,127,138,202),endTime=518)

assMatrix<-array(0,dim=c(8,8,10))
assMatrix[,,4]<-amg2
assMatrix[,,9]<-1
ng2t4 <-nbdaData(label="ng2t4",assMatrix=assMatrix, asoc_ilv=asoc2, multi_ilv=asoc2, orderAcq=oag2t4, idname=c("Tope","Jam","Pod","Meg","Rosie","Dot","Charlie","Biscuit"),timeAcq=c(19,21,22,22,24,44,51,90),endTime=1213)

assMatrix<-array(0,dim=c(8,8,10))
assMatrix[,,5]<-amg2
assMatrix[,,10]<-1
ng2t5 <-nbdaData(label="ng2t5",assMatrix=assMatrix, asoc_ilv=asoc2, multi_ilv=asoc2, orderAcq=oag2t5, idname=c("Tope","Jam","Pod","Meg","Rosie","Dot","Charlie","Biscuit"),timeAcq=c(21,23,26,27,49,51,53,57),endTime=654)

#Repeat for group 3 (Tamar)
sex3<-c(1,0,1,0,0,0,1,1,0,0,1,0)
age3<-c(10,9,5,5,3,3,1,1,1,1,0.75,0.75)
sex3<-cbind(sex3)
age3<-cbind(age3)
asoc3<-c("sex3","age3")
amg3<-matrix(data=c(0.0, 0.46, 0.44, 0.41, 0.44, 0.41, 0.36, 0.48, 0.39, 0.38, 0.40, 0.41,	
			  0.46, 0.0, 0.49, 0.52, 0.56, 0.57, 0.43, 0.51, 0.46, 0.47, 0.34, 0.32, 
 			  0.44, 0.49, 0.0, 0.56, 0.47, 0.44, 0.49, 0.56, 0.51, 0.49, 0.36, 0.37, 
			  0.41, 0.52, 0.56, 0.0, 0.54, 0.53, 0.48, 0.59, 0.48, 0.52, 0.36, 0.38, 
			  0.44, 0.56, 0.47, 0.54, 0.0, 0.62, 0.47, 0.51, 0.44, 0.46, 0.34, 0.34, 
			  0.41, 0.57, 0.44, 0.53, 0.62, 0.0, 0.42, 0.46, 0.43, 0.41, 0.29, 0.28, 
			  0.36, 0.43, 0.49, 0.48, 0.47, 0.42, 0.0, 0.57, 0.56, 0.49, 0.36, 0.34, 
			  0.48, 0.51, 0.56, 0.59, 0.51, 0.46, 0.57, 0.0, 0.56, 0.57, 0.39, 0.42, 
			  0.39, 0.46, 0.51, 0.48, 0.44, 0.43, 0.56, 0.56, 0.0, 0.56, 0.32, 0.34,
   			  0.38, 0.47, 0.49, 0.52, 0.46, 0.41, 0.49, 0.57, 0.56, 0.0, 0.38, 0.41, 
			  0.40, 0.34, 0.36, 0.36, 0.34, 0.29, 0.36, 0.39, 0.32, 0.38, 0.0, 0.58, 
  			  0.41, 0.32, 0.37, 0.38, 0.34, 0.28, 0.34, 0.42, 0.34, 0.41, 0.58, 0.0), nrow=12)
oag3t1<-c(4,9,6,2,1)
oag3t2<-c(1,8,4,6,9,12,5,3,7)
oag3t3<-c(6,9,8,12,1,4,7,3,5,2,10,11)
oag3t4<-c(1,4,5,7,6,9,2,12,8,3,10,11)
oag3t5<-c(3,1,4,6,9,10,5,7,11,8,2,12)

dim(amg3)

assMatrix<-array(0,dim=c(12,12,10))
assMatrix[,,1]<-amg3
assMatrix[,,6]<-1
ng3t1 <-nbdaData(label="ng3t1",assMatrix=assMatrix, asoc_ilv=asoc3, multi_ilv=asoc3,orderAcq=oag3t1, idname=c("Leah","Feet","India","Chai","Cassia","Cameron","Hazel","Daisy","Harry","Dougie","Rani","Khan"),timeAcq=c(3,5,15,28,40),endTime=73)

assMatrix<-array(0,dim=c(12,12,10))
assMatrix[,,2]<-amg3
assMatrix[,,7]<-1
ng3t2 <-nbdaData(label="ng3t2",assMatrix=assMatrix, asoc_ilv=asoc3, multi_ilv=asoc3, orderAcq=oag3t2, idname=c("Leah","Feet","India","Chai","Cassia","Cameron","Hazel","Daisy","Harry","Dougie","Rani","Khan"),timeAcq=c(3,5,8,8,9,9,28,59,77),endTime=108)

assMatrix<-array(0,dim=c(12,12,10))
assMatrix[,,3]<-amg3
assMatrix[,,8]<-1
ng3t3 <-nbdaData(label="ng3t3",assMatrix=assMatrix, asoc_ilv=asoc3, multi_ilv=asoc3, orderAcq=oag3t3, idname=c("Leah","Feet","India","Chai","Cassia","Cameron","Hazel","Daisy","Harry","Dougie","Rani","Khan"),timeAcq=c(2,18,19,19,21,21,21,22,30,53,55,138),endTime=270)

assMatrix<-array(0,dim=c(12,12,10))
assMatrix[,,4]<-amg3
assMatrix[,,9]<-1
ng3t4 <-nbdaData(label="ng3t4",assMatrix=assMatrix, asoc_ilv=asoc3, multi_ilv=asoc3, orderAcq=oag3t4, idname=c("Leah","Feet","India","Chai","Cassia","Cameron","Hazel","Daisy","Harry","Dougie","Rani","Khan"),timeAcq=c(4,4,4,4,4,8,10,11,12,17,57,85),endTime=193)

assMatrix<-array(0,dim=c(12,12,10))
assMatrix[,,5]<-amg3
assMatrix[,,10]<-1
ng3t5 <-nbdaData(label="ng3t5",assMatrix=assMatrix, asoc_ilv=asoc3, multi_ilv=asoc3, orderAcq=oag3t5, idname=c("Leah","Feet","India","Chai","Cassia","Cameron","Hazel","Daisy","Harry","Dougie","Rani","Khan"),timeAcq=c(5,8,8,8,9,10,11,14,14,16,21,148),endTime=263)

#Create model vectors
constraintsVectMatrix<-rbind(
##Social network models
#Task 1 
##models fit to group networks where social transmission is constrained to be the same across tasks
c(0,0,0,0,0,1,1,1,1,1,0,0,0,0), #No ILVs
c(0,0,0,0,0,1,1,1,1,1,0,0,2,0), #Additive - Sex
c(0,0,0,0,0,1,1,1,1,1,0,0,0,2), #Additive - Age
c(0,0,0,0,0,1,1,1,1,1,0,0,2,3), #Additive - Sex and Age
c(0,0,0,0,0,1,1,1,1,1,2,0,0,0), #Multiplicative - Sex
c(0,0,0,0,0,1,1,1,1,1,0,2,0,0), #Multiplicative - Age
c(0,0,0,0,0,1,1,1,1,1,2,3,0,0), #Multiplicative - Sex and Age
##models fit to group networks where social transmission is constrained to be the different across tasks
c(0,0,0,0,0,1,2,3,4,5,0,0,0,0),
c(0,0,0,0,0,1,2,3,4,5,0,0,6,0),
c(0,0,0,0,0,1,2,3,4,5,0,0,0,6),
c(0,0,0,0,0,1,2,3,4,5,0,0,6,7),
c(0,0,0,0,0,1,2,3,4,5,6,0,0,0),
c(0,0,0,0,0,1,2,3,4,5,0,6,0,0),
c(0,0,0,0,0,1,2,3,4,5,6,7,0,0),
##models fit to social networks where social transmission is constrained to be the same across tasks
c(1,1,1,1,1,0,0,0,0,0,0,0,0,0),
c(1,1,1,1,1,0,0,0,0,0,0,0,2,0),
c(1,1,1,1,1,0,0,0,0,0,0,0,0,2),
c(1,1,1,1,1,0,0,0,0,0,0,0,2,3),
c(1,1,1,1,1,0,0,0,0,0,2,0,0,0),
c(1,1,1,1,1,0,0,0,0,0,0,2,0,0),
c(1,1,1,1,1,0,0,0,0,0,2,3,0,0),
####models fit to social networks where social transmission is constrained to be the different across tasks
c(1,2,3,4,5,0,0,0,0,0,0,0,0,0), 
c(1,2,3,4,5,0,0,0,0,0,0,0,6,0),
c(1,2,3,4,5,0,0,0,0,0,0,0,0,6),
c(1,2,3,4,5,0,0,0,0,0,0,0,6,7),
c(1,2,3,4,5,0,0,0,0,0,6,0,0,0),
c(1,2,3,4,5,0,0,0,0,0,0,6,0,0),
c(1,2,3,4,5,0,0,0,0,0,6,7,0,0),
##asocial models
c(0,0,0,0,0,0,0,0,0,0,0,0,0,0),
c(0,0,0,0,0,0,0,0,0,0,0,0,1,0),
c(0,0,0,0,0,0,0,0,0,0,0,0,0,1),
c(0,0,0,0,0,0,0,0,0,0,0,0,1,2))

#Fit constant and gamma baseline rates of asocial learning to models 
constraintsVectMatrix2<-rbind(constraintsVectMatrix,constraintsVectMatrix)
baselineVect<-rep(c("constant","gamma"),each=dim(constraintsVectMatrix)[1])

#Create AIC table
aicTable_tada<-tadaAICtable(nbdadata=c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5","ng2t1", "ng2t2", "ng2t3", "ng2t4", "ng2t5","ng3t1", "ng3t2", "ng3t3", "ng3t4", "ng3t5"),
                            constraintsVectMatrix= constraintsVectMatrix2,baselineVect=baselineVect, iterations=1000)

print(aicTable_tada)

#Get support for each type of model, with each baseline rate, fit to each network type
typeByNetworksSupport(aicTable_tada)

write.csv(aicTable_tada@printTable,file="aicTable_tada.csv")

#Can also get support for the different network combinations across all model types and baselines
networksSupport(aicTable_tada)
#                         support numberOfModels
#0:0:0:0:0:0:0:0:0:0 1.540617e-13              8  asocial learning
#0:0:0:0:0:1:1:1:1:1 1.232553e-04             14  social learning through the group network where rate of social transmission is the same across all task types
#0:0:0:0:0:1:2:3:4:5 7.513868e-02             14  social learning through the group network where rate of social transmission is the different across all task types
#1:1:1:1:1:0:0:0:0:0 7.381589e-04             14  social learning through the social network where rate of social transmission is the same across all task types
#1:2:3:4:5:0:0:0:0:0 9.239999e-01             14  social learning through the social network where rate of social transmission is the different across all task types

#Because all the variables have been specified in the data objects, new constrained data objects with only the variables
#in the best model need to be created. 
#The easiest way to do this is to look at which model is top (56) and then use row 9 of the constraintsVectMatrix to constrain the data objects:
ng1t1_bestModel<-constrainedNBDAdata(ng1t1,constraintsVect =constraintsVectMatrix2[56,])
ng1t2_bestModel<-constrainedNBDAdata(ng1t2,constraintsVect =constraintsVectMatrix2[56,])
ng1t3_bestModel<-constrainedNBDAdata(ng1t3,constraintsVect =constraintsVectMatrix2[56,])
ng1t4_bestModel<-constrainedNBDAdata(ng1t4,constraintsVect =constraintsVectMatrix2[56,])
ng1t5_bestModel<-constrainedNBDAdata(ng1t5,constraintsVect =constraintsVectMatrix2[56,])
ng2t1_bestModel<-constrainedNBDAdata(ng2t1,constraintsVect =constraintsVectMatrix2[56,])
ng2t2_bestModel<-constrainedNBDAdata(ng2t2,constraintsVect =constraintsVectMatrix2[56,])
ng2t3_bestModel<-constrainedNBDAdata(ng2t3,constraintsVect =constraintsVectMatrix2[56,])
ng2t4_bestModel<-constrainedNBDAdata(ng2t4,constraintsVect =constraintsVectMatrix2[56,])
ng2t5_bestModel<-constrainedNBDAdata(ng2t5,constraintsVect =constraintsVectMatrix2[56,])
ng3t1_bestModel<-constrainedNBDAdata(ng3t1,constraintsVect =constraintsVectMatrix2[56,])
ng3t2_bestModel<-constrainedNBDAdata(ng3t2,constraintsVect =constraintsVectMatrix2[56,])
ng3t3_bestModel<-constrainedNBDAdata(ng3t3,constraintsVect =constraintsVectMatrix2[56,])
ng3t4_bestModel<-constrainedNBDAdata(ng3t4,constraintsVect =constraintsVectMatrix2[56,])
ng3t5_bestModel<-constrainedNBDAdata(ng3t5,constraintsVect =constraintsVectMatrix2[56,])

bestModelNew<-tadaFit(list(ng1t1_bestModel,ng1t2_bestModel,ng1t3_bestModel,ng1t4_bestModel,ng1t5_bestModel,
				   ng2t1_bestModel,ng2t2_bestModel,ng2t3_bestModel,ng2t4_bestModel,ng2t5_bestModel,
				   ng3t1_bestModel,ng3t2_bestModel,ng3t3_bestModel,ng3t4_bestModel,ng3t5_bestModel),
			    baseline="gamma")
bestModelNew
bestModelNew@aicc
print(aicTable_tada)[1,]

#A neat way to extract the fitted parameters is:
data.frame(
Variable=bestModelNew@varNames,
MLE=bestModelNew@outputPar)

#                 Variable          MLE
#1         Scale (1/rate): 2072.6557292
#2                   Shape    0.5996747
#3 1 Social transmission 1    1.0881123 #Social transmission rate in Task 1
#4 2 Social transmission 2    2.6084788 #Social transmission rate in Task 2
#5 3 Social transmission 3    3.2617435 #Social transmission rate in Task 3
#6 4 Social transmission 4    6.0562472 #Social transmission rate in Task 4
#7 5 Social transmission 5    8.3298699 #Social transmission rate in Task 5
#8 6 Social= asocial: age1    0.1160667 #Effect of age on learning rate


#The second number on the left indicates which parameter numbers to ask for below when ascertaining 95% confidence intervals for s parameters.
#So paramter 1 is s1 and so on:

#Plot profile likelihoods for s1 to find range for lower and upper CIs  
plotProfLik(which=1,model=bestModelNew,range=c(0,20),resolution=10)
plotProfLik(which=1,model=bestModelNew,range=c(0,10),resolution=30)
#Ascertain lower and upper 95% CIs
profLikCI(which = 1,model=bestModelNew,upperRange = c(4,6)) #In this case only had to constrain upper range as lower CI cannot be below zero for s parameters
#Lower CI   Upper CI 
#[1] 0.00 5.133396   

#Reapeat for parameters 2-6
plotProfLik(which=2,model=bestModelNew,range=c(0,60),resolution=10)
plotProfLik(which=2,model=bestModelNew,range=c(0,10),resolution=30)
profLikCI(which = 2,model=bestModelNew,upperRange = c(7,9),lowerRange=c(0,2))
#Lower CI   Upper CI 
#[2]  0.5450606 7.7512870   

plotProfLik(which=3,model=bestModelNew,range=c(0,50),resolution=10)
plotProfLik(which=3,model=bestModelNew,range=c(0,20),resolution=30)
profLikCI(which = 3,model=bestModelNew,upperRange = c(7,10),lowerRange=c(0,3))
# Lower CI   Upper CI 
#[3] 0.8008696 9.4143231 

plotProfLik(which=4,model=bestModelNew,range=c(0,50),resolution=10)
plotProfLik(which=4,model=bestModelNew,range=c(0,20),resolution=30)
profLikCI(which = 4,model=bestModelNew,upperRange = c(15,18),lowerRange = c(0,4))
# Lower CI   Upper CI 
#[4] 1.992449 15.520164  

plotProfLik(which=5,model=bestModelNew,range=c(0,60),resolution=10)
plotProfLik(which=5,model=bestModelNew,range=c(0,10),resolution=30)
plotProfLik(which=5,model=bestModelNew,range=c(20,23),resolution=30)
profLikCI(which = 5,model=bestModelNew,upperRange = c(31,33),lowerRange = c(2,4))
#Lower CI Upper CI 
#[5] 2.607495 31.000082

exp(0.1160667)
#[Age] 1.123071
plotProfLik(which=6,model=bestModelNew,range=c(0,0.5),resolution=50)
profLikCI(which = 6,model=bestModelNew,upperRange = c(0.1,0.2),lowerRange = c(0,0.1))
#Lower CI         Upper CI 
#[6] 0.06172406   0.16503709
exp(0.06172406)
exp(0.16503709)
#[6] 1.063669     1.179437

#Get the proportion of events estimated to due to social transmission
nbdaPropSolveByST.byevent(model=bestModelNew)
nbdaPropSolveByST(model=bestModelNew)

#As it stands these figures are the proportion of all events explained by each network
#e.g. the proportion of acquisition events across all tasks explained by the network for task 1.
#Instead we want the proportion of task 1 acquisition events explained by the task 1 network etc.
#We therefore make the following adjustment:

#Divide P(Network 1) by the number of acquisition events for task 1 (excluding the innovators - the first otters in each group to interact with task 1 apparatus)/
#the total number of aquisition events across all task types (excluding the innovators)
0.06856/(12/96)
#[1] 0.54848

#Repeat for tasks 2 to 5
0.14739/(18/96)
#[2] 0.78608
        
0.19091/(22/96)
#[3] 0.8330618

0.20661/(22/96)
#[4] 0.9015709

0.21199/(22/96)
#[5] 0.9250473

#Ascertain lower 95% confidence interval for the number of acquisition events that were due to social learning for task 1
ng1t1_lower<-constrainedNBDAdata(ng1t1_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng1t2_lower<-constrainedNBDAdata(ng1t2_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng1t3_lower<-constrainedNBDAdata(ng1t3_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng1t4_lower<-constrainedNBDAdata(ng1t4_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng1t5_lower<-constrainedNBDAdata(ng1t5_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng2t1_lower<-constrainedNBDAdata(ng2t1_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng2t2_lower<-constrainedNBDAdata(ng2t2_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng2t3_lower<-constrainedNBDAdata(ng2t3_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng2t4_lower<-constrainedNBDAdata(ng2t4_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng2t5_lower<-constrainedNBDAdata(ng2t5_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng3t1_lower<-constrainedNBDAdata(ng3t1_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng3t2_lower<-constrainedNBDAdata(ng3t2_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng3t3_lower<-constrainedNBDAdata(ng3t3_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng3t4_lower<-constrainedNBDAdata(ng3t4_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
ng3t5_lower<-constrainedNBDAdata(ng3t5_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(0,0,0,0,0,0))
tadabestGamma_lower<-tadaFit(nbdadata=c("ng1t1_lower", "ng1t2_lower", "ng1t3_lower", "ng1t4_lower", "ng1t5_lower",
						    "ng2t1_lower", "ng1t2_lower", "ng2t3_lower", "ng2t4_lower", "ng2t5_lower",
 						    "ng3t1_lower", "ng3t2_lower", "ng3t3_lower", "ng3t4_lower", "ng3t5_lower"),baseline="gamma",iterations=10000)
tadabestGamma_lower@outputPar
nbdaPropSolveByST.byevent(model=tadabestGamma_lower)
nbdaPropSolveByST(model=tadabestGamma_lower)
0.00000/(12/96)
#[1] 0

#Ascertain higher 95% confidence interval for the number of acquisition events that were due to social learning for task 1
ng1t1_higher<-constrainedNBDAdata(ng1t1_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng1t2_higher<-constrainedNBDAdata(ng1t2_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng1t3_higher<-constrainedNBDAdata(ng1t3_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng1t4_higher<-constrainedNBDAdata(ng1t4_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng1t5_higher<-constrainedNBDAdata(ng1t5_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng2t1_higher<-constrainedNBDAdata(ng2t1_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng2t2_higher<-constrainedNBDAdata(ng2t2_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng2t3_higher<-constrainedNBDAdata(ng2t3_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng2t4_higher<-constrainedNBDAdata(ng2t4_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng2t5_higher<-constrainedNBDAdata(ng2t5_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng3t1_higher<-constrainedNBDAdata(ng3t1_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng3t2_higher<-constrainedNBDAdata(ng3t2_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng3t3_higher<-constrainedNBDAdata(ng3t3_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng3t4_higher<-constrainedNBDAdata(ng3t4_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
ng3t5_higher<-constrainedNBDAdata(ng3t5_bestModel,constraintsVect=c(0,1,2,3,4,5),offsetVect = c(5.133396,0,0,0,0,0))
tadabestGamma_higher<-tadaFit(nbdadata=c("ng1t1_higher", "ng1t2_higher", "ng1t3_higher", "ng1t4_higher", "ng1t5_higher",
						     "ng2t1_higher", "ng2t2_higher", "ng2t3_higher", "ng2t4_higher", "ng2t5_higher",
						     "ng3t1_higher", "ng3t2_higher", "ng3t3_higher", "ng3t4_higher", "ng3t5_higher"),baseline="gamma",iterations=10000)
tadabestGamma_higher@outputPar
nbdaPropSolveByST.byevent(model=tadabestGamma_higher)
nbdaPropSolveByST(model=tadabestGamma_higher)
0.10528 /(12/96)
#[1] 0.84224

#Repeat for tasks 2 to 5
ng1t1_lower2<-constrainedNBDAdata(ng1t1_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng1t2_lower2<-constrainedNBDAdata(ng1t2_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng1t3_lower2<-constrainedNBDAdata(ng1t3_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng1t4_lower2<-constrainedNBDAdata(ng1t4_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng1t5_lower2<-constrainedNBDAdata(ng1t5_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng2t1_lower2<-constrainedNBDAdata(ng2t1_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng2t2_lower2<-constrainedNBDAdata(ng2t2_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng2t3_lower2<-constrainedNBDAdata(ng2t3_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng2t4_lower2<-constrainedNBDAdata(ng2t4_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng2t5_lower2<-constrainedNBDAdata(ng2t5_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng3t1_lower2<-constrainedNBDAdata(ng3t1_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng3t2_lower2<-constrainedNBDAdata(ng3t2_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng3t3_lower2<-constrainedNBDAdata(ng3t3_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng3t4_lower2<-constrainedNBDAdata(ng3t4_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
ng3t5_lower2<-constrainedNBDAdata(ng3t5_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,0.5450606,0,0,0,0))
tadabestGamma_lower2<-tadaFit(nbdadata=c("ng1t1_lower2", "ng1t2_lower2", "ng1t3_lower2", "ng1t4_lower2", "ng1t5_lower2",
						     "ng2t1_lower2", "ng1t2_lower2", "ng2t3_lower2", "ng2t4_lower2", "ng2t5_lower2",
 						     "ng3t1_lower2", "ng3t2_lower2", "ng3t3_lower2", "ng3t4_lower2", "ng3t5_lower2"),baseline="gamma",iterations=10000)
tadabestGamma_lower2@outputPar
nbdaPropSolveByST.byevent(model=tadabestGamma_lower2)
nbdaPropSolveByST(model=tadabestGamma_lower2)
0.07406/(18/96)
#[2] 0.3949867
ng1t1_higher2<-constrainedNBDAdata(ng1t1_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng1t2_higher2<-constrainedNBDAdata(ng1t2_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng1t3_higher2<-constrainedNBDAdata(ng1t3_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng1t4_higher2<-constrainedNBDAdata(ng1t4_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng1t5_higher2<-constrainedNBDAdata(ng1t5_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng2t1_higher2<-constrainedNBDAdata(ng2t1_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng2t2_higher2<-constrainedNBDAdata(ng2t2_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng2t3_higher2<-constrainedNBDAdata(ng2t3_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng2t4_higher2<-constrainedNBDAdata(ng2t4_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng2t5_higher2<-constrainedNBDAdata(ng2t5_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng3t1_higher2<-constrainedNBDAdata(ng3t1_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng3t2_higher2<-constrainedNBDAdata(ng3t2_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng3t3_higher2<-constrainedNBDAdata(ng3t3_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng3t4_higher2<-constrainedNBDAdata(ng3t4_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
ng3t5_higher2<-constrainedNBDAdata(ng3t5_bestModel,constraintsVect=c(1,0,2,3,4,5),offsetVect = c(0,7.7512870,0,0,0,0))
tadabestGamma_higher2<-tadaFit(nbdadata=c("ng1t1_higher2", "ng1t2_higher2", "ng1t3_higher2", "ng1t4_higher2", "ng1t5_higher2",
						      "ng2t1_higher2", "ng1t2_higher2", "ng2t3_higher2", "ng2t4_higher2", "ng2t5_higher2",
 						      "ng3t1_higher2", "ng3t2_higher2", "ng3t3_higher2", "ng3t4_higher2", "ng3t5_higher2"),baseline="gamma",iterations=10000)
tadabestGamma_higher2@outputPar
nbdaPropSolveByST.byevent(model=tadabestGamma_higher2)
nbdaPropSolveByST(model=tadabestGamma_higher2)
0.15381/(18/96) 
#[2] 0.96555

ng1t1_lower3<-constrainedNBDAdata(ng1t1_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng1t2_lower3<-constrainedNBDAdata(ng1t2_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng1t3_lower3<-constrainedNBDAdata(ng1t3_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng1t4_lower3<-constrainedNBDAdata(ng1t4_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng1t5_lower3<-constrainedNBDAdata(ng1t5_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng2t1_lower3<-constrainedNBDAdata(ng2t1_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng2t2_lower3<-constrainedNBDAdata(ng2t2_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng2t3_lower3<-constrainedNBDAdata(ng2t3_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng2t4_lower3<-constrainedNBDAdata(ng2t4_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng2t5_lower3<-constrainedNBDAdata(ng2t5_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng3t1_lower3<-constrainedNBDAdata(ng3t1_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng3t2_lower3<-constrainedNBDAdata(ng3t2_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng3t3_lower3<-constrainedNBDAdata(ng3t3_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng3t4_lower3<-constrainedNBDAdata(ng3t4_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
ng3t5_lower3<-constrainedNBDAdata(ng3t5_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,0.8008696,0,0,0))
tadabestGamma_lower3<-tadaFit(nbdadata=c("ng1t1_lower3", "ng1t2_lower3", "ng1t3_lower3", "ng1t4_lower3", "ng1t5_lower3",
						     "ng2t1_lower3", "ng1t2_lower3", "ng2t3_lower3", "ng2t4_lower3", "ng2t5_lower3",
 						     "ng3t1_lower3", "ng3t2_lower3", "ng3t3_lower3", "ng3t4_lower3", "ng3t5_lower3"),baseline="gamma",iterations=10000)
tadabestGamma_lower3@outputPar
nbdaPropSolveByST.byevent(model=tadabestGamma_lower3)
nbdaPropSolveByST(model=tadabestGamma_lower3)
0.13704/(22/96)
#[3] 0.5979927
ng1t1_higher3<-constrainedNBDAdata(ng1t1_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng1t2_higher3<-constrainedNBDAdata(ng1t2_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng1t3_higher3<-constrainedNBDAdata(ng1t3_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng1t4_higher3<-constrainedNBDAdata(ng1t4_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng1t5_higher3<-constrainedNBDAdata(ng1t5_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng2t1_higher3<-constrainedNBDAdata(ng2t1_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng2t2_higher3<-constrainedNBDAdata(ng2t2_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng2t3_higher3<-constrainedNBDAdata(ng2t3_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng2t4_higher3<-constrainedNBDAdata(ng2t4_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng2t5_higher3<-constrainedNBDAdata(ng2t5_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng3t1_higher3<-constrainedNBDAdata(ng3t1_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng3t2_higher3<-constrainedNBDAdata(ng3t2_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng3t3_higher3<-constrainedNBDAdata(ng3t3_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng3t4_higher3<-constrainedNBDAdata(ng3t4_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
ng3t5_higher3<-constrainedNBDAdata(ng3t5_bestModel,constraintsVect=c(1,2,0,3,4,5),offsetVect = c(0,0,9.4143231,0,0,0))
tadabestGamma_higher3<-tadaFit(nbdadata=c("ng1t1_higher3", "ng1t2_higher3", "ng1t3_higher3", "ng1t4_higher3", "ng1t5_higher3",
						      "ng2t1_higher3", "ng1t2_higher3", "ng2t3_higher3", "ng2t4_higher3", "ng2t5_higher3",
 						      "ng3t1_higher3", "ng3t2_higher3", "ng3t3_higher3", "ng3t4_higher3", "ng3t5_higher3"),baseline="gamma",iterations=10000)
tadabestGamma_higher3@outputPar
nbdaPropSolveByST.byevent(model=tadabestGamma_higher3)
nbdaPropSolveByST(model=tadabestGamma_higher3)
0.21792/(22/96)
#[3] 0.9509236

ng1t1_lower4<-constrainedNBDAdata(ng1t1_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng1t2_lower4<-constrainedNBDAdata(ng1t2_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng1t3_lower4<-constrainedNBDAdata(ng1t3_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng1t4_lower4<-constrainedNBDAdata(ng1t4_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng1t5_lower4<-constrainedNBDAdata(ng1t5_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng2t1_lower4<-constrainedNBDAdata(ng2t1_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng2t2_lower4<-constrainedNBDAdata(ng2t2_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng2t3_lower4<-constrainedNBDAdata(ng2t3_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng2t4_lower4<-constrainedNBDAdata(ng2t4_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng2t5_lower4<-constrainedNBDAdata(ng2t5_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng3t1_lower4<-constrainedNBDAdata(ng3t1_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng3t2_lower4<-constrainedNBDAdata(ng3t2_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng3t3_lower4<-constrainedNBDAdata(ng3t3_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng3t4_lower4<-constrainedNBDAdata(ng3t4_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
ng3t5_lower4<-constrainedNBDAdata(ng3t5_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,1.992449,0,0))
tadabestGamma_lower4<-tadaFit(nbdadata=c("ng1t1_lower4", "ng1t2_lower4", "ng1t3_lower4", "ng1t4_lower4", "ng1t5_lower4",
						     "ng2t1_lower4", "ng1t2_lower4", "ng2t3_lower4", "ng2t4_lower4", "ng2t5_lower4",
 						     "ng3t1_lower4", "ng3t2_lower4", "ng3t3_lower4", "ng3t4_lower4", "ng3t5_lower4"),baseline="gamma",iterations=10000)
tadabestGamma_lower4@outputPar
nbdaPropSolveByST.byevent(model=tadabestGamma_lower4)
nbdaPropSolveByST(model=tadabestGamma_lower4)
0.17915/(22/96)
#[4] 0.7817455

ng1t1_higher4<-constrainedNBDAdata(ng1t1_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng1t2_higher4<-constrainedNBDAdata(ng1t2_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng1t3_higher4<-constrainedNBDAdata(ng1t3_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng1t4_higher4<-constrainedNBDAdata(ng1t4_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng1t5_higher4<-constrainedNBDAdata(ng1t5_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng2t1_higher4<-constrainedNBDAdata(ng2t1_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng2t2_higher4<-constrainedNBDAdata(ng2t2_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng2t3_higher4<-constrainedNBDAdata(ng2t3_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng2t4_higher4<-constrainedNBDAdata(ng2t4_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng2t5_higher4<-constrainedNBDAdata(ng2t5_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng3t1_higher4<-constrainedNBDAdata(ng3t1_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng3t2_higher4<-constrainedNBDAdata(ng3t2_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng3t3_higher4<-constrainedNBDAdata(ng3t3_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng3t4_higher4<-constrainedNBDAdata(ng3t4_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
ng3t5_higher4<-constrainedNBDAdata(ng3t5_bestModel,constraintsVect=c(1,2,3,0,4,5),offsetVect = c(0,0,0,15.520164,0,0))
tadabestGamma_higher4<-tadaFit(nbdadata=c("ng1t1_higher4", "ng1t2_higher4", "ng1t3_higher4", "ng1t4_higher4", "ng1t5_higher4",
						      "ng2t1_higher4", "ng1t2_higher4", "ng2t3_higher4", "ng2t4_higher4", "ng2t5_higher4",
 						      "ng3t1_higher4", "ng3t2_higher4", "ng3t3_higher4", "ng3t4_higher4", "ng3t5_higher4"),baseline="gamma",iterations=10000)
tadabestGamma_higher4@outputPar
nbdaPropSolveByST.byevent(model=tadabestGamma_higher4)
nbdaPropSolveByST(model=tadabestGamma_higher4)
0.22415/(22/96)
#[4] 0.9781091

ng1t1_lower5<-constrainedNBDAdata(ng1t1_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng1t2_lower5<-constrainedNBDAdata(ng1t2_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng1t3_lower5<-constrainedNBDAdata(ng1t3_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng1t4_lower5<-constrainedNBDAdata(ng1t4_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng1t5_lower5<-constrainedNBDAdata(ng1t5_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng2t1_lower5<-constrainedNBDAdata(ng2t1_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng2t2_lower5<-constrainedNBDAdata(ng2t2_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng2t3_lower5<-constrainedNBDAdata(ng2t3_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng2t4_lower5<-constrainedNBDAdata(ng2t4_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng2t5_lower5<-constrainedNBDAdata(ng2t5_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng3t1_lower5<-constrainedNBDAdata(ng3t1_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng3t2_lower5<-constrainedNBDAdata(ng3t2_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng3t3_lower5<-constrainedNBDAdata(ng3t3_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng3t4_lower5<-constrainedNBDAdata(ng3t4_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
ng3t5_lower5<-constrainedNBDAdata(ng3t5_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,2.607495,0))
tadabestGamma_lower5<-tadaFit(nbdadata=c("ng1t1_lower5", "ng1t2_lower5", "ng1t3_lower5", "ng1t4_lower5", "ng1t5_lower5",
						     "ng2t1_lower5", "ng1t2_lower5", "ng2t3_lower5", "ng2t4_lower5", "ng2t5_lower5",
 						     "ng3t1_lower5", "ng3t2_lower5", "ng3t3_lower5", "ng3t4_lower5", "ng3t5_lower5"),baseline="gamma",iterations=10000)
tadabestGamma_lower5@outputPar
nbdaPropSolveByST.byevent(model=tadabestGamma_lower5)
nbdaPropSolveByST(model=tadabestGamma_lower5)
0.18854/(22/96)
#[5] 0.82272

ng1t1_higher5<-constrainedNBDAdata(ng1t1_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng1t2_higher5<-constrainedNBDAdata(ng1t2_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng1t3_higher5<-constrainedNBDAdata(ng1t3_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng1t4_higher5<-constrainedNBDAdata(ng1t4_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng1t5_higher5<-constrainedNBDAdata(ng1t5_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng2t1_higher5<-constrainedNBDAdata(ng2t1_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng2t2_higher5<-constrainedNBDAdata(ng2t2_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng2t3_higher5<-constrainedNBDAdata(ng2t3_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng2t4_higher5<-constrainedNBDAdata(ng2t4_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng2t5_higher5<-constrainedNBDAdata(ng2t5_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng3t1_higher5<-constrainedNBDAdata(ng3t1_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng3t2_higher5<-constrainedNBDAdata(ng3t2_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng3t3_higher5<-constrainedNBDAdata(ng3t3_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng3t4_higher5<-constrainedNBDAdata(ng3t4_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
ng3t5_higher5<-constrainedNBDAdata(ng3t5_bestModel,constraintsVect=c(1,2,3,4,0,5),offsetVect = c(0,0,0,0,31.000082,0))
tadabestGamma_higher5<-tadaFit(nbdadata=c("ng1t1_higher5", "ng1t2_higher5", "ng1t3_higher5", "ng1t4_higher5", "ng1t5_higher5",
						      "ng2t1_higher5", "ng1t2_higher5", "ng2t3_higher5", "ng2t4_higher5", "ng2t5_higher5",
 						      "ng3t1_higher5", "ng3t2_higher5", "ng3t3_higher5", "ng3t4_higher5", "ng3t5_higher5"),baseline="gamma",iterations=10000)
tadabestGamma_higher5@outputPar
nbdaPropSolveByST.byevent(model=tadabestGamma_higher5)
nbdaPropSolveByST(model=tadabestGamma_higher5)
0.22889/(22/96)
#[5] 0.9987927 



