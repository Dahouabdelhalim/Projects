#Load the NBDA package as follows
library(NBDA)

#Set working directory to where the data files are located on your computer



#############################################################################
# TUTORIAL 6.1
# MODELLING MULTIPLE DIFFUSIONS USING OADA
# 4 diffusions
# 1 static network
# 1 time constant ILV
#############################################################################

#Read in the social network and order of acquisition vector as shown in Tutorial 1.
socNet1<-as.matrix(read.csv(file="multiDiffSocNet1.csv"))
socNet1<-array(socNet1,dim=c(30,30,1))
oa1<-c(16,9,20,2,14,28,8,15,19,30,3,25,22,18,24,11,29,27,6,4,26,21,10,7,5,13,1,17,12,23)
ILVdata<-read.csv(file="multiDiffILV1.csv")

#Then extract the ILV- age
age<-cbind(ILVdata$age)
#age in years
#Each ILV needs to be stored as a column matrix; this is to allow extension to time-varying ILVs (see below)

#I am going to standardize age
stAge<-(age-mean(age))/sd(age)
#Since age is centred on 0, and males=0, females=1, the baseline asocial rate is a male of mean age
#Therefore s is estimated relative to this baseline

asoc1<-c("stAge")
#We create a vector of the names of the variables to be included in the analysis

#Now we create an nbdaData object as before, but specifying the ILVs to be included

#Or the "unconstrained" model
nbdaDataMulti1<-nbdaData(label="MultiDiffusion1",assMatrix=socNet1,orderAcq = oa1,asoc_ilv = asoc1,int_ilv=asoc1)


#Now do the same to create separate objects for the other 3 diffusions
socNet2<-as.matrix(read.csv(file="multiDiffSocNet2.csv"))
socNet2<-array(socNet2,dim=c(30,30,1))
oa2<-c(26,29,30,8,19,21,22,3,14,12,11,1,17,28,5,9,15,7,6,25,4,13,20,23,2,16,27,18,24,10)
ILVdata2<-read.csv(file="multiDiffILV2.csv")
age2<-cbind(ILVdata2$age)
stAge2<-(age-mean(age2))/sd(age2)
asoc2<-c("stAge2")
nbdaDataMulti2<-nbdaData(label="MultiDiffusion2",assMatrix=socNet2,orderAcq = oa2,asoc_ilv = asoc2,int_ilv=asoc2)

socNet3<-as.matrix(read.csv(file="multiDiffSocNet3.csv"))
socNet3<-array(socNet3,dim=c(30,30,1))
oa3<-c(20,3,2,1,12,11,10,16,26,30,24,5,29,17,13,21,28,4,18,22,19,15,7,14,27,8,9,6,25,23)
ILVdata3<-read.csv(file="multiDiffILV3.csv")
age3<-cbind(ILVdata3$age)
stAge3<-(age-mean(age3))/sd(age3)
asoc3<-c("stAge3")
nbdaDataMulti3<-nbdaData(label="MultiDiffusion3",assMatrix=socNet3,orderAcq = oa3,asoc_ilv = asoc3,int_ilv=asoc3)

socNet4<-as.matrix(read.csv(file="multiDiffSocNet4.csv"))
socNet4<-array(socNet4,dim=c(30,30,1))
oa4<-c(19,14,13,15,18,24,12,23,8,30,5,11,7,29,1,16,2,6,28,4,17,9,22,21,26,25,20,10,3,27)
ILVdata4<-read.csv(file="multiDiffILV4.csv")
age4<-cbind(ILVdata4$age)
stAge4<-(age-mean(age4))/sd(age4)
asoc4<-c("stAge4")
nbdaDataMulti4<-nbdaData(label="MultiDiffusion4",assMatrix=socNet4,orderAcq = oa4,asoc_ilv = asoc4,int_ilv=asoc4)

#Instead of specifying a single nbdaData object, we specify a list containing all of the diffusions we wish to include:
multiDiffModel1<-oadaFit(list(nbdaDataMulti1,nbdaDataMulti2,nbdaDataMulti3,nbdaDataMulti4))
data.frame(Variable=multiDiffModel1@varNames,MLE=multiDiffModel1@outputPar,SE=multiDiffModel1@se)
#              Variable         MLE        SE
#1 1 Social transmission 1  2.94289456 2.1577225
#2        2 Asocial: stAge -0.26368738 0.3249724
#3         3 Social: stAge  0.08629645 0.1244382

multiDiffModel1_asocial<-oadaFit(list(nbdaDataMulti1,nbdaDataMulti2,nbdaDataMulti3,nbdaDataMulti4),type="asocial")
multiDiffModel1@aicc
#[1] 581.0339
multiDiffModel1_asocial@aicc
#[1] 599.1968

#Further inference can then proceed as in previous tutorials

#############################################################################
# TUTORIAL 6.2
# MODELLING MULTIPLE DIFFUSIONS USING TADA
# 4 diffusions
# 1 static network
# 1 time constant ILV + group ILVs
#############################################################################

#We may prefer to fit a TADA model for the reasons explained in the main text. To do this, we need the times of acqusition,
#which we can enter as a vector for each diffusion
ta1<-c(2069,2080,2083,2095,2105,2112,2117,2118,2122,2129,2130,2131,2132,2136,2138,2140,2145,2150,2151,2153,2154,2158,2158,2161,2166,2168,2177,2180,2194,2236)
ta2<-c(964,971,985,988,990,998,999,1000,1000,1001,1001,1002,1002,1002,1004,1006,1006,1008,1011,1015,1017,1020,1022,1023,1025,1027,1029,1037,1039,1079)
ta3<-c(306,1008,1796,1848,2046,2223,2402,2621,3281,3935,4428,4548,5729,6251,6270,6764,7231,7327,8009,8010,8168,9037,9497,10237,10249,11354,13448,14275,16586,16797)
ta4<-c(2268,2495,2967,3016,7988,8157,9716,9856,9878,10773,10829,12550,13957,16205,18075,18262,19073,23988,25781,26089,26216,26944,27698,27985,28515,28968,30447,30451,31126,31264)

#We can then create new nbdaData objects with times included
nbdaDataMulti1t<-nbdaData(label="MultiDiffusion1",assMatrix=socNet1,orderAcq = oa1,asoc_ilv = asoc1,int_ilv=asoc1,timeAcq = ta1,endTime=max(ta1)+1)
nbdaDataMulti2t<-nbdaData(label="MultiDiffusion2",assMatrix=socNet2,orderAcq = oa2,asoc_ilv = asoc2,int_ilv=asoc2,timeAcq = ta2,endTime=max(ta2)+1)
nbdaDataMulti3t<-nbdaData(label="MultiDiffusion3",assMatrix=socNet3,orderAcq = oa3,asoc_ilv = asoc3,int_ilv=asoc3,timeAcq = ta3,endTime=max(ta3)+1)
nbdaDataMulti4t<-nbdaData(label="MultiDiffusion4",assMatrix=socNet4,orderAcq = oa4,asoc_ilv = asoc4,int_ilv=asoc4,timeAcq = ta4,endTime=max(ta4)+1)


multiDiffModel1t<-tadaFit(list(nbdaDataMulti1t,nbdaDataMulti2t,nbdaDataMulti3t,nbdaDataMulti4t))
data.frame(Variable=multiDiffModel1t@varNames,MLE=multiDiffModel1t@outputPar,SE=multiDiffModel1t@se)
#               Variable           MLE           SE
#1         Scale (1/rate): 32164.5098863 1.181200e+04
#2 1 Social transmission 1     1.6883211 7.580722e-01
#3        2 Asocial: stAge    -0.5137037 2.929822e-01
#4         3 Social: stAge     0.1777465 1.178132e-01

multiDiffModel1_asocialt<-tadaFit(list(nbdaDataMulti1t,nbdaDataMulti2t,nbdaDataMulti3t,nbdaDataMulti4t),type="asocial",gradient=F)
multiDiffModel1t@aicc
#[1] 2308.022
multiDiffModel1_asocialt@aicc
#[1] 2369.141

#A huge difference in AICc here. However, if there were group differences in asocial learning rate, this might cause a spurious positive result for
#social transmission. In some cases, we might be able to justify assuming that such differences are implausible.
#If not, we would want to allow for the possibility that different groups learned asocially at different rates.
#we need to create a set of 3 (4-1) indicator variables to make up a group factor:

#Indicator for group 2
group2_1<-cbind(rep(0,30))
group2_2<-cbind(rep(1,30))
group2_3<-cbind(rep(0,30))
group2_4<-cbind(rep(0,30))

#Indicator for group 3
group3_1<-cbind(rep(0,30))
group3_2<-cbind(rep(0,30))
group3_3<-cbind(rep(1,30))
group3_4<-cbind(rep(0,30))

#Indicator for group 4
group4_1<-cbind(rep(0,30))
group4_2<-cbind(rep(0,30))
group4_3<-cbind(rep(0,30))
group4_4<-cbind(rep(1,30))

asoc1_group<-c("stAge","group2_1","group3_1","group4_1")
asoc2_group<-c("stAge2","group2_2","group3_2","group4_2")
asoc3_group<-c("stAge3","group2_3","group3_3","group4_3")
asoc4_group<-c("stAge4","group2_4","group3_4","group4_4")

#We can then recreate the nbdaData objects with group variable added to asoc_ilv
#If we added the group variable to the int_ilv slot, we could test for differences among groups in social transmission rate
#But we will suggest an alternative way to do this below
nbdaDataMulti1t<-nbdaData(label="MultiDiffusion1",assMatrix=socNet1,orderAcq = oa1,asoc_ilv = asoc1_group,int_ilv=asoc1,timeAcq = ta1,endTime=max(ta1)+1)
nbdaDataMulti2t<-nbdaData(label="MultiDiffusion2",assMatrix=socNet2,orderAcq = oa2,asoc_ilv = asoc2_group,int_ilv=asoc2,timeAcq = ta2,endTime=max(ta2)+1)
nbdaDataMulti3t<-nbdaData(label="MultiDiffusion3",assMatrix=socNet3,orderAcq = oa3,asoc_ilv = asoc3_group,int_ilv=asoc3,timeAcq = ta3,endTime=max(ta3)+1)
nbdaDataMulti4t<-nbdaData(label="MultiDiffusion4",assMatrix=socNet4,orderAcq = oa4,asoc_ilv = asoc4_group,int_ilv=asoc4,timeAcq = ta4,endTime=max(ta4)+1)

multiDiffModel1tg<-tadaFit(list(nbdaDataMulti1t,nbdaDataMulti2t,nbdaDataMulti3t,nbdaDataMulti4t))
data.frame(Variable=multiDiffModel1tg@varNames,MLE=multiDiffModel1tg@outputPar,SE=multiDiffModel1tg@se)

#             Variable           MLE           SE
#1         Scale (1/rate): 2706.81197398 637.37535976
#2 1 Social transmission 1    0.07500999   0.02627594
#3       2 Asocial: stAge1   -0.04110105   0.14469662
#4     3 Asocial: group2_1    0.88489896   0.30669228
#5     4 Asocial: group3_1   -1.58149187   0.40669331
#6     5 Asocial: group4_1   -3.52206298   0.66550921
#7         6 Social: stAge    0.11403643   0.17346317

multiDiffModel1tg_asocial<-tadaFit(list(nbdaDataMulti1t,nbdaDataMulti2t,nbdaDataMulti3t,nbdaDataMulti4t),type="asocial",gradient=F)
multiDiffModel1tg@aicc
#[1] 2218.423
multiDiffModel1tg_asocial@aicc
#[1] 2243.665

#We can now see that the evidence for social transmission, and the esimate of s is reduced once we account for potential group differences
#in asocial learning rate
#It is neither necessary nor possible to include group differences in a standard OADA (see exercise 6.1 above) since groups are not assumed to have the
#same baseline rate function; this means that group differences are absorbed into the baseline rate function for each group

#We also need to allow for the possibility that social transmission is a spurious result of an increasing baseline function, using weibull and/or gamma
#baseline rate functions:

multiDiffModel1tg_weibull<-tadaFit(list(nbdaDataMulti1t,nbdaDataMulti2t,nbdaDataMulti3t,nbdaDataMulti4t),baseline="weibull")
data.frame(Variable=multiDiffModel1tg_weibull@varNames,MLE=multiDiffModel1tg_weibull@outputPar,SE=multiDiffModel1tg_weibull@se)
#                 Variable           MLE  SE
#1         Scale (1/rate):  8.176185e+13 NaN
#2                   Shape  1.302955e-01 NaN
#3 1 Social transmission 1  9.224943e+01 NaN
#4       2 Asocial: stAge1 -7.152392e-01 NaN
#5     3 Asocial: group2_1  4.645889e-01 NaN
#6     4 Asocial: group3_1  5.977118e-01 NaN
#7     5 Asocial: group4_1  8.548826e-02 NaN
#8         6 Social: stAge  2.100482e-01 NaN

#(no SEs can be derived in this case so NaNs are returned)
#In fact, the shape is <1 indicating a decreasing baseline rate, which was potentially obscuring social transmission in the constant baseline model

#fit asocial model for comparison
multiDiffModel1tg_weibull_asocial<-tadaFit(list(nbdaDataMulti1t,nbdaDataMulti2t,nbdaDataMulti3t,nbdaDataMulti4t),baseline="weibull",type="asocial",gradient=F)

multiDiffModel1tg_gamma<-tadaFit(list(nbdaDataMulti1t,nbdaDataMulti2t,nbdaDataMulti3t,nbdaDataMulti4t),baseline="gamma")
data.frame(Variable=multiDiffModel1tg_gamma@varNames,MLE=multiDiffModel1tg_gamma@outputPar,SE=multiDiffModel1tg_gamma@se)
#Variable           MLE           SE
#1         Scale (1/rate): 72.5449262908 1.162253e+01
#2                   Shape 28.5156358622 4.361002e+00
#3 1 Social transmission 1  0.0007663599 5.022673e-04
#4       2 Asocial: stAge1  0.0236731875 1.085577e-01
#5     3 Asocial: group2_1  7.7968752485 1.026750e+00
#6     4 Asocial: group3_1 -3.8346202673 2.666178e-01
#7     5 Asocial: group4_1 -5.6046110620 4.379389e-01
#8         6 Social: stAge -0.0369562901 3.828953e-01

#The gamma model reaches the opposite conclusion: that there is an increasing baseline rate (shape >1) and no social transmission

#fit asocial model for comparison
multiDiffModel1tg_gamma_asocial<-tadaFit(list(nbdaDataMulti1t,nbdaDataMulti2t,nbdaDataMulti3t,nbdaDataMulti4t),baseline="gamma",type="asocial")

#Compare all AICs
multiDiffModel1tg@aicc
#2218.423
multiDiffModel1tg_asocial@aicc
#2243.665
multiDiffModel1tg_weibull@aicc
#2169.994
multiDiffModel1tg_weibull_asocial@aicc
#2136.367
multiDiffModel1tg_gamma@aicc
#2041.317
multiDiffModel1tg_gamma_asocial@aicc
#2039.465

#This suggest that the gamma baseline function with asocial learning only is favoured
#However, given the wildly different results depending on our assumptions about the baseline function, it suggests the
#time course of events is dominating the analysis. In such circmstances, we would be inclined to switch to the OADA which makes inferences
#only on the order of spread.

#However, as noted in the main text, OADA does not take into account diffusions taking off at different times in different diffusions.
#A compromise between OADA and TADA is the stratified OADA, which we turn to in the next exercise:

#############################################################################
# TUTORIAL 6.3
# MODELLING MULTIPLE DIFFUSIONS USING STRATIFIED OADA
# 4 diffusions
# 1 static network
# 1 time constant ILV + group ILVs
#############################################################################

#This simply involves modelling the 4 diffiusions as a single diffusion with zero connections between those in different groups
#So the first stage is to build a big network covering all groups. The dimensions will be 120 x 120 x 1, but we can calculate this
#automatically as follows
c((dim(socNet1)+dim(socNet2)+dim(socNet3)+dim(socNet4))[1:2],1)
#And build an array of 0s of the appropriate size
bigSocNet<-array(0,dim=c((dim(socNet1)+dim(socNet2)+dim(socNet3)+dim(socNet4))[1:2],1))
dim(bigSocNet)
#[1] 120 120   1

#Now we can fill in the appropriate parts of the network as follows
bigSocNet[1:30,1:30,]<-socNet1
bigSocNet[31:60,31:60,]<-socNet2
bigSocNet[61:90,61:90,]<-socNet3
bigSocNet[91:120,91:120,]<-socNet4

#Now we need to get the order of acquisition across each group
#First we need to change the identities of individuals in group 2-4 so we have 31-60 in group 2:
oa2_new<-oa2+30
range(oa2_new)
#61-90 in group 3
oa3_new<-oa3+60
range(oa3_new)
#91-120 in group 4
oa4_new<-oa4+90
range(oa4_new)

#Now let us combine the oa vectors into one
oa_combined_unordered<-c(oa1,oa2_new,oa3_new,oa4_new)
#but this is not yet ordered across diffusions. To do this we need to combined the ta vectors into one:
ta_combined<-c(ta1,ta2,ta3,ta4)
#Now we can order oa_unordered by ta_combined
oa_combined<-oa_combined_unordered[order(ta_combined)]

#Now we need to combine the ILVs for the groups

stAgeAll<-rbind(stAge,stAge2,stAge3,stAge4)

#Note that in all these analyses we have standardized age within each group; depending on the variable and circumstance,
#it may make more sense to standardize ILVs across groups

group2All<-rbind(group2_1,group2_2,group2_3,group2_4)
group3All<-rbind(group3_1,group3_2,group3_3,group3_4)
group4All<-rbind(group4_1,group4_2,group4_3,group4_4)

asocAll<-c("stAgeAll")
asocAll_group<-c("stAgeAll","group2All","group3All","group4All")

#Now we can make a single nbdaData object for all diffusions:
nbdaDataMultiAll<-nbdaData(label="MultiDiffusionAll",assMatrix=bigSocNet,orderAcq = oa_combined,asoc_ilv = asocAll_group,int_ilv=asocAll)

#And fit the model the same way we would a standard OADA:
stratOADAModel<-oadaFit(nbdaDataMultiAll)
data.frame(Variable=stratOADAModel@varNames,MLE=stratOADAModel@outputPar,SE=stratOADAModel@se)

#            Variable         MLE          SE
#1 1 Social transmission 1  0.00180374 0.002851674
#2     2 Asocial: stAgeAll -0.09075625 0.108566126
#3    3 Asocial: group2All  4.01542653 0.763192945
#4    4 Asocial: group3All -2.77795403 0.504445792
#5    5 Asocial: group4All -4.52478958 0.604461563
#6      6 Social: stAgeAll -0.25408653 0.462389941

stratOADAModel_asocial<-oadaFit(nbdaDataMultiAll,type="asocial")
stratOADAModel@aicc
#[1] 718.6968
stratOADAModel_asocial@aicc
#[1] 715.4256

#In this case, we have no support for social transmission using a stratified OADA but very strong support in the normal OADA. We suspect this situation
#will be unusual for real data. If this result were obtained, it could be the result of a false assumption that the shape of the baseline rate
#is the same for all diffusions.
#In fact, the scenario these data were simulated from is a little different as we will see in the next exercise.


#############################################################################
# TUTORIAL 6.4
# TESTING FOR A DIFFERENCE IN S AMONG DIFFUSIONS
# 4 diffusions
# 1 static network
# 1 time constant ILV
#############################################################################

#What if we wish to test if the rate of social transmission, per unit connection, is different in different networks?
#For example, let us assume that diffusions 1 and 2 were conducted under one experimental condition (A), and diffusions 3 and 4 were conducted under another
#condition (B), and we wish to test if it affected the rate of transmission.
#The way to model this is to break the network down into two, one representing social transmission in condition A, and another in condition B.

#We will first see how to do this in a standard OADA, but it would be done in the same way for a TADA:

#For each group, we create a 2-network array of 30 x 30 x 2. The first network is for condition A and the second is for condition B.
#Therefore, network 2 if full of zeroes for diffusions 1 and 2, and network 2 is full of zeroes for diffusions 3 and 4.

newSocNet1<-array(0,dim=c(30,30,2))
newSocNet2<-array(0,dim=c(30,30,2))
newSocNet3<-array(0,dim=c(30,30,2))
newSocNet4<-array(0,dim=c(30,30,2))

newSocNet1[,,1]<-socNet1
newSocNet2[,,1]<-socNet2
newSocNet3[,,2]<-socNet3
newSocNet4[,,2]<-socNet4

#Now create new nbdaData objects

nbdaDataMulti1_2nets<-nbdaData(label="MultiDiffusion1_2nets",assMatrix=newSocNet1,orderAcq = oa1,asoc_ilv = asoc1,int_ilv=asoc1)
nbdaDataMulti2_2nets<-nbdaData(label="MultiDiffusion2_2nets",assMatrix=newSocNet2,orderAcq = oa2,asoc_ilv = asoc2,int_ilv=asoc2)
nbdaDataMulti3_2nets<-nbdaData(label="MultiDiffusion3_2nets",assMatrix=newSocNet3,orderAcq = oa3,asoc_ilv = asoc3,int_ilv=asoc3)
nbdaDataMulti4_2nets<-nbdaData(label="MultiDiffusion4_2nets",assMatrix=newSocNet4,orderAcq = oa4,asoc_ilv = asoc4,int_ilv=asoc4)

#and fit the model
multiDiffModel1_2nets<-oadaFit(list(nbdaDataMulti1_2nets,nbdaDataMulti2_2nets,nbdaDataMulti3_2nets,nbdaDataMulti4_2nets))
data.frame(Variable=multiDiffModel1_2nets@varNames,MLE=multiDiffModel1_2nets@outputPar,SE=multiDiffModel1_2nets@se)
#             Variable         MLE         SE
#1 1 Social transmission 1  9.49088144 11.7038043
#2 2 Social transmission 2  0.97543719  1.4395802
#3        3 Asocial: stAge -0.36288324  0.4242249
#4         4 Social: stAge  0.04118782  0.1364657

#Now we can see that there may be a big difference in s between conditions (let us call these parameters sA and sB)
#Let us test for a difference between sA and sB, by fitting a model in which they are constrained to be the same

#Create the constrained nbdaData objects:
nbdaDataMulti1_sAequalsB<-constrainedNBDAdata(nbdaDataMulti1_2nets,constraintsVect = c(1,1,2,3))
nbdaDataMulti2_sAequalsB<-constrainedNBDAdata(nbdaDataMulti2_2nets,constraintsVect = c(1,1,2,3))
nbdaDataMulti3_sAequalsB<-constrainedNBDAdata(nbdaDataMulti3_2nets,constraintsVect = c(1,1,2,3))
nbdaDataMulti4_sAequalsB<-constrainedNBDAdata(nbdaDataMulti4_2nets,constraintsVect = c(1,1,2,3))
#Fit the model
multiDiffModel1_sAequalsB<-oadaFit(list(nbdaDataMulti1_sAequalsB,nbdaDataMulti2_sAequalsB,nbdaDataMulti3_sAequalsB,nbdaDataMulti4_sAequalsB))
data.frame(Variable=multiDiffModel1_sAequalsB@varNames,MLE=multiDiffModel1_sAequalsB@outputPar,SE=multiDiffModel1_sAequalsB@se)
#              Variable         MLE        SE
#1 1 Social transmission 1  2.94289456 2.1577225
#2        2 Asocial: stAge -0.26368738 0.3249724
#3         3 Social: stAge  0.08629645 0.1244382

#Note here we have arrived at the same model we fitted in Exercise 6.1

#Compare AICcs
multiDiffModel1_2nets@aicc
#[1] 579.4531
multiDiffModel1_sAequalsB@aicc
#[1] 581.0339
exp(0.5*(multiDiffModel1_sAequalsB@aicc-multiDiffModel1_2nets@aicc))
#[1] 2.204332
#2.2x more support for a difference in rates between conditions

#LRT:
2*(multiDiffModel1_sAequalsB@loglik-multiDiffModel1_2nets@loglik)
#[1] 3.721778
pchisq(2*(multiDiffModel1_sAequalsB@loglik-multiDiffModel1_2nets@loglik),df=1,lower.tail = F)
#[1] 0.05370714
#Some fairly weak evidence of a difference.

#We can also get a 95% C.I. for sA-sB
plotProfLik(which=1,model=multiDiffModel1_2nets,range=c(0,500),constraintsVect = c(1,1,2,3))
#upper limit 200-300, lower limit is <0 so we will deal with this below
profLikCI(which=1,model=multiDiffModel1_2nets,upperRange =c(200,300),constraintsVect = c(1,1,2,3))
#Lower CI Upper CI
#0.0000 232.2918

#Repeat in the sB-sA direction:
plotProfLik(which=2,model=multiDiffModel1_2nets,range=c(0,0.5),constraintsVect = c(1,1,2,3))
#between 0.1 and 0.2
profLikCI(which=2,model=multiDiffModel1_2nets,upperRange =c(0.1,0.2),constraintsVect = c(1,1,2,3))
#Lower CI  Upper CI
#0.0000000 0.1380924

#So the confidence interval for sA-sB is -0.14 to 232.29

#Further inference can proceed as in previous tutorials

#Let us now look at how to do the same thing using a stratified OADA
#Here we need a single 120x120x2 array, with connections for diffusions 1 and 2 entered into the 1st network
#and connections for diffusions 3 and 4 entered into the 2nd network

bigSocNet_2nets<-array(0,dim=c(120,120,2))

bigSocNet_2nets[1:30,1:30,1]<-socNet1
bigSocNet_2nets[31:60,31:60,1]<-socNet2
bigSocNet_2nets[61:90,61:90,2]<-socNet3
bigSocNet_2nets[91:120,91:120,2]<-socNet4

#Now we can make a single nbdaData object for all diffusions:
nbdaDataMultiAll_2nets<-nbdaData(label="MultiDiffusionAll_2nets",assMatrix=bigSocNet_2nets,orderAcq = oa_combined,asoc_ilv = asocAll_group,int_ilv=asocAll)
#And fit the model the same way we would a standard OADA:
stratOADAModel_2nets<-oadaFit(nbdaDataMultiAll_2nets)
data.frame(Variable=stratOADAModel_2nets@varNames,MLE=stratOADAModel_2nets@outputPar,SE=stratOADAModel_2nets@se)
#             Variable          MLE          SE
#1 1 Social transmission 1 160.00493890 188.3402929
#2 2 Social transmission 2   0.18877431   0.3038777
#3     3 Asocial: stAgeAll   0.03687398   0.1761472
#4    4 Asocial: group2All   2.73828124   1.4705340
#5    5 Asocial: group3All   1.52388911   1.0879288
#6    6 Asocial: group4All  -0.29845588   1.1571600
#7      7 Social: stAgeAll  -0.09463353   0.1465292

#Fit the sA=sB constrained model
nbdaDataMulti1_sAequalsB<-constrainedNBDAdata(nbdaDataMultiAll_2nets,constraintsVect = c(1,1,2,3,4,5,6))
#Fit the model
stratOADAModel_sAequalsB<-oadaFit(nbdaDataMulti1_sAequalsB)

#Compare AICcs
stratOADAModel_2nets@aicc
#[1] 665.0402
stratOADAModel_sAequalsB@aicc
#[1] 719.2048
exp(0.5*(stratOADAModel_sAequalsB@aicc-stratOADAModel_2nets@aicc))
#[1] 577691928395
#Overwhelming support for a difference in rates between conditions

#LRT:
2*(stratOADAModel_sAequalsB@loglik-stratOADAModel_2nets@loglik)
#[1] 56.42125
pchisq(2*(stratOADAModel_sAequalsB@loglik-stratOADAModel_2nets@loglik),df=1,lower.tail = F)
#[1] 5.849477e-14
#Very strong evidence of a difference.

#We could also get a 95% C.I. for sA-sB, though it is best to reparameterize so that group 2 is the baseline as shown in
#Tutorial 2.

#Further inference can proceed as in previous tutorials



#############################################################################
# TUTORIAL 6.5
# INCLUDING A GROUP NETWORK IN A cTADA AND STRATIFIED OADA
# 4 diffusions
# 1 static network
#############################################################################


#In a multiple diffusion analysis using TADA or stratified OADA, comparing a network-based model of social
#learning to an asocial model does not test whether the diffusion follows the network within each group
#(see main text for details)

#Here we show how to include a group network in the analysis
#For this tutorial, we will analyse two simulated datasets using the same social networks as tutorials 6.2-6.4 above.
#In the first dataset (Dataset A) we have social transmission following the social network
#In the second datset (Dataset B) we have social transmission following a different social network to the one we have measured
#(i.e. have randomized the connections within groups). This represents the case where our social network does not provide
#a good approximation to the true pathways of social transmission within each group.

#We will look at cTADA first for both datasets (we will assume a constant baseline rate to simplify things)


#Read in the social networks
socNet1<-as.matrix(read.csv(file="multiDiffSocNet1.csv"))
socNet1<-array(socNet1,dim=c(30,30,1))

socNet2<-as.matrix(read.csv(file="multiDiffSocNet2.csv"))
socNet2<-array(socNet2,dim=c(30,30,1))

socNet3<-as.matrix(read.csv(file="multiDiffSocNet3.csv"))
socNet3<-array(socNet3,dim=c(30,30,1))

socNet4<-as.matrix(read.csv(file="multiDiffSocNet4.csv"))
socNet4<-array(socNet4,dim=c(30,30,1))


#DATA SET A: First, we will look at data set A where we have simulated the diffusion following the network

oa1_simFollow<-c(17,21,8,14,28,5,4,16,20,2,7,15,3,30,9,19,1,18,11,12,25,10,22,29,6,26,24,13,23,27)
oa2_simFollow<-c(2,19,25,1,22,14,7,6,21,12,30,9,26,11,28,17,23,5,18,13,20,27,15,24,3,4,29,16,8,10)
oa3_simFollow<-c(8,9,3,1,28,15,29,20,17,11,16,18,10,27,30,23,26,19,24,4,5,12,7,21,14,13,2,6,22,25)
oa4_simFollow<-c(22,1,8,14,9,19,15,2,7,6,11,17,13,25,24,26,30,18,21,10,20,5,27,4,3,12,29,23,16,28)

ta1_simFollow<-c(2,4,13,49,78,107,109,111,114,129,151,190,213,221,234,239,256,258,271,278,287,301,302,304,306,317,330,345,352,353)
ta2_simFollow<-c(4,8,19,20,34,34,40,66,97,113,120,128,128,137,164,179,183,201,205,205,206,239,244,245,249,266,267,273,290,292)
ta3_simFollow<-c(13,14,20,23,34,61,74,92,108,110,128,138,170,174,182,185,186,195,196,211,212,218,219,226,259,277,281,292,297,302)
ta4_simFollow<-c(7,8,12,18,23,26,27,28,31,36,48,51,83,86,104,142,161,172,180,185,189,198,199,211,217,218,239,239,248,249)


#We can create the NBDA objects for cTADA in the usual way
nbdaDataMulti1t_simFollow<-nbdaData(label="MultiDiffusion1_simFollow",assMatrix=socNet1,orderAcq = oa1_simFollow,timeAcq = ta1_simFollow,endTime=max(ta1_simFollow)+1)
nbdaDataMulti2t_simFollow<-nbdaData(label="MultiDiffusion2_simFollow",assMatrix=socNet2,orderAcq = oa2_simFollow,timeAcq = ta2_simFollow,endTime=max(ta2_simFollow)+1)
nbdaDataMulti3t_simFollow<-nbdaData(label="MultiDiffusion3_simFollow",assMatrix=socNet3,orderAcq = oa3_simFollow,timeAcq = ta3_simFollow,endTime=max(ta3_simFollow)+1)
nbdaDataMulti4t_simFollow<-nbdaData(label="MultiDiffusion4_simFollow",assMatrix=socNet4,orderAcq = oa4_simFollow,timeAcq = ta4_simFollow,endTime=max(ta4_simFollow)+1)


#Then fit the model
multiDiffModel1tg_simFollow<-tadaFit(list(nbdaDataMulti1t_simFollow,nbdaDataMulti2t_simFollow,nbdaDataMulti3t_simFollow,nbdaDataMulti4t_simFollow))
data.frame(Variable=multiDiffModel1tg_simFollow@varNames,MLE=multiDiffModel1tg_simFollow@outputPar,SE=multiDiffModel1tg_simFollow@se)

#And an asocial model
multiDiffModel1tg_simFollow_asocial<-tadaFit(list(nbdaDataMulti1t_simFollow,nbdaDataMulti2t_simFollow,nbdaDataMulti3t_simFollow,nbdaDataMulti4t_simFollow),type="asocial")
data.frame(Variable=multiDiffModel1tg_simFollow_asocial@varNames,MLE=multiDiffModel1tg_simFollow_asocial@outputPar,SE=multiDiffModel1tg_simFollow_asocial@se)

#Now we can fit another model with a homogeneous network:

#First create networks full of 1s the same size as the social network for each group:
groupNet1<-array(1,dim=dim(socNet1))
groupNet2<-array(1,dim=dim(socNet2))
groupNet3<-array(1,dim=dim(socNet3))
groupNet4<-array(1,dim=dim(socNet4))

#And then create nbdaData objects for each group with the group network:
nbdaDataMulti1t_simFollow_grNet<-nbdaData(label="MultiDiffusion1_simFollow",assMatrix=groupNet1,orderAcq = oa1_simFollow,timeAcq = ta1_simFollow,endTime=max(ta1_simFollow)+1)
nbdaDataMulti2t_simFollow_grNet<-nbdaData(label="MultiDiffusion2_simFollow",assMatrix=groupNet2,orderAcq = oa2_simFollow,timeAcq = ta2_simFollow,endTime=max(ta2_simFollow)+1)
nbdaDataMulti3t_simFollow_grNet<-nbdaData(label="MultiDiffusion3_simFollow",assMatrix=groupNet3,orderAcq = oa3_simFollow,timeAcq = ta3_simFollow,endTime=max(ta3_simFollow)+1)
nbdaDataMulti4t_simFollow_grNet<-nbdaData(label="MultiDiffusion4_simFollow",assMatrix=groupNet4,orderAcq = oa4_simFollow,timeAcq = ta4_simFollow,endTime=max(ta4_simFollow)+1)

#Now fit the model with the group network
multiDiffModel1tg_simFollow_grNet<-tadaFit(list(nbdaDataMulti1t_simFollow_grNet,nbdaDataMulti2t_simFollow_grNet,nbdaDataMulti3t_simFollow_grNet,nbdaDataMulti4t_simFollow_grNet))
data.frame(Variable=multiDiffModel1tg_simFollow_grNet@varNames,MLE=multiDiffModel1tg_simFollow_grNet@outputPar,SE=multiDiffModel1tg_simFollow_grNet@se)


#And now we can compare AICs
multiDiffModel1tg_simFollow@aicc
multiDiffModel1tg_simFollow_grNet@aicc
multiDiffModel1tg_simFollow_asocial@aicc

#The model with the social network is strongly supported over the group network model, and the asocial model
#So we have evidence that social transmission follows the network within each group



#DATA SET B: Now we will look at data set B where we have simulated the diffusion not following the network

oa1_simNotFollow<-c(6,17,28,1,5,15,23,2,26,18,9,8,25,14,29,19,13,20,12,27,21,22,30,10,4,3,16,11,24,7)
oa2_simNotFollow<-c(7,30,6,10,15,27,3,26,20,16,12,2,17,19,14,28,11,29,21,22,5,24,13,23,8,1,9,18,4,25)
oa3_simNotFollow<-c(7,22,20,6,10,17,5,9,14,18,1,11,19,16,13,3,26,28,24,21,29,12,15,25,8,2,4,30,27,23)
oa4_simNotFollow<-c(24,21,13,27,9,17,15,10,19,23,29,25,7,8,12,11,1,5,22,2,20,16,30,3,4,26,28,18,6,14)

ta1_simNotFollow<-c(1,6,12,24,39,49,74,96,102,103,109,121,134,143,148,152,165,170,184,193,195,214,219,221,223,235,247,284,290,332)
ta2_simNotFollow<-c(1,28,39,45,64,68,80,91,100,115,128,148,154,165,165,167,174,182,196,208,216,225,227,231,236,251,262,263,267,272)
ta3_simNotFollow<-c(4,10,11,25,38,43,50,64,91,93,95,96,105,110,115,117,118,121,121,121,150,156,156,159,165,187,204,204,211,226)
ta4_simNotFollow<-c(1,21,26,35,45,70,110,133,154,170,190,190,200,212,218,221,247,248,256,262,263,275,290,311,312,328,343,344,361,368)


#Create the nbdaData objects as usual for cTADA
nbdaDataMulti1t_simNotFollow<-nbdaData(label="MultiDiffusion1_simNotFollow",assMatrix=socNet1,orderAcq = oa1_simNotFollow,timeAcq = ta1_simNotFollow,endTime=max(ta1_simNotFollow)+1)
nbdaDataMulti2t_simNotFollow<-nbdaData(label="MultiDiffusion2_simNotFollow",assMatrix=socNet2,orderAcq = oa2_simNotFollow,timeAcq = ta2_simNotFollow,endTime=max(ta2_simNotFollow)+1)
nbdaDataMulti3t_simNotFollow<-nbdaData(label="MultiDiffusion3_simNotFollow",assMatrix=socNet3,orderAcq = oa3_simNotFollow,timeAcq = ta3_simNotFollow,endTime=max(ta3_simNotFollow)+1)
nbdaDataMulti4t_simNotFollow<-nbdaData(label="MultiDiffusion4_simNotFollow",assMatrix=socNet4,orderAcq = oa4_simNotFollow,timeAcq = ta4_simNotFollow,endTime=max(ta4_simNotFollow)+1)


#Fit the social network model
multiDiffModel1tg_simNotFollow<-tadaFit(list(nbdaDataMulti1t_simNotFollow,nbdaDataMulti2t_simNotFollow,nbdaDataMulti3t_simNotFollow,nbdaDataMulti4t_simNotFollow))
data.frame(Variable=multiDiffModel1tg_simNotFollow@varNames,MLE=multiDiffModel1tg_simNotFollow@outputPar,SE=multiDiffModel1tg_simNotFollow@se)

#And the asocial model
multiDiffModel1tg_simNotFollow_asocial<-tadaFit(list(nbdaDataMulti1t_simNotFollow,nbdaDataMulti2t_simNotFollow,nbdaDataMulti3t_simNotFollow,nbdaDataMulti4t_simNotFollow),type="asocial")
data.frame(Variable=multiDiffModel1tg_simNotFollow_asocial@varNames,MLE=multiDiffModel1tg_simNotFollow_asocial@outputPar,SE=multiDiffModel1tg_simNotFollow_asocial@se)


#Now the group network model

#As above, create networks full of 1s the same size as the social network for each group:
groupNet1<-array(1,dim=dim(socNet1))
groupNet2<-array(1,dim=dim(socNet2))
groupNet3<-array(1,dim=dim(socNet3))
groupNet4<-array(1,dim=dim(socNet4))

#And then create nbdaData objects for each group with the group network:
nbdaDataMulti1t_simNotFollow_grNet<-nbdaData(label="MultiDiffusion1_simFollow_grNet",assMatrix=groupNet1,orderAcq = oa1_simNotFollow,timeAcq = ta1_simNotFollow,endTime=max(ta1_simNotFollow)+1)
nbdaDataMulti2t_simNotFollow_grNet<-nbdaData(label="MultiDiffusion2_simFollow_grNet",assMatrix=groupNet2,orderAcq = oa2_simNotFollow,timeAcq = ta2_simNotFollow,endTime=max(ta2_simNotFollow)+1)
nbdaDataMulti3t_simNotFollow_grNet<-nbdaData(label="MultiDiffusion3_simFollow_grNet",assMatrix=groupNet3,orderAcq = oa3_simNotFollow,timeAcq = ta3_simNotFollow,endTime=max(ta3_simNotFollow)+1)
nbdaDataMulti4t_simNotFollow_grNet<-nbdaData(label="MultiDiffusion4_simFollow_grNet",assMatrix=groupNet4,orderAcq = oa4_simNotFollow,timeAcq = ta4_simNotFollow,endTime=max(ta4_simNotFollow)+1)

#Now fit the group network model
multiDiffModel1tg_simNotFollow_grNet<-tadaFit(list(nbdaDataMulti1t_simNotFollow_grNet,nbdaDataMulti2t_simNotFollow_grNet,nbdaDataMulti3t_simNotFollow_grNet,nbdaDataMulti4t_simNotFollow_grNet))
data.frame(Variable=multiDiffModel1tg_simNotFollow_grNet@varNames,MLE=multiDiffModel1tg_simNotFollow_grNet@outputPar,SE=multiDiffModel1tg_simNotFollow_grNet@se)

#And compare AICcs
multiDiffModel1tg_simNotFollow@aicc
multiDiffModel1tg_simNotFollow_grNet@aicc
multiDiffModel1tg_simNotFollow_asocial@aicc

#This time the group network model is favoured over the social network model, showing a lack of evidence that the diffusion follows the network
#within each group.

#However, it is strongly favoured over the asocial model. Could this be taken as evidence of social transmission that does not follow the network?
#An alternative possibility is that groups differ in their rate of asocial learning (see above). If this can be ruled out as unlikely a priori
#for the target behvaviour then we suggest that this result could be taken as evidence of social transmission.
#Otherwise, we could include ILVs representing groups to allow for the possibility of different asocial learning rates in the same way as for
#tutorial 6.2 above:

#Indicator for group 2
group2_1<-cbind(rep(0,30))
group2_2<-cbind(rep(1,30))
group2_3<-cbind(rep(0,30))
group2_4<-cbind(rep(0,30))

#Indicator for group 3
group3_1<-cbind(rep(0,30))
group3_2<-cbind(rep(0,30))
group3_3<-cbind(rep(1,30))
group3_4<-cbind(rep(0,30))

#Indicator for group 4
group4_1<-cbind(rep(0,30))
group4_2<-cbind(rep(0,30))
group4_3<-cbind(rep(0,30))
group4_4<-cbind(rep(1,30))

asoc1_group<-c("group2_1","group3_1","group4_1")
asoc2_group<-c("group2_2","group3_2","group4_2")
asoc3_group<-c("group2_3","group3_3","group4_3")
asoc4_group<-c("group2_4","group3_4","group4_4")


#Now re-fit the group network model with group included as an ILV:
nbdaDataMulti1t_simNotFollow_grNet_withGrILV<-nbdaData(label="MultiDiffusion1_simFollow_grNet",assMatrix=groupNet1,orderAcq = oa1_simNotFollow,timeAcq = ta1_simNotFollow,asoc_ilv=asoc1_group,endTime=max(ta1_simNotFollow)+1)
nbdaDataMulti2t_simNotFollow_grNet_withGrILV<-nbdaData(label="MultiDiffusion2_simFollow_grNet",assMatrix=groupNet2,orderAcq = oa2_simNotFollow,timeAcq = ta2_simNotFollow,asoc_ilv=asoc2_group,endTime=max(ta2_simNotFollow)+1)
nbdaDataMulti3t_simNotFollow_grNet_withGrILV<-nbdaData(label="MultiDiffusion3_simFollow_grNet",assMatrix=groupNet3,orderAcq = oa3_simNotFollow,timeAcq = ta3_simNotFollow,asoc_ilv=asoc3_group,endTime=max(ta3_simNotFollow)+1)
nbdaDataMulti4t_simNotFollow_grNet_withGrILV<-nbdaData(label="MultiDiffusion4_simFollow_grNet",assMatrix=groupNet4,orderAcq = oa4_simNotFollow,timeAcq = ta4_simNotFollow,asoc_ilv=asoc4_group,endTime=max(ta4_simNotFollow)+1)

#Now re-fit the group network model
multiDiffModel1tg_simNotFollow_grNet_withGrILV<-tadaFit(list(nbdaDataMulti1t_simNotFollow_grNet_withGrILV,nbdaDataMulti2t_simNotFollow_grNet_withGrILV,nbdaDataMulti3t_simNotFollow_grNet_withGrILV,nbdaDataMulti4t_simNotFollow_grNet_withGrILV))
data.frame(Variable=multiDiffModel1tg_simNotFollow_grNet_withGrILV@varNames,MLE=multiDiffModel1tg_simNotFollow_grNet_withGrILV@outputPar,SE=multiDiffModel1tg_simNotFollow_grNet_withGrILV@se)
#And an asocial model indcluding the group ILV
multiDiffModel1tg_simNotFollow_asocial_withGrILV<-tadaFit(list(nbdaDataMulti1t_simNotFollow_grNet_withGrILV,nbdaDataMulti2t_simNotFollow_grNet_withGrILV,nbdaDataMulti3t_simNotFollow_grNet_withGrILV,nbdaDataMulti4t_simNotFollow_grNet_withGrILV),
                                                          type="asocial")
data.frame(Variable=multiDiffModel1tg_simNotFollow_asocial_withGrILV@varNames,MLE=multiDiffModel1tg_simNotFollow_asocial_withGrILV@outputPar,SE=multiDiffModel1tg_simNotFollow_asocial_withGrILV@se)

#And compare using AICc
multiDiffModel1tg_simNotFollow_grNet_withGrILV@aicc
multiDiffModel1tg_simNotFollow_asocial_withGrILV@aicc

#We can see that the group network model is still favoured even when we allow for the possibility of differences among groups in the asocial rate of learning


#STRATIFIED OADA
#DATASET A (following the network)

#The first stage is to build a big network covering all groups. The dimensions will be 120 x 120 x 1, but we can calculate this
#automatically as follows
c((dim(socNet1)+dim(socNet2)+dim(socNet3)+dim(socNet4))[1:2],1)
#And build an array of 0s of the appropriate size
bigSocNet<-array(0,dim=c((dim(socNet1)+dim(socNet2)+dim(socNet3)+dim(socNet4))[1:2],1))
dim(bigSocNet)
#[1] 120 120   1

#Now we can fill in the appropriate parts of the network as follows
bigSocNet[1:30,1:30,]<-socNet1
bigSocNet[31:60,31:60,]<-socNet2
bigSocNet[61:90,61:90,]<-socNet3
bigSocNet[91:120,91:120,]<-socNet4

#Now we need to get the order of acquisition across each group
#First we need to change the identities of individuals in group 2-4 so we have 31-60 in group 2:
oa2_new<-oa2_simFollow+30
range(oa2_new)
#61-90 in group 3
oa3_new<-oa3_simFollow+60
range(oa3_new)
#91-120 in group 4
oa4_new<-oa4_simFollow+90
range(oa4_new)

#Now let us combine the oa vectors into one
oa_combined_unordered<-c(oa1_simFollow,oa2_new,oa3_new,oa4_new)
#but this is not yet ordered across diffusions. To do this we need to combined the ta vectors into one:
ta_combined<-c(ta1_simFollow,ta2_simFollow,ta3_simFollow,ta4_simFollow)
#Now we can order oa_unordered by ta_combined
oa_combined_Follow<-oa_combined_unordered[order(ta_combined)]


#Now we can make a single nbdaData object for all diffusions:
nbdaDataMultiAll_simFollow<-nbdaData(label="MultiDiffusionAll_simFollow",assMatrix=bigSocNet,orderAcq = oa_combined_Follow)

#Fit the social network model
stratOADAModel_simFollow<-oadaFit(nbdaDataMultiAll_simFollow)
data.frame(Variable=stratOADAModel_simFollow@varNames,MLE=stratOADAModel_simFollow@outputPar,SE=stratOADAModel_simFollow@se)

#Fit the asocial model
stratOADAModel_simFollow_asocial<-oadaFit(nbdaDataMultiAll_simFollow,type="asocial")
data.frame(Variable=stratOADAModel_simFollow_asocial@varNames,MLE=stratOADAModel_simFollow_asocial@outputPar,SE=stratOADAModel_simFollow_asocial@se)
#This is a model with no variables


#Next we create the group network: a single big network but with 1s among individuals in the same group and 0s among individuals
#in different groups
bigSocNet_grNet<-array(0,dim=c(120,120,1))


#Now we can fill in the appropriate parts of the network as follows
bigSocNet_grNet[1:30,1:30,]<-1
bigSocNet_grNet[31:60,31:60,]<-1
bigSocNet_grNet[61:90,61:90,]<-1
bigSocNet_grNet[91:120,91:120,]<-1

#Create the nbdaData object for the group network model
nbdaDataMultiAll_simFollow_grNet<-nbdaData(label="MultiDiffusionAll_simFollow_grNet",assMatrix=bigSocNet_grNet,orderAcq = oa_combined_Follow)

#And fit the group network model to the dataset:
stratOADAModel_simFollow_grNet<-oadaFit(nbdaDataMultiAll_simFollow_grNet)
data.frame(Variable=stratOADAModel_simFollow_grNet@varNames,MLE=stratOADAModel_simFollow_grNet@outputPar,SE=stratOADAModel_simFollow_grNet@se)


#Now we can compare all models using AICc
stratOADAModel_simFollow@aicc
stratOADAModel_simFollow_grNet@aicc
stratOADAModel_simFollow_asocial@aicc

#As with the cTADA, the model with the social network is strongly supported over the group network model, and the asocial model
#So we have evidence that social transmission follows the network within each group


#DATASET B (not following the network)


#Now we need to get the order of acquisition across each group
#First we need to change the identities of individuals in group 2-4 so we have 31-60 in group 2:
oa2_new<-oa2_simNotFollow+30
range(oa2_new)
#61-90 in group 3
oa3_new<-oa3_simNotFollow+60
range(oa3_new)
#91-120 in group 4
oa4_new<-oa4_simNotFollow+90
range(oa4_new)

#Now let us combine the oa vectors into one
oa_combined_unordered<-c(oa1_simNotFollow,oa2_new,oa3_new,oa4_new)
#but this is not yet ordered across diffusions. To do this we need to combined the ta vectors into one:
ta_combined<-c(ta1_simNotFollow,ta2_simNotFollow,ta3_simNotFollow,ta4_simNotFollow)
#Now we can order oa_unordered by ta_combined
oa_combined_simNotFollow<-oa_combined_unordered[order(ta_combined)]


#Now we can make a single nbdaData object for all diffusions:
nbdaDataMultiAll_simNotFollow<-nbdaData(label="MultiDiffusionAll_simNotFollow",assMatrix=bigSocNet,orderAcq = oa_combined_simNotFollow)

#And fit the social networks model:
stratOADAModel_simNotFollow<-oadaFit(nbdaDataMultiAll_simNotFollow)
data.frame(Variable=stratOADAModel_simNotFollow@varNames,MLE=stratOADAModel_simNotFollow@outputPar,SE=stratOADAModel_simNotFollow@se)

#And the asocial model
stratOADAModel_simNotFollow_asocial<-oadaFit(nbdaDataMultiAll_simNotFollow,type="asocial")
data.frame(Variable=stratOADAModel_simNotFollow_asocial@varNames,MLE=stratOADAModel_simNotFollow_asocial@outputPar,SE=stratOADAModel_simNotFollow_asocial@se)
#which has no variables in this case


#Now create the group network as for Data set A
bigSocNet_grNet<-array(0,dim=c(120,120,1))


#Now we can fill in the appropriate parts of the network as follows
bigSocNet_grNet[1:30,1:30,]<-1
bigSocNet_grNet[31:60,31:60,]<-1
bigSocNet_grNet[61:90,61:90,]<-1
bigSocNet_grNet[91:120,91:120,]<-1

#And make an nbdaData object for the group network
nbdaDataMultiAll_simNotFollow_grNet<-nbdaData(label="MultiDiffusionAll_simNotFollow_grNet",assMatrix=bigSocNet_grNet,orderAcq = oa_combined_simNotFollow)

#And fit the group network model to the dataset:
stratOADAModel_simNotFollow_grNet<-oadaFit(nbdaDataMultiAll_simNotFollow_grNet)
data.frame(Variable=stratOADAModel_simNotFollow_grNet@varNames,MLE=stratOADAModel_simNotFollow_grNet@outputPar,SE=stratOADAModel_simNotFollow_grNet@se)

#And compare AICcs again
stratOADAModel_simNotFollow@aicc
stratOADAModel_simNotFollow_grNet@aicc
stratOADAModel_simNotFollow_asocial@aicc

#As with the cTADA, for data set B the group network model is favoured over the social network model, showing a lack of evidence that the diffusion follows the network
#within each group.

#Again, it is strongly favoured over the asocial model. Could this be taken as evidence of social transmission that does not follow the network?
#An alternative possibility is that groups differ in their rate of asocial learning (see above). If this can be ruled out as unlikely a priori
#for the target behvaviour then we suggest that this result could be taken as evidence of social transmission.
#Otherwise, we could include ILVs representing groups to allow for the possibility of different asocial learning rates as for tutorial 6.3:

group2All<-rbind(group2_1,group2_2,group2_3,group2_4)
group3All<-rbind(group3_1,group3_2,group3_3,group3_4)
group4All<-rbind(group4_1,group4_2,group4_3,group4_4)

asocAll_group<-c("group2All","group3All","group4All")

#Now we can make a single nbdaData object for all diffusions:
nbdaDataMultiAll_simNotFollow_grNet__withGrILV<-nbdaData(label="MultiDiffusionAll",assMatrix=bigSocNet_grNet,orderAcq = oa_combined_simNotFollow,asoc_ilv = asocAll_group)


#Now re-fit the group network model with group included as an ILV:
stratOADAModel_simNotFollow_grNet_withGrILV<-oadaFit(nbdaDataMultiAll_simNotFollow_grNet__withGrILV)
data.frame(Variable=stratOADAModel_simNotFollow_grNet_withGrILV@varNames,MLE=stratOADAModel_simNotFollow_grNet_withGrILV@outputPar,SE=stratOADAModel_simNotFollow_grNet_withGrILV@se)

#And an asocial model including the group ILV
stratOADAModel_simNotFollow_asocial_withGrILV<-oadaFit(nbdaDataMultiAll_simNotFollow_grNet__withGrILV,type="asocial")
data.frame(Variable=stratOADAModel_simNotFollow_asocial_withGrILV@varNames,MLE=stratOADAModel_simNotFollow_asocial_withGrILV@outputPar,SE=stratOADAModel_simNotFollow_asocial_withGrILV@se)

#And compare using AICc
stratOADAModel_simNotFollow_grNet_withGrILV@aicc
stratOADAModel_simNotFollow_asocial_withGrILV@aicc

#This time we can see that there is little evidence of social transmission once we have allowed for the possibility of group differences in asocial rate of learning
#In fact, s is estimated at 0. This suggests that the evidence found in the cTADA was a result of the different time course of events, so the robustness of the finding
#may depend on whether it remains when other baselines are used, or if we can a priori assume a constant baseline rate.

#An alternative approach is to compare a model with a group network and no group ILVs, to an asocial model with the group ILVs, i.e.

stratOADAModel_simNotFollow_grNet@aicc
stratOADAModel_simNotFollow_asocial_withGrILV@aicc

#In this case, the asocial model is still favoured. This is despite the data being generated from a model of social transmission, showing that it may be
#difficult to infer social transmission from a group network alone unless we can a priori rule out differences among groups in asocial learning rate.


#############################################################################
# TUTORIAL 6.6
# INCLUDING A SPACE-USE NETWORK
# 1 diffusion
# 1 static network
#############################################################################

#Another situation, similar to the case for a group network in tutorial 6.5, might require a space-use matrix as a null model, if
#different subsets of the population inhabit or use different areas. Only those individuals that use the same areas can have a non-zero connection.
#Therefore, a positive result (evidence for s>0) does not necesarily indicate that the diffusion has followed the provided social network provided, per se,
#or at least only to the extent to which the social network indicates which individuals utilise the same space.

#For example, take a population of birds for whom we have a network giving the frequency with which each dyad engages in allo-preening.
#We fit an NBDA model and find evidence for social transmission. This does not necessarily indicate that preening interactions are an important indicator
#of who learns from whom, since only individuals that inhabit the same space are able to preen one another, and the positive result may only reflect this.

#We can therefore fit an additional model which reflects shared space usage. In many cases, this situation is directly analogous to the group network
#above, but now instead of '1 = same group' and '0= different group', we have '1 = use same space' and '0= use different space'. Some subtle differences may arise if, for 
#example, A shares space with B, and B with C, but B and C do not overlap in space use- in which case dyads do not divide sharply into groups (one
#may wish to use % overlap instead of a binary network in such cases). However, the way in which the two types of network could be included in the 
#analysis remain the same, and the way in which the results are interpreted remain parallel.

#For this reason, we use the use DATASET A from tutorial 6.5 above, but interpret in the context of a space-use network.
#We use an OADA (which becomes equivalent to the stratified OADA above), however the same principles could be used in a TADA.

#Here the bigSocNet built in tutorial 6.5 is the single social network collected by the researcher, in this case for a population of
#120 individuals:
bigSocNet

#And the order of acquisition is 
oa_combined_Follow

#We can make an nbdaData object for the diffusion in the usual way
nbdaDataSpaceUseExample_socNet<-nbdaData(label="nbdaDataSpaceUseExample_socNet",assMatrix=bigSocNet,orderAcq = oa_combined_Follow)

#Fit the social network model
OADAModelSpaceUseExample<-oadaFit(nbdaDataSpaceUseExample_socNet)
data.frame(Variable=OADAModelSpaceUseExample@varNames,MLE=OADAModelSpaceUseExample@outputPar,SE=OADAModelSpaceUseExample@se)

#Fit the asocial model
OADAModelSpaceUseExample_asocial<-oadaFit(nbdaDataSpaceUseExample_socNet,type="asocial")
data.frame(Variable=OADAModelSpaceUseExample_asocial@varNames,MLE=OADAModelSpaceUseExample_asocial@outputPar,SE=OADAModelSpaceUseExample_asocial@se)
#This is a model with no variables


#Next we create the space use network. In this case, we know that individuals 1-30 inhabit the same area, as do 31-60, 61-90 and 91-120.
spaceUseNet<-array(0,dim=c(120,120,1))
spaceUseNet[1:30,1:30,]<-1
spaceUseNet[31:60,31:60,]<-1
spaceUseNet[61:90,61:90,]<-1
spaceUseNet[91:120,91:120,]<-1

#Create the nbdaData object for the space use network model
OADAModelSpaceUseExample_suNet<-nbdaData(label="OADAModelSpaceUseExample_suNet",assMatrix=spaceUseNet,orderAcq = oa_combined_Follow)

#And fit the space use network model to the dataset:
OADAModelSpaceUseExample_suNet<-oadaFit(OADAModelSpaceUseExample_suNet)
data.frame(Variable=OADAModelSpaceUseExample_suNet@varNames,MLE=OADAModelSpaceUseExample_suNet@outputPar,SE=OADAModelSpaceUseExample_suNet@se)


#Now we can compare all models using AICc
OADAModelSpaceUseExample@aicc
OADAModelSpaceUseExample_asocial@aicc
OADAModelSpaceUseExample_suNet@aicc

#The model with the social network is strongly supported over the space use network model, and the asocial model
#So we have evidence that social transmission follows the preening network more closely than the space-use network, suggesting that the
#frequency of preening interactions do have an impact on the likelihood that members of a dyad learn from one another.



