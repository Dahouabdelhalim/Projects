#Load the NBDA package as follows
library(NBDA)

#Set working directory to where the data files are located on your computer



#############################################################################
# TUTORIAL 2.1
# ADDING TIME-CONSTANT ILVs TO AN OADA MODEL
# 1 diffusion
# 1 static network
# 2 Time constant ILVs
#############################################################################

#Read in the social network and order of acquisition vector as shown in tutorial 1.
socNet1<-as.matrix(read.csv(file="exampleStaticSocNet.csv"))
socNet1<-array(socNet1,dim=c(30,30,1))
oa1<-c(26,29,30,8,19,21,22,3,14,12,11,1,17,28,5,9,15,7,6,25,4,13,27,18,20,24,23,16,2,10)


#In our made-up case study, we have the sex and age of all 30 individuals
#We will read in the csv file into a dataframe
ILVdata<-read.csv(file="exampleTimeConstantILVs.csv")

#Then extract each ILV
female<-cbind(ILVdata$female)
#females=1 males=0
age<-cbind(ILVdata$age)
#age in years
#Each ILV needs to be stored as a column matrix; this is to allow extension to time-varying ILVs (see below)

#I am going to standardize age
stAge<-(age-mean(age))/sd(age)
#Since age is centred on 0, and males=0, females=1, the baseline asocial rate is a male of mean age
#Therefore s is estimated relative to this baseline

asoc<-c("female","stAge")
#We create a vector of the names of the variables to be included in the analysis

#Now we create an nbdaData object as before, but specifying the ILVs to be included
#asoc_ilv indicates the ILVs assumed to affect asocial learning rate
#int_ilv indicates the ILVs assumed to affect social learning rate
#multi_ilv indicates the ILVs assumed to affect asocial and social learning rate the same amount (multiplicative model)

#Therefore, if we wanted to fit an "additive" model we would specify:
nbdaData2_add<-nbdaData(label="ExampleDiffusion2",assMatrix=socNet1,orderAcq = oa1,asoc_ilv = asoc)
#Then fit the model:
model2_add<-oadaFit(nbdaData2_add)
#And get the output:
data.frame(Variable=model2_add@varNames,MLE=model2_add@outputPar,SE=model2_add@se)

#             Variable           MLE           SE
#1 1 Social transmission 1  05121.827376 1.176533e+07
#2       2 Asocial: female     11.540250 3.858018e+01
#3        3 Asocial: stAge     -1.027096 1.034692e+00


#If we wanted to fit a "multiplicative" model we would specify:
nbdaData2_multi<-nbdaData(label="ExampleDiffusion2",assMatrix=socNet1,orderAcq = oa1,multi_ilv = asoc)
model2_multi<-oadaFit(nbdaData2_multi)
data.frame(Variable=model2_multi@varNames,MLE=model2_multi@outputPar,SE=model2_multi@se)

#                 Variable       MLE        SE
#1   1 Social transmission 1 7.0160870 9.3076887
#2 2 Social= asocial: female 0.3364630 0.3941035
#3  3 Social= asocial: stAge 0.1345857 0.2041686

#Or the "unconstrained" model
nbdaData2_uc<-nbdaData(label="ExampleDiffusion2",assMatrix=socNet1,orderAcq = oa1,asoc_ilv = asoc,int_ilv=asoc)
model2_uc<-oadaFit(nbdaData2_uc)
data.frame(Variable=model2_uc@varNames,MLE=model2_uc@outputPar,SE=model2_uc@se)

#              Variable           MLE           SE
#1 1 Social transmission 1  1.490592e+05 5.186920e+06
#2       2 Asocial: female  1.093709e+01 3.481992e+01
#3        3 Asocial: stAge -1.173364e+00 9.928235e-01
#4        4 Social: female -4.946070e-02 5.183760e-01
#5         5 Social: stAge  2.985279e-01 2.386751e-01

#Note that one does not need to specify the same set of ILVs in asoc_ilv int_ilv and multi_ilv
#So a model can be fitted in which some variables affect only asocial learning, some only social learning,
#some affect both the same amount, and some affect social and asocial learning differently.

#In practise, there will be little reason to decide, a priori, which ILVs to put in which slot.
#Our preferred approach is to include in multi_ilv only ILVs for which there is a logical reason to believe it will
#affect asocial and social learning the same amount, put all other ILVs in both the asoc_ilv and int_ilv
#slots, and then perform multi-model inferencing (see Tutorial 7)

#For simplicity in the tutorial, we will take the approach of choosing the best of the 3 models based on AICc

model2_add@aicc
model2_multi@aicc
model2_uc@aicc

#We see the additive model is favoured so we will proceed with this model as an example
#(this is not a surprise to us, the data were simulated from a model in which the ILVs
#affect only asocial learning). The output is as follows:

data.frame(Variable=model2_add@varNames,MLE=model2_add@outputPar,SE=model2_add@se)
#            Variable           MLE           SE
#1 1 Social transmission 1 305121.827376 1.176533e+07
#2       2 Asocial: female     11.540250 3.858018e+01
#3        3 Asocial: stAge     -1.027096 1.034692e+00

#And we can compare with an asocial model containing the same ILVs:
model2_asocial<-oadaFit(nbdaData2_add, type="asocial")
model2_add@aicc
#[1] 144.2358
model2_asocial@aicc
#[1] 153.6641

#From looking at the MLEs for the parameters we can see that s is estimated to be very large.
#Indeed if we look at the profile log-likelihood plot for s:

plotProfLik(which=1,model=model2_add,range=c(0,50),resolution=20)
plotProfLik(which=1,model=model2_add,range=c(0,200),resolution=20)
plotProfLik(which=1,model=model2_add,range=c(0,1000),resolution=20)
#We can see it appears to level out as s tend to infinity

#However, this may well be an artifact of the asocial baseline chosen- we can see that
#females are estimated to be much faster than males at asocial learning, and males are set
#as the baseline. This means s is being estimated relative to a very small baseline rate of
#asocial learning
#We can reparameterise the model so that females (of mean age) are the baseline.

male<-1-female
asoc2<-c("male","stAge")
#Now males=1 and females=0

nbdaData3_add<-nbdaData(label="ExampleDiffusion2_reparam",assMatrix=socNet1,orderAcq = oa1,asoc_ilv = asoc2)
#Then fit the model:
model3_add<-oadaFit(nbdaData3_add)
#And get the output:
data.frame(Variable=model3_add@varNames,MLE=model3_add@outputPar,SE=model3_add@se)

#               Variable        MLE           SE
#1 1 Social transmission 1   2.970893     3.922358
#2         2 Asocial: male -19.839490 11618.172259
#3        3 Asocial: stAge  -1.026818     1.034753

#Now we get a large negative coefficient for males, and a much easier to interpret MLE of s=2.97.

#We can also see that the AICc is fractionally better for model 3, showing that the optimum has been
#found more precisely
model2_add@aicc
model3_add@aicc

#Note that the two models specified are the same, just parameterized differently- so you may even
#see the same AICc here.

#You may also think we have somehow magically reduced the importance of social transmission, given the
#much lower estimation of s! But this is not the case, s is merely estimated reletive to a much higher
#baseline rate of asocial learning. This can be seen by comparing %ST for the two models:

nbdaPropSolveByST(model=model2_add)
#P(Network 1)  P(S offset)
#0.90487      0.00000
nbdaPropSolveByST(model=model3_add)
#P(Network 1)  P(S offset)
#0.90491      0.00000

#About the same in each case (and the minor discrepancy is just due to model2_add not quite finding the
#optimum)

#Now we also get a profile log-likelihood for s we can work with:
plotProfLik(which=1,model=model3_add,range=c(0,10),resolution=20)
plotProfLik(which=1,model=model3_add,range=c(0,110),resolution=20)

#We can see the lower limit is between 0 and 2, and the upper from 90-110

profLikCI(which=1,model=model3_add,upperRange = c(90,110),lowerRange = c(0,2))
#Lower CI    Upper CI
#0.4020731 101.4168944


#We can again get the %ST corresponding to the upper and lower points of the 95% CI. However, before we
#do so, we need to look at the constrainedNBDAdata function- so we will come back to this after we have
#looked at the estimated ILV effects and 95% C.I.s for those.


#What about confidence intervals for the coefficients for the ILVs?
#We could get Wald 95% C.I.s by taking MLE +/- 1.96x SE

data.frame(Variable=model3_add@varNames,MLE=model3_add@outputPar,SE=model3_add@se,
           WaldLower=model3_add@outputPar-1.96*model3_add@se,
           WaldUpper=model3_add@outputPar+1.96*model3_add@se)

#We have already seen why the Wald C.I.s are misleading for s
#For age the Wald 95% C.I.s look like they *could* be reasonable
#But for the male effect it looks highly suspect, suggesting that a very large difference in either
#direction is plausible
#The reason for this is that asymmetrical profile log-likelihoods can also arise for ILVs.
#But we can use the same approach to getting profile likelihood C.I.s for the other parameters in the model


#AGE

#We will start with age since that is easier to interpret in this case
#Since age is the third parameter in the model output, we set which=3
plotProfLik(which=3,model=model3_add,range=c(-5,1),resolution=20)

#We can see a little asymmetry here but it is not too bad, suggesting the Wald 95% C.I.s might be
#good in this case, but let us get the profile likelihood intervals anyway. We can see the lower limit is
#around -4 to -3 and the upper limit between 0 and 1

profLikCI(which=3,model=model3_add,lowerRange=c(-4,-3),upperRange = c(0,1))

#Lower CI   Upper CI
#-3.9371721  0.6454884

#So even this mild asymmetry has made a bit of difference when compared to the Wald intervals.
#The first thing to note is that 0 is well within the 95% C.I. so there is not much evidence that age affects
#asocial learning rate.

#Nonetheless, let us examine the back-transformed effect and C.I.
#The MLE for age is -1.027
exp(-1.027)
#[1] 0.3580796
#So the rate of asocial learning decreases by an estimated factor of x0.36 for an increase of 1 S.D. in age, with
#back transformed 95% C.I. of
exp(-3.9371721)
exp(0.6454884)
#0.020 - 1.9x

#So we cannot rule out a large negative effect or a moderate positive effect

#What if we prefered to interpret effect sizes per year, rather than per SD? We simply divide the coefficient by the
#SD for age (original variable, not the standardized version)

exp(-1.027/sd(age))
#[1] 0.6168412
#So the rate of asocial learning decreases by an estimated factor of x0.62 for an increase of 1 year of age, with
#back transformed 95% C.I. of

exp(-3.9371721/sd(age))
exp(0.6454884/sd(age))
#0.157 - 1.35x



#since male is the second parameter in the model output, we set which=2
plotProfLik(which=2,model=model3_add,range=c(-20,5),resolution=20)

#We can see the profile log-likelihood is indeed highly asymmetrical here too, explaining the large SE.
#We can also see there is something a bit odd going on at the left side of the plot, but first let us zoom in and
#find the upper limit

plotProfLik(which=2,model=model3_add,range=c(-5,5),resolution=20)
#We can see the upper limit is somewhere around 0, between -1 and 1. What about the lower (most negative) limit?

plotProfLik(which=2,model=model3_add,range=c(-30,5),resolution=20)

#We can see the profile log-likelihood seems to be going to an asymptote at around 68.6, and then suddenly jumps up
#to a new value at 69.1 and stays there.
#In reality, what is happening here is that the profile log-likelihood *is* approaching an asymptote at
#~68.6. But when the coefficient for male =-15, we are saying that males are
exp(-15)
#[1] 3.059023e-07
#times slower than females at asocial learning. This tiny number causes a computational error meaning the
#profile log-likelihood is not found accurately when the coefficient is less than -14


#What has happened here is only 2 individuals ever learned the behaviour when they had a connection of 0
#to informed individuals, and they were both female. Since this diffusion follows the network, by the logic
#underlying OADA it is plausible that *only* females ever learn asocially, and males always learn by social
#transmission.
#Therefore a value of -Inf is plausible for the male ILV.

#In such cases, all we can do is look for a plausible upper limit for the effect, which we know is between
#-1 and 1

profLikCI(which=2,model=model3_add,upperRange = c(-1,1))
#Lower CI     Upper CI
#NA         -0.001126415

#So the 95% CI only just excludes zero, and there is not huge evidence that females are faster at
#asocial learning than males. We back-transform the effect as follows:

exp(-0.001126415)
#[1] 0.9988742

#This gives us the upper limit of the ratio (male asocial learning rate)/(female asocial learning rate)
#You may find it easier to report by reversing the sign so we are reporting the faster/slower category:

exp(0.001126415)
#[1] 1.001127

#We conclude that females are at least 1.001x faster than males at asocial learning.
#It does not look like a groundbreaking conclusion when stated like this! However, put another way, we
#have found it is *implausible* that males are faster at asocial learning than females, which might be
#enough to cast doubt on a previously favoured hypothesis for that species/ context.

#If you are lucky, all your ILVs will yield clear upper and lower limits making intepretation much easier!


#CONSTRAINING MODELS

#There are a number of reasons we might want to fit a constrained version of a model, e.g.

#1. Constraining one or more parameters to zero
#2. Constraining two or more parameters to have the same value

#Perhaps the most obvious reason to do this is to test a null hypothesis. Imagine we want to get a model
#in which the effect of age is zero. The way to fit such a model is to first create a constrained version
#of the data using the constrainedNBDAdata function

dropAgeData_add<-constrainedNBDAdata(nbdaData3_add,constraintsVect = c(1,2,0))
#Then fit a model to that data
dropAgeModel_add<-oadaFit(dropAgeData_add)
#And get the output:
data.frame(Variable=dropAgeModel_add@varNames,MLE=dropAgeModel_add@outputPar,SE=dropAgeModel_add@se)
#              Variable        MLE           SE
#1 1 Social transmission 1   2.475542     2.978714
#2         2 Asocial: male -19.544137 11683.910820

#We can see that age has been dropped

#The key argument is the constraintsVect. This determines which of the parameters are constrained and how.
#Any variables that have a value of 0 are constrained to be =0 (dropped from the model)
#Any variables that have the same non-zero value are constrained to have the same effect
#Only variables of the same type can be validly constrained, i.e. s parameters; asoc_ilv; int_ilv, multi_ilv

#We can then compare AICc for the two models
model3_add@aicc
dropAgeModel_add@aicc
#showing that AICc is improved by dropping age.

#We can also conduct a likelihood ratio test (LRT)
#Test statistic
2*(dropAgeModel_add@loglik-model3_add@loglik)
#[1] 1.343349

#There are 3 parameters in model3_add, and 2 in dropAgeModel_add, so we have 1 d.f.
pchisq(2*(dropAgeModel_add@loglik-model3_add@loglik),df=1,lower.tail=F)
#[1] 0.2464442
#p= 0.246


#Let us also test for an effect of sex (male)
dropMaleData_add<-constrainedNBDAdata(nbdaData3_add,constraintsVect = c(1,0,2))
dropMaleModel_add<-oadaFit(dropMaleData_add)
data.frame(Variable=dropMaleModel_add@varNames,MLE=dropMaleModel_add@outputPar,SE=dropMaleModel_add@se)
#              Variable        MLE        SE
#1 1 Social transmission 1  6.8853229 9.5740967
#2        2 Asocial: stAge -0.6355158 0.8473629

#We can then compare AICc for the two models
model3_add@aicc
dropMaleModel_add@aicc
#showing that AICc is not improved by dropping sex.

#We can also conduct a likelihood ratio test (LRT)
#Test statistic
2*(dropMaleModel_add@loglik-model3_add@loglik)
#[1] 3.844285

#There are 3 parameters in model3_add, and 2 in dropMaleModel_add, so we have 1 d.f.
pchisq(2*(dropMaleModel_add@loglik-model3_add@loglik),df=1,lower.tail=F)
#[1] 0.04991579
#p= 0.0499
#Reasonable evidence that females are faster asocial learners than males but not overwhelming


#OTHER USES OF constrainedNBDAdata

#If we wanted to fit a model in which a particular parameter is set to a specific value
#we can use the offsetVect argument
#This adds an offset term to the model multiplied by the numbers given in offsetVect
#So here offsetVect=c(0,0,2) would add an offset of age x 2
#If we fit a model with the effect of age constrained to be zero in constraintsVect = c(1,2,0)
#but with offsetVect=c(0,0,2), we are constraining the coefficient of age to be 2

ageTwoData_add<-constrainedNBDAdata(nbdaData3_add,constraintsVect = c(1,2,0),offsetVect = c(0,0,2))
ageTwoModel_add<-oadaFit(ageTwoData_add)

#In this example we set the effect of male to -3 and the effect of age to 2:

ageTwoMaleMinusThreeData_add<-constrainedNBDAdata(nbdaData3_add,constraintsVect = c(1,0,0),offsetVect = c(0,-3,2))
ageTwoMaleMinusThreeModel_add<-oadaFit(ageTwoMaleMinusThreeData_add)

#We can also constrain parameters of the same type to be equal, e.g.

ageAndSexEqualData_add<-constrainedNBDAdata(nbdaData3_add,constraintsVect = c(1,2,2))
ageAndSexEqualModel_add<-oadaFit(ageAndSexEqualData_add)

#constrains the effect of age and sex to be the same. There is no reason we would want to do this
#for age and sex but we will see examples later where it makes sense.
#Or we can fit a model with a specified difference between two parameters

ageMinusSex1Data_add<-constrainedNBDAdata(nbdaData3_add,constraintsVect = c(1,2,2), offsetVect=c(0,0,1))
ageMinusSex1Model_add<-oadaFit(ageMinusSex1Data_add)

#This fits a model in which coefficient for age - coefficient for sex = 1
#Again there is little reason to do this for these two parameters

#In models with more than one network, we can also constrain s parameters to be equal to zero, equal to another number,
#constrain s parameters to be the same or have a specified difference between them
#However, if all s parameters are constrained to 0 this will trip an error:

sZeroData_add<-constrainedNBDAdata(nbdaData3_add,constraintsVect = c(0,1,2))

#So if we wish to constrain all s parameters to be 0, or constrained to some value, we have to fit an asocial model.
#To fit a model with s constrained to 0.4020731 (the lower bound of the 95% C.I.) we set offsetVect=c(0.4020731,0,0)

sLowerData_add<-constrainedNBDAdata(nbdaData3_add,constraintsVect = c(1,2,3),offsetVect=c(0.4020731,0,0))

#Then fit an asocial model
sLowerModel_add<-oadaFit(sLowerData_add, type="asocial")


#We are now in a position to return to the problem we had above, of estimating %ST for the upper and lower bounds
#of the 95% C.I.
#Lower CI    Upper CI
#0.4020731 101.4168944

#In order to estimate %ST, we need to input values for every parameter in the model. So we need to find the
#values of the other parameters that correspond to the upper and lower limits of s.
#We do this by fitting models in which s is constrained to its upper and lower limits:

#Lower limit first:
sLowerData_add<-constrainedNBDAdata(nbdaData3_add,constraintsVect = c(1,2,3),offsetVect=c(0.4020731,0,0))
sLowerModel_add<-oadaFit(sLowerData_add, type="asocial")
#We can recover the values of the other two parameters as follows:
sLowerModel_add@outputPar
#And input these to get %ST
nbdaPropSolveByST(par=c(0.4020731,sLowerModel_add@outputPar),nbdadata = nbdaData3_add)
#P(Network 1)  P(S offset)
#0.75866      0.00000

#Alternatively, we can set model= constrained model to get the same result
nbdaPropSolveByST(model=sLowerModel_add)
#P(S offset)
#0.75866
#Note now %ST is assigned to the offset for s.

#Upper limit next
sUpperData_add<-constrainedNBDAdata(nbdaData3_add,constraintsVect = c(1,2,3),offsetVect=c(101.4168944,0,0))
sUpperModel_add<-oadaFit(sUpperData_add, type="asocial")
#We can recover the values of the other two parameters as follows:
sUpperModel_add@outputPar
#And input these to get %ST
nbdaPropSolveByST(par=c(101.4168944,sUpperModel_add@outputPar),nbdadata = nbdaData3_add)
#P(Network 1)  P(S offset)
#0.95632      0.00000

#Alternatively, we can set model= constrained model to get the same result
nbdaPropSolveByST(model=sUpperModel_add)
#P(S offset)
#0.95632

#So in this case we estimate that 75.9 - 95.6 % of events were by social transmission


#############################################################################
# TUTORIAL 2.2
# ADDING TIME-VARYING ILVs TO AN OADA MODEL
# 1 diffusion
# 1 static network
# 1 Time constant ILV
# 1 Time varying ILV
#############################################################################

#Time-varying ILVs are easily added to an OADA
#It might take a bit of work to set up the nbdaData object, but once this is done the analysis proceeds as above.

#If one or more of our ILVs is time-varying, we need to set up an array for ALL of our ILVs specifying their
#values for every individual for every acquistion event.

#For example, let us assume that individuals 1-5 and 16-20 were undergoing some kind of experimental manipulation for events 1-15
#and 6-10 and 21-30 were undergoing the same experimental manipulation for events 16-30.
#Let us set up a matrix with rows = number of individuals, columns = number of acquistion event (in this case, 30 for both)

treatment<-matrix(NA,nrow=30,ncol=30)
treatment[,1]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,2]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,3]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,4]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,5]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,6]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,7]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,8]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,9]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,10]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,11]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,12]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,13]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,14]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,15]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))

treatment[,16]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,17]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,18]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,19]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,20]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,21]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,22]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,23]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,24]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,25]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,26]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,27]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,28]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,29]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,30]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment

#let us also assume we want to include the male variable from above
male
#We have to input this as a time varying ILV as well, even though it does not change
#However it is easy to create this using the byrow=F argument:
maleTV<-matrix(male,nrow=30,ncol=30,byrow = F)
maleTV

asocTV<-c("treatment","maleTV")

#We can then create the nbdaData object for the unconstrained model, specifying asocialTreatment="timevarying"

nbdaData3_TVadd<-nbdaData(label="ExampleDiffusion2",assMatrix=socNet1,orderAcq = oa1,asoc_ilv = asocTV,int_ilv = asocTV,asocialTreatment = "timevarying")
#Fit the model
model_socialTV<-oadaFit(nbdaData3_TVadd)
#Display the output
data.frame(Variable=model_socialTV@varNames,MLE=model_socialTV@outputPar,SE=model_socialTV@se)




