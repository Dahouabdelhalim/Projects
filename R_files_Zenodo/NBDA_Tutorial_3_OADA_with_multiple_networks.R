#Load the NBDA package as follows
library(NBDA)

#Set working directory to where the data files are located on your computer



#############################################################################
# TUTORIAL 3.1
# INCLUDING MULTIPLE NETWORKS IN AN OADA MODEL
# 1 diffusion
# 2 static networks
# 2 Time constant ILVs
#############################################################################

#Read in the 2 social networks and order of acquisition vector as shown in Tutorial 1
socNet1<-as.matrix(read.csv(file="exampleStaticSocNet.csv"))
socNet2<-as.matrix(read.csv(file="exampleStaticSocNet2.csv"))
oa1<-c(26,29,30,8,19,21,22,3,14,12,11,1,17,28,5,9,15,7,6,25,4,13,27,18,20,24,23,16,2,10)

#Load and process the ILVs as shown in Tutorial 2
ILVdata<-read.csv(file="exampleTimeConstantILVs.csv")
female<-cbind(ILVdata$female)
male<-1-female
age<-cbind(ILVdata$age)
stAge<-(age-mean(age))/sd(age)
asoc<-c("male","stAge")

#In a multi-network NBDA, we need to combine our networks into an array
#If all of our networks are static (i.e. do not change over time), this is a three-dimensional array of size 
#no. individuals x no.individuals x no.networks
socNets<-array(NA,dim=c(30,30,2))
#Then slot the networks into the array:
socNets[,,1]<-socNet1
socNets[,,2]<-socNet2

#ASIDE:
#If we have dynamic (time-varying) networks we need to create a four dimensional array of size
#no. individuals x no.individuals x no.networks x number of time periods
#and provide an assMatrixIndex vector as shown in Tutorial 1.

#If only one of the networks in the model is dynamic, we need to treat the other static networks as if they were dynamic too
#and enter the same connection values into each time period.

#Then we go on to create our nbdaData object as before. Let us assume that we are interested in the unconstrained model:
nbdaData_multiNet<-nbdaData(label="ExampleDiffusion2",assMatrix=socNets,orderAcq = oa1,asoc_ilv = asoc,int_ilv = asoc )
#Then fit the model:
model1_multiNet<-oadaFit(nbdaData_multiNet)
#And get the output:
data.frame(Variable=model1_multiNet@varNames,MLE=model1_multiNet@outputPar,SE=model1_multiNet@se)

#              Variable          MLE           SE
#1 1 Social transmission 1   2.52517170 4.590037e+00
#2 2 Social transmission 2   0.00000000 1.105857e+00
#3         3 Asocial: male -19.32342498 1.112089e+04
#4        4 Asocial: stAge  -1.17332662 1.174882e+00
#5          5 Social: male   0.04954004 5.366979e-01
#6         6 Social: stAge   0.29852561 2.864485e-01

#You will notice that we have a second s parameter in the model. The s parameters correspond to the order the networks
#are entered into the socNets array. So the first s parameter,s1, labelled "1 Social transmission 1"  corresponds to 
#the network in socNets[,,1] and the second s parameter,s2, labelled "2 Social transmission 2" to the network in 
#socNets[,,2].

#Inference regarding the ILVs can proceed as it did in Tutorial 2, so here we will focus on inference about the s parameters.

#It looks like we have some evidence of social transmission following network 1 but not network 2, at first sight.
#Let us fit some constrained models to test some null hypotheses.


#First let us test the hypothesis s1=0.
#We first create the nbdaData object with the constraint s1=0:
s1equals0_data<-constrainedNBDAdata(nbdaData_multiNet,constraintsVect = c(0,1,2,3,4,5))
#Here the first parameter, s1, is constrained to be 0 and all other parameters are unconstrained
#Then fit the model:
s1equals0_model<-oadaFit(s1equals0_data)

#We can then compare AICcs
model1_multiNet@aicc
#[1] 151.2664
s1equals0_model@aicc
#[1] 157.4304
#The model with s1>0 is favoured but by how much?
exp(0.5*(s1equals0_model@aicc-model1_multiNet@aicc))
#[1] 21.80218
#21x more support for a model in which there is social transmission following network 1

#We can also conduct a likelihood ratio test (LRT)
#Test statistic
2*(s1equals0_model@loglik-model1_multiNet@loglik)
#The difference in number of parameters is 1, so df=1
pchisq(2*(s1equals0_model@loglik-model1_multiNet@loglik),df=1,lower.tail = F)
#[1] 0.002271373
#Strong evidence for social transmission following network 1



#Now we can do the same for the hypothesis s2=0
s2equals0_data<-constrainedNBDAdata(nbdaData_multiNet,constraintsVect = c(1,0,2,3,4,5))
s2equals0_model<-oadaFit(s2equals0_data)
model1_multiNet@aicc
#[1] 151.2664
s2equals0_model@aicc
#[1] 148.1142
#The model with s2=0 is favoured but by how much?
exp(0.5*(model1_multiNet@aicc-s2equals0_model@aicc))
#[1] 4.835995

#Now the LRT
#Test statistic
2*(s2equals0_model@loglik-model1_multiNet@loglik)
#The difference in number of parameters is 1, so df=1
pchisq(2*(s2equals0_model@loglik-model1_multiNet@loglik),df=1,lower.tail = F)
#[1] 0.9999791
#No evidence for social transmission following network 2


#This does not necessarily mean that we have strong evidence for s1>s2. It could be that s2 has very wide
#confidence intervals.
#We need to test the hypothesis s1=s2 separately, as follows
s1equalss2_data<-constrainedNBDAdata(nbdaData_multiNet,constraintsVect = c(1,1,2,3,4,5))
#In the constraintsVect, parameter 1 and 2 have the same number, meaning they are constrained to have the
#same value, i.e. s1=s2 as required
#Fit the model:
s1equalss2_model<-oadaFit(s1equalss2_data)
#compare AICcs:
model1_multiNet@aicc
s1equalss2_model@aicc
exp(0.5*(s1equalss2_model@aicc-model1_multiNet@aicc))
#[1] 1.130765

#LRT:
#Test statistic
2*(s1equalss2_model@loglik-model1_multiNet@loglik)
#The difference in number of parameters is 1, so df=1
pchisq(2*(s1equalss2_model@loglik-model1_multiNet@loglik),df=1,lower.tail = F)
#[1] 0.06527703
#Only weak evidence that s1>s2.

#So overall we have:
#1. strong evidence for social transmission following network 1. 
#2. no evidence for social transmission following network 2.
#3. weak evidence that the social transmission following network 2, if it exists, is weaker than the social 
#transmission following network 1

#In general terms: it is often tempting, when we find strong evidence of effect A, and no evidence of effect B
#to conclude that we have strong evidence that effect A > effect B.
#But this is a logical error--remember that "no evidence of an effect" does not equate to "strong evidence of no effect".
#This made-up example illustrates this point nicely.

#We would strongly advise users of NBDA to present confidence intervals (C.I.s) for effect sizes. These can
#be obtained for s1 and s2 in the same manner as in Tutorial 2:

#for s1, which=1
plotProfLik(which=1,model=model1_multiNet, range=c(0,10),resolution=20)
#lower limit between 0 and 1
plotProfLik(which=1,model=model1_multiNet, range=c(0,150),resolution=20)
#upper between 100 and 130
profLikCI(which=1,model=model1_multiNet,lowerRange = c(0,1),upperRange = c(100,130))
#Lower CI     Upper CI 
#0.01910204 108.60099419 

#for s2, which=2
plotProfLik(which=2,model=model1_multiNet, range=c(0,10),resolution=20)
#s2=0 is below the line so we know the lower limit is 0, since s2 can only take values greater than or equal to zero.
#The upper limit is between 7 and 9
profLikCI(which=2,model=model1_multiNet,upperRange = c(7,9))
#Lower CI Upper CI 
#0.00000  7.63635  

#So we have a clearer picture already as to why we saw the pattern of significance in the LRTs
#There is quite a lot of overlap in the 95% C.I.s.
#However, do not fall into the trap of thinking that because the C.I.s for two parameters overlap there is NOT
#evidence for a difference between them!
#Instead one should aim to get the confidence interval for the difference, in this case for s1-s2.

#We can do this using the constraintsVect argument. By setting this to constrain s1=s2, and by setting
#which=1, we are obtaining the profile log-likelihood for the s1-s2. If we set which=2 we would get the profile
#log-likelihood for s2-s1, as we will see below.

#In general, if we want the C.I.s for the difference between 2 parameters of the same type, we constrain them to be the
#same in constraintsVect and specify one of them as "which".
#So this same approach can be used to get C.I.s for the difference in the effects of 2 ILVs on asocial learning, on
#social learning, or the effects of two ILVs that are both assumed to have a multiplicative effect.

plotProfLik(which=1,model=model1_multiNet, range=c(0,120),resolution=20,constraintsVect = c(1,1,2,3,4,5))
#The upper limit is between 100 and 120. We can find this as follows:
profLikCI(which=1,model=model1_multiNet,upperRange = c(100,120),constraintsVect = c(1,1,2,3,4,5))
#Lower CI Upper CI 
#0.000  108.601 

#The lower limit has been identified as 0 since 0 is below the line. But this is not the case since s1-s2 can take values
#less than 0. We can find the lower limit by repeating the process for s2-s1, and reversing the sign of the upper limit found:

plotProfLik(which=2,model=model1_multiNet, range=c(0,2),resolution=20,constraintsVect = c(1,1,2,3,4,5))
profLikCI(which=2,model=model1_multiNet,upperRange = c(0,0.5),constraintsVect = c(1,1,2,3,4,5))
#Lower CI  Upper CI 
#0.0000000 0.1351511 

#So the 95% C.I. for s1-s2 is -0.13 to 108.6:
#It is plausible that s1 is a great deal larger than s2, but it is implausible that s2 is anything but a tiny 
#amount larger than s1.

#You might wonder why we did not simply extend the range for s1-s2 to include negative values thus:
plotProfLik(which=1,model=model1_multiNet, range=c(-1,120),resolution=20,constraintsVect = c(1,1,2,3,4,5))
#You can see that we get an error in this case.
#This is because doing things this way allows s1 to get pushed into negative values during optimization, which trips errors.
#In some cases, you might get lucky and it works, but it is better to rely on the 2 step process above if the C.I. includes
#zero.
#There is no problem in setting a range with negative values for difference in the ILVs, just for the s parameters.


#We can also get %ST for the proportion of events that ocurred via social transmission for each network:

nbdaPropSolveByST(model=model1_multiNet)
#P(Network 1) P(Network 2)  P(S offset) 
#0.88848      0.00000      0.00000 

#If we want to get %ST corresponding to the upper and lower limits of C.I.s we can do so using the same procedure as in previous tutorials:



#s1 
#Lower CI     Upper CI 
#0.01910204 108.60099419 

#lower limit
s1AtLower_data<-constrainedNBDAdata(nbdaData_multiNet,constraintsVect = c(0,1,2,3,4,5),offsetVect = c(0.01910204,0,0,0,0,0))
s1AtLower_model<-oadaFit(s1AtLower_data)
nbdaPropSolveByST(model=s1AtLower_model)
#P(Network 1)  P(S offset) 
#0.0000       0.3313 
#Note that network 1 is now moved to the offset, so we have 33.1% at the lower limit for s1

#upper limit
s1AtUpper_data<-constrainedNBDAdata(nbdaData_multiNet,constraintsVect = c(0,1,2,3,4,5),offsetVect = c(108.60099419,0,0,0,0,0))
s1AtUpper_model<-oadaFit(s1AtUpper_data)
nbdaPropSolveByST(model=s1AtUpper_model)
#P(Network 1)  P(S offset) 
#0.00000      0.95568 

#So a plausible range for the % of events occurring by social transmission through network 1 is 33.1 - 95.6%



#s2
#Lower CI Upper CI 
#0.00000  7.63635  

#lower limit is 0 so must correspond to 0%

#upper limit
s2AtUpper_data<-constrainedNBDAdata(nbdaData_multiNet,constraintsVect = c(1,0,2,3,4,5),offsetVect = c(0,7.63635,0,0,0,0))
s2AtUpper_model<-oadaFit(s2AtUpper_data)
nbdaPropSolveByST(model=s2AtUpper_model)
#P(Network 1)  P(S offset) 
#0.7760       0.1717 

#So a plausible range for the % of events occurring by social transmission through network 2 is 0 - 17.2%

#Aside: why do these ranges for %ST not overlap when the 95% C.I.s for s do?
#Recall that %ST takes into account the connections of the network, so this can occur if one network has more/stronger connections than the other.
#In this case the reason is more subtle. Network 2 is actually just a randomized version of network 1, so they have the same number and
#strength of connections. %ST also depends on the pathway that the diffusion takes through the network, i.e. whether it passes through
#parts of the network with strong or weak connections.

