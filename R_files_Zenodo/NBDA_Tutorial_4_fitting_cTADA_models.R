#Load the NBDA package as follows
library(NBDA)

#Set working directory to where the data files are located on your computer



#############################################################################
# TUTORIAL 4.1
# FITTING A cTADA MODEL
# 1 diffusion
# 1 static network
# 1 Time constant ILV
#############################################################################

#Read in the social network and order of acquisition vector as shown in Tutorial 1.
socNet1<-as.matrix(read.csv(file="exampleStaticSocNet.csv"))
socNet1<-array(socNet1,dim=c(30,30,1))
oa1<-c(26,29,30,8,19,21,22,3,14,12,11,1,17,28,5,9,15,7,6,25,4,13,27,18,20,24,23,16,2,10)

#Load and process the ILVs as shown in Tutorial 2: we will just include age to keep things simple
ILVdata<-read.csv(file="exampleTimeConstantILVs.csv")
age<-cbind(ILVdata$age)
stAge<-(age-mean(age))/sd(age)
asoc<-c("stAge")

#Now in addition, we need a vector giving the times at which each individual acquired the behaviour:
ta1<-c(234,252,262,266,273,296,298,310,313,326,332,334,334,337,338,340,343,367,374,376,377,393,402,405,407,407,435,472,499,567)

#We also need to specify a time at which the diffusion ends. If the diffusion is complete--i.e. everyone in the population learns the behaviour
#you can input any value greater than the time at which the last individual learned and it makes no difference.
nbdaData_cTADA<-nbdaData(label="ExampleDiffusion_cTADA",assMatrix=socNet1,orderAcq = oa1,timeAcq = ta1, endTime = 568,asoc_ilv = asoc, int_ilv=asoc)

#We can still choose to fit an OADA model to this object, and ignore the time data:
model_oada<-oadaFit(nbdaData_cTADA)
model_oada_asocial<-oadaFit(nbdaData_cTADA,type="asocial")

#We fit a cTADA model using tadaFit:
model_constant<-tadaFit(nbdaData_cTADA)
#And get the output:
data.frame(Variable=model_constant@varNames,MLE=model_constant@outputPar,SE=model_constant@se)
#              Variable          MLE           SE
#1         Scale (1/rate): 4734.4023556 4016.0905054
#2 1 Social transmission 1   16.8826352   14.9552133
#3        2 Asocial: stAge   -0.7520922    0.8683823
#4         3 Social: stAge    0.3035350    0.2151248

#By default a constant baseline rate is fitted, so we have a single parameter lambda0 or rate.
#However, the optimization algorithm works better if we specify the value of 1/rate- this is the first
#parameter estimated in the table above. This can be thought of as the estimated mean time to learn asocially,
#at the baseline asocial level (in this case, this is age set to its mean)

#The numbering of the variables starts after the baseline function parameters, so we have
#"1 Social transmission 1" . This means when we are adding constraints and offsets the baseline parameters are ignored as we will see below.

#Let us next fit an asocial model with a constant baseline rate:
model_constant_asocial<-tadaFit(nbdaData_cTADA,type="asocial")
data.frame(Variable=model_constant_asocial@varNames,MLE=model_constant_asocial@outputPar,SE=model_constant_asocial@se)
#       Variable           MLE         SE
#1  Scale (1/rate): 355.000154854 64.6994285
#2 1 Asocial: stAge  -0.005259915  0.1859458

#Note that not only is the s parameter missing, but also the effect of age on social learning- since this parameter is meaningless when s=0


#And compare AICcs
model_constant@aicc
#[1] 349.1325
model_constant_asocial@aicc
#[1] 416.8778

exp(0.5*(model_constant_asocial@aicc-model_constant@aicc))
#[1] 5.136959e+14

#Far more support for social model, and a much bigger difference than for OADA:
model_oada@aicc
model_oada_asocial@aicc
#Showing how we can get far more power from a cTADA if the assumptions are true.

#However, maybe these results are simply the result of an increasing baseline rate?
#We can fit two models which allow for an increasing baseline rate, the gamma and Weibull models:

model_gamma<-tadaFit(nbdaData_cTADA, baseline = "gamma")
#This model can take a long time to fit relative to the other two options
data.frame(Variable=model_gamma@varNames,MLE=model_gamma@outputPar,SE=model_gamma@se)
#              Variable        MLE         SE
#1         Scale (1/rate): 31.8749336 23.3470376
#2                   Shape 13.2441307  7.5801080
#3 1 Social transmission 1  0.2216704  0.3800373
#4        2 Asocial: stAge -0.7376322  0.5477921
#5         3 Social: stAge  1.3152773  0.7394257

#Here we have another baseline rate parameter, shape. If this is >1 it corresponds to an increasing baseline function
#and <1 to a decreasing baseline function.

#Let us fit the asocial model:
model_gamma_asocial<-tadaFit(nbdaData_cTADA, baseline = "gamma", type="asocial")

#And compare AICcs
model_gamma@aicc
#[1] 346.1326
model_gamma_asocial@aicc
#[1] 347.2452
exp(0.5*(model_gamma_asocial@aicc-model_gamma@aicc))
#[1] 1.744256

#So if we allow for an increasing baseline rate, the support for social learning is greatly weakened.


#We will also fit the Weibull model:
model_weibull<-tadaFit(nbdaData_cTADA, baseline = "weibull")
data.frame(Variable=model_weibull@varNames,MLE=model_weibull@outputPar,SE=model_weibull@se)
#             Variable           MLE           SE
#1         Scale (1/rate):  4.037055e+02 20.158151117
#2                   Shape  4.763533e+00  0.662451787
#3 1 Social transmission 1  8.777549e-04  0.001389945
#4        2 Asocial: stAge -2.598330e-01  0.290078859
#5         3 Social: stAge  4.562159e+00  1.062500396

#Again we have a shape parameter, which is >1 for an increasing baseline rate and <1 for a decreasing baseline rate

#Let us fit the asocial model:
model_weibull_asocial<-tadaFit(nbdaData_cTADA, baseline = "weibull", type="asocial")

#And compare AICcs
model_weibull@aicc
#[1] 352.2198
model_weibull_asocial@aicc
#[1] 352.8123
exp(0.5*(model_weibull_asocial@aicc-model_weibull@aicc))
#[1] 1.344825

#Again, there is much less evidence for social learning if we assume a Weibull baseline rate

#In terms of AICc, the gamma model is favoured, so let us use this for further inference
#Inference proceeds as for the OADA in the previous tutorials. The only difference here is that we ignore the baseline parameters
#when constraining models, and specifying "which" parameter when finding C.I.s

#So to get the 95% C.I. for s:
#s is the first parameter, not the 3rd. This is why it is labelled "1 Social transmission" in the table above
plotProfLik(which=1, model=model_gamma, range=c(0,10),resolution=20)
#upper limit between 8 and 10. Difficult to see the lower limit, so we can zoom in:
plotProfLik(which=1, model=model_gamma, range=c(0,1),resolution=20)
#between 0 and 0.1
profLikCI(which=1, model=model_gamma, lowerRange = c(0,0.1),upperRange=c(8,10))
#Lower CI    Upper CI
#0.003444185 8.717068774


#Note that we can often get warnings like
#Warning messages:
#  1: In dgamma(solveTimes, shape = shape, rate = rate) : NaNs produced
#  2: In pgamma(solveTimes, shape = shape, rate = rate, lower = FALSE) :
#  NaNs produced

#when fitting a TADA with non-constant baseline rate. This just means the optimization algorithm is "wandering" into regions
#of parameter space for which the values of shape and rate are unrealistic and return NaN. So long as the model has fitted this
#is not a concern.

#In cases where the results seem to be heavily influenced by whether we assume the baseline rate is constant or not
#it suggests that the analysis is dominated by the time course of events rather than the pattern of diffusion through
#the network. Here, we can see that we seem to have less power in the cTADA then the OADA, once we have allowed for an increasing
#baseline rate, since s is estimated very low, as is %ST:

nbdaPropSolveByST(model=model_gamma)
#P(Network 1)  P(S offset)
#0.44494      0.00000

nbdaPropSolveByST(model=model_oada)
#P(Network 1)  P(S offset)
#0.84249      0.00000

#In such cases, we think there is a good case for switching to OADA which is only sensitive to the pattern of spread through the
#network.




#############################################################################
# TUTORIAL 4.2
# ADDING TIME-VARYING ILVs TO A cTADA MODEL
# 1 diffusion
# 1 static network
# 1 Time constant ILV
# 1 Time varying ILV
#############################################################################

#Time-varying ILVs are easily added to a cTADA
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

#let us also assume we want to include the male variable we created in Tutorial 2.1

female<-ILVdata$female
male<-1-female
#Now males=1 and females=0

male
#We have to input this as a time varying ILV as well, even though it does not change
#However it is easy to create this using the byrow=F argument:
maleTV<-matrix(male,nrow=30,ncol=30,byrow = F)
maleTV

asocTV<-c("treatment","maleTV")

#Now in addition, we need a vector giving the times at which each individual acquired the behaviour:
oa2<-c(6,16,9,11,2,8,14,21,26,3,15,18,23,12,25,27,5,1,4,28,29,30,24,20,10,13,7,17,22,19)
ta2<-c(40,53,75,91,92,135,144,146,158,178,181,191,196,197,206,208,212,212,217,235,235,237,241,243,244,250,261,265,265,267)
#And since all individuals learn the endTime can be set to the maximum +1 =268

#We can then create the nbdaData object for the unconstrained model, specifying asocialTreatment="timevarying"

nbdaData3_TV<-nbdaData(label="ExampleDiffusion2",assMatrix=socNet1,orderAcq = oa2,timeAcq = ta2, endTime = 268,
                          asoc_ilv = asocTV,int_ilv = asocTV,asocialTreatment = "timevarying")
#Fit the model
model_socialTV<-tadaFit(nbdaData3_TV)
#Display the output
data.frame(Variable=model_socialTV@varNames,MLE=model_socialTV@outputPar,SE=model_socialTV@se)

#############################################################################
# TUTORIAL 4.3
# COMPARISON TO A HOMOGENEOUS NETWORK IN A cTADA MODEL
# 1 diffusion
# 1 static network
# 1 Time constant ILV
#############################################################################

#A positive result in a cTADA (evidence for s>0) does not necessarily constitute strong evidence that the diffusion follows the specific network 
#provided (see main text).
#To test for evidence of this, we can fit an additional model with a homogeneous network and compare its fit to the model with the provided social
#network.

#Fit the constant baseline model from Tutorial 4.1:
model_constant<-tadaFit(nbdaData_cTADA)
#And get the output:
data.frame(Variable=model_constant@varNames,MLE=model_constant@outputPar,SE=model_constant@se)

#A homogeneous network can be created as an array the same size as socNet1 but with 1s in every cell:
homoNet1<-array(1,dim=dim(socNet1))
#We can then create a second nbdaData object with the homogeneous network in place of the 
nbdaData_cTADA_homog<-nbdaData(label="ExampleDiffusion_cTADA",assMatrix=homoNet1,orderAcq = oa1,timeAcq = ta1, endTime = 568,asoc_ilv = asoc, int_ilv=asoc)
#And then we can fit a model to this nbdaData object:
model_constant_homog<-tadaFit(nbdaData_cTADA_homog)
#And get the output as follows:
data.frame(Variable=model_constant_homog@varNames,MLE=model_constant_homog@outputPar,SE=model_constant_homog@se)

model_constant@aicc
model_constant_homog@aicc

#Here we can see the homogeneous network receives less support than the network model, showing evidence that the diffusion follows the network.

#We can also fit models with the social and homogeneous networks with gamma and weibull baseline rates:

model_gamma<-tadaFit(nbdaData_cTADA, baseline = "gamma")
model_gamma_homog<-tadaFit(nbdaData_cTADA_homog, baseline = "gamma")
model_weibull<-tadaFit(nbdaData_cTADA, baseline = "weibull")
model_weibull_homog<-tadaFit(nbdaData_cTADA_homog, baseline = "weibull")

model_constant@aicc
model_constant_homog@aicc
model_gamma@aicc
model_gamma_homog@aicc
model_weibull@aicc
model_weibull_homog@aicc

#Now we can see the model with the gamma baseline rate and social network is still the favoured model, with

exp(0.5*(model_gamma_homog@aicc-model_gamma@aicc))

#6.9 times more support for social transmission following the social network than for the homogeneous network.