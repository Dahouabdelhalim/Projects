#Load the NBDA package as follows
library(NBDA)

#Set working directory to where the data files are located on your computer



#############################################################################
# TUTORIAL 5.1
# FITTING A dTADA MODEL
# 1 diffusion
# 1 static network
# 1 Time constant ILV
#############################################################################

#It may be that you already have your data collected in discrete time periods. However, as explained in the main text, there are circumstances
#under which you may wish to convert continuous TADA data into discrete TADA data, so we will start with that.

#We will first load in the data we used in Tutorial 4 with the cTADA:

#######################################################################
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

#Create the nbdaData object for comparison
nbdaData_cTADA<-nbdaData(label="ExampleDiffusion_cTADA",assMatrix=socNet1,orderAcq = oa1,timeAcq = ta1, endTime = 568,asoc_ilv = asoc, int_ilv=asoc)


#######################################################################

#Now let us assume that we have decided to have time periods of length 10
length_of_time_period<-10
#We need to convert the ta1 into the time period in which each individual learned.
#We can do this by using the integer divide %/% which divides two numbers but rounds DOWN, and adding 1, since anyone learned in time 0 to 10
#learned in time period 1, not time period 0.
tp1<-ta1%/%length_of_time_period+1
tp1

#The endTime is now the final time period, since this is a complete diffusion we can set this to the period in which the last individual learned: 57

#We then use the dTADAData function to create the dTADAData object
nbdaData_dTADA<-dTADAData(label="ExampleDiffusion_dTADA",assMatrix=socNet1,orderAcq = oa1,timeAcq = tp1, endTime = max(tp1),asoc_ilv = asoc, int_ilv=asoc)

#We can see this is a different class of object to an nbdaData object
class(nbdaData_dTADA)
#[1] "dTADAData"
#attr(,"package")
#[1] "NBDA"
class(nbdaData_cTADA)
#[1] "nbdaData"
#attr(,"package")
#[1] "NBDA"
#So it can only be used to fit dTADA models, but it also means the tadaFit function can detect whether to fit a cTADA or dTADA model:
model_dTADA_constant<-tadaFit(nbdaData_dTADA)
data.frame(Variable=model_dTADA_constant@varNames,MLE=model_dTADA_constant@outputPar,SE=model_dTADA_constant@se)
#                 Variable         MLE          SE
#1         Scale (1/rate): 459.1536420 390.7832751
#2 1 Social transmission 1  17.6903472  15.7640528
#3        2 Asocial: stAge  -0.8104670   0.8734262
#4         3 Social: stAge   0.3262354   0.2195108

#We get the same outputs from a dTADA model which can be interpreted in the same way as for cTADA (see tutorial 4)
#We can see the estimates are different to the cTADA but gives us a similar message.

model_cTADA_constant<-tadaFit(nbdaData_cTADA)
data.frame(Variable=model_cTADA_constant@varNames,MLE=model_cTADA_constant@outputPar,SE=model_cTADA_constant@se)
#Variable          MLE           SE
#1         Scale (1/rate): 4734.4023556 4016.0905054
#2 1 Social transmission 1   16.8826352   14.9552133
#3        2 Asocial: stAge   -0.7520922    0.8683823
#4         3 Social: stAge    0.3035350    0.2151248


#The main difference is in scale, since we have reduced 10 time units into 1, so scale is about 10x smaller here. If required, we can tweak
#the dTADA model so the scale is appropriate to the real units of time, by specifying the length of each time period
#as the vector dTimePeriods; by default, this is set to rep(1,endTime), so we can change this to rep(10,max(tp1))
nbdaData_dTADA2<-dTADAData(label="ExampleDiffusion_cTADA",assMatrix=socNet1,orderAcq = oa1,timeAcq = tp1, endTime = max(tp1),
                           dTimePeriods=rep(10,max(tp1)),asoc_ilv = asoc, int_ilv=asoc)
model_dTADA_constant2<-tadaFit(nbdaData_dTADA2)
data.frame(Variable=model_dTADA_constant2@varNames,MLE=model_dTADA_constant2@outputPar,SE=model_dTADA_constant2@se)
#Variable          MLE           SE
#1         Scale (1/rate): 4591.5342657 3907.8295243
#2 1 Social transmission 1   17.6903369   15.7640387
#3        2 Asocial: stAge   -0.8104671    0.8734260
#4         3 Social: stAge    0.3262362    0.2195108

#Nothing has changed except the scale, which is now similar to the cTADA


#The fits will converge on the same estimates as the cTADA as time periods get shorter--as illustrated if we set the time periods
#to length 0.1. We can do this by multiplying the times by 10, and setting the length of time period to 1

length_of_time_period<-1
tp2<-(ta1*10)%/%length_of_time_period+1
tp2


nbdaData_dTADA3<-dTADAData(label="ExampleDiffusion_cTADA",assMatrix=socNet1,orderAcq = oa1,timeAcq = tp2, endTime = max(tp2),
                              dTimePeriods=rep(0.1,max(tp2)),asoc_ilv = asoc, int_ilv=asoc)
#Warning this model takes a LONG time to fit--you can take the results below on trust if you prefer not to wait!
model_dTADA_constant3<-tadaFit(nbdaData_dTADA3)
data.frame(Variable=model_dTADA_constant3@varNames,MLE=model_dTADA_constant3@outputPar,SE=model_dTADA_constant3@se)

#Variable          MLE           SE
#1         Scale (1/rate): 4686.8958713 3978.4947503
#2 1 Social transmission 1   16.6814100   14.8005503
#3        2 Asocial: stAge   -0.7715070    0.8698590
#4         3 Social: stAge    0.3071391    0.2161572

#You can see that the fit is now very close to the cTADA. However, you probably also noticed it takes a lot longer to fit the model.
#Therefore, it is better to use the cTADA than a dTADA with very small time periods if you have precise and accurate times of acquisition.

#We can fit models with non-constant baselines in the same way as for cTADA, e.g.
model_dTADA_gamma<-tadaFit(nbdaData_dTADA2,baseline="gamma")
data.frame(Variable=model_dTADA_gamma@varNames,MLE=model_dTADA_gamma@outputPar,SE=model_dTADA_gamma@se)
#               Variable        MLE         SE
#1         Scale (1/rate): 26.6724367 15.9161497
#2                   Shape 15.2014027  7.3921090
#3 1 Social transmission 1  0.1406955  0.2430626
#4        2 Asocial: stAge -0.6429907  0.4794006
#5         3 Social: stAge  1.5460867  0.8375252

#compare to results for cTADA:
model_cTADA_gamma<-tadaFit(nbdaData_cTADA,baseline="gamma")
data.frame(Variable=model_cTADA_gamma@varNames,MLE=model_cTADA_gamma@outputPar,SE=model_cTADA_gamma@se)
#              Variable        MLE         SE
#1         Scale (1/rate): 31.8749336 23.3470376
#2                   Shape 13.2441307  7.5801080
#3 1 Social transmission 1  0.2216704  0.3800373
#4        2 Asocial: stAge -0.7376322  0.5477921
#5         3 Social: stAge  1.3152773  0.7394257


model_dTADA_weibull<-tadaFit(nbdaData_dTADA2,baseline="weibull")
data.frame(Variable=model_dTADA_weibull@varNames,MLE=model_dTADA_weibull@outputPar,SE=model_dTADA_weibull@se)
#              Variable         MLE         SE
#1         Scale (1/rate): 469.1442840 84.4030608
#2                   Shape   4.1952996  0.9327357
#3 1 Social transmission 1   0.2727547  0.4564938
#4        2 Asocial: stAge  -0.8275726  0.5776173
#5         3 Social: stAge   1.4412530  0.7349971

#compare to results for cTADA:
model_cTADA_weibull<-tadaFit(nbdaData_cTADA,baseline="weibull")
data.frame(Variable=model_cTADA_weibull@varNames,MLE=model_cTADA_weibull@outputPar,SE=model_cTADA_weibull@se)
#              Variable           MLE           SE
#1         Scale (1/rate):  4.037055e+02 20.158151117
#2                   Shape  4.763533e+00  0.662451787
#3 1 Social transmission 1  8.777549e-04  0.001389945
#4        2 Asocial: stAge -2.598330e-01  0.290078859
#5         3 Social: stAge  4.562159e+00  1.062500396

#Note that we cannot compare AICc between a cTADA and dTADA model, since they are fitted to a different response variable.


#So if we have a static network or networks, and time-constant ILVs, fitting a dTADA is usually just a case of coverting the times of acquisition into
#discrete times (if they are not already in that format), using the dTADAData function to create the data object, then proceeding as usual.

#Minor complications arise if

#1. We have time periods of unequal length-
# this is solved by providing a vector of lengths of each time period to the dTimePeriods argument, in the dTADAData function

#2. We have a dynamic network. As for OADA and cTADA, we provide a 4-dimensional array (see Tutorial 1), and a vector assMatrixIndex. However, here assMatrixIndex
# is of length equal to the number of time periods (endTime), and specifies which time section of the network is appropriate for each time period.

#3. We have time-varying ILVs. We proceed as for OADA and cTADA (see Tutorial 2), except instead of providing values of the ILVs for each acquisition event, we
#provide them for each time period (example below).



#############################################################################
# TUTORIAL 5.2
# ADDING TIME-VARYING ILVs TO A dTADA MODEL
# 1 diffusion
# 1 static network
# 1 Time constant ILV
# 1 Time varying ILV
#############################################################################

#Time-varying ILVs are easily added to a dTADA
#It might take a bit of work to set up the dTADAData object, but once this is done the analysis proceeds as above.

#NOTE the meaning of the structure is a little different to a continuous TADA or an OADA
#If one or more of our ILVs is time-varying, we need to set up an array for ALL of our ILVs specifying their
#values for every individual for every acquistion event.

#For example, let us assume that individuals 1-5 and 16-20 were undergoing some kind of experimental manipulation during
#TIME PERIODS 1-30 and 6-10 and 21-30 were undergoing the same experimental manipulation for TIME PERDIODS 31-57.
#Let us set up a matrix with rows = number of individuals, columns = number of time periods (in this case, 30 x 57)

treatment<-matrix(NA,nrow=30,ncol=57)
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

treatment[,16]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,17]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,18]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,19]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,20]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,21]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,22]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,23]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,24]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,25]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,26]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,27]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,28]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,29]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))
treatment[,30]<-c(rep(1,5),rep(0,10),rep(1,5),rep(0,10))

treatment[,31]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,32]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,33]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,34]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,35]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,36]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,37]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,38]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,39]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,40]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,41]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,42]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,43]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,44]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,45]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))

treatment[,46]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,47]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,48]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,49]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,50]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,51]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,52]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,53]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,54]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,55]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,56]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))
treatment[,57]<-c(rep(0,5),rep(1,10),rep(0,5),rep(1,10))





treatment

#let us also assume we want to include the stAge variable from above
stAge
#We have to input this as a time varying ILV as well, even though it does not change
#However it is easy to create this using the byrow=F argument:
stAgeTV<-matrix(stAge,nrow=30,ncol=57,byrow = F)
stAgeTV

asocTV<-c("treatment","stAgeTV")

#Now in addition, will try the same oa1 and tp1 vectors from above
oa1
tp1
#And we will set the maximum time steps to 60

#We can then create the dTADAData object for the unconstrained model, specifying asocialTreatment="timevarying"
#We will just stick to an additive model here

nbdaData_dTADA_TVadd<-dTADAData(label="ExampleDiffusion2_dTADA",assMatrix=socNet1,orderAcq = oa1,timeAcq = tp1,
                                endTime = 57,asoc_ilv = asocTV,asocialTreatment = "timevarying" )

#Fit the model
model_dTADA_TVadd<-tadaFit(nbdaData_dTADA_TVadd)
#Display the output
data.frame(Variable=model_dTADA_TVadd@varNames,MLE=model_dTADA_TVadd@outputPar,SE=model_dTADA_TVadd@se)



