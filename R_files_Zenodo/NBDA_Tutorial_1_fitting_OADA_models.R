#To load the NBDA package, you first need to install the package "devtools" in the usual way
#Next load it up as follows:
library(devtools)
#Then download and install my NBDA package from GitHub:
devtools::install_github("whoppitt/NBDA")
#And load it as follows
library(NBDA)

#If you cannot get devtools working, you can download the NBDA package from gitHub 
#manually and install it as source code. Got to https://github.com/whoppitt/nbda
#and click on "Clone or download" then "Download ZIP"
#Next install the package as follows specifying the file path to where you saved the zip file
#install.packages("C:/User1/Downloads/NBDA-master.zip", repos=NULL, type="source")


#Set working directory to where the data files are located on your computer


#############################################################################
# TUTORIAL 1.1
# FITTING A BASIC OADA MODEL
# 1 diffusion
# 1 static network
# 0 ILVs
#############################################################################

#Read in the csv file containing the social network, converting it to a matrix
socNet1<-as.matrix(read.csv(file="exampleStaticSocNet.csv"))

#Convert the matrix to a three dimensional array- this is because the NBDA package
#is designed to work with multiple networks as well as a single network
socNet1<-array(socNet1,dim=c(30,30,1))
#The network needs to be arranged such that row N contains the incoming connections for
#individual N.

#Enter a vector giving the order in which individuals learned the target behaviour
#This corresponds to the individuals' positions in the social network matrix
oa1<-c(26,29,30,8,19,21,22,3,14,12,11,1,17,28,5,9,15,7,6,25,4,13,27,18,20,24,23,16,2,10)
#e.g. the first individual to learn here is in the 26th row and column of socNet1

#Create an nbdaData object containing the data we need to fit an OADA model
#The label is simply a string of text you can use to remind yourself what data is 
#stored in this object
#assMatrix stands for association matrix, since these are most commonly used, but this can be any type of social network
nbdaData1<-nbdaData(label="ExampleDiffusion1",assMatrix=socNet1,orderAcq = oa1)

#We can now fit an OADA model using oadaFit, here we store it in an object named model_social
model_social<-oadaFit(nbdaData1)

#The maximum likelihood estimates (MLEs) for the parameter(s) is stored in the @outputPar slot
model_social@outputPar
#The standard errors for the parameter(s) is stored in the @se slot
model_social@se
#But we can get a neat printout of the model fit as follows

data.frame(Variable=model_social@varNames,MLE=model_social@outputPar,SE=model_social@se)
#giving us:

#               Variable      MLE       SE
#1 1 Social transmission 1 5.614655 7.551307

#So we can see that s has been estimated at 5.6, but the SE looks large in comparison at
#7.55. Does this mean there is not good evidence for social transmission? Not necessarily
#as we shall see.

#The value estimated for s might be difficult to interpret, depending on the network
#used. In such cases, we can obtain an estimate of the % of events that occured by social
#transmission as opposed to asocial learning (%ST) as follows:

nbdaPropSolveByST(model=model_social)
#P(Network 1)  P(S offset) 
#0.90783      0.00000 

#This tells us that the estimated value for s corresponds to 90.7% (the function returns a 
#proportion so multiply by 100 to get %ST)
#P(Network 1) stands for the proportion of events that were a result of transmission through 
#network 1 (we only have one network)
#P(S offset) can be ignored for now

#By default the calculation of %ST excludes the first learning event (innovation) for an unseeded
#diffusion, since we know this had to be asocial learning. If we want to include them we can do so:

nbdaPropSolveByST(model=model_social, exclude.innovations = F)
#P(Network 1)  P(S offset) 
#0.87757      0.00000 


#Let us fit an asocial model (with s=0), by specifying type="asocial"
model_asocial<-oadaFit(nbdaData1,type="asocial")

#And then compare the social and asocial models using 

model_social@aicc
model_asocial@aicc
model_asocial@aicc-model_social@aicc

#So the social model is favoured by 5.35 AICc units. This means:
exp(0.5*(model_asocial@aicc-model_social@aicc))
#[1] 14.50248
#the social model is 14.5x more likely to be the best K-L model, out of the two.
#Or we can say the social model has 14.5x more support than the asocial model.

#We can also conduct a likelihood ratio test (LRT) for social transmission
#The @loglik slot contains the -log-likelihood- i.e. minus the log-likelihood
#So we can get the test statistic as double the difference in -log-likelihood as follows:
2*(model_asocial@loglik-model_social@loglik)
#[1] 7.491497

#There is 1 parameter in model_social, and 0 in model_asocial, so we have 1 d.f.
pchisq(2*(model_asocial@loglik-model_social@loglik),df=1,lower.tail=F)
#[1] 0.006199102
#p= 0.0062, strong evidence of an effect consistent with social transmission


#We can get 95% confidence intervals (C.I.s) for the parameter by first plotting the
#profile log-likelihood function. This is the -log-likelihood for a specified
#value of the parameter, when all other parameters in the model have been optimized.
#In this case there are no other parameters, so the profile log-likelihood is the
#same as the -log-likelihood.
#We specify which=1 because we are interested in the first parameter in the model
#as listed in the model output above. We specify the name of the model, and also the range
#we want to plot over, and resolution determines how many points will be plotted
plotProfLik(which=1,model=model_social,range=c(0,10),resolution=20)

#Any values of s for which the profile log-likelihood is above the dotted line would be
#REJECTED in a likelihood ratio test (LRT), and therefore are OUTSIDE the 95% C.I.
#Therefore we can get the 95% C.I. by finding the crossing points

#Here we can see that there is one crossing point between 0 and 2, but it looks like we
#have to zoom out a bit to see the upper limit. We can make the range wider to do this
plotProfLik(which=1,model=model_social,range=c(0,150),resolution=20)

#The upper limit seems to be between 120 and 150.

#Before we move on to find these points, take a moment to note the asymmetry in the
#profile log-likelihood. We have quite a lot of certainty about the lower limit of
#s, but little certainty about the upper value.
#It is this uncertainty about the upper value that led to a high SE above, and thus
#the SE failed to quantify the strength of evidence against the null hypothesis, s=0.

#So we need to find the cross over points to get the 95% C.I. We use the profLikCI
#function, specifying the upperRange and lowerRange to search in.

profLikCI(which=1,model=model_social,upperRange=c(120,150),lowerRange=c(0,2))

#Lower CI   Upper CI 
#0.300734 135.946785 

#s=0 is not included in the 95% C.I. so there is at least reasonable evidence for social transmission
#or at least a statistical effect consistent with social transmission.
#Note we can obtain C.I.s for a different level of confidence by setting, e.g. conf=0.99 in the 
#plotProfLik and profLikCI functions

#We can get an estimate of %ST corresponding to the upper and lower limits of the
#95% C.I. as follows.
#Instead of specifying the model, we specify the parameter values and the name of the nbdaData object

nbdaPropSolveByST(par=0.300734,nbdadata=nbdaData1)
#P(Network 1)  P(S offset) 
#0.47783      0.00000 
nbdaPropSolveByST(par=135.946785,nbdadata=nbdaData1)
#P(Network 1)  P(S offset) 
#0.96294      0.00000

#So between 47.8 - 97.3% of events are estimated to have occurred by social transmission.


#############################################################################
# TUTORIAL 1.2
# ADDING SEEDED DEMONSTRATORS TO AN OADA MODEL
# 1 diffusion
# 1 static network
# Seeded demonstrators
#############################################################################

#Adding seeded demonstrators to a diffusion is straightfoward: simply create a vector showing who started
#the diffusion informed (1) or naive, and input it to the demons argument in nbdaData

#Let us imagine the first 3 individuals to learn were in fact trained and seeded by the
#experimenter

oa1<-c(26,29,30,8,19,21,22,3,14,12,11,1,17,28,5,9,15,7,6,25,4,13,27,18,20,24,23,16,2,10)

#So let create a demons vector full of zeroes
demons<-rep(0,30)
#and slot in 1s for 26,29 and 30
demons[26]<-demons[29]<-demons[30]<-1
demons

#we can then remove them from the order of acquistion:
oa2<-c(8,19,21,22,3,14,12,11,1,17,28,5,9,15,7,6,25,4,13,27,18,20,24,23,16,2,10)

#create the nbdaData object
nbdaData1_seeded<-nbdaData(label="ExampleDiffusion1",assMatrix=socNet1,orderAcq = oa2, demons=demons)
#and fit the model
model_seeded<-oadaFit(nbdaData1_seeded)

#Note that models with seeded demonstrators cannot be compared to models with those same individuals
#included in the order of acqusition using AICc or LRTs, since they are being fitted to different
#data (different orders of acquisition)


#############################################################################
# TUTORIAL 1.3
# ADDING TRANSMISSION WEIGHTS TO AN OADA MODEL
# 1 diffusion
# 1 static network
# Transmission weights
#############################################################################

#Imagine we have times of acquisition (in mins) as follows:
ta1<-c(64,106,114,123,132,217,229,231,331,361,365,462,494,509,724,728,755,799,917,992,1044,1045,1120,1277,1298,1347,1415,1470,1492,1596)

#And that the diffusion ended when the final individual learned (1596 mins)

#We observe the following number of performances by each individual
performances<-c(17,14,13,20,14,13,13,15,24,10,14,11,8,8,11,9,5,15,4,10,8,7,7,1,2,5,4,1,2,0)

#The average rate (per min) at which each individual performed the behaviour once they were informed is calculated thus
transWeight<-performances/(1596-ta1)

#The last individual gets 0/0 = NaN so we can replace this with a zero
transWeight[30]<-0

#We can mutliply by 60, to make the weights per hour, making them more managable numbers:

transWeight<-transWeight*60
transWeight

#Now we can include them by specifying them in the nbdaData function
nbdaData1_weights<-nbdaData(label="ExampleDiffusion1",assMatrix=socNet1,orderAcq = oa1, weights = transWeight)
model_transWeights<-oadaFit(nbdaData1_weights)

data.frame(Variable=model_transWeights@varNames,MLE=model_transWeights@outputPar,SE=model_transWeights@se)
#                Variable     MLE       SE
#1 1 Social transmission 1 4.92103 6.795499

#We can compare weighted and unweighted models using AICc if they are fitted to the same order of acquisition
model_transWeights@aicc
#[1] 146.0437
model_social@aicc
#[1] 143.9678

#In this case the unweighted model is favoured


#############################################################################
# TUTORIAL 1.4
# USING A DYNAMIC NETWORK IN OADA
# 1 diffusion
# 1 dynamic network
#############################################################################

#Imagine we believe the social network changed sometime between the 20th and 21st acquisition event.
#Let us load in the networks for each time period
socNet2a<-as.matrix(read.csv(file="exampleDynamicSocNetA.csv"))
socNet2b<-as.matrix(read.csv(file="exampleDynamicSocNetB.csv"))

#We need to combine these in a 4 dimensional array, with the 4th dimension for time periods

#Create the empty array
socNet2<-array(NA,dim=c(30,30,1,2))
#Slot in the network for each time period
socNet2[,,,1]<-socNet2a
socNet2[,,,2]<-socNet2b

#Now we need to create a vector specifying which time period corresponds to which acquisition event
#Here the first 20 events occured in the first period and the next 10 events occured in the second time period
assMatrixIndex<-c(rep(1,20),rep(2,10))
#Now we enter the 4 dimensional network and assMatrixIndex as follows
nbdaData1_weights<-nbdaData(label="ExampleDiffusionDynNet",assMatrix=socNet2,orderAcq = oa1,assMatrixIndex =assMatrixIndex )
model_dynamic<-oadaFit(nbdaData1_weights)
data.frame(Variable=model_dynamic@varNames,MLE=model_dynamic@outputPar,SE=model_dynamic@se)
#              Variable      MLE       SE
#1 1 Social transmission 1 6.022794 7.911388

#Models with a dynamic network can be compared to static network models if they are fitted to the same order of acquisition
model_dynamic@aicc
#[1] 142.9997
model_social@aicc
#[1] 143.9678

#Here the dynamic model fits slightly better than the static network above

#Note we can combine dynamic networks, transmission weights and seeded demonstrators in the same model if required.


#############################################################################
# TUTORIAL 1.5
# USING THE PRESENCE MATRIX IN OADA
# 1 diffusion
# 1 static network
#############################################################################

#In some cases, we may have individuals entering and/or leaving the population during the course of the diffusion.
#We need to provide this information when constructing the nbdaData object since if a naive individual is absent for
#an event, then it was not available to learn for that event, and if an informed individual is absent for an event, it
#is not available to be learned from for that event.

#We can provide this information using the presenceMatrix argument in the nbdaData function.

#We will use the same diffusion data as for tutorial 1.1:

socNet1<-as.matrix(read.csv(file="exampleStaticSocNet.csv"))
socNet1<-array(socNet1,dim=c(30,30,1))
oa1<-c(26,29,30,8,19,21,22,3,14,12,11,1,17,28,5,9,15,7,6,25,4,13,27,18,20,24,23,16,2,10)

#And then let us add the following information:

#individual 30 left the population between events 9 and 10
#individual 20 entered the population between events 14 and 15
#individual 1 entered the population between events 4 and 5 and left again between events 19 and 20

#The presenceMatrix is (from the help file):
#an optional binary matrix specifying who was present in the diffusion for each event- set to 1s by default. 
#Number of rows to match the number of individuals, and the number of columns to match the number of events 
#(length of orderAcq). 1 denotes that an individual was present in the diffusion for a given event, 
#0 denotes that an individual was absent, and so could neither learn nor transmit the behaviour to others 
#for that event.

#Since
length(oa1)
#is 30 we can start by creating a 30 (individuals) x 30 (events) matrix of 1s:

PresMat1<-matrix(1, nrow=30, ncol=30)

#Then we replace elements with a 0 for the individuals that were missing for specific events:
#individual 30:
PresMat1[30,10:30]<-0
#individual 20:
PresMat1[20,1:14]<-0
#individual 1:
PresMat1[1,1:4]<-0
PresMat1[1,20:30]<-0
PresMat1



#Create an nbdaData object with the presenceMatrix specified
nbdaData1_PresMat1<-nbdaData(label="ExampleDiffusion1_PresMat1",assMatrix=socNet1,orderAcq = oa1, presenceMatrix = PresMat1)

#and fit the OADA model:
model_social_PresMat1<-oadaFit(nbdaData1_PresMat1)

#with the output:
data.frame(Variable=model_social_PresMat1@varNames,MLE=model_social_PresMat1@outputPar,SE=model_social_PresMat1@se)

#which we can see is different to the model which does not include the presence matrix:
data.frame(Variable=model_social@varNames,MLE=model_social@outputPar,SE=model_social@se)
#Though in this case it makes little difference to the estimate of s

#Inference then proceeds as in tutorial 1.1




