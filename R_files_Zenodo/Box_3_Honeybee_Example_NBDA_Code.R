#To load the NBDA package, you first need to install the package "devtools" in the usual way
#Next load it up as follows:
library(devtools)
#Then download and install the NBDA package from GitHub:
devtools::install_github("whoppitt/NBDA")
#And load it as follows
library(NBDA)

#Set working directory to where the data files are located on your computer


#############################################################################
# Box 3
# TESTING FOR SOCIAL TRANSMISSION ACROSS MULTIPLE PATHWAYS
# 1 diffusion
# 3 networks (time-varying and static alternatives)
# 0 ILVs
#############################################################################

#Honeybees gather information from experienced foragers in the hive before leaving to search for the advertised food source in the field
#Here, foragers were searching for a scented artificial feeding station to which a subset of individuals had previously been trained.
#Dances provide spatial information about a foraging site
#Trophallaxis involves receiving nectar from a successful forager; thus, it can transmit scent, taste, and nectar quality information.
#Antennation allows for detecting the scent of a food source on another forager.
#These latter two pathways may thereby allow bees to learn the scent of the target foraging patch and improve the likelihood of successfully finding it
#when searching in the field.
#These network were treated as 'observation networks', in the sense that these interactions advertised the availability of a specific feeding site and
#provided information relevant to successfully locating it.
#To capture the time course of these interactions, we allow the networks to update during the diffusion
#Networks updated when the next forager to discover the feeding station last left the hive before its arrival at the station.
#In one case, a forager discovered the feeder after another bee, but had left before that individual.
#To avoid having the network "rewind" in time, this first individual used the update time of the latter
#This resulted in 15 points at which the networks updated, indicated by Index
#Read in the .csv files containing the social networks, converting them to a matrix

DanceNetwork.Index1<-as.matrix(read.csv(file="DanceNetwork_Index1.csv"))
DanceNetwork.Index2<-as.matrix(read.csv(file="DanceNetwork_Index2.csv"))
DanceNetwork.Index3<-as.matrix(read.csv(file="DanceNetwork_Index3.csv"))
DanceNetwork.Index4<-as.matrix(read.csv(file="DanceNetwork_Index4.csv"))
DanceNetwork.Index5<-as.matrix(read.csv(file="DanceNetwork_Index5.csv"))
DanceNetwork.Index6<-as.matrix(read.csv(file="DanceNetwork_Index6.csv"))
DanceNetwork.Index7<-as.matrix(read.csv(file="DanceNetwork_Index7.csv"))
DanceNetwork.Index8<-as.matrix(read.csv(file="DanceNetwork_Index8.csv"))
DanceNetwork.Index9<-as.matrix(read.csv(file="DanceNetwork_Index9.csv"))
DanceNetwork.Index10<-as.matrix(read.csv(file="DanceNetwork_Index10.csv"))
DanceNetwork.Index11<-as.matrix(read.csv(file="DanceNetwork_Index11.csv"))
DanceNetwork.Index12<-as.matrix(read.csv(file="DanceNetwork_Index12.csv"))
DanceNetwork.Index13<-as.matrix(read.csv(file="DanceNetwork_Index13.csv"))
DanceNetwork.Index14<-as.matrix(read.csv(file="DanceNetwork_Index14.csv"))
DanceNetwork.Index15<-as.matrix(read.csv(file="DanceNetwork_Index15.csv"))

TrophallaxisNetwork.Index1<-as.matrix(read.csv(file="TrophallaxisNetwork_Index1.csv"))
TrophallaxisNetwork.Index2<-as.matrix(read.csv(file="TrophallaxisNetwork_Index2.csv"))
TrophallaxisNetwork.Index3<-as.matrix(read.csv(file="TrophallaxisNetwork_Index3.csv"))
TrophallaxisNetwork.Index4<-as.matrix(read.csv(file="TrophallaxisNetwork_Index4.csv"))
TrophallaxisNetwork.Index5<-as.matrix(read.csv(file="TrophallaxisNetwork_Index5.csv"))
TrophallaxisNetwork.Index6<-as.matrix(read.csv(file="TrophallaxisNetwork_Index6.csv"))
TrophallaxisNetwork.Index7<-as.matrix(read.csv(file="TrophallaxisNetwork_Index7.csv"))
TrophallaxisNetwork.Index8<-as.matrix(read.csv(file="TrophallaxisNetwork_Index8.csv"))
TrophallaxisNetwork.Index9<-as.matrix(read.csv(file="TrophallaxisNetwork_Index9.csv"))
TrophallaxisNetwork.Index10<-as.matrix(read.csv(file="TrophallaxisNetwork_Index10.csv"))
TrophallaxisNetwork.Index11<-as.matrix(read.csv(file="TrophallaxisNetwork_Index11.csv"))
TrophallaxisNetwork.Index12<-as.matrix(read.csv(file="TrophallaxisNetwork_Index12.csv"))
TrophallaxisNetwork.Index13<-as.matrix(read.csv(file="TrophallaxisNetwork_Index13.csv"))
TrophallaxisNetwork.Index14<-as.matrix(read.csv(file="TrophallaxisNetwork_Index14.csv"))
TrophallaxisNetwork.Index15<-as.matrix(read.csv(file="TrophallaxisNetwork_Index15.csv"))

AntennationNetwork.Index1<-as.matrix(read.csv(file="AntennationNetwork_Index1.csv"))
AntennationNetwork.Index2<-as.matrix(read.csv(file="AntennationNetwork_Index2.csv"))
AntennationNetwork.Index3<-as.matrix(read.csv(file="AntennationNetwork_Index3.csv"))
AntennationNetwork.Index4<-as.matrix(read.csv(file="AntennationNetwork_Index4.csv"))
AntennationNetwork.Index5<-as.matrix(read.csv(file="AntennationNetwork_Index5.csv"))
AntennationNetwork.Index6<-as.matrix(read.csv(file="AntennationNetwork_Index6.csv"))
AntennationNetwork.Index7<-as.matrix(read.csv(file="AntennationNetwork_Index7.csv"))
AntennationNetwork.Index8<-as.matrix(read.csv(file="AntennationNetwork_Index8.csv"))
AntennationNetwork.Index9<-as.matrix(read.csv(file="AntennationNetwork_Index9.csv"))
AntennationNetwork.Index10<-as.matrix(read.csv(file="AntennationNetwork_Index10.csv"))
AntennationNetwork.Index11<-as.matrix(read.csv(file="AntennationNetwork_Index11.csv"))
AntennationNetwork.Index12<-as.matrix(read.csv(file="AntennationNetwork_Index12.csv"))
AntennationNetwork.Index13<-as.matrix(read.csv(file="AntennationNetwork_Index13.csv"))
AntennationNetwork.Index14<-as.matrix(read.csv(file="AntennationNetwork_Index14.csv"))
AntennationNetwork.Index15<-as.matrix(read.csv(file="AntennationNetwork_Index15.csv"))

#We need to combine these in a 4-dimensional array
#Dimensions 1 and 2 indicate the number of rows and columns, and are equal to the number of individuals in the diffusion
#Dimension 3 indicates the number of network types
#Dimension 4 indicates the number of time periods
#Create the empty array
assMatrix.TimeVarying<-array(NA,dim=c(35,35,3,15))

#Slot in the network for each time period
assMatrix.TimeVarying[,,1,1]<-DanceNetwork.Index1
assMatrix.TimeVarying[,,2,1]<-TrophallaxisNetwork.Index1
assMatrix.TimeVarying[,,3,1]<-AntennationNetwork.Index1
assMatrix.TimeVarying[,,1,2]<-DanceNetwork.Index2
assMatrix.TimeVarying[,,2,2]<-TrophallaxisNetwork.Index2
assMatrix.TimeVarying[,,3,2]<-AntennationNetwork.Index2
assMatrix.TimeVarying[,,1,3]<-DanceNetwork.Index3
assMatrix.TimeVarying[,,2,3]<-TrophallaxisNetwork.Index3
assMatrix.TimeVarying[,,3,3]<-AntennationNetwork.Index3
assMatrix.TimeVarying[,,1,4]<-DanceNetwork.Index4
assMatrix.TimeVarying[,,2,4]<-TrophallaxisNetwork.Index4
assMatrix.TimeVarying[,,3,4]<-AntennationNetwork.Index4
assMatrix.TimeVarying[,,1,5]<-DanceNetwork.Index5
assMatrix.TimeVarying[,,2,5]<-TrophallaxisNetwork.Index5
assMatrix.TimeVarying[,,3,5]<-AntennationNetwork.Index5
assMatrix.TimeVarying[,,1,6]<-DanceNetwork.Index6
assMatrix.TimeVarying[,,2,6]<-TrophallaxisNetwork.Index6
assMatrix.TimeVarying[,,3,6]<-AntennationNetwork.Index6
assMatrix.TimeVarying[,,1,7]<-DanceNetwork.Index7
assMatrix.TimeVarying[,,2,7]<-TrophallaxisNetwork.Index7
assMatrix.TimeVarying[,,3,7]<-AntennationNetwork.Index7
assMatrix.TimeVarying[,,1,8]<-DanceNetwork.Index8
assMatrix.TimeVarying[,,2,8]<-TrophallaxisNetwork.Index8
assMatrix.TimeVarying[,,3,8]<-AntennationNetwork.Index8
assMatrix.TimeVarying[,,1,9]<-DanceNetwork.Index9
assMatrix.TimeVarying[,,2,9]<-TrophallaxisNetwork.Index9
assMatrix.TimeVarying[,,3,9]<-AntennationNetwork.Index9
assMatrix.TimeVarying[,,1,10]<-DanceNetwork.Index10
assMatrix.TimeVarying[,,2,10]<-TrophallaxisNetwork.Index10
assMatrix.TimeVarying[,,3,10]<-AntennationNetwork.Index10
assMatrix.TimeVarying[,,1,11]<-DanceNetwork.Index11
assMatrix.TimeVarying[,,2,11]<-TrophallaxisNetwork.Index11
assMatrix.TimeVarying[,,3,11]<-AntennationNetwork.Index11
assMatrix.TimeVarying[,,1,12]<-DanceNetwork.Index12
assMatrix.TimeVarying[,,2,12]<-TrophallaxisNetwork.Index12
assMatrix.TimeVarying[,,3,12]<-AntennationNetwork.Index12
assMatrix.TimeVarying[,,1,13]<-DanceNetwork.Index13
assMatrix.TimeVarying[,,2,13]<-TrophallaxisNetwork.Index13
assMatrix.TimeVarying[,,3,13]<-AntennationNetwork.Index13
assMatrix.TimeVarying[,,1,14]<-DanceNetwork.Index14
assMatrix.TimeVarying[,,2,14]<-TrophallaxisNetwork.Index14
assMatrix.TimeVarying[,,3,14]<-AntennationNetwork.Index14
assMatrix.TimeVarying[,,1,15]<-DanceNetwork.Index15
assMatrix.TimeVarying[,,2,15]<-TrophallaxisNetwork.Index15
assMatrix.TimeVarying[,,3,15]<-AntennationNetwork.Index15

#Now we need to create a vector specifying which time period corresponds to which acquisition event
#Here, each event occurred in its own time period except for events 9 and 10 that occurred during the same period (Index 9)
assMatrixIndex<-c(1,2,3,4,5,6,7,8,9,9,10,11,12,13,14,15)

#Enter a vector giving the order in which honeybee recruits arrived at the target feeding station
#These values correspond to individuals' position (row and column) in the social network matrix
order.Recruitment<-c(18,14,29,22,13,32,25,35,31,20,33,17,34,15,24,21)

#In this experiment, some bees had previously been trained to the target feeding site in order to recruit others during the trial
#We indicate which individuals in the diffusion were these previously trained individuals with the demons argument
#A value of 1 indicates a trained, informed demonstrator, while a value of 0 indicates a naive, potential recruit
#The order in the demons statement corresponds to an individual's position in the social network
#Thus, individuals 1, 2, 4-12, and 27 were our trained demonstrators
demonBees<-c(1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0)

#create the nbdaData object
label<-"Honeybee Recruitment Example"
HoneybeeRecruitmentOADA<-nbdaData(label=label, assMatrix=assMatrix.TimeVarying, orderAcq=order.Recruitment, assMatrixIndex=assMatrixIndex, demons=demonBees)

#Fit a model in which only the dance following network is included
#First, set up a new NBDA object that constrains s=0 for both the trophallaxis and antennation networks
RecruitOADA1<-constrainedNBDAdata(HoneybeeRecruitmentOADA, constraintsVect=c(1,0,0))
#Fit the model
M1<-oadaFit(RecruitOADA1)

#Fit a model in which s is estimated separately for all three networks
RecruitOADA2<-constrainedNBDAdata(HoneybeeRecruitmentOADA, constraintsVect=c(1,2,3))
#Fit the model
M2<-oadaFit(RecruitOADA2)

#Fit a model in which s is constrained to be equal across all three networks
RecruitOADA3<-constrainedNBDAdata(HoneybeeRecruitmentOADA, constraintsVect=c(1,1,1))
#Fit the model
M3<-oadaFit(RecruitOADA3)

#Fit an asocial model in which s is constrained to 0 for all three networks
M4<-oadaFit(HoneybeeRecruitmentOADA, type="asocial")

#We can also construct the static equivalent of our dynamic observation networks
#Each network includes all interactions across the complete trial between informed foragers and individuals that have yet to discover the feeding station
#Read in the csv files containing the social networks, converting them to a matrix
DanceNetwork.Static<-as.matrix(read.csv(file="DanceNetwork_Static.csv"))
TrophallaxisNetwork.Static<-as.matrix(read.csv(file="TrophallaxisNetwork_Static.csv"))
AntennationNetwork.Static<-as.matrix(read.csv(file="AntennationNetwork_Static.csv"))

#Combine these in a 3-dimensional array
#The fourth dimension for time period is no longer needed, as all networks in the model are static
assMatrix.Static<-array(NA,dim=c(35,35,3))
assMatrix.Static[,,1]<-DanceNetwork.Static
assMatrix.Static[,,2]<-TrophallaxisNetwork.Static
assMatrix.Static[,,3]<-AntennationNetwork.Static

order.Recruitment<-c(18,14,29,22,13,32,25,35,31,20,33,17,34,15,24,21)

demonBees<-c(1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0)

#create the NBDA object
label<-"Honeybee Recruitment Example_Static Network"
HoneybeeRecruitmentOADA.StaticNetwork<-nbdaData(label=label, assMatrix=assMatrix.Static, orderAcq=order.Recruitment, demons=demonBees)

#Fit a model that includes the static dance following network only; s=0 for both trophallaxis and antennation networks
RecruitOADA.Static<-constrainedNBDAdata(HoneybeeRecruitmentOADA.StaticNetwork, constraintsVect=c(1,0,0))
#Fit the model
M5<-oadaFit(RecruitOADA.Static)

#Code to create Box 3 Table 1
#Model number corresponds to the model labels used above
model.number<-c(1,2,3,4,5)
#Model type refers to either: a social model using time-varying (TV) networks, a social model using static networks, or an asocial model
model.type<-c("TV Net", "TV Net", "TV Net", "Asocial", "Static Net")
#Obtain the -log-likelihood from each model
model.log<-c(M1@loglik, M2@loglik, M3@loglik, M4@loglik, M5@loglik)
#The number of variables refers to the number of fitted parameters
#For example, in M1 only 1 parameter is fitted: the social transmission rate for the dance following network
#In M2, 3 parameters are fitted: a separate transmission rate, s, for all three networks
number.variables<-c(1,3,1,0,1)
#Obtain the AICc for each model
aicc<-c(M1@aicc,M2@aicc,M3@aicc,M4@aicc,M5@aicc)
#Calculate deltaAICc for each model (i.e. subtract its AICc from the lowest AICc in the model set)
DeltaAicc<-c(M1@aicc-M1@aicc,M2@aicc-M1@aicc,M3@aicc-M1@aicc,M4@aicc-M1@aicc,M5@aicc-M1@aicc)
#Obtain model likelihoods
model.likelihood<-c(exp((-0.5)*DeltaAicc))
#Transform model likelihoods into Akaike weights
AICcWeight<-c(model.likelihood/sum(model.likelihood))

data.frame(Model=model.number, Model.Type=model.type, Log.likelihood=model.log, K=number.variables, AICc=aicc, deltaAICc=DeltaAicc, Akaike.Weight=AICcWeight)

#giving us:

#  Model Model.Type -Log.likelihood K     AICc deltaAICc Akaike.Weight
#1     1     TV Net       30.11308 1 62.51188  0.000000  9.422541e-01
#2     2     TV Net       30.11308 3 68.22616  5.714286  5.411612e-02
#3     3     TV Net       35.67486 1 73.63544 11.123562  3.620073e-03
#4     4    Asocial       43.08151 0 86.16303 23.651151  6.892643e-06
#5     5 Static Net       42.85108 1 87.98787 25.475990  2.767752e-06

#We can compare support for fitting a separate s parameter to all three networks versus constraining s to be equal across them
#This comparison involves models 2 and 3
#As model 3 is nested within model 2, we can employ a likelihood ratio test
#Test statistic
2*(M3@loglik-M2@loglik)
#[1] 11.12356

#There are 3 parameters in M2, and 1 in M3, so we have 2 d.f.
pchisq(2*(M3@loglik-M2@loglik),df=2,lower.tail=F)
#[1] 0.003841928
#p= 0.004, strong evidence in favor of separate social transmission rates across these networks

#We can also compare these models on the basis of AIC
M3@aicc-M2@aicc
#M2 is favoured by 5.41 AICc units. This means:
exp(0.5*(M3@aicc-M2@aicc))
#[1] 14.94891
#Fitting a separate transmission rate to each network has 14.9x more support than constraining transmission rates to be equal across them

M2@outputPar
#In Model 2, the rate of transmission is estimated as essentially zero in the trophallaxis and antennation networks
#The dance following network has a huge social transmission rate
#Indeed, the best supported model (M1) from the five considered here includes only the dance following network
#This is shown by it having the lowest AICc value and an Akaike weight of 0.94
#This latter value indicates that there is little uncertainty as to the best model in this particular model set

#The time-varying dance following network was strongly favoured over a static equivalent
#This can be seen by comparing AICc or Akaike weights
M5@aicc-M1@aicc
exp(0.5*(M5@aicc-M1@aicc))
#[1] 340440.3
#The model that includes the dynamic dance network received 340,440x as much support as the model that included the static network!

#Similarly, M1 is strongly supported (136,704x as much) over the asocial model (M4) that assumed that bees discovered the feeding station entirely through independent search
M4@aicc-M1@aicc
exp(0.5*(M4@aicc-M1@aicc))
#[1] 136704.3

#The estimated social transmission effect is very large for our best model
M1@outputPar
#[1] 99443294

#This is probably due to the order of recruitment following the dance following network extremely closely; see Section 6.4 of the main text
#Indeed, we may not be able to get an upper bound to this effect
plotProfLik(which=1,model=M1,range=c(0,10),resolution=20)
plotProfLik(which=1,model=M1,range=c(0,100),resolution=20)
plotProfLik(which=1,model=M1,range=c(0,1000),resolution=20)
#We can see it appears to level out as s tend to infinity

#We can at least get the lower bound
profLikCI(which=1,model=M1,lowerRange=c(0,2))

# Lower CI  Upper CI 
#0.4450421        NA 

#Let's convert these values into %ST, i.e. the percentage of recruitment events that occurred due to social transmission (via the waggle dance)
nbdaPropSolveByST(par=M1@outputPar,nbdadata = RecruitOADA1)

#P(Network 1)  P(S offset) 
#           1            0 

#So an estimated 100% of the 16 recruitments occurred through social transmission

#The lower bound
nbdaPropSolveByST(par=c(0.4450421), nbdadata=RecruitOADA1)

#P(Network 1)  P(S offset) 
#     0.91213      0.00000 

#This tells us that at least 91% of the recruitment events occurred through dance-mediate social transmission
