#Load the NBDA package as follows
library(NBDA)

#Set working directory to where the data files are located on your computer



#############################################################################
# TUTORIAL 7.1
# MULTI-MODEL INFERENCE IN AN OADA: UNCONSTRAINED MODEL
# 1 diffusion
# 1 static network
# 2 time constant ILVs
#############################################################################


#Read in the social network and order of acquisition vector as shown in Tutorials 1 and 2
socNet1<-as.matrix(read.csv(file="exampleStaticSocNet.csv"))
socNet1<-array(socNet1,dim=c(30,30,1))
oa1<-c(26,29,30,8,19,21,22,3,14,12,11,1,17,28,5,9,15,7,6,25,4,13,27,18,20,24,23,16,2,10)

ILVdata<-read.csv(file="exampleTimeConstantILVs.csv")
#Then extract each ILV
female<-cbind(ILVdata$female)
male<-1-female
#females=0 males=1
age<-cbind(ILVdata$age)
#age in years
#I am going to standardize age
stAge<-(age-mean(age))/sd(age)

asoc<-c("male","stAge")

#Create an object for the "unconstrained" model
nbdaData2_uc<-nbdaData(label="ExampleDiffusion2",assMatrix=socNet1,orderAcq = oa1,asoc_ilv = asoc,int_ilv=asoc)
#Fit and display the model
model2_uc<-oadaFit(nbdaData2_uc)
data.frame(Variable=model2_uc@varNames,MLE=model2_uc@outputPar,SE=model2_uc@se)

#We are not sure which of these variables to include in our model, and decide to use a multi-model inference approach.
#The first stage is to define the set of models we are interested in. Recall that we can specify a model simpler than the one determined by the nbdaData
#object by using constrainedNBDAdata, to cut down the data, and then fit the model to the constrained data. Inside the constrainedNBDAdata function, we
#specify a constraintsVect giving the way in which the data/model is constrained.
#Any parameters with a 0 are constrained to be 0 (omitted from the model), and any parameters with the same non-zero number are constrained to be the same,
#with the restriction that only parameters of the same type can be constrained to be the same: i.e. only 2 or more s parameters, 2 or more effects on asocial
#learning; 2 or more effects on social learning; or 2 or more multiplicative effects.
#e.g. here constraintsVect=c(1,2,0,0,3), would be a social model with an effect of sex on asocial learning and an effect of age on social learning.

#To specify a set of models, we construct a constraintsVectMatrix, with each row giving a constraintsVect and specifying a different model in the set.

constraintsVectMatrix<-rbind(
  #The full model
  c(1,2,3,4,5),
  #Social models with only ILVs affecting asocial learning
  c(1,2,3,0,0),
  c(1,2,0,0,0),
  c(1,0,2,0,0),
  #Social models with only ILVs affecting social learning
  c(1,0,0,2,3),
  c(1,0,0,2,0),
  c(1,0,0,0,2),
  #Social models with ILVs affecting both asocial and social learning
  c(1,2,0,3,0),
  c(1,2,0,0,3),
  c(1,2,0,3,4),
  c(1,0,2,3,0),
  c(1,0,2,0,3),
  c(1,0,2,3,4),
  c(1,2,3,4,0),
  c(1,2,3,0,4),
  #Social model with no ILVs
  c(1,0,0,0,0),
  #Asocial models with all combinations of ILVs
  c(0,0,0,0,0),
  c(0,1,0,0,0),
  c(0,0,1,0,0),
  c(0,1,2,0,0)
  #Remember asocial models cannot have an effect of ILVs on social learning.
)

#You may find it easier to set out the models in a more systematic way to ensure you do not miss any, as follows:
constraintsVectMatrix<-rbind(
  #social models
  c(1,0,0,0,0),
  c(1,0,0,0,2),
  c(1,0,0,2,0),
  c(1,0,0,2,3),
  c(1,0,2,0,0),
  c(1,0,2,0,3),
  c(1,0,2,3,0),
  c(1,0,2,3,4),
  c(1,2,0,0,0),
  c(1,2,0,0,3),
  c(1,2,0,3,0),
  c(1,2,0,3,4),
  c(1,2,3,0,0),
  c(1,2,3,0,4),
  c(1,2,3,4,0),
  c(1,2,3,4,5),
  #asocial models
  c(0,0,0,0,0),
  c(0,1,0,0,0),
  c(0,1,2,0,0),
  c(0,0,1,0,0))


#or you could write code to generate all combinations, or create the table in Excel and import it into R

#To fit the set of models we use
modelSet_uc<-oadaAICtable(nbdadata = nbdaData2_uc,constraintsVectMatrix =constraintsVectMatrix)

#We can print out the set of models ordered by AICc as follows:
print(modelSet_uc)
#This shows us the model number (row in the original constraintsVectMatrix)
#Model type- asocial, unconstrained if ILVs are affecting social learning or additive if they are affecting only asocial learning, NA if there are no ILVs
#netCombo we will come to below
#baseline is NA because this is an OADA
#The next columns beginning with CONS are the constraints vector defining that model
#The next columns beginning with OFF are the offset vector used to define the model (here not used)
#convergence says whether the optimisation algorithm "thinks" it has converged
#loglik the -log-likelihood for the model
#next are the parameter estimates; these are set to zero if the parameter is constrained to zero
#next are the SEs
#then AIC and AICc; AICc is used by default
#deltaAICc is the difference in AICc from the best model, the relative support for the model compared to the best model, and the Akaike weight

#Here we can see we do not have a clearly favoured best model, so it makes sense to use multi-model inferencing

#The first thing we can do is get support for the network
networksSupport(modelSet_uc)
#which tells us that network 1 has 98.5% support compared to 1.5% support for no network (asocial)
#But we can also see that the number of models is not equal (see main text) so this is not a good way to compare support for social and asocial models
typeSupport(modelSet_uc)
#We can also break down support by model type, but again we have unequal model numbers so this is not so useful here

#More useful is the total support for each variable:
variableSupport(modelSet_uc)
#Showing some support (60%) for an effect of sex on asocial learning but little support for other effects of ILVs
#s1 also has high support but is present in more models than it is absent as we saw above, so this could be misleading

#We can get model averaged estimated as follows:
modelAverageEstimates(modelSet_uc)

#If we prefer, we can get model weighted medians instead (see main text):
modelAverageEstimates(modelSet_uc,averageType = "median")

#We can get unconditional SEs:
unconditionalStdErr(modelSet_uc)

#A neat way to tie this up together is:

rbind(
  support=variableSupport(modelSet_uc),
  MAE=modelAverageEstimates(modelSet_uc),
  USE=unconditionalStdErr(modelSet_uc))

#             s1     ASOC:male ASOC:stAge SOCIAL:male SOCIAL:stAge
#support  0.984571  5.962426e-01  0.3492146  0.22726889   0.35722149
#MAE      3.968586 -1.132901e+01 -0.3361343 -0.03219899   0.09769119
#USE     39.012369  7.255577e+07  0.5473529  0.05297392   0.03940915

#We could use the USE to generate Wald CIs, but we saw this can be misleading in an NBDA where profile likelihoods are asymmetrical
#We already know that the profile likelihoods are asymmetrical for the 2 parameters with support here from a previous tutorial.

#One option would be to obtain the profile likelihood 95% CI for every model in which an s parameter is present and check how robust these are to model
#selection uncertainty. However, we have seen in previous tutorials that finding the 95% CIs is an interactive process, with the user finding the points
#at which the profile -log-likelihood crosses the maximum -log-likelihood +1.92 cutoff line.

#That said, it is quite easy to find the lower limit for an s parameter, since we know it must be between 0 and the MLE for s in that model.
#So we can use an automated procedure to find the lower limit for s in each model and the estimate of propST (%ST) corresponding to this. It is the lower
#limit for s that is likely to be of most interest, since this tells us the strength of evidence for social transmission.

#We can obtain the lower limit of an s parameter (determined by which if there are more than 1) for a given confidence level
#(conf) in all models containing that parameter as follows:

lowerLimitsByModel<-multiModelLowerLimits(which=1,aicTable = modelSet_uc,conf=0.95)
lowerLimitsByModel

#   model netCombo    lowerCI  propST deltaAICc akaikeWeight adjAkWeight cumulAdjAkWeight
#1      9        1 0.38994144 0.75010 0.0000000   0.19047275  0.19345760        0.1934576
#2      1        1 0.30073259 0.47783 0.8673898   0.12344728  0.12538179        0.3188394
#3     13        1 0.40209279 0.75866 1.1352833   0.10797152  0.10966351        0.4285029
#4     10        1 0.40866399 0.75008 1.2875510   0.10005636  0.10162432        0.5301272
#5      2        1 0.26908549 0.46200 1.9136852   0.07316129  0.07430779        0.6044350
#6     14        1 0.14529140 0.45126 2.1230682   0.06588923  0.06692176        0.6713568
#7     11        1 0.09396559 0.61701 2.4186965   0.05683549  0.05772614        0.7290829
#8      5        1 0.40372469 0.54502 2.5009361   0.05454582  0.05540059        0.7844835
#9      3        1 0.35226547 0.52724 2.8398628   0.04604309  0.04676462        0.8312481
#10     6        1 0.05757081 0.29719 3.0352408   0.04175790  0.04241228        0.8736604
#11    15        1 0.12429051 0.64278 3.7513244   0.02919061  0.02964804        0.9033084
#12    12        1 0.04061835 0.56999 3.9587446   0.02631494  0.02672731        0.9300358
#13     4        1 0.23745183 0.48312 4.2344673   0.02292609  0.02328536        0.9533211
#14     7        1 0.47633772 0.58845 4.6267354   0.01884300  0.01913828        0.9724594
#15    16        1 0.01909823 0.33129 5.0137907   0.01552752  0.01577085        0.9882302
#16     8        1 0.01176647 0.22252 5.5990505   0.01158816  0.01176976        1.0000000

#This function also returns the propST corresponding to each lower limit.
#deltaAICc is the difference in AICc from the best model in the full set provided to the function
#akaikeWeight is the corresponding Akaike weight from the full model set.
#adjAkWeight is the adjusted Akaike weight conditional on the parameter being in the model (so they sum to one in this reduced set)
#cumulAdjAkWeight is the cumulative adjusted Akaike weight.

#We can see that in the top few models, the lower limit for s and propST is consistently above 0.3 and 45%, but the effect gets plausibly smaller as we move down the table.
#Across all models, the lower 95% CI for s goes as low as 0.011 and propST as low as 0.223 (22.3%).
#But seeing this at least adds to our confidence that we have good evidence for social transmission, since 0 is outside the 95% CI in all models.

#But these lowest values for lowerCI and propST are in the most implausible models, so these are highly conservative lower limits
#for the effect of social transmission.

#We can obtain a model averaged lower-limit for propST as follows:

sum(lowerLimitsByModel$propST*lowerLimitsByModel$adjAkWeight)
#[1] 0.5966962

#############################################################################
# TUTORIAL 7.2
# MULTI-MODEL INFERENCE IN AN OADA: ADDITIVE VERSUS MULTIPLICATIVE VERSUS ASOCIAL MODELS
# 1 diffusion
# 1 static network
# 2 time constant ILVs
#############################################################################

#An alternative is to consider only additive versus multiplicative versus asocial models (See main text and the 
#additional guidance in the supporting information for dicussion of pros and cons).
#To do this we need to build a different starting nbdaData object:

nbdaData2_addVmulti<-nbdaData(label="ExampleDiffusion2",assMatrix=socNet1,orderAcq = oa1,asoc_ilv = asoc,multi_ilv=asoc)
#Here providing asoc as the argument to asoc_ilv and multi_ilv (not int_ilv)

#A full additive model is therefore given by constraintsVect= c(1,2,3,0,0), whereas a full multiplicative model is given by c(1,0,0,2,3)
#We can specify constrainstVectMatrix to give the 3 required model subsets as follows:


constraintsVectMatrix<-rbind(
  #no ILVs
  c(1,0,0,0,0),
  #additive models
  c(1,2,0,0,0),
  c(1,2,3,0,0),
  c(1,0,2,0,0),
  #multipicative models
  c(1,0,0,2,0),
  c(1,0,0,2,3),
  c(1,0,0,0,2),
  #asocial models
  c(0,0,0,0,0),
  c(0,1,0,0,0),
  c(0,1,2,0,0),
  c(0,0,1,0,0))


modelSet_addVmulti<-oadaAICtable(nbdadata = nbdaData2_addVmulti,constraintsVectMatrix =constraintsVectMatrix)

#We can print out the set of models ordered by AICc as follows:
print(modelSet_addVmulti)

networksSupport(modelSet_addVmulti)
#which tells us that network 1 has 97.5% support compared to 2.5% support for no network (asocial)
#But we can also see that the number of models is not equal (see main text)


typeSupport(modelSet_addVmulti)
#We can also break down support by model type. This is much more useful in this case, since we can do a fair three-way model comparison.
#The noILVs model can be added to both additive and multiplicative sets, since they reduce to the same model when there are no ILVs
#Giving us a fair 4 model comparison
#additive
0.56151439+0.19637216
#multiplicative
0.21757004+0.19637216

#Both are much better supported than asocial learning, with additive being favoured slightly over multiplicative

#We can get the total support for each variable
variableSupport(modelSet_addVmulti)

#or we can condition on a particular model type:
variableSupport(modelSet_addVmulti,typeFilter = "additive")
#               s1 ASOC:male ASOC:stAge A&S:male A&S:stAge
#support 0.9581212 0.8208817  0.4517089        0         0

#This gives the support across all models consistent with the additive assumption--i.e. additive, noILVs and asocial models.
#This has equal numbers of models with each parameter present and absent (assuming we included all combinations above), so it is a fair measure of support for each

#We can get model averaged estimates as follows. Here it probably also makes sense to condition on model type
modelAverageEstimates(modelSet_addVmulti,typeFilter = "additive")

#We can get unconditional SEs:
unconditionalStdErr(modelSet_addVmulti,typeFilter = "additive")

#A neat way to tie this up together is:

rbind(
  support=variableSupport(modelSet_addVmulti,typeFilter = "additive"),
  MAE=modelAverageEstimates(modelSet_addVmulti,typeFilter = "additive"),
  USE=unconditionalStdErr(modelSet_addVmulti,typeFilter = "additive"))

#              s1     ASOC:male ASOC:stAge A&S:male A&S:stAge
#support  0.9581212  8.208817e-01  0.4517089        0         0
#MAE      3.7834965 -1.192426e+01 -0.2960734        0         0
#USE     34.4544549  8.249478e+07  0.5035605        0         0


#If the multiplicative model had been favoured we would have obtained an output like this:
rbind(
  support=variableSupport(modelSet_addVmulti,typeFilter = "multiplicative"),
  MAE=modelAverageEstimates(modelSet_addVmulti,typeFilter = "multiplicative"),
  USE=unconditionalStdErr(modelSet_addVmulti,typeFilter = "multiplicative"))

#                s1   ASOC:male ASOC:stAge    A&S:male  A&S:stAge
#support  0.8986285  0.02617684 0.02562927  0.55018591 0.49376378
#MAE      5.7268541 -0.11157355 0.04073732 -0.11006709 0.04108839
#USE     64.5971940  0.07622397 0.01619253  0.07421112 0.01558622


#The MAEs and USEs in the ASOC (asocial) columns includes the asocial models whereas those in the A&S (asocial and social) columns do not.
#So if we want MAEs including the asocial models we use the ASOC columns, whereas if we want MAEs not including the asocial models we use the A&S columns

#The support in the ASOC columns is for models containing ONLY an effect of the ILV on asocial learning
#whereas the support in the A&S columns is for models containing the same effect on both social and asocial learning
#The total support for the variable, conditional on the multiplicative model can be obtained by adding these two columns together.


#We can examine the plausible lower limit for s across models as before

lowerLimitsByModel_addVmulti<-multiModelLowerLimits(which=1,aicTable = modelSet_addVmulti,conf=0.95)
lowerLimitsByModel_addVmulti

#  model netCombo   lowerCI  propST deltaAICc akaikeWeight adjAkWeight cumulAdjAkWeight
#1     2        1 0.3899414 0.75010 0.0000000   0.30299205  0.31061561        0.3106156
#2     1        1 0.3007326 0.47783 0.8673898   0.19637216  0.20131307        0.5119287
#3     3        1 0.4020928 0.75866 1.1352833   0.17175429  0.17607579        0.6880045
#4     5        1 0.3960483 0.53915 2.2570050   0.09802319  0.10048955        0.7884940
#5     4        1 0.4037247 0.54502 2.5009361   0.08676805  0.08895122        0.8774452
#6     7        1 0.3878685 0.53454 2.5571639   0.08436264  0.08648528        0.9639305
#7     6        1 0.4667553 0.57501 4.3062178   0.03518422  0.03606949        1.0000000

#We can obtain a model averaged lower-limit for propST as follows:

sum(lowerLimitsByModel_addVmulti$propST*lowerLimitsByModel_addVmulti$adjAkWeight)

#############################################################################
# TUTORIAL 7.3
# MULTI-MODEL INFERENCE IN A TADA: UNCONSTRAINED MODEL
# 1 diffusion
# 1 static network
# 2 time constant ILVs
#############################################################################

#Read in the social network and order of acquisition vector as shown in Tutorials 1 and 2
socNet1<-as.matrix(read.csv(file="exampleStaticSocNet.csv"))
socNet1<-array(socNet1,dim=c(30,30,1))
oa1<-c(26,29,30,8,19,21,22,3,14,12,11,1,17,28,5,9,15,7,6,25,4,13,27,18,20,24,23,16,2,10)

ILVdata<-read.csv(file="exampleTimeConstantILVs.csv")
#Then extract each ILV
female<-cbind(ILVdata$female)
male<-1-female
#females=0 males=1
age<-cbind(ILVdata$age)
#age in years
#I am going to standardize age
stAge<-(age-mean(age))/sd(age)

asoc<-c("male","stAge")

#Now in addition, we need a vector giving the times at which each individual acquired the behaviour:
ta1<-c(234,252,262,266,273,296,298,310,313,326,332,334,334,337,338,340,343,367,374,376,377,393,402,405,407,407,435,472,499,567)

#Create an object for the "unconstrained" model
nbdaData_tada<-nbdaData(label="ExampleDiffusion2",assMatrix=socNet1,orderAcq = oa1,asoc_ilv = asoc,int_ilv=asoc,timeAcq = ta1,endTime = 568)
#Fit and display the model (constant baseline)
model_tada<-tadaFit(nbdaData_tada)
data.frame(Variable=model_tada@varNames,MLE=model_tada@outputPar,SE=model_tada@se)

#                 Variable           MLE  SE
#1         Scale (1/rate):  2616.1989966 NaN
#2 1 Social transmission 1    10.3941207 NaN
#3         2 Asocial: male -6352.3570347 NaN
#4        3 Asocial: stAge    -1.1416307 NaN
#5          4 Social: male    -0.1782846 NaN
#6         5 Social: stAge     0.2767679 NaN


#Now for the model set, let us consider all possible sets of ILVs in an unconstrained and an asocial model:

constraintsVectMatrix<-rbind(
  #social models
  c(1,0,0,0,0),
  c(1,0,0,0,2),
  c(1,0,0,2,0),
  c(1,0,0,2,3),
  c(1,0,2,0,0),
  c(1,0,2,0,3),
  c(1,0,2,3,0),
  c(1,0,2,3,4),
  c(1,2,0,0,0),
  c(1,2,0,0,3),
  c(1,2,0,3,0),
  c(1,2,0,3,4),
  c(1,2,3,0,0),
  c(1,2,3,0,4),
  c(1,2,3,4,0),
  c(1,2,3,4,5),
  #asocial models
  c(0,0,0,0,0),
  c(0,1,0,0,0),
  c(0,1,2,0,0),
  c(0,0,1,0,0))

#However, in a TADA we likely want to consider different baseline functions too. So we want to consider all possible sets of ILVs in an
#unconstrained and an asocial model for the constant, weibull and gamma baseline functions
#So we want to replicate the constraintsVectMartix above 3 times:
constraintsVectMatrixAll<-rbind(constraintsVectMatrix,constraintsVectMatrix,constraintsVectMatrix)

#And we also need to provide a vector saying what baseline function each model will have:

baselineVect<-rep(c("constant","weibull","gamma"),each=dim(constraintsVectMatrix)[1])

#Check:
cbind(baselineVect,constraintsVectMatrixAll)

#To fit the set of models, we use
modelSet_tada<-tadaAICtable(nbdadata = nbdaData_tada,constraintsVectMatrix =constraintsVectMatrixAll,baselineVect = baselineVect)
#Note that progress slows as we get to the gamma models since these take longer to fit

#Aside: for both oadaAICtable and tadaAICtable, we can use multiple cores on the computer to speed things up. Most modern computers have at least 4 cores, and 8 "virtual cores".
#So you can probably use 6 and still have computing power left to do other things on your computer while you wait; just specify cores=6. If you have more cores, you can specify more.
#You do not get a progress bar when using multiple cores
modelSet_tada<-tadaAICtable(nbdadata = nbdaData_tada,constraintsVectMatrix =constraintsVectMatrixAll,baselineVect = baselineVect,cores=6)

print(modelSet_tada)
#We can see that the table now also specifies a baseline function for each model, and values for the baseline parameters too
#We can get the support for each baseline as follows:
baselineSupport(modelSet_tada)
#             support numberOfModels
#constant 0.3370271             20
#gamma    0.5212694             20
#weibull  0.1417035             20

#We can then get support for each combination of networks across all models
networksSupport(modelSet_tada)
#      support numberOfModels
#0 0.207637             12
#1 0.792363             48

#But since we are using the unconstrained model, there are more social than asocial models, making this of limited use
#We can get MAE and USEs as follows

rbind(
support=variableSupport(modelSet_tada),
MAE=modelAverageEstimates(modelSet_tada),
USE=unconditionalStdErr(modelSet_tada))

#            s1    ASOC:male ASOC:stAge SOCIAL:male SOCIAL:stAge
#support 0.792363    0.4467811  0.3138500 0.184076914    0.3825203
#MAE     4.621754 -953.5084687 -0.2286202 0.002706989    0.5644324
#USE          NaN          NaN        NaN         NaN          NaN

#Support for the ILVs is also useful here because there are equal numbers of models with and without each ILV
#Notice that the USE are all NaN; this is because SEs could not be derived for all models. A practical solution is to replace these with the
#model-weighted average SE across all other models, and recalculate the USE:

rbind(
  support=variableSupport(modelSet_tada),
  MAE=modelAverageEstimates(modelSet_tada),
  USE=unconditionalStdErr(modelSet_tada,nanReplace = T))

#             s1     ASOC:male ASOC:stAge SOCIAL:male SOCIAL:stAge
#support  0.792363  4.467811e-01  0.3138500 0.184076914    0.3825203
#MAE      4.621754 -9.535085e+02 -0.2286202 0.002706989    0.5644324
#USE     76.480108  3.780749e+06  0.2776675 0.649489562    4.7175325

#However, we recommend that these figures only be used if a few models with low Akaike weight are missing SEs. If we look at the SEs for the individual models:
print(modelSet_tada)
#we can see that a lot are missing in this case, especially for s.

#We can also obtain these results conditional on specific baseline functions, to check whether the choice of baseline influences our conclusions

rbind(
  support=variableSupport(modelSet_tada,baselineFilter = "constant"),
  MAE=modelAverageEstimates(modelSet_tada,baselineFilter = "constant"),
  USE=unconditionalStdErr(modelSet_tada,baselineFilter = "constant"))

#             s1     ASOC:male ASOC:stAge SOCIAL:male SOCIAL:stAge
#support  1.0000     0.5527769  0.2982583  0.24843070    0.3750632
#MAE     11.2341 -2494.9170202 -0.2742417 -0.06992442    0.1052773
#USE         NaN           NaN        NaN         NaN          NaN

rbind(
  support=variableSupport(modelSet_tada,baselineFilter = "gamma"),
  MAE=modelAverageEstimates(modelSet_tada,baselineFilter = "gamma"),
  USE=unconditionalStdErr(modelSet_tada,baselineFilter = "gamma"))

#              s1    ASOC:male ASOC:stAge SOCIAL:male SOCIAL:stAge
#support 0.6240646    0.3421164  0.3232841  0.13227921    0.3979534
#MAE     0.6746116 -170.2258982 -0.1849467  0.06348583    0.9245151
#USE           NaN          NaN        NaN         NaN          NaN

rbind(
  support=variableSupport(modelSet_tada,baselineFilter = "weibull"),
  MAE=modelAverageEstimates(modelSet_tada,baselineFilter = "weibull"),
  USE=unconditionalStdErr(modelSet_tada,baselineFilter = "weibull"))

#            s1    ASOC:male ASOC:stAge SOCIAL:male SOCIAL:stAge
#support 0.917621    0.5797001  0.3162287   0.2215607    0.3434843
#MAE     3.414911 -168.8053590 -0.2807712  -0.0481273    0.3318878
#USE          NaN          NaN        NaN         NaN          NaN

#The ILVs seem robust to choice of baseline function, but the support for s changes a lot. While the support for s is difficult to interpret in the
#unconstrained model (see above), the fact that support for s, and the MAE for s varies so much depending on baseline function suggests that our conclusions
#about social transmission are unlikely to be robust to choice of baseline function.

#We can look into this further by fitting different model sets for each baseline function:

modelSet_constant<-tadaAICtable(nbdadata = nbdaData_tada,constraintsVectMatrix =constraintsVectMatrix,baselineVect = rep("constant",dim(constraintsVectMatrix)[1]),cores=6)
modelSet_weibull<-tadaAICtable(nbdadata = nbdaData_tada,constraintsVectMatrix =constraintsVectMatrix,baselineVect = rep("weibull",dim(constraintsVectMatrix)[1]),cores=6)
modelSet_gamma<-tadaAICtable(nbdadata = nbdaData_tada,constraintsVectMatrix =constraintsVectMatrix,baselineVect = rep("gamma",dim(constraintsVectMatrix)[1]),cores=6)

lowerLimitsByModel_constant<-multiModelLowerLimits(which=1,aicTable = modelSet_constant,conf=0.95)
lowerLimitsByModel_constant

#   model netCombo  lowerCI  propST deltaAICc akaikeWeight adjAkWeight cumulAdjAkWeight
#1      9        1 1.902186 0.88778 0.0000000   0.18039315  0.18039315        0.1803931
#2      1        1 3.934956 0.88586 0.5274254   0.13857697  0.13857697        0.3189701
#3     10        1 1.959740 0.89041 0.9975866   0.10954609  0.10954609        0.4285162
#4      2        1 3.998064 0.88533 1.2956095   0.09438044  0.09438044        0.5228966
#5     13        1 2.096165 0.89415 1.5402517   0.08351388  0.08351388        0.6064105
#6     11        1 2.007292 0.89102 2.1435039   0.06176808  0.06176808        0.6681786
#7      3        1 4.413323 0.89317 2.2568543   0.05836471  0.05836471        0.7265433
#8      5        1 4.197446 0.88987 2.3641868   0.05531506  0.05531506        0.7818584
#9     14        1 2.148816 0.89240 2.5143788   0.05131326  0.05131326        0.8331716
#10     6        1 4.212441 0.88056 3.0873795   0.03853045  0.03853045        0.8717021
#11     4        1 4.187817 0.89167 3.5951715   0.02989086  0.02989086        0.9015930
#12    12        1 1.882475 0.88800 3.6680434   0.02882136  0.02882136        0.9304143
#13    15        1 2.253014 0.89797 3.8948059   0.02573201  0.02573201        0.9561463
#14     7        1 4.741689 0.89654 4.2754043   0.02127295  0.02127295        0.9774193
#15    16        1 2.065459 0.88914 5.4729671   0.01168908  0.01168908        0.9891084
#16     8        1 4.430636 0.88660 5.6142865   0.01089164  0.01089164        1.0000000

lowerLimitsByModel_weibull<-multiModelLowerLimits(which=1,aicTable = modelSet_weibull,conf=0.95)
lowerLimitsByModel_weibull

#   model netCombo      lowerCI  propST deltaAICc akaikeWeight adjAkWeight cumulAdjAkWeight
#1      9        1 3.574105e-01 0.74106 0.0000000  0.173878178 0.189488003        0.1894880
#2      1        1 2.560752e-01 0.44195 0.5142289  0.134456465 0.146527226        0.3360152
#3     10        1 3.324885e-01 0.73927 0.8028994  0.116385180 0.126833601        0.4628488
#4     13        1 3.886976e-01 0.75474 1.6705158  0.075421850 0.082192808        0.5450416
#5     14        1 3.310529e-01 0.74626 1.9146074  0.066756470 0.072749498        0.6177911
#6      2        1 1.362141e-02 0.15350 2.0695996  0.061778482 0.067324614        0.6851158
#7      3        1 7.863683e-01 0.63223 2.3111190  0.054750990 0.059666232        0.7447820
#8     11        1 2.730730e-01 0.71327 2.4739420  0.050470244 0.055001184        0.7997832
#9      5        1 4.425570e-01 0.56332 2.5822441  0.047809908 0.052102018        0.8518852
#11    16        1 1.098705e-01 0.65163 3.6883515  0.027499812 0.029968594        0.8818538
#12    12        1 5.843642e-02 0.59234 3.8764400  0.025031497 0.027278687        0.9091325
#13     8        1 5.779235e-02 0.26650 4.1389313  0.021952681 0.023923471        0.9330559
#14    15        1 3.475977e-01 0.74324 4.3594770  0.019660596 0.021425616        0.9544816
#15     6        1 7.914711e-05 0.13646 4.3683302  0.019573759 0.021330983        0.9758125
#17     7        1 7.253772e-01 0.61514 4.5708758  0.017688541 0.019276521        0.9950891
#20     4        1 0.000000e+00 0.00000 7.3057191  0.004506385 0.004910944        1.0000000


lowerLimitsByModel_gamma<-multiModelLowerLimits(which=1,aicTable = modelSet_gamma,conf=0.95)
lowerLimitsByModel_gamma

#   model netCombo      lowerCI  propST deltaAICc akaikeWeight adjAkWeight cumulAdjAkWeight
#2      2        1 7.572302e-05 0.06099  1.188017  0.117891048 0.188908394        0.1889084
#3      6        1 3.427088e-03 0.15313  1.296982  0.111639878 0.178891531        0.3677999
#5      9        1 0.000000e+00 0.00000  2.020740  0.077742290 0.124574099        0.4923740
#7     10        1 2.481921e-01 0.53180  2.759314  0.053737559 0.086108964        0.5784830
#8      1        1 0.000000e+00 0.00000  3.008375  0.047445510 0.076026596        0.6545096
#9      4        1 5.880874e-05 0.11791  3.441307  0.038210716 0.061228779        0.7157384
#10    14        1 7.376884e-01 0.81218  3.629206  0.034784304 0.055738303        0.7714767
#11    13        1 0.000000e+00 0.00000  3.756859  0.032633511 0.052291877        0.8237685
#12     8        1 6.941501e-05 0.13515  4.392990  0.023742683 0.038045231        0.8618138
#14     3        1 0.000000e+00 0.00000  4.782076  0.019545233 0.031319244        0.8931330
#15    11        1 0.000000e+00 0.00000  4.797658  0.019393551 0.031076189        0.9242092
#16     5        1 0.000000e+00 0.00000  5.193476  0.015911331 0.025496286        0.9497055
#17    12        1 1.790286e-04 0.51753  5.877753  0.011301021 0.018108734        0.9678142
#18    15        1 5.556320e-05 0.51670  6.704368  0.007475157 0.011978178        0.9797924
#19    16        1 2.808353e-04 0.51776  6.939440  0.006646227 0.010649901        0.9904423
#20     7        1 0.000000e+00 0.00000  7.155848  0.005964619 0.009557695        1.0000000


#So if we assume a constant baseline rate, we get strong evidence of a strong social transmission effect.
#If we assume a weibull baseline function, we get at least reasonable evidence of social transmission in all models
#since in all but the weakest model the 95% CI does not include 0, but in some models the effect is plausibly quite weak (low
#propST).
#If we assume a gamma baseline function, many models have 95% CI including 0, or very low plausible effects for social transmission.

#In cases like this, where the conclusions are strongly dependent on the assumptions about baseline rate,
#it shows that the analysis is dominated by the time course of events rather than the pattern of spread through the network.
#Unless there is a good a priori reason to believe that the baseline rate is constant, we would recommend switching to OADA in
#such cases, which ignores the time course of events and is only sensitive to the pattern of spread through the network.


#############################################################################
# TUTORIAL 7.4
# MULTI-MODEL INFERENCE IN AN OADA WITH MULTIPLE NETWORKS
# 1 diffusion
# 2 static networks
# 2 time constant ILVs
#############################################################################


#Read in the 2 social networks and order of acquisition vector as shown in Tutorial 3
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
#If all of our networks are static (do not change over time) this is a three-dimensional array of size
#no. individuals x no.individuals x no.networks
socNets<-array(NA,dim=c(30,30,2))
#Then slot the networks into the array:
socNets[,,1]<-socNet1
socNets[,,2]<-socNet2

#Then we go on to create our nbdaData object as before. Let us assume we are interested in the unconstrained model:
nbdaData_multiNet<-nbdaData(label="ExampleDiffusion2",assMatrix=socNets,orderAcq = oa1,asoc_ilv = asoc,int_ilv = asoc )
#Then fit the model:
model1_multiNet<-oadaFit(nbdaData_multiNet)
#And get the output:
data.frame(Variable=model1_multiNet@varNames,MLE=model1_multiNet@outputPar,SE=model1_multiNet@se)

#Let us assume we have 4 competing hypotheses about social transmission.
#1. Transmission through network 1 only (s2=0)
#2. Transmission through network 2 only (s1=0)
#3. Transmission through network 1 and network 2 at equal rates per unit connection (s1=s2)
#4. Transmission through network 1 and network 2 at different rate (no constraint)

#Each of these can be represented by different constraints between the s parameters, which we refer to as a network combination or netcombo:
#1. 1 0
#2. 0 1
#3. 1 1
#4. 1 2

#We can fit models with each netcombo and compare the fit. However, we are unsure which ILVs to include, so we want to use multi-model inference
#We need to set up a constraintsVectMatrix that considers each netcombo but also every combination of ILVs

constraintsVectMatrix<-rbind(
  #netcombo 1 0
  c(1,0,0,0,0,0),
  c(1,0,0,0,0,2),
  c(1,0,0,0,2,0),
  c(1,0,0,0,2,3),
  c(1,0,0,2,0,0),
  c(1,0,0,2,0,3),
  c(1,0,0,2,3,0),
  c(1,0,0,2,3,4),
  c(1,0,2,0,0,0),
  c(1,0,2,0,0,3),
  c(1,0,2,0,3,0),
  c(1,0,2,0,3,4),
  c(1,0,2,3,0,0),
  c(1,0,2,3,0,4),
  c(1,0,2,3,4,0),
  c(1,0,2,3,4,5),
  
  #netcombo 0 1
  c(0,1,0,0,0,0),
  c(0,1,0,0,0,2),
  c(0,1,0,0,2,0),
  c(0,1,0,0,2,3),
  c(0,1,0,2,0,0),
  c(0,1,0,2,0,3),
  c(0,1,0,2,3,0),
  c(0,1,0,2,3,4),
  c(0,1,2,0,0,0),
  c(0,1,2,0,0,3),
  c(0,1,2,0,3,0),
  c(0,1,2,0,3,4),
  c(0,1,2,3,0,0),
  c(0,1,2,3,0,4),
  c(0,1,2,3,4,0),
  c(0,1,2,3,4,5),
  
  #netcombo 1 1
  c(1,1,0,0,0,0),
  c(1,1,0,0,0,2),
  c(1,1,0,0,2,0),
  c(1,1,0,0,2,3),
  c(1,1,0,2,0,0),
  c(1,1,0,2,0,3),
  c(1,1,0,2,3,0),
  c(1,1,0,2,3,4),
  c(1,1,2,0,0,0),
  c(1,1,2,0,0,3),
  c(1,1,2,0,3,0),
  c(1,1,2,0,3,4),
  c(1,1,2,3,0,0),
  c(1,1,2,3,0,4),
  c(1,1,2,3,4,0),
  c(1,1,2,3,4,5),
  
  #netcombo 1 2
  c(1,2,0,0,0,0),
  c(1,2,0,0,0,3),
  c(1,2,0,0,3,0),
  c(1,2,0,0,3,4),
  c(1,2,0,3,0,0),
  c(1,2,0,3,0,4),
  c(1,2,0,3,4,0),
  c(1,2,0,3,4,5),
  c(1,2,3,0,0,0),
  c(1,2,3,0,0,4),
  c(1,2,3,0,4,0),
  c(1,2,3,0,4,5),
  c(1,2,3,4,0,0),
  c(1,2,3,4,0,5),
  c(1,2,3,4,5,0),
  c(1,2,3,4,5,6),
  
  #netcombo 0 0 
  c(0,0,0,0,0,0),
  c(0,0,1,0,0,0),
  c(0,0,0,1,0,0),
  c(0,0,1,2,0,0)
  )

modelSet_multiNet<-oadaAICtable(nbdaData_multiNet,constraintsVectMatrix = constraintsVectMatrix)

networksSupport(modelSet_multiNet)

#support numberOfModels
#0:0 0.01065492              4
#0:1 0.00818532             16
#1:0 0.67992442             16
#1:1 0.11457221             16
#1:2 0.18666313             16

#Note that we have equal numbers of models for each network combo save for asocial learning (0 0), so we can compare the support for different models
#of social transmission. However, we recommend examining the CIs for s parameters to judge the evidence for social transmission versus asocial learning (see main text).
#We can see strongest support (68%) for social transmission only through network 1, and second most (18.7%) for transmission of different strengths through each network

#We can also get MAEs etc as before:

rbind(
  support=variableSupport(modelSet_multiNet),
  MAE=modelAverageEstimates(modelSet_multiNet),
  USE=unconditionalStdErr(modelSet_multiNet))

#              s1        s2   ASOC:male ASOC:stAge SOCIAL:male SOCIAL:stAge
#support 0.9811598 0.3094207   0.5935977  0.3502696  0.22537376    0.3445204
#MAE     3.6470815 0.1433855 -11.2021988 -0.3394773 -0.02000315    0.1031443
#USE           NaN       NaN         NaN  0.5557081         NaN          NaN

#The missing USEs are due to some models in which the SEs could not be derived. If we look at
print(modelSet_multiNet)
#We can see that the SEs are present in all the models with high weight, so it seems reasonable to get an approximate USE by replacing NAs and NaNs with the weighted mean across
#all other models:

rbind(
  support=variableSupport(modelSet_multiNet),
  MAE=modelAverageEstimates(modelSet_multiNet),
  USE=unconditionalStdErr(modelSet_multiNet,nanReplace = T))

#              s1        s2     ASOC:male    ASOC:stAge SOCIAL:male SOCIAL:stAge
#support  0.9811598 0.3094207  5.935977e-01  0.3502696  0.22537376    0.3445204
#MAE      3.6470815 0.1433855 -1.120220e+01 -0.3394773 -0.02000315    0.1031443
#USE     36.0285332 1.0059692  7.310227e+07  0.5557081  0.14394848    0.1284004

#Averaged across models we can see that s1 is estimated to be considerably greater than s2.
#We would recommend presenting the 95% CI for s1 and s2 in the best model in which they are present
#and the 95% CI for s1-s2 in the best model in which they are both present and unconstrained.

#We can check the sensitivity of the 95% CIs for s1 to model selection uncertainty as before:

lowerLimitsByModel<-multiModelLowerLimits(which=1,aicTable = modelSet_multiNet,conf=0.95)
lowerLimitsByModel

#   model netCombo      lowerCI  propST deltaAICc akaikeWeight  adjAkWeight cumulAdjAkWeight
#1      9      1:0 3.899414e-01 0.53571 0.0000000 0.1315365430 0.1340623086        0.1340623
#2      1      1:0 3.007326e-01 0.47783 0.8673898 0.0852501433 0.0868871171        0.2209494
#3     13      1:0 4.020928e-01 0.54412 1.1352833 0.0745628988 0.0759946561        0.2969441
#4     10      1:0 4.086640e-01 0.54993 1.2875510 0.0690968545 0.0704236527        0.3673677
#5      2      1:0 2.690855e-01 0.45299 1.9136852 0.0505236774 0.0514938334        0.4188616
#6     14      1:0 1.452914e-01 0.38811 2.1230682 0.0455017414 0.0463754662        0.4652370
#7     11      1:0 9.396559e-02 0.30242 2.4186965 0.0392494122 0.0400030796        0.5052401
#8     57      1:2 3.743582e-01 0.73650 2.4786325 0.0380906356 0.0388220522        0.5440622
#9      5      1:0 4.037247e-01 0.54338 2.5009361 0.0376682150 0.0383915203        0.5824537
#10     3      1:0 3.522655e-01 0.51316 2.8398628 0.0317964081 0.0324069629        0.6148606
#11     6      1:0 5.757081e-02 0.23110 3.0352408 0.0288371449 0.0293908759        0.6442515
#12    49      1:2 3.007326e-01 0.47783 3.1689771 0.0269719152 0.0274898301        0.6717414
#13    41      1:1 1.573432e-01 0.41185 3.5213567 0.0226148409 0.0230490912        0.6947904
#14    15      1:0 1.242905e-01 0.34397 3.7513244 0.0201584287 0.0205455110        0.7153360
#15    61      1:2 3.639542e-01 0.74839 3.8122064 0.0195540321 0.0199295088        0.7352655
#16    12      1:0 4.061835e-02 0.26161 3.9587446 0.0181725536 0.0185215031        0.7537870
#17    58      1:2 4.066797e-01 0.74767 3.9644741 0.0181205685 0.0184685197        0.7722555
#18    45      1:1 1.856705e-01 0.50556 4.0029245 0.0177755243 0.0181168501        0.7903723
#19    33      1:1 1.274416e-02 0.04379 4.1942219 0.0161540980 0.0164642891        0.8068366
#20     4      1:0 2.374518e-01 0.43650 4.2344673 0.0158322823 0.0161362939        0.8229729
#21    50      1:2 2.690855e-01 0.46200 4.3923177 0.0146307554 0.0149116953        0.8378846
#22     7      1:0 4.763377e-01 0.58550 4.6267354 0.0130125853 0.0132624530        0.8511471
#23    53      1:2 4.031556e-01 0.54471 4.9795686 0.0109080429 0.0111174992        0.8622646
#24    16      1:0 1.909823e-02 0.24849 5.0137907 0.0107229829 0.0109288857        0.8731935
#25    62      1:2 1.452915e-01 0.45126 5.0230682 0.0106733566 0.0108783065        0.8840718
#26    59      1:2 9.396390e-02 0.61701 5.0956195 0.0102931120 0.0104907604        0.8945625
#27    51      1:2 3.522655e-01 0.52724 5.3184952 0.0092076724 0.0093844782        0.9039470
#28    42      1:1 1.565209e-01 0.41558 5.5865516 0.0080527119 0.0082073401        0.9121543
#29     8      1:0 1.176647e-02 0.19174 5.5990505 0.0080025438 0.0081562087        0.9203105
#30    54      1:2 5.757082e-02 0.29719 5.7121639 0.0075625072 0.0077077225        0.9280183
#31    37      1:1 5.059767e-02 0.14921 5.7482011 0.0074274614 0.0075700836        0.9355884
#32    46      1:1 7.603493e-02 0.37244 5.9165026 0.0068280110 0.0069591225        0.9425475
#33    43      1:1 2.559226e-03 0.52095 5.9343480 0.0067673575 0.0068973043        0.9494448
#34    34      1:1 1.081990e-02 0.15035 6.1186879 0.0061714929 0.0062899979        0.9557348
#36    35      1:1 9.875455e-03 0.15061 6.4670343 0.0051849882 0.0052845504        0.9610193
#37    47      1:1 2.528128e-02 0.55149 6.5904721 0.0048746518 0.0049682549        0.9659876
#38    63      1:2 1.242834e-01 0.64277 6.6513244 0.0047285684 0.0048193664        0.9708070
#39    60      1:2 4.062878e-02 0.57001 6.8586465 0.0042629502 0.0043448074        0.9751518
#40    52      1:2 2.374518e-01 0.48312 6.9113904 0.0041519973 0.0042317240        0.9793835
#41    38      1:1 4.534872e-05 0.11161 7.2380830 0.0035262790 0.0035939906        0.9829775
#42    55      1:2 4.753915e-01 0.58807 7.3036585 0.0034125351 0.0034780626        0.9864555
#43    44      1:1 6.649519e-05 0.51734 8.0635604 0.0023338155 0.0023786295        0.9888342
#44    64      1:2 1.909815e-02 0.33128 8.1659646 0.0022173270 0.0022599041        0.9910941
#45    39      1:1 4.755244e-02 0.22975 8.1991606 0.0021808276 0.0022227038        0.9933168
#49    48      1:1 6.307494e-05 0.28074 8.4117525 0.0019609096 0.0019985630        0.9953153
#50    56      1:2 1.176656e-02 0.22252 8.4990505 0.0018771590 0.0019132042        0.9972285
#51    36      1:1 3.834547e-02 0.25710 8.5969567 0.0017874792 0.0018218024        0.9990503
#53    40      1:1 7.916526e-05 0.16341 9.8999189 0.0009317649 0.0009496567        1.0000000

#All 95% CI are above zero, adding to our confidence that we have evidence for social transmission through network 1


#We can obtain a model averaged lower-limit for propST as follows:

sum(lowerLimitsByModel$propST*lowerLimitsByModel$adjAkWeight)
#[1] 0.5473627

#However, note that this includes models with netCombo= 1 1- i.e. one in which the value of s1 is constrained to be = s2. 
#Our preferred approach is to refit the model set including only those models in which the value of s1 is unconstrained, and re-run the exercise:


constraintsVectMatrixs1<-rbind(
  #netcombo 1 0
  c(1,0,0,0,0,0),
  c(1,0,0,0,0,2),
  c(1,0,0,0,2,0),
  c(1,0,0,0,2,3),
  c(1,0,0,2,0,0),
  c(1,0,0,2,0,3),
  c(1,0,0,2,3,0),
  c(1,0,0,2,3,4),
  c(1,0,2,0,0,0),
  c(1,0,2,0,0,3),
  c(1,0,2,0,3,0),
  c(1,0,2,0,3,4),
  c(1,0,2,3,0,0),
  c(1,0,2,3,0,4),
  c(1,0,2,3,4,0),
  c(1,0,2,3,4,5),

  #netcombo 1 2
  c(1,2,0,0,0,0),
  c(1,2,0,0,0,3),
  c(1,2,0,0,3,0),
  c(1,2,0,0,3,4),
  c(1,2,0,3,0,0),
  c(1,2,0,3,0,4),
  c(1,2,0,3,4,0),
  c(1,2,0,3,4,5),
  c(1,2,3,0,0,0),
  c(1,2,3,0,0,4),
  c(1,2,3,0,4,0),
  c(1,2,3,0,4,5),
  c(1,2,3,4,0,0),
  c(1,2,3,4,0,5),
  c(1,2,3,4,5,0),
  c(1,2,3,4,5,6))

modelSet_s1<-oadaAICtable(nbdadata = nbdaData_multiNet,constraintsVectMatrix =constraintsVectMatrixs1)
lowerLimitsByModel_s1<-multiModelLowerLimits(which=1,aicTable = modelSet_s1,conf=0.95)
lowerLimitsByModel_s1

sum(lowerLimitsByModel_s1$propST*lowerLimitsByModel_s1$adjAkWeight)
#[1] 0.5857381

#Giving a model-averaged lower limit of 58.6% of events occurring as a result of social transmission via network 1


#We can now do the same for s2:

constraintsVectMatrixs2<-rbind(
  #netcombo 0 1
  c(0,1,0,0,0,0),
  c(0,1,0,0,0,2),
  c(0,1,0,0,2,0),
  c(0,1,0,0,2,3),
  c(0,1,0,2,0,0),
  c(0,1,0,2,0,3),
  c(0,1,0,2,3,0),
  c(0,1,0,2,3,4),
  c(0,1,2,0,0,0),
  c(0,1,2,0,0,3),
  c(0,1,2,0,3,0),
  c(0,1,2,0,3,4),
  c(0,1,2,3,0,0),
  c(0,1,2,3,0,4),
  c(0,1,2,3,4,0),
  c(0,1,2,3,4,5),

  #netcombo 1 2
  c(1,2,0,0,0,0),
  c(1,2,0,0,0,3),
  c(1,2,0,0,3,0),
  c(1,2,0,0,3,4),
  c(1,2,0,3,0,0),
  c(1,2,0,3,0,4),
  c(1,2,0,3,4,0),
  c(1,2,0,3,4,5),
  c(1,2,3,0,0,0),
  c(1,2,3,0,0,4),
  c(1,2,3,0,4,0),
  c(1,2,3,0,4,5),
  c(1,2,3,4,0,0),
  c(1,2,3,4,0,5),
  c(1,2,3,4,5,0),
  c(1,2,3,4,5,6))

modelSet_s2<-oadaAICtable(nbdadata = nbdaData_multiNet,constraintsVectMatrix =constraintsVectMatrixs2)

lowerLimitsByModel_s2<-multiModelLowerLimits(which=2,aicTable = modelSet_s2,conf=0.95)
lowerLimitsByModel_s2

#model netCombo    lowerCI  propST  deltaAICc akaikeWeight  adjAkWeight cumulAdjAkWeight
#1     25      1:2 0.00000000 0.00000  0.0000000 0.1954885179 0.1954885179        0.1954885
#2     17      1:2 0.00000000 0.00000  0.6903446 0.1384250918 0.1384250918        0.3339136
#3     29      1:2 0.00000000 0.00000  1.3335739 0.1003550795 0.1003550795        0.4342687
#4     26      1:2 0.00000000 0.00000  1.4858416 0.0929982664 0.0929982664        0.5272670
#5     18      1:2 0.00000000 0.00000  1.9136852 0.0750878699 0.0750878699        0.6023548
#6     21      1:2 0.00000000 0.00000  2.5009361 0.0559821884 0.0559821884        0.6583370
#7     30      1:2 0.00000000 0.00000  2.5444357 0.0547777329 0.0547777329        0.7131147
#8     27      1:2 0.00000000 0.00000  2.6169871 0.0528262440 0.0528262440        0.7659410
#9     19      1:2 0.00000000 0.00000  2.8398627 0.0472555577 0.0472555577        0.8131965
#10    22      1:2 0.00000000 0.00000  3.2335314 0.0388122512 0.0388122512        0.8520088
#11    31      1:2 0.00000000 0.00000  4.1726919 0.0242679289 0.0242679289        0.8762767
#12    28      1:2 0.00000000 0.00000  4.3800140 0.0218782859 0.0218782859        0.8981550
#13    20      1:2 0.00000000 0.00000  4.4327579 0.0213088542 0.0213088542        0.9194639
#14    23      1:2 0.00000000 0.00000  4.8250260 0.0175137910 0.0175137910        0.9369777
#15    32      1:2 0.00000000 0.00000  5.6873321 0.0113797515 0.0113797515        0.9483574
#16     1      0:1 0.00000000 0.00000  5.8518172 0.0104813034 0.0104813034        0.9588387
#17    24      1:2 0.00000000 0.00000  6.0204180 0.0096339435 0.0096339435        0.9684727
#18     2      0:1 0.00000000 0.00000  7.3398351 0.0049807681 0.0049807681        0.9734534
#19    13      0:1 0.00000000 0.00000  7.9931050 0.0035928621 0.0035928621        0.9770463
#20     9      0:1 0.00000000 0.00000  8.0761826 0.0034466762 0.0034466762        0.9804930
#21     5      0:1 0.00000000 0.00000  8.1326620 0.0033507045 0.0033507045        0.9838437
#22     3      0:1 0.00000000 0.00000  8.1453457 0.0033295221 0.0033295221        0.9871732
#23     4      0:1 0.00000000 0.00000  8.9206772 0.0022595440 0.0022595440        0.9894327
#24     6      0:1 0.00000000 0.00000  9.3040168 0.0018654331 0.0018654331        0.9912982
#25    10      0:1 0.00000000 0.00000  9.6072278 0.0016030174 0.0016030174        0.9929012
#26    11      0:1 0.00000000 0.00000  9.8584312 0.0014138070 0.0014138070        0.9943150
#27    15      0:1 0.00000000 0.00000  9.9529192 0.0013485663 0.0013485663        0.9956636
#28    14      0:1 0.00000000 0.00000 10.3194716 0.0011227329 0.0011227329        0.9967863
#29     7      0:1 0.00000000 0.00000 10.6018689 0.0009748879 0.0009748879        0.9977612
#30     8      0:1 0.00000000 0.00000 10.7016820 0.0009274287 0.0009274287        0.9986886
#31    12      0:1 0.00000000 0.00000 11.0238742 0.0007894368 0.0007894368        0.9994780
#32    16      0:1 0.07159997 0.44518 11.8513525 0.0005219546 0.0005219546        1.0000000

#In this case, for all models except the worst one the 95% CI for s2 contain 0, so we can be sure we do not have good evidence for social
#transmission via network 2.

sum(lowerLimitsByModel_s2$propST*lowerLimitsByModel_s2$adjAkWeight)
#[1] 0.0002323638

#############################################################################
# TUTORIAL 7.5
# MULTI-MODEL INFERENCE IN A cTADA: INCLUDING A HOMOGENEOUS NETWORK
# 1 diffusion
# 1 static network
# 2 time constant ILVs
#############################################################################

# We can include multiple networks in a TADA in the same way as described for an OADA in tutorial 7.4.
# Furthermore, in a TADA, even when we only have one social network of interest, we may wish to include
# a homogeneous network in our multi-model procedure (see tutorial 4.3).
# If we have multiple diffusions and use a cTADA we may also want to include a group network (see tutorial
# 6.5)
# In this tutorial we show how to include a homogeneous network in a multi-model inference procedure for
# a cTADA of a single diffusion. A group network could be included in a cTADA/ stratified OADA of multiple
# diffusions in the same way.


#We use the same data as for tutoral 7.3:
#########################################
#Read in the social network and order of acquisition vector as shown in Tutorials 1 and 2
socNet1<-as.matrix(read.csv(file="exampleStaticSocNet.csv"))
socNet1<-array(socNet1,dim=c(30,30,1))
oa1<-c(26,29,30,8,19,21,22,3,14,12,11,1,17,28,5,9,15,7,6,25,4,13,27,18,20,24,23,16,2,10)

ILVdata<-read.csv(file="exampleTimeConstantILVs.csv")
#Then extract each ILV
female<-cbind(ILVdata$female)
male<-1-female
#females=0 males=1
age<-cbind(ILVdata$age)
#age in years
#I am going to standardize age
stAge<-(age-mean(age))/sd(age)

asoc<-c("male","stAge")

#Now in addition, we need a vector giving the times at which each individual acquired the behaviour:
ta1<-c(234,252,262,266,273,296,298,310,313,326,332,334,334,337,338,340,343,367,374,376,377,393,402,405,407,407,435,472,499,567)
#######################################

#Now we add in the homogeneous network as a second network:
dim(socNet1)
socNets<-array(NA,dim=c(30,30,2))
socNets[,,1]<-socNet1
socNets[,,2]<-1


#Create an object for the "unconstrained" model
nbdaData_tada<-nbdaData(label="ExampleDiffusion2",assMatrix=socNets,orderAcq = oa1,asoc_ilv = asoc,int_ilv=asoc,timeAcq = ta1,endTime = 568)

#Now we have an object with two networks- the first is the social network, the second the homogeneous network

#Now for the model set, we consider all possible sets of ILVs in an unconstrained model for both social and
#homogeneous networks and an asocial model too:

constraintsVectMatrix<-rbind(
  #social network models (1 in first slot, 0 in second)
  c(1,0,0,0,0,0),
  c(1,0,0,0,0,2),
  c(1,0,0,0,2,0),
  c(1,0,0,0,2,3),
  c(1,0,0,2,0,0),
  c(1,0,0,2,0,3),
  c(1,0,0,2,3,0),
  c(1,0,0,2,3,4),
  c(1,0,2,0,0,0),
  c(1,0,2,0,0,3),
  c(1,0,2,0,3,0),
  c(1,0,2,0,3,4),
  c(1,0,2,3,0,0),
  c(1,0,2,3,0,4),
  c(1,0,2,3,4,0),
  c(1,0,2,3,4,5),
  #homogeneous network models (0 in first slot, 1 in second)
  c(0,1,0,0,0,0),
  c(0,1,0,0,0,2),
  c(0,1,0,0,2,0),
  c(0,1,0,0,2,3),
  c(0,1,0,2,0,0),
  c(0,1,0,2,0,3),
  c(0,1,0,2,3,0),
  c(0,1,0,2,3,4),
  c(0,1,2,0,0,0),
  c(0,1,2,0,0,3),
  c(0,1,2,0,3,0),
  c(0,1,2,0,3,4),
  c(0,1,2,3,0,0),
  c(0,1,2,3,0,4),
  c(0,1,2,3,4,0),
  c(0,1,2,3,4,5),
  #asocial models
  c(0,0,0,0,0,0),
  c(0,0,1,0,0,0),
  c(0,0,1,2,0,0),
  c(0,0,0,1,0,0))

#Again, in a TADA we likely want to consider different baseline functions too. 
#To speed things up we will just consider the constant and weibull baseline functions here,
#so we want to replicate the constraintsVectMartix above 2 times:
constraintsVectMatrixAll<-rbind(constraintsVectMatrix,constraintsVectMatrix)

#And we also need to provide a vector saying what baseline function each model will have.


baselineVect<-rep(c("constant","weibull"),each=dim(constraintsVectMatrix)[1])

#Check:
cbind(baselineVect,constraintsVectMatrixAll)

#To fit the set of models, we use
modelSet_tada<-tadaAICtable(nbdadata = nbdaData_tada,constraintsVectMatrix =constraintsVectMatrixAll,baselineVect = baselineVect)


#We can then get support for each combination of networks across all models
networksSupport(modelSet_tada)
#       support numberOfModels
#0:0 0.02347081              8
#0:1 0.03745207             32
#1:0 0.93907712             32

#Here we can see that the social network (1:0= 93.9%) has far more support than the homogeneous network (0:1= 3.7%) and the number of
#models is the same making this a fair comparison. The ratio of support is
0.93907712/0.03745207 
#One could report this result as something like:
# "25.1x more support for the social network than the homogeneous network, providing evidence that the diffusion follows the
# social network."

#One could then proceed to make inferences as in tutorial 7.3 above. Given the different scales of the networks, it would be reasonable
#to refit the model set to just include the social network before calculating model-averaged estimates. However, if the support for
#the homogeneous network is tiny (e.g. <0.1%) then it is unlikely to affect the model-averaged estimates at all.
