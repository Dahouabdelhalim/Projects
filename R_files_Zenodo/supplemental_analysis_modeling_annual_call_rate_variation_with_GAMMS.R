##  Load data and packages  ###############
library(lme4)
library(lmerTest)
library(effects)
library(car)

load("multiplicata calls with temp corrections")

#fix population issue: BIP-1 and BIP-2 are ~.1m apart and should be considered same population
#because they are already factors and I don't want to fuck around with the factor levels, 
#naming them both "BIP-1" for purposes of this analysis.
data.m[data.m$pop=="BIP-2",1]<-"BIP-1"
data.m<-droplevels(data.m)
levels(data.m$pop) #worked--21 levels (populations) with all BIP collapsed into BIP-1. 
######################################

## A.  Confirm effect of year-of-recording 

##############################
#We know from ongoing work that call rates have been changing over years, want to control for 
#this variation.  First confirm that call rates change over time:

year.mod<-lmer("call.rate.tempcor ~ year + (1|pop)",
               data = data.m)
summary(year.mod) 
anova(year.mod, ddf = "Kenward-Roger")#sig increase with year on temp corrected calls

#Thus, we should include year as a fixed effect in subsequent models.  

#Modeling with a GAM:
library(mgcv)
year.mod.gam<-gamm(call.rate.tempcor ~ s(year) + s(pop, bs = "re"), data = data.m)
summary(year.mod.gam$lme) #get parameter estimates
summary(year.mod.gam$gam) #get sig. smooth terms. significant effect of year as expected

#compare models:
AIC(year.mod) #2864.174
AIC(year.mod.gam$lme) #2861.653 

#gam has smaller AIC but delta AIC is 2.5 so it doesn't matter much...
summary(year.mod.gam$gam)$r.sq #0.2921821
library(MuMIn)
r.squaredGLMM(year.mod) #0.288393 #basically the same R-squared.
 
#GAM doesn't offer a much better fit than linear but let's check to make sure the results are the 
#same either way
######################################


## B.  Do call rates differ in sympatry and allopatry? 

##############################
#Assign pops to sympatry or allopatry.
#unique(data.m$pop)
sym.pops<-c("Guy Miller", "410", "Shrimp", "Sulphur Draw", "Javelina", "Stateline Rd", "Mesquite", "Zent", "Sky Ranch")
allo.pops<-c("Acacia", "Dead Cow", "Peach Orchard", "Horseshoe", "BIP-1", "Crissal", "Observatory", "Crater", "Oneal Ditch", 
"Pampas", "Silver Creek", "SkeletonCanyon")


topy<-rep(NA, times = nrow(data.m))
for (i in 1:length(topy)) {
  pop.var<-data.m[i,1]
  if(pop.var %in% sym.pops) {topy[i]<-"sym"}
  else if (pop.var %in% allo.pops ) {topy[i]<- "allo"}
  else {topy[i]<-NA}
}
data.m<-cbind.data.frame(data.m, topy)

topy.mod<-lmer(call.rate.tempcor ~ year + topy + (1|pop), data = data.m)
summary(topy.mod)
anova(topy.mod, ddf = "Kenward-Roger")
#sym and allo marginally diff in call rate (p = 0.09)

#checking with gam:
topy.mod.gam<-gamm(call.rate.tempcor ~ s(year) + topy + s(pop, bs = "re"), data = data.m)
summary(topy.mod.gam$lme) #get estimates and SE
summary(topy.mod.gam$gam) #get sig. of smooth terms
anova(topy.mod.gam$gam) #get p and F of parametric terms

#same results qualitatively: sig effect of year, marginal effect of population type.


######################################

## C.  Effect of Elevation on call rates 

######################################
year.elev.mod<-lmer(call.rate.tempcor ~ elevation + year + (1|pop), data = data.m)
summary(year.elev.mod)
anova(year.elev.mod, ddf = "Kenward-Roger")
#year, elev both sig

#checking with gamm:
year.elev.gam<-gamm(call.rate.tempcor ~ elevation + s(year) + s(pop, bs = "re"), data = data.m)
summary(year.elev.gam$lme)$tTable #get parameter estimates

summary(year.elev.gam$gam) #get p and F of smooth terms: sig. effect of year 
anova(year.elev.gam$gam) #get p and F of parametric terms
######################################

## E.  Does Female preference predict male call rate across populations?

######################################
df.choice<-read.csv(file = "multiplicata_fvs_choice_updated2020.csv", header = T)

#Clean and organize choice data, and combine with calls df: #########
#Only interested in high-water testing condition for this analysis, removing low-water testing columns:
df.choice<-df.choice[,c(1:17)]

#
#Remove hybrids/likely hybrids/any frog with a boss or slight boss (S. bombifrons characters) 
#as per Karin's visual inspection:

#1746_LF1 = 1721_lf1 (duplicated test. Same female,, she changed cages and therefore cage number.)  
#Remove 1746_LF1 test, she chose "S" in both tests so doesn't matter which one is dropped, just going w/ most recent.
df.choice<-df.choice[df.choice$ID!="1746_lf1",]
#Likely hybrids or hybrids:
df.choice<-df.choice[df.choice$ID!="7803_lf3_lf4",]
df.choice<-df.choice[df.choice$ID!="7812_lf2",]
df.choice<-df.choice[df.choice$ID!="7813_lf1",]
#Frogs Karin classified as SM with boss or slight boss:
df.choice<-df.choice[df.choice$ID!="1721_LF1",]
df.choice<-df.choice[df.choice$ID!="1757_lf1",]
df.choice<-df.choice[df.choice$ID!="1759_lf3",]
df.choice<-df.choice[df.choice$ID!="12102",]
df.choice<-df.choice[df.choice$ID!="10279",]
df.choice<-df.choice[df.choice$ID!="11416",]

#Excluding 'no choice' females from analysis:
  df.choice<-df.choice[(df.choice$mate_choice_hi!="NC")&is.na(df.choice$mate_choice_hi)==F,]
#re-factor choices so slow = 0, fast = 1:
df.choice$mate_choice_hi<-factor(df.choice$mate_choice_hi, levels = c("S","F"))
#ensure populations are treated as factors
df.choice$pop<-factor(df.choice$pop)
df.choice<-droplevels(df.choice)

#Fixing population nameing issues:
#PeachOrchard and PeachOrchardRd2 are the same population:
which(df.choice$pop == "PeachOrchardRd2") #row 84
df.choice[84,4]<-"PeachOrchard"
#Javelina and Lazy River are the same population:
which(df.choice$pop == "Lazy River") 
df.choice[c(129:133),4]<-"Javelina"

#Get number choosing fast and slow for each population and add to df.choice:
#add prop.choose.f column
pop<-unique(df.choice$pop)
empty<-rep(NA, times = length(unique(df.choice$pop)))
N.chose.f<-empty
N.chose.s<-empty

for (i in 1:length(pop)) {
  pop.var<-pop[i]
  subset<-df.choice[df.choice$pop==pop.var,13]
  N.chose.f[i]<-length(subset[subset=="F"])
  N.chose.s[i]<-length(subset[subset=="S"])
}
bypopdf<-cbind.data.frame(pop,N.chose.f, N.chose.s)

#make sure population names are transferrable between the data frames
cbind.data.frame(unique(data.m$pop),unique(data.m$pop) %in% df.choice$pop)

#pops in calls data that are not in choice data: Skeleton Canyon, Shrimp, Stateline Rd, 
#Crissal, Silver Creek.  These will be left as NA values for the append matrix.

#fixing pops with different naming conventions between the datasets:
bypopdf$pop<-as.character(bypopdf$pop)
bypopdf$pop[25]<-"Oneal Ditch"  #Oneal Ditch = Oneal Rd
bypopdf$pop[15]<-"Peach Orchard"  #PeachOrchard = Peach Orchard

#combine choice and calls df:
append.df<-data.frame(matrix(NA, nrow = nrow(data.m), ncol = 2))
colnames(append.df)<-names(bypopdf)[c(2:3)]
for (i in 1:nrow(data.m)){
  pop.var <-data.m[i,1]
  if(pop.var=="BIP-1") {   #BIP1 and BIP2 are same population
    pop.var<-"BIP1"
  }
  if(pop.var=="BIP-2") {   #BIP1 and BIP2 are same population
    pop.var<-"BIP1"
  }
  if(pop.var %in% bypopdf$pop == F) {
    next
  }
  subset<-bypopdf[bypopdf$pop==pop.var,c(2,3)]
  append.df[i,]<-subset
}
data.m.choice<-cbind.data.frame(data.m,append.df)
data.m.choice$prop.chose.f<-(data.m.choice$N.chose.f/(data.m.choice$N.chose.f+data.m.choice$N.chose.s))
data.m.choice$N.females<-(data.m.choice$N.chose.f+data.m.choice$N.chose.s) 
############


#Test whether female preferences correlate with male call rates: 

#remove pops where we have no female choice data:
data.m.choice2<-droplevels(data.m.choice[is.na(data.m.choice$N.females)==F,])

#restrict analysis to populations with at least 4 females sampled
data.m.choice.sub<-data.m.choice2[data.m.choice2$N.females>3,]


#Now, see if preferences predict call rates.  
choice.mod1.gam<-gamm(call.rate.tempcor ~ car::logit(prop.chose.f) + s(year) + elevation + s(pop, bs = "re"), data = data.m.choice.sub)
choice.mod2.gam<-gamm(call.rate.tempcor ~ car::logit(prop.chose.f) + s(year) + s(pop, bs = "re"), data = data.m.choice.sub)

summary(choice.mod1.gam$lme) #get parameter estimates
summary(choice.mod2.gam$lme)

summary(choice.mod1.gam$gam)# get smoothed term p and F
summary(choice.mod2.gam$gam)# year sig in both models

anova(choice.mod1.gam$gam) #get parametric terms p and F
anova(choice.mod2.gam$gam) # female preference NS in both models.

