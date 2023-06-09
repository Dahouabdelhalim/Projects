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

#reviewer suggested analysis with only allopatric populations
#making a datafframe of only allpatric populations, then we will run all 
#subsequent analyses with all pops and allotopic only, and compare results.  
data.m.allo<-data.m[data.m$topy == "allo",] 

######################################

## A.  Confirm effect of year-of-recording 

##############################
#We know from ongoing work that call rates have been changing over years, want to control for 
#this variation.  First confirm that call rates change over time:
year.mod<-lmer("call.rate.tempcor ~ year + (1|pop)",
                 data = data.m.allo)
summary(year.mod)
anova(year.mod, ddf = "Kenward-Roger")

#Thus, we should include year as a fixed effect in subsequent models.  

#skipping B., comparison of call rates in sympatry and allopatry because, for an allopatry-only
#dataset, this is not possible

##############################

## C.  Effect of male genotype on call rates 

##############################
#could sympatric calls be faster than allopatric due to introgression from fast-calling 
#S. bombifrons?

#Test whether pure-M (males scoring 5) and introgressed (males scoring from 4 to 5) have 
#different call rates:
data.pure.m.allo<-data.m.allo[(is.na(data.m.allo$male.genotype)==F & data.m.allo$male.genotype>4.99),]
data.introg.m.allo<-data.m.allo[(is.na(data.m.allo$male.genotype)==F & data.m.allo$male.genotype<4.99),]

t.test(x=data.pure.m.allo$call.rate.tempcor, y=data.introg.m.allo$call.rate.tempcor) #not sig diff 
#no difference between pure m and introgressed males in call rate for allopatric pops. 


######################################

## D.  Effect of Elevation on call rates 

######################################
year.elev.mod<-lmer(call.rate.tempcor ~ elevation + year + (1|pop), data = data.m.allo)
summary(year.elev.mod)
anova(year.elev.mod, ddf = "Kenward-Roger")
#year sig in both, elevation marginal for allopatry-only analysis. In addition to cutting the sample size
#by looking at allopatry only we are also sampling across a narrower range of elevations
#since sympatry is low-elevation and allopatry is high-elevation.

######################################

## E.  Does Female preference predict male call rate across populations?

######################################
######################################
df.choice<-read.csv(file = "multiplicata_fvs_choice_updated2020.csv", header = T)

#Clean and organize choice data, and combine with calls df: #########
#Only interested in high-water testing condition for this analysis, removing low-water testing columns:
df.choice<-df.choice[,c(1:17)]

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

#restrict analysis to allopatric populations only:
data.m.choice.sub.allo<-data.m.choice.sub[data.m.choice.sub$topy=="allo",]
#choice data for pops in both datasets:
test2<-unique(data.m.choice.sub.allo[,c(1,29:33)]) #7 pops, N=110 females
summary(test2$N.females)

#Now, see if preferences predict call rates.  
choice.mod1<-lmer(call.rate.tempcor ~ car::logit(prop.chose.f)+elevation + year + (1|pop), 
                   weights = N.females, data = data.m.choice.sub.allo)
#singular fit---model is overfitted.  
library(rstanarm)
set.seed(42069)
choice.mod1<-stan_lmer(call.rate.tempcor ~ car::logit(prop.chose.f)+elevation + year + (1|pop), 
                   weights = N.females, data = data.m.choice.sub.allo, adapt_delta = 0.99)
#specifying adapt_delta= 0.99 above because using default of 0.8 results in 2 divergent transitions after warmup.

summary(choice.mod1, pars = c("car::logit(prop.chose.f)", "elevation", "year"), 
        probs = c(0.025, 0.975), digits = 2)

choice.mod2<-stan_lmer(call.rate.tempcor ~ car::logit(prop.chose.f) + year + (1|pop), 
                   weights = N.females, data = data.m.choice.sub.allo)

summary(choice.mod2, pars = c("car::logit(prop.chose.f)", "year"), 
        probs = c(0.025, 0.975), digits = 2)

