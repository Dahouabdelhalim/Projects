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
#We know from work in review that call rates have been changing over years, want to control for 
#this variation.  First confirm that call rates change over time:
year.mod<-lmer("call.rate.tempcor ~ year + (1|pop)",
               data = data.m)
summary(year.mod) 
anova(year.mod, ddf = "Kenward-Roger")#sig increase with year on temp corrected calls

#Thus, we should include year as a fixed effect in subsequent models of call rate.  

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

##################
# Fig. 2##
##################
require(ggplot2)
require(ggbeeswarm)
fig2_start<-ggplot(data=data.m, 
             mapping = aes(x=topy, y=call.rate.tempcor))+
 geom_quasirandom(alpha = 0.4)+ #plots points with horizonta jitter to reduce overplotting while showing distribution info
  theme_classic()+
  scale_x_discrete(name="Population type", breaks=c("allo", "sym"), labels=c("allopatry", "sympatry"))+
  ylab(label = "Call rate (calls/min)")

#now get model predicted values using bootstrapping via bootMer
  #first make a prediction dataframe (supplies values of fixed effects for which to generate model predictions)
pred_df<-expand.grid(year = 1996:2018, pop = levels(data.m$pop))
get_topy <- function(x) {
  if (x %in% sym.pops) {
    return("sym")
  } else if (x %in% allo.pops) {
    return("allo")
  }
  return(NA)
}
pred_df$topy<-as.factor(vapply(X = pred_df$pop, FUN = get_topy, FUN.VALUE = character(1)))
set.seed(42069)
booted.sym<-bootMer(x = topy.mod, FUN = function(x)(predict(x, newdata = pred_df[pred_df$topy == "sym",], re.form = NULL)),
                nsim = 1000)
booted.allo<-bootMer(x = topy.mod, FUN = function(x)(predict(x, newdata = pred_df[pred_df$topy == "allo",], re.form = NULL)),
                    nsim = 1000)
booteddf1<-data.frame(boot.cr = booted.allo$t[,1], topy = "allo")
booteddf2<-data.frame(boot.cr = booted.sym$t[,1], topy = "sym")
booteddf<-rbind.data.frame(booteddf1, booteddf2)

fig2<-fig2_start + 
  geom_boxplot(data = booteddf, aes(x = topy, y = boot.cr), alpha = 0.5)

#svg output
svg(filename="C:/Users/Gina/Dropbox/Work/Manuscripts/calls_prefs_mismatch/revision/figures/fig2.svg", 
    width=6, 
    height=4, 
    pointsize=12, 
    antialias = c("subpixel")
)
fig2
dev.off()
 
######################################

## C.  Effect of male genotype on call rates 

##############################
#could sympatric calls be faster than allopatric due to introgression from fast-calling 
#S. bombifrons?

summary(data.m$male.genotype) #298 NAs, most males above 4.82.  

#Test whether pure-M (males scoring 5) and introgressed (males scoring from 4 to 5) have 
#different call rates:
data.pure.m<-data.m[(is.na(data.m$male.genotype)==F & data.m$male.genotype>4.99),]
data.introg.m<-data.m[(is.na(data.m$male.genotype)==F & data.m$male.genotype<4.99),]

t.test(x=data.pure.m$call.rate.tempcor, y=data.introg.m$call.rate.tempcor) #not sig diff 
#mean of pure-m call rate is 35.43, mean of introgressed is 35.21
sd(data.pure.m$call.rate.tempcor)/sqrt(length(data.pure.m$call.rate.tempcor))
sd(data.introg.m$call.rate.tempcor)/sqrt(length(data.introg.m$call.rate.tempcor))

gt.mod<-lmer(call.rate.tempcor ~ male.genotype + year + (1|pop), data = data.m)
summary(gt.mod)
anova(gt.mod, ddf = "Kenward-Roger")

######################################

## D.  Effect of Elevation on call rates 

######################################
year.elev.mod<-lmer(call.rate.tempcor ~ elevation + year + (1|pop), data = data.m)
summary(year.elev.mod)
anova(year.elev.mod, ddf = "Kenward-Roger")
#year, elev both sig

############
## fig1
###########
require(sjPlot)
fig1<-sjPlot::plot_model(model = year.elev.mod, #creates glmer prediction plot as a ggplot object
                         type = "pred",
                         terms = c("elevation"),  #plotting only 1 fixed effect at a time
                         pred.type = "fe",  #it's a fixed effect
                         ci.lvl = 0.95,
                         title = "")+
  theme_bw()+
  theme(axis.line = element_line(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())+

  ylab(label = "Call rate (calls/min)")+
  xlab(label = "Elevation (m)")+
  geom_point(data = data.m, mapping = aes(x=elevation, y=call.rate.tempcor), color = "grey34")


#svg output
svg(filename="C:/Users/Gina/Dropbox/Work/Manuscripts/calls_prefs_mismatch/revision/figures/fig1.svg", 
    width=6, 
    height=4, 
    pointsize=12, 
    antialias = c("subpixel")
)
fig1
dev.off()

######################################

## E.  Does Female preference predict male call rate across populations?

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
#of females choosing fast in each population 

#remove pops where we have no female choice data:
data.m.choice2<-droplevels(data.m.choice[is.na(data.m.choice$N.females)==F,])

#restrict analysis to populations with at least 4 females sampled
data.m.choice.sub<-data.m.choice2[data.m.choice2$N.females>3,]

#choice data for pops in both datasets:
test2<-unique(data.m.choice.sub[,c(1,29:33)]) #14 pops, N=169 females
summary(test2$N.females)

#test for female preferences in sympatry and allopatry:
sym<-test2[test2$topy=="sym",]
allo<-test2[test2$topy=="allo",]
sum(sym$N.chose.f) #30
sum(sym$N.chose.s) #29
#test for preferences in sympatry using exact binomial test:
binom.test(x = c(30,29)) #no sig pref in sypatry: p = 1
#test for preferences in sympatry using exact binomial test:
sum(allo$N.chose.f) #70
sum(allo$N.chose.s) #40
binom.test(x = c(70,40)) #sig pref for fast in allopatry, p = 0.005

#test whether sample size predicts binomial outcomes:
choice.sample.mod<-glm(formula = cbind(N.chose.f, N.chose.s) ~ N.females, data = test2, family = binomial)
summary(choice.sample.mod)
anova(choice.sample.mod, test = "Chisq")

library(psych)
print(corr.test(x=test2$prop.chose.f, y=test2$N.females, method = "spearman"), short = FALSE)
#spearman correlation coefficient: r = 0.12, p = 0.7


#Now, see if preferences predict call rates.  
choice.mod1<-lmer(call.rate.tempcor ~ car::logit(prop.chose.f) + year + elevation + (1|pop), weights = N.females, data = data.m.choice.sub)
choice.mod2<-lmer(call.rate.tempcor ~ car::logit(prop.chose.f) + year + (1|pop), weights = N.females, data = data.m.choice.sub)

anova(choice.mod1, ddf = "Kenward-Roger") #no effect of female preference
anova(choice.mod2, ddf = "Kenward-Roger") #no effect of female preference

summary(choice.mod1)
summary(choice.mod2)

#model for plotting:
data.m.choice.sub$logitpref<-car::logit(data.m.choice.sub$prop.chose.f)
choice.mod.3<-lmer(call.rate.tempcor ~ logitpref + year + (1|pop), weights = N.females, data = data.m.choice.sub)


#############
## fig 3
#############
fig3<-sjPlot::plot_model(model = choice.mod.3, #creates glmer prediction plot as a ggplot object
                         type = "pred",
                         terms = c("logitpref"),  #plotting only 1 fixed effect at a time
                         pred.type = "fe",  #it's a fixed effect
                         ci.lvl = 0.95,
                         title = "")+
  theme_classic()+
  scale_y_continuous(name = "Male call rate (calls/min)", breaks = c(20,30,40,50), labels = c("20", "30", "40", "50"))+
  expand_limits(y=c(20,50))+
  xlab(label = "Log odds of female preference for fast calls")+
  geom_point(data = data.m.choice.sub, mapping = aes(x=logitpref, y=call.rate.tempcor), color = "grey34")

#svg output
svg(filename="C:/Users/Gina/Dropbox/Work/Manuscripts/calls_prefs_mismatch/revision/figures/fig3.svg", 
    width=6, 
    height=4, 
    pointsize=12, 
    antialias = c("subpixel")
)
fig3
dev.off()

################################################

###Fig S1####
data.m.gtknown<-data.m[is.na(data.m$male.genotype)==F,]

summary(data.m.gtknown$topy) #74 allo, 154 sym
summary(data.m$topy) #284 allo, 242 sym

#Let's visualize the differences in datasets:
sym.gt<-data.m.gtknown[data.m.gtknown$topy=="sym",]
allo.gt<-data.m.gtknown[data.m.gtknown$topy=="allo",]
sym.all<-data.m[data.m$topy=="sym",]
allo.all<-data.m[data.m$topy=="allo",]
boxdf1<-data.frame(allo.all$call.rate.tempcor)
boxdf1$group<-"allo.all"
names(boxdf1)
boxdf2<-data.frame(allo.gt$call.rate.tempcor)
boxdf2$group<-"allo.gt"
boxdf3<-data.frame(sym.all$call.rate.tempcor)
boxdf3$group<-"sym.all"
boxdf4<-data.frame(sym.gt$call.rate.tempcor)
boxdf4$group<-"sym.gt"
names(boxdf1)[1]<-"call.rate.tempcor"
names(boxdf2)[1]<-"call.rate.tempcor"
names(boxdf3)[1]<-"call.rate.tempcor"
names(boxdf4)[1]<-"call.rate.tempcor"
boxdf<-rbind.data.frame(boxdf1, boxdf2, boxdf3, boxdf4)

namevec<-c("allopatry\\nall", "allopatry\\ngenotyped", "sympatry\\nall", "sympatry\\ngenotyped")

fig.s1a<-ggplot(data = boxdf, aes(x = group, y = call.rate.tempcor))+
  geom_boxplot()+
  theme_classic()+
  labs(y = "Call rate (calls/min)", x = "Dataset")+
  scale_x_discrete(labels = namevec)

#plot overlapping histograms of call rates from data.m and data.m.gtknown
#First need to set up a data frame with both sets of data
histdf.all<-data.frame(data.m$call.rate.tempcor)
histdf.all$dataset<-"all males"
names(histdf.all)<-c("call.rate.tempcor", "dataset")
histdf.gtknown<-data.frame(data.m.gtknown$call.rate.tempcor)
histdf.gtknown$dataset<-"genotpyed males"
names(histdf.gtknown)<-c("call.rate.tempcor", "dataset")
histdf<-rbind.data.frame(histdf.all, histdf.gtknown)

#now plot a normalized overlapping histogram:
fig.s1b<-ggplot(histdf, aes (x = call.rate.tempcor, fill = dataset))+
  geom_histogram(aes(y=..density..), alpha = 0.5, position = "identity", binwidth = 1)+ #makes it a normalized (proportional) histogram
  theme_classic()+
  labs(x = "Call rate (calls/min)", y = "Normalized count")+
  theme(legend.position = "top", legend.title = element_blank())

fig.s1<-ggarrange(fig.s1b, fig.s1a, nrow = 1, ncol = 2, labels = "auto")

#export figure
#svg output
svg(filename="C:/Users/Gina/Dropbox/Work/Manuscripts/calls_prefs_mismatch/revision/figures/figS1.svg", 
    width=6, 
    height=4, 
    pointsize=12, 
    antialias = c("subpixel")
)
fig.s1
dev.off()







######## meta questions for methods etc.###############
#sampling years/ages of females used for preference data?
#first trim down preference dataset to only those females used in this analysis
#(those females from populations where we also had male call data)

calls.pops<-unique(data.m.choice.sub$pop) #pops used in calls-prefs analysis
females.used<-df.choice[df.choice$pop%in%calls.pops,] #index choice dataset by those pops

#get collection range of females:
summary(females.used$yr_collected) #2008-2017
#get years of testing range of females.  First convert date of testing and collection
#into R date format.
females.used$date_hi<-as.Date(females.used$date_hi, format = "%Y%m%d")
females.used$yr_tested<-format(females.used$date_hi, "%Y")
females.used$yr_collected<-as.character(females.used$yr_collected)
#calculate how long they had been in the colony at testing (our best proxy for age...)
females.used$age.at.test<-as.numeric(females.used$yr_tested)-as.numeric(females.used$yr_collected)

#test whether female age or testing year predicts preference for individual females:
age.indiv.mod<-glmer(formula = mate_choice_hi ~ age.at.test + (1|pop), 
                   family = binomial, data = females.used, na.action = na.fail)
null.indiv.mod<-glmer(formula = mate_choice_hi ~ (1|pop),
                      family = binomial, data = females.used)
library(lmerTest)
summary(age.indiv.mod)
anova(age.indiv.mod, null.indiv.mod) #no sig effect of age at testing

tstyr.indiv.mod<-glmer(formula = mate_choice_hi ~ yr_tested + (1|pop),
                        family = binomial, data = females.used)
summary(tstyr.indiv.mod)
anova(tstyr.indiv.mod, null.indiv.mod) #no sig effect of year of testing

