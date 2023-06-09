#Contact author E.H. Schultheis (schulth5@msu.edu) for any questions on analysis
#Autor A.E. Berardi added phylogenetic analysis 2.15.15

######-----IMPORT DATA------######

damagedata<-read.csv("damage.csv",header=T) #import data
damagedata$collectionyear =as.factor(damagedata$collectionyear)
summary(damagedata)
names(damagedata)
str(damagedata)

damagedata.na <- na.omit(damagedata) #remove rows with missing data
str(damagedata.na)

###----------------------------------------------------------------------------###
###-----------Means & SE for Figure 1, Figure 2, Supplemental Figure 1---------###
###----------------------------------------------------------------------------###
library(vegan)
library(reshape)

se <- function(x) {
  s <- sd(x, na.rm=TRUE)
  n=length(x)
  s/sqrt(n)
} #function defining standard error for melt and cast

se.prop <- function(x) {
  p=(sum(x))/length(x)
  sqrt(p*(1-p)/length(x))
} #function defining standard error for BINOMIAL (y/n) data for melt and cast

#---------melt and cast by species---------#
reshape <- melt(damagedata.na, id=c("plot", "species", "name", "family", "type", "status", "collectionyear"), measured=c("propleaf", "propbranchesbrowsed"))
str(reshape) #take a look at the data

#cast by all predictor variables to get means and se for species
species.means=cast(reshape, species + family + type + status + collectionyear ~ variable, fun.aggregate=mean, na.rm=FALSE)
species.se=cast(reshape, species + family + type + status + collectionyear ~ variable, fun.aggregate=se)
str(species.means)
str(species.se)

######-----CREATE NEW DATA FILE------######

fix(species.means) #make a new data file that is the mean for each species separated out by each year. This data file will be used for all ANOVA analysis.
write.table(species.means, "speciesmean.csv", quote=FALSE, sep=",", row.names=FALSE)

#---------re-melt using species means from previous cast---------#
year.status.melt <- melt(species.means, id=c("species", "family", "type", "status", "collectionyear"), measured=c("propleaf", "propbranchesbrowsed")) 
str(year.status.melt)

#cast by status and year
#mean and standard error data used for Figure 1
year.means1=cast(year.status.melt, status + collectionyear ~ variable, fun.aggregate=mean, na.rm=TRUE)
year.se1=cast(year.status.melt, status + collectionyear ~ variable, fun.aggregate=se)
fix(year.means1) #fix lets you see the actual data to extract for figures
fix(year.se1)

#cast by family and year
#mean and standard error data used for Supplemental Figure 1
year.means2=cast(year.status.melt, family + collectionyear ~ variable, fun.aggregate=mean, na.rm=TRUE)
year.se2=cast(year.status.melt, family + collectionyear ~ variable, fun.aggregate=se)
fix(year.means2)
fix(year.se2)

#cast by family, status, and year
#mean and standard error data used to show interactions between status and family for Figure 2
year.means3=cast(year.status.melt, family + status + collectionyear ~ variable, fun.aggregate=mean, na.rm=TRUE)
year.se3=cast(year.status.melt, family + status + collectionyear ~ variable, fun.aggregate=se)
fix(year.means3)
fix(year.se3)

######-----CREATE NEW DATA FILE------######

#cast by species only, data to be used for all logistic ANCOVA analysis
time.means=cast(year.status.melt, species + status + family ~ variable, fun.aggregate=mean, na.rm=TRUE)
fix(time.means) #each species averaged across all three years
write.table(time.means, "timespeciesmean_whole.csv", quote=FALSE, sep=",", row.names=FALSE) #There are a few key differences between this file that is created, and the "timespeciesmean.csv" file that is used later on in all logistic ANCOVA analysis. First, native species are removed because the analysis will just be performed on introduced species (exotic and invasive). Second, the year introduced to Michigan (year), residence time (time), and number of counties spread (county) data columns will be added to look at whether enemy damage is dynamic with time or spread.

###----------------------------------------------------------------------------###
###-ANOVA & Tukey tests for Table 1, Figure 1, Figure 2, Supplemental Figure 1-###
###----------------------------------------------------------------------------###
#Table 1 is now a supplemental table
library(multcomp) #for tukey contrasts
library(agricolae) #for tukey HSD.test
library(car) #use Anova function with lm to add type III error

######-----IMPORT DATA------######

speciesmeans<-read.csv("speciesmean.csv",header=T) #import new data file, the one created from melting and casting earlier
speciesmeans$collectionyear = as.factor(speciesmeans$collectionyear)
summary(speciesmeans)
names(speciesmeans)
str(speciesmeans)

#subset species mean data by year
data.2011 <- subset(speciesmeans, collectionyear=="2011")
data.2012 <- subset(speciesmeans, collectionyear=="2012")
data.2013 <- subset(speciesmeans, collectionyear=="2013")
#Give list of species present in each year of data (to prune the tree)
summary(data.2011$species)
summary(data.2012$species)
summary(data.2013$species)
#rename variables: native = 1, exotic = 2, invasive = 3

#subset dataset by family for 2012 - for Tukey interaction analysis
#will no longer 
aster.2012 <- subset(data.2012, family=="Asteraceae")
fab.2012 <- subset(data.2012, family=="Fabaceae")
poa.2012 <- subset(data.2012, family=="Poaceae")

######-----HERBIVORY------######

#2011 Herbivory- can't do any interaction terms without getting singularities. There is no replication within family to make it work.
model_herbivory_2011 = aov(propleaf ~ status + family, data=data.2011)
summary(model_herbivory_2011)
Anova(model_herbivory_2011, type="III") #Type III error
posthoc <- HSD.test(model_herbivory_2011, "status", std.err) #Tukey test
posthoc

#2012 Herbivory
model_herbivory_2012 = aov(propleaf ~ status * family, data=data.2012)
summary(model_herbivory_2012)
Anova(model_herbivory_2012, type="III") #Type III error
posthoc <- HSD.test(model_herbivory_2012, "status", std.err) #Tukey test
posthoc
posthoc <- HSD.test(model_herbivory_2012, "family", std.err) #Tukey test
posthoc

#to get Tukey letters for the status x family interaction in 2012 for herbivory
model.aster.2012 = aov(propleaf ~ status, data=aster.2012) #Asteraceae
posthoc <- HSD.test(model.aster.2012, "status", std.err)
posthoc
model.fab.2012 = aov(propleaf ~ status, data=fab.2012) #Fabaceae
posthoc2 <- HSD.test(model.fab.2012, "status", std.err)
posthoc2
model.poa.2012 = aov(propleaf ~ status, data=poa.2012) #Poaceae
posthoc3 <- HSD.test(model.poa.2012, "status", std.err)
posthoc3

#2013 Herbivory
model_herbivory_2013 = aov(propleaf ~ status * family, data=data.2013)
summary(model_herbivory_2013)
Anova(model_herbivory_2013, type="III") #Interaction term not significant for this year, so take it out of the model
model_herbivory_2013b = aov(propleaf ~ status + family, data=data.2013)
summary(model_herbivory_2013b)
Anova(model_herbivory_2013b, type="III")
posthoc <- HSD.test(model_herbivory_2013b, "family", std.err)
posthoc

######-----BROWSING-----######

#2011 Browsing - can't do any interaction terms without getting singularities. There is no replication within family to make it work.
model_browsing_2011 = aov(propbranchesbrowsed ~ status + family, data=data.2011)
summary(model_browsing_2011)
Anova(model_browsing_2011, type="III") #Type III error

#2012 Browsing
model_browsing_2012 = aov(propbranchesbrowsed ~ status * family, data=data.2012)
summary(model_browsing_2012)
Anova(model_browsing_2012, type="III") #Interaction term not significant for this year, so take it out of the model
model_browsing_2012b = aov(propbranchesbrowsed ~ status + family, data=data.2012)
summary(model_browsing_2012b)
Anova(model_browsing_2012b, type="III")

#2013 Browsing
model_browsing_2013 = aov(propbranchesbrowsed ~ status * family, data=data.2013)
summary(model_browsing_2013)
Anova(model_browsing_2013, type="III") #Interaction term not significant for this year, so take it out of the model
model_browsing_2013b = aov(propbranchesbrowsed ~ status + family, data=data.2013)
summary(model_browsing_2013b)
Anova(model_browsing_2013b, type="III") 
posthoc <- HSD.test(model_browsing_2013b, "family") #Tukey test for Family
posthoc

###----------------------------------------------------------------------------###
###-----------Logistic ANCOVA analysis for Table 2 and Figures 3 & 4-----------###
###----------------------------------------------------------------------------###
#Table 2 is now a supplemental table
library(multcomp) #for tukey contrasts
library(agricolae) #for tukey
library(car) 

######-----IMPORT DATA------######

timespeciesmean<-read.csv("timespeciesmean.csv",header=T) #This data is the herbivory and browsing for each species, averaged across years. Just depicts exotic and invasive species because natives do not have a date of introduciton or area of spread.
#Just like the "timespeciesmean_whole.csv" data file created by melting and casting earlier, except native species have been removed and time and spread data added as new data columns.
summary(timespeciesmean)
names(timespeciesmean)
str(timespeciesmean)

#get keep list for species present in ANCOVA data to prune tree
summary(timespeciesmean$species)

#subset data into exotic and invasive species seperate
data.exotic <- subset(timespeciesmean, status=="exotic")
data.invasive <- subset(timespeciesmean, status=="invasive")

#relationships are non-linear
#use glm and logit
#http://www.umass.edu/landeco/teaching/ecodata/schedule/methods1.pdf (page 20)

######-----TEST CORRELATION BETWEEN COUNTY & TIME------######

model_x <-lm(county ~  time, data=timespeciesmean) #r squared is around 0.5, pretty tight correlation between variables. If the two explanatory variables (time and county) are to be interpreted seperately in a multiple regression, there must be no collinearity or correlation between any of the independent variables. Solution is to analyze each seperate (Underwood p. 439, Experiments in Ecology). 
Anova(model_x)
summary(model_x) #(r = 0.70, P < 0.001)

model_x1 <-lm(county ~  uscounty, data=timespeciesmean)
Anova(model_x1)
summary(model_x1) #(r = 0.72, P < 0.001), spread in MI correlated with spread in the US

model_x2 <-lm(county ~  statecounty, data=timespeciesmean)
Anova(model_x2)
summary(model_x2) #(r = 0.86, P < 0.001), spread in MI correlated with spread in nearby 5 states

######-----HERBIVORY & MI SPREAD------###### Suplemental Table 2, Figure 4

model_herbivory_spread <-glm(herbivory ~  county * status, family=quasibinomial(link=logit), data=timespeciesmean) #County is the predictor variable for number of counties occupied by an introduced species. Status is the predictor variable that indicates whether an introduced species is exotic or invasive.
summary(model_herbivory_spread)
anova(model_herbivory_spread, test = "Chisq") #county and interaction significant

#Is the relationship significant for exotics? YES
model_ca <-glm(herbivory ~  county, family=quasibinomial(link=logit), data=data.exotic)
summary(model_ca)
anova(model_ca, test="Chisq") #exotic significant

# pseudo R^2 in the book Extending the Linear Model with R, Julian J. Faraway (p. 59).
#http://stats.stackexchange.com/questions/11676/pseudo-r-squared-formula-for-glms
1 - (model_ca$deviance/model_ca$null.deviance)

#Is the relationship significant for invasives? NO
model_caa <-glm(herbivory ~  county, family=quasibinomial(link=logit), data= data.invasive)
summary(model_caa)
anova(model_caa, test="Chisq") #invasive not significant

#Figure 4a
plot(timespeciesmean$herbivory ~ timespeciesmean$county, pch=(1), col="black", xlab="Number of MI Counties", ylab="Proportion Leaf Herbivory", main="")
points(timespeciesmean$herbivory[timespeciesmean$status =="exotic"] ~ timespeciesmean $county[timespeciesmean$status=="exotic"], pch=(16), col = "gray")
points(timespeciesmean $herbivory[timespeciesmean $status =="invasive"] ~ timespeciesmean $county[timespeciesmean $status=="invasive"], pch=(16), col = "black")
model_a2e <-glm(herbivory ~  county, family=quasibinomial(link=logit), data= data.exotic) #need to take out status or the curve function will not work
curve(predict(model_a2e,data.frame(county =x),type="resp"),add=TRUE, lwd=3, lty=1, col = "gray") 
legend("topleft", c("Exotic", "Invasive"), pch = 21, col="black", pt.bg=c("gray","black"), bg="white")

######----HERBIVORY & MI, WI, IL, IN, OH SPREAD-----###### Suplemental Table 2, Figure 4b

model_herbivory_spread2 <-glm(herbivory ~  statecounty * status, family=quasibinomial(link=logit), data=timespeciesmean) #County is the predictor variable for number of counties occupied by an introduced species. Status is the predictor variable that indicates whether an introduced species is exotic or invasive.
summary(model_herbivory_spread2)
anova(model_herbivory_spread2, test = "Chisq") 

# pseudo R^2 
1 - (model_herbivory_spread2$deviance/model_herbivory_spread2$null.deviance)

#Potential figure of 5 state county spread and herbivory
plot(timespeciesmean$herbivory ~ timespeciesmean$statecounty, pch=(1), col="black", xlab="Five State Counties", ylab="Proportion Leaf Herbivory", main="")
points(timespeciesmean$herbivory[timespeciesmean$status =="exotic"] ~ timespeciesmean $statecounty[timespeciesmean$status=="exotic"], pch=(16), col = "gray")
points(timespeciesmean $herbivory[timespeciesmean $status =="invasive"] ~ timespeciesmean $statecounty[timespeciesmean $status=="invasive"], pch=(16), col = "black")
model_a2b <-glm(herbivory ~  statecounty, family=quasibinomial(link=logit), data= timespeciesmean) #need to take out status or the curve function will not work
curve(predict(model_a2b,data.frame(statecounty =x),type="resp"), add=TRUE, lwd=3, lty=1, col="gray") #logistic regression line
curve(predict(model_a2b,data.frame(statecounty =x),type="resp"), add=TRUE, lwd=3, lty=2, col="black") #logistic regression line
legend("topleft", c("Exotic", "Invasive"), pch = 21, col="black", pt.bg=c("gray","black"), bg="white")

######-----HERBIVORY & US SPREAD------###### Suplemental Table 2, Figure 4c

model_herbivory_usspread <-glm(herbivory ~  uscounty * status, family=quasibinomial(link=logit), data=timespeciesmean) #County is the predictor variable for number of counties occupied by an introduced species. Status is the predictor variable that indicates whether an introduced species is exotic or invasive.
summary(model_herbivory_usspread)
anova(model_herbivory_usspread, test = "Chisq") 
#county is highly significant, p < 0.001. No interaction.
#Same pattern as MI data, just much stronger. With increasing US counties, proportion leaf herbivory increases

# pseudo R^2 
1 - (model_herbivory_usspread$deviance/model_herbivory_usspread$null.deviance)

#Potential figure of US county spread and herbivory
plot(timespeciesmean$herbivory ~ timespeciesmean$uscounty, pch=(1), col="black", xlab="Number of US Counties", ylab="Proportion Leaf Herbivory", main="")
points(timespeciesmean$herbivory[timespeciesmean$status =="exotic"] ~ timespeciesmean $uscounty[timespeciesmean$status=="exotic"], pch=(16), col = "gray")
points(timespeciesmean $herbivory[timespeciesmean $status =="invasive"] ~ timespeciesmean $uscounty[timespeciesmean $status=="invasive"], pch=(16), col = "black")
model_a2b <-glm(herbivory ~  uscounty, family=quasibinomial(link=logit), data= timespeciesmean) #need to take out status or the curve function will not work
curve(predict(model_a2b,data.frame(uscounty=x),type="resp"), add=TRUE, lwd=3, lty=1, col="gray") #logistic regression line
curve(predict(model_a2b,data.frame(uscounty=x),type="resp"), add=TRUE, lwd=3, lty=2, col="black") #logistic regression line
legend("topleft", c("Exotic", "Invasive"), pch = 21, col="black", pt.bg=c("gray","black"), bg="white")

######-----HERBIVORY & RESIDENCE TIME------###### Supplemental Table 2, Figure 3a

model_a22 <-glm(herbivory ~  time * status, family=quasibinomial(link=logit), data= timespeciesmean) #Time is the predictor variable for residence time, or number of years a species has been in Michigan, according to herbarium records. Status is the predictor variable that indicates whether an introduced species is exotic or invasive.
summary(model_a22)
anova(model_a22, test = "Chisq")

#Is the relationship significant for exotics? YES
model_a3a <-glm(herbivory ~  time, family=quasibinomial(link=logit), data=data.exotic)
summary(model_a3a)
anova(model_a3a, test="Chisq") #exotic significant

# pseudo R^2 
1 - (model_a3a$deviance/model_a3a$null.deviance)

#Is the relationship significant for invasives? NO
model_a3aa <-glm(herbivory ~  time, family=quasibinomial(link=logit), data= data.invasive)
summary(model_a3aa)
anova(model_a3aa, test="Chisq") #invasive not significant

#Figure 3a
plot(timespeciesmean$herbivory ~ timespeciesmean$time, pch=(1), col="black", xlab="Number of Years in MI", ylab="Proportion Leaf Herbivory", main="")
points(timespeciesmean$herbivory[timespeciesmean$status =="exotic"] ~ timespeciesmean $time[timespeciesmean$status=="exotic"], pch=(16), col = "gray")
points(timespeciesmean $herbivory[timespeciesmean $status =="invasive"] ~ timespeciesmean $time[timespeciesmean $status=="invasive"], pch=(16), col = "black")
model_a2e <-glm(herbivory ~  time, family=quasibinomial(link=logit), data= data.exotic) #need to take out status or the curve function will not work
curve(predict(model_a2e,data.frame(time=x),type="resp"),add=TRUE, lwd=3, lty=1, col = "gray") 
legend("topleft", c("Exotic", "Invasive"), pch = 21, col="black", pt.bg=c("gray","black"), bg="white")

######-----BROWSING & MI SPREAD------###### Supplemental Table 2, Figure 3c

model_a3 <-glm(browsing ~  county * status, family=quasibinomial(link=logit), data= timespeciesmean)
summary(model_a3)
anova(model_a3, test="Chisq")

#Is the relationship significant for exotics? YES
model_a3a <-glm(browsing ~  county, family=quasibinomial(link=logit), data=data.exotic)
summary(model_a3a)
anova(model_a3a, test="Chisq") #exotic significant

# pseudo R^2 
1 - (model_a3a$deviance/model_a3a$null.deviance)

#Is the relationship significant for invasives? NO
model_a3aa <-glm(browsing ~  county, family=quasibinomial(link=logit), data= data.invasive)
summary(model_a3aa)
anova(model_a3aa, test="Chisq") #invasive not significant

#Figure 4c
plot(timespeciesmean$browsing ~ timespeciesmean$county, pch=(1), col="black", xlab="Number of MI Counties", ylab="Proportion Branches Browsed", main="")
points(data.exotic $browsing ~ data.exotic $county, pch=(16), col = "gray")
model_a2e <-glm(browsing ~  county, family=quasibinomial(link=logit), data= data.exotic) #need to take out status or the curve function will not work
curve(predict(model_a2e,data.frame(county=x),type="resp"),add=TRUE, lwd=3, lty=1, col = "gray") 
points(timespeciesmean $browsing[timespeciesmean $status =="invasive"] ~ timespeciesmean $county[timespeciesmean $status=="invasive"], pch=(16), col = "black")
legend("top", c("Exotic", "Invasive"), pch = 21, col="black", pt.bg=c("gray","black"), bg="white")

######-----BROWSING & MI, WI, IL, IN, OH SPREAD------######

model_browsing_spread2 <-glm(browsing ~  statecounty * status, family=quasibinomial(link=logit), data=timespeciesmean) #County is the predictor variable for number of counties occupied by an introduced species. Status is the predictor variable that indicates whether an introduced species is exotic or invasive.
summary(model_browsing_spread2)
anova(model_browsing_spread2, test = "Chisq") 
# spread marginally significant (p=0.099), and interaction significant (p<0.001)

# pseudo R^2 
1 - (model_browsing_spread2$deviance/model_browsing_spread2$null.deviance)

#Is the relationship significant for exotics? YES
model_eb <-glm(browsing ~  statecounty, family=quasibinomial(link=logit), data=data.exotic)
summary(model_eb)
anova(model_eb, test="Chisq") #exotic significant

#Is the relationship significant for invasives? YES
model_ib <-glm(browsing ~  statecounty, family=quasibinomial(link=logit), data= data.invasive)
summary(model_ib)
anova(model_ib, test="Chisq") #invasive significant

#Potential figure of 5 state county spread and herbivory
plot(timespeciesmean$browsing ~ timespeciesmean$statecounty, pch=(1), col="black", xlab="Number of Counties in 5 States", ylab="Proportion Branches Browsed", main="")
points(timespeciesmean$browsing[timespeciesmean$status =="exotic"] ~ timespeciesmean $statecounty[timespeciesmean$status=="exotic"], pch=(16), col = "gray")
points(timespeciesmean $browsing[timespeciesmean $status =="invasive"] ~ timespeciesmean $statecounty[timespeciesmean $status=="invasive"], pch=(16), col = "black")
model_1 <-glm(browsing ~  statecounty, family=quasibinomial(link=logit), data= data.exotic) #need to take out status or the curve function will not work
curve(predict(model_1,data.frame(statecounty =x),type="resp"),add=TRUE, lwd=3, lty=1, col = "gray") 
model_2 <-glm(browsing ~  statecounty, family=quasibinomial(link=logit), data= data.invasive) #need to take out status or the curve function will not work
curve(predict(model_2,data.frame(statecounty=x),type="resp"),add=TRUE, lwd=3, lty=1, col = "black") 
legend("top", c("Exotic", "Invasive"), pch = 21, col="black", pt.bg=c("gray","black"), bg="white")

######-----BROWSING & US SPREAD------######

model_usa3 <-glm(browsing ~  uscounty * status, family=quasibinomial(link=logit), data= timespeciesmean)
summary(model_usa3)
anova(model_usa3, test="Chisq") #interaction significant

# pseudo R^2 
1 - (model_usa3$deviance/model_usa3$null.deviance)

#Is the relationship significant for exotics? YES
model_eb <-glm(browsing ~  uscounty, family=quasibinomial(link=logit), data=data.exotic)
summary(model_eb)
anova(model_eb, test="Chisq") #exotic significant

#Is the relationship significant for invasives? YES
model_ib <-glm(browsing ~  uscounty, family=quasibinomial(link=logit), data= data.invasive)
summary(model_ib)
anova(model_ib, test="Chisq") #invasive significant

#Potential figure showing relationship between browsing and spread in the US
plot(timespeciesmean$browsing ~ timespeciesmean$uscounty, pch=(1), col="black", xlab="Number of US Counties", ylab="Proportion Branches Browsed", main="")
points(data.exotic $browsing ~ data.exotic$uscounty, pch=(16), col = "gray")
points(timespeciesmean $browsing[timespeciesmean $status =="invasive"] ~ timespeciesmean $uscounty[timespeciesmean $status=="invasive"], pch=(16), col = "black")
model_1 <-glm(browsing ~  uscounty, family=quasibinomial(link=logit), data= data.exotic)
curve(predict(model_1,data.frame(uscounty =x),type="resp"),add=TRUE, lwd=3, lty=1, col = "gray") 
model_2 <-glm(browsing ~  uscounty, family=quasibinomial(link=logit), data= data.invasive) 
curve(predict(model_2,data.frame(uscounty =x),type="resp"),add=TRUE, lwd=3, lty=1, col = "black") 
legend("top", c("Exotic", "Invasive"), pch = 21, col="black", pt.bg=c("gray","black"), bg="white")

######-----BROWSING & TIME------###### Table 2, Figure 3d

model_a2 <-glm(browsing ~  time * status, family=quasibinomial(link=logit), data= timespeciesmean)
summary(model_a2)
anova(model_a2, test="Chisq")

# pseudo R^2 
1 - (model_a2$deviance/model_a2$null.deviance)

#Figure 3d
plot(timespeciesmean$browsing ~ timespeciesmean$time, pch=(1), col="black", xlab="Number of Years in MI", ylab="Proportion Branches Browsed", main="")
points(timespeciesmean$browsing[timespeciesmean$status =="exotic"] ~ timespeciesmean $time[timespeciesmean$status=="exotic"], pch=(16), col = "gray")
points(timespeciesmean $browsing[timespeciesmean $status =="invasive"] ~ timespeciesmean $time[timespeciesmean $status=="invasive"], pch=(16), col = "black")
legend("topright", c("Exotic", "Invasive"), pch = 21, col="black", pt.bg=c("gray","black"), bg="white")

###----------------------------------------------------------------------------###
###--------------------------OTHER STATISTICAL TESTS---------------------------###
###----------------------------------------------------------------------------###

######-----Is there a correlation between herbivory and browsing?------######

model_herbrow <-lm(herbivory ~  browsing, data=timespeciesmean)
summary(model_herbrow)
Anova(model_herbrow, type="III") #There was no relationship between the amount of damage a species received from insect herbivores and that received from mammalian browsers (R squared = 0.02, P = 0.46).

plot(timespeciesmean$herbivory ~ timespeciesmean$browsing, pch=(1), col="black", xlab="Proportion Branches Browsed", ylab="Proportion Leaf Herbivory", main="")
points(timespeciesmean$herbivory ~ timespeciesmean $browsing, pch=(16), col = "black")

######-----Supplemental Figure 2, boxplots showing data for each species------######

#Supplemental Figure 2a
speciesmeans$species_order<-factor(speciesmeans$species, levels=c("MELAL", "TRIRE", "CICIN","CORVA","CENST","MELOF","POTAR","EUPPE","TRIPR","TRIHY","AMCAN","MEDLU","HELAU","LOTCO","HESMA","DESCA","POTRE","DAUCA","POAPR","SOLCA","CONCA","CORTR","LESCU","SYMPI","CORLA","LESCA","HELFL","BROIN","AGRRE","BROHO","COSBI","GAIPU","POANE","PHLPR","POACO","BROKA","ACHMI","LEUVU","CACAR","CORTI","ERIAN","PANVI","SOLRI","POATR","SOLGR","SORNU","POTAG","CORPA","COSSU","TEPVI","SCHSC","CENCY","LACSA","SPOHE","TAROF")) #species are sorted in order of decending means, averaged across the 3 years. Boxplots show median and 1st and 3rd quartiles (explains why they don't look like they're perfectly descending)
boxplot(propleaf ~ species_order, data= speciesmeans, las=2, notch=FALSE, horizontal=FALSE,
 col=(c("gray16","ivory3","ivory3","gray16","gray16","gray16","ivory3","white","gray16","gray16","white","ivory3","ivory3","gray16","gray16","white","ivory3","ivory3","gray16","white","white","white","gray16","white","white","white","ivory3","gray16","gray16","ivory3","ivory3","ivory3","white","ivory3","gray16","white","white","ivory3","ivory3","ivory3","white","white","white","ivory3","white","white","white","white","ivory3","white","white","ivory3","ivory3","white","ivory3")),
  main="All Years Together", xlab="", ylab="Proportion Leaf Area Removed By Herbivory")
legend("topright", c("Native", "Exotic", "Invasive"), pch = 22, col="black", pt.bg=c("white","ivory3","gray16"), bg="white")

#Supplemental Figure 2b
speciesmeans$species_order<-factor(speciesmeans$species, levels=c("MELAL", "COSSU", "CENCY","COSBI","MELOF","GAIPU","SOLGR","CORTR","CENST","SOLCA","CORLA","EUPPE","LOTCO","SYMPI","CORPA","BROHO","LESCU","DESCA","ERIAN","POANE","HELAU","CORTI","TRIPR","LESCA","SOLRI","CACAR","MEDLU","DAUCA","POTAR","POTRE","SORNU","ACHMI","AGRRE","AMCAN","BROIN","BROKA","CICIN","CONCA","CORVA","HELFL","HESMA","LACSA","LEUVU","PANVI","PHLPR","POACO","POAPR","POATR","POTAG","SCHSC","SPOHE","TAROF","TEPVI","TRIHY","TRIRE")) #species are sorted in order of decending means, averaged across the 3 years. Boxplots show median and 1st and 3rd quartiles (explains why they don't look like they're perfectly descending)
boxplot(propbranchesbrowsed ~ species_order, data=speciesmeans, las=2, notch=FALSE, 
 col=(c("gray16","ivory3","ivory3","ivory3","gray16","ivory3","white","white","gray16","white","white","white","gray16","white","white","ivory3","gray16","white","white","white","ivory3","ivory3","ivory3","white","white","ivory3","ivory3","ivory3","ivory3","ivory3","white")),
  main="All Years Together", xlab="", ylab="Proportion Branches with Browsing Damage")
legend("topright", c("Native", "Exotic", "Invasive"), pch = 22, col="black", pt.bg=c("white", "ivory3","gray16"), bg="white")