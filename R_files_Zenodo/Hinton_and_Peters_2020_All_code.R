library(bipartite)
##read in Network data (plant by sp matrix .csv, EX: Full_network_plant_by_species_matrix)
network<-read.csv(file.choose(),row.names = 1, header=TRUE)

sum(network)#number of interactions/links
#Creates a matrix of the data for use in several functions including the modularity analysis
nw<-sortweb(network)#Sort web to minimize interaction crossings
web<-as.matrix(nw)


######## Observed Modularity Analysis ################
#Used DIRTLPAwb+ (Beckett 2016) modularity algorithm
res1<-metaComputeModules(web, method = "Beckett", deleteOriginalFiles = TRUE, N = 100)
res1@likelihood#Extract the observed modularity value (Q-obs)
##Plot the module web to show the number of modules, the module organization i.e. which pollinators and plants are grouped together in modules.
plotModuleWeb(res1, weighted = TRUE, displayAlabels = TRUE, displayBlabels = TRUE, labsize = 0.5, square.border = "white")
#weighted c- and z-values of each plant species from our observed. For plants:level="lower" and for pollinators: level="higher"
webcz1<-czvalues(res1, weighted = TRUE, level = "lower")
webcz1
#Weighted values are calculated based on species strength, not degree. See package details for more detail.

######## Null Model Modularity Comparison and Critical Thresholds #########

#Create 100 null models for a null model comparison. This output will be used to standardize the observed modularity (Z-score) and calculate the Adjusted Modularity (Q-adj) values for the observed networks (Dormann and Strauss 2014, Carstensen et al. 2016, Watts et al. 2016).
nulls <- nullmodel(web, N=100, method=3)
modules.nulls <- sapply(nulls, metaComputeModules)

#Z-score to assess if the network is significantly modularity; Z > 2 means the network is significantly modular (Carstensen et al. 2016, Saunders and Rader 2019)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
(z <- (res1@likelihood - mean(like.nulls))/sd(like.nulls))

#Adjusted modularity score, which considers network size and sampling intensity (Dormann and Strauss 2014, Watts et al. 2016)
Qadj<-res1@likelihood - mean(like.nulls)
Qadj

#Calculates the weighted c and z values for all previvously created null models to objectively define critical thresholds using 95 quantiles (Watts et al. 2016)
cz.nulls <- sapply(modules.nulls, czvalues, weighted=TRUE,level="lower")
c.nulls <- as.data.frame(unlist(cz.nulls[1,]))
colnames(c.nulls)[1] <- "cval"
c.crit <- quantile(c.nulls$cval,probs=c(0.95))#95% quantiles of your null c values which creates your c critical threshold
z.nulls <- as.data.frame(unlist(cz.nulls[2,]))
colnames(z.nulls)[1] <- "zval"
z.nulls <- na.omit(z.nulls)
z.crit <- quantile(z.nulls$zval,probs=c(0.95))#95% quantiles of your null z values which creates your z critical threshold

##Review critical thresholds derived from null models
c.crit
z.crit

######## Network Plot Codes#############
##Plot c z scores

cz<-read.csv(file.choose(), header =TRUE)#read in cz_values .csv, EX: Full_network_cz_values.csv
s<-subset(cz, c>0.65, z>1.79)#Subset using the critical thresholds, reference above values obtained from "c.crit" and "z.crit"
with(cz,
plot(cz$c,cz$z, 
     type = "p",
     pch = 16,
     cex = 2.5,
     col = c("#006600", "#00FF00")[cz$flowering],
     ylab = "z, within module participation",
     xlab = "c, among module connectivity",
     xlim = c(0, 1.0), ylim = c(-1.3,2.5)))
text(s[,3:4], labels=s[,1], pos=3, offset = 1)
abline(v=0.65) #c.crit
abline(h=1.79) #z.crit

legend("bottomright", 
       legend = c("Continuous", "Brief"), 
       col = c("#00FF00", "#006600"), 
       pch = 16, 
       pt.cex = 2,
       cex = 1.05,
       bty = "n", 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))


#NW plot
plot<-plotweb(web, method ="cca", text.rot = 90, arrow =  "down.center", col.interaction = "grey", col.high = "orange", col.low = c("#00FF00","#00FF00","#00FF00","#00FF00","#00FF00", "#00FF00", "#006600", "#006600", "#006600", "#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600"), y.lim = c(-1,2.5), bor.col.high = "orange", bor.col.low = c("#00FF00","#00FF00","#00FF00","#00FF00","#00FF00", "#00FF00", "#006600", "#006600", "#006600", "#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600","#006600"), bor.col.interaction = "grey")#color coded for the full network. Light colors = continuous phenology.

#species level attributes, used Full_network_plant_by_sp_matrix.csv data file 
specieslevel(web, index = c("normalised degree","betweenness","closeness","species strength","degree"), level="lower")


######## Sampling Completeness #########
#For further methodology descriptions, see Chacoff et al. 2012
library(vegan)

#Measuring sampling sufficiency in the pollinator community
#read in "Specpool_poll_sample_by_sp_matrix.csv" in the sampling sufficiency folder
poll<-read.csv(file.choose(), header=TRUE, row.names=1)

#rarefaction curves for accumulated pollinator species
poll.accum<-specaccum(poll, method = "rarefaction", conditioned =TRUE, gamma = "chao", w = NULL)

#estimated asymptotic richness of pollinator species
est.poll<-specpool(poll, smallsample = TRUE)
est.poll


#measuring sampling completeness for interactions, by sample
#read in "Specpool_interaction_sample_by_interaction_matrix.csv"
int<-read.csv(file.choose(), header=TRUE, row.names=1)

#rarefaction curve for accumulated interactions
int.accum<-specaccum(int, method = "rarefaction",conditioned =TRUE, gamma = "chao", w = NULL)

#estimated asymptotic richness of interactions, we used the chao estimator results
est.int<-specpool(int,smallsample=TRUE)
est.int


###Accumulation curve plots, with asymptotic richness ablines
#for pollinator richness (all plants pooled)
#observed richness = 171, Expected=247 +- 27.8
# ~69% of pollinator community captured
plot(poll.accum, add = FALSE, random = FALSE, ci = 2,
     ci.type = "polygon", col = "black", lty = 1, lwd=1,
     ci.col = "orange", ci.lty = 0,
     xvar = "individuals",
     ylim = c(0,300))
abline(h=c(219.2,247,274.8), col=c("grey50","black","grey50"), lty=c(2,1,2), lwd=c(1.5,2,1.5))
#for unique links/interactions (all plants pooled)
#Observed int richness=506, Expected=921 +- 73.3
# ~55% interactions captured
plot(int.accum, add = FALSE, random = FALSE, ci = 2,
     ci.type = "polygon", col = "black", lty = 1,
     ci.col = "#00FF00", ci.lty = 0,
     xvar = "individuals",
     ylim = c(0,1050))
abline(h=c(847.7,921,994.3), col=c("grey50","black","grey50"), lty=c(2,1,2), lwd=c(1.5,2,1.5))


#Measuring sampling completeness for each plant species separately
#read in "Specpool_env.csv"
env<-read.csv(file.choose(), header = TRUE, row.names = 1)
#estimated asymptotic richness by plant species
pool <- with(env, specpool(poll, Plant))
pool


####### GLMM and LMM Code, GGplot Code #######
##Packages
library(lme4)
library(Rmisc)
library(ggplot2)
library(MASS)
library(dplyr)
library(car)
library(gridExtra)
#The following analysis only considers shrubs with a continuous flowering phenology
#read in "Model_data_shrub_only.csv"
shrub<-read.csv(file.choose(), header = TRUE)

#Test for fit
attach(shrub)
plant1<-as.factor(plant)
site1<-as.factor(site)

model.poisson <- glmer(rich~log(florbund+0.5)+(1|site1), family=poisson)

list(residual.deviance=deviance(model.poisson), residual.degrees.of.freedom =df.residual(model.poisson), chisq.p.value=pchisq(deviance(model.poisson),df.residual(model.poisson),lower=F)) 

model.nbin <- glmer.nb(rich~log(florbund+0.5)+(1|site1), data=shrub)

list(residual.deviance=deviance(model.nbin), residual.degrees.of.freedom =df.residual(model.nbin), chisq.p.value=pchisq(deviance(model.nbin),df.residual(model.nbin),lower=F)) 

model.poisson <- glmer(rich~plant1+(1|site1), family=poisson)

list(residual.deviance=deviance(model.poisson), residual.degrees.of.freedom =df.residual(model.poisson), chisq.p.value=pchisq(deviance(model.poisson),df.residual(model.poisson),lower=F)) 

model.nbin <- glmer.nb(rich~plant1+(1|site1), data=shrub)

list(residual.deviance=deviance(model.nbin), residual.degrees.of.freedom =df.residual(model.nbin), chisq.p.value=pchisq(deviance(model.nbin),df.residual(model.nbin),lower=F)) 


#Testing pollinator richness, modelling pollinator richness as the repsonse variable, shrub species identity as the explanatory variable, and site as the random effect. 
model1<- glmer.nb(rich~plant+(1|site), data=shrub)  
model2<- glmer.nb(rich~1+(1|site), data=shrub)
anova(model1,model2)


library(multcomp)
###Tukey's posthoc comparison witha  single-step adjustement to identify which shrub species are statistically different from one another
posthoc <- glht(model1, linfct=mcp(plant="Tukey"))

mcs <- summary(posthoc, test=adjusted("single-step"))

cld(mcs, level=0.05, decreasing =TRUE)


#Testing pollinator abundance; modelling pollinator abundance as the repsonse variable, shrub species identity as the explanatory variable, and site as the random effect. 
model1<- lmer(log(abund+0.5)~plant+(1|site), data=shrub)  
model2<- lmer(log(abund+0.5)~1+(1|site), data=shrub)
anova(model1,model2)

###Tukey's posthoc comparison
posthoc <- glht(model1, linfct=mcp(plant="Tukey"))

mcs <- summary(posthoc, test=adjusted("single-step"))

cld(mcs, level=0.05, decreasing =TRUE)


##Testing pollinator richness; modelling pollinator richness as the repsonse variable, floral abundance as the explanatory variable, and site as the random effect. 
model1<- glmer.nb(rich~log(florbund+0.5)+(1|site), data=shrub)  
model2<- glmer.nb(rich~1+(1|site), data=shrub)
anova(model1,model2)


##Testing pollinator abundance; modelling pollinator abundance as the repsonse variable, floral abundance as the explanatory variable, and site as the random effect. 
model1<- lmer(log(abund+0.5)~log(florbund)+(1|site), data=shrub)  
model2<- lmer(log(abund+0.5)~1+(1|site), data=shrub)
anova(model1,model2)

##Plotting Code
library (multcomp)

#Average pollinator abundance per plant species
Data2<-summarySE(data=shrub, "abund", groupvars="plant", conf.interval=0.95)
Data2

Table1<-as.table(Data2$abund)

rownames(Table1)=Data2$plant

#Average pollinator species richness per plant species
Data3<- summarySE(data=shrub, "rich", groupvars="plant", conf.interval=0.95)
Data3

Table1<-as.table(Data3$rich)
Table1
rownames(Table1)=Data3$plant

#####Arranging more than one ggplot in a window
###Requires minor modifications to axis titles
p<-ggplot(Data2, aes(x=plant, y=abund, ymax=45,ymin=0.0))+geom_bar(stat="identity", fill="green2", colour="black", width=0.7)+geom_errorbar(aes(ymax=abund+se,ymin=abund-se), width=0.2, size=0.5, color="black")+ylab("Pollinator Abundance")+xlab("")

p2<-ggplot(Data3, aes(x=plant, y=rich, ymax=10.0,ymin=0.0))+geom_bar(stat="identity", fill="green2", colour="black", width=0.7)+geom_errorbar(aes(ymax=rich+se,ymin=rich-se), width=0.2, size=0.5, color="black")+xlab("Plant Species")+ylab("Pollinator Richness")

p3<-ggplot(shrub, aes(x=log(florbund+0.5), y=log(abund+0.5)))+ geom_point()+stat_smooth(method = lm, se=TRUE)+xlab("Pollinator Abundance (log)")+ylab("")

p4<-ggplot(shrub, aes(log(florbund+0.5), rich)) + geom_point()+stat_smooth(method = lm, se=TRUE)+xlab("Floral Abundance")+ylab("")


grid.arrange(p,p3,p2,p4)
#Final edits to figures were completed in Adobe Illustrator