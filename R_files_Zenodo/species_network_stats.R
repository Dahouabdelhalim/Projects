setwd("C:/Users/llf44/OneDrive/Documents/ACADEMIA/Graduate School/R Files/Cornell/bee_pathogens/2015-FINAL/science/verified")
spnetR<-read.csv("beenetspeciesR.csv")
spnetUR<-read.csv("beenetspeciesUR.csv")

library(lme4)
library(car)

#Unresolved network
spnetUR$total<-spnetUR$tryp+spnetUR$neo+spnetUR$N.c+spnetUR$N.b
spnetUR$gen<-pmin(spnetUR$total, 1)

# General pathogen prevalence
## first checking variance inflation, if less than two all factors can be included in the same model
vif(glmer(gen ~ mod.group + degree + betweenness + (1|site) + (1|species), data=spnetUR, family=binomial)) # mod.group really increases VIF
vif(glmer(gen ~  degree + betweenness + (1|site) + (1|species), data=spnetUR, family=binomial)) # looks good (<2), so will run complete model

gen.null<-glmer(gen ~ 1 + (1|site) + (1|species), data=spnetUR, family=binomial)

##module group
genmod.group<-glmer(gen ~ mod.group + (1|site) + (1|species), data=spnetUR, family=binomial)
anova(genmod.group, gen.null)

#degree and betweenness centrality
genfull<-glmer(gen ~ degree + betweenness + (1|site) + (1|species), data=spnetUR, family=binomial)
anova(gen.null,genfull)

#Resolved network
spnetR$total<-spnetR$tryp+spnetR$neo+spnetR$N.c+spnetR$N.b
spnetR$gen<-pmin(spnetR$total, 1)

# General pathogen prevalence
## first checking variance inflation, if less than two all factors can be included in the same model
vif(glmer(gen ~ mod.group + degree + betweenness + (1|site) + (1|species), data=spnetR, family=binomial)) # mod.group really increases VIF
vif(glmer(gen ~  degree + betweenness + (1|site) + (1|species), data=spnetR, family=binomial)) # looks good (<2), so will run complete model

gen.null<-glmer(gen ~ 1 + (1|site) + (1|species), data=spnetR, family=binomial)

##module group
genmod.group<-glmer(gen ~ mod.group + (1|site) + (1|species), data=spnetR, family=binomial)
anova(genmod.group, gen.null)

#degree and betweenness centrality
genfull<-glmer(gen ~ degree + betweenness + (1|site) + (1|species), data=spnetR, family=binomial)
anova(gen.null,genfull)
