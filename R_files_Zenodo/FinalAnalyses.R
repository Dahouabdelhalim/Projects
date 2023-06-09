##Read and check the data
data <- read.csv(file.choose(), na.strings=c(""))
summary(data)
head(data)

#Calculate sum of ep chicks and brood size by configuration
tapply(data$EGPChicks, data$Configuration, sum, na.rm=TRUE)
tapply(data$BroodSize, data$Configuration, sum, na.rm=TRUE)


###CONFIGURATION DIFFERENCES FOR TABLE

wilcox.test(data$RelatednessQ~data$Configuration)
tapply(data$RelatednessQ, data$Configuration, median, na.rm=TRUE)
tapply(data$RelatednessQ, data$Configuration, quantile, na.rm=TRUE)

wilcox.test(data$GroupSize~data$Configuration)
tapply(data$GroupSize, data$Configuration, median, na.rm=TRUE)
tapply(data$GroupSize, data$Configuration, quantile, na.rm=TRUE)

wilcox.test(data$NeighbourhoodDensity~data$Configuration)
tapply(data$NeighbourhoodDensity, data$Configuration, median, na.rm=TRUE)
tapply(data$NeighbourhoodDensity, data$Configuration, quantile, na.rm=TRUE)

wilcox.test(data$ImmNeigh~data$Configuration)
tapply(data$ImmNeigh, data$Configuration, median, na.rm=TRUE)
tapply(data$ImmNeigh, data$Configuration, quantile, na.rm=TRUE)

wilcox.test(data$AvgNestDist~data$Configuration)
tapply(data$AvgNestDist, data$Configuration, median, na.rm=TRUE)
tapply(data$AvgNestDist, data$Configuration, quantile, na.rm=TRUE)

wilcox.test(data$BroodSize~data$Configuration)
tapply(data$BroodSize, data$Configuration, median, na.rm=TRUE)
tapply(data$BroodSize, data$Configuration, quantile, na.rm=TRUE)

wilcox.test(data$MotherMass~data$Configuration)
tapply(data$MotherMass, data$Configuration, median, na.rm=TRUE)
tapply(data$MotherMass, data$Configuration, quantile, na.rm=TRUE)

wilcox.test(data$FatherMass~data$Configuration)
tapply(data$FatherMass, data$Configuration, median, na.rm=TRUE)
tapply(data$FatherMass, data$Configuration, quantile, na.rm=TRUE)



##POPULATION DIFFERENCES 
##GENETIC
hist(data$RelatednessQ)
wilcox.test(data$RelatednessQ~data$Population) 
tapply(data$RelatednessQ, data$Population, median, na.rm=TRUE)
tapply(data$RelatednessQ, data$Population, quantile, na.rm=TRUE)

##From site data - allelic richness
wilcox.test(sitedat$AllelicRichness~sitedat$Population)
tapply(sitedat$AllelicRichness, sitedat$Population, median, na.rm=TRUE)
tapply(sitedat$AllelicRichness, sitedat$Population, quantile, na.rm=TRUE)

##From site data - genetic diversity
wilcox.test(sitedat$GeneticDiversity~sitedat$Population)
tapply(sitedat$GeneticDiversity, sitedat$Population, median, na.rm=TRUE)
tapply(sitedat$GeneticDiversity, sitedat$Population, quantile, na.rm=TRUE)




#### 

###FINAL GLM MODEL OF EGP
m1full <- glm (formula = cbind(EGPChicks, BroodSize-EGPChicks) ~ Configuration+RelatednessQ+Population+PairGroup, family = binomial, data=data) 
summary(m1full)
par(mfrow=c(2,2))
plot(m1full)

#Check for population differences in configuration effect (interaction)
m1x <- glm (formula = cbind(EGPChicks, BroodSize-EGPChicks) ~ Configuration*Population+RelatednessQ+PairGroup, family = binomial, data=data) 
summary(m1x)
par(mfrow=c(2,2))
plot(m1x)

#Substitute Number of Neighbouring territories for Configuration
m1d <- glm (formula = cbind(EGPChicks, BroodSize-EGPChicks) ~ ImmNeighbours+RelatednessQ+Population+PairGroup, family = binomial, data=data) 
summary(m1d)
par(mfrow=c(2,2))
plot(m1d)

#Substitute Nearest nest distance for Configuration
m1n <- glm (formula = cbind(EGPChicks, BroodSize-EGPChicks) ~ AvgNestDist+RelatednessQ+Population+PairGroup, family = binomial, data=data) 
summary(m1n)
par(mfrow=c(2,2))
plot(m1n)

