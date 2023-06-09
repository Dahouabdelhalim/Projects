##########################################################################
#
#	Script for meta-analyses of West Indian extinct mammal dates
#
#	Files required in same folder
#	data: fossil.csv archeo.csv
#	jags models: model.txt model0.txt model2.txt model3.txt model4.txt
#
##########################################################################
library(R2jags)
library(nlme)
rm(list = ls())

#read data
arch<-read.csv("archeo.csv", fileEncoding="latin1")
foss<-read.csv("fossil.csv", fileEncoding="latin1")
arc1<-subset(arch, (!is.na(arch$hage)))
arc1$lon<-arc1$lon*-1
arc2<-aggregate(hage ~ Island, data = arc1, max)
arc3 <- merge(arc1, arc2, by = c("Island", "hage"))
fos1<-subset(foss, Extinct=="Y")
fos1<-subset(fos1, synonym=="N")
fos2<-aggregate(age ~ island+Species, data = fos1, min)
fos3 <- merge(fos1, fos2, by = c("island", "Species", "age"), all.x=F, all.y=T)
fos3$island2<-ifelse(fos3$island=="Cuba", ifelse(fos3$age < 5000, "Cuba_a", "Cuba_b"), ifelse(fos3$island=="Hispaniola", ifelse(fos3$age < 5000, "Hispaniola_a", "Hispaniola_b"), as.character(fos3$island)))
alld<-merge(fos3, arc3, by.x="island", by.y="Island")
alld<-subset(alld, age<12000)
alld$diff<-alld$hage-alld$age
alld$island<-factor(alld$island)
alld$island2<-factor(alld$island2)

#hierarchical maximum likelihood models
m0<-lme(diff~1, alld, random = ~1|island, method="ML")
m1<-lme(diff~Direct.Date.of.mammal, alld, random = ~1|island, method="ML")
m2<-lme(diff~Direct.Date.of.mammal + Archaeo.or.Palaeo, alld, random = ~1|island, method="ML")
m2b<-lme(diff~Archaeo.or.Palaeo, alld, random = ~1|island, method="ML")
m3<-lme(diff~Direct.Date.of.mammal, alld, random = ~1|island2, method="ML")

#write out ML results
sink("ml_model_comparison.txt")
print(anova(m0, m1, m2, m2b, m3))
sink()
sink("ml_best_model.txt")
print("Summary and fixed coefficients")
print(summary(m3))
print("All coefficients")
print(coef(m3))
sink()

#write out coefficients from ML model
table_m<-as.data.frame(cbind(coef(m3)[,1], c(as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Abaco"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Anguilla"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Antigua"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Barbados"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Barbuda"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Cuba_a"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Cuba_b"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Desirade"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Grenada"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Guadeloupe"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Hispaniola_a"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Hispaniola_b"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island=="Jamaica"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island=="Marie Galante"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Martinique"))[1]), 
as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Montserrat"))[1]), 
as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Nevis"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island=="Puerto Rico"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="Saba"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="St. Eustatius"))[1]), as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="St. Kitts"))[1]), 
as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="St. Lucia"))[1]),
as.numeric(VarCorr(m3)[3])/sqrt(dim(subset(alld, island2=="St. Martin"))[1]))))
colnames(table_m)<-c("mean", "standard error")
rownames(table_m)<-unique(alld$island2)
table_m$min<-table_m[,1]-table_m[,2]
table_m$max<-table_m[,1]+table_m[,2]
write.csv(table_m, "ml_island_effect.csv")

#prepare data for Bayesian analyses
y<-alld$diff
m<-as.numeric(as.factor(alld$Direct.Date.of.mammal))
d<-as.numeric(as.factor(alld$Archaeo.or.Palaeo))
island<-as.factor(as.numeric(alld$island))
J<-length(levels(island))
N<-length(y)

#initial model, see .txt file for description
car0.data<-list("N", "y", "J", "island")
car0.parameters<-c("a", "sigma.y", "e.y", "mu.a", "sigma.a", "e.a")
car0.inits<-function(){list(a=rnorm(J), sigma.y=runif(1), mu.a=rnorm(1), sigma.a=runif(1))}
car.0<-jags(car0.data, car0.inits, car0.parameters, "model0.txt", n.chains=4, n.iter=5000)

#see .txt file for description
car.data<-list("N", "y", "J", "island", "m")
car.parameters<-c("a", "b0", "sigma.y", "e.y", "mu.a", "sigma.a", "e.a")
car.inits<-function(){list(a=rnorm(J), sigma.y=runif(1), mu.a=rnorm(1), sigma.a=runif(1), b0=rnorm(3) )}
car.1<-jags(car.data, car.inits, car.parameters, "model.txt", n.chains=4, n.iter=5000)

#see .txt file for description
car2.data<-list("N", "y", "J", "island", "m", "d")
car2.parameters<-c("a", "b0", "b1", "sigma.y", "e.y", "mu.a", "sigma.a", "e.a")
car2.inits<-function(){list(a=rnorm(J), sigma.y=runif(1), mu.a=rnorm(1), sigma.a=runif(1), b0=rnorm(3), b1=rnorm(3) )}
car.2<-jags(car2.data, car2.inits, car2.parameters, "model2.txt", n.chains=4, n.iter=5000)

#see .txt file for description
car3.data<-list("N", "y", "J", "island", "m", "d")
car3.parameters<-c("a", "b0", "b1", "sigma.y", "e.y", "mu.a", "sigma.a", "e.a")
car3.inits<-function(){list(a=rnorm(J), sigma.y=runif(1), mu.a=rnorm(1), sigma.a=runif(J), b0=rnorm(3), b1=rnorm(3) )}
car.3<-jags(car3.data, car3.inits, car3.parameters, "model3.txt", n.chains=4, n.iter=5000)
island<-as.factor(as.numeric(alld$island2))
J<-length(levels(island))

#see .txt file for description
car4.data<-list("N", "y", "J", "island", "m")
car4.parameters<-c("a", "b0", "sigma.y", "e.y", "mu.a", "sigma.a", "e.a")
car4.inits<-function(){list(a=rnorm(J), sigma.y=runif(1), mu.a=rnorm(1), sigma.a=runif(J), b0=rnorm(3))}
car.4<-jags(car4.data, car4.inits, car4.parameters, "model4.txt", n.chains=4, n.iter=5000)

#see .txt file for description
car5.inits<-function(){list(a=rnorm(J), sigma.y=runif(1), mu.a=rnorm(1), sigma.a=runif(1), b0=rnorm(3))}
car.5<-jags(car4.data, car5.inits, car4.parameters, "model.txt", n.chains=4, n.iter=5000)

#print out Bayesian models
pdf("bayes_models.pdf")
plot(car.0)
plot(car.1)
plot(car.2)
plot(car.3)
plot(car.4)
plot(car.5)
dev.off()

#estimate R2 of Bayesian models
attach.jags(car.0)
rsquared0.y<-1-mean(apply(e.y,1,var))/var(na.omit(y))
rsquared0.a<-1-mean(apply(e.a,1,var))/mean(apply(a,1,var))
lambda0.y<-1-var(apply (e.y, 2, mean))/mean(apply (e.y, 1, var))
lambda0.a<-1- var(apply (e.a, 2, mean))/mean(apply(e.a, 1,var))
detach.jags()
attach.jags(car.1)
rsquared.y<-1-mean(apply(e.y,1,var))/var(na.omit(y))
rsquared.a<-1-mean(apply(e.a,1,var))/mean(apply(a,1,var))
lambda.y<-1-var(apply (e.y, 2, mean))/mean(apply (e.y, 1, var))
lambda.a<-1- var(apply (e.a, 2, mean))/mean(apply(e.a, 1,var))
detach.jags()
attach.jags(car.2)
rsquared2.y<-1-mean(apply(e.y,1,var))/var(na.omit(y))
rsquared2.a<-1-mean(apply(e.a,1,var))/mean(apply(a,1,var))
lambda2.y<-1-var(apply (e.y, 2, mean))/mean(apply (e.y, 1, var))
lambda2.a<-1- var(apply (e.a, 2, mean))/mean(apply(e.a, 1,var))
detach.jags()
attach.jags(car.3)
rsquared3.y<-1-mean(apply(e.y,1,var))/var(na.omit(y))
rsquared3.a<-1-mean(apply(e.a,1,var))/mean(apply(a,1,var))
lambda3.y<-1-var(apply (e.y, 2, mean))/mean(apply (e.y, 1, var))
lambda3.a<-1- var(apply (e.a, 2, mean))/mean(apply(e.a, 1,var))
detach.jags()
attach.jags(car.4)
rsquared4.y<-1-mean(apply(e.y,1,var))/var(na.omit(y))
rsquared4.a<-1-mean(apply(e.a,1,var))/mean(apply(a,1,var))
lambda4.y<-1-var(apply (e.y, 2, mean))/mean(apply (e.y, 1, var))
lambda4.a<-1- var(apply (e.a, 2, mean))/mean(apply(e.a, 1,var))
detach.jags()
attach.jags(car.5)
rsquared5.y<-1-mean(apply(e.y,1,var))/var(na.omit(y))
rsquared5.a<-1-mean(apply(e.a,1,var))/mean(apply(a,1,var))
lambda5.y<-1-var(apply (e.y, 2, mean))/mean(apply (e.y, 1, var))
lambda5.a<-1- var(apply (e.a, 2, mean))/mean(apply(e.a, 1,var))
detach.jags()

#print out Bayesian modle results
sink("bayes_model_results.txt")
print("Simplest model no predictors")
print(car.0)
print("Simpler model direct/indirect date")
print(car.1)
print("Complex model direct/indirect date and archeo vs. paleo")
print(car.2)
print("Complex direct/indirect date and archeo vs. paleo, heteroskedastic model")
print(car.3)
print("Simpler model direct/indirect date, two waves, heteroskedastic model")
print(car.4)
print("Simpler model direct/indirect date, two waves")
print(car.5)
sink()
sink("bayes_model_r2.txt")
print("First row r2, second row lambda, simplest no predictors")
print(round(c(rsquared0.y, rsquared0.a),2))
print(round(c(lambda0.y, lambda0.a),2))
print("First row r2, second row lambda, simple direct/indirect date")
print(round(c(rsquared.y, rsquared.a),2))
print(round(c(lambda.y, lambda.a),2))
print("First row r2, second row lambda, complex direct/indirect date and archeo vs. paleo")
print(round(c(rsquared2.y, rsquared2.a),2))
print(round(c(lambda2.y, lambda2.a),2))
print("First row r2, second row lambda, complex direct/indirect date and archeo vs. paleo, heteroskedastic")
print(round(c(rsquared3.y, rsquared3.a),2))
print(round(c(lambda3.y, lambda3.a),2))
print("First row r2, second row lambda, simpler model direct/indirect date, heteroskedastic, two waves")
print(round(c(rsquared4.y, rsquared4.a),2))
print(round(c(lambda4.y, lambda4.a),2))
print("First row r2, second row lambda, simpler model direct/indirect date, two waves")
print(round(c(rsquared5.y, rsquared5.a),2))
print(round(c(lambda5.y, lambda5.a),2))
sink()
save.image("dates.Rdata")
