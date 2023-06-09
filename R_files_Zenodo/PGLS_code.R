library(ape)
library(nlme)
library(robustHD)
data<-read.csv("Main data.csv",row.names=1)
flytree<-read.tree("Phylogeny data.txt")

data$log.sperm.length<-log10(data$sperm.length)
data$log.seminal.receptacle.length<-log10(data$seminal.receptacle.length)
data$log.male.thorax.length<-log10(data$male.thorax.length)
data$log.female.thorax.length<-log10(data$female.thorax.length)

standardize(data$sperm.metabolic.rate)->data$sperm.metabolic.rate.s
standardize(data$female.mating.frequency)->data$female.mating.frequency.s
standardize(data$sperm.length)->data$sperm.length.s
standardize(data$log.sperm.length)->data$log.sperm.length.s
standardize(data$log.seminal.receptacle.length)->data$log.seminal.receptacle.length.s
standardize(data$log.male.thorax.length)->data$log.male.thorax.length.s
standardize(data$log.female.thorax.length)->data$log.female.thorax.length.s

#Exclude D. hydei
data.nh<-data[-1,]
flytree.nh<-drop.tip(flytree, c("Drosophila_hydei"), trim.internal = TRUE, subtree = FALSE,root.edge = 0, rooted = is.rooted(flytree), collapse.singles = TRUE, interactive = FALSE)

data.nh$log.sperm.length<-log10(data.nh$sperm.length)
data.nh$log.seminal.receptacle.length<-log10(data.nh$seminal.receptacle.length)
data.nh$log.male.thorax.length<-log10(data.nh$male.thorax.length)
data.nh$log.female.thorax.length<-log10(data.nh$female.thorax.length)

standardize(data.nh$sperm.metabolic.rate)->data.nh$sperm.metabolic.rate.s
standardize(data.nh$female.mating.frequency)->data.nh$female.mating.frequency.s
standardize(data.nh$sperm.length)->data.nh$sperm.length.s
standardize(data.nh$log.sperm.length)->data.nh$log.sperm.length.s
standardize(data.nh$log.seminal.receptacle.length)->data.nh$log.seminal.receptacle.length.s
standardize(data.nh$log.male.thorax.length)->data.nh$log.male.thorax.length.s
standardize(data.nh$log.female.thorax.length)->data.nh$log.female.thorax.length.s

#Table 1

#Model 1 (fig 4)
lambda <- seq(0, 1, length.out = 500)
gls(female.mating.frequency.s~sperm.metabolic.rate.s, correlation = corPagel(value=lambda, phy=flytree, fixed=TRUE), data, method="ML")->pgls1
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s,correlation = corPagel(value = lambda, phy = flytree, fixed = TRUE), data = data, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model 1.nh (fig 4)
gls(female.mating.frequency.s~sperm.metabolic.rate.s, correlation = corPagel(value=lambda, phy=flytree.nh, fixed=TRUE), data.nh, method="ML")->pgls1.nh
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s,correlation = corPagel(value = lambda, phy = flytree.nh, fixed = TRUE), data = data.nh, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model 2
gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree, fixed=TRUE), data, method="ML")->pgls2
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.sperm.length.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree, fixed = TRUE), data = data, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model 2.nh
gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nh, fixed=TRUE), data.nh, method="ML")->pgls2.nh
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate+log.sperm.length.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nh, fixed = TRUE), data = data.nh, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model 3
gls(female.mating.frequency.s~log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree, fixed=TRUE), data, method="ML")->pgls3
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~log.sperm.length.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree, fixed = TRUE), data = data, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model 3.nh
gls(female.mating.frequency.s~log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nh, fixed=TRUE), data.nh, method="ML")->pgls3.nh
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~log.sperm.length.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nh, fixed = TRUE), data = data.nh, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model 4 (fig 5)
gls(sperm.metabolic.rate.s~log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree, fixed=TRUE), data, method="ML")->pgls4
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(sperm.metabolic.rate.s~log.sperm.length.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree, fixed = TRUE), data = data, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model 4.nh (Fig 5)
gls(sperm.metabolic.rate.s~log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nh, fixed=TRUE), data.nh, method="ML")->pgls4.nh
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(sperm.metabolic.rate.s~log.sperm.length.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nh, fixed = TRUE), data = data.nh, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually near 0.5; redo with lambda starting value set to 0.5
gls(sperm.metabolic.rate.s~log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=0.5, phy = flytree.nh, fixed=FALSE),data=data.nh, method="ML")->pgls4.1.nh

#Table S1

#Model S1
gls(female.mating.frequency.s~log.seminal.receptacle.length.s+log.female.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree, fixed=TRUE), data, method="ML")->pglsS1
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~log.seminal.receptacle.length.s+log.female.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree, fixed = TRUE), data = data, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model S1.nh
gls(female.mating.frequency.s~log.seminal.receptacle.length.s+log.female.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nh, fixed=TRUE), data.nh, method="ML")->pglsS1.nh
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~log.seminal.receptacle.length.s+log.female.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nh, fixed = TRUE), data = data.nh, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model S2
gls(log.sperm.length.s~log.seminal.receptacle.length.s+log.male.thorax.length.s+log.female.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree, fixed=TRUE), data, method="ML")->pglsS2
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(log.sperm.length.s~log.seminal.receptacle.length.s+log.male.thorax.length.s+log.female.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree, fixed = TRUE), data = data, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually 1; redo with lambda = 1
gls(log.sperm.length.s~log.seminal.receptacle.length.s+log.male.thorax.length.s+log.female.thorax.length.s, correlation = corPagel(value=1, phy = flytree, fixed=TRUE),data=data, method="ML")->pglsS2.1

#Model S2.nh
gls(log.sperm.length.s~log.seminal.receptacle.length.s+log.male.thorax.length.s+log.female.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nh, fixed=TRUE), data.nh, method="ML")->pglsS2.nh
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(log.sperm.length.s~log.seminal.receptacle.length.s+log.male.thorax.length.s+log.female.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nh, fixed = TRUE), data = data.nh, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually 1; redo with lambda = 1
gls(log.sperm.length.s~log.seminal.receptacle.length.s+log.male.thorax.length.s+log.female.thorax.length.s, correlation = corPagel(value=1, phy = flytree.nh, fixed=TRUE),data=data.nh, method="ML")->pglsS2.nh.1

#Table S2

#Model S3
gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree, fixed=TRUE), data, method="ML")->pglsS3
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree, fixed = TRUE), data = data, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model S3.nh
gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nh, fixed=TRUE), data.nh, method="ML")->pglsS3.nh
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nh, fixed = TRUE), data = data.nh, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Calculate adjusted variance
raw.data<-read.csv("Raw data.csv")
a2_mean = tapply(raw.data$sperm.metabolic.rate, raw.data$species, mean)
a2_N = tapply(raw.data$sperm.metabolic.rate, raw.data$species, length)
a2_var = tapply(raw.data$sperm.metabolic.rate, raw.data$species, var)
a2_cv = sqrt(a2_var)/a2_mean
a2_cv.pool = sum(a2_cv * a2_N, na.rm = T)/sum(a2_N, na.rm = T)
a2_var.adj = (a2_cv.pool * a2_mean)^2
#This is estimated.sperm.metabolic.rate.variance in "Main data.csv"

#Model S4 (Model 1 accounting for within-species variance)
gls(female.mating.frequency.s~sperm.metabolic.rate.s, correlation = corPagel(value=lambda, phy=flytree, fixed=TRUE), data, method="ML",weights = varFixed(~estimated.sperm.metabolic.rate.variance))->pglsS4
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s,correlation = corPagel(value = lambda, phy = flytree, fixed = TRUE),weights = varFixed(~estimated.sperm.metabolic.rate.variance), data = data, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model S4.nh (Model 1.nh accounting for within-species variance)
gls(female.mating.frequency.s~sperm.metabolic.rate.s, correlation = corPagel(value=lambda, phy=flytree.nh, fixed=TRUE), data.nh, method="ML",weights = varFixed(~estimated.sperm.metabolic.rate.variance))->pglsS4.nh
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s,correlation = corPagel(value = lambda, phy = flytree.nh, fixed = TRUE),weights = varFixed(~estimated.sperm.metabolic.rate.variance), data = data.nh, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model S5
gls(sperm.metabolic.rate.s~sperm.length.s, correlation = corPagel(value=lambda, phy=flytree, fixed=TRUE), data, method="ML")->pglsS5
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(sperm.metabolic.rate.s~sperm.length.s,correlation = corPagel(value = lambda, phy = flytree, fixed = TRUE), data = data, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually near 1; redo with lambda = 1
gls(sperm.metabolic.rate.s~sperm.length.s, correlation = corPagel(value=1, phy = flytree, fixed=FALSE),data=data, method="ML")->pglsS5.1

#Model S5.nh
gls(sperm.metabolic.rate.s~sperm.length.s, correlation = corPagel(value=lambda, phy=flytree.nh, fixed=TRUE), data.nh, method="ML")->pglsS5.nh
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(sperm.metabolic.rate~sperm.length.s,correlation = corPagel(value = lambda, phy = flytree.nh, fixed = TRUE), data = data.nh, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually near 1; redo with lambda = 1
gls(sperm.metabolic.rate.s~sperm.length.s, correlation = corPagel(value=1, phy = flytree.nh, fixed=FALSE),data=data.nh, method="ML")->pglsS5.nh.1

#Redo excluding D. affinis
data.na<-data[-11,]
flytree.na<-drop.tip(flytree, c("Drosophila_affinis"), trim.internal = TRUE, subtree = FALSE,root.edge = 0, rooted = is.rooted(flytree), collapse.singles = TRUE, interactive = FALSE)

data.na$log.sperm.length<-log10(data.na$sperm.length)
data.na$log.seminal.receptacle.length<-log10(data.na$seminal.receptacle.length)
data.na$log.male.thorax.length<-log10(data.na$male.thorax.length)
data.na$log.female.thorax.length<-log10(data.na$female.thorax.length)

standardize(data.na$sperm.metabolic.rate)->data.na$sperm.metabolic.rate.s
standardize(data.na$female.mating.frequency)->data.na$female.mating.frequency.s
standardize(data.na$sperm.length)->data.na$sperm.length.s
standardize(data.na$log.sperm.length)->data.na$log.sperm.length.s
standardize(data.na$log.seminal.receptacle.length)->data.na$log.seminal.receptacle.length.s
standardize(data.na$log.male.thorax.length)->data.na$log.male.thorax.length.s
standardize(data.na$log.female.thorax.length)->data.na$log.female.thorax.length.s

data.nha<-data[-c(1,11),]
flytree.nha<-drop.tip(flytree, c("Drosophila_hydei","Drosophila_affinis"), trim.internal = TRUE, subtree = FALSE,root.edge = 0, rooted = is.rooted(flytree), collapse.singles = TRUE, interactive = FALSE)

data.nha$log.sperm.length<-log10(data.nha$sperm.length)
data.nha$log.seminal.receptacle.length<-log10(data.nha$seminal.receptacle.length)
data.nha$log.male.thorax.length<-log10(data.nha$male.thorax.length)
data.nha$log.female.thorax.length<-log10(data.nha$female.thorax.length)

standardize(data.nha$sperm.metabolic.rate)->data.nha$sperm.metabolic.rate.s
standardize(data.nha$female.mating.frequency)->data.nha$female.mating.frequency.s
standardize(data.nha$sperm.length)->data.nha$sperm.length.s
standardize(data.nha$log.sperm.length)->data.nha$log.sperm.length.s
standardize(data.nha$log.seminal.receptacle.length)->data.nha$log.seminal.receptacle.length.s
standardize(data.nha$log.male.thorax.length)->data.nha$log.male.thorax.length.s
standardize(data.nha$log.female.thorax.length)->data.nha$log.female.thorax.length.s

#Table S3

#Model 1.na (fig 4)
gls(female.mating.frequency.s~sperm.metabolic.rate.s, correlation = corPagel(value=lambda, phy=flytree.na, fixed=TRUE), data.na, method="ML")->pgls1.na
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s,correlation = corPagel(value = lambda, phy = flytree.na, fixed = TRUE), data = data.na, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model 1.nha (fig 4)
gls(female.mating.frequency.s~sperm.metabolic.rate.s, correlation = corPagel(value=lambda, phy=flytree.nha, fixed=TRUE), data.nha, method="ML")->pgls1.nha
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s,correlation = corPagel(value = lambda, phy = flytree.nha, fixed = TRUE), data = data.nha, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model 2.na
gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.na, fixed=TRUE), data.na, method="ML")->pgls2.na
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.sperm.length.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.na, fixed = TRUE), data = data.na, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model 2.nha
gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nha, fixed=TRUE), data.nha, method="ML")->pgls2.nha
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.sperm.length.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nha, fixed = TRUE), data = data.nha, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model 3.na
gls(female.mating.frequency.s~log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.na, fixed=TRUE), data.na, method="ML")->pgls3.na
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~log.sperm.length.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.na, fixed = TRUE), data = data.na, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model 3.nha
gls(female.mating.frequency.s~log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nha, fixed=TRUE), data.nha, method="ML")->pgls3.nha
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~log.sperm.length.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nha, fixed = TRUE), data = data.nha, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model 4.na (fig 5)
gls(sperm.metabolic.rate.s~log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.na, fixed=TRUE), data.na, method="ML")->pgls4.na
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(sperm.metabolic.rate.s~log.sperm.length.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.na, fixed = TRUE), data = data.na, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually near 0.5; redo with lambda = 0.5
gls(sperm.metabolic.rate.s~log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=0.5, phy = flytree.na, fixed=FALSE),data=data.na, method="ML")->pgls4.na.1

#Model 4.nha (fig 5)
gls(sperm.metabolic.rate.s~log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nha, fixed=TRUE), data.nha, method="ML")->pgls4.nha
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(sperm.metabolic.rate.s~log.sperm.length.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nha, fixed = TRUE), data = data.nha, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually near 0.5; redo with lambda = 0.5
gls(sperm.metabolic.rate.s~log.sperm.length.s+log.male.thorax.length.s, correlation = corPagel(value=0.5, phy = flytree.nha, fixed=FALSE),data=data.nha, method="ML")->pgls4.nha.1

#Table S5

#Model S1.na
gls(female.mating.frequency.s~log.seminal.receptacle.length.s+log.female.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.na, fixed=TRUE), data.na, method="ML")->pglsS1.na
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~log.seminal.receptacle.length.s+log.female.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.na, fixed = TRUE), data = data.na, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model S1.nha
gls(female.mating.frequency.s~log.seminal.receptacle.length.s+log.female.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nha, fixed=TRUE), data.nha, method="ML")->pglsS1.nha
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~log.seminal.receptacle.length.s+log.female.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nha, fixed = TRUE), data = data.nha, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model S2.na
gls(log.sperm.length.s~log.seminal.receptacle.length.s+log.male.thorax.length.s+log.female.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.na, fixed=TRUE), data.na, method="ML")->pglsS2.na
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(log.sperm.length.s~log.seminal.receptacle.length.s+log.male.thorax.length.s+log.female.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.na, fixed = TRUE), data = data.na, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually 1; redo with lambda = 1
gls(log.sperm.length.s~log.seminal.receptacle.length.s+log.male.thorax.length.s+log.female.thorax.length.s, correlation = corPagel(value=1, phy = flytree.na, fixed=TRUE),data=data.na, method="ML")->pglsS2.na.1

#Model S2.nha
gls(log.sperm.length.s~log.seminal.receptacle.length.s+log.male.thorax.length.s+log.female.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nha, fixed=TRUE), data.nha, method="ML")->pglsS2.nha
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(log.sperm.length.s~log.seminal.receptacle.length.s+log.male.thorax.length.s+log.female.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nha, fixed = TRUE), data = data.nha, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually 1; redo with lambda = 1
gls(log.sperm.length.s~log.seminal.receptacle.length.s+log.male.thorax.length.s+log.female.thorax.length.s, correlation = corPagel(value=1, phy = flytree.nha, fixed=TRUE),data=data.nha, method="ML")->pglsS2.nha.1

#Table S6

#Model S3.na
gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.na, fixed=TRUE), data.na, method="ML")->pglsS3.na
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.na, fixed = TRUE), data = data.na, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model S3.nha
gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nha, fixed=TRUE), data.nha, method="ML")->pglsS3.nha
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s+log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nha, fixed = TRUE), data = data.nha, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model S4.na (Model 1.na accounting for within-species variance)
gls(female.mating.frequency.s~sperm.metabolic.rate.s, correlation = corPagel(value=lambda, phy=flytree.na, fixed=TRUE), data.na, method="ML",weights = varFixed(~estimated.sperm.metabolic.rate.variance))->pglsS4.na
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s,correlation = corPagel(value = lambda, phy = flytree.na, fixed = TRUE),weights = varFixed(~estimated.sperm.metabolic.rate.variance), data = data.na, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model S4.nha (Model 1.nha accounting for within-species variance)
gls(female.mating.frequency.s~sperm.metabolic.rate.s, correlation = corPagel(value=lambda, phy=flytree.nha, fixed=TRUE), data.nha, method="ML",weights = varFixed(~estimated.sperm.metabolic.rate.variance))->pglsS4.nha
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(female.mating.frequency.s~sperm.metabolic.rate.s,correlation = corPagel(value = lambda, phy = flytree.nha, fixed = TRUE),weights = varFixed(~estimated.sperm.metabolic.rate.variance), data = data.nha, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model S5.na
gls(sperm.metabolic.rate.s~sperm.length.s, correlation = corPagel(value=lambda, phy=flytree.nh, fixed=TRUE), data.na, method="ML")->pglsS5.na
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(sperm.metabolic.rate.s~sperm.length.s,correlation = corPagel(value = lambda, phy = flytree.na, fixed = TRUE), data = data.na, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))

#Model S5.nha
gls(sperm.metabolic.rate.s~sperm.length.s, correlation = corPagel(value=lambda, phy=flytree.nh, fixed=TRUE), data.nha, method="ML")->pglsS5.nha
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(sperm.metabolic.rate.s~sperm.length.s,correlation = corPagel(value = lambda, phy = flytree.nha, fixed = TRUE), data = data.nha, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually near 1; redo with lambda = 1
gls(sperm.metabolic.rate.s~sperm.length.s, correlation = corPagel(value=1, phy = flytree.nha, fixed=FALSE),data=data.nha, method="ML")->pglsS5.nha.1

#asymptotic test for the equality of coefficients of variation
library(cvequality)
with(raw.data, asymptotic_test(sperm.metabolic.rate, species))
