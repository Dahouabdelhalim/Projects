library(ape)
library(phylopath)
data<-read.csv("Main data.csv",row.names=1)
flytree<-read.tree("Phylogeny data.txt")

data$log.sperm.length<-log10(data$sperm.length)
data$log.male.thorax.length<-log10(data$male.thorax.length)
standardize(data$log.sperm.length)->data$log.sperm.length.s
standardize(data$log.male.thorax.length)->data$log.male.thorax.length.s

lambda <- seq(0, 1, length.out = 500)
gls(log.sperm.length.s~log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree, fixed=TRUE), data, method="ML")->pgls.res
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(log.sperm.length.s~log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree, fixed = TRUE), data = data, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually 1; redo with lambda = 1
gls(log.sperm.length.s~log.male.thorax.length.s, correlation = corPagel(value=1, phy = flytree, fixed=TRUE),data=data, method="ML")->pgls.res.1
data$residual.sperm.length<-resid(pgls.res.1)

models <- define_model_set(
     one   = c(female.mating.frequency ~ sperm.metabolic.rate, sperm.metabolic.rate ~ residual.sperm.length),
     two   = c(female.mating.frequency~ sperm.metabolic.rate, female.mating.frequency ~ residual.sperm.length),
     three   = c(residual.sperm.length ~ female.mating.frequency, female.mating.frequency ~sperm.metabolic.rate),
     four   = c(sperm.metabolic.rate ~residual.sperm.length, sperm.metabolic.rate ~ female.mating.frequency),
     five   = c(female.mating.frequency~ residual.sperm.length, sperm.metabolic.rate ~ female.mating.frequency),
     six = c(female.mating.frequency ~ residual.sperm.length, sperm.metabolic.rate ~ residual.sperm.length))
result <- phylo_path(models, data = data, tree = flytree, model = 'lambda')

#Excluding D. hydei
flytree.nh<-drop.tip(flytree, c("Drosophila_hydei"), trim.internal = TRUE, subtree = FALSE,root.edge = 0, rooted = is.rooted(flytree), collapse.singles = TRUE, interactive = FALSE)
data.nh<-data[-1,]

data.nh$log.sperm.length<-log10(data.nh$sperm.length)
data.nh$log.male.thorax.length<-log10(data.nh$male.thorax.length)
standardize(data.nh$log.sperm.length)->data.nh$log.sperm.length.s
standardize(data.nh$log.male.thorax.length)->data.nh$log.male.thorax.length.s

gls(log.sperm.length.s~log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nh, fixed=TRUE), data.nh, method="ML")->pgls.res.nh
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(log.sperm.length.s~log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nh, fixed = TRUE), data = data.nh, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually 1; redo with lambda = 1
gls(log.sperm.length.s~log.male.thorax.length.s, correlation = corPagel(value=1, phy = flytree.nh, fixed=TRUE),data=data.nh, method="ML")->pgls.res.nh.1

data.nh$residual.sperm.length<-resid(pgls.res.nh.1)
result.nh <- phylo_path(models, data = data.nh, tree = flytree.nh, model = 'lambda')

#Table 2
summary(result)
summary(result.nh)

#Fig. 6
average(result)
average(result.nh)

#Excluding D. affinis
flytree.na<-drop.tip(flytree, c("Drosophila_affinis"), trim.internal = TRUE, subtree = FALSE,root.edge = 0, rooted = is.rooted(flytree), collapse.singles = TRUE, interactive = FALSE)
data.na<-data[-11,]

data.na$log.sperm.length<-log10(data.na$sperm.length)
data.na$log.male.thorax.length<-log10(data.na$male.thorax.length)
standardize(data.na$log.sperm.length)->data.na$log.sperm.length.s
standardize(data.na$log.male.thorax.length)->data.na$log.male.thorax.length.s

flytree.nha<-drop.tip(flytree, c("Drosophila_hydei","Drosophila_affinis"), trim.internal = TRUE, subtree = FALSE,root.edge = 0, rooted = is.rooted(flytree), collapse.singles = TRUE, interactive = FALSE)
data.nha<-data[-c(1,11),]

data.nha$log.sperm.length<-log10(data.nha$sperm.length)
data.nha$log.male.thorax.length<-log10(data.nha$male.thorax.length)
standardize(data.nha$log.sperm.length)->data.nha$log.sperm.length.s
standardize(data.nha$log.male.thorax.length)->data.nha$log.male.thorax.length.s

gls(log.sperm.length.s~log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.na, fixed=TRUE), data.na, method="ML")->pgls.res.na
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(log.sperm.length.s~log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.na, fixed = TRUE), data = data.na, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually 1; redo with lambda = 1
gls(log.sperm.length.s~log.male.thorax.length.s, correlation = corPagel(value=1, phy = flytree.na, fixed=TRUE),data=data.na, method="ML")->pgls.res.na.1
data.na$residual.sperm.length<-resid(pgls.res.na.1)
result.na <- phylo_path(models, data = data.na, tree = flytree.na, model = 'lambda')

gls(log.sperm.length.s~log.male.thorax.length.s, correlation = corPagel(value=lambda, phy=flytree.nha, fixed=TRUE), data.nha, method="ML")->pgls.res.nha
#Check that estimated lambda has highest log likelihood
lik <- sapply(lambda, function(lambda) logLik(gls(log.sperm.length.s~log.male.thorax.length.s,correlation = corPagel(value = lambda, phy = flytree.nha, fixed = TRUE), data = data.nha, method="ML")))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ", lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
#Lambda with highest log likelihood is actually 1; redo with lambda = 1
gls(log.sperm.length.s~log.male.thorax.length.s, correlation = corPagel(value=1, phy = flytree.nha, fixed=TRUE),data=data.nha, method="ML")->pgls.res.nha.1
data.nha$residual.sperm.length<-resid(pgls.res.nha.1)
result.nha <- phylo_path(models, data = data.nha, tree = flytree.nha, model = 'lambda')

#Table S4
summary(result.na)
summary(result.nha)

