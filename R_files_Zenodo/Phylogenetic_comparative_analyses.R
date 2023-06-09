library(ape)
library(geiger)
library(phytools)
library(diversitree)
library(nlme)
library(qpcR)
library(reshape)
library(calibrate)

#####################################################################
# Read in tree
hyper_tree <- read.tree("/Users/portik/Dropbox/Hyperoliids_SSD/Data/Hyperoliid-tree.tre")

# check labels
hyper_tree$tip.label

# plot the tree
plot(hyper_tree, cex = 0.5)


#####################################################################
# Read in body size and SSD data
SSD_data <- read.delim("/Users/portik/Dropbox/Hyperoliids_SSD/Data/Body-size-data.txt", header = TRUE, sep = "\\t")

SSD_frame <- data.frame(SSD_data[,2:13])
rownames(SSD_frame) <- SSD_data[,1]
SSD_frame

#####################################################################
# Nice function for matching tips and sorting 
treeresults <- treedata(hyper_tree, SSD_frame, sort=TRUE, warnings=TRUE)
treeresults

# analysis tree is the pruned phy
analysis_tree <- treeresults$phy

# grab matched data (sorted)
matched_data <- treeresults$data
matched_data

#####################################################################
# Create stand-alone variables so we don't always have to write matched_data$VARIABLE below

SUL_m <- as.numeric(matched_data[,1])
names(SUL_m) <- rownames(matched_data)

SUL_f <- as.numeric(matched_data[,2])
names(SUL_f) <- rownames(matched_data)

log_SUL_m <- as.numeric(matched_data[,3])
names(log_SUL_m) <- rownames(matched_data)

log_SUL_f <- as.numeric(matched_data[,4])
names(log_SUL_f) <- rownames(matched_data)

SSD <- as.numeric(matched_data[,5])
names(SSD) <- rownames(matched_data)

Subfamily <- as.factor(matched_data[,7])
names(Subfamily) <- rownames(matched_data)

Hyperoliinae <- as.factor(matched_data[,8])
names(Hyperoliinae) <- rownames(matched_data)

Hyperolius <- as.factor(matched_data[,9])
names(Hyperolius) <- rownames(matched_data)

Dichromatism <- as.factor(matched_data[,10])
names(Dichromatism) <- rownames(matched_data)

Ecology <- as.factor(matched_data[,11])
names(Ecology) <- rownames(matched_data)

Sexchange <- as.factor(matched_data[,12])
names(Sexchange) <- rownames(matched_data)

#####################################################################
# Summary of data
summary(SUL_m)
summary(SUL_f)
summary(log_SUL_m)
summary(log_SUL_f)
summary(SSD)
summary(Subfamily)
summary(Hyperoliinae)
summary(Hyperolius)
summary(Dichromatism)
summary(Ecology)
summary(Sexchange)

#####################################################################
# Quick data visualization on ancestral states using fastanc function

# male body size
contMap(analysis_tree, SUL_m, res=300, fsize=c(0.5,0.5), lwd = 2, lims = c(15,65), outline = FALSE)
contMap(analysis_tree, log_SUL_m, res=300, fsize=c(0.5,0.5), lwd = 2, outline = FALSE)

# female body size
contMap(analysis_tree, SUL_f, res=300, fsize=c(0.5,0.5), lwd = 2, lims = c(15,65))
contMap(analysis_tree, log_SUL_f, res=300, fsize=c(0.5,0.5), lwd = 2, outline = FALSE)

# SSD
contMap(analysis_tree, SSD, res=300, fsize=c(0.5,0.5), lwd = 2, lims=c(-0.38,0.38))

# change color scheme
ssdcontmap <- contMap(analysis_tree, SSD, res=300, fsize=c(0.5,0.5), lwd = 2)
n<-length(ssdcontmap$cols)
ssdcontmap$cols[1:n]<-colorRampPalette(c("purple","darkgrey","green"), space="Lab")(n)

# plot with new color ramp
plot(ssdcontmap, res=300, fsize=c(0.5,0.5), lwd = 3, outline = FALSE)


#####################################################################
# Exploratory signal measures using lambda and K
# Blomberg's K (Blomberg, S. P., and T. Garland, Jr. 2002. Journal of Evolutionary Biology 15:899-910.), 
# it is the ratio of the observed statistic to that expected under Brownian motion thus K= 1 suggests 
# perfect Brownian motion, K < 1 weaker signal than expected under BM or stronger if K > 1 (i.e. more 
# closely related species are even more similar in their traits than you expect from Brownian motion)
# If K>1 then variance tends to be among clades; while if K<1 then variance is within clades 
# (with BM as reference). The variance on K for a given process is quite large.

# With phylosig() you can specify method (presently either "K" or "lambda") and optionally return a 
# P-value for either the randomization test of Blomberg et al. or a likelihood ratio test against the 
# null hypothesis that lambda=0.

# However, whenever additional factors, unrelated to the phylogenetic history, have an impact on trait 
# evolution, the inﬂuence of the phylogeny needs to be down-weighted. The coeﬃcient k deﬁnes this weight
# and is ﬁtted to observed data such that it scales the Brownian phylogenetic covariances down to the 
# actually observed ones. In other words, k is the transformation of the phylogeny that ensures the best
# ﬁt of trait data to a BM model. Pagel’s k can adopt values larger than one (traits of related species 
# are more similar than expected under BM) but in practice the upper limit is restricted because the oﬀ-
# diagonal elements in the variance–covariance matrix cannot be larger than the diagonal elements 
# (Freckleton, Harvey & Pagel 2002).

# male body size
log_SUL_m.lambda <- phylosig(analysis_tree, log_SUL_m, method="lambda", test=T)
log_SUL_m.K <- phylosig(analysis_tree, log_SUL_m, method="K", nsim=10000, test=T)
log_SUL_m.lambda
log_SUL_m.K

# female body size
log_SUL_f.lambda <- phylosig(analysis_tree, log_SUL_f, method="lambda", test=T)
log_SUL_f.K <- phylosig(analysis_tree, log_SUL_f, method="K", nsim=10000, test=T)
log_SUL_f.lambda
log_SUL_f.K

# SSD
SSD.lambda <- phylosig(analysis_tree, SSD, method="lambda", test=T)
SSD.K <- phylosig(analysis_tree, SSD, method="K", nsim=1000, test=T)
SSD.lambda
SSD.K


#####################################################################
# PGLS

# formula for lm: A typical model has the form Response ~ Terms where response is 
# the (numeric) response vector and terms is a series of terms which specifies a 
# linear predictor for response. A terms specification of the form first + second 
# indicates all the terms in first together with all the terms in second with 
# duplicates removed. (ie lm(y ~ x))


########## Brownian Motion Model ##############

### Examine predictors of male body size

# Brownian motion, fit model where Ecology is a predictor of male SUL 
result.pureBM<-gls(SUL_m ~ Ecology, data=SSD_frame, correlation=corBrownian(phy=analysis_tree),method="ML")
summary(result.pureBM)
plot(SUL_m ~ Ecology, xlab="Ecology", ylab="Male SUL", pch=21, bg="grey", cex = 3)

# Brownian motion, fit model where dichromatism is a predictor of male SUL 
result.pureBM<-gls(SUL_m ~ Dichromatism, data=SSD_frame, correlation=corBrownian(phy=analysis_tree),method="ML")
summary(result.pureBM)
plot(SUL_m ~ Dichromatism, xlab="Color", ylab="Male SUL", pch=21, bg="grey", cex = 3)

# Brownian motion, fit model where sexchange is a predictor of male SUL 
result.pureBM<-gls(SUL_m ~ Sexchange, data=SSD_frame, correlation=corBrownian(phy=analysis_tree),method="ML")
summary(result.pureBM)
plot(SUL_m ~ Sexchange, xlab="Sex Change", ylab="Male SUL", pch=21, bg="grey", cex = 3)


### Examine predictors of female body size

# Brownian motion, fit model where Ecology is a predictor of female SUL 
result.pureBM<-gls(SUL_f ~ Ecology, data=SSD_frame, correlation=corBrownian(phy=analysis_tree),method="ML")
summary(result.pureBM)
plot(SUL_f ~ Ecology, xlab="Ecology", ylab="female SUL", pch=21, bg="grey", cex = 3)

# Brownian motion, fit model where dichromatism is a predictor of female SUL 
result.pureBM<-gls(SUL_f ~ Dichromatism, data=SSD_frame, correlation=corBrownian(phy=analysis_tree),method="ML")
summary(result.pureBM)
plot(SUL_f ~ Dichromatism, xlab="Coloration", ylab="female SUL", pch=21, bg="grey", cex = 3)

# Brownian motion, fit model where sexchange is a predictor of female SUL 
result.pureBM<-gls(SUL_f ~ Sexchange, data=SSD_frame, correlation=corBrownian(phy=analysis_tree),method="ML")
summary(result.pureBM)
plot(SUL_f ~ Sexchange, xlab="Sex Change", ylab="female SUL", pch=21, bg="grey", cex = 3)


### Examine predictors of SSDi

# Brownian motion, fit model where Ecology is a predictor of SSDi 
result.pureBM<-gls(SSD_1 ~ Ecology, data=SSD_frame, correlation=corBrownian(phy=analysis_tree),method="ML")
summary(result.pureBM)
plot(SSD ~ Ecology, xlab="Ecology", ylab="SSDi", pch=21, bg="grey", cex = 3, ylim=c(-0.2,0.4))
 
# Brownian motion, fit model where Dichromatism is a predictor of SSDi 
result.pureBM<-gls(SSD_1 ~ Dichromatism, data=SSD_frame, correlation=corBrownian(phy=analysis_tree),method="ML")
summary(result.pureBM)
plot(SSD ~ Dichromatism, xlab="Coloration", ylab="SSDi", pch=21, bg="grey", cex = 3, ylim=c(-0.2,0.4))

# Brownian motion, fit model where sexchange is a predictor of SSDi 
result.pureBM<-gls(SSD_1 ~ Sexchange, data=SSD_frame, correlation=corBrownian(phy=analysis_tree),method="ML")
summary(result.pureBM)
plot(SSD ~ Sexchange, xlab="Sex Change", ylab="SSDi", pch=21, bg="grey", cex = 3, ylim=c(-0.2,0.4))



########## Brownian Motion + Lambda Model ##############

### Examine predictors of male body size

# Pagel, fit model where Ecology is a predictor of male SUL 
result.pagel<-gls(SUL_m ~ Ecology, data=SSD_frame, correlation=corPagel(1,phy=analysis_tree),method="ML")
summary(result.pagel)
plot(SUL_m ~ Ecology, xlab="Ecology", ylab="male SUL", pch=21, bg="grey", cex = 3)

# Pagel, fit model where dichromatism is a predictor of male SUL 
result.pagel<-gls(SUL_m ~ Dichromatism, data=SSD_frame, correlation=corPagel(1,phy=analysis_tree),method="ML")
summary(result.pagel)
plot(SUL_m ~ Dichromatism, xlab="Coloration", ylab="male SUL", pch=21, bg="grey", cex = 3)

# Pagel, fit model where Sexchange is a predictor of male SUL 
result.pagel<-gls(SUL_m ~ Sexchange, data=SSD_frame, correlation=corPagel(1,phy=analysis_tree),method="ML")
summary(result.pagel)
plot(SUL_m ~ Sexchange, xlab="Sex Change", ylab="male SUL", pch=21, bg="grey", cex = 3)


### Examine predictors of female body size

# Pagel, fit model where Ecology is a predictor of female SUL 
result.pagel<-gls(SUL_f ~ Ecology, data=SSD_frame, correlation=corPagel(1,phy=analysis_tree),method="ML")
summary(result.pagel)
plot(SUL_f ~ Ecology, xlab="Ecology", ylab="female SUL", pch=21, bg="grey", cex = 3)

# Pagel, fit model where dichromatism is a predictor of female SUL 
result.pagel<-gls(SUL_f ~ Dichromatism, data=SSD_frame, correlation=corPagel(1,phy=analysis_tree),method="ML")
summary(result.pagel)
plot(SUL_f ~ Dichromatism, xlab="Coloration", ylab="female SUL", pch=21, bg="grey", cex = 3)

# Pagel, fit model where Sexchange is a predictor of female SUL 
result.pagel<-gls(SUL_f ~ Sexchange, data=SSD_frame, correlation=corPagel(1,phy=analysis_tree),method="ML")
summary(result.pagel)
plot(SUL_f ~ Sexchange, xlab="Sex Change", ylab="female SUL", pch=21, bg="grey", cex = 3)


### Examine predictors of SSDi

# Pagel, fit model where Ecology is a predictor of SSDi 
result.pagel<-gls(SSD_1 ~ Ecology, data=SSD_frame, correlation=corPagel(1,phy=analysis_tree),method="ML")
summary(result.pagel)
plot(SSD ~ Ecology, xlab="Ecology", ylab="SSDi", pch=21, bg="grey", cex = 3, ylim=c(-0.2,0.4))

# Pagel, fit model where Dichromatism is a predictor of SSDi 
result.pagel<-gls(SSD_1 ~ Dichromatism, data=SSD_frame, correlation=corPagel(1,phy=analysis_tree),method="ML")
summary(result.pagel)
plot(SSD ~ Dichromatism, xlab="Coloration", ylab="SSDi", pch=21, bg="grey", cex = 3, ylim=c(-0.2,0.4))

# Pagel, fit model where Sexchange is a predictor of SSDi 
result.pagel<-gls(SSD_1 ~ Sexchange, data=SSD_frame, correlation=corPagel(1,phy=analysis_tree),method="ML")
summary(result.pagel)
plot(SSD ~ Sexchange, xlab="Sex Change", ylab="SSDi", pch=21, bg="grey", cex = 3, ylim=c(-0.2,0.4))


########## OU Model ##############

### Examine predictors of male body size

# OU, fit model where Ecology is a predictor of male SUL 
result.OU<-gls(SUL_m ~ Ecology, data=SSD_frame, correlation=corMartins(0.001,phy=analysis_tree),method="ML")
summary(result.OU)
plot(SUL_m ~ Ecology, xlab="Ecology", ylab="male SUL", pch=21, bg="grey", cex = 3)

# OU, fit model where dichromatism is a predictor of male SUL 
result.OU<-gls(SUL_m ~ Dichromatism, data=SSD_frame, correlation=corMartins(0.001,phy=analysis_tree),method="ML")
summary(result.OU)
plot(SUL_m ~ Dichromatism, xlab="Coloration", ylab="male SUL", pch=21, bg="grey", cex = 3)

# OU, fit model where Sexchange is a predictor of male SUL 
result.OU<-gls(SUL_m ~ Sexchange, data=SSD_frame, correlation=corMartins(0.001,phy=analysis_tree),method="ML")
summary(result.OU)
plot(SUL_m ~ Sexchange, xlab="Sex Change", ylab="male SUL", pch=21, bg="grey", cex = 3)


### Examine predictors of female body size

# OU, fit model where Ecology is a predictor of female SUL 
result.OU<-gls(SUL_f ~ Ecology, data=SSD_frame, correlation=corMartins(0.001,phy=analysis_tree),method="ML")
summary(result.OU)
plot(SUL_f ~ Ecology, xlab="Ecology", ylab="female SUL", pch=21, bg="grey", cex = 3)

# OU, fit model where dichromatism is a predictor of female SUL 
result.OU<-gls(SUL_f ~ Dichromatism, data=SSD_frame, correlation=corMartins(0.001,phy=analysis_tree),method="ML")
summary(result.OU)
plot(SUL_f ~ Dichromatism, xlab="Coloration", ylab="female SUL", pch=21, bg="grey", cex = 3)

# OU, fit model where Sexchange is a predictor of female SUL 
result.OU<-gls(SUL_f ~ Sexchange, data=SSD_frame, correlation=corMartins(0.001,phy=analysis_tree),method="ML")
summary(result.OU)
plot(SUL_f ~ Sexchange, xlab="Sex Change", ylab="female SUL", pch=21, bg="grey", cex = 3)


### Examine predictors of SSDi

# OU, fit model where Ecology is a predictor of SSDi 
result.OU<-gls(SSD_1 ~ Ecology, data=SSD_frame, correlation=corMartins(0.01,phy=analysis_tree),method="ML")
summary(result.OU)
plot(SSD ~ Ecology, xlab="Ecology", ylab="SSDi", pch=21, bg="grey", cex = 3, ylim=c(-0.2,0.4))

# OU, fit model where Dichromatism is a predictor of SSDi 
result.OU<-gls(SSD_1 ~ Dichromatism, data=SSD_frame, correlation=corMartins(0.01,phy=analysis_tree),method="ML")
summary(result.OU)
plot(SSD ~ Dichromatism, xlab="Coloration", ylab="SSDi", pch=21, bg="grey", cex = 3, ylim=c(-0.2,0.4))

# OU, fit model where Sexchange is a predictor of SSDi 
result.OU<-gls(SSD_1 ~ Sexchange, data=SSD_frame, correlation=corMartins(0.01,phy=analysis_tree),method="ML")
summary(result.OU)
plot(SSD ~ Sexchange, xlab="Sex Change", ylab="SSDi", pch=21, bg="grey", cex = 3, ylim=c(-0.2,0.4))


#####################################################################
# Phylogenetic reduced major axis regressions
# phyl.RMA(x, y, tree, method="BM", lambda=NULL, fixed=FALSE, h0=1.0)

# First perform on all hyperoliid body sizes with females on x axis and males on y
femresult.RMA <- phyl.RMA(log_SUL_f, log_SUL_m, analysis_tree, method="BM", lambda=NULL, fixed=FALSE)
femresult.RMA$test
femresult.RMA$RMA.beta

plot(y=log_SUL_m, x=log_SUL_f, ylab="log SUL males", xlab="log SUL females", pch=21, bg="grey", cex=2.5, xlim=c(2.8, 4.2), ylim=c(2.8, 4.2))
abline(femresult.RMA$RMA.beta, lwd=3)
# plot the isometry line
abline(0,1, lwd=3, lty=3)


# Next perform on all hyperoliid body sizes with males on x axis and females on y
maleresult.RMA <- phyl.RMA(log_SUL_m, log_SUL_f, analysis_tree, method="BM", lambda=NULL, fixed=FALSE)
maleresult.RMA$test
maleresult.RMA$RMA.beta

plot(x=log_SUL_m, y=log_SUL_f, xlab="log SUL males", ylab="log SUL females", pch=21, bg="grey", cex=2.5, xlim=c(2.8, 4.2), ylim=c(2.8, 4.2))
abline(maleresult.RMA$RMA.beta, lwd=3)
#a nd the isometry line
abline(0,1, lwd=3, lty=3)

# notice they are the same slope!



#####################################################################
# Perform PGLS on subset of hyperoliid data for sex change species
#####################################################################

# Extract data for hyperoliioids
# this is required to test whether sex change is a predictor
# of body sizes or SSDi within this sub-group

# read in tree
hyper_tree <- read.tree("/Users/portik/Dropbox/Hyperoliids_SSD/Data/Hyperoliid-tree.tre")

plot(hyper_tree, cex = 0.5)

# get node number for mrca of group
mrca(hyper_tree)["Cryptothylax_greshoffii", "Hyperolius_argus"]
hyp_tree <- extract.clade(hyper_tree, 142, root.edge = 0)

# check that this pruned correctly
plot(hyp_tree, cex = 0.5)

# read in data subset
SSD_data1 <- read.delim("/Users/portik/Dropbox/Hyperoliids_SSD/Data/Hyperolius_Sexchange.txt", header = TRUE, sep = "\\t")

# create new data frame
SSD_frame1 <- data.frame(SSD_data1[,2:13])
rownames(SSD_frame1) <- SSD_data1[,1]
SSD_frame1

# check labels match
hyp_results <- treedata(hyp_tree, SSD_frame1, sort=TRUE, warnings=TRUE)

# assign tree and data
hyp_analysis <- hyp_results$phy
hyp_morpho <- hyp_results$data

# assign names to variables
hyp_SUL_m <- as.numeric(hyp_morpho[,1])
names(hyp_SUL_m) <- rownames(hyp_morpho)
hyp_SUL_f <- as.numeric(hyp_morpho[,2])
names(hyp_SUL_f) <- rownames(hyp_morpho)
hyp_log_SUL_m <- as.numeric(hyp_morpho[,3])
names(hyp_log_SUL_m) <- rownames(hyp_morpho)
hyp_log_SUL_f <- as.numeric(hyp_morpho[,4])
names(hyp_log_SUL_f) <- rownames(hyp_morpho)
hyp_SSD <- as.numeric(hyp_morpho[,5])
names(hyp_SSD) <- rownames(hyp_morpho)
hyp_Hyperolius <- as.factor(hyp_morpho[,9])
names(hyp_Hyperolius) <- rownames(hyp_morpho)
hyp_Dichromatism <- as.factor(hyp_morpho[,10])
names(hyp_Dichromatism) <- rownames(hyp_morpho)
hyp_Sexchange <- as.factor(hyp_morpho[,12])
names(hyp_Sexchange) <- rownames(hyp_morpho)

# get summary of each variable
summary(hyp_SUL_m)
summary(hyp_SUL_f)
summary(hyp_log_SUL_m)
summary(hyp_log_SUL_f)
summary(hyp_SSD)
summary(hyp_Hyperolius)
summary(hyp_Dichromatism)
summary(hyp_Sexchange)


##########  Perform PGLS on data subset ########## 

########## Brownian Motion Model

#Brownian motion, fit model where sexchange is a predictor of SSDi 
result.pureBM<-gls(SSD_1 ~ Sex, data=SSD_frame1, correlation=corBrownian(phy=hyp_analysis),method="ML")
summary(result.pureBM)

#Brownian motion, fit model where sexchange is a predictor of female SUL 
result.pureBM<-gls(SUL_m ~ Sex, data=SSD_frame1, correlation=corBrownian(phy=hyp_analysis),method="ML")
summary(result.pureBM)

#Brownian motion, fit model where sexchange is a predictor of male SUL 
result.pureBM<-gls(SUL_f ~ Sex, data=SSD_frame1, correlation=corBrownian(phy=hyp_analysis),method="ML")
summary(result.pureBM)


########## OU Model

#OU, fit model where Sexchange is a predictor of SSDi 
result.OU<-gls(SSD_1 ~ Sex, data=SSD_frame1, correlation=corMartins(0.01,phy=hyp_analysis),method="ML")
summary(result.OU)

#OU, fit model where Sexchange is a predictor of male SUL 
result.OU<-gls(SUL_m ~ Sex, data=SSD_frame1, correlation=corMartins(0.01,phy=hyp_analysis),method="ML")
summary(result.OU)

#OU, fit model where Sexchange is a predictor of female SUL 
result.OU<-gls(SUL_f ~ Sex, data=SSD_frame1, correlation=corMartins(0.001,phy=hyp_analysis),method="ML")
summary(result.OU)


########## Brownian Motion + Lambda Model

#Pagel, fit model where Sexchange is a predictor of SSDi 
result.pagel<-gls(SSD_1 ~ Sex, data=SSD_frame1, correlation=corPagel(1,phy=hyp_analysis),method="ML")
summary(result.pagel)

#Pagel, fit model where Sexchange is a predictor of male SUL 
result.pagel<-gls(SUL_m ~ Sex, data=SSD_frame1, correlation=corPagel(1,phy=hyp_analysis),method="ML")
summary(result.pagel)

#Pagel, fit model where Sexchange is a predictor of female SUL 
result.pagel<-gls(SUL_f ~ Sex, data=SSD_frame1, correlation=corPagel(1,phy=hyp_analysis),method="ML")
summary(result.pagel)



#####################################################################
# Perform phyRMA on dichromatic and monochromatic species
#####################################################################

#####################################################################
# Read in data for dichromatic species
dichro_data <- read.delim("/Users/portik/Dropbox/Hyperoliids_SSD/Data/Hyperolius_Dichro.txt", header = TRUE, sep = "\\t")

dichro_frame <- data.frame(dichro_data[,2:11])
rownames(dichro_frame) <- dichro_data[,1]
dichro_frame

# Nice function for matching tips and sorting 
dichrotreeresults <- treedata(hyper_tree, dichro_frame, sort=TRUE, warnings=TRUE)
dichrotreeresults

# analysis tree is the pruned phy
dichro_tree <- dichrotreeresults$phy

# grab matched data (sorted)
dichro_matched_data <- dichrotreeresults$data
dichro_matched_data

# Create stand-alone variables
dichro_log_SUL_m <- as.numeric(dichro_matched_data[,3])
names(dichro_log_SUL_m) <- rownames(dichro_matched_data)

dichro_log_SUL_f <- as.numeric(dichro_matched_data[,4])
names(dichro_log_SUL_f) <- rownames(dichro_matched_data)

dichro_SSD <- as.numeric(dichro_matched_data[,5])
names(dichro_SSD) <- rownames(dichro_matched_data)


#####################################################################
# Read in data for monochromatic species
mono_data <- read.delim("/Users/portik/Dropbox/Hyperoliids_SSD/Data/Hyperolius_Mono.txt", header = TRUE, sep = "\\t")

mono_frame <- data.frame(mono_data[,2:11])
rownames(mono_frame) <- mono_data[,1]
mono_frame

# Nice function for matching tips and sorting 
monotreeresults <- treedata(hyper_tree, mono_frame, sort=TRUE, warnings=TRUE)
monotreeresults

# analysis tree is the pruned phy
mono_tree <- monotreeresults$phy

# grab matched data (sorted)
mono_matched_data <- monotreeresults$data
mono_matched_data

# Create stand-alone variables 
mono_log_SUL_m <- as.numeric(mono_matched_data[,3])
names(mono_log_SUL_m) <- rownames(mono_matched_data)

mono_log_SUL_f <- as.numeric(mono_matched_data[,4])
names(mono_log_SUL_f) <- rownames(mono_matched_data)

mono_SSD <- as.numeric(mono_matched_data[,5])
names(mono_SSD) <- rownames(mono_matched_data)

#####################################################################
# now perform phyRMAs for each of the subsets

# First perform on all dichromatic body sizes with females on x axis and males on y
dichro.RMA <- phyl.RMA(dichro_log_SUL_f, dichro_log_SUL_m, dichro_tree, method="BM", lambda=NULL, fixed=FALSE)
dichro.RMA$test
dichro.RMA$RMA.beta

plot(y=dichro_log_SUL_m, x=dichro_log_SUL_f, ylab="log SUL males", xlab="log SUL females", pch=21, bg="green", cex=2.5, xlim=c(2.8, 4.2), ylim=c(2.8, 4.2))
abline(dichro.RMA$RMA.beta, lwd=3)
#and the isometry line
abline(0,1, lwd=3, lty=3)


# Then perform on all monochromatic body sizes with females on x axis and males on y
mono.RMA <- phyl.RMA(mono_log_SUL_f, mono_log_SUL_m, mono_tree, method="BM", lambda=NULL, fixed=FALSE)
mono.RMA$test
mono.RMA$RMA.beta

plot(y=mono_log_SUL_m, x=mono_log_SUL_f, ylab="log SUL males", xlab="log SUL females", pch=21, bg="darkgrey", cex=2.5, xlim=c(2.8, 4.2), ylim=c(2.8, 4.2))
abline(mono.RMA$RMA.beta, lwd=3)
# and the isometry line
abline(0,1, lwd=3, lty=3)
