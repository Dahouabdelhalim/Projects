library(xlsx)
library(lmerTest)
library(lme4)
library(erer)
library(nlme)
library(lme4)
library(PolyPatEx)
library(qdapTools)
library(car)
###Define functions

id.tr <- function(x){
	x <- x
}

mean.scale<-function (vec){
	(1/mean(vec))* vec
}


c.factor <- function(...){
  as.factor(do.call(c, lapply(list(...), as.character)))
}

zw.prep1 <- function(w,transform){
	#w <- w[!is.na(w$ft),]
	w$rf <- mean.scale(transform(w$ft))
	for (i in 5: 16){
		w[,i] <- scale(w[,i])
	}
	return(w)
}

#Prepare data for Ancova
Anc.prep1 <- function(w,x,y,z,tr){
	zero1 <-(c.factor(w[is.na(w$ft),"ID"],x[is.na(x$ft),"ID"]))
	zero2 <-(c.factor(y[is.na(y$ft),"ID"],z[is.na(z$ft),"ID"]))
	w <- w[!is.element(w$ID,zero1),]
	x <- x[!is.element(x$ID,zero1),]
	y <- y[!is.element(y$ID,zero2),]
	z <- z[!is.element(z$ID,zero2),]
	return(rbind(zw.prep1(w,tr),zw.prep1(x,tr),zw.prep1(y,tr),zw.prep1(z,tr)))
	}


#Fit five ancova models
int <- function(anc.dat){
intn <- lmer(rf ~  (1|ID), REML=FALSE, data= anc.dat)
int0 <- lmer(rf ~  Sex + E + Year + (1|ID), REML=FALSE, data=anc.dat)
int1 <- lmer(rf ~ (TL + TW.F + PW + PL + LB + STE + ANE + NF + HT) + Sex + E + Year + (1|ID), REML=FALSE, data=anc.dat)
int2 <- lmer(rf ~ (TL + TW.F + PW + PL + LB + STE + ANE + NF+ HT)  + Sex + E + Year + (TL + TW.F + PW + PL + LB + STE + ANE + NF+ HT)*E + (1|ID),REML=FALSE, data=anc.dat)
int3 <- lmer(rf ~ (TL + TW.F + PW + PL + LB + STE + ANE + NF+ HT)  + Sex + E + Year + (TL + TW.F + PW + PL + LB + STE + ANE + NF + HT)*Sex + (1|ID) ,REML=FALSE, data=anc.dat)
int4 <- lmer(rf ~ (TL + TW.F + PW + PL + LB + STE + ANE + NF+ HT) +(TL + TW.F + PW + PL + LB + STE + ANE + NF+ HT)*E + (TL + TW.F + PW + PL + LB + STE + ANE + NF + HT)*Sex + (1|ID) + Year + Sex + E, REML=FALSE,data=anc.dat)
int5 <- lmer(rf ~ (TL + TW.F + PW + PL + LB + STE + ANE + NF+ HT) + (TL + TW.F + PW + PL + LB + STE + ANE + NF+ HT)*E + (TL + TW.F + PW + PL + LB + STE + ANE + NF + HT)*Sex + (TL + TW.F + PW + PL + LB + STE + ANE + NF+ HT)*Sex*E +  + Year + Sex + E + (1|ID),REML=FALSE, data=anc.dat)
rt <- list(summary(intn), summary(int0), summary(int1),summary(int2), summary(int3), summary(int4), summary(int5))
out <- as.data.frame(matrix(0,nrow=7,ncol=5),row.names=c("M_Null","M0","M1","M2","M3","M4","M5"))
names(out) <- c("AIC", "BIC","logLik","deviance","df.resid")
for (i in 1:7) {
	out[i,] <- rt[[i]]$AICtab
}
return(out)}

#Fit linear regression
myp <- function(x)
{
lmer(ft/mean(ft) ~ scale(TL) + scale(TW.F)  + scale(PL)  + scale(PW) + scale(LB) + scale(STE) + scale(ANE) +  scale(NF) + scale(HT) + (1| Year) ,data = x)
}

#Permutation p value
perm <- function(x)
{
	x[sample(1:length(x),length(x) )]
}

permdata <- function(dat){
d <- dat
colnames(d)<-colnames(dat)
d$ft <- perm(dat$ft)
d
}
permP <- function(data,g, nrep) {
	fobs=Anova(g(data),type=3)[,1]
	results=matrix(0,nrow=length(fobs),ncol=nrep)
	for (i in 1:nrep){
	results[,i]<- Anova(g(permdata(data)),type=3)[,1]}
	permutationpvalues <- matrix(0,nrow=length(fobs),ncol=1)
	for (i in 1:length(fobs))
	{
	permutationpvalues[i] <- sum(results[i,] > fobs[i])/nrep}
	permutationpvalues}
	
#False discovery rate
fdr <- function (q,pvals) {
	
	cutoff <- apply(as.matrix(order(pvals,decreasing=TRUE)),1, function(x){(x*q)/length(pvals)})
	
	fdrvec <- pvals > cutoff
	apply(fdrvec,1,function(x){ifelse(x == TRUE,"","*")})

	
}

#Extract male fertility data from PolyPatEx output

sire.count <- function(x)
{
x$Mother <- factor(x$Mother,levels=levels(x$potentialFather))
x$Self <- x$Mother == x$potentialFather
x$minmis <- 0
for (i in 1:dim(x)[1]) {
	x[i,"minmis"] <- min(x[x$Progeny==x$Progeny[i],"VLTotal"]-x[x$Progeny==x$Progeny[i],"FLCount"])
}
x$incl <- 0
for (i in 1:dim(x)[1]) {
	x[i,"incl"] <- (x[i,"VLTotal"]- x[i,"FLCount"] - x[i,"minmis"])
}
x <- subset(x,incl==0)

x$repeats <- 0
for (i in 1:dim(x)[1]) {
	x[i,"repeats"] <- dim(x[(x$Progeny==x$Progeny[i]),])[1]
}
return(x)
}

male.fer <- function(x) {
	x <- x[x$repeats==1,]
	s <- sort(unique(x$potentialFather))
	
	mw_ns <- matrix(0,nrow=length(s),ncol=2)
	mw_ns <- as.data.frame(mw_ns)
	mw_ns[,1] <- s
	for (i in 1:length(s)) {
	mw_ns[i,2] <- length(x[(x$potentialFather==s[i])&(x$Self==FALSE | x$Self=="NA"),"repeats"])
	}
	return(mw_ns)
}

####Read in data
#Set working directory here:
wd="myfolder"
setwd(wd)

traitnames <- c("Tube length","Tube width","Petal length","Petal width","Fringe number","Stigama exsertion","Anther exsertion","Flower number","Display height")

d121m <- read.csv("2012_early_male.csv")
d121f <- read.csv("2012_early_female.csv")
d122m <- read.csv("2012_late_male.csv")
d122f <- read.csv("2012_late_female.csv")
d131m <- read.csv("2013_early_male.csv")
d131f <- read.csv("2013_early_female.csv")
d132m <- read.csv("2013_late_male.csv")
d132f <- read.csv("2013_late_female.csv")


### paternity analysis
pDataFile <- "genotype_2012_early.csv"
pData <- inputData(pDataFile,numLoci=8, ploidy=4,dataType="phenotype",dioecious=FALSE,selfCompatible=TRUE, lociMin = 1, matMismatches=1)
PPE121 <- phenotPPE(pData)
pf121  <- potentialFatherIDs(PPE121,mismatches=3,VLTMin=5)
pf121 <- sire.count(pf121)
mw121 <- male.fer(pf121)

pDataFile <- "genotype_2013_early.csv"
pData <- inputData(pDataFile,numLoci=8, ploidy=4,dataType="phenotype",dioecious=FALSE,selfCompatible=TRUE, lociMin = 1, matMismatches=1)
PPE131 <- phenotPPE(pData)
pf131  <- potentialFatherIDs(PPE131,mismatches=3,VLTMin=5)
pf131 <- sire.count(pf131)
mw131 <- male.fer(pf131)


pDataFile <- "genotype_2012_late.csv"
pData <- inputData(pDataFile,numLoci=8, ploidy=4,dataType="phenotype",dioecious=FALSE,selfCompatible=TRUE, lociMin = 1, matMismatches=1)
PPE122 <- phenotPPE(pData)
pf122  <- potentialFatherIDs(PPE122,mismatches=3,VLTMin=5)
pf122 <- sire.count(pf122)
mw122 <- male.fer(pf122)



pDataFile <- "genotype_2013_late.csv"
pData <- inputData(pDataFile,numLoci=8, ploidy=4,dataType="phenotype",dioecious=FALSE,selfCompatible=TRUE, lociMin = 1, matMismatches=1)
PPE132 <- phenotPPE(pData)
pf132  <- potentialFatherIDs(PPE132,mismatches=3,VLTMin=5)
pf132 <- sire.count(pf132)
mw132 <- male.fer(pf132)


d121m$ft <- lookup(as.factor(d121m$ID), mw121)
d121m$ft[is.na(d121m$ft)] = 0

d131m$ft <- lookup(as.factor(d131m$ID), mw131)
d131m$ft[is.na(d131m$ft)] = 0

d122m$ft <- lookup(as.factor(d122m$ID), mw122)
d122m$ft[is.na(d122m$ft)] = 0

d132m$ft <- lookup(as.factor(d132m$ID), mw132)
d132m$ft[is.na(d132m$ft)] = 0

d121f$rf = d121f$ft/mean(d121f$ft)
d131f$rf = d131f$ft/mean(d131f$ft)
d121m$rf = d121m$ft/mean(d121m$ft)
d131m$rf = d131m$ft/mean(d131m$ft)
d122f$rf = d122f$ft/mean(d122f$ft)
d132f$rf = d132f$ft/mean(d132f$ft)
d122m$rf = d122m$ft/mean(d122m$ft)
d132m$rf = d132m$ft/mean(d132m$ft)


########ANCOVA

p12 <- Anc.prep1 (d121m,d121f,d122m,d122f,id.tr)
p13 <- Anc.prep1 (d131m,d131f,d132m,d132f,id.tr)

models <- int(rbind(p12,p13))
write.csv(round(models,1),"aic.csv")

##Anova

out <- lmer(rf ~ (TL + TW.F + PL + PW + LB + STE + ANE + NF+ HT)  + Sex + E +Year + (TL + TW.F + PL + PW + LB + STE + ANE + NF+ HT)*Sex + (1|ID) ,REML=FALSE, data=rbind(p12,p13))
write.csv(round(anova(out),3),"ancova.csv")


###Pool data across years

dat_F_E=rbind(d121f,d131f)
dat_M_E=rbind(d121m,d131m)
dat_F_L=rbind(d122f,d132f)
dat_M_L=rbind(d122m,d132m)


###Selection analyses
#Pooled across years

nrep=1000

pvals_M_E=permP(dat_M_E,myp,nrep)
sresults_M_E = cbind(summary(myp(dat_M_E))$coefficients[2:10,1], pvals_M_E[2:10])

pvals_M_L=permP(dat_M_L,myp,nrep)
sresults_M_L = cbind(summary(myp(dat_M_L))$coefficients[2:10,1], pvals_M_L[2:10])


pvals_F_E =permP(dat_F_E,myp,nrep)
sresults_F_E = cbind(summary(myp(dat_F_E))$coefficients[2:10,1], pvals_F_E[2:10])

pvals_F_L =permP(dat_F_L,myp,nrep)
sresults_F_L = cbind(summary(myp(dat_F_L))$coefficients[2:10,1], pvals_F_L[2:10])

selectionresults_pooled_season <- cbind(sresults_M_E, sresults_F_E, sresults_M_L, sresults_F_L)
colnames(selectionresults_pooled_season) <- rep(c("beta","Pmu"),4)
rownames(selectionresults_pooled_season) <- traitnames

write.csv(round(selectionresults_pooled_season,3),"selectionresults_pooled_season.csv")


#Pooled across year & period + Quadratic terms

mypn <- function(x)
{
lmer(ft ~ scale(TL) + scale(TW.F)  + scale(PL)  + scale(PW) + scale(LB) + scale(STE) + scale(ANE) +  scale(NF) + scale(HT)  + I(scale(TL)^2) + I(scale(TW.F)^2)  + I(scale(PW)^2) + I(scale(PL)^2) + I(scale(LB)^2) + I(scale(STE)^2) + I(scale(ANE)^2) + I(scale(NF)^2) + I(scale(HT)^2) + (1|Year) + (1|E) ,data = x)
}


# results_f=myp(rbind(dat_F_E,dat_F_L))
# results_m=myp(rbind(dat_M_E,dat_M_L))

results_f=mypn(rbind(dat_F_E,dat_F_L))
results_m=mypn(rbind(dat_M_E,dat_M_L))

anova(results_f)
anova(results_m)

summary(results_f)
summary(results_m)

nrep=1000
pvals_F=permP(rbind(dat_F_E,dat_F_L),mypn,nrep)
sresults_F = cbind(round(summary(mypn(rbind(dat_F_E,dat_F_L)))$coefficients[2:19,1],3), paste(pvals_F,fdr(.10,pvals_F),sep=""))

pvals_M=permP(rbind(dat_M_E,dat_M_L),mypn, nrep)
sresults_M = cbind(round(summary(mypn(rbind(dat_M_E,dat_M_L)))$coefficients[2:19,1],3), paste(pvals_M,fdr(.1,pvals_M),sep=""))


selectionresults_pooled <- cbind(sresults_M, sresults_F)

write.csv(selectionresults_pooled,"selectionresults_pooled.csv")


#####Oviposition preferences

eggdat <- read.csv("eggdat.csv")
# ft = "number of eggs in flower"

results=lmer(ft/mean(ft) ~ scale(TL) + scale(TW.F)  + scale(PL)  + scale(PW) + scale(LB) + scale(STE) + scale(ANE) + (1|Year) +  (1| ID), data = eggdat, na.action=na.omit)
summary(results)$coefficients


lmegg <- function(x)
{
lmer(ft/mean(ft) ~ scale(TL) + scale(TW.F)  + scale(PL)  + scale(PW) + scale(LB) + scale(STE) + scale(ANE) + (1|Year) +  (1| ID), data = x, na.action=na.omit)}

anova(lmegg(eggdat),type=3)[,5]

nrep=1000
pval <- permP(eggdat,lmegg,nrep)
write.csv(cbind(round(summary(results)$coefficients[2:8,1],3), pval[2:8]),"oviposition_out.csv")

