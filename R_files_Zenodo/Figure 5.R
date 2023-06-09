### This code produces Figure 5 in the manuscript: 
### Burgess SC, Powell J, Bueno M. Dispersal, kin aggregation, and the fitness consequences of not spreading sibling larvae.
# Code modified Jun 2022 to increase gap between panels
# Any comments or error reporting, please contact Scott Burgess. sburgess@bio.fsu.edu

# R version 4.0.5 (2021-03-31) -- "Shake and Throw"

# Load required libraries
library(lme4) # version 1.1-27.1 


# Import data
# set.wd() # set working directory, or file path, before importing data
dat <- read.csv("Figure 5 data.csv")

# Prepare data
dat$sib.total <- as.factor(dat$sib.total)
dat$sib.group.size <- as.factor(dat$sib.group.size)
dat$total.group.size <- as.factor(dat$total.group.size)
dat$mother.ID <- as.factor(dat$mother.ID)
dat$sheet <- as.factor(dat$sheet)
# head(dat)
# summary(dat)

# get focal individuals only
datF <- dat[dat$Ind=="F",] 
# head(datF)

# How many individuals?
dim(dat[dat$age=="10",])[1] # 280 individuals total
dim(datF[datF$age=="10",])[1] # 120 focal individuals total

# How many groups?
length(unique(dat$sheet)) # 88 groups

# How many mothers?
length(unique(dat$mother.ID)) # 4 families
with(datF[datF$age=="10",],table(mother.ID)) # 30 settlers from each maternal family

# When were data collected
unique(dat$age.d) # 10, 24, 38 days since settlement

# Subset data
# Age 10 days
d10 <- datF[datF$age.d=="10",]
d10$sib.total <- factor(d10$sib.total,levels=c("11","15","55"))
# Age 24 days
d24 <- datF[datF$age.d=="24",]
d24$sib.total <- factor(d24$sib.total,levels=c("11","15","55"))
# Age 38 days
d38 <- datF[datF$age.d=="38",]
d38$sib.total <- factor(d38$sib.total,levels=c("11","15","55"))


## Fit models

## Growth 
# = Size (# bifurcations) at 24 days
m1s <- glmer(bifurcations ~ sib.total + (1|mother.ID), data=d24, family="poisson")
m2s <- glmer(bifurcations ~ 1 + (1|mother.ID), data=d24, family="poisson")
round(exp(fixef(m1s)[1]),2) 
round(exp(confint(m1s)[2,]),2) 
# At 24 days, colonies averaged 2.99 (2.07 - 4.14) bifurcations at low density

exp(fixef(m1s)[3])
exp(confint(m1s)[4,]) # (0.39 - 0.88, 95% CI)
1 - exp(confint(m1s)[4,]) # (12 - 60, 95% CI)
anova(m1s,m2s,test="Chisq") # X2=7.27, df=2, p=0.026)
# The combination of high density and high relatedness reduced the number of bifurcations by 41% (12 – 60, 95% CI) compared to that at low density (Χ^2=7.27, df=2, p=0.026).


## Survival
# 10 days
m1.10 <- glmer(survival ~ sib.total + (1|mother.ID), data=d10, family="binomial",na.action=na.exclude)
m2.10 <- glmer(survival ~ 1 + (1|mother.ID), data=d10, family="binomial",na.action=na.exclude)
anova(m1.10,m2.10,test="Chisq") # Survival at 10 days, X2=2.01, df=2, p=0.37
plogis(fixef(m2.10)) # 75% survival
plogis(confint(m2.10)[2,]) # (52 - 90, 95% CI)

## 24 days
m1.24 <- glmer(survival ~ sib.total + (1|mother.ID), data=d24, family="binomial",na.action=na.exclude)
m2.24 <- glmer(survival ~ 1 + (1|mother.ID), data=d24, family="binomial",na.action=na.exclude)
anova(m1.24,m2.24,test="Chisq") # Survival at 24 days, X2=0.34, df=2, p=0.85
plogis(fixef(m2.24)) # 55% survival
plogis(confint(m2.24)[2,]) # (43 - 0.68, 95% CI)

## 38 days
m1.38 <- glmer(survival ~ sib.total + (1|mother.ID), data=d38, family="binomial",na.action=na.exclude)
m2.38 <- glmer(survival ~ 1 + (1|mother.ID), data=d38, family="binomial",na.action=na.exclude)
anova(m1.38,m2.38,test="Chisq") # Survival at 38 days, X2=6.11, df=2, p=0.047
plogis(fixef(m1.38)[1]) # 45% survival at low density
plogis(confint(m1.38)[2,]) # (25 - 0.66, 95% CI)

d38$sib.total <- factor(d38$sib.total,levels=c("15","55","11"))
m <- glmer(survival ~ sib.total + (1|mother.ID), data=d38, family="binomial",na.action=na.exclude)
exp(fixef(m)[2]) # Odds of survival changed by factor of 0.31 
exp(confint(m)[3,]) # (0.12 - 0.80, 95% CI) 
1-exp(fixef(m)[2]) # Odds of survival were were 69% lower
1-exp(confint(m)[3,]) # (20 - 88%, 95%CI)
d38$sib.total <- factor(d38$sib.total,levels=c("11","15","55"))




# Prepare model predictions for plot
pred.frame <- with(d24, expand.grid(sib.total=levels(sib.total))) 
X <- model.matrix(~sib.total,data=pred.frame) 
pred <- data.frame(pred.frame,mean=(X%*%fixef(m1s)))
V <- vcov(m1s)
# Fixed effects uncertainty only
pred.var1 <- diag(X %*% V %*% t(X)) # XVX^T
# diag(X %*% tcrossprod(V,X)) # just another way

# Attach to dataframe, calculate standard errors and confidence intervals 
predictionsSize <- data.frame(pred,pred.se1=sqrt(pred.var1))
predictionsSize <- with(predictionsSize, data.frame(predictionsSize,
	p1.lo = mean-1.96*pred.se1,
	p1.hi = mean+1.96*pred.se1))


pred.frame <- with(d38, expand.grid(sib.total=levels(sib.total))) 
X <- model.matrix(~sib.total,data=pred.frame) 
pred <- data.frame(pred.frame,mean=(X%*%fixef(m1.38)))
V <- vcov(m1.38)
# Fixed effects uncertainty only
pred.var1 <- diag(X %*% V %*% t(X)) # XVX^T
# Fixed effects uncertainty + random effects variance
pred.var2 <- pred.var1 + as.numeric(VarCorr(m1.38)$mother.ID[1]) 

# Attach to dataframe, calculate standard errors and confidence intervals
predictions <- data.frame(pred,pred.se1=sqrt(pred.var1))
predictions <- with(predictions, data.frame(predictions,
	p1.lo = mean-1.96*pred.se1,
	p1.hi = mean+1.96*pred.se1))



# Make Figure 5
x.vec <- 1:3
cx <- 1.5
mID <- unique(dat$mother.ID)

xlabs=c("1 sib alone\\n\\nlow density\\n","1 sib +\\n4 unrelated\\nhigh density\\nlow relatedness","5 sibs\\n\\nhigh density\\nhigh relatedness")

quartz(width=10,height=5.2)
par(mar=c(2,3.5,2,2),oma=c(3,2,0,0),mfrow=c(1,2))

# a)
plot(x.vec,exp(predictionsSize$mean),type="n",ylim=c(0,5),xlim=c(0.9,3.1),ylab="",xlab="",xaxt="n",yaxt="n",bty="l")
axis(side=1,at=c(1,1.9,2.8),labels=xlabs,tick=F,line=2,cex.axis=1)
axis(side=2,at=seq(0,5,1),labels=seq(0,5,1),las=1,cex.axis=1.4)
mtext("Size (# bifurcations at 24 days)",side=2,line=2.5,cex=cx)

set.seed(99)
with(d24,points(jitter(as.numeric(sib.total),0.5), bifurcations,
		pch=19,
		cex=1.2,
		# col=adjustcolor("grey",alpha.f=0.6)
		col="grey"
		))

for(i in 1:length(mID)){
		y.vec <- predictionsSize$mean + ranef(m1s)$mother.ID[i,]
		lines(x.vec,exp(y.vec),col="grey",lwd=2)	
}
points(x.vec,exp(predictionsSize$mean),pch=19,cex=2)
segments(x.vec,exp(predictionsSize$p1.lo),x.vec,exp(predictionsSize$p1.hi),lwd=2)
mtext("a)",side=3,cex=cx,adj=0)

# b)
plot(x.vec,predictions$mean,type="n",ylim=c(0,1),xlim=c(0.9,3.1),ylab="",xlab="",xaxt="n",yaxt="n",bty="l")
axis(side=1,at=c(1,1.9,2.8),labels=xlabs,tick=F,line=2,cex.axis=1)
axis(side=2,at=seq(0,1,0.1),labels=seq(0,1,0.1),las=1,cex.axis=1.4)
mtext("Probability of survival (38 days)",side=2,line=3.5,cex=cx,adj=0.8)

for(i in 1:length(mID)){
		y.vec <- predictions$mean + ranef(m1.38)$mother.ID[i,]
		lines(x.vec,plogis(y.vec),col="grey",lwd=2)	
}
points(x.vec,plogis(predictions$mean),pch=19,cex=2)
segments(x.vec,plogis(predictions$p1.lo),x.vec,plogis(predictions$p1.hi),lwd=2)
mtext("b)",side=3,cex=cx,adj=0)


