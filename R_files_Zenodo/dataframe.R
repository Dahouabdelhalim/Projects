
# Script to reproduce the data used in:
#
# Cowling BJ, Ng DMW, Ip DKM, Liao Q, Lam WWT, Wu JT, et al.
# Community Psychological and Behavioral Responses through the First Wave
# of the 2009 Influenza A(H1N1) Pandemic in Hong Kong
# JID 2010; 202(6):867-876.




require(Design)

survey <- read.csv("H1N1_surveys.csv")
subset <- survey[survey$wave>2&survey$wave!=6,]

# exclude fever+cough OR fever+sore throat
subset$sick <- 1*((subset$sm1_1==1&subset$sm1_5==1)|(subset$sm1_1==1&subset$sm1_9==1))
table (subset$sick, exclude=NULL)

subset <- subset[subset$sick==0&!is.na(subset$sick),]

anxiety <- cbind(subset$ax1_1, subset$ax1_2, subset$ax1_3, subset$ax1_4, subset$ax1_5, subset$ax1_6, subset$ax1_7, subset$ax1_8, subset$ax1_9, subset$ax1_10)

withNA <- apply(anxiety, 1, function(x) ifelse(sum(is.na(x))>=1,1,0))
sum(withNA)

subset$bf3_1[!is.na(subset$bf3_1)&subset$bf3_1==3] <- NA
subset$bf3_3[!is.na(subset$bf3_3)&subset$bf3_3==3] <- NA
subset$bf3_4[!is.na(subset$bf3_4)&subset$bf3_4==3] <- NA
subset$pm2[!is.na(subset$pm2)&subset$pm2==5] <- NA
subset$pm3[!is.na(subset$pm3)&subset$pm3==5] <- NA
subset$pm3a[!is.na(subset$pm3a)&subset$pm3a==5] <- NA
subset$pm4[!is.na(subset$pm4)&subset$pm4==5] <- NA
subset$pm5[!is.na(subset$pm5)&subset$pm5==5] <- NA
subset$pm7[!is.na(subset$pm7)&subset$pm7==5] <- NA
subset$pm7b[!is.na(subset$pm7b)&subset$pm7b==5] <- NA
subset$pm10_1[!is.na(subset$pm10_1)&subset$pm10_1==4] <- NA
subset$pm10_2[!is.na(subset$pm10_2)&subset$pm10_2==4] <- NA
subset$pm10_3[!is.na(subset$pm10_3)&subset$pm10_3==4] <- NA
subset$pm10_5[!is.na(subset$pm10_5)&subset$pm10_5==4] <- NA
subset$pm10_6[!is.na(subset$pm10_6)&subset$pm10_6==4] <- NA

# Multiple imputation
set.seed(12345)

subset.i <- aregImpute( ~ I(pm3)+I(pm5)+sex+I(age_gp)+I(edu)+I(ax1_1)+I(ax1_2)+I(ax1_3)+I(ax1_4)+I(ax1_5)+I(ax1_6)+I(ax1_7)+I(ax1_8)+I(ax1_9)+I(ax1_10)+I(ph1)+I(bf1)+I(bf2)+I(pm2)+I(pm3a)+I(pm4)+I(pm7)+I(pm7b)+I(pm10_1)+I(pm10_2)+I(pm10_3)+I(pm10_5)+I(pm10_6)+I(bf5)+I(bf4)+I(bf3_1)+I(bf3_3)+I(bf3_4)+wave,
data=subset, n.impute=10)

subset.nomiss <- list(subset, subset, subset, subset, subset, subset, subset, subset, subset, subset)
for(i in 1:10){
    subset.nomiss[[i]]$pm3[is.na(subset.nomiss[[i]]$pm3)] <- subset.i$imputed$pm3[,i]
    subset.nomiss[[i]]$pm3a[is.na(subset.nomiss[[i]]$pm3a)] <- subset.i$imputed$pm3a[,i]
    subset.nomiss[[i]]$pm4[is.na(subset.nomiss[[i]]$pm4)] <- subset.i$imputed$pm4[,i]
    subset.nomiss[[i]]$pm5[is.na(subset.nomiss[[i]]$pm5)] <- subset.i$imputed$pm5[,i]
    subset.nomiss[[i]]$pm7[is.na(subset.nomiss[[i]]$pm7)] <- subset.i$imputed$pm7[,i]
    subset.nomiss[[i]]$pm7b[is.na(subset.nomiss[[i]]$pm7b)] <- subset.i$imputed$pm7b[,i]
    subset.nomiss[[i]]$sex[is.na(subset.nomiss[[i]]$sex)] <- subset.i$imputed$sex[,i]
    subset.nomiss[[i]]$age_gp[is.na(subset.nomiss[[i]]$age_gp)] <- subset.i$imputed$age_gp[,i]
    subset.nomiss[[i]]$edu[is.na(subset.nomiss[[i]]$edu)] <- subset.i$imputed$edu[,i]
    subset.nomiss[[i]]$ax1_1[is.na(subset.nomiss[[i]]$ax1_1)] <- subset.i$imputed$ax1_1[,i]
    subset.nomiss[[i]]$ax1_2[is.na(subset.nomiss[[i]]$ax1_2)] <- subset.i$imputed$ax1_2[,i]
    subset.nomiss[[i]]$ax1_3[is.na(subset.nomiss[[i]]$ax1_3)] <- subset.i$imputed$ax1_3[,i]
    subset.nomiss[[i]]$ax1_4[is.na(subset.nomiss[[i]]$ax1_4)] <- subset.i$imputed$ax1_4[,i]
    subset.nomiss[[i]]$ax1_5[is.na(subset.nomiss[[i]]$ax1_5)] <- subset.i$imputed$ax1_5[,i]
    subset.nomiss[[i]]$ax1_6[is.na(subset.nomiss[[i]]$ax1_6)] <- subset.i$imputed$ax1_6[,i]
    subset.nomiss[[i]]$ax1_7[is.na(subset.nomiss[[i]]$ax1_7)] <- subset.i$imputed$ax1_7[,i]
    subset.nomiss[[i]]$ax1_8[is.na(subset.nomiss[[i]]$ax1_8)] <- subset.i$imputed$ax1_8[,i]
    subset.nomiss[[i]]$ax1_9[is.na(subset.nomiss[[i]]$ax1_9)] <- subset.i$imputed$ax1_9[,i]
    subset.nomiss[[i]]$ax1_10[is.na(subset.nomiss[[i]]$ax1_10)] <- subset.i$imputed$ax1_10[,i]
    subset.nomiss[[i]]$ph1[is.na(subset.nomiss[[i]]$ph1)] <- subset.i$imputed$ph1[,i]
    subset.nomiss[[i]]$bf1[is.na(subset.nomiss[[i]]$bf1)] <- subset.i$imputed$bf1[,i]
    subset.nomiss[[i]]$bf2[is.na(subset.nomiss[[i]]$bf2)] <- subset.i$imputed$bf2[,i]
    subset.nomiss[[i]]$pm10_1[is.na(subset.nomiss[[i]]$pm10_1)] <- subset.i$imputed$pm10_1[,i]
    subset.nomiss[[i]]$pm10_2[is.na(subset.nomiss[[i]]$pm10_2)] <- subset.i$imputed$pm10_2[,i]
    subset.nomiss[[i]]$pm10_3[is.na(subset.nomiss[[i]]$pm10_3)] <- subset.i$imputed$pm10_3[,i]
    subset.nomiss[[i]]$pm10_5[is.na(subset.nomiss[[i]]$pm10_5)] <- subset.i$imputed$pm10_5[,i]
    subset.nomiss[[i]]$pm10_6[is.na(subset.nomiss[[i]]$pm10_6)] <- subset.i$imputed$pm10_6[,i]
    subset.nomiss[[i]]$bf5[is.na(subset.nomiss[[i]]$bf5)] <- subset.i$imputed$bf5[,i]
    subset.nomiss[[i]]$bf4[is.na(subset.nomiss[[i]]$bf4)] <- subset.i$imputed$bf4[,i]
    subset.nomiss[[i]]$bf3_1[is.na(subset.nomiss[[i]]$bf3_1)] <- subset.i$imputed$bf3_1[,i]
    subset.nomiss[[i]]$bf3_3[is.na(subset.nomiss[[i]]$bf3_3)] <- subset.i$imputed$bf3_3[,i]
    subset.nomiss[[i]]$bf3_4[is.na(subset.nomiss[[i]]$bf3_4)] <- subset.i$imputed$bf3_4[,i]

}

for (i in 1:10){
# anxiety
subset.nomiss[[i]]$ax1_1 <- 5-subset.nomiss[[i]]$ax1_1
subset.nomiss[[i]]$ax1_2 <- 5-subset.nomiss[[i]]$ax1_2
subset.nomiss[[i]]$ax1_3 <- 5-subset.nomiss[[i]]$ax1_3
subset.nomiss[[i]]$ax1_4 <- 5-subset.nomiss[[i]]$ax1_4
subset.nomiss[[i]]$ax1_5 <- 5-subset.nomiss[[i]]$ax1_5

subset.nomiss[[i]]$ax1_score <- (subset.nomiss[[i]]$ax1_1+subset.nomiss[[i]]$ax1_2+subset.nomiss[[i]]$ax1_3+subset.nomiss[[i]]$ax1_4+subset.nomiss[[i]]$ax1_5+subset.nomiss[[i]]$ax1_6+subset.nomiss[[i]]$ax1_7+subset.nomiss[[i]]$ax1_8+subset.nomiss[[i]]$ax1_9+subset.nomiss[[i]]$ax1_10)/10

subset.nomiss[[i]]$ax1_score[subset.nomiss[[i]]$ax1_score<2] <- 1
subset.nomiss[[i]]$ax1_score[subset.nomiss[[i]]$ax1_score>=2&subset.nomiss[[i]]$ax1_score<2.5] <- 2
subset.nomiss[[i]]$ax1_score[subset.nomiss[[i]]$ax1_score>=2.5&subset.nomiss[[i]]$ax1_score<=4] <- 3

# sex ref:M
subset.nomiss[[i]]$sex[subset.nomiss[[i]]$sex==1] <- 0
subset.nomiss[[i]]$sex[subset.nomiss[[i]]$sex==2] <- 1

# age ref: 35-44
subset.nomiss[[i]]$age_gp[subset.nomiss[[i]]$age_gp==3] <- 0

# education ref:primary
subset.nomiss[[i]]$edu[subset.nomiss[[i]]$edu<=3] <- 1
subset.nomiss[[i]]$edu[subset.nomiss[[i]]$edu<=6&subset.nomiss[[i]]$edu>=4] <- 2
subset.nomiss[[i]]$edu[subset.nomiss[[i]]$edu==7] <- 3

# perceived health
subset.nomiss[[i]]$ph1[subset.nomiss[[i]]$ph1==3] <- 0

# absolute susceptibility
subset.nomiss[[i]]$bf1[subset.nomiss[[i]]$bf1==4] <- 0
subset.nomiss[[i]]$bf1[subset.nomiss[[i]]$bf1==7] <- 6

# relative susceptibility
subset.nomiss[[i]]$bf2[subset.nomiss[[i]]$bf2==4] <- 0
subset.nomiss[[i]]$bf2[subset.nomiss[[i]]$bf2==7] <- 6

# severity vs SARS
subset.nomiss[[i]]$bf5[subset.nomiss[[i]]$bf5==3] <- 0

# how worry if developed ILI tomorrow
subset.nomiss[[i]]$bf4[subset.nomiss[[i]]$bf4==4] <- 0

# knowledge
subset.nomiss[[i]]$bf3_1[subset.nomiss[[i]]$bf3_1==2] <- 0
subset.nomiss[[i]]$bf3_3[subset.nomiss[[i]]$bf3_3==2] <- 0
subset.nomiss[[i]]$bf3_4[subset.nomiss[[i]]$bf3_4==2] <- 0

# cover mouth
subset.nomiss[[i]]$pm2[subset.nomiss[[i]]$pm2<=2] <- 1
subset.nomiss[[i]]$pm2[subset.nomiss[[i]]$pm2>=3&subset.nomiss[[i]]$pm2<=4] <- 0

# handwashing after sneezing
subset.nomiss[[i]]$pm3[subset.nomiss[[i]]$pm3<=2] <- 1
subset.nomiss[[i]]$pm3[subset.nomiss[[i]]$pm3>=3&subset.nomiss[[i]]$pm3<=4] <- 0

# handwashing after home
subset.nomiss[[i]]$pm3a[subset.nomiss[[i]]$pm3a<=2] <- 1
subset.nomiss[[i]]$pm3a[subset.nomiss[[i]]$pm3a>=3&subset.nomiss[[i]]$pm3a<=4] <- 0

# handwashing - use liquid soup
subset.nomiss[[i]]$pm4[subset.nomiss[[i]]$pm4<=2] <- 1
subset.nomiss[[i]]$pm4[subset.nomiss[[i]]$pm4>=3&subset.nomiss[[i]]$pm4<=4] <- 0

# wear mask
subset.nomiss[[i]]$pm5[subset.nomiss[[i]]$pm5<=2] <- 1
subset.nomiss[[i]]$pm5[subset.nomiss[[i]]$pm5>=3&subset.nomiss[[i]]$pm5<=4] <- 0

# adopt any preventive measures when touching common objects (pm7)
subset.nomiss[[i]]$pm7[subset.nomiss[[i]]$pm7<=2] <- 1
subset.nomiss[[i]]$pm7[subset.nomiss[[i]]$pm7>=3&subset.nomiss[[i]]$pm7<=4] <- 0

# wash hands after toucning common objects (pm7b)
subset.nomiss[[i]]$pm7b[subset.nomiss[[i]]$pm7b<=2] <- 1
subset.nomiss[[i]]$pm7b[subset.nomiss[[i]]$pm7b>=3&subset.nomiss[[i]]$pm7b<=4] <- 0

# distancing
subset.nomiss[[i]]$pm10_1[subset.nomiss[[i]]$pm10_1>=2&subset.nomiss[[i]]$pm10_1<=3] <- 0
subset.nomiss[[i]]$pm10_2[subset.nomiss[[i]]$pm10_2>=2&subset.nomiss[[i]]$pm10_2<=3] <- 0
subset.nomiss[[i]]$pm10_3[subset.nomiss[[i]]$pm10_3>=2&subset.nomiss[[i]]$pm10_3<=3] <- 0
subset.nomiss[[i]]$pm10_5[subset.nomiss[[i]]$pm10_5>=2&subset.nomiss[[i]]$pm10_5<=3] <- 0
subset.nomiss[[i]]$pm10_6[subset.nomiss[[i]]$pm10_6>=2&subset.nomiss[[i]]$pm10_6<=3] <- 0
}

combine.mi <- function(model, n.impute){
	betas <- matrix(c(model[[1]][[4]]$fixed, model[[2]][[4]]$fixed, model[[3]][[4]]$fixed,
           model[[4]][[4]]$fixed, model[[5]][[4]]$fixed,model[[6]][[4]]$fixed, model[[7]][[4]]$fixed,
           model[[8]][[4]]$fixed, model[[9]][[4]]$fixed, model[[10]][[4]]$fixed),
           byrow=FALSE, ncol=n.impute)  # beta
	vars <- matrix(c(diag(model[[1]][[5]]), diag(model[[2]][[5]]), diag(model[[3]][[5]]),
          diag(model[[4]][[5]]), diag(model[[5]][[5]]),diag(model[[6]][[5]]), diag(model[[7]][[5]]),
          diag(model[[8]][[5]]), diag(model[[9]][[5]]), diag(model[[10]][[5]])),
          byrow=FALSE, ncol=n.impute)
	coef.names <- names(model[[1]][[4]]$fixed)
	mean.coefs <- rowMeans(betas)
	Ubar <- rowMeans(vars)
	B <- rowSums((betas - mean.coefs)*(betas-mean.coefs) /
		(n.impute - 1))
	T <- (1 + 1/n.impute) * B + Ubar
	degf <- (n.impute - 1)*(1 + Ubar / ((1 + 1/n.impute)*B))*
		(1 + Ubar / ((1 + 1/n.impute)*B))
	data.frame(OR = exp(mean.coefs),
		lowerCI = exp(mean.coefs - qt(0.975, df=degf)*sqrt(T))
		,upperCI = exp(mean.coefs + qt(0.975, df=degf)*sqrt(T))
		#,p.value = 2*(1 - pt(abs(mean.coefs)/sqrt(T), df=degf))
    ,row.names=coef.names)
}

# End of script
