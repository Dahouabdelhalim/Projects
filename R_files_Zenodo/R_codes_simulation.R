library(clusterGeneration)
library(evolqg)




##Generation of variance/covariance (V/CV) matrix with 300 dimensions (Grabowski and Porto, 2017)
ran.gen.300=genPositiveDefMat("c-vine", dim=300, rangeVar=c(0.5, 0.6))


###Mean and SD of ICV sampling sample size 11-150 with r2 value from 0.02 to 0.7
##Means of ICV resampling with 11-150 sample sizes
ICVm.300 = function(p){
decom=svd(ran.gen.300$Sigma)
recom=c(p, rep(1, 299))* decom$d
ran.gen.de= decom$u %*% diag(recom) %*% t(decom$v)
@p = multiplication number to the first eigenvalue. 

Sigma = ran.gen.de
random_data = mvrnorm(n = 10000, rep(0, 300), Sigma, empirical=FALSE)

random.sam = function(n,t,r){
as.numeric(lapply(1:r,function(i) CalcICV(cov(random_data[sample(1:10000, n), sample(1:300, t)]))))
}
  @n = number of individuals resampled
  @t = number of traits resampled
  @r = number of resampling
  
samplem=sapply(1:140, function(i) mean(random.sam(seq(from=11, to=150, by=1)[i], 10, 1000)))
return(samplem)
}

vary.r1.m.300 = sapply(1:5, function(i) ICVm.300(seq(from=15, to=75, by=12)[i]))
vary.r2.m.300 = sapply(1:5, function(i) ICVm.300(seq(from=78, to=178, by=18)[i]))
vary.r3.m.300 = sapply(1:5, function(i) ICVm.300(seq(from=175, to=295, by=30)[i]))
vary.r4.m.300 = sapply(1:5, function(i) ICVm.300(seq(from=340, to=590, by=60)[i]))
vary.r5.m.300 = sapply(1:5, function(i) ICVm.300(seq(from=665, to=1165, by=110)[i]))
vary.r6.m.300 = sapply(1:5, function(i) ICVm.300(seq(from=1200, to=2300, by=220)[i]))
vary.r7.m.300 = sapply(1:5, function(i) ICVm.300(seq(from=2450, to=4950, by=500)[i]))

vary.m.300=cbind(vary.r1.m.300, vary.r2.m.300, vary.r3.m.300, vary.r4.m.300, vary.r5.m.300, vary.r6.m.300, vary.r7.m.300)

boxplot(vary.m.300, ylab="Mean ICV", xlab="r2", 
        names=c(varyr2.300, varyr3.300, varyr4.300, varyr5.300, varyr6.300, varyr7.300, varyr8.300), cex.axis=1.5, cex.lab=1.5)



##SD of means of ICV resampling with 11-150 sample sizes
sd.m.300=function(x){
sd(vary.m.300[,x])
}

vary.m.sd.300=as.numeric(lapply(1:35, function(i) sd.m.300(i)))

vary.r2.300=c(varyr2.300, varyr3.300, varyr4.300, varyr5.300, varyr6.300, varyr7.300, varyr8.300)

plot(vary.r2.300, vary.m.sd.300, ylab="Standard deviation of Mean ICV", xlab="r2", 
     cex.axis=1.5, cex.lab=1.5)


###Calculation of means of ICV calculated by resampling with 11-150 sample sizes based on specific parameter r2 value
#matrix decomposition
decom = svd(ran.gen.300$Sigma)
recom=c(p, rep(1, 299))* decom$d
ran.gen.de= decom $u %*% diag(recom) %*% t(decom$v)
Sigma = ran.gen.de
@p = multiplication number to the first eigenvalue.

#Generation of multivariate 10000 multivariate normal population based on parameter V/CV matrix 
random_data = mvrnorm(n = 10000, rep(0, 300), Sigma, empirical=FALSE)

# ICV calculation with 10 traits resampling out of 300 traits with varying sample sizes from multivariate normal population
random.sam = function(n,t,r){
as.numeric(lapply(1:r,function(i) CalcICV(cov(random_data[sample(1:10000, n), sample(1:300, t)]))))
}
@n = number of individuals resampled
@t = number of traits resampled
@r = number of resampling

samplem =sapply(1:140, function(i) mean(random.sam(seq(from=11, to=150, by=1)[i], 10, 1000)))


###Calculation of means of ICV calculated by resampling method with 11-150 sample sizes when r2=0.05 
#matrix decomposition
decom1=svd(ran.gen.300$Sigma)
recom1=c(32, rep(1, 299))* decom1 $d
ran.gen.de1= decom1 $u %*% diag(recom1) %*% t(decom1 $v)

#Generation of multivariate 10000 multivariate normal population based on parameter V/CV matrix 
ran.sam1= function(x){
Sigma1=ran.gen.de1
mvrnorm(n = x, rep(0, 300), Sigma1, empirical=FALSE)
}

random_data1=ran.sam1(10000)

# ICV calculation with 10 traits resampling out of 300 traits with varying sample sizes from multivariate normal population
random.sam1 = function(n,t,r){
as.numeric(lapply(1:r,function(i) CalcICV(cov(random_data1[sample(1:10000, n), sample(1:300, t)]))))
}

samplem1=sapply(1:140, function(i) random.sam1(seq(from=11, to=150, by=1)[i], 10, 1000))

ss=seq(from=11, to=150, by=1)

boxplot(samplem1, names=ss, ylab="ICV", xlab="sample size", cex.axis=1.5, cex.lab=1.5)



###Calculation of means of ICV calculated by resampling method with 11-150 sample sizes when r2=0.08.
#matrix decomposition
decom2=svd(ran.gen.300$Sigma)
recom2=c(50.5, rep(1, 299))* decom2 $d
ran.gen.de2= decom2 $u %*% diag(recom2) %*% t(decom2 $v)

#Generation of multivariate 10000 multivariate normal population based on parameter V/CV matrix 
ran.sam2= function(x){
Sigma2=ran.gen.de2
mvrnorm(n = x, rep(0, 300), Sigma2, empirical=FALSE)
}

random_data2=ran.sam2(10000)

# ICV calculation with 10 traits resampling out of 300 traits with varying sample sizes from multivariate normal population
random.sam2 = function(n,t,r){
as.numeric(lapply(1:r,function(i) CalcICV(cov(random_data2[sample(1:10000, n), sample(1:300, t)]))))
}

Samplem2=sapply(1:140, function(i) random.sam2(seq(from=11, to=150, by=1)[i], 10, 1000))

ss=seq(from=11, to=150, by=1)

boxplot(samplem2, names=ss, ylab="ICV", xlab="sample size", cex.axis=1.5, cex.lab=1.5)



###Calculation of means of ICV calculated by resampling method with 11-150 sample sizes when r2=0.12
#matrix decomposition
decom3=svd(ran.gen.300$Sigma)
recom3=c(78, rep(1, 299))* decom3 $d
ran.gen.de3= decom3 $u %*% diag(recom3) %*% t(decom3 $v)

#Generation of multivariate 10000 multivariate normal population based on parameter V/CV matrix 
ran.sam3= function(x){
Sigma3=ran.gen.de3
mvrnorm(n = x, rep(0, 300), Sigma3, empirical=FALSE)
}

random_data3=ran.sam3(10000)

# ICV calculation with 10 traits resampling out of 300 traits with varying sample sizes from multivariate normal population
random.sam3 = function(n,t,r){
as.numeric(lapply(1:r,function(i) CalcICV(cov(random_data3[sample(1:10000, n), sample(1:300, t)]))))
}

samplem3=sapply(1:140, function(i) random.sam3(seq(from=11, to=150, by=1)[i], 10, 1000))

ss=seq(from=11, to=150, by=1)

boxplot(samplem3, names=ss, ylab="ICV", xlab="sample size", cex.axis=1.5, cex.lab=1.5)



###Calculation of means of ICV calculated by resampling method with 11-150 sample sizes when r2=0.2 #matrix decomposition
decom4=svd(ran.gen.300$Sigma)
recom4=c(150, rep(1, 299))* decom4 $d
ran.gen.de4= decom4 $u %*% diag(recom4) %*% t(decom4 $v)

#Generation of multivariate 10000 multivariate normal population based on parameter V/CV matrix 
ran.sam4= function(x){
Sigma4=ran.gen.de4
mvrnorm(n = x, rep(0, 300), Sigma4, empirical=FALSE)
}

random_data4=ran.sam4(10000)

# ICV calculation with 10 traits resampling out of 300 traits with varying sample sizes from multivariate normal population
random.sam4 = function(n,t,r){
as.numeric(lapply(1:r,function(i) CalcICV(cov(random_data4[sample(1:10000, n), sample(1:300, t)]))))
}

samplem4=sapply(1:140, function(i) random.sam4(seq(from=11, to=150, by=1)[i], 10, 1000))

ss=seq(from=11, to=150, by=1)

boxplot(samplem4, names=ss, ylab="ICV", xlab="sample size", cex.axis=1.5, cex.lab=1.5)


###Calculation of means of ICV calculated by resampling method with 11-150 sample sizes when r2=0.35 
#matrix decomposition
decom5=svd(ran.gen.300$Sigma)
recom5=c(385, rep(1, 299))* decom5 $d
ran.gen.de5= decom5 $u %*% diag(recom5) %*% t(decom5 $v)

#Generation of multivariate 10000 multivariate normal population based on parameter V/CV matrix 
ran.sam5= function(x){
Sigma5=ran.gen.de5
mvrnorm(n = x, rep(0, 300), Sigma5, empirical=FALSE)
}

random_data5=ran.sam5(10000)

# ICV calculation with 10 traits resampling out of 300 traits with varying sample sizes from multivariate normal population
random.sam5 = function(n,t,r){
as.numeric(lapply(1:r,function(i) CalcICV(cov(random_data5[sample(1:10000, n), sample(1:300, t)]))))
}

samplem5=sapply(1:140, function(i) random.sam5(seq(from=11, to=150, by=1)[i], 10, 1000))

ss=seq(from=11, to=150, by=1)

boxplot(samplem5, names=ss, ylab="ICV", xlab="sample size", cex.axis=1.5, cex.lab=1.5)



###Effect of resampled trait numbers and the number of total traits on ICV
##300 of total traits and 30 sample size.
trait.300.sam.30 = function(p){
decom300.30=svd(ran.gen.300$Sigma)
recom300.30=c(p, rep(1, 299))* decom300.30$d
ran.gen.de300.30= decom300.30$u %*% diag(recom300.30) %*% t(decom300.30$v)
Sigma300.30 = ran.gen.de300.30
random_data300.30 = mvrnorm(n = 10000, rep(0, 300), Sigma300.30, empirical=FALSE)
random.sam300.30 = function(n,t,r){
as.numeric(lapply(1:r,function(i) CalcICV(cov(random_data300.30[sample(1:10000, n), sample(1:300, t)]))))
}
samplem300.30=sapply(1:20, function(i) mean(random.sam300.30 (30, seq(from=10, to=29, by=1)[i], 1000)))
return(samplem300.30)
}

#r2=0.1
t300s30r0.1 = trait.300.sam.30(63.5) 
#r2=0.5
t300s30r0.5 = trait.300.sam.30(926)


##300 of total traits and 100 sample size.
trait.300.sam.100 = function(p){
decom300.100=svd(ran.gen.300$Sigma)
recom300.100=c(p, rep(1, 299))* decom300.100$d
ran.gen.de300.100= decom300.100$u %*% diag(recom300.100) %*% t(decom300.100$v)
  
Sigma300.100 = ran.gen.de300.100
random_data300.100 = mvrnorm(n = 10000, rep(0, 300), Sigma300.100, empirical=FALSE)
random.sam300.100 = function(n,t,r){
as.numeric(lapply(1:r,function(i) CalcICV(cov(random_data300.100 [sample(1:10000, n), sample(1:300, t)]))))
}
samplem300.100=sapply(1:20, function(i) mean(random.sam300.100 (100, seq(from=10, to=29, by=1)[i], 1000)))
return(samplem300.100)
}

#r2=0.1
t300s100r0.1 = trait.300.sam.100(63.5)
#r2=0.5
t300s100r0.5 = trait.300.sam.100(926)


##150 of total traits and 30 sample size.
ran.gen.150=genPositiveDefMat("c-vine", dim=150 , rangeVar=c(0.5, 0.6))

trait.150.sam.30 = function(p){
decom150.30=svd(ran.gen.150$Sigma)
recom150.30=c(p, rep(1, 149))* decom150.30$d
ran.gen.de150.30= decom150.30$u %*% diag(recom150.30) %*% t(decom150.30$v)
Sigma150.30 = ran.gen.de150.30
random_data150.30 = mvrnorm(n = 10000, rep(0, 150), Sigma150.30, empirical=FALSE)
random.sam150.30 = function(n,t,r){
as.numeric(lapply(1:r,function(i) CalcICV(cov(random_data150.30 [sample(1:10000, n), sample(1:150, t)]))))
}
samplem150.30=sapply(1:20, function(i) mean(random.sam150.30 (30, seq(from=10, to=29, by=1)[i], 1000)))
return(samplem)
}

#r2=0.1
t150s30r0.1 = trait.150.sam.30(30.6)
#r2=0.5
t150s30r0.5 = trait.150.sam.30(420)


##150 of total traits and 100 sample size.
trait.150.sam.100 = function(p){
decom150.100=svd(ran.gen.150$Sigma)
recom150.100=c(p, rep(1, 149))* decom150.100$d
ran.gen.de150.100= decom150.100$u %*% diag(recom150.100) %*% t(decom150.100$v)
Sigma150.100 = ran.gen.de150.100
random_data150.100 = mvrnorm(n = 10000, rep(0, 150), Sigma150.100, empirical=FALSE)
random.sam150.100 = function(n,t,r){
as.numeric(lapply(1:r,function(i) CalcICV(cov(random_data150.100 [sample(1:10000, n), sample(1:150, t)]))))
}
samplem150.100=sapply(1:20, function(i) mean(random.sam150.100 (100, seq(from=10, to=29, by=1)[i], 1000)))
return(samplem150.100)
}

#r2=0.1
t150s100r0.1 = trait.150.sam.100(30.6)
#r2=0.5
t150s100r0.5 = trait.150.sam.100(420)

##75 of total traits and 30 sample size.
ran.gen.75=genPositiveDefMat("c-vine", dim=75 , rangeVar=c(0.5, 0.6))

trait.75.sam.30 = function(p){
decom75.30=svd(ran.gen.75$Sigma)
recom75.30=c(p, rep(1, 74))* decom75.30$d
ran.gen.de75.30= decom75.30$u %*% diag(recom75.30) %*% t(decom75.30$v)
Sigma75.30 = ran.gen.de75.30
random_data75.30 = mvrnorm(n = 10000, rep(0, 75), Sigma75.30, empirical=FALSE)
random.sam75.30 = function(n,t,r){
as.numeric(lapply(1:r,function(i) CalcICV(cov(random_data75.30 [sample(1:10000, n), sample(1:75, t)]))))
}
samplem75.30=sapply(1:20, function(i) mean(random.sam75.30 (30, seq(from=10, to=29, by=1)[i], 1000)))
return(samplem75.30)
}

#r2=0.1
t75s30r0.1 = trait.75.sam.30(13.7)
#r2=0.5
t75s30r0.5 = trait.75.sam.30(164)


##75 of total traits and 100 sample size.
trait.75.sam.100 = function(p){
decom75.100=svd(ran.gen.75$Sigma)
recom75.100=c(p, rep(1, 74))* decom75.100$d
ran.gen.de75.100= decom75.100$u %*% diag(recom75.100) %*% t(decom75.100$v)

Sigma75.100 = ran.gen.de75.100
random_data75.100 = mvrnorm(n = 10000, rep(0, 75), Sigma75.100, empirical=FALSE)
random.sam75.100 = function(n,t,r){
as.numeric(lapply(1:r,function(i) CalcICV(cov(random_data75.100 [sample(1:10000, n), sample(1:75, t)]))))
}
samplem75.100=sapply(1:20, function(i) mean(random.sam75.100 (100, seq(from=10, to=29, by=1)[i], 1000)))
return(samplem75.100)
}

#r2=0.1
t75s100r0.1 = trait.75.sam.100(13.7)
#r2=0.5
t75s100r0.5 = trait.75.sam.100(164)



###calculation of smallest eigenvalue and skewness of log-eigenvalue distribution when r2=0.1
skewICV=function(t){
ran.gen.skew =genPositiveDefMat("c-vine", dim=t , rangeVar=c(0.5, 0.6))
  
decom.skew=svd(ran.gen.skew $Sigma)
recom.skew =c(60, rep(1, t-1))* decom.skew $d
ran.gen.de.skew = decom.skew $u %*% diag(recom.skew) %*% t(decom.skew $v)
  
ICVcal = CalcICV(ran.gen.de.skew)
  
skew = skewness(log(eigen(ran.gen.skew $Sigma)$values))
  
small = eigen(ran.gen.de.skew)$values[t]
  
total = c(ICVcal, skew, small)
  
return(total)
}
@t=number of traits


fil.test = sapply(1:1000, function(i) skewICV(50))


skewness = plot(fil.test[1,], fil.test[2,], xlab="ICV", ylab="Skewness of log-eigenvalue distribution")

eigen.small = plot(fil.test[1,], fil.test[3,], xlab="ICV", ylab="Smallest eigenvalue")


fil.data=data.frame(fil.test)

head(fil.data)

summary(lm(fil.data[1,] ~ fil.data[2,] + fil.data[3,], data=fil.data))

as.vector(fil.data[1,])
