library(AICcmodavg)
library(nls.multstart)
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)
library(MuMIn)
library(nlme)
library(plotrix)
library(lattice)
library(visreg)
library(rotl)
library(ape)
library(metafor)
library(phylosignal)
library(adephylo)
library(phylobase)

###################################################################################################
#Pre-treatment/calculation of plasticity rates, including fitting of broken stick and exponential decay
Raw.data <- read.delim("C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/Raw data.txt")

#standardize time variable to hours
Raw.data$time_h <- ifelse(is.na(Raw.data$time_d), Raw.data$time_h, Raw.data$time_d*24)

#add number of observations used for estimating rates for each experiment
n <- as.data.frame(table(Raw.data$Experiment))
names(n) <- c("Experiment", "nr.obs")
Raw.data <- merge(Raw.data,n)             

#Calculate Dt
#Use first measurement of Traitvalue as z0, use max measurement
#(for upregulation of thermal tolerance) or min measurement (for downregulation of 
#thermal tolerance) as z.inf

Exp <- vector()
z0 <- vector()
zinf <- vector()
for (i in (1:max(Raw.data$Experiment))) {
temp <- Raw.data[Raw.data$Experiment == i,]
Exp[i] <- i
z0[i] <- mean(temp$Traitvalue[temp$time_h == 0]) #for experiments with several measurements at time 0
zinf[i] <- ifelse(temp$regulation[1] == "UP", max(temp$Traitvalue), min(temp$Traitvalue))
}
temp <- as.data.frame(cbind(Exp, z0, zinf))
names(temp) <- c("Experiment", "z0", "zinf")
Raw.data <- merge(temp, Raw.data)
Raw.data$Dt <- (Raw.data$Traitvalue-Raw.data$zinf)/(Raw.data$z0-Raw.data$zinf)
#View(Data)


#fit broken stick and exponential decline models to each experiment and extract model characteristics
sigma.linear <- vector()
b <- vector()
b.SE <- vector()
b.P <- vector()
sigma.exp <- vector()
lambda <- vector()
lambda.SE <-vector()
lambda.P <- vector()
Exp <- vector() 
slope.at.end <- vector() #slope of the fitted exponential decay function at final measurement 


for (i in (1:max(Raw.data$Experiment))) {
       temp <- Raw.data[Raw.data$Experiment == i,]
#broken stick can often not be fitted if Dt = 0 for second observation (a can then take on any value between these two observations, with equally good fit)
if (temp$Dt[2] > 0)  {     mod1 <- nls_multstart(Dt ~ ifelse(time_h < a, 1-(1/a)*time_h, 0),
                     data = temp,
                     iter = 500,
                     start_lower = c(a=0.01),
                     start_upper = c(a=max(temp$time_h)*2),
                     supp_errors = 'Y',
                     convergence_count = FALSE,
                     na.action = na.omit)
sigma.linear[i] <- sigma(mod1)
b[i] <- 1/summary(mod1)$coefficients[1, 1] #this is |slope| of linear declining part, i.e. linear lambda
b.SE[i] <- summary(mod1)$coefficients[1, 2] 
b.P[i] <- summary(mod1)$coefficients[1, 4] 
}
  mod2 <- nls_multstart(Dt~exp(-lambda*time_h), data = temp, 
                        iter = 500,
                        start_lower = c(lambda=0),
                        start_upper = c(lambda=1), 
                        supp_errors = 'Y',
                        convergence_count = FALSE,
                        na.action = na.omit)
sigma.exp[i] <- sigma(mod2)
lambda[i] <- summary(mod2)$coefficients[1, 1]
lambda.SE[i] <- summary(mod2)$coefficients[1, 2]
lambda.P[i] <- summary(mod2)$coefficients[1, 4]
slope.at.end[i] <- -lambda[i]*exp(-lambda[i]*max(temp$time_h)) #derivative of function at final measurement

mod3 <- lm(Dt~1, data=temp)
Exp[i] <- temp$Experiment[1]
}

#create new data frame Model.estimates containing model parameters
Model.estimates <- as.data.frame(cbind(Exp, sigma.linear, b,b.SE,b.P, sigma.exp, lambda, lambda.SE, lambda.P, slope.at.end))
names(Model.estimates) <- c("Experiment","sigma.linear", "b", "b.SE", "b.P", "sigma.exp", "lambda", "lambda.SE","lambda.P","slope.at.end")
Model.estimates$sigma.diff <- Model.estimates$sigma.linear - Model.estimates$sigma.exp
Model.estimates$best.model <- ifelse(Model.estimates$sigma.exp<Model.estimates$sigma.linear,"Exp","Linear")

#################################################################
# Fig. S5
#################################################################
#plot data from each experiment and fitted lines for broken stick and exponential decay as estimated above 
par(mfrow = c(5,4), mai=c(0.4,0.4,0.4,0.4))
for (i in (1:max(Raw.data$Experiment))) {
  temp = Raw.data[Raw.data$Experiment ==i,]
  lambda <- Model.estimates$lambda[Model.estimates$Experiment == i]
  b <- Model.estimates$b[Model.estimates$Experiment == i]
  best <- Model.estimates$best.model[Model.estimates$Experiment == i]
  sigma <- Model.estimates$sigma.exp[Model.estimates$Experiment == i] #deltaAIC null-exp
  plot(temp$time_h,temp$Dt, xlab = "", ylab = "")  
  title(xlab = "Acclimation duration (h)", line = 2)
  title(ylab = "Acclimation response (D)", line = 2)
  text(max(temp$time_h)-max(temp$time_h)*0.2, 0.9*max(temp$Dt), labels =temp$Experiment[1])
  text(max(temp$time_h)-max(temp$time_h)*0.2, 0.75*max(temp$Dt), labels =best)
  pred.time <- seq(0,max(temp$time_h), by=max(temp$time_h)/100)
  pred.Dt <- exp(-lambda*pred.time)
  points(pred.time, pred.Dt, type = "l")
  a <- 1/b
  pred.Dt <- ifelse(pred.time < a, 1-b*pred.time, 0)
  points(pred.time, pred.Dt, type = "l", lty=2)
}


#add number of observations for each experiment to Model.estimates
for (i in (1:max(Raw.data$Experiment))) {
  temp <- Raw.data[Raw.data$Experiment == i,]
  Model.estimates$nr.obs[i] <- length(temp$Traitvalue)
}

#180 of the 290 experiments for which both models could be fitted
#best described by exponential decay
table(Model.estimates$best.model)

#but effect of number of observations
#################################################################
# Fig. S6
#################################################################
par(mfrow=c(1,1))
dev.off()
fig<-ggplot(Model.estimates, aes(x=as.factor(nr.obs), y=sigma.diff)) +
  geom_boxplot(lwd=1.5)
fig <- fig + theme_classic(base_size= 20)+xlab("Number of observations in experiment") + ylab("Residual standard error difference")
fig <- fig + geom_hline(yintercept = 0, lwd=2)
fig

#remove columns from Raw.data containing time specific values and 
#merge with Model.estimates. Resulting Final.data to be used for
#statistical analyses
str(Raw.data)
temp <- Raw.data[,c(-30, -31, -32, -37)]
#merge model outputs with raw data file
Final.data <- merge(temp, Model.estimates)
Final.data <- unique(Final.data)
#calculate a single mass variable, depending on data availability
Final.data$mass <- Final.data$Mass.g.from.paper
Final.data$mass <- ifelse(is.na(Final.data$mass), Final.data$Med.mass.g, Final.data$mass)
Final.data$mass <- ifelse(is.na(Final.data$mass), (Final.data$Min.adultmass.g + Final.data$Max.adultmass.g)/2,Final.data$mass)

#one experiment gives species name not recognised by onetree, changed to synonymouos here
Final.data$species <- ifelse(Final.data$species == "Mucrosomia caeca", "Cryptopygus caecus",Final.data$species)

#find number of species with data for each class
str(Final.data)
temp <- Final.data[,c(9, 23)]
temp <- unique(temp)
table(temp$Class)


write.table(Final.data, "C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/Final data.txt", sep=",")

###################################################################################################
###################################################################################################
#Explore estimated lambda vs. SE of estimates#
###################################################################################################
###################################################################################################

Final.data <- read.table("C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/Final data.txt", sep=",")

#################################################################
# Fig. S4
#################################################################

Final.data$variance.level <- cut(Final.data$lambda.SE,c(0,0.01, 0.02, 0.03, 100), labels=c("0-0.01", "0.01-0.02", "0.02-0.03", ">0.03"))
par(mfrow=c(1,1))
plot(as.factor(Final.data$variance.level), Final.data$lambda, xlab = "SE of estimated lambda", ylab="Estimated lambda")
text(1,2, round(with(Final.data, median(lambda[lambda.SE<0.01])), digits = 3), cex=1.1)
text(0.93,1.75,"N =", cex=1.1)
text(1.08,1.75, with(Final.data, length(lambda[lambda.SE<0.01])), cex=1.1)

text(2,2, round(with(Final.data, median(lambda[lambda.SE<0.02 & lambda.SE>=0.01 ])), digits = 3), cex=1.1)
text(1.93,1.75,"N =", cex=1.1)
text(2.08,1.75, with(Final.data, length(lambda[lambda.SE<0.02 & lambda.SE>=0.01 ])), cex=1.1)

text(3,2, round(with(Final.data, median(lambda[lambda.SE<0.03 & lambda.SE>=0.02])), digits = 3), cex=1.1)
text(2.93,1.75,"N =", cex=1.1)
text(3.08,1.75, with(Final.data, length(lambda[lambda.SE<0.03 & lambda.SE>=0.02])), cex=1.1)

text(4,2, round(with(Final.data, median(lambda[lambda.SE>=0.03])), digits = 3), cex=1.1)
text(3.93,1.75,"N =", cex=1.1)
text(4.08,1.75, with(Final.data, length(lambda[lambda.SE>=0.03])), cex=1.1)


#################################################################
# Fig. S3
#################################################################
#plot of number of experiments included vs SE threshold for
#exclusion

threshold <- seq(from=0.0001,to=10, length.out= 10000)
nr.exps <- vector()
for (i in (1:length(threshold))) {
  t <- threshold[i]
  temp <- Final.data[Final.data$lambda.SE < t,]
  nr.exps[i] <- length(temp$lambda.SE)
}
temp <- as.data.frame(cbind(threshold, nr.exps))
par(mfrow=c(2,1))
plot(temp$threshold, temp$nr.exps, type = "l", xlab = "SE of estimated lambda", ylab="Cumulative number of experiments")
plot(temp$threshold, temp$nr.exps, type = "l", xlim = c(0,0.1), xlab = "SE of estimated lambda", ylab="Cumulative number of experiments")
abline(v=0.02, lty="dashed")
abline(v=0.01, lty="dashed")

length(Final.data$Experiment)-1
temp <- Final.data[Final.data$lambda.SE<0.01,]
length(temp$Experiment)-1
temp <- Final.data[Final.data$lambda.SE<0.02,]
length(temp$Experiment)-1
###################################################################################################
###################################################################################################
#check if model with phylogeny fits better than one without
#Only include experiments where body mass data are available
#and lambda.SE < 0.01#
###################################################################################################
###################################################################################################

Final.data <- read.table("C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/Final data.txt", sep=",")
#remove data without body mass data
Final.data <- Final.data[!is.na(Final.data$mass),]
Final.data <- Final.data[Final.data$lambda.SE<0.01,]
#one negative estimate, remove this
Final.data <- Final.data[Final.data$lambda>0,]
#six experiments lacked information on acclimation temperature, exclude these
Final.data <- Final.data[!is.na(Final.data$acctemp2),]
#Thecostraca and Turbellaria only 1 and 2 species, respectively. Remove these
Final.data <- Final.data[Final.data$Class != "Thecostraca" & Final.data$Class != "Turbellaria",] 



#build phylogeny
species.list <- Final.data$species
species.list <- unique(species.list)
taxa <- tnrs_match_names(names = species.list)
taxa$search_string <-   gsub("(^[[:alpha:]])", "\\\\U\\\\1", taxa$search_string, perl=TRUE)
taxon_map <- structure(taxa$search_string, names = taxa$unique_name)
tree <- tol_induced_subtree(ott_ids = ott_id(taxa)[is_in_tree(ott_id(taxa))])
otl_tips <- strip_ott_ids(tree$tip.label, remove_underscores = TRUE)
tree$tip.label <- taxon_map[ otl_tips ]
tree$node.label <- NULL
Final.data <- Final.data[Final.data$species %in% tree$tip.label, ]

# converting the non-ultrametric tree to ultrametric tree 
tree_sub <- compute.brlen(tree, power = 1)

# plotting the ultrametric tree
plot(tree_sub, cex = .8, label.offset = .1, no.margin = TRUE)

# computing the variance-covariance matrix from the phylogeny tree
A <- vcv(tree_sub)

# creating a variable to distinguish the phylogenetic component from the 
# non-phlylogenetic component
Final.data$phylo <- Final.data$species

#calculate sampling variance from standard errors
Final.data$sampvar <- Final.data$lambda.SE^2

Final.data$id <- 1:nrow(Final.data)

# fitting models with different random structure
mod1 <- rma.mv(lambda, sampvar,mods = ~ Class +mass+slope.at.end + acctemp2+measurement.type.simplified, random = list(~ 1 | study, ~ 1 | species, ~ 1 | phylo, ~1|id),
               R = list(phylo = A), data = Final.data, sparse = TRUE, method = "REML")

mod2 <- rma.mv(lambda, sampvar,mods = ~ Class +mass+slope.at.end + acctemp2+measurement.type.simplified, random = list(~ 1 | study, ~ 1 | species, ~1|id),
               data = Final.data, sparse = TRUE, method = "REML")

AICc(mod1,mod2)

###################################################################################################
###################################################################################################
#proceed to find best fixed structure in model without phylogeny
###################################################################################################
###################################################################################################

###################################################################################################
#First with only experiments where body mass data are available
#, including only experiments with lambda.SE < 0.01#
###################################################################################################
Final.data <- read.table("C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/Final data.txt", sep=",")
#remove data without body mass data
Final.data <- Final.data[!is.na(Final.data$mass),]
Final.data <- Final.data[Final.data$lambda.SE<0.01,]
#one negative estimate, remove this
Final.data <- Final.data[Final.data$lambda>0,]
#six experiments lacked information on acclimation temperature, exclude these
Final.data <- Final.data[!is.na(Final.data$acctemp2),]
#Thecostraca and Turbellaria only 1 and 2 species, respectively. Remove these
Final.data <- Final.data[Final.data$Class != "Thecostraca" & Final.data$Class != "Turbellaria",] 

#descriptives of this data set
length(Final.data$Experiment)
str(Final.data)
temp <- Final.data[,c(9,23)]
temp <- unique(temp)
length(temp$species)
str(temp)
table(temp$Class)
temp <- Final.data[,c(9,23,39)]
means <- as.data.frame(as.table(tapply(temp$lambda, temp$species, mean)))
mean(means$Freq)
sd(means$Freq)
min(means$Freq)
max(means$Freq)

#calculate sampling variance from standard errors
Final.data$sampvar <- Final.data$lambda.SE^2

#add observation level random effect variable
Final.data$id <- 1:nrow(Final.data)

#to make MuMIn work with metafor models
eval(metafor:::.MuMIn)

# fitting full model
mod <- rma.mv(lambda, sampvar,mods = ~ Class +mass+slope.at.end + acctemp2+measurement.type.simplified, random = list(~ 1 | study, ~ 1 | species, ~1|id),
                data = Final.data, sparse = TRUE, method = "ML")

#check residual distributions
res <- residuals(mod)
hist(res)
Final.data$Class <- as.factor(Final.data$Class)
plot(Final.data$Class,res)
plot(Final.data$acctemp2,res)
plot(Final.data$slope.at.end,res) #some pattern here
#but results remain similar when removing 
#slope at end < -0.002 and <-0.001 (supplement)

##################for Table 1
#compare model fits
res <- dredge(mod, trace=2)
subset(res, delta <= 20, recalc.weights = FALSE)

# #multimodel inference
# summary(model.avg(res))
# 
# #relative importance values for the predictors 
# sw(res)

#refit best model with REML to get parameter estimates
mod <- rma.mv(lambda, sampvar,mods = ~ Class +    slope.at.end + acctemp2, random = list(~ 1 | study, ~ 1 | species, ~1|id),
               data = Final.data, sparse = TRUE, method = "REML")
summary(mod)



###################################################################################################
#Previous analysis show no effect of body mass. Re-analyse including experiments without body size data
#, including only experiments with lambda.SE < 0.01#
###################################################################################################
Final.data <- read.table("C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/Final data.txt", sep=",")
Final.data <- Final.data[Final.data$lambda.SE<0.01,]
#one negative estimate, remove this
Final.data <- Final.data[Final.data$lambda>0,]
#six experiments lacked information on acclimation temperature, exclude these
Final.data <- Final.data[!is.na(Final.data$acctemp2),]
#Thecostraca and Turbellaria only 1 and 2 species, respectively. Remove these
Final.data <- Final.data[Final.data$Class != "Thecostraca" & Final.data$Class != "Turbellaria",] 

#descriptives of this data set
length(Final.data$Experiment)
str(Final.data)
temp <- Final.data[,c(9,23)]
temp <- unique(temp)
length(temp$species)
table(temp$Class)
temp <- Final.data[,c(9,23,39)]
means <- as.data.frame(as.table(tapply(temp$lambda, temp$species, mean)))
mean(means$Freq)
sd(means$Freq)
min(means$Freq)
max(means$Freq)


#calculate sampling variance from standard errors
Final.data$sampvar <- Final.data$lambda.SE^2

#add observation level random effect variable
Final.data$id <- 1:nrow(Final.data)

#to make MuMIn work with metafor models
eval(metafor:::.MuMIn)

# fitting full model
mod <- rma.mv(lambda, sampvar,mods = ~ Class +    slope.at.end + acctemp2+measurement.type.simplified, random = list(~ 1 | study, ~ 1 | species, ~1|id),
               data = Final.data, sparse = TRUE, method = "ML")

#check residual distributions
res <- residuals(mod)
hist(res)
Final.data$Class <- as.factor(Final.data$Class)
plot(Final.data$Class,res)
plot(Final.data$acctemp2,res)
plot(Final.data$slope.at.end,res) #some pattern here
#but results remain similar when removing 
#slope at end < -0.002 and <-0.001 (supplement)

##################for Table 1
#compare model fits
res <- dredge(mod, trace=2)
subset(res, delta <= 20, recalc.weights = FALSE)

# #multimodel inference
# summary(model.avg(res))
# 
# #relative importance values for the predictors 
# sw(res)

##################for Table 2
#refit best model with REML to get parameter estimates
mod <- rma.mv(lambda, sampvar,mods = ~ Class-1 +    slope.at.end + acctemp2, random = list(~ 1 | study, ~ 1 | species, ~1|id),
              data = Final.data, sparse = TRUE, method = "REML")
summary(mod)
#mod$I2

# calculate I2, based on https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
W <- diag(1/mod$vi)
X <- model.matrix(mod)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2 <- 100 * sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))

#how much of the total variance can be attributed to different levels (study, species, id)
I2s <- 100 * mod$sigma2 / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))

#pseudo R2
mod1<-rma.mv(lambda, sampvar,mods = ~ Class-1 +    slope.at.end + acctemp2, random = list(~ 1 | study, ~ 1 | species, ~1|id),
             data = Final.data, sparse = TRUE, method = "REML")
mod2<-rma.mv(lambda, sampvar,random = list(~ 1 | study, ~ 1 | species, ~1|id),
             data = Final.data, sparse = TRUE, method = "REML")

###pseudo-R2
(sum(mod2$sigma2) - sum(mod1$sigma2)) / sum(mod2$sigma2)

###################################################################################################
#repeat above but including only experiments with lambda.SE < 0.02#
###################################################################################################
Final.data <- read.table("C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/Final data.txt", sep=",")
Final.data <- Final.data[Final.data$lambda.SE<0.02,]
#one negative estimate, remove this
Final.data <- Final.data[Final.data$lambda>0,]
#six experiments lacked information on acclimation temperature, exclude these
Final.data <- Final.data[!is.na(Final.data$acctemp2),]
#Thecostraca and Turbellaria only 1 and 2 species, respectively. Remove these
Final.data <- Final.data[Final.data$Class != "Thecostraca" & Final.data$Class != "Turbellaria",] 

#descriptives of this data set
length(Final.data$Experiment)
str(Final.data)
temp <- Final.data[,c(9,23)]
temp <- unique(temp)
length(temp$species)
table(temp$Class)
temp <- Final.data[,c(9,23,39)]
means <- as.data.frame(as.table(tapply(temp$lambda, temp$species, mean)))
mean(means$Freq)
sd(means$Freq)
min(means$Freq)
max(means$Freq)


#calculate sampling variance from standard errors
Final.data$sampvar <- Final.data$lambda.SE^2

#add observation level random effect variable
Final.data$id <- 1:nrow(Final.data)

#to make MuMIn work with metafor models
eval(metafor:::.MuMIn)

# fitting full model
mod <- rma.mv(lambda, sampvar,mods = ~ Class +    slope.at.end + acctemp2+measurement.type.simplified, random = list(~ 1 | study, ~ 1 | species, ~1|id),
              data = Final.data, sparse = TRUE, method = "ML")

#check residual distributions
res <- residuals(mod)
hist(res)
Final.data$Class <- as.factor(Final.data$Class)
plot(Final.data$Class,res)
plot(Final.data$acctemp2,res)
plot(Final.data$slope.at.end,res) #some pattern here
#but results remain similar when removing 
#slope at end < -0.002 and <-0.001 (supplement)



##################for Table S1
#compare model fits
res <- dredge(mod, trace=2)
subset(res, delta <= 20, recalc.weights = FALSE)

# #multimodel inference
# summary(model.avg(res))
# 
# #relative importance values for the predictors 
# sw(res)


#refit best model with REML to get parameter estimates
mod <- rma.mv(lambda, sampvar,mods = ~ Class +    slope.at.end + acctemp2, random = list(~ 1 | study, ~ 1 | species, ~1|id),
              data = Final.data, sparse = TRUE, method = "REML")
summary(mod)



###################################################################################################
# Re-analyse including experiments without body size data
#, including only experiments with lambda.SE < 0.01 & slope at final measuremen > -0.002#
###################################################################################################
Final.data <- read.table("C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/Final data.txt", sep=",")
Final.data <- Final.data[Final.data$lambda.SE<0.01,]
#one negative estimate, remove this
Final.data <- Final.data[Final.data$lambda>0,]
#six experiments lacked information on acclimation temperature, exclude these
Final.data <- Final.data[!is.na(Final.data$acctemp2),]
#Thecostraca and Turbellaria only 1 and 2 species, respectively. Remove these
Final.data <- Final.data[Final.data$Class != "Thecostraca" & Final.data$Class != "Turbellaria",] 
#final slope criterion
Final.data <- Final.data[Final.data$slope.at.end > -0.002,] 



#descriptives of this data set
length(Final.data$Experiment)
str(Final.data)
temp <- Final.data[,c(9,23)]
temp <- unique(temp)
length(temp$species)
table(temp$Class)
temp <- Final.data[,c(9,23,39)]
means <- as.data.frame(as.table(tapply(temp$lambda, temp$species, mean)))
mean(means$Freq)
sd(means$Freq)
min(means$Freq)
max(means$Freq)


#calculate sampling variance from standard errors
Final.data$sampvar <- Final.data$lambda.SE^2

#add observation level random effect variable
Final.data$id <- 1:nrow(Final.data)

#to make MuMIn work with metafor models
eval(metafor:::.MuMIn)

# fitting full model
mod <- rma.mv(lambda, sampvar,mods = ~ Class +    slope.at.end + acctemp2+measurement.type.simplified, random = list(~ 1 | study, ~ 1 | species, ~1|id),
              data = Final.data, sparse = TRUE, method = "ML")

#check residual distributions
res <- residuals(mod)
hist(res)
Final.data$Class <- as.factor(Final.data$Class)
plot(Final.data$Class,res)
plot(Final.data$acctemp2,res)
plot(Final.data$slope.at.end,res)

##################for Table S3
#compare model fits
res <- dredge(mod, trace=2)
subset(res, delta <= 20, recalc.weights = FALSE)

# #multimodel inference
# summary(model.avg(res))
# 
# #relative importance values for the predictors 
# sw(res)

##################for Table S5
#refit best model with REML to get parameter estimates
mod <- rma.mv(lambda, sampvar,mods = ~ Class-1 +    slope.at.end + acctemp2+measurement.type.simplified-1, random = list(~ 1 | study, ~ 1 | species, ~1|id),
              data = Final.data, sparse = TRUE, method = "REML")
summary(mod)


# calculate I2, based on https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
W <- diag(1/mod$vi)
X <- model.matrix(mod)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2 <- 100 * sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))

#how much of the total variance can be attributed to different levels (study, species, id)
I2s <- 100 * mod$sigma2 / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))

#pseudo R2
mod1<-rma.mv(lambda, sampvar,mods = ~ Class-1 +    slope.at.end + acctemp2+measurement.type.simplified-1, random = list(~ 1 | study, ~ 1 | species, ~1|id),
             data = Final.data, sparse = TRUE, method = "REML")
mod2<-rma.mv(lambda, sampvar,random = list(~ 1 | study, ~ 1 | species, ~1|id),
             data = Final.data, sparse = TRUE, method = "REML")

###pseudo-R2
(sum(mod2$sigma2) - sum(mod1$sigma2)) / sum(mod2$sigma2)


###################################################################################################
# Re-analyse including experiments without body size data
#, including only experiments with lambda.SE < 0.01 & slope at final measuremen > -0.001#
###################################################################################################
Final.data <- read.table("C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/Final data.txt", sep=",")
Final.data <- Final.data[Final.data$lambda.SE<0.01,]
#one negative estimate, remove this
Final.data <- Final.data[Final.data$lambda>0,]
#six experiments lacked information on acclimation temperature, exclude these
Final.data <- Final.data[!is.na(Final.data$acctemp2),]
#Thecostraca and Turbellaria only 1 and 2 species, respectively. Remove these
Final.data <- Final.data[Final.data$Class != "Thecostraca" & Final.data$Class != "Turbellaria",] 
#final slope criterion
Final.data <- Final.data[Final.data$slope.at.end > -0.001,] 



#descriptives of this data set
length(Final.data$Experiment)
str(Final.data)
temp <- Final.data[,c(9,23)]
temp <- unique(temp)
length(temp$species)
table(temp$Class)
temp <- Final.data[,c(9,23,39)]
means <- as.data.frame(as.table(tapply(temp$lambda, temp$species, mean)))
mean(means$Freq)
sd(means$Freq)
min(means$Freq)
max(means$Freq)


#calculate sampling variance from standard errors
Final.data$sampvar <- Final.data$lambda.SE^2

#add observation level random effect variable
Final.data$id <- 1:nrow(Final.data)

#to make MuMIn work with metafor models
eval(metafor:::.MuMIn)

# fitting full model
mod <- rma.mv(lambda, sampvar,mods = ~ Class +    slope.at.end + acctemp2+measurement.type.simplified, random = list(~ 1 | study, ~ 1 | species, ~1|id),
              data = Final.data, sparse = TRUE, method = "ML")

#check residual distributions
res <- residuals(mod)
hist(res)
Final.data$Class <- as.factor(Final.data$Class)
plot(Final.data$Class,res)
plot(Final.data$acctemp2,res)
plot(Final.data$slope.at.end,res)


##################for Table S4
#compare model fits
res <- dredge(mod, trace=2)
subset(res, delta <= 20, recalc.weights = FALSE)

# #multimodel inference
# summary(model.avg(res))
# 
# #relative importance values for the predictors 
# sw(res)

##################for Table S6
#refit best model with REML to get parameter estimates
mod <- rma.mv(lambda, sampvar,mods = ~ Class-1 +    slope.at.end + acctemp2, random = list(~ 1 | study, ~ 1 | species, ~1|id),
              data = Final.data, sparse = TRUE, method = "REML")
summary(mod)



# calculate I2, based on https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
W <- diag(1/mod$vi)
X <- model.matrix(mod)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2 <- 100 * sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))

#how much of the total variance can be attributed to different levels (study, species, id)
I2s <- 100 * mod$sigma2 / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))

#pseudo R2
mod1<-rma.mv(lambda, sampvar,mods = ~ Class-1 +    slope.at.end + acctemp2, random = list(~ 1 | study, ~ 1 | species, ~1|id),
             data = Final.data, sparse = TRUE, method = "REML")
mod2<-rma.mv(lambda, sampvar,random = list(~ 1 | study, ~ 1 | species, ~1|id),
             data = Final.data, sparse = TRUE, method = "REML")

###pseudo-R2
(sum(mod2$sigma2) - sum(mod1$sigma2)) / sum(mod2$sigma2)




##########################
#create output for table S2
Final.data <- read.table("C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/Final data.txt", sep=",")
str(Final.data)
temp <- Final.data[,c(1,6,7,9,11,12,23,39,40,45)]
temp$mass <- ifelse(is.na(temp$mass) ==T, "No", "Yes")
temp <- unique(temp)
write.table(temp, "C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/for supplementary table.txt", sep=",")
str(temp)

#############################





###################
# Figure 2
###################

Final.data <- read.table("C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/Final data.txt", sep=",")
Final.data <- Final.data[Final.data$lambda.SE<0.01,]
#one negative estimate, remove this
Final.data <- Final.data[Final.data$lambda>0,]
#six experiments lacked information on acclimation temperature, exclude these
Final.data <- Final.data[!is.na(Final.data$acctemp2),]
#Thecostraca and Turbellaria only 1 and 2 species, respectively. Remove these
Final.data <- Final.data[Final.data$Class != "Thecostraca" & Final.data$Class != "Turbellaria",] 

Final.data$Class[Final.data$Class == "Amphibia"] <-"Amphibians"
Final.data$Class[Final.data$Class == "Reptilia"] <-"Reptiles"
Final.data$Class[Final.data$Class == "Insecta"] <-"Insects"
Final.data$Class[Final.data$Class == "Malacostraca"] <-"Crustaceans"
Final.data$Class[Final.data$Class == "Osteichtyes"] <-"Fishes"

#calculate mean lambda for each species
temp <- as.data.frame(as.table(tapply(Final.data$lambda,Final.data$species, mean)))
names(temp) <- c("species", "mean.lambda")

#merge with class information
temp2 <- Final.data[,c(9,23)]
temp2 <- unique(temp2)
temp3 <- merge(temp, temp2)

par(mfrow=c(1,1))
temp3$Class <- factor(temp3$Class, levels = c("Amphibians", "Reptiles", "Insects", "Crustaceans", "Fishes"))

p1<-ggplot(temp3, aes(x=Class, y=mean.lambda, fill=Class)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"))+
  geom_boxplot(lwd=1.5)
p1 <- p1 + theme_classic(base_size= 20)+ ylab("Rate of acclimation (lamda)")
p1 <- p1 + theme(legend.position = "none")
p1 <- p1 + theme(axis.title.x = element_blank())
p1

###################
# Figure 3
###################
Final.data <- read.table("C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/Final data.txt", sep=",")
Final.data <- Final.data[Final.data$lambda.SE<0.01,]
#one negative estimate, remove this
Final.data <- Final.data[Final.data$lambda>0,]
#six experiments lacked information on acclimation temperature, exclude these
Final.data <- Final.data[!is.na(Final.data$acctemp2),]
#Thecostraca and Turbellaria only 1 and 2 species, respectively. Remove these
Final.data <- Final.data[Final.data$Class != "Thecostraca" & Final.data$Class != "Turbellaria",] 
Final.data$Class[Final.data$Class == "Amphibia"] <-"Amphibians"
Final.data$Class[Final.data$Class == "Reptilia"] <-"Reptiles"
Final.data$Class[Final.data$Class == "Insecta"] <-"Insects"
Final.data$Class[Final.data$Class == "Malacostraca"] <-"Crustaceans"
Final.data$Class[Final.data$Class == "Osteichtyes"] <-"Fishes"
#calculate sampling variance from standard errors
Final.data$sampvar <- Final.data$lambda.SE^2

Final.data$id <- 1:nrow(Final.data)

mod1 <- rma.mv(lambda, sampvar,mods = ~ Class-1 +    slope.at.end + acctemp2, random = list(~ 1 | study, ~ 1 | species, ~1|id),
               data = Final.data, sparse = TRUE, method = "REML")

amph <- Final.data[Final.data$Class== "Amphibians",]
crus <- Final.data[Final.data$Class== "Crustaceans",]
fish <- Final.data[Final.data$Class== "Fishes",]
inse <- Final.data[Final.data$Class== "Insects",]
rept <- Final.data[Final.data$Class== "Reptiles",]

p <- ggplot(Final.data, aes(x = acctemp2, y = lambda, color = Class) ) +
  scale_color_manual(values=c("#E69F00","#F0E442","#0072B2","#009E73",   "#56B4E9"))+
  geom_point(size = 3) + 
geom_segment(aes(x = min(amph$acctemp2), xend = max(amph$acctemp2), y=mod1$b[1] + 0.0006*min(amph$acctemp2),yend=mod1$b[1] + 0.0006*max(amph$acctemp2)), colour="#E69F00",size=2)+
geom_segment(aes(x = min(crus$acctemp2), xend = max(crus$acctemp2), y=mod1$b[2] + 0.0006*min(crus$acctemp2),yend=mod1$b[2] + 0.0006*max(crus$acctemp2)), colour="#F0E442",size=2)+
geom_segment(aes(x = min(fish$acctemp2), xend = max(fish$acctemp2), y=mod1$b[3] + 0.0006*min(fish$acctemp2),yend=mod1$b[3] + 0.0006*max(fish$acctemp2)), colour="#0072B2",size=2)+
geom_segment(aes(x = min(inse$acctemp2), xend = max(inse$acctemp2), y=mod1$b[4] + 0.0006*min(inse$acctemp2),yend=mod1$b[4] + 0.0006*max(inse$acctemp2)), colour="#009E73",size=2)+
geom_segment(aes(x = min(rept$acctemp2), xend = max(rept$acctemp2), y=mod1$b[5] + 0.0006*min(rept$acctemp2),yend=mod1$b[5] + 0.0006*max(rept$acctemp2)), colour="#56B4E9",size=2)
  p <- p + theme_classic(base_size= 20)+xlab("Acclimation temperature") + ylab("Rate of acclimation (lamda)")
p


###################
# Figure S7
###################
Final.data <- read.table("C:/Daglig sync/work done in Texas/rate of thermal acclimation/resubmission ecology letters/Final data.txt", sep=",")
Final.data <- Final.data[Final.data$lambda.SE<0.01,]
Final.data <- Final.data[Final.data$Class != "Thecostraca" & Final.data$Class != "Turbellaria",] #only 1-2 species from these classes
Final.data$Class[Final.data$Class == "Amphibia"] <-"Amphibians"
Final.data$Class[Final.data$Class == "Reptilia"] <-"Reptiles"
Final.data$Class[Final.data$Class == "Insecta"] <-"Insects"
Final.data$Class[Final.data$Class == "Malacostraca"] <-"Crustaceans"
Final.data$Class[Final.data$Class == "Osteichtyes"] <-"Fishes"

xyplot(Final.data$lambda~Final.data$slope.at.end|Final.data$Class, xlab = "Slope at end", ylab="Lambda")

#############################


###################
# Figure S2
###################
par(mfrow=c(2,2))
lambda <- 0.1
time <- 0:100
trait <- 37 + 3*exp(-lambda*time)
plot(time, trait, xlab = "Time", ylab = "CTmax", xlim = c(0,100),ylim = c(37,40))
text(98,40,"A")
time2 <- 0:30
trait.observed <- 37 + 3*exp(-lambda*time2)
plot(time2, trait.observed, xlab = "Time", ylab = "CTmax", xlim = c(0,30),ylim = c(37,40))
text(29,40,"C")

time.obs2 <- 0:100
trait.obs2 <- 37 + 3*exp(-lambda*time.obs2)
Dt2 <- (trait.obs2-min(trait.obs2))/(max(trait.obs2)-min(trait.obs2))
plot(time.obs2, Dt2, xlab = "Time", ylab = "D", xlim = c(0,100))
text(98,1,"B")
temp2 <- as.data.frame(cbind(time.obs2, Dt2))
mod <- nls_multstart(Dt2~exp(-lambda*time.obs2), data = temp2, 
                     iter = 250,
                     start_lower = c(lambda=0),
                     start_upper = c(lambda=1),
                     supp_errors = 'Y',
                     convergence_count = 100,
                     na.action = na.omit,
                     lower = c(lambda=0))
est.lambda2 <- summary(mod)$coefficients[1, 1]

time.obs1 <- 0:30
trait.obs1 <- 37 + 3*exp(-lambda*time.obs1)
Dt1 <- (trait.obs1-min(trait.obs1))/(max(trait.obs1)-min(trait.obs1))
plot(time.obs1, Dt1, xlab = "Time", ylab = "D", xlim = c(0,30))
text(29,1,"D")
temp1 <- as.data.frame(cbind(time.obs1, Dt1))
mod <- nls_multstart(Dt1~exp(-lambda*time.obs1), data = temp1, 
                     iter = 250,
                     start_lower = c(lambda=0),
                     start_upper = c(lambda=1),
                     supp_errors = 'Y',
                     convergence_count = 100,
                     na.action = na.omit,
                     lower = c(lambda=0))
est.lambda1 <- summary(mod)$coefficients[1, 1]


