library(psych)
library(FactoMineR)
library(lme4)
library(car)
library(effects)
library(vegan)
library(MCMCglmm)
library(export)
library(ape)

########
# Data #
########
data <- read.table("CHC data.txt", header=T)
data <- na.omit(data)

pedigree <- read.table("pedigree.txt", header=T)

n = length(data[,6:length(data)])

perc.data <- data
perc.data[,6:length(data)]  <- (data[,6:length(data)]/apply(data[,6:length(data)], 1, sum))*100

lr.data <- data
lr.data[,6:length(data)] <- log(data[,6:length(data)]/apply(data[,6:length(data)], 1, geometric.mean))


#########################
# Mean profiles and PCA #
#########################
input <- lr.data

mean.profiles <- data.frame(matrix(NA, nrow=n, ncol=3))
row.names(mean.profiles) <- colnames(input)[6:length(input)]
colnames(mean.profiles) <- c("female","male","nymph")
for (i in 1:n){
  mean.profiles[i,] <- round(tapply(input[,i+5], input$gender, mean, na.rm=T),2)
}

pca <- prcomp(input[,6:length(input)], scale.=F, tol=0.5)
pca$sdev^2
summary(pca)
plot(pca$x[,2]~pca$x[,1], pch=21, bg=input$gender, xlab="PC1 (51.4 %)", ylab="PC2 (19.9 %)")
abline(h=0, lty=2)
abline(v=0, lty=2)

write.csv(cbind(data, pca$x), file="pca-individuals.csv")
write.csv(pca$rotation, file="pca-variables.csv")


######################
# Pairwise distances #
######################
relation <- read.table("Pairwise relationships.txt", header=T, row.names=1)
relatedness <- read.table("Rlat SS0625.txt", header=T, row.names=1)

nymph.data <- subset(data, age=="nymph")
str(nymph.data)
nymph.perc.data <- nymph.data
nymph.perc.data[,6:length(nymph.perc.data)] = (nymph.perc.data[,6:length(nymph.perc.data)]/apply(nymph.perc.data[,6:length(nymph.perc.data)], 1, sum))*100

nymph.lr.data <- nymph.data
nymph.lr.data[,6:length(nymph.lr.data)] = log(nymph.lr.data[,6:length(nymph.lr.data)]/apply(nymph.lr.data[,6:length(nymph.lr.data)], 1, geometric.mean))

input <- nymph.lr.data

## all compounds together
dist <- as.matrix(dist(input[,6:length(input)], diag=T, upper=T))
results <- tapply(dist, as.matrix(relation), mean)
results.sd <- tapply(dist, as.matrix(relation), sd)
write.csv(rbind(results, results.sd), file="pair_mean_all_compounds.csv")
mantel.test(dist, relatedness, nperm=1000, graph=T)

# 5 most heritable compounds are: 4, 8, 15, 19, 22
names(input[,c(4, 8, 15, 19, 22)+5])
dist <- as.matrix(dist(input[,c(4, 8, 15, 19, 22)+5], diag=T, upper=T))
results <- tapply(dist, as.matrix(relation), mean)
results.sd <- tapply(dist, as.matrix(relation), sd)
write.csv(rbind(results, results.sd), file="pair_mean_5_most_heritable.csv")
mantel.test(dist, relatedness, nperm=1000, graph=T)

# 5  least heritable compounds are: 5, 10, 17, 21, 25
names(input[,c(5, 10, 17, 21, 25)+5])
dist <- as.matrix(dist(input[,c(5, 10, 17, 21, 25)+5], diag=T, upper=T))
results <- tapply(dist, as.matrix(relation), mean)
results.sd <- tapply(dist, as.matrix(relation), sd)
write.csv(rbind(results, results.sd), file="pair_mean_5_least_heritable.csv")
mantel.test(dist, relatedness, nperm=1000, graph=T)


##########################
# Heritability estimates #
##########################
input <- lr.data
pca <- prcomp(input[,6:length(input)], scale.=F, tol=0.5)
input <- cbind(input, pca$x)

n <- length(input[,6:length(input)])

df <- data.frame(matrix(NA, ncol=3, nrow=n))
row.names(df) <- colnames(input)[6:length(input)]
colnames(df) <- c("estimate", "lower", "upper")

posterior.heritability <-  do.call("list", replicate(n, df, simplify=FALSE))
names(posterior.heritability) <- row.names(df)

genetic.correlation <-  do.call("list", replicate(n, df, simplify=FALSE))
names(genetic.correlation) <- row.names(df)

maternal.correlation <-  do.call("list", replicate(n, df, simplify=FALSE))
names(maternal.correlation) <- row.names(df)

for (i in 1:n){
  for (j in 1:n){
    temp.data <- input[,c(1:5,i+5,j+5)]
    colnames(temp.data)[6:7] <- c("CHC.1", "CHC.2")
    
    phen.var <- matrix(0, nrow=2, ncol=2)
    diag(phen.var) <- diag(var(cbind(temp.data$CHC.1, temp.data$CHC.2)))
    
    prior <- list(G=list(G1=list(V=phen.var/4,n=2),
                         G2=list(V=phen.var/4,n=2),
                         G3=list(V=phen.var/4,n=2)),
                  R=list(V=phen.var/4,n=2))
    
    model.bi <- MCMCglmm(cbind(CHC.1, CHC.2) ~ trait-1 + trait:age,
                         random= ~ us(trait):animal + us(trait):mother + us(trait):strain,
                         rcov= ~ us(trait):units,
                         family=c("gaussian","gaussian"),
                         pedigree=pedigree, data=temp.data,
                         nitt=55000, thin=50, burnin=5000,
                         prior=prior, verbose=T)
    
    posterior.heritability.CHC.1 <- model.bi$VCV[,"traitCHC.1:traitCHC.1.animal"]/(model.bi$VCV[,"traitCHC.1:traitCHC.1.animal"]+model.bi$VCV[,"traitCHC.1:traitCHC.1.strain"]+model.bi$VCV[,"traitCHC.1:traitCHC.1.mother"]+model.bi$VCV[,"traitCHC.1:traitCHC.1.units"])
    posterior.heritability[[i]][j,1] <- round(posterior.mode(posterior.heritability.CHC.1),3)
    posterior.heritability[[i]][j,c(2:3)] <- round(HPDinterval(posterior.heritability.CHC.1, 0.95),3)

    posterior.heritability.CHC.2 <- model.bi$VCV[,"traitCHC.2:traitCHC.2.animal"]/(model.bi$VCV[,"traitCHC.2:traitCHC.2.animal"]+model.bi$VCV[,"traitCHC.2:traitCHC.2.strain"]+model.bi$VCV[,"traitCHC.2:traitCHC.2.mother"]+model.bi$VCV[,"traitCHC.2:traitCHC.2.units"])
    posterior.heritability[[j]][i,1] <- round(posterior.mode(posterior.heritability.CHC.2),3)
    posterior.heritability[[j]][i,c(2:3)] <- round(HPDinterval(posterior.heritability.CHC.2, 0.95),3)

    genetic.correlation.CHC.1.CHC.2 <- model.bi$VCV[,"traitCHC.1:traitCHC.2.animal"]/sqrt(model.bi$VCV[,"traitCHC.1:traitCHC.1.animal"]*model.bi$VCV[,"traitCHC.2:traitCHC.2.animal"])
    genetic.correlation[[i]][j,1] <- round(posterior.mode(genetic.correlation.CHC.1.CHC.2),3)
    genetic.correlation[[i]][j,c(2:3)] <- round(HPDinterval(genetic.correlation.CHC.1.CHC.2, 0.95),3)
    genetic.correlation[[j]][i,1] <- round(posterior.mode(genetic.correlation.CHC.1.CHC.2),3)
    genetic.correlation[[j]][i,c(2:3)] <- round(HPDinterval(genetic.correlation.CHC.1.CHC.2, 0.95),3)

    maternal.correlation.CHC.1.CHC.2 <- model.bi$VCV[,"traitCHC.1:traitCHC.2.mother"]/sqrt(model.bi$VCV[,"traitCHC.1:traitCHC.1.mother"]*model.bi$VCV[,"traitCHC.2:traitCHC.2.mother"])
    maternal.correlation[[i]][j,1] <- round(posterior.mode(maternal.correlation.CHC.1.CHC.2),3)
    maternal.correlation[[i]][j,c(2:3)] <- round(HPDinterval(maternal.correlation.CHC.1.CHC.2, 0.95),3)
    maternal.correlation[[j]][i,1] <- round(posterior.mode(maternal.correlation.CHC.1.CHC.2),3)
    maternal.correlation[[j]][i,c(2:3)] <- round(HPDinterval(maternal.correlation.CHC.1.CHC.2, 0.95),3)
}}

posterior.heritability.mean <- df
for (p in 1:n){
  posterior.heritability.mean[p,] <- apply(posterior.heritability[[p]],2,mean)
}

posterior.heritability.mean <- posterior.heritability.mean
genetic.correlation <- genetic.correlation
maternal.correlation <- maternal.correlation

write.csv(posterior.heritability.mean, file="posterior heritability lr-transformed unscaled.csv")
write.csv(genetic.correlation, file="genetic correlation lr-transformed unscaled.csv")
write.csv(maternal.correlation, file="maternal correlation lr-transformed unscaled.csv")

genetic.correlation.matrix <- matrix(NA, n, n)
for (i in 1:n){
  for (j in 1:n){
    genetic.correlation.matrix[i,j] <- genetic.correlation[[i]][j,1]
  }
}

maternal.correlation.matrix <- matrix(NA, n, n)
for (i in 1:n){
  for (j in 1:n){
    maternal.correlation.matrix[i,j] <- maternal.correlation[[i]][j,1]
  }
}

figure.data <- data.frame(matrix(NA, ncol=n, nrow=n))
colnames(figure.data) <- colnames(input[6:length(input)])
row.names(figure.data) <- colnames(input[6:length(input)])

figure.data[lower.tri(figure.data)] <- genetic.correlation.matrix[lower.tri(genetic.correlation.matrix)]
figure.data[upper.tri(figure.data)] <- maternal.correlation.matrix[upper.tri(maternal.correlation.matrix)]
diag(figure.data) <- posterior.heritability.mean[,1]

write.csv(figure.data, file="figure data lr-transformed unscaled.csv")


# Chain length effect #
chain.length <- c(27, 27, 27, 27, 27, 27, 28, 27, 28, 28, 28, 29, 29, 29, 29, 29, 29, 29, 29, 30, 30, 31, 31, 31, 32)
cl.data <- cbind(chain.length, genetic.correlation.matrix[1:25,])
#cl.data <- cbind(chain.length, maternal.correlation.matrix[1:25,]) 

fit <- lm(n.C27.estimate ~ chain.length , data=cl.data)
summary(fit)$coef[2,c(1,4)]

fit <- lm(nC29.estimate ~ chain.length , data=cl.data)
summary(fit)$coef[2,c(1,4)]

fit <- lm(X10_12MeC32.estimate ~ chain.length , data=cl.data)
summary(fit)$coef[2,c(1,4)]


############################
## Supplementary Material ##
############################

## Univariate animal model ##
input <- perc.data
input[,6:30] <- scale(input[,6:30])

df <- data.frame(matrix(NA, ncol=3, nrow=25))
row.names(df) <- colnames(input)[6:length(input)]
colnames(df) <- c("estimate", "lower", "upper")
fixed.results <- list(model1.2=df,
                      model1.3=df,
                      model1.4=df)

random.results <- data.frame(matrix(NA, ncol=4, nrow=25))
row.names(random.results) <- colnames(input)[6:length(input)]
colnames(random.results) <- c("model1.1", "model1.2", "model1.3", "model1.4")

df <- data.frame(matrix(NA, ncol=3, nrow=25))
row.names(df) <- colnames(input)[6:length(input)]
colnames(df) <- c("estimate", "lower", "upper")
posterior.heritability <- list(model1.1=df,
                               model1.2=df,
                               model1.3=df,
                               model1.4=df)

for (x in 1:25){
  temp.data <- input[,c(1:5,x+5)]
  colnames(temp.data)[6] <- "CHC"
  
  pheno.var <- var(temp.data$CHC, na.rm=TRUE)
  
  prior1.1 <- list(G=list(G1=list(V=matrix(pheno.var/2),n=1)),
                   R=list(V=matrix(pheno.var/2),n=1))
  model1.1 <- MCMCglmm(CHC ~ 1, random= ~ animal, pedigree=pedigree, data=temp.data, nitt=60000, thin=50, burnin=10000, prior=prior1.1)
  
  prior1.2 <- list(G=list(G1=list(V=matrix(pheno.var/2),n=1)),
                   R=list(V=matrix(pheno.var/2),n=1))
  model1.2 <- MCMCglmm(CHC ~ age, random= ~ animal, pedigree=pedigree, data=temp.data, nitt=60000, thin=50, burnin=10000, prior=prior1.2)
  
  
  prior1.3 <- list(G=list(G1=list(V=matrix(pheno.var/3),n=1), 
                          G2=list(V=matrix(pheno.var/3),n=1)),
                   R=list(V=matrix(pheno.var/3),n=1))
  model1.3 <- MCMCglmm(CHC ~ age, random= ~ animal + mother, pedigree=pedigree, data=temp.data, nitt=60000, thin=50, burnin=10000, prior=prior1.3)
  
  prior1.4 <- list(G=list(G1=list(V=matrix(pheno.var/4),n=1), 
                          G2=list(V=matrix(pheno.var/4),n=1),
                          G3=list(V=matrix(pheno.var/4),n=1)),
                   R=list(V=matrix(pheno.var/4),n=1))
  model1.4 <- MCMCglmm(CHC ~ age, random= ~ animal + mother + strain, pedigree=pedigree, data=temp.data, nitt=60000, thin=50, burnin=10000, prior=prior1.4)
  
  random.results[x,1] <- round(model1.1$DIC,1)
  
  fixed.results$model1.2[x,1] <- round(posterior.mode(model1.2$Sol)["agenymph"],3)
  fixed.results$model1.2[x,c(2:3)] <- round(HPDinterval(model1.2$Sol)["agenymph",],3)
  random.results[x,2] <- round(model1.2$DIC,1)
  
  fixed.results$model1.3[x,1] <- round(posterior.mode(model1.3$Sol)["agenymph"],3)
  fixed.results$model1.3[x,c(2:3)] <- round(HPDinterval(model1.3$Sol)["agenymph",],3)
  random.results[x,3] <- round(model1.3$DIC,1)
  
  fixed.results$model1.4[x,1] <- round(posterior.mode(model1.4$Sol)["agenymph"],3)
  fixed.results$model1.4[x,c(2:3)] <- round(HPDinterval(model1.4$Sol)["agenymph",],3)
  random.results[x,4] <- round(model1.4$DIC,1)
  
  posterior.heritability1.1 <- model1.1$VCV[,"animal"]/(model1.1$VCV[,"animal"] + model1.1$VCV[,"units"])
  posterior.heritability$model1.1[x,1] <- round(posterior.mode(posterior.heritability1.1),3)
  posterior.heritability$model1.1[x,c(2:3)] <- round(HPDinterval(posterior.heritability1.1, 0.95),3)
  
  posterior.heritability1.2 <- model1.2$VCV[,"animal"]/(model1.2$VCV[,"animal"] + model1.2$VCV[,"units"])
  posterior.heritability$model1.2[x,1] <- round(posterior.mode(posterior.heritability1.2),3)
  posterior.heritability$model1.2[x,c(2:3)] <- round(HPDinterval(posterior.heritability1.2, 0.95),3)
  
  posterior.heritability1.3 <- model1.3$VCV[,"animal"]/(model1.3$VCV[,"animal"]+model1.3$VCV[,"mother"]+model1.3$VCV[,"units"])
  posterior.heritability$model1.3[x,1] <- round(posterior.mode(posterior.heritability1.3),3)
  posterior.heritability$model1.3[x,c(2:3)] <- round(HPDinterval(posterior.heritability1.3, 0.95),3)
  
  posterior.heritability1.4 <- model1.4$VCV[,"animal"]/(model1.4$VCV[,"animal"]+model1.4$VCV[,"mother"]+model1.4$VCV[,"strain"]+model1.4$VCV[,"units"])
  posterior.heritability$model1.4[x,1] <- round(posterior.mode(posterior.heritability1.4),3)
  posterior.heritability$model1.4[x,c(2:3)] <- round(HPDinterval(posterior.heritability1.4, 0.95),3)
}

write.csv(fixed.results, file="fixed results percentages.csv")
write.csv(random.results, file="random results percentages.csv")
write.csv(posterior.heritability, file="posterior heritability percentages.csv")


## Paternal half-sib covariance ##
phc.data <- read.table("CHC data FS HS.txt", header=T, row.names=1)

perc.phc.data <- phc.data
perc.phc.data[,3:length(perc.phc.data)] <- (perc.phc.data[,3:length(perc.phc.data)]/apply(perc.phc.data[,3:length(perc.phc.data)], 1, sum))*100

lr.phc.data <- phc.data
lr.phc.data[,3:length(lr.phc.data)] <- log(lr.phc.data[,3:length(lr.phc.data)]/apply(lr.phc.data[,3:length(lr.phc.data)], 1, geometric.mean))

input <- lr.phc.data

H2.reml.bootstrap <- function(x, sire, dam, reps) {
  reml <- lmer(x ~ (1|sire/dam), data=input)
  varcomps <- as.data.frame(VarCorr(reml))[,"vcov"]
  relvarcomps <- varcomps/sum(varcomps)
  H2 <- round(4*relvarcomps[2],5)
  newout <- c(among.s=round(varcomps[2],4),
              among.d=round(varcomps[1],4),
              within.d=round(varcomps[3],4),
              H2=H2)
  estimates <- matrix(numeric(reps), nrow=reps, ncol=4)
  data <- as.data.frame(cbind(x, sire, dam))
  dam.numbers <- length(levels(dam))  
  for (i in 1:reps){ 
    boot <- data.frame()
    dam.sample <- sample(dam.numbers,dam.numbers,replace=T)
    for (j in 1:length(dam.sample)){
      boot <- rbind(boot,data[data$dam==dam.sample[j],])
    }
    remlb <- lmer(x ~ (1|sire/dam), data=boot)
    varcompsb <- as.data.frame(VarCorr(remlb))[,"vcov"]
    relvarcompsb <- varcompsb/sum(varcompsb)
    H2b <- 4*relvarcompsb[2]
    estimates[i,] <- cbind(varcompsb[2],varcompsb[1],varcompsb[3],H2b) 
  }				
  est <- mean(estimates[,4])
  se <- sd(estimates[,4])/sqrt(length(estimates[,4]))
  up <- as.numeric(quantile(estimates[,4], probs=0.95))
  low <- as.numeric(quantile(estimates[,4], probs=0.05))
  results <- data.frame(Estimate=est, SE=se, Lower95CL=low, Upper95CL=up)
  return(results)
}

heritabilities <- data.frame(Compound=character(), H2=numeric(), SE=numeric(), Lower95CL=numeric(), Upper95CL=numeric())
for (p in 1:(length(input)-2)){
  res <- H2.reml.bootstrap(input[,p+2], input$sire, input$dam, 1000)
  result.per.compound <- data.frame(Compound=colnames(input)[p+2], H2=round(res[1],3), SE=round(res[2],3), Lower95CL=round(res[3],3), Upper95CL=round(res[4],3))
  heritabilities <- rbind(heritabilities, result.per.compound)
}


## Parent-offspring regression ##
data <- read.table("CHC data.txt", header=T)

mother <- data[1:36, 6:30]; row.names(mother) <- data$animal[1:36]
mother <- droplevels(mother)

father <- data[37:63, 6:30]; row.names(father) <- data$animal[37:63]
father <- rbind(father[1:9,], father)
father <- droplevels(father)

young <- data[64:351, 6:30]; row.names(young) <- data$animal[64:351]
young <- data.frame(father=rep(row.names(father), each=8), mother=rep(row.names(mother), each=8), young)
young <- droplevels(young)


perc.young <- young; perc.young[,3:27] <- (perc.young[,3:27]/apply(perc.young[,3:27], 1, sum))*100
perc.father <- father; perc.father <- (perc.father/apply(perc.father, 1, sum))*100
perc.mother <- mother; perc.mother <- (perc.mother/apply(perc.mother, 1, sum))*100

lr.young <- young; lr.young[,3:27] <- log(lr.young[,3:27]/apply(lr.young[,3:27], 1, geometric.mean))
lr.father <- father; lr.father <- log(lr.father/apply(lr.father, 1, geometric.mean))
lr.mother <- mother; lr.mother <- log(lr.mother/apply(lr.mother, 1, geometric.mean))


young.input <- lr.young
father.input <- lr.father
mother.input <- lr.mother

# mother
reps <- 1000
regression.reps <- data.frame(matrix(NA, nrow=reps, ncol=25))
colnames(regression.reps) <- colnames(mother.input)
for (h in 1:reps){
  mother.sample <- sample(row.names(mother.input), nrow(mother.input), replace=T)
  young.sample <- data.frame()
  for (j in 1:length(mother.sample)){
    young.sample = rbind(young.sample, young.input[young.input$mother==mother.sample[j],])
  }
  mother.sample <- mother.input[match(young.sample$mother, row.names(mother.input)),]
  
  regression <- data.frame(matrix(NA, nrow=25, ncol=2))
  colnames(regression) <- c("intercept", "slope")
  for (i in 1:25){
    fit <- lm(young.sample[,i+2] ~ mother.sample[,i])
    regression[i,] <- round(coef(fit),3)
    row.names(regression)[i] <- colnames(mother.sample)[i]
  }
  regression.reps[h,] <- regression[,2]
}
est <- apply(regression.reps, 2, mean)
se <- apply(regression.reps, 2, function(x) sd(x)/sqrt(length(x)))
up <- apply(regression.reps, 2, function(x) quantile(x, probs=0.95))
low <- apply(regression.reps, 2, function(x) quantile(x, probs=0.05))
results <- data.frame(Estimate=est, SE=se, Lower95CL=low, Upper95CL=up)

# father
reps <- 1000
regression.reps <- data.frame(matrix(NA, nrow=reps, ncol=25))
colnames(regression.reps) <- colnames(father.input)
for (h in 1:reps){
  father.sample <- sample(row.names(father.input), nrow(father.input), replace=T)
  young.sample <- data.frame()
  for (j in 1:length(father.sample)){
    young.sample = rbind(young.sample, young.input[young.input$father==father.sample[j],])
  }
  father.sample <- father.input[match(young.sample$father, row.names(father.input)),]
  
  regression <- data.frame(matrix(NA, nrow=25, ncol=2))
  colnames(regression) <- c("intercept", "slope")
  for (i in 1:25){
    fit <- lm(young.sample[,i+2] ~ father.sample[,i])
    regression[i,] <- round(coef(fit),3)
    row.names(regression)[i] <- colnames(father.sample)[i]
  }
  regression.reps[h,] <- regression[,2]
}
est <- apply(regression.reps, 2, mean)
se <- apply(regression.reps, 2, function(x) sd(x)/sqrt(length(x)))
up <- apply(regression.reps, 2, function(x) quantile(x, probs=0.95))
low <- apply(regression.reps, 2, function(x) quantile(x, probs=0.05))
results <- data.frame(Estimate=est, SE=se, Lower95CL=low, Upper95CL=up)

# midparent
reps <- 1000
regression.reps <- data.frame(matrix(NA, nrow=reps, ncol=25))
colnames(regression.reps) <- colnames(mother.input)
for (h in 1:reps){
  mother.sample <- sample(row.names(mother.input), nrow(mother.input), replace=T)
  young.sample <- data.frame()
  for (j in 1:length(mother.sample)){
    young.sample = rbind(young.sample, young.input[young.input$mother==mother.sample[j],])
  }
  mother.sample <- mother.input[match(young.sample$mother, row.names(mother.input)),]
  father.sample <- father.input[match(young.sample$mother, row.names(mother.input)),]
  
  regression <- data.frame(matrix(NA, nrow=25, ncol=2))
  colnames(regression) <- c("intercept", "slope")
  for (i in 1:25){
    fit <- lm(young.sample[,i+2] ~ as.vector((mother.sample[,i]+father.sample[,i])/2))
    regression[i,] <- round(coef(fit),3)
    row.names(regression)[i] <- colnames(mother.sample)[i]
  }
  regression.reps[h,] <- regression[,2]
}
est <- apply(regression.reps, 2, mean)
se <- apply(regression.reps, 2, function(x) sd(x)/sqrt(length(x)))
up <- apply(regression.reps, 2, function(x) quantile(x, probs=0.95))
low <- apply(regression.reps, 2, function(x) quantile(x, probs=0.05))
results <- data.frame(Estimate=est, SE=se, Lower95CL=low, Upper95CL=up)

