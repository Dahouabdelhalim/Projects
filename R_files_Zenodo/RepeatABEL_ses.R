library(RepeatABEL)

#version 24 May 2016
get.SEh2 <- function(formula.FixedEffects = y ~ 1, genabel.data, phenotype.data, id.name = "id", GWAS.output, GRM = NULL) {
  if (class(genabel.data) != "gwaa.data") 
    stop("The input of genabel.data is not a GenABEL object")
  if (is.null(genabel.data@phdata$id)) 
    stop("IDs not given as id in the phdata list")
  if (is.null(GRM)) GRM <- compute.GRM(genabel.data)
  trait <- all.vars(formula.FixedEffects)[1]
  y.all <- phenotype.data[, names(phenotype.data) %in% trait]
  phenotype.data <- phenotype.data[!is.na(y.all), ]
  id1 <- phenotype.data[, names(phenotype.data) %in% id.name]
  id2 <- genabel.data@phdata$id
  test1 <- id1 %in% id2
  test2 <- id2 %in% id1
  genabel.data <- genabel.data[test2, ]
  phenotype.data <- phenotype.data[test1, ]
  id1 <- phenotype.data[, names(phenotype.data) %in% id.name]
  id2 <- genabel.data@phdata$id
  N = length(id2)
  n = length(id1)
  indx <- numeric(n)
  for (i in 1:N) {
    indx <- indx + i * (id1 %in% id2[i])
  }
  Z.indx <- diag(N)[indx, ]
  y <- phenotype.data[, names(phenotype.data) %in% trait]
  X <- model.matrix(formula.FixedEffects, data = phenotype.data)
  eig <- eigen(GRM)
  non_zero.eigenvalues <- eig$values > (1e-06)
  eig$values[!non_zero.eigenvalues] <- 0
  #print("GRM ready")
  Z.GRM <- (eig$vectors %*% diag(sqrt(eig$values)))[indx, ]
  Z <- cbind(Z.GRM, Z.indx)
  mod1b <- GWAS.output@call$hglm
  V <- constructV(Z=Z,RandC = c(ncol(Z.GRM), ncol(Z.indx)), ratio=mod1b$varRanef/mod1b$varFix)
  V <- mod1b$varFix*V #Bug fixed by LRN 24 May 2016
  tr <- function(X) sum(diag(X))
  invV <- solve(V)
  P <- invV - invV%*%X%*%solve( t(X)%*%invV%*%X )%*%t(X)%*%invV
  b <- c( mod1b$varRanef[2]+mod1b$varFix, -mod1b$varRanef[1], -mod1b$varRanef[1]) / ( (sum(mod1b$varRanef) + mod1b$varFix)^2 )
  row1 <- cbind( tr( P%*%tcrossprod(Z.GRM)%*%P%*%tcrossprod(Z.GRM)), tr( P%*%tcrossprod(Z.GRM)%*%P%*%tcrossprod(Z.indx)), tr( P%*%tcrossprod(Z.GRM)%*%P))
  row2 <- cbind( tr( P%*%tcrossprod(Z.indx)%*%P%*%tcrossprod(Z.GRM)), tr( P%*%tcrossprod(Z.indx)%*%P%*%tcrossprod(Z.indx)), tr( P%*%tcrossprod(Z.indx)%*%P))
  row3 <- cbind( tr( P%*%P%*%tcrossprod(Z.GRM)), tr( P%*%P%*%tcrossprod(Z.indx)), tr( P%*%P))
  C <- 0.5*rbind(row1, row2, row3)
  if (!isSymmetric(C)) warning("Non-symmetric matrix constructed")
  b <- matrix(b,1,3)
  #print(b)
  #print(round(cov2cor(solve(C)),3))
  h2 <- mod1b$varRanef[1]/(sum(mod1b$varRanef)+mod1b$varFix)
  return( c(h2, sqrt( b%*%solve(C)%*%t(b) ) ) ) 
}

###########################################################################################################################
################
#Example 1
data(Phen.Data) #Phenotype data with repeated observations
data(gen.data) #GenABEL object including IDs and marker genotypes
GWAS1 <- rGLS(y ~ age + sex, genabel.data = gen.data, phenotype.data = Phen.Data)
plot(GWAS1, main="")
summary(GWAS1)
#Summary for variance component estimation without SNP effects
summary(GWAS1@call$hglm)


h2.SE <- get.SEh2(formula.FixedEffects = y ~ age + sex, genabel.data=gen.data, phenotype.data=Phen.Data, GWAS.output=GWAS1) 
cat("The estimated heritability is ", h2.SE[1], "with an approximate SE of ", h2.SE[2], ".","\\n" )
###########
#Example 2
#Simulate 4 observations per individual with a heritability of 0.5
set.seed(1234)
Phen.Sim <- simulate_PhenData(y ~ age, genabel.data=gen.data,
                              n.obs=rep(4, nids(gen.data)), SNP.eff=0, SNP.nr=1000, VC=c(1,0,10))
GWAS1 <- rGLS(y ~ age, genabel.data = gen.data, phenotype.data = Phen.Sim)
plot(GWAS1, main="Simulated Data Results")

h2.SE <- get.SEh2(formula.FixedEffects = y ~ age, genabel.data=gen.data, phenotype.data=Phen.Sim, GWAS.output=GWAS1) 
cat("The estimated heritability is ", h2.SE[1], "with an approximate SE of ", h2.SE[2], ".","\\n" )
###########


