
library("GenABEL")
library("RepeatABEL")


get.SEh2_genabel <- function(formula.FixedEffects= y~1, genabel.data, Va, Ve, GRM=NULL) {
#Version corrected 8 Sept 2016 by LRN
  if (is.null(GRM)) GRM = compute.GRM(genabel.data)
  #Check IDs in genotypes and phenotypes
  if (any(genabel.data@gtdata@idnames!=genabel.data@phdata$id)) stop("IDs of phenotypes and genotypes not in the same order")
  trait <- all.vars(formula.FixedEffects)[1]
  phenotype.data <- genabel.data@phdata
  y.all <- phenotype.data[, names(phenotype.data) %in% trait]
  n.all = length(y.all)
  X.all <- model.matrix(formula.FixedEffects, data = phenotype.data)
  test.missing <- !is.na(y.all) 
  y <- y.all[test.missing]
  n <- length(y)
  X <- X.all #BUG FIXED
  Z <- diag(n.all)[test.missing, ]
  ZGZ <- Z%*%GRM%*%t(Z) #BUG FIXED
  V <- ZGZ*Va + diag(n)*Ve
  tr <- function(X) sum(diag(X))
  invV <- solve(V)
  P <- invV - invV%*%X%*%solve( t(X)%*%invV%*%X )%*%t(X)%*%invV
  b <- c( Ve, -Va) / ( (Va + Ve)^2 )
  row1 <- cbind( tr( P%*%ZGZ%*%P%*%ZGZ), tr( P%*%ZGZ%*%P))
  row2 <- cbind( tr( P%*%P%*%ZGZ),  tr( P%*%P))
  C <- 0.5*rbind(row1, row2)
  if (!isSymmetric(C)) warning("Non-symmetric matrix constructed")
  b <- matrix(b,1,2)
  h2 <- Va/(Va + Ve)
  return( c(h2, sqrt( b%*%solve(C)%*%t(b) ) ) ) 
}


gwaa0.2 <- load.gwaa.data(pheno,geno,force=TRUE,makemap=F)

# When we have the VCEs, the heritability and SE are computed as
K = compute.GRM(gwaa0.2) #Computations the same as K = ibs(gen.data, w="freq") but returns a symmetric matrix


h2_SE <- get.SEh2_genabel(phenotype ~1, gwaa0.2, Va, Ve, GRM=K)

cat("Estimated heritability", h2_SE[1], "with SE", h2_SE[2], "\\n")

