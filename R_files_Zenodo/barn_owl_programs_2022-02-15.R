
# procomp2: Estimate components of projection matrix (enhanced version from package 'lmf') ----

procomp2 <- function (formula, data)
{
  # Keep data as x
  x <- data.frame(data)
  # Keep call
  mf <- match.call(expand.dots = FALSE)
  # Match arguments with call
  m <- match(c("formula", "data"), names(mf), 0L)
  # Organize call
  mf <- mf[c(1L, m)]
  # Extract variables in formula
  alv <- all.vars(mf$formula)
  # Extract survival
  survival <- x[, alv[2]]
  # Check
  if(!all(unique(survival[!is.na(survival)]) %in% c(1, 0)))
    stop("Survival vector should only have elements 0 or 1. Check order of variables in formula")
  # Estimate projection matrix components
  setNames(stats:::aggregate.formula(formula, data = x, FUN = mean), c(all.vars(mf$formula[[3]]), all.vars(mf$formula[[2]])))
}



# promat2: Set-up of projection matrix (enhanced version from package 'lmf') ----

promat2 <- function (pc)
{
  # Extract number of age classes
  nage <- nrow(pc)
  # Set up diagonal matrix with survivals
  lcJ <- diag(pc[, 3], nrow = nage, ncol = nage)
  # If only one age class (i.e. all individuals of equal age)
  if (nage == 1) {
    # Set up 'projection matrix'
    l <- matrix(rbind(pc[, 2], lcJ), ncol = nage)
  }
  # If more than one age class
  else {
    # Set up projection matrix
    l <- matrix(rbind(pc[, 2], lcJ[-nage, ]), ncol = nage)
    # Insert survival for final age
    l[nage, nage] <- lcJ[nage, nage]
  }
  ret <- list(age = as.vector(pc[, 1]), l = l)
  ret
}



# eigenl2: Calculate lambda, u and v ----

eigenl2 <- function (pm) 
{
  
  if(is.list(pm)) {
    # Split pm
    age <- pm[[1]]
    l <- pm[[2]] 
  }
  else {
    age <- 1:ncol(pm)
    l <- pm
  }
  
  # Get lambda, u and v
  if (dim(l)[1] == dim(l)[2]) {
    ret <- list(age = age)
    ret$lambda <- as.numeric(eigen(l, only.values = TRUE)$values[1])
    ret$u <- abs(eigen(l)$vectors[, 1])
    ret$v <- abs(eigen(t(l))$vectors[, 1])
    ret$u <- ret$u/sum(ret$u)
    ret$v <- ret$v/sum(ret$u * ret$v)
  }
  else {
    ret <- list(lambda = colSums(l))
    ret$u <- 1
    ret$v <- 1
  }
  ret
}



# ind.rep.value: Calculate individual reproductive values ----
ind.rep.value <- function (fecundity, survival, age, v)
{

  # Arguments:
  # - fecundity, survival, age: Vectors with individual records of fec, surv and age
  # - v: Data frame with all age classes and reproductive values for each age class in descending order
  
  # Repeate the last value of v
  v <- rbind(v, tail(v, 1))
  # Now v is structured such that we always can choose v for the next age class in the individual reproductive value calculation below
  
  # Set up v2, a vector with individual records of the reproductive value in the next age class after the one in which an individual is now
  v2 <- left_join(data.frame(age = age), data.frame(age = head(v$age, -1), v = tail(v$v, -1)), by = "age")
  
  # Calculate individual reproductive value
  wj <- v$v[1] * fecundity + v2$v * survival
  
  # Return
  wj
  
}

#




# scalevar: Scale variable ----
scalevar <- function(x, center = TRUE, scale = "sd", na.rm = TRUE){
  if(center){
    x_uc <- x
    x <- x-mean(x, na.rm = na.rm)
  }
  if(scale == "sd"){
    x <- x/sd(x, na.rm = na.rm)
  }
  if(scale == "mean"){
    x <- x/mean(x_uc, na.rm = na.rm)
  }
  x
}

#



# se: Standard error ----

# Create function se() to calculate standard error (SE).
se <- function(x,
               na.rm = FALSE,
               ...)
{
  if(na.rm)
    x <- x[!is.na(x)]
  return(sqrt(var(x) / length(x)))
}
#



# pval_stars: Stars from p-values ----

pval_stars <- function(pval) {
  
  # Arguments:
  # - pval: vector or matrix with p-values
  
  # Value
  # Returns the object with stars for p-values
  
  # Create stars
  pval[pval < 0.001] <- "***"
  pval[pval < 0.01 & !pval %in% c("***", NA)] <- "**"
  pval[pval < 0.05 & !pval %in% c("***", "**", NA)] <- "*"
  pval[pval <= 0.1 & !pval %in% c("***", "**", "*", NA)] <- "."
  pval[pval > 0.1 & !pval %in% c("***", "**", "*", ".", NA)] <- "n.s."
  pval[is.na(pval)] <- ""
  
  # Return object
  pval
  
}
#


