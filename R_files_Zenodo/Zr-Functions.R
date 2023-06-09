Z.Vrel <- function(x,phy = NULL){ 
  require(geomorph)
  if(length(dim(x))==3){ 
    x <- two.d.array(x)
  }
  n <- nrow(x)
  if (!is.null(phy)) {
    phy.parts <- geomorph:::phylo.mat(x, phy)
    Ptrans <- phy.parts$D.mat %*% (diag(n) - matrix(1, n) %*% 
                                     crossprod(matrix(1, n), phy.parts$invC)/sum(phy.parts$invC))
    x <- Ptrans %*% x
  }
  eig.obs <- eigen(cov(x))$values
  eig.obs <- eig.obs[which(zapsmall(eig.obs)>0)] 
  p <- length(eig.obs)
  VNull <- (p+2)/((p*(n-1))+2)
  Vtrans <- 2*VNull-1
  ZN <- 0.5*log((1+Vtrans) / (1-Vtrans))
  Re.obs <- var(eig.obs) / (mean(eig.obs)^2*p)
  Re.obs.trans <- ((2*Re.obs)-1)
  if(Re.obs.trans==1){Re.obs.trans=0.999}
  if(Re.obs.trans==-1){Re.obs.trans=-0.999}
  Z.obs <- 0.5*log((1+Re.obs.trans) / (1-Re.obs.trans))
  ZR <- Z.obs + abs(ZN)
  ZR.var <- 1/(n-3)
  out <- list(Re.obs = Re.obs, Z.obs = Z.obs, ZR = ZR, ZR.var = ZR.var)
  class(out) <- "rel.eig"
  out
}


compare.ZR <- function(...,two.tailed = TRUE){
  dots <- list(...)
  tails <- if(two.tailed) 2 else 1
  if(length(dots) < 2) stop("At least two objects of class rel.eig are needed")
  is.rel.eig <- function(x) class(x) == "rel.eig"
  list.check <- sapply(1:length(dots), function(j) any(is.rel.eig(dots[[j]])))
  if(any(list.check == FALSE)) stop("Not all objects are class rel.eig")
  k <- length(list.check)
  list.names <- as.list(substitute(list(...)))[-1L]
  k.combn <- combn(k,2)
  list.re.obs <- sapply(1:k, function(j) dots[[j]]$Re.obs)
  list.ZRs <- sapply(1:k, function(j) dots[[j]]$ZR)
  list.vars <- sapply(1:k, function(j) dots[[j]]$ZR.var) 
  z12 <- sapply(1:ncol(k.combn), function(j){
    a <- k.combn[1,j]; b <- k.combn[2,j]
    r1 <- list.ZRs[a]; r2 <- list.ZRs[b]; var1 <- list.vars[a]; var2 <- list.vars[b]
    abs(r1-r2)/sqrt( var1+var2)
  })
  z12.p <- sapply(1:length(z12), function(j) pnorm(abs(z12[[j]]), lower.tail = FALSE) * tails)
  d <- rep(0,k); names(d) <- list.names
  D <-dist(d)
  z12.pw <- p12.pw <- D
  for(i in 1:length(z12)) z12.pw[i] <-z12[i]
  for(i in 1:length(z12)) p12.pw[i] <-z12.p[i]
  names(list.ZRs) <-list.names
  pairwise.z <- as.matrix(z12.pw)
  pairwise.P <- as.matrix(p12.pw)
  diag(pairwise.P) <- 1
  
  out <- list(RelEig.obs = list.re.obs, sample.ZR = list.ZRs, 
              Z.var = list.vars,
              pairwise.z = pairwise.z,
              pairwise.P = pairwise.P)
  class(out) <- "compare.rel.eig"
  out
}

print.compare.r.eig <- function(x,...){
  z <- x$sample.z
  z.pw <- x$pairwise.z
  p <- x$pairwise.P
  cat("\\nEffect sizes\\n\\n")
  print(z)
  cat("\\nEffect sizes for pairwise differences in rel.eig effect size\\n\\n")
  print(z.pw)
  cat("\\nP-values\\n\\n")
  print(p)
  invisible(x)
}

summary.compare.r.eig <- function(object, ...) print.compare.r.eig(object,...)