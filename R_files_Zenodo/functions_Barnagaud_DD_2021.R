###################################################################
# Tests for phylogenetic and traits clustering vs overdispersion  #
###################################################################

TPQE <- function(df, dis, nrep = 999, popw = NULL, alter = "two-sided"){
  
  if(is.null(popw))
    popw <- apply(df, 2, sum) / sum(df)
  
  fun1 <- function(df1, dis){
    stati <- divc(df1, dis)
    statimean <- sum(popw * stati[, 1])
    vtot <- apply(sweep(sweep(df1, 2, apply(df1, 2, sum), "/"), 2, popw, "*"), 1, sum)
    statpop <- divc(cbind.data.frame(vtot), dis)[, 1]
    
    targstat <- (statpop - statimean) / statpop
    
    return(targstat)
    
  }
  
  valobs <- fun1(df, dis)
  valsim <- sapply(1:nrep, function(i) fun1(df[sample(1:nrow(df)), ], dis))
  
  test1 <- as.randtest(valsim, valobs, alter = alter)
  test1$call <- "TPQE"
  return(test1)
  
}

################################
# Test for phylogenetic signal #
################################

rtest.decdiv <- function(phy, freq, dis = NULL, nrep = 99, vranking = "complexity", ties.method = "average", option = 1:3, optiontest = NULL, tol = 1e-08)
{
  
  #*******************************************************************************#
  #                         Checking of the parameters                            #
  #*******************************************************************************#
  
  if(!is.vector(freq)) stop("freq must be a unique vector")
  if (!is.numeric(nrep) | nrep <= 1) 
    stop("Non convenient nrep")
  if(sum(freq) < tol) stop("empty sample")
  if(any(freq < -tol)) stop("negative values in df")
  
  #*******************************************************************************#
  #                               Basic notations                                 #
  #*******************************************************************************#
  
  freq[freq < tol] <- 0
  freq <- freq / sum(freq)
  
  nsp <- length(phy$leaves)
  nnodes <- length(phy$nodes)
  if(is.null(dis))
    dis <- as.dist(sqrt(2*matrix(1, nsp, nsp) - diag(rep(1, nsp))))
  
  #*******************************************************************************#
  #                               Node ranking                                    #
  #*******************************************************************************#
  
  complexity <- function(phy){   
    
    matno <- as.data.frame(matrix(0, nnodes, nnodes))
    rownames(matno) <- names(phy$nodes)
    names(matno) <- names(phy$nodes)
    pathnodes <- phy$path[-(1:nsp)]
    for(i in 1:nnodes){
      matno[pathnodes[[i]], i] <- 1
    }
    listno <- lapply(1:nnodes, function(i) names(matno)[matno[i, ] > 0])
    names(listno) <- names(phy$nodes)
    nbdes <- cbind.data.frame(lapply(phy$parts, function(x) prod(1:length(x))))
    compl <- lapply(listno, function(x) prod(nbdes[x]))
    compltab <- cbind.data.frame(compl)
    compltab <- cbind.data.frame(t(compltab))
    names(compltab) <- "complexity"
    return(compltab)
    
  }
  
  droot <- function(phy){
    roottab <- cbind.data.frame(phy$droot[-(1:nsp)])
    names(roottab) <- "droot"
    return(roottab)
  }
  
  if(is.numeric(vranking)){
    vrank <- as.data.frame(rank(vranking, ties.method = ties.method))
    names(vrank) <- "free"
  }
  else
    vrank <- sapply(vranking, function(x) rank(get(x)(phy), ties.method = ties.method))
  
  if(!any(option == 3))
    r1 <- length(option)
  else
    r1 <- length(option) + length(vranking) - 1
  
  #*******************************************************************************#
  #                       Field observations                                      #
  #*******************************************************************************#
  
  vobs <- decdiv(phy, freq, dis, tol = 1e-08)
  
  #*******************************************************************************#
  #                       Statistics for the four tests                           #
  #*******************************************************************************#    
  
  namnodes <- rownames(vobs)
  stat1 <- function(v){
    v <- v/sum(v)
    return(max(v))
  }
  stat2 <- function(v){
    v <- v/sum(v)
    fun1 <- function(m){
      return(abs(sum(sort(v)[1:m]) - m/nnodes))
    }
    return(max(unlist(lapply(1:nnodes, fun1))))
  }
  stat3 <- function(v){
    # this statistics has been sightly changed because the consideration of ties in Ollier et al.
    # was not explained, althought ties always happen with such a methodology. 
    funstat3 <- function(vrank1){
      v <- v/sum(v)
      return(sum(rank(vrank1, ties.method = ties.method)*v)/nnodes)
    }
    return(apply(vrank, 2, funstat3))
  }
  
  methods <- c("stat1", "stat2", "stat3")[option]
  
  #*******************************************************************************#
  #                      Statistics on field observations                         #
  #*******************************************************************************#
  
  statobs <- unlist(sapply(methods, function(x) get(x)(vobs[, 1])))
  
  #*******************************************************************************#
  #                              Permutation scheme                               #
  #*******************************************************************************#     
  
  funperm <- function(i){
    e <- sample(1:nsp)
    vtheo <- decdiv(phy, freq[e], as.dist(as.matrix(dis)[e, e]), tol = tol)
    stattheo <- unlist(sapply(methods, function(x) get(x)(vtheo[, 1])))
    return(stattheo)
  }
  
  tabsimu <- as.data.frame(t(cbind.data.frame(lapply(1:nrep, funperm))))
  rownames(tabsimu) <- paste("t", 1:nrep, sep="")
  if(r1 == 2 & methods[1] == "stat3")
    names(tabsimu) <- paste("stat3", names(tabsimu), sep=".")
  
  #*******************************************************************************#
  #                                     End                                       #
  #*******************************************************************************# 
  
  optiondefault <- c("greater", "greater", "two-sided", "two-sided", "two-sided")
  names(optiondefault) <- c("stat1", "stat2", "stat3.complexity", "stat3.droot", "stat3.free")
  
  if(r1 == 1)
  {
    if(!is.null(optiontest))
      return(as.randtest(obs = statobs, sim = tabsimu[, 1], alter = optiontest, call = "rtest.decdiv"))
    else
      return(as.randtest(tabsimu[, 1], statobs, alter = optiondefault[names(tabsimu)], call = "rtest.decdiv"))
  }
  
  if(!is.null(optiontest))
    return(as.krandtest(obs = statobs, sim = tabsimu, alter = optiontest, call = "rtest.decdiv"))
  else
    return(as.krandtest(obs = statobs, sim = tabsimu, alter = optiondefault[names(tabsimu)], 
                        call = "rtest.decdiv"))
  
}

decdiv <- function(phy, df, dis = NULL, tol = 1e-08){
  
  if(is.vector(df)){
    df <- cbind.data.frame(df)
  }
  if(!is.data.frame(df)) stop("df should be a data frame")
  if(any(apply(df, 2, sum)<tol)) stop("null column in df")
  if(any(df < -tol)) stop("negative values in df")
  df[df < tol] <- 0
  df <- as.data.frame(apply(df, 2, function(x) x/sum(x)))
  
  disc2 <- function(samples, dis = NULL, structures = NULL, tol = 1e-08) 
  {
    if (!inherits(samples, "data.frame")) 
      stop("Non convenient samples")
    if (any(samples < 0)) 
      stop("Negative value in samples")
    if (any(apply(samples, 2, sum) < 1e-16)) 
      stop("Empty samples")
    if (!is.null(dis)) {
      if (!inherits(dis, "dist")) 
        stop("Object of class 'dist' expected for distance")
      if (!is.euclid(dis)) 
        warning("Euclidean property is expected for distance")
      dis <- as.matrix(dis)
      if (nrow(samples) != nrow(dis)) 
        stop("Non convenient samples")
    }
    if (is.null(dis)) 
      dis <- (matrix(1, nrow(samples), nrow(samples)) - diag(rep(1, 
                                                                 nrow(samples)))) * sqrt(2)
    if (!is.null(structures)) {
      if (!inherits(structures, "data.frame")) 
        stop("Non convenient structures")
      m <- match(apply(structures, 2, function(x) length(x)), 
                 ncol(samples), 0)
      if (length(m[m == 1]) != ncol(structures)) 
        stop("Non convenient structures")
      m <- match(tapply(1:ncol(structures), as.factor(1:ncol(structures)), 
                        function(x) is.factor(structures[, x])), TRUE, 0)
      if (length(m[m == 1]) != ncol(structures)) 
        stop("Non convenient structures")
    }
    Structutil <- function(dp2, Np, unit) {
      if (!is.null(unit)) {
        modunit <- model.matrix(~-1 + unit)
        sumcol <- apply(Np, 2, sum)
        Ng <- modunit * sumcol
        lesnoms <- levels(unit)
      }
      else {
        Ng <- as.matrix(Np)
        lesnoms <- colnames(Np)
      }
      sumcol <- apply(Ng, 2, sum)
      Lg <- t(t(Ng)/sumcol)
      colnames(Lg) <- lesnoms
      Pg <- as.matrix(apply(Ng, 2, sum)/nbhaplotypes)
      rownames(Pg) <- lesnoms
      deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*% 
                                  dp2 %*% x))
      ug <- matrix(1, ncol(Lg), 1)
      dg2 <- t(Lg) %*% dp2 %*% Lg - 1/2 * (deltag %*% t(ug) + 
                                             ug %*% t(deltag))
      colnames(dg2) <- lesnoms
      rownames(dg2) <- lesnoms
      return(list(dg2 = dg2, Ng = Ng, Pg = Pg))
    }
    Diss <- function(dis, nbhaplotypes, samples, structures) {
      structutil <- list(0)
      structutil[[1]] <- Structutil(dp2 = dis, Np = samples, 
                                    NULL)
      
      ###
      diss <- list(as.dist(structutil[[1]]$dg2))
      fun1 <- function(x){
        y <- x
        y[y<tol] <- 0
        return(y)
      }
      diss <- lapply(diss, fun1)
      diss <- lapply(diss, function(x) sqrt(2*x))
      ###
      
      if (!is.null(structures)) {
        for (i in 1:length(structures)) {
          structutil[[i + 1]] <- Structutil(structutil[[1]]$dg2, 
                                            structutil[[1]]$Ng, structures[, i])
        }
        ###
        diss <- c(diss, tapply(1:length(structures), factor(1:length(structures)),
                               function(x) (as.dist(structutil[[x + 1]]$dg2))))
        diss <- lapply(diss, fun1)
        diss <- lapply(diss, function(x) sqrt(2*x))
        ###
      }
      return(diss)
    }
    nbhaplotypes <- sum(samples)
    diss <- Diss(dis^2, nbhaplotypes, samples, structures)
    names(diss) <- c("samples", names(structures))
    if (!is.null(structures)) {
      return(diss)
    }
    return(diss$samples)
  }
  
  decdivV <- function(freq){
    nsp <- length(phy$leaves)
    nnodes <- length(phy$nodes)
    matno <- as.data.frame(matrix(0, nnodes, nsp))
    rownames(matno) <- names(phy$nodes)
    names(matno) <- names(phy$leaves)
    for(i in 1:nsp){
      matno[phy$path[[i]][-length(phy$path[[i]])], i] <- 1
    }
    matfr <- as.matrix(matno) %*% diag(freq)
    matfr2 <- as.data.frame(t(matfr))
    divno <- divc(matfr2, dis)
    matfr3 <- cbind.data.frame(matfr2, diag(freq))
    names(matfr3) <- c(names(matfr2), names(phy$leaves))
    matfr4 <- matfr3[, apply(matfr3, 2, sum)!=0]
    if(ncol(matfr4)==0) stop("only one species considered")
    discno <- disc2(matfr4, dis, tol = tol)
    lambdano <- apply(matfr4, 2, sum)
    prdist <- diag(lambdano)%*%as.matrix(discno^2/2)%*%diag(lambdano)
    colnames(prdist) <- rownames(prdist) <- names(matfr4)
    fun1 <- function(x){
      x <- x[apply(matfr3[x], 2, sum)!=0]
      if(length(x) == 1) return(0)
      else return(sum(prdist[x, x])/2)
    }
    res <- unlist(lapply(phy$parts, fun1))
    lambdano <- apply(matfr3, 2, sum)
    lambdano[lambdano < tol] <- 1
    res <- res * 1/as.vector(lambdano)[1:nnodes]
    return(res)
  }
  return(apply(df, 2, decdivV))
}
