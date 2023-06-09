popsim <- function(datpop,nj,popsize,boot,nrep,npop,nloci) {
  stats <- list()
  shannon <- matrix(0,nrep,npop)
  ehet <- matrix(0,nrep,npop)
  nallele <- matrix(0,nrep,npop)
  for (t in 1:nrep) {
    if (t%%100==0) {
      print(t)
    }
    for (j in 1:npop) {
      whichrows <- sample.int(nj[j],popsize,replace=boot)
      tmp <- datpop[[j]][whichrows,3:ncol(datpop[[j]])]
      for (k in 1:nloci) {
        rng <- c(2*(k-1)+1,2*k)
        alleles <- c(as.matrix(tmp[,rng]))
        prb <- as.matrix(table(alleles))
        prb <- prb/sum(prb)
        shannon[t,j] <- shannon[t,j] + sum(-prb*log(prb))
        ehet[t,j] <- ehet[t,j] + sum(prb^2)
        nallele[t,j] <- length(prb) + nallele[t,j]
      }
    }
  }
  stats[[1]] <- shannon
  stats[[2]] <- ehet
  stats[[3]] <- nallele
  return(stats)
}