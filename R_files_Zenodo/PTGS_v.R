# Functions for PTGS computations

# Copyright (C) 2012-2017 Juuso Parkkinen, Pekka Kohonen, Samuel Kaski & Roland GrafstrÃ¶m
# All rights reserved.
# 
# This program is open source software; you can redistribute it and/or modify it under the terms of the FreeBSD License (keep this notice): http://en.wikipedia.org/wiki/BSD_licenses
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# Write ranked gene lists for each experiment
WriteRankedGenes <- function(data, input.folder) {
  
  message("Writing ranked genes for ", nrow(data)," experiments to folder ", input.folder, "...", appendLF=FALSE)
  if (!is.null(rownames(data)))
    exp.names <- rownames(data)
  else
    exp.names <- paste0("experiment", 1:nrow(data))
  for (i in 1:nrow(data)) {
    ranked.dat <- sort(data[i,])
    filename <- file.path(input.folder, paste0(exp.names[i],".rnk"))
    write(c("GeneSymbols", "LogRatios"), ncol=2, sep="\\t", file=filename)
    write.table(ranked.dat, row.names=T, col.names=F, quote=F, sep="\\t", file=filename, append=T)  
  }
  message("DONE") 
}

# Estimate component probabilities for a new sample
EstimateComponents <- function(sparse.counts, alpha, phi, Niter) {
  
  C <- nrow(phi)
  Ndata <- length(sparse.counts)
  
  z <- sample(1:C, Ndata, replace=TRUE)
  n <- rep(0, C)
  temp <- table(z)
  n[as.numeric(names(temp))] <- temp
  logmeanprob <- rep(NA, Niter)
  for (i in 1:Niter) {
    if (i %% 10 == 0) cat(i,"..")
    for (d in 1:Ndata) {
      # Get word and subtract count
      n[z[d]] <- n[z[d]] -1
      # Compute probs for components
      probs <- rep(0, C)
      for (c in 1:C) 
        probs[c] <- (n[c] + alpha)*phi[c, sparse.counts[d]]
      
      # Sample new component and add count
      z[d] <- multinom.single(probs)
      n[z[d]] <- n[z[d]] +1
      logmeanprob[i] <- log(probs[z[d]]/sum(probs))
    }
    
  }
  cat("DONE\\n")
  #  theta <- (n + alpha)/ sum(n + C*alpha)
  theta <- (n + alpha)/ sum(n + alpha)
  return(list(z=z, n=n, theta=theta, logmeanprob=logmeanprob))
}

multinom.single <- function(prob) {
  cs <- cumsum(prob)
  which.max(runif(1) <= cs/cs[length(cs)])
}


# Get sparse counts from count matrix
GetSparseCounts <- function(counts) {
  
  sparse.counts <- c()
  if (max(counts)>0) {
    for (m in 0:(max(counts)-1)) {
      temp <- which((counts-m)>0)
      sparse.counts <- c(sparse.counts, temp)
    }
  }
  return(sparse.counts)
}

# Convert FDR q-values to counts
ConvertFDRtoCounts <- function(fdr.values, genesets, max.count) {
  
  counts <- c()
  for (direction in c("pos", "neg")) {
    
    # Convert FDR q-values to counts
    fdr.raw <- fdr.values[[direction]]
    fdr.df <- droplevels(subset(fdr.raw, NAME %in% genesets))
    counts.temp <- round(-log2(fdr.df$FDR.q.val)) -1
    counts.temp[which(counts.temp==-1 | counts.temp==-Inf)] <- 0
    counts.temp[which(counts.temp==Inf)] <- max.count
    
    # Match to genests and add to results
    counts.res <- rep(0, length(genesets))
    counts.res[match(fdr.df$NAME, genesets)] <- counts.temp
    counts <- c(counts, counts.res)
  } 
  return(counts)
}







