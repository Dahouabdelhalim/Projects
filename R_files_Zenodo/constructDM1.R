# construct a DM (Decision Maker)

constructDM = function(assessments, weights_, DMx, L, U, quantiles) {

  Nexperts = length(assessments)      
  N = NROW(assessments[[1]])
  Nquantiles=ncol(assessments[[1]])
  
  numSamplePoints = 50000

  # get x coordinates if needed
  if (TRUE == is.null(DMx)) {
    DMx = list()
    for (q in 1:N) {
      x = list()
      for (e in 1:Nexperts) {
        for (quantile in 1:Nquantiles) {
          x = c(x, assessments[[e]][q, quantile])
        }
      }
      DMx[[q]] = c(L[q], unique(sort(unlist(x))), U[q])
    }
  }

  tmpDMy = list()
  isItemWeights = NCOL(weights_) > 1

  for (q in 1:N) {
    tmpDMy[[q]] = vector()
    for (x in DMx[[q]][c(-1,-length(DMx[[q]]))]) {
      y = 0
      for (e in 1:Nexperts) {
        tmp = assessments[[e]][q,]
        tmp = unname(tmp)
        if (isItemWeights) {
          weight = weights_[q, e]
        } else {
          weight = weights_[e]
        }
        if(any(diff(t(tmp))<=0))browser()
        y = y + unname(weight) * approx(c(L[q], tmp, U[q]), c(0, quantiles, 1), x, n=numSamplePoints)$y
      }
      tmpDMy[[q]] = c(tmpDMy[[q]], y)
    }
  }

  # add boundary values to y
  for (q in 1:N) {
    tmpDMy[[q]] = c(0, tmpDMy[[q]], 1)
  }

  # get assessments from DMs
  # use approx with x and y swapped to input the y
  DM = matrix(0, nrow=N, ncol=Nquantiles)
  for (q in 1:N) {
    for (i in 1:Nquantiles) {
      quantile = quantiles[i]
      if(any(diff(t(tmpDMy[[q]]))<=0))browser()
      DM[q, i] = approx(tmpDMy[[q]], DMx[[q]], quantile, n=numSamplePoints)$y
     }
  }

  return (DM)
}


