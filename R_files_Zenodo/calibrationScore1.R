# this function calculates the calibration score and returns it

calculateCalibrationScore = function (assessments, realizations) {

  Nexperts = length(assessments)      
  Ncal = length(which(!is.na(realizations)))
  Nquantiles=ncol(assessments[[1]])

  # check realization against assessments using bucketing
  df = matrix(0, nrow=Nexperts, ncol=6)
    
    for (e in 1:Nexperts) {
      for (q in 1:Ncal) {
        if (realizations[q] <= assessments[[e]][q,][1]) {
          df[e, 1] = df[e, 1] + 1
        } else if (realizations[q] <= assessments[[e]][q,][2]) {  
          df[e, 2] = df[e, 2] + 1
        } else if (realizations[q] <= assessments[[e]][q,][3]) {
          df[e, 3] = df[e, 3] + 1
        } else if (realizations[q] <= assessments[[e]][q,][4]) {
          df[e, 4] = df[e, 4] + 1
        } else if (realizations[q] <= assessments[[e]][q,][5]) {
          df[e, 5] = df[e, 5] + 1
        } else {
          df[e, 6] = df[e, 6] + 1
        }
      }
    }
    
    # empirical probability vector, normalized
    se = list(array(0, dim=Nexperts))
    for (e in 1:Nexperts) {
      se[[e]] = c(as.vector(df[e,]))
      se[[e]] = se[[e]] / sum(se[[e]])
    }
    
    # theoretical probability vector
    p = c(0.05, 0.2, 0.25, 0.25, 0.2, 0.05)
    
    l = array(0, dim=Nexperts)
    for (e in 1:Nexperts) {
      for (i in 1:6) {
        if (se[[e]][i] > 1e-6) {
          l[e] = l[e] + se[[e]][i] * log(se[[e]][i] / p[i])
        }
      }
    }
    
    calibrationScores = 1 - pchisq(2 * Ncal * l, df=5)
  

  #print(calibrationScores)
  return(calibrationScores)
}

