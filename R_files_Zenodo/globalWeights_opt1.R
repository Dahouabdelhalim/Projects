#performance-based weights function that would provide the weights given a set of calibration questions 
#corresponding to the best performing (combined score) DM

perfWeights_opt<-function(assessments,realizations)
{
  
  # check all alpha cutoff values for
  bestAlpha = 0.0
  bestCalScore = 0.0
  bestInfoScore = 0.0
  bestCombScore = 0.0
  bestWeights = 0
  alpha = 0.0   # cutoff value
  
  calibrationScores = calculateCalibrationScore(assessments,realizations)
  
  tmp = calculateInformationScore(assessments, realizations, k=0.1, bounds=NULL)
  informationScores = tmp[[1]][which(!is.na(realizations)),]
  L = tmp[[2]][which(!is.na(realizations))]
  U = tmp[[3]][which(!is.na(realizations))]
  
  informationScoresCalAvg = colMeans(informationScores) 
  
  # combined score
  combinedScores = calibrationScores * informationScoresCalAvg
  
  orderedCalScores = unique(sort(unname(calibrationScores)))
  numAlphas = length(orderedCalScores)
  
  for (a in 1:numAlphas) {
    alpha = orderedCalScores[a]
    expertsForDM = calibrationScores >= alpha
    tmpWeights = combinedScores * (calibrationScores >= alpha)
    tmpWeights = tmpWeights / sum(tmpWeights)

    # create optimized DM using best Alpha
    tmpDM = constructDM(assessments[expertsForDM], tmpWeights[expertsForDM], NULL, L, U, quantiles)
    
    # check DM scores
    tmpCalScore = calculateCalibrationScore(list(tmpDM), realizations)
    tmpInfoScore = mean(calculateInformationScore(list(tmpDM), realizations, k=0.1, bounds=NULL, L, U)[[1]])
    tmpCombScore = tmpCalScore * tmpInfoScore
    
    if (tmpCombScore > bestCombScore) {
      bestAlpha = alpha
      bestCalScore = tmpCalScore
      bestInfoScore = tmpInfoScore
      bestWeights = tmpWeights
      bestCombScore = tmpCombScore
      bestDM = tmpDM
    }
  }
  
  return(bestWeights)
}