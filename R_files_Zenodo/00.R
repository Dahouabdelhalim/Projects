################################################################################
#                              Define functions                                #
#                                                                              #
#                                 Eli Strauss                                  #
#                                                                              #
#                                                                              #
#                               February 2019                                  #
################################################################################


elo_deviance <- function(intx, id, k = 20){
  if(nrow(intx) < 1){stop('No aggressions supplied')}
  intx <- intx[order(intx$date, intx$time),]
  intx$order <- 1:nrow(intx)
  #remove interactions with no rank
  intx <- intx[!is.na(intx$recip.rank) & !is.na(intx$agg.rank),]
  #remove interactions between ids with same rank
  intx <- intx[intx$agg.rank != intx$recip.rank,]
  if(nrow(intx) < 1){warning('No interactions with known, non-equal rank'); return(NA)}
  
  we_1 <- intx[intx$recip.rank - intx$agg.rank > 0, c('date', 'aggressor', 'recip', 'order')]
  we_2 <- intx[intx$recip.rank - intx$agg.rank < 0,c('date', 'recip', 'aggressor', 'order')]
  
  names(we_1) <- c('date', 'exp.winner', 'exp.loser', 'order')
  names(we_2) <- c('date','exp.winner', 'exp.loser', 'order')
  wins_expected <- rbind(we_1, we_2)
  wins_expected <- wins_expected[order(wins_expected$order),]
  
  #ensure that the focal individual is first in the matrix
  idents <- unique(c(intx$aggressor, intx$recip))
  idents <- c(idents[which(idents == id)], idents[-which(idents == id)])
  ##get elo scores for all individuals, although only valid for focal individual
  elo_observed <- aniDom::elo_scores(winners = intx$aggressor, losers = intx$recip, return.trajectories = T, identities = idents, K = k)
  elo_expected <- aniDom::elo_scores(winners = wins_expected$exp.winner, losers = wins_expected$exp.loser, return.trajectories = T, identities = idents, K = k)
  return(elo_observed[1,length(elo_observed[1,])] - elo_expected[1,length(elo_expected[1,])])
}

elo_obs <- function(intx, id, k = 20){
  if(nrow(intx) < 1){stop('No aggressions supplied')}
  intx <- intx[order(intx$date, intx$time),]
  intx$order <- 1:nrow(intx)
  #remove interactions with no rank
  intx <- intx[!is.na(intx$recip.rank) & !is.na(intx$agg.rank),]
  #remove interactions between ids with same rank
  intx <- intx[intx$agg.rank != intx$recip.rank,]
  if(nrow(intx) < 1){warning('No interactions with known, non-equal rank'); return(NA)}
  
  we_1 <- intx[intx$recip.rank - intx$agg.rank > 0, c('date', 'aggressor', 'recip', 'order')]
  we_2 <- intx[intx$recip.rank - intx$agg.rank < 0,c('date', 'recip', 'aggressor', 'order')]
  
  names(we_1) <- c('date', 'exp.winner', 'exp.loser', 'order')
  names(we_2) <- c('date','exp.winner', 'exp.loser', 'order')
  wins_expected <- rbind(we_1, we_2)
  wins_expected <- wins_expected[order(wins_expected$order),]
  
  #ensure that the focal individual is first in the matrix
  idents <- unique(c(intx$aggressor, intx$recip))
  idents <- c(idents[which(idents == id)], idents[-which(idents == id)])
  ##get elo scores for all individuals, although only valid for focal individual
  elo_observed <- aniDom::elo_scores(winners = intx$aggressor, losers = intx$recip, return.trajectories = T, identities = idents, K = k)
  elo_expected <- aniDom::elo_scores(winners = wins_expected$exp.winner, losers = wins_expected$exp.loser, return.trajectories = T, identities = idents, K = k)
  return(elo_observed[1,length(elo_observed[1,])])
}


elo_deviance_traj <- function(intx, id, k = 20){
  if(nrow(intx) < 1){stop('No aggressions supplied')}
  intx <- intx[order(intx$date, intx$time),]
  intx$order <- 1:nrow(intx)
  #remove interactions with no rank
  intx <- intx[!is.na(intx$recip.rank) & !is.na(intx$agg.rank),]
  #remove interactions between ids with same rank
  intx <- intx[intx$agg.rank != intx$recip.rank,]
  if(nrow(intx) < 1){warning('No interactions with known, non-equal rank'); return(NA)}
  
  we_1 <- intx[intx$recip.rank - intx$agg.rank > 0, c('date', 'aggressor', 'recip', 'order')]
  we_2 <- intx[intx$recip.rank - intx$agg.rank < 0,c('date', 'recip', 'aggressor', 'order')]
  
  names(we_1) <- c('date', 'exp.winner', 'exp.loser', 'order')
  names(we_2) <- c('date','exp.winner', 'exp.loser', 'order')
  wins_expected <- rbind(we_1, we_2)
  wins_expected <- wins_expected[order(wins_expected$order),]
  
  #ensure that the focal individual is first in the matrix
  idents <- unique(c(intx$aggressor, intx$recip))
  idents <- c(idents[which(idents == id)], idents[-which(idents == id)])
  ##get elo scores for all individuals, although only valid for focal individual
  elo_observed <- aniDom::elo_scores(winners = intx$aggressor, losers = intx$recip, return.trajectories = T, identities = idents, K = k)
  elo_expected <- aniDom::elo_scores(winners = wins_expected$exp.winner, losers = wins_expected$exp.loser, return.trajectories = T, identities = idents, K = k)
  return(elo_observed[1,] - elo_expected[1,])
}


### Function for producing tables of Cox hazards models
cox.tab <- function(mod, fixef.labs = NULL){
  if(is.null(fixef.labs))
    fixef.labs <- names(fixef(mod))
  hr <- jstable::coxme.display(mod, dec = 3)
  hr <- as.data.frame(hr$table)
  
  s <- summary(mod)
  n = mod$n[2]
  
  if(length(fixef.labs) == 1){
    hr.name <- ''
  }else{
    hr.name <- 'adj. '
  }
  df <- data.frame(Predictor = fixef.labs,
                         Estimate = round(fixef(mod), 3),
                         SE = round(diag(sqrt(vcov(mod))), 3),
                         "Hazard Ratio\\n(95% CI)" = hr[,paste0(hr.name,"HR(95%CI)")],
                         "P value" = hr[,paste0(hr.name,"P value")], check.names = FALSE, row.names = NULL)
              
        footnote = paste0("n = ", n, ';',
              '       Random effect of clan (variance) = ', round(VarCorr(mod)$clan, 4))
  
  # ran.df <- data.frame("Group" = names(ranef(mod)),
  #                      Variable = 'Intercept',
  #                      "Std Dev" = round(sqrt(VarCorr(mod)$clan), 3),
  #                      Variance = round(VarCorr(mod)$clan, 3),check.names = FALSE, row.names = NULL)
  #                  
  return(list(df, footnote))
}


glm.tab <- function(mod, fixef.labs = NULL){
  if(is.null(fixef.labs))
    fixef.labs <- names(fixef(mod))
  s <- summary(mod)
  n = nrow(mod@frame)
  
  df <- data.frame(Predictor = fixef.labs,
                   Estimate = round(fixef(mod), 3),
                   SE = round(diag(sqrt(vcov(mod))), 3),
                   "P value" = round(s$coefficients[,4], 4), check.names = FALSE, row.names = NULL)
  
  df[,'P value'][s$coefficients[,4] < 0.0001] <- '< 0.0001'
  
  
  footnote = paste0("n = ", n, ';',
                    '       Random effect of clan (variance) = ', round(VarCorr(mod)$clan[1], 4))
  
      
  return(list(df, footnote))
}

