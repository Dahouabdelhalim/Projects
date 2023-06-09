################################################################################ 

# SCRIPT 3: ANCESTRAL STATE ESTIMATION, PHYLOGENETIC SIGNAL ESTIMATES AND
# MODEL-FIT ESTIMATES.

# AUTHOR: Jeremy D. Wilson

# ARTICLE: Wilson et al. 2022 - Chronogram or phylogram for ancestral state
# estimation? Model-fit statistics indicate the branch lengths underlying a
# binary characterâ€™s evolution

################################################################################

# LOAD PACKAGES & IMPORT DATA---------------------------------------------------


require(tidyverse)
require(caper)
require(expm)
require(phytools)
require(geiger)
load("trees&characters.Rdata")


# ASR OF MARKOV CHARACTERS------------------------------------------------------


new.ace.discrete <- function(Data, Phy, Model){
  reconstruction <- ace(x = Data, phy = Phy, type = "discrete", model = Model)
  reconstruction$AIC <- AIC(reconstruction)
  return(reconstruction)
}


# ER on chronogram.

data$chrono_ER <- purrr::map2(data$tip_states,
                              data$chronogram,
                              new.ace.discrete,
                              "ER")


# ARD on chronogram.

data$chrono_ARD <- purrr::map2(data$tip_states,
                               data$chronogram,
                               new.ace.discrete,
                               "ARD")


# ER on phylogram.

data$phylo_ER <- purrr::map2(data$tip_states,
                             data$phylogram,
                             new.ace.discrete,
                             "ER")


# ARD on phylogram.

data$phylo_ARD <- purrr::map2(data$tip_states,
                              data$phylogram,
                              new.ace.discrete,
                              "ARD")


data <- hoist(data,
              chrono_ER,
              chrono_ER_loglik = "loglik",
              chrono_ER_AIC = "AIC",
              chrono_ER_marginals = "lik.anc",
              chrono_ER_rates = "rates"
)


data <- hoist(data,
              chrono_ARD,
              chrono_ARD_loglik = "loglik",
              chrono_ARD_AIC = "AIC",
              chrono_ARD_marginals = "lik.anc",
              chrono_ARD_rates = "rates"
)


data <- hoist(data,
              phylo_ER,
              phylo_ER_loglik = "loglik",
              phylo_ER_AIC = "AIC",
              phylo_ER_marginals = "lik.anc",
              phylo_ER_rates = "rates"
)


data <- hoist(data,
              phylo_ARD,
              phylo_ARD_loglik = "loglik",
              phylo_ARD_AIC = "AIC",
              phylo_ARD_marginals = "lik.anc",
              phylo_ARD_rates = "rates"
)


# ASRs OF HIDDEN RATES CHARACTERS------------------------------------------------


# ER on chronogram.

data$HRM_chrono_ER <- purrr::map2(data$HRM_tip_states,
                                  data$chronogram,
                                  new.ace.discrete,
                                  "ER")


# ARD on chronogram.

data$HRM_chrono_ARD <- purrr::map2(data$HRM_tip_states,
                                   data$chronogram,
                                   new.ace.discrete,
                                   "ARD")


# ER on phylogram.

data$HRM_phylo_ER <- purrr::map2(data$HRM_tip_states,
                                 data$phylogram,
                                 new.ace.discrete,
                                 "ER")


# ARD on phylogram.

data$HRM_phylo_ARD <- purrr::map2(data$HRM_tip_states,
                                  data$phylogram,
                                  new.ace.discrete,
                                  "ARD")


data <- hoist(data,
              HRM_chrono_ER,
              HRM_chrono_ER_loglik = "loglik",
              HRM_chrono_ER_AIC = "AIC",
              HRM_chrono_ER_marginals = "lik.anc",
              HRM_chrono_ER_rates = "rates"
)

data <- hoist(data,
              HRM_chrono_ARD,
              HRM_chrono_ARD_loglik = "loglik",
              HRM_chrono_ARD_AIC = "AIC",
              HRM_chrono_ARD_marginals = "lik.anc",
              HRM_chrono_ARD_rates = "rates"
              
)

data <- hoist(data,
              HRM_phylo_ER,
              HRM_phylo_ER_loglik = "loglik",
              HRM_phylo_ER_AIC = "AIC",
              HRM_phylo_ER_marginals = "lik.anc",
              HRM_phylo_ER_rates = "rates"
)

data <- hoist(data,
              HRM_phylo_ARD,
              HRM_phylo_ARD_loglik = "loglik",
              HRM_phylo_ARD_AIC = "AIC",
              HRM_phylo_ARD_marginals = "lik.anc",
              HRM_phylo_ARD_rates = "rates"
)


# ASRs OF AMPLIFIED HIDDEN RATES CHARACTERS-------------------------------------


# ER on chronogram.

data$HRMs_chrono_ER <- purrr::map2(data$HRMs_tip_states,
                                   data$chronogram,
                                   new.ace.discrete,
                                   "ER")


# ARD on chronogram.

data$HRMs_chrono_ARD <- purrr::map2(data$HRMs_tip_states,
                                    data$chronogram,
                                    new.ace.discrete,
                                    "ARD")


# ER on phylogram.

data$HRMs_phylo_ER <- purrr::map2(data$HRMs_tip_states,
                                  data$phylogram,
                                  new.ace.discrete,
                                  "ER")


# ARD on phylogram.

data$HRMs_phylo_ARD <- purrr::map2(data$HRMs_tip_states,
                                   data$phylogram,
                                   new.ace.discrete,
                                   "ARD")


data <- hoist(data,
              HRMs_chrono_ER,
              HRMs_chrono_ER_loglik = "loglik",
              HRMs_chrono_ER_AIC = "AIC",
              HRMs_chrono_ER_marginals = "lik.anc",
              HRMs_chrono_ER_rates = "rates"
)

data <- hoist(data,
              HRMs_chrono_ARD,
              HRMs_chrono_ARD_loglik = "loglik",
              HRMs_chrono_ARD_AIC = "AIC",
              HRMs_chrono_ARD_marginals = "lik.anc",
              HRMs_chrono_ARD_rates = "rates"
              
)

data <- hoist(data,
              HRMs_phylo_ER,
              HRMs_phylo_ER_loglik = "loglik",
              HRMs_phylo_ER_AIC = "AIC",
              HRMs_phylo_ER_marginals = "lik.anc",
              HRMs_phylo_ER_rates = "rates"
)

data <- hoist(data,
              HRMs_phylo_ARD,
              HRMs_phylo_ARD_loglik = "loglik",
              HRMs_phylo_ARD_AIC = "AIC",
              HRMs_phylo_ARD_marginals = "lik.anc",
              HRMs_phylo_ARD_rates = "rates"
)


# GENERATING AICc VALUES FOR ALL MODEL/TREE COMBOS------------------------------


AICc_mod <- function(Phy, Loglik, Mod){
  n_obs = length(Phy$tip.label)
  if(Mod == "ER"){
    DF <- 1 
  } else {
    DF <- 2
  }
  AICc <- -2 * Loglik + 2 * DF * (n_obs/(n_obs - DF - 1))
  return(AICc)
}


data$chrono_ER_AICc <- purrr::pmap_dbl(list(data$chronogram,
                                            data$chrono_ER_loglik,
                                            "ER"),
                                       AICc_mod)
data$chrono_ARD_AICc <- purrr::pmap_dbl(list(data$chronogram,
                                             data$chrono_ARD_loglik,
                                             "ARD"),
                                        AICc_mod)
data$phylo_ER_AICc <- purrr::pmap_dbl(list(data$chronogram,
                                           data$phylo_ER_loglik,
                                           "ER"),
                                      AICc_mod)
data$phylo_ARD_AICc <- purrr::pmap_dbl(list(data$chronogram,
                                            data$phylo_ARD_loglik,
                                            "ARD"),
                                       AICc_mod)
data$HRM_chrono_ER_AICc <- purrr::pmap_dbl(list(data$chronogram,
                                                data$HRM_chrono_ER_loglik,
                                                "ER"),
                                           AICc_mod)
data$HRM_chrono_ARD_AICc <- purrr::pmap_dbl(list(data$chronogram,
                                                 data$HRM_chrono_ARD_loglik,
                                                 "ARD"),
                                            AICc_mod)
data$HRM_phylo_ER_AICc <- purrr::pmap_dbl(list(data$chronogram,
                                               data$HRM_phylo_ER_loglik,
                                               "ER"),
                                          AICc_mod)
data$HRM_phylo_ARD_AICc <- purrr::pmap_dbl(list(data$chronogram,
                                                data$HRM_phylo_ARD_loglik,
                                                "ARD"),
                                           AICc_mod)
data$HRMs_chrono_ER_AICc <- purrr::pmap_dbl(list(data$chronogram,
                                                 data$HRMs_chrono_ER_loglik,
                                                 "ER"),
                                            AICc_mod)
data$HRMs_chrono_ARD_AICc <- purrr::pmap_dbl(list(data$chronogram,
                                                  data$HRMs_chrono_ARD_loglik,
                                                  "ARD"),
                                             AICc_mod)
data$HRMs_phylo_ER_AICc <- purrr::pmap_dbl(list(data$chronogram,
                                                data$HRMs_phylo_ER_loglik,
                                                "ER"),
                                           AICc_mod)
data$HRMs_phylo_ARD_AICc <- purrr::pmap_dbl(list(data$chronogram,
                                                 data$HRMs_phylo_ARD_loglik,
                                                 "ARD"),
                                            AICc_mod)


# GENERATING BIC VALUES FOR ALL MODEL/TREE COMBOS------------------------------


BIC_mod <- function(Phy, Loglik, Mod){
  n_obs = length(Phy$tip.label)
  if(Mod == "ER"){
    DF <- 1 
  } else {
    DF <- 2
  }
  BIC <- -2 * Loglik + DF * log(n_obs)
  return(BIC)
}

data$chrono_ER_BIC <- purrr::pmap_dbl(list(data$chronogram,
                                           data$chrono_ER_loglik,
                                           "ER"),
                                      BIC_mod)
data$chrono_ARD_BIC <- purrr::pmap_dbl(list(data$chronogram,
                                            data$chrono_ARD_loglik,
                                            "ARD"),
                                       BIC_mod)
data$phylo_ER_BIC <- purrr::pmap_dbl(list(data$chronogram,
                                          data$phylo_ER_loglik,
                                          "ER"),
                                     BIC_mod)
data$phylo_ARD_BIC <- purrr::pmap_dbl(list(data$chronogram,
                                           data$phylo_ARD_loglik,
                                           "ARD"),
                                      BIC_mod)
data$HRM_chrono_ER_BIC <- purrr::pmap_dbl(list(data$chronogram,
                                               data$HRM_chrono_ER_loglik,
                                               "ER"),
                                          BIC_mod)
data$HRM_chrono_ARD_BIC <- purrr::pmap_dbl(list(data$chronogram,
                                                data$HRM_chrono_ARD_loglik,
                                                "ARD"),
                                           BIC_mod)
data$HRM_phylo_ER_BIC <- purrr::pmap_dbl(list(data$chronogram,
                                              data$HRM_phylo_ER_loglik,
                                              "ER"),
                                         BIC_mod)
data$HRM_phylo_ARD_BIC <- purrr::pmap_dbl(list(data$chronogram,
                                               data$HRM_phylo_ARD_loglik,
                                               "ARD"),
                                          BIC_mod)
data$HRMs_chrono_ER_BIC <- purrr::pmap_dbl(list(data$chronogram,
                                                data$HRMs_chrono_ER_loglik,
                                                "ER"),
                                           BIC_mod)
data$HRMs_chrono_ARD_BIC <- purrr::pmap_dbl(list(data$chronogram,
                                                 data$HRMs_chrono_ARD_loglik,
                                                 "ARD"),
                                            BIC_mod)
data$HRMs_phylo_ER_BIC <- purrr::pmap_dbl(list(data$chronogram,
                                               data$HRMs_phylo_ER_loglik,
                                               "ER"),
                                          BIC_mod)
data$HRMs_phylo_ARD_BIC <- purrr::pmap_dbl(list(data$chronogram,
                                                data$HRMs_phylo_ARD_loglik,
                                                "ARD"),
                                           BIC_mod)


# CALCULATE FRITZ D FOR ALL CHARACTERS ON BOTH TREES----------------------------


fritz.d <- function(Data, Phy){
  d <- phylo.d(data.frame(names = names(Data),
                          states = Data),
               Phy,
               names.col = names,
               binvar = states
  )
  return(d$DEstimate)
}


data$Fritz_D_chronogram <- map2_dbl(data$tip_states,
                                    data$chronogram,
                                    fritz.d)
data$Fritz_D_phylogram <- map2_dbl(data$tip_states,
                                   data$phylogram,
                                   fritz.d)
data$HRM_Fritz_D_chronogram <- map2_dbl(data$HRM_tip_states,
                                        data$chronogram,
                                        fritz.d)
data$HRM_Fritz_D_phylogram <- map2_dbl(data$HRM_tip_states,
                                       data$phylogram,
                                       fritz.d)
data$HRMs_Fritz_D_chronogram <- map2_dbl(data$HRMs_tip_states,
                                         data$chronogram,
                                         fritz.d)
data$HRMs_Fritz_D_phylogram <- map2_dbl(data$HRMs_tip_states,
                                        data$phylogram,
                                        fritz.d)


# CALCULATE PAGEL'S LAMBDA FOR ALL CHARACTERS ON BOTH TREES---------------------


lambda <- function(Data, Phy){
  result <- fitDiscrete(Phy,
                        Data,
                        model = "ARD",
                        transform = "lambda"
  )
  return(result$opt$lambda)
}

data$lambda_chronogram <- map2_dbl(data$tip_states,
                                   data$chronogram,
                                   lambda)
data$lambda_phylogram <- map2_dbl(data$tip_states,
                                  data$phylogram,
                                  lambda)
data$HRM_lambda_chronogram <- map2_dbl(data$HRM_tip_states,
                                       data$chronogram,
                                       lambda)
data$HRM_lambda_phylogram <- map2_dbl(data$HRM_tip_states,
                                      data$phylogram,
                                      lambda)
data$HRMs_lambda_chronogram <- map2_dbl(data$HRMs_tip_states,
                                        data$chronogram,
                                        lambda)
data$HRMs_lambda_phylogram <- map2_dbl(data$HRMs_tip_states,
                                       data$phylogram,
                                       lambda)



# CALCULATE BORGES DELTA FOR ALL TREES, USING MARGINALS FROM ARD MODEL----------

# Functions needed to calculate Borge's delta, taken and modified from "code.R" 
# by mrborges23 -https://github.com/mrborges23/delta_statistic

nentropy <- function(prob) {
  k <- ncol(prob) # number of states
  prob[prob > 1 / k] <- prob[prob > 1 / k] / (1 - k) - 1 / (1 - k) # state entropies
  tent <- apply(prob, 1, sum) # node entropy
  
  # correct absolute 0/1
  tent[tent == 0] <- tent[tent == 0] + runif(1, 0, 1) / 10000
  tent[tent == 1] <- tent[tent == 1] - runif(1, 0, 1) / 10000
  
  return(tent)
}

lpalpha <- function(a, b, x, l0) { # log posterior alpha
  N <- length(x)
  lp <- N * (lgamma(a + b) - lgamma(a)) - a * (l0 - sum(log(x)))
  return(lp)
}

lpbeta <- function(a, b, x, l0) { # log posterior beta
  N <- length(x)
  lp <- N * (lgamma(a + b) - lgamma(b)) - b * (l0 - sum(log(1 - x)))
  return(lp)
}

mhalpha <- function(a, b, x, l0, se) { # metropolis hastings alpha
  a0 <- a
  a1 <- exp(rnorm(1, log(a0), se))
  
  r <- min(1, exp(lpalpha(a1, b, x, l0) - lpalpha(a0, b, x, l0)))
  
  while (is.na(r) == T) {
    a1 <- exp(rnorm(1, log(a0), se))
    r <- min(1, exp(lpalpha(a1, b, x, l0) - lpalpha(a0, b, x, l0)))
  }
  
  if (runif(1) < r) {
    return(a1)
  } else {
    return(a0)
  }
}

mhbeta <- function(a, b, x, l0, se) { # metropolis hastings beta
  b0 <- b
  b1 <- exp(rnorm(1, log(b0), se))
  
  r <- min(1, exp(lpbeta(a, b1, x, l0) - lpbeta(a, b0, x, l0)))
  
  while (is.na(r) == T) {
    b1 <- exp(rnorm(1, log(b0), se))
    r <- min(1, exp(lpbeta(a, b1, x, l0) - lpbeta(a, b0, x, l0)))
  }
  
  if (runif(1) < r) {
    return(b1)
  } else {
    return(b0)
  }
}

emcmc <- function(alpha, beta, x, l0, se, sim, thin, burn) {
  usim <- seq(burn, sim, thin)
  gibbs <- matrix(NA, ncol = 2, nrow = length(usim))
  p <- 1
  
  for (i in 1:sim) {
    alpha <- mhalpha(alpha, beta, x, l0, se)
    beta <- mhbeta(alpha, beta, x, l0, se)
    
    if (i == usim[p]) {
      gibbs[p, ] <- c(alpha, beta)
      p <- p + 1
    }
  }
  return(gibbs)
}

delta <- function(marginals, lambda0, se, sim, thin, burn) {
  ar <- marginals
  x <- nentropy(ar)
  mc1 <- emcmc(rexp(1), rexp(1), x, lambda0, se, sim, thin, burn)
  mc2 <- emcmc(rexp(1), rexp(1), x, lambda0, se, sim, thin, burn)
  mchain <- rbind(mc1, mc2)
  delta <- mean(mchain[, 2] / mchain[, 1])
  return(delta)
}

delta.fast <- function(Marginals) {
  lambda0 <- 0.1 # rate parameter of the proposal
  se <- 0.5 # standard deviation of the proposal
  sim <- 10000 # number of iterations
  thin <- 10 # we kept only each 10th iterate
  burn <- 100 # 100 iterates are burned-in
  Borgesdelta <- delta(Marginals, lambda0, se, sim, thin, burn)
  return(Borgesdelta)
}


data$Borges_D_chrono <- map_dbl(data$chrono_ARD_marginals,
                                delta.fast)

data$Borges_D_phylo <- map_dbl(data$phylo_ARD_marginals,
                               delta.fast)

data$HRM_Borges_D_chrono <- map_dbl(data$HRM_chrono_ARD_marginals,
                                    delta.fast)

data$HRM_Borges_D_phylo <- map_dbl(data$HRM_phylo_ARD_marginals,
                                   delta.fast)

data$HRMs_Borges_D_chrono <- map_dbl(data$HRMs_chrono_ARD_marginals,
                                     delta.fast)

data$HRMs_Borges_D_phylo <- map_dbl(data$HRMs_phylo_ARD_marginals,
                                    delta.fast)


# SAVE DATA---------------------------------------------------------------------


save(data, file = "trees&characters&metrics.RData")