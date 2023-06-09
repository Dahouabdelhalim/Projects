loa_maker <- function(sample, bias,V_bias,logs2,V_logs2, rho) {
  n_sample <- length(unique(sample))
  n_ES <- length(sample)

  ##bias and tau2 estimation
  V_list_b = impute_covariance_matrix(vi = V_bias,  #known correlation vector
                                      cluster = sample,   #CitationID, #study ID
                                      r = rho) #assumed correlation

  MVmodel_b = rma.mv(yi = bias, #effect size
                     V = V_list_b, #variance
                     random = ~1 | sample, #nesting structure: highest / lowest
                     test= "t", #use t-tests
                     method="REML") #estimate variances using REML

  tau2_b = summary(MVmodel_b)$sigma2

  mvcf_b <- coef_test(MVmodel_b,#estimation model above
                      cluster=sample, #define cluster IDs
                      vcov = "CR2") #estimation method (CR2 is best)

  mean_b <- mvcf_b$beta
  se_b <- mvcf_b$SE

  ##sd2 estimation
  V_list_s = impute_covariance_matrix(vi = V_logs2,  #known correlation vector
                                      cluster = sample,   #CitationID, #study ID
                                      r = rho) #assumed correlation

  MVmodel_s = rma.mv(yi = logs2, #effect size
                     V = V_list_s, #variance
                     random = ~1 | sample, #nesting structure: highest / lowest
                     test= "t", #use t-tests
                     method="REML") #estimate variances using REML

  mvcf_s <- coef_test(MVmodel_s,#estimation model above
                      cluster=sample, #define cluster IDs
                      vcov = "CR2") #estimation method (CR2 is best)

  mean_logs2 <- mvcf_s$beta
  se_logs2 <- mvcf_s$SE

  ##### put it together
  bias_row <- c(n_sample,n_ES, mean_b, se_b, tau2_b)
  logs2_row <- c(mean_logs2, se_logs2)

  bias_mean <- bias_row[3]
  tau2_est <- bias_row[5]
  sd2_est <- exp(logs2_row[1])

  LOA_L <- bias_mean - 1.96*sqrt(sd2_est + tau2_est)
  LOA_U <- bias_mean + 1.96*sqrt(sd2_est + tau2_est)

  m <- bias_row[1]
  k <- bias_row[2]
  tcrit <- qt(1-.05/2,m-1)

  B1 <- sd2_est^2/(sd2_est + tau2_est)
  B2 <- 1/(sd2_est + tau2_est)
  V_logT2 <- 2/sum((V_bias + tau2_est)^(-2))

  #V_LOA_mod <- bias_row[4]^2 + B1*logs2_row[4] + B2*V_logT2
  V_LOA_rve <- bias_row[4]^2 + B1*logs2_row[2]^2 + B2*V_logT2

  #CI_L_mod <- LOA_L - tcrit*sqrt(V_LOA_mod)
  #CI_U_mod <- LOA_U + tcrit*sqrt(V_LOA_mod)

  CI_L_rve <- LOA_L - tcrit*sqrt(V_LOA_rve)
  CI_U_rve <- LOA_U + tcrit*sqrt(V_LOA_rve)

  #c(m, k, rho, bias_mean, sd2_est, tau2_est, LOA_L, LOA_U, CI_L_rve, CI_U_rve)
  c("nsample" = m, "nES" = k, "imputed_corr" = rho, "mean_bias" = bias_mean, "s2_est" = sd2_est, "tau2_est" = tau2_est, "LOA_L" = LOA_L, "LOA_U" = LOA_U, "CI_L" = CI_L_rve, "CI_U" = CI_U_rve)

}

