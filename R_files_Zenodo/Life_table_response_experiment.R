##########################################
###   LIFE TABLE RESPONSE EXPERIMENT   ###
##########################################

# Whole-island ------------------------------------------------------------

# creating treatment and m-prime matrices ---------------------------------

horse_matrix <- 
  function(demographic_rates){
    # Define horse life stages
    stages <- c("F_foal", "F_1st_yr","F_2nd_yr", "F_3rd_yr", "F_Adt", "M_foal", "M_1st_yr","M_2nd_yr", "M_3rd_yr", "M_Adt")
    
    # Build the 10x10 matrix
    result <- 
      matrix(c(
        # top row
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        # second row
        demographic_rates$F_foal_surv,0, 0, 0, 0, 0, 0, 0, 0, 0,
        # third row
        0, demographic_rates$F_yearling_surv, 0, 0, 0, 0, 0, 0, 0, 0,
        # fourth row
        0, 0, demographic_rates$F_2yr_surv, 0, 0, 0, 0, 0, 0, 0,
        #fifth row
        0, 0, 0, demographic_rates$F_3yr_surv, demographic_rates$F_adult_surv, 0, 0, 0, 0, 0,
        # sixth row
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        # seventh row
        0, 0, 0, 0, 0, demographic_rates$M_foal_surv, 0, 0, 0, 0,
        # eighth row
        0, 0, 0, 0, 0, 0, demographic_rates$M_yearling_surv, 0, 0, 0,
        # ninth row
        0, 0, 0, 0, 0, 0, 0, demographic_rates$M_2yr_surv, 0, 0,
        #tenth row
        0, 0, 0, 0, 0, 0, 0, 0, demographic_rates$M_3yr_surv, demographic_rates$M_adult_surv),
        nrow = length(stages), byrow = TRUE,
        dimnames = list(stages, stages))
    result
  }

#constructing treatment matrix based on the raw demographic rates

surv_rates <- data.frame(F_foal_surv = 0.796,
                         F_yearling_surv = 0.896,
                         F_2yr_surv = 0.817,
                         F_3yr_surv = 0.963,
                         F_adult_surv =0.879,
                         M_foal_surv = 0.785,
                         M_yearling_surv = 0.936,
                         M_2yr_surv = 0.881,
                         M_3yr_surv = 0.975,
                         M_adult_surv = 0.927)

horse_treat <- list(F_foal_surv = surv_rates$F_foal_surv,
                    F_yearling_surv = surv_rates$F_yearling_surv,
                    F_2yr_surv = surv_rates$F_2yr_surv,
                    F_3yr_surv = surv_rates$F_3yr_surv,
                    F_adult_surv = surv_rates$F_adult_surv,
                    M_foal_surv = surv_rates$M_foal_surv,
                    M_yearling_surv = surv_rates$M_yearling_surv,
                    M_2yr_surv = surv_rates$M_2yr_surv,
                    M_3yr_surv = surv_rates$M_3yr_surv,
                    M_adult_surv = surv_rates$M_adult_surv,
                    h = 1.714,
                    k = 1,
                    FSR = 0.5)

#now creating the m-prime matrix

horse_mprime_male <- list(F_foal_surv = (surv_rates$F_foal_surv + surv_rates$M_foal_surv)/2,
                          F_yearling_surv = (surv_rates$F_yearling_surv + surv_rates$M_yearling_surv)/2,
                          F_2yr_surv = (surv_rates$F_2yr_surv + surv_rates$M_2yr_surv)/2,
                          F_3yr_surv = (surv_rates$F_3yr_surv + surv_rates$M_3yr_surv)/2,
                          F_adult_surv = (surv_rates$F_adult_surv + surv_rates$M_adult_surv)/2,
                          M_foal_surv = (surv_rates$M_foal_surv + surv_rates$M_foal_surv)/2,
                          M_yearling_surv = (surv_rates$M_yearling_surv + surv_rates$M_yearling_surv)/2,
                          M_2yr_surv = (surv_rates$M_2yr_surv + surv_rates$M_2yr_surv)/2,
                          M_3yr_surv = (surv_rates$M_3yr_surv + surv_rates$M_3yr_surv)/2,
                          M_adult_surv = (surv_rates$M_adult_surv + surv_rates$M_adult_surv)/2,
                          h = (1.714*2)/2,
                          k = 1,
                          FSR = 0.5)

horse_treatment_matrix <- horse_matrix(horse_treat)

horse_mprime_matrix <- horse_matrix(horse_mprime_male)

# Estimating the ASR for treatment and m-prime ----------------------------

matrix_ASR <-
  function(M, n = rep(10, nrow(M)), h = 1, k = 1, 
           iterations = 1000, FSR = 0.5, plot = FALSE){
    # Number of stages in matrix
    x <- length(n)
    # Number of time steps to simulate
    t <- iterations
    # an empty t by x matrix to store the stage distributions
    stage <- matrix(numeric(x * t), nrow = x)
    # an empty t vector to store the population sizes
    pop <- numeric(t)
    # for loop that goes through each of t time steps
    for (i in 1:t) {
      # stage distribution at time t
      stage[,i] <- n
      # population size at time t
      pop[i] <- sum(n)
      # number of male adults at time t
      M2 <- stage[10, i]
      # number of female adults at time t
      F2 <- stage[5, i]
      # Female freq-dep fecundity of female foals
      M[1,x/2]        <- min(0.5*0.58, k*M2/(M2 + ((0.58*F2)/h)))*(1-FSR)
      # Female freq-dep fecundity of male foals
      M[(x/8)*5,x/2]  <- min(0.5*0.58, k*M2/(M2 + ((0.58*F2)/h)))*FSR
      # Male freq-dep fecundity of female foals
      M[1,x]          <- min(0.5*(0.58*F2)/M2, k*(0.58*F2)/(M2 + ((0.58*F2)/h)))*(1-FSR)
      # Male freq-dep fecundity of male foals
      M[(x/10)*6,x]    <- min(0.5*(0.58*F2)/M2, k*(0.58*F2)/(M2 + ((0.58*F2)/h)))*FSR
      # define the new n (i.e., new stage distribution at time t)
      n <- M %*% n
      # define rownames of stage matrix
      rownames(stage) <- rownames(M)
      # define colnames of stage matrix
      colnames(stage) <- 0:(t - 1)
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, t]
    }
    # calc ASR as the proportion of the adult stable stage class that is male
    ASR <- stable.stage[x]/(stable.stage[x/2] + stable.stage[x])
    
    if(plot)
    {
      # plot distribution to assure that it is not chaotic
      matplot(rownames(t(stage)), t(stage), type='l', lwd=2, las=1)
    }
    # make a list of results
    pop.proj <- list(ASR = ASR,
                     stable.stage = stable.stage,
                     stage.vectors = stage,
                     SSD_M2 = stable.stage[10],
                     SSD_F2 = stable.stage[5])
    # print the list as output to the function
    pop.proj
  }

#treatment first
horse_treatment_ASR_analysis <- matrix_ASR(M = horse_treatment_matrix, h = horse_treat$h, FSR = horse_treat$FSR, iterations = 1000)

horse_treatment_ASR <- horse_treatment_ASR_analysis$ASR

horse_treatment_ASR

#now the mprime
horse_mprime_ASR_analysis <- matrix_ASR(M = horse_mprime_matrix, h = horse_mprime_male$h, FSR = horse_mprime_male$FSR, iterations = 1000)

horse_mprime_ASR <- horse_mprime_ASR_analysis$ASR

horse_mprime_ASR

# Sensitivity analysis ----------------------------------------------------

matrix_structure <- expression(
  # top row of matrix
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  # second row of matrix
  F_foal_surv, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  # third row of matrix
  0, F_1st_yr_surv, 0, 0, 0, 0, 0, 0, 0, 0,
  # fourth row of matrix
  0, 0, F_2yr_surv, 0, 0, 0, 0, 0, 0, 0,
  #fifth row of matrix
  0, 0, 0, F_3yr_surv, F_Adt_surv, 0, 0, 0, 0, 0,
  # sixth row
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  # seventh row
  0, 0, 0, 0, 0, M_foal_surv, 0, 0, 0, 0,
  # eighth row
  0, 0, 0, 0, 0, 0, M_1st_yr_surv, 0, 0, 0,
  # ninth row
  0, 0, 0, 0, 0, 0, 0, M_2yr_surv, 0, 0,
  # tenth row
  0, 0, 0, 0, 0, 0, 0, 0, M_3yr_surv, M_Adt_surv)

sensitivity_analysis <-
  function(vital_rates, matrix_str, h, k, FSR, niter = 1000, ASR){
    
    # make a list of all parameters
    vr <- list(F_foal_surv = vital_rates$F_foal_surv,
               F_1st_yr_surv = vital_rates$F_yearling_surv,
               F_2yr_surv = vital_rates$F_2yr_surv,
               F_3yr_surv = vital_rates$F_3yr_surv,
               F_Adt_surv = vital_rates$F_adult_surv,
               M_foal_surv = vital_rates$M_foal_surv,
               M_1st_yr_surv = vital_rates$M_yearling_surv,
               M_2yr_surv = vital_rates$M_2yr_surv,
               M_3yr_surv = vital_rates$M_3yr_surv,
               M_Adt_surv = vital_rates$M_adult_surv)
    # number of stages in the matrix
    no_stages <- sqrt(length(matrix_str))
    
    # Define horse
    stages <- c("F_foal", "F_1st_yr", "F_2yr", "F_3yr", "F_Adt", "M_foal", "M_1st_yr", "M_2yr", "M_3yr", "M_Adt")
    
    # an empty t by x matrix
    stage <- matrix(numeric(no_stages * niter), nrow = no_stages)
    
    # an empty t vector to store the population sizes
    pop <- numeric(niter)
    
    # dataframe to store the perturbation results
    ASR_pert_results <-
      data.frame(parameter = c("F_foal_surv", "F_1st_yr_surv","F_2yr_surv", "F_3yr_surv", "F_Adt_surv",
                               "M_foal_surv", "M_1st_yr_surv", "M_2yr_surv", "M_3yr_surv", "M_Adt_surv",
                               "h", "FSR"),
                 sensitivities = numeric(12),
                 elasticities = numeric(12))
    
    # specifiy how many survival rates there are
    n <- length(vr)
    
    # create vectors of perturbations to test on parameters of the matrix model
    vr_nums <- seq(0, 1, 0.01) # prop changes in survival and FSR (i.e., between 0 and 1)
    h_nums <- seq(0, 2, 0.02) # proportional changes in h index (i.e., between 0 and 2)
    
    # create empty dataframes to store the perturbation results for ASR
    vr_pert_ASR <- matrix(numeric(n * length(vr_nums)),
                          ncol = n, dimnames = list(vr_nums, names(vr)))
    h_pert_ASR <- matrix(numeric(length(h_nums)),
                         ncol = 1, dimnames = list(h_nums, "h"))
    # k_pert_ASR <- matrix(numeric(length(k_nums)),
    #                      ncol = 1, dimnames = list(k_nums, "k"))
    FSR_pert_ASR <- matrix(numeric(length(vr_nums)),
                           ncol = 1, dimnames = list(vr_nums, "FSR"))
    
    # perturbation of vital rates survival rates
    for (g in 1:n) # pick a column (i.e., a variable)
    {
      vr2 <- vr # reset the vital rates to the original
      for (i in 1:length(vr_nums)) # pick a perturbation level
      {
        vr2[[g]] <- vr_nums[i] # specify the vital rate with the new perturbation level
        A <- matrix(sapply(matrix_str, eval, vr2, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        # reset the starting stage distribution for simulation (all with 10 individuals)
        m <- rep(10, no_stages) 
        for (j in 1:niter) { # project the matrix through t iteration
          # stage distribution at time t
          stage[,j] <- m
          # population size at time t
          pop[j] <- sum(m)
          # number of male adults at time t
          M2 <- stage[10, j]
          # number of female adults at time t
          F2 <- stage[5, j]
          # Female freq-dep fecundity of female foals
          A[1,no_stages/2]        <- min((0.5*0.58), k*M2/(M2 + ((F2*0.58)/h)))*(1-FSR)
          # Female freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.58), k*M2/(M2 + ((F2*0.58)/h)))*FSR
          # Male freq-dep fecundity of female foals
          A[1,no_stages]          <- min(0.5*(0.58*F2)/M2, k*(0.58*F2)/(M2 + ((0.58*F2)/h)))*(1-FSR)
          # Male freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.58*F2)/M2, k*(0.58*F2)/(M2 + ((0.58*F2)/h)))*FSR
          # define the new n (i.e., new stage distribution at time t)
          m <- A %*% m
        }
        # define rownames of stage matrix
        rownames(stage) <- rownames(A)
        # define colnames of stage matrix
        colnames(stage) <- 0:(niter - 1)
        # calculate the proportional stable stage distribution
        stage <- apply(stage, 2, function(x) x/sum(x))
        # define stable stage as the last stage
        stable.stage <- stage[, niter]
        # calc ASR as the proportion of the adult stable stage class that is male
        vr_pert_ASR[i, g] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                        stable.stage[no_stages])
      }
      # get the spline function of ASR
      spl_ASR <- smooth.spline(vr_pert_ASR[,g] ~ rownames(vr_pert_ASR))
      # estimate the slope of the tangent of the spline at the vital rate
      ASR_pert_results[g, 2] <- predict(spl_ASR, x=vr[[g]], deriv=1)$y
      # re-scale sensitivity into elasticity
      ASR_pert_results[g, 3] <- vr[[g]]/ASR * ASR_pert_results[g, 2]
    }
    # perturbation of the h index parameter
    for (i in (1:length(h_nums))) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                  nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                  dimnames = list(stages, stages)) # build the matrix with the new value
      # reset the starting stage distribution for simulation (all with 10 individuals)
      m <- rep(10, no_stages) 
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        # population size at time t
        pop[j] <- sum(m)
        # number of male adults at time t
        M2 <- stage[10, j]
        # number of female adults at time t
        F2 <- stage[5, j]
        # Female freq-dep fecundity of female foals
        A[1,no_stages/2]        <- min((0.5*0.58), k*M2/(M2 + ((0.58*F2)/h_nums[i])))*(1-FSR)
        # Female freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages/2] <- min((0.5*0.58), k*M2/(M2 + ((0.58*F2)/h_nums[i])))*FSR
        # Male freq-dep fecundity of female foals
        A[1,no_stages]          <- min(0.5*(0.58*F2)/M2, k*(0.58*F2)/(M2 + ((0.58*F2)/h_nums[i])))*(1-FSR)
        # Male freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages] <- min(0.5*(0.58*F2)/M2, k*(0.58*F2)/(M2 + ((0.58*F2)/h_nums[i])))*FSR
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% m
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      # calc ASR as the proportion of the adult stable stage class that is male
      h_pert_ASR[i,] <- 
        stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
    }
    # get the spline function of ASR
    spl_ASR <- smooth.spline(h_pert_ASR[, 1] ~ rownames(h_pert_ASR))
    # estimate the slope of the tangent of the spline at the vital rate
    ASR_pert_results[n+1, 2] <- predict(spl_ASR, x=h, deriv=1)$y
    # re-scale sensitivity into elasticity
    ASR_pert_results[n+1, 3] <- h/ASR * ASR_pert_results[n+1, 2]
    # perturbation of FSR
    for (i in 1:length(vr_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                  nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                  dimnames = list(stages, stages)) # build the matrix with the new value
      # reset the starting stage distribution for simulation (all with 10 individuals)
      m <- rep(10, no_stages) 
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        # population size at time t
        pop[j] <- sum(m)
        # number of male adults at time t
        M2 <- stage[10, j]
        # number of female adults at time t
        F2 <- stage[5, j]
        # Female freq-dep fecundity of female foals
        A[1,no_stages/2]        <- min((0.5*0.58), k*M2/(M2 + ((0.58*F2)/h)))*(1-vr_nums[i])
        # Female freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages/2] <- min((0.5*0.58), k*M2/(M2 + ((0.58*F2)/h)))*vr_nums[i]
        # Male freq-dep fecundity of female foals
        A[1,no_stages]          <-min(0.5*(0.58*F2)/M2, k*(0.64*F2)/(M2 + ((0.58*F2)/h)))*(1-vr_nums[i])
        # Male freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages] <- min(0.5*(0.58*F2)/M2, k*(0.58*F2)/(M2 + ((0.58*F2)/h)))*vr_nums[i]
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% m
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      # calc ASR as the proportion of the adult stable stage class that is male
      FSR_pert_ASR[i,] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                     stable.stage[no_stages])
    }
    # get the spline function of ASR
    spl_ASR <- smooth.spline(FSR_pert_ASR[,1] ~ rownames(FSR_pert_ASR))
    # estimate the slope of the tangent of the spline at the vital rate    
    ASR_pert_results[n+2, 2] <- predict(spl_ASR, x=FSR, deriv=1)$y
    # re-scale sensitivity into elasticity
    ASR_pert_results[n+2, 3] <- FSR/ASR * ASR_pert_results[n+3, 2]
    
    result <- list(ASR_pert_results = ASR_pert_results)
  }


horse_treat_sensitivity_analysis <- sensitivity_analysis(vital_rates = horse_treat, matrix_str = matrix_structure, h = horse_treat$h, k = horse_treat$k, FSR = horse_treat$FSR, niter = 1000, ASR = horse_treatment_ASR)

horse_mprime_sensitivity_analysis_male <- sensitivity_analysis(vital_rates = horse_mprime_male, matrix_str = matrix_structure, h = horse_mprime_male$h, k = horse_mprime_male$k, FSR = horse_mprime_male$FSR, niter = 1000, ASR = horse_mprime_ASR)

# Life table response experiment ------------------------------------------

LTRE_analysis <-
  function(Mprime_sens, vital_rates, sex){
    
    # make empty dataframes to stroe LTRE results for ASR
    LTRE_ASR <-
      data.frame(parameter = c("Foal survival", "Yearling survival", "2-year-old survival", "3-year-old survival",
                               "Adult survival", "Foal sex ratio",
                               "Mating system"),
                 contribution = numeric(7))
    
    # run a for loop to extract the parameter contributions
    if (sex == "male")
      # male rates scenario
      for(i in 1:nrow(LTRE_ASR))
      {
        LTRE_ASR[i, 2] <-
          ifelse(i < 6, (vital_rates[[i]] - vital_rates[[i + 5]]) * Mprime_sens$ASR_pert_results$sensitivities[i], #deviations in survival
                 ifelse(i == 6, (vital_rates[[13]] - 0.5) * Mprime_sens$ASR_pert_results$sensitivities[12], #deviations in offspring sex ratio
                        (vital_rates[[11]] - vital_rates[[11]]) * Mprime_sens$ASR_pert_results$sensitivities[11])) #and then deviations in the mating system
      }
    
    else
      # female rates scenario
      for(i in 1:nrow(LTRE_ASR))
      {
        LTRE_ASR[i, 2] <-
          ifelse(i < 6, (vital_rates[[i + 5]] - vital_rates[[i]]) *
                   Mprime_sens$ASR_pert_results$sensitivities[i + 5],
                 ifelse(i == 5, (vital_rates[[13]] - 0.5) *
                          Mprime_sens$ASR_pert_results$sensitivities[12],
                        (vital_rates[[11]] - vital_rates[[11]]) * 
                          Mprime_sens$ASR_pert_results$sensitivities[11]))
      }
    
    LTRE_ASR$parameter <- factor(LTRE_ASR$parameter, levels = c("Adult survival",
                                                                "2-year-old survival",
                                                                "3-year-old survival",
                                                                "Yearling survival",
                                                                "Foal survival",
                                                                "Foal sex ratio",
                                                                "Mating system"))
    
    LTRE_ASR$model <- "ASR"
    LTRE_df <- LTRE_ASR
    LTRE_df 
  }

horse_LTRE_male <- LTRE_analysis(Mprime_sens = horse_mprime_sensitivity_analysis_male, vital_rates = horse_treat, sex = "male")

horse_LTRE_male

LTRE_contributions_check <- 
  function(M_matrix_vital_rates, Mprime_sensitivities, M_matrix_ASR, scenario) {
    if (scenario == "male")
      contribution_sum <- 
        sum(
          (M_matrix_vital_rates[[1]] - M_matrix_vital_rates[[6]]) * 
            Mprime_sensitivities[1],
          (M_matrix_vital_rates[[2]] - M_matrix_vital_rates[[7]]) * 
            Mprime_sensitivities[2],
          (M_matrix_vital_rates[[3]] - M_matrix_vital_rates[[8]]) * 
            Mprime_sensitivities[3],
          (M_matrix_vital_rates[[4]] - M_matrix_vital_rates[[9]]) * 
            Mprime_sensitivities[4],
          (M_matrix_vital_rates[[5]] - M_matrix_vital_rates[[10]]) * 
            Mprime_sensitivities[5],
          (M_matrix_vital_rates[[6]] - M_matrix_vital_rates[[6]]) * 
            Mprime_sensitivities[6],
          (M_matrix_vital_rates[[7]] - M_matrix_vital_rates[[7]]) * 
            Mprime_sensitivities[7],
          (M_matrix_vital_rates[[8]] - M_matrix_vital_rates[[8]]) * 
            Mprime_sensitivities[8],
          (M_matrix_vital_rates[[9]] - M_matrix_vital_rates[[9]]) * 
            Mprime_sensitivities[9],
          (M_matrix_vital_rates[[10]] - M_matrix_vital_rates[[10]]) * 
            Mprime_sensitivities[10],
          (M_matrix_vital_rates[[13]] - 0.5) *
            Mprime_sensitivities[12],
          (M_matrix_vital_rates[[11]] - M_matrix_vital_rates[[11]]) * Mprime_sensitivities[11]
        )
    
    else
      contribution_sum <- 
        sum(
          (M_matrix_vital_rates[[5]] - M_matrix_vital_rates[[1]]) * 
            Mprime_sensitivities[5],
          (M_matrix_vital_rates[[6]] - M_matrix_vital_rates[[2]]) * 
            Mprime_sensitivities[6],
          (M_matrix_vital_rates[[7]] - M_matrix_vital_rates[[3]]) * 
            Mprime_sensitivities[7],
          (M_matrix_vital_rates[[8]] - M_matrix_vital_rates[[4]]) * 
            Mprime_sensitivities[8],
          (M_matrix_vital_rates[[1]] - M_matrix_vital_rates[[1]]) * 
            Mprime_sensitivities[1],
          (M_matrix_vital_rates[[2]] - M_matrix_vital_rates[[2]]) * 
            Mprime_sensitivities[2],
          (M_matrix_vital_rates[[3]] - M_matrix_vital_rates[[3]]) * 
            Mprime_sensitivities[3],
          (M_matrix_vital_rates[[4]] - M_matrix_vital_rates[[4]]) * 
            Mprime_sensitivities[4],
          (M_matrix_vital_rates[[13]] - 0.5) *
            Mprime_sensitivities[12],
          (M_matrix_vital_rates[[11]] - M_matrix_vital_rates[[11]]) * Mprime_sensitivities[11]
        )
    ASR_bias <- abs(M_matrix_ASR - 0.5)
    absolute_difference <- abs(ASR_bias) - abs(contribution_sum)
    
    return(list(contribution_sum = as.vector(contribution_sum), 
                ASR_bias = as.vector(ASR_bias), 
                absolute_difference = as.vector(absolute_difference)))
  }


LTRE_contributions_check(
  M_matrix_vital_rates = horse_treat, 
  M_matrix_ASR = horse_treatment_ASR, 
  Mprime_sensitivities = 
    horse_mprime_sensitivity_analysis_male$ASR_pert_results$sensitivities, 
  scenario = "male")

# Western subdivision -----------------------------------------------------

# creating matrices -------------------------------------------------------

horse_matrix <- 
  function(demographic_rates){
    # Define horse life stages
    stages <- c("F_foal", "F_1st_yr","F_2nd_yr", "F_3rd_yr", "F_Adt", "M_foal", "M_1st_yr","M_2nd_yr", "M_3rd_yr", "M_Adt")
    
    # Build the 10x10 matrix
    result <- 
      matrix(c(
        # top row
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        # second row
        demographic_rates$F_foal_surv,0, 0, 0, 0, 0, 0, 0, 0, 0,
        # third row
        0, demographic_rates$F_yearling_surv, 0, 0, 0, 0, 0, 0, 0, 0,
        # fourth row
        0, 0, demographic_rates$F_2yr_surv, 0, 0, 0, 0, 0, 0, 0,
        #fifth row
        0, 0, 0, demographic_rates$F_3yr_surv, demographic_rates$F_adult_surv, 0, 0, 0, 0, 0,
        # sixth row
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        # seventh row
        0, 0, 0, 0, 0, demographic_rates$M_foal_surv, 0, 0, 0, 0,
        # eighth row
        0, 0, 0, 0, 0, 0, demographic_rates$M_yearling_surv, 0, 0, 0,
        # ninth row
        0, 0, 0, 0, 0, 0, 0, demographic_rates$M_2yr_surv, 0, 0,
        #tenth row
        0, 0, 0, 0, 0, 0, 0, 0, demographic_rates$M_3yr_surv, demographic_rates$M_adult_surv),
        nrow = length(stages), byrow = TRUE,
        dimnames = list(stages, stages))
    result
  }

#constructing matrices

surv_rates_treat <- data.frame(F_foal = 0.793,
                               F_yearling = 0.891,
                               F_2yr = 0.807,
                               F_3yr = 0.961,
                               F_adult =0.870,
                               M_foal = 0.751,
                               M_yearling = 0.917,
                               M_2yr = 0.861,
                               M_3yr = 0.988,
                               M_adult = 0.920)


imm.rates_treat <- data.frame(F_foal = 0.025,
                              F_yearling = 0.014,
                              F_2yr = 0.028,
                              F_3yr = 0.025,
                              F_adult = 0.035,
                              M_foal = 0.023,
                              M_yearling = 0.042,
                              M_2yr = 0.027,
                              M_3yr = 0.088,
                              M_adult = 0.037)

emm.rates_treat <- data.frame(F_foal = 0.042,
                              F_yearling = 0.048,
                              F_2yr = 0.038,
                              F_3yr = 0.068,
                              F_adult = 0.040,
                              M_foal = 0.082,
                              M_yearling = 0.086,
                              M_2yr = 0.075,
                              M_3yr = 0.054,
                              M_adult = 0.038)

surv_rates_mprime <- data.frame(F_foal = (surv_rates_treat$F_foal + surv_rates_treat$M_foal)/2,
                                F_yearling = (surv_rates_treat$F_yearling + surv_rates_treat$M_yearling)/2,
                                F_2yr = (surv_rates_treat$F_2yr + surv_rates_treat$M_2yr)/2,
                                F_3yr = (surv_rates_treat$F_3yr + surv_rates_treat$M_3yr)/2,
                                F_adult = (surv_rates_treat$F_adult + surv_rates_treat$M_adult)/2,
                                M_foal = (surv_rates_treat$M_foal*2)/2,
                                M_yearling = (surv_rates_treat$M_yearling*2)/2,
                                M_2yr = (surv_rates_treat$M_2yr + surv_rates_treat$M_2yr)/2,
                                M_3yr = (surv_rates_treat$F_3yr + surv_rates_treat$M_3yr)/2,
                                M_adult = (surv_rates_treat$M_adult*2)/2)

imm.rates_mprime <- data.frame(F_foal = (imm.rates_treat$F_foal + imm.rates_treat$M_foal)/2,
                               F_yearling = (imm.rates_treat$F_yearling + imm.rates_treat$M_yearling)/2,
                               F_2yr = (imm.rates_treat$F_2yr + imm.rates_treat$M_2yr)/2,
                               F_3yr = (imm.rates_treat$F_3yr + imm.rates_treat$M_3yr)/2,
                               F_adult = (imm.rates_treat$F_adult + imm.rates_treat$M_adult)/2,
                               M_foal = (imm.rates_treat$M_foal*2)/2,
                               M_yearling = (imm.rates_treat$M_yearling*2)/2,
                               M_2yr = (imm.rates_treat$M_2yr + imm.rates_treat$M_2yr)/2,
                               M_3yr = (imm.rates_treat$M_3yr + imm.rates_treat$M_3yr)/2,
                               M_adult = (imm.rates_treat$M_adult*2)/2)

emm.rates_mprime <- data.frame(F_foal = (emm.rates_treat$F_foal + emm.rates_treat$M_foal)/2,
                               F_yearling = (emm.rates_treat$F_yearling + emm.rates_treat$M_yearling)/2,
                               F_2yr = (emm.rates_treat$F_2yr + emm.rates_treat$M_2yr)/2,
                               F_3yr = (emm.rates_treat$F_3yr + emm.rates_treat$M_3yr)/2,
                               F_adult = (emm.rates_treat$F_adult + emm.rates_treat$M_adult)/2,
                               M_foal = (emm.rates_treat$M_foal*2)/2,
                               M_yearling = (emm.rates_treat$M_yearling*2)/2,
                               M_2yr = (emm.rates_treat$M_2yr + emm.rates_treat$M_2yr)/2,
                               M_3yr = (emm.rates_treat$M_3yr + emm.rates_treat$M_3yr)/2,
                               M_adult = (emm.rates_treat$M_adult*2)/2)

horse_treat <- list(F_foal_surv = surv_rates_treat$F_foal*(imm.rates_treat$F_foal + (1-emm.rates_treat$F_foal)),
                    F_yearling_surv = surv_rates_treat$F_yearling*(imm.rates_treat$F_yearling + (1-emm.rates_treat$F_yearling)),
                    F_2yr_surv = surv_rates_treat$F_2yr*(imm.rates_treat$F_2yr + (1-emm.rates_treat$F_2yr)),
                    F_3yr_surv = surv_rates_treat$F_3yr*(imm.rates_treat$F_3yr + (1-emm.rates_treat$F_3yr)),
                    F_adult_surv = surv_rates_treat$F_adult*(imm.rates_treat$F_adult + (1-emm.rates_treat$F_adult)),
                    M_foal_surv = surv_rates_treat$M_foal*(imm.rates_treat$M_foal + (1-emm.rates_treat$M_foal)),
                    M_yearling_surv = surv_rates_treat$M_yearling*(imm.rates_treat$M_yearling + (1-emm.rates_treat$M_yearling)),
                    M_2yr_surv = surv_rates_treat$M_2yr*(imm.rates_treat$M_2yr + (1-emm.rates_treat$M_2yr)),
                    M_3yr_surv = surv_rates_treat$M_3yr*(imm.rates_treat$M_3yr + (1-emm.rates_treat$M_3yr)), 
                    M_adult_surv = surv_rates_treat$M_adult*(imm.rates_treat$M_adult + (1-emm.rates_treat$M_adult)),
                    h = 1.818,
                    k = 1,
                    FSR = 0.5)

horse_mprime_male <- list(F_foal_surv = ((surv_rates_treat$F_foal*(imm.rates_treat$F_foal + (1-emm.rates_treat$F_foal))) + (surv_rates_treat$M_foal*(imm.rates_treat$M_foal + (1-emm.rates_treat$M_foal))))/2,
                          F_yearling_surv = ((surv_rates_treat$F_yearling*(imm.rates_treat$F_yearling + (1-emm.rates_treat$F_yearling))) + (surv_rates_treat$M_yearling*(imm.rates_treat$M_yearling + (1-emm.rates_treat$M_yearling))))/2,
                          F_2yr_surv = ((surv_rates_treat$F_2yr*(imm.rates_treat$F_2yr + (1-emm.rates_treat$F_2yr))) + (surv_rates_treat$M_2yr*(imm.rates_treat$M_2yr + (1-emm.rates_treat$M_2yr))))/2,
                          F_3yr_surv = ((surv_rates_treat$F_3yr*(imm.rates_treat$F_3yr + (1-emm.rates_treat$F_3yr))) + (surv_rates_treat$M_3yr*(imm.rates_treat$M_3yr + (1-emm.rates_treat$M_3yr))))/2,
                          F_adult_surv = ((surv_rates_treat$F_adult*(imm.rates_treat$F_adult + (1-emm.rates_treat$F_adult))) + (surv_rates_treat$M_adult*(imm.rates_treat$M_adult + (1-emm.rates_treat$M_adult))))/2,
                          M_foal_surv = ((surv_rates_treat$M_foal*(imm.rates_treat$M_foal + (1-emm.rates_treat$M_foal))) + (surv_rates_treat$M_foal*(imm.rates_treat$M_foal + (1-emm.rates_treat$M_foal))))/2,
                          M_yearling_surv = ((surv_rates_treat$M_yearling*(imm.rates_treat$M_yearling + (1-emm.rates_treat$M_yearling))) + (surv_rates_treat$M_yearling*(imm.rates_treat$M_yearling + (1-emm.rates_treat$M_yearling))))/2,
                          M_2yr_surv = ((surv_rates_treat$M_2yr*(imm.rates_treat$M_2yr + (1-emm.rates_treat$M_2yr))) + (surv_rates_treat$M_2yr*(imm.rates_treat$M_2yr + (1-emm.rates_treat$M_2yr))))/2,
                          M_3yr_surv = ((surv_rates_treat$M_3yr*(imm.rates_treat$M_3yr + (1-emm.rates_treat$M_3yr))) + (surv_rates_treat$M_3yr*(imm.rates_treat$M_3yr + (1-emm.rates_treat$M_3yr))))/2,
                          M_adult_surv = ((surv_rates_treat$M_adult*(imm.rates_treat$M_adult + (1-emm.rates_treat$M_adult))) + (surv_rates_treat$M_adult*(imm.rates_treat$M_adult + (1-emm.rates_treat$M_adult))))/2,
                          h = (1.818*2)/2,
                          k = 1,
                          FSR = 0.5)

horse_treatment_matrix <- horse_matrix(horse_treat)

horse_mprime_matrix <- horse_matrix(horse_mprime_male)

# Estimating the ASR ------------------------------------------------------

matrix_ASR <-
  function(M, n = rep(10, nrow(M)), h = 1, k = 1, 
           iterations = 1000, FSR = 0.5, plot = FALSE){
    # Number of stages in matrix
    x <- length(n)
    # Number of time steps to simulate
    t <- iterations
    # an empty t by x matrix to store the stage distributions
    stage <- matrix(numeric(x * t), nrow = x)
    # an empty t vector to store the population sizes
    pop <- numeric(t)
    # for loop that goes through each of t time steps
    for (i in 1:t) {
      # stage distribution at time t
      stage[,i] <- n
      # population size at time t
      pop[i] <- sum(n)
      # number of male adults at time t
      M2 <- stage[10, i]
      # number of female adults at time t
      F2 <- stage[5, i]
      # Female freq-dep fecundity of female foals
      M[1,x/2]        <- min(0.5*0.60, k*M2/(M2 + ((0.60*F2)/h)))*(1-FSR)
      # Female freq-dep fecundity of male foals
      M[(x/10)*6,x/2]  <- min(0.5*0.60, k*M2/(M2 + ((0.60*F2)/h)))*FSR
      # Male freq-dep fecundity of female foals
      M[1,x]          <- min(0.5*(0.60*F2)/M2, k*(0.60*F2)/(M2 + ((0.60*F2)/h)))*(1-FSR)
      # Male freq-dep fecundity of male foals
      M[(x/10)*6,x]    <- min(0.5*(0.60*F2)/M2, k*(0.60*F2)/(M2 + ((0.60*F2)/h)))*FSR
      
      # define the new n (i.e., new stage distribution at time t)
      n <- M %*% n
      # define rownames of stage matrix
      rownames(stage) <- rownames(M)
      # define colnames of stage matrix
      colnames(stage) <- 0:(t - 1)
      # calculate the proportional stable stage distribution - starts completely even as we start with 10 inds in each
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, t]
    }
    # calc ASR as the proportion of the adult stable stage class that is male
    ASR <- stable.stage[x]/(stable.stage[x/2] + stable.stage[x])
    
    if(plot)
    {
      # plot distrubution to assure that it is not chaotic
      matplot(rownames(t(stage)), t(stage), type='l', lwd=2, las=1)
    }
    # make a list of results
    pop.proj <- list(ASR = ASR,
                     stable.stage = stable.stage,
                     stage.vectors = stage,
                     SSD_M2 = stable.stage[10],
                     SSD_F2 = stable.stage[5])
    # print the list as output to the function
    pop.proj
  }

#treatment first
horse_treatment_ASR_analysis <- 
  matrix_ASR(M = horse_treatment_matrix, h = horse_treat$h, FSR = horse_treat$FSR, 
             iterations = 1000, plot = T)
horse_treatment_ASR <- horse_treatment_ASR_analysis$ASR
horse_treatment_ASR

#now the mprime matrix
horse_mprime_ASR_analysis <- 
  matrix_ASR(M = horse_mprime_matrix, h = horse_mprime_male$h, FSR = horse_mprime_male$FSR, 
             iterations = 1000, plot = T)
horse_mprime_ASR <- horse_mprime_ASR_analysis$ASR
horse_mprime_ASR

# Sensitivity analysis ----------------------------------------------------

matrix_structure <- expression(
  # top row of matrix
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  # second row of matrix
  F_foal, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  # third row of matrix
  0, F_1st_yr, 0, 0, 0, 0, 0, 0, 0, 0,
  # fourth row of matrix
  0, 0, F_2yr, 0, 0, 0, 0, 0, 0, 0,
  # fifth row
  0, 0, 0, F_3yr, F_Adt, 0, 0, 0, 0, 0,
  # sixth row
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  # seventh row
  0, 0, 0, 0, 0, M_foal, 0, 0, 0, 0,
  # eighth row
  0, 0, 0, 0, 0, 0, M_1st_yr, 0, 0, 0,
  #ninth row
  0, 0, 0, 0, 0, 0, 0, M_2yr, 0, 0,
  # tenth row
  0, 0, 0, 0, 0, 0, 0, 0, M_3yr, M_Adt)

sensitivity_analysis <-
  function(vital_rates, matrix_str, h, k, FSR, niter = 1000, ASR, imm_rates, emm_rates){
    
    # make a list of all parameters
    vr <- list(F_foal = vital_rates$F_foal,
               F_1st_yr = vital_rates$F_yearling,
               F_2yr = vital_rates$F_2yr,
               F_3yr = vital_rates$F_3yr,
               F_Adt = vital_rates$F_adult,
               M_foal = vital_rates$M_foal,
               M_1st_yr = vital_rates$M_yearling,
               M_2yr = vital_rates$M_2yr,
               M_3yr = vital_rates$M_3yr,
               M_Adt = vital_rates$M_adult)
    
    emm <- list(F_foal = emm_rates$F_foal,
                F_1st_yr = emm_rates$F_yearling,
                F_2yr = emm_rates$F_2yr,
                F_3yr = emm_rates$F_3yr,
                F_Adt = emm_rates$F_adult,
                M_foal = emm_rates$M_foal,
                M_1st_yr = emm_rates$M_yearling,
                M_2yr = emm_rates$M_yearling,
                M_3yr = emm_rates$M_yearling,
                M_Adt = emm_rates$M_adult)
    
    imm <- list(F_foal = imm_rates$F_foal,
                F_1st_yr = imm_rates$F_yearling,
                F_2yr = imm_rates$F_2yr,
                F_3yr = imm_rates$F_3yr,
                F_Adt = imm_rates$F_adult,
                M_foal = imm_rates$M_foal,
                M_1st_yr = imm_rates$M_yearling,
                M_2yr = imm_rates$M_yearling,
                M_3yr = imm_rates$M_yearling,
                M_Adt = imm_rates$M_adult)
    # number of stages in the matrix
    no_stages <- sqrt(length(matrix_str))
    
    # Define horse life stages
    stages <- c("F_foal", "F_1st_yr", "F_2yr", "F_3yr", "F_Adt", "M_foal", "M_1st_yr",  "M_2yr", "M_3yr", "M_Adt")
    
    # an empty t by x matrix
    stage <- matrix(numeric(no_stages * niter), nrow = no_stages)
    
    # an empty t vector to store the population sizes
    pop <- numeric(niter)
    
    # dataframe to store the perturbation results
    ASR_pert_results <-
      data.frame(parameter = c("F_foal_surv", "F_1st_yr_surv", "F_2yr_surv","F_3yr_surv", "F_Adt_surv",
                               "M_foal_surv", "M_1st_yr_surv", "M_2yr_surv", "M_3yr_surv", "M_Adt_surv",
                               "F_foal_imm", "F_1st_yr_imm", "F_2yr_imm","F_3yr_imm", "F_Adt_imm",
                               "M_foal_imm", "M_1st_yr_imm", "M_2yr_imm","M_3yr_imm", "M_Adt_imm", "F_foal_emm", "F_1st_yr_emm", "F_2yr_emm", "F_3yr_emm", "F_Adt_emm",
                               "M_foal_emm", "M_1st_yr_emm", "M_2yr_emm", "M_3yr_emm", "M_Adt_emm","h", "FSR"),
                 sensitivities = numeric(32),
                 elasticities = numeric(32))
    
    # specifiy how many survival rates there are
    n <- length(vr)
    
    # create vectors of perturbations to test on parameters of the matrix model
    vr_nums <- seq(0, 1, 0.01) # changes in survival and FSR (i.e., between 0 and 1)
    h_nums <- seq(0, 2, 0.02) # changes in h index (i.e., between 0 and 2)
    emm_nums <- seq(0, 1, 0.01) #changes in emigration
    imm_nums <- seq(0, 1, 0.01) #changes in immigration
    
    # create empty dataframes to store the perturbation results for ASR
    vr_pert_ASR <- matrix(numeric(n * length(vr_nums)),
                          ncol = n, dimnames = list(vr_nums, names(vr)))
    h_pert_ASR <- matrix(numeric(length(h_nums)),
                         ncol = 1, dimnames = list(h_nums, "h"))
    # k_pert_ASR <- matrix(numeric(length(k_nums)),
    #                      ncol = 1, dimnames = list(k_nums, "k"))
    FSR_pert_ASR <- matrix(numeric(length(vr_nums)),
                           ncol = 1, dimnames = list(vr_nums, "FSR"))
    emm_pert_ASR <- matrix(numeric(n * length(emm_nums)),
                           ncol = n, dimnames = list(emm_nums, names(emm)))
    imm_pert_ASR <- matrix(numeric(n * length(imm_nums)),
                           ncol = n, dimnames = list(imm_nums, names(imm)))
    # perturbation of vital rates survival rates
    for (g in 1:n) # pick a column (i.e., a variable)
    {
      vr2 <- vr # reset the vital rates to the original
      for (i in 1:length(vr_nums)) # pick a perturbation level
      {
        vr2[[g]] <- vr_nums[i] # specify the vital rate with the new perturbation level
        A <- matrix(sapply(matrix_str, eval, vr2, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        I <- matrix(sapply(matrix_str, eval, imm, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        E <- matrix(sapply(matrix_str, eval, emm, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        Mov.M <- I + (1-E)
        # reset the starting stage distribution for simulation (all with 10 individuals)
        A <- A*Mov.M
        m <- rep(10, no_stages) 
        for (j in 1:niter) { # project the matrix through t iteration
          # stage distribution at time t
          stage[,j] <- m
          # population size at time t
          pop[j] <- sum(m)
          
          # number of male adults at time t
          M2 <- stage[10, j]
          # number of female adults at time t
          F2 <- stage[5, j]
          # Female freq-dep fecundity of female foals
          A[1,no_stages/2]        <- min((0.5*0.60), k*M2/(M2 + ((F2*0.60)/h)))*(1-FSR)
          # Female freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.60), k*M2/(M2 + ((F2*0.60)/h)))*FSR
          # Male freq-dep fecundity of female foals
          A[1,no_stages]          <- min(0.5*(0.60*F2)/M2, k*(0.60*F2)/(M2 + ((0.60*F2)/h)))*(1-FSR)
          # Male freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.60*F2)/M2, k*(0.60*F2)/(M2 + ((0.60*F2)/h)))*FSR
          # define the new n (i.e., new stage distribution at time t)
          m <- A %*% m
        }
        # define rownames of stage matrix
        rownames(stage) <- rownames(A)
        # define colnames of stage matrix
        colnames(stage) <- 0:(niter - 1)
        # calculate the proportional stable stage distribution
        stage <- apply(stage, 2, function(x) x/sum(x))
        # define stable stage as the last stage
        stable.stage <- stage[, niter]
        # calc ASR as the proportion of the adult stable stage class that is male
        vr_pert_ASR[i, g] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                        stable.stage[no_stages])
      }
      # get the spline function of ASR
      spl_ASR <- smooth.spline(vr_pert_ASR[,g] ~ rownames(vr_pert_ASR))
      # estimate the slope of the tangent of the spline at the vital rate
      ASR_pert_results[g, 2] <- predict(spl_ASR, x=vr[[g]], deriv=1)$y
      # re-scale sensitivity into elasticity
      ASR_pert_results[g, 3] <- vr[[g]]/ASR * ASR_pert_results[g, 2]
    }
    # perturbation of the h index parameter
    for (i in (1:length(h_nums))) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                  nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                  dimnames = list(stages, stages)) # build the matrix with the new value
      # reset the starting stage distribution for simulation (all with 10 individuals)
      
      A <- A*Mov.M
      
      m <- rep(10, no_stages) 
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        # population size at time t
        pop[j] <- sum(m)
        
        # number of male adults at time t
        M2 <- stage[10, j]
        # number of female adults at time t
        F2 <- stage[5, j]
        # Female freq-dep fecundity of female foals
        A[1,no_stages/2]        <- min((0.5*0.60), k*M2/(M2 + ((F2*0.60)/h_nums[i])))*(1-FSR)
        # Female freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.60), k*M2/(M2 + ((F2*0.60)/h_nums[i])))*FSR
        # Male freq-dep fecundity of female foals
        A[1,no_stages]          <- min(0.5*(0.60*F2)/M2, k*(0.60*F2)/(M2 + ((0.60*F2)/h_nums[i])))*(1-FSR)
        # Male freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.60*F2)/M2, k*(0.60*F2)/(M2 + ((0.60*F2)/h_nums[i])))*FSR
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% m
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      # calc ASR as the proportion of the adult stable stage class that is male
      h_pert_ASR[i,] <- 
        stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
    }
    # get the spline function of ASR
    spl_ASR <- smooth.spline(h_pert_ASR[, 1] ~ rownames(h_pert_ASR))
    # estimate the slope of the tangent of the spline at the vital rate
    ASR_pert_results[(n*3)+1, 2] <- predict(spl_ASR, x=h, deriv=1)$y
    # re-scale sensitivity into elasticity
    ASR_pert_results[(n*3)+1, 3] <- h/ASR * ASR_pert_results[(n*3)+1, 2]
    
    # perturbation of FSR
    for (i in 1:length(vr_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                  nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                  dimnames = list(stages, stages)) # build the matrix with the new value
      
      A <- A*Mov.M
      
      # reset the starting stage distribution for simulation (all with 10 individuals)
      m <- rep(10, no_stages) 
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        # population size at time t
        pop[j] <- sum(m)
        
        # number of male adults at time t
        M2 <- stage[10, j]
        # number of female adults at time t
        F2 <- stage[5, j]
        # Female freq-dep fecundity of female foals
        A[1,no_stages/2]        <- min((0.5*0.60), k*M2/(M2 + ((F2*0.60)/h)))*(1-vr_nums[i])
        # Female freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.60), k*M2/(M2 + ((F2*0.60)/h)))*vr_nums[i]
        # Male freq-dep fecundity of female foals
        A[1,no_stages]          <- min(0.5*(0.60*F2)/M2, k*(0.60*F2)/(M2 + ((0.60*F2)/h)))*(1-vr_nums[i])
        # Male freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.60*F2)/M2, k*(0.60*F2)/(M2 + ((0.60*F2)/h)))*vr_nums[i]
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% m
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      # calc ASR as the proportion of the adult stable stage class that is male
      FSR_pert_ASR[i,] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                     stable.stage[no_stages])
    }
    # get the spline function of ASR
    spl_ASR <- smooth.spline(FSR_pert_ASR[,1] ~ rownames(FSR_pert_ASR))
    # estimate the slope of the tangent of the spline at the vital rate    
    ASR_pert_results[(n*3)+2, 2] <- predict(spl_ASR, x=FSR, deriv=1)$y
    # re-scale sensitivity into elasticity
    ASR_pert_results[(n*3)+2, 3] <- FSR/ASR * ASR_pert_results[(n*3)+2, 2]
    
    #Now perturbing immigration
    for (g in 1:n) # pick a column (i.e., a variable)
    {
      imm2 <- imm # reset the vital rates to the original
      for (i in 1:length(imm_nums)) # pick a perturbation level
      {
        imm2[[g]] <- imm_nums[i] # specify the vital rate with the new perturbation level
        A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        I <- matrix(sapply(matrix_str, eval, imm2, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        E <- matrix(sapply(matrix_str, eval, emm, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        # reset the starting stage distribution for simulation (all with 10 individuals)
        Mov.M <- I + (1-E)
        
        A <- A*Mov.M
        
        m <- rep(10, no_stages) 
        for (j in 1:niter) { # project the matrix through t iteration
          # stage distribution at time t
          stage[,j] <- m
          # population size at time t
          pop[j] <- sum(m)
          
          # number of male adults at time t
          M2 <- stage[10, j]
          # number of female adults at time t
          F2 <- stage[5, j]
          # Female freq-dep fecundity of female foals
          A[1,no_stages/2]        <- min((0.5*0.60), k*M2/(M2 + ((F2*0.60)/h)))*(1-FSR)
          # Female freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.60), k*M2/(M2 + ((F2*0.60)/h)))*FSR
          # Male freq-dep fecundity of female foals
          A[1,no_stages]          <- min(0.5*(0.60*F2)/M2, k*(0.60*F2)/(M2 + ((0.60*F2)/h)))*(1-FSR)
          # Male freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.60*F2)/M2, k*(0.60*F2)/(M2 + ((0.60*F2)/h)))*FSR
          # define the new n (i.e., new stage distribution at time t)
          m <- A %*% m
        }
        # define rownames of stage matrix
        rownames(stage) <- rownames(A)
        # define colnames of stage matrix
        colnames(stage) <- 0:(niter - 1)
        # calculate the proportional stable stage distribution
        stage <- apply(stage, 2, function(x) x/sum(x))
        # define stable stage as the last stage
        stable.stage <- stage[, niter]
        # calc ASR as the proportion of the adult stable stage class that is male
        imm_pert_ASR[i, g] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                         stable.stage[no_stages])
      }
      # get the spline function of ASR
      spl_ASR <- smooth.spline(imm_pert_ASR[,g] ~ rownames(imm_pert_ASR))
      # estimate the slope of the tangent of the spline at the vital rate
      ASR_pert_results[n+g, 2] <- predict(spl_ASR, x=imm[[g]], deriv=1)$y
      # re-scale sensitivity into elasticity
      ASR_pert_results[n+g, 3] <- imm[[g]]/ASR * ASR_pert_results[n+g, 2]
    }
    
    #Now perturbing emmigration
    for (g in 1:n) # pick a column (i.e., a variable)
    {
      emm2 <- emm # reset the vital rates to the original
      for (i in 1:length(emm_nums)) # pick a perturbation level
      {
        emm2[[g]] <- emm_nums[i] # specify the vital rate with the new perturbation level
        A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        I <- matrix(sapply(matrix_str, eval, imm, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        E <- matrix(sapply(matrix_str, eval, emm2, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        # reset the starting stage distribution for simulation (all with 10 individuals)
        Mov.M <- I + (1-E)
        
        A <- A*Mov.M
        
        m <- rep(10, no_stages) 
        for (j in 1:niter) { # project the matrix through t iteration
          # stage distribution at time t
          stage[,j] <- m
          # population size at time t
          pop[j] <- sum(m)
          
          # number of male adults at time t
          M2 <- stage[10, j]
          # number of female adults at time t
          F2 <- stage[5, j]
          # Female freq-dep fecundity of female foals
          A[1,no_stages/2]        <- min((0.5*0.60), k*M2/(M2 + ((F2*0.60)/h)))*(1-FSR)
          # Female freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.60), k*M2/(M2 + ((F2*0.60)/h)))*FSR
          # Male freq-dep fecundity of female foals
          A[1,no_stages]          <- min(0.5*(0.60*F2)/M2, k*(0.60*F2)/(M2 + ((0.60*F2)/h)))*(1-FSR)
          # Male freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.60*F2)/M2, k*(0.60*F2)/(M2 + ((0.60*F2)/h)))*FSR
          # define the new n (i.e., new stage distribution at time t)
          m <- A %*% m
        }
        # define rownames of stage matrix
        rownames(stage) <- rownames(A)
        # define colnames of stage matrix
        colnames(stage) <- 0:(niter - 1)
        # calculate the proportional stable stage distribution
        stage <- apply(stage, 2, function(x) x/sum(x))
        # define stable stage as the last stage
        stable.stage <- stage[, niter]
        # calc ASR as the proportion of the adult stable stage class that is male
        emm_pert_ASR[i, g] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                         stable.stage[no_stages])
      }
      # get the spline function of ASR
      spl_ASR <- smooth.spline(emm_pert_ASR[,g] ~ rownames(emm_pert_ASR))
      # estimate the slope of the tangent of the spline at the vital rate
      ASR_pert_results[(n*2)+g, 2] <- predict(spl_ASR, x=emm[[g]], deriv=1)$y
      # re-scale sensitivity into elasticity
      ASR_pert_results[(n*2)+g, 3] <- emm[[g]]/ASR * ASR_pert_results[(n*2)+g, 2]
    }
    
    result <- list(ASR_pert_results = ASR_pert_results)
  }

horse_treat_sensitivity_analysis <- sensitivity_analysis(vital_rates = surv_rates_treat, matrix_str = matrix_structure, h = horse_treat$h, k = horse_treat$k, FSR = horse_treat$FSR, niter = 1000, ASR = horse_treatment_ASR, imm_rates = imm.rates_treat, emm_rates = emm.rates_treat)

horse_mprime_sensitivity_analysis_male <- sensitivity_analysis(vital_rates = surv_rates_mprime, matrix_str = matrix_structure, h = horse_mprime_male$h, k = horse_mprime_male$k, FSR = horse_mprime_male$FSR, niter = 1000, ASR = horse_mprime_ASR, imm_rates = imm.rates_mprime, emm_rates = emm.rates_mprime)

# Life table response experiment ------------------------------------------

LTRE_analysis <-
  function(Mprime_sens, surv_rates, imm_rates, emm_rates, sex, treat_h, treat_FSR){
    
    # make empty dataframes to store LTRE results for ASR
    LTRE_ASR <-
      data.frame(parameter = c("Foal survival", "Yearling survival", "2-year-old survival", "3-year-old survival",
                               "Adult survival", "Foal immigration", "Yearling immigration", "2-year-old immigration", "3-year-old immigration", "Adult immigration", "Foal emmigration", "Yearling emmigration", "2-year-old emmigration", "3-year-old emmigration", "Adult emmigration","Foal sex ratio",
                               "Mating system"),
                 contribution = numeric(17))
    
    # run a for loop to extract the parameter contributions
    if (sex == "male")
      # male rates scenario
      for(i in 1:nrow(LTRE_ASR))
      {
        LTRE_ASR[i, 2] <-
          ifelse(i < 6, (surv_rates[[i]] - surv_rates[[i + 5]]) * Mprime_sens$ASR_pert_results$sensitivities[i], #deviations in survival
                 ifelse(i >5 & i <11, (imm_rates[[i - 5]] - imm_rates[[i]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 5], #deviations in immigration
                        ifelse(i > 10 & i < 16, (emm_rates[[i - 10]] - emm_rates[[i - 5]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 10], #deviations in emigration
                               ifelse(i == 16, (treat_FSR - 0.5) * Mprime_sens$ASR_pert_results$sensitivities[i + 16], #then looking at the effect of deviations in offspring sex ratio
                                      (treat_h - treat_h) * Mprime_sens$ASR_pert_results$sensitivities[i + 14])))) #and then deviations in the mating system
      }
    
    else
      # female rates scenario
      for(i in 1:nrow(LTRE_ASR))
      {
        LTRE_ASR[i, 2] <-
          ifelse(i < 6, (surv_rates[[i + 5]] - surv_rates[[i]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 5],
                 ifelse(i > 5 & i < 11, (imm_rates[[i]] - imm_rates[[i - 5]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 19],
                        ifelse(i > 10 & i < 16, (emm_rates[[i - 5]] - emm_rates[[i - 10]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 15],
                               ifelse(i == 16, (treat_FSR - 0.5) * Mprime_sens$ASR_pert_results$sensitivities[i + 16], 
                                      treat_h - treat_h * Mprime_sens$ASR_pert_results$sensitivities[i + 14]))))
      }
    
    LTRE_ASR$parameter <- factor(LTRE_ASR$parameter, levels = c("Adult survival",
                                                                "3-year-old survival",
                                                                "2-year-old survival",
                                                                "Yearling survival",
                                                                "Foal survival",
                                                                "Adult immigration",
                                                                "3-year-old immigration",
                                                                "2-year-old immigration",
                                                                "Yearling immigration",
                                                                "Foal immigration",
                                                                "Adult emmigration",
                                                                "3-year-old emmigration",
                                                                "2-year-old emmigration",
                                                                "Yearling emmigration",
                                                                "Foal emmigration",
                                                                "Foal sex ratio",
                                                                "Mating system"))
    
    LTRE_ASR$model <- "ASR"
    LTRE_df <- LTRE_ASR
    LTRE_df 
  }


LTRE_west <- LTRE_analysis(Mprime_sens = horse_mprime_sensitivity_analysis_male, surv_rates = surv_rates_treat, imm_rates = imm.rates_treat, emm_rates = emm.rates_treat, treat_FSR = horse_treat$FSR, treat_h = horse_treat$h, sex = "male")

LTRE_west

LTRE_contributions_check <- 
  function(surv_rates, imm_rates, emm_rates, treat_FSR, treat_h, Mprime_sensitivities, M_matrix_ASR, scenario) {
    if (scenario == "male")
      contribution_sum <- 
        sum(
          (surv_rates[[1]] - surv_rates[[6]]) * 
            Mprime_sensitivities[1],
          (surv_rates[[2]] - surv_rates[[7]]) * 
            Mprime_sensitivities[2],
          (surv_rates[[3]] - surv_rates[[8]]) * 
            Mprime_sensitivities[3],
          (surv_rates[[4]] - surv_rates[[9]]) * 
            Mprime_sensitivities[4],
          (surv_rates[[5]] - surv_rates[[10]]) * 
            Mprime_sensitivities[5],
          (surv_rates[[6]] - surv_rates[[6]]) * 
            Mprime_sensitivities[6],
          (surv_rates[[7]] - surv_rates[[7]]) * 
            Mprime_sensitivities[7],
          (surv_rates[[8]] - surv_rates[[8]]) * 
            Mprime_sensitivities[8],
          (surv_rates[[9]] - surv_rates[[9]]) * 
            Mprime_sensitivities[9],
          (surv_rates[[10]] - surv_rates[[10]]) * 
            Mprime_sensitivities[10],
          (imm_rates[[1]] - imm_rates[[6]]) * 
            Mprime_sensitivities[11],
          (imm_rates[[2]] - imm_rates[[7]]) * 
            Mprime_sensitivities[12],
          (imm_rates[[3]] - imm_rates[[8]]) * 
            Mprime_sensitivities[13],
          (imm_rates[[4]] - imm_rates[[9]]) * 
            Mprime_sensitivities[14],
          (imm_rates[[5]] - imm_rates[[10]]) * 
            Mprime_sensitivities[15],
          (imm_rates[[6]] - imm_rates[[6]]) * 
            Mprime_sensitivities[16],
          (imm_rates[[7]] - imm_rates[[7]]) * 
            Mprime_sensitivities[17],
          (imm_rates[[8]] - imm_rates[[8]]) * 
            Mprime_sensitivities[18],
          (imm_rates[[9]] - imm_rates[[10]]) * 
            Mprime_sensitivities[19],
          (imm_rates[[10]] - imm_rates[[10]]) * 
            Mprime_sensitivities[20],
          (emm_rates[[1]] - emm_rates[[6]]) * 
            Mprime_sensitivities[21],
          (emm_rates[[2]] - emm_rates[[7]]) * 
            Mprime_sensitivities[22],
          (emm_rates[[3]] - emm_rates[[8]]) * 
            Mprime_sensitivities[23],
          (emm_rates[[4]] - emm_rates[[9]]) * 
            Mprime_sensitivities[24],
          (emm_rates[[5]] - emm_rates[[10]]) * 
            Mprime_sensitivities[25],
          (emm_rates[[6]] - emm_rates[[6]]) * 
            Mprime_sensitivities[26],
          (emm_rates[[7]] - emm_rates[[7]]) * 
            Mprime_sensitivities[27],
          (emm_rates[[8]] - emm_rates[[8]]) * 
            Mprime_sensitivities[28],
          (emm_rates[[9]] - emm_rates[[9]]) * 
            Mprime_sensitivities[29],
          (emm_rates[[10]] - emm_rates[[10]]) * 
            Mprime_sensitivities[30],
          (treat_h - treat_h) *
            Mprime_sensitivities[31],
          (treat_FSR - 0.5) * Mprime_sensitivities[32]
        )
    
    else
      contribution_sum <- 
        sum(
          (surv_rates[[6]] - surv_rates[[1]]) * 
            Mprime_sensitivities[6],
          (surv_rates[[7]] - surv_rates[[2]]) * 
            Mprime_sensitivities[7],
          (surv_rates[[8]] - surv_rates[[3]]) * 
            Mprime_sensitivities[8],
          (surv_rates[[9]] - surv_rates[[4]]) * 
            Mprime_sensitivities[9],
          (surv_rates[[10]] - surv_rates[[5]]) * 
            Mprime_sensitivities[10],
          (surv_rates[[1]] - surv_rates[[1]]) * 
            Mprime_sensitivities[1],
          (surv_rates[[2]] - surv_rates[[2]]) * 
            Mprime_sensitivities[2],
          (surv_rates[[3]] - surv_rates[[3]]) * 
            Mprime_sensitivities[3],
          (surv_rates[[4]] - surv_rates[[4]]) * 
            Mprime_sensitivities[4],
          (surv_rates[[5]] - surv_rates[[5]]) * 
            Mprime_sensitivities[5],
          
          (imm_rates[[6]] - imm_rates[[1]]) * 
            Mprime_sensitivities[16],
          (imm_rates[[7]] - imm_rates[[2]]) * 
            Mprime_sensitivities[17],
          (imm_rates[[8]] - imm_rates[[3]]) * 
            Mprime_sensitivities[18],
          (imm_rates[[9]] - imm_rates[[4]]) * 
            Mprime_sensitivities[19],
          (imm_rates[[10]] - imm_rates[[1]]) * 
            Mprime_sensitivities[20],
          (imm_rates[[1]] - imm_rates[[1]]) * 
            Mprime_sensitivities[11],
          (imm_rates[[2]] - imm_rates[[2]]) * 
            Mprime_sensitivities[12],
          (imm_rates[[3]] - imm_rates[[3]]) * 
            Mprime_sensitivities[13],
          (imm_rates[[4]] - imm_rates[[4]]) * 
            Mprime_sensitivities[14],
          (imm_rates[[5]] - imm_rates[[5]]) * 
            Mprime_sensitivities[15],
          
          (emm_rates[[6]] - emm_rates[[1]]) * 
            Mprime_sensitivities[26],
          (emm_rates[[7]] - emm_rates[[2]]) * 
            Mprime_sensitivities[27],
          (emm_rates[[8]] - emm_rates[[3]]) * 
            Mprime_sensitivities[28],
          (emm_rates[[9]] - emm_rates[[4]]) * 
            Mprime_sensitivities[29],
          (emm_rates[[10]] - emm_rates[[5]]) * 
            Mprime_sensitivities[30],
          (emm_rates[[1]] - emm_rates[[1]]) * 
            Mprime_sensitivities[21],
          (emm_rates[[2]] - emm_rates[[2]]) * 
            Mprime_sensitivities[22],
          (emm_rates[[3]] - emm_rates[[3]]) * 
            Mprime_sensitivities[23],
          (emm_rates[[4]] - emm_rates[[4]]) * 
            Mprime_sensitivities[24],
          (emm_rates[[5]] - emm_rates[[5]]) * 
            Mprime_sensitivities[25],
          (treat_h - treat_h) *
            Mprime_sensitivities[31],
          (treat_FSR - 0.5) * Mprime_sensitivities[32]
        )
    
    ASR_bias <- abs(M_matrix_ASR - 0.5)
    absolute_difference <- abs(ASR_bias) - abs(contribution_sum)
    
    return(list(contribution_sum = as.vector(contribution_sum), 
                ASR_bias = as.vector(ASR_bias), 
                absolute_difference = as.vector(absolute_difference)))
  }


LTRE_contributions_check(surv_rates = surv_rates_treat, imm_rates = imm.rates_treat, emm_rates = emm.rates_treat, treat_h = horse_treat$h, treat_FSR = horse_treat$FSR, M_matrix_ASR = horse_treatment_ASR, Mprime_sensitivities = horse_mprime_sensitivity_analysis_male$ASR_pert_results$sensitivities, scenario = "male")

# Central subdivision -----------------------------------------------------

# creating matrices -------------------------------------------------------

horse_matrix <- 
  function(demographic_rates){
    # Define horse life stages
    stages <- c("F_foal", "F_1st_yr","F_2nd_yr", "F_3rd_yr", "F_Adt", "M_foal", "M_1st_yr","M_2nd_yr", "M_3rd_yr", "M_Adt")
    
    # Build the 10x10 matrix
    result <- 
      matrix(c(
        # top row
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        # second row
        demographic_rates$F_foal_surv,0, 0, 0, 0, 0, 0, 0, 0, 0,
        # third row
        0, demographic_rates$F_yearling_surv, 0, 0, 0, 0, 0, 0, 0, 0,
        # fourth row
        0, 0, demographic_rates$F_2yr_surv, 0, 0, 0, 0, 0, 0, 0,
        #fifth row
        0, 0, 0, demographic_rates$F_3yr_surv, demographic_rates$F_adult_surv, 0, 0, 0, 0, 0,
        # sixth row
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        # seventh row
        0, 0, 0, 0, 0, demographic_rates$M_foal_surv, 0, 0, 0, 0,
        # eighth row
        0, 0, 0, 0, 0, 0, demographic_rates$M_yearling_surv, 0, 0, 0,
        # ninth row
        0, 0, 0, 0, 0, 0, 0, demographic_rates$M_2yr_surv, 0, 0,
        #tenth row
        0, 0, 0, 0, 0, 0, 0, 0, demographic_rates$M_3yr_surv, demographic_rates$M_adult_surv),
        nrow = length(stages), byrow = TRUE,
        dimnames = list(stages, stages))
    result
  }

surv_rates_treat <- data.frame(F_foal = 0.782,
                               F_yearling = 0.871,
                               F_2yr = 0.870,
                               F_3yr = 0.975,
                               F_adult =0.873,
                               M_foal = 0.792,
                               M_yearling = 0.930,
                               M_2yr = 0.899,
                               M_3yr = 0.975, 
                               M_adult = 0.925)

imm.rates_treat <- data.frame(F_foal = 0.086,
                              F_yearling = 0.105,
                              F_2yr = 0.147,
                              F_3yr = 0.120,
                              F_adult = 0.074,
                              M_foal = 0.123,
                              M_yearling = 0.136,
                              M_2yr = 0.069, 
                              M_3yr = 0.103, 
                              M_adult = 0.066)

emm.rates_treat <- data.frame(F_foal = 0.063,
                              F_yearling = 0.091,
                              F_2yr = 0.062,
                              F_3yr = 0.055,
                              F_adult = 0.081,
                              M_foal = 0.063,
                              M_yearling = 0.101,
                              M_2yr = 0.079, 
                              M_3yr = 0.101, 
                              M_adult = 0.072)

surv_rates_mprime <- data.frame(F_foal = (surv_rates_treat$F_foal + surv_rates_treat$M_foal)/2,
                                F_yearling = (surv_rates_treat$F_yearling + surv_rates_treat$M_yearling)/2,
                                F_2yr = (surv_rates_treat$F_2yr + surv_rates_treat$M_2yr)/2,
                                F_3yr = (surv_rates_treat$F_3yr + surv_rates_treat$M_3yr)/2,
                                F_adult = (surv_rates_treat$F_adult + surv_rates_treat$M_adult)/2,
                                M_foal = (surv_rates_treat$M_foal*2)/2,
                                M_yearling = (surv_rates_treat$M_yearling*2)/2,
                                M_2yr = (surv_rates_treat$M_2yr*2)/2, 
                                M_3yr  = (surv_rates_treat$M_3yr*2)/2, 
                                M_adult = (surv_rates_treat$M_adult*2)/2)

imm.rates_mprime <- data.frame(F_foal = (imm.rates_treat$F_foal + imm.rates_treat$M_foal)/2,
                               F_yearling = (imm.rates_treat$F_yearling + imm.rates_treat$M_yearling)/2,
                               F_2yr = (imm.rates_treat$F_2yr + imm.rates_treat$M_2yr)/2,
                               F_3yr = (imm.rates_treat$F_3yr + imm.rates_treat$M_3yr)/2,
                               F_adult = (imm.rates_treat$F_adult + imm.rates_treat$M_adult)/2,
                               M_foal = (imm.rates_treat$M_foal*2)/2,
                               M_yearling = (imm.rates_treat$M_yearling*2)/2,
                               M_2yr = (imm.rates_treat$M_2yr*2)/2, 
                               M_3yr = (imm.rates_treat$M_3yr*2)/2, 
                               M_adult = (imm.rates_treat$M_adult*2)/2)

emm.rates_mprime <- data.frame(F_foal = (emm.rates_treat$F_foal + emm.rates_treat$M_foal)/2,
                               F_yearling = (emm.rates_treat$F_yearling + emm.rates_treat$M_yearling)/2,
                               F_2yr = (emm.rates_treat$F_2yr + emm.rates_treat$M_2yr)/2,
                               F_3yr = (emm.rates_treat$F_3yr + emm.rates_treat$M_3yr)/2,
                               F_adult = (emm.rates_treat$F_adult + emm.rates_treat$M_adult)/2,
                               M_foal = (emm.rates_treat$M_foal*2)/2,
                               M_yearling = (emm.rates_treat$M_yearling*2)/2,
                               M_2yr = (emm.rates_treat$M_2yr*2)/2, 
                               M_3yr = (emm.rates_treat$M_3yr*2)/2, 
                               M_adult = (emm.rates_treat$M_adult*2)/2)

horse_treat <- list(F_foal_surv = surv_rates_treat$F_foal*(imm.rates_treat$F_foal + (1-emm.rates_treat$F_foal)),
                    F_yearling_surv = surv_rates_treat$F_yearling*(imm.rates_treat$F_yearling + (1-emm.rates_treat$F_yearling)),
                    F_2yr_surv = surv_rates_treat$F_2yr*(imm.rates_treat$F_2yr + (1-emm.rates_treat$F_2yr)),
                    F_3yr_surv = surv_rates_treat$F_3yr*(imm.rates_treat$F_3yr + (1-emm.rates_treat$F_3yr)),
                    F_adult_surv = surv_rates_treat$F_adult*(imm.rates_treat$F_adult + (1-emm.rates_treat$F_adult)),
                    M_foal_surv = surv_rates_treat$M_foal*(imm.rates_treat$M_foal + (1-emm.rates_treat$M_foal)),
                    M_yearling_surv = surv_rates_treat$M_yearling*(imm.rates_treat$M_yearling + (1-emm.rates_treat$M_yearling)),
                    M_2yr_surv = surv_rates_treat$M_2yr*(imm.rates_treat$M_2yr + (1-emm.rates_treat$M_2yr)), 
                    M_3yr_surv = surv_rates_treat$M_3yr*(imm.rates_treat$M_3yr + (1-emm.rates_treat$M_3yr)), 
                    M_adult_surv = surv_rates_treat$M_adult*(imm.rates_treat$M_adult + (1-emm.rates_treat$M_adult)),
                    h = 1.709,
                    k = 1,
                    FSR = 0.5)

horse_mprime_male <- list(F_foal_surv = ((surv_rates_treat$F_foal*(imm.rates_treat$F_foal + (1-emm.rates_treat$F_foal))) + (surv_rates_treat$M_foal*(imm.rates_treat$M_foal + (1-emm.rates_treat$M_foal))))/2,
                          F_yearling_surv = ((surv_rates_treat$F_yearling*(imm.rates_treat$F_yearling + (1-emm.rates_treat$F_yearling))) + (surv_rates_treat$M_yearling*(imm.rates_treat$M_yearling + (1-emm.rates_treat$M_yearling))))/2,
                          F_2yr_surv = ((surv_rates_treat$F_2yr*(imm.rates_treat$F_2yr + (1-emm.rates_treat$F_2yr))) + (surv_rates_treat$M_2yr*(imm.rates_treat$M_2yr + (1-emm.rates_treat$M_2yr))))/2,
                          F_3yr_surv = ((surv_rates_treat$F_3yr*(imm.rates_treat$F_3yr + (1-emm.rates_treat$F_3yr))) + (surv_rates_treat$M_3yr*(imm.rates_treat$M_3yr + (1-emm.rates_treat$M_3yr))))/2,
                          F_adult_surv = ((surv_rates_treat$F_adult*(imm.rates_treat$F_adult + (1-emm.rates_treat$F_adult))) + (surv_rates_treat$M_adult*(imm.rates_treat$M_adult + (1-emm.rates_treat$M_adult))))/2,
                          M_foal_surv = ((surv_rates_treat$M_foal*(imm.rates_treat$M_foal + (1-emm.rates_treat$M_foal))) + (surv_rates_treat$M_foal*(imm.rates_treat$M_foal + (1-emm.rates_treat$M_foal))))/2,
                          M_yearling_surv = ((surv_rates_treat$M_yearling*(imm.rates_treat$M_yearling + (1-emm.rates_treat$M_yearling))) + (surv_rates_treat$M_yearling*(imm.rates_treat$M_yearling + (1-emm.rates_treat$M_yearling))))/2,
                          M_2yr_surv = ((surv_rates_treat$M_2yr*(imm.rates_treat$M_2yr + (1-emm.rates_treat$M_2yr))) + (surv_rates_treat$M_2yr*(imm.rates_treat$M_2yr + (1-emm.rates_treat$M_2yr))))/2,
                          M_3yr_surv = ((surv_rates_treat$M_3yr*(imm.rates_treat$M_3yr + (1-emm.rates_treat$M_3yr))) + (surv_rates_treat$M_3yr*(imm.rates_treat$M_3yr + (1-emm.rates_treat$M_3yr))))/2,
                          M_adult_surv = ((surv_rates_treat$M_adult*(imm.rates_treat$M_adult + (1-emm.rates_treat$M_adult))) + (surv_rates_treat$M_adult*(imm.rates_treat$M_adult + (1-emm.rates_treat$M_adult))))/2,
                          h = (1.709*2)/2,
                          k = 1,
                          FSR = 0.5)

horse_treatment_matrix <- horse_matrix(horse_treat)

horse_mprime_matrix <- horse_matrix(horse_mprime_male)

# Estimating the ASR ------------------------------------------------------

matrix_ASR <-
  function(M, n = rep(10, nrow(M)), h = 1, k = 1, 
           iterations = 1000, FSR = 0.5, plot = FALSE){
    # Number of stages in matrix
    x <- length(n)
    # Number of time steps to simulate
    t <- iterations
    # an empty t by x matrix to store the stage distributions
    stage <- matrix(numeric(x * t), nrow = x)
    # an empty t vector to store the population sizes
    pop <- numeric(t)
    # for loop that goes through each of t time steps
    for (i in 1:t) {
      # stage distribution at time t
      stage[,i] <- n
      # population size at time t
      pop[i] <- sum(n)
      # number of male adults at time t
      M2 <- stage[10, i]
      # number of female adults at time t
      F2 <- stage[5, i]
      # Female freq-dep fecundity of female foals 
      M[1,x/2]        <- min(0.5*0.59, k*M2/(M2 + ((0.59*F2)/h)))*(1-FSR)
      # Female freq-dep fecundity of male foals
      M[(x/10)*6,x/2]  <- min(0.5*0.59, k*M2/(M2 + ((0.59*F2)/h)))*FSR
      # Male freq-dep fecundity of female foals
      M[1,x]          <- min(0.5*(0.59*F2)/M2, k*(0.59*F2)/(M2 + ((0.59*F2)/h)))*(1-FSR)
      # Male freq-dep fecundity of male foals
      M[(x/10)*6,x]    <- min(0.5*(0.59*F2)/M2, k*(0.59*F2)/(M2 + ((0.59*F2)/h)))*FSR
      # define the new n (i.e., new stage distribution at time t)
      n <- M %*% n
      # define rownames of stage matrix
      rownames(stage) <- rownames(M)
      # define colnames of stage matrix
      colnames(stage) <- 0:(t - 1)
      # calculate the proportional stable stage distribution - starts completely even as we start with 10 inds in each
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, t]
    }
    # calc ASR as the proportion of the adult stable stage class that is male
    ASR <- stable.stage[x]/(stable.stage[x/2] + stable.stage[x])
    
    if(plot)
    {
      # plot distrubution to assure that it is not chaotic
      matplot(rownames(t(stage)), t(stage), type='l', lwd=2, las=1)
    }
    # make a list of results
    pop.proj <- list(ASR = ASR,
                     stable.stage = stable.stage,
                     stage.vectors = stage,
                     SSD_M2 = stable.stage[10],
                     SSD_F2 = stable.stage[5])
    # print the list as output to the function
    pop.proj
  }

#treatment first
horse_treatment_ASR_analysis <- 
  matrix_ASR(M = horse_treatment_matrix, h = horse_treat$h, FSR = horse_treat$FSR, 
             iterations = 1000, plot = T)
horse_treatment_ASR <- horse_treatment_ASR_analysis$ASR
horse_treatment_ASR

#now the mprime matrix
horse_mprime_ASR_analysis <- 
  matrix_ASR(M = horse_mprime_matrix, h = horse_mprime_male$h, FSR = horse_mprime_male$FSR, 
             iterations = 1000, plot = T)
horse_mprime_ASR <- horse_mprime_ASR_analysis$ASR
horse_mprime_ASR

# Sensitivity analysis ----------------------------------------------------

matrix_structure <- expression(
  # top row of matrix
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  # second row of matrix
  F_foal, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  # third row of matrix
  0, F_1st_yr, 0, 0, 0, 0, 0, 0, 0, 0,
  # fourth row of matrix
  0, 0, F_2yr, 0, 0, 0, 0, 0, 0, 0,
  # fifth row
  0, 0, 0, F_3yr, F_Adt, 0, 0, 0, 0, 0,
  # sixth row
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  # seventh row
  0, 0, 0, 0, 0, M_foal, 0, 0, 0, 0,
  # eighth row
  0, 0, 0, 0, 0, 0, M_1st_yr, 0, 0, 0,
  #ninth row
  0, 0, 0, 0, 0, 0, 0, M_2yr, 0, 0,
  # tenth row
  0, 0, 0, 0, 0, 0, 0, 0, M_3yr, M_Adt)

sensitivity_analysis <-
  function(vital_rates, matrix_str, h, k, FSR, niter = 1000, ASR, imm_rates, emm_rates){
    
    # make a list of all parameters
    vr <- list(F_foal = vital_rates$F_foal,
               F_1st_yr = vital_rates$F_yearling,
               F_2yr = vital_rates$F_2yr,
               F_3yr = vital_rates$F_3yr,
               F_Adt = vital_rates$F_adult,
               M_foal = vital_rates$M_foal,
               M_1st_yr = vital_rates$M_yearling,
               M_2yr = vital_rates$M_2yr,
               M_3yr = vital_rates$M_3yr,
               M_Adt = vital_rates$M_adult)
    
    emm <- list(F_foal = emm_rates$F_foal,
                F_1st_yr = emm_rates$F_yearling,
                F_2yr = emm_rates$F_2yr,
                F_3yr = emm_rates$F_3yr,
                F_Adt = emm_rates$F_adult,
                M_foal = emm_rates$M_foal,
                M_1st_yr = emm_rates$M_yearling,
                M_2yr = emm_rates$M_yearling,
                M_3yr = emm_rates$M_yearling,
                M_Adt = emm_rates$M_adult)
    
    imm <- list(F_foal = imm_rates$F_foal,
                F_1st_yr = imm_rates$F_yearling,
                F_2yr = imm_rates$F_2yr,
                F_3yr = imm_rates$F_3yr,
                F_Adt = imm_rates$F_adult,
                M_foal = imm_rates$M_foal,
                M_1st_yr = imm_rates$M_yearling,
                M_2yr = imm_rates$M_yearling,
                M_3yr = imm_rates$M_yearling,
                M_Adt = imm_rates$M_adult)
    # number of stages in the matrix
    no_stages <- sqrt(length(matrix_str))
    
    # Define horse life stages
    stages <- c("F_foal", "F_1st_yr", "F_2yr", "F_3yr", "F_Adt", "M_foal", "M_1st_yr",  "M_2yr", "M_3yr", "M_Adt")
    
    # an empty t by x matrix
    stage <- matrix(numeric(no_stages * niter), nrow = no_stages)
    
    # an empty t vector to store the population sizes
    pop <- numeric(niter)
    
    # dataframe to store the perturbation results
    ASR_pert_results <-
      data.frame(parameter = c("F_foal_surv", "F_1st_yr_surv", "F_2yr_surv","F_3yr_surv", "F_Adt_surv",
                               "M_foal_surv", "M_1st_yr_surv", "M_2yr_surv", "M_3yr_surv", "M_Adt_surv",
                               "F_foal_imm", "F_1st_yr_imm", "F_2yr_imm","F_3yr_imm", "F_Adt_imm",
                               "M_foal_imm", "M_1st_yr_imm", "M_2yr_imm","M_3yr_imm", "M_Adt_imm", "F_foal_emm", "F_1st_yr_emm", "F_2yr_emm", "F_3yr_emm", "F_Adt_emm",
                               "M_foal_emm", "M_1st_yr_emm", "M_2yr_emm", "M_3yr_emm", "M_Adt_emm","h", "FSR"),
                 sensitivities = numeric(32),
                 elasticities = numeric(32))
    
    # specifiy how many survival rates there are
    n <- length(vr)
    
    # create vectors of perturbations to test on parameters of the matrix model
    vr_nums <- seq(0, 1, 0.01) # changes in survival and FSR (i.e., between 0 and 1)
    h_nums <- seq(0, 2, 0.02) # changes in h index (i.e., between 0 and 2)
    emm_nums <- seq(0, 1, 0.01) #changes in emigration
    imm_nums <- seq(0, 1, 0.01) #changes in immigration

    # create empty dataframes to store the perturbation results for ASR
    vr_pert_ASR <- matrix(numeric(n * length(vr_nums)),
                          ncol = n, dimnames = list(vr_nums, names(vr)))
    h_pert_ASR <- matrix(numeric(length(h_nums)),
                         ncol = 1, dimnames = list(h_nums, "h"))
    FSR_pert_ASR <- matrix(numeric(length(vr_nums)),
                           ncol = 1, dimnames = list(vr_nums, "FSR"))
    emm_pert_ASR <- matrix(numeric(n * length(emm_nums)),
                           ncol = n, dimnames = list(emm_nums, names(emm)))
    imm_pert_ASR <- matrix(numeric(n * length(imm_nums)),
                           ncol = n, dimnames = list(imm_nums, names(imm)))
    # perturbation of vital rates survival rates
    for (g in 1:n) # pick a column (i.e., a variable)
    {
      vr2 <- vr # reset the vital rates to the original
      for (i in 1:length(vr_nums)) # pick a perturbation level
      {
        vr2[[g]] <- vr_nums[i] # specify the vital rate with the new perturbation level
        A <- matrix(sapply(matrix_str, eval, vr2, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        I <- matrix(sapply(matrix_str, eval, imm, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        E <- matrix(sapply(matrix_str, eval, emm, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        
        Mov.M <- I + (1-E)
        # reset the starting stage distribution for simulation (all with 10 individuals)
        A <- A*Mov.M
        m <- rep(10, no_stages) 
        for (j in 1:niter) { # project the matrix through t iteration
          # stage distribution at time t
          stage[,j] <- m
          # population size at time t
          pop[j] <- sum(m)
          
          # number of male adults at time t
          M2 <- stage[10, j]
          # number of female adults at time t
          F2 <- stage[5, j]
          # Female freq-dep fecundity of female foals
          A[1,no_stages/2]        <- min((0.5*0.59), k*M2/(M2 + ((F2*0.59)/h)))*(1-FSR)
          # Female freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.59), k*M2/(M2 + ((F2*0.59)/h)))*FSR
          # Male freq-dep fecundity of female foals
          A[1,no_stages]          <- min(0.5*(0.59*F2)/M2, k*(0.59*F2)/(M2 + ((0.59*F2)/h)))*(1-FSR)
          # Male freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.59*F2)/M2, k*(0.59*F2)/(M2 + ((0.59*F2)/h)))*FSR
          # define the new n (i.e., new stage distribution at time t)
          m <- A %*% m
        }
        # define rownames of stage matrix
        rownames(stage) <- rownames(A)
        # define colnames of stage matrix
        colnames(stage) <- 0:(niter - 1)
        # calculate the proportional stable stage distribution
        stage <- apply(stage, 2, function(x) x/sum(x))
        # define stable stage as the last stage
        stable.stage <- stage[, niter]
        # calc ASR as the proportion of the adult stable stage class that is male
        vr_pert_ASR[i, g] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                        stable.stage[no_stages])
      }
      # get the spline function of ASR
      spl_ASR <- smooth.spline(vr_pert_ASR[,g] ~ rownames(vr_pert_ASR))
      # estimate the slope of the tangent of the spline at the vital rate
      ASR_pert_results[g, 2] <- predict(spl_ASR, x=vr[[g]], deriv=1)$y
      # re-scale sensitivity into elasticity
      ASR_pert_results[g, 3] <- vr[[g]]/ASR * ASR_pert_results[g, 2]
    }
    # perturbation of the h index parameter
    for (i in (1:length(h_nums))) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                  nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                  dimnames = list(stages, stages)) # build the matrix with the new value
      # reset the starting stage distribution for simulation (all with 10 individuals)
      
      A <- A*Mov.M
      
      m <- rep(10, no_stages) 
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        # population size at time t
        pop[j] <- sum(m)
        
        # number of male adults at time t
        M2 <- stage[10, j]
        # number of female adults at time t
        F2 <- stage[5, j]
        # Female freq-dep fecundity of female foals
        A[1,no_stages/2]        <- min((0.5*0.59), k*M2/(M2 + ((F2*0.59)/h_nums[i])))*(1-FSR)
        # Female freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.59), k*M2/(M2 + ((F2*0.59)/h_nums[i])))*FSR
        # Male freq-dep fecundity of female foals
        A[1,no_stages]          <- min(0.5*(0.59*F2)/M2, k*(0.59*F2)/(M2 + ((0.59*F2)/h_nums[i])))*(1-FSR)
        # Male freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.59*F2)/M2, k*(0.59*F2)/(M2 + ((0.59*F2)/h_nums[i])))*FSR
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% m
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      # calc ASR as the proportion of the adult stable stage class that is male
      h_pert_ASR[i,] <- 
        stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
    }
    # get the spline function of ASR
    spl_ASR <- smooth.spline(h_pert_ASR[, 1] ~ rownames(h_pert_ASR))
    # estimate the slope of the tangent of the spline at the vital rate
    ASR_pert_results[(n*3)+1, 2] <- predict(spl_ASR, x=h, deriv=1)$y
    # re-scale sensitivity into elasticity
    ASR_pert_results[(n*3)+1, 3] <- h/ASR * ASR_pert_results[(n*3)+1, 2]
    
    # perturbation of FSR
    for (i in 1:length(vr_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                  nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                  dimnames = list(stages, stages)) # build the matrix with the new value
      
      A <- A*Mov.M
      
      # reset the starting stage distribution for simulation (all with 10 individuals)
      m <- rep(10, no_stages) 
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        # population size at time t
        pop[j] <- sum(m)
        
        # number of male adults at time t
        M2 <- stage[10, j]
        # number of female adults at time t
        F2 <- stage[5, j]
        # Female freq-dep fecundity of female foals
        A[1,no_stages/2]        <- min((0.5*0.59), k*M2/(M2 + ((F2*0.59)/h)))*(1-vr_nums[i])
        # Female freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.59), k*M2/(M2 + ((F2*0.59)/h)))*vr_nums[i]
        # Male freq-dep fecundity of female foals
        A[1,no_stages]          <- min(0.5*(0.59*F2)/M2, k*(0.59*F2)/(M2 + ((0.59*F2)/h)))*(1-vr_nums[i])
        # Male freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.59*F2)/M2, k*(0.59*F2)/(M2 + ((0.59*F2)/h)))*vr_nums[i]
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% m
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      # calc ASR as the proportion of the adult stable stage class that is male
      FSR_pert_ASR[i,] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                     stable.stage[no_stages])
    }
    # get the spline function of ASR
    spl_ASR <- smooth.spline(FSR_pert_ASR[,1] ~ rownames(FSR_pert_ASR))
    # estimate the slope of the tangent of the spline at the vital rate    
    ASR_pert_results[(n*3)+2, 2] <- predict(spl_ASR, x=FSR, deriv=1)$y
    # re-scale sensitivity into elasticity
    ASR_pert_results[(n*3)+2, 3] <- FSR/ASR * ASR_pert_results[(n*3)+2, 2]
    
    #Now perturbing immigration
    for (g in 1:n) # pick a column (i.e., a variable)
    {
      imm2 <- imm # reset the vital rates to the original
      for (i in 1:length(imm_nums)) # pick a perturbation level
      {
        imm2[[g]] <- imm_nums[i] # specify the vital rate with the new perturbation level
        A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        I <- matrix(sapply(matrix_str, eval, imm2, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        E <- matrix(sapply(matrix_str, eval, emm, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        # reset the starting stage distribution for simulation (all with 10 individuals)
        Mov.M <- I + (1-E)
        
        A <- A*Mov.M
        
        m <- rep(10, no_stages) 
        for (j in 1:niter) { # project the matrix through t iteration
          # stage distribution at time t
          stage[,j] <- m
          # population size at time t
          pop[j] <- sum(m)
          
          # number of male adults at time t
          M2 <- stage[10, j]
          # number of female adults at time t
          F2 <- stage[5, j]
          # Female freq-dep fecundity of female foals
          A[1,no_stages/2]        <- min((0.5*0.59), k*M2/(M2 + ((F2*0.59)/h)))*(1-FSR)
          # Female freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.59), k*M2/(M2 + ((F2*0.59)/h)))*FSR
          # Male freq-dep fecundity of female foals
          A[1,no_stages]          <- min(0.5*(0.59*F2)/M2, k*(0.59*F2)/(M2 + ((0.59*F2)/h)))*(1-FSR)
          # Male freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.59*F2)/M2, k*(0.59*F2)/(M2 + ((0.59*F2)/h)))*FSR
          # define the new n (i.e., new stage distribution at time t)
          m <- A %*% m
        }
        # define rownames of stage matrix
        rownames(stage) <- rownames(A)
        # define colnames of stage matrix
        colnames(stage) <- 0:(niter - 1)
        # calculate the proportional stable stage distribution
        stage <- apply(stage, 2, function(x) x/sum(x))
        # define stable stage as the last stage
        stable.stage <- stage[, niter]
        # calc ASR as the proportion of the adult stable stage class that is male
        imm_pert_ASR[i, g] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                         stable.stage[no_stages])
      }
      # get the spline function of ASR
      spl_ASR <- smooth.spline(imm_pert_ASR[,g] ~ rownames(imm_pert_ASR))
      # estimate the slope of the tangent of the spline at the vital rate
      ASR_pert_results[n+g, 2] <- predict(spl_ASR, x=imm[[g]], deriv=1)$y
      # re-scale sensitivity into elasticity
      ASR_pert_results[n+g, 3] <- imm[[g]]/ASR * ASR_pert_results[n+g, 2]
    }
    
    #Now perturbing emmigration
    for (g in 1:n) # pick a column (i.e., a variable)
    {
      emm2 <- emm # reset the vital rates to the original
      for (i in 1:length(emm_nums)) # pick a perturbation level
      {
        emm2[[g]] <- emm_nums[i] # specify the vital rate with the new perturbation level
        A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        I <- matrix(sapply(matrix_str, eval, imm, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        E <- matrix(sapply(matrix_str, eval, emm2, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        # reset the starting stage distribution for simulation (all with 10 individuals)
        Mov.M <- I + (1-E)
        
        A <- A*Mov.M
        
        m <- rep(10, no_stages) 
        for (j in 1:niter) { # project the matrix through t iteration
          # stage distribution at time t
          stage[,j] <- m
          # population size at time t
          pop[j] <- sum(m)
          
          # number of male adults at time t
          M2 <- stage[10, j]
          # number of female adults at time t
          F2 <- stage[5, j]
          # Female freq-dep fecundity of female foals
          A[1,no_stages/2]        <- min((0.5*0.59), k*M2/(M2 + ((F2*0.59)/h)))*(1-FSR)
          # Female freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.59), k*M2/(M2 + ((F2*0.59)/h)))*FSR
          # Male freq-dep fecundity of female foals
          A[1,no_stages]          <- min(0.5*(0.59*F2)/M2, k*(0.59*F2)/(M2 + ((0.59*F2)/h)))*(1-FSR)
          # Male freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.59*F2)/M2, k*(0.59*F2)/(M2 + ((0.59*F2)/h)))*FSR
          # define the new n (i.e., new stage distribution at time t)
          m <- A %*% m
        }
        # define rownames of stage matrix
        rownames(stage) <- rownames(A)
        # define colnames of stage matrix
        colnames(stage) <- 0:(niter - 1)
        # calculate the proportional stable stage distribution
        stage <- apply(stage, 2, function(x) x/sum(x))
        # define stable stage as the last stage
        stable.stage <- stage[, niter]
        # calc ASR as the proportion of the adult stable stage class that is male
        emm_pert_ASR[i, g] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                         stable.stage[no_stages])
      }
      # get the spline function of ASR
      spl_ASR <- smooth.spline(emm_pert_ASR[,g] ~ rownames(emm_pert_ASR))
      # estimate the slope of the tangent of the spline at the vital rate
      ASR_pert_results[(n*2)+g, 2] <- predict(spl_ASR, x=emm[[g]], deriv=1)$y
      # re-scale sensitivity into elasticity
      ASR_pert_results[(n*2)+g, 3] <- emm[[g]]/ASR * ASR_pert_results[(n*2)+g, 2]
    }
    
    result <- list(ASR_pert_results = ASR_pert_results)
  }

horse_treat_sensitivity_analysis <- sensitivity_analysis(vital_rates = surv_rates_treat, matrix_str = matrix_structure, h = horse_treat$h, k = horse_treat$k, FSR = horse_treat$FSR, niter = 1000, ASR = horse_treatment_ASR, imm_rates = imm.rates_treat, emm_rates = emm.rates_treat)

horse_mprime_sensitivity_analysis_male <- sensitivity_analysis(vital_rates = surv_rates_mprime, matrix_str = matrix_structure, h = horse_mprime_male$h, k = horse_mprime_male$k, FSR = horse_mprime_male$FSR, niter = 1000, ASR = horse_mprime_ASR, imm_rates = imm.rates_mprime, emm_rates = emm.rates_mprime)

# Life table response experiment ------------------------------------------

LTRE_analysis <-
  function(Mprime_sens, surv_rates, imm_rates, emm_rates, sex, treat_h, treat_FSR){
    
    # make empty dataframes to stroe LTRE results for ASR
    LTRE_ASR <-
      data.frame(parameter = c("Foal survival", "Yearling survival", "2-year-old survival", "3-year-old survival",
                               "Adult survival", "Foal immigration", "Yearling immigration", "2-year-old immigration", "3-year-old immigration", "Adult immigration", "Foal emmigration", "Yearling emmigration", "2-year-old emmigration", "3-year-old emmigration", "Adult emmigration","Foal sex ratio",
                               "Mating system"),
                 contribution = numeric(17))
    
    # run a for loop to extract the parameter contributions
    if (sex == "male")
      # male rates scenario
      for(i in 1:nrow(LTRE_ASR))
      {
        LTRE_ASR[i, 2] <-
          ifelse(i < 6, (surv_rates[[i]] - surv_rates[[i + 5]]) * Mprime_sens$ASR_pert_results$sensitivities[i], #deviations in survival
                 ifelse(i >5 & i <11, (imm_rates[[i - 5]] - imm_rates[[i]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 5], #deviations in immigration
                        ifelse(i > 10 & i < 16, (emm_rates[[i - 10]] - emm_rates[[i - 5]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 10], #deviations in emigration
                               ifelse(i == 16, (treat_FSR - 0.5) * Mprime_sens$ASR_pert_results$sensitivities[i + 16], #deviations in foal sex ratio
                                      (treat_h - treat_h) * Mprime_sens$ASR_pert_results$sensitivities[i + 14])))) #and then deviations in the mating system
      }
    
    else
      # female rates scenario
      for(i in 1:nrow(LTRE_ASR))
      {
        LTRE_ASR[i, 2] <-
          ifelse(i < 6, (surv_rates[[i + 5]] - surv_rates[[i]]) *
                   Mprime_sens$ASR_pert_results$sensitivities[i + 5],
                 ifelse(i > 5 & i < 11, (imm_rates[[i]] - imm_rates[[i - 5]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 19],
                        ifelse(i > 10 & i < 16, (emm_rates[[i - 5]] - emm_rates[[i - 10]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 15],
                               ifelse(i == 16, (treat_FSR - 0.5) * Mprime_sens$ASR_pert_results$sensitivities[i + 16], 
                                      (treat_h - treat_h) * Mprime_sens$ASR_pert_results$sensitivities[i + 14]))))
      }
    
    LTRE_ASR$parameter <- factor(LTRE_ASR$parameter, levels = c("Adult survival",
                                                                "3-year-old survival",
                                                                "2-year-old survival",
                                                                "Yearling survival",
                                                                "Foal survival",
                                                                "Adult immigration",
                                                                "3-year-old immigration",
                                                                "2-year-old immigration",
                                                                "Yearling immigration",
                                                                "Foal immigration",
                                                                "Adult emmigration",
                                                                "3-year-old emmigration",
                                                                "2-year-old emmigration",
                                                                "Yearling emmigration",
                                                                "Foal emmigration",
                                                                "Foal sex ratio",
                                                                "Mating system"))
    
    LTRE_ASR$model <- "ASR"
    LTRE_df <- LTRE_ASR
    LTRE_df 
  }


LTRE_central <- LTRE_analysis(Mprime_sens = horse_mprime_sensitivity_analysis_male, surv_rates = surv_rates_treat, imm_rates = imm.rates_treat, emm_rates = emm.rates_treat, treat_FSR = horse_treat$FSR, treat_h = horse_treat$h, sex = "male")

LTRE_central

LTRE_contributions_check <- 
  function(surv_rates, imm_rates, emm_rates, treat_FSR, treat_h, Mprime_sensitivities, M_matrix_ASR, scenario) {
    if (scenario == "male")
      contribution_sum <- 
        sum(
          (surv_rates[[1]] - surv_rates[[6]]) * 
            Mprime_sensitivities[1],
          (surv_rates[[2]] - surv_rates[[7]]) * 
            Mprime_sensitivities[2],
          (surv_rates[[3]] - surv_rates[[8]]) * 
            Mprime_sensitivities[3],
          (surv_rates[[4]] - surv_rates[[9]]) * 
            Mprime_sensitivities[4],
          (surv_rates[[5]] - surv_rates[[10]]) * 
            Mprime_sensitivities[5],
          (surv_rates[[6]] - surv_rates[[6]]) * 
            Mprime_sensitivities[6],
          (surv_rates[[7]] - surv_rates[[7]]) * 
            Mprime_sensitivities[7],
          (surv_rates[[8]] - surv_rates[[8]]) * 
            Mprime_sensitivities[8],
          (surv_rates[[9]] - surv_rates[[9]]) * 
            Mprime_sensitivities[9],
          (surv_rates[[10]] - surv_rates[[10]]) * 
            Mprime_sensitivities[10],
          (imm_rates[[1]] - imm_rates[[6]]) * 
            Mprime_sensitivities[11],
          (imm_rates[[2]] - imm_rates[[7]]) * 
            Mprime_sensitivities[12],
          (imm_rates[[3]] - imm_rates[[8]]) * 
            Mprime_sensitivities[13],
          (imm_rates[[4]] - imm_rates[[9]]) * 
            Mprime_sensitivities[14],
          (imm_rates[[5]] - imm_rates[[10]]) * 
            Mprime_sensitivities[15],
          (imm_rates[[6]] - imm_rates[[6]]) * 
            Mprime_sensitivities[16],
          (imm_rates[[7]] - imm_rates[[7]]) * 
            Mprime_sensitivities[17],
          (imm_rates[[8]] - imm_rates[[8]]) * 
            Mprime_sensitivities[18],
          (imm_rates[[9]] - imm_rates[[10]]) * 
            Mprime_sensitivities[19],
          (imm_rates[[10]] - imm_rates[[10]]) * 
            Mprime_sensitivities[20],
          (emm_rates[[1]] - emm_rates[[6]]) * 
            Mprime_sensitivities[21],
          (emm_rates[[2]] - emm_rates[[7]]) * 
            Mprime_sensitivities[22],
          (emm_rates[[3]] - emm_rates[[8]]) * 
            Mprime_sensitivities[23],
          (emm_rates[[4]] - emm_rates[[9]]) * 
            Mprime_sensitivities[24],
          (emm_rates[[5]] - emm_rates[[10]]) * 
            Mprime_sensitivities[25],
          (emm_rates[[6]] - emm_rates[[6]]) * 
            Mprime_sensitivities[26],
          (emm_rates[[7]] - emm_rates[[7]]) * 
            Mprime_sensitivities[27],
          (emm_rates[[8]] - emm_rates[[8]]) * 
            Mprime_sensitivities[28],
          (emm_rates[[9]] - emm_rates[[9]]) * 
            Mprime_sensitivities[29],
          (emm_rates[[10]] - emm_rates[[10]]) * 
            Mprime_sensitivities[30],
          (treat_h - treat_h) *
            Mprime_sensitivities[31],
          (treat_FSR - 0.5) * Mprime_sensitivities[32]
        )
    
    else
      contribution_sum <- 
        sum(
          (surv_rates[[6]] - surv_rates[[1]]) * 
            Mprime_sensitivities[6],
          (surv_rates[[7]] - surv_rates[[2]]) * 
            Mprime_sensitivities[7],
          (surv_rates[[8]] - surv_rates[[3]]) * 
            Mprime_sensitivities[8],
          (surv_rates[[9]] - surv_rates[[4]]) * 
            Mprime_sensitivities[9],
          (surv_rates[[10]] - surv_rates[[5]]) * 
            Mprime_sensitivities[10],
          (surv_rates[[1]] - surv_rates[[1]]) * 
            Mprime_sensitivities[1],
          (surv_rates[[2]] - surv_rates[[2]]) * 
            Mprime_sensitivities[2],
          (surv_rates[[3]] - surv_rates[[3]]) * 
            Mprime_sensitivities[3],
          (surv_rates[[4]] - surv_rates[[4]]) * 
            Mprime_sensitivities[4],
          (surv_rates[[5]] - surv_rates[[5]]) * 
            Mprime_sensitivities[5],
          
          (imm_rates[[6]] - imm_rates[[1]]) * 
            Mprime_sensitivities[16],
          (imm_rates[[7]] - imm_rates[[2]]) * 
            Mprime_sensitivities[17],
          (imm_rates[[8]] - imm_rates[[3]]) * 
            Mprime_sensitivities[18],
          (imm_rates[[9]] - imm_rates[[4]]) * 
            Mprime_sensitivities[19],
          (imm_rates[[10]] - imm_rates[[1]]) * 
            Mprime_sensitivities[20],
          (imm_rates[[1]] - imm_rates[[1]]) * 
            Mprime_sensitivities[11],
          (imm_rates[[2]] - imm_rates[[2]]) * 
            Mprime_sensitivities[12],
          (imm_rates[[3]] - imm_rates[[3]]) * 
            Mprime_sensitivities[13],
          (imm_rates[[4]] - imm_rates[[4]]) * 
            Mprime_sensitivities[14],
          (imm_rates[[5]] - imm_rates[[5]]) * 
            Mprime_sensitivities[15],
          
          (emm_rates[[6]] - emm_rates[[1]]) * 
            Mprime_sensitivities[26],
          (emm_rates[[7]] - emm_rates[[2]]) * 
            Mprime_sensitivities[27],
          (emm_rates[[8]] - emm_rates[[3]]) * 
            Mprime_sensitivities[28],
          (emm_rates[[9]] - emm_rates[[4]]) * 
            Mprime_sensitivities[29],
          (emm_rates[[10]] - emm_rates[[5]]) * 
            Mprime_sensitivities[30],
          (emm_rates[[1]] - emm_rates[[1]]) * 
            Mprime_sensitivities[21],
          (emm_rates[[2]] - emm_rates[[2]]) * 
            Mprime_sensitivities[22],
          (emm_rates[[3]] - emm_rates[[3]]) * 
            Mprime_sensitivities[23],
          (emm_rates[[4]] - emm_rates[[4]]) * 
            Mprime_sensitivities[24],
          (emm_rates[[5]] - emm_rates[[5]]) * 
            Mprime_sensitivities[25],
          (treat_h - treat_h) *
            Mprime_sensitivities[31],
          (treat_FSR - 0.5) * Mprime_sensitivities[32]
        )
    
    ASR_bias <- abs(M_matrix_ASR - 0.5)
    absolute_difference <- abs(ASR_bias) - abs(contribution_sum)
    
    return(list(contribution_sum = as.vector(contribution_sum), 
                ASR_bias = as.vector(ASR_bias), 
                absolute_difference = as.vector(absolute_difference)))
  }

LTRE_contributions_check(surv_rates = surv_rates_treat, imm_rates = imm.rates_treat, emm_rates = emm.rates_treat, treat_h = horse_treat$h, treat_FSR = horse_treat$FSR, M_matrix_ASR = horse_treatment_ASR, Mprime_sensitivities = horse_mprime_sensitivity_analysis_male$ASR_pert_results$sensitivities, scenario = "male")

# Eastern subdivision -----------------------------------------------------

# creating matrices -------------------------------------------------------

horse_matrix <- 
  function(demographic_rates){
    # Define horse life stages
    stages <- c("F_foal", "F_1st_yr","F_2nd_yr", "F_3rd_yr", "F_Adt", "M_foal", "M_1st_yr","M_2nd_yr", "M_3rd_yr", "M_Adt")
    
    # Build the 10x10 matrix
    result <- 
      matrix(c(
        # top row
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        # second row
        demographic_rates$F_foal_surv,0, 0, 0, 0, 0, 0, 0, 0, 0,
        # third row
        0, demographic_rates$F_yearling_surv, 0, 0, 0, 0, 0, 0, 0, 0,
        # fourth row
        0, 0, demographic_rates$F_2yr_surv, 0, 0, 0, 0, 0, 0, 0,
        #fifth row
        0, 0, 0, demographic_rates$F_3yr_surv, demographic_rates$F_adult_surv, 0, 0, 0, 0, 0,
        # sixth row
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        # seventh row
        0, 0, 0, 0, 0, demographic_rates$M_foal_surv, 0, 0, 0, 0,
        # eighth row
        0, 0, 0, 0, 0, 0, demographic_rates$M_yearling_surv, 0, 0, 0,
        # ninth row
        0, 0, 0, 0, 0, 0, 0, demographic_rates$M_2yr_surv, 0, 0,
        #tenth row
        0, 0, 0, 0, 0, 0, 0, 0, demographic_rates$M_3yr_surv, demographic_rates$M_adult_surv),
        nrow = length(stages), byrow = TRUE,
        dimnames = list(stages, stages))
    result
  }

surv_rates_treat <- data.frame(F_foal = 0.762,
                               F_yearling = 0.886,
                               F_2yr = 0.804,
                               F_3yr = 0.954,
                               F_adult =0.854,
                               M_foal = 0.738,
                               M_yearling = 0.934,
                               M_2yr = 0.855, 
                               M_3yr = 0.981, 
                               M_adult = 0.910)

imm.rates_treat <- data.frame(F_foal = 0.016,
                              F_yearling = 0.050,
                              F_2yr = 0.053,
                              F_3yr = 0.059,
                              F_adult = 0.058,
                              M_foal = 0.050,
                              M_yearling = 0.109,
                              M_2yr = 0.112, 
                              M_3yr = 0.068, 
                              M_adult = 0.059)

emm.rates_treat <- data.frame(F_foal = 0.031,
                              F_yearling = 0.060,
                              F_2yr = 0.072,
                              F_3yr = 0.073,
                              F_adult = 0.042,
                              M_foal = 0.023,
                              M_yearling = 0.059,
                              M_2yr = 0.049, 
                              M_3yr = 0.099, 
                              M_adult = 0.048)

surv_rates_mprime <- data.frame(F_foal = (surv_rates_treat$F_foal + surv_rates_treat$M_foal)/2,
                                F_yearling = (surv_rates_treat$F_yearling + surv_rates_treat$M_yearling)/2,
                                F_2yr = (surv_rates_treat$F_2yr + surv_rates_treat$M_2yr)/2,
                                F_3yr = (surv_rates_treat$F_3yr + surv_rates_treat$M_3yr)/2,
                                F_adult = (surv_rates_treat$F_adult + surv_rates_treat$M_adult)/2,
                                M_foal = (surv_rates_treat$M_foal*2)/2,
                                M_yearling = (surv_rates_treat$M_yearling*2)/2,
                                M_2yr = (surv_rates_treat$M_2yr*2)/2, 
                                M_3yr  = (surv_rates_treat$M_3yr*2)/2, 
                                M_adult = (surv_rates_treat$M_adult*2)/2)

imm.rates_mprime <- data.frame(F_foal = (imm.rates_treat$F_foal + imm.rates_treat$M_foal)/2,
                               F_yearling = (imm.rates_treat$F_yearling + imm.rates_treat$M_yearling)/2,
                               F_2yr = (imm.rates_treat$F_2yr + imm.rates_treat$M_2yr)/2,
                               F_3yr = (imm.rates_treat$F_3yr + imm.rates_treat$M_3yr)/2,
                               F_adult = (imm.rates_treat$F_adult + imm.rates_treat$M_adult)/2,
                               M_foal = (imm.rates_treat$M_foal*2)/2,
                               M_yearling = (imm.rates_treat$M_yearling*2)/2,
                               M_2yr = (imm.rates_treat$M_2yr*2)/2, 
                               M_3yr = (imm.rates_treat$M_3yr*2)/2, 
                               M_adult = (imm.rates_treat$M_adult*2)/2)

emm.rates_mprime <- data.frame(F_foal = (emm.rates_treat$F_foal + emm.rates_treat$M_foal)/2,
                               F_yearling = (emm.rates_treat$F_yearling + emm.rates_treat$M_yearling)/2,
                               F_2yr = (emm.rates_treat$F_2yr + emm.rates_treat$M_2yr)/2,
                               F_3yr = (emm.rates_treat$F_3yr + emm.rates_treat$M_3yr)/2,
                               F_adult = (emm.rates_treat$F_adult + emm.rates_treat$M_adult)/2,
                               M_foal = (emm.rates_treat$M_foal*2)/2,
                               M_yearling = (emm.rates_treat$M_yearling*2)/2,
                               M_2yr = (emm.rates_treat$M_2yr*2)/2, 
                               M_3yr = (emm.rates_treat$M_3yr*2)/2, 
                               M_adult = (emm.rates_treat$M_adult*2)/2)

horse_treat <- list(F_foal_surv = surv_rates_treat$F_foal*(imm.rates_treat$F_foal + (1-emm.rates_treat$F_foal)),
                    F_yearling_surv = surv_rates_treat$F_yearling*(imm.rates_treat$F_yearling + (1-emm.rates_treat$F_yearling)),
                    F_2yr_surv = surv_rates_treat$F_2yr*(imm.rates_treat$F_2yr + (1-emm.rates_treat$F_2yr)),
                    F_3yr_surv = surv_rates_treat$F_3yr*(imm.rates_treat$F_3yr + (1-emm.rates_treat$F_3yr)),
                    F_adult_surv = surv_rates_treat$F_adult*(imm.rates_treat$F_adult + (1-emm.rates_treat$F_adult)),
                    M_foal_surv = surv_rates_treat$M_foal*(imm.rates_treat$M_foal + (1-emm.rates_treat$M_foal)),
                    M_yearling_surv = surv_rates_treat$M_yearling*(imm.rates_treat$M_yearling + (1-emm.rates_treat$M_yearling)),
                    M_2yr_surv = surv_rates_treat$M_2yr*(imm.rates_treat$M_2yr + (1-emm.rates_treat$M_2yr)), 
                    M_3yr_surv = surv_rates_treat$M_3yr*(imm.rates_treat$M_3yr + (1-emm.rates_treat$M_3yr)), 
                    M_adult_surv = surv_rates_treat$M_adult*(imm.rates_treat$M_adult + (1-emm.rates_treat$M_adult)),
                    h = 1.590,
                    k = 1,
                    FSR = 0.5)

horse_mprime_male <- list(F_foal_surv = ((surv_rates_treat$F_foal*(imm.rates_treat$F_foal + (1-emm.rates_treat$F_foal))) + (surv_rates_treat$M_foal*(imm.rates_treat$M_foal + (1-emm.rates_treat$M_foal))))/2,
                          F_yearling_surv = ((surv_rates_treat$F_yearling*(imm.rates_treat$F_yearling + (1-emm.rates_treat$F_yearling))) + (surv_rates_treat$M_yearling*(imm.rates_treat$M_yearling + (1-emm.rates_treat$M_yearling))))/2,
                          F_2yr_surv = ((surv_rates_treat$F_2yr*(imm.rates_treat$F_2yr + (1-emm.rates_treat$F_2yr))) + (surv_rates_treat$M_2yr*(imm.rates_treat$M_2yr + (1-emm.rates_treat$M_2yr))))/2,
                          F_3yr_surv = ((surv_rates_treat$F_3yr*(imm.rates_treat$F_3yr + (1-emm.rates_treat$F_3yr))) + (surv_rates_treat$M_3yr*(imm.rates_treat$M_3yr + (1-emm.rates_treat$M_3yr))))/2,
                          F_adult_surv = ((surv_rates_treat$F_adult*(imm.rates_treat$F_adult + (1-emm.rates_treat$F_adult))) + (surv_rates_treat$M_adult*(imm.rates_treat$M_adult + (1-emm.rates_treat$M_adult))))/2,
                          M_foal_surv = ((surv_rates_treat$M_foal*(imm.rates_treat$M_foal + (1-emm.rates_treat$M_foal))) + (surv_rates_treat$M_foal*(imm.rates_treat$M_foal + (1-emm.rates_treat$M_foal))))/2,
                          M_yearling_surv = ((surv_rates_treat$M_yearling*(imm.rates_treat$M_yearling + (1-emm.rates_treat$M_yearling))) + (surv_rates_treat$M_yearling*(imm.rates_treat$M_yearling + (1-emm.rates_treat$M_yearling))))/2,
                          M_2yr_surv = ((surv_rates_treat$M_2yr*(imm.rates_treat$M_2yr + (1-emm.rates_treat$M_2yr))) + (surv_rates_treat$M_2yr*(imm.rates_treat$M_2yr + (1-emm.rates_treat$M_2yr))))/2,
                          M_3yr_surv = ((surv_rates_treat$M_3yr*(imm.rates_treat$M_3yr + (1-emm.rates_treat$M_3yr))) + (surv_rates_treat$M_3yr*(imm.rates_treat$M_3yr + (1-emm.rates_treat$M_3yr))))/2,
                          M_adult_surv = ((surv_rates_treat$M_adult*(imm.rates_treat$M_adult + (1-emm.rates_treat$M_adult))) + (surv_rates_treat$M_adult*(imm.rates_treat$M_adult + (1-emm.rates_treat$M_adult))))/2,
                          h = (1.590*2)/2,
                          k = 1,
                          FSR = 0.5)

horse_treatment_matrix <- horse_matrix(horse_treat)

horse_mprime_matrix <- horse_matrix(horse_mprime_male)

# Estimating the ASR ------------------------------------------------------

matrix_ASR <-
  function(M, n = rep(10, nrow(M)), h = 1, k = 1, 
           iterations = 1000, FSR = 0.5, plot = FALSE){
    # Number of stages in matrix
    x <- length(n)
    # Number of time steps to simulate
    t <- iterations
    # an empty t by x matrix to store the stage distributions
    stage <- matrix(numeric(x * t), nrow = x)
    # an empty t vector to store the population sizes
    pop <- numeric(t)
    # for loop that goes through each of t time steps
    for (i in 1:t) {
      # stage distribution at time t
      stage[,i] <- n
      # population size at time t
      pop[i] <- sum(n)
      # number of male adults at time t
      M2 <- stage[10, i]
      # number of female adults at time t
      F2 <- stage[5, i]
      # Female freq-dep fecundity of female foals
      M[1,x/2]        <- min(0.5*0.61, k*M2/(M2 + ((0.61*F2)/h)))*(1-FSR)
      # Female freq-dep fecundity of male foals
      M[(x/10)*6,x/2]  <- min(0.5*0.61, k*M2/(M2 + ((0.61*F2)/h)))*FSR
      # Male freq-dep fecundity of female foals
      M[1,x]          <- min(((0.5*(0.61*F2))/M2), k*(0.61*F2)/(M2 + ((0.61*F2)/h)))*(1-FSR)
      # Male freq-dep fecundity of male foals
      M[(x/10)*6,x]    <- min(((0.5*(0.61*F2))/M2), k*(0.61*F2)/(M2 + ((0.61*F2)/h)))*FSR
      # define the new n (i.e., new stage distribution at time t)
      n <- M %*% n
      # define rownames of stage matrix
      rownames(stage) <- rownames(M)
      # define colnames of stage matrix
      colnames(stage) <- 0:(t - 1)
      # calculate the proportional stable stage distribution - starts completely even as we start with 10 inds in each
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, t]
    }
    # calc ASR as the proportion of the adult stable stage class that is male
    ASR <- stable.stage[x]/(stable.stage[x/2] + stable.stage[x])
    
    if(plot)
    {
      # plot distrubution to assure that it is not chaotic
      matplot(rownames(t(stage)), t(stage), type='l', lwd=2, las=1)
    }
    # make a list of results
    pop.proj <- list(ASR = ASR,
                     stable.stage = stable.stage,
                     stage.vectors = stage,
                     SSD_M2 = stable.stage[10],
                     SSD_F2 = stable.stage[5])
    # print the list as output to the function
    pop.proj
  }

#treatment first
horse_treatment_ASR_analysis <- 
  matrix_ASR(M = horse_treatment_matrix, h = horse_treat$h, FSR = horse_treat$FSR, 
             iterations = 1000, plot = T)
horse_treatment_ASR <- horse_treatment_ASR_analysis$ASR
horse_treatment_ASR

#now the mprime matrix
horse_mprime_ASR_analysis <- 
  matrix_ASR(M = horse_mprime_matrix, h = horse_mprime_male$h, FSR = horse_mprime_male$FSR, 
             iterations = 1000, plot = T)
horse_mprime_ASR <- horse_mprime_ASR_analysis$ASR
horse_mprime_ASR

# Sensitivity analysis ----------------------------------------------------

matrix_structure <- expression(
  # top row of matrix
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  # second row of matrix
  F_foal, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  # third row of matrix
  0, F_1st_yr, 0, 0, 0, 0, 0, 0, 0, 0,
  # fourth row of matrix
  0, 0, F_2yr, 0, 0, 0, 0, 0, 0, 0,
  # fifth row
  0, 0, 0, F_3yr, F_Adt, 0, 0, 0, 0, 0,
  # sixth row
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  # seventh row
  0, 0, 0, 0, 0, M_foal, 0, 0, 0, 0,
  # eighth row
  0, 0, 0, 0, 0, 0, M_1st_yr, 0, 0, 0,
  #ninth row
  0, 0, 0, 0, 0, 0, 0, M_2yr, 0, 0,
  # tenth row
  0, 0, 0, 0, 0, 0, 0, 0, M_3yr, M_Adt)

sensitivity_analysis <-
  function(vital_rates, matrix_str, h, k, FSR, niter = 1000, ASR, imm_rates, emm_rates){
    
    # make a list of all parameters
    vr <- list(F_foal = vital_rates$F_foal,
               F_1st_yr = vital_rates$F_yearling,
               F_2yr = vital_rates$F_2yr,
               F_3yr = vital_rates$F_3yr,
               F_Adt = vital_rates$F_adult,
               M_foal = vital_rates$M_foal,
               M_1st_yr = vital_rates$M_yearling,
               M_2yr = vital_rates$M_2yr,
               M_3yr = vital_rates$M_3yr,
               M_Adt = vital_rates$M_adult)
    
    emm <- list(F_foal = emm_rates$F_foal,
                F_1st_yr = emm_rates$F_yearling,
                F_2yr = emm_rates$F_2yr,
                F_3yr = emm_rates$F_3yr,
                F_Adt = emm_rates$F_adult,
                M_foal = emm_rates$M_foal,
                M_1st_yr = emm_rates$M_yearling,
                M_2yr = emm_rates$M_yearling,
                M_3yr = emm_rates$M_yearling,
                M_Adt = emm_rates$M_adult)
    
    imm <- list(F_foal = imm_rates$F_foal,
                F_1st_yr = imm_rates$F_yearling,
                F_2yr = imm_rates$F_2yr,
                F_3yr = imm_rates$F_3yr,
                F_Adt = imm_rates$F_adult,
                M_foal = imm_rates$M_foal,
                M_1st_yr = imm_rates$M_yearling,
                M_2yr = imm_rates$M_yearling,
                M_3yr = imm_rates$M_yearling,
                M_Adt = imm_rates$M_adult)
    # number of stages in the matrix
    no_stages <- sqrt(length(matrix_str))
    
    # Define horse life stages
    stages <- c("F_foal", "F_1st_yr", "F_2yr", "F_3yr", "F_Adt", "M_foal", "M_1st_yr",  "M_2yr", "M_3yr", "M_Adt")
    
    # an empty t by x matrix
    stage <- matrix(numeric(no_stages * niter), nrow = no_stages)
    
    # an empty t vector to store the population sizes
    pop <- numeric(niter)
    
    # dataframe to store the perturbation results
    ASR_pert_results <-
      data.frame(parameter = c("F_foal_surv", "F_1st_yr_surv", "F_2yr_surv","F_3yr_surv", "F_Adt_surv",
                               "M_foal_surv", "M_1st_yr_surv", "M_2yr_surv", "M_3yr_surv", "M_Adt_surv",
                               "F_foal_imm", "F_1st_yr_imm", "F_2yr_imm","F_3yr_imm", "F_Adt_imm",
                               "M_foal_imm", "M_1st_yr_imm", "M_2yr_imm","M_3yr_imm", "M_Adt_imm", "F_foal_emm", "F_1st_yr_emm", "F_2yr_emm", "F_3yr_emm", "F_Adt_emm",
                               "M_foal_emm", "M_1st_yr_emm", "M_2yr_emm", "M_3yr_emm", "M_Adt_emm","h", "FSR"),
                 sensitivities = numeric(32),
                 elasticities = numeric(32))
    
    # specifiy how many survival rates there are
    n <- length(vr)
    
    # create vectors of perturbations to test on parameters of the matrix model
    vr_nums <- seq(0, 1, 0.01) # changes in survival and FSR (i.e., between 0 and 1)
    h_nums <- seq(0, 2, 0.02) # changes in h index (i.e., between 0 and 2)
    emm_nums <- seq(0, 1, 0.01) #changes in emigration
    imm_nums <- seq(0, 1, 0.01) #changes in immigration

    # create empty dataframes to store the perturbation results for ASR
    vr_pert_ASR <- matrix(numeric(n * length(vr_nums)),
                          ncol = n, dimnames = list(vr_nums, names(vr)))
    h_pert_ASR <- matrix(numeric(length(h_nums)),
                         ncol = 1, dimnames = list(h_nums, "h"))
    FSR_pert_ASR <- matrix(numeric(length(vr_nums)),
                           ncol = 1, dimnames = list(vr_nums, "FSR"))
    emm_pert_ASR <- matrix(numeric(n * length(emm_nums)),
                           ncol = n, dimnames = list(emm_nums, names(emm)))
    imm_pert_ASR <- matrix(numeric(n * length(imm_nums)),
                           ncol = n, dimnames = list(imm_nums, names(imm)))
    # perturbation of vital rates survival rates
    for (g in 1:n) # pick a column (i.e., a variable)
    {
      vr2 <- vr # reset the vital rates to the original
      for (i in 1:length(vr_nums)) # pick a perturbation level
      {
        vr2[[g]] <- vr_nums[i] # specify the vital rate with the new perturbation level
        A <- matrix(sapply(matrix_str, eval, vr2, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        I <- matrix(sapply(matrix_str, eval, imm, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        E <- matrix(sapply(matrix_str, eval, emm, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        Mov.M <- I + (1-E)
        # reset the starting stage distribution for simulation (all with 10 individuals)

        A <- A*Mov.M
        
        m <- rep(10, no_stages) 
        for (j in 1:niter) { # project the matrix through t iteration
          # stage distribution at time t
          stage[,j] <- m
          # population size at time t
          pop[j] <- sum(m)
          
          # number of male adults at time t
          M2 <- stage[10, j]
          # number of female adults at time t
          F2 <- stage[5, j]
          # Female freq-dep fecundity of female foals
          A[1,no_stages/2]        <- min((0.5*0.61), k*M2/(M2 + ((F2*0.61)/h)))*(1-FSR)
          # Female freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.61), k*M2/(M2 + ((F2*0.61)/h)))*FSR
          # Male freq-dep fecundity of female foals
          A[1,no_stages]          <- min(0.5*(0.61*F2)/M2, k*(0.61*F2)/(M2 + ((0.61*F2)/h)))*(1-FSR)
          # Male freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.61*F2)/M2, k*(0.61*F2)/(M2 + ((0.61*F2)/h)))*FSR
          # define the new n (i.e., new stage distribution at time t)
          m <- A %*% m
        }
        # define rownames of stage matrix
        rownames(stage) <- rownames(A)
        # define colnames of stage matrix
        colnames(stage) <- 0:(niter - 1)
        # calculate the proportional stable stage distribution
        stage <- apply(stage, 2, function(x) x/sum(x))
        # define stable stage as the last stage
        stable.stage <- stage[, niter]
        # calc ASR as the proportion of the adult stable stage class that is male
        vr_pert_ASR[i, g] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                        stable.stage[no_stages])
      }
      # get the spline function of ASR
      spl_ASR <- smooth.spline(vr_pert_ASR[,g] ~ rownames(vr_pert_ASR))
      # estimate the slope of the tangent of the spline at the vital rate
      ASR_pert_results[g, 2] <- predict(spl_ASR, x=vr[[g]], deriv=1)$y
      # re-scale sensitivity into elasticity
      ASR_pert_results[g, 3] <- vr[[g]]/ASR * ASR_pert_results[g, 2]
    }
    # perturbation of the h index parameter
    for (i in (1:length(h_nums))) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                  nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                  dimnames = list(stages, stages)) # build the matrix with the new value
      # reset the starting stage distribution for simulation (all with 10 individuals)
      
      A <- A*Mov.M
      
      m <- rep(10, no_stages) 
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        # population size at time t
        pop[j] <- sum(m)
        
        # number of male adults at time t
        M2 <- stage[10, j]
        # number of female adults at time t
        F2 <- stage[5, j]
        # Female freq-dep fecundity of female foals
        A[1,no_stages/2]        <- min((0.5*0.61), k*M2/(M2 + ((F2*0.61)/h_nums[i])))*(1-FSR)
        # Female freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.61), k*M2/(M2 + ((F2*0.61)/h_nums[i])))*FSR
        # Male freq-dep fecundity of female foals
        A[1,no_stages]          <- min(0.5*(0.61*F2)/M2, k*(0.61*F2)/(M2 + ((0.61*F2)/h_nums[i])))*(1-FSR)
        # Male freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.61*F2)/M2, k*(0.61*F2)/(M2 + ((0.61*F2)/h_nums[i])))*FSR
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% m
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      # calc ASR as the proportion of the adult stable stage class that is male
      h_pert_ASR[i,] <- 
        stable.stage[no_stages]/(stable.stage[no_stages/2] + stable.stage[no_stages])
    }
    # get the spline function of ASR
    spl_ASR <- smooth.spline(h_pert_ASR[, 1] ~ rownames(h_pert_ASR))
    # estimate the slope of the tangent of the spline at the vital rate
    ASR_pert_results[(n*3)+1, 2] <- predict(spl_ASR, x=h, deriv=1)$y
    # re-scale sensitivity into elasticity
    ASR_pert_results[(n*3)+1, 3] <- h/ASR * ASR_pert_results[(n*3)+1, 2]
    
    # perturbation of FSR
    for (i in 1:length(vr_nums)) # pick a perturbation level
    {
      A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                  nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                  dimnames = list(stages, stages)) # build the matrix with the new value
      
      A <- A*Mov.M
      
      # reset the starting stage distribution for simulation (all with 10 individuals)
      m <- rep(10, no_stages) 
      for (j in 1:niter) { # project the matrix through t iteration
        # stage distribution at time t
        stage[,j] <- m
        # population size at time t
        pop[j] <- sum(m)
        
        # number of male adults at time t
        M2 <- stage[10, j]
        # number of female adults at time t
        F2 <- stage[5, j]
        # Female freq-dep fecundity of female foals
        A[1,no_stages/2]        <- min((0.5*0.61), k*M2/(M2 + ((F2*0.61)/h)))*(1-vr_nums[i])
        # Female freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.61), k*M2/(M2 + ((F2*0.61)/h)))*vr_nums[i]
        # Male freq-dep fecundity of female foals
        A[1,no_stages]          <- min(0.5*(0.61*F2)/M2, k*(0.61*F2)/(M2 + ((0.61*F2)/h)))*(1-vr_nums[i])
        # Male freq-dep fecundity of male foals
        A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.61*F2)/M2, k*(0.61*F2)/(M2 + ((0.61*F2)/h)))*vr_nums[i]
        # define the new n (i.e., new stage distribution at time t)
        m <- A %*% m
      }
      # define rownames of stage matrix
      rownames(stage) <- rownames(A)
      # define colnames of stage matrix
      colnames(stage) <- 0:(niter - 1)
      # calculate the proportional stable stage distribution
      stage <- apply(stage, 2, function(x) x/sum(x))
      # define stable stage as the last stage
      stable.stage <- stage[, niter]
      # calc ASR as the proportion of the adult stable stage class that is male
      FSR_pert_ASR[i,] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                     stable.stage[no_stages])
    }
    # get the spline function of ASR
    spl_ASR <- smooth.spline(FSR_pert_ASR[,1] ~ rownames(FSR_pert_ASR))
    # estimate the slope of the tangent of the spline at the vital rate    
    ASR_pert_results[(n*3)+2, 2] <- predict(spl_ASR, x=FSR, deriv=1)$y
    # re-scale sensitivity into elasticity
    ASR_pert_results[(n*3)+2, 3] <- FSR/ASR * ASR_pert_results[(n*3)+2, 2]
    
    #Now perturbing immigration
    for (g in 1:n) # pick a column (i.e., a variable)
    {
      imm2 <- imm # reset the vital rates to the original
      for (i in 1:length(imm_nums)) # pick a perturbation level
      {
        imm2[[g]] <- imm_nums[i] # specify the vital rate with the new perturbation level
        A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        I <- matrix(sapply(matrix_str, eval, imm2, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        E <- matrix(sapply(matrix_str, eval, emm, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        # reset the starting stage distribution for simulation (all with 10 individuals)
        Mov.M <- I + (1-E)
        
        A <- A*Mov.M
        
        m <- rep(10, no_stages) 
        for (j in 1:niter) { # project the matrix through t iteration
          # stage distribution at time t
          stage[,j] <- m
          # population size at time t
          pop[j] <- sum(m)
          
          # number of male adults at time t
          M2 <- stage[10, j]
          # number of female adults at time t
          F2 <- stage[5, j]
          # Female freq-dep fecundity of female foals
          A[1,no_stages/2]        <- min((0.5*0.61), k*M2/(M2 + ((F2*0.61)/h)))*(1-FSR)
          # Female freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.61), k*M2/(M2 + ((F2*0.61)/h)))*FSR
          # Male freq-dep fecundity of female foals
          A[1,no_stages]          <- min(0.5*(0.61*F2)/M2, k*(0.61*F2)/(M2 + ((0.61*F2)/h)))*(1-FSR)
          # Male freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.61*F2)/M2, k*(0.61*F2)/(M2 + ((0.61*F2)/h)))*FSR
          # define the new n (i.e., new stage distribution at time t)
          m <- A %*% m
        }
        # define rownames of stage matrix
        rownames(stage) <- rownames(A)
        # define colnames of stage matrix
        colnames(stage) <- 0:(niter - 1)
        # calculate the proportional stable stage distribution
        stage <- apply(stage, 2, function(x) x/sum(x))
        # define stable stage as the last stage
        stable.stage <- stage[, niter]
        # calc ASR as the proportion of the adult stable stage class that is male
        imm_pert_ASR[i, g] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                         stable.stage[no_stages])
      }
      # get the spline function of ASR
      spl_ASR <- smooth.spline(imm_pert_ASR[,g] ~ rownames(imm_pert_ASR))
      # estimate the slope of the tangent of the spline at the vital rate
      ASR_pert_results[n+g, 2] <- predict(spl_ASR, x=imm[[g]], deriv=1)$y
      # re-scale sensitivity into elasticity
      ASR_pert_results[n+g, 3] <- imm[[g]]/ASR * ASR_pert_results[n+g, 2]
    }
    
    #Now perturbing emmigration
    for (g in 1:n) # pick a column (i.e., a variable)
    {
      emm2 <- emm # reset the vital rates to the original
      for (i in 1:length(emm_nums)) # pick a perturbation level
      {
        emm2[[g]] <- emm_nums[i] # specify the vital rate with the new perturbation level
        A <- matrix(sapply(matrix_str, eval, vr, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        I <- matrix(sapply(matrix_str, eval, imm, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        E <- matrix(sapply(matrix_str, eval, emm2, NULL), 
                    nrow = sqrt(length(matrix_str)), byrow=TRUE, 
                    dimnames = list(stages, stages)) # build the matrix with the new value
        
        # reset the starting stage distribution for simulation (all with 10 individuals)
        Mov.M <- I + (1-E)
        
        A <- A*Mov.M
        
        m <- rep(10, no_stages) 
        for (j in 1:niter) { # project the matrix through t iteration
          # stage distribution at time t
          stage[,j] <- m
          # population size at time t
          pop[j] <- sum(m)
          
          # number of male adults at time t
          M2 <- stage[10, j]
          # number of female adults at time t
          F2 <- stage[5, j]
          # Female freq-dep fecundity of female foals
          A[1,no_stages/2]        <- min((0.5*0.61), k*M2/(M2 + ((F2*0.61)/h)))*(1-FSR)
          # Female freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages/2]  <- min((0.5*0.61), k*M2/(M2 + ((F2*0.61)/h)))*FSR
          # Male freq-dep fecundity of female foals
          A[1,no_stages]          <- min(0.5*(0.61*F2)/M2, k*(0.61*F2)/(M2 + ((0.61*F2)/h)))*(1-FSR)
          # Male freq-dep fecundity of male foals
          A[(no_stages/10)*6,no_stages]    <- min(0.5*(0.61*F2)/M2, k*(0.61*F2)/(M2 + ((0.61*F2)/h)))*FSR
          # define the new n (i.e., new stage distribution at time t)
          m <- A %*% m
        }
        # define rownames of stage matrix
        rownames(stage) <- rownames(A)
        # define colnames of stage matrix
        colnames(stage) <- 0:(niter - 1)
        # calculate the proportional stable stage distribution
        stage <- apply(stage, 2, function(x) x/sum(x))
        # define stable stage as the last stage
        stable.stage <- stage[, niter]
        # calc ASR as the proportion of the adult stable stage class that is male
        emm_pert_ASR[i, g] <- stable.stage[no_stages]/(stable.stage[no_stages/2] + 
                                                         stable.stage[no_stages])
      }
      # get the spline function of ASR
      spl_ASR <- smooth.spline(emm_pert_ASR[,g] ~ rownames(emm_pert_ASR))
      # estimate the slope of the tangent of the spline at the vital rate
      ASR_pert_results[(n*2)+g, 2] <- predict(spl_ASR, x=emm[[g]], deriv=1)$y
      # re-scale sensitivity into elasticity
      ASR_pert_results[(n*2)+g, 3] <- emm[[g]]/ASR * ASR_pert_results[(n*2)+g, 2]
    }
    
    result <- list(ASR_pert_results = ASR_pert_results)
  }

horse_treat_sensitivity_analysis <- sensitivity_analysis(vital_rates = surv_rates_treat, matrix_str = matrix_structure, h = horse_treat$h, k = horse_treat$k, FSR = horse_treat$FSR, niter = 1000, ASR = horse_treatment_ASR, imm_rates = imm.rates_treat, emm_rates = emm.rates_treat)

horse_mprime_sensitivity_analysis_male <- sensitivity_analysis(vital_rates = surv_rates_mprime, matrix_str = matrix_structure, h = horse_mprime_male$h, k = horse_mprime_male$k, FSR = horse_mprime_male$FSR, niter = 1000, ASR = horse_mprime_ASR, imm_rates = imm.rates_mprime, emm_rates = emm.rates_mprime)

# Life table response experiment ------------------------------------------

LTRE_analysis <-
  function(Mprime_sens, surv_rates, imm_rates, emm_rates, sex, treat_h, treat_FSR){
    
    # make empty dataframes to stroe LTRE results for ASR
    LTRE_ASR <-
      data.frame(parameter = c("Foal survival", "Yearling survival", "2-year-old survival", "3-year-old survival",
                               "Adult survival", "Foal immigration", "Yearling immigration", "2-year-old immigration", "3-year-old immigration", "Adult immigration", "Foal emmigration", "Yearling emmigration", "2-year-old emmigration", "3-year-old emmigration", "Adult emmigration","Foal sex ratio",
                               "Mating system"),
                 contribution = numeric(17))
    
    # run a for loop to extract the parameter contributions
    if (sex == "male")
      # male rates scenario
      for(i in 1:nrow(LTRE_ASR))
      {
        LTRE_ASR[i, 2] <-
          ifelse(i < 6, (surv_rates[[i]] - surv_rates[[i + 5]]) * Mprime_sens$ASR_pert_results$sensitivities[i], #deviations in survival
                 ifelse(i >5 & i <11, (imm_rates[[i - 5]] - imm_rates[[i]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 5], #deviations in immigration
                        ifelse(i > 10 & i < 16, (emm_rates[[i - 10]] - emm_rates[[i - 5]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 10], #deviations in emigration
                               ifelse(i == 16, (treat_FSR - 0.5) * Mprime_sens$ASR_pert_results$sensitivities[i + 16], #deviations in foal sex ratio
                                      (treat_h - treat_h) * Mprime_sens$ASR_pert_results$sensitivities[i + 14])))) #and then deviations in the mating system
      }
    
    else
      # female rates scenario
      for(i in 1:nrow(LTRE_ASR))
      {
        LTRE_ASR[i, 2] <-
          ifelse(i < 6, (surv_rates[[i + 5]] - surv_rates[[i]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 5],
                 ifelse(i > 5 & i < 11, (imm_rates[[i]] - imm_rates[[i - 5]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 19],
                        ifelse(i > 10 & i < 16, (emm_rates[[i - 5]] - emm_rates[[i - 10]]) * Mprime_sens$ASR_pert_results$sensitivities[i + 15],
                               ifelse(i == 16, (treat_FSR - 0.5) * Mprime_sens$ASR_pert_results$sensitivities[i + 16], 
                                      (treat_h - treat_h) * Mprime_sens$ASR_pert_results$sensitivities[i + 14]))))
      }
    
    LTRE_ASR$parameter <- factor(LTRE_ASR$parameter, levels = c("Adult survival",
                                                                "3-year-old survival",
                                                                "2-year-old survival",
                                                                "Yearling survival",
                                                                "Foal survival",
                                                                "Adult immigration",
                                                                "3-year-old immigration",
                                                                "2-year-old immigration",
                                                                "Yearling immigration",
                                                                "Foal immigration",
                                                                "Adult emmigration",
                                                                "3-year-old emmigration",
                                                                "2-year-old emmigration",
                                                                "Yearling emmigration",
                                                                "Foal emmigration",
                                                                "Foal sex ratio",
                                                                "Mating system"))
    
    LTRE_ASR$model <- "ASR"
    LTRE_df <- LTRE_ASR
    LTRE_df 
  }

LTRE_east <- LTRE_analysis(Mprime_sens = horse_mprime_sensitivity_analysis_male, surv_rates = surv_rates_treat, imm_rates = imm.rates_treat, emm_rates = emm.rates_treat, treat_FSR = horse_treat$FSR, treat_h = horse_treat$h, sex = "male")

LTRE_east

LTRE_contributions_check <- 
  function(surv_rates, imm_rates, emm_rates, treat_FSR, treat_h, Mprime_sensitivities, M_matrix_ASR, scenario) {
    if (scenario == "male")
      contribution_sum <- 
        sum(
          (surv_rates[[1]] - surv_rates[[6]]) * 
            Mprime_sensitivities[1],
          (surv_rates[[2]] - surv_rates[[7]]) * 
            Mprime_sensitivities[2],
          (surv_rates[[3]] - surv_rates[[8]]) * 
            Mprime_sensitivities[3],
          (surv_rates[[4]] - surv_rates[[9]]) * 
            Mprime_sensitivities[4],
          (surv_rates[[5]] - surv_rates[[10]]) * 
            Mprime_sensitivities[5],
          (surv_rates[[6]] - surv_rates[[6]]) * 
            Mprime_sensitivities[6],
          (surv_rates[[7]] - surv_rates[[7]]) * 
            Mprime_sensitivities[7],
          (surv_rates[[8]] - surv_rates[[8]]) * 
            Mprime_sensitivities[8],
          (surv_rates[[9]] - surv_rates[[9]]) * 
            Mprime_sensitivities[9],
          (surv_rates[[10]] - surv_rates[[10]]) * 
            Mprime_sensitivities[10],
          (imm_rates[[1]] - imm_rates[[6]]) * 
            Mprime_sensitivities[11],
          (imm_rates[[2]] - imm_rates[[7]]) * 
            Mprime_sensitivities[12],
          (imm_rates[[3]] - imm_rates[[8]]) * 
            Mprime_sensitivities[13],
          (imm_rates[[4]] - imm_rates[[9]]) * 
            Mprime_sensitivities[14],
          (imm_rates[[5]] - imm_rates[[10]]) * 
            Mprime_sensitivities[15],
          (imm_rates[[6]] - imm_rates[[6]]) * 
            Mprime_sensitivities[16],
          (imm_rates[[7]] - imm_rates[[7]]) * 
            Mprime_sensitivities[17],
          (imm_rates[[8]] - imm_rates[[8]]) * 
            Mprime_sensitivities[18],
          (imm_rates[[9]] - imm_rates[[10]]) * 
            Mprime_sensitivities[19],
          (imm_rates[[10]] - imm_rates[[10]]) * 
            Mprime_sensitivities[20],
          (emm_rates[[1]] - emm_rates[[6]]) * 
            Mprime_sensitivities[21],
          (emm_rates[[2]] - emm_rates[[7]]) * 
            Mprime_sensitivities[22],
          (emm_rates[[3]] - emm_rates[[8]]) * 
            Mprime_sensitivities[23],
          (emm_rates[[4]] - emm_rates[[9]]) * 
            Mprime_sensitivities[24],
          (emm_rates[[5]] - emm_rates[[10]]) * 
            Mprime_sensitivities[25],
          (emm_rates[[6]] - emm_rates[[6]]) * 
            Mprime_sensitivities[26],
          (emm_rates[[7]] - emm_rates[[7]]) * 
            Mprime_sensitivities[27],
          (emm_rates[[8]] - emm_rates[[8]]) * 
            Mprime_sensitivities[28],
          (emm_rates[[9]] - emm_rates[[9]]) * 
            Mprime_sensitivities[29],
          (emm_rates[[10]] - emm_rates[[10]]) * 
            Mprime_sensitivities[30],
          (treat_h - treat_h) *
            Mprime_sensitivities[31],
          (treat_FSR - 0.5) * Mprime_sensitivities[32]
        )
    
    else
      contribution_sum <- 
        sum(
          (surv_rates[[6]] - surv_rates[[1]]) * 
            Mprime_sensitivities[6],
          (surv_rates[[7]] - surv_rates[[2]]) * 
            Mprime_sensitivities[7],
          (surv_rates[[8]] - surv_rates[[3]]) * 
            Mprime_sensitivities[8],
          (surv_rates[[9]] - surv_rates[[4]]) * 
            Mprime_sensitivities[9],
          (surv_rates[[10]] - surv_rates[[5]]) * 
            Mprime_sensitivities[10],
          (surv_rates[[1]] - surv_rates[[1]]) * 
            Mprime_sensitivities[1],
          (surv_rates[[2]] - surv_rates[[2]]) * 
            Mprime_sensitivities[2],
          (surv_rates[[3]] - surv_rates[[3]]) * 
            Mprime_sensitivities[3],
          (surv_rates[[4]] - surv_rates[[4]]) * 
            Mprime_sensitivities[4],
          (surv_rates[[5]] - surv_rates[[5]]) * 
            Mprime_sensitivities[5],
          
          (imm_rates[[6]] - imm_rates[[1]]) * 
            Mprime_sensitivities[16],
          (imm_rates[[7]] - imm_rates[[2]]) * 
            Mprime_sensitivities[17],
          (imm_rates[[8]] - imm_rates[[3]]) * 
            Mprime_sensitivities[18],
          (imm_rates[[9]] - imm_rates[[4]]) * 
            Mprime_sensitivities[19],
          (imm_rates[[10]] - imm_rates[[1]]) * 
            Mprime_sensitivities[20],
          (imm_rates[[1]] - imm_rates[[1]]) * 
            Mprime_sensitivities[11],
          (imm_rates[[2]] - imm_rates[[2]]) * 
            Mprime_sensitivities[12],
          (imm_rates[[3]] - imm_rates[[3]]) * 
            Mprime_sensitivities[13],
          (imm_rates[[4]] - imm_rates[[4]]) * 
            Mprime_sensitivities[14],
          (imm_rates[[5]] - imm_rates[[5]]) * 
            Mprime_sensitivities[15],
          
          (emm_rates[[6]] - emm_rates[[1]]) * 
            Mprime_sensitivities[26],
          (emm_rates[[7]] - emm_rates[[2]]) * 
            Mprime_sensitivities[27],
          (emm_rates[[8]] - emm_rates[[3]]) * 
            Mprime_sensitivities[28],
          (emm_rates[[9]] - emm_rates[[4]]) * 
            Mprime_sensitivities[29],
          (emm_rates[[10]] - emm_rates[[5]]) * 
            Mprime_sensitivities[30],
          (emm_rates[[1]] - emm_rates[[1]]) * 
            Mprime_sensitivities[21],
          (emm_rates[[2]] - emm_rates[[2]]) * 
            Mprime_sensitivities[22],
          (emm_rates[[3]] - emm_rates[[3]]) * 
            Mprime_sensitivities[23],
          (emm_rates[[4]] - emm_rates[[4]]) * 
            Mprime_sensitivities[24],
          (emm_rates[[5]] - emm_rates[[5]]) * 
            Mprime_sensitivities[25],
          (treat_h - treat_h) *
            Mprime_sensitivities[31],
          (treat_FSR - 0.5) * Mprime_sensitivities[32]
        )
    
    ASR_bias <- abs(M_matrix_ASR - 0.5)
    absolute_difference <- abs(ASR_bias) - abs(contribution_sum)
    
    return(list(contribution_sum = as.vector(contribution_sum), 
                ASR_bias = as.vector(ASR_bias), 
                absolute_difference = as.vector(absolute_difference)))
  }

LTRE_contributions_check(surv_rates = surv_rates_treat, imm_rates = imm.rates_treat, emm_rates = emm.rates_treat, treat_h = horse_treat$h, treat_FSR = horse_treat$FSR, M_matrix_ASR = horse_treatment_ASR, Mprime_sensitivities = horse_mprime_sensitivity_analysis_male$ASR_pert_results$sensitivities, scenario = "male")

