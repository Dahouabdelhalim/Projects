
model{
  # Community-level (hyper) priors --------------------------------------------
  tau <- pow(2.25, -2)
  # # Intercepts:
  # Intial occupancy
  occ_int_location ~ dnorm(0, tau)
  occ_int_sigma ~ dt(0, tau, 1) T(0, )
  occ_int_scale <- pow(occ_int_sigma, -2)
  # Detection
  det_int_location ~ dnorm(0, tau)
  det_int_sigma ~ dt(0, tau, 1) T(0, )
  det_int_scale <- pow(det_int_sigma, -2)
  # Survival and colonization
  surv_int_sigma ~ dt(0, tau, 1) T(0, )
  colon_int_sigma ~ dt(0, tau, 1) T(0, )
  for (time in 1:(n_time_periods - 1)) {
    surv_int_location[time] ~ dnorm(0, tau)
    surv_int_scale[time] <- pow(surv_int_sigma, -2)
    colon_int_location[time] ~ dnorm(0, tau)
    colon_int_scale[time] <- pow(colon_int_sigma, -2)
  }
  # # Coefficients:
  # Initial occupancy
  for (cov in 1:n_occ_covs) {
    occ_coef_sigma[cov] ~ dt(0, tau, 1) T(0, )
    occ_coef_location[cov] ~ dnorm(0, tau)
    occ_coef_scale[cov] <- pow(occ_coef_sigma[cov], -2)
  }
  # Detection
  for (cov in 1:n_det_covs) {
    det_coef_sigma[cov] ~ dt(0, tau, 1) T(0, )
    det_coef_location[cov] ~ dnorm(0, tau)
    det_coef_scale[cov] <- pow(det_coef_sigma[cov], -2)
  }
  # Survival and colonization
  for (cov in 1:n_surv_covs) {
    surv_coef_sigma[cov] ~ dt(0, tau, 1) T(0, )
  }
  for (cov in 1:n_colon_covs) {
    colon_coef_sigma[cov] ~ dt(0, tau, 1) T(0, )
  }
  for (time in 1:(n_time_periods - 1)) {
    for (cov in 1:n_surv_covs) {
      surv_coef_location[time, cov] ~ dnorm(0, tau)
      surv_coef_scale[time, cov] <- pow(surv_coef_sigma[cov], -2)
    }
    for (cov in 1:n_colon_covs) {
      colon_coef_location[time, cov] ~ dnorm(0, tau)
      colon_coef_scale[time, cov] <- pow(colon_coef_sigma[cov], -2)
    }
  }
  # Likelihood loop -----------------------------------------------------------
  for (species in 1:n_species) {
    # # Species-level intercept and coefficient priors:
    # Occupancy
    occ_int[species] ~ dnorm(occ_int_location, occ_int_scale)
    for (cov in 1:n_occ_covs) {
      occ_coef[species, cov] ~ dnorm(occ_coef_location[cov],
                                     occ_coef_scale[cov])
    }
    # Survival and colonization
    for (time in 1:(n_time_periods - 1)) {
      surv_int[species, time] ~ dnorm(surv_int_location[time],
                                      surv_int_scale[time])
      for (cov in 1:n_surv_covs) {
        surv_coef[species, time, cov] ~ dnorm(surv_coef_location[time, cov],
                                              surv_coef_scale[time, cov])
      }
      colon_int[species, time] ~ dnorm(colon_int_location[time],
                                       colon_int_scale[time])
      for (cov in 1:n_colon_covs) {
        colon_coef[species, time, cov] ~ dnorm(colon_coef_location[time, cov],
                                               colon_coef_scale[time, cov])
      }
    }
    for (site in 1:n_sites) {
      # # Regressions:
      # Initial occupancy state (time = 1)
      logit(occ[species, site, 1]) <- occ_int[species] +
        inprod(occ_coef[species, ], occ_covs[site, ])
      mean_incidence[species, site, 1] <- occ[species, site, 1]
      # Model of changes in occupancy state for time=2, ..., n
      for (time in 1:(n_time_periods - 1)) {
        logit(surv[species, site, time]) <- surv_int[species, time] +
          inprod(surv_coef[species, time, ], surv_covs[site, time, ])
        logit(colon[species, site, time]) <- colon_int[species, time] +
          inprod(colon_coef[species, time, ], colon_covs[site, time, ])
        mean_incidence[species, site, time + 1] <- surv[species, site, time] *
          incidence[species, site, time] + colon[species, site, time] *
          (1 - incidence[species, site, time])
      }
    }
    # # Detection loops:
    # Species-level intercept and coefficient priors
    det_int[species] ~ dnorm(det_int_location, det_int_scale)
    for (cov in 1:n_det_covs) {
      det_coef[species, cov] ~ dnorm(det_coef_location[cov],
                                     det_coef_scale[cov])
    }
    for (site in 1:n_sites) {
      for (time in 1:n_time_periods) {
        for (visit in 1:visits_by_site[site, time]) {
          # Regression	
          logit(det[species, visit, site, time]) <- det_int[species] +
            inprod(det_coef[species, 1:(n_det_covs - n_site_det_covs)], det_covs[visit, site, time, ]) +
            inprod(det_coef[species, (n_det_covs - n_site_det_covs+1):n_det_covs], site_det_covs[site, time, ])
          mean_det_matrix[species, visit, site, time] <-
            det[species, visit, site, time] * incidence[species, site, time]
          # Generate a likelihood for the MCMC sampler by linking the
          # model to the field data
          det_matrix[species, visit, site, time] ~
            dbern(mean_det_matrix[species, visit, site, time])
          # Simulate the detection matrix that is predicted by the model
          # (used to calculate PPC)
          det_matrix_simulated[species, visit, site, time] ~
            dbern(mean_det_matrix[species, visit, site, time])
        }
        # Partially-latent incidence matrix
        incidence[species, site, time] ~
          dbern(mean_incidence[species, site, time])
      }
    }
  }
}

