model{
det_ind <- n_det_covs - n_sl_det_covs
# Community-level Intercept Priors
# Community-level prior shape parameters
tau <- 1 / (2.25 ^ 2)
occ_inter_sig ~ dt(0, tau, 1) T(0, )
det_inter_sig ~ dt(0, tau, 1) T(0, )
surv_inter_sig ~ dt(0, tau, 1) T(0, )
colon_inter_sig ~ dt(0, tau, 1) T(0, )

# Prior for occupancy because only the initial occupancy probability is estimated
# by the model and any additional occupancy values are derived, only one prior is
# specified.
occ_inter_loc ~ dnorm(0, tau)
occ_inter_scale <- pow(occ_inter_sig, -2)

# Prior for detection probability
det_inter_loc ~ dnorm(0, tau)
det_inter_scale <- pow(det_inter_sig, -2)

# Priors for persistence and colonization
for (time in 1:(n_time_periods - 1)) {
  surv_inter_loc[time] ~ dnorm(0, tau)
  surv_inter_scale[time] <- pow(surv_inter_sig, -2)
  colon_inter_loc[time] ~ dnorm(0, tau)
  colon_inter_scale[time] <- pow(colon_inter_sig, -2)
}

# Community-level priors: occupancy coefficients
for (cov in 1:n_occ_covs) {
  occ_coef_sig[cov] ~ dt(0, tau, 1) T(0, )
  occ_coef_loc[cov] ~ dnorm(0, tau)
  occ_coef_scale[cov] <- pow(occ_coef_sig[cov], -2)
}
# Community-level priors: detection coefficients
for (cov in 1:n_det_covs) {
  det_coef_sig[cov] ~ dt(0, tau, 1) T(0, )
  det_coef_loc[cov] ~ dnorm(0, tau)
  det_coef_scale[cov] <- pow(det_coef_sig[cov], -2)
}

# Community-level priors: persistence and survival coefficients
for (cov in 1:n_colon_covs) {
  colon_coef_sig[cov] ~ dt(0, tau, 1) T(0, )
}
for (cov in 1:n_surv_covs) {
  surv_coef_sig[cov] ~ dt(0, tau, 1) T(0, )
}
for (time in 1:(n_time_periods - 1)) {
  for (cov in 1:n_colon_covs) {
    colon_coef_loc[time, cov] ~ dnorm(0, tau)
    colon_coef_scale[time, cov] <- pow(colon_coef_sig[cov], -2)
  }
  for (cov in 1:n_surv_covs) {
    surv_coef_loc[time, cov] ~ dnorm(0, tau)
    surv_coef_scale[time, cov] <- pow(surv_coef_sig[cov], -2)
  }
}
# Beginning of model
for (sp in 1:n_spp) {
  occ_inter[sp] ~ dnorm(occ_inter_loc, occ_inter_scale)
  for (cov in 1:n_occ_covs) {
    occ_coef[sp, cov] ~ dnorm(occ_coef_loc[cov], occ_coef_scale[cov])
  }
  for (time in 1:(n_time_periods - 1)) {
    colon_inter[sp, time] ~ dnorm(colon_inter_loc[time],
                                   colon_inter_scale[time])
    for (cov in 1:n_colon_covs) {
      colon_coef[sp, time, cov] ~ dnorm(colon_coef_loc[time, cov],
                                         colon_coef_scale[time, cov])
    }
    surv_inter[sp, time] ~ dnorm(surv_inter_loc[time], surv_inter_scale[time])
    for (cov in 1:n_surv_covs) {
      surv_coef[sp, time, cov] ~ dnorm(surv_coef_loc[time, cov], 
                                        surv_coef_scale[time, cov])
    }
  }
  # Initial occupancy state at time = 1
  for (site in 1:n_sites) {
    logit(occ[sp, site, 1]) <- occ_inter[sp] +
      inprod(occ_coef[sp, ], occ_covs[site, ])
    mu_incid[sp, site, 1] <- occ[sp, site, 1]
    incid[sp, site, 1] ~ dbern(mu_incid[sp, site, 1])
    # Model of changes in occ state for time=2, ... n
    for (time in 1:(n_time_periods - 1)) {
      logit(colon[sp, site, time]) <- colon_inter[sp, time] +
        inprod(colon_coef[sp, time, ], colon_covs[site, time, ])
      logit(surv[sp, site, time]) <- surv_inter[sp, time] + 
        inprod(surv_coef[sp, time, ], surv_covs[site, time, ])
      # Occupancy after the first time step is a derived value, and though done
      # here, it is not necessary to estimate within the model. This simplifies
      # how the code appears when detection is calculated.
      occ[sp, site, time + 1] <- surv[sp, site, time] *
        occ[sp, site, time] + colon[sp, site, time] *
        (1 - occ[sp, site, time])
      mu_incid[sp, site, time + 1] <- surv[sp, site, time] *
        incid[sp, site, time] + colon[sp, site, time] *
        (1 - incid[sp, site, time])
      incid[sp, site, time + 1] ~ dbern(mu_incid[sp, site, time + 1])
    }
  }
  # Detection
  det_inter[sp] ~ dnorm(det_inter_loc, det_inter_scale)
  for (cov in 1:n_det_covs) {
    det_coef[sp, cov] ~ dnorm(det_coef_loc[cov], det_coef_scale[cov])
  }
  for (site in 1:n_sites) {
    for (time in 1:n_time_periods) {
      for (visit in 1:visits_by_site[site, time]) {
        logit(det[sp, visit, site, time]) <- det_inter[sp] + 
          inprod(det_coef[sp, 1:det_ind], 
                 det_covs[visit, site, time, ]) +
          inprod(det_coef[sp, (det_ind + 1):n_det_covs],
                 sl_det_covs[site, time, ])
        mu_eh[sp, visit, site, time] <- det[sp, visit, site, time] *
          incid[sp, site, time]
        eh[sp, visit, site, time] ~ dbern(mu_eh[sp, visit, site, time])
        # Simulate the detection matrix predicted by the model
        eh_sim[sp, visit, site, time] ~ dbern(mu_eh[sp, visit, site, time])
      }
    }
  }
}
}
