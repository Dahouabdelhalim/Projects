model {
# Detection model ----
for(st in 1:n_sites) {
  for(ti in 1:n_time_periods) {
    for(vt in 1:visits_by_site[st, ti]) {
      eh[vt, st, ti] ~ dcat(dpm[incid[st, ti], , vt, st, ti])
    }
  }
  # Latent state model ----
  # for first season
  incid[st, 1] ~ dcat(occ_st[ , st])
  for(ti in 2:n_time_periods){
    incid[st, ti] ~ dcat(trans[incid[st, ti - 1], , st, ti - 1] )
  }
}

# Build the vectors/matrices for the above models.
for (st in 1:n_sites) {
  # Initial occupancy states (occ_st) ----
  # Model assumes initial occupancy is independent of the incidence of the
  #  other species.
  # All occupancy states at a site sum to 1.
  occ_st[1, st] <- (1 - occ[1, st]) * (1 - occ[2, st]) #U
  occ_st[2, st] <-      occ[1, st]  * (1 - occ[2, st]) #A
  occ_st[3, st] <- (1 - occ[1, st]) *      occ[2, st]  #B
  occ_st[4, st] <-      occ[1, st]  *      occ[2, st]  #AB
  # Transition probability matrix (trans) ----
  # The indexes on the probabilities (colonization (col), presistence (prs),
  #   continued absence (abs), extinction (ext)) correspond to the following
  #   condintional states: (1) A|B, (2) A|b, (3) B|A, (4) B|a. Where lowercase
  #   denotes the absence of the species in the next time step (t + 1).
  # Rows sum to 1.
  for (ti in 1:(n_time_periods - 1)) {
    # Create new variables to make filling of matrix tidier
    # Abs == continued absence of a species at a site
    abs[1:4, st, ti] <- 1 - col[1:4, st, ti]
    # Ext == extinction of a species at a site
    ext[1:4, st, ti] <- 1 - prs[1:4, st, ti]
    # Row 1 == initial state unoccupied
    trans[1, 1, st, ti] <- abs[2, st, ti] * abs[4, st, ti]
    trans[1, 2, st, ti] <- col[2, st, ti] * abs[3, st, ti]
    trans[1, 3, st, ti] <- abs[1, st, ti] * col[4, st, ti]
    trans[1, 4, st, ti] <- col[1, st, ti] * col[3, st, ti]
    # Row 2 == initial state A only
    trans[2, 1, st, ti] <- ext[2, st, ti] * abs[4, st, ti]
    trans[2, 2, st, ti] <- prs[2, st, ti] * abs[3, st, ti]
    trans[2, 3, st, ti] <- ext[1, st, ti] * col[4, st, ti]
    trans[2, 4, st, ti] <- prs[1, st, ti] * col[3, st, ti]
    # Row 3 == initial state B only
    trans[3, 1, st, ti] <- abs[2, st, ti] * ext[4, st, ti]
    trans[3, 2, st, ti] <- col[2, st, ti] * ext[3, st, ti]
    trans[3, 3, st, ti] <- abs[1, st, ti] * prs[4, st, ti]
    trans[3, 4, st, ti] <- col[1, st, ti] * prs[3, st, ti]
    # Row 4 == initial state A and B
    trans[4, 1, st, ti] <- ext[2, st, ti] * ext[4, st, ti]
    trans[4, 2, st, ti] <- prs[2, st, ti] * ext[3, st, ti]
    trans[4, 3, st, ti] <- ext[1, st, ti] * prs[4, st, ti]
    trans[4, 4, st, ti] <- prs[1, st, ti] * prs[3, st, ti]
  }
  # Detection probability matrix ----
  for (ti in 1:n_time_periods) {
    for (vt in 1:visits_by_site[st, ti]) {
      # Model assumes detection of one species is independent of the incidence of
      #   the other species. Rows sum to 1.
      # Row 1 == occ_st[1]: A and B are absent
      dpm[1, 1, vt, st, ti] <- 1
      dpm[1, 2, vt, st, ti] <- 0
      dpm[1, 3, vt, st, ti] <- 0
      dpm[1, 4, vt, st, ti] <- 0
      # Row 2 == occ_st[2]: A is present. B is absent.
      dpm[2, 1, vt, st, ti] <- (1 - det[1, vt, st, ti])
      dpm[2, 2, vt, st, ti] <- det[1, vt, st, ti]
      dpm[2, 3, vt, st, ti] <- 0
      dpm[2, 4, vt, st, ti] <- 0
      # Row 3 == occ_st[3]: A is absent. B is present.
      dpm[3, 1, vt, st, ti] <- 1 - det[2, vt, st, ti]
      dpm[3, 2, vt, st, ti] <- 0
      dpm[3, 3, vt, st, ti] <- det[2, vt, st, ti]
      dpm[3, 4, vt, st, ti] <- 0
      # Row 4 == occ_st[4]: A and B are both present.
      dpm[4, 1, vt, st, ti] <- (1 - det[1, vt, st, ti]) * (1 - det[2, vt, st, ti])
      dpm[4, 2, vt, st, ti] <-      det[1, vt, st, ti]  * (1 - det[2, vt, st, ti])
      dpm[4, 3, vt, st, ti] <- (1 - det[1, vt, st, ti]) *      det[2, vt, st, ti]
      dpm[4, 4, vt, st, ti] <-      det[1, vt, st, ti]  *      det[2, vt, st, ti]
    }
  }
}

# Linear predictors ------------------------------------------------------------
# Logit-links for occupancy, colonization, persistence, and detection
for (st in 1:n_sites) {
  for (sp in 1:n_spp) {
    logit(occ[sp, st]) <- occ_inter[sp] + inprod(occ_coef[sp, ],
                                                 occ_covs[st, ])
  }
  for (ti in 1:(n_time_periods - 1)) {
    # Colonization
    # A|B
    logit(col[1, st, ti]) <- col_inter[1, ti] + inx_col[ti]
    # A|b
    logit(col[2, st, ti]) <- col_inter[1, ti]
    # B|A
    logit(col[3, st, ti]) <- col_inter[2, ti] + inx_col[ti]
    # B|a
    logit(col[4, st, ti]) <- col_inter[2, ti]
    # Persistence
    # A|B
    logit(prs[1, st, ti]) <- prs_inter[1, ti] + inx_prs[ti]
    # A|b
    logit(prs[2, st, ti]) <- prs_inter[1, ti]
    # B|A
    logit(prs[3, st, ti]) <- prs_inter[2, ti] + inx_prs[ti]
    # B|a
    logit(prs[4, st, ti]) <- prs_inter[2, ti]
  }
}
# Covariate index to tidy logit-link
det_ind <- n_det_covs - n_sl_det_covs
for (sp in 1:n_spp) {
  for (st in 1:n_sites) {
    for (ti in 1:n_time_periods) {
      for (vt in 1:visits_by_site[st, ti]) {
        logit(det[sp, vt, st, ti]) <- det_inter[sp] +
          inprod(det_coef[sp, 1:det_ind],
                 det_covs[vt, st, ti, ]) +
          inprod(det_coef[sp, (det_ind + 1):n_det_covs], sl_det_covs[st, ti, ])
      }
    }
  }
}

# Priors -----------------------------------------------------------------------
# Shape parameter
tau <- 1 / (2.25 ^ 2)
# Species interactions
for (ti in 1:(n_time_periods - 1)) {
  inx_col[ti] ~ dnorm(0, tau)
  inx_prs[ti] ~ dnorm(0, tau)
}
for (sp in 1:n_spp){
  occ_inter[sp] ~ dnorm(0, tau)
  for (cov in 1:n_occ_covs) {
    occ_coef[sp, cov] ~ dnorm(0, tau)
  }
  det_inter[sp] ~ dnorm(0, tau)
  for (cov in 1:n_det_covs) {
    det_coef[sp, cov] ~ dnorm(0, tau)
  }
  for (ti in 1:(n_time_periods - 1)) {
    col_inter[sp, ti] ~ dnorm(0, tau)
    prs_inter[sp, ti] ~ dnorm(0, tau)
  }
}
}


