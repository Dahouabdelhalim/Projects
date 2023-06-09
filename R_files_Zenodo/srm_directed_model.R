model_code_0 <- "
data{
    int N;

    int N_i__j_kt_ID;
    int N_ij_kt_ID;
    int N_i2j_kt_ID;
    // int N_k_ID; // commented out
    
    int y[N];
    int denominator[N];
    
    int ikt_ID[N];
    int jkt_ID[N];
    
    int ijkt_ID[N];
    int i2jkt_ID[N];
    
    // int k_ID[N]; // commented out

}

parameters{

    real a;

    ////////// Variances and covariances
    
    // Node-level
    vector <lower=0> [2] sigma_i__j_kt_ID;
    cholesky_factor_corr[2] L_Rho_i__j_kt_ID;
    
    // Symmetric dyads
    real <lower=0> sigma_ij_kt_ID;
    
    // Directed dyads
    real <lower=0> sigma_i2j_kt_ID;
    
    // Group
    // real <lower=0> sigma_k_ID; // commented out

    /////////// Effects
    matrix[2 , N_i__j_kt_ID] z_i__j_kt_ID; // Node-level varying intercepts
    vector[N_ij_kt_ID] v_ij_kt_ID; // Symmetric dyad effects
    vector[N_i2j_kt_ID] v_i2j_kt_ID; // Directed dyad effects
    // vector[N_k_ID] v_k_ID; // Group effects // commented out

}

transformed parameters{

    // Transformed effects
    matrix [N_i__j_kt_ID,2] v_i__j_kt_ID;
    
    // Calculate effects from Cholesky matrices
    v_i__j_kt_ID = (diag_pre_multiply ( sigma_i__j_kt_ID, L_Rho_i__j_kt_ID ) * z_i__j_kt_ID)';

}

model{

    vector[N] p;

    //priors for fixed effects
    a ~ normal( -4, 5 );

    // priors on varying effects
    to_vector(z_i__j_kt_ID) ~ normal( 0,1 );
    to_vector(v_ij_kt_ID) ~ normal( 0,1 );
    to_vector(v_i2j_kt_ID) ~ normal( 0,1 );
    // to_vector(v_k_ID) ~ normal( 0,1 ); // commented out
    
    // Sigmas and covariance
    sigma_i__j_kt_ID ~ exponential( 1 );
    L_Rho_i__j_kt_ID ~ lkj_corr_cholesky( 4 );
    
    sigma_ij_kt_ID ~ exponential( 1 );
    sigma_i2j_kt_ID ~ exponential( 1 );
    // sigma_k_ID ~ exponential( 1 ); // commented out
    

// LIKELIHOODS

    // 
    for ( i in 1:N ) {
        p[i] = a
           + v_i__j_kt_ID [ikt_ID[i],1]
           + v_i__j_kt_ID [jkt_ID[i],2]
           + v_ij_kt_ID [ijkt_ID[i]] * sigma_ij_kt_ID
           + v_i2j_kt_ID [i2jkt_ID[i]] * sigma_i2j_kt_ID
           // + v_k_ID [k_ID[i]] * sigma_k_ID // commented out
           ;
          }

          
          // update targets
          y ~ binomial_logit( denominator, p);

}
          
generated quantities{
    matrix[2,2] Rho_i__j_kt_ID;

    real<lower=0> variance_i_kt;
    real<lower=0> variance_j_kt;
    
    Rho_i__j_kt_ID = L_Rho_i__j_kt_ID * L_Rho_i__j_kt_ID';

    variance_i_kt = sigma_i__j_kt_ID[1]^2;
    variance_j_kt = sigma_i__j_kt_ID[2]^2;
}
"
