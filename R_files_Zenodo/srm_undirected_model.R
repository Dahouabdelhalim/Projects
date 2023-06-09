model_code_0 <- "
data{
    int N;

    int N_i__j_kt_ID;
    int N_ij_kt_ID;
    // int N_k_ID; // commented out
    
    int y[N];
    int denominator[N];
    
    int ikt_ID[N];
    int jkt_ID[N];
    
    int ijkt_ID[N];

    // int k_ID[N]; // commented out

}

parameters{

    real a;

    ////////// Variances and covariances
    
    // Node-level
    real <lower=0> sigma_i__j_kt_ID;
    
    // Symmetric dyads
    real <lower=0> sigma_ij_kt_ID;

    // Group
    // real <lower=0> sigma_k_ID; // Commented out

    /////////// Effects
    vector[N_i__j_kt_ID] v_i__j_kt_ID; // Node-level varying intercepts
    vector[N_ij_kt_ID] v_ij_kt_ID; // Symmetric dyad effects
    // vector[N_k_ID] v_k_ID; // Group effects // Commented out

}


model{

    vector[N] p;

    //priors for fixed effects
    a ~ normal( -4, 5 );

    // priors on varying effects
    to_vector(v_i__j_kt_ID) ~ normal( 0,1 );
    to_vector(v_ij_kt_ID) ~ normal( 0,1 );
    // to_vector(v_k_ID) ~ normal( 0,1 ); // Commented out
    
    // Sigmas and covariance
    sigma_i__j_kt_ID ~ exponential( 1 );
    sigma_ij_kt_ID ~ exponential( 1 );
    // sigma_k_ID ~ exponential( 1 ); // Commented out
    

// LIKELIHOODS

    // Undirected outcome
    for ( i in 1:N ) {
        p[i] = a
           + v_i__j_kt_ID [ikt_ID[i]] * sigma_i__j_kt_ID
           + v_i__j_kt_ID [jkt_ID[i]] * sigma_i__j_kt_ID
           + v_ij_kt_ID [ijkt_ID[i]] * sigma_ij_kt_ID
           // + v_k_ID [k_ID[i]] * sigma_k_ID // Commented out
           ;
          }

          
          // update targets
          y ~ binomial_logit( denominator, p);

}
"
