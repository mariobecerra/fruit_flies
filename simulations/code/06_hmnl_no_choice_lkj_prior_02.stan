data {
  int<lower=2> J; // number of alternatives/outcomes
  int<lower=2> S; // number of choice sets
  int<lower=1> n_var; // of covariates
  int<lower=1> n_obs; // of observations
  int<lower=1> n_exp; // Number of experiments
  vector[n_obs] Ny; // counts of the response variable
  matrix[n_obs, n_var] X; // attribute matrix
  int start[S]; // the starting observation for each decision/choice set
  int end[S]; // the ending observation for each decision/choice set
  int experiment[S]; // mapping for each experiment
}



parameters {
  vector[n_var] beta;
  vector[n_var] beta_exp[n_exp];
  cholesky_factor_corr[n_var] L_Omega; // // cholesky of correlation matrix
  vector<lower=0>[n_var] tau;
}



transformed parameters {
  cholesky_factor_cov[n_var] L; // cholesky factors
  L = diag_pre_multiply(tau, L_Omega); // Returns the product of the diagonal matrix formed from the vector tau and the matrix L_Omega, i.e., diag_matrix(tau) * L_Omega.
  
}



model {
  tau ~ cauchy(0, 2.5); 
  L_Omega ~ lkj_corr_cholesky(1);
  
  beta ~ normal(3, 1);
  
  for(e in 1:n_exp){
    beta_exp[e] ~ multi_normal_cholesky(beta, L);
  }
  
  for(s in 1:S){
    vector[J] utilities;
    utilities = X[start[s]:end[s]]*beta_exp[experiment[s]];
    target += sum(log_softmax(utilities) .* Ny[start[s]:end[s]]);
  }
}



generated quantities {
  matrix[n_var, n_var] Sigma;
  corr_matrix[n_var] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  // multiply_lower_tri_self_transpose(L_Omega) returns the product of L_Omega by its transpose, i.e., L_Omega * t(L_Omega).
  
  Sigma = quad_form_diag(Omega, tau);
  // quad_form_diag(matrix Omega, vector tau) returns the quadratic form, i.e., diag_matrix(tau) * Omega * diag_matrix(tau).
  // Then:
  // Sigma = diag_matrix(tau) * (L_Omega * t(L_Omega)) * diag_matrix(tau)
  // Sigma = diag_matrix(tau) * Omega * diag_matrix(tau)
}

