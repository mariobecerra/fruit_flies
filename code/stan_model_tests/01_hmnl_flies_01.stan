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
  corr_matrix[n_var] Omega;
  vector<lower=0>[n_var] tau;
}

transformed parameters {
  cov_matrix[n_var] Sigma = quad_form_diag(Omega, tau);
  // quad_form_diag(matrix Omega, vector tau) returns the quadratic form, i.e., diag_matrix(tau) * Omega * diag_matrix(tau).
  // Sigma = diag_matrix(tau) * Omega * diag_matrix(tau)
}

model {
  
  // tau ~ cauchy(0, 2.5); 
  tau ~ normal(1, 3); // Seems to be a tight prior that only led to 2 divergent transitions
  Omega ~ lkj_corr(1);
  
  beta ~ normal(0, 1);
  
  for(e in 1:n_exp){
    beta_exp[e] ~ multi_normal(beta, Sigma);
  }
  
  for(s in 1:S){
    vector[J] utilities;
    utilities = X[start[s]:end[s]]*beta_exp[experiment[s]];
    target += sum(log_softmax(utilities) .* Ny[start[s]:end[s]]);
  }
  
}
