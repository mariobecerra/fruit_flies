data {
  int<lower=2> n_alt; // number of alternatives/outcomes
  int<lower=1> n_var; // of covariates
  int<lower=1> n_obs; // of observations
  int<lower=1> n_exp; // Number of experiments
  int<lower=1> n_betas_level_2; // Number of level 2 betas
  int<lower=1> n_betas_level_3; // Number of level 3 betas
  vector[n_obs] Ny; // counts of the response variable
  matrix[n_obs, n_var] X; // attribute matrix
  int start[n_betas_level_3]; // the starting observation for each image
  int end[n_betas_level_3]; // the ending observation for each image
  int experiment[n_betas_level_3]; // mapping for each experiment
  int choice_set[n_betas_level_3]; // mapping for each choice set
  int image[n_betas_level_3]; // mapping for each experiment
}

parameters {
  vector[n_var] beta_level_0;
  vector[n_var] beta_level_1[n_exp];
  vector[n_var] beta_level_2[n_betas_level_2];
  vector[n_var] beta_level_3[n_betas_level_3];
  corr_matrix[n_var] Omega_0;
  vector<lower=0>[n_var] tau_0;
  corr_matrix[n_var] Omega_1;
  vector<lower=0>[n_var] tau_1;
  corr_matrix[n_var] Omega_2;
  vector<lower=0>[n_var] tau_2;
}

transformed parameters {
  cov_matrix[n_var] Sigma_0 = quad_form_diag(Omega_0, tau_0);
  cov_matrix[n_var] Sigma_1 = quad_form_diag(Omega_1, tau_1);
  cov_matrix[n_var] Sigma_2 = quad_form_diag(Omega_2, tau_2);
  // quad_form_diag(matrix Omega, vector tau) returns the quadratic form, i.e., diag_matrix(tau) * Omega * diag_matrix(tau).
  // Sigma = diag_matrix(tau) * Omega * diag_matrix(tau)
}

model {
  
  // tau ~ cauchy(0, 2.5); 
  tau_0 ~ normal(1, 3); 
  Omega_0 ~ lkj_corr(1);
  
  tau_1 ~ normal(1, 3); 
  Omega_1 ~ lkj_corr(1);
  
  tau_2 ~ normal(1, 3); 
  Omega_2 ~ lkj_corr(1);
  
  beta_level_0 ~ normal(0, 1);
  
  
  for(i in 1:n_exp){
    beta_level_1[i] ~ multi_normal(beta_level_0, Sigma_0);
  }
  
  for(i in 1:n_betas_level_2){
    beta_level_2[i] ~ multi_normal(beta_level_1[experiment[i]], Sigma_1);
  }
  
  // for(i in 1:n_betas_level_3){
    //   beta_level_3[i] ~ multi_normal(beta_level_2, Sigma_2);
    //       
    //   vector[n_alt] utilities;
    //   utilities = X[start[i]:end[i]]*beta_level_3[image[i]];
    //   target += sum(log_softmax(utilities) .* Ny[start[i]:end[i]]);
    // }
    
    for(i in 1:n_betas_level_3){
      beta_level_3[i] ~ multi_normal(beta_level_2[choice_set[i]], Sigma_2);
    }
    
    for(i in 1:n_betas_level_3){
      vector[n_alt] utilities;
      utilities = X[start[i]:end[i]]*beta_level_3[image[i]];
      target += sum(log_softmax(utilities) .* Ny[start[i]:end[i]]);
    }
    
    
}

