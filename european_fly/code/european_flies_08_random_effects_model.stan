data {
  int<lower=2> n_alt; // number of alternatives/outcomes
  int<lower=1> n_var; // number of covariates
  int<lower=1> n_obs; // number of observations
  int<lower=1> n_exp; // Number of experiments
  // int<lower=1> n_betas_level_2; // Number of level 2 betas
  vector[n_obs] Ny; // counts of the response variable
  matrix[n_obs, n_var] X; // attribute matrix
  // int experiment[n_betas_level_2]; // mapping for each experiment
  // int choice_set[n_betas_level_2]; // mapping for each choice set
  int n_images;
  int start[n_images]; // the starting observation for each image
  int end[n_images]; // the ending observation for each image
  // int experiment_index_2[n_images];
  int exp_index[n_images];
  
  int n_mixture_cols;
  int n_extra_vars;
  
  matrix[n_obs, n_exp] Z; 
}

parameters {
  vector[n_mixture_cols+1] beta_level_0;
  
  vector[n_extra_vars] alpha;
  
  corr_matrix[n_exp] Omega_level_0;
  vector<lower=0>[n_exp] tau_level_0;
  
  vector[n_exp] u;
  
}

transformed parameters {
  
  cov_matrix[n_exp] Sigma_level_0 = quad_form_diag(Omega_level_0, tau_level_0);
  
  vector[n_alt*n_images] utilities_all;
  
  // Save the utilities of each observation in a vector
  for(i in 1:n_images){
    // J dimension vector
    utilities_all[start[i]:end[i]] = 
       X[start[i]:end[i], 1:(n_mixture_cols+1)] * beta_level_0 + // mixture and process variables
       X[start[i]:end[i], (n_mixture_cols+2):n_var] * alpha + // is_right
       Z[start[i]:end[i], 1:n_exp] * u; // forgot how to get the whole row of a matrix in Stan
       
  }
  
}

model {
  
  
  beta_level_0 ~ normal(0, 2);
  
  alpha ~ normal(0, 2);
  

  tau_level_0 ~ normal(1, 0.5);   // tau_level_1 ~ normal(0.5, 0.1);   // beta_level_0 ~ normal(0, 2);  /
  Omega_level_0 ~ lkj_corr(5);
  u ~ multi_normal(rep_vector(0, n_exp), Sigma_level_0);
  
  
  for(i in 1:n_images){
    target += sum(log_softmax(utilities_all[start[i]:end[i]]) .* Ny[start[i]:end[i]]);
  }
  
}



