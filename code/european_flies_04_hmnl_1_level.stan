data {
  int<lower=2> n_alt; // number of alternatives/outcomes
  int<lower=1> n_var; // of covariates
  int<lower=1> n_obs; // of observations
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
}

parameters {
  vector[n_mixture_cols+1] beta_level_0;
  vector[n_mixture_cols+1] beta_level_1[n_exp];
  // vector[n_mixture_cols+1] beta_level_2[n_betas_level_2];
  corr_matrix[n_mixture_cols+1] Omega_level_0;
  vector<lower=0>[n_mixture_cols+1] tau_level_0;
  corr_matrix[n_mixture_cols+1] Omega_level_1;
  vector<lower=0>[n_mixture_cols+1] tau_level_1;
  
  vector[n_extra_vars] alpha;
}

transformed parameters {
  cov_matrix[n_mixture_cols+1] Sigma_level_0 = quad_form_diag(Omega_level_0, tau_level_0);
  cov_matrix[n_mixture_cols+1] Sigma_level_1 = quad_form_diag(Omega_level_1, tau_level_1);
  // quad_form_diag(matrix Omega, vector tau) returns the quadratic form, i.e., diag_matrix(tau) * Omega * diag_matrix(tau).
  // Sigma = diag_matrix(tau) * Omega * diag_matrix(tau)
  
  vector[n_alt*n_images] utilities_all;
  
  // Save the utilities of each observation in a vector
  for(i in 1:n_images){
    utilities_all[start[i]:end[i]] = X[start[i]:end[i], 1:(n_mixture_cols+1)]*beta_level_1[exp_index[i]] + X[start[i]:end[i], (n_mixture_cols+2):n_var]*alpha;
  }
  
}

model {
  
  // vector[n_alt*n_images] utilities_all;
  // matrix[n_images, n_alt] utilities_all;
  
  // In LKJ:
  // if eta = 1, then the density is uniform over correlation matrices of a given order K (the number of row/column).
  // if eta > 1, the correlation values in correlation matrices are going to centered around 0. higher eta indicate no correlations (converge to identity correlation matrix).
  // https://yingqijing.medium.com/lkj-correlation-distribution-in-stan-29927b69e9be
  
  // tau ~ cauchy(0, 2.5); 
  
  // tau_level_0 ~ normal(0.77, 0.5); 
  tau_level_0 ~ normal(1, 0.01);   // tau_level_1 ~ normal(0.5, 0.1);   // beta_level_0 ~ normal(0, 2);  /
  Omega_level_0 ~ lkj_corr(100);
  
  
  tau_level_1 ~ normal(0.5, 0.1); 
  Omega_level_1 ~ lkj_corr(10);
  
  
  // beta_level_0 ~ normal(0, 2);
  beta_level_0 ~ normal([-2.99, -1.02, -0.25, -0.63, -2.54, -2.19, 5.41, 8.31, 2.64, 1.13, 15.71, 9.36, 3.25, 2.00, 3.65, 5.79, 2.75, -0.68, 7.72, -6.42, -4.48, 11.68, 8.30, 3.36, -7.40, 9.02, 6.10, 1.14, 3.09, 1.70, -0.21, -0.93, -0.32, -0.79, -0.04, 4.79], 2);
  
  alpha ~ normal(0, 2);
  

  for(i in 1:n_exp){
    beta_level_1[i] ~ multi_normal(beta_level_0, Sigma_level_0);
  }
  
  // for(i in 1:n_betas_level_2){
  //   beta_level_2[i] ~ multi_normal(beta_level_1[experiment[i]], Sigma_level_1);
  // }
  
  
  for(i in 1:n_images){
    target += sum(log_softmax(utilities_all[start[i]:end[i]]) .* Ny[start[i]:end[i]]);
  }
  
}



