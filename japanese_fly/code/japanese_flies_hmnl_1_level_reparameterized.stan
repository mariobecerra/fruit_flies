data {
  int<lower=2> n_alt; // number of alternatives/outcomes
  int<lower=1> n_var; // of covariates
  int<lower=1> n_obs; // of observations
  int<lower=1> n_exp; // Number of experiments
  vector[n_obs] Ny; // counts of the response variable
  matrix[n_obs, n_var] X; // attribute matrix
  int n_images; // Number of images
  int start[n_images]; // the starting observation for each image
  int end[n_images]; // the ending observation for each image
  int exp_index[n_images];
  
  // priors
  real tau_level_0_mean;
  real<lower=0> tau_level_0_sd;
  real L_Omega_level_0_param;
  vector[n_var] beta_level_0_mean;
  vector[n_var] beta_level_0_sd;
}

parameters {
  vector[n_var] beta_level_0;
  matrix[n_exp, n_var] z;
  // corr_matrix[n_var] Omega_level_0;
  cholesky_factor_corr[n_var] L_Omega_level_0;
  vector<lower=0>[n_var] tau_level_0;
}

transformed parameters {
  // reparametrization of beta
  matrix[n_exp, n_var] beta_level_1 = rep_matrix(beta_level_0', n_exp) + z*diag_pre_multiply(tau_level_0, L_Omega_level_0);
}

model {
  vector[n_alt] utilities;

  to_vector(z) ~ normal(0, 1);
  // In LKJ:
  // if eta = 1, then the density is uniform over correlation matrices of a given order K (the number of row/column).
  // if eta > 1, the correlation values in correlation matrices are going to centered around 0. higher eta indicate no correlations (converge to identity correlation matrix).
  // https://yingqijing.medium.com/lkj-correlation-distribution-in-stan-29927b69e9be
  
  tau_level_0 ~ normal(tau_level_0_mean, tau_level_0_sd); 
  L_Omega_level_0 ~ lkj_corr_cholesky(L_Omega_level_0_param);
  
  beta_level_0 ~ normal(beta_level_0_mean, beta_level_0_sd);
  

  for(i in 1:n_images){
    // print("loop iteration: ", i);
    utilities = X[start[i]:end[i], 1:n_var]*beta_level_1[exp_index[i]]';
    target += sum(log_softmax(utilities) .* Ny[start[i]:end[i]]);
  }
}


generated quantities{
  matrix[n_var, n_var] Sigma_level_0 = diag_pre_multiply(tau_level_0, L_Omega_level_0) * diag_pre_multiply(tau_level_0, L_Omega_level_0)';
  
  // // Save log-likelihood for LOO package
  //   vector[n_images] log_lik;
  // 
  //   for(i in 1:n_images){
  //       log_lik[i] = sum(log_softmax(X[start[i]:end[i]]*beta_level_1[exp_index[i]]) .* Ny[start[i]:end[i]]);
  //     }
  
}
 
 
