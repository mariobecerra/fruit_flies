data {
  int<lower=2> n_alt; // number of alternatives/outcomes
  int<lower=1> n_var; // of covariates
  int<lower=1> n_obs; // of observations
  int<lower=1> n_exp; // Number of experiments
  vector[n_obs] Ny; // counts of the response variable
  matrix[n_obs, n_var] X; // attribute matrix
  int n_images;
  int start[n_images]; // the starting observation for each image
  int end[n_images]; // the ending observation for each image
  int exp_index[n_images];

}

parameters {
  vector[n_var] beta_level_0;
  vector[n_var] beta_level_1[n_exp];
  corr_matrix[n_var] Omega_level_0;
  vector<lower=0>[n_var] tau_level_0;
}

transformed parameters {
  cov_matrix[n_var] Sigma_level_0 = quad_form_diag(Omega_level_0, tau_level_0);
  
}

model {
  vector[n_alt] utilities;

  
  // In LKJ:
  // if eta = 1, then the density is uniform over correlation matrices of a given order K (the number of row/column).
  // if eta > 1, the correlation values in correlation matrices are going to centered around 0. higher eta indicate no correlations (converge to identity correlation matrix).
  // https://yingqijing.medium.com/lkj-correlation-distribution-in-stan-29927b69e9be
  
  tau_level_0 ~ normal(1, 0.5); 
  Omega_level_0 ~ lkj_corr(5);
  
  beta_level_0 ~ normal(0, 2);
  
  for(i in 1:n_exp){
    beta_level_1[i] ~ multi_normal(beta_level_0, Sigma_level_0);
  }
  

  for(i in 1:n_images){
    // utilities = X[start[i]:end[i], 1:(n_mixture_cols+1)]*beta_level_1[exp_index[i]] + X[start[i]:end[i], (n_mixture_cols+2):n_var]*alpha;
    utilities = X[start[i]:end[i], 1:n_var]*beta_level_1[exp_index[i]];
    
    target += sum(log_softmax(utilities) .* Ny[start[i]:end[i]]);
  }
  
  
}

 
// // Save log-likelihood for LOO package
// generated quantities{
//   vector[n_images] log_lik;
// 
//   for(i in 1:n_images){
//       log_lik[i] = sum(log_softmax(X[start[i]:end[i]]*beta_level_1[exp_index[i]]) .* Ny[start[i]:end[i]]);
//     }
// }
// 

