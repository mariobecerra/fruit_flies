// The input data
data {
  int<lower=2> J; // number of alternatives/outcomes
  int<lower=2> S; // number of choice sets
  int<lower=1> n_var; // of covariates
  int<lower=1> n_obs; // of observations
  int<lower=1> n_exp; // Number of experiments
  vector[n_obs] Ny; // counts of the response variable -> a vector of length 'n_obs'
  matrix[n_obs, n_var] X; // attribute matrix
  int start[S]; // the starting observation for each decision/choice set
  int end[S]; // the ending observation for each decision/choice set
  int experiment[S]; // mapping for each experiment
}

// The parameters to be estimated.
parameters {
  vector[n_var] beta;
  vector[n_var] beta_exp[n_exp];
  cov_matrix[n_var] Sigma;
}

// The model to be estimated.
model {
  matrix[n_var, n_var] Vprior;
  //Vprior = (n_var+1) * diag_matrix(rep_vector(1, n_var));
  // A more informative prior
  Vprior = (0.9) * diag_matrix(rep_vector(1, n_var));
  
  beta ~ normal(0, 10);
  
  // IW becomes more informative when degrees of freedom increase and values in the scale matrix become smaller
  // (from https://mplus.sites.uu.nl/wp-content/uploads/sites/24/2012/07/mplusiw.pdf)
  // Sigma ~ inv_wishart((n_var + 1), Vprior);
  Sigma ~ inv_wishart((n_var + 1), Vprior);
  
  for(e in 1:n_exp){
    beta_exp[e] ~ multi_normal(beta, Sigma);
  }
  
  
  for(s in 1:S){
    vector[J] utilities;
    utilities = X[start[s]:end[s]]*beta_exp[experiment[s]];
    target += sum(log_softmax(utilities) .* Ny[start[s]:end[s]]); // elementwaise multiplication
  }

}

