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
  matrix[n_var, n_var] Sigma; // known variance
}

parameters {
  vector[n_var] beta;
  vector[n_var] beta_exp[n_exp];
}


model {
  
  beta ~ normal(3, 1);
  
  for(e in 1:n_exp){
    beta_exp[e] ~ multi_normal(beta, Sigma);
  }
  
  for(s in 1:S){
    vector[J] utilities;
    utilities = X[start[s]:end[s]]*beta_exp[experiment[s]];
    target += sum(log_softmax(utilities) .* Ny[start[s]:end[s]]);
  }
  
}
