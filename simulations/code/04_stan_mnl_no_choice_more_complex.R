library(MASS)
library(tidyverse)
library(mlogit)
library(opdesmixr)
library(rstan)
# library(here)

source("utils.R")


q = 5
J = 2
S = 25

# beta_vec_01_mixture = c(2, 4, 1, 3)
beta_vec_01_mixture = c(2, 2, 2, 2)
beta_no_choice = 1

sd_individual = 0.00001
Sigma = matrix(0.0, ncol = length(beta_vec_01_mixture) + 1, nrow = length(beta_vec_01_mixture) + 1)
diag(Sigma) = sd_individual





# Create an I-optimal design
design_array_i_opt_01_opdesmixr = opdesmixr::mnl_mixture_coord_exch(
  q = q, J = J, S = S, 
  beta = beta_vec_01_mixture, seed = 2021, order = 1, opt_crit = "I", transform_beta = F,
  n_random_starts = 1, n_cores = 1, max_it = 2)


# Simulate data for 30 heterogeneous respondents
simulated_data_i_opt_01 = lapply(1:30, function(r){
  print(r)
  
  real_beta_r = mvrnorm(1, c(beta_vec_01_mixture, beta_no_choice), Sigma)
  
  S = dim(design_array_i_opt_01_opdesmixr$X)[3]
  
  sim_data_r = simulate_mixture_choice_data_no_choice(
    design_array_i_opt_01_opdesmixr$X, 
    beta_vector_mixture = real_beta_r[1:(q-1)],
    beta_no_choice = real_beta_r[q],
    order = 1, 
    append_model_matrix = T)
  
  out = sim_data_r %>% 
    mutate(point = choice_set) %>% 
    mutate(
      choice_set = as.character(as.numeric(choice_set) + S*(r - 1)),
      resp = r
    ) %>% 
    bind_cols(
      matrix(rep(real_beta_r, (J+1)*S), ncol = length(real_beta_r), byrow = T) %>% 
        as_tibble() %>% 
        set_names(paste0("real_beta_resp_", 1:ncol(.)))
    )
  
  return(out)
  
}) %>% 
  bind_rows()





# Transform data to be used in mlogit
sim_data_mlogit_i_opt_01 <- dfidx(
  simulated_data_i_opt_01 %>% 
    select(starts_with("comp"), no_choice_ind, choice_set, choice), 
  idx = c("choice_set"))


# Fit a frequentist multinomial logit with mlogit
model_i_opt_01 = mlogit(choice ~ -1 + 
                          comp_1 + comp_2 + comp_3 + comp_4 + no_choice_ind,
                        sim_data_mlogit_i_opt_01)
model_i_opt_01








X_stan = lapply(split(simulated_data_i_opt_01, simulated_data_i_opt_01$choice_set), function(x){
  x %>% 
    select(starts_with("model_mat_col")) %>% 
    as.matrix()
})


Y_stan = sapply(split(simulated_data_i_opt_01, simulated_data_i_opt_01$choice_set), function(x){
  choice_vec = x %>% 
    pull(choice)
  
  return(which(choice_vec == 1))
})



stan_data_01 <- list(J = J+1, N = length(Y_stan), 
                     n_var = length(beta_vec_01_mixture) + 1, 
                     Y = Y_stan, X = X_stan)


m_01 = "
data {
    int<lower=2> J; // number of alternatives/outcomes
    int<lower=1> N; // of observations
    int<lower=1> n_var; // of covariates
    int<lower=0,upper=J> Y[N];
    matrix[J,n_var] X[N];
}

parameters {
    vector[n_var] beta;  // attribute effects 
}

model {
    for (i in 1:N)
        Y[i] ~ categorical_logit(X[i]*beta);
}
"

model_stan_01 <- stan(model_code = m_01, data = stan_data_01, iter = 2000, chains = 4, seed = 2023)
model_stan_01









# Grouped counts
simulated_data_i_opt_02 = simulated_data_i_opt_01 %>% 
  group_by(point, across(starts_with("model_mat"))) %>% 
  summarize(Ny = sum(choice)) %>% 
  rename(choice_set = point) %>% 
  ungroup()

X_stan_02 = simulated_data_i_opt_02 %>% 
  select(starts_with("model_mat")) %>% 
  as.matrix()

index_dataframe = lapply(unique(simulated_data_i_opt_02$choice_set), function(i){
  indices_i = which(simulated_data_i_opt_02$choice_set == i)
  out = tibble(
    choice_set = i, start = indices_i[1], end = indices_i[length(indices_i)]
  )
}) %>% 
  bind_rows()


stan_data_02 <- list(
  J = J+1, 
  S = max(index_dataframe$choice_set), 
  n_var = ncol(X_stan_02), 
  n_obs = nrow(X_stan_02),
  Ny = simulated_data_i_opt_02$Ny,
  X = X_stan_02,
  start = index_dataframe$start, # the starting observation for each decision/choice set
  end = index_dataframe$end  # the ending observation for each decision/choice set
  )


m_02 = "
data {
    int<lower=2> J; // number of alternatives/outcomes
    int<lower=2> S; // number of choice sets
    int<lower=1> n_var; // of covariates
    int<lower=1> n_obs; // of observations
    vector[n_obs] Ny; // counts of the response variable
    matrix[n_obs, n_var] X; // attribute matrix
    int start[S]; // the starting observation for each decision/choice set
    int end[S]; // the ending observation for each decision/choice set
}

parameters {
    vector[n_var] beta;  // attribute effects 
}

model {
    
    for(s in 1:S){
    
      vector[J] utilities;
      
      utilities = X[start[s]:end[s]]*beta;
      
      target += sum(log_softmax(utilities) .* Ny[start[s]:end[s]]);
    }
  
}
"

model_stan_02 <- stan(model_code = m_02, data = stan_data_02, iter = 2000, chains = 4, seed = 2023)
model_stan_02


