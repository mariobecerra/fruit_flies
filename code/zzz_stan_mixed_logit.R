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

sd_individual = 0.5
Sigma = matrix(0.0, ncol = length(beta_vec_01), nrow = length(beta_vec_01))
diag(Sigma) = sd_individual


beta_vec_01 = rep(5, 4)


# Create an I-optimal design
design_array_i_opt_01_opdesmixr = opdesmixr::mnl_mixture_coord_exch(
  q = q, J = J, S = S, 
  beta = beta_vec_01, seed = 2021, order = 1, opt_crit = "I", transform_beta = F,
  n_random_starts = 1, n_cores = 1, max_it = 2)


# Repeat with 30 respondents
simulated_data_i_opt_01 = lapply(1:30, function(r){
  print(r)
  
  real_beta_r = mvrnorm(1, beta_vec_01, Sigma)
  
  S = dim(design_array_i_opt_01_opdesmixr$X)[3]
  
  sim_data_r = simulate_mixture_choice_data(
    design_array_i_opt_01_opdesmixr$X, 
    beta = real_beta_r, 
    order = 1, 
    append_model_matrix = F)
  
  out = design_array_i_opt_01_opdesmixr$X %>%
    mnl_design_array_to_dataframe() %>% 
    set_names(c(paste0("comp_", 1:(ncol(.)-1)), "choice_set")) %>% 
    mutate(point = choice_set) %>% 
    mutate(
      choice_set = as.character(as.numeric(choice_set) + S*(r - 1)),
      resp = r
    ) %>% 
    mutate(choice = sim_data_r$choice) %>% 
    bind_cols(
      matrix(rep(real_beta_r, 2*S), ncol = length(real_beta_r), byrow = T) %>% 
        as_tibble() %>% 
        set_names(paste0("real_beta_resp_", 1:ncol(.)))
    )
  
  return(out)
  
}) %>% 
  bind_rows()





# Transform data to be used in mlogit
sim_data_mlogit_i_opt_01 <- dfidx(
  simulated_data_i_opt_01 %>% 
    select(starts_with("comp"), choice_set, choice), 
  idx = c("choice_set"))


# Fit a multinomial logit
model_i_opt_01 = mlogit(choice ~ -1 + 
                          comp_1 + comp_2 + comp_3 + comp_4,
                        sim_data_mlogit_i_opt_01)
model_i_opt_01


# Transform data for a MNL in Stan

X_stan_01 = lapply(split(simulated_data_i_opt_01, simulated_data_i_opt_01$choice_set), function(x){
  x %>% 
    select(starts_with("comp")) %>% 
    select(-comp_5) %>% # remove last component
    as.matrix()
})


Y_stan_01 = sapply(split(simulated_data_i_opt_01, simulated_data_i_opt_01$choice_set), function(x){
  first_opt = x %>% 
    slice(1) %>% 
    pull(choice)
  
  return(as.numeric(first_opt == 0) + 1)
})



stan_data_01 <- list(J = J, N = length(Y_stan_01), n_var = length(beta_vec_01), Y = Y_stan_01, X = X_stan_01)


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









### Do mixed logit here


