library(tidyverse)
library(mlogit)
library(opdesmixr)
library(rstan)
# library(here)

source("utils.R")


q = 5
J = 2
S = 60


beta_vec_01 = rep(5, 4)


# Create an I-optimal design
design_array_i_opt_01_opdesmixr = opdesmixr::mnl_mixture_coord_exch(
  q = q, J = J, S = S, 
  beta = beta_vec_01, seed = 2021, order = 1, opt_crit = "I", transform_beta = F,
  n_random_starts = 1, n_cores = 1, max_it = 2)


# Repeat with 50 respondents
design_array_i_opt_01_df = lapply(1:50, function(r){
  print(r)
  
  S = dim(design_array_i_opt_01_opdesmixr$X)[3]
  
  out = design_array_i_opt_01_opdesmixr$X %>%
    mnl_design_array_to_dataframe() %>% 
    mutate(point = choice_set) %>% 
    mutate(
      choice_set = as.character(as.numeric(choice_set) + S*(r - 1)),
      resp = r
    )
  
  return(out)
  
}) %>% 
  bind_rows()


design_array_i_opt_01 = mnl_design_dataframe_to_array(design_array_i_opt_01_df %>% select(-resp, -point))



simulated_data_i_opt_01 = simulate_mixture_choice_data(design_array_i_opt_01, beta = beta_vec_01, order = 1, append_model_matrix = F)



# Transform data to be used in mlogit
sim_data_mlogit_i_opt_01 <- dfidx(
  simulated_data_i_opt_01 %>% 
    select(starts_with("comp"), choice_set, choice), 
  idx = c("choice_set"))

model_i_opt_01 = mlogit(choice ~ -1 + 
                          comp_1 + comp_2 + comp_3 + comp_4,
                        sim_data_mlogit_i_opt_01)
model_i_opt_01



X_stan = lapply(split(simulated_data_i_opt_01, simulated_data_i_opt_01$choice_set), function(x){
  x %>% 
    select(starts_with("comp")) %>% 
    select(-comp_5) %>% # remove last component
    as.matrix()
})


Y_stan = sapply(split(simulated_data_i_opt_01, simulated_data_i_opt_01$choice_set), function(x){
  first_opt = x %>% 
     slice(1) %>% 
    pull(choice)
  
  return(as.numeric(first_opt == 0) + 1)
})



stan_data_01 <- list(J = J, N = length(Y_stan), n_var = 4, Y = Y_stan, X = X_stan)




m1 = "
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

p1 <- stan(model_code = m1, data = stan_data_01, iter = 2000, chains = 4, seed = 2023)





