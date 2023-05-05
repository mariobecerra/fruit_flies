library(rstan)
library(tidyverse)
library(ggtern)
library(here)

simulated_data = readRDS(here("out/sims_out/model_2_levels_sim_simulated_data_tibble_2023-05-03.rds"))
model_stan = readRDS(here("out/sims_out/model_2_levels_sim_stan_object_2023-05-04.rds"))



create_lattice_design = function(n_var = 3, n_levels = 5){
  # Example: create_lattice_design(3, 10)
  
  # create_lattice_design(q, 20) %>% 
  #   ggtern::ggtern(ggtern::aes(x1, x2, x3)) +
  #   geom_point(shape = "x", size = 4) +
  #   theme_minimal() +
  #   ggtern::theme_nomask()
  
  lattice_df = lapply(1:n_var, function(i){
    tibble(x = seq(0, 1, by = (1/n_levels))) %>% 
      set_names(paste0("x", i))
  }) %>% 
    bind_cols() %>% 
    expand.grid() %>% 
    as_tibble() %>% 
    slice(which(abs(apply(., 1, sum) - 1) < 1e-12))
  
  return(lattice_df)
}




lattice_design_3_ing = create_lattice_design(n_var = 3, n_levels = 20)

ggtern(data = lattice_design_3_ing) +
  geom_point(aes(x1, x2, x3)) +
  theme_bw() +
  theme_nomask()


tibble_to_predict = expand_grid(
  lattice_design_3_ing, 
  var_1 = c(1, 0),
  var_2= c(1, 0),
  no_choice = 0
  )

model_matrix_to_predict = tibble_to_predict %>% 
  select(x1, x2, no_choice, var_1, var_2) %>% 
  as.matrix()


beta_level_0_draws = rstan::extract(model_stan, "beta_level_0")
Sigma_level_0_draws = rstan::extract(model_stan, "Sigma_level_0")

Sigma_level_1_draws = rstan::extract(model_stan, "Sigma_level_1")

alpha_draws = rstan::extract(model_stan, "alpha")

n_draws = nrow(beta_level_0_draws[[1]])


# Proper predictive
utilities_post = matrix(0.0, ncol = n_draws, nrow = nrow(model_matrix_to_predict))

for(d in 1:ncol(utilities_post)){
  
  # I don't know if I have to sample from a Normal here, I think I don't because there's no variance.
  alpha_vec_d = alpha_draws[[1]][d, ]
  
  
  Sigma_level_0_d = Sigma_level_0_draws[[1]][d, , ]
  beta_level_0_d = beta_level_0_draws[[1]][d, ]
  beta_level_1 = mvrnorm(1, beta_level_0_d, Sigma_level_0_d)
  
  Sigma_level_1 = Sigma_level_1_draws[[1]][d, , ]
  beta_level_2 = mvrnorm(1, beta_level_1, Sigma_level_1)
  
  param_vec = c(beta_level_2, alpha_vec_d)
  
  utilities_d = model_matrix_to_predict %*% param_vec
  
  utilities_post[, d] = utilities_d
  
}
  

quantiles_utilities = t(apply(utilities_post, 1, function(x) c(quantile(x, probs = c(0.1, 0.5, 0.9)), sd = sd(x)) ))



lattice_utilities = tibble_to_predict %>% 
  mutate(order_original = 1:n()) %>% 
  bind_cols(quantiles_utilities %>% 
              as_tibble() %>% 
              set_names(c("p10", "p50", "p90", "sd"))) %>% 
  group_by(var_1, var_2) %>% 
  arrange(desc(p50)) %>% 
  ungroup() %>% 
  mutate(ix1 = 1:n())










# plug-in

utilities_plugin = matrix(0.0, ncol = n_draws, nrow = nrow(model_matrix_to_predict))

for(d in 1:ncol(utilities_plugin)){
  
  alpha_vec_d = alpha_draws[[1]][d, ]
  beta_level_0_d = beta_level_0_draws[[1]][d, ]
  
  param_vec = c(beta_level_0_d, alpha_vec_d)
  
  utilities_d = model_matrix_to_predict %*% param_vec
  
  utilities_plugin[, d] = utilities_d
  
}


quantiles_utilities_plugin = t(apply(utilities_plugin, 1, function(x) c(quantile(x, probs = c(0.1, 0.5, 0.9)), sd = sd(x)) ))



lattice_utilities_2 = lattice_utilities %>% 
  arrange(order_original) %>% 
  bind_cols(quantiles_utilities_plugin %>% 
              as_tibble() %>% 
              set_names(c("p10_pl", "p50_pl", "p90_pl", "sd_pl"))) %>% 
  group_by(var_1, var_2) %>% 
  arrange(desc(p50_pl)) %>% 
  ungroup() %>% 
  mutate(ix2 = 1:n())

lattice_utilities_2 %>% 
  filter(ix1 != ix2)


lattice_utilities_2 %>% 
  group_by(var_1, var_2) %>% 
  top_n(1, p50_pl) %>% 
  as.data.frame()



utilities_all_mean = get_posterior_mean(model_stan, "utilities_all")
utilities_all_mean[, ncol(utilities_all_mean)] %>% head()

simulated_data %>% 
  mutate(utility_model = utilities_all_mean[, ncol(utilities_all_mean)]) %>% 
  filter(no_choice_ind == 0) %>% 
  select(comp_1, comp_2, comp_3, other_vars_1, other_vars_2, utility_model) %>% 
  distinct() %>% 
  group_by(other_vars_1, other_vars_2) %>% 
  top_n(1, utility_model) %>% 
  as.data.frame()

  

