library(rstan)
library(MASS)
library(tidyverse)
library(ggtern)
library(here)

# counts_european = read_csv(here("out/counts_choice_format.csv")) %>% 
#   filter(experiment != 1) %>% 
#   mutate(experiment = experiment - 1)

model_stan = readRDS(here("out/european_flies_hmnl_1_level_stan_object.rds"))



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



mixture_variable_names_no_UV = c('R','O','Y','G','B','P')
mixture_variable_names = c(mixture_variable_names_no_UV, 'UV')

lattice_design_7_ing = create_lattice_design(n_var = 7, n_levels = 5) %>% 
  set_names(mixture_variable_names)




tibble_to_predict = expand_grid(
  lattice_design_7_ing, 
  intensity = seq(-1, 2, by = 0.2),
  is_right = c(1, 0),
  is_daytime= c(1, 0),
  no_choice = 0
)






mixture_pairwise_interaction_names = c(
  'R*O','R*Y','R*G','R*B','R*P','R*UV',
  'O*Y','O*G','O*B','O*P','O*UV',
  'Y*G','Y*B','Y*P','Y*UV',
  'G*B','G*P','G*UV',
  'B*P','B*UV',
  'P*UV')

mixture_intensity_interaction_names = c('R*intensity','O*intensity','Y*intensity','G*intensity',
                                        'B*intensity','P*intensity','UV*intensity')


X_other_vars_nc_to_predict <- tibble_to_predict %>% 
  select(all_of(c('no_choice','is_right', 'is_daytime'))) %>%
  as.data.frame()

X_main_to_predict <- tibble_to_predict %>% 
  select(all_of(c(mixture_variable_names,'intensity'))) %>%
  as.data.frame()


X_color_pairwise_to_predict <- combn(mixture_variable_names, 2, function(x) X_main_to_predict[x[1]]*X_main_to_predict[x[2]], simplify = FALSE) %>%
  bind_cols() %>%
  set_names(mixture_pairwise_interaction_names)


X_color_intensity_to_predict <- X_main_to_predict %>%
  select(c(mixture_variable_names,'intensity')) %>%
  mutate(across(c(mixture_variable_names,'intensity'), ~.*intensity)) %>%
  set_names(c(mixture_intensity_interaction_names,'intensity^2'))


# X_stan <- X_main %>% 
#   select(-UV) %>% 
#   as.matrix()



X_to_predict <- X_main_to_predict %>% 
  select(-UV, -intensity) %>% 
  bind_cols(X_color_pairwise_to_predict) %>% 
  bind_cols(X_color_intensity_to_predict) %>%
  bind_cols(X_other_vars_nc_to_predict) %>% 
  as.matrix()




beta_level_0_draws = rstan::extract(model_stan, "beta_level_0")
Sigma_level_0_draws = rstan::extract(model_stan, "Sigma_level_0")

Sigma_level_1_draws = rstan::extract(model_stan, "Sigma_level_1")

alpha_draws = rstan::extract(model_stan, "alpha")

n_draws = nrow(beta_level_0_draws[[1]])


# Proper predictive
utilities_post = matrix(0.0, ncol = n_draws, nrow = nrow(X_to_predict))

for(d in 1:ncol(utilities_post)){
  
  # I don't know if I have to sample from a Normal here, I think I don't because there's no variance.
  alpha_vec_d = alpha_draws[[1]][d, ]
  
  
  Sigma_level_0_d = Sigma_level_0_draws[[1]][d, , ]
  beta_level_0_d = beta_level_0_draws[[1]][d, ]
  beta_level_1 = mvrnorm(1, beta_level_0_d, Sigma_level_0_d)
  
  param_vec = c(beta_level_1, alpha_vec_d)
  
  utilities_d = X_to_predict %*% param_vec
  
  utilities_post[, d] = utilities_d
  
}


quantiles_utilities = t(apply(utilities_post, 1, function(x) c(quantile(x, probs = c(0.1, 0.5, 0.9)), sd = sd(x)) ))



lattice_utilities = tibble_to_predict %>% 
  mutate(order_original = 1:n()) %>% 
  bind_cols(quantiles_utilities %>% 
              as_tibble() %>% 
              set_names(c("p10", "p50", "p90", "sd"))) %>% 
  group_by(is_right, is_daytime) %>% 
  arrange(desc(p50)) %>% 
  ungroup() %>% 
  mutate(ix1 = 1:n())










# plug-in

utilities_plugin = matrix(0.0, ncol = n_draws, nrow = nrow(X_to_predict))

for(d in 1:ncol(utilities_plugin)){
  
  alpha_vec_d = alpha_draws[[1]][d, ]
  beta_level_0_d = beta_level_0_draws[[1]][d, ]
  
  param_vec = c(beta_level_0_d, alpha_vec_d)
  
  utilities_d = X_to_predict %*% param_vec
  
  utilities_plugin[, d] = utilities_d
  
}


quantiles_utilities_plugin = t(apply(utilities_plugin, 1, function(x) c(quantile(x, probs = c(0.1, 0.5, 0.9)), sd = sd(x)) ))



lattice_utilities_2 = lattice_utilities %>% 
  arrange(order_original) %>% 
  bind_cols(quantiles_utilities_plugin %>% 
              as_tibble() %>% 
              set_names(c("p10_pl", "p50_pl", "p90_pl", "sd_pl"))) %>% 
  group_by(is_right, is_daytime) %>% 
  arrange(desc(p50_pl)) %>% 
  ungroup() %>% 
  mutate(ix2 = 1:n())

lattice_utilities_2 %>% 
  filter(ix1 != ix2)


lattice_utilities_2 %>% 
  group_by(is_right, is_daytime) %>% 
  top_n(1, p50_pl) %>% 
  as.data.frame()


lattice_utilities_2 %>% 
  group_by(is_right, is_daytime) %>% 
  top_n(3, p50) %>% 
  as.data.frame()




lattice_utilities_2 %>% 
  filter(intensity > 0) %>% 
  group_by(is_right, is_daytime) %>% 
  top_n(3, p50) %>% 
  as.data.frame()
