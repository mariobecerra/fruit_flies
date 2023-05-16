library(rstan)
library(MASS)
library(tidyverse)
library(ggtern)
library(here)

source(here("code/european_flies_06_hmnl_1_level_utils.R"))

# counts_european = read_csv(here("out/counts_choice_format.csv")) %>% 
#   filter(experiment != 1) %>% 
#   mutate(experiment = experiment - 1)

model_stan = readRDS(here("out/european_flies_hmnl_1_level_stan_object_noisdaytime.rds"))



lattice_design_7_ing = create_lattice_design(n_var = 7, n_levels = 5) %>% 
  set_names(c('R','O','Y','G','B','P', 'UV'))




df_to_predict = expand_grid(
  lattice_design_7_ing, 
  intensity = seq(-1, 2, by = 0.2),
  is_right = c(1, 0),
  no_choice = 0
) %>% 
  as.data.frame()


X_to_predict_list <- create_model_matrix_second_order_scheffe(df_to_predict)





betas_level_0_summary = summary(model_stan, pars = c("beta_level_0", "alpha"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_to_predict_list$X),
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F))


betas_level_0_summary %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Hierarchical model", subtitle = "Betas level 0")






names_betas_level_1 = X_to_predict_list$names_betas_level_1


betas_level_1_summary_0 = summary(model_stan, pars = c("beta_level_1"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "2.5%", "10%", "90%", "97.5%")

n_experiments = nrow(betas_level_1_summary_0)/length(names_betas_level_1)

betas_level_1_summary = betas_level_1_summary_0 %>% 
  mutate(
    exp = paste0(
      "Exp ",
      rep(1:n_experiments, each = length(names_betas_level_1))
    ),
    variable = rep(names_betas_level_1, n_experiments)
  ) %>% 
  group_by(exp) %>% 
  mutate(ix = 1:n()) %>%
  ungroup() %>% 
  # mutate(variable = paste0(exp, ", ", var)) %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = T))



betas_level_1_summary %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  facet_wrap(~exp) +
  coord_flip() +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  ggtitle("Hierarchical model", subtitle = "Betas level 1")



betas_level_1_summary %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = F)) %>% 
  ggplot(aes(x = variable, color = exp)) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.6), size = 0.5) +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), size = 0.3, position = position_dodge(width = 0.6)) +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), size = 0.7, position = position_dodge(width = 0.6)) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  geom_vline(xintercept = 1:36 + 0.5, size = 0.3, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle("Hierarchical model", subtitle = "Betas level 1")







beta_level_0_draws = rstan::extract(model_stan, "beta_level_0")
Sigma_level_0_draws = rstan::extract(model_stan, "Sigma_level_0")

# Sigma_level_1_draws = rstan::extract(model_stan, "Sigma_level_1")

alpha_draws = rstan::extract(model_stan, "alpha")

n_draws = nrow(beta_level_0_draws[[1]])





# Proper predictive
utilities_post = matrix(0.0, ncol = n_draws, nrow = nrow(X_to_predict_list$X))

for(d in 1:ncol(utilities_post)){
  
  # I don't know if I have to sample from a Normal here, I think I don't because there's no variance.
  alpha_vec_d = alpha_draws[[1]][d, ]
  
  
  Sigma_level_0_d = Sigma_level_0_draws[[1]][d, , ]
  beta_level_0_d = beta_level_0_draws[[1]][d, ]
  beta_level_1 = mvrnorm(1, beta_level_0_d, Sigma_level_0_d)
  
  param_vec = c(beta_level_1, alpha_vec_d)
  
  utilities_d = X_to_predict_list$X %*% param_vec
  
  utilities_post[, d] = utilities_d
  
}


quantiles_utilities = t(apply(utilities_post, 1, function(x) c(quantile(x, probs = c(0.1, 0.5, 0.9)), sd = sd(x)) ))
rm(utilities_post)


lattice_utilities = df_to_predict %>% 
  mutate(order_original = 1:n()) %>% 
  bind_cols(quantiles_utilities %>% 
              as_tibble() %>% 
              set_names(c("p10", "p50", "p90", "sd"))) %>% 
  group_by(is_right) %>% 
  arrange(desc(p50)) %>% 
  ungroup() %>% 
  mutate(ix1 = 1:n())


lattice_utilities %>% 
  group_by(is_right) %>% 
  top_n(3, p50) %>% 
  as.data.frame()


lattice_utilities %>% 
  filter(intensity > 0) %>% 
  group_by(is_right) %>% 
  top_n(3, p50) %>% 
  as.data.frame()



















lattice_design_3_ing = create_lattice_design(n_var = 3, n_levels = 25) %>%
  set_names(c("R", "O", "UV"))



df_to_predict_2 = lattice_design_3_ing %>% 
  mutate(
    Y = 0, G = 0, B = 0, P = 0
  ) %>% 
  select(R, O, Y, G, B, P, UV) %>% 
  expand_grid(
    .,
    intensity = seq(-1, 2, by = 0.05), 
    is_right = c(1, 0),
    no_choice = 0
  )



X_to_predict_list_2 <- create_model_matrix_second_order_scheffe(df_to_predict_2)



# Proper predictive
utilities_post_2 = matrix(0.0, ncol = n_draws, nrow = nrow(X_to_predict_list_2$X))

for(d in 1:ncol(utilities_post_2)){
  
  # I don't know if I have to sample from a Normal here, I think I don't because there's no variance.
  alpha_vec_d = alpha_draws[[1]][d, ]
  
  
  Sigma_level_0_d = Sigma_level_0_draws[[1]][d, , ]
  beta_level_0_d = beta_level_0_draws[[1]][d, ]
  beta_level_1 = mvrnorm(1, beta_level_0_d, Sigma_level_0_d)
  
  param_vec = c(beta_level_1, alpha_vec_d)
  
  utilities_d = X_to_predict_list_2$X %*% param_vec
  
  utilities_post_2[, d] = utilities_d
  
}


quantiles_utilities_2 = t(apply(utilities_post_2, 1, function(x) c(quantile(x, probs = c(0.1, 0.5, 0.9)), sd = sd(x)) ))
rm(utilities_post_2)


lattice_utilities_2 = df_to_predict_2 %>% 
  mutate(order_original = 1:n()) %>% 
  bind_cols(quantiles_utilities_2 %>% 
              as_tibble() %>% 
              set_names(c("p10", "p50", "p90", "sd"))) %>% 
  group_by(is_right) %>% 
  arrange(desc(p50)) %>% 
  ungroup() %>% 
  mutate(ix1 = 1:n())


lattice_utilities_2 %>% 
  group_by(is_right) %>% 
  top_n(3, p50) %>% 
  as.data.frame()


lattice_utilities_2 %>% 
  filter(intensity > 0) %>% 
  group_by(is_right) %>% 
  top_n(3, p50) %>% 
  as.data.frame()






lattice_utilities_2 %>% 
  group_by(is_right) %>% 
  top_n(30, p50) %>% 
  as.data.frame() %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point() +
  geom_linerange(aes(ymin = p10, ymax = p90)) +
  theme_bw()






lattice_utilities_2 %>% 
  group_by(is_right) %>% 
  top_n(30, p50) %>% 
  bind_rows(
    lattice_utilities %>% 
      group_by(is_right) %>% 
      top_n(-30, p50)
  ) %>% 
  as.data.frame() %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point() +
  geom_linerange(aes(ymin = p10, ymax = p90)) +
  theme_bw()











# More granular around the best values predicted in previous step
lattice_design_3_ing_2 = create_lattice_design(n_var = 3, n_levels = 50, limits = cbind(c(0, 0.3), c(0.2, 0.4), c(0.4, 0.3))) %>% 
  set_names(c("R", "O", "UV"))


df_to_predict_2_2 = lattice_design_3_ing_2 %>% 
  mutate(
    Y = 0, G = 0, B = 0, P = 0
  ) %>% 
  select(R, O, Y, G, B, P, UV) %>% 
  expand_grid(
    .,
    intensity = seq(-1, 2, by = 0.05), 
    is_right = c(1, 0),
    no_choice = 0
  )



X_to_predict_list_2_2 <- create_model_matrix_second_order_scheffe(df_to_predict_2_2)



# Proper predictive
utilities_post_2_2 = matrix(0.0, ncol = n_draws, nrow = nrow(X_to_predict_list_2_2$X))

for(d in 1:ncol(utilities_post_2_2)){
  
  # I don't know if I have to sample from a Normal here, I think I don't because there's no variance.
  alpha_vec_d = alpha_draws[[1]][d, ]
  
  
  Sigma_level_0_d = Sigma_level_0_draws[[1]][d, , ]
  beta_level_0_d = beta_level_0_draws[[1]][d, ]
  beta_level_1 = mvrnorm(1, beta_level_0_d, Sigma_level_0_d)
  
  param_vec = c(beta_level_1, alpha_vec_d)
  
  utilities_d = X_to_predict_list_2_2$X %*% param_vec
  
  utilities_post_2_2[, d] = utilities_d
  
}


quantiles_utilities_2_2 = t(apply(utilities_post_2_2, 1, function(x) c(quantile(x, probs = c(0.1, 0.5, 0.9)), sd = sd(x)) ))
rm(utilities_post_2_2)


lattice_utilities_2_2 = df_to_predict_2_2 %>% 
  mutate(order_original = 1:n()) %>% 
  bind_cols(quantiles_utilities_2_2 %>% 
              as_tibble() %>% 
              set_names(c("p10", "p50", "p90", "sd"))) %>% 
  group_by(is_right) %>% 
  arrange(desc(p50)) %>% 
  ungroup() %>% 
  mutate(ix1 = 1:n())


lattice_utilities_2_2 %>% 
  group_by(is_right) %>% 
  top_n(3, p50) %>% 
  as.data.frame()

lattice_utilities_2_2 %>% 
  group_by(is_right) %>% 
  top_n(30, p50) %>% 
  as.data.frame() %>% 
  arrange(desc(p50))


lattice_utilities_2_2 %>% 
  filter(intensity > 0) %>% 
  group_by(is_right) %>% 
  top_n(3, p50) %>% 
  as.data.frame()






lattice_utilities_2_2 %>% 
  group_by(is_right) %>% 
  top_n(30, p50) %>% 
  as.data.frame() %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point() +
  geom_linerange(aes(ymin = p10, ymax = p90)) +
  theme_bw()






lattice_utilities_2_2 %>% 
  group_by(is_right) %>% 
  top_n(30, p50) %>% 
  bind_rows(
    lattice_utilities %>% 
      group_by(is_right) %>% 
      top_n(-30, p50)
  ) %>% 
  as.data.frame() %>% 
  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, y = p50)) +
  geom_point() +
  geom_linerange(aes(ymin = p10, ymax = p90)) +
  theme_bw()















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



lattice_utilities_bind_cols = lattice_utilities %>% 
  arrange(order_original) %>% 
  bind_cols(quantiles_utilities_plugin %>% 
              as_tibble() %>% 
              set_names(c("p10_pl", "p50_pl", "p90_pl", "sd_pl"))) %>% 
  group_by(is_right) %>% 
  arrange(desc(p50_pl)) %>% 
  ungroup() %>% 
  mutate(ix2 = 1:n())

lattice_utilities_bind_cols %>% 
  filter(ix1 != ix2)


lattice_utilities_bind_cols %>% 
  group_by(is_right) %>% 
  top_n(1, p50_pl) %>% 
  as.data.frame()


lattice_utilities_bind_cols %>% 
  group_by(is_right) %>% 
  top_n(3, p50) %>% 
  as.data.frame()




lattice_utilities_bind_cols %>% 
  filter(intensity > 0) %>% 
  group_by(is_right) %>% 
  top_n(3, p50) %>% 
  as.data.frame()

