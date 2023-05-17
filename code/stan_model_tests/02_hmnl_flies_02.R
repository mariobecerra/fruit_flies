library(mlogit)
library(rstan)
library(tidyverse)
library(here)

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)


design_mapping = read_csv("out/design_mapping_european.csv")
counts_european = read_csv(here("out/counts_european_choice_format.csv"))

max_counts = counts_european %>% 
  group_by(experiment, choice_set, folder, alternative) %>% 
  summarize(count = max(count)) %>% 
  ungroup %>% 
  filter(alternative != "3_none")




count_data_flies = max_counts %>% 
  full_join(
    counts_european %>% 
      select(experiment, alternative, choice_set, folder, R, O, Y, G, B, P, UV, intensity, no_choice) %>% 
      distinct()
  ) %>% 
  arrange(experiment, choice_set, folder, alternative) %>% 
  group_by(experiment, choice_set, folder) %>% 
  mutate(no_choice_count = 80 - sum(count, na.rm = T)) %>% 
  mutate(count = ifelse(is.na(count), no_choice_count, count)) %>% 
  select(-no_choice_count) %>% 
  ungroup()






unique_experiment_choice_set_01 = count_data_flies %>% 
  select(experiment, choice_set) %>% 
  distinct()

index_dataframe_01 = lapply(1:nrow(unique_experiment_choice_set_01), function(i){
  exp_i = unique_experiment_choice_set_01$experiment[i]
  cs_i =  unique_experiment_choice_set_01$choice_set[i]
  
  indices_i = which(count_data_flies$choice_set == cs_i & count_data_flies$experiment == exp_i)
  out = tibble(
    choice_set = cs_i, 
    start = indices_i[1], 
    end = indices_i[length(indices_i)],
    exp = exp_i
  )
}) %>% 
  bind_rows()

X_stan_01 = count_data_flies %>% 
  select(R, O, Y, G, B, P, intensity, no_choice) %>% 
  as.matrix()


stan_data_01 <- list(
  J = nrow(count_data_flies)/max(count_data_flies$choice_set), 
  S = max(count_data_flies$choice_set), 
  n_exp = max(count_data_flies$experiment),
  n_var = ncol(X_stan_01), 
  n_obs = nrow(X_stan_01),
  Ny = count_data_flies$count,
  X = X_stan_01,
  start = index_dataframe_01$start, # the starting observation for each decision/choice set
  end = index_dataframe_01$end,  # the ending observation for each decision/choice set
  experiment = index_dataframe_01$exp # mapping for each experiment
)

model_stan_01 <- stan(
  file = here("code/02_hmnl_flies_max_counts.stan"),
  data = stan_data_01, 
  iter = 10000,
  warmup = 8000,
  # iter = 300, 
  # warmup = 100,
  chains = 4, 
  seed = 2023, 
  cores = 4)

model_stan_01


summary(model_stan_01, pars = c("beta"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_01, pars = c("beta_exp"), probs = c(0.1, 0.5, 0.9))$summary


summary(model_stan_01, pars = c("beta"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_stan_01)) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")





summary(model_stan_01, pars = c("beta_exp"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = paste0("exp_", rep(1:max(count_data_flies$experiment), each = ncol(X_stan_01)), "_", colnames(X_stan_01))) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")


# Same plot as before but in different order
summary(model_stan_01, pars = c("beta_exp"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = paste0(colnames(X_stan_01), "_exp_", rep(1:max(count_data_flies$experiment), each = ncol(X_stan_01)))) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")


plot(model_stan_01, plotfun="trace", pars=("beta"))
plot(model_stan_01, plotfun="trace", pars=("beta_exp"))













X_main <- count_data_flies %>%
  select(c('R','O','Y','G','B','P','UV','intensity','no_choice')) %>%
  as.data.frame()

X_color_pairwise <- combn(c('R','O','Y','G','B','P','UV'), 2, function(x) X_main[x[1]]*X_main[x[2]], simplify = FALSE) %>%
  bind_cols() %>%
  set_names(c('R*O','R*Y','R*G','R*B','R*P','R*UV',
              'O*Y','O*G','O*B','O*P','O*UV',
              'Y*G','Y*B','Y*P','Y*UV',
              'G*B','G*P','G*UV',
              'B*P','B*UV',
              'P*UV'))

X_color_intensity <- X_main %>%
  select(c('R','O','Y','G','B','P','UV','intensity')) %>%
  mutate(across(c('R','O','Y','G','B','P','UV','intensity'), ~.*intensity)) %>%
  set_names(c('R*intensity','O*intensity','Y*intensity','G*intensity',
              'B*intensity','P*intensity','UV*intensity','intensity^2'))








X_stan_02 <- X_main %>% 
  select(-UV, -intensity) %>% 
  bind_cols(X_color_pairwise) %>% 
  # bind_cols(X_color_intensity) %>% 
  as.matrix()





stan_data_02 <- list(
  J = nrow(count_data_flies)/max(count_data_flies$choice_set), 
  S = max(count_data_flies$choice_set), 
  n_exp = max(count_data_flies$experiment),
  n_var = ncol(X_stan_02), 
  n_obs = nrow(X_stan_02),
  Ny = count_data_flies$count,
  X = X_stan_02,
  start = index_dataframe_01$start, # the starting observation for each decision/choice set
  end = index_dataframe_01$end,  # the ending observation for each decision/choice set
  experiment = index_dataframe_01$exp # mapping for each experiment
)

# 5 minutes for 500 iterations
model_stan_02 <- stan(
  file = here("code/nc_model_iw_prior.stan"),
  data = stan_data_02, 
  iter = 1000,
  warmup = 500,
  # iter = 300, 
  # warmup = 100,
  chains = 4, 
  seed = 2023, 
  cores = 4)

model_stan_02


summary(model_stan_02, pars = c("beta"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_02, pars = c("beta_exp"), probs = c(0.1, 0.5, 0.9))$summary


summary(model_stan_02, pars = c("beta"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_stan_02)) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")





summary(model_stan_02, pars = c("beta_exp"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = paste0("exp_", rep(1:max(count_data_flies$experiment), each = ncol(X_stan_02)), "_", colnames(X_stan_02))) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")


# Same plot as before but in different order
summary(model_stan_02, pars = c("beta_exp"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = paste0(colnames(X_stan_02), "_exp_", rep(1:max(count_data_flies$experiment), each = ncol(X_stan_02)))) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")


# plot(model_stan_02, plotfun="trace", pars=("beta"))
# plot(model_stan_02, plotfun="trace", pars=("beta_exp"))


















X_stan_03 <- X_main %>% 
  select(-UV, -intensity) %>% 
  bind_cols(X_color_pairwise) %>% 
  bind_cols(X_color_intensity) %>% 
  as.matrix()





stan_data_03 <- list(
  J = nrow(count_data_flies)/max(count_data_flies$choice_set), 
  S = max(count_data_flies$choice_set), 
  n_exp = max(count_data_flies$experiment),
  n_var = ncol(X_stan_03), 
  n_obs = nrow(X_stan_03),
  Ny = count_data_flies$count,
  X = X_stan_03,
  start = index_dataframe_01$start, # the starting observation for each decision/choice set
  end = index_dataframe_01$end,  # the ending observation for each decision/choice set
  experiment = index_dataframe_01$exp # mapping for each experiment
)

model_stan_03 <- stan(
  file = here("code/nc_model_iw_prior.stan"),
  data = stan_data_03, 
  iter = 5000,
  warmup = 3000,
  # iter = 300, 
  # warmup = 100,
  chains = 4, 
  seed = 2023, 
  cores = 4)

model_stan_03


