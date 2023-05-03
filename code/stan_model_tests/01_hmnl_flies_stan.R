library(mlogit)
library(rstan)
library(tidyverse)
library(here)

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)


count_data_01 = read_csv(here("data/eur_max.csv")) %>% 
  mutate(experiment = 1) %>% 
  rename(choice_set = choiceset,
         no_choice = D)


unique_experiment_choice_set_01 = count_data_01 %>% 
  select(experiment, choice_set) %>% 
  distinct()

index_dataframe_01 = lapply(1:nrow(unique_experiment_choice_set_01), function(i){
  exp_i = unique_experiment_choice_set_01$experiment[i]
  cs_i =  unique_experiment_choice_set_01$choice_set[i]
  
  indices_i = which(count_data_01$choice_set == cs_i & count_data_01$experiment == exp_i)
  out = tibble(
    choice_set = cs_i, 
    start = indices_i[1], 
    end = indices_i[length(indices_i)],
    exp = exp_i
  )
}) %>% 
  bind_rows()

X_stan_01 = count_data_01 %>% 
  select(R, O, Y, G, B, P, intensity) %>% 
  as.matrix()


stan_data_01 <- list(
  J = nrow(count_data_01)/max(count_data_01$choice_set), 
  S = max(count_data_01$choice_set), 
  n_exp = 1,
  n_var = ncol(X_stan_01), 
  n_obs = nrow(X_stan_01),
  Ny = count_data_01$max,
  X = X_stan_01,
  start = index_dataframe_01$start, # the starting observation for each decision/choice set
  end = index_dataframe_01$end,  # the ending observation for each decision/choice set
  experiment = index_dataframe_01$exp # mapping for each experiment
)

model_stan_01 <- stan(
  file = here("code/01_hmnl_flies_01.stan"),
  data = stan_data_01, 
  iter = 10000, 
  warmup = 8000,
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
  mutate(variable = colnames(X_stan_01)) %>% 
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









count_data_02 = read_csv(here("data/eur_max.csv")) %>% 
  mutate(experiment = 1) %>% 
  rename(choice_set = choiceset,
         no_choice = D) %>% 
  filter(no_choice == 0)

unique_experiment_choice_set_02 = count_data_02 %>% 
  select(experiment, choice_set) %>% 
  distinct()

index_dataframe_02 = lapply(1:nrow(unique_experiment_choice_set_02), function(i){
  exp_i = unique_experiment_choice_set_02$experiment[i]
  cs_i =  unique_experiment_choice_set_02$choice_set[i]
  
  indices_i = which(count_data_02$choice_set == cs_i & count_data_02$experiment == exp_i)
  out = tibble(
    choice_set = cs_i, 
    start = indices_i[1], 
    end = indices_i[length(indices_i)],
    exp = exp_i
  )
}) %>% 
  bind_rows()

X_stan_02 = count_data_02 %>% 
  select(R, O, Y, G, B, P, intensity) %>% 
  as.matrix()


stan_data_02 <- list(
  J = nrow(count_data_02)/max(count_data_02$choice_set), 
  S = max(count_data_02$choice_set), 
  n_exp = 1,
  n_var = ncol(X_stan_02), 
  n_obs = nrow(X_stan_02),
  Ny = count_data_02$max,
  X = X_stan_02,
  start = index_dataframe_02$start, # the starting observation for each decision/choice set
  end = index_dataframe_02$end,  # the ending observation for each decision/choice set
  experiment = index_dataframe_02$exp # mapping for each experiment
)

model_stan_02 <- stan(
  file = here("code/01_hmnl_flies_01.stan"),
  data = stan_data_02, 
  iter = 10000, 
  warmup = 8000,
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
  mutate(variable = colnames(X_stan_02)) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")




plot(model_stan_02, plotfun="trace", pars=("beta"))
plot(model_stan_02, plotfun="trace", pars=("beta_exp"))
