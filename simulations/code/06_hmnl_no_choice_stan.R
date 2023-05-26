library(MASS)
library(tidyverse)
library(mlogit)
library(opdesmixr)
library(rstan)
# library(here)

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)

source("utils.R")



q = 5
J = 2
S = 25

n_flies_per_experiment = 80

# First four are mixture coefficients, last is no-choice
beta = c(3, 3, 3, 3, 4)


n_experiments = 8
sd_experiments = 2
Sigma = matrix(0.0, ncol = length(beta), nrow = length(beta))
diag(Sigma) = sd_experiments





# Create an I-optimal design
design_array_i_opt = opdesmixr::mnl_mixture_coord_exch(
  q = q, J = J, S = S, 
  beta = beta[1:(q-1)], seed = 2021, order = 1, opt_crit = "I", transform_beta = F,
  n_random_starts = 1, n_cores = 1, max_it = 2)



simulated_data = lapply(1:n_experiments, function(e){
  
  beta_exp = mvrnorm(1, beta, Sigma)
  
  # Simulate data for homogeneous flies in each experiment
  simulated_data_i_opt_01 = lapply(1:n_flies_per_experiment, function(r){
    
    print(paste("Experiment", e, ", Respondent", r))
    
    S = dim(design_array_i_opt$X)[3]
    
    sim_data_r = simulate_mixture_choice_data_no_choice(
      design_array_i_opt$X, 
      beta_vector_mixture = beta_exp[1:(q-1)],
      beta_no_choice = beta_exp[q],
      order = 1, 
      append_model_matrix = T)
    
    out = sim_data_r %>% 
      mutate(point = choice_set) %>% 
      mutate(
        choice_set = as.numeric(choice_set) + S*(r - 1),
        resp = r
      ) %>% 
      bind_cols(
        matrix(rep(beta_exp, (J+1)*S), ncol = length(beta_exp), byrow = T) %>% 
          as_tibble() %>% 
          set_names(paste0("real_beta_resp_", 1:ncol(.)))
      ) %>% 
      select(-utilities, -probs)
    
    return(out)
    
  }) %>% 
    bind_rows() %>% 
    mutate(experiment = e)
}) %>% 
  bind_rows()



# Get the real betas of each experiment
real_betas = simulated_data %>% 
  select(starts_with("real_beta")) %>% 
  distinct() %>% 
  as.matrix() %>% 
  t() %>% 
  as.vector()



# Grouped counts
simulated_data_counts = simulated_data %>% 
  group_by(experiment, point, across(starts_with("model_mat"))) %>% 
  summarize(Ny = sum(choice)) %>% 
  rename(choice_set = point) %>% 
  ungroup()





X_stan = simulated_data_counts %>% 
  select(starts_with("model_mat")) %>% 
  as.matrix()


unique_experiment_choice_set = simulated_data_counts %>% 
  select(experiment, choice_set) %>% 
  distinct()

index_dataframe = lapply(1:nrow(unique_experiment_choice_set), function(i){
  exp_i = unique_experiment_choice_set$experiment[i]
  cs_i =  unique_experiment_choice_set$choice_set[i]
  
  indices_i = which(simulated_data_counts$choice_set == cs_i & simulated_data_counts$experiment == exp_i)
  out = tibble(
    choice_set = cs_i, 
    start = indices_i[1], 
    end = indices_i[length(indices_i)],
    exp = exp_i
  )
}) %>% 
  bind_rows()


stan_data <- list(
  J = J+1, 
  S = nrow(index_dataframe), 
  n_exp = n_experiments,
  n_var = ncol(X_stan), 
  n_obs = nrow(X_stan),
  Ny = simulated_data_counts$Ny,
  X = X_stan,
  start = index_dataframe$start, # the starting observation for each decision/choice set
  end = index_dataframe$end,  # the ending observation for each decision/choice set
  experiment = index_dataframe$exp, # mapping for each experiment
  Sigma = Sigma # For first example with known variance
)




########################################################
########################################################
# Using a known Sigma variance matrix
########################################################
########################################################


# Around 40 seconds
model_text_known_sigma = paste(readLines("06_hmnl_no_choice_known_Sigma.stan"), collapse = "\n")
model_stan_known_sigma <- stan(
  model_code = model_text_known_sigma, 
  data = stan_data, 
  iter = 3000, 
  warmup = 2000,
  chains = 4, 
  seed = 2023, 
  cores = 4)

model_stan_known_sigma


summary(model_stan_known_sigma, pars = c("beta"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_known_sigma, pars = c("beta_exp"), probs = c(0.1, 0.5, 0.9))$summary


summary(model_stan_known_sigma, pars = c("beta_exp"), probs = c(0.1, 0.5, 0.9))$summary %>% 
  as.data.frame() %>% 
  select("mean", "10%", "90%") %>% 
  rownames_to_column() %>% 
  mutate(real_betas = real_betas) %>% 
  ggplot(aes(x = rowname)) +
  geom_point(aes(y = real_betas), color = "red", shape = "x", size = 4) +
  geom_point(aes(y = mean), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "dark grey") +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")



summary(model_stan_known_sigma, pars = c("beta"), probs = c(0.1, 0.5, 0.9))$summary %>% 
  as.data.frame() %>% 
  select("mean", "10%", "90%") %>% 
  rownames_to_column() %>% 
  mutate(real_betas = beta) %>% 
  ggplot(aes(x = rowname)) +
  geom_point(aes(y = real_betas), color = "red", shape = "x", size = 4) +
  geom_point(aes(y = mean), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "dark grey") +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")



plot(model_stan_known_sigma, plotfun="trace", pars=("beta"))







########################################################
########################################################
# Using an LKJ prior 
########################################################
########################################################



model_text_lkj = paste(readLines("06_hmnl_no_choice_lkj_prior_01.stan"), collapse = "\n")
# # Another way to specify the model thaht might be more efficient, but less clear (and different prior)
# model_text_lkj = paste(readLines("06_hmnl_no_choice_lkj_prior_02.stan"), collapse = "\n")
# Around 1.5 minutes
model_stan_lkj <- stan(
  model_code = model_text_lkj, 
  data = stan_data, 
  iter = 3000, 
  warmup = 2000,
  chains = 4, 
  seed = 2023, 
  cores = 4)

model_stan_lkj


summary(model_stan_lkj, pars = c("beta"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_lkj, pars = c("beta_exp"), probs = c(0.1, 0.5, 0.9))$summary


summary(model_stan_lkj, pars = c("beta_exp"), probs = c(0.1, 0.5, 0.9))$summary %>% 
  as.data.frame() %>% 
  select("mean", "10%", "90%") %>% 
  rownames_to_column() %>% 
  mutate(real_betas = real_betas) %>% 
  ggplot(aes(x = rowname)) +
  geom_point(aes(y = real_betas), color = "red", shape = "x", size = 4) +
  geom_point(aes(y = mean), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "dark grey") +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")



summary(model_stan_lkj, pars = c("beta"), probs = c(0.1, 0.5, 0.9))$summary %>% 
  as.data.frame() %>% 
  select("mean", "10%", "90%") %>% 
  rownames_to_column() %>% 
  mutate(real_betas = beta) %>% 
  ggplot(aes(x = rowname)) +
  geom_point(aes(y = real_betas), color = "red", shape = "x", size = 4) +
  geom_point(aes(y = mean), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "dark grey") +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")



plot(model_stan_lkj, plotfun="trace", pars=("beta"))
plot(model_stan_lkj, plotfun="trace", pars=("beta_exp"))


matrix(as.data.frame(summary(model_stan_lkj, pars = c("Sigma"), probs = c(0.1, 0.5, 0.9))$summary)$mean, 
       ncol = stan_data$n_var)





########################################################
########################################################
# Using an Inverse-Wishart prior for covariance matrix
########################################################
########################################################

model_text_iw = paste(readLines("06_hmnl_no_choice_iw_prior.stan"), collapse = "\n")

# model_stan_iw <- stan(model_code = model_text, data = stan_data, iter = 200, chains = 4, seed = 2023)
model_stan_iw <- stan(model_code = model_text_iw, data = stan_data, iter = 2000, chains = 4, seed = 2023, cores = 4)
model_stan_iw


plot(model_stan_iw, pars = "beta")
plot(model_stan_iw, pars = "beta_exp")


summary(model_stan_iw, pars = c("beta"), probs = c(0.1, 0.5, 0.9))$summary


summary(model_stan_iw, pars = c("beta_exp"), probs = c(0.1, 0.5, 0.9))$summary

summary(model_stan_iw, pars = c("beta_exp"), probs = c(0.1, 0.5, 0.9))$summary %>% 
  as.data.frame() %>% 
  select("mean", "10%", "90%") %>% 
  rownames_to_column() %>% 
  mutate(real_betas = real_betas) %>% 
  ggplot(aes(x = rowname)) +
  geom_point(aes(y = real_betas), color = "red", shape = "x", size = 4) +
  geom_point(aes(y = mean), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "dark grey") +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")



summary(model_stan_iw, pars = c("beta"), probs = c(0.1, 0.5, 0.9))$summary %>% 
  as.data.frame() %>% 
  select("mean", "10%", "90%") %>% 
  rownames_to_column() %>% 
  mutate(real_betas = beta) %>% 
  ggplot(aes(x = rowname)) +
  geom_point(aes(y = real_betas), color = "red", shape = "x", size = 4) +
  geom_point(aes(y = mean), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "dark grey") +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")




