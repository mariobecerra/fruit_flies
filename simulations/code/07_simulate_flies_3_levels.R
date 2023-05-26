library(MASS)
library(tidyverse)
library(opdesmixr)
library(rstan)
library(here)

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)

source(here("simulation_code/utils.R"))



# q = 7
# J = 2
# S = 100
# n_images_per_choice_set = 250
# 
# n_experiments = 8
# n_flies_per_experiment = 80
# 
# # First are mixture coefficients, last is no-choice
# beta_level_0 = c(3, 3, 3, 3, 3, 3, 4)




q = 3
J = 2
S = 10
n_images_per_choice_set = 10

n_experiments = 2
n_flies_per_experiment = 80

# First are mixture coefficients, last is no-choice
beta_level_0 = c(3, 3, 4)


# Variation between experiments
Sigma_level_0 = diag(0.6, ncol = length(beta_level_0), nrow = length(beta_level_0))

# Variation between choice sets within same experiment
Sigma_level_1 = diag(0.11, ncol = length(beta_level_0), nrow = length(beta_level_0))

# Variation between images within same choice set
Sigma_level_2 = diag(0.15, ncol = length(beta_level_0), nrow = length(beta_level_0))



# Create an I-optimal design
design_array_i_opt = opdesmixr::mnl_mixture_coord_exch(
  q = q, J = J, S = S,  no_choice = T,
  beta = beta_level_0, seed = 2021, order = 1, opt_crit = "I", transform_beta = F,
  n_random_starts = 1, n_cores = 1, max_it = 4)



# Tibble with n_images_per_choice_set*S*n_experiments*(J+1) rows
simulated_data = lapply(1:n_experiments, function(e){
  cat(paste("\nExperiment", e, "of", n_experiments, "\n"))
  
  # Experiment level beta
  beta_level_1 = mvrnorm(1, beta_level_0, Sigma_level_0)
  
  
  out_exp_e = lapply(1:S, function(s){
    
    cat(paste("\tchoice set", s, "of", S, "\n"))
    
    beta_level_2 = mvrnorm(1, beta_level_1, Sigma_level_1)
    
    
    out_s = lapply(1:n_images_per_choice_set, function(im){
      
      beta_level_3 = mvrnorm(1, beta_level_2, Sigma_level_2)
      
      sim_data_im = simulate_mixture_choice_data_no_choice_many_respondents(
        design_array = design_array_i_opt$X[, , s, drop = F], 
        beta_vector_mixture = beta_level_3[1:(q-1)],
        beta_no_choice = beta_level_3[q],
        n_resp = n_flies_per_experiment,
        order = 1, 
        append_model_matrix = T
      )
      
      
      out = sim_data_im %>% 
        mutate(point = choice_set) %>% 
        mutate(
          choice_set = s
        ) %>% 
        select(-utilities, -probs) %>% 
        bind_cols(
          matrix(rep(beta_level_3, (J+1)), ncol = length(beta_level_3), byrow = T) %>% 
            as_tibble() %>% 
            set_names(paste0("real_beta_level_3_", 1:ncol(.)))
        ) %>% 
        bind_cols(
          matrix(rep(beta_level_2, (J+1)), ncol = length(beta_level_2), byrow = T) %>% 
            as_tibble() %>% 
            set_names(paste0("real_beta_level_2_", 1:ncol(.)))
        ) %>% 
        bind_cols(
          matrix(rep(beta_level_1, (J+1)), ncol = length(beta_level_1), byrow = T) %>% 
            as_tibble() %>% 
            set_names(paste0("real_beta_level_1_", 1:ncol(.)))
        ) %>% 
        mutate(image = im)
      
      return(out)
      
    }) %>% 
      bind_rows() %>% 
      mutate(s = s)
    
    return(out_s)
    
    
  }) %>% 
    bind_rows() %>% 
    mutate(experiment = e)
  
  return(out_exp_e)
}) %>% 
  bind_rows()








X_stan = simulated_data %>% 
  select(starts_with("model_mat")) %>% 
  as.matrix()


# unique_experiment_choice_set = simulated_data %>% 
#   select(experiment, choice_set) %>% 
#   distinct()
# 
# index_dataframe = lapply(1:nrow(unique_experiment_choice_set), function(i){
#   exp_i = unique_experiment_choice_set$experiment[i]
#   cs_i =  unique_experiment_choice_set$choice_set[i]
#   
#   indices_i = which(simulated_data$choice_set == cs_i & simulated_data$experiment == exp_i)
#   out = tibble(
#     choice_set = cs_i, 
#     start = indices_i[1], 
#     end = indices_i[length(indices_i)],
#     exp = exp_i
#   )
# }) %>% 
#   bind_rows()






unique_experiment_choice_set_image = simulated_data %>% 
  select(experiment, choice_set, image) %>% 
  distinct()


index_dataframe_image = lapply(1:nrow(unique_experiment_choice_set_image), function(i){
  exp_i = unique_experiment_choice_set_image$experiment[i]
  cs_i =  unique_experiment_choice_set_image$choice_set[i]
  image_i =  unique_experiment_choice_set_image$image[i]
  
  indices_i = which(simulated_data$choice_set == cs_i & simulated_data$experiment == exp_i & simulated_data$image == image_i)
  out = tibble(
    exp = exp_i,
    choice_set = cs_i, 
    image = image_i,
    start = indices_i[1], 
    end = indices_i[length(indices_i)]
  )
}) %>% 
  bind_rows()


n_betas_level_2 = unique_experiment_choice_set_image %>% 
  select(experiment, choice_set) %>% 
  distinct() %>% 
  nrow()



stan_data <- list(
  n_alt = J+1, 
  n_exp = n_experiments,
  n_betas_level_2 = n_betas_level_2,
  n_betas_level_3 = nrow(index_dataframe_image), # same as nrow(X_stan)/(J+1),
  n_var = ncol(X_stan), 
  n_obs = nrow(X_stan),
  Ny = simulated_data$Ny,
  X = X_stan,
  start = index_dataframe_image$start, # the starting observation for each image
  end = index_dataframe_image$end,  # the ending observation for each image
  image = index_dataframe_image$image,
  experiment = index_dataframe_image$exp,
  choice_set = index_dataframe_image$choice_set
)







# 300 seconds (5 minutes) for 500 iterations
# With these little iterations the model doesn't recover well
# Maybe it's lack of iterations, maybe it's lack of data.

# I also think that the specification of the model in Stan is wrong.
# I think I need to change the choice set indices so that there's a unique index for every choice set across all experiments, like in the model with two levels in files 12_blah_blah

model_stan_lkj <- stan(
  file = here("simulation_code/07_hmnl_flies_3_levels.stan"),
  data = stan_data, 
  seed = 2023,
  # iter = 50, chains = 1, cores = 1)
  iter = 500,  warmup = 100, chains = 2, cores = 2)
  # iter = 3000,
  # warmup = 2000,
  # chains = 4,
  # seed = 2023,
  # cores = 4)

model_stan_lkj


summary(model_stan_lkj, pars = c("beta_level_0"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_lkj, pars = c("beta_level_1"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_lkj, pars = c("beta_level_2"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_lkj, pars = c("beta_level_3"), probs = c(0.1, 0.5, 0.9))$summary




summary(model_stan_lkj, pars = c("beta_level_0"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_lkj, pars = c("beta_level_1"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_lkj, pars = c("beta_level_2"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_lkj, pars = c("beta_level_3"), probs = c(0.1, 0.5, 0.9))$summary


summary(model_stan_lkj, pars = c("beta_level_0"), probs = c(0.1, 0.5, 0.9))$summary %>% 
  as.data.frame() %>% 
  select("mean", "10%", "90%") %>% 
  rownames_to_column() %>% 
  mutate(real_betas = beta_level_0) %>% 
  ggplot(aes(x = rowname)) +
  geom_point(aes(y = real_betas), color = "red", shape = "x", size = 4) +
  geom_point(aes(y = mean), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "dark grey") +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")





# Get the real betas of each experiment
real_betas_level_1 = simulated_data %>% 
  select(starts_with("real_beta_level_1")) %>% 
  distinct() %>% 
  as.matrix() %>% 
  t() %>% 
  as.vector()

summary(model_stan_lkj, pars = c("beta_level_1"), probs = c(0.1, 0.5, 0.9))$summary %>% 
  as.data.frame() %>% 
  select("mean", "10%", "90%") %>% 
  rownames_to_column() %>% 
  mutate(real_betas = real_betas_level_1) %>% 
  ggplot(aes(x = rowname)) +
  geom_point(aes(y = real_betas_level_1), color = "red", shape = "x", size = 4) +
  geom_point(aes(y = mean), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "dark grey") +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")



# Get the real betas of each choice set
real_betas_level_2 = simulated_data %>% 
  select(starts_with("real_beta_level_2")) %>% 
  distinct() %>% 
  as.matrix() %>% 
  t() %>% 
  as.vector()



summary(model_stan_lkj, pars = c("beta_level_2"), probs = c(0.1, 0.5, 0.9))$summary %>% 
  as.data.frame() %>% 
  select("mean", "10%", "90%") %>% 
  rownames_to_column() %>% 
  mutate(real_betas = real_betas_level_2) %>% 
  ggplot(aes(x = rowname)) +
  geom_point(aes(y = real_betas_level_2), color = "red", shape = "x", size = 4) +
  geom_point(aes(y = mean), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "dark grey") +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")




# Get the real betas of each choice set
real_betas_level_3 = simulated_data %>% 
  select(starts_with("real_beta_level_3")) %>% 
  distinct() %>% 
  as.matrix() %>% 
  t() %>% 
  as.vector()



summary(model_stan_lkj, pars = c("beta_level_3"), probs = c(0.1, 0.5, 0.9))$summary %>% 
  as.data.frame() %>% 
  select("mean", "10%", "90%") %>% 
  rownames_to_column() %>% 
  mutate(real_betas = real_betas_level_3) %>% 
  ggplot(aes(x = rowname)) +
  geom_point(aes(y = real_betas_level_3), color = "red", shape = "x", size = 4) +
  geom_point(aes(y = mean), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "dark grey") +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")




plot(model_stan_lkj, plotfun="trace", pars=("beta_level_0"))
plot(model_stan_lkj, plotfun="trace", pars=("beta_level_1"))
plot(model_stan_lkj, plotfun="trace", pars=("beta_level_2"))


matrix(as.data.frame(summary(model_stan_lkj, pars = c("Sigma_level_0"), probs = c(0.1, 0.5, 0.9))$summary)$mean, 
       ncol = stan_data$n_var)






