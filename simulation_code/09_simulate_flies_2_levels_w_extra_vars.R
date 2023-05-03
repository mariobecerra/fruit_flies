library(MASS)
library(opdesmixr)
library(rstan)
library(forcats)
library(tidyverse)
library(here)

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)

source(here("simulation_code/utils.R"))

sims_out_folder = here("out/sims_out/")
dir.create(sims_out_folder, showWarnings = F)

date_today = Sys.Date()


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
S = 30
n_images_per_choice_set = 3

n_experiments = 6
n_flies_per_experiment = 80

# First are mixture coefficients, last is no-choice
beta_level_0 = c(2, 3, 4)


# Variation between experiments
Sigma_level_0 = diag(0.6, ncol = length(beta_level_0), nrow = length(beta_level_0))

# Variation between choice sets within same experiment
Sigma_level_1 = diag(0.15, ncol = length(beta_level_0), nrow = length(beta_level_0))

# Variation between images within same choice set
Sigma_level_2 = diag(0.01, ncol = length(beta_level_0), nrow = length(beta_level_0))



# Create an I-optimal design
design_array_i_opt = opdesmixr::mnl_mixture_coord_exch(
  q = q, J = J, S = S,  no_choice = T,
  beta = beta_level_0, seed = 2021, order = 1, opt_crit = "I", transform_beta = F,
  n_random_starts = 1, n_cores = 1, max_it = 10)



alpha_parameters = c(1, 1)

# Tibble with n_images_per_choice_set*S*n_experiments*(J+1) rows
simulated_data = lapply(1:n_experiments, function(e){
  cat(paste("\nExperiment", e, "of", n_experiments, "\n"))
  
  # Experiment level beta
  beta_level_1 = mvrnorm(1, beta_level_0, Sigma_level_0)
  
  z1 = rep(c(rep(1, floor(S/4)), rep(0, ceiling(S/4))), 2)
  z2 = c(rep(1, floor(S/2)), rep(0, ceiling(S/2)))
  
  out_exp_e = lapply(1:S, function(s){
    
    cat(paste("\tchoice set", s, "of", S, "\n"))
    
    beta_level_2 = mvrnorm(1, beta_level_1, Sigma_level_1)
    
    
    out_s = lapply(1:n_images_per_choice_set, function(im){
      
      # beta_level_3 = mvrnorm(1, beta_level_2, Sigma_level_2)
      
      sim_data_im = simulate_mixture_choice_data_no_choice_many_respondents_with_covariates(
        design_array = design_array_i_opt$X[, , s, drop = F], 
        # beta_vector_mixture = beta_level_3[1:(q-1)],
        # beta_no_choice = beta_level_3[q],
        beta_vector_mixture = beta_level_2[1:(q-1)],
        beta_no_choice = beta_level_2[q],
        extra_covariates_beta = alpha_parameters,
        extra_covariates_model_matrix = t(replicate(J, c(z1[s], z2[s]))),
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
        # bind_cols(
        #   matrix(rep(beta_level_3, (J+1)), ncol = length(beta_level_3), byrow = T) %>% 
        #     as_tibble() %>% 
        #     set_names(paste0("real_beta_level_3_", 1:ncol(.)))
        # ) %>% 
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


saveRDS(simulated_data, paste0(sims_out_folder, "model_2_levels_sim_simulated_data_tibble_", date_today, ".rds"))






model_mat_1 = simulated_data %>% 
  select(starts_with("model_mat"))
  
# Reorder columns so that first columns are mixture relates, then no_choice, and then the extra variables
n_mixture_cols = ncol(model_mat_1)-length(alpha_parameters)-1
order_columns_aux_1 = c(1:n_mixture_cols, ncol(model_mat_1))
order_columns_aux_2 = setdiff(1:ncol(model_mat_1), order_columns_aux_1)
n_extra_vars = length(order_columns_aux_2)
order_columns = c(order_columns_aux_1, order_columns_aux_2)

X_stan = model_mat_1 %>% 
  select(order_columns) %>% 
  as.matrix()


unique_experiment_choice_set = simulated_data %>%
  select(experiment, choice_set) %>%
  distinct()

index_dataframe = lapply(1:nrow(unique_experiment_choice_set), function(i){
  exp_i = unique_experiment_choice_set$experiment[i]
  cs_i =  unique_experiment_choice_set$choice_set[i]
  
  indices_i = which(simulated_data$choice_set == cs_i & simulated_data$experiment == exp_i)
  out = tibble(
    choice_set = cs_i,
    start = indices_i[1],
    end = indices_i[length(indices_i)],
    exp = exp_i
  )
}) %>%
  bind_rows()


n_betas_level_2 = nrow(index_dataframe)



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
  bind_rows() %>% 
  mutate(choice_set_2 = as.integer(fct(paste(exp, choice_set))))




# 
# n_betas_level_2 = unique_experiment_choice_set_image %>% 
#   select(experiment, choice_set) %>% 
#   distinct() %>% 
#   nrow()



stan_data <- list(
  n_alt = J+1, 
  n_exp = n_experiments,
  n_betas_level_2 = n_betas_level_2,
  n_var = ncol(X_stan), 
  n_obs = nrow(X_stan),
  Ny = simulated_data$Ny,
  X = X_stan,
  experiment = index_dataframe$exp,
  choice_set = index_dataframe$choice_set,
  
  start = index_dataframe_image$start, # the starting observation for each image
  end = index_dataframe_image$end,  # the ending observation for each image
  n_images = nrow(index_dataframe_image),
  # experiment_index_2 = index_dataframe_image$exp,
  cs_index_2 = index_dataframe_image$choice_set_2,
  
  n_mixture_cols = n_mixture_cols,
  n_extra_vars = n_extra_vars
)





# 3 minutes for 1500 iterations and 6 experiments, 50 choice sets and 5 images per choice set
# 15 hours for 4000 iterations, 8 experiments, 120 choice sets per experiment, and 25 images per choice set
model_stan_lkj <- stan(
  file = here("simulation_code/09_hmnl_flies_2_levels_w_extra_vars.stan"),
  data = stan_data, 
  seed = 2023,
  # iter = 50, chains = 1, cores = 1)
  # iter = 4000,  warmup = 3000, chains = 6, cores = 6)
  iter = 1500,  warmup = 1000, chains = 4, cores = 4)

model_stan_lkj


saveRDS(model_stan_lkj, paste0(sims_out_folder, "model_2_levels_sim_stan_object_", date_today, ".rds"))


summary(model_stan_lkj, pars = c("alpha"), probs = c(0.1, 0.5, 0.9))$summary


summary(model_stan_lkj, pars = c("beta_level_0"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_lkj, pars = c("beta_level_1"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_lkj, pars = c("beta_level_2"), probs = c(0.1, 0.5, 0.9))$summary



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
  geom_point(aes(y = mean), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "dark grey") +
  geom_point(aes(y = real_betas), color = "red", shape = "x", size = 4) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")



summary(model_stan_lkj, pars = c("beta_level_2"), probs = c(0.1, 0.5, 0.9))$summary %>% 
  as.data.frame() %>% 
  select("mean", "10%", "90%") %>% 
  rownames_to_column() %>% 
  mutate(real_betas = real_betas_level_2) %>% 
  arrange(desc(mean)) %>% 
  mutate(ix = 1:n()) %>% 
  ggplot(aes(x = ix)) +
  geom_point(aes(y = mean), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "dark grey") +
  geom_point(aes(y = real_betas), color = "red", shape = "x", size = 4) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")





plot(model_stan_lkj, plotfun="trace", pars=("beta_level_0"))
# plot(model_stan_lkj, plotfun="trace", pars=("beta_level_1"))
# plot(model_stan_lkj, plotfun="trace", pars=("beta_level_2"))


matrix(as.data.frame(summary(model_stan_lkj, pars = c("Sigma_level_0"), probs = c(0.5))$summary)$`50%`, 
       ncol = stan_data$n_mixture_cols+1)

matrix(as.data.frame(summary(model_stan_lkj, pars = c("Sigma_level_1"), probs = c(0.5))$summary)$`50%`, 
       ncol = stan_data$n_mixture_cols+1)






mean_beta_level_1 = summary(model_stan_lkj, pars = c("beta_level_1"), probs = c(0.1, 0.5, 0.9))$summary %>% 
  as.data.frame() %>% 
  select("mean", "10%", "90%") %>% 
  pull("mean")


summary(model_stan_lkj, pars = c("beta_level_2"), probs = c(0.1, 0.5, 0.9))$summary %>% 
  as.data.frame() %>% 
  select("50%", "10%", "90%", "mean") %>% 
  rownames_to_column() %>% 
  mutate(real_betas = real_betas_level_2) %>% 
  ggplot +
  geom_point(aes(x = real_betas_level_2, y = mean)) + 
  geom_hline(yintercept = mean_beta_level_1) +
  geom_abline(slope = 1) +
  theme_bw()






