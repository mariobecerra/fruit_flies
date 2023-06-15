# Fit 3 models and compare:
# 1) simple MNL model with 30 images per choice set
# 2) simple MNL model with the maximum number of flies per choice set
# 3) Hierarchical MNL model with 30 images per choice set

library(MASS)
library(rstan)
library(opdesmixr)
library(tidyverse)
library(here)



init_fun <- function() {
  out_list = list(
    # tau_level_0 = rep(1, n_mixture_cols+1), 
    # Omega_level_0 = diag(1, ncol = n_mixture_cols+1, nrow = n_mixture_cols+1)
    # tau_level_0 = rnorm(n_mixture_cols+1, mean = 1, sd = 0.1), 
    tau_level_0 = rgamma(n_mixture_cols+1, shape = 4, scale = 0.25),
    Omega_level_0 = rlkjcorr(K = n_mixture_cols+1, eta = 30),
    L_Omega_level_0 = chol(rlkjcorr(K = n_mixture_cols+1, eta = 30))
  )
  
  return(out_list)
}

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)

source(here("japanese_fly/code/japanese_flies_hmnl_1_level_utils.R"))


counts_japanese = read_csv(here("japanese_fly/out/counts_japanese_choice_format.csv"))


counts_japanese_filtered = counts_japanese %>% 
  group_by(experiment, choice_set) %>% 
  mutate(
    min_time = min(date_time), 
    max_time = max(date_time)
  ) %>%  
  mutate(elapsed_mins = difftime(max_time ,min_time, "mins")) %>% 
  filter(date_time > min_time + 0.25*elapsed_mins) %>%  # remove the first 25% of the data (first 5 minutes)
  ungroup() %>% 
  as.data.frame()


n_images_per_cs = counts_japanese_filtered %>%
  group_by(experiment, folder, choice_set) %>%
  summarize(n_images = length(unique(image)))


n_images_per_cs_to_sample = 30
# Take subset of images
set.seed(2023)
images_indices = counts_japanese_filtered %>%
  select(experiment, folder, choice_set, image) %>% 
  distinct() %>% 
  group_by(experiment, folder, choice_set) %>% 
  slice_sample(n = n_images_per_cs_to_sample) %>% 
  arrange(experiment, choice_set, image) %>% 
  ungroup()


counts_japanese_sample <- images_indices %>% 
  left_join(counts_japanese_filtered) %>% 
  as.data.frame() %>% 
  mutate(time_day = as.numeric(substr(as.character(date_time), 12, 13))) %>% 
  mutate(
    is_right = ifelse(alternative == "2_right", 1, 0),
    is_daytime = ifelse(time_day >= 9 & time_day <= 20, 1, 0) # 1 = 9-20h, 0 = 21-23h + 0-8h
  ) 





unique_experiment_choice_set = counts_japanese_sample %>%
  select(experiment, folder) %>%
  distinct()


index_dataframe = lapply(1:nrow(unique_experiment_choice_set), function(i){
  exp_i = unique_experiment_choice_set$experiment[i]
  cs_i =  unique_experiment_choice_set$folder[i]
  
  indices_i = which(counts_japanese_sample$folder == cs_i & counts_japanese_sample$experiment == exp_i)
  out = tibble(
    choice_set = cs_i,
    start = indices_i[1],
    end = indices_i[length(indices_i)],
    exp = exp_i
  )
}) %>%
  bind_rows()


unique_experiment_choice_set_image = counts_japanese_sample %>%
  select(experiment, folder, image) %>%
  distinct()

index_dataframe_image = lapply(1:nrow(unique_experiment_choice_set_image), function(i){
  exp_i = unique_experiment_choice_set_image$experiment[i]
  cs_i =  unique_experiment_choice_set_image$folder[i]
  image_i =  unique_experiment_choice_set_image$image[i]
  
  indices_i = which(counts_japanese_sample$folder == cs_i & counts_japanese_sample$experiment == exp_i & counts_japanese_sample$image == image_i)
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



df_to_fit <- counts_japanese_sample %>% 
  select(all_of(c('R','O','Y','G','B','P','UV','intensity', 'no_choice'))) %>%
  as.data.frame()


X_stan_list_01 = create_model_matrix_second_order_scheffe(df_to_fit)



n_experiments = max(counts_japanese_sample$experiment)
n_mixture_cols = dim(X_stan_list_01$X)[2] - 1











# Simple MNL model
model_text_01 = "
data {
    int<lower=2> J; // number of alternatives/outcomes
    int<lower=2> S; // number of choice sets
    int<lower=1> n_var; // of covariates
    int<lower=1> n_obs; // of observations
    vector[n_obs] Ny; // counts of the response variable
    matrix[n_obs, n_var] X; // attribute matrix
    int nrow_image_index_df;
    int start[nrow_image_index_df]; // the starting observation for each image
    int end[nrow_image_index_df]; // the ending observation for each image
    int n_images;
}

parameters {
    vector[n_var] beta;  // attribute effects 
}


model {
    vector[J] utilities;
    
    beta ~ normal(0, 5);
    
    for(i in 1:n_images){
      utilities = X[start[i]:end[i]]*beta;
      target += sum(log_softmax(utilities) .* Ny[start[i]:end[i]]);
      
    }
}

// // Save log-likelihood for LOO package
// generated quantities{
//   vector[n_images] log_lik;
//   
//   for(i in 1:n_images){
//       log_lik[i] = sum(log_softmax(X[start[i]:end[i]]*beta) .* Ny[start[i]:end[i]]);
//     }
// }

"



stan_data_01 <- list(
  J = 3, 
  S = max(counts_japanese_sample$choice_set), 
  n_var = ncol(X_stan_list_01$X), 
  n_obs = nrow(X_stan_list_01$X),
  Ny = counts_japanese_sample$count,
  X = X_stan_list_01$X,
  nrow_image_index_df = nrow(index_dataframe_image),
  start = index_dataframe_image$start, # the starting observation for each decision/choice set
  end = index_dataframe_image$end,  # the ending observation for each decision/choice set
  n_images = nrow(index_dataframe_image)
)

model_01_filename = here(paste0("japanese_fly/out/japanese_mnl_model_experiments_1_to_", n_experiments, "_stan_object_", n_images_per_cs_to_sample, "images.rds"))

if(file.exists(model_01_filename)){
  cat("File already exists")
  model_stan_01 = readRDS(model_01_filename)
} else{
  
  (t1 = Sys.time())
  # 26 minutes for 500 iterations and 4 chains in 4 cores.
  model_stan_01 <- stan(
    model_code = model_text_01, 
    data = stan_data_01, 
    iter = 500,  warmup = 350, chains = 4, cores = 4,
    # iter = 1500, warmup = 1000,
    # iter = 500, warmup = 200,
    seed = 2023)
  (t2 = Sys.time())
  t2 - t1
  
  saveRDS(model_stan_01, model_01_filename)
}


model_stan_01


betas_model_01_summary = summary(model_stan_01, pars = c("beta"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  mutate(variable = colnames(X_stan_list_01$X),
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F))

betas_model_01_summary %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  theme_bw() +
  xlab("Beta") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



betas_model_01_summary %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), color = "dark grey") +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), color = "black", size = 1) +
  theme_bw() +
  xlab("Beta") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

















max_counts = counts_japanese %>% 
  group_by(experiment, choice_set, folder, alternative) %>% 
  summarize(count = max(count)) %>% 
  ungroup() %>% 
  filter(alternative != "3_none")




count_data_flies_max = max_counts %>% 
  full_join(
    counts_japanese %>% 
      select(experiment, alternative, choice_set, folder, R, O, Y, G, B, P, UV, intensity, no_choice) %>% 
      distinct()
  ) %>% 
  arrange(experiment, choice_set, folder, alternative) %>% 
  group_by(experiment, choice_set, folder) %>% 
  mutate(no_choice_count = 80 - sum(count, na.rm = T)) %>% 
  mutate(count = ifelse(is.na(count), no_choice_count, count)) %>% 
  select(-no_choice_count) %>% 
  ungroup()


X_stan_list_02 = create_model_matrix_second_order_scheffe(
  count_data_flies_max %>% 
    select(all_of(c('R','O','Y','G','B','P','UV','intensity','no_choice')))
)



unique_experiment_choice_set_max = count_data_flies_max %>%
  select(experiment, choice_set) %>%
  distinct()

index_dataframe_cs = lapply(1:nrow(unique_experiment_choice_set_max), function(i){
  exp_i = unique_experiment_choice_set_max$experiment[i]
  cs_i =  unique_experiment_choice_set_max$choice_set[i]
  
  indices_i = which(count_data_flies_max$choice_set == cs_i & count_data_flies_max$experiment == exp_i)
  out = tibble(
    choice_set = cs_i,
    start = indices_i[1],
    end = indices_i[length(indices_i)],
    exp = exp_i
  )
}) %>%
  bind_rows() 


stan_data_02 <- list(
  J = 3,
  S = nrow(index_dataframe_cs),
  n_var = ncol(X_stan_list_02$X),
  n_obs = nrow(X_stan_list_02$X),
  Ny = count_data_flies_max$count,
  X = X_stan_list_02$X,
  start = index_dataframe_cs$start, # the starting observation for each decision/choice set
  end = index_dataframe_cs$end  # the ending observation for each decision/choice set
)






model_text_02 = "
data {
    int<lower=2> J; // number of alternatives/outcomes
    int<lower=2> S; // number of choice sets
    int<lower=1> n_var; // of covariates
    int<lower=1> n_obs; // of observations
    vector[n_obs] Ny; // counts of the response variable
    matrix[n_obs, n_var] X; // attribute matrix
    int start[S]; // the starting observation for each decision/choice set
    int end[S]; // the ending observation for each decision/choice set
}

parameters {
    vector[n_var] beta;  // attribute effects 
}

model {
    vector[J] utilities;
    
    for(s in 1:S){
      utilities = X[start[s]:end[s]]*beta;
      target += sum(log_softmax(utilities) .* Ny[start[s]:end[s]]);
    }
}


// Save log-likelihood for LOO package
generated quantities{
 vector[S] log_lik;
 
 vector[J] utilities;
 
 for(s in 1:S){
    utilities = X[start[s]:end[s]]*beta;
    log_lik[s] = sum(log_softmax(utilities) .* Ny[start[s]:end[s]]);
   }
}

"


(t1 = Sys.time())
# 2 minutes for 1000 iterations and 4 chains in 4 cores
model_stan_02 <- stan(
  model_code = model_text_02, 
  data = stan_data_02, 
  warmup = 850,
  iter = 1000, 
  chains = 4, 
  seed = 2023,
  cores = 4)
(t2 = Sys.time())
t2 - t1

model_stan_02




betas_model_02_summary = summary(model_stan_02, pars = c("beta"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  mutate(variable = colnames(X_stan_list_02$X),
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F))

betas_model_02_summary %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  theme_bw() +
  xlab("Beta") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



betas_model_02_summary %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), color = "dark grey") +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), color = "black", size = 1) +
  theme_bw() +
  xlab("Beta") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))













betas_model_01_summary %>% 
  mutate(model = "model_1") %>% 
  bind_rows(
    betas_model_02_summary %>% 
      mutate(model = "model_2")
  ) %>% 
  ggplot(aes(x = variable, color = model)) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 1, position = position_dodge(width = 0.5)) +
  theme_bw() +
  xlab("Beta") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))










stan_data_03 <- list(
  n_alt = 3, 
  n_exp = n_experiments,
  n_var = ncol(X_stan_list_01$X), 
  n_obs = nrow(X_stan_list_01$X),
  Ny = counts_japanese_sample$count,
  X = X_stan_list_01$X,
  start = index_dataframe_image$start, # the starting observation for each image
  end = index_dataframe_image$end,  # the ending observation for each image
  n_images = nrow(index_dataframe_image),
  exp_index = index_dataframe_image$exp,
  
  # prior
  tau_level_0_mean = 1,
  tau_level_0_sd = 0.5,
  L_Omega_level_0_param = 5,
  beta_level_0_mean = rep(0, ncol(X_stan_list_01$X)),
  beta_level_0_sd = rep(2, ncol(X_stan_list_01$X))
)


# With 4 experiments and 50 images
# Gradient evaluation took 0.054098 seconds
# 000 transitions using 10 leapfrog steps per transition would take 540.98 seconds

# With 4 experiments and 20 images:
# Gradient evaluation took 0.035647 seconds
# 1000 transitions using 10 leapfrog steps per transition would take 356.47 seconds.
# 17 minutes for 50 iterations
#
# 3 hours with 200 iterations
# There were 400 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
# The largest R-hat is 1.21, indicating chains have not mixed.
# Running the chains for more iterations may help
#
# 4.5 hours with 500 iterations (350 warmup)
# There were 179 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10

model_03_filename = here(paste0("japanese_fly/out/japanese_hmnl_model_experiments_1_to_", n_experiments, "_stan_object_", n_images_per_cs_to_sample, "images.rds"))

if(file.exists(model_03_filename)){
  cat("File already exists")
  model_stan_03 = readRDS(model_03_filename)
} else{
  Sys.time()
  model_stan_03 <- stan(
    file = here("japanese_fly/code/japanese_flies_hmnl_1_level_reparameterized.stan"),
    data = stan_data,
    seed = 2023,
    include = F, pars = c("z", "L_Omega_level_0", "tau_level_0"),
    # iter = 1500,  warmup = 1000, chains = 4, cores = 4,
    # iter = 2500,  warmup = 2000, chains = 4, cores = 4,
    # iter = 3000,  warmup = 2000, chains = 6, cores = 6,
    # iter = 500,  chains = 4, cores = 4,
    # iter = 5, chains = 1,
    # iter = 50, chains = 4, cores = 4,
    iter = 500,  warmup = 350, chains = 4, cores = 4,
    init = init_fun
  )
  Sys.time()
  saveRDS(model_stan_03, model_03_filename)
}

model_stan_03


betas_level_0_summary_03 = summary(model_stan_03, pars = c("beta_level_0"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_stan_list_01$X),
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F))



betas_level_0_summary_03 %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))







betas_model_01_summary %>% 
  select(variable, mean, sd) %>% 
  mutate(model = "model_1") %>% 
  bind_rows(
    betas_model_02_summary %>% 
      select(variable, mean, sd) %>% 
      mutate(model = "model_2")
  ) %>% 
  bind_rows(
    betas_level_0_summary_03 %>% 
      select(variable, mean, sd) %>% 
      mutate(model = "model_3")
  ) %>% 
  ggplot(aes(x = variable, color = model)) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 1, position = position_dodge(width = 0.5)) +
  theme_bw() +
  xlab("Beta") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))




beta_level_0_draws_hmnl = rstan::extract(model_stan_03, "beta_level_0")
beta_level_1_draws_hmnl = rstan::extract(model_stan_03, "beta_level_1")
Sigma_level_0_draws_hmnl = rstan::extract(model_stan_03, "Sigma_level_0")
n_draws = nrow(beta_level_0_draws_hmnl[[1]])

#########################################
####### Predictions in data used to fit
#########################################


predict_n_flies_original_data = function(stan_data, data_used_to_fit, beta_level_1_draws_hmnl, n_flies_opt = 80, J = 3){
  n_draws = dim(beta_level_1_draws_hmnl[[1]])[1]
  n_exp = dim(beta_level_1_draws_hmnl[[1]])[2]
    

  out_all_experiments = lapply(1:n_exp, function(ex){
    
    print(paste("Experiment", ex))
    
    exp_indices = which(data_used_to_fit$experiment == ex)
    n_obs_ex = length(exp_indices)
    
    utilities_post = matrix(0.0, ncol = n_draws, nrow = n_obs_ex)
    probs_post = array(0.0, c(n_draws, J, n_obs_ex/J))
    N_flies_post = array(0.0, dim(probs_post))
    
    for(d in 1:n_draws){
      
      beta_level_1_d_ex = beta_level_1_draws_hmnl[[1]][d, ex, ]
      
      
      utilities_d = stan_data$X[exp_indices, ] %*% beta_level_1_d_ex
    
      utilities_post[, d] = utilities_d
      
      
      choice_sets = rep(1:(n_obs_ex/J), each = J)
      
      probs_d = sapply(split(utilities_d, choice_sets), function(x){
        exp_x = exp(x)
        exp_x/sum(exp_x)
      })
      
      for(i in 1:ncol(probs_d)){
        N_flies_post[d, , i] = as.numeric(rmultinom(n = 1, size = n_flies_opt, prob = probs_d[, i]))
      }
      
      
    }
    
    
    # Predictions for number of flies
    mean_N_flies_post = apply(N_flies_post, c(2, 3), mean)
    median_N_flies_post = apply(N_flies_post, c(2, 3), median)
    p10_N_flies_post = apply(N_flies_post, c(2, 3), quantile, 0.1)
    p90_N_flies_post = apply(N_flies_post, c(2, 3), quantile, 0.9)
    
    design_df_preds = data.frame(
      experiment = ex,
      choice_set = choice_sets,
      Ny = stan_data$Ny[exp_indices]
    ) %>% 
      mutate(
        mean_N_flies = as.numeric(mean_N_flies_post),
        median_N_flies = as.numeric(median_N_flies_post),
        p10_N_flies = as.numeric(p10_N_flies_post),
        p90_N_flies = as.numeric(p90_N_flies_post)
      )
    
    return(
      list(
        utilities_post = utilities_post,
        probs_post = probs_post,
        N_flies_post = N_flies_post,
        design_df_preds = design_df_preds
      )
    )
    
  })
  
}




back_predictions = predict_n_flies_original_data(stan_data_03, counts_japanese_sample, beta_level_1_draws_hmnl, n_flies_opt = 80, J = 3)

glimpse(back_predictions)

back_preds_df = lapply(1:length(back_predictions), function(i){
  back_predictions[[i]]$design_df_preds
}) %>% 
  bind_rows() %>% 
  mutate(no_choice = counts_japanese_sample$no_choice)

back_preds_df %>% 
  ggplot(aes(x = Ny)) +
  geom_jitter(aes(y = median_N_flies), size = 0.3) +
  facet_wrap(~experiment) +
  theme_bw()


back_preds_df %>% 
  filter(no_choice == 0) %>% 
  ggplot(aes(x = Ny)) +
  geom_jitter(aes(y = median_N_flies), size = 0.3) +
  facet_wrap(~experiment) +
  theme_bw()


back_preds_df %>% 
  filter(no_choice == 1) %>% 
  ggplot(aes(x = Ny)) +
  geom_jitter(aes(y = median_N_flies), size = 0.3) +
  facet_wrap(~experiment) +
  theme_bw()






back_preds_df %>% 
  ggplot(aes(x = Ny)) +
  geom_jitter(aes(y = median_N_flies, color = as.character(experiment)), size = 0.3, alpha = 0.6) +
  theme_bw()


back_preds_df %>% 
  filter(no_choice == 0) %>% 
  ggplot(aes(x = Ny)) +
  geom_jitter(aes(y = median_N_flies, color = as.character(experiment)), size = 0.7, alpha = 0.6) +
  theme_bw()


back_preds_df %>% 
  filter(no_choice == 1) %>% 
  ggplot(aes(x = Ny)) +
  geom_jitter(aes(y = median_N_flies, color = as.character(experiment)), size = 0.7, alpha = 0.6) +
  theme_bw()




#########################################
####### Predictions in new data
#########################################



predict_n_flies_design = function(
  design_df, beta_level_0_draws_hmnl, Sigma_level_0_draws_hmnl, J = 3, n_flies_opt = 80
  ){
  
  correct_names_nc = c('R','O','Y','G','B','P', 'UV', 'intensity', 'no_choice', 'choice_set')
  correct_names_no_nc = c('R','O','Y','G','B','P', 'UV', 'intensity', 'choice_set')
  
  all_eq_1 = all.equal(names(design_df), correct_names_no_nc)
  all_eq_2 = all.equal(names(design_df), correct_names_nc)
  
  if(is.character(all_eq_1) & is.character(all_eq_2)){
    stop("Variables in design_df must be either c(", paste(correct_names_no_nc, collapse = ", "), ")",
         " or c(", paste(correct_names_nc, collapse = ", "), ")")
  }
  
  if(is.character(all_eq_1)) {
    no_choice_included = T
  } else{
    no_choice_included = F
  }
  
  n_draws = nrow(beta_level_0_draws_hmnl[[1]])
  
  
  # If the provided data frame doesn't include the "no_choice", include it
  if(!no_choice_included){
    design_df = design_df %>% 
      mutate(
        no_choice = 0,
        alt = rep(1:(J-1), nrow(.)/(J-1))
      ) %>% 
      bind_rows(
        as.data.frame(matrix(0.0, ncol = ncol(design_df), nrow = nrow(design_df)/2)) %>% 
          set_names(correct_names_no_nc) %>% 
          mutate(no_choice = 1, alt = J, choice_set = 1:nrow(.))
      ) %>% 
      arrange(choice_set, alt) %>% 
      select(-alt)
  }
  
  
  
  X_to_predict_list <- create_model_matrix_second_order_scheffe(
    design_df %>% 
      select(-choice_set)
  )
  
  
  utilities_post = matrix(0.0, ncol = n_draws, nrow = nrow(X_to_predict_list$X))
  probs_post = array(0.0, c(n_draws, J, nrow(X_to_predict_list$X)/J))
  N_flies_post = array(0.0, dim(probs_post))
  
  for(d in 1:ncol(utilities_post)){
    
    Sigma_level_0_d = Sigma_level_0_draws_hmnl[[1]][d, , ]
    beta_level_0_d = beta_level_0_draws_hmnl[[1]][d, ]
    beta_level_1 = mvrnorm(1, beta_level_0_d, Sigma_level_0_d)
    
    utilities_d = X_to_predict_list$X %*% beta_level_1
    
    utilities_post[, d] = utilities_d
    
    probs_d = sapply(split(utilities_d, design_df$choice_set), function(x){
      exp_x = exp(x)
      exp_x/sum(exp_x)
    })
    
    for(i in 1:ncol(probs_d)){
      N_flies_post[d, , i] = as.numeric(rmultinom(n = 1, size = n_flies_opt, prob = probs_d[, i]))
    }
    
    probs_post[d, , ] = probs_d
    
  }
  
  
  
  # Predictions for number of flies
  mean_N_flies_post = apply(N_flies_post, c(2, 3), mean)
  median_N_flies_post = apply(N_flies_post, c(2, 3), median)
  p10_N_flies_post = apply(N_flies_post, c(2, 3), quantile, 0.1)
  p90_N_flies_post = apply(N_flies_post, c(2, 3), quantile, 0.9)
  
  
  design_df_preds = design_df %>% 
    mutate(
      mean_N_flies = as.numeric(mean_N_flies_post),
      median_N_flies = as.numeric(median_N_flies_post),
      p10_N_flies = as.numeric(p10_N_flies_post),
      p90_N_flies = as.numeric(p90_N_flies_post)
    )
  
  return(
    list(
      utilities_post = utilities_post,
      probs_post = probs_post,
      N_flies_post = N_flies_post,
      design_df_preds = design_df_preds
    )
  )
  
}


# If design 2 was a new experiment, these would be the predictions
design_df_exp2 = counts_japanese_sample %>% 
  filter(experiment == 2) %>% 
  select(R, O, Y, G, B, P, UV, intensity, no_choice, choice_set, alternative) %>%
  rename(cs_orig = choice_set) %>% 
  mutate(choice_set = as.integer(as.factor(cs_orig))) %>% 
  distinct()


preds_design_exp2 = predict_n_flies_design(
  design_df = design_df_exp2 %>% select(-cs_orig, -alternative), 
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl, 
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl, 
  J = 3, 
  n_flies_opt = 80
)


counts_and_preds_exp2 = design_df_exp2 %>% 
  bind_cols(
    preds_design_exp2$design_df_preds %>% 
      select(mean_N_flies, median_N_flies, p10_N_flies, p90_N_flies)
  ) %>% 
  select(-choice_set, -R, -O, -Y, -G, -B, -P, -UV, -intensity, -no_choice) %>% 
  rename(choice_set = cs_orig) %>% 
  full_join(
    counts_japanese_sample %>% 
      filter(experiment == 2)
  ) %>% 
  arrange(choice_set, image, alternative)

sum(counts_and_preds_exp2$count >= counts_and_preds_exp2$p10_N_flies & counts_and_preds_exp2$count <= counts_and_preds_exp2$p90_N_flies)

mean(counts_and_preds_exp2$count >= counts_and_preds_exp2$p10_N_flies & counts_and_preds_exp2$count <= counts_and_preds_exp2$p90_N_flies)




counts_and_preds_exp2 %>% 
  ggplot(aes(x = count)) +
  geom_jitter(aes(y = median_N_flies))


counts_and_preds_exp2 %>% 
  filter(no_choice == 0) %>% 
  ggplot(aes(x = count)) +
  geom_jitter(aes(y = median_N_flies))


counts_and_preds_exp2 %>% 
  filter(no_choice == 1) %>% 
  ggplot(aes(x = count)) +
  geom_jitter(aes(y = median_N_flies))








#### 4th I-optimal design

i_opt_designs_list = readRDS(here("japanese_fly/out/japanese_flies_4th_i_optimal_design.rds"))

i_opts = sapply(seq_along(i_opt_designs_list), function(i) i_opt_designs_list[[i]]$opt_crit_value)

i_opt_design = i_opt_designs_list[[which.min(i_opts)]]

i_opt_design_df = mnl_design_array_to_dataframe(i_opt_design$X) %>% 
  set_names(c("R", "O", "Y", "G", "B", "P", "UV", "intensity", "cs")) %>% 
  mutate(choice_set = as.numeric(cs)) %>% 
  as.data.frame()



preds_i_opt_design = predict_n_flies_design(
  design_df = i_opt_design_df %>% select(-cs), 
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl, 
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl, 
  J = 3, 
  n_flies_opt = 80
)



i_opt_design_df_preds = preds_i_opt_design$design_df_preds

head(i_opt_design_df_preds)




#### Lattice (need to check)

lattice_design_7_ing = create_lattice_design(n_var = 7, n_levels = 5) %>% 
  set_names(c('R','O','Y','G','B','P', 'UV'))




df_to_predict_01 = expand_grid(
  lattice_design_7_ing, 
  intensity = seq(-1, 2, by = 0.2),
  no_choice = 0
) %>% 
  as.data.frame()




predict_utilities_design = function(
  df_to_predict, beta_level_0_draws_hmnl, Sigma_level_0_draws_hmnl,
  include_all_utilities = T
  ){
  X_to_predict_list <- create_model_matrix_second_order_scheffe(df_to_predict)
  
  # Proper predictive
  utilities_post = matrix(0.0, ncol = n_draws, nrow = nrow(X_to_predict_list$X))
  
  for(d in 1:ncol(utilities_post)){
    Sigma_level_0_d = Sigma_level_0_draws_hmnl[[1]][d, , ]
    beta_level_0_d = beta_level_0_draws_hmnl[[1]][d, ]
    beta_level_1 = mvrnorm(1, beta_level_0_d, Sigma_level_0_d)
    
    utilities_d = X_to_predict_list$X %*% beta_level_1
    
    utilities_post[, d] = utilities_d
  }
  
  
  quantiles_utilities = t(apply(utilities_post, 1, function(x) c(quantile(x, probs = c(0.1, 0.5, 0.9)), sd = sd(x)) ))
  
  
  out_utilities_summary = df_to_predict %>% 
    mutate(order_original = 1:n()) %>% 
    bind_cols(quantiles_utilities %>% 
                as_tibble() %>% 
                set_names(c("p10", "p50", "p90", "sd")))
  
  
  if(include_all_utilities){
    out = list(
      utilities_summary = out_utilities_summary,
      utilities_post = utilities_post
    )
  }else{
    out = list(
      utilities_summary = out_utilities_summary
    )
  }
  
}



predicted_utilities_df_01 = predict_utilities_design(
  df_to_predict = df_to_predict_01, 
  beta_level_0_draws_hmnl = beta_level_0_draws_hmnl, 
  Sigma_level_0_draws_hmnl = Sigma_level_0_draws_hmnl,
  include_all_utilities = T
)



predicted_utilities_df_01$utilities_summary %>% 
  arrange(desc(p50))











library("loo")


log_lik_2 <- extract_log_lik(model_stan_02, merge_chains = FALSE)
r_eff_2 <- relative_eff(exp(log_lik_2), cores = 4) 
loo_2 <- loo(log_lik_2, r_eff = r_eff_2, cores = 4)
print(loo_2)





