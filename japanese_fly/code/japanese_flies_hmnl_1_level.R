library(rstan)
library(tidyverse)
library(here)

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)



create_model_matrix_second_order_scheffe = function(df_to_convert){
  
  correct_var_names = c("R", "O", "Y", "G", "B", "P", "UV", "intensity", "no_choice")
  # correct_var_names = c("R", "O", "Y", "G", "B", "P", "UV", "intensity", "is_right", "is_daytime", "no_choice")
  
  if(!all.equal(names(df_to_convert), correct_var_names)) {
    stop("Names of df_to_convert must be c(", paste(correct_var_names, collapse = ", "), ")")
  }
  
  
  
  mixture_variable_names_no_UV = c('R','O','Y','G','B','P')
  mixture_variable_names = c(mixture_variable_names_no_UV, 'UV')
  
  mixture_pairwise_interaction_names = c(
    'R*O','R*Y','R*G','R*B','R*P','R*UV',
    'O*Y','O*G','O*B','O*P','O*UV',
    'Y*G','Y*B','Y*P','Y*UV',
    'G*B','G*P','G*UV',
    'B*P','B*UV',
    'P*UV')
  
  mixture_intensity_interaction_names = c('R*intensity','O*intensity','Y*intensity','G*intensity',
                                          'B*intensity','P*intensity','UV*intensity')
  
  
  X_other_vars_nc <- df_to_convert %>% 
    # dplyr::select(all_of(c('no_choice','is_right', 'is_daytime'))) %>%
    dplyr::select(all_of(c('no_choice'))) %>%
    as.data.frame()
  
  X_main <- df_to_convert %>% 
    dplyr::select(all_of(c(mixture_variable_names,'intensity'))) %>%
    as.data.frame()
  
  
  # X_color_pairwise <- combn(mixture_variable_names, 2, function(x) X_main[x[1]]*X_main[x[2]], simplify = FALSE) %>%
  #   bind_cols() %>%
  #   set_names(mixture_pairwise_interaction_names)
  
  X_color_pairwise = combn(seq_along(mixture_variable_names), 2, function(x) X_main[, x[1]]*X_main[, x[2]], simplify = T) %>% 
    as.data.frame() %>% 
    set_names(mixture_pairwise_interaction_names)
  
  
  X_color_intensity <- X_main %>%
    dplyr::select(all_of(c(mixture_variable_names,'intensity'))) %>%
    mutate(across(c(mixture_variable_names,'intensity'), ~.*intensity)) %>%
    set_names(c(mixture_intensity_interaction_names,'intensity^2'))
  
  
  X <- X_main %>% 
    dplyr::select(-UV, -intensity) %>% 
    bind_cols(X_color_pairwise) %>% 
    bind_cols(X_color_intensity) %>%
    bind_cols(X_other_vars_nc) %>% 
    as.matrix()
  
  names_betas_level_1 = c(mixture_variable_names_no_UV, mixture_pairwise_interaction_names, mixture_intensity_interaction_names,'intensity^2', "no_choice")
  
  return(list(
    X = X,
    mixture_variable_names_no_UV = mixture_variable_names_no_UV,
    mixture_variable_names = mixture_variable_names,
    mixture_pairwise_interaction_names = mixture_pairwise_interaction_names,
    mixture_intensity_interaction_names = mixture_intensity_interaction_names,
    names_betas_level_1 = names_betas_level_1
  ))
}


rlkjcorr = function(K, eta = 1)  {
  stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
  stopifnot(eta > 0)
  
  alpha <- eta + (K - 2)/2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, K, K)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (K > 2) 
    for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  return(crossprod(R))
}







counts = read_csv(here("out/counts_japanese_choice_format.csv")) %>% 
  filter(experiment <= 2)


# counts %>% group_by(choice_set) %>% summarize(min_time = min(date_time), max_time = max(date_time)) %>% mutate(elapsed_mins = as.numeric(difftime(max_time ,min_time, "mins"))) 



count_data_flies_all = counts %>% 
  group_by(choice_set) %>% 
  mutate(
    min_time = min(date_time), 
    max_time = max(date_time)
  ) %>%  
  mutate(elapsed_mins = difftime(max_time ,min_time, "mins")) %>% 
  filter(date_time > min_time + 0.25*elapsed_mins) %>%  # remove the first 25% of the data (first 5 minutes)
  ungroup()



n_images_per_cs = count_data_flies_all %>% 
  group_by(experiment, folder, choice_set) %>% 
  summarize(n_images = max(image))


# Take subset of images
set.seed(2023)
images_indices = count_data_flies_all %>%
  left_join(n_images_per_cs) %>% 
  filter(n_images > 1) %>% 
  select(experiment, folder, choice_set, image) %>% 
  distinct() %>% 
  group_by(experiment, folder, choice_set) %>% 
  slice_sample(n = 20) %>% 
  arrange(experiment, choice_set, image) %>% 
  ungroup()


count_data_flies <- images_indices %>% 
  left_join(count_data_flies_all) %>% 
  as.data.frame() %>% 
  mutate(time_day = as.numeric(substr(as.character(date_time), 12, 13))) %>% 
  mutate(
    is_right = ifelse(alternative == "2_right", 1, 0),
    is_daytime = ifelse(time_day >= 9 & time_day <= 20, 1, 0) # 1 = 9-20h, 0 = 21-23h + 0-8h
  ) 





X_stan_list_01 = create_model_matrix_second_order_scheffe(
  count_data_flies %>% 
    select(all_of(c('R','O','Y','G','B','P','UV','intensity','no_choice')))
)


X_stan_01 <- X_stan_list_01$X


unique_experiment_choice_set_image = count_data_flies %>%
  select(experiment, folder, image) %>%
  distinct()

index_dataframe_image = lapply(1:nrow(unique_experiment_choice_set_image), function(i){
  exp_i = unique_experiment_choice_set_image$experiment[i]
  cs_i =  unique_experiment_choice_set_image$folder[i]
  image_i =  unique_experiment_choice_set_image$image[i]
  
  indices_i = which(count_data_flies$folder == cs_i & count_data_flies$experiment == exp_i & count_data_flies$image == image_i)
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

// Save log-likelihood for LOO package
generated quantities{
  vector[n_images] log_lik;
  
  for(i in 1:n_images){
      log_lik[i] = sum(log_softmax(X[start[i]:end[i]]*beta) .* Ny[start[i]:end[i]]);
    }
}



"



stan_data_01 <- list(
  J = 3, 
  S = max(count_data_flies$choice_set), 
  n_var = ncol(X_stan_01), 
  n_obs = nrow(X_stan_01),
  Ny = count_data_flies$count,
  X = X_stan_01,
  nrow_image_index_df = nrow(index_dataframe_image),
  start = index_dataframe_image$start, # the starting observation for each decision/choice set
  end = index_dataframe_image$end,  # the ending observation for each decision/choice set
  n_images = nrow(index_dataframe_image)
)




# 5 minutes for 500 iterations and 4 chains in 4 cores.
# 9 minutes for 1500 iterations and 4 chains in 4 cores.
model_stan_01 <- stan(
  model_code = model_text_01, 
  data = stan_data_01, 
  iter = 1500, warmup = 1000,
  # iter = 500, warmup = 200,
  chains = 4, 
  seed = 2023,
  cores = 4)

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



# To use as prior in design
model_01_summary = summary(model_stan_01, pars = c("beta"), probs = c(0.5))$summary %>% 
  as.data.frame() %>% 
  mutate(variable = colnames(X_stan_01)) %>% 
  select(variable, mean, sd) 


saveRDS(model_stan_01, here("out/japanese_stan_object_mnl_model_2nd_experiment_20_images_per_cs.rds"))
saveRDS(model_01_summary, here("out/japanese_model_mnl_model_2nd_experiment_20_images_per_cs.rds"))















max_counts = count_data_flies_all %>% 
  group_by(experiment, choice_set, folder, alternative) %>% 
  summarize(count = max(count)) %>% 
  ungroup() %>% 
  filter(alternative != "3_none")




count_data_flies_max = max_counts %>% 
  full_join(
    counts %>% 
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


X_stan_02 <- X_stan_list_02$X


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
  n_var = ncol(X_stan_02),
  n_obs = nrow(X_stan_02),
  Ny = count_data_flies_max$count,
  X = X_stan_02,
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
    
    for(s in 1:S){
    
      vector[J] utilities;
      
      utilities = X[start[s]:end[s]]*beta;
      
      target += sum(log_softmax(utilities) .* Ny[start[s]:end[s]]);
    }
  
}
"

# 60 seconds for 5000 iterations and 4 chains in 4 cores
model_stan_02 <- stan(
  model_code = model_text_02, 
  data = stan_data_02, 
  warmup = 2000,
  iter = 5000, 
  chains = 4, 
  seed = 2023,
  cores = 4)

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




# To use as prior in design
model_02_summary = summary(model_stan_02, pars = c("beta"), probs = c(0.5))$summary %>% 
  as.data.frame() %>% 
  mutate(variable = colnames(X_stan_02)) %>% 
  select(variable, mean, sd) 


saveRDS(model_stan_02, here("out/japanese_stan_object_mnl_model_2nd_experiment_max_counts.rds"))
saveRDS(model_02_summary, here("out/japanese_model_mnl_model_2nd_experiment__max_counts.rds"))











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
  n_exp = max(count_data_flies$experiment),
  n_var = ncol(X_stan_01), 
  n_obs = nrow(X_stan_01),
  Ny = count_data_flies$count,
  X = X_stan_01,
  
  start = index_dataframe_image$start, # the starting observation for each image
  end = index_dataframe_image$end,  # the ending observation for each image
  n_images = nrow(index_dataframe_image),
  exp_index = index_dataframe_image$exp
)



# Stan model ------------

init_fun <- function() {
  out_list = list(
    tau_level_0 = rgamma(ncol(X_stan_01), shape = 4, scale = 0.25),
    Omega_level_0 = rlkjcorr(K = ncol(X_stan_01), eta = 30)
  )
  
  return(out_list)
}


Sys.time()
# 5 minutes for 100 iterations and 4 chains in 4 cores.
# 1 hour for 1500 iterations and 4 chains in 4 cores, 217 divergent transitions.
# 2 hours with 3500 iter (2500 warmup), 1004 divergent transitions after warmup. 
# 3 hours with 5000 iter (3000 warmup), 1237 divergent transitions after warmup.
model_stan_03 <- stan(
  file = here("code/japanese_flies_05_hmnl_for_next_experiment.stan"),
  data = stan_data_03,
  seed = 2023,
  # iter = 1500,  warmup = 1000, chains = 4, cores = 4,
  iter = 5000,  warmup = 3000, chains = 4, cores = 4,
  # iter = 50, chains = 1,
  init = init_fun,
  save_warmup = F
)

Sys.time()

saveRDS(model_stan_03, here("out/japanese_stan_object_hmnl_model_2nd_experiment.rds"))

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







betas_model_01_summary_03 %>% 
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






betas_level_1_summary_03 = summary(model_stan_03, pars = c("beta_level_1"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(
    exp = paste0(
      "Exp ",
      rep(1:max(count_data_flies$experiment), each = length(colnames(X_stan_list_01$X)))
    ),
    variable = rep(colnames(X_stan_list_01$X), max(count_data_flies$experiment))
  ) %>% 
  group_by(exp) %>% 
  mutate(ix = 1:n()) %>%
  ungroup() %>% 
  # mutate(variable = paste0(exp, ", ", var)) %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = T)) 



betas_level_1_summary_03 %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = F)) %>% 
  ggplot(aes(x = variable, color = exp)) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.6), size = 0.9) +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), size = 0.6, position = position_dodge(width = 0.6)) +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), size = 0.9, position = position_dodge(width = 0.6)) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  geom_vline(xintercept = 1:36 + 0.5, size = 0.3, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())












library("rstanarm")
library("bayesplot")
library("loo")


log_lik_1 <- extract_log_lik(model_stan_01, merge_chains = FALSE)
r_eff_1 <- relative_eff(exp(log_lik_1), cores = 4) 
loo_1 <- loo(log_lik_1, r_eff = r_eff_1, cores = 4)
print(loo_1)




log_lik_3 <- extract_log_lik(model_stan_03, merge_chains = FALSE)
r_eff_3 <- relative_eff(exp(log_lik_3), cores = 4) 
loo_3 <- loo(log_lik_3, r_eff = r_eff_3, cores = 4)
print(loo_3)


comp <- loo_compare(loo_1, loo_3)

