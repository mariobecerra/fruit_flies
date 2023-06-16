library(rstan)
library(tidyverse)
library(here)

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)


# Check that experiment number is correct!!!!!!!!!!!
counts = read_csv(here("japanese_fly/out/counts_japanese_choice_format.csv")) %>% 
  filter(experiment <= 8)




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






count_data_flies_all = counts %>% 
  group_by(choice_set) %>% 
  mutate(
    min_time = min(date_time), 
    max_time = max(date_time)
  ) %>%  
  mutate(elapsed_mins = difftime(max_time ,min_time, "mins")) %>% 
  filter(date_time > min_time + 0.25*elapsed_mins) %>%  # remove the first 25% of the data (first 5 minutes)
  ungroup()



# # Take subset of images
# set.seed(2023)
# images_indices = count_data_flies_all %>%
#   select(experiment, folder, choice_set, image) %>% 
#   distinct() %>% 
#   group_by(experiment, folder, choice_set) %>% 
#   slice_sample(n = 20) %>% 
#   arrange(experiment, choice_set, image) %>% 
#   ungroup()
# 
# 
# count_data_flies <- images_indices %>% 
#   left_join(count_data_flies_all) %>% 
#   as.data.frame() %>% 
#   mutate(time_day = as.numeric(substr(as.character(date_time), 12, 13))) %>% 
#   mutate(
#     is_right = ifelse(alternative == "2_right", 1, 0),
#     is_daytime = ifelse(time_day >= 9 & time_day <= 20, 1, 0) # 1 = 9-20h, 0 = 21-23h + 0-8h
#   ) 





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


X_stan_list_01 = create_model_matrix_second_order_scheffe(
  count_data_flies_max %>% 
    select(all_of(c('R','O','Y','G','B','P','UV','intensity','no_choice')))
)


X_stan_01 <- X_stan_list_01$X


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


stan_data_01 <- list(
  J = 3,
  S = nrow(index_dataframe_cs),
  n_var = ncol(X_stan_01),
  n_obs = nrow(X_stan_01),
  Ny = count_data_flies_max$count,
  X = X_stan_01,
  start = index_dataframe_cs$start, # the starting observation for each decision/choice set
  end = index_dataframe_cs$end  # the ending observation for each decision/choice set
)






model_text_01 = "
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

# 110 seconds for 5000 iterations and 4 chains in 4 cores
model_stan_01 <- stan(
  model_code = model_text_01, 
  data = stan_data_01, 
  warmup = 2000,
  iter = 5000, 
  chains = 4, 
  seed = 2023,
  cores = 4)

model_stan_01

# saveRDS(model_stan_01, here("japanese_fly/out/japanese_max_mnl_model_experiments_1_to_8_stan_object.rds"))


betas_model_01_summary = summary(model_stan_01, pars = c("beta"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  mutate(variable = colnames(X_stan_list_01$X),
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F))

saveRDS(betas_model_01_summary, here("japanese_fly/out/japanese_max_mnl_model_experiments_1_to_8_betas_summary.rds"))


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





