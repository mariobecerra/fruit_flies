library(tidyverse)
library(rstan)
library(here)

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)



counts_european = read_csv(here("out/counts_european_choice_format.csv")) %>% 
  filter(experiment != 1) %>% 
  mutate(experiment = experiment - 1)

# Take subset of images
set.seed(2023)
images_indices = counts_european %>%
  select(experiment, folder, choice_set, image) %>% 
  distinct() %>% 
  group_by(experiment, folder, choice_set) %>% 
  slice_sample(n = 10) %>% 
  arrange(experiment, choice_set, image) %>% 
  ungroup()

# There are some choice sets with only one image. Check that.
data_dmela_sample <- images_indices %>% 
  left_join(counts_european) %>% 
  as.data.frame() %>% 
  mutate(time_day = substr(as.character(date_time), 12, 13)) %>% 
  mutate(
    is_right = ifelse(alternative == "2_right", 1, 0),
    is_daytime = ifelse(time_day >= 9 & time_day <= 20, 1, 0) # 1 = 9-20h, 0 = 21-23h + 0-8h
  ) 
  



# Stan data ------------
unique_experiment_choice_set = data_dmela_sample %>%
  select(experiment, folder) %>%
  distinct()


index_dataframe = lapply(1:nrow(unique_experiment_choice_set), function(i){
  exp_i = unique_experiment_choice_set$experiment[i]
  cs_i =  unique_experiment_choice_set$folder[i]
  
  indices_i = which(data_dmela_sample$folder == cs_i & data_dmela_sample$experiment == exp_i)
  out = tibble(
    choice_set = cs_i,
    start = indices_i[1],
    end = indices_i[length(indices_i)],
    exp = exp_i
  )
}) %>%
  bind_rows()

n_betas_level_2 = nrow(index_dataframe)

unique_experiment_choice_set_image = data_dmela_sample %>%
  select(experiment, folder, image) %>%
  distinct()

index_dataframe_image = lapply(1:nrow(unique_experiment_choice_set_image), function(i){
  exp_i = unique_experiment_choice_set_image$experiment[i]
  cs_i =  unique_experiment_choice_set_image$folder[i]
  image_i =  unique_experiment_choice_set_image$image[i]
  
  indices_i = which(data_dmela_sample$folder == cs_i & data_dmela_sample$experiment == exp_i & data_dmela_sample$image == image_i)
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


X_main <- data_dmela_sample %>% ungroup() %>%
  select(c('R','O','Y','G','B','P','UV','intensity','no_choice','is_right', 'is_daytime')) %>%
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


# X_stan <- X_main %>% 
#   select(-UV) %>% 
#   as.matrix()



X_stan <- X_main %>% 
  select(-UV) %>% 
  bind_cols(X_color_pairwise) %>% 
  bind_cols(X_color_intensity) %>%
  as.matrix()


n_experiments = max(data_dmela_sample$experiment)
n_extra_vars = 2
n_mixture_cols = dim(X_stan)[2] - n_extra_vars - 1



stan_data <- list(
  n_alt = 3, 
  n_exp = n_experiments,
  n_betas_level_2 = n_betas_level_2,
  n_var = ncol(X_stan), 
  n_obs = nrow(X_stan),
  Ny = data_dmela_sample$Ny,
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



# Stan model ------------
model_stan_lkj <- stan(
  file = here("code/european_flies_06_hmnl_2_levels.stan"),
  data = stan_data, 
  seed = 2023,
  iter = 1500,  warmup = 1000, chains = 1, cores = 4)
# Chain 1:  Elapsed Time: 2617.93 seconds (Warm-up)
# Chain 1:                1365.87 seconds (Sampling)
# Chain 1:                3983.8  seconds (Total)


summary(model_stan_lkj, pars = c("alpha"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_lkj, pars = c("beta_level_0"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_lkj, pars = c("beta_level_1"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan_lkj, pars = c("beta_level_2"), probs = c(0.1, 0.5, 0.9))$summary



plot(model_stan_lkj, plotfun="trace", pars=("beta_level_0"))
plot(model_stan_lkj, plotfun="trace", pars=("beta_level_1"))
plot(model_stan_lkj, plotfun="trace", pars=("beta_level_2"))









