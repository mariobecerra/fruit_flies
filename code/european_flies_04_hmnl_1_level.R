library(tidyverse)
library(rstan)
library(here)

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)



counts_european = read_csv(here("out/counts_choice_format.csv")) %>% 
  filter(experiment != 1) %>% 
  mutate(experiment = experiment - 1)

# There are some choice sets with only one image. Check that.
n_images_per_cs = counts_european %>% 
  group_by(experiment, folder, choice_set) %>% 
  summarize(n_images = max(image))

# Take subset of images
set.seed(2023)
images_indices = counts_european %>%
  left_join(n_images_per_cs) %>% 
  filter(n_images > 1) %>% 
  select(experiment, folder, choice_set, image) %>% 
  distinct() %>% 
  group_by(experiment, folder, choice_set) %>% 
  slice_sample(n = 20) %>% 
  arrange(experiment, choice_set, image) %>% 
  ungroup()


data_dmela_sample <- images_indices %>% 
  left_join(counts_european) %>% 
  as.data.frame() %>% 
  mutate(time_day = as.numeric(substr(as.character(date_time), 12, 13))) %>% 
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

# n_betas_level_2 = nrow(index_dataframe)

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


X_other_vars_nc <- data_dmela_sample %>% 
  select(c('no_choice','is_right', 'is_daytime')) %>%
  as.data.frame()

X_main <- data_dmela_sample %>% ungroup() %>%
  select(c(mixture_variable_names,'intensity')) %>%
  as.data.frame()


X_color_pairwise <- combn(mixture_variable_names, 2, function(x) X_main[x[1]]*X_main[x[2]], simplify = FALSE) %>%
  bind_cols() %>%
  set_names(mixture_pairwise_interaction_names)


X_color_intensity <- X_main %>%
  select(c(mixture_variable_names,'intensity')) %>%
  mutate(across(c(mixture_variable_names,'intensity'), ~.*intensity)) %>%
  set_names(c(mixture_intensity_interaction_names,'intensity^2'))


# X_stan <- X_main %>% 
#   select(-UV) %>% 
#   as.matrix()



X_stan <- X_main %>% 
  select(-UV, -intensity) %>% 
  bind_cols(X_color_pairwise) %>% 
  bind_cols(X_color_intensity) %>%
  bind_cols(X_other_vars_nc) %>% 
  as.matrix()


n_experiments = max(data_dmela_sample$experiment)
n_extra_vars = 2
n_mixture_cols = dim(X_stan)[2] - n_extra_vars - 1



stan_data <- list(
  n_alt = 3, 
  n_exp = n_experiments,
  # n_betas_level_2 = n_betas_level_2,
  n_var = ncol(X_stan), 
  n_obs = nrow(X_stan),
  Ny = data_dmela_sample$count,
  X = X_stan,
  experiment = index_dataframe$exp,
  choice_set = index_dataframe$choice_set,
  
  start = index_dataframe_image$start, # the starting observation for each image
  end = index_dataframe_image$end,  # the ending observation for each image
  n_images = nrow(index_dataframe_image),
  # experiment_index_2 = index_dataframe_image$exp,
  # cs_index_2 = index_dataframe_image$choice_set_2,
  exp_index = index_dataframe_image$exp,
  
  n_mixture_cols = n_mixture_cols,
  n_extra_vars = n_extra_vars
)



# Stan model ------------

init_list = list(
  list(tau_level_0 = rep(1, n_mixture_cols+1), 
       Omega_level_0 = diag(1, ncol = n_mixture_cols+1, nrow = n_mixture_cols+1),
       tau_level_1 = rep(1, n_mixture_cols+1), 
       Omega_level_1 = diag(1, ncol = n_mixture_cols+1, nrow = n_mixture_cols+1)),
  list(tau_level_0 = rep(0.8, n_mixture_cols+1), 
       Omega_level_0 = diag(0.8, ncol = n_mixture_cols+1, nrow = n_mixture_cols+1),
       tau_level_1 = rep(0.8, n_mixture_cols+1), 
       Omega_level_1 = diag(0.8, ncol = n_mixture_cols+1, nrow = n_mixture_cols+1)),
  list(tau_level_0 = rep(0.9, n_mixture_cols+1), 
       Omega_level_0 = diag(0.9, ncol = n_mixture_cols+1, nrow = n_mixture_cols+1),
       tau_level_1 = rep(0.9, n_mixture_cols+1), 
       Omega_level_1 = diag(0.9, ncol = n_mixture_cols+1, nrow = n_mixture_cols+1)),
  list(tau_level_0 = rep(0.5, n_mixture_cols+1), 
       Omega_level_0 = diag(0.5, ncol = n_mixture_cols+1, nrow = n_mixture_cols+1),
       tau_level_1 = rep(0.5, n_mixture_cols+1), 
       Omega_level_1 = diag(0.5, ncol = n_mixture_cols+1, nrow = n_mixture_cols+1))
)

init_fun <- function() {
  out_list = list(tau_level_0 = rep(1, n_mixture_cols+1), 
       Omega_level_0 = diag(1, ncol = n_mixture_cols+1, nrow = n_mixture_cols+1),
       tau_level_1 = rep(1, n_mixture_cols+1), 
       Omega_level_1 = diag(1, ncol = n_mixture_cols+1, nrow = n_mixture_cols+1))
  
  return(out_list)
}

# 2 hours with 1500 iterations
model_stan <- stan(
  file = here("code/european_flies_04_hmnl_1_level.stan"),
  data = stan_data,
  seed = 2023,
  iter = 1500,  warmup = 1000, chains = 4, cores = 4,
  init = init_fun
)

saveRDS(model_stan, here("out/european_flies_hmnl_1_level_stan_object.rds"))

# model_stan <- stan(
#   file = here("code/european_flies_05_hmnl_1_level.stan"),
#   data = stan_data, 
#   seed = 2023,
#   iter = 1500,  warmup = 1000, chains = 4, cores = 4)






summary(model_stan, pars = c("alpha"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan, pars = c("beta_level_0"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan, pars = c("beta_level_1"), probs = c(0.1, 0.5, 0.9))$summary



plot(model_stan, plotfun="trace", pars=("alpha"))
plot(model_stan, plotfun="trace", pars=("beta_level_0"))
plot(model_stan, plotfun="trace", pars=("beta_level_1"))



prior_params_summary = tibble(prior_mean = c(-2.99, -1.02, -0.25, -0.63, -2.54, -2.19, 5.41, 8.31, 2.64, 1.13, 15.71, 9.36, 3.25, 2.00, 3.65, 5.79, 2.75, -0.68, 7.72, -6.42, -4.48, 11.68, 8.30, 3.36, -7.40, 9.02, 6.10, 1.14, 3.09, 1.70, -0.21, -0.93, -0.32, -0.79, -0.04, 4.79, 0, 0)) %>%
  mutate(prior_sd = 2)


betas_level_0_summary = summary(model_stan, pars = c("beta_level_0", "alpha"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_stan),
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = T)) %>% 
  bind_cols(prior_params_summary)

betas_level_0_summary %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value")



betas_level_0_summary %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = F)) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



# betas_level_0_summary %>% 
#   mutate(variable = fct_reorder(variable, ix, .desc = F)) %>% 
#   ggplot(aes(x = variable)) +
#   geom_point(aes(y = mean), color = "blue", position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), color = "blue", position = position_dodge(width = 0.5)) +
#   geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), color = "blue", size = 1, position = position_dodge(width = 0.5)) +
#   geom_point(aes(y = prior_mean), color = "red", position = position_dodge(width = -0.5)) +
#   geom_linerange(aes(ymin = prior_mean - prior_sd, ymax = prior_mean + prior_sd), color = "red", position = position_dodge(width = -0.5)) +
#   geom_linerange(aes(ymin = prior_mean - 2*prior_sd, ymax = prior_mean + 2*prior_sd), color = "red", size = 1, position = position_dodge(width = -0.5)) +
#   coord_flip() +
#   theme_bw() +
#   xlab("Parameter") +
#   ylab("Value")



betas_level_0_summary %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = F)) %>% 
  select(variable, mean, sd) %>% 
  mutate(dist = "posterior") %>% 
  bind_rows(
    betas_level_0_summary %>% 
      select(variable, mean = prior_mean, sd = prior_sd) %>% 
      mutate(dist = "prior")
  ) %>% 
  ggplot(aes(x = variable, color = dist)) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 0.9, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), size = 0.5, position = position_dodge(width = 0.5)) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  ggtitle("Beta level 0", subtitle = "Prior and posterior") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))




names_betas_level_1 = c(mixture_variable_names_no_UV, colnames(X_color_pairwise), colnames(X_color_intensity), "no_choice")


betas_level_1_summary = summary(model_stan, pars = c("beta_level_1"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "2.5%", "10%", "90%", "97.5%") %>% 
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
  ylab("Value")



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
        panel.grid.minor.x = element_blank())



Sigma_level_0_posterior_median = matrix(as.data.frame(summary(model_stan, pars = c("Sigma_level_0"), probs = c(0.5))$summary)$`50%`, 
       ncol = stan_data$n_mixture_cols+1)

Sigma_level_1_posterior_median = matrix(as.data.frame(summary(model_stan, pars = c("Sigma_level_1"), probs = c(0.5))$summary)$`50%`, 
       ncol = stan_data$n_mixture_cols+1)


diag(Sigma_level_0_posterior_median)
diag(Sigma_level_1_posterior_median)



utilities_all_mean = get_posterior_mean(model_stan, "utilities_all")

data_dmela_sample %>% 
  mutate(utility_model = utilities_all_mean[, ncol(utilities_all_mean)]) %>% 
  filter(no_choice == 0) %>% 
  select(all_of(c('R','O','Y','G','B','P', 'UV', 'intensity', colnames(X_other_vars_nc))), utility_model) %>% 
  distinct() %>% 
  group_by(across(all_of(colnames(X_other_vars_nc)))) %>% 
  top_n(1, utility_model) %>% 
  as.data.frame()



