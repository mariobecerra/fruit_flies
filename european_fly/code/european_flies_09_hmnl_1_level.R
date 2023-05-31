library(tidyverse)
library(rstan)
library(here)

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)

source(here("european_fly/code/european_flies_09_hmnl_1_level_utils.R"))

counts_european = read_csv(here("european_fly/out/counts_european_choice_format.csv")) %>% 
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
  slice_sample(n = 50) %>% 
  arrange(experiment, choice_set, image) %>% 
  ungroup()


counts_european_sample <- images_indices %>% 
  left_join(counts_european) %>% 
  as.data.frame() %>% 
  mutate(time_day = as.numeric(substr(as.character(date_time), 12, 13))) %>% 
  mutate(
    is_right = ifelse(alternative == "2_right", 1, 0),
    is_daytime = ifelse(time_day >= 9 & time_day <= 20, 1, 0) # 1 = 9-20h, 0 = 21-23h + 0-8h
  ) 


# Stan data ------------
unique_experiment_choice_set = counts_european_sample %>%
  select(experiment, folder) %>%
  distinct()


index_dataframe = lapply(1:nrow(unique_experiment_choice_set), function(i){
  exp_i = unique_experiment_choice_set$experiment[i]
  cs_i =  unique_experiment_choice_set$folder[i]
  
  indices_i = which(counts_european_sample$folder == cs_i & counts_european_sample$experiment == exp_i)
  out = tibble(
    choice_set = cs_i,
    start = indices_i[1],
    end = indices_i[length(indices_i)],
    exp = exp_i
  )
}) %>%
  bind_rows()


unique_experiment_choice_set_image = counts_european_sample %>%
  select(experiment, folder, image) %>%
  distinct()

index_dataframe_image = lapply(1:nrow(unique_experiment_choice_set_image), function(i){
  exp_i = unique_experiment_choice_set_image$experiment[i]
  cs_i =  unique_experiment_choice_set_image$folder[i]
  image_i =  unique_experiment_choice_set_image$image[i]
  
  indices_i = which(counts_european_sample$folder == cs_i & counts_european_sample$experiment == exp_i & counts_european_sample$image == image_i)
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



df_to_fit <- counts_european_sample %>% 
  select(all_of(c('R','O','Y','G','B','P','UV','intensity', 'no_choice'))) %>%
  as.data.frame()


X_stan_list = create_model_matrix_second_order_scheffe(df_to_fit)



n_experiments = max(counts_european_sample$experiment)
n_mixture_cols = dim(X_stan_list$X)[2] - 1



stan_data <- list(
  n_alt = 3, 
  n_exp = n_experiments,
  n_var = ncol(X_stan_list$X), 
  n_obs = nrow(X_stan_list$X),
  Ny = counts_european_sample$count,
  X = X_stan_list$X,
  experiment = index_dataframe$exp,
  choice_set = index_dataframe$choice_set,
  
  start = index_dataframe_image$start, # the starting observation for each image
  end = index_dataframe_image$end,  # the ending observation for each image
  n_images = nrow(index_dataframe_image),
  exp_index = index_dataframe_image$exp,
  
  n_mixture_cols = n_mixture_cols
)





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






Sys.time()
# 1158 seconds with 100 iterations and 10 images
# 1000 transitions using 10 leapfrog steps per transition would take 176 to 210 seconds


# 27 hours with 5000 iterations and 50 images (4000 warmup), 6 chains in 6 cores.
# 10 divergent transitions
# Although it was running at the same time of another model with 6 cores, so might have been slower because of that.
# Chain 1: Gradient evaluation took 0.080809 seconds
# Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 808.09 seconds.
# 
# Chain 2: Gradient evaluation took 0.103633 seconds
# Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1036.33 seconds.
# 
# Chain 3: Gradient evaluation took 0.070737 seconds
# Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 707.37 seconds.
# 
# Chain 4: Gradient evaluation took 0.071202 seconds
# Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 712.02 seconds.
# 
# Chain 5: Gradient evaluation took 0.115092 seconds
# Chain 5: 1000 transitions using 10 leapfrog steps per transition would take 1150.92 seconds.
# 
# Chain 6: Gradient evaluation took 0.105459 seconds
# Chain 6: 1000 transitions using 10 leapfrog steps per transition would take 1054.59 seconds.
model_stan_01 <- stan(
  file = here("european_fly/code/european_flies_09_hmnl_1_level_01.stan"),
  data = stan_data,
  seed = 2023,
  iter = 5000,  warmup = 4000, chains = 6, cores = 6,
  save_warmup = F,
  # iter = 2500,  warmup = 2000, chains = 4, cores = 4,
  # iter = 3000,  warmup = 2000, chains = 6, cores = 6,
  # iter = 100, chains = 6, cores = 6,
  # iter = 25, chains = 1,
  init = init_fun
)

Sys.time()
saveRDS(model_stan_01, "european_fly/out/09_hmnl_1_level_01_stan_object_5000_iter_50_images.rds")
model_stan_01








# 534 seconds with 50 iterations
# Gradient evaluation took 0.0181 seconds
# 1000 transitions using 10 leapfrog steps per transition would take 181 to 185 seconds


# ??? hours with 5000 iterations and 50 images (4000 warmup), 6 chains in 6 cores.
# 10 divergent transitions
# Chain 1: Gradient evaluation took 0.048389 seconds
# Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 483.89 seconds.
# 
# Chain 2: Gradient evaluation took 0.05093 seconds
# Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 509.3 seconds.
# 
# Chain 3: Gradient evaluation took 0.051514 seconds
# Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 515.14 seconds.
# 
# Chain 4: Gradient evaluation took 0.064938 seconds
# Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 649.38 seconds.
# 
# Chain 5: Gradient evaluation took 0.065467 seconds
# Chain 5: 1000 transitions using 10 leapfrog steps per transition would take 654.67 seconds.
# 
# Chain 6: Gradient evaluation took 0.065957 seconds
# Chain 6: 1000 transitions using 10 leapfrog steps per transition would take 659.57 seconds.

model_stan_03 <- stan(
  file = here("european_fly/code/european_flies_09_hmnl_1_level_03.stan"),
  data = stan_data,
  seed = 2023,
  iter = 5000,  warmup = 4000, chains = 6, cores = 6,
  save_warmup = F,
  # iter = 2500,  warmup = 2000, chains = 4, cores = 4,
  # iter = 3000,  warmup = 2000, chains = 6, cores = 6,
  # iter = 100, chains = 6, cores = 6,
  # iter = 25, chains = 1,
  init = init_fun
)

Sys.time()
saveRDS(model_stan_03, "european_fly/out/09_hmnl_1_level_03_stan_object_5000_iter_50_images.rds")
Sys.time()

model_stan_03













# 1153 seconds with 100 iterations
# Gradient evaluation took 0.0209 seconds
# 1000 transitions using 10 leapfrog steps per transition would take 209 to 220 seconds
model_stan_02 <- stan(
  file = here("european_fly/code/european_flies_09_hmnl_1_level_02.stan"),
  data = stan_data,
  seed = 2023,
  # iter = 1500,  warmup = 1000, chains = 4, cores = 4,
  # iter = 2500,  warmup = 2000, chains = 4, cores = 4,
  # iter = 3000,  warmup = 2000, chains = 6, cores = 6,
  iter = 100,  chains = 4, cores = 4,
  # iter = 25, chains = 1,
  init = init_fun
)

model_stan_02








































betas_level_0_summary_01 = summary(model_stan_01, pars = c("beta_level_0"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_stan_list$X),
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = T))

betas_level_0_summary_01 %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value")



betas_level_0_summary_01 %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = F)) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))




names_betas_level_1 = X_stan_list$names_betas_level_1


betas_level_1_summary_01 = summary(model_stan_01, pars = c("beta_level_1"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
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
  mutate(variable = fct_reorder(variable, ix, .desc = F)) 



betas_level_1_summary_01 %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  facet_wrap(~exp) +
  coord_flip() +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value")



betas_level_1_summary_01 %>% 
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



Sigma_level_0_posterior_median_01 = matrix(as.data.frame(summary(model_stan_01, pars = c("Sigma_level_0"), probs = c(0.5))$summary)$`50%`, 
                                           ncol = stan_data$n_mixture_cols+1)



diag(Sigma_level_0_posterior_median_01)
























betas_level_0_summary_02 = summary(model_stan_02, pars = c("beta_level_0"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_stan_list$X),
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = T))

betas_level_0_summary_02 %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value")



betas_level_0_summary_02 %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = F)) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



names_betas_level_1 = X_stan_list$names_betas_level_1


betas_level_1_summary_02 = summary(model_stan_02, pars = c("beta_level_1"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
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
  mutate(variable = fct_reorder(variable, ix, .desc = F)) 



betas_level_1_summary_02 %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  facet_wrap(~exp) +
  coord_flip() +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value")



betas_level_1_summary_02 %>% 
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



Sigma_level_0_posterior_median_02 = matrix(as.data.frame(summary(model_stan_02, pars = c("Sigma_level_0"), probs = c(0.5))$summary)$`50%`, 
                                           ncol = stan_data$n_mixture_cols+1)



diag(Sigma_level_0_posterior_median_02)

















betas_level_0_summary_03 = summary(model_stan_03, pars = c("beta_level_0"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_stan_list$X),
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F)) 

betas_level_0_summary_03 %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value")



betas_level_0_summary_03 %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = F)) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



names_betas_level_1 = X_stan_list$names_betas_level_1


betas_level_1_summary_03 = summary(model_stan_03, pars = c("beta_level_1"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
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
  mutate(variable = fct_reorder(variable, ix, .desc = F)) 



betas_level_1_summary_03 %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  facet_wrap(~exp) +
  coord_flip() +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value")



betas_level_1_summary_03 %>% 
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



Sigma_level_0_posterior_median_03 = matrix(as.data.frame(summary(model_stan_02, pars = c("Sigma_level_0"), probs = c(0.5))$summary)$`50%`, 
                                           ncol = stan_data$n_mixture_cols+1)



diag(Sigma_level_0_posterior_median_03)















betas_level_0_summary_01 %>% 
  mutate(model = "m_01") %>% 
  bind_rows(
    betas_level_0_summary_02 %>% 
      mutate(model = "m_02")
  ) %>% 
  bind_rows(
    betas_level_0_summary_03 %>% 
      mutate(model = "m_03")
  ) %>% 
  ggplot(aes(x = variable, color = model)) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 0.9, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), size = 0.5, position = position_dodge(width = 0.5)) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
