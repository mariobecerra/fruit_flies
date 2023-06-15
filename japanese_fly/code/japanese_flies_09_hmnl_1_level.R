library(tidyverse)
library(rstan)
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

n_experiments = max(counts_japanese_sample$experiment)

model_01_dataframe_filename = here(paste0("japanese_fly/out/japanese_hmnl_model_experiments_1_to_", n_experiments, "_dataframe_", n_images_per_cs_to_sample, "images.rds"))

saveRDS(counts_japanese_sample, model_01_dataframe_filename)


# Stan data ------------
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


X_stan_list = create_model_matrix_second_order_scheffe(df_to_fit)




n_mixture_cols = dim(X_stan_list$X)[2] - 1



european_fly_betas_level_0_mean_sd = read_csv(here("japanese_fly/data/european_fly_betas_level_0_mean_sd.csv"))

prior_tibble = tibble(
  variable = colnames(X_stan_list$X)
) %>% 
  left_join(european_fly_betas_level_0_mean_sd)


stan_data <- list(
  n_alt = 3, 
  n_exp = n_experiments,
  n_var = ncol(X_stan_list$X), 
  n_obs = nrow(X_stan_list$X),
  Ny = counts_japanese_sample$count,
  X = X_stan_list$X,
  start = index_dataframe_image$start, # the starting observation for each image
  end = index_dataframe_image$end,  # the ending observation for each image
  n_images = nrow(index_dataframe_image),
  exp_index = index_dataframe_image$exp,
  
  # prior
  tau_level_0_mean = 1,
  tau_level_0_sd = 0.5,
  L_Omega_level_0_param = 5,
  # beta_level_0_mean = rep(0, ncol(X_stan_list$X)),
  # beta_level_0_sd = rep(2, ncol(X_stan_list$X))
  beta_level_0_mean = prior_tibble$mean,
  
  # 3 times the estimated SD as a prior
  beta_level_0_sd = 3*prior_tibble$sd
)













# With 4 experiments and 50 images
# Gradient evaluation took 0.054098 seconds
# 1000 transitions using 10 leapfrog steps per transition would take 540.98 seconds

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



# With 7 experiments and 30 images per choice set:
# Gradient evaluation took 0.056564 seconds
# 1000 transitions using 10 leapfrog steps per transition would take 565.64 seconds.
# 10 hours for 600 iterations (300 warmup)
# No warnings about treedepth or divergent transitions

model_01_filename = here(paste0("japanese_fly/out/japanese_hmnl_model_experiments_1_to_", n_experiments, "_stan_object_", n_images_per_cs_to_sample, "images.rds"))

if(file.exists(model_01_filename)){
  cat("File already exists")
  model_stan_01 = readRDS(model_01_filename)
} else{
  Sys.time()
  model_stan_01 <- stan(
    file = here("japanese_fly/code/japanese_flies_hmnl_1_level_reparameterized.stan"),
    data = stan_data,
    seed = 2023,
    include = F, 
    pars = c("z", "L_Omega_level_0", "tau_level_0"),
    # control = list(max_treedepth = 10),
    # iter = 1500,  warmup = 1000, chains = 4, cores = 4,
    # iter = 2500,  warmup = 2000, chains = 4, cores = 4,
    # iter = 3000,  warmup = 2000, chains = 6, cores = 6,
    # iter = 500,  chains = 4, cores = 4,
    # iter = 5, chains = 1,
    # iter = 50, chains = 4, cores = 4,
    iter = 600,  warmup = 300, chains = 4, cores = 4,
    init = init_fun
  )
  Sys.time()
  saveRDS(model_stan_01, model_01_filename)
}

model_stan_01










betas_level_0_summary_01 = summary(model_stan_01, pars = c("beta_level_0"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_stan_list$X),
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F)) 




betas_level_0_summary_01 %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = F)) %>% 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle(paste0("Japanese fly hyper betas, experiments 1 to ", n_experiments))




names_betas_level_1 = X_stan_list$names_betas_level_1


betas_level_1_summary_01 = summary(model_stan_01, pars = c("beta_level_1"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean","sd", "2.5%", "10%", "90%", "97.5%") %>% 
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
  geom_hline(yintercept = 0) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  facet_wrap(~exp) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  ggtitle(paste0("Japanese fly betas of experiments 1 to ", n_experiments)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



betas_level_1_summary_01 %>% 
  ggplot(aes(x = variable, color = exp)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.6), size = 0.9) +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), size = 0.8, position = position_dodge(width = 0.6)) +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), size = 1.2, position = position_dodge(width = 0.6)) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  geom_vline(xintercept = 1:36 + 0.5, size = 0.3, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle(paste0("Japanese fly betas of experiments 1 to ", n_experiments))


betas_level_1_summary_01 %>% 
  ggplot(aes(x = variable, color = exp)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.6), size = 0.9) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), size = 0.8, position = position_dodge(width = 0.6)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 1.2, position = position_dodge(width = 0.6)) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  geom_vline(xintercept = 1:36 + 0.5, size = 0.3, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle(paste0("Japanese fly betas of experiments 1 to ", n_experiments))



betas_level_0_summary_01 %>% 
  mutate(exp = "All experiments", line_thickness_1 = 2, line_thickness_2 = 1) %>% 
  bind_rows(betas_level_1_summary_01 %>% 
              mutate(line_thickness_1 = 0.6, line_thickness_2 = 0.4)) %>% 
  ggplot(aes(x = variable, group = exp)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.7), size = 0.9) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd, size = line_thickness_2), 
                 position = position_dodge(width = 0.7)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd, size = line_thickness_1),
                 position = position_dodge(width = 0.7)) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  geom_vline(xintercept = 1:36 + 0.5, size = 0.3, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle(paste0("Japanese fly betas of experiments 1 to ", n_experiments)) +
  scale_size_continuous(range = c(0.5, 2.5))



betas_level_0_summary_01 %>% 
  mutate(exp = "All experiments", line_thickness_1 = 2, line_thickness_2 = 1) %>% 
  bind_rows(betas_level_1_summary_01 %>% 
              mutate(line_thickness_1 = 0.8, line_thickness_2 = 0.5)) %>% 
  ggplot(aes(x = variable, color = exp)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.7), size = 0.9) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd, size = line_thickness_2), 
                 position = position_dodge(width = 0.7), show.legend = FALSE) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd, size = line_thickness_1),
                 position = position_dodge(width = 0.7), show.legend = FALSE) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  geom_vline(xintercept = 1:36 + 0.5, size = 0.3, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle(paste0("Japanese fly betas of experiments 1 to ", n_experiments)) +
  scale_size_continuous(range = c(0.5, 2)) +
  scale_color_manual(values = c("red", "#381532", "#4b1b42", "#5d2252", "#702963",
                                "#833074", "#953784", "#a83e95"))



Sigma_level_0_posterior_median_01 = matrix(
  as.data.frame(summary(model_stan_01, pars = c("Sigma_level_0"), probs = c(0.5))$summary)$`50%`, 
  ncol = stan_data$n_var)



diag(Sigma_level_0_posterior_median_01)

























