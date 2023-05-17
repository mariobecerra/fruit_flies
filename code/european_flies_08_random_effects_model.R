library(tidyverse)
library(rstan)
library(here)

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)

source(here("code/european_flies_06_hmnl_1_level_utils.R"))

counts_european = read_csv(here("out/counts_european_choice_format.csv")) %>% 
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

saveRDS(counts_european_sample, here("out/european_flies_random_effects_sample.rds"))


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
  select(all_of(c('R','O','Y','G','B','P','UV','intensity','is_right', 'no_choice'))) %>%
  as.data.frame()


X_stan_list = create_model_matrix_second_order_scheffe(df_to_fit)



n_experiments = max(counts_european_sample$experiment)
# n_extra_vars = 2
n_extra_vars = 1
n_mixture_cols = dim(X_stan_list$X)[2] - n_extra_vars - 1


Z = matrix(0.0, nrow = nrow(X_stan_list$X), ncol = n_experiments)
for(i in 1:nrow(index_dataframe_image)){
  Z[index_dataframe_image$start[i]:index_dataframe_image$end[i], index_dataframe_image$exp[i]] = 1
}


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
  
  n_mixture_cols = n_mixture_cols,
  n_extra_vars = n_extra_vars,
  Z = Z
)



# Stan model ------------

init_fun <- function() {
  out_list = list(
    # tau_level_0 = rep(1, n_mixture_cols+1), 
    # Omega_level_0 = diag(1, ncol = n_mixture_cols+1, nrow = n_mixture_cols+1)
    # tau_level_0 = rnorm(n_mixture_cols+1, mean = 1, sd = 0.1), 
    tau_level_0 = rgamma(n_experiments, shape = 4, scale = 0.25),
    Omega_level_0 = rlkjcorr(K = n_experiments, eta = 30)
  )
  
  return(out_list)
}


Sys.time()

# about 1 hour with 1000 iterations and 5 images per cs, 312 divergent transitions
# 5 hours with 4000 iterations and 20 images, 423 divergent transitions
# 8 hours with 6000 iterations, 6 cores, 50 images per cs, 1468 divergent transitions. Ran out of RAM or something because got error of not being able to allocate vector oZ 6 Gb.
model_stan <- stan(
  file = here("code/european_flies_08_random_effects_model.stan"),
  data = stan_data,
  seed = 2023,
  # iter = 1000,  warmup = 500, chains = 4, cores = 4,
  # iter = 2500,  warmup = 2000, chains = 4, cores = 4,
  iter = 4000,  warmup = 3000, chains = 6, cores = 6,
  # iter = 25, chains = 1,
  init = init_fun,
  save_warmup = F 
)

Sys.time()

saveRDS(model_stan, here("out/european_flies_random_effects_model_stan_object.rds"))

Sys.time()

# model_stan <- stan(
#   file = here("code/european_flies_05_hmnl_1_level.stan"),
#   data = stan_data, 
#   seed = 2023,
#   iter = 1500,  warmup = 1000, chains = 4, cores = 4)






summary(model_stan, pars = c("alpha"), probs = c(0.1, 0.5, 0.9))$summary
summary(model_stan, pars = c("beta_level_0"), probs = c(0.1, 0.5, 0.9))$summary
# summary(model_stan, pars = c("beta_level_1"), probs = c(0.1, 0.5, 0.9))$summary



plot(model_stan, plotfun="trace", pars=("alpha"))
plot(model_stan, plotfun="trace", pars=("beta_level_0"))
# plot(model_stan, plotfun="trace", pars=("beta_level_1"))



betas_level_0_summary = summary(model_stan, pars = c("beta_level_0", "alpha"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_stan_list$X),
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = T))

betas_level_0_summary %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  coord_flip() +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  ggtitle("Random effects model")



betas_level_0_summary %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = F)) %>% 
  ggplot(aes(x = variable)) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Random effects model")



Sigma_level_0_posterior_median = matrix(as.data.frame(summary(model_stan, pars = c("Sigma_level_0"), probs = c(0.5))$summary)$`50%`, 
                                        ncol = stan_data$n_exp)


Sigma_level_0_posterior_median
diag(Sigma_level_0_posterior_median)



utilities_all_mean = get_posterior_mean(model_stan, "utilities_all")

counts_european_sample %>% 
  mutate(utility_model = utilities_all_mean[, ncol(utilities_all_mean)]) %>% 
  filter(no_choice == 0) %>% 
  select(all_of(names(df_to_fit)), utility_model) %>% 
  distinct() %>% 
  group_by(across(all_of(c('no_choice','is_right')))) %>% 
  top_n(3, utility_model) %>% 
  as.data.frame() %>% 
  ungroup() %>% 
  arrange(desc(utility_model))





# counts_by_image = counts_european_sample %>% 
#   filter(no_choice != 1) %>% 
#   group_by(experiment, image, choice_set, date_time, is_daytime) %>% 
#   summarize(n_flies = sum(count))
# 
# 
# counts_by_image %>% 
#   ggplot(aes(date_time, n_flies, color = is_daytime)) +
#   geom_point(size = 0.6) +
#   geom_line(size = 0.4) +
#   facet_wrap(~experiment, scales = "free") +
#   theme_bw()
# 
# 
# 
# counts_by_image %>% 
#   mutate(time_day = as.numeric(substr(as.character(date_time), 12, 13))) %>% 
#   group_by(experiment,time_day) %>% 
#   summarize(avg_n_flies = sum(n_flies)/n()) %>% 
#   ggplot(aes(time_day, avg_n_flies)) +
#   geom_point(size = 0.6) +
#   geom_line(size = 0.4) +
#   facet_wrap(~experiment, scales = "free") +
#   theme_bw()
# 
