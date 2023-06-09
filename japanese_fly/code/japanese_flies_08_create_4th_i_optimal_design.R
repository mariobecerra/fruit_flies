library(tidyverse)
library(opdesmixr)
library(rstan)
library(here)

var_names = c("R", "O", "Y", "G", "B", "P", "R*O", "R*Y", "R*G", "R*B", "R*P", "R*UV", "O*Y", "O*G", "O*B", "O*P", "O*UV", "Y*G", "Y*B", "Y*P", "Y*UV", "G*B", "G*P", "G*UV", "B*P", "B*UV", "P*UV", "R*intensity", "O*intensity", "Y*intensity", "G*intensity", "B*intensity", "P*intensity", "UV*intensity", "intensity^2", "no_choice")


model_01_summary = readRDS(here("japanese_fly/out/japanese_max_mnl_model_experiments_1_to_7_betas_summary.rds"))

beta_values_flies_01 = tibble(
  var = var_names
) %>% 
  left_join(
    model_01_summary, by = c("var" = "variable")
  )


length_beta_flies_01 = nrow(beta_values_flies_01)

q_flies_01 = 7
J_flies_01 = 2
S_flies_01 = 60
order_flies_01 = 2
n_pv_flies_01 = 1
no_choice_flies_01 = T

set.seed(2023)

n_draws_flies_01 = 128

beta_flies_01_draws = get_halton_draws(
  beta = beta_values_flies_01$mean, 
  sd = beta_values_flies_01$sd, 
  n_draws = n_draws_flies_01)

# beta_flies_01_draws = get_correlated_halton_draws(
#   beta = beta_values_flies_01$mean, 
#   sigma = diag(beta_values_flies_01$sd^2), 
#   n_draws = 120)
# 
# beta_flies_01_draws = matrix(0.0, nrow = n_draws_flies_01, ncol = length_beta_flies_01)
# for(i in 1:length_beta_flies_01){
#   beta_flies_01_draws[, i] = as.numeric(
#     get_halton_draws(
#       beta = beta_values_flies_01$mean[i], 
#       sd = beta_values_flies_01$sd[i],
#       ndraws = n_draws_flies_01)
#   )
# }


# To see how different the means and SDs from the draws are from the theoretical values
data.frame(
  variable = beta_values_flies_01$var,
  real_mean = beta_values_flies_01$mean,
  mean_draws = apply(beta_flies_01_draws, 2, mean),
  real_sd = beta_values_flies_01$sd,
  sd_draws = apply(beta_flies_01_draws, 2, sd)
) %>% 
  mutate(
    mean_ratio = mean_draws/real_mean,
    sd_ratio = sd_draws/real_sd
  )


# Plot the real betas versus the draws
data.frame(
  variable = beta_values_flies_01$var,
  mean = beta_values_flies_01$mean,
  sd = beta_values_flies_01$sd,
  type = "real",
  ix = 1:nrow(beta_values_flies_01)
) %>% 
  bind_rows(
    data.frame(
      variable = beta_values_flies_01$var,
      mean = apply(beta_flies_01_draws, 2, mean),
      sd = apply(beta_flies_01_draws, 2, sd),
      type = "draws",
      ix = 1:nrow(beta_values_flies_01)
    )
  ) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F)) %>% 
  ggplot(aes(x = variable, color = type)) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), size = 0.6, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 1.2, position = position_dodge(width = 0.5)) +
  theme_bw() +
  xlab("Beta") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))




i_opt_design_filename_flies_01 = paste0(here("japanese_fly/out/japanese_flies_4th_i_optimal_design.rds"))

# Around 20 minutes per iteration
# Hence, 15 iterations would take around 5 hours
# And, 20 iterations would take around 6.7 hours
# 6 initial designs with 25 iterations on 4 cores (should have been 6) on the HP machine took 9.3 hours
if(file.exists(i_opt_design_filename_flies_01)){
  cat("I_B optimal design already exists.\n")
  i_opt_design_flies_01 = readRDS(i_opt_design_filename_flies_01)
} else{
  cat("Doing I_B optimal design for insect experiment.\n")
  (t1I = Sys.time())
  i_opt_design_flies_01 =  mnl_mixture_coord_exch(
    n_random_starts = 6,
    q = q_flies_01, 
    J = J_flies_01, 
    S = S_flies_01, 
    n_pv = n_pv_flies_01,
    order = order_flies_01, 
    beta = beta_flies_01_draws,
    transform_beta = F,
    opt_method = "B",
    opt_crit = "I",
    max_it = 25,
    verbose = 1,
    plot_designs = F,
    seed = 2023,
    n_cores = 4,
    save_all_designs = T,
    no_choice = T
  )
  (t2I = Sys.time())
  t2I - t1I
  
  saveRDS(i_opt_design_flies_01, i_opt_design_filename_flies_01)
}


