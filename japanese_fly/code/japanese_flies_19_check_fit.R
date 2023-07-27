library(rstan)
library(MASS)
library(tidyverse)
library(bayesplot)
library(here)

source(here("japanese_fly/code/japanese_flies_hmnl_1_level_utils.R"))

# Update to correct one
# model_stan = readRDS(here("japanese_fly/out/japanese_hmnl_model_experiments_1_to_4_stan_object_30images.rds"))
model_stan = readRDS(here("japanese_fly/out/japanese_hmnl_model_experiments_1_to_9_stan_object_30images.rds"))

data_model = readRDS(here("japanese_fly/out/japanese_hmnl_model_experiments_1_to_9_dataframe_30images.rds"))

beta_level_0_draws_hmnl = rstan::extract(model_stan, "beta_level_0")
beta_level_1_draws_hmnl = rstan::extract(model_stan, "beta_level_1")
Sigma_level_0_draws_hmnl = rstan::extract(model_stan, "Sigma_level_0")
n_draws = nrow(beta_level_0_draws_hmnl[[1]])



names_betas_level_0 = c("R", "O", "Y", "G", "B", "P", "R*O", "R*Y", "R*G", "R*B", "R*P", "R*UV", "O*Y", "O*G", "O*B", "O*P", "O*UV", "Y*G", "Y*B", "Y*P", "Y*UV", "G*B", "G*P", "G*UV", "B*P", "B*UV", "P*UV", "R*intensity", "O*intensity", "Y*intensity", "G*intensity", "B*intensity", "P*intensity", "UV*intensity", "intensity^2", "no_choice")

betas_level_0_summary = summary(model_stan, pars = c("beta_level_0"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = names_betas_level_0,
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = T))



n_experiments = dim(beta_level_1_draws_hmnl[[1]])[2]

betas_level_1_summary = summary(model_stan, pars = c("beta_level_1"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean","sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(
    exp = paste0(
      "Exp ",
      rep(1:n_experiments, each = length(names_betas_level_0))
    ),
    variable = rep(names_betas_level_0, n_experiments)
  ) %>% 
  group_by(exp) %>% 
  mutate(ix = 1:n()) %>%
  ungroup() %>% 
  # mutate(variable = paste0(exp, ", ", var)) %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = F))






create_level_1_betas_hmnl = function(
  beta_level_0_draws_hmnl, Sigma_level_0_draws_hmnl, n_draws_per_posterior_sample = 1,
  seed = NULL){
  
  n_draws = nrow(beta_level_0_draws_hmnl[[1]])
  
  # Proper predictive
  betas_level_1_hmnl = matrix(0.0, nrow = n_draws_per_posterior_sample*n_draws, ncol = ncol(beta_level_0_draws_hmnl[[1]]))
  
  if(!is.null(seed)) set.seed(seed = seed)
  
  counter = 1
  for(d in 1:n_draws){
    Sigma_level_0_d = Sigma_level_0_draws_hmnl[[1]][d, , ]
    beta_level_0_d = beta_level_0_draws_hmnl[[1]][d, ]
    beta_level_1_d = mvrnorm(n_draws_per_posterior_sample, beta_level_0_d, Sigma_level_0_d)
    
    betas_level_1_hmnl[counter:(counter + n_draws_per_posterior_sample - 1), ] = beta_level_1_d
    
    counter = counter + n_draws_per_posterior_sample
  }
  
  return(betas_level_1_hmnl)
  
}









#########################################
####### Predictions in data used to fit
#########################################


predict_n_flies_original_data = function(X, Ny, data_used_to_fit, beta_level_1_draws_hmnl, n_flies_opt = 80, J = 3){
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
      
      
      utilities_d = X[exp_indices, ] %*% beta_level_1_d_ex
      
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
      Ny = Ny[exp_indices]
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




df_to_fit <- data_model %>% 
  select(all_of(c('R','O','Y','G','B','P','UV','intensity', 'no_choice'))) %>%
  as.data.frame()


X_stan_list_01 = create_model_matrix_second_order_scheffe(df_to_fit)




back_predictions = predict_n_flies_original_data(X_stan_list_01$X, data_model$count, data_model, beta_level_1_draws_hmnl, n_flies_opt = 80, J = 3)

glimpse(back_predictions)













n_y_rep_pp_check = 1000
n_y_rep_pp_check_ix = sample(1:dim(back_predictions[[1]]$N_flies_post)[1], size = n_y_rep_pp_check)
pp_check_list = lapply(1:(3*length(back_predictions)), function(x) x)
counter = 1
for(j in 1:3){
  
  if(j == 1) side = "Left"
  if(j == 2) side = "Right"
  if(j == 3) side = "No choice"
  
  for(i in 1:length(back_predictions)){
    
    out = pp_check(
      back_predictions[[i]]$design_df_preds$Ny[seq(j, length(back_predictions[[i]]$design_df_preds$Ny), by = 3)],
      back_predictions[[i]]$N_flies_post[n_y_rep_pp_check_ix, j, ], 
      ppc_dens_overlay
    ) +
      ggtitle(paste0("Experiment ", i), subtitle = side) +
      theme_bw()
    
    pp_check_list[[counter]] = out
    counter = counter + 1
  }
}

# cowplot::plot_grid(plotlist = pp_check_list, ncol = 3)
# ppc_dens_overlay_plots = gridExtra::grid.arrange(grobs = pp_check_list, ncol = length(back_predictions))
ppc_dens_overlay_plots = gridExtra::arrangeGrob(grobs = pp_check_list, ncol = length(back_predictions))


ggsave(
  here("japanese_fly/out/plots/fly_counts_ppc_dens_overlay.png"), 
  plot = ppc_dens_overlay_plots, 
  width = 70, height = 20, units = "cm", dpi = 300
  )

rm(pp_check_list)
rm(ppc_dens_overlay_plots)




pp_check_list_2 = lapply(1:(3*length(back_predictions)), function(x) x)
counter = 1
for(j in 1:3){
  
  if(j == 1) side = "Left"
  if(j == 2) side = "Right"
  if(j == 3) side = "No choice"
  
  for(i in 1:length(back_predictions)){
    
    y_rep_counts_exp_i = lapply(1:dim(back_predictions[[i]]$N_flies_post)[1], function(k){
      table(back_predictions[[i]]$N_flies_post[k, j, ]) %>% 
        as.data.frame() %>% 
        mutate(sim = k, Var1 = as.integer(Var1))
    }) %>% 
      bind_rows()
    
    y_observed_counts_exp_i = table(back_predictions[[i]]$design_df_preds$Ny[seq(j, length(back_predictions[[i]]$design_df_preds$Ny), by = 3)]) %>% 
      as.data.frame() %>% 
      mutate(group = 1, Var1 = as.integer(Var1))
    
    out = y_rep_counts_exp_i %>% 
      ggplot() +
      geom_line(
        aes(x = Var1, y = Freq, group = sim),
        color = "gray", size = 0.3, alpha = 0.2) +
      theme_bw() +
      xlab("Fly counts") +
      ylab("Frequency") +
      geom_line(
        data = y_observed_counts_exp_i,
        aes(x = Var1, y = Freq, group = group),
        color = "black", size = 0.8, alpha = 1,
        inherit.aes = F
      ) +
      ggtitle(paste0("Experiment ", i), subtitle = side)
    
    pp_check_list_2[[counter]] = out
    counter = counter + 1
  }
}


ppc_lines_overlay_plots = gridExtra::arrangeGrob(grobs = pp_check_list_2, ncol = length(back_predictions))


ggsave(
  here("japanese_fly/out/plots/fly_counts_ppc_lines_overlay.png"), 
  plot = ppc_lines_overlay_plots, 
  width = 70, height = 20, units = "cm", dpi = 300
)

rm(pp_check_list_2)














pp_check_list_3 = lapply(1:(3*length(back_predictions)), function(x) x)
counter = 1
for(j in 1:3){
  
  if(j == 1) side = "Left"
  if(j == 2) side = "Right"
  if(j == 3) side = "No choice"
  
  for(i in 1:length(back_predictions)){
    
    out = pp_check(
      back_predictions[[i]]$design_df_preds$Ny[seq(j, length(back_predictions[[i]]$design_df_preds$Ny), by = 3)],
      back_predictions[[i]]$N_flies_post[, j, ], 
      ppc_bars
    ) +
      ggtitle(paste0("Experiment ", i), subtitle = side) +
      theme_bw()
    
    pp_check_list_3[[counter]] = out
    counter = counter + 1
  }
}

ppc_bars_plots = gridExtra::arrangeGrob(grobs = pp_check_list_3, ncol = length(back_predictions))


ggsave(
  here("japanese_fly/out/plots/fly_counts_ppc_bars.png"), 
  plot = ppc_bars_plots, 
  width = 130, height = 50, units = "cm", dpi = 300, limitsize = F
)

rm(pp_check_list_3)





